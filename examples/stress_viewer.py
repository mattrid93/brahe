#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import platform
import subprocess
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

import brahe

SECONDS_PER_DAY = 86400.0
PLANET_ID = 0
MOON_IDS = (1, 2, 3, 4)

PLANET_MU = 126_686_534.0
PLANET_RADIUS = 71_492.0
PLANET_SOI_RADIUS = 48_000_000.0

DEFAULT_RADIUS = 180_000.0
DEFAULT_THETA = 0.0
DEFAULT_SPEED = 30.0
DEFAULT_FAILURE_LOG = "examples/stress_failures.jsonl"
MIN_PREVIEW_SAMPLES_PER_SEGMENT = 120
SAMPLES_PER_REVOLUTION = 80
MAX_PREVIEW_SAMPLES_PER_SEGMENT = 20_000
PERIAPSIS_SAMPLES_PER_REVOLUTION = 260
HYPERBOLIC_PERIAPSIS_SAMPLES = 520

MOON_DEFS = (
    {
        "mu": 3_200.0,
        "radius": 1_400.0,
        "soi_radius": 24_000.0,
        "orbit_radius": 310_000.0,
        "phase": 0.20,
        "rate_scale": 1.00,
    },
    {
        "mu": 12_500.0,
        "radius": 2_600.0,
        "soi_radius": 48_000.0,
        "orbit_radius": 520_000.0,
        "phase": 1.70,
        "rate_scale": 1.08,
    },
    {
        "mu": 38_000.0,
        "radius": 4_400.0,
        "soi_radius": 82_000.0,
        "orbit_radius": 860_000.0,
        "phase": 3.20,
        "rate_scale": 0.93,
    },
    {
        "mu": 95_000.0,
        "radius": 7_200.0,
        "soi_radius": 130_000.0,
        "orbit_radius": 1_380_000.0,
        "phase": 4.60,
        "rate_scale": 1.14,
    },
)

PLOT_LIMIT = 1.2 * max(spec["orbit_radius"] for spec in MOON_DEFS)


def current_git_commit():
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        )
    except Exception:
        return None
    return result.stdout.strip()


def make_stress_system():
    planet = brahe.BodyDef()
    planet.id = PLANET_ID
    planet.parent_id = brahe.InvalidBody
    planet.mu = PLANET_MU
    planet.radius = PLANET_RADIUS
    planet.soi_radius = PLANET_SOI_RADIUS

    builder = brahe.BodySystemBuilder()
    builder.add_body(planet)

    for body_id, spec in zip(MOON_IDS, MOON_DEFS):
        moon = brahe.BodyDef()
        moon.id = body_id
        moon.parent_id = PLANET_ID
        moon.mu = spec["mu"]
        moon.radius = spec["radius"]
        moon.soi_radius = spec["soi_radius"]
        moon.orbit_radius = spec["orbit_radius"]
        moon.angular_rate = (
            math.sqrt(PLANET_MU / moon.orbit_radius**3) * spec["rate_scale"]
        )
        moon.phase_at_epoch = spec["phase"]
        builder.add_body(moon)

    status, system = builder.build()
    if status != brahe.SolveStatus.Ok:
        raise RuntimeError(f"failed to build stress system: {status!r}")
    return system


def initial_state(radius, theta, speed):
    c = math.cos(theta)
    s = math.sin(theta)
    return brahe.State2(
        brahe.Vec2(radius * c, radius * s),
        brahe.Vec2(-speed * s, speed * c),
    )


def make_conditions(args, radius, theta, speed, end_time):
    return {
        "radius": radius,
        "theta": theta,
        "speed": speed,
        "end_time_days": end_time,
        "max_segments": args.max_segments,
        "samples": args.samples,
        "moons": [
            {
                "id": body_id,
                "mu": spec["mu"],
                "radius": spec["radius"],
                "soi_radius": spec["soi_radius"],
                "orbit_radius": spec["orbit_radius"],
                "phase": spec["phase"],
                "rate_scale": spec["rate_scale"],
                "angular_rate": math.sqrt(PLANET_MU / spec["orbit_radius"]**3)
                * spec["rate_scale"],
            }
            for body_id, spec in zip(MOON_IDS, MOON_DEFS)
        ],
    }


def build_request(radius, theta, speed, end_time, max_segments):
    req = brahe.PreviewRequest()
    req.central_body = PLANET_ID
    req.start_time = 0.0
    req.end_time = end_time * SECONDS_PER_DAY
    req.max_segments = max_segments
    req.initial_state = initial_state(radius, theta, speed)
    return req


def build_request_from_conditions(conditions):
    return build_request(
        conditions["radius"],
        conditions["theta"],
        conditions["speed"],
        conditions["end_time_days"],
        conditions["max_segments"],
    )


def build_trajectory(system, req):
    tolerances = brahe.Tolerances()
    tolerances.root_epsilon = 1e-3
    tolerances.max_event_refine_iterations = 96
    builder = brahe.TrajectoryBuilder(system, tolerances)
    status, trajectory = builder.build_preview(req)
    if status not in {brahe.SolveStatus.Ok, brahe.SolveStatus.CapacityExceeded}:
        raise RuntimeError(f"trajectory build failed: {status!r}")
    return status, trajectory


def segment_records(trajectory):
    return [
        {
            "central_body": segment.central_body,
            "start_time": segment.start_time,
            "end_time": segment.end_time,
            "end_reason": segment.end_reason.name,
            "initial_state": {
                "rx": segment.initial_state.r.x,
                "ry": segment.initial_state.r.y,
                "vx": segment.initial_state.v.x,
                "vy": segment.initial_state.v.y,
            },
        }
        for segment in trajectory.segments
    ]


def request_record(req):
    return {
        "central_body": req.central_body,
        "start_time": req.start_time,
        "end_time": req.end_time,
        "max_segments": req.max_segments,
        "initial_state": {
            "rx": req.initial_state.r.x,
            "ry": req.initial_state.r.y,
            "vx": req.initial_state.v.x,
            "vy": req.initial_state.v.y,
        },
    }


def log_failure(path, conditions, req=None, status=None, trajectory=None, exc=None):
    if path is None:
        return

    record = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "git_commit": current_git_commit(),
        "conditions": conditions,
    }
    if req is not None:
        record["request"] = request_record(req)
    if status is not None:
        record["status"] = status.name
    if trajectory is not None:
        record["segments"] = segment_records(trajectory)
    if exc is not None:
        record["exception"] = {
            "type": type(exc).__name__,
            "message": str(exc),
            "traceback": traceback.format_exc(),
        }

    log_path = Path(path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as f:
        f.write(json.dumps(record, sort_keys=True) + "\n")


def format_status(status, trajectory):
    if not trajectory.segments:
        return f"{status.name}: no segments"
    last = trajectory.segments[-1]
    return (
        f"{status.name}: {len(trajectory.segments)} segment(s), "
        f"last event {last.end_reason.name} at t={last.end_time / SECONDS_PER_DAY:.3f} d"
    )


def local_circular_speed(radius):
    return math.sqrt(PLANET_MU / radius)


def linear_times(start_time, end_time, count):
    count = max(1, count)
    if count == 1:
        return [start_time]
    return [
        start_time + i * (end_time - start_time) / (count - 1)
        for i in range(count)
    ]


def solve_eccentric_anomaly(mean_anomaly, eccentricity):
    reduced = math.fmod(mean_anomaly, 2.0 * math.pi)
    if reduced > math.pi:
        reduced -= 2.0 * math.pi
    if reduced <= -math.pi:
        reduced += 2.0 * math.pi
    offset = mean_anomaly - reduced
    if reduced == 0.0:
        return offset

    sign = -1.0 if reduced < 0.0 else 1.0
    target = abs(reduced)
    lo = 0.0
    hi = math.pi
    estimate = min(max(target + eccentricity * math.sin(target), lo), hi)
    for _ in range(64):
        residual = estimate - eccentricity * math.sin(estimate) - target
        if abs(residual) < 1e-14 * max(1.0, target):
            return offset + sign * estimate
        if residual > 0.0:
            hi = estimate
        else:
            lo = estimate
        derivative = 1.0 - eccentricity * math.cos(estimate)
        if abs(derivative) < 1e-30:
            estimate = 0.5 * (lo + hi)
            continue
        next_estimate = estimate - residual / derivative
        if not (lo < next_estimate < hi) or not math.isfinite(next_estimate):
            next_estimate = 0.5 * (lo + hi)
        estimate = next_estimate
    return offset + sign * estimate


def solve_hyperbolic_anomaly(mean_anomaly, eccentricity):
    if mean_anomaly == 0.0:
        return 0.0

    sign = -1.0 if mean_anomaly < 0.0 else 1.0
    target = abs(mean_anomaly)
    lo = 0.0
    hi = max(1.0, math.log(2.0 * target / eccentricity + 1.0))

    def residual(value):
        return eccentricity * math.sinh(value) - value - target

    while residual(hi) < 0.0:
        hi *= 2.0
        if not math.isfinite(hi) or hi > 1000.0:
            return sign * hi

    estimate = 0.5 * (lo + hi)
    for _ in range(64):
        value = residual(estimate)
        if abs(value) < 1e-14 * max(1.0, target):
            return sign * estimate
        if value > 0.0:
            hi = estimate
        else:
            lo = estimate
        derivative = eccentricity * math.cosh(estimate) - 1.0
        if abs(derivative) < 1e-30:
            estimate = 0.5 * (lo + hi)
            continue
        next_estimate = estimate - value / derivative
        if not (lo < next_estimate < hi) or not math.isfinite(next_estimate):
            next_estimate = 0.5 * (lo + hi)
        estimate = next_estimate
    return sign * estimate


def preview_segment_times(system, segment, base_samples):
    base_count = max(1, base_samples)
    body = system.get_body(segment.central_body)
    if body is None:
        return linear_times(segment.start_time, segment.end_time, base_count)

    status, elements = brahe.TwoBody.to_elements(body.mu, segment.initial_state)
    if status != brahe.SolveStatus.Ok:
        return linear_times(segment.start_time, segment.end_time, base_count)

    a = elements.semi_major_axis
    if not math.isfinite(a):
        return linear_times(segment.start_time, segment.end_time, base_count)

    if elements.type == brahe.ConicType.Hyperbola:
        if a >= 0.0:
            return linear_times(segment.start_time, segment.end_time, base_count)
        duration = max(0.0, segment.end_time - segment.start_time)
        mean_motion = math.sqrt(body.mu / (-(a**3)))
        sign_h = 1.0 if elements.angular_momentum_z >= 0.0 else -1.0
        mean_start = elements.mean_anomaly
        mean_end = mean_start + sign_h * mean_motion * duration
        anomaly_start = elements.hyperbolic_anomaly
        anomaly_end = solve_hyperbolic_anomaly(mean_end, elements.eccentricity)
        count = min(
            max(base_count, MIN_PREVIEW_SAMPLES_PER_SEGMENT,
                HYPERBOLIC_PERIAPSIS_SAMPLES),
            MAX_PREVIEW_SAMPLES_PER_SEGMENT,
        )
        if count <= 1:
            return [segment.start_time]
        times = []
        for i in range(count):
            anomaly = (
                anomaly_start
                + i * (anomaly_end - anomaly_start) / (count - 1)
            )
            mean = elements.eccentricity * math.sinh(anomaly) - anomaly
            dt = sign_h * (mean - mean_start) / mean_motion
            times.append(min(max(segment.start_time + dt, segment.start_time),
                             segment.end_time))
        return times

    if elements.type != brahe.ConicType.Ellipse or a <= 0.0:
        return linear_times(segment.start_time, segment.end_time, base_count)

    period = 2.0 * math.pi * math.sqrt(a**3 / body.mu)
    duration = max(0.0, segment.end_time - segment.start_time)
    revolutions = duration / period if period > 0.0 else 0.0
    adaptive = int(math.ceil(revolutions * SAMPLES_PER_REVOLUTION)) + 1
    count = min(
        max(base_count, MIN_PREVIEW_SAMPLES_PER_SEGMENT, adaptive),
        MAX_PREVIEW_SAMPLES_PER_SEGMENT,
    )
    if count <= 1:
        return [segment.start_time]

    e = elements.eccentricity
    if e < 0.2:
        return linear_times(segment.start_time, segment.end_time, count)

    periapsis_count = min(
        max(
            count,
            int(math.ceil(max(1.0, revolutions) * PERIAPSIS_SAMPLES_PER_REVOLUTION))
            + 1,
        ),
        MAX_PREVIEW_SAMPLES_PER_SEGMENT,
    )
    n = 2.0 * math.pi / period
    sign_h = 1.0 if elements.angular_momentum_z >= 0.0 else -1.0
    mean_start = elements.mean_anomaly
    mean_end = mean_start + sign_h * n * duration
    eccentric_start = elements.eccentric_anomaly
    eccentric_end = solve_eccentric_anomaly(mean_end, e)
    times = []
    for i in range(periapsis_count):
        eccentric = (
            eccentric_start
            + i * (eccentric_end - eccentric_start) / (periapsis_count - 1)
        )
        mean = eccentric - e * math.sin(eccentric)
        dt = sign_h * (mean - mean_start) / n if n > 0.0 else 0.0
        times.append(min(max(segment.start_time + dt, segment.start_time),
                         segment.end_time))
    return times


def preview_sample_count(system, segment, base_samples):
    return len(preview_segment_times(system, segment, base_samples))


def plot_stress_trajectory(system, trajectory, ax, base_samples, show_events=True):
    for segment_index, segment in enumerate(trajectory.segments):
        times = preview_segment_times(system, segment, base_samples)
        xs = []
        ys = []
        body = system.get_body(segment.central_body)
        if body is None:
            continue
        for t in times:
            status, state = brahe.TwoBody.propagate(
                body.mu, segment.initial_state, t - segment.start_time
            )
            if status != brahe.SolveStatus.Ok:
                raise RuntimeError(f"propagation failed with status {status!r}")
            body_state = system.state_in_root_frame(segment.central_body, t)
            xs.append(state.r.x + body_state.r.x)
            ys.append(state.r.y + body_state.r.y)

        ax.plot(xs, ys, label=f"segment {segment_index}")

        if show_events and xs:
            ax.plot(
                [xs[-1]],
                [ys[-1]],
                marker="x",
                linestyle="",
                label=f"{segment_index}: {segment.end_reason.name}",
            )

    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel("x (km)")
    ax.set_ylabel("y (km)")


def open_interactive_plot(args):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
    from matplotlib.widgets import Button, Slider

    fig, ax = plt.subplots(figsize=(9, 7))
    fig.subplots_adjust(left=0.10, right=0.96, bottom=0.34, top=0.92)

    status_text = fig.text(0.10, 0.95, "", ha="left", va="center")
    time_text = fig.text(0.10, 0.925, "", ha="left", va="center")
    sample_count = max(8, args.samples)

    slider_axes = {
        "radius": fig.add_axes((0.18, 0.22, 0.70, 0.025)),
        "theta": fig.add_axes((0.18, 0.18, 0.70, 0.025)),
        "speed": fig.add_axes((0.18, 0.14, 0.70, 0.025)),
        "end_time": fig.add_axes((0.18, 0.10, 0.70, 0.025)),
    }
    sliders = {
        "radius": Slider(
            slider_axes["radius"],
            "radius km",
            PLANET_RADIUS + 1_000.0,
            1_600_000.0,
            valinit=args.radius,
        ),
        "theta": Slider(
            slider_axes["theta"], "theta rad", 0.0, 2.0 * math.pi, valinit=args.theta
        ),
        "speed": Slider(
            slider_axes["speed"], "tangent km/s", 1.0, 65.0, valinit=args.speed
        ),
        "end_time": Slider(
            slider_axes["end_time"], "end d", 1.0, 400.0, valinit=args.end_time
        ),
    }

    play_ax = fig.add_axes((0.68, 0.015, 0.10, 0.04))
    play_button = Button(play_ax, "Play")
    reset_ax = fig.add_axes((0.80, 0.015, 0.10, 0.04))
    reset_button = Button(reset_ax, "Reset")

    animation = {
        "system": None,
        "samples": None,
        "frame": 0,
        "playing": False,
        "trail": None,
        "craft": None,
        "time_marker": None,
        "body_artists": [],
    }
    timer = fig.canvas.new_timer(interval=args.animation_interval_ms)

    def stop_animation():
        animation["playing"] = False
        timer.stop()
        play_button.label.set_text("Play")

    def clear_animation_artists():
        for key in ("trail", "craft", "time_marker"):
            artist = animation.get(key)
            if artist is not None:
                artist.remove()
                animation[key] = None
        for _, marker, soi_circle, label in animation["body_artists"]:
            for artist in (marker, soi_circle, label):
                artist.remove()
        animation["body_artists"] = []
        time_text.set_text("")

    def update_animation_frame(frame):
        samples = animation["samples"]
        system = animation["system"]
        if not samples or system is None:
            return

        frame = max(0, min(frame, len(samples["t"]) - 1))
        animation["frame"] = frame

        xs = samples["x"][: frame + 1]
        ys = samples["y"][: frame + 1]
        t = samples["t"][frame]

        animation["trail"].set_data(xs, ys)
        animation["craft"].set_data([samples["x"][frame]], [samples["y"][frame]])
        time_label = f"t = {t / SECONDS_PER_DAY:.3f} d"
        animation["time_marker"].set_text(time_label)
        time_text.set_text(time_label)

        for body_id, marker, soi_circle, label in animation["body_artists"]:
            body_state = system.state_in_root_frame(body_id, t)
            x = body_state.r.x
            y = body_state.r.y
            marker.set_data([x], [y])
            soi_circle.center = (x, y)
            label.set_position((x, y))

        fig.canvas.draw_idle()

    def setup_animation_overlay(system, trajectory):
        clear_animation_artists()
        frame_count = sample_count * max(1, len(trajectory.segments))
        samples = brahe.sample_trajectory_uniform(
            trajectory, system, samples=frame_count, frame="root"
        )
        animation["system"] = system
        animation["samples"] = samples
        animation["frame"] = 0

        (trail,) = ax.plot([], [], color="black", linewidth=2.0, label="craft trail")
        (craft,) = ax.plot(
            [], [], marker="o", color="red", markersize=6, linestyle="", label="craft"
        )
        time_marker = ax.text(
            0.02,
            0.98,
            "",
            transform=ax.transAxes,
            ha="left",
            va="top",
            bbox={"facecolor": "white", "alpha": 0.75, "edgecolor": "none"},
        )
        animation["trail"] = trail
        animation["craft"] = craft
        animation["time_marker"] = time_marker

        body_artists = []
        start_time = samples["t"][0] if samples["t"] else 0.0
        for body_id in system.body_ids():
            body = system.get_body(body_id)
            if body is None:
                continue
            body_state = system.state_in_root_frame(body_id, start_time)
            x = body_state.r.x
            y = body_state.r.y
            (marker,) = ax.plot(
                [x], [y], marker="o", linestyle="", markersize=5, label=f"body {body_id}"
            )
            soi_circle = Circle(
                (x, y), body.soi_radius, fill=False, linewidth=0.8, alpha=0.18
            )
            ax.add_patch(soi_circle)
            label = ax.annotate(str(body_id), (x, y), textcoords="offset points",
                                xytext=(4, 4))
            body_artists.append((body_id, marker, soi_circle, label))
        animation["body_artists"] = body_artists
        update_animation_frame(0)

    def advance_animation():
        if not animation["playing"]:
            return
        samples = animation["samples"]
        if not samples:
            stop_animation()
            return
        next_frame = animation["frame"] + 1
        if next_frame >= len(samples["t"]):
            stop_animation()
            return
        update_animation_frame(next_frame)

    timer.add_callback(advance_animation)

    def redraw(_=None):
        stop_animation()
        clear_animation_artists()
        conditions = make_conditions(
            args,
            sliders["radius"].val,
            sliders["theta"].val,
            sliders["speed"].val,
            sliders["end_time"].val,
        )
        req = None
        system = None
        status = None
        trajectory = None
        try:
            system = make_stress_system()
            req = build_request_from_conditions(conditions)
            status, trajectory = build_trajectory(system, req)

            ax.clear()
            plot_stress_trajectory(system, trajectory, ax, sample_count)
            setup_animation_overlay(system, trajectory)
            ax.set_xlim(-PLOT_LIMIT, PLOT_LIMIT)
            ax.set_ylim(-PLOT_LIMIT, PLOT_LIMIT)
            vcirc = local_circular_speed(conditions["radius"])
            ax.set_title("Multi-moon propagation stress viewer")
            ax.legend(loc="upper right", fontsize="small")
            status_text.set_text(
                f"{format_status(status, trajectory)} | v/vcirc={conditions['speed'] / vcirc:.3f}"
            )
        except Exception as exc:
            log_failure(
                args.failure_log,
                conditions,
                req=req,
                status=status,
                trajectory=trajectory,
                exc=exc,
            )
            ax.clear()
            if system is not None:
                try:
                    brahe.plot_body_system(system, ax=ax)
                except Exception:
                    pass
            ax.set_title("Stress trajectory build failed")
            status_text.set_text(
                f"{type(exc).__name__}: {exc} "
                f"(logged to {args.failure_log})"
            )
        fig.canvas.draw_idle()

    def toggle_play(_=None):
        samples = animation["samples"]
        if not samples:
            return
        if animation["playing"]:
            stop_animation()
            return
        if animation["frame"] >= len(samples["t"]) - 1:
            update_animation_frame(0)
        animation["playing"] = True
        play_button.label.set_text("Pause")
        timer.start()

    def reset(_=None):
        stop_animation()
        for slider in sliders.values():
            slider.reset()

    for slider in sliders.values():
        slider.on_changed(redraw)
    play_button.on_clicked(toggle_play)
    reset_button.on_clicked(reset)

    redraw()
    plt.show()


def run_conditions(conditions, failure_log=None):
    req = None
    trajectory = None
    status = None
    try:
        system = make_stress_system()
        req = build_request_from_conditions(conditions)
        status, trajectory = build_trajectory(system, req)
        print(format_status(status, trajectory))
        for i, segment in enumerate(trajectory.segments):
            print(
                f"  {i}: body={segment.central_body} "
                f"event={segment.end_reason.name} "
                f"t={segment.end_time / SECONDS_PER_DAY:.6f} d"
            )
    except Exception as exc:
        log_failure(
            failure_log,
            conditions,
            req=req,
            status=status,
            trajectory=trajectory,
            exc=exc,
        )
        raise


def smoke_check(args):
    conditions = make_conditions(
        args,
        args.radius,
        args.theta,
        args.speed,
        args.end_time,
    )
    run_conditions(conditions, failure_log=args.failure_log)


def load_failure(path, index):
    with Path(path).open("r", encoding="utf-8") as f:
        records = [json.loads(line) for line in f if line.strip()]
    if not records:
        raise RuntimeError(f"failure log is empty: {path}")
    return records[index]


def replay_failure(args):
    record = load_failure(args.replay_failure, args.failure_index)
    run_conditions(record["conditions"], failure_log=args.failure_log)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Open a multi-moon Brahe trajectory stress viewer."
    )
    parser.add_argument("--radius", type=float, default=DEFAULT_RADIUS)
    parser.add_argument("--theta", type=float, default=DEFAULT_THETA)
    parser.add_argument("--speed", type=float, default=DEFAULT_SPEED)
    parser.add_argument("--end-time", type=float, default=175.0,
                        help="preview duration in days")
    parser.add_argument("--max-segments", type=int, default=24)
    parser.add_argument("--samples", type=int, default=220)
    parser.add_argument("--animation-interval-ms", type=int, default=40)
    parser.add_argument("--failure-log", default=DEFAULT_FAILURE_LOG)
    parser.add_argument(
        "--replay-failure",
        help="replay conditions from a JSONL failure log instead of using CLI sliders",
    )
    parser.add_argument("--failure-index", type=int, default=-1)
    parser.add_argument(
        "--smoke",
        action="store_true",
        help="build one trajectory and exit without opening a plot window",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    if args.replay_failure:
        replay_failure(args)
    elif args.smoke:
        smoke_check(args)
    else:
        open_interactive_plot(args)


if __name__ == "__main__":
    main()
