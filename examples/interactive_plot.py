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
EARTH_ID = 0
MOON_ID = 1

EARTH_MU = 398600.4418       # km^3 / s^2
EARTH_RADIUS = 6378.137      # km
EARTH_SOI_RADIUS = 925000.0  # km

MOON_MU = 4902.800066        # km^3 / s^2
MOON_RADIUS = 1737.4         # km
MOON_ORBIT_RADIUS = 384400.0 # km
MOON_PHASE_DEG = 130.0

LOW_EARTH_ORBIT_RADIUS = 6678.0  # km, roughly 300 km altitude
TLI_SPEED = 10.85                # km/s
DEFAULT_FAILURE_LOG = "examples/interactive_failures.jsonl"


def make_demo_system(moon_mu=MOON_MU, moon_orbit_radius=MOON_ORBIT_RADIUS):
    earth = brahe.BodyDef()
    earth.id = EARTH_ID
    earth.parent_id = brahe.InvalidBody
    earth.mu = EARTH_MU
    earth.radius = EARTH_RADIUS
    earth.soi_radius = EARTH_SOI_RADIUS

    moon = brahe.BodyDef()
    moon.id = MOON_ID
    moon.parent_id = EARTH_ID
    moon.mu = moon_mu
    moon.radius = MOON_RADIUS
    moon.orbit_radius = moon_orbit_radius
    moon.soi_radius = moon_orbit_radius * (moon_mu / EARTH_MU) ** 0.4
    moon.angular_rate = math.sqrt(EARTH_MU / moon_orbit_radius**3)
    moon.phase_at_epoch = math.radians(MOON_PHASE_DEG)

    builder = brahe.BodySystemBuilder()
    builder.add_body(earth)
    builder.add_body(moon)
    status, system = builder.build()
    if status != brahe.SolveStatus.Ok:
        raise RuntimeError(f"failed to build demo system: {status!r}")
    return system


def moon_derived_values(moon_mu, moon_orbit_radius):
    if moon_mu <= 0.0 or moon_orbit_radius <= 0.0:
        return {
            "moon_soi_radius": None,
            "moon_angular_rate": None,
        }
    return {
        "moon_soi_radius": moon_orbit_radius * (moon_mu / EARTH_MU) ** 0.4,
        "moon_angular_rate": math.sqrt(EARTH_MU / moon_orbit_radius**3),
    }


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


def make_conditions(args, x0, y0, vx0, vy0, end_time, moon_mu, moon_orbit_radius):
    conditions = {
        "x0": x0,
        "y0": y0,
        "vx0": vx0,
        "vy0": vy0,
        "end_time_days": end_time,
        "max_segments": args.max_segments,
        "samples": args.samples,
        "moon_mu": moon_mu,
        "moon_orbit_radius": moon_orbit_radius,
    }
    conditions.update(moon_derived_values(moon_mu, moon_orbit_radius))
    return conditions


def build_request(x0, y0, vx0, vy0, end_time, max_segments):
    req = brahe.PreviewRequest()
    req.central_body = EARTH_ID
    req.start_time = 0.0
    req.end_time = end_time * SECONDS_PER_DAY
    req.max_segments = max_segments
    req.initial_state = brahe.State2(brahe.Vec2(x0, y0), brahe.Vec2(vx0, vy0))
    return req


def build_request_from_conditions(conditions):
    return build_request(
        conditions["x0"],
        conditions["y0"],
        conditions["vx0"],
        conditions["vy0"],
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


def open_interactive_plot(args):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
    from matplotlib.widgets import Button, Slider

    fig, ax = plt.subplots(figsize=(9, 7))
    fig.subplots_adjust(left=0.10, right=0.96, bottom=0.42, top=0.92)

    status_text = fig.text(0.10, 0.95, "", ha="left", va="center")
    time_text = fig.text(0.10, 0.925, "", ha="left", va="center")
    sample_count = max(8, args.samples)

    slider_axes = {
        "x0": fig.add_axes((0.18, 0.32, 0.70, 0.025)),
        "y0": fig.add_axes((0.18, 0.28, 0.70, 0.025)),
        "vx0": fig.add_axes((0.18, 0.24, 0.70, 0.025)),
        "vy0": fig.add_axes((0.18, 0.20, 0.70, 0.025)),
        "end_time": fig.add_axes((0.18, 0.16, 0.70, 0.025)),
        "moon_orbit_radius": fig.add_axes((0.18, 0.12, 0.70, 0.025)),
        "moon_mu": fig.add_axes((0.18, 0.08, 0.70, 0.025)),
    }
    sliders = {
        "x0": Slider(slider_axes["x0"], "x0 km", 6400.0, 9000.0, valinit=args.x0),
        "y0": Slider(slider_axes["y0"], "y0 km", -2000.0, 2000.0, valinit=args.y0),
        "vx0": Slider(slider_axes["vx0"], "vx0 km/s", -1.0, 1.0, valinit=args.vx0),
        "vy0": Slider(slider_axes["vy0"], "vy0 km/s", 10.5, 11.2, valinit=args.vy0),
        "end_time": Slider(
            slider_axes["end_time"], "end d", 3.0, 14.0, valinit=args.end_time
        ),
        "moon_orbit_radius": Slider(
            slider_axes["moon_orbit_radius"],
            "Moon orbit km",
            320000.0,
            430000.0,
            valinit=args.moon_orbit_radius,
        ),
        "moon_mu": Slider(
            slider_axes["moon_mu"], "Moon mu", 3000.0, 7000.0, valinit=args.moon_mu
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
        for artists in animation["body_artists"]:
            for artist in artists:
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
        samples = brahe.sample_trajectory(
            trajectory, system, samples_per_segment=sample_count, frame="root"
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
            sliders["x0"].val,
            sliders["y0"].val,
            sliders["vx0"].val,
            sliders["vy0"].val,
            sliders["end_time"].val,
            sliders["moon_mu"].val,
            sliders["moon_orbit_radius"].val,
        )
        req = None
        system = None
        try:
            system = make_demo_system(
                moon_mu=conditions["moon_mu"],
                moon_orbit_radius=conditions["moon_orbit_radius"],
            )
            req = build_request_from_conditions(conditions)
            status, trajectory = build_trajectory(system, req)

            ax.clear()
            brahe.plot_trajectory(
                system,
                trajectory,
                ax=ax,
                samples_per_segment=sample_count,
                show_bodies=False,
                show_events=True,
            )
            setup_animation_overlay(system, trajectory)
            ax.set_title("Earth-Moon free-return-style trajectory preview")
            ax.legend(loc="upper right", fontsize="small")
            status_text.set_text(format_status(status, trajectory))
        except Exception as exc:
            log_failure(args.failure_log, conditions, req=req, exc=exc)
            ax.clear()
            if system is not None:
                try:
                    brahe.plot_body_system(system, ax=ax)
                except Exception:
                    pass
            ax.set_title("Trajectory build failed")
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


def smoke_check(args):
    conditions = make_conditions(
        args,
        args.x0,
        args.y0,
        args.vx0,
        args.vy0,
        args.end_time,
        args.moon_mu,
        args.moon_orbit_radius,
    )
    run_conditions(conditions, failure_log=args.failure_log)


def load_failure(path, index):
    with Path(path).open("r", encoding="utf-8") as f:
        records = [json.loads(line) for line in f if line.strip()]
    if not records:
        raise RuntimeError(f"failure log is empty: {path}")
    return records[index]


def run_conditions(conditions, failure_log=None):
    req = None
    trajectory = None
    status = None
    try:
        system = make_demo_system(
            moon_mu=conditions["moon_mu"],
            moon_orbit_radius=conditions["moon_orbit_radius"],
        )
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


def replay_failure(args):
    record = load_failure(args.replay_failure, args.failure_index)
    run_conditions(record["conditions"], failure_log=args.failure_log)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Open an interactive Brahe trajectory plot."
    )
    parser.add_argument("--x0", type=float, default=LOW_EARTH_ORBIT_RADIUS)
    parser.add_argument("--y0", type=float, default=0.0)
    parser.add_argument("--vx0", type=float, default=0.0)
    parser.add_argument("--vy0", type=float, default=TLI_SPEED)
    parser.add_argument("--end-time", type=float, default=7.0,
                        help="preview duration in days")
    parser.add_argument("--moon-orbit-radius", type=float, default=MOON_ORBIT_RADIUS)
    parser.add_argument("--moon-mu", type=float, default=MOON_MU)
    parser.add_argument("--max-segments", type=int, default=8)
    parser.add_argument("--samples", type=int, default=240)
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
