#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math

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


def build_request(x0, y0, vx0, vy0, end_time, max_segments):
    req = brahe.PreviewRequest()
    req.central_body = EARTH_ID
    req.start_time = 0.0
    req.end_time = end_time * SECONDS_PER_DAY
    req.max_segments = max_segments
    req.initial_state = brahe.State2(brahe.Vec2(x0, y0), brahe.Vec2(vx0, vy0))
    return req


def build_trajectory(system, req):
    tolerances = brahe.Tolerances()
    tolerances.root_epsilon = 1e-3
    tolerances.max_event_refine_iterations = 96
    builder = brahe.TrajectoryBuilder(system, tolerances)
    status, trajectory = builder.build_preview(req)
    if status not in {brahe.SolveStatus.Ok, brahe.SolveStatus.CapacityExceeded}:
        raise RuntimeError(f"trajectory build failed: {status!r}")
    return status, trajectory


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
    from matplotlib.widgets import Button, Slider

    fig, ax = plt.subplots(figsize=(9, 7))
    fig.subplots_adjust(left=0.10, right=0.96, bottom=0.42, top=0.92)

    status_text = fig.text(0.10, 0.95, "", ha="left", va="center")
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

    reset_ax = fig.add_axes((0.80, 0.015, 0.10, 0.04))
    reset_button = Button(reset_ax, "Reset")

    def redraw(_=None):
        system = make_demo_system(
            moon_mu=sliders["moon_mu"].val,
            moon_orbit_radius=sliders["moon_orbit_radius"].val,
        )
        req = build_request(
            sliders["x0"].val,
            sliders["y0"].val,
            sliders["vx0"].val,
            sliders["vy0"].val,
            sliders["end_time"].val,
            args.max_segments,
        )
        status, trajectory = build_trajectory(system, req)

        ax.clear()
        brahe.plot_trajectory(
            system,
            trajectory,
            ax=ax,
            samples_per_segment=sample_count,
            show_bodies=True,
            show_events=True,
        )
        ax.set_title("Earth-Moon free-return-style trajectory preview")
        ax.legend(loc="upper right", fontsize="small")
        status_text.set_text(format_status(status, trajectory))
        fig.canvas.draw_idle()

    def reset(_=None):
        for slider in sliders.values():
            slider.reset()

    for slider in sliders.values():
        slider.on_changed(redraw)
    reset_button.on_clicked(reset)

    redraw()
    plt.show()


def smoke_check(args):
    system = make_demo_system(
        moon_mu=args.moon_mu,
        moon_orbit_radius=args.moon_orbit_radius,
    )
    req = build_request(args.x0, args.y0, args.vx0, args.vy0, args.end_time,
                        args.max_segments)
    status, trajectory = build_trajectory(system, req)
    print(format_status(status, trajectory))


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
    parser.add_argument(
        "--smoke",
        action="store_true",
        help="build one trajectory and exit without opening a plot window",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    if args.smoke:
        smoke_check(args)
    else:
        open_interactive_plot(args)


if __name__ == "__main__":
    main()
