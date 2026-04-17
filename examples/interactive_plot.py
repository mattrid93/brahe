#!/usr/bin/env python3
from __future__ import annotations

import argparse

import brahe


def make_demo_system():
    sun = brahe.BodyDef()
    sun.id = 0
    sun.parent_id = brahe.InvalidBody
    sun.mu = 1000.0
    sun.radius = 1.0
    sun.soi_radius = 200.0

    moon = brahe.BodyDef()
    moon.id = 1
    moon.parent_id = 0
    moon.mu = 1.0
    moon.radius = 0.35
    moon.soi_radius = 4.0
    moon.orbit_radius = 28.0
    moon.angular_rate = 0.18
    moon.phase_at_epoch = 0.45

    builder = brahe.BodySystemBuilder()
    builder.add_body(sun)
    builder.add_body(moon)
    status, system = builder.build()
    if status != brahe.SolveStatus.Ok:
        raise RuntimeError(f"failed to build demo system: {status!r}")
    return system


def build_request(x0, y0, vx0, vy0, end_time, max_segments):
    req = brahe.PreviewRequest()
    req.central_body = 0
    req.start_time = 0.0
    req.end_time = end_time
    req.max_segments = max_segments
    req.initial_state = brahe.State2(brahe.Vec2(x0, y0), brahe.Vec2(vx0, vy0))
    return req


def build_trajectory(system, req):
    builder = brahe.TrajectoryBuilder(system)
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
        f"last event {last.end_reason.name} at t={last.end_time:.3f}"
    )


def open_interactive_plot(args):
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Button, Slider

    system = make_demo_system()

    fig, ax = plt.subplots(figsize=(9, 7))
    fig.subplots_adjust(left=0.10, right=0.96, bottom=0.34, top=0.92)

    status_text = fig.text(0.10, 0.95, "", ha="left", va="center")
    sample_count = max(8, args.samples)

    slider_axes = {
        "x0": fig.add_axes((0.16, 0.24, 0.72, 0.025)),
        "y0": fig.add_axes((0.16, 0.20, 0.72, 0.025)),
        "vx0": fig.add_axes((0.16, 0.16, 0.72, 0.025)),
        "vy0": fig.add_axes((0.16, 0.12, 0.72, 0.025)),
        "end_time": fig.add_axes((0.16, 0.08, 0.72, 0.025)),
    }
    sliders = {
        "x0": Slider(slider_axes["x0"], "x0", 5.0, 70.0, valinit=args.x0),
        "y0": Slider(slider_axes["y0"], "y0", -35.0, 35.0, valinit=args.y0),
        "vx0": Slider(slider_axes["vx0"], "vx0", -15.0, 15.0, valinit=args.vx0),
        "vy0": Slider(slider_axes["vy0"], "vy0", -15.0, 15.0, valinit=args.vy0),
        "end_time": Slider(
            slider_axes["end_time"], "end", 1.0, 30.0, valinit=args.end_time
        ),
    }

    reset_ax = fig.add_axes((0.80, 0.015, 0.10, 0.04))
    reset_button = Button(reset_ax, "Reset")

    def redraw(_=None):
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
        ax.set_title("Brahe trajectory preview")
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
    system = make_demo_system()
    req = build_request(args.x0, args.y0, args.vx0, args.vy0, args.end_time,
                        args.max_segments)
    status, trajectory = build_trajectory(system, req)
    print(format_status(status, trajectory))


def parse_args():
    parser = argparse.ArgumentParser(
        description="Open an interactive Brahe trajectory plot."
    )
    parser.add_argument("--x0", type=float, default=12.0)
    parser.add_argument("--y0", type=float, default=-9.0)
    parser.add_argument("--vx0", type=float, default=4.5)
    parser.add_argument("--vy0", type=float, default=7.0)
    parser.add_argument("--end-time", type=float, default=14.0)
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
