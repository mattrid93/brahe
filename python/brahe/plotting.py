from __future__ import annotations

from collections import defaultdict
from math import isfinite

from .sampling import sample_segment, sample_trajectory


def _get_pyplot():
    try:
        import matplotlib.pyplot as plt
        from matplotlib.patches import Circle
    except ImportError as exc:
        raise ImportError(
            "brahe plotting helpers require matplotlib; install with "
            "'pip install brahe[plot]'"
        ) from exc
    return plt, Circle


def _figure_axes(ax):
    plt, _ = _get_pyplot()
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure
    return fig, ax


def _body_ids(system, body_ids):
    if body_ids is not None:
        return list(body_ids)
    if hasattr(system, "body_ids"):
        return list(system.body_ids())
    return [system.root_id()]


def _body_position(system, body_id, t):
    state = system.state_in_root_frame(body_id, t)
    return state.r.x, state.r.y


def _finite_positive(value):
    return isfinite(value) and value > 0.0


def plot_body_system(
    system,
    t: float = 0.0,
    ax=None,
    body_ids=None,
    show_soi: bool = True,
    show_radius: bool = True,
    label: bool = True,
):
    """Plot body positions, radii, and SOI circles in the root frame."""
    _, Circle = _get_pyplot()
    fig, ax = _figure_axes(ax)

    for body_id in _body_ids(system, body_ids):
        body = system.get_body(body_id)
        if body is None:
            continue
        x, y = _body_position(system, body_id, t)

        ax.plot([x], [y], marker="o", linestyle="", label=f"body {body_id}")
        if label:
            ax.annotate(str(body_id), (x, y), textcoords="offset points",
                        xytext=(4, 4))

        if show_soi and _finite_positive(body.soi_radius):
            ax.add_patch(
                Circle((x, y), body.soi_radius, fill=False, linewidth=0.8,
                       alpha=0.25)
            )
        if show_radius and _finite_positive(body.radius):
            ax.add_patch(
                Circle((x, y), body.radius, fill=False, linewidth=1.2,
                       alpha=0.8)
            )

    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    return fig, ax


def plot_trajectory(
    system,
    trajectory,
    ax=None,
    frame: str = "root",
    samples_per_segment: int = 200,
    show_bodies: bool = True,
    show_events: bool = True,
    body_ids=None,
):
    """Plot a sampled trajectory, optionally overlaying bodies and events."""
    fig, ax = _figure_axes(ax)

    if show_bodies:
        t = trajectory.segments[0].start_time if trajectory.segments else 0.0
        plot_body_system(system, t=t, ax=ax, body_ids=body_ids)

    samples = sample_trajectory(
        trajectory, system, samples_per_segment=samples_per_segment, frame=frame
    )
    by_segment = defaultdict(lambda: {"x": [], "y": []})
    for x, y, segment_index in zip(
        samples["x"], samples["y"], samples["segment_index"]
    ):
        by_segment[segment_index]["x"].append(x)
        by_segment[segment_index]["y"].append(y)

    for segment_index in sorted(by_segment):
        segment = by_segment[segment_index]
        ax.plot(segment["x"], segment["y"], label=f"segment {segment_index}")

    if show_events:
        for segment_index, segment in enumerate(trajectory.segments):
            sampled = sample_segment(segment, system, samples=2, frame=frame)
            if not sampled["x"]:
                continue
            ax.plot(
                [sampled["x"][-1]],
                [sampled["y"][-1]],
                marker="x",
                linestyle="",
                label=f"{segment_index}: {segment.end_reason.name}",
            )

    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    return fig, ax
