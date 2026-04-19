from __future__ import annotations

from ._brahe import EventType, SolveStatus, TwoBody


def _body_mu(system, body_id: int) -> float:
    body = system.get_body(body_id)
    if body is None:
        raise ValueError(f"unknown body id {body_id}")
    return body.mu


def _linspace(start: float, stop: float, count: int) -> list[float]:
    if count <= 1:
        return [start]
    step = (stop - start) / float(count - 1)
    return [start + i * step for i in range(count)]


def sample_segment(segment, system, samples: int = 200, frame: str = "central"):
    """Sample one trajectory segment.

    Returns a dict of Python lists. If NumPy is available, callers can convert
    the values directly with ``numpy.asarray``.
    """
    if samples <= 0:
        raise ValueError("samples must be positive")
    if frame not in {"central", "root"}:
        raise ValueError("frame must be 'central' or 'root'")

    mu = _body_mu(system, segment.central_body)
    out = {
        "t": [],
        "x": [],
        "y": [],
        "vx": [],
        "vy": [],
        "central_body": [],
    }
    for t in _linspace(segment.start_time, segment.end_time, samples):
        status, state = TwoBody.propagate(
            mu, segment.initial_state, t - segment.start_time
        )
        if status != SolveStatus.Ok:
            raise RuntimeError(f"propagation failed with status {status!r}")

        if frame == "root":
            body_state = system.state_in_root_frame(segment.central_body, t)
            state.r.x += body_state.r.x
            state.r.y += body_state.r.y
            state.v.x += body_state.v.x
            state.v.y += body_state.v.y

        out["t"].append(t)
        out["x"].append(state.r.x)
        out["y"].append(state.r.y)
        out["vx"].append(state.v.x)
        out["vy"].append(state.v.y)
        out["central_body"].append(segment.central_body)
    return out


def _empty_trajectory_samples():
    return {
        "t": [],
        "x": [],
        "y": [],
        "vx": [],
        "vy": [],
        "central_body": [],
        "segment_index": [],
        "event_type": [],
    }


def _sample_segment_at_time(segment, system, t: float, frame: str):
    mu = _body_mu(system, segment.central_body)
    status, state = TwoBody.propagate(mu, segment.initial_state, t - segment.start_time)
    if status != SolveStatus.Ok:
        raise RuntimeError(f"propagation failed with status {status!r}")

    if frame == "root":
        body_state = system.state_in_root_frame(segment.central_body, t)
        state.r.x += body_state.r.x
        state.r.y += body_state.r.y
        state.v.x += body_state.v.x
        state.v.y += body_state.v.y

    return state


def sample_trajectory(trajectory, system, samples_per_segment: int = 200,
                      frame: str = "root"):
    """Sample every segment in a trajectory into plotting-friendly lists."""
    if samples_per_segment <= 0:
        raise ValueError("samples_per_segment must be positive")

    out = _empty_trajectory_samples()
    for i, segment in enumerate(trajectory.segments):
        sampled = sample_segment(segment, system, samples_per_segment, frame)
        for key in ("t", "x", "y", "vx", "vy", "central_body"):
            out[key].extend(sampled[key])
        out["segment_index"].extend([i] * len(sampled["t"]))
        out["event_type"].extend([segment.end_reason] * len(sampled["t"]))
    return out


def sample_trajectory_uniform(trajectory, system, samples: int = 200,
                              frame: str = "root"):
    """Sample a trajectory at uniformly spaced absolute times.

    Unlike :func:`sample_trajectory`, this allocates samples across the whole
    trajectory duration rather than giving every segment the same count. It is
    intended for animation playback where each frame should advance simulation
    time by a consistent amount.
    """
    if samples <= 0:
        raise ValueError("samples must be positive")
    if frame not in {"central", "root"}:
        raise ValueError("frame must be 'central' or 'root'")
    if not trajectory.segments:
        return _empty_trajectory_samples()

    start_time = trajectory.segments[0].start_time
    end_time = trajectory.segments[-1].end_time
    times = _linspace(start_time, end_time, samples)
    out = _empty_trajectory_samples()

    segment_index = 0
    last_index = len(trajectory.segments) - 1
    for t in times:
        while (segment_index < last_index and
               t > trajectory.segments[segment_index].end_time):
            segment_index += 1

        segment = trajectory.segments[segment_index]
        state = _sample_segment_at_time(segment, system, t, frame)
        out["t"].append(t)
        out["x"].append(state.r.x)
        out["y"].append(state.r.y)
        out["vx"].append(state.v.x)
        out["vy"].append(state.v.y)
        out["central_body"].append(segment.central_body)
        out["segment_index"].append(segment_index)
        out["event_type"].append(segment.end_reason)

    return out
