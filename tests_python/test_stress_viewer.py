import importlib.util
import math
import subprocess
import sys
from pathlib import Path

import brahe


def load_stress_viewer():
    path = Path(__file__).resolve().parents[1] / "examples" / "stress_viewer.py"
    spec = importlib.util.spec_from_file_location("stress_viewer", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_stress_benchmark_smoke():
    path = Path(__file__).resolve().parents[1] / "examples" / "stress_benchmark.py"
    result = subprocess.run(
        [
            sys.executable,
            str(path),
            "--samples",
            "3",
            "--seed",
            "7",
            "--worst",
            "1",
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    assert "Brahe stress preview benchmark" in result.stdout
    assert "samples: 3" in result.stdout
    assert "Timing" in result.stdout


def test_stress_worst_cases_benchmark_smoke():
    path = (
        Path(__file__).resolve().parents[1]
        / "examples"
        / "stress_worst_cases_benchmark.py"
    )
    result = subprocess.run(
        [sys.executable, str(path), "--repeats", "1"],
        check=True,
        capture_output=True,
        text=True,
    )

    assert "Brahe worst-case stress benchmark" in result.stdout
    assert "source_index=97" in result.stdout
    assert "overall" in result.stdout


def test_stress_system_default_request_builds():
    viewer = load_stress_viewer()
    system = viewer.make_stress_system()
    req = viewer.build_request(
        viewer.DEFAULT_RADIUS,
        viewer.DEFAULT_THETA,
        viewer.DEFAULT_SPEED,
        35.0,
        24,
    )

    status, trajectory = viewer.build_trajectory(system, req)

    assert status in {brahe.SolveStatus.Ok, brahe.SolveStatus.CapacityExceeded}
    assert trajectory.segments


def test_stress_initial_state_uses_tangential_velocity():
    viewer = load_stress_viewer()
    radius = 250_000.0
    theta = 1.25
    speed = 18.0

    state = viewer.initial_state(radius, theta, speed)

    assert math.isclose(math.hypot(state.r.x, state.r.y), radius)
    assert math.isclose(math.hypot(state.v.x, state.v.y), speed)
    assert abs(state.r.x * state.v.x + state.r.y * state.v.y) < 1e-9


def test_stress_preview_sampling_increases_for_many_revolutions():
    viewer = load_stress_viewer()
    system = viewer.make_stress_system()
    req = viewer.build_request(
        viewer.DEFAULT_RADIUS,
        viewer.DEFAULT_THETA,
        viewer.local_circular_speed(viewer.DEFAULT_RADIUS),
        175.0,
        24,
    )
    status, trajectory = viewer.build_trajectory(system, req)
    assert status == brahe.SolveStatus.Ok
    assert len(trajectory.segments) == 1

    count = viewer.preview_sample_count(system, trajectory.segments[0], 220)

    assert count > 10_000


def test_stress_preview_sampling_clusters_near_elliptic_periapsis():
    viewer = load_stress_viewer()
    system = viewer.make_stress_system()
    rp = 180_000.0
    ra = 1_200_000.0
    semi_major = 0.5 * (rp + ra)
    periapsis_speed = math.sqrt(
        viewer.PLANET_MU * (2.0 / rp - 1.0 / semi_major)
    )
    period = 2.0 * math.pi * math.sqrt(semi_major**3 / viewer.PLANET_MU)

    segment = brahe.Segment()
    segment.central_body = viewer.PLANET_ID
    segment.start_time = 0.0
    segment.end_time = period
    segment.end_reason = brahe.EventType.TimeLimit
    segment.initial_state = viewer.initial_state(rp, 0.0, periapsis_speed)

    times = viewer.preview_segment_times(system, segment, 220)
    deltas = [b - a for a, b in zip(times, times[1:])]

    assert len(times) >= viewer.PERIAPSIS_SAMPLES_PER_REVOLUTION
    assert min(deltas[:10]) < max(deltas) * 0.25


def test_stress_preview_sampling_clusters_near_hyperbolic_periapsis():
    viewer = load_stress_viewer()
    system = viewer.make_stress_system()
    rp = 180_000.0
    periapsis_speed = math.sqrt(2.0 * viewer.PLANET_MU / rp) * 1.08

    segment = brahe.Segment()
    segment.central_body = viewer.PLANET_ID
    segment.start_time = 0.0
    segment.end_time = 8.0 * viewer.SECONDS_PER_DAY
    segment.end_reason = brahe.EventType.TimeLimit
    segment.initial_state = viewer.initial_state(rp, 0.0, periapsis_speed)

    times = viewer.preview_segment_times(system, segment, 220)
    deltas = [b - a for a, b in zip(times, times[1:])]

    assert len(times) >= viewer.HYPERBOLIC_PERIAPSIS_SAMPLES
    assert min(deltas[:10]) < max(deltas) * 0.25
