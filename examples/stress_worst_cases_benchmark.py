#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import statistics
import time

import brahe

SECONDS_PER_DAY = 86400.0
PLANET_ID = 0
PLANET_MU = 126_686_534.0
PLANET_RADIUS = 71_492.0
PLANET_SOI_RADIUS = 48_000_000.0

MOONS = (
    (1, 3_200.0, 1_400.0, 24_000.0, 310_000.0, 0.20, 1.00),
    (2, 12_500.0, 2_600.0, 48_000.0, 520_000.0, 1.70, 1.08),
    (3, 38_000.0, 4_400.0, 82_000.0, 860_000.0, 3.20, 0.93),
    (4, 95_000.0, 7_200.0, 130_000.0, 1_380_000.0, 4.60, 1.14),
)

CASES = (
    {
        "source_index": 97,
        "radius": 242_598.819617,
        "theta": 1.352099192,
        "speed": 21.177021862,
    },
    {
        "source_index": 48,
        "radius": 238_626.355304,
        "theta": 5.653726381,
        "speed": 17.630564246,
    },
    {
        "source_index": 14,
        "radius": 257_152.380723,
        "theta": 2.090385500,
        "speed": 24.020838391,
    },
)


def percentile(sorted_values, fraction):
    if len(sorted_values) == 1:
        return sorted_values[0]
    position = fraction * (len(sorted_values) - 1)
    lo = int(math.floor(position))
    hi = int(math.ceil(position))
    if lo == hi:
        return sorted_values[lo]
    weight = position - lo
    return sorted_values[lo] * (1.0 - weight) + sorted_values[hi] * weight


def make_system():
    planet = brahe.BodyDef()
    planet.id = PLANET_ID
    planet.parent_id = brahe.InvalidBody
    planet.mu = PLANET_MU
    planet.radius = PLANET_RADIUS
    planet.soi_radius = PLANET_SOI_RADIUS

    builder = brahe.BodySystemBuilder()
    builder.add_body(planet)
    for body_id, mu, radius, soi_radius, orbit_radius, phase, rate_scale in MOONS:
        moon = brahe.BodyDef()
        moon.id = body_id
        moon.parent_id = PLANET_ID
        moon.mu = mu
        moon.radius = radius
        moon.soi_radius = soi_radius
        moon.orbit_radius = orbit_radius
        moon.angular_rate = math.sqrt(PLANET_MU / orbit_radius**3) * rate_scale
        moon.phase_at_epoch = phase
        builder.add_body(moon)

    status, system = builder.build()
    if status != brahe.SolveStatus.Ok:
        raise RuntimeError(f"failed to build system: {status!r}")
    return system


def initial_state(radius, theta, speed):
    c = math.cos(theta)
    s = math.sin(theta)
    return brahe.State2(
        brahe.Vec2(radius * c, radius * s),
        brahe.Vec2(-speed * s, speed * c),
    )


def make_request(case, end_time_days, max_segments):
    req = brahe.PreviewRequest()
    req.central_body = PLANET_ID
    req.start_time = 0.0
    req.end_time = end_time_days * SECONDS_PER_DAY
    req.max_segments = max_segments
    req.initial_state = initial_state(case["radius"], case["theta"], case["speed"])
    return req


def make_builder(system):
    tolerances = brahe.Tolerances()
    tolerances.root_epsilon = 1e-3
    tolerances.max_event_refine_iterations = 96
    return brahe.TrajectoryBuilder(system, tolerances)


def run_case(builder, request, repeats):
    timings = []
    last_status = None
    last_trajectory = None
    for _ in range(repeats):
        start = time.perf_counter()
        status, trajectory = builder.build_preview(request)
        timings.append(time.perf_counter() - start)
        last_status = status
        last_trajectory = trajectory
    return timings, last_status, last_trajectory


def print_case_summary(case, timings, status, trajectory):
    sorted_timings = sorted(timings)
    total = sum(timings)
    final = trajectory.segments[-1].end_reason.name if trajectory.segments else "None"
    print(f"case source_index={case['source_index']}")
    print(
        "  conditions: "
        f"radius={case['radius']:.6f}, theta={case['theta']:.9f}, "
        f"speed={case['speed']:.9f}"
    )
    print(
        f"  result: status={status.name}, segments={len(trajectory.segments)}, "
        f"final={final}"
    )
    print(f"  total_ms: {total * 1000.0:.6f}")
    print(f"  mean_ms: {statistics.fmean(timings) * 1000.0:.6f}")
    print(f"  median_ms: {percentile(sorted_timings, 0.50) * 1000.0:.6f}")
    print(f"  p95_ms: {percentile(sorted_timings, 0.95) * 1000.0:.6f}")
    print(f"  min_ms: {min(timings) * 1000.0:.6f}")
    print(f"  max_ms: {max(timings) * 1000.0:.6f}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Benchmark the three worst sampled stress-preview cases."
    )
    parser.add_argument("--repeats", type=int, default=30)
    parser.add_argument("--end-time-days", type=float, default=300.0)
    parser.add_argument("--max-segments", type=int, default=50)
    return parser.parse_args()


def main():
    args = parse_args()
    if args.repeats <= 0:
        raise SystemExit("--repeats must be positive")

    system = make_system()
    builder = make_builder(system)

    all_timings = []
    print("Brahe worst-case stress benchmark")
    print(f"  repeats: {args.repeats}")
    print(f"  end_time_days: {args.end_time_days}")
    print(f"  max_segments: {args.max_segments}")
    print()

    for case in CASES:
        request = make_request(case, args.end_time_days, args.max_segments)
        timings, status, trajectory = run_case(builder, request, args.repeats)
        all_timings.extend(timings)
        print_case_summary(case, timings, status, trajectory)
        print()

    sorted_timings = sorted(all_timings)
    print("overall")
    print(f"  total_runs: {len(all_timings)}")
    print(f"  total_ms: {sum(all_timings) * 1000.0:.6f}")
    print(f"  mean_ms: {statistics.fmean(all_timings) * 1000.0:.6f}")
    print(f"  median_ms: {percentile(sorted_timings, 0.50) * 1000.0:.6f}")
    print(f"  p95_ms: {percentile(sorted_timings, 0.95) * 1000.0:.6f}")
    print(f"  max_ms: {max(all_timings) * 1000.0:.6f}")


if __name__ == "__main__":
    main()
