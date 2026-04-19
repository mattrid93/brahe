#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib.util
import math
import random
import statistics
import sys
import time
from collections import Counter
from pathlib import Path

import brahe


def load_stress_viewer():
    path = Path(__file__).resolve().parent / "stress_viewer.py"
    spec = importlib.util.spec_from_file_location("stress_viewer", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def percentile(sorted_values, fraction):
    if not sorted_values:
        return 0.0
    if len(sorted_values) == 1:
        return sorted_values[0]
    position = fraction * (len(sorted_values) - 1)
    lo = int(math.floor(position))
    hi = int(math.ceil(position))
    if lo == hi:
        return sorted_values[lo]
    weight = position - lo
    return sorted_values[lo] * (1.0 - weight) + sorted_values[hi] * weight


def make_requests(viewer, count, seed, end_time_days, max_segments):
    rng = random.Random(seed)
    radius_min = viewer.PLANET_RADIUS + 1_000.0
    radius_max = 1_600_000.0
    requests = []
    conditions = []

    for _ in range(count):
        radius = rng.uniform(radius_min, radius_max)
        theta = rng.uniform(0.0, 2.0 * math.pi)
        circular_speed = viewer.local_circular_speed(radius)
        speed = rng.uniform(0.0, 1.5 * circular_speed)

        req = viewer.build_request(
            radius,
            theta,
            speed,
            end_time_days,
            max_segments,
        )
        requests.append(req)
        conditions.append((radius, theta, speed, circular_speed))

    return requests, conditions


def run_benchmark(args):
    viewer = load_stress_viewer()
    system = viewer.make_stress_system()

    tolerances = brahe.Tolerances()
    tolerances.root_epsilon = 1e-3
    tolerances.max_event_refine_iterations = 96
    builder = brahe.TrajectoryBuilder(system, tolerances)

    requests, conditions = make_requests(
        viewer,
        args.samples,
        args.seed,
        args.end_time_days,
        args.max_segments,
    )

    timings = []
    status_counts = Counter()
    segment_counts = Counter()
    event_counts = Counter()
    worst = []

    wall_start = time.perf_counter()
    for i, req in enumerate(requests):
        start = time.perf_counter()
        status, trajectory = builder.build_preview(req)
        elapsed = time.perf_counter() - start

        timings.append(elapsed)
        status_counts[status.name] += 1
        segment_counts[len(trajectory.segments)] += 1
        if trajectory.segments:
            event_counts[trajectory.segments[-1].end_reason.name] += 1

        if len(worst) < args.worst:
            worst.append((elapsed, i, status, trajectory))
            worst.sort(reverse=True, key=lambda item: item[0])
        elif args.worst > 0 and elapsed > worst[-1][0]:
            worst[-1] = (elapsed, i, status, trajectory)
            worst.sort(reverse=True, key=lambda item: item[0])

        if args.progress and (i + 1) % args.progress == 0:
            print(f"completed {i + 1}/{args.samples}", file=sys.stderr)

    wall_elapsed = time.perf_counter() - wall_start
    sorted_timings = sorted(timings)
    total_build_time = sum(timings)

    print("Brahe stress preview benchmark")
    print(f"  samples: {args.samples}")
    print(f"  seed: {args.seed}")
    print(f"  end_time_days: {args.end_time_days}")
    print(f"  max_segments: {args.max_segments}")
    print("  radius_range_km: "
          f"[{viewer.PLANET_RADIUS + 1_000.0:.1f}, 1600000.0]")
    print("  speed_range: [0, 1.5 * local circular speed]")
    print()
    print("Timing")
    print(f"  wall_time_s: {wall_elapsed:.6f}")
    print(f"  total_build_time_s: {total_build_time:.6f}")
    print(f"  previews_per_second: {args.samples / total_build_time:.3f}")
    print(f"  mean_ms: {statistics.fmean(timings) * 1000.0:.6f}")
    print(f"  median_ms: {percentile(sorted_timings, 0.50) * 1000.0:.6f}")
    print(f"  p90_ms: {percentile(sorted_timings, 0.90) * 1000.0:.6f}")
    print(f"  p95_ms: {percentile(sorted_timings, 0.95) * 1000.0:.6f}")
    print(f"  p99_ms: {percentile(sorted_timings, 0.99) * 1000.0:.6f}")
    print(f"  max_ms: {max(timings) * 1000.0:.6f}")
    print()
    print("Statuses")
    for name, count in status_counts.most_common():
        print(f"  {name}: {count}")
    print()
    print("Final events")
    for name, count in event_counts.most_common():
        print(f"  {name}: {count}")
    print()
    print("Segment counts")
    for segment_count, count in sorted(segment_counts.items()):
        print(f"  {segment_count}: {count}")

    if args.worst > 0:
        print()
        print(f"Worst {len(worst)} cases")
        for elapsed, index, status, trajectory in worst:
            radius, theta, speed, circular_speed = conditions[index]
            final = trajectory.segments[-1].end_reason.name if trajectory.segments else "None"
            print(
                f"  #{index}: {elapsed * 1000.0:.6f} ms, "
                f"status={status.name}, segments={len(trajectory.segments)}, "
                f"final={final}, radius={radius:.6f}, theta={theta:.9f}, "
                f"speed={speed:.9f}, speed/vcirc={speed / circular_speed:.6f}"
            )


def parse_args():
    parser = argparse.ArgumentParser(
        description="Benchmark stress-system trajectory preview construction."
    )
    parser.add_argument("--samples", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--end-time-days", type=float, default=300.0)
    parser.add_argument("--max-segments", type=int, default=24)
    parser.add_argument("--worst", type=int, default=5)
    parser.add_argument(
        "--progress",
        type=int,
        default=0,
        help="print progress every N previews to stderr; 0 disables progress",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    if args.samples <= 0:
        raise SystemExit("--samples must be positive")
    if args.max_segments <= 0:
        raise SystemExit("--max-segments must be positive")
    run_benchmark(args)


if __name__ == "__main__":
    main()
