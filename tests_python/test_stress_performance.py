import importlib.util
import math
import statistics
import time
from collections import Counter
from pathlib import Path

import brahe


SAMPLES = 100
SEED = 1
END_TIME_DAYS = 300.0
MAX_SEGMENTS = 50

# Local performance budget for this development machine. Representative runs
# before adding this test were: total ~= 0.32s, mean ~= 3.3ms, p95 ~= 15ms,
# max ~= 22ms.
MAX_TOTAL_SECONDS = 0.45
MAX_MEAN_MS = 4.5
MAX_P95_MS = 25.0
MAX_SINGLE_PREVIEW_MS = 35.0


def load_stress_benchmark():
    path = Path(__file__).resolve().parents[1] / "examples" / "stress_benchmark.py"
    spec = importlib.util.spec_from_file_location("stress_benchmark", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def percentile(sorted_values, fraction):
    position = fraction * (len(sorted_values) - 1)
    lo = int(math.floor(position))
    hi = int(math.ceil(position))
    if lo == hi:
        return sorted_values[lo]
    weight = position - lo
    return sorted_values[lo] * (1.0 - weight) + sorted_values[hi] * weight


def test_stress_preview_build_performance_local_machine():
    benchmark = load_stress_benchmark()
    viewer = benchmark.load_stress_viewer()
    system = viewer.make_stress_system()

    tolerances = brahe.Tolerances()
    tolerances.root_epsilon = 1e-3
    tolerances.max_event_refine_iterations = 96
    builder = brahe.TrajectoryBuilder(system, tolerances)

    requests, _ = benchmark.make_requests(
        viewer,
        SAMPLES,
        SEED,
        END_TIME_DAYS,
        MAX_SEGMENTS,
    )

    timings = []
    status_counts = Counter()
    segment_counts = Counter()

    for req in requests:
        start = time.perf_counter()
        status, trajectory = builder.build_preview(req)
        timings.append(time.perf_counter() - start)
        status_counts[status.name] += 1
        segment_counts[len(trajectory.segments)] += 1

    sorted_timings = sorted(timings)
    total = sum(timings)
    mean_ms = statistics.fmean(timings) * 1000.0
    p95_ms = percentile(sorted_timings, 0.95) * 1000.0
    max_ms = max(timings) * 1000.0
    summary = (
        f"total={total:.6f}s mean={mean_ms:.3f}ms p95={p95_ms:.3f}ms "
        f"max={max_ms:.3f}ms statuses={dict(status_counts)} "
        f"segments={dict(sorted(segment_counts.items()))}"
    )

    assert total < MAX_TOTAL_SECONDS, summary
    assert mean_ms < MAX_MEAN_MS, summary
    assert p95_ms < MAX_P95_MS, summary
    assert max_ms < MAX_SINGLE_PREVIEW_MS, summary
