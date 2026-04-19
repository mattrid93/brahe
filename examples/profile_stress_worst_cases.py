#!/usr/bin/env python3
from __future__ import annotations

import argparse
import cProfile
import pstats
import runpy
import sys
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Profile the deterministic worst-case stress benchmark."
    )
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--sort", default="cumulative")
    parser.add_argument("--limit", type=int, default=30)
    return parser.parse_args()


def main():
    args = parse_args()
    benchmark = Path(__file__).resolve().parent / "stress_worst_cases_benchmark.py"
    profiler = cProfile.Profile()
    old_argv = sys.argv
    try:
        sys.argv = [str(benchmark), "--repeats", str(args.repeats)]
        profiler.enable()
        runpy.run_path(str(benchmark), run_name="__main__")
        profiler.disable()
    finally:
        sys.argv = old_argv

    stats = pstats.Stats(profiler)
    stats.strip_dirs().sort_stats(args.sort).print_stats(args.limit)


if __name__ == "__main__":
    main()
