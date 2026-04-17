import importlib.util
from pathlib import Path

import matplotlib

import brahe


matplotlib.use("Agg")


def load_interactive_example():
    path = Path(__file__).resolve().parents[1] / "examples" / "interactive_plot.py"
    spec = importlib.util.spec_from_file_location("interactive_plot", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_first_logged_failure_plots_without_propagation_error():
    example = load_interactive_example()
    conditions = {
        "end_time_days": 7.0,
        "max_segments": 8,
        "moon_mu": 4902.800066,
        "moon_orbit_radius": 384400.0,
        "samples": 240,
        "vx0": 0.0,
        "vy0": 10.85,
        "x0": 6678.0,
        "y0": 561.9047619047624,
    }

    system = example.make_demo_system(
        moon_mu=conditions["moon_mu"],
        moon_orbit_radius=conditions["moon_orbit_radius"],
    )
    req = example.build_request_from_conditions(conditions)
    status, trajectory = example.build_trajectory(system, req)
    assert status == brahe.SolveStatus.Ok

    fig, ax = brahe.plot_trajectory(
        system,
        trajectory,
        samples_per_segment=conditions["samples"],
        show_bodies=True,
        show_events=True,
    )
    try:
        assert fig is ax.figure
    finally:
        import matplotlib.pyplot as plt

        plt.close(fig)
