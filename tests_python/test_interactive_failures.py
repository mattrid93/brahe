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


def assert_conditions_build_and_plot(conditions):
    example = load_interactive_example()
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


def test_first_logged_failure_plots_without_propagation_error():
    assert_conditions_build_and_plot(
        {
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
    )


def test_second_logged_failure_plots_without_propagation_error():
    assert_conditions_build_and_plot(
        {
            "end_time_days": 7.0,
            "max_segments": 8,
            "moon_mu": 4902.800066,
            "moon_orbit_radius": 384400.0,
            "samples": 240,
            "vx0": 0.0,
            "vy0": 10.85,
            "x0": 6750.793650793651,
            "y0": 0.0,
        }
    )


def test_third_logged_failure_builds_without_numerical_failure():
    assert_conditions_build_and_plot(
        {
            "end_time_days": 7.0,
            "max_segments": 8,
            "moon_mu": 4902.800066,
            "moon_orbit_radius": 384400.0,
            "samples": 240,
            "vx0": 0.0,
            "vy0": 10.85,
            "x0": 6678.0,
            "y0": 726.9841269841272,
        }
    )


def test_default_free_return_impact_is_stable_when_end_time_changes():
    example = load_interactive_example()
    system = example.make_demo_system()
    expected_impact_day = 6.74580583190028

    for end_time_days in (6.75, 6.8, 7.0, 7.5, 8.0, 9.0, 10.0, 12.0, 14.0):
        req = example.build_request(
            example.LOW_EARTH_ORBIT_RADIUS,
            0.0,
            0.0,
            example.TLI_SPEED,
            end_time_days,
            8,
        )
        status, trajectory = example.build_trajectory(system, req)

        assert status == brahe.SolveStatus.Ok
        assert trajectory.segments[-1].end_reason == brahe.EventType.Impact
        assert (
            abs(trajectory.segments[-1].end_time / example.SECONDS_PER_DAY
                - expected_impact_day)
            < 1e-5
        )
