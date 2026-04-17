import gc
import unittest

import brahe


def make_body_system():
    sun = brahe.BodyDef()
    sun.id = 0
    sun.parent_id = brahe.InvalidBody
    sun.mu = 1000.0
    sun.radius = 1.0
    sun.soi_radius = 1000.0

    moon = brahe.BodyDef()
    moon.id = 1
    moon.parent_id = 0
    moon.mu = 1.0
    moon.radius = 0.2
    moon.soi_radius = 2.0
    moon.orbit_radius = 20.0
    moon.angular_rate = 0.01

    builder = brahe.BodySystemBuilder()
    builder.add_body(sun)
    builder.add_body(moon)
    status, system = builder.build()
    assert status == brahe.SolveStatus.Ok
    return system


class BindingsTests(unittest.TestCase):
    def test_two_body_propagates_state(self):
        state = brahe.State2(brahe.Vec2(10.0, 0.0), brahe.Vec2(0.0, 10.0))

        status, out = brahe.TwoBody.propagate(1000.0, state, 0.1)

        self.assertEqual(status, brahe.SolveStatus.Ok)
        self.assertNotEqual(out.r.y, 0.0)

    def test_body_system_builder_returns_system(self):
        system = make_body_system()

        self.assertEqual(system.root_id(), 0)
        self.assertEqual(system.body_count(), 2)
        self.assertIsNotNone(system.get_body(1))
        self.assertEqual(system.children(0), [1])
        self.assertEqual(system.body_ids(), [0, 1])

    def test_event_detector_returns_status_and_event(self):
        system = make_body_system()
        req = brahe.EventSearchRequest()
        req.central_body = 0
        req.start_time = 0.0
        req.time_limit = 1.0
        req.initial_state = brahe.State2(
            brahe.Vec2(10.0, 0.0), brahe.Vec2(0.0, 10.0)
        )

        detector = brahe.EventDetector(system)
        status, event, diagnostics = detector.find_next_event(req, diagnostics=True)

        self.assertEqual(status, brahe.SolveStatus.Ok)
        self.assertEqual(event.type, brahe.EventType.TimeLimit)
        self.assertGreater(diagnostics.scan_steps, 0)

    def test_trajectory_builder_survives_system_reference_drop(self):
        system = make_body_system()
        builder = brahe.TrajectoryBuilder(system)
        del system
        gc.collect()

        req = brahe.PreviewRequest()
        req.central_body = 0
        req.start_time = 0.0
        req.end_time = 1.0
        req.max_segments = 4
        req.initial_state = brahe.State2(
            brahe.Vec2(10.0, 0.0), brahe.Vec2(0.0, 10.0)
        )

        status, trajectory = builder.build_preview(req)

        self.assertEqual(status, brahe.SolveStatus.Ok)
        self.assertEqual(len(trajectory.segments), 1)

    def test_sampling_returns_plot_friendly_lists(self):
        system = make_body_system()
        builder = brahe.TrajectoryBuilder(system)

        req = brahe.PreviewRequest()
        req.central_body = 0
        req.start_time = 0.0
        req.end_time = 1.0
        req.max_segments = 4
        req.initial_state = brahe.State2(
            brahe.Vec2(10.0, 0.0), brahe.Vec2(0.0, 10.0)
        )
        status, trajectory = builder.build_preview(req)
        self.assertEqual(status, brahe.SolveStatus.Ok)

        samples = brahe.sample_trajectory(trajectory, system, samples_per_segment=5)

        self.assertEqual(len(samples["t"]), 5)
        self.assertEqual(len(samples["x"]), 5)
        self.assertEqual(samples["segment_index"], [0, 0, 0, 0, 0])


class PlottingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import matplotlib

        matplotlib.use("Agg")

    def test_plot_body_system_returns_figure_and_axes(self):
        import matplotlib.pyplot as plt

        system = make_body_system()

        fig, ax = brahe.plot_body_system(system, t=0.0)

        self.assertIs(fig, ax.figure)
        self.assertGreaterEqual(len(ax.lines), 2)
        self.assertGreaterEqual(len(ax.patches), 2)
        plt.close(fig)

    def test_plot_trajectory_returns_figure_and_axes(self):
        import matplotlib.pyplot as plt

        system = make_body_system()
        builder = brahe.TrajectoryBuilder(system)

        req = brahe.PreviewRequest()
        req.central_body = 0
        req.start_time = 0.0
        req.end_time = 1.0
        req.max_segments = 4
        req.initial_state = brahe.State2(
            brahe.Vec2(10.0, 0.0), brahe.Vec2(0.0, 10.0)
        )
        status, trajectory = builder.build_preview(req)
        self.assertEqual(status, brahe.SolveStatus.Ok)

        fig, ax = brahe.plot_trajectory(
            system, trajectory, samples_per_segment=8, show_events=True
        )

        self.assertIs(fig, ax.figure)
        self.assertGreaterEqual(len(ax.lines), 4)
        plt.close(fig)


if __name__ == "__main__":
    unittest.main()
