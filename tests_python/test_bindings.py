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
            brahe.Vec2(0.0, 19.0), brahe.Vec2(-7.25, 0.0)
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

    def test_uniform_sampling_uses_even_times_across_segments(self):
        system = make_body_system()
        trajectory = brahe.Trajectory()

        first = brahe.Segment()
        first.central_body = 0
        first.start_time = 0.0
        first.end_time = 2.0
        first.end_reason = brahe.EventType.SoiEntry
        first.initial_state = brahe.State2(
            brahe.Vec2(10.0, 0.0), brahe.Vec2(0.0, 10.0)
        )

        second = brahe.Segment()
        second.central_body = 0
        second.start_time = 2.0
        second.end_time = 10.0
        second.end_reason = brahe.EventType.TimeLimit
        second.initial_state = brahe.State2(
            brahe.Vec2(10.0, 0.0), brahe.Vec2(0.0, 10.0)
        )

        trajectory.segments = [first, second]

        samples = brahe.sample_trajectory_uniform(trajectory, system, samples=6)

        self.assertEqual(samples["t"], [0.0, 2.0, 4.0, 6.0, 8.0, 10.0])
        self.assertEqual(samples["segment_index"], [0, 0, 1, 1, 1, 1])


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

    def test_body_artist_cleanup_removes_only_artists(self):
        import matplotlib.pyplot as plt
        from matplotlib.patches import Circle

        fig, ax = plt.subplots()
        marker, = ax.plot([0.0], [0.0])
        soi_circle = Circle((0.0, 0.0), 1.0)
        ax.add_patch(soi_circle)
        label = ax.annotate("1", (0.0, 0.0))
        body_artists = [(1, marker, soi_circle, label)]

        for _, marker, soi_circle, label in body_artists:
            for artist in (marker, soi_circle, label):
                artist.remove()

        self.assertEqual(len(ax.lines), 0)
        self.assertEqual(len(ax.patches), 0)
        self.assertEqual(len(ax.texts), 0)
        plt.close(fig)


if __name__ == "__main__":
    unittest.main()
