#include <catch2/catch_test_macros.hpp>

#include "brahe/body_system.h"
#include "brahe/trajectory.h"
#include "brahe/trajectory_builder.h"
#include "brahe/vec2.h"

#include <cmath>
#include <random>

using namespace brahe;

// --- Section 7.10: Fuzz tests ---

namespace {

struct SystemBuild {
    BodySystem sys;
    bool ok = false;
};

// Deterministically-seeded random two- or three-body system. Returns ok=false
// if the random draw produced constraints the builder rejects.
SystemBuild random_system(std::mt19937_64& rng) {
    std::uniform_int_distribution<int> n_children(1, 2);
    std::uniform_real_distribution<double> orbit(30.0, 80.0);
    std::uniform_real_distribution<double> soi(2.0, 8.0);
    std::uniform_real_distribution<double> rad(0.1, 1.0);
    std::uniform_real_distribution<double> rate(-0.1, 0.1);
    std::uniform_real_distribution<double> phase(-M_PI, M_PI);

    BodySystemBuilder b;
    BodyDef sun{};
    sun.id = 1;
    sun.parent_id = InvalidBody;
    sun.mu = 1.0;
    sun.radius = 0.1;
    sun.soi_radius = 1e6;
    b.add_body(sun);

    int nc = n_children(rng);
    for (int i = 0; i < nc; ++i) {
        BodyDef m{};
        m.id = static_cast<BodyId>(2 + i);
        m.parent_id = 1;
        m.mu = 0.001;
        m.radius = rad(rng);
        double orb = orbit(rng);
        double s = soi(rng);
        // Constraints: soi > radius, soi < orbit_radius.
        if (s <= m.radius) s = m.radius + 0.5;
        if (s >= orb) s = orb * 0.5;
        m.soi_radius = s;
        m.orbit_radius = orb;
        m.angular_rate = rate(rng);
        m.phase_at_epoch = phase(rng);
        b.add_body(m);
    }

    SystemBuild out;
    out.ok = (b.build(out.sys) == SolveStatus::Ok);
    return out;
}

bool is_finite_state(const State2& s) {
    return std::isfinite(s.r.x) && std::isfinite(s.r.y) &&
           std::isfinite(s.v.x) && std::isfinite(s.v.y);
}

bool is_finite_segment(const Segment& s) {
    return std::isfinite(s.start_time) && std::isfinite(s.end_time) &&
           is_finite_state(s.initial_state);
}

}  // namespace

TEST_CASE("RandomValidSystemsAndStates_NoCrashesNoNaNs", "[phase4][fuzz]") {
    std::mt19937_64 rng(0xF00DULL);
    std::uniform_real_distribution<double> r_dist(-80.0, 80.0);
    std::uniform_real_distribution<double> v_dist(-2.0, 2.0);

    for (int trial = 0; trial < 60; ++trial) {
        SystemBuild sb = random_system(rng);
        if (!sb.ok) continue;

        PreviewRequest req;
        req.central_body = 1;
        req.start_time = 0.0;
        req.initial_state = State2{{r_dist(rng), r_dist(rng)}, {v_dist(rng), v_dist(rng)}};
        if (length_squared(req.initial_state.r) < 1.0) continue;
        req.end_time = 20.0;
        req.max_segments = 16;

        Trajectory out;
        TrajectoryBuilder tb(sb.sys);
        (void)tb.build_preview(req, out);  // must not crash
        for (const Segment& s : out.segments) {
            REQUIRE(is_finite_segment(s));
        }
    }
}

TEST_CASE("RandomNestedParentChildTransitions_NoImmediateRepatchLoops",
          "[phase4][fuzz]") {
    // Any patch-ending segment must be followed by a segment with non-trivial
    // duration. The suppression rule forbids degenerate repatch oscillations.
    std::mt19937_64 rng(0xBADA55ULL);
    std::uniform_real_distribution<double> r_dist(20.0, 70.0);
    std::uniform_real_distribution<double> v_dist(-3.0, 3.0);

    Tolerances tol{};
    for (int trial = 0; trial < 60; ++trial) {
        SystemBuild sb = random_system(rng);
        if (!sb.ok) continue;

        PreviewRequest req;
        req.central_body = 1;
        req.start_time = 0.0;
        req.initial_state =
            State2{{r_dist(rng), r_dist(rng) - 45.0}, {v_dist(rng), v_dist(rng)}};
        if (length_squared(req.initial_state.r) < 1.0) continue;
        req.end_time = 30.0;
        req.max_segments = 32;

        Trajectory out;
        TrajectoryBuilder tb(sb.sys, tol);
        (void)tb.build_preview(req, out);

        for (size_t i = 1; i < out.segments.size(); ++i) {
            EventType prev = out.segments[i - 1].end_reason;
            if (prev == EventType::SoiEntry || prev == EventType::SoiExit) {
                double dur = out.segments[i].end_time - out.segments[i].start_time;
                REQUIRE(dur > tol.time_epsilon);
            }
        }
    }
}

TEST_CASE("RandomPreviewBuilds_AllSegmentBodiesAreValid", "[phase4][fuzz]") {
    std::mt19937_64 rng(0xB00B1EULL);
    std::uniform_real_distribution<double> r_dist(-70.0, 70.0);
    std::uniform_real_distribution<double> v_dist(-2.0, 2.0);

    for (int trial = 0; trial < 60; ++trial) {
        SystemBuild sb = random_system(rng);
        if (!sb.ok) continue;

        PreviewRequest req;
        req.central_body = 1;
        req.start_time = 0.0;
        req.initial_state = State2{{r_dist(rng), r_dist(rng)}, {v_dist(rng), v_dist(rng)}};
        if (length_squared(req.initial_state.r) < 1.0) continue;
        req.end_time = 20.0;
        req.max_segments = 16;

        Trajectory out;
        TrajectoryBuilder tb(sb.sys);
        (void)tb.build_preview(req, out);

        for (const Segment& s : out.segments) {
            REQUIRE(sb.sys.get_body(s.central_body) != nullptr);
        }
    }
}

TEST_CASE("RandomPreviewBuilds_AllSegmentTimesAreOrdered", "[phase4][fuzz]") {
    std::mt19937_64 rng(0x1234ABCDULL);
    std::uniform_real_distribution<double> r_dist(-70.0, 70.0);
    std::uniform_real_distribution<double> v_dist(-2.0, 2.0);

    Tolerances tol{};
    for (int trial = 0; trial < 60; ++trial) {
        SystemBuild sb = random_system(rng);
        if (!sb.ok) continue;

        PreviewRequest req;
        req.central_body = 1;
        req.start_time = 0.0;
        req.initial_state = State2{{r_dist(rng), r_dist(rng)}, {v_dist(rng), v_dist(rng)}};
        if (length_squared(req.initial_state.r) < 1.0) continue;
        req.end_time = 20.0;
        req.max_segments = 16;

        Trajectory out;
        TrajectoryBuilder tb(sb.sys, tol);
        (void)tb.build_preview(req, out);

        for (const Segment& s : out.segments) {
            REQUIRE(s.end_time >= s.start_time);
        }
        for (size_t i = 1; i < out.segments.size(); ++i) {
            REQUIRE(out.segments[i].start_time >= out.segments[i - 1].start_time - tol.time_epsilon);
            // Adjacent segment boundaries meet.
            REQUIRE(std::abs(out.segments[i].start_time -
                             out.segments[i - 1].end_time) <= tol.time_epsilon);
        }
    }
}
