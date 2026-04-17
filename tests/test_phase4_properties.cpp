#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "brahe/body_system.h"
#include "brahe/patcher.h"
#include "brahe/trajectory.h"
#include "brahe/trajectory_builder.h"
#include "brahe/two_body.h"
#include "brahe/vec2.h"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <random>

using namespace brahe;
using Catch::Matchers::WithinAbs;

// --- Section 7.9: Property tests ---

namespace {

constexpr BodyId kSun = 1;
constexpr BodyId kMoon = 2;

BodySystem sun_with_moon(double phase = 0.0, double rate = 0.08) {
    BodySystemBuilder b;
    BodyDef sun{};
    sun.id = kSun;
    sun.parent_id = InvalidBody;
    sun.mu = 1.0;
    sun.radius = 0.1;
    sun.soi_radius = 1e6;
    b.add_body(sun);

    BodyDef moon{};
    moon.id = kMoon;
    moon.parent_id = kSun;
    moon.mu = 0.001;
    moon.radius = 0.2;
    moon.soi_radius = 5.0;
    moon.orbit_radius = 50.0;
    moon.angular_rate = rate;
    moon.phase_at_epoch = phase;
    b.add_body(moon);

    BodySystem sys;
    (void)b.build(sys);
    return sys;
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

TEST_CASE("PatchRoundTrip_PreservesAbsoluteInertialState", "[phase4][property][patch]") {
    BodySystem sys = sun_with_moon(0.3, 0.05);
    Patcher patcher(sys);
    Tolerances tol{};

    std::mt19937_64 rng(0xC0FFEEULL);
    std::uniform_real_distribution<double> r_dist(-60.0, 60.0);
    std::uniform_real_distribution<double> v_dist(-2.0, 2.0);
    std::uniform_real_distribution<double> t_dist(0.0, 50.0);

    for (int trial = 0; trial < 50; ++trial) {
        double t = t_dist(rng);
        State2 sc_parent{{r_dist(rng), r_dist(rng)}, {v_dist(rng), v_dist(rng)}};

        State2 in_child{};
        REQUIRE(patcher.patch_parent_to_child(t, kSun, kMoon, sc_parent, in_child) ==
                SolveStatus::Ok);

        State2 back_to_parent{};
        REQUIRE(patcher.patch_child_to_parent(t, kMoon, in_child, back_to_parent) ==
                SolveStatus::Ok);

        REQUIRE_THAT(back_to_parent.r.x, WithinAbs(sc_parent.r.x, tol.position_epsilon));
        REQUIRE_THAT(back_to_parent.r.y, WithinAbs(sc_parent.r.y, tol.position_epsilon));
        REQUIRE_THAT(back_to_parent.v.x, WithinAbs(sc_parent.v.x, tol.velocity_epsilon));
        REQUIRE_THAT(back_to_parent.v.y, WithinAbs(sc_parent.v.y, tol.velocity_epsilon));
    }
}

TEST_CASE("PreviewChain_EachNextSegmentBeginsFromPreviousTerminalEventState",
          "[phase4][property][chain]") {
    // For SoiEntry boundaries: segment[i+1].initial_state should be the
    // parent-frame propagated end of segment[i] minus child ephemeris at t_patch.
    // We approximate the parent-frame end by propagating segment[i] to its end_time.
    BodySystem sys = sun_with_moon();
    Tolerances tol{};
    TrajectoryBuilder tb(sys, tol);

    PreviewRequest req;
    req.central_body = kSun;
    req.start_time = 0.0;
    req.initial_state = State2{{30.0, 0.5}, {5.0, 0.0}};
    req.end_time = 40.0;
    req.max_segments = 16;

    Trajectory out;
    REQUIRE(tb.build_preview(req, out) == SolveStatus::Ok);

    for (size_t i = 0; i + 1 < out.segments.size(); ++i) {
        const Segment& cur = out.segments[i];
        const Segment& next = out.segments[i + 1];

        const BodyDef* cb = sys.get_body(cur.central_body);
        REQUIRE(cb != nullptr);

        // Propagate cur to end_time.
        State2 end_state{};
        double dt = cur.end_time - cur.start_time;
        SolveStatus ps = TwoBody::propagate(cb->mu, cur.initial_state, dt, end_state);
        REQUIRE(ps == SolveStatus::Ok);
        REQUIRE(is_finite_state(end_state));

        if (cur.end_reason == EventType::SoiEntry) {
            // Reconstruct the parent-frame state from the next segment's
            // initial state (which is in the child frame).
            Vec2 r_ch = sys.position_in_parent(next.central_body, cur.end_time);
            Vec2 v_ch = sys.velocity_in_parent(next.central_body, cur.end_time);
            Vec2 r_recon = next.initial_state.r + r_ch;
            Vec2 v_recon = next.initial_state.v + v_ch;

            REQUIRE_THAT(r_recon.x, WithinAbs(end_state.r.x, tol.position_epsilon));
            REQUIRE_THAT(r_recon.y, WithinAbs(end_state.r.y, tol.position_epsilon));
            REQUIRE_THAT(v_recon.x, WithinAbs(end_state.v.x, tol.velocity_epsilon));
            REQUIRE_THAT(v_recon.y, WithinAbs(end_state.v.y, tol.velocity_epsilon));
        } else if (cur.end_reason == EventType::SoiExit) {
            // Reconstruct child-frame end from next parent-frame initial.
            Vec2 r_ch = sys.position_in_parent(cur.central_body, cur.end_time);
            Vec2 v_ch = sys.velocity_in_parent(cur.central_body, cur.end_time);
            Vec2 r_recon_child = next.initial_state.r - r_ch;
            Vec2 v_recon_child = next.initial_state.v - v_ch;

            REQUIRE_THAT(r_recon_child.x, WithinAbs(end_state.r.x, tol.position_epsilon));
            REQUIRE_THAT(r_recon_child.y, WithinAbs(end_state.r.y, tol.position_epsilon));
            REQUIRE_THAT(v_recon_child.x, WithinAbs(end_state.v.x, tol.velocity_epsilon));
            REQUIRE_THAT(v_recon_child.y, WithinAbs(end_state.v.y, tol.velocity_epsilon));
        }
    }
}

TEST_CASE("PreviewChain_FinalSegmentAlwaysEndsAtImpactOrTimeLimitOrPatchEvent",
          "[phase4][property][chain]") {
    BodySystem sys = sun_with_moon();
    TrajectoryBuilder tb(sys);

    // A few varied starts.
    struct Case {
        State2 s0;
        double t_end;
    };
    Case cases[] = {
        {State2{{30.0, 0.5}, {5.0, 0.0}}, 40.0},
        {State2{{10.0, 0.0}, {0.0, 0.31622776601683794}}, 5.0},  // bound circular
        {State2{{30.0, 0.001}, {5.0, 0.0}}, 50.0},
    };

    for (const Case& c : cases) {
        PreviewRequest req;
        req.central_body = kSun;
        req.start_time = 0.0;
        req.initial_state = c.s0;
        req.end_time = c.t_end;
        req.max_segments = 16;

        Trajectory out;
        SolveStatus status = tb.build_preview(req, out);
        REQUIRE(!out.segments.empty());

        const EventType final_reason = out.segments.back().end_reason;
        if (status == SolveStatus::Ok) {
            REQUIRE((final_reason == EventType::Impact ||
                     final_reason == EventType::TimeLimit ||
                     final_reason == EventType::SoiExit));
        } else {
            // CapacityExceeded -> final reason is a patch event.
            REQUIRE(status == SolveStatus::CapacityExceeded);
            REQUIRE((final_reason == EventType::SoiEntry ||
                     final_reason == EventType::SoiExit));
        }
    }
}

TEST_CASE("PreviewChain_NoNaNsInAnySegmentEndpoints", "[phase4][property][chain]") {
    BodySystem sys = sun_with_moon();
    TrajectoryBuilder tb(sys);

    std::mt19937_64 rng(0xBEEFULL);
    std::uniform_real_distribution<double> r_dist(-80.0, 80.0);
    std::uniform_real_distribution<double> v_dist(-1.5, 1.5);

    for (int trial = 0; trial < 40; ++trial) {
        PreviewRequest req;
        req.central_body = kSun;
        req.start_time = 0.0;
        req.initial_state = State2{{r_dist(rng), r_dist(rng)}, {v_dist(rng), v_dist(rng)}};
        // Reject initial r ~ 0 (event detector rejects those).
        if (length_squared(req.initial_state.r) < 1.0) continue;
        req.end_time = 30.0;
        req.max_segments = 16;

        Trajectory out;
        (void)tb.build_preview(req, out);
        for (const Segment& s : out.segments) {
            REQUIRE(is_finite_segment(s));
        }
    }
}

TEST_CASE("DynamicAndFixedPreviewAgreeUpToCapacity", "[phase4][property][chain]") {
    BodySystem sys = sun_with_moon(0.0, 0.0);
    TrajectoryBuilder tb(sys);

    PreviewRequest req;
    req.central_body = kSun;
    req.start_time = 0.0;
    req.initial_state = State2{{30.0, 0.5}, {5.0, 0.0}};
    req.end_time = 40.0;
    req.max_segments = 16;

    Trajectory dyn;
    REQUIRE(tb.build_preview(req, dyn) == SolveStatus::Ok);

    // Tight fixed-cap that still fits the full chain.
    TrajectoryFixed<16> fxd;
    REQUIRE(tb.build_preview_fixed(req, fxd) == SolveStatus::Ok);

    REQUIRE(dyn.segments.size() == fxd.count);
    for (size_t i = 0; i < dyn.segments.size(); ++i) {
        REQUIRE(std::memcmp(&dyn.segments[i], &fxd.segments[i], sizeof(Segment)) == 0);
    }

    // And a strict cap that truncates early.
    PreviewRequest req_small = req;
    req_small.max_segments = 1;
    TrajectoryFixed<1> small;
    SolveStatus status = tb.build_preview_fixed(req_small, small);
    REQUIRE(status == SolveStatus::CapacityExceeded);
    REQUIRE(small.count == 1);
    REQUIRE(std::memcmp(&small.segments[0], &dyn.segments[0], sizeof(Segment)) == 0);
}
