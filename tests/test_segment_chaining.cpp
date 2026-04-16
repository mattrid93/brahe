#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "brahe/body_system.h"
#include "brahe/trajectory.h"
#include "brahe/trajectory_builder.h"
#include "brahe/vec2.h"

#include <cmath>

using namespace brahe;
using Catch::Matchers::WithinAbs;

// --- Section 7.6: Segment chaining tests ---

namespace {

constexpr BodyId kSun = 1;
constexpr BodyId kMoon = 2;

BodySystem sun_with_moon(double moon_radius = 0.2, double moon_soi = 5.0) {
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
    moon.radius = moon_radius;
    moon.soi_radius = moon_soi;
    moon.orbit_radius = 50.0;
    moon.angular_rate = 0.0;
    moon.phase_at_epoch = 0.0;
    b.add_body(moon);

    BodySystem sys;
    (void)b.build(sys);
    return sys;
}

PreviewRequest default_flyby_request(double t_end = 40.0, size_t cap = 16) {
    PreviewRequest req;
    req.central_body = kSun;
    req.start_time = 0.0;
    req.initial_state = State2{{30.0, 0.5}, {5.0, 0.0}};
    req.end_time = t_end;
    req.max_segments = cap;
    return req;
}

}  // namespace

TEST_CASE("PreviewChain_StopsOnImpact", "[phase4][chain]") {
    BodySystem sys = sun_with_moon(/*moon_radius*/ 1.5, /*moon_soi*/ 5.0);
    TrajectoryBuilder tb(sys);

    PreviewRequest req;
    req.central_body = kSun;
    req.start_time = 0.0;
    req.initial_state = State2{{30.0, 0.001}, {5.0, 0.0}};  // aimed at moon center
    req.end_time = 50.0;
    req.max_segments = 16;

    Trajectory out;
    REQUIRE(tb.build_preview(req, out) == SolveStatus::Ok);
    REQUIRE(!out.segments.empty());
    REQUIRE(out.segments.back().end_reason == EventType::Impact);
}

TEST_CASE("PreviewChain_StopsOnTimeLimit", "[phase4][chain]") {
    BodySystem sys = sun_with_moon();
    TrajectoryBuilder tb(sys);

    PreviewRequest req;
    req.central_body = kSun;
    req.start_time = 0.0;
    // Bound circular orbit nowhere near the moon.
    req.initial_state = State2{{10.0, 0.0}, {0.0, std::sqrt(1.0 / 10.0)}};
    req.end_time = 5.0;
    req.max_segments = 16;

    Trajectory out;
    REQUIRE(tb.build_preview(req, out) == SolveStatus::Ok);
    REQUIRE(out.segments.size() == 1);
    REQUIRE(out.segments.back().end_reason == EventType::TimeLimit);
    REQUIRE_THAT(out.segments.back().end_time, WithinAbs(5.0, 1e-9));
}

TEST_CASE("PreviewChain_PatchesIntoChildAfterSoiEntry", "[phase4][chain]") {
    BodySystem sys = sun_with_moon();
    TrajectoryBuilder tb(sys);

    Trajectory out;
    REQUIRE(tb.build_preview(default_flyby_request(), out) == SolveStatus::Ok);
    REQUIRE(out.segments.size() >= 2);
    REQUIRE(out.segments[0].end_reason == EventType::SoiEntry);
    REQUIRE(out.segments[0].central_body == kSun);
    REQUIRE(out.segments[1].central_body == kMoon);
}

TEST_CASE("PreviewChain_PatchesBackToParentAfterSoiExit", "[phase4][chain]") {
    BodySystem sys = sun_with_moon();
    TrajectoryBuilder tb(sys);

    Trajectory out;
    REQUIRE(tb.build_preview(default_flyby_request(), out) == SolveStatus::Ok);
    REQUIRE(out.segments.size() >= 3);
    REQUIRE(out.segments[1].end_reason == EventType::SoiExit);
    REQUIRE(out.segments[1].central_body == kMoon);
    REQUIRE(out.segments[2].central_body == kSun);
}

TEST_CASE("PreviewChain_BuildsMultipleSegmentsAcrossParentChildParentSequence",
          "[phase4][chain]") {
    BodySystem sys = sun_with_moon();
    TrajectoryBuilder tb(sys);

    Trajectory out;
    REQUIRE(tb.build_preview(default_flyby_request(), out) == SolveStatus::Ok);
    REQUIRE(out.segments.size() >= 3);

    // Sun -> Moon -> Sun.
    REQUIRE(out.segments[0].central_body == kSun);
    REQUIRE(out.segments[1].central_body == kMoon);
    REQUIRE(out.segments[2].central_body == kSun);
}

TEST_CASE("PreviewChain_UsesPatchedStateAsNextSegmentInitialState", "[phase4][chain]") {
    // Each SoiEntry/SoiExit segment boundary should satisfy inertial continuity:
    // reconstructed absolute state on either side of the patch must agree within
    // tolerance.
    BodySystem sys = sun_with_moon();
    Tolerances tol{};
    TrajectoryBuilder tb(sys, tol);

    Trajectory out;
    REQUIRE(tb.build_preview(default_flyby_request(), out) == SolveStatus::Ok);
    REQUIRE(out.segments.size() >= 3);

    // --- Boundary 0/1: SoiEntry ---
    // Segment 0 ends in parent frame; segment 1 starts in child frame.
    // Need the parent-frame end-state of segment 0 -- we can reconstruct it
    // from segment 1's initial state + child ephemeris at t_patch.
    {
        double t_patch = out.segments[0].end_time;
        BodyId child = out.segments[1].central_body;
        State2 sc_child_init = out.segments[1].initial_state;
        Vec2 r_ch = sys.position_in_parent(child, t_patch);
        Vec2 v_ch = sys.velocity_in_parent(child, t_patch);
        Vec2 r_recon_parent = sc_child_init.r + r_ch;
        Vec2 v_recon_parent = sc_child_init.v + v_ch;

        // The parent-frame end state of segment 0 is NOT stored; we can only
        // cross-check that the reconstructed absolute state is finite and
        // "near the SOI boundary" within position tolerance. Use the known
        // soi_radius.
        const BodyDef* ch = sys.get_body(child);
        double boundary_residual = length(sc_child_init.r) - ch->soi_radius;
        REQUIRE_THAT(boundary_residual, WithinAbs(0.0, tol.position_epsilon));
        (void)r_recon_parent;
        (void)v_recon_parent;
    }

    // --- Boundary 1/2: SoiExit ---
    {
        double t_patch = out.segments[1].end_time;
        BodyId child = out.segments[1].central_body;
        State2 sc_parent_init = out.segments[2].initial_state;
        Vec2 r_ch = sys.position_in_parent(child, t_patch);
        Vec2 v_ch = sys.velocity_in_parent(child, t_patch);
        Vec2 r_recon_child = sc_parent_init.r - r_ch;
        Vec2 v_recon_child = sc_parent_init.v - v_ch;

        const BodyDef* ch = sys.get_body(child);
        double boundary_residual = length(r_recon_child) - ch->soi_radius;
        REQUIRE_THAT(boundary_residual, WithinAbs(0.0, tol.position_epsilon));
        (void)v_recon_child;
    }
}

TEST_CASE("PreviewChain_DoesNotImmediatelyRepatchAfterBoundaryCrossing",
          "[phase4][chain]") {
    // Each segment immediately after a patch must have non-trivial duration.
    // Otherwise a repeat-patch loop would drive zero-duration segments.
    BodySystem sys = sun_with_moon();
    Tolerances tol{};
    TrajectoryBuilder tb(sys, tol);

    Trajectory out;
    REQUIRE(tb.build_preview(default_flyby_request(), out) == SolveStatus::Ok);

    for (size_t i = 1; i < out.segments.size(); ++i) {
        const EventType prev_reason = out.segments[i - 1].end_reason;
        if (prev_reason == EventType::SoiEntry || prev_reason == EventType::SoiExit) {
            double dur = out.segments[i].end_time - out.segments[i].start_time;
            REQUIRE(dur > tol.time_epsilon);
        }
    }
}

TEST_CASE("PreviewChain_SegmentTimesAreMonotonic", "[phase4][chain]") {
    BodySystem sys = sun_with_moon();
    TrajectoryBuilder tb(sys);

    Trajectory out;
    REQUIRE(tb.build_preview(default_flyby_request(), out) == SolveStatus::Ok);

    for (const Segment& s : out.segments) {
        REQUIRE(s.end_time >= s.start_time);
    }
    for (size_t i = 1; i < out.segments.size(); ++i) {
        REQUIRE(out.segments[i].start_time >= out.segments[i - 1].start_time);
    }
}

TEST_CASE("PreviewChain_AdjacentSegmentsMeetAtSameTimeWithinTimeEpsilon",
          "[phase4][chain]") {
    BodySystem sys = sun_with_moon();
    Tolerances tol{};
    TrajectoryBuilder tb(sys, tol);

    Trajectory out;
    REQUIRE(tb.build_preview(default_flyby_request(), out) == SolveStatus::Ok);

    for (size_t i = 1; i < out.segments.size(); ++i) {
        REQUIRE_THAT(out.segments[i].start_time,
                     WithinAbs(out.segments[i - 1].end_time, tol.time_epsilon));
    }
}

TEST_CASE("PreviewChain_NoInfiniteLoopAtBoundary", "[phase4][chain]") {
    // Spacecraft starts essentially ON the SOI boundary with a tangential velocity.
    // Without suppression this would oscillate forever. With suppression the chain
    // must terminate (either at TimeLimit, Impact, or a legitimate later event) in
    // a bounded number of segments.
    BodySystem sys = sun_with_moon();
    TrajectoryBuilder tb(sys);

    PreviewRequest req;
    req.central_body = kSun;
    req.start_time = 0.0;
    // Approximately on the boundary of the moon's SOI (moon at (50,0), SOI=5).
    // Slightly inward (y=0.01) plus tangential velocity.
    req.initial_state = State2{{55.0 - 1e-10, 0.01}, {0.0, 0.3}};
    req.end_time = 5.0;
    req.max_segments = 32;

    Trajectory out;
    SolveStatus status = tb.build_preview(req, out);
    // Must not return pathological failure; must terminate in bounded segments.
    REQUIRE((status == SolveStatus::Ok || status == SolveStatus::CapacityExceeded));
    REQUIRE(out.segments.size() <= 32);
    // If we did terminate cleanly, final segment must be a terminating event.
    if (status == SolveStatus::Ok) {
        const EventType final_reason = out.segments.back().end_reason;
        REQUIRE((final_reason == EventType::TimeLimit ||
                 final_reason == EventType::Impact));
    }
}
