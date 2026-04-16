#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "brahe/body_system.h"
#include "brahe/event_detector.h"
#include "brahe/vec2.h"

#include <cmath>

using namespace brahe;
using Catch::Matchers::WithinAbs;

// --- Radial (h = 0) event detection tests ---
//
// The EventDetector degrades gracefully when TwoBody::propagate fails, so
// radial trajectories currently get reported as TimeLimit regardless of what
// they would physically do. Once the propagator supports the radial case,
// these tests pin down the physically correct event classifications.

namespace {

BodySystem build_root_only(double mu, double radius, double soi) {
    BodySystemBuilder b;
    BodyDef root{};
    root.id = 1;
    root.parent_id = InvalidBody;
    root.mu = mu;
    root.radius = radius;
    root.soi_radius = soi;
    b.add_body(root);
    BodySystem sys;
    (void)b.build(sys);
    return sys;
}

}  // namespace

TEST_CASE("FindNextEvent_RadialInboundAboveBody_DetectsImpact",
          "[phase3][radial]") {
    // Spacecraft aimed straight at the central body. With v=(-1, 0) from r=(10, 0)
    // and mu=1, eps = 0.4 (hyperbolic), so it falls directly into the body.
    BodySystem sys = build_root_only(1.0, 0.5, 100.0);
    EventDetector det(sys);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    req.initial_state = State2{{10.0, 0.0}, {-1.0, 0.0}};
    req.time_limit = 50.0;

    PredictedEvent event;
    REQUIRE(det.find_next_event(req, event) == SolveStatus::Ok);
    REQUIRE(event.type == EventType::Impact);
    REQUIRE(event.from_body == 1);
    REQUIRE(event.time > req.start_time);
    REQUIRE(event.time < req.time_limit);
    REQUIRE_THAT(length(event.state.r), WithinAbs(0.5, 1e-4));
}

TEST_CASE("FindNextEvent_RadialOutboundEscape_DetectsSoiExit",
          "[phase3][radial]") {
    // Purely radial hyperbolic escape: r=(2, 0), v=(1.5, 0), mu=1.
    // eps = 1.125 - 0.5 = 0.625 > 0, so spacecraft flies out to the SOI.
    BodySystem sys = build_root_only(1.0, 0.5, 50.0);
    EventDetector det(sys);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    req.initial_state = State2{{2.0, 0.0}, {1.5, 0.0}};
    req.time_limit = 500.0;

    PredictedEvent event;
    REQUIRE(det.find_next_event(req, event) == SolveStatus::Ok);
    REQUIRE(event.type == EventType::SoiExit);
    REQUIRE(event.from_body == 1);
    REQUIRE(event.time > req.start_time);
    REQUIRE_THAT(length(event.state.r), WithinAbs(50.0, 1e-3));
}

TEST_CASE("FindNextEvent_RadialBoundOutbound_ReturnsTimeLimitIfBelowTurningPoint",
          "[phase3][radial]") {
    // Bound radial: r=(10, 0), v=(0.3, 0), mu=1 -> eps=-0.055, turning point ~18.18.
    // Body radius 0.5 and SOI 100. Over a short horizon that stays above the body
    // and below the turning point, no event should fire; expect TimeLimit.
    BodySystem sys = build_root_only(1.0, 0.5, 100.0);
    EventDetector det(sys);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    req.initial_state = State2{{10.0, 0.0}, {0.3, 0.0}};
    req.time_limit = 5.0;

    PredictedEvent event;
    REQUIRE(det.find_next_event(req, event) == SolveStatus::Ok);
    REQUIRE(event.type == EventType::TimeLimit);
    REQUIRE_THAT(event.time, WithinAbs(req.time_limit, 1e-9));
    // End state should still be on the +x axis, moved outward.
    REQUIRE_THAT(event.state.r.y, WithinAbs(0.0, 1e-6));
    REQUIRE(event.state.r.x > 10.0);
}

TEST_CASE("FindNextEvent_RadialInboundOffAxis_DetectsImpact",
          "[phase3][radial]") {
    // Radial along a non-axis direction: r and v are parallel but pointing inward
    // along theta = 0.7 rad. Impact still expected; final state should still be
    // colinear with the initial radial line.
    BodySystem sys = build_root_only(1.0, 0.5, 100.0);
    EventDetector det(sys);

    double theta = 0.7;
    Vec2 dir = {std::cos(theta), std::sin(theta)};
    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    req.initial_state =
        State2{{10.0 * dir.x, 10.0 * dir.y}, {-1.0 * dir.x, -1.0 * dir.y}};
    req.time_limit = 50.0;

    PredictedEvent event;
    REQUIRE(det.find_next_event(req, event) == SolveStatus::Ok);
    REQUIRE(event.type == EventType::Impact);
    REQUIRE_THAT(length(event.state.r), WithinAbs(0.5, 1e-4));
    // Impact point should still be on the original radial line.
    double cross = event.state.r.x * dir.y - event.state.r.y * dir.x;
    REQUIRE_THAT(cross, WithinAbs(0.0, 1e-6));
}
