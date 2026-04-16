#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "brahe/body_system.h"
#include "brahe/event_detector.h"
#include "brahe/vec2.h"

#include <cmath>

using namespace brahe;
using Catch::Matchers::WithinAbs;

// --- Section 7.5: End-to-end detector tests: impact ---

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

TEST_CASE("FindNextEvent_ReturnsImpactForInboundCollision", "[phase3][impact]") {
    // mu=1, radius=1. Sub-circular tangential orbit -> periapsis < radius -> impact.
    BodySystem sys = build_root_only(1.0, 1.0, 100.0);
    EventDetector det(sys);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    req.initial_state = State2{{3.0, 0.0}, {0.0, 0.3}};
    req.time_limit = 50.0;

    PredictedEvent event;
    REQUIRE(det.find_next_event(req, event) == SolveStatus::Ok);
    REQUIRE(event.type == EventType::Impact);
    REQUIRE(event.from_body == 1);
    // State at event time should be at body radius.
    REQUIRE_THAT(length(event.state.r), WithinAbs(1.0, 1e-5));
}

TEST_CASE("FindNextEvent_ReturnsTimeLimitWhenNoImpactOccursBeforeLimit", "[phase3][impact]") {
    // Circular orbit well outside body radius.
    BodySystem sys = build_root_only(1.0, 1.0, 100.0);
    EventDetector det(sys);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    // Circular at r=5: v = sqrt(1/5).
    req.initial_state = State2{{5.0, 0.0}, {0.0, std::sqrt(1.0 / 5.0)}};
    req.time_limit = 20.0;

    PredictedEvent event;
    REQUIRE(det.find_next_event(req, event) == SolveStatus::Ok);
    REQUIRE(event.type == EventType::TimeLimit);
    REQUIRE_THAT(event.time, WithinAbs(20.0, 1e-9));
}

TEST_CASE("FindNextEvent_ImpactAtStartTimeIsDetected", "[phase3][impact]") {
    // Spacecraft starts at exactly body radius.
    BodySystem sys = build_root_only(1.0, 1.0, 100.0);
    EventDetector det(sys);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 7.5;
    req.initial_state = State2{{1.0, 0.0}, {0.0, 0.5}};  // |r| = 1.0 = body radius
    req.time_limit = 20.0;

    PredictedEvent event;
    REQUIRE(det.find_next_event(req, event) == SolveStatus::Ok);
    REQUIRE(event.type == EventType::Impact);
    REQUIRE_THAT(event.time, WithinAbs(7.5, 1e-9));
}

TEST_CASE("FindNextEvent_PrefersEarlierImpactOverLaterSoiExit", "[phase3][impact]") {
    // Highly eccentric orbit: periapsis < body_radius AND apoapsis > soi_radius
    // (so both events are possible in principle), but the spacecraft starts
    // inbound so impact occurs first.
    BodySystem sys = build_root_only(1.0, 1.0, 8.0);
    EventDetector det(sys);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    // At (4, 0), v=(0, 0.3): E = 0.045 - 0.25 = -0.205, a ~= 2.44
    // h = 4*0.3 = 1.2, p = 1.44, e = sqrt(1 - 1.44/2.44) = sqrt(0.41) ~= 0.64
    // rp ~= 0.88 < 1 (impact), ra ~= 4.0 (well within SOI)
    // Moving tangentially, but v_r just after start goes negative -> inward -> impact.
    req.initial_state = State2{{4.0, 0.0}, {-0.1, 0.3}};  // slight inward radial component
    req.time_limit = 50.0;

    PredictedEvent event;
    REQUIRE(det.find_next_event(req, event) == SolveStatus::Ok);
    REQUIRE(event.type == EventType::Impact);
}
