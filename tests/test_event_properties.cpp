#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "brahe/body_system.h"
#include "brahe/event_detector.h"
#include "brahe/two_body.h"
#include "brahe/vec2.h"

#include <cmath>

using namespace brahe;
using Catch::Matchers::WithinAbs;

// --- Section 7.11: Property tests ---

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

BodySystem build_with_child(double child_soi, double orbit_radius) {
    BodySystemBuilder b;
    BodyDef root{};
    root.id = 1;
    root.parent_id = InvalidBody;
    root.mu = 1.0;
    root.radius = 0.5;
    root.soi_radius = 1000.0;
    b.add_body(root);

    BodyDef child{};
    child.id = 2;
    child.parent_id = 1;
    child.mu = 0.001;
    child.radius = 0.1;
    child.soi_radius = child_soi;
    child.orbit_radius = orbit_radius;
    child.angular_rate = 0.0;
    child.phase_at_epoch = 0.0;
    b.add_body(child);

    BodySystem sys;
    (void)b.build(sys);
    return sys;
}

}  // namespace

TEST_CASE("Detector_ReturnedEventTimeIsNeverEarlierThanStartTimeMinusEpsilon",
          "[phase3][property]") {
    BodySystem sys = build_root_only(1.0, 1.0, 100.0);
    EventDetector det(sys);
    double eps = 1e-6;

    struct Case {
        double start_time;
        State2 state;
        double time_limit;
    };
    Case cases[] = {
        {0.0, {{3.0, 0.0}, {0.0, 0.3}}, 50.0},  // impact
        {10.0, {{5.0, 0.0}, {0.0, std::sqrt(0.2)}}, 30.0},  // no event -> time limit
        {0.0, {{3.0, 0.0}, {0.0, 2.0}}, 100.0},  // hyperbolic exit
    };

    for (const auto& c : cases) {
        EventSearchRequest req;
        req.central_body = 1;
        req.start_time = c.start_time;
        req.initial_state = c.state;
        req.time_limit = c.time_limit;
        PredictedEvent ev;
        REQUIRE(det.find_next_event(req, ev) == SolveStatus::Ok);
        REQUIRE(ev.time >= c.start_time - eps);
    }
}

TEST_CASE("Detector_ReturnedEventTimeIsNeverLaterThanTimeLimitPlusEpsilonUnlessFailure",
          "[phase3][property]") {
    BodySystem sys = build_root_only(1.0, 1.0, 100.0);
    EventDetector det(sys);
    double eps = 1e-6;

    struct Case {
        double start_time;
        State2 state;
        double time_limit;
    };
    Case cases[] = {
        {0.0, {{3.0, 0.0}, {0.0, 0.3}}, 50.0},
        {0.0, {{5.0, 0.0}, {0.0, std::sqrt(0.2)}}, 100.0},
        {0.0, {{3.0, 0.0}, {0.0, 2.0}}, 100.0},
    };

    for (const auto& c : cases) {
        EventSearchRequest req;
        req.central_body = 1;
        req.start_time = c.start_time;
        req.initial_state = c.state;
        req.time_limit = c.time_limit;
        PredictedEvent ev;
        REQUIRE(det.find_next_event(req, ev) == SolveStatus::Ok);
        REQUIRE(ev.time <= c.time_limit + eps);
    }
}

TEST_CASE("Detector_ReturnedStateMatchesDirectPropagationToEventTime", "[phase3][property]") {
    BodySystem sys = build_root_only(1.0, 1.0, 100.0);
    EventDetector det(sys);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    req.initial_state = State2{{3.0, 0.0}, {0.0, 2.0}};  // hyperbolic
    req.time_limit = 100.0;

    PredictedEvent ev;
    REQUIRE(det.find_next_event(req, ev) == SolveStatus::Ok);

    // Direct propagation from start to event.time must agree with returned state.
    State2 direct;
    REQUIRE(TwoBody::propagate(1.0, req.initial_state, ev.time - req.start_time, direct) ==
            SolveStatus::Ok);
    REQUIRE_THAT(ev.state.r.x, WithinAbs(direct.r.x, 1e-5));
    REQUIRE_THAT(ev.state.r.y, WithinAbs(direct.r.y, 1e-5));
    REQUIRE_THAT(ev.state.v.x, WithinAbs(direct.v.x, 1e-5));
    REQUIRE_THAT(ev.state.v.y, WithinAbs(direct.v.y, 1e-5));
}

TEST_CASE("Detector_EarliestEventIsInvariantAcrossRepeatedRunsOnSamePlatform",
          "[phase3][property]") {
    BodySystem sys = build_with_child(5.0, 50.0);
    EventDetector det(sys);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    req.initial_state = State2{{30.0, 0.1}, {10.0, 0.0}};
    req.time_limit = 10.0;

    PredictedEvent e1, e2, e3;
    REQUIRE(det.find_next_event(req, e1) == SolveStatus::Ok);
    REQUIRE(det.find_next_event(req, e2) == SolveStatus::Ok);
    REQUIRE(det.find_next_event(req, e3) == SolveStatus::Ok);

    // Bit-identical across runs.
    REQUIRE(e1.time == e2.time);
    REQUIRE(e2.time == e3.time);
    REQUIRE(e1.type == e2.type);
    REQUIRE(e2.type == e3.type);
    REQUIRE(e1.to_body == e2.to_body);
    REQUIRE(e2.to_body == e3.to_body);
    REQUIRE(e1.state.r.x == e2.state.r.x);
    REQUIRE(e2.state.r.x == e3.state.r.x);
}

TEST_CASE("Detector_GrazingDetectionNeverReportsMissAsImpactWhenMinimumExceedsRootEpsilon",
          "[phase3][property]") {
    BodySystem sys = build_root_only(1.0, 1.0, 100.0);
    EventDetector det(sys);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    // Orbit whose periapsis is slightly above body radius -- close but clear miss.
    // At r=3, v_y = 0.58 gives an elliptic orbit with rp ~= 1.13 (above radius=1).
    req.initial_state = State2{{3.0, 0.0}, {0.0, 0.58}};
    req.time_limit = 30.0;

    PredictedEvent ev;
    REQUIRE(det.find_next_event(req, ev) == SolveStatus::Ok);
    // No impact should be reported because periapsis is above radius by well over root_epsilon.
    REQUIRE(ev.type != EventType::Impact);
}
