#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "brahe/body_system.h"
#include "brahe/event_detector.h"
#include "brahe/vec2.h"

#include <cmath>

using namespace brahe;
using Catch::Matchers::WithinAbs;

// --- Section 7.6: End-to-end detector tests: SOI exit ---

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

TEST_CASE("FindNextEvent_ReturnsSoiExitForOutboundTrajectory", "[phase3][soi_exit]") {
    // Hyperbolic outbound: will exit the SOI.
    BodySystem sys = build_root_only(1.0, 0.1, 10.0);
    EventDetector det(sys);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    // v_esc at r=2 is 1.0; use 1.5 perpendicular -> clearly hyperbolic.
    req.initial_state = State2{{2.0, 0.0}, {0.0, 1.5}};
    req.time_limit = 100.0;

    PredictedEvent event;
    REQUIRE(det.find_next_event(req, event) == SolveStatus::Ok);
    REQUIRE(event.type == EventType::SoiExit);
    REQUIRE(event.from_body == 1);
    REQUIRE_THAT(length(event.state.r), WithinAbs(10.0, 1e-4));
}

TEST_CASE("FindNextEvent_DoesNotReturnSoiExitForBoundOrbitInsideLimit", "[phase3][soi_exit]") {
    // Circular orbit well inside SOI.
    BodySystem sys = build_root_only(1.0, 0.1, 100.0);
    EventDetector det(sys);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    req.initial_state = State2{{5.0, 0.0}, {0.0, std::sqrt(1.0 / 5.0)}};  // circular at r=5
    req.time_limit = 50.0;

    PredictedEvent event;
    REQUIRE(det.find_next_event(req, event) == SolveStatus::Ok);
    REQUIRE(event.type != EventType::SoiExit);
    REQUIRE(event.type == EventType::TimeLimit);
}

TEST_CASE("FindNextEvent_SoiExitAtStartHandledPerBoundaryPolicy", "[phase3][soi_exit]") {
    // Policy (spec 6.3): outward radial trend at boundary returns SoiExit immediately.
    // Inward radial trend at boundary does NOT return spurious SoiExit.
    BodySystem sys = build_root_only(1.0, 0.1, 10.0);
    EventDetector det(sys);

    // Case A: at boundary, moving outward -> immediate exit.
    {
        EventSearchRequest req;
        req.central_body = 1;
        req.start_time = 0.0;
        req.initial_state = State2{{10.0, 0.0}, {1.5, 0.0}};  // radially outward at boundary
        req.time_limit = 10.0;

        PredictedEvent event;
        REQUIRE(det.find_next_event(req, event) == SolveStatus::Ok);
        REQUIRE(event.type == EventType::SoiExit);
        REQUIRE_THAT(event.time, WithinAbs(0.0, 1e-9));
    }

    // Case B: at boundary, moving inward -> no spurious immediate exit.
    {
        EventSearchRequest req;
        req.central_body = 1;
        req.start_time = 0.0;
        req.initial_state = State2{{10.0, 0.0}, {-0.5, 0.0}};  // radially inward at boundary
        req.time_limit = 5.0;

        PredictedEvent event;
        REQUIRE(det.find_next_event(req, event) == SolveStatus::Ok);
        REQUIRE(event.time > 0.0);  // event must be strictly after start_time
    }
}

TEST_CASE("FindNextEvent_ReturnsEarlierOfImpactAndSoiExit", "[phase3][soi_exit]") {
    // Spacecraft on a trajectory where SOI exit would occur, but impact happens first.
    BodySystem sys = build_root_only(1.0, 0.5, 5.0);
    EventDetector det(sys);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    // Sub-circular at r=2, tangential velocity v=0.3 (v_circ at r=2 is ~0.707).
    // This creates an elliptic orbit with rp < 0.5 (impact occurs before reaching SOI).
    req.initial_state = State2{{2.0, 0.0}, {0.0, 0.3}};
    req.time_limit = 100.0;

    PredictedEvent event;
    REQUIRE(det.find_next_event(req, event) == SolveStatus::Ok);
    REQUIRE(event.type == EventType::Impact);
}
