#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "brahe/body_system.h"
#include "brahe/event_detector.h"

#include <cmath>

using namespace brahe;
using Catch::Matchers::WithinAbs;

// --- Section 7.4: Coarse-scan / candidate-generation tests ---
// Exercised through the full detector using scenarios that specifically
// stress scan bracketing and adaptive step sizing.

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

BodySystem build_root_with_stationary_child(double mu_root, double radius_root, double soi_root,
                                            BodyId child_id, double child_radius, double child_soi,
                                            double child_orbit_radius) {
    BodySystemBuilder b;
    BodyDef root{};
    root.id = 1;
    root.parent_id = InvalidBody;
    root.mu = mu_root;
    root.radius = radius_root;
    root.soi_radius = soi_root;
    b.add_body(root);

    BodyDef child{};
    child.id = child_id;
    child.parent_id = 1;
    child.mu = 0.001;
    child.radius = child_radius;
    child.soi_radius = child_soi;
    child.orbit_radius = child_orbit_radius;
    child.angular_rate = 0.0;
    child.phase_at_epoch = 0.0;  // child at (+orbit_radius, 0)
    b.add_body(child);

    BodySystem sys;
    (void)b.build(sys);
    return sys;
}

}  // namespace

TEST_CASE("CoarseScan_BracketsImpactEventOnSimpleInboundTrajectory", "[phase3][scan]") {
    // mu=1, radius=1, soi=100. Sub-circular tangential velocity drives periapsis below radius.
    BodySystem sys = build_root_only(1.0, 1.0, 100.0);
    EventDetector det(sys);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    // At r=3, v_circ = sqrt(1/3) ~= 0.577. Use tangential v=0.3 (sub-circular) -> impact.
    req.initial_state = State2{{3.0, 0.0}, {0.0, 0.3}};
    req.time_limit = 50.0;

    PredictedEvent event;
    REQUIRE(det.find_next_event(req, event) == SolveStatus::Ok);
    REQUIRE(event.type == EventType::Impact);
    REQUIRE(event.time > 0.0);
    REQUIRE(event.time < 50.0);
}

TEST_CASE("CoarseScan_BracketsSoiExitOnSimpleOutboundTrajectory", "[phase3][scan]") {
    // Hyperbolic outbound: v > escape velocity.
    BodySystem sys = build_root_only(1.0, 0.1, 10.0);
    EventDetector det(sys);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    // At r=2, v_esc = sqrt(2*1/2) = 1.0. v=1.5 perpendicular -> hyperbolic outbound.
    req.initial_state = State2{{2.0, 0.0}, {0.0, 1.5}};
    req.time_limit = 100.0;

    PredictedEvent event;
    REQUIRE(det.find_next_event(req, event) == SolveStatus::Ok);
    REQUIRE(event.type == EventType::SoiExit);
    REQUIRE(event.time > 0.0);
    REQUIRE(event.time < 100.0);
}

TEST_CASE("CoarseScan_BracketsChildEntryOnCrossingTrajectory", "[phase3][scan]") {
    // Root mu=1, child stationary at (50, 0) with SOI=5.
    BodySystem sys =
        build_root_with_stationary_child(1.0, 0.5, 500.0, 2, 0.1, 5.0, 50.0);
    EventDetector det(sys);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    // Near-straight line toward child: high velocity minimizes gravitational curvature.
    req.initial_state = State2{{30.0, 0.1}, {10.0, 0.0}};
    req.time_limit = 20.0;

    PredictedEvent event;
    REQUIRE(det.find_next_event(req, event) == SolveStatus::Ok);
    REQUIRE(event.type == EventType::SoiEntry);
    REQUIRE(event.to_body == 2);
}

TEST_CASE("CoarseScan_FindsGrazingCandidateWithoutSignChange", "[phase3][scan]") {
    // Trajectory passes tangent to child SOI (distance of closest approach ~= child_soi).
    BodySystem sys =
        build_root_with_stationary_child(1.0, 0.5, 500.0, 2, 0.1, 5.0, 50.0);
    EventDetector det(sys);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    // Spacecraft at (30, 4.9), v=(10, 0). Closest approach to child(50,0) is ~4.9 < 5.0.
    // The distance function drops below zero only briefly, producing a near-tangent dip.
    req.initial_state = State2{{30.0, 4.9}, {10.0, 0.0}};
    req.time_limit = 10.0;

    PredictedEvent event;
    REQUIRE(det.find_next_event(req, event) == SolveStatus::Ok);
    REQUIRE(event.type == EventType::SoiEntry);
    REQUIRE(event.to_body == 2);
}

TEST_CASE("CoarseScan_FindsChildEntryWhenCoarseEndpointsStayOutside", "[phase3][scan]") {
    // This is the true no-sign-change grazing case: the spacecraft dips just inside the
    // child SOI between coarse samples, but both coarse endpoints remain outside.
    BodySystem sys =
        build_root_with_stationary_child(1.0, 0.5, 1000.0, 2, 0.1, 5.0, 50.0);
    EventDetector det(sys);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    req.initial_state = State2{{47.5, 4.999}, {100.0, 0.0}};
    req.time_limit = 1.0;

    PredictedEvent event;
    REQUIRE(det.find_next_event(req, event) == SolveStatus::Ok);
    REQUIRE(event.type == EventType::SoiEntry);
    REQUIRE(event.to_body == 2);
    REQUIRE(event.time > 0.0);
    REQUIRE(event.time < 0.05);
}

TEST_CASE("CoarseScan_UsesAdaptiveStepSmallEnoughToNotSkipSmallChildSoi", "[phase3][scan]") {
    // Very small child SOI. Without adaptive stepping, a naive scan could skip it.
    BodySystem sys =
        build_root_with_stationary_child(1.0, 0.5, 1000.0, 2, 0.05, 0.5, 50.0);
    EventDetector det(sys);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    // Near-straight line heading directly at child center (50, 0).
    // Crossing time through SOI ~= 2*0.5 / 10 = 0.1, much smaller than default scan step.
    req.initial_state = State2{{30.0, 0.01}, {10.0, 0.0}};
    req.time_limit = 5.0;

    PredictedEvent event;
    REQUIRE(det.find_next_event(req, event) == SolveStatus::Ok);
    REQUIRE(event.type == EventType::SoiEntry);
    REQUIRE(event.to_body == 2);
}

TEST_CASE("CoarseScan_RespectsMinimumStepFloor", "[phase3][scan]") {
    // Confirm scan completes in bounded time with very small child SOI and high speed.
    // A missing floor could collapse step size to pathological values and hang.
    BodySystem sys =
        build_root_with_stationary_child(1.0, 0.5, 1000.0, 2, 0.01, 0.1, 50.0);
    EventDetector det(sys);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    // Spacecraft nowhere near the child. Should find no event before time_limit.
    req.initial_state = State2{{10.0, 0.0}, {0.0, 0.3}};  // bound orbit, doesn't reach child
    req.time_limit = 5.0;

    PredictedEvent event;
    SolveStatus s = det.find_next_event(req, event);
    REQUIRE(s == SolveStatus::Ok);  // completes without hanging
    REQUIRE(std::isfinite(event.time));
}
