#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "brahe/body_system.h"
#include "brahe/event_detector.h"
#include "brahe/vec2.h"

#include <cmath>

using namespace brahe;
using Catch::Matchers::WithinAbs;

// --- Section 7.9: Start-state and boundary tests ---
// Spec 6.3 rules:
//   Impact at start:          always return Impact at start_time
//   SOI exit at boundary:     only return if outward radial trend or outside by >root_epsilon
//   Child entry at boundary:  only return if inward trend or inside by >root_epsilon

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

BodySystem build_with_stationary_child(BodyId child_id, double child_soi, double orbit_radius) {
    BodySystemBuilder b;
    BodyDef root{};
    root.id = 1;
    root.parent_id = InvalidBody;
    root.mu = 1.0;
    root.radius = 0.5;
    root.soi_radius = 1000.0;
    b.add_body(root);

    BodyDef child{};
    child.id = child_id;
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

TEST_CASE("Detector_StartExactlyInsideImpactRadiusReturnsImpactImmediately",
          "[phase3][boundary]") {
    BodySystem sys = build_root_only(1.0, 1.0, 100.0);
    EventDetector det(sys);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 2.5;
    req.initial_state = State2{{0.5, 0.0}, {0.0, 0.1}};  // inside body radius
    req.time_limit = 20.0;

    PredictedEvent event;
    REQUIRE(det.find_next_event(req, event) == SolveStatus::Ok);
    REQUIRE(event.type == EventType::Impact);
    REQUIRE_THAT(event.time, WithinAbs(2.5, 1e-9));
}

TEST_CASE("Detector_StartExactlyOnChildBoundaryWithOutwardMotionDoesNotSpuriouslyEnter",
          "[phase3][boundary]") {
    // Child at (+50, 0), SOI=5. Spacecraft at (45, 0) on inner boundary, moving outward from
    // the child (-x direction increases distance to child at (50,0)).
    BodySystem sys = build_with_stationary_child(2, 5.0, 50.0);
    EventDetector det(sys);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    // Distance to (50,0) is 5 (on boundary). Moving away from child -> outward trend.
    req.initial_state = State2{{45.0, 0.0}, {-1.0, 0.0}};
    req.time_limit = 5.0;

    PredictedEvent event;
    REQUIRE(det.find_next_event(req, event) == SolveStatus::Ok);
    // Should NOT return an immediate (t=0) child entry.
    if (event.type == EventType::SoiEntry) {
        REQUIRE(event.time > 1e-9);
    }
}

TEST_CASE("Detector_StartExactlyOnChildBoundaryWithInwardMotionReturnsEntry",
          "[phase3][boundary]") {
    // Distance to child at (50,0) is 5. Moving toward child -> inward trend -> entry at start.
    BodySystem sys = build_with_stationary_child(2, 5.0, 50.0);
    EventDetector det(sys);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    req.initial_state = State2{{45.0, 0.0}, {1.0, 0.0}};  // moving +x toward child
    req.time_limit = 5.0;

    PredictedEvent event;
    REQUIRE(det.find_next_event(req, event) == SolveStatus::Ok);
    REQUIRE(event.type == EventType::SoiEntry);
    REQUIRE(event.to_body == 2);
    REQUIRE_THAT(event.time, WithinAbs(0.0, 1e-9));
}

TEST_CASE("Detector_StartExactlyOnSoiBoundaryWithOutwardMotionReturnsExit",
          "[phase3][boundary]") {
    BodySystem sys = build_root_only(1.0, 0.1, 10.0);
    EventDetector det(sys);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    req.initial_state = State2{{10.0, 0.0}, {1.5, 0.0}};  // radially outward at boundary
    req.time_limit = 5.0;

    PredictedEvent event;
    REQUIRE(det.find_next_event(req, event) == SolveStatus::Ok);
    REQUIRE(event.type == EventType::SoiExit);
    REQUIRE_THAT(event.time, WithinAbs(0.0, 1e-9));
}

TEST_CASE("Detector_StartExactlyOnSoiBoundaryWithInwardMotionDoesNotSpuriouslyExit",
          "[phase3][boundary]") {
    BodySystem sys = build_root_only(1.0, 0.1, 10.0);
    EventDetector det(sys);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    req.initial_state = State2{{10.0, 0.0}, {-0.5, 0.0}};  // radially inward at boundary
    req.time_limit = 5.0;

    PredictedEvent event;
    REQUIRE(det.find_next_event(req, event) == SolveStatus::Ok);
    // Should NOT return an immediate (t=0) SOI exit.
    if (event.type == EventType::SoiExit) {
        REQUIRE(event.time > 1e-9);
    }
}
