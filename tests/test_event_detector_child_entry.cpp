#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "brahe/body_system.h"
#include "brahe/event_detector.h"
#include "brahe/vec2.h"

#include <cmath>

using namespace brahe;
using Catch::Matchers::WithinAbs;

// --- Section 7.7: End-to-end detector tests: child SOI entry ---

namespace {

// One stationary child at (+orbit_radius, 0).
BodySystem build_single_stationary_child(BodyId child_id, double child_soi, double orbit_radius) {
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

// Two stationary children at different phases (same orbit_radius).
BodySystem build_two_stationary_children(BodyId id_a, double phase_a, BodyId id_b, double phase_b,
                                         double child_soi, double orbit_radius) {
    BodySystemBuilder b;
    BodyDef root{};
    root.id = 1;
    root.parent_id = InvalidBody;
    root.mu = 1.0;
    root.radius = 0.5;
    root.soi_radius = 1000.0;
    b.add_body(root);

    BodyDef ca{};
    ca.id = id_a;
    ca.parent_id = 1;
    ca.mu = 0.001;
    ca.radius = 0.1;
    ca.soi_radius = child_soi;
    ca.orbit_radius = orbit_radius;
    ca.angular_rate = 0.0;
    ca.phase_at_epoch = phase_a;
    b.add_body(ca);

    BodyDef cb{};
    cb.id = id_b;
    cb.parent_id = 1;
    cb.mu = 0.001;
    cb.radius = 0.1;
    cb.soi_radius = child_soi;
    cb.orbit_radius = orbit_radius;
    cb.angular_rate = 0.0;
    cb.phase_at_epoch = phase_b;
    b.add_body(cb);

    BodySystem sys;
    (void)b.build(sys);
    return sys;
}

}  // namespace

TEST_CASE("FindNextEvent_ReturnsChildSoiEntryWhenTrajectoryIntersectsChildSoi",
          "[phase3][child_entry]") {
    // Child at (50, 0) with SOI=5. Near-straight-line trajectory heading toward it.
    BodySystem sys = build_single_stationary_child(2, 5.0, 50.0);
    EventDetector det(sys);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    req.initial_state = State2{{30.0, 0.1}, {10.0, 0.0}};
    req.time_limit = 10.0;

    PredictedEvent event;
    REQUIRE(det.find_next_event(req, event) == SolveStatus::Ok);
    REQUIRE(event.type == EventType::SoiEntry);
    REQUIRE(event.from_body == 1);
    REQUIRE(event.to_body == 2);
    REQUIRE(event.time > 0.0);
}

TEST_CASE("FindNextEvent_ReturnsEarliestChildEntryAmongMultipleChildren",
          "[phase3][child_entry]") {
    // Two children at same orbit_radius but different phases. Spacecraft heads toward
    // the one at phase=0 (the nearer one).
    BodySystem sys = build_two_stationary_children(2, 0.0, 3, M_PI, 5.0, 50.0);
    // child id=2 at (+50, 0); child id=3 at (-50, 0).
    EventDetector det(sys);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    req.initial_state = State2{{30.0, 0.1}, {10.0, 0.0}};  // heading to (+50, 0)
    req.time_limit = 10.0;

    PredictedEvent event;
    REQUIRE(det.find_next_event(req, event) == SolveStatus::Ok);
    REQUIRE(event.type == EventType::SoiEntry);
    REQUIRE(event.to_body == 2);
}

TEST_CASE("FindNextEvent_UsesLowestBodyIdWhenSiblingEntriesTieWithinTimeEpsilon",
          "[phase3][child_entry]") {
    // Two non-overlapping children on the +X axis. High spacecraft speed makes their
    // entry times differ by less than time_epsilon; lowest BodyId should win.
    BodySystemBuilder b;
    BodyDef root{};
    root.id = 1;
    root.parent_id = InvalidBody;
    root.mu = 1.0;
    root.radius = 0.5;
    root.soi_radius = 1000.0;
    b.add_body(root);

    BodyDef child_a{};
    child_a.id = 10;
    child_a.parent_id = 1;
    child_a.mu = 0.001;
    child_a.radius = 0.1;
    child_a.soi_radius = 1.0;
    child_a.orbit_radius = 50.0;
    child_a.angular_rate = 0.0;
    child_a.phase_at_epoch = 0.0;
    b.add_body(child_a);

    BodyDef child_b = child_a;
    child_b.id = 20;
    child_b.orbit_radius = 52.1;
    b.add_body(child_b);

    BodySystem sys;
    REQUIRE(b.build(sys) == SolveStatus::Ok);
    EventDetector det(sys);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    req.initial_state = State2{{30.0, 0.0}, {3000000.0, 0.0}};
    req.time_limit = 0.00001;

    PredictedEvent event;
    REQUIRE(det.find_next_event(req, event) == SolveStatus::Ok);
    REQUIRE(event.type == EventType::SoiEntry);
    REQUIRE(event.to_body == 10);  // lowest BodyId among tied entries
}

TEST_CASE("FindNextEvent_DoesNotReportFalseChildEntryForNearMiss", "[phase3][child_entry]") {
    // Trajectory passes the child with closest approach well outside child's SOI.
    BodySystem sys = build_single_stationary_child(2, 5.0, 50.0);
    EventDetector det(sys);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    // Spacecraft at (30, 20), v=(10, 0). Closest approach to (50, 0) has y-distance 20 >> 5.
    req.initial_state = State2{{30.0, 20.0}, {10.0, 0.0}};
    req.time_limit = 10.0;

    PredictedEvent event;
    REQUIRE(det.find_next_event(req, event) == SolveStatus::Ok);
    REQUIRE(event.type != EventType::SoiEntry);
}

TEST_CASE("FindNextEvent_DetectsGrazingChildEntry", "[phase3][child_entry]") {
    // Closest approach is just barely inside the child SOI -- distance function dips
    // very briefly below zero and returns. A pure sign-change scan could miss the bracket.
    BodySystem sys = build_single_stationary_child(2, 5.0, 50.0);
    EventDetector det(sys);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    // y-offset 4.9 gives closest-approach ~4.9 when x = 50: inside SOI=5 by a tiny margin.
    req.initial_state = State2{{30.0, 4.9}, {10.0, 0.0}};
    req.time_limit = 10.0;

    PredictedEvent event;
    REQUIRE(det.find_next_event(req, event) == SolveStatus::Ok);
    REQUIRE(event.type == EventType::SoiEntry);
    REQUIRE(event.to_body == 2);
}
