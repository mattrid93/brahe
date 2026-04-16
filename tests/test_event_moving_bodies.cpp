#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "brahe/body_system.h"
#include "brahe/event_detector.h"
#include "brahe/vec2.h"

#include <cmath>

using namespace brahe;
using Catch::Matchers::WithinAbs;

// --- Section 7.13: Moving-body system tests ---
// All prior phase-3 tests use angular_rate = 0.0 for child bodies. These tests
// exercise the time-varying child ephemeris path -- `BodySystem::position_in_parent`
// is explicitly time-dependent, so a detector that uses a snapshot of the child
// position would pass stationary-child tests while failing these.
//
// Shared geometry:
//   Root: mu=1.0, radius=0.5, soi_radius=1000.0
//   Spacecraft ellipse: rp=5, ra=15 => a=10, e=0.5, period T = 2*pi*sqrt(1000) ~= 198.69
//     initial state: (5, 0), v=(0, sqrt(0.3)) ~= (0, 0.5477) -- periapsis, CCW
//     apoapsis at (-15, 0) at t = T/2 ~= 99.346
//   Moon circular orbit at r=15, soi=4 (reaches from r=11 to r=19)

namespace {

constexpr double kScVp = 0.5477225575051661;  // sqrt(0.3), tangential speed at rp=5, mu=1
constexpr double kApoapsisTime = 99.3459;     // T/2 for the above ellipse

BodySystem build_with_one_moon(double orbit_radius, double soi, double angular_rate,
                               double phase_at_epoch) {
    BodySystemBuilder b;
    BodyDef root{};
    root.id = 1;
    root.parent_id = InvalidBody;
    root.mu = 1.0;
    root.radius = 0.5;
    root.soi_radius = 1000.0;
    b.add_body(root);

    BodyDef moon{};
    moon.id = 2;
    moon.parent_id = 1;
    moon.mu = 0.001;
    moon.radius = 0.3;
    moon.soi_radius = soi;
    moon.orbit_radius = orbit_radius;
    moon.angular_rate = angular_rate;
    moon.phase_at_epoch = phase_at_epoch;
    b.add_body(moon);

    BodySystem sys;
    (void)b.build(sys);
    return sys;
}

EventSearchRequest make_default_request() {
    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    req.initial_state = State2{{5.0, 0.0}, {0.0, kScVp}};
    req.time_limit = 110.0;  // a bit beyond apoapsis
    return req;
}

}  // namespace

TEST_CASE("Detector_MovingMoonRotatesIntoPath_EntryOnlyOccursWhenMoonIsMoving",
          "[phase3][moving]") {
    // Moon phase_at_epoch = -0.507 places a stationary moon at ~(13.15, -7.28) which
    // is entirely outside the spacecraft's elliptic orbit (x in [-15, +5]) -- no encounter.
    // With angular_rate = 0.1, the same phase at epoch causes the moon to arrive at
    // (-15, 0) at t ~= 99.35, coinciding with spacecraft apoapsis -- a rendezvous.
    double phase = -0.507;

    SECTION("stationary moon -- no entry") {
        BodySystem sys = build_with_one_moon(15.0, 4.0, /*rate=*/0.0, phase);
        EventDetector det(sys);

        PredictedEvent ev;
        REQUIRE(det.find_next_event(make_default_request(), ev) == SolveStatus::Ok);
        REQUIRE(ev.type != EventType::SoiEntry);
    }

    SECTION("moving moon -- SoiEntry") {
        BodySystem sys = build_with_one_moon(15.0, 4.0, /*rate=*/0.1, phase);
        EventDetector det(sys);

        PredictedEvent ev;
        REQUIRE(det.find_next_event(make_default_request(), ev) == SolveStatus::Ok);
        REQUIRE(ev.type == EventType::SoiEntry);
        REQUIRE(ev.to_body == 2);
        // Rendezvous occurs close to spacecraft apoapsis; entry happens a few seconds
        // before apoapsis since the moon SOI has radius 4.
        REQUIRE(ev.time > 85.0);
        REQUIRE(ev.time < kApoapsisTime + 1.0);
    }
}

TEST_CASE("Detector_MovingMoon_EventStateMatchesTimeEvaluatedMoonPosition",
          "[phase3][moving]") {
    // This is the critical property test: at the returned event time, the distance
    // between the reported spacecraft state and the moon's time-dependent position
    // must equal the moon's SOI radius (within root_epsilon). A detector that
    // snapshotted the moon position at t=0 would fail this test.
    double moon_rate = 0.1;
    double moon_soi = 4.0;
    double orbit_radius = 15.0;
    BodySystem sys = build_with_one_moon(orbit_radius, moon_soi, moon_rate, -0.507);
    EventDetector det(sys);

    PredictedEvent ev;
    REQUIRE(det.find_next_event(make_default_request(), ev) == SolveStatus::Ok);
    REQUIRE(ev.type == EventType::SoiEntry);
    REQUIRE(ev.to_body == 2);

    // Evaluate moon position at the event time and check it is exactly SOI away.
    Vec2 moon_at_event = sys.position_in_parent(2, ev.time);
    Vec2 delta = ev.state.r - moon_at_event;
    double dist = length(delta);
    REQUIRE_THAT(dist, WithinAbs(moon_soi, 1e-3));
}

TEST_CASE("Detector_StationaryVsMovingMoon_ProduceDifferentEventTimesForSameTrajectory",
          "[phase3][moving]") {
    // Stationary moon at phase=pi sits at (-15, 0); entry (SOI=4) occurs around t~77.5.
    // Moving moon rate=0.1 phase=-0.507 places moon at (-15, 0) at the *moment of apoapsis*;
    // entry occurs later (t~96). Both are SoiEntry events, but the times must differ.
    EventSearchRequest req = make_default_request();
    req.time_limit = 150.0;

    PredictedEvent ev_static;
    {
        BodySystem sys = build_with_one_moon(15.0, 4.0, /*rate=*/0.0, /*phase=*/M_PI);
        EventDetector det(sys);
        REQUIRE(det.find_next_event(req, ev_static) == SolveStatus::Ok);
        REQUIRE(ev_static.type == EventType::SoiEntry);
        REQUIRE(ev_static.to_body == 2);
    }

    PredictedEvent ev_moving;
    {
        BodySystem sys = build_with_one_moon(15.0, 4.0, /*rate=*/0.1, /*phase=*/-0.507);
        EventDetector det(sys);
        REQUIRE(det.find_next_event(req, ev_moving) == SolveStatus::Ok);
        REQUIRE(ev_moving.type == EventType::SoiEntry);
        REQUIRE(ev_moving.to_body == 2);
    }

    // Static case: entry near t=77.5 (symmetric arc across apoapsis at 99.35).
    // Moving case: entry near t=96 (moon "catches up" near apoapsis).
    REQUIRE(ev_moving.time > ev_static.time + 10.0);
}

TEST_CASE("Detector_MovingMoon_DeterminismAcrossRepeatedRuns", "[phase3][moving]") {
    // Bit-identical output across repeated calls when the system has a moving child.
    BodySystem sys = build_with_one_moon(15.0, 4.0, 0.1, -0.507);
    EventDetector det(sys);
    EventSearchRequest req = make_default_request();

    PredictedEvent e1, e2, e3;
    REQUIRE(det.find_next_event(req, e1) == SolveStatus::Ok);
    REQUIRE(det.find_next_event(req, e2) == SolveStatus::Ok);
    REQUIRE(det.find_next_event(req, e3) == SolveStatus::Ok);

    REQUIRE(e1.time == e2.time);
    REQUIRE(e2.time == e3.time);
    REQUIRE(e1.type == e2.type);
    REQUIRE(e2.type == e3.type);
    REQUIRE(e1.to_body == e2.to_body);
    REQUIRE(e2.to_body == e3.to_body);
    REQUIRE(e1.state.r.x == e2.state.r.x);
    REQUIRE(e2.state.r.x == e3.state.r.x);
    REQUIRE(e1.state.r.y == e2.state.r.y);
    REQUIRE(e1.state.v.x == e2.state.v.x);
    REQUIRE(e1.state.v.y == e2.state.v.y);
}

TEST_CASE("Detector_RetrogradeMoon_EntryStillDetected", "[phase3][moving]") {
    // Negative angular_rate (retrograde moon) with a phase chosen so a rendezvous
    // exists in the search window. Because the moon now sweeps clockwise while the
    // spacecraft moves CCW, the first encounter generally occurs earlier than the
    // prograde case -- the two meet head-on rather than the moon having to catch up.
    // We only assert that an entry is detected, that it falls strictly inside the
    // time window, and that the reported state matches the time-evaluated moon
    // position (the critical moving-ephemeris property).
    double rate = -0.1;
    double phase = std::fmod(M_PI - rate * kApoapsisTime, 2.0 * M_PI);
    BodySystem sys = build_with_one_moon(15.0, 4.0, rate, phase);
    EventDetector det(sys);

    EventSearchRequest req = make_default_request();
    PredictedEvent ev;
    REQUIRE(det.find_next_event(req, ev) == SolveStatus::Ok);
    REQUIRE(ev.type == EventType::SoiEntry);
    REQUIRE(ev.to_body == 2);
    REQUIRE(ev.time > req.start_time);
    REQUIRE(ev.time <= req.time_limit);

    // Property: state at event time is exactly moon_soi away from the moon's
    // time-evaluated position (not its t=0 snapshot).
    Vec2 moon_at_event = sys.position_in_parent(2, ev.time);
    double dist = length(ev.state.r - moon_at_event);
    REQUIRE_THAT(dist, WithinAbs(4.0, 1e-3));
}

TEST_CASE("Detector_TwoOrbitingMoons_CorrectMoonReportedByEarliestEntry",
          "[phase3][moving]") {
    // Two moons at the same orbit_radius and rate, different phases. The one that
    // rotates into rendezvous first should be reported.
    //
    // Moon A (id=2): phase chosen so it reaches (-15, 0) at t ~= 60 -- an "early" encounter.
    //   phase_A + 0.1*60 = pi => phase_A = pi - 6 ~= -2.858
    // Moon B (id=3): phase chosen so it reaches (-15, 0) at t ~= 95 -- a "late" encounter.
    //   phase_B + 0.1*95 = pi => phase_B = pi - 9.5 ~= -6.358 (~= -0.075 mod 2*pi)
    //
    // Caveat: moon A's rendezvous point at t=60 must actually be near the spacecraft.
    // Spacecraft at t=60 is *not* at (-15, 0); it is somewhere between perihelion and
    // apoapsis. So moon A at (-15, 0) at t=60 is not automatically an encounter.
    //
    // Instead we calibrate using the apoapsis time more carefully:
    // Let moon A rendezvous with spacecraft near spacecraft-apoapsis at t_A = 99.
    // Let moon B rendezvous with spacecraft near spacecraft-apoapsis at t_B = 100.
    // Tiny phase differences -- but detector must still pick moon A (lower ID already,
    // so this tests arbitration plus moving-body ephemeris).
    BodySystemBuilder b;
    BodyDef root{};
    root.id = 1;
    root.parent_id = InvalidBody;
    root.mu = 1.0;
    root.radius = 0.5;
    root.soi_radius = 1000.0;
    b.add_body(root);

    // Moon A arrives at (-15, 0) slightly before spacecraft apoapsis.
    BodyDef moon_a{};
    moon_a.id = 2;
    moon_a.parent_id = 1;
    moon_a.mu = 0.001;
    moon_a.radius = 0.3;
    moon_a.soi_radius = 4.0;
    moon_a.orbit_radius = 15.0;
    moon_a.angular_rate = 0.1;
    moon_a.phase_at_epoch = M_PI - 0.1 * (kApoapsisTime - 2.0);  // arrives at t~97.35
    b.add_body(moon_a);

    // Moon B arrives at (-15, 0) slightly after spacecraft apoapsis -- so encounter
    // should be later than moon A.
    BodyDef moon_b{};
    moon_b.id = 3;
    moon_b.parent_id = 1;
    moon_b.mu = 0.001;
    moon_b.radius = 0.3;
    moon_b.soi_radius = 4.0;
    moon_b.orbit_radius = 15.0;
    moon_b.angular_rate = 0.1;
    moon_b.phase_at_epoch = M_PI - 0.1 * (kApoapsisTime + 5.0);  // arrives at t~104.35
    b.add_body(moon_b);

    BodySystem sys;
    REQUIRE(b.build(sys) == SolveStatus::Ok);
    EventDetector det(sys);

    EventSearchRequest req = make_default_request();
    req.time_limit = 110.0;

    PredictedEvent ev;
    REQUIRE(det.find_next_event(req, ev) == SolveStatus::Ok);
    REQUIRE(ev.type == EventType::SoiEntry);
    // Moon A is the earlier encounter.
    REQUIRE(ev.to_body == 2);
}
