#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "brahe/body_system.h"
#include "brahe/event_detector.h"
#include "brahe/patcher.h"
#include "brahe/vec2.h"

#include <cmath>

using namespace brahe;
using Catch::Matchers::WithinAbs;

// --- Section 7.3: State continuity across patch boundaries ---

namespace {

BodySystem make_sun_moving_planet() {
    BodySystemBuilder b;
    BodyDef sun{};
    sun.id = 1;
    sun.parent_id = InvalidBody;
    sun.mu = 1.0;
    sun.radius = 0.5;
    sun.soi_radius = 1000.0;
    b.add_body(sun);

    BodyDef planet{};
    planet.id = 2;
    planet.parent_id = 1;
    planet.mu = 0.01;
    planet.radius = 0.2;
    planet.soi_radius = 5.0;
    planet.orbit_radius = 50.0;
    planet.angular_rate = 0.08;
    planet.phase_at_epoch = 0.1;
    b.add_body(planet);

    BodySystem sys;
    (void)b.build(sys);
    return sys;
}

}  // namespace

TEST_CASE("PatchBoundary_PositionContinuityWithinPositionEpsilon",
          "[phase4][continuity][patch]") {
    BodySystem sys = make_sun_moving_planet();
    Tolerances tol{};
    EventDetector detector(sys, tol);
    Patcher patcher(sys, tol);

    // Hyperbolic flyby heading toward the planet.
    // Start in Sun frame well outside the planet's SOI, aimed to graze it.
    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    req.initial_state = State2{{30.0, 4.5}, {5.0, 0.0}};
    req.time_limit = 20.0;

    PredictedEvent ev{};
    REQUIRE(detector.find_next_event(req, ev) == SolveStatus::Ok);
    REQUIRE(ev.type == EventType::SoiEntry);
    REQUIRE(ev.to_body == 2);

    // Patch into the child frame, then reconstruct the absolute parent-frame
    // position and verify it matches the detector's reported state.
    State2 sc_child{};
    REQUIRE(patcher.patch_parent_to_child(ev.time, ev.from_body, ev.to_body, ev.state,
                                          sc_child) == SolveStatus::Ok);

    Vec2 r_pl = sys.position_in_parent(ev.to_body, ev.time);
    Vec2 r_recon = sc_child.r + r_pl;

    REQUIRE_THAT(r_recon.x, WithinAbs(ev.state.r.x, tol.position_epsilon));
    REQUIRE_THAT(r_recon.y, WithinAbs(ev.state.r.y, tol.position_epsilon));
}

TEST_CASE("PatchBoundary_VelocityContinuityWithinVelocityEpsilon",
          "[phase4][continuity][patch]") {
    BodySystem sys = make_sun_moving_planet();
    Tolerances tol{};
    EventDetector detector(sys, tol);
    Patcher patcher(sys, tol);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    req.initial_state = State2{{30.0, 4.5}, {5.0, 0.0}};
    req.time_limit = 20.0;

    PredictedEvent ev{};
    REQUIRE(detector.find_next_event(req, ev) == SolveStatus::Ok);
    REQUIRE(ev.type == EventType::SoiEntry);

    State2 sc_child{};
    REQUIRE(patcher.patch_parent_to_child(ev.time, ev.from_body, ev.to_body, ev.state,
                                          sc_child) == SolveStatus::Ok);

    Vec2 v_pl = sys.velocity_in_parent(ev.to_body, ev.time);
    Vec2 v_recon = sc_child.v + v_pl;

    REQUIRE_THAT(v_recon.x, WithinAbs(ev.state.v.x, tol.velocity_epsilon));
    REQUIRE_THAT(v_recon.y, WithinAbs(ev.state.v.y, tol.velocity_epsilon));
}

TEST_CASE("PatchBoundary_ReverseTransitionPreservesAbsoluteState",
          "[phase4][continuity][patch]") {
    // Complement: exercise child -> parent with a realistic state constructed
    // at an SOI boundary in the child frame.
    BodySystem sys = make_sun_moving_planet();
    Tolerances tol{};
    Patcher patcher(sys, tol);

    const double t = 3.7;
    // A spacecraft state at the planet's SOI boundary, moving outward.
    State2 sc_child{{5.0, 0.0}, {0.2, 0.3}};

    State2 sc_parent{};
    REQUIRE(patcher.patch_child_to_parent(t, 2, sc_child, sc_parent) == SolveStatus::Ok);

    // Round-trip back to child frame via the inverse patch.
    State2 sc_child_rt{};
    REQUIRE(patcher.patch_parent_to_child(t, 1, 2, sc_parent, sc_child_rt) ==
            SolveStatus::Ok);

    REQUIRE_THAT(sc_child_rt.r.x, WithinAbs(sc_child.r.x, tol.position_epsilon));
    REQUIRE_THAT(sc_child_rt.r.y, WithinAbs(sc_child.r.y, tol.position_epsilon));
    REQUIRE_THAT(sc_child_rt.v.x, WithinAbs(sc_child.v.x, tol.velocity_epsilon));
    REQUIRE_THAT(sc_child_rt.v.y, WithinAbs(sc_child.v.y, tol.velocity_epsilon));
}
