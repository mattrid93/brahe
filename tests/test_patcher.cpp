#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "brahe/body_system.h"
#include "brahe/patcher.h"
#include "brahe/vec2.h"

#include <cmath>
#include <limits>

using namespace brahe;
using Catch::Matchers::WithinAbs;

// --- Section 7.2: Patch-math unit tests ---

namespace {

constexpr BodyId kSun = 1;
constexpr BodyId kPlanet = 2;
constexpr BodyId kOtherPlanet = 3;
constexpr BodyId kMoon = 4;  // child of kPlanet

// Sun with a stationary planet at (+50, 0). Chosen so the child ephemeris
// evaluates to an easily-verifiable absolute state at every t.
BodySystem make_sun_stationary_planet() {
    BodySystemBuilder b;
    BodyDef sun{};
    sun.id = kSun;
    sun.parent_id = InvalidBody;
    sun.mu = 1.0;
    sun.radius = 0.5;
    sun.soi_radius = 1000.0;
    b.add_body(sun);

    BodyDef planet{};
    planet.id = kPlanet;
    planet.parent_id = kSun;
    planet.mu = 0.01;
    planet.radius = 0.2;
    planet.soi_radius = 5.0;
    planet.orbit_radius = 50.0;
    planet.angular_rate = 0.0;  // stationary
    planet.phase_at_epoch = 0.0;
    b.add_body(planet);

    BodySystem sys;
    (void)b.build(sys);
    return sys;
}

// Sun with a moving planet (angular_rate != 0).
BodySystem make_sun_moving_planet() {
    BodySystemBuilder b;
    BodyDef sun{};
    sun.id = kSun;
    sun.parent_id = InvalidBody;
    sun.mu = 1.0;
    sun.radius = 0.5;
    sun.soi_radius = 1000.0;
    b.add_body(sun);

    BodyDef planet{};
    planet.id = kPlanet;
    planet.parent_id = kSun;
    planet.mu = 0.01;
    planet.radius = 0.2;
    planet.soi_radius = 5.0;
    planet.orbit_radius = 50.0;
    planet.angular_rate = 0.1;
    planet.phase_at_epoch = 0.25;
    b.add_body(planet);

    BodySystem sys;
    (void)b.build(sys);
    return sys;
}

// Sun with two planets (for invalid parent/child relationship tests).
BodySystem make_sun_two_planets() {
    BodySystemBuilder b;
    BodyDef sun{};
    sun.id = kSun;
    sun.parent_id = InvalidBody;
    sun.mu = 1.0;
    sun.radius = 0.5;
    sun.soi_radius = 1000.0;
    b.add_body(sun);

    BodyDef p2{};
    p2.id = kPlanet;
    p2.parent_id = kSun;
    p2.mu = 0.01;
    p2.radius = 0.2;
    p2.soi_radius = 5.0;
    p2.orbit_radius = 50.0;
    p2.angular_rate = 0.0;
    p2.phase_at_epoch = 0.0;
    b.add_body(p2);

    BodyDef p3{};
    p3.id = kOtherPlanet;
    p3.parent_id = kSun;
    p3.mu = 0.01;
    p3.radius = 0.2;
    p3.soi_radius = 5.0;
    p3.orbit_radius = 80.0;
    p3.angular_rate = 0.0;
    p3.phase_at_epoch = 0.0;
    b.add_body(p3);

    BodySystem sys;
    (void)b.build(sys);
    return sys;
}

// Nested hierarchy: sun -> planet -> moon.
BodySystem make_sun_planet_moon() {
    BodySystemBuilder b;
    BodyDef sun{};
    sun.id = kSun;
    sun.parent_id = InvalidBody;
    sun.mu = 1.0;
    sun.radius = 0.5;
    sun.soi_radius = 1000.0;
    b.add_body(sun);

    BodyDef planet{};
    planet.id = kPlanet;
    planet.parent_id = kSun;
    planet.mu = 0.01;
    planet.radius = 0.2;
    planet.soi_radius = 5.0;
    planet.orbit_radius = 50.0;
    planet.angular_rate = 0.0;
    planet.phase_at_epoch = 0.0;
    b.add_body(planet);

    BodyDef moon{};
    moon.id = kMoon;
    moon.parent_id = kPlanet;
    moon.mu = 0.0001;
    moon.radius = 0.05;
    moon.soi_radius = 0.3;
    moon.orbit_radius = 2.0;
    moon.angular_rate = 0.0;
    moon.phase_at_epoch = 0.0;
    b.add_body(moon);

    BodySystem sys;
    (void)b.build(sys);
    return sys;
}

bool is_finite_state(const State2& s) {
    return std::isfinite(s.r.x) && std::isfinite(s.r.y) &&
           std::isfinite(s.v.x) && std::isfinite(s.v.y);
}

}  // namespace

// --- Parent -> child ---

TEST_CASE("PatchParentToChild_SubtractsChildPositionAndVelocity",
          "[phase4][patch][parent_to_child]") {
    // Stationary planet at (50, 0), v=(0,0).
    BodySystem sys = make_sun_stationary_planet();
    Patcher patcher(sys);

    State2 sc_parent{{52.0, 3.0}, {0.5, 1.2}};
    State2 sc_child{};
    REQUIRE(patcher.patch_parent_to_child(0.0, kSun, kPlanet, sc_parent, sc_child) ==
            SolveStatus::Ok);

    // Exact subtraction: (52,3) - (50,0) = (2,3); (0.5,1.2) - (0,0) = (0.5,1.2).
    REQUIRE_THAT(sc_child.r.x, WithinAbs(2.0, 1e-15));
    REQUIRE_THAT(sc_child.r.y, WithinAbs(3.0, 1e-15));
    REQUIRE_THAT(sc_child.v.x, WithinAbs(0.5, 1e-15));
    REQUIRE_THAT(sc_child.v.y, WithinAbs(1.2, 1e-15));
}

TEST_CASE("PatchParentToChild_SubtractsMovingChildStateAtEventTime",
          "[phase4][patch][parent_to_child]") {
    // Moving planet: at a specific t we can compute its state analytically.
    BodySystem sys = make_sun_moving_planet();
    Patcher patcher(sys);

    const double t = 2.0;
    Vec2 r_pl = sys.position_in_parent(kPlanet, t);
    Vec2 v_pl = sys.velocity_in_parent(kPlanet, t);

    // Arbitrary spacecraft state in parent frame at time t.
    State2 sc_parent{{r_pl.x + 1.5, r_pl.y - 2.25}, {v_pl.x + 0.3, v_pl.y + 0.7}};
    State2 sc_child{};
    REQUIRE(patcher.patch_parent_to_child(t, kSun, kPlanet, sc_parent, sc_child) ==
            SolveStatus::Ok);

    REQUIRE_THAT(sc_child.r.x, WithinAbs(1.5, 1e-12));
    REQUIRE_THAT(sc_child.r.y, WithinAbs(-2.25, 1e-12));
    REQUIRE_THAT(sc_child.v.x, WithinAbs(0.3, 1e-12));
    REQUIRE_THAT(sc_child.v.y, WithinAbs(0.7, 1e-12));
}

TEST_CASE("PatchParentToChild_PreservesAbsoluteInertialState",
          "[phase4][patch][parent_to_child][continuity]") {
    BodySystem sys = make_sun_moving_planet();
    Patcher patcher(sys);

    const double t = 4.3;
    State2 sc_parent_original{{30.0, 20.0}, {-0.2, 0.9}};
    State2 sc_child{};
    REQUIRE(patcher.patch_parent_to_child(t, kSun, kPlanet, sc_parent_original, sc_child) ==
            SolveStatus::Ok);

    // Reconstruct absolute parent-frame state by adding the child ephemeris at t.
    Vec2 r_pl = sys.position_in_parent(kPlanet, t);
    Vec2 v_pl = sys.velocity_in_parent(kPlanet, t);
    State2 sc_parent_reconstructed{sc_child.r + r_pl, sc_child.v + v_pl};

    REQUIRE_THAT(sc_parent_reconstructed.r.x, WithinAbs(sc_parent_original.r.x, 1e-12));
    REQUIRE_THAT(sc_parent_reconstructed.r.y, WithinAbs(sc_parent_original.r.y, 1e-12));
    REQUIRE_THAT(sc_parent_reconstructed.v.x, WithinAbs(sc_parent_original.v.x, 1e-12));
    REQUIRE_THAT(sc_parent_reconstructed.v.y, WithinAbs(sc_parent_original.v.y, 1e-12));
}

// --- Child -> parent ---

TEST_CASE("PatchChildToParent_AddsChildPositionAndVelocity",
          "[phase4][patch][child_to_parent]") {
    BodySystem sys = make_sun_moving_planet();
    Patcher patcher(sys);

    const double t = 1.7;
    Vec2 r_pl = sys.position_in_parent(kPlanet, t);
    Vec2 v_pl = sys.velocity_in_parent(kPlanet, t);

    State2 sc_child{{3.0, -1.5}, {0.8, -0.4}};
    State2 sc_parent{};
    REQUIRE(patcher.patch_child_to_parent(t, kPlanet, sc_child, sc_parent) ==
            SolveStatus::Ok);

    REQUIRE_THAT(sc_parent.r.x, WithinAbs(3.0 + r_pl.x, 1e-12));
    REQUIRE_THAT(sc_parent.r.y, WithinAbs(-1.5 + r_pl.y, 1e-12));
    REQUIRE_THAT(sc_parent.v.x, WithinAbs(0.8 + v_pl.x, 1e-12));
    REQUIRE_THAT(sc_parent.v.y, WithinAbs(-0.4 + v_pl.y, 1e-12));
}

TEST_CASE("PatchChildToParent_PreservesAbsoluteInertialState",
          "[phase4][patch][child_to_parent][continuity]") {
    BodySystem sys = make_sun_moving_planet();
    Patcher patcher(sys);

    const double t = 3.14;
    State2 sc_child_original{{2.2, -0.5}, {-0.4, 0.6}};
    State2 sc_parent{};
    REQUIRE(patcher.patch_child_to_parent(t, kPlanet, sc_child_original, sc_parent) ==
            SolveStatus::Ok);

    // Reconstruct child-frame state by subtracting the child ephemeris at t.
    Vec2 r_pl = sys.position_in_parent(kPlanet, t);
    Vec2 v_pl = sys.velocity_in_parent(kPlanet, t);
    State2 sc_child_reconstructed{sc_parent.r - r_pl, sc_parent.v - v_pl};

    REQUIRE_THAT(sc_child_reconstructed.r.x, WithinAbs(sc_child_original.r.x, 1e-12));
    REQUIRE_THAT(sc_child_reconstructed.r.y, WithinAbs(sc_child_original.r.y, 1e-12));
    REQUIRE_THAT(sc_child_reconstructed.v.x, WithinAbs(sc_child_original.v.x, 1e-12));
    REQUIRE_THAT(sc_child_reconstructed.v.y, WithinAbs(sc_child_original.v.y, 1e-12));
}

// --- Round-trip ---

TEST_CASE("PatchParentToChildThenBackToParent_PreservesStateWithinTolerance",
          "[phase4][patch][roundtrip]") {
    BodySystem sys = make_sun_moving_planet();
    Patcher patcher(sys);

    const double t = 2.718;
    State2 sc_parent_original{{45.0, 12.5}, {0.1, -0.3}};

    State2 sc_child{};
    REQUIRE(patcher.patch_parent_to_child(t, kSun, kPlanet, sc_parent_original, sc_child) ==
            SolveStatus::Ok);

    State2 sc_parent_roundtrip{};
    REQUIRE(patcher.patch_child_to_parent(t, kPlanet, sc_child, sc_parent_roundtrip) ==
            SolveStatus::Ok);

    REQUIRE_THAT(sc_parent_roundtrip.r.x, WithinAbs(sc_parent_original.r.x, 1e-12));
    REQUIRE_THAT(sc_parent_roundtrip.r.y, WithinAbs(sc_parent_original.r.y, 1e-12));
    REQUIRE_THAT(sc_parent_roundtrip.v.x, WithinAbs(sc_parent_original.v.x, 1e-12));
    REQUIRE_THAT(sc_parent_roundtrip.v.y, WithinAbs(sc_parent_original.v.y, 1e-12));
}

TEST_CASE("PatchChildToParentThenBackToChild_PreservesStateWithinTolerance",
          "[phase4][patch][roundtrip]") {
    BodySystem sys = make_sun_moving_planet();
    Patcher patcher(sys);

    const double t = 5.55;
    State2 sc_child_original{{2.5, -1.8}, {0.07, 0.22}};

    State2 sc_parent{};
    REQUIRE(patcher.patch_child_to_parent(t, kPlanet, sc_child_original, sc_parent) ==
            SolveStatus::Ok);

    State2 sc_child_roundtrip{};
    REQUIRE(patcher.patch_parent_to_child(t, kSun, kPlanet, sc_parent, sc_child_roundtrip) ==
            SolveStatus::Ok);

    REQUIRE_THAT(sc_child_roundtrip.r.x, WithinAbs(sc_child_original.r.x, 1e-12));
    REQUIRE_THAT(sc_child_roundtrip.r.y, WithinAbs(sc_child_original.r.y, 1e-12));
    REQUIRE_THAT(sc_child_roundtrip.v.x, WithinAbs(sc_child_original.v.x, 1e-12));
    REQUIRE_THAT(sc_child_roundtrip.v.y, WithinAbs(sc_child_original.v.y, 1e-12));
}

// --- Numeric guards ---

TEST_CASE("PatchFunctions_RejectUnknownBodies", "[phase4][patch][invalid]") {
    BodySystem sys = make_sun_stationary_planet();
    Patcher patcher(sys);

    State2 sc{{1.0, 2.0}, {0.3, 0.4}};
    State2 out{};

    // Unknown parent.
    REQUIRE(patcher.patch_parent_to_child(0.0, 999, kPlanet, sc, out) ==
            SolveStatus::InvalidInput);
    // Unknown child.
    REQUIRE(patcher.patch_parent_to_child(0.0, kSun, 999, sc, out) ==
            SolveStatus::InvalidInput);
    // InvalidBody sentinels.
    REQUIRE(patcher.patch_parent_to_child(0.0, InvalidBody, kPlanet, sc, out) ==
            SolveStatus::InvalidInput);
    REQUIRE(patcher.patch_parent_to_child(0.0, kSun, InvalidBody, sc, out) ==
            SolveStatus::InvalidInput);

    // Unknown child for child->parent.
    REQUIRE(patcher.patch_child_to_parent(0.0, 999, sc, out) == SolveStatus::InvalidInput);
    REQUIRE(patcher.patch_child_to_parent(0.0, InvalidBody, sc, out) ==
            SolveStatus::InvalidInput);
}

TEST_CASE("PatchFunctions_RejectInvalidParentChildRelationship", "[phase4][patch][invalid]") {
    // parent_to_child where child_body is not actually a child of parent_body.
    BodySystem sys = make_sun_planet_moon();
    Patcher patcher(sys);

    State2 sc{{1.0, 0.0}, {0.0, 1.0}};
    State2 out{};

    // kMoon is a grandchild of kSun, not a direct child.
    REQUIRE(patcher.patch_parent_to_child(0.0, kSun, kMoon, sc, out) ==
            SolveStatus::InvalidInput);

    // Two siblings: kPlanet is not a parent of kOtherPlanet, and vice versa.
    BodySystem sys2 = make_sun_two_planets();
    Patcher patcher2(sys2);
    REQUIRE(patcher2.patch_parent_to_child(0.0, kPlanet, kOtherPlanet, sc, out) ==
            SolveStatus::InvalidInput);

    // child_to_parent for the root body: root has no parent.
    REQUIRE(patcher.patch_child_to_parent(0.0, kSun, sc, out) == SolveStatus::InvalidInput);
}

TEST_CASE("PatchFunctions_RejectNonFiniteTimeOrState", "[phase4][patch][invalid]") {
    BodySystem sys = make_sun_moving_planet();
    Patcher patcher(sys);

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();

    State2 good{{10.0, 0.0}, {0.0, 1.0}};
    State2 nan_r{{nan, 0.0}, {0.0, 1.0}};
    State2 nan_v{{10.0, 0.0}, {0.0, nan}};
    State2 inf_r{{inf, 0.0}, {0.0, 1.0}};
    State2 out{};

    REQUIRE(patcher.patch_parent_to_child(nan, kSun, kPlanet, good, out) ==
            SolveStatus::InvalidInput);
    REQUIRE(patcher.patch_parent_to_child(inf, kSun, kPlanet, good, out) ==
            SolveStatus::InvalidInput);
    REQUIRE(patcher.patch_parent_to_child(0.0, kSun, kPlanet, nan_r, out) ==
            SolveStatus::InvalidInput);
    REQUIRE(patcher.patch_parent_to_child(0.0, kSun, kPlanet, nan_v, out) ==
            SolveStatus::InvalidInput);
    REQUIRE(patcher.patch_parent_to_child(0.0, kSun, kPlanet, inf_r, out) ==
            SolveStatus::InvalidInput);

    REQUIRE(patcher.patch_child_to_parent(nan, kPlanet, good, out) ==
            SolveStatus::InvalidInput);
    REQUIRE(patcher.patch_child_to_parent(0.0, kPlanet, nan_r, out) ==
            SolveStatus::InvalidInput);
}

TEST_CASE("PatchFunctions_DoNotEmitNaNOnOrdinaryInvalidInput", "[phase4][patch][invalid]") {
    BodySystem sys = make_sun_moving_planet();
    Patcher patcher(sys);

    const double nan = std::numeric_limits<double>::quiet_NaN();
    State2 sc{{10.0, 0.0}, {0.0, 1.0}};
    State2 bad{{nan, nan}, {nan, nan}};

    // Pre-initialise output so that failure modes can't produce NaN by leaving it
    // uninitialised. Any deliberate write must also be finite.
    State2 out{{0.0, 0.0}, {0.0, 0.0}};
    (void)patcher.patch_parent_to_child(0.0, 999, kPlanet, sc, out);
    REQUIRE(is_finite_state(out));

    out = State2{{0.0, 0.0}, {0.0, 0.0}};
    (void)patcher.patch_parent_to_child(nan, kSun, kPlanet, sc, out);
    REQUIRE(is_finite_state(out));

    out = State2{{0.0, 0.0}, {0.0, 0.0}};
    (void)patcher.patch_parent_to_child(0.0, kSun, kPlanet, bad, out);
    REQUIRE(is_finite_state(out));

    out = State2{{0.0, 0.0}, {0.0, 0.0}};
    (void)patcher.patch_child_to_parent(0.0, 999, sc, out);
    REQUIRE(is_finite_state(out));
}
