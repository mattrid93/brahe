#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "brahe/body_system.h"
#include "brahe/vec2.h"

#include <cmath>

using namespace brahe;
using Catch::Matchers::WithinAbs;

// star at origin, planet orbit_radius=100 phase=0, moon orbit_radius=10 phase=pi/2
static BodySystem make_transform_system() {
    BodySystemBuilder b;

    BodyDef star;
    star.id = 0;
    star.parent_id = InvalidBody;
    star.mu = 1.0e10;
    star.radius = 50.0;
    star.soi_radius = 1e9;

    BodyDef planet;
    planet.id = 1;
    planet.parent_id = 0;
    planet.mu = 1.0e6;
    planet.radius = 5.0;
    planet.soi_radius = 40.0;
    planet.orbit_radius = 100.0;
    planet.angular_rate = 0.1;
    planet.phase_at_epoch = 0.0;

    BodyDef moon;
    moon.id = 2;
    moon.parent_id = 1;
    moon.mu = 1.0e3;
    moon.radius = 0.5;
    moon.soi_radius = 3.0;
    moon.orbit_radius = 10.0;
    moon.angular_rate = 1.0;
    moon.phase_at_epoch = M_PI / 2.0;

    b.add_body(star);
    b.add_body(planet);
    b.add_body(moon);

    BodySystem sys;
    b.build(sys);
    return sys;
}

// --- E. Frame transform tests ---

TEST_CASE("StateInAncestorFrame_PlanetToRoot_MatchesPositionInParent", "[transforms]") {
    BodySystem sys = make_transform_system();
    double t = 0.0;

    State2 out;
    REQUIRE(sys.state_in_ancestor_frame(1, 0, t, out) == SolveStatus::Ok);

    Vec2 pos = sys.position_in_parent(1, t);
    Vec2 vel = sys.velocity_in_parent(1, t);

    REQUIRE_THAT(out.r.x, WithinAbs(pos.x, 1e-12));
    REQUIRE_THAT(out.r.y, WithinAbs(pos.y, 1e-12));
    REQUIRE_THAT(out.v.x, WithinAbs(vel.x, 1e-12));
    REQUIRE_THAT(out.v.y, WithinAbs(vel.y, 1e-12));
}

TEST_CASE("StateInAncestorFrame_MoonToPlanet_MatchesDirectEphemeris", "[transforms]") {
    BodySystem sys = make_transform_system();
    double t = 0.0;

    State2 out;
    REQUIRE(sys.state_in_ancestor_frame(2, 1, t, out) == SolveStatus::Ok);

    Vec2 pos = sys.position_in_parent(2, t);
    Vec2 vel = sys.velocity_in_parent(2, t);

    REQUIRE_THAT(out.r.x, WithinAbs(pos.x, 1e-12));
    REQUIRE_THAT(out.r.y, WithinAbs(pos.y, 1e-12));
    REQUIRE_THAT(out.v.x, WithinAbs(vel.x, 1e-12));
    REQUIRE_THAT(out.v.y, WithinAbs(vel.y, 1e-12));
}

TEST_CASE("StateInAncestorFrame_MoonToRoot_EqualsPlanetPlusMoonState", "[transforms]") {
    BodySystem sys = make_transform_system();
    double t = 0.0;

    // Moon in root frame
    State2 moon_root;
    REQUIRE(sys.state_in_ancestor_frame(2, 0, t, moon_root) == SolveStatus::Ok);

    // Planet pos/vel in root frame + moon pos/vel in planet frame
    Vec2 planet_pos = sys.position_in_parent(1, t);
    Vec2 planet_vel = sys.velocity_in_parent(1, t);
    Vec2 moon_pos = sys.position_in_parent(2, t);
    Vec2 moon_vel = sys.velocity_in_parent(2, t);

    REQUIRE_THAT(moon_root.r.x, WithinAbs(planet_pos.x + moon_pos.x, 1e-12));
    REQUIRE_THAT(moon_root.r.y, WithinAbs(planet_pos.y + moon_pos.y, 1e-12));
    REQUIRE_THAT(moon_root.v.x, WithinAbs(planet_vel.x + moon_vel.x, 1e-12));
    REQUIRE_THAT(moon_root.v.y, WithinAbs(planet_vel.y + moon_vel.y, 1e-12));
}

TEST_CASE("StateInRootFrame_EqualsStateInAncestorFrameWithRoot", "[transforms]") {
    BodySystem sys = make_transform_system();
    double t = 2.5;

    State2 via_ancestor;
    REQUIRE(sys.state_in_ancestor_frame(2, 0, t, via_ancestor) == SolveStatus::Ok);

    State2 via_root = sys.state_in_root_frame(2, t);

    REQUIRE(via_ancestor.r.x == via_root.r.x);
    REQUIRE(via_ancestor.r.y == via_root.r.y);
    REQUIRE(via_ancestor.v.x == via_root.v.x);
    REQUIRE(via_ancestor.v.y == via_root.v.y);
}

TEST_CASE("Transforms_AreConsistentAcrossMultipleTimes", "[transforms]") {
    BodySystem sys = make_transform_system();

    for (double t : {0.0, 1.0, 10.0, -3.0}) {
        // moon in root = planet in root + moon in planet
        State2 moon_root;
        REQUIRE(sys.state_in_ancestor_frame(2, 0, t, moon_root) == SolveStatus::Ok);

        Vec2 pp = sys.position_in_parent(1, t);
        Vec2 pv = sys.velocity_in_parent(1, t);
        Vec2 mp = sys.position_in_parent(2, t);
        Vec2 mv = sys.velocity_in_parent(2, t);

        REQUIRE_THAT(moon_root.r.x, WithinAbs(pp.x + mp.x, 1e-10));
        REQUIRE_THAT(moon_root.r.y, WithinAbs(pp.y + mp.y, 1e-10));
        REQUIRE_THAT(moon_root.v.x, WithinAbs(pv.x + mv.x, 1e-10));
        REQUIRE_THAT(moon_root.v.y, WithinAbs(pv.y + mv.y, 1e-10));
    }
}

TEST_CASE("AncestorAccumulation_DoesNotDependOnBuildOrder", "[transforms]") {
    auto build_with_order = [](std::vector<int> order) {
        BodyDef star;
        star.id = 0;
        star.parent_id = InvalidBody;
        star.mu = 1.0e10;
        star.radius = 50.0;
        star.soi_radius = 1e9;

        BodyDef planet;
        planet.id = 1;
        planet.parent_id = 0;
        planet.mu = 1.0e6;
        planet.radius = 5.0;
        planet.soi_radius = 40.0;
        planet.orbit_radius = 100.0;
        planet.angular_rate = 0.1;
        planet.phase_at_epoch = 0.0;

        BodyDef moon;
        moon.id = 2;
        moon.parent_id = 1;
        moon.mu = 1.0e3;
        moon.radius = 0.5;
        moon.soi_radius = 3.0;
        moon.orbit_radius = 10.0;
        moon.angular_rate = 1.0;
        moon.phase_at_epoch = M_PI / 2.0;

        BodyDef defs[] = {star, planet, moon};
        BodySystemBuilder b;
        for (int i : order) b.add_body(defs[i]);

        BodySystem sys;
        b.build(sys);
        return sys;
    };

    BodySystem a = build_with_order({0, 1, 2});
    BodySystem b = build_with_order({2, 1, 0});

    double t = 7.3;
    State2 sa = a.state_in_root_frame(2, t);
    State2 sb = b.state_in_root_frame(2, t);

    REQUIRE(sa.r.x == sb.r.x);
    REQUIRE(sa.r.y == sb.r.y);
    REQUIRE(sa.v.x == sb.v.x);
    REQUIRE(sa.v.y == sb.v.y);
}

// --- F. Numeric guardrail tests ---

TEST_CASE("Vec2_Normalize_RejectsNearZeroVector", "[numeric]") {
    Vec2 tiny = {1e-20, 1e-20};
    Vec2 out;
    REQUIRE_FALSE(normalize(tiny, out));
    REQUIRE(out.x == 0.0);
    REQUIRE(out.y == 0.0);

    Vec2 zero = {0.0, 0.0};
    REQUIRE_FALSE(normalize(zero, out));
}

TEST_CASE("Vec2_LengthSquared_AndDot_AreExactForSimpleInputs", "[numeric]") {
    Vec2 a = {3.0, 4.0};
    REQUIRE(length_squared(a) == 25.0);
    REQUIRE(dot(a, a) == 25.0);

    Vec2 b = {1.0, 0.0};
    REQUIRE(dot(a, b) == 3.0);
}

TEST_CASE("StateQueries_DoNotEmitNaN_OnValidInputs", "[numeric]") {
    BodySystem sys = make_transform_system();

    for (BodyId id : {BodyId(0), BodyId(1), BodyId(2)}) {
        for (double t : {0.0, 1.0, -1.0, 1e6, -1e6}) {
            Vec2 pos = sys.position_in_parent(id, t);
            Vec2 vel = sys.velocity_in_parent(id, t);
            REQUIRE(is_finite(pos));
            REQUIRE(is_finite(vel));

            State2 root = sys.state_in_root_frame(id, t);
            REQUIRE(is_finite(root.r));
            REQUIRE(is_finite(root.v));
        }
    }
}
