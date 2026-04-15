#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "brahe/body_system.h"
#include "brahe/vec2.h"

#include <cmath>

using namespace brahe;
using Catch::Matchers::WithinAbs;

// Build a system with known orbit params for analytic verification
// planet: orbit_radius=10, angular_rate=2, phase_at_epoch=0.25
static BodySystem make_test_system() {
    BodySystemBuilder b;

    BodyDef star;
    star.id = 0;
    star.parent_id = InvalidBody;
    star.mu = 1.0e10;
    star.radius = 100.0;
    star.soi_radius = 1e9;

    BodyDef planet;
    planet.id = 1;
    planet.parent_id = 0;
    planet.mu = 1.0e6;
    planet.radius = 1.0;
    planet.soi_radius = 2.0;
    planet.orbit_radius = 10.0;
    planet.angular_rate = 2.0;
    planet.phase_at_epoch = 0.25;

    b.add_body(star);
    b.add_body(planet);

    BodySystem sys;
    b.build(sys);
    return sys;
}

// --- D. Circular ephemeris tests ---

TEST_CASE("PositionInParent_AtEpoch_MatchesAnalyticFormula", "[ephemeris]") {
    BodySystem sys = make_test_system();
    Vec2 pos = sys.position_in_parent(1, 0.0);

    double r = 10.0;
    double phase = 0.25;
    REQUIRE_THAT(pos.x, WithinAbs(r * std::cos(phase), 1e-12));
    REQUIRE_THAT(pos.y, WithinAbs(r * std::sin(phase), 1e-12));
}

TEST_CASE("VelocityInParent_AtEpoch_MatchesAnalyticFormula", "[ephemeris]") {
    BodySystem sys = make_test_system();
    Vec2 vel = sys.velocity_in_parent(1, 0.0);

    double r = 10.0;
    double w = 2.0;
    double phase = 0.25;
    REQUIRE_THAT(vel.x, WithinAbs(r * w * (-std::sin(phase)), 1e-12));
    REQUIRE_THAT(vel.y, WithinAbs(r * w * std::cos(phase), 1e-12));
}

TEST_CASE("PositionInParent_AtPositiveTime_MatchesAnalyticFormula", "[ephemeris]") {
    BodySystem sys = make_test_system();
    double t = 3.7;
    Vec2 pos = sys.position_in_parent(1, t);

    double r = 10.0;
    double w = 2.0;
    double phase = 0.25;
    double theta = phase + w * t;
    REQUIRE_THAT(pos.x, WithinAbs(r * std::cos(theta), 1e-12));
    REQUIRE_THAT(pos.y, WithinAbs(r * std::sin(theta), 1e-12));
}

TEST_CASE("VelocityInParent_AtPositiveTime_MatchesAnalyticFormula", "[ephemeris]") {
    BodySystem sys = make_test_system();
    double t = 3.7;
    Vec2 vel = sys.velocity_in_parent(1, t);

    double r = 10.0;
    double w = 2.0;
    double phase = 0.25;
    double theta = phase + w * t;
    REQUIRE_THAT(vel.x, WithinAbs(r * w * (-std::sin(theta)), 1e-12));
    REQUIRE_THAT(vel.y, WithinAbs(r * w * std::cos(theta), 1e-12));
}

TEST_CASE("RootBody_PositionAndVelocity_AreAlwaysZero", "[ephemeris]") {
    BodySystem sys = make_test_system();

    for (double t : {-100.0, 0.0, 1.0, 999.0}) {
        Vec2 pos = sys.position_in_parent(0, t);
        Vec2 vel = sys.velocity_in_parent(0, t);
        REQUIRE(pos.x == 0.0);
        REQUIRE(pos.y == 0.0);
        REQUIRE(vel.x == 0.0);
        REQUIRE(vel.y == 0.0);
    }
}

TEST_CASE("NegativeTime_EphemerisStillMatchesAnalyticFormula", "[ephemeris]") {
    BodySystem sys = make_test_system();
    double t = -2.5;
    Vec2 pos = sys.position_in_parent(1, t);

    double r = 10.0;
    double w = 2.0;
    double phase = 0.25;
    double theta = phase + w * t;
    REQUIRE_THAT(pos.x, WithinAbs(r * std::cos(theta), 1e-12));
    REQUIRE_THAT(pos.y, WithinAbs(r * std::sin(theta), 1e-12));
}

TEST_CASE("ZeroAngularRate_GivesConstantPositionAndZeroVelocity", "[ephemeris]") {
    BodySystemBuilder b;

    BodyDef star;
    star.id = 0;
    star.parent_id = InvalidBody;
    star.mu = 1.0e10;
    star.radius = 100.0;
    star.soi_radius = 1e9;

    BodyDef planet;
    planet.id = 1;
    planet.parent_id = 0;
    planet.mu = 1.0;
    planet.radius = 1.0;
    planet.soi_radius = 2.0;
    planet.orbit_radius = 50.0;
    planet.angular_rate = 0.0;
    planet.phase_at_epoch = 0.5;

    b.add_body(star);
    b.add_body(planet);

    BodySystem sys;
    REQUIRE(b.build(sys) == SolveStatus::Ok);

    double expected_x = 50.0 * std::cos(0.5);
    double expected_y = 50.0 * std::sin(0.5);

    for (double t : {0.0, 1.0, 100.0, -50.0}) {
        Vec2 pos = sys.position_in_parent(1, t);
        REQUIRE_THAT(pos.x, WithinAbs(expected_x, 1e-12));
        REQUIRE_THAT(pos.y, WithinAbs(expected_y, 1e-12));

        Vec2 vel = sys.velocity_in_parent(1, t);
        REQUIRE_THAT(vel.x, WithinAbs(0.0, 1e-15));
        REQUIRE_THAT(vel.y, WithinAbs(0.0, 1e-15));
    }
}
