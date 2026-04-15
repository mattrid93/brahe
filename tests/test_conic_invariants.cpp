#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "brahe/two_body.h"
#include "brahe/vec2.h"

#include <cmath>
#include <limits>
#include <type_traits>

using namespace brahe;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

// --- A. Header and API compile tests ---

TEST_CASE("TwoBody_HeaderCompiles", "[phase2][api]") {
    // Instantiate types in evaluated context to verify layout
    State2 s{};
    ConicElements2D e{};
    ConicState cs{};
    (void)s;
    (void)e;
    (void)cs;
    SUCCEED();
}

TEST_CASE("ConicElements2D_IsTriviallyCopyable", "[phase2][api]") {
    static_assert(std::is_trivially_copyable_v<ConicElements2D>);
    SUCCEED();
}

TEST_CASE("ConicElements2D_IsStandardLayout", "[phase2][api]") {
    static_assert(std::is_standard_layout_v<ConicElements2D>);
    SUCCEED();
}

TEST_CASE("TwoBody_PublicFunctions_AreStaticAndCallable", "[phase2][api]") {
    // Verify all public static functions exist with expected signatures
    double mu = 1.0;
    State2 s{{1, 0}, {0, 1}};
    ConicType ct;
    ConicElements2D el;
    State2 out;

    (void)TwoBody::classify(mu, s, ct);
    (void)TwoBody::to_elements(mu, s, el);
    (void)TwoBody::from_elements(el, out);
    (void)TwoBody::propagate(mu, s, 0.0, out);
    (void)TwoBody::specific_energy(mu, s);
    (void)TwoBody::specific_angular_momentum_z(s);
    (void)TwoBody::eccentricity(mu, s);
    (void)TwoBody::semi_major_axis(mu, s);
    (void)TwoBody::semi_latus_rectum(mu, s);
    (void)TwoBody::periapsis_radius(mu, s);
    (void)TwoBody::apoapsis_radius(mu, s);
    SUCCEED();
}

// --- B. Invariant tests ---

// Circular unit orbit: mu=1, r=(1,0), v=(0,1)
// eps=-0.5, h=1, e=0, a=1, p=1, rp=1, ra=1

TEST_CASE("SpecificEnergy_CircularUnitOrbit_IsMinusHalf", "[phase2][invariants]") {
    double mu = 1.0;
    State2 s{{1, 0}, {0, 1}};
    REQUIRE_THAT(TwoBody::specific_energy(mu, s), WithinAbs(-0.5, 1e-12));
}

TEST_CASE("AngularMomentumZ_CircularUnitOrbit_IsOne", "[phase2][invariants]") {
    State2 s{{1, 0}, {0, 1}};
    REQUIRE_THAT(TwoBody::specific_angular_momentum_z(s), WithinAbs(1.0, 1e-12));
}

TEST_CASE("Eccentricity_CircularUnitOrbit_IsZero", "[phase2][invariants]") {
    double mu = 1.0;
    State2 s{{1, 0}, {0, 1}};
    REQUIRE_THAT(TwoBody::eccentricity(mu, s), WithinAbs(0.0, 1e-10));
}

TEST_CASE("SemiMajorAxis_CircularUnitOrbit_IsOne", "[phase2][invariants]") {
    double mu = 1.0;
    State2 s{{1, 0}, {0, 1}};
    REQUIRE_THAT(TwoBody::semi_major_axis(mu, s), WithinAbs(1.0, 1e-12));
}

TEST_CASE("SemiLatusRectum_CircularUnitOrbit_IsOne", "[phase2][invariants]") {
    double mu = 1.0;
    State2 s{{1, 0}, {0, 1}};
    REQUIRE_THAT(TwoBody::semi_latus_rectum(mu, s), WithinAbs(1.0, 1e-12));
}

TEST_CASE("Periapsis_CircularUnitOrbit_IsOne", "[phase2][invariants]") {
    double mu = 1.0;
    State2 s{{1, 0}, {0, 1}};
    REQUIRE_THAT(TwoBody::periapsis_radius(mu, s), WithinAbs(1.0, 1e-12));
}

TEST_CASE("Apoapsis_CircularUnitOrbit_IsOne", "[phase2][invariants]") {
    double mu = 1.0;
    State2 s{{1, 0}, {0, 1}};
    REQUIRE_THAT(TwoBody::apoapsis_radius(mu, s), WithinAbs(1.0, 1e-12));
}

// Elliptic: mu=1, a=2, e=0.5, rp=1, vp=sqrt(1.5), state at periapsis
TEST_CASE("Invariants_AtEllipticPeriapsis_MatchExpectedValues", "[phase2][invariants]") {
    double mu = 1.0;
    double a = 2.0;
    double e = 0.5;
    double rp = a * (1.0 - e); // 1.0
    double vp = std::sqrt(mu * (1.0 + e) / (a * (1.0 - e))); // sqrt(1.5)
    State2 s{{rp, 0}, {0, vp}};

    REQUIRE_THAT(TwoBody::specific_energy(mu, s), WithinAbs(-mu / (2.0 * a), 1e-12));
    REQUIRE_THAT(TwoBody::specific_angular_momentum_z(s), WithinAbs(rp * vp, 1e-12));
    REQUIRE_THAT(TwoBody::eccentricity(mu, s), WithinAbs(e, 1e-10));
    REQUIRE_THAT(TwoBody::semi_major_axis(mu, s), WithinAbs(a, 1e-10));
    REQUIRE_THAT(TwoBody::periapsis_radius(mu, s), WithinAbs(rp, 1e-10));

    double ra = a * (1.0 + e); // 3.0
    REQUIRE_THAT(TwoBody::apoapsis_radius(mu, s), WithinAbs(ra, 1e-10));
}

// Hyperbolic: mu=1, r=(1,0), v=(0,2)
// eps=1, h=2, p=4, e=3, rp=1, ra=infinity
TEST_CASE("Invariants_HyperbolicCase_MatchExpectedValues", "[phase2][invariants]") {
    double mu = 1.0;
    State2 s{{1, 0}, {0, 2}};

    REQUIRE_THAT(TwoBody::specific_energy(mu, s), WithinAbs(1.0, 1e-12));
    REQUIRE_THAT(TwoBody::specific_angular_momentum_z(s), WithinAbs(2.0, 1e-12));
    REQUIRE_THAT(TwoBody::semi_latus_rectum(mu, s), WithinAbs(4.0, 1e-12));
    REQUIRE_THAT(TwoBody::eccentricity(mu, s), WithinAbs(3.0, 1e-10));
    REQUIRE_THAT(TwoBody::periapsis_radius(mu, s), WithinAbs(1.0, 1e-10));
    REQUIRE(TwoBody::apoapsis_radius(mu, s) == std::numeric_limits<double>::infinity());
}
