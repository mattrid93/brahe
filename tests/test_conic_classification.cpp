#include <catch2/catch_test_macros.hpp>

#include "brahe/two_body.h"

#include <cmath>
#include <limits>

using namespace brahe;

TEST_CASE("Classify_CircularOrbit_AsEllipse", "[phase2][classify]") {
    double mu = 1.0;
    State2 s{{1, 0}, {0, 1}};
    ConicType ct;
    REQUIRE(TwoBody::classify(mu, s, ct) == SolveStatus::Ok);
    REQUIRE(ct == ConicType::Ellipse);
}

TEST_CASE("Classify_EllipticOrbit_AsEllipse", "[phase2][classify]") {
    double mu = 1.0;
    // e = 0.5 at periapsis
    double a = 2.0, e = 0.5;
    double vp = std::sqrt(mu * (1.0 + e) / (a * (1.0 - e)));
    State2 s{{1, 0}, {0, vp}};
    ConicType ct;
    REQUIRE(TwoBody::classify(mu, s, ct) == SolveStatus::Ok);
    REQUIRE(ct == ConicType::Ellipse);
}

TEST_CASE("Classify_HyperbolicOrbit_AsHyperbola", "[phase2][classify]") {
    double mu = 1.0;
    State2 s{{1, 0}, {0, 2}};
    ConicType ct;
    REQUIRE(TwoBody::classify(mu, s, ct) == SolveStatus::Ok);
    REQUIRE(ct == ConicType::Hyperbola);
}

TEST_CASE("Classify_EccentricityWithinBandOfOne_AsParabolaLike", "[phase2][classify]") {
    // Construct a state with e very close to 1.
    // At periapsis with r=(rp,0), v=(0,vp):
    //   e = rp*vp^2/mu - 1
    // For e=1 exactly: vp = sqrt(2*mu/rp)
    // Use rp=1, mu=1 => vp = sqrt(2) gives e=1 exactly
    double mu = 1.0;
    double vp = std::sqrt(2.0); // e = 1.0 exactly
    State2 s{{1, 0}, {0, vp}};
    ConicType ct;
    REQUIRE(TwoBody::classify(mu, s, ct) == SolveStatus::Ok);
    REQUIRE(ct == ConicType::ParabolaLike);
}

TEST_CASE("Classify_EccentricityJustBelowBand_AsEllipse", "[phase2][classify]") {
    // Need e = 1 - band - small_offset to be clearly elliptic
    // band = 1e-6, so target e ~ 0.999998
    // At periapsis: e = rp*vp^2/mu - 1
    // vp^2 = mu*(1+e)/rp => vp = sqrt(1 + 0.999998) = sqrt(1.999998)
    double mu = 1.0;
    double e_target = 1.0 - 1.5e-6; // clearly below band
    double vp = std::sqrt(mu * (1.0 + e_target));
    State2 s{{1, 0}, {0, vp}};
    ConicType ct;
    REQUIRE(TwoBody::classify(mu, s, ct) == SolveStatus::Ok);
    REQUIRE(ct == ConicType::Ellipse);
}

TEST_CASE("Classify_EccentricityJustAboveBand_AsHyperbola", "[phase2][classify]") {
    double mu = 1.0;
    double e_target = 1.0 + 1.5e-6; // clearly above band
    double vp = std::sqrt(mu * (1.0 + e_target));
    State2 s{{1, 0}, {0, vp}};
    ConicType ct;
    REQUIRE(TwoBody::classify(mu, s, ct) == SolveStatus::Ok);
    REQUIRE(ct == ConicType::Hyperbola);
}

TEST_CASE("Classify_InvalidMu_Fails", "[phase2][classify]") {
    State2 s{{1, 0}, {0, 1}};
    ConicType ct;
    REQUIRE(TwoBody::classify(0.0, s, ct) == SolveStatus::InvalidInput);
    REQUIRE(TwoBody::classify(-1.0, s, ct) == SolveStatus::InvalidInput);
}

TEST_CASE("Classify_ZeroRadius_Fails", "[phase2][classify]") {
    double mu = 1.0;
    State2 s{{0, 0}, {0, 1}};
    ConicType ct;
    REQUIRE(TwoBody::classify(mu, s, ct) == SolveStatus::InvalidInput);
}

TEST_CASE("Classify_NonFiniteInput_Fails", "[phase2][classify]") {
    double mu = 1.0;
    double inf = std::numeric_limits<double>::infinity();
    double nan = std::numeric_limits<double>::quiet_NaN();

    ConicType ct;
    REQUIRE(TwoBody::classify(mu, State2{{inf, 0}, {0, 1}}, ct) == SolveStatus::InvalidInput);
    REQUIRE(TwoBody::classify(mu, State2{{nan, 0}, {0, 1}}, ct) == SolveStatus::InvalidInput);
    REQUIRE(TwoBody::classify(mu, State2{{1, 0}, {0, inf}}, ct) == SolveStatus::InvalidInput);
}
