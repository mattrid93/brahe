#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "brahe/two_body.h"
#include "brahe/vec2.h"

#include <cmath>

using namespace brahe;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

// --- F. Propagation tests ---

// F1. Basic propagation tests

TEST_CASE("Propagate_CircularOrbit_ByZeroDt_ReturnsSameState", "[phase2][propagate]") {
    double mu = 1.0;
    State2 s{{1, 0}, {0, 1}};
    State2 out;
    REQUIRE(TwoBody::propagate(mu, s, 0.0, out) == SolveStatus::Ok);
    REQUIRE_THAT(out.r.x, WithinAbs(1.0, 1e-12));
    REQUIRE_THAT(out.r.y, WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(out.v.x, WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(out.v.y, WithinAbs(1.0, 1e-12));
}

TEST_CASE("Propagate_CircularOrbit_ByQuarterPeriod_RotatesStateCorrectly",
          "[phase2][propagate]") {
    // mu=1, a=1 => period = 2*pi, quarter = pi/2
    double mu = 1.0;
    State2 s{{1, 0}, {0, 1}};
    double dt = M_PI / 2.0;
    State2 out;
    REQUIRE(TwoBody::propagate(mu, s, dt, out) == SolveStatus::Ok);
    REQUIRE_THAT(out.r.x, WithinAbs(0.0, 1e-9));
    REQUIRE_THAT(out.r.y, WithinAbs(1.0, 1e-9));
    REQUIRE_THAT(out.v.x, WithinAbs(-1.0, 1e-9));
    REQUIRE_THAT(out.v.y, WithinAbs(0.0, 1e-9));
}

TEST_CASE("Propagate_CircularOrbit_ByHalfPeriod_RotatesStateCorrectly", "[phase2][propagate]") {
    double mu = 1.0;
    State2 s{{1, 0}, {0, 1}};
    double dt = M_PI;
    State2 out;
    REQUIRE(TwoBody::propagate(mu, s, dt, out) == SolveStatus::Ok);
    REQUIRE_THAT(out.r.x, WithinAbs(-1.0, 1e-9));
    REQUIRE_THAT(out.r.y, WithinAbs(0.0, 1e-9));
    REQUIRE_THAT(out.v.x, WithinAbs(0.0, 1e-9));
    REQUIRE_THAT(out.v.y, WithinAbs(-1.0, 1e-9));
}

TEST_CASE("Propagate_EllipticOrbit_OneFullPeriod_ReturnsOriginalState", "[phase2][propagate]") {
    double mu = 1.0;
    double a = 2.0, e = 0.5;
    double vp = std::sqrt(mu * (1.0 + e) / (a * (1.0 - e)));
    State2 s{{1, 0}, {0, vp}};
    double period = 2.0 * M_PI * std::sqrt(a * a * a / mu);

    State2 out;
    REQUIRE(TwoBody::propagate(mu, s, period, out) == SolveStatus::Ok);
    REQUIRE_THAT(out.r.x, WithinAbs(s.r.x, 1e-8));
    REQUIRE_THAT(out.r.y, WithinAbs(s.r.y, 1e-8));
    REQUIRE_THAT(out.v.x, WithinAbs(s.v.x, 1e-8));
    REQUIRE_THAT(out.v.y, WithinAbs(s.v.y, 1e-8));
}

TEST_CASE("Propagate_EllipticOrbit_TwoStepMatchesSingleStep", "[phase2][propagate]") {
    double mu = 1.0;
    double a = 2.0, e = 0.5;
    double vp = std::sqrt(mu * (1.0 + e) / (a * (1.0 - e)));
    State2 s{{1, 0}, {0, vp}};

    double dt1 = 1.5, dt2 = 2.3;

    // Two-step
    State2 mid, two_step;
    REQUIRE(TwoBody::propagate(mu, s, dt1, mid) == SolveStatus::Ok);
    REQUIRE(TwoBody::propagate(mu, mid, dt2, two_step) == SolveStatus::Ok);

    // Single step
    State2 single_step;
    REQUIRE(TwoBody::propagate(mu, s, dt1 + dt2, single_step) == SolveStatus::Ok);

    REQUIRE_THAT(two_step.r.x, WithinAbs(single_step.r.x, 1e-8));
    REQUIRE_THAT(two_step.r.y, WithinAbs(single_step.r.y, 1e-8));
    REQUIRE_THAT(two_step.v.x, WithinAbs(single_step.v.x, 1e-8));
    REQUIRE_THAT(two_step.v.y, WithinAbs(single_step.v.y, 1e-8));
}

TEST_CASE("Propagate_HyperbolicOrbit_PreservesEnergyAndAngularMomentum", "[phase2][propagate]") {
    double mu = 1.0;
    State2 s{{1, 0}, {0, 2}};

    double eps0 = TwoBody::specific_energy(mu, s);
    double h0 = TwoBody::specific_angular_momentum_z(s);

    State2 out;
    REQUIRE(TwoBody::propagate(mu, s, 5.0, out) == SolveStatus::Ok);

    double eps1 = TwoBody::specific_energy(mu, out);
    double h1 = TwoBody::specific_angular_momentum_z(out);

    REQUIRE_THAT(eps1, WithinAbs(eps0, 1e-10));
    REQUIRE_THAT(h1, WithinAbs(h0, 1e-10));
}

TEST_CASE("Propagate_ParabolaLikeCase_RemainsFinite", "[phase2][propagate]") {
    double mu = 1.0;
    double vp = std::sqrt(2.0);
    State2 s{{1, 0}, {0, vp}};

    State2 out;
    REQUIRE(TwoBody::propagate(mu, s, 1.0, out) == SolveStatus::Ok);
    REQUIRE(std::isfinite(out.r.x));
    REQUIRE(std::isfinite(out.r.y));
    REQUIRE(std::isfinite(out.v.x));
    REQUIRE(std::isfinite(out.v.y));

    // Energy should be approximately zero (parabolic)
    double eps = TwoBody::specific_energy(mu, out);
    REQUIRE_THAT(eps, WithinAbs(0.0, 1e-8));
}

// F2. High-eccentricity and near-degenerate tests

TEST_CASE("Propagate_HighEccentricityEllipse_DoesNotExplodeNearPeriapsis", "[phase2][propagate]") {
    double mu = 1.0;
    double a = 100.0, e = 0.999;
    double rp = a * (1.0 - e); // 0.1
    double vp = std::sqrt(mu * (1.0 + e) / (a * (1.0 - e)));
    State2 s{{rp, 0}, {0, vp}};

    double period = 2.0 * M_PI * std::sqrt(a * a * a / mu);

    State2 out;
    REQUIRE(TwoBody::propagate(mu, s, period, out) == SolveStatus::Ok);
    REQUIRE_THAT(out.r.x, WithinAbs(s.r.x, 1e-4));
    REQUIRE_THAT(out.r.y, WithinAbs(s.r.y, 1e-4));
}

TEST_CASE("Propagate_NearParabolicBelowBand_RemainsStable", "[phase2][propagate]") {
    double mu = 1.0;
    double e_target = 1.0 - 1.5e-6; // just below band -> Ellipse
    double vp = std::sqrt(mu * (1.0 + e_target));
    State2 s{{1, 0}, {0, vp}};

    State2 out;
    REQUIRE(TwoBody::propagate(mu, s, 10.0, out) == SolveStatus::Ok);
    REQUIRE(std::isfinite(out.r.x));
    REQUIRE(std::isfinite(out.r.y));

    double eps0 = TwoBody::specific_energy(mu, s);
    double eps1 = TwoBody::specific_energy(mu, out);
    REQUIRE_THAT(eps1, WithinAbs(eps0, 1e-6));
}

TEST_CASE("Propagate_NearParabolicAboveBand_RemainsStable", "[phase2][propagate]") {
    double mu = 1.0;
    double e_target = 1.0 + 1.5e-6; // just above band -> Hyperbola
    double vp = std::sqrt(mu * (1.0 + e_target));
    State2 s{{1, 0}, {0, vp}};

    State2 out;
    REQUIRE(TwoBody::propagate(mu, s, 10.0, out) == SolveStatus::Ok);
    REQUIRE(std::isfinite(out.r.x));
    REQUIRE(std::isfinite(out.r.y));

    double eps0 = TwoBody::specific_energy(mu, s);
    double eps1 = TwoBody::specific_energy(mu, out);
    REQUIRE_THAT(eps1, WithinAbs(eps0, 1e-6));
}

TEST_CASE("Propagate_VerySmallDt_ProducesConsistentState", "[phase2][propagate]") {
    double mu = 1.0;
    State2 s{{1, 0}, {0, 1}};

    State2 out;
    REQUIRE(TwoBody::propagate(mu, s, 1e-10, out) == SolveStatus::Ok);

    // Should barely move from initial state
    REQUIRE_THAT(out.r.x, WithinAbs(s.r.x, 1e-6));
    REQUIRE_THAT(out.r.y, WithinAbs(s.r.y, 1e-6));
    REQUIRE_THAT(out.v.x, WithinAbs(s.v.x, 1e-6));
    REQUIRE_THAT(out.v.y, WithinAbs(s.v.y, 1e-6));
}

TEST_CASE("Propagate_LargeDtEllipse_MultiplePeriodsStillConsistent", "[phase2][propagate]") {
    double mu = 1.0;
    double a = 2.0, e = 0.3;
    double rp = a * (1.0 - e);
    double vp = std::sqrt(mu * (1.0 + e) / (a * (1.0 - e)));
    State2 s{{rp, 0}, {0, vp}};

    double period = 2.0 * M_PI * std::sqrt(a * a * a / mu);
    double dt = 10.0 * period; // 10 full orbits

    State2 out;
    REQUIRE(TwoBody::propagate(mu, s, dt, out) == SolveStatus::Ok);

    // After exactly 10 periods, should return to start
    REQUIRE_THAT(out.r.x, WithinAbs(s.r.x, 1e-5));
    REQUIRE_THAT(out.r.y, WithinAbs(s.r.y, 1e-5));
}

// F3. Directionality and physical sanity tests

TEST_CASE("Propagate_ForwardTime_DoesNotReverseAngularMomentumSign", "[phase2][propagate]") {
    double mu = 1.0;
    State2 s{{1, 0}, {0, 1}};
    double h0 = TwoBody::specific_angular_momentum_z(s);
    REQUIRE(h0 > 0.0);

    State2 out;
    REQUIRE(TwoBody::propagate(mu, s, 2.0, out) == SolveStatus::Ok);
    double h1 = TwoBody::specific_angular_momentum_z(out);
    REQUIRE_THAT(h1, WithinAbs(h0, 1e-10));
}

TEST_CASE("Propagate_PericenterState_MovesAwayFromPericenterAfterPositiveDt",
          "[phase2][propagate]") {
    double mu = 1.0;
    double a = 2.0, e = 0.5;
    double rp = a * (1.0 - e);
    double vp = std::sqrt(mu * (1.0 + e) / (a * (1.0 - e)));
    State2 s{{rp, 0}, {0, vp}};

    State2 out;
    REQUIRE(TwoBody::propagate(mu, s, 0.5, out) == SolveStatus::Ok);

    double r_after = length(out.r);
    REQUIRE(r_after > rp + 1e-6);
}

TEST_CASE("Propagate_HyperbolicOutboundState_RemainsOutboundForSufficientlyLargeDt",
          "[phase2][propagate]") {
    double mu = 1.0;
    // Start outbound on hyperbola
    State2 s{{2, 1}, {0.5, 1.5}};
    double r0 = length(s.r);

    State2 out;
    REQUIRE(TwoBody::propagate(mu, s, 100.0, out) == SolveStatus::Ok);
    double r1 = length(out.r);
    REQUIRE(r1 > r0);
}
