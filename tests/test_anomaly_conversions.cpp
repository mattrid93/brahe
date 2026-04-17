#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "brahe/two_body.h"

#include <cmath>
#include <limits>

using namespace brahe;
using namespace brahe::detail;
using Catch::Matchers::WithinAbs;

static constexpr double kTol = 1e-12;

// --- E. Anomaly conversion tests ---

// Ellipse tests

TEST_CASE("TrueToEccentricToTrue_RoundTrip_Ellipse", "[phase2][anomaly]") {
    double e = 0.5;
    for (double nu : {0.0, M_PI / 6, M_PI / 4, M_PI / 2, M_PI * 0.9, -M_PI / 3, -M_PI + 0.01}) {
        double E = true_to_eccentric_anomaly(nu, e);
        double nu_back = eccentric_to_true_anomaly(E, e);
        REQUIRE_THAT(nu_back, WithinAbs(nu, kTol));
    }
}

TEST_CASE("MeanToEccentricToMean_RoundTrip_Ellipse", "[phase2][anomaly]") {
    double e = 0.6;
    for (double M : {0.0, 0.5, 1.0, M_PI / 2, M_PI, -1.0, -M_PI + 0.01}) {
        double E;
        REQUIRE(mean_to_eccentric_anomaly(M, e, E) == SolveStatus::Ok);
        double M_back = eccentric_to_mean_anomaly(E, e);
        REQUIRE_THAT(M_back, WithinAbs(M, 1e-10));
    }
}

TEST_CASE("MeanToEccentric_ConvergesForHighEccentricitySmallMeanAnomaly", "[phase2][anomaly]") {
    double e = 0.9937767804941645;
    double M = 0.05128749566606613;

    double E;
    REQUIRE(mean_to_eccentric_anomaly(M, e, E) == SolveStatus::Ok);
    double M_back = eccentric_to_mean_anomaly(E, e);
    REQUIRE_THAT(M_back, WithinAbs(M, 1e-10));
}

TEST_CASE("EccentricToMean_MonotonicOverPrincipalRange", "[phase2][anomaly]") {
    double e = 0.3;
    double prev_M = eccentric_to_mean_anomaly(-M_PI, e);
    int steps = 100;
    for (int i = 1; i <= steps; ++i) {
        double E = -M_PI + 2.0 * M_PI * i / steps;
        double M = eccentric_to_mean_anomaly(E, e);
        REQUIRE(M > prev_M);
        prev_M = M;
    }
}

// Hyperbola tests

TEST_CASE("TrueToHyperbolicToTrue_RoundTrip_Hyperbola", "[phase2][anomaly]") {
    double e = 2.0;
    double nu_max = std::acos(-1.0 / e) - 0.01; // just inside asymptote
    for (double nu : {0.0, 0.5, 1.0, nu_max, -0.5, -1.0, -nu_max}) {
        double H = true_to_hyperbolic_anomaly(nu, e);
        double nu_back = hyperbolic_to_true_anomaly(H, e);
        REQUIRE_THAT(nu_back, WithinAbs(nu, kTol));
    }
}

TEST_CASE("MeanToHyperbolicToMean_RoundTrip_Hyperbola", "[phase2][anomaly]") {
    double e = 1.5;
    for (double Mh : {0.0, 0.5, 2.0, 10.0, -0.5, -2.0, -10.0}) {
        double H;
        REQUIRE(mean_to_hyperbolic_anomaly(Mh, e, H) == SolveStatus::Ok);
        double Mh_back = hyperbolic_to_mean_anomaly(H, e);
        REQUIRE_THAT(Mh_back, WithinAbs(Mh, 1e-10));
    }
}

TEST_CASE("MeanToHyperbolic_ConvergesNearParabolicUnitMeanAnomaly", "[phase2][anomaly]") {
    double e = 1.0196305827783696;
    double Mh = 0.9999999596267415;

    double H;
    REQUIRE(mean_to_hyperbolic_anomaly(Mh, e, H) == SolveStatus::Ok);
    double Mh_back = hyperbolic_to_mean_anomaly(H, e);
    REQUIRE_THAT(Mh_back, WithinAbs(Mh, 1e-10));
}

TEST_CASE("HyperbolicConversions_RemainFiniteNearAsymptote", "[phase2][anomaly]") {
    double e = 1.5;
    double nu_max = std::acos(-1.0 / e);
    // Approach asymptote closely but don't hit it
    double nu_near = nu_max - 1e-6;
    double H = true_to_hyperbolic_anomaly(nu_near, e);
    REQUIRE(std::isfinite(H));

    double Mh = hyperbolic_to_mean_anomaly(H, e);
    REQUIRE(std::isfinite(Mh));

    double H_back;
    REQUIRE(mean_to_hyperbolic_anomaly(Mh, e, H_back) == SolveStatus::Ok);
    REQUIRE(std::isfinite(H_back));
}

// Guardrail tests

TEST_CASE("InverseTrigInputs_AreClamped_NotNaN", "[phase2][anomaly]") {
    // safe_acos should never return NaN
    REQUIRE(std::isfinite(safe_acos(1.0)));
    REQUIRE(std::isfinite(safe_acos(-1.0)));
    REQUIRE(std::isfinite(safe_acos(1.0 + 1e-15)));  // slightly out of range
    REQUIRE(std::isfinite(safe_acos(-1.0 - 1e-15)));
    REQUIRE(std::isfinite(safe_acos(2.0)));            // well out of range

    // safe_asin similarly
    REQUIRE(std::isfinite(safe_asin(1.0)));
    REQUIRE(std::isfinite(safe_asin(-1.0)));
    REQUIRE(std::isfinite(safe_asin(1.0 + 1e-15)));
    REQUIRE(std::isfinite(safe_asin(-1.0 - 1e-15)));
}

TEST_CASE("TinyNegativeDiscriminant_IsTreatedAsZero", "[phase2][anomaly]") {
    // For a circular orbit (e ≈ 0), computing eccentric anomaly from true anomaly
    // should not produce NaN even when intermediate discriminants are slightly negative
    double e = 1e-14;
    double nu = 1.0;
    double E = true_to_eccentric_anomaly(nu, e);
    REQUIRE(std::isfinite(E));

    double nu_back = eccentric_to_true_anomaly(E, e);
    REQUIRE(std::isfinite(nu_back));
    REQUIRE_THAT(nu_back, WithinAbs(nu, 1e-8));
}
