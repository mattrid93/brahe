#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "brahe/two_body.h"
#include "brahe/vec2.h"

#include <cmath>
#include <limits>

using namespace brahe;
using Catch::Matchers::WithinAbs;

static constexpr double kTol = 1e-10;

// --- D. Cartesian <-> conic conversion tests ---

TEST_CASE("ToElements_CircularOrbit_ProducesExpectedCanonicalValues", "[phase2][conversion]") {
    double mu = 1.0;
    State2 s{{1, 0}, {0, 1}};
    ConicElements2D el;
    REQUIRE(TwoBody::to_elements(mu, s, el) == SolveStatus::Ok);

    REQUIRE(el.type == ConicType::Ellipse);
    REQUIRE_THAT(el.eccentricity, WithinAbs(0.0, kTol));
    REQUIRE_THAT(el.semi_latus_rectum, WithinAbs(1.0, kTol));
    REQUIRE_THAT(el.semi_major_axis, WithinAbs(1.0, kTol));

    // Circular convention: argument_of_periapsis = 0, true_anomaly = atan2(y, x)
    REQUIRE_THAT(el.argument_of_periapsis, WithinAbs(0.0, kTol));
    REQUIRE_THAT(el.true_anomaly, WithinAbs(0.0, kTol)); // atan2(0, 1) = 0
}

TEST_CASE("ToElements_EllipticPeriapsis_ProducesExpectedValues", "[phase2][conversion]") {
    double mu = 1.0;
    double a = 2.0, e = 0.5;
    double vp = std::sqrt(mu * (1.0 + e) / (a * (1.0 - e)));
    State2 s{{1, 0}, {0, vp}};
    ConicElements2D el;
    REQUIRE(TwoBody::to_elements(mu, s, el) == SolveStatus::Ok);

    REQUIRE(el.type == ConicType::Ellipse);
    REQUIRE_THAT(el.eccentricity, WithinAbs(0.5, kTol));
    REQUIRE_THAT(el.semi_major_axis, WithinAbs(2.0, kTol));
    REQUIRE_THAT(el.true_anomaly, WithinAbs(0.0, kTol));       // at periapsis
    REQUIRE_THAT(el.mean_anomaly, WithinAbs(0.0, kTol));       // at periapsis
    REQUIRE_THAT(el.eccentric_anomaly, WithinAbs(0.0, kTol));  // at periapsis
}

TEST_CASE("ToElements_HyperbolicPeriapsis_ProducesExpectedValues", "[phase2][conversion]") {
    double mu = 1.0;
    State2 s{{1, 0}, {0, 2}};
    ConicElements2D el;
    REQUIRE(TwoBody::to_elements(mu, s, el) == SolveStatus::Ok);

    REQUIRE(el.type == ConicType::Hyperbola);
    REQUIRE_THAT(el.eccentricity, WithinAbs(3.0, kTol));
    REQUIRE_THAT(el.true_anomaly, WithinAbs(0.0, kTol));         // at periapsis
    REQUIRE_THAT(el.hyperbolic_anomaly, WithinAbs(0.0, kTol));   // at periapsis
}

TEST_CASE("FromElements_Ellipse_ReconstructsOriginalState", "[phase2][conversion]") {
    double mu = 1.0;
    double a = 2.0, e = 0.5;
    double vp = std::sqrt(mu * (1.0 + e) / (a * (1.0 - e)));
    State2 original{{1, 0}, {0, vp}};

    ConicElements2D el;
    REQUIRE(TwoBody::to_elements(mu, original, el) == SolveStatus::Ok);

    State2 reconstructed;
    REQUIRE(TwoBody::from_elements(el, reconstructed) == SolveStatus::Ok);

    REQUIRE_THAT(reconstructed.r.x, WithinAbs(original.r.x, kTol));
    REQUIRE_THAT(reconstructed.r.y, WithinAbs(original.r.y, kTol));
    REQUIRE_THAT(reconstructed.v.x, WithinAbs(original.v.x, kTol));
    REQUIRE_THAT(reconstructed.v.y, WithinAbs(original.v.y, kTol));
}

TEST_CASE("FromElements_Hyperbola_ReconstructsOriginalState", "[phase2][conversion]") {
    double mu = 1.0;
    State2 original{{1, 0}, {0, 2}};

    ConicElements2D el;
    REQUIRE(TwoBody::to_elements(mu, original, el) == SolveStatus::Ok);

    State2 reconstructed;
    REQUIRE(TwoBody::from_elements(el, reconstructed) == SolveStatus::Ok);

    REQUIRE_THAT(reconstructed.r.x, WithinAbs(original.r.x, kTol));
    REQUIRE_THAT(reconstructed.r.y, WithinAbs(original.r.y, kTol));
    REQUIRE_THAT(reconstructed.v.x, WithinAbs(original.v.x, kTol));
    REQUIRE_THAT(reconstructed.v.y, WithinAbs(original.v.y, kTol));
}

TEST_CASE("RoundTrip_CartesianToElementsToCartesian_Ellipse_MatchesOriginal",
          "[phase2][roundtrip]") {
    // Non-axis-aligned elliptic state
    double mu = 3.986e5; // Earth-like
    State2 original{{7000, 1000}, {-0.5, 7.5}};

    ConicElements2D el;
    REQUIRE(TwoBody::to_elements(mu, original, el) == SolveStatus::Ok);
    REQUIRE(el.type == ConicType::Ellipse);

    State2 reconstructed;
    REQUIRE(TwoBody::from_elements(el, reconstructed) == SolveStatus::Ok);

    REQUIRE_THAT(reconstructed.r.x, WithinAbs(original.r.x, 1e-6));
    REQUIRE_THAT(reconstructed.r.y, WithinAbs(original.r.y, 1e-6));
    REQUIRE_THAT(reconstructed.v.x, WithinAbs(original.v.x, 1e-9));
    REQUIRE_THAT(reconstructed.v.y, WithinAbs(original.v.y, 1e-9));
}

TEST_CASE("RoundTrip_CartesianToElementsToCartesian_Hyperbola_MatchesOriginal",
          "[phase2][roundtrip]") {
    double mu = 3.986e5;
    State2 original{{7000, 1000}, {-2.0, 12.0}};

    ConicElements2D el;
    REQUIRE(TwoBody::to_elements(mu, original, el) == SolveStatus::Ok);
    REQUIRE(el.type == ConicType::Hyperbola);

    State2 reconstructed;
    REQUIRE(TwoBody::from_elements(el, reconstructed) == SolveStatus::Ok);

    REQUIRE_THAT(reconstructed.r.x, WithinAbs(original.r.x, 1e-6));
    REQUIRE_THAT(reconstructed.r.y, WithinAbs(original.r.y, 1e-6));
    REQUIRE_THAT(reconstructed.v.x, WithinAbs(original.v.x, 1e-9));
    REQUIRE_THAT(reconstructed.v.y, WithinAbs(original.v.y, 1e-9));
}

TEST_CASE("RoundTrip_CartesianToElementsToCartesian_ParabolaLike_RemainsFiniteAndConsistent",
          "[phase2][roundtrip]") {
    // Near-parabolic: e ≈ 1
    double mu = 1.0;
    double vp = std::sqrt(2.0); // exactly parabolic at r=1
    State2 original{{1, 0}, {0, vp}};

    ConicElements2D el;
    REQUIRE(TwoBody::to_elements(mu, original, el) == SolveStatus::Ok);
    REQUIRE(el.type == ConicType::ParabolaLike);

    // Elements should be finite
    REQUIRE(std::isfinite(el.eccentricity));
    REQUIRE(std::isfinite(el.semi_latus_rectum));
    REQUIRE(std::isfinite(el.true_anomaly));

    State2 reconstructed;
    REQUIRE(TwoBody::from_elements(el, reconstructed) == SolveStatus::Ok);

    // Loose invariant comparison: radius and speed should match
    double r_orig = length(original.r);
    double r_recon = length(reconstructed.r);
    REQUIRE_THAT(r_recon, WithinAbs(r_orig, 1e-8));

    double v_orig = length(original.v);
    double v_recon = length(reconstructed.v);
    REQUIRE_THAT(v_recon, WithinAbs(v_orig, 1e-8));
}

TEST_CASE("FromElements_InvalidInput_FailsWithoutNaN", "[phase2][conversion]") {
    State2 out;

    // Negative mu
    ConicElements2D bad_mu{};
    bad_mu.mu = -1.0;
    bad_mu.eccentricity = 0.5;
    bad_mu.semi_latus_rectum = 1.0;
    REQUIRE(TwoBody::from_elements(bad_mu, out) == SolveStatus::InvalidInput);

    // Negative semi_latus_rectum
    ConicElements2D bad_p{};
    bad_p.mu = 1.0;
    bad_p.eccentricity = 0.5;
    bad_p.semi_latus_rectum = -1.0;
    REQUIRE(TwoBody::from_elements(bad_p, out) == SolveStatus::InvalidInput);

    // Negative eccentricity
    ConicElements2D bad_e{};
    bad_e.mu = 1.0;
    bad_e.eccentricity = -0.5;
    bad_e.semi_latus_rectum = 1.0;
    REQUIRE(TwoBody::from_elements(bad_e, out) == SolveStatus::InvalidInput);
}
