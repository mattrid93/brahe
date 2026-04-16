#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "brahe/two_body.h"
#include "brahe/vec2.h"

#include <cmath>

using namespace brahe;
using Catch::Matchers::WithinAbs;

// --- Radial (h = 0) propagation tests ---
//
// Purely radial trajectories have angular momentum h = r x v = 0, which makes
// the conic-elements parameterisation r(nu) = p / (1 + e cos(nu)) degenerate:
// p = 0, e collapses to 1 regardless of energy, and the classic Kepler
// propagation pipeline (to_elements -> Kepler solve -> from_elements) breaks.
// These tests pin down the correct behaviour so that a radial Kepler solver
// can be added without regressing the normal code paths.

TEST_CASE("Propagate_RadialOutbound_BoundCase_RemainsFiniteAndOnAxis",
          "[two_body][radial]") {
    // mu=1, r=(10, 0), v=(0.3, 0) -> eps = 0.045 - 0.1 = -0.055 (bound radial).
    // Turning point at r = -mu/eps ~= 18.18. A short propagation must stay
    // strictly on the +x axis (y stays 0) with r between 10 and 18.18.
    double mu = 1.0;
    State2 s{{10.0, 0.0}, {0.3, 0.0}};

    State2 out;
    REQUIRE(TwoBody::propagate(mu, s, 1.0, out) == SolveStatus::Ok);
    REQUIRE(std::isfinite(out.r.x));
    REQUIRE(std::isfinite(out.r.y));
    REQUIRE(std::isfinite(out.v.x));
    REQUIRE(std::isfinite(out.v.y));
    REQUIRE_THAT(out.r.y, WithinAbs(0.0, 1e-9));
    REQUIRE_THAT(out.v.y, WithinAbs(0.0, 1e-9));
    REQUIRE(out.r.x > 10.0);       // moved outward
    REQUIRE(out.r.x < 18.19);      // did not cross turning point
}

TEST_CASE("Propagate_RadialOutbound_BoundCase_ConservesEnergy",
          "[two_body][radial]") {
    double mu = 1.0;
    State2 s{{10.0, 0.0}, {0.3, 0.0}};
    double eps0 = TwoBody::specific_energy(mu, s);

    State2 out;
    REQUIRE(TwoBody::propagate(mu, s, 1.0, out) == SolveStatus::Ok);
    double eps1 = TwoBody::specific_energy(mu, out);
    REQUIRE_THAT(eps1, WithinAbs(eps0, 1e-10));
}

TEST_CASE("Propagate_RadialOutbound_AngularMomentumStaysZero",
          "[two_body][radial]") {
    double mu = 1.0;
    State2 s{{10.0, 0.0}, {0.3, 0.0}};

    State2 out;
    REQUIRE(TwoBody::propagate(mu, s, 1.0, out) == SolveStatus::Ok);
    double h1 = TwoBody::specific_angular_momentum_z(out);
    REQUIRE_THAT(h1, WithinAbs(0.0, 1e-12));
}

TEST_CASE("Propagate_RadialParabolic_MatchesClosedFormSolution",
          "[two_body][radial]") {
    // Outbound radial parabolic: v = sqrt(2*mu/r_0), eps = 0.
    // Closed form: r(t) = (r_0^{3/2} + (3/2) sqrt(2 mu) t)^{2/3}.
    // At mu=1, r_0=1, t=1: r = (1 + 1.5*sqrt(2))^{2/3} = 3.12132^{2/3} ~= 2.13436.
    double mu = 1.0;
    double r_0 = 1.0;
    double v_0 = std::sqrt(2.0 * mu / r_0);
    State2 s{{r_0, 0.0}, {v_0, 0.0}};

    double dt = 1.0;
    double r_expected_3half = std::pow(r_0, 1.5) + 1.5 * std::sqrt(2.0 * mu) * dt;
    double r_expected = std::pow(r_expected_3half, 2.0 / 3.0);
    double v_expected = std::sqrt(2.0 * mu / r_expected);

    State2 out;
    REQUIRE(TwoBody::propagate(mu, s, dt, out) == SolveStatus::Ok);
    REQUIRE_THAT(out.r.x, WithinAbs(r_expected, 1e-9));
    REQUIRE_THAT(out.r.y, WithinAbs(0.0, 1e-9));
    REQUIRE_THAT(out.v.x, WithinAbs(v_expected, 1e-9));
    REQUIRE_THAT(out.v.y, WithinAbs(0.0, 1e-9));
}

TEST_CASE("Propagate_RadialHyperbolic_OutboundGrowsMonotonically",
          "[two_body][radial]") {
    // Radial hyperbolic escape: v > sqrt(2*mu/r_0) means eps > 0.
    // r must increase monotonically; energy and zero angular momentum preserved.
    double mu = 1.0;
    State2 s{{1.0, 0.0}, {2.0, 0.0}};  // eps = 2 - 1 = 1 > 0
    double eps0 = TwoBody::specific_energy(mu, s);

    State2 out1, out2;
    REQUIRE(TwoBody::propagate(mu, s, 1.0, out1) == SolveStatus::Ok);
    REQUIRE(TwoBody::propagate(mu, s, 2.0, out2) == SolveStatus::Ok);

    REQUIRE(out1.r.x > s.r.x);
    REQUIRE(out2.r.x > out1.r.x);
    REQUIRE_THAT(out1.r.y, WithinAbs(0.0, 1e-9));
    REQUIRE_THAT(out2.r.y, WithinAbs(0.0, 1e-9));

    double eps1 = TwoBody::specific_energy(mu, out1);
    double eps2 = TwoBody::specific_energy(mu, out2);
    REQUIRE_THAT(eps1, WithinAbs(eps0, 1e-10));
    REQUIRE_THAT(eps2, WithinAbs(eps0, 1e-10));
}

TEST_CASE("Propagate_RadialInbound_ShortDt_MovesTowardOrigin",
          "[two_body][radial]") {
    // Inbound radial, short dt so we remain well above the singularity at r=0.
    double mu = 1.0;
    State2 s{{10.0, 0.0}, {-0.2, 0.0}};
    double eps0 = TwoBody::specific_energy(mu, s);

    State2 out;
    REQUIRE(TwoBody::propagate(mu, s, 0.5, out) == SolveStatus::Ok);
    REQUIRE(std::isfinite(out.r.x));
    REQUIRE_THAT(out.r.y, WithinAbs(0.0, 1e-9));
    REQUIRE(out.r.x < s.r.x);   // moved inward
    REQUIRE(out.r.x > 0.0);      // did not pass origin

    double eps1 = TwoBody::specific_energy(mu, out);
    REQUIRE_THAT(eps1, WithinAbs(eps0, 1e-10));
}

TEST_CASE("Propagate_RadialOffAxis_StateStaysOnRadialLine",
          "[two_body][radial]") {
    // Radial along a non-axis direction: r parallel to v but not on the x-axis.
    // Position and velocity should remain parallel to the initial radial line.
    double mu = 1.0;
    double theta = 0.7;  // arbitrary direction
    Vec2 dir = {std::cos(theta), std::sin(theta)};
    State2 s{{5.0 * dir.x, 5.0 * dir.y}, {0.4 * dir.x, 0.4 * dir.y}};

    State2 out;
    REQUIRE(TwoBody::propagate(mu, s, 1.5, out) == SolveStatus::Ok);
    REQUIRE(std::isfinite(out.r.x));
    REQUIRE(std::isfinite(out.r.y));

    // out.r must be along dir (cross product with dir is zero)
    double cross = out.r.x * dir.y - out.r.y * dir.x;
    REQUIRE_THAT(cross, WithinAbs(0.0, 1e-9));
    double v_cross = out.v.x * dir.y - out.v.y * dir.x;
    REQUIRE_THAT(v_cross, WithinAbs(0.0, 1e-9));
}

TEST_CASE("Propagate_RadialOutbound_ZeroDt_ReturnsSameState",
          "[two_body][radial]") {
    // Radial case with dt == 0 must still return the initial state exactly.
    double mu = 1.0;
    State2 s{{7.5, 0.0}, {0.2, 0.0}};

    State2 out;
    REQUIRE(TwoBody::propagate(mu, s, 0.0, out) == SolveStatus::Ok);
    REQUIRE(out.r.x == s.r.x);
    REQUIRE(out.r.y == s.r.y);
    REQUIRE(out.v.x == s.v.x);
    REQUIRE(out.v.y == s.v.y);
}
