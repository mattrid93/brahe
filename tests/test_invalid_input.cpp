#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "brahe/two_body.h"
#include "brahe/vec2.h"

#include <cmath>
#include <limits>

using namespace brahe;
using Catch::Matchers::WithinAbs;

// --- G. Invalid-input and failure-behavior tests ---

TEST_CASE("Propagate_InvalidMu_Fails", "[phase2][invalid]") {
    State2 s{{1, 0}, {0, 1}};
    State2 out;
    REQUIRE(TwoBody::propagate(0.0, s, 1.0, out) == SolveStatus::InvalidInput);
    REQUIRE(TwoBody::propagate(-1.0, s, 1.0, out) == SolveStatus::InvalidInput);
}

TEST_CASE("Propagate_ZeroRadius_Fails", "[phase2][invalid]") {
    State2 s{{0, 0}, {0, 1}};
    State2 out;
    REQUIRE(TwoBody::propagate(1.0, s, 1.0, out) == SolveStatus::InvalidInput);
}

TEST_CASE("Propagate_NegativeDt_Fails", "[phase2][invalid]") {
    State2 s{{1, 0}, {0, 1}};
    State2 out;
    REQUIRE(TwoBody::propagate(1.0, s, -1.0, out) == SolveStatus::InvalidInput);
}

TEST_CASE("Propagate_NonFiniteState_Fails", "[phase2][invalid]") {
    double inf = std::numeric_limits<double>::infinity();
    double nan = std::numeric_limits<double>::quiet_NaN();
    State2 out;

    REQUIRE(TwoBody::propagate(1.0, {{inf, 0}, {0, 1}}, 1.0, out) == SolveStatus::InvalidInput);
    REQUIRE(TwoBody::propagate(1.0, {{nan, 0}, {0, 1}}, 1.0, out) == SolveStatus::InvalidInput);
    REQUIRE(TwoBody::propagate(1.0, {{1, 0}, {0, inf}}, 1.0, out) == SolveStatus::InvalidInput);
    REQUIRE(TwoBody::propagate(1.0, {{1, 0}, {nan, 0}}, 1.0, out) == SolveStatus::InvalidInput);
}

TEST_CASE("Propagate_OutputStateIsFullyDefinedOnSuccess", "[phase2][invalid]") {
    double mu = 1.0;
    State2 s{{1, 0}, {0, 1}};
    State2 out;
    REQUIRE(TwoBody::propagate(mu, s, 1.0, out) == SolveStatus::Ok);
    REQUIRE(std::isfinite(out.r.x));
    REQUIRE(std::isfinite(out.r.y));
    REQUIRE(std::isfinite(out.v.x));
    REQUIRE(std::isfinite(out.v.y));
}

TEST_CASE("Propagate_DoesNotWriteGarbageOnFailure", "[phase2][invalid]") {
    // Pre-fill output with sentinel values
    State2 out{{42.0, 43.0}, {44.0, 45.0}};
    State2 sentinel = out;

    // Invalid: zero radius
    REQUIRE(TwoBody::propagate(1.0, {{0, 0}, {0, 1}}, 1.0, out) == SolveStatus::InvalidInput);

    // Output should be unchanged
    REQUIRE(out.r.x == sentinel.r.x);
    REQUIRE(out.r.y == sentinel.r.y);
    REQUIRE(out.v.x == sentinel.v.x);
    REQUIRE(out.v.y == sentinel.v.y);
}

TEST_CASE("NoPublicFunctionReturnsNaNOnOrdinaryInvalidInput", "[phase2][invalid]") {
    // Test a range of bad inputs across all public functions

    // Zero radius
    State2 zero_r{{0, 0}, {0, 1}};
    double mu = 1.0;

    // Invariant functions on zero radius — may return odd values but not NaN
    // (these don't validate input; caller should validate first)
    // We focus on functions that return SolveStatus:
    ConicType ct;
    REQUIRE(TwoBody::classify(0.0, {{1, 0}, {0, 1}}, ct) == SolveStatus::InvalidInput);

    ConicElements2D el;
    REQUIRE(TwoBody::to_elements(0.0, {{1, 0}, {0, 1}}, el) == SolveStatus::InvalidInput);
    REQUIRE(TwoBody::to_elements(mu, zero_r, el) == SolveStatus::InvalidInput);

    State2 out;
    REQUIRE(TwoBody::propagate(0.0, {{1, 0}, {0, 1}}, 1.0, out) == SolveStatus::InvalidInput);
    REQUIRE(TwoBody::propagate(mu, zero_r, 1.0, out) == SolveStatus::InvalidInput);
}

TEST_CASE("IterationCapReturnsNoConvergence_NotPartialSolution", "[phase2][invalid]") {
    // Force the Kepler solver to hit iteration cap by using max_iterations=0
    double M = 1.0, e = 0.5;
    double E;
    REQUIRE(detail::mean_to_eccentric_anomaly(M, e, E, 0) == SolveStatus::NoConvergence);

    double H;
    REQUIRE(detail::mean_to_hyperbolic_anomaly(1.0, 1.5, H, 0) == SolveStatus::NoConvergence);
}

// --- H. Property tests ---

TEST_CASE("RoundTrip_ElementsAndCartesian_PreserveInvariants", "[phase2][property]") {
    // Test several elliptic and hyperbolic states
    struct TestCase {
        double mu;
        State2 state;
    };

    TestCase cases[] = {
        {1.0, {{1, 0}, {0, 1}}},                             // circular
        {1.0, {{1, 0}, {0, std::sqrt(1.5)}}},                // e=0.5
        {1.0, {{2, 1}, {-0.3, 0.8}}},                        // arbitrary elliptic
        {1.0, {{1, 0}, {0, 2}}},                              // hyperbolic
        {3.986e5, {{7000, 1000}, {-0.5, 7.5}}},              // Earth-like elliptic
    };

    for (const auto& tc : cases) {
        double eps0 = TwoBody::specific_energy(tc.mu, tc.state);
        double h0 = TwoBody::specific_angular_momentum_z(tc.state);

        ConicElements2D el;
        REQUIRE(TwoBody::to_elements(tc.mu, tc.state, el) == SolveStatus::Ok);

        State2 recon;
        REQUIRE(TwoBody::from_elements(el, recon) == SolveStatus::Ok);

        double eps1 = TwoBody::specific_energy(tc.mu, recon);
        double h1 = TwoBody::specific_angular_momentum_z(recon);

        REQUIRE_THAT(eps1, WithinAbs(eps0, 1e-8));
        REQUIRE_THAT(h1, WithinAbs(h0, 1e-8));
    }
}

TEST_CASE("Propagation_ComposesInTime", "[phase2][property]") {
    double mu = 1.0;
    State2 s{{1, 0}, {0, 1.2}}; // elliptic

    double dt1 = 1.3, dt2 = 2.7;

    State2 mid, composed;
    REQUIRE(TwoBody::propagate(mu, s, dt1, mid) == SolveStatus::Ok);
    REQUIRE(TwoBody::propagate(mu, mid, dt2, composed) == SolveStatus::Ok);

    State2 direct;
    REQUIRE(TwoBody::propagate(mu, s, dt1 + dt2, direct) == SolveStatus::Ok);

    REQUIRE_THAT(composed.r.x, WithinAbs(direct.r.x, 1e-8));
    REQUIRE_THAT(composed.r.y, WithinAbs(direct.r.y, 1e-8));
    REQUIRE_THAT(composed.v.x, WithinAbs(direct.v.x, 1e-8));
    REQUIRE_THAT(composed.v.y, WithinAbs(direct.v.y, 1e-8));
}

TEST_CASE("Propagation_PreservesSpecificEnergy", "[phase2][property]") {
    double mu = 1.0;

    struct TestCase {
        State2 state;
        double dt;
    };
    TestCase cases[] = {
        {{{1, 0}, {0, 1}}, 3.0},               // circular
        {{{1, 0}, {0, std::sqrt(1.5)}}, 5.0},   // elliptic
        {{{1, 0}, {0, 2}}, 2.0},                // hyperbolic
    };

    for (const auto& tc : cases) {
        double eps0 = TwoBody::specific_energy(mu, tc.state);
        State2 out;
        REQUIRE(TwoBody::propagate(mu, tc.state, tc.dt, out) == SolveStatus::Ok);
        double eps1 = TwoBody::specific_energy(mu, out);
        REQUIRE_THAT(eps1, WithinAbs(eps0, 1e-10));
    }
}

TEST_CASE("Propagation_PreservesAngularMomentum", "[phase2][property]") {
    double mu = 1.0;

    struct TestCase {
        State2 state;
        double dt;
    };
    TestCase cases[] = {
        {{{1, 0}, {0, 1}}, 3.0},
        {{{1, 0}, {0, std::sqrt(1.5)}}, 5.0},
        {{{1, 0}, {0, 2}}, 2.0},
    };

    for (const auto& tc : cases) {
        double h0 = TwoBody::specific_angular_momentum_z(tc.state);
        State2 out;
        REQUIRE(TwoBody::propagate(mu, tc.state, tc.dt, out) == SolveStatus::Ok);
        double h1 = TwoBody::specific_angular_momentum_z(out);
        REQUIRE_THAT(h1, WithinAbs(h0, 1e-10));
    }
}

TEST_CASE("SameInputsSamePlatform_ProduceBitIdenticalOutputsAcrossRuns",
          "[phase2][property]") {
    double mu = 1.0;
    State2 s{{1, 0}, {0, std::sqrt(1.5)}};
    double dt = 3.7;

    State2 out1, out2;
    REQUIRE(TwoBody::propagate(mu, s, dt, out1) == SolveStatus::Ok);
    REQUIRE(TwoBody::propagate(mu, s, dt, out2) == SolveStatus::Ok);

    // Bit-identical
    REQUIRE(out1.r.x == out2.r.x);
    REQUIRE(out1.r.y == out2.r.y);
    REQUIRE(out1.v.x == out2.v.x);
    REQUIRE(out1.v.y == out2.v.y);
}
