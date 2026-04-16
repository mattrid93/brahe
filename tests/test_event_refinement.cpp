#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "brahe/event_detector.h"

#include <cmath>

using namespace brahe;
using Catch::Matchers::WithinAbs;

// --- Section 7.3: Root-refinement tests ---

// Sign-change refinement

TEST_CASE("RefineRoot_FindsLinearSignChangeRoot", "[phase3][refine]") {
    // f(t) = t - 3. Root at t=3.
    auto f = [](double t) { return t - 3.0; };
    double out = -1.0;
    double f_lo = f(0.0);
    double f_hi = f(10.0);
    REQUIRE(detail::refine_root_bisection(f, 0.0, 10.0, f_lo, f_hi, 1e-12, 200, out) ==
            SolveStatus::Ok);
    REQUIRE_THAT(out, WithinAbs(3.0, 1e-10));
}

TEST_CASE("RefineRoot_FindsQuadraticRootWithinTolerance", "[phase3][refine]") {
    // f(t) = (t - 2)(t - 5). On [1, 3]: f(1)=4 > 0, f(3)=-2 < 0. Root at t=2.
    auto f = [](double t) { return (t - 2.0) * (t - 5.0); };
    double out = -1.0;
    double f_lo = f(1.0);
    double f_hi = f(3.0);
    REQUIRE(detail::refine_root_bisection(f, 1.0, 3.0, f_lo, f_hi, 1e-12, 200, out) ==
            SolveStatus::Ok);
    REQUIRE_THAT(out, WithinAbs(2.0, 1e-8));
}

TEST_CASE("RefineRoot_StopsWithinRootEpsilon", "[phase3][refine]") {
    auto f = [](double t) { return t - 7.0; };
    double tol = 1e-6;
    double out = -1.0;
    REQUIRE(detail::refine_root_bisection(f, 0.0, 10.0, f(0.0), f(10.0), tol, 200, out) ==
            SolveStatus::Ok);
    REQUIRE(std::abs(out - 7.0) < tol);
}

TEST_CASE("RefineRoot_ReturnsNoConvergenceOnIterationCap", "[phase3][refine]") {
    // Tight tolerance + very few iterations must return NoConvergence.
    auto f = [](double t) { return t - 5.0; };
    double out = -1.0;
    REQUIRE(detail::refine_root_bisection(f, 0.0, 10.0, f(0.0), f(10.0), 1e-15, 3, out) ==
            SolveStatus::NoConvergence);
    REQUIRE(std::isfinite(out));
}

// Grazing / minimum refinement

TEST_CASE("RefineMinimum_DetectsTangentialContactForParabola", "[phase3][refine]") {
    // f(t) = (t - 3)^2. Minimum at t=3 with value 0.
    auto f = [](double t) { return (t - 3.0) * (t - 3.0); };
    double out_t = -1.0;
    double out_f = -1.0;
    REQUIRE(detail::refine_minimum(f, 0.0, 6.0, 1e-10, 200, out_t, out_f) == SolveStatus::Ok);
    REQUIRE_THAT(out_t, WithinAbs(3.0, 1e-5));
    REQUIRE_THAT(out_f, WithinAbs(0.0, 1e-8));
}

TEST_CASE("RefineMinimum_RejectsNearMissOutsideRootEpsilon", "[phase3][refine]") {
    // f(t) = (t - 3)^2 + 0.1. Min at t=3 with value 0.1 (well above root_epsilon).
    auto f = [](double t) { return (t - 3.0) * (t - 3.0) + 0.1; };
    double out_t = -1.0;
    double out_f = -1.0;
    REQUIRE(detail::refine_minimum(f, 0.0, 6.0, 1e-10, 200, out_t, out_f) == SolveStatus::Ok);
    REQUIRE(out_f > 1e-6);  // well above root_epsilon so caller can reject
    REQUIRE_THAT(out_f, WithinAbs(0.1, 1e-6));
}

TEST_CASE("RefineMinimum_DoesNotEmitNaNOnFlatFunction", "[phase3][refine]") {
    // Constant function: no unique minimum. Must not emit NaN.
    auto f = [](double /*t*/) { return 5.0; };
    double out_t = -1.0;
    double out_f = -1.0;
    auto status = detail::refine_minimum(f, 0.0, 10.0, 1e-10, 200, out_t, out_f);
    (void)status;  // status can be Ok or NoConvergence, but outputs must be finite
    REQUIRE(std::isfinite(out_t));
    REQUIRE(std::isfinite(out_f));
    REQUIRE_THAT(out_f, WithinAbs(5.0, 1e-12));
    REQUIRE(out_t >= 0.0);
    REQUIRE(out_t <= 10.0);
}
