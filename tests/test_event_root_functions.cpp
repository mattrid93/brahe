#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "brahe/event_detector.h"

#include <cmath>

using namespace brahe;
using Catch::Matchers::WithinAbs;

// --- Section 7.2: Root-function unit tests ---

// Impact root tests

TEST_CASE("ImpactFunction_IsNegativeInsideRadius", "[phase3][root]") {
    REQUIRE(detail::impact_function({0.5, 0.0}, 1.0) < 0.0);
    REQUIRE(detail::impact_function({0.3, 0.4}, 1.0) < 0.0);  // |r|=0.5
}

TEST_CASE("ImpactFunction_IsZeroOnRadius", "[phase3][root]") {
    REQUIRE_THAT(detail::impact_function({1.0, 0.0}, 1.0), WithinAbs(0.0, 1e-15));
    REQUIRE_THAT(detail::impact_function({0.6, 0.8}, 1.0), WithinAbs(0.0, 1e-15));  // |r|=1
}

TEST_CASE("ImpactFunction_IsPositiveOutsideRadius", "[phase3][root]") {
    REQUIRE(detail::impact_function({2.0, 0.0}, 1.0) > 0.0);
    REQUIRE(detail::impact_function({3.0, 4.0}, 1.0) > 0.0);  // |r|=5
}

// SOI exit root tests

TEST_CASE("SoiExitFunction_IsNegativeInsideSoi", "[phase3][root]") {
    REQUIRE(detail::soi_exit_function({5.0, 0.0}, 10.0) < 0.0);
    REQUIRE(detail::soi_exit_function({3.0, 4.0}, 10.0) < 0.0);
}

TEST_CASE("SoiExitFunction_IsZeroOnSoiBoundary", "[phase3][root]") {
    REQUIRE_THAT(detail::soi_exit_function({10.0, 0.0}, 10.0), WithinAbs(0.0, 1e-15));
    REQUIRE_THAT(detail::soi_exit_function({6.0, 8.0}, 10.0), WithinAbs(0.0, 1e-15));
}

TEST_CASE("SoiExitFunction_IsPositiveOutsideSoi", "[phase3][root]") {
    REQUIRE(detail::soi_exit_function({15.0, 0.0}, 10.0) > 0.0);
    REQUIRE(detail::soi_exit_function({0.0, 20.0}, 10.0) > 0.0);
}

// Child entry root tests

TEST_CASE("ChildEntryFunction_IsPositiveOutsideChildSoi", "[phase3][root]") {
    // spacecraft at (100, 0), child at (50, 0), child SOI=5. Distance=50, outside.
    REQUIRE(detail::child_entry_function({100.0, 0.0}, {50.0, 0.0}, 5.0) > 0.0);
}

TEST_CASE("ChildEntryFunction_IsZeroOnChildSoiBoundary", "[phase3][root]") {
    // spacecraft at (55, 0), child at (50, 0). Distance=5, on boundary.
    REQUIRE_THAT(detail::child_entry_function({55.0, 0.0}, {50.0, 0.0}, 5.0),
                 WithinAbs(0.0, 1e-15));
}

TEST_CASE("ChildEntryFunction_IsNegativeInsideChildSoi", "[phase3][root]") {
    // spacecraft at (52, 0), child at (50, 0). Distance=2, inside SOI.
    REQUIRE(detail::child_entry_function({52.0, 0.0}, {50.0, 0.0}, 5.0) < 0.0);
}

// Numeric guard tests

TEST_CASE("RootFunctions_DoNotEmitNaN_ForFiniteInputs", "[phase3][root]") {
    REQUIRE(std::isfinite(detail::impact_function({0.0, 0.0}, 1.0)));
    REQUIRE(std::isfinite(detail::soi_exit_function({0.0, 0.0}, 1.0)));
    REQUIRE(std::isfinite(detail::child_entry_function({0.0, 0.0}, {0.0, 0.0}, 1.0)));

    // Typical orbital values
    REQUIRE(std::isfinite(detail::impact_function({7000.0, 1000.0}, 6378.0)));
    REQUIRE(std::isfinite(detail::child_entry_function({7e8, 0.0}, {3.84e8, 0.0}, 6.6e7)));
}

TEST_CASE("RootFunctions_HandleVerySmallPositiveDistances", "[phase3][root]") {
    REQUIRE(std::isfinite(detail::impact_function({1e-15, 0.0}, 1.0)));
    REQUIRE(std::isfinite(detail::soi_exit_function({1e-300, 1e-300}, 1.0)));
    REQUIRE(std::isfinite(detail::child_entry_function({1e-15, 0.0}, {0.0, 0.0}, 1.0)));
}

TEST_CASE("RootFunctions_HandleBoundaryValuesWithinRootEpsilon", "[phase3][root]") {
    double root_eps = 1e-6;
    double r_just_outside = 1.0 + 0.5 * root_eps;
    double r_just_inside = 1.0 - 0.5 * root_eps;

    double f_out = detail::impact_function({r_just_outside, 0.0}, 1.0);
    double f_in = detail::impact_function({r_just_inside, 0.0}, 1.0);

    REQUIRE(std::isfinite(f_out));
    REQUIRE(std::isfinite(f_in));
    REQUIRE(std::abs(f_out) < root_eps);
    REQUIRE(std::abs(f_in) < root_eps);
}
