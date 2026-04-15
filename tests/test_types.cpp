#include <catch2/catch_test_macros.hpp>

#include "brahe/types.h"

#include <type_traits>

using namespace brahe;

// --- A. Core layout and type tests ---

TEST_CASE("BodyId_InvalidBody_IsAllOnes", "[types]") {
    REQUIRE(InvalidBody == ~BodyId{0});
    REQUIRE(InvalidBody == 0xFFFFFFFF);
}

TEST_CASE("Vec2_State2_BodyDef_AreTriviallyCopyable", "[types]") {
    static_assert(std::is_trivially_copyable_v<Vec2>);
    static_assert(std::is_trivially_copyable_v<State2>);
    static_assert(std::is_trivially_copyable_v<BodyDef>);
    SUCCEED();
}

TEST_CASE("Vec2_State2_BodyDef_AreStandardLayout", "[types]") {
    static_assert(std::is_standard_layout_v<Vec2>);
    static_assert(std::is_standard_layout_v<State2>);
    static_assert(std::is_standard_layout_v<BodyDef>);
    SUCCEED();
}

TEST_CASE("DefaultTolerances_MatchSpec", "[types]") {
    Tolerances tol;
    REQUIRE(tol.position_epsilon == 1e-3);
    REQUIRE(tol.velocity_epsilon == 1e-6);
    REQUIRE(tol.angle_epsilon == 1e-9);
    REQUIRE(tol.time_epsilon == 1e-6);
    REQUIRE(tol.root_epsilon == 1e-6);
    REQUIRE(tol.lambert_residual_epsilon == 1e-10);
    REQUIRE(tol.parabolic_eccentricity_band == 1e-6);
}
