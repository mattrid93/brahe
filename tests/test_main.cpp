#include <catch2/catch_test_macros.hpp>

#include "brahe/brahe.h"

TEST_CASE("add returns correct sum", "[add]") {
    REQUIRE(brahe::add(2, 3) == 5);
    REQUIRE(brahe::add(-1, 1) == 0);
    REQUIRE(brahe::add(0, 0) == 0);
}
