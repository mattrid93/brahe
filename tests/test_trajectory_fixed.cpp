#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "brahe/body_system.h"
#include "brahe/trajectory.h"
#include "brahe/trajectory_builder.h"

#include <array>
#include <cstdint>
#include <cstring>

using namespace brahe;

// --- Section 7.7: Fixed-capacity tests ---

namespace {

constexpr BodyId kSun = 1;
constexpr BodyId kMoon = 2;

BodySystem sun_with_moon() {
    BodySystemBuilder b;
    BodyDef sun{};
    sun.id = kSun;
    sun.parent_id = InvalidBody;
    sun.mu = 1.0;
    sun.radius = 0.1;
    sun.soi_radius = 1e6;
    b.add_body(sun);

    BodyDef moon{};
    moon.id = kMoon;
    moon.parent_id = kSun;
    moon.mu = 0.001;
    moon.radius = 0.2;
    moon.soi_radius = 5.0;
    moon.orbit_radius = 50.0;
    moon.angular_rate = 0.0;
    moon.phase_at_epoch = 0.0;
    b.add_body(moon);

    BodySystem sys;
    (void)b.build(sys);
    return sys;
}

PreviewRequest flyby(size_t cap) {
    PreviewRequest req;
    req.central_body = kSun;
    req.start_time = 0.0;
    req.initial_state = State2{{30.0, 0.5}, {5.0, 0.0}};
    req.end_time = 40.0;
    req.max_segments = cap;
    return req;
}

}  // namespace

TEST_CASE("FixedCapacity_SucceedsWhenSegmentCountFits", "[phase4][fixed_capacity]") {
    // Flyby produces three segments; capacity 8 is generous.
    BodySystem sys = sun_with_moon();
    TrajectoryBuilder tb(sys);

    TrajectoryFixed<8> out;
    SolveStatus status = tb.build_preview_fixed(flyby(8), out);
    REQUIRE(status == SolveStatus::Ok);
    REQUIRE(out.count >= 3);
    REQUIRE(out.count <= 8);
}

TEST_CASE("FixedCapacity_ReturnsCapacityFailureWhenSegmentCountWouldOverflow",
          "[phase4][fixed_capacity]") {
    // Flyby needs at least 3 segments; capacity 2 forces CapacityExceeded.
    BodySystem sys = sun_with_moon();
    TrajectoryBuilder tb(sys);

    TrajectoryFixed<2> out;
    SolveStatus status = tb.build_preview_fixed(flyby(2), out);
    REQUIRE(status == SolveStatus::CapacityExceeded);
    REQUIRE(out.count <= 2);
}

TEST_CASE("FixedCapacity_DoesNotWritePastArrayBounds", "[phase4][fixed_capacity]") {
    // Surround the TrajectoryFixed with guard regions of a known bit pattern
    // and verify the guards are untouched after a chain that would overflow
    // the fixed buffer.
    struct Guarded {
        std::array<uint64_t, 16> front;
        TrajectoryFixed<2> fixed;
        std::array<uint64_t, 16> back;
    };

    Guarded g;
    constexpr uint64_t kGuard = 0xDEADBEEFCAFEBABEULL;
    for (auto& w : g.front) w = kGuard;
    for (auto& w : g.back) w = kGuard;

    BodySystem sys = sun_with_moon();
    TrajectoryBuilder tb(sys);

    SolveStatus status = tb.build_preview_fixed(flyby(2), g.fixed);
    REQUIRE(status == SolveStatus::CapacityExceeded);
    REQUIRE(g.fixed.count <= 2);

    for (auto w : g.front) REQUIRE(w == kGuard);
    for (auto w : g.back) REQUIRE(w == kGuard);
}

TEST_CASE("FixedCapacity_PreservesValidPrefixOnFailure", "[phase4][fixed_capacity]") {
    // Build twice: once with ample capacity (reference), once with a tight cap.
    // The cap'd prefix must match the first `count` segments of the reference.
    BodySystem sys = sun_with_moon();
    TrajectoryBuilder tb(sys);

    TrajectoryFixed<8> ref;
    REQUIRE(tb.build_preview_fixed(flyby(8), ref) == SolveStatus::Ok);

    TrajectoryFixed<2> capped;
    REQUIRE(tb.build_preview_fixed(flyby(2), capped) == SolveStatus::CapacityExceeded);

    REQUIRE(capped.count <= ref.count);
    for (size_t i = 0; i < capped.count; ++i) {
        // Byte-identical segment prefix.
        REQUIRE(std::memcmp(&capped.segments[i], &ref.segments[i], sizeof(Segment)) == 0);
    }
}
