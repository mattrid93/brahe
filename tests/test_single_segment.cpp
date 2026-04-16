#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "brahe/body_system.h"
#include "brahe/trajectory.h"
#include "brahe/trajectory_builder.h"
#include "brahe/vec2.h"

#include <cmath>
#include <cstring>

using namespace brahe;
using Catch::Matchers::WithinAbs;

// --- Section 7.5: Segment emission behavior ---

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

}  // namespace

TEST_CASE("SingleSegment_EndsAtImpactEvent", "[phase4][single_segment]") {
    BodySystem sys = sun_with_moon();
    TrajectoryBuilder tb(sys);

    PreviewRequest req;
    req.central_body = kSun;
    req.start_time = 0.0;
    // Aimed directly into the sun: slight y-offset avoids a zero-h radial case.
    req.initial_state = State2{{5.0, 0.01}, {-3.0, 0.0}};
    req.end_time = 50.0;
    req.max_segments = 1;

    Trajectory out;
    SolveStatus status = tb.build_preview(req, out);
    REQUIRE((status == SolveStatus::Ok || status == SolveStatus::CapacityExceeded));
    REQUIRE(out.segments.size() == 1);
    REQUIRE(out.segments[0].end_reason == EventType::Impact);
    REQUIRE(out.segments[0].central_body == kSun);
}

TEST_CASE("SingleSegment_EndsAtSoiEntryEvent", "[phase4][single_segment]") {
    BodySystem sys = sun_with_moon();
    TrajectoryBuilder tb(sys);

    PreviewRequest req;
    req.central_body = kSun;
    req.start_time = 0.0;
    req.initial_state = State2{{30.0, 0.5}, {5.0, 0.0}};
    req.end_time = 20.0;
    req.max_segments = 1;  // cap forces chain to stop after first emitted segment

    Trajectory out;
    SolveStatus status = tb.build_preview(req, out);
    // Chain wants to patch to child next; capacity denies it.
    REQUIRE(status == SolveStatus::CapacityExceeded);
    REQUIRE(out.segments.size() == 1);
    REQUIRE(out.segments[0].end_reason == EventType::SoiEntry);
    REQUIRE(out.segments[0].central_body == kSun);
}

TEST_CASE("SingleSegment_EndsAtSoiExitEvent", "[phase4][single_segment]") {
    // Start directly in the moon's frame, heading outward. First event is SoiExit.
    BodySystem sys = sun_with_moon();
    TrajectoryBuilder tb(sys);

    PreviewRequest req;
    req.central_body = kMoon;
    req.start_time = 0.0;
    req.initial_state = State2{{4.0, 0.0}, {1.0, 0.0}};
    req.end_time = 20.0;
    req.max_segments = 1;

    Trajectory out;
    SolveStatus status = tb.build_preview(req, out);
    REQUIRE(status == SolveStatus::CapacityExceeded);
    REQUIRE(out.segments.size() == 1);
    REQUIRE(out.segments[0].end_reason == EventType::SoiExit);
    REQUIRE(out.segments[0].central_body == kMoon);
}

TEST_CASE("SingleSegment_EndsAtTimeLimitWhenNoEarlierEventOccurs",
          "[phase4][single_segment]") {
    BodySystem sys = sun_with_moon();
    TrajectoryBuilder tb(sys);

    PreviewRequest req;
    req.central_body = kSun;
    req.start_time = 0.0;
    // Circular orbit at r=10, well inside sun's SOI and away from the moon at (50,0).
    req.initial_state = State2{{10.0, 0.0}, {0.0, std::sqrt(1.0 / 10.0)}};
    req.end_time = 5.0;
    req.max_segments = 1;

    Trajectory out;
    SolveStatus status = tb.build_preview(req, out);
    REQUIRE(status == SolveStatus::Ok);
    REQUIRE(out.segments.size() == 1);
    REQUIRE(out.segments[0].end_reason == EventType::TimeLimit);
    REQUIRE(out.segments[0].central_body == kSun);
}

TEST_CASE("SingleSegment_StartAndEndTimesAreCorrect", "[phase4][single_segment]") {
    BodySystem sys = sun_with_moon();
    TrajectoryBuilder tb(sys);
    Tolerances tol{};

    const double t0 = 12.34;
    PreviewRequest req;
    req.central_body = kSun;
    req.start_time = t0;
    req.initial_state = State2{{10.0, 0.0}, {0.0, std::sqrt(1.0 / 10.0)}};
    req.end_time = t0 + 7.5;
    req.max_segments = 4;

    Trajectory out;
    REQUIRE(tb.build_preview(req, out) == SolveStatus::Ok);
    REQUIRE(!out.segments.empty());
    REQUIRE_THAT(out.segments[0].start_time, WithinAbs(t0, tol.time_epsilon));
    REQUIRE(out.segments.back().end_time > t0);
    REQUIRE_THAT(out.segments.back().end_time, WithinAbs(t0 + 7.5, tol.time_epsilon));
}

TEST_CASE("SingleSegment_StoresInitialStateExactly", "[phase4][single_segment]") {
    // The first segment's initial_state must match the request bitwise -- no
    // silent renormalisation or projection.
    BodySystem sys = sun_with_moon();
    TrajectoryBuilder tb(sys);

    State2 exact{{10.0, 0.0}, {0.0, std::sqrt(1.0 / 10.0)}};
    PreviewRequest req;
    req.central_body = kSun;
    req.start_time = 1.0;
    req.initial_state = exact;
    req.end_time = 2.0;
    req.max_segments = 1;

    Trajectory out;
    REQUIRE(tb.build_preview(req, out) == SolveStatus::Ok);
    REQUIRE(out.segments.size() == 1);

    // Compare the POD representations byte-for-byte.
    REQUIRE(std::memcmp(&out.segments[0].initial_state, &exact, sizeof(State2)) == 0);
}
