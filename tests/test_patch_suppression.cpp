#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "brahe/body_system.h"
#include "brahe/trajectory.h"
#include "brahe/trajectory_builder.h"
#include "brahe/vec2.h"

#include <cmath>

using namespace brahe;
using Catch::Matchers::WithinAbs;

// --- Section 7.4: Re-entry suppression (observed through preview chains) ---

namespace {

constexpr BodyId kSun = 1;
constexpr BodyId kMoonA = 2;
constexpr BodyId kMoonB = 3;

// Sun + single stationary moon. Planet radius is much smaller than SOI so a
// grazing hyperbolic trajectory enters -> transits -> exits cleanly.
BodySystem sun_with_one_stationary_moon(double moon_soi, double moon_radius,
                                        double moon_orbit_radius) {
    BodySystemBuilder b;
    BodyDef sun{};
    sun.id = kSun;
    sun.parent_id = InvalidBody;
    sun.mu = 1.0;
    sun.radius = 0.1;
    sun.soi_radius = 1e6;
    b.add_body(sun);

    BodyDef moon{};
    moon.id = kMoonA;
    moon.parent_id = kSun;
    moon.mu = 0.001;
    moon.radius = moon_radius;
    moon.soi_radius = moon_soi;
    moon.orbit_radius = moon_orbit_radius;
    moon.angular_rate = 0.0;
    moon.phase_at_epoch = 0.0;
    b.add_body(moon);

    BodySystem sys;
    (void)b.build(sys);
    return sys;
}

// Sun + two stationary moons on the +x axis. Spacecraft can be aimed to
// traverse moon A then moon B in sequence.
BodySystem sun_with_two_stationary_moons(double soi, double radius,
                                         double orbit_a, double orbit_b) {
    BodySystemBuilder b;
    BodyDef sun{};
    sun.id = kSun;
    sun.parent_id = InvalidBody;
    sun.mu = 1.0;
    sun.radius = 0.1;
    sun.soi_radius = 1e6;
    b.add_body(sun);

    BodyDef ma{};
    ma.id = kMoonA;
    ma.parent_id = kSun;
    ma.mu = 0.001;
    ma.radius = radius;
    ma.soi_radius = soi;
    ma.orbit_radius = orbit_a;
    ma.angular_rate = 0.0;
    ma.phase_at_epoch = 0.0;
    b.add_body(ma);

    BodyDef mb{};
    mb.id = kMoonB;
    mb.parent_id = kSun;
    mb.mu = 0.001;
    mb.radius = radius;
    mb.soi_radius = soi;
    mb.orbit_radius = orbit_b;
    mb.angular_rate = 0.0;
    mb.phase_at_epoch = 0.0;
    b.add_body(mb);

    BodySystem sys;
    (void)b.build(sys);
    return sys;
}

// A small helper that builds a preview chain and sanity-checks completion.
SolveStatus build_chain(const BodySystem& sys, const Tolerances& tol,
                        BodyId central_body, double t0, const State2& s0, double t_end,
                        size_t cap, Trajectory& out) {
    TrajectoryBuilder tb(sys, tol);
    PreviewRequest req;
    req.central_body = central_body;
    req.start_time = t0;
    req.initial_state = s0;
    req.end_time = t_end;
    req.max_segments = cap;
    return tb.build_preview(req, out);
}

}  // namespace

TEST_CASE("Suppression_ActivatesImmediatelyAfterParentToChildPatch",
          "[phase4][suppression]") {
    // Hyperbolic flyby heading through moon A's SOI. After SoiEntry, the child
    // segment must have non-trivial duration: it must not immediately terminate
    // at SoiExit of the same boundary (which would indicate suppression failure).
    BodySystem sys = sun_with_one_stationary_moon(5.0, 0.2, 50.0);
    Tolerances tol{};

    Trajectory chain;
    REQUIRE(build_chain(sys, tol, kSun, 0.0, State2{{30.0, 0.5}, {5.0, 0.0}}, 40.0, 16,
                        chain) == SolveStatus::Ok);

    REQUIRE(chain.segments.size() >= 2);
    REQUIRE(chain.segments[0].end_reason == EventType::SoiEntry);

    // The second (child-frame) segment must have meaningful duration.
    const Segment& child_seg = chain.segments[1];
    REQUIRE(child_seg.central_body == kMoonA);
    REQUIRE((child_seg.end_time - child_seg.start_time) > tol.time_epsilon);
}

TEST_CASE("Suppression_ActivatesImmediatelyAfterChildToParentPatch",
          "[phase4][suppression]") {
    // After the child-frame transit and SoiExit, the new parent-frame segment
    // must also have non-trivial duration -- we must not immediately re-enter
    // the moon we just exited.
    BodySystem sys = sun_with_one_stationary_moon(5.0, 0.2, 50.0);
    Tolerances tol{};

    Trajectory chain;
    REQUIRE(build_chain(sys, tol, kSun, 0.0, State2{{30.0, 0.5}, {5.0, 0.0}}, 40.0, 16,
                        chain) == SolveStatus::Ok);

    REQUIRE(chain.segments.size() >= 3);
    REQUIRE(chain.segments[1].end_reason == EventType::SoiExit);

    const Segment& post_exit = chain.segments[2];
    REQUIRE(post_exit.central_body == kSun);
    REQUIRE((post_exit.end_time - post_exit.start_time) > tol.time_epsilon);
}

TEST_CASE("Suppression_BlocksOnlyTheCrossedBoundary", "[phase4][suppression]") {
    // After entering moon A and transiting to the other side, the spacecraft
    // continues toward moon B. Moon B's SOI entry must not be suppressed -- it
    // is a different boundary. The chain should therefore include a second
    // SoiEntry event for moon B.
    BodySystem sys = sun_with_two_stationary_moons(2.0, 0.2, 50.0, 70.0);
    Tolerances tol{};

    Trajectory chain;
    REQUIRE(build_chain(sys, tol, kSun, 0.0, State2{{47.0, 0.5}, {5.0, 0.0}}, 20.0, 16,
                        chain) == SolveStatus::Ok);

    // Expected sequence: SoiEntry(A), SoiExit(A), SoiEntry(B), SoiExit(B), ...
    bool saw_entry_a = false;
    bool saw_exit_a = false;
    bool saw_entry_b = false;
    for (const Segment& s : chain.segments) {
        if (s.end_reason == EventType::SoiEntry && s.central_body == kSun) {
            // Look at the next segment's central body for which SOI we entered.
            size_t idx = static_cast<size_t>(&s - chain.segments.data());
            if (idx + 1 < chain.segments.size()) {
                BodyId entered = chain.segments[idx + 1].central_body;
                if (entered == kMoonA) saw_entry_a = true;
                if (entered == kMoonB) saw_entry_b = true;
            }
        }
        if (s.end_reason == EventType::SoiExit && s.central_body == kMoonA) saw_exit_a = true;
    }
    REQUIRE(saw_entry_a);
    REQUIRE(saw_exit_a);
    REQUIRE(saw_entry_b);
}

TEST_CASE("Suppression_ExpiresAfterTimeEpsilonEvenIfStillNearBoundary",
          "[phase4][suppression]") {
    // A very short time_epsilon lets the suppression lapse quickly. The chain
    // must still reach the natural SOI exit on the far side of the moon,
    // producing a properly-terminated SoiExit segment -- not a stalled chain.
    BodySystem sys = sun_with_one_stationary_moon(5.0, 0.2, 50.0);
    Tolerances tol{};
    tol.time_epsilon = 1e-9;

    Trajectory chain;
    REQUIRE(build_chain(sys, tol, kSun, 0.0, State2{{30.0, 0.5}, {5.0, 0.0}}, 40.0, 16,
                        chain) == SolveStatus::Ok);

    // Look for a child-frame segment that cleanly terminates at SoiExit.
    bool found_clean_exit = false;
    for (const Segment& s : chain.segments) {
        if (s.central_body == kMoonA && s.end_reason == EventType::SoiExit) {
            REQUIRE((s.end_time - s.start_time) > tol.time_epsilon);
            found_clean_exit = true;
        }
    }
    REQUIRE(found_clean_exit);
}

TEST_CASE("Suppression_ExpiresAfterMovingAwayByPositionThresholdEvenBeforeTimeEpsilon",
          "[phase4][suppression]") {
    // Wide position tolerance vs tight time tolerance: the spacecraft moves
    // many position-epsilons away during transit, so suppression must have
    // lifted by the time it approaches the far boundary. We observe this as a
    // successfully emitted SoiExit for the crossed moon.
    BodySystem sys = sun_with_one_stationary_moon(5.0, 0.2, 50.0);
    Tolerances tol{};
    tol.position_epsilon = 1e-6;

    Trajectory chain;
    REQUIRE(build_chain(sys, tol, kSun, 0.0, State2{{30.0, 0.5}, {5.0, 0.0}}, 40.0, 16,
                        chain) == SolveStatus::Ok);

    bool found_clean_exit = false;
    for (const Segment& s : chain.segments) {
        if (s.central_body == kMoonA && s.end_reason == EventType::SoiExit) {
            found_clean_exit = true;
        }
    }
    REQUIRE(found_clean_exit);
}

TEST_CASE("Suppression_DoesNotBlockImpactDetection", "[phase4][suppression]") {
    // Spacecraft enters moon SOI and is aimed directly at the body. Impact
    // must be detected even though it occurs inside the just-entered SOI.
    // Use a larger moon radius so the direct hit is unambiguous.
    BodySystem sys = sun_with_one_stationary_moon(5.0, 1.5, 50.0);
    Tolerances tol{};

    Trajectory chain;
    REQUIRE(build_chain(sys, tol, kSun, 0.0, State2{{30.0, 0.001}, {5.0, 0.0}}, 20.0, 16,
                        chain) == SolveStatus::Ok);

    REQUIRE(chain.segments.size() >= 2);
    REQUIRE(chain.segments.back().end_reason == EventType::Impact);
    REQUIRE(chain.segments.back().central_body == kMoonA);
}

TEST_CASE("Suppression_DoesNotBlockUnrelatedSiblingChildEntry", "[phase4][suppression]") {
    // After exiting moon A, heading straight into moon B. Moon B's SoiEntry
    // must not be suppressed -- it's a completely different boundary.
    BodySystem sys = sun_with_two_stationary_moons(2.0, 0.2, 50.0, 70.0);
    Tolerances tol{};

    Trajectory chain;
    REQUIRE(build_chain(sys, tol, kSun, 0.0, State2{{47.0, 0.5}, {5.0, 0.0}}, 20.0, 16,
                        chain) == SolveStatus::Ok);

    // Must find a Moon B segment (only possible if entry wasn't suppressed).
    bool saw_moon_b = false;
    for (const Segment& s : chain.segments) {
        if (s.central_body == kMoonB) saw_moon_b = true;
    }
    REQUIRE(saw_moon_b);
}
