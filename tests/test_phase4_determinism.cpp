#include <catch2/catch_test_macros.hpp>

#include "brahe/body_system.h"
#include "brahe/trajectory.h"
#include "brahe/trajectory_builder.h"

#include <cmath>
#include <cstring>

using namespace brahe;

// --- Section 7.8: Determinism tests ---

namespace {

constexpr BodyId kSun = 1;
constexpr BodyId kMoonA = 2;
constexpr BodyId kMoonB = 3;

BodyDef make_sun() {
    BodyDef d{};
    d.id = kSun;
    d.parent_id = InvalidBody;
    d.mu = 1.0;
    d.radius = 0.1;
    d.soi_radius = 1e6;
    return d;
}

BodyDef make_moon(BodyId id, double orbit_radius, double phase = 0.0) {
    BodyDef d{};
    d.id = id;
    d.parent_id = kSun;
    d.mu = 0.001;
    d.radius = 0.2;
    d.soi_radius = 5.0;
    d.orbit_radius = orbit_radius;
    d.angular_rate = 0.0;
    d.phase_at_epoch = phase;
    return d;
}

BodySystem build_with_order(const BodyDef* defs, size_t n) {
    BodySystemBuilder b;
    for (size_t i = 0; i < n; ++i) b.add_body(defs[i]);
    BodySystem sys;
    (void)b.build(sys);
    return sys;
}

bool segments_bit_identical(const Trajectory& a, const Trajectory& b) {
    if (a.segments.size() != b.segments.size()) return false;
    for (size_t i = 0; i < a.segments.size(); ++i) {
        if (std::memcmp(&a.segments[i], &b.segments[i], sizeof(Segment)) != 0) return false;
    }
    return true;
}

}  // namespace

TEST_CASE("BuildPreview_RepeatedRunsProduceBitIdenticalSegmentsOnSamePlatform",
          "[phase4][determinism]") {
    BodyDef defs[] = {make_sun(), make_moon(kMoonA, 50.0)};
    BodySystem sys = build_with_order(defs, 2);

    TrajectoryBuilder tb(sys);
    PreviewRequest req;
    req.central_body = kSun;
    req.start_time = 0.0;
    req.initial_state = State2{{30.0, 0.5}, {5.0, 0.0}};
    req.end_time = 40.0;
    req.max_segments = 16;

    Trajectory a, b, c;
    REQUIRE(tb.build_preview(req, a) == SolveStatus::Ok);
    REQUIRE(tb.build_preview(req, b) == SolveStatus::Ok);
    REQUIRE(tb.build_preview(req, c) == SolveStatus::Ok);

    REQUIRE(segments_bit_identical(a, b));
    REQUIRE(segments_bit_identical(a, c));
}

TEST_CASE("BuildPreview_IsIndependentOfBodyInsertionOrder", "[phase4][determinism]") {
    // Same bodies, different insertion orders. BodySystemBuilder should
    // canonicalise by BodyId -> outputs must be identical.
    BodyDef forward[] = {make_sun(), make_moon(kMoonA, 50.0), make_moon(kMoonB, 70.0)};
    BodyDef reverse[] = {make_moon(kMoonB, 70.0), make_moon(kMoonA, 50.0), make_sun()};

    BodySystem sys_fwd = build_with_order(forward, 3);
    BodySystem sys_rev = build_with_order(reverse, 3);

    PreviewRequest req;
    req.central_body = kSun;
    req.start_time = 0.0;
    req.initial_state = State2{{47.0, 0.01}, {5.0, 0.0}};
    req.end_time = 20.0;
    req.max_segments = 16;

    Trajectory out_fwd, out_rev;
    REQUIRE(TrajectoryBuilder(sys_fwd).build_preview(req, out_fwd) == SolveStatus::Ok);
    REQUIRE(TrajectoryBuilder(sys_rev).build_preview(req, out_rev) == SolveStatus::Ok);

    REQUIRE(segments_bit_identical(out_fwd, out_rev));
}

TEST_CASE("ChildEntryTieBreak_RemainsStableAcrossRepeatedRuns", "[phase4][determinism]") {
    // Two non-overlapping children on the +X axis. At high spacecraft speed, their
    // entry events fall within time_epsilon. Phase 3 resolves this by lowest BodyId;
    // Phase 4 must preserve that choice identically run-to-run.
    BodyDef defs[] = {make_sun(), make_moon(10, 50.0, 0.0), make_moon(20, 60.1, 0.0)};
    BodySystem sys = build_with_order(defs, 3);
    TrajectoryBuilder tb(sys);

    PreviewRequest req;
    req.central_body = kSun;
    req.start_time = 0.0;
    req.initial_state = State2{{30.0, 0.0}, {20000000.0, 0.0}};
    req.end_time = 0.00001;
    req.max_segments = 8;

    Trajectory a, b;
    SolveStatus sa = tb.build_preview(req, a);
    SolveStatus sb = tb.build_preview(req, b);
    REQUIRE(sa == sb);
    REQUIRE(segments_bit_identical(a, b));

    // The first patch must pick BodyId 10 (lower of 10/20).
    REQUIRE(a.segments.size() >= 2);
    REQUIRE(a.segments[0].end_reason == EventType::SoiEntry);
    REQUIRE(a.segments[1].central_body == 10);
}

TEST_CASE("SuppressionBehavior_IsStableAcrossRepeatedRuns", "[phase4][determinism]") {
    // A near-boundary start that exercises suppression. Its outcome must be
    // bit-identical across repeated runs.
    BodyDef defs[] = {make_sun(), make_moon(kMoonA, 50.0)};
    BodySystem sys = build_with_order(defs, 2);
    TrajectoryBuilder tb(sys);

    PreviewRequest req;
    req.central_body = kSun;
    req.start_time = 0.0;
    // Essentially on the SOI boundary (SOI=5 around (50,0)): x = 55 - 1e-10.
    req.initial_state = State2{{55.0 - 1e-10, 0.01}, {0.0, 0.3}};
    req.end_time = 5.0;
    req.max_segments = 16;

    Trajectory a, b;
    SolveStatus sa = tb.build_preview(req, a);
    SolveStatus sb = tb.build_preview(req, b);
    REQUIRE(sa == sb);
    REQUIRE(segments_bit_identical(a, b));
}

TEST_CASE("DynamicAndFixedPreviewAgreeOnPrefix", "[phase4][determinism]") {
    // The dynamic and fixed-capacity paths share one implementation; their
    // emitted prefixes must match byte-for-byte.
    BodyDef defs[] = {make_sun(), make_moon(kMoonA, 50.0)};
    BodySystem sys = build_with_order(defs, 2);
    TrajectoryBuilder tb(sys);

    PreviewRequest req;
    req.central_body = kSun;
    req.start_time = 0.0;
    req.initial_state = State2{{30.0, 0.5}, {5.0, 0.0}};
    req.end_time = 40.0;
    req.max_segments = 16;

    Trajectory dyn;
    TrajectoryFixed<16> fxd;
    REQUIRE(tb.build_preview(req, dyn) == SolveStatus::Ok);
    REQUIRE(tb.build_preview_fixed(req, fxd) == SolveStatus::Ok);

    REQUIRE(dyn.segments.size() == fxd.count);
    for (size_t i = 0; i < dyn.segments.size(); ++i) {
        REQUIRE(std::memcmp(&dyn.segments[i], &fxd.segments[i], sizeof(Segment)) == 0);
    }
}
