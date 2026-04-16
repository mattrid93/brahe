#include <catch2/catch_test_macros.hpp>

#include "brahe/body_system.h"
#include "brahe/patcher.h"
#include "brahe/trajectory.h"
#include "brahe/trajectory_builder.h"
#include "brahe/types.h"

#include <type_traits>

using namespace brahe;

// --- Section 7.1: Compile-only API tests ---

namespace {

BodySystem make_minimal_system() {
    BodySystemBuilder builder;
    BodyDef root{};
    root.id = 1;
    root.parent_id = InvalidBody;
    root.mu = 1.0;
    root.radius = 0.1;
    root.soi_radius = 1000.0;
    builder.add_body(root);

    BodyDef child{};
    child.id = 2;
    child.parent_id = 1;
    child.mu = 0.001;
    child.radius = 0.05;
    child.soi_radius = 5.0;
    child.orbit_radius = 50.0;
    child.angular_rate = 0.0;
    child.phase_at_epoch = 0.0;
    builder.add_body(child);

    BodySystem sys;
    (void)builder.build(sys);
    return sys;
}

}  // namespace

TEST_CASE("TrajectoryBuilder_CanBeConstructedFromBodySystem", "[phase4][api]") {
    BodySystem sys = make_minimal_system();
    TrajectoryBuilder builder(sys);
    (void)builder;
    SUCCEED();
}

TEST_CASE("BuildPreview_AcceptsPreviewRequestAndDynamicTrajectory", "[phase4][api]") {
    BodySystem sys = make_minimal_system();
    TrajectoryBuilder builder(sys);

    PreviewRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    req.initial_state = State2{{5.0, 0.0}, {0.0, 1.0}};
    req.end_time = 10.0;
    req.max_segments = 8;

    Trajectory out;
    SolveStatus status = builder.build_preview(req, out);
    (void)status;
    SUCCEED();
}

TEST_CASE("BuildPreviewFixed_AcceptsPreviewRequestAndFixedTrajectory", "[phase4][api]") {
    BodySystem sys = make_minimal_system();
    TrajectoryBuilder builder(sys);

    PreviewRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    req.initial_state = State2{{5.0, 0.0}, {0.0, 1.0}};
    req.end_time = 10.0;
    req.max_segments = 4;

    TrajectoryFixed<4> out;
    SolveStatus status = builder.build_preview_fixed(req, out);
    (void)status;
    SUCCEED();
}

TEST_CASE("Patcher_CanBeConstructedFromBodySystem", "[phase4][api]") {
    BodySystem sys = make_minimal_system();
    Patcher patcher(sys);
    (void)patcher;
    SUCCEED();
}

TEST_CASE("PatchParentToChild_HasStableSignature", "[phase4][api]") {
    BodySystem sys = make_minimal_system();
    Patcher patcher(sys);

    State2 sc_parent{{50.5, 0.0}, {0.0, 1.0}};
    State2 sc_child{};
    SolveStatus status = patcher.patch_parent_to_child(0.0, 1, 2, sc_parent, sc_child);
    (void)status;
    SUCCEED();
}

TEST_CASE("PatchChildToParent_HasStableSignature", "[phase4][api]") {
    BodySystem sys = make_minimal_system();
    Patcher patcher(sys);

    State2 sc_child{{0.5, 0.0}, {0.0, 1.0}};
    State2 sc_parent{};
    SolveStatus status = patcher.patch_child_to_parent(0.0, 2, sc_child, sc_parent);
    (void)status;
    SUCCEED();
}

TEST_CASE("Segment_IsTriviallyCopyable", "[phase4][api]") {
    STATIC_REQUIRE(std::is_trivially_copyable_v<Segment>);
}

TEST_CASE("Segment_IsStandardLayout", "[phase4][api]") {
    STATIC_REQUIRE(std::is_standard_layout_v<Segment>);
}

TEST_CASE("TrajectoryEvent_IsTriviallyCopyable", "[phase4][api]") {
    STATIC_REQUIRE(std::is_trivially_copyable_v<TrajectoryEvent>);
}

TEST_CASE("TrajectoryEvent_IsStandardLayout", "[phase4][api]") {
    STATIC_REQUIRE(std::is_standard_layout_v<TrajectoryEvent>);
}
