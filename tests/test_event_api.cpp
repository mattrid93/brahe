#include <catch2/catch_test_macros.hpp>

#include "brahe/body_system.h"
#include "brahe/event_detector.h"
#include "brahe/types.h"

#include <type_traits>

using namespace brahe;

// --- Section 7.1: Compile-only API tests ---

namespace {

BodySystem make_single_body_system() {
    BodySystemBuilder builder;
    BodyDef root{};
    root.id = 1;
    root.parent_id = InvalidBody;
    root.mu = 1.0;
    root.radius = 1.0;
    root.soi_radius = 100.0;
    builder.add_body(root);
    BodySystem sys;
    (void)builder.build(sys);
    return sys;
}

}  // namespace

TEST_CASE("EventDetector_CanBeConstructedFromBodySystem", "[phase3][api]") {
    BodySystem sys = make_single_body_system();
    EventDetector detector(sys);
    (void)detector;
    SUCCEED();
}

TEST_CASE("FindNextEvent_AcceptsSearchRequestAndOutputEvent", "[phase3][api]") {
    BodySystem sys = make_single_body_system();
    EventDetector detector(sys);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    req.initial_state = State2{{5.0, 0.0}, {0.0, 1.0}};
    req.time_limit = 10.0;

    PredictedEvent event;
    SolveStatus status = detector.find_next_event(req, event);
    (void)status;
    SUCCEED();
}

TEST_CASE("PredictedEvent_IsTriviallyCopyable", "[phase3][api]") {
    STATIC_REQUIRE(std::is_trivially_copyable_v<PredictedEvent>);
}

TEST_CASE("PredictedEvent_IsStandardLayout", "[phase3][api]") {
    STATIC_REQUIRE(std::is_standard_layout_v<PredictedEvent>);
}

TEST_CASE("EventSearchRequest_IsTriviallyCopyable", "[phase3][api]") {
    STATIC_REQUIRE(std::is_trivially_copyable_v<EventSearchRequest>);
}

TEST_CASE("EventSearchRequest_IsStandardLayout", "[phase3][api]") {
    STATIC_REQUIRE(std::is_standard_layout_v<EventSearchRequest>);
}
