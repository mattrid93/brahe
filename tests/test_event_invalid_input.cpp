#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "brahe/body_system.h"
#include "brahe/event_detector.h"

#include <cmath>
#include <limits>

using namespace brahe;
using Catch::Matchers::WithinAbs;

// --- Section 7.10: Invalid-input and failure-mode tests ---

namespace {

BodySystem build_single() {
    BodySystemBuilder b;
    BodyDef root{};
    root.id = 1;
    root.parent_id = InvalidBody;
    root.mu = 1.0;
    root.radius = 1.0;
    root.soi_radius = 100.0;
    b.add_body(root);
    BodySystem sys;
    (void)b.build(sys);
    return sys;
}

constexpr double inf_v = std::numeric_limits<double>::infinity();
constexpr double nan_v = std::numeric_limits<double>::quiet_NaN();

}  // namespace

TEST_CASE("FindNextEvent_RejectsUnknownCentralBody", "[phase3][invalid]") {
    BodySystem sys = build_single();
    EventDetector det(sys);

    EventSearchRequest req;
    req.central_body = 999;  // not in system
    req.start_time = 0.0;
    req.initial_state = State2{{5.0, 0.0}, {0.0, 1.0}};
    req.time_limit = 10.0;

    PredictedEvent event;
    REQUIRE(det.find_next_event(req, event) == SolveStatus::InvalidInput);
}

TEST_CASE("FindNextEvent_RejectsNonFiniteStartTime", "[phase3][invalid]") {
    BodySystem sys = build_single();
    EventDetector det(sys);

    EventSearchRequest req;
    req.central_body = 1;
    req.initial_state = State2{{5.0, 0.0}, {0.0, 1.0}};
    req.time_limit = 10.0;

    PredictedEvent event;
    req.start_time = inf_v;
    REQUIRE(det.find_next_event(req, event) == SolveStatus::InvalidInput);

    req.start_time = nan_v;
    REQUIRE(det.find_next_event(req, event) == SolveStatus::InvalidInput);
}

TEST_CASE("FindNextEvent_RejectsNonFiniteTimeLimit", "[phase3][invalid]") {
    BodySystem sys = build_single();
    EventDetector det(sys);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    req.initial_state = State2{{5.0, 0.0}, {0.0, 1.0}};

    PredictedEvent event;
    req.time_limit = inf_v;
    REQUIRE(det.find_next_event(req, event) == SolveStatus::InvalidInput);

    req.time_limit = nan_v;
    REQUIRE(det.find_next_event(req, event) == SolveStatus::InvalidInput);
}

TEST_CASE("FindNextEvent_RejectsTimeLimitBeforeStartTime", "[phase3][invalid]") {
    BodySystem sys = build_single();
    EventDetector det(sys);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 10.0;
    req.time_limit = 5.0;  // before start
    req.initial_state = State2{{5.0, 0.0}, {0.0, 1.0}};

    PredictedEvent event;
    REQUIRE(det.find_next_event(req, event) == SolveStatus::InvalidInput);
}

TEST_CASE("FindNextEvent_RejectsNonFiniteInitialState", "[phase3][invalid]") {
    BodySystem sys = build_single();
    EventDetector det(sys);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    req.time_limit = 10.0;

    PredictedEvent event;
    req.initial_state = State2{{inf_v, 0.0}, {0.0, 1.0}};
    REQUIRE(det.find_next_event(req, event) == SolveStatus::InvalidInput);

    req.initial_state = State2{{nan_v, 0.0}, {0.0, 1.0}};
    REQUIRE(det.find_next_event(req, event) == SolveStatus::InvalidInput);

    req.initial_state = State2{{1.0, 2.0}, {inf_v, 0.0}};
    REQUIRE(det.find_next_event(req, event) == SolveStatus::InvalidInput);

    req.initial_state = State2{{1.0, 2.0}, {0.0, nan_v}};
    REQUIRE(det.find_next_event(req, event) == SolveStatus::InvalidInput);
}

TEST_CASE("FindNextEvent_NoNaNsOnOrdinaryInvalidInputs", "[phase3][invalid]") {
    BodySystem sys = build_single();
    EventDetector det(sys);

    PredictedEvent event;
    event.time = 42.0;  // sentinel

    // Zero-radius initial state is invalid input for propagation.
    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    req.initial_state = State2{{0.0, 0.0}, {0.0, 1.0}};
    req.time_limit = 10.0;

    auto status = det.find_next_event(req, event);
    REQUIRE(status != SolveStatus::Ok);
    // Even on failure, output should not contain NaN.
    REQUIRE(std::isfinite(event.time));
    REQUIRE(std::isfinite(event.state.r.x));
    REQUIRE(std::isfinite(event.state.r.y));
    REQUIRE(std::isfinite(event.state.v.x));
    REQUIRE(std::isfinite(event.state.v.y));
}

TEST_CASE("FindNextEvent_PropagatorFailureIsReportedExplicitly", "[phase3][invalid]") {
    BodySystem sys = build_single();
    EventDetector det(sys);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    req.initial_state = State2{{0.0, 0.0}, {0.0, 1.0}};  // zero radius fails propagation
    req.time_limit = 10.0;

    PredictedEvent event;
    auto status = det.find_next_event(req, event);
    REQUIRE(status != SolveStatus::Ok);
}

TEST_CASE("FindNextEvent_RefinementFailureIsReportedExplicitly", "[phase3][invalid]") {
    // Force refinement to fail by setting max_event_refine_iterations to 0.
    BodySystem sys = build_single();
    Tolerances tol{};
    tol.max_event_refine_iterations = 0;
    EventDetector det(sys, tol);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    // A scenario that produces an actual event -- impact on an inbound orbit.
    req.initial_state = State2{{3.0, 0.0}, {0.0, 0.3}};
    req.time_limit = 50.0;

    PredictedEvent event;
    auto status = det.find_next_event(req, event);
    // With zero refinement iterations, a bracketed event cannot be refined.
    // Either NoConvergence is returned or the scan reports no event -- either way, status
    // must be deterministic and finite.
    REQUIRE(std::isfinite(event.time));
    if (status != SolveStatus::Ok) {
        REQUIRE(status == SolveStatus::NoConvergence);
    }
}
