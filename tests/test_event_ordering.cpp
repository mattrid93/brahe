#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "brahe/body_system.h"
#include "brahe/event_detector.h"

#include <cmath>

using namespace brahe;
using Catch::Matchers::WithinAbs;

// --- Section 7.8: Event-ordering tests ---
// Priority (spec 9.5.1, phase3 6.4):
//   Impact > SoiEntry > SoiExit > Burn > TimeLimit
// With lowest-BodyId tiebreak for sibling SoiEntries.

namespace {

PredictedEvent make_event(EventType t, double time, BodyId to_body = InvalidBody) {
    PredictedEvent e;
    e.type = t;
    e.time = time;
    e.from_body = 1;
    e.to_body = to_body;
    return e;
}

}  // namespace

TEST_CASE("EventOrdering_ImpactBeatsChildEntryWithinTimeEpsilon", "[phase3][ordering]") {
    double eps = 1e-6;
    auto impact = make_event(EventType::Impact, 1.0);
    auto entry = make_event(EventType::SoiEntry, 1.0 + 0.5 * eps, 2);
    REQUIRE(detail::event_precedes(impact, entry, eps));
    REQUIRE_FALSE(detail::event_precedes(entry, impact, eps));
}

TEST_CASE("EventOrdering_ChildEntryBeatsSoiExitWithinTimeEpsilon", "[phase3][ordering]") {
    double eps = 1e-6;
    auto entry = make_event(EventType::SoiEntry, 5.0, 3);
    auto exit = make_event(EventType::SoiExit, 5.0 + 0.5 * eps);
    REQUIRE(detail::event_precedes(entry, exit, eps));
    REQUIRE_FALSE(detail::event_precedes(exit, entry, eps));
}

TEST_CASE("EventOrdering_SoiExitBeatsTimeLimitWithinTimeEpsilon", "[phase3][ordering]") {
    double eps = 1e-6;
    auto exit = make_event(EventType::SoiExit, 10.0);
    auto limit = make_event(EventType::TimeLimit, 10.0 + 0.5 * eps);
    REQUIRE(detail::event_precedes(exit, limit, eps));
    REQUIRE_FALSE(detail::event_precedes(limit, exit, eps));
}

TEST_CASE("EventOrdering_MultipleSiblingEntriesUseLowestBodyId", "[phase3][ordering]") {
    double eps = 1e-6;
    auto ea = make_event(EventType::SoiEntry, 3.0, 10);
    auto eb = make_event(EventType::SoiEntry, 3.0 + 0.5 * eps, 20);
    REQUIRE(detail::event_precedes(ea, eb, eps));
    REQUIRE_FALSE(detail::event_precedes(eb, ea, eps));

    // Reversed-order input with same body IDs must still select the lower BodyId.
    auto ec = make_event(EventType::SoiEntry, 3.0, 50);
    auto ed = make_event(EventType::SoiEntry, 3.0 + 0.5 * eps, 30);
    REQUIRE(detail::event_precedes(ed, ec, eps));
    REQUIRE_FALSE(detail::event_precedes(ec, ed, eps));
}

TEST_CASE("EventOrdering_StableAcrossRepeatedRuns", "[phase3][ordering]") {
    // Full detector: same inputs must produce bit-identical events.
    BodySystemBuilder b;
    BodyDef root{};
    root.id = 1;
    root.parent_id = InvalidBody;
    root.mu = 1.0;
    root.radius = 1.0;
    root.soi_radius = 100.0;
    b.add_body(root);
    BodySystem sys;
    REQUIRE(b.build(sys) == SolveStatus::Ok);

    EventDetector det(sys);

    EventSearchRequest req;
    req.central_body = 1;
    req.start_time = 0.0;
    req.initial_state = State2{{3.0, 0.0}, {0.0, 0.3}};  // elliptic with rp<1 -> impact
    req.time_limit = 50.0;

    PredictedEvent e1, e2;
    REQUIRE(det.find_next_event(req, e1) == SolveStatus::Ok);
    REQUIRE(det.find_next_event(req, e2) == SolveStatus::Ok);

    REQUIRE(e1.type == e2.type);
    REQUIRE(e1.time == e2.time);  // bit-identical
    REQUIRE(e1.from_body == e2.from_body);
    REQUIRE(e1.to_body == e2.to_body);
    REQUIRE(e1.state.r.x == e2.state.r.x);
    REQUIRE(e1.state.r.y == e2.state.r.y);
    REQUIRE(e1.state.v.x == e2.state.v.x);
    REQUIRE(e1.state.v.y == e2.state.v.y);
}
