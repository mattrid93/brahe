#include "brahe/event_detector.h"

#include "brahe/two_body.h"
#include "brahe/vec2.h"

#include <algorithm>
#include <cmath>

namespace brahe {

// --- EventDetector ---

EventDetector::EventDetector(const BodySystem& bodies, const Tolerances& tolerances)
    : bodies_(bodies), tolerances_(tolerances) {}

SolveStatus EventDetector::find_next_event(const EventSearchRequest& /*req*/,
                                           PredictedEvent& out_event) const {
    out_event = PredictedEvent{};
    return SolveStatus::InvalidInput;  // stub
}

// --- detail helpers ---

namespace detail {

bool event_precedes(const PredictedEvent& a, const PredictedEvent& b, double time_epsilon) {
    // If times differ by more than epsilon, earlier wins
    if (a.time < b.time - time_epsilon) return true;
    if (b.time < a.time - time_epsilon) return false;

    // Within epsilon: use type priority (lower enum value = higher priority)
    // Impact(1) > SoiEntry(2) > SoiExit(3) > Burn(4) > TimeLimit(5)
    if (a.type != b.type) return static_cast<int>(a.type) < static_cast<int>(b.type);

    // Same type within epsilon: lowest BodyId wins (for sibling SOI entries)
    return a.to_body < b.to_body;
}

}  // namespace detail

}  // namespace brahe
