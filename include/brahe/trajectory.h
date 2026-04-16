#pragma once

#include "brahe/types.h"

#include <array>
#include <cstddef>
#include <vector>

namespace brahe {

// --- Segment ---
// A single-central-body arc with a starting state and a terminating event.
// Segment is valid over [start_time, end_time] propagated from initial_state
// in the central_body inertial frame.
struct Segment {
    BodyId central_body = InvalidBody;
    double start_time = 0.0;
    State2 initial_state = State2{};
    double end_time = 0.0;
    EventType end_reason = EventType::None;
};

// --- Trajectory event ---
// Describes a terminating event on a segment, optionally including the post-patch
// state in the new central-body frame (for SoiEntry / SoiExit transitions).
struct TrajectoryEvent {
    EventType type = EventType::None;
    double time = 0.0;
    BodyId from_body = InvalidBody;
    BodyId to_body = InvalidBody;  // InvalidBody if not applicable
    State2 state_before = State2{};
    State2 state_after = State2{};  // post-patch state for SOI transitions
};

// --- Heap-backed chain output ---
struct Trajectory {
    std::vector<Segment> segments;
};

// --- Fixed-capacity chain output ---
// 'count' is the number of valid segments; the tail beyond 'count' is unspecified.
template <size_t MaxSegments>
struct TrajectoryFixed {
    std::array<Segment, MaxSegments> segments = {};
    size_t count = 0;
};

// --- Preview-build request ---
// No burns in Phase 4; the request is narrow.
struct PreviewRequest {
    BodyId central_body = InvalidBody;
    double start_time = 0.0;
    State2 initial_state = State2{};
    double end_time = 0.0;
    size_t max_segments = 0;
};

}  // namespace brahe
