#include "brahe/trajectory_builder.h"

#include "brahe/event_detector.h"
#include "brahe/patcher.h"
#include "brahe/two_body.h"
#include "brahe/vec2.h"

#include <algorithm>
#include <cmath>

namespace brahe {

namespace {

bool is_finite_state(const State2& s) {
    return std::isfinite(s.r.x) && std::isfinite(s.r.y) &&
           std::isfinite(s.v.x) && std::isfinite(s.v.y);
}

// Internal post-patch suppression record. Prevents immediate re-detection of
// the SOI boundary just crossed due to floating-point noise in the patched
// state. See PHASE_4_PLAN.md section 6.6-6.7.
struct SuppressionState {
    EventType suppressed_event_type = EventType::None;  // SoiEntry or SoiExit
    BodyId suppressed_body = InvalidBody;
    double patch_time = 0.0;
    double boundary_radius = 0.0;
    bool active = false;
};

// Return true if `event` corresponds to the boundary that suppression is
// currently blocking. `central_body` is the central body the search is
// running in at this moment.
bool event_matches_suppressed_boundary(const PredictedEvent& event, BodyId central_body,
                                       const SuppressionState& s) {
    if (!s.active) return false;
    if (s.suppressed_event_type == EventType::SoiExit &&
        event.type == EventType::SoiExit &&
        central_body == s.suppressed_body) {
        return true;
    }
    if (s.suppressed_event_type == EventType::SoiEntry &&
        event.type == EventType::SoiEntry &&
        event.to_body == s.suppressed_body) {
        return true;
    }
    return false;
}

// Suppression is no longer in effect if EITHER the time or distance criterion
// has been met (PHASE_4_PLAN 6.7).
bool suppression_has_expired(const PredictedEvent& event, const BodySystem& bodies,
                             BodyId central_body, const SuppressionState& s,
                             const Tolerances& tol) {
    // Time criterion.
    double dt = event.time - s.patch_time;
    if (dt >= tol.time_epsilon) return true;

    // Distance criterion: |distance to suppressed boundary| >= 10 * pos_eps.
    double dist_from_boundary = 0.0;
    if (s.suppressed_event_type == EventType::SoiExit) {
        // Currently inside child frame; suppressed body IS the central body.
        // distance-from-boundary = |r_sc| - soi.
        double r = length(event.state.r);
        dist_from_boundary = std::abs(r - s.boundary_radius);
    } else {
        // Currently in parent frame; suppressed body is a child of central_body.
        Vec2 r_child = bodies.position_in_parent(s.suppressed_body, event.time);
        Vec2 d = event.state.r - r_child;
        double dist = length(d);
        dist_from_boundary = std::abs(dist - s.boundary_radius);
        (void)central_body;
    }
    if (dist_from_boundary >= 10.0 * tol.position_epsilon) return true;

    return false;
}

// Find the next non-suppressed event. If the detector returns an event that
// matches the currently suppressed boundary AND the suppression has not yet
// expired, advance the search request past the suppression window (propagate
// state by 2*time_epsilon) and try again. The advancement is bounded to a few
// retries to prevent infinite loops on pathological inputs.
SolveStatus find_next_event_respecting_suppression(
    const EventDetector& detector, const BodySystem& bodies,
    EventSearchRequest request, const SuppressionState& suppression,
    const Tolerances& tol, PredictedEvent& out_event) {
    constexpr int kMaxRetries = 4;
    for (int retry = 0; retry <= kMaxRetries; ++retry) {
        SolveStatus ds = detector.find_next_event(request, out_event);
        if (ds != SolveStatus::Ok) return ds;

        if (!event_matches_suppressed_boundary(out_event, request.central_body,
                                               suppression)) {
            return SolveStatus::Ok;
        }
        if (suppression_has_expired(out_event, bodies, request.central_body, suppression,
                                    tol)) {
            return SolveStatus::Ok;
        }

        // Event is suppressed -- advance the search past the suppression window.
        double new_start = suppression.patch_time + 2.0 * tol.time_epsilon;
        if (new_start >= request.time_limit) {
            // Suppression window extends past the horizon -- synthesize a
            // TimeLimit event at the horizon end. Use a best-effort propagated
            // final state (fall back to the current initial_state if the
            // propagator fails on a degenerate case).
            const BodyDef* cb = bodies.get_body(request.central_body);
            State2 final_state = request.initial_state;
            if (cb != nullptr) {
                State2 tmp{};
                if (TwoBody::propagate(cb->mu, request.initial_state,
                                       request.time_limit - request.start_time,
                                       tmp) == SolveStatus::Ok &&
                    is_finite_state(tmp)) {
                    final_state = tmp;
                }
            }
            out_event.type = EventType::TimeLimit;
            out_event.time = request.time_limit;
            out_event.from_body = request.central_body;
            out_event.to_body = InvalidBody;
            out_event.state = final_state;
            return SolveStatus::Ok;
        }

        // Propagate the search state forward to new_start.
        const BodyDef* cb = bodies.get_body(request.central_body);
        if (cb == nullptr) return SolveStatus::InvalidInput;
        State2 advanced{};
        SolveStatus ps = TwoBody::propagate(cb->mu, request.initial_state,
                                            new_start - request.start_time, advanced);
        if (ps != SolveStatus::Ok || !is_finite_state(advanced)) {
            return (ps != SolveStatus::Ok) ? ps : SolveStatus::NumericalFailure;
        }
        request.start_time = new_start;
        request.initial_state = advanced;
    }
    return SolveStatus::NumericalFailure;  // exceeded bounded retries
}

}  // namespace

TrajectoryBuilder::TrajectoryBuilder(const BodySystem& bodies, const Tolerances& tolerances)
    : bodies_(bodies), tolerances_(tolerances) {}

SolveStatus TrajectoryBuilder::build_preview(const PreviewRequest& req,
                                              Trajectory& out_trajectory) const {
    out_trajectory.segments.clear();

    // Bound-free input path: write directly into the output vector's storage
    // while still sharing the fixed-size append helper.
    if (req.max_segments == 0) return SolveStatus::InvalidInput;
    out_trajectory.segments.resize(req.max_segments);
    size_t count = 0;
    SolveStatus status =
        build_preview_into(req, out_trajectory.segments.data(), req.max_segments, count);
    out_trajectory.segments.resize(count);
    return status;
}

SolveStatus TrajectoryBuilder::build_preview_into(const PreviewRequest& req,
                                                   Segment* out_segments, size_t capacity,
                                                   size_t& out_count) const {
    out_count = 0;

    // --- Input validation ---
    if (!std::isfinite(req.start_time) || !std::isfinite(req.end_time)) {
        return SolveStatus::InvalidInput;
    }
    if (req.end_time < req.start_time) return SolveStatus::InvalidInput;
    if (!is_finite_state(req.initial_state)) return SolveStatus::InvalidInput;
    if (bodies_.get_body(req.central_body) == nullptr) return SolveStatus::InvalidInput;
    if (req.max_segments == 0) return SolveStatus::InvalidInput;

    const size_t max_seg = std::min(capacity, req.max_segments);

    // Empty horizon: no segments (PHASE_4_PLAN 10.3, preferred policy).
    if (req.end_time <= req.start_time + tolerances_.time_epsilon) {
        return SolveStatus::Ok;
    }

    EventDetector detector(bodies_, tolerances_);
    Patcher patcher(bodies_, tolerances_);

    BodyId cur_cb = req.central_body;
    double cur_t = req.start_time;
    State2 cur_s = req.initial_state;
    SuppressionState suppression{};

    // Bound the outer loop defensively. A healthy chain terminates on Impact,
    // TimeLimit, or CapacityExceeded; this cap only guards against pathological
    // internal bugs that would otherwise spin forever.
    const size_t kOuterLoopCap = max_seg + 16;
    for (size_t iter = 0; iter < kOuterLoopCap; ++iter) {
        if (out_count >= max_seg) {
            return SolveStatus::CapacityExceeded;
        }

        EventSearchRequest er;
        er.central_body = cur_cb;
        er.start_time = cur_t;
        er.initial_state = cur_s;
        er.time_limit = req.end_time;

        PredictedEvent event{};
        SolveStatus ds = find_next_event_respecting_suppression(
            detector, bodies_, er, suppression, tolerances_, event);
        if (ds != SolveStatus::Ok) return ds;

        // Emit one segment for the time window we just searched.
        Segment& seg = out_segments[out_count];
        seg.central_body = cur_cb;
        seg.start_time = cur_t;
        seg.initial_state = cur_s;
        seg.end_time = event.time;
        seg.end_reason = event.type;
        ++out_count;

        // Terminating events end the chain.
        if (event.type == EventType::Impact || event.type == EventType::TimeLimit) {
            return SolveStatus::Ok;
        }
        if (event.type == EventType::None) {
            // Detector returned no event but Ok -- shouldn't happen.
            return SolveStatus::NumericalFailure;
        }

        // Patch for the next segment.
        if (event.type == EventType::SoiEntry) {
            const BodyDef* ch = bodies_.get_body(event.to_body);
            if (ch == nullptr) return SolveStatus::NumericalFailure;
            State2 patched{};
            SolveStatus ps = patcher.patch_parent_to_child(event.time, cur_cb,
                                                            event.to_body, event.state,
                                                            patched);
            if (ps != SolveStatus::Ok) return ps;

            cur_cb = event.to_body;
            cur_t = event.time;
            cur_s = patched;
            suppression = SuppressionState{EventType::SoiExit, event.to_body, event.time,
                                            ch->soi_radius, true};
        } else if (event.type == EventType::SoiExit) {
            const BodyDef* cb_def = bodies_.get_body(cur_cb);
            if (cb_def == nullptr) {
                return SolveStatus::NumericalFailure;
            }
            if (cb_def->parent_id == InvalidBody) {
                return SolveStatus::Ok;
            }
            State2 patched{};
            SolveStatus ps = patcher.patch_child_to_parent(event.time, cur_cb, event.state,
                                                           patched);
            if (ps != SolveStatus::Ok) return ps;

            BodyId exited_child = cur_cb;
            double exited_soi = cb_def->soi_radius;
            cur_cb = cb_def->parent_id;
            cur_t = event.time;
            cur_s = patched;
            suppression = SuppressionState{EventType::SoiEntry, exited_child, event.time,
                                            exited_soi, true};
        } else {
            return SolveStatus::NumericalFailure;
        }
    }
    return SolveStatus::NumericalFailure;
}

}  // namespace brahe
