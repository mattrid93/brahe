#include "brahe/patcher.h"

#include "brahe/vec2.h"

#include <cmath>

namespace brahe {

namespace {

bool is_finite_state(const State2& s) {
    return std::isfinite(s.r.x) && std::isfinite(s.r.y) &&
           std::isfinite(s.v.x) && std::isfinite(s.v.y);
}

}  // namespace

Patcher::Patcher(const BodySystem& bodies, const Tolerances& tolerances)
    : bodies_(bodies), tolerances_(tolerances) {}

SolveStatus Patcher::patch_parent_to_child(double event_time, BodyId parent_body,
                                            BodyId child_body,
                                            const State2& spacecraft_in_parent_frame,
                                            State2& out_spacecraft_in_child_frame) const {
    if (!std::isfinite(event_time)) return SolveStatus::InvalidInput;
    if (!is_finite_state(spacecraft_in_parent_frame)) return SolveStatus::InvalidInput;

    const BodyDef* parent = bodies_.get_body(parent_body);
    const BodyDef* child = bodies_.get_body(child_body);
    if (parent == nullptr || child == nullptr) return SolveStatus::InvalidInput;
    // child must be a *direct* child of parent_body.
    if (child->parent_id != parent_body) return SolveStatus::InvalidInput;

    Vec2 r_child = bodies_.position_in_parent(child_body, event_time);
    Vec2 v_child = bodies_.velocity_in_parent(child_body, event_time);

    out_spacecraft_in_child_frame.r = spacecraft_in_parent_frame.r - r_child;
    out_spacecraft_in_child_frame.v = spacecraft_in_parent_frame.v - v_child;
    return SolveStatus::Ok;
}

SolveStatus Patcher::patch_child_to_parent(double event_time, BodyId child_body,
                                            const State2& spacecraft_in_child_frame,
                                            State2& out_spacecraft_in_parent_frame) const {
    if (!std::isfinite(event_time)) return SolveStatus::InvalidInput;
    if (!is_finite_state(spacecraft_in_child_frame)) return SolveStatus::InvalidInput;

    const BodyDef* child = bodies_.get_body(child_body);
    if (child == nullptr) return SolveStatus::InvalidInput;
    // Root body has no parent -- child->parent patch is undefined.
    if (child->parent_id == InvalidBody) return SolveStatus::InvalidInput;
    // Parent must exist in the body system.
    if (bodies_.get_body(child->parent_id) == nullptr) return SolveStatus::InvalidInput;

    Vec2 r_child = bodies_.position_in_parent(child_body, event_time);
    Vec2 v_child = bodies_.velocity_in_parent(child_body, event_time);

    out_spacecraft_in_parent_frame.r = spacecraft_in_child_frame.r + r_child;
    out_spacecraft_in_parent_frame.v = spacecraft_in_child_frame.v + v_child;
    return SolveStatus::Ok;
}

}  // namespace brahe
