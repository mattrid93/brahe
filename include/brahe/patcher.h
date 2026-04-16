#pragma once

#include "brahe/body_system.h"
#include "brahe/types.h"

namespace brahe {

// Frame-change math for SOI patch transitions.
//
// parent_to_child: given a spacecraft state in the parent-centered inertial
// frame at event_time, produce its state in the child-centered inertial frame.
// The child must be a direct child of parent_body in the body hierarchy.
//
// child_to_parent: given a spacecraft state in the child-centered inertial
// frame at event_time, produce its state in the parent-centered inertial
// frame. The parent is inferred from child_body's parent link.
//
// Both preserve inertial position and velocity (no rotation — 2D inertial
// frames are parallel). Inputs must be finite. Unknown bodies, invalid
// parent/child relationships, non-finite states, or non-finite event_time
// yield InvalidInput and leave the output untouched.
class Patcher {
public:
    explicit Patcher(const BodySystem& bodies, const Tolerances& tolerances = {});

    SolveStatus patch_parent_to_child(double event_time,
                                      BodyId parent_body,
                                      BodyId child_body,
                                      const State2& spacecraft_in_parent_frame,
                                      State2& out_spacecraft_in_child_frame) const;

    SolveStatus patch_child_to_parent(double event_time,
                                      BodyId child_body,
                                      const State2& spacecraft_in_child_frame,
                                      State2& out_spacecraft_in_parent_frame) const;

private:
    const BodySystem& bodies_;
    Tolerances tolerances_;
};

}  // namespace brahe
