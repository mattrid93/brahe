#pragma once

#include "brahe/body_system.h"
#include "brahe/trajectory.h"
#include "brahe/types.h"

#include <cstddef>

namespace brahe {

// TrajectoryBuilder assembles a chain of conic segments forward in time by
// repeatedly asking the event detector for the next event, terminating the
// current segment at that event, and (for SoiEntry / SoiExit) patching into
// the appropriate new central-body frame. See PHASE_4_PLAN.md section 6 for
// the full behavioral contract.
//
// Termination:
//   - Impact or TimeLimit end the chain.
//   - Reaching max_segments yields CapacityExceeded (valid prefix preserved).
//   - InvalidInput or downstream failures short-circuit with an explicit
//     SolveStatus and leave the output cleared of trailing unspecified data.
class TrajectoryBuilder {
public:
    explicit TrajectoryBuilder(const BodySystem& bodies, const Tolerances& tolerances = {});

    // Heap-backed output. out_trajectory.segments is cleared before building.
    SolveStatus build_preview(const PreviewRequest& req, Trajectory& out_trajectory) const;

    // Fixed-capacity output. out_trajectory.count is set to the prefix length
    // that was actually written. On CapacityExceeded, the prefix is valid but
    // truncated — the final segment's end_reason reflects the last terminating
    // event processed before the cap was hit.
    template <size_t MaxSegments>
    SolveStatus build_preview_fixed(const PreviewRequest& req,
                                    TrajectoryFixed<MaxSegments>& out_trajectory) const {
        out_trajectory.count = 0;
        return build_preview_into(req, out_trajectory.segments.data(), MaxSegments,
                                  out_trajectory.count);
    }

private:
    // Shared append-into-buffer implementation. out_count is always set to the
    // number of valid segments actually written (the prefix length), even on
    // failure returns.
    SolveStatus build_preview_into(const PreviewRequest& req, Segment* out_segments,
                                   size_t capacity, size_t& out_count) const;

    const BodySystem& bodies_;
    Tolerances tolerances_;
};

}  // namespace brahe
