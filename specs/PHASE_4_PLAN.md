# Patched Conics Library — Phase 4 Spec
## Patching and segment chaining

This document defines **Phase 4** of the implementation plan for the 2D patched-conics library.

Phase 4 covers:

- `Segment`, `Trajectory`, and `TrajectoryFixed`
- parent → child patch generation
- child → parent patch generation
- state continuity checks across patch boundaries
- post-patch re-entry suppression
- segment chaining
- preview-chain construction from repeated event detection + patching

This phase does **not** include impulsive burns, `rebuild_from(T)`, Lambert solving, transfer planning, or multiplayer snapshot handling. Those belong to later phases.

---

## 1. Development approach

For C++, use this sequence:

1. **write the public headers first**
2. **write tests against those headers**
3. **write implementation to satisfy the tests**

That is still TDD. The headers define:

- ownership and mutability
- output layout and serialization expectations
- fixed-capacity behavior
- what is and is not stored per segment
- how patch state is represented
- how preview building terminates
- what happens on capacity exhaustion or invalid geometry

This phase is particularly sensitive to API shape because it becomes the core substrate for:

- preview rendering
- maneuver editing
- runtime query mode
- later snapshot reconciliation

The correct order is therefore:

1. public API/header draft
2. compile-only API tests
3. patch math/unit tests
4. segment/preview behavior tests
5. property/fuzz tests
6. implementation
7. determinism hardening

---

## 2. Goals of Phase 4

By the end of this phase, the library shall be able to:

- represent a trajectory as an ordered chain of conic segments
- terminate a segment at a detected event
- generate a new segment by patching between parent and child frames
- preserve inertial position and velocity continuity across patch boundaries within configured tolerances
- suppress immediate re-detection of the SOI boundary that was just crossed
- repeatedly chain segments until impact, time limit, or segment cap
- support both heap-backed and fixed-capacity trajectory outputs
- operate deterministically on the same platform for identical inputs

Phase 4 should end with a standalone preview-chain builder that can take:

- an initial central body
- an initial state
- a start time
- a horizon / end time
- a max-segment limit

and produce the future conic chain with patch events applied.

---

## 3. Dependencies and assumed completed work

Phase 4 assumes the prior phases are complete.

Required from earlier phases:

### From Phase 1
- POD core types and enums
- immutable `BodySystem`
- circular ephemerides
- frame transforms
- deterministic body ordering
- tolerance defaults

### From Phase 2
- two-body propagation
- conic classification
- state-at-time propagation in a central-body frame

### From Phase 3
- event detection kernel
- impact detection
- SOI exit detection
- child SOI entry detection
- deterministic event arbitration
- bounded refinement
- grazing handling

If the event detector does not return enough information to distinguish:

- current central body
- target child body for SOI entry
- event time
- spacecraft state at event time

then stop and refine Phase 3 first. Phase 4 depends on that data being available and well-defined.

---

## 4. Scope of this phase

### In scope
- public trajectory/segment types
- patch generation math
- parent→child and child→parent frame changes
- end reason / event recording
- re-entry suppression state
- segment append logic
- preview chain construction to a horizon
- fixed-capacity and heap-backed trajectory outputs
- deterministic termination and tie behavior

### Out of scope
- burns
- `rebuild_from(T)`
- editing an existing chain
- Lambert
- transfer planning
- snapshot reconciliation
- optimization or planner heuristics
- UI concerns beyond data layout

---

## 5. Public API to lock before tests

The exact API can vary, but Phase 4 should freeze something close to the following.

## 5.1 Segment and trajectory types

```cpp
struct Segment {
    BodyId central_body;
    double start_time;
    State2 initial_state;
    double end_time;
    EventType end_reason;
};

struct TrajectoryEvent {
    EventType type;
    double time;
    BodyId from_body;
    BodyId to_body;     // InvalidBody if not applicable
    State2 state_before;
    State2 state_after; // post-patch state for SOI transitions
};

struct Trajectory {
    std::vector<Segment> segments;
};

template <size_t MaxSegments>
struct TrajectoryFixed {
    std::array<Segment, MaxSegments> segments;
    size_t count;
};
```

This is already consistent with the main spec. Phase 4 may need a small amount of extra metadata during construction, but that metadata should stay internal if possible.

## 5.2 Preview-build request

Recommended:

```cpp
struct PreviewRequest {
    BodyId central_body;
    double start_time;
    State2 initial_state;
    double end_time;
    size_t max_segments;
};
```

This phase should not include burns, so the request remains narrow.

## 5.3 Preview-build API

Recommended dynamic and fixed-capacity entry points:

```cpp
class TrajectoryBuilder {
public:
    TrajectoryBuilder(const BodySystem& bodies, const Tolerances& tolerances = {});

    SolveStatus build_preview(const PreviewRequest& req,
                              Trajectory& out_trajectory) const;

    template <size_t MaxSegments>
    SolveStatus build_preview_fixed(const PreviewRequest& req,
                                    TrajectoryFixed<MaxSegments>& out_trajectory) const;
};
```

## 5.4 Patch-generation API

Patch generation should be testable independently from full preview chaining.

Recommended public or internal-testable helper surface:

```cpp
class Patcher {
public:
    Patcher(const BodySystem& bodies, const Tolerances& tolerances = {});

    SolveStatus patch_parent_to_child(double event_time,
                                      BodyId parent_body,
                                      BodyId child_body,
                                      const State2& spacecraft_in_parent_frame,
                                      State2& out_spacecraft_in_child_frame) const;

    SolveStatus patch_child_to_parent(double event_time,
                                      BodyId child_body,
                                      const State2& spacecraft_in_child_frame,
                                      State2& out_spacecraft_in_parent_frame) const;
};
```

This API is valuable because it isolates the frame-transition math from the chain builder.

## 5.5 Optional internal suppression metadata

The spec requires post-patch re-entry suppression, but the public segment layout does not include suppression metadata. That is acceptable. Suppression can remain internal to the chain builder.

Recommended internal struct:

```cpp
struct SuppressionState {
    EventType suppressed_event_type;   // SoiEntry or SoiExit
    BodyId suppressed_body;
    double patch_time;
    double boundary_radius;
    bool active;
};
```

Only the builder needs this.

---

## 6. Behavioral rules to specify before coding

These rules should be written down before implementation.

## 6.1 Segment meaning

A `Segment` represents:

- one central body
- one initial state in that central-body inertial frame
- one start time
- one end time
- one end reason

A segment is valid over:

```cpp
[start_time, end_time]
```

in the sense that:
- the initial state is defined at `start_time`
- propagation within the segment frame is valid until the terminating event at `end_time`

`end_reason` must be one of:
- `Impact`
- `SoiEntry`
- `SoiExit`
- `TimeLimit`

Phase 4 does not use `Burn`.

## 6.2 Segment continuity rules

### Across ordinary within-frame propagation
Propagation inside one segment is continuous in position and velocity by definition.

### Across patch boundaries
At a parent↔child or child↔parent patch:

- inertial position must be continuous within `position_epsilon`
- inertial velocity must be continuous within `velocity_epsilon`

Acceleration is allowed to be discontinuous.

This must be checked in tests, not just assumed from the formulas.

## 6.3 Parent → child patch rule

At event time `t_e`:

1. compute spacecraft absolute state in the parent inertial frame
2. compute child absolute state in the same parent inertial frame
3. subtract child state from spacecraft state
4. resulting state is the new segment initial state in the child-centered inertial frame

Formula:

```cpp
r_sc_child = r_sc_parent_abs - r_child_parent_abs
v_sc_child = v_sc_parent_abs - v_child_parent_abs
```

The new segment central body becomes `child_body`.

## 6.4 Child → parent patch rule

At event time `t_e`:

1. compute spacecraft relative state in the child frame
2. compute child absolute state in the parent inertial frame
3. add child state to spacecraft state
4. resulting state is the new segment initial state in the parent-centered inertial frame

Formula:

```cpp
r_sc_parent = r_sc_child + r_child_parent_abs
v_sc_parent = v_sc_child + v_child_parent_abs
```

The new segment central body becomes the child's parent.

## 6.5 Preview-chain construction rule

A preview chain is built by repeating:

1. start from current segment initial state
2. ask event detector for earliest event before preview end time
3. terminate current segment at that event time
4. if the event is:
   - `Impact`: stop
   - `TimeLimit`: stop
   - `SoiEntry`: patch parent→child and continue
   - `SoiExit`: patch child→parent and continue
5. repeat until stop condition or segment cap

## 6.6 Re-entry suppression rule

After a patch, the new segment shall suppress detection of the same SOI boundary just crossed until:

- at least `time_epsilon` has elapsed **or**
- the spacecraft has moved at least `position_epsilon * 10` away from that boundary,

whichever occurs first.

This rule is required to prevent immediate re-patch loops caused by floating-point subtraction and exact-boundary starts.

Important clarification:

- after parent→child patch, suppress the child SOI boundary that was just entered
- after child→parent patch, suppress the child SOI boundary that was just exited

Suppression applies only to the crossed boundary. It must not suppress:
- impact detection
- unrelated sibling child entries
- unrelated future SOI exits/entries

## 6.7 Suppression deactivation rule

Suppression becomes inactive as soon as either condition is met:

```cpp
(current_time - patch_time) >= time_epsilon
```

or

```cpp
abs(distance_to_suppressed_boundary) >= position_epsilon * 10
```

where `distance_to_suppressed_boundary` is the signed boundary function value in the new frame.

## 6.8 Capacity behavior

For `TrajectoryFixed<N>`:

- if the preview would require more than `N` segments, return an explicit failure status
- do not write beyond capacity
- the written prefix must remain valid and deterministic

Prefer adding a dedicated status:

```cpp
enum class SolveStatus {
    Ok,
    InvalidInput,
    NoConvergence,
    NoSolution,
    NumericalFailure,
    CapacityExceeded
};
```

If the enum is already frozen, document which existing status is used for segment-cap exhaustion.

## 6.9 Time ordering invariants

For every output chain:

- segment start times are monotonic nondecreasing
- each segment has `end_time >= start_time`
- for consecutive segments, next `start_time == previous end_time` within `time_epsilon`
- no infinite loops of zero-duration segments are allowed

A zero-duration final segment should generally be avoided. If one can occur from exact-boundary/time-limit coincidence, define the policy explicitly and test it.

Recommended:
- do not emit a trailing zero-duration continuation segment after a terminating event

---

## 7. TDD plan

Write headers first, then all of the following tests before implementation.

## 7.1 Compile-only API tests

These should be added first.

1. **`TrajectoryBuilder_CanBeConstructedFromBodySystem`**
2. **`BuildPreview_AcceptsPreviewRequestAndDynamicTrajectory`**
3. **`BuildPreviewFixed_AcceptsPreviewRequestAndFixedTrajectory`**
4. **`Patcher_CanBeConstructedFromBodySystem`**
5. **`PatchParentToChild_HasStableSignature`**
6. **`PatchChildToParent_HasStableSignature`**
7. **`Segment_IsTriviallyCopyable`**
8. **`Segment_IsStandardLayout`**
9. **`TrajectoryEvent_IsTriviallyCopyable`**
10. **`TrajectoryEvent_IsStandardLayout`**

These do not need runtime assertions.

---

## 7.2 Patch-math unit tests

Write these before any preview-builder logic.

Use simple axis-aligned systems so expected vectors are obvious.

### Parent → child patch tests

11. **`PatchParentToChild_SubtractsChildPositionAndVelocity`**
   - Use a known child state and known spacecraft parent-frame state.
   - Verify exact vector subtraction.

12. **`PatchParentToChild_ChangesCentralBodyToChildSemantically`**
   - This may be tested through builder output rather than the patch helper itself.

13. **`PatchParentToChild_PreservesAbsoluteInertialState`**
   - Convert the patched child-frame state back into parent-frame absolute coordinates.
   - Compare with original spacecraft absolute parent-frame state.

### Child → parent patch tests

14. **`PatchChildToParent_AddsChildPositionAndVelocity`**
15. **`PatchChildToParent_PreservesAbsoluteInertialState`**

### Round-trip tests

16. **`PatchParentToChildThenBackToParent_PreservesStateWithinTolerance`**
17. **`PatchChildToParentThenBackToChild_PreservesStateWithinTolerance`**

### Numeric guard tests

18. **`PatchFunctions_RejectUnknownBodies`**
19. **`PatchFunctions_RejectInvalidParentChildRelationship`**
20. **`PatchFunctions_RejectNonFiniteTimeOrState`**
21. **`PatchFunctions_DoNotEmitNaNOnOrdinaryInvalidInput`**

---

## 7.3 State continuity tests at patch boundaries

These are core correctness tests and must be explicit.

22. **`PatchBoundary_PositionContinuityWithinPositionEpsilon`**
23. **`PatchBoundary_VelocityContinuityWithinVelocityEpsilon`**

Use a realistic parent/child ephemeris plus a detected event-time state.

These tests should compare:
- spacecraft absolute inertial state before patch
- spacecraft absolute inertial state reconstructed from the post-patch state

---

## 7.4 Re-entry suppression unit tests

These should be written before preview-chain implementation if suppression logic is factored cleanly.

24. **`Suppression_ActivatesImmediatelyAfterParentToChildPatch`**
25. **`Suppression_ActivatesImmediatelyAfterChildToParentPatch`**
26. **`Suppression_BlocksOnlyTheCrossedBoundary`**
27. **`Suppression_ExpiresAfterTimeEpsilonEvenIfStillNearBoundary`**
28. **`Suppression_ExpiresAfterMovingAwayByPositionThresholdEvenBeforeTimeEpsilon`**
29. **`Suppression_DoesNotBlockImpactDetection`**
30. **`Suppression_DoesNotBlockUnrelatedSiblingChildEntry`**

If suppression is purely internal, test it through preview scenarios.

---

## 7.5 Segment construction tests

These validate segment emission behavior separate from long chains.

31. **`SingleSegment_EndsAtImpactEvent`**
32. **`SingleSegment_EndsAtSoiEntryEvent`**
33. **`SingleSegment_EndsAtSoiExitEvent`**
34. **`SingleSegment_EndsAtTimeLimitWhenNoEarlierEventOccurs`**
35. **`SingleSegment_StartAndEndTimesAreCorrect`**
36. **`SingleSegment_StoresInitialStateExactly`**

---

## 7.6 Segment chaining tests

These are the central Phase 4 tests.

37. **`PreviewChain_StopsOnImpact`**
38. **`PreviewChain_StopsOnTimeLimit`**
39. **`PreviewChain_PatchesIntoChildAfterSoiEntry`**
40. **`PreviewChain_PatchesBackToParentAfterSoiExit`**
41. **`PreviewChain_BuildsMultipleSegmentsAcrossParentChildParentSequence`**
42. **`PreviewChain_UsesPatchedStateAsNextSegmentInitialState`**
43. **`PreviewChain_DoesNotImmediatelyRepatchAfterBoundaryCrossing`**
44. **`PreviewChain_SegmentTimesAreMonotonic`**
45. **`PreviewChain_AdjacentSegmentsMeetAtSameTimeWithinTimeEpsilon`**
46. **`PreviewChain_NoInfiniteLoopAtBoundary`**

---

## 7.7 Fixed-capacity tests

47. **`FixedCapacity_SucceedsWhenSegmentCountFits`**
48. **`FixedCapacity_ReturnsCapacityFailureWhenSegmentCountWouldOverflow`**
49. **`FixedCapacity_DoesNotWritePastArrayBounds`**
50. **`FixedCapacity_PreservesValidPrefixOnFailure`**

---

## 7.8 Determinism tests

51. **`BuildPreview_RepeatedRunsProduceBitIdenticalSegmentsOnSamePlatform`**
52. **`BuildPreview_IsIndependentOfBodyInsertionOrder`**
53. **`ChildEntryTieBreak_RemainsStableAcrossRepeatedRuns`**
54. **`SuppressionBehavior_IsStableAcrossRepeatedRuns`**

These should use checked binary comparisons of POD segment outputs where possible.

---

## 7.9 Property tests

55. **`PatchRoundTrip_PreservesAbsoluteInertialState`**
56. **`PreviewChain_EachNextSegmentBeginsFromPreviousTerminalEventState`**
57. **`PreviewChain_FinalSegmentAlwaysEndsAtImpactOrTimeLimitOrPatchEvent`**
58. **`PreviewChain_NoNaNsInAnySegmentEndpoints`**
59. **`DynamicAndFixedPreviewAgreeUpToCapacity`**

---

## 7.10 Fuzz tests

60. **`RandomValidSystemsAndStates_NoCrashesNoNaNs`**
61. **`RandomNestedParentChildTransitions_NoImmediateRepatchLoops`**
62. **`RandomPreviewBuilds_AllSegmentBodiesAreValid`**
63. **`RandomPreviewBuilds_AllSegmentTimesAreOrdered`**

---

## 8. Implementation plan after tests exist

Implement in thin slices.

## 8.1 Slice 1: headers only

Write or update:

- `trajectory.h`
- `patcher.h`
- `trajectory_builder.h`

Lock:
- public types
- status returns
- fixed-capacity semantics
- preview request shape

Add the compile-only tests first.

## 8.2 Slice 2: patch math only

Implement:
- parent→child patch
- child→parent patch

Pass tests 11–21.

## 8.3 Slice 3: continuity validation

Implement or at least test helper logic that reconstructs inertial state before/after patch and verifies continuity.

Pass tests 22–23.

## 8.4 Slice 4: suppression helper

Implement internal suppression state and gating logic.

Pass tests 24–30.

## 8.5 Slice 5: single-segment builder

Implement logic that:
- asks event detector for next event
- emits one segment ending at that event

Pass tests 31–36.

## 8.6 Slice 6: multi-segment chaining

Implement repeated patch-and-continue behavior.

Pass tests 37–46.

## 8.7 Slice 7: fixed-capacity output

Implement `TrajectoryFixed<N>` path and overflow handling.

Pass tests 47–50.

## 8.8 Slice 8: determinism/property/fuzz hardening

Pass tests 51–63.

---

## 9. Internal implementation recommendations

These are implementation recommendations, not required public API commitments.

## 9.1 Separate patch math from chain orchestration

Keep these concerns isolated:

- event detection
- patch state generation
- suppression logic
- segment append logic
- full preview loop

This makes the tests smaller and the failure modes clearer.

## 9.2 Use one internal “segment step” operation

Recommended internal helper:

```cpp
struct SegmentStepResult {
    Segment segment;
    bool has_followup_segment;
    BodyId next_central_body;
    State2 next_initial_state;
    double next_start_time;
};
```

Or equivalent.

This should represent:
- one emitted segment
- whether patching continues
- the next segment seed if it does

This helps keep the preview loop simple.

## 9.3 Centralize suppression checks

Do not scatter suppression logic across multiple event paths.

Recommended internal helper:

```cpp
bool is_event_suppressed(const PredictedEvent& event,
                         const SuppressionState& suppression,
                         double current_time,
                         const State2& current_state,
                         const Tolerances& tol);
```

Every patch-induced suppression decision should pass through one function.

## 9.4 Avoid zero-length segment loops

After every chained step, assert one of:

- time advanced by more than `time_epsilon`, or
- the crossed-boundary event is now suppressed, or
- the builder returns explicit failure

This is essential to prevent oscillation at exact SOI boundaries.

## 9.5 Prefer event detector to stay authoritative for event timing

The preview builder should not recompute event times in ad hoc ways. It should:

- trust the event detector for event time
- use patch math only at that returned time
- not adjust event times except through explicit epsilon-policy handling

---

## 10. Edge cases to lock before coding

These should not be left vague.

## 10.1 Entry into child SOI exactly at preview start

Recommended:
- if the detector returns `SoiEntry` at `start_time`, emit a segment ending at `start_time` only if doing so is required to preserve event semantics
- otherwise patch immediately and start the first emitted segment in the child frame

Preferred policy:
- avoid emitting a meaningless zero-duration pre-patch segment
- patch immediately and begin the child-centered segment at `start_time`

This policy should be stated and tested.

## 10.2 Exit from child SOI exactly at preview start

Use the same policy:
- patch immediately to the parent frame
- do not emit a zero-duration pre-exit segment

## 10.3 Time limit equal to start time

Recommended:
- emit no segments and return `Ok`, or
- emit one zero-duration segment ending at `TimeLimit`

Preferred policy:
- emit no segments for empty horizon
- keep output empty and deterministic

Lock and test this.

## 10.4 Impact exactly at patch boundary time

Event ordering from Phase 3 should already resolve this:
- impact beats SOI entry/exit within `time_epsilon`

The chain builder must not override that.

## 10.5 Overlapping sibling SOIs

If the detector resolves simultaneous sibling entry by lowest `BodyId`, the chain builder must accept that result exactly and patch to that child only.

---

## 11. Failure handling requirements

Every public Phase 4 API must:

- return explicit `SolveStatus`
- fully initialize written outputs
- never expose uninitialized segment contents
- never emit NaNs on ordinary invalid input
- never loop indefinitely

Required failure categories:

- invalid preview request
- invalid parent/child relationship during patching
- event detector failure
- numerical failure during propagation or patch reconstruction
- capacity exhaustion

Output rules:

### Dynamic trajectory output
Recommended:
- clear the output before build
- on failure, preserve the valid prefix already written
- document that the prefix is valid up to the returned failure

### Fixed-capacity output
- `count` must always reflect the valid prefix length
- the array beyond `count` is unspecified
- no out-of-bounds writes are allowed

---

## 12. Determinism requirements specific to Phase 4

Phase 4 adds several determinism hazards:

- repeated parent/child transitions
- suppression activation/deactivation boundaries
- segment-cap termination
- fixed vs dynamic output path drift
- boundary-time zero-length segment policies

The implementation shall avoid:

- unordered child traversal
- unstable segment append order
- inconsistent epsilon checks for chaining
- hidden dependence on vector capacity growth or pointer values
- special-case suppression logic distributed across multiple code paths

Required same-platform properties:

- same body system + same request + same tolerances -> bit-identical segment outputs
- fixed and dynamic trajectories agree on overlapping prefix semantics
- patch round-trips are stable across repeated runs

Recommended tests:
- build the same preview twice and compare raw bytes of segments
- build with body insertion orders permuted and compare outputs
- compare fixed-capacity and dynamic outputs where the full chain fits

---

## 13. Suggested file layout for this phase

```text
include/patched_conics/trajectory.h
include/patched_conics/patcher.h
include/patched_conics/trajectory_builder.h

src/patcher.cpp
src/trajectory_builder.cpp

tests/test_trajectory_api.cpp
tests/test_patcher.cpp
tests/test_patch_continuity.cpp
tests/test_patch_suppression.cpp
tests/test_single_segment.cpp
tests/test_segment_chaining.cpp
tests/test_trajectory_fixed.cpp
tests/test_phase4_determinism.cpp
tests/test_phase4_properties.cpp
tests/test_phase4_fuzz.cpp
```

If the codebase is still small, fewer source files are fine. The important part is conceptual separation.

---

## 14. Definition of done for Phase 4

Phase 4 is complete when:

- the public patching/trajectory API is header-defined and stable
- tests exist before implementation for all major patching and chaining behaviors
- parent→child and child→parent patching work correctly
- inertial continuity is preserved across patch boundaries within tolerance
- post-patch suppression prevents immediate re-patch loops
- preview chains can be built through multiple SOI transitions
- fixed-capacity and dynamic outputs both behave correctly
- invalid inputs return explicit failure without crashes or NaNs
- repeated runs on the same platform produce bit-identical outputs

---

## 15. Recommended execution order summary

1. draft `trajectory.h`, `patcher.h`, and `trajectory_builder.h`
2. add compile-only API tests
3. write patch-math unit tests
4. write patch continuity tests
5. write suppression tests
6. write single-segment builder tests
7. write multi-segment chaining tests
8. write fixed-capacity tests
9. write determinism/property/fuzz tests
10. implement patch math
11. implement suppression helper
12. implement single-segment emission
13. implement repeated patch chaining
14. implement fixed-capacity path
15. harden edge cases and determinism

This phase should end with a robust trajectory-chaining core that later phases can extend with burns, rebuild, and planning.
