# Patched Conics Library — Phase 3 Spec
## Event detection kernel

This document defines **Phase 3** of the implementation plan for the 2D patched-conics library.

Phase 3 covers the **event detection kernel** only:

- impact root functions
- SOI exit root functions
- child SOI entry search
- coarse scan + bounded refinement
- grazing handling
- deterministic event ordering and tie-breaks

This phase does **not** include patch application, segment chaining, impulsive burns, rebuild logic, Lambert solving, or transfer planning. Those belong to later phases.

---

## 1. Development approach

For C++, the right order here is:

1. **write the public headers first**
2. **write tests against those headers**
3. **write implementation to satisfy the tests**

That is still TDD. The headers are the API contract. Writing them first is useful because event detection is one of the most interface-sensitive parts of the library:

- it determines which inputs are explicit
- it defines what failure looks like
- it fixes deterministic tie-break rules
- it prevents implementation details from leaking into callers

The sequence for this phase should be:

1. public API/header draft
2. compile-only API tests
3. unit tests for root functions and refinement
4. scenario tests for SOI entry/exit/impact ordering
5. property and fuzz tests
6. implementation
7. determinism hardening

---

## 2. Goals of Phase 3

By the end of this phase, the library shall be able to:

- predict the earliest event after a given start time for a spacecraft moving on a single two-body conic around a current central body
- detect impact with the current central body
- detect exit from the current body's SOI
- detect entry into any child body's SOI
- detect grazing encounters that do not produce a sign change
- refine candidate event times with bounded scalar root solving
- apply deterministic event ordering when multiple events occur within `time_epsilon`
- return explicit status on invalid input or numerical failure
- do all of the above without emitting NaNs on ordinary invalid inputs

This phase produces the event detector that later phases will use for patching, preview-chain building, burns, and reconciliation.

---

## 3. Dependencies and assumed completed work

Phase 3 assumes **Phase 1** and **Phase 2** are complete.

Required from earlier phases:

### From Phase 1
- POD core types and enums
- immutable `BodySystem`
- body hierarchy and child lookup
- circular ephemerides
- frame transforms
- tolerance defaults

### From Phase 2
- two-body propagation
- conic classification
- state propagation from `(mu, initial_state, dt) -> State2`
- safe low-level vector math
- any anomaly stepping helpers that are already available

If Phase 2 did not expose a reusable propagation API, stop and fix that first. The event detector must depend on a stable propagate function, not on internal propagation details.

---

## 4. Scope of this phase

### In scope
- event root function definitions
- event candidate search over time
- coarse stepping / bracketing
- bounded root refinement
- local-minimum checks for grazing events
- child iteration and sibling tie-breaks
- deterministic earliest-event selection
- event result structs and failure statuses

### Out of scope
- applying the event to create a new segment
- patch state generation
- post-patch re-entry suppression
- burns
- preview chain building
- `rebuild_from(T)`
- Lambert
- transfer planning

---

## 5. Public API to lock before tests

The exact API can vary, but the event detector needs a clear and narrow public surface.

## 5.1 Event result types

```cpp
enum class EventType {
    None,
    Impact,
    SoiEntry,
    SoiExit,
    Burn,
    TimeLimit
};

struct PredictedEvent {
    EventType type;
    double time;
    BodyId from_body;
    BodyId to_body;      // child body for SoiEntry, parent for SoiExit, InvalidBody otherwise
    State2 state;        // spacecraft state in current central-body frame at event time
};
```

For this phase, `state` should mean the spacecraft state **just before** the event in the current central-body frame. Do not mix pre-event and post-patch state in this phase.

## 5.2 Detector request type

Recommended:

```cpp
struct EventSearchRequest {
    BodyId central_body;
    double start_time;
    State2 initial_state;    // state relative to central_body at start_time
    double time_limit;       // inclusive upper search bound
};
```

The detector should search for the earliest event in `[start_time, time_limit]`.

## 5.3 Detector API

Recommended:

```cpp
class EventDetector {
public:
    EventDetector(const BodySystem& bodies, const Tolerances& tolerances = {});

    SolveStatus find_next_event(const EventSearchRequest& req,
                                PredictedEvent& out_event) const;
};
```

## 5.4 Optional lower-level testable helpers

These helpers do not need to be public API, but they should exist as separable units inside the implementation and may be exposed in an internal header for tests.

Recommended internal helpers:

```cpp
double impact_function(const Vec2& r, double body_radius);
double soi_exit_function(const Vec2& r, double soi_radius);
double child_entry_function(const Vec2& spacecraft_abs,
                            const Vec2& child_abs,
                            double child_soi_radius);
```

Or time-parametric versions:

```cpp
SolveStatus eval_impact_f(double t, const Context&, double& out_f);
SolveStatus eval_soi_exit_f(double t, const Context&, double& out_f);
SolveStatus eval_child_entry_f(double t, const Context&, BodyId child, double& out_f);
```

The important point is separation. Root evaluation, bracketing, refinement, and candidate ordering should not be fused into one monolithic function.

---

## 6. Behavioral rules to specify before coding

These rules must be fixed before implementation.

## 6.1 Event types and meaning

### Impact
Occurs when the spacecraft distance from the current central body reaches the body's physical radius:

```cpp
|r_rel(t)| - body.radius = 0
```

### SOI exit
Occurs when the spacecraft distance from the current central body reaches the current body's SOI radius from the inside:

```cpp
|r_rel(t)| - body.soi_radius = 0
```

### Child SOI entry
Occurs when the spacecraft distance to a child body center reaches that child body's SOI radius from the outside:

```cpp
|r_sc_abs(t) - r_child_abs(t)| - child.soi_radius = 0
```

where both spacecraft and child states are expressed in the same inertial frame of the current central body's parent.

## 6.2 Search interval

The detector searches for the earliest event with:

- `t >= start_time`
- `t <= time_limit + time_epsilon`

If no event is found in that interval:
- return `Ok`
- emit `EventType::TimeLimit` or `EventType::None`

Pick one and lock it.

Recommended:
- `find_next_event` returns `Ok`
- `out_event.type = EventType::TimeLimit`
- `out_event.time = time_limit`

That makes the detector usable by later preview builders without special casing “no event found”.

## 6.3 Start-on-boundary behavior

A major ambiguity is what to do if the initial state is already exactly on a boundary at `start_time`.

Recommended rules:

- **Impact at start**: return `Impact` at `start_time`
- **SOI exit at start**: only return it if the radial trend is outward or if the state is outside by more than `root_epsilon`
- **Child SOI entry at start**: only return it if the trend is inward or if the state is inside by more than `root_epsilon`

This avoids spurious immediate events due to tiny roundoff, while still treating genuine start-state collisions as real.

Alternative rule:
- reject all start-on-boundary states as invalid

Do **not** choose that. The library will need to handle boundary states robustly.

## 6.4 Deterministic event ordering

If multiple events occur within `time_epsilon`, priority is:

1. impact
2. SOI entry into child
3. SOI exit
4. time limit

Among multiple sibling SOI entries within `time_epsilon`, choose the **lowest `BodyId`**.

This is part of the determinism contract and must be enforced centrally in exactly one place.

## 6.5 Bounded refinement

All root refinement must use:

- explicit tolerance: `root_epsilon`
- explicit max iteration count: `max_event_refine_iterations`
- explicit failure return on iteration-cap hit

The detector must never return a “best effort” partial refined root after iteration cap exhaustion.

## 6.6 Grazing encounters

The detector must not rely on sign changes alone.

For child SOI entry and possibly for impact/exit in edge cases, the distance-to-boundary function may touch zero tangentially without changing sign.

The detector shall therefore:

- step coarsely forward
- evaluate both sign changes and local minima
- refine candidate windows where:
  - a sign change occurs, or
  - a local minimum lies within `root_epsilon` of zero

## 6.7 Adaptive coarse stepping

Coarse scan step must be bounded so small child SOIs cannot be skipped.

The scan step should be bounded by a function like:

```cpp
dt_scan <= min_child_soi_radius / max(|v_rel|, v_floor)
```

with a lower bound:

```cpp
dt_scan >= time_epsilon * 100
```

The exact formula can vary, but the behavior must satisfy:

- steps shrink when relative speed is high
- steps shrink when the smallest child SOI is small
- steps do not collapse to pathological tiny values in ordinary cases

If anomaly-based stepping is available from Phase 2, prefer it near periapsis for high-eccentricity cases.

---

## 7. Test-driven development plan

Write the headers first, then the tests below before implementation.

## 7.1 Compile-only API tests

These ensure the API is usable and stable.

1. **`EventDetector_CanBeConstructedFromBodySystem`**
2. **`FindNextEvent_AcceptsSearchRequestAndOutputEvent`**
3. **`PredictedEvent_IsTriviallyCopyable`**
4. **`PredictedEvent_IsStandardLayout`**
5. **`EventSearchRequest_IsTriviallyCopyable`**
6. **`EventSearchRequest_IsStandardLayout`**

These should compile before any implementation exists.

---

## 7.2 Root-function unit tests

These should be written before the detector implementation.

### Impact root tests

7. **`ImpactFunction_IsNegativeInsideRadius`**
8. **`ImpactFunction_IsZeroOnRadius`**
9. **`ImpactFunction_IsPositiveOutsideRadius`**

### SOI exit root tests

10. **`SoiExitFunction_IsNegativeInsideSoi`**
11. **`SoiExitFunction_IsZeroOnSoiBoundary`**
12. **`SoiExitFunction_IsPositiveOutsideSoi`**

### Child entry root tests

13. **`ChildEntryFunction_IsPositiveOutsideChildSoi`**
14. **`ChildEntryFunction_IsZeroOnChildSoiBoundary`**
15. **`ChildEntryFunction_IsNegativeInsideChildSoi`**

### Numeric guard tests

16. **`RootFunctions_DoNotEmitNaN_ForFiniteInputs`**
17. **`RootFunctions_HandleVerySmallPositiveDistances`**
18. **`RootFunctions_HandleBoundaryValuesWithinRootEpsilon`**

These tests are small and local. They pin down sign conventions before any search logic exists.

---

## 7.3 Root refinement tests

The refinement logic should be tested independently of the full detector.

Use synthetic scalar functions where possible, not only orbital scenarios.

### Sign-change refinement tests

19. **`RefineRoot_FindsLinearSignChangeRoot`**
20. **`RefineRoot_FindsQuadraticRootWithinTolerance`**
21. **`RefineRoot_StopsWithinRootEpsilon`**
22. **`RefineRoot_ReturnsNoConvergenceOnIterationCap`**

### Grazing/minimum refinement tests

23. **`RefineMinimum_DetectsTangentialContactForParabola`**
   - Example scalar function: `f(t) = (t - 3)^2`

24. **`RefineMinimum_RejectsNearMissOutsideRootEpsilon`**
25. **`RefineMinimum_DoesNotEmitNaNOnFlatFunction`**

These tests matter because grazing handling is one of the core requirements of the spec.

---

## 7.4 Coarse scan / candidate generation tests

These verify that the scan brackets events correctly before refinement.

26. **`CoarseScan_BracketsImpactEventOnSimpleInboundTrajectory`**
27. **`CoarseScan_BracketsSoiExitOnSimpleOutboundTrajectory`**
28. **`CoarseScan_BracketsChildEntryOnCrossingTrajectory`**
29. **`CoarseScan_FindsGrazingCandidateWithoutSignChange`**
30. **`CoarseScan_UsesAdaptiveStepSmallEnoughToNotSkipSmallChildSoi`**
31. **`CoarseScan_RespectsMinimumStepFloor`**

Where possible, make these tests white-box at the helper level rather than full end-to-end detector tests. That keeps failures localized.

---

## 7.5 End-to-end detector tests: impact

Use simple single-body systems where no children exist.

32. **`FindNextEvent_ReturnsImpactForInboundCollision`**
33. **`FindNextEvent_ReturnsTimeLimitWhenNoImpactOccursBeforeLimit`**
34. **`FindNextEvent_ImpactAtStartTimeIsDetected`**
35. **`FindNextEvent_PrefersEarlierImpactOverLaterSoiExit`**

Use scenarios with analytically obvious behavior.

---

## 7.6 End-to-end detector tests: SOI exit

Use a body with a defined SOI and no children.

36. **`FindNextEvent_ReturnsSoiExitForOutboundTrajectory`**
37. **`FindNextEvent_DoesNotReturnSoiExitForBoundOrbitInsideLimit`**
38. **`FindNextEvent_SoiExitAtStartHandledPerBoundaryPolicy`**
39. **`FindNextEvent_ReturnsEarlierOfImpactAndSoiExit`**

---

## 7.7 End-to-end detector tests: child SOI entry

Use a parent body with at least one child.

40. **`FindNextEvent_ReturnsChildSoiEntryWhenTrajectoryIntersectsChildSoi`**
41. **`FindNextEvent_ReturnsEarliestChildEntryAmongMultipleChildren`**
42. **`FindNextEvent_UsesLowestBodyIdWhenSiblingEntriesTieWithinTimeEpsilon`**
43. **`FindNextEvent_DoesNotReportFalseChildEntryForNearMiss`**
44. **`FindNextEvent_DetectsGrazingChildEntry`**

These are core determinism tests.

---

## 7.8 Event ordering tests

These must exist before implementation because they define game-visible semantics.

45. **`EventOrdering_ImpactBeatsChildEntryWithinTimeEpsilon`**
46. **`EventOrdering_ChildEntryBeatsSoiExitWithinTimeEpsilon`**
47. **`EventOrdering_SoiExitBeatsTimeLimitWithinTimeEpsilon`**
48. **`EventOrdering_MultipleSiblingEntriesUseLowestBodyId`**
49. **`EventOrdering_StableAcrossRepeatedRuns`**

---

## 7.9 Start-state and boundary tests

50. **`Detector_StartExactlyInsideImpactRadiusReturnsImpactImmediately`**
51. **`Detector_StartExactlyOnChildBoundaryWithOutwardMotionDoesNotSpuriouslyEnter`**
52. **`Detector_StartExactlyOnChildBoundaryWithInwardMotionReturnsEntry`**
53. **`Detector_StartExactlyOnSoiBoundaryWithOutwardMotionReturnsExit`**
54. **`Detector_StartExactlyOnSoiBoundaryWithInwardMotionDoesNotSpuriouslyExit`**

These tests are where many event kernels fail.

---

## 7.10 Invalid-input and failure-mode tests

55. **`FindNextEvent_RejectsUnknownCentralBody`**
56. **`FindNextEvent_RejectsNonFiniteStartTime`**
57. **`FindNextEvent_RejectsNonFiniteTimeLimit`**
58. **`FindNextEvent_RejectsTimeLimitBeforeStartTime`**
59. **`FindNextEvent_RejectsNonFiniteInitialState`**
60. **`FindNextEvent_NoNaNsOnOrdinaryInvalidInputs`**
61. **`FindNextEvent_PropagatorFailureIsReportedExplicitly`**
62. **`FindNextEvent_RefinementFailureIsReportedExplicitly`**

---

## 7.11 Property tests

63. **`Detector_ReturnedEventTimeIsNeverEarlierThanStartTimeMinusEpsilon`**
64. **`Detector_ReturnedEventTimeIsNeverLaterThanTimeLimitPlusEpsilonUnlessFailure`**
65. **`Detector_ReturnedStateMatchesDirectPropagationToEventTime`**
66. **`Detector_EarliestEventIsInvariantAcrossRepeatedRunsOnSamePlatform`**
67. **`Detector_GrazingDetectionNeverReportsMissAsImpactWhenMinimumExceedsRootEpsilon`**

---

## 7.12 Fuzz tests

68. **`RandomValidSystemsAndStates_NoCrashesNoNaNs`**
69. **`RandomHighEccentricityCases_NoSkippedEventsWithinLimit`**
70. **`RandomSiblingChildSystems_EventOrderingAlwaysReturnsValidBodyIds`**
71. **`RandomBoundaryStarts_DoNotInfiniteLoop`**

---

## 8. Implementation plan after tests exist

Implement this phase in small slices.

## 8.1 Slice 1: headers only

Write or update:

- `event_detector.h`
- maybe `event_types.h`
- internal helper declarations for root functions and refinement

Add the compile-only tests first.

## 8.2 Slice 2: scalar root functions

Implement:

- impact distance-to-boundary
- SOI exit distance-to-boundary
- child entry distance-to-boundary

Pass tests 7–18.

## 8.3 Slice 3: bounded scalar refinement

Implement reusable bounded refinement for:

- sign-change roots
- local minima / tangent detection

Pass tests 19–25.

## 8.4 Slice 4: coarse scan helpers

Implement scan logic that:

- samples forward from `start_time`
- produces candidate brackets
- respects adaptive step sizing
- tracks minima for grazing detection

Pass tests 26–31.

## 8.5 Slice 5: end-to-end impact and SOI exit detection

Integrate propagation + scan + refinement for:
- impact
- SOI exit

Pass tests 32–39.

## 8.6 Slice 6: child SOI entry detection

Add child-iteration logic and parent-frame child ephemeris queries.

Pass tests 40–44.

## 8.7 Slice 7: deterministic event arbitration

Implement one centralized arbitration function that:

- merges near-simultaneous events within `time_epsilon`
- applies type priority
- applies lowest-`BodyId` tie-break among siblings

Pass tests 45–49.

## 8.8 Slice 8: boundary handling and failure hardening

Pass tests 50–71.

---

## 9. Internal implementation recommendations

These are implementation recommendations, not public API requirements.

## 9.1 Keep helper layers separate

Do not merge all of this into one giant function.

Separate:

- root evaluation
- scalar refinement
- coarse scan
- candidate generation
- event arbitration

That will make the TDD process much cleaner.

## 9.2 Use one canonical event comparison function

Recommended internal helper:

```cpp
bool event_precedes(const CandidateEvent& a,
                    const CandidateEvent& b,
                    const Tolerances& tol);
```

This function should encapsulate:

- earlier time wins if separated by more than `time_epsilon`
- otherwise priority order wins
- otherwise lower `BodyId` wins where applicable

All event-order decisions should go through this function.

## 9.3 Use parent-frame accumulation consistently for child entry

For child SOI entry in a current central-body frame:

- propagate spacecraft state in central-body frame to time `t`
- convert spacecraft state to the parent inertial frame if needed
- query child body state in the same frame
- subtract and evaluate distance

Do not mix coordinate frames inside root evaluation.

## 9.4 Clamp and sanitize all fragile scalar operations

Required:
- clamp inverse trig inputs to `[-1, 1]`
- treat tiny negative discriminants within epsilon as zero
- reject non-finite states up front
- never propagate NaNs into event comparison logic

## 9.5 Prevent infinite loops

The event detector itself should not loop indefinitely if:
- propagation repeatedly fails
- scan step collapses
- candidate refinement oscillates

Every loop must have an explicit bound.

---

## 10. Edge cases to decide before coding

These should not be left to implementation taste.

## 10.1 Time-limit semantics

Recommended:
- no event before horizon returns `TimeLimit` at exactly `time_limit`

That is more useful than `None` for later segment builders.

## 10.2 Starting outside current SOI

This can happen from bad input or partially reconstructed state.

Recommended:
- if outside by more than `root_epsilon`, return `SoiExit` at `start_time`

This makes the detector robust and simplifies later recovery code.

## 10.3 Starting inside a child SOI

Recommended:
- if already inside by more than `root_epsilon`, return `SoiEntry` at `start_time` for the lowest-`BodyId` qualifying child

This is a recovery-friendly behavior and fits the patching model.

## 10.4 Simultaneous entry into overlapping authored child SOIs

The authored data may be bad even if the tree is valid.

Recommended:
- if multiple children qualify within `time_epsilon`, choose lowest `BodyId`
- document that overlapping sibling SOIs are legal but ambiguous, resolved deterministically

## 10.5 Very short time limits

If `time_limit - start_time <= time_epsilon`, do not force a scan.
Check start-state boundary events, otherwise return `TimeLimit`.

---

## 11. Failure handling requirements

Every public event-detection API must:

- return an explicit `SolveStatus`
- leave outputs fully initialized
- never emit NaNs for ordinary invalid inputs
- never return partially refined event times as success if convergence failed

Required failure categories:
- invalid input
- propagation failure
- refinement no-convergence
- numerical failure

If there is no event before the time limit, that is **not** a failure.

---

## 12. Determinism requirements specific to Phase 3

Phase 3 is one of the most determinism-sensitive parts of the library.

The implementation shall avoid:

- unordered child traversal
- unstable sorting of candidate events
- floating-point comparisons done inconsistently in different code paths
- “best effort” convergence exits
- hidden dependence on insertion order of bodies

Required same-platform properties:

- same `BodySystem`, same request, same tolerances, same binary -> bit-identical event result
- tie cases resolve identically across repeated runs
- no platform-local nondeterminism from iteration order or uninitialized data

Recommended determinism tests:

- repeat the same event search many times in one process and compare raw bytes
- build the same body system from different insertion orders and verify event results are identical

---

## 13. Suggested file layout for this phase

```text
include/patched_conics/event_detector.h

src/event_detector.cpp
src/event_root_functions.cpp
src/event_refinement.cpp

tests/test_event_api.cpp
tests/test_event_root_functions.cpp
tests/test_event_refinement.cpp
tests/test_event_coarse_scan.cpp
tests/test_event_detector_impact.cpp
tests/test_event_detector_soi_exit.cpp
tests/test_event_detector_child_entry.cpp
tests/test_event_ordering.cpp
tests/test_event_boundaries.cpp
tests/test_event_invalid_input.cpp
tests/test_event_properties.cpp
tests/test_event_fuzz.cpp
```

If the codebase is still small, `src/event_detector.cpp` and one internal header are enough. The important part is conceptual separation, not file count.

---

## 14. Definition of done for Phase 3

Phase 3 is complete when:

- the event-detector API is header-defined and stable
- tests exist before implementation for all major event behaviors
- impact, SOI exit, and child SOI entry are detected correctly
- grazing encounters are detected without relying on sign changes alone
- refinement is bounded and failure is explicit
- deterministic event ordering and sibling tie-breaks are enforced
- no NaNs or crashes occur on ordinary invalid inputs
- same-platform repeated runs produce bit-identical results

---

## 15. Recommended execution order summary

1. draft `event_detector.h`
2. add compile-only API tests
3. write root-function unit tests
4. write scalar refinement tests
5. write coarse-scan tests
6. write end-to-end impact tests
7. write end-to-end SOI exit tests
8. write end-to-end child-entry tests
9. write event-ordering tests
10. write boundary and invalid-input tests
11. implement root functions
12. implement bounded refinement
13. implement coarse scan
14. integrate impact and SOI exit
15. integrate child-entry search
16. implement deterministic arbitration
17. run property and fuzz tests
18. harden determinism and edge-case behavior

This phase should end with a standalone event kernel that later phases can call without redesign.
