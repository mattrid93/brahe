# Patched Conics Library — Phase 1 Detailed Plan

## Scope

Build the minimum immutable system needed to answer:
- what bodies exist
- how they are related
- where a body is at time `t`
- what a body’s velocity is at time `t`
- what a body’s state is in an ancestor or root frame
- whether the body registry is valid and deterministic

This phase should produce:
- `BodyId`, `InvalidBody`
- `Vec2`, `State2`
- enums and `Tolerances`
- `BodyDef`
- `BodySystemBuilder`
- immutable `BodySystem`
- deterministic child lists ordered by ascending `BodyId`
- cached `parent`, `depth`, `root`, `children` ranges
- `position_in_parent`
- `velocity_in_parent`
- `state_in_ancestor_frame`
- `state_in_root_frame`

## Design decisions to lock now

- **Canonical body order:** sort internal storage by ascending `BodyId` during build. All child lists are also ascending `BodyId`. This avoids insertion-order dependence.
- **Immutable runtime object:** `BodySystemBuilder` is mutable; `BodySystem` is frozen after `build`.
- **No heap in query path:** allocation is allowed in builder/build; state queries must allocate nothing.
- **Ancestor queries only:** reject invalid ancestor requests explicitly.
- **Numeric behavior:** use `double` only; no SIMD or approximate trig.

## TDD plan

Write these tests first. Do not write implementation until the test list exists.

### A. Core layout and type tests

1. **`BodyId_InvalidBody_IsAllOnes`**
   - Assert `InvalidBody == ~BodyId{0}`.

2. **`Vec2_State2_BodyDef_AreTriviallyCopyable`**
   - `static_assert(std::is_trivially_copyable_v<...>)`.

3. **`Vec2_State2_BodyDef_AreStandardLayout`**
   - `static_assert(std::is_standard_layout_v<...>)`.

4. **`EnumUnderlyingTypes_AreStable`**
   - Assert exact underlying type choice if you fix one.
   - If you do not fix it, omit this test.

5. **`DefaultTolerances_MatchSpec`**
   - Assert the exact numeric defaults from §13.

### B. Builder and hierarchy validation tests

6. **`Build_Fails_OnDuplicateBodyId`**
   - Two bodies with same `BodyId`.
   - Expect `InvalidInput`.

7. **`Build_Fails_WhenParentMissing`**
   - Non-root body references nonexistent parent.
   - Expect `InvalidInput`.

8. **`Build_Fails_WhenNoRootExists`**
   - All bodies have parents.
   - Expect `InvalidInput`.

9. **`Build_Fails_WhenMultipleRootsExist`**
   - Two bodies with `parent_id == InvalidBody`.
   - Expect `InvalidInput`.

10. **`Build_Fails_OnCycle`**
    - A→B→C→A.
    - Expect `InvalidInput`.

11. **`Build_Fails_WhenNonRootOrbitRadiusNotPositive`**
    - `orbit_radius <= 0`.
    - Expect `InvalidInput`.

12. **`Build_Fails_WhenSoiNotLargerThanRadius`**
    - `soi_radius <= radius`.
    - Expect `InvalidInput`.

13. **`Build_Fails_WhenChildSoiDoesNotFitOrbit`**
    - `soi_radius >= orbit_radius`.
    - Expect `InvalidInput`.

14. **`Build_Succeeds_OnValidStarPlanetMoonTree`**
    - Simple three-body system.
    - Expect `Ok`.

15. **`GetBody_ReturnsNullForUnknownId`**
    - Query missing ID.
    - Expect `nullptr`.

16. **`GetBody_ReturnsStableDefinitionForKnownId`**
    - Built system returns exact authored values.

### C. Deterministic ordering and cached topology tests

17. **`Build_CanonicalizesStorageIndependentOfInsertionOrder`**
    - Add same bodies in two different orders.
    - After build, queries and traversal-visible behavior must match bit-for-bit.

18. **`SiblingTraversalOrder_IsAscendingBodyId`**
    - Same parent, unsorted insertion order.
    - Check child iteration order indirectly through a test-only accessor or query behavior.

19. **`DepthCache_IsCorrectForNestedTree`**
    - root depth 0, child depth 1, grandchild depth 2.

20. **`AncestorQuery_RejectsNonAncestor`**
    - Ask for state of moon in sibling frame.
    - Expect `InvalidInput` or equivalent explicit failure path.

### D. Circular ephemeris tests

Use a tiny deterministic system with known angles:
- `orbit_radius = 10`
- `angular_rate = 2`
- `phase_at_epoch = 0.25`

21. **`PositionInParent_AtEpoch_MatchesAnalyticFormula`**
   - Compare against `r*[cos(phase), sin(phase)]`.

22. **`VelocityInParent_AtEpoch_MatchesAnalyticFormula`**
   - Compare against `r*w*[-sin(phase), cos(phase)]`.

23. **`PositionInParent_AtPositiveTime_MatchesAnalyticFormula`**

24. **`VelocityInParent_AtPositiveTime_MatchesAnalyticFormula`**

25. **`RootBody_PositionAndVelocity_AreAlwaysZero`**
   - Any `t`.

26. **`NegativeTime_EphemerisStillMatchesAnalyticFormula`**
   - Since time is a signed `double`.

27. **`ZeroAngularRate_GivesConstantPositionAndZeroVelocity`**
   - Only if you allow `angular_rate == 0`.
   - If not allowed by policy, replace with validation failure test.

### E. Frame transform tests

Use:
- star at origin
- planet orbit radius 100, phase 0
- moon orbit radius 10, phase π/2
- choose `t = 0`

28. **`StateInAncestorFrame_PlanetToRoot_MatchesPositionInParent`**

29. **`StateInAncestorFrame_MoonToPlanet_MatchesDirectEphemeris`**

30. **`StateInAncestorFrame_MoonToRoot_EqualsPlanetPlusMoonState`**
   - Position and velocity both.

31. **`StateInRootFrame_EqualsStateInAncestorFrameWithRoot`**

32. **`Transforms_AreConsistentAcrossMultipleTimes`**
   - A small table of `t` values: `0`, `1`, `10`, `-3`.

33. **`AncestorAccumulation_DoesNotDependOnBuildOrder`**
   - Same system, different insertion orders, same answers bit-for-bit.

### F. Numeric guardrail tests for phase 1

34. **`Vec2_Normalize_RejectsNearZeroVector`**
   - Required for later burn-frame/event code.

35. **`Vec2_LengthSquared_AndDot_AreExactForSimpleInputs`**
   - Basic arithmetic sanity.

36. **`StateQueries_DoNotEmitNaN_OnValidInputs`**
   - Smoke test across a few bodies and times.

## Implementation sequence

Only after the tests exist.

### 1. Minimal public headers
Implement declarations only:
- ids, enums, POD structs
- `Tolerances`
- `BodySystemBuilder`
- `BodySystem`

Run compile-only tests.

### 2. `Vec2` primitives
Implement only what phase 1 tests need:
- add/subtract
- scalar multiply
- dot
- length / length_squared
- finite check
- safe normalize with explicit failure

Run the core numeric tests.

### 3. Builder storage and canonicalization
Implement builder as append-only.
In `build()`:
- copy defs into temporary array
- sort by `BodyId`
- detect duplicates
- map ids to indices with deterministic structure
- compute root count
- validate parent links
- compute depths
- detect cycles

Run hierarchy tests.

### 4. Immutable `BodySystem`
Store:
- sorted body defs
- parent index per body
- depth per body
- contiguous child ranges in ascending `BodyId`
- root index

Expose `get_body()`.

Run topology tests.

### 5. Ephemeris functions
Implement:
- root body special case
- `theta = phase + angular_rate * t`
- analytic `position_in_parent`
- analytic `velocity_in_parent`

Run ephemeris tests.

### 6. Ancestor/root transforms
Implement iterative accumulation from body up to ancestor:
- sum positions
- sum velocities
- reject invalid ancestor relation

Run transform tests.

### 7. Refactor for query-path allocation-free behavior
Check that query functions use no temporary heap objects.
Keep all scratch on stack.

### 8. Lock behavior with regression tests
Add one end-to-end star/planet/moon test that exercises:
- build
- lookup
- ephemeris
- root transform
- insertion-order independence

## Definition of done

Phase 1 is done when:
- all phase-1 tests pass
- `BodySystem::build/validate` rejects every invalid hierarchy in the spec
- body state queries match the analytic circular model
- ancestor/root transforms are correct
- query results are insertion-order independent
- all relevant types remain trivially serializable
- valid inputs produce no NaNs

## Expected code artifacts

- `include/patched_conics/types.h`
- `include/patched_conics/vec2.h`
- `include/patched_conics/body_system.h`
- `src/body_system.cpp`
- `tests/test_types.cpp`
- `tests/test_body_builder.cpp`
- `tests/test_ephemeris.cpp`
- `tests/test_frame_transforms.cpp`
