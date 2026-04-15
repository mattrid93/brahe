# Patched Conics Library — Phase 2 Detailed Spec

## Scope

Phase 2 implements the **two-body conic math layer** for 2D patched-conic simulation. This phase sits between the immutable world model from Phase 1 and the event/patch/trajectory logic planned for later phases.

The purpose of Phase 2 is to provide deterministic, bounded, analytically grounded operations on a spacecraft state relative to a single central body:

- classify the conic from Cartesian state
- compute core invariants
- convert between Cartesian state and 2D conic representation
- propagate analytically from one time to another
- expose anomaly conversion helpers
- handle ellipse, hyperbola, and parabolic-like edge cases without NaNs

This phase does **not** include:

- SOI event detection
- patch generation
- trajectory segment chaining
- burns
- Lambert
- multiplayer synchronization logic

Those depend on Phase 2 outputs, but are not implemented here.

---

## Goals

Phase 2 is complete when the library can, for any valid `(mu, r, v)` in a body-centered inertial frame:

- compute specific orbital energy
- compute eccentricity magnitude
- classify the orbit as ellipse / parabolic-like / hyperbola
- compute periapsis radius
- compute apoapsis radius where defined
- convert Cartesian state into a canonical 2D conic representation
- reconstruct Cartesian state from that representation
- propagate analytically by `dt >= 0`
- do all of the above deterministically on the same platform
- return explicit status on invalid input or numerical failure
- never emit NaNs on ordinary invalid input

---

## Recommended development order

For C++, the right order here is:

1. **Write the public headers first**
2. **Write the tests against those headers**
3. **Implement the smallest amount of code needed to pass the tests**
4. **Refine internals without changing the public contract unless tests prove the API is wrong**

That is the preferred workflow for this phase.

### Why headers first is a good idea here

For this phase, the API surface is a large part of the design. The main risk is not syntax; it is locking the wrong conceptual model. A header-first pass helps answer these questions before implementation details take over:

- what state is public and what is derived
- what exactly is considered invalid input
- which helpers are public versus internal
- what outputs are required for downstream phases
- how failure is represented
- whether later code can stay allocation-free and deterministic

The headers for Phase 2 should be treated as a **design draft**, not as production code. They exist to make the tests concrete.

### Recommended workflow rule

Do **not** write `.cpp` implementation until:

- the Phase 2 public headers compile
- the test list exists
- the expected status behavior is written down
- the canonical formulas and numerical conventions are written down

---

## Phase 2 deliverables

Phase 2 should produce:

- public header for conic math API
- internal header(s) for anomaly solving and numerical helpers
- unit tests for invariants, conversion, classification, propagation, and failure handling
- optional property tests for round-trip and time-composition behavior
- implementation of deterministic two-body analytic propagation in 2D

Expected file set:

- `include/patched_conics/conics.h`
- `include/patched_conics/two_body.h`
- `src/two_body.cpp`
- `src/conic_elements.cpp` or merged into `two_body.cpp`
- `src/anomaly_solvers.cpp` or internal-only helper section
- `tests/test_conic_classification.cpp`
- `tests/test_conic_invariants.cpp`
- `tests/test_conic_roundtrip.cpp`
- `tests/test_two_body_propagation.cpp`
- `tests/test_anomaly_conversions.cpp`
- `tests/test_invalid_input.cpp`

---

## Public API design draft

The Phase 1 spec already defines a minimal public surface:

```cpp
struct ConicState {
    BodyId central_body;
    double epoch;
    State2 state_at_epoch;
    ConicType type;
};

class TwoBody {
public:
    static SolveStatus propagate(double mu, const State2& initial, double dt, State2& out_state);
    static SolveStatus classify(double mu, const State2& state, ConicType& out_type);
    static double specific_energy(double mu, const State2& state);
    static double eccentricity(double mu, const State2& state);
    static double periapsis_radius(double mu, const State2& state);
    static double apoapsis_radius(double mu, const State2& state);  // infinity for open orbits
};
```

That is enough for the minimal spec, but Phase 2 needs more explicit structure to support robust TDD.

### Recommended header additions

Add a public 2D conic representation and explicit conversion functions.

```cpp
struct ConicElements2D {
    ConicType type;

    double mu;
    double semi_major_axis;      // negative for hyperbola, undefined-sign handling for parabola-like
    double eccentricity;
    double semi_latus_rectum;
    double argument_of_periapsis;
    double true_anomaly;

    double mean_anomaly;         // elliptic if defined
    double eccentric_anomaly;    // elliptic if defined
    double hyperbolic_anomaly;   // hyperbolic if defined
};

class TwoBody {
public:
    static SolveStatus classify(double mu, const State2& state, ConicType& out_type);
    static SolveStatus to_elements(double mu, const State2& state, ConicElements2D& out);
    static SolveStatus from_elements(const ConicElements2D& elements, State2& out_state);
    static SolveStatus propagate(double mu, const State2& initial, double dt, State2& out_state);

    static double specific_energy(double mu, const State2& state);
    static double specific_angular_momentum_z(const State2& state);
    static double eccentricity(double mu, const State2& state);
    static double semi_major_axis(double mu, const State2& state);
    static double semi_latus_rectum(double mu, const State2& state);
    static double periapsis_radius(double mu, const State2& state);
    static double apoapsis_radius(double mu, const State2& state);
};
```

### Why expose `to_elements` / `from_elements`

Because Phase 2 must test more than just propagation. It must lock:

- canonical meaning of the orbit representation
- round-trip correctness
- anomaly definitions
- behavior on near-degenerate states

Without explicit conversion functions, those tests become indirect and fragile.

### What should stay internal

These helpers are useful but should stay internal unless a downstream requirement appears:

- Kepler equation solvers
- universal-variable Stumpff helpers
- branch-selection helpers for near-parabolic propagation
- anomaly normalization helpers
- clamped inverse-trig helpers
- discriminant cleanup helpers

---

## Canonical mathematical model for Phase 2

All math in this phase is body-centered and 2D.

Given:

- `mu > 0`
- `r = (x, y)`
- `v = (vx, vy)`
- `|r| > 0`

The motion is governed by:

```text
r¨ = -mu * r / |r|^3
```

### Required invariants

Implement and test these first:

- radius: `r = |position|`
- speed squared: `v² = dot(v, v)`
- specific angular momentum scalar:
  - `h = x * vy - y * vx`
- specific orbital energy:
  - `eps = v² / 2 - mu / r`
- eccentricity vector:
  - `e_vec = ((v² - mu/r) * r - dot(r, v) * v * ???)`

Use the standard 2D equivalent of the 3D formula, implemented carefully. A safer implementation is:

```text
e_vec = (1/mu) * [ (v² - mu/r) * r - dot(r, v) * v ]
```

where `r` and `v` are vectors.

Then:

- `e = |e_vec|`
- semi-latus rectum:
  - `p = h² / mu`
- semi-major axis:
  - `a = -mu / (2 * eps)` when `eps` is not zero within tolerance

### Conic classification rules

Use the spec tolerance:

```cpp
parabolic_eccentricity_band = 1e-6
```

Required classification:

- `e < 1 - band` → `Ellipse`
- `|e - 1| <= band` → `ParabolaLike`
- `e > 1 + band` → `Hyperbola`

Do **not** classify based solely on specific energy. Use eccentricity-based classification because the spec makes the parabolic band part of the public determinism contract.

### Required geometric outputs

- periapsis radius:
  - `rp = p / (1 + e)`
- apoapsis radius:
  - ellipse: `ra = p / (1 - e)`
  - parabola-like / hyperbola: `+infinity`

---

## Propagation strategy

Phase 2 needs a deterministic propagation algorithm that is stable for:

- ellipses
- hyperbolas
- near-parabolic states

### Recommended first implementation

Use a **universal-variable propagator** with deterministic bounded iteration.

Reason:

- one framework can handle ellipse, hyperbola, and near-parabolic cases
- it avoids excessive branch fragmentation
- it composes well with TDD because the same input/output shape applies across conic types

### Acceptable alternative

A branch-based implementation is also acceptable:

- ellipse via eccentric anomaly and Kepler’s equation
- hyperbola via hyperbolic anomaly and hyperbolic Kepler equation
- parabola-like via special handling or universal-variable fallback

If branch-based propagation is used, the near-parabolic branch must still be explicitly tested and must not devolve into unstable `1/(1-e)` behavior.

### Propagation contract

`TwoBody::propagate(mu, initial, dt, out)` shall:

- require `mu > 0`
- require finite input state
- require `|r| > position_epsilon`
- require `dt >= 0`
- return failure on invalid input
- preserve orbital invariants to within test tolerance for valid inputs
- return `out_state` fully initialized on success
- never leave `out_state` partially uninitialized on failure

---

## Numerical policy for Phase 2

Phase 2 must implement the numerical rules from the main spec and make them testable.

### Required rules

- clamp `acos` / `asin` inputs to `[-1, 1]`
- treat tiny negative discriminants within epsilon as zero
- reject normalization of near-zero vectors
- never return NaN on ordinary invalid input
- bounded iteration counts are failure conditions, not partial-success exits
- no tolerance-dependent hidden behavior outside `Tolerances`

### Additional Phase 2 guidance

- use `std::isfinite` checks on all externally supplied doubles
- normalize angles only where the API needs canonical values
- do not normalize away physically meaningful sign information for hyperbolic anomalies
- explicitly document whether `true_anomaly` is returned in `(-pi, pi]` or `[0, 2pi)`

Recommended canonical choice:

- `argument_of_periapsis`: `(-pi, pi]`
- `true_anomaly`: `(-pi, pi]`

Lock this in tests.

---

## TDD strategy for Phase 2

The test plan should be written before implementation. The most effective structure is:

1. **header/API compile tests**
2. **invariant tests**
3. **classification tests**
4. **conversion round-trip tests**
5. **anomaly conversion tests**
6. **propagation tests**
7. **invalid-input and numerical-guardrail tests**
8. **property tests**

That order is intentional. It lets the core formulas settle before propagation complexity appears.

---

## Tests to write before implementation

## A. Header and API compile tests

These are first. Their job is to force the API shape to exist.

1. **`TwoBody_HeaderCompiles`**
   - Include public headers only.
   - Instantiate `State2`, `ConicElements2D`, and call signatures in unevaluated context.

2. **`ConicElements2D_IsTriviallyCopyable`**
   - `static_assert`.

3. **`ConicElements2D_IsStandardLayout`**
   - `static_assert`.

4. **`TwoBody_PublicFunctions_AreStaticAndCallable`**
   - Build-only test.

These tests may seem shallow, but they are useful in C++ because they lock the API before implementation drift starts.

---

## B. Invariant tests

These should be the first runtime tests.

Use hand-authored states with easy expected values.

### Test fixture suggestions

Use `mu = 1` for simple algebraic cases.

#### Circular orbit fixture

```text
r = (1, 0)
v = (0, 1)
```

Expected:

- `r = 1`
- `v² = 1`
- `h = 1`
- `eps = -0.5`
- `e = 0`
- `a = 1`
- `p = 1`
- `rp = 1`
- `ra = 1`

### Tests

5. **`SpecificEnergy_CircularUnitOrbit_IsMinusHalf`**

6. **`AngularMomentumZ_CircularUnitOrbit_IsOne`**

7. **`Eccentricity_CircularUnitOrbit_IsZero`**

8. **`SemiMajorAxis_CircularUnitOrbit_IsOne`**

9. **`SemiLatusRectum_CircularUnitOrbit_IsOne`**

10. **`Periapsis_CircularUnitOrbit_IsOne`**

11. **`Apoapsis_CircularUnitOrbit_IsOne`**

### Elliptic non-circular fixture
Choose a known periapsis state:

```text
mu = 1
a = 2
e = 0.5
rp = a(1-e) = 1
vp = sqrt(mu * (1+e)/(a(1-e))) = sqrt(1.5)
r = (1, 0)
v = (0, sqrt(1.5))
```

12. **`Invariants_AtEllipticPeriapsis_MatchExpectedValues`**

### Hyperbolic fixture
Choose:

```text
mu = 1
r = (1, 0)
v = (0, 2)
```

Expected:

- `eps = 1`
- `h = 2`
- `p = 4`
- `e = 3`
- `rp = 1`
- `ra = infinity`

13. **`Invariants_HyperbolicCase_MatchExpectedValues`**

---

## C. Conic classification tests

These lock the public classification semantics.

14. **`Classify_CircularOrbit_AsEllipse`**

15. **`Classify_EllipticOrbit_AsEllipse`**

16. **`Classify_HyperbolicOrbit_AsHyperbola`**

17. **`Classify_EccentricityWithinBandOfOne_AsParabolaLike`**
   - Construct or mock a state numerically near `e = 1`.

18. **`Classify_EccentricityJustBelowBand_AsEllipse`**

19. **`Classify_EccentricityJustAboveBand_AsHyperbola`**

20. **`Classify_InvalidMu_Fails`**

21. **`Classify_ZeroRadius_Fails`**

22. **`Classify_NonFiniteInput_Fails`**

---

## D. Cartesian ↔ conic conversion tests

These are essential. Write them before propagation.

### Required semantics

`to_elements()` should output a canonical representation for valid states.

`from_elements()` should reconstruct a state consistent with those elements.

Because some angular quantities wrap, compare:

- radius and speed vectors directly after round-trip
- angle outputs using wrapped-angle comparison helper

### Tests

23. **`ToElements_CircularOrbit_ProducesExpectedCanonicalValues`**
   - Circular case is degenerate in periapsis direction.
   - The spec should define a canonical convention.

### Required circular convention

For `e <= position tolerance equivalent`, choose:

- `argument_of_periapsis = atan2(r.y, r.x)` at epoch state
- `true_anomaly = 0`

or

- `argument_of_periapsis = 0`
- `true_anomaly = atan2(r.y, r.x)`

Pick one and test it. Do not leave circular handling implicit.

Recommended choice:

- `argument_of_periapsis = 0`
- `true_anomaly = polar_angle(position)`

This is simpler and avoids rotating the whole frame to a fake periapsis.

24. **`ToElements_EllipticPeriapsis_ProducesExpectedValues`**

25. **`ToElements_HyperbolicPeriapsis_ProducesExpectedValues`**

26. **`FromElements_Ellipse_ReconstructsOriginalState`**

27. **`FromElements_Hyperbola_ReconstructsOriginalState`**

28. **`RoundTrip_CartesianToElementsToCartesian_Ellipse_MatchesOriginal`**

29. **`RoundTrip_CartesianToElementsToCartesian_Hyperbola_MatchesOriginal`**

30. **`RoundTrip_CartesianToElementsToCartesian_ParabolaLike_RemainsFiniteAndConsistent`**
   - Compare loose invariants, not exact anomaly values.

31. **`FromElements_InvalidInput_FailsWithoutNaN`**
   - Negative `mu`, negative `p`, nonsensical eccentricity combinations.

---

## E. Anomaly conversion tests

Even if anomaly helpers remain internal, write tests either directly against internal test hooks or indirectly through `to_elements` and `from_elements`.

### Ellipse tests

32. **`TrueToEccentricToTrue_RoundTrip_Ellipse`**
   - Multiple values, including near `0`, `pi/2`, `pi`, `-pi + eps`.

33. **`MeanToEccentricToMean_RoundTrip_Ellipse`**

34. **`EccentricToMean_MonotonicOverPrincipalRange`**

### Hyperbola tests

35. **`TrueToHyperbolicToTrue_RoundTrip_Hyperbola`**

36. **`MeanToHyperbolicToMean_RoundTrip_Hyperbola`**

37. **`HyperbolicConversions_RemainFiniteNearAsymptote`**

### Guardrail tests

38. **`InverseTrigInputs_AreClamped_NotNaN`**

39. **`TinyNegativeDiscriminant_IsTreatedAsZero`**

---

## F. Propagation tests

Write these only after the invariant and conversion tests are specified.

### F1. Basic propagation tests

40. **`Propagate_CircularOrbit_ByZeroDt_ReturnsSameState`**

41. **`Propagate_CircularOrbit_ByQuarterPeriod_RotatesStateCorrectly`**
   - Use the unit circular orbit.
   - `period = 2pi`
   - at `dt = pi/2`, expect `r ≈ (0,1)`, `v ≈ (-1,0)`.

42. **`Propagate_CircularOrbit_ByHalfPeriod_RotatesStateCorrectly`**

43. **`Propagate_EllipticOrbit_OneFullPeriod_ReturnsOriginalState`**
   - Compare within tolerance.

44. **`Propagate_EllipticOrbit_TwoStepMatchesSingleStep`**
   - propagate `0→t1`, then `t1→t2`
   - compare with direct `0→t2`

45. **`Propagate_HyperbolicOrbit_PreservesEnergyAndAngularMomentum`**

46. **`Propagate_ParabolaLikeCase_RemainsFinite`**

### F2. High-eccentricity and near-degenerate tests

47. **`Propagate_HighEccentricityEllipse_DoesNotExplodeNearPeriapsis`**
   - Example `e = 0.999`.

48. **`Propagate_NearParabolicBelowBand_RemainsStable`**

49. **`Propagate_NearParabolicAboveBand_RemainsStable`**

50. **`Propagate_VerySmallDt_ProducesConsistentState`**

51. **`Propagate_LargeDtEllipse_MultiplePeriodsStillConsistent`**

### F3. Directionality and physical sanity tests

52. **`Propagate_ForwardTime_DoesNotReverseAngularMomentumSign`**

53. **`Propagate_PericenterState_MovesAwayFromPericenterAfterPositiveDt`**
   - For a non-circular ellipse.

54. **`Propagate_HyperbolicOutboundState_RemainsOutboundForSufficientlyLargeDt`**

---

## G. Invalid-input and failure-behavior tests

These are mandatory.

55. **`Propagate_InvalidMu_Fails`**

56. **`Propagate_ZeroRadius_Fails`**

57. **`Propagate_NegativeDt_Fails`**

58. **`Propagate_NonFiniteState_Fails`**

59. **`Propagate_OutputStateIsFullyDefinedOnSuccess`**

60. **`Propagate_DoesNotWriteGarbageOnFailure`**
   - Pre-fill `out_state` with sentinel values and verify documented behavior.
   - Recommended behavior: leave unchanged on failure.

61. **`NoPublicFunctionReturnsNaNOnOrdinaryInvalidInput`**
   - Classification, conversion, propagation.

62. **`IterationCapReturnsNoConvergence_NotPartialSolution`**
   - If solver internals can be forced into cap through a test hook.

---

## H. Property tests

These are strongly recommended in this phase.

63. **`RoundTrip_ElementsAndCartesian_PreserveInvariants`**
   - Random valid ellipses/hyperbolas.

64. **`Propagation_ComposesInTime`**
   - `prop(dt1)` then `prop(dt2)` equals `prop(dt1 + dt2)` within tolerance.

65. **`Propagation_PreservesSpecificEnergy`**

66. **`Propagation_PreservesAngularMomentum`**

67. **`SameInputsSamePlatform_ProduceBitIdenticalOutputsAcrossRuns`**
   - Especially important for deterministic math code.

---

## Suggested implementation sequence after tests exist

## 1. Draft and freeze the headers

Write only declarations first:

- `ConicElements2D`
- `TwoBody`
- any minimal numeric helper declarations that tests require

Compile. Do not implement behavior yet.

### Header review checklist

Before writing tests, check:

- are the outputs sufficient for later event detection and Lambert work?
- are failure states explicit?
- are any public fields redundant or underdefined?
- is circular-orbit canonicalization specified?
- is near-parabolic behavior represented clearly enough?

If the answer is no, fix the headers first.

---

## 2. Implement invariant helpers first

Functions to implement first:

- `specific_energy`
- `specific_angular_momentum_z`
- `eccentricity`
- `semi_latus_rectum`
- `semi_major_axis`
- `periapsis_radius`
- `apoapsis_radius`
- `classify`

These should pass tests B and C before any conversion or propagation code exists.

---

## 3. Implement `to_elements`

This is the next milestone. It forces all orbital geometry conventions to become concrete.

### Required behavior

- compute `e_vec`
- compute `e`
- compute `p`
- compute `a` where meaningful
- compute `argument_of_periapsis`
- compute `true_anomaly`
- apply circular canonicalization rule
- classify the conic
- avoid NaNs in all angle computations

Pass conversion half-tests before moving on.

---

## 4. Implement `from_elements`

Construct Cartesian state from canonical conic data.

Recommended formula path:

- perifocal position:
  - `r_pf = [ p cos(nu)/(1+e cos(nu)), p sin(nu)/(1+e cos(nu)) ]`
- perifocal velocity:
  - `v_pf = sqrt(mu/p) * [ -sin(nu), e + cos(nu) ]`
- rotate by `argument_of_periapsis`

Pass round-trip tests.

---

## 5. Implement anomaly helpers and solvers

Only now add the machinery for:

- true ↔ eccentric anomaly
- true ↔ hyperbolic anomaly
- mean ↔ eccentric anomaly via bounded solve
- mean ↔ hyperbolic anomaly via bounded solve

These can remain internal. They should be testable either directly or via internal test-only exposure.

---

## 6. Implement propagation

At this point choose one of:

- universal-variable propagator
- branch-based conic-specific propagation

### Recommendation

Use universal variables if you are comfortable implementing and testing them. Otherwise:

- implement ellipse branch first
- implement hyperbola branch second
- implement parabola-like fallback third

The branch-based path is easier to debug in C++ if orbital mechanics is not already familiar.

### Practical recommendation for this project

Given the library’s constraints and the stated uncertainty about C++, the safest plan is:

- **keep the public API small**
- **write the headers first**
- **write the tests next**
- **implement branch-based propagation first**
- **switch to universal variables later only if the parabola-like branch becomes too awkward**

That minimizes implementation risk.

---

## 7. Add failure-path hardening

After propagation works on nominal cases, explicitly harden:

- invalid input paths
- iteration-cap behavior
- non-finite intermediate checks
- unchanged output-on-failure contract

Do not treat this as cleanup. It is part of the core Phase 2 deliverable.

---

## Canonical conventions to decide and lock in tests

These must be specified before implementation.

### 1. Circular orbit angle convention

Recommended:

- if `e <= circular_epsilon`, set `argument_of_periapsis = 0`
- set `true_anomaly = atan2(y, x)`

This makes circular states representable without inventing a noisy periapsis direction.

### 2. Signed angular momentum convention

Use:

```text
h = x * vy - y * vx
```

Positive means CCW. Negative means CW.

Phase 2 does not need to reject negative `h`; it should represent both directions.

### 3. Angle range convention

Recommended:

- all public angular outputs in `(-pi, pi]`

### 4. Failure-write convention

Recommended:

- output structs remain unchanged on failure

This is easy to test and safer for engine code.

### 5. Near-parabolic `a` convention

For `ParabolaLike`, `a` becomes numerically unstable. Choose one explicit rule:

- still compute `a = -mu/(2eps)` if finite
- or set `a = infinity`
- or leave it unspecified and document that callers must not use it

Recommended:

- preserve computed finite `a` if numerically valid
- document that callers must not branch on `a` for `ParabolaLike`

---

## Phase 2 acceptance criteria

Phase 2 is complete when:

- public headers exist and compile cleanly
- all Phase 2 unit tests described above exist
- invariants match known analytic expectations
- classification follows the eccentricity-band rule from the spec
- Cartesian ↔ conic conversion works for ellipse and hyperbola
- circular handling follows a documented canonical convention
- analytic propagation passes closed-orbit and composition tests
- parabolic-like cases remain finite and explicitly classified
- invalid input returns explicit failure and does not produce NaNs
- same-platform repeated runs are bit-stable

---

## Risks and likely trouble spots

These are the areas most likely to cause churn.

### 1. Circular and near-circular angle ambiguity

This is the biggest representational issue. Solve it in the spec and tests, not in code.

### 2. Near-parabolic instability

Do not let `1 - e` denominators creep into public behavior without tests around the band edges.

### 3. Hyperbolic anomaly sign errors

These are common. Write explicit signed round-trip tests.

### 4. Hidden partial-success behavior

Iteration caps must fail cleanly. Never return approximate partial states under a success status.

### 5. API drift during implementation

This is why the headers should exist first.

---

## Recommended next concrete action

Write these artifacts in this exact order:

1. `include/patched_conics/two_body.h`
2. `include/patched_conics/conics.h`
3. compile-only API tests
4. invariant tests
5. classification tests
6. conversion tests
7. propagation tests
8. implementation

That is the cleanest TDD path for this phase.
