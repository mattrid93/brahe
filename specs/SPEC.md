# Patched Conics Library Spec

**Target:** game engine integration
**Scope:** 2D only, fixed circular orbits, deterministic-enough, fast, stable, no n-body integration
**Primary use cases:** trajectory preview, maneuver planning, transfer generation, encounter prediction
**Multiplayer model:** server-authoritative simulation with client-side deterministic preview (hybrid)

---

## 1. Goals

The library shall provide:

- deterministic 2D patched-conic trajectory propagation (same-platform bit-stable, cross-platform best-effort)
- circular-orbit ephemerides for stars, planets, and moons
- sphere-of-influence patching between parent and child bodies
- impulsive maneuvers
- Lambert transfer solving between two positions in the same central-body frame
- encounter prediction and orbit-segment chaining
- stable runtime behavior suitable for per-frame or on-demand planning
- serializable state suitable for authoritative server snapshots

The library shall **not** provide:

- n-body integration
- perturbations
- non-circular body orbits
- 3D dynamics
- atmospheric drag
- low-thrust propagation
- relativistic or real-world time-scale handling
- bit-exact cross-architecture determinism (see §11.2)

---

## 2. Design constraints

- 2D only
- double precision throughout
- single continuous simulation time scalar (see §4.3 for range)
- all bodies on fixed circular orbits
- piecewise two-body propagation
- instantaneous SOI transitions
- instantaneous burns
- same-platform deterministic; cross-platform agreement is best-effort
- no heap allocation in hot-path propagation unless explicitly configured
- physics path must avoid unordered container traversal, parallelism, and pointer-address-dependent ordering

---

## 3. World model

### 3.1 Body hierarchy

Bodies form a tree.

Example:

```
star
├── planet A
│   ├── moon A1
│   └── moon A2
├── planet B
└── planet C
```

Each body orbits exactly one parent, except the root body.

### 3.2 Body assumptions

Each non-root body has:

- fixed circular orbit around parent
- fixed orbital radius
- fixed angular rate
- fixed phase at epoch 0
- fixed SOI radius
- fixed physical radius

No body-body interactions are simulated beyond patched-conic frame switching.

### 3.3 Immutability

After `BodySystem::validate()` succeeds, the body hierarchy shall be immutable. Runtime authoring is supported via a separate builder type that produces an immutable `BodySystem`. This permits caching of derived data (depth, ancestor chains) and is required for the determinism guarantees in §11.

---

## 4. Units and conventions

The library shall enforce one consistent unit system, selected by integrator at setup.

Recommended default:

- distance: meters or kilometers
- time: seconds
- angle: radians
- gravitational parameter: distance³ / time²
- velocity: distance / time

### 4.1 Coordinate conventions

2D vectors are Cartesian: `Vec2 = { x, y }`

Angle convention:

- 0 rad along +X
- positive rotation is counterclockwise (CCW)

Burn frame conventions:

- **prograde** = unit vector along velocity
- **retrograde** = -prograde
- **radial-out** = unit vector along position (from central body)
- **radial-in** = -radial-out
- **normal** = prograde rotated +90° (CCW)
- **anti-normal** = -normal

### 4.2 Reference frames

Supported frames:

- **Local inertial frame of body B**: origin at body center, non-rotating
- **Parent inertial frame**: origin at parent center, non-rotating

The root body frame is the global inertial frame. All trajectory segments are defined relative to a single central body.

### 4.3 Time

Simulation time is a single `double` scalar. Supported range is approximately ±10⁹ seconds (~30 years) from epoch 0 with sub-microsecond resolution. Sessions exceeding this range should periodically rebase the epoch (server-coordinated for multiplayer).

---

## 5. Core mathematical model

For each segment, the spacecraft follows a pure two-body conic around the current central body.

Supported conic types:

- ellipse
- near-parabolic (within tolerance band, see §13)
- hyperbola

Given state relative to central body `r(t), v(t)`, the segment is propagated analytically under:

```
r¨ = -mu * r / |r|^3
```

A trajectory is a chain of segments:

```
Segment[0], Segment[1], ..., Segment[n]
```

Each segment ends due to one of:

- SOI exit
- SOI entry into child body
- impact with central body
- user burn
- explicit end time / planning horizon

---

## 6. Circular-orbit ephemeris model

For body B orbiting parent P:

Inputs: `orbit_radius`, `angular_rate`, `phase_at_epoch`

At simulation time t:

```
theta(t)       = phase_at_epoch + angular_rate * t
pos_B_in_P(t)  = orbit_radius * [cos(theta), sin(theta)]
vel_B_in_P(t)  = orbit_radius * angular_rate * [-sin(theta), cos(theta)]
```

The root body position is always (0,0) and velocity (0,0) in global inertial frame.

The library shall support recursive frame accumulation so any body state can be queried in any ancestor frame.

---

## 7. Sphere of influence model

Each body has a fixed SOI radius. The library shall treat SOI as a hard boundary.

### 7.1 Transition rules

- if trajectory in parent-body frame intersects child SOI, patch into child-centered frame
- if trajectory in child-centered frame exceeds child SOI, patch back to parent-centered frame

### 7.2 Patch event state

- compute spacecraft absolute inertial state at event time
- compute target body absolute inertial state at event time
- subtract to obtain new relative state in new central-body frame

Position and velocity shall be continuous across patch boundaries in inertial space to within `position_epsilon` and `velocity_epsilon`. Acceleration is allowed to be discontinuous.

### 7.3 Re-entry suppression

After a patch, the new segment shall suppress detection of the SOI boundary it just crossed for a duration of `time_epsilon` OR until the spacecraft moves at least `position_epsilon * 10` away from that boundary, whichever occurs first. This prevents immediate re-patch loops caused by floating-point roundoff in the patch subtraction.

---

## 8. Lambert scope

The library shall include a 2D single-revolution Lambert solver.

### 8.1 Required use

- generate transfer between two points around the same central body
- support both rotation directions
- support elliptic and hyperbolic transfers where physically valid

The solver is **not** required to support:

- multi-revolution solutions
- 3D transfers
- universal optimization over multiple bodies

### 8.2 Lambert inputs

```cpp
struct LambertRequest {
    double mu;
    Vec2 r1;
    Vec2 r2;
    double dt;
    LambertDirection direction;  // CCW or CW
};
```

**Branch parameterization rationale:** in 2D, "short-way / long-way" and "prograde / retrograde" are not independent — once a rotation direction is chosen, the transfer angle is uniquely determined by the geometry. The API uses rotation direction (CCW/CW) as the single branch selector to avoid ambiguity, particularly when the transfer angle is near π.

### 8.3 Lambert outputs

```cpp
struct LambertSolution {
    SolveStatus status;
    Vec2 v1;
    Vec2 v2;
    ConicType conic_type;
    int iterations;       // diagnostic only
    double residual;      // diagnostic only
};
```

### 8.4 Lambert validity rules

The solver shall return failure for:

- `|r1| == 0` or `|r2| == 0`
- `dt <= 0`
- degenerate geometry outside numerical tolerance
- non-convergence
- physically invalid direction request

### 8.5 Lambert robustness requirements

The solver shall explicitly handle:

- nearly equal position vectors
- transfer angle near 0
- transfer angle near π (collinear `r1`, `r2`): in 2D the solution is uniquely determined by the rotation direction; the solver shall dispatch to the directional branch automatically rather than failing
- very short TOF
- very long TOF within single-rev regime
- hyperbolic high-energy solutions
- floating-point drift in inverse trig inputs

The solver shall never emit NaNs on ordinary invalid input. It must return a failure status instead.

### 8.6 Lambert convergence

Convergence shall be defined by a residual threshold (`lambert_residual_epsilon`). The iteration cap (`max_lambert_iterations`) is a failure signal, not a normal exit path. Hitting the cap returns `NoConvergence`, never a "best-effort" partial solution. This is required for determinism.

### 8.7 Lambert implementation recommendation

Recommended first implementation:

- a 2D-adapted universal-variable or Izzo-style single-rev solver
- direction-explicit API
- bounded iteration count
- residual-based convergence check
- fallback bracketing if Newton-style iteration fails

---

## 9. Functional requirements

### 9.1 Body registry

The library shall provide immutable body definitions (see §3.3).

```cpp
struct BodyDef {
    BodyId id;
    BodyId parent_id;          // InvalidBody for root
    double mu;
    double radius;
    double soi_radius;
    double orbit_radius;       // ignored for root
    double angular_rate;       // ignored for root
    double phase_at_epoch;     // ignored for root
};

constexpr BodyId InvalidBody = ~BodyId{0};
```

Validation requirements:

- tree validity (acyclic, single root)
- unique IDs
- parent must exist (or be `InvalidBody` for root)
- orbit radius > 0 for non-root
- soi radius > radius
- soi radius < orbit radius (child SOI must fit inside parent's domain near the child)

### 9.2 State queries

The library shall provide:

- body position in parent frame at time t
- body velocity in parent frame at time t
- body position/velocity in any ancestor frame at time t
- spacecraft state in current segment frame at time t

### 9.3 Conic conversion utilities

Required conversions:

- Cartesian state → orbital invariants
- Cartesian state → conic elements in 2D
- conic elements → Cartesian state
- true anomaly / eccentric anomaly / hyperbolic anomaly conversions
- periapsis, apoapsis, semi-major axis, eccentricity, mean motion where defined

`apoapsis_radius` returns `std::numeric_limits<double>::infinity()` for open orbits (parabolic or hyperbolic).

### 9.4 Propagation

Given `(r0, v0, t0, mu)`, the library shall propagate analytically to any `t >= t0`.

Outputs: `r(t)`, `v(t)`, optional anomaly / orbital phase.

The propagator shall support elliptical, hyperbolic, and near-parabolic (within tolerance) orbits.

### 9.5 Event detection

The propagator shall predict the earliest event after a given start time among:

- impact with central body radius
- SOI exit
- SOI entry into any child body of current central body
- user-defined time limit

#### 9.5.1 Event ordering

If multiple events occur within `time_epsilon`, the priority order is:

1. impact
2. SOI entry into child (tiebreaker among multiple sibling SOI entries: lowest `BodyId`)
3. SOI exit
4. time limit

The deterministic tiebreaker is required for multiplayer agreement.

#### 9.5.2 Event detection method

The detector performs a two-stage scan: a coarse forward step to bracket candidates, then a scalar root-solve to refine. The implementation shall handle:

- **Sign-change brackets**: standard root in `f(t) = |r_rel(t)| - boundary_radius`
- **Grazing encounters**: the distance function may be tangent to zero without sign change. The scan must also track local minima of the distance-to-boundary function and refine candidates where the minimum is within `root_epsilon` of zero.
- **Adaptive step size**: the coarse scan step shall be bounded by `min(child_soi_radius) / |v_rel|` across all candidate child bodies, so that no SOI smaller than the smallest child can be skipped. A floor of `time_epsilon * 100` prevents pathological tiny steps.
- **High-eccentricity orbits**: anomaly-based stepping is preferred over time-based stepping near periapsis to avoid skipping events.

Functions to solve:

- impact: `|r_rel(t)| - body.radius = 0`
- SOI exit: `|r_rel(t)| - body.soi_radius = 0`
- child SOI entry: `|r_sc_abs(t) - r_child_abs(t)| - child.soi_radius = 0`

Root refinement uses bounded iteration with explicit tolerance. Hitting the iteration cap returns a failure status, not a partial result.

### 9.6 Patch generation

#### Parent → child patch

At event time `t_e`:

1. compute spacecraft absolute state in parent inertial frame
2. compute child absolute state in parent inertial frame
3. subtract child state from spacecraft state
4. new segment central body = child
5. apply re-entry suppression (§7.3)

#### Child → parent patch

At event time `t_e`:

1. compute spacecraft relative state in child frame
2. compute child absolute state in parent frame
3. add child state to spacecraft state
4. new segment central body = parent
5. apply re-entry suppression (§7.3)

### 9.7 Maneuver support

```cpp
struct ImpulsiveBurn {
    double time;
    Vec2 delta_v;
    BurnFrame frame;
};

enum class BurnFrame {
    Inertial,            // delta_v in current segment's inertial frame
    ProgradeRadial,      // x = prograde, y = radial-out
    ProgradeNormal,      // x = prograde, y = normal (+90° CCW from prograde)
};
```

At burn time:

1. propagate current segment to burn time
2. transform delta_v into segment inertial frame per `BurnFrame`
3. add delta_v
4. terminate current segment
5. create new segment from updated state
6. recompute downstream patches

### 9.8 Trajectory recomputation

The API shall expose `rebuild_from(time T)` as a first-class operation. This is the hot path for both burn editing in planning UI and authoritative state reconciliation in multiplayer (§11.5).

### 9.9 Transfer planning

The library shall support Lambert-driven transfer generation between two bodies orbiting the same parent.

Input:

- departure body
- arrival body
- departure time
- arrival time
- Lambert direction

Steps:

1. query departure body absolute position at t0 in common parent frame
2. query arrival body absolute position at t1 in same frame
3. solve Lambert in parent frame
4. compute departure v_inf relative to departure body
5. compute arrival v_inf relative to arrival body
6. return transfer solution record

Conversion of v_inf to parking-orbit escape/capture burns is intentionally outside the core library.

---

## 10. API surface

### 10.1 Public types

```cpp
using BodyId = uint32_t;
constexpr BodyId InvalidBody = ~BodyId{0};

struct Vec2 { double x, y; };
struct State2 { Vec2 r, v; };

enum class ConicType { Ellipse, ParabolaLike, Hyperbola };

enum class EventType {
    None, Impact, SoiEntry, SoiExit, Burn, TimeLimit
};

enum class SolveStatus {
    Ok, InvalidInput, NoConvergence, NoSolution, NumericalFailure
};

enum class LambertDirection { CCW, CW };
enum class BurnFrame { Inertial, ProgradeRadial, ProgradeNormal };
```

### 10.2 Body system API

```cpp
class BodySystemBuilder {
public:
    void add_body(const BodyDef& def);
    SolveStatus build(BodySystem& out) const;
};

class BodySystem {
public:
    SolveStatus validate() const;
    const BodyDef* get_body(BodyId id) const;
    Vec2 position_in_parent(BodyId id, double t) const;
    Vec2 velocity_in_parent(BodyId id, double t) const;
    State2 state_in_ancestor_frame(BodyId body, BodyId ancestor, double t) const;
    State2 state_in_root_frame(BodyId body, double t) const;
};
```

### 10.3 Conic API

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

### 10.4 Lambert API

```cpp
class LambertSolver {
public:
    LambertSolution solve(const LambertRequest& req) const;
};
```

### 10.5 Event API

```cpp
struct TrajectoryEvent {
    EventType type;
    double time;
    BodyId from_body;
    BodyId to_body;     // InvalidBody if not applicable
    State2 state_before;
    State2 state_after; // for burns / patches
};
```

### 10.6 Segment and trajectory API

```cpp
struct Segment {
    BodyId central_body;
    double start_time;
    State2 initial_state;
    double end_time;
    EventType end_reason;
};

// Heap-allocated variant for planning
struct Trajectory {
    std::vector<Segment> segments;
};

// Fixed-capacity variant for hot path / no-heap usage
template <size_t MaxSegments>
struct TrajectoryFixed {
    std::array<Segment, MaxSegments> segments;
    size_t count;
};
```

### 10.7 Planner API

```cpp
struct TransferPlanRequest {
    BodyId departure_body;
    BodyId arrival_body;
    double departure_time;
    double arrival_time;
    LambertDirection direction;
};

struct TransferPlanResult {
    SolveStatus status;
    BodyId common_central_body;
    State2 transfer_departure_state;
    State2 transfer_arrival_state;
    Vec2 v_inf_departure;
    Vec2 v_inf_arrival;
    LambertSolution lambert;
};

class Planner {
public:
    TransferPlanResult solve_transfer(const TransferPlanRequest& req) const;
};
```

### 10.8 Serialization

All state structs (`State2`, `Segment`, `Trajectory`, `TrajectoryEvent`, `BodyDef`) shall be trivially serializable — POD layout with no pointers, no internal references, no padding-dependent semantics. This is required for authoritative server snapshots (§11).

---

## 11. Engine integration requirements

### 11.1 Runtime modes

The library shall support two runtime modes:

1. **Preview mode** — build full segment chain until time limit or event limit. Used for map display, maneuver planning, and client-side prediction.
2. **Lightweight query mode** — propagate only current segment or next event. Used for flight runtime and server tick processing.

### 11.2 Determinism

The library targets **same-platform bit-stable determinism** with **best-effort cross-platform agreement**. This matches the hybrid multiplayer model (§11.5).

#### Same-platform guarantees (required)

For two runs on the same binary, same CPU architecture, and same compiler flags, given identical inputs, all public APIs shall produce bit-identical outputs. This requires:

- fixed iteration order (no unordered containers in the physics path)
- no parallelism in the physics path
- no pointer-address-dependent ordering
- residual-based convergence everywhere; iteration caps return failure, never partial results
- fixed default tolerance constants (§13)
- deterministic tiebreakers for simultaneous events (§9.5.1)

#### Cross-platform behavior (best-effort)

Across different CPU architectures, compilers, or `libm` implementations, results may differ in the last few bits of mantissa due to:

- transcendental function precision (`sin`, `cos`, `sqrt`, `atan2`)
- FMA contraction
- x87 vs SSE vs NEON

These differences are acceptable for the hybrid multiplayer model because:

- the server is the source of truth for actual simulation state
- client-side previews are cosmetic predictions that get reconciled against server snapshots
- visible divergence in preview rendering across platforms is bounded and acceptable

Build flags for maximum cross-platform agreement (recommended but not enforced):

- `-fno-fast-math`
- `-ffp-contract=off`
- strict IEEE 754 mode
- avoid compiler intrinsics with implementation-defined precision

### 11.3 Performance targets

For a typical solar system of fewer than 100 bodies:

- body state query: O(depth)
- segment propagation: O(1)
- event search for current body children: O(num_children × event_cost)
- preview trajectory chain generation: bounded by configured max segments

Recommended defaults:

- `max_segments = 64`
- `max_event_refine_iterations = 32`
- `max_lambert_iterations = 32`

### 11.4 Memory

- stack-only single-query usage supported via `TrajectoryFixed<N>`
- optional caller-provided scratch buffers
- no mandatory global mutable state
- `std::vector`-based `Trajectory` is opt-in for planning code paths

### 11.5 Multiplayer integration model (hybrid)

The library is designed for the following multiplayer architecture:

**Server:**

- runs the authoritative simulation in lightweight query mode
- advances spacecraft state tick-by-tick
- detects SOI transitions and impacts authoritatively
- executes scheduled burns
- broadcasts state snapshots periodically (recommended: low rate during cruise, higher rate near burns and SOI transitions)

**Client:**

- receives authoritative state snapshots
- runs preview mode locally for trajectory visualization and planning UI
- predicts current state forward from last snapshot for responsive rendering
- on snapshot receipt: discards local segments from snapshot time forward, applies authoritative state, calls `rebuild_from(snapshot_time)` to regenerate the preview

**Determinism contract:**

- Server simulation is the single source of truth
- Same-platform determinism guarantees that a given client running multiple previews from the same snapshot produces stable visualization
- Cross-platform preview agreement is best-effort; the server's execution result is what actually applies
- Trajectory previews shown to the player may differ by sub-pixel amounts across client platforms, which is cosmetic

**Serialization requirements:**

- snapshots transmit current `State2`, current central `BodyId`, simulation time, and any pending scheduled burns
- the body system is established at session start and not transmitted per-snapshot
- segment chains are not transmitted; clients regenerate them locally

---

## 12. Failure handling

The library shall never crash or return partially uninitialized outputs for ordinary invalid game input. Every public solve/propagate/planning API shall return explicit status.

Failure categories:

- invalid body hierarchy
- invalid units or constants
- invalid initial state
- impossible event geometry
- Lambert no-solution
- numerical non-convergence
- max segment limit exceeded

Optional debug hooks:

- warning callback
- trace callback for planner decisions

---

## 13. Numerical tolerances

```cpp
struct Tolerances {
    double position_epsilon          = 1e-3;    // distance units
    double velocity_epsilon          = 1e-6;    // velocity units
    double angle_epsilon             = 1e-9;    // radians
    double time_epsilon              = 1e-6;    // seconds
    double root_epsilon              = 1e-6;    // for event refinement
    double lambert_residual_epsilon  = 1e-10;
    double parabolic_eccentricity_band = 1e-6;  // |e - 1| within this is ParabolaLike
};
```

Default values are part of the determinism contract — overriding them changes results.

Required behaviors:

- clamp `acos` and `asin` inputs to [-1, 1]
- treat tiny negative discriminants within epsilon as zero
- merge events within `time_epsilon`
- reject near-zero vector normalization

---

## 14. Simplifications specific to this game scope

Because all body orbits are fixed circles:

- ephemerides are analytic and exact within model
- no time-scale conversions are needed
- no frame ambiguity beyond parent-child inertial transforms
- SOI sizes can be authored rather than computed physically
- Lambert only needs same-central-body transfers

This allows aggressive simplification.

---

## 15. Known model limitations

The following artifacts are expected and acceptable:

- acceleration discontinuity at SOI boundary
- patched flybys are approximate
- no third-body perturbations near SOI edges
- no libration-point behavior
- long-term trajectory drift from reality if compared to n-body physics
- transfer feasibility may differ from full simulation if engine later adds perturbations
- sub-bit cross-platform divergence in client previews (cosmetic only under hybrid multiplayer model)

These are acceptable because the library is intended for gameplay and planning, not high-fidelity mission analysis.

---

## 16. Test requirements

### 16.1 Unit tests

Required coverage:

- circular ephemeris position/velocity against analytic expectation
- state/frame transforms across hierarchy
- ellipse propagation closed-orbit consistency
- hyperbola propagation asymptotic behavior
- SOI entry/exit detection (including grazing cases)
- impact detection
- patch state continuity (within `position_epsilon` / `velocity_epsilon`)
- post-patch re-entry suppression (no immediate re-patch loops)
- Lambert nominal cases
- Lambert near-degenerate cases (transfer angle near 0 and π)
- invalid-input handling (no NaNs, explicit failure status)

### 16.2 Property tests

Recommended:

- patch parent → child → parent without burn preserves absolute inertial state
- propagation to t1, then to t2, matches direct propagation to t2
- body-relative and ancestor-frame transformations are inverse-consistent
- same-platform: identical inputs produce bit-identical outputs across runs

### 16.3 Fuzz tests

- random valid body hierarchies + random initial states + random burn sequences
- assert: no NaNs, no crashes, all status codes valid, round-trip patch invariants hold

### 16.4 Golden scenario tests

Fixed scenarios with checked-in expected outputs:

- planet-to-planet transfer around star
- moon injection from planet orbit
- moon flyby and escape
- impact trajectory
- short-way vs long-way Lambert comparison (CCW vs CW direction)
- multi-segment preview chain through nested SOIs

### 16.5 Multiplayer integration tests

- server-side simulation produces snapshots; client rebuilds preview from snapshot and matches server's next snapshot to within tolerance
- snapshot reconciliation: inject artificial divergence on client, verify `rebuild_from` correctly resyncs
- platform matrix: run golden scenarios on x86-64, ARM64, and document observed cross-platform deltas

---

## 17. Recommended implementation phases

### Phase 1: core math
- Vec2
- body registry + builder
- circular ephemerides
- frame transforms
- two-body propagation
- conic classification

### Phase 2: events and patching
- impact detection
- SOI exit
- child SOI entry (with grazing handling and adaptive stepping)
- segment chaining
- post-patch re-entry suppression

### Phase 3: maneuvers
- impulsive burns
- downstream recomputation
- `rebuild_from(T)` operation
- preview trajectory builder

### Phase 4: Lambert
- single-rev solver
- direction selection
- transfer planner

### Phase 5: multiplayer integration
- serialization layout audit
- snapshot/reconcile reference implementation
- cross-platform test matrix
- determinism property tests

### Phase 6: engine polish
- visualization helpers
- debug traces
- tuning tolerances
- scenario tests

---

## 18. Optional extensions

Not required for v1, but compatible with this architecture:

- parking orbit helpers
- burn solver for circularize/escape/capture
- porkchop grid generation
- closest-approach search
- swingby turning-angle utilities
- authored transfer assists
- non-circular body ephemerides
- 3D upgrade path
- fixed-point math variant (if cross-platform bit-exactness ever becomes required)

---

## 19. Minimal v1 acceptance criteria

v1 is complete when the library can:

- define a star-planet-moon system via the builder API
- propagate a spacecraft conic around any body
- detect SOI entry/exit and impact (including grazing cases)
- patch correctly between frames without re-patch loops
- apply an impulsive burn
- build a preview chain of future segments
- solve a same-central-body Lambert transfer
- generate planet-to-planet transfer velocities in the parent frame
- serialize and restore spacecraft state for snapshot-based multiplayer
- rebuild a preview chain from an injected authoritative snapshot
- produce bit-identical results across runs on the same platform
- do all of the above without NaNs or crashes on invalid input

---

## 20. Example usage flows

### A. Preview current trajectory (client)

1. spacecraft has current body, state, time
2. build segment
3. detect earliest event
4. patch or terminate
5. repeat until horizon or max segments
6. render segments

### B. Plan a transfer

1. pick departure body and arrival body
2. choose departure and arrival times
3. query both body states in common parent frame
4. solve Lambert
5. compute departure v_inf
6. convert to burn outside core if needed
7. inject burn into preview chain

### C. Multiplayer snapshot reconciliation (client)

1. receive authoritative snapshot `(state, central_body, time)` from server
2. discard local segments with `start_time >= snapshot.time`
3. install snapshot state as the initial state of a new segment
4. call `rebuild_from(snapshot.time)` to regenerate forward preview
5. continue local prediction from new chain

### D. Server tick

1. for each spacecraft: advance current segment to tick time
2. check for events (impact, SOI transition, scheduled burn)
3. apply event if reached
4. on snapshot interval: serialize spacecraft state, broadcast to clients

---

## 21. Non-goals for v1

The following shall be deferred:

- automatic multi-body search
- gravity-assist search
- low-thrust
- optimizer integration
- multi-revolution Lambert
- continuous collision with arbitrary shapes
- station-keeping or thrust arcs
- bit-exact cross-platform determinism
- fixed-point math
