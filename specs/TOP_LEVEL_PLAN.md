# Patched Conics Library — Top-Level Plan

## 1. Foundation and world model
- Core POD types and enums
- `Vec2` / `State2` math primitives
- `Tolerances`
- `BodyDef`, `BodySystemBuilder`, immutable `BodySystem`
- body validation
- circular ephemerides
- ancestor/root frame transforms
- deterministic body ordering and cached hierarchy metadata

## 2. Two-body conic math
- Cartesian ↔ invariants / 2D conic elements
- conic classification
- anomaly conversion helpers
- analytic propagation for ellipse, hyperbola, parabola-like cases
- numerical guardrails

## 3. Event detection kernel
- impact root functions
- SOI exit root functions
- child SOI entry search
- coarse scan + bounded refinement
- grazing handling
- deterministic event ordering and tie-breaks

## 4. Patching and segment chaining
- `Segment`, `Trajectory`, `TrajectoryFixed`
- parent→child and child→parent patch generation
- state continuity checks
- post-patch re-entry suppression
- preview chain builder

## 5. Maneuvers and recomputation
- `ImpulsiveBurn`
- burn-frame transforms
- segment termination on burn
- downstream recomputation
- `rebuild_from(T)`

## 6. Lambert and transfer planning
- 2D single-rev Lambert solver
- CCW/CW branch handling
- bounded convergence / fallback bracketing
- `Planner::solve_transfer`

## 7. Serialization and determinism hardening
- POD layout audit
- fixed traversal rules
- tolerance centralization
- same-platform determinism tests
- no-heap hot-path verification for fixed-capacity mode

## 8. Scenario, fuzz, and multiplayer integration
- golden trajectory cases
- fuzz/property tests
- snapshot/reconcile reference flow
- platform-delta documentation
