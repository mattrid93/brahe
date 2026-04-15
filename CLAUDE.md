# Brahe

2D patched-conics orbital mechanics library for game engine integration. C++17, CMake.

## Build

```sh
cmake --preset default
cmake --build --preset default
```

## Test

```sh
ctest --preset default
```

## Project layout

- `include/brahe/types.h` — POD types, enums, BodyId, Vec2, State2, Tolerances, BodyDef
- `include/brahe/vec2.h` — inline Vec2 math operations
- `include/brahe/body_system.h` — BodySystemBuilder (mutable) and BodySystem (immutable)
- `include/brahe/conics.h` — ConicElements2D, ConicState (POD)
- `include/brahe/two_body.h` — TwoBody class (invariants, classification, conversion, propagation), detail anomaly helpers
- `src/body_system.cpp` — builder validation, ephemeris, frame transforms
- `src/two_body.cpp` — conic math, anomaly solvers, branch-based propagation
- `src/main.cpp` — demo executable
- `tests/` — Catch2 tests (102 total across 10 files)
- `specs/` — spec and phase plans

## Conventions

- C++17, double precision throughout
- Formatting: `.clang-format` (Google-based, 4-space indent, 100 col)
- Static analysis: `.clang-tidy`
- Warnings: `-Wall -Wextra -Wpedantic`
- All physics types are trivially copyable and standard layout (POD)
- Body storage is canonically sorted by ascending BodyId
- Query-path functions do not heap-allocate
- Angles in (-pi, pi], circular convention: omega=0 nu=atan2(y,x)
- Kepler solvers: bounded Newton iteration, cap = failure not partial result
- Convergence tolerance scales with problem magnitude for numerical stability
