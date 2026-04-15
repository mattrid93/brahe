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
- `src/body_system.cpp` — builder validation, ephemeris, frame transforms
- `src/main.cpp` — demo executable
- `tests/` — Catch2 tests (test_types, test_body_builder, test_ephemeris, test_frame_transforms)
- `specs/` — spec and phase plans

## Conventions

- C++17, double precision throughout
- Formatting: `.clang-format` (Google-based, 4-space indent, 100 col)
- Static analysis: `.clang-tidy`
- Warnings: `-Wall -Wextra -Wpedantic`
- All physics types are trivially copyable and standard layout (POD)
- Body storage is canonically sorted by ascending BodyId
- Query-path functions do not heap-allocate
