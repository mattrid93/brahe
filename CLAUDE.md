# Brahe

C++ project using CMake.

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

- `include/brahe/` — public headers
- `src/` — library and executable sources
- `tests/` — Catch2 tests

## Conventions

- C++17
- Formatting: `.clang-format` (Google-based, 4-space indent, 100 col)
- Static analysis: `.clang-tidy`
- Warnings: `-Wall -Wextra -Wpedantic`
