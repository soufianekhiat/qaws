# Contributing to Qaws

Thank you for your interest in contributing to Qaws.

## Building

Qaws uses **Sharpmake** (primary) or **CMake** (backup).

```bash
# CMake
cmake -B build -DQAWS_SCALAR_IS_FLOAT=ON
cmake --build build --config Release
ctest --build-config Release --output-on-failure

# Sharpmake
sharpmake /sources(@'buildsystem/sharpmake/src/main.sharpmake.cs')
```

### Build options

| Option | Default | Effect |
|---|:---:|---|
| `QAWS_SCALAR_IS_FLOAT` | `ON` | `float` (32-bit) or `double` (64-bit) |
| `QAWS_ENABLE_SIMD` | `OFF` | AVX2 / SSE2 / NEON batch eval |
| `QAWS_ENABLE_WARNINGS` | `ON` | Compiler warnings |
| `QAWS_WARNINGS_AS_ERRORS` | `ON`(f64) / `OFF`(f32) | `-Werror` / `/WX` |

## Code style

- C11 standard, no compiler extensions (`C_EXTENSIONS NO`)
- Public API prefixed with `qaws_`
- All fallible functions return `qaws_status`
- Tab indentation in source files
- No external dependencies — C11 standard library only (`<math.h>`, `<stdlib.h>`, `<string.h>`)

## Testing

Tests live in `tests/` and follow the `XX_name.c` naming convention. Each file exports a single `test_XX_name_main(void)` function that returns 0 on success or 1 on failure. The `test_runner.c` file calls all test mains sequentially.

To add a new test:

1. Create `tests/XX_name.c` (pick the next available index)
2. Include `test_common.h`
3. Write test functions using `TEST_ASSERT` / `TEST_ASSERT_STATUS`
4. Export a `test_XX_name_main(void)` function that calls them
5. Add a forward declaration and registry entry in `test_runner.c`
6. Add the file to `tests/CMakeLists.txt`

Run the suite:

```bash
cmake --build build --target check
```

## Adding a curve family

1. Create `src/qaws/qaws_<family>.h` with the descriptor struct and `qaws_curve_create_<family>`
2. Create `src/qaws/qaws_<family>.c` implementing creation and the eval vtable
3. Add `QAWS_CURVE_KIND_<FAMILY>` to `qaws_types.h`
4. Include the header in `qaws.h`
5. Add the `.c` file to `QAWS_SOURCES` in `CMakeLists.txt`
6. Add the `.h` to the install list in `CMakeLists.txt`
7. Write tests and documentation

## Pull requests

- One logical change per PR
- Tests must pass for both f32 and f64 configurations
- New public API functions need documentation in `docs/api/`

## License

By contributing you agree that your contributions will be licensed under the [MIT License](LICENSE).
