# Building Qaws

Sharpmake is the primary build system (Visual Studio 2022, Windows). CMake
mirrors the public API for other platforms, installs, packaging and CI. For
builds without either, see [integration.md](integration.md).

## Both precisions in one library

The library holds every precision-bearing function and type twice, as
`name_f32` and `name_f64` ([precision.md](precision.md)). Both build systems
compile each `.c` file through two generated wrappers, one per precision,
never directly:

```c
/* generated: <source> as f32 */
#undef QAWS_SCALAR_IS_FLOAT
#define QAWS_SCALAR_IS_FLOAT 1
#include "<absolute path to the source>"
```

The precision-free sources (the status strings, the memory module with the
default, thread and scratch allocators and the arena) define their symbols in
the f64 pass only, so the allocator state is one for both precisions.
Wrappers go to `<build>/qaws_precision/` (CMake) or `tmp/qaws_precision/`
(Sharpmake).

`QAWS_SCALAR_IS_FLOAT` (default 1) then only selects what a caller's
unsuffixed names, `qaws_scalar` and the scalar macros mean. A translation unit
that wants both sets side by side includes `qaws_both.h`.

After adding or renaming a public declaration, regenerate
`src/qaws/qaws_names.h` with `tools/qaws_api_tool` (see CONTRIBUTING.md);
the CTest test `qaws_api_names` fails while it is stale.

## Sharpmake (primary)

### Prerequisites

- Visual Studio 2022
- .NET SDK 6.0 or later
- Git submodules initialized (`extern/Sharpmake`)

### Bootstrap

```bat
git submodule update --init --recursive
buildsystem\bootstrap.bat
```

`bootstrap.bat` builds Sharpmake from `extern/Sharpmake`, copies its binaries
to `buildsystem/sharpmake/` and runs `generate_projects.bat`, which writes
`Qaws_vs2022_win64.sln` at the repository root. Regenerate after adding or
removing sources:

```bat
buildsystem\generate_projects.bat
```

The project sources are `buildsystem/sharpmake/src/` (`main.cs`, `common.cs`,
`QawsLib.cs`, `QawsTests.cs`, `QawsAmalgam.cs`).

### Projects

| Project | Output | Contents |
|---|---|---|
| `Qaws` | static library | every `.c` under `src/qaws/` (subdirectories included), through the per-precision wrappers |
| `QawsTests` | console executable in `working_dir/` | every `.c` under `tests/` except `debug_*.c` |
| `QawsAmalgam` | static library | the generated `qaws.c`, when `cmake` is on PATH (see [integration.md](integration.md)) |

### Configurations

Projects are generated for every combination of optimization (Debug,
Release), scalar (Float, Double) and SIMD (Off, On), named
`{Optimization}_{f32|f64}[_simd]`, for example `Debug_f32`,
`Release_f64_simd`. The solution lists `Debug_f32`, `Debug_f64`,
`Release_f32`, `Release_f64`.

- The scalar only sets `QAWS_SCALAR_IS_FLOAT` for callers (the tests): the
  library holds both precisions in every configuration.
- f32 executables get the suffix `_f32`, so both precisions share
  `working_dir/` without overwriting each other.
- Every configuration compiles as C (`/TC /std:c11`) with `/W4 /WX` and
  `/fp:precise`.
- SIMD On defines `QAWS_SIMD_AVX2=1`, adds `/arch:AVX2` and compiles
  `src/qaws/batch/`; SIMD Off excludes `batch/`.
- Intermediates go to `tmp/<project>/<configuration>/`, libraries to
  `tmp/lib/`.

## CMake

### Basic build

```bash
cmake -B build
cmake --build build --config Release
ctest --test-dir build --build-config Release --output-on-failure
```

CMake 3.20 or later. The project version is 3.0.0 (`QAWS_VERSION_STRING`
in `qaws.h`).

### Options

| Option | Default | Description |
|---|---|---|
| `QAWS_SCALAR_IS_FLOAT` | `ON` | What callers' unsuffixed names mean: `ON` f32, `OFF` f64 (an interface definition of the target). The library holds both. |
| `QAWS_ENABLE_WARNINGS` | `ON` | `/W4` or `-Wall -Wextra -Wpedantic` |
| `QAWS_WARNINGS_AS_ERRORS` | `ON` with f64, `OFF` with f32 | `/WX` or `-Werror` |
| `QAWS_DISABLE_STRICT_ALIASING` | `ON` | `-fno-strict-aliasing` on GCC and Clang |
| `QAWS_ENABLE_SIMD` | `OFF` | Probe AVX2, SSE2, NEON (in that order), define the first found (`QAWS_SIMD_AVX2`, `_SSE2`, `_NEON`) and compile `batch/` |
| `QAWS_ENABLE_EXACT` | `ON` | Compile the exact kernel (`src/qaws/exact/`) and define `QAWS_ENABLE_EXACT=1` publicly |
| `QAWS_UNITY_BUILD` | `OFF` | One translation unit per precision |
| `QAWS_AMALGAMATION` | `OFF` | Generate and build the `qaws.h` + `qaws.c` amalgamation |
| `QAWS_TEST_AMALGAMATION` | `OFF` | Link the tests against the amalgamation (implies `QAWS_AMALGAMATION`) |
| `QAWS_AMALGAMATION_DIR` | `<build>/amalgam` | Where the generated pair is written |
| `QAWS_BUILD_EXAMPLES` | `OFF` | Build the C showcases in `examples/` into `working_dir/` |
| `QAWS_DEBUG_POSTFIX` | `d` | Appended to the library name in Debug (`qawsd`) |
| `QAWS_DEBUG_PREFIX` | empty | Prepended to the library name in Debug |
| `QAWS_INSTALL` | `ON` | Install rules |
| `QAWS_INSTALL_PKGCONFIG` | `ON` | Install `qaws.pc` |
| `BUILD_TESTING` | `ON` | Build `qaws_tests` and register the CTest tests |

The library target is `qaws` (alias `Qaws::qaws`), C11 without extensions,
position-independent; it links `m` on non-Windows platforms. With
`BUILD_SHARED_LIBS` on Windows every symbol is exported.

### Tests

```bash
cmake --build build --target check   # build qaws_tests, then ctest --output-on-failure
ctest --test-dir build -R "^14"      # suites whose name starts with 14
ctest --test-dir build -j 4          # several suites at once
```

CTest registers one test per suite, named after its file (`85_shatter`,
`142_curve_view`, ...), plus `qaws_api_convention` and `qaws_api_names`
(the naming tool's checks). `qaws_tests` is written to `working_dir/` and
runs from there; run by hand it runs every suite in turn and reports each
one's time. Environment variables:

| Variable | Effect |
|---|---|
| `QAWS_TEST=text` | Run only the suites whose name contains `text` |
| `QAWS_TEST==name` | Run only the suite named `name` |
| `QAWS_TEST_HEAVY=1` | Also run the slow checks (large meshes, robustness sweeps, the smallest tolerances) |
| `QAWS_TEST_SCRATCH=1` | Give the main thread's scratch an arena that must be empty between suites |

Suites write SVG and OBJ files to `working_dir/tests_output/` for visual
checks.

### Installing

```bash
cmake --install build --prefix /usr/local
```

Installs the library, the public headers (`include/qaws*.h`, including
`qaws_names.h`, `qaws_both.h`, `qaws_precision.h` and the bootstrap headers
`qaws_hlsl.h`, `qaws_glsl.h`, `qaws_halide.h`), the core headers
(`include/core/`), the CMake package (`lib/cmake/Qaws/`, version
compatibility SameMajorVersion) and `lib/pkgconfig/qaws.pc`.

The install list does not contain `qaws_exact.h` and `qaws_motion.h`, which
`qaws.h` includes, nor the core headers `qaws_frame_core.h`,
`qaws_curve_record_core.h` and `qaws_surface_record_core.h`. Use the source
tree's `src/qaws` as include path when you need them.

```cmake
find_package(Qaws REQUIRED)
target_link_libraries(myapp PRIVATE Qaws::qaws)
```

## Shader and Halide backends

The C build produces the runtime library for C and C++ callers. HLSL, GLSL
and Halide users do not link it: they include a bootstrap header and the
`core/` headers they need. The backend is chosen at include time
([backends.md](backends.md)).

## Compilers

| Compiler | Notes |
|---|---|
| MSVC 2022 | Sharpmake and CMake; `/W4 /WX` |
| GCC, Clang | CMake; CI builds Linux and macOS |

The library requires C11.
