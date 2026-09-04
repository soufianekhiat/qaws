# Building Qaws

Qaws has two build systems: CMake (cross-platform) and Sharpmake (Visual Studio).

## CMake

### Basic build

```bash
cd qaws
mkdir build && cd build
cmake ..
cmake --build .
```

### Configuration options

| Option | Default | Description |
|---|---|---|
| `QAWS_SCALAR_IS_FLOAT` | `ON` | Use `float` (f32). Set to `OFF` for `double` (f64). |
| `QAWS_ENABLE_WARNINGS` | `ON` | Enable compiler warnings |
| `QAWS_WARNINGS_AS_ERRORS` | `ON` (f64) / `OFF` (f32) | Treat warnings as errors |
| `QAWS_DISABLE_STRICT_ALIASING` | `ON` | Pass `-fno-strict-aliasing` on GCC/Clang |
| `QAWS_ENABLE_SIMD` | `OFF` | Enable SIMD batch evaluation. Auto-detects the best available ISA at configure time (AVX2, SSE2, or NEON). When enabled, the `batch/` source files are compiled into the library. |
| `QAWS_DEBUG_POSTFIX` | `"d"` | String appended to the library name in Debug builds (e.g. `qawsd`). Set to `""` to disable. |
| `QAWS_DEBUG_PREFIX` | `""` | String prepended to the library name in Debug builds (e.g. `dbg_qaws`). |
| `QAWS_INSTALL` | `ON` | Enable install rules |
| `QAWS_INSTALL_PKGCONFIG` | `ON` | Install the `qaws.pc` pkg-config file |
| `BUILD_TESTING` | `ON` | Build test suite |
| `QAWS_UNITY_BUILD` | `OFF` | Compile the library as a single translation unit |
| `QAWS_AMALGAMATION` | `OFF` | Generate and build the `qaws.h` + `qaws.c` amalgamation. Requires Python 3 |
| `QAWS_TEST_AMALGAMATION` | `OFF` | Run the test suite against the amalgamated library. Implies `QAWS_AMALGAMATION` |
| `QAWS_AMALGAMATION_DIR` | `<build>/amalgam` | Where the generated `qaws.h`/`qaws.c` are written |

See [integration.md](integration.md) for the single-file and no-build-system
paths, including the `amalgamate` target and the Xcode/iOS notes.

### f32 vs f64 build

```bash
# Float (default)
cmake .. -DQAWS_SCALAR_IS_FLOAT=ON

# Double
cmake .. -DQAWS_SCALAR_IS_FLOAT=OFF
```

### SIMD build

```bash
cmake .. -DQAWS_ENABLE_SIMD=ON
```

When `QAWS_ENABLE_SIMD` is `ON`, CMake probes the compiler for AVX2, SSE2, and
NEON support (in that order) and enables the first one it finds. The
corresponding define (`QAWS_SIMD_AVX2`, `QAWS_SIMD_SSE2`, or `QAWS_SIMD_NEON`)
is set automatically and the appropriate compiler flags are added (e.g.
`/arch:AVX2` on MSVC, `-mavx2 -mfma` on GCC/Clang).

### Running tests

```bash
cmake --build . --target check
# or
ctest --output-on-failure
# or run directly
./tests/Debug/qaws_tests      # Windows
./tests/qaws_tests             # Linux/macOS
```

Tests produce SVG files in `working_dir/tests_output/` for visual verification of curve shapes.

### Installing

```bash
cmake --install . --prefix /usr/local
```

Installs:
- Library: `lib/libqaws.a` (or `.so` / `.dll`)
- Public headers: `include/qaws*.h` (all public API headers)
- Core algorithm headers: `include/core/qaws_*_core.h` (backend-agnostic evaluation kernels)
- Bootstrap headers: `include/qaws_hlsl.h`, `include/qaws_glsl.h`, `include/qaws_halide.h`
- CMake config: `lib/cmake/Qaws/`
- pkg-config: `lib/pkgconfig/qaws.pc` (when `QAWS_INSTALL_PKGCONFIG` is `ON`)

### Using from another CMake project

```cmake
find_package(Qaws REQUIRED)
target_link_libraries(myapp PRIVATE Qaws::qaws)
```

## Multi-backend note

Qaws is a multi-backend library. The C build (CMake or Sharpmake) produces the
runtime library used by C/C++ callers. The evaluation algorithms live in
backend-agnostic `core/*.h` headers that contain pure arithmetic with no
platform dependencies. HLSL, GLSL, and Halide users do **not** link the C
library -- they include the corresponding bootstrap header (`qaws_hlsl.h`,
`qaws_glsl.h`, or `qaws_halide.h`) which pulls in the core headers directly.
The backend is therefore selected at include time, not at build time.

## Sharpmake (Visual Studio)

### Prerequisites

- Visual Studio 2022
- .NET SDK 6.0 or later
- Git submodules initialized

### Bootstrap

```bash
git submodule update --init --recursive
cd buildsystem
bootstrap.bat
```

This will:
1. Build Sharpmake from source
2. Copy Sharpmake binaries to `buildsystem/sharpmake/`
3. Generate Visual Studio solution

### Opening the solution

After bootstrap, open `Qaws_vs2022_win64.sln` in the repo root.

### Configuration names

Sharpmake generates configurations for every combination of `Optimization`,
`ScalarType` (Float / Double), and `SimdMode` (On / Off). The resulting
configuration names follow the pattern:

```
{Optimization}_{scalar}[_simd]
```

For example: `Debug_f32`, `Release_f64`, `Debug_f32_simd`, `Release_f64_simd`.

When `SimdMode` is `On`, the build defines `QAWS_SIMD_AVX2=1`, adds `/arch:AVX2`,
and compiles the `batch/` source files. When `SimdMode` is `Off`, the `batch/`
sources are excluded via a build regex filter.

### Regenerating projects

```bash
cd buildsystem
generate_projects.bat
```

### Amalgamation project

Generation also adds a `QawsAmalgam` static library built from the generated
single-file `qaws.c`. It produces the SIMD and non-SIMD variants up front,
compiles whichever matches each configuration's `SimdMode`, and regenerates
before every build. It is skipped with a warning if no Python 3 interpreter is
found. See [integration.md](integration.md).

## Compiler compatibility

| Compiler | Status |
|---|---|
| MSVC 2022 (cl) | Tested, /W4 clean |
| GCC 11+ | Supported |
| Clang 13+ | Supported |

The library requires C11 (`_Static_assert`, anonymous structs, etc.).
