# Integrating Qaws Without a Build System

Qaws is dependency-free C11 and can be compiled as part of any target. Three
ways in, from least to most self-contained:

1. [Add the sources to your project](#1-add-the-sources-to-your-project).
2. [One translation unit per precision](#2-one-translation-unit-per-precision).
3. [Amalgamation](#3-amalgamation): one `qaws.h` and one `qaws.c`.

See [building.md](building.md) for the Sharpmake and CMake builds.

## The one rule: compile each source twice

The library holds both precisions ([precision.md](precision.md)): every
`.c` file is compiled once with `QAWS_SCALAR_IS_FLOAT` 1 (the `_f32` names)
and once with 0 (the `_f64` names). The precision-free sources
(`qaws_status.c`, `internal/qaws_internal_memory.c`,
`internal/qaws_internal_cpu.c`) guard their bodies and define them in the 0
pass only, so compiling them twice is harmless. Compile each source through a
small wrapper per precision:

```c
/* qaws_bezier_f32.c */
#undef QAWS_SCALAR_IS_FLOAT
#define QAWS_SCALAR_IS_FLOAT 1
#include "path/to/src/qaws/qaws_bezier.c"
```

and the same with `0` for `qaws_bezier_f64.c`. This is what the CMake and
Sharpmake builds generate. Your own code then defines `QAWS_SCALAR_IS_FLOAT`
(default 1) only to choose what its unsuffixed names mean; it does not change
the library.

## 1. Add the sources to your project

| Setting | Value | Notes |
|---|---|---|
| Sources | wrappers for every `.c` under `src/qaws/`, `src/qaws/internal/` and `src/qaws/exact/` | add `src/qaws/batch/*.c` for the SIMD kernels |
| Header search path | `src/qaws` | headers include each other relatively |
| Language standard | C11 | |
| Exact kernel | define `QAWS_ENABLE_EXACT=1` | as the Sharpmake and CMake builds do; the 2D path code then decides orientations exactly |
| Recommended flag | `-fno-strict-aliasing` | GCC and Clang; CMake sets it by default |
| Link | `m` | POSIX only; nothing on Apple or Windows |

### Xcode and iOS

Add the wrappers to the target, set **Header Search Paths** to `src/qaws`,
and set the caller's `QAWS_SCALAR_IS_FLOAT` under **Preprocessor Macros**.
`libm` is part of `libSystem` on Apple platforms.

### SIMD

`src/qaws/simd/qaws_simd.h` selects its ISA from macros the compiler
predefines (`__AVX2__`, `__SSE2__` / `_M_X64`, `__ARM_NEON`), not from the
`QAWS_SIMD_*` defines the build systems set. Clang predefines `__ARM_NEON`
for arm64, so Apple silicon and iOS get the NEON path by adding
`src/qaws/batch/*.c`. On x86-64 SSE2 is always available; AVX2 needs
`/arch:AVX2` or `-mavx2 -mfma`.

The `batch/` kernels have no public header and no public entry point calls
them; leaving `batch/` out changes nothing in the public API.

## 2. One translation unit per precision

Write one unity file per precision, each including every source after
setting the precision:

```c
/* qaws_unity_f32.c */
#undef QAWS_SCALAR_IS_FLOAT
#define QAWS_SCALAR_IS_FLOAT 1
#include "qaws_status.c"
#include "qaws_curve.c"
/* ... every other .c ... */
```

Both precisions cannot share one unit: each would define every static
function twice. CMake does the same grouping with:

```bash
cmake -B build -DQAWS_UNITY_BUILD=ON
```

## 3. Amalgamation

```bash
cmake -DQAWS_OUTPUT_DIR=amalgam -P buildsystem/amalgamate.cmake                 # without batch/
cmake -DQAWS_OUTPUT_DIR=amalgam -DQAWS_SIMD=ON -P buildsystem/amalgamate.cmake  # with batch/
```

This writes `amalgam/qaws.h` (every public header inlined) and
`amalgam/qaws.c` (the internal headers and every source). The generator
inlines quoted includes recursively, emits each file at most once, and leaves
angle-bracket includes in place (several sit inside ISA detection blocks).
The bootstrap headers `qaws_hlsl.h`, `qaws_glsl.h`, `qaws_halide.h` are left
out.

Compile `qaws.c` twice, once per precision, as the banner of the generated
files says:

```bash
cc -std=c11 -c -DQAWS_SCALAR_IS_FLOAT=1 -Iamalgam amalgam/qaws.c -o qaws_f32.o
cc -std=c11 -c -DQAWS_SCALAR_IS_FLOAT=0 -Iamalgam amalgam/qaws.c -o qaws_f64.o
```

Including `qaws.h` declares the precision `QAWS_SCALAR_IS_FLOAT` selects
(default 1); include it again under the other value, or use the per-precision
guards as `qaws_both.h` does, for both. `examples/amalgam_standalone.c` is a
minimal consumer.

### From CMake

```bash
cmake -B build && cmake --build build --target amalgamate   # generate only
cmake -B build -DQAWS_AMALGAMATION=ON                       # also build qaws_amalgam (Qaws::amalgam)
cmake -B build -DQAWS_TEST_AMALGAMATION=ON                  # and link the tests against it
```

The `qaws_amalgam` target compiles `qaws.c` through one wrapper per
precision. `QAWS_AMALGAMATION` follows `QAWS_ENABLE_SIMD`. Generating needs
CMake (`cmake -P`); consuming the pair does not.

### From Sharpmake

Generation adds a `QawsAmalgam` static library next to `Qaws` and
`QawsTests`. Both variants (`tmp/amalgam/simd/` and `tmp/amalgam/nosimd/`)
are generated up front; each configuration compiles the one matching its SIMD
mode, through one wrapper per precision, and regenerates it before every
build. Without `cmake` on PATH the project is skipped with a warning.

## What CI checks

`.github/workflows/permutations.yml` builds and tests every combination of
OS (Linux, macOS, Windows), caller precision (f32, f64), SIMD (on, off) and
layout (per-file, unity, amalgamation with the tests linked against it).
