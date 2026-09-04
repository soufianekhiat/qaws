# Integrating Qaws Without CMake

Qaws is dependency-free C11. It needs no build system of its own, and nothing
stops you from compiling it as part of whatever target you already have. There
are three ways in, from least to most self-contained.

1. [Drop the tree into your project](#1-drop-the-tree-into-your-project) --
   add the `.c` files to an existing Xcode / Visual Studio / Makefile target.
2. [Single translation unit](#2-single-translation-unit) -- compile the whole
   library as one `.c` for faster builds and better inlining.
3. [Amalgamation](#3-amalgamation) -- generate one `qaws.h` and one `qaws.c`
   and forget the rest of the repository exists.

See [building.md](building.md) for the CMake and Sharpmake paths.

## 1. Drop the tree into your project

Add every `.c` file under `src/qaws/` to your target, add `src/qaws` to the
header search path, and compile as C11. That is the whole integration.

| Setting | Value | Notes |
|---|---|---|
| Sources | `src/qaws/**/*.c` | 57 files, plus 4 more under `batch/` for SIMD |
| Header search path | `src/qaws` | Headers include each other relatively |
| Language standard | C11 | `_Static_assert`, anonymous structs |
| Required define | `QAWS_SCALAR_IS_FLOAT` | `1` for `float`, `0` for `double` |
| Recommended flag | `-fno-strict-aliasing` | GCC/Clang; CMake sets it by default |
| Link | `m` | POSIX only. Nothing to add on Apple or Windows |

`QAWS_SCALAR_IS_FLOAT` defaults to `1` if you do not define it, but define it
explicitly: it changes the ABI of every struct in the public API, so it must
match everywhere it is seen.

### Xcode and iOS

Drag `src/qaws` into the project, let Xcode add the `.c` files to your target,
and set **Header Search Paths** to the `src/qaws` directory. Set
`QAWS_SCALAR_IS_FLOAT=1` under **Preprocessor Macros**. Nothing needs linking:
`libm` is part of `libSystem` on Apple platforms.

**SIMD needs no flags there.** `src/qaws/simd/qaws_simd.h` selects its ISA from
macros the compiler already predefines -- `__AVX2__`, `__SSE2__`, `__ARM_NEON`
-- not from the `QAWS_SIMD_*` defines that CMake sets. Clang predefines
`__ARM_NEON` for every arm64 target, so on Apple silicon and iOS you get the
NEON path simply by adding `src/qaws/batch/*.c` to the target. The
`QAWS_SIMD_*` defines only exist so the CMake and Sharpmake builds can report
which ISA they probed; they do not drive the dispatch.

On x86-64 the same header falls back to SSE2 unconditionally (`_M_X64` /
`__SSE2__` are always defined there) and only reaches the AVX2 path when you
actually pass `/arch:AVX2` or `-mavx2 -mfma`.

If you leave `batch/` out, everything still builds; you just lose the batch
evaluation entry points.

## 2. Single translation unit

Compiling the library as one translation unit is faster than 57 separate ones
and lets the compiler inline across what used to be file boundaries. Either
write the unity file yourself:

```c
/* qaws_unity.c */
#include "qaws_status.c"
#include "qaws_curve.c"
/* ... every other .c under src/qaws/ ... */
```

or let CMake do it:

```bash
cmake -B build -DQAWS_UNITY_BUILD=ON
```

The tree is kept free of the symbol collisions that normally make this fail,
and CI compiles this layout on every supported platform.

## 3. Amalgamation

The most self-contained option: collapse the entire library into one header and
one implementation file.

```bash
python buildsystem/amalgamate.py --output-dir amalgam        # scalar only
python buildsystem/amalgamate.py --output-dir amalgam --simd # with batch/
```

That writes `amalgam/qaws.h` (~79 KiB) and `amalgam/qaws.c` (~1 MiB). Those two
files are the whole library. Compile them with no include paths, no library and
no build system:

```bash
cc -std=c11 -DQAWS_SCALAR_IS_FLOAT=1 -Iamalgam \
   examples/amalgam_standalone.c amalgam/qaws.c -lm -o smoke
```

```bat
cl /std:c11 /TC /DQAWS_SCALAR_IS_FLOAT=1 /Iamalgam ^
   examples\amalgam_standalone.c amalgam\qaws.c /Fe:smoke.exe
```

`examples/amalgam_standalone.c` is a complete worked example. Drag `qaws.h` and
`qaws.c` into an Xcode project and it behaves like any other pair of files.

The generator inlines quoted includes recursively and emits each file once,
which matches the include-guard semantics of the original tree. Angle-bracket
includes are left where they are, because several of them sit inside ISA
detection blocks and have to stay conditional.

### From CMake

```bash
# Just generate the pair, build nothing.
cmake -B build && cmake --build build --target amalgamate

# Build a static library from the generated qaws.c.
cmake -B build -DQAWS_AMALGAMATION=ON

# ... and run the full test suite against it.
cmake -B build -DQAWS_TEST_AMALGAMATION=ON && ctest --test-dir build
```

| Option | Default | Description |
|---|---|---|
| `QAWS_UNITY_BUILD` | `OFF` | Compile the library as one translation unit |
| `QAWS_AMALGAMATION` | `OFF` | Generate and build `qaws.h` + `qaws.c` |
| `QAWS_TEST_AMALGAMATION` | `OFF` | Point the test suite at the amalgamated library (implies the above) |
| `QAWS_AMALGAMATION_DIR` | `<build>/amalgam` | Where the generated pair is written |

`QAWS_AMALGAMATION` follows `QAWS_ENABLE_SIMD`: with SIMD on, the `batch/`
sources are folded in too. Generating requires a Python 3 interpreter;
consuming the result does not.

### From Sharpmake

`generate_projects.bat` adds a `QawsAmalgam` static library project alongside
`Qaws` and `QawsTests`. It generates both variants up front, compiles the one
matching each configuration's `SimdMode`, and regenerates on every build so the
pair cannot go stale after an edit under `src/qaws/`. If no Python 3
interpreter is found the project is skipped with a warning and the rest of the
solution generates as before.

## What CI checks

`.github/workflows/permutations.yml` builds every combination of OS (Linux,
macOS, Windows) x scalar type (f32, f64) x SIMD (on, off) x layout (per-file,
unity, amalgamation), and separately compiles and runs
`examples/amalgam_standalone.c` against a freshly generated amalgamation on all
three platforms. The f32 amalgamation is uploaded as a build artifact, so you
can download `qaws.h` and `qaws.c` without running the generator at all.
