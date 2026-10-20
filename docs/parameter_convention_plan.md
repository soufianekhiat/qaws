# Parameter convention: aligning Qaws with Alwan, Suwar and Abaad

Alwan's parameter convention (v2.0, `CONTRIBUTING.md` there) is shared by
Suwar and Abaad: outputs first, a buffer's stride next to it, the context
last. Qaws grew the opposite order: the context or object first, then
inputs, then outputs. This plan moves all of Qaws to the shared convention.

## The rules, as Qaws will apply them

| # | Rule | Source |
|---|---|---|
| R1 | The context (`qaws_diff_context`, any future `qaws_ctx`) is the last parameter. | Alwan R1 |
| R2 | Each `*_stride` follows the buffer it strides. | Alwan R2 |
| R3 | Outputs come before inputs. | Alwan R3 |
| R4 | Sizes (`count`, `*_count`, `dim`, `width`, `height`) come after the buffer block. | Alwan R4 |
| R5 | Pure value maths (return by value, no pointers) takes no context. | Alwan R5 |

Alwan's rules cover buffers. Qaws also has objects (curves, solvers, trees),
so these follow Abaad's existing practice:

- **Order of a full signature:** plain outputs; then summaries or reports
  (`qaws_repulsive_report*`); then in-out buffers (positions updated in
  place); then the object handle (`qaws_repulsive_curves_solver*`); then
  input buffers, each followed by its stride; then sizes; then the
  descriptor or options struct; then enums and scalars; the context last.
  - Abaad: `abaad_kdtree_knn_f32(indices, dist2, found, tree, ...)`.
  - Abaad: `abaad_lm_solve_f32(summary, params, param_count, residuals, ..., options, ctx)`.
- **Constructors:** the new object first:
  `qaws_repulsive_curves_create(&curves, &desc)`. This matches Abaad's
  `abaad_kdtree_build_f32(&tree, points, count, dim, ctx)`. Destructors
  are unchanged: they take one argument.
- **Descriptors:** settings structs where zero means the default, as now.
  They are knobs: they go after the buffers and sizes.
- **Names:** an output parameter is `out_*` or `*_out`; a bulk input is its
  plain name or `*_in`. The survey keys on these names and on `const`,
  as Alwan's does.
- **Shader cores** (`src/qaws/core/*_core.h`): `QAWS_OUT` marks outputs,
  because HLSL and GLSL drop `const`. Outputs move first there too, so the
  cores and their C callers read the same way. The cores are nearly
  compliant already (3 of 166 flagged).

## Scope (survey of the headers at this commit)

| Headers | Declarations | Compliant | R1 | R3 | R4 |
|---|---|---|---|---|---|
| Public, `src/qaws/*.h` (81 files) | 587 | 132 | 48 | 451 | 98 |
| Internal, `src/qaws/internal/*.h` | 130 | 58 | 17 | 65 | 21 |
| Shader cores, `src/qaws/core/*.h` (15) | 166 | 163 | 0 | 3 | 0 |

No parameter is named `*_stride` today, so R2 applies to no current
declaration. Nearly every violation is R3: an output after its inputs. The
heaviest headers are `qaws_inspect.h` (43), `qaws_repulsive.h` (44),
`qaws_surface_region.h` (37), `qaws_diff.h` (29 R3, 15 R1) and
`qaws_exact.h` (25). All the R1 violations are in the differentiation
headers, whose `qaws_diff_context` comes first.

Callers to update with each header: the tests (`tests/*.c`), the examples
and showcases (`examples/**/*.c`), the benchmarks, the HLSL/GLSL samples
that mirror the cores, the amalgamation, and every code snippet in
`docs/api/*.md`. About 218 C files call the API.

## How

1. **The survey, as a tool.** Qaws uses no Python, so the survey is a C
   program, `tools/qaws_api_survey.c`, built by CMake and Sharpmake. It
   parses the public, internal and core headers and checks R1 to R4
   using `const`, `QAWS_OUT` and the names. A ctest test fails on any
   violation not on an explicit, shrinking allow-list. It lands first,
   with the whole current list allowed; each migrated header then leaves
   the list. The build expects 100% at the end, as Alwan's does.
2. **One module at a time.** For each header group (repulsive, diff,
   inspect, exact, surfaces, regions, paths, 2D operations, evaluation and
   the small curve and surface families):
   - reorder its declarations and definitions;
   - update every caller in the same commit;
   - rerun that module's tests in both precisions (MSVC) and under GCC
     strict on Linux;
   - remove it from the allow-list.
   Each commit is `[refactor] parameter convention: <module>`.
   - The signatures change, so a missed call site fails to compile. Swaps
     between two pointers of the same type would not, so each commit also
     runs that module's full tests.
   - Swaps between same-typed pointers are where silent bugs would hide.
     Where a function has two same-typed buffers (`out, in` of
     `qaws_scalar*`), its tests must check the values, not just the status.
3. **Docs and version.** `CONTRIBUTING.md` gains the rules, citing Alwan's.
   The API docs' snippets follow each module. The change is breaking, so it
   carries a major version bump, noted in `CHANGELOG.md`. No compatibility
   layer: the old order would only survive in macros nobody can type-check.
4. **The order of the modules** (step 3 below):
   1. repulsive (the current branch's module, and the most recently
      designed);
   2. the differentiation headers (the only R1 violations);
   3. then by size: inspect, surface regions, exact, shell regions, paths;
   4. then the rest.

## Decisions

- **Both precisions in one binary.** As in Alwan, Suwar and Abaad, every
  scalar-bearing type and function exists as `_f32` and `_f64`, both in
  the same library. The build-time choice (`QAWS_SCALAR_IS_FLOAT`) no
  longer selects the library; it only selects which precision the
  unsuffixed names mean in a caller's translation unit.
- **A context only if needed.** Qaws has no global state today, and none
  of its allocators, executors or differentiation contexts need one, so
  no `qaws_ctx` is introduced unless the migration finds a need.
- **A dedicated branch**, `parameter_convention`, cut from `repulsive`.

## Both precisions: the design

Suwar includes an `_impl.inc` twice per source, under an f32 and an f64
setup. Qaws has about 85 sources written against one `qaws_scalar`, so it
does the same at the scale of the whole library instead:

- **Compiled twice.** Every precision-bearing source is compiled once with
  `QAWS_SCALAR_IS_FLOAT=1` and once with `=0`, into one library (two CMake
  object libraries, one Sharpmake configuration each). The few sources
  that do not depend on the scalar type are compiled once: the status
  strings, the integer clipper and offsetter (`qaws_clip64.c`,
  `qaws_offset64.c`) and their headers. The library has no global data
  (`nm` finds none), so nothing splits into two copies.
- **Renamed.** A generated header, `qaws_names.h`, maps every
  precision-bearing name to its suffixed form for the current precision:
  `#define qaws_curve_eval qaws_curve_eval_f32`, and the types alike
  (`qaws_vec2` to `qaws_vec2_f32`, `qaws_curve` to `qaws_curve_f32`).
  - It is re-includable: including it again under the other precision
    flips every macro.
  - About 780 external symbols (124 of them internal) and the scalar types
    and opaque handles are mapped.
  - A missed symbol cannot go unnoticed: the two object sets then define
    it twice, and the link fails.
- **Headers in two parts.** Each public header keeps a once-only part for
  what does not depend on the precision (enums, constants), and a part
  guarded per precision (`QAWS_<NAME>_F32`, `QAWS_<NAME>_F64`) for the
  rest. A caller can therefore include it under both precisions.
- **For callers.** Including a header gives the explicitly suffixed names
  of its precision, plus the unsuffixed names as aliases for the
  precision the caller chose (`QAWS_SCALAR_IS_FLOAT`, default 1, as now).
  `qaws_both.h` includes everything under both precisions, for a
  translation unit that wants `_f32` and `_f64` side by side. Existing
  callers therefore keep compiling with only the convention's reordering.
- **Shader cores** (`core/*_core.h`) are static inline, compiled by
  whoever includes them, and need no renaming.
- **The scalar.** `qaws_f32` and `qaws_f64` are the scalar types;
  `qaws_scalar` is the alias of the current precision.

The tool generating `qaws_names.h` is the same C program as the survey
(`tools/qaws_api_tool.c`): it reads the headers' declarations, so the map
and the convention check agree on what the API is.

As built (steps 1 and 2):

- **The wrappers** are generated by the build systems, not committed:
  - CMake writes `qaws_precision/<name>_f32.c` and `<name>_f64.c` into the
    build tree at configure time, and `QawsLib.cs` does the same into
    `tmp/qaws_precision/`;
  - each wrapper is three lines: `#undef QAWS_SCALAR_IS_FLOAT`, its value,
    `#include` of the source;
  - the sources themselves are never compiled directly. With one source
    compiled directly and also included by its f32 wrapper, MSBuild did
    not rebuild the f32 object when the source changed, so stale objects
    were linked. With two wrappers, touching a source rebuilds exactly its
    two objects.
- **The status strings** (`qaws_status.c`) are defined by the f64 pass
  only, which compiles them once in the library and in an amalgamation
  compiled twice.
- **`qaws_status.h` includes `qaws_platform.h`.** Most headers include the
  status first, so the names are active before anything is declared.
- **MSVC:** `/wd4305 /wd4244` for the library, since its f32 objects
  narrow double literals; f32 builds were not `/WX` before.
- **The amalgamation** is compiled through the same two wrappers. Two
  fixes were needed:
  - its generator now inlines include lines that end with a comment;
  - four static names collided between sources (`sr_job`, `rm_cross`,
    `sm_d3`, `SM_CORNER`; already broken on `repulsive`).
  It builds with `-Werror`, holds both precisions, and passes the suites
  run against it.
- **Checks:** two tests, `qaws_api_convention` (the survey with
  `tools/qaws_api_allow.txt`, 744 entries at the start) and
  `qaws_api_names` (`qaws_names.h` current). Test 121 uses both precisions
  in one translation unit, from an f32 and an f64 caller.
- **Results:**
  - the library holds 779 functions in each precision, plus
    `qaws_status_to_string`;
  - GCC builds it with `-Werror` and MSVC with `/W4 /WX`;
  - the full test suite passes against it on Linux in both callers'
    precisions, except `108_shell_trimmed` in f64, which fails the same
    way on `repulsive` before this work;
  - the repulsive suites (112 to 121) pass on MSVC in both.

## The order of the work

1. The tool (survey and names), with the current violations allowed.
2. Both precisions: the build, `qaws_names.h`, the headers' two parts. All
   tests then run against one library, in each precision.
3. The convention, one module at a time, as above.
4. `CONTRIBUTING.md`, the docs' snippets, the version and the changelog.
5. Back to `repulsive`: the Repulsive Surfaces figures.

Steps 1 to 4 are done (3.0.0). How step 3 went:

- **The order of each declaration.** `qaws_api_tool survey -order` sorts a
  declaration's parameters by rank (stably) and prints the result as old
  indices. A perl script, not part of the tree, applied it to every call,
  declaration and definition. It parses the arguments with balanced
  brackets and keeps the separators. It also covers the `_f32` and `_f64`
  spellings and the docs, C++ and shader files. Calls inside comments or
  strings with other arguments were left alone and fixed by hand where
  they described an order.
- **The classifier's rules**, refined module by module:
  - a struct is data when a public header gives it a body, and an opaque
    handle otherwise;
  - a pointer is an input only when it is const at every level (`T const**`
    hands back a pointer, so it is an output);
  - function pointers are callbacks;
  - an array bound may hold a `*`;
  - an array of descriptors is a buffer, and its count follows it;
  - the allocator and executor go last.
- **The commits:** repulsive, differentiation, curves, surfaces, 2D
  operations, regions, exact, internal headers, shader cores. Each family
  passed the whole f64 suite on Linux (108 fails as on `repulsive`), and
  every target built with MSVC. The cores' HLSL kernels compile with
  `dxc` (cs_6_0), as before; three entries already failed with `fxc`
  before this work.
- **The result:** 851 of 851 declarations in the convention, and the
  allow-list is empty.

## Merged into meshing

`meshing` took the convention's order first (its own API included: the
meshing headers, views, allocators, arenas, the repulsive additions), then
both precisions as designed above. What differs from `parameter_convention`:

- **The allocator exists once.** `meshing` keeps process-wide state: the
  default allocator (atomics), the thread and scratch allocators
  (thread-local), and the arena, all in `internal/qaws_internal_memory.c`.
  Compiled per precision, it would hold two defaults and two thread
  allocators, and an allocator installed for one precision would not serve
  the other. So:
  - `qaws_allocator.h` (new, once-only, precision-free) holds the
    descriptor, `QAWS_MEMORY_ALIGNMENT`, `qaws_allocator_extended`, the
    default, thread and scratch settings and the arena; `qaws_types.h` and
    `qaws_alloc.h` include it, and `qaws_alloc.h` keeps only the `_ex`
    creators, per precision;
  - `internal/qaws_internal_memory.h` is once-only too, and now also holds
    `qaws_internal_alloc`, `_dealloc`, `_calloc`, `_scratch` and
    `_scratch_release` (moved from `qaws_internal_curve.c`: they only
    forward to the memory module);
  - `qaws_internal_memory.c` is compiled by the f64 pass only, like
    `qaws_status.c` (skipped for f32 by `CMakeLists.txt` and `QawsLib.cs`,
    `#if !QAWS_SCALAR_IS_FLOAT` for the twice-compiled amalgamation).
  The library exports each of these names once, unsuffixed; test 155
  installs one counting thread allocator and sees an f32 and an f64 curve
  both allocate from it and give every block back.
- **No other process-wide state.** Every other file-scope variable is a
  const table. `rv_trace` (the CVT remesher's test probe) exists once per
  precision: the tests set the one of their own precision.
- **The mesh module's internals** (`src/qaws/mesh/*.h`) name their
  functions with their own prefixes (`hm_`, `cf_`, `mf_`, ...), and
  `vm_corners`, `vm_facets` and `rv_trace` are extern variables. The names
  tool maps them too when given after `-any`; both CTest checks read the
  mesh headers with the internal ones.
- **The test** is `155_precision_both` (121 is `121_mesh` here).
- **The clipper, offsetter and exact kernel** are compiled twice as on
  `parameter_convention`: they hold no state, so the copies are only code.
