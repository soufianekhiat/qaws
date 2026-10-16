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
4. **The order of the modules:**
   1. repulsive (the current branch's module, and the most recently
      designed);
   2. the differentiation headers (the only R1 violations);
   3. then by size: inspect, surface regions, exact, shell regions, paths;
   4. then the rest.

## Open questions

- **Precision forms.** Alwan, Suwar and Abaad expose `_f32` and `_f64`
  forms of every numerical operation. Qaws picks one precision at build
  time (`QAWS_SCALAR_IS_FLOAT`). Moving to both forms is a separate,
  larger change (every type and entry point twice), so it is not part of
  this plan unless decided otherwise.
- **A context object.** Qaws has no general `qaws_ctx`: allocators,
  executors and differentiation contexts are passed in descriptors or
  separately. Alwan, Suwar and Abaad carry an allocator and thread pool in
  their contexts. R1 applies to whatever context exists; introducing a
  `qaws_ctx` would be its own decision.
- **The branch.** Commits go to the current branch (`repulsive`); a
  separate branch for the migration is the user's call.
