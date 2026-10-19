# API Reference: qaws_alloc.h

Custom allocator support for curve and surface creation.

---

## Overview

Every allocation of the library goes through one place (`src/qaws/internal/qaws_internal_memory.c`, the only caller of `malloc`, `realloc` and `free`). Memory comes from, in order:

1. the allocator passed to an `_ex` function;
2. the calling thread's allocator, `qaws_set_thread_allocator`;
3. the process default, `qaws_set_default_allocator`;
4. `malloc`.

This covers objects, results, and all internal scratch: operations, sampling, traversal, inspection, clipping, meshing, batches and repulsive solvers.

- **Ownership.** Each block carries a 16-byte header naming the allocator that gave it. Freeing and growing always go back to that allocator, whatever the settings at the time and whichever thread releases it.
- **Allocator requirements.** An allocator must return 16-byte aligned memory and outlive its blocks.
- **Setting the default.** Set the default before calls run. Set a thread's allocator on that thread.

```c
void qaws_set_default_allocator(qaws_allocator const* allocator);   /* NULL: malloc */
void qaws_set_thread_allocator(qaws_allocator const* allocator);    /* this thread; NULL: the default */
```

## Caller memory: qaws_arena

```c
typedef struct qaws_arena {
    qaws_allocator allocator;   /* pass &arena.allocator */
    unsigned char* base;
    size_t size, used, peak, last;
} qaws_arena;

void qaws_arena_init(qaws_arena* arena, void* memory, size_t size);
void qaws_arena_reset(qaws_arena* arena);
```

An allocator over memory the host provides: a mapped file, pinned or GPU-visible memory, or a frame arena.

- **Allocation.** Blocks are cut in order and 16-byte aligned.
- **Freeing.** A free releases only the latest block, so a growing array reuses its space. Everything else is released at `qaws_arena_reset`. `peak` tells how much a workload needed.
- **When full.** Allocation fails, and calls return `QAWS_STATUS_ALLOCATION_FAILURE`. Nothing falls back to the heap.
- **Threads.** One thread at a time. For workers, use one arena per thread with `qaws_set_thread_allocator(&arena.allocator)`.

Test 154 checks this:

- a counting default allocator sees and gets back every block of a mixed workload;
- the thread allocator takes precedence over the default;
- blocks return to their own allocator;
- the workload runs inside a 2 MB caller buffer (peak 386 KB);
- a 2 KB buffer fails cleanly.

## _ex creation functions

The `_ex` variants of creation functions behave identically to their non-`_ex` counterparts but take their memory from the provided `qaws_allocator`. Pass NULL for `allocator` to use the order above.

Curves and surfaces created with a custom allocator MUST be destroyed with `qaws_curve_destroy` or `qaws_surface_destroy` (which read the stored allocator). The allocator must remain valid until the object is destroyed.

---

## qaws_curve_create_bezier_ex

```c
qaws_status qaws_curve_create_bezier_ex(
    qaws_bezier_desc const* desc,
    qaws_allocator const* allocator,
    qaws_curve** out_curve);
```

Allocator-aware variant of `qaws_curve_create_bezier`.

**Parameters:**
- `desc` -- Bezier curve descriptor.
- `allocator` -- Custom allocator, or NULL for default.
- `out_curve` -- Receives the created curve.

**Returns:** `QAWS_STATUS_OK` on success.

---

## qaws_curve_create_hermite_ex

```c
qaws_status qaws_curve_create_hermite_ex(
    qaws_hermite_desc const* desc,
    qaws_allocator const* allocator,
    qaws_curve** out_curve);
```

Allocator-aware variant of `qaws_curve_create_hermite`.

**Parameters:**
- `desc` -- Hermite curve descriptor.
- `allocator` -- Custom allocator, or NULL for default.
- `out_curve` -- Receives the created curve.

**Returns:** `QAWS_STATUS_OK` on success.

---

## qaws_curve_create_catmull_rom_ex

```c
qaws_status qaws_curve_create_catmull_rom_ex(
    qaws_catmull_rom_desc const* desc,
    qaws_allocator const* allocator,
    qaws_curve** out_curve);
```

Allocator-aware variant of `qaws_curve_create_catmull_rom`.

**Parameters:**
- `desc` -- Catmull-Rom curve descriptor.
- `allocator` -- Custom allocator, or NULL for default.
- `out_curve` -- Receives the created curve.

**Returns:** `QAWS_STATUS_OK` on success.

---

## qaws_curve_create_bspline_ex

```c
qaws_status qaws_curve_create_bspline_ex(
    qaws_bspline_desc const* desc,
    qaws_allocator const* allocator,
    qaws_curve** out_curve);
```

Allocator-aware variant of `qaws_curve_create_bspline`.

**Parameters:**
- `desc` -- B-Spline curve descriptor.
- `allocator` -- Custom allocator, or NULL for default.
- `out_curve` -- Receives the created curve.

**Returns:** `QAWS_STATUS_OK` on success.

---

## qaws_curve_create_nurbs_ex

```c
qaws_status qaws_curve_create_nurbs_ex(
    qaws_nurbs_desc const* desc,
    qaws_allocator const* allocator,
    qaws_curve** out_curve);
```

Allocator-aware variant of `qaws_curve_create_nurbs`.

**Parameters:**
- `desc` -- NURBS curve descriptor.
- `allocator` -- Custom allocator, or NULL for default.
- `out_curve` -- Receives the created curve.

**Returns:** `QAWS_STATUS_OK` on success.

---

## qaws_curve_create_trajectory_ex

```c
qaws_status qaws_curve_create_trajectory_ex(
    qaws_trajectory_desc const* desc,
    qaws_allocator const* allocator,
    qaws_curve** out_curve);
```

Allocator-aware variant of `qaws_curve_create_trajectory`.

**Parameters:**
- `desc` -- Trajectory curve descriptor.
- `allocator` -- Custom allocator, or NULL for default.
- `out_curve` -- Receives the created curve.

**Returns:** `QAWS_STATUS_OK` on success.

---

## The other curve families

```c
qaws_status qaws_curve_create_rational_bezier_ex(qaws_rational_bezier_desc const* desc, qaws_allocator const* allocator, qaws_curve** out_curve);
qaws_status qaws_curve_create_yuksel_ex(qaws_yuksel_desc const* desc, qaws_allocator const* allocator, qaws_curve** out_curve);
qaws_status qaws_curve_create_arc_ex(qaws_arc_desc const* desc, qaws_allocator const* allocator, qaws_curve** out_curve);
qaws_status qaws_curve_create_polynomial_ex(qaws_polynomial_desc const* desc, qaws_allocator const* allocator, qaws_curve** out_curve);
qaws_status qaws_curve_create_clothoid_ex(qaws_clothoid_desc const* desc, qaws_allocator const* allocator, qaws_curve** out_curve);
qaws_status qaws_curve_create_composite_ex(qaws_composite_desc const* desc, qaws_allocator const* allocator, qaws_curve** out_curve);
qaws_status qaws_curve_create_subdivision_ex(qaws_subdivision_desc const* desc, qaws_allocator const* allocator, qaws_curve** out_curve);
```

Each is the allocator-aware variant of the constructor without `_ex`,
which calls it with NULL. Every allocation the constructor makes goes
through the allocator: the curve, its stored data, and its scratch (the
subdivision's refinement levels). A composite owns its segments only on
success; the segments keep the allocators they were made with.

**Checked by test 146.** Every curve and surface constructor here runs with
an allocator that fails at its k-th allocation, for every k. Each failure
returns `QAWS_STATUS_ALLOCATION_FAILURE` with no object and nothing left
allocated. A success, once destroyed, leaves nothing either.

---

## qaws_surface_create_bezier_ex

```c
qaws_status qaws_surface_create_bezier_ex(
    qaws_surface_bezier_desc const* desc,
    qaws_allocator const* allocator,
    qaws_surface** out_surface);
```

Allocator-aware variant of `qaws_surface_create_bezier`.

**Parameters:**
- `desc` -- Bezier surface descriptor.
- `allocator` -- Custom allocator, or NULL for default.
- `out_surface` -- Receives the created surface.

**Returns:** `QAWS_STATUS_OK` on success.

---

## qaws_surface_create_bspline_ex

```c
qaws_status qaws_surface_create_bspline_ex(
    qaws_surface_bspline_desc const* desc,
    qaws_allocator const* allocator,
    qaws_surface** out_surface);
```

Allocator-aware variant of `qaws_surface_create_bspline`.

**Parameters:**
- `desc` -- B-Spline surface descriptor.
- `allocator` -- Custom allocator, or NULL for default.
- `out_surface` -- Receives the created surface.

**Returns:** `QAWS_STATUS_OK` on success.

---

## qaws_surface_create_nurbs_ex

```c
qaws_status qaws_surface_create_nurbs_ex(
    qaws_surface_nurbs_desc const* desc,
    qaws_allocator const* allocator,
    qaws_surface** out_surface);
```

Allocator-aware variant of `qaws_surface_create_nurbs`.

**Parameters:**
- `desc` -- NURBS surface descriptor.
- `allocator` -- Custom allocator, or NULL for default.
- `out_surface` -- Receives the created surface.

**Returns:** `QAWS_STATUS_OK` on success.

---

## The other surface families

```c
qaws_status qaws_surface_create_bilinear_ex(qaws_surface_bilinear_desc const* desc, qaws_allocator const* allocator, qaws_surface** out_surface);
qaws_status qaws_surface_create_biquadratic_ex(qaws_surface_biquadratic_desc const* desc, qaws_allocator const* allocator, qaws_surface** out_surface);
qaws_status qaws_surface_create_coons_ex(qaws_surface_coons_desc const* desc, qaws_allocator const* allocator, qaws_surface** out_surface);
qaws_status qaws_surface_create_extrusion_ex(qaws_surface_extrusion_desc const* desc, qaws_allocator const* allocator, qaws_surface** out_surface);
qaws_status qaws_surface_create_offset_ex(qaws_surface_offset_desc const* desc, qaws_allocator const* allocator, qaws_surface** out_surface);
qaws_status qaws_surface_create_pipe_ex(qaws_surface_pipe_desc const* desc, qaws_allocator const* allocator, qaws_surface** out_surface);
qaws_status qaws_surface_create_revolution_ex(qaws_surface_revolution_desc const* desc, qaws_allocator const* allocator, qaws_surface** out_surface);
qaws_status qaws_surface_create_ruled_ex(qaws_surface_ruled_desc const* desc, qaws_allocator const* allocator, qaws_surface** out_surface);
qaws_status qaws_surface_create_swept_ex(qaws_surface_swept_desc const* desc, qaws_allocator const* allocator, qaws_surface** out_surface);
qaws_status qaws_surface_create_loft_ex(qaws_surface_loft_desc const* desc, qaws_allocator const* allocator, qaws_surface** out_surface);
qaws_status qaws_surface_create_gordon_ex(qaws_surface_gordon_desc const* desc, qaws_allocator const* allocator, qaws_surface** out_surface);
qaws_status qaws_surface_create_trimmed_ex(qaws_surface_trim_desc const* desc, qaws_allocator const* allocator, qaws_surface** out_surface);
qaws_status qaws_surface_create_tspline_ex(qaws_surface_tspline_desc const* desc, qaws_allocator const* allocator, qaws_surface** out_surface);
qaws_status qaws_surface_create_subdiv_ex(qaws_surface_subdiv_desc const* desc, qaws_allocator const* allocator, qaws_surface** out_surface);
qaws_status qaws_surface_create_fillet_ex(qaws_surface_fillet_desc const* desc, qaws_allocator const* allocator, qaws_surface** out_surface);
qaws_status qaws_surface_create_fillet_with_trims_ex(qaws_surface_fillet_desc const* desc, qaws_allocator const* allocator,
    qaws_surface** out_fillet, qaws_curve** out_trim_a, qaws_curve** out_trim_b);
```

As for curves, the plain constructors call these with NULL, and every
allocation the constructor makes goes through the allocator. That covers
the surface, its stored data, and construction scratch, such as subdivision
meshes and the fillet's intersection points. The fillet's trim curves come
from it too. Two things are not counted: the intersection search the fillet
runs, which allocates its own scratch, and the curves a surface borrows.

**Evaluation.** The loft and Gordon evaluators need scratch per section
curve. Up to 32 sections it is on the stack; beyond that it comes from the
surface's own allocator, so evaluation never calls malloc. A host that
evaluates a large loft or Gordon surface from several threads needs an
allocator that is safe across threads. Trimmed surfaces test containment
without allocating.

Test 146 covers all of these with an allocator failing at each allocation.

---

## Derived curves

An operation that makes curves from a source curve gives its results the
source curve's allocator: reverse, extract, split, join, conversion, degree
elevation and reduction, merging a chain (the first curve's allocator),
transforms, and the arc-length and length-matched wrappers. With a host
allocator on the source, everything derived from it lives in that
allocator, and is destroyed with `qaws_curve_destroy` as usual.

Two kinds of curve keep the default allocator for now: curves built from
raw arrays (polyline and offset outputs, clipping, skeletons, stacks blending
several curves), and the scratch an operation frees before it returns.
Test 146 checks reverse, extract, split, degree elevation, B-spline to
NURBS, and both wrappers. It injects a failure at each allocation and
checks that the results are counted in the source's allocator.
