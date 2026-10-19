# API Reference: qaws_alloc.h

Where the library's memory comes from, and how hosts give it their own.

---

## Overview

Every allocation of the library goes through one place: `src/qaws/internal/qaws_internal_memory.c` is the only caller of the C allocator. Memory comes from, in order:

1. the allocator passed to an `_ex` function;
2. the calling thread's allocator, `qaws_set_thread_allocator`;
3. the process default, `qaws_set_default_allocator`;
4. the C allocator, aligned to `QAWS_MEMORY_ALIGNMENT` (`_aligned_malloc` on MSVC).

This covers objects, results and all internal scratch: operations, sampling, traversal, inspection, clipping, meshing, batches and repulsive solvers.

## The descriptor

```c
typedef struct qaws_allocator {
    void* (*alloc)(unsigned long size, void* user_data);              /* legacy, or qaws_allocator_extended */
    void  (*dealloc)(void* ptr, void* user_data);                      /* legacy */
    void* user_data;
    void* (*alloc_sized)(size_t size, size_t alignment, void* user_data);
    void  (*dealloc_sized)(void* ptr, size_t size, size_t alignment, void* user_data);
    int   (*resize)(void* ptr, size_t old_size, size_t new_size, size_t alignment, void* user_data);
} qaws_allocator;

#define QAWS_MEMORY_ALIGNMENT 16
void* qaws_allocator_extended(unsigned long size, void* user_data);   /* the marker */
```

The library recognizes two kinds of descriptor:

- **Legacy.** Existing code sets `alloc`, `dealloc` and `user_data`. The other members are never read, so descriptors filled field by field on the stack keep working, and so do binaries built against the old three-member layout. Sizes are `unsigned long`, which is 32 bits on Windows. A larger request fails without calling `alloc`; it is never truncated.
- **Extended.** Set `alloc = qaws_allocator_extended` (a marker that never allocates) and the sized members. Sizes are `size_t`, with the alignment passed. `dealloc_sized` receives the size and alignment it was allocated with. `resize`, if set, grows or shrinks a block in place and returns non-zero, or returns 0 to let the library allocate, copy and free.

Both kinds have the same alignment rule:

- Blocks must be aligned to `QAWS_MEMORY_ALIGNMENT` (16). That is enough for every type the library keeps on the heap; its aligned SIMD loads only read stack arrays.
- A misaligned block is given back to the allocator and the call fails with `QAWS_STATUS_ALLOCATION_FAILURE`. Debug builds also assert.

## Lifetimes

- **The descriptor** can be short-lived. Settings and `_ex` arguments are copied: callbacks and `user_data`. Each block keeps the copy it was made with, and so does every curve, surface or prepared object that allocates again later (derived curves, evaluation scratch). A descriptor can therefore be a local that goes out of scope before the objects it made (test 154).
- **The state `user_data` points to** must outlive every block of that allocator. Copying the callbacks does not extend it.
- **Blocks** are freed and grown by the allocator that gave them, whatever the settings at that time and whichever thread releases them (test 154: four workers' blocks freed on another thread under another allocator).
- **Results** live as long as their objects and do not depend on any call's scratch. Scratch is meant to be freed before a call returns. Test 154 checks this for a mixed workload, both completed and failing for lack of memory inside a 2 KB arena (no block left held). It is not checked for every failure path of every operation.

## Settings

```c
qaws_status qaws_set_default_allocator(qaws_allocator const* allocator);   /* NULL: the C allocator */
void qaws_set_thread_allocator(qaws_allocator const* allocator);            /* this thread; NULL: the default */
void qaws_get_thread_allocator(qaws_allocator* out_allocator);              /* zeroed: none */
void qaws_exchange_thread_allocator(qaws_allocator const* allocator, qaws_allocator* out_previous);
```

- **The default is write-once in effect.** It can be set or cleared until the library first allocates from it. After that, `qaws_set_default_allocator` returns `QAWS_STATUS_UNSUPPORTED_OPERATION` and changes nothing. Setting and sealing are synchronized with atomics. Once sealed, reading the default is a plain load, so allocations take no lock.
- **A thread's allocator** applies to that thread only, and the library never changes it. To nest, install with `qaws_exchange_thread_allocator(&mine, &previous)` and restore with `qaws_set_thread_allocator(&previous)` on every return path: success, failure and cancellation.

```c
qaws_allocator previous;
qaws_exchange_thread_allocator(&frame_arena.allocator, &previous);
status = qaws_curve_offset_2d(curve, d, 1, out, cap, &count);   /* every return path below restores */
qaws_set_thread_allocator(&previous);
```

## Caller memory: qaws_arena

```c
typedef struct qaws_arena {
    qaws_allocator allocator;      /* pass &arena.allocator */
    unsigned char* base;
    size_t size;
    size_t used;                   /* up to the end of the latest block */
    size_t peak;                   /* most used since init */
    size_t last;                   /* the latest block's offset */
    size_t live;                   /* blocks held */
    size_t failures;               /* requests refused for room */
    size_t invalid_frees;          /* frees of blocks it does not hold */
    unsigned int generation;       /* resets so far */
} qaws_arena;

void qaws_arena_init(qaws_arena* arena, void* memory, size_t size);
qaws_status qaws_arena_reset(qaws_arena* arena);
```

An extended allocator over memory the host provides: a mapped file, pinned memory, or a frame arena.

- **CPU access.** The memory must be CPU-accessible and coherent for the CPU. Putting the library's data in GPU-visible memory does not run any of its algorithms on a GPU.
- **Allocation.** Blocks are cut in order and aligned. The latest block grows and shrinks in place, which is how growing arrays stay compact. A freed latest block is released at once; other freed blocks are released at reset.
- **Hard bound.** When full, allocation fails and the call returns `QAWS_STATUS_ALLOCATION_FAILURE`. Nothing falls back to the heap. All offsets are computed without overflow, including requests near `SIZE_MAX`.
- **Reset.** `qaws_arena_reset` refuses (`QAWS_STATUS_INVALID_ARGUMENT`) while blocks are held (`live > 0`). Destroy every arena-backed object and finish every call first. Freeing a block the arena does not hold (another arena, an earlier generation, or a double free) is counted in `invalid_frees`, and debug builds assert.
- **Threads.** One thread at a time. For workers, use one arena per thread with `qaws_exchange_thread_allocator`.

Arena-backed objects never move to the heap. To keep an object beyond a reset, create it with another allocator (the `_ex` argument) while scratch uses the arena (the thread setting).

## _ex creation functions

The `_ex` variants of creation functions behave identically to their non-`_ex` counterparts but take their memory from the provided `qaws_allocator`. Pass NULL for `allocator` to use the order above.

Curves and surfaces created with a custom allocator are destroyed with `qaws_curve_destroy` or `qaws_surface_destroy` as usual; they keep a copy of the descriptor, so only the state behind `user_data` must remain valid until then.

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
