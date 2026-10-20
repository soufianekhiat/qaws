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
void qaws_exchange_thread_allocator(qaws_allocator* out_previous, qaws_allocator const* allocator);
```

- **The default is write-once in effect.** It can be set or cleared until the library first allocates from it. After that, `qaws_set_default_allocator` returns `QAWS_STATUS_UNSUPPORTED_OPERATION` and changes nothing. Setting and sealing are synchronized with atomics. Once sealed, reading the default is a plain load, so allocations take no lock.
- **A thread's allocator** applies to that thread only, and the library never changes it. To nest, install with `qaws_exchange_thread_allocator(&mine, &previous)` and restore with `qaws_set_thread_allocator(&previous)` on every return path: success, failure and cancellation.

```c
qaws_allocator previous;
qaws_exchange_thread_allocator(&previous, &frame_arena.allocator);
status = qaws_curve_offset_2d(out, &count, curve, cap, d, 1);   /* every return path below restores */
qaws_set_thread_allocator(&previous);
```

- **Setting the default against the first allocation.** Sealing and setting go through one compare-and-swap state, and the descriptor is copied before the state that publishes it is stored (sequentially consistent atomics, `Interlocked*` on MSVC). An allocation therefore sees either no default or the whole copy, never a partly written one. Test 154 starts 24 fresh processes in which seven workers allocate while thread 0 sets the default: at once, while a worker enters its first allocation, or after it. Each run ends with the default installed or refused as sealed, the allocator's callbacks always see their own `user_data`, and every block returns to the allocator that gave it.

## Scratch apart from results

```c
void qaws_exchange_thread_scratch_allocator(qaws_allocator* out_previous, qaws_allocator const* allocator);   /* NULL: none */
```

With a scratch allocator installed on a thread, the temporaries of adaptive sampling, tessellation, meshing, clipping and the repulsive solvers come from it: blocks allocated and freed within one call. Results and objects still follow the order above (`_ex` argument, thread allocator, default). A frame arena can therefore serve as scratch and be reset after every call while the results stay on the heap.

```c
qaws_allocator previous;
qaws_exchange_thread_scratch_allocator(&previous, &frame_arena.allocator);
status = qaws_mesh_generate_surface(&mesh, NULL, surface, &desc);   /* the mesh: heap; its temporaries: the arena */
qaws_exchange_thread_scratch_allocator(NULL, &previous);
qaws_arena_reset(&frame_arena);                                    /* OK: the arena holds nothing after the call */
```

- **Which allocations count as scratch.** A local pointer allocated and freed in the same function and never stored elsewhere (386 sites in those modules). Everything else, including every block a result keeps, uses the persistent order, so a result never holds scratch memory.
- **Coverage.** Other operations take their temporaries from the persistent order as before. Without a scratch allocator, nothing changes.
- **Tests.** Test 154 runs sampling, offsets, clipping, tessellation, meshing and a repulsive solve with an arena as scratch and a counting heap as the allocator. After each call the arena holds no block and resets, and the kept results (curves, a clip result, a mesh, a solver) are read after the resets and then returned to the heap. With `QAWS_TEST_SCRATCH=1`, the test runner installs a scratch arena for every suite and fails a suite that leaves a block in it.

## Sizes and overflow

- **Allocation sizes** are computed with checked arithmetic (`qaws_internal_mul2/3/4`, `qaws_internal_add2`). A product or sum that does not fit `size_t` becomes `SIZE_MAX`, which every allocation refuses before calling an allocator. The call returns `QAWS_STATUS_ALLOCATION_FAILURE`: the request cannot be represented, so it is reported the same way as memory the allocator does not have. A size never wraps into a short block that the library would then overrun.
- **Element counts** that the library indexes with `unsigned int` (control points `u × v` of a surface, indices, list capacities) are refused at the public entry with the function's count status (`QAWS_STATUS_INVALID_CONTROL_POINT_COUNT` and similar) when their product does not fit. Growing lists double up to `UINT_MAX` and fail beyond it.
- **Order.** Invalid dimensions and counts are refused before any allocation. Allocation failures happen before any result is published: a failed call leaves its output pointers `NULL` and frees what it allocated (test 154, every part in arenas from 256 bytes up, 0 blocks held and 0 results published on failure).
- **The legacy interface** counts in `unsigned long`. A request above `ULONG_MAX` (4 GB on Windows) fails without calling `alloc`.

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
- **Resize.** `qaws_internal_realloc` grows or shrinks the latest block in place. Any other block is copied into a new block, which the arena records as the latest, and the old one is marked freed. If a reallocation fails, the old block stays held and unchanged (test 154).
- **Threads.** One thread at a time: allocation, free, reset, and the destruction of arena-backed objects. Destroying an arena-backed object on another thread while the owner resets is a race, and that use is unsupported. Reset refuses while any block is held, so a reset that does not race cannot invalidate a live object or result. For workers, use one arena per thread with `qaws_exchange_thread_allocator` (test 154 runs four workers this way, with arenas that are large enough and arenas that fail).

Arena-backed objects never move to the heap. To keep an object beyond a reset, create it with another allocator (the `_ex` argument) while scratch uses the arena (the thread setting).

## _ex creation functions

The `_ex` variants of creation functions behave identically to their non-`_ex` counterparts but take their memory from the provided `qaws_allocator`. Pass NULL for `allocator` to use the order above.

Curves and surfaces created with a custom allocator are destroyed with `qaws_curve_destroy` or `qaws_surface_destroy` as usual; they keep a copy of the descriptor, so only the state behind `user_data` must remain valid until then.

---

## qaws_curve_create_bezier_ex

```c
qaws_status qaws_curve_create_bezier_ex(
    qaws_curve** out_curve,
    qaws_bezier_desc const* desc,
    qaws_allocator const* allocator);
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
    qaws_curve** out_curve,
    qaws_hermite_desc const* desc,
    qaws_allocator const* allocator);
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
    qaws_curve** out_curve,
    qaws_catmull_rom_desc const* desc,
    qaws_allocator const* allocator);
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
    qaws_curve** out_curve,
    qaws_bspline_desc const* desc,
    qaws_allocator const* allocator);
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
    qaws_curve** out_curve,
    qaws_nurbs_desc const* desc,
    qaws_allocator const* allocator);
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
    qaws_curve** out_curve,
    qaws_trajectory_desc const* desc,
    qaws_allocator const* allocator);
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
qaws_status qaws_curve_create_rational_bezier_ex(qaws_curve** out_curve, qaws_rational_bezier_desc const* desc, qaws_allocator const* allocator);
qaws_status qaws_curve_create_yuksel_ex(qaws_curve** out_curve, qaws_yuksel_desc const* desc, qaws_allocator const* allocator);
qaws_status qaws_curve_create_arc_ex(qaws_curve** out_curve, qaws_arc_desc const* desc, qaws_allocator const* allocator);
qaws_status qaws_curve_create_polynomial_ex(qaws_curve** out_curve, qaws_polynomial_desc const* desc, qaws_allocator const* allocator);
qaws_status qaws_curve_create_clothoid_ex(qaws_curve** out_curve, qaws_clothoid_desc const* desc, qaws_allocator const* allocator);
qaws_status qaws_curve_create_composite_ex(qaws_curve** out_curve, qaws_composite_desc const* desc, qaws_allocator const* allocator);
qaws_status qaws_curve_create_subdivision_ex(qaws_curve** out_curve, qaws_subdivision_desc const* desc, qaws_allocator const* allocator);
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
    qaws_surface** out_surface,
    qaws_surface_bezier_desc const* desc,
    qaws_allocator const* allocator);
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
    qaws_surface** out_surface,
    qaws_surface_bspline_desc const* desc,
    qaws_allocator const* allocator);
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
    qaws_surface** out_surface,
    qaws_surface_nurbs_desc const* desc,
    qaws_allocator const* allocator);
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
qaws_status qaws_surface_create_bilinear_ex(qaws_surface** out_surface, qaws_surface_bilinear_desc const* desc, qaws_allocator const* allocator);
qaws_status qaws_surface_create_biquadratic_ex(qaws_surface** out_surface, qaws_surface_biquadratic_desc const* desc, qaws_allocator const* allocator);
qaws_status qaws_surface_create_coons_ex(qaws_surface** out_surface, qaws_surface_coons_desc const* desc, qaws_allocator const* allocator);
qaws_status qaws_surface_create_extrusion_ex(qaws_surface** out_surface, qaws_surface_extrusion_desc const* desc, qaws_allocator const* allocator);
qaws_status qaws_surface_create_offset_ex(qaws_surface** out_surface, qaws_surface_offset_desc const* desc, qaws_allocator const* allocator);
qaws_status qaws_surface_create_pipe_ex(qaws_surface** out_surface, qaws_surface_pipe_desc const* desc, qaws_allocator const* allocator);
qaws_status qaws_surface_create_revolution_ex(qaws_surface** out_surface, qaws_surface_revolution_desc const* desc, qaws_allocator const* allocator);
qaws_status qaws_surface_create_ruled_ex(qaws_surface** out_surface, qaws_surface_ruled_desc const* desc, qaws_allocator const* allocator);
qaws_status qaws_surface_create_swept_ex(qaws_surface** out_surface, qaws_surface_swept_desc const* desc, qaws_allocator const* allocator);
qaws_status qaws_surface_create_loft_ex(qaws_surface** out_surface, qaws_surface_loft_desc const* desc, qaws_allocator const* allocator);
qaws_status qaws_surface_create_gordon_ex(qaws_surface** out_surface, qaws_surface_gordon_desc const* desc, qaws_allocator const* allocator);
qaws_status qaws_surface_create_trimmed_ex(qaws_surface** out_surface, qaws_surface_trim_desc const* desc, qaws_allocator const* allocator);
qaws_status qaws_surface_create_tspline_ex(qaws_surface** out_surface, qaws_surface_tspline_desc const* desc, qaws_allocator const* allocator);
qaws_status qaws_surface_create_subdiv_ex(qaws_surface** out_surface, qaws_surface_subdiv_desc const* desc, qaws_allocator const* allocator);
qaws_status qaws_surface_create_fillet_ex(qaws_surface** out_surface, qaws_surface_fillet_desc const* desc, qaws_allocator const* allocator);
qaws_status qaws_surface_create_fillet_with_trims_ex(qaws_surface** out_fillet, qaws_curve** out_trim_a,
    qaws_curve** out_trim_b, qaws_surface_fillet_desc const* desc, qaws_allocator const* allocator);
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
