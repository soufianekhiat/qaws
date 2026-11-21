# API Reference: qaws_alloc.h, qaws_allocator.h

Where the library's memory comes from, and how hosts give it their own.

- `qaws_allocator.h` is precision-free, declared once and defined once in
  the library: the allocator descriptor, the process default, thread
  allocators, scratch allocators, kept memory and arenas are the same for
  the f32 and the f64 functions.
- `qaws_alloc.h` declares the allocator-aware `_ex` creation functions of
  every curve and surface family and includes `qaws_allocator.h`. These
  names exist as `name_f32` and `name_f64`; this page uses the unsuffixed
  names (see [precision.md](precision.md)).

See [memory.md](../memory.md) for the memory model.

---

## Where memory comes from

Without an explicit allocator (a NULL argument to the `_ex` functions, and
all internal scratch and results), memory comes from, in order:

1. the calling thread's allocator (`qaws_set_thread_allocator`), if set;
2. the process default (`qaws_set_default_allocator`), if set;
3. the C allocator, aligned to `QAWS_MEMORY_ALIGNMENT`.

An allocator passed to an `_ex` function takes precedence over all three
for that object. When an allocation returns NULL, the call fails with
`QAWS_STATUS_ALLOCATION_FAILURE`.

## The descriptor

```c
typedef struct qaws_allocator {
    void* (*alloc)(unsigned long size, void* user_data);
    void  (*dealloc)(void* ptr, void* user_data);
    void* user_data;
    void* (*alloc_sized)(size_t size, size_t alignment, void* user_data);
    void  (*dealloc_sized)(void* ptr, size_t size, size_t alignment, void* user_data);
    int   (*resize)(void* ptr, size_t old_size, size_t new_size, size_t alignment, void* user_data);
} qaws_allocator;

#define QAWS_MEMORY_ALIGNMENT 16

void* qaws_allocator_extended(unsigned long size, void* user_data);
```

Two interfaces:

- **Legacy.** Set `alloc`, `dealloc` and `user_data` only; the other
  members are never read. Sizes are `unsigned long` (32 bits on Windows):
  a larger request fails without calling `alloc`.
- **Extended.** Set `alloc = qaws_allocator_extended` (a marker that never
  allocates when called) and the sized members. `alloc_sized` receives the
  size and alignment; `dealloc_sized` receives the size and alignment the
  block was allocated with. `resize` is optional: it grows or shrinks a
  block in place and returns non-zero, or returns 0 to let the library
  allocate, copy and free.

Blocks must be aligned to `QAWS_MEMORY_ALIGNMENT` (16). A misaligned block
is given back and the call fails with `QAWS_STATUS_ALLOCATION_FAILURE`.

## Lifetimes

- **The descriptor** may be a short-lived local. Every setting and every
  `_ex` argument is copied (callbacks and `user_data`) into every block the
  library allocates and into the objects that allocate later.
- **The state `user_data` points to** must outlive every block the
  allocator gave.
- **Blocks** are freed and grown by the copy they were made with, whatever
  the settings at that time and on whichever thread.

## Settings

```c
qaws_status qaws_set_default_allocator(qaws_allocator const* allocator);   /* NULL: the C allocator */
void qaws_set_thread_allocator(qaws_allocator const* allocator);            /* this thread; NULL: the default */
void qaws_get_thread_allocator(qaws_allocator* out_allocator);              /* zeroed: none */
void qaws_exchange_thread_allocator(qaws_allocator* out_previous, qaws_allocator const* allocator);
```

- **`qaws_set_default_allocator`** is write-once in effect: the default can
  be set (or cleared with NULL) until the library first allocates from it.
  After that the call returns `QAWS_STATUS_UNSUPPORTED_OPERATION` and
  nothing changes. Reading the sealed default is a plain load, so
  allocations take no lock.
- **`qaws_set_thread_allocator`** applies to the calling thread only; the
  library never changes it.
- **`qaws_get_thread_allocator`** copies the thread's allocator (zeroed:
  none) into `out_allocator`.
- **`qaws_exchange_thread_allocator`** installs `allocator` and returns the
  one it replaces in `out_previous` (zeroed: none; `out_previous` may be
  NULL), so a host can nest settings and restore them on every return path
  with `qaws_set_thread_allocator(&previous)`.

```c
qaws_allocator previous;
qaws_exchange_thread_allocator(&previous, &frame_arena.allocator);
status = qaws_curve_offset_2d(out, &count, curve, cap, d, 1);
qaws_set_thread_allocator(&previous);   /* on every return path */
```

## Scratch apart from results

```c
void qaws_exchange_thread_scratch_allocator(qaws_allocator* out_previous, qaws_allocator const* allocator);
```

Installs a scratch allocator on the calling thread (NULL: none) and returns
the one it replaces (zeroed: none; `out_previous` may be NULL). With one
installed, the temporaries of adaptive sampling, tessellation, meshing,
clipping and the repulsive solvers (blocks allocated and released within
one call) come from it. Results and objects still come from the `_ex`
argument, the thread allocator or the default, so a frame arena can serve
as scratch and be reset after each call while results stay. Other
operations take their scratch from the order above. Restore it on every
return path.

```c
qaws_allocator previous;
qaws_exchange_thread_scratch_allocator(&previous, &frame_arena.allocator);
status = qaws_curve_sample_2d(points, &count, curve, capacity, &desc);
qaws_exchange_thread_scratch_allocator(NULL, &previous);
qaws_arena_reset(&frame_arena);   /* the arena holds nothing after the call */
```

## Memory kept for reuse

```c
void qaws_set_kept_memory_limit(size_t bytes);   /* 0: none kept; default 16 MB */
void qaws_release_kept_memory(void);
```

Blocks of the C allocator (no custom allocator in effect) from 4 KB to
1 MB that the library frees are kept, up to the limit in all, and handed to
the next requests of their size, so repeated calls do not take their memory
from the system each time. Such blocks are made at their size class's
bytes, a quarter octave apart. The limit is process-wide and may be changed
at any time; lowering it frees the kept blocks past it.
`qaws_release_kept_memory` frees every kept block and leaves the limit.
Both calls are thread-safe. Custom allocators are not affected.

## Caller memory: qaws_arena

```c
typedef struct qaws_arena {
    qaws_allocator allocator;      /* pass &arena.allocator wherever an allocator is taken */
    unsigned char* base;
    size_t size;
    size_t used;                   /* bytes up to the end of the latest block */
    size_t peak;                   /* most bytes used since init */
    size_t last;                   /* the latest block's header offset ((size_t)-1: none) */
    size_t live;                   /* blocks held */
    size_t failures;               /* allocations refused for room */
    size_t invalid_frees;          /* frees of blocks it does not hold */
    unsigned int generation;       /* resets so far */
} qaws_arena;

void qaws_arena_init(qaws_arena* arena, void* memory, size_t size);
qaws_status qaws_arena_reset(qaws_arena* arena);
```

An extended allocator over caller memory: a mapped file, pinned memory, a
host frame arena.

- **`qaws_arena_init`** sets up the arena over `memory` (aligned up to 16
  bytes inside it) and fills `allocator`.
- **CPU access.** The memory must be CPU-accessible and coherent for the
  CPU; placing the library's data there does not move any computation to a
  GPU.
- **Allocation.** Blocks are cut in order, aligned to
  `QAWS_MEMORY_ALIGNMENT`, each after a small header. The latest block
  grows and shrinks in place. A freed block is marked, and the freed blocks
  at the top give their space back, so blocks released in last-in-first-out
  order (nested scratch) are reclaimed at once; others wait for the blocks
  above them, or the reset.
- **Hard bound.** When full, allocation fails and calls return
  `QAWS_STATUS_ALLOCATION_FAILURE`; nothing falls back to the heap.
- **`qaws_arena_reset`** refuses (`QAWS_STATUS_INVALID_ARGUMENT`) while
  blocks are held (`live > 0`): destroy every arena-backed object and
  finish every call first. Freeing a block the arena does not hold (another
  arena, an earlier generation, twice) is counted in `invalid_frees` and
  asserts in debug builds.
- **Threads.** One thread at a time. For workers, use one arena per thread
  with `qaws_set_thread_allocator` or `qaws_exchange_thread_allocator`.

Objects and results allocated in an arena live in it. To keep an object
beyond a reset, create it with another allocator (the `_ex` argument) while
scratch uses the arena.

## _ex creation functions

Each behaves as its constructor without `_ex` (which calls it with NULL)
and routes every allocation the constructor makes through `allocator`:
the object, its stored data, and construction scratch. Pass NULL to use the
order above. Objects keep a copy of the descriptor and are destroyed with
`qaws_curve_destroy` or `qaws_surface_destroy` as usual; the state behind
`user_data` must stay valid until then. On failure, `*out_curve` /
`*out_surface` is not set to an object and nothing stays allocated.

### Curves

| Function | Descriptor ([family page](curve.md)) |
|---|---|
| `qaws_status qaws_curve_create_bezier_ex(qaws_curve** out_curve, qaws_bezier_desc const* desc, qaws_allocator const* allocator);` | [bezier.md](bezier.md) |
| `qaws_status qaws_curve_create_hermite_ex(qaws_curve** out_curve, qaws_hermite_desc const* desc, qaws_allocator const* allocator);` | [hermite.md](hermite.md) |
| `qaws_status qaws_curve_create_catmull_rom_ex(qaws_curve** out_curve, qaws_catmull_rom_desc const* desc, qaws_allocator const* allocator);` | [catmull_rom.md](catmull_rom.md) |
| `qaws_status qaws_curve_create_bspline_ex(qaws_curve** out_curve, qaws_bspline_desc const* desc, qaws_allocator const* allocator);` | [bspline.md](bspline.md) |
| `qaws_status qaws_curve_create_nurbs_ex(qaws_curve** out_curve, qaws_nurbs_desc const* desc, qaws_allocator const* allocator);` | [nurbs.md](nurbs.md) |
| `qaws_status qaws_curve_create_trajectory_ex(qaws_curve** out_curve, qaws_trajectory_desc const* desc, qaws_allocator const* allocator);` | [trajectory.md](trajectory.md) |
| `qaws_status qaws_curve_create_rational_bezier_ex(qaws_curve** out_curve, qaws_rational_bezier_desc const* desc, qaws_allocator const* allocator);` | [rational_bezier.md](rational_bezier.md) |
| `qaws_status qaws_curve_create_yuksel_ex(qaws_curve** out_curve, qaws_yuksel_desc const* desc, qaws_allocator const* allocator);` | [yuksel.md](yuksel.md) |
| `qaws_status qaws_curve_create_arc_ex(qaws_curve** out_curve, qaws_arc_desc const* desc, qaws_allocator const* allocator);` | [arc.md](arc.md) |
| `qaws_status qaws_curve_create_polynomial_ex(qaws_curve** out_curve, qaws_polynomial_desc const* desc, qaws_allocator const* allocator);` | [polynomial.md](polynomial.md) |
| `qaws_status qaws_curve_create_clothoid_ex(qaws_curve** out_curve, qaws_clothoid_desc const* desc, qaws_allocator const* allocator);` | [clothoid.md](clothoid.md) |
| `qaws_status qaws_curve_create_composite_ex(qaws_curve** out_curve, qaws_composite_desc const* desc, qaws_allocator const* allocator);` | [composite.md](composite.md) |
| `qaws_status qaws_curve_create_subdivision_ex(qaws_curve** out_curve, qaws_subdivision_desc const* desc, qaws_allocator const* allocator);` | [subdivision.md](subdivision.md) |

A composite owns its segments only on success; the segments keep the
allocators they were made with. The subdivision curve's refinement levels
come from the allocator.

### Surfaces

| Function | Descriptor |
|---|---|
| `qaws_status qaws_surface_create_bezier_ex(qaws_surface** out_surface, qaws_surface_bezier_desc const* desc, qaws_allocator const* allocator);` | [surface_families.md](surface_families.md) |
| `qaws_status qaws_surface_create_bspline_ex(qaws_surface** out_surface, qaws_surface_bspline_desc const* desc, qaws_allocator const* allocator);` | [surface_families.md](surface_families.md) |
| `qaws_status qaws_surface_create_nurbs_ex(qaws_surface** out_surface, qaws_surface_nurbs_desc const* desc, qaws_allocator const* allocator);` | [surface_families.md](surface_families.md) |
| `qaws_status qaws_surface_create_bilinear_ex(qaws_surface** out_surface, qaws_surface_bilinear_desc const* desc, qaws_allocator const* allocator);` | [surface_families.md](surface_families.md) |
| `qaws_status qaws_surface_create_biquadratic_ex(qaws_surface** out_surface, qaws_surface_biquadratic_desc const* desc, qaws_allocator const* allocator);` | [surface_families.md](surface_families.md) |
| `qaws_status qaws_surface_create_coons_ex(qaws_surface** out_surface, qaws_surface_coons_desc const* desc, qaws_allocator const* allocator);` | [surface_families.md](surface_families.md) |
| `qaws_status qaws_surface_create_extrusion_ex(qaws_surface** out_surface, qaws_surface_extrusion_desc const* desc, qaws_allocator const* allocator);` | [surface_families.md](surface_families.md) |
| `qaws_status qaws_surface_create_offset_ex(qaws_surface** out_surface, qaws_surface_offset_desc const* desc, qaws_allocator const* allocator);` | [surface_families.md](surface_families.md) |
| `qaws_status qaws_surface_create_pipe_ex(qaws_surface** out_surface, qaws_surface_pipe_desc const* desc, qaws_allocator const* allocator);` | [surface_families.md](surface_families.md) |
| `qaws_status qaws_surface_create_revolution_ex(qaws_surface** out_surface, qaws_surface_revolution_desc const* desc, qaws_allocator const* allocator);` | [surface_families.md](surface_families.md) |
| `qaws_status qaws_surface_create_ruled_ex(qaws_surface** out_surface, qaws_surface_ruled_desc const* desc, qaws_allocator const* allocator);` | [surface_families.md](surface_families.md) |
| `qaws_status qaws_surface_create_swept_ex(qaws_surface** out_surface, qaws_surface_swept_desc const* desc, qaws_allocator const* allocator);` | [surface_families.md](surface_families.md) |
| `qaws_status qaws_surface_create_loft_ex(qaws_surface** out_surface, qaws_surface_loft_desc const* desc, qaws_allocator const* allocator);` | [surface_families.md](surface_families.md) |
| `qaws_status qaws_surface_create_gordon_ex(qaws_surface** out_surface, qaws_surface_gordon_desc const* desc, qaws_allocator const* allocator);` | [surface_families.md](surface_families.md) |
| `qaws_status qaws_surface_create_trimmed_ex(qaws_surface** out_surface, qaws_surface_trim_desc const* desc, qaws_allocator const* allocator);` | [surface_trim.md](surface_trim.md) |
| `qaws_status qaws_surface_create_tspline_ex(qaws_surface** out_surface, qaws_surface_tspline_desc const* desc, qaws_allocator const* allocator);` | [surface_tspline.md](surface_tspline.md) |
| `qaws_status qaws_surface_create_subdiv_ex(qaws_surface** out_surface, qaws_surface_subdiv_desc const* desc, qaws_allocator const* allocator);` | [surface_subdiv.md](surface_subdiv.md) |

The curves a surface borrows are not counted. Loft and Gordon evaluation
needs scratch per section curve: on the stack up to 32 sections, beyond
that from the surface's own allocator, so a host that evaluates a large
loft or Gordon surface from several threads needs an allocator that is
safe across threads.

## Derived curves

An operation that makes curves from a source curve gives its results the
source curve's allocator: reverse, extract, split, join, conversion, degree
elevation and reduction, merging a chain (the first curve's allocator),
transforms, and the arc-length and length-matched wrappers. Curves built
from raw arrays (polyline and offset outputs, clipping, skeletons, stacks
blending several curves) take the default order instead.
