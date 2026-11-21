# Memory and Allocation

This page explains where qaws gets its memory, who owns each object, and how to run without the heap. API pages:
[alloc](api/alloc.md) (`qaws_alloc.h` and `qaws_allocator.h`) and [inline](api/inline.md).

The allocator, its settings and the arena are precision-free. They are declared and defined once, so the same default
allocator, thread allocators and arenas serve both the f32 and the f64 functions (see [precision.md](precision.md)).

## Where memory comes from

Every allocation goes through one internal layer. For each block, the allocator is chosen in this order:

1. The allocator passed explicitly to an `_ex` function or a `prepare` function, if it is not NULL.
2. The calling thread's allocator (`qaws_set_thread_allocator`), if one is set.
3. The process default (`qaws_set_default_allocator`), if one is set.
4. The C allocator, aligned to `QAWS_MEMORY_ALIGNMENT` (16).

This order covers objects, results and internal scratch alike. One exception applies: with a thread **scratch**
allocator installed (`qaws_exchange_thread_scratch_allocator`), the temporaries of adaptive sampling, tessellation,
meshing, clipping and the repulsive solvers come from it. Those are blocks allocated and released within one call.
Results still follow the order above, so a frame arena can serve as scratch and be reset after each call while the
results stay.

A failed allocation makes the call return `QAWS_STATUS_ALLOCATION_FAILURE`. Size products are checked, and a size that
would overflow is refused instead of wrapping.

## The allocator descriptor

A `qaws_allocator` has two forms:

- **Legacy:** `alloc(size, user_data)` and `dealloc(ptr, user_data)`. Sizes are `unsigned long`, which is 32 bits on
  Windows, so larger requests fail.
- **Sized (extended):** set `alloc = qaws_allocator_extended` as the marker, then fill in
  `alloc_sized(size, alignment, user)`, `dealloc_sized(ptr, size, alignment, user)` and, optionally,
  `resize(...)`. `resize` grows or shrinks a block in place, or returns 0 to let the library allocate, copy and free.

Blocks must be 16-byte aligned. A misaligned block is given back and the call fails.

**Lifetimes.** The library copies the descriptor (its callbacks and `user_data`) into every block it allocates and into
every object that allocates later. The descriptor itself can therefore be a short-lived local. What `user_data` points
to must outlive every block of that allocator. A block is always freed and grown by the copy it was made with,
whatever the settings are at that time and on whichever thread.

## Settings

- **Process default.** `qaws_set_default_allocator` is write-once in effect. It can be set, or cleared with NULL,
  until the library first allocates from it. After that it returns `UNSUPPORTED_OPERATION` and nothing changes. Reading
  it is a plain load, so allocations take no lock.
- **Thread allocator.** `qaws_set_thread_allocator`, `qaws_get_thread_allocator` and
  `qaws_exchange_thread_allocator` apply to the calling thread only. Exchange returns the previous setting, so a host
  can nest settings and restore them on every return path. The library never changes a thread's setting itself.
- **Kept memory.** When no custom allocator is in effect, blocks of the C allocator from 4 KB to 1 MB are kept for
  reuse when they are freed, up to a process-wide limit (16 MB by default). This saves repeated calls from taking pages
  from the system each time. `qaws_set_kept_memory_limit` changes the limit at any time (0 keeps nothing), and
  `qaws_release_kept_memory` frees every kept block. Both are thread-safe.

## Arenas

A `qaws_arena` is an allocator over caller memory, such as a mapped file, pinned memory or a frame buffer. The memory
must be CPU-accessible. Pass `&arena.allocator` anywhere an allocator is accepted.

- **Allocation.** Blocks are cut in order, each after a small header. The latest block grows and shrinks in place.
- **Freeing.** Freed blocks at the top give their space back at once, so last-in-first-out scratch is reclaimed
  immediately. Other freed blocks wait for the blocks above them, or for the reset.
- **Hard bound.** When the arena is full, allocation fails, and nothing falls back to the heap.
- **One thread at a time.** For parallel work, give each worker its own arena through `qaws_set_thread_allocator`.
- **Reset.** `qaws_arena_reset` refuses (`INVALID_ARGUMENT`) while blocks are still held (`live > 0`). Destroy every
  arena-backed object first.
- **Statistics.** The struct exposes `used`, `peak`, `live`, `failures`, `invalid_frees` and `generation`. A free of a
  block the arena does not hold is counted in `invalid_frees` and asserts in debug builds.

## Ownership and destroy rules

| Object | Made by | Released by | Notes |
|---|---|---|---|
| curve | `qaws_curve_create_*`, `_ex`, operations | `qaws_curve_destroy` | freed by the allocator stored in it; that allocator must stay valid until then |
| surface | `qaws_surface_create_*`, `_ex` | `qaws_surface_destroy` | same |
| prepared curve / surface | `qaws_*_prepare` | `qaws_*_prepared_destroy` | `prepare_into` objects live in caller memory; their destroy does nothing |
| curve / surface set, seams, shatter, clip result, B-rep shell, traversal | their `_create` | their `_destroy` | results inside (pieces, path views) belong to the object |
| inline curve | `qaws_curve_init_*_inline` | nothing | never call `qaws_curve_destroy` on it |

Further rules:

- **Borrowing.** Many objects borrow other objects instead of owning them: arc-length wrappers borrow their source,
  trimmed surfaces their base and trim curves, B-rep shells their surfaces and curves, sets their inputs, and borrowed
  curves the host's arrays. Whatever is borrowed must outlive the borrower.
- **Derived curves.** A curve derived from one source (split, extract, reverse, join, conversions) is allocated with
  the source's allocator.
- **NULL is safe.** Every `destroy` function accepts NULL.
- **No destroy during use.** Never destroy an object while another thread is using it.

## Inline curves: no heap at all

A `qaws_curve_inline` is a 512-byte union that holds a complete single-span Bezier or polynomial curve (2D or 3D,
degree up to 7) in its own storage. `qaws_curve_init_bezier_inline` and `qaws_curve_init_polynomial_inline` fill it
without allocating. `qaws_curve_inline_get_curve` returns a `qaws_curve*` that works with every evaluation and
inspection function.

The struct is self-contained, so it can be copied with `memcpy` or plain assignment. After copying, call
`qaws_curve_inline_get_curve(&copy)` once to re-point the copy at its own storage, before other threads use it. A
pointer taken from the original keeps reading the original.

## Running without allocation

There are three ways to avoid heap allocation:

- **No allocation at all.** Inline curves, `qaws_curve_prepare_into` / `qaws_surface_prepare_into` with prepared-object
  evaluation, and the `qaws_prepare.h` functions (caller-allocated outputs) never allocate.
- **No allocation during evaluation.** Most evaluation of an existing object does not allocate. When it needs
  temporaries, they come from a small stack buffer or the object's own allocator.
- **A bounded budget.** To put everything else under a fixed budget, install an arena as the thread allocator for the
  duration of the work.
