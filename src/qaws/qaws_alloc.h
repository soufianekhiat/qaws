#ifndef QAWS_ALLOC_H
#define QAWS_ALLOC_H

#include "qaws_types.h"
#include "qaws_status.h"
#include <stddef.h>
#include "qaws_bezier.h"
#include "qaws_hermite.h"
#include "qaws_catmull_rom.h"
#include "qaws_bspline.h"
#include "qaws_nurbs.h"
#include "qaws_trajectory.h"
#include "qaws_rational_bezier.h"
#include "qaws_yuksel.h"
#include "qaws_arc.h"
#include "qaws_polynomial.h"
#include "qaws_clothoid.h"
#include "qaws_composite.h"
#include "qaws_subdivision.h"
#include "qaws_surface_bezier.h"
#include "qaws_surface_bspline.h"
#include "qaws_surface_nurbs.h"
#include "qaws_surface_bilinear.h"
#include "qaws_surface_biquadratic.h"
#include "qaws_surface_coons.h"
#include "qaws_surface_extrusion.h"
#include "qaws_surface_offset.h"
#include "qaws_surface_pipe.h"
#include "qaws_surface_revolution.h"
#include "qaws_surface_ruled.h"
#include "qaws_surface_swept.h"
#include "qaws_surface_loft.h"
#include "qaws_surface_gordon.h"
#include "qaws_surface_trim.h"
#include "qaws_surface_tspline.h"
#include "qaws_surface_subdiv.h"
#include "qaws_surface_fillet.h"

/*
 * Allocator-aware creation functions (_ex variants).
 *
 * These behave identically to their non-_ex counterparts but route
 * all internal allocations through the provided qaws_allocator.
 * Pass NULL for allocator to use the default malloc/free.
 *
 * Curves created with a custom allocator MUST be destroyed with
 * qaws_curve_destroy (which reads the stored allocator). The
 * allocator must remain valid until the curve is destroyed.
 */

/*
 * Where the library's memory comes from. Without an explicit allocator (a NULL argument to the _ex functions, and all
 * internal scratch and results), memory comes from the calling thread's allocator if set, else the process default if
 * set, else the C allocator (aligned to QAWS_MEMORY_ALIGNMENT).
 *
 * Lifetimes: every setting and every _ex argument is copied (callbacks and user_data). A descriptor may be a local that
 * goes out of scope; each block keeps the copy it was made with and is freed and grown by it, whatever the settings at
 * that time and on whichever thread. What user_data points to must outlive every block of that allocator.
 *
 * The default is write-once in effect: it can be set (or cleared with NULL) until the library first allocates from it;
 * after that the call returns QAWS_STATUS_UNSUPPORTED_OPERATION and nothing changes. Reading it then is a plain load,
 * so allocations take no lock.
 *
 * A thread's allocator applies to that thread only. qaws_exchange_thread_allocator installs one and returns the one it
 * replaces (zeroed: none) so a host can nest settings and restore them, on every return path, with
 * qaws_set_thread_allocator(&previous). The library never changes it.
 */
qaws_status qaws_set_default_allocator(qaws_allocator const* allocator);   /* NULL: the C allocator */
void qaws_set_thread_allocator(qaws_allocator const* allocator);            /* this thread; NULL: the default */
void qaws_get_thread_allocator(qaws_allocator* out_allocator);              /* zeroed: none */
void qaws_exchange_thread_allocator(qaws_allocator const* allocator, qaws_allocator* out_previous);

/*
 * An allocator over caller memory (sized interface): a mapped file, pinned memory, a host frame arena. The memory must
 * be CPU-accessible and coherent for the CPU; placing the library's data there does not move any computation to a GPU.
 *
 * Blocks are cut in order, aligned to QAWS_MEMORY_ALIGNMENT, each after a small header. The latest block grows and shrinks
 * in place. A freed block is marked, and the freed blocks at the top give their space back, so blocks released in any
 * last-in-first-out order (nested scratch) are reclaimed at once; others wait for the blocks above them, or the reset. It is a hard bound: when full, allocation fails and
 * calls return QAWS_STATUS_ALLOCATION_FAILURE; nothing falls back to the heap. One thread at a time (one arena per
 * worker with qaws_set_thread_allocator). Pass &arena.allocator wherever an allocator is taken.
 *
 * Lifetime: objects and results allocated in an arena live in it. qaws_arena_reset refuses (QAWS_STATUS_INVALID_ARGUMENT)
 * while blocks are still held (live > 0): destroy every arena-backed object and finish every call first. Freeing a block
 * the arena does not hold (another arena, an earlier generation, twice) is counted in invalid_frees and asserts in
 * debug builds.
 */
typedef struct qaws_arena
{
	qaws_allocator allocator;
	unsigned char* base;
	size_t size;
	size_t used;                        /* bytes up to the end of the latest block */
	size_t peak;                        /* most bytes used since init */
	size_t last;                        /* the latest block's header offset ((size_t)-1: none) */
	size_t live;                        /* blocks held */
	size_t failures;                    /* allocations refused for room */
	size_t invalid_frees;
	unsigned int generation;            /* resets so far */
} qaws_arena;

void qaws_arena_init(qaws_arena* arena, void* memory, size_t size);
qaws_status qaws_arena_reset(qaws_arena* arena);

/* --- Curve creation _ex variants --- */

qaws_status qaws_curve_create_bezier_ex(
	qaws_bezier_desc const* desc,
	qaws_allocator const* allocator,
	qaws_curve** out_curve);

qaws_status qaws_curve_create_hermite_ex(
	qaws_hermite_desc const* desc,
	qaws_allocator const* allocator,
	qaws_curve** out_curve);

qaws_status qaws_curve_create_catmull_rom_ex(
	qaws_catmull_rom_desc const* desc,
	qaws_allocator const* allocator,
	qaws_curve** out_curve);

qaws_status qaws_curve_create_bspline_ex(
	qaws_bspline_desc const* desc,
	qaws_allocator const* allocator,
	qaws_curve** out_curve);

qaws_status qaws_curve_create_nurbs_ex(
	qaws_nurbs_desc const* desc,
	qaws_allocator const* allocator,
	qaws_curve** out_curve);

qaws_status qaws_curve_create_trajectory_ex(
	qaws_trajectory_desc const* desc,
	qaws_allocator const* allocator,
	qaws_curve** out_curve);

qaws_status qaws_curve_create_rational_bezier_ex(
	qaws_rational_bezier_desc const* desc,
	qaws_allocator const* allocator,
	qaws_curve** out_curve);

qaws_status qaws_curve_create_yuksel_ex(
	qaws_yuksel_desc const* desc,
	qaws_allocator const* allocator,
	qaws_curve** out_curve);

qaws_status qaws_curve_create_arc_ex(
	qaws_arc_desc const* desc,
	qaws_allocator const* allocator,
	qaws_curve** out_curve);

qaws_status qaws_curve_create_polynomial_ex(
	qaws_polynomial_desc const* desc,
	qaws_allocator const* allocator,
	qaws_curve** out_curve);

qaws_status qaws_curve_create_clothoid_ex(
	qaws_clothoid_desc const* desc,
	qaws_allocator const* allocator,
	qaws_curve** out_curve);

qaws_status qaws_curve_create_composite_ex(
	qaws_composite_desc const* desc,
	qaws_allocator const* allocator,
	qaws_curve** out_curve);

qaws_status qaws_curve_create_subdivision_ex(
	qaws_subdivision_desc const* desc,
	qaws_allocator const* allocator,
	qaws_curve** out_curve);

/* --- Surface creation _ex variants --- */

qaws_status qaws_surface_create_bezier_ex(
	qaws_surface_bezier_desc const* desc,
	qaws_allocator const* allocator,
	qaws_surface** out_surface);

qaws_status qaws_surface_create_bspline_ex(
	qaws_surface_bspline_desc const* desc,
	qaws_allocator const* allocator,
	qaws_surface** out_surface);

qaws_status qaws_surface_create_nurbs_ex(
	qaws_surface_nurbs_desc const* desc,
	qaws_allocator const* allocator,
	qaws_surface** out_surface);

qaws_status qaws_surface_create_bilinear_ex(
	qaws_surface_bilinear_desc const* desc,
	qaws_allocator const* allocator,
	qaws_surface** out_surface);

qaws_status qaws_surface_create_biquadratic_ex(
	qaws_surface_biquadratic_desc const* desc,
	qaws_allocator const* allocator,
	qaws_surface** out_surface);

qaws_status qaws_surface_create_coons_ex(
	qaws_surface_coons_desc const* desc,
	qaws_allocator const* allocator,
	qaws_surface** out_surface);

qaws_status qaws_surface_create_extrusion_ex(
	qaws_surface_extrusion_desc const* desc,
	qaws_allocator const* allocator,
	qaws_surface** out_surface);

qaws_status qaws_surface_create_offset_ex(
	qaws_surface_offset_desc const* desc,
	qaws_allocator const* allocator,
	qaws_surface** out_surface);

qaws_status qaws_surface_create_pipe_ex(
	qaws_surface_pipe_desc const* desc,
	qaws_allocator const* allocator,
	qaws_surface** out_surface);

qaws_status qaws_surface_create_revolution_ex(
	qaws_surface_revolution_desc const* desc,
	qaws_allocator const* allocator,
	qaws_surface** out_surface);

qaws_status qaws_surface_create_ruled_ex(
	qaws_surface_ruled_desc const* desc,
	qaws_allocator const* allocator,
	qaws_surface** out_surface);

qaws_status qaws_surface_create_swept_ex(
	qaws_surface_swept_desc const* desc,
	qaws_allocator const* allocator,
	qaws_surface** out_surface);

qaws_status qaws_surface_create_loft_ex(
	qaws_surface_loft_desc const* desc,
	qaws_allocator const* allocator,
	qaws_surface** out_surface);

qaws_status qaws_surface_create_gordon_ex(
	qaws_surface_gordon_desc const* desc,
	qaws_allocator const* allocator,
	qaws_surface** out_surface);

qaws_status qaws_surface_create_trimmed_ex(
	qaws_surface_trim_desc const* desc,
	qaws_allocator const* allocator,
	qaws_surface** out_surface);

qaws_status qaws_surface_create_tspline_ex(
	qaws_surface_tspline_desc const* desc,
	qaws_allocator const* allocator,
	qaws_surface** out_surface);

qaws_status qaws_surface_create_subdiv_ex(
	qaws_surface_subdiv_desc const* desc,
	qaws_allocator const* allocator,
	qaws_surface** out_surface);

qaws_status qaws_surface_create_fillet_ex(
	qaws_surface_fillet_desc const* desc,
	qaws_allocator const* allocator,
	qaws_surface** out_surface);

/* the fillet with its trim curves on both surfaces; all three from the allocator */
qaws_status qaws_surface_create_fillet_with_trims_ex(
	qaws_surface_fillet_desc const* desc,
	qaws_allocator const* allocator,
	qaws_surface** out_fillet,
	qaws_curve** out_trim_a,
	qaws_curve** out_trim_b);

#endif /* QAWS_ALLOC_H */
