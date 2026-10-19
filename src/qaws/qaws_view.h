#ifndef QAWS_VIEW_H
#define QAWS_VIEW_H

#include "qaws_types.h"
#include "qaws_status.h"
#include <stddef.h>

/*
 * Borrowed geometry: curves evaluated from host memory without copying it.
 *
 *   qaws_array_view      a typed strided array the host owns: f32 or f64,
 *                        packed (any stride, interleaved records included)
 *                        or planar (one pointer per component).
 *   qaws_curve_view      a B-spline, NURBS or Bezier curve described by array
 *                        views: points, optional weights, optional knots.
 *   qaws_curve_prepared  derived evaluation data (the knots in the computing
 *                        precision, the span table), owned by the prepared
 *                        object; the points and weights stay borrowed.
 *                        Immutable: any number of threads may evaluate it.
 *   borrowed qaws_curve  qaws_curve_create_borrowed: a qaws_curve over a
 *                        view (kind QAWS_CURVE_KIND_VIEW), usable by every
 *                        operation that evaluates curves (sampling, arc
 *                        length, closest points, frames, ...); extract,
 *                        split and reverse return owned curves (made through
 *                        an owned copy). Operations that read a family's own
 *                        arrays (control point and knot inspection,
 *                        transforms, derivatives with respect to control
 *                        points) refuse it: the host holds those arrays.
 *   owned copy           qaws_curve_create_from_view: an owned B-spline or
 *                        NURBS curve from a view.
 *
 * Precision: views hold f32 or f64 data whatever the library's qaws_scalar;
 * values are converted one at a time as they are read (no converted copy).
 *
 * Lifetime contract:
 *   - qaws never writes to, frees or keeps host memory a view points to,
 *     except as stated here for prepared objects and borrowed curves.
 *   - Functions taking a view borrow it for the call only.
 *   - A prepared object or a borrowed curve borrows the view's arrays for
 *     its whole life: they must stay valid (not freed, not moved) until it is
 *     destroyed, and must not change while an evaluation with it runs.
 *   - Between evaluations the host may change what the object does not
 *     depend on: qaws_curve_prepared_dependencies says what it does. A curve
 *     prepared from knots depends on the knots (and the layout: pointers,
 *     strides, counts, types, degree), not on the points' or weights'
 *     values, which may be rewritten in place (an animation, a solver
 *     writing control points) without preparing again. A change to anything
 *     depended on requires destroying and preparing again.
 *   - The prepared object's own memory comes from the given allocator (or
 *     malloc), or from caller memory (qaws_curve_prepare_into: no
 *     allocation at all). Evaluation never allocates.
 */

typedef enum qaws_element_type
{
	QAWS_ELEMENT_F32 = 1,
	QAWS_ELEMENT_F64 = 2
} qaws_element_type;

/* Packed (data non-NULL): element i's component c at data + i * stride + c * size. Planar (data NULL): at planes[c] + i * stride.
   size: 4 (f32) or 8 (f64); stride 0 means tight (packed: components * size; planar: size). */
typedef struct qaws_array_view
{
	void const* data;
	void const* planes[4];
	qaws_element_type type;
	unsigned int components;            /* 1 to 4 */
	unsigned int count;
	size_t stride;                      /* bytes */
} qaws_array_view;

/* The same, writable: outputs */
typedef struct qaws_array_view_mut
{
	void* data;
	void* planes[4];
	qaws_element_type type;
	unsigned int components;
	unsigned int count;
	size_t stride;
} qaws_array_view_mut;

/* QAWS_STATUS_INVALID_ARGUMENT unless: a known type; 1 to 4 components; when count > 0, non-NULL pointers aligned to the scalar's
   size and a stride that is a multiple of it, at least an element (packed) or a value (planar). */
qaws_status qaws_array_view_validate(qaws_array_view const* view);

/* QAWS_STATUS_OUT_OF_RANGE unless every byte the view reads lies in [base, base + bytes) (a host buffer's extent) */
qaws_status qaws_array_view_check_extent(qaws_array_view const* view, void const* base, size_t bytes);

#define QAWS_VIEW_MAX_DEGREE 31

typedef struct qaws_curve_view
{
	qaws_dimension dimension;           /* 2 or 3 */
	unsigned int degree;                /* 1 to QAWS_VIEW_MAX_DEGREE */
	qaws_array_view points;             /* n points: dimension components, or dimension + 1 homogeneous (w x, w y, [w z,] w) */
	qaws_array_view weights;            /* count 0: none; else n weights > 0 (with Cartesian points) */
	qaws_array_view knots;              /* count 0: clamped uniform 0, 1, ..., n - degree, as qaws_curve_create_bspline (n = degree + 1:
	                                       a Bezier curve on [0, 1]); else n + degree + 1 non-decreasing values, the domain
	                                       [knots[degree], knots[n]] */
} qaws_curve_view;

/* The view's arrays and their values: counts, layouts, knots non-decreasing over a non-empty domain, weights (or homogeneous
   w) positive and finite. Every function taking a view checks this. */
qaws_status qaws_curve_view_validate(qaws_curve_view const* view);

typedef struct qaws_curve_prepared qaws_curve_prepared;

/* what a prepared object reads from host memory and so depends on */
#define QAWS_DEPENDS_POINTS  1u         /* the points' values */
#define QAWS_DEPENDS_WEIGHTS 2u         /* the weights' values */
#define QAWS_DEPENDS_KNOTS   4u         /* the knots' values */
#define QAWS_DEPENDS_LAYOUT  8u         /* pointers, strides, counts, types, the degree */

qaws_status qaws_curve_prepare(qaws_curve_view const* view, qaws_allocator const* allocator, qaws_curve_prepared** out_prepared);

/* bytes (and the alignment, 16) qaws_curve_prepare_into needs for this view; 0 when the view is invalid */
size_t qaws_curve_prepared_size(qaws_curve_view const* view);

/* into caller memory (no allocation): memory 16-byte aligned, at least qaws_curve_prepared_size bytes; the object lives there
   until the caller reuses it (qaws_curve_prepared_destroy is a no-op for it, and may be called) */
qaws_status qaws_curve_prepare_into(qaws_curve_view const* view, void* memory, size_t bytes, qaws_curve_prepared** out_prepared);

void qaws_curve_prepared_destroy(qaws_curve_prepared* prepared);

/* QAWS_DEPENDS_* flags: what the host must not change without preparing again (the points' and weights' values never) */
unsigned int qaws_curve_prepared_dependencies(qaws_curve_prepared const* prepared);

qaws_range qaws_curve_prepared_get_range(qaws_curve_prepared const* prepared);

/* One parameter (clamped to the domain); derivatives with respect to it, as qaws_curve_evaluate_*; thread-safe, no allocation */
qaws_status qaws_curve_prepared_evaluate_2d(qaws_curve_prepared const* prepared, qaws_scalar parameter, unsigned int eval_flags,
	qaws_eval_result_2d* out_result);
qaws_status qaws_curve_prepared_evaluate_3d(qaws_curve_prepared const* prepared, qaws_scalar parameter, unsigned int eval_flags,
	qaws_eval_result_3d* out_result);

/* Outputs of a range evaluation: each written when it has a pointer (data or planes[0]), dimension components, at least `end`
   elements, f32 or f64. */
typedef struct qaws_curve_outputs
{
	qaws_array_view_mut position, d1, d2, d3;
} qaws_curve_outputs;

/* Parameters begin to end - 1 of a view of one component, each evaluated into the outputs' same elements. Disjoint ranges write
   disjoint elements: a host's workers may evaluate ranges of one call at once. No allocation; the first failing status, the
   outputs written up to it. */
qaws_status qaws_curve_prepared_evaluate_range(qaws_curve_prepared const* prepared, qaws_array_view const* parameters, unsigned int begin,
	unsigned int end, qaws_curve_outputs const* outputs);

/* A qaws_curve over the view without copying its points, weights or knots (a prepared object inside, from the allocator); the same
   contract as a prepared object for as long as the curve lives. Kind QAWS_CURVE_KIND_VIEW. */
qaws_status qaws_curve_create_borrowed(qaws_curve_view const* view, qaws_allocator const* allocator, qaws_curve** out_curve);

/* An owned copy: a B-spline curve (polynomial views) or a NURBS curve (weights or homogeneous points, divided by w), with the
   view's knots (or the clamped uniform ones), converted to qaws_scalar; the view is borrowed for the call only. */
qaws_status qaws_curve_create_from_view(qaws_curve_view const* view, qaws_allocator const* allocator, qaws_curve** out_curve);

#endif /* QAWS_VIEW_H */
