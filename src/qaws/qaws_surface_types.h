#ifndef QAWS_SURFACE_TYPES_H
#define QAWS_SURFACE_TYPES_H

#include "qaws_types.h"
#include "qaws_status.h"

/*
 * Thread safety: qaws_surface objects are immutable after creation.
 * Multiple threads may concurrently evaluate or sample the same surface.
 */

typedef enum qaws_surface_kind
{
	QAWS_SURFACE_KIND_INVALID = 0,
	QAWS_SURFACE_KIND_BEZIER,
	QAWS_SURFACE_KIND_BSPLINE,
	QAWS_SURFACE_KIND_NURBS,
	QAWS_SURFACE_KIND_SWEPT,
	QAWS_SURFACE_KIND_RULED,
	QAWS_SURFACE_KIND_COONS,
	QAWS_SURFACE_KIND_EXTRUSION,
	QAWS_SURFACE_KIND_REVOLUTION,
	QAWS_SURFACE_KIND_PIPE,
	QAWS_SURFACE_KIND_BILINEAR,
	QAWS_SURFACE_KIND_BIQUADRATIC,
	QAWS_SURFACE_KIND_LOFT,
	QAWS_SURFACE_KIND_GORDON,
	QAWS_SURFACE_KIND_OFFSET,
	QAWS_SURFACE_KIND_TRIMMED
} qaws_surface_kind;

typedef enum qaws_surface_eval_flags
{
	QAWS_SURFACE_EVAL_NONE     = 0,
	QAWS_SURFACE_EVAL_POSITION = 1 << 0,
	QAWS_SURFACE_EVAL_DU       = 1 << 1,  /* partial derivative w.r.t. u */
	QAWS_SURFACE_EVAL_DV       = 1 << 2,  /* partial derivative w.r.t. v */
	QAWS_SURFACE_EVAL_NORMAL   = 1 << 3,  /* unit surface normal (du x dv) */
	QAWS_SURFACE_EVAL_DUU      = 1 << 4,  /* second partial d^2/du^2 */
	QAWS_SURFACE_EVAL_DUV      = 1 << 5,  /* mixed partial d^2/dudv */
	QAWS_SURFACE_EVAL_DVV      = 1 << 6   /* second partial d^2/dv^2 */
} qaws_surface_eval_flags;

typedef struct qaws_surface_eval_result
{
	qaws_vec3 position;
	qaws_vec3 du;       /* dS/du */
	qaws_vec3 dv;       /* dS/dv */
	qaws_vec3 normal;   /* unit normal */
	qaws_vec3 duu;      /* d^2S/du^2 */
	qaws_vec3 duv;      /* d^2S/dudv */
	qaws_vec3 dvv;      /* d^2S/dv^2 */
	unsigned int valid_flags;
} qaws_surface_eval_result;

typedef struct qaws_surface qaws_surface;

#endif /* QAWS_SURFACE_TYPES_H */
