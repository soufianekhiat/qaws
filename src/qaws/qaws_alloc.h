#ifndef QAWS_ALLOC_H
#define QAWS_ALLOC_H

#include "qaws_types.h"
#include "qaws_status.h"
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
