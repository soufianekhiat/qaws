#ifndef QAWS_DIFF_FUNCTIONALS_H
#define QAWS_DIFF_FUNCTIONALS_H

#include "qaws_diff_types.h"
#include "qaws_surface_types.h"

/*
 * Integral functionals of curves and surfaces with exact derivatives of
 * their quadrature:
 *
 *   value, tangent, tangent2  along one parameter direction
 *   gradient                  accumulated (+=) into parameter views
 *   hvp                       Hessian-vector product (+=), direct for
 *                             families linear in their fields
 *                             (QAWS_CAP_LINEAR); others return
 *                             QAWS_STATUS_UNSUPPORTED_OPERATION and the
 *                             caller composes HVPs from tangents/adjoints.
 *
 * Curves integrate span by span with Gauss-Legendre points
 * (quadrature = points per span, 0 = 6). Surfaces integrate a grid of
 * cells with Gauss points (quadrature = cells per direction, 0 = 8;
 * 4 x 4 points per cell).
 */

typedef enum qaws_curve_functional
{
	QAWS_FUNCTIONAL_LENGTH = 0,            /* integral of |C'| dt */
	QAWS_FUNCTIONAL_BENDING,               /* integral of |C''|^2 dt */
	QAWS_FUNCTIONAL_CURVATURE_SQUARED      /* integral of kappa^2 ds */
} qaws_curve_functional;

typedef enum qaws_surface_functional
{
	QAWS_FUNCTIONAL_AREA = 0,              /* integral of |Su x Sv| du dv */
	QAWS_FUNCTIONAL_THIN_PLATE,            /* integral of |Suu|^2 + 2|Suv|^2 + |Svv|^2 */
	QAWS_FUNCTIONAL_WILLMORE               /* integral of H^2 dA */
} qaws_surface_functional;

qaws_status qaws_curve_functional_eval(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_curve_functional functional,
	unsigned int quadrature,
	qaws_diff_views const* direction,
	qaws_scalar* out_value,
	qaws_scalar* out_tangent,
	qaws_scalar* out_tangent2);

qaws_status qaws_curve_functional_gradient(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_curve_functional functional,
	unsigned int quadrature,
	qaws_diff_views* gradient,
	qaws_scalar* out_value);

qaws_status qaws_curve_functional_hvp(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_curve_functional functional,
	unsigned int quadrature,
	qaws_diff_views const* direction,
	qaws_diff_views* out_hv);

qaws_status qaws_surface_functional_eval(
	qaws_diff_context const* ctx,
	qaws_surface const* surface,
	qaws_surface_functional functional,
	unsigned int quadrature,
	qaws_diff_views const* direction,
	qaws_scalar* out_value,
	qaws_scalar* out_tangent,
	qaws_scalar* out_tangent2);

qaws_status qaws_surface_functional_gradient(
	qaws_diff_context const* ctx,
	qaws_surface const* surface,
	qaws_surface_functional functional,
	unsigned int quadrature,
	qaws_diff_views* gradient,
	qaws_scalar* out_value);

qaws_status qaws_surface_functional_hvp(
	qaws_diff_context const* ctx,
	qaws_surface const* surface,
	qaws_surface_functional functional,
	unsigned int quadrature,
	qaws_diff_views const* direction,
	qaws_diff_views* out_hv);

#endif /* QAWS_DIFF_FUNCTIONALS_H */
