#ifndef QAWS_INTERNAL_DIFF_H
#define QAWS_INTERNAL_DIFF_H

#include "../qaws_diff_types.h"
#include "../qaws_surface_types.h"
#include "qaws_internal_types.h"
#include "qaws_internal_surface.h"

/*
 * Differential rules attached to a curve family.
 *
 * Linear families (output linear in every field) only implement
 * linear_support; the generic engine in qaws_diff.c derives primal,
 * tangent, tangent2 and adjoint from the exact basis weights.
 * Non-linear families implement tangent_span / adjoint_span directly.
 */
struct qaws_curve_diff_vtable
{
	unsigned int capabilities;   /* qaws_diff_capability bits */
	qaws_diff_class diff_class;  /* with respect to the evaluation coordinate */

	/* Fill out[0..cap) and return the total number of fields. */
	unsigned int (*describe_fields)(
		qaws_curve const* curve,
		qaws_field_desc* out,
		unsigned int capacity);

	/* Read-only access to the primal storage of a field. */
	qaws_status (*primal_field)(
		qaws_curve const* curve,
		qaws_diff_field field,
		qaws_scalar const** out_data,
		unsigned int* out_count,
		unsigned int* out_components);

	/* Exact basis weights for derivative orders 0..order with respect to
	   the global parameter, evaluated at span-local coordinate local_t. */
	qaws_status (*linear_support)(
		qaws_curve const* curve,
		unsigned int span_index,
		qaws_scalar local_t,
		unsigned int order,
		qaws_local_support* out);
};

struct qaws_surface_diff_vtable
{
	unsigned int capabilities;
	qaws_diff_class diff_class;

	unsigned int (*describe_fields)(
		qaws_surface const* surface,
		qaws_field_desc* out,
		unsigned int capacity);

	qaws_status (*primal_field)(
		qaws_surface const* surface,
		qaws_diff_field field,
		qaws_scalar const** out_data,
		unsigned int* out_count,
		unsigned int* out_components);

	qaws_status (*linear_support)(
		qaws_surface const* surface,
		qaws_scalar u,
		qaws_scalar v,
		unsigned int order,
		qaws_surface_local_support* out);
};

/* Fill one field descriptor. */
QAWS_INLINE qaws_field_desc qaws_internal_field_desc(
	qaws_diff_field field,
	qaws_value_type value_type,
	unsigned int count,
	qaws_value_domain domain,
	qaws_value_constraint constraint,
	qaws_diff_class diff_class,
	unsigned int capabilities)
{
	qaws_field_desc d;
	d.field = field;
	d.value_type = value_type;
	d.count = count;
	d.domain = domain;
	d.constraint = constraint;
	d.default_lower = -QAWS_DIFF_UNBOUNDED;
	d.default_upper = QAWS_DIFF_UNBOUNDED;
	if (constraint == QAWS_CONSTRAINT_POSITIVE || constraint == QAWS_CONSTRAINT_NON_NEGATIVE)
		d.default_lower = 0;
	if (constraint == QAWS_CONSTRAINT_UNIT_INTERVAL)
	{
		d.default_lower = 0;
		d.default_upper = 1;
	}
	d.diff_class = diff_class;
	d.capabilities = capabilities;
	return d;
}

/* B-spline basis derivatives for orders 0..k (k may exceed degree; higher
   rows are zero). out has (k+1) rows of stride (degree+1). */
void qaws_internal_bspline_basis_derivs_any(
	qaws_scalar const* knots,
	unsigned int knot_count,
	unsigned int degree,
	unsigned int span,
	qaws_scalar t,
	unsigned int k,
	qaws_scalar* out_ders);

/* Bernstein basis derivatives of degree n at t for orders 0..k.
   out[r * (n+1) + i] = d^r B_{i,n}(t) / dt^r. */
void qaws_internal_bernstein_derivs(
	unsigned int n,
	qaws_scalar t,
	unsigned int k,
	qaws_scalar* out);

#endif /* QAWS_INTERNAL_DIFF_H */
