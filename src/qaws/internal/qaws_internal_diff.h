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
		qaws_surface_support* out);

	/* Primal jet up to third order for non-linear families (NULL when the
	   jet follows from linear_support). */
	qaws_status (*eval_jet)(
		qaws_surface const* surface,
		qaws_scalar u,
		qaws_scalar v,
		unsigned int channels,
		qaws_surface_jet* out);
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

/* ------------------------------------------------------------------ */
/*  Shared adjoint accumulation engine                                */
/* ------------------------------------------------------------------ */

/* One parameter adjoint contribution of one sample. */
typedef struct qaws_diff_entry
{
	qaws_diff_field field;
	unsigned int element;
	qaws_scalar g[3];
} qaws_diff_entry;

/*
 * Produce the parameter contributions of `sample`. first_pass is non-zero
 * exactly once per sample per accumulate call: that is when coordinate
 * adjoints and report notes must be written.
 */
typedef qaws_status (*qaws_diff_collect_fn)(
	void const* user,
	unsigned int sample,
	int first_pass,
	qaws_diff_entry* entries,
	unsigned int capacity,
	unsigned int* out_count);

/* Accumulate (+=) every sample's contributions into views with the
   strategy selected by ctx (scatter, tiled or gather). */
qaws_status qaws_internal_diff_accumulate(
	qaws_diff_context const* ctx,
	unsigned int components,
	unsigned int sample_count,
	unsigned int entry_capacity,
	qaws_diff_collect_fn collect,
	void const* user,
	qaws_diff_views* views);

void qaws_internal_diff_report_note(
	qaws_diff_context const* ctx,
	qaws_diff_class diff_class,
	qaws_diff_validity validity,
	unsigned int frozen,
	unsigned int index);

unsigned int qaws_internal_view_stride(qaws_field_view const* v);
int qaws_internal_view_element_active(qaws_field_view const* v, unsigned int e);
int qaws_internal_view_component_active(qaws_field_view const* v, unsigned int c);

/* Masked read of one tangent element (zeros when inactive). */
void qaws_internal_view_read(qaws_field_view const* v, unsigned int e,
	unsigned int components, qaws_scalar* out);

/* Masked add into one adjoint element. */
void qaws_internal_view_add(qaws_field_view* v, unsigned int e,
	unsigned int components, qaws_scalar const* g);

/* Components of a field: 1 for scalar fields (weights, knots, radius...), else dim. */
unsigned int qaws_internal_field_components(qaws_diff_field field, unsigned int dim);

/* Every view with data must match qaws_internal_field_components. */
qaws_status qaws_internal_check_views(qaws_diff_views const* views, unsigned int components);

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
