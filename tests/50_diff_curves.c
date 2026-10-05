/*
 * Test 50: Differentiable curve evaluation
 *
 * Families (Bezier, B-spline, Hermite, polynomial, NURBS, rational Bezier):
 *   - primal jets match qaws_curve_evaluate_*
 *   - parameter tangents match finite differences
 *   - coordinate tangents match the next spatial derivative
 *   - adjoint identity over batches, all accumulation strategies agree
 *   - second tangent matches the finite difference of the tangent
 *   - activity masks, local support, support index, report, schema
 *   - composition with a dual kernel (unit tangent)
 */

#include "test_diff.h"

#define MAX_PARAMS 64
#define SAMPLE_COUNT 9

/* Weights and knots are scalar fields; every other field has one value per dimension. */
#define FIELD_COMPS(f, r) (((f)->fields[r] == QAWS_FIELD_WEIGHTS || (f)->fields[r] == QAWS_FIELD_KNOTS) ? 1u : (f)->dim)

/* ------------------------------------------------------------------ */
/*  Family fixtures                                                   */
/* ------------------------------------------------------------------ */

typedef struct family
{
	char const* name;
	unsigned int dim;
	unsigned int field_count;
	qaws_diff_field fields[3];
	unsigned int counts[3];
	qaws_scalar params[3][MAX_PARAMS];
	qaws_status (*create)(struct family const* f, qaws_scalar const* const* params, qaws_curve** out);
	qaws_scalar t_min, t_max;
	qaws_scalar knots[16];
	unsigned int knot_count;
	qaws_parameterization cr_param;
	int cr_closed;
	int nonlinear;              /* fields enter non-linearly (no QAWS_CAP_LINEAR) */
} family;

static qaws_status create_bezier(family const* f, qaws_scalar const* const* p, qaws_curve** out)
{
	qaws_bezier_desc d;
	d.dimension = (qaws_dimension)f->dim;
	d.degree = f->counts[0] - 1;
	d.control_points = p[0];
	d.control_point_count = f->counts[0];
	return qaws_curve_create_bezier(&d, out);
}

static qaws_status create_bspline(family const* f, qaws_scalar const* const* p, qaws_curve** out)
{
	qaws_bspline_desc d;
	memset(&d, 0, sizeof(d));
	d.dimension = (qaws_dimension)f->dim;
	d.degree = f->knot_count - f->counts[0] - 1;
	d.control_points = p[0];
	d.control_point_count = f->counts[0];
	d.knots = f->knots;
	d.knot_count = f->knot_count;
	return qaws_curve_create_bspline(&d, out);
}

static qaws_status create_bspline_knots(family const* f, qaws_scalar const* const* p, qaws_curve** out)
{
	qaws_bspline_desc d;
	memset(&d, 0, sizeof(d));
	d.dimension = (qaws_dimension)f->dim;
	d.degree = f->counts[1] - f->counts[0] - 1;
	d.control_points = p[0];
	d.control_point_count = f->counts[0];
	d.knots = p[1];
	d.knot_count = f->counts[1];
	return qaws_curve_create_bspline(&d, out);
}

static qaws_status create_hermite(family const* f, qaws_scalar const* const* p, qaws_curve** out)
{
	qaws_hermite_desc d;
	d.dimension = (qaws_dimension)f->dim;
	d.degree = 3;
	d.points = p[0];
	d.derivatives = p[1];
	d.point_count = f->counts[0];
	d.derivative_count = f->counts[1];
	return qaws_curve_create_hermite(&d, out);
}

static qaws_status create_polynomial(family const* f, qaws_scalar const* const* p, qaws_curve** out)
{
	qaws_polynomial_desc d;
	d.dimension = (qaws_dimension)f->dim;
	d.degree = f->counts[0] - 1;
	d.coefficients = p[0];
	d.coefficient_count = f->counts[0];
	d.t_min = f->t_min;
	d.t_max = f->t_max;
	return qaws_curve_create_polynomial(&d, out);
}

static qaws_status create_nurbs(family const* f, qaws_scalar const* const* p, qaws_curve** out)
{
	qaws_nurbs_desc d;
	memset(&d, 0, sizeof(d));
	d.dimension = (qaws_dimension)f->dim;
	d.degree = f->knot_count - f->counts[0] - 1;
	d.control_points = p[0];
	d.control_point_count = f->counts[0];
	d.knots = f->knots;
	d.knot_count = f->knot_count;
	d.weights = p[1];
	d.weight_count = f->counts[1];
	return qaws_curve_create_nurbs(&d, out);
}

static qaws_status create_nurbs_knots(family const* f, qaws_scalar const* const* p, qaws_curve** out)
{
	qaws_nurbs_desc d;
	memset(&d, 0, sizeof(d));
	d.dimension = (qaws_dimension)f->dim;
	d.degree = f->counts[2] - f->counts[0] - 1;
	d.control_points = p[0];
	d.control_point_count = f->counts[0];
	d.weights = p[1];
	d.weight_count = f->counts[1];
	d.knots = p[2];
	d.knot_count = f->counts[2];
	return qaws_curve_create_nurbs(&d, out);
}

static qaws_status create_catmull_rom(family const* f, qaws_scalar const* const* p, qaws_curve** out)
{
	qaws_catmull_rom_desc d;
	memset(&d, 0, sizeof(d));
	d.dimension = (qaws_dimension)f->dim;
	d.control_points = p[0];
	d.control_point_count = f->counts[0];
	d.parameterization = f->cr_param;
	d.closed = f->cr_closed;
	return qaws_curve_create_catmull_rom(&d, out);
}

static qaws_status create_yuksel(family const* f, qaws_scalar const* const* p, qaws_curve** out)
{
	qaws_yuksel_desc d;
	memset(&d, 0, sizeof(d));
	d.dimension = (qaws_dimension)f->dim;
	d.control_points = p[0];
	d.control_point_count = f->counts[0];
	d.mode = QAWS_YUKSEL_MODE_BEZIER;
	d.closed = f->cr_closed;
	return qaws_curve_create_yuksel(&d, out);
}

static qaws_status create_rational_bezier(family const* f, qaws_scalar const* const* p, qaws_curve** out)
{
	qaws_rational_bezier_desc d;
	d.dimension = (qaws_dimension)f->dim;
	d.degree = f->counts[0] - 1;
	d.control_points = p[0];
	d.control_point_count = f->counts[0];
	d.weights = p[1];
	d.weight_count = f->counts[1];
	return qaws_curve_create_rational_bezier(&d, out);
}

static void make_families(family* fams, unsigned int* count)
{
	unsigned int i;
	family* f;
	memset(fams, 0, sizeof(family) * 13);
	diff_seed(2024);

	f = &fams[0];
	f->name = "bezier";
	f->dim = 3;
	f->field_count = 1;
	f->fields[0] = QAWS_FIELD_CONTROL_POINTS;
	f->counts[0] = 6;
	diff_rand_fill(f->params[0], 6 * 3);
	f->create = create_bezier;
	f->t_min = 0; f->t_max = 1;

	f = &fams[1];
	f->name = "bspline";
	f->dim = 3;
	f->field_count = 1;
	f->fields[0] = QAWS_FIELD_CONTROL_POINTS;
	f->counts[0] = 7;
	diff_rand_fill(f->params[0], 7 * 3);
	f->create = create_bspline;
	{
		static qaws_scalar const k[11] = { 0, 0, 0, 0, 0.5f, 1.25f, 2, 3, 3, 3, 3 };
		for (i = 0; i < 11; i++) f->knots[i] = k[i];
		f->knot_count = 11;
	}
	f->t_min = 0; f->t_max = 3;

	f = &fams[2];
	f->name = "hermite";
	f->dim = 2;
	f->field_count = 2;
	f->fields[0] = QAWS_FIELD_POINTS;
	f->fields[1] = QAWS_FIELD_DERIVATIVES;
	f->counts[0] = 4;
	f->counts[1] = 4;
	diff_rand_fill(f->params[0], 4 * 2);
	diff_rand_fill(f->params[1], 4 * 2);
	f->create = create_hermite;
	f->t_min = 0; f->t_max = 3;

	f = &fams[3];
	f->name = "polynomial";
	f->dim = 3;
	f->field_count = 1;
	f->fields[0] = QAWS_FIELD_COEFFICIENTS;
	f->counts[0] = 5;
	diff_rand_fill(f->params[0], 5 * 3);
	f->create = create_polynomial;
	f->t_min = -1; f->t_max = 2;

	f = &fams[4];
	f->name = "nurbs";
	f->dim = 3;
	f->field_count = 2;
	f->fields[0] = QAWS_FIELD_CONTROL_POINTS;
	f->fields[1] = QAWS_FIELD_WEIGHTS;
	f->counts[0] = 7;
	f->counts[1] = 7;
	diff_rand_fill(f->params[0], 7 * 3);
	for (i = 0; i < 7; i++)
		f->params[1][i] = (qaws_scalar)1.25 + (qaws_scalar)0.75 * diff_rand();
	f->create = create_nurbs;
	{
		static qaws_scalar const k[11] = { 0, 0, 0, 0, 0.5f, 1.25f, 2, 3, 3, 3, 3 };
		for (i = 0; i < 11; i++) f->knots[i] = k[i];
		f->knot_count = 11;
	}
	f->t_min = 0; f->t_max = 3;

	f = &fams[5];
	f->name = "rational_bezier";
	f->dim = 2;
	f->field_count = 2;
	f->fields[0] = QAWS_FIELD_CONTROL_POINTS;
	f->fields[1] = QAWS_FIELD_WEIGHTS;
	f->counts[0] = 5;
	f->counts[1] = 5;
	diff_rand_fill(f->params[0], 5 * 2);
	for (i = 0; i < 5; i++)
		f->params[1][i] = (qaws_scalar)1.25 + (qaws_scalar)0.75 * diff_rand();
	f->create = create_rational_bezier;
	f->t_min = 0; f->t_max = 1;


	/* Knots as a differentiable field: unclamped so every knot can move
	   both ways, samples kept away from knots (d3 jumps there). */
	f = &fams[6];
	f->name = "bspline_knots";
	f->dim = 2;
	f->field_count = 2;
	f->fields[0] = QAWS_FIELD_CONTROL_POINTS;
	f->fields[1] = QAWS_FIELD_KNOTS;
	f->counts[0] = 7;
	f->counts[1] = 11;
	diff_rand_fill(f->params[0], 7 * 2);
	{
		static qaws_scalar const k[11] = { 0.1f, 0.6f, 1.0f, 1.5f, 2.13f, 2.58f, 3.03f, 3.5f, 3.9f, 4.4f, 4.8f };
		for (i = 0; i < 11; i++) f->params[1][i] = k[i];
	}
	f->create = create_bspline_knots;
	f->t_min = 1.5f; f->t_max = 3.5f;


	/* NURBS with weights and knots: the rational quotient carries the
	   knot terms (3D, unclamped, samples away from knots). */
	f = &fams[7];
	f->name = "nurbs_knots";
	f->dim = 3;
	f->field_count = 3;
	f->fields[0] = QAWS_FIELD_CONTROL_POINTS;
	f->fields[1] = QAWS_FIELD_WEIGHTS;
	f->fields[2] = QAWS_FIELD_KNOTS;
	f->counts[0] = 7;
	f->counts[1] = 7;
	f->counts[2] = 11;
	diff_rand_fill(f->params[0], 7 * 3);
	for (i = 0; i < 7; i++)
		f->params[1][i] = (qaws_scalar)1.25 + (qaws_scalar)0.75 * diff_rand();
	{
		static qaws_scalar const k[11] = { 0.1f, 0.6f, 1.0f, 1.5f, 2.13f, 2.58f, 3.03f, 3.5f, 3.9f, 4.4f, 4.8f };
		for (i = 0; i < 11; i++) f->params[2][i] = k[i];
	}
	f->create = create_nurbs_knots;
	f->t_min = 1.5f; f->t_max = 3.5f;


	/* Catmull-Rom: knot intervals |P_i+1 - P_i|^alpha make centripetal
	   and chordal curves non-linear in their points. */
	f = &fams[8];
	f->name = "catmull_rom_centripetal";
	f->dim = 2;
	f->field_count = 1;
	f->fields[0] = QAWS_FIELD_POINTS;
	f->counts[0] = 6;
	diff_rand_fill(f->params[0], 6 * 2);
	for (i = 0; i < 6; i++)
		f->params[0][2 * i] += (qaws_scalar)(1.5 * i);
	f->create = create_catmull_rom;
	f->cr_param = QAWS_PARAMETERIZATION_CENTRIPETAL;
	f->nonlinear = 1;
	f->t_min = 0; f->t_max = 3;

	f = &fams[9];
	f->name = "catmull_rom_chordal_closed";
	f->dim = 3;
	f->field_count = 1;
	f->fields[0] = QAWS_FIELD_POINTS;
	f->counts[0] = 5;
	for (i = 0; i < 5; i++)
	{
		f->params[0][3 * i + 0] = (qaws_scalar)(2.0 * cos(1.2566 * i)) + (qaws_scalar)0.2 * diff_rand();
		f->params[0][3 * i + 1] = (qaws_scalar)(2.0 * sin(1.2566 * i)) + (qaws_scalar)0.2 * diff_rand();
		f->params[0][3 * i + 2] = (qaws_scalar)0.5 * diff_rand();
	}
	f->create = create_catmull_rom;
	f->cr_param = QAWS_PARAMETERIZATION_CHORDAL;
	f->cr_closed = 1;
	f->nonlinear = 1;
	f->t_min = 0; f->t_max = 5;

	f = &fams[10];
	f->name = "catmull_rom_uniform";
	f->dim = 3;
	f->field_count = 1;
	f->fields[0] = QAWS_FIELD_POINTS;
	f->counts[0] = 6;
	diff_rand_fill(f->params[0], 6 * 3);
	f->create = create_catmull_rom;
	f->cr_param = QAWS_PARAMETERIZATION_UNIFORM;
	f->nonlinear = 1;
	f->t_min = 0; f->t_max = 3;


	/* Yuksel C2 interpolating splines (Bezier mode): the sub-curve
	   parameter is the root of a cubic in the points. */
	f = &fams[11];
	f->name = "yuksel_open";
	f->dim = 2;
	f->field_count = 1;
	f->fields[0] = QAWS_FIELD_POINTS;
	f->counts[0] = 6;
	diff_rand_fill(f->params[0], 6 * 2);
	for (i = 0; i < 6; i++)
		f->params[0][2 * i] += (qaws_scalar)(1.4 * i);
	f->create = create_yuksel;
	f->nonlinear = 1;
	f->t_min = 0; f->t_max = 5;

	f = &fams[12];
	f->name = "yuksel_closed";
	f->dim = 3;
	f->field_count = 1;
	f->fields[0] = QAWS_FIELD_POINTS;
	f->counts[0] = 5;
	for (i = 0; i < 5; i++)
	{
		f->params[0][3 * i + 0] = (qaws_scalar)(2.0 * cos(1.2566 * i)) + (qaws_scalar)0.3 * diff_rand();
		f->params[0][3 * i + 1] = (qaws_scalar)(2.0 * sin(1.2566 * i)) + (qaws_scalar)0.3 * diff_rand();
		f->params[0][3 * i + 2] = (qaws_scalar)0.5 * diff_rand();
	}
	f->create = create_yuksel;
	f->cr_closed = 1;
	f->nonlinear = 1;
	f->t_min = 0; f->t_max = 5;

	*count = 13;
}

static qaws_curve* family_curve(family const* f, qaws_scalar const (*params)[MAX_PARAMS])
{
	qaws_scalar const* p[3];
	qaws_curve* c = NULL;
	p[0] = params[0];
	p[1] = params[1];
	p[2] = params[2];
	if (f->create(f, p, &c) != QAWS_STATUS_OK)
		return NULL;
	return c;
}

/* Curve built from params + h * dir. */
static qaws_curve* family_curve_shifted(family const* f, qaws_scalar const (*dir)[MAX_PARAMS], double h)
{
	qaws_scalar p[3][MAX_PARAMS];
	unsigned int r, i;
	for (r = 0; r < f->field_count; r++)
		for (i = 0; i < f->counts[r] * FIELD_COMPS(f, r); i++)
			p[r][i] = (qaws_scalar)(f->params[r][i] + h * dir[r][i]);
	return family_curve(f, (qaws_scalar const (*)[MAX_PARAMS])p);
}

static void family_views(family const* f, qaws_scalar (*storage)[MAX_PARAMS],
	qaws_field_view* views_storage, qaws_diff_views* views)
{
	unsigned int r;
	for (r = 0; r < f->field_count; r++)
		views_storage[r] = qaws_field_view_make(f->fields[r], storage[r], f->counts[r], FIELD_COMPS(f, r));
	views->fields = views_storage;
	views->field_count = f->field_count;
	views->children = NULL;
	views->child_count = 0;
}

static qaws_scalar sample_t(family const* f, unsigned int i)
{
	return f->t_min + (f->t_max - f->t_min) * ((qaws_scalar)i + (qaws_scalar)0.37) / (qaws_scalar)SAMPLE_COUNT;
}

/* Evaluate any family into a 3D jet (z = 0 for 2D). */
static qaws_status eval_jet(family const* f, qaws_curve const* c, qaws_scalar t,
	qaws_scalar t_dot, qaws_diff_views const* tangent_views,
	qaws_curve_jet_3d* primal, qaws_curve_jet_3d* tangent, qaws_curve_jet_3d* tangent2)
{
	qaws_status st;
	unsigned int k;
	if (f->dim == 3)
	{
		if (tangent2)
			return qaws_curve_eval_batch_tangent2_3d(NULL, c, &t, &t_dot, 1, 0xF, tangent_views, primal, tangent, tangent2);
		return qaws_curve_eval_tangent_3d(NULL, c, t, t_dot, 0xF, tangent_views, primal, tangent);
	}
	else
	{
		qaws_curve_jet_2d p2, t2, tt2;
		if (tangent2)
			st = qaws_curve_eval_batch_tangent2_2d(NULL, c, &t, &t_dot, 1, 0xF, tangent_views, &p2, &t2, &tt2);
		else
			st = qaws_curve_eval_tangent_2d(NULL, c, t, t_dot, 0xF, tangent_views, &p2, &t2);
		for (k = 0; k < 4; k++)
		{
			if (primal) { primal->d[k].x = p2.d[k].x; primal->d[k].y = p2.d[k].y; primal->d[k].z = 0; }
			if (tangent) { tangent->d[k].x = t2.d[k].x; tangent->d[k].y = t2.d[k].y; tangent->d[k].z = 0; }
			if (tangent2) { tangent2->d[k].x = tt2.d[k].x; tangent2->d[k].y = tt2.d[k].y; tangent2->d[k].z = 0; }
		}
		return st;
	}
}

static int jet_close(qaws_curve_jet_3d const* a, qaws_curve_jet_3d const* b, double tol)
{
	unsigned int k;
	for (k = 0; k < 4; k++)
		if (!diff_close(a->d[k].x, b->d[k].x, tol) || !diff_close(a->d[k].y, b->d[k].y, tol) ||
		    !diff_close(a->d[k].z, b->d[k].z, tol))
			return 0;
	return 1;
}

/* ------------------------------------------------------------------ */
/*  Per family checks                                                 */
/* ------------------------------------------------------------------ */

static void check_primal(family const* f, qaws_curve const* c)
{
	unsigned int i;
	int ok = 1;
	for (i = 0; i < SAMPLE_COUNT; i++)
	{
		qaws_scalar t = sample_t(f, i);
		qaws_curve_jet_3d p, tg;
		qaws_curve_jet_3d ref;
		if (eval_jet(f, c, t, 0, NULL, &p, &tg, NULL) != QAWS_STATUS_OK) { ok = 0; break; }
		memset(&ref, 0, sizeof(ref));
		if (f->dim == 3)
		{
			qaws_eval_result_3d r;
			qaws_curve_evaluate_3d(c, t, 0xF, &r);
			ref.d[0] = r.position; ref.d[1] = r.d1; ref.d[2] = r.d2; ref.d[3] = r.d3;
		}
		else
		{
			qaws_eval_result_2d r;
			qaws_curve_evaluate_2d(c, t, 0xF, &r);
			ref.d[0].x = r.position.x; ref.d[0].y = r.position.y;
			ref.d[1].x = r.d1.x; ref.d[1].y = r.d1.y;
			ref.d[2].x = r.d2.x; ref.d[2].y = r.d2.y;
			ref.d[3].x = r.d3.x; ref.d[3].y = r.d3.y;
		}
		if (!jet_close(&p, &ref, 1e-4))
			ok = 0;
	}
	printf("    %s: primal jets match evaluate: %s\n", f->name, ok ? "yes" : "NO");
	TEST_ASSERT(ok, "primal jet matches qaws_curve_evaluate");
}

static void check_param_tangent_fd(family const* f, qaws_curve const* c)
{
	qaws_scalar dir[3][MAX_PARAMS];
	qaws_field_view vs[3];
	qaws_diff_views views;
	unsigned int r, i, k;
	int ok = 1;
	double h = DIFF_FD_STEP;
	qaws_curve *cp, *cm;

	for (r = 0; r < f->field_count; r++)
		diff_rand_fill(dir[r], f->counts[r] * FIELD_COMPS(f, r));
	family_views(f, dir, vs, &views);

	cp = family_curve_shifted(f, (qaws_scalar const (*)[MAX_PARAMS])dir, h);
	cm = family_curve_shifted(f, (qaws_scalar const (*)[MAX_PARAMS])dir, -h);
	for (i = 0; i < SAMPLE_COUNT && cp && cm; i++)
	{
		qaws_scalar t = sample_t(f, i);
		qaws_curve_jet_3d p, tg, pp, pm, dummy;
		eval_jet(f, c, t, 0, &views, &p, &tg, NULL);
		eval_jet(f, cp, t, 0, NULL, &pp, &dummy, NULL);
		eval_jet(f, cm, t, 0, NULL, &pm, &dummy, NULL);
		for (k = 0; k < 4; k++)
		{
			double fx = (pp.d[k].x - pm.d[k].x) / (2 * h);
			double fy = (pp.d[k].y - pm.d[k].y) / (2 * h);
			double fz = (pp.d[k].z - pm.d[k].z) / (2 * h);
			if (!diff_close(fx, tg.d[k].x, DIFF_TOL * 50) || !diff_close(fy, tg.d[k].y, DIFF_TOL * 50) ||
			    !diff_close(fz, tg.d[k].z, DIFF_TOL * 50))
				ok = 0;
		}
	}
	if (cp) qaws_curve_destroy(cp);
	if (cm) qaws_curve_destroy(cm);
	printf("    %s: parameter tangent matches finite differences: %s\n", f->name, ok ? "yes" : "NO");
	TEST_ASSERT(ok, "parameter tangent matches finite differences");
}

static void check_coordinate_tangent(family const* f, qaws_curve const* c)
{
	unsigned int i, k;
	int ok = 1;
	for (i = 0; i < SAMPLE_COUNT; i++)
	{
		qaws_scalar t = sample_t(f, i);
		qaws_curve_jet_3d p, tg;
		eval_jet(f, c, t, 1, NULL, &p, &tg, NULL);
		for (k = 0; k < 3; k++)
			if (!diff_close(tg.d[k].x, p.d[k + 1].x, 1e-4) || !diff_close(tg.d[k].y, p.d[k + 1].y, 1e-4) ||
			    !diff_close(tg.d[k].z, p.d[k + 1].z, 1e-4))
				ok = 0;
	}
	TEST_ASSERT(ok, "coordinate tangent equals next spatial derivative");
}

/* Run the batch adjoint for one family with a given accumulation strategy. */
static qaws_status run_adjoint(family const* f, qaws_curve const* c, qaws_diff_accumulation acc,
	unsigned int tile, qaws_scalar const* ts, void const* ybar,
	qaws_scalar (*pbar)[MAX_PARAMS], qaws_scalar* tbar)
{
	qaws_diff_context ctx;
	qaws_field_view vs[3];
	qaws_diff_views views;
	unsigned int r;

	qaws_diff_context_init(&ctx);
	ctx.accumulation = acc;
	ctx.tile_size = tile;
	for (r = 0; r < 3; r++)
		memset(pbar[r], 0, sizeof(qaws_scalar) * MAX_PARAMS);
	memset(tbar, 0, sizeof(qaws_scalar) * SAMPLE_COUNT);
	family_views(f, pbar, vs, &views);
	if (f->dim == 3)
		return qaws_curve_eval_batch_adjoint_3d(&ctx, c, ts, SAMPLE_COUNT, 0xF,
			(qaws_curve_jet_3d const*)ybar, &views, tbar);
	return qaws_curve_eval_batch_adjoint_2d(&ctx, c, ts, SAMPLE_COUNT, 0xF,
		(qaws_curve_jet_2d const*)ybar, &views, tbar);
}

static void check_adjoint_identity(family const* f, qaws_curve const* c)
{
	qaws_scalar ts[SAMPLE_COUNT], tdots[SAMPLE_COUNT], tbar[SAMPLE_COUNT];
	qaws_scalar dir[3][MAX_PARAMS], pbar[3][MAX_PARAMS];
	qaws_curve_jet_3d ybar3[SAMPLE_COUNT], tan3[SAMPLE_COUNT];
	qaws_curve_jet_2d ybar2[SAMPLE_COUNT], tan2[SAMPLE_COUNT];
	qaws_field_view vs[3];
	qaws_diff_views views;
	unsigned int i, r;
	double lhs = 0, rhs = 0;
	qaws_status st;

	memset(dir, 0, sizeof(dir));
	for (i = 0; i < SAMPLE_COUNT; i++)
	{
		ts[i] = sample_t(f, i);
		tdots[i] = diff_rand();
		diff_rand_jet3(&ybar3[i]);
		diff_rand_jet2(&ybar2[i]);
	}
	for (r = 0; r < f->field_count; r++)
		diff_rand_fill(dir[r], f->counts[r] * FIELD_COMPS(f, r));
	family_views(f, dir, vs, &views);

	if (f->dim == 3)
		st = qaws_curve_eval_batch_tangent_3d(NULL, c, ts, tdots, SAMPLE_COUNT, 0xF, &views, NULL, tan3);
	else
		st = qaws_curve_eval_batch_tangent_2d(NULL, c, ts, tdots, SAMPLE_COUNT, 0xF, &views, NULL, tan2);
	TEST_ASSERT_STATUS(st);

	for (i = 0; i < SAMPLE_COUNT; i++)
		lhs += (f->dim == 3) ? diff_jet3_dot(&ybar3[i], &tan3[i], 0xF) : diff_jet2_dot(&ybar2[i], &tan2[i], 0xF);

	st = run_adjoint(f, c, QAWS_ACCUMULATE_SCATTER, 0, ts,
		(f->dim == 3) ? (void const*)ybar3 : (void const*)ybar2, pbar, tbar);
	TEST_ASSERT_STATUS(st);
	for (r = 0; r < f->field_count; r++)
		rhs += diff_dot(pbar[r], dir[r], f->counts[r] * FIELD_COMPS(f, r));
	rhs += diff_dot(tbar, tdots, SAMPLE_COUNT);

	printf("    %s: <ybar, J xdot> = %.9g  <J^T ybar, xdot> = %.9g\n", f->name, lhs, rhs);
	TEST_ASSERT(diff_close(lhs, rhs, DIFF_TOL), "adjoint identity over a batch");

	/* All accumulation strategies agree. */
	{
		qaws_scalar pbar_t[3][MAX_PARAMS], pbar_g[3][MAX_PARAMS];
		qaws_scalar tbar_t[SAMPLE_COUNT], tbar_g[SAMPLE_COUNT];
		int ok = 1;
		unsigned int n;
		st = run_adjoint(f, c, QAWS_ACCUMULATE_TILED, 2, ts,
			(f->dim == 3) ? (void const*)ybar3 : (void const*)ybar2, pbar_t, tbar_t);
		TEST_ASSERT_STATUS(st);
		st = run_adjoint(f, c, QAWS_ACCUMULATE_GATHER, 0, ts,
			(f->dim == 3) ? (void const*)ybar3 : (void const*)ybar2, pbar_g, tbar_g);
		TEST_ASSERT_STATUS(st);
		for (r = 0; r < f->field_count; r++)
			for (n = 0; n < f->counts[r] * FIELD_COMPS(f, r); n++)
				if (!diff_close(pbar[r][n], pbar_t[r][n], DIFF_TOL) || !diff_close(pbar[r][n], pbar_g[r][n], DIFF_TOL))
					ok = 0;
		for (n = 0; n < SAMPLE_COUNT; n++)
			if (!diff_close(tbar[n], tbar_t[n], DIFF_TOL) || !diff_close(tbar[n], tbar_g[n], DIFF_TOL))
				ok = 0;
		TEST_ASSERT(ok, "scatter, tiled and gather accumulation agree");
	}
}

static void check_tangent2(family const* f, qaws_curve const* c)
{
	qaws_scalar dir[3][MAX_PARAMS];
	qaws_field_view vs[3];
	qaws_diff_views views;
	unsigned int r, i, k;
	int ok = 1;
	double h = DIFF_FD_STEP;

	memset(dir, 0, sizeof(dir));
	for (r = 0; r < f->field_count; r++)
		diff_rand_fill(dir[r], f->counts[r] * FIELD_COMPS(f, r));
	family_views(f, dir, vs, &views);

	for (i = 0; i < SAMPLE_COUNT; i++)
	{
		qaws_scalar t = sample_t(f, i);
		qaws_scalar tdot = (qaws_scalar)0.7;
		qaws_curve_jet_3d p, tg, tt, pp, tgp, pm, tgm;
		qaws_curve *cp, *cm;
		if (eval_jet(f, c, t, tdot, &views, &p, &tg, &tt) != QAWS_STATUS_OK) { ok = 0; break; }
		cp = family_curve_shifted(f, (qaws_scalar const (*)[MAX_PARAMS])dir, h);
		cm = family_curve_shifted(f, (qaws_scalar const (*)[MAX_PARAMS])dir, -h);
		eval_jet(f, cp, (qaws_scalar)(t + h * tdot), tdot, &views, &pp, &tgp, NULL);
		eval_jet(f, cm, (qaws_scalar)(t - h * tdot), tdot, &views, &pm, &tgm, NULL);
		for (k = 0; k < 4; k++)
		{
			double fx = (tgp.d[k].x - tgm.d[k].x) / (2 * h);
			double fy = (tgp.d[k].y - tgm.d[k].y) / (2 * h);
			double fz = (tgp.d[k].z - tgm.d[k].z) / (2 * h);
			/* relative to the vector: small components of large jets carry its rounding */
			double scale = fabs(tt.d[k].x) + fabs(tt.d[k].y) + fabs(tt.d[k].z);
			if (!diff_close(fx / (1 + scale), tt.d[k].x / (1 + scale), DIFF_TOL * 50) ||
			    !diff_close(fy / (1 + scale), tt.d[k].y / (1 + scale), DIFF_TOL * 50) ||
			    !diff_close(fz / (1 + scale), tt.d[k].z / (1 + scale), DIFF_TOL * 50))
				ok = 0;
		}
		qaws_curve_destroy(cp);
		qaws_curve_destroy(cm);
	}
	printf("    %s: second tangent matches finite differences: %s\n", f->name, ok ? "yes" : "NO");
	TEST_ASSERT(ok, "second tangent matches finite difference of the tangent");
}

static void test_curve_families(void)
{
	family fams[13];
	unsigned int n, i;
	make_families(fams, &n);
	for (i = 0; i < n; i++)
	{
		qaws_curve* c = family_curve(&fams[i], (qaws_scalar const (*)[MAX_PARAMS])fams[i].params);
		TEST_ASSERT(c != NULL, "family curve created");
		if (!c)
			continue;
		{
			unsigned int caps = qaws_curve_get_diff_capabilities(c);
			int rational = fams[i].field_count > 1 && fams[i].fields[1] == QAWS_FIELD_WEIGHTS;
			TEST_ASSERT((caps & (QAWS_CAP_TANGENT | QAWS_CAP_ADJOINT | QAWS_CAP_TANGENT2)) ==
				(QAWS_CAP_TANGENT | QAWS_CAP_ADJOINT | QAWS_CAP_TANGENT2), "family capabilities");
			TEST_ASSERT(((caps & QAWS_CAP_LINEAR) != 0) == !(rational || fams[i].nonlinear), "linear flag only for linear families");
		}
		check_primal(&fams[i], c);
		check_param_tangent_fd(&fams[i], c);
		check_coordinate_tangent(&fams[i], c);
		check_adjoint_identity(&fams[i], c);
		check_tangent2(&fams[i], c);
		qaws_curve_destroy(c);
	}
}

/* ------------------------------------------------------------------ */
/*  Masks, support, report, schema, composition                       */
/* ------------------------------------------------------------------ */

static qaws_curve* make_bspline(void)
{
	static qaws_scalar const cps[7 * 3] = {
		0, 0, 0,  1, 2, 0,  2, -1, 1,  3, 1, 0,  4, 0, -1,  5, 2, 1,  6, 0, 0 };
	static qaws_scalar const knots[11] = { 0, 0, 0, 0, 0.5f, 1.25f, 2, 3, 3, 3, 3 };
	qaws_bspline_desc d;
	qaws_curve* c = NULL;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_3D;
	d.degree = 3;
	d.control_points = cps;
	d.control_point_count = 7;
	d.knots = knots;
	d.knot_count = 11;
	qaws_curve_create_bspline(&d, &c);
	return c;
}

static void test_masks(void)
{
	qaws_curve* c = make_bspline();
	qaws_scalar dir[21], dir_zeroed[21], pbar[21];
	unsigned char active[7] = { 1, 1, 0, 1, 1, 1, 1 };
	qaws_field_view v;
	qaws_diff_views views;
	qaws_curve_jet_3d p, tg_masked, tg_ref, ybar;
	qaws_scalar t = (qaws_scalar)0.8;
	unsigned int i;
	int untouched = 1;

	diff_seed(5);
	diff_rand_fill(dir, 21);
	memcpy(dir_zeroed, dir, sizeof(dir));
	dir_zeroed[6] = dir_zeroed[7] = dir_zeroed[8] = 0;

	views.fields = &v;
	views.field_count = 1;
	views.children = NULL;
	views.child_count = 0;

	v = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, dir, 7, 3);
	v.active = active;
	qaws_curve_eval_tangent_3d(NULL, c, t, 0, 0xF, &views, &p, &tg_masked);

	v = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, dir_zeroed, 7, 3);
	qaws_curve_eval_tangent_3d(NULL, c, t, 0, 0xF, &views, &p, &tg_ref);
	TEST_ASSERT(jet_close(&tg_masked, &tg_ref, 1e-6), "inactive element contributes no tangent");

	memset(pbar, 0, sizeof(pbar));
	v = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, pbar, 7, 3);
	v.active = active;
	v.component_mask = 1u | 4u; /* x and z only */
	diff_rand_jet3(&ybar);
	qaws_curve_eval_adjoint_3d(NULL, c, t, 0xF, &ybar, &views, NULL);
	for (i = 0; i < 3; i++)
		if (pbar[6 + i] != 0)
			untouched = 0;
	for (i = 0; i < 7; i++)
		if (pbar[i * 3 + 1] != 0)
			untouched = 0;
	TEST_ASSERT(untouched, "inactive elements and components receive no adjoint");
	TEST_ASSERT(pbar[0] != 0 || pbar[3] != 0, "active elements receive adjoint");

	qaws_curve_destroy(c);
}

static void test_support(void)
{
	qaws_curve* c = make_bspline();
	qaws_local_support s;
	qaws_scalar ts[5] = { 0.1f, 0.6f, 1.0f, 2.2f, 2.9f };
	unsigned int offsets[8], samples[64], entries = 0, i, j, e;
	int partition = 1, consistent = 1;

	for (i = 0; i < 5; i++)
	{
		qaws_scalar sum = 0;
		TEST_ASSERT_STATUS(qaws_curve_local_support(c, ts[i], 1, &s));
		for (j = 0; j < s.ranges[0].count; j++)
			sum += s.weights[0][0][j];
		if (!approx_eq_loose(sum, 1))
			partition = 0;
	}
	TEST_ASSERT(partition, "basis weights form a partition of unity");
	TEST_ASSERT(s.kind == QAWS_SUPPORT_LOCAL && s.ranges[0].count == 4, "cubic B-spline support is 4 local points");

	TEST_ASSERT(qaws_curve_build_support_index(c, ts, 5, QAWS_FIELD_CONTROL_POINTS,
		offsets, 8, NULL, 0, &entries) == QAWS_STATUS_BUFFER_TOO_SMALL, "size query");
	TEST_ASSERT(entries == 20, "support index entry count");
	TEST_ASSERT_STATUS(qaws_curve_build_support_index(c, ts, 5, QAWS_FIELD_CONTROL_POINTS,
		offsets, 8, samples, 64, &entries));

	for (i = 0; i < 5; i++)
	{
		qaws_curve_local_support(c, ts[i], 0, &s);
		for (j = 0; j < s.ranges[0].count; j++)
		{
			int found = 0;
			e = s.ranges[0].first + j;
			for (unsigned int q = offsets[e]; q < offsets[e + 1]; q++)
				if (samples[q] == i)
					found = 1;
			if (!found)
				consistent = 0;
		}
	}
	TEST_ASSERT(consistent, "support index lists every influenced sample");
	qaws_curve_destroy(c);
}

static void test_report(void)
{
	qaws_curve* c = make_bspline();
	qaws_diff_context ctx;
	qaws_diff_report report;
	qaws_curve_jet_3d p, tg;

	qaws_diff_context_init(&ctx);
	qaws_diff_report_reset(&report);
	ctx.report = &report;

	qaws_curve_eval_tangent_3d(&ctx, c, (qaws_scalar)0.8, 0, 0xF, NULL, &p, &tg);
	TEST_ASSERT(report.validity == QAWS_DIFF_VALID, "parameter-only tangent is valid");

	qaws_curve_eval_tangent_3d(&ctx, c, (qaws_scalar)1.25, 1, 0xF, NULL, &p, &tg);
	TEST_ASSERT(report.validity == QAWS_DIFF_AT_BOUNDARY, "coordinate tangent at a knot is one-sided");
	TEST_ASSERT(report.worst_index == 0, "worst index recorded");
	TEST_ASSERT((report.frozen_used & QAWS_FREEZE_SPAN) != 0, "span selection recorded as frozen state");
	TEST_ASSERT(report.diff_class == QAWS_DIFF_PIECEWISE_SMOOTH, "B-spline class is piecewise smooth");
	TEST_ASSERT(report.evaluation_count == 2, "evaluation count");
	qaws_curve_destroy(c);
}

static void test_schema(void)
{
	qaws_curve* c = make_bspline();
	qaws_field_desc fields[4];
	unsigned int n = 0;
	qaws_scalar cps[21];
	unsigned int got = 0;

	TEST_ASSERT_STATUS(qaws_curve_describe_fields(c, fields, 4, &n));
	TEST_ASSERT(n == 2, "B-spline exposes control points and knots");
	TEST_ASSERT(fields[0].field == QAWS_FIELD_CONTROL_POINTS && fields[0].domain == QAWS_DOMAIN_POSITION &&
		fields[0].count == 7 && fields[0].value_type == QAWS_VALUE_VEC3, "control point schema");
	TEST_ASSERT(fields[1].field == QAWS_FIELD_KNOTS && fields[1].constraint == QAWS_CONSTRAINT_MONOTONIC &&
		(fields[1].capabilities & QAWS_CAP_ADJOINT) && !(fields[1].capabilities & QAWS_CAP_LINEAR),
		"knots are monotonic and differentiable (non-linear)");
	TEST_ASSERT(qaws_curve_describe_fields(c, fields, 1, &n) == QAWS_STATUS_BUFFER_TOO_SMALL && n == 2,
		"describe reports required capacity");

	TEST_ASSERT_STATUS(qaws_curve_read_field(c, QAWS_FIELD_CONTROL_POINTS, cps, 21, &got));
	TEST_ASSERT(got == 21 && cps[3] == 1 && cps[4] == 2, "read_field returns primal values");
	TEST_ASSERT(qaws_curve_get_coordinate_kind(c) == QAWS_COORDINATE_PARAMETRIC, "parametric coordinate");
	qaws_curve_destroy(c);

	{
		/* Not yet differentiable families report it explicitly. */
		qaws_clothoid_desc d;
		qaws_curve* cl = NULL;
		qaws_curve_jet_2d p, tg;
		memset(&d, 0, sizeof(d));
		d.start_curvature = (qaws_scalar)0.1;
		d.end_curvature = (qaws_scalar)0.8;
		d.length = (qaws_scalar)2.0;
		qaws_curve_create_clothoid(&d, &cl);
		TEST_ASSERT(qaws_curve_get_diff_capabilities(cl) == 0, "no capabilities without rules");
		TEST_ASSERT(qaws_curve_get_diff_class(cl) == QAWS_DIFF_UNSUPPORTED, "unsupported class");
		TEST_ASSERT(qaws_curve_eval_tangent_2d(NULL, cl, (qaws_scalar)0.5, 0, 1, NULL, &p, &tg)
			== QAWS_STATUS_UNSUPPORTED_OPERATION, "tangent refused, never approximated");
		qaws_curve_destroy(cl);
	}
}

/* Unit tangent T = normalize(C'(t)): curve rule composed with a dual kernel. */
static void test_unit_tangent_composition(void)
{
	qaws_curve* c = make_bspline();
	qaws_scalar dir[21], pbar[21];
	qaws_field_view v;
	qaws_diff_views views;
	unsigned int i;
	int ok = 1;

	views.fields = &v;
	views.field_count = 1;
	views.children = NULL;
	views.child_count = 0;
	diff_seed(77);

	for (i = 0; i < 6; i++)
	{
		qaws_scalar t = (qaws_scalar)0.2 + (qaws_scalar)0.45 * (qaws_scalar)i;
		qaws_scalar tdot = diff_rand(), tbar = 0;
		qaws_curve_jet_3d p, tg, ybar;
		qaws_dual3 T;
		qaws_vec3 Tbar = diff_rand_vec3();
		double lhs, rhs;

		diff_rand_fill(dir, 21);
		v = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, dir, 7, 3);
		qaws_curve_eval_tangent_3d(NULL, c, t, tdot, QAWS_EVAL_FLAG_D1, &views, &p, &tg);
		T = qaws_dual3_normalize(qaws_dual3_make(p.d[1], tg.d[1], qaws_v3_zero()));
		lhs = qaws_v3_dot(Tbar, T.t);

		memset(&ybar, 0, sizeof(ybar));
		ybar.d[1] = qaws_normalize_adjoint(p.d[1], Tbar);
		memset(pbar, 0, sizeof(pbar));
		v = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, pbar, 7, 3);
		qaws_curve_eval_adjoint_3d(NULL, c, t, QAWS_EVAL_FLAG_D1, &ybar, &views, &tbar);
		rhs = diff_dot(pbar, dir, 21) + (double)tbar * tdot;

		if (!diff_close(lhs, rhs, DIFF_TOL))
			ok = 0;
	}
	TEST_ASSERT(ok, "unit tangent composition satisfies the adjoint identity");
	qaws_curve_destroy(c);
}

/* ------------------------------------------------------------------ */
/*  Composite curves: rules chain into the segments (children)        */
/* ------------------------------------------------------------------ */

typedef struct composite_params
{
	qaws_scalar bez[12];   /* cubic Bezier, 4 points */
	qaws_scalar nrb[21];   /* NURBS, 7 points */
	qaws_scalar w[7];
} composite_params;

static qaws_scalar const g_comp_knots[11] = { 0, 0, 0, 0, 0.5f, 1.25f, 2, 3, 3, 3, 3 };

static qaws_curve* composite_build(composite_params const* p)
{
	qaws_bezier_desc bd;
	qaws_nurbs_desc nd;
	qaws_composite_desc cd;
	qaws_curve* segs[2] = { NULL, NULL };
	qaws_curve* c = NULL;

	bd.dimension = QAWS_DIMENSION_3D;
	bd.degree = 3;
	bd.control_points = p->bez;
	bd.control_point_count = 4;
	qaws_curve_create_bezier(&bd, &segs[0]);

	memset(&nd, 0, sizeof(nd));
	nd.dimension = QAWS_DIMENSION_3D;
	nd.degree = 3;
	nd.control_points = p->nrb;
	nd.control_point_count = 7;
	nd.knots = g_comp_knots;
	nd.knot_count = 11;
	nd.weights = p->w;
	nd.weight_count = 7;
	qaws_curve_create_nurbs(&nd, &segs[1]);

	cd.dimension = QAWS_DIMENSION_3D;
	cd.segments = segs;
	cd.segment_count = 2;
	qaws_curve_create_composite(&cd, &c);
	return c;
}

typedef struct composite_storage
{
	composite_params data;
	qaws_field_view v0, v1[2];
	qaws_diff_views child[2], views;
} composite_storage;

static void composite_bind(composite_storage* s)
{
	s->v0 = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, s->data.bez, 4, 3);
	s->v1[0] = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, s->data.nrb, 7, 3);
	s->v1[1] = qaws_field_view_make(QAWS_FIELD_WEIGHTS, s->data.w, 7, 1);
	s->child[0].fields = &s->v0;
	s->child[0].field_count = 1;
	s->child[0].children = NULL;
	s->child[0].child_count = 0;
	s->child[1].fields = s->v1;
	s->child[1].field_count = 2;
	s->child[1].children = NULL;
	s->child[1].child_count = 0;
	s->views.fields = NULL;
	s->views.field_count = 0;
	s->views.children = s->child;
	s->views.child_count = 2;
}

static double composite_dot(composite_params const* a, composite_params const* b)
{
	return diff_dot(a->bez, b->bez, 12) + diff_dot(a->nrb, b->nrb, 21) + diff_dot(a->w, b->w, 7);
}

static void test_composite(void)
{
	composite_params base, pp, pm;
	qaws_curve* c;
	unsigned int i, n, k;
	int ok_jet = 1, ok_fd = 1, ok_t2 = 1, ok_adj = 1;
	double h = DIFF_FD_STEP;
	qaws_diff_child kids[2];
	unsigned int nk = 0;

	diff_seed(321);
	diff_rand_fill(base.bez, 12);
	diff_rand_fill(base.nrb, 21);
	for (n = 0; n < 7; n++)
		base.w[n] = (qaws_scalar)1.2 + (qaws_scalar)0.5 * diff_rand();
	c = composite_build(&base);
	TEST_ASSERT(c != NULL, "composite built");
	TEST_ASSERT_STATUS(qaws_curve_diff_children(c, kids, 2, &nk));
	TEST_ASSERT(nk == 2 && kids[0].curve && kids[1].curve, "composite exposes its segments as children");

	for (i = 0; i < 8; i++)
	{
		qaws_scalar t = (qaws_scalar)0.13 + (qaws_scalar)0.23 * (qaws_scalar)i;
		qaws_scalar tdot = diff_rand(), tbar = 0;
		composite_storage dir, bar;
		qaws_curve_jet_3d p, tg, tt, ybar;
		qaws_eval_result_3d r;
		qaws_curve *cp, *cm;
		double lhs = 0, rhs;

		qaws_curve_evaluate_3d(c, t, 0xF, &r);
		memset(&dir, 0, sizeof(dir));
		diff_rand_fill(dir.data.bez, 12);
		diff_rand_fill(dir.data.nrb, 21);
		diff_rand_fill(dir.data.w, 7);
		composite_bind(&dir);
		TEST_ASSERT_STATUS(qaws_curve_eval_batch_tangent2_3d(NULL, c, &t, &tdot, 1, 0xF, &dir.views, &p, &tg, &tt));
		if (!diff_close(p.d[0].x, r.position.x, 1e-5) || !diff_close(p.d[1].y, r.d1.y, 1e-5) ||
		    !diff_close(p.d[2].z, r.d2.z, 1e-5) || !diff_close(p.d[3].x, r.d3.x, 1e-4))
			ok_jet = 0;

		for (n = 0; n < 12; n++) { pp.bez[n] = (qaws_scalar)(base.bez[n] + h * dir.data.bez[n]); pm.bez[n] = (qaws_scalar)(base.bez[n] - h * dir.data.bez[n]); }
		for (n = 0; n < 21; n++) { pp.nrb[n] = (qaws_scalar)(base.nrb[n] + h * dir.data.nrb[n]); pm.nrb[n] = (qaws_scalar)(base.nrb[n] - h * dir.data.nrb[n]); }
		for (n = 0; n < 7; n++) { pp.w[n] = (qaws_scalar)(base.w[n] + h * dir.data.w[n]); pm.w[n] = (qaws_scalar)(base.w[n] - h * dir.data.w[n]); }
		cp = composite_build(&pp);
		cm = composite_build(&pm);
		{
			qaws_curve_jet_3d jp, jm, tp, tm;
			qaws_curve_eval_tangent_3d(NULL, cp, (qaws_scalar)(t + h * tdot), tdot, 0xF, &dir.views, &jp, &tp);
			qaws_curve_eval_tangent_3d(NULL, cm, (qaws_scalar)(t - h * tdot), tdot, 0xF, &dir.views, &jm, &tm);
			for (k = 0; k < 4; k++)
			{
				if (!diff_close((jp.d[k].x - jm.d[k].x) / (2 * h), tg.d[k].x, DIFF_TOL * 100) ||
				    !diff_close((jp.d[k].z - jm.d[k].z) / (2 * h), tg.d[k].z, DIFF_TOL * 100))
					ok_fd = 0;
				if (!diff_close((tp.d[k].y - tm.d[k].y) / (2 * h), tt.d[k].y, DIFF_TOL * 300))
					ok_t2 = 0;
			}
		}
		qaws_curve_destroy(cp);
		qaws_curve_destroy(cm);

		diff_rand_jet3(&ybar);
		lhs = diff_jet3_dot(&ybar, &tg, 0xF);
		memset(&bar, 0, sizeof(bar));
		composite_bind(&bar);
		TEST_ASSERT_STATUS(qaws_curve_eval_adjoint_3d(NULL, c, t, 0xF, &ybar, &bar.views, &tbar));
		rhs = composite_dot(&bar.data, &dir.data) + (double)tbar * tdot;
		if (!diff_close(lhs, rhs, DIFF_TOL * 10))
			ok_adj = 0;
	}
	printf("    composite: jets %s, tangent %s, second tangent %s, adjoint %s\n",
		ok_jet ? "ok" : "NO", ok_fd ? "ok" : "NO", ok_t2 ? "ok" : "NO", ok_adj ? "ok" : "NO");
	TEST_ASSERT(ok_jet, "composite jets match evaluate");
	TEST_ASSERT(ok_fd, "composite tangent matches finite differences through segments");
	TEST_ASSERT(ok_t2, "composite second tangent matches finite differences");
	TEST_ASSERT(ok_adj, "composite adjoint identity across Bezier and NURBS segments");
	qaws_curve_destroy(c);
}

int test_50_diff_curves_main(void)
{
	g_pass = 0;
	g_fail = 0;

	printf("Test 50: Differentiable curve evaluation\n");
	test_curve_families();
	test_masks();
	test_support();
	test_report();
	test_schema();
	test_unit_tangent_composition();
	test_composite();

	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
