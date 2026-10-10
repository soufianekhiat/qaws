/*
 * 2D paths and regions: construction, transforms, measures, point location,
 * and Clipper2's polyline utilities.
 */

#include "qaws_path.h"
#include "qaws_curve.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "qaws_operations.h"
#include "qaws_bezier.h"
#include "qaws_rational_bezier.h"
#include "qaws_bspline.h"
#include "qaws_nurbs.h"
#include "qaws_polynomial.h"
#include "qaws_arc.h"
#include "qaws_clothoid.h"
#include "qaws_composite.h"
#include "qaws_platform.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_kinds.h"
#include "internal/qaws_internal_flatten.h"
#if defined(QAWS_ENABLE_EXACT) && QAWS_ENABLE_EXACT
#include "qaws_exact.h"
#endif
#include <stdlib.h>
#include <string.h>
#include <math.h>

#if QAWS_SCALAR_IS_FLOAT
#define PT_ON_REL   1e-5
#define PT_AREA_REL 1e-6
#else
#define PT_ON_REL   1e-10
#define PT_AREA_REL 1e-14
#endif

#define PT_PI 3.14159265358979323846

int qaws_fill_rule_inside(qaws_fill_rule rule, int w)
{
	switch (rule)
	{
	case QAWS_FILL_EVEN_ODD: return (w & 1) != 0;
	case QAWS_FILL_NON_ZERO: return w != 0;
	case QAWS_FILL_POSITIVE: return w > 0;
	case QAWS_FILL_NEGATIVE: return w < 0;
	default: return 0;
	}
}

static int pt_check_path(qaws_path_2d const* p)
{
	unsigned int i;
	if (!p || (p->curve_count && !p->curves))
		return 0;
	for (i = 0; i < p->curve_count; i++)
		if (!p->curves[i] || p->curves[i]->dimension != QAWS_DIMENSION_2D)
			return 0;
	return 1;
}

static qaws_status pt_eval(qaws_curve const* c, double t, unsigned int flags, double* p, double* d)
{
	qaws_eval_result_2d r;
	qaws_status s = qaws_curve_evaluate_2d(c, (qaws_scalar)t, flags, &r);
	p[0] = r.position.x; p[1] = r.position.y;
	if (d) { d[0] = r.d1.x; d[1] = r.d1.y; }
	return s;
}

/* ======================================================================== */
/*  Construction                                                            */
/* ======================================================================== */

qaws_status qaws_curve_create_polyline_2d(qaws_scalar const* points, unsigned int count, int closed, qaws_curve** out_curve)
{
	unsigned int n, i;
	qaws_scalar *cp, *kn;
	qaws_bspline_desc d;
	qaws_status s;
	if (!points || !out_curve || count < 2)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_curve = NULL;
	n = count + (closed ? 1u : 0u);
	cp = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (2 * n + n + 2));
	if (!cp)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	kn = cp + 2 * n;
	memcpy(cp, points, sizeof(qaws_scalar) * 2 * count);
	if (closed)
	{
		cp[2 * count] = points[0];
		cp[2 * count + 1] = points[1];
	}
	kn[0] = 0;
	for (i = 0; i < n; i++)
		kn[i + 1] = (qaws_scalar)i;
	kn[n + 1] = (qaws_scalar)(n - 1);
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D;
	d.degree = 1;
	d.control_points = cp;
	d.control_point_count = n;
	d.knots = kn;
	d.knot_count = n + 2;
	s = qaws_curve_create_bspline(&d, out_curve);
	free(cp);
	return s;
}

qaws_status qaws_curve_create_ellipse_2d(qaws_scalar cx, qaws_scalar cy, qaws_scalar rx, qaws_scalar ry,
	qaws_scalar rotation, qaws_curve** out_curve)
{
	static double const ux[9] = { 1, 1, 0, -1, -1, -1, 0, 1, 1 }, uy[9] = { 0, 1, 1, 1, 0, -1, -1, -1, 0 };
	qaws_scalar cp[18], w[9], kn[12] = { 0, 0, 0, (qaws_scalar)0.25, (qaws_scalar)0.25, (qaws_scalar)0.5,
		(qaws_scalar)0.5, (qaws_scalar)0.75, (qaws_scalar)0.75, 1, 1, 1 };
	double c = cos((double)rotation), s = sin((double)rotation);
	qaws_nurbs_desc d;
	unsigned int i;
	if (!out_curve || !(rx > 0) || !(ry > 0))
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < 9; i++)
	{
		double x = rx * ux[i], y = ry * uy[i];
		cp[2 * i] = (qaws_scalar)(cx + c * x - s * y);
		cp[2 * i + 1] = (qaws_scalar)(cy + s * x + c * y);
		w[i] = (qaws_scalar)(i % 2 ? sqrt(0.5) : 1.0);
	}
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D;
	d.degree = 2;
	d.control_points = cp;
	d.control_point_count = 9;
	d.knots = kn;
	d.knot_count = 12;
	d.weights = w;
	d.weight_count = 9;
	return qaws_curve_create_nurbs(&d, out_curve);
}

/* ======================================================================== */
/*  Transforms                                                              */
/* ======================================================================== */

static void pt_map(qaws_scalar const m[6], qaws_scalar const* in, qaws_scalar* out, unsigned int n)
{
	unsigned int i;
	for (i = 0; i < n; i++)
	{
		qaws_scalar x = in[2 * i], y = in[2 * i + 1];
		out[2 * i] = m[0] * x + m[1] * y + m[2];
		out[2 * i + 1] = m[3] * x + m[4] * y + m[5];
	}
}

/* A similarity: m's linear part is s R or s R F (uniform scale, rotation,
   optional reflection). Writes the scale, the rotation angle and det sign. */
static int pt_similarity(qaws_scalar const m[6], double* scale, double* angle, int* flip)
{
	double a = m[0], b = m[1], c = m[3], d = m[4];
	double det = a * d - b * c, n2 = a * a + c * c, tol = 1e-12 * (n2 + b * b + d * d);
	if (det > 0)
	{
		/* [a b; c d] = s [cos -sin; sin cos] */
		if (fabs(a - d) > sqrt(tol) || fabs(b + c) > sqrt(tol))
			return 0;
		*flip = 0;
	}
	else if (det < 0)
	{
		/* s [cos sin; sin -cos] */
		if (fabs(a + d) > sqrt(tol) || fabs(b - c) > sqrt(tol))
			return 0;
		*flip = 1;
	}
	else
		return 0;
	*scale = sqrt(n2);
	*angle = atan2(c, a);
	return 1;
}

/* an arc curve as an exact NURBS: pieces of at most a quarter turn */
static qaws_status pt_arc_to_nurbs(qaws_curve const* curve, qaws_curve** out)
{
	qaws_arc_impl const* impl = (qaws_arc_impl const*)curve->impl;
	unsigned int i, k, n = 0, np = 0;
	qaws_scalar *cp, *w, *kn;
	qaws_nurbs_desc d;
	qaws_status s;
	for (i = 0; i < impl->segment_count; i++)
	{
		double sweep = fabs((double)(impl->segments[i].angle_end - impl->segments[i].angle_start));
		n += (unsigned int)ceil(sweep / (PT_PI / 2) - 1e-9) > 0 ? (unsigned int)ceil(sweep / (PT_PI / 2) - 1e-9) : 1;
	}
	cp = (qaws_scalar*)malloc(sizeof(qaws_scalar) * ((2 * n + 1) * 3 + 2 * n + 4));
	if (!cp)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	w = cp + 2 * (2 * n + 1);
	kn = w + (2 * n + 1);
	for (i = 0; i < impl->segment_count; i++)
	{
		qaws_arc_segment const* g = &impl->segments[i];
		double a0 = g->angle_start, sweep = g->angle_end - g->angle_start;
		unsigned int m = (unsigned int)ceil(fabs(sweep) / (PT_PI / 2) - 1e-9), j;
		if (m == 0) m = 1;
		for (j = 0; j < m; j++)
		{
			double t0 = a0 + sweep * j / m, t1 = a0 + sweep * (j + 1) / m, h = (t1 - t0) / 2, tm = t0 + h;
			if (np == 0)
			{
				cp[0] = (qaws_scalar)(g->center[0] + g->radius * cos(t0));
				cp[1] = (qaws_scalar)(g->center[1] + g->radius * sin(t0));
				w[0] = 1;
				np = 1;
			}
			cp[2 * np] = (qaws_scalar)(g->center[0] + g->radius * cos(tm) / cos(h));
			cp[2 * np + 1] = (qaws_scalar)(g->center[1] + g->radius * sin(tm) / cos(h));
			w[np] = (qaws_scalar)cos(h);
			cp[2 * np + 2] = (qaws_scalar)(g->center[0] + g->radius * cos(t1));
			cp[2 * np + 3] = (qaws_scalar)(g->center[1] + g->radius * sin(t1));
			w[np + 1] = 1;
			np += 2;
		}
	}
	/* knots: the arc-length parameter of the source at each piece end */
	k = 0;
	kn[k++] = curve->parameter_range.min_value;
	kn[k++] = curve->parameter_range.min_value;
	kn[k++] = curve->parameter_range.min_value;
	{
		double acc = curve->parameter_range.min_value;
		for (i = 0; i < impl->segment_count; i++)
		{
			qaws_arc_segment const* g = &impl->segments[i];
			double len = g->radius * fabs((double)(g->angle_end - g->angle_start));
			unsigned int m = (unsigned int)ceil(len / g->radius / (PT_PI / 2) - 1e-9), j;
			if (m == 0) m = 1;
			for (j = 0; j < m; j++)
			{
				acc += len / m;
				kn[k++] = (qaws_scalar)acc;
				kn[k++] = (qaws_scalar)acc;
			}
		}
		kn[k - 1] = kn[k - 2] = curve->parameter_range.max_value;
		kn[k++] = curve->parameter_range.max_value;
	}
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D;
	d.degree = 2;
	d.control_points = cp;
	d.control_point_count = np;
	d.knots = kn;
	d.knot_count = k;
	d.weights = w;
	d.weight_count = np;
	s = qaws_curve_create_nurbs(&d, out);
	free(cp);
	return s;
}

qaws_status qaws_curve_transform_2d(qaws_curve const* curve, qaws_scalar const m[6], qaws_curve** out_curve)
{
	qaws_status s;
	double scale, angle;
	int flip;
	if (!curve || !m || !out_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_curve = NULL;
	if (curve->dimension != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;

	switch (curve->kind)
	{
	case QAWS_CURVE_KIND_BEZIER:
	{
		qaws_bezier_impl const* impl = (qaws_bezier_impl const*)curve->impl;
		qaws_scalar* cp = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * impl->control_point_count);
		qaws_bezier_desc d;
		if (!cp) return QAWS_STATUS_ALLOCATION_FAILURE;
		pt_map(m, impl->control_points, cp, impl->control_point_count);
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = curve->degree;
		d.control_points = cp; d.control_point_count = impl->control_point_count;
		s = qaws_curve_create_bezier(&d, out_curve);
		free(cp);
		return s;
	}
	case QAWS_CURVE_KIND_RATIONAL_BEZIER:
	{
		qaws_rational_bezier_impl const* impl = (qaws_rational_bezier_impl const*)curve->impl;
		qaws_scalar* cp = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * impl->control_point_count);
		qaws_rational_bezier_desc d;
		if (!cp) return QAWS_STATUS_ALLOCATION_FAILURE;
		pt_map(m, impl->control_points, cp, impl->control_point_count);
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = curve->degree;
		d.control_points = cp; d.control_point_count = impl->control_point_count;
		d.weights = impl->weights; d.weight_count = impl->control_point_count;
		s = qaws_curve_create_rational_bezier(&d, out_curve);
		free(cp);
		return s;
	}
	case QAWS_CURVE_KIND_BSPLINE:
	{
		qaws_bspline_impl const* impl = (qaws_bspline_impl const*)curve->impl;
		qaws_scalar* cp = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * impl->control_point_count);
		qaws_bspline_desc d;
		if (!cp) return QAWS_STATUS_ALLOCATION_FAILURE;
		pt_map(m, impl->control_points, cp, impl->control_point_count);
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = curve->degree;
		d.control_points = cp; d.control_point_count = impl->control_point_count;
		d.knots = impl->knots; d.knot_count = impl->knot_count;
		s = qaws_curve_create_bspline(&d, out_curve);
		free(cp);
		return s;
	}
	case QAWS_CURVE_KIND_NURBS:
	{
		qaws_nurbs_impl const* impl = (qaws_nurbs_impl const*)curve->impl;
		qaws_scalar* cp = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * impl->control_point_count);
		qaws_nurbs_desc d;
		if (!cp) return QAWS_STATUS_ALLOCATION_FAILURE;
		pt_map(m, impl->control_points, cp, impl->control_point_count);
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = curve->degree;
		d.control_points = cp; d.control_point_count = impl->control_point_count;
		d.knots = impl->knots; d.knot_count = impl->knot_count;
		d.weights = impl->weights; d.weight_count = impl->weight_count;
		s = qaws_curve_create_nurbs(&d, out_curve);
		free(cp);
		return s;
	}
	case QAWS_CURVE_KIND_POLYNOMIAL:
	{
		/* linear part on every coefficient, the translation on c0 only */
		qaws_polynomial_impl const* impl = (qaws_polynomial_impl const*)curve->impl;
		unsigned int n = impl->coefficient_count, i;
		qaws_scalar* co = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * n);
		qaws_polynomial_desc d;
		if (!co) return QAWS_STATUS_ALLOCATION_FAILURE;
		for (i = 0; i < n; i++)
		{
			qaws_scalar x = impl->coefficients[2 * i], y = impl->coefficients[2 * i + 1];
			co[2 * i] = m[0] * x + m[1] * y + (i == 0 ? m[2] : 0);
			co[2 * i + 1] = m[3] * x + m[4] * y + (i == 0 ? m[5] : 0);
		}
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = curve->degree;
		d.coefficients = co; d.coefficient_count = n;
		d.t_min = curve->parameter_range.min_value; d.t_max = curve->parameter_range.max_value;
		s = qaws_curve_create_polynomial(&d, out_curve);
		free(co);
		return s;
	}
	case QAWS_CURVE_KIND_ARC:
	{
		qaws_arc_impl const* impl = (qaws_arc_impl const*)curve->impl;
		qaws_arc_segment* seg;
		qaws_arc_desc d;
		unsigned int i;
		if (!pt_similarity(m, &scale, &angle, &flip))
		{
			qaws_curve* nurbs = NULL;
			s = pt_arc_to_nurbs(curve, &nurbs);
			if (s == QAWS_STATUS_OK)
				s = qaws_curve_transform_2d(nurbs, m, out_curve);
			qaws_curve_destroy(nurbs);
			return s;
		}
		seg = (qaws_arc_segment*)malloc(sizeof(qaws_arc_segment) * impl->segment_count);
		if (!seg) return QAWS_STATUS_ALLOCATION_FAILURE;
		for (i = 0; i < impl->segment_count; i++)
		{
			qaws_arc_segment g = impl->segments[i];
			qaws_scalar c[2];
			pt_map(m, g.center, c, 1);
			g.center[0] = c[0]; g.center[1] = c[1];
			g.radius = (qaws_scalar)(g.radius * scale);
			if (flip)
			{
				/* reflection about the axis at angle/2: theta -> angle - theta */
				g.angle_start = (qaws_scalar)(angle - impl->segments[i].angle_start);
				g.angle_end = (qaws_scalar)(angle - impl->segments[i].angle_end);
			}
			else
			{
				g.angle_start = (qaws_scalar)(impl->segments[i].angle_start + angle);
				g.angle_end = (qaws_scalar)(impl->segments[i].angle_end + angle);
			}
			seg[i] = g;
		}
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.segments = seg; d.segment_count = impl->segment_count;
		s = qaws_curve_create_arc(&d, out_curve);
		free(seg);
		return s;
	}
	case QAWS_CURVE_KIND_CLOTHOID:
	{
		qaws_clothoid_impl const* impl = (qaws_clothoid_impl const*)curve->impl;
		qaws_clothoid_desc d;
		qaws_scalar o[2], src[2];
		if (!pt_similarity(m, &scale, &angle, &flip))
		{
			qaws_curve* fit = NULL;
			s = qaws_internal_curve_fit_cubic(curve, curve->parameter_range.min_value, curve->parameter_range.max_value, &fit);
			if (s == QAWS_STATUS_OK)
				s = qaws_curve_transform_2d(fit, m, out_curve);
			qaws_curve_destroy(fit);
			return s;
		}
		src[0] = impl->origin_x; src[1] = impl->origin_y;
		pt_map(m, src, o, 1);
		memset(&d, 0, sizeof(d));
		d.origin_x = o[0]; d.origin_y = o[1];
		d.start_angle = (qaws_scalar)(flip ? angle - impl->start_angle : impl->start_angle + angle);
		d.start_curvature = (qaws_scalar)((flip ? -impl->kappa_0 : impl->kappa_0) / scale);
		d.end_curvature = (qaws_scalar)((flip ? -impl->kappa_1 : impl->kappa_1) / scale);
		d.length = (qaws_scalar)(impl->length * scale);
		return qaws_curve_create_clothoid(&d, out_curve);
	}
	case QAWS_CURVE_KIND_COMPOSITE:
	{
		qaws_composite_impl const* impl = (qaws_composite_impl const*)curve->impl;
		qaws_curve** seg = (qaws_curve**)malloc(sizeof(qaws_curve*) * impl->segment_count);
		qaws_composite_desc d;
		unsigned int i, made = 0;
		if (!seg) return QAWS_STATUS_ALLOCATION_FAILURE;
		s = QAWS_STATUS_OK;
		for (i = 0; i < impl->segment_count && s == QAWS_STATUS_OK; i++)
			if ((s = qaws_curve_transform_2d(impl->segments[i], m, &seg[i])) == QAWS_STATUS_OK)
				made++;
		if (s == QAWS_STATUS_OK)
		{
			memset(&d, 0, sizeof(d));
			d.dimension = QAWS_DIMENSION_2D; d.segments = seg; d.segment_count = impl->segment_count;
			s = qaws_curve_create_composite(&d, out_curve);
		}
		if (s != QAWS_STATUS_OK)
			for (i = 0; i < made; i++) qaws_curve_destroy(seg[i]);
		free(seg);
		return s;
	}
	default:
	{
		/* Hermite, Catmull-Rom, trajectory, subdivision, Yuksel, wrappers:
		   their exact extracted form, then that */
		qaws_curve* e = NULL;
		s = qaws_curve_extract(curve, curve->parameter_range.min_value, curve->parameter_range.max_value, &e);
		if (s == QAWS_STATUS_OK)
			s = qaws_curve_transform_2d(e, m, out_curve);
		qaws_curve_destroy(e);
		return s;
	}
	}
}

/* ======================================================================== */
/*  Area                                                                    */
/* ======================================================================== */

/* 8- and 16-point Gauss-Legendre nodes / weights on [-1, 1] (half sets) */
static double const pt_g8x[4] = { 0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363 };
static double const pt_g8w[4] = { 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763 };
static double const pt_g16x[8] = { 0.0950125098376374, 0.2816035507792589, 0.4580167776572274, 0.6178762444026438,
	0.7554044083550030, 0.8656312023878318, 0.9445750230732326, 0.9894009349916499 };
static double const pt_g16w[8] = { 0.1894506104550685, 0.1826034150449236, 0.1691565193950025, 0.1495959888165767,
	0.1246289712555339, 0.0951585116824928, 0.0622535239386479, 0.0271524594117541 };

/* integral of x dy over [a, b] with n-point Gauss-Legendre (n = 8 or 16) */
static qaws_status pt_gauss(qaws_curve const* c, double a, double b, int sixteen, double* out)
{
	double mid = 0.5 * (a + b), half = 0.5 * (b - a), sum = 0.0;
	unsigned int i, n = sixteen ? 8 : 4;
	for (i = 0; i < n; i++)
	{
		double xi = sixteen ? pt_g16x[i] : pt_g8x[i], wi = sixteen ? pt_g16w[i] : pt_g8w[i];
		int sg;
		for (sg = -1; sg <= 1; sg += 2)
		{
			double p[2], d[2];
			qaws_status s = pt_eval(c, mid + sg * half * xi, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, p, d);
			if (s != QAWS_STATUS_OK) return s;
			sum += wi * p[0] * d[1];
		}
	}
	*out = sum * half;
	return QAWS_STATUS_OK;
}

static qaws_status pt_area_piece(qaws_curve const* c, double a, double b, double tol, unsigned int depth, double* out)
{
	double g8 = 0.0, g16 = 0.0, l = 0.0, r = 0.0;
	qaws_status s = pt_gauss(c, a, b, 0, &g8);
	if (s == QAWS_STATUS_OK) s = pt_gauss(c, a, b, 1, &g16);
	if (s != QAWS_STATUS_OK) return s;
	if (fabs(g16 - g8) <= tol || depth >= 24)
	{
		*out = g16;
		return QAWS_STATUS_OK;
	}
	s = pt_area_piece(c, a, 0.5 * (a + b), tol * 0.5, depth + 1, &l);
	if (s == QAWS_STATUS_OK) s = pt_area_piece(c, 0.5 * (a + b), b, tol * 0.5, depth + 1, &r);
	*out = l + r;
	return s;
}

/* rough size of a path, from its bounds */
static double pt_extent(qaws_path_2d const* path)
{
	qaws_vec2 lo, hi;
	if (qaws_path_compute_bounds_2d(path, &lo, &hi) != QAWS_STATUS_OK)
		return 1.0;
	return sqrt((double)((hi.x - lo.x) * (hi.x - lo.x) + (hi.y - lo.y) * (hi.y - lo.y)));
}

qaws_status qaws_path_compute_area_2d(qaws_path_2d const* path, qaws_scalar* out_area)
{
	double area = 0.0, ext, tol;
	unsigned int i, s;
	if (!pt_check_path(path) || !out_area)
		return QAWS_STATUS_INVALID_ARGUMENT;
	ext = pt_extent(path);
	tol = ext * ext * PT_AREA_REL;
	for (i = 0; i < path->curve_count; i++)
	{
		qaws_curve const* c = path->curves[i];
		for (s = 0; s < c->span_count; s++)
		{
			double a = c->span_boundaries[s], b = c->span_boundaries[s + 1], v;
			qaws_status st;
			if (!(b > a)) continue;
			st = pt_area_piece(c, a, b, tol / (c->span_count + 1), 0, &v);
			if (st != QAWS_STATUS_OK) return st;
			area += v;
		}
	}
	/* an open path is closed by the straight line back to its start */
	if (!path->closed && path->curve_count)
	{
		qaws_curve const* f = path->curves[0];
		qaws_curve const* l = path->curves[path->curve_count - 1];
		double p0[2], p1[2];
		pt_eval(f, f->parameter_range.min_value, QAWS_EVAL_FLAG_POSITION, p0, NULL);
		pt_eval(l, l->parameter_range.max_value, QAWS_EVAL_FLAG_POSITION, p1, NULL);
		area += 0.5 * (p1[0] + p0[0]) * (p0[1] - p1[1]);
	}
	*out_area = (qaws_scalar)area;
	return QAWS_STATUS_OK;
}

qaws_status qaws_region_compute_area_2d(qaws_path_2d const* paths, unsigned int path_count, qaws_scalar* out_area)
{
	double sum = 0.0;
	unsigned int i;
	if ((!paths && path_count) || !out_area)
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < path_count; i++)
	{
		qaws_scalar a;
		qaws_status s = qaws_path_compute_area_2d(&paths[i], &a);
		if (s != QAWS_STATUS_OK) return s;
		sum += a;
	}
	*out_area = (qaws_scalar)sum;
	return QAWS_STATUS_OK;
}

int qaws_path_is_positive_2d(qaws_path_2d const* path)
{
	qaws_scalar a = 0;
	return qaws_path_compute_area_2d(path, &a) == QAWS_STATUS_OK && a >= 0;
}

qaws_status qaws_path_compute_length_2d(qaws_path_2d const* path, qaws_scalar* out_length)
{
	qaws_scalar sum = 0;
	unsigned int i;
	if (!pt_check_path(path) || !out_length)
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < path->curve_count; i++)
	{
		qaws_curve const* c = path->curves[i];
		qaws_scalar l;
		qaws_status s = qaws_curve_compute_arc_length(c, c->parameter_range.min_value, c->parameter_range.max_value, &l);
		if (s != QAWS_STATUS_OK) return s;
		sum += l;
	}
	*out_length = sum;
	return QAWS_STATUS_OK;
}

/* ======================================================================== */
/*  Bounds                                                                  */
/* ======================================================================== */

/* where coordinate k's derivative vanishes on [a, b]: 32 samples bracket the
   sign changes, bisection takes each to rounding */
static qaws_status pt_extrema(qaws_curve const* c, double a, double b, double* lo, double* hi)
{
	unsigned int const n = 32;
	double prev_t = a, prev_d[2], p[2], d[2];
	unsigned int i, k;
	qaws_status s = pt_eval(c, a, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, p, prev_d);
	if (s != QAWS_STATUS_OK) return s;
	for (i = 1; i <= n; i++)
	{
		double t = a + (b - a) * i / n;
		s = pt_eval(c, t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, p, d);
		if (s != QAWS_STATUS_OK) return s;
		for (k = 0; k < 2; k++)
		{
			if ((prev_d[k] < 0) != (d[k] < 0))
			{
				double l = prev_t, r = t, dl = prev_d[k], q[2], dq[2];
				unsigned int it;
				for (it = 0; it < 60; it++)
				{
					double m = 0.5 * (l + r);
					if (m <= l || m >= r) break;
					pt_eval(c, m, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, q, dq);
					if ((dq[k] < 0) == (dl < 0)) { l = m; dl = dq[k]; } else r = m;
				}
				pt_eval(c, 0.5 * (l + r), QAWS_EVAL_FLAG_POSITION, q, NULL);
				if (q[k] < lo[k]) lo[k] = q[k];
				if (q[k] > hi[k]) hi[k] = q[k];
			}
			prev_d[k] = d[k];
		}
		prev_t = t;
	}
	return QAWS_STATUS_OK;
}

qaws_status qaws_path_compute_bounds_2d(qaws_path_2d const* path, qaws_vec2* out_min, qaws_vec2* out_max)
{
	double lo[2] = { 1e300, 1e300 }, hi[2] = { -1e300, -1e300 };
	unsigned int i, s, k;
	if (!pt_check_path(path) || !out_min || !out_max || path->curve_count == 0)
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < path->curve_count; i++)
	{
		qaws_curve const* c = path->curves[i];
		for (s = 0; s <= c->span_count; s++)
		{
			double p[2];
			qaws_status st = pt_eval(c, c->span_boundaries[s], QAWS_EVAL_FLAG_POSITION, p, NULL);
			if (st != QAWS_STATUS_OK) return st;
			for (k = 0; k < 2; k++)
			{
				if (p[k] < lo[k]) lo[k] = p[k];
				if (p[k] > hi[k]) hi[k] = p[k];
			}
			if (s < c->span_count && c->span_boundaries[s + 1] > c->span_boundaries[s])
			{
				st = pt_extrema(c, c->span_boundaries[s], c->span_boundaries[s + 1], lo, hi);
				if (st != QAWS_STATUS_OK) return st;
			}
		}
	}
	out_min->x = (qaws_scalar)lo[0]; out_min->y = (qaws_scalar)lo[1];
	out_max->x = (qaws_scalar)hi[0]; out_max->y = (qaws_scalar)hi[1];
	return QAWS_STATUS_OK;
}

/* ======================================================================== */
/*  Point location                                                          */
/* ======================================================================== */

typedef struct pt_seg
{
	qaws_curve const* c;
	double t0, t1, p0[2], p1[2], r;
	unsigned int path;
} pt_seg;

typedef struct pt_segs
{
	pt_seg* s;
	unsigned int n, cap;
} pt_segs;

static int pt_push(pt_segs* q, pt_seg const* g)
{
	if (q->n == q->cap)
	{
		unsigned int cap = q->cap ? 2 * q->cap : 256;
		pt_seg* s = (pt_seg*)realloc(q->s, sizeof(pt_seg) * cap);
		if (!s) return 0;
		q->s = s;
		q->cap = cap;
	}
	q->s[q->n++] = *g;
	return 1;
}

static double pt_seg_dist(double const* q, double const* a, double const* b)
{
	double ab[2] = { b[0] - a[0], b[1] - a[1] }, aq[2] = { q[0] - a[0], q[1] - a[1] };
	double l2 = ab[0] * ab[0] + ab[1] * ab[1], u = l2 > 0 ? (aq[0] * ab[0] + aq[1] * ab[1]) / l2 : 0;
	if (u < 0) u = 0;
	if (u > 1) u = 1;
	return hypot(aq[0] - u * ab[0], aq[1] - u * ab[1]);
}

/* deviation of the curve on [t0, t1] from its chord, from three inner points */
static double pt_dev(pt_seg const* g)
{
	double dev = 0.0;
	unsigned int i;
	for (i = 1; i <= 3; i++)
	{
		double p[2], d;
		pt_eval(g->c, g->t0 + (g->t1 - g->t0) * i / 4, QAWS_EVAL_FLAG_POSITION, p, NULL);
		d = pt_seg_dist(p, g->p0, g->p1);
		if (d > dev) dev = d;
	}
	return 2 * dev;
}

/* closest distance from q to the curve of segment g (Gauss-Newton on the
   segment's parameter range, from the chord's nearest point) */
static double pt_true_dist(pt_seg const* g, double const* q)
{
	double ab[2] = { g->p1[0] - g->p0[0], g->p1[1] - g->p0[1] }, l2 = ab[0] * ab[0] + ab[1] * ab[1];
	double u = l2 > 0 ? ((q[0] - g->p0[0]) * ab[0] + (q[1] - g->p0[1]) * ab[1]) / l2 : 0, t, p[2], d[2];
	unsigned int it;
	if (u < 0) u = 0;
	if (u > 1) u = 1;
	t = g->t0 + u * (g->t1 - g->t0);
	for (it = 0; it < 30; it++)
	{
		double dd, step;
		pt_eval(g->c, t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, p, d);
		dd = d[0] * d[0] + d[1] * d[1];
		if (!(dd > 0)) break;
		step = ((p[0] - q[0]) * d[0] + (p[1] - q[1]) * d[1]) / dd;
		t -= step;
		if (t < g->t0) t = g->t0;
		if (t > g->t1) t = g->t1;
		if (fabs(step) <= fabs(g->t1 - g->t0) * 1e-15) break;
	}
	pt_eval(g->c, t, QAWS_EVAL_FLAG_POSITION, p, NULL);
	return hypot(p[0] - q[0], p[1] - q[1]);
}

/* every curve of every path flattened into segments */
static qaws_status pt_flatten(qaws_path_2d const* paths, unsigned int path_count, double flat, pt_segs* out)
{
	unsigned int i, j, k;
	for (i = 0; i < path_count; i++)
		for (j = 0; j < paths[i].curve_count; j++)
		{
			qaws_flat_seg* fs = NULL;
			unsigned int n = 0, cap = 0;
			qaws_status s = qaws_internal_flatten_curve(paths[i].curves[j], 2, (qaws_scalar)flat, i, &fs, &n, &cap);
			if (s != QAWS_STATUS_OK) { free(fs); return s; }
			for (k = 0; k < n; k++)
			{
				pt_seg g;
				g.c = paths[i].curves[j];
				g.path = i;
				g.t0 = fs[k].t0; g.t1 = fs[k].t1;
				g.p0[0] = fs[k].p0[0]; g.p0[1] = fs[k].p0[1];
				g.p1[0] = fs[k].p1[0]; g.p1[1] = fs[k].p1[1];
				g.r = fs[k].r;
				if (!pt_push(out, &g)) { free(fs); return QAWS_STATUS_ALLOCATION_FAILURE; }
			}
			free(fs);
		}
	return QAWS_STATUS_OK;
}

/* winding of the segments [from, to) around q; segments within four times
   their deviation of q are split until the chords cannot cross q, so the
   polygon of chords winds around q as the curves do */
static int pt_wind(pt_segs* q, unsigned int from, unsigned int to, double const* pt, double* angle)
{
	unsigned int i;
	for (i = from; i < to; i++)
	{
		pt_seg g = q->s[i];
		unsigned int depth = 0;
		/* local stack of pieces */
		pt_seg stack[64];
		unsigned int top = 0;
		stack[top++] = g;
		while (top)
		{
			pt_seg h = stack[--top];
			double dist = pt_seg_dist(pt, h.p0, h.p1);
			if (dist <= 4 * h.r && top + 2 <= 64 && fabs(h.t1 - h.t0) > 1e-15 * (fabs(h.t0) + fabs(h.t1) + 1))
			{
				pt_seg a = h, b = h;
				double tm = 0.5 * (h.t0 + h.t1), pm[2];
				pt_eval(h.c, tm, QAWS_EVAL_FLAG_POSITION, pm, NULL);
				a.t1 = tm; a.p1[0] = pm[0]; a.p1[1] = pm[1];
				b.t0 = tm; b.p0[0] = pm[0]; b.p0[1] = pm[1];
				a.r = pt_dev(&a);
				b.r = pt_dev(&b);
				stack[top++] = b;
				stack[top++] = a;
				depth++;
				continue;
			}
			{
				double ax = h.p0[0] - pt[0], ay = h.p0[1] - pt[1], bx = h.p1[0] - pt[0], by = h.p1[1] - pt[1];
				*angle += atan2(ax * by - ay * bx, ax * bx + ay * by);
			}
		}
		(void)depth;
	}
	return 1;
}

qaws_status qaws_region_locate_point_2d(qaws_path_2d const* paths, unsigned int path_count, qaws_fill_rule fill_rule,
	qaws_vec2 point, qaws_scalar tolerance, qaws_point_location* out_location)
{
	pt_segs q;
	double ext = 0.0, tol, pt[2], lo[2] = { 1e300, 1e300 }, hi[2] = { -1e300, -1e300 };
	int winding = 0;
	unsigned int i, start = 0;
	qaws_status s;
	if ((!paths && path_count) || !out_location)
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < path_count; i++)
		if (!pt_check_path(&paths[i]))
			return QAWS_STATUS_INVALID_ARGUMENT;
	*out_location = QAWS_POINT_OUTSIDE;
	memset(&q, 0, sizeof(q));
	/* extent from a first coarse flattening */
	for (i = 0; i < path_count; i++)
	{
		qaws_vec2 a, b;
		if (paths[i].curve_count && qaws_path_compute_bounds_2d(&paths[i], &a, &b) == QAWS_STATUS_OK)
		{
			if (a.x < lo[0]) lo[0] = a.x;
			if (a.y < lo[1]) lo[1] = a.y;
			if (b.x > hi[0]) hi[0] = b.x;
			if (b.y > hi[1]) hi[1] = b.y;
		}
	}
	if (hi[0] >= lo[0])
		ext = hypot(hi[0] - lo[0], hi[1] - lo[1]);
	if (!(ext > 0)) ext = 1.0;
	tol = tolerance > 0 ? (double)tolerance : ext * PT_ON_REL;
	pt[0] = point.x; pt[1] = point.y;
	s = pt_flatten(paths, path_count, ext / 1024, &q);
	if (s != QAWS_STATUS_OK) { free(q.s); return s; }
	/* on a boundary? */
	for (i = 0; i < q.n; i++)
		if (pt_seg_dist(pt, q.s[i].p0, q.s[i].p1) <= q.s[i].r + tol && pt_true_dist(&q.s[i], pt) <= tol)
		{
			*out_location = QAWS_POINT_ON;
			free(q.s);
			return QAWS_STATUS_OK;
		}
	/* winding, path by path: a path's segments are contiguous */
	for (i = 0; i < path_count; i++)
	{
		unsigned int n = 0;
		double angle = 0.0;
		while (start + n < q.n && q.s[start + n].path == i)
			n++;
		pt_wind(&q, start, start + n, pt, &angle);
		/* an open path is closed by its chord back to the start */
		if (!paths[i].closed && n)
		{
			double ax = q.s[start + n - 1].p1[0] - pt[0], ay = q.s[start + n - 1].p1[1] - pt[1];
			double bx = q.s[start].p0[0] - pt[0], by = q.s[start].p0[1] - pt[1];
			angle += atan2(ax * by - ay * bx, ax * bx + ay * by);
		}
		winding += (int)floor(angle / (2 * PT_PI) + 0.5);
		start += n;
	}
	free(q.s);
	*out_location = qaws_fill_rule_inside(fill_rule, winding) ? QAWS_POINT_INSIDE : QAWS_POINT_OUTSIDE;
	return QAWS_STATUS_OK;
}

qaws_status qaws_path_compute_winding_2d(qaws_path_2d const* path, qaws_vec2 point, int* out_winding)
{
	qaws_point_location loc;
	qaws_status s;
	int w;
	if (!out_winding)
		return QAWS_STATUS_INVALID_ARGUMENT;
	/* the non-zero rule tells inside from outside; the sign and size come
	   from the even-odd and positive rules only up to |w| <= 1, so compute
	   the angle sum directly */
	s = qaws_region_locate_point_2d(path, 1, QAWS_FILL_NON_ZERO, point, 0, &loc);
	if (s != QAWS_STATUS_OK)
		return s;
	if (loc == QAWS_POINT_ON)
		return QAWS_STATUS_INVALID_ARGUMENT;
	{
		pt_segs q;
		double angle = 0.0, pt[2] = { point.x, point.y };
		memset(&q, 0, sizeof(q));
		s = pt_flatten(path, 1, pt_extent(path) / 1024, &q);
		if (s == QAWS_STATUS_OK)
		{
			pt_wind(&q, 0, q.n, pt, &angle);
			if (!path->closed && q.n)
			{
				double ax = q.s[q.n - 1].p1[0] - pt[0], ay = q.s[q.n - 1].p1[1] - pt[1];
				double bx = q.s[0].p0[0] - pt[0], by = q.s[0].p0[1] - pt[1];
				angle += atan2(ax * by - ay * bx, ax * bx + ay * by);
			}
		}
		free(q.s);
		w = (int)floor(angle / (2 * PT_PI) + 0.5);
	}
	*out_winding = w;
	return s;
}

/* ======================================================================== */
/*  Polylines                                                               */
/* ======================================================================== */

/* sign of (b - a) x (c - a): exact in f64 builds with the exact kernel */
static int pt_orient(qaws_scalar const* a, qaws_scalar const* b, qaws_scalar const* c)
{
#if defined(QAWS_ENABLE_EXACT) && QAWS_ENABLE_EXACT
	double da[2] = { a[0], a[1] }, db[2] = { b[0], b[1] }, dc[2] = { c[0], c[1] };
	qaws_exact_sign sg = QAWS_EXACT_ZERO;
	/* orient2d(b, c, a) = (b - a) x (c - a) */
	if (qaws_exact_orient2d(db, dc, da, &sg, NULL) == QAWS_STATUS_OK)
		return (int)sg;
#endif
	{
		double v = ((double)b[0] - a[0]) * ((double)c[1] - a[1]) - ((double)b[1] - a[1]) * ((double)c[0] - a[0]);
		return v > 0 ? 1 : (v < 0 ? -1 : 0);
	}
}

static int pt_same(qaws_scalar const* a, qaws_scalar const* b)
{
	return a[0] == b[0] && a[1] == b[1];
}

qaws_status qaws_polyline_strip_duplicates_2d(qaws_scalar const* points, unsigned int count, int closed,
	qaws_scalar* out, unsigned int* out_count)
{
	unsigned int i, n = 0;
	if ((!points && count) || !out || !out_count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < count; i++)
	{
		qaws_scalar x = points[2 * i], y = points[2 * i + 1];
		if (n && out[2 * n - 2] == x && out[2 * n - 1] == y)
			continue;
		out[2 * n] = x;
		out[2 * n + 1] = y;
		n++;
	}
	if (closed)
		while (n > 1 && pt_same(&out[2 * n - 2], &out[0]))
			n--;
	*out_count = n;
	return QAWS_STATUS_OK;
}

qaws_status qaws_polyline_trim_collinear_2d(qaws_scalar const* points, unsigned int count, int closed,
	qaws_scalar* out, unsigned int* out_count)
{
	qaws_scalar* w;
	unsigned int n, i, changed = 1;
	qaws_status s;
	if ((!points && count) || !out || !out_count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	w = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * (count ? count : 1));
	if (!w)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	s = qaws_polyline_strip_duplicates_2d(points, count, closed, w, &n);
	/* drop a point collinear with its neighbours (a spike too: collinear,
	   with the path turning back) until none is left */
	while (changed && n >= 3)
	{
		unsigned int m = 0, first = closed ? 0 : 1, last = closed ? n : n - 1;
		changed = 0;
		if (!closed)
		{
			out[0] = w[0]; out[1] = w[1];
			m = 1;
		}
		for (i = first; i < last; i++)
		{
			qaws_scalar const* p = &w[2 * ((i + n - 1) % n)];
			qaws_scalar const* c = &w[2 * i];
			qaws_scalar const* q = &w[2 * ((i + 1) % n)];
			if (m && !closed)
				p = &out[2 * m - 2];
			if (pt_orient(p, c, q) == 0)
			{
				changed = 1;
				continue;
			}
			out[2 * m] = c[0];
			out[2 * m + 1] = c[1];
			m++;
		}
		if (!closed)
		{
			out[2 * m] = w[2 * n - 2]; out[2 * m + 1] = w[2 * n - 1];
			m++;
		}
		memcpy(w, out, sizeof(qaws_scalar) * 2 * m);
		n = m;
		if (closed && changed)
		{
			/* strip again: removing a spike can leave two equal neighbours */
			qaws_polyline_strip_duplicates_2d(w, n, 1, w, &n);
		}
	}
	if (n < (closed ? 3u : 2u))
		n = 0;
	memcpy(out, w, sizeof(qaws_scalar) * 2 * n);
	free(w);
	*out_count = n;
	return s;
}

static double pt_perp2(qaws_scalar const* p, qaws_scalar const* a, qaws_scalar const* b)
{
	/* squared distance of p from the line through a and b */
	double dx = (double)b[0] - a[0], dy = (double)b[1] - a[1];
	double cr = ((double)p[0] - a[0]) * dy - ((double)p[1] - a[1]) * dx, l2 = dx * dx + dy * dy;
	if (!(l2 > 0))
		return ((double)p[0] - a[0]) * ((double)p[0] - a[0]) + ((double)p[1] - a[1]) * ((double)p[1] - a[1]);
	return cr * cr / l2;
}

qaws_status qaws_polyline_simplify_2d(qaws_scalar const* points, unsigned int count, int closed, qaws_scalar epsilon,
	qaws_scalar* out, unsigned int* out_count)
{
	unsigned int *prev, *next, i, alive, n = count;
	unsigned char* gone;
	double *dist, eps2 = (double)epsilon * epsilon;
	if ((!points && count) || !out || !out_count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (count < 4)
	{
		memmove(out, points, sizeof(qaws_scalar) * 2 * count);
		*out_count = count;
		return QAWS_STATUS_OK;
	}
	prev = (unsigned int*)malloc(sizeof(unsigned int) * 2 * n);
	dist = (double*)malloc(sizeof(double) * n);
	gone = (unsigned char*)calloc(n, 1);
	if (!prev || !dist || !gone)
	{
		free(prev); free(dist); free(gone);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	next = prev + n;
	for (i = 0; i < n; i++)
	{
		prev[i] = (i + n - 1) % n;
		next[i] = (i + 1) % n;
	}
	for (i = 0; i < n; i++)
		dist[i] = (!closed && (i == 0 || i == n - 1)) ? 1e300
			: pt_perp2(&points[2 * i], &points[2 * prev[i]], &points[2 * next[i]]);
	alive = n;
	/* remove the nearest point while it is within epsilon (n^2, as Clipper's
	   result does not depend on the order of ties) */
	while (alive > (closed ? 3u : 2u))
	{
		unsigned int best = n;
		double bd = eps2;
		for (i = 0; i < n; i++)
			if (!gone[i] && dist[i] <= bd)
			{
				if (best == n || dist[i] < bd)
				{
					best = i;
					bd = dist[i];
				}
			}
		if (best == n)
			break;
		gone[best] = 1;
		alive--;
		next[prev[best]] = next[best];
		prev[next[best]] = prev[best];
		{
			unsigned int a = prev[best], b = next[best];
			if (closed || a != 0)
				dist[a] = pt_perp2(&points[2 * a], &points[2 * prev[a]], &points[2 * next[a]]);
			if (closed || b != n - 1)
				dist[b] = pt_perp2(&points[2 * b], &points[2 * prev[b]], &points[2 * next[b]]);
		}
	}
	{
		unsigned int m = 0;
		for (i = 0; i < n; i++)
			if (!gone[i])
			{
				qaws_scalar x = points[2 * i], y = points[2 * i + 1];
				out[2 * m] = x;
				out[2 * m + 1] = y;
				m++;
			}
		*out_count = m;
	}
	free(prev); free(dist); free(gone);
	return QAWS_STATUS_OK;
}

static void pt_rdp(qaws_scalar const* p, unsigned int a, unsigned int b, double eps2, unsigned char* keep)
{
	unsigned int i, far = a;
	double best = -1.0;
	if (b <= a + 1)
		return;
	for (i = a + 1; i < b; i++)
	{
		double d = pt_perp2(&p[2 * i], &p[2 * a], &p[2 * b]);
		if (d > best) { best = d; far = i; }
	}
	if (best > eps2)
	{
		keep[far] = 1;
		pt_rdp(p, a, far, eps2, keep);
		pt_rdp(p, far, b, eps2, keep);
	}
}

qaws_status qaws_polyline_rdp_2d(qaws_scalar const* points, unsigned int count, qaws_scalar epsilon,
	qaws_scalar* out, unsigned int* out_count)
{
	unsigned char* keep;
	unsigned int i, m = 0;
	if ((!points && count) || !out || !out_count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (count < 3)
	{
		memmove(out, points, sizeof(qaws_scalar) * 2 * count);
		*out_count = count;
		return QAWS_STATUS_OK;
	}
	keep = (unsigned char*)calloc(count, 1);
	if (!keep)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	keep[0] = keep[count - 1] = 1;
	pt_rdp(points, 0, count - 1, (double)epsilon * epsilon, keep);
	for (i = 0; i < count; i++)
		if (keep[i])
		{
			qaws_scalar x = points[2 * i], y = points[2 * i + 1];
			out[2 * m] = x;
			out[2 * m + 1] = y;
			m++;
		}
	free(keep);
	*out_count = m;
	return QAWS_STATUS_OK;
}
