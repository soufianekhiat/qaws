/*
 * Test 142: borrowed curves (qaws_view.h)
 *
 *   - parity: prepared and borrowed curves against the owned B-spline and
 *     NURBS curves on the same data (non-uniform knots with a double knot,
 *     clamped uniform default knots, weights and homogeneous points),
 *     position to D3;
 *   - layouts: packed f64, interleaved f32 records with a stride, planar;
 *     f32 data gives what f64 data holding the same values gives;
 *   - lifetime: points rewritten in host memory between evaluations, without
 *     preparing again, evaluate as an owned curve made from the new points;
 *     dependencies name the knots and the layout only;
 *   - range evaluation into f32 and f64 outputs, and from four host threads
 *     over disjoint ranges of one call, bit for bit the serial answers;
 *   - allocation: an allocator failing at each allocation in turn (prepare,
 *     borrowed curve, owned copy), every allocation released; prepare_into
 *     in caller memory (size, alignment, too small);
 *   - validation and extents; operations on a borrowed curve (arc length,
 *     extract, split, reverse) as on the owned curve.
 */

#include "test_common.h"
#include "test_threads.h"
#include "qaws_view.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define CV_N 9
#define CV_DEG 3

#if QAWS_SCALAR_IS_FLOAT
#define CV_TOL 1e-4
#else
#define CV_TOL 1e-11
#endif

static double const cv_knots[CV_N + CV_DEG + 1] = { 0, 0, 0, 0, 0.5, 1.25, 1.25, 2, 3.5, 4, 4, 4, 4 };

static void cv_points(double* p, unsigned int n, double phase)
{
	unsigned int i;
	for (i = 0; i < n; i++)
	{
		p[i * 3 + 0] = cos(0.9 * i + phase) * (1 + 0.1 * i);
		p[i * 3 + 1] = sin(1.3 * i + phase) * 0.8;
		p[i * 3 + 2] = 0.25 * i - 0.05 * i * i;
	}
}

static qaws_array_view cv_packed(void const* data, qaws_element_type t, unsigned int comps, unsigned int n, size_t stride)
{
	qaws_array_view v;
	memset(&v, 0, sizeof v);
	v.data = data;
	v.type = t;
	v.components = comps;
	v.count = n;
	v.stride = stride;
	return v;
}

/* largest difference of the orders both results flag, relative to 1 + |value| */
static double cv_diff3(qaws_eval_result_3d const* a, qaws_eval_result_3d const* b)
{
	qaws_vec3 const* va[4] = { &a->position, &a->d1, &a->d2, &a->d3 };
	qaws_vec3 const* vb[4] = { &b->position, &b->d1, &b->d2, &b->d3 };
	double worst = 0;
	unsigned int r;
	for (r = 0; r < 4; r++)
		if (a->valid_flags & b->valid_flags & (1u << r))
		{
			double e = fabs((double)va[r]->x - vb[r]->x) + fabs((double)va[r]->y - vb[r]->y) + fabs((double)va[r]->z - vb[r]->z);
			double s = 1 + fabs((double)vb[r]->x) + fabs((double)vb[r]->y) + fabs((double)vb[r]->z);
			worst = e / s > worst ? e / s : worst;
		}
	return worst;
}

static double cv_diff2(qaws_eval_result_2d const* a, qaws_eval_result_2d const* b)
{
	qaws_vec2 const* va[4] = { &a->position, &a->d1, &a->d2, &a->d3 };
	qaws_vec2 const* vb[4] = { &b->position, &b->d1, &b->d2, &b->d3 };
	double worst = 0;
	unsigned int r;
	for (r = 0; r < 4; r++)
		if (a->valid_flags & b->valid_flags & (1u << r))
		{
			double e = fabs((double)va[r]->x - vb[r]->x) + fabs((double)va[r]->y - vb[r]->y);
			double s = 1 + fabs((double)vb[r]->x) + fabs((double)vb[r]->y);
			worst = e / s > worst ? e / s : worst;
		}
	return worst;
}

#define CV_ALL (QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3)

static qaws_curve* cv_owned_bspline(double const* p, unsigned int n, double const* knots)
{
	qaws_scalar cp[CV_N * 3], k[CV_N + CV_DEG + 1];
	qaws_bspline_desc d;
	qaws_curve* c = NULL;
	unsigned int i;
	for (i = 0; i < n * 3; i++)
		cp[i] = (qaws_scalar)p[i];
	for (i = 0; knots && i < n + CV_DEG + 1; i++)
		k[i] = (qaws_scalar)knots[i];
	memset(&d, 0, sizeof d);
	d.dimension = QAWS_DIMENSION_3D;
	d.degree = CV_DEG;
	d.control_points = cp;
	d.control_point_count = n;
	d.knots = knots ? k : NULL;
	d.knot_count = knots ? n + CV_DEG + 1 : 0;
	qaws_curve_create_bspline(&d, &c);
	return c;
}

/* worst difference over 41 parameters of the domain, prepared against owned (3D) */
static double cv_compare3(qaws_curve_prepared const* p, qaws_curve const* owned)
{
	qaws_range r = qaws_curve_get_parameter_range(owned);
	double worst = 0, e;
	unsigned int i;
	for (i = 0; i <= 40; i++)
	{
		qaws_scalar t = r.min_value + (r.max_value - r.min_value) * (qaws_scalar)i / 40;
		qaws_eval_result_3d a, b;
		if (qaws_curve_prepared_evaluate_3d(p, t, CV_ALL, &a) != QAWS_STATUS_OK || qaws_curve_evaluate_3d(owned, t, CV_ALL, &b) != QAWS_STATUS_OK
			|| (b.valid_flags & ~a.valid_flags))
			return 1e30;
		e = cv_diff3(&a, &b);
		worst = e > worst ? e : worst;
	}
	return worst;
}

static void test_parity_bspline(void)
{
	double p[CV_N * 3];
	qaws_curve_view v;
	qaws_curve_prepared* pr = NULL;
	qaws_curve *owned, *owned_u, *borrowed = NULL;
	qaws_range ro, rb;
	double e, eu, eb = 0;
	unsigned int i;
	char msg[256];
	cv_points(p, CV_N, 0.3);
	memset(&v, 0, sizeof v);
	v.dimension = QAWS_DIMENSION_3D;
	v.degree = CV_DEG;
	v.points = cv_packed(p, QAWS_ELEMENT_F64, 3, CV_N, 0);
	v.knots = cv_packed(cv_knots, QAWS_ELEMENT_F64, 1, CV_N + CV_DEG + 1, 0);
	owned = cv_owned_bspline(p, CV_N, cv_knots);
	owned_u = cv_owned_bspline(p, CV_N, NULL);
	TEST_ASSERT(owned && owned_u, "owned B-splines created");
	if (!owned || !owned_u)
		return;
	TEST_ASSERT(qaws_curve_prepare(&v, NULL, &pr) == QAWS_STATUS_OK, "prepare: packed f64, non-uniform knots");
	e = cv_compare3(pr, owned);
	TEST_ASSERT(qaws_curve_prepared_dependencies(pr) == (QAWS_DEPENDS_KNOTS | QAWS_DEPENDS_LAYOUT), "dependencies: knots and layout");
	qaws_curve_prepared_destroy(pr);
	v.knots.count = 0;
	TEST_ASSERT(qaws_curve_prepare(&v, NULL, &pr) == QAWS_STATUS_OK, "prepare: default knots");
	eu = cv_compare3(pr, owned_u);
	TEST_ASSERT(qaws_curve_prepared_dependencies(pr) == QAWS_DEPENDS_LAYOUT, "dependencies without knots: layout");
	qaws_curve_prepared_destroy(pr);
	v.knots.count = CV_N + CV_DEG + 1;
	TEST_ASSERT(qaws_curve_create_borrowed(&v, NULL, &borrowed) == QAWS_STATUS_OK, "borrowed curve created");
	if (borrowed)
	{
		ro = qaws_curve_get_parameter_range(owned);
		rb = qaws_curve_get_parameter_range(borrowed);
		TEST_ASSERT(qaws_curve_get_kind(borrowed) == QAWS_CURVE_KIND_VIEW && qaws_curve_get_span_count(borrowed) == qaws_curve_get_span_count(owned)
			&& ro.min_value == rb.min_value && ro.max_value == rb.max_value && !qaws_curve_is_rational(borrowed)
			&& qaws_curve_get_continuity(borrowed) == QAWS_CONTINUITY_C1,
			"borrowed curve: kind, spans, range, polynomial, C1 at the double knot");
		for (i = 0; i <= 40; i++)
		{
			qaws_scalar t = ro.min_value + (ro.max_value - ro.min_value) * (qaws_scalar)i / 40;
			qaws_eval_result_3d a, b;
			qaws_curve_evaluate_3d(borrowed, t, CV_ALL, &a);
			qaws_curve_evaluate_3d(owned, t, CV_ALL, &b);
			e = cv_diff3(&a, &b) > e ? cv_diff3(&a, &b) : e;
			eb = cv_diff3(&a, &b) > eb ? cv_diff3(&a, &b) : eb;
		}
		qaws_curve_destroy(borrowed);
	}
	sprintf(msg, "B-spline parity, position to D3: prepared %.2g, default knots %.2g, borrowed curve %.2g", e, eu, eb);
	TEST_ASSERT(e <= CV_TOL && eu <= CV_TOL && eb <= CV_TOL, msg);
	printf("  %s\n", msg);
	qaws_curve_destroy(owned);
	qaws_curve_destroy(owned_u);
}

static void test_parity_nurbs(void)
{
	/* a 2D quadratic NURBS: three quarter circles and a bump, weights 1 and sqrt(1/2) */
	double const s = sqrt(0.5);
	double xy[7 * 2] = { 1, 0, 1, 1, 0, 1, -1, 1, -1, 0, -1, -1, 0, -1.5 };
	double w[7] = { 1, s, 1, s, 1, 0.7, 1 };
	double hom[7 * 3];
	double kn[10] = { 0, 0, 0, 1, 1, 2, 2, 3, 3, 3 };
	qaws_scalar cps[14], ws[7], ks[10];
	qaws_nurbs_desc nd;
	qaws_curve* owned = NULL;
	qaws_curve_view v;
	qaws_curve_prepared *pw = NULL, *ph = NULL;
	double ew = 0, eh = 0, e;
	unsigned int i;
	char msg[256];
	for (i = 0; i < 14; i++)
		cps[i] = (qaws_scalar)xy[i];
	for (i = 0; i < 7; i++)
	{
		ws[i] = (qaws_scalar)w[i];
		hom[i * 3 + 0] = xy[i * 2] * w[i];
		hom[i * 3 + 1] = xy[i * 2 + 1] * w[i];
		hom[i * 3 + 2] = w[i];
	}
	for (i = 0; i < 10; i++)
		ks[i] = (qaws_scalar)kn[i];
	memset(&nd, 0, sizeof nd);
	nd.dimension = QAWS_DIMENSION_2D;
	nd.degree = 2;
	nd.control_points = cps;
	nd.control_point_count = 7;
	nd.knots = ks;
	nd.knot_count = 10;
	nd.weights = ws;
	nd.weight_count = 7;
	TEST_ASSERT(qaws_curve_create_nurbs(&nd, &owned) == QAWS_STATUS_OK, "owned NURBS created");
	memset(&v, 0, sizeof v);
	v.dimension = QAWS_DIMENSION_2D;
	v.degree = 2;
	v.points = cv_packed(xy, QAWS_ELEMENT_F64, 2, 7, 0);
	v.weights = cv_packed(w, QAWS_ELEMENT_F64, 1, 7, 0);
	v.knots = cv_packed(kn, QAWS_ELEMENT_F64, 1, 10, 0);
	TEST_ASSERT(qaws_curve_prepare(&v, NULL, &pw) == QAWS_STATUS_OK, "prepare: weights");
	v.points = cv_packed(hom, QAWS_ELEMENT_F64, 3, 7, 0);
	memset(&v.weights, 0, sizeof v.weights);
	TEST_ASSERT(qaws_curve_prepare(&v, NULL, &ph) == QAWS_STATUS_OK, "prepare: homogeneous points");
	if (!owned || !pw || !ph)
		return;
	for (i = 0; i <= 60; i++)
	{
		qaws_scalar t = (qaws_scalar)(3.0 * i / 60);
		qaws_eval_result_2d a, b, c;
		qaws_curve_prepared_evaluate_2d(pw, t, CV_ALL, &a);
		qaws_curve_prepared_evaluate_2d(ph, t, CV_ALL, &b);
		qaws_curve_evaluate_2d(owned, t, CV_ALL, &c);
		e = cv_diff2(&a, &c);
		ew = e > ew ? e : ew;
		e = cv_diff2(&b, &c);
		eh = e > eh ? e : eh;
	}
	sprintf(msg, "NURBS parity, position to D3: weights %.2g, homogeneous points %.2g", ew, eh);
	TEST_ASSERT(ew <= CV_TOL && eh <= CV_TOL * 10, msg);
	printf("  %s\n", msg);
	qaws_curve_prepared_destroy(pw);
	qaws_curve_prepared_destroy(ph);
	qaws_curve_destroy(owned);
}

/* host records: a float position among other fields */
typedef struct cv_record
{
	int id;
	float position[3];
	float color[4];
} cv_record;

static void test_layouts(void)
{
	double p[CV_N * 3], pf[CV_N * 3], px[CV_N], py[CV_N], pz[CV_N];
	cv_record rec[CV_N];
	float kf[CV_N + CV_DEG + 1];
	qaws_curve_view v;
	qaws_curve_prepared *a = NULL, *b = NULL, *c = NULL;
	unsigned int i, diff = 0;
	cv_points(p, CV_N, 1.1);
	for (i = 0; i < CV_N; i++)
	{
		memset(&rec[i], 0, sizeof rec[i]);
		rec[i].id = (int)i;
		rec[i].position[0] = (float)p[i * 3];
		rec[i].position[1] = (float)p[i * 3 + 1];
		rec[i].position[2] = (float)p[i * 3 + 2];
		pf[i * 3] = px[i] = rec[i].position[0];
		pf[i * 3 + 1] = py[i] = rec[i].position[1];
		pf[i * 3 + 2] = pz[i] = rec[i].position[2];
	}
	for (i = 0; i < CV_N + CV_DEG + 1; i++)
		kf[i] = (float)cv_knots[i];
	memset(&v, 0, sizeof v);
	v.dimension = QAWS_DIMENSION_3D;
	v.degree = CV_DEG;
	/* f32 records, f32 knots */
	v.points = cv_packed(rec[0].position, QAWS_ELEMENT_F32, 3, CV_N, sizeof(cv_record));
	v.knots = cv_packed(kf, QAWS_ELEMENT_F32, 1, CV_N + CV_DEG + 1, 0);
	TEST_ASSERT(qaws_curve_prepare(&v, NULL, &a) == QAWS_STATUS_OK, "prepare: interleaved f32 records, f32 knots");
	/* the same values as packed f64 */
	v.points = cv_packed(pf, QAWS_ELEMENT_F64, 3, CV_N, 0);
	v.knots = cv_packed(cv_knots, QAWS_ELEMENT_F64, 1, CV_N + CV_DEG + 1, 0);
	TEST_ASSERT(qaws_curve_prepare(&v, NULL, &b) == QAWS_STATUS_OK, "prepare: packed f64");
	/* planar f64 */
	memset(&v.points, 0, sizeof v.points);
	v.points.planes[0] = px;
	v.points.planes[1] = py;
	v.points.planes[2] = pz;
	v.points.type = QAWS_ELEMENT_F64;
	v.points.components = 3;
	v.points.count = CV_N;
	TEST_ASSERT(qaws_curve_prepare(&v, NULL, &c) == QAWS_STATUS_OK, "prepare: planar f64");
	if (!a || !b || !c)
		return;
	for (i = 0; i <= 50; i++)
	{
		qaws_scalar t = (qaws_scalar)(4.0 * i / 50);
		qaws_eval_result_3d ra, rb, rc;
		qaws_curve_prepared_evaluate_3d(a, t, CV_ALL, &ra);
		qaws_curve_prepared_evaluate_3d(b, t, CV_ALL, &rb);
		qaws_curve_prepared_evaluate_3d(c, t, CV_ALL, &rc);
		diff += memcmp(&ra, &rb, sizeof ra) != 0 || memcmp(&rb, &rc, sizeof rb) != 0;
	}
	TEST_ASSERT(diff == 0, "f32 records, packed f64 and planar f64 holding the same values evaluate bit for bit alike");
	qaws_curve_prepared_destroy(a);
	qaws_curve_prepared_destroy(b);
	qaws_curve_prepared_destroy(c);
}

static void test_mutation(void)
{
	double p[CV_N * 3];
	qaws_curve_view v;
	qaws_curve_prepared* pr = NULL;
	qaws_curve* owned;
	double e0, e1;
	char msg[200];
	cv_points(p, CV_N, 0.0);
	memset(&v, 0, sizeof v);
	v.dimension = QAWS_DIMENSION_3D;
	v.degree = CV_DEG;
	v.points = cv_packed(p, QAWS_ELEMENT_F64, 3, CV_N, 0);
	v.knots = cv_packed(cv_knots, QAWS_ELEMENT_F64, 1, CV_N + CV_DEG + 1, 0);
	if (qaws_curve_prepare(&v, NULL, &pr) != QAWS_STATUS_OK)
	{
		TEST_ASSERT(0, "mutation: prepared");
		return;
	}
	owned = cv_owned_bspline(p, CV_N, cv_knots);
	e0 = cv_compare3(pr, owned);
	qaws_curve_destroy(owned);
	/* the host rewrites its points in place (an animation frame); the prepared object is kept */
	cv_points(p, CV_N, 2.0);
	owned = cv_owned_bspline(p, CV_N, cv_knots);
	e1 = cv_compare3(pr, owned);
	qaws_curve_destroy(owned);
	sprintf(msg, "points rewritten in place without preparing again: before %.2g, after %.2g", e0, e1);
	TEST_ASSERT(e0 <= CV_TOL && e1 <= CV_TOL, msg);
	qaws_curve_prepared_destroy(pr);
}

#define CV_RANGE 4096
#define CV_THREADS 4

typedef struct cv_job
{
	qaws_curve_prepared const* p;
	qaws_array_view const* params;
	qaws_curve_outputs const* out;
	unsigned int begin, end;
	qaws_status status;
} cv_job;

static void cv_range_job(void* a)
{
	cv_job* j = (cv_job*)a;
	unsigned int r;
	for (r = 0; r < 4; r++)
		j->status |= qaws_curve_prepared_evaluate_range(j->p, j->params, j->begin, j->end, j->out);
}

static void test_ranges(void)
{
	static double params[CV_RANGE], serial[CV_RANGE * 3], d1s[CV_RANGE * 3], pos[CV_RANGE * 3], d1[CV_RANGE * 3];
	static float pos32[CV_RANGE * 3];
	double p[CV_N * 3];
	qaws_curve_view v;
	qaws_curve_prepared* pr = NULL;
	qaws_array_view pv;
	qaws_curve_outputs o;
	cv_job job[CV_THREADS];
	void* args[CV_THREADS];
	unsigned int i, bad = 0, bad32 = 0, st = 0;
	char msg[256];
	cv_points(p, CV_N, 0.7);
	for (i = 0; i < CV_RANGE; i++)
		params[i] = 4.0 * i / (CV_RANGE - 1);
	memset(&v, 0, sizeof v);
	v.dimension = QAWS_DIMENSION_3D;
	v.degree = CV_DEG;
	v.points = cv_packed(p, QAWS_ELEMENT_F64, 3, CV_N, 0);
	v.knots = cv_packed(cv_knots, QAWS_ELEMENT_F64, 1, CV_N + CV_DEG + 1, 0);
	if (qaws_curve_prepare(&v, NULL, &pr) != QAWS_STATUS_OK)
	{
		TEST_ASSERT(0, "ranges: prepared");
		return;
	}
	pv = cv_packed(params, QAWS_ELEMENT_F64, 1, CV_RANGE, 0);
	/* serial: position and D1 in one call, against single evaluations */
	memset(&o, 0, sizeof o);
	o.position.data = serial;
	o.position.type = QAWS_ELEMENT_F64;
	o.position.components = 3;
	o.position.count = CV_RANGE;
	o.d1 = o.position;
	o.d1.data = d1s;
	TEST_ASSERT(qaws_curve_prepared_evaluate_range(pr, &pv, 0, CV_RANGE, &o) == QAWS_STATUS_OK, "range evaluation");
	for (i = 0; i < CV_RANGE; i += 7)
	{
		qaws_eval_result_3d r;
		qaws_curve_prepared_evaluate_3d(pr, (qaws_scalar)params[i], QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &r);
		bad += (double)r.position.x != serial[i * 3] || (double)r.position.y != serial[i * 3 + 1] || (double)r.position.z != serial[i * 3 + 2]
			|| (double)r.d1.x != d1s[i * 3] || (double)r.d1.z != d1s[i * 3 + 2];
	}
	TEST_ASSERT(bad == 0, "range evaluation: the single evaluations' values");
	/* f32 output */
	memset(&o, 0, sizeof o);
	o.position.data = pos32;
	o.position.type = QAWS_ELEMENT_F32;
	o.position.components = 3;
	o.position.count = CV_RANGE;
	TEST_ASSERT(qaws_curve_prepared_evaluate_range(pr, &pv, 0, CV_RANGE, &o) == QAWS_STATUS_OK, "range evaluation into f32");
	for (i = 0; i < CV_RANGE * 3; i++)
		bad32 += pos32[i] != (float)serial[i];
	TEST_ASSERT(bad32 == 0, "f32 outputs: the f64 values rounded");
	/* four threads, disjoint ranges of one call, four times each */
	memset(&o, 0, sizeof o);
	o.position.data = pos;
	o.position.type = QAWS_ELEMENT_F64;
	o.position.components = 3;
	o.position.count = CV_RANGE;
	o.d1 = o.position;
	o.d1.data = d1;
	for (i = 0; i < CV_THREADS; i++)
	{
		job[i].p = pr;
		job[i].params = &pv;
		job[i].out = &o;
		job[i].begin = CV_RANGE * i / CV_THREADS;
		job[i].end = CV_RANGE * (i + 1) / CV_THREADS;
		job[i].status = QAWS_STATUS_OK;
		args[i] = &job[i];
	}
	test_run_threads(CV_THREADS, cv_range_job, args);
	for (i = 0; i < CV_THREADS; i++)
		st |= (unsigned int)job[i].status;
	sprintf(msg, "%u threads over disjoint ranges of %u parameters: bit for bit the serial outputs", CV_THREADS, CV_RANGE);
	TEST_ASSERT(st == 0 && memcmp(pos, serial, sizeof pos) == 0 && memcmp(d1, d1s, sizeof d1) == 0, msg);
	/* refusals */
	o.d1.count = 10;
	TEST_ASSERT(qaws_curve_prepared_evaluate_range(pr, &pv, 0, CV_RANGE, &o) == QAWS_STATUS_BUFFER_TOO_SMALL, "range: short output refused");
	o.d1.count = CV_RANGE;
	o.d1.components = 2;
	TEST_ASSERT(qaws_curve_prepared_evaluate_range(pr, &pv, 0, CV_RANGE, &o) == QAWS_STATUS_INVALID_DIMENSION, "range: output of 2 components refused");
	qaws_curve_prepared_destroy(pr);
}

/* an allocator failing at allocation number `fail_at` (counting from 0), counting what is live */
typedef struct cv_alloc_state
{
	int calls, fail_at, live;
} cv_alloc_state;

static void* cv_alloc(unsigned long size, void* u)
{
	cv_alloc_state* s = (cv_alloc_state*)u;
	if (s->calls++ == s->fail_at)
		return NULL;
	s->live++;
	return malloc(size);
}

static void cv_dealloc(void* ptr, void* u)
{
	if (ptr)
		((cv_alloc_state*)u)->live--;
	free(ptr);
}

static void test_allocation(void)
{
	double p[CV_N * 3];
	qaws_curve_view v;
	cv_alloc_state s;
	qaws_allocator al;
	unsigned int f, which, leaks = 0, wrong = 0, ok[3] = { 0, 0, 0 };
	unsigned char* mem;
	size_t bytes;
	qaws_curve_prepared *pi = NULL, *pa = NULL;
	char msg[256];
	cv_points(p, CV_N, 0.2);
	memset(&v, 0, sizeof v);
	v.dimension = QAWS_DIMENSION_3D;
	v.degree = CV_DEG;
	v.points = cv_packed(p, QAWS_ELEMENT_F64, 3, CV_N, 0);
	v.knots = cv_packed(cv_knots, QAWS_ELEMENT_F64, 1, CV_N + CV_DEG + 1, 0);
	al.alloc = cv_alloc;
	al.dealloc = cv_dealloc;
	al.user_data = &s;
	for (which = 0; which < 3; which++)
		for (f = 0; f < 16 && !ok[which]; f++)
		{
			qaws_curve_prepared* pr = NULL;
			qaws_curve* c = NULL;
			qaws_status st;
			s.calls = 0;
			s.fail_at = (int)f;
			s.live = 0;
			st = which == 0 ? qaws_curve_prepare(&v, &al, &pr) : which == 1 ? qaws_curve_create_borrowed(&v, &al, &c)
				: qaws_curve_create_from_view(&v, &al, &c);
			if (st == QAWS_STATUS_OK)
			{
				ok[which] = f + 1;
				qaws_curve_prepared_destroy(pr);
				qaws_curve_destroy(c);
			}
			else
				wrong += st != QAWS_STATUS_ALLOCATION_FAILURE || pr || c;
			leaks += s.live != 0;
		}
	sprintf(msg, "allocation failures: prepare, borrowed, copy succeed once %u, %u, %u allocations pass; %u wrong statuses, %u leaks",
		ok[0] - 1, ok[1] - 1, ok[2] - 1, wrong, leaks);
	TEST_ASSERT(ok[0] == 2 && ok[1] && ok[2] && wrong == 0 && leaks == 0, msg);
	/* caller memory */
	bytes = qaws_curve_prepared_size(&v);
	mem = (unsigned char*)malloc(bytes + 32);
	{
		unsigned char* at = (unsigned char*)(((size_t)mem + 15) & ~(size_t)15);
		TEST_ASSERT(bytes > 0 && qaws_curve_prepare_into(&v, at, bytes - 1, &pi) == QAWS_STATUS_BUFFER_TOO_SMALL && !pi, "prepare_into: one byte short refused");
		TEST_ASSERT(qaws_curve_prepare_into(&v, at + 8, bytes, &pi) == QAWS_STATUS_INVALID_ARGUMENT, "prepare_into: misaligned memory refused");
		TEST_ASSERT(qaws_curve_prepare_into(&v, at, bytes, &pi) == QAWS_STATUS_OK, "prepare_into: in caller memory");
		TEST_ASSERT(qaws_curve_prepare(&v, NULL, &pa) == QAWS_STATUS_OK, "prepare: allocated");
		if (pi && pa)
		{
			qaws_eval_result_3d a, b;
			qaws_curve_prepared_evaluate_3d(pi, (qaws_scalar)1.7, CV_ALL, &a);
			qaws_curve_prepared_evaluate_3d(pa, (qaws_scalar)1.7, CV_ALL, &b);
			TEST_ASSERT(memcmp(&a, &b, sizeof a) == 0, "prepare_into evaluates as prepare");
		}
		qaws_curve_prepared_destroy(pi);
		qaws_curve_prepared_destroy(pa);
	}
	free(mem);
}

static void test_validation(void)
{
	double p[CV_N * 3], kn[CV_N + CV_DEG + 1], w[CV_N];
	qaws_curve_view v;
	qaws_curve_prepared* pr = (qaws_curve_prepared*)&v;
	unsigned int i;
	cv_points(p, CV_N, 0.0);
	memcpy(kn, cv_knots, sizeof kn);
	for (i = 0; i < CV_N; i++)
		w[i] = 1;
	memset(&v, 0, sizeof v);
	v.dimension = QAWS_DIMENSION_3D;
	v.degree = CV_DEG;
	v.points = cv_packed(p, QAWS_ELEMENT_F64, 3, CV_N, 0);
	v.knots = cv_packed(kn, QAWS_ELEMENT_F64, 1, CV_N + CV_DEG + 1, 0);
	TEST_ASSERT(qaws_curve_view_validate(&v) == QAWS_STATUS_OK, "valid view");
	kn[5] = 0.1;
	TEST_ASSERT(qaws_curve_view_validate(&v) == QAWS_STATUS_INVALID_KNOT_VECTOR, "decreasing knots refused");
	TEST_ASSERT(qaws_curve_prepare(&v, NULL, &pr) == QAWS_STATUS_INVALID_KNOT_VECTOR && !pr, "prepare refuses, output cleared");
	kn[5] = cv_knots[5];
	v.knots.count = 12;
	TEST_ASSERT(qaws_curve_view_validate(&v) == QAWS_STATUS_INVALID_KNOT_COUNT, "knot count refused");
	v.knots.count = CV_N + CV_DEG + 1;
	v.degree = 32;
	TEST_ASSERT(qaws_curve_view_validate(&v) == QAWS_STATUS_INVALID_DEGREE, "degree 32 refused");
	v.degree = CV_DEG;
	v.weights = cv_packed(w, QAWS_ELEMENT_F64, 1, CV_N, 0);
	w[3] = 0;
	TEST_ASSERT(qaws_curve_view_validate(&v) == QAWS_STATUS_INVALID_ARGUMENT, "zero weight refused");
	w[3] = 1;
	v.weights.count = CV_N - 1;
	TEST_ASSERT(qaws_curve_view_validate(&v) == QAWS_STATUS_INVALID_WEIGHT_COUNT, "weight count refused");
	memset(&v.weights, 0, sizeof v.weights);
	p[4] = NAN;
	TEST_ASSERT(qaws_curve_view_validate(&v) == QAWS_STATUS_INVALID_ARGUMENT, "NaN point refused");
	p[4] = 0;
	v.points.stride = 12;
	TEST_ASSERT(qaws_curve_view_validate(&v) == QAWS_STATUS_INVALID_ARGUMENT, "stride shorter than a point refused");
	v.points.stride = 0;
	TEST_ASSERT(qaws_array_view_check_extent(&v.points, p, sizeof p) == QAWS_STATUS_OK, "extent: inside the host buffer");
	TEST_ASSERT(qaws_array_view_check_extent(&v.points, p, sizeof p - 1) == QAWS_STATUS_OUT_OF_RANGE, "extent: one byte past refused");
	v.points.stride = 32;
	TEST_ASSERT(qaws_array_view_check_extent(&v.points, p, sizeof p) == QAWS_STATUS_OUT_OF_RANGE, "extent: a stride reading past refused");
}

static void test_operations(void)
{
	double p[CV_N * 3];
	qaws_curve_view v;
	qaws_curve *owned, *borrowed = NULL, *eo = NULL, *eb = NULL, *l = NULL, *r = NULL, *rv = NULL;
	qaws_scalar lo = 0, lb = 0;
	double e = 0;
	unsigned int i;
	char msg[256];
	cv_points(p, CV_N, 0.9);
	memset(&v, 0, sizeof v);
	v.dimension = QAWS_DIMENSION_3D;
	v.degree = CV_DEG;
	v.points = cv_packed(p, QAWS_ELEMENT_F64, 3, CV_N, 0);
	v.knots = cv_packed(cv_knots, QAWS_ELEMENT_F64, 1, CV_N + CV_DEG + 1, 0);
	owned = cv_owned_bspline(p, CV_N, cv_knots);
	if (!owned || qaws_curve_create_borrowed(&v, NULL, &borrowed) != QAWS_STATUS_OK)
	{
		TEST_ASSERT(0, "operations: curves created");
		return;
	}
	qaws_curve_compute_arc_length(owned, 0, 4, &lo);
	qaws_curve_compute_arc_length(borrowed, 0, 4, &lb);
	sprintf(msg, "arc length: owned %.12g, borrowed %.12g", (double)lo, (double)lb);
	TEST_ASSERT(fabs((double)lo - lb) <= CV_TOL * (1 + fabs((double)lo)), msg);
	TEST_ASSERT(qaws_curve_extract(owned, (qaws_scalar)0.3, (qaws_scalar)3.1, &eo) == QAWS_STATUS_OK
		&& qaws_curve_extract(borrowed, (qaws_scalar)0.3, (qaws_scalar)3.1, &eb) == QAWS_STATUS_OK, "extract: owned and borrowed");
	TEST_ASSERT(qaws_curve_split(borrowed, (qaws_scalar)1.25, &l, &r) == QAWS_STATUS_OK, "split of a borrowed curve");
	TEST_ASSERT(qaws_curve_reverse(borrowed, &rv) == QAWS_STATUS_OK, "reverse of a borrowed curve");
	if (eo && eb && l && r && rv)
	{
		for (i = 0; i <= 20; i++)
		{
			qaws_scalar t = (qaws_scalar)(0.3 + 2.8 * i / 20);
			qaws_eval_result_3d a, b, c;
			qaws_curve_evaluate_3d(eo, t, QAWS_EVAL_FLAG_POSITION, &a);
			qaws_curve_evaluate_3d(eb, t, QAWS_EVAL_FLAG_POSITION, &b);
			e = cv_diff3(&a, &b) > e ? cv_diff3(&a, &b) : e;
			qaws_curve_evaluate_3d(owned, t, QAWS_EVAL_FLAG_POSITION, &c);
			qaws_curve_evaluate_3d(t < (qaws_scalar)1.25 ? l : r, t, QAWS_EVAL_FLAG_POSITION, &a);
			e = cv_diff3(&a, &c) > e ? cv_diff3(&a, &c) : e;
			qaws_curve_evaluate_3d(rv, 4 - t, QAWS_EVAL_FLAG_POSITION, &a);
			e = cv_diff3(&a, &c) > e ? cv_diff3(&a, &c) : e;
		}
		sprintf(msg, "extract, split and reverse of the borrowed curve as of the owned one: %.2g", e);
		TEST_ASSERT(e <= CV_TOL, msg);
		printf("  %s\n", msg);
	}
	qaws_curve_destroy(eo);
	qaws_curve_destroy(eb);
	qaws_curve_destroy(l);
	qaws_curve_destroy(r);
	qaws_curve_destroy(rv);
	qaws_curve_destroy(borrowed);
	qaws_curve_destroy(owned);
}

/* worst difference over 61 parameters, the prepared view and the borrowed curve against an owned curve (3D) */
static double cv_family_diff(qaws_curve_view const* v, qaws_curve const* owned)
{
	qaws_curve_prepared* p = NULL;
	qaws_curve* b = NULL;
	qaws_range r = qaws_curve_get_parameter_range(owned);
	double worst = 0, e;
	unsigned int i;
	if (qaws_curve_prepare(v, NULL, &p) != QAWS_STATUS_OK || qaws_curve_create_borrowed(v, NULL, &b) != QAWS_STATUS_OK)
	{
		qaws_curve_prepared_destroy(p);
		return 1e30;
	}
	for (i = 0; i <= 60; i++)
	{
		qaws_scalar t = r.min_value + (r.max_value - r.min_value) * (qaws_scalar)i / 60;
		qaws_eval_result_3d a, c, o;
		qaws_curve_prepared_evaluate_3d(p, t, CV_ALL, &a);
		qaws_curve_evaluate_3d(b, t, CV_ALL, &c);
		qaws_curve_evaluate_3d(owned, t, CV_ALL, &o);
		e = cv_diff3(&a, &o) > cv_diff3(&c, &o) ? cv_diff3(&a, &o) : cv_diff3(&c, &o);
		worst = e > worst ? e : worst;
		if (o.valid_flags & ~a.valid_flags)
			worst = 1e30;
	}
	qaws_curve_prepared_destroy(p);
	qaws_curve_destroy(b);
	return worst;
}

static void test_families(void)
{
	double p[8 * 3], tg[8 * 3], cf[5 * 3], rng[2] = { -1, 2 };
	qaws_scalar ps[8 * 3], ts[8 * 3], cs[5 * 3];
	qaws_curve_view v;
	qaws_curve* owned = NULL;
	double eh, eh2, ecr[4], ep;
	unsigned int i, k;
	char msg[300];
	cv_points(p, 8, 0.4);
	for (i = 0; i < 8 * 3; i++)
	{
		tg[i] = 0.7 * cos(0.37 * i + 1.0);
		ps[i] = (qaws_scalar)p[i];
		ts[i] = (qaws_scalar)tg[i];
	}
	for (i = 0; i < 5 * 3; i++)
		cs[i] = (qaws_scalar)(cf[i] = 0.4 * sin(1.3 * i + 0.1));
	/* Hermite */
	{
		qaws_hermite_desc d;
		memset(&d, 0, sizeof d);
		d.dimension = QAWS_DIMENSION_3D;
		d.degree = 3;
		d.points = ps;
		d.derivatives = ts;
		d.point_count = 8;
		d.derivative_count = 8;
		qaws_curve_create_hermite(&d, &owned);
		memset(&v, 0, sizeof v);
		v.family = QAWS_VIEW_HERMITE;
		v.dimension = QAWS_DIMENSION_3D;
		v.degree = 3;
		v.points = cv_packed(p, QAWS_ELEMENT_F64, 3, 8, 0);
		v.tangents = cv_packed(tg, QAWS_ELEMENT_F64, 3, 8, 0);
		eh = owned ? cv_family_diff(&v, owned) : 1e30;
		qaws_curve_destroy(owned);
		owned = NULL;
		/* the tangents rewritten in place: a prepared object follows them */
		for (i = 0; i < 8 * 3; i++)
			ts[i] = (qaws_scalar)(tg[i] = -0.5 * tg[i] + 0.1);
		qaws_curve_create_hermite(&d, &owned);
		eh2 = owned ? cv_family_diff(&v, owned) : 1e30;
		qaws_curve_destroy(owned);
		owned = NULL;
	}
	/* Catmull-Rom: open uniform, chordal, centripetal; closed centripetal */
	for (k = 0; k < 4; k++)
	{
		qaws_catmull_rom_desc d;
		memset(&d, 0, sizeof d);
		d.dimension = QAWS_DIMENSION_3D;
		d.control_points = ps;
		d.control_point_count = 8;
		d.parameterization = k == 0 ? QAWS_PARAMETERIZATION_UNIFORM : k == 1 ? QAWS_PARAMETERIZATION_CHORDAL : QAWS_PARAMETERIZATION_CENTRIPETAL;
		d.closed = k == 3;
		for (i = 0; i < 8 * 3; i++)
			ps[i] = (qaws_scalar)p[i];
		qaws_curve_create_catmull_rom(&d, &owned);
		memset(&v, 0, sizeof v);
		v.family = QAWS_VIEW_CATMULL_ROM;
		v.dimension = QAWS_DIMENSION_3D;
		v.degree = 3;
		v.points = cv_packed(p, QAWS_ELEMENT_F64, 3, 8, 0);
		v.parameterization = d.parameterization;
		v.closed = d.closed;
		ecr[k] = owned ? cv_family_diff(&v, owned) : 1e30;
		qaws_curve_destroy(owned);
		owned = NULL;
	}
	/* polynomial of degree 4 over [-1, 2] */
	{
		qaws_polynomial_desc d;
		memset(&d, 0, sizeof d);
		d.dimension = QAWS_DIMENSION_3D;
		d.degree = 4;
		d.coefficients = cs;
		d.coefficient_count = 5;
		d.t_min = -1;
		d.t_max = 2;
		qaws_curve_create_polynomial(&d, &owned);
		memset(&v, 0, sizeof v);
		v.family = QAWS_VIEW_POLYNOMIAL;
		v.dimension = QAWS_DIMENSION_3D;
		v.degree = 4;
		v.points = cv_packed(cf, QAWS_ELEMENT_F64, 3, 5, 0);
		v.knots = cv_packed(rng, QAWS_ELEMENT_F64, 1, 2, 0);
		ep = owned ? cv_family_diff(&v, owned) : 1e30;
		qaws_curve_destroy(owned);
		owned = NULL;
	}
	sprintf(msg, "families as views against the owned curves: Hermite %.2g (tangents rewritten %.2g), Catmull-Rom uniform %.2g, chordal %.2g, "
		"centripetal %.2g, closed %.2g, polynomial %.2g", eh, eh2, ecr[0], ecr[1], ecr[2], ecr[3], ep);
	TEST_ASSERT(eh <= CV_TOL && eh2 <= CV_TOL && ecr[0] <= CV_TOL && ecr[1] <= CV_TOL && ecr[2] <= CV_TOL && ecr[3] <= CV_TOL && ep <= CV_TOL, msg);
	printf("  %s\n", msg);
	/* refusals */
	memset(&v, 0, sizeof v);
	v.family = QAWS_VIEW_HERMITE;
	v.dimension = QAWS_DIMENSION_3D;
	v.degree = 3;
	v.points = cv_packed(p, QAWS_ELEMENT_F64, 3, 8, 0);
	v.tangents = cv_packed(tg, QAWS_ELEMENT_F64, 3, 7, 0);
	TEST_ASSERT(qaws_curve_view_validate(&v) == QAWS_STATUS_INVALID_CONTROL_POINT_COUNT, "Hermite: tangents not one per point refused");
	v.tangents.count = 8;
	v.degree = 2;
	TEST_ASSERT(qaws_curve_view_validate(&v) == QAWS_STATUS_INVALID_DEGREE, "Hermite: degree 2 refused");
	v.family = QAWS_VIEW_CATMULL_ROM;
	v.degree = 3;
	v.points.count = 3;
	TEST_ASSERT(qaws_curve_view_validate(&v) == QAWS_STATUS_INVALID_CONTROL_POINT_COUNT, "Catmull-Rom: 3 points open refused");
	v.closed = 1;
	TEST_ASSERT(qaws_curve_view_validate(&v) == QAWS_STATUS_OK, "Catmull-Rom: 3 points closed accepted");
	v.family = QAWS_VIEW_POLYNOMIAL;
	v.closed = 0;
	v.degree = 2;
	rng[0] = 3;
	v.knots = cv_packed(rng, QAWS_ELEMENT_F64, 1, 2, 0);
	TEST_ASSERT(qaws_curve_view_validate(&v) == QAWS_STATUS_INVALID_PARAMETER_RANGE, "polynomial: an empty range refused");
	/* operations through owned copies: extract of a borrowed Catmull-Rom curve as of the owned */
	{
		qaws_catmull_rom_desc d;
		qaws_curve *b = NULL, *eo = NULL, *eb = NULL;
		double e = 0;
		memset(&d, 0, sizeof d);
		d.dimension = QAWS_DIMENSION_3D;
		d.control_points = ps;
		d.control_point_count = 8;
		d.parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL;
		qaws_curve_create_catmull_rom(&d, &owned);
		memset(&v, 0, sizeof v);
		v.family = QAWS_VIEW_CATMULL_ROM;
		v.dimension = QAWS_DIMENSION_3D;
		v.degree = 3;
		v.points = cv_packed(p, QAWS_ELEMENT_F64, 3, 8, 0);
		v.parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL;
		if (owned && qaws_curve_create_borrowed(&v, NULL, &b) == QAWS_STATUS_OK && qaws_curve_extract(owned, (qaws_scalar)0.6, (qaws_scalar)3.7, &eo) == QAWS_STATUS_OK
			&& qaws_curve_extract(b, (qaws_scalar)0.6, (qaws_scalar)3.7, &eb) == QAWS_STATUS_OK)
		{
			for (i = 0; i <= 20; i++)
			{
				qaws_scalar t = (qaws_scalar)(0.6 + 3.1 * i / 20);
				qaws_eval_result_3d a, c;
				qaws_curve_evaluate_3d(eo, t, QAWS_EVAL_FLAG_POSITION, &a);
				qaws_curve_evaluate_3d(eb, t, QAWS_EVAL_FLAG_POSITION, &c);
				e = cv_diff3(&a, &c) > e ? cv_diff3(&a, &c) : e;
			}
			sprintf(msg, "extract of a borrowed Catmull-Rom curve as of the owned: %.2g", e);
			TEST_ASSERT(e <= CV_TOL, msg);
		}
		else
			TEST_ASSERT(0, "Catmull-Rom extracts");
		qaws_curve_destroy(eo);
		qaws_curve_destroy(eb);
		qaws_curve_destroy(b);
		qaws_curve_destroy(owned);
	}
}

int test_142_curve_view_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 142: borrowed curves\n");
	test_parity_bspline();
	test_parity_nurbs();
	test_layouts();
	test_mutation();
	test_ranges();
	test_allocation();
	test_validation();
	test_operations();
	test_families();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
