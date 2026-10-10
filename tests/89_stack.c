/*
 * Test 89: 2.5D stacks of parallel contours (qaws_stack)
 *
 *   - a cone: a circle at z 0, a point at z 4; the section at z 1 is the
 *     circle of radius 1.5 exactly (matched contours, the point promoted)
 *   - a sphere from 17 circles at uneven heights: the cubic blend follows
 *     sqrt(1 - z^2) closer than the linear one
 *   - a peanut splitting into two blobs between levels (different contour
 *     counts): the distance-field section lies between the two levels
 *   - outside the stack: empty; at a level: that level's region
 *   - cone x sphere: one level per height of either stack, each the
 *     intersection of the sections, the same serially and on an executor
 */

#include "test_common.h"
#include <math.h>
#include <string.h>

#define STT_PI 3.14159265358979323846
/* areas agree to 1e-9 in double, to float rounding in single */
#define STT_TOL(ref) (sizeof(qaws_scalar) == 4 ? 2e-6 * (fabs(ref) + 1) : 1e-9)

static double stt_area(qaws_clip_result const* r)
{
	unsigned int i, n = qaws_clip_result_get_path_count(r);
	double s = 0.0;
	for (i = 0; i < n; i++)
	{
		qaws_path_2d p;
		qaws_scalar a = 0;
		qaws_clip_result_get_path(r, i, &p);
		qaws_path_compute_area_2d(&p, &a);
		s += a;
	}
	return s;
}

typedef struct stt_stack
{
	qaws_curve* curves[64];
	qaws_curve const* views[64];
	qaws_path_2d paths[64];
	qaws_stack_level levels[64];
	unsigned int n;
	qaws_stack_2d st;
} stt_stack;

static void stt_add(stt_stack* s, double z, qaws_curve* c)
{
	unsigned int i = s->n++;
	s->curves[i] = c;
	s->views[i] = c;
	s->paths[i].curves = &s->views[i];
	s->paths[i].curve_count = 1;
	s->paths[i].closed = 1;
	s->levels[i].z = (qaws_scalar)z;
	s->levels[i].paths = &s->paths[i];
	s->levels[i].path_count = 1;
}

static void stt_done(stt_stack* s, qaws_stack_interp interp)
{
	memset(&s->st, 0, sizeof(s->st));
	s->st.levels = s->levels;
	s->st.level_count = s->n;
	s->st.fill_rule = QAWS_FILL_NON_ZERO;
	s->st.interp = interp;
}

static void stt_free(stt_stack* s)
{
	unsigned int i;
	for (i = 0; i < s->n; i++) qaws_curve_destroy(s->curves[i]);
}

/* a point contour: a closed polyline of two equal points */
static qaws_curve* stt_point(double x, double y)
{
	qaws_scalar p[4];
	qaws_curve* c = NULL;
	p[0] = p[2] = (qaws_scalar)x;
	p[1] = p[3] = (qaws_scalar)y;
	qaws_curve_create_polyline_2d(p, 2, 1, &c);
	return c;
}

static qaws_curve* stt_circle(double cx, double cy, double r)
{
	qaws_curve* c = NULL;
	qaws_curve_create_ellipse_2d((qaws_scalar)cx, (qaws_scalar)cy, (qaws_scalar)r, (qaws_scalar)r, 0, &c);
	return c;
}

static void test_cone(void)
{
	stt_stack s;
	qaws_clip_result* r = NULL;
	char msg[256];
	double a;
	memset(&s, 0, sizeof(s));
	stt_add(&s, 0, stt_circle(0, 0, 2));
	stt_add(&s, 4, stt_point(0, 0));
	stt_done(&s, QAWS_STACK_LINEAR);
	qaws_stack_section_2d(&s.st, 1, &r);
	a = r ? stt_area(r) : 0;
	sprintf(msg, "cone: the section at z 1 is the circle of radius 1.5 (area %.12g vs %.12g)", a, STT_PI * 2.25);
	TEST_ASSERT(r && fabs(a - STT_PI * 2.25) < STT_TOL(STT_PI * 2.25), msg);
	if (r)
	{
		qaws_path_2d p;
		qaws_clip_result_get_path(r, 0, &p);
		TEST_ASSERT(p.curve_count >= 1 && qaws_curve_get_kind(p.curves[0]) == QAWS_CURVE_KIND_NURBS, "the cone's section is a NURBS circle, not a fit");
	}
	qaws_clip_result_destroy(r);
	r = NULL;
	qaws_stack_section_2d(&s.st, 5, &r);
	TEST_ASSERT(r && qaws_clip_result_get_path_count(r) == 0, "above the apex: empty");
	qaws_clip_result_destroy(r);
	r = NULL;
	qaws_stack_section_2d(&s.st, 0, &r);
	TEST_ASSERT(r && fabs(stt_area(r) - 4 * STT_PI) < STT_TOL(4 * STT_PI), "at a level: that level's region");
	qaws_clip_result_destroy(r);
	stt_free(&s);
}

static void test_sphere(void)
{
	stt_stack s;
	unsigned int i;
	double el = 0, ec = 0;
	char msg[256];
	memset(&s, 0, sizeof(s));
	/* 17 circles at heights cos(theta): denser near the equator's edge */
	for (i = 0; i <= 16; i++)
	{
		double th = STT_PI * (16 - i) / 16, z = cos(th), r = sin(th);
		stt_add(&s, z, r > 1e-9 ? stt_circle(0, 0, r) : stt_point(0, 0));
	}
	/* sections between levels: radius against sqrt(1 - z^2) */
	for (i = 0; i < 2; i++)
	{
		unsigned int k;
		stt_done(&s, i ? QAWS_STACK_CUBIC : QAWS_STACK_LINEAR);
		for (k = 2; k < 14; k++)   /* away from the poles, where sqrt(1 - z^2) has no slope a blend can follow */
		{
			double z = 0.5 * ((double)s.levels[k].z + (double)s.levels[k + 1].z), want = sqrt(1 - z * z), got, e;
			qaws_clip_result* r = NULL;
			qaws_stack_section_2d(&s.st, (qaws_scalar)z, &r);
			got = r ? sqrt(stt_area(r) / STT_PI) : 0;
			e = fabs(got - want);
			if (i) { if (e > ec) ec = e; } else { if (e > el) el = e; }
			qaws_clip_result_destroy(r);
		}
	}
	sprintf(msg, "sphere from 17 circles: cubic blend closer to sqrt(1 - z^2) than linear away from the poles (worst radius error %.2e vs %.2e)", ec, el);
	TEST_ASSERT(ec < el && ec < 0.25 * el, msg);
	stt_free(&s);
}

static void test_split(void)
{
	/* a long ellipse at z 0 splitting into two circles at z 1 */
	qaws_curve* e = NULL;
	qaws_curve *c1 = stt_circle(-2, 0, 1), *c2 = stt_circle(2, 0, 1);
	qaws_curve const* v[3];
	qaws_path_2d p[3];
	qaws_stack_level lv[2];
	qaws_stack_2d st;
	qaws_clip_result* r = NULL;
	double a;
	char msg[256];
	qaws_curve_create_ellipse_2d(0, 0, 3, (qaws_scalar)1.2, 0, &e);
	v[0] = e; v[1] = c1; v[2] = c2;
	p[0].curves = &v[0]; p[0].curve_count = 1; p[0].closed = 1;
	p[1].curves = &v[1]; p[1].curve_count = 1; p[1].closed = 1;
	p[2].curves = &v[2]; p[2].curve_count = 1; p[2].closed = 1;
	lv[0].z = 0; lv[0].paths = &p[0]; lv[0].path_count = 1;
	lv[1].z = 1; lv[1].paths = &p[1]; lv[1].path_count = 2;
	memset(&st, 0, sizeof(st));
	st.levels = lv; st.level_count = 2; st.fill_rule = QAWS_FILL_NON_ZERO;
	qaws_stack_section_2d(&st, (qaws_scalar)0.9, &r);
	a = r ? stt_area(r) : 0;
	sprintf(msg, "ellipse splitting into two circles: the distance section at z 0.9 has 2 parts (%u) and an area between the levels' (%.6g in [%.6g, %.6g])",
		r ? qaws_clip_result_get_path_count(r) : 0, a, 2 * STT_PI, STT_PI * 3.6);
	TEST_ASSERT(r && qaws_clip_result_get_path_count(r) == 2 && a > 2 * STT_PI * 0.9 && a < STT_PI * 3.6, msg);
	qaws_clip_result_destroy(r);
	qaws_curve_destroy(e); qaws_curve_destroy(c1); qaws_curve_destroy(c2);
}

/* a serial executor that runs chunks backwards (results must not care) */
static void stt_backwards(void* user, unsigned int count, qaws_batch_task_fn task, void* ctx)
{
	unsigned int i;
	(void)user;
	for (i = count; i-- > 0;)
		task(ctx, i, i + 1);
}

static void test_boolean(void)
{
	stt_stack cone, sphere;
	qaws_stack_result *r = NULL, *r2 = NULL;
	qaws_batch_executor ex;
	unsigned int i, n, same = 1;
	char msg[256];
	memset(&cone, 0, sizeof(cone));
	memset(&sphere, 0, sizeof(sphere));
	stt_add(&cone, 0, stt_circle(0, 0, 1.5));
	stt_add(&cone, 2, stt_point(0, 0));
	stt_done(&cone, QAWS_STACK_LINEAR);
	for (i = 0; i <= 12; i++)
	{
		double th = STT_PI * (12 - i) / 12, z = 1 + cos(th), rr = sin(th);
		stt_add(&sphere, z, rr > 1e-9 ? stt_circle(0.5, 0, rr) : stt_point(0.5, 0));
	}
	stt_done(&sphere, QAWS_STACK_LINEAR);
	qaws_stack_boolean_2d(QAWS_CLIP_INTERSECTION, &cone.st, &sphere.st, NULL, &r);
	memset(&ex, 0, sizeof(ex));
	ex.parallel_for = stt_backwards;
	qaws_stack_boolean_2d(QAWS_CLIP_INTERSECTION, &cone.st, &sphere.st, &ex, &r2);
	n = r ? qaws_stack_result_get_level_count(r) : 0;
	for (i = 0; r && r2 && i < n; i++)
		same = same && fabs(stt_area(qaws_stack_result_get_level(r, i)) - stt_area(qaws_stack_result_get_level(r2, i))) == 0;
	sprintf(msg, "cone x sphere: one level per height of either stack (%u: 13 sphere heights, the cone's 0 and 2 among them), the same on an executor (%d)", n, same);
	TEST_ASSERT(n == 13 && same, msg);
	if (r)
	{
		/* at z = 1 (a sphere level): cone radius 0.75 at (0, 0), sphere radius 1 at (0.5, 0) */
		for (i = 0; i < n; i++)
			if (fabs((double)qaws_stack_result_get_z(r, i) - 1) < 1e-12)
			{
				double r1 = 0.75, r2 = 1, dd = 0.5;
				double lens = r1 * r1 * acos((dd * dd + r1 * r1 - r2 * r2) / (2 * dd * r1)) + r2 * r2 * acos((dd * dd + r2 * r2 - r1 * r1) / (2 * dd * r2))
					- 0.5 * sqrt((-dd + r1 + r2) * (dd + r1 - r2) * (dd - r1 + r2) * (dd + r1 + r2));
				double a = stt_area(qaws_stack_result_get_level(r, i));
				sprintf(msg, "at z 1 the intersection is the lens of the cone's r 0.75 disc and the sphere's r 1 disc (%.10g vs %.10g)", a, lens);
				TEST_ASSERT(fabs(a - lens) < STT_TOL(lens), msg);
			}
	}
	qaws_stack_result_destroy(r);
	qaws_stack_result_destroy(r2);
	stt_free(&cone);
	stt_free(&sphere);
}

int test_89_stack_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 89: 2.5D stacks of parallel contours\n");
	test_cone();
	test_sphere();
	test_split();
	test_boolean();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
