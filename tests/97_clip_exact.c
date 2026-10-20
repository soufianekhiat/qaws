/*
 * Test 97: exact Boolean operations on curves (qaws_clip_exact)
 *
 *   - two squares, every clip type: areas exact, the float engine's paths
 *   - squares sharing an edge (an overlap merged) and touching at a corner
 *   - a parabola touching a square's side at an exact point (curvature
 *     decides the order round the vertex): union, intersection, difference
 *   - three segments through one rational point (1/3, 1/3): one vertex
 *   - a square with a hole and an island in it: three paths nested 0, 1, 2
 *   - a pentagram, even-odd and non-zero
 *   - crossing ellipses (NURBS conics) with a duplicate: the float engine's
 *     area and paths
 */

#include "test_common.h"
#include "qaws_clip_exact.h"
#include <math.h>
#include <string.h>

typedef struct cxt_scene
{
	qaws_curve* curves[64];
	qaws_curve const* views[64];
	qaws_path_2d paths[32];
	unsigned int nc, np;
} cxt_scene;

static qaws_path_2d* cxt_path(cxt_scene* s, qaws_curve* c)
{
	qaws_path_2d* p = &s->paths[s->np++];
	s->curves[s->nc] = c;
	s->views[s->nc] = c;
	p->curves = &s->views[s->nc];
	p->curve_count = 1;
	p->closed = 1;
	s->nc++;
	return p;
}

static qaws_curve* cxt_poly(double const* xy, unsigned int n)
{
	qaws_scalar p[64];
	qaws_curve* c = NULL;
	unsigned int i;
	for (i = 0; i < 2 * n; i++)
		p[i] = (qaws_scalar)xy[i];
	qaws_curve_create_polyline_2d(&c, p, n, 1);
	return c;
}

static qaws_curve* cxt_rect(double x0, double y0, double x1, double y1)
{
	double xy[8];
	xy[0] = x0; xy[1] = y0; xy[2] = x1; xy[3] = y0; xy[4] = x1; xy[5] = y1; xy[6] = x0; xy[7] = y1;
	return cxt_poly(xy, 4);
}

static void cxt_free(cxt_scene* s)
{
	unsigned int i;
	for (i = 0; i < s->nc; i++)
		qaws_curve_destroy(s->curves[i]);
}

static double cxt_area(qaws_clip_result const* r)
{
	unsigned int i, n = r ? qaws_clip_result_get_path_count(r) : 0;
	double a = 0;
	for (i = 0; i < n; i++)
	{
		qaws_path_2d p;
		qaws_scalar s = 0;
		qaws_clip_result_get_path(&p, r, i);
		qaws_path_compute_area_2d(&s, &p);
		a += (double)s;
	}
	return a;
}

typedef struct cxt_run
{
	qaws_status st, fst;
	unsigned int paths, fpaths;
	double area, farea;
	qaws_clip_exact_stats stats;
	qaws_clip_result* r;
} cxt_run;

/* the exact engine, and the float one on the same input */
static void cxt_do(cxt_run* o, qaws_path_2d const* subj, unsigned int ns, qaws_path_2d const* clip, unsigned int nc, qaws_clip_type ct, qaws_fill_rule fr)
{
	qaws_clip_exact_desc d;
	qaws_clip_desc fd;
	qaws_clip_result* fr_res = NULL;
	memset(&d, 0, sizeof(d));
	d.subjects = subj;
	d.subject_count = ns;
	d.clips = clip;
	d.clip_count = nc;
	d.clip_type = ct;
	d.fill_rule = fr;
	o->r = NULL;
	o->st = qaws_clip_exact_execute(&o->r, &o->stats, &d);
	o->paths = o->r ? qaws_clip_result_get_path_count(o->r) : 0;
	o->area = cxt_area(o->r);
	memset(&fd, 0, sizeof(fd));
	fd.subjects = subj;
	fd.subject_count = ns;
	fd.clips = clip;
	fd.clip_count = nc;
	fd.clip_type = ct;
	fd.fill_rule = fr;
	o->fst = qaws_clip_execute(&fr_res, &fd);
	o->fpaths = fr_res ? qaws_clip_result_get_path_count(fr_res) : 0;
	o->farea = cxt_area(fr_res);
	qaws_clip_result_destroy(fr_res);
}

static void cxt_msg(char* msg, char const* name, cxt_run const* o)
{
	sprintf(msg, "%s: status %d, %u paths, area %.12g (float engine: status %d, %u paths, area %.12g); %u vertices, %u edges, %u faces, %u parts, %u overlaps, %u tangencies, %u identity tests, %u ray retries",
		name, (int)o->st, o->paths, o->area, (int)o->fst, o->fpaths, o->farea, o->stats.vertex_count, o->stats.edge_count, o->stats.face_count,
		o->stats.component_count, o->stats.overlap_count, o->stats.tangency_count, o->stats.identity_tests, o->stats.ray_retries);
}

static void test_squares(void)
{
	static char const* names[4] = { "intersection", "union", "difference", "xor" };
	static double const want[4] = { 1, 7, 3, 6 };
	cxt_scene s;
	qaws_path_2d *a, *b;
	unsigned int k;
	char msg[600];
	memset(&s, 0, sizeof(s));
	a = cxt_path(&s, cxt_rect(0, 0, 2, 2));
	b = cxt_path(&s, cxt_rect(1, 1, 3, 3));
	for (k = 0; k < 4; k++)
	{
		cxt_run o;
		char name[64];
		cxt_do(&o, a, 1, b, 1, (qaws_clip_type)(QAWS_CLIP_INTERSECTION + k), QAWS_FILL_NON_ZERO);
		sprintf(name, "two squares, %s", names[k]);
		cxt_msg(msg, name, &o);
		TEST_ASSERT(o.st == QAWS_STATUS_OK && fabs(o.area - want[k]) < 1e-12 && o.paths == o.fpaths, msg);
		qaws_clip_result_destroy(o.r);
	}
	cxt_free(&s);
}

static void test_shared(void)
{
	char msg[600];
	cxt_scene s;
	qaws_path_2d *a, *b, *c;
	cxt_run o;
	memset(&s, 0, sizeof(s));
	a = cxt_path(&s, cxt_rect(0, 0, 1, 1));
	b = cxt_path(&s, cxt_rect(1, 0, 2, 1));
	c = cxt_path(&s, cxt_rect(1, 1, 2, 2));
	(void)c;
	/* a and b share the edge x = 1: one square of area 2 */
	cxt_do(&o, a, 1, b, 1, QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO);
	cxt_msg(msg, "squares sharing an edge, union", &o);
	TEST_ASSERT(o.st == QAWS_STATUS_OK && o.paths == 1 && fabs(o.area - 2) < 1e-12 && o.stats.overlap_count >= 1, msg);
	qaws_clip_result_destroy(o.r);
	/* a and c touch at the corner (1, 1): two paths, area 2 */
	cxt_do(&o, a, 1, &s.paths[2], 1, QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO);
	cxt_msg(msg, "squares touching at a corner, union", &o);
	TEST_ASSERT(o.st == QAWS_STATUS_OK && o.paths == 2 && fabs(o.area - 2) < 1e-12, msg);
	qaws_clip_result_destroy(o.r);
	cxt_do(&o, a, 1, &s.paths[2], 1, QAWS_CLIP_INTERSECTION, QAWS_FILL_NON_ZERO);
	cxt_msg(msg, "squares touching at a corner, intersection", &o);
	TEST_ASSERT(o.st == QAWS_STATUS_OK && o.paths == 0 && fabs(o.area) < 1e-12, msg);
	qaws_clip_result_destroy(o.r);
	cxt_free(&s);
}

static void test_tangent(void)
{
	/* the lens above y = (x - 1)^2 and below y = 1 (a parabola closed by a line), in the square [0, 2]^2: it touches the bottom side at (1, 0) */
	static double const par[6] = { 0, 1, 1, -1, 2, 1 };
	static char const* names[3] = { "union", "intersection", "difference" };
	static qaws_clip_type const ct[3] = { QAWS_CLIP_UNION, QAWS_CLIP_INTERSECTION, QAWS_CLIP_DIFFERENCE };
	double const want[3] = { 4, 4.0 / 3.0, 8.0 / 3.0 };
	char msg[600];
	cxt_scene s;
	qaws_path_2d *sq, *lens;
	qaws_curve* pc = NULL;
	qaws_bezier_desc bd;
	qaws_scalar cp[6];
	unsigned int k;
	memset(&s, 0, sizeof(s));
	sq = cxt_path(&s, cxt_rect(0, 0, 2, 2));
	for (k = 0; k < 6; k++) cp[k] = (qaws_scalar)par[k];
	memset(&bd, 0, sizeof(bd));
	bd.dimension = QAWS_DIMENSION_2D;
	bd.degree = 2;
	bd.control_points = cp;
	bd.control_point_count = 3;
	qaws_curve_create_bezier(&pc, &bd);
	lens = cxt_path(&s, pc);   /* the path closes with the line back from (2, 1) to (0, 1) */
	for (k = 0; k < 3; k++)
	{
		cxt_run o;
		char name[80];
		cxt_do(&o, sq, 1, lens, 1, ct[k], QAWS_FILL_NON_ZERO);
		sprintf(name, "parabola touching a square's side, %s", names[k]);
		cxt_msg(msg, name, &o);
		TEST_ASSERT(o.st == QAWS_STATUS_OK && fabs(o.area - want[k]) < (sizeof(qaws_scalar) == 4 ? 1e-6 : 1e-9) && o.stats.tangency_count >= 1, msg);
		qaws_clip_result_destroy(o.r);
	}
	cxt_free(&s);
}

static void test_concurrent(void)
{
	/* y = x, y = (1 - x) / 2 and y = 1 - 2x all pass through (1/3, 1/3) */
	static double const t1[6] = { 0, 0, 1, 1, -1, 1 };
	static double const t2[6] = { 1, 0, 0, 0.5, 1, 1 };
	static double const t3[6] = { 0, 1, 0.5, 0, 0, 0 };
	char msg[600];
	cxt_scene s;
	qaws_path_2d* p;
	cxt_run o;
	memset(&s, 0, sizeof(s));
	p = cxt_path(&s, cxt_poly(t1, 3));
	cxt_path(&s, cxt_poly(t2, 3));
	cxt_path(&s, cxt_poly(t3, 3));
	cxt_do(&o, p, 3, NULL, 0, QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO);
	cxt_msg(msg, "three triangles with edges through (1/3, 1/3), union", &o);
	TEST_ASSERT(o.st == QAWS_STATUS_OK && fabs(o.area - o.farea) < 1e-9 && o.paths == o.fpaths, msg);
	qaws_clip_result_destroy(o.r);
	cxt_do(&o, p, 3, NULL, 0, QAWS_CLIP_UNION, QAWS_FILL_EVEN_ODD);
	cxt_msg(msg, "three triangles with edges through (1/3, 1/3), even-odd", &o);
	TEST_ASSERT(o.st == QAWS_STATUS_OK && fabs(o.area - o.farea) < 1e-9 && o.paths == o.fpaths, msg);
	qaws_clip_result_destroy(o.r);
	cxt_free(&s);
}

static void test_nested(void)
{
	char msg[600];
	cxt_scene s;
	qaws_path_2d* p;
	cxt_run o;
	qaws_curve* hole;
	double hxy[8] = { 1, 1, 1, 5, 5, 5, 5, 1 };   /* clockwise */
	memset(&s, 0, sizeof(s));
	p = cxt_path(&s, cxt_rect(0, 0, 6, 6));
	hole = cxt_poly(hxy, 4);
	cxt_path(&s, hole);
	cxt_path(&s, cxt_rect(2, 2, 4, 4));
	cxt_do(&o, p, 3, NULL, 0, QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO);
	cxt_msg(msg, "a square with a hole and an island", &o);
	{
		unsigned int i, depth[3] = { 9, 9, 9 };
		int ok = o.st == QAWS_STATUS_OK && o.paths == 3;
		for (i = 0; ok && i < 3; i++)
			depth[i] = qaws_clip_result_get_depth(o.r, i);
		ok = ok && fabs(o.area - (36 - 16 + 4)) < 1e-12 && depth[0] + depth[1] + depth[2] == 3;
		TEST_ASSERT(ok, msg);
	}
	qaws_clip_result_destroy(o.r);
	cxt_free(&s);
}

static void test_star(void)
{
	char msg[600];
	cxt_scene s;
	qaws_path_2d* p;
	double xy[10];
	unsigned int k;
	cxt_run o;
	memset(&s, 0, sizeof(s));
	for (k = 0; k < 5; k++)
	{
		/* a pentagram on dyadic points (roughly round) */
		double a = 3.14159265358979 / 2 + k * 4 * 3.14159265358979 / 5;
		xy[2 * k] = floor(cos(a) * 1024) / 256;
		xy[2 * k + 1] = floor(sin(a) * 1024) / 256;
	}
	p = cxt_path(&s, cxt_poly(xy, 5));
	cxt_do(&o, p, 1, NULL, 0, QAWS_CLIP_UNION, QAWS_FILL_EVEN_ODD);
	cxt_msg(msg, "pentagram, even-odd", &o);
	TEST_ASSERT(o.st == QAWS_STATUS_OK && fabs(o.area - o.farea) < 1e-9 && o.paths == o.fpaths && o.paths == 5, msg);
	qaws_clip_result_destroy(o.r);
	cxt_do(&o, p, 1, NULL, 0, QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO);
	cxt_msg(msg, "pentagram, non-zero", &o);
	TEST_ASSERT(o.st == QAWS_STATUS_OK && fabs(o.area - o.farea) < 1e-9 && o.paths == o.fpaths && o.paths == 1, msg);
	qaws_clip_result_destroy(o.r);
	cxt_free(&s);
}

static void test_ellipses(void)
{
	char msg[600];
	cxt_scene s;
	qaws_path_2d* p;
	qaws_curve* e = NULL;
	unsigned int k;
	cxt_run o;
	memset(&s, 0, sizeof(s));
	for (k = 0; k < 4; k++)
	{
		e = NULL;
		qaws_curve_create_ellipse_2d(&e, (qaws_scalar)(0.75 * k), (qaws_scalar)(0.25 * k), (qaws_scalar)1.5, (qaws_scalar)1.0, (qaws_scalar)(0.3 * k));
		cxt_path(&s, e);
	}
	/* the second one again */
	e = NULL;
	qaws_curve_create_ellipse_2d(&e, (qaws_scalar)0.75, (qaws_scalar)0.25, (qaws_scalar)1.5, (qaws_scalar)1.0, (qaws_scalar)0.3);
	cxt_path(&s, e);
	p = &s.paths[0];
	cxt_do(&o, p, 4, NULL, 0, QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO);
	cxt_msg(msg, "four crossing ellipses, union", &o);
	TEST_ASSERT(o.st == QAWS_STATUS_OK && fabs(o.area - o.farea) < 1e-6 * fabs(o.farea) && o.paths == o.fpaths, msg);
	qaws_clip_result_destroy(o.r);
	cxt_do(&o, p, 5, NULL, 0, QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO);
	cxt_msg(msg, "four crossing ellipses and a duplicate, union", &o);
	TEST_ASSERT(o.st == QAWS_STATUS_OK && fabs(o.area - o.farea) < 1e-6 * fabs(o.farea) && o.stats.overlap_count >= 1, msg);
	qaws_clip_result_destroy(o.r);
	cxt_do(&o, p, 2, &s.paths[2], 2, QAWS_CLIP_XOR, QAWS_FILL_NON_ZERO);
	cxt_msg(msg, "two ellipses xor two ellipses", &o);
	TEST_ASSERT(o.st == QAWS_STATUS_OK && fabs(o.area - o.farea) < 1e-6 * fabs(o.farea) && o.paths == o.fpaths, msg);
	qaws_clip_result_destroy(o.r);
	cxt_free(&s);
}

int test_97_clip_exact_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 97: exact Boolean operations on curves\n");
	test_squares();
	test_shared();
	test_tangent();
	test_concurrent();
	test_nested();
	test_star();
	test_ellipses();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
