/*
 * Test 160: quality reports (qaws_*_measure_elements, qaws_quality_summarize)
 *
 *   - 2D: an equilateral triangle (angles 60, shape, aspect and scaled
 *     Jacobian 1), a right isosceles one (90 and 45 degrees, its known
 *     shape and aspect), the equilateral turned clockwise (size, shape and
 *     Jacobian negative, angles past 180), a square (Jacobian 1, angles 90,
 *     shape not applicable);
 *   - surfaces: the equilateral triangle tilted in space (shape 1), a flat
 *     triangle (shape 0, aspect infinite), a perturbed quad grid against
 *     qaws_mesh_3d_measure_shape (the same smallest angle and Jacobian);
 *   - volumes: a regular tetrahedron (shape, aspect, Jacobian 1, dihedral
 *     70.53), a unit cube (Jacobian 1, dihedral 90, volume 1), a regular
 *     prism and a pyramid of equal edges (Jacobian 1); a sheared block
 *     against qaws_volume_mesh_measure_quality (the same smallest
 *     Jacobian);
 *   - the summary: statistics, the median, a histogram over a given span
 *     and over the values' own, the worst elements in order for both
 *     directions, failing counts, NaN skipped;
 *   - refusals: a measure out of range, NULL outputs.
 */

#include "test_common.h"
#include "qaws_quality_report.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

static int qt_near(double a, double b, double tol)
{
	return fabs(a - b) <= tol * (1 + fabs(b));
}

static double qt_tol(void)
{
	return sizeof(qaws_scalar) == 4 ? 1e-5 : 1e-11;
}

static void test_2d(void)
{
	double h = sqrt(3.0) / 2;
	qaws_scalar x[] = { 0, 0, 1, 0, (qaws_scalar)0.5, (qaws_scalar)h, 2, 0, 3, 0, 2, 1, 4, 0, 5, 0, (qaws_scalar)4.5, (qaws_scalar)h, 6, 0, 7, 0, 7, 1, 6, 1 };
	unsigned int f[] = { 0, 1, 2, 3, 4, 5, 6, 8, 7, 9, 10, 11, 12 };
	unsigned int fs[] = { 3, 3, 3, 4 };
	qaws_scalar v[4];
	qaws_mesh_2d_desc d;
	qaws_mesh_2d* m = NULL;
	memset(&d, 0, sizeof d);
	d.positions = x;
	d.vertex_count = 13;
	d.face_vertices = f;
	d.face_sizes = fs;
	d.face_count = 4;
	/* (the third triangle, 4 0, 4.5 h, 5 0, is clockwise) */
	TEST_ASSERT(qaws_mesh_2d_create(&m, &d) == QAWS_STATUS_OK, "2D: mesh made");
	if (!m)
		return;
	TEST_ASSERT(qaws_mesh_2d_get_triangle_count(m) == 3 && qaws_mesh_2d_get_quad_count(m) == 1, "2D: three triangles and a quad");
	qaws_mesh_2d_measure_elements(v, m, QAWS_MEASURE_MIN_ANGLE);
	TEST_ASSERT(qt_near(v[0], 60, qt_tol()) && qt_near(v[1], 45, qt_tol()) && qt_near(v[3], 90, qt_tol()), "2D: smallest angles 60, 45, 90");
	qaws_mesh_2d_measure_elements(v, m, QAWS_MEASURE_MAX_ANGLE);
	TEST_ASSERT(qt_near(v[0], 60, qt_tol()) && qt_near(v[1], 90, qt_tol()) && v[2] > 180 && qt_near(v[3], 90, qt_tol()),
		"2D: largest angles 60, 90, past 180 when turned, 90");
	qaws_mesh_2d_measure_elements(v, m, QAWS_MEASURE_SHAPE);
	TEST_ASSERT(qt_near(v[0], 1, qt_tol()) && qt_near(v[1], sqrt(3.0) / 2, qt_tol()) && qt_near(v[2], -1, qt_tol()) && v[3] != v[3],
		"2D: shape 1, sqrt 3 / 2 for the right isosceles, -1 turned, NaN for the quad");
	qaws_mesh_2d_measure_elements(v, m, QAWS_MEASURE_ASPECT_RATIO);
	TEST_ASSERT(qt_near(v[0], 1, qt_tol()) && qt_near(v[1], 1 / (2 * (sqrt(2.0) - 1)), qt_tol()), "2D: aspect 1 and 1 / (2 (sqrt 2 - 1))");
	qaws_mesh_2d_measure_elements(v, m, QAWS_MEASURE_SCALED_JACOBIAN);
	TEST_ASSERT(qt_near(v[0], 1, qt_tol()) && qt_near(v[1], 1 / sqrt(2.0) / h, qt_tol()) && v[2] < 0 && qt_near(v[3], 1, qt_tol()),
		"2D: scaled Jacobian 1, sin 45 / sin 60, negative turned, 1 for the square");
	qaws_mesh_2d_measure_elements(v, m, QAWS_MEASURE_SIZE);
	TEST_ASSERT(qt_near(v[0], h / 2, qt_tol()) && qt_near(v[1], 0.5, qt_tol()) && qt_near(v[2], -h / 2, qt_tol()) && qt_near(v[3], 1, qt_tol()),
		"2D: areas, the turned one negative");
	qaws_mesh_2d_measure_elements(v, m, QAWS_MEASURE_EDGE_RATIO);
	TEST_ASSERT(qt_near(v[0], 1, qt_tol()) && qt_near(v[1], sqrt(2.0), qt_tol()) && qt_near(v[3], 1, qt_tol()), "2D: edge ratios");
	qaws_mesh_2d_destroy(m);
}

static void test_surface(void)
{
	/* the equilateral triangle tilted (rotated about x by 40 degrees, about z by 25) and a flat one */
	double h = sqrt(3.0) / 2, ca = cos(0.698), sa = sin(0.698), cb = cos(0.436), sb = sin(0.436), p[3][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 0.5, h, 0 } };
	qaws_scalar x[18], v[2];
	unsigned int f[6] = { 0, 1, 2, 3, 4, 5 }, i;
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	for (i = 0; i < 3; i++)
	{
		double y = p[i][1] * ca, z = p[i][1] * sa, xx = p[i][0];
		x[3 * i] = (qaws_scalar)(xx * cb - y * sb);
		x[3 * i + 1] = (qaws_scalar)(xx * sb + y * cb);
		x[3 * i + 2] = (qaws_scalar)z;
	}
	x[9] = 2; x[10] = 0; x[11] = 0;
	x[12] = 3; x[13] = 0; x[14] = 0;
	x[15] = 4; x[16] = 0; x[17] = 0;
	memset(&d, 0, sizeof d);
	d.positions = x;
	d.vertex_count = 6;
	d.face_vertices = f;
	d.face_count = 2;
	TEST_ASSERT(qaws_mesh_3d_create(&m, &d) == QAWS_STATUS_OK, "surface: mesh made");
	if (!m)
		return;
	qaws_mesh_3d_measure_elements(v, m, QAWS_MEASURE_SHAPE);
	TEST_ASSERT(qt_near(v[0], 1, 1e-5) && fabs(v[1]) < 1e-12, "surface: shape 1 tilted, 0 flat");
	qaws_mesh_3d_measure_elements(v, m, QAWS_MEASURE_ASPECT_RATIO);
	TEST_ASSERT(qt_near(v[0], 1, 1e-5) && v[1] > 1e30, "surface: aspect 1 tilted, infinite flat");
	qaws_mesh_3d_destroy(m);
	/* a perturbed 8 x 8 quad grid on a wavy surface against measure_shape */
	{
		unsigned int n = 8, nv = 81, j, k = 0, fq[256];
		qaws_scalar g[243], vals[64];
		qaws_shape_quality sq;
		double lo = HUGE_VAL, jl = HUGE_VAL;
		unsigned int s = 7;
		for (j = 0; j <= n; j++)
			for (i = 0; i <= n; i++)
			{
				s = s * 1664525u + 1013904223u;
				g[3 * (j * 9 + i)] = (qaws_scalar)(i + 0.25 * ((s >> 9) / 8388608.0 - 0.5));
				s = s * 1664525u + 1013904223u;
				g[3 * (j * 9 + i) + 1] = (qaws_scalar)(j + 0.25 * ((s >> 9) / 8388608.0 - 0.5));
				g[3 * (j * 9 + i) + 2] = (qaws_scalar)(0.3 * sin(0.7 * i) * cos(0.5 * j));
			}
		for (j = 0; j < n; j++)
			for (i = 0; i < n; i++)
			{
				unsigned int a = j * 9 + i;
				fq[k++] = a; fq[k++] = a + 1; fq[k++] = a + 10; fq[k++] = a + 9;
			}
		{
			unsigned int sizes[64];
			for (i = 0; i < 64; i++)
				sizes[i] = 4;
			d.positions = g;
			d.vertex_count = nv;
			d.face_vertices = fq;
			d.face_sizes = sizes;
			d.face_count = 64;
			TEST_ASSERT(qaws_mesh_3d_create(&m, &d) == QAWS_STATUS_OK, "surface: quad grid made");
		}
		qaws_mesh_3d_measure_shape(&sq, NULL, m);
		qaws_mesh_3d_measure_elements(vals, m, QAWS_MEASURE_MIN_ANGLE);
		for (i = 0; i < 64; i++)
			lo = vals[i] < lo ? vals[i] : lo;
		qaws_mesh_3d_measure_elements(vals, m, QAWS_MEASURE_SCALED_JACOBIAN);
		for (i = 0; i < 64; i++)
			jl = vals[i] < jl ? vals[i] : jl;
		printf("  quad grid: smallest angle %.4f (measure_shape %.4f), Jacobian %.4f (%.4f)\n", lo, (double)sq.min_angle, jl, (double)sq.min_jacobian);
		TEST_ASSERT(qt_near(lo, sq.min_angle, 1e-6) && qt_near(jl, sq.min_jacobian, 1e-6), "surface: as qaws_mesh_3d_measure_shape's smallest angle and Jacobian");
		qaws_mesh_3d_destroy(m);
	}
}

static qaws_volume_mesh* qt_volume(qaws_scalar const* x, unsigned int nv, unsigned int const* cv, qaws_cell_type const* ty, unsigned int nc)
{
	qaws_volume_mesh_desc d;
	qaws_volume_mesh* m = NULL;
	memset(&d, 0, sizeof d);
	d.positions = x;
	d.vertex_count = nv;
	d.cell_vertices = cv;
	d.cell_types = ty;
	d.cell_count = nc;
	qaws_volume_mesh_create(&m, &d);
	return m;
}

static void test_volume(void)
{
	double s3 = sqrt(3.0), s2 = sqrt(2.0);
	/* a regular tetrahedron of edge 1, a unit cube, a prism (equilateral base, height 1), a pyramid of equal edges: apart */
	qaws_scalar x[] = {
		0, 0, 0, 1, 0, 0, (qaws_scalar)0.5, (qaws_scalar)(s3 / 2), 0, (qaws_scalar)0.5, (qaws_scalar)(s3 / 6), (qaws_scalar)sqrt(2.0 / 3),
		3, 0, 0, 4, 0, 0, 4, 1, 0, 3, 1, 0, 3, 0, 1, 4, 0, 1, 4, 1, 1, 3, 1, 1,
		6, 0, 0, 7, 0, 0, (qaws_scalar)6.5, (qaws_scalar)(s3 / 2), 0, 6, 0, 1, 7, 0, 1, (qaws_scalar)6.5, (qaws_scalar)(s3 / 2), 1,
		9, 0, 0, 10, 0, 0, 10, 1, 0, 9, 1, 0, (qaws_scalar)9.5, (qaws_scalar)0.5, (qaws_scalar)(1 / s2) };
	unsigned int cv[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22 };
	qaws_cell_type ty[] = { QAWS_CELL_TETRAHEDRON, QAWS_CELL_HEXAHEDRON, QAWS_CELL_PRISM, QAWS_CELL_PYRAMID };
	qaws_scalar v[4];
	qaws_volume_mesh* m = qt_volume(x, 23, cv, ty, 4);
	TEST_ASSERT(m != NULL, "volume: cells made");
	if (!m)
		return;
	qaws_volume_mesh_measure_elements(v, m, QAWS_MEASURE_SCALED_JACOBIAN);
	TEST_ASSERT(qt_near(v[0], 1, 1e-5) && qt_near(v[1], 1, qt_tol()) && qt_near(v[2], 1, 1e-5) && qt_near(v[3], 1, 1e-5),
		"volume: scaled Jacobian 1 for the regular tetrahedron, the cube, the prism, the pyramid");
	qaws_volume_mesh_measure_elements(v, m, QAWS_MEASURE_SHAPE);
	TEST_ASSERT(qt_near(v[0], 1, 1e-5) && v[1] != v[1], "volume: shape 1 for the tetrahedron, NaN for the cube");
	qaws_volume_mesh_measure_elements(v, m, QAWS_MEASURE_ASPECT_RATIO);
	TEST_ASSERT(qt_near(v[0], 1, 1e-5), "volume: aspect 1 for the tetrahedron");
	qaws_volume_mesh_measure_elements(v, m, QAWS_MEASURE_MIN_ANGLE);
	TEST_ASSERT(qt_near(v[0], acos(1.0 / 3) * 180 / 3.14159265358979323846, 1e-5) && qt_near(v[1], 90, 1e-5) && qt_near(v[2], 60, 1e-5),
		"volume: dihedral 70.53 (tetrahedron), 90 (cube), 60 (prism)");
	qaws_volume_mesh_measure_elements(v, m, QAWS_MEASURE_MAX_ANGLE);
	TEST_ASSERT(qt_near(v[1], 90, 1e-5) && qt_near(v[2], 90, 1e-5), "volume: largest dihedral 90 (cube, prism)");
	qaws_volume_mesh_measure_elements(v, m, QAWS_MEASURE_SIZE);
	TEST_ASSERT(qt_near(v[0], 1 / (6 * s2), 1e-5) && qt_near(v[1], 1, qt_tol()) && qt_near(v[2], s3 / 4, 1e-5) && qt_near(v[3], 1 / (3 * s2), 1e-5),
		"volume: the volumes");
	qaws_volume_mesh_destroy(m);
	/* a sheared 2 x 2 x 2 block of hexahedra against measure_quality */
	{
		qaws_scalar g[81];
		unsigned int c[64], i, j, k, n = 0;
		qaws_cell_type t[8];
		qaws_volume_quality q;
		qaws_scalar vals[8];
		double lo = HUGE_VAL;
		for (k = 0; k < 3; k++)
			for (j = 0; j < 3; j++)
				for (i = 0; i < 3; i++)
				{
					unsigned int id = (k * 3 + j) * 3 + i;
					g[3 * id] = (qaws_scalar)(i + 0.3 * j + 0.1 * k * k);
					g[3 * id + 1] = (qaws_scalar)(j + 0.05 * i * k);
					g[3 * id + 2] = (qaws_scalar)(k + 0.1 * i);
				}
		for (k = 0; k < 2; k++)
			for (j = 0; j < 2; j++)
				for (i = 0; i < 2; i++, n++)
				{
					unsigned int a = (k * 3 + j) * 3 + i;
					c[8 * n] = a; c[8 * n + 1] = a + 1; c[8 * n + 2] = a + 4; c[8 * n + 3] = a + 3;
					c[8 * n + 4] = a + 9; c[8 * n + 5] = a + 10; c[8 * n + 6] = a + 13; c[8 * n + 7] = a + 12;
					t[n] = QAWS_CELL_HEXAHEDRON;
				}
		m = qt_volume(g, 27, c, t, 8);
		TEST_ASSERT(m != NULL, "volume: sheared block made");
		if (!m)
			return;
		qaws_volume_mesh_measure_quality(&q, m);
		qaws_volume_mesh_measure_elements(vals, m, QAWS_MEASURE_SCALED_JACOBIAN);
		for (i = 0; i < 8; i++)
			lo = vals[i] < lo ? vals[i] : lo;
		printf("  sheared block: smallest scaled Jacobian %.6f (measure_quality %.6f)\n", lo, (double)q.min_scaled_jacobian);
		TEST_ASSERT(qt_near(lo, q.min_scaled_jacobian, 1e-5), "volume: as qaws_volume_mesh_measure_quality's smallest scaled Jacobian");
		qaws_volume_mesh_destroy(m);
	}
}

static void test_summary(void)
{
	qaws_scalar v[9] = { 5, 1, 9, (qaws_scalar)NAN, 3, 7, 2, 8, 4 };
	unsigned int hist[4], worst[3];
	qaws_quality_summary s;
	qaws_quality_summary_desc d;
	memset(&d, 0, sizeof d);
	d.measure = QAWS_MEASURE_MIN_ANGLE;
	d.check = 1;
	d.threshold = (qaws_scalar)3.5;
	d.bins = 4;
	d.worst_count = 3;
	TEST_ASSERT(qaws_quality_summarize(&s, hist, worst, v, 9, &d) == QAWS_STATUS_OK, "summary: ran");
	TEST_ASSERT(s.count == 8 && s.skipped == 1 && s.min == 1 && s.max == 9 && s.mean == (qaws_scalar)4.875 && s.median == (qaws_scalar)4.5,
		"summary: count, NaN skipped, min, max, mean, median");
	TEST_ASSERT(s.failing == 3, "summary: three below the threshold (lower is worse)");
	TEST_ASSERT(s.worst_found == 3 && worst[0] == 1 && worst[1] == 6 && worst[2] == 4, "summary: the worst three, the smallest first");
	/* span 1..9 in 4 bins: [1,3) 1 2, [3,5) 3 4, [5,7) 5, [7,9] 7 8 9 */
	TEST_ASSERT(hist[0] == 2 && hist[1] == 2 && hist[2] == 1 && hist[3] == 3, "summary: the histogram over the values' span");
	d.measure = QAWS_MEASURE_ASPECT_RATIO;
	d.threshold = 7;
	d.range[0] = 0;
	d.range[1] = 4;
	qaws_quality_summarize(&s, hist, worst, v, 9, &d);
	TEST_ASSERT(s.failing == 2 && worst[0] == 2 && worst[1] == 7 && worst[2] == 5, "summary: higher is worse for aspect: two above 7, the largest first");
	TEST_ASSERT(hist[0] == 0 && hist[1] == 1 && hist[2] == 1 && hist[3] == 6 && s.range[1] == 4, "summary: a given span, values past it in the end bins");
	d.bins = 0;
	d.worst_count = 0;
	TEST_ASSERT(qaws_quality_summarize(&s, NULL, NULL, v, 0, &d) == QAWS_STATUS_OK && s.count == 0, "summary: no values");
}

static void test_refusals(void)
{
	qaws_scalar x[9] = { 0, 0, 0, 1, 0, 0, 0, 1, 0 }, v[1];
	unsigned int f[3] = { 0, 1, 2 };
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	qaws_quality_summary s;
	qaws_quality_summary_desc sd;
	memset(&d, 0, sizeof d);
	d.positions = x;
	d.vertex_count = 3;
	d.face_vertices = f;
	d.face_count = 1;
	qaws_mesh_3d_create(&m, &d);
	TEST_ASSERT(qaws_mesh_3d_measure_elements(v, m, (qaws_element_measure)99) == QAWS_STATUS_INVALID_ARGUMENT, "refusals: a measure out of range");
	TEST_ASSERT(qaws_mesh_3d_measure_elements(NULL, m, QAWS_MEASURE_SIZE) == QAWS_STATUS_INVALID_ARGUMENT, "refusals: no output");
	memset(&sd, 0, sizeof sd);
	sd.bins = 4;
	TEST_ASSERT(qaws_quality_summarize(&s, NULL, NULL, v, 1, &sd) == QAWS_STATUS_INVALID_ARGUMENT, "refusals: bins without a histogram");
	qaws_mesh_3d_destroy(m);
}

int test_160_quality_report_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 160: quality reports\n");
	test_2d();
	test_surface();
	test_volume();
	test_summary();
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
