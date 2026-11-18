/*
 * Test 85: seams and shatter of curved surfaces (qaws_shatter.h)
 *
 *   - two square sheets crossing at right angles: one open seam along their
 *     line of crossing, its length, every point on both sheets;
 *   - a saddle (a bilinear patch) and a plane above its middle: two open
 *     seams (hyperbola branches), every point on both surfaces;
 *   - A's own seams: two crossing sheets given as one operand (one A-A
 *     seam), a cube of six faces (none: they only share borders);
 *   - shatter: a square sheet by a cube (six bilinear faces, a solid): two
 *     pieces, the one inside the cube a quarter of the sheet; the saddle by
 *     the plane: three pieces, not told inside or outside;
 *   - refusals.
 */

#include "test_common.h"
#include <math.h>
#include <string.h>

#if QAWS_SCALAR_IS_FLOAT
#define SHT_TOL 2e-4
#else
#define SHT_TOL 1e-8
#endif

static qaws_vec3 sht_v(double x, double y, double z)
{
	qaws_vec3 v;
	v.x = (qaws_scalar)x;
	v.y = (qaws_scalar)y;
	v.z = (qaws_scalar)z;
	return v;
}

static qaws_surface* sht_bilinear(qaws_vec3 p00, qaws_vec3 p10, qaws_vec3 p01, qaws_vec3 p11)
{
	qaws_surface_bilinear_desc d;
	qaws_surface* s = NULL;
	d.p00 = p00;
	d.p10 = p10;
	d.p01 = p01;
	d.p11 = p11;
	qaws_surface_create_bilinear(&s, &d);
	return s;
}

/* the six faces of the box [lo, hi], normals outward (S_u x S_v) */
static void sht_box(double const* lo, double const* hi, qaws_surface** f)
{
	double x0 = lo[0], y0 = lo[1], z0 = lo[2], x1 = hi[0], y1 = hi[1], z1 = hi[2];
	f[0] = sht_bilinear(sht_v(x0, y0, z0), sht_v(x0, y1, z0), sht_v(x1, y0, z0), sht_v(x1, y1, z0));   /* z = z0, -z */
	f[1] = sht_bilinear(sht_v(x0, y0, z1), sht_v(x1, y0, z1), sht_v(x0, y1, z1), sht_v(x1, y1, z1));   /* z = z1, +z */
	f[2] = sht_bilinear(sht_v(x0, y0, z0), sht_v(x1, y0, z0), sht_v(x0, y0, z1), sht_v(x1, y0, z1));   /* y = y0, -y */
	f[3] = sht_bilinear(sht_v(x0, y1, z0), sht_v(x0, y1, z1), sht_v(x1, y1, z0), sht_v(x1, y1, z1));   /* y = y1, +y */
	f[4] = sht_bilinear(sht_v(x0, y0, z0), sht_v(x0, y0, z1), sht_v(x0, y1, z0), sht_v(x0, y1, z1));   /* x = x0, -x */
	f[5] = sht_bilinear(sht_v(x1, y0, z0), sht_v(x1, y1, z0), sht_v(x1, y0, z1), sht_v(x1, y1, z1));   /* x = x1, +x */
}

/* the largest gap between the two surfaces at a seam's points, and its length */
static void sht_check_seam(qaws_surface_seams const* S, unsigned int i, qaws_surface const* const* a, qaws_surface const* const* b, double* gap,
	double* len, qaws_seam_info* info)
{
	qaws_seam_point* p;
	unsigned int j;
	qaws_seam_info in;
	qaws_surface_seams_get(&in, NULL, S, i);
	p = (qaws_seam_point*)malloc(sizeof(qaws_seam_point) * (in.point_count + 1));
	qaws_surface_seams_get(NULL, p, S, i);
	*gap = 0;
	*len = 0;
	for (j = 0; j < in.point_count; j++)
	{
		qaws_surface_eval_result r1, r2;
		qaws_surface const* s1 = in.operand1 ? b[in.surface1] : a[in.surface1];
		qaws_surface const* s2 = in.operand2 ? b[in.surface2] : a[in.surface2];
		double dx, dy, dz;
		qaws_surface_evaluate(&r1, s1, p[j].u1, p[j].v1, QAWS_EVAL_FLAG_POSITION);
		qaws_surface_evaluate(&r2, s2, p[j].u2, p[j].v2, QAWS_EVAL_FLAG_POSITION);
		dx = (double)r1.position.x - (double)r2.position.x;
		dy = (double)r1.position.y - (double)r2.position.y;
		dz = (double)r1.position.z - (double)r2.position.z;
		*gap = fmax(*gap, sqrt(dx * dx + dy * dy + dz * dz));
		if (j + 1 < in.point_count + (in.closed ? 1u : 0u))
		{
			qaws_seam_point const* q = &p[(j + 1) % in.point_count];
			dx = (double)q->position.x - (double)p[j].position.x;
			dy = (double)q->position.y - (double)p[j].position.y;
			dz = (double)q->position.z - (double)p[j].position.z;
			*len += sqrt(dx * dx + dy * dy + dz * dz);
		}
	}
	free(p);
	*info = in;
}

static void test_seams(void)
{
	qaws_surface *sh = sht_bilinear(sht_v(-1, -1, 0), sht_v(1, -1, 0), sht_v(-1, 1, 0), sht_v(1, 1, 0));
	qaws_surface *up = sht_bilinear(sht_v(0.13, -1.3, -1), sht_v(0.13, 1.3, -1), sht_v(0.13, -1.3, 1), sht_v(0.13, 1.3, 1));
	qaws_surface *saddle = sht_bilinear(sht_v(-1, -1, 0.5), sht_v(1, -1, -0.5), sht_v(-1, 1, -0.5), sht_v(1, 1, 0.5));
	qaws_surface *plane = sht_bilinear(sht_v(-2, -2, 0.2), sht_v(2, -2, 0.2), sht_v(-2, 2, 0.2), sht_v(2, 2, 0.2));
	qaws_surface* cube[6];
	double const lo[3] = { -0.5, -0.5, -0.5 }, hi[3] = { 0.5, 0.5, 0.5 };
	qaws_surface_seam_desc d;
	qaws_surface_seams* S = NULL;
	qaws_seam_info in;
	double gap, len;
	unsigned int i;
	char msg[256];
	sht_box(lo, hi, cube);
	/* two sheets crossing along x = 0.13, y in [-1, 1] */
	memset(&d, 0, sizeof d);
	d.a = (qaws_surface const* const*)&sh;
	d.a_count = 1;
	d.b = (qaws_surface const* const*)&up;
	d.b_count = 1;
	d.kinds = QAWS_SEAM_A_B;
	TEST_ASSERT_STATUS(qaws_surface_seams_create(&S, &d));
	TEST_ASSERT(qaws_surface_seams_get_count(S) == 1, "two sheets: one seam");
	if (qaws_surface_seams_get_count(S) == 1)
	{
		sht_check_seam(S, 0, d.a, d.b, &gap, &len, &in);
		sprintf(msg, "two sheets: open %d, kind %u, length %.12f (2), gap %.3g", !in.closed, in.kind, len, gap);
		printf("  %s\n", msg);
		TEST_ASSERT(!in.closed && in.kind == QAWS_SEAM_A_B && in.operand1 == 0 && in.operand2 == 1 && fabs(len - 2) < 0.02 && gap < SHT_TOL, msg);
	}
	qaws_surface_seams_destroy(S);
	/* a saddle and a plane: two hyperbola branches */
	d.a = (qaws_surface const* const*)&saddle;
	d.b = (qaws_surface const* const*)&plane;
	S = NULL;
	TEST_ASSERT_STATUS(qaws_surface_seams_create(&S, &d));
	sprintf(msg, "saddle and plane: %u seams (2)", qaws_surface_seams_get_count(S));
	printf("  %s\n", msg);
	TEST_ASSERT(qaws_surface_seams_get_count(S) == 2, msg);
	for (i = 0; i < qaws_surface_seams_get_count(S); i++)
	{
		sht_check_seam(S, i, d.a, d.b, &gap, &len, &in);
		sprintf(msg, "saddle and plane, seam %u: %u points, open %d, length %.6f, gap %.3g", i, in.point_count, !in.closed, len, gap);
		printf("  %s\n", msg);
		TEST_ASSERT(!in.closed && in.point_count > 10 && gap < SHT_TOL, msg);
	}
	qaws_surface_seams_destroy(S);
	/* A's own: two sheets as one operand cross; a cube's faces only share borders */
	{
		qaws_surface* two[2];
		two[0] = sh;
		two[1] = up;
		memset(&d, 0, sizeof d);
		d.a = (qaws_surface const* const*)two;
		d.a_count = 2;
		d.kinds = QAWS_SEAM_A_A;
		S = NULL;
		TEST_ASSERT_STATUS(qaws_surface_seams_create(&S, &d));
		TEST_ASSERT(qaws_surface_seams_get_count(S) == 1, "two sheets as one operand: one A-A seam");
		if (qaws_surface_seams_get_count(S) == 1)
		{
			sht_check_seam(S, 0, d.a, d.a, &gap, &len, &in);
			TEST_ASSERT(in.kind == QAWS_SEAM_A_A && fabs(len - 2) < 0.02 && gap < SHT_TOL, "two sheets as one operand: the A-A seam on both");
		}
		qaws_surface_seams_destroy(S);
		d.a = (qaws_surface const* const*)cube;
		d.a_count = 6;
		S = NULL;
		TEST_ASSERT_STATUS(qaws_surface_seams_create(&S, &d));
		sprintf(msg, "a cube's faces: %u A-A seams (0)", qaws_surface_seams_get_count(S));
		printf("  %s\n", msg);
		TEST_ASSERT(qaws_surface_seams_get_count(S) == 0, msg);
		qaws_surface_seams_destroy(S);
	}
	/* refusals */
	S = (qaws_surface_seams*)1;
	memset(&d, 0, sizeof d);
	TEST_ASSERT(qaws_surface_seams_create(&S, &d) == QAWS_STATUS_INVALID_ARGUMENT && S == NULL, "seams: no surface refused");
	d.a = (qaws_surface const* const*)&sh;
	d.a_count = 1;
	TEST_ASSERT(qaws_surface_seams_create(&S, &d) == QAWS_STATUS_INVALID_ARGUMENT, "seams: no kind refused");
	qaws_surface_destroy(sh);
	qaws_surface_destroy(up);
	qaws_surface_destroy(saddle);
	qaws_surface_destroy(plane);
	for (i = 0; i < 6; i++)
		qaws_surface_destroy(cube[i]);
}

/* the share of a trimmed piece's (u, v) domain inside its trim, sampled on a grid */
static double sht_trim_share(qaws_surface const* piece)
{
	unsigned int i, j, in_count = 0, n = 100;
	for (j = 0; j < n; j++)
		for (i = 0; i < n; i++)
		{
			int in = 0;
			qaws_surface_trim_contains(&in, piece, (qaws_scalar)((i + 0.5) / n), (qaws_scalar)((j + 0.5) / n));
			in_count += in != 0;
		}
	return (double)in_count / (n * n);
}

static void test_shatter(void)
{
	qaws_surface* sh = sht_bilinear(sht_v(-1, -1, 0.05), sht_v(1, -1, 0.05), sht_v(-1, 1, 0.05), sht_v(1, 1, 0.05));
	qaws_surface* saddle = sht_bilinear(sht_v(-1, -1, 0.5), sht_v(1, -1, -0.5), sht_v(-1, 1, -0.5), sht_v(1, 1, 0.5));
	qaws_surface* plane = sht_bilinear(sht_v(-2, -2, 0.2), sht_v(2, -2, 0.2), sht_v(-2, 2, 0.2), sht_v(2, 2, 0.2));
	qaws_surface* cube[6];
	double const lo[3] = { -0.5, -0.5, -0.5 }, hi[3] = { 0.5, 0.5, 0.5 };
	qaws_surface_seam_desc d;
	qaws_surface_shatter* H = NULL;
	unsigned int i, n;
	double share_in = 0, share_out = 0, total = 0;
	int bad = 0;
	char msg[256];
	sht_box(lo, hi, cube);
	memset(&d, 0, sizeof d);
	d.a = (qaws_surface const* const*)&sh;
	d.a_count = 1;
	d.b = (qaws_surface const* const*)cube;
	d.b_count = 6;
	TEST_ASSERT_STATUS(qaws_surface_shatter_create(&H, &d, 1));
	n = qaws_surface_shatter_get_count(H);
	for (i = 0; i < n; i++)
	{
		qaws_surface const* p = NULL;
		unsigned int from = 99;
		int inside = -2;
		double s;
		qaws_surface_shatter_get(&p, &from, &inside, H, i);
		s = sht_trim_share(p);
		bad += from != 0 || inside < 0;
		if (inside == 1)
			share_in += s;
		else
			share_out += s;
		total += s;
	}
	sprintf(msg, "sheet by a cube: %u pieces (2), inside %.4f (0.25), outside %.4f (0.75)", n, share_in, share_out);
	printf("  %s\n", msg);
	TEST_ASSERT(n == 2 && !bad && fabs(share_in - 0.25) < 0.01 && fabs(share_out - 0.75) < 0.01, msg);
	qaws_surface_shatter_destroy(H);
	/* the saddle by the plane: three pieces, not told */
	d.a = (qaws_surface const* const*)&saddle;
	d.b = (qaws_surface const* const*)&plane;
	d.b_count = 1;
	H = NULL;
	TEST_ASSERT_STATUS(qaws_surface_shatter_create(&H, &d, 0));
	n = qaws_surface_shatter_get_count(H);
	total = 0;
	bad = 0;
	for (i = 0; i < n; i++)
	{
		qaws_surface const* p = NULL;
		int inside = 0;
		qaws_surface_shatter_get(&p, NULL, &inside, H, i);
		total += sht_trim_share(p);
		bad += inside != -1;
	}
	sprintf(msg, "saddle by a plane: %u pieces (3), their shares summing to %.4f (1)", n, total);
	printf("  %s\n", msg);
	TEST_ASSERT(n == 3 && !bad && fabs(total - 1) < 0.01, msg);
	qaws_surface_shatter_destroy(H);
	/* refusals: B empty */
	d.b_count = 0;
	H = (qaws_surface_shatter*)1;
	TEST_ASSERT(qaws_surface_shatter_create(&H, &d, 0) == QAWS_STATUS_INVALID_ARGUMENT && H == NULL, "shatter: no B refused");
	qaws_surface_destroy(sh);
	qaws_surface_destroy(saddle);
	qaws_surface_destroy(plane);
	for (i = 0; i < 6; i++)
		qaws_surface_destroy(cube[i]);
}

int test_85_shatter_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 85: seams and shatter of curved surfaces\n");
	test_seams();
	test_shatter();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
