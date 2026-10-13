/*
 * Test 107: curves and regions across the faces of a shell
 * (qaws_shell_region.h)
 *
 *   - a cube of six flat faces: twelve shared edges; a straightest ray over
 *     an edge against the unfolded line; geodesic circles over one edge and
 *     over two (near a vertex, not round it) against pi r^2 and 2 pi r
 *     (the faces unfold flat about the edges); a disk round a corner against
 *     the unfolded cone of angle 3 pi / 2; two circles over edges by
 *     inclusion and exclusion
 *   - a sphere as two hemispheres meeting at the equator: a circle over
 *     the equator against the cap 2 pi R^2 (1 - cos(r / R)) and its
 *     perimeter; a great circle over the equator twice, back at its start
 *   - five squares round a saddle vertex (angle 5 pi / 2): a disk through it
 *     against the plane disk plus the sector behind it, (pi / 4)(r - rho0)^2
 *   - six squares round a vertex (a monkey saddle, angle 3 pi): a disk
 *     against pi r^2 + (pi / 2)(r - rho0)^2; the same six as hypars (far
 *     corners pushed off the plane, the faces curved): the disk at two
 *     tolerances, warped either way alike, larger than the flat one
 */

#include "test_common.h"
#include "qaws_shell_region.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define SHT_PI 3.14159265358979323846
#define SHT_F32 (sizeof(qaws_scalar) == 4)

static qaws_surface* sht_face(double const* p00, double const* p10, double const* p01, double const* p11)
{
	qaws_surface_bilinear_desc d;
	qaws_surface* s = NULL;
	memset(&d, 0, sizeof(d));
	d.p00.x = (qaws_scalar)p00[0]; d.p00.y = (qaws_scalar)p00[1]; d.p00.z = (qaws_scalar)p00[2];
	d.p10.x = (qaws_scalar)p10[0]; d.p10.y = (qaws_scalar)p10[1]; d.p10.z = (qaws_scalar)p10[2];
	d.p01.x = (qaws_scalar)p01[0]; d.p01.y = (qaws_scalar)p01[1]; d.p01.z = (qaws_scalar)p01[2];
	d.p11.x = (qaws_scalar)p11[0]; d.p11.y = (qaws_scalar)p11[1]; d.p11.z = (qaws_scalar)p11[2];
	qaws_surface_create_bilinear(&d, &s);
	return s;
}

/* the cube [-1, 1]^3, faces oriented outward: +z, -z, +x, -x, +y, -y */
static void sht_cube(qaws_surface** f)
{
	static double const P[8][3] = { { -1, -1, -1 }, { 1, -1, -1 }, { -1, 1, -1 }, { 1, 1, -1 }, { -1, -1, 1 }, { 1, -1, 1 }, { -1, 1, 1 }, { 1, 1, 1 } };
	f[0] = sht_face(P[4], P[5], P[6], P[7]);
	f[1] = sht_face(P[0], P[2], P[1], P[3]);
	f[2] = sht_face(P[1], P[3], P[5], P[7]);
	f[3] = sht_face(P[0], P[4], P[2], P[6]);
	f[4] = sht_face(P[2], P[6], P[3], P[7]);
	f[5] = sht_face(P[0], P[1], P[4], P[5]);
}

static qaws_shell_point sht_pt(unsigned int face, double u, double v)
{
	qaws_shell_point p;
	p.face = face; p.uv.x = (qaws_scalar)u; p.uv.y = (qaws_scalar)v;
	return p;
}

static double sht_circle_area(qaws_shell_atlas const* a, qaws_shell_point c, double r, double* per, unsigned int* runs, qaws_status* st)
{
	qaws_shell_curve* k = NULL;
	qaws_shell_region* g = NULL;
	qaws_scalar A = 0, L = 0;
	*st = qaws_shell_curve_create_geodesic_circle(a, c, (qaws_scalar)r, 0, &k);
	if (runs) *runs = k ? qaws_shell_curve_get_run_count(k) : 0;
	if (*st == QAWS_STATUS_OK) *st = qaws_shell_region_create(a, (qaws_shell_curve const* const*)&k, 1, &g);
	if (*st == QAWS_STATUS_OK) *st = qaws_shell_region_get_area(g, &A);
	if (*st == QAWS_STATUS_OK && per) *st = qaws_shell_region_get_perimeter(g, &L);
	if (per) *per = (double)L;
	qaws_shell_region_destroy(g);
	qaws_shell_curve_destroy(k);
	return (double)A;
}

static void test_cube(void)
{
	qaws_surface* f[6];
	qaws_shell_atlas* a = NULL;
	qaws_status st;
	double tol = SHT_F32 ? 2e-4 : 1e-9;
	unsigned int i;
	char msg[300];
	sht_cube(f);
	st = qaws_shell_atlas_create((qaws_surface const* const*)f, 6, 0, &a);
	sprintf(msg, "cube: the atlas of six faces, %u shared edges", a ? qaws_shell_atlas_get_edge_count(a) : 0);
	TEST_ASSERT(st == QAWS_STATUS_OK && qaws_shell_atlas_get_edge_count(a) == 12, msg);
	/* a ray from (0.3, -0.2, 1) at 30 degrees from +x, run 1.5: over the
	   edge x = 1 and on down the face, unfolded a straight line */
	{
		qaws_shell_point q;
		qaws_vec2 d;
		qaws_vec3 vel;
		unsigned int cr = 0;
		qaws_surface_eval_result e;
		double x0 = 0.3, y0 = -0.2, th = SHT_PI / 6, L = 1.5, tx = (1 - x0) / cos(th), ex, ey, ez, err = 1e300;
		d.x = (qaws_scalar)cos(th); d.y = (qaws_scalar)sin(th);
		st = qaws_shell_exp(a, sht_pt(0, (x0 + 1) / 2, (y0 + 1) / 2), d, (qaws_scalar)L, &q, &vel, &cr);
		/* unfolded: the rest of the run goes down x = 1 */
		ex = 1; ey = y0 + L * sin(th); ez = 1 - (L - tx) * cos(th);
		if (st == QAWS_STATUS_OK)
		{
			qaws_surface_evaluate(f[q.face], q.uv.x, q.uv.y, QAWS_SURFACE_EVAL_POSITION, &e);
			err = fabs(e.position.x - ex) + fabs(e.position.y - ey) + fabs(e.position.z - ez);
		}
		sprintf(msg, "cube: a straightest ray over an edge (%u crossing) against the unfolded line within %.2g", cr, err);
		TEST_ASSERT(st == QAWS_STATUS_OK && cr == 1 && q.face == 2 && err <= tol, msg);
	}
	/* circles: over one edge, over two near a vertex (not round it) */
	{
		double per1, per2, a1, a2, r1 = 0.8, r2 = 0.5;
		unsigned int n1, n2;
		qaws_status s1, s2;
		a1 = sht_circle_area(a, sht_pt(0, 0.75, 0.5), r1, &per1, &n1, &s1);
		a2 = sht_circle_area(a, sht_pt(0, 0.8, 0.8), r2, &per2, &n2, &s2);
		sprintf(msg, "cube: geodesic circles over one edge (%u runs) and two edges by a vertex (%u runs): areas against pi r^2 %.2g, %.2g; perimeters against 2 pi r %.2g, %.2g",
			n1, n2, fabs(a1 - SHT_PI * r1 * r1), fabs(a2 - SHT_PI * r2 * r2), fabs(per1 - 2 * SHT_PI * r1), fabs(per2 - 2 * SHT_PI * r2));
		TEST_ASSERT(!s1 && !s2 && n1 == 2 && n2 == 4 && fabs(a1 - SHT_PI * r1 * r1) <= tol && fabs(a2 - SHT_PI * r2 * r2) <= tol &&
			fabs(per1 - 2 * SHT_PI * r1) <= tol && fabs(per2 - 2 * SHT_PI * r2) <= tol, msg);
	}
	/* two circles over edges: union + intersection = both */
	{
		qaws_shell_curve *k1 = NULL, *k2 = NULL;
		qaws_shell_region *g1 = NULL, *g2 = NULL, *gu = NULL, *gi = NULL;
		qaws_scalar A1 = 0, A2 = 0, AU = 0, AI = 0;
		qaws_status s1;
		s1 = qaws_shell_curve_create_geodesic_circle(a, sht_pt(0, 0.75, 0.5), (qaws_scalar)0.7, 0, &k1);
		if (!s1) s1 = qaws_shell_curve_create_geodesic_circle(a, sht_pt(2, 0.6, 0.8), (qaws_scalar)0.6, 0, &k2);
		if (!s1) s1 = qaws_shell_region_create(a, (qaws_shell_curve const* const*)&k1, 1, &g1);
		if (!s1) s1 = qaws_shell_region_create(a, (qaws_shell_curve const* const*)&k2, 1, &g2);
		if (!s1) s1 = qaws_shell_region_clip(g1, g2, QAWS_CLIP_UNION, &gu);
		if (!s1) s1 = qaws_shell_region_clip(g1, g2, QAWS_CLIP_INTERSECTION, &gi);
		if (!s1)
		{
			qaws_shell_region_get_area(g1, &A1); qaws_shell_region_get_area(g2, &A2);
			qaws_shell_region_get_area(gu, &AU); qaws_shell_region_get_area(gi, &AI);
		}
		sprintf(msg, "cube: two circles over edges: union + intersection - both %.2g (intersection %.4f)", (double)(AU + AI - A1 - A2), (double)AI);
		TEST_ASSERT(!s1 && fabs((double)(AU + AI - A1 - A2)) <= tol && AI > 0.05, msg);
		qaws_shell_region_destroy(g1); qaws_shell_region_destroy(g2); qaws_shell_region_destroy(gu); qaws_shell_region_destroy(gi);
		qaws_shell_curve_destroy(k1); qaws_shell_curve_destroy(k2);
	}
	/* a disk round a corner (a cone of angle 3 pi / 2): unfolded, the plane
	   disk less its part in the missing quarter behind the apex,
	   int 1/2 R(phi)^2 over [3 pi / 4, 5 pi / 4], R the distance from the
	   apex to the circle that way (Simpson) */
	{
		double t = 0.2, rho0 = t * sqrt(2.0), r = 0.6, back = 0, A0 = 0;
		unsigned int n = 20000, q, nf = 0;
		qaws_shell_region* g = NULL;
		qaws_scalar A = 0;
		qaws_status s1;
		for (q = 0; q <= n; q++)
		{
			double phi = 0.75 * SHT_PI + 0.5 * SHT_PI * q / n, R = rho0 * cos(phi) + sqrt(r * r - rho0 * rho0 * sin(phi) * sin(phi));
			back += (q == 0 || q == n ? 1 : (q % 2 ? 4 : 2)) * 0.5 * R * R;
		}
		back *= 0.5 * SHT_PI / n / 3;
		A0 = SHT_PI * r * r - back;
		s1 = qaws_shell_region_create_geodesic_disk(a, sht_pt(0, 1 - t / 2, 1 - t / 2), (qaws_scalar)r, 0, &g);
		if (s1 == QAWS_STATUS_OK)
		{
			qaws_shell_region_get_area(g, &A);
			for (q = 0; q < 6; q++) if (qaws_shell_region_get_face(g, q)) nf++;
		}
		sprintf(msg, "cube: a geodesic disk round a corner (on %u faces) against the unfolded cone's %.2g", nf, fabs((double)A - A0));
		TEST_ASSERT(s1 == QAWS_STATUS_OK && nf == 3 && fabs((double)A - A0) <= (SHT_F32 ? 2e-4 : 1e-9), msg);
		qaws_shell_region_destroy(g);
	}
	qaws_shell_atlas_destroy(a);
	for (i = 0; i < 6; i++) qaws_surface_destroy(f[i]);
}

static qaws_surface* sht_hemi(double a0, double a1, qaws_curve** prof)
{
	qaws_arc_segment g;
	qaws_arc_desc ad;
	qaws_surface_revolution_desc rd;
	qaws_surface* s = NULL;
	memset(&g, 0, sizeof(g));
	g.radius = 1; g.angle_start = (qaws_scalar)a0; g.angle_end = (qaws_scalar)a1;
	memset(&ad, 0, sizeof(ad));
	ad.dimension = QAWS_DIMENSION_2D; ad.segments = &g; ad.segment_count = 1;
	qaws_curve_create_arc(&ad, prof);
	memset(&rd, 0, sizeof(rd));
	rd.profile = *prof; rd.axis_direction.z = 1;
	qaws_surface_create_revolution(&rd, &s);
	return s;
}

static void test_hemispheres(void)
{
	qaws_curve *p1 = NULL, *p2 = NULL;
	qaws_surface* h[2];
	qaws_shell_atlas* a = NULL;
	qaws_status st, s1;
	double tol = SHT_F32 ? 2e-4 : 1e-9, r = 0.6, A, per;
	unsigned int n = 0;
	char msg[300];
	h[0] = sht_hemi(-SHT_PI / 2, 0, &p1);
	h[1] = sht_hemi(0, SHT_PI / 2, &p2);
	st = qaws_shell_atlas_create((qaws_surface const* const*)h, 2, 0, &a);
	A = sht_circle_area(a, sht_pt(1, 0.3, 0.1), r, &per, &n, &s1);
	sprintf(msg, "hemispheres: one edge (the equator); a circle over it (%u runs) against the cap %.2g and 2 pi sin r %.2g",
		n, fabs(A - 2 * SHT_PI * (1 - cos(r))), fabs(per - 2 * SHT_PI * sin(r)));
	TEST_ASSERT(st == QAWS_STATUS_OK && qaws_shell_atlas_get_edge_count(a) == 1 && !s1 && n == 2 &&
		fabs(A - 2 * SHT_PI * (1 - cos(r))) <= tol && fabs(per - 2 * SHT_PI * sin(r)) <= tol, msg);
	/* a great circle from the equator heading north-east: back after 2 pi */
	{
		qaws_shell_point q;
		qaws_vec2 d;
		qaws_vec3 vel;
		unsigned int cr = 0;
		double du, dv;
		d.x = (qaws_scalar)0.1; d.y = (qaws_scalar)0.2;
		st = qaws_shell_exp(a, sht_pt(1, 0.1, 0.02), d, (qaws_scalar)(2 * SHT_PI), &q, &vel, &cr);
		du = fabs(q.uv.x - 0.1); du = fmin(du, 1 - du);
		dv = fabs(q.uv.y - 0.02);
		sprintf(msg, "hemispheres: a great circle over the equator (%u crossings) back at its start within %.2g", cr, du + dv);
		TEST_ASSERT(st == QAWS_STATUS_OK && cr == 2 && q.face == 1 && du + dv <= tol, msg);
	}
	qaws_shell_atlas_destroy(a);
	qaws_surface_destroy(h[0]); qaws_surface_destroy(h[1]);
	qaws_curve_destroy(p1); qaws_curve_destroy(p2);
}

/* five unit squares round the origin (three of the floor, two walls from
   it): a saddle vertex of angle 5 pi / 2. A disk round it through a point
   at rho0 from it: the plane disk, and behind the vertex the sector of
   the angle left over (pi / 2) and radius r - rho0, reached through it */
static void test_saddle(void)
{
	static double const S[5][4][3] = {
		{ { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 1, 1, 0 } },
		{ { -1, 0, 0 }, { 0, 0, 0 }, { -1, 1, 0 }, { 0, 1, 0 } },
		{ { -1, -1, 0 }, { 0, -1, 0 }, { -1, 0, 0 }, { 0, 0, 0 } },
		{ { 0, -1, -1 }, { 0, 0, -1 }, { 0, -1, 0 }, { 0, 0, 0 } },
		{ { 0, 0, -1 }, { 1, 0, -1 }, { 0, 0, 0 }, { 1, 0, 0 } } };
	qaws_surface* f[5];
	qaws_shell_atlas* a = NULL;
	qaws_shell_region* g = NULL;
	qaws_scalar A = 0;
	double rho0 = 0.2 * sqrt(2.0), r = 0.5, A0 = SHT_PI * r * r + 0.25 * SHT_PI * (r - rho0) * (r - rho0);
	unsigned int i, nf = 0;
	qaws_status st, s1 = QAWS_STATUS_INVALID_ARGUMENT;
	char msg[300];
	for (i = 0; i < 5; i++) f[i] = sht_face(S[i][0], S[i][1], S[i][2], S[i][3]);
	st = qaws_shell_atlas_create((qaws_surface const* const*)f, 5, 0, &a);
	if (st == QAWS_STATUS_OK)
		s1 = qaws_shell_region_create_geodesic_disk(a, sht_pt(0, 0.2, 0.2), (qaws_scalar)r, 0, &g);
	if (s1 == QAWS_STATUS_OK)
	{
		qaws_shell_region_get_area(g, &A);
		for (i = 0; i < 5; i++) if (qaws_shell_region_get_face(g, i)) nf++;
	}
	sprintf(msg, "saddle (five squares round a vertex): a disk through it (on %u faces) against the plane disk and the sector behind %.2g",
		nf, fabs((double)A - A0));
	TEST_ASSERT(st == QAWS_STATUS_OK && s1 == QAWS_STATUS_OK && nf == 5 && fabs((double)A - A0) <= (SHT_F32 ? 2e-4 : 1e-9), msg);
	qaws_shell_region_destroy(g);
	qaws_shell_atlas_destroy(a);
	for (i = 0; i < 5; i++) qaws_surface_destroy(f[i]);
}

/* six faces round the origin, face i spanned by axis directions d[i] and
   d[i + 1] (the monkey saddle), the far corner pushed by w along the normal */
static void sht_monkey(qaws_surface** f, double w)
{
	static double const D[6][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { -1, 0, 0 }, { 0, -1, 0 }, { 0, 0, -1 } };
	static double const O[3] = { 0, 0, 0 };
	unsigned int i, k;
	for (i = 0; i < 6; i++)
	{
		double const* a = D[i];
		double const* b = D[(i + 1) % 6];
		double n[3] = { a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0] }, p11[3];
		for (k = 0; k < 3; k++) p11[k] = a[k] + b[k] + w * n[k];
		f[i] = sht_face(O, a, b, p11);
	}
}

static double sht_monkey_disk(double w, double tolerance, unsigned int* faces)
{
	qaws_surface* f[6];
	qaws_shell_atlas* a = NULL;
	qaws_shell_region* g = NULL;
	qaws_scalar A = -1;
	unsigned int i;
	sht_monkey(f, w);
	*faces = 0;
	if (qaws_shell_atlas_create((qaws_surface const* const*)f, 6, 0, &a) == QAWS_STATUS_OK
		&& qaws_shell_region_create_geodesic_disk(a, sht_pt(0, 0.3, 0.3), (qaws_scalar)0.6, (qaws_scalar)tolerance, &g) == QAWS_STATUS_OK)
	{
		qaws_shell_region_get_area(g, &A);
		for (i = 0; i < 6; i++) if (qaws_shell_region_get_face(g, i)) (*faces)++;
	}
	qaws_shell_region_destroy(g);
	qaws_shell_atlas_destroy(a);
	for (i = 0; i < 6; i++) qaws_surface_destroy(f[i]);
	return (double)A;
}

static void test_monkey_saddle(void)
{
	double rho0 = 0.3 * sqrt(2.0), r = 0.6, A0 = SHT_PI * r * r + 0.5 * SHT_PI * (r - rho0) * (r - rho0);
	double flat, fine, coarse, other;
	unsigned int n0, n1, n2, n3;
	char msg[300];
	flat = sht_monkey_disk(0, 0, &n0);
	sprintf(msg, "monkey saddle (six squares, angle 3 pi): a disk through the vertex (on %u faces) against pi r^2 + (pi/2)(r - rho)^2 %.2g",
		n0, fabs(flat - A0));
	TEST_ASSERT(n0 == 6 && fabs(flat - A0) <= (SHT_F32 ? 2e-4 : 1e-9), msg);
	fine = sht_monkey_disk(0.4, 0, &n1);
	coarse = sht_monkey_disk(0.4, SHT_F32 ? 1e-3 : 1e-6, &n2);
	other = sht_monkey_disk(-0.4, 0, &n3);
	sprintf(msg, "monkey saddle of hypars: the disk %.9f (flat %.9f) at two tolerances %.2g, warped the other way %.2g",
		fine, flat, fabs(fine - coarse), fabs(fine - other));
	TEST_ASSERT(n1 == 6 && n2 == 6 && n3 == 6 && fine > flat + 1e-3 && fabs(fine - coarse) <= (SHT_F32 ? 2e-3 : 1e-6)
		&& fabs(fine - other) <= (SHT_F32 ? 2e-4 : 1e-9), msg);
}

int test_107_shell_region_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 107: curves and regions across the faces of a shell\n");
	test_cube();
	test_hemispheres();
	test_saddle();
	test_monkey_saddle();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
