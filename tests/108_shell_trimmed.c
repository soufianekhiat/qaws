/*
 * Test 108: shells of trimmed faces (qaws_shell_atlas_create_trimmed)
 *
 *   - a plane split by a circle into two trimmed faces of different
 *     surfaces (the square with a hole, the disk; the circle parametrised
 *     differently in each): one curved edge; a ray over the disk and back
 *     out, a straight line; a disk across the edge, pi r^2 and 2 pi r, each
 *     face's part against the circles' lens
 *   - a plate with a round hole and a tube hanging from it (the hole's trim
 *     against the tube's untrimmed top side): straightest rays into the
 *     tube against the helices of the unrolled cylinder
 *   - a cube of planes trimmed to its faces: twelve edges; circles over one
 *     edge and over two by a vertex against pi r^2 and 2 pi r; a disk round
 *     a corner against the untrimmed cube's
 *   - T-junctions: a side meeting two faces side by side, flat and as walls:
 *     split where they meet; disks round the T point against pi r^2, 2 pi r
 *   - a cylinder trimmed by its rims (loops wrapping round) and a hole on the
 *     seam, filled by a disk of the same surface: disks across the hole and
 *     the seam against pi r^2 and 2 pi r (the cylinder unrolled), a ray over
 *     both against the unrolled line
 *   - a tube capped by one disk or by two half disks (T points on the rim of
 *     a closed face): disks over the rim alike either way and round the tube
 *   - a torus with a hole where its seams cross, filled by a disk of the same
 *     surface: disks over the hole and a ray against the whole torus
 *   - shortest paths on a plate with a round hole against the taut string
 *   - shortest paths round two and three holes in a row against the closed
 *     form; from a neighbouring face round a hole against the unfolded one
 *   - holes far smaller than a trace step stop the rays through them
 *   - an open loop refused; a centre off its trimmed face refused
 */

#include "test_common.h"
#include "qaws_shell_region.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define STR_PI 3.14159265358979323846
#define STR_F32 (sizeof(qaws_scalar) == 4)

static qaws_surface* str_plane(double const* p00, double const* p10, double const* p01, double const* p11)
{
	qaws_surface_bilinear_desc d;
	qaws_surface* s = NULL;
	memset(&d, 0, sizeof(d));
	d.p00.x = (qaws_scalar)p00[0]; d.p00.y = (qaws_scalar)p00[1]; d.p00.z = (qaws_scalar)p00[2];
	d.p10.x = (qaws_scalar)p10[0]; d.p10.y = (qaws_scalar)p10[1]; d.p10.z = (qaws_scalar)p10[2];
	d.p01.x = (qaws_scalar)p01[0]; d.p01.y = (qaws_scalar)p01[1]; d.p01.z = (qaws_scalar)p01[2];
	d.p11.x = (qaws_scalar)p11[0]; d.p11.y = (qaws_scalar)p11[1]; d.p11.z = (qaws_scalar)p11[2];
	qaws_surface_create_bilinear(&s, &d);
	return s;
}

static qaws_curve* str_line(double x0, double y0, double x1, double y1)
{
	qaws_bezier_desc bd;
	qaws_scalar cp[4];
	qaws_curve* c = NULL;
	cp[0] = (qaws_scalar)x0; cp[1] = (qaws_scalar)y0; cp[2] = (qaws_scalar)x1; cp[3] = (qaws_scalar)y1;
	memset(&bd, 0, sizeof(bd));
	bd.dimension = QAWS_DIMENSION_2D; bd.degree = 1; bd.control_points = cp; bd.control_point_count = 2;
	qaws_curve_create_bezier(&c, &bd);
	return c;
}

static qaws_curve* str_arc(double cx, double cy, double r, double a0, double a1)
{
	qaws_arc_segment g;
	qaws_arc_desc ad;
	qaws_curve* c = NULL;
	memset(&g, 0, sizeof(g));
	g.center[0] = (qaws_scalar)cx; g.center[1] = (qaws_scalar)cy; g.radius = (qaws_scalar)r;
	g.angle_start = (qaws_scalar)a0; g.angle_end = (qaws_scalar)a1;
	memset(&ad, 0, sizeof(ad));
	ad.dimension = QAWS_DIMENSION_2D; ad.segments = &g; ad.segment_count = 1;
	qaws_curve_create_arc(&c, &ad);
	return c;
}

/* the square [lo, hi]^2 of (u, v), counter-clockwise */
static void str_square(qaws_curve** c, double lo, double hi)
{
	c[0] = str_line(lo, lo, hi, lo);
	c[1] = str_line(hi, lo, hi, hi);
	c[2] = str_line(hi, hi, lo, hi);
	c[3] = str_line(lo, hi, lo, lo);
}

static qaws_shell_point str_pt(unsigned int face, double u, double v)
{
	qaws_shell_point p;
	p.face = face; p.uv.x = (qaws_scalar)u; p.uv.y = (qaws_scalar)v;
	return p;
}

static qaws_vec3 str_at(qaws_shell_atlas const* a, qaws_shell_point p)
{
	qaws_surface_eval_result r;
	memset(&r, 0, sizeof(r));
	qaws_surface_evaluate(&r, qaws_shell_atlas_get_face(a, p.face), p.uv.x, p.uv.y, QAWS_SURFACE_EVAL_POSITION);
	return r.position;
}

/* the plane z = 0 as (x, y) = (-2 + 4u, -2 + 4v) */
static double const STR_P00[3] = { -2, -2, 0 }, STR_P10[3] = { 2, -2, 0 }, STR_P01[3] = { -2, 2, 0 }, STR_P11[3] = { 2, 2, 0 };

static void test_split_plane(void)
{
	/* the disk's own surface: [-1, 1]^2 with u along y and v along -x */
	static double const Q00[3] = { 1, -1, 0 }, Q10[3] = { 1, 1, 0 }, Q01[3] = { -1, -1, 0 }, Q11[3] = { -1, 1, 0 };
	qaws_surface *sa = str_plane(STR_P00, STR_P10, STR_P01, STR_P11), *sb = str_plane(Q00, Q10, Q01, Q11);
	qaws_curve *oc[4], *hole, *disk;
	qaws_curve const *ol[4], *hl[1], *dl[1];
	qaws_path_2d pa[2], pb[1];
	qaws_shell_face_desc fd[2];
	qaws_shell_atlas* a = NULL;
	qaws_status st;
	double tol = STR_F32 ? 1e-4 : 1e-9;
	unsigned int i;
	char msg[300];
	str_square(oc, 0.125, 0.875);
	hole = str_arc(0.5, 0.5, 0.2, 0, 2 * STR_PI);
	disk = str_arc(0.5, 0.5, 0.4, 1.0, 1.0 + 2 * STR_PI);
	for (i = 0; i < 4; i++) ol[i] = oc[i];
	hl[0] = hole; dl[0] = disk;
	pa[0].curves = ol; pa[0].curve_count = 4; pa[0].closed = 1;
	pa[1].curves = hl; pa[1].curve_count = 1; pa[1].closed = 1;
	pb[0].curves = dl; pb[0].curve_count = 1; pb[0].closed = 1;
	fd[0].surface = sa; fd[0].trims = pa; fd[0].trim_count = 2;
	fd[1].surface = sb; fd[1].trims = pb; fd[1].trim_count = 1;
	st = qaws_shell_atlas_create_trimmed(&a, fd, 2, 0);
	sprintf(msg, "split plane: the square with a hole and the disk, %u shared edges (the circle, split where each face's circle starts)", a ? qaws_shell_atlas_get_edge_count(a) : 0);
	TEST_ASSERT(st == QAWS_STATUS_OK && qaws_shell_atlas_get_edge_count(a) == 2, msg);
	if (st == QAWS_STATUS_OK)
	{
		/* a ray from (-1.2, 0.1) over the disk and back out: a straight line */
		qaws_shell_point q;
		qaws_vec2 d;
		qaws_vec3 e;
		unsigned int cr = 0;
		double dx = 1, dy = 0.3, n = hypot(dx, dy), L = 2, off;
		d.x = (qaws_scalar)(dx / n); d.y = (qaws_scalar)(dy / n);
		st = qaws_shell_exp(&q, NULL, &cr, a, str_pt(0, (-1.2 + 2) / 4, (0.1 + 2) / 4), d, (qaws_scalar)L);
		e = str_at(a, q);
		off = hypot((double)e.x - (-1.2 + L * dx / n), (double)e.y - (0.1 + L * dy / n)) + fabs((double)e.z);
		sprintf(msg, "split plane: a ray over the disk and out (%u crossings) against the straight line %.2g", cr, off);
		TEST_ASSERT(st == QAWS_STATUS_OK && cr == 2 && q.face == 0 && off <= (STR_F32 ? 2e-6 : 1e-12), msg);
	}
	if (st == QAWS_STATUS_OK)
	{
		/* a disk from (0.9, 0.2), r 0.5, over the circle */
		qaws_shell_region* g = NULL;
		qaws_scalar A = 0, P = 0, a0 = 0, a1 = 0;
		double dd = hypot(0.9, 0.2), r1 = 0.5, r2 = 0.8, lens;
		lens = r1 * r1 * acos((dd * dd + r1 * r1 - r2 * r2) / (2 * dd * r1)) + r2 * r2 * acos((dd * dd + r2 * r2 - r1 * r1) / (2 * dd * r2))
			- 0.5 * sqrt((-dd + r1 + r2) * (dd + r1 - r2) * (dd - r1 + r2) * (dd + r1 + r2));
		st = qaws_shell_region_create_geodesic_disk(&g, a, str_pt(0, (0.9 + 2) / 4, (0.2 + 2) / 4), (qaws_scalar)r1, 0);
		if (st == QAWS_STATUS_OK)
		{
			qaws_shell_region_get_area(&A, g);
			qaws_shell_region_get_perimeter(&P, g);
			if (qaws_shell_region_get_face(g, 0)) qaws_surface_region_get_area(&a0, qaws_shell_region_get_face(g, 0));
			if (qaws_shell_region_get_face(g, 1)) qaws_surface_region_get_area(&a1, qaws_shell_region_get_face(g, 1));
		}
		sprintf(msg, "split plane: a disk over the curved edge against pi r^2 %.2g and 2 pi r %.2g; the disk's lens %.2g, the rest %.2g",
			fabs((double)A - STR_PI * r1 * r1), fabs((double)P - 2 * STR_PI * r1), fabs((double)a1 - lens), fabs((double)a0 - (STR_PI * r1 * r1 - lens)));
		TEST_ASSERT(st == QAWS_STATUS_OK && fabs((double)A - STR_PI * r1 * r1) <= tol && fabs((double)P - 2 * STR_PI * r1) <= tol
			&& fabs((double)a1 - lens) <= tol && fabs((double)a0 - (STR_PI * r1 * r1 - lens)) <= tol, msg);
		qaws_shell_region_destroy(g);
	}
	qaws_shell_atlas_destroy(a);
	qaws_surface_destroy(sa); qaws_surface_destroy(sb);
	for (i = 0; i < 4; i++) qaws_curve_destroy(oc[i]);
	qaws_curve_destroy(hole); qaws_curve_destroy(disk);
}

static void test_tube(void)
{
	qaws_surface *sa = str_plane(STR_P00, STR_P10, STR_P01, STR_P11), *tube = NULL;
	qaws_curve *prof = str_line(0.5, 0, 0.5, -1.5), *hole = str_arc(0.5, 0.5, 0.125, 0, 2 * STR_PI), *oc[4];
	qaws_curve const *ol[4], *hl[1];
	qaws_surface_revolution_desc rd;
	qaws_path_2d pa[2];
	qaws_shell_face_desc fd[2];
	qaws_shell_atlas* a = NULL;
	qaws_status st;
	unsigned int i;
	char msg[300];
	memset(&rd, 0, sizeof(rd));
	rd.profile = prof; rd.axis_direction.z = 1;
	qaws_surface_create_revolution(&tube, &rd);
	str_square(oc, 0.125, 0.875);
	for (i = 0; i < 4; i++) ol[i] = oc[i];
	hl[0] = hole;
	pa[0].curves = ol; pa[0].curve_count = 4; pa[0].closed = 1;
	pa[1].curves = hl; pa[1].curve_count = 1; pa[1].closed = 1;
	fd[0].surface = sa; fd[0].trims = pa; fd[0].trim_count = 2;
	fd[1].surface = tube; fd[1].trims = NULL; fd[1].trim_count = 0;
	st = qaws_shell_atlas_create_trimmed(&a, fd, 2, 0);
	sprintf(msg, "plate and tube: the hole's trim against the tube's top, %u shared edge", a ? qaws_shell_atlas_get_edge_count(a) : 0);
	TEST_ASSERT(st == QAWS_STATUS_OK && qaws_shell_atlas_get_edge_count(a) == 1, msg);
	if (st == QAWS_STATUS_OK)
	{
		/* rays from (1.2, 0.3) into the hole: on the tube, helices (the
		   components along the rim and down it kept, the cylinder unrolled) */
		static double const dirs[3][2] = { { -1, -0.1 }, { -1, -0.5 }, { -1, 0.15 } };
		double worst = 0;
		unsigned int z, crs = 0;
		for (z = 0; z < 3 && st == QAWS_STATUS_OK; z++)
		{
			qaws_shell_point q;
			qaws_vec2 d;
			qaws_vec3 e;
			unsigned int cr = 0;
			double R = 0.5, x0 = 1.2, y0 = 0.3, n = hypot(dirs[z][0], dirs[z][1]), dx = dirs[z][0] / n, dy = dirs[z][1] / n, L = 1.6;
			double b = x0 * dx + y0 * dy, c = x0 * x0 + y0 * y0 - R * R, t = -b - sqrt(b * b - c), px = x0 + t * dx, py = y0 + t * dy;
			double al = dx * (py / R) + dy * (-px / R), be = -(dx * px + dy * py) / R, l = L - t, ph = atan2(py, px) - al * l / R, off;
			d.x = (qaws_scalar)dx; d.y = (qaws_scalar)dy;
			st = qaws_shell_exp(&q, NULL, &cr, a, str_pt(0, (x0 + 2) / 4, (y0 + 2) / 4), d, (qaws_scalar)L);
			e = str_at(a, q);
			off = sqrt(((double)e.x - R * cos(ph)) * ((double)e.x - R * cos(ph)) + ((double)e.y - R * sin(ph)) * ((double)e.y - R * sin(ph))
				+ ((double)e.z + be * l) * ((double)e.z + be * l));
			if (off > worst) worst = off;
			if (q.face == 1 && cr == 1) crs++;
		}
		sprintf(msg, "plate and tube: three rays into the tube (%u on it after one crossing) against the helices within %.2g", crs, worst);
		TEST_ASSERT(st == QAWS_STATUS_OK && crs == 3 && worst <= (STR_F32 ? 2e-6 : 1e-12), msg);
	}
	qaws_shell_atlas_destroy(a);
	qaws_surface_destroy(sa); qaws_surface_destroy(tube);
	qaws_curve_destroy(prof); qaws_curve_destroy(hole);
	for (i = 0; i < 4; i++) qaws_curve_destroy(oc[i]);
}

/* the cube [-1, 1]^3, faces oriented outward: +z, -z, +x, -x, +y, -y; each
   a plane half as large again, trimmed to its square when trim is set */
static void str_cube(qaws_surface** f, double grow)
{
	static double const P[8][3] = { { -1, -1, -1 }, { 1, -1, -1 }, { -1, 1, -1 }, { 1, 1, -1 }, { -1, -1, 1 }, { 1, -1, 1 }, { -1, 1, 1 }, { 1, 1, 1 } };
	static int const C[6][4] = { { 4, 5, 6, 7 }, { 0, 2, 1, 3 }, { 1, 3, 5, 7 }, { 0, 4, 2, 6 }, { 2, 6, 3, 7 }, { 0, 1, 4, 5 } };
	unsigned int i, k, m;
	for (i = 0; i < 6; i++)
	{
		double q[4][3], c[3];
		for (m = 0; m < 3; m++) c[m] = 0.25 * (P[C[i][0]][m] + P[C[i][1]][m] + P[C[i][2]][m] + P[C[i][3]][m]);
		for (k = 0; k < 4; k++)
			for (m = 0; m < 3; m++) q[k][m] = c[m] + grow * (P[C[i][k]][m] - c[m]);
		f[i] = str_plane(q[0], q[1], q[2], q[3]);
	}
}

static double str_circle(qaws_shell_atlas const* a, qaws_shell_point c, double r, double* per, qaws_status* st)
{
	qaws_shell_curve* k = NULL;
	qaws_shell_region* g = NULL;
	qaws_scalar A = 0, L = 0;
	*st = qaws_shell_curve_create_geodesic_circle(&k, a, c, (qaws_scalar)r, 0);
	if (*st == QAWS_STATUS_OK) *st = qaws_shell_region_create(&g, a, (qaws_shell_curve const* const*)&k, 1);
	if (*st == QAWS_STATUS_OK) *st = qaws_shell_region_get_area(&A, g);
	if (*st == QAWS_STATUS_OK) *st = qaws_shell_region_get_perimeter(&L, g);
	*per = (double)L;
	qaws_shell_region_destroy(g);
	qaws_shell_curve_destroy(k);
	return (double)A;
}

static void test_trimmed_cube(void)
{
	qaws_surface *f[6], *g[6];
	qaws_curve* sq[6][4];
	qaws_curve const* sl[6][4];
	qaws_path_2d loops[6];
	qaws_shell_face_desc fd[6];
	qaws_shell_atlas *a = NULL, *b = NULL;
	qaws_status st, sb;
	double tol = STR_F32 ? 2e-4 : 1e-9, lo = 1.0 / 6, span = 2.0 / 3;
	unsigned int i, k;
	char msg[300];
	str_cube(f, 1.5);
	str_cube(g, 1);
	for (i = 0; i < 6; i++)
	{
		str_square(sq[i], lo, 1 - lo);
		for (k = 0; k < 4; k++) sl[i][k] = sq[i][k];
		loops[i].curves = sl[i]; loops[i].curve_count = 4; loops[i].closed = 1;
		fd[i].surface = f[i]; fd[i].trims = &loops[i]; fd[i].trim_count = 1;
	}
	st = qaws_shell_atlas_create_trimmed(&a, fd, 6, 0);
	sb = qaws_shell_atlas_create(&b, (qaws_surface const* const*)g, 6, 0);
	sprintf(msg, "trimmed cube: planes trimmed to the faces, %u shared edges", a ? qaws_shell_atlas_get_edge_count(a) : 0);
	TEST_ASSERT(st == QAWS_STATUS_OK && sb == QAWS_STATUS_OK && qaws_shell_atlas_get_edge_count(a) == 12, msg);
	if (st == QAWS_STATUS_OK)
	{
		/* (a face's (u, v) of the cube at lo + span (u, v) of the plane) */
		double per1 = 0, per2 = 0, a1, a2, r1 = 0.8, r2 = 0.5;
		qaws_status s1, s2;
		a1 = str_circle(a, str_pt(0, lo + span * 0.75, lo + span * 0.5), r1, &per1, &s1);
		a2 = str_circle(a, str_pt(0, lo + span * 0.8, lo + span * 0.8), r2, &per2, &s2);
		sprintf(msg, "trimmed cube: circles over one edge and over two by a vertex: areas against pi r^2 %.2g, %.2g; perimeters against 2 pi r %.2g, %.2g",
			fabs(a1 - STR_PI * r1 * r1), fabs(a2 - STR_PI * r2 * r2), fabs(per1 - 2 * STR_PI * r1), fabs(per2 - 2 * STR_PI * r2));
		TEST_ASSERT(s1 == QAWS_STATUS_OK && s2 == QAWS_STATUS_OK && fabs(a1 - STR_PI * r1 * r1) <= tol && fabs(a2 - STR_PI * r2 * r2) <= tol &&
			fabs(per1 - 2 * STR_PI * r1) <= tol && fabs(per2 - 2 * STR_PI * r2) <= tol, msg);
	}
	if (st == QAWS_STATUS_OK && sb == QAWS_STATUS_OK)
	{
		/* a disk round a corner: the trimmed cube's against the plain one's */
		qaws_shell_region *ga = NULL, *gb = NULL;
		qaws_scalar A = 0, B = 0;
		qaws_status s1, s2;
		double t = 0.3;
		s1 = qaws_shell_region_create_geodesic_disk(&ga, a, str_pt(0, lo + span * (1 - t / 2), lo + span * (1 - t / 2)), (qaws_scalar)0.75, 0);
		s2 = qaws_shell_region_create_geodesic_disk(&gb, b, str_pt(0, 1 - t / 2, 1 - t / 2), (qaws_scalar)0.75, 0);
		if (s1 == QAWS_STATUS_OK) qaws_shell_region_get_area(&A, ga);
		if (s2 == QAWS_STATUS_OK) qaws_shell_region_get_area(&B, gb);
		sprintf(msg, "trimmed cube: a disk round a corner %.9f against the untrimmed cube's %.9f (%.2g)", (double)A, (double)B, fabs((double)A - B));
		TEST_ASSERT(s1 == QAWS_STATUS_OK && s2 == QAWS_STATUS_OK && fabs((double)A - B) <= (STR_F32 ? 2e-4 : 1e-8) && A > 1.5, msg);
		qaws_shell_region_destroy(ga);
		qaws_shell_region_destroy(gb);
	}
	qaws_shell_atlas_destroy(a);
	qaws_shell_atlas_destroy(b);
	for (i = 0; i < 6; i++)
	{
		qaws_surface_destroy(f[i]); qaws_surface_destroy(g[i]);
		for (k = 0; k < 4; k++) qaws_curve_destroy(sq[i][k]);
	}
}

/* a T-junction: the side y = 0 of B (x in [0, 2], y in [-1, 0]) meets two
   unit squares, side by side, flat beyond it or standing up from it as
   walls (the angle round the T point 2 pi either way: unfolded, a plane) */
static void test_tee(int walls)
{
	double B00[3] = { 0, -1, 0 }, B10[3] = { 2, -1, 0 }, B01[3] = { 0, 0, 0 }, B11[3] = { 2, 0, 0 };
	double A00[3] = { 0, 0, 0 }, A10[3] = { 1, 0, 0 }, A01[3] = { 0, 1, 0 }, A11[3] = { 1, 1, 0 };
	double C00[3] = { 1, 0, 0 }, C10[3] = { 2, 0, 0 }, C01[3] = { 1, 1, 0 }, C11[3] = { 2, 1, 0 };
	static double const cs[3][3] = { { 1.0, -0.3, 0.5 }, { 0.7, -0.25, 0.5 }, { 0.6, -0.4, 0.55 } };
	qaws_surface* f[3];
	qaws_shell_atlas* a = NULL;
	qaws_status st;
	double worst = 0, tol = STR_F32 ? 2e-4 : 1e-9;
	unsigned int i, z, ok = 0;
	char msg[300];
	if (walls) { A01[1] = 0; A01[2] = 1; A11[1] = 0; A11[2] = 1; C01[1] = 0; C01[2] = 1; C11[1] = 0; C11[2] = 1; }
	f[0] = str_plane(B00, B10, B01, B11);
	f[1] = str_plane(A00, A10, A01, A11);
	f[2] = str_plane(C00, C10, C01, C11);
	st = qaws_shell_atlas_create(&a, (qaws_surface const* const*)f, 3, 0);
	sprintf(msg, "T-junction (%s): a side meeting two faces' sides, split there: %u shared edges", walls ? "walls" : "flat",
		a ? qaws_shell_atlas_get_edge_count(a) : 0);
	TEST_ASSERT(st == QAWS_STATUS_OK && qaws_shell_atlas_get_edge_count(a) == 3, msg);
	for (z = 0; z < 3 && st == QAWS_STATUS_OK; z++)
	{
		qaws_shell_region* g = NULL;
		qaws_scalar A = 0, P = 0;
		double r = cs[z][2], e;
		if (qaws_shell_region_create_geodesic_disk(&g, a, str_pt(0, cs[z][0] / 2, cs[z][1] + 1), (qaws_scalar)r, 0) == QAWS_STATUS_OK)
		{
			qaws_shell_region_get_area(&A, g);
			qaws_shell_region_get_perimeter(&P, g);
			e = fabs((double)A - STR_PI * r * r) + fabs((double)P - 2 * STR_PI * r);
			if (e > worst) worst = e;
			if (e <= tol) ok++;
		}
		else
			worst = 1e300;
		qaws_shell_region_destroy(g);
	}
	sprintf(msg, "T-junction (%s): three disks (two round the T point) against pi r^2 and 2 pi r within %.2g", walls ? "walls" : "flat", worst);
	TEST_ASSERT(st == QAWS_STATUS_OK && ok == 3, msg);
	if (walls && st == QAWS_STATUS_OK)
	{
		/* a disk on a wall just touching the edge below it: on that wall only */
		qaws_shell_region* g = NULL;
		qaws_scalar A = 0;
		unsigned int nf = 0;
		qaws_status s1 = qaws_shell_region_create_geodesic_disk(&g, a, str_pt(2, 0.5, 0.45), (qaws_scalar)0.45, 0);
		if (s1 == QAWS_STATUS_OK)
		{
			qaws_shell_region_get_area(&A, g);
			for (i = 0; i < 3; i++) if (qaws_shell_region_get_face(g, i)) nf++;
		}
		sprintf(msg, "T-junction (walls): a disk touching the edge below it, on %u face, against pi r^2 %.2g", nf, fabs((double)A - STR_PI * 0.45 * 0.45));
		TEST_ASSERT(s1 == QAWS_STATUS_OK && nf == 1 && fabs((double)A - STR_PI * 0.45 * 0.45) <= tol, msg);
		qaws_shell_region_destroy(g);
	}
	qaws_shell_atlas_destroy(a);
	for (i = 0; i < 3; i++) qaws_surface_destroy(f[i]);
}

/* the cylinder of radius 1/2 and height pi (its (u, v) isotropic: pi to the
   unit), trimmed by its rims (lines wrapping round once) and a hole of
   radius 0.1 centred on the seam; the hole filled by a disk of the same
   surface whose circle starts elsewhere in another copy */
static void test_cylinder(void)
{
	qaws_curve *prof = str_line(0.5, 0, 0.5, -STR_PI), *bot = str_line(0, 0, 1, 0), *top = str_line(1, 1, 0, 1);
	qaws_curve *hole = str_arc(0, 0.5, 0.1, 0, 2 * STR_PI), *disk = str_arc(1, 0.5, 0.1, 2, 2 + 2 * STR_PI);
	qaws_curve const *l1[1], *l2[1], *l3[1], *l4[1];
	qaws_surface* s = NULL;
	qaws_surface_revolution_desc rd;
	qaws_path_2d pa[3], pb[1];
	qaws_shell_face_desc fd[2];
	qaws_shell_atlas* a = NULL;
	qaws_status st;
	static double const cs[3][3] = { { 0.85, 0.5, 0.3 }, { 0.9, 0.3, 0.25 }, { 0.5, 0.5, 0.4 } };
	double worst = 0, tol = STR_F32 ? 2e-4 : 1e-8;
	unsigned int z, ok = 0;
	char msg[300];
	memset(&rd, 0, sizeof(rd));
	rd.profile = prof; rd.axis_direction.z = 1;
	qaws_surface_create_revolution(&s, &rd);
	l1[0] = bot; l2[0] = top; l3[0] = hole; l4[0] = disk;
	pa[0].curves = l1; pa[0].curve_count = 1; pa[0].closed = 1;
	pa[1].curves = l2; pa[1].curve_count = 1; pa[1].closed = 1;
	pa[2].curves = l3; pa[2].curve_count = 1; pa[2].closed = 1;
	pb[0].curves = l4; pb[0].curve_count = 1; pb[0].closed = 1;
	fd[0].surface = s; fd[0].trims = pa; fd[0].trim_count = 3;
	fd[1].surface = s; fd[1].trims = pb; fd[1].trim_count = 1;
	st = qaws_shell_atlas_create_trimmed(&a, fd, 2, 0);
	sprintf(msg, "cylinder: trimmed by its rims and a hole on the seam, the hole's disk: %u shared edges", a ? qaws_shell_atlas_get_edge_count(a) : 0);
	TEST_ASSERT(st == QAWS_STATUS_OK && qaws_shell_atlas_get_edge_count(a) == 2, msg);
	for (z = 0; z < 3 && st == QAWS_STATUS_OK; z++)
	{
		qaws_shell_region* g = NULL;
		qaws_scalar A = 0, P = 0;
		double r = cs[z][2], e = 1e300;
		if (qaws_shell_region_create_geodesic_disk(&g, a, str_pt(0, cs[z][0], cs[z][1]), (qaws_scalar)r, 0) == QAWS_STATUS_OK)
		{
			qaws_shell_region_get_area(&A, g);
			qaws_shell_region_get_perimeter(&P, g);
			e = fabs((double)A - STR_PI * r * r) + fabs((double)P - 2 * STR_PI * r);
		}
		if (e > worst) worst = e;
		if (e <= tol) ok++;
		qaws_shell_region_destroy(g);
	}
	sprintf(msg, "cylinder: three disks (one over the hole and the seam) against pi r^2 and 2 pi r within %.2g", worst);
	TEST_ASSERT(st == QAWS_STATUS_OK && ok == 3, msg);
	if (st == QAWS_STATUS_OK)
	{
		/* a ray from (0.7, 0.45) over the hole's disk and the seam: unrolled, a line */
		qaws_shell_point q;
		qaws_vec2 d;
		qaws_vec3 e;
		unsigned int cr = 0;
		double dx = 1, dy = 0.15, n = hypot(dx, dy), L = 1.6, X = STR_PI * 0.7 + L * dx / n, Y = STR_PI * 0.45 + L * dy / n, off;
		double ex = 0.5 * cos(2 * X), ey = 0.5 * sin(2 * X), ez = -Y;
		d.x = (qaws_scalar)(dx / n); d.y = (qaws_scalar)(dy / n);
		st = qaws_shell_exp(&q, NULL, &cr, a, str_pt(0, 0.7, 0.45), d, (qaws_scalar)L);
		e = str_at(a, q);
		off = sqrt(((double)e.x - ex) * ((double)e.x - ex) + ((double)e.y - ey) * ((double)e.y - ey) + ((double)e.z - ez) * ((double)e.z - ez));
		sprintf(msg, "cylinder: a ray over the hole's disk and the seam (%u crossings) against the unrolled line %.2g", cr, off);
		TEST_ASSERT(st == QAWS_STATUS_OK && cr == 2 && off <= (STR_F32 ? 2e-6 : 1e-12), msg);
	}
	qaws_shell_atlas_destroy(a);
	qaws_surface_destroy(s);
	qaws_curve_destroy(prof); qaws_curve_destroy(bot); qaws_curve_destroy(top); qaws_curve_destroy(hole); qaws_curve_destroy(disk);
}

/* a tube (closed in u) capped by one disk, or by two half disks whose
   diameter meets the tube's rim at two T points (the rim, a side along a
   closed face, split there): the same surface either way */
static void test_capped_tube(void)
{
	static double const C00[3] = { -1, -1, 0 }, C10[3] = { -1, 1, 0 }, C01[3] = { 1, -1, 0 }, C11[3] = { 1, 1, 0 };
	static double const cs[3][2] = { { 0.0, 0.15 }, { 0.25, 0.15 }, { 0.6, 0.2 } };
	qaws_curve *prof = str_line(0.5, 0, 0.5, -1.5), *full = str_arc(0.5, 0.5, 0.25, STR_PI / 2, 5 * STR_PI / 2);
	qaws_curve *h1 = str_arc(0.5, 0.5, 0.25, -STR_PI / 2, STR_PI / 2), *d1 = str_line(0.5, 0.75, 0.5, 0.25);
	qaws_curve *h2 = str_arc(0.5, 0.5, 0.25, STR_PI / 2, 3 * STR_PI / 2), *d2 = str_line(0.5, 0.25, 0.5, 0.75);
	qaws_surface *tube = NULL, *cap = str_plane(C00, C10, C01, C11);
	qaws_surface_revolution_desc rd;
	qaws_curve const *lf[1], *la[2], *lb[2];
	qaws_path_2d pf[1], pa[1], pb[1];
	qaws_shell_face_desc A[2], B[3];
	qaws_shell_atlas *aa = NULL, *ab = NULL;
	qaws_status s1, s2;
	double worst = 0, ref = -1, tol = STR_F32 ? 2e-4 : 1e-9;
	unsigned int z, ok = 0;
	char msg[300];
	memset(&rd, 0, sizeof(rd));
	rd.profile = prof; rd.axis_direction.z = 1;
	qaws_surface_create_revolution(&tube, &rd);
	lf[0] = full; la[0] = h1; la[1] = d1; lb[0] = h2; lb[1] = d2;
	pf[0].curves = lf; pf[0].curve_count = 1; pf[0].closed = 1;
	pa[0].curves = la; pa[0].curve_count = 2; pa[0].closed = 1;
	pb[0].curves = lb; pb[0].curve_count = 2; pb[0].closed = 1;
	memset(A, 0, sizeof(A));
	memset(B, 0, sizeof(B));
	A[0].surface = tube; A[1].surface = cap; A[1].trims = pf; A[1].trim_count = 1;
	B[0].surface = tube; B[1].surface = cap; B[1].trims = pa; B[1].trim_count = 1; B[2].surface = cap; B[2].trims = pb; B[2].trim_count = 1;
	s1 = qaws_shell_atlas_create_trimmed(&aa, A, 2, 0);
	s2 = qaws_shell_atlas_create_trimmed(&ab, B, 3, 0);
	sprintf(msg, "capped tube: one cap, %u shared edge; two half caps meeting the rim at T points, %u shared edges",
		aa ? qaws_shell_atlas_get_edge_count(aa) : 0, ab ? qaws_shell_atlas_get_edge_count(ab) : 0);
	TEST_ASSERT(s1 == QAWS_STATUS_OK && s2 == QAWS_STATUS_OK && qaws_shell_atlas_get_edge_count(aa) == 1 && qaws_shell_atlas_get_edge_count(ab) == 3, msg);
	/* disks over the rim (one centred on the seam right below a T point):
	   the same either way, and by symmetry the same wherever round */
	/* (by default the two disks against the first: the third, one cap against two again, under QAWS_TEST_HEAVY) */
	for (z = 0; z < (test_heavy() ? 3u : 2u) && s1 == QAWS_STATUS_OK && s2 == QAWS_STATUS_OK; z++)
	{
		qaws_shell_region *g = NULL, *h = NULL;
		qaws_scalar Ag = 0, Ah = 0, Pg = 0, Ph = 0;
		double e;
		if (qaws_shell_region_create_geodesic_disk(&g, aa, str_pt(0, cs[z][0], cs[z][1]), (qaws_scalar)0.5, 0) == QAWS_STATUS_OK &&
			qaws_shell_region_create_geodesic_disk(&h, ab, str_pt(0, cs[z][0], cs[z][1]), (qaws_scalar)0.5, 0) == QAWS_STATUS_OK)
		{
			qaws_shell_region_get_area(&Ag, g); qaws_shell_region_get_perimeter(&Pg, g);
			qaws_shell_region_get_area(&Ah, h); qaws_shell_region_get_perimeter(&Ph, h);
			if (z == 0) ref = (double)Ag + (double)Pg;
			e = fabs((double)Ag - Ah) + fabs((double)Pg - Ph);
			if (z < 2) e += fabs((double)Ag + (double)Pg - ref);
			if (e > worst) worst = e;
			if (e <= tol) ok++;
		}
		else
			worst = 1e300;
		qaws_shell_region_destroy(g);
		qaws_shell_region_destroy(h);
	}
	sprintf(msg, "capped tube: disks over the rim (one on the seam below a T point) alike with one cap or two, round the tube alike, within %.2g", worst);
	TEST_ASSERT(ok == (test_heavy() ? 3u : 2u), msg);
	qaws_shell_atlas_destroy(aa);
	qaws_shell_atlas_destroy(ab);
	qaws_surface_destroy(tube); qaws_surface_destroy(cap);
	qaws_curve_destroy(prof); qaws_curve_destroy(full); qaws_curve_destroy(h1); qaws_curve_destroy(d1); qaws_curve_destroy(h2); qaws_curve_destroy(d2);
}

/* a torus (closed both ways) with a hole where its seams cross, filled by a
   disk of the same surface: the loops as given (the face on their left);
   disks over the hole and rays against the whole torus's */
static void test_torus(void)
{
	static double const cs[3][3] = { { 0.15, 0.1, 0.3 }, { -0.13, 0.08, 0.25 }, { 0.03, -0.17, 0.3 } };
	qaws_curve *prof = str_arc(1.0, 0, 0.35, 0, 2 * STR_PI), *hole = str_arc(0, 0, 0.12, 2 * STR_PI, 0), *disk = str_arc(0, 0, 0.12, 1.0, 1.0 + 2 * STR_PI);
	qaws_surface* s = NULL;
	qaws_surface_revolution_desc rd;
	qaws_curve const *lh[1], *ld[1];
	qaws_path_2d ph[1], pd[1];
	qaws_shell_face_desc A[1], B[2];
	qaws_shell_atlas *aa = NULL, *ab = NULL;
	qaws_status s1, s2;
	double worst = 0, tol = STR_F32 ? 2e-4 : 1e-9;
	unsigned int z, ok = 0;
	char msg[300];
	memset(&rd, 0, sizeof(rd));
	rd.profile = prof; rd.axis_direction.z = 1;
	qaws_surface_create_revolution(&s, &rd);
	lh[0] = hole; ld[0] = disk;
	ph[0].curves = lh; ph[0].curve_count = 1; ph[0].closed = 1;
	pd[0].curves = ld; pd[0].curve_count = 1; pd[0].closed = 1;
	memset(A, 0, sizeof(A));
	memset(B, 0, sizeof(B));
	A[0].surface = s;
	B[0].surface = s; B[0].trims = ph; B[0].trim_count = 1;
	B[1].surface = s; B[1].trims = pd; B[1].trim_count = 1;
	s1 = qaws_shell_atlas_create_trimmed(&aa, A, 1, 0);
	s2 = qaws_shell_atlas_create_trimmed(&ab, B, 2, 0);
	sprintf(msg, "torus: whole, %u edges; with a hole on both seams and its disk, %u shared edges",
		aa ? qaws_shell_atlas_get_edge_count(aa) : 0, ab ? qaws_shell_atlas_get_edge_count(ab) : 0);
	TEST_ASSERT(s1 == QAWS_STATUS_OK && s2 == QAWS_STATUS_OK && qaws_shell_atlas_get_edge_count(aa) == 0 && qaws_shell_atlas_get_edge_count(ab) == 2, msg);
	for (z = 0; z < 3 && s1 == QAWS_STATUS_OK && s2 == QAWS_STATUS_OK; z++)
	{
		qaws_shell_region *g = NULL, *h = NULL;
		qaws_scalar Ag = 0, Ah = 0, Pg = 0, Ph = 0;
		double e = 1e300;
		if (qaws_shell_region_create_geodesic_disk(&g, aa, str_pt(0, cs[z][0], cs[z][1]), (qaws_scalar)cs[z][2], 0) == QAWS_STATUS_OK &&
			qaws_shell_region_create_geodesic_disk(&h, ab, str_pt(0, cs[z][0], cs[z][1]), (qaws_scalar)cs[z][2], 0) == QAWS_STATUS_OK)
		{
			qaws_shell_region_get_area(&Ag, g); qaws_shell_region_get_perimeter(&Pg, g);
			qaws_shell_region_get_area(&Ah, h); qaws_shell_region_get_perimeter(&Ph, h);
			e = fabs((double)Ag - Ah) + fabs((double)Pg - Ph);
		}
		if (e > worst) worst = e;
		if (e <= tol) ok++;
		qaws_shell_region_destroy(g);
		qaws_shell_region_destroy(h);
	}
	sprintf(msg, "torus: three disks over the hole (across both seams) against the whole torus's within %.2g", worst);
	TEST_ASSERT(ok == 3, msg);
	if (s1 == QAWS_STATUS_OK && s2 == QAWS_STATUS_OK)
	{
		/* a ray over the hole: as on the whole torus */
		qaws_shell_point qa, qb;
		qaws_vec2 d;
		qaws_vec3 ea, eb;
		unsigned int ca = 0, cb = 0;
		qaws_status t1, t2;
		double off;
		d.x = (qaws_scalar)-0.03; d.y = (qaws_scalar)-1;
		t1 = qaws_shell_exp(&qa, NULL, &ca, aa, str_pt(0, 0.08, 0.6), d, (qaws_scalar)1.2);
		t2 = qaws_shell_exp(&qb, NULL, &cb, ab, str_pt(0, 0.08, 0.6), d, (qaws_scalar)1.2);
		ea = str_at(aa, qa); eb = str_at(ab, qb);
		off = sqrt(((double)ea.x - eb.x) * ((double)ea.x - eb.x) + ((double)ea.y - eb.y) * ((double)ea.y - eb.y) + ((double)ea.z - eb.z) * ((double)ea.z - eb.z));
		sprintf(msg, "torus: a ray over the hole (%u crossings) against the whole torus's %.2g", cb, off);
		TEST_ASSERT(t1 == QAWS_STATUS_OK && t2 == QAWS_STATUS_OK && cb == 1 && off <= (STR_F32 ? 2e-6 : 1e-12), msg);
	}
	qaws_shell_atlas_destroy(aa);
	qaws_shell_atlas_destroy(ab);
	qaws_surface_destroy(s);
	qaws_curve_destroy(prof); qaws_curve_destroy(hole); qaws_curve_destroy(disk);
}

/* the taut length from a to b round the circle of radius R at the origin */
static double str_taut(double ax, double ay, double bx, double by, double R)
{
	double da = hypot(ax, ay), db = hypot(bx, by), d = fabs(atan2(ay, ax) - atan2(by, bx)), cr;
	if (d > STR_PI) d = 2 * STR_PI - d;
	cr = fabs(ax * by - ay * bx) / hypot(bx - ax, by - ay);
	if (cr >= R || (ax * (bx - ax) + ay * (by - ay)) * (bx * (bx - ax) + by * (by - ay)) > 0)
		return hypot(bx - ax, by - ay);
	return sqrt(da * da - R * R) + sqrt(db * db - R * R) + R * (d - acos(R / da) - acos(R / db));
}

static void test_shortest_path(void)
{
	/* a plate with a hole of radius 0.5 at the origin */
	static double const pts[4][4] = { { -1.0, 0.05, 1.0, -0.1 }, { -1.2, 0.3, 0.9, 0.2 }, { 0.0, -1.0, 0.1, 1.2 }, { -1.0, 0.6, 1.0, 0.6 } };
	qaws_surface* sa = str_plane(STR_P00, STR_P10, STR_P01, STR_P11);
	qaws_curve *oc[4], *hole = str_arc(0.5, 0.5, 0.125, 0, 2 * STR_PI);
	qaws_curve const *ol[4], *hl[1];
	qaws_path_2d pa[2];
	qaws_shell_face_desc fd;
	qaws_shell_atlas* a = NULL;
	qaws_status st;
	double worst = 0, tol = STR_F32 ? 1e-5 : 1e-9;
	unsigned int i, z, bad = 0;
	char msg[200];
	str_square(oc, 0.125, 0.875);
	for (i = 0; i < 4; i++) ol[i] = oc[i];
	hl[0] = hole;
	pa[0].curves = ol; pa[0].curve_count = 4; pa[0].closed = 1;
	pa[1].curves = hl; pa[1].curve_count = 1; pa[1].closed = 1;
	fd.surface = sa; fd.trims = pa; fd.trim_count = 2;
	st = qaws_shell_atlas_create_trimmed(&a, &fd, 1, 0);
	for (z = 0; z < 4 && st == QAWS_STATUS_OK; z++)
	{
		qaws_shell_curve* c = NULL;
		qaws_scalar L = 0;
		double off;
		if (qaws_shell_curve_create_shortest_path(&c, a, str_pt(0, (pts[z][0] + 2) / 4, (pts[z][1] + 2) / 4), str_pt(0, (pts[z][2] + 2) / 4, (pts[z][3] + 2) / 4), 0) != QAWS_STATUS_OK)
		{
			bad++;
			continue;
		}
		qaws_shell_curve_get_length(&L, c);
		off = fabs((double)L - str_taut(pts[z][0], pts[z][1], pts[z][2], pts[z][3], 0.5));
		if (off > worst) worst = off;
		qaws_shell_curve_destroy(c);
	}
	sprintf(msg, "shortest paths round a hole (three round it, one clear of it): %u failed, against the taut string %.2g", bad, worst);
	TEST_ASSERT(st == QAWS_STATUS_OK && bad == 0 && worst <= tol, msg);
	qaws_shell_atlas_destroy(a);
	qaws_surface_destroy(sa);
	for (i = 0; i < 4; i++) qaws_curve_destroy(oc[i]);
	qaws_curve_destroy(hole);
}

/* a plate with holes of radius r (space) at (xs[i], 0), y = 0 their line */
static qaws_shell_atlas* str_holes(qaws_surface* sa, qaws_curve** oc, qaws_curve** hc, double const* xs, unsigned int n, double r)
{
	qaws_curve const *ol[4], *hl[4];
	qaws_path_2d pa[5];
	qaws_shell_face_desc fd;
	qaws_shell_atlas* a = NULL;
	unsigned int i;
	str_square(oc, 0.125, 0.875);
	for (i = 0; i < 4; i++) ol[i] = oc[i];
	pa[0].curves = ol; pa[0].curve_count = 4; pa[0].closed = 1;
	for (i = 0; i < n; i++)
	{
		hc[i] = str_arc((xs[i] + 2) / 4, 0.5, r / 4, 0, 2 * STR_PI);
		hl[i] = hc[i];
		pa[1 + i].curves = &hl[i]; pa[1 + i].curve_count = 1; pa[1 + i].closed = 1;
	}
	fd.surface = sa; fd.trims = pa; fd.trim_count = 1 + n;
	qaws_shell_atlas_create_trimmed(&a, &fd, 1, 0);
	return a;
}

static void test_shortest_holes(void)
{
	/* holes in a row, the ends on their line beyond them: over their tops,
	   the tangents from the ends, the arcs to the tops and the segments
	   between the tops */
	static double const x2[2] = { -0.6, 0.6 }, x3[3] = { -0.8, 0, 0.8 };
	double worst = 0, tol = STR_F32 ? 1e-5 : 1e-9;
	unsigned int z, i, bad = 0;
	char msg[200];
	for (z = 0; z < 2; z++)
	{
		qaws_surface* sa = str_plane(STR_P00, STR_P10, STR_P01, STR_P11);
		qaws_curve *oc[4], *hc[3];
		unsigned int n = z ? 3 : 2;
		double r = z ? 0.25 : 0.3, e = z ? 1.5 : 1.4, d = e - (z ? 0.8 : 0.6), span = z ? 1.6 : 1.2, ref, L0 = 0;
		qaws_shell_atlas* a = str_holes(sa, oc, hc, z ? x3 : x2, n, r);
		qaws_shell_curve* c = NULL;
		ref = 2 * sqrt(d * d - r * r) + 2 * r * (STR_PI / 2 - acos(r / d)) + span;
		if (a && qaws_shell_curve_create_shortest_path(&c, a, str_pt(0, (-e + 2) / 4, 0.5), str_pt(0, (e + 2) / 4, 0.5), 0) == QAWS_STATUS_OK)
		{
			qaws_scalar L = 0;
			qaws_shell_curve_get_length(&L, c);
			L0 = (double)L;
			if (fabs(L0 - ref) > worst) worst = fabs(L0 - ref);
		}
		else
			bad++;
		qaws_shell_curve_destroy(c);
		qaws_shell_atlas_destroy(a);
		qaws_surface_destroy(sa);
		for (i = 0; i < 4; i++) qaws_curve_destroy(oc[i]);
		for (i = 0; i < n; i++) qaws_curve_destroy(hc[i]);
	}
	sprintf(msg, "shortest paths round two and three holes in a row: %u failed, against the closed form %.2g", bad, worst);
	TEST_ASSERT(bad == 0 && worst <= tol, msg);
}

static void test_shortest_other_face(void)
{
	/* an end on a neighbouring face: a plate x in [-1, 1] with a hole of
	   radius 0.3 at the origin and a plate x in [-3, -1] beside it, in the
	   same plane (unfolded, the taut string of the plane) */
	static double const A00[3] = { -1.5, -1.5, 0 }, A10[3] = { 1.5, -1.5, 0 }, A01[3] = { -1.5, 1.5, 0 }, A11[3] = { 1.5, 1.5, 0 };
	static double const B00[3] = { -3, -1, 0 }, B10[3] = { -1, -1, 0 }, B01[3] = { -3, 1, 0 }, B11[3] = { -1, 1, 0 };
	qaws_surface *sa = str_plane(A00, A10, A01, A11), *sb = str_plane(B00, B10, B01, B11);
	qaws_curve *oc[4], *hole = str_arc(0.5, 0.5, 0.1, 0, 2 * STR_PI);
	qaws_curve const *ol[4], *hl[1];
	qaws_path_2d pa[2];
	qaws_shell_face_desc fd[2];
	qaws_shell_atlas* a = NULL;
	qaws_shell_curve* c = NULL;
	qaws_scalar L = 0;
	qaws_status st;
	double ref = str_taut(-1.5, 0.05, 0.8, -0.1, 0.3), tol = STR_F32 ? 1e-4 : 1e-6;
	unsigned int i;
	char msg[200];
	str_square(oc, 1.0 / 6, 5.0 / 6);
	for (i = 0; i < 4; i++) ol[i] = oc[i];
	hl[0] = hole;
	pa[0].curves = ol; pa[0].curve_count = 4; pa[0].closed = 1;
	pa[1].curves = hl; pa[1].curve_count = 1; pa[1].closed = 1;
	fd[0].surface = sa; fd[0].trims = pa; fd[0].trim_count = 2;
	fd[1].surface = sb; fd[1].trims = NULL; fd[1].trim_count = 0;
	st = qaws_shell_atlas_create_trimmed(&a, fd, 2, 0);
	if (st == QAWS_STATUS_OK)
		st = qaws_shell_curve_create_shortest_path(&c, a, str_pt(1, 0.75, 0.525), str_pt(0, (0.8 + 1.5) / 3, (-0.1 + 1.5) / 3), 0);
	if (st == QAWS_STATUS_OK) qaws_shell_curve_get_length(&L, c);
	sprintf(msg, "shortest path from a neighbouring face round the hole (%d): against the unfolded taut string %.2g", (int)st, fabs((double)L - ref));
	TEST_ASSERT(st == QAWS_STATUS_OK && fabs((double)L - ref) <= tol, msg);
	qaws_shell_curve_destroy(c);
	qaws_shell_atlas_destroy(a);
	qaws_surface_destroy(sa); qaws_surface_destroy(sb);
	for (i = 0; i < 4; i++) qaws_curve_destroy(oc[i]);
	qaws_curve_destroy(hole);
}

static void test_thin_holes(void)
{
	/* holes far smaller than a trace step: a ray through each stops at it,
	   one passing beside it does not */
	static double const rs[3] = { 0.004, 0.0015, 0.0005 };
	unsigned int z, i, stopped = 0, passed = 0;
	char msg[200];
	for (z = 0; z < 3; z++)
	{
		qaws_surface* sa = str_plane(STR_P00, STR_P10, STR_P01, STR_P11);
		qaws_curve *oc[4], *hole = str_arc(0.5, 0.5, rs[z], 0, 2 * STR_PI);
		qaws_curve const *ol[4], *hl[1];
		qaws_path_2d pa[2];
		qaws_shell_face_desc fd;
		qaws_shell_atlas* a = NULL;
		qaws_shell_point q;
		qaws_vec2 d;
		str_square(oc, 0.125, 0.875);
		for (i = 0; i < 4; i++) ol[i] = oc[i];
		hl[0] = hole;
		pa[0].curves = ol; pa[0].curve_count = 4; pa[0].closed = 1;
		pa[1].curves = hl; pa[1].curve_count = 1; pa[1].closed = 1;
		fd.surface = sa; fd.trims = pa; fd.trim_count = 2;
		if (qaws_shell_atlas_create_trimmed(&a, &fd, 1, 0) == QAWS_STATUS_OK)
		{
			d.x = 1; d.y = 0;
			stopped += qaws_shell_exp(&q, NULL, NULL, a, str_pt(0, 0.3, 0.5), d, (qaws_scalar)1.6) == QAWS_STATUS_OUT_OF_RANGE;
			d.y = (qaws_scalar)(2.5 * rs[z] / 0.2);
			passed += qaws_shell_exp(&q, NULL, NULL, a, str_pt(0, 0.3, 0.5), d, (qaws_scalar)1.6) == QAWS_STATUS_OK;
		}
		qaws_shell_atlas_destroy(a);
		qaws_surface_destroy(sa);
		for (i = 0; i < 4; i++) qaws_curve_destroy(oc[i]);
		qaws_curve_destroy(hole);
	}
	sprintf(msg, "thin holes (0.016, 0.006 and 0.002 across, a step 0.2): rays through them stopped %u/3, beside them passing %u/3", stopped, passed);
	TEST_ASSERT(stopped == 3 && passed == 3, msg);
}

static void test_refusals(void)
{
	qaws_surface* sa = str_plane(STR_P00, STR_P10, STR_P01, STR_P11);
	qaws_curve *open3[3], *prof = str_arc(1.0, 0, 0.3, 0, 2 * STR_PI), *loop = str_arc(0.5, 0.5, 0.2, 0, 2 * STR_PI);
	qaws_curve const *ol[3], *cl[1];
	qaws_path_2d p, pc;
	qaws_shell_face_desc fd;
	qaws_shell_atlas* a = NULL;
	qaws_surface* tube = NULL;
	qaws_surface_revolution_desc rd;
	qaws_status s1, s2;
	unsigned int i;
	char msg[200];
	open3[0] = str_line(0.2, 0.2, 0.8, 0.2); open3[1] = str_line(0.8, 0.2, 0.8, 0.8); open3[2] = str_line(0.8, 0.8, 0.3, 0.2);
	for (i = 0; i < 3; i++) ol[i] = open3[i];
	p.curves = ol; p.curve_count = 3; p.closed = 1;
	fd.surface = sa; fd.trims = &p; fd.trim_count = 1;
	s1 = qaws_shell_atlas_create_trimmed(&a, &fd, 1, 0);
	qaws_shell_atlas_destroy(a);
	a = NULL;
	memset(&rd, 0, sizeof(rd));
	rd.profile = prof; rd.axis_direction.z = 1;
	qaws_surface_create_revolution(&tube, &rd);
	cl[0] = loop;
	pc.curves = cl; pc.curve_count = 1; pc.closed = 1;
	fd.surface = tube; fd.trims = &pc; fd.trim_count = 1;
	s2 = qaws_shell_atlas_create_trimmed(&a, &fd, 1, 0);
	{
		/* a disk of a torus: a centre off it refused */
		qaws_shell_region* g = NULL;
		qaws_status s3 = s2 == QAWS_STATUS_OK ? qaws_shell_region_create_geodesic_disk(&g, a, str_pt(0, 0.1, 0.1), (qaws_scalar)0.1, 0) : QAWS_STATUS_OK;
		qaws_shell_region_destroy(g);
		qaws_shell_atlas_destroy(a);
		sprintf(msg, "refused: a loop that does not close (%d); a disk of a torus taken (%d), a centre off it (%d)", (int)s1, (int)s2, (int)s3);
		TEST_ASSERT(s1 == QAWS_STATUS_INVALID_ARGUMENT && s2 == QAWS_STATUS_OK && s3 == QAWS_STATUS_INVALID_ARGUMENT, msg);
	}
	qaws_surface_destroy(sa); qaws_surface_destroy(tube);
	for (i = 0; i < 3; i++) qaws_curve_destroy(open3[i]);
	qaws_curve_destroy(prof); qaws_curve_destroy(loop);
}

int test_108_shell_trimmed_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 108: shells of trimmed faces\n");
	test_split_plane();
	test_tube();
	test_trimmed_cube();
	test_tee(0);
	test_tee(1);
	test_cylinder();
	test_capped_tube();
	test_torus();
	test_shortest_path();
	test_shortest_holes();
	test_shortest_other_face();
	test_thin_holes();
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
