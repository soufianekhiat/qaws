/*
 * Test 111: sampling, meshing and offsets across the faces of a shell
 * (qaws_shell_region_sample, qaws_shell_region_triangulate, qaws_shell_region_offset)
 *
 *   - a geodesic disk round a cube's corner (three faces): its mesh's area
 *     against the region's; no cracks where the faces meet (the edges used
 *     by one triangle add up to the region's perimeter) and no edge in more
 *     than two triangles
 *   - samples of the disk: all inside, shared out by the faces' areas
 *   - offsets: disks across edges and round a corner grown by d against the
 *     disks of r + d, and shrunk (no vertex inside) against r - d
 */

#include "test_common.h"
#include "qaws_shell_region.h"
#include "qaws_surface_bilinear.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

static qaws_surface* smt_plane(double const* p00, double const* p10, double const* p01, double const* p11)
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

typedef struct smt_edge
{
	unsigned int a, b;
} smt_edge;

static int smt_cmp(void const* x, void const* y)
{
	smt_edge const *p = (smt_edge const*)x, *q = (smt_edge const*)y;
	if (p->a != q->a) return p->a < q->a ? -1 : 1;
	return p->b < q->b ? -1 : p->b > q->b;
}

static double smt_d3(qaws_vec3 a, qaws_vec3 b)
{
	return sqrt(((double)a.x - b.x) * ((double)a.x - b.x) + ((double)a.y - b.y) * ((double)a.y - b.y) + ((double)a.z - b.z) * ((double)a.z - b.z));
}

static void test_cube_corner(void)
{
	static double const P[8][3] = { { -1, -1, -1 }, { 1, -1, -1 }, { -1, 1, -1 }, { 1, 1, -1 }, { -1, -1, 1 }, { 1, -1, 1 }, { -1, 1, 1 }, { 1, 1, 1 } };
	static int const C[6][4] = { { 4, 5, 6, 7 }, { 0, 2, 1, 3 }, { 1, 3, 5, 7 }, { 0, 4, 2, 6 }, { 2, 6, 3, 7 }, { 0, 1, 4, 5 } };
	qaws_surface* f[6];
	qaws_shell_atlas* a = NULL;
	qaws_shell_region* g = NULL;
	qaws_shell_mesh* m = NULL;
	qaws_shell_point c;
	qaws_scalar A = 0, Pr = 0;
	qaws_status st;
	unsigned int i, nt = 0, ne = 0, over = 0, faces = 0;
	double area = 0, blen = 0;
	smt_edge* E = NULL;
	char msg[300];
	for (i = 0; i < 6; i++) f[i] = smt_plane(P[C[i][0]], P[C[i][1]], P[C[i][2]], P[C[i][3]]);
	qaws_shell_atlas_create((qaws_surface const* const*)f, 6, 0, &a);
	c.face = 0; c.uv.x = (qaws_scalar)0.85; c.uv.y = (qaws_scalar)0.85;
	st = qaws_shell_region_create_geodesic_disk(a, c, (qaws_scalar)0.75, 0, &g);
	if (st == QAWS_STATUS_OK)
	{
		qaws_shell_region_get_area(g, &A);
		qaws_shell_region_get_perimeter(g, &Pr);
		for (i = 0; i < 6; i++) faces += qaws_shell_region_get_face(g, i) != NULL;
		st = qaws_shell_region_triangulate(g, NULL, &m);
	}
	if (st == QAWS_STATUS_OK)
	{
		nt = qaws_shell_mesh_get_triangle_count(m);
		E = (smt_edge*)malloc(sizeof(smt_edge) * 3 * nt + 1);
		for (i = 0; E && i < nt; i++)
		{
			unsigned int t[3], k;
			qaws_vec3 x[3];
			double u[3], v[3], cx, cy, cz;
			qaws_shell_mesh_get_triangle(m, i, t);
			for (k = 0; k < 3; k++) qaws_shell_mesh_get_vertex(m, t[k], NULL, &x[k]);
			u[0] = x[1].x - x[0].x; u[1] = x[1].y - x[0].y; u[2] = x[1].z - x[0].z;
			v[0] = x[2].x - x[0].x; v[1] = x[2].y - x[0].y; v[2] = x[2].z - x[0].z;
			cx = u[1] * v[2] - u[2] * v[1]; cy = u[2] * v[0] - u[0] * v[2]; cz = u[0] * v[1] - u[1] * v[0];
			area += 0.5 * sqrt(cx * cx + cy * cy + cz * cz);
			for (k = 0; k < 3; k++)
			{
				unsigned int p = t[k], q = t[(k + 1) % 3];
				E[ne].a = p < q ? p : q; E[ne].b = p < q ? q : p;
				ne++;
			}
		}
		if (E)
		{
			qsort(E, ne, sizeof(smt_edge), smt_cmp);
			for (i = 0; i < ne; )
			{
				unsigned int j = i;
				while (j < ne && E[j].a == E[i].a && E[j].b == E[i].b) j++;
				if (j - i == 1)
				{
					qaws_vec3 x, y;
					qaws_shell_mesh_get_vertex(m, E[i].a, NULL, &x);
					qaws_shell_mesh_get_vertex(m, E[i].b, NULL, &y);
					blen += smt_d3(x, y);
				}
				if (j - i > 2)
					over++;
				i = j;
			}
		}
	}
	sprintf(msg, "shell mesh of a disk round a cube's corner (%u faces, %u triangles): area %.2g off, boundary %.2g off the perimeter, %u edges in > 2 triangles (status %d)",
		faces, nt, fabs(area - (double)A), fabs(blen - (double)Pr), over, (int)st);
	TEST_ASSERT(st == QAWS_STATUS_OK && faces == 3 && fabs(area - (double)A) <= 2e-3 * (double)A && fabs(blen - (double)Pr) <= 1e-3 * (double)Pr && over == 0, msg);
	{
		enum { N = 1500 };
		static qaws_shell_point pts[N];
		unsigned int per[6] = { 0, 0, 0, 0, 0, 0 }, out = 0, off = 0;
		st = g ? qaws_shell_region_sample(g, NULL, N, pts, NULL) : QAWS_STATUS_INVALID_ARGUMENT;
		for (i = 0; i < N && st == QAWS_STATUS_OK; i++)
		{
			int in = 0;
			per[pts[i].face]++;
			qaws_shell_region_contains(g, pts[i], &in);
			out += in == 0;
		}
		for (i = 0; i < 6 && st == QAWS_STATUS_OK; i++)
		{
			qaws_scalar fa = 0;
			if (qaws_shell_region_get_face(g, i)) qaws_surface_region_get_area(qaws_shell_region_get_face(g, i), &fa);
			if (fabs(per[i] - N * (double)fa / (double)A) > 1) off++;
		}
		sprintf(msg, "shell samples of the disk: %u outside, %u faces off their share by area (status %d)", out, off, (int)st);
		TEST_ASSERT(st == QAWS_STATUS_OK && out == 0 && off == 0, msg);
	}
	free(E);
	qaws_shell_mesh_destroy(m);
	qaws_shell_region_destroy(g);
	qaws_shell_atlas_destroy(a);
	for (i = 0; i < 6; i++) qaws_surface_destroy(f[i]);
}

static void test_offsets(void)
{
	/* a geodesic disk's offset by d is the disk of radius r + d (the points
	   within d of it), and shrunk the disk of radius r - d where no vertex
	   is inside (a cone point leaves points behind it farther from the
	   boundary): two disks across edges and one round a cube's corner */
	static double const P[8][3] = { { -1, -1, -1 }, { 1, -1, -1 }, { -1, 1, -1 }, { 1, 1, -1 }, { -1, -1, 1 }, { 1, -1, 1 }, { -1, 1, 1 }, { 1, 1, 1 } };
	static int const C[6][4] = { { 4, 5, 6, 7 }, { 0, 2, 1, 3 }, { 1, 3, 5, 7 }, { 0, 4, 2, 6 }, { 2, 6, 3, 7 }, { 0, 1, 4, 5 } };
	static double const cases[3][4] = { { 0.5, 0.5, 0.4, 1 }, { 0.7, 0.4, 0.5, 1 }, { 0.85, 0.85, 0.75, 0 } };
	qaws_surface* f[6];
	qaws_shell_atlas* a = NULL;
	unsigned int i, z, bad = 0;
	double worst = 0, tol = sizeof(qaws_scalar) == 4 ? 5e-5 : 1e-8;
	char msg[300];
	for (i = 0; i < 6; i++) f[i] = smt_plane(P[C[i][0]], P[C[i][1]], P[C[i][2]], P[C[i][3]]);
	qaws_shell_atlas_create((qaws_surface const* const*)f, 6, 0, &a);
	for (z = 0; z < 3; z++)
	{
		qaws_shell_point c;
		qaws_shell_region *g = NULL, *gr = NULL, *sh = NULL, *big = NULL, *small = NULL;
		qaws_scalar A1 = 0, A2 = 0, B1 = 0, B2 = 0;
		double r = cases[z][2];
		c.face = 0; c.uv.x = (qaws_scalar)cases[z][0]; c.uv.y = (qaws_scalar)cases[z][1];
		qaws_shell_region_create_geodesic_disk(a, c, (qaws_scalar)r, 0, &g);
		qaws_shell_region_create_geodesic_disk(a, c, (qaws_scalar)(r + 0.25), 0, &big);
		if (cases[z][3]) qaws_shell_region_create_geodesic_disk(a, c, (qaws_scalar)(r - 0.25), 0, &small);
		if (!g || !big || qaws_shell_region_offset(g, (qaws_scalar)0.25, 0, &gr) != QAWS_STATUS_OK)
			bad++;
		if (cases[z][3] && (!small || qaws_shell_region_offset(g, (qaws_scalar)-0.25, 0, &sh) != QAWS_STATUS_OK))
			bad++;
		if (gr) { qaws_shell_region_get_area(gr, &A1); qaws_shell_region_get_area(big, &B1); worst = fmax(worst, fabs((double)A1 - B1)); }
		if (sh) { qaws_shell_region_get_area(sh, &A2); qaws_shell_region_get_area(small, &B2); worst = fmax(worst, fabs((double)A2 - B2)); }
		qaws_shell_region_destroy(g); qaws_shell_region_destroy(gr); qaws_shell_region_destroy(sh);
		qaws_shell_region_destroy(big); qaws_shell_region_destroy(small);
	}
	sprintf(msg, "shell offsets: disks grown (two across edges, one round a corner) and shrunk against the disks of r +- d: %u failed, %.2g off", bad, worst);
	TEST_ASSERT(bad == 0 && worst <= tol, msg);
	qaws_shell_atlas_destroy(a);
	for (i = 0; i < 6; i++) qaws_surface_destroy(f[i]);
}

int test_111_shell_mesh_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 111: sampling and meshing on shells\n");
	test_cube_corner();
	test_offsets();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
