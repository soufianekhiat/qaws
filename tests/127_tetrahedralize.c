/*
 * Test 127: volume meshes and tetrahedralization (DELAUNAY)
 *
 *   - the container: a unit cube as one hexahedron (scaled Jacobian 1,
 *     volume 1, its border six quads) and as six tetrahedra round its
 *     diagonal (facets matched pairwise, the border twelve triangles, closed);
 *     an inverted cell refused by validation; a facet of three cells refused
 *     at creation;
 *   - a cube's surface (six quads): its corners' tetrahedra (the border's
 *     coplanar ties flipped), then filled; the volume exactly 1;
 *   - a NURBS torus's generated surface filled: valid, its volume exactly
 *     the surface's, its border the surface, the angles bounded;
 *   - a box generated as a shell of six faces, filled: the volume exact;
 *   - refusals: an open surface, the reserved methods.
 */

#include "test_common.h"
#include "qaws_tetrahedralize.h"
#include "qaws_mesh_generate.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

static qaws_scalar const tv_cube_x[24] = { 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1 };

static qaws_mesh_3d* tv_cube_surface(void)
{
	static unsigned int const f[24] = { 0, 3, 2, 1, 4, 5, 6, 7, 0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6, 3, 0, 4, 7 };
	static unsigned int const fs[6] = { 4, 4, 4, 4, 4, 4 };
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	memset(&d, 0, sizeof d);
	d.positions = tv_cube_x;
	d.vertex_count = 8;
	d.face_vertices = f;
	d.face_sizes = fs;
	d.face_count = 6;
	qaws_mesh_3d_create(&d, &m);
	return m;
}

/* the volume a closed surface encloses */
static double tv_enclosed(qaws_mesh_3d const* m)
{
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), nf = qaws_mesh_3d_get_face_count(m), f, q, k;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	unsigned int* fv = (unsigned int*)malloc(sizeof(unsigned int) * qaws_mesh_3d_get_corner_count(m));
	unsigned int* fs = (unsigned int*)malloc(sizeof(unsigned int) * nf);
	double vol = 0;
	qaws_mesh_3d_get_positions(m, x);
	qaws_mesh_3d_get_faces(m, fv, fs);
	for (f = 0, q = 0; f < nf; q += fs[f], f++)
		for (k = 1; k + 1 < fs[f]; k++)
		{
			qaws_scalar const *a = x + 3 * fv[q], *b = x + 3 * fv[q + k], *c = x + 3 * fv[q + k + 1];
			vol += ((double)a[0] * ((double)b[1] * c[2] - (double)b[2] * c[1]) + (double)a[1] * ((double)b[2] * c[0] - (double)b[0] * c[2])
				+ (double)a[2] * ((double)b[0] * c[1] - (double)b[1] * c[0])) / 6;
		}
	free(x); free(fv); free(fs);
	return vol;
}

static void test_container(void)
{
	static unsigned int const hex[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
	/* the six tetrahedra round the diagonal 0-6, one per order of the axes */
	static unsigned int const path[6][2] = { { 1, 2 }, { 1, 5 }, { 3, 2 }, { 3, 7 }, { 4, 5 }, { 4, 7 } };
	qaws_volume_mesh_desc d;
	qaws_volume_mesh* m = NULL;
	qaws_volume_quality q;
	qaws_mesh_3d* b = NULL;
	qaws_mesh_topology top;
	qaws_cell_type ht = QAWS_CELL_HEXAHEDRON;
	unsigned int tets[24], adj[24], i, nborder = 0, sym = 1;
	qaws_status s1, s2;
	char msg[300];
	memset(&d, 0, sizeof d);
	d.positions = tv_cube_x;
	d.vertex_count = 8;
	d.cell_vertices = hex;
	d.cell_types = &ht;
	d.cell_count = 1;
	s1 = qaws_volume_mesh_create(&d, &m);
	if (s1 == QAWS_STATUS_OK)
	{
		qaws_volume_mesh_measure_quality(m, &q);
		s2 = qaws_volume_mesh_get_boundary(m, &b);
		qaws_mesh_3d_get_topology(b, &top);
		sprintf(msg, "a hexahedron: validate %s, scaled Jacobian %.3f, volume %.6f, border %u faces (%u quads), closed %d", qaws_status_to_string(qaws_volume_mesh_validate(m)),
			(double)q.min_scaled_jacobian, (double)q.volume, qaws_mesh_3d_get_face_count(b), qaws_mesh_3d_get_face_count(b), top.closed);
		TEST_ASSERT(s2 == QAWS_STATUS_OK && qaws_volume_mesh_validate(m) == QAWS_STATUS_OK && q.hexahedra == 1 && fabs((double)q.min_scaled_jacobian - 1) < 1e-12
			&& fabs((double)q.volume - 1) < 1e-12 && qaws_mesh_3d_get_face_count(b) == 6 && top.closed && top.genus == 0, msg);
		qaws_mesh_3d_destroy(b);
		b = NULL;
	}
	else
		TEST_ASSERT(0, "a hexahedron: created");
	qaws_volume_mesh_destroy(m);
	m = NULL;
	for (i = 0; i < 6; i++)
	{
		qaws_scalar const *p0 = tv_cube_x, *p1 = tv_cube_x + 3 * path[i][0], *p2 = tv_cube_x + 3 * path[i][1], *p3 = tv_cube_x + 18;
		double u[3] = { p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2] }, v[3] = { p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2] };
		double w[3] = { p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2] };
		double det = u[0] * (v[1] * w[2] - v[2] * w[1]) + u[1] * (v[2] * w[0] - v[0] * w[2]) + u[2] * (v[0] * w[1] - v[1] * w[0]);
		tets[4 * i] = 0;
		tets[4 * i + 1] = det > 0 ? path[i][0] : path[i][1];
		tets[4 * i + 2] = det > 0 ? path[i][1] : path[i][0];
		tets[4 * i + 3] = 6;
	}
	d.cell_vertices = tets;
	d.cell_types = NULL;
	d.cell_count = 6;
	s1 = qaws_volume_mesh_create(&d, &m);
	if (s1 == QAWS_STATUS_OK)
	{
		qaws_volume_mesh_get_cell_adjacency(m, adj);
		for (i = 0; i < 24; i++)
		{
			unsigned int o = adj[i], k, back = 0;
			if (o == QAWS_MESH_NONE)
			{
				nborder++;
				continue;
			}
			for (k = 0; k < 4; k++)
				back |= adj[4 * o + k] == i / 4;
			sym &= back;
		}
		qaws_volume_mesh_measure_quality(m, &q);
		qaws_volume_mesh_get_boundary(m, &b);
		qaws_mesh_3d_get_topology(b, &top);
		sprintf(msg, "six tetrahedra: validate %s, %u border facets, adjacency symmetric %u, volume %.6f, border %u faces closed %d",
			qaws_status_to_string(qaws_volume_mesh_validate(m)), nborder, sym, (double)q.volume, qaws_mesh_3d_get_face_count(b), top.closed);
		TEST_ASSERT(qaws_volume_mesh_validate(m) == QAWS_STATUS_OK && qaws_volume_mesh_get_facet_count(m) == 24 && nborder == 12 && sym
			&& fabs((double)q.volume - 1) < 1e-12 && qaws_mesh_3d_get_face_count(b) == 12 && top.closed && top.genus == 0, msg);
		qaws_mesh_3d_destroy(b);
	}
	else
		TEST_ASSERT(0, "six tetrahedra: created");
	qaws_volume_mesh_destroy(m);
	m = NULL;
	/* one cell turned over; one facet in three cells */
	{
		unsigned int bad[4] = { 0, 2, 1, 4 }, three[12] = { 0, 1, 2, 4, 0, 2, 1, 5, 0, 1, 2, 6 };
		qaws_status s3, s4;
		d.cell_vertices = bad;
		d.cell_count = 1;
		s3 = qaws_volume_mesh_create(&d, &m);
		s4 = m ? qaws_volume_mesh_validate(m) : QAWS_STATUS_OK;
		qaws_volume_mesh_destroy(m);
		m = NULL;
		d.cell_vertices = three;
		d.cell_count = 3;
		s1 = qaws_volume_mesh_create(&d, &m);
		sprintf(msg, "refusals: an inverted cell validates %s, a facet of three cells %s", qaws_status_to_string(s4), qaws_status_to_string(s1));
		TEST_ASSERT(s3 == QAWS_STATUS_OK && s4 == QAWS_STATUS_INVALID_TOPOLOGY && s1 == QAWS_STATUS_INVALID_TOPOLOGY && !m, msg);
	}
}

/* a tetrahedralization's checks: valid, its volume the surface's, its border closed with the given genus and vertex count */
static void tv_check(char const* name, qaws_mesh_3d const* surface, qaws_volume_mesh const* v, qaws_status st, unsigned int genus, double min_dihedral)
{
	qaws_volume_quality q;
	qaws_mesh_3d* b = NULL;
	qaws_mesh_topology top;
	double encl = tv_enclosed(surface);
	char msg[400];
	if (st != QAWS_STATUS_OK || !v)
	{
		sprintf(msg, "%s: tetrahedralized (%s)", name, qaws_status_to_string(st));
		TEST_ASSERT(0, msg);
		return;
	}
	qaws_volume_mesh_measure_quality(v, &q);
	qaws_volume_mesh_get_boundary(v, &b);
	qaws_mesh_3d_get_topology(b, &top);
	sprintf(msg, "%s: %u tetrahedra, validate %s, volume %.12f for %.12f, border %u vertices (surface %u) closed %d genus %u, dihedral %.1f..%.1f",
		name, q.tetrahedra, qaws_status_to_string(qaws_volume_mesh_validate(v)), (double)q.volume, encl, qaws_mesh_3d_get_vertex_count(b),
		qaws_mesh_3d_get_vertex_count(surface), top.closed, top.genus, (double)q.min_dihedral, (double)q.max_dihedral);
	TEST_ASSERT(qaws_volume_mesh_validate(v) == QAWS_STATUS_OK && q.tetrahedra == qaws_volume_mesh_get_cell_count(v)
		&& fabs((double)q.volume - encl) < 1e-10 * fabs(encl) && qaws_mesh_3d_get_vertex_count(b) >= qaws_mesh_3d_get_vertex_count(surface) && top.closed
		&& top.genus == genus && (double)q.min_dihedral >= min_dihedral, msg);
	qaws_mesh_3d_destroy(b);
}

static void test_cube(void)
{
	qaws_mesh_3d* s = tv_cube_surface();
	qaws_tetrahedralize_desc d;
	qaws_volume_mesh* v = NULL;
	qaws_mesh_report rep;
	qaws_status st;
	memset(&d, 0, sizeof d);
	d.border_only = 1;
	st = qaws_volume_mesh_tetrahedralize(s, &d, &v, &rep);
	tv_check("a cube from its corners", s, v, st, 0, 30);
	TEST_ASSERT(!v || (qaws_volume_mesh_get_vertex_count(v) == 8 && rep.splits == 0), "a cube from its corners: no vertex added");
	qaws_volume_mesh_destroy(v);
	v = NULL;
	d.border_only = 0;
	d.target_length = (qaws_scalar)0.25;
	st = qaws_volume_mesh_tetrahedralize(s, &d, &v, &rep);
	tv_check("a cube filled", s, v, st, 0, 5);
	TEST_ASSERT(!v || (qaws_volume_mesh_get_vertex_count(v) > 30 && rep.moves > 20), "a cube filled: points inside");
	qaws_volume_mesh_destroy(v);
	qaws_mesh_3d_destroy(s);
}

static qaws_surface* tv_torus(double R, double r)
{
	static double const cx[9] = { 1, 1, 0, -1, -1, -1, 0, 1, 1 }, cy[9] = { 0, 1, 1, 1, 0, -1, -1, -1, 0 };
	static qaws_scalar const kn[12] = { 0, 0, 0, 0.25f, 0.25f, 0.5f, 0.5f, 0.75f, 0.75f, 1, 1, 1 };
	qaws_vec3 cp[81];
	qaws_scalar w[81];
	qaws_surface_nurbs_desc d;
	qaws_surface* s = NULL;
	unsigned int i, j;
	for (i = 0; i < 9; i++)
		for (j = 0; j < 9; j++)
		{
			double rad = R + r * cx[j];
			cp[i * 9 + j].x = (qaws_scalar)(rad * cx[i]);
			cp[i * 9 + j].y = (qaws_scalar)(rad * cy[i]);
			cp[i * 9 + j].z = (qaws_scalar)(r * cy[j]);
			w[i * 9 + j] = (qaws_scalar)((i & 1 ? sqrt(0.5) : 1) * (j & 1 ? sqrt(0.5) : 1));
		}
	memset(&d, 0, sizeof d);
	d.u_degree = 2;
	d.v_degree = 2;
	d.control_points = cp;
	d.u_point_count = 9;
	d.v_point_count = 9;
	d.weights = w;
	d.u_knots = kn;
	d.u_knot_count = 12;
	d.v_knots = kn;
	d.v_knot_count = 12;
	qaws_surface_create_nurbs(&d, &s);
	return s;
}

static void test_torus(void)
{
	qaws_surface* t = tv_torus(2, 0.7);
	qaws_generate_desc g;
	qaws_tetrahedralize_desc d;
	qaws_mesh_3d* s = NULL;
	qaws_volume_mesh* v = NULL;
	qaws_status st;
	memset(&g, 0, sizeof g);
	g.target_length = (qaws_scalar)0.25;
	st = qaws_mesh_generate_surface(t, &g, &s, NULL);
	memset(&d, 0, sizeof d);
	if (st == QAWS_STATUS_OK)
		st = qaws_volume_mesh_tetrahedralize(s, &d, &v, NULL);
	tv_check("a torus filled", s, v, st, 1, 8);
	if (v)
	{
		qaws_volume_quality q;
		qaws_volume_mesh_measure_quality(v, &q);
		TEST_ASSERT((double)q.max_radius_edge < 2.5 && (double)q.max_dihedral < 170, "a torus filled: circumradius over shortest edge below 2.5, dihedral below 170");
	}
	qaws_volume_mesh_destroy(v);
	qaws_mesh_3d_destroy(s);
	qaws_surface_destroy(t);
}

static void test_box(void)
{
	static int const q[6][4][3] = {
		{ { -1, -1, -1 }, { -1, 1, -1 }, { 1, -1, -1 }, { 1, 1, -1 } }, { { -1, -1, 1 }, { 1, -1, 1 }, { -1, 1, 1 }, { 1, 1, 1 } },
		{ { -1, -1, -1 }, { 1, -1, -1 }, { -1, -1, 1 }, { 1, -1, 1 } }, { { -1, 1, -1 }, { -1, 1, 1 }, { 1, 1, -1 }, { 1, 1, 1 } },
		{ { -1, -1, -1 }, { -1, -1, 1 }, { -1, 1, -1 }, { -1, 1, 1 } }, { { 1, -1, -1 }, { 1, 1, -1 }, { 1, -1, 1 }, { 1, 1, 1 } } };
	static double const h[3] = { 1, 0.7, 0.5 };
	qaws_surface* f[6];
	qaws_shell_atlas* atlas = NULL;
	qaws_generate_desc g;
	qaws_tetrahedralize_desc d;
	qaws_mesh_3d* s = NULL;
	qaws_volume_mesh* v = NULL;
	qaws_status st;
	unsigned int i, k;
	for (i = 0; i < 6; i++)
	{
		qaws_surface_bilinear_desc bd;
		qaws_vec3* p[4];
		p[0] = &bd.p00; p[1] = &bd.p10; p[2] = &bd.p01; p[3] = &bd.p11;
		for (k = 0; k < 4; k++)
		{
			p[k]->x = (qaws_scalar)(h[0] * q[i][k][0]);
			p[k]->y = (qaws_scalar)(h[1] * q[i][k][1]);
			p[k]->z = (qaws_scalar)(h[2] * q[i][k][2]);
		}
		f[i] = NULL;
		qaws_surface_create_bilinear(&bd, &f[i]);
	}
	st = qaws_shell_atlas_create((qaws_surface const* const*)f, 6, 0, &atlas);
	memset(&g, 0, sizeof g);
	g.target_length = (qaws_scalar)0.15;
	if (st == QAWS_STATUS_OK)
		st = qaws_mesh_generate_shell(atlas, &g, &s, NULL);
	memset(&d, 0, sizeof d);
	if (st == QAWS_STATUS_OK)
		st = qaws_volume_mesh_tetrahedralize(s, &d, &v, NULL);
	tv_check("a box filled", s, v, st, 0, 5);
	TEST_ASSERT(!v || fabs(tv_enclosed(s) - 2.8) < 1e-12, "a box filled: its surface encloses 2 x 1.4 x 1");
	qaws_volume_mesh_destroy(v);
	qaws_mesh_3d_destroy(s);
	qaws_shell_atlas_destroy(atlas);
	for (i = 0; i < 6; i++)
		qaws_surface_destroy(f[i]);
}

static void test_refusals(void)
{
	static qaws_scalar const p[12] = { 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0 };
	static unsigned int const f[4] = { 0, 1, 2, 3 };
	qaws_mesh_3d_desc md;
	qaws_mesh_3d *open = NULL, *cube = tv_cube_surface();
	qaws_tetrahedralize_desc d;
	qaws_volume_mesh* v = NULL;
	qaws_status s1, s2, s3;
	char msg[200];
	memset(&md, 0, sizeof md);
	md.positions = p;
	md.vertex_count = 4;
	md.face_vertices = f;
	md.face_count = 1;
	{
		unsigned int fs = 4;
		md.face_sizes = &fs;
		qaws_mesh_3d_create(&md, &open);
	}
	memset(&d, 0, sizeof d);
	s1 = qaws_volume_mesh_tetrahedralize(open, &d, &v, NULL);
	d.method = QAWS_TETS_ENVELOPE;
	s2 = qaws_volume_mesh_tetrahedralize(cube, &d, &v, NULL);
	s3 = qaws_volume_mesh_tetrahedralize(cube, NULL, &v, NULL);
	sprintf(msg, "refusals: an open surface %s, ENVELOPE %s, no description %s", qaws_status_to_string(s1), qaws_status_to_string(s2), qaws_status_to_string(s3));
	TEST_ASSERT(s1 == QAWS_STATUS_INVALID_TOPOLOGY && s2 == QAWS_STATUS_UNSUPPORTED_OPERATION && s3 == QAWS_STATUS_INVALID_ARGUMENT && !v, msg);
	qaws_mesh_3d_destroy(open);
	qaws_mesh_3d_destroy(cube);
}

int test_127_tetrahedralize_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 127: volume meshes and tetrahedralization\n");
	test_container();
	test_cube();
	test_torus();
	test_box();
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
