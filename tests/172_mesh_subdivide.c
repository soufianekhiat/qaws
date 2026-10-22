/*
 * Test 172: subdivision of meshes (qaws_mesh_3d_subdivide), small meshes
 *
 *   - Catmull-Clark on the cube [-1, 1]^3 of quads: 26 vertices, 24 quads,
 *     a corner at 5/9 (1, 1, 1), an edge's point at (3/4, 3/4, 0); a
 *     tetrahedron: 14 vertices, 12 quads, a sphere's topology;
 *   - Loop on the octahedron: a vertex of valence 4 at 1 - 4 beta = 0.515625,
 *     an edge's point at (3/8, 3/8, 0); with every edge a crease the
 *     octahedron itself;
 *   - the cube's edges creases: two levels keep the cube (volume 8, every
 *     vertex on it), 48 crease edges, the corners kept; semi-sharp, the
 *     volume grows with the sharpness;
 *   - a flat sheet of 2 x 2 quads: its corners stay, its border stays on it;
 *   - attributes: a vertex value linear in the position follows it exactly,
 *     a normalized one stays unit, corners' values linear in each face,
 *     faces' ids their parent's, a nearest vertex label kept;
 *   - refusals: Loop on quads, a negative sharpness, an unknown scheme.
 */

#include "test_common.h"
#include "qaws_mesh_subdivide.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

static unsigned int const sd_cube_faces[24] = { 0, 2, 3, 1, 4, 5, 7, 6, 0, 1, 5, 4, 2, 6, 7, 3, 0, 4, 6, 2, 1, 3, 7, 5 };

static qaws_mesh_3d* sd_make(qaws_scalar const* x, unsigned int nv, unsigned int const* fv, unsigned int const* fs, unsigned int nf)
{
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	memset(&d, 0, sizeof d);
	d.positions = x;
	d.vertex_count = nv;
	d.face_vertices = fv;
	d.face_sizes = fs;
	d.face_count = nf;
	qaws_mesh_3d_create(&m, &d);
	return m;
}

static qaws_mesh_3d* sd_cube(void)
{
	qaws_scalar x[24];
	unsigned int fs[6] = { 4, 4, 4, 4, 4, 4 }, i;
	for (i = 0; i < 8; i++)
	{
		x[3 * i] = (qaws_scalar)(i & 1 ? 1 : -1);
		x[3 * i + 1] = (qaws_scalar)(i & 2 ? 1 : -1);
		x[3 * i + 2] = (qaws_scalar)(i & 4 ? 1 : -1);
	}
	return sd_make(x, 8, sd_cube_faces, fs, 6);
}

static qaws_mesh_3d* sd_octahedron(void)
{
	static qaws_scalar const x[18] = { 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1 };
	static unsigned int const t[24] = { 0, 2, 4, 2, 1, 4, 1, 3, 4, 3, 0, 4, 2, 0, 5, 1, 2, 5, 3, 1, 5, 0, 3, 5 };
	return sd_make(x, 6, t, NULL, 8);
}

static qaws_mesh_3d* sd_run(qaws_mesh_3d const* m, qaws_mesh_subdivide_scheme scheme, unsigned int levels, qaws_scalar const* sharpness)
{
	qaws_mesh_subdivide_desc d;
	qaws_mesh_3d* o = NULL;
	memset(&d, 0, sizeof d);
	d.scheme = scheme;
	d.levels = levels;
	d.sharpness = sharpness;
	qaws_mesh_3d_subdivide(&o, m, &d);
	return o;
}

/* the enclosed volume, over the faces' fan triangles */
static double sd_volume(qaws_mesh_3d const* m)
{
	unsigned int nt = qaws_mesh_3d_get_triangle_count(m), i;
	double v = 0;
	for (i = 0; i < nt; i++)
	{
		unsigned int t[3];
		qaws_vec3 a, b, c;
		qaws_mesh_3d_get_triangle(t, m, i);
		a = qaws_mesh_3d_get_vertex(m, t[0]);
		b = qaws_mesh_3d_get_vertex(m, t[1]);
		c = qaws_mesh_3d_get_vertex(m, t[2]);
		v += ((double)a.x * ((double)b.y * c.z - (double)b.z * c.y) - (double)a.y * ((double)b.x * c.z - (double)b.z * c.x) +
			(double)a.z * ((double)b.x * c.y - (double)b.y * c.x)) / 6;
	}
	return v;
}

/* the distance from p to the nearest vertex */
static double sd_nearest(qaws_mesh_3d const* m, double x, double y, double z)
{
	unsigned int n = qaws_mesh_3d_get_vertex_count(m), i;
	double best = 1e30;
	for (i = 0; i < n; i++)
	{
		qaws_vec3 p = qaws_mesh_3d_get_vertex(m, i);
		double d = sqrt(((double)p.x - x) * ((double)p.x - x) + ((double)p.y - y) * ((double)p.y - y) + ((double)p.z - z) * ((double)p.z - z));
		if (d < best)
			best = d;
	}
	return best;
}

static void test_counts(void)
{
	qaws_mesh_3d *c = sd_cube(), *o;
	qaws_mesh_quality q;
	qaws_mesh_topology tp;
	qaws_vec3 p;
	char msg[300];
	memset(&q, 0, sizeof q);
	memset(&tp, 0, sizeof tp);
	o = sd_run(c, QAWS_MESH_SUBDIVIDE_CATMULL_CLARK, 1, NULL);
	TEST_ASSERT(o != NULL, "Catmull-Clark on the cube");
	if (o)
	{
		qaws_mesh_3d_measure_quality(&q, o);
		qaws_mesh_3d_get_topology(&tp, o);
		p = qaws_mesh_3d_get_vertex(o, 7);
		sprintf(msg, "Catmull-Clark on the cube: %u vertices, %u quads of %u faces, closed %d; corner 7 at (%.12f %.12f %.12f), 5/9; an edge's point (3/4, 3/4, 0) %.1e away",
			qaws_mesh_3d_get_vertex_count(o), q.quads, qaws_mesh_3d_get_face_count(o), tp.closed, (double)p.x, (double)p.y, (double)p.z,
			sd_nearest(o, 0.75, 0.75, 0));
		printf("  %s\n", msg);
		TEST_ASSERT(qaws_mesh_3d_get_vertex_count(o) == 26 && q.quads == 24 && qaws_mesh_3d_get_face_count(o) == 24 && tp.closed &&
			fabs((double)p.x - 5.0 / 9) < 1e-12 && fabs((double)p.y - 5.0 / 9) < 1e-12 && fabs((double)p.z - 5.0 / 9) < 1e-12 &&
			sd_nearest(o, 0.75, 0.75, 0) < 1e-12 && qaws_mesh_3d_validate(o) == QAWS_STATUS_OK, msg);
		qaws_mesh_3d_destroy(o);
	}
	qaws_mesh_3d_destroy(c);
	{
		static qaws_scalar const x[12] = { 1, 1, 1, 1, -1, -1, -1, 1, -1, -1, -1, 1 };
		static unsigned int const t[12] = { 0, 1, 2, 0, 3, 1, 0, 2, 3, 1, 3, 2 };
		qaws_mesh_3d* tet = sd_make(x, 4, t, NULL, 4);
		o = sd_run(tet, QAWS_MESH_SUBDIVIDE_CATMULL_CLARK, 1, NULL);
		if (o)
		{
			qaws_mesh_3d_get_topology(&tp, o);
			qaws_mesh_3d_measure_quality(&q, o);
		}
		sprintf(msg, "a tetrahedron: %u vertices, %u quads, Euler %d, closed %d", o ? qaws_mesh_3d_get_vertex_count(o) : 0, o ? q.quads : 0,
			o ? tp.euler_characteristic : 0, o ? tp.closed : 0);
		printf("  %s\n", msg);
		TEST_ASSERT(o && qaws_mesh_3d_get_vertex_count(o) == 14 && q.quads == 12 && tp.euler_characteristic == 2 && tp.closed, msg);
		qaws_mesh_3d_destroy(o);
		qaws_mesh_3d_destroy(tet);
	}
}

static void test_loop(void)
{
	qaws_mesh_3d *oc = sd_octahedron(), *o, *s;
	qaws_vec3 p;
	double v = 0;
	char msg[300];
	o = sd_run(oc, QAWS_MESH_SUBDIVIDE_LOOP, 1, NULL);
	TEST_ASSERT(o != NULL, "Loop on the octahedron");
	if (!o)
		return;
	p = qaws_mesh_3d_get_vertex(o, 0);
	qaws_mesh_3d_detect_features(oc, (qaws_scalar)0.5, 0, 0);
	s = sd_run(oc, QAWS_MESH_SUBDIVIDE_LOOP, 2, NULL);
	if (s)
		v = sd_volume(s);
	sprintf(msg, "Loop on the octahedron: %u vertices, %u triangles; vertex 0 at %.12f (0.515625); an edge's point (3/8, 3/8, 0) %.1e away; creased, two levels: volume %.12f (4/3), %u faces",
		qaws_mesh_3d_get_vertex_count(o), qaws_mesh_3d_get_face_count(o), (double)p.x, sd_nearest(o, 0.375, 0.375, 0), v, s ? qaws_mesh_3d_get_face_count(s) : 0);
	printf("  %s\n", msg);
	TEST_ASSERT(qaws_mesh_3d_get_vertex_count(o) == 18 && qaws_mesh_3d_get_face_count(o) == 32 && fabs((double)p.x - 0.515625) < 1e-12 &&
		fabs((double)p.y) < 1e-12 && sd_nearest(o, 0.375, 0.375, 0) < 1e-12 && s && fabs(v - 4.0 / 3) < 1e-12 && qaws_mesh_3d_get_face_count(s) == 128, msg);
	qaws_mesh_3d_destroy(o);
	qaws_mesh_3d_destroy(s);
	qaws_mesh_3d_destroy(oc);
}

static void test_creases(void)
{
	qaws_mesh_3d *c = sd_cube(), *o;
	double v = 0, off = 0, vs[5];
	unsigned int i, corners = 0, n = 0;
	qaws_scalar sh[12];
	static double const level[5] = { 0, 0.5, 1, 2, -1 };
	char msg[400];
	qaws_mesh_3d_detect_features(c, (qaws_scalar)0.5, 0, 0);
	o = sd_run(c, QAWS_MESH_SUBDIVIDE_CATMULL_CLARK, 2, NULL);
	if (o)
	{
		v = sd_volume(o);
		n = qaws_mesh_3d_get_vertex_count(o);
		for (i = 0; i < n; i++)
		{
			qaws_vec3 p = qaws_mesh_3d_get_vertex(o, i);
			double m = fmax(fabs((double)p.x), fmax(fabs((double)p.y), fabs((double)p.z)));
			off = fmax(off, fabs(m - 1));
		}
		for (i = 0; i < 8; i++) corners += qaws_mesh_3d_is_feature_corner(o, i) != 0;
	}
	sprintf(msg, "the cube's edges creases, two levels: volume %.12f (8), vertices off the cube %.1e, %u crease edges (48), corners kept %u of 8",
		v, off, o ? qaws_mesh_3d_get_feature_edge_count(o) : 0, corners);
	printf("  %s\n", msg);
	TEST_ASSERT(o && fabs(v - 8) < 1e-12 && off < 1e-12 && qaws_mesh_3d_get_feature_edge_count(o) == 48 && corners == 8, msg);
	qaws_mesh_3d_destroy(o);
	for (i = 0; i < 5; i++)
	{
		unsigned int j;
		for (j = 0; j < 12; j++) sh[j] = (qaws_scalar)level[i];
		o = sd_run(c, QAWS_MESH_SUBDIVIDE_CATMULL_CLARK, 3, level[i] < 0 ? NULL : sh);
		vs[i] = o ? sd_volume(o) : 0;
		qaws_mesh_3d_destroy(o);
	}
	sprintf(msg, "semi-sharp, three levels: volumes %.5f (smooth), %.5f (0.5), %.5f (1), %.5f (2), %.5f (sharp)", vs[0], vs[1], vs[2], vs[3], vs[4]);
	printf("  %s\n", msg);
	TEST_ASSERT(vs[0] > 1 && vs[0] < vs[1] && vs[1] < vs[2] && vs[2] < vs[3] && vs[3] < vs[4] && fabs(vs[4] - 8) < 1e-12, msg);
	qaws_mesh_3d_destroy(c);
}

static void test_sheet(void)
{
	qaws_scalar x[27];
	unsigned int fv[16] = { 0, 1, 4, 3, 1, 2, 5, 4, 3, 4, 7, 6, 4, 5, 8, 7 }, fs[4] = { 4, 4, 4, 4 }, i, n = 0, off = 0;
	qaws_mesh_3d *s, *o;
	double fixed = 0, flat = 0;
	char msg[300];
	for (i = 0; i < 9; i++)
	{
		x[3 * i] = (qaws_scalar)(i % 3);
		x[3 * i + 1] = (qaws_scalar)(i / 3);
		x[3 * i + 2] = (qaws_scalar)(0.3 * (i == 4));
	}
	s = sd_make(x, 9, fv, fs, 4);
	o = sd_run(s, QAWS_MESH_SUBDIVIDE_CATMULL_CLARK, 2, NULL);
	if (o)
	{
		static unsigned int const cs[4] = { 0, 2, 6, 8 };
		n = qaws_mesh_3d_get_vertex_count(o);
		for (i = 0; i < 4; i++)
		{
			qaws_vec3 p = qaws_mesh_3d_get_vertex(o, cs[i]);
			fixed = fmax(fixed, fabs((double)p.x - (double)x[3 * cs[i]]) + fabs((double)p.y - (double)x[3 * cs[i] + 1]) + fabs((double)p.z));
		}
		for (i = 0; i < n; i++)
			if (qaws_mesh_3d_is_boundary_vertex(o, i))
			{
				qaws_vec3 p = qaws_mesh_3d_get_vertex(o, i);
				double dx = fmin(fabs((double)p.x), fabs((double)p.x - 2)), dy = fmin(fabs((double)p.y), fabs((double)p.y - 2));
				flat = fmax(flat, fmin(dx, dy) + fabs((double)p.z));
				off++;
			}
	}
	sprintf(msg, "a sheet of 2 x 2 quads, two levels: %u vertices, corners moved %.1e, %u border vertices off the border %.1e", n, fixed, off, flat);
	printf("  %s\n", msg);
	TEST_ASSERT(o && n == 81 && fixed < 1e-12 && off == 32 && flat < 1e-12, msg);
	qaws_mesh_3d_destroy(o);
	qaws_mesh_3d_destroy(s);
}

static void test_attributes(void)
{
	qaws_mesh_3d *c = sd_cube(), *o;
	qaws_scalar vx[8], vn[24], cuv[48], fid[6], lab[8];
	qaws_attribute_desc d;
	unsigned int i, ids[5], nv = 0, nf = 0, nc = 0, bad_id = 0, bad_lab = 0;
	double ex = 0, en = 0, euv = 0;
	char msg[400];
	for (i = 0; i < 8; i++)
	{
		qaws_vec3 p = qaws_mesh_3d_get_vertex(c, i);
		double l = sqrt(3.0);
		vx[i] = (qaws_scalar)(2 * (double)p.x + (double)p.y - 1);
		vn[3 * i] = (qaws_scalar)(p.x / l);
		vn[3 * i + 1] = (qaws_scalar)(p.y / l);
		vn[3 * i + 2] = (qaws_scalar)(p.z / l);
		lab[i] = (qaws_scalar)i;
	}
	for (i = 0; i < 24; i++)
	{
		cuv[2 * i] = (qaws_scalar)(i % 4 == 1 || i % 4 == 2);
		cuv[2 * i + 1] = (qaws_scalar)(i % 4 >= 2) + (qaws_scalar)(i / 4);
	}
	for (i = 0; i < 6; i++) fid[i] = (qaws_scalar)(10 + i);
	memset(&d, 0, sizeof d);
	d.name = "linear";
	d.element = QAWS_MESH_VERTEX;
	d.components = 1;
	d.rule = QAWS_ATTRIBUTE_LINEAR;
	qaws_mesh_3d_add_attribute(&ids[0], c, vx, &d);
	d.name = "normal";
	d.components = 3;
	d.rule = QAWS_ATTRIBUTE_NORMALIZED;
	qaws_mesh_3d_add_attribute(&ids[1], c, vn, &d);
	d.name = "uv";
	d.element = QAWS_MESH_CORNER;
	d.components = 2;
	d.rule = QAWS_ATTRIBUTE_LINEAR;
	qaws_mesh_3d_add_attribute(&ids[2], c, cuv, &d);
	d.name = "id";
	d.element = QAWS_MESH_FACE;
	d.components = 1;
	d.rule = QAWS_ATTRIBUTE_NEAREST;
	qaws_mesh_3d_add_attribute(&ids[3], c, fid, &d);
	d.name = "label";
	d.element = QAWS_MESH_VERTEX;
	qaws_mesh_3d_add_attribute(&ids[4], c, lab, &d);
	o = sd_run(c, QAWS_MESH_SUBDIVIDE_CATMULL_CLARK, 2, NULL);
	if (o)
	{
		qaws_scalar *a, *b, *u, *f, *l;
		unsigned int id[5];
		nv = qaws_mesh_3d_get_vertex_count(o);
		nf = qaws_mesh_3d_get_face_count(o);
		nc = qaws_mesh_3d_get_corner_count(o);
		for (i = 0; i < 5; i++)
		{
			static char const* const names[5] = { "linear", "normal", "uv", "id", "label" };
			qaws_mesh_3d_find_attribute(&id[i], o, names[i]);
		}
		a = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv);
		b = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
		u = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * nc);
		f = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nf);
		l = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv);
		qaws_mesh_3d_get_attribute(a, o, id[0]);
		qaws_mesh_3d_get_attribute(b, o, id[1]);
		qaws_mesh_3d_get_attribute(u, o, id[2]);
		qaws_mesh_3d_get_attribute(f, o, id[3]);
		qaws_mesh_3d_get_attribute(l, o, id[4]);
		for (i = 0; i < nv; i++)
		{
			qaws_vec3 p = qaws_mesh_3d_get_vertex(o, i);
			ex = fmax(ex, fabs((double)a[i] - (2 * (double)p.x + (double)p.y - 1)));
			en = fmax(en, fabs(sqrt((double)b[3 * i] * b[3 * i] + (double)b[3 * i + 1] * b[3 * i + 1] + (double)b[3 * i + 2] * b[3 * i + 2]) - 1));
			bad_lab += i < 8 && (unsigned int)l[i] != i;
		}
		for (i = 0; i < nf; i++) bad_id += (unsigned int)f[i] != 10 + i / 16;
		/* face 0's first grandchild: corner 0's uv (0, 0), then its quarter edge, centre, quarter edge of a unit square */
		euv = fabs((double)u[0]) + fabs((double)u[1]) + fabs((double)u[2] - 0.25) + fabs((double)u[3]) + fabs((double)u[4] - 0.25) + fabs((double)u[5] - 0.25) +
			fabs((double)u[6]) + fabs((double)u[7] - 0.25);
		free(a);
		free(b);
		free(u);
		free(f);
		free(l);
	}
	sprintf(msg, "attributes, two levels (%u vertices, %u faces): linear off %.1e, normalized off unit %.1e, corners' uv off %.1e, faces' ids wrong %u, labels wrong %u",
		nv, nf, ex, en, euv, bad_id, bad_lab);
	printf("  %s\n", msg);
	TEST_ASSERT(o && ex < 1e-12 && en < 1e-12 && euv < 1e-12 && bad_id == 0 && bad_lab == 0 && qaws_mesh_3d_get_attribute_count(o) == 5, msg);
	qaws_mesh_3d_destroy(o);
	qaws_mesh_3d_destroy(c);
}

static void test_refusals(void)
{
	qaws_mesh_3d *c = sd_cube(), *o = NULL;
	qaws_mesh_subdivide_desc d;
	qaws_scalar sh[12];
	qaws_status s0, s1, s2;
	char msg[300];
	unsigned int i;
	memset(&d, 0, sizeof d);
	d.scheme = QAWS_MESH_SUBDIVIDE_LOOP;
	s0 = qaws_mesh_3d_subdivide(&o, c, &d);
	qaws_mesh_3d_detect_features(c, (qaws_scalar)0.5, 0, 0);
	for (i = 0; i < 12; i++) sh[i] = (qaws_scalar)(i == 5 ? -1 : 1);
	d.scheme = QAWS_MESH_SUBDIVIDE_CATMULL_CLARK;
	d.sharpness = sh;
	s1 = qaws_mesh_3d_subdivide(&o, c, &d);
	d.sharpness = NULL;
	d.scheme = (qaws_mesh_subdivide_scheme)5;
	s2 = qaws_mesh_3d_subdivide(&o, c, &d);
	sprintf(msg, "refused: Loop on quads (%d), a negative sharpness (%d), an unknown scheme (%d)", (int)s0, (int)s1, (int)s2);
	printf("  %s\n", msg);
	TEST_ASSERT(s0 == QAWS_STATUS_INVALID_ARGUMENT && s1 == QAWS_STATUS_INVALID_ARGUMENT && s2 == QAWS_STATUS_INVALID_ARGUMENT && !o, msg);
	qaws_mesh_3d_destroy(c);
}

int test_172_mesh_subdivide_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 172: subdivision of meshes\n");
	test_counts();
	test_loop();
	test_creases();
	test_sheet();
	test_attributes();
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
