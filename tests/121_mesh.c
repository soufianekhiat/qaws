/*
 * Test 121: polygon meshes (qaws_mesh_3d over the corner store)
 *
 *   - building: an icosahedron, a cube of quads, a torus and a disk: counts,
 *     topology (components, boundaries, Euler characteristic, genus), the
 *     store's invariants;
 *   - input that is not a manifold split, never refused: a bowtie (one
 *     vertex, two fans), a fin (three faces on an edge), a neighbour turned
 *     the other way; copies traced to their source vertices;
 *   - the dual: the cube's face adjacency, its dual the octahedron;
 *   - random edits (splits, collapses, flips, face splits) on a sphere, the
 *     invariants after each, the Euler characteristic kept; a linear vertex
 *     attribute and a linear corner attribute kept exact (the positions at
 *     midpoints and centroids), a face attribute inherited;
 *   - refusals: a tetrahedron's edge (it would fold), a feature corner, a
 *     feature edge's flip, a boundary flip; nothing changed;
 *   - features: the cube's 12 creases and 8 corners, a split crease's halves;
 *   - footpoints set, kept and given to a new vertex;
 *   - quality of the icosahedron: 60 degrees, aspect 1, valence 5.
 */

#include "test_common.h"
#include "qaws_mesh.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

static qaws_mesh_3d* tm_make(qaws_scalar const* x, unsigned int nv, unsigned int const* f, unsigned int const* fs, unsigned int nf)
{
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	memset(&d, 0, sizeof d);
	d.positions = x;
	d.vertex_count = nv;
	d.face_vertices = f;
	d.face_sizes = fs;
	d.face_count = nf;
	if (qaws_mesh_3d_create(&m, &d) != QAWS_STATUS_OK)
		return NULL;
	return m;
}

static qaws_mesh_3d* tm_icosahedron(void)
{
	static unsigned int const f[60] = { 0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11, 1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8,
		3, 9, 4, 3, 4, 2, 3, 2, 6, 3, 6, 8, 3, 8, 9, 4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7, 9, 8, 1 };
	double t = (1 + sqrt(5.0)) / 2;
	qaws_scalar x[36];
	double const v[36] = { -1, t, 0, 1, t, 0, -1, -t, 0, 1, -t, 0, 0, -1, t, 0, 1, t, 0, -1, -t, 0, 1, -t, t, 0, -1, t, 0, 1, -t, 0, -1, -t, 0, 1 };
	unsigned int i;
	for (i = 0; i < 36; i++)
		x[i] = (qaws_scalar)(v[i] / sqrt(1 + t * t));
	return tm_make(x, 12, f, NULL, 20);
}

static qaws_mesh_3d* tm_cube(void)
{
	static qaws_scalar const x[24] = { -1, -1, -1, 1, -1, -1, 1, 1, -1, -1, 1, -1, -1, -1, 1, 1, -1, 1, 1, 1, 1, -1, 1, 1 };
	static unsigned int const f[24] = { 0, 3, 2, 1, 4, 5, 6, 7, 0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6, 3, 0, 4, 7 };
	static unsigned int const fs[6] = { 4, 4, 4, 4, 4, 4 };
	return tm_make(x, 8, f, fs, 6);
}

/* an nu x nv grid of quads split in triangles: a torus (closed both ways) or a disk-like patch */
static qaws_mesh_3d* tm_grid(unsigned int nu, unsigned int nv, int torus)
{
	unsigned int cu = torus ? nu : nu + 1, cv = torus ? nv : nv + 1, i, j, k = 0, n = cu * cv;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n);
	unsigned int* f = (unsigned int*)malloc(sizeof(unsigned int) * 6 * nu * nv);
	qaws_mesh_3d* m;
	for (j = 0; j < cv; j++)
		for (i = 0; i < cu; i++)
		{
			double u = 2 * M_PI * i / nu, v = 2 * M_PI * j / nv;
			qaws_scalar* p = x + 3 * (j * cu + i);
			if (torus)
			{
				p[0] = (qaws_scalar)((2 + 0.7 * cos(v)) * cos(u));
				p[1] = (qaws_scalar)((2 + 0.7 * cos(v)) * sin(u));
				p[2] = (qaws_scalar)(0.7 * sin(v));
			}
			else
			{
				p[0] = (qaws_scalar)i;
				p[1] = (qaws_scalar)j;
				p[2] = (qaws_scalar)(0.1 * sin((double)i) * cos((double)j));
			}
		}
	for (j = 0; j < nv; j++)
		for (i = 0; i < nu; i++)
		{
			unsigned int a = j * cu + i, b = j * cu + (i + 1) % cu, c = ((j + 1) % cv) * cu + (i + 1) % cu, d = ((j + 1) % cv) * cu + i;
			f[k++] = a; f[k++] = b; f[k++] = c;
			f[k++] = a; f[k++] = c; f[k++] = d;
		}
	m = tm_make(x, n, f, NULL, 2 * nu * nv);
	free(x);
	free(f);
	return m;
}

static void test_build(void)
{
	qaws_mesh_3d *ico = tm_icosahedron(), *cube = tm_cube(), *torus = tm_grid(12, 8, 1), *disk = tm_grid(6, 5, 0);
	qaws_mesh_topology t[4];
	char msg[400];
	int ok = ico && cube && torus && disk;
	TEST_ASSERT(ok, "build: the four meshes");
	if (!ok)
		goto done;
	qaws_mesh_3d_get_topology(&t[0], ico);
	qaws_mesh_3d_get_topology(&t[1], cube);
	qaws_mesh_3d_get_topology(&t[2], torus);
	qaws_mesh_3d_get_topology(&t[3], disk);
	sprintf(msg, "icosahedron: %u vertices, %u edges, %u faces, chi %d, genus %u, closed %d, valid %d",
		qaws_mesh_3d_get_vertex_count(ico), qaws_mesh_3d_get_edge_count(ico), qaws_mesh_3d_get_face_count(ico),
		t[0].euler_characteristic, t[0].genus, t[0].closed, qaws_mesh_3d_validate(ico) == QAWS_STATUS_OK);
	TEST_ASSERT(qaws_mesh_3d_get_vertex_count(ico) == 12 && qaws_mesh_3d_get_edge_count(ico) == 30 && t[0].euler_characteristic == 2
		&& t[0].genus == 0 && t[0].closed && t[0].components == 1 && qaws_mesh_3d_validate(ico) == QAWS_STATUS_OK, msg);
	sprintf(msg, "cube of quads: %u vertices, %u edges, %u faces, %u corners, %u fan triangles, chi %d, valid %d",
		qaws_mesh_3d_get_vertex_count(cube), qaws_mesh_3d_get_edge_count(cube), qaws_mesh_3d_get_face_count(cube),
		qaws_mesh_3d_get_corner_count(cube), qaws_mesh_3d_get_triangle_count(cube), t[1].euler_characteristic,
		qaws_mesh_3d_validate(cube) == QAWS_STATUS_OK);
	TEST_ASSERT(qaws_mesh_3d_get_edge_count(cube) == 12 && qaws_mesh_3d_get_corner_count(cube) == 24
		&& qaws_mesh_3d_get_triangle_count(cube) == 12 && t[1].euler_characteristic == 2 && qaws_mesh_3d_validate(cube) == QAWS_STATUS_OK, msg);
	sprintf(msg, "torus 12 x 8: chi %d, genus %u, closed %d; disk 6 x 5: chi %d, %u boundary loop, genus %u; both valid %d",
		t[2].euler_characteristic, t[2].genus, t[2].closed, t[3].euler_characteristic, t[3].boundary_loops, t[3].genus,
		qaws_mesh_3d_validate(torus) == QAWS_STATUS_OK && qaws_mesh_3d_validate(disk) == QAWS_STATUS_OK);
	TEST_ASSERT(t[2].euler_characteristic == 0 && t[2].genus == 1 && t[2].closed && t[3].euler_characteristic == 1
		&& t[3].boundary_loops == 1 && t[3].genus == 0 && qaws_mesh_3d_validate(torus) == QAWS_STATUS_OK
		&& qaws_mesh_3d_validate(disk) == QAWS_STATUS_OK, msg);
	{
		/* a fan-triangulated face read back: the cube's first face */
		unsigned int tri[3], fv[24], fs[6];
		qaws_mesh_3d_get_faces(fv, fs, cube);
		qaws_mesh_3d_get_triangle(tri, cube, 1);
		sprintf(msg, "cube: face 0 (%u %u %u %u), its second fan triangle (%u %u %u)", fv[0], fv[1], fv[2], fv[3], tri[0], tri[1], tri[2]);
		TEST_ASSERT(fs[0] == 4 && tri[0] == fv[0] && tri[1] == fv[2] && tri[2] == fv[3], msg);
	}
	{
		/* bad input refused: a face with a repeated vertex, an index out of range */
		static qaws_scalar const x[9] = { 0, 0, 0, 1, 0, 0, 0, 1, 0 };
		static unsigned int const rep[3] = { 0, 1, 1 }, out[3] = { 0, 1, 3 };
		qaws_mesh_3d *a = tm_make(x, 3, rep, NULL, 1), *b = tm_make(x, 3, out, NULL, 1);
		TEST_ASSERT(!a && !b, "build refusals: a repeated vertex, an index out of range");
		qaws_mesh_3d_destroy(a);
		qaws_mesh_3d_destroy(b);
	}
done:
	qaws_mesh_3d_destroy(ico);
	qaws_mesh_3d_destroy(cube);
	qaws_mesh_3d_destroy(torus);
	qaws_mesh_3d_destroy(disk);
}

static void test_non_manifold(void)
{
	static qaws_scalar const x[15] = { 0, 0, 0, 1, 0, 0, 0, 1, 0, -1, 0, 0, 0, -1, 0 };
	static unsigned int const bowtie[6] = { 0, 1, 2, 0, 3, 4 };
	static unsigned int const fin[9] = { 0, 1, 2, 1, 0, 4, 0, 1, 3 };
	static unsigned int const turned[6] = { 0, 1, 2, 0, 1, 4 };
	qaws_mesh_3d *a = tm_make(x, 5, bowtie, NULL, 2), *b = tm_make(x, 5, fin, NULL, 3), *c = tm_make(x, 5, turned, NULL, 2);
	qaws_mesh_topology ta, tb;
	char msg[400];
	int ok = a && b && c;
	TEST_ASSERT(ok, "non-manifold input accepted");
	if (!ok)
		goto done;
	qaws_mesh_3d_get_topology(&ta, a);
	qaws_mesh_3d_get_topology(&tb, b);
	sprintf(msg, "bowtie: %u vertices (copy of %u), %u components, valid %d; fin: %u vertices (copies of %u, %u), valid %d; turned neighbour: %u vertices, valid %d",
		qaws_mesh_3d_get_vertex_count(a), qaws_mesh_3d_get_vertex_source(a, 5), ta.components, qaws_mesh_3d_validate(a) == QAWS_STATUS_OK,
		qaws_mesh_3d_get_vertex_count(b), qaws_mesh_3d_get_vertex_source(b, 5), qaws_mesh_3d_get_vertex_source(b, 6),
		qaws_mesh_3d_validate(b) == QAWS_STATUS_OK, qaws_mesh_3d_get_vertex_count(c), qaws_mesh_3d_validate(c) == QAWS_STATUS_OK);
	TEST_ASSERT(qaws_mesh_3d_get_vertex_count(a) == 6 && qaws_mesh_3d_get_vertex_source(a, 5) == 0 && ta.components == 2
		&& qaws_mesh_3d_validate(a) == QAWS_STATUS_OK && qaws_mesh_3d_get_vertex_count(b) == 7
		&& qaws_mesh_3d_get_vertex_source(b, 5) == 0 && qaws_mesh_3d_get_vertex_source(b, 6) == 1 && tb.components == 2
		&& qaws_mesh_3d_validate(b) == QAWS_STATUS_OK && qaws_mesh_3d_get_vertex_count(c) == 7 && qaws_mesh_3d_validate(c) == QAWS_STATUS_OK, msg);
done:
	qaws_mesh_3d_destroy(a);
	qaws_mesh_3d_destroy(b);
	qaws_mesh_3d_destroy(c);
}

static void test_dual(void)
{
	qaws_mesh_3d *cube = tm_cube(), *dual = NULL;
	unsigned int adj[24], fv[24], fs[6], i, j, ok = 1;
	qaws_mesh_topology t;
	char msg[300];
	if (!cube)
	{
		TEST_ASSERT(0, "dual: the cube");
		return;
	}
	qaws_mesh_3d_get_face_adjacency(adj, cube);
	qaws_mesh_3d_get_faces(fv, fs, cube);
	/* each face's neighbour across each edge shares that edge, reversed */
	for (i = 0; i < 24; i++)
	{
		unsigned int f = i / 4, a = fv[i], b = fv[f * 4 + (i % 4 + 1) % 4], g = adj[i], hit = 0;
		if (g == QAWS_MESH_NONE || g == f)
		{
			ok = 0;
			continue;
		}
		for (j = 0; j < 4; j++)
			hit |= fv[g * 4 + j] == b && fv[g * 4 + (j + 1) % 4] == a;
		ok &= hit;
	}
	TEST_ASSERT(ok, "dual: the cube's face adjacency, each neighbour on the reversed edge");
	TEST_ASSERT(qaws_mesh_3d_create_dual(&dual, cube) == QAWS_STATUS_OK && dual, "dual: created");
	if (dual)
	{
		qaws_mesh_3d_get_topology(&t, dual);
		sprintf(msg, "the cube's dual: %u vertices, %u faces (%u triangles), %u edges, chi %d, valid %d", qaws_mesh_3d_get_vertex_count(dual),
			qaws_mesh_3d_get_face_count(dual), qaws_mesh_3d_get_triangle_count(dual), qaws_mesh_3d_get_edge_count(dual), t.euler_characteristic,
			qaws_mesh_3d_validate(dual) == QAWS_STATUS_OK);
		TEST_ASSERT(qaws_mesh_3d_get_vertex_count(dual) == 6 && qaws_mesh_3d_get_face_count(dual) == 8 && qaws_mesh_3d_get_triangle_count(dual) == 8
			&& qaws_mesh_3d_get_edge_count(dual) == 12 && t.euler_characteristic == 2 && t.closed && qaws_mesh_3d_validate(dual) == QAWS_STATUS_OK, msg);
	}
	qaws_mesh_3d_destroy(dual);
	qaws_mesh_3d_destroy(cube);
}

static unsigned int tm_rand(unsigned int* s)
{
	*s = *s * 1664525u + 1013904223u;
	return *s >> 8;
}

static double tm_f(double const* p) { return p[0] + 2 * p[1] - 3 * p[2]; }
static double tm_g(double const* p) { return 2 * p[0] - p[1] + 0.5 * p[2]; }

/* the attributes against their linear functions: vertex values, corner values */
static void tm_attr_error(qaws_mesh_3d* m, unsigned int va, unsigned int ca, double* out_v, double* out_c)
{
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), nc = qaws_mesh_3d_get_corner_count(m), nf = qaws_mesh_3d_get_face_count(m), i;
	qaws_scalar *x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv), *vv = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv),
		*cc = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nc);
	unsigned int *fv = (unsigned int*)malloc(sizeof(unsigned int) * nc), *fs = (unsigned int*)malloc(sizeof(unsigned int) * nf);
	*out_v = *out_c = 0;
	qaws_mesh_3d_get_positions(x, m);
	qaws_mesh_3d_get_attribute(vv, m, va);
	qaws_mesh_3d_get_faces(fv, fs, m);
	for (i = 0; i < nv; i++)
	{
		double p[3] = { x[3 * i], x[3 * i + 1], x[3 * i + 2] };
		*out_v = fmax(*out_v, fabs(vv[i] - tm_f(p)));
	}
	if (ca != QAWS_MESH_NONE)
	{
		qaws_mesh_3d_get_attribute(cc, m, ca);
		for (i = 0; i < nc; i++)
		{
			double p[3] = { x[3 * fv[i]], x[3 * fv[i] + 1], x[3 * fv[i] + 2] };
			*out_c = fmax(*out_c, fabs(cc[i] - tm_g(p)));
		}
	}
	free(x); free(vv); free(cc); free(fv); free(fs);
}

/* attributes from the linear functions, on mesh m */
static void tm_add_attrs(qaws_mesh_3d* m, unsigned int* va, unsigned int* ca, unsigned int* fa)
{
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), nc = qaws_mesh_3d_get_corner_count(m), nf = qaws_mesh_3d_get_face_count(m), i;
	qaws_scalar *x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv), *vv = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv),
		*cc = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nc), *ff = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nf);
	unsigned int *fv = (unsigned int*)malloc(sizeof(unsigned int) * nc), *fs = (unsigned int*)malloc(sizeof(unsigned int) * nf);
	qaws_attribute_desc d;
	qaws_mesh_3d_get_positions(x, m);
	qaws_mesh_3d_get_faces(fv, fs, m);
	for (i = 0; i < nv; i++)
	{
		double p[3] = { x[3 * i], x[3 * i + 1], x[3 * i + 2] };
		vv[i] = (qaws_scalar)tm_f(p);
	}
	for (i = 0; i < nc; i++)
	{
		double p[3] = { x[3 * fv[i]], x[3 * fv[i] + 1], x[3 * fv[i] + 2] };
		cc[i] = (qaws_scalar)tm_g(p);
	}
	for (i = 0; i < nf; i++)
		ff[i] = (qaws_scalar)(i % 3);
	memset(&d, 0, sizeof d);
	d.name = "f"; d.element = QAWS_MESH_VERTEX; d.components = 1; d.rule = QAWS_ATTRIBUTE_LINEAR;
	qaws_mesh_3d_add_attribute(va, m, vv, &d);
	d.name = "g"; d.element = QAWS_MESH_CORNER;
	qaws_mesh_3d_add_attribute(ca, m, cc, &d);
	d.name = "material"; d.element = QAWS_MESH_FACE; d.rule = QAWS_ATTRIBUTE_NEAREST;
	qaws_mesh_3d_add_attribute(fa, m, ff, &d);
	free(x); free(vv); free(cc); free(ff); free(fv); free(fs);
}

/* n random edits; with collapses or not (corner attributes are kept by a collapse, not re-interpolated) */
static void tm_edit(qaws_mesh_3d* m, unsigned int n, int collapses, unsigned int seed, unsigned int* out_ops, int* out_valid, int* out_chi)
{
	unsigned int it, s = seed, ops[4] = { 0, 0, 0, 0 };
	qaws_mesh_topology t0, t1;
	qaws_mesh_3d_get_topology(&t0, m);
	*out_valid = 1;
	for (it = 0; it < n && *out_valid; it++)
	{
		unsigned int ne = qaws_mesh_3d_get_edge_count(m), *e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne), k = tm_rand(&s) % ne;
		unsigned int a, b, op = tm_rand(&s) % (collapses ? 4 : 3), v;
		qaws_status st = QAWS_STATUS_OK;
		qaws_vec3 pa, pb;
		qaws_scalar mid[3];
		qaws_mesh_3d_get_edges(e, m);
		a = e[2 * k];
		b = e[2 * k + 1];
		free(e);
		pa = qaws_mesh_3d_get_vertex(m, a);
		pb = qaws_mesh_3d_get_vertex(m, b);
		mid[0] = (pa.x + pb.x) / 2; mid[1] = (pa.y + pb.y) / 2; mid[2] = (pa.z + pb.z) / 2;
		if (qaws_mesh_3d_get_vertex_count(m) > 400)
			op = collapses ? 3 : 1;   /* (keep it small) */
		if (op == 0)
			st = qaws_mesh_3d_split_edge(&v, m, mid, NULL, a, b);
		else if (op == 1)
			st = qaws_mesh_3d_can_flip_edge(m, a, b) ? qaws_mesh_3d_flip_edge(m, a, b) : QAWS_STATUS_INVALID_TOPOLOGY;
		else if (op == 2)
			st = qaws_mesh_3d_split_face(&v, m, NULL, NULL, tm_rand(&s) % qaws_mesh_3d_get_face_count(m));
		else
			st = qaws_mesh_3d_can_collapse_edge(m, a, b) ? qaws_mesh_3d_collapse_edge(m, a, b, mid, NULL) : QAWS_STATUS_INVALID_TOPOLOGY;
		if (st == QAWS_STATUS_OK)
			ops[op]++;
		else if (st != QAWS_STATUS_INVALID_TOPOLOGY)
			*out_valid = 0;
		if (qaws_mesh_3d_validate(m) != QAWS_STATUS_OK)
			*out_valid = 0;
	}
	qaws_mesh_3d_get_topology(&t1, m);
	*out_chi = t1.euler_characteristic == t0.euler_characteristic && t1.closed == t0.closed;
	memcpy(out_ops, ops, sizeof ops);
}

static void test_edits(void)
{
	qaws_mesh_3d* m = tm_icosahedron();
	unsigned int va, ca, fa, ops[4], nf, i;
	int valid, chi, materials = 1;
	double ev, ec;
	qaws_scalar* ff;
	char msg[500];
	if (!m)
	{
		TEST_ASSERT(0, "edits: the icosahedron");
		return;
	}
	tm_add_attrs(m, &va, &ca, &fa);
	tm_edit(m, 600, 0, 7, ops, &valid, &chi);
	tm_attr_error(m, va, ca, &ev, &ec);
	nf = qaws_mesh_3d_get_face_count(m);
	ff = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nf);
	qaws_mesh_3d_get_attribute(ff, m, fa);
	for (i = 0; i < nf; i++)
		materials &= ff[i] == 0 || ff[i] == 1 || ff[i] == 2;
	free(ff);
	sprintf(msg, "600 edits without collapses (%u splits, %u flips, %u face splits): valid after each %d, chi kept %d, %u vertices; vertex attribute error %.2g, corner attribute error %.2g, materials inherited %d",
		ops[0], ops[1], ops[2], valid, chi, qaws_mesh_3d_get_vertex_count(m), ev, ec, materials);
	TEST_ASSERT(valid && chi && ops[0] > 50 && ops[1] > 20 && ops[2] > 50 && ev < 1e-4 && ec < 1e-4 && materials, msg);
	tm_edit(m, 3000, 1, 11, ops, &valid, &chi);
	tm_attr_error(m, va, QAWS_MESH_NONE, &ev, &ec);
	sprintf(msg, "3000 edits with collapses (%u splits, %u flips, %u face splits, %u collapses): valid after each %d, chi kept %d, %u vertices; vertex attribute error %.2g",
		ops[0], ops[1], ops[2], ops[3], valid, chi, qaws_mesh_3d_get_vertex_count(m), ev);
	TEST_ASSERT(valid && chi && ops[3] > 200 && ev < 1e-4, msg);
	qaws_mesh_3d_destroy(m);
	/* the same on a torus with a border (a cut grid) */
	m = tm_grid(10, 6, 0);
	if (m)
	{
		tm_add_attrs(m, &va, &ca, &fa);
		tm_edit(m, 2000, 1, 5, ops, &valid, &chi);
		tm_attr_error(m, va, QAWS_MESH_NONE, &ev, &ec);
		sprintf(msg, "2000 edits on an open grid (%u splits, %u flips, %u face splits, %u collapses): valid after each %d, chi and border kept %d, vertex attribute error %.2g",
			ops[0], ops[1], ops[2], ops[3], valid, chi, ev);
		TEST_ASSERT(valid && chi && ops[3] > 100 && ev < 1e-4, msg);
	}
	qaws_mesh_3d_destroy(m);
}

static void test_refusals_and_features(void)
{
	static qaws_scalar const tx[12] = { 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1 };
	static unsigned int const tf[12] = { 0, 2, 1, 0, 1, 3, 0, 3, 2, 1, 2, 3 };
	qaws_mesh_3d *tet = tm_make(tx, 4, tf, NULL, 4), *cube = tm_cube(), *disk = tm_grid(3, 3, 0);
	unsigned int f, v, nfe, corners = 0, e[64];
	qaws_status s1, s2, s3, s4;
	char msg[400];
	if (!tet || !cube || !disk)
	{
		TEST_ASSERT(0, "refusals: the meshes");
		goto done;
	}
	s1 = qaws_mesh_3d_collapse_edge(tet, 0, 1, NULL, NULL);
	s2 = qaws_mesh_3d_collapse_edge(tet, 0, 0, NULL, NULL);
	s3 = qaws_mesh_3d_flip_edge(disk, 0, 1);   /* (a border edge) */
	sprintf(msg, "refusals: a tetrahedron's edge collapse %s (still 4 vertices, valid %d), a non-edge %s, a border flip %s",
		qaws_status_to_string(s1), qaws_mesh_3d_validate(tet) == QAWS_STATUS_OK, qaws_status_to_string(s2), qaws_status_to_string(s3));
	TEST_ASSERT(s1 == QAWS_STATUS_INVALID_TOPOLOGY && qaws_mesh_3d_get_vertex_count(tet) == 4 && s2 == QAWS_STATUS_INVALID_ARGUMENT
		&& s3 == QAWS_STATUS_INVALID_TOPOLOGY, msg);
	/* the cube's creases: triangulated by face splits, then detected */
	for (f = 0; f < 6; f++)
		qaws_mesh_3d_split_face(NULL, cube, NULL, NULL, f);
	qaws_mesh_3d_detect_features(cube, (qaws_scalar)(30 * M_PI / 180), 0, 0);
	nfe = qaws_mesh_3d_get_feature_edge_count(cube);
	for (v = 0; v < qaws_mesh_3d_get_vertex_count(cube); v++)
		corners += qaws_mesh_3d_is_feature_corner(cube, v);
	qaws_mesh_3d_get_feature_edges(e, cube);
	s1 = qaws_mesh_3d_collapse_edge(cube, e[0], e[1], NULL, NULL);   /* (a corner into a corner) */
	s2 = qaws_mesh_3d_flip_edge(cube, e[0], e[1]);
	s3 = qaws_mesh_3d_collapse_edge(cube, 8, 0, NULL, NULL);          /* (face 0's centre into its corner 0: allowed) */
	sprintf(msg, "cube features: %u creases, %u corners; refused: collapsing a corner %s, flipping a crease %s; a face centre into a corner %s",
		nfe, corners, qaws_status_to_string(s1), qaws_status_to_string(s2), qaws_status_to_string(s3));
	TEST_ASSERT(nfe == 12 && corners == 8 && s1 == QAWS_STATUS_INVALID_TOPOLOGY && s2 == QAWS_STATUS_INVALID_TOPOLOGY
		&& s3 == QAWS_STATUS_OK && qaws_mesh_3d_validate(cube) == QAWS_STATUS_OK, msg);
	qaws_mesh_3d_get_feature_edges(e, cube);
	s4 = qaws_mesh_3d_split_edge(&v, cube, NULL, NULL, e[0], e[1]);
	sprintf(msg, "a crease split: %s, %u creases now, the new vertex a corner %d (no), valid %d", qaws_status_to_string(s4),
		qaws_mesh_3d_get_feature_edge_count(cube), qaws_mesh_3d_is_feature_corner(cube, v), qaws_mesh_3d_validate(cube) == QAWS_STATUS_OK);
	TEST_ASSERT(s4 == QAWS_STATUS_OK && qaws_mesh_3d_get_feature_edge_count(cube) == 13 && !qaws_mesh_3d_is_feature_corner(cube, v)
		&& qaws_mesh_3d_validate(cube) == QAWS_STATUS_OK, msg);
done:
	qaws_mesh_3d_destroy(tet);
	qaws_mesh_3d_destroy(cube);
	qaws_mesh_3d_destroy(disk);
}

static void test_footpoints_and_quality(void)
{
	qaws_mesh_3d* m = tm_icosahedron();
	qaws_mesh_location loc[13], got[13], one;
	qaws_mesh_quality q;
	unsigned int i, v = 0, ok = 1;
	char msg[400];
	if (!m)
	{
		TEST_ASSERT(0, "footpoints: the icosahedron");
		return;
	}
	qaws_mesh_3d_measure_quality(&q, m);
	sprintf(msg, "icosahedron quality: angles %.6g to %.6g, aspect %.6g to %.6g, valence 5 at %u vertices, %u irregular", (double)q.min_angle,
		(double)q.max_angle, (double)q.mean_aspect, (double)q.max_aspect, q.valence[2], q.irregular);
	TEST_ASSERT(fabs(q.min_angle - 60) < 1e-3 && fabs(q.max_angle - 60) < 1e-3 && fabs(q.max_aspect - 1) < 1e-4 && q.valence[2] == 12
		&& q.irregular == 12 && q.triangles == 20, msg);
	TEST_ASSERT(!qaws_mesh_3d_has_footpoints(m) && qaws_mesh_3d_get_footpoints(got, m) == QAWS_STATUS_INVALID_ARGUMENT, "footpoints: none at first");
	for (i = 0; i < 12; i++)
	{
		loc[i].patch = i;
		loc[i].a = (qaws_scalar)(0.1 * i);
		loc[i].b = (qaws_scalar)0.25;
	}
	qaws_mesh_3d_set_footpoints(m, loc);
	one.patch = 99; one.a = (qaws_scalar)0.5; one.b = (qaws_scalar)0.5;
	qaws_mesh_3d_split_edge(&v, m, NULL, &one, 0, 11);
	qaws_mesh_3d_get_footpoints(got, m);
	for (i = 0; i < 12; i++)
		ok &= got[i].patch == loc[i].patch && got[i].a == loc[i].a && got[i].b == loc[i].b;
	sprintf(msg, "footpoints kept through a split (%u), the new vertex %u given its own (patch %u)", ok, v, got[v].patch);
	TEST_ASSERT(ok && v == 12 && got[12].patch == 99, msg);
	qaws_mesh_3d_destroy(m);
}

int test_121_mesh_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 121: polygon meshes over the corner store\n");
	test_build();
	test_non_manifold();
	test_dual();
	test_edits();
	test_refusals_and_features();
	test_footpoints_and_quality();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
