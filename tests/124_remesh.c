/*
 * Test 124: remeshing (qaws_mesh_remesh, LOCAL and FIELD_ALIGNED)
 *
 *   - a coarse icosphere with its normals remeshed on its PN reference: a
 *     valid closed mesh of genus 0, every vertex its footpoint evaluated,
 *     near the sphere, edges near the target, angles above 25 degrees,
 *     mostly of valence 6, the faces' centres near the reference;
 *   - a cube of quads with its creases: every vertex on the cube, its 8
 *     corners kept, the crease edges along the cube's 12 edges (their
 *     length 24), a linear vertex attribute exact at every vertex, a corner
 *     attribute with seams on the creases (a value per side) kept per side;
 *   - an open bumpy square: one border loop, the border on the square's
 *     border, its 4 corners kept;
 *   - adaptive: a curvature field on a bumped sphere, shorter edges on the
 *     bump;
 *   - FIELD_ALIGNED on the sphere and on the cube with its creases: valid
 *     closed meshes of genus 0, footpoints, near the surface, edges near the
 *     target, mostly of valence 6 (the hexagonal lattice), the cube's
 *     corners kept and its creases along its edges;
 *   - refusals (FIELD_ALIGNED with a metric, a size field on another reference).
 */

#include "test_common.h"
#include "qaws_remesh.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

static qaws_mesh_3d* rt_make(double const* x, unsigned int nv, unsigned int const* f, unsigned int const* fs, unsigned int nf)
{
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	qaws_scalar* p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	unsigned int i;
	for (i = 0; i < 3 * nv; i++)
		p[i] = (qaws_scalar)x[i];
	memset(&d, 0, sizeof d);
	d.positions = p;
	d.vertex_count = nv;
	d.face_vertices = f;
	d.face_sizes = fs;
	d.face_count = nf;
	qaws_mesh_3d_create(&m, &d);
	free(p);
	return m;
}

/* the icosahedron split once, on the unit sphere, with its radial normals */
static qaws_mesh_3d* rt_icosphere(qaws_scalar* normals)
{
	static unsigned int const f0[60] = { 0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11, 1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8,
		3, 9, 4, 3, 4, 2, 3, 2, 6, 3, 6, 8, 3, 8, 9, 4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7, 9, 8, 1 };
	double g = (1 + sqrt(5.0)) / 2, x[3 * 42];
	double const v0[36] = { -1, g, 0, 1, g, 0, -1, -g, 0, 1, -g, 0, 0, -1, g, 0, 1, g, 0, -1, -g, 0, 1, -g, g, 0, -1, g, 0, 1, -g, 0, -1, -g, 0, 1 };
	unsigned int f[240], keys[3 * 30], nk = 0, nv = 12, nf = 0, i, j, q, k;
	memcpy(x, v0, sizeof v0);
	for (i = 0; i < 20; i++)
	{
		unsigned int e[3];
		for (j = 0; j < 3; j++)
		{
			unsigned int a = f0[3 * i + j], b = f0[3 * i + (j + 1) % 3], lo = a < b ? a : b, hi = a < b ? b : a;
			e[j] = 0xFFFFFFFFu;
			for (q = 0; q < nk; q++)
				if (keys[3 * q] == lo && keys[3 * q + 1] == hi)
					e[j] = keys[3 * q + 2];
			if (e[j] == 0xFFFFFFFFu)
			{
				for (k = 0; k < 3; k++)
					x[3 * nv + k] = 0.5 * (x[3 * a + k] + x[3 * b + k]);
				keys[3 * nk] = lo; keys[3 * nk + 1] = hi; keys[3 * nk + 2] = nv;
				e[j] = nv++;
				nk++;
			}
		}
		{
			unsigned int a = f0[3 * i], b = f0[3 * i + 1], c = f0[3 * i + 2];
			unsigned int t[12] = { a, e[0], e[2], e[0], b, e[1], e[2], e[1], c, e[0], e[1], e[2] };
			memcpy(f + 3 * nf, t, sizeof t);
			nf += 4;
		}
	}
	for (i = 0; i < nv; i++)
	{
		double r = sqrt(x[3 * i] * x[3 * i] + x[3 * i + 1] * x[3 * i + 1] + x[3 * i + 2] * x[3 * i + 2]);
		for (k = 0; k < 3; k++)
		{
			x[3 * i + k] /= r;
			if (normals)
				normals[3 * i + k] = (qaws_scalar)x[3 * i + k];
		}
	}
	return rt_make(x, nv, f, NULL, nf);
}

static qaws_mesh_3d* rt_remesh(qaws_mesh_3d const* m, qaws_mesh_reference const* r, double L, qaws_size_field const* sf, qaws_mesh_report* rep,
	qaws_status* st)
{
	qaws_remesh_desc d;
	qaws_mesh_3d* out = NULL;
	memset(&d, 0, sizeof d);
	d.reference = r;
	d.target_length = (qaws_scalar)L;
	d.size = sf;
	*st = qaws_mesh_remesh(&out, rep, m, &d);
	return out;
}

static void test_sphere(void)
{
	qaws_scalar nrm[3 * 42];
	qaws_mesh_3d *m = rt_icosphere(nrm), *out = NULL;
	qaws_mesh_reference_desc rd;
	qaws_mesh_reference* r = NULL;
	qaws_mesh_report rep = { 0 };
	qaws_mesh_topology t;
	qaws_status st = QAWS_STATUS_INTERNAL_ERROR;
	unsigned int nv = 0, i, v6 = 0, interior = 0;
	double foot = 0, off = 0;
	char msg[500];
	memset(&rd, 0, sizeof rd);
	rd.normals = nrm;
	qaws_mesh_reference_create(&r, m, &rd);
	if (r)
		out = rt_remesh(m, r, 0.12, NULL, &rep, &st);
	if (!out)
	{
		TEST_ASSERT(0, "sphere: remeshed");
		goto done;
	}
	nv = qaws_mesh_3d_get_vertex_count(out);
	{
		qaws_mesh_location* l = (qaws_mesh_location*)malloc(sizeof(qaws_mesh_location) * nv);
		qaws_scalar *x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv), *y = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
		qaws_mesh_3d_get_positions(x, out);
		qaws_mesh_3d_get_footpoints(l, out);
		qaws_mesh_reference_evaluate(y, NULL, NULL, r, l, nv);
		for (i = 0; i < nv; i++)
		{
			unsigned int ring[64], n = 0;
			foot = fmax(foot, fabs(x[3 * i] - y[3 * i]) + fabs(x[3 * i + 1] - y[3 * i + 1]) + fabs(x[3 * i + 2] - y[3 * i + 2]));
			off = fmax(off, fabs(sqrt(x[3 * i] * x[3 * i] + x[3 * i + 1] * x[3 * i + 1] + x[3 * i + 2] * x[3 * i + 2]) - 1));
			qaws_mesh_3d_get_vertex_ring(ring, &n, out, 64, i);
			interior++;
			v6 += n == 6;
		}
		free(l); free(x); free(y);
	}
	qaws_mesh_3d_get_topology(&t, out);
	sprintf(msg, "sphere remeshed to 0.12: %u vertices, valid %d, closed %d, genus %u; footpoints %.2g, off the sphere %.2g; edges %.3g..%.3g mean %.3g, angles %.3g..%.3g, valence 6 at %.0f%%, deviation %.2g",
		nv, qaws_mesh_3d_validate(out) == QAWS_STATUS_OK, t.closed, t.genus, foot, off, (double)rep.quality.min_edge, (double)rep.quality.max_edge,
		(double)rep.quality.mean_edge, (double)rep.quality.min_angle, (double)rep.quality.max_angle, 100.0 * v6 / interior, (double)rep.max_deviation);
	TEST_ASSERT(st == QAWS_STATUS_OK && qaws_mesh_3d_validate(out) == QAWS_STATUS_OK && t.closed && t.genus == 0 && t.components == 1 && foot < 1e-9
		&& off < 0.01 && fabs(rep.quality.mean_edge - 0.12) < 0.025 && rep.quality.min_angle > 25 && v6 > 0.7 * interior && rep.max_deviation < 0.01, msg);
done:
	qaws_mesh_reference_destroy(r);
	qaws_mesh_3d_destroy(m);
	qaws_mesh_3d_destroy(out);
}

static void test_cube(void)
{
	static double const x[24] = { -1, -1, -1, 1, -1, -1, 1, 1, -1, -1, 1, -1, -1, -1, 1, 1, -1, 1, 1, 1, 1, -1, 1, 1 };
	static unsigned int const f[24] = { 0, 3, 2, 1, 4, 5, 6, 7, 0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6, 3, 0, 4, 7 };
	static unsigned int const fs[6] = { 4, 4, 4, 4, 4, 4 };
	qaws_mesh_3d *m = rt_make(x, 8, f, fs, 6), *out = NULL;
	qaws_mesh_report rep;
	qaws_status st = QAWS_STATUS_INTERNAL_ERROR;
	qaws_attribute_desc ad;
	qaws_scalar vals[24];
	unsigned int i, va = 0, ca = 0, nv, nfe, corners = 0, *fe, ok_seam = 1;
	double onc = 0, flen = 0, fon = 0, ev = 0;
	char msg[500];
	for (i = 0; i < 8; i++)
		vals[i] = (qaws_scalar)(x[3 * i] + 2 * x[3 * i + 1] + 3 * x[3 * i + 2]);
	memset(&ad, 0, sizeof ad);
	ad.name = "f"; ad.element = QAWS_MESH_VERTEX; ad.components = 1;
	qaws_mesh_3d_add_attribute(&va, m, vals, &ad);
	for (i = 0; i < 24; i++)
		vals[i] = (qaws_scalar)(i / 4);   /* (each face's corners: its index) */
	ad.name = "side"; ad.element = QAWS_MESH_CORNER; ad.rule = QAWS_ATTRIBUTE_NEAREST;
	qaws_mesh_3d_add_attribute(&ca, m, vals, &ad);
	qaws_mesh_3d_detect_features(m, (qaws_scalar)(M_PI / 6), 0, 0);
	out = rt_remesh(m, NULL, 0.25, NULL, &rep, &st);
	if (!out)
	{
		TEST_ASSERT(0, "cube: remeshed");
		goto done;
	}
	nv = qaws_mesh_3d_get_vertex_count(out);
	nfe = qaws_mesh_3d_get_feature_edge_count(out);
	{
		qaws_scalar *p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv), *fv = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv);
		unsigned int found[8] = { 0 }, a2, c2, nc = qaws_mesh_3d_get_corner_count(out), nf = qaws_mesh_3d_get_face_count(out), k, q;
		unsigned int *fvx = (unsigned int*)malloc(sizeof(unsigned int) * nc), *fsz = (unsigned int*)malloc(sizeof(unsigned int) * nf);
		qaws_scalar* side = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nc);
		fe = (unsigned int*)malloc(sizeof(unsigned int) * 2 * nfe);
		qaws_mesh_3d_get_positions(p, out);
		qaws_mesh_3d_find_attribute(&a2, out, "f");
		qaws_mesh_3d_find_attribute(&c2, out, "side");
		qaws_mesh_3d_get_attribute(fv, out, a2);
		for (i = 0; i < nv; i++)
		{
			double mx = fmax(fabs(p[3 * i]), fmax(fabs(p[3 * i + 1]), fabs(p[3 * i + 2])));
			onc = fmax(onc, fabs(mx - 1));
			ev = fmax(ev, fabs(fv[i] - (p[3 * i] + 2 * p[3 * i + 1] + 3 * p[3 * i + 2])));
			for (k = 0; k < 8; k++)
				if (fabs(p[3 * i] - x[3 * k]) + fabs(p[3 * i + 1] - x[3 * k + 1]) + fabs(p[3 * i + 2] - x[3 * k + 2]) < 1e-9)
					found[k] = 1;
			corners += qaws_mesh_3d_is_feature_corner(out, i);
		}
		qaws_mesh_3d_get_feature_edges(fe, out);
		for (i = 0; i < nfe; i++)
		{
			qaws_scalar const *a = p + 3 * fe[2 * i], *b = p + 3 * fe[2 * i + 1];
			unsigned int on = 0;
			flen += sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
			for (k = 0; k < 3; k++)
				on += fabs(fabs(a[k]) - 1) < 1e-9 && fabs(fabs(b[k]) - 1) < 1e-9;
			fon += on < 2;
		}
		/* each face's corners: the side its centre is on */
		qaws_mesh_3d_get_faces(fvx, fsz, out);
		qaws_mesh_3d_get_attribute(side, out, c2);
		for (i = 0, q = 0; i < nf; q += fsz[i], i++)
		{
			double cen[3] = { 0, 0, 0 };
			unsigned int j, want = 0;
			static unsigned int const axis_side[3][2] = { { 5, 3 }, { 2, 4 }, { 0, 1 } };   /* (-x: face 5, +x: 3; -y: 2, +y: 4; -z: 0, +z: 1) */
			for (j = 0; j < fsz[i]; j++)
				for (k = 0; k < 3; k++)
					cen[k] += p[3 * fvx[q + j] + k] / fsz[i];
			for (k = 0; k < 3; k++)
				if (fabs(cen[k]) > 1 - 1e-6)
					want = axis_side[k][cen[k] > 0];
			for (j = 0; j < fsz[i]; j++)
				ok_seam &= (unsigned int)side[q + j] == want;
		}
		corners = corners >= 8 && found[0] && found[1] && found[2] && found[3] && found[4] && found[5] && found[6] && found[7];
		free(p); free(fv); free(fe); free(fvx); free(fsz); free(side);
	}
	sprintf(msg, "cube remeshed to 0.25: %u vertices, valid %d, on the cube %.2g, its corners kept %u, %u crease edges of total length %.12g (off the cube's edges %.0f), linear attribute %.2g, sides per corner %u, min angle %.3g",
		nv, qaws_mesh_3d_validate(out) == QAWS_STATUS_OK, onc, corners, nfe, flen, fon, ev, ok_seam, (double)rep.quality.min_angle);
	TEST_ASSERT(st == QAWS_STATUS_OK && qaws_mesh_3d_validate(out) == QAWS_STATUS_OK && onc < TEST_BOUND(1e-9, 32) && corners && fabs(flen - 24) < TEST_BOUND(1e-9, 32 * 24) && fon == 0
		&& ev < TEST_BOUND(1e-9, 32) && ok_seam && rep.quality.min_angle > 20, msg);
done:
	qaws_mesh_3d_destroy(m);
	qaws_mesh_3d_destroy(out);
}

static void test_open_square(void)
{
	unsigned int n = 6, i, j, k = 0, nb = 0, okb = 1, found = 0;
	double x[3 * 49];
	unsigned int f[6 * 36];
	qaws_mesh_3d *m, *out = NULL;
	qaws_mesh_report rep;
	qaws_mesh_topology t;
	qaws_status st = QAWS_STATUS_INTERNAL_ERROR;
	char msg[300];
	for (j = 0; j <= n; j++)
		for (i = 0; i <= n; i++)
		{
			double u = (double)i / n, v = (double)j / n;
			x[3 * (j * (n + 1) + i)] = u;
			x[3 * (j * (n + 1) + i) + 1] = v;
			x[3 * (j * (n + 1) + i) + 2] = 0.15 * sin(3 * u) * cos(2 * v);
		}
	for (j = 0; j < n; j++)
		for (i = 0; i < n; i++)
		{
			unsigned int a = j * (n + 1) + i;
			f[k++] = a; f[k++] = a + 1; f[k++] = a + n + 2;
			f[k++] = a; f[k++] = a + n + 2; f[k++] = a + n + 1;
		}
	m = rt_make(x, (n + 1) * (n + 1), f, NULL, 2 * n * n);
	out = rt_remesh(m, NULL, 0.06, NULL, &rep, &st);
	if (!out)
	{
		TEST_ASSERT(0, "open square: remeshed");
		goto done;
	}
	qaws_mesh_3d_get_topology(&t, out);
	{
		unsigned int nv = qaws_mesh_3d_get_vertex_count(out);
		qaws_scalar* p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
		qaws_mesh_3d_get_positions(p, out);
		for (i = 0; i < nv; i++)
		{
			double u = p[3 * i], v = p[3 * i + 1];
			int corner = (fabs(u) < 1e-9 || fabs(u - 1) < 1e-9) && (fabs(v) < 1e-9 || fabs(v - 1) < 1e-9);
			found += corner;
			if (!qaws_mesh_3d_is_boundary_vertex(out, i))
				continue;
			nb++;
			okb &= fabs(u) < 1e-9 || fabs(u - 1) < 1e-9 || fabs(v) < 1e-9 || fabs(v - 1) < 1e-9;
		}
		free(p);
	}
	sprintf(msg, "an open bumpy square remeshed to 0.06: valid %d, %u border loop, %u border vertices all on the square's border %u, its 4 corners kept (%u), min angle %.3g",
		qaws_mesh_3d_validate(out) == QAWS_STATUS_OK, t.boundary_loops, nb, okb, found, (double)rep.quality.min_angle);
	TEST_ASSERT(st == QAWS_STATUS_OK && qaws_mesh_3d_validate(out) == QAWS_STATUS_OK && t.boundary_loops == 1 && okb && found == 4 && nb > 40
		&& rep.quality.min_angle > 20, msg);
done:
	qaws_mesh_3d_destroy(m);
	qaws_mesh_3d_destroy(out);
}

/* FIELD_ALIGNED: valence, angles, the edges' spread; (on the cube) corners and crease length */
static void rt_aligned_case(qaws_mesh_3d const* m, qaws_mesh_reference const* r, double L, char const* name, int cube)
{
	qaws_remesh_desc d;
	qaws_mesh_3d* out = NULL;
	qaws_mesh_report rep;
	qaws_mesh_topology t;
	qaws_status st;
	unsigned int nv = 0, i, v6 = 0, nin = 0, corners = 0;
	double foot = 0, off = 0, crease = 0;
	char msg[600];
	memset(&d, 0, sizeof d);
	memset(&rep, 0, sizeof rep);
	d.method = QAWS_REMESH_FIELD_ALIGNED;
	d.reference = r;
	d.target_length = (qaws_scalar)L;
	st = qaws_mesh_remesh(&out, &rep, m, &d);
	if (!out)
	{
		sprintf(msg, "FIELD_ALIGNED, %s: %s", name, qaws_status_to_string(st));
		TEST_ASSERT(0, msg);
		return;
	}
	nv = qaws_mesh_3d_get_vertex_count(out);
	{
		qaws_mesh_location* l = (qaws_mesh_location*)malloc(sizeof(qaws_mesh_location) * nv);
		qaws_scalar *x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv), *y = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
		unsigned int nfe = qaws_mesh_3d_get_feature_edge_count(out), *fe = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (nfe + 1));
		unsigned char* onf = (unsigned char*)calloc(nv + 1, 1);
		qaws_mesh_3d_get_feature_edges(fe, out);
		for (i = 0; i < 2 * nfe; i++)
			onf[fe[i]] = 1;
		qaws_mesh_3d_get_positions(x, out);
		qaws_mesh_3d_get_footpoints(l, out);
		qaws_mesh_reference_evaluate(y, NULL, NULL, r, l, nv);
		for (i = 0; i < nv; i++)
		{
			unsigned int ring[64], n = 0;
			foot = fmax(foot, fabs(x[3 * i] - y[3 * i]) + fabs(x[3 * i + 1] - y[3 * i + 1]) + fabs(x[3 * i + 2] - y[3 * i + 2]));
			if (cube)
				off = fmax(off, fabs(fmax(fabs(x[3 * i]), fmax(fabs(x[3 * i + 1]), fabs(x[3 * i + 2]))) - 1));
			else
				off = fmax(off, fabs(sqrt(x[3 * i] * x[3 * i] + x[3 * i + 1] * x[3 * i + 1] + x[3 * i + 2] * x[3 * i + 2]) - 1));
			corners += qaws_mesh_3d_is_feature_corner(out, i);
			if (onf[i])
				continue;
			qaws_mesh_3d_get_vertex_ring(ring, &n, out, 64, i);
			nin++;
			v6 += n == 6;
		}
		free(l); free(x); free(y); free(fe); free(onf);
	}
	if (cube)
	{
		unsigned int ne = qaws_mesh_3d_get_feature_edge_count(out), e;
		unsigned int* fe = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (ne + 1));
		qaws_mesh_3d_get_feature_edges(fe, out);
		for (e = 0; e < ne; e++)
		{
			qaws_vec3 a = qaws_mesh_3d_get_vertex(out, fe[2 * e]), b = qaws_mesh_3d_get_vertex(out, fe[2 * e + 1]);
			crease += sqrt((double)((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z)));
		}
		free(fe);
	}
	qaws_mesh_3d_get_topology(&t, out);
	sprintf(msg, "FIELD_ALIGNED, %s to %.2g: %u vertices, valid %d, closed %d, genus %u; footpoints %.2g, off %.2g; edges %.3g..%.3g mean %.3g, "
		"angles %.3g..%.3g, valence 6 at %.0f%% of %u smooth vertices; corners %u, crease length %.3g",
		name, L, nv, qaws_mesh_3d_validate(out) == QAWS_STATUS_OK, t.closed, t.genus, foot, off, (double)rep.quality.min_edge,
		(double)rep.quality.max_edge, (double)rep.quality.mean_edge, (double)rep.quality.min_angle, (double)rep.quality.max_angle,
		nin ? 100.0 * v6 / nin : 0, nin, corners, crease);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && qaws_mesh_3d_validate(out) == QAWS_STATUS_OK && t.closed && t.genus == 0 && t.components == 1 && foot < 1e-9
		&& off < 0.02 && fabs(rep.quality.mean_edge - L) < 0.25 * L && v6 > (cube ? 0.7 : 0.85) * nin && (!cube || (corners == 8 && fabs(crease - 24) < 0.05)), msg);
	qaws_mesh_3d_destroy(out);
}

static void test_field_aligned(void)
{
	static double const x[24] = { -1, -1, -1, 1, -1, -1, 1, 1, -1, -1, 1, -1, -1, -1, 1, 1, -1, 1, 1, 1, 1, -1, 1, 1 };
	static unsigned int const f[24] = { 0, 3, 2, 1, 4, 5, 6, 7, 0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6, 3, 0, 4, 7 };
	static unsigned int const fs[6] = { 4, 4, 4, 4, 4, 4 };
	qaws_scalar nrm[3 * 42];
	qaws_mesh_3d *m = rt_icosphere(nrm), *cube = rt_make(x, 8, f, fs, 6);
	qaws_mesh_reference_desc rd;
	qaws_mesh_reference *r = NULL, *rc = NULL;
	memset(&rd, 0, sizeof rd);
	rd.normals = nrm;
	qaws_mesh_reference_create(&r, m, &rd);
	if (r)
		rt_aligned_case(m, r, 0.15, "sphere", 0);
	if (cube)
	{
		qaws_mesh_3d_detect_features(cube, (qaws_scalar)(3.14159265358979 / 6), 0, 0);
		qaws_mesh_reference_create(&rc, cube, NULL);
	}
	if (rc)
		rt_aligned_case(cube, rc, 0.25, "cube with creases", 1);
	qaws_mesh_reference_destroy(r);
	qaws_mesh_reference_destroy(rc);
	qaws_mesh_3d_destroy(m);
	qaws_mesh_3d_destroy(cube);
}

static void test_adaptive_and_refusals(void)
{
	qaws_mesh_3d *m = rt_icosphere(NULL), *out = NULL, *other = rt_icosphere(NULL);
	qaws_mesh_reference *r = NULL, *ro = NULL;
	qaws_curvature_size_desc cd;
	qaws_size_field *sf = NULL, *so = NULL;
	qaws_mesh_report rep;
	qaws_status st = QAWS_STATUS_INTERNAL_ERROR, s1, s2;
	unsigned int nv, ne, i, nbump = 0, nfar = 0, *e;
	double lb = 0, lf = 0;
	qaws_scalar* p;
	char msg[400];
	nv = qaws_mesh_3d_get_vertex_count(m);
	p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	qaws_mesh_3d_get_positions(p, m);
	for (i = 0; i < nv; i++)
	{
		double d = p[3 * i] * p[3 * i] + p[3 * i + 1] * p[3 * i + 1] + (p[3 * i + 2] - 1) * (p[3 * i + 2] - 1), s = 1 + 0.6 * exp(-6 * d);
		p[3 * i] = (qaws_scalar)(p[3 * i] * s);
		p[3 * i + 1] = (qaws_scalar)(p[3 * i + 1] * s);
		p[3 * i + 2] = (qaws_scalar)(p[3 * i + 2] * s);
	}
	qaws_mesh_3d_set_positions(m, p);
	free(p);
	qaws_mesh_reference_create(&r, m, NULL);
	memset(&cd, 0, sizeof cd);
	cd.angle_tolerance = (qaws_scalar)0.2;
	cd.min_length = (qaws_scalar)0.02;
	cd.max_length = (qaws_scalar)0.4;
	cd.gradation = (qaws_scalar)0.4;
	if (r)
		qaws_size_field_create_curvature(&sf, r, &cd);
	out = sf ? rt_remesh(m, r, 0, sf, &rep, &st) : NULL;
	if (out)
	{
		nv = qaws_mesh_3d_get_vertex_count(out);
		ne = qaws_mesh_3d_get_edge_count(out);
		p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
		e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne);
		qaws_mesh_3d_get_positions(p, out);
		qaws_mesh_3d_get_edges(e, out);
		for (i = 0; i < ne; i++)
		{
			qaws_scalar const *a = p + 3 * e[2 * i], *b = p + 3 * e[2 * i + 1];
			double l = sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
			if (a[2] > 1.25 && b[2] > 1.25)
			{
				lb += l;
				nbump++;
			}
			else if (a[2] < -0.3 && b[2] < -0.3)
			{
				lf += l;
				nfar++;
			}
		}
		free(p);
		free(e);
	}
	sprintf(msg, "adaptive by curvature on a bumped sphere: valid %d, mean edge on the bump %.3g (%u), far from it %.3g (%u)",
		out && qaws_mesh_3d_validate(out) == QAWS_STATUS_OK, nbump ? lb / nbump : 0, nbump, nfar ? lf / nfar : 0, nfar);
	TEST_ASSERT(out && st == QAWS_STATUS_OK && qaws_mesh_3d_validate(out) == QAWS_STATUS_OK && nbump > 10 && nfar > 10 && lb / nbump < 0.6 * lf / nfar, msg);
	/* refusals: FIELD_ALIGNED with a metric, a field on another reference */
	{
		qaws_remesh_desc d;
		qaws_mesh_3d* bad = NULL;
		memset(&d, 0, sizeof d);
		d.method = QAWS_REMESH_FIELD_ALIGNED;
		d.metric = (qaws_metric_field const*)&d;   /* (any: refused before it is read) */
		s1 = qaws_mesh_remesh(&bad, NULL, m, &d);
		qaws_mesh_reference_create(&ro, other, NULL);
		if (ro)
			qaws_size_field_create_curvature(&so, ro, &cd);
		memset(&d, 0, sizeof d);
		d.reference = r;
		d.size = so;
		s2 = so ? qaws_mesh_remesh(&bad, NULL, m, &d) : QAWS_STATUS_OK;
		sprintf(msg, "remesh refusals: FIELD_ALIGNED with a metric %s, a size field on another reference %s", qaws_status_to_string(s1), qaws_status_to_string(s2));
		TEST_ASSERT(s1 == QAWS_STATUS_UNSUPPORTED_OPERATION && s2 == QAWS_STATUS_INVALID_ARGUMENT && !bad, msg);
	}
	qaws_size_field_destroy(sf);
	qaws_size_field_destroy(so);
	qaws_mesh_reference_destroy(r);
	qaws_mesh_reference_destroy(ro);
	qaws_mesh_3d_destroy(m);
	qaws_mesh_3d_destroy(other);
	qaws_mesh_3d_destroy(out);
}

int test_124_remesh_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 124: remeshing\n");
	test_sphere();
	test_cube();
	test_open_square();
	test_field_aligned();
	test_adaptive_and_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
