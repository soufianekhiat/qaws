/*
 * Test 125: cross fields and quadrangulation (INSTANT)
 *
 *   - cross fields: on an icosphere the indices sum to 4 chi = 8 quarter
 *     turns, on a torus to 0 (Poincare-Hopf); the field along the border of
 *     an open patch;
 *   - a torus quadrangulated: a valid closed mesh of genus 1, mostly quads,
 *     every vertex its footpoint evaluated and near the torus, edges near the
 *     target; all quads with `pure`;
 *   - a can with its creases, all quads: closed, every face outward, the
 *     crease edges on the rims;
 *   - a cube with a linear vertex attribute: transferred exactly;
 *   - refusals.
 */

#include "test_common.h"
#include "qaws_quadrangulate.h"
#include "qaws_mesh_field.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

static qaws_mesh_3d* qt_make(double const* x, unsigned int nv, unsigned int const* f, unsigned int const* fs, unsigned int nf)
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

static qaws_mesh_3d* qt_torus(unsigned int nu, unsigned int nv, double R, double r)
{
	double* x = (double*)malloc(sizeof(double) * 3 * nu * nv);
	unsigned int* f = (unsigned int*)malloc(sizeof(unsigned int) * 6 * nu * nv), i, j, k = 0;
	qaws_mesh_3d* m;
	for (j = 0; j < nv; j++)
		for (i = 0; i < nu; i++)
		{
			double u = 2 * M_PI * i / nu, v = 2 * M_PI * (j + 0.5 * (i & 1)) / nv;
			x[3 * (j * nu + i)] = (R + r * cos(v)) * cos(u);
			x[3 * (j * nu + i) + 1] = (R + r * cos(v)) * sin(u);
			x[3 * (j * nu + i) + 2] = r * sin(v);
		}
	for (j = 0; j < nv; j++)
		for (i = 0; i < nu; i++)
		{
			unsigned int a = j * nu + i, b = j * nu + (i + 1) % nu, c = ((j + 1) % nv) * nu + (i + 1) % nu, d = ((j + 1) % nv) * nu + i;
			f[k++] = a; f[k++] = b; f[k++] = c;
			f[k++] = a; f[k++] = c; f[k++] = d;
		}
	m = qt_make(x, nu * nv, f, NULL, 2 * nu * nv);
	free(x);
	free(f);
	return m;
}

/* the icosahedron split twice, on the unit sphere */
static qaws_mesh_3d* qt_sphere(void)
{
	static unsigned int const f0[60] = { 0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11, 1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8,
		3, 9, 4, 3, 4, 2, 3, 2, 6, 3, 6, 8, 3, 8, 9, 4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7, 9, 8, 1 };
	double g = (1 + sqrt(5.0)) / 2, x[3 * 162];
	double const v0[36] = { -1, g, 0, 1, g, 0, -1, -g, 0, 1, -g, 0, 0, -1, g, 0, 1, g, 0, -1, -g, 0, 1, -g, g, 0, -1, g, 0, 1, -g, 0, -1, -g, 0, 1 };
	unsigned int f[3 * 320], tmp[3 * 320], nv = 12, nf = 20, l, i, j, q, k;
	memcpy(x, v0, sizeof v0);
	memcpy(f, f0, sizeof f0);
	for (l = 0; l < 2; l++)
	{
		unsigned int keys[3 * 480], nk = 0, onf = nf;
		memcpy(tmp, f, sizeof(unsigned int) * 3 * nf);
		nf = 0;
		for (i = 0; i < onf; i++)
		{
			unsigned int e[3];
			for (j = 0; j < 3; j++)
			{
				unsigned int a = tmp[3 * i + j], b = tmp[3 * i + (j + 1) % 3], lo = a < b ? a : b, hi = a < b ? b : a;
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
				unsigned int a = tmp[3 * i], b = tmp[3 * i + 1], c = tmp[3 * i + 2];
				unsigned int t[12] = { a, e[0], e[2], e[0], b, e[1], e[2], e[1], c, e[0], e[1], e[2] };
				memcpy(f + 3 * nf, t, sizeof t);
				nf += 4;
			}
		}
	}
	for (i = 0; i < nv; i++)
	{
		double r = sqrt(x[3 * i] * x[3 * i] + x[3 * i + 1] * x[3 * i + 1] + x[3 * i + 2] * x[3 * i + 2]);
		for (k = 0; k < 3; k++)
			x[3 * i + k] /= r;
	}
	return qt_make(x, nv, f, NULL, nf);
}

static void test_fields(void)
{
	qaws_mesh_3d *s = qt_sphere(), *t = qt_torus(30, 12, 2, 0.75);
	unsigned int i, nfs = qaws_mesh_3d_get_face_count(s), nft = qaws_mesh_3d_get_face_count(t), nv = qaws_mesh_3d_get_vertex_count(s);
	int *is = (int*)malloc(sizeof(int) * nfs), *it = (int*)malloc(sizeof(int) * nft), ss = 0, st = 0;
	qaws_scalar* dir = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	double unit = 0, tang = 0;
	char msg[300];
	qaws_scalar* p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	qaws_mesh_3d_direction_field(dir, is, s, NULL);
	qaws_mesh_3d_direction_field(NULL, it, t, NULL);
	qaws_mesh_3d_get_vertex_normals(p, s);   /* (tangent to the normals the field uses) */
	for (i = 0; i < nfs; i++)
		ss += is[i];
	for (i = 0; i < nft; i++)
		st += it[i];
	for (i = 0; i < nv; i++)
	{
		double l = sqrt(dir[3 * i] * dir[3 * i] + dir[3 * i + 1] * dir[3 * i + 1] + dir[3 * i + 2] * dir[3 * i + 2]);
		unit = fmax(unit, fabs(l - 1));
		tang = fmax(tang, fabs(dir[3 * i] * p[3 * i] + dir[3 * i + 1] * p[3 * i + 1] + dir[3 * i + 2] * p[3 * i + 2]));
	}
	sprintf(msg, "cross fields: indices on an icosphere sum to %d quarter turns (4 chi = 8), on a torus to %d (0); unit %.2g, tangent %.2g",
		ss, st, unit, tang);
	TEST_ASSERT(ss == 8 && st == 0 && unit < 1e-9 && tang < 1e-9, msg);
	free(is); free(it); free(dir); free(p);
	qaws_mesh_3d_destroy(s);
	qaws_mesh_3d_destroy(t);
}

static qaws_mesh_3d* qt_quads(qaws_mesh_3d const* m, double L, int pure, qaws_mesh_report* rep, qaws_status* st)
{
	qaws_quadrangulate_desc d;
	qaws_mesh_3d* out = NULL;
	memset(&d, 0, sizeof d);
	d.target_length = (qaws_scalar)L;
	d.pure = pure;
	*st = qaws_mesh_quadrangulate(&out, rep, m, &d);
	return out;
}

static void test_torus(void)
{
	qaws_mesh_3d *t = qt_torus(30, 12, 2, 0.75), *q = NULL, *p = NULL;
	qaws_mesh_reference* r = NULL;
	qaws_mesh_report rep, rep2;
	qaws_mesh_topology top;
	qaws_status st = QAWS_STATUS_INTERNAL_ERROR, st2 = QAWS_STATUS_INTERNAL_ERROR;
	unsigned int nv = 0, i;
	double foot = 0, off = 0;
	char msg[500];
	double const L = test_heavy() ? 0.3 : 0.45;   /* (by default 0.45: the quadrangulations 2.25 times cheaper; the edge bound follows) */
	q = qt_quads(t, L, 0, &rep, &st);
	p = qt_quads(t, L, 1, &rep2, &st2);
	qaws_mesh_reference_create(&r, t, NULL);
	if (!q || !p || !r)
	{
		TEST_ASSERT(0, "torus: quadrangulated");
		goto done;
	}
	qaws_mesh_3d_get_topology(&top, q);
	nv = qaws_mesh_3d_get_vertex_count(q);
	{
		qaws_mesh_location* l = (qaws_mesh_location*)malloc(sizeof(qaws_mesh_location) * nv);
		qaws_scalar *x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv), *y = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
		qaws_mesh_3d_get_positions(x, q);
		qaws_mesh_3d_get_footpoints(l, q);
		qaws_mesh_reference_evaluate(y, NULL, NULL, r, l, nv);
		for (i = 0; i < nv; i++)
		{
			double rho = sqrt(x[3 * i] * x[3 * i] + x[3 * i + 1] * x[3 * i + 1]), dz = x[3 * i + 2];
			foot = fmax(foot, fabs(x[3 * i] - y[3 * i]) + fabs(x[3 * i + 1] - y[3 * i + 1]) + fabs(x[3 * i + 2] - y[3 * i + 2]));
			off = fmax(off, fabs(sqrt((rho - 2) * (rho - 2) + dz * dz) - 0.75));
		}
		free(l); free(x); free(y);
	}
	sprintf(msg, "torus quadrangulated to %.2f: valid %d, closed %d, genus %u; %u quads of %u faces; footpoints %.2g, off the torus %.2g; mean edge %.3g; pure: %u quads of %u, valid %d",
		L, qaws_mesh_3d_validate(q) == QAWS_STATUS_OK, top.closed, top.genus, rep.quality.quads, qaws_mesh_3d_get_face_count(q), foot, off,
		(double)rep.quality.mean_edge, rep2.quality.quads, qaws_mesh_3d_get_face_count(p), qaws_mesh_3d_validate(p) == QAWS_STATUS_OK);
	TEST_ASSERT(st == QAWS_STATUS_OK && st2 == QAWS_STATUS_OK && qaws_mesh_3d_validate(q) == QAWS_STATUS_OK && top.closed && top.genus == 1
		&& rep.quality.quads > 0.85 * qaws_mesh_3d_get_face_count(q) && foot < 1e-9 && off < 0.02 && fabs(rep.quality.mean_edge - L) < 0.08 * L / 0.3
		&& rep2.quality.quads == qaws_mesh_3d_get_face_count(p) && qaws_mesh_3d_validate(p) == QAWS_STATUS_OK, msg);
done:
	qaws_mesh_reference_destroy(r);
	qaws_mesh_3d_destroy(t);
	qaws_mesh_3d_destroy(q);
	qaws_mesh_3d_destroy(p);
}

static void test_can_and_cube(void)
{
	unsigned int n = 8, i, k = 0, nf, f, c0, bad = 0, nfe, onrim = 1;
	double x[3 * 18];
	unsigned int fv[3 * 32];
	qaws_mesh_3d *can, *q = NULL, *cube = NULL, *qc = NULL;
	qaws_mesh_report rep;
	qaws_mesh_topology top;
	qaws_status st = QAWS_STATUS_INTERNAL_ERROR, st2 = QAWS_STATUS_INTERNAL_ERROR;
	char msg[500];
	double const L = test_heavy() ? 0.16 : 0.24, Lc = test_heavy() ? 0.3 : 0.5;   /* (by default coarser, more than 40 crease edges still) */
	for (i = 0; i < n; i++)
	{
		double a = 2 * M_PI * i / n;
		x[3 * i] = cos(a); x[3 * i + 1] = sin(a); x[3 * i + 2] = -0.8;
		x[3 * (n + i)] = cos(a); x[3 * (n + i) + 1] = sin(a); x[3 * (n + i) + 2] = 0.8;
	}
	x[6 * n] = 0; x[6 * n + 1] = 0; x[6 * n + 2] = -0.8;
	x[6 * n + 3] = 0; x[6 * n + 4] = 0; x[6 * n + 5] = 0.8;
	for (i = 0; i < n; i++)
	{
		unsigned int a = i, b = (i + 1) % n;
		fv[k++] = a; fv[k++] = b; fv[k++] = n + b;
		fv[k++] = a; fv[k++] = n + b; fv[k++] = n + a;
		fv[k++] = 2 * n; fv[k++] = b; fv[k++] = a;
		fv[k++] = 2 * n + 1; fv[k++] = n + a; fv[k++] = n + b;
	}
	can = qt_make(x, 2 * n + 2, fv, NULL, 4 * n);
	qaws_mesh_3d_detect_features(can, (qaws_scalar)(M_PI / 3), 0, 0);   /* (the rims at 90 degrees; the sides meet at 45) */
	q = qt_quads(can, L, 1, &rep, &st);
	if (!q)
	{
		TEST_ASSERT(0, "can: quadrangulated");
		goto done;
	}
	qaws_mesh_3d_get_topology(&top, q);
	nf = qaws_mesh_3d_get_face_count(q);
	nfe = qaws_mesh_3d_get_feature_edge_count(q);
	{
		unsigned int nvq = qaws_mesh_3d_get_vertex_count(q), *ff = (unsigned int*)malloc(sizeof(unsigned int) * qaws_mesh_3d_get_corner_count(q));
		unsigned int *fs = (unsigned int*)malloc(sizeof(unsigned int) * nf), *fe = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (nfe + 1));
		qaws_scalar* p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nvq);
		qaws_mesh_3d_get_positions(p, q);
		qaws_mesh_3d_get_faces(ff, fs, q);
		for (f = 0, c0 = 0; f < nf; c0 += fs[f], f++)
		{
			double cen[3] = { 0, 0, 0 }, nn[3] = { 0, 0, 0 }, out[3];
			unsigned int j;
			for (j = 0; j < fs[f]; j++)
			{
				qaws_scalar const *a = p + 3 * ff[c0 + j], *b = p + 3 * ff[c0 + (j + 1) % fs[f]];
				nn[0] += (a[1] - b[1]) * (a[2] + b[2]); nn[1] += (a[2] - b[2]) * (a[0] + b[0]); nn[2] += (a[0] - b[0]) * (a[1] + b[1]);
				cen[0] += a[0] / fs[f]; cen[1] += a[1] / fs[f]; cen[2] += a[2] / fs[f];
			}
			if (fabs(cen[2]) > 0.79) { out[0] = 0; out[1] = 0; out[2] = cen[2] > 0 ? 1 : -1; }
			else { out[0] = cen[0]; out[1] = cen[1]; out[2] = 0; }
			bad += nn[0] * out[0] + nn[1] * out[1] + nn[2] * out[2] <= 0;
		}
		qaws_mesh_3d_get_feature_edges(fe, q);
		for (i = 0; i < 2 * nfe; i++)
		{
			onrim &= fabs(fabs(p[3 * fe[i] + 2]) - 0.8) < 1e-6;
		}
		free(ff); free(fs); free(fe); free(p);
	}
	sprintf(msg, "can with creases, all quads: %u quads of %u faces, valid %d, closed %d, faces turned in %u, %u crease edges all on the rims %u",
		rep.quality.quads, nf, qaws_mesh_3d_validate(q) == QAWS_STATUS_OK, top.closed, bad, nfe, onrim);
	TEST_ASSERT(st == QAWS_STATUS_OK && rep.quality.quads == nf && qaws_mesh_3d_validate(q) == QAWS_STATUS_OK && top.closed && bad == 0 && nfe > 40 && onrim, msg);
	/* a cube's linear attribute, transferred */
	{
		static double const cx[24] = { -1, -1, -1, 1, -1, -1, 1, 1, -1, -1, 1, -1, -1, -1, 1, 1, -1, 1, 1, 1, 1, -1, 1, 1 };
		static unsigned int const cf[24] = { 0, 3, 2, 1, 4, 5, 6, 7, 0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6, 3, 0, 4, 7 };
		static unsigned int const cs[6] = { 4, 4, 4, 4, 4, 4 };
		qaws_scalar vals[8];
		qaws_attribute_desc ad;
		unsigned int a, nvq;
		double ev = 0;
		cube = qt_make(cx, 8, cf, cs, 6);
		for (i = 0; i < 8; i++)
			vals[i] = (qaws_scalar)(cx[3 * i] - 2 * cx[3 * i + 1] + 0.5 * cx[3 * i + 2]);
		memset(&ad, 0, sizeof ad);
		ad.name = "f"; ad.element = QAWS_MESH_VERTEX; ad.components = 1;
		qaws_mesh_3d_add_attribute(&a, cube, vals, &ad);
		qaws_mesh_3d_detect_features(cube, (qaws_scalar)(M_PI / 6), 0, 0);
		qc = qt_quads(cube, Lc, 1, &rep, &st2);
		if (qc && qaws_mesh_3d_find_attribute(&a, qc, "f") == QAWS_STATUS_OK)
		{
			qaws_scalar *p, *v;
			nvq = qaws_mesh_3d_get_vertex_count(qc);
			p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nvq);
			v = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nvq);
			qaws_mesh_3d_get_positions(p, qc);
			qaws_mesh_3d_get_attribute(v, qc, a);
			for (i = 0; i < nvq; i++)
				ev = fmax(ev, fabs(v[i] - (p[3 * i] - 2 * p[3 * i + 1] + 0.5 * p[3 * i + 2])));
			free(p);
			free(v);
		}
		sprintf(msg, "cube with a linear attribute, all quads: %u quads, attribute error %.2g", qc ? rep.quality.quads : 0, ev);
		TEST_ASSERT(st2 == QAWS_STATUS_OK && qc && ev < 1e-9, msg);
	}
done:
	qaws_mesh_3d_destroy(can);
	qaws_mesh_3d_destroy(q);
	qaws_mesh_3d_destroy(cube);
	qaws_mesh_3d_destroy(qc);
}

static void test_refusals(void)
{
	qaws_mesh_3d* t = qt_torus(12, 6, 2, 0.75);
	qaws_quadrangulate_desc d;
	qaws_direction_field_desc dd;
	qaws_mesh_3d* out = NULL;
	qaws_status s1, s2;
	char msg[200];
	memset(&d, 0, sizeof d);
	d.method = QAWS_QUADS_LAYOUT;
	s1 = qaws_mesh_quadrangulate(&out, NULL, t, &d);
	memset(&dd, 0, sizeof dd);
	dd.symmetry = 6;
	s2 = qaws_mesh_3d_direction_field(NULL, NULL, t, &dd);
	sprintf(msg, "refusals: LAYOUT %s, a 6-symmetric field %s", qaws_status_to_string(s1), qaws_status_to_string(s2));
	TEST_ASSERT(s1 == QAWS_STATUS_UNSUPPORTED_OPERATION && !out && s2 == QAWS_STATUS_UNSUPPORTED_OPERATION, msg);
	qaws_mesh_3d_destroy(t);
}

int test_125_quadrangulate_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 125: cross fields and quadrangulation\n");
	test_fields();
	test_torus();
	test_can_and_cube();
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
