/*
 * Test 122: mesh references (qaws_mesh_reference)
 *
 *   - a PN triangle against Mathematica's from the paper's definitions
 *     (reference/122_pn_triangle.wls): position, first and second
 *     derivatives, the normal field and its derivatives, the curvatures;
 *   - an icosphere with its radial normals: the vertices interpolated, the
 *     patches meeting along their edges (positions and normals), PN far
 *     closer to the sphere than the flat triangles, curvatures near 1;
 *     PHONG (Boubekeur and Alexa) the same: vertices, edges met, about twice
 *     closer than the flat triangles (its shape factor 3/4 lifts an edge's
 *     middle by 3/16 |d|^2 where the sphere needs 1/8: past it by half the
 *     flat sag);
 *   - a cylinder: curvatures 1 / r and 0, the first direction round it;
 *   - a cube with its creases detected: every PN point on the cube;
 *   - SUBDIVISION against qaws_mesh_3d_subdivide (Loop, 2 levels, limit):
 *     on an icosahedron, an icosphere, a cube and the cube with its creases,
 *     each vertex's location the same limit point, every limit point on the
 *     reference (projected), its limit normal the reference's (smooth ones);
 *   - projection against brute force (LINEAR and PN), and with a hint the
 *     sheet the hint is on, not a nearer one;
 *   - walking: straight on a plane across triangles, and round a prism
 *     against the prism unrolled;
 *   - attributes at locations; refusals.
 */

#include "test_common.h"
#include "qaws_mesh_reference.h"
#include "qaws_mesh_subdivide.h"
#include "mesh/qaws_mesh_internal.h"
#include "reference/122_pn_triangle.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

typedef struct ts_mesh
{
	double* x;
	unsigned int nv, cap;
	unsigned int* f;
	unsigned int nf, capf;
} ts_mesh;

static unsigned int ts_vert(ts_mesh* m, double x, double y, double z)
{
	if (m->nv == m->cap)
	{
		m->cap = m->cap ? 2 * m->cap : 64;
		m->x = (double*)realloc(m->x, sizeof(double) * 3 * m->cap);
	}
	m->x[3 * m->nv] = x; m->x[3 * m->nv + 1] = y; m->x[3 * m->nv + 2] = z;
	return m->nv++;
}

static void ts_tri(ts_mesh* m, unsigned int a, unsigned int b, unsigned int c)
{
	if (m->nf == m->capf)
	{
		m->capf = m->capf ? 2 * m->capf : 64;
		m->f = (unsigned int*)realloc(m->f, sizeof(unsigned int) * 3 * m->capf);
	}
	m->f[3 * m->nf] = a; m->f[3 * m->nf + 1] = b; m->f[3 * m->nf + 2] = c;
	m->nf++;
}

static qaws_mesh_3d* ts_make(ts_mesh const* t)
{
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	qaws_scalar* p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * t->nv);
	unsigned int i;
	for (i = 0; i < 3 * t->nv; i++)
		p[i] = (qaws_scalar)t->x[i];
	memset(&d, 0, sizeof d);
	d.positions = p;
	d.vertex_count = t->nv;
	d.face_vertices = t->f;
	d.face_count = t->nf;
	qaws_mesh_3d_create(&m, &d);
	free(p);
	return m;
}

/* the midpoint vertex of edge (a, b), made once (pairs kept in a list) */
static unsigned int ts_mid(ts_mesh* m, unsigned int* keys, unsigned int* nk, unsigned int a, unsigned int b)
{
	unsigned int lo = a < b ? a : b, hi = a < b ? b : a, i, v;
	for (i = 0; i < *nk; i++)
		if (keys[3 * i] == lo && keys[3 * i + 1] == hi)
			return keys[3 * i + 2];
	v = ts_vert(m, (m->x[3 * a] + m->x[3 * b]) / 2, (m->x[3 * a + 1] + m->x[3 * b + 1]) / 2, (m->x[3 * a + 2] + m->x[3 * b + 2]) / 2);
	keys[3 * *nk] = lo; keys[3 * *nk + 1] = hi; keys[3 * *nk + 2] = v;
	(*nk)++;
	return v;
}

static void ts_icosphere(ts_mesh* m, unsigned int level)
{
	static unsigned int const f[60] = { 0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11, 1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8,
		3, 9, 4, 3, 4, 2, 3, 2, 6, 3, 6, 8, 3, 8, 9, 4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7, 9, 8, 1 };
	double g = (1 + sqrt(5.0)) / 2;
	double const v[36] = { -1, g, 0, 1, g, 0, -1, -g, 0, 1, -g, 0, 0, -1, g, 0, 1, g, 0, -1, -g, 0, 1, -g, g, 0, -1, g, 0, 1, -g, 0, -1, -g, 0, 1 };
	unsigned int i, l;
	memset(m, 0, sizeof *m);
	for (i = 0; i < 12; i++)
		ts_vert(m, v[3 * i], v[3 * i + 1], v[3 * i + 2]);
	for (i = 0; i < 20; i++)
		ts_tri(m, f[3 * i], f[3 * i + 1], f[3 * i + 2]);
	for (l = 0; l < level; l++)
	{
		unsigned int nf = m->nf, *old = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nf), *keys = (unsigned int*)malloc(sizeof(unsigned int) * 9 * nf), nk = 0;
		memcpy(old, m->f, sizeof(unsigned int) * 3 * nf);
		m->nf = 0;
		for (i = 0; i < nf; i++)
		{
			unsigned int a = old[3 * i], b = old[3 * i + 1], c = old[3 * i + 2];
			unsigned int ab = ts_mid(m, keys, &nk, a, b), bc = ts_mid(m, keys, &nk, b, c), ca = ts_mid(m, keys, &nk, c, a);
			ts_tri(m, a, ab, ca); ts_tri(m, ab, b, bc); ts_tri(m, ca, bc, c); ts_tri(m, ab, bc, ca);
		}
		free(old);
		free(keys);
	}
	for (i = 0; i < m->nv; i++)
	{
		double* p = m->x + 3 * i, r = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
		p[0] /= r; p[1] /= r; p[2] /= r;
	}
}

/* a prism of n sides, radius R (at its vertices), rings at heights 0..h */
static void ts_prism(ts_mesh* m, unsigned int n, unsigned int rings, double R, double h)
{
	unsigned int i, j;
	memset(m, 0, sizeof *m);
	for (j = 0; j < rings; j++)
		for (i = 0; i < n; i++)
			ts_vert(m, R * cos(2 * M_PI * i / n), R * sin(2 * M_PI * i / n), h * j / (rings - 1));
	for (j = 0; j + 1 < rings; j++)
		for (i = 0; i < n; i++)
		{
			unsigned int a = j * n + i, b = j * n + (i + 1) % n, c = (j + 1) * n + (i + 1) % n, d = (j + 1) * n + i;
			ts_tri(m, a, b, c);
			ts_tri(m, a, c, d);
		}
}

static void ts_free(ts_mesh* m)
{
	free(m->x);
	free(m->f);
}

static qaws_mesh_reference* ts_ref(qaws_mesh_3d const* m, qaws_interpolant kind, qaws_scalar const* normals)
{
	qaws_mesh_reference_desc d;
	qaws_mesh_reference* r = NULL;
	memset(&d, 0, sizeof d);
	d.interpolant = kind;
	d.normals = normals;
	qaws_mesh_reference_create(&r, m, &d);
	return r;
}

static unsigned int ts_rand(unsigned int* s)
{
	*s = *s * 1664525u + 1013904223u;
	return *s >> 8;
}

static double ts_u(unsigned int* s)
{
	return (ts_rand(s) & 0xFFFFFF) / (double)0x1000000;
}

static qaws_mesh_location ts_loc(unsigned int* s, unsigned int nt)
{
	qaws_mesh_location l;
	double a = ts_u(s), b = ts_u(s);
	if (a + b > 1)
	{
		a = 1 - a;
		b = 1 - b;
	}
	l.patch = ts_rand(s) % nt;
	l.a = (qaws_scalar)a;
	l.b = (qaws_scalar)b;
	return l;
}

static void test_pn_against_mathematica(void)
{
	mi_tri t;
	double const *P[3], *N[3];
	double emax = 0, scale = 0;
	unsigned int i, k, q;
	char msg[300];
	for (k = 0; k < 3; k++)
	{
		P[k] = g122_P + 3 * k;
		N[k] = g122_N + 3 * k;
	}
	mi_build(&t, QAWS_INTERPOLANT_PN, P, N, NULL);
	for (i = 0; i < 5; i++)
	{
		double v[27], k1, k2, d1[3], d2[3];
		double const* want = g122_rows + 29 * i;
		mi_eval(&t, g122_pts[2 * i], g122_pts[2 * i + 1], v, v + 3, v + 6, v + 9, v + 12, v + 15, v + 18, v + 21, v + 24);
		mi_curvature(v + 3, v + 6, v + 18, v + 21, v + 24, &k1, &k2, d1, d2);
		for (q = 0; q < 27; q++)
		{
			emax = fmax(emax, fabs(v[q] - want[q]));
			scale = fmax(scale, fabs(want[q]));
		}
		emax = fmax(emax, fmax(fabs(k1 - want[27]), fabs(k2 - want[28])));
	}
	sprintf(msg, "PN triangle against Mathematica at 5 points: position, derivatives to second order, normal field and its derivatives, curvatures: error %.2g (of %.3g)",
		emax, scale);
	TEST_ASSERT(emax < 1e-12 * fmax(1, scale), msg);
}

static void test_sphere(void)
{
	ts_mesh t;
	qaws_mesh_3d* m;
	qaws_mesh_reference *pn, *lin, *ph, *ng, *g1, *refs[4];
	qaws_scalar* nrm;
	unsigned int i, s = 5, nt, okv = 1;
	double vmax = 0, cmax = 0, nmax = 0, epn = 0, elin = 0, kmax = 0, eph = 0, eng = 0, eg1 = 0, gjump[4] = { 0, 0, 0, 0 };
	unsigned int r;
	char msg[400];
	ts_icosphere(&t, 2);
	m = ts_make(&t);
	nrm = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * t.nv);
	for (i = 0; i < 3 * t.nv; i++)
		nrm[i] = (qaws_scalar)t.x[i];
	pn = ts_ref(m, QAWS_INTERPOLANT_PN, nrm);
	lin = ts_ref(m, QAWS_INTERPOLANT_LINEAR, nrm);
	ph = ts_ref(m, QAWS_INTERPOLANT_PHONG, nrm);
	ng = ts_ref(m, QAWS_INTERPOLANT_NAGATA, nrm);
	g1 = ts_ref(m, QAWS_INTERPOLANT_G1, nrm);
	refs[0] = pn; refs[1] = ph; refs[2] = ng; refs[3] = g1;
	if (!m || !pn || !lin || !ph || !ng || !g1)
	{
		TEST_ASSERT(0, "sphere: the references");
		goto done;
	}
	nt = qaws_mesh_reference_get_triangle_count(pn);
	/* the vertices interpolated (PN, PHONG, NAGATA, G1) */
	for (r = 0; r < 4; r++)
	for (i = 0; i < t.nv; i++)
	{
		qaws_mesh_location l;
		qaws_scalar x[3];
		okv &= qaws_mesh_reference_get_vertex_location(&l, refs[r], i) == QAWS_STATUS_OK;
		qaws_mesh_reference_evaluate(x, NULL, NULL, refs[r], &l, 1);
		vmax = fmax(vmax, fabs(x[0] - t.x[3 * i]) + fabs(x[1] - t.x[3 * i + 1]) + fabs(x[2] - t.x[3 * i + 2]));
	}
	/* the patches meet along their edges: each triangle edge's points from both sides */
	for (r = 0; r < 4; r++)
	for (i = 0; i < 3 * nt; i++)
	{
		qaws_mesh_reference const* rr = refs[r];
		unsigned int nb = rr->ta[i], tt = i / 3, e = i % 3, j;
		if (nb == HM_NONE)
			continue;
		for (j = 1; j < 4; j++)
		{
			double lam = j / 4.0, b1[3] = { 0, 0, 0 }, b2[3] = { 0, 0, 0 }, x1[3], x2[3], n1[3], n2[3], a1[3], c1[3], a2[3], c2[3], g[2][3];
			unsigned int q, k;
			unsigned int u = nb / 3, ue = nb % 3;
			b1[e] = 1 - lam; b1[(e + 1) % 3] = lam;
			b2[ue] = lam; b2[(ue + 1) % 3] = 1 - lam;
			mi_eval(&rr->patch[tt], b1[1], b1[2], x1, a1, c1, NULL, NULL, NULL, n1, NULL, NULL);
			mi_eval(&rr->patch[u], b2[1], b2[2], x2, a2, c2, NULL, NULL, NULL, n2, NULL, NULL);
			/* the surface's own normals, xa x xb, from both sides */
			for (q = 0; q < 2; q++)
			{
				double const *A = q ? a2 : a1, *B = q ? c2 : c1;
				double l;
				g[q][0] = A[1] * B[2] - A[2] * B[1];
				g[q][1] = A[2] * B[0] - A[0] * B[2];
				g[q][2] = A[0] * B[1] - A[1] * B[0];
				l = sqrt(g[q][0] * g[q][0] + g[q][1] * g[q][1] + g[q][2] * g[q][2]);
				for (k = 0; k < 3; k++)
					g[q][k] /= l;
			}
			gjump[r] = fmax(gjump[r], fabs(g[0][0] - g[1][0]) + fabs(g[0][1] - g[1][1]) + fabs(g[0][2] - g[1][2]));
			cmax = fmax(cmax, fabs(x1[0] - x2[0]) + fabs(x1[1] - x2[1]) + fabs(x1[2] - x2[2]));
			nmax = fmax(nmax, fabs(n1[0] - n2[0]) + fabs(n1[1] - n2[1]) + fabs(n1[2] - n2[2]));
		}
	}
	sprintf(msg, "icosphere, PN, PHONG, NAGATA and G1: vertices interpolated (%.2g), patches meet along edges (positions %.2g, normals %.2g)", vmax, cmax, nmax);
	TEST_ASSERT(okv && vmax < TEST_BOUND(1e-12, 32) && cmax < TEST_BOUND(1e-12, 32) && nmax < TEST_BOUND(1e-12, 32), msg);
	sprintf(msg, "icosphere: the surfaces' own normals across edges jump by PN %.3g, PHONG %.3g, NAGATA %.3g, G1 %.2g (G1)", gjump[0], gjump[1],
		gjump[2], gjump[3]);
	printf("  %s\n", msg);
	TEST_ASSERT(gjump[3] < TEST_BOUND(1e-10, 32) && gjump[0] > 1e-4, msg);
	/* closer to the sphere than the flat triangles; curvatures near 1 */
	for (i = 0; i < 2000; i++)
	{
		qaws_mesh_location l = ts_loc(&s, nt);
		qaws_scalar x[3], y[3], k[8];
		qaws_mesh_reference_evaluate(x, NULL, k, pn, &l, 1);
		qaws_mesh_reference_evaluate(y, NULL, NULL, lin, &l, 1);
		{
			qaws_scalar z[3];
			qaws_mesh_reference_evaluate(z, NULL, NULL, ph, &l, 1);
			eph = fmax(eph, fabs(sqrt(z[0] * z[0] + z[1] * z[1] + z[2] * z[2]) - 1));
			qaws_mesh_reference_evaluate(z, NULL, NULL, ng, &l, 1);
			eng = fmax(eng, fabs(sqrt(z[0] * z[0] + z[1] * z[1] + z[2] * z[2]) - 1));
			qaws_mesh_reference_evaluate(z, NULL, NULL, g1, &l, 1);
			eg1 = fmax(eg1, fabs(sqrt(z[0] * z[0] + z[1] * z[1] + z[2] * z[2]) - 1));
		}
		epn = fmax(epn, fabs(sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]) - 1));
		elin = fmax(elin, fabs(sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]) - 1));
		kmax = fmax(kmax, fmax(fabs(k[0] - 1), fabs(k[1] - 1)));
	}
	sprintf(msg, "icosphere of %u triangles: off the sphere PN %.3g, PHONG %.3g, NAGATA %.3g, G1 %.3g, flat %.3g (%.0fx, %.1fx, %.0fx, %.0fx); PN curvatures off 1 by %.3g", nt, epn,
		eph, eng, eg1, elin, elin / epn, elin / eph, elin / eng, elin / eg1, kmax);
	printf("  %s\n", msg);
	TEST_ASSERT(epn < elin / 10 && eph < elin / 1.5 && eph > elin / 3 && eng < elin / 10 && eg1 < elin / 10 && kmax < 0.08, msg);
done:
	qaws_mesh_reference_destroy(pn);
	qaws_mesh_reference_destroy(lin);
	qaws_mesh_reference_destroy(ph);
	qaws_mesh_reference_destroy(ng);
	qaws_mesh_reference_destroy(g1);
	qaws_mesh_3d_destroy(m);
	free(nrm);
	ts_free(&t);
}

static void test_cylinder_and_cube(void)
{
	ts_mesh t;
	qaws_mesh_3d *m, *cube = NULL;
	qaws_mesh_reference *r = NULL, *rc = NULL;
	qaws_scalar* nrm;
	unsigned int i, s = 9, nt;
	double k1e = 0, k2e = 0, dax = 0, off = 0, onface = 0;
	char msg[400];
	ts_prism(&t, 48, 6, 1.5, 3);
	m = ts_make(&t);
	nrm = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * t.nv);
	for (i = 0; i < t.nv; i++)
	{
		nrm[3 * i] = (qaws_scalar)(t.x[3 * i] / 1.5);
		nrm[3 * i + 1] = (qaws_scalar)(t.x[3 * i + 1] / 1.5);
		nrm[3 * i + 2] = 0;
	}
	r = ts_ref(m, QAWS_INTERPOLANT_PN, nrm);
	if (!r)
	{
		TEST_ASSERT(0, "cylinder: the reference");
		goto done;
	}
	nt = qaws_mesh_reference_get_triangle_count(r);
	for (i = 0; i < 1000; i++)
	{
		qaws_mesh_location l = ts_loc(&s, nt);
		qaws_scalar k[8];
		qaws_mesh_reference_evaluate(NULL, NULL, k, r, &l, 1);
		k1e = fmax(k1e, fabs(k[0] - 1 / 1.5));
		k2e = fmax(k2e, fabs(k[1]));
		dax = fmax(dax, fabs(k[4]));   /* (d1's z: round the cylinder, not along it) */
	}
	sprintf(msg, "cylinder r 1.5: k1 off 1/r by %.3g, k2 off 0 by %.3g, the first direction's axial part %.3g", k1e, k2e, dax);
	TEST_ASSERT(k1e < 0.01 && k2e < 0.01 && dax < 0.05, msg);
	/* a cube, its creases detected: PN stays on the cube */
	{
		static qaws_scalar const x[24] = { -1, -1, -1, 1, -1, -1, 1, 1, -1, -1, 1, -1, -1, -1, 1, 1, -1, 1, 1, 1, 1, -1, 1, 1 };
		static unsigned int const f[24] = { 0, 3, 2, 1, 4, 5, 6, 7, 0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6, 3, 0, 4, 7 };
		static unsigned int const fs[6] = { 4, 4, 4, 4, 4, 4 };
		qaws_mesh_3d_desc d;
		memset(&d, 0, sizeof d);
		d.positions = x; d.vertex_count = 8; d.face_vertices = f; d.face_sizes = fs; d.face_count = 6;
		qaws_mesh_3d_create(&cube, &d);
		for (i = 0; cube && i < 6; i++)
			qaws_mesh_3d_split_face(NULL, cube, NULL, NULL, i);
		if (cube)
			qaws_mesh_3d_detect_features(cube, (qaws_scalar)(M_PI / 6), 0, 0);
		rc = cube ? ts_ref(cube, QAWS_INTERPOLANT_PN, NULL) : NULL;
	}
	if (rc)
	{
		nt = qaws_mesh_reference_get_triangle_count(rc);
		for (i = 0; i < 2000; i++)
		{
			qaws_mesh_location l = ts_loc(&s, nt);
			qaws_scalar x[3];
			double mx;
			qaws_mesh_reference_evaluate(x, NULL, NULL, rc, &l, 1);
			mx = fmax(fabs(x[0]), fmax(fabs(x[1]), fabs(x[2])));
			off = fmax(off, mx - 1);
			onface = fmax(onface, 1 - mx);
		}
	}
	sprintf(msg, "cube with creases, PN with split normals: points beyond the cube %.2g, inside its faces %.2g", off, onface);
	TEST_ASSERT(rc && off < 1e-12 && onface < 1e-12, msg);
done:
	qaws_mesh_reference_destroy(r);
	qaws_mesh_reference_destroy(rc);
	qaws_mesh_3d_destroy(m);
	qaws_mesh_3d_destroy(cube);
	free(nrm);
	ts_free(&t);
}

static void test_projection(void)
{
	ts_mesh t;
	qaws_mesh_3d* m;
	qaws_mesh_reference *pn, *lin;
	unsigned int i, s = 13, nt, kind;
	double emax[2] = { 0, 0 }, smax = 0;
	char msg[400];
	ts_icosphere(&t, 1);
	/* bumps: a non-convex shape */
	for (i = 0; i < t.nv; i++)
	{
		double* p = t.x + 3 * i, f = 1 + 0.15 * sin(3 * p[0]) * cos(2 * p[1]);
		p[0] *= f; p[1] *= f; p[2] *= f;
	}
	m = ts_make(&t);
	pn = ts_ref(m, QAWS_INTERPOLANT_PN, NULL);
	lin = ts_ref(m, QAWS_INTERPOLANT_LINEAR, NULL);
	if (!pn || !lin)
	{
		TEST_ASSERT(0, "projection: the references");
		goto done;
	}
	nt = qaws_mesh_reference_get_triangle_count(pn);
	for (kind = 0; kind < 2; kind++)
	{
		qaws_mesh_reference const* r = kind ? lin : pn;
		for (i = 0; i < 300; i++)
		{
			qaws_scalar p[3], d;
			qaws_mesh_location l;
			double q[3], best = HUGE_VAL, a, b;
			unsigned int tt;
			q[0] = 2.6 * ts_u(&s) - 1.3; q[1] = 2.6 * ts_u(&s) - 1.3; q[2] = 2.6 * ts_u(&s) - 1.3;
			p[0] = (qaws_scalar)q[0]; p[1] = (qaws_scalar)q[1]; p[2] = (qaws_scalar)q[2];
			qaws_mesh_reference_project(&l, &d, r, p, NULL, 1);
			for (tt = 0; tt < nt; tt++)
				best = fmin(best, mr_closest_on_patch(r, tt, q, &a, &b));
			emax[kind] = fmax(emax[kind], fabs(d - sqrt(best)));
			smax = fmax(smax, sqrt(best));
		}
	}
	sprintf(msg, "projection of 300 points against brute force: PN %.2g, LINEAR %.2g (distances up to %.3g)", emax[0], emax[1], smax);
	TEST_ASSERT(emax[0] < TEST_BOUND(1e-9, 32) && emax[1] < TEST_BOUND(1e-12, 32), msg);
	/* two sheets: the hint's, not the nearer one */
	{
		static qaws_scalar const x[24] = { 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0.1f, 1, 0, 0.1f, 1, 1, 0.1f, 0, 1, 0.1f };
		static unsigned int const f[12] = { 0, 1, 2, 0, 2, 3, 4, 5, 6, 4, 6, 7 };
		qaws_mesh_3d_desc d;
		qaws_mesh_3d* two = NULL;
		qaws_mesh_reference* r = NULL;
		qaws_scalar p[3] = { 0.5f, 0.4f, 0.06f }, d0 = 0, d1 = 0;
		qaws_mesh_location l0 = { 0 }, l1 = { 0 }, hint;
		memset(&d, 0, sizeof d);
		d.positions = x; d.vertex_count = 8; d.face_vertices = f; d.face_count = 4;
		qaws_mesh_3d_create(&two, &d);
		r = two ? ts_ref(two, QAWS_INTERPOLANT_PN, NULL) : NULL;
		hint.patch = 0; hint.a = (qaws_scalar)0.2; hint.b = (qaws_scalar)0.2;
		if (r)
		{
			qaws_mesh_reference_project(&l0, &d0, r, p, NULL, 1);
			qaws_mesh_reference_project(&l1, &d1, r, p, &hint, 1);
		}
		sprintf(msg, "two sheets 0.1 apart, a point 0.06 above the first: unhinted to triangle %u at %.4g, hinted on the first to triangle %u at %.4g",
			l0.patch, (double)d0, l1.patch, (double)d1);
		TEST_ASSERT(r && l0.patch >= 2 && fabs(d0 - 0.04) < 1e-6 && l1.patch < 2 && fabs(d1 - 0.06) < 1e-6, msg);
		qaws_mesh_reference_destroy(r);
		qaws_mesh_3d_destroy(two);
	}
done:
	qaws_mesh_reference_destroy(pn);
	qaws_mesh_reference_destroy(lin);
	qaws_mesh_3d_destroy(m);
	ts_free(&t);
}

static void test_walk(void)
{
	ts_mesh t;
	qaws_mesh_3d *m, *pm = NULL;
	qaws_mesh_reference *r = NULL, *rp = NULL;
	unsigned int i, j, n = 12;
	double e1 = 0, e2 = 0;
	char msg[400];
	/* a plane of 10 x 10 squares */
	memset(&t, 0, sizeof t);
	for (j = 0; j <= 10; j++)
		for (i = 0; i <= 10; i++)
			ts_vert(&t, i * 0.1, j * 0.1, 0);
	for (j = 0; j < 10; j++)
		for (i = 0; i < 10; i++)
		{
			unsigned int a = j * 11 + i;
			ts_tri(&t, a, a + 1, a + 12);
			ts_tri(&t, a, a + 12, a + 11);
		}
	m = ts_make(&t);
	r = m ? ts_ref(m, QAWS_INTERPOLANT_LINEAR, NULL) : NULL;
	if (r)
	{
		qaws_mesh_location from, to;
		qaws_scalar p[3] = { 0.43f, 0.51f, 0 }, x[3];
		for (i = 0; i < 20; i++)
		{
			double ang = 2 * M_PI * i / 20, len = 0.05 + 0.02 * i;
			qaws_scalar w[3];
			w[0] = (qaws_scalar)(len * cos(ang)); w[1] = (qaws_scalar)(len * sin(ang)); w[2] = 0;
			qaws_mesh_reference_project(&from, NULL, r, p, NULL, 1);
			qaws_mesh_reference_walk(&to, r, w, from);
			qaws_mesh_reference_evaluate(x, NULL, NULL, r, &to, 1);
			e1 = fmax(e1, fabs(x[0] - (p[0] + w[0])) + fabs(x[1] - (p[1] + w[1])));
		}
	}
	sprintf(msg, "walking straight on a plane across triangles: 20 directions, error %.2g", e1);
	TEST_ASSERT(r && e1 < TEST_BOUND(1e-12, 32), msg);
	/* round a 12-sided prism: against the prism unrolled */
	ts_free(&t);
	ts_prism(&t, n, 5, 1, 4);
	pm = ts_make(&t);
	rp = pm ? ts_ref(pm, QAWS_INTERPOLANT_LINEAR, NULL) : NULL;
	if (rp)
	{
		double side = 2 * sin(M_PI / n);
		for (i = 0; i < 12; i++)
		{
			double th = 0.15 + 0.07 * i, L = 2.5 + 0.3 * i, z0 = 1.0, p0 = 0.5 * side, pp, zz, want[3];
			double c0[3] = { 1, 0, 0 }, c1[3] = { cos(2 * M_PI / n), sin(2 * M_PI / n), 0 }, u[3];
			qaws_scalar start[3], w[3], x[3];
			qaws_mesh_location from, to;
			unsigned int k, face;
			for (k = 0; k < 3; k++)
				u[k] = (c1[k] - c0[k]) / side;
			start[0] = (qaws_scalar)((c0[0] + c1[0]) / 2); start[1] = (qaws_scalar)((c0[1] + c1[1]) / 2); start[2] = (qaws_scalar)z0;
			for (k = 0; k < 3; k++)
				w[k] = (qaws_scalar)(L * (cos(th) * u[k] + (k == 2 ? sin(th) * 0.4 : 0)));
			{
				/* the walked length is |w|: unrolled, its horizontal and vertical parts */
				double h = L * cos(th), v = L * sin(th) * 0.4, len = sqrt(h * h + v * v);
				qaws_mesh_reference_project(&from, NULL, rp, start, NULL, 1);
				qaws_mesh_reference_walk(&to, rp, w, from);
				qaws_mesh_reference_evaluate(x, NULL, NULL, rp, &to, 1);
				(void)len;
				pp = p0 + h;
				zz = z0 + v;
			}
			face = (unsigned int)floor(pp / side) % n;
			{
				double f = fmod(pp, side) / side, a0 = 2 * M_PI * face / n, a1 = 2 * M_PI * (face + 1) / n;
				want[0] = (1 - f) * cos(a0) + f * cos(a1);
				want[1] = (1 - f) * sin(a0) + f * sin(a1);
				want[2] = zz;
			}
			e2 = fmax(e2, fabs(x[0] - want[0]) + fabs(x[1] - want[1]) + fabs(x[2] - want[2]));
		}
	}
	sprintf(msg, "walking round a 12-sided prism, 12 helical paths over several faces: against the prism unrolled %.2g", e2);
	TEST_ASSERT(rp && e2 < TEST_BOUND(1e-9, 32), msg);
	qaws_mesh_reference_destroy(r);
	qaws_mesh_reference_destroy(rp);
	qaws_mesh_3d_destroy(m);
	qaws_mesh_3d_destroy(pm);
	ts_free(&t);
}

static void test_attributes_and_refusals(void)
{
	ts_mesh t;
	qaws_mesh_3d* m;
	qaws_mesh_reference* r = NULL;
	qaws_attribute_desc d;
	unsigned int i, s = 21, va = 0, fa = 0, nt;
	double ev = 0, ef = 0;
	qaws_scalar* vals;
	char msg[300];
	ts_icosphere(&t, 1);
	m = ts_make(&t);
	vals = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (t.nv > t.nf ? t.nv : t.nf));
	for (i = 0; i < t.nv; i++)
		vals[i] = (qaws_scalar)(t.x[3 * i] - 2 * t.x[3 * i + 1]);
	memset(&d, 0, sizeof d);
	d.name = "f"; d.element = QAWS_MESH_VERTEX; d.components = 1;
	qaws_mesh_3d_add_attribute(&va, m, vals, &d);
	for (i = 0; i < t.nf; i++)
		vals[i] = (qaws_scalar)i;
	d.name = "id"; d.element = QAWS_MESH_FACE; d.rule = QAWS_ATTRIBUTE_NEAREST;
	qaws_mesh_3d_add_attribute(&fa, m, vals, &d);
	r = ts_ref(m, QAWS_INTERPOLANT_LINEAR, NULL);
	nt = r ? qaws_mesh_reference_get_triangle_count(r) : 0;
	for (i = 0; r && i < 500; i++)
	{
		qaws_mesh_location l = ts_loc(&s, nt);
		qaws_scalar x[3], v, id;
		qaws_mesh_reference_evaluate(x, NULL, NULL, r, &l, 1);
		qaws_mesh_reference_evaluate_attribute(&v, r, &l, 1, va);
		qaws_mesh_reference_evaluate_attribute(&id, r, &l, 1, fa);
		ev = fmax(ev, fabs(v - (x[0] - 2 * x[1])));
		ef = fmax(ef, fabs(id - (double)l.patch));
	}
	sprintf(msg, "attributes at 500 locations of the flat reference: a linear vertex attribute %.2g, the face ids %.2g", ev, ef);
	TEST_ASSERT(r && ev < TEST_BOUND(1e-12, 32) && ef == 0, msg);
	{
		qaws_mesh_reference_desc dd;
		qaws_mesh_reference* bad = NULL;
		qaws_mesh_location l;
		qaws_status s1, s2;
		memset(&dd, 0, sizeof dd);
		dd.interpolant = QAWS_INTERPOLANT_MLS;
		s1 = qaws_mesh_reference_create(&bad, m, &dd);
		l.patch = 100000; l.a = 0; l.b = 0;
		s2 = r ? qaws_mesh_reference_evaluate(NULL, NULL, NULL, r, &l, 1) : QAWS_STATUS_OK;
		sprintf(msg, "refusals: a reserved interpolant %s, a location off the triangles %s", qaws_status_to_string(s1), qaws_status_to_string(s2));
		TEST_ASSERT(s1 == QAWS_STATUS_UNSUPPORTED_OPERATION && !bad && s2 == QAWS_STATUS_INVALID_ARGUMENT, msg);
	}
	qaws_mesh_reference_destroy(r);
	qaws_mesh_3d_destroy(m);
	free(vals);
	ts_free(&t);
}

/* SUBDIVISION: the reference's vertices against the Loop limit points (vertex i stays vertex i), every limit point projected
   on it, the limit normals (no creases) against its normals there */
static void ts_subdivision_case(qaws_mesh_3d const* m, char const* name, int smooth)
{
	qaws_mesh_reference* r = ts_ref(m, QAWS_INTERPOLANT_SUBDIVISION, NULL);
	qaws_mesh_3d* L = NULL;
	qaws_mesh_subdivide_desc sd;
	unsigned int i, nv = qaws_mesh_3d_get_vertex_count(m), nl, id;
	double ev = 0, ep = 0, en = 0;
	qaws_scalar* ln = NULL;
	char msg[300];
	memset(&sd, 0, sizeof sd);
	sd.scheme = QAWS_MESH_SUBDIVIDE_LOOP;
	sd.levels = 2;
	sd.limit = 1;
	sd.limit_normals = smooth;
	qaws_mesh_3d_subdivide(&L, m, &sd);
	if (!r || !L)
	{
		TEST_ASSERT(0, name);
		goto done;
	}
	nl = qaws_mesh_3d_get_vertex_count(L);
	if (smooth)
	{
		ln = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nl);
		qaws_mesh_3d_find_attribute(&id, L, "limit_normal");
		qaws_mesh_3d_get_attribute(ln, L, id);
	}
	for (i = 0; i < nv; i++)
	{
		qaws_mesh_location l;
		qaws_scalar x[3];
		qaws_vec3 y = qaws_mesh_3d_get_vertex(L, i);
		qaws_mesh_reference_get_vertex_location(&l, r, i);
		qaws_mesh_reference_evaluate(x, NULL, NULL, r, &l, 1);
		ev = fmax(ev, fabs(x[0] - y.x) + fabs(x[1] - y.y) + fabs(x[2] - y.z));
	}
	for (i = 0; i < nl; i++)
	{
		qaws_vec3 y = qaws_mesh_3d_get_vertex(L, i);
		qaws_scalar q[3], d, n[3];
		qaws_mesh_location l;
		q[0] = y.x; q[1] = y.y; q[2] = y.z;
		qaws_mesh_reference_project(&l, &d, r, q, NULL, 1);
		ep = fmax(ep, (double)d);
		if (smooth)
		{
			qaws_mesh_reference_evaluate(NULL, n, NULL, r, &l, 1);
			en = fmax(en, fabs(n[0] - ln[3 * i]) + fabs(n[1] - ln[3 * i + 1]) + fabs(n[2] - ln[3 * i + 2]));
		}
	}
	sprintf(msg, "SUBDIVISION, %s: vertices at their limit points %.2g, %u limit points off the reference %.2g, limit normals %.2g", name, ev, nl,
		ep, en);
	printf("  %s\n", msg);
	TEST_ASSERT(ev < TEST_BOUND(1e-12, 32) && ep < TEST_BOUND(1e-9, 32) && en < TEST_BOUND(1e-12, 32), msg);
done:
	free(ln);
	qaws_mesh_reference_destroy(r);
	qaws_mesh_3d_destroy(L);
}

static void test_subdivision(void)
{
	static qaws_scalar const x[24] = { -1, -1, -1, 1, -1, -1, 1, 1, -1, -1, 1, -1, -1, -1, 1, 1, -1, 1, 1, 1, 1, -1, 1, 1 };
	static unsigned int const f[24] = { 0, 3, 2, 1, 4, 5, 6, 7, 0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6, 3, 0, 4, 7 };
	static unsigned int const fs[6] = { 4, 4, 4, 4, 4, 4 };
	ts_mesh t;
	qaws_mesh_3d *m, *cube = NULL;
	qaws_mesh_3d_desc d;
	unsigned int i, lv;
	for (lv = 0; lv < 2; lv++)
	{
		ts_icosphere(&t, lv);
		m = ts_make(&t);
		ts_subdivision_case(m, lv ? "icosphere" : "icosahedron", 1);
		qaws_mesh_3d_destroy(m);
		ts_free(&t);
	}
	memset(&d, 0, sizeof d);
	d.positions = x; d.vertex_count = 8; d.face_vertices = f; d.face_sizes = fs; d.face_count = 6;
	qaws_mesh_3d_create(&cube, &d);
	for (i = 0; cube && i < 6; i++)
		qaws_mesh_3d_split_face(NULL, cube, NULL, NULL, i);
	ts_subdivision_case(cube, "cube", 1);
	if (cube)
		qaws_mesh_3d_detect_features(cube, (qaws_scalar)(M_PI / 6), 0, 0);
	ts_subdivision_case(cube, "cube with creases", 0);
	qaws_mesh_3d_destroy(cube);
}

int test_122_mesh_reference_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 122: mesh references\n");
	test_pn_against_mathematica();
	test_sphere();
	test_cylinder_and_cube();
	test_subdivision();
	test_projection();
	test_walk();
	test_attributes_and_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
