/*
 * Test 159: spatial queries on meshes (qaws_mesh_bvh_*)
 *
 *   - a closed sphere of triangles: closest points and rays against a
 *     brute-force loop over every triangle (the same distances, the hit's
 *     weights giving its point); the winding number against the exact sum
 *     (1 inside, 0 outside, the dipole expansion within 1e-3); signed
 *     distances negative inside; the same sphere with a cap removed still
 *     about 1 at its centre; refit after scaling equal to a new build; an
 *     executor running the chunks backwards gives the same results;
 *   - watertight rays: a plane of triangles shot through every vertex and
 *     every edge's midpoint, no ray slipping through;
 *   - a 2D mesh: points located in their triangle (distance 0, weights
 *     giving the point); its border: winding 1 inside, 0 outside, rays from
 *     inside hitting it, signed distances;
 *   - a volume of hexahedra: points located in their cell, a point outside
 *     at its distance from the block;
 *   - refusals: a bad dimension, tetrahedra in 2D, an index out of range,
 *     a position not finite, rays and winding numbers where they do not
 *     apply.
 */

#include "test_common.h"
#include "qaws_mesh_bvh.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define BT_PI 3.14159265358979323846

static unsigned int g_seed = 12345u;
static double bt_rand(void)
{
	g_seed = g_seed * 1664525u + 1013904223u;
	return (g_seed >> 8) / 16777216.0;
}

/* a UV sphere of radius 1, outward, na around and nb from pole to pole; positions and triangles allocated */
static void bt_sphere(unsigned int na, unsigned int nb, qaws_scalar** ox, unsigned int* onv, unsigned int** ot, unsigned int* ont)
{
	unsigned int nv = 2 + na * (nb - 1), nt = 2 * na * (nb - 1), i, j, t = 0;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	unsigned int* f = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nt);
	x[0] = 0; x[1] = 0; x[2] = 1;
	x[3] = 0; x[4] = 0; x[5] = -1;
	for (j = 1; j < nb; j++)
		for (i = 0; i < na; i++)
		{
			double th = BT_PI * j / nb, ph = 2 * BT_PI * i / na;
			unsigned int v = 2 + (j - 1) * na + i;
			x[3 * v] = (qaws_scalar)(sin(th) * cos(ph));
			x[3 * v + 1] = (qaws_scalar)(sin(th) * sin(ph));
			x[3 * v + 2] = (qaws_scalar)cos(th);
		}
	for (i = 0; i < na; i++)
	{
		unsigned int a = 2 + i, b = 2 + (i + 1) % na;
		f[t++] = 0; f[t++] = a; f[t++] = b;
		a = 2 + (nb - 2) * na + i;
		b = 2 + (nb - 2) * na + (i + 1) % na;
		f[t++] = 1; f[t++] = b; f[t++] = a;
	}
	for (j = 1; j + 1 < nb; j++)
		for (i = 0; i < na; i++)
		{
			unsigned int a = 2 + (j - 1) * na + i, b = 2 + (j - 1) * na + (i + 1) % na, c = b + na, d = a + na;
			f[t++] = a; f[t++] = d; f[t++] = c;
			f[t++] = a; f[t++] = c; f[t++] = b;
		}
	*ox = x;
	*onv = nv;
	*ot = f;
	*ont = t / 3;
}

/* brute force: squared distance from p to triangle abc by dense barycentric search refined (Ericson in the library; here a plain projection check) */
static double bt_tri_d2(double const* p, double const* a, double const* b, double const* c)
{
	/* minimize over the triangle by projecting onto the plane, else the three edges */
	double e1[3], e2[3], n[3], w[3], d, best = HUGE_VAL;
	unsigned int k, e;
	double const* v[3] = { a, b, c };
	for (k = 0; k < 3; k++)
	{
		e1[k] = b[k] - a[k];
		e2[k] = c[k] - a[k];
	}
	n[0] = e1[1] * e2[2] - e1[2] * e2[1];
	n[1] = e1[2] * e2[0] - e1[0] * e2[2];
	n[2] = e1[0] * e2[1] - e1[1] * e2[0];
	{
		double nn = n[0] * n[0] + n[1] * n[1] + n[2] * n[2], q[3], s, l1, l2, l0, area[3];
		if (nn > 0)
		{
			s = ((p[0] - a[0]) * n[0] + (p[1] - a[1]) * n[1] + (p[2] - a[2]) * n[2]) / nn;
			for (k = 0; k < 3; k++)
				q[k] = p[k] - s * n[k];
			for (e = 0; e < 3; e++)
			{
				double const *u = v[(e + 1) % 3], *t = v[(e + 2) % 3];
				double x1[3], x2[3], cr[3];
				for (k = 0; k < 3; k++)
				{
					x1[k] = u[k] - q[k];
					x2[k] = t[k] - q[k];
				}
				cr[0] = x1[1] * x2[2] - x1[2] * x2[1];
				cr[1] = x1[2] * x2[0] - x1[0] * x2[2];
				cr[2] = x1[0] * x2[1] - x1[1] * x2[0];
				area[e] = (cr[0] * n[0] + cr[1] * n[1] + cr[2] * n[2]) / nn;
			}
			l0 = area[0];
			l1 = area[1];
			l2 = area[2];
			if (l0 >= 0 && l1 >= 0 && l2 >= 0)
				return s * s * nn;
		}
	}
	for (e = 0; e < 3; e++)
	{
		double const *u = v[e], *t = v[(e + 1) % 3];
		double ed[3], l = 0, tt = 0;
		for (k = 0; k < 3; k++)
		{
			ed[k] = t[k] - u[k];
			l += ed[k] * ed[k];
			tt += (p[k] - u[k]) * ed[k];
		}
		tt = l > 0 ? tt / l : 0;
		tt = tt < 0 ? 0 : tt > 1 ? 1 : tt;
		d = 0;
		for (k = 0; k < 3; k++)
			d += (u[k] + tt * ed[k] - p[k]) * (u[k] + tt * ed[k] - p[k]);
		if (d < best)
			best = d;
		(void)w;
	}
	return best;
}

static double bt_solid(double const* p, double const* a, double const* b, double const* c)
{
	double A[3], B[3], C[3], la, lb, lc, num, den;
	unsigned int k;
	for (k = 0; k < 3; k++)
	{
		A[k] = a[k] - p[k];
		B[k] = b[k] - p[k];
		C[k] = c[k] - p[k];
	}
	la = sqrt(A[0] * A[0] + A[1] * A[1] + A[2] * A[2]);
	lb = sqrt(B[0] * B[0] + B[1] * B[1] + B[2] * B[2]);
	lc = sqrt(C[0] * C[0] + C[1] * C[1] + C[2] * C[2]);
	num = A[0] * (B[1] * C[2] - B[2] * C[1]) + A[1] * (B[2] * C[0] - B[0] * C[2]) + A[2] * (B[0] * C[1] - B[1] * C[0]);
	den = la * lb * lc + (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]) * lc + (B[0] * C[0] + B[1] * C[1] + B[2] * C[2]) * la
		+ (C[0] * A[0] + C[1] * A[1] + C[2] * A[2]) * lb;
	return 2 * atan2(num, den) / (4 * BT_PI);
}

static void bt_get(double* o, qaws_scalar const* x, unsigned int v)
{
	o[0] = x[3 * v];
	o[1] = x[3 * v + 1];
	o[2] = x[3 * v + 2];
}

/* runs the chunks backwards, one at a time */
static void bt_backwards(void* user, unsigned int count, qaws_batch_task_fn task, void* ctx)
{
	(void)user;
	while (count--)
		task(ctx, count, count + 1);
}

static double bt_tol(void)
{
	return sizeof(qaws_scalar) == 4 ? 2e-5 : 1e-11;
}

static void test_sphere(void)
{
	qaws_scalar *x, *pts, *dir, *wn, *wn6, *sd, *x2;
	unsigned int nv, nt, *f, i, t, n = 400, bad_d = 0, bad_w = 0, bad_r = 0, bad_wn = 0, bad_sd = 0, bad_ex = 0, bad_rf = 0;
	double worst_wn = 0, worst_wn6 = 0;
	qaws_mesh_bvh_desc d;
	qaws_mesh_bvh *b = NULL, *b2 = NULL;
	qaws_mesh_hit *h, *h2;
	qaws_batch_executor ex = { bt_backwards, NULL };
	bt_sphere(48, 24, &x, &nv, &f, &nt);
	memset(&d, 0, sizeof d);
	d.positions = x;
	d.vertex_count = nv;
	d.dimension = 3;
	d.simplices = f;
	d.simplex_size = 3;
	d.simplex_count = nt;
	TEST_ASSERT(qaws_mesh_bvh_create(&b, &d) == QAWS_STATUS_OK && b, "sphere: built");
	pts = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n);
	dir = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n);
	wn = (qaws_scalar*)malloc(sizeof(qaws_scalar) * n);
	wn6 = (qaws_scalar*)malloc(sizeof(qaws_scalar) * n);
	sd = (qaws_scalar*)malloc(sizeof(qaws_scalar) * n);
	h = (qaws_mesh_hit*)malloc(sizeof(qaws_mesh_hit) * n);
	h2 = (qaws_mesh_hit*)malloc(sizeof(qaws_mesh_hit) * n);
	x2 = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	for (i = 0; i < 3 * n; i++)
	{
		pts[i] = (qaws_scalar)(4 * bt_rand() - 2);
		dir[i] = (qaws_scalar)(2 * bt_rand() - 1);
	}
	/* closest against brute force */
	TEST_ASSERT(qaws_mesh_bvh_closest(h, b, pts, n, 0, NULL) == QAWS_STATUS_OK, "sphere: closest ran");
	for (i = 0; i < n; i++)
	{
		double p[3] = { pts[3 * i], pts[3 * i + 1], pts[3 * i + 2] }, best = HUGE_VAL, q[3] = { 0, 0, 0 };
		unsigned int k, j;
		for (t = 0; t < nt; t++)
		{
			double a[3], bb[3], c[3], dd;
			bt_get(a, x, f[3 * t]);
			bt_get(bb, x, f[3 * t + 1]);
			bt_get(c, x, f[3 * t + 2]);
			dd = bt_tri_d2(p, a, bb, c);
			if (dd < best)
				best = dd;
		}
		if (fabs(sqrt(best) - h[i].distance) > bt_tol() * (1 + sqrt(best)))
		{
			bad_d++;
		}
		for (j = 0; j < 3; j++)
			for (k = 0; k < 3; k++)
				q[k] += h[i].weights[j] * x[3 * h[i].vertices[j] + k];
		for (k = 0; k < 3; k++)
			if (fabs(q[k] - h[i].point[k]) > bt_tol() * 4)
				bad_w++;
		if (h[i].element != h[i].simplex || h[i].vertices[3] != QAWS_MESH_NONE)
			bad_w++;
	}
	TEST_ASSERT(bad_d == 0, "sphere: every closest distance equals the brute force's");
	TEST_ASSERT(bad_w == 0, "sphere: every hit's weights give its point");
	/* rays against brute force (Moller-Trumbore in double, both sides) */
	TEST_ASSERT(qaws_mesh_bvh_ray(h, b, pts, dir, n, 0, 0, NULL) == QAWS_STATUS_OK, "sphere: rays ran");
	for (i = 0; i < n; i++)
	{
		double o[3] = { pts[3 * i], pts[3 * i + 1], pts[3 * i + 2] }, dv[3] = { dir[3 * i], dir[3 * i + 1], dir[3 * i + 2] }, best = HUGE_VAL;
		for (t = 0; t < nt; t++)
		{
			double a[3], bb[3], c[3], e1[3], e2[3], pv[3], tv[3], qv[3], det, u, v, tt;
			unsigned int k;
			bt_get(a, x, f[3 * t]);
			bt_get(bb, x, f[3 * t + 1]);
			bt_get(c, x, f[3 * t + 2]);
			for (k = 0; k < 3; k++)
			{
				e1[k] = bb[k] - a[k];
				e2[k] = c[k] - a[k];
				tv[k] = o[k] - a[k];
			}
			pv[0] = dv[1] * e2[2] - dv[2] * e2[1];
			pv[1] = dv[2] * e2[0] - dv[0] * e2[2];
			pv[2] = dv[0] * e2[1] - dv[1] * e2[0];
			det = e1[0] * pv[0] + e1[1] * pv[1] + e1[2] * pv[2];
			if (det == 0)
				continue;
			u = (tv[0] * pv[0] + tv[1] * pv[1] + tv[2] * pv[2]) / det;
			qv[0] = tv[1] * e1[2] - tv[2] * e1[1];
			qv[1] = tv[2] * e1[0] - tv[0] * e1[2];
			qv[2] = tv[0] * e1[1] - tv[1] * e1[0];
			v = (dv[0] * qv[0] + dv[1] * qv[1] + dv[2] * qv[2]) / det;
			tt = (e2[0] * qv[0] + e2[1] * qv[1] + e2[2] * qv[2]) / det;
			if (u >= 0 && v >= 0 && u + v <= 1 && tt >= 0 && tt < best)
				best = tt;
		}
		if ((best == HUGE_VAL) != (h[i].element == QAWS_MESH_NONE) || (best != HUGE_VAL && fabs(best - h[i].distance) > bt_tol() * 10 * (1 + best)))
			bad_r++;
	}
	TEST_ASSERT(bad_r == 0, "sphere: every ray's first hit equals the brute force's");
	/* winding numbers against the exact sum, at the default accuracy and at 6 */
	TEST_ASSERT(qaws_mesh_bvh_winding_number(wn, b, pts, n, NULL) == QAWS_STATUS_OK, "sphere: winding ran");
	d.winding_accuracy = 6;
	TEST_ASSERT(qaws_mesh_bvh_create(&b2, &d) == QAWS_STATUS_OK && qaws_mesh_bvh_winding_number(wn6, b2, pts, n, NULL) == QAWS_STATUS_OK, "sphere: winding at accuracy 6 ran");
	d.winding_accuracy = 0;
	qaws_mesh_bvh_destroy(b2);
	b2 = NULL;
	TEST_ASSERT(qaws_mesh_bvh_signed_distance(sd, h2, b, pts, n, NULL) == QAWS_STATUS_OK, "sphere: signed distances ran");
	TEST_ASSERT(qaws_mesh_bvh_closest(h, b, pts, n, 0, NULL) == QAWS_STATUS_OK, "sphere: closest again");
	for (i = 0; i < n; i++)
	{
		double p[3] = { pts[3 * i], pts[3 * i + 1], pts[3 * i + 2] }, w = 0, r = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
		for (t = 0; t < nt; t++)
		{
			double a[3], bb[3], c[3];
			bt_get(a, x, f[3 * t]);
			bt_get(bb, x, f[3 * t + 1]);
			bt_get(c, x, f[3 * t + 2]);
			w += bt_solid(p, a, bb, c);
		}
		if (fabs(w - wn[i]) > worst_wn)
			worst_wn = fabs(w - wn[i]);
		if (fabs(w - wn6[i]) > worst_wn6)
			worst_wn6 = fabs(w - wn6[i]);
		if ((r < 0.95 && fabs(wn[i] - 1) > 5e-3) || (r > 1.05 && fabs(wn[i]) > 5e-3))
			bad_wn++;
		if (fabs(fabs(sd[i]) - h[i].distance) > bt_tol() || (r < 0.95 && !(sd[i] < 0)) || (r > 1.05 && !(sd[i] > 0)) || h2[i].simplex != h[i].simplex)
			bad_sd++;
	}
	printf("  winding: worst difference from the exact sum %.2e (accuracy 2), %.2e (accuracy 6)\n", worst_wn, worst_wn6);
	TEST_ASSERT(worst_wn < 5e-3 && worst_wn6 < 1e-4, "sphere: winding numbers within 5e-3 of the exact sum by default, 1e-4 at accuracy 6");
	TEST_ASSERT(bad_wn == 0, "sphere: winding numbers 1 inside, 0 outside");
	TEST_ASSERT(bad_sd == 0, "sphere: signed distances negative inside, their size the closest distance");
	/* the executor: chunks backwards, the same results */
	TEST_ASSERT(qaws_mesh_bvh_closest(h2, b, pts, n, 0, &ex) == QAWS_STATUS_OK, "sphere: closest through an executor");
	for (i = 0; i < n; i++)
		if (h2[i].simplex != h[i].simplex || h2[i].distance != h[i].distance)
			bad_ex++;
	TEST_ASSERT(bad_ex == 0, "sphere: the executor's results equal the calling thread's");
	/* refit after scaling by 2 against a new build */
	for (i = 0; i < 3 * nv; i++)
		x2[i] = 2 * x[i];
	TEST_ASSERT(qaws_mesh_bvh_refit(b, x2) == QAWS_STATUS_OK, "refit: ran");
	d.positions = x2;
	TEST_ASSERT(qaws_mesh_bvh_create(&b2, &d) == QAWS_STATUS_OK, "refit: a new build");
	qaws_mesh_bvh_closest(h, b, pts, n, 0, NULL);
	qaws_mesh_bvh_closest(h2, b2, pts, n, 0, NULL);
	qaws_mesh_bvh_winding_number(wn, b, pts, n, NULL);
	qaws_mesh_bvh_winding_number(sd, b2, pts, n, NULL);
	for (i = 0; i < n; i++)
		if (fabs(h[i].distance - h2[i].distance) > bt_tol() || fabs(wn[i] - sd[i]) > 2e-3)
			bad_rf++;
	TEST_ASSERT(bad_rf == 0, "refit: distances and winding numbers as a new build's");
	qaws_mesh_bvh_destroy(b2);
	b2 = NULL;
	/* a cap removed (the triangles round the north pole): open, the centre still about 1 */
	d.positions = x;
	d.simplices = f + 3 * 48;
	d.simplex_count = nt - 48;
	TEST_ASSERT(qaws_mesh_bvh_create(&b2, &d) == QAWS_STATUS_OK, "open sphere: built");
	{
		qaws_scalar c[6] = { 0, 0, 0, 0, 0, (qaws_scalar)0.995 }, w[2];
		double cap = 0;
		qaws_mesh_bvh_winding_number(w, b2, c, 2, NULL);
		for (t = 0; t < 48; t++)
		{
			double a[3], bb[3], cc[3], o[3] = { 0, 0, 0 };
			bt_get(a, x, f[3 * t]);
			bt_get(bb, x, f[3 * t + 1]);
			bt_get(cc, x, f[3 * t + 2]);
			cap += bt_solid(o, a, bb, cc);
		}
		printf("  open sphere: centre %.4f (1 less the cap's %.4f), under the hole %.3f\n", (double)w[0], cap, (double)w[1]);
		TEST_ASSERT(fabs(w[0] - (1 - cap)) < 1e-3 && w[0] > 0.9, "open sphere: the centre's winding number 1 less the cap's solid angle");
		TEST_ASSERT(w[1] > 0.05 && w[1] < 0.95, "open sphere: a point under the hole between in and out");
	}
	qaws_mesh_bvh_destroy(b2);
	qaws_mesh_bvh_destroy(b);
	free(x);
	free(x2);
	free(f);
	free(pts);
	free(dir);
	free(wn);
	free(wn6);
	free(sd);
	free(h);
	free(h2);
}

static void test_watertight(void)
{
	/* a 16 x 16 grid in z = 0, rays straight down through every vertex and edge midpoint */
	unsigned int n = 16, nv = (n + 1) * (n + 1), i, j, t = 0, r = 0, miss = 0, nr = (2 * n + 1) * (2 * n + 1);
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	unsigned int* f = (unsigned int*)malloc(sizeof(unsigned int) * 6 * n * n);
	qaws_scalar *o = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nr), *dv = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nr);
	qaws_mesh_hit* h = (qaws_mesh_hit*)malloc(sizeof(qaws_mesh_hit) * nr);
	qaws_mesh_bvh_desc d;
	qaws_mesh_bvh* b = NULL;
	for (j = 0; j <= n; j++)
		for (i = 0; i <= n; i++)
		{
			/* a slight irregular warp so the grid's coordinates are not round */
			x[3 * (j * (n + 1) + i)] = (qaws_scalar)(i * 0.1 + 0.013 * sin(1.7 * j));
			x[3 * (j * (n + 1) + i) + 1] = (qaws_scalar)(j * 0.1 + 0.011 * cos(2.3 * i));
			x[3 * (j * (n + 1) + i) + 2] = 0;
		}
	for (j = 0; j < n; j++)
		for (i = 0; i < n; i++)
		{
			unsigned int a = j * (n + 1) + i, b2 = a + 1, c = a + n + 2, e = a + n + 1;
			f[t++] = a; f[t++] = b2; f[t++] = c;
			f[t++] = a; f[t++] = c; f[t++] = e;
		}
	for (j = 0; j <= 2 * n; j++)
		for (i = 0; i <= 2 * n; i++, r++)
		{
			/* vertices, horizontal and vertical edge midpoints and the diagonals' (all on edges) */
			unsigned int i0 = i / 2, j0 = j / 2, i1 = (i + 1) / 2, j1 = (j + 1) / 2, k;
			for (k = 0; k < 2; k++)
				o[3 * r + k] = (qaws_scalar)(0.5 * (x[3 * (j0 * (n + 1) + i0) + k] + x[3 * (j1 * (n + 1) + i1) + k]));
			o[3 * r + 2] = 1;
			dv[3 * r] = 0;
			dv[3 * r + 1] = 0;
			dv[3 * r + 2] = -1;
		}
	memset(&d, 0, sizeof d);
	d.positions = x;
	d.vertex_count = nv;
	d.dimension = 3;
	d.simplices = f;
	d.simplex_size = 3;
	d.simplex_count = 2 * n * n;
	TEST_ASSERT(qaws_mesh_bvh_create(&b, &d) == QAWS_STATUS_OK, "watertight: built");
	qaws_mesh_bvh_ray(h, b, o, dv, nr, 0, 0, NULL);
	/* (the outer border's points may round off the mesh: not counted) */
	for (r = 0; r < nr; r++)
		if (r % (2 * n + 1) && r % (2 * n + 1) != 2 * n && r > 2 * n && r < nr - 2 * n - 1
			&& (h[r].element == QAWS_MESH_NONE || fabs(h[r].distance - 1) > 1e-6))
		{
			if (miss < 3)
			miss++;
		}
	printf("  watertight: %u rays through vertices and edges, %u missed\n", nr, miss);
	TEST_ASSERT(miss == 0, "watertight: no ray through a vertex or an edge slips between the triangles");
	qaws_mesh_bvh_destroy(b);
	free(x);
	free(f);
	free(o);
	free(dv);
	free(h);
}

static void test_2d(void)
{
	/* a 10 x 10 grid of the unit square in triangles */
	unsigned int n = 10, nv = (n + 1) * (n + 1), i, j, t = 0, bad = 0, bad_b = 0, np = 200;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * nv);
	unsigned int* f = (unsigned int*)malloc(sizeof(unsigned int) * 6 * n * n);
	qaws_scalar pts[400], wn[200], sd[200], dir[400];
	qaws_mesh_hit h[200];
	qaws_mesh_2d_desc md;
	qaws_mesh_2d* m = NULL;
	qaws_mesh_bvh *b = NULL, *bb = NULL;
	for (j = 0; j <= n; j++)
		for (i = 0; i <= n; i++)
		{
			x[2 * (j * (n + 1) + i)] = (qaws_scalar)((double)i / n);
			x[2 * (j * (n + 1) + i) + 1] = (qaws_scalar)((double)j / n);
		}
	for (j = 0; j < n; j++)
		for (i = 0; i < n; i++)
		{
			unsigned int a = j * (n + 1) + i, b2 = a + 1, c = a + n + 2, e = a + n + 1;
			f[t++] = a; f[t++] = b2; f[t++] = c;
			f[t++] = a; f[t++] = c; f[t++] = e;
		}
	memset(&md, 0, sizeof md);
	md.positions = x;
	md.vertex_count = nv;
	md.face_vertices = f;
	md.face_count = 2 * n * n;
	TEST_ASSERT(qaws_mesh_2d_create(&m, &md) == QAWS_STATUS_OK, "2D: mesh made");
	TEST_ASSERT(qaws_mesh_bvh_create_mesh_2d(&b, m, 0, NULL) == QAWS_STATUS_OK, "2D: faces' BVH built");
	TEST_ASSERT(qaws_mesh_bvh_create_mesh_2d(&bb, m, 1, NULL) == QAWS_STATUS_OK, "2D: border's BVH built");
	for (i = 0; i < np; i++)
	{
		pts[2 * i] = (qaws_scalar)(1.4 * bt_rand() - 0.2);
		pts[2 * i + 1] = (qaws_scalar)(1.4 * bt_rand() - 0.2);
	}
	qaws_mesh_bvh_closest(h, b, pts, np, 0, NULL);
	for (i = 0; i < np; i++)
	{
		double px = pts[2 * i], py = pts[2 * i + 1], dx = px < 0 ? -px : px > 1 ? px - 1 : 0, dy = py < 0 ? -py : py > 1 ? py - 1 : 0;
		double q[2] = { 0, 0 };
		unsigned int k;
		for (k = 0; k < 3; k++)
		{
			q[0] += h[i].weights[k] * x[2 * h[i].vertices[k]];
			q[1] += h[i].weights[k] * x[2 * h[i].vertices[k] + 1];
		}
		if (fabs(h[i].distance - sqrt(dx * dx + dy * dy)) > 1e-5 || fabs(q[0] - h[i].point[0]) > 1e-5 || fabs(q[1] - h[i].point[1]) > 1e-5)
		{
			bad++;
		}
		if (!dx && !dy)
		{
			/* located: the element is the triangle holding the point */
			unsigned int ci = (unsigned int)(px * n), cj = (unsigned int)(py * n);
			ci = ci >= n ? n - 1 : ci;
			cj = cj >= n ? n - 1 : cj;
			if (h[i].element / 2 != cj * n + ci && h[i].distance == 0)
			{
				/* on a cell's border it may be the neighbour's: accept when the point is on its edge */
				double fx = px * n - ci, fy = py * n - cj;
				if (fx > 1e-6 && fx < 1 - 1e-6 && fy > 1e-6 && fy < 1 - 1e-6)
					bad++;
			}
		}
	}
	TEST_ASSERT(bad == 0, "2D: points located in their triangle, outside at their distance from the square");
	qaws_mesh_bvh_winding_number(wn, bb, pts, np, NULL);
	qaws_mesh_bvh_signed_distance(sd, NULL, bb, pts, np, NULL);
	for (i = 0; i < np; i++)
	{
		double px = pts[2 * i], py = pts[2 * i + 1];
		int in = px > 0 && px < 1 && py > 0 && py < 1;
		double dist = in ? fmin(fmin(px, 1 - px), fmin(py, 1 - py)) : h[i].distance;
		if ((in && fabs(wn[i] - 1) > 1e-2) || (!in && fabs(wn[i]) > 1e-2) || fabs(fabs(sd[i]) - dist) > 1e-5 || (in != (sd[i] < 0)))
		{
			bad_b++;
		}
	}
	TEST_ASSERT(bad_b == 0, "2D border: winding 1 inside and 0 outside within 1e-2, signed distances");
	/* rays from the centre in every direction hit the border at the square's distance */
	for (i = 0; i < 100; i++)
	{
		double a = 2 * BT_PI * (i + 0.5) / 100;
		pts[2 * i] = (qaws_scalar)0.5;
		pts[2 * i + 1] = (qaws_scalar)0.5;
		dir[2 * i] = (qaws_scalar)cos(a);
		dir[2 * i + 1] = (qaws_scalar)sin(a);
	}
	qaws_mesh_bvh_ray(h, bb, pts, dir, 100, 0, 0, NULL);
	for (i = 0, bad = 0; i < 100; i++)
	{
		double a = 2 * BT_PI * (i + 0.5) / 100, tt = fmin(0.5 / fabs(cos(a)), 0.5 / fabs(sin(a)));
		if (h[i].element == QAWS_MESH_NONE || fabs(h[i].distance - tt) > 1e-5)
			bad++;
	}
	TEST_ASSERT(bad == 0, "2D border: rays from the centre hit it at the square's distance");
	qaws_mesh_bvh_destroy(b);
	qaws_mesh_bvh_destroy(bb);
	qaws_mesh_2d_destroy(m);
	free(x);
	free(f);
}

static void test_volume(void)
{
	/* a 3 x 3 x 3 block of unit hexahedra */
	unsigned int n = 3, nv = 64, i, j, k, c = 0, bad = 0;
	qaws_scalar x[192], pts[300];
	unsigned int cv[27 * 8];
	qaws_cell_type ty[27];
	qaws_mesh_hit h[100];
	qaws_volume_mesh_desc vd;
	qaws_volume_mesh* v = NULL;
	qaws_mesh_bvh* b = NULL;
	for (k = 0; k <= n; k++)
		for (j = 0; j <= n; j++)
			for (i = 0; i <= n; i++)
			{
				unsigned int id = (k * 4 + j) * 4 + i;
				x[3 * id] = (qaws_scalar)i;
				x[3 * id + 1] = (qaws_scalar)j;
				x[3 * id + 2] = (qaws_scalar)k;
			}
	for (k = 0; k < n; k++)
		for (j = 0; j < n; j++)
			for (i = 0; i < n; i++, c++)
			{
				unsigned int a = (k * 4 + j) * 4 + i;
				cv[8 * c] = a;
				cv[8 * c + 1] = a + 1;
				cv[8 * c + 2] = a + 5;
				cv[8 * c + 3] = a + 4;
				cv[8 * c + 4] = a + 16;
				cv[8 * c + 5] = a + 17;
				cv[8 * c + 6] = a + 21;
				cv[8 * c + 7] = a + 20;
				ty[c] = QAWS_CELL_HEXAHEDRON;
			}
	memset(&vd, 0, sizeof vd);
	vd.positions = x;
	vd.vertex_count = nv;
	vd.cell_vertices = cv;
	vd.cell_types = ty;
	vd.cell_count = 27;
	TEST_ASSERT(qaws_volume_mesh_create(&v, &vd) == QAWS_STATUS_OK, "volume: block made");
	TEST_ASSERT(qaws_mesh_bvh_create_volume(&b, v, NULL) == QAWS_STATUS_OK, "volume: BVH built");
	for (i = 0; i < 99; i++)
		pts[i] = (qaws_scalar)(0.02 + 2.96 * bt_rand());
	pts[99] = 0;
	qaws_mesh_bvh_closest(h, b, pts, 33, 0, NULL);
	for (i = 0; i < 33; i++)
	{
		unsigned int ci = (unsigned int)pts[3 * i], cj = (unsigned int)pts[3 * i + 1], ck = (unsigned int)pts[3 * i + 2], q, m;
		double p[3] = { 0, 0, 0 };
		for (q = 0; q < 4; q++)
			for (m = 0; m < 3; m++)
				p[m] += h[i].weights[q] * x[3 * h[i].vertices[q] + m];
		if (h[i].distance != 0 || h[i].element != (ck * 3 + cj) * 3 + ci || fabs(p[0] - pts[3 * i]) > 1e-5 || fabs(p[1] - pts[3 * i + 1]) > 1e-5
			|| fabs(p[2] - pts[3 * i + 2]) > 1e-5)
			bad++;
	}
	TEST_ASSERT(bad == 0, "volume: points located in their cell, the weights giving them");
	{
		qaws_scalar o[3] = { 4, (qaws_scalar)1.5, (qaws_scalar)1.5 };
		qaws_mesh_hit ho;
		qaws_mesh_bvh_closest(&ho, b, o, 1, 0, NULL);
		TEST_ASSERT(fabs(ho.distance - 1) < 1e-6 && ho.point[0] == 3, "volume: a point outside at its distance from the block");
		TEST_ASSERT(qaws_mesh_bvh_closest(&ho, b, o, 1, (qaws_scalar)0.5, NULL) == QAWS_STATUS_OK && ho.element == QAWS_MESH_NONE,
			"volume: nothing within a smaller max_distance");
		TEST_ASSERT(qaws_mesh_bvh_ray(&ho, b, o, o, 1, 0, 0, NULL) == QAWS_STATUS_UNSUPPORTED_OPERATION, "refusals: rays on tetrahedra");
		TEST_ASSERT(qaws_mesh_bvh_winding_number(o, b, o, 1, NULL) == QAWS_STATUS_UNSUPPORTED_OPERATION, "refusals: winding on tetrahedra");
	}
	qaws_mesh_bvh_destroy(b);
	qaws_volume_mesh_destroy(v);
}

static void test_refusals(void)
{
	qaws_scalar x[9] = { 0, 0, 0, 1, 0, 0, 0, 1, 0 };
	unsigned int t[4] = { 0, 1, 2, 0 };
	qaws_mesh_bvh_desc d;
	qaws_mesh_bvh* b = NULL;
	memset(&d, 0, sizeof d);
	d.positions = x;
	d.vertex_count = 3;
	d.dimension = 4;
	d.simplices = t;
	d.simplex_size = 3;
	d.simplex_count = 1;
	TEST_ASSERT(qaws_mesh_bvh_create(&b, &d) == QAWS_STATUS_INVALID_ARGUMENT && !b, "refusals: dimension 4");
	d.dimension = 2;
	d.simplex_size = 4;
	TEST_ASSERT(qaws_mesh_bvh_create(&b, &d) == QAWS_STATUS_INVALID_ARGUMENT && !b, "refusals: tetrahedra in 2D");
	d.dimension = 3;
	d.simplex_size = 3;
	t[2] = 7;
	TEST_ASSERT(qaws_mesh_bvh_create(&b, &d) == QAWS_STATUS_INVALID_ARGUMENT && !b, "refusals: an index out of range");
	t[2] = 2;
	x[4] = (qaws_scalar)NAN;
	TEST_ASSERT(qaws_mesh_bvh_create(&b, &d) == QAWS_STATUS_INVALID_ARGUMENT && !b, "refusals: a position not finite");
	x[4] = 0;
	d.dimension = 2;
	d.vertex_count = 4;
	TEST_ASSERT(qaws_mesh_bvh_create(&b, &d) == QAWS_STATUS_OK, "refusals: 2D triangles taken");
	{
		qaws_scalar w, p[2] = { 0, 0 };
		TEST_ASSERT(qaws_mesh_bvh_winding_number(&w, b, p, 1, NULL) == QAWS_STATUS_UNSUPPORTED_OPERATION, "refusals: winding on 2D triangles");
	}
	qaws_mesh_bvh_destroy(b);
}

static void test_scale(void)
{
	qaws_scalar *x, *p;
	unsigned int nv, nt, *f, i, n = 100000;
	qaws_mesh_bvh_desc d;
	qaws_mesh_bvh* b = NULL;
	qaws_mesh_hit* h;
	clock_t c0, c1, c2;
	bt_sphere(400, 250, &x, &nv, &f, &nt);
	p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n);
	h = (qaws_mesh_hit*)malloc(sizeof(qaws_mesh_hit) * n);
	/* points near the surface (within 0.1 of it), as remeshing and repair ask */
	for (i = 0; i < n; i++)
	{
		double u = 2 * bt_rand() - 1, a = 2 * BT_PI * bt_rand(), s = sqrt(1 - u * u), r = 0.9 + 0.2 * bt_rand();
		p[3 * i] = (qaws_scalar)(r * s * cos(a));
		p[3 * i + 1] = (qaws_scalar)(r * s * sin(a));
		p[3 * i + 2] = (qaws_scalar)(r * u);
	}
	memset(&d, 0, sizeof d);
	d.positions = x;
	d.vertex_count = nv;
	d.dimension = 3;
	d.simplices = f;
	d.simplex_size = 3;
	d.simplex_count = nt;
	c0 = clock();
	TEST_ASSERT(qaws_mesh_bvh_create(&b, &d) == QAWS_STATUS_OK, "scale: built");
	c1 = clock();
	qaws_mesh_bvh_closest(h, b, p, n, 0, NULL);
	c2 = clock();
	printf("  scale: %u triangles built in %.3f s, %u closest points near the surface in %.3f s\n", nt, (double)(c1 - c0) / CLOCKS_PER_SEC, n,
		(double)(c2 - c1) / CLOCKS_PER_SEC);
	qaws_mesh_bvh_destroy(b);
	free(x);
	free(f);
	free(p);
	free(h);
}

int test_159_mesh_bvh_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 159: spatial queries on meshes\n");
	test_sphere();
	test_watertight();
	test_2d();
	test_volume();
	test_refusals();
	test_scale();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
