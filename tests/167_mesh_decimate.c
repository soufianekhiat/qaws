/*
 * Test 167: decimation by quadric error (qaws_mesh_3d_decimate), small meshes
 *
 *   - a unit icosphere (1280 triangles) down to 320: the count reached, still
 *     closed, its certified deviation from the input within 2% of the radius;
 *     an error bound stops it with every collapse under the bound;
 *   - a flat grid: down to a handful of triangles, its area and border
 *     unchanged (the plane costs nothing); with the border locked, its
 *     border vertices stay, unmoved;
 *   - a cube in 6 x 6 grids per side with its creases as features: down to
 *     a few triangles per side, its volume exact;
 *   - silhouettes: a sphere rippled round z decimated seen along z or from an
 *     eye on z keeps that outline far better than without; 4 x 4 octahedral
 *     views from above put the faces on the upper half (the lower one
 *     unseen);
 *   - creases found by angle on a cube without features; importance 0 keeps
 *     vertices; a tight bound on the faces' turn stops earlier;
 *   - refusals: no mesh, a negative error bound.
 */

#include "test_common.h"
#include "qaws_mesh_decimate.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

static qaws_mesh_3d* dm_icosphere(unsigned int levels)
{
	static int const ic[60] = { 0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11, 1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8, 3, 9, 4, 3, 4, 2,
		3, 2, 6, 3, 6, 8, 3, 8, 9, 4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7, 9, 8, 1 };
	double p = (1 + sqrt(5.0)) / 2;
	double const v0[36] = { -1, p, 0, 1, p, 0, -1, -p, 0, 1, -p, 0, 0, -1, p, 0, 1, p, 0, -1, -p, 0, 1, -p, p, 0, -1, p, 0, 1, -p, 0, -1, -p, 0, 1 };
	unsigned int nt = 20, nv = 12, i, j, lv, nm, cap = 20, (*mid)[3];
	double* x;
	unsigned int *t, *t2;
	qaws_scalar* xs;
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	for (lv = 0; lv < levels; lv++)
		cap *= 4;
	x = (double*)malloc(sizeof(double) * 3 * (cap / 2 + 2));
	t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * cap);
	t2 = (unsigned int*)malloc(sizeof(unsigned int) * 3 * cap);
	mid = (unsigned int(*)[3])malloc(sizeof(unsigned int) * 3 * 2 * cap);
	for (i = 0; i < 12; i++)
	{
		double l = sqrt(v0[3 * i] * v0[3 * i] + v0[3 * i + 1] * v0[3 * i + 1] + v0[3 * i + 2] * v0[3 * i + 2]);
		for (j = 0; j < 3; j++)
			x[3 * i + j] = v0[3 * i + j] / l;
	}
	for (i = 0; i < 60; i++)
		t[i] = (unsigned int)ic[i];
	for (lv = 0; lv < levels; lv++)
	{
		nm = 0;
		for (i = 0; i < nt; i++)
		{
			unsigned int c[3], e;
			for (e = 0; e < 3; e++)
			{
				unsigned int a = t[3 * i + e], b = t[3 * i + (e + 1) % 3], lo = a < b ? a : b, hi = a < b ? b : a, s;
				double l;
				for (s = 0; s < nm && !(mid[s][0] == lo && mid[s][1] == hi); s++)
					;
				if (s == nm)
				{
					mid[nm][0] = lo;
					mid[nm][1] = hi;
					mid[nm++][2] = nv;
					for (j = 0; j < 3; j++)
						x[3 * nv + j] = 0.5 * (x[3 * a + j] + x[3 * b + j]);
					l = sqrt(x[3 * nv] * x[3 * nv] + x[3 * nv + 1] * x[3 * nv + 1] + x[3 * nv + 2] * x[3 * nv + 2]);
					for (j = 0; j < 3; j++)
						x[3 * nv + j] /= l;
					nv++;
				}
				c[e] = mid[s][2];
			}
			{
				unsigned int const q[12] = { t[3 * i], c[0], c[2], c[0], t[3 * i + 1], c[1], c[2], c[1], t[3 * i + 2], c[0], c[1], c[2] };
				memcpy(t2 + 12 * i, q, sizeof q);
			}
		}
		nt *= 4;
		memcpy(t, t2, sizeof(unsigned int) * 3 * nt);
	}
	xs = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	for (i = 0; i < 3 * nv; i++)
		xs[i] = (qaws_scalar)x[i];
	memset(&d, 0, sizeof d);
	d.positions = xs;
	d.vertex_count = nv;
	d.face_vertices = t;
	d.face_count = nt;
	qaws_mesh_3d_create(&m, &d);
	free(x); free(t); free(t2); free(mid); free(xs);
	return m;
}

static qaws_mesh_3d* dm_grid(unsigned int n)
{
	unsigned int nv = (n + 1) * (n + 1), i, j, k = 0;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 6 * n * n);
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	for (j = 0; j <= n; j++)
		for (i = 0; i <= n; i++)
		{
			x[3 * (j * (n + 1) + i)] = (qaws_scalar)((double)i / n);
			x[3 * (j * (n + 1) + i) + 1] = (qaws_scalar)((double)j / n);
			x[3 * (j * (n + 1) + i) + 2] = 0;
		}
	for (j = 0; j < n; j++)
		for (i = 0; i < n; i++)
		{
			unsigned int v = j * (n + 1) + i;
			t[k++] = v; t[k++] = v + 1; t[k++] = v + n + 2;
			t[k++] = v; t[k++] = v + n + 2; t[k++] = v + n + 1;
		}
	memset(&d, 0, sizeof d);
	d.positions = x;
	d.vertex_count = nv;
	d.face_vertices = t;
	d.face_count = 2 * n * n;
	qaws_mesh_3d_create(&m, &d);
	free(x);
	free(t);
	return m;
}

/* the unit cube's surface in n x n grids, vertices shared along its edges (by their integer coordinates) */
static qaws_mesh_3d* dm_cube(unsigned int n)
{
	unsigned int s = n + 1, nv = 0, nt = 0, f, i, j, k;
	unsigned int* id = (unsigned int*)malloc(sizeof(unsigned int) * s * s * s);
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * s * s * s);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 36 * n * n);
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	for (i = 0; i < s * s * s; i++)
		id[i] = 0xFFFFFFFFu;
	for (f = 0; f < 6; f++)
	{
		/* side f: axis f / 2 fixed at 0 or n, the other two running so the normal points out */
		unsigned int ax = f / 2, u = (ax + 1) % 3, w = (ax + 2) % 3, hi = f % 2;
		for (j = 0; j < n; j++)
			for (i = 0; i < n; i++)
			{
				unsigned int q[4], c;
				for (c = 0; c < 4; c++)
				{
					unsigned int g[3], key;
					g[ax] = hi ? n : 0;
					g[u] = i + (c == 1 || c == 2);
					g[w] = j + (c >= 2);
					key = (g[0] * s + g[1]) * s + g[2];
					if (id[key] == 0xFFFFFFFFu)
					{
						for (k = 0; k < 3; k++)
							x[3 * nv + k] = (qaws_scalar)((double)g[k] / n);
						id[key] = nv++;
					}
					q[c] = id[key];
				}
				if (hi)
				{
					t[3 * nt] = q[0]; t[3 * nt + 1] = q[1]; t[3 * nt + 2] = q[2]; nt++;
					t[3 * nt] = q[0]; t[3 * nt + 1] = q[2]; t[3 * nt + 2] = q[3]; nt++;
				}
				else
				{
					t[3 * nt] = q[0]; t[3 * nt + 1] = q[2]; t[3 * nt + 2] = q[1]; nt++;
					t[3 * nt] = q[0]; t[3 * nt + 1] = q[3]; t[3 * nt + 2] = q[2]; nt++;
				}
			}
	}
	memset(&d, 0, sizeof d);
	d.positions = x;
	d.vertex_count = nv;
	d.face_vertices = t;
	d.face_count = nt;
	qaws_mesh_3d_create(&m, &d);
	free(id); free(x); free(t);
	return m;
}

static double dm_volume(qaws_mesh_3d const* m, double* area)
{
	unsigned int nt = qaws_mesh_3d_get_triangle_count(m), i;
	double v = 0;
	*area = 0;
	for (i = 0; i < nt; i++)
	{
		unsigned int c[3];
		qaws_vec3 a, b, e;
		double u[3], w[3], n[3];
		qaws_mesh_3d_get_triangle(c, m, i);
		a = qaws_mesh_3d_get_vertex(m, c[0]);
		b = qaws_mesh_3d_get_vertex(m, c[1]);
		e = qaws_mesh_3d_get_vertex(m, c[2]);
		v += ((double)a.x * (b.y * e.z - b.z * e.y) + (double)a.y * (b.z * e.x - b.x * e.z) + (double)a.z * (b.x * e.y - b.y * e.x)) / 6;
		u[0] = b.x - a.x; u[1] = b.y - a.y; u[2] = b.z - a.z;
		w[0] = e.x - a.x; w[1] = e.y - a.y; w[2] = e.z - a.z;
		n[0] = u[1] * w[2] - u[2] * w[1];
		n[1] = u[2] * w[0] - u[0] * w[2];
		n[2] = u[0] * w[1] - u[1] * w[0];
		*area += 0.5 * sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
	}
	return v;
}

static void test_sphere(void)
{
	qaws_mesh_3d* m = dm_icosphere(3);
	qaws_decimate_desc d;
	qaws_decimate_report r;
	qaws_mesh_topology top;
	qaws_status st;
	char msg[256];
	memset(&d, 0, sizeof d);
	d.target_faces = 320;
	d.flags = QAWS_DECIMATE_MEASURE;
	st = qaws_mesh_3d_decimate(&r, m, &d);
	qaws_mesh_3d_get_topology(&top, m);
	sprintf(msg, "icosphere %u -> %u faces (%u collapses): closed %d, genus %u, largest error %.4f, deviation from the input in [%.4f, %.4f]",
		r.faces_before, r.faces_after, r.collapses, top.closed, top.genus, (double)r.max_error, (double)r.deviation.hausdorff[0],
		(double)r.deviation.hausdorff[1]);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && r.faces_after <= 320 && r.faces_after >= 300 && top.closed && top.genus == 0
		&& r.deviation.hausdorff[1] < 0.02, msg);
	qaws_mesh_3d_destroy(m);
	m = dm_icosphere(3);
	memset(&d, 0, sizeof d);
	d.max_error = (qaws_scalar)0.004;
	st = qaws_mesh_3d_decimate(&r, m, &d);
	sprintf(msg, "icosphere with an error bound 0.004: %u -> %u faces, the largest error %.4f", r.faces_before, r.faces_after, (double)r.max_error);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && r.max_error <= (qaws_scalar)0.004 && r.faces_after < r.faces_before && r.faces_after > 100, msg);
	qaws_mesh_3d_destroy(m);
}

static void test_plane(void)
{
	unsigned int n = 16, k;
	char msg[256];
	for (k = 0; k < 2; k++)
	{
		qaws_mesh_3d* m = dm_grid(n);
		qaws_decimate_desc d;
		qaws_decimate_report r;
		double area, lo[2] = { 1e300, 1e300 }, hi[2] = { -1e300, -1e300 };
		unsigned int nv, i, nb = 0, ok = 1;
		qaws_status st;
		memset(&d, 0, sizeof d);
		d.target_faces = 2;
		d.max_error = (qaws_scalar)1e-9;
		d.flags = k ? QAWS_DECIMATE_LOCK_BORDER : 0;
		st = qaws_mesh_3d_decimate(&r, m, &d);
		dm_volume(m, &area);
		nv = qaws_mesh_3d_get_vertex_count(m);
		for (i = 0; i < nv; i++)
		{
			qaws_vec3 p = qaws_mesh_3d_get_vertex(m, i);
			lo[0] = fmin(lo[0], p.x);
			lo[1] = fmin(lo[1], p.y);
			hi[0] = fmax(hi[0], p.x);
			hi[1] = fmax(hi[1], p.y);
			ok &= fabs(p.z) < 1e-12;
			if (qaws_mesh_3d_is_boundary_vertex(m, i))
			{
				double gx = p.x * n, gy = p.y * n;
				nb++;
				ok &= fabs(gx - floor(gx + 0.5)) < 1e-12 && fabs(gy - floor(gy + 0.5)) < 1e-12;   /* (on the grid: unmoved) */
			}
		}
		sprintf(msg, "flat grid %s: %u -> %u faces, area %.12f, box [%g, %g] x [%g, %g], %u border vertices", k ? "with its border locked" : "",
			r.faces_before, r.faces_after, area, lo[0], hi[0], lo[1], hi[1], nb);
		printf("  %s\n", msg);
		TEST_ASSERT(st == QAWS_STATUS_OK && ok && fabs(area - 1) < 1e-12 && lo[0] == 0 && lo[1] == 0 && hi[0] == 1 && hi[1] == 1
			&& (k ? nb == 4 * n && r.faces_after < r.faces_before : r.faces_after <= 8), msg);
		qaws_mesh_3d_destroy(m);
	}
}

static void test_cube(void)
{
	qaws_mesh_3d* m = dm_cube(6);
	qaws_decimate_desc d;
	qaws_decimate_report r;
	double area, vol;
	qaws_status st;
	char msg[200];
	qaws_mesh_3d_detect_features(m, (qaws_scalar)0.5, (qaws_scalar)0.5, 0);
	memset(&d, 0, sizeof d);
	d.target_faces = 12;
	st = qaws_mesh_3d_decimate(&r, m, &d);
	vol = dm_volume(m, &area);
	sprintf(msg, "a cube in 6 x 6 grids, its creases features: %u -> %u faces, volume %.12f, area %.12f", r.faces_before, r.faces_after, vol, area);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && r.faces_after <= 48 && fabs(vol - 1) < 1e-12 && fabs(area - 6) < 1e-12, msg);
	qaws_mesh_3d_destroy(m);
}

/* the icosphere rippled round z (its outline along z wavy) */
static qaws_mesh_3d* dm_rippled(void)
{
	qaws_mesh_3d* m = dm_icosphere(3);
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), i;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	qaws_mesh_3d_get_positions(x, m);
	for (i = 0; i < nv; i++)
	{
		double r = 1 + 0.08 * sin(6 * atan2((double)x[3 * i + 1], (double)x[3 * i]));
		x[3 * i] = (qaws_scalar)(r * x[3 * i]);
		x[3 * i + 1] = (qaws_scalar)(r * x[3 * i + 1]);
		x[3 * i + 2] = (qaws_scalar)(r * x[3 * i + 2]);
	}
	qaws_mesh_3d_set_positions(m, x);
	free(x);
	return m;
}

/* the outline seen along axis ax: the largest gap between the two meshes' support functions in the other two axes' plane */
static double dm_outline_gap(qaws_mesh_3d const* a, qaws_mesh_3d const* b, unsigned int ax)
{
	unsigned int s, i, u = (ax + 1) % 3, w = (ax + 2) % 3;
	double gap = 0;
	for (s = 0; s < 180; s++)
	{
		double th = 2 * 3.14159265358979323846 * s / 180, h[2] = { -1e300, -1e300 };
		unsigned int k;
		for (k = 0; k < 2; k++)
		{
			qaws_mesh_3d const* m = k ? b : a;
			for (i = 0; i < qaws_mesh_3d_get_vertex_count(m); i++)
			{
				qaws_vec3 p = qaws_mesh_3d_get_vertex(m, i);
				double c[3] = { p.x, p.y, p.z }, d = c[u] * cos(th) + c[w] * sin(th);
				h[k] = d > h[k] ? d : h[k];
			}
		}
		gap = fmax(gap, fabs(h[0] - h[1]));
	}
	return gap;
}

/* the outline gap seen along direction d (support functions over 90 angles in the plane across d) */
static double dm_outline_gap_dir(qaws_mesh_3d const* a, qaws_mesh_3d const* b, double const* d)
{
	double u[3], w[3], l, gap = 0;
	unsigned int s, i, k;
	u[0] = fabs(d[0]) < 0.9 ? 0 : -d[2];
	u[1] = fabs(d[0]) < 0.9 ? -d[2] : 0;
	u[2] = fabs(d[0]) < 0.9 ? d[1] : d[0];
	l = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
	for (k = 0; k < 3; k++)
		u[k] /= l;
	w[0] = d[1] * u[2] - d[2] * u[1];
	w[1] = d[2] * u[0] - d[0] * u[2];
	w[2] = d[0] * u[1] - d[1] * u[0];
	for (s = 0; s < 90; s++)
	{
		double th = 2 * 3.14159265358979323846 * s / 90, h2[2] = { -1e300, -1e300 }, e[3];
		for (k = 0; k < 3; k++)
			e[k] = u[k] * cos(th) + w[k] * sin(th);
		for (k = 0; k < 2; k++)
		{
			qaws_mesh_3d const* m = k ? b : a;
			for (i = 0; i < qaws_mesh_3d_get_vertex_count(m); i++)
			{
				qaws_vec3 p = qaws_mesh_3d_get_vertex(m, i);
				double v = p.x * e[0] + p.y * e[1] + p.z * e[2];
				h2[k] = v > h2[k] ? v : h2[k];
			}
		}
		gap = fmax(gap, fabs(h2[0] - h2[1]));
	}
	return gap;
}

/* the mean outline gap over the 26 directions of a cube's faces, edges and corners: those with z > 0 (side 1), z < 0 (-1), or all (0) */
static double dm_outline_mean(qaws_mesh_3d const* a, qaws_mesh_3d const* b, int side)
{
	int x, y, z;
	double s = 0;
	unsigned int n = 0;
	for (x = -1; x <= 1; x++)
		for (y = -1; y <= 1; y++)
			for (z = -1; z <= 1; z++)
				if ((x || y || z) && (side == 0 || side * z > 0))
				{
					double d[3] = { x, y, z }, l = sqrt((double)(x * x + y * y + z * z));
					d[0] /= l; d[1] /= l; d[2] /= l;
					s += dm_outline_gap_dir(a, b, d);
					n++;
				}
	return s / n;
}

static void test_silhouette(void)
{
	qaws_mesh_3d* in = dm_rippled();
	unsigned int k;
	double gap[4], up[2] = { 0, 0 };
	unsigned int above[2] = { 0, 0 };
	char msg[600];
	qaws_scalar const zdir[3] = { 0, 0, 1 }, eye[3] = { 0, 0, 6 };
	for (k = 0; k < 4; k++)
	{
		/* 0: plain, 1: along z, 2: from an eye on z, 3: octahedral 4 x 4 */
		qaws_mesh_3d* m = dm_rippled();
		qaws_decimate_desc d;
		qaws_decimate_report r;
		qaws_status st;
		memset(&d, 0, sizeof d);
		d.target_faces = 200;
		if (k == 1 || k == 2)
		{
			d.views = k == 1 ? zdir : eye;
			d.view_count = 1;
			d.flags = k == 2 ? QAWS_DECIMATE_VIEW_POINTS : 0;
		}
		if (k == 3)
		{
			d.octahedral_views = 4;
			d.flags = QAWS_DECIMATE_HEMISPHERE;
		}
		st = qaws_mesh_3d_decimate(&r, m, &d);
		gap[k] = dm_outline_gap(in, m, 2);
		if (k == 0 || k == 3)
		{
			/* the outlines all round, and where the faces went: above (seen) or below */
			unsigned int f, c3[3];
			up[k ? 1 : 0] = dm_outline_mean(in, m, 0);
			for (f = 0; f < qaws_mesh_3d_get_triangle_count(m); f++)
			{
				qaws_mesh_3d_get_triangle(c3, m, f);
				above[k ? 1 : 0] += qaws_mesh_3d_get_vertex(m, c3[0]).z + qaws_mesh_3d_get_vertex(m, c3[1]).z + qaws_mesh_3d_get_vertex(m, c3[2]).z > 0;
			}
		}
		TEST_ASSERT(st == QAWS_STATUS_OK && r.faces_after <= 200 && (k == 0 || r.silhouette_edges > 0), "rippled sphere decimated");
		qaws_mesh_3d_destroy(m);
	}
	sprintf(msg, "rippled sphere to 200 faces, its outline along z off by %.4f (plain), %.4f (seen along z), %.4f (from an eye on z); the three axes' "
		"outlines all round (26 directions) off by %.4f on average plainly, %.4f with 4 x 4 octahedral views from above; faces in the upper half "
		"%u of 200 plainly, %u with those views", gap[0], gap[1], gap[2], up[0], up[1], above[0], above[1]);
	printf("  %s\n", msg);
	TEST_ASSERT(gap[1] < 0.7 * gap[0] && gap[2] < 0.7 * gap[0] && above[1] > above[0] + 20 && up[1] < up[0], msg);
	qaws_mesh_3d_destroy(in);
}

static void test_options(void)
{
	char msg[240];
	{
		/* creases found by angle on a cube without features */
		qaws_mesh_3d* m = dm_cube(6);
		qaws_decimate_desc d;
		qaws_decimate_report r;
		double area, vol;
		qaws_status st;
		memset(&d, 0, sizeof d);
		d.target_faces = 12;
		d.crease_angle = (qaws_scalar)0.5;
		st = qaws_mesh_3d_decimate(&r, m, &d);
		vol = dm_volume(m, &area);
		sprintf(msg, "a cube without features, creases by angle: %u creases found, %u -> %u faces, %u feature edges left, volume %.12f", r.creases,
			r.faces_before, r.faces_after, qaws_mesh_3d_get_feature_edge_count(m), vol);
		printf("  %s\n", msg);
		TEST_ASSERT(st == QAWS_STATUS_OK && r.creases == 72 && r.faces_after <= 16 && qaws_mesh_3d_get_feature_edge_count(m) >= 12
			&& fabs(vol - 1) < 1e-12, msg);
		qaws_mesh_3d_destroy(m);
	}
	{
		/* importance 0 on the cap z > 0.5: its vertices stay where they were */
		qaws_mesh_3d* m = dm_icosphere(3);
		unsigned int nv = qaws_mesh_3d_get_vertex_count(m), i, j, kept = 0, locked = 0, nv2;
		qaws_scalar* imp = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv);
		qaws_scalar* x0 = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
		qaws_decimate_desc d;
		qaws_decimate_report r;
		qaws_mesh_3d_get_positions(x0, m);
		for (i = 0; i < nv; i++)
		{
			imp[i] = (qaws_scalar)(x0[3 * i + 2] > 0.5 ? 0 : 1);
			locked += x0[3 * i + 2] > 0.5;
		}
		memset(&d, 0, sizeof d);
		d.target_faces = 200;
		d.importance = imp;
		qaws_mesh_3d_decimate(&r, m, &d);
		nv2 = qaws_mesh_3d_get_vertex_count(m);
		for (i = 0; i < nv; i++)
			if (x0[3 * i + 2] > 0.5)
				for (j = 0; j < nv2; j++)
				{
					qaws_vec3 p = qaws_mesh_3d_get_vertex(m, j);
					if (p.x == x0[3 * i] && p.y == x0[3 * i + 1] && p.z == x0[3 * i + 2])
					{
						kept++;
						break;
					}
				}
		sprintf(msg, "importance 0 on a cap: %u of its %u vertices still there, unmoved; %u -> %u faces", kept, locked, r.faces_before, r.faces_after);
		printf("  %s\n", msg);
		TEST_ASSERT(kept == locked && r.faces_after < r.faces_before, msg);
		qaws_mesh_3d_destroy(m);
		free(imp);
		free(x0);
	}
	{
		/* a tight bound on the faces' turn stops earlier */
		qaws_mesh_3d* a = dm_rippled();
		qaws_mesh_3d* b = dm_rippled();
		qaws_decimate_desc d;
		qaws_decimate_report ra, rb;
		memset(&d, 0, sizeof d);
		d.target_faces = 60;
		qaws_mesh_3d_decimate(&ra, a, &d);
		d.max_normal_angle = (qaws_scalar)0.15;
		qaws_mesh_3d_decimate(&rb, b, &d);
		sprintf(msg, "faces turning at most 0.15 rad: stops at %u faces (%u with the default bound)", rb.faces_after, ra.faces_after);
		printf("  %s\n", msg);
		TEST_ASSERT(rb.faces_after > ra.faces_after, msg);
		qaws_mesh_3d_destroy(a);
		qaws_mesh_3d_destroy(b);
	}
}

/* a UV seam: the grid's corners' UVs cut apart along x = 1/2 (the right half's moved by 10); decimated, no face crosses it */
static void test_seams(void)
{
	unsigned int n = 16, nf, nc, f, i, cross = 0, ends = 0, id;
	qaws_mesh_3d* m = dm_grid(n);
	unsigned int* fv;
	qaws_scalar* uv;
	qaws_attribute_desc ad;
	qaws_decimate_desc d;
	qaws_decimate_report r;
	qaws_status st;
	char msg[300];
	nf = qaws_mesh_3d_get_face_count(m);
	nc = qaws_mesh_3d_get_corner_count(m);
	fv = (unsigned int*)malloc(sizeof(unsigned int) * nc);
	uv = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * nc);
	qaws_mesh_3d_get_faces(fv, NULL, m);
	for (f = 0; f < nf; f++)
	{
		double cx = 0;
		for (i = 0; i < 3; i++) cx += (double)qaws_mesh_3d_get_vertex(m, fv[3 * f + i]).x / 3;
		for (i = 0; i < 3; i++)
		{
			qaws_vec3 q = qaws_mesh_3d_get_vertex(m, fv[3 * f + i]);
			uv[2 * (3 * f + i)] = (qaws_scalar)((double)q.x + (cx > 0.5 ? 10 : 0));
			uv[2 * (3 * f + i) + 1] = q.y;
		}
	}
	memset(&ad, 0, sizeof ad);
	ad.name = "uv";
	ad.element = QAWS_MESH_CORNER;
	ad.components = 2;
	ad.rule = QAWS_ATTRIBUTE_LINEAR;
	qaws_mesh_3d_add_attribute(&id, m, uv, &ad);
	memset(&d, 0, sizeof d);
	d.target_faces = 32;
	memset(&r, 0, sizeof r);
	st = qaws_mesh_3d_decimate(&r, m, &d);
	nf = qaws_mesh_3d_get_face_count(m);
	for (f = 0; f < nf; f++)
	{
		unsigned int t[3], left = 0, right = 0;
		qaws_mesh_3d_get_triangle(t, m, f);
		for (i = 0; i < 3; i++)
		{
			double x = (double)qaws_mesh_3d_get_vertex(m, t[i]).x;
			left += x < 0.5 - 1e-12;
			right += x > 0.5 + 1e-12;
		}
		cross += left && right;
	}
	for (i = 0; i < qaws_mesh_3d_get_vertex_count(m); i++)
	{
		qaws_vec3 q = qaws_mesh_3d_get_vertex(m, i);
		ends += fabs((double)q.x - 0.5) < 1e-12 && (q.y == 0 || q.y == 1);
	}
	sprintf(msg, "a UV seam across a grid: %u seam edges made features, %u -> %u faces, faces across the seam %u, its two ends kept %u",
		r.seams, r.faces_before, r.faces_after, cross, ends);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && r.seams == n && r.faces_after <= 40 && cross == 0 && ends == 2, msg);
	free(fv);
	free(uv);
	qaws_mesh_3d_destroy(m);
}

static void test_refusals(void)
{
	qaws_mesh_3d* m = dm_grid(2);
	qaws_decimate_desc d;
	memset(&d, 0, sizeof d);
	TEST_ASSERT(qaws_mesh_3d_decimate(NULL, NULL, NULL) == QAWS_STATUS_INVALID_ARGUMENT, "no mesh refused");
	d.max_error = (qaws_scalar)-1;
	TEST_ASSERT(qaws_mesh_3d_decimate(NULL, m, &d) == QAWS_STATUS_INVALID_ARGUMENT, "a negative error bound refused");
	qaws_mesh_3d_destroy(m);
}

int test_167_mesh_decimate_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 167: decimation\n");
	test_sphere();
	test_plane();
	test_cube();
	test_silhouette();
	test_options();
	test_seams();
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
