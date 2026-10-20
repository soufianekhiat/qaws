/*
 * Test 138: sizing by local feature size (qaws_size_field_create_feature_size)
 *
 *   - a block with a fin 0.02 thick (test 132's slab, held by a block),
 *     remeshed by LOCAL with the field at 3 cells across: the fin's top at
 *     a third of its thickness, at least 2 edges across its thin sides,
 *     angles in 25..120, a valid closed genus 0 mesh, the block's far end
 *     back at the largest size; sizes graded (no edge's sizes further apart
 *     than the gradation allows); built in seconds;
 *   - two plates 0.03 apart: the faces across the gap sized by the gap,
 *     their outer faces by their thickness;
 *   - a sphere of radius 1: twice the radius over the cells, or the
 *     largest size when smaller; combined with the curvature field by min;
 *   - refusals.
 */

#include "test_common.h"
#include "qaws_mesh_field.h"
#include "qaws_remesh.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

static qaws_mesh_3d* fz_make(double const* x, unsigned int nv, unsigned int const* f, unsigned int const* fs, unsigned int nf)
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

/* the profile in (x, z), counter-clockwise: a block [-0.5, -0.2] x [-0.2, 0.2] and a fin [-0.2, 0.15] x [-0.01, 0.01], extruded
   over y in [-0.15, 0.15] */
static double const fz_profile[16] = { -0.5, -0.2, -0.2, -0.2, -0.2, -0.01, 0.15, -0.01, 0.15, 0.01, -0.2, 0.01, -0.2, 0.2, -0.5, 0.2 };

static qaws_mesh_3d* fz_fin(void)
{
	double x[48];
	unsigned int f[32 + 20], fs[12], nf = 0, k, i, n = 0;
	static unsigned int const hex[6] = { 0, 1, 2, 5, 6, 7 }, fin[4] = { 2, 3, 4, 5 };
	for (k = 0; k < 16; k++)
	{
		x[3 * k] = fz_profile[2 * (k % 8)];
		x[3 * k + 1] = k < 8 ? -0.15 : 0.15;
		x[3 * k + 2] = fz_profile[2 * (k % 8) + 1];
	}
	for (k = 0; k < 8; k++)
	{
		f[n++] = k; f[n++] = 8 + k; f[n++] = 8 + (k + 1) % 8; f[n++] = (k + 1) % 8;
		fs[nf++] = 4;
	}
	for (i = 0; i < 6; i++)
		f[n++] = hex[i];
	fs[nf++] = 6;
	for (i = 0; i < 4; i++)
		f[n++] = fin[i];
	fs[nf++] = 4;
	for (i = 0; i < 6; i++)
		f[n++] = 8 + hex[5 - i];
	fs[nf++] = 6;
	for (i = 0; i < 4; i++)
		f[n++] = 8 + fin[3 - i];
	fs[nf++] = 4;
	return fz_make(x, 16, f, fs, nf);
}

/* boxes [lo, hi], one component each */
static qaws_mesh_3d* fz_boxes(double const* lo, double const* hi, unsigned int count)
{
	static unsigned int const f0[24] = { 0, 3, 2, 1, 4, 5, 6, 7, 0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6, 3, 0, 4, 7 };
	double x[3 * 8 * 4];
	unsigned int f[24 * 4], fs[6 * 4], b, k;
	for (b = 0; b < count; b++)
	{
		double const *l = lo + 3 * b, *h = hi + 3 * b;
		for (k = 0; k < 8; k++)
		{
			x[24 * b + 3 * k] = (k == 1 || k == 2 || k == 5 || k == 6) ? h[0] : l[0];
			x[24 * b + 3 * k + 1] = (k == 2 || k == 3 || k == 6 || k == 7) ? h[1] : l[1];
			x[24 * b + 3 * k + 2] = k >= 4 ? h[2] : l[2];
		}
		for (k = 0; k < 24; k++)
			f[24 * b + k] = 8 * b + f0[k];
		for (k = 0; k < 6; k++)
			fs[6 * b + k] = 4;
	}
	return fz_make(x, 8 * count, f, fs, 6 * count);
}

/* the icosahedron split twice on a sphere of radius `radius`, its normals radial */
static qaws_mesh_3d* fz_sphere(double radius, qaws_scalar* normals)
{
	static unsigned int const f0[60] = { 0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11, 1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8,
		3, 9, 4, 3, 4, 2, 3, 2, 6, 3, 6, 8, 3, 8, 9, 4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7, 9, 8, 1 };
	double g = (1 + sqrt(5.0)) / 2, x[3 * 162];
	double const v0[36] = { -1, g, 0, 1, g, 0, -1, -g, 0, 1, -g, 0, 0, -1, g, 0, 1, g, 0, -1, -g, 0, 1, -g, g, 0, -1, g, 0, 1, -g, 0, -1, -g, 0, 1 };
	unsigned int f[3 * 320], old[3 * 320], keys[3 * 480], nv = 12, nf = 20, l, i, j, q, k;
	memcpy(x, v0, sizeof v0);
	memcpy(f, f0, sizeof f0);
	for (l = 0; l < 2; l++)
	{
		unsigned int nk = 0, onf = nf;
		memcpy(old, f, sizeof(unsigned int) * 3 * nf);
		nf = 0;
		for (i = 0; i < onf; i++)
		{
			unsigned int e[3];
			for (j = 0; j < 3; j++)
			{
				unsigned int a = old[3 * i + j], b = old[3 * i + (j + 1) % 3], lo = a < b ? a : b, hi = a < b ? b : a;
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
				unsigned int a = old[3 * i], b = old[3 * i + 1], c = old[3 * i + 2];
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
		{
			x[3 * i + k] /= r;
			normals[3 * i + k] = (qaws_scalar)x[3 * i + k];
			x[3 * i + k] *= radius;
		}
	}
	return fz_make(x, nv, f, NULL, nf);
}

static qaws_size_field* fz_field(qaws_mesh_reference const* r, double cells, double lo, double hi, double grad, double* seconds)
{
	qaws_feature_size_desc d;
	qaws_size_field* sf = NULL;
	clock_t t0 = clock();
	memset(&d, 0, sizeof d);
	d.cells_across = (qaws_scalar)cells;
	d.min_length = (qaws_scalar)lo;
	d.max_length = (qaws_scalar)hi;
	d.gradation = (qaws_scalar)grad;
	qaws_size_field_create_feature_size(&sf, r, &d);
	if (seconds)
		*seconds = (double)(clock() - t0) / CLOCKS_PER_SEC;
	return sf;
}

/* the field at a point, through its closest location on the reference */
static double fz_at(qaws_mesh_reference const* r, qaws_size_field const* sf, double x, double y, double z)
{
	qaws_scalar p[3], s = 0;
	qaws_mesh_location l;
	p[0] = (qaws_scalar)x; p[1] = (qaws_scalar)y; p[2] = (qaws_scalar)z;
	qaws_mesh_reference_project(&l, NULL, r, p, NULL, 1);
	qaws_size_field_evaluate(&s, sf, &l, p, 1);
	return (double)s;
}

/* LOCAL at L or by sf, its angles held in min_angle..110 degrees (0: no bounds) */
static qaws_mesh_3d* fz_remesh(qaws_mesh_3d const* m, qaws_mesh_reference const* r, qaws_size_field const* sf, double L, double min_angle,
	qaws_mesh_report* rep, double* seconds)
{
	qaws_remesh_desc d;
	qaws_mesh_3d* out = NULL;
	clock_t t0 = clock();
	memset(&d, 0, sizeof d);
	d.method = QAWS_REMESH_LOCAL;
	d.reference = r;
	d.target_length = (qaws_scalar)L;
	d.size = sf;
	d.min_angle = (qaws_scalar)(min_angle * M_PI / 180);
	d.max_angle = (qaws_scalar)(min_angle > 0 ? 110 * M_PI / 180 : 0);
	if (qaws_mesh_remesh(&out, rep, m, &d) != QAWS_STATUS_OK)
		out = NULL;
	*seconds = (double)(clock() - t0) / CLOCKS_PER_SEC;
	return out;
}

/* a result's triangles whose centres pass `keep`: their mean edge length, the largest z extent */
typedef struct fz_stats
{
	unsigned int count;
	double mean_edge, max_dz;
} fz_stats;

static void fz_measure(qaws_mesh_3d const* m, int (*keep)(double const* c), fz_stats* s)
{
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), nt = qaws_mesh_3d_get_triangle_count(m), i, k;
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (nt + 1));
	qaws_scalar* p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (nv + 1));
	double sum = 0;
	memset(s, 0, sizeof *s);
	qaws_mesh_3d_get_triangles(t, m);
	qaws_mesh_3d_get_positions(p, m);
	for (i = 0; i < nt; i++)
	{
		double c[3] = { 0, 0, 0 }, zl = HUGE_VAL, zh = -HUGE_VAL;
		for (k = 0; k < 3; k++)
		{
			qaws_scalar const* v = p + 3 * t[3 * i + k];
			c[0] += v[0] / 3.0; c[1] += v[1] / 3.0; c[2] += v[2] / 3.0;
			zl = fmin(zl, (double)v[2]);
			zh = fmax(zh, (double)v[2]);
		}
		if (!keep(c))
			continue;
		for (k = 0; k < 3; k++)
		{
			qaws_scalar const *a = p + 3 * t[3 * i + k], *b = p + 3 * t[3 * i + (k + 1) % 3];
			sum += sqrt((double)((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2])));
		}
		s->count++;
		s->max_dz = fmax(s->max_dz, zh - zl);
	}
	s->mean_edge = s->count ? sum / (3.0 * s->count) : 0;
	free(t);
	free(p);
}

/* the fin's thin sides (y = +-0.15 and x = 0.15 within the fin), its top, the block's far end */
static int fz_fin_side(double const* c)
{
	return c[0] > -0.2 && fabs(c[2]) < 0.0101 && (fabs(fabs(c[1]) - 0.15) < 1e-4 || c[0] > 0.15 - 1e-4);
}

static int fz_fin_top(double const* c)
{
	return c[0] > -0.15 && c[0] < 0.1 && fabs(c[1]) < 0.1 && fabs(fabs(c[2]) - 0.01) < 1e-4;
}

static int fz_far_end(double const* c)
{
	return c[0] < -0.5 + 1e-4;
}

/* the largest growth of the field along the result's edges, over the gradation */
static double fz_growth(qaws_mesh_3d const* m, qaws_size_field const* sf, double beta)
{
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), ne = qaws_mesh_3d_get_edge_count(m), i;
	unsigned int* e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (ne + 1));
	qaws_scalar *p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (nv + 1)), *h = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (nv + 1));
	qaws_mesh_location* l = (qaws_mesh_location*)malloc(sizeof(qaws_mesh_location) * (nv + 1));
	double worst = 0;
	qaws_mesh_3d_get_edges(e, m);
	qaws_mesh_3d_get_positions(p, m);
	qaws_mesh_3d_get_footpoints(l, m);
	qaws_size_field_evaluate(h, sf, l, p, nv);
	for (i = 0; i < ne; i++)
	{
		qaws_scalar const *a = p + 3 * e[2 * i], *b = p + 3 * e[2 * i + 1];
		double d = sqrt((double)((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2])));
		if (d > 0)
			worst = fmax(worst, fabs((double)(h[e[2 * i]] - h[e[2 * i + 1]])) / (beta * d));
	}
	free(e); free(p); free(h); free(l);
	return worst;
}

static void test_fin(void)
{
	qaws_mesh_3d *m = fz_fin(), *out = NULL, *uni = NULL;
	qaws_mesh_reference* r = NULL;
	qaws_size_field* sf = NULL;
	qaws_mesh_report rep, urep;
	qaws_mesh_topology t;
	fz_stats side, uside, top, far;
	double build = 0, sec = 0, usec = 0, at_top, at_far, growth;
	char msg[700];
	qaws_mesh_3d_detect_features(m, (qaws_scalar)(M_PI / 6), 0, 0);
	qaws_mesh_reference_create(&r, m, NULL);
	sf = r ? fz_field(r, 3, 0.002, 0.05, 0.3, &build) : NULL;
	if (!sf)
	{
		TEST_ASSERT(0, "fin: field created");
		goto done;
	}
	at_top = fz_at(r, sf, 0, 0, 0.01);
	at_far = fz_at(r, sf, -0.5, 0, 0);
	printf("  the fin's field: %.5f on its top (0.02 / 3 = %.5f), %.5f at the block's far end, built in %.2f s\n", at_top, 0.02 / 3, at_far,
		build);
	TEST_ASSERT(fabs(at_top - 0.02 / 3) < 0.03 * 0.02 / 3 && fabs(at_far - 0.05) < 1e-9 && build < 5, "fin: the field's sizes, its time");
	out = fz_remesh(m, r, sf, 0.05, 30, &rep, &sec);
	uni = fz_remesh(m, r, NULL, 0.05, 30, &urep, &usec);
	if (!out || !uni)
	{
		TEST_ASSERT(0, "fin: remeshed");
		goto done;
	}
	qaws_mesh_3d_get_topology(&t, out);
	fz_measure(out, fz_fin_side, &side);
	fz_measure(uni, fz_fin_side, &uside);
	fz_measure(out, fz_fin_top, &top);
	fz_measure(out, fz_far_end, &far);
	growth = fz_growth(out, sf, 0.3);
	/* a triangle on a thin side spanning its 0.02: an edge across it alone */
	sprintf(msg, "the fin by LOCAL: %u vertices in %.2f s, valid %d, closed %d, genus %u, angles %.1f..%.1f; thin sides: %u triangles spanning "
		"%.4f at most (uniform at 0.05: %u, %.4f, angles %.1f..%.1f); top: mean edge %.5f; far end: mean edge %.4f; growth along edges %.2f of "
		"the gradation", qaws_mesh_3d_get_vertex_count(out), sec, qaws_mesh_3d_validate(out) == QAWS_STATUS_OK, t.closed, t.genus,
		(double)rep.quality.min_angle, (double)rep.quality.max_angle, side.count, side.max_dz, uside.count, uside.max_dz,
		(double)urep.quality.min_angle, (double)urep.quality.max_angle, top.mean_edge, far.mean_edge, growth);
	printf("  %s\n", msg);
	TEST_ASSERT(qaws_mesh_3d_validate(out) == QAWS_STATUS_OK && t.closed && t.genus == 0 && t.components == 1, msg);
	TEST_ASSERT(side.count > 0 && side.max_dz < 0.02 * (1 - 1e-3) && uside.max_dz > 0.02 * (1 - 1e-3),
		"fin: at least cells_across - 1 edges across the thin sides (one alone with uniform sizes)");
	TEST_ASSERT(rep.quality.min_angle >= 25 && rep.quality.max_angle <= 120, "fin: angles in 25..120");
	TEST_ASSERT(fabs(top.mean_edge - 0.02 / 3) < 0.2 * 0.02 / 3 && far.count > 0 && far.mean_edge > 0.75 * 0.05, "fin: sizes on the top, at the far end");
	TEST_ASSERT(growth <= 1.25, "fin: graded");
done:
	qaws_mesh_3d_destroy(out);
	qaws_mesh_3d_destroy(uni);
	qaws_size_field_destroy(sf);
	qaws_mesh_reference_destroy(r);
	qaws_mesh_3d_destroy(m);
}

static int fz_gap_face(double const* c)
{
	return fabs(c[0]) < 0.1 && fabs(c[1]) < 0.1 && fabs(fabs(c[2]) - 0.015) < 1e-4;
}

static int fz_outer_face(double const* c)
{
	return fabs(c[0]) < 0.1 && fabs(c[1]) < 0.1 && fabs(fabs(c[2]) - 0.08) < 1e-4;
}

static void test_gap(void)
{
	static double const lo[6] = { -0.15, -0.15, -0.08, -0.15, -0.15, 0.015 }, hi[6] = { 0.15, 0.15, -0.015, 0.15, 0.15, 0.08 };
	qaws_mesh_3d *m = fz_boxes(lo, hi, 2), *out = NULL;
	qaws_mesh_reference* r = NULL;
	qaws_size_field* sf = NULL;
	qaws_mesh_report rep;
	fz_stats gap, outer;
	double sec = 0, at_gap, at_outer;
	char msg[400];
	qaws_mesh_3d_detect_features(m, (qaws_scalar)(M_PI / 6), 0, 0);
	qaws_mesh_reference_create(&r, m, NULL);
	sf = r ? fz_field(r, 3, 0.002, 0.05, 0.3, NULL) : NULL;
	if (!sf)
	{
		TEST_ASSERT(0, "gap: field created");
		goto done;
	}
	at_gap = fz_at(r, sf, 0, 0, 0.015);
	at_outer = fz_at(r, sf, 0, 0, 0.08);
	out = fz_remesh(m, r, sf, 0.05, 30, &rep, &sec);
	if (!out)
	{
		TEST_ASSERT(0, "gap: remeshed");
		goto done;
	}
	fz_measure(out, fz_gap_face, &gap);
	fz_measure(out, fz_outer_face, &outer);
	sprintf(msg, "two plates 0.065 thick, 0.03 apart: the field %.5f across the gap (0.01), %.5f outside (0.065 / 3 = %.5f); by LOCAL in %.2f s: "
		"mean edge %.5f on the gap's faces, %.5f outside, angles %.1f..%.1f", at_gap, at_outer, 0.065 / 3, sec, gap.mean_edge, outer.mean_edge,
		(double)rep.quality.min_angle, (double)rep.quality.max_angle);
	printf("  %s\n", msg);
	TEST_ASSERT(fabs(at_gap - 0.01) < 0.03 * 0.01 && fabs(at_outer - 0.065 / 3) < 0.03 * 0.065 / 3, msg);
	TEST_ASSERT(gap.count > 0 && fabs(gap.mean_edge - 0.01) < 0.2 * 0.01 && fabs(outer.mean_edge - 0.065 / 3) < 0.2 * 0.065 / 3
		&& qaws_mesh_3d_validate(out) == QAWS_STATUS_OK, "gap: the faces across sized by the gap");
done:
	qaws_mesh_3d_destroy(out);
	qaws_size_field_destroy(sf);
	qaws_mesh_reference_destroy(r);
	qaws_mesh_3d_destroy(m);
}

static void test_sphere(void)
{
	qaws_scalar normals[3 * 162], *pos = NULL, *h = NULL, *hc = NULL, *hm = NULL, *hk = NULL;
	qaws_mesh_3d* m = fz_sphere(1, normals);
	qaws_mesh_reference_desc rd;
	qaws_mesh_reference* r = NULL;
	qaws_size_field *sf = NULL, *capped = NULL, *curv = NULL, *both = NULL;
	qaws_curvature_size_desc cd;
	qaws_mesh_location* l = NULL;
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), i;
	double lo = HUGE_VAL, up = 0, clo = HUGE_VAL, cup = 0, worst = 0;
	char msg[300];
	memset(&rd, 0, sizeof rd);
	rd.normals = normals;
	qaws_mesh_reference_create(&r, m, &rd);
	sf = r ? fz_field(r, 3, 0.01, 10, 0, NULL) : NULL;
	capped = r ? fz_field(r, 3, 0.01, 0.3, 0.5, NULL) : NULL;
	memset(&cd, 0, sizeof cd);
	cd.angle_tolerance = (qaws_scalar)0.2;
	cd.max_length = 10;
	if (r)
		qaws_size_field_create_curvature(&curv, r, &cd);
	if (sf && curv)
		qaws_size_field_create_min(&both, sf, curv);
	if (!sf || !capped || !both)
	{
		TEST_ASSERT(0, "sphere: fields created");
		goto done;
	}
	pos = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	h = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv);
	hc = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv);
	hm = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv);
	hk = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv);
	l = (qaws_mesh_location*)malloc(sizeof(qaws_mesh_location) * nv);
	qaws_mesh_3d_get_positions(pos, m);
	for (i = 0; i < nv; i++)
		qaws_mesh_reference_get_vertex_location(&l[i], r, i);
	qaws_size_field_evaluate(h, sf, l, pos, nv);
	qaws_size_field_evaluate(hc, capped, l, pos, nv);
	qaws_size_field_evaluate(hm, both, l, pos, nv);
	qaws_size_field_evaluate(hk, curv, l, pos, nv);
	for (i = 0; i < nv; i++)
	{
		lo = fmin(lo, (double)h[i]);
		up = fmax(up, (double)h[i]);
		clo = fmin(clo, (double)hc[i]);
		cup = fmax(cup, (double)hc[i]);
		worst = fmax(worst, fabs((double)hm[i] - fmin((double)h[i], (double)hk[i])));
	}
	sprintf(msg, "a sphere of radius 1 at 3 cells across: %.4f..%.4f (2 / 3); capped at 0.3: %.4f..%.4f; min with the curvature field off by %.1e",
		lo, up, clo, cup, worst);
	printf("  %s\n", msg);
	TEST_ASSERT(lo > 0.97 * 2 / 3.0 && up < 1.03 * 2 / 3.0 && fabs(clo - 0.3) < TEST_BOUND(1e-9, 32) && fabs(cup - 0.3) < TEST_BOUND(1e-9, 32) && worst < TEST_BOUND(1e-12, 32)
		&& qaws_size_field_get_reference(both) == r, msg);
done:
	free(pos); free(h); free(hc); free(hm); free(hk); free(l);
	qaws_size_field_destroy(both);
	qaws_size_field_destroy(curv);
	qaws_size_field_destroy(capped);
	qaws_size_field_destroy(sf);
	qaws_mesh_reference_destroy(r);
	qaws_mesh_3d_destroy(m);
}

static void test_refusals(void)
{
	static double const lo[3] = { 0, 0, 0 }, hi[3] = { 1, 1, 1 };
	qaws_mesh_3d* m = fz_boxes(lo, hi, 1);
	qaws_mesh_reference* r = NULL;
	qaws_feature_size_desc d;
	qaws_size_field* sf = NULL;
	qaws_status a, b, c, e;
	qaws_mesh_reference_create(&r, m, NULL);
	memset(&d, 0, sizeof d);
	d.cells_across = -1;
	a = qaws_size_field_create_feature_size(&sf, r, &d);
	d.cells_across = 3;
	d.min_length = (qaws_scalar)0.5;
	d.max_length = (qaws_scalar)0.1;
	b = qaws_size_field_create_feature_size(&sf, r, &d);
	c = qaws_size_field_create_feature_size(&sf, NULL, NULL);
	e = qaws_size_field_create_feature_size(&sf, r, NULL);
	TEST_ASSERT(a == QAWS_STATUS_INVALID_ARGUMENT && b == QAWS_STATUS_INVALID_ARGUMENT && c == QAWS_STATUS_INVALID_ARGUMENT && e == QAWS_STATUS_OK
		&& sf, "refusals; the defaults");
	qaws_size_field_destroy(sf);
	qaws_mesh_reference_destroy(r);
	qaws_mesh_3d_destroy(m);
}

int test_138_feature_size_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 138: sizing by local feature size\n");
	test_fin();
	test_gap();
	test_sphere();
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
