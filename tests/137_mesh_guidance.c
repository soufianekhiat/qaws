/*
 * Test 137: user guidance of direction fields and quads
 *
 *   - density paint: a sphere painted fine on its upper half, quadrangulated
 *     all quads: the mean edge on each side within 20% of its painted size,
 *     every quad's scaled Jacobian positive;
 *   - directions: a stroke round the equator at 30 degrees to it, pulled
 *     and held, through both direction fields (LOCAL, SMOOTHEST): the field
 *     within a few degrees of it on the stroke; quads aligned there;
 *   - edge loops: a parallel of the sphere given as guide edges: every
 *     segment (sampled) within a small distance of the quads' edges, far
 *     closer than without;
 *   - pinned vertices: quad vertices where they are, flagged corners; kept
 *     by LOCAL and CVT remeshing;
 *   - placed cones: eight +1 at a cube's corner directions on the sphere
 *     (SMOOTHEST: 8 singular faces, each round its vertex; seven refused by
 *     Poincare-Hopf; LOCAL refused), quads a cube-sphere with exactly 8
 *     valence-3 vertices on them (INSTANT, GLOBAL); a +1 / -1 pair on a
 *     square (the field, and quads with valence 3 and 5 on them);
 *   - attractors: a free cone of the smoothest field drawn to a point and
 *     pushed off one (SMOOTHEST, LOCAL); quads' irregular vertex drawn;
 *   - refusals: bad indices, directions on a foreign reference.
 */

#include "test_common.h"
#include "qaws_quadrangulate.h"
#include "qaws_remesh.h"
#include "qaws_mesh_field.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

/* a latitude-longitude sphere: the poles, then nlat - 1 parallels of nlon vertices (vertex 1 + (j - 1) nlon + i) */
static qaws_mesh_3d* gd_sphere(unsigned int nlat, unsigned int nlon, double** out_x, unsigned int* out_nv)
{
	unsigned int nv = 2 + (nlat - 1) * nlon, nf = 2 * nlon * (nlat - 1), i, j, k = 0;
	double* x = (double*)malloc(sizeof(double) * 3 * nv);
	unsigned int* f = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nf);
	qaws_scalar* p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	x[0] = 0; x[1] = 0; x[2] = 1;
	x[3 * (nv - 1)] = 0; x[3 * (nv - 1) + 1] = 0; x[3 * (nv - 1) + 2] = -1;
	for (j = 1; j < nlat; j++)
		for (i = 0; i < nlon; i++)
		{
			double th = M_PI * j / nlat, ph = 2 * M_PI * (i + 0.5 * (j & 1)) / nlon;
			unsigned int v = 1 + (j - 1) * nlon + i;
			x[3 * v] = sin(th) * cos(ph);
			x[3 * v + 1] = sin(th) * sin(ph);
			x[3 * v + 2] = cos(th);
		}
	for (i = 0; i < nlon; i++)
	{
		unsigned int a = 1 + i, b = 1 + (i + 1) % nlon, c = 1 + (nlat - 2) * nlon + i, e = 1 + (nlat - 2) * nlon + (i + 1) % nlon;
		f[k++] = 0; f[k++] = a; f[k++] = b;
		f[k++] = nv - 1; f[k++] = e; f[k++] = c;
	}
	for (j = 1; j + 1 < nlat; j++)
		for (i = 0; i < nlon; i++)
		{
			unsigned int a = 1 + (j - 1) * nlon + i, b = 1 + (j - 1) * nlon + (i + 1) % nlon, c = 1 + j * nlon + i, e = 1 + j * nlon + (i + 1) % nlon;
			if (j & 1)
			{
				f[k++] = a; f[k++] = c; f[k++] = e;
				f[k++] = a; f[k++] = e; f[k++] = b;
			}
			else
			{
				f[k++] = a; f[k++] = c; f[k++] = b;
				f[k++] = b; f[k++] = c; f[k++] = e;
			}
		}
	for (i = 0; i < 3 * nv; i++)
		p[i] = (qaws_scalar)x[i];
	memset(&d, 0, sizeof d);
	d.positions = p;
	d.vertex_count = nv;
	d.face_vertices = f;
	d.face_count = nf;
	qaws_mesh_3d_create(&d, &m);
	free(p);
	free(f);
	*out_x = x;
	*out_nv = nv;
	return m;
}

/* the stroke's direction at a point of the unit sphere: 30 degrees from the parallel toward the north */
static void gd_stroke(double const* p, double* d)
{
	double r = sqrt(p[0] * p[0] + p[1] * p[1]), e[3], nth[3];
	unsigned int k;
	e[0] = -p[1] / r; e[1] = p[0] / r; e[2] = 0;
	nth[0] = -p[0] * p[2] / r; nth[1] = -p[1] * p[2] / r; nth[2] = r;
	for (k = 0; k < 3; k++)
		d[k] = cos(M_PI / 6) * e[k] + sin(M_PI / 6) * nth[k];
}

/* the angle between tangent q and the cross of d at the point of unit normal n, in degrees, 0 to 45 */
static double gd_cross_angle(double const* q, double const* d, double const* n)
{
	double t[3], a, b, g;
	t[0] = n[1] * d[2] - n[2] * d[1];
	t[1] = n[2] * d[0] - n[0] * d[2];
	t[2] = n[0] * d[1] - n[1] * d[0];
	a = q[0] * d[0] + q[1] * d[1] + q[2] * d[2];
	b = q[0] * t[0] + q[1] * t[1] + q[2] * t[2];
	g = fmod(fabs(atan2(b, a)) * 180 / M_PI, 90);
	return g > 45 ? 90 - g : g;
}

/* a quad's smallest scaled Jacobian about its Newell normal */
static double gd_quad_sj(double const* p, unsigned int const* v)
{
	double n[3] = { 0, 0, 0 }, ln, worst = 1;
	unsigned int k, i;
	for (k = 0; k < 4; k++)
	{
		double const *a = p + 3 * v[k], *b = p + 3 * v[(k + 1) & 3];
		n[0] += (a[1] - b[1]) * (a[2] + b[2]);
		n[1] += (a[2] - b[2]) * (a[0] + b[0]);
		n[2] += (a[0] - b[0]) * (a[1] + b[1]);
	}
	ln = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
	for (k = 0; k < 4; k++)
	{
		double const *o = p + 3 * v[k], *a = p + 3 * v[(k + 1) & 3], *b = p + 3 * v[(k + 3) & 3];
		double e1[3], e2[3], l1, l2;
		for (i = 0; i < 3; i++)
		{
			e1[i] = a[i] - o[i];
			e2[i] = b[i] - o[i];
		}
		l1 = sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
		l2 = sqrt(e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2]);
		worst = fmin(worst, l1 > 0 && l2 > 0 && ln > 0 ? ((e1[1] * e2[2] - e1[2] * e2[1]) * n[0] + (e1[2] * e2[0] - e1[0] * e2[2]) * n[1]
			+ (e1[0] * e2[1] - e1[1] * e2[0]) * n[2]) / (ln * l1 * l2) : -1);
	}
	return worst;
}

/* a mesh's arrays as doubles (caller frees x, fv, fs, e) */
typedef struct gd_arrays
{
	double* x;
	unsigned int nv, nf, ne, *fv, *fs, *e;
} gd_arrays;

static void gd_get(qaws_mesh_3d const* m, gd_arrays* a)
{
	qaws_scalar* p;
	unsigned int i;
	a->nv = qaws_mesh_3d_get_vertex_count(m);
	a->nf = qaws_mesh_3d_get_face_count(m);
	a->ne = qaws_mesh_3d_get_edge_count(m);
	a->x = (double*)malloc(sizeof(double) * 3 * (a->nv + 1));
	a->fv = (unsigned int*)malloc(sizeof(unsigned int) * (qaws_mesh_3d_get_corner_count(m) + 1));
	a->fs = (unsigned int*)malloc(sizeof(unsigned int) * (a->nf + 1));
	a->e = (unsigned int*)malloc(sizeof(unsigned int) * (2 * a->ne + 2));
	p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (a->nv + 1));
	qaws_mesh_3d_get_positions(m, p);
	for (i = 0; i < 3 * a->nv; i++)
		a->x[i] = (double)p[i];
	free(p);
	qaws_mesh_3d_get_faces(m, a->fv, a->fs);
	qaws_mesh_3d_get_edges(m, a->e);
}

static void gd_free(gd_arrays* a)
{
	free(a->x); free(a->fv); free(a->fs); free(a->e);
}

/* all quads, the smallest scaled Jacobian */
static double gd_quads(gd_arrays const* a, unsigned int* nonquads)
{
	unsigned int f, c0;
	double s = 1;
	*nonquads = 0;
	for (f = 0, c0 = 0; f < a->nf; c0 += a->fs[f], f++)
	{
		if (a->fs[f] != 4)
		{
			(*nonquads)++;
			continue;
		}
		s = fmin(s, gd_quad_sj(a->x, a->fv + c0));
	}
	return s;
}

static qaws_mesh_3d* gd_quadrangulate(qaws_mesh_3d const* m, qaws_quadrangulate_desc* d, qaws_mesh_report* rep, qaws_status* st)
{
	qaws_mesh_3d* out = NULL;
	*st = qaws_mesh_quadrangulate(m, d, &out, rep);
	return out;
}

static void test_density(void)
{
	double* x0;
	unsigned int nv0, i;
	qaws_mesh_3d *s = gd_sphere(24, 48, &x0, &nv0), *q = NULL;
	qaws_mesh_reference* r = NULL;
	qaws_size_field* sf = NULL;
	qaws_scalar* sizes = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv0);
	qaws_quadrangulate_desc d;
	qaws_mesh_report rep;
	qaws_status st = QAWS_STATUS_INTERNAL_ERROR;
	gd_arrays a;
	double sum[2] = { 0, 0 }, sj = -1;
	unsigned int cnt[2] = { 0, 0 }, nonq = 1;
	char msg[400];
	memset(&a, 0, sizeof a);
	for (i = 0; i < nv0; i++)
		sizes[i] = (qaws_scalar)(x0[3 * i + 2] > 0 ? 0.08 : 0.16);
	qaws_mesh_reference_create(s, NULL, &r);
	qaws_size_field_create_vertices(r, sizes, &sf);
	memset(&d, 0, sizeof d);
	d.reference = r;
	d.size = sf;
	d.pure = 1;
	q = gd_quadrangulate(s, &d, &rep, &st);
	if (q)
	{
		gd_get(q, &a);
		for (i = 0; i < a.ne; i++)
		{
			double const *p = a.x + 3 * a.e[2 * i], *u = a.x + 3 * a.e[2 * i + 1];
			double z = 0.5 * (p[2] + u[2]), l = sqrt((p[0] - u[0]) * (p[0] - u[0]) + (p[1] - u[1]) * (p[1] - u[1]) + (p[2] - u[2]) * (p[2] - u[2]));
			if (fabs(z) < 0.3)
				continue;   /* (the transition) */
			sum[z < 0] += l;
			cnt[z < 0]++;
		}
		sj = gd_quads(&a, &nonq);
	}
	sprintf(msg, "density paint 0.08 above, 0.16 below: %s, mean edge above %.4f (%u), below %.4f (%u), %u faces not quads, SJ min %.3f, %u singular",
		qaws_status_to_string(st), cnt[0] ? sum[0] / cnt[0] : 0, cnt[0], cnt[1] ? sum[1] / cnt[1] : 0, cnt[1], nonq, sj, rep.quality.irregular);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && cnt[0] && cnt[1] && fabs(sum[0] / cnt[0] - 0.08) < 0.2 * 0.08 && fabs(sum[1] / cnt[1] - 0.16) < 0.2 * 0.16
		&& nonq == 0 && sj > 0, msg);
	gd_free(&a);
	qaws_mesh_3d_destroy(q);
	qaws_size_field_destroy(sf);
	qaws_mesh_reference_destroy(r);
	qaws_mesh_3d_destroy(s);
	free(sizes);
	free(x0);
}

static void test_directions(void)
{
	double* x0;
	unsigned int nv0, i, k, j;
	qaws_mesh_3d* s = gd_sphere(24, 48, &x0, &nv0);
	qaws_scalar *g = (qaws_scalar*)calloc(3 * nv0, sizeof(qaws_scalar)), *w = (qaws_scalar*)calloc(nv0, sizeof(qaws_scalar));
	qaws_scalar* q = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv0);
	double mean[4], worst[4], before = 0;
	double devs[2] = { 90, 0 };
	unsigned int ns = 0;
	char msg[500];
	for (i = 0; i < nv0; i++)
		if (fabs(x0[3 * i + 2]) < 0.2)
		{
			double d[3];
			gd_stroke(x0 + 3 * i, d);
			for (k = 0; k < 3; k++)
				g[3 * i + k] = (qaws_scalar)d[k];
			ns++;
		}
	/* LOCAL and SMOOTHEST, pulled (weight 1) and held (-1) */
	for (j = 0; j < 5; j++)
	{
		qaws_direction_field_desc dd;
		memset(&dd, 0, sizeof dd);
		dd.method = j & 2 ? QAWS_DIRECTION_FIELD_SMOOTHEST : QAWS_DIRECTION_FIELD_LOCAL;
		if (j < 4)
		{
			for (i = 0; i < nv0; i++)
				w[i] = (qaws_scalar)(g[3 * i] != 0 || g[3 * i + 1] != 0 || g[3 * i + 2] != 0 ? (j & 1 ? -1 : 1) : 0);
			dd.guide_directions = g;
			dd.guide_weights = w;
		}
		qaws_mesh_3d_direction_field(s, &dd, q, NULL);
		mean[j & 3] = 0;
		worst[j & 3] = 0;
		for (i = 0; i < nv0; i++)
		{
			double d[3], qq[3];
			if (fabs(x0[3 * i + 2]) >= 0.2)
				continue;
			gd_stroke(x0 + 3 * i, d);
			for (k = 0; k < 3; k++)
				qq[k] = (double)q[3 * i + k];
			if (j == 4)
				before += gd_cross_angle(qq, d, x0 + 3 * i) / ns;
			else
			{
				double a = gd_cross_angle(qq, d, x0 + 3 * i);
				mean[j] += a / ns;
				worst[j] = fmax(worst[j], a);
			}
		}
	}
	sprintf(msg, "a stroke at 30 degrees on %u vertices: LOCAL pulled mean %.2f max %.2f, held %.2f %.2f; SMOOTHEST pulled %.2f %.2f, held %.2f %.2f; unguided mean %.1f",
		ns, mean[0], worst[0], mean[1], worst[1], mean[2], worst[2], mean[3], worst[3], before);
	printf("  %s\n", msg);
	TEST_ASSERT(mean[0] < 3 && worst[0] < 8 && mean[1] < 0.5 && worst[1] < 1 && mean[2] < 3 && worst[2] < 8 && mean[3] < 0.5 && worst[3] < 1
		&& before > 5, msg);
	/* quads: the edges on the stroke along it */
	for (j = 0; j < 2; j++)
	{
		qaws_quadrangulate_desc d;
		qaws_mesh_report rep;
		qaws_status st;
		qaws_mesh_3d* m;
		gd_arrays a;
		double dev = 0;
		unsigned int n = 0;
		for (i = 0; i < nv0; i++)
			w[i] = (qaws_scalar)(g[3 * i] != 0 || g[3 * i + 1] != 0 || g[3 * i + 2] != 0 ? 1 : 0);
		memset(&d, 0, sizeof d);
		d.target_length = (qaws_scalar)0.12;
		d.pure = 1;
		if (j == 0)
		{
			d.guide_directions = g;
			d.guide_weights = w;
		}
		m = gd_quadrangulate(s, &d, &rep, &st);
		if (!m)
		{
			TEST_ASSERT(0, "stroke: quadrangulated");
			continue;
		}
		gd_get(m, &a);
		for (i = 0; i < a.ne; i++)
		{
			double const *p = a.x + 3 * a.e[2 * i], *u = a.x + 3 * a.e[2 * i + 1];
			double mid[3], t[3], dd[3], l, dn;
			for (k = 0; k < 3; k++)
			{
				mid[k] = 0.5 * (p[k] + u[k]);
				t[k] = u[k] - p[k];
			}
			l = sqrt(mid[0] * mid[0] + mid[1] * mid[1] + mid[2] * mid[2]);
			for (k = 0; k < 3; k++)
				mid[k] /= l;
			if (fabs(mid[2]) > 0.15)
				continue;
			dn = t[0] * mid[0] + t[1] * mid[1] + t[2] * mid[2];
			for (k = 0; k < 3; k++)
				t[k] -= dn * mid[k];
			gd_stroke(mid, dd);
			dev += gd_cross_angle(t, dd, mid);
			n++;
		}
		dev = n ? dev / n : 90;
		sprintf(msg, "quads %s the stroke: %s, %u edges on it, mean angle to its crosses %.2f degrees, %u faces, %u singular",
			j ? "without" : "with", qaws_status_to_string(st), n, dev, a.nf, rep.quality.irregular);
		printf("  %s\n", msg);
		devs[j] = dev;
		TEST_ASSERT(st == QAWS_STATUS_OK && n > 50, msg);
		gd_free(&a);
		qaws_mesh_3d_destroy(m);
	}
	sprintf(msg, "quads on the stroke: %.2f degrees off it with it, %.2f without", devs[0], devs[1]);
	TEST_ASSERT(devs[0] < 10 && devs[0] < 0.6 * devs[1], msg);
	free(g); free(w); free(q); free(x0);
	qaws_mesh_3d_destroy(s);
}

/* the distance from p to the nearest edge of a */
static double gd_edge_distance(gd_arrays const* a, double const* p)
{
	unsigned int i, k;
	double best = HUGE_VAL;
	for (i = 0; i < a->ne; i++)
	{
		double const *u = a->x + 3 * a->e[2 * i], *v = a->x + 3 * a->e[2 * i + 1];
		double uv[3], t, d2 = 0, l2 = 0;
		for (k = 0; k < 3; k++)
		{
			uv[k] = v[k] - u[k];
			l2 += uv[k] * uv[k];
		}
		t = l2 > 0 ? ((p[0] - u[0]) * uv[0] + (p[1] - u[1]) * uv[1] + (p[2] - u[2]) * uv[2]) / l2 : 0;
		t = t < 0 ? 0 : t > 1 ? 1 : t;
		for (k = 0; k < 3; k++)
			d2 += (u[k] + t * uv[k] - p[k]) * (u[k] + t * uv[k] - p[k]);
		best = fmin(best, d2);
	}
	return sqrt(best);
}

static void test_edge_loops(void)
{
	double* x0;
	unsigned int nv0, i, j, nlon = 48, ring = 8;   /* (the parallel at 60 degrees from the north pole) */
	qaws_mesh_3d* s = gd_sphere(24, nlon, &x0, &nv0);
	unsigned int* ge = (unsigned int*)malloc(sizeof(unsigned int) * 2 * nlon);
	double worst[2] = { 0, 0 }, mean[2] = { 0, 0 }, L = 0.12;
	char msg[400];
	for (i = 0; i < nlon; i++)
	{
		ge[2 * i] = 1 + (ring - 1) * nlon + i;
		ge[2 * i + 1] = 1 + (ring - 1) * nlon + (i + 1) % nlon;
	}
	for (j = 0; j < 2; j++)
	{
		qaws_quadrangulate_desc d;
		qaws_mesh_report rep;
		qaws_status st;
		qaws_mesh_3d* m;
		gd_arrays a;
		unsigned int nq;
		double sj;
		memset(&d, 0, sizeof d);
		d.target_length = (qaws_scalar)L;
		d.pure = 1;
		if (j == 0)
		{
			d.guide_edges = ge;
			d.guide_edge_count = nlon;
		}
		m = gd_quadrangulate(s, &d, &rep, &st);
		if (!m)
		{
			TEST_ASSERT(0, "edge loop: quadrangulated");
			continue;
		}
		gd_get(m, &a);
		sj = gd_quads(&a, &nq);
		/* each segment sampled at a quarter, a half and three quarters, on the sphere */
		for (i = 0; i < nlon; i++)
		{
			double const *pa = x0 + 3 * ge[2 * i], *pb = x0 + 3 * ge[2 * i + 1];
			unsigned int k, t;
			for (t = 1; t < 4; t++)
			{
				double p[3], l, dd;
				for (k = 0; k < 3; k++)
					p[k] = pa[k] + 0.25 * t * (pb[k] - pa[k]);
				l = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
				for (k = 0; k < 3; k++)
					p[k] /= l;
				dd = gd_edge_distance(&a, p);
				worst[j] = fmax(worst[j], dd);
				mean[j] += dd / (3 * nlon);
			}
		}
		sprintf(msg, "edge loop on a parallel %s guide edges: %s, segments to the quads' edges mean %.4f max %.4f (target %.2f), %u faces not quads, SJ min %.3f, %u singular",
			j ? "without" : "with", qaws_status_to_string(st), mean[j], worst[j], L, nq, sj, rep.quality.irregular);
		printf("  %s\n", msg);
		TEST_ASSERT(st == QAWS_STATUS_OK && nq == 0 && sj > 0 && (j || worst[0] < 0.08 * L), msg);
		gd_free(&a);
		qaws_mesh_3d_destroy(m);
	}
	sprintf(msg, "edge loop: max distance %.4f with the guide against %.4f without", worst[0], worst[1]);
	TEST_ASSERT(worst[1] > 4 * worst[0], msg);
	free(ge); free(x0);
	qaws_mesh_3d_destroy(s);
}

static void test_pinned(void)
{
	double* x0;
	unsigned int nv0, i, j, k;
	qaws_mesh_3d* s = gd_sphere(24, 48, &x0, &nv0);
	unsigned int pins[4] = { 101, 377, 610, 845 };
	double worst[3] = { 0, 0, 0 };
	unsigned int corner[3] = { 0, 0, 0 };
	qaws_status sts[3];
	char msg[400];
	for (j = 0; j < 3; j++)
	{
		qaws_mesh_3d* m = NULL;
		gd_arrays a;
		if (j == 0)
		{
			qaws_quadrangulate_desc d;
			memset(&d, 0, sizeof d);
			d.target_length = (qaws_scalar)0.12;
			d.pure = 1;
			d.pinned = pins;
			d.pinned_count = 4;
			sts[j] = qaws_mesh_quadrangulate(s, &d, &m, NULL);
		}
		else
		{
			qaws_remesh_desc d;
			memset(&d, 0, sizeof d);
			d.method = j == 1 ? QAWS_REMESH_LOCAL : QAWS_REMESH_CVT;
			d.target_length = (qaws_scalar)0.09;
			d.pinned = pins;
			d.pinned_count = 4;
			sts[j] = qaws_mesh_remesh(s, &d, &m, NULL);
		}
		if (!m)
		{
			worst[j] = HUGE_VAL;
			continue;
		}
		gd_get(m, &a);
		for (i = 0; i < 4; i++)
		{
			double best = HUGE_VAL;
			unsigned int bv = 0;
			for (k = 0; k < a.nv; k++)
			{
				double const* p = a.x + 3 * k;
				double d2 = (p[0] - x0[3 * pins[i]]) * (p[0] - x0[3 * pins[i]]) + (p[1] - x0[3 * pins[i] + 1]) * (p[1] - x0[3 * pins[i] + 1])
					+ (p[2] - x0[3 * pins[i] + 2]) * (p[2] - x0[3 * pins[i] + 2]);
				if (d2 < best)
				{
					best = d2;
					bv = k;
				}
			}
			worst[j] = fmax(worst[j], sqrt(best));
			corner[j] += qaws_mesh_3d_is_feature_corner(m, bv) != 0;
		}
		gd_free(&a);
		qaws_mesh_3d_destroy(m);
	}
	sprintf(msg, "pinned: quads %s, nearest vertex %.2g, corners %u of 4; LOCAL %s %.2g %u; CVT %s %.2g %u", qaws_status_to_string(sts[0]), worst[0],
		corner[0], qaws_status_to_string(sts[1]), worst[1], corner[1], qaws_status_to_string(sts[2]), worst[2], corner[2]);
	printf("  %s\n", msg);
	TEST_ASSERT(sts[0] == QAWS_STATUS_OK && sts[1] == QAWS_STATUS_OK && sts[2] == QAWS_STATUS_OK && worst[0] < 1e-12 && worst[1] < 1e-12
		&& worst[2] < 1e-12 && corner[0] == 4 && corner[1] == 4 && corner[2] == 4, msg);
	free(x0);
	qaws_mesh_3d_destroy(s);
}


/* a square [-1, 1]^2 at z = 0, an (n + 1)^2 grid split along alternating diagonals (vertex j (n + 1) + i) */
static qaws_mesh_3d* gd_square(unsigned int n, double** out_x, unsigned int* out_nv)
{
	unsigned int nv = (n + 1) * (n + 1), nf = 2 * n * n, i, j, k = 0;
	double* x = (double*)malloc(sizeof(double) * 3 * nv);
	unsigned int* f = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nf);
	qaws_scalar* p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	for (j = 0; j <= n; j++)
		for (i = 0; i <= n; i++)
		{
			x[3 * (j * (n + 1) + i)] = -1 + 2.0 * i / n;
			x[3 * (j * (n + 1) + i) + 1] = -1 + 2.0 * j / n;
			x[3 * (j * (n + 1) + i) + 2] = 0;
		}
	for (j = 0; j < n; j++)
		for (i = 0; i < n; i++)
		{
			unsigned int a = j * (n + 1) + i, b = a + 1, c = a + n + 1, e = c + 1;
			if ((i + j) & 1)
			{
				f[k++] = a; f[k++] = b; f[k++] = e;
				f[k++] = a; f[k++] = e; f[k++] = c;
			}
			else
			{
				f[k++] = a; f[k++] = b; f[k++] = c;
				f[k++] = b; f[k++] = e; f[k++] = c;
			}
		}
	for (i = 0; i < 3 * nv; i++)
		p[i] = (qaws_scalar)x[i];
	memset(&d, 0, sizeof d);
	d.positions = p;
	d.vertex_count = nv;
	d.face_vertices = f;
	d.face_count = nf;
	qaws_mesh_3d_create(&d, &m);
	free(p);
	free(f);
	*out_x = x;
	*out_nv = nv;
	return m;
}

/* the vertex of x (nv of them) nearest the point p */
static unsigned int gd_nearest(double const* x, unsigned int nv, double const* p)
{
	unsigned int i, best = 0;
	double bd = HUGE_VAL;
	for (i = 0; i < nv; i++)
	{
		double a = x[3 * i] - p[0], b = x[3 * i + 1] - p[1], c = x[3 * i + 2] - p[2];
		if (a * a + b * b + c * c < bd)
		{
			bd = a * a + b * b + c * c;
			best = i;
		}
	}
	return best;
}

/* a field's singular faces: their count; at_v: how many have index want and v among their vertices */
static unsigned int gd_singular(gd_arrays const* a, int const* idx, unsigned int v, int want, unsigned int* at_v)
{
	unsigned int f, c0, n = 0, i;
	*at_v = 0;
	for (f = 0, c0 = 0; f < a->nf; c0 += a->fs[f], f++)
	{
		if (!idx[f])
			continue;
		n++;
		for (i = 0; i < a->fs[f]; i++)
			*at_v += a->fv[c0 + i] == v && idx[f] == want;
	}
	return n;
}

/* the singular face whose centroid is nearest p: its distance, its centroid (may be NULL) */
static double gd_cone_distance(gd_arrays const* a, int const* idx, double const* p, double* out_c)
{
	unsigned int f, c0, i, k;
	double best = HUGE_VAL;
	for (f = 0, c0 = 0; f < a->nf; c0 += a->fs[f], f++)
	{
		double c[3] = { 0, 0, 0 }, d2 = 0;
		if (!idx[f])
			continue;
		for (i = 0; i < a->fs[f]; i++)
			for (k = 0; k < 3; k++)
				c[k] += a->x[3 * a->fv[c0 + i] + k] / a->fs[f];
		for (k = 0; k < 3; k++)
			d2 += (c[k] - p[k]) * (c[k] - p[k]);
		if (sqrt(d2) < best)
		{
			best = sqrt(d2);
			if (out_c)
				memcpy(out_c, c, sizeof c);
		}
	}
	return best;
}

/* the quads' vertices off the border whose valence is not 4 (counted; the distance from to to the nearest, when to is given,
   lowered into out_near); per cone (3 each) the distance to the nearest vertex
   (the worst kept) and that vertex's valence */
static unsigned int gd_irregular(gd_arrays const* a, double const* cones, unsigned int ncone, double* worst, unsigned int* val, double const* to,
	double* out_near)
{
	unsigned int *deg = (unsigned int*)calloc(a->nv, sizeof(unsigned int)), *edges = (unsigned int*)calloc(a->nv, sizeof(unsigned int));
	unsigned int i, j, k, n = 0, c0, f;
	/* faces and edges round each vertex: as many off the border, one edge more on it */
	for (f = 0, c0 = 0; f < a->nf; c0 += a->fs[f], f++)
		for (i = 0; i < a->fs[f]; i++)
			deg[a->fv[c0 + i]]++;
	for (i = 0; i < a->ne; i++)
	{
		edges[a->e[2 * i]]++;
		edges[a->e[2 * i + 1]]++;
	}
	for (i = 0; i < a->nv; i++)
		if (edges[i] == deg[i] && deg[i] != 4)
		{
			n++;
			if (to)
				*out_near = fmin(*out_near, sqrt((a->x[3 * i] - to[0]) * (a->x[3 * i] - to[0]) + (a->x[3 * i + 1] - to[1]) * (a->x[3 * i + 1] - to[1])
					+ (a->x[3 * i + 2] - to[2]) * (a->x[3 * i + 2] - to[2])));
		}
	*worst = 0;
	for (j = 0; j < ncone; j++)
	{
		unsigned int bv = gd_nearest(a->x, a->nv, cones + 3 * j);
		double d2 = 0;
		for (k = 0; k < 3; k++)
			d2 += (a->x[3 * bv + k] - cones[3 * j + k]) * (a->x[3 * bv + k] - cones[3 * j + k]);
		*worst = fmax(*worst, sqrt(d2));
		val[j] = deg[bv];
	}
	free(deg);
	free(edges);
	return n;
}

static void test_cones(void)
{
	double *x0, *xs, cube[24], pair[6], worst, sj;
	unsigned int nv0, nvs, i, k, j, cv[8], pv[2], n = 0, at, at2, val[8], corners, nonq;
	int ci[8] = { 1, 1, 1, 1, 1, 1, 1, 1 }, pi[2] = { 1, -1 }, *idx;
	qaws_mesh_3d *s = gd_sphere(24, 48, &x0, &nv0), *sq = gd_square(24, &xs, &nvs), *q = NULL;
	qaws_direction_field_desc dd;
	qaws_quadrangulate_desc qd;
	qaws_mesh_report rep;
	qaws_status st, s2, s3;
	gd_arrays a;
	char msg[600];
	/* the sphere: eight +1 cones at a cube's corner directions */
	for (i = 0; i < 8; i++)
	{
		double d[3];
		d[0] = (i & 1 ? 1 : -1) / sqrt(3.0);
		d[1] = (i & 2 ? 1 : -1) / sqrt(3.0);
		d[2] = (i & 4 ? 1 : -1) / sqrt(3.0);
		cv[i] = gd_nearest(x0, nv0, d);
		for (k = 0; k < 3; k++)
			cube[3 * i + k] = x0[3 * cv[i] + k];
	}
	gd_get(s, &a);
	idx = (int*)malloc(sizeof(int) * a.nf);
	memset(&dd, 0, sizeof dd);
	dd.method = QAWS_DIRECTION_FIELD_SMOOTHEST;
	dd.cone_vertices = cv;
	dd.cone_indices = ci;
	dd.cone_count = 8;
	st = qaws_mesh_3d_direction_field(s, &dd, NULL, idx);
	for (i = 0, at = 0; i < 8; i++)
	{
		n = gd_singular(&a, idx, cv[i], 1, &at2);
		at += at2 == 1;
	}
	dd.cone_count = 7;
	s2 = qaws_mesh_3d_direction_field(s, &dd, NULL, idx);
	dd.cone_count = 8;
	dd.method = QAWS_DIRECTION_FIELD_LOCAL;
	s3 = qaws_mesh_3d_direction_field(s, &dd, NULL, idx);
	sprintf(msg, "sphere, 8 cones placed (SMOOTHEST): %s, %u singular faces, %u of 8 cones +1 on a face round their vertex; 7 cones %s; LOCAL %s",
		qaws_status_to_string(st), n, at, qaws_status_to_string(s2), qaws_status_to_string(s3));
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && n == 8 && at == 8 && s2 == QAWS_STATUS_INVALID_ARGUMENT && s3 == QAWS_STATUS_UNSUPPORTED_OPERATION, msg);
	free(idx);
	gd_free(&a);
	/* quads: a cube-sphere, its eight valence-3 vertices exactly on the cones (INSTANT, GLOBAL) */
	for (j = 0; j < 2; j++)
	{
		memset(&qd, 0, sizeof qd);
		qd.method = j ? QAWS_QUADS_GLOBAL : QAWS_QUADS_INSTANT;
		qd.target_length = (qaws_scalar)0.15;
		qd.pure = 1;
		qd.cone_vertices = cv;
		qd.cone_indices = ci;
		qd.cone_count = 8;
		memset(&rep, 0, sizeof rep);
		q = gd_quadrangulate(s, &qd, &rep, &st);
		memset(&a, 0, sizeof a);
		n = 0;
		worst = HUGE_VAL;
		at = 0;
		corners = 0;
		nonq = 1;
		sj = -1;
		if (q)
		{
			gd_get(q, &a);
			n = gd_irregular(&a, cube, 8, &worst, val, NULL, NULL);
			sj = gd_quads(&a, &nonq);
			for (i = 0; i < 8; i++)
				at += val[i] == 3;
			for (i = 0; i < a.nv; i++)
				corners += qaws_mesh_3d_is_feature_corner(q, i) != 0;
			sprintf(msg, "cube-sphere %s: %s, %u faces (%u not quads), SJ min %.3f, %u irregular; cones: nearest vertex %.2g, %u of 8 of valence 3; %u corners; map %u",
				j ? "GLOBAL" : "INSTANT", qaws_status_to_string(st), a.nf, nonq, sj, n, worst, at, corners, (unsigned int)rep.iterations);
		}
		else
			sprintf(msg, "cube-sphere %s: %s", j ? "GLOBAL" : "INSTANT", qaws_status_to_string(st));
		printf("  %s\n", msg);
		TEST_ASSERT(st == QAWS_STATUS_OK && n == 8 && worst < 1e-12 && at == 8 && corners == 0 && nonq == 0 && sj > 0.3, msg);
		gd_free(&a);
		qaws_mesh_3d_destroy(q);
	}
	memset(&qd, 0, sizeof qd);
	qd.cone_vertices = cv;
	qd.cone_indices = ci;
	qd.cone_count = 7;
	q = gd_quadrangulate(s, &qd, &rep, &st);
	sprintf(msg, "quads with 7 cones on a sphere: %s", qaws_status_to_string(st));
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_INVALID_ARGUMENT && !q, msg);
	/* the square: a +1 and -1 pair, the border's corners take the rest */
	{
		double da[3] = { -0.5, 0.1, 0 }, db[3] = { 0.5, -0.1, 0 };
		unsigned int ia, ib;
		pv[0] = gd_nearest(xs, nvs, da);
		pv[1] = gd_nearest(xs, nvs, db);
		for (k = 0; k < 3; k++)
		{
			pair[k] = xs[3 * pv[0] + k];
			pair[3 + k] = xs[3 * pv[1] + k];
		}
		gd_get(sq, &a);
		idx = (int*)malloc(sizeof(int) * a.nf);
		memset(&dd, 0, sizeof dd);
		dd.method = QAWS_DIRECTION_FIELD_SMOOTHEST;
		dd.cone_vertices = pv;
		dd.cone_indices = pi;
		dd.cone_count = 2;
		st = qaws_mesh_3d_direction_field(sq, &dd, NULL, idx);
		n = gd_singular(&a, idx, pv[0], 1, &ia);
		gd_singular(&a, idx, pv[1], -1, &ib);
		sprintf(msg, "square, a +1 and -1 pair placed (SMOOTHEST): %s, %u singular faces, +1 at its vertex %u, -1 at its %u", qaws_status_to_string(st), n,
			ia, ib);
		printf("  %s\n", msg);
		TEST_ASSERT(st == QAWS_STATUS_OK && n == 2 && ia == 1 && ib == 1, msg);
		free(idx);
		gd_free(&a);
		for (j = 0; j < 2; j++)
		{
			memset(&qd, 0, sizeof qd);
			qd.method = j ? QAWS_QUADS_GLOBAL : QAWS_QUADS_INSTANT;
			qd.target_length = (qaws_scalar)0.12;
			qd.pure = 1;
			qd.cone_vertices = pv;
			qd.cone_indices = pi;
			qd.cone_count = 2;
			memset(&rep, 0, sizeof rep);
			q = gd_quadrangulate(sq, &qd, &rep, &st);
			memset(&a, 0, sizeof a);
			n = 0;
			worst = HUGE_VAL;
			val[0] = val[1] = 0;
			if (q)
			{
				gd_get(q, &a);
				n = gd_irregular(&a, pair, 2, &worst, val, NULL, NULL);
			}
			sprintf(msg, "square quads %s with the pair: %s, %u faces, %u irregular inside; cones: nearest vertex %.2g, valences %u and %u; map %u",
				j ? "GLOBAL" : "INSTANT", qaws_status_to_string(st), a.nf, n, worst, val[0], val[1], (unsigned int)rep.iterations);
			printf("  %s\n", msg);
			TEST_ASSERT(st == QAWS_STATUS_OK && n == 2 && worst < 1e-12 && val[0] == 3 && val[1] == 5, msg);
			gd_free(&a);
			qaws_mesh_3d_destroy(q);
		}
	}
	qaws_mesh_3d_destroy(s);
	qaws_mesh_3d_destroy(sq);
	free(x0);
	free(xs);
}

/* attractors: a free cone of the sphere's smoothest field drawn to a point, and pushed off one */
static void test_attractors(void)
{
	double* x0;
	unsigned int nv0, j;
	qaws_mesh_3d* s = gd_sphere(24, 48, &x0, &nv0);
	qaws_direction_field_desc dd;
	gd_arrays a;
	int* idx;
	double dist[6], p[3], c0[3] = { 0, 0, 1 };
	unsigned int ns[6], tmp;
	qaws_scalar ap[3], str;
	char msg[500];
	gd_get(s, &a);
	idx = (int*)malloc(sizeof(int) * a.nf);
	/* p off the free cones (the direction (1, 2, 3)); the repeller on the free cone nearest it */
	p[0] = 1 / sqrt(14.0);
	p[1] = 2 / sqrt(14.0);
	p[2] = 3 / sqrt(14.0);
	for (j = 0; j < 6; j++)
	{
		double const* at = j % 3 == 2 ? c0 : p;
		memset(&dd, 0, sizeof dd);
		dd.method = j < 3 ? QAWS_DIRECTION_FIELD_SMOOTHEST : QAWS_DIRECTION_FIELD_LOCAL;
		str = (qaws_scalar)((j % 3 == 0 ? 0 : j % 3 == 1 ? 1 : -1) * (j < 3 ? 1 : 3));   /* (SMOOTHEST in cones, LOCAL in weights) */
		ap[0] = (qaws_scalar)at[0];
		ap[1] = (qaws_scalar)at[1];
		ap[2] = (qaws_scalar)at[2];
		dd.attractors = ap;
		dd.attractor_strengths = &str;
		dd.attractor_count = str != 0;
		dd.attractor_radius = (qaws_scalar)0.3;
		qaws_mesh_3d_direction_field(s, &dd, NULL, idx);
		dist[j] = gd_cone_distance(&a, idx, at, j % 3 == 0 ? c0 : NULL);
		ns[j] = gd_singular(&a, idx, 0, 0, &tmp);
	}
	sprintf(msg, "attractors on the sphere (radius 0.3), the nearest cone: SMOOTHEST to p %.3f free, %.3f drawn (strength 1), to a free cone %.3f pushed off it (-1); LOCAL %.3f, %.3f (3), %.3f (-3); singular faces %u %u %u, %u %u %u",
		dist[0], dist[1], dist[2], dist[3], dist[4], dist[5], ns[0], ns[1], ns[2], ns[3], ns[4], ns[5]);
	printf("  %s\n", msg);
	TEST_ASSERT(dist[1] < 0.1 && dist[1] < 0.5 * dist[0] && dist[2] > 0.2 && dist[4] < dist[3] && dist[5] > 0.2 && ns[0] == 8 && ns[1] == 8 && ns[2] == 8, msg);
	free(idx);
	gd_free(&a);
	{
		/* quads at 0.15: an irregular vertex drawn to a point (strength 1, reach 0.1, about the lattice's spacing) */
		double near[2] = { HUGE_VAL, HUGE_VAL }, worst, pq[3] = { 0.62, -0.45, 0.64 }, l = sqrt(0.62 * 0.62 + 0.45 * 0.45 + 0.64 * 0.64);
		unsigned int irr[2] = { 0, 0 }, val[1], k;
		qaws_status st[2];
		for (k = 0; k < 3; k++)
		{
			pq[k] /= l;
			ap[k] = (qaws_scalar)pq[k];
		}
		str = 1;
		for (j = 0; j < 2; j++)
		{
			qaws_quadrangulate_desc qd;
			qaws_mesh_3d* q;
			memset(&qd, 0, sizeof qd);
			qd.target_length = (qaws_scalar)0.15;
			qd.pure = 1;
			qd.attractors = ap;
			qd.attractor_strengths = &str;
			qd.attractor_count = j;
			qd.attractor_radius = (qaws_scalar)0.1;
			q = gd_quadrangulate(s, &qd, NULL, &st[j]);
			if (!q)
				continue;
			gd_get(q, &a);
			irr[j] = gd_irregular(&a, pq, 1, &worst, val, pq, &near[j]);
			gd_free(&a);
			qaws_mesh_3d_destroy(q);
		}
		sprintf(msg, "quads with an attractor (strength 1, reach 0.1): %s, the nearest irregular vertex %.3f (%u irregular), %.3f without (%u)",
			qaws_status_to_string(st[1]), near[1], irr[1], near[0], irr[0]);
		printf("  %s\n", msg);
		TEST_ASSERT(st[0] == QAWS_STATUS_OK && st[1] == QAWS_STATUS_OK && near[1] < 0.1 && near[0] > 0.3 && irr[1] == 8, msg);
	}
	qaws_mesh_3d_destroy(s);
	free(x0);
}

static void test_refusals(void)
{
	double *x0, *x1;
	unsigned int nv0, nv1, bad[2] = { 0, 100000 };
	qaws_mesh_3d *s = gd_sphere(6, 8, &x0, &nv0), *t = gd_sphere(7, 9, &x1, &nv1), *out = NULL;
	qaws_mesh_reference* rt = NULL;
	qaws_scalar* g = (qaws_scalar*)calloc(3 * nv0, sizeof(qaws_scalar));
	qaws_quadrangulate_desc d;
	qaws_remesh_desc rd;
	qaws_status s1, s2, s3, s4;
	char msg[300];
	memset(&d, 0, sizeof d);
	d.pinned = bad + 1;
	d.pinned_count = 1;
	s1 = qaws_mesh_quadrangulate(s, &d, &out, NULL);
	memset(&d, 0, sizeof d);
	d.guide_edges = bad;
	d.guide_edge_count = 1;
	s2 = qaws_mesh_quadrangulate(s, &d, &out, NULL);
	qaws_mesh_reference_create(t, NULL, &rt);
	memset(&d, 0, sizeof d);
	d.reference = rt;
	d.guide_directions = g;
	s3 = qaws_mesh_quadrangulate(s, &d, &out, NULL);
	memset(&rd, 0, sizeof rd);
	rd.pinned = bad + 1;
	rd.pinned_count = 1;
	s4 = qaws_mesh_remesh(s, &rd, &out, NULL);
	sprintf(msg, "refusals: a pin past the vertices %s, a guide edge past them %s, directions on a foreign reference %s, remeshing's pin %s",
		qaws_status_to_string(s1), qaws_status_to_string(s2), qaws_status_to_string(s3), qaws_status_to_string(s4));
	TEST_ASSERT(s1 == QAWS_STATUS_INVALID_ARGUMENT && s2 == QAWS_STATUS_INVALID_ARGUMENT && s3 == QAWS_STATUS_INVALID_ARGUMENT
		&& s4 == QAWS_STATUS_INVALID_ARGUMENT && !out, msg);
	qaws_mesh_reference_destroy(rt);
	qaws_mesh_3d_destroy(s);
	qaws_mesh_3d_destroy(t);
	free(g); free(x0); free(x1);
}

int test_137_mesh_guidance_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 137: user guidance\n");
	test_refusals();
	test_pinned();
	test_cones();
	test_attractors();
	test_directions();
	test_edge_loops();
	test_density();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
