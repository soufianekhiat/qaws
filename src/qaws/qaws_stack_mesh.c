/*
 * The surface of a stack of contours as a 3D triangle mesh: sections
 * sampled on their curves, contours paired between heights and stitched,
 * fans at tips, caps where the topology changes and at the ends.
 */

#include "qaws_stack.h"
#include "qaws_triangulate.h"
#include "qaws_curve.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "qaws_platform.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_parallel.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define SM_NONE 0xFFFFFFFFu
#define SM_DENSE 32             /* samples per span to measure arc length */

struct qaws_mesh_3d
{
	double* v;                  /* 3 per vertex */
	unsigned int nv, capv;
	unsigned int* t;            /* 3 per triangle */
	unsigned int nt, capt;
	int failed;
};

static unsigned int sm_vert(qaws_mesh_3d* m, double x, double y, double z)
{
	if (m->nv == m->capv)
	{
		unsigned int nc = m->capv ? 2 * m->capv : 1024;
		double* g = (double*)realloc(m->v, sizeof(double) * 3 * nc);
		if (!g) { m->failed = 1; return SM_NONE; }
		m->v = g;
		m->capv = nc;
	}
	m->v[3 * m->nv] = x; m->v[3 * m->nv + 1] = y; m->v[3 * m->nv + 2] = z;
	return m->nv++;
}

static void sm_tri(qaws_mesh_3d* m, unsigned int a, unsigned int b, unsigned int c)
{
	if (a == b || b == c || a == c) return;
	if (m->nt == m->capt)
	{
		unsigned int nc = m->capt ? 2 * m->capt : 2048;
		unsigned int* g = (unsigned int*)realloc(m->t, sizeof(unsigned int) * 3 * nc);
		if (!g) { m->failed = 1; return; }
		m->t = g;
		m->capt = nc;
	}
	m->t[3 * m->nt] = a; m->t[3 * m->nt + 1] = b; m->t[3 * m->nt + 2] = c;
	m->nt++;
}

/* ======================================================================== */
/*  Sections and their contours                                             */
/* ======================================================================== */

typedef struct sm_contour
{
	qaws_path_2d path;
	double length, area;        /* signed area: counter-clockwise outer, clockwise hole */
	double c[2], lo[2], hi[2];
	unsigned int n;             /* samples (the chain's) */
	unsigned int first;         /* its first vertex in the mesh */
	unsigned int up, down;      /* paired contour index in the next / previous section */
	unsigned int chain;
} sm_contour;

typedef struct sm_section
{
	double z;
	qaws_clip_result* r;
	sm_contour* c;
	unsigned int nc;
} sm_section;

typedef struct sm_job
{
	qaws_stack_2d const* st;
	sm_section* s;
} sm_job;

static qaws_status sm_section_job(void* ctx, unsigned int chunk, unsigned int begin, unsigned int end)
{
	sm_job* j = (sm_job*)ctx;
	unsigned int i;
	(void)chunk;
	for (i = begin; i < end; i++)
	{
		qaws_status s = qaws_stack_section_2d(j->st, (qaws_scalar)j->s[i].z, &j->s[i].r);
		if (s != QAWS_STATUS_OK) return s;
	}
	return QAWS_STATUS_OK;
}

static void sm_eval(qaws_curve const* c, double t, double* p)
{
	qaws_eval_result_2d r;
	qaws_curve_evaluate_2d(c, (qaws_scalar)t, QAWS_EVAL_FLAG_POSITION, &r);
	p[0] = r.position.x;
	p[1] = r.position.y;
}

/* a contour's dense polyline: (curve, parameter, cumulative length) per point */
typedef struct sm_dense
{
	unsigned int* curve;
	double* t;
	double* s;
	double* p;
	unsigned int n;
} sm_dense;

static int sm_densify(qaws_path_2d const* p, sm_dense* d)
{
	unsigned int k, total = 1, n = 0;
	for (k = 0; k < p->curve_count; k++) total += SM_DENSE * p->curves[k]->span_count;
	d->curve = (unsigned int*)malloc(sizeof(unsigned int) * total);
	d->t = (double*)malloc(sizeof(double) * total);
	d->s = (double*)malloc(sizeof(double) * total);
	d->p = (double*)malloc(sizeof(double) * 2 * total);
	if (!d->curve || !d->t || !d->s || !d->p) return 0;
	for (k = 0; k < p->curve_count; k++)
	{
		qaws_curve const* c = p->curves[k];
		unsigned int m = SM_DENSE * c->span_count, j;
		double t0 = c->parameter_range.min_value, t1 = c->parameter_range.max_value;
		for (j = (k ? 1u : 0u); j <= m; j++)
		{
			double t = t0 + (t1 - t0) * j / m;
			d->curve[n] = k;
			d->t[n] = t;
			sm_eval(c, t, d->p + 2 * n);
			d->s[n] = n ? d->s[n - 1] + hypot(d->p[2 * n] - d->p[2 * n - 2], d->p[2 * n + 1] - d->p[2 * n - 1]) : 0;
			n++;
		}
	}
	d->n = n;
	return 1;
}

static void sm_dense_free(sm_dense* d)
{
	free(d->curve); free(d->t); free(d->s); free(d->p);
	memset(d, 0, sizeof(*d));
}

static int sm_measure(sm_contour* c)
{
	sm_dense d;
	unsigned int i;
	double a = 0, cx = 0, cy = 0;
	memset(&d, 0, sizeof(d));
	if (!sm_densify(&c->path, &d)) { sm_dense_free(&d); return 0; }
	c->length = d.n ? d.s[d.n - 1] : 0;
	c->lo[0] = c->lo[1] = 1e300;
	c->hi[0] = c->hi[1] = -1e300;
	for (i = 0; i < d.n; i++)
	{
		double const* p = d.p + 2 * i;
		double const* q = d.p + 2 * ((i + 1) % d.n);
		double w = p[0] * q[1] - q[0] * p[1];
		a += w;
		cx += (p[0] + q[0]) * w;
		cy += (p[1] + q[1]) * w;
		if (p[0] < c->lo[0]) c->lo[0] = p[0];
		if (p[1] < c->lo[1]) c->lo[1] = p[1];
		if (p[0] > c->hi[0]) c->hi[0] = p[0];
		if (p[1] > c->hi[1]) c->hi[1] = p[1];
	}
	c->area = 0.5 * a;
	if (fabs(a) > 0) { c->c[0] = cx / (3 * a); c->c[1] = cy / (3 * a); }
	else { c->c[0] = 0.5 * (c->lo[0] + c->hi[0]); c->c[1] = 0.5 * (c->lo[1] + c->hi[1]); }
	sm_dense_free(&d);
	return 1;
}

/* the contour's n samples by arc length, on its curves, into the mesh */
static int sm_emit(qaws_mesh_3d* m, sm_contour* c, double z)
{
	sm_dense d;
	unsigned int i, j = 0;
	memset(&d, 0, sizeof(d));
	if (!sm_densify(&c->path, &d)) { sm_dense_free(&d); return 0; }
	c->first = m->nv;
	for (i = 0; i < c->n; i++)
	{
		double target = c->length * i / c->n, p[2], t, f;
		while (j + 2 < d.n && d.s[j + 1] < target) j++;
		f = d.s[j + 1] > d.s[j] ? (target - d.s[j]) / (d.s[j + 1] - d.s[j]) : 0;
		if (d.curve[j] == d.curve[j + 1])
		{
			t = d.t[j] + f * (d.t[j + 1] - d.t[j]);
			sm_eval(c->path.curves[d.curve[j]], t, p);
		}
		else
		{
			p[0] = d.p[2 * j] + f * (d.p[2 * j + 2] - d.p[2 * j]);
			p[1] = d.p[2 * j + 1] + f * (d.p[2 * j + 3] - d.p[2 * j + 1]);
		}
		if (sm_vert(m, p[0], p[1], z) == SM_NONE) { sm_dense_free(&d); return 0; }
	}
	sm_dense_free(&d);
	return 1;
}

/* ======================================================================== */
/*  Pairing, bands, fans, caps                                              */
/* ======================================================================== */

/* contours of sections a and b paired: same orientation, overlapping
   boxes, mutually nearest centres */
static void sm_pair(sm_section* a, sm_section* b)
{
	unsigned int i, j;
	for (i = 0; i < a->nc; i++)
	{
		sm_contour* p = &a->c[i];
		unsigned int best = SM_NONE;
		double bd = 1e300;
		for (j = 0; j < b->nc; j++)
		{
			sm_contour const* q = &b->c[j];
			double d;
			if ((p->area > 0) != (q->area > 0)) continue;
			if (p->hi[0] < q->lo[0] || q->hi[0] < p->lo[0] || p->hi[1] < q->lo[1] || q->hi[1] < p->lo[1]) continue;
			d = hypot(p->c[0] - q->c[0], p->c[1] - q->c[1]);
			if (d < bd) { bd = d; best = j; }
		}
		if (best == SM_NONE) continue;
		/* mutual: no contour of a nearer to it */
		{
			sm_contour const* q = &b->c[best];
			unsigned int k;
			for (k = 0; k < a->nc; k++)
			{
				sm_contour const* r = &a->c[k];
				if (k == i || (r->area > 0) != (q->area > 0)) continue;
				if (r->hi[0] < q->lo[0] || q->hi[0] < r->lo[0] || r->hi[1] < q->lo[1] || q->hi[1] < r->lo[1]) continue;
				if (hypot(r->c[0] - q->c[0], r->c[1] - q->c[1]) < bd) break;
			}
			if (k < a->nc) continue;
		}
		p->up = best;
		b->c[best].down = i;
	}
}

static double sm_d2(qaws_mesh_3d const* m, unsigned int a, unsigned int b)
{
	double dx = m->v[3 * a] - m->v[3 * b], dy = m->v[3 * a + 1] - m->v[3 * b + 1], dz = m->v[3 * a + 2] - m->v[3 * b + 2];
	return dx * dx + dy * dy + dz * dz;
}

/* the band between paired contours p (below) and q (above), n samples each */
static void sm_band(qaws_mesh_3d* m, sm_contour const* p, sm_contour const* q)
{
	unsigned int n = p->n, s, best = 0, i;
	double bd = 1e300;
	/* the alignment: q rotated to the least squared distance to p */
	for (s = 0; s < n; s++)
	{
		double acc = 0;
		unsigned int step = n > 256 ? n / 256 : 1;
		for (i = 0; i < n; i += step)
		{
			double dx = m->v[3 * (p->first + i)] - m->v[3 * (q->first + (i + s) % n)];
			double dy = m->v[3 * (p->first + i) + 1] - m->v[3 * (q->first + (i + s) % n) + 1];
			acc += dx * dx + dy * dy;
			if (acc >= bd) break;
		}
		if (acc < bd) { bd = acc; best = s; }
	}
	for (i = 0; i < n; i++)
	{
		unsigned int a = p->first + i, b = p->first + (i + 1) % n;
		unsigned int c = q->first + (i + 1 + best) % n, d = q->first + (i + best) % n;
		/* travel a -> b with the outside on the right, up toward d, c */
		if (sm_d2(m, a, c) <= sm_d2(m, b, d))
		{
			sm_tri(m, a, b, c);
			sm_tri(m, a, c, d);
		}
		else
		{
			sm_tri(m, a, b, d);
			sm_tri(m, b, c, d);
		}
	}
}

/* a fan from contour p to its centre at height z (up: the tip above p;
   the triangles keep the contour's travel, so holes need no turn) */
static void sm_fan(qaws_mesh_3d* m, sm_contour const* p, double z, int up)
{
	unsigned int apex = sm_vert(m, p->c[0], p->c[1], z), i;
	if (apex == SM_NONE) return;
	for (i = 0; i < p->n; i++)
	{
		unsigned int a = p->first + i, b = p->first + (i + 1) % p->n;
		if (up) sm_tri(m, a, b, apex);
		else sm_tri(m, b, a, apex);
	}
}

/* flat faces over the contours cs (their samples' rings) at their height:
   facing up when up, down otherwise; the triangulation uses the samples as
   its only vertices */
static int sm_cap(qaws_mesh_3d* m, sm_contour const* const* cs, unsigned int n, int up)
{
	qaws_curve** rings = (qaws_curve**)calloc(n + 1, sizeof(qaws_curve*));
	qaws_curve const** views = (qaws_curve const**)malloc(sizeof(qaws_curve*) * (n + 1));
	qaws_path_2d* paths = (qaws_path_2d*)malloc(sizeof(qaws_path_2d) * (n + 1));
	qaws_triangulate_desc d;
	qaws_mesh_2d* mesh = NULL;
	unsigned int i, k, total = 0;
	int ok = 1;
	if (!rings || !views || !paths) ok = 0;
	for (i = 0; ok && i < n; i++)
	{
		qaws_scalar* xy = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * cs[i]->n);
		if (!xy) { ok = 0; break; }
		for (k = 0; k < cs[i]->n; k++)
		{
			xy[2 * k] = (qaws_scalar)m->v[3 * (cs[i]->first + k)];
			xy[2 * k + 1] = (qaws_scalar)m->v[3 * (cs[i]->first + k) + 1];
		}
		if (qaws_curve_create_polyline_2d(xy, cs[i]->n, 1, &rings[i]) != QAWS_STATUS_OK) ok = 0;
		free(xy);
		views[i] = rings[i];
		paths[i].curves = &views[i];
		paths[i].curve_count = 1;
		paths[i].closed = 1;
		total += cs[i]->n;
	}
	if (ok)
	{
		memset(&d, 0, sizeof(d));
		d.paths = paths;
		d.path_count = n;
		d.fill_rule = QAWS_FILL_EVEN_ODD;
		d.tolerance = (qaws_scalar)1e300;
		if (qaws_triangulate_region_2d(&d, &mesh) != QAWS_STATUS_OK) ok = 0;
	}
	if (ok && mesh)
	{
		unsigned int nv = qaws_mesh_2d_get_vertex_count(mesh), nt = qaws_mesh_2d_get_triangle_count(mesh);
		unsigned int* map = (unsigned int*)malloc(sizeof(unsigned int) * (nv + 1));
		double z = m->v[3 * cs[0]->first + 2];
		if (!map) ok = 0;
		/* each vertex back to its sample (the same point up to the lattice) */
		for (i = 0; ok && i < nv; i++)
		{
			qaws_vec2 q = qaws_mesh_2d_get_vertex(mesh, i);
			unsigned int best = SM_NONE, c, s;
			double bd = 1e300;
			for (c = 0; c < n; c++)
				for (s = 0; s < cs[c]->n; s++)
				{
					unsigned int v = cs[c]->first + s;
					double dx = m->v[3 * v] - (double)q.x, dy = m->v[3 * v + 1] - (double)q.y, dd = dx * dx + dy * dy;
					if (dd < bd) { bd = dd; best = v; }
				}
			/* a point the triangulation made (where rings cross): its own vertex */
			if (best == SM_NONE || bd > 1e-18 * (1 + fabs((double)q.x) + fabs((double)q.y)))
				best = sm_vert(m, q.x, q.y, z);
			map[i] = best;
		}
		for (i = 0; ok && i < nt; i++)
		{
			unsigned int t[3];
			qaws_mesh_2d_get_triangle(mesh, i, t);
			if (up) sm_tri(m, map[t[0]], map[t[1]], map[t[2]]);
			else sm_tri(m, map[t[0]], map[t[2]], map[t[1]]);
		}
		free(map);
	}
	(void)total;
	qaws_mesh_2d_destroy(mesh);
	for (i = 0; rings && i < n; i++) qaws_curve_destroy(rings[i]);
	free(rings); free((void*)views); free(paths);
	return ok && !m->failed;
}

/* ======================================================================== */
/*  Public                                                                  */
/* ======================================================================== */

static unsigned int sm_find(unsigned int* parent, unsigned int x)
{
	while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
	return x;
}

qaws_status qaws_stack_mesh_3d(qaws_stack_mesh_desc const* desc, qaws_mesh_3d** out)
{
	qaws_stack_2d const* st;
	qaws_mesh_3d* m;
	sm_section* s = NULL;
	unsigned int ns = 0, i, j, total = 0, per;
	unsigned int *base = NULL, *parent = NULL, *want = NULL;
	double spacing, ext = 0;
	sm_job job;
	qaws_status status = QAWS_STATUS_OK;
	if (!desc || !out || !desc->stack) return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	st = desc->stack;
	m = (qaws_mesh_3d*)calloc(1, sizeof(qaws_mesh_3d));
	if (!m) return QAWS_STATUS_ALLOCATION_FAILURE;
	if (st->level_count == 0) { *out = m; return QAWS_STATUS_OK; }
	/* the heights: every level, and the sections between */
	per = desc->sections + 1;
	s = (sm_section*)calloc((st->level_count - 1) * per + 1, sizeof(sm_section));
	if (!s) { free(m); return QAWS_STATUS_ALLOCATION_FAILURE; }
	for (i = 0; i < st->level_count; i++)
	{
		s[ns++].z = (double)st->levels[i].z;
		if (i + 1 < st->level_count)
			for (j = 1; j < per; j++)
				s[ns++].z = (double)st->levels[i].z + ((double)st->levels[i + 1].z - (double)st->levels[i].z) * j / per;
	}
	job.st = st;
	job.s = s;
	status = qaws_internal_parallel(desc->executor, ns, 1, sm_section_job, &job);
	/* contours, measured */
	for (i = 0; status == QAWS_STATUS_OK && i < ns; i++)
	{
		unsigned int nc = s[i].r ? qaws_clip_result_get_path_count(s[i].r) : 0;
		s[i].c = (sm_contour*)calloc(nc + 1, sizeof(sm_contour));
		if (!s[i].c) { status = QAWS_STATUS_ALLOCATION_FAILURE; break; }
		s[i].nc = nc;
		for (j = 0; j < nc; j++)
		{
			sm_contour* c = &s[i].c[j];
			qaws_clip_result_get_path(s[i].r, j, &c->path);
			c->up = c->down = SM_NONE;
			if (!sm_measure(c)) { status = QAWS_STATUS_ALLOCATION_FAILURE; break; }
			if (c->hi[0] - c->lo[0] > ext) ext = c->hi[0] - c->lo[0];
			if (c->hi[1] - c->lo[1] > ext) ext = c->hi[1] - c->lo[1];
		}
	}
	spacing = desc->spacing > 0 ? (double)desc->spacing : (ext > 0 ? ext / 64 : 1);
	for (i = 0; status == QAWS_STATUS_OK && i + 1 < ns; i++)
		sm_pair(&s[i], &s[i + 1]);
	/* chains of paired contours share the largest sample count among them */
	if (status == QAWS_STATUS_OK)
	{
		base = (unsigned int*)malloc(sizeof(unsigned int) * (ns + 1));
		if (!base) status = QAWS_STATUS_ALLOCATION_FAILURE;
		for (i = 0; base && i < ns; i++) { base[i] = total; total += s[i].nc; }
		parent = (unsigned int*)malloc(sizeof(unsigned int) * (total + 1));
		want = (unsigned int*)calloc(total + 1, sizeof(unsigned int));
		if (!parent || !want) status = QAWS_STATUS_ALLOCATION_FAILURE;
	}
	if (status == QAWS_STATUS_OK)
	{
		for (i = 0; i < total; i++) parent[i] = i;
		for (i = 0; i + 1 < ns; i++)
			for (j = 0; j < s[i].nc; j++)
				if (s[i].c[j].up != SM_NONE)
				{
					unsigned int a = sm_find(parent, base[i] + j), b = sm_find(parent, base[i + 1] + s[i].c[j].up);
					if (a != b) parent[a] = b;
				}
		for (i = 0; i < ns; i++)
			for (j = 0; j < s[i].nc; j++)
			{
				unsigned int r = sm_find(parent, base[i] + j);
				unsigned int w = (unsigned int)ceil(s[i].c[j].length / spacing);
				if (w < 8) w = 8;
				if (w > want[r]) want[r] = w;
			}
		for (i = 0; i < ns && status == QAWS_STATUS_OK; i++)
			for (j = 0; j < s[i].nc; j++)
			{
				s[i].c[j].n = want[sm_find(parent, base[i] + j)];
				if (!sm_emit(m, &s[i].c[j], s[i].z)) { status = QAWS_STATUS_ALLOCATION_FAILURE; break; }
			}
	}
	/* bands, fans and caps */
	for (i = 0; status == QAWS_STATUS_OK && i < ns; i++)
	{
		for (j = 0; j < s[i].nc; j++)
		{
			sm_contour const* c = &s[i].c[j];
			int hole = c->area < 0;
			if (i + 1 < ns)
			{
				if (c->up != SM_NONE)
					sm_band(m, c, &s[i + 1].c[c->up]);
				else if (s[i + 1].nc == 0)
					sm_fan(m, c, s[i + 1].z, 1);
				else
				{
					sm_contour const* one = c;
					if (!sm_cap(m, &one, 1, !hole)) status = QAWS_STATUS_ALLOCATION_FAILURE;
				}
			}
			if (i > 0 && c->down == SM_NONE)
			{
				if (s[i - 1].nc == 0)
					sm_fan(m, c, s[i - 1].z, 0);
				else
				{
					sm_contour const* one = c;
					if (!sm_cap(m, &one, 1, hole)) status = QAWS_STATUS_ALLOCATION_FAILURE;
				}
			}
		}
	}
	if (status == QAWS_STATUS_OK && desc->caps)
	{
		unsigned int e;
		for (e = 0; e < 2 && status == QAWS_STATUS_OK; e++)
		{
			sm_section const* sec = &s[e ? ns - 1 : 0];
			sm_contour const** cs;
			if (!sec->nc || (e == 1 && ns == 1)) continue;
			cs = (sm_contour const**)malloc(sizeof(sm_contour*) * sec->nc);
			if (!cs) { status = QAWS_STATUS_ALLOCATION_FAILURE; break; }
			for (j = 0; j < sec->nc; j++) cs[j] = &sec->c[j];
			if (!sm_cap(m, cs, sec->nc, e == 1)) status = QAWS_STATUS_ALLOCATION_FAILURE;
			free((void*)cs);
		}
	}
	if (m->failed && status == QAWS_STATUS_OK) status = QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; s && i < ns; i++)
	{
		qaws_clip_result_destroy(s[i].r);
		free(s[i].c);
	}
	free(s); free(base); free(parent); free(want);
	if (status != QAWS_STATUS_OK) { qaws_mesh_3d_destroy(m); return status; }
	*out = m;
	return QAWS_STATUS_OK;
}

void qaws_mesh_3d_destroy(qaws_mesh_3d* m)
{
	if (!m) return;
	free(m->v);
	free(m->t);
	free(m);
}

unsigned int qaws_mesh_3d_get_vertex_count(qaws_mesh_3d const* m) { return m ? m->nv : 0; }

qaws_vec3 qaws_mesh_3d_get_vertex(qaws_mesh_3d const* m, unsigned int i)
{
	qaws_vec3 r;
	r.x = (qaws_scalar)m->v[3 * i];
	r.y = (qaws_scalar)m->v[3 * i + 1];
	r.z = (qaws_scalar)m->v[3 * i + 2];
	return r;
}

unsigned int qaws_mesh_3d_get_triangle_count(qaws_mesh_3d const* m) { return m ? m->nt : 0; }

void qaws_mesh_3d_get_triangle(qaws_mesh_3d const* m, unsigned int t, unsigned int out[3])
{
	out[0] = m->t[3 * t];
	out[1] = m->t[3 * t + 1];
	out[2] = m->t[3 * t + 2];
}
