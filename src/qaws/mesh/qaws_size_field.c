/*
 * Size fields: uniform, from a reference's curvature with tolerances and
 * gradation, by local feature size (medial balls on dense samples), per
 * reference vertex, by callback, the smaller of two.
 */

#include "qaws_mesh_internal.h"
#include "../internal/qaws_internal_memory.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

enum { SF_UNIFORM, SF_VERTICES, SF_SAMPLES, SF_CALLBACK, SF_MIN, SF_SCALED };

struct qaws_size_field
{
	int kind;
	double length;
	qaws_mesh_reference const* ref;
	double* vs;                     /* per vertex of the reference's mesh; SF_SAMPLES: per sample */
	unsigned int *lv, *first;       /* SF_SAMPLES: per triangle its grid's level and its first sample */
	qaws_scalar (*fn)(void* user, qaws_scalar const* point);
	void* user;
	qaws_size_field const *a, *b;
};

/* a triangle's barycentric grid of level n: sample (i, j) at (b1, b2) = (i, j) / n, i + j <= n, row j after rows 0..j-1 */
#define FS_AT(n, i, j) ((j) * ((n) + 1) - (j) * ((j) - 1) / 2 + (i))

/* the samples' sizes interpolated over the grid's cell holding the location */
static double fs_interp(unsigned int const* lv, unsigned int const* first, double const* vs, unsigned int t, double a, double b)
{
	unsigned int n = lv[t], i, j;
	double const* h = vs + first[t];
	double u = a * n, v = b * n, fu, fv;
	u = u < 0 ? 0 : u > n ? n : u;
	v = v < 0 ? 0 : v > n ? n : v;
	if (u + v > n)
	{
		double s = n / (u + v);
		u *= s;
		v *= s;
	}
	i = (unsigned int)u;
	if (i > n - 1)
		i = n - 1;
	j = (unsigned int)v;
	if (j > n - 1 - i)
		j = n - 1 - i;
	fu = u - i;
	fv = v - j;
	if (fu + fv <= 1)
		return (1 - fu - fv) * h[FS_AT(n, i, j)] + fu * h[FS_AT(n, i + 1, j)] + fv * h[FS_AT(n, i, j + 1)];
	if (i + j + 2 <= n)
		return (fu + fv - 1) * h[FS_AT(n, i + 1, j + 1)] + (1 - fu) * h[FS_AT(n, i, j + 1)] + (1 - fv) * h[FS_AT(n, i + 1, j)];
	return (fu * h[FS_AT(n, i + 1, j)] + fv * h[FS_AT(n, i, j + 1)]) / (fu + fv);   /* (rounding past the long side) */
}

double sf_eval(qaws_size_field const* f, qaws_mesh_location const* l, double const* x)
{
	switch (f->kind)
	{
	case SF_UNIFORM:
		return f->length;
	case SF_VERTICES:
	{
		unsigned int const* tv = f->ref->tv + 3 * (size_t)l->patch;
		double b1 = (double)l->a, b2 = (double)l->b;
		return (1 - b1 - b2) * f->vs[tv[0]] + b1 * f->vs[tv[1]] + b2 * f->vs[tv[2]];
	}
	case SF_SAMPLES:
		return fs_interp(f->lv, f->first, f->vs, l->patch, (double)l->a, (double)l->b);
	case SF_CALLBACK:
	{
		qaws_scalar p[3];
		p[0] = (qaws_scalar)x[0]; p[1] = (qaws_scalar)x[1]; p[2] = (qaws_scalar)x[2];
		return (double)f->fn(f->user, p);
	}
	case SF_SCALED:
		return f->length * sf_eval(f->a, l, x);
	default:
	{
		double u = sf_eval(f->a, l, x), v = sf_eval(f->b, l, x);
		return u < v ? u : v;
	}
	}
}

static qaws_size_field* sf_new(int kind)
{
	qaws_size_field* f = (qaws_size_field*)qaws_internal_zalloc(1, sizeof *f);
	if (f)
		f->kind = kind;
	return f;
}

qaws_status qaws_size_field_create_uniform(qaws_size_field** out, qaws_scalar length)
{
	if (!out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	if (!(length > 0))
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = sf_new(SF_UNIFORM);
	if (!*out)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	(*out)->length = (double)length;
	return QAWS_STATUS_OK;
}

qaws_status qaws_size_field_create_callback(qaws_size_field** out, qaws_scalar (*size)(void* user, qaws_scalar const* point), void* user)
{
	if (!out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	if (!size)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = sf_new(SF_CALLBACK);
	if (!*out)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	(*out)->fn = size;
	(*out)->user = user;
	return QAWS_STATUS_OK;
}

qaws_status qaws_size_field_create_min(qaws_size_field** out, qaws_size_field const* a, qaws_size_field const* b)
{
	qaws_mesh_reference const *ra, *rb;
	if (!out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	if (!a || !b)
		return QAWS_STATUS_INVALID_ARGUMENT;
	ra = qaws_size_field_get_reference(a);
	rb = qaws_size_field_get_reference(b);
	if (ra && rb && ra != rb)
		return QAWS_STATUS_INVALID_ARGUMENT;   /* (one reference's locations for both) */
	*out = sf_new(SF_MIN);
	if (!*out)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	(*out)->a = a;
	(*out)->b = b;
	return QAWS_STATUS_OK;
}

static qaws_status sf_vertices(qaws_mesh_reference const* r, double* vs, qaws_size_field** out)
{
	*out = sf_new(SF_VERTICES);
	if (!*out)
	{
		qaws_internal_free(vs);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	(*out)->ref = r;
	(*out)->vs = vs;
	return QAWS_STATUS_OK;
}

qaws_status qaws_size_field_create_vertices(qaws_size_field** out, qaws_mesh_reference const* r, qaws_scalar const* sizes)
{
	unsigned int nv, i;
	double* vs;
	if (!out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	if (!r || !sizes)
		return QAWS_STATUS_INVALID_ARGUMENT;
	nv = r->mesh->hm.nv;
	for (i = 0; i < nv; i++)
		if (!(sizes[i] > 0))
			return QAWS_STATUS_INVALID_ARGUMENT;
	vs = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), (nv + 1)));
	if (!vs)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < nv; i++)
		vs[i] = (double)sizes[i];
	return sf_vertices(r, vs, out);
}

static unsigned int sf_hm_ring(void* m, unsigned int u, unsigned int* out, unsigned int cap)
{
	return hm_ring((qaws_hm const*)m, u, out, cap);
}

/* sizes no faster growing than beta per unit of distance along a graph's
   edges (n points x, the neighbours of u from nbr: at most cap of them
   kept): Dijkstra from every point at once, h(v) <= h(u) + beta |uv| */
static int sf_gradation(unsigned int n, double const* x, double* h, double beta, unsigned int (*nbr)(void* user, unsigned int u, unsigned int* out,
	unsigned int cap), void* user)
{
	unsigned int *heap = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (n + 1))), *pos = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (n + 1)));
	unsigned int count = n, i, ring[256];
	if (!heap || !pos)
	{
		qaws_internal_free(heap);
		qaws_internal_free(pos);
		return 0;
	}
	for (i = 0; i < n; i++)
	{
		heap[i] = i;
		pos[i] = i;
	}
	/* heapify by h */
#define SF_LESS(x, y) (h[heap[x]] < h[heap[y]])
#define SF_SWAP(x, y) do { unsigned int t_ = heap[x]; heap[x] = heap[y]; heap[y] = t_; pos[heap[x]] = x; pos[heap[y]] = y; } while (0)
	for (i = n / 2; i-- > 0;)
	{
		unsigned int j = i;
		for (;;)
		{
			unsigned int l = 2 * j + 1, s = j;
			if (l < count && SF_LESS(l, s)) s = l;
			if (l + 1 < count && SF_LESS(l + 1, s)) s = l + 1;
			if (s == j) break;
			SF_SWAP(s, j);
			j = s;
		}
	}
	while (count)
	{
		unsigned int u = heap[0], k, nr;
		SF_SWAP(0, count - 1);
		pos[u] = HM_NONE;
		count--;
		{
			unsigned int j = 0;
			for (;;)
			{
				unsigned int l = 2 * j + 1, s = j;
				if (l < count && SF_LESS(l, s)) s = l;
				if (l + 1 < count && SF_LESS(l + 1, s)) s = l + 1;
				if (s == j) break;
				SF_SWAP(s, j);
				j = s;
			}
		}
		nr = nbr(user, u, ring, 256);
		for (k = 0; k < nr && k < 256; k++)
		{
			unsigned int w = ring[k], j;
			double d = 0, cand;
			for (j = 0; j < 3; j++)
				d += (x[3 * (size_t)u + j] - x[3 * (size_t)w + j]) * (x[3 * (size_t)u + j] - x[3 * (size_t)w + j]);
			cand = h[u] + beta * sqrt(d);
			if (pos[w] == HM_NONE || cand >= h[w])
				continue;
			h[w] = cand;
			for (j = pos[w]; j > 0 && SF_LESS(j, (j - 1) / 2); j = (j - 1) / 2)
				SF_SWAP(j, (j - 1) / 2);
		}
	}
#undef SF_LESS
#undef SF_SWAP
	qaws_internal_free(heap);
	qaws_internal_free(pos);
	return 1;
}

qaws_status qaws_size_field_create_curvature(qaws_size_field** out, qaws_mesh_reference const* r, qaws_curvature_size_desc const* desc)
{
	qaws_hm const* m;
	double *kmax, *vs, lo, hi, e, th;
	unsigned int t, i, k;
	if (!out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	if (!r || !desc || desc->chord_tolerance < 0 || desc->angle_tolerance < 0 || desc->min_length < 0 || desc->max_length < 0
		|| desc->gradation < 0 || (desc->chord_tolerance == 0 && desc->angle_tolerance == 0))
		return QAWS_STATUS_INVALID_ARGUMENT;
	m = &r->mesh->hm;
	lo = desc->min_length > 0 ? (double)desc->min_length : r->extent / 1000;
	hi = desc->max_length > 0 ? (double)desc->max_length : r->extent / 10;
	if (lo > hi)
		return QAWS_STATUS_INVALID_ARGUMENT;
	e = (double)desc->chord_tolerance;
	th = (double)desc->angle_tolerance;
	kmax = (double*)qaws_internal_scratch_zalloc((size_t)m->nv + 1, sizeof(double));
	vs = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), ((size_t)m->nv + 1)));
	if (!kmax || !vs)
	{
		qaws_internal_free(kmax);
		qaws_internal_free(vs);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	/* the largest curvature near each vertex: at its corners and at its triangles' centres */
	for (t = 0; t < r->nt; t++)
	{
		static double const at[4][2] = { { 0, 0 }, { 1, 0 }, { 0, 1 }, { 1.0 / 3, 1.0 / 3 } };
		double kc = 0;
		for (i = 0; i < 4; i++)
		{
			double xa[3], xb[3], n[3], na[3], nb[3], k1, k2, d1[3], d2[3], kk;
			mi_eval(&r->patch[t], at[i][0], at[i][1], NULL, xa, xb, NULL, NULL, NULL, n, na, nb);
			mi_curvature(xa, xb, n, na, nb, &k1, &k2, d1, d2);
			kk = fmax(fabs(k1), fabs(k2));
			if (i < 3)
				kmax[r->tv[3 * t + i]] = fmax(kmax[r->tv[3 * t + i]], kk);
			else
				kc = kk;
		}
		for (k = 0; k < 3; k++)
			kmax[r->tv[3 * t + k]] = fmax(kmax[r->tv[3 * t + k]], kc);
	}
	for (i = 0; i < m->nv; i++)
	{
		double L = hi;
		if (kmax[i] > 0)
		{
			double R = 1 / kmax[i];
			if (e > 0)
				L = fmin(L, e < R ? 2 * sqrt(e * (2 * R - e)) : 2 * R);
			if (th > 0)
				L = fmin(L, 2 * R * sin(0.5 * fmin(th, 3.14159265358979323846)));
		}
		vs[i] = L < lo ? lo : L > hi ? hi : L;
	}
	qaws_internal_free(kmax);
	if (desc->gradation > 0 && !sf_gradation(m->nv, m->x, vs, (double)desc->gradation, sf_hm_ring, (void*)m))
	{
		qaws_internal_free(vs);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	return sf_vertices(r, vs, out);
}

/* ---- sizing by local feature size --------------------------------------
 *
 * Samples on a barycentric grid per triangle; per sample and side the
 * largest ball tangent there that no opposite sheet enters (Ma, Bae, Choi,
 * Kim 2012: the ball shrunk through the nearest opposite sample inside it,
 * then its contact carried onto the reference by projection); the size is
 * the ball's diameter over cells_across. A first pass at half the largest
 * size, a second with a sample per size where the first's are smaller (its
 * balls only on the triangles the first found thin, the others at the
 * largest size until graded).
 */

#define FS_MAX_SAMPLES (1u << 18)

/* a k-d tree over points with normals, in place: each range's median on its widest axis at its middle, the range's box
   and its normals' box kept there */
typedef struct fs_tree
{
	double const *x, *nrm;
	unsigned int* perm;
	double* box;                    /* per range, at its middle: lo (3), hi (3), the normals' lo (3), hi (3) */
	unsigned int n;
} fs_tree;

static void fs_tree_build(fs_tree* t, unsigned int lo, unsigned int hi)
{
	unsigned int mid = lo + (hi - lo) / 2, ax = 0, i, k;
	double *bl = t->box + 12 * (size_t)mid, *bh = bl + 3;
	int a = (int)lo, b = (int)hi - 1;
	if (lo >= hi)
		return;
	for (k = 0; k < 3; k++)
	{
		bl[k] = bl[6 + k] = HUGE_VAL;
		bh[k] = bh[6 + k] = -HUGE_VAL;
	}
	for (i = lo; i < hi; i++)
		for (k = 0; k < 3; k++)
		{
			double v = t->x[3 * (size_t)t->perm[i] + k], u = t->nrm[3 * (size_t)t->perm[i] + k];
			bl[k] = fmin(bl[k], v);
			bh[k] = fmax(bh[k], v);
			bl[6 + k] = fmin(bl[6 + k], u);
			bh[6 + k] = fmax(bh[6 + k], u);
		}
	for (k = 1; k < 3; k++)
		if (bh[k] - bl[k] > bh[ax] - bl[ax])
			ax = k;
	/* quickselect (Hoare) of the median along ax */
	while (a < b)
	{
		double pv = t->x[3 * (size_t)t->perm[(a + b) / 2] + ax];
		int p = a, q = b;
		while (p <= q)
		{
			while (t->x[3 * (size_t)t->perm[p] + ax] < pv)
				p++;
			while (t->x[3 * (size_t)t->perm[q] + ax] > pv)
				q--;
			if (p <= q)
			{
				unsigned int s = t->perm[p];
				t->perm[p] = t->perm[q];
				t->perm[q] = s;
				p++;
				q--;
			}
		}
		if ((int)mid <= q)
			b = q;
		else if ((int)mid >= p)
			a = p;
		else
			break;
	}
	fs_tree_build(t, lo, mid);
	fs_tree_build(t, mid + 1, hi);
}

/* the squared distance from c to range [lo, hi)'s box (HUGE_VAL when empty) */
static double fs_box_d2(fs_tree const* t, unsigned int lo, unsigned int hi, double const* c)
{
	double const* b = t->box + 12 * (size_t)(lo + (hi - lo) / 2);
	double d2 = 0;
	unsigned int k;
	if (lo >= hi)
		return HUGE_VAL;
	for (k = 0; k < 3; k++)
	{
		double e = c[k] < b[k] ? b[k] - c[k] : c[k] > b[3 + k] ? c[k] - b[3 + k] : 0;
		d2 += e * e;
	}
	return d2;
}

/* the nearest accepted point to c closer than sqrt(best2) */
typedef struct fs_near
{
	fs_tree const* t;
	double c[3], best2;
	unsigned int at;                /* HM_NONE: none */
	int (*accept)(void* user, unsigned int i);
	int (*may)(void* user, double const* box);   /* 0: no point of a range with this box (12) is accepted */
	void* user;
} fs_near;

static void fs_nearest(fs_near* q, unsigned int lo, unsigned int hi)
{
	unsigned int mid = lo + (hi - lo) / 2, i, ax = 0, k;
	double const *p, *bx = q->t->box + 12 * (size_t)mid;
	double d2;
	if (fs_box_d2(q->t, lo, hi, q->c) >= q->best2 || !q->may(q->user, q->t->box + 12 * (size_t)mid))
		return;
	for (k = 1; k < 3; k++)   /* (the split's axis) */
		if (bx[3 + k] - bx[k] > bx[3 + ax] - bx[ax])
			ax = k;
	i = q->t->perm[mid];
	p = q->t->x + 3 * (size_t)i;
	d2 = (p[0] - q->c[0]) * (p[0] - q->c[0]) + (p[1] - q->c[1]) * (p[1] - q->c[1]) + (p[2] - q->c[2]) * (p[2] - q->c[2]);
	if (d2 < q->best2 && q->accept(q->user, i))
	{
		q->best2 = d2;
		q->at = i;
	}
	if (q->c[ax] <= p[ax])
	{
		fs_nearest(q, lo, mid);
		fs_nearest(q, mid + 1, hi);
	}
	else
	{
		fs_nearest(q, mid + 1, hi);
		fs_nearest(q, lo, mid);
	}
}

/* the points within sqrt(r2) of c but self, at most cap */
typedef struct fs_within
{
	fs_tree const* t;
	double c[3], r2;
	unsigned int *out, n, cap, self;
} fs_within;

static void fs_gather(fs_within* q, unsigned int lo, unsigned int hi)
{
	unsigned int mid = lo + (hi - lo) / 2, i;
	double const* p;
	if (q->n >= q->cap || fs_box_d2(q->t, lo, hi, q->c) > q->r2)
		return;
	i = q->t->perm[mid];
	p = q->t->x + 3 * (size_t)i;
	if (i != q->self && (p[0] - q->c[0]) * (p[0] - q->c[0]) + (p[1] - q->c[1]) * (p[1] - q->c[1]) + (p[2] - q->c[2]) * (p[2] - q->c[2]) <= q->r2)
		q->out[q->n++] = i;
	fs_gather(q, lo, mid);
	fs_gather(q, mid + 1, hi);
}

typedef struct fs_work
{
	qaws_mesh_reference const* r;
	unsigned int *lv, *first, *tri; /* per triangle its level and first sample; per sample its triangle */
	double *x, *nrm, *h;            /* per sample: position, unit normal, size */
	double *emax, *spacing;         /* per triangle: its longest side, its samples' spacing */
	unsigned int ns;
	fs_tree tree;
} fs_work;

/* neighbours for the gradation: the samples within 1.5 times the spacing of u's triangle */
static unsigned int fs_nbr(void* user, unsigned int u, unsigned int* out, unsigned int cap)
{
	fs_work const* w = (fs_work const*)user;
	fs_within q;
	double rad = 1.5 * w->spacing[w->tri[u]];
	q.t = &w->tree;
	memcpy(q.c, w->x + 3 * (size_t)u, sizeof q.c);
	q.r2 = rad * rad;
	q.out = out;
	q.n = 0;
	q.cap = cap;
	q.self = u;
	fs_gather(&q, 0, w->tree.n);
	return q.n;
}

/* the ball tangent at p on side d (unit) */
typedef struct fs_ball
{
	double const *x, *nrm;
	double p[3], np[3], d[3];
	double eps2;                    /* points closer than this (squared) are one: copies of a sample on two triangles */
} fs_ball;

/* the radius of the ball through q, or HUGE_VAL when q is no opposite sheet: behind p, or its normal within 60 degrees of
   perpendicular to p's (a face across a crease), or parallel to it (oriented alike: the same sheet, or a sheet of an open
   surface oriented the other way) and q more than 45 degrees off the ball's axis (contacts under 90 degrees apart) */
static double fs_radius(fs_ball const* b, double const* q, double const* nq)
{
	double dq[3], w, rr, c;
	unsigned int k;
	for (k = 0; k < 3; k++)
		dq[k] = q[k] - b->p[k];
	w = dq[0] * b->d[0] + dq[1] * b->d[1] + dq[2] * b->d[2];
	rr = dq[0] * dq[0] + dq[1] * dq[1] + dq[2] * dq[2];
	c = b->np[0] * nq[0] + b->np[1] * nq[1] + b->np[2] * nq[2];
	if (w > 0 && rr > b->eps2 && (c <= -0.5 || (c >= 0.5 && 2 * w * w >= rr)))
		return rr / (2 * w);
	return HUGE_VAL;
}

static int fs_accept(void* user, unsigned int i)
{
	fs_ball const* b = (fs_ball const*)user;
	return fs_radius(b, b->x + 3 * (size_t)i, b->nrm + 3 * (size_t)i) < HUGE_VAL;
}

/* may a range hold an opposite sheet: some of its box in front of the sample, its normals not all near perpendicular */
static int fs_may(void* user, double const* box)
{
	fs_ball const* b = (fs_ball const*)user;
	double front = 0, up = 0, down = 0;
	unsigned int k;
	for (k = 0; k < 3; k++)
	{
		double u = b->np[k] * box[6 + k], v = b->np[k] * box[9 + k];
		front += b->d[k] * ((b->d[k] > 0 ? box[3 + k] : box[k]) - b->p[k]);
		up += fmax(u, v);
		down += fmin(u, v);
	}
	return front > 0 && (up >= 0.5 || down <= -0.5);
}

/* the largest ball tangent at b's sample on its side that no opposite sheet enters, up to R; seed: how far to look for an
   opposite sample near the sample itself (beyond the samples' spacing a sheet has samples inside the balls) */
static double fs_side(fs_work const* w, fs_ball* b, double R, double seed)
{
	fs_near q;
	unsigned int hint = HM_NONE, at = HM_NONE, it, k;
	double r = R;
	q.t = &w->tree;
	q.accept = fs_accept;
	q.may = fs_may;
	q.user = b;
	/* the seed: the nearest opposite sample, for a sheet closer than the samples' spacing (none of its samples inside the
	   balls, too far off to the side) */
	memcpy(q.c, b->p, sizeof q.c);
	q.best2 = seed * seed;
	q.at = HM_NONE;
	fs_nearest(&q, 0, w->tree.n);
	if (q.at != HM_NONE)
	{
		hint = w->tri[at = q.at];
		r = fmin(R, fs_radius(b, w->x + 3 * (size_t)q.at, w->nrm + 3 * (size_t)q.at));
	}
	/* shrinking: the nearest opposite sample inside the ball, the ball through it; done within a hundredth (a sphere's
	   samples all lie about as far from its centre) */
	for (it = 0; it < 64; it++)
	{
		double rq;
		for (k = 0; k < 3; k++)
			q.c[k] = b->p[k] + r * b->d[k];
		q.best2 = r * r * (1 - 2e-2);
		q.at = HM_NONE;
		fs_nearest(&q, 0, w->tree.n);
		if (q.at == HM_NONE)
			break;
		rq = fs_radius(b, w->x + 3 * (size_t)q.at, w->nrm + 3 * (size_t)q.at);
		if (!(rq < r))
			break;
		r = rq;
		hint = w->tri[at = q.at];
	}
	/* the contact carried from the samples onto the opposite sheet: the closest point to the centre on the contact's triangle,
	   or across the side it lands on; done within a thousandth. Needless where the samples there are close beside the ball
	   (an offset s errs by about s^2 / 4r) or the contact sits on the ball's axis (within a hundredth of its depth squared) */
	if (at != HM_NONE)
	{
		double dq[3], dw, rr;
		for (k = 0; k < 3; k++)
			dq[k] = w->x[3 * (size_t)at + k] - b->p[k];
		dw = dq[0] * b->d[0] + dq[1] * b->d[1] + dq[2] * b->d[2];
		rr = dq[0] * dq[0] + dq[1] * dq[1] + dq[2] * dq[2];
		if (rr - dw * dw <= 1e-2 * dw * dw)
			hint = HM_NONE;
	}
	for (it = 0; hint != HM_NONE && w->spacing[hint] > 0.1 * r && it < 4; it++)
	{
		double c[3], y[3], ny[3], ry, a0, b0, best, ba, bb;
		unsigned int e, bt = hint;
		for (k = 0; k < 3; k++)
			c[k] = b->p[k] + r * b->d[k];
		best = mr_closest_on_patch(w->r, hint, c, &a0, &b0);
		ba = a0;
		bb = b0;
		for (e = 0; e < 3; e++)
		{
			unsigned int t = w->r->ta[3 * (size_t)hint + e];
			double a, bc, d;
			/* side e from corner e to e + 1: b = 0, a + b = 1, a = 0 */
			if (t == HM_NONE || (e == 0 ? b0 : e == 1 ? 1 - a0 - b0 : a0) > 1e-9)
				continue;
			d = mr_closest_on_patch(w->r, t / 3, c, &a, &bc);
			if (d < best)
			{
				best = d;
				bt = t / 3;
				ba = a;
				bb = bc;
			}
		}
		mi_eval(&w->r->patch[bt], ba, bb, y, NULL, NULL, NULL, NULL, NULL, ny, NULL, NULL);
		ry = fs_radius(b, y, ny);
		if (!(ry < r))
			break;
		k = ry > r * (1 - 1e-3);
		r = ry;
		hint = bt;
		if (k)
			break;
	}
	return r;
}

/* per triangle the level giving samples at most target[t] apart, the count held to FS_MAX_SAMPLES by widening every target
   alike; the samples' count */
static size_t fs_levels(fs_work* w, double const* target)
{
	unsigned int nt = w->r->nt, t;
	double scale = 1;
	for (;;)
	{
		size_t total = 0;
		int coarsest = 1;
		for (t = 0; t < nt; t++)
		{
			double q = w->emax[t] / (target[t] * scale);
			unsigned int n = q <= 1 ? 1 : q >= 1024 ? 1024 : (unsigned int)ceil(q);
			w->lv[t] = n;
			w->first[t] = (unsigned int)total;
			total += (size_t)(n + 1) * (n + 2) / 2;
			coarsest &= n == 1;
		}
		if (total <= FS_MAX_SAMPLES || coarsest)
			return total;
		scale *= 1.1 * sqrt((double)total / FS_MAX_SAMPLES);
	}
}

static void fs_free_samples(fs_work* w)
{
	qaws_internal_free(w->x);
	qaws_internal_free(w->nrm);
	qaws_internal_free(w->h);
	qaws_internal_free(w->tri);
	qaws_internal_free(w->tree.perm);
	qaws_internal_free(w->tree.box);
	w->x = w->nrm = w->h = NULL;
	w->tri = w->tree.perm = NULL;
	w->tree.box = NULL;
}

/* the samples of w's levels, their sizes graded; 0 on allocation failure. thin (per triangle): set where a size is under hi
   when first, else read: the balls only there, hi elsewhere */
static int fs_pass(fs_work* w, size_t ns, double cells, double lo, double hi, double beta, unsigned char* thin, int first)
{
	qaws_mesh_reference const* r = w->r;
	unsigned int t, i, j;
	double R = 0.5 * cells * hi, smax = 0;
	w->ns = (unsigned int)ns;
	w->x = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (ns + 1)));
	w->nrm = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (ns + 1)));
	w->h = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), (ns + 1)));
	w->tri = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (ns + 1)));
	w->tree.perm = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (ns + 1)));
	w->tree.box = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 12, (ns + 1)));
	if (!w->x || !w->nrm || !w->h || !w->tri || !w->tree.perm || !w->tree.box)
		return 0;
	for (t = 0; t < r->nt; t++)
	{
		unsigned int n = w->lv[t];
		w->spacing[t] = w->emax[t] / n;
		smax = fmax(smax, w->spacing[t]);
		for (j = 0; j <= n; j++)
			for (i = 0; i + j <= n; i++)
			{
				size_t k = w->first[t] + FS_AT(n, i, j);
				mi_eval(&r->patch[t], (double)i / n, (double)j / n, w->x + 3 * k, NULL, NULL, NULL, NULL, NULL, w->nrm + 3 * k, NULL, NULL);
				w->tri[k] = t;
			}
	}
	w->tree.x = w->x;
	w->tree.nrm = w->nrm;
	w->tree.n = w->ns;
	for (i = 0; i < w->ns; i++)
		w->tree.perm[i] = i;
	fs_tree_build(&w->tree, 0, w->ns);
	for (i = 0; i < w->ns; i++)
	{
		fs_ball b;
		double rad = R;
		int side;
		if (!first && !thin[w->tri[i]])
		{
			w->h[i] = hi;
			continue;
		}
		b.x = w->x;
		b.nrm = w->nrm;
		b.eps2 = 1e-14 * r->extent * r->extent;
		memcpy(b.p, w->x + 3 * (size_t)i, sizeof b.p);
		memcpy(b.np, w->nrm + 3 * (size_t)i, sizeof b.np);
		for (side = -1; side <= 1; side += 2)   /* inward (thickness) and outward (gaps), for outward normals */
		{
			unsigned int k;
			for (k = 0; k < 3; k++)
				b.d[k] = side * b.np[k];
			rad = fmin(rad, fs_side(w, &b, rad, fmin(2 * R, 1.5 * smax)));
		}
		w->h[i] = fmin(hi, fmax(lo, 2 * rad / cells));
		if (first && w->h[i] < hi)
			thin[w->tri[i]] = 1;
	}
	return sf_gradation(w->ns, w->x, w->h, beta, fs_nbr, w);
}

qaws_status qaws_size_field_create_feature_size(qaws_size_field** out, qaws_mesh_reference const* r, qaws_feature_size_desc const* desc)
{
	qaws_feature_size_desc dd;
	fs_work w;
	double lo, hi, cells, beta, *target = NULL;
	unsigned int* prev = NULL;
	unsigned char* thin = NULL;
	unsigned int nt, t, i;
	size_t ns;
	int ok = 0, finer = 0;
	if (!out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	if (!desc)
	{
		memset(&dd, 0, sizeof dd);
		desc = &dd;
	}
	if (!r || !(desc->cells_across >= 0) || !(desc->min_length >= 0) || !(desc->max_length >= 0) || !(desc->gradation >= 0))
		return QAWS_STATUS_INVALID_ARGUMENT;
	nt = r->nt;
	cells = desc->cells_across > 0 ? (double)desc->cells_across : 3;
	lo = desc->min_length > 0 ? (double)desc->min_length : r->extent / 1000;
	hi = desc->max_length > 0 ? (double)desc->max_length : r->extent / 10;
	if (nt && (!(lo > 0) || lo > hi))
		return QAWS_STATUS_INVALID_ARGUMENT;
	/* unlimited: only samples at one place (shared corners and sides) joined */
	beta = desc->gradation > 0 ? (double)desc->gradation : hi / (1e-6 * r->extent + 1e-300);
	memset(&w, 0, sizeof w);
	w.r = r;
	w.lv = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (nt + 1)));
	w.first = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (nt + 1)));
	w.emax = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), (nt + 1)));
	w.spacing = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), (nt + 1)));
	target = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), (nt + 1)));
	prev = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (nt + 1)));
	thin = (unsigned char*)qaws_internal_scratch_zalloc(nt + 1, 1);
	if (!w.lv || !w.first || !w.emax || !w.spacing || !target || !prev || !thin)
		goto done;
	for (t = 0; t < nt; t++)
	{
		mi_tri const* pt = &r->patch[t];
		double e = 0;
		for (i = 0; i < 3; i++)
		{
			double const *a = pt->p[i], *b = pt->p[(i + 1) % 3];
			e = fmax(e, sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2])));
		}
		w.emax[t] = e;
		target[t] = 0.5 * hi;
	}
	ns = fs_levels(&w, target);
	if (!fs_pass(&w, ns, cells, lo, hi, beta, thin, 1))
		goto done;
	/* a sample per size where the sizes are smaller than the first pass's spacing */
	for (t = 0; t < nt; t++)
	{
		unsigned int n = w.lv[t], k, m = (n + 1) * (n + 2) / 2;
		double hm = hi;
		for (k = 0; k < m; k++)
			hm = fmin(hm, w.h[w.first[t] + k]);
		target[t] = fmin(0.5 * hi, hm);
	}
	memcpy(prev, w.lv, sizeof(unsigned int) * nt);
	ns = fs_levels(&w, target);   /* (the same levels: the same first samples, the first pass kept) */
	for (t = 0; t < nt; t++)
		finer |= w.lv[t] != prev[t];
	if (finer)
	{
		fs_free_samples(&w);
		if (!fs_pass(&w, ns, cells, lo, hi, beta, thin, 0))
			goto done;
	}
	*out = sf_new(SF_SAMPLES);
	if (*out)
	{
		(*out)->ref = r;
		(*out)->vs = w.h;
		(*out)->lv = w.lv;
		(*out)->first = w.first;
		w.h = NULL;
		w.lv = w.first = NULL;
		ok = 1;
	}
done:
	fs_free_samples(&w);
	qaws_internal_free(w.lv);
	qaws_internal_free(w.first);
	qaws_internal_free(w.emax);
	qaws_internal_free(w.spacing);
	qaws_internal_free(target);
	qaws_internal_free(prev);
	qaws_internal_free(thin);
	return ok ? QAWS_STATUS_OK : QAWS_STATUS_ALLOCATION_FAILURE;
}

void qaws_size_field_destroy(qaws_size_field* f)
{
	if (!f)
		return;
	qaws_internal_free(f->vs);
	qaws_internal_free(f->lv);
	qaws_internal_free(f->first);
	qaws_internal_free(f);
}

qaws_mesh_reference const* qaws_size_field_get_reference(qaws_size_field const* f)
{
	qaws_mesh_reference const* r;
	if (!f)
		return NULL;
	if (f->kind == SF_VERTICES || f->kind == SF_SAMPLES)
		return f->ref;
	if (f->kind == SF_SCALED)
		return qaws_size_field_get_reference(f->a);
	if (f->kind != SF_MIN)
		return NULL;
	r = qaws_size_field_get_reference(f->a);
	return r ? r : qaws_size_field_get_reference(f->b);
}

qaws_status qaws_size_field_evaluate(qaws_scalar* out, qaws_size_field const* f, qaws_mesh_location const* loc, qaws_scalar const* points,
	unsigned int count)
{
	qaws_mesh_reference const* r;
	unsigned int i;
	if (!f || (count && (!points || !out)))
		return QAWS_STATUS_INVALID_ARGUMENT;
	r = qaws_size_field_get_reference(f);
	if (r && count && !loc)
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; r && i < count; i++)
		if (loc[i].patch >= r->nt)
			return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < count; i++)
	{
		double x[3];
		x[0] = (double)points[3 * (size_t)i];
		x[1] = (double)points[3 * (size_t)i + 1];
		x[2] = (double)points[3 * (size_t)i + 2];
		out[i] = (qaws_scalar)sf_eval(f, loc ? &loc[i] : NULL, x);
	}
	return QAWS_STATUS_OK;
}

qaws_status sf_create_scaled(qaws_size_field const* a, double factor, qaws_size_field** out)
{
	*out = sf_new(SF_SCALED);
	if (!*out)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	(*out)->a = a;
	(*out)->length = factor;
	return QAWS_STATUS_OK;
}

/* ---- metric fields ------------------------------------------------------ */

void mf_eval(qaws_metric_field const* f, qaws_mesh_location const* l, double* M)
{
	unsigned int const* tv = f->ref->tv + 3 * (size_t)l->patch;
	double b1 = (double)l->a, b2 = (double)l->b, b0 = 1 - b1 - b2;
	unsigned int k;
	for (k = 0; k < 6; k++)
		M[k] = b0 * f->m[6 * (size_t)tv[0] + k] + b1 * f->m[6 * (size_t)tv[1] + k] + b2 * f->m[6 * (size_t)tv[2] + k];
}

void mf_basis(double const* n, double* t1, double* t2)
{
	double a[3] = { 0, 0, 0 }, l;
	a[fabs(n[0]) < fabs(n[1]) ? (fabs(n[0]) < fabs(n[2]) ? 0 : 2) : (fabs(n[1]) < fabs(n[2]) ? 1 : 2)] = 1;
	t1[0] = n[1] * a[2] - n[2] * a[1];
	t1[1] = n[2] * a[0] - n[0] * a[2];
	t1[2] = n[0] * a[1] - n[1] * a[0];
	l = sqrt(t1[0] * t1[0] + t1[1] * t1[1] + t1[2] * t1[2]);
	t1[0] /= l; t1[1] /= l; t1[2] /= l;
	t2[0] = n[1] * t1[2] - n[2] * t1[1];
	t2[1] = n[2] * t1[0] - n[0] * t1[2];
	t2[2] = n[0] * t1[1] - n[1] * t1[0];
}

double mf_lift(double* M, double const* g, double const* t1, double const* t2, double const* n, double lo, double hi)
{
	double mean = 0.5 * (g[0] + g[2]), r = sqrt(0.25 * (g[0] - g[2]) * (g[0] - g[2]) + g[1] * g[1]), th = 0.5 * atan2(2 * g[1], g[0] - g[2]);
	double mu[2], e[2][3], c = cos(th), s = sin(th), nn;
	unsigned int k, i, j, at = 0;
	mu[0] = fmin(fmax(fabs(mean + r), lo), hi);
	mu[1] = fmin(fmax(fabs(mean - r), lo), hi);
	for (k = 0; k < 3; k++)
	{
		e[0][k] = c * t1[k] + s * t2[k];
		e[1][k] = -s * t1[k] + c * t2[k];
	}
	nn = fmin(mu[0], mu[1]);
	for (i = 0; i < 3; i++)
		for (j = i; j < 3; j++)
			M[at++] = mu[0] * e[0][i] * e[0][j] + mu[1] * e[1][i] * e[1][j] + nn * n[i] * n[j];
	return mu[0] * mu[1];
}

static qaws_status mf_new(qaws_metric_field** out, qaws_mesh_reference const* r, double* m)
{
	*out = (qaws_metric_field*)qaws_internal_zalloc(1, sizeof **out);
	if (!*out)
	{
		qaws_internal_free(m);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	(*out)->ref = r;
	(*out)->m = m;
	return QAWS_STATUS_OK;
}

qaws_status qaws_metric_field_create_vertices(qaws_metric_field** out, qaws_mesh_reference const* r, qaws_scalar const* metrics)
{
	unsigned int nv, i, k;
	double* m;
	if (!out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	if (!r || !metrics)
		return QAWS_STATUS_INVALID_ARGUMENT;
	nv = r->mesh->hm.nv;
	m = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), 6 * (size_t)nv + 1));
	if (!m)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < nv; i++)
	{
		double* M = m + 6 * (size_t)i;
		for (k = 0; k < 6; k++)
			M[k] = (double)metrics[6 * (size_t)i + k];
		/* positive definite: the leading minors positive (Sylvester) */
		if (!(M[0] > 0) || !(M[0] * M[3] - M[1] * M[1] > 0)
			|| !(M[0] * (M[3] * M[5] - M[4] * M[4]) - M[1] * (M[1] * M[5] - M[4] * M[2]) + M[2] * (M[1] * M[4] - M[3] * M[2]) > 0))
		{
			qaws_internal_free(m);
			return QAWS_STATUS_INVALID_ARGUMENT;
		}
	}
	return mf_new(out, r, m);
}

/* vertex v's two rings (v aside) into out, at most cap */
static unsigned int mf_rings(qaws_hm const* h, unsigned int v, unsigned int* out, unsigned int cap)
{
	unsigned int one[64], two[64], n1 = hm_ring(h, v, one, 64), n = 0, i, j, k;
	n1 = n1 < 64 ? n1 : 64;
	for (i = 0; i < n1 && n < cap; i++)
		out[n++] = one[i];
	for (i = 0; i < n1; i++)
	{
		unsigned int n2 = hm_ring(h, one[i], two, 64);
		n2 = n2 < 64 ? n2 : 64;
		for (j = 0; j < n2 && n < cap; j++)
		{
			for (k = 0; k < n && out[k] != two[j]; k++)
				;
			if (k == n && two[j] != v)
				out[n++] = two[j];
		}
	}
	return n;
}

/* the least-squares Hessian (hxx, hxy, hyy) of the values round v in its tangent plane (t1, t2): the fit of
   f(p) - f(v) = g . d + d^T H d / 2 over its two rings, d the offset in the plane; 0 when they do not fix it */
static void mf_hessian(qaws_hm const* h, unsigned int v, qaws_scalar const* f, double const* t1, double const* t2, double* H)
{
	unsigned int ring[256], n = mf_rings(h, v, ring, 256), i, j, k, p;
	double A[5][6], s = 0;
	H[0] = H[1] = H[2] = 0;
	if (n < 5)
		return;
	for (i = 0; i < n; i++)
	{
		double const *a = h->x + 3 * (size_t)v, *b = h->x + 3 * (size_t)ring[i];
		s += sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
	}
	s /= n;   /* (the offsets in units of their mean, for the conditioning) */
	if (!(s > 0))
		return;
	memset(A, 0, sizeof A);
	for (i = 0; i < n; i++)
	{
		double d[3], x, y, row[6];
		for (k = 0; k < 3; k++)
			d[k] = h->x[3 * (size_t)ring[i] + k] - h->x[3 * (size_t)v + k];
		x = (d[0] * t1[0] + d[1] * t1[1] + d[2] * t1[2]) / s;
		y = (d[0] * t2[0] + d[1] * t2[1] + d[2] * t2[2]) / s;
		row[0] = x; row[1] = y; row[2] = 0.5 * x * x; row[3] = x * y; row[4] = 0.5 * y * y;
		row[5] = (double)f[ring[i]] - (double)f[v];
		for (j = 0; j < 5; j++)
			for (k = 0; k < 6; k++)
				A[j][k] += row[j] * row[k];
	}
	/* Gaussian elimination, partial pivoting */
	for (j = 0; j < 5; j++)
	{
		double big = 0;
		p = j;
		for (i = j; i < 5; i++)
			if (fabs(A[i][j]) > big)
			{
				big = fabs(A[i][j]);
				p = i;
			}
		if (big < 1e-12 * (fabs(A[0][0]) + fabs(A[1][1]) + fabs(A[2][2])))
			return;
		for (k = 0; k < 6; k++)
		{
			double t = A[j][k];
			A[j][k] = A[p][k];
			A[p][k] = t;
		}
		for (i = j + 1; i < 5; i++)
		{
			double q = A[i][j] / A[j][j];
			for (k = j; k < 6; k++)
				A[i][k] -= q * A[j][k];
		}
	}
	for (j = 5; j-- > 0;)
	{
		double t = A[j][5];
		for (k = j + 1; k < 5; k++)
			t -= A[j][k] * A[k][5];
		A[j][5] = t / A[j][j];
	}
	H[0] = A[2][5] / (s * s);
	H[1] = A[3][5] / (s * s);
	H[2] = A[4][5] / (s * s);
}

qaws_status qaws_metric_field_create_hessian(qaws_metric_field** out, qaws_mesh_reference const* r, qaws_scalar const* values,
	qaws_hessian_metric_desc const* desc)
{
	unsigned int nv, i;
	double *m, lo, hi, c;
	if (!out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	if (!r || !values || !desc || !(desc->error > 0) || desc->min_length < 0 || desc->max_length < 0)
		return QAWS_STATUS_INVALID_ARGUMENT;
	lo = desc->min_length > 0 ? (double)desc->min_length : 1e-3 * r->extent;
	hi = desc->max_length > 0 ? (double)desc->max_length : 0.1 * r->extent;
	if (!(lo <= hi) || !(lo > 0))
		return QAWS_STATUS_INVALID_ARGUMENT;
	nv = r->mesh->hm.nv;
	for (i = 0; i < nv; i++)
		if (!isfinite((double)values[i]))
			return QAWS_STATUS_INVALID_ARGUMENT;
	m = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), 6 * (size_t)nv + 1));
	if (!m)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	c = (2.0 / 9) / (double)desc->error;
	for (i = 0; i < nv; i++)
	{
		qaws_mesh_location l;
		double x[3], n[3], t1[3], t2[3], H[3] = { 0, 0, 0 };
		if (qaws_mesh_reference_get_vertex_location(&l, r, i) == QAWS_STATUS_OK)
		{
			mr_eval(r, &l, x, n);
			mf_basis(n, t1, t2);
			mf_hessian(&r->mesh->hm, i, values, t1, t2, H);
		}
		else
		{
			n[0] = 0; n[1] = 0; n[2] = 1;
			mf_basis(n, t1, t2);
		}
		H[0] *= c; H[1] *= c; H[2] *= c;
		mf_lift(m + 6 * (size_t)i, H, t1, t2, n, 1 / (hi * hi), 1 / (lo * lo));
	}
	return mf_new(out, r, m);
}

void qaws_metric_field_destroy(qaws_metric_field* f)
{
	if (!f)
		return;
	qaws_internal_free(f->m);
	qaws_internal_free(f);
}

qaws_status qaws_metric_field_evaluate(qaws_scalar* out, qaws_metric_field const* f, qaws_mesh_location const* loc, unsigned int count)
{
	unsigned int i, k;
	if (!f || (count && (!out || !loc)))
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < count; i++)
		if (loc[i].patch >= f->ref->nt)
			return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < count; i++)
	{
		double M[6];
		mf_eval(f, &loc[i], M);
		for (k = 0; k < 6; k++)
			out[6 * (size_t)i + k] = (qaws_scalar)M[k];
	}
	return QAWS_STATUS_OK;
}

/* ---- metric intersection and gradation ---------------------------------- */

static void mf_full(double const* s, double (*A)[3])
{
	A[0][0] = s[0]; A[0][1] = A[1][0] = s[1]; A[0][2] = A[2][0] = s[2];
	A[1][1] = s[3]; A[1][2] = A[2][1] = s[4]; A[2][2] = s[5];
}

/* the symmetric 3 x 3 A diagonalized by Jacobi rotations: A = Q diag(d) Q^T */
static void mf_jacobi(double (*A)[3], double (*Q)[3], double* d)
{
	unsigned int it, p, q, k;
	for (p = 0; p < 3; p++)
		for (q = 0; q < 3; q++)
			Q[p][q] = p == q;
	for (it = 0; it < 50; it++)
	{
		double off = fabs(A[0][1]) + fabs(A[0][2]) + fabs(A[1][2]);
		if (off < 1e-15 * (fabs(A[0][0]) + fabs(A[1][1]) + fabs(A[2][2])))
			break;
		for (p = 0; p < 2; p++)
			for (q = p + 1; q < 3; q++)
			{
				double th, t, c, s;
				if (A[p][q] == 0)
					continue;
				th = (A[q][q] - A[p][p]) / (2 * A[p][q]);
				t = (th >= 0 ? 1 : -1) / (fabs(th) + sqrt(th * th + 1));
				c = 1 / sqrt(t * t + 1);
				s = t * c;
				for (k = 0; k < 3; k++)
				{
					double a = A[k][p], b = A[k][q];
					A[k][p] = c * a - s * b;
					A[k][q] = s * a + c * b;
				}
				for (k = 0; k < 3; k++)
				{
					double a = A[p][k], b = A[q][k];
					A[p][k] = c * a - s * b;
					A[q][k] = s * a + c * b;
				}
				for (k = 0; k < 3; k++)
				{
					double a = Q[k][p], b = Q[k][q];
					Q[k][p] = c * a - s * b;
					Q[k][q] = s * a + c * b;
				}
			}
	}
	for (k = 0; k < 3; k++)
		d[k] = A[k][k];
}

/* the intersection of metrics a and b (six values each): the largest unit ball inside both, by their simultaneous reduction
   (a = L L^T, b seen as L^-1 b L^-T = Q diag(c) Q^T, the result L Q diag(max(1, c)) Q^T L^T); returns whether b was finer
   than a somewhere (by more than rounding) */
static int mf_intersect(double const* a, double const* b, double* out)
{
	double A[3][3], B[3][3], L[3][3] = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } }, X[3][3], C[3][3], Q[3][3], c[3], R[3][3];
	unsigned int i, j, k, at = 0;
	int finer = 0;
	mf_full(a, A);
	mf_full(b, B);
	/* Cholesky of a */
	for (j = 0; j < 3; j++)
	{
		double s = A[j][j];
		for (k = 0; k < j; k++)
			s -= L[j][k] * L[j][k];
		if (!(s > 0))
		{
			memcpy(out, b, 6 * sizeof(double));
			return 1;
		}
		L[j][j] = sqrt(s);
		for (i = j + 1; i < 3; i++)
		{
			double t = A[i][j];
			for (k = 0; k < j; k++)
				t -= L[i][k] * L[j][k];
			L[i][j] = t / L[j][j];
		}
	}
	/* X = L^-1 B, then C = L^-1 X^T (= L^-1 B L^-T, symmetric) */
	for (j = 0; j < 3; j++)
		for (i = 0; i < 3; i++)
		{
			double t = B[i][j];
			for (k = 0; k < i; k++)
				t -= L[i][k] * X[k][j];
			X[i][j] = t / L[i][i];
		}
	for (j = 0; j < 3; j++)
		for (i = 0; i < 3; i++)
		{
			double t = X[j][i];
			for (k = 0; k < i; k++)
				t -= L[i][k] * C[k][j];
			C[i][j] = t / L[i][i];
		}
	mf_jacobi(C, Q, c);
	for (k = 0; k < 3; k++)
	{
		finer |= c[k] > 1 + 1e-9;
		c[k] = c[k] > 1 ? c[k] : 1;
	}
	/* R = L Q, out = R diag(c) R^T */
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
		{
			R[i][j] = 0;
			for (k = 0; k <= i; k++)
				R[i][j] += L[i][k] * Q[k][j];
		}
	for (i = 0; i < 3; i++)
		for (j = i; j < 3; j++)
			out[at++] = R[i][0] * c[0] * R[j][0] + R[i][1] * c[1] * R[j][1] + R[i][2] * c[2] * R[j][2];
	return finer;
}

qaws_status qaws_metric_field_create_intersection(qaws_metric_field** out, qaws_metric_field const* const* fields, unsigned int count)
{
	unsigned int nv, i, k;
	double* m;
	if (!out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	if (!fields || !count || !fields[0])
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (k = 1; k < count; k++)
		if (!fields[k] || fields[k]->ref != fields[0]->ref)
			return QAWS_STATUS_INVALID_ARGUMENT;
	nv = fields[0]->ref->mesh->hm.nv;
	m = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), 6 * (size_t)nv + 1));
	if (!m)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	memcpy(m, fields[0]->m, sizeof(double) * 6 * nv);
	for (k = 1; k < count; k++)
		for (i = 0; i < nv; i++)
		{
			double r[6];
			mf_intersect(m + 6 * (size_t)i, fields[k]->m + 6 * (size_t)i, r);
			memcpy(m + 6 * (size_t)i, r, sizeof r);
		}
	return mf_new(out, fields[0]->ref, m);
}

qaws_status qaws_metric_field_create_graded(qaws_metric_field** out, qaws_metric_field const* field, qaws_scalar gradation)
{
	qaws_hm const* h;
	unsigned int nv, c, pass, k;
	double* m, lb;
	if (!out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	if (!field || !(gradation > 1))
		return QAWS_STATUS_INVALID_ARGUMENT;
	h = &field->ref->mesh->hm;
	nv = h->nv;
	lb = log((double)gradation);
	m = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), 6 * (size_t)nv + 1));
	if (!m)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	memcpy(m, field->m, sizeof(double) * 6 * nv);
	/* each edge's ends: a vertex's metric spread to the other, its sizes grown along the edge, intersected there, until still */
	for (pass = 0; pass < 64; pass++)
	{
		int changed = 0;
		for (c = 0; c < h->nc; c++)
		{
			unsigned int s;
			if (h->cf[c] == HM_NONE || (h->co[c] != HM_NONE && h->co[c] < c))
				continue;
			for (s = 0; s < 2; s++)
			{
				unsigned int i = s ? hm_to(h, c) : h->cv[c], j = s ? h->cv[c] : hm_to(h, c);
				double const *Mi = m + 6 * (size_t)i, *xi = h->x + 3 * (size_t)i, *xj = h->x + 3 * (size_t)j;
				double e[3] = { xj[0] - xi[0], xj[1] - xi[1], xj[2] - xi[2] }, l, eta, G[6], r[6];
				l = sqrt(Mi[0] * e[0] * e[0] + 2 * Mi[1] * e[0] * e[1] + 2 * Mi[2] * e[0] * e[2] + Mi[3] * e[1] * e[1] + 2 * Mi[4] * e[1] * e[2]
					+ Mi[5] * e[2] * e[2]);
				eta = 1 / ((1 + l * lb) * (1 + l * lb));
				for (k = 0; k < 6; k++)
					G[k] = eta * Mi[k];
				if (mf_intersect(m + 6 * (size_t)j, G, r))
				{
					memcpy(m + 6 * (size_t)j, r, sizeof r);
					changed = 1;
				}
			}
		}
		if (!changed)
			break;
	}
	return mf_new(out, field->ref, m);
}
