/*
 * Size fields: uniform, from a reference's curvature with tolerances and
 * gradation, per reference vertex, by callback, the smaller of two.
 */

#include "qaws_mesh_internal.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

enum { SF_UNIFORM, SF_VERTICES, SF_CALLBACK, SF_MIN };

struct qaws_size_field
{
	int kind;
	double length;
	qaws_mesh_reference const* ref;
	double* vs;                     /* per vertex of the reference's mesh */
	qaws_scalar (*fn)(void* user, qaws_scalar const* point);
	void* user;
	qaws_size_field const *a, *b;
};

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
	case SF_CALLBACK:
	{
		qaws_scalar p[3];
		p[0] = (qaws_scalar)x[0]; p[1] = (qaws_scalar)x[1]; p[2] = (qaws_scalar)x[2];
		return (double)f->fn(f->user, p);
	}
	default:
	{
		double u = sf_eval(f->a, l, x), v = sf_eval(f->b, l, x);
		return u < v ? u : v;
	}
	}
}

static qaws_size_field* sf_new(int kind)
{
	qaws_size_field* f = (qaws_size_field*)calloc(1, sizeof *f);
	if (f)
		f->kind = kind;
	return f;
}

qaws_status qaws_size_field_create_uniform(qaws_scalar length, qaws_size_field** out)
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

qaws_status qaws_size_field_create_callback(qaws_scalar (*size)(void* user, qaws_scalar const* point), void* user, qaws_size_field** out)
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

qaws_status qaws_size_field_create_min(qaws_size_field const* a, qaws_size_field const* b, qaws_size_field** out)
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
		free(vs);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	(*out)->ref = r;
	(*out)->vs = vs;
	return QAWS_STATUS_OK;
}

qaws_status qaws_size_field_create_vertices(qaws_mesh_reference const* r, qaws_scalar const* sizes, qaws_size_field** out)
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
	vs = (double*)malloc(sizeof(double) * (nv + 1));
	if (!vs)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < nv; i++)
		vs[i] = (double)sizes[i];
	return sf_vertices(r, vs, out);
}

/* sizes no faster growing than beta per unit of distance along the mesh's
   edges: Dijkstra from every vertex at once, h(v) <= h(u) + beta |uv| */
static int sf_gradation(qaws_hm const* m, double* h, double beta)
{
	unsigned int n = m->nv, *heap = (unsigned int*)malloc(sizeof(unsigned int) * (n + 1)), *pos = (unsigned int*)malloc(sizeof(unsigned int) * (n + 1));
	unsigned int count = n, i, ring[256];
	if (!heap || !pos)
	{
		free(heap);
		free(pos);
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
		nr = hm_ring(m, u, ring, 256);
		for (k = 0; k < nr && k < 256; k++)
		{
			unsigned int w = ring[k], j;
			double d = 0, cand;
			for (j = 0; j < 3; j++)
				d += (m->x[3 * (size_t)u + j] - m->x[3 * (size_t)w + j]) * (m->x[3 * (size_t)u + j] - m->x[3 * (size_t)w + j]);
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
	free(heap);
	free(pos);
	return 1;
}

qaws_status qaws_size_field_create_curvature(qaws_mesh_reference const* r, qaws_curvature_size_desc const* desc, qaws_size_field** out)
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
	kmax = (double*)calloc((size_t)m->nv + 1, sizeof(double));
	vs = (double*)malloc(sizeof(double) * ((size_t)m->nv + 1));
	if (!kmax || !vs)
	{
		free(kmax);
		free(vs);
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
	free(kmax);
	if (desc->gradation > 0 && !sf_gradation(m, vs, (double)desc->gradation))
	{
		free(vs);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	return sf_vertices(r, vs, out);
}

void qaws_size_field_destroy(qaws_size_field* f)
{
	if (!f)
		return;
	free(f->vs);
	free(f);
}

qaws_mesh_reference const* qaws_size_field_get_reference(qaws_size_field const* f)
{
	qaws_mesh_reference const* r;
	if (!f)
		return NULL;
	if (f->kind == SF_VERTICES)
		return f->ref;
	if (f->kind != SF_MIN)
		return NULL;
	r = qaws_size_field_get_reference(f->a);
	return r ? r : qaws_size_field_get_reference(f->b);
}

qaws_status qaws_size_field_evaluate(qaws_size_field const* f, qaws_mesh_location const* loc, qaws_scalar const* points, unsigned int count,
	qaws_scalar* out)
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
