/*
 * Booleans carrying the operands' attributes (qaws_mesh_boolean.h): each
 * result corner's values interpolated by their rule within its source face.
 */

#include "../internal/qaws_internal_memory.h"
#include "../qaws_mesh_boolean.h"
#include "qaws_mesh_internal.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* an operand mesh as arrays */
typedef struct cy_part
{
	qaws_mesh_3d const* m;
	unsigned int nv, nf, nc, base;   /* base: its first face in the soup */
	qaws_scalar* x;
	unsigned int *fv, *fs, *fo;      /* faces' vertices, sizes, first corners */
} cy_part;

/* an attribute carried: its desc and each part's index of it (QAWS_MESH_NONE: absent there) and values */
typedef struct cy_attr
{
	qaws_attribute_desc d;
	unsigned int* idx;
	qaws_scalar** val;
} cy_attr;

/* a polygon's unit normal (Newell) */
static void cy_normal(cy_part const* p, unsigned int f, double* n)
{
	unsigned int i, s = p->fs[f], k;
	double l;
	n[0] = n[1] = n[2] = 0;
	for (i = 0; i < s; i++)
	{
		qaws_scalar const *a = p->x + 3 * (size_t)p->fv[p->fo[f] + i], *b = p->x + 3 * (size_t)p->fv[p->fo[f] + (i + 1) % s];
		n[0] += ((double)a[1] - (double)b[1]) * ((double)a[2] + (double)b[2]);
		n[1] += ((double)a[2] - (double)b[2]) * ((double)a[0] + (double)b[0]);
		n[2] += ((double)a[0] - (double)b[0]) * ((double)a[1] + (double)b[1]);
	}
	l = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
	for (k = 0; l > 0 && k < 3; k++)
		n[k] /= l;
}

/* p's place in polygon f: the fan triangle (0, i, i + 1) holding it best, its three corners' ranks r and weights w
   (clamped to the triangle: a point off the face takes the nearest values, never extrapolated ones) */
static void cy_locate(cy_part const* p, unsigned int f, double const* q, unsigned int* r, double* w)
{
	unsigned int s = p->fs[f], i, k;
	double best = -HUGE_VAL;
	r[0] = r[1] = r[2] = 0;
	w[0] = 1;
	w[1] = w[2] = 0;
	for (i = 1; i + 1 < s; i++)
	{
		qaws_scalar const *a = p->x + 3 * (size_t)p->fv[p->fo[f]], *b = p->x + 3 * (size_t)p->fv[p->fo[f] + i], *c = p->x + 3 * (size_t)p->fv[p->fo[f] + i + 1];
		double e0[3], e1[3], e2[3], d00 = 0, d01 = 0, d11 = 0, d20 = 0, d21 = 0, den, v, u, m;
		for (k = 0; k < 3; k++)
		{
			e0[k] = (double)b[k] - (double)a[k];
			e1[k] = (double)c[k] - (double)a[k];
			e2[k] = q[k] - (double)a[k];
			d00 += e0[k] * e0[k];
			d01 += e0[k] * e1[k];
			d11 += e1[k] * e1[k];
			d20 += e2[k] * e0[k];
			d21 += e2[k] * e1[k];
		}
		den = d00 * d11 - d01 * d01;
		if (!(den > 1e-300))
			continue;
		v = (d11 * d20 - d01 * d21) / den;
		u = (d00 * d21 - d01 * d20) / den;
		m = fmin(1 - v - u, fmin(v, u));
		if (m > best)
		{
			double t;
			best = m;
			r[0] = 0;
			r[1] = i;
			r[2] = i + 1;
			w[0] = fmax(1 - v - u, 0);
			w[1] = fmax(v, 0);
			w[2] = fmax(u, 0);
			t = w[0] + w[1] + w[2];
			for (k = 0; k < 3; k++)
				w[k] /= t;
		}
	}
}

/* three values mixed by weights in their rule's space */
static void cy_mix(double* out, qaws_scalar const* const* v, double const* w, unsigned int n, unsigned int rule)
{
	unsigned int i, k;
	double x[16];
	memset(out, 0, sizeof(double) * n);
	if (rule == QAWS_ATTRIBUTE_NEAREST)
	{
		i = w[1] > w[0] ? (w[2] > w[1] ? 2 : 1) : (w[2] > w[0] ? 2 : 0);
		for (k = 0; k < n; k++)
			out[k] = (double)v[i][k];
		return;
	}
	for (i = 0; i < 3; i++)
	{
		for (k = 0; k < n; k++)
			x[k] = (double)v[i][k];
		hm_space_in(x, n, rule);
		for (k = 0; k < n; k++)
			out[k] += w[i] * x[k];
	}
	hm_space_out(out, n, rule);
}

static void cy_part_free(cy_part* p)
{
	qaws_internal_free(p->x);
	qaws_internal_free(p->fv);
	qaws_internal_free(p->fs);
	qaws_internal_free(p->fo);
}

qaws_status cy_boolean(qaws_mesh_3d** out, qaws_boolean_report* report, qaws_mesh_3d const* const* a, unsigned int na,
	qaws_mesh_3d const* const* b, unsigned int nb, qaws_boolean_op op, qaws_intersect_desc const* desc)
{
	unsigned int np = na + nb, i, j, k, nattr = 0, capattr = 0, nt = 0, src_at;
	cy_part* P = NULL;
	cy_attr* A = NULL;
	qaws_scalar *sx[2] = { NULL, NULL }, *rx = NULL, *rsrc = NULL, *ov = NULL;
	unsigned int *sfv[2] = { NULL, NULL }, *sfs[2] = { NULL, NULL }, *rt = NULL;
	unsigned char* turned = NULL;
	qaws_mesh_3d_desc d[2];
	qaws_status st = QAWS_STATUS_OK;
	if (!out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	if (!na || (nb && !b))
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < np; i++)
		if (!(i < na ? a[i] : b[i - na]))
			return QAWS_STATUS_INVALID_ARGUMENT;
	P = (cy_part*)qaws_internal_zalloc(np, sizeof *P);
	if (!P)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	/* the parts as arrays, A's then B's faces (the soups the Boolean's "source" counts in) */
	for (i = 0, k = 0; i < np; i++)
	{
		cy_part* p = &P[i];
		p->m = i < na ? a[i] : b[i - na];
		p->nv = qaws_mesh_3d_get_vertex_count(p->m);
		p->nf = qaws_mesh_3d_get_face_count(p->m);
		p->nc = qaws_mesh_3d_get_corner_count(p->m);
		p->base = k;
		k += p->nf;
		p->x = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)p->nv + 1));
		p->fv = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)p->nc + 1));
		p->fs = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)p->nf + 1));
		p->fo = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)p->nf + 1));
		if (!p->x || !p->fv || !p->fs || !p->fo)
			goto oom;
		qaws_mesh_3d_get_positions(p->x, p->m);
		qaws_mesh_3d_get_faces(p->fv, p->fs, p->m);
		for (j = 0, p->fo[0] = 0; j + 1 < p->nf; j++)
			p->fo[j + 1] = p->fo[j] + p->fs[j];
	}
	/* the attributes, matched by name (unnamed ones never match; an input's "source" gives way to the result's) */
	for (i = 0; i < np; i++)
	{
		unsigned int count = qaws_mesh_3d_get_attribute_count(P[i].m);
		for (j = 0; j < count; j++)
		{
			qaws_attribute_desc ad;
			cy_attr* t = NULL;
			qaws_mesh_3d_get_attribute_desc(&ad, P[i].m, j);
			if (ad.name && !strcmp(ad.name, "source"))
				continue;
			for (k = 0; ad.name && k < nattr && !t; k++)
				if (A[k].d.name && !strcmp(A[k].d.name, ad.name))
					t = &A[k];
			if (t && (t->d.element != ad.element || t->d.components != ad.components || t->d.rule != ad.rule))
			{
				st = QAWS_STATUS_INVALID_ARGUMENT;
				goto done;
			}
			if (!t)
			{
				if (nattr == capattr)
				{
					unsigned int cap = capattr ? 2 * capattr : 8;
					cy_attr* g = (cy_attr*)qaws_internal_realloc(A, sizeof(cy_attr) * cap);
					if (!g)
						goto oom;
					A = g;
					capattr = cap;
				}
				t = &A[nattr++];
				t->d = ad;
				t->idx = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * np);
				t->val = (qaws_scalar**)qaws_internal_zalloc(np, sizeof(qaws_scalar*));
				if (!t->idx || !t->val)
					goto oom;
				for (k = 0; k < np; k++)
					t->idx[k] = QAWS_MESH_NONE;
			}
			if (t->idx[i] == QAWS_MESH_NONE)
			{
				size_t count_el = ad.element == QAWS_MESH_VERTEX ? P[i].nv : ad.element == QAWS_MESH_CORNER ? P[i].nc : P[i].nf;
				t->idx[i] = j;
				t->val[i] = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul3(sizeof(qaws_scalar), ad.components, count_el + 1));
				if (!t->val[i])
					goto oom;
				qaws_mesh_3d_get_attribute(t->val[i], P[i].m, j);
			}
		}
	}
	/* the two soups */
	for (j = 0; j < 2; j++)
	{
		unsigned int lo = j ? na : 0, hi = j ? np : na, nv = 0, nf = 0, nc = 0, q;
		memset(&d[j], 0, sizeof d[j]);
		if (lo == hi)
			continue;
		for (i = lo; i < hi; i++)
		{
			nv += P[i].nv;
			nf += P[i].nf;
			nc += P[i].nc;
		}
		sx[j] = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)nv + 1));
		sfv[j] = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nc + 1));
		sfs[j] = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nf + 1));
		if (!sx[j] || !sfv[j] || !sfs[j])
			goto oom;
		for (i = lo, nv = nf = nc = 0; i < hi; i++)
		{
			memcpy(sx[j] + 3 * (size_t)nv, P[i].x, sizeof(qaws_scalar) * 3 * P[i].nv);
			for (q = 0; q < P[i].nc; q++)
				sfv[j][nc + q] = P[i].fv[q] + nv;
			memcpy(sfs[j] + nf, P[i].fs, sizeof(unsigned int) * P[i].nf);
			nv += P[i].nv;
			nf += P[i].nf;
			nc += P[i].nc;
		}
		d[j].positions = sx[j];
		d[j].vertex_count = nv;
		d[j].face_vertices = sfv[j];
		d[j].face_sizes = sfs[j];
		d[j].face_count = nf;
	}
	st = sx_boolean(out, report, &d[0], nb ? &d[1] : NULL, op, desc, nattr > 0, nattr ? &turned : NULL);
	if (st != QAWS_STATUS_OK || !*out || !nattr)
		goto done;
	/* each result corner's values within its source face */
	nt = qaws_mesh_3d_get_face_count(*out);
	rx = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)qaws_mesh_3d_get_vertex_count(*out) + 1));
	rt = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)nt + 1));
	rsrc = (qaws_scalar*)qaws_internal_malloc(sizeof(qaws_scalar) * ((size_t)nt + 1));
	ov = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 48, (size_t)nt + 1));   /* (3 corners of 16 components at most) */
	if (!rx || !rt || !rsrc || !ov || qaws_mesh_3d_find_attribute(&src_at, *out, "source") != QAWS_STATUS_OK)
		goto oom;
	qaws_mesh_3d_get_positions(rx, *out);
	qaws_mesh_3d_get_faces(rt, NULL, *out);
	qaws_mesh_3d_get_attribute(rsrc, *out, src_at);
	for (k = 0; k < nattr && st == QAWS_STATUS_OK; k++)
	{
		cy_attr const* t = &A[k];
		unsigned int n = t->d.components, f, id;
		qaws_attribute_desc od = t->d;
		for (f = 0; f < nt; f++)
		{
			unsigned int s = (unsigned int)rsrc[f], pi = 0, sf, c;
			cy_part const* p;
			double nrm[3], sign = t->d.rule == QAWS_ATTRIBUTE_NORMAL && turned[f] ? -1 : 1;
			while (pi + 1 < np && P[pi + 1].base <= s)
				pi++;
			p = &P[pi];
			sf = s - p->base;
			if (t->idx[pi] == QAWS_MESH_NONE && t->d.rule == QAWS_ATTRIBUTE_NORMAL)
				cy_normal(p, sf, nrm);   /* (absent on this operand: its face's normal) */
			for (c = 0; c < (t->d.element == QAWS_MESH_FACE ? 1u : 3u); c++)
			{
				qaws_scalar* o = ov + (size_t)n * (t->d.element == QAWS_MESH_FACE ? f : 3 * (size_t)f + c);
				double y[16];
				if (t->idx[pi] == QAWS_MESH_NONE)
				{
					for (j = 0; j < n; j++)
						y[j] = t->d.rule == QAWS_ATTRIBUTE_NORMAL ? nrm[j] : 0;
				}
				else if (t->d.element == QAWS_MESH_FACE)
				{
					for (j = 0; j < n; j++)
						y[j] = (double)t->val[pi][(size_t)n * sf + j];
				}
				else
				{
					double q[3], w[3];
					unsigned int r[3];
					qaws_scalar const* v[3];
					for (j = 0; j < 3; j++)
						q[j] = (double)rx[3 * (size_t)rt[3 * (size_t)f + c] + j];
					cy_locate(p, sf, q, r, w);
					for (j = 0; j < 3; j++)
					{
						unsigned int pc = p->fo[sf] + r[j];
						v[j] = t->val[pi] + (size_t)n * (t->d.element == QAWS_MESH_VERTEX ? p->fv[pc] : pc);
					}
					cy_mix(y, v, w, n, t->d.rule);
				}
				for (j = 0; j < n; j++)
					o[j] = (qaws_scalar)(sign * y[j]);
			}
		}
		/* vertex attributes per corner: each side of a cut keeps its own */
		od.element = t->d.element == QAWS_MESH_FACE ? QAWS_MESH_FACE : QAWS_MESH_CORNER;
		st = qaws_mesh_3d_add_attribute(&id, *out, ov, &od);
	}
	goto done;
oom:
	st = QAWS_STATUS_ALLOCATION_FAILURE;
done:
	if (st != QAWS_STATUS_OK && out && *out)
	{
		qaws_mesh_3d_destroy(*out);
		*out = NULL;
	}
	for (i = 0; i < np; i++)
		cy_part_free(&P[i]);
	for (k = 0; k < nattr; k++)
	{
		for (i = 0; A[k].val && i < np; i++)
			qaws_internal_free(A[k].val[i]);
		qaws_internal_free(A[k].val);
		qaws_internal_free(A[k].idx);
	}
	qaws_internal_free(A);
	qaws_internal_free(P);
	for (j = 0; j < 2; j++)
	{
		qaws_internal_free(sx[j]);
		qaws_internal_free(sfv[j]);
		qaws_internal_free(sfs[j]);
	}
	qaws_internal_free(rx); qaws_internal_free(rt); qaws_internal_free(rsrc); qaws_internal_free(ov); qaws_internal_free(turned);
	return st;
}

qaws_status qaws_mesh_3d_boolean(qaws_mesh_3d** out_mesh, qaws_boolean_report* report, qaws_mesh_3d const* a, qaws_mesh_3d const* b,
	qaws_boolean_op op, qaws_intersect_desc const* desc)
{
	if (out_mesh)
		*out_mesh = NULL;
	if (!a)
		return QAWS_STATUS_INVALID_ARGUMENT;
	return cy_boolean(out_mesh, report, &a, 1, b ? &b : NULL, b ? 1 : 0, op, desc);
}
