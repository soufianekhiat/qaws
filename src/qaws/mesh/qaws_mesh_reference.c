/*
 * Mesh references: a copied mesh, its faces fanned into triangles with
 * their adjacency, a patch per triangle (PN or LINEAR), a grid over the
 * patches' control nets; evaluation, projection, walking, attributes.
 */

#include "qaws_mesh_internal.h"
#include "../internal/qaws_internal_memory.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

static double mr_dot(double const* a, double const* b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static double mr_unit(double* v)
{
	double l = sqrt(mr_dot(v, v));
	if (l > 0)
	{
		v[0] /= l; v[1] /= l; v[2] /= l;
	}
	return l;
}

static void mr_cross(double const* a, double const* b, double* out)
{
	out[0] = a[1] * b[2] - a[2] * b[1];
	out[1] = a[2] * b[0] - a[0] * b[2];
	out[2] = a[0] * b[1] - a[1] * b[0];
}

/* the angle at corner c of its face */
static double mr_corner_angle(qaws_hm const* m, unsigned int c)
{
	double const *p = m->x + 3 * (size_t)m->cv[c], *a = m->x + 3 * (size_t)m->cv[hm_prev(m, c)], *b = m->x + 3 * (size_t)hm_to(m, c);
	double u[3], w[3], d;
	unsigned int k;
	for (k = 0; k < 3; k++)
	{
		u[k] = a[k] - p[k];
		w[k] = b[k] - p[k];
	}
	if (mr_unit(u) == 0 || mr_unit(w) == 0)
		return 0;
	d = mr_dot(u, w);
	return acos(d > 1 ? 1 : d < -1 ? -1 : d);
}

/* Corner normals: the faces' normals round each vertex weighted by their
   angle, one normal per run of faces between feature edges. */
static void mr_corner_normals(qaws_hm const* m, double* cn, double* fnrm)
{
	unsigned int f, v;
	for (f = 0; f < m->nf; f++)
	{
		hm_face_normal(m, f, fnrm + 3 * (size_t)f);
		mr_unit(fnrm + 3 * (size_t)f);
	}
	for (v = 0; v < m->nv; v++)
	{
		unsigned int start = m->vc[v], k, first, guard = 0;
		if (start == HM_NONE)
			continue;
		/* an interior vertex on features: start right after a feature edge */
		if (!hm_is_border_vertex(m, v))
			for (k = start; guard++ < 1000;)
			{
				if (m->cflag[k] & HM_FEATURE)
				{
					/* k's edge is a feature: the run starts at the corner past it, clockwise */
					start = hm_cw(m, k) != HM_NONE ? k : start;
					break;
				}
				k = hm_round(m, v, k);
				if (k == HM_NONE)
					break;
			}
		/* runs counter-clockwise from start: a feature edge between corners k and ccw(k) is prev(k)'s edge */
		first = start;
		k = first;
		while (k != HM_NONE)
		{
			double acc[3] = { 0, 0, 0 };
			unsigned int run = k, q, n;
			/* gather the run */
			for (q = k;;)
			{
				double w = mr_corner_angle(m, q);
				unsigned int nx;
				acc[0] += w * fnrm[3 * (size_t)m->cf[q]];
				acc[1] += w * fnrm[3 * (size_t)m->cf[q] + 1];
				acc[2] += w * fnrm[3 * (size_t)m->cf[q] + 2];
				nx = hm_ccw(m, q);
				if (nx == HM_NONE || nx == first || (m->cflag[hm_prev(m, q)] & HM_FEATURE))
				{
					n = nx;
					break;
				}
				q = nx;
			}
			if (mr_unit(acc) == 0)
				memcpy(acc, fnrm + 3 * (size_t)m->cf[run], sizeof acc);
			/* write it to the run's corners */
			for (q = run;;)
			{
				unsigned int nx;
				memcpy(cn + 3 * (size_t)q, acc, sizeof acc);
				nx = hm_ccw(m, q);
				if (nx == HM_NONE || nx == first || (m->cflag[hm_prev(m, q)] & HM_FEATURE))
					break;
				q = nx;
			}
			k = n == first ? HM_NONE : n;
		}
	}
}

/* per vertex on a feature line (two feature or border edges, not a corner,
   turning by less than 60 degrees): the line's tangent there, from its two
   neighbours along the line; 0 elsewhere */
static void mr_crease_tangents(qaws_hm const* m, double* tan, unsigned char* has)
{
	unsigned int v, c;
	for (v = 0; v < m->nv; v++)
	{
		unsigned int w[2], n = 0, last = HM_NONE;
		has[v] = 0;
		if (m->vc[v] == HM_NONE || (m->vflag[v] & HM_CORNER))
			continue;
		for (c = m->vc[v]; c != HM_NONE; c = hm_round(m, v, c))
		{
			int line = (m->cflag[c] & HM_FEATURE) || m->co[c] == HM_NONE;
			if (line && n < 2)
				w[n] = hm_to(m, c);
			n += line;
			last = c;
		}
		if (last != HM_NONE && hm_is_border_vertex(m, v))
		{
			/* (the border edge into v is a line edge) */
			if (n < 2)
				w[n] = m->cv[hm_prev(m, last)];
			n++;
		}
		if (n != 2)
			continue;
		{
			double const *p = m->x + 3 * (size_t)v, *a = m->x + 3 * (size_t)w[0], *b = m->x + 3 * (size_t)w[1];
			double u[3], s[3];
			unsigned int k;
			for (k = 0; k < 3; k++)
			{
				u[k] = p[k] - a[k];
				s[k] = b[k] - p[k];
			}
			if (mr_dot(u, s) < 0.5 * sqrt(mr_dot(u, u) * mr_dot(s, s)))
				continue;   /* (turning by more than 60 degrees: a corner) */
		}
		tan[3 * (size_t)v] = m->x[3 * (size_t)w[1]] - m->x[3 * (size_t)w[0]];
		tan[3 * (size_t)v + 1] = m->x[3 * (size_t)w[1] + 1] - m->x[3 * (size_t)w[0] + 1];
		tan[3 * (size_t)v + 2] = m->x[3 * (size_t)w[1] + 2] - m->x[3 * (size_t)w[0] + 2];
		has[v] = mr_unit(tan + 3 * (size_t)v) > 0;
	}
}

void qaws_mesh_reference_destroy(qaws_mesh_reference* r)
{
	if (!r)
		return;
	qaws_mesh_3d_destroy(r->mesh);
	qaws_internal_free(r->tv); qaws_internal_free(r->tc); qaws_internal_free(r->tf); qaws_internal_free(r->ta); qaws_internal_free(r->tcrease); qaws_internal_free(r->treg);
	qaws_internal_free(r->patch); qaws_internal_free(r->box); qaws_internal_free(r->vloc);
	qaws_internal_grid_destroy(r->grid);
	qaws_internal_free(r->seamk);
	qaws_internal_free(r->loop);
	qaws_internal_free(r->loopi);
	qaws_internal_free(r->mls);
	qaws_internal_free(r->mlsi);
	qaws_internal_free(r->mlsd);
	qaws_internal_free(r);
}

/* SUBDIVISION: triangle t's neighbourhood, every triangle round its corners (each corner's fan walked both ways across ta):
   triangles tl, local vertices vl (t's corners first), local triangles lt, crease and border edges cl; 0 past MI_LOOP_MAX */
static int mr_loop_gather(qaws_mesh_reference const* r, unsigned int t, unsigned int* tl, unsigned int* vl, unsigned int* lt, unsigned int* cl,
	unsigned int* cnt)
{
	unsigned int ntl = 1, nvl = 0, ncl = 0, k, d, i, j, q;
	tl[0] = t;
	for (k = 0; k < 3; k++)
	{
		unsigned int v = r->tv[3 * t + k];
		for (d = 0; d < 2; d++)
		{
			unsigned int u = t, e = d ? (k + 2) % 3 : k, steps = 0, closed = 0;
			for (;;)
			{
				unsigned int nb = r->ta[3 * u + e], u2, ue, kk;
				if (nb == HM_NONE)
					break;
				u2 = nb / 3;
				ue = nb % 3;
				if (u2 == t)
				{
					closed = 1;
					break;
				}
				for (i = 0; i < ntl && tl[i] != u2; i++)
					;
				if (i == ntl)
				{
					if (ntl == MI_LOOP_MAX)
						return 0;
					tl[ntl++] = u2;
				}
				for (kk = 0; kk < 3 && r->tv[3 * u2 + kk] != v; kk++)
					;
				if (kk == 3 || ++steps > MI_LOOP_MAX)
					break;
				e = kk == ue ? (kk + 2) % 3 : kk;
				u = u2;
			}
			if (closed)
				break;
		}
	}
	for (i = 0; i < ntl; i++)
		for (k = 0; k < 3; k++)
		{
			unsigned int v = r->tv[3 * tl[i] + k];
			for (j = 0; j < nvl && vl[j] != v; j++)
				;
			if (j == nvl)
			{
				if (nvl == MI_LOOP_MAX)
					return 0;
				vl[nvl++] = v;
			}
			lt[3 * i + k] = j;
		}
	for (i = 0; i < ntl; i++)
		for (k = 0; k < 3; k++)
			if (r->tcrease[3 * tl[i] + k])
			{
				unsigned int a = lt[3 * i + k], b = lt[3 * i + (k + 1) % 3];
				for (q = 0; q < ncl && !((cl[2 * q] == a && cl[2 * q + 1] == b) || (cl[2 * q] == b && cl[2 * q + 1] == a)); q++)
					;
				if (q == ncl)
				{
					if (ncl == MI_LOOP_MAX)
						return 0;
					cl[2 * ncl] = a;
					cl[2 * ncl + 1] = b;
					ncl++;
				}
			}
	/* the corners fixed: a feature corner, a vertex of one face */
	cnt[3] = 0;
	for (k = 0; k < 3; k++)
	{
		unsigned int v = r->tv[3 * t + k], f0 = HM_NONE, one = 1;
		for (i = 0; i < ntl; i++)
			if (r->tv[3 * tl[i]] == v || r->tv[3 * tl[i] + 1] == v || r->tv[3 * tl[i] + 2] == v)
			{
				if (f0 == HM_NONE)
					f0 = r->tf[tl[i]];
				else if (r->tf[tl[i]] != f0)
					one = 0;
			}
		if (r->mesh->hm.vflag && (r->mesh->hm.vflag[v] & HM_CORNER))
			cnt[3] |= 1u << k;
		if (one)
			cnt[3] |= 8u << k;
	}
	cnt[0] = nvl;
	cnt[1] = ntl;
	cnt[2] = ncl;
	return 1;
}

/* MLS: the vertices' normals (the given ones, or their corners' averaged: creases rounded) and kernel radii (scale times the
   mean length of their edges); per triangle its candidates, walked from its corners across edges while a vertex's kernel
   reaches the triangle's box widened by its longest side */
static qaws_status mr_mls_setup(qaws_mesh_reference* r, double const* cn, qaws_mesh_reference_desc const* desc)
{
	qaws_hm const* m = &r->mesh->hm;
	unsigned int nv = m->nv, i, k, t, *off = NULL, *adj = NULL, *cnt = NULL, *stamp = NULL, *queue = NULL, cap = 0, used = 0;
	double *nrm, *rad, *len = NULL, scale = desc && desc->mls_radius > 0 ? (double)desc->mls_radius : 2.5;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	r->mlsd = (double*)qaws_internal_zalloc(4 * (size_t)nv + 1, sizeof(double));
	r->mls = (mi_mls*)qaws_internal_zalloc((size_t)r->nt + 1, sizeof(mi_mls));
	off = (unsigned int*)qaws_internal_scratch_zalloc((size_t)nv + 2, sizeof(unsigned int));
	cnt = (unsigned int*)qaws_internal_scratch_zalloc((size_t)nv + 1, sizeof(unsigned int));
	stamp = (unsigned int*)qaws_internal_scratch_zalloc((size_t)nv + 1, sizeof(unsigned int));
	queue = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)nv + 1));
	len = (double*)qaws_internal_scratch_zalloc((size_t)nv + 1, sizeof(double));
	adj = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 6, (size_t)r->nt + 1));
	if (!r->mlsd || !r->mls || !off || !cnt || !stamp || !queue || !len || !adj)
		goto done;
	nrm = r->mlsd;
	rad = r->mlsd + 3 * (size_t)nv;
	for (i = 0; i < m->nc; i++)
		if (m->cf[i] != HM_NONE)
			for (k = 0; k < 3; k++)
				nrm[3 * (size_t)m->cv[i] + k] += desc && desc->normals ? (double)desc->normals[3 * (size_t)m->cv[i] + k] : cn[3 * (size_t)i + k];
	for (i = 0; i < nv; i++)
		mr_unit(nrm + 3 * (size_t)i);
	/* edges per vertex (each triangle side from both ends: inner sides twice, the mean the same) */
	for (t = 0; t < r->nt; t++)
		for (k = 0; k < 3; k++)
		{
			unsigned int a = r->tv[3 * t + k], b = r->tv[3 * t + (k + 1) % 3];
			double const *pa = m->x + 3 * (size_t)a, *pb = m->x + 3 * (size_t)b;
			double l = sqrt((pa[0] - pb[0]) * (pa[0] - pb[0]) + (pa[1] - pb[1]) * (pa[1] - pb[1]) + (pa[2] - pb[2]) * (pa[2] - pb[2]));
			len[a] += l;
			len[b] += l;
			off[a + 1]++;
			off[b + 1]++;
		}
	for (i = 0; i < nv; i++)
		rad[i] = off[i + 1] ? scale * len[i] / off[i + 1] : 0;
	for (i = 0; i < nv; i++)
		off[i + 1] += off[i];
	for (t = 0; t < r->nt; t++)
		for (k = 0; k < 3; k++)
		{
			unsigned int a = r->tv[3 * t + k], b = r->tv[3 * t + (k + 1) % 3];
			adj[off[a] + cnt[a]++] = b;
			adj[off[b] + cnt[b]++] = a;
		}
	for (t = 0; t < r->nt; t++)
	{
		double lo[3], hi[3], L = 0;
		unsigned int head = 0, tail = 0, q;
		for (k = 0; k < 3; k++)
		{
			unsigned int a = r->tv[3 * t + k], b = r->tv[3 * t + (k + 1) % 3];
			double const *pa = m->x + 3 * (size_t)a, *pb = m->x + 3 * (size_t)b;
			L = fmax(L, sqrt((pa[0] - pb[0]) * (pa[0] - pb[0]) + (pa[1] - pb[1]) * (pa[1] - pb[1]) + (pa[2] - pb[2]) * (pa[2] - pb[2])));
		}
		for (k = 0; k < 3; k++)
		{
			lo[k] = fmin(m->x[3 * (size_t)r->tv[3 * t] + k], fmin(m->x[3 * (size_t)r->tv[3 * t + 1] + k], m->x[3 * (size_t)r->tv[3 * t + 2] + k])) - L;
			hi[k] = fmax(m->x[3 * (size_t)r->tv[3 * t] + k], fmax(m->x[3 * (size_t)r->tv[3 * t + 1] + k], m->x[3 * (size_t)r->tv[3 * t + 2] + k])) + L;
		}
		for (k = 0; k < 3; k++)
			if (stamp[r->tv[3 * t + k]] != t + 1)
			{
				stamp[r->tv[3 * t + k]] = t + 1;
				queue[tail++] = r->tv[3 * t + k];
			}
		while (head < tail)
		{
			unsigned int v = queue[head++];
			for (q = off[v]; q < off[v + 1]; q++)
			{
				unsigned int w = adj[q];
				double d2 = 0;
				if (stamp[w] == t + 1)
					continue;
				stamp[w] = t + 1;
				for (k = 0; k < 3; k++)
				{
					double c = m->x[3 * (size_t)w + k], e = c < lo[k] ? lo[k] - c : c > hi[k] ? c - hi[k] : 0;
					d2 += e * e;
				}
				if (d2 < rad[w] * rad[w])
					queue[tail++] = w;
			}
		}
		if (used + tail > cap)
		{
			unsigned int nc2 = 2 * (used + tail) + 64, *g = (unsigned int*)qaws_internal_realloc(r->mlsi, qaws_internal_mul2(sizeof(unsigned int), nc2));
			if (!g)
				goto done;
			r->mlsi = g;
			cap = nc2;
		}
		memcpy(r->mlsi + used, queue, sizeof(unsigned int) * tail);
		r->mls[t].ncand = tail;
		r->mls[t].cand = NULL;
		r->mls[t].x = m->x;
		r->mls[t].nrm = nrm;
		r->mls[t].rad = rad;
		for (k = 0; k < 3; k++)
			r->mls[t].v[k] = r->tv[3 * t + k];
		/* (the pool's offset for now: the pointers once it stops moving) */
		r->mls[t].cand = (unsigned int const*)(size_t)used;
		used += tail;
	}
	for (t = 0; t < r->nt; t++)
		r->mls[t].cand = r->mlsi + (size_t)r->mls[t].cand;
	st = QAWS_STATUS_OK;
done:
	qaws_internal_free(off); qaws_internal_free(cnt); qaws_internal_free(stamp); qaws_internal_free(queue); qaws_internal_free(len);
	qaws_internal_free(adj);
	return st;
}

/* SUBDIVISION: every triangle's neighbourhood into r->loop (two passes: sizes, then the pool) */
static qaws_status mr_loop_rings(qaws_mesh_reference* r)
{
	unsigned int buf[MI_LOOP_MAX * 7], cnt[4], t, pass;
	size_t total = 0, at = 0;
	r->loop = (mi_loop*)qaws_internal_zalloc((size_t)r->nt + 1, sizeof(mi_loop));
	if (!r->loop)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (pass = 0; pass < 2; pass++)
	{
		if (pass && !(r->loopi = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), total + 1))))
			return QAWS_STATUS_ALLOCATION_FAILURE;
		for (t = 0; t < r->nt; t++)
		{
			unsigned int *tl = buf, *vl = buf + MI_LOOP_MAX, *lt = buf + 2 * MI_LOOP_MAX, *cl = buf + 5 * MI_LOOP_MAX;
			mi_loop* L = &r->loop[t];
			if (!mr_loop_gather(r, t, tl, vl, lt, cl, cnt))
				return QAWS_STATUS_UNSUPPORTED_OPERATION;   /* (a vertex of more than about 50 triangles) */
			if (!pass)
			{
				total += cnt[0] + 3 * (size_t)cnt[1] + 2 * (size_t)cnt[2];
				continue;
			}
			memcpy(r->loopi + at, vl, sizeof(unsigned int) * cnt[0]);
			memcpy(r->loopi + at + cnt[0], lt, sizeof(unsigned int) * 3 * cnt[1]);
			memcpy(r->loopi + at + cnt[0] + 3 * cnt[1], cl, sizeof(unsigned int) * 2 * cnt[2]);
			L->x = r->mesh->hm.x;
			L->v = r->loopi + at;
			L->t = L->v + cnt[0];
			L->c = L->t + 3 * cnt[1];
			L->nv = cnt[0];
			L->nt = cnt[1];
			L->nc = cnt[2];
			L->fixed = cnt[3];
			at += cnt[0] + 3 * (size_t)cnt[1] + 2 * (size_t)cnt[2];
		}
	}
	return QAWS_STATUS_OK;
}

int mr_seam(qaws_hm const* h, unsigned int c, unsigned int o)
{
	unsigned int a, k;
	for (a = 0; a < h->na; a++)
	{
		hm_attr const* at = &h->attr[a];
		if (at->element == QAWS_MESH_CORNER)
		{
			double const *p0 = hm_attr_at(h, a, c), *p1 = hm_attr_at(h, a, hm_next(h, o)), *q0 = hm_attr_at(h, a, hm_next(h, c)), *q1 = hm_attr_at(h, a, o);
			for (k = 0; k < at->n; k++)
				if (fabs(p0[k] - p1[k]) > 1e-12 * (1 + fabs(p0[k])) || fabs(q0[k] - q1[k]) > 1e-12 * (1 + fabs(q0[k])))
					return 1;
		}
		else if (at->element == QAWS_MESH_FACE)
		{
			double const *p0 = hm_attr_at(h, a, h->cf[c]), *p1 = hm_attr_at(h, a, h->cf[o]);
			for (k = 0; k < at->n; k++)
				if (fabs(p0[k] - p1[k]) > 1e-12 * (1 + fabs(p0[k])))
					return 1;
		}
		else if (at->element == QAWS_MESH_EDGE)
			for (k = 0; k < at->n; k++)
				if (hm_attr_at(h, a, c)[k] != 0)
					return 1;
	}
	return 0;
}

static int mr_key_cmp(void const* p, void const* q)
{
	unsigned long long a = *(unsigned long long const*)p, b = *(unsigned long long const*)q;
	return a < b ? -1 : a > b;
}

/* the copy's seams made features (creases), their edges remembered */
static qaws_status mr_seams(qaws_mesh_reference* r)
{
	qaws_hm const* h = &r->mesh->hm;
	unsigned int c, n = 0, nfe = qaws_mesh_3d_get_feature_edge_count(r->mesh);
	unsigned int* fe;
	qaws_status st;
	for (c = 0; c < h->nc; c++)
		n += h->cf[c] != HM_NONE && h->co[c] != HM_NONE && c < h->co[c] && !(h->cflag[c] & HM_FEATURE) && mr_seam(h, c, h->co[c]);
	if (!n)
		return QAWS_STATUS_OK;
	fe = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 2, (size_t)nfe + n + 1));
	r->seamk = (unsigned long long*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned long long), (size_t)n + 1));
	if (!fe || !r->seamk)
	{
		qaws_internal_free(fe);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	qaws_mesh_3d_get_feature_edges(fe, r->mesh);
	for (c = 0; c < h->nc; c++)
	{
		unsigned int o = h->co[c];
		unsigned long long a, b;
		if (h->cf[c] == HM_NONE || o == HM_NONE || o < c || (h->cflag[c] & HM_FEATURE) || !mr_seam(h, c, o))
			continue;
		fe[2 * (size_t)nfe] = h->cv[c];
		fe[2 * (size_t)nfe + 1] = h->cv[o];
		nfe++;
		a = h->cv[c];
		b = h->cv[o];
		r->seamk[r->nseam++] = a < b ? (a << 32 | b) : (b << 32 | a);
	}
	qsort(r->seamk, r->nseam, sizeof(unsigned long long), mr_key_cmp);
	st = qaws_mesh_3d_set_feature_edges(r->mesh, fe, nfe);
	qaws_internal_free(fe);
	return st;
}

qaws_status qaws_mesh_reference_create(qaws_mesh_reference** out, qaws_mesh_3d const* mesh, qaws_mesh_reference_desc const* desc)
{
	qaws_mesh_reference* r;
	qaws_hm* m;
	qaws_interpolant kind = desc ? desc->interpolant : QAWS_INTERPOLANT_PN;
	unsigned int f, t, i, k, nt = 0, *ce2t = NULL;
	double *cn = NULL, *fnrm = NULL, *tan = NULL, lo[3] = { 1e300, 1e300, 1e300 }, hi[3] = { -1e300, -1e300, -1e300 };
	unsigned char* has = NULL;
	qaws_status st;
	if (!out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	if (!mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if ((unsigned int)kind > QAWS_INTERPOLANT_MLS)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (kind != QAWS_INTERPOLANT_PN && kind != QAWS_INTERPOLANT_LINEAR && kind != QAWS_INTERPOLANT_PHONG && kind != QAWS_INTERPOLANT_NAGATA && kind != QAWS_INTERPOLANT_G1 &&
		kind != QAWS_INTERPOLANT_SUBDIVISION && kind != QAWS_INTERPOLANT_MLS)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	r = (qaws_mesh_reference*)qaws_internal_zalloc(1, sizeof *r);
	if (!r)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	r->kind = kind;
	st = qaws_mesh_3d_clone(&r->mesh, mesh);
	if (st != QAWS_STATUS_OK)
	{
		qaws_internal_free(r);
		return st;
	}
	if (desc && desc->seams && (st = mr_seams(r)) != QAWS_STATUS_OK)
	{
		qaws_mesh_reference_destroy(r);
		return st;
	}
	m = &r->mesh->hm;
	for (f = 0; f < m->nf; f++)
		nt += m->fn[f] - 2;
	r->nt = nt;
	r->tv = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (nt + 1)));
	r->tc = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (nt + 1)));
	r->tf = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (nt + 1)));
	r->ta = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (nt + 1)));
	r->tcrease = (unsigned char*)qaws_internal_zalloc(3 * (size_t)nt + 1, 1);
	r->treg = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (nt + 1)));
	r->patch = (mi_tri*)qaws_internal_malloc(qaws_internal_mul2(sizeof(mi_tri), (nt + 1)));
	r->box = (qaws_bp_box*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_bp_box), (nt + 1)));
	r->vloc = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (m->nv + 1)));
	ce2t = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (m->nc + 1)));
	cn = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 3, (m->nc + 1)));
	fnrm = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 3, (m->nf + 1)));
	tan = (double*)qaws_internal_scratch_zalloc(3 * (size_t)m->nv + 1, sizeof(double));
	has = (unsigned char*)qaws_internal_scratch_zalloc((size_t)m->nv + 1, 1);
	if (!r->tv || !r->tc || !r->tf || !r->ta || !r->tcrease || !r->treg || !r->patch || !r->box || !r->vloc || !ce2t || !cn || !fnrm || !tan || !has)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	/* triangles: face f fanned from its first corner; each face edge's triangle */
	for (f = 0, t = 0; f < m->nf; f++)
	{
		unsigned int c0 = m->fc[f], s = m->fn[f], j;
		for (j = 0; j + 2 < s; j++, t++)
		{
			r->tc[3 * t] = c0;
			r->tc[3 * t + 1] = c0 + j + 1;
			r->tc[3 * t + 2] = c0 + j + 2;
			for (k = 0; k < 3; k++)
				r->tv[3 * t + k] = m->cv[r->tc[3 * t + k]];
			r->tf[t] = f;
			/* fan diagonals: edge 2 of triangle j against edge 0 of j + 1 */
			r->ta[3 * t] = j > 0 ? 3 * (t - 1) + 2 : HM_NONE;
			r->ta[3 * t + 1] = HM_NONE;
			r->ta[3 * t + 2] = j + 3 < s ? 3 * (t + 1) : HM_NONE;
		}
		for (j = 0; j < s; j++)
			ce2t[c0 + j] = j == 0 ? 3 * (t - (s - 2)) : j <= s - 2 ? 3 * (t - (s - 2) + j - 1) + 1 : 3 * (t - 1) + 2;
	}
	/* across the faces' edges */
	for (f = 0; f < m->nf; f++)
		for (i = m->fc[f]; i < m->fc[f] + m->fn[f]; i++)
		{
			unsigned int e = ce2t[i];
			r->ta[e] = m->co[i] == HM_NONE ? HM_NONE : ce2t[m->co[i]];
			r->tcrease[e] = (m->cflag[i] & HM_FEATURE) != 0 || m->co[i] == HM_NONE;   /* (the border: a crease) */
		}
	/* normals per corner slot */
	if (desc && desc->normals)
	{
		for (i = 0; i < m->nc; i++)
			if (m->cf[i] != HM_NONE)
			{
				for (k = 0; k < 3; k++)
					cn[3 * (size_t)i + k] = (double)desc->normals[3 * (size_t)m->cv[i] + k];
				mr_unit(cn + 3 * (size_t)i);
			}
	}
	else
		mr_corner_normals(m, cn, fnrm);
	mr_crease_tangents(m, tan, has);
	if (kind == QAWS_INTERPOLANT_SUBDIVISION && (st = mr_loop_rings(r)) != QAWS_STATUS_OK)
		goto done;
	if (kind == QAWS_INTERPOLANT_MLS && (st = mr_mls_setup(r, cn, desc)) != QAWS_STATUS_OK)
		goto done;
	/* patches and their boxes */
	for (t = 0; t < nt; t++)
	{
		double const *P[3], *N[3], *crease[3];
		double T[3][6];
		unsigned int e, q;
		for (k = 0; k < 3; k++)
		{
			P[k] = m->x + 3 * (size_t)r->tv[3 * t + k];
			N[k] = cn + 3 * (size_t)r->tc[3 * t + k];
		}
		for (e = 0; e < 3; e++)
		{
			crease[e] = NULL;
			if (!r->tcrease[3 * t + e])
				continue;
			/* each end: its crease line's tangent, or the chord at a corner or a line's end */
			for (q = 0; q < 2; q++)
			{
				unsigned int v = r->tv[3 * t + (e + q) % 3], w = r->tv[3 * t + (e + 1 - q) % 3];
				if (has[v])
					memcpy(T[e] + 3 * q, tan + 3 * (size_t)v, sizeof(double) * 3);
				else
				{
					for (k = 0; k < 3; k++)
						T[e][3 * q + k] = m->x[3 * (size_t)w + k] - m->x[3 * (size_t)v + k];
					mr_unit(T[e] + 3 * q);
				}
			}
			crease[e] = T[e];
		}
		r->patch[t].loop = r->loop ? &r->loop[t] : NULL;
		r->patch[t].mls = r->mls ? &r->mls[t] : NULL;
		mi_build(&r->patch[t], kind, P, N, crease);
		for (k = 0; k < 3; k++)
		{
			r->box[t].lo[k] = 1e300;
			r->box[t].hi[k] = -1e300;
		}
		/* (the Loop limit: within its neighbourhood's hull) */
		for (q = 0; q < (r->loop ? r->loop[t].nv : r->patch[t].gregory ? 16u : 10u); q++)
			for (k = 0; k < 3; k++)
			{
				double v = r->loop ? m->x[3 * (size_t)r->loop[t].v[q] + k] : q < 10 ? r->patch[t].p[q][k] : r->patch[t].n[q - 10][k];
				if (v < r->box[t].lo[k]) r->box[t].lo[k] = v;
				if (v > r->box[t].hi[k]) r->box[t].hi[k] = v;
				if (v < lo[k]) lo[k] = v;
				if (v > hi[k]) hi[k] = v;
			}
		if (r->mls)
		{
			/* (MLS: points sampled 4 x 4 over the triangle, widened by a tenth of its longest side) */
			unsigned int ia, ib;
			double L = 0;
			for (ia = 0; ia <= 4; ia++)
				for (ib = 0; ia + ib <= 4; ib++)
				{
					double y[3];
					mi_eval(&r->patch[t], ia / 4.0, ib / 4.0, y, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
					for (k = 0; k < 3; k++)
					{
						r->box[t].lo[k] = fmin(r->box[t].lo[k], y[k]);
						r->box[t].hi[k] = fmax(r->box[t].hi[k], y[k]);
					}
				}
			for (k = 0; k < 3; k++)
			{
				double const *A = m->x + 3 * (size_t)r->tv[3 * t + k], *B = m->x + 3 * (size_t)r->tv[3 * t + (k + 1) % 3];
				L = fmax(L, sqrt((A[0] - B[0]) * (A[0] - B[0]) + (A[1] - B[1]) * (A[1] - B[1]) + (A[2] - B[2]) * (A[2] - B[2])));
			}
			for (k = 0; k < 3; k++)
			{
				r->box[t].lo[k] -= 0.1 * L;
				r->box[t].hi[k] += 0.1 * L;
				lo[k] = fmin(lo[k], r->box[t].lo[k]);
				hi[k] = fmax(hi[k], r->box[t].hi[k]);
			}
		}
	}
	r->extent = nt ? sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) + (hi[2] - lo[2]) * (hi[2] - lo[2])) : 1;
	for (i = 0; i < m->nv; i++)
		r->vloc[i] = HM_NONE;
	for (t = nt; t > 0; t--)
		for (k = 0; k < 3; k++)
			r->vloc[r->tv[3 * (t - 1) + k]] = 3 * (t - 1) + k;
	/* regions: triangles connected across edges that are not creases */
	for (t = 0; t < nt; t++)
		r->treg[t] = HM_NONE;
	for (t = 0; t < nt; t++)
	{
		unsigned int head = 0, tail = 0, *queue = ce2t;   /* (ce2t is free now: nc >= nt) */
		if (r->treg[t] != HM_NONE)
			continue;
		r->treg[t] = r->nreg;
		queue[tail++] = t;
		while (head < tail)
		{
			unsigned int u = queue[head++], e;
			for (e = 0; e < 3; e++)
			{
				unsigned int nb = r->ta[3 * u + e];
				if (nb == HM_NONE || r->tcrease[3 * u + e] || r->treg[nb / 3] != HM_NONE)
					continue;
				r->treg[nb / 3] = r->nreg;
				queue[tail++] = nb / 3;
			}
		}
		r->nreg++;
	}
	st = nt ? qaws_internal_grid_create(&r->grid, r->box, nt, 3) : QAWS_STATUS_OK;
done:
	qaws_internal_free(ce2t); qaws_internal_free(cn); qaws_internal_free(fnrm); qaws_internal_free(tan); qaws_internal_free(has);
	if (st != QAWS_STATUS_OK)
	{
		qaws_mesh_reference_destroy(r);
		return st;
	}
	*out = r;
	return QAWS_STATUS_OK;
}

unsigned int qaws_mesh_reference_get_triangle_count(qaws_mesh_reference const* r)
{
	return r ? r->nt : 0;
}

qaws_status qaws_mesh_reference_get_triangle(unsigned int out[3], qaws_mesh_reference const* r, unsigned int t)
{
	if (!r || t >= r->nt || !out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	out[0] = r->tv[3 * t];
	out[1] = r->tv[3 * t + 1];
	out[2] = r->tv[3 * t + 2];
	return QAWS_STATUS_OK;
}

qaws_mesh_3d const* qaws_mesh_reference_get_mesh(qaws_mesh_reference const* r)
{
	return r ? r->mesh : NULL;
}

qaws_status qaws_mesh_reference_get_vertex_location(qaws_mesh_location* out, qaws_mesh_reference const* r, unsigned int v)
{
	unsigned int l;
	if (!r || !out || v >= r->mesh->hm.nv || r->vloc[v] == HM_NONE)
		return QAWS_STATUS_INVALID_ARGUMENT;
	l = r->vloc[v];
	out->patch = l / 3;
	out->a = (qaws_scalar)(l % 3 == 1);
	out->b = (qaws_scalar)(l % 3 == 2);
	return QAWS_STATUS_OK;
}

/* ======================================================================== */
/*  Evaluation                                                              */
/* ======================================================================== */

qaws_status qaws_mesh_reference_evaluate(qaws_scalar* out_x, qaws_scalar* out_n, qaws_scalar* out_k,
	qaws_mesh_reference const* r, qaws_mesh_location const* loc, unsigned int count)
{
	unsigned int i, k;
	if (!r || (count && !loc))
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < count; i++)
		if (loc[i].patch >= r->nt)
			return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < count; i++)
	{
		double x[3], xa[3], xb[3], n[3], na[3], nb[3], k1, k2, d1[3], d2[3];
		mi_eval(&r->patch[loc[i].patch], (double)loc[i].a, (double)loc[i].b, x, xa, xb, NULL, NULL, NULL, n, na, nb);
		for (k = 0; k < 3; k++)
		{
			if (out_x) out_x[3 * (size_t)i + k] = (qaws_scalar)x[k];
			if (out_n) out_n[3 * (size_t)i + k] = (qaws_scalar)n[k];
		}
		if (out_k)
		{
			qaws_scalar* o = out_k + 8 * (size_t)i;
			mi_curvature(xa, xb, n, na, nb, &k1, &k2, d1, d2);
			o[0] = (qaws_scalar)k1;
			o[1] = (qaws_scalar)k2;
			for (k = 0; k < 3; k++)
			{
				o[2 + k] = (qaws_scalar)d1[k];
				o[5 + k] = (qaws_scalar)d2[k];
			}
		}
	}
	return QAWS_STATUS_OK;
}

/* ======================================================================== */
/*  Projection                                                              */
/* ======================================================================== */

/* the closest point of triangle (A, B, C) to p, as (a, b) on B and C (Ericson) */
static void mr_closest_flat(double const* A, double const* B, double const* C, double const* p, double* oa, double* ob)
{
	double ab[3], ac[3], ap[3], bp[3], cp[3], d1, d2, d3, d4, d5, d6, va, vb, vc, v, w, den;
	int k;
	for (k = 0; k < 3; k++)
	{
		ab[k] = B[k] - A[k];
		ac[k] = C[k] - A[k];
		ap[k] = p[k] - A[k];
		bp[k] = p[k] - B[k];
		cp[k] = p[k] - C[k];
	}
	d1 = mr_dot(ab, ap); d2 = mr_dot(ac, ap);
	if (d1 <= 0 && d2 <= 0) { *oa = 0; *ob = 0; return; }
	d3 = mr_dot(ab, bp); d4 = mr_dot(ac, bp);
	if (d3 >= 0 && d4 <= d3) { *oa = 1; *ob = 0; return; }
	vc = d1 * d4 - d3 * d2;
	if (vc <= 0 && d1 >= 0 && d3 <= 0) { *oa = d1 / (d1 - d3); *ob = 0; return; }
	d5 = mr_dot(ab, cp); d6 = mr_dot(ac, cp);
	if (d6 >= 0 && d5 <= d6) { *oa = 0; *ob = 1; return; }
	vb = d5 * d2 - d1 * d6;
	if (vb <= 0 && d2 >= 0 && d6 <= 0) { *oa = 0; *ob = d2 / (d2 - d6); return; }
	va = d3 * d6 - d5 * d4;
	if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0)
	{
		w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
		*oa = 1 - w;
		*ob = w;
		return;
	}
	den = 1 / (va + vb + vc);
	v = vb * den;
	w = vc * den;
	*oa = v;
	*ob = w;
}

/* (a, b) into the triangle a >= 0, b >= 0, a + b <= 1 */
static void mr_clamp(double* a, double* b)
{
	if (*a < 0) *a = 0;
	if (*b < 0) *b = 0;
	if (*a + *b > 1)
	{
		double s = *a + *b;
		*a /= s;
		*b /= s;
	}
}

double mr_closest_on_patch(qaws_mesh_reference const* r, unsigned int t, double const* p, double* out_a, double* out_b)
{
	mi_tri const* P = &r->patch[t];
	double a, b, x[3], xa[3], xb[3], d, best;
	int it, k;
	mr_closest_flat(P->p[0], P->p[1], P->p[2], p, &a, &b);
	if (r->kind == QAWS_INTERPOLANT_LINEAR)
	{
		mi_eval(P, a, b, x, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
		*out_a = a;
		*out_b = b;
		return (x[0] - p[0]) * (x[0] - p[0]) + (x[1] - p[1]) * (x[1] - p[1]) + (x[2] - p[2]) * (x[2] - p[2]);
	}
	/* projected Gauss-Newton on |S(a, b) - p|^2, the step halved while it does not improve */
	mi_eval(P, a, b, x, xa, xb, NULL, NULL, NULL, NULL, NULL, NULL);
	best = (x[0] - p[0]) * (x[0] - p[0]) + (x[1] - p[1]) * (x[1] - p[1]) + (x[2] - p[2]) * (x[2] - p[2]);
	for (it = 0; it < 24; it++)
	{
		double rr[3], g0, g1, h00, h01, h11, det, da, db, s, na, nb;
		int improved = 0;
		for (k = 0; k < 3; k++)
			rr[k] = x[k] - p[k];
		g0 = mr_dot(xa, rr); g1 = mr_dot(xb, rr);
		h00 = mr_dot(xa, xa); h01 = mr_dot(xa, xb); h11 = mr_dot(xb, xb);
		det = h00 * h11 - h01 * h01;
		if (det <= 1e-300)
			break;
		da = -(h11 * g0 - h01 * g1) / det;
		db = -(h00 * g1 - h01 * g0) / det;
		for (s = 1; s > 1e-4; s *= 0.5)
		{
			double y[3];
			na = a + s * da;
			nb = b + s * db;
			mr_clamp(&na, &nb);
			mi_eval(P, na, nb, y, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
			d = (y[0] - p[0]) * (y[0] - p[0]) + (y[1] - p[1]) * (y[1] - p[1]) + (y[2] - p[2]) * (y[2] - p[2]);
			if (d < best)
			{
				improved = fabs(na - a) + fabs(nb - b) > 1e-14;
				a = na;
				b = nb;
				best = d;
				break;
			}
		}
		if (!improved)
			break;
		mi_eval(P, a, b, x, xa, xb, NULL, NULL, NULL, NULL, NULL, NULL);
	}
	*out_a = a;
	*out_b = b;
	return best;
}

/* the squared distance from p to a box */
static double mr_box_d2(qaws_bp_box const* b, double const* p)
{
	double s = 0;
	int k;
	for (k = 0; k < 3; k++)
	{
		double d = p[k] < b->lo[k] ? b->lo[k] - p[k] : p[k] > b->hi[k] ? p[k] - b->hi[k] : 0;
		s += d * d;
	}
	return s;
}

typedef struct mr_query
{
	qaws_mesh_reference const* r;
	double p[3];
	double best;                    /* squared */
	qaws_mesh_location loc;
	mr_scratch* s;
	int crease;                     /* only the creases and the border */
	unsigned int region;            /* only this region's triangles; HM_NONE: any */
} mr_query;

/* the closest point of the edge curve e of patch t (from corner e to e + 1) to p: (a, b) and the squared distance */
static double mr_closest_on_edge(qaws_mesh_reference const* r, unsigned int t, unsigned int e, double const* p, double* oa, double* ob)
{
	static double const a0[3] = { 0, 1, 0 }, b0[3] = { 0, 0, 1 }, da[3] = { 1, -1, 0 }, db[3] = { 0, 1, -1 };
	mi_tri const* P = &r->patch[t];
	double s, best = HUGE_VAL, bs = 0;
	int it, k;
	/* start from the chord's closest point */
	{
		double const *A = P->p[e], *B = P->p[(e + 1) % 3];
		double d[3], w = 0, l = 0;
		for (k = 0; k < 3; k++)
		{
			d[k] = B[k] - A[k];
			w += (p[k] - A[k]) * d[k];
			l += d[k] * d[k];
		}
		s = l > 0 ? w / l : 0;
		s = s < 0 ? 0 : s > 1 ? 1 : s;
	}
	for (it = 0; it < 20; it++)
	{
		double x[3], xa[3], xb[3], c1[3], g = 0, h = 0, dist = 0, ns;
		mi_eval(P, a0[e] + s * da[e], b0[e] + s * db[e], x, xa, xb, NULL, NULL, NULL, NULL, NULL, NULL);
		for (k = 0; k < 3; k++)
		{
			c1[k] = da[e] * xa[k] + db[e] * xb[k];
			g += c1[k] * (x[k] - p[k]);
			h += c1[k] * c1[k];
			dist += (x[k] - p[k]) * (x[k] - p[k]);
		}
		if (dist < best)
		{
			best = dist;
			bs = s;
		}
		if (h <= 0)
			break;
		ns = s - g / h;
		ns = ns < 0 ? 0 : ns > 1 ? 1 : ns;
		if (fabs(ns - s) < 1e-14)
			break;
		s = ns;
	}
	*oa = a0[e] + bs * da[e];
	*ob = b0[e] + bs * db[e];
	return best;
}

static void mr_try(mr_query* q, unsigned int t)
{
	double a, b, d;
	unsigned int e;
	if (q->s->tried[t] == q->s->gen)
		return;
	q->s->tried[t] = q->s->gen;
	if (q->region != HM_NONE && q->r->treg[t] != q->region)
		return;
	if (mr_box_d2(&q->r->box[t], q->p) >= q->best)
		return;
	if (!q->crease)
	{
		d = mr_closest_on_patch(q->r, t, q->p, &a, &b);
		if (d < q->best)
		{
			q->best = d;
			q->loc.patch = t;
			q->loc.a = (qaws_scalar)a;
			q->loc.b = (qaws_scalar)b;
		}
		return;
	}
	for (e = 0; e < 3; e++)
		if (q->r->tcrease[3 * t + e] || q->r->ta[3 * t + e] == HM_NONE)
		{
			d = mr_closest_on_edge(q->r, t, e, q->p, &a, &b);
			if (d < q->best)
			{
				q->best = d;
				q->loc.patch = t;
				q->loc.a = (qaws_scalar)a;
				q->loc.b = (qaws_scalar)b;
			}
		}
}

static void mr_visit(void* user, unsigned int item)
{
	mr_try((mr_query*)user, item);
}

int mr_scratch_init(mr_scratch* s, unsigned int nt)
{
	s->tried = (unsigned int*)qaws_internal_zalloc((size_t)nt + 1, sizeof(unsigned int));
	s->queued = (unsigned int*)qaws_internal_zalloc((size_t)nt + 1, sizeof(unsigned int));
	s->queue = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), ((size_t)nt + 1)));
	s->gen = 0;
	if (!s->tried || !s->queued || !s->queue)
	{
		mr_scratch_free(s);
		return 0;
	}
	return 1;
}

void mr_scratch_free(mr_scratch* s)
{
	qaws_internal_free(s->tried);
	qaws_internal_free(s->queued);
	qaws_internal_free(s->queue);
	memset(s, 0, sizeof *s);
}

static double mr_query_run(qaws_mesh_reference const* r, double const* p, unsigned int hint, mr_scratch* s, int crease, unsigned int region,
	qaws_mesh_location* out)
{
	mr_query q;
	q.r = r;
	q.p[0] = p[0]; q.p[1] = p[1]; q.p[2] = p[2];
	q.best = HUGE_VAL;
	q.loc.patch = HM_NONE;
	q.loc.a = q.loc.b = 0;
	q.s = s;
	q.crease = crease;
	q.region = region;
	if (++s->gen == 0)
	{
		/* (the stamps wrapped round: cleared once in 2^32 queries) */
		memset(s->tried, 0, sizeof(unsigned int) * r->nt);
		memset(s->queued, 0, sizeof(unsigned int) * r->nt);
		s->gen = 1;
	}
	if (hint < r->nt && (region == HM_NONE || r->treg[hint] == region))
	{
		/* the triangles reachable from the hint's across edges (within the region), each box within twice the best distance so
		   far (a path to the closest point may pass a little farther than the best found before it) */
		unsigned int head = 0, tail = 0;
		s->queue[tail++] = hint;
		s->queued[hint] = s->gen;
		while (head < tail)
		{
			unsigned int t = s->queue[head++], e;
			mr_try(&q, t);
			for (e = 0; e < 3; e++)
			{
				unsigned int u = r->ta[3 * t + e];
				if (u == HM_NONE || (region != HM_NONE && r->treg[u / 3] != region))
					continue;
				u /= 3;
				if (s->queued[u] != s->gen && mr_box_d2(&r->box[u], q.p) < 4 * q.best)
				{
					s->queued[u] = s->gen;
					s->queue[tail++] = u;
				}
			}
		}
	}
	else
	{
		unsigned int ring;
		for (ring = 0;; ring++)
		{
			double lb = qaws_internal_grid_ring(r->grid, q.p, ring, mr_visit, &q);
			if (lb == HUGE_VAL || lb * lb > q.best)
				break;
		}
	}
	*out = q.loc;
	return q.best;
}

double mr_project(qaws_mesh_reference const* r, double const* p, unsigned int hint, mr_scratch* s, qaws_mesh_location* out)
{
	return mr_query_run(r, p, hint, s, 0, HM_NONE, out);
}

double mr_project_region(qaws_mesh_reference const* r, double const* p, unsigned int hint, unsigned int region, mr_scratch* s,
	qaws_mesh_location* out)
{
	return mr_query_run(r, p, hint, s, 0, region, out);
}

double mr_project_crease(qaws_mesh_reference const* r, double const* p, unsigned int hint, mr_scratch* s, qaws_mesh_location* out)
{
	return mr_query_run(r, p, hint, s, 1, HM_NONE, out);
}

void mr_eval(qaws_mesh_reference const* r, qaws_mesh_location const* l, double* x, double* n)
{
	mi_eval(&r->patch[l->patch], (double)l->a, (double)l->b, x, NULL, NULL, NULL, NULL, NULL, n, NULL, NULL);
}

qaws_status qaws_mesh_reference_project(qaws_mesh_location* out_loc, qaws_scalar* out_dist, qaws_mesh_reference const* r,
	qaws_scalar const* points, qaws_mesh_location const* hints, unsigned int count)
{
	mr_scratch s;
	unsigned int i;
	if (!r || (count && (!points || !out_loc)))
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (!r->nt)
		return count ? QAWS_STATUS_INVALID_ARGUMENT : QAWS_STATUS_OK;
	if (!mr_scratch_init(&s, r->nt))
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < count; i++)
	{
		double p[3], d;
		p[0] = (double)points[3 * (size_t)i];
		p[1] = (double)points[3 * (size_t)i + 1];
		p[2] = (double)points[3 * (size_t)i + 2];
		d = mr_project(r, p, hints && hints[i].patch < r->nt ? hints[i].patch : HM_NONE, &s, &out_loc[i]);
		if (out_dist)
			out_dist[i] = (qaws_scalar)sqrt(d);
	}
	mr_scratch_free(&s);
	return QAWS_STATUS_OK;
}

/* ======================================================================== */
/*  Walking                                                                 */
/* ======================================================================== */

int mr_walk(qaws_mesh_reference const* r, qaws_mesh_location from, double const* tangent, int stop_at_creases, qaws_mesh_location* out_to)
{
	double w[3], dir[3], rem, b[3];
	unsigned int t, k, guard, stalls = 0;
	int stopped = 0;
	t = from.patch;
	b[1] = (double)from.a;
	b[2] = (double)from.b;
	b[0] = 1 - b[1] - b[2];
	for (k = 0; k < 3; k++)
		w[k] = tangent[k];
	rem = sqrt(mr_dot(w, w));
	for (guard = 0; rem > 0 && guard < 100000; guard++)
	{
		double const *P0 = r->patch[t].p[0], *P1 = r->patch[t].p[1], *P2 = r->patch[t].p[2];
		double e1[3], e2[3], nrm[3], g00, g01, g11, det, r0, r1, da, db, d[3], s = HUGE_VAL;
		unsigned int hit = HM_NONE, i;
		for (k = 0; k < 3; k++)
		{
			e1[k] = P1[k] - P0[k];
			e2[k] = P2[k] - P0[k];
		}
		mr_cross(e1, e2, nrm);
		if (mr_unit(nrm) == 0)
			break;
		/* the direction in this triangle's plane */
		{
			double wn = mr_dot(w, nrm);
			for (k = 0; k < 3; k++)
				dir[k] = w[k] - wn * nrm[k];
			if (mr_unit(dir) == 0)
				break;
		}
		g00 = mr_dot(e1, e1); g01 = mr_dot(e1, e2); g11 = mr_dot(e2, e2);
		det = g00 * g11 - g01 * g01;
		r0 = mr_dot(dir, e1); r1 = mr_dot(dir, e2);
		da = (g11 * r0 - g01 * r1) / det;
		db = (g00 * r1 - g01 * r0) / det;
		d[0] = -da - db; d[1] = da; d[2] = db;
		for (i = 0; i < 3; i++)
			if (d[i] < 0 && -b[i] / d[i] < s)
			{
				s = -b[i] / d[i];
				hit = i;
			}
		if (s < 0)
			s = 0;
		/* through a vertex the steps can vanish: give up after a few */
		if (s < 1e-14 * r->extent && ++stalls > 16)
			break;
		if (rem <= s)
		{
			for (i = 0; i < 3; i++)
				b[i] += rem * d[i];
			rem = 0;
			break;
		}
		for (i = 0; i < 3; i++)
			b[i] += s * d[i];
		b[hit] = 0;
		rem -= s;
		{
			/* across the edge opposite corner `hit`: edge (hit + 1) % 3, from corner e0 to e1 */
			unsigned int e = (hit + 1) % 3, c0 = e, c1 = (e + 1) % 3, nb = r->ta[3 * t + e], u, ue, j;
			double const *A = r->patch[t].p[c0], *B = r->patch[t].p[c1];
			double ax[3], mt[3], mu[3], n2[3], lam, de, dm, tw[3];
			if (nb == HM_NONE)
				break;   /* the border */
			if (stop_at_creases && r->tcrease[3 * t + e])
			{
				/* just short of the crease, on this side */
				b[hit] = 1e-6;
				b[c0] -= 0.5e-6;
				b[c1] -= 0.5e-6;
				stopped = 1;
				break;
			}
			u = nb / 3;
			ue = nb % 3;
			for (k = 0; k < 3; k++)
				ax[k] = B[k] - A[k];
			mr_unit(ax);
			/* out of t across the edge, and into u */
			mr_cross(ax, nrm, mt);
			mr_unit(mt);
			{
				double const *Q0 = r->patch[u].p[0], *Q1 = r->patch[u].p[1], *Q2 = r->patch[u].p[2];
				double f1[3], f2[3];
				for (k = 0; k < 3; k++)
				{
					f1[k] = Q1[k] - Q0[k];
					f2[k] = Q2[k] - Q0[k];
				}
				mr_cross(f1, f2, n2);
				mr_unit(n2);
			}
			mr_cross(ax, n2, mu);
			mr_unit(mu);
			de = mr_dot(dir, ax);
			dm = mr_dot(dir, mt);
			for (k = 0; k < 3; k++)
				tw[k] = de * ax[k] + dm * mu[k];
			/* the point on the edge: lam of the way from A to B, then u's barycentrics (its edge ue runs B to A) */
			lam = b[c1];
			memset(b, 0, sizeof b);
			b[ue] = lam;
			b[(ue + 1) % 3] = 1 - lam;
			for (j = 0; j < 3; j++)
				w[j] = tw[j] * rem;
			t = u;
		}
	}
	out_to->patch = t;
	out_to->a = (qaws_scalar)(b[1] < 0 ? 0 : b[1]);
	out_to->b = (qaws_scalar)(b[2] < 0 ? 0 : b[2]);
	return stopped;
}

qaws_status qaws_mesh_reference_walk(qaws_mesh_location* out_to, qaws_mesh_reference const* r, qaws_scalar const* tangent,
	qaws_mesh_location from)
{
	double w[3];
	if (!r || !tangent || !out_to || from.patch >= r->nt)
		return QAWS_STATUS_INVALID_ARGUMENT;
	w[0] = (double)tangent[0];
	w[1] = (double)tangent[1];
	w[2] = (double)tangent[2];
	mr_walk(r, from, w, 0, out_to);
	return QAWS_STATUS_OK;
}

/* ======================================================================== */
/*  Attributes                                                              */
/* ======================================================================== */

qaws_status qaws_mesh_reference_evaluate_attribute(qaws_scalar* out, qaws_mesh_reference const* r, qaws_mesh_location const* loc,
	unsigned int count, unsigned int attribute)
{
	qaws_hm const* m;
	hm_attr const* at;
	unsigned int i, k;
	if (!r || attribute >= r->mesh->hm.na || (count && (!loc || !out)))
		return QAWS_STATUS_INVALID_ARGUMENT;
	m = &r->mesh->hm;
	at = &m->attr[attribute];
	if (at->element == QAWS_MESH_EDGE)
		return QAWS_STATUS_INVALID_ARGUMENT;   /* (an edge attribute has no value inside a face) */
	for (i = 0; i < count; i++)
		if (loc[i].patch >= r->nt)
			return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < count; i++)
	{
		unsigned int t = loc[i].patch, j;
		double w[3], v[16] = { 0 }, best = -1, r0[16];
		w[1] = (double)loc[i].a;
		w[2] = (double)loc[i].b;
		w[0] = 1 - w[1] - w[2];
		if (at->element == QAWS_MESH_FACE)
			memcpy(v, hm_attr_at(m, attribute, r->tf[t]), sizeof(double) * at->n);
		else
			for (j = 0; j < 3; j++)
			{
				double s[16];
				memcpy(s, hm_attr_at(m, attribute, at->element == QAWS_MESH_VERTEX ? r->tv[3 * t + j] : r->tc[3 * t + j]), sizeof(double) * at->n);
				if (at->rule != QAWS_ATTRIBUTE_NEAREST)
					hm_space_in(s, at->n, at->rule);
				if (j == 0)
					memcpy(r0, s, sizeof(double) * at->n);
				hm_align(s, r0, at->n, at->rule);
				if (at->rule == QAWS_ATTRIBUTE_NEAREST)
				{
					if (w[j] > best)
					{
						best = w[j];
						memcpy(v, s, sizeof(double) * at->n);
					}
				}
				else
					for (k = 0; k < at->n; k++)
						v[k] += w[j] * s[k];
			}
		if (at->element != QAWS_MESH_FACE && at->rule != QAWS_ATTRIBUTE_NEAREST)
			hm_space_out(v, at->n, at->rule);
		for (k = 0; k < at->n; k++)
			out[(size_t)at->n * i + k] = (qaws_scalar)v[k];
	}
	return QAWS_STATUS_OK;
}
