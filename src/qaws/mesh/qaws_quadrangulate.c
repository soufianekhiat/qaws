/*
 * Quadrangulation, method INSTANT (Instant Meshes with QuadriFlow's
 * consistency): the input remeshed densely on the reference; a cross field
 * (twice: free, and drawn to the principal curvature directions; the one
 * with fewer singularities kept) and a lattice position field on it; each
 * dense edge's step between lattice points rounded, then mended until the
 * steps close round every face (a min-cost flow, qaws_quad_flow.h: the
 * cones on lattice points, folds undone, positions refitted); the quad mesh
 * read off the steps (lattice points merged, the dense faces' half cells
 * paired across their diagonals into quads, holes filled, doublets
 * removed), for all quads at the target spacing or else split from twice
 * it, relaxed on the reference (creases and the border along themselves);
 * footpoints and attributes as remeshing.
 */

#include "qaws_mesh_internal.h"
#include "../internal/qaws_internal_memory.h"
#include "../qaws_quadrangulate.h"
#include "../qaws_remesh.h"
#include "qaws_cross_field.h"
#include "qaws_quad_flow.h"
#include "qaws_quad_improve.h"
#include "qaws_quad_global.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

static double qd_dot(double const* a, double const* b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static unsigned int qd_find(unsigned int* uf, unsigned int i)
{
	while (uf[i] != i)
	{
		uf[i] = uf[uf[i]];
		i = uf[i];
	}
	return i;
}

/* the faces round which the cross field turns (its singularities) */
static unsigned int qd_singular(qaws_hm const* dm, mf_level const* L)
{
	unsigned int f, i, n = 0;
	for (f = 0; f < dm->nf; f++)
	{
		int t = 0;
		for (i = 0; i < 3; i++)
		{
			unsigned int a = dm->cv[dm->fc[f] + i], b = dm->cv[dm->fc[f] + (i + 1) % 3];
			t += mf_rotation(L->q + 3 * (size_t)a, L->nrm + 3 * (size_t)a, L->q + 3 * (size_t)b, L->nrm + 3 * (size_t)b);
		}
		n += (t & 3) != 0;
	}
	return n;
}

/* Doublets removed (Instant Meshes' clean-up): a vertex off the features in two quads only, which share both its edges
   (a, v, b, x and b, v, a, y), goes and the two become one (a, y, b, x); its two neighbours lose the extra edge. nv: the
   vertices the faces index. */
static qaws_status qd_doublets(unsigned int* fv, unsigned int* fs, unsigned int* nf, unsigned int* nc, unsigned int nv, unsigned char const* kind)
{
	unsigned int *first = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (nv + 2))), *inc = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (*nc + 1)));
	unsigned int* at = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (*nf + 1)));
	unsigned char* mark = (unsigned char*)qaws_internal_scratch_malloc((size_t)*nf + 1);   /* (1: changed this pass, 2: joined away) */
	unsigned int f, k, v, changed;
	if (!first || !inc || !at || !mark)
	{
		qaws_internal_free(first); qaws_internal_free(inc); qaws_internal_free(at); qaws_internal_free(mark);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	do
	{
		unsigned int c0, w, wc;
		changed = 0;
		memset(first, 0, sizeof(unsigned int) * (nv + 2));
		memset(mark, 0, *nf + 1);
		/* each vertex's faces: inc[first[v] .. first[v + 1]) */
		for (f = 0, c0 = 0; f < *nf; c0 += fs[f], f++)
		{
			at[f] = c0;
			for (k = 0; k < fs[f]; k++)
				first[fv[c0 + k] + 2]++;
		}
		for (v = 0; v < nv; v++)
			first[v + 2] += first[v + 1];
		for (f = 0; f < *nf; f++)
			for (k = 0; k < fs[f]; k++)
				inc[first[fv[at[f] + k] + 1]++] = f;
		for (v = 0; v < nv; v++)
		{
			unsigned int f1, f2, i, j, *p, *q, a, b, x, y;
			if (first[v + 1] - first[v] != 2 || kind[v])
				continue;
			f1 = inc[first[v]];
			f2 = inc[first[v] + 1];
			if (f1 == f2 || fs[f1] != 4 || fs[f2] != 4 || mark[f1] || mark[f2])
				continue;
			p = fv + at[f1];
			q = fv + at[f2];
			for (i = 0; p[i] != v; i++) {}
			for (j = 0; q[j] != v; j++) {}
			a = p[(i + 3) & 3];
			b = p[(i + 1) & 3];
			x = p[(i + 2) & 3];
			y = q[(j + 2) & 3];
			if (q[(j + 1) & 3] != a || q[(j + 3) & 3] != b || x == y)
				continue;
			p[0] = a;
			p[1] = y;
			p[2] = b;
			p[3] = x;
			mark[f1] = 1;
			mark[f2] = 2;
			changed = 1;
		}
		/* the joined-away faces dropped */
		for (f = 0, w = 0, wc = 0; f < *nf; f++)
			if (mark[f] != 2)
			{
				memmove(fv + wc, fv + at[f], sizeof(unsigned int) * fs[f]);
				fs[w++] = fs[f];
				wc += fs[f];
			}
		*nf = w;
		*nc = wc;
	} while (changed);
	qaws_internal_free(first); qaws_internal_free(inc); qaws_internal_free(at); qaws_internal_free(mark);
	return QAWS_STATUS_OK;
}

/* where an output vertex goes on the reference: the creases for one made of
   feature vertices, else the region of its first member's footpoint */
static void qd_place(qaws_mesh_reference const* r, mr_scratch* s, double const* p, qaws_mesh_location hint, int feature, double* x,
	qaws_mesh_location* l)
{
	double d = HUGE_VAL;
	if (feature)
		d = mr_project_crease(r, p, hint.patch, s, l);
	if (d == HUGE_VAL || l->patch == HM_NONE)
	{
		unsigned int reg = hint.patch < r->nt ? r->treg[hint.patch] : HM_NONE;
		mr_project_region(r, p, hint.patch, reg, s, l);
	}
	if (l->patch == HM_NONE)
		mr_project(r, p, HM_NONE, s, l);
	mr_eval(r, l, x, NULL);
}

void qd_out_free(qd_out* o)
{
	qaws_internal_free(o->x); qaws_internal_free(o->foot); qaws_internal_free(o->kind); qaws_internal_free(o->fv); qaws_internal_free(o->fs);
	memset(o, 0, sizeof *o);
}

/* The mesh read off consistent steps (QuadriFlow's extraction): dense vertices joined by zero steps are one lattice point
   (a cluster, placed at its members' refitted lattice points); a dense face whose steps close and span half a cell (area
   1/2, no step longer than one) is that half cell, and two half cells sharing their diagonal are its quad. The topology
   comes from the steps alone, so a cell sheared or folded in space still comes out whole. Half cells left alone stay
   triangles; faces over a fold (negative area), a step of two or a cone off the lattice leave holes. */
static qaws_status qd_extract(qaws_hm const* dm, mf_hier const* h, int const* T, unsigned char const* K, qaws_mesh_reference const* r,
	mr_scratch* s, qd_out* out)
{
	mf_level const* L = &h->lv[0];
	unsigned int n = dm->nv, i, c, f, k, ncl = 0, nt = 0, mask;
	unsigned int *uf = NULL, *cl = NULL, *member = NULL, *tri = NULL, *hash = NULL, *used = NULL;
	double *cx = NULL, *cw = NULL, *gx = NULL;
	unsigned char *ckind = NULL, *paired = NULL;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	memset(out, 0, sizeof *out);
	uf = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (n + 1)));
	cl = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (n + 1)));
	member = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (n + 1)));
	cx = (double*)qaws_internal_scratch_zalloc(3 * (size_t)n + 3, sizeof(double));
	cw = (double*)qaws_internal_scratch_zalloc((size_t)n + 1, sizeof(double));
	gx = (double*)qaws_internal_scratch_zalloc(4 * (size_t)n + 4, sizeof(double));   /* (per cluster its guide members' lattice points summed, their count) */
	ckind = (unsigned char*)qaws_internal_scratch_zalloc((size_t)n + 1, 1);
	tri = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (4 * (size_t)dm->nf + 4)));
	paired = (unsigned char*)qaws_internal_scratch_zalloc((size_t)dm->nf + 1, 1);
	if (!uf || !cl || !member || !cx || !cw || !gx || !ckind || !tri || !paired)
		goto done;
	/* the clusters */
	for (i = 0; i < n; i++)
		uf[i] = i;
	for (c = 0; c < dm->nc; c++)
		if (dm->cf[c] != HM_NONE && (dm->co[c] == HM_NONE || c < dm->co[c]) && !T[2 * (size_t)c] && !T[2 * (size_t)c + 1])
		{
			unsigned int ra = qd_find(uf, dm->cv[c]), rb = qd_find(uf, hm_to(dm, c));
			if (ra != rb)
				uf[ra] = rb;
		}
	for (i = 0; i < n; i++)
		cl[i] = HM_NONE;
	for (i = 0; i < n; i++)
	{
		unsigned int root = qd_find(uf, i), q;
		if (cl[root] == HM_NONE)
		{
			member[ncl] = i;
			cl[root] = ncl++;
		}
		q = cl[i] = cl[root];
		for (k = 0; k < 3; k++)
		{
			cx[3 * (size_t)q + k] += L->o[3 * (size_t)i + k];
		}
		cw[q] += 1;
		if (L->con[i] & MF_PINNED)
		{
			ckind[q] = 2;
			member[q] = i;
		}
		else if (L->con[i] & MF_GUIDE)
		{
			/* on a guide line (3, below the creases): placed at its guide members' lattice point */
			for (k = 0; k < 3; k++)
				gx[4 * (size_t)q + k] += L->o[3 * (size_t)i + k];
			gx[4 * (size_t)q + 3] += 1;
			if (!ckind[q])
			{
				ckind[q] = 3;
				member[q] = i;
			}
		}
		else if ((L->con[i] & MF_ON_LINE) && (ckind[q] < 1 || ckind[q] == 3))
			ckind[q] = 1;
	}
	for (i = 0; i < ncl; i++)
		for (k = 0; k < 3; k++)
			cx[3 * (size_t)i + k] /= cw[i];
	/* the half cells: per face its three clusters and which of its edges is the cell's diagonal */
	for (f = 0; f < dm->nf; f++)
	{
		unsigned int c0 = dm->fc[f], a, b, d = 3;
		int e[6], j, i0 = qf_image(dm, f, T, K, e);
		if (i0 < 0)
			continue;
		for (j = 0; j < 3; j++)
		{
			if (abs(e[2 * j]) > 1 || abs(e[2 * j + 1]) > 1)
				break;
			if (abs(e[2 * j]) == 1 && abs(e[2 * j + 1]) == 1)
				d = (unsigned int)j;
		}
		if (j < 3 || d == 3 || e[0] * e[3] - e[1] * e[2] != 1)
			continue;
		a = cl[dm->cv[c0 + i0]];
		b = cl[dm->cv[c0 + (i0 + 1) % 3]];
		c = cl[dm->cv[c0 + (i0 + 2) % 3]];
		if (a == b || b == c || a == c)
			continue;
		tri[4 * (size_t)nt] = a;
		tri[4 * (size_t)nt + 1] = b;
		tri[4 * (size_t)nt + 2] = c;
		tri[4 * (size_t)nt + 3] = d;
		nt++;
	}
	/* each diagonal once in a hash (a second half cell on the same side of it overlaps the first: dropped) */
	mask = qaws_internal_grow(1, qaws_internal_mul2(2, (size_t)nt + 1), 1);
	if (!mask || (mask & (mask - 1)))
		goto done;   /* (a power of two: the probe masks with it) */
	hash = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), mask));
	out->fv = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (4 * (size_t)nt + 1)));
	out->fs = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), ((size_t)nt + 1)));
	if (!hash || !out->fv || !out->fs)
		goto done;
	for (i = 0; i < mask; i++)
		hash[i] = HM_NONE;
	mask -= 1;
	for (i = 0; i < nt; i++)
	{
		unsigned int const* t = tri + 4 * (size_t)i;
		unsigned int u = t[t[3]], v = t[(t[3] + 1) % 3], j;
		for (j = (u * 2654435761u ^ v * 40503u) & mask; hash[j] != HM_NONE; j = (j + 1) & mask)
		{
			unsigned int const* o = tri + 4 * (size_t)hash[j];
			if (o[o[3]] == u && o[(o[3] + 1) % 3] == v)
				break;
		}
		if (hash[j] != HM_NONE)
			paired[i] = 2;
		else
			hash[j] = i;
	}
	/* quads: a half cell and the one across its diagonal; then the half cells left, as triangles */
	for (i = 0; i < nt; i++)
	{
		unsigned int const* t = tri + 4 * (size_t)i;
		unsigned int u = t[t[3]], v = t[(t[3] + 1) % 3], w = t[(t[3] + 2) % 3], j;
		if (paired[i])
			continue;
		for (j = (v * 2654435761u ^ u * 40503u) & mask; hash[j] != HM_NONE; j = (j + 1) & mask)
		{
			unsigned int const* o = tri + 4 * (size_t)hash[j];
			if (o[o[3]] == v && o[(o[3] + 1) % 3] == u)
				break;
		}
		if (hash[j] == HM_NONE || paired[hash[j]])
			continue;
		/* (a partner on the same three clusters is a pillow, a fold: both dropped, qd_fill closes the hole if one is left) */
		if (tri[4 * (size_t)hash[j] + (tri[4 * (size_t)hash[j] + 3] + 2) % 3] == w)
		{
			paired[i] = paired[hash[j]] = 2;
			continue;
		}
		{
			unsigned int const* o = tri + 4 * (size_t)hash[j];
			unsigned int* q = out->fv + out->nc;
			q[0] = v;
			q[1] = w;
			q[2] = u;
			q[3] = o[(o[3] + 2) % 3];
			out->fs[out->nf++] = 4;
			out->nc += 4;
			paired[i] = paired[hash[j]] = 1;
		}
	}
	for (i = 0; i < nt; i++)
		if (!paired[i])
		{
			memcpy(out->fv + out->nc, tri + 4 * (size_t)i, sizeof(unsigned int) * 3);
			out->fs[out->nf++] = 3;
			out->nc += 3;
		}
	if (qd_doublets(out->fv, out->fs, &out->nf, &out->nc, ncl, ckind) != QAWS_STATUS_OK)
		goto done;
	/* the vertices the faces use, placed on the reference */
	used = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (ncl + 1)));
	if (!used)
		goto done;
	for (i = 0; i < ncl; i++)
		used[i] = HM_NONE;
	for (i = 0; i < out->nc; i++)
		if (used[out->fv[i]] == HM_NONE)
			used[out->fv[i]] = out->nv++;
	out->x = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (out->nv + 1)));
	out->foot = (qaws_mesh_location*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_mesh_location), (out->nv + 1)));
	out->kind = (unsigned char*)qaws_internal_malloc(out->nv + 1);
	if (!out->x || !out->foot || !out->kind)
		goto done;
	for (i = 0; i < ncl; i++)
	{
		unsigned int q = used[i];
		if (q == HM_NONE)
			continue;
		if (ckind[i] == 2)
		{
			/* a corner: where its pinned vertex is */
			memcpy(out->x + 3 * (size_t)q, dm->x + 3 * (size_t)member[i], sizeof(double) * 3);
			out->foot[q] = dm->foot[member[i]];
		}
		else if (ckind[i] == 3)
		{
			double g[3];
			for (k = 0; k < 3; k++)
				g[k] = gx[4 * (size_t)i + k] / gx[4 * (size_t)i + 3];
			qd_place(r, s, g, dm->foot[member[i]], 0, out->x + 3 * (size_t)q, &out->foot[q]);
		}
		else
			qd_place(r, s, cx + 3 * (size_t)i, dm->foot[member[i]], ckind[i] != 0, out->x + 3 * (size_t)q, &out->foot[q]);
		out->kind[q] = ckind[i];
	}
	for (i = 0; i < out->nc; i++)
		out->fv[i] = used[out->fv[i]];
	st = QAWS_STATUS_OK;
done:
	qaws_internal_free(uf); qaws_internal_free(cl); qaws_internal_free(member); qaws_internal_free(tri); qaws_internal_free(hash); qaws_internal_free(used);
	qaws_internal_free(cx); qaws_internal_free(cw); qaws_internal_free(gx); qaws_internal_free(ckind); qaws_internal_free(paired);
	if (st != QAWS_STATUS_OK)
		qd_out_free(out);
	return st;
}

/* Small holes the extraction left (a loop of at most 8 border edges) closed by a face each; a loop round a lattice point
   with no vertex (8 edges, or up to 10 round a cone on a lattice point off its face's vertices: cone holds those faces'
   centres) by a vertex there: joined to every other loop vertex (6, 8, 10 edges: quads round a vertex of valence 3, 4
   or 5), or to all of them (a fan, its triangles left to pair with the half cells across). */
static qaws_status qd_fill(qd_out* o, qaws_mesh_reference const* r, mr_scratch* s, double const* cone, unsigned int ncone, double rho)
{
	unsigned int f, k, c0, mask, *h, *start, nadd = 0, added = 0, *nf, *nfs, nv0 = o->nv, open = 0;
	qaws_status st = QAWS_STATUS_OK;
	for (k = 0; k < 3 * r->nt && !open; k++)
		open = r->ta[k] == HM_NONE;
	mask = qaws_internal_grow(1, qaws_internal_mul2(2, (size_t)o->nc + 1), 1);
	if (!mask || (mask & (mask - 1)))
		return QAWS_STATUS_ALLOCATION_FAILURE;   /* (a power of two: the probe masks with it) */
	h = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 2, mask));
	start = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (o->nv + 1)));
	if (!h || !start)
	{
		qaws_internal_free(h);
		qaws_internal_free(start);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (k = 0; k < mask; k++)
		h[2 * k] = HM_NONE;
	mask -= 1;
	/* the faces' directed edges, in a hash */
	for (f = 0, c0 = 0; f < o->nf; c0 += o->fs[f], f++)
		for (k = 0; k < o->fs[f]; k++)
		{
			unsigned int a = o->fv[c0 + k], b = o->fv[c0 + (k + 1) % o->fs[f]], j;
			for (j = (a * 2654435761u ^ b * 40503u) & mask; h[2 * j] != HM_NONE; j = (j + 1) & mask) {}
			h[2 * j] = a;
			h[2 * j + 1] = b;
		}
	/* a border edge a -> b (no b -> a): the hole goes b -> a, then on from a; start[a] = b. A vertex on two holes: none */
	for (k = 0; k < o->nv; k++)
		start[k] = HM_NONE;
	for (f = 0, c0 = 0; f < o->nf; c0 += o->fs[f], f++)
		for (k = 0; k < o->fs[f]; k++)
		{
			unsigned int a = o->fv[c0 + k], b = o->fv[c0 + (k + 1) % o->fs[f]], j, rev = 0;
			for (j = (b * 2654435761u ^ a * 40503u) & mask; h[2 * j] != HM_NONE && !rev; j = (j + 1) & mask)
				rev = h[2 * j] == b && h[2 * j + 1] == a;
			if (rev)
				continue;
			start[b] = start[b] == HM_NONE ? a : HM_NONE - 1;   /* (HM_NONE - 1: ambiguous) */
		}
	nf = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (o->nc + 8 * (size_t)o->nv + 1)));
	nfs = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (o->nf + o->nv + 1)));
	{
		/* (room for a centre per loop, the loops of three or more sharing no vertex) */
		size_t nmax = (size_t)o->nv + o->nv / 3 + 1;
		double* x = (double*)qaws_internal_realloc(o->x, qaws_internal_mul3(sizeof(double), 3, nmax));
		qaws_mesh_location* foot;
		unsigned char* kind;
		if (x)
			o->x = x;
		foot = x ? (qaws_mesh_location*)qaws_internal_realloc(o->foot, qaws_internal_mul2(sizeof(qaws_mesh_location), nmax)) : NULL;
		if (foot)
			o->foot = foot;
		kind = foot ? (unsigned char*)qaws_internal_realloc(o->kind, nmax) : NULL;
		if (kind)
			o->kind = kind;
		if (!nf || !nfs || !kind)
		{
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			goto done;
		}
	}
	memcpy(nf, o->fv, sizeof(unsigned int) * o->nc);
	memcpy(nfs, o->fs, sizeof(unsigned int) * o->nf);
	for (k = 0; k < nv0; k++)
	{
		unsigned int loop[11], len = 0, v = k, ok = 0, j, q;
		if (start[k] >= HM_NONE - 1)
			continue;
		while (len < 11)
		{
			unsigned int nx = start[v];
			if (nx >= HM_NONE - 1)
				break;
			loop[len++] = v;
			v = nx;
			if (v == k)
			{
				ok = 1;
				break;
			}
		}
		for (j = 0, q = 0; j < len; j++)
		{
			start[loop[j]] = HM_NONE;   /* (each hole once) */
			q += o->kind[loop[j]] != 0;
		}
		if (!ok || len < 3 || len > 10 || (open && q == len))
			continue;   /* (on an open reference a loop all on its lines is taken for its border, left open) */
		if (len == 8 || ncone)
		{
			double m[3] = { 0, 0, 0 }, near = HUGE_VAL, sum[2] = { 0, 0 };
			for (j = 0; j < len; j++)
				for (q = 0; q < 3; q++)
					m[q] += o->x[3 * (size_t)loop[j] + q] / len;
			for (j = 0; j < ncone && len != 8; j++)
			{
				double d2 = 0;
				for (q = 0; q < 3; q++)
					d2 += (cone[3 * (size_t)j + q] - m[q]) * (cone[3 * (size_t)j + q] - m[q]);
				near = d2 < near ? d2 : near;
			}
			if (len == 8 || near < rho * rho)
			{
				unsigned int c = o->nv++, off;
				/* the spokes to the loop's nearer half (the lattice neighbours, not the cells' far corners) */
				for (j = 0; j < len; j++)
				{
					double d2 = 0;
					for (q = 0; q < 3; q++)
						d2 += (o->x[3 * (size_t)loop[j] + q] - m[q]) * (o->x[3 * (size_t)loop[j] + q] - m[q]);
					sum[j & 1] += sqrt(d2);
				}
				off = sum[1] < sum[0];
				qd_place(r, s, m, o->foot[loop[0]], 0, o->x + 3 * (size_t)c, &o->foot[c]);
				o->kind[c] = 0;
				for (j = 0; (len & 1) || len == 4 ? j < len : 0; j++)
				{
					/* (a loop of the cells' diagonals, or of some: a fan, its triangles then paired with the half cells across) */
					unsigned int* t = nf + o->nc + added;
					t[0] = loop[j];
					t[1] = loop[(j + 1) % len];
					t[2] = c;
					nfs[o->nf + nadd++] = 3;
					added += 3;
				}
				for (j = 0; (len & 1) || len == 4 ? 0 : j < len / 2; j++)
				{
					unsigned int* t = nf + o->nc + added;
					t[0] = loop[(off + 2 * j) % len];
					t[1] = loop[(off + 2 * j + 1) % len];
					t[2] = loop[(off + 2 * j + 2) % len];
					t[3] = c;
					nfs[o->nf + nadd++] = 4;
					added += 4;
				}
				continue;
			}
		}
		if (len > 8)
			continue;
		memcpy(nf + o->nc + added, loop, sizeof(unsigned int) * len);
		nfs[o->nf + nadd++] = len;
		added += len;
	}
	qaws_internal_free(o->fv);
	qaws_internal_free(o->fs);
	o->fv = nf;
	o->fs = nfs;
	o->nf += nadd;
	o->nc += added;
	nf = nfs = NULL;
done:
	qaws_internal_free(nf);
	qaws_internal_free(nfs);
	qaws_internal_free(h);
	qaws_internal_free(start);
	return st;
}

/* every face split into quads round its centre through its edges' midpoints */
static qaws_status qd_split(qd_out* o, qaws_mesh_reference const* r, mr_scratch* s)
{
	unsigned int f, k, q, nv2, *mid, *keys, c0, nf2 = 0, mask;
	qd_out n;
	memset(&n, 0, sizeof n);
	mask = qaws_internal_grow(1, qaws_internal_mul2(2, (size_t)o->nc + 1), 1);
	if (!mask || (mask & (mask - 1)))
		return QAWS_STATUS_ALLOCATION_FAILURE;   /* (a power of two: the probe masks with it) */
	mid = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (o->nc + 1)));
	keys = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, mask));   /* (a hash: low end, high end, midpoint) */
	n.x = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)o->nv + o->nc + o->nf + 1));
	n.foot = (qaws_mesh_location*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_mesh_location), (size_t)o->nv + o->nc + o->nf + 1));
	n.kind = (unsigned char*)qaws_internal_zalloc((size_t)o->nv + o->nc + o->nf + 1, 1);
	n.fv = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 4, (o->nc + 1)));
	n.fs = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (o->nc + 1)));
	if (!mid || !keys || !n.x || !n.foot || !n.kind || !n.fv || !n.fs)
	{
		qaws_internal_free(mid);
		qaws_internal_free(keys);
		qd_out_free(&n);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (k = 0; k < mask; k++)
		keys[3 * k] = HM_NONE;
	mask -= 1;
	memcpy(n.x, o->x, sizeof(double) * 3 * o->nv);
	memcpy(n.foot, o->foot, sizeof(qaws_mesh_location) * o->nv);
	memcpy(n.kind, o->kind, o->nv);
	nv2 = o->nv;
	for (f = 0, c0 = 0; f < o->nf; c0 += o->fs[f], f++)
	{
		double p[3] = { 0, 0, 0 };
		unsigned int c, reg = HM_NONE, hint = HM_NONE;
		qaws_mesh_location l;
		/* the face's region: from a vertex off the features, else under its centre */
		for (k = 0; k < o->fs[f]; k++)
			for (q = 0; q < 3; q++)
				p[q] += o->x[3 * (size_t)o->fv[c0 + k] + q] / o->fs[f];
		for (k = 0; k < o->fs[f] && hint == HM_NONE; k++)
			if (!o->kind[o->fv[c0 + k]] && o->foot[o->fv[c0 + k]].patch < r->nt)
				hint = o->foot[o->fv[c0 + k]].patch;
		if (hint == HM_NONE)
		{
			mr_project(r, p, HM_NONE, s, &l);
			hint = l.patch;
		}
		reg = hint < r->nt ? r->treg[hint] : HM_NONE;
		/* its edges' midpoints, made once: on the crease for an edge along one, in the region otherwise */
		for (k = 0; k < o->fs[f]; k++)
		{
			unsigned int a = o->fv[c0 + k], b = o->fv[c0 + (k + 1) % o->fs[f]], j, found = HM_NONE;
			unsigned int lo = a < b ? a : b, hi = a < b ? b : a;
			for (j = (lo * 2654435761u ^ hi * 40503u) & mask; keys[3 * j] != HM_NONE; j = (j + 1) & mask)
				if (keys[3 * j] == lo && keys[3 * j + 1] == hi)
				{
					found = keys[3 * j + 2];
					break;
				}
			if (found == HM_NONE)
			{
				double m2[3], d2 = 0;
				int along = 0;
				for (q = 0; q < 3; q++)
				{
					m2[q] = 0.5 * (o->x[3 * (size_t)a + q] + o->x[3 * (size_t)b + q]);
					d2 += (o->x[3 * (size_t)a + q] - o->x[3 * (size_t)b + q]) * (o->x[3 * (size_t)a + q] - o->x[3 * (size_t)b + q]);
				}
				if (o->kind[a] && o->kind[b] && o->kind[a] != 3 && o->kind[b] != 3)
					along = mr_project_crease(r, m2, o->foot[a].patch, s, &l) < 0.01 * d2;
				if (along)
					mr_eval(r, &l, n.x + 3 * (size_t)nv2, NULL);
				else
				{
					mr_project_region(r, m2, hint, reg, s, &l);
					if (l.patch == HM_NONE)
						mr_project(r, m2, HM_NONE, s, &l);
					mr_eval(r, &l, n.x + 3 * (size_t)nv2, NULL);
				}
				n.foot[nv2] = l;
				n.kind[nv2] = (unsigned char)(o->kind[a] == 3 && o->kind[b] == 3 ? 3 : along);   /* (along a guide line: on it) */
				keys[3 * j] = lo; keys[3 * j + 1] = hi; keys[3 * j + 2] = nv2;
				found = nv2++;
			}
			mid[c0 + k] = found;
		}
		/* its centre, in the region, and its quads */
		mr_project_region(r, p, hint, reg, s, &l);
		if (l.patch == HM_NONE)
			mr_project(r, p, HM_NONE, s, &l);
		mr_eval(r, &l, n.x + 3 * (size_t)nv2, NULL);
		n.foot[nv2] = l;
		c = nv2++;
		for (k = 0; k < o->fs[f]; k++)
		{
			unsigned int* t = n.fv + 4 * (size_t)nf2;
			t[0] = o->fv[c0 + k];
			t[1] = mid[c0 + k];
			t[2] = c;
			t[3] = mid[c0 + (k + o->fs[f] - 1) % o->fs[f]];
			n.fs[nf2++] = 4;
		}
	}
	n.nv = nv2;
	n.nf = nf2;
	n.nc = 4 * nf2;
	qaws_internal_free(mid);
	qaws_internal_free(keys);
	qd_out_free(o);
	*o = n;
	return QAWS_STATUS_OK;
}

/* triangles that share a side merged into a quad each (a half cell and its neighbour across a fold, a hole's filling).
   shortcut: pairs found by a scan over all triangle pairs, fine for the few left; a hash of their sides if many are */
static qaws_status qd_pair_triangles(qd_out* o)
{
	unsigned int f, g, k, j, c0, nf = 0, nc = 0;
	unsigned int *at = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), ((size_t)o->nf + 1))), *mate = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), ((size_t)o->nf + 1)));
	unsigned int *nfv = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), ((size_t)o->nc + 1))), *nfs = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), ((size_t)o->nf + 1)));
	unsigned char* side = (unsigned char*)qaws_internal_scratch_malloc((size_t)o->nf + 1);
	if (!at || !mate || !nfv || !nfs || !side)
	{
		qaws_internal_free(at); qaws_internal_free(mate); qaws_internal_free(nfv); qaws_internal_free(nfs); qaws_internal_free(side);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (f = 0, c0 = 0; f < o->nf; c0 += o->fs[f], f++)
	{
		at[f] = c0;
		mate[f] = HM_NONE;
	}
	for (f = 0; f < o->nf; f++)
		for (k = 0; k < 3 && o->fs[f] == 3 && mate[f] == HM_NONE; k++)
			for (g = f + 1; g < o->nf && mate[f] == HM_NONE; g++)
			{
				unsigned int const *t = o->fv + at[f], *u = o->fv + at[g];
				if (o->fs[g] != 3 || mate[g] != HM_NONE)
					continue;
				for (j = 0; j < 3; j++)
					if (u[j] == t[(k + 1) % 3] && u[(j + 1) % 3] == t[k] && u[(j + 2) % 3] != t[(k + 2) % 3])
					{
						mate[f] = g;
						mate[g] = f;
						side[f] = (unsigned char)k;
						side[g] = (unsigned char)j;
						break;
					}
			}
	for (f = 0; f < o->nf; f++)
	{
		unsigned int const* t = o->fv + at[f];
		if (mate[f] == HM_NONE)
		{
			memcpy(nfv + nc, t, sizeof(unsigned int) * o->fs[f]);
			nfs[nf++] = o->fs[f];
			nc += o->fs[f];
		}
		else if (f < mate[f])
		{
			unsigned int const* u = o->fv + at[mate[f]];
			k = side[f];
			j = side[mate[f]];
			nfv[nc] = t[(k + 1) % 3];
			nfv[nc + 1] = t[(k + 2) % 3];
			nfv[nc + 2] = t[k];
			nfv[nc + 3] = u[(j + 2) % 3];
			nfs[nf++] = 4;
			nc += 4;
		}
	}
	qaws_internal_free(o->fv);
	qaws_internal_free(o->fs);
	o->fv = nfv;
	o->fs = nfs;
	o->nf = nf;
	o->nc = nc;
	qaws_internal_free(at); qaws_internal_free(mate); qaws_internal_free(side);
	return QAWS_STATUS_OK;
}

/* the quad of triangles t (side k shared) and u (side j): its corners from t's far end of the side round to u's apex */
static void qd_quad_of(unsigned int const* t, unsigned int k, unsigned int const* u, unsigned int j, unsigned int* q)
{
	q[0] = t[(k + 1) % 3];
	q[1] = t[(k + 2) % 3];
	q[2] = t[k];
	q[3] = u[(j + 2) % 3];
}

/* a quad's worth: the smallest sine of its corners (0 when it is not convex in its mean plane), times the cosine of the
   fold between its triangles (none past 40 degrees) */
static double qd_quad_worth(double const* x, unsigned int const* q)
{
	double n[3] = { 0, 0, 0 }, w = HUGE_VAL, nl, f0[3], f1[3], c0, l0, l1;
	unsigned int i, k;
	for (i = 0; i < 4; i++)
	{
		double const *a = x + 3 * (size_t)q[i], *b = x + 3 * (size_t)q[(i + 1) % 4];
		n[0] += (a[1] - b[1]) * (a[2] + b[2]);
		n[1] += (a[2] - b[2]) * (a[0] + b[0]);
		n[2] += (a[0] - b[0]) * (a[1] + b[1]);
	}
	nl = sqrt(qd_dot(n, n));
	if (!(nl > 0))
		return 0;
	for (i = 0; i < 4; i++)
	{
		double const *p0 = x + 3 * (size_t)q[i], *p1 = x + 3 * (size_t)q[(i + 1) % 4], *p3 = x + 3 * (size_t)q[(i + 3) % 4];
		double e1[3], e2[3], cr[3], s;
		for (k = 0; k < 3; k++) { e1[k] = p1[k] - p0[k]; e2[k] = p3[k] - p0[k]; }
		cr[0] = e1[1] * e2[2] - e1[2] * e2[1];
		cr[1] = e1[2] * e2[0] - e1[0] * e2[2];
		cr[2] = e1[0] * e2[1] - e1[1] * e2[0];
		s = qd_dot(cr, n) / nl;
		l0 = sqrt(qd_dot(e1, e1) * qd_dot(e2, e2));
		if (!(s > 0) || !(l0 > 0))
			return 0;
		w = fmin(w, s / l0);
	}
	/* the fold: the normals of (q0 q1 q2) and (q0 q2 q3) */
	{
		double const *a = x + 3 * (size_t)q[0], *b = x + 3 * (size_t)q[1], *c = x + 3 * (size_t)q[2], *d = x + 3 * (size_t)q[3];
		double u[3], v[3], ww[3];
		for (k = 0; k < 3; k++) { u[k] = b[k] - a[k]; v[k] = c[k] - a[k]; ww[k] = d[k] - a[k]; }
		f0[0] = u[1] * v[2] - u[2] * v[1]; f0[1] = u[2] * v[0] - u[0] * v[2]; f0[2] = u[0] * v[1] - u[1] * v[0];
		f1[0] = v[1] * ww[2] - v[2] * ww[1]; f1[1] = v[2] * ww[0] - v[0] * ww[2]; f1[2] = v[0] * ww[1] - v[1] * ww[0];
		l0 = sqrt(qd_dot(f0, f0));
		l1 = sqrt(qd_dot(f1, f1));
		c0 = l0 > 0 && l1 > 0 ? qd_dot(f0, f1) / (l0 * l1) : -1;
	}
	return c0 < 0.766 ? 0 : w * c0;
}

/* PAIRING and SPLIT: the input's own faces (no dense remesh, no fields). PAIRING pairs its triangles across their shared
   sides (off feature lines) into quads: a greedy matching by the quads' worth, then augmenting paths (an unpaired triangle
   takes a paired neighbour, whose old partner takes the next, until one takes an unpaired triangle); with pure,
   the result split round each face's centre (all quads at half the size). SPLIT splits every face. The vertices are the
   input's, on the reference; their kinds from its features and border. */
static qaws_status qd_direct(qaws_mesh_3d const* mesh, qaws_quadrangulate_desc const* desc, qaws_mesh_reference const* r, mr_scratch* s, qd_out* o)
{
	qaws_hm const* m = &mesh->hm;
	unsigned int v, f, k, nc = 0, *mate = NULL, *side = NULL, *nfv = NULL, *nfs = NULL, nf2 = 0, nc2 = 0;
	double* wt = NULL;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	memset(o, 0, sizeof *o);
	for (f = 0; f < m->nf; f++) nc += m->fn[f];
	o->x = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)m->nv + 1));
	o->foot = (qaws_mesh_location*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_mesh_location), (size_t)m->nv + 1));
	o->kind = (unsigned char*)qaws_internal_zalloc((size_t)m->nv + 1, 1);
	o->fv = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)nc + 1));
	o->fs = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)m->nf + 1));
	if (!o->x || !o->foot || !o->kind || !o->fv || !o->fs)
		goto done;
	o->nv = m->nv;
	memcpy(o->x, m->x, sizeof(double) * 3 * m->nv);
	for (v = 0; v < m->nv; v++)
	{
		mr_project(r, o->x + 3 * (size_t)v, HM_NONE, s, &o->foot[v]);
		o->kind[v] = (unsigned char)((m->vflag[v] & HM_CORNER) ? 2 : (hm_is_border_vertex(m, v) || hm_feature_degree(m, v)) ? 1 : 0);
	}
	for (f = 0; f < m->nf; f++)
	{
		for (k = 0; k < m->fn[f]; k++) o->fv[o->nc + k] = m->cv[m->fc[f] + k];
		o->fs[o->nf++] = m->fn[f];
		o->nc += m->fn[f];
	}
	if (desc->method == QAWS_QUADS_PAIRING)
	{
		/* the candidates: each triangle's sides shared with a triangle, off the features, with a quad worth having */
		unsigned int changed = 1, pass;
		mate = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)m->nf + 1));
		side = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)m->nf + 1));
		wt = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)m->nf + 1));
		nfv = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)nc + 1));
		nfs = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)m->nf + 1));
		if (!mate || !side || !wt || !nfv || !nfs)
			goto done;
		for (f = 0; f < m->nf; f++)
		{
			mate[f] = HM_NONE;
			for (k = 0; k < 3; k++)
			{
				unsigned int c = m->fc[f] + k, oc = m->co[c], g, q[4], j;
				wt[3 * (size_t)f + k] = 0;
				if (m->fn[f] != 3 || oc == HM_NONE || (m->cflag[c] & HM_FEATURE))
					continue;
				g = m->cf[oc];
				if (m->fn[g] != 3)
					continue;
				j = oc - m->fc[g];
				qd_quad_of(m->cv + m->fc[f], k, m->cv + m->fc[g], j, q);
				wt[3 * (size_t)f + k] = qd_quad_worth(m->x, q);
			}
		}
		/* greedy: the worthiest pairs first (a pass per tenth of worth, best sides first within it) */
		for (pass = 10; pass-- > 0;)
			for (f = 0; f < m->nf; f++)
			{
				unsigned int bk = 3;
				double bw = 0.1 * pass;
				if (mate[f] != HM_NONE)
					continue;
				for (k = 0; k < 3; k++)
				{
					unsigned int g;
					if (!(wt[3 * (size_t)f + k] > bw) || wt[3 * (size_t)f + k] < 0.05)
						continue;
					g = m->cf[m->co[m->fc[f] + k]];
					if (mate[g] != HM_NONE)
						continue;
					bw = wt[3 * (size_t)f + k];
					bk = k;
				}
				if (bk < 3)
				{
					unsigned int oc = m->co[m->fc[f] + bk], g = m->cf[oc];
					mate[f] = g;
					mate[g] = f;
					side[f] = bk;
					side[g] = oc - m->fc[g];
				}
			}
		/* augmenting paths (Edmonds' search, blossoms not shrunk): from each unpaired triangle, breadth first along unpaired
		   sides to a neighbour and on along its pair; reaching an unpaired triangle, the pairs along the path are swapped,
		   one triangle more paired */
		{
			unsigned int *par = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)m->nf + 1));
			unsigned int *stamp = (unsigned int*)qaws_internal_scratch_zalloc((size_t)m->nf + 1, sizeof(unsigned int));
			unsigned int *queue = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)m->nf + 1));
			unsigned int round = 0;
			if (!par || !stamp || !queue)
			{
				qaws_internal_free(par);
				qaws_internal_free(stamp);
				qaws_internal_free(queue);
				goto done;
			}
			while (changed)
			{
				changed = 0;
				for (f = 0; f < m->nf; f++)
				{
					unsigned int qh = 0, qt = 0, end = HM_NONE, endk = 0;
					if (m->fn[f] != 3 || mate[f] != HM_NONE)
						continue;
					round++;
					stamp[f] = round;
					queue[qt++] = f;
					while (qh < qt && end == HM_NONE && qt < 4096)
					{
						unsigned int e = queue[qh++], ka;
						for (ka = 0; ka < 3 && end == HM_NONE; ka++)
						{
							unsigned int a2, b2;
							if (wt[3 * (size_t)e + ka] < 0.05)
								continue;
							a2 = m->cf[m->co[m->fc[e] + ka]];
							if (stamp[a2] == round || a2 == mate[e])
								continue;
							stamp[a2] = round;
							par[a2] = e;   /* (a2 reached from e along an unpaired side) */
							if (mate[a2] == HM_NONE)
							{
								end = a2;
								endk = ka;
								break;
							}
							b2 = mate[a2];
							if (stamp[b2] == round)
								continue;
							stamp[b2] = round;
							queue[qt++] = b2;
						}
					}
					if (end == HM_NONE)
						continue;
					/* along the path back: each (par, node) paired, its node's old pair freed toward the start */
					{
						unsigned int node = end, from;
						(void)endk;
						for (;;)
						{
							unsigned int old, kk;
							from = par[node];
							old = mate[from];
							for (kk = 0; kk < 3 && m->cf[m->co[m->fc[from] + kk]] != node; kk++) {}
							mate[from] = node;
							mate[node] = from;
							side[from] = kk;
							side[node] = m->co[m->fc[from] + kk] - m->fc[node];
							if (old == HM_NONE)
								break;
							node = old;
						}
					}
					changed = 1;
				}
			}
			qaws_internal_free(par);
			qaws_internal_free(stamp);
			qaws_internal_free(queue);
		}
		for (f = 0; f < m->nf; f++)
		{
			unsigned int const* t = m->cv + m->fc[f];
			if (mate[f] == HM_NONE)
			{
				memcpy(nfv + nc2, t, sizeof(unsigned int) * m->fn[f]);
				nfs[nf2++] = m->fn[f];
				nc2 += m->fn[f];
			}
			else if (f < mate[f])
			{
				qd_quad_of(t, side[f], m->cv + m->fc[mate[f]], side[mate[f]], nfv + nc2);
				nfs[nf2++] = 4;
				nc2 += 4;
			}
		}
		qaws_internal_free(o->fv);
		qaws_internal_free(o->fs);
		o->fv = nfv;
		o->fs = nfs;
		o->nf = nf2;
		o->nc = nc2;
		nfv = nfs = NULL;
	}
	st = desc->method == QAWS_QUADS_SPLIT || desc->pure ? qd_split(o, r, s) : QAWS_STATUS_OK;
done:
	qaws_internal_free(mate);
	qaws_internal_free(side);
	qaws_internal_free(wt);
	qaws_internal_free(nfv);
	qaws_internal_free(nfs);
	return st;
}

/* a vertex where the quads' grid is not regular: a feature corner, a feature line's or the border's (the separatrices end
   there), or one inside whose valence is not 4 or round which a face is not a quad */
static int qd_layout_stop(qaws_hm const* m, unsigned int v)
{
	unsigned int c;
	if ((m->vflag[v] & HM_CORNER) || hm_is_border_vertex(m, v) || hm_feature_degree(m, v) || hm_valence(m, v) != 4)
		return 1;
	for (c = m->vc[v]; c != HM_NONE; c = hm_round(m, v, c))
		if (m->fn[m->cf[c]] != 4)
			return 1;
	return 0;
}

static void qd_layout_cut(qaws_hm const* m, unsigned char* cut, unsigned int a, unsigned int b)
{
	unsigned int c = hm_find_edge(m, a, b);
	if (c == HM_NONE)
		return;
	cut[c] = 1;
	if (m->co[c] != HM_NONE)
		cut[m->co[c]] = 1;
}

/* a separatrix from v through its neighbour w, straight on through regular vertices (marking them), until a vertex
   already on one or where the grid is not regular */
static void qd_layout_trace(qaws_hm const* m, unsigned char* cut, unsigned char* on, unsigned int v, unsigned int w)
{
	unsigned int prev = v, cur = w, steps;
	qd_layout_cut(m, cut, prev, cur);
	for (steps = 0; steps < m->nv; steps++)
	{
		unsigned int r2[4], j;
		if (on[cur] || qd_layout_stop(m, cur) || hm_ring(m, cur, r2, 4) != 4)
			break;
		on[cur] = 2;
		for (j = 0; j < 4 && r2[j] != prev; j++) {}
		if (j == 4)
			break;
		prev = cur;
		cur = r2[(j + 2) % 4];
		qd_layout_cut(m, cut, prev, cur);
	}
}

/* LAYOUT: the quads parted into patches by their motorcycle graph: from each vertex where the grid is not regular (a
   valence other than 4 inside or on a feature line, 3 on the border; a feature corner) a separatrix along each of its
   edges off the border and the features, straight on through regular vertices (the edge across), until it reaches the
   border, a feature line, a singular vertex or one an earlier separatrix passed. The patches: faces joined across edges
   no separatrix, feature line or border parts, each a grid of quads (a rectangle; or periodic where no singularity
   bounds it: an annulus between two creases, a torus); their indices a face attribute "patch" (NEAREST); out_count, how
   many. */
static qaws_status qd_patches(qaws_mesh_3d* out, unsigned int* out_count)
{
	qaws_hm* m = &out->hm;
	unsigned char *cut = (unsigned char*)qaws_internal_scratch_zalloc((size_t)m->nc + 1, 1), *on = (unsigned char*)qaws_internal_scratch_zalloc((size_t)m->nv + 1, 1);
	unsigned int *uf = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)m->nf + 1)), v, f, c, n = 0, id;
	qaws_scalar* val = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(qaws_scalar), (size_t)m->nf + 1));
	qaws_attribute_desc d;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (!cut || !on || !uf || !val)
		goto done;
	/* the grid's singular vertices: a feature corner, or a valence not the regular one where the vertex is (3 on the border,
	   4 inside or on a feature line), or a face round it not a quad */
	for (v = 0; v < m->nv; v++)
	{
		unsigned int want = hm_is_border_vertex(m, v) ? 3 : 4, cc;
		if (m->vc[v] == HM_NONE)
			continue;
		on[v] = (unsigned char)((m->vflag[v] & HM_CORNER) || hm_valence(m, v) != want);
		for (cc = m->vc[v]; !on[v] && cc != HM_NONE; cc = hm_round(m, v, cc))
			on[v] = (unsigned char)(m->fn[m->cf[cc]] != 4);
	}
	for (v = 0; v < m->nv; v++)
	{
		unsigned int ring[64], nr, i;
		if (on[v] != 1)
			continue;
		nr = hm_ring(m, v, ring, 64);
		for (i = 0; i < nr && i < 64; i++)
		{
			unsigned int ce = hm_find_edge(m, v, ring[i]);
			if (ce == HM_NONE || m->co[ce] == HM_NONE || (m->cflag[ce] & HM_FEATURE))
				continue;   /* (border and feature edges part patches already) */
			qd_layout_trace(m, cut, on, v, ring[i]);
		}
	}
	for (f = 0; f < m->nf; f++) uf[f] = f;
	for (c = 0; c < m->nc; c++)
	{
		unsigned int a, b2;
		if (m->cf[c] == HM_NONE || m->co[c] == HM_NONE || cut[c] || (m->cflag[c] & HM_FEATURE))
			continue;
		a = qd_find(uf, m->cf[c]);
		b2 = qd_find(uf, m->cf[m->co[c]]);
		if (a != b2)
			uf[a] = b2;
	}
	for (f = 0; f < m->nf; f++) val[f] = (qaws_scalar)-1;
	for (f = 0; f < m->nf; f++)
	{
		unsigned int r = qd_find(uf, f);
		if (val[r] < 0)
			val[r] = (qaws_scalar)n++;
		val[f] = val[r];
	}
	memset(&d, 0, sizeof d);
	d.name = "patch";
	d.element = QAWS_MESH_FACE;
	d.components = 1;
	d.rule = QAWS_ATTRIBUTE_NEAREST;
	st = qaws_mesh_3d_add_attribute(&id, out, val, &d);
	*out_count = n;
done:
	qaws_internal_free(cut);
	qaws_internal_free(on);
	qaws_internal_free(uf);
	qaws_internal_free(val);
	return st;
}

/* The quad mesh at lattice spacing rho from the fields of the dense mesh dm (level 0 of h): the steps rounded, made
   consistent (qaws_quad_flow.h), the positions refitted to them, the mesh read off, its holes filled (the cones on
   lattice points their centres). T and K: room per corner of dm. */
static qaws_status qd_lattice(qaws_hm const* dm, mf_hier* h, double rho, int* T, unsigned char* K, qaws_mesh_reference const* r,
	mr_scratch* s, qd_out* o)
{
	mf_level* L = &h->lv[0];
	unsigned int f, k, odd, ncone = 0;
	double *Tc = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), (4 * (size_t)dm->nc + 4))), *cone = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), (3 * (size_t)dm->nf + 3)));
	unsigned char* fr = (unsigned char*)qaws_internal_scratch_malloc((size_t)dm->nc + 1);
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (!Tc || !cone || !fr)
		goto done;
	qf_steps(dm, L, rho, T, K, Tc);
	if (qf_consistent(dm, L, T, K, Tc, &odd, NULL) == HM_NONE)
		goto done;
	qf_refit(dm, L, rho, T, K, 10);
	for (f = 0; f < dm->nf; f++)
	{
		int S[2], e[6], t = qf_closure(dm, f, T, K, fr, S);
		if (!(t & 1) || ((S[0] + S[1]) & 1) || qf_image(dm, f, T, K, e) >= 0)
			continue;   /* (a cone off the lattice, or on a vertex: no hole) */
		for (k = 0; k < 3; k++)
			cone[3 * (size_t)ncone + k] = (L->x[3 * (size_t)dm->cv[dm->fc[f]] + k] + L->x[3 * (size_t)dm->cv[dm->fc[f] + 1] + k]
				+ L->x[3 * (size_t)dm->cv[dm->fc[f] + 2] + k]) / 3;
		ncone++;
	}
	st = qd_extract(dm, h, T, K, r, s, o);
	if (st == QAWS_STATUS_OK)
		st = qd_fill(o, r, s, cone, ncone, rho);
	if (st == QAWS_STATUS_OK)
		st = qd_pair_triangles(o);
done:
	qaws_internal_free(Tc);
	qaws_internal_free(cone);
	qaws_internal_free(fr);
	return st;
}

/* a few steps of each free vertex toward its neighbours' centre, along the reference */
static void qd_relax(qaws_hm* m, qaws_mesh_reference const* r, mr_scratch* s, unsigned char const* kind, unsigned int steps)
{
	unsigned int it, v, k, ring[64];
	for (it = 0; it < steps; it++)
		for (v = 0; v < m->nv; v++)
		{
			unsigned int n, i, c;
			double cen[3] = { 0, 0, 0 }, d[3], x[3], nrm[3], dn;
			qaws_mesh_location l;
			int bad = 0;
			if (m->vc[v] == HM_NONE || kind[v] == 2)
				continue;
			if (kind[v] == 1)
			{
				/* on a crease or the border: along it, halfway to the middle of its two neighbours there */
				unsigned int nb[2], nn = 0;
				for (c = m->vc[v]; c != HM_NONE && nn <= 2; c = hm_round(m, v, c))
				{
					unsigned int p = hm_prev(m, c), w[2], j;
					w[0] = (m->cflag[c] & HM_FEATURE) || m->co[c] == HM_NONE ? hm_to(m, c) : HM_NONE;
					w[1] = (m->cflag[p] & HM_FEATURE) || m->co[p] == HM_NONE ? m->cv[p] : HM_NONE;
					for (j = 0; j < 2; j++)
						if (w[j] != HM_NONE && (nn == 0 || nb[0] != w[j]) && (nn < 2 || nb[1] != w[j]))
						{
							if (nn < 2)
								nb[nn] = w[j];
							nn++;
						}
				}
				if (nn != 2)
					continue;
				for (k = 0; k < 3; k++)
					cen[k] = 0.5 * m->x[3 * (size_t)v + k] + 0.25 * (m->x[3 * (size_t)nb[0] + k] + m->x[3 * (size_t)nb[1] + k]);
				if (mr_project_crease(r, cen, m->foot[v].patch, s, &l) == HUGE_VAL || l.patch == HM_NONE)
					continue;
				mr_eval(r, &l, x, NULL);
				mr_eval(r, &m->foot[v], NULL, nrm);
			}
			else
			{
				/* inside: halfway to its neighbours' centre, walking on the reference */
				if (hm_is_border_vertex(m, v))
					continue;
				n = hm_ring(m, v, ring, 64);
				if (n == 0 || n > 64)
					continue;
				for (i = 0; i < n; i++)
					for (k = 0; k < 3; k++)
						cen[k] += m->x[3 * (size_t)ring[i] + k] / n;
				mr_eval(r, &m->foot[v], x, nrm);
				for (k = 0; k < 3; k++)
					d[k] = 0.5 * (cen[k] - m->x[3 * (size_t)v + k]);
				dn = qd_dot(d, nrm);
				for (k = 0; k < 3; k++)
					d[k] -= dn * nrm[k];
				if (mr_walk(r, m->foot[v], d, 1, &l))
					continue;
				mr_eval(r, &l, x, NULL);
			}
			/* no face round v turned against the surface (one turned already may be turned back) */
			for (c = m->vc[v]; c != HM_NONE && !bad; c = hm_round(m, v, c))
			{
				double n0[3], n1[3], keep[3];
				hm_face_normal(m, m->cf[c], n0);
				memcpy(keep, m->x + 3 * (size_t)v, sizeof keep);
				memcpy(m->x + 3 * (size_t)v, x, sizeof keep);
				hm_face_normal(m, m->cf[c], n1);
				memcpy(m->x + 3 * (size_t)v, keep, sizeof keep);
				bad = qd_dot(n1, nrm) <= 0 && qd_dot(n0, nrm) > 0;
			}
			if (bad)
				continue;
			memcpy(m->x + 3 * (size_t)v, x, sizeof x);
			m->foot[v] = l;
		}
}

/* The guidance carried onto the dense mesh dm (level 0 of h, normals nrm). Density: each vertex's size at its footpoint,
   the lattice's spacing (rho becomes the smallest, sc the rest over it). Directions: each dense vertex's from the input
   vertices of its footpoint's triangle (crosses matched to the first), a pull weighed by the barycentrics (fading to 0
   at the vertices without one), held where the nearest corner holds. Guide edges: the dense vertices within half a dense
   edge of a segment held along it and on it (moved onto it, MF_GUIDE). out_soft: the pulls (3 per vertex, weight times
   direction) for SMOOTHEST, NULL when none. */
static qaws_status qd_guide(qaws_quadrangulate_desc const* desc, qaws_mesh_3d const* mesh, qaws_mesh_reference const* r, qaws_hm const* dm,
	double const* nrm, mf_hier* h, double* rho, double** out_soft, int** out_cone, double** out_vw)
{
	mf_level* L = &h->lv[0];
	unsigned int n = dm->nv, v, i, k, e, c, any = 0, *band = NULL;
	double *dir = NULL, *w = NULL, *sz = NULL;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	*out_soft = NULL;
	*out_cone = NULL;
	*out_vw = NULL;
	if (desc->size)
	{
		double lo = HUGE_VAL;
		sz = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), (n + 1)));
		if (!sz)
			return st;
		for (v = 0; v < n; v++)
		{
			sz[v] = sf_eval(desc->size, &dm->foot[v], dm->x + 3 * (size_t)v);
			if (!(sz[v] > 0) || sz[v] == HUGE_VAL)
			{
				qaws_internal_free(sz);
				return QAWS_STATUS_INVALID_ARGUMENT;
			}
			lo = sz[v] < lo ? sz[v] : lo;
		}
		*rho = n ? lo : *rho;
		for (v = 0; v < n; v++)
			sz[v] /= *rho;
		k = (unsigned int)mf_scales(h, sz);
		qaws_internal_free(sz);
		if (!k)
			return st;
	}
	if (desc->cone_count || desc->attractor_count)
	{
		/* cones: on the dense vertex the remesh kept at each (pinned there); attractors: the weights at the dense vertices */
		double lo[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL }, hi[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL }, rad = (double)desc->attractor_radius;
		int* cone = desc->cone_count ? (int*)qaws_internal_zalloc((size_t)n + 1, sizeof(int)) : NULL;
		double* vw = desc->attractor_count ? (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 2, (n + 1))) : NULL;
		if ((desc->cone_count && !cone) || (desc->attractor_count && !vw))
		{
			qaws_internal_free(cone);
			qaws_internal_free(vw);
			return st;
		}
		*out_cone = cone;
		*out_vw = vw;
		/* (with a border the flow keeps each cone's star, the border taking what cannot close: the cone stays on its pinned
		   vertex's lattice point. On a closed surface the cones must take it, so they may move a step: their vertices
		   are not pinned, the quad vertex is moved onto them after) */
		for (c = 0; c < dm->nc && (dm->cf[c] == HM_NONE || dm->co[c] != HM_NONE); c++) {}
		for (e = 0; e < desc->cone_count; e++)
		{
			double const* p = mesh->hm.x + 3 * (size_t)desc->cone_vertices[e];
			double best = HUGE_VAL;
			unsigned int bv = 0;
			for (v = 0; v < n; v++)
			{
				double d2 = 0;
				for (k = 0; k < 3; k++)
					d2 += (dm->x[3 * (size_t)v + k] - p[k]) * (dm->x[3 * (size_t)v + k] - p[k]);
				if (d2 < best)
				{
					best = d2;
					bv = v;
				}
			}
			cone[bv] += desc->cone_indices[e];
			for (i = 0; i < desc->pinned_count && desc->pinned[i] != desc->cone_vertices[e]; i++) {}
			L->con[bv] = (unsigned char)(c == dm->nc && i == desc->pinned_count && L->con[bv] == MF_PINNED ? MF_CONE : L->con[bv] | MF_CONE);
		}
		for (v = 0; v < mesh->hm.nv; v++)
			for (k = 0; k < 3; k++)
			{
				lo[k] = fmin(lo[k], mesh->hm.x[3 * (size_t)v + k]);
				hi[k] = fmax(hi[k], mesh->hm.x[3 * (size_t)v + k]);
			}
		if (rad == 0)
			rad = 0.1 * sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) + (hi[2] - lo[2]) * (hi[2] - lo[2]));
		if (vw)
		{
			mf_attract(dm->x, L->area, n, desc->attractors, desc->attractor_strengths, desc->attractor_count, rad > 0 ? rad : 1, vw);
			if (!mf_weights(h, vw))
				return st;
		}
	}
	if (!desc->guide_directions && !desc->guide_edge_count)
		return QAWS_STATUS_OK;
	dir = (double*)qaws_internal_zalloc(3 * (size_t)n + 3, sizeof(double));
	w = (double*)qaws_internal_scratch_zalloc((size_t)n + 1, sizeof(double));
	band = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (n + 1)));
	if (!dir || !w || !band)
		goto done;
	for (v = 0; desc->guide_directions && v < n; v++)
	{
		qaws_mesh_location const* l = &dm->foot[v];
		unsigned int const* tv = r->tv + 3 * (size_t)l->patch;
		double bc[3], ref[3], ws = 0, best = -1;
		int have = 0, hard = 0;
		bc[0] = 1 - (double)l->a - (double)l->b;
		bc[1] = (double)l->a;
		bc[2] = (double)l->b;
		for (i = 0; i < 3; i++)
		{
			double wi = desc->guide_weights ? (double)desc->guide_weights[tv[i]] : 1, d[3], a[3], b[3], dn;
			int held = !(wi > 0 && wi < HUGE_VAL);
			for (k = 0; k < 3; k++)
				d[k] = (double)desc->guide_directions[3 * (size_t)tv[i] + k];
			dn = qd_dot(d, nrm + 3 * (size_t)v);
			for (k = 0; k < 3; k++)
				d[k] -= dn * nrm[3 * (size_t)v + k];
			dn = sqrt(qd_dot(d, d));
			if (bc[i] > best)
			{
				best = bc[i];
				hard = 0;
			}
			if (wi == 0 || wi != wi || !(dn > 0))
				continue;
			hard = bc[i] == best && held;
			for (k = 0; k < 3; k++)
				d[k] /= dn;
			if (!have)
				memcpy(ref, d, sizeof ref);
			have = 1;
			mf_compat_orientation(ref, nrm + 3 * (size_t)v, d, nrm + 3 * (size_t)v, a, b);
			for (k = 0; k < 3; k++)
				dir[3 * (size_t)v + k] += bc[i] * (held ? 1 : wi) * b[k];
			ws += bc[i] * (held ? 1 : wi);
		}
		w[v] = hard ? -1 : ws;
		any |= w[v] != 0;
	}
	/* guide edges: the nearest segment within half a dense edge (a third of the lattice's spacing), held along it.
	   shortcut: every vertex against every segment (with a box test); a grid of the segments when guides get long */
	for (v = 0; v < n; v++)
	{
		double const* x = dm->x + 3 * (size_t)v;
		double tol = 0.5 * mf_pair_scale(L, v, v, *rho) / 3, bd = tol * tol, bp[3] = { 0, 0, 0 }, bdir[3] = { 0, 0, 0 };
		band[v] = 0;
		for (e = 0; e < desc->guide_edge_count; e++)
		{
			double const *pa = mesh->hm.x + 3 * (size_t)desc->guide_edges[2 * e], *pb = mesh->hm.x + 3 * (size_t)desc->guide_edges[2 * e + 1];
			double ab[3], t, d2 = 0, p[3];
			for (k = 0; k < 3; k++)
				if (x[k] < fmin(pa[k], pb[k]) - tol || x[k] > fmax(pa[k], pb[k]) + tol)
					break;
			if (k < 3)
				continue;
			for (k = 0; k < 3; k++)
				ab[k] = pb[k] - pa[k];
			t = qd_dot(ab, ab) > 0 ? ((x[0] - pa[0]) * ab[0] + (x[1] - pa[1]) * ab[1] + (x[2] - pa[2]) * ab[2]) / qd_dot(ab, ab) : 0;
			t = t < 0 ? 0 : t > 1 ? 1 : t;
			for (k = 0; k < 3; k++)
			{
				p[k] = pa[k] + t * ab[k];
				d2 += (p[k] - x[k]) * (p[k] - x[k]);
			}
			if (d2 <= bd && qd_dot(ab, ab) > 0)
			{
				bd = d2;
				memcpy(bp, p, sizeof bp);
				memcpy(bdir, ab, sizeof bdir);
				band[v] = 1;
			}
		}
		if (!band[v] || (L->con[v] & MF_ALONG))
		{
			band[v] = 0;
			continue;   /* (a feature line wins) */
		}
		memcpy(dir + 3 * (size_t)v, bdir, sizeof bdir);
		w[v] = -1;
		memcpy(L->x + 3 * (size_t)v, bp, sizeof bp);
		any = 1;
	}
	if (!any)
	{
		st = QAWS_STATUS_OK;
		goto done;
	}
	mf_guide(h, dir, w);
	for (v = 0; v < n; v++)
		if (band[v])
			L->con[v] |= MF_ON_LINE | MF_GUIDE;
	/* the pulls, as mf_guide left them */
	for (v = 0; v < n && !(L->con[v] & MF_PREFER); v++) {}
	if (v < n)
	{
		*out_soft = dir;
		for (v = 0; v < 3 * n; v++)
			dir[v] = (L->con[v / 3] & MF_PREFER) ? L->cq[v] : 0;
		dir = NULL;
	}
	st = QAWS_STATUS_OK;
done:
	qaws_internal_free(dir);
	qaws_internal_free(w);
	qaws_internal_free(band);
	return st;
}

/* each placed cone's lattice point on its dense vertex, its neighbours' on it too: the steps round the cone's star 0, so
   the cone sits on that point whatever turn its face makes (the point a rotation about it keeps) */
static void qd_cone_points(mf_hier* h, int const* cone)
{
	mf_level* L = &h->lv[0];
	unsigned int v, k;
	for (v = 0; cone && v < L->n; v++)
		if (cone[v])
		{
			memcpy(L->o + 3 * (size_t)v, L->x + 3 * (size_t)v, sizeof(double) * 3);
			for (k = L->adj0[v]; k < L->adj0[v + 1]; k++)
				memcpy(L->o + 3 * (size_t)L->adj[k], L->x + 3 * (size_t)v, sizeof(double) * 3);
		}
}

/* qaws_mesh_quadrangulate; map (may be NULL): GLOBAL's integer-grid map kept (left empty when it fell back) */
static qaws_status qd_run(qaws_mesh_3d const* mesh, qaws_quadrangulate_desc const* desc, qaws_mesh_3d** out_mesh, qaws_mesh_report* out_report,
	qg_map* map)
{
	qaws_mesh_reference* own = NULL;
	qaws_mesh_reference const* r;
	qaws_mesh_3d *dense = NULL, *out = NULL;
	qaws_remesh_desc rd;
	qaws_mesh_report rep;
	mf_hier h;
	mr_scratch s;
	qd_out o;
	qaws_mesh_location* fl = NULL;
	double rho, rf, *nrm = NULL, *soft = NULL, *vw = NULL;
	int* dcone = NULL;
	unsigned int* pins = NULL;
	qaws_size_field* dsize = NULL;
	int* T = NULL;
	unsigned char* K = NULL;
	unsigned int v, f, c, k;
	int direct;
	qaws_status st;
	if (!out_mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_mesh = NULL;
	if (!mesh || !desc || (unsigned int)desc->method > QAWS_QUADS_LAYOUT || desc->target_length < 0)
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (v = 0; v < desc->pinned_count; v++)
		if (!desc->pinned || desc->pinned[v] >= mesh->hm.nv)
			return QAWS_STATUS_INVALID_ARGUMENT;
	for (v = 0; v < 2 * desc->guide_edge_count; v++)
		if (!desc->guide_edges || desc->guide_edges[v] >= mesh->hm.nv)
			return QAWS_STATUS_INVALID_ARGUMENT;
	for (v = 0; v < desc->cone_count; v++)
		if (!desc->cone_vertices || !desc->cone_indices || desc->cone_vertices[v] >= mesh->hm.nv)
			return QAWS_STATUS_INVALID_ARGUMENT;
	if ((desc->attractor_count && (!desc->attractors || !desc->attractor_strengths)) || !(desc->attractor_radius >= 0))
		return QAWS_STATUS_INVALID_ARGUMENT;
	memset(&h, 0, sizeof h);
	memset(&s, 0, sizeof s);
	memset(&o, 0, sizeof o);
	memset(&rep, 0, sizeof rep);
	r = desc->reference;
	if (!r)
	{
		st = qaws_mesh_reference_create(&own, mesh, NULL);
		if (st != QAWS_STATUS_OK)
			return st;
		r = own;
	}
	if (desc->guide_directions && r->mesh->hm.nv != mesh->hm.nv)
	{
		qaws_mesh_reference_destroy(own);
		return QAWS_STATUS_INVALID_ARGUMENT;   /* (directions per input vertex, carried through the reference's triangles) */
	}
	rho = (double)desc->target_length;
	if (rho == 0)
	{
		qaws_mesh_quality q;
		qaws_mesh_3d_measure_quality(&q, mesh);
		rho = q.mean_edge > 0 ? 2 * (double)q.mean_edge : 1;
	}
	rf = rho;
	direct = desc->method == QAWS_QUADS_PAIRING || desc->method == QAWS_QUADS_SPLIT;
	if (direct)
	{
		/* PAIRING and SPLIT: the input's own faces */
		if (!mr_scratch_init(&s, r->nt))
		{
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			goto done;
		}
		st = qd_direct(mesh, desc, r, &s, &o);
		if (st != QAWS_STATUS_OK)
			goto done;
		goto output;
	}
	/* 1. a dense sampling of the surface */
	memset(&rd, 0, sizeof rd);
	rd.reference = r;
	rd.target_length = (qaws_scalar)(rf / 3);
	rd.iterations = 6;
	if (desc->cone_count)
	{
		/* the cones' vertices kept by the dense remesh too */
		pins = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), ((size_t)desc->pinned_count + desc->cone_count + 1)));
		if (!pins)
		{
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			goto done;
		}
		if (desc->pinned_count)
			memcpy(pins, desc->pinned, sizeof(unsigned int) * desc->pinned_count);
		memcpy(pins + desc->pinned_count, desc->cone_vertices, sizeof(unsigned int) * desc->cone_count);
	}
	rd.pinned = pins ? pins : desc->pinned;
	rd.pinned_count = desc->pinned_count + desc->cone_count;
	st = desc->size ? sf_create_scaled(desc->size, 1.0 / 3, &dsize) : QAWS_STATUS_OK;
	rd.size = dsize;
	if (st == QAWS_STATUS_OK)
		st = qaws_mesh_remesh(&dense, NULL, mesh, &rd);
	if (st != QAWS_STATUS_OK)
		goto done;
	/* 2. the fields, with the reference's normals */
	nrm = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (dense->hm.nv + 1)));
	if (!nrm || !mr_scratch_init(&s, r->nt))
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	for (v = 0; v < dense->hm.nv; v++)
		mr_eval(r, &dense->hm.foot[v], NULL, nrm + 3 * (size_t)v);
	st = mf_build(&h, &dense->hm, nrm);
	if (st != QAWS_STATUS_OK)
		goto done;
	st = qd_guide(desc, mesh, r, &dense->hm, nrm, &h, &rho, &soft, &dcone, &vw);
	if (st != QAWS_STATUS_OK)
		goto done;
	rf = rho;
	{
		/* the cross field twice from the direction of largest curvature: free (Instant Meshes, along the features), and drawn
		   to it on every level, as strongly as the surface bends one way more than the other (1 round a cylinder
		   or a torus, none where it bends alike or hardly at the lattice's scale); the one with fewer singularities kept
		   (a closed tube wants the pull, a patch held by its border not). Then twice globally (Knoppel et al. 2013, one
		   factor): the smoothest field, kept when it has no more singularities than the best so far, and the one drawn to
		   the same direction, kept when it has fewer; the singularities sit where the whole surface wants them, not where
		   a coarse start left them. Held by lines, either is first relaxed at the lattice's scale (Ginzburg-Landau), which
		   parts the zeros the relaxed energy lets coincide: a disc held by its border gets four cones well inside it (an
		   O-grid), neither one at its centre nor four on its rim. Any failure there keeps the field there is */
		double* init = (double*)qaws_internal_scratch_zalloc(3 * (size_t)dense->hm.nv + 3, sizeof(double));
		double *keep = NULL, *cand = NULL;
		unsigned int pulled = 0, s0, l, t, placed = 0;
		qaws_status fst = QAWS_STATUS_OK;
		size_t at, total;
		cf_system* sys = NULL;
		if (!init)
		{
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			goto done;
		}
		for (v = 0; v < dense->hm.nv; v++)
		{
			qaws_mesh_location const* l0 = &dense->hm.foot[v];
			double xa[3], xb[3], n[3], na[3], nb[3], k1, k2, d1[3], d2[3], an;
			mi_eval(&r->patch[l0->patch], (double)l0->a, (double)l0->b, NULL, xa, xb, NULL, NULL, NULL, n, na, nb);
			mi_curvature(xa, xb, n, na, nb, &k1, &k2, d1, d2);
			an = (fabs(k1) + fabs(k2)) * rf > 1e-3 ? (k1 - k2) / (fabs(k1) + fabs(k2)) : 0;
			if (an > 0.15)
			{
				for (k = 0; k < 3; k++)
					init[3 * (size_t)v + k] = d1[k] * an;
				pulled++;
			}
			if (soft && (soft[3 * (size_t)v] != 0 || soft[3 * (size_t)v + 1] != 0 || soft[3 * (size_t)v + 2] != 0))
				memcpy(init + 3 * (size_t)v, soft + 3 * (size_t)v, sizeof(double) * 3);   /* (the user's pull, not the curvature's) */
		}
		mf_orientations(&h, desc->sweeps ? desc->sweeps : 8, init, 0);
		for (total = 0, l = 0; l < h.nl; l++)
			total += 3 * (size_t)h.lv[l].n;
		keep = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), (total + 1)));
		cand = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), (3 * (size_t)dense->hm.nv + 3)));
		if (!keep || !cand)
		{
			qaws_internal_free(init);
			qaws_internal_free(keep);
			qaws_internal_free(cand);
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			goto done;
		}
		for (at = 0, l = 0; l < h.nl; at += 3 * (size_t)h.lv[l].n, l++)
			memcpy(keep + at, h.lv[l].q, sizeof(double) * 3 * h.lv[l].n);
		s0 = qd_singular(&dense->hm, &h.lv[0]);
		/* each candidate: kept (into keep) when it has fewer singularities (the smoothest: no more), else undone */
		for (t = 0; t < 3; t++)
		{
			unsigned int s1;
			int better;
			if (t == 0)
			{
				if (!pulled || dcone)
					continue;
				mf_orientations(&h, desc->sweeps ? desc->sweeps : 8, init, 1);
			}
			else
			{
				if (t == 1 && ((fst = cf_create(&dense->hm, nrm, &h.lv[0], 4, &sys)) != QAWS_STATUS_OK
					|| (fst = cf_steer(sys, vw ? vw + dense->hm.nv : NULL, dcone)) != QAWS_STATUS_OK
					|| (soft && (fst = cf_guide(sys, soft)) != QAWS_STATUS_OK)))
					break;
				if (t == 2 && !pulled)
					break;
				if ((fst = t == 1 ? cf_smoothest(sys, cand) : cf_aligned(sys, init, cand)) != QAWS_STATUS_OK)
					continue;
				if ((fst = cf_relax(sys, rf, 100, cand)) != QAWS_STATUS_OK)
					continue;
				memcpy(h.lv[0].q, cand, sizeof(double) * 3 * dense->hm.nv);
			}
			s1 = qd_singular(&dense->hm, &h.lv[0]);
			better = t == 1 ? s1 <= s0 || dcone : s1 < s0;   /* (placed cones: the field with them, whatever it costs) */
			placed |= better && t;
			if (better && t)
				cf_restrict(&h);
			s0 = better ? s1 : s0;
			for (at = 0, l = 0; l < h.nl; at += 3 * (size_t)h.lv[l].n, l++)
				if (better)
					memcpy(keep + at, h.lv[l].q, sizeof(double) * 3 * h.lv[l].n);
				else
					memcpy(h.lv[l].q, keep + at, sizeof(double) * 3 * h.lv[l].n);
		}
		cf_free(sys);
		qaws_internal_free(keep);
		qaws_internal_free(cand);
		qaws_internal_free(init);
		if (dcone && !placed)
		{
			st = fst != QAWS_STATUS_OK ? fst : QAWS_STATUS_INTERNAL_ERROR;   /* (the cones refused: Poincare-Hopf) */
			goto done;
		}
	}
	mf_positions(&h, rf, desc->sweeps ? desc->sweeps : 8);
	qd_cone_points(&h, dcone);
	/* 3. the steps between lattice points made to close round every face, the mesh read off them; all quads at the
	   target when that leaves no other face, else (a triangle or a polygon left) the same at twice the target, every face
	   split; holes left are not looked for, the split would not close them either */
	T = (int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(int), (2 * (size_t)dense->hm.nc + 2)));
	K = (unsigned char*)qaws_internal_scratch_malloc((size_t)dense->hm.nc + 1);
	if (!T || !K)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	st = qd_lattice(&dense->hm, &h, rf, T, K, r, &s, &o);
	if (st == QAWS_STATUS_OK && (desc->method == QAWS_QUADS_GLOBAL || desc->method == QAWS_QUADS_LAYOUT))
	{
		/* GLOBAL: the integer-grid map of the same steps, its quads placed on the reference; INSTANT's kept when it gives
		   none */
		qd_out g;
		st = qg_quads(&dense->hm, &h.lv[0], rf, T, K, &g, map);
		if (st == QAWS_STATUS_OK)
		{
			for (v = 0; v < g.nv; v++)
				if (g.kind[v] != 2)
				{
					double p[3];
					memcpy(p, g.x + 3 * (size_t)v, sizeof p);
					qd_place(r, &s, p, g.foot[v], g.kind[v] != 0, g.x + 3 * (size_t)v, &g.foot[v]);
				}
			qd_out_free(&o);
			o = g;
			rep.iterations = 1;
		}
		else if (st == QAWS_STATUS_INVALID_TOPOLOGY)
			st = QAWS_STATUS_OK;
	}
	for (f = 0; st == QAWS_STATUS_OK && desc->pure && !rep.iterations && f < o.nf && o.fs[f] == 4; f++) {}
	if (st == QAWS_STATUS_OK && desc->pure && !rep.iterations && f < o.nf)
	{
		qd_out o2;
		memset(&o2, 0, sizeof o2);
		mf_positions(&h, 2 * rho, desc->sweeps ? desc->sweeps : 8);
		qd_cone_points(&h, dcone);
		st = qd_lattice(&dense->hm, &h, 2 * rho, T, K, r, &s, &o2);
		/* (nothing at twice a target about the shape's own size: the first mesh split instead) */
		if (st == QAWS_STATUS_OK && o2.nf)
		{
			qd_out_free(&o);
			o = o2;
		}
		else
			qd_out_free(&o2);
		if (st == QAWS_STATUS_OK)
			st = qd_split(&o, r, &s);
	}
	if (st != QAWS_STATUS_OK)
		goto done;
output:
	st = hm_mesh_3d_create(o.x, o.nv, o.fv, o.fs, o.nf, &out);
	if (st != QAWS_STATUS_OK)
		goto done;
	if (!hm_foot_enable(&out->hm))
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	{
		/* copies split off at creation take their sources' footpoints and kinds */
		unsigned char* kind = (unsigned char*)qaws_internal_scratch_malloc(out->hm.nv + 1);
		if (!kind)
		{
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			goto done;
		}
		for (v = 0; v < out->hm.nv; v++)
		{
			unsigned int src = out->hm.vsrc[v] < o.nv ? out->hm.vsrc[v] : 0;
			out->hm.foot[v] = o.foot[src];
			kind[v] = o.kind[src];
			if (kind[v] == 2)
				out->hm.vflag[v] |= HM_CORNER;
		}
		/* each cone's quad vertex: the nearest within the spacing of its valence (4 - index), else the nearest; moved onto
		   the cone's vertex and held there like a corner (not flagged one) by the relaxation and the finishing */
		for (c = 0; !direct && c < desc->cone_count; c++)
		{
			double const* p = mesh->hm.x + 3 * (size_t)desc->cone_vertices[c];
			double best[2] = { HUGE_VAL, HUGE_VAL };
			unsigned int bv[2] = { HM_NONE, HM_NONE }, want = (unsigned int)(4 - desc->cone_indices[c]);
			for (v = 0; v < out->hm.nv; v++)
			{
				double d2 = 0;
				if (out->hm.vc[v] == HM_NONE)
					continue;
				for (k = 0; k < 3; k++)
					d2 += (out->hm.x[3 * (size_t)v + k] - p[k]) * (out->hm.x[3 * (size_t)v + k] - p[k]);
				if (kind[v] == 2 && d2 == 0)
				{
					for (k = 0; k < desc->pinned_count && desc->pinned[k] != desc->cone_vertices[c]; k++) {}
					if (k == desc->pinned_count)
						out->hm.vflag[v] = (unsigned char)(out->hm.vflag[v] & ~HM_CORNER);
					break;   /* (pinned there: the dense vertex held its lattice point) */
				}
				for (f = 0; kind[v] != 2 && f < 2; f++)
					if (d2 < best[f] && (f || (hm_valence(&out->hm, v) == want && !hm_is_border_vertex(&out->hm, v) && d2 < rf * rf)))
					{
						best[f] = d2;
						bv[f] = v;
					}
			}
			v = v < out->hm.nv ? HM_NONE : bv[0] != HM_NONE ? bv[0] : bv[1];
			if (v == HM_NONE)
				continue;
			qd_place(r, &s, p, out->hm.foot[v], 0, out->hm.x + 3 * (size_t)v, &out->hm.foot[v]);
			kind[v] = 2;
		}
		/* border vertices near the reference's border but off it (a lattice point by the border whose cluster had no
		   border vertex, or a midpoint next to one): onto it */
		for (v = 0; v < out->hm.nv; v++)
		{
			qaws_mesh_location l;
			if (kind[v] || !hm_is_border_vertex(&out->hm, v))
				continue;
			if (mr_project_crease(r, out->hm.x + 3 * (size_t)v, out->hm.foot[v].patch, &s, &l) < 0.25 * rf * rf && l.patch != HM_NONE)
			{
				mr_eval(r, &l, out->hm.x + 3 * (size_t)v, NULL);
				out->hm.foot[v] = l;
				kind[v] = 1;
			}
		}
		/* feature edges: between feature vertices, along a crease */
		for (c = 0; c < out->hm.nc; c++)
		{
			unsigned int a, b;
			double mid[3], d2 = 0;
			qaws_mesh_location l;
			if (out->hm.cf[c] == HM_NONE)
				continue;
			a = out->hm.cv[c];
			b = hm_to(&out->hm, c);
			if (!kind[a] || !kind[b])
				continue;
			for (k = 0; k < 3; k++)
			{
				mid[k] = 0.5 * (out->hm.x[3 * (size_t)a + k] + out->hm.x[3 * (size_t)b + k]);
				d2 += (out->hm.x[3 * (size_t)a + k] - out->hm.x[3 * (size_t)b + k]) * (out->hm.x[3 * (size_t)a + k] - out->hm.x[3 * (size_t)b + k]);
			}
			if (mr_project_crease(r, mid, out->hm.foot[a].patch, &s, &l) < 0.01 * d2)
				out->hm.cflag[c] |= HM_FEATURE;
		}
		for (v = 0; v < out->hm.nv; v++)
			kind[v] = kind[v] == 3 ? 2 : kind[v];   /* (guide lines held by the relaxation and the finishing, not corners) */
		if (!direct)
		{
			qd_relax(&out->hm, r, &s, kind, 10);
			if (qi_improve(&out->hm, r, &s, kind, rho, 4) == HM_NONE)
			{
				qaws_internal_free(kind);
				st = QAWS_STATUS_ALLOCATION_FAILURE;
				goto done;
			}
		}
		qaws_internal_free(kind);
	}
	/* 4. attributes and the report */
	fl = (qaws_mesh_location*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_mesh_location), (out->hm.nf + 1)));
	if (!fl)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	for (f = 0; f < out->hm.nf; f++)
	{
		double cen[3] = { 0, 0, 0 }, d;
		qaws_mesh_location hint = out->hm.foot[out->hm.cv[out->hm.fc[f]]];
		for (c = out->hm.fc[f]; c < out->hm.fc[f] + out->hm.fn[f]; c++)
			for (k = 0; k < 3; k++)
				cen[k] += out->hm.x[3 * (size_t)out->hm.cv[c] + k] / out->hm.fn[f];
		d = sqrt(mr_project(r, cen, hint.patch, &s, &fl[f]));
		rep.max_deviation = (qaws_scalar)(d > (double)rep.max_deviation ? d : (double)rep.max_deviation);
		rep.mean_deviation += (qaws_scalar)(d / out->hm.nf);
	}
	st = mx_transfer(&out->hm, r, fl);
	if (st == QAWS_STATUS_OK && desc->method == QAWS_QUADS_LAYOUT)
		st = qd_patches(out, &rep.splits);
	if (st == QAWS_STATUS_OK)
		qaws_mesh_3d_measure_quality(&rep.quality, out);
done:
	qaws_internal_free(fl);
	qaws_internal_free(nrm);
	qaws_internal_free(soft);
	qaws_internal_free(vw);
	qaws_internal_free(dcone);
	qaws_internal_free(pins);
	qaws_size_field_destroy(dsize);
	qaws_internal_free(T);
	qaws_internal_free(K);
	qd_out_free(&o);
	mf_free(&h);
	mr_scratch_free(&s);
	qaws_mesh_3d_destroy(dense);
	qaws_mesh_reference_destroy(own);
	if (st != QAWS_STATUS_OK)
	{
		qaws_mesh_3d_destroy(out);
		return st;
	}
	if (out_report)
		*out_report = rep;
	*out_mesh = out;
	return QAWS_STATUS_OK;
}

qaws_status qaws_mesh_quadrangulate(qaws_mesh_3d** out_mesh, qaws_mesh_report* out_report, qaws_mesh_3d const* mesh, qaws_quadrangulate_desc const* desc)
{
	return qd_run(mesh, desc, out_mesh, out_report, NULL);
}

qaws_status qg_parametrization(qaws_mesh_3d const* mesh, qaws_quadrangulate_desc const* desc, qg_map* out)
{
	qaws_mesh_3d* q = NULL;
	qaws_status st;
	memset(out, 0, sizeof *out);
	if (!desc || desc->method != QAWS_QUADS_GLOBAL)
		return QAWS_STATUS_INVALID_ARGUMENT;
	st = qd_run(mesh, desc, &q, NULL, out);
	qaws_mesh_3d_destroy(q);
	if (st == QAWS_STATUS_OK && !out->nt)
		st = QAWS_STATUS_INVALID_TOPOLOGY;
	if (st != QAWS_STATUS_OK)
		qg_map_free(out);
	return st;
}
