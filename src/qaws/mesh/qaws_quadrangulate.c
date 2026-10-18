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
#include "../qaws_quadrangulate.h"
#include "../qaws_remesh.h"
#include "qaws_cross_field.h"
#include "qaws_quad_flow.h"
#include "qaws_quad_improve.h"
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
	unsigned int *first = (unsigned int*)malloc(sizeof(unsigned int) * (nv + 2)), *inc = (unsigned int*)malloc(sizeof(unsigned int) * (*nc + 1));
	unsigned int* at = (unsigned int*)malloc(sizeof(unsigned int) * (*nf + 1));
	unsigned char* mark = (unsigned char*)malloc(*nf + 1);   /* (1: changed this pass, 2: joined away) */
	unsigned int f, k, v, changed;
	if (!first || !inc || !at || !mark)
	{
		free(first); free(inc); free(at); free(mark);
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
	free(first); free(inc); free(at); free(mark);
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

/* The quad mesh read off the fields of the dense mesh dm (level 0 of h):
   vertices x (3 each) with footpoints, kind (2: pinned, 1: on a line), faces. */
typedef struct qd_out
{
	double* x;
	qaws_mesh_location* foot;
	unsigned char* kind;
	unsigned int nv;
	unsigned int *fv, *fs, nf, nc;
} qd_out;

static void qd_out_free(qd_out* o)
{
	free(o->x); free(o->foot); free(o->kind); free(o->fv); free(o->fs);
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
	double *cx = NULL, *cw = NULL;
	unsigned char *ckind = NULL, *paired = NULL;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	memset(out, 0, sizeof *out);
	uf = (unsigned int*)malloc(sizeof(unsigned int) * (n + 1));
	cl = (unsigned int*)malloc(sizeof(unsigned int) * (n + 1));
	member = (unsigned int*)malloc(sizeof(unsigned int) * (n + 1));
	cx = (double*)calloc(3 * (size_t)n + 3, sizeof(double));
	cw = (double*)calloc((size_t)n + 1, sizeof(double));
	ckind = (unsigned char*)calloc((size_t)n + 1, 1);
	tri = (unsigned int*)malloc(sizeof(unsigned int) * (4 * (size_t)dm->nf + 4));
	paired = (unsigned char*)calloc((size_t)dm->nf + 1, 1);
	if (!uf || !cl || !member || !cx || !cw || !ckind || !tri || !paired)
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
		else if ((L->con[i] & MF_ON_LINE) && ckind[q] < 1)
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
	for (mask = 1; mask < 2 * nt + 2; mask <<= 1) {}
	hash = (unsigned int*)malloc(sizeof(unsigned int) * mask);
	out->fv = (unsigned int*)malloc(sizeof(unsigned int) * (4 * (size_t)nt + 1));
	out->fs = (unsigned int*)malloc(sizeof(unsigned int) * ((size_t)nt + 1));
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
	used = (unsigned int*)malloc(sizeof(unsigned int) * (ncl + 1));
	if (!used)
		goto done;
	for (i = 0; i < ncl; i++)
		used[i] = HM_NONE;
	for (i = 0; i < out->nc; i++)
		if (used[out->fv[i]] == HM_NONE)
			used[out->fv[i]] = out->nv++;
	out->x = (double*)malloc(sizeof(double) * 3 * (out->nv + 1));
	out->foot = (qaws_mesh_location*)malloc(sizeof(qaws_mesh_location) * (out->nv + 1));
	out->kind = (unsigned char*)malloc(out->nv + 1);
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
		else
			qd_place(r, s, cx + 3 * (size_t)i, dm->foot[member[i]], ckind[i] != 0, out->x + 3 * (size_t)q, &out->foot[q]);
		out->kind[q] = ckind[i];
	}
	for (i = 0; i < out->nc; i++)
		out->fv[i] = used[out->fv[i]];
	st = QAWS_STATUS_OK;
done:
	free(uf); free(cl); free(member); free(tri); free(hash); free(used);
	free(cx); free(cw); free(ckind); free(paired);
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
	unsigned int f, k, c0, mask, *h, *start, nadd = 0, added = 0, *nf, *nfs, nv0 = o->nv;
	qaws_status st = QAWS_STATUS_OK;
	for (mask = 1; mask < 2 * o->nc + 2; mask <<= 1) {}
	h = (unsigned int*)malloc(sizeof(unsigned int) * 2 * mask);
	start = (unsigned int*)malloc(sizeof(unsigned int) * (o->nv + 1));
	if (!h || !start)
	{
		free(h);
		free(start);
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
	nf = (unsigned int*)malloc(sizeof(unsigned int) * (o->nc + 8 * (size_t)o->nv + 1));
	nfs = (unsigned int*)malloc(sizeof(unsigned int) * (o->nf + o->nv + 1));
	{
		/* (room for a centre per loop of six or more, the loops sharing no vertex) */
		size_t nmax = (size_t)o->nv + o->nv / 6 + 1;
		double* x = (double*)realloc(o->x, sizeof(double) * 3 * nmax);
		qaws_mesh_location* foot;
		unsigned char* kind;
		if (x)
			o->x = x;
		foot = x ? (qaws_mesh_location*)realloc(o->foot, sizeof(qaws_mesh_location) * nmax) : NULL;
		if (foot)
			o->foot = foot;
		kind = foot ? (unsigned char*)realloc(o->kind, nmax) : NULL;
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
		for (j = 0; j < len; j++)
			start[loop[j]] = HM_NONE;   /* (each hole once) */
		if (!ok || len < 3 || len > 10)
			continue;
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
	free(o->fv);
	free(o->fs);
	o->fv = nf;
	o->fs = nfs;
	o->nf += nadd;
	o->nc += added;
	nf = nfs = NULL;
done:
	free(nf);
	free(nfs);
	free(h);
	free(start);
	return st;
}

/* every face split into quads round its centre through its edges' midpoints */
static qaws_status qd_split(qd_out* o, qaws_mesh_reference const* r, mr_scratch* s)
{
	unsigned int f, k, q, nv2, *mid, *keys, c0, nf2 = 0, mask;
	qd_out n;
	memset(&n, 0, sizeof n);
	for (mask = 1; mask < 2 * o->nc + 2; mask <<= 1) {}
	mid = (unsigned int*)malloc(sizeof(unsigned int) * (o->nc + 1));
	keys = (unsigned int*)malloc(sizeof(unsigned int) * 3 * mask);   /* (a hash: low end, high end, midpoint) */
	n.x = (double*)malloc(sizeof(double) * 3 * (o->nv + o->nc + o->nf + 1));
	n.foot = (qaws_mesh_location*)malloc(sizeof(qaws_mesh_location) * (o->nv + o->nc + o->nf + 1));
	n.kind = (unsigned char*)calloc((size_t)o->nv + o->nc + o->nf + 1, 1);
	n.fv = (unsigned int*)malloc(sizeof(unsigned int) * 4 * (o->nc + 1));
	n.fs = (unsigned int*)malloc(sizeof(unsigned int) * (o->nc + 1));
	if (!mid || !keys || !n.x || !n.foot || !n.kind || !n.fv || !n.fs)
	{
		free(mid);
		free(keys);
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
				if (o->kind[a] && o->kind[b])
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
				n.kind[nv2] = (unsigned char)along;
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
	free(mid);
	free(keys);
	qd_out_free(o);
	*o = n;
	return QAWS_STATUS_OK;
}

/* triangles that share a side merged into a quad each (a half cell and its neighbour across a fold, a hole's filling).
   shortcut: pairs found by a scan over all triangle pairs, fine for the few left; a hash of their sides if many are */
static qaws_status qd_pair_triangles(qd_out* o)
{
	unsigned int f, g, k, j, c0, nf = 0, nc = 0;
	unsigned int *at = (unsigned int*)malloc(sizeof(unsigned int) * ((size_t)o->nf + 1)), *mate = (unsigned int*)malloc(sizeof(unsigned int) * ((size_t)o->nf + 1));
	unsigned int *nfv = (unsigned int*)malloc(sizeof(unsigned int) * ((size_t)o->nc + 1)), *nfs = (unsigned int*)malloc(sizeof(unsigned int) * ((size_t)o->nf + 1));
	unsigned char* side = (unsigned char*)malloc((size_t)o->nf + 1);
	if (!at || !mate || !nfv || !nfs || !side)
	{
		free(at); free(mate); free(nfv); free(nfs); free(side);
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
					if (u[j] == t[(k + 1) % 3] && u[(j + 1) % 3] == t[k])
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
	free(o->fv);
	free(o->fs);
	o->fv = nfv;
	o->fs = nfs;
	o->nf = nf;
	o->nc = nc;
	free(at); free(mate); free(side);
	return QAWS_STATUS_OK;
}

/* The quad mesh at lattice spacing rho from the fields of the dense mesh dm (level 0 of h): the steps rounded, made
   consistent (qaws_quad_flow.h), the positions refitted to them, the mesh read off, its holes filled (the cones on
   lattice points their centres). T and K: room per corner of dm. */
static qaws_status qd_lattice(qaws_hm const* dm, mf_hier* h, double rho, int* T, unsigned char* K, qaws_mesh_reference const* r,
	mr_scratch* s, qd_out* o)
{
	mf_level* L = &h->lv[0];
	unsigned int f, k, odd, ncone = 0;
	double *Tc = (double*)malloc(sizeof(double) * (4 * (size_t)dm->nc + 4)), *cone = (double*)malloc(sizeof(double) * (3 * (size_t)dm->nf + 3));
	unsigned char* fr = (unsigned char*)malloc((size_t)dm->nc + 1);
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
	free(Tc);
	free(cone);
	free(fr);
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

qaws_status qaws_mesh_quadrangulate(qaws_mesh_3d const* mesh, qaws_quadrangulate_desc const* desc, qaws_mesh_3d** out_mesh, qaws_mesh_report* out_report)
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
	double rho, rf, *nrm = NULL;
	int* T = NULL;
	unsigned char* K = NULL;
	unsigned int v, f, c, k;
	qaws_status st;
	if (!out_mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_mesh = NULL;
	if (!mesh || !desc || (unsigned int)desc->method > QAWS_QUADS_LAYOUT || desc->target_length < 0)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (desc->method != QAWS_QUADS_INSTANT)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	memset(&h, 0, sizeof h);
	memset(&s, 0, sizeof s);
	memset(&o, 0, sizeof o);
	memset(&rep, 0, sizeof rep);
	r = desc->reference;
	if (!r)
	{
		st = qaws_mesh_reference_create(mesh, NULL, &own);
		if (st != QAWS_STATUS_OK)
			return st;
		r = own;
	}
	rho = (double)desc->target_length;
	if (rho == 0)
	{
		qaws_mesh_quality q;
		qaws_mesh_3d_measure_quality(mesh, &q);
		rho = q.mean_edge > 0 ? 2 * (double)q.mean_edge : 1;
	}
	rf = rho;
	/* 1. a dense sampling of the surface */
	memset(&rd, 0, sizeof rd);
	rd.reference = r;
	rd.target_length = (qaws_scalar)(rf / 3);
	rd.iterations = 6;
	st = qaws_mesh_remesh(mesh, &rd, &dense, NULL);
	if (st != QAWS_STATUS_OK)
		goto done;
	/* 2. the fields, with the reference's normals */
	nrm = (double*)malloc(sizeof(double) * 3 * (dense->hm.nv + 1));
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
		double* init = (double*)calloc(3 * (size_t)dense->hm.nv + 3, sizeof(double));
		double *keep = NULL, *cand = NULL;
		unsigned int pulled = 0, s0, l, t;
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
		}
		mf_orientations(&h, desc->sweeps ? desc->sweeps : 8, init, 0);
		for (total = 0, l = 0; l < h.nl; l++)
			total += 3 * (size_t)h.lv[l].n;
		keep = (double*)malloc(sizeof(double) * (total + 1));
		cand = (double*)malloc(sizeof(double) * (3 * (size_t)dense->hm.nv + 3));
		if (!keep || !cand)
		{
			free(init);
			free(keep);
			free(cand);
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
				if (!pulled)
					continue;
				mf_orientations(&h, desc->sweeps ? desc->sweeps : 8, init, 1);
			}
			else
			{
				if (t == 1 && cf_create(&dense->hm, nrm, &h.lv[0], 4, &sys) != QAWS_STATUS_OK)
					break;
				if (t == 2 && !pulled)
					break;
				if ((t == 1 ? cf_smoothest(sys, cand) : cf_aligned(sys, init, cand)) != QAWS_STATUS_OK)
					continue;
				if (cf_relax(sys, rf, 100, cand) != QAWS_STATUS_OK)
					continue;
				memcpy(h.lv[0].q, cand, sizeof(double) * 3 * dense->hm.nv);
			}
			s1 = qd_singular(&dense->hm, &h.lv[0]);
			better = t == 1 ? s1 <= s0 : s1 < s0;
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
		free(keep);
		free(cand);
		free(init);
	}
	mf_positions(&h, rf, desc->sweeps ? desc->sweeps : 8);
	/* 3. the steps between lattice points made to close round every face, the mesh read off them; all quads at the
	   target when that leaves no other face, else (a hole or a cone off its lattice point left) the same at twice the
	   target, every face split */
	T = (int*)malloc(sizeof(int) * (2 * (size_t)dense->hm.nc + 2));
	K = (unsigned char*)malloc((size_t)dense->hm.nc + 1);
	if (!T || !K)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	st = qd_lattice(&dense->hm, &h, rf, T, K, r, &s, &o);
	for (f = 0; st == QAWS_STATUS_OK && desc->pure && f < o.nf && o.fs[f] == 4; f++) {}
	if (st == QAWS_STATUS_OK && desc->pure && f < o.nf)
	{
		qd_out_free(&o);
		mf_positions(&h, 2 * rho, desc->sweeps ? desc->sweeps : 8);
		st = qd_lattice(&dense->hm, &h, 2 * rho, T, K, r, &s, &o);
		if (st == QAWS_STATUS_OK)
			st = qd_split(&o, r, &s);
	}
	if (st != QAWS_STATUS_OK)
		goto done;
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
		unsigned char* kind = (unsigned char*)malloc(out->hm.nv + 1);
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
		qd_relax(&out->hm, r, &s, kind, 30);
		free(kind);
	}
	/* 4. attributes and the report */
	fl = (qaws_mesh_location*)malloc(sizeof(qaws_mesh_location) * (out->hm.nf + 1));
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
	if (st == QAWS_STATUS_OK)
		qaws_mesh_3d_measure_quality(out, &rep.quality);
done:
	free(fl);
	free(nrm);
	free(T);
	free(K);
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
