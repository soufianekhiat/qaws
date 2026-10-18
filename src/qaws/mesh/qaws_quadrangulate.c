/*
 * Quadrangulation, method INSTANT (Instant Meshes): the input remeshed
 * densely on the reference, a cross field and a lattice position field on
 * it, the quad mesh read off the fields (coinciding lattice points merged,
 * neighbours one step apart joined, faces traced round each vertex's
 * neighbours in turn), optionally split into all quads, relaxed on the
 * reference; footpoints and attributes as remeshing.
 *
 * shortcut: the read-off keeps the position field's local inconsistencies
 * (some triangles and pentagons, extra irregular vertices); QuadriFlow's
 * global consistency step (min-cost flow on the integer offsets) removes
 * them, the GLOBAL method to come.
 */

#include "qaws_mesh_internal.h"
#include "../qaws_quadrangulate.h"
#include "../qaws_remesh.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

static double qd_dot(double const* a, double const* b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void qd_cross(double const* a, double const* b, double* o)
{
	o[0] = a[1] * b[2] - a[2] * b[1];
	o[1] = a[2] * b[0] - a[0] * b[2];
	o[2] = a[0] * b[1] - a[1] * b[0];
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

typedef struct qd_pair
{
	unsigned int a, b;
} qd_pair;

static int qd_pair_cmp(void const* x, void const* y)
{
	qd_pair const *p = (qd_pair const*)x, *q = (qd_pair const*)y;
	if (p->a != q->a)
		return p->a < q->a ? -1 : 1;
	return p->b < q->b ? -1 : p->b > q->b ? 1 : 0;
}

/* a neighbour seen round its vertex: the angle, in the vertex's tangent frame */
typedef struct qd_spoke
{
	double ang;
	unsigned int to;
} qd_spoke;

static int qd_spoke_cmp(void const* x, void const* y)
{
	double a = ((qd_spoke const*)x)->ang, b = ((qd_spoke const*)y)->ang;
	return a < b ? -1 : a > b ? 1 : 0;
}

typedef struct qd_near
{
	unsigned int* uf;
	unsigned int const* member;
	double const *cx, *cn;
	double h2;
} qd_near;

static qaws_status qd_near_visit(void* user, unsigned int i, unsigned int j)
{
	qd_near const* w = (qd_near const*)user;
	double const *a = w->cx + 3 * (size_t)i, *b = w->cx + 3 * (size_t)j;
	double d[3];
	d[0] = a[0] - b[0];
	d[1] = a[1] - b[1];
	d[2] = a[2] - b[2];
	if (qd_dot(d, d) < w->h2 && qd_dot(w->cn + 3 * (size_t)i, w->cn + 3 * (size_t)j) > 0)
	{
		unsigned int ra = qd_find(w->uf, w->member[i]), rb = qd_find(w->uf, w->member[j]);
		if (ra != rb)
			w->uf[ra] = rb;
	}
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

static qaws_status qd_extract(qaws_hm const* dm, mf_hier const* h, double rho, qaws_mesh_reference const* r, mr_scratch* s, qd_out* out)
{
	mf_level const* L = &h->lv[0];
	unsigned int n = dm->nv, i, c, np = 0, ncl = 0, k, changed, pass;
	unsigned int *uf = NULL, *cl = NULL, *deg = NULL, *first = NULL, *adj = NULL, *member = NULL, *used = NULL, *face = NULL;
	double *cx = NULL, *cn = NULL, *cw = NULL;
	qd_pair* pairs = NULL;
	qd_spoke* spokes = NULL;
	unsigned char *ckind = NULL, *dirused = NULL;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	memset(out, 0, sizeof *out);
	uf = (unsigned int*)malloc(sizeof(unsigned int) * (n + 1));
	pairs = (qd_pair*)malloc(sizeof(qd_pair) * (dm->nc + 1));
	if (!uf || !pairs)
		goto done;
	for (i = 0; i < n; i++)
		uf[i] = i;
	/* each dense edge: the same lattice point (merged), one step apart (joined), or neither */
	for (c = 0; c < dm->nc; c++)
	{
		unsigned int a, b;
		double qa[3], qb[3], ob[3], d[3], t[3], dx, dy;
		if (dm->cf[c] == HM_NONE || (dm->co[c] != HM_NONE && dm->co[c] < c))
			continue;
		a = dm->cv[c];
		b = hm_to(dm, c);
		/* the two vertices' own lattice points, their offset in steps of a's lattice, b's unfolded into a's tangent plane (across a crease) */
		mf_compat_orientation(L->q + 3 * (size_t)a, L->nrm + 3 * (size_t)a, L->q + 3 * (size_t)b, L->nrm + 3 * (size_t)b, qa, qb);
		mf_unfold(L->x + 3 * (size_t)a, L->nrm + 3 * (size_t)a, L->x + 3 * (size_t)b, L->nrm + 3 * (size_t)b, L->o + 3 * (size_t)b, ob);
		for (k = 0; k < 3; k++)
			d[k] = ob[k] - L->o[3 * (size_t)a + k];
		qd_cross(L->nrm + 3 * (size_t)a, qa, t);
		dx = floor(qd_dot(d, qa) / rho + 0.5);
		dy = floor(qd_dot(d, t) / rho + 0.5);
		if (dx == 0 && dy == 0)
		{
			unsigned int ra = qd_find(uf, a), rb = qd_find(uf, b);
			if (ra != rb)
				uf[ra] = rb;
		}
		else if (fabs(dx) + fabs(dy) == 1)
		{
			pairs[np].a = a;
			pairs[np].b = b;
			np++;
		}
	}
	/* the clusters: their lattice points averaged; then clusters closer than a third of a step, alike in normal, merged
	   (one lattice point the dense edges did not tie: a vertex near a cell's centre rounding to a corner no neighbour did) */
	cl = (unsigned int*)malloc(sizeof(unsigned int) * (n + 1));
	member = (unsigned int*)malloc(sizeof(unsigned int) * (n + 1));
	cx = (double*)malloc(sizeof(double) * (3 * (size_t)n + 3));
	cn = (double*)malloc(sizeof(double) * (3 * (size_t)n + 3));
	cw = (double*)malloc(sizeof(double) * ((size_t)n + 1));
	ckind = (unsigned char*)malloc((size_t)n + 1);
	if (!cl || !member || !cx || !cn || !cw || !ckind)
		goto done;
	for (pass = 0; pass < 2; pass++)
	{
		ncl = 0;
		for (i = 0; i < n; i++)
			cl[i] = HM_NONE;
		for (i = 0; i < n; i++)
		{
			unsigned int root = qd_find(uf, i);
			if (cl[root] == HM_NONE)
			{
				member[ncl] = i;
				cl[root] = ncl++;
			}
			cl[i] = cl[root];
		}
		memset(cx, 0, sizeof(double) * 3 * ncl);
		memset(cn, 0, sizeof(double) * 3 * ncl);
		memset(cw, 0, sizeof(double) * ncl);
		memset(ckind, 0, ncl);
		for (i = 0; i < n; i++)
		{
			unsigned int q = cl[i];
			for (k = 0; k < 3; k++)
			{
				cx[3 * (size_t)q + k] += L->o[3 * (size_t)i + k];
				cn[3 * (size_t)q + k] += L->nrm[3 * (size_t)i + k];
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
		{
			double l = sqrt(qd_dot(cn + 3 * (size_t)i, cn + 3 * (size_t)i));
			for (k = 0; k < 3; k++)
			{
				cx[3 * (size_t)i + k] /= cw[i];
				cn[3 * (size_t)i + k] /= l > 0 ? l : 1;
			}
		}
		if (pass == 0)
		{
			qd_near w;
			qaws_bp_box* box = (qaws_bp_box*)malloc(sizeof(qaws_bp_box) * (ncl + 1));
			if (!box)
				goto done;
			for (i = 0; i < ncl; i++)
				for (k = 0; k < 3; k++)
				{
					box[i].lo[k] = cx[3 * (size_t)i + k] - 0.15 * rho;
					box[i].hi[k] = cx[3 * (size_t)i + k] + 0.15 * rho;
				}
			w.uf = uf;
			w.member = member;
			w.cx = cx;
			w.cn = cn;
			w.h2 = 0.09 * rho * rho;
			k = qaws_internal_broadphase(box, ncl, 3, NULL, qd_near_visit, &w, NULL) == QAWS_STATUS_OK;
			free(box);
			if (!k)
				goto done;
		}
	}
	/* the clusters' edges, once each */
	for (k = 0, i = 0; i < np; i++)
	{
		unsigned int a = cl[pairs[i].a], b = cl[pairs[i].b];
		if (a == b)
			continue;
		pairs[k].a = a < b ? a : b;
		pairs[k].b = a < b ? b : a;
		k++;
	}
	np = k;
	qsort(pairs, np, sizeof(qd_pair), qd_pair_cmp);
	for (k = 0, i = 0; i < np; i++)
		if (i == 0 || pairs[i].a != pairs[k - 1].a || pairs[i].b != pairs[k - 1].b)
			pairs[k++] = pairs[i];
	np = k;
	/* vertices on fewer than two edges go, until none */
	deg = (unsigned int*)calloc((size_t)ncl + 1, sizeof(unsigned int));
	if (!deg)
		goto done;
	do
	{
		changed = 0;
		memset(deg, 0, sizeof(unsigned int) * (ncl + 1));
		for (i = 0; i < np; i++)
		{
			deg[pairs[i].a]++;
			deg[pairs[i].b]++;
		}
		for (k = 0, i = 0; i < np; i++)
			if (deg[pairs[i].a] >= 2 && deg[pairs[i].b] >= 2)
				pairs[k++] = pairs[i];
			else
				changed = 1;
		np = k;
	} while (changed);
	/* each vertex's neighbours in turn round its normal */
	first = (unsigned int*)calloc((size_t)ncl + 2, sizeof(unsigned int));
	adj = (unsigned int*)malloc(sizeof(unsigned int) * (2 * (size_t)np + 1));
	spokes = (qd_spoke*)malloc(sizeof(qd_spoke) * (2 * (size_t)np + 1));
	dirused = (unsigned char*)calloc(2 * (size_t)np + 1, 1);
	if (!first || !adj || !spokes || !dirused)
		goto done;
	for (i = 0; i < np; i++)
	{
		first[pairs[i].a + 1]++;
		first[pairs[i].b + 1]++;
	}
	for (i = 0; i < ncl; i++)
		first[i + 1] += first[i];
	{
		unsigned int* fill = (unsigned int*)malloc(sizeof(unsigned int) * (ncl + 1));
		if (!fill)
			goto done;
		memcpy(fill, first, sizeof(unsigned int) * (ncl + 1));
		for (i = 0; i < np; i++)
		{
			adj[fill[pairs[i].a]++] = pairs[i].b;
			adj[fill[pairs[i].b]++] = pairs[i].a;
		}
		free(fill);
	}
	for (i = 0; i < ncl; i++)
	{
		double const *p = cx + 3 * (size_t)i, *nn = cn + 3 * (size_t)i;
		double u[3], w[3], a[3] = { 1, 0, 0 };
		unsigned int j;
		if (fabs(nn[0]) > 0.9)
		{
			a[0] = 0;
			a[1] = 1;
		}
		qd_cross(nn, a, u);
		{
			double l = sqrt(qd_dot(u, u));
			for (k = 0; k < 3; k++)
				u[k] /= l;
		}
		qd_cross(nn, u, w);
		for (j = first[i]; j < first[i + 1]; j++)
		{
			double d[3];
			for (k = 0; k < 3; k++)
				d[k] = cx[3 * (size_t)adj[j] + k] - p[k];
			spokes[j].ang = atan2(qd_dot(d, w), qd_dot(d, u));
			spokes[j].to = adj[j];
		}
		qsort(spokes + first[i], first[i + 1] - first[i], sizeof(qd_spoke), qd_spoke_cmp);
		for (j = first[i]; j < first[i + 1]; j++)
			adj[j] = spokes[j].to;
	}
	/* faces: from a to b, then the neighbour of b just before a counter-clockwise */
	out->fv = (unsigned int*)malloc(sizeof(unsigned int) * (2 * (size_t)np + 1));
	out->fs = (unsigned int*)malloc(sizeof(unsigned int) * (2 * (size_t)np + 1));
	face = (unsigned int*)malloc(sizeof(unsigned int) * 8);
	if (!out->fv || !out->fs || !face)
		goto done;
	for (i = 0; i < ncl; i++)
	{
		unsigned int e;
		for (e = first[i]; e < first[i + 1]; e++)
		{
			unsigned int a = i, b = adj[e], len = 0, ok = 1, cur = e, guard;
			if (dirused[e])
				continue;
			for (guard = 0; guard < 8; guard++)
			{
				unsigned int j, at = HM_NONE, nx;
				if (len == 8)
				{
					ok = 0;
					break;
				}
				face[len++] = a;
				dirused[cur] = 1;
				/* at b: where a is in its list, then the one before */
				for (j = first[b]; j < first[b + 1]; j++)
					if (adj[j] == a)
						at = j;
				if (at == HM_NONE)
				{
					ok = 0;
					break;
				}
				nx = at == first[b] ? first[b + 1] - 1 : at - 1;
				a = b;
				b = adj[nx];
				cur = nx;
				if (a == i && cur == e)
					break;
				if (dirused[cur])
				{
					ok = a == i && cur == e;
					break;
				}
			}
			if (!ok || guard == 8 || len < 3 || len > 6)
				continue;
			memcpy(out->fv + out->nc, face, sizeof(unsigned int) * len);
			out->fs[out->nf++] = len;
			out->nc += len;
		}
	}
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
	free(uf); free(cl); free(deg); free(first); free(adj); free(member); free(used); free(face);
	free(cx); free(cn); free(cw); free(pairs); free(spokes); free(ckind); free(dirused);
	if (st != QAWS_STATUS_OK)
		qd_out_free(out);
	return st;
}

/* small holes the tracing left (a loop of at most 8 border edges) closed by a face each */
static qaws_status qd_fill(qd_out* o)
{
	unsigned int f, k, c0, mask, *h, *start, nadd = 0, added = 0, *nf, *nfs;
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
	if (!nf || !nfs)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	memcpy(nf, o->fv, sizeof(unsigned int) * o->nc);
	memcpy(nfs, o->fs, sizeof(unsigned int) * o->nf);
	for (k = 0; k < o->nv; k++)
	{
		unsigned int loop[9], len = 0, v = k, ok = 0, j;
		if (start[k] >= HM_NONE - 1)
			continue;
		while (len < 9)
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
		if (!ok || len < 3 || len > 8)
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

/* a few steps of each free vertex toward its neighbours' centre, along the reference */
static void qd_relax(qaws_hm* m, qaws_mesh_reference const* r, unsigned char const* kind, unsigned int steps)
{
	unsigned int it, v, k, ring[64];
	for (it = 0; it < steps; it++)
		for (v = 0; v < m->nv; v++)
		{
			unsigned int n, i, c;
			double cen[3] = { 0, 0, 0 }, d[3], x[3], nrm[3], dn;
			qaws_mesh_location l;
			int bad = 0;
			if (kind[v] || hm_is_border_vertex(m, v) || m->vc[v] == HM_NONE)
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
	rf = desc->pure ? 2 * rho : rho;
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
		/* the start: the direction of largest curvature where the surface bends one way more than the other */
		double* init = (double*)calloc(3 * (size_t)dense->hm.nv + 3, sizeof(double));
		if (!init)
		{
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			goto done;
		}
		for (v = 0; v < dense->hm.nv; v++)
		{
			qaws_mesh_location const* l = &dense->hm.foot[v];
			double xa[3], xb[3], n[3], na[3], nb[3], k1, k2, d1[3], d2[3], big;
			mi_eval(&r->patch[l->patch], (double)l->a, (double)l->b, NULL, xa, xb, NULL, NULL, NULL, n, na, nb);
			mi_curvature(xa, xb, n, na, nb, &k1, &k2, d1, d2);
			big = fmax(fabs(k1), fabs(k2));
			if (big > 0 && k1 - k2 > 0.15 * big)
				memcpy(init + 3 * (size_t)v, d1, sizeof d1);
		}
		mf_orientations(&h, desc->sweeps ? desc->sweeps : 8, init);
		free(init);
	}
	mf_positions(&h, rf, desc->sweeps ? desc->sweeps : 8);
	/* 3. the mesh read off them */
	st = qd_extract(&dense->hm, &h, rf, r, &s, &o);
	if (st == QAWS_STATUS_OK)
		st = qd_fill(&o);
	if (st == QAWS_STATUS_OK && desc->pure)
		st = qd_split(&o, r, &s);
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
		qd_relax(&out->hm, r, kind, 3);
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
