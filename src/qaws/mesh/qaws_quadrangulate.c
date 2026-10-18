/*
 * Quadrangulation, method INSTANT (Instant Meshes with QuadriFlow's
 * consistency): the input remeshed densely on the reference; a cross field
 * (twice: free, and drawn to the principal curvature directions; the one
 * with fewer singularities kept) and a lattice position field on it; each
 * dense edge's step between lattice points rounded, then mended until the
 * steps close round every face (dislocations cancelled in pairs along the
 * cheapest paths); the quad mesh read off the steps (lattice points merged,
 * neighbours one step apart joined, faces traced round each vertex's
 * neighbours in turn, doublets removed), optionally split into all quads,
 * relaxed on the reference (creases and the border along themselves);
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

/* the lattice coordinates c turned by k quarter turns counter-clockwise */
static void qd_turn(int const* c, int k, int* o)
{
	int x = c[0], y = c[1];
	switch (k & 3)
	{
	case 0: o[0] = x; o[1] = y; break;
	case 1: o[0] = -y; o[1] = x; break;
	case 2: o[0] = -x; o[1] = -y; break;
	default: o[0] = y; o[1] = -x; break;
	}
}

/* Per corner of the dense mesh (the edge from its vertex a to the next, b): T, the step from a's lattice point to b's in
   whole steps of a's frame (q, n x q); K, the quarter turns taking b's frame onto a's (coordinates in b's frame are those
   in a's turned by K). The two corners of an edge are kept opposite. */
static void qd_steps(qaws_hm const* dm, mf_level const* L, double rho, int* T, unsigned char* K)
{
	unsigned int c, k;
	for (c = 0; c < dm->nc; c++)
	{
		unsigned int a, b, o = dm->co[c];
		double ob[3], d[3], t[3];
		if (dm->cf[c] == HM_NONE || (o != HM_NONE && o < c))
			continue;
		a = dm->cv[c];
		b = hm_to(dm, c);
		K[c] = (unsigned char)mf_rotation(L->q + 3 * (size_t)a, L->nrm + 3 * (size_t)a, L->q + 3 * (size_t)b, L->nrm + 3 * (size_t)b);
		mf_unfold(L->x + 3 * (size_t)a, L->nrm + 3 * (size_t)a, L->x + 3 * (size_t)b, L->nrm + 3 * (size_t)b, L->o + 3 * (size_t)b, ob);
		for (k = 0; k < 3; k++)
			d[k] = ob[k] - L->o[3 * (size_t)a + k];
		qd_cross(L->nrm + 3 * (size_t)a, L->q + 3 * (size_t)a, t);
		T[2 * (size_t)c] = (int)floor(qd_dot(d, L->q + 3 * (size_t)a) / rho + 0.5);
		T[2 * (size_t)c + 1] = (int)floor(qd_dot(d, t) / rho + 0.5);
		if (o != HM_NONE)
		{
			int m[2];
			m[0] = -T[2 * (size_t)c];
			m[1] = -T[2 * (size_t)c + 1];
			K[o] = (unsigned char)((4 - K[c]) & 3);
			qd_turn(m, K[c], T + 2 * (size_t)o);
		}
	}
}

/* face f's steps summed in its first vertex's frame (zero: the lattice closes round it); fr per corner: the quarter turns
   from its vertex's frame to the face's. Nonzero return: the frames turn round f (a singularity of the cross field). */
static int qd_closure(qaws_hm const* dm, unsigned int f, int const* T, unsigned char const* K, unsigned char* fr, int* S)
{
	unsigned int c0 = dm->fc[f], i;
	int turn = 0, s[2];
	S[0] = S[1] = 0;
	for (i = 0; i < 3; i++)
	{
		fr[c0 + i] = (unsigned char)(turn & 3);
		qd_turn(T + 2 * (size_t)(c0 + i), turn, s);
		S[0] += s[0];
		S[1] += s[1];
		turn -= K[c0 + i];
	}
	return (turn & 3) != 0;
}

/* the lattice axis (0: q, 1: n x q of its vertex) corner c's edge runs along */
static int qd_axis(qaws_hm const* dm, mf_level const* L, unsigned int c)
{
	unsigned int a = dm->cv[c], b = hm_to(dm, c), k;
	double d[3], t[3];
	for (k = 0; k < 3; k++)
		d[k] = L->x[3 * (size_t)b + k] - L->x[3 * (size_t)a + k];
	qd_cross(L->nrm + 3 * (size_t)a, L->q + 3 * (size_t)a, t);
	return fabs(qd_dot(d, L->q + 3 * (size_t)a)) >= fabs(qd_dot(d, t)) ? 0 : 1;
}

/* a binary heap of (cost, node), cheapest on top */
typedef struct qd_item
{
	double d;
	unsigned int n;
} qd_item;

static void qd_push(qd_item* h, unsigned int* n, double d, unsigned int v)
{
	unsigned int i = (*n)++;
	while (i > 0 && h[(i - 1) / 2].d > d)
	{
		h[i] = h[(i - 1) / 2];
		i = (i - 1) / 2;
	}
	h[i].d = d;
	h[i].n = v;
}

static qd_item qd_pop(qd_item* h, unsigned int* n)
{
	qd_item top = h[0], last = h[--(*n)];
	unsigned int i = 0, j;
	while ((j = 2 * i + 1) < *n)
	{
		if (j + 1 < *n && h[j + 1].d < h[j].d)
			j++;
		if (h[j].d >= last.d)
			break;
		h[i] = h[j];
		i = j;
	}
	h[i] = last;
	return top;
}

/* the cost of an edge's step t becoming n: lattice points merged cheapest, a lattice line between the two next, a step
   turned or made diagonal (the cells shear) dear, a step of two (a lattice point with no vertex) a last resort */
static double qd_cost(int const* t, int const* n)
{
	int a0 = abs(t[0]) + abs(t[1]), a1 = abs(n[0]) + abs(n[1]);
	if (abs(n[0]) > 1 || abs(n[1]) > 1)
		return 50;
	if (a1 < a0)
		return 1;
	return a0 == 0 ? 2 : 6;
}

/* corner c's edge takes the unit step u (in the frame of the path's first face; tr: its turn into c's face): -u off c's
   face's sum, +u onto the face across */
static void qd_step(qaws_hm const* dm, int* T, unsigned char const* K, unsigned char const* fr, unsigned char const* tr, int const* u, unsigned int c)
{
	unsigned int o = dm->co[c];
	int uP[2], dd[2], m[2];
	qd_turn(u, tr[dm->cf[c]], uP);
	qd_turn(uP, -(int)fr[c], dd);
	T[2 * (size_t)c] -= dd[0];
	T[2 * (size_t)c + 1] -= dd[1];
	if (o != HM_NONE)
	{
		qd_turn(dd, K[c], m);
		T[2 * (size_t)o] += m[0];
		T[2 * (size_t)o + 1] += m[1];
	}
}

/* QuadriFlow's consistency (Huang, Zhou, Funkhouser, Niessner 2018): a face whose steps do not close holds a dislocation
   (an extra half row of the lattice ends there). Dislocations cancelled in pairs, or at the border, along shortest paths
   across the faces: each edge crossed takes the unit step, the faces between keep their sums. Feature and border edges,
   lattice lines, take steps along themselves only. Returns the dislocations left (round singular faces, or alone).
   shortcut: nearest partners in turn and paths around singular faces stand in for its min-cost flow; the flow when
   the extra rows should go where the stretch is least. */
static unsigned int qd_consistent(qaws_hm const* dm, mf_level const* L, int* T, unsigned char const* K)
{
	unsigned int nf = dm->nf, f0, f, run = 0, left = 0;
	unsigned int* from = (unsigned int*)malloc(sizeof(unsigned int) * (nf + 1));
	unsigned int* stamp = (unsigned int*)calloc((size_t)nf + 1, sizeof(unsigned int));
	unsigned char *fr = (unsigned char*)malloc((size_t)dm->nc + 1), *tr = (unsigned char*)malloc((size_t)nf + 1), *sing = (unsigned char*)malloc((size_t)nf + 1);
	int* S = (int*)malloc(sizeof(int) * (2 * (size_t)nf + 2));
	double* dist = (double*)malloc(sizeof(double) * (nf + 1));
	qd_item* heap = (qd_item*)malloc(sizeof(qd_item) * (3 * (size_t)nf + dm->nc + 8));   /* (each face pushed at most thrice, each exit once) */
	if (!from || !stamp || !fr || !tr || !sing || !S || !dist || !heap)
	{
		left = HM_NONE;
		goto done;
	}
	for (f = 0; f < nf; f++)
		sing[f] = (unsigned char)qd_closure(dm, f, T, K, fr, S + 2 * (size_t)f);
	for (f0 = 0; f0 < nf; f0++)
		while (!sing[f0] && (S[2 * (size_t)f0] || S[2 * (size_t)f0 + 1]))
		{
			int u[2], ax = abs(S[2 * (size_t)f0]) >= abs(S[2 * (size_t)f0 + 1]) ? 0 : 1;
			unsigned int hn = 0, hit = HM_NONE, exitc = HM_NONE, B;
			u[0] = u[1] = 0;
			u[ax] = S[2 * (size_t)f0 + ax] > 0 ? 1 : -1;
			stamp[f0] = ++run;
			dist[f0] = 0;
			tr[f0] = 0;
			from[f0] = HM_NONE;
			qd_push(heap, &hn, 0, f0);
			/* the cheapest path from f0 to a face holding the opposite step, or out through the border */
			while (hn && hit == HM_NONE)
			{
				qd_item it = qd_pop(heap, &hn);
				unsigned int A = it.n, c;
				int uA[2];
				if (A >= nf)
				{
					/* out through the border (along it: checked when pushed) */
					exitc = A - nf;
					hit = dm->cf[exitc];
					break;
				}
				if (it.d > dist[A])
					continue;
				qd_turn(u, tr[A], uA);
				if (A != f0 && S[2 * (size_t)A] * uA[0] + S[2 * (size_t)A + 1] * uA[1] < 0)
				{
					hit = A;
					break;
				}
				for (c = dm->fc[A]; c < dm->fc[A] + 3; c++)
				{
					unsigned int o = dm->co[c];
					int ua[2], nt[2];
					double w;
					qd_turn(uA, -(int)fr[c], ua);
					if ((o == HM_NONE || ((dm->cflag[c] | dm->cflag[o]) & HM_FEATURE)) && ua[qd_axis(dm, L, c) ^ 1] != 0)
						continue;
					nt[0] = T[2 * (size_t)c] - ua[0];
					nt[1] = T[2 * (size_t)c + 1] - ua[1];
					w = it.d + qd_cost(T + 2 * (size_t)c, nt);
					if (o == HM_NONE)
					{
						qd_push(heap, &hn, w, nf + c);
						continue;
					}
					B = dm->cf[o];
					if (sing[B] || (stamp[B] == run && dist[B] <= w))
						continue;
					stamp[B] = run;
					dist[B] = w;
					tr[B] = (unsigned char)((tr[A] - fr[c] + K[c] + fr[o]) & 3);
					from[B] = c;
					qd_push(heap, &hn, w, B);
				}
			}
			if (hit == HM_NONE)
				break;   /* (no partner: left) */
			/* each edge on the way takes the step: -u off the face before it, +u onto the one after (the border: none) */
			if (exitc != HM_NONE)
				qd_step(dm, T, K, fr, tr, u, exitc);
			for (B = hit; from[B] != HM_NONE; B = dm->cf[from[B]])
				qd_step(dm, T, K, fr, tr, u, from[B]);
			qd_closure(dm, f0, T, K, fr, S + 2 * (size_t)f0);
			qd_closure(dm, hit, T, K, fr, S + 2 * (size_t)hit);
		}
	for (f = 0; f < nf; f++)
		left += !sing[f] && (S[2 * (size_t)f] || S[2 * (size_t)f + 1]);
done:
	free(from); free(stamp); free(fr); free(tr); free(sing); free(S); free(dist); free(heap);
	return left;
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

static qaws_status qd_extract(qaws_hm const* dm, mf_hier const* h, double rho, int const* T, qaws_mesh_reference const* r, mr_scratch* s,
	qd_out* out)
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
		int dx, dy;
		if (dm->cf[c] == HM_NONE || (dm->co[c] != HM_NONE && dm->co[c] < c))
			continue;
		a = dm->cv[c];
		b = hm_to(dm, c);
		dx = T[2 * (size_t)c];
		dy = T[2 * (size_t)c + 1];
		if (dx == 0 && dy == 0)
		{
			unsigned int ra = qd_find(uf, a), rb = qd_find(uf, b);
			if (ra != rb)
				uf[ra] = rb;
		}
		else if (abs(dx) + abs(dy) == 1)
		{
			pairs[np].a = a;
			pairs[np].b = b;
			np++;
		}
	}
	/* the clusters: their members' positions averaged (the mended steps may join members of different lattice points);
	   then clusters closer than a third of a step, alike in normal, merged
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
				cx[3 * (size_t)q + k] += L->x[3 * (size_t)i + k];
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
		/* the cross field twice from the direction of largest curvature: free (Instant Meshes, along the features), and drawn
		   to it on every level, as strongly as the surface bends one way more than the other (1 round a cylinder
		   or a torus, none where it bends alike or hardly at the lattice's scale); the one with fewer singularities kept
		   (a closed tube wants the pull, a patch held by its border not) */
		double* init = (double*)calloc(3 * (size_t)dense->hm.nv + 3, sizeof(double));
		double* keep = NULL;
		unsigned int pulled = 0, s0, l;
		size_t at;
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
		if (pulled)
		{
			for (at = 0, l = 0; l < h.nl; l++)
				at += 3 * (size_t)h.lv[l].n;
			keep = (double*)malloc(sizeof(double) * (at + 1));
			if (!keep)
			{
				free(init);
				st = QAWS_STATUS_ALLOCATION_FAILURE;
				goto done;
			}
			for (at = 0, l = 0; l < h.nl; at += 3 * (size_t)h.lv[l].n, l++)
				memcpy(keep + at, h.lv[l].q, sizeof(double) * 3 * h.lv[l].n);
			s0 = qd_singular(&dense->hm, &h.lv[0]);
			mf_orientations(&h, desc->sweeps ? desc->sweeps : 8, init, 1);
			if (qd_singular(&dense->hm, &h.lv[0]) >= s0)
				for (at = 0, l = 0; l < h.nl; at += 3 * (size_t)h.lv[l].n, l++)
					memcpy(h.lv[l].q, keep + at, sizeof(double) * 3 * h.lv[l].n);
			free(keep);
		}
		free(init);
	}
	mf_positions(&h, rf, desc->sweeps ? desc->sweeps : 8);
	/* 3. the steps between lattice points made to close round every face, the mesh read off them */
	T = (int*)malloc(sizeof(int) * (2 * (size_t)dense->hm.nc + 2));
	K = (unsigned char*)malloc((size_t)dense->hm.nc + 1);
	if (!T || !K)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	qd_steps(&dense->hm, &h.lv[0], rf, T, K);
	if (qd_consistent(&dense->hm, &h.lv[0], T, K) == HM_NONE)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	st = qd_extract(&dense->hm, &h, rf, T, r, &s, &o);
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
