/*
 * Local remeshing on a reference (Botsch and Kobbelt 2004, with the
 * refinements of pmp, geometry-central and Mmgs): long edges split at their
 * midpoint on the reference, short ones collapsed onto theirs, flips toward
 * the valence each vertex's sectors call for, vertices relaxed to the
 * size-weighted centroid of their ring by walking on the reference, needles
 * collapsed and caps flipped, Delaunay flips to finish; on demand an
 * angle-bounded pass (Wang, Yan et al. 2019) and a distance guard (Mmgs's
 * hausd). Feature corners fixed, feature lines and the border slid along;
 * footpoints kept throughout; the reference mesh's attributes transferred
 * at the end.
 */

#include "qaws_mesh_internal.h"
#include "../internal/qaws_internal_memory.h"
#include "../qaws_remesh.h"
#include "qaws_remesh_cvt.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define RM_PI 3.14159265358979323846
#define RM_DEG (RM_PI / 180)

typedef struct rmsh
{
	qaws_hm* m;                     /* the working mesh (the output's store) */
	qaws_mesh_reference const* r;
	qaws_size_field const* size;
	double L;                       /* the uniform target */
	unsigned int sa;                /* the size attribute (per vertex) */
	unsigned int ra;                /* the region attribute (per face): the reference's triangles between creases the face is on */
	mr_scratch s;
	qaws_mesh_report rep;
	unsigned int* ep;               /* the edges as vertex pairs at a pass's start (the Delaunay pass's queue) */
	unsigned int nep, cep;
	double tol;                     /* the distance guard (0: none) */
	qaws_size_field* tsize;         /* with a tolerance: the sizes its chord gap allows where the surface bends (under the target) */
	int strict;                     /* after the loop: the guard samples quarter points too, and holds where faces strayed already */
	double lo, hi;                  /* the angle bounds in radians (0: none), set for the angle pass */
	int bounded;                    /* operations held to the bounds (in the angle pass) */
	double moved;                   /* the longest move since it was cleared */
	unsigned int* chg;              /* the vertices the angle pass changed */
	unsigned int nchg, cchg;
	qaws_metric_field const* metric;   /* anisotropic: lengths and angles in it (NULL: none) */
	unsigned int ma;                /* with a metric: per vertex, its metric in its tangent plane (six) */
	double M[6];                    /* with a metric: the one the current operation measures in (rm_use) */
	int ready, failed;
} rmsh;

static double rm_d2(double const* a, double const* b)
{
	return (a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]);
}

static double* rm_x(qaws_hm const* m, unsigned int v)
{
	return m->x + 3 * (size_t)v;
}

static double rm_size_of(rmsh const* w, unsigned int v)
{
	return hm_attr_at(w->m, w->sa, v)[0];
}

/* u^T M t (M NULL: u . t) */
static double rm_mdot(double const* M, double const* u, double const* t)
{
	if (!M)
		return u[0] * t[0] + u[1] * t[1] + u[2] * t[2];
	return u[0] * (M[0] * t[0] + M[1] * t[1] + M[2] * t[2]) + u[1] * (M[1] * t[0] + M[3] * t[1] + M[4] * t[2])
		+ u[2] * (M[2] * t[0] + M[4] * t[1] + M[5] * t[2]);
}

/* the metric the current operation measures in, NULL without one */
static double const* rm_M(rmsh const* w)
{
	return w->metric ? w->M : NULL;
}

/* the current operation's metric: the mean of a's and b's */
static void rm_use(rmsh* w, unsigned int a, unsigned int b)
{
	unsigned int k;
	for (k = 0; w->metric && k < 6; k++)
		w->M[k] = 0.5 * (hm_attr_at(w->m, w->ma, a)[k] + hm_attr_at(w->m, w->ma, b)[k]);
}

/* the target length between vertices a and b; with a metric, along direction e: its length over its length in the mean of
   a's and b's metrics */
static double rm_target_along(rmsh const* w, unsigned int a, unsigned int b, double const* e)
{
	if (w->metric)
	{
		double M[6], q;
		unsigned int k;
		for (k = 0; k < 6; k++)
			M[k] = 0.5 * (hm_attr_at(w->m, w->ma, a)[k] + hm_attr_at(w->m, w->ma, b)[k]);
		q = rm_mdot(M, e, e);
		if (q > 0)
			return sqrt(rm_mdot(NULL, e, e) / q);
	}
	return 0.5 * (rm_size_of(w, a) + rm_size_of(w, b));
}

/* the target length of edge (a, b) */
static double rm_target(rmsh const* w, unsigned int a, unsigned int b)
{
	double e[3];
	unsigned int k;
	for (k = 0; k < 3; k++)
		e[k] = rm_x(w->m, b)[k] - rm_x(w->m, a)[k];
	return rm_target_along(w, a, b, e);
}

/* the size at v from its footpoint; with a metric, also v's metric in its tangent plane (capped by the sizes when a size
   field or tolerance gives them), the size its mean (so the relaxation's weights are metric areas) */
static void rm_set_size(rmsh* w, unsigned int v)
{
	double h = w->size ? sf_eval(w->size, &w->m->foot[v], rm_x(w->m, v)) : w->L;
	if (w->tsize)
		h = fmin(h, sf_eval(w->tsize, &w->m->foot[v], rm_x(w->m, v)));
	if (w->metric)
	{
		double M[6], x[3], n[3], t1[3], t2[3], g[3], ln;
		mf_eval(w->metric, &w->m->foot[v], M);
		mr_eval(w->r, &w->m->foot[v], x, n);
		ln = sqrt(rm_mdot(NULL, n, n));
		if (ln > 0)
		{
			n[0] /= ln; n[1] /= ln; n[2] /= ln;
			mf_basis(n, t1, t2);
			g[0] = rm_mdot(M, t1, t1);
			g[1] = rm_mdot(M, t1, t2);
			g[2] = rm_mdot(M, t2, t2);
			h = pow(mf_lift(hm_attr_at(w->m, w->ma, v), g, t1, t2, n, w->size || w->tsize ? 1 / (h * h) : 0, HUGE_VAL), -0.25);
		}
		else
			memcpy(hm_attr_at(w->m, w->ma, v), M, sizeof M);
	}
	hm_attr_at(w->m, w->sa, v)[0] = h;
}

/* a triangle's shape: 4 sqrt(3) area over the sum of its squared edges (1 equilateral, 0 flat); in metric M (NULL: none) */
static double rm_quality(double const* M, double const* A, double const* B, double const* C)
{
	double u[3], t[3], n[3], s = 0;
	unsigned int k;
	if (M)
	{
		double uu, tt, ut, r[3];
		for (k = 0; k < 3; k++)
		{
			u[k] = B[k] - A[k];
			t[k] = C[k] - A[k];
			r[k] = C[k] - B[k];
		}
		uu = rm_mdot(M, u, u);
		tt = rm_mdot(M, t, t);
		ut = rm_mdot(M, u, t);
		s = uu + tt + rm_mdot(M, r, r);
		return s > 0 ? 2 * sqrt(3.0) * sqrt(fmax(uu * tt - ut * ut, 0)) / s : 0;
	}
	for (k = 0; k < 3; k++)
	{
		u[k] = B[k] - A[k];
		t[k] = C[k] - A[k];
		s += u[k] * u[k] + t[k] * t[k] + (C[k] - B[k]) * (C[k] - B[k]);
	}
	n[0] = u[1] * t[2] - u[2] * t[1]; n[1] = u[2] * t[0] - u[0] * t[2]; n[2] = u[0] * t[1] - u[1] * t[0];
	return s > 0 ? 2 * sqrt(3.0) * sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]) / s : 0;
}

/* the cosine of the angle at P of triangle (P, A, B), in metric M (NULL: none) */
static double rm_cos(double const* M, double const* P, double const* A, double const* B)
{
	double u[3], t[3], uu, tt, ut;
	unsigned int k;
	for (k = 0; k < 3; k++)
	{
		u[k] = A[k] - P[k];
		t[k] = B[k] - P[k];
	}
	uu = rm_mdot(M, u, u);
	tt = rm_mdot(M, t, t);
	ut = rm_mdot(M, u, t);
	return uu > 0 && tt > 0 ? ut / sqrt(uu * tt) : 1;
}

/* the angle at P of triangle (P, A, B), in metric M (NULL: none) */
static double rm_angle(double const* M, double const* P, double const* A, double const* B)
{
	double c = rm_cos(M, P, A, B);
	return acos(c < -1 ? -1 : c > 1 ? 1 : c);
}

/* the angle at corner c of its triangle; with a metric, in the mean of its corners' */
static double rm_corner(rmsh const* w, unsigned int c)
{
	qaws_hm const* m = w->m;
	unsigned int p = m->cv[c], a = hm_to(m, c), b = m->cv[hm_prev(m, c)], k;
	double M[6];
	for (k = 0; w->metric && k < 6; k++)
		M[k] = (hm_attr_at(m, w->ma, p)[k] + hm_attr_at(m, w->ma, a)[k] + hm_attr_at(m, w->ma, b)[k]) / 3;
	return rm_angle(w->metric ? M : NULL, rm_x(m, p), rm_x(m, a), rm_x(m, b));
}

/* edge c's length; with a metric, in the current operation's */
static double rm_len(rmsh const* w, unsigned int c)
{
	double e[3];
	unsigned int k;
	if (!w->metric)
		return hm_length(w->m, c);
	for (k = 0; k < 3; k++)
		e[k] = rm_x(w->m, hm_to(w->m, c))[k] - rm_x(w->m, w->m->cv[c])[k];
	return sqrt(rm_mdot(w->M, e, e));
}

/* face f's longest edge */
static unsigned int rm_longest(qaws_hm const* m, unsigned int f)
{
	unsigned int c, e = m->fc[f];
	for (c = e + 1; c < m->fc[f] + 3; c++)
		if (hm_length(m, c) > hm_length(m, e))
			e = c;
	return e;
}

/* a triangle's smallest angle */
static double rm_tri_min(double const* M, double const* A, double const* B, double const* C)
{
	return fmin(rm_angle(M, A, B, C), fmin(rm_angle(M, B, C, A), rm_angle(M, C, A, B)));
}

#define RM_FAIR 0.3   /* below this a triangle is poor: no operation may leave the poorest one round it poorer */

/* would moving v to y turn over a face round v (faces holding `skip` aside), or
   leave the poorest face poorer than before and than fair? */
static int rm_folds(rmsh const* w, unsigned int v, double const* y, unsigned int skip)
{
	qaws_hm const* m = w->m;
	unsigned int c;
	double qb = 1, qa = 1;
	for (c = m->vc[v]; c != HM_NONE; c = hm_round(m, v, c))
	{
		unsigned int a = hm_to(m, c), b = m->cv[hm_prev(m, c)], k;
		double const *p = rm_x(m, v), *A = rm_x(m, a), *B = rm_x(m, b);
		double u[3], t[3], n0[3], n1[3], d0, d1, dd;
		if (a == skip || b == skip)
			continue;
		for (k = 0; k < 3; k++) { u[k] = A[k] - p[k]; t[k] = B[k] - p[k]; }
		n0[0] = u[1] * t[2] - u[2] * t[1]; n0[1] = u[2] * t[0] - u[0] * t[2]; n0[2] = u[0] * t[1] - u[1] * t[0];
		for (k = 0; k < 3; k++) { u[k] = A[k] - y[k]; t[k] = B[k] - y[k]; }
		n1[0] = u[1] * t[2] - u[2] * t[1]; n1[1] = u[2] * t[0] - u[0] * t[2]; n1[2] = u[0] * t[1] - u[1] * t[0];
		d0 = n0[0] * n0[0] + n0[1] * n0[1] + n0[2] * n0[2];
		d1 = n1[0] * n1[0] + n1[1] * n1[1] + n1[2] * n1[2];
		dd = n0[0] * n1[0] + n0[1] * n1[1] + n0[2] * n1[2];
		if (dd <= 0.2 * sqrt(d0 * d1) || d1 < 1e-12 * d0)
			return 1;   /* turned by more than about 78 degrees, or flattened */
		qb = fmin(qb, rm_quality(rm_M(w), p, A, B));
		qa = fmin(qa, rm_quality(rm_M(w), y, A, B));
	}
	return qa < qb && qa < RM_FAIR;
}

/* the region of face f */
static unsigned int rm_region(rmsh const* w, unsigned int f)
{
	return (unsigned int)hm_attr_at(w->m, w->ra, f)[0];
}

/* an end's footpoint triangle in the region as a hint (HM_NONE: none is) */
static unsigned int rm_hint(rmsh const* w, unsigned int a, unsigned int b, unsigned int region)
{
	unsigned int ta = w->m->foot[a].patch, tb = w->m->foot[b].patch;
	if (ta < w->r->nt && w->r->treg[ta] == region)
		return ta;
	return tb < w->r->nt && w->r->treg[tb] == region ? tb : HM_NONE;
}

/* the edges as vertex pairs, taken at a pass's start: an edit moves the corners of the faces it
   touches past the end of the store, so walking the corners would skip their edges */
static int rm_snapshot(rmsh* w)
{
	qaws_hm const* m = w->m;
	unsigned int c;
	if (m->nc > w->cep)
	{
		unsigned int* g = (unsigned int*)qaws_internal_realloc(w->ep, qaws_internal_mul3(sizeof(unsigned int), 2, m->nc));
		if (!g)
		{
			w->failed = 1;
			return 0;
		}
		w->ep = g;
		w->cep = m->nc;
	}
	w->nep = 0;
	for (c = 0; c < m->nc; c++)
		if (m->cf[c] != HM_NONE && (m->co[c] == HM_NONE || c < m->co[c]))
		{
			w->ep[2 * w->nep] = m->cv[c];
			w->ep[2 * w->nep + 1] = hm_to(m, c);
			w->nep++;
		}
	return 1;
}

/* the corner of snapshot edge i now, or HM_NONE when it is gone */
static unsigned int rm_edge(rmsh const* w, unsigned int i)
{
	unsigned int a = w->ep[2 * i], b = w->ep[2 * i + 1];
	if ((w->m->vflag[a] | w->m->vflag[b]) & HM_DELETED)
		return HM_NONE;
	return hm_find_edge(w->m, a, b);
}

/* an edge (a, b) pushed on the snapshot as on a stack */
static int rm_push(rmsh* w, unsigned int a, unsigned int b)
{
	if (w->nep == w->cep)
	{
		unsigned int* g = (unsigned int*)qaws_internal_realloc(w->ep, qaws_internal_mul3(sizeof(unsigned int), 2, (2 * (size_t)w->cep + 64)));
		if (!g)
		{
			w->failed = 1;
			return 0;
		}
		w->ep = g;
		w->cep = 2 * w->cep + 64;
	}
	w->ep[2 * w->nep] = a;
	w->ep[2 * w->nep + 1] = b;
	w->nep++;
	return 1;
}

/* a vertex the angle pass changed, for its relaxation */
static void rm_touch(rmsh* w, unsigned int v)
{
	if (!w->bounded)
		return;
	if (w->nchg == w->cchg)
	{
		unsigned int* g = (unsigned int*)qaws_internal_realloc(w->chg, qaws_internal_mul2(sizeof(unsigned int), (2 * (size_t)w->cchg + 64)));
		if (!g)
		{
			w->failed = 1;
			return;
		}
		w->chg = g;
		w->cchg = 2 * w->cchg + 64;
	}
	w->chg[w->nchg++] = v;
}

/* ---- angle bounds ------------------------------------------------------- */

/* the angle of the sector holding corner c: its vertex's corners between two feature (or border) edges, summed, their count
   in faces (may be NULL); 0 when the vertex is on no feature */
static double rm_sector(qaws_hm const* m, unsigned int c, unsigned int* faces)
{
	unsigned int k = c, n, f = 0;
	double th = 0;
	if (faces)
		*faces = 0;
	if (!hm_feature_degree(m, m->cv[c]))
		return 0;
	for (n = 0; n < 256 && !(m->cflag[k] & HM_FEATURE); n++)
	{
		unsigned int b = hm_cw(m, k);
		if (b == HM_NONE || b == c)
			break;
		k = b;
	}
	for (n = 0; n < 256; n++)
	{
		th += rm_angle(NULL, rm_x(m, m->cv[k]), rm_x(m, hm_to(m, k)), rm_x(m, m->cv[hm_prev(m, k)]));
		f++;
		if (m->cflag[hm_prev(m, k)] & HM_FEATURE)
			break;
		k = hm_ccw(m, k);
		if (k == HM_NONE)
			break;
	}
	if (faces)
		*faces = f;
	return th;
}

/* is corner c exempt from the bounds: in a sector of angle theta where no whole number k of corners has theta / k inside them */
static int rm_exempt(rmsh const* w, unsigned int c)
{
	double th, k;
	if (w->lo <= 0 || (th = rm_sector(w->m, c, NULL)) <= 0)
		return 0;
	k = floor(th / w->lo);
	return k < 1 || (w->hi > 0 && th / k > w->hi);
}

/* how far an angle is outside the bounds (negative: its margin inside) */
static double rm_viol(rmsh const* w, double a)
{
	double v = -RM_PI;
	if (w->lo > 0)
		v = w->lo - a;
	if (w->hi > 0 && a - w->hi > v)
		v = a - w->hi;
	return v;
}

typedef struct rm_stat
{
	double amin, amax;              /* over every corner */
	double worst;                   /* the largest violation of the bounds over the corners not exempt */
	unsigned int bad;               /* corners not exempt outside the bounds */
	double sq;                      /* the squared distances of every angle from 60 degrees */
	unsigned int n;
} rm_stat;

static void rm_stat_init(rm_stat* s)
{
	s->amin = RM_PI;
	s->amax = 0;
	s->worst = -RM_PI;
	s->bad = 0;
	s->sq = 0;
	s->n = 0;
}

static void rm_stat_add(rmsh const* w, rm_stat* s, double a, int exempt)
{
	double v;
	s->amin = fmin(s->amin, a);
	s->amax = fmax(s->amax, a);
	s->sq += (a - 60 * RM_DEG) * (a - 60 * RM_DEG);
	s->n++;
	if (exempt || (w->lo <= 0 && w->hi <= 0))
		return;
	v = rm_viol(w, a);
	s->worst = fmax(s->worst, v);
	s->bad += v > 0;
}

/* the angle pass's rule: the worst violation in the star falls (or every corner stays inside), no more corners are outside,
   and none that was inside leaves (left) */
static int rm_better(rm_stat const* b, rm_stat const* a, unsigned int left)
{
	return !left && a->bad <= b->bad && (a->worst < b->worst - 1e-12 || (a->worst <= 0 && b->worst <= 0));
}

/* the corners of the faces round v (and round u, for the collapse of edge (v, u): its two faces dropped) before (b) and after
   (a) v (and u) move to y; returns how many corners inside the bounds would leave them */
static unsigned int rm_star(rmsh const* w, unsigned int v, unsigned int u, double const* y, rm_stat* b, rm_stat* a)
{
	qaws_hm const* m = w->m;
	unsigned int i, c, left = 0, bounded = w->lo > 0 || w->hi > 0;
	rm_stat_init(b);
	rm_stat_init(a);
	for (i = 0; i < 2; i++)
	{
		unsigned int s = i ? u : v, o = i ? v : u;
		if (s == HM_NONE)
			continue;
		for (c = m->vc[s]; c != HM_NONE; c = hm_round(m, s, c))
		{
			unsigned int cc[3], k;
			double const *p0[3], *p1[3];
			int gone;
			cc[0] = c;
			cc[1] = hm_next(m, c);
			cc[2] = hm_prev(m, c);
			gone = o != HM_NONE && (m->cv[cc[1]] == o || m->cv[cc[2]] == o);
			if (gone && i)
				continue;   /* (counted from v) */
			for (k = 0; k < 3; k++)
			{
				p0[k] = rm_x(m, m->cv[cc[k]]);
				p1[k] = k ? p0[k] : y;
			}
			for (k = 0; k < 3; k++)
			{
				int ex = bounded && rm_exempt(w, cc[k]);
				double a0 = rm_angle(rm_M(w), p0[k], p0[(k + 1) % 3], p0[(k + 2) % 3]), a1;
				rm_stat_add(w, b, a0, ex);
				if (gone)
					continue;
				a1 = rm_angle(rm_M(w), p1[k], p1[(k + 1) % 3], p1[(k + 2) % 3]);
				rm_stat_add(w, a, a1, ex);
				left += bounded && !ex && rm_viol(w, a0) <= 0 && rm_viol(w, a1) > 0;
			}
		}
	}
	return left;
}

/* the corners of the faces round the vertices vs (each face once) against the bounds */
static void rm_patch_stat(rmsh const* w, unsigned int const* vs, unsigned int n, rm_stat* s)
{
	qaws_hm const* m = w->m;
	unsigned int fs[256], nf = 0, i, j, c, k;
	rm_stat_init(s);
	for (i = 0; i < n; i++)
		for (c = m->vc[vs[i]]; c != HM_NONE; c = hm_round(m, vs[i], c))
		{
			unsigned int f = m->cf[c];
			for (j = 0; j < nf && fs[j] != f; j++)
			{
			}
			if (j < nf || nf == 256)
				continue;
			fs[nf++] = f;
			for (k = m->fc[f]; k < m->fc[f] + 3; k++)
				rm_stat_add(w, s, rm_corner(w, k), rm_exempt(w, k));
		}
}

/* ---- the distance guard -------------------------------------------------- */

/* is a point of patch t within squared distance t2 of p: Gauss-Newton from (a, b), done once a point is near enough (any
   point of the reference bounds the distance); 1 yes, 0 no (a minimum inside the patch farther), -1 unknown (the patch's
   border reached: the closest point may lie beyond it). shortcut: a minimum inside the patch is taken as the distance, so a
   nearer point on another patch is missed and the operation refused (never wrongly kept); upgrade to the region search
   when references fold back within a tolerance */
static int rm_near(qaws_mesh_reference const* r, unsigned int t, double const* p, double a, double b, double t2)
{
	unsigned int it, k;
	for (it = 0; it < 6; it++)
	{
		double x[3], xa[3], xb[3], rr[3], g0 = 0, g1 = 0, h00 = 0, h01 = 0, h11 = 0, det, da, db, na, nb;
		mi_eval(&r->patch[t], a, b, x, xa, xb, NULL, NULL, NULL, NULL, NULL, NULL);
		for (k = 0; k < 3; k++)
		{
			rr[k] = x[k] - p[k];
			g0 += xa[k] * rr[k]; g1 += xb[k] * rr[k];
			h00 += xa[k] * xa[k]; h01 += xa[k] * xb[k]; h11 += xb[k] * xb[k];
		}
		if (rr[0] * rr[0] + rr[1] * rr[1] + rr[2] * rr[2] <= t2)
			return 1;
		det = h00 * h11 - h01 * h01;
		if (det <= 1e-300)
			return -1;
		da = -(h11 * g0 - h01 * g1) / det;
		db = -(h00 * g1 - h01 * g0) / det;
		na = a + da;
		nb = b + db;
		if (na < 0 || nb < 0 || na + nb > 1)
			return -1;
		if (fabs(da) + fabs(db) < 1e-12)
			return 0;
		a = na;
		b = nb;
	}
	return -1;
}

/* is every sample of the face (A, B, C) of region reg, its corners' footpoints LA, LB, LC, within the tolerance of the
   reference: its centre, its edges' midpoints (and their quarter points once strict). shortcut: samples only, a face bulging
   between them passes (dense sampling of the blob at 2e-3 finds 2.38e-3); upgrade to 130's certified deviation */
static int rm_within(rmsh* w, double const* A, double const* B, double const* C, qaws_mesh_location const* LA, qaws_mesh_location const* LB,
	qaws_mesh_location const* LC, unsigned int reg)
{
	static double const at[3] = { 0.5, 0.25, 0.75 };
	double const* P[3];
	qaws_mesh_location const* L[3];
	double s[3], t2 = w->tol * w->tol, e2 = 0;
	qaws_mesh_location l;
	unsigned int i, j, k;
	int one;
	if (w->tol <= 0)
		return 1;
	P[0] = A; P[1] = B; P[2] = C;
	L[0] = LA; L[1] = LB; L[2] = LC;
	for (i = 0; i < 3; i++)
		e2 = fmax(e2, rm_d2(P[i], P[(i + 1) % 3]));
	/* every point of a face is within its longest edge over sqrt 3 of a corner, and the corners are on the reference */
	if (e2 <= 3 * t2)
		return 1;
	one = LB->patch == LA->patch && LC->patch == LA->patch;
	for (i = 0; i < 4; i++)
		for (j = 0; j < (i < 3 && w->strict ? 3u : 1u); j++)
		{
			double t = at[j];
			int near = 0, h;
			for (k = 0; k < 3; k++)
				s[k] = i < 3 ? (1 - t) * P[i][k] + t * P[(i + 1) % 3][k] : (A[k] + B[k] + C[k]) / 3;
			/* from the sample's parameters interpolated when the corners share a patch, else from each corner's */
			if (one && LA->patch < w->r->nt)
			{
				double a = i < 3 ? (1 - t) * (double)L[i]->a + t * (double)L[(i + 1) % 3]->a : ((double)LA->a + (double)LB->a + (double)LC->a) / 3;
				double b = i < 3 ? (1 - t) * (double)L[i]->b + t * (double)L[(i + 1) % 3]->b : ((double)LA->b + (double)LB->b + (double)LC->b) / 3;
				near = rm_near(w->r, LA->patch, s, a, b, t2);
			}
			else
				for (k = 0; k < 3 && near <= 0; k++)
					if (L[k]->patch < w->r->nt && (k == 0 || L[k]->patch != L[0]->patch) && (k < 2 || L[k]->patch != L[1]->patch))
					{
						h = rm_near(w->r, L[k]->patch, s, (double)L[k]->a, (double)L[k]->b, t2);
						near = h > 0 ? 1 : h < 0 ? -1 : near;
					}
			if (near > 0)
				continue;
			if (!near || mr_project_region(w->r, s, LA->patch, reg, &w->s, &l) > t2)
				return 0;
		}
	return 1;
}

/* the edge to split when face f strays past the tolerance: its longest (Rivara's bisection: no angle falls under half the
   smallest there was); HM_NONE when it does not stray */
static unsigned int rm_face_stray(rmsh* w, unsigned int f)
{
	qaws_hm const* m = w->m;
	unsigned int c = m->fc[f], a = m->cv[c], b = m->cv[c + 1], d = m->cv[c + 2];
	return rm_within(w, rm_x(m, a), rm_x(m, b), rm_x(m, d), &m->foot[a], &m->foot[b], &m->foot[d], rm_region(w, f)) ? HM_NONE : rm_longest(m, f);
}

/* does the guard let an operation leave a straying face: in the loop, when one of the faces it replaces (fa, fb, or those
   round v and u) strays already (the splits refine those; held, a mesh still far from the tolerance could not be cleaned up) */
static int rm_lenient(rmsh* w, unsigned int fa, unsigned int fb, unsigned int v, unsigned int u)
{
	qaws_hm const* m = w->m;
	unsigned int i, c;
	if (w->strict)
		return 0;
	if ((fa != HM_NONE && rm_face_stray(w, fa) != HM_NONE) || (fb != HM_NONE && rm_face_stray(w, fb) != HM_NONE))
		return 1;
	for (i = 0; i < 2; i++)
	{
		unsigned int s = i ? u : v;
		if (s != HM_NONE)
			for (c = m->vc[s]; c != HM_NONE; c = hm_round(m, s, c))
				if (rm_face_stray(w, m->cf[c]) != HM_NONE)
					return 1;
	}
	return 0;
}

/* would the faces round v (and u, as in rm_star) stay within the tolerance once moved to y, its footpoint ly (or may they
   stray, as rm_lenient says) */
static int rm_star_within(rmsh* w, unsigned int v, unsigned int u, double const* y, qaws_mesh_location const* ly)
{
	qaws_hm const* m = w->m;
	unsigned int i, c;
	if (w->tol <= 0)
		return 1;
	for (i = 0; i < 2; i++)
	{
		unsigned int s = i ? u : v, o = i ? v : u;
		if (s == HM_NONE)
			continue;
		for (c = m->vc[s]; c != HM_NONE; c = hm_round(m, s, c))
		{
			unsigned int A = hm_to(m, c), B = m->cv[hm_prev(m, c)];
			if (A == o || B == o)
				continue;
			if (!rm_within(w, y, rm_x(m, A), rm_x(m, B), ly, &m->foot[A], &m->foot[B], rm_region(w, m->cf[c])))
				return rm_lenient(w, HM_NONE, HM_NONE, v, u);
		}
	}
	return 1;
}

/* ---- splits ------------------------------------------------------------- */

/* c's edge's midpoint on the reference: on the creases for a feature edge, in the edge's region otherwise (0: none found) */
static int rm_mid(rmsh* w, unsigned int c, double* x, qaws_mesh_location* l)
{
	qaws_hm const* m = w->m;
	unsigned int a = m->cv[c], b = hm_to(m, c), k;
	double mid[3], d;
	for (k = 0; k < 3; k++)
		mid[k] = 0.5 * (rm_x(m, a)[k] + rm_x(m, b)[k]);
	if (m->cflag[c] & HM_FEATURE)
		d = mr_project_crease(w->r, mid, m->foot[a].patch, &w->s, l);
	else
	{
		unsigned int reg = rm_region(w, m->cf[c]);
		d = mr_project_region(w->r, mid, rm_hint(w, a, b, reg), reg, &w->s, l);
	}
	if (d == HUGE_VAL || l->patch == HM_NONE)
		return 0;
	mr_eval(w->r, l, x, NULL);
	return 1;
}

/* is the midpoint of edge (a, b) surely within the tolerance (its ends on one patch, a point there near enough) */
static int rm_near_mid(rmsh const* w, unsigned int a, unsigned int b)
{
	qaws_mesh_location const *la = &w->m->foot[a], *lb = &w->m->foot[b];
	double mid[3];
	unsigned int k;
	if (la->patch != lb->patch || la->patch >= w->r->nt)
		return 0;
	for (k = 0; k < 3; k++)
		mid[k] = 0.5 * (rm_x(w->m, a)[k] + rm_x(w->m, b)[k]);
	return rm_near(w->r, la->patch, mid, 0.5 * ((double)la->a + (double)lb->a), 0.5 * ((double)la->b + (double)lb->b), w->tol * w->tol) > 0;
}

/* c's edge split at x, its footpoint l */
static unsigned int rm_split_at(rmsh* w, unsigned int c, double const* x, qaws_mesh_location const* l)
{
	unsigned int v = hm_split_edge(w->m, c, 0.5, x, 1);
	if (v == HM_NONE)
	{
		w->failed = 1;
		return HM_NONE;
	}
	w->m->foot[v] = *l;
	rm_set_size(w, v);
	w->rep.splits++;
	return v;
}

static void rm_split(rmsh* w)
{
	qaws_hm* m = w->m;
	unsigned int i;
	for (i = 0; (i > 0 || rm_snapshot(w)) && i < w->nep && !w->failed; i++)
	{
		unsigned int a, b, k, c = rm_edge(w, i);
		double mid[3], x[3];
		qaws_mesh_location l;
		int longer;
		if (c == HM_NONE)
			continue;
		a = m->cv[c];
		b = hm_to(m, c);
		longer = rm_d2(rm_x(m, a), rm_x(m, b)) > (16.0 / 9) * rm_target(w, a, b) * rm_target(w, a, b);
		if (!longer && (w->tol <= 0 || rm_near_mid(w, a, b)))
			continue;
		if (!rm_mid(w, c, x, &l))
			continue;
		/* with a tolerance, also an edge whose midpoint strays farther from the reference, whatever its length */
		for (k = 0; k < 3; k++)
			mid[k] = 0.5 * (rm_x(m, a)[k] + rm_x(m, b)[k]);
		if (!longer && rm_d2(x, mid) <= w->tol * w->tol)
			continue;
		rm_split_at(w, c, x, &l);
	}
}

/* ---- collapses ---------------------------------------------------------- */

static int rm_along_crease(rmsh* w, unsigned int p, unsigned int q);

/* where a collapse of c's edge puts the kept end: the edge's midpoint on the reference when both ends are off the features,
   on the crease when the edge runs along a feature line between two of its inner vertices; 0: the kept end stays */
static int rm_place(rmsh* w, unsigned int c, double* y, qaws_mesh_location* l)
{
	qaws_hm const* m = w->m;
	unsigned int a = m->cv[c], b = hm_to(m, c), fa = hm_feature_degree(m, a), fb = hm_feature_degree(m, b);
	if (((fa || fb) && (!(m->cflag[c] & HM_FEATURE) || fa != 2 || fb != 2)) || ((m->vflag[a] | m->vflag[b]) & HM_CORNER))
		return 0;   /* (a corner stays, a pinned one off the lines too) */
	return rm_mid(w, c, y, l);
}

/* may x, an end of an edge, go into the other end y, which moves to dst (its footpoint ldst)? (the link condition and the
   features checked by the caller) */
static int rm_can_remove(rmsh* w, unsigned int x, unsigned int y, double const* dst, qaws_mesh_location const* ldst)
{
	qaws_hm const* m = w->m;
	unsigned int ring[256], n, i, j, k, left;
	double e[3];
	int moved = rm_d2(dst, rm_x(m, y)) > 0;
	rm_stat b, a;
	for (j = 0; j < 1u + (unsigned int)moved; j++)
	{
		n = hm_ring(m, j ? y : x, ring, 256);
		if (n > 256)
			return 0;
		for (i = 0; i < n; i++)
		{
			double hi;
			if (ring[i] == y || ring[i] == x)
				continue;
			for (k = 0; k < 3; k++)
				e[k] = rm_x(m, ring[i])[k] - dst[k];
			hi = (4.0 / 3) * rm_target_along(w, y, ring[i], e);
			if (rm_d2(dst, rm_x(m, ring[i])) > hi * hi)
				return 0;
			/* no new edge along a crease over its own vertices (but x's feature edge, which goes on as y's) */
			if (!j && hm_find_edge(m, y, ring[i]) == HM_NONE && !(m->cflag[hm_find_edge(m, x, ring[i])] & HM_FEATURE) && rm_along_crease(w, y, ring[i]))
				return 0;
		}
	}
	if (rm_folds(w, x, dst, y) || (moved && rm_folds(w, y, dst, x)))
		return 0;
	/* no corner past 120 degrees, or past the largest there was (Zheng and Lv 2025) */
	left = rm_star(w, x, y, dst, &b, &a);
	if (a.amax > fmax(120 * RM_DEG, b.amax) || (w->bounded && !rm_better(&b, &a, left)))
		return 0;
	return rm_star_within(w, x, y, dst, ldst);
}

/* would x's and y's rings stay within 4/3 of their targets from dst, give or take slack? (rm_can_remove's first test, run on
   the plain midpoint before projecting it: most short edges fail it, and the projection is most of a collapse's cost) */
static int rm_reach(rmsh* w, unsigned int x, unsigned int y, double const* dst, double slack)
{
	qaws_hm const* m = w->m;
	unsigned int ring[256], n, i, j, k;
	double e[3];
	for (j = 0; j < 2; j++)
	{
		n = hm_ring(m, j ? y : x, ring, 256);
		if (n > 256)
			return 0;
		for (i = 0; i < n; i++)
		{
			double hi;
			if (ring[i] == y || ring[i] == x)
				continue;
			for (k = 0; k < 3; k++)
				e[k] = rm_x(m, ring[i])[k] - dst[k];
			hi = (4.0 / 3) * rm_target_along(w, y, ring[i], e) + slack;
			if (rm_d2(dst, rm_x(m, ring[i])) > hi * hi)
				return 0;
		}
	}
	return 1;
}

/* x merged into y: onto the placed point when the guards allow it there, else at y */
static int rm_try_collapse(rmsh* w, unsigned int c, unsigned int x, unsigned int y)
{
	qaws_hm* m = w->m;
	qaws_mesh_location l;
	double p[3];
	unsigned int k;
	int placed;
	if ((m->vflag[x] & HM_CORNER) || (hm_feature_degree(m, x) && !(m->cflag[c] & HM_FEATURE)) || !hm_can_collapse(m, c, x))
		return 0;
	rm_use(w, x, y);
	for (k = 0; k < 3; k++)
		p[k] = 0.5 * (rm_x(m, x)[k] + rm_x(m, y)[k]);
	l = m->foot[y];
	/* (the projection moves the midpoint by about the edge's sagitta: a tenth of the edge covers it) */
	placed = rm_reach(w, x, y, p, 0.1 * sqrt(rm_d2(rm_x(m, x), rm_x(m, y)))) && rm_place(w, c, p, &l);
	if (placed && !rm_can_remove(w, x, y, p, &l))
		placed = 0;
	if (!placed)
	{
		for (k = 0; k < 3; k++)
			p[k] = rm_x(m, y)[k];
		l = m->foot[y];
		if (!rm_can_remove(w, x, y, p, &l))
			return 0;
	}
	if (!hm_collapse(m, c, x, placed ? 0.5 : 1, p))
	{
		w->failed = 1;
		return 0;
	}
	m->foot[y] = l;
	if (placed)
		rm_set_size(w, y);
	w->rep.collapses++;
	rm_touch(w, y);
	return 1;
}

/* which end of c's edge goes: not the one on a feature, else the one of lower valence */
static void rm_pick(qaws_hm const* m, unsigned int c, unsigned int* x, unsigned int* y)
{
	unsigned int a = m->cv[c], b = hm_to(m, c), fa = hm_feature_degree(m, a), fb = hm_feature_degree(m, b);
	if (fa > fb || (!fa && !fb && hm_valence(m, a) > hm_valence(m, b)))
	{
		*x = b;
		*y = a;
	}
	else
	{
		*x = a;
		*y = b;
	}
}

static void rm_collapse(rmsh* w)
{
	qaws_hm* m = w->m;
	unsigned int i;
	for (i = 0; (i > 0 || rm_snapshot(w)) && i < w->nep && !w->failed; i++)
	{
		unsigned int a, b, x, y, c = rm_edge(w, i);
		double lo;
		if (c == HM_NONE)
			continue;
		a = m->cv[c];
		b = hm_to(m, c);
		lo = 0.8 * rm_target(w, a, b);
		if (rm_d2(rm_x(m, a), rm_x(m, b)) >= lo * lo)
			continue;
		rm_pick(m, c, &x, &y);
		if (!rm_try_collapse(w, c, x, y) && !w->failed)
			rm_try_collapse(w, c, y, x);
	}
}

/* needles (a corner under 15 degrees): the face's shortest edge collapsed as rm_collapse would */
static void rm_needles(rmsh* w)
{
	qaws_hm* m = w->m;
	unsigned int i, s;
	for (i = 0; (i > 0 || rm_snapshot(w)) && i < w->nep && !w->failed; i++)
	{
		unsigned int x, y, c = rm_edge(w, i);
		if (c == HM_NONE)
			continue;
		for (s = 0; s < 2; s++)
		{
			unsigned int e = s ? m->co[c] : c;
			double l;
			if (e == HM_NONE)
				continue;
			rm_use(w, m->cv[c], hm_to(m, c));
			l = rm_len(w, e);
			if (l > rm_len(w, hm_next(m, e)) || l > rm_len(w, hm_prev(m, e)) || rm_corner(w, hm_prev(m, e)) >= 15 * RM_DEG)
				continue;
			rm_pick(m, c, &x, &y);
			if (!rm_try_collapse(w, c, x, y) && !w->failed)
				rm_try_collapse(w, c, y, x);
			break;
		}
	}
}

/* ---- flips -------------------------------------------------------------- */

/* how far the valence at corner c is from ideal: inside, the vertex's valence over 6; on features and the border, the faces
   of the sector holding c over the whole number of 60-degree corners nearest its angle (at least 1): so 6 on a smooth crease
   and 4 on a straight border, per side */
static int rm_off(qaws_hm const* m, unsigned int c)
{
	unsigned int faces, ideal;
	double th;
	if (!hm_feature_degree(m, m->cv[c]))
		return (int)hm_valence(m, m->cv[c]) - 6;
	th = rm_sector(m, c, &faces);
	ideal = (unsigned int)floor(th / (60 * RM_DEG) + 0.5 + 1e-6);
	return (int)faces - (int)(ideal ? ideal : 1);
}

/* would a new edge from p to q run along a crease (both ends on features, its
   middle on one): a chord over a feature line's own vertices? */
static int rm_along_crease(rmsh* w, unsigned int p, unsigned int q)
{
	qaws_hm const* m = w->m;
	double mid[3], d2 = 0;
	qaws_mesh_location l;
	unsigned int k;
	if (!hm_feature_degree(m, p) || !hm_feature_degree(m, q))
		return 0;
	for (k = 0; k < 3; k++)
	{
		mid[k] = 0.5 * (rm_x(m, p)[k] + rm_x(m, q)[k]);
		d2 += (rm_x(m, p)[k] - rm_x(m, q)[k]) * (rm_x(m, p)[k] - rm_x(m, q)[k]);
	}
	return mr_project_crease(w->r, mid, m->foot[p].patch, &w->s, &l) < 1e-6 * d2;
}

/* may c's edge flip: inside and off the features, the new edge off the creases, the new pair (a, q, p) and (q, b, p) facing
   the old pair's way, the poorer of them no poorer than before or than fair, both within the tolerance */
static int rm_flip_ok(rmsh* w, unsigned int c)
{
	qaws_hm const* m = w->m;
	unsigned int d = m->co[c], a, b, p, q, k;
	double const *A, *B, *P, *Q;
	double n1[3], n2[3], n0[3], g[3], u[3], t[3];
	if (d == HM_NONE || (m->cflag[c] & HM_FEATURE) || !hm_can_flip(m, c))
		return 0;
	a = m->cv[c]; b = hm_to(m, c); p = m->cv[hm_prev(m, c)]; q = m->cv[hm_prev(m, d)];
	if (rm_along_crease(w, p, q))
		return 0;
	A = rm_x(m, a); B = rm_x(m, b); P = rm_x(m, p); Q = rm_x(m, q);
	for (k = 0; k < 3; k++) { u[k] = Q[k] - A[k]; t[k] = P[k] - A[k]; }
	n1[0] = u[1] * t[2] - u[2] * t[1]; n1[1] = u[2] * t[0] - u[0] * t[2]; n1[2] = u[0] * t[1] - u[1] * t[0];
	for (k = 0; k < 3; k++) { u[k] = B[k] - Q[k]; t[k] = P[k] - Q[k]; }
	n2[0] = u[1] * t[2] - u[2] * t[1]; n2[1] = u[2] * t[0] - u[0] * t[2]; n2[2] = u[0] * t[1] - u[1] * t[0];
	hm_face_normal(m, m->cf[c], n0);
	hm_face_normal(m, m->cf[d], g);
	for (k = 0; k < 3; k++)
		n0[k] += g[k];
	if (n1[0] * n0[0] + n1[1] * n0[1] + n1[2] * n0[2] <= 0 || n2[0] * n0[0] + n2[1] * n0[1] + n2[2] * n0[2] <= 0
		|| n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2] <= 0)
		return 0;
	{
		double qo = fmin(rm_quality(rm_M(w), A, B, P), rm_quality(rm_M(w), B, A, Q)), qn = fmin(rm_quality(rm_M(w), A, Q, P), rm_quality(rm_M(w), Q, B, P));
		if (qn < qo && qn < RM_FAIR)
			return 0;
	}
	if (w->tol > 0)
	{
		unsigned int reg = rm_region(w, m->cf[c]);
		if ((!rm_within(w, A, Q, P, &m->foot[a], &m->foot[q], &m->foot[p], reg) || !rm_within(w, Q, B, P, &m->foot[q], &m->foot[b], &m->foot[p], reg))
			&& !rm_lenient(w, m->cf[c], m->cf[d], HM_NONE, HM_NONE))
			return 0;
	}
	return 1;
}

/* the bounds' verdict on flipping c's edge: the new corners take the exemption of the old corner at their vertex */
static int rm_flip_better(rmsh const* w, unsigned int c)
{
	qaws_hm const* m = w->m;
	unsigned int d = m->co[c];
	double const *A = rm_x(m, m->cv[c]), *B = rm_x(m, hm_to(m, c)), *P = rm_x(m, m->cv[hm_prev(m, c)]), *Q = rm_x(m, m->cv[hm_prev(m, d)]);
	int ea = rm_exempt(w, c), eb = rm_exempt(w, hm_next(m, c)), ep = rm_exempt(w, hm_prev(m, c)), eq = rm_exempt(w, hm_prev(m, d));
	rm_stat b, a;
	rm_stat_init(&b);
	rm_stat_init(&a);
	rm_stat_add(w, &b, rm_angle(rm_M(w), A, B, P), ea); rm_stat_add(w, &b, rm_angle(rm_M(w), B, P, A), eb); rm_stat_add(w, &b, rm_angle(rm_M(w), P, A, B), ep);
	rm_stat_add(w, &b, rm_angle(rm_M(w), B, A, Q), eb); rm_stat_add(w, &b, rm_angle(rm_M(w), A, Q, B), ea); rm_stat_add(w, &b, rm_angle(rm_M(w), Q, B, A), eq);
	rm_stat_add(w, &a, rm_angle(rm_M(w), A, Q, P), ea); rm_stat_add(w, &a, rm_angle(rm_M(w), Q, P, A), eq); rm_stat_add(w, &a, rm_angle(rm_M(w), P, A, Q), ep);
	rm_stat_add(w, &a, rm_angle(rm_M(w), Q, B, P), eq); rm_stat_add(w, &a, rm_angle(rm_M(w), B, P, Q), eb); rm_stat_add(w, &a, rm_angle(rm_M(w), P, Q, B), ep);
	return rm_better(&b, &a, 0);
}

/* flips toward the ideal valences, the new pair's smallest angle no smaller than the old pair's or than 30 degrees */
static void rm_flip(rmsh* w)
{
	qaws_hm* m = w->m;
	unsigned int i;
	for (i = 0; (i > 0 || rm_snapshot(w)) && i < w->nep && !w->failed; i++)
	{
		unsigned int d, a, b, p, q, c = rm_edge(w, i);
		int va, vb, vp, vq, before, after;
		double const *A, *B, *P, *Q;
		if (c == HM_NONE || (d = m->co[c]) == HM_NONE || (m->cflag[c] & HM_FEATURE))
			continue;
		a = m->cv[c]; b = hm_to(m, c); p = m->cv[hm_prev(m, c)]; q = m->cv[hm_prev(m, d)];
		if (p == q)
			continue;
		rm_use(w, a, b);
		va = rm_off(m, c);
		vb = rm_off(m, hm_next(m, c));
		vp = rm_off(m, hm_prev(m, c));
		vq = rm_off(m, hm_prev(m, d));
		before = abs(va) + abs(vb) + abs(vp) + abs(vq);
		after = abs(va - 1) + abs(vb - 1) + abs(vp + 1) + abs(vq + 1);
		if (after >= before)
			continue;
		A = rm_x(m, a); B = rm_x(m, b); P = rm_x(m, p); Q = rm_x(m, q);
		if (fmin(rm_tri_min(rm_M(w), A, Q, P), rm_tri_min(rm_M(w), Q, B, P)) < fmin(fmin(rm_tri_min(rm_M(w), A, B, P), rm_tri_min(rm_M(w), B, A, Q)), 30 * RM_DEG))
			continue;
		if (!rm_flip_ok(w, c) || (w->bounded && !rm_flip_better(w, c)))
			continue;
		if (!hm_flip(m, c))
		{
			w->failed = 1;
			break;
		}
		w->rep.flips++;
		rm_touch(w, a); rm_touch(w, b); rm_touch(w, p); rm_touch(w, q);
	}
}

/* caps (an angle over 120 degrees): the edge across it flipped when that makes the pair better (as pmp's remove_caps) */
static void rm_caps(rmsh* w)
{
	qaws_hm* m = w->m;
	unsigned int i;
	double lim = cos(120 * RM_DEG);
	for (i = 0; (i > 0 || rm_snapshot(w)) && i < w->nep && !w->failed; i++)
	{
		unsigned int d, c = rm_edge(w, i);
		double const *A, *B, *P, *Q;
		if (c == HM_NONE || (d = m->co[c]) == HM_NONE || (m->cflag[c] & HM_FEATURE) || !hm_can_flip(m, c)
			|| rm_along_crease(w, m->cv[hm_prev(m, c)], m->cv[hm_prev(m, d)]))
			continue;
		rm_use(w, m->cv[c], hm_to(m, c));
		A = rm_x(m, m->cv[c]); B = rm_x(m, hm_to(m, c)); P = rm_x(m, m->cv[hm_prev(m, c)]); Q = rm_x(m, m->cv[hm_prev(m, d)]);
		if (rm_cos(rm_M(w), P, A, B) > lim && rm_cos(rm_M(w), Q, B, A) > lim)
			continue;
		if (fmin(rm_quality(rm_M(w), A, Q, P), rm_quality(rm_M(w), Q, B, P)) <= fmin(rm_quality(rm_M(w), A, B, P), rm_quality(rm_M(w), B, A, Q)))
			continue;
		if (w->tol > 0)
		{
			unsigned int reg = rm_region(w, m->cf[c]);
			if ((!rm_within(w, A, Q, P, &m->foot[m->cv[c]], &m->foot[m->cv[hm_prev(m, d)]], &m->foot[m->cv[hm_prev(m, c)]], reg)
				|| !rm_within(w, Q, B, P, &m->foot[m->cv[hm_prev(m, d)]], &m->foot[hm_to(m, c)], &m->foot[m->cv[hm_prev(m, c)]], reg))
				&& !rm_lenient(w, m->cf[c], m->cf[d], HM_NONE, HM_NONE))
				continue;
		}
		if (!hm_flip(m, c))
		{
			w->failed = 1;
			break;
		}
		w->rep.flips++;
	}
}

/* Delaunay flips to finish (Lawson's queue): an inner edge off the features flipped while the angles across it sum past pi,
   under the flips' guards, each flip queueing the four edges round it */
static void rm_delaunay(rmsh* w)
{
	qaws_hm* m = w->m;
	unsigned int n, flips = 0, cap = m->nv > 0xFFFFFFFFu / 100 ? 0xFFFFFFFFu : 100 * m->nv;   /* a flip budget: saturate, never wrap */
	if (!rm_snapshot(w))
		return;
	for (n = w->nep; n > 0 && flips < cap && !w->failed;)
	{
		unsigned int a = w->ep[2 * n - 2], b = w->ep[2 * n - 1], c, d, p, q;
		n--;
		c = hm_find_edge(m, a, b);
		if (c == HM_NONE || (d = m->co[c]) == HM_NONE || (m->cflag[c] & HM_FEATURE))
			continue;
		a = m->cv[c]; b = hm_to(m, c); p = m->cv[hm_prev(m, c)]; q = m->cv[hm_prev(m, d)];
		rm_use(w, a, b);
		if (rm_angle(rm_M(w), rm_x(m, p), rm_x(m, a), rm_x(m, b)) + rm_angle(rm_M(w), rm_x(m, q), rm_x(m, b), rm_x(m, a)) <= RM_PI + 1e-12 || !rm_flip_ok(w, c))
			continue;
		if (!hm_flip(m, c))
		{
			w->failed = 1;
			break;
		}
		flips++;
		w->nep = n;
		if (!rm_push(w, a, p) || !rm_push(w, p, b) || !rm_push(w, b, q) || !rm_push(w, q, a))
			break;
		n = w->nep;
	}
	w->rep.flips += flips;
}

/* ---- tangential relaxation ---------------------------------------------- */

/* v relaxed: inside, to the size-weighted area centroid of its ring (each face's centroid weighted by its area over its mean
   size squared, so a graded mesh keeps its grading); on a feature line, to the size-weighted centre of its two line edges
   (each weighted by its length over its size); a full step on the reference, half when that folds or a guard refuses */
static int rm_smooth_vertex(rmsh* w, unsigned int v)
{
	qaws_hm* m = w->m;
	unsigned int fd, k, half;
	double d[3], x[3], nrm[3], dd;
	qaws_mesh_location l;
	if ((m->vflag[v] & (HM_DELETED | HM_CORNER)) || m->vc[v] == HM_NONE)
		return 0;
	fd = hm_feature_degree(m, v);
	if (fd != 0 && fd != 2)
		return 0;
	rm_use(w, v, v);
	if (fd == 0)
	{
		unsigned int c;
		double cen[3] = { 0, 0, 0 }, W = 0, dn;
		for (c = m->vc[v]; c != HM_NONE; c = hm_round(m, v, c))
		{
			unsigned int a = hm_to(m, c), b = m->cv[hm_prev(m, c)];
			double const *P = rm_x(m, v), *A = rm_x(m, a), *B = rm_x(m, b);
			double u[3], t[3], n[3], h = (rm_size_of(w, v) + rm_size_of(w, a) + rm_size_of(w, b)) / 3, wt;
			for (k = 0; k < 3; k++) { u[k] = A[k] - P[k]; t[k] = B[k] - P[k]; }
			n[0] = u[1] * t[2] - u[2] * t[1]; n[1] = u[2] * t[0] - u[0] * t[2]; n[2] = u[0] * t[1] - u[1] * t[0];
			wt = 0.5 * sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]) / (h * h);
			for (k = 0; k < 3; k++)
				cen[k] += wt * (P[k] + A[k] + B[k]) / 3;
			W += wt;
		}
		if (W <= 0)
			return 0;
		mr_eval(w->r, &m->foot[v], x, nrm);
		for (k = 0; k < 3; k++)
			d[k] = cen[k] / W - rm_x(m, v)[k];
		dn = d[0] * nrm[0] + d[1] * nrm[1] + d[2] * nrm[2];
		for (k = 0; k < 3; k++)
			d[k] -= dn * nrm[k];
	}
	else
	{
		/* along the feature line: its two line neighbours f, the move along the chord between them */
		unsigned int f[2], nf = 0, c, last = HM_NONE;
		double t[3], dt, wt[2], cen[3] = { 0, 0, 0 };
		for (c = m->vc[v]; c != HM_NONE; c = hm_round(m, v, c))
		{
			if ((m->cflag[c] & HM_FEATURE) && nf < 2)
				f[nf++] = hm_to(m, c);
			last = c;
		}
		if (nf < 2 && last != HM_NONE && hm_is_border_vertex(m, v) && (m->cflag[hm_prev(m, last)] & HM_FEATURE))
			f[nf++] = m->cv[hm_prev(m, last)];
		if (nf != 2)
			return 0;
		for (c = 0; c < 2; c++)
			wt[c] = sqrt(rm_d2(rm_x(m, v), rm_x(m, f[c]))) / rm_target(w, v, f[c]);
		if (wt[0] + wt[1] <= 0)
			return 0;
		for (k = 0; k < 3; k++)
		{
			t[k] = rm_x(m, f[1])[k] - rm_x(m, f[0])[k];
			cen[k] = (wt[0] * (rm_x(m, v)[k] + rm_x(m, f[0])[k]) + wt[1] * (rm_x(m, v)[k] + rm_x(m, f[1])[k])) / (2 * (wt[0] + wt[1]));
		}
		dt = t[0] * t[0] + t[1] * t[1] + t[2] * t[2];
		if (dt <= 0)
			return 0;
		dt = ((cen[0] - rm_x(m, v)[0]) * t[0] + (cen[1] - rm_x(m, v)[1]) * t[1] + (cen[2] - rm_x(m, v)[2]) * t[2]) / dt;
		for (k = 0; k < 3; k++)
			d[k] = dt * t[k];
	}
	dd = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
	if (dd <= 1e-24 * w->L * w->L)
		return 0;
	for (half = 0; half < 2; half++)
	{
		if (half)
			for (k = 0; k < 3; k++)
				d[k] *= 0.5;
		if (fd == 0)
			mr_walk(w->r, m->foot[v], d, 1, &l);   /* (stopped short of a crease: kept at its end) */
		else
		{
			double y[3];
			for (k = 0; k < 3; k++)
				y[k] = rm_x(m, v)[k] + d[k];
			if (mr_project_crease(w->r, y, m->foot[v].patch, &w->s, &l) == HUGE_VAL || l.patch == HM_NONE)
				return 0;
		}
		mr_eval(w->r, &l, x, NULL);
		if (rm_folds(w, v, x, HM_NONE) || !rm_star_within(w, v, HM_NONE, x, &l))
			continue;
		if (w->bounded)
		{
			rm_stat b, a;
			unsigned int left = rm_star(w, v, HM_NONE, x, &b, &a);
			if (!rm_better(&b, &a, left))
				continue;
		}
		w->moved = fmax(w->moved, sqrt(rm_d2(x, rm_x(m, v))));
		for (k = 0; k < 3; k++)
			rm_x(m, v)[k] = x[k];
		m->foot[v] = l;
		rm_set_size(w, v);
		w->rep.moves++;
		return 1;
	}
	return 0;
}

static void rm_smooth(rmsh* w)
{
	unsigned int v;
	for (v = 0; v < w->m->nv; v++)
		rm_smooth_vertex(w, v);
}

/* ---- the distance guard made good ---------------------------------------- */

/* every face within the tolerance: the longest edge of a face whose samples stray split, round after round, each round
   finished by Delaunay flips; then needles collapsed and a relaxation sweep (all under the guard) */
static void rm_tolerance(rmsh* w)
{
	qaws_hm* m = w->m;
	unsigned int round, f, n, e;
	if (w->tol <= 0)
		return;
	for (round = 0; round < 30 && !w->failed; round++)
	{
		unsigned int nf = m->nf;
		for (f = 0, n = 0; f < nf && !w->failed; f++)
		{
			double x[3];
			qaws_mesh_location l;
			if (m->fflag[f] & HM_DELETED)
				continue;
			if ((e = rm_face_stray(w, f)) == HM_NONE)
				continue;
			if (rm_mid(w, e, x, &l) && rm_split_at(w, e, x, &l) != HM_NONE)
				n++;
		}
		hm_compact(m);
		if (!n)
			break;
		rm_delaunay(w);
	}
	rm_needles(w);
	hm_compact(m);
	rm_smooth(w);
}

/* ---- the angle-bounded pass (Wang, Yan et al. 2019) ---------------------- */

typedef struct rm_bad
{
	double e;                       /* how far its worst corner is out */
	unsigned int v[3];              /* the face, by its vertices */
} rm_bad;

static int rm_bad_cmp(void const* a, void const* b)
{
	double d = ((rm_bad const*)b)->e - ((rm_bad const*)a)->e;
	return d > 0 ? 1 : d < 0 ? -1 : 0;
}

/* face f's worst corner not exempt above the upper bound (large) or under the lower, and by how much; HM_NONE: none */
static unsigned int rm_worst(rmsh const* w, unsigned int f, int large, double* out)
{
	qaws_hm const* m = w->m;
	unsigned int c, best = HM_NONE;
	double e = 0;
	if ((large && w->hi <= 0) || (!large && w->lo <= 0))
		return HM_NONE;
	for (c = m->fc[f]; c < m->fc[f] + 3; c++)
	{
		double a = rm_corner(w, c), x = large ? a - w->hi : w->lo - a;
		if (x > e && !rm_exempt(w, c))
		{
			e = x;
			best = c;
		}
	}
	*out = e;
	return best;
}

/* the faces with a corner out of the bounds (large: above), worst first */
static unsigned int rm_gather(rmsh* w, int large, rm_bad** list, unsigned int* cap)
{
	qaws_hm const* m = w->m;
	unsigned int f, n = 0, k;
	double e;
	for (f = 0; f < m->nf; f++)
	{
		if ((m->fflag[f] & HM_DELETED) || rm_worst(w, f, large, &e) == HM_NONE)
			continue;
		if (n == *cap)
		{
			unsigned int nc = qaws_internal_grow(*cap, (size_t)n + 1, 64);
			rm_bad* g = nc ? (rm_bad*)qaws_internal_realloc(*list, qaws_internal_mul2(sizeof(rm_bad), nc)) : NULL;
			if (!g)
			{
				w->failed = 1;
				break;
			}
			*list = g;
			*cap = nc;
		}
		(*list)[n].e = e;
		for (k = 0; k < 3; k++)
			(*list)[n].v[k] = m->cv[m->fc[f] + k];
		n++;
	}
	qsort(*list, n, sizeof(rm_bad), rm_bad_cmp);
	return n;
}

/* the face a listed one still is, or HM_NONE */
static unsigned int rm_face_of(qaws_hm const* m, rm_bad const* t)
{
	unsigned int c;
	if ((m->vflag[t->v[0]] | m->vflag[t->v[1]] | m->vflag[t->v[2]]) & HM_DELETED)
		return HM_NONE;
	c = hm_find(m, t->v[0], t->v[1]);
	return c != HM_NONE && m->cv[hm_prev(m, c)] == t->v[2] ? m->cf[c] : HM_NONE;
}

/* the worst violation over every corner not exempt, and how many are out */
static double rm_census(rmsh const* w, unsigned int* bad)
{
	qaws_hm const* m = w->m;
	unsigned int f, c;
	double worst = -RM_PI;
	*bad = 0;
	for (f = 0; f < m->nf; f++)
		if (!(m->fflag[f] & HM_DELETED))
			for (c = m->fc[f]; c < m->fc[f] + 3; c++)
			{
				double v = rm_viol(w, rm_corner(w, c));
				if (v > 0 && !rm_exempt(w, c))
				{
					worst = fmax(worst, v);
					(*bad)++;
				}
			}
	return worst;
}

static int rm_relocate(rmsh* w, unsigned int v);

/* v and its ring relaxed freely, then moved toward the bounds */
static void rm_settle(rmsh* w, unsigned int v)
{
	unsigned int ring[32], nr = hm_ring(w->m, v, ring, 32), j, k;
	int bd = w->bounded;
	nr = nr > 32 ? 32 : nr;
	w->bounded = 0;
	for (j = 0; j < 3; j++)
	{
		rm_smooth_vertex(w, v);
		for (k = 0; k < nr; k++)
			rm_smooth_vertex(w, ring[k]);
	}
	rm_relocate(w, v);
	for (k = 0; k < nr; k++)
		rm_relocate(w, ring[k]);
	rm_relocate(w, v);
	w->bounded = bd;
}

/* how a patch stands for the split's choice: its worst violation of the bounds first, then its corners out of them, then the
   root mean square of its angles' distance from 60 degrees */
static double rm_score(rm_stat const* s)
{
	return 1e4 * fmax(s->worst, 0) + s->bad + (s->n ? sqrt(s->sq / s->n) : 0);
}

/* the large-angle split: e's edge split at its midpoint on the reference; then (off the features) no flip or the flip of one
   of the sides round the new vertex away from the apex (the large corner's vertex), each tried with the new vertex and its
   ring settled, the best kept; all of it undone unless the angles round it got better, or (need) the apex's sector has too
   few faces to be inside the bounds at all: then the face it gains is kept for later rounds to polish */
static int rm_split_large(rmsh* w, unsigned int e, unsigned int apex, int need)
{
	qaws_hm* m = w->m;
	unsigned int vs[16], nvs = 0, j, k, v, a = m->cv[e], c, side[2], opt[8][2], no = 1, o, best = 0, far = HM_NONE;
	double x[3], sx[16][3], ssz[16], bscore = HUGE_VAL;
	qaws_mesh_location l, sl[16];
	qaws_mesh_report rep = w->rep;
	rm_stat sb, sa;
	rm_use(w, a, hm_to(m, e));
	/* the region judged: the faces round the quad's corners and round the corners across its sides */
	side[0] = e;
	side[1] = m->co[e];
	for (j = 0; j < 2; j++)
	{
		unsigned int s = side[j], t;
		if (s == HM_NONE)
			continue;
		vs[nvs++] = m->cv[s];
		vs[nvs++] = m->cv[hm_prev(m, s)];
		for (t = hm_next(m, s); t != s; t = hm_next(m, t))
			if (m->co[t] != HM_NONE)
				vs[nvs++] = m->cv[hm_prev(m, m->co[t])];
	}
	vs[nvs++] = hm_to(m, e);
	rm_patch_stat(w, vs, nvs, &sb);
	if (!rm_mid(w, e, x, &l) || (v = rm_split_at(w, e, x, &l)) == HM_NONE)
		return 0;
	vs[nvs++] = v;
	for (k = 0; k < nvs; k++)
	{
		memcpy(sx[k], rm_x(m, vs[k]), sizeof sx[k]);
		sl[k] = m->foot[vs[k]];
		ssz[k] = rm_size_of(w, vs[k]);
	}
	/* the options: none, or a side round v off the apex (none on a feature: its halves stay) */
	for (c = m->vc[v]; c != HM_NONE && no < 8; c = hm_round(m, v, c))
	{
		unsigned int t = hm_next(m, c);
		if (!hm_feature_degree(m, v) && m->cv[t] != apex && hm_to(m, t) != apex && rm_flip_ok(w, t))
		{
			opt[no][0] = m->cv[t];
			opt[no][1] = hm_to(m, t);
			no++;
		}
	}
	/* each tried (flipped, settled, judged, undone), then the best done again */
	rm_stat_init(&sa);
	sa.worst = HUGE_VAL;
	for (o = 0; o <= no; o++)
	{
		unsigned int oo = o < no ? o : best;
		double score;
		far = HM_NONE;
		if (oo)
		{
			c = hm_find(m, opt[oo][0], opt[oo][1]);
			if (c == HM_NONE || m->co[c] == HM_NONE)
				continue;
			far = m->cv[hm_prev(m, m->co[c])];
			if (!hm_flip(m, c))
			{
				w->failed = 1;
				return 1;
			}
		}
		rm_settle(w, v);
		rm_patch_stat(w, vs, nvs, &sa);
		if (o == no)
			break;
		score = rm_score(&sa);
		if (score < bscore)
		{
			bscore = score;
			best = oo;
		}
		for (k = 0; k < nvs; k++)
		{
			memcpy(rm_x(m, vs[k]), sx[k], sizeof sx[k]);
			m->foot[vs[k]] = sl[k];
			hm_attr_at(m, w->sa, vs[k])[0] = ssz[k];
		}
		if (far != HM_NONE && ((c = hm_find_edge(m, v, far)) == HM_NONE || !hm_flip(m, c)))
		{
			w->failed = 1;
			return 1;
		}
	}
	if ((need || rm_better(&sb, &sa, 0)) && rm_star_within(w, v, HM_NONE, rm_x(m, v), &m->foot[v]))
	{
		unsigned int ring[32], nr = hm_ring(m, v, ring, 32);
		rm_touch(w, v);
		for (k = 0; k < nr && k < 32; k++)
			rm_touch(w, ring[k]);
		return 1;
	}
	/* undone: the places back, the flip back, v merged into a */
	for (k = 0; k < nvs; k++)
	{
		memcpy(rm_x(m, vs[k]), sx[k], sizeof sx[k]);
		m->foot[vs[k]] = sl[k];
		hm_attr_at(m, w->sa, vs[k])[0] = ssz[k];
	}
	if (far != HM_NONE && ((c = hm_find_edge(m, v, far)) == HM_NONE || !hm_flip(m, c)))
	{
		w->failed = 1;
		return 1;
	}
	c = hm_find_edge(m, v, a);
	if (c != HM_NONE && hm_can_collapse(m, c, v))
	{
		memcpy(x, rm_x(m, a), sizeof x);
		if (!hm_collapse(m, c, v, 1, x))
		{
			w->failed = 1;
			return 1;
		}
		w->rep = rep;
	}
	return 0;
}

/* v moved by a pattern search on the reference: 16 directions in its tangent plane (2 along its line on a feature), the step
   halved from a quarter of its size, toward the smallest worst violation of the bounds in its star (for a small angle, its
   largest smallest angle). shortcut: a pattern search, not Wang et al.'s exact feasible region; upgrade if the pass stalls
   on stars whose feasible region is a sliver the 16 directions miss (8 left the cube at 0.15 one corner at 86.45) */
static int rm_relocate(rmsh* w, unsigned int v)
{
	qaws_hm* m = w->m;
	unsigned int fd, it, k, dirs, f[2] = { HM_NONE, HM_NONE }, nf = 0;
	double step, cur, cx[3], t1[3], t2[3] = { 0, 0, 0 }, nrm[3], q[3], bx[3] = { 0, 0, 0 };
	qaws_mesh_location cl, bl = { 0, 0, 0 };
	rm_stat b0, a0;
	int found = 0;
	if ((m->vflag[v] & (HM_DELETED | HM_CORNER)) || m->vc[v] == HM_NONE)
		return 0;
	rm_use(w, v, v);
	fd = hm_feature_degree(m, v);
	if (fd != 0 && fd != 2)
		return 0;
	if (fd == 2)
	{
		unsigned int c, last = HM_NONE;
		for (c = m->vc[v]; c != HM_NONE; c = hm_round(m, v, c))
		{
			if ((m->cflag[c] & HM_FEATURE) && nf < 2)
				f[nf++] = hm_to(m, c);
			last = c;
		}
		if (nf < 2 && last != HM_NONE && hm_is_border_vertex(m, v) && (m->cflag[hm_prev(m, last)] & HM_FEATURE))
			f[nf++] = m->cv[hm_prev(m, last)];
		if (nf != 2)
			return 0;
	}
	rm_star(w, v, HM_NONE, rm_x(m, v), &b0, &a0);
	cur = b0.worst;
	memcpy(cx, rm_x(m, v), sizeof cx);
	cl = m->foot[v];
	dirs = fd ? 2 : 16;
	for (it = 0, step = 0.25 * rm_size_of(w, v); it < 10; it++)
	{
		int improved = 0;
		if (fd)
		{
			double tt = 0;
			for (k = 0; k < 3; k++)
			{
				t1[k] = rm_x(m, f[1])[k] - rm_x(m, f[0])[k];
				tt += t1[k] * t1[k];
			}
			for (k = 0; k < 3; k++)
				t1[k] /= sqrt(tt);
		}
		else
		{
			double nn;
			mr_eval(w->r, &cl, q, nrm);
			k = fabs(nrm[0]) < 0.6 ? 0 : fabs(nrm[1]) < 0.6 ? 1 : 2;
			t1[0] = t1[1] = t1[2] = 0;
			t1[k] = 1;
			nn = t1[0] * nrm[0] + t1[1] * nrm[1] + t1[2] * nrm[2];
			for (k = 0; k < 3; k++)
				t1[k] -= nn * nrm[k];
			nn = sqrt(t1[0] * t1[0] + t1[1] * t1[1] + t1[2] * t1[2]);
			for (k = 0; k < 3; k++)
				t1[k] /= nn;
			t2[0] = nrm[1] * t1[2] - nrm[2] * t1[1]; t2[1] = nrm[2] * t1[0] - nrm[0] * t1[2]; t2[2] = nrm[0] * t1[1] - nrm[1] * t1[0];
		}
		for (k = 0; k < dirs; k++)
		{
			double d[3], y[3], ang = 2 * RM_PI * k / dirs;
			qaws_mesh_location l;
			rm_stat b, a;
			unsigned int j, left;
			for (j = 0; j < 3; j++)
				d[j] = step * (cos(ang) * t1[j] + (fd ? 0 : sin(ang) * t2[j]));
			if (fd)
			{
				for (j = 0; j < 3; j++)
					y[j] = cx[j] + d[j];
				if (mr_project_crease(w->r, y, cl.patch, &w->s, &l) == HUGE_VAL || l.patch == HM_NONE)
					continue;
			}
			else
				mr_walk(w->r, cl, d, 1, &l);
			mr_eval(w->r, &l, y, NULL);
			if (rm_folds(w, v, y, HM_NONE))
				continue;
			left = rm_star(w, v, HM_NONE, y, &b, &a);
			if (left || a.bad > b.bad || a.worst >= cur - 1e-12 || !rm_star_within(w, v, HM_NONE, y, &l))
				continue;
			cur = a.worst;
			memcpy(bx, y, sizeof bx);
			bl = l;
			improved = 1;
		}
		if (improved)
		{
			memcpy(cx, bx, sizeof cx);
			cl = bl;
			found = 1;
		}
		else
			step *= 0.5;
	}
	if (!found)
		return 0;
	memcpy(rm_x(m, v), cx, sizeof cx);
	m->foot[v] = cl;
	rm_set_size(w, v);
	w->rep.moves++;
	rm_touch(w, v);
	return 1;
}

/* is a listed face gone, or its corners now inside the bound (large: the upper) */
static int rm_fixed(rmsh const* w, rm_bad const* t, int large)
{
	unsigned int f = rm_face_of(w->m, t);
	double e;
	return f == HM_NONE || rm_worst(w, f, large, &e) == HM_NONE;
}

/* the face's corners moved in turn (the vertex at corner c first, the vertex across the face's longest edge last) until it is
   fixed; 1 when it is */
static int rm_relocate_face(rmsh* w, rm_bad const* t, unsigned int c, int large)
{
	qaws_hm const* m = w->m;
	unsigned int e = rm_longest(m, m->cf[c]), vs[4], i;
	vs[0] = m->cv[c];
	vs[1] = m->cv[hm_next(m, c)];
	vs[2] = m->cv[hm_prev(m, c)];
	vs[3] = m->co[e] != HM_NONE ? m->cv[hm_prev(m, m->co[e])] : HM_NONE;
	for (i = 0; i < 4 && !w->failed; i++)
		if (vs[i] != HM_NONE && rm_relocate(w, vs[i]) && rm_fixed(w, t, large))
			return 1;
	return 0;
}

/* c's edge flipped under the bounds' rule; 1 when it was */
static int rm_flip_bounded(rmsh* w, unsigned int c)
{
	qaws_hm* m = w->m;
	rm_use(w, m->cv[c], hm_to(m, c));
	if (!rm_flip_ok(w, c) || !rm_flip_better(w, c))
		return 0;
	rm_touch(w, m->cv[c]); rm_touch(w, hm_to(m, c)); rm_touch(w, m->cv[hm_prev(m, c)]); rm_touch(w, m->cv[hm_prev(m, m->co[c])]);
	if (!hm_flip(m, c))
	{
		w->failed = 1;
		return 0;
	}
	w->rep.flips++;
	return 1;
}

/* a face with a corner past the upper bound: the edge across the corner flipped (cutting the corner in two), else its
   corners moved, else its longest edge split, else the edge across the corner's vertex in a neighbouring face of its sector;
   1 when the face got fixed or changed */
static int rm_large(rmsh* w, rm_bad const* t)
{
	qaws_hm* m = w->m;
	unsigned int f = rm_face_of(m, t), apex, nf;
	double ex, th;
	int need;
	if (f == HM_NONE || (apex = rm_worst(w, f, 1, &ex)) == HM_NONE)
		return 0;
	if (rm_flip_bounded(w, hm_next(m, apex)) || w->failed || rm_relocate_face(w, t, apex, 1) || w->failed)
		return 1;
	if ((f = rm_face_of(m, t)) == HM_NONE || (apex = rm_worst(w, f, 1, &ex)) == HM_NONE)
		return 1;
	/* (a sector with fewer faces than its angle over the upper bound cannot be inside: it must gain one) */
	th = rm_sector(m, apex, &nf);
	need = th > 0 && nf < (unsigned int)ceil(th / w->hi - 1e-9);
	if (rm_split_large(w, rm_longest(m, f), m->cv[apex], need && rm_longest(m, f) == hm_next(m, apex)) || w->failed || (f = rm_face_of(m, t)) == HM_NONE
		|| (apex = rm_worst(w, f, 1, &ex)) == HM_NONE)
		return 1;
	/* else the edge across the apex in a neighbouring face of its sector (a sector of too few faces gets one more either way) */
	{
		unsigned int v = m->cv[apex], i, e[2][2], n = 0, c;
		for (i = 0; i < 2; i++)
		{
			unsigned int s = i ? hm_cw(m, apex) : hm_ccw(m, apex);
			if (s != HM_NONE && !(m->cflag[i ? apex : hm_prev(m, apex)] & HM_FEATURE))
			{
				e[n][0] = m->cv[hm_next(m, s)];
				e[n][1] = m->cv[hm_prev(m, s)];
				n++;
			}
		}
		/* (by vertices: a split undone leaves new corner slots) */
		for (i = 0; i < n && !w->failed; i++)
			if ((c = hm_find(m, e[i][0], e[i][1])) != HM_NONE && rm_split_large(w, c, v, need))
				return 1;
	}
	return 0;
}

/* a face with a corner under the lower bound: the edge across it collapsed, else the face's longest edge flipped, else its
   corners moved; 1 when the face got fixed or changed */
static int rm_small(rmsh* w, rm_bad const* t)
{
	qaws_hm* m = w->m;
	unsigned int f = rm_face_of(m, t), c0, e, x, y;
	double ex;
	if (f == HM_NONE || (c0 = rm_worst(w, f, 0, &ex)) == HM_NONE)
		return 0;
	e = hm_next(m, c0);
	rm_pick(m, e, &x, &y);
	if (rm_try_collapse(w, e, x, y) || w->failed || rm_try_collapse(w, e, y, x) || w->failed || rm_flip_bounded(w, rm_longest(m, f)) || w->failed)
		return 1;
	return rm_relocate_face(w, t, c0, 0);
}

/* the changed vertices' two-rings relaxed, each move held to the bounds */
static void rm_smooth_changed(rmsh* w)
{
	qaws_hm* m = w->m;
	unsigned char* mk = (unsigned char*)qaws_internal_scratch_zalloc((size_t)m->nv + 1, 1);
	unsigned int *list = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), ((size_t)m->nv + 1))), n = 0, i, k, depth, start = 0;
	if (!mk || !list)
	{
		qaws_internal_free(mk);
		qaws_internal_free(list);
		w->failed = 1;
		return;
	}
	for (i = 0; i < w->nchg; i++)
	{
		unsigned int v = w->chg[i];
		if (v < m->nv && !(m->vflag[v] & HM_DELETED) && !mk[v])
		{
			mk[v] = 1;
			list[n++] = v;
		}
	}
	for (depth = 0; depth < 2; depth++)
	{
		unsigned int end = n;
		for (i = start; i < end; i++)
		{
			unsigned int ring[64], nr = hm_ring(m, list[i], ring, 64);
			for (k = 0; k < nr && k < 64; k++)
				if (!mk[ring[k]])
				{
					mk[ring[k]] = 1;
					list[n++] = ring[k];
				}
		}
		start = end;
	}
	w->nchg = 0;
	for (i = 0; i < n; i++)
		rm_smooth_vertex(w, list[i]);
	qaws_internal_free(mk);
	qaws_internal_free(list);
}

/* rounds of: the worst large angles split away, valence flips and relaxation, the worst small angles collapsed, flipped or
   moved away, valence flips and relaxation; each operation kept only when its star gets better; until every corner not
   exempt is inside the bounds, or the worst and the count stop falling for 3 rounds */
static void rm_angles(rmsh* w)
{
	rm_bad* list = NULL;
	unsigned int cap = 0, n, i, k, done, tries, round, stall = 0, nbest = 0xFFFFFFFFu, bad;
	double best = HUGE_VAL, worst;
	w->bounded = 1;
	for (round = 0; round < 50 && !w->failed; round++)
	{
		hm_compact(w->m);
		worst = rm_census(w, &bad);
		if (!bad)
			break;
		/* (progress: the worst a tenth of a degree better, or a twentieth fewer out) */
		if (worst < best - 0.1 * RM_DEG || bad + nbest / 20 < nbest)
		{
			stall = 0;
			best = fmin(best, worst);
			nbest = bad < nbest ? bad : nbest;
		}
		else if (++stall >= 3)
			break;
		k = bad / 5 > 1 ? bad / 5 : 1;
		/* (k faces fixed in each half, worst first, out of 2k + 16 tried: the ones that resist do not hold the others
		   up, nor cost more than their share; with 2k alone the last few stalled, 13 left out on the cube at 0.6) */
		n = rm_gather(w, 1, &list, &cap);
		for (i = 0, done = 0, tries = 0; i < n && done < k && tries < 2 * k + 16 && !w->failed; i++)
			if (!rm_fixed(w, &list[i], 1))
			{
				tries++;
				done += (unsigned int)rm_large(w, &list[i]);
			}
		rm_flip(w);
		rm_smooth_changed(w);
		hm_compact(w->m);
		n = rm_gather(w, 0, &list, &cap);
		for (i = 0, done = 0, tries = 0; i < n && done < k && tries < 2 * k + 16 && !w->failed; i++)
			if (!rm_fixed(w, &list[i], 0))
			{
				tries++;
				done += (unsigned int)rm_small(w, &list[i]);
			}
		rm_flip(w);
		rm_smooth_changed(w);
	}
	w->bounded = 0;
	w->nchg = 0;
	qaws_internal_free(list);
}

/* after the loop: the distance guard made good on every face (and strict from here on), the Delaunay flips, the
   angle-bounded pass */
static void rm_post(rmsh* w, qaws_remesh_desc const* desc)
{
	if (!w->ready || w->failed)
		return;
	hm_compact(w->m);
	w->strict = 1;
	rm_tolerance(w);
	rm_delaunay(w);
	hm_compact(w->m);
	w->lo = desc->min_angle > 0 ? (double)desc->min_angle : 0;
	w->hi = desc->max_angle > 0 ? (double)desc->max_angle : 0;
	if (w->lo > 0 || w->hi > 0)
		rm_angles(w);
}

/* ---- set-up and the end ------------------------------------------------- */

/* the working mesh: the input's faces fanned, its features and the border as features, footpoints */
static qaws_status rm_begin(rmsh* w, qaws_mesh_3d const* in, qaws_mesh_3d** out)
{
	unsigned int nv = qaws_mesh_3d_get_vertex_count(in), nt = qaws_mesh_3d_get_triangle_count(in), nfe = qaws_mesh_3d_get_feature_edge_count(in), i, k;
	unsigned int *tri = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (3 * (size_t)nt + 1))), *fe = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (2 * (size_t)nfe + 1)));
	qaws_hm* m;
	qaws_status st;
	if (!tri || !fe)
	{
		qaws_internal_free(tri);
		qaws_internal_free(fe);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	qaws_mesh_3d_get_triangles(tri, in);
	qaws_mesh_3d_get_feature_edges(fe, in);
	st = hm_mesh_3d_create(in->hm.x, nv, tri, NULL, nt, out);
	qaws_internal_free(tri);
	if (st != QAWS_STATUS_OK)
	{
		qaws_internal_free(fe);
		return st;
	}
	m = &(*out)->hm;
	w->m = m;
	for (i = 0; i < nfe; i++)
	{
		unsigned int c = hm_find_edge(m, fe[2 * i], fe[2 * i + 1]);
		if (c == HM_NONE)
			continue;
		m->cflag[c] |= HM_FEATURE;
		if (m->co[c] != HM_NONE)
			m->cflag[m->co[c]] |= HM_FEATURE;
	}
	qaws_internal_free(fe);
	for (i = 0; i < m->nc; i++)
		if (m->cf[i] != HM_NONE && m->co[i] == HM_NONE)
			m->cflag[i] |= HM_FEATURE;   /* (the border) */
	/* corners: the input's, where feature lines meet or end, and where one turns by more than 60 degrees */
	for (i = 0; i < m->nv; i++)
	{
		unsigned int fd = hm_feature_degree(m, i), f[2], nf = 0, c, last = HM_NONE;
		if ((i < in->hm.nv && (in->hm.vflag[i] & HM_CORNER)) || (fd != 0 && fd != 2))
		{
			m->vflag[i] |= HM_CORNER;
			continue;
		}
		if (fd != 2)
			continue;
		for (c = m->vc[i]; c != HM_NONE; c = hm_round(m, i, c))
		{
			if ((m->cflag[c] & HM_FEATURE) && nf < 2)
				f[nf++] = hm_to(m, c);
			last = c;
		}
		if (nf < 2 && last != HM_NONE && hm_is_border_vertex(m, i) && (m->cflag[hm_prev(m, last)] & HM_FEATURE))
			f[nf++] = m->cv[hm_prev(m, last)];
		if (nf == 2)
		{
			double u[3], v[3], uu = 0, vv = 0, uv = 0;
			for (k = 0; k < 3; k++)
			{
				u[k] = rm_x(m, i)[k] - rm_x(m, f[0])[k];
				v[k] = rm_x(m, f[1])[k] - rm_x(m, i)[k];
				uu += u[k] * u[k];
				vv += v[k] * v[k];
				uv += u[k] * v[k];
			}
			if (uv < 0.5 * sqrt(uu * vv))
				m->vflag[i] |= HM_CORNER;
		}
	}
	if (!hm_foot_enable(m) || !mr_scratch_init(&w->s, w->r->nt))
		return QAWS_STATUS_ALLOCATION_FAILURE;
	w->sa = hm_attr_add(m, "qaws.size", QAWS_MESH_VERTEX, QAWS_ATTRIBUTE_LINEAR, 1);
	w->ra = hm_attr_add(m, "qaws.region", QAWS_MESH_FACE, QAWS_ATTRIBUTE_NEAREST, 1);
	w->ma = w->metric ? hm_attr_add(m, "qaws.metric", QAWS_MESH_VERTEX, QAWS_ATTRIBUTE_LINEAR, 6) : 0;
	if (w->sa == HM_NONE || w->ra == HM_NONE || w->ma == HM_NONE)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	/* each face's region: its own triangle's when the reference is the input's fan, else under its centre */
	for (i = 0; i < m->nf; i++)
	{
		unsigned int reg = HM_NONE;
		if (w->r->mesh->hm.nv == in->hm.nv && w->r->nt == m->nf && w->r->tv[3 * i] == m->cv[m->fc[i]] && w->r->tv[3 * i + 1] == m->cv[m->fc[i] + 1]
			&& w->r->tv[3 * i + 2] == m->cv[m->fc[i] + 2])
			reg = w->r->treg[i];
		else
		{
			double cen[3] = { 0, 0, 0 };
			qaws_mesh_location l;
			unsigned int c;
			for (c = m->fc[i]; c < m->fc[i] + 3; c++)
				for (k = 0; k < 3; k++)
					cen[k] += rm_x(m, m->cv[c])[k] / 3;
			mr_project(w->r, cen, HM_NONE, &w->s, &l);
			reg = l.patch < w->r->nt ? w->r->treg[l.patch] : 0;
		}
		hm_attr_at(m, w->ra, i)[0] = reg;
	}
	/* footpoints: the input's when they are on this reference (they evaluate to the vertex: a generated mesh's are on its
	   surfaces, its faces and (u, v)), else projected (feature vertices onto the creases) */
	for (i = 0; i < m->nv; i++)
	{
		qaws_mesh_location l;
		double x[3];
		unsigned int src = m->vsrc[i];
		if (in->hm.foot && src < in->hm.nv && in->hm.foot[src].patch < w->r->nt
			&& (mr_eval(w->r, &in->hm.foot[src], x, NULL), fabs(x[0] - rm_x(m, i)[0]) + fabs(x[1] - rm_x(m, i)[1]) + fabs(x[2] - rm_x(m, i)[2]) <= 1e-6 * w->r->extent))
			l = in->hm.foot[src];
		else if (hm_feature_degree(m, i) && mr_project_crease(w->r, rm_x(m, i), HM_NONE, &w->s, &l) != HUGE_VAL && l.patch != HM_NONE)
		{
		}
		else if (m->vc[i] != HM_NONE)
			mr_project_region(w->r, rm_x(m, i), HM_NONE, rm_region(w, m->cf[m->vc[i]]), &w->s, &l);
		else
			mr_project(w->r, rm_x(m, i), HM_NONE, &w->s, &l);
		if (l.patch == HM_NONE)
			return QAWS_STATUS_INVALID_ARGUMENT;
		mr_eval(w->r, &l, x, NULL);
		for (k = 0; k < 3; k++)
			rm_x(m, i)[k] = x[k];
		m->foot[i] = l;
		rm_set_size(w, i);
	}
	w->ready = 1;
	return QAWS_STATUS_OK;
}

/* barycentric coordinates of p in triangle t's flat corners, clamped into it */
static void rm_bary(qaws_mesh_reference const* r, unsigned int t, double const* p, double* a, double* b)
{
	double const *P0 = r->patch[t].p[0], *P1 = r->patch[t].p[1], *P2 = r->patch[t].p[2];
	double d00 = 0, d01 = 0, d11 = 0, d20 = 0, d21 = 0, den, u, v;
	unsigned int k;
	for (k = 0; k < 3; k++)
	{
		double e1 = P1[k] - P0[k], e2 = P2[k] - P0[k], q = p[k] - P0[k];
		d00 += e1 * e1; d01 += e1 * e2; d11 += e2 * e2; d20 += q * e1; d21 += q * e2;
	}
	den = d00 * d11 - d01 * d01;
	u = den > 0 ? (d11 * d20 - d01 * d21) / den : 1.0 / 3;
	v = den > 0 ? (d00 * d21 - d01 * d20) / den : 1.0 / 3;
	if (u < 0) u = 0;
	if (v < 0) v = 0;
	if (u + v > 1)
	{
		double s = u + v;
		u /= s;
		v /= s;
	}
	*a = u;
	*b = v;
}

/* the reference mesh's attributes onto m by its footpoints and its faces' centre locations fl */
qaws_status mx_transfer(qaws_hm* m, qaws_mesh_reference const* r, qaws_mesh_location const* fl)
{
	qaws_hm const* src = &r->mesh->hm;
	unsigned int f, k, v, c;
	for (k = 0; k < src->na; k++)
	{
		hm_attr const* at = &src->attr[k];
		unsigned int a = hm_attr_add(m, at->name, at->element, at->rule, at->n);
		qaws_scalar vals[16];
		if (a == HM_NONE)
		{
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}
		if (at->element == QAWS_MESH_VERTEX)
			for (v = 0; v < m->nv; v++)
			{
				unsigned int q;
				qaws_mesh_reference_evaluate_attribute(vals, r, &m->foot[v], 1, k);
				for (q = 0; q < at->n; q++)
					hm_attr_at(m, a, v)[q] = (double)vals[q];
			}
		else
			for (f = 0; f < m->nf; f++)
			{
				unsigned int q;
				if (fl[f].patch == HM_NONE)
					continue;
				if (at->element == QAWS_MESH_FACE)
				{
					qaws_mesh_reference_evaluate_attribute(vals, r, &fl[f], 1, k);
					for (q = 0; q < at->n; q++)
						hm_attr_at(m, a, f)[q] = (double)vals[q];
					continue;
				}
				/* corners: the vertex in the reference triangle under the face's centre */
				for (c = m->fc[f]; c < m->fc[f] + m->fn[f]; c++)
				{
					qaws_mesh_location l;
					double ba, bb;
					rm_bary(r, fl[f].patch, rm_x(m, m->cv[c]), &ba, &bb);
					l.patch = fl[f].patch;
					l.a = (qaws_scalar)ba;
					l.b = (qaws_scalar)bb;
					qaws_mesh_reference_evaluate_attribute(vals, r, &l, 1, k);
					for (q = 0; q < at->n; q++)
						hm_attr_at(m, a, c)[q] = (double)vals[q];
				}
			}
	}
	return QAWS_STATUS_OK;
}

/* the reference mesh's attributes onto the result; the faces' deviations when report */
static qaws_status rm_end(rmsh* w, int report)
{
	qaws_hm* m = w->m;
	qaws_mesh_location* fl;
	unsigned int f, k, c;
	double dmax = 0, dsum = 0;
	int need = report;
	hm_compact(m);
	fl = (qaws_mesh_location*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_mesh_location), (m->nf + 1)));
	if (!fl)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	/* each face's centre on the reference, from its first vertex's side: the deviation reported, and where face and corner
	   attributes are read; neither wanted (INSTANT's dense mesh), none projected */
	for (k = 0; k < w->r->mesh->hm.na && !need; k++)
		need = w->r->mesh->hm.attr[k].element != QAWS_MESH_VERTEX;
	for (f = 0; !need && f < m->nf; f++)
		fl[f].patch = HM_NONE;
	for (f = 0; need && f < m->nf; f++)
	{
		double cen[3] = { 0, 0, 0 }, d;
		for (c = m->fc[f]; c < m->fc[f] + m->fn[f]; c++)
			for (k = 0; k < 3; k++)
				cen[k] += rm_x(m, m->cv[c])[k] / m->fn[f];
		{
			unsigned int reg = rm_region(w, f), h = m->foot[m->cv[m->fc[f]]].patch;
			d = sqrt(mr_project_region(w->r, cen, h < w->r->nt && w->r->treg[h] == reg ? h : HM_NONE, reg, &w->s, &fl[f]));
		}
		dmax = d > dmax ? d : dmax;
		dsum += d;
	}
	if (w->metric)
		hm_attr_remove(m, w->ma);   /* (added last: its removal leaves the others where they are) */
	hm_attr_remove(m, w->ra > w->sa ? w->ra : w->sa);
	hm_attr_remove(m, w->ra > w->sa ? w->sa : w->ra);
	w->rep.max_deviation = (qaws_scalar)dmax;
	w->rep.mean_deviation = (qaws_scalar)(m->nf ? dsum / m->nf : 0);
	{
		qaws_status st = mx_transfer(m, w->r, fl);
		qaws_internal_free(fl);
		return st;
	}
}

qaws_status qaws_mesh_remesh(qaws_mesh_3d** out_mesh, qaws_mesh_report* out_report, qaws_mesh_3d const* mesh, qaws_remesh_desc const* desc)
{
	rmsh w;
	qaws_mesh_reference* own = NULL;
	qaws_mesh_3d* out = NULL;
	qaws_status st;
	unsigned int it, iterations, relax;
	if (!out_mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_mesh = NULL;
	if (!mesh || !desc || (unsigned int)desc->method > QAWS_REMESH_FIELD_ALIGNED || desc->target_length < 0)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (desc->method == QAWS_REMESH_FIELD_ALIGNED || (desc->metric && desc->method == QAWS_REMESH_CVT))
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	/* the angle bounds about 60 degrees, below one another; the tolerance not negative */
	if (desc->min_angle < 0 || desc->max_angle < 0 || desc->tolerance < 0 || (desc->min_angle > 0 && (double)desc->min_angle >= RM_PI / 3)
		|| (desc->max_angle > 0 && (double)desc->max_angle <= RM_PI / 3) || (desc->min_angle > 0 && desc->max_angle > 0 && desc->min_angle >= desc->max_angle))
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (it = 0; it < desc->pinned_count; it++)
		if (!desc->pinned || desc->pinned[it] >= mesh->hm.nv)
			return QAWS_STATUS_INVALID_ARGUMENT;
	memset(&w, 0, sizeof w);
	w.r = desc->reference;
	w.tol = desc->tolerance > 0 ? (double)desc->tolerance : 0;
	if (!w.r)
	{
		st = qaws_mesh_reference_create(&own, mesh, NULL);
		if (st != QAWS_STATUS_OK)
			return st;
		w.r = own;
	}
	if (desc->size && qaws_size_field_get_reference(desc->size) && qaws_size_field_get_reference(desc->size) != w.r)
	{
		qaws_mesh_reference_destroy(own);
		return QAWS_STATUS_INVALID_ARGUMENT;   /* (a field sampled on another reference) */
	}
	if (desc->metric && desc->metric->ref != w.r)
	{
		qaws_mesh_reference_destroy(own);
		return QAWS_STATUS_INVALID_ARGUMENT;
	}
	w.metric = desc->metric;
	if (!w.r->nt)
	{
		qaws_mesh_reference_destroy(own);
		return qaws_mesh_3d_clone(out_mesh, mesh);
	}
	w.size = desc->size;
	w.L = desc->target_length > 0 ? (double)desc->target_length : 0;
	if (w.L == 0)
	{
		qaws_mesh_quality q;
		qaws_mesh_3d_measure_quality(&q, mesh);
		w.L = q.mean_edge > 0 ? (double)q.mean_edge : 1;
	}
	/* a tolerance shrinks the sizes where the surface bends, as Mmgs's hausd does: an edge's chord gap at most 3/4 of it */
	if (w.tol > 0)
	{
		qaws_curvature_size_desc cd;
		memset(&cd, 0, sizeof cd);
		cd.chord_tolerance = (qaws_scalar)(0.75 * w.tol);
		cd.max_length = (qaws_scalar)w.r->extent;
		st = qaws_size_field_create_curvature(&w.tsize, w.r, &cd);
		if (st != QAWS_STATUS_OK)
		{
			qaws_mesh_reference_destroy(own);
			return st;
		}
	}
	/* (three relaxation sweeps in the last third only for graded or bounded meshes: INSTANT's dense mesh does better with one) */
	relax = desc->size || desc->metric || w.tol > 0 || desc->min_angle > 0 || desc->max_angle > 0 ? 3 : 1;
	iterations = desc->method == QAWS_REMESH_CVT ? 10 : desc->iterations ? desc->iterations : desc->metric ? 30 : 10;
	st = rm_begin(&w, mesh, &out);
	/* pinned vertices: corners (never moved or collapsed away) */
	for (it = 0; st == QAWS_STATUS_OK && it < desc->pinned_count; it++)
		w.m->vflag[desc->pinned[it]] |= HM_CORNER;
	for (it = 0; st == QAWS_STATUS_OK && it < iterations && !w.failed; it++)
	{
		unsigned int ops = w.rep.splits + w.rep.collapses + w.rep.flips, sweep, flips = ~0u;
		rm_split(&w);
		hm_compact(w.m);
		rm_collapse(&w);
		hm_compact(w.m);
		/* (flip passes while they find flips, two at most) */
		for (sweep = 0; sweep < 2 && flips != w.rep.flips; sweep++)
		{
			flips = w.rep.flips;
			rm_flip(&w);
		}
		hm_compact(w.m);
		w.moved = 0;
		for (sweep = 3 * it >= 2 * iterations ? relax : 1; sweep > 0; sweep--)
			rm_smooth(&w);
		rm_caps(&w);
		rm_needles(&w);
		hm_compact(w.m);
		w.rep.iterations++;
		/* converged: no edit, no move farther than a thousandth of the target */
		if (ops == w.rep.splits + w.rep.collapses + w.rep.flips && w.moved <= 1e-3 * w.L)
			break;
	}
	if (st == QAWS_STATUS_OK && !w.failed && desc->method == QAWS_REMESH_CVT && !rv_optimize(w.m, w.r, w.size, w.L, w.sa, w.ra, &w.s, desc->iterations ? desc->iterations : 30, &w.rep))
		w.failed = 1;
	rm_post(&w, desc);
	if (w.failed && st == QAWS_STATUS_OK)
		st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (st == QAWS_STATUS_OK)
		st = rm_end(&w, out_report != NULL);
	mr_scratch_free(&w.s);
	qaws_internal_free(w.ep);
	qaws_internal_free(w.chg);
	qaws_size_field_destroy(w.tsize);
	qaws_mesh_reference_destroy(own);
	if (st != QAWS_STATUS_OK)
	{
		qaws_mesh_3d_destroy(out);
		return st;
	}
	if (out_report)
	{
		qaws_mesh_3d_measure_quality(&w.rep.quality, out);
		*out_report = w.rep;
	}
	*out_mesh = out;
	return QAWS_STATUS_OK;
}
