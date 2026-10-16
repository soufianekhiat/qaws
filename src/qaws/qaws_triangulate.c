/*
 * Triangulating regions of curves: boundary sampled on the true curves,
 * constrained Delaunay triangulation on a lattice with exact predicates,
 * Ruppert refinement with boundary splits on the curves, remeshing.
 */

#include "qaws_triangulate.h"
#include "qaws_clip.h"
#include "qaws_curve_batch.h"
#include "qaws_curve.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "qaws_platform.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_wide.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define TM_NONE 0xFFFFFFFFu
#define TM_PI 3.14159265358979323846
#define TM_LATTICE_BITS 48      /* the region's extent spans 2^48 lattice steps */
#define TM_SUPER 3              /* the first three vertices are the super triangle */

typedef struct tm_vert
{
	double p[2];                /* on the lattice */
	int64_t q[2];               /* lattice coordinates */
	unsigned int curve;         /* boundary vertex: its curve, else TM_NONE */
	double t;                   /* its parameter */
	unsigned int tri;           /* a triangle using it */
} tm_vert;

typedef struct tm_tri
{
	unsigned int v[3];          /* counter-clockwise */
	unsigned int n[3];          /* neighbour across the edge opposite v[i] */
	unsigned int s[3];          /* that edge's segment + 1, 0 = free */
	unsigned char alive, in, mark;
} tm_tri;

typedef struct tm_seg
{
	unsigned int a, b;          /* the region on the left from a to b */
	unsigned int curve;
	double t0, t1;              /* parameters of a and b on the curve */
	unsigned char alive;
	unsigned char kind;          /* 0 boundary, 1 constraint (inside: no side change across it) */
} tm_seg;

struct qaws_mesh_2d
{
	qaws_clip_result* region;   /* owns the boundary curves */
	qaws_clip_result* cons;     /* owns the constraint pieces */
	unsigned int nbcurves;      /* curves [0, nbcurves) bound the region, the rest are constraints */
	qaws_curve const** curves;
	unsigned int ncurves;
	tm_vert* v;
	unsigned int nv, capv;
	tm_tri* t;
	unsigned int nt, capt;
	unsigned int* dead;         /* free triangle slots */
	unsigned int ndead, capdead;
	tm_seg* sg;
	unsigned int ns, caps;
	double origin[2], step, ext, tol;
	unsigned int last;          /* where walks start */
	int failed;
	/* the output, compacted */
	double* ov;
	unsigned int onv;
	unsigned int* otri;
	unsigned int ont;
	unsigned int* ovc;          /* output vertex: its boundary curve or TM_NONE */
	double* ovt;                /* and its parameter there */
	unsigned int* obe;          /* feature edges: the boundary (region on the left), then the constraints */
	unsigned int* obc;          /* their curves */
	double* obt;                /* their parameter intervals on them */
	unsigned char* obk;         /* 0 boundary, 1 constraint */
	unsigned int onbe, once;
	unsigned int* oquad;        /* quadrilaterals (after quadrangulation) */
	unsigned int onq;
	/* cavity scratch */
	unsigned int* stack;
	unsigned int capstack;
	unsigned int* cav;
	unsigned int ncav, capcav;
};

/* ======================================================================== */
/*  Lattice and exact predicates                                            */
/* ======================================================================== */

static void tm_snap(qaws_mesh_2d const* m, double x, double y, tm_vert* v)
{
	double a = floor((x - m->origin[0]) / m->step + 0.5), b = floor((y - m->origin[1]) / m->step + 0.5);
	v->q[0] = (int64_t)a;
	v->q[1] = (int64_t)b;
	v->p[0] = m->origin[0] + a * m->step;
	v->p[1] = m->origin[1] + b * m->step;
}

/* sign of (b - a) x (c - a): positive counter-clockwise */
static int tm_orient(qaws_mesh_2d const* m, unsigned int a, unsigned int b, unsigned int c)
{
	int64_t const* A = m->v[a].q;
	int64_t const* B = m->v[b].q;
	int64_t const* C = m->v[c].q;
	int64_t bx = B[0] - A[0], by = B[1] - A[1], cx = C[0] - A[0], cy = C[1] - A[1];
	double l = (double)bx * (double)cy, r = (double)by * (double)cx, d = l - r;
	if (fabs(d) > 4e-16 * (fabs(l) + fabs(r)))
		return d > 0 ? 1 : -1;
	{
		qaws_wide w = qaws_wide_cross(bx, by, cx, cy);
		return qaws_wide_sign(&w);
	}
}

/* sign of the in-circle determinant: positive when d lies inside the circle
   through a, b, c (counter-clockwise) */
static int tm_incircle(qaws_mesh_2d const* m, unsigned int a, unsigned int b, unsigned int c, unsigned int d)
{
	int64_t const* D = m->v[d].q;
	int64_t adx = m->v[a].q[0] - D[0], ady = m->v[a].q[1] - D[1];
	int64_t bdx = m->v[b].q[0] - D[0], bdy = m->v[b].q[1] - D[1];
	int64_t cdx = m->v[c].q[0] - D[0], cdy = m->v[c].q[1] - D[1];
	double fadx = (double)adx, fady = (double)ady, fbdx = (double)bdx, fbdy = (double)bdy, fcdx = (double)cdx, fcdy = (double)cdy;
	double bc = fbdx * fcdy - fcdx * fbdy, ca = fcdx * fady - fadx * fcdy, ab = fadx * fbdy - fbdx * fady;
	double al = fadx * fadx + fady * fady, bl = fbdx * fbdx + fbdy * fbdy, cl = fcdx * fcdx + fcdy * fcdy;
	double det = al * bc + bl * ca + cl * ab;
	double perm = al * (fabs(fbdx * fcdy) + fabs(fcdx * fbdy)) + bl * (fabs(fcdx * fady) + fabs(fadx * fcdy)) + cl * (fabs(fadx * fbdy) + fabs(fbdx * fady));
	if (fabs(det) > 4e-15 * perm)
		return det > 0 ? 1 : -1;
	{
		qaws_wide wbc = qaws_wide_cross(bdx, bdy, cdx, cdy), wca = qaws_wide_cross(cdx, cdy, adx, ady), wab = qaws_wide_cross(adx, ady, bdx, bdy);
		qaws_wide wal = qaws_wide_dot(adx, ady, adx, ady), wbl = qaws_wide_dot(bdx, bdy, bdx, bdy), wcl = qaws_wide_dot(cdx, cdy, cdx, cdy);
		qaws_wide x = qaws_wide_mul(&wal, &wbc), y = qaws_wide_mul(&wbl, &wca), z = qaws_wide_mul(&wcl, &wab);
		x = qaws_wide_add(x, &y);
		x = qaws_wide_add(x, &z);
		return qaws_wide_sign(&x);
	}
}

/* ======================================================================== */
/*  Storage                                                                 */
/* ======================================================================== */

static unsigned int tm_vert_add(qaws_mesh_2d* m, double x, double y, unsigned int curve, double t)
{
	tm_vert* v;
	if (m->nv == m->capv)
	{
		unsigned int nc = m->capv ? 2 * m->capv : 1024;
		tm_vert* g = (tm_vert*)realloc(m->v, sizeof(tm_vert) * nc);
		if (!g) { m->failed = 1; return TM_NONE; }
		m->v = g;
		m->capv = nc;
	}
	v = &m->v[m->nv];
	memset(v, 0, sizeof(*v));
	tm_snap(m, x, y, v);
	v->curve = curve;
	v->t = t;
	v->tri = TM_NONE;
	return m->nv++;
}

static unsigned int tm_tri_new(qaws_mesh_2d* m)
{
	unsigned int i;
	if (m->ndead)
		i = m->dead[--m->ndead];
	else
	{
		if (m->nt == m->capt)
		{
			unsigned int nc = m->capt ? 2 * m->capt : 2048;
			tm_tri* g = (tm_tri*)realloc(m->t, sizeof(tm_tri) * nc);
			if (!g) { m->failed = 1; return TM_NONE; }
			m->t = g;
			m->capt = nc;
		}
		i = m->nt++;
	}
	memset(&m->t[i], 0, sizeof(tm_tri));
	m->t[i].alive = 1;
	m->t[i].n[0] = m->t[i].n[1] = m->t[i].n[2] = TM_NONE;
	return i;
}

static void tm_tri_kill(qaws_mesh_2d* m, unsigned int i)
{
	m->t[i].alive = 0;
	if (m->ndead == m->capdead)
	{
		unsigned int nc = m->capdead ? 2 * m->capdead : 256;
		unsigned int* g = (unsigned int*)realloc(m->dead, sizeof(unsigned int) * nc);
		if (!g) { m->failed = 1; return; }
		m->dead = g;
		m->capdead = nc;
	}
	m->dead[m->ndead++] = i;
}

static unsigned int tm_seg_add(qaws_mesh_2d* m, unsigned int a, unsigned int b, unsigned int curve, double t0, double t1)
{
	if (m->ns == m->caps)
	{
		unsigned int nc = m->caps ? 2 * m->caps : 256;
		tm_seg* g = (tm_seg*)realloc(m->sg, sizeof(tm_seg) * nc);
		if (!g) { m->failed = 1; return TM_NONE; }
		m->sg = g;
		m->caps = nc;
	}
	m->sg[m->ns].a = a; m->sg[m->ns].b = b; m->sg[m->ns].curve = curve;
	m->sg[m->ns].t0 = t0; m->sg[m->ns].t1 = t1; m->sg[m->ns].alive = 1; m->sg[m->ns].kind = 0;
	return m->ns++;
}

static int tm_push(unsigned int** buf, unsigned int* n, unsigned int* cap, unsigned int x)
{
	if (*n == *cap)
	{
		unsigned int nc = *cap ? 2 * *cap : 64;
		unsigned int* g = (unsigned int*)realloc(*buf, sizeof(unsigned int) * nc);
		if (!g) return 0;
		*buf = g;
		*cap = nc;
	}
	(*buf)[(*n)++] = x;
	return 1;
}

/* the edge index of t whose neighbour is u */
static unsigned int tm_back(qaws_mesh_2d const* m, unsigned int t, unsigned int u)
{
	unsigned int k;
	for (k = 0; k < 3; k++)
		if (m->t[t].n[k] == u) return k;
	return TM_NONE;
}

/* ======================================================================== */
/*  Location and insertion                                                  */
/* ======================================================================== */

/* the triangle containing vertex p (walking); *on receives the edge p lies
   on (TM_NONE if inside), *dup an existing vertex at p (TM_NONE if none).
   With stop_at_seg, a walk that must cross a boundary edge stops there and
   returns TM_NONE with *seg_tri / *seg_edge set. */
static unsigned int tm_locate(qaws_mesh_2d* m, unsigned int p, unsigned int start, int stop_at_seg,
	unsigned int* dup, unsigned int* seg_tri, unsigned int* seg_edge)
{
	unsigned int t = start, steps = 0, r = 0;
	*dup = TM_NONE;
	if (t == TM_NONE || t >= m->nt || !m->t[t].alive) t = m->last;
	if (t == TM_NONE || t >= m->nt || !m->t[t].alive)
	{
		for (t = 0; t < m->nt && !m->t[t].alive; t++) {}
	}
	while (steps++ < 4 * m->nt + 16)
	{
		tm_tri const* T = &m->t[t];
		unsigned int k, moved = 0, zeros = 0, zk = 0;
		r = (r + 1) % 3;
		for (k = 0; k < 3; k++)
		{
			unsigned int e = (k + r) % 3, a = T->v[(e + 1) % 3], b = T->v[(e + 2) % 3];
			int o = tm_orient(m, a, b, p);
			if (o < 0)
			{
				if (stop_at_seg && T->s[e])
				{
					*seg_tri = t;
					*seg_edge = e;
					return TM_NONE;
				}
				if (T->n[e] == TM_NONE) return TM_NONE;
				t = T->n[e];
				moved = 1;
				break;
			}
			if (o == 0) { zeros++; zk = e; }
		}
		if (moved) continue;
		if (zeros >= 2)
		{
			/* on two edges: at their common vertex */
			for (k = 0; k < 3; k++)
				if (m->v[T->v[k]].q[0] == m->v[p].q[0] && m->v[T->v[k]].q[1] == m->v[p].q[1])
					*dup = T->v[k];
		}
		(void)zk;
		return t;
	}
	return TM_NONE;
}

/* Bowyer-Watson: the triangles whose circumcircles hold p, grown from t0
   without crossing boundary edges, replaced by the fan round p. Returns 0 on
   failure. The new triangles take the side flag `in` of the cavity's first
   triangle. */
static int tm_insert(qaws_mesh_2d* m, unsigned int p, unsigned int t0)
{
	unsigned int ns = 0, i, k, nb = 0, nnew, base_in = m->t[t0].in;
	unsigned int* bnd = NULL;   /* pairs: cavity triangle, edge */
	unsigned int capb = 0;
	m->ncav = 0;
	if (!tm_push(&m->stack, &ns, &m->capstack, t0)) { m->failed = 1; return 0; }
	m->t[t0].mark = 1;
	if (!tm_push(&m->cav, &m->ncav, &m->capcav, t0)) { m->failed = 1; return 0; }
	while (ns)
	{
		unsigned int t = m->stack[--ns];
		for (k = 0; k < 3; k++)
		{
			unsigned int u = m->t[t].n[k];
			if (u == TM_NONE || m->t[u].mark || m->t[t].s[k]) continue;
			if (tm_incircle(m, m->t[u].v[0], m->t[u].v[1], m->t[u].v[2], p) > 0)
			{
				m->t[u].mark = 1;
				if (!tm_push(&m->stack, &ns, &m->capstack, u) || !tm_push(&m->cav, &m->ncav, &m->capcav, u)) { m->failed = 1; return 0; }
			}
		}
	}
	/* the cavity's boundary, oriented counter-clockwise round p */
	for (i = 0; i < m->ncav; i++)
	{
		unsigned int t = m->cav[i];
		for (k = 0; k < 3; k++)
		{
			unsigned int u = m->t[t].n[k];
			if (u != TM_NONE && m->t[u].mark) continue;
			/* a boundary edge that p lies on or behind would make a flat or
			   inverted triangle: the cavity is not star-shaped from p */
			if (tm_orient(m, m->t[t].v[(k + 1) % 3], m->t[t].v[(k + 2) % 3], p) <= 0)
			{
				unsigned int j;
				for (j = 0; j < m->ncav; j++) m->t[m->cav[j]].mark = 0;
				free(bnd);
				return 0;
			}
			if (!tm_push(&bnd, &nb, &capb, t) || !tm_push(&bnd, &nb, &capb, k)) { m->failed = 1; free(bnd); return 0; }
		}
	}
	nnew = nb / 2;
	{
		unsigned int* nt = (unsigned int*)malloc(sizeof(unsigned int) * (nnew + 1));
		unsigned int* ea = (unsigned int*)malloc(sizeof(unsigned int) * (nnew + 1));
		unsigned int* eb = (unsigned int*)malloc(sizeof(unsigned int) * (nnew + 1));
		unsigned int* on = (unsigned int*)malloc(sizeof(unsigned int) * (nnew + 1));
		unsigned int* os = (unsigned int*)malloc(sizeof(unsigned int) * (nnew + 1));
		unsigned int* ob = (unsigned int*)malloc(sizeof(unsigned int) * (nnew + 1));
		if (!nt || !ea || !eb || !on || !os || !ob) { free(nt); free(ea); free(eb); free(on); free(os); free(ob); free(bnd); m->failed = 1; return 0; }
		for (i = 0; i < nnew; i++)
		{
			unsigned int t = bnd[2 * i], e = bnd[2 * i + 1];
			ea[i] = m->t[t].v[(e + 1) % 3];
			eb[i] = m->t[t].v[(e + 2) % 3];
			on[i] = m->t[t].n[e];
			os[i] = m->t[t].s[e];
			/* the outside neighbour's edge back, while the cavity is alive (its
			   slots are reused at once) */
			ob[i] = on[i] == TM_NONE ? 0 : tm_back(m, on[i], t);
		}
		for (i = 0; i < m->ncav; i++)
		{
			m->t[m->cav[i]].mark = 0;
			tm_tri_kill(m, m->cav[i]);
		}
		for (i = 0; i < nnew; i++)
		{
			unsigned int T = tm_tri_new(m);
			if (T == TM_NONE) break;
			nt[i] = T;
			m->t[T].v[0] = p; m->t[T].v[1] = ea[i]; m->t[T].v[2] = eb[i];
			m->t[T].s[0] = os[i];
			m->t[T].in = (unsigned char)base_in;
			m->t[T].n[0] = on[i];
			if (on[i] != TM_NONE) m->t[on[i]].n[ob[i]] = T;
			m->v[ea[i]].tri = T;
			m->v[eb[i]].tri = T;
		}
		if (!m->failed)
		{
			/* the fan's own edges: (b_i, p) meets the triangle starting at b_i */
			for (i = 0; i < nnew; i++)
			{
				unsigned int j;
				for (j = 0; j < nnew; j++)
					if (ea[j] == eb[i]) { m->t[nt[i]].n[1] = nt[j]; m->t[nt[j]].n[2] = nt[i]; break; }
			}
			m->v[p].tri = nt[0];
			m->last = nt[0];
		}
		free(nt); free(ea); free(eb); free(on); free(os); free(ob);
	}
	free(bnd);
	return !m->failed;
}

/* insert vertex p (located from start); returns the vertex now at p (an
   existing one for a duplicate), TM_NONE on failure */
static unsigned int tm_add_point(qaws_mesh_2d* m, unsigned int p, unsigned int start)
{
	unsigned int dup, st, se;
	unsigned int t = tm_locate(m, p, start, 0, &dup, &st, &se);
	if (t == TM_NONE) return TM_NONE;
	if (dup != TM_NONE) return dup;
	if (!tm_insert(m, p, t)) return TM_NONE;
	return p;
}

/* ======================================================================== */
/*  Boundary edges                                                          */
/* ======================================================================== */

/* the triangle and edge index holding edge a -> b (with a, b in its
   counter-clockwise order), by the fan of a */
static int tm_find_edge(qaws_mesh_2d const* m, unsigned int a, unsigned int b, unsigned int* out_t, unsigned int* out_k)
{
	unsigned int t0 = m->v[a].tri, t = t0, guard = 0;
	if (t == TM_NONE || !m->t[t].alive)
	{
		/* a stale hint: search (rare) */
		for (t = 0; t < m->nt; t++)
			if (m->t[t].alive && (m->t[t].v[0] == a || m->t[t].v[1] == a || m->t[t].v[2] == a)) break;
		if (t == m->nt) return 0;
		t0 = t;
	}
	/* round a one way, then the other (a boundary of the super triangle
	   cannot be met: a is inside it) */
	do
	{
		unsigned int i, k;
		for (i = 0; i < 3; i++) if (m->t[t].v[i] == a) break;
		if (i == 3) return 0;
		/* edge a -> v[i+1] is opposite v[i+2] */
		if (m->t[t].v[(i + 1) % 3] == b) { *out_t = t; *out_k = (i + 2) % 3; return 1; }
		k = (i + 1) % 3;      /* across edge (v[i+2], a) is opposite v[i+1] */
		t = m->t[t].n[k];
		if (t == TM_NONE) break;
	} while (t != t0 && ++guard < 1000000);
	return 0;
}

static void tm_mark_seg(qaws_mesh_2d* m, unsigned int t, unsigned int k, unsigned int seg)
{
	unsigned int u = m->t[t].n[k];
	m->t[t].s[k] = seg + 1;
	if (u != TM_NONE)
	{
		unsigned int j = tm_back(m, u, t);
		if (j != TM_NONE) m->t[u].s[j] = seg + 1;
	}
}

/* link new triangles nt[0..nn) to each other and to the cavity's boundary
   records rec (a, b, outer, outer edge, segment + 1: edge a -> b seen from
   inside) */
static void tm_relink(qaws_mesh_2d* m, unsigned int const* nt, unsigned int nn, unsigned int const* rec, unsigned int nrec)
{
	unsigned int i, k, r, j;
	for (i = 0; i < nn; i++)
	{
		tm_tri* T = &m->t[nt[i]];
		for (k = 0; k < 3; k++)
		{
			unsigned int va = T->v[(k + 1) % 3], vb = T->v[(k + 2) % 3], found = 0;
			m->v[va].tri = nt[i];
			for (r = 0; r < nrec && !found; r++)
				if (rec[5 * r] == va && rec[5 * r + 1] == vb)
				{
					T->n[k] = rec[5 * r + 2];
					T->s[k] = rec[5 * r + 4];
					if (rec[5 * r + 2] != TM_NONE) m->t[rec[5 * r + 2]].n[rec[5 * r + 3]] = nt[i];
					found = 1;
				}
			for (j = 0; j < nn && !found; j++)
			{
				unsigned int q;
				if (j == i) continue;
				for (q = 0; q < 3; q++)
					if (m->t[nt[j]].v[(q + 1) % 3] == vb && m->t[nt[j]].v[(q + 2) % 3] == va)
					{
						T->n[k] = nt[j];
						found = 1;
						break;
					}
			}
		}
	}
}

/* the constrained Delaunay triangulation of the pseudo-polygon first, the
   chain ch[0..n) (left of first -> last), last (Anglada) */
static int tm_pseudo(qaws_mesh_2d* m, unsigned int first, unsigned int last, unsigned int const* ch, unsigned int n,
	unsigned int** out, unsigned int* nout, unsigned int* cap)
{
	unsigned int c = 0, i, T;
	if (n == 0) return 1;
	for (i = 1; i < n; i++)
		if (tm_incircle(m, first, last, ch[c], ch[i]) > 0) c = i;
	T = tm_tri_new(m);
	if (T == TM_NONE) return 0;
	m->t[T].v[0] = first; m->t[T].v[1] = last; m->t[T].v[2] = ch[c];
	if (!tm_push(out, nout, cap, T)) return 0;
	return tm_pseudo(m, first, ch[c], ch, c, out, nout, cap) && tm_pseudo(m, ch[c], last, ch + c + 1, n - c - 1, out, nout, cap);
}

static int tm_recover(qaws_mesh_2d* m, unsigned int seg);

/* split segment seg at vertex x lying on it, both halves recovered */
static int tm_split_at(qaws_mesh_2d* m, unsigned int seg, unsigned int x)
{
	tm_seg s = m->sg[seg];
	double la = hypot(m->v[x].p[0] - m->v[s.a].p[0], m->v[x].p[1] - m->v[s.a].p[1]);
	double lb = hypot(m->v[s.b].p[0] - m->v[s.a].p[0], m->v[s.b].p[1] - m->v[s.a].p[1]);
	double tx = s.t0 + (lb > 0 ? la / lb : 0.5) * (s.t1 - s.t0);
	unsigned int s2;
	m->sg[seg].b = x;
	m->sg[seg].t1 = tx;
	s2 = tm_seg_add(m, x, s.b, s.curve, tx, s.t1);
	if (s2 == TM_NONE) return 0;
	m->sg[s2].kind = s.kind;
	return tm_recover(m, seg) && tm_recover(m, s2);
}

/* segments s and o cross: a vertex at the crossing of their chords, the
   boundary edge of o freed, both split there and recovered */
static int tm_cross_split(qaws_mesh_2d* m, unsigned int s, unsigned int o)
{
	tm_seg S = m->sg[s], O = m->sg[o];
	double const* A = m->v[S.a].p;
	double const* B = m->v[S.b].p;
	double const* P = m->v[O.a].p;
	double const* Q = m->v[O.b].p;
	double d = (B[0] - A[0]) * (Q[1] - P[1]) - (B[1] - A[1]) * (Q[0] - P[0]), u, w, x[2];
	unsigned int t, k, v, vx, s2, o2;
	if (!(fabs(d) > 0)) return 0;
	u = ((P[0] - A[0]) * (Q[1] - P[1]) - (P[1] - A[1]) * (Q[0] - P[0])) / d;
	w = ((P[0] - A[0]) * (B[1] - A[1]) - (P[1] - A[1]) * (B[0] - A[0])) / d;
	if (u < 0) u = 0;
	if (u > 1) u = 1;
	if (w < 0) w = 0;
	if (w > 1) w = 1;
	x[0] = A[0] + u * (B[0] - A[0]);
	x[1] = A[1] + u * (B[1] - A[1]);
	/* onto both curves: Newton on C_s(ts) = C_o(to) from the chords' crossing
	   (refinement splits next to it on the curves, which must stay on their
	   sides of each other) */
	{
		qaws_curve const* cs = m->curves[S.curve];
		qaws_curve const* co = m->curves[O.curve];
		double ts = S.t0 + u * (S.t1 - S.t0), to = O.t0 + w * (O.t1 - O.t0);
		double slo = S.t0 < S.t1 ? S.t0 : S.t1, shi = S.t0 < S.t1 ? S.t1 : S.t0;
		double olo = O.t0 < O.t1 ? O.t0 : O.t1, ohi = O.t0 < O.t1 ? O.t1 : O.t0;
		unsigned int it;
		for (it = 0; it < 30; it++)
		{
			qaws_eval_result_2d es, eo;
			double fx, fy, det, dts, dto;
			qaws_curve_evaluate_2d(&es, cs, (qaws_scalar)ts, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1);
			qaws_curve_evaluate_2d(&eo, co, (qaws_scalar)to, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1);
			fx = (double)es.position.x - eo.position.x;
			fy = (double)es.position.y - eo.position.y;
			/* [Cs' -Co'] (dts, dto) = -f */
			det = (double)es.d1.x * -(double)eo.d1.y - (double)es.d1.y * -(double)eo.d1.x;
			if (!(fabs(det) > 0)) break;
			dts = (-fx * -(double)eo.d1.y - -fy * -(double)eo.d1.x) / det;
			dto = ((double)es.d1.x * -fy - (double)es.d1.y * -fx) / det;
			ts += dts; to += dto;
			if (ts < slo) ts = slo;
			if (ts > shi) ts = shi;
			if (to < olo) to = olo;
			if (to > ohi) to = ohi;
			if (fabs(dts) + fabs(dto) < 1e-15 * (1 + fabs(ts) + fabs(to))) break;
		}
		{
			qaws_eval_result_2d es;
			qaws_curve_evaluate_2d(&es, cs, (qaws_scalar)ts, QAWS_EVAL_FLAG_POSITION);
			if (hypot((double)es.position.x - x[0], (double)es.position.y - x[1]) < 4 * hypot(B[0] - A[0], B[1] - A[1]) + 4 * hypot(Q[0] - P[0], Q[1] - P[1]))
			{
				x[0] = es.position.x; x[1] = es.position.y;
				u = S.t1 != S.t0 ? (ts - S.t0) / (S.t1 - S.t0) : 0.5;
				w = O.t1 != O.t0 ? (to - O.t0) / (O.t1 - O.t0) : 0.5;
			}
		}
	}
	/* the crossed edge freed so the point can go in across it */
	if (tm_find_edge(m, O.a, O.b, &t, &k) || tm_find_edge(m, O.b, O.a, &t, &k))
	{
		unsigned int nb = m->t[t].n[k];
		m->t[t].s[k] = 0;
		if (nb != TM_NONE) { unsigned int j = tm_back(m, nb, t); if (j != TM_NONE) m->t[nb].s[j] = 0; }
	}
	/* the vertex takes the constraint's curve (a boundary one if o bounds the region) */
	v = tm_vert_add(m, x[0], x[1], O.kind == 0 ? O.curve : S.curve, O.kind == 0 ? O.t0 + w * (O.t1 - O.t0) : S.t0 + u * (S.t1 - S.t0));
	if (v == TM_NONE) return 0;
	vx = tm_add_point(m, v, m->v[S.a].tri);
	if (vx == TM_NONE) return 0;
	if (vx != v) m->nv--;
	/* split both (a crossing at an end splits only the other) */
	if (vx != S.a && vx != S.b)
	{
		m->sg[s].b = vx; m->sg[s].t1 = S.t0 + u * (S.t1 - S.t0);
		s2 = tm_seg_add(m, vx, S.b, S.curve, m->sg[s].t1, S.t1);
		if (s2 == TM_NONE) return 0;
		m->sg[s2].kind = S.kind;
	}
	else
		s2 = TM_NONE;
	if (vx != O.a && vx != O.b)
	{
		m->sg[o].b = vx; m->sg[o].t1 = O.t0 + w * (O.t1 - O.t0);
		o2 = tm_seg_add(m, vx, O.b, O.curve, m->sg[o].t1, O.t1);
		if (o2 == TM_NONE) return 0;
		m->sg[o2].kind = O.kind;
	}
	else
		o2 = TM_NONE;
	return tm_recover(m, o) && (o2 == TM_NONE || tm_recover(m, o2)) && tm_recover(m, s) && (s2 == TM_NONE || tm_recover(m, s2));
}

/* make a -> b an edge and mark it segment seg: the triangles it crosses
   removed, the two sides retriangulated (Anglada); a vertex exactly on it
   splits it there */
static int tm_recover(qaws_mesh_2d* m, unsigned int seg)
{
	unsigned int a = m->sg[seg].a, b = m->sg[seg].b, t, k;
	unsigned int *cross = NULL, *L = NULL, *R = NULL, *rec = NULL, *nt = NULL;
	unsigned int nc = 0, cc = 0, nl = 0, cl = 0, nr = 0, cr = 0, nrec = 0, crec = 0, nn = 0, cn = 0, i;
	unsigned int t0, tt, start = TM_NONE, e = 0, guard = 0;
	int ok = 0;
	if (tm_find_edge(m, a, b, &t, &k) || tm_find_edge(m, b, a, &t, &k)) { tm_mark_seg(m, t, k, seg); return 1; }
	/* the triangle round a that a -> b leaves through */
	t0 = m->v[a].tri;
	if (t0 == TM_NONE || !m->t[t0].alive) return 0;
	tt = t0;
	do
	{
		unsigned int x, y;
		int ox, oy;
		for (i = 0; i < 3; i++) if (m->t[tt].v[i] == a) break;
		if (i == 3) return 0;
		x = m->t[tt].v[(i + 1) % 3]; y = m->t[tt].v[(i + 2) % 3];
		ox = tm_orient(m, a, b, x); oy = tm_orient(m, a, b, y);
		if (ox == 0 && x >= TM_SUPER && (m->v[x].q[0] - m->v[a].q[0]) * (double)(m->v[b].q[0] - m->v[a].q[0]) + (m->v[x].q[1] - m->v[a].q[1]) * (double)(m->v[b].q[1] - m->v[a].q[1]) > 0)
			return tm_split_at(m, seg, x);
		if (ox < 0 && oy > 0) { start = tt; e = i; break; }
		tt = m->t[tt].n[(i + 1) % 3];
	} while (tt != TM_NONE && tt != t0 && ++guard < 1000000);
	if (start == TM_NONE) return 0;
	/* walk to b: the crossed triangles, the vertices left and right */
	tt = start;
	if (!tm_push(&cross, &nc, &cc, tt)) goto done;
	if (!tm_push(&R, &nr, &cr, m->t[tt].v[(e + 1) % 3]) || !tm_push(&L, &nl, &cl, m->t[tt].v[(e + 2) % 3])) goto done;
	for (guard = 0; guard < 4 * m->nt + 16; guard++)
	{
		unsigned int u, j, z, lv, rv, p, q;
		int o;
		if (m->t[tt].s[e])
		{
			/* another segment crosses this one: constraints meet there (a new
			   vertex at the crossing of the two chords, both split); two
			   boundary edges never cross (the region was cleaned) */
			unsigned int other = m->t[tt].s[e] - 1;
			if (m->sg[other].kind == 0 && m->sg[seg].kind == 0) goto done;
			free(cross); free(L); free(R);
			return tm_cross_split(m, seg, other);
		}
		u = m->t[tt].n[e];
		if (u == TM_NONE) goto done;
		j = tm_back(m, u, tt);
		z = m->t[u].v[j];
		if (!tm_push(&cross, &nc, &cc, u)) goto done;
		if (z == b) { ok = 1; break; }
		o = tm_orient(m, a, b, z);
		if (o == 0)
		{
			free(cross); free(L); free(R);
			return tm_split_at(m, seg, z);
		}
		p = m->t[u].v[(j + 1) % 3]; q = m->t[u].v[(j + 2) % 3];
		if (tm_orient(m, a, b, p) > 0) { lv = p; rv = q; } else { lv = q; rv = p; }
		(void)rv;
		if (o > 0)
		{
			if (!tm_push(&L, &nl, &cl, z)) goto done;
			/* out through (rv, z): opposite lv */
			for (e = 0; e < 3; e++) if (m->t[u].v[e] == lv) break;
		}
		else
		{
			if (!tm_push(&R, &nr, &cr, z)) goto done;
			for (e = 0; e < 3; e++) if (m->t[u].v[e] == rv) break;
		}
		tt = u;
	}
	if (!ok) goto done;
	ok = 0;
	/* the cavity's boundary records, then the crossed triangles removed */
	for (i = 0; i < nc; i++) m->t[cross[i]].mark = 1;
	for (i = 0; i < nc; i++)
	{
		tm_tri const* T = &m->t[cross[i]];
		for (k = 0; k < 3; k++)
		{
			unsigned int u = T->n[k];
			if (u != TM_NONE && m->t[u].mark) continue;
			if (!tm_push(&rec, &nrec, &crec, T->v[(k + 1) % 3]) || !tm_push(&rec, &nrec, &crec, T->v[(k + 2) % 3])
				|| !tm_push(&rec, &nrec, &crec, u) || !tm_push(&rec, &nrec, &crec, u == TM_NONE ? 0 : tm_back(m, u, cross[i]))
				|| !tm_push(&rec, &nrec, &crec, T->s[k])) goto done;
		}
	}
	for (i = 0; i < nc; i++) { m->t[cross[i]].mark = 0; tm_tri_kill(m, cross[i]); }
	/* the two sides: a, L, b on the left; b, R reversed, a on the right */
	for (i = 0; i < nr / 2; i++) { unsigned int w = R[i]; R[i] = R[nr - 1 - i]; R[nr - 1 - i] = w; }
	if (!tm_pseudo(m, a, b, L, nl, &nt, &nn, &cn) || !tm_pseudo(m, b, a, R, nr, &nt, &nn, &cn)) goto done;
	tm_relink(m, nt, nn, rec, nrec / 5);
	if (!tm_find_edge(m, a, b, &t, &k) && !tm_find_edge(m, b, a, &t, &k)) goto done;
	tm_mark_seg(m, t, k, seg);
	m->last = nt[0];
	ok = 1;
done:
	free(cross); free(L); free(R); free(rec); free(nt);
	if (!ok && !m->failed) return 0;
	return ok;
}

/* inside = odd number of boundary edges crossed from the super triangle */
static int tm_classify(qaws_mesh_2d* m)
{
	unsigned int* st = NULL;
	unsigned int ns = 0, cap = 0, i;
	for (i = 0; i < m->nt; i++) m->t[i].mark = 0;
	for (i = 0; i < m->nt; i++)
	{
		tm_tri* T = &m->t[i];
		if (!T->alive || T->mark) continue;
		if (T->v[0] < TM_SUPER || T->v[1] < TM_SUPER || T->v[2] < TM_SUPER)
		{
			unsigned int k;
			T->mark = 1;
			T->in = 0;
			if (!tm_push(&st, &ns, &cap, i)) { free(st); return 0; }
			while (ns)
			{
				unsigned int c = st[--ns];
				for (k = 0; k < 3; k++)
				{
					unsigned int u = m->t[c].n[k];
					if (u == TM_NONE || m->t[u].mark) continue;
					m->t[u].mark = 1;
					m->t[u].in = (unsigned char)(m->t[c].in ^ (m->t[c].s[k] && m->sg[m->t[c].s[k] - 1].kind == 0));
					if (!tm_push(&st, &ns, &cap, u)) { free(st); return 0; }
				}
			}
		}
	}
	for (i = 0; i < m->nt; i++) m->t[i].mark = 0;
	free(st);
	return 1;
}

/* ======================================================================== */
/*  Boundary sampling                                                       */
/* ======================================================================== */

static void tm_eval(qaws_curve const* c, double t, double* p)
{
	qaws_eval_result_2d r;
	qaws_curve_evaluate_2d(&r, c, (qaws_scalar)t, QAWS_EVAL_FLAG_POSITION);
	p[0] = r.position.x;
	p[1] = r.position.y;
}

static double tm_seg_dist(double const* q, double const* a, double const* b)
{
	double dx = b[0] - a[0], dy = b[1] - a[1], l2 = dx * dx + dy * dy, u;
	u = l2 > 0 ? ((q[0] - a[0]) * dx + (q[1] - a[1]) * dy) / l2 : 0;
	if (u < 0) u = 0;
	if (u > 1) u = 1;
	return hypot(q[0] - a[0] - u * dx, q[1] - a[1] - u * dy);
}

typedef struct tm_samples
{
	double* t;
	unsigned int n, cap;
} tm_samples;

static int tm_sample_push(tm_samples* s, double t)
{
	if (s->n == s->cap)
	{
		unsigned int nc = s->cap ? 2 * s->cap : 64;
		double* g = (double*)realloc(s->t, sizeof(double) * nc);
		if (!g) return 0;
		s->t = g;
		s->cap = nc;
	}
	s->t[s->n++] = t;
	return 1;
}

/* the parameters after t0 up to t1 (t1 included) where the chords stay
   within the tolerance (and the size field) */
static int tm_sample(qaws_mesh_2d const* m, qaws_curve const* c, double t0, double const* p0, double t1, double const* p1,
	qaws_triangulate_desc const* d, tm_samples* s, unsigned int depth)
{
	double tm = 0.5 * (t0 + t1), pa[2], pm[2], pb[2], dev, len = hypot(p1[0] - p0[0], p1[1] - p0[1]);
	int split;
	tm_eval(c, 0.5 * (t0 + tm), pa);
	tm_eval(c, tm, pm);
	tm_eval(c, 0.5 * (tm + t1), pb);
	dev = tm_seg_dist(pm, p0, p1);
	if (tm_seg_dist(pa, p0, p1) > dev) dev = tm_seg_dist(pa, p0, p1);
	if (tm_seg_dist(pb, p0, p1) > dev) dev = tm_seg_dist(pb, p0, p1);
	/* a closed curve (its chord a point) is split before its deviation counts */
	split = dev > m->tol || (depth < 2 && len <= 1e-9 * m->ext);
	if (!split && d->size_fn)
	{
		qaws_vec2 q;
		q.x = (qaws_scalar)pm[0]; q.y = (qaws_scalar)pm[1];
		split = len > (double)d->size_fn(d->size_user, q);
	}
	if (split && depth < 30 && (len > 4 * m->step || depth < 2))
	{
		return tm_sample(m, c, t0, p0, tm, pm, d, s, depth + 1) && tm_sample(m, c, tm, pm, t1, p1, d, s, depth + 1);
	}
	return tm_sample_push(s, t1);
}

/* a sample forced onto a boundary curve: where a constraint piece ends on
   it (the two share the point, so they share the vertex) */
typedef struct tm_force
{
	unsigned int curve;
	double t, p[2];
} tm_force;

static void tm_point(qaws_curve const* c, unsigned int ci, double t, tm_force const* f, unsigned int nf, double* q)
{
	unsigned int i;
	for (i = 0; i < nf; i++)
		if (f[i].curve == ci && f[i].t == t) { q[0] = f[i].p[0]; q[1] = f[i].p[1]; return; }
	tm_eval(c, t, q);
}

static int tm_cmp_double(void const* a, void const* b)
{
	double x = *(double const*)a, y = *(double const*)b;
	return x < y ? -1 : (x > y ? 1 : 0);
}

/* the sample parameters of curve c (index ci) after its start, its end
   last: span by span (corners of polylines, knots of splines) and through
   the forced parameters */
static int tm_curve_samples(qaws_mesh_2d const* m, qaws_triangulate_desc const* d, qaws_curve const* c, unsigned int ci,
	tm_force const* f, unsigned int nf, tm_samples* s)
{
	double t0 = c->parameter_range.min_value, t1 = c->parameter_range.max_value, *br, pa[2], ta;
	unsigned int nb = 0, i, sp;
	br = (double*)malloc(sizeof(double) * (c->span_count + nf + 2));
	if (!br) return 0;
	for (sp = 1; sp < c->span_count; sp++)
		if (c->span_boundaries[sp] > t0 && c->span_boundaries[sp] < t1) br[nb++] = c->span_boundaries[sp];
	for (i = 0; i < nf; i++)
		if (f[i].curve == ci && f[i].t > t0 && f[i].t < t1) br[nb++] = f[i].t;
	qsort(br, nb, sizeof(double), tm_cmp_double);
	br[nb++] = t1;
	s->n = 0;
	tm_point(c, ci, t0, f, nf, pa);
	ta = t0;
	for (i = 0; i < nb; i++)
	{
		double pb[2];
		if (!(br[i] > ta)) continue;
		tm_point(c, ci, br[i], f, nf, pb);
		if (!tm_sample(m, c, ta, pa, br[i], pb, d, s, 0)) { free(br); return 0; }
		ta = br[i];
		pa[0] = pb[0]; pa[1] = pb[1];
	}
	free(br);
	return 1;
}

/* the vertices along the curves of path p (curve indices from ci0) and the
   segments of the kind between them; a closed path's last sample is its
   first */
static int tm_chain(qaws_mesh_2d* m, qaws_triangulate_desc const* d, qaws_path_2d const* p, unsigned int ci0, int closed,
	unsigned char kind, tm_force const* f, unsigned int nf, tm_samples* s)
{
	unsigned int first = TM_NONE, prev = TM_NONE, prev_curve = 0, k, j;
	double prev_t = 0;
	for (k = 0; k < p->curve_count; k++)
	{
		qaws_curve const* c = p->curves[k];
		unsigned int ci = ci0 + k;
		double t0 = c->parameter_range.min_value;
		if (!tm_curve_samples(m, d, c, ci, f, nf, s)) return 0;
		if (prev == TM_NONE)
		{
			double q[2];
			unsigned int v, w;
			tm_point(c, ci, t0, f, nf, q);
			v = tm_vert_add(m, q[0], q[1], ci, t0);
			if (v == TM_NONE) return 0;
			w = tm_add_point(m, v, m->last);
			if (w == TM_NONE) return 0;
			if (w != v) m->nv--;
			first = prev = w;
		}
		else if (kind == 0)
		{
			/* the junction starts this curve */
			m->v[prev].curve = ci;
			m->v[prev].t = t0;
		}
		prev_curve = ci;
		prev_t = t0;
		for (j = 0; j < s->n; j++)
		{
			unsigned int v, w;
			if (closed && j + 1 == s->n && k + 1 == p->curve_count)
				w = first;
			else
			{
				double q[2];
				tm_point(c, ci, s->t[j], f, nf, q);
				v = tm_vert_add(m, q[0], q[1], ci, s->t[j]);
				if (v == TM_NONE) return 0;
				w = tm_add_point(m, v, m->v[prev].tri);
				if (w == TM_NONE) return 0;
				if (w != v) m->nv--;
			}
			if (w != prev)
			{
				unsigned int sg = tm_seg_add(m, prev, w, prev_curve, prev_t, s->t[j]);
				if (sg == TM_NONE) return 0;
				m->sg[sg].kind = kind;
				prev = w;
				prev_t = s->t[j];
			}
		}
	}
	return 1;
}

/* ======================================================================== */
/*  Refinement                                                              */
/* ======================================================================== */

static double tm_len(qaws_mesh_2d const* m, unsigned int a, unsigned int b)
{
	return hypot(m->v[b].p[0] - m->v[a].p[0], m->v[b].p[1] - m->v[a].p[1]);
}

/* a vertex inside the diametral circle of segment a b */
static int tm_encroaches(qaws_mesh_2d const* m, double const* p, unsigned int a, unsigned int b)
{
	double const* A = m->v[a].p;
	double const* B = m->v[b].p;
	return (p[0] - A[0]) * (p[0] - B[0]) + (p[1] - A[1]) * (p[1] - B[1]) < 0;
}

/* split segment seg at its parameter midpoint on its curve */
static int tm_split_seg(qaws_mesh_2d* m, unsigned int seg)
{
	tm_seg s = m->sg[seg];
	double tmid = 0.5 * (s.t0 + s.t1), p[2];
	unsigned int t, k, v, w, s1, s2;
	tm_eval(m->curves[s.curve], tmid, p);
	/* the boundary edge freed, then the point inserted */
	if (tm_find_edge(m, s.a, s.b, &t, &k) || tm_find_edge(m, s.b, s.a, &t, &k))
	{
		unsigned int u = m->t[t].n[k];
		m->t[t].s[k] = 0;
		if (u != TM_NONE) { unsigned int j = tm_back(m, u, t); if (j != TM_NONE) m->t[u].s[j] = 0; }
	}
	v = tm_vert_add(m, p[0], p[1], s.curve, tmid);
	if (v == TM_NONE) return 0;
	if ((m->v[v].q[0] == m->v[s.a].q[0] && m->v[v].q[1] == m->v[s.a].q[1]) || (m->v[v].q[0] == m->v[s.b].q[0] && m->v[v].q[1] == m->v[s.b].q[1]))
	{
		/* below the lattice: keep the segment */
		m->nv--;
		return tm_recover(m, seg) ? 2 : 0;
	}
	w = tm_add_point(m, v, m->v[s.a].tri);
	if (w == TM_NONE) return 0;
	m->sg[seg].alive = 0;
	s1 = tm_seg_add(m, s.a, w, s.curve, s.t0, tmid);
	s2 = tm_seg_add(m, w, s.b, s.curve, tmid, s.t1);
	if (s1 == TM_NONE || s2 == TM_NONE) return 0;
	m->sg[s1].kind = m->sg[s2].kind = s.kind;
	if (!tm_recover(m, s1) || !tm_recover(m, s2)) return 0;
	return tm_classify(m);
}

/* the circumcentre of triangle t */
static void tm_circumcentre(qaws_mesh_2d const* m, unsigned int t, double* c)
{
	double const* A = m->v[m->t[t].v[0]].p;
	double const* B = m->v[m->t[t].v[1]].p;
	double const* C = m->v[m->t[t].v[2]].p;
	double bx = B[0] - A[0], by = B[1] - A[1], cx = C[0] - A[0], cy = C[1] - A[1];
	double d = 2 * (bx * cy - by * cx), b2 = bx * bx + by * by, c2 = cx * cx + cy * cy;
	c[0] = A[0] + (cy * b2 - by * c2) / d;
	c[1] = A[1] + (bx * c2 - cx * b2) / d;
}

typedef struct tm_quality
{
	double sin_min;             /* sin of the angle bound */
	double max_area;
	qaws_size_fn size_fn;
	void* size_user;
	double hmin;                /* edges below it are not refined further */
} tm_quality;

/* 1: bad by angle, 2: bad by size; 0 good */
static int tm_bad(qaws_mesh_2d const* m, unsigned int t, tm_quality const* q)
{
	unsigned int const* V = m->t[t].v;
	double l0 = tm_len(m, V[1], V[2]), l1 = tm_len(m, V[2], V[0]), l2 = tm_len(m, V[0], V[1]);
	double lmin = l0 < l1 ? (l0 < l2 ? l0 : l2) : (l1 < l2 ? l1 : l2);
	double lmax = l0 > l1 ? (l0 > l2 ? l0 : l2) : (l1 > l2 ? l1 : l2);
	double A = 0.5 * fabs((m->v[V[1]].p[0] - m->v[V[0]].p[0]) * (m->v[V[2]].p[1] - m->v[V[0]].p[1]) - (m->v[V[1]].p[1] - m->v[V[0]].p[1]) * (m->v[V[2]].p[0] - m->v[V[0]].p[0]));
	if (q->max_area > 0 && A > q->max_area) return 2;
	if (q->size_fn)
	{
		qaws_vec2 c;
		c.x = (qaws_scalar)((m->v[V[0]].p[0] + m->v[V[1]].p[0] + m->v[V[2]].p[0]) / 3);
		c.y = (qaws_scalar)((m->v[V[0]].p[1] + m->v[V[1]].p[1] + m->v[V[2]].p[1]) / 3);
		if (lmax > (double)q->size_fn(q->size_user, c)) return 2;
	}
	if (q->sin_min > 0 && lmin > q->hmin)
	{
		/* the smallest angle is opposite the shortest edge: sin = lmin / 2R,
		   R = l0 l1 l2 / 4A */
		double s = A > 0 ? lmin * 2 * A / (l0 * l1 * l2) : 0;
		if (s < q->sin_min)
		{
			/* the input's own angle: both edges at the small angle are boundary
			   edges (a sharp corner of the input; no point can open it) */
			unsigned int opp = l0 == lmin ? 0u : (l1 == lmin ? 1u : 2u);
			if (m->t[t].s[(opp + 1) % 3] && m->t[t].s[(opp + 2) % 3])
				return 0;
			return 1;
		}
	}
	return 0;
}

static int tm_refine(qaws_mesh_2d* m, tm_quality const* q, unsigned int max_vertices)
{
	int changed = 1;
	unsigned int pass = 0, i;
	/* encroached segments first */
	while (changed && pass++ < 64 && m->nv < max_vertices)
	{
		changed = 0;
		for (i = 0; i < m->ns && m->nv < max_vertices; i++)
		{
			unsigned int t, k, u;
			if (!m->sg[i].alive || tm_len(m, m->sg[i].a, m->sg[i].b) < q->hmin) continue;
			if (!tm_find_edge(m, m->sg[i].a, m->sg[i].b, &t, &k) && !tm_find_edge(m, m->sg[i].b, m->sg[i].a, &t, &k)) continue;
			u = m->t[t].n[k];
			if (tm_encroaches(m, m->v[m->t[t].v[k]].p, m->sg[i].a, m->sg[i].b)
				|| (u != TM_NONE && tm_encroaches(m, m->v[m->t[u].v[tm_back(m, u, t)]].p, m->sg[i].a, m->sg[i].b)))
			{
				int r = tm_split_seg(m, i);
				if (!r) return 0;
				changed |= r == 1;
			}
		}
	}
	changed = 1;
	pass = 0;
	while (changed && pass++ < 1000 && m->nv < max_vertices)
	{
		changed = 0;
		for (i = 0; i < m->nt && m->nv < max_vertices; i++)
		{
			double c[2];
			unsigned int dup, st = TM_NONE, se = TM_NONE, loc, v, k, enc = TM_NONE;
			if (!m->t[i].alive || !m->t[i].in || !tm_bad(m, i, q)) continue;
			tm_circumcentre(m, i, c);
			if (!(fabs(c[0]) < 1e300 && fabs(c[1]) < 1e300)) continue;
			v = tm_vert_add(m, c[0], c[1], TM_NONE, 0);
			if (v == TM_NONE) return 0;
			loc = tm_locate(m, v, i, 1, &dup, &st, &se);
			if (loc == TM_NONE)
			{
				/* the walk met the boundary: that segment is encroached */
				m->nv--;
				if (st != TM_NONE && m->t[st].s[se])
				{
					int r = tm_split_seg(m, m->t[st].s[se] - 1);
					if (!r) return 0;
					changed |= r == 1;
				}
				continue;
			}
			if (dup != TM_NONE || !m->t[loc].in) { m->nv--; continue; }
			/* boundary edges round the cavity that the point encroaches */
			{
				unsigned int ns = 0, j, kk;
				m->ncav = 0;
				tm_push(&m->stack, &ns, &m->capstack, loc);
				tm_push(&m->cav, &m->ncav, &m->capcav, loc);
				m->t[loc].mark = 1;
				while (ns && enc == TM_NONE)
				{
					unsigned int t = m->stack[--ns];
					for (kk = 0; kk < 3; kk++)
					{
						unsigned int u = m->t[t].n[kk];
						if (m->t[t].s[kk])
						{
							tm_seg const* s = &m->sg[m->t[t].s[kk] - 1];
							if (tm_encroaches(m, m->v[v].p, s->a, s->b)) { enc = m->t[t].s[kk] - 1; break; }
							continue;
						}
						if (u == TM_NONE || m->t[u].mark) continue;
						if (tm_incircle(m, m->t[u].v[0], m->t[u].v[1], m->t[u].v[2], v) > 0)
						{
							m->t[u].mark = 1;
							tm_push(&m->stack, &ns, &m->capstack, u);
							tm_push(&m->cav, &m->ncav, &m->capcav, u);
						}
					}
				}
				for (j = 0; j < m->ncav; j++) m->t[m->cav[j]].mark = 0;
			}
			if (enc != TM_NONE)
			{
				int r;
				m->nv--;
				if (tm_len(m, m->sg[enc].a, m->sg[enc].b) < q->hmin) continue;
				r = tm_split_seg(m, enc);
				if (!r) return 0;
				changed |= r == 1;
				continue;
			}
			if (!tm_insert(m, v, loc)) { m->nv--; if (m->failed) return 0; continue; }
			changed = 1;
			(void)k;
		}
	}
	return !m->failed;
}

/* ======================================================================== */
/*  Output                                                                  */
/* ======================================================================== */

/* is edge k of inside triangle i a feature edge to output: a boundary edge
   (the other side outside) or a constraint edge (from the triangle where it
   runs the segment's way, so once) */
static int tm_feature(qaws_mesh_2d const* m, unsigned int i, unsigned int k, unsigned int kind)
{
	unsigned int s = m->t[i].s[k], n = m->t[i].n[k];
	tm_seg const* sg;
	if (!s) return 0;
	sg = &m->sg[s - 1];
	if (sg->kind != kind) return 0;
	if (kind == 0) return n == TM_NONE || !m->t[n].in;
	return sg->a == m->t[i].v[(k + 1) % 3];
}

static int tm_compact(qaws_mesh_2d* m)
{
	unsigned int* map = (unsigned int*)malloc(sizeof(unsigned int) * (m->nv + 1));
	unsigned int i, k, kind, nt = 0, nb = 0;
	free(m->ov); free(m->otri); free(m->ovc); free(m->ovt); free(m->obe); free(m->obc); free(m->obt); free(m->obk); free(m->oquad);
	m->ov = NULL; m->otri = NULL; m->ovc = NULL; m->ovt = NULL; m->obe = NULL; m->obc = NULL; m->obt = NULL; m->obk = NULL; m->oquad = NULL;
	m->onv = m->ont = m->onbe = m->once = m->onq = 0;
	if (!map) return 0;
	for (i = 0; i < m->nv; i++) map[i] = TM_NONE;
	for (i = 0; i < m->nt; i++)
		if (m->t[i].alive && m->t[i].in)
		{
			nt++;
			for (k = 0; k < 3; k++)
				if (tm_feature(m, i, k, 0) || tm_feature(m, i, k, 1)) nb++;
		}
	m->otri = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (nt + 1));
	m->ov = (double*)malloc(sizeof(double) * 2 * (m->nv + 1));
	m->ovc = (unsigned int*)malloc(sizeof(unsigned int) * (m->nv + 1));
	m->ovt = (double*)malloc(sizeof(double) * (m->nv + 1));
	m->obe = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (nb + 1));
	m->obc = (unsigned int*)malloc(sizeof(unsigned int) * (nb + 1));
	m->obt = (double*)malloc(sizeof(double) * 2 * (nb + 1));
	m->obk = (unsigned char*)malloc(nb + 1);
	if (!m->otri || !m->ov || !m->ovc || !m->ovt || !m->obe || !m->obc || !m->obt || !m->obk) { free(map); return 0; }
	for (i = 0; i < m->nt; i++)
	{
		if (!m->t[i].alive || !m->t[i].in) continue;
		for (k = 0; k < 3; k++)
		{
			unsigned int v = m->t[i].v[k];
			if (map[v] == TM_NONE)
			{
				map[v] = m->onv;
				m->ov[2 * m->onv] = m->v[v].p[0];
				m->ov[2 * m->onv + 1] = m->v[v].p[1];
				m->ovc[m->onv] = m->v[v].curve;
				m->ovt[m->onv] = m->v[v].t;
				m->onv++;
			}
			m->otri[3 * m->ont + k] = map[v];
		}
		m->ont++;
	}
	/* feature edges: the boundary first, then the constraints */
	for (kind = 0; kind < 2; kind++)
		for (i = 0; i < m->nt; i++)
		{
			if (!m->t[i].alive || !m->t[i].in) continue;
			for (k = 0; k < 3; k++)
				if (tm_feature(m, i, k, kind))
				{
					unsigned int e = m->onbe + m->once;
					tm_seg const* sg = &m->sg[m->t[i].s[k] - 1];
					int same = sg->a == m->t[i].v[(k + 1) % 3];
					m->obe[2 * e] = map[m->t[i].v[(k + 1) % 3]];
					m->obe[2 * e + 1] = map[m->t[i].v[(k + 2) % 3]];
					m->obc[e] = sg->curve;
					m->obt[2 * e] = same ? sg->t0 : sg->t1;
					m->obt[2 * e + 1] = same ? sg->t1 : sg->t0;
					m->obk[e] = (unsigned char)kind;
					if (kind) m->once++; else m->onbe++;
				}
		}
	free(map);
	return 1;
}

/* ======================================================================== */
/*  Public                                                                  */
/* ======================================================================== */

void qaws_mesh_2d_destroy(qaws_mesh_2d* m)
{
	if (!m) return;
	qaws_clip_result_destroy(m->region);
	qaws_clip_result_destroy(m->cons);
	free((void*)m->curves);
	free(m->v); free(m->t); free(m->dead); free(m->sg);
	free(m->ov); free(m->otri); free(m->ovc); free(m->ovt); free(m->obe); free(m->obc); free(m->obt); free(m->obk); free(m->oquad);
	free(m->stack); free(m->cav);
	free(m);
}

qaws_status qaws_triangulate_region_2d(qaws_mesh_2d** out, qaws_triangulate_desc const* d)
{
	qaws_mesh_2d* m;
	qaws_clip_desc cd;
	qaws_status st;
	unsigned int i, k, np, ncons = 0, ncc = 0, nforced = 0;
	tm_force* forced = NULL;
	double lo[2] = { 1e300, 1e300 }, hi[2] = { -1e300, -1e300 };
	if (!d || !out || (d->path_count && !d->paths) || (d->constraint_count && !d->constraints)) return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	m = (qaws_mesh_2d*)calloc(1, sizeof(qaws_mesh_2d));
	if (!m) return QAWS_STATUS_ALLOCATION_FAILURE;
	m->last = TM_NONE;
	/* the region cleaned */
	memset(&cd, 0, sizeof(cd));
	cd.subjects = d->paths;
	cd.subject_count = d->path_count;
	cd.clip_type = QAWS_CLIP_UNION;
	cd.fill_rule = d->fill_rule;
	st = qaws_clip_execute(&m->region, &cd);
	if (st != QAWS_STATUS_OK) { qaws_mesh_2d_destroy(m); return st; }
	np = qaws_clip_result_get_path_count(m->region);
	for (i = 0; i < np; i++)
	{
		qaws_path_2d p;
		qaws_vec2 a, b;
		qaws_clip_result_get_path(&p, m->region, i);
		m->ncurves += p.curve_count;
		if (qaws_path_compute_bounds_2d(&a, &b, &p) == QAWS_STATUS_OK)
		{
			if (a.x < lo[0]) lo[0] = a.x;
			if (a.y < lo[1]) lo[1] = a.y;
			if (b.x > hi[0]) hi[0] = b.x;
			if (b.y > hi[1]) hi[1] = b.y;
		}
	}
	/* the constraints: their pieces inside the region */
	if (d->constraint_count)
	{
		qaws_clip_desc kd;
		qaws_path_2d* rp = (qaws_path_2d*)malloc(sizeof(qaws_path_2d) * (np + 1));
		if (!rp) { qaws_mesh_2d_destroy(m); return QAWS_STATUS_ALLOCATION_FAILURE; }
		for (i = 0; i < np; i++) qaws_clip_result_get_path(&rp[i], m->region, i);
		memset(&kd, 0, sizeof(kd));
		kd.open_subjects = d->constraints;
		kd.open_subject_count = d->constraint_count;
		kd.clips = rp;
		kd.clip_count = np;
		kd.clip_type = QAWS_CLIP_INTERSECTION;
		kd.fill_rule = QAWS_FILL_NON_ZERO;
		st = qaws_clip_execute(&m->cons, &kd);
		free(rp);
		if (st != QAWS_STATUS_OK) { qaws_mesh_2d_destroy(m); return st; }
		ncons = qaws_clip_result_get_open_path_count(m->cons);
		for (i = 0; i < ncons; i++)
		{
			qaws_path_2d p;
			qaws_clip_result_get_open_path(&p, m->cons, i);
			ncc += p.curve_count;
		}
	}
	m->nbcurves = m->ncurves;
	m->ncurves += ncc;
	m->curves = (qaws_curve const**)malloc(sizeof(qaws_curve*) * (m->ncurves + 1));
	if (!m->curves) { qaws_mesh_2d_destroy(m); return QAWS_STATUS_ALLOCATION_FAILURE; }
	{
		unsigned int ci = 0;
		for (i = 0; i < np; i++)
		{
			qaws_path_2d p;
			qaws_clip_result_get_path(&p, m->region, i);
			for (k = 0; k < p.curve_count; k++) m->curves[ci++] = p.curves[k];
		}
		for (i = 0; i < ncons; i++)
		{
			qaws_path_2d p;
			qaws_clip_result_get_open_path(&p, m->cons, i);
			for (k = 0; k < p.curve_count; k++) m->curves[ci++] = p.curves[k];
		}
	}
	if (!np || !(hi[0] >= lo[0]))
	{
		*out = m;
		return tm_compact(m) ? QAWS_STATUS_OK : QAWS_STATUS_ALLOCATION_FAILURE;
	}
	m->ext = hi[0] - lo[0] > hi[1] - lo[1] ? hi[0] - lo[0] : hi[1] - lo[1];
	if (!(m->ext > 0)) m->ext = 1;
	m->tol = d->tolerance > 0 ? (double)d->tolerance : 1e-3 * m->ext;
	/* the lattice: a power of two step, the origin on it */
	m->step = ldexp(1.0, (int)ceil(log2(m->ext)) - TM_LATTICE_BITS);
	m->origin[0] = floor(0.5 * (lo[0] + hi[0]) / m->step) * m->step;
	m->origin[1] = floor(0.5 * (lo[1] + hi[1]) / m->step) * m->step;
	/* the super triangle round lattice box [-2^47, 2^47]^2 */
	{
		double S = ldexp(m->step, 50);
		unsigned int a = tm_vert_add(m, m->origin[0] - S, m->origin[1] - 0.5 * S, TM_NONE, 0);
		unsigned int b = tm_vert_add(m, m->origin[0] + S, m->origin[1] - 0.5 * S, TM_NONE, 0);
		unsigned int c = tm_vert_add(m, m->origin[0], m->origin[1] + S, TM_NONE, 0);
		unsigned int t = tm_tri_new(m);
		if (m->failed) { qaws_mesh_2d_destroy(m); return QAWS_STATUS_ALLOCATION_FAILURE; }
		m->t[t].v[0] = a; m->t[t].v[1] = b; m->t[t].v[2] = c;
		m->v[a].tri = m->v[b].tri = m->v[c].tri = t;
		m->last = t;
	}
	/* where constraint pieces end on the boundary: forced boundary samples */
	if (ncons && m->nbcurves)
	{
		qaws_scalar* pts = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 4 * ncons);
		qaws_closest_point* cl = (qaws_closest_point*)malloc(sizeof(qaws_closest_point) * 2 * ncons);
		qaws_closest_desc qd;
		forced = (tm_force*)malloc(sizeof(tm_force) * 2 * ncons);
		if (!pts || !cl || !forced) { free(pts); free(cl); free(forced); qaws_mesh_2d_destroy(m); return QAWS_STATUS_ALLOCATION_FAILURE; }
		for (i = 0; i < ncons; i++)
		{
			qaws_path_2d p;
			double a[2], b[2];
			qaws_clip_result_get_open_path(&p, m->cons, i);
			tm_eval(p.curves[0], p.curves[0]->parameter_range.min_value, a);
			tm_eval(p.curves[p.curve_count - 1], p.curves[p.curve_count - 1]->parameter_range.max_value, b);
			pts[4 * i] = (qaws_scalar)a[0]; pts[4 * i + 1] = (qaws_scalar)a[1];
			pts[4 * i + 2] = (qaws_scalar)b[0]; pts[4 * i + 3] = (qaws_scalar)b[1];
		}
		memset(&qd, 0, sizeof(qd));
		qd.curves = m->curves;
		qd.curve_count = m->nbcurves;
		qd.points = pts;
		qd.point_count = 2 * ncons;
		qd.max_distance = (qaws_scalar)(sizeof(qaws_scalar) == 4 ? 1e-4 * m->ext : 1e-8 * m->ext);
		if (qaws_curve_batch_find_closest(cl, NULL, &qd) == QAWS_STATUS_OK)
			for (i = 0; i < 2 * ncons; i++)
				if (cl[i].curve != QAWS_CURVE_BATCH_NONE)
				{
					forced[nforced].curve = cl[i].curve;
					forced[nforced].t = (double)cl[i].parameter;
					forced[nforced].p[0] = pts[2 * i];
					forced[nforced].p[1] = pts[2 * i + 1];
					nforced++;
				}
		free(pts);
		free(cl);
	}
	/* the boundary loops, then the constraint pieces (their ends on the
	   boundary meet the forced samples: the same points) */
	{
		unsigned int ci = 0;
		tm_samples s;
		memset(&s, 0, sizeof(s));
		for (i = 0; i < np && !m->failed; i++)
		{
			qaws_path_2d p;
			qaws_clip_result_get_path(&p, m->region, i);
			if (!tm_chain(m, d, &p, ci, 1, 0, forced, nforced, &s)) m->failed = 1;
			ci += p.curve_count;
		}
		for (i = 0; i < ncons && !m->failed; i++)
		{
			qaws_path_2d p;
			qaws_clip_result_get_open_path(&p, m->cons, i);
			if (!tm_chain(m, d, &p, ci, 0, 1, forced, nforced, &s)) m->failed = 1;
			ci += p.curve_count;
		}
		free(s.t);
	}
	free(forced);
	if (m->failed) { qaws_mesh_2d_destroy(m); return QAWS_STATUS_ALLOCATION_FAILURE; }
	for (i = 0; i < m->ns; i++)
		if (m->sg[i].alive && !tm_recover(m, i)) { qaws_mesh_2d_destroy(m); return QAWS_STATUS_DEGENERATE_CURVE; }
	if (!tm_classify(m)) { qaws_mesh_2d_destroy(m); return QAWS_STATUS_ALLOCATION_FAILURE; }
	if (d->min_angle > 0 || d->max_area > 0 || d->size_fn)
	{
		tm_quality q;
		double ang = (double)d->min_angle;
		if (ang > 33) ang = 33;
		q.sin_min = ang > 0 ? sin(ang * TM_PI / 180) : 0;
		q.max_area = (double)d->max_area;
		q.size_fn = d->size_fn;
		q.size_user = d->size_user;
		q.hmin = 1e-6 * m->ext;
		if (!tm_refine(m, &q, d->max_vertices ? d->max_vertices : 4000000u)) { qaws_mesh_2d_destroy(m); return QAWS_STATUS_ALLOCATION_FAILURE; }
	}
	if (!tm_compact(m)) { qaws_mesh_2d_destroy(m); return QAWS_STATUS_ALLOCATION_FAILURE; }
	*out = m;
	return QAWS_STATUS_OK;
}

unsigned int qaws_mesh_2d_get_vertex_count(qaws_mesh_2d const* m) { return m ? m->onv : 0; }

qaws_vec2 qaws_mesh_2d_get_vertex(qaws_mesh_2d const* m, unsigned int i)
{
	qaws_vec2 r;
	r.x = (qaws_scalar)m->ov[2 * i];
	r.y = (qaws_scalar)m->ov[2 * i + 1];
	return r;
}

unsigned int qaws_mesh_2d_get_triangle_count(qaws_mesh_2d const* m) { return m ? m->ont : 0; }

void qaws_mesh_2d_get_triangle(unsigned int out[3], qaws_mesh_2d const* m, unsigned int t)
{
	out[0] = m->otri[3 * t];
	out[1] = m->otri[3 * t + 1];
	out[2] = m->otri[3 * t + 2];
}

int qaws_mesh_2d_get_vertex_source(unsigned int* curve, qaws_scalar* t, qaws_mesh_2d const* m, unsigned int i)
{
	if (m->ovc[i] == TM_NONE) return 0;
	if (curve) *curve = m->ovc[i];
	if (t) *t = (qaws_scalar)m->ovt[i];
	return 1;
}

unsigned int qaws_mesh_2d_get_boundary_curve_count(qaws_mesh_2d const* m) { return m ? m->ncurves : 0; }
qaws_curve const* qaws_mesh_2d_get_boundary_curve(qaws_mesh_2d const* m, unsigned int c) { return m->curves[c]; }
unsigned int qaws_mesh_2d_get_boundary_edge_count(qaws_mesh_2d const* m) { return m ? m->onbe : 0; }

void qaws_mesh_2d_get_boundary_edge(unsigned int out[2], qaws_mesh_2d const* m, unsigned int e)
{
	out[0] = m->obe[2 * e];
	out[1] = m->obe[2 * e + 1];
}

/* ======================================================================== */
/*  Remeshing (Botsch-Kobbelt) on the output mesh                           */
/* ======================================================================== */

#define RM_MAX_VERTICES 4000000u   /* splitting stops there whatever the size field asks */

typedef struct rm_ctx
{
	qaws_mesh_2d* m;
	double* p;                  /* positions */
	unsigned int* vc;           /* boundary curve or TM_NONE */
	double* vt;                 /* parameter on it */
	unsigned char* corner;      /* never moved nor removed */
	unsigned int nv, capv;
	unsigned int* tri;          /* 3 per triangle, TM_NONE first = dead */
	unsigned int nt, capt;
	unsigned int* be;           /* boundary edges: a, b (region on the left) */
	unsigned int* bc;
	double* bt;                 /* t at a, t at b on curve bc */
	unsigned char* bk;          /* 0 boundary, 1 constraint */
	unsigned int nb, capb;
	double len;
	qaws_size_fn size_fn;
	void* size_user;
	/* per pass */
	uint64_t* ekey;             /* sorted: (min << 32 | max) per triangle edge */
	unsigned int* eref;         /* triangle * 3 + edge, in ekey's order */
	unsigned int ne;
	uint64_t* bkey;             /* sorted boundary keys */
	unsigned int* bref;
	unsigned int nbk;           /* their count when built (edges added since are not in it) */
	unsigned char* used;
	int failed;
} rm_ctx;

static uint64_t rm_key(unsigned int a, unsigned int b)
{
	return a < b ? ((uint64_t)a << 32) | b : ((uint64_t)b << 32) | a;
}

static double rm_target(rm_ctx const* x, double const* p)
{
	if (x->size_fn)
	{
		qaws_vec2 q;
		double v;
		q.x = (qaws_scalar)p[0]; q.y = (qaws_scalar)p[1];
		v = (double)x->size_fn(x->size_user, q);
		/* a size field at or below zero there: the fallback length */
		if (v > 0) return v;
	}
	return x->len;
}

/* twice the signed area of a, b, c (a's position replaced by pa if given) */
static double tri_cross(rm_ctx const* x, unsigned int a, unsigned int b, unsigned int c, double const* pa)
{
	double const* A = pa ? pa : x->p + 2 * a;
	double const* B = x->p + 2 * b;
	double const* C = x->p + 2 * c;
	return (B[0] - A[0]) * (C[1] - A[1]) - (B[1] - A[1]) * (C[0] - A[0]);
}

static unsigned int rm_vert(rm_ctx* x, double px, double py, unsigned int curve, double t)
{
	if (x->nv == x->capv)
	{
		unsigned int nc = x->capv ? 2 * x->capv : 256;
		double* p = (double*)realloc(x->p, sizeof(double) * 2 * nc);
		unsigned int* vc;
		double* vt;
		unsigned char* co;
		if (!p) { x->failed = 1; return TM_NONE; }
		x->p = p;
		vc = (unsigned int*)realloc(x->vc, sizeof(unsigned int) * nc);
		if (!vc) { x->failed = 1; return TM_NONE; }
		x->vc = vc;
		vt = (double*)realloc(x->vt, sizeof(double) * nc);
		if (!vt) { x->failed = 1; return TM_NONE; }
		x->vt = vt;
		co = (unsigned char*)realloc(x->corner, nc);
		if (!co) { x->failed = 1; return TM_NONE; }
		x->corner = co;
		x->capv = nc;
	}
	x->p[2 * x->nv] = px; x->p[2 * x->nv + 1] = py;
	x->vc[x->nv] = curve; x->vt[x->nv] = t; x->corner[x->nv] = 0;
	return x->nv++;
}

static unsigned int rm_tri(rm_ctx* x, unsigned int a, unsigned int b, unsigned int c)
{
	if (x->nt == x->capt)
	{
		unsigned int nc = x->capt ? 2 * x->capt : 512;
		unsigned int* g = (unsigned int*)realloc(x->tri, sizeof(unsigned int) * 3 * nc);
		unsigned char* u;
		if (!g) { x->failed = 1; return TM_NONE; }
		x->tri = g;
		u = (unsigned char*)realloc(x->used, nc);
		if (!u) { x->failed = 1; return TM_NONE; }
		memset(u + x->capt, 0, nc - x->capt);
		x->used = u;
		x->capt = nc;
	}
	x->tri[3 * x->nt] = a; x->tri[3 * x->nt + 1] = b; x->tri[3 * x->nt + 2] = c;
	x->used[x->nt] = 1;
	return x->nt++;
}

static unsigned int rm_bedge(rm_ctx* x, unsigned int a, unsigned int b, unsigned int c, double t0, double t1, unsigned char kind)
{
	if (x->nb == x->capb)
	{
		unsigned int nc = x->capb ? 2 * x->capb : 256;
		unsigned int* be = (unsigned int*)realloc(x->be, sizeof(unsigned int) * 2 * nc);
		unsigned int* bc;
		double* bt;
		if (!be) { x->failed = 1; return TM_NONE; }
		x->be = be;
		bc = (unsigned int*)realloc(x->bc, sizeof(unsigned int) * nc);
		if (!bc) { x->failed = 1; return TM_NONE; }
		x->bc = bc;
		bt = (double*)realloc(x->bt, sizeof(double) * 2 * nc);
		if (!bt) { x->failed = 1; return TM_NONE; }
		x->bt = bt;
		{
			unsigned char* bk = (unsigned char*)realloc(x->bk, nc);
			if (!bk) { x->failed = 1; return TM_NONE; }
			x->bk = bk;
		}
		x->capb = nc;
	}
	x->be[2 * x->nb] = a; x->be[2 * x->nb + 1] = b; x->bc[x->nb] = c;
	x->bt[2 * x->nb] = t0; x->bt[2 * x->nb + 1] = t1;
	x->bk[x->nb] = kind;
	return x->nb++;
}

/* keys with a payload, sorted (heap sort: no recursion, no allocation) */
static void rm_sift(uint64_t* k, unsigned int* r, unsigned int root, unsigned int n)
{
	for (;;)
	{
		unsigned int c = 2 * root + 1;
		uint64_t t;
		unsigned int u;
		if (c >= n) break;
		if (c + 1 < n && k[c + 1] > k[c]) c++;
		if (k[root] >= k[c]) break;
		t = k[root]; u = r[root]; k[root] = k[c]; r[root] = r[c]; k[c] = t; r[c] = u;
		root = c;
	}
}

static void rm_sort(uint64_t* k, unsigned int* r, unsigned int n)
{
	unsigned int start, end;
	if (n < 2) return;
	for (start = n / 2; start-- > 0;)
		rm_sift(k, r, start, n);
	for (end = n - 1; end > 0; end--)
	{
		uint64_t t = k[0];
		unsigned int u = r[0];
		k[0] = k[end]; r[0] = r[end]; k[end] = t; r[end] = u;
		rm_sift(k, r, 0, end);
	}
}

/* the edge and boundary maps of the live triangles */
static int rm_maps(rm_ctx* x)
{
	unsigned int i, k, n = 0;
	free(x->ekey); free(x->eref); free(x->bkey); free(x->bref);
	x->ekey = (uint64_t*)malloc(sizeof(uint64_t) * (3 * x->nt + 1));
	x->eref = (unsigned int*)malloc(sizeof(unsigned int) * (3 * x->nt + 1));
	x->bkey = (uint64_t*)malloc(sizeof(uint64_t) * (x->nb + 1));
	x->bref = (unsigned int*)malloc(sizeof(unsigned int) * (x->nb + 1));
	if (!x->ekey || !x->eref || !x->bkey || !x->bref) { x->failed = 1; return 0; }
	for (i = 0; i < x->nt; i++)
	{
		if (x->tri[3 * i] == TM_NONE) continue;
		for (k = 0; k < 3; k++)
		{
			x->ekey[n] = rm_key(x->tri[3 * i + (k + 1) % 3], x->tri[3 * i + (k + 2) % 3]);
			x->eref[n] = 3 * i + k;
			n++;
		}
	}
	x->ne = n;
	rm_sort(x->ekey, x->eref, n);
	for (i = 0; i < x->nb; i++)
	{
		x->bkey[i] = rm_key(x->be[2 * i], x->be[2 * i + 1]);
		x->bref[i] = i;
	}
	rm_sort(x->bkey, x->bref, x->nb);
	x->nbk = x->nb;
	return 1;
}

static unsigned int rm_find(uint64_t const* k, unsigned int n, uint64_t key)
{
	unsigned int lo = 0, hi = n;
	while (lo < hi)
	{
		unsigned int mid = (lo + hi) / 2;
		if (k[mid] < key) lo = mid + 1; else hi = mid;
	}
	return lo < n && k[lo] == key ? lo : TM_NONE;
}

/* the boundary edge a b (either way), TM_NONE if free */
static unsigned int rm_boundary(rm_ctx const* x, unsigned int a, unsigned int b)
{
	unsigned int i = rm_find(x->bkey, x->nbk, rm_key(a, b));
	return i == TM_NONE ? TM_NONE : x->bref[i];
}

/* the live triangle refs (tri * 3 + edge) on edge a b; their count */
static unsigned int rm_edge_tris(rm_ctx const* x, unsigned int a, unsigned int b, unsigned int* r)
{
	uint64_t key = rm_key(a, b);
	unsigned int i = rm_find(x->ekey, x->ne, key), n = 0;
	if (i == TM_NONE) return 0;
	while (i < x->ne && x->ekey[i] == key && n < 2)
	{
		/* every triangle that had the edge when the map was built (a caller
		   skips the edge if one of them changed since: its used flag) */
		if (x->tri[3 * (x->eref[i] / 3)] != TM_NONE)
			r[n++] = x->eref[i];
		i++;
	}
	return n;
}

static int rm_split(rm_ctx* x)
{
	unsigned int i, ne, splits = 0;
	if (!rm_maps(x)) return -1;
	ne = x->ne;
	memset(x->used, 0, x->capt);
	for (i = 0; i < ne && !x->failed && x->nv < RM_MAX_VERTICES; i++)
	{
		unsigned int r[2], n, a, b, t1, k1, bi, mv, j;
		double mid[2], tm = 0, len;
		if (i > 0 && x->ekey[i] == x->ekey[i - 1]) continue;
		t1 = x->eref[i] / 3; k1 = x->eref[i] % 3;
		if (x->tri[3 * t1] == TM_NONE || x->used[t1]) continue;
		a = x->tri[3 * t1 + (k1 + 1) % 3]; b = x->tri[3 * t1 + (k1 + 2) % 3];
		n = rm_edge_tris(x, a, b, r);
		if (!n) continue;
		for (j = 0; j < n; j++) if (x->used[r[j] / 3]) break;
		if (j < n) continue;
		mid[0] = 0.5 * (x->p[2 * a] + x->p[2 * b]); mid[1] = 0.5 * (x->p[2 * a + 1] + x->p[2 * b + 1]);
		len = hypot(x->p[2 * b] - x->p[2 * a], x->p[2 * b + 1] - x->p[2 * a + 1]);
		if (!(len > 4.0 / 3.0 * rm_target(x, mid))) continue;
		bi = rm_boundary(x, a, b);
		if (bi != TM_NONE)
		{
			tm = 0.5 * (x->bt[2 * bi] + x->bt[2 * bi + 1]);
			tm_eval(x->m->curves[x->bc[bi]], tm, mid);
		}
		/* the new triangles must stay counter-clockwise */
		for (j = 0; j < n; j++)
		{
			unsigned int t = r[j] / 3, k = r[j] % 3, c = x->tri[3 * t + k], ea = x->tri[3 * t + (k + 1) % 3], eb = x->tri[3 * t + (k + 2) % 3];
			if (!(tri_cross(x, TM_NONE, c, ea, mid) > 0) || !(tri_cross(x, TM_NONE, eb, c, mid) > 0)) break;
		}
		if (j < n) continue;
		mv = rm_vert(x, mid[0], mid[1], bi != TM_NONE ? x->bc[bi] : TM_NONE, tm);
		if (mv == TM_NONE) return -1;
		for (j = 0; j < n; j++)
		{
			unsigned int t = r[j] / 3, k = r[j] % 3, c = x->tri[3 * t + k], ea = x->tri[3 * t + (k + 1) % 3], eb = x->tri[3 * t + (k + 2) % 3];
			x->tri[3 * t] = c; x->tri[3 * t + 1] = ea; x->tri[3 * t + 2] = mv;
			x->used[t] = 1;
			if (rm_tri(x, c, mv, eb) == TM_NONE) return -1;
		}
		if (bi != TM_NONE)
		{
			unsigned int bb = x->be[2 * bi + 1], bcv = x->bc[bi];
			double tb = x->bt[2 * bi + 1];
			x->be[2 * bi + 1] = mv; x->bt[2 * bi + 1] = tm;
			if (rm_bedge(x, mv, bb, bcv, tm, tb, x->bk[bi]) == TM_NONE) return -1;
		}
		splits++;
	}
	return x->failed ? -1 : (int)splits;
}

/* vertex stars: the triangles of v at star[first[v] .. first[v + 1]) */
static int rm_stars(rm_ctx* x, unsigned int** first, unsigned int** star)
{
	unsigned int i, k;
	unsigned int* f = (unsigned int*)calloc(x->nv + 2, sizeof(unsigned int));
	unsigned int* s = (unsigned int*)malloc(sizeof(unsigned int) * (3 * x->nt + 1));
	unsigned int* fill = (unsigned int*)malloc(sizeof(unsigned int) * (x->nv + 1));
	if (!f || !s || !fill) { free(f); free(s); free(fill); x->failed = 1; return 0; }
	for (i = 0; i < x->nt; i++)
		if (x->tri[3 * i] != TM_NONE)
			for (k = 0; k < 3; k++) f[x->tri[3 * i + k] + 1]++;
	for (i = 0; i < x->nv; i++) f[i + 1] += f[i];
	memcpy(fill, f, sizeof(unsigned int) * (x->nv + 1));
	for (i = 0; i < x->nt; i++)
		if (x->tri[3 * i] != TM_NONE)
			for (k = 0; k < 3; k++) s[fill[x->tri[3 * i + k]]++] = i;
	free(fill);
	*first = f;
	*star = s;
	return 1;
}

static int rm_has(rm_ctx const* x, unsigned int t, unsigned int v)
{
	return x->tri[3 * t] != TM_NONE && (x->tri[3 * t] == v || x->tri[3 * t + 1] == v || x->tri[3 * t + 2] == v);
}

static int rm_collapse(rm_ctx* x)
{
	unsigned int *first = NULL, *star = NULL, i, ne, done = 0;
	unsigned int* bin = NULL;     /* per vertex: the boundary edge ending there */
	unsigned int* bout = NULL;    /* and the one starting there */
	unsigned char* lock = NULL;
	if (!rm_maps(x) || !rm_stars(x, &first, &star)) return -1;
	bin = (unsigned int*)malloc(sizeof(unsigned int) * (x->nv + 1));
	bout = (unsigned int*)malloc(sizeof(unsigned int) * (x->nv + 1));
	lock = (unsigned char*)calloc(x->nv + 1, 1);
	if (!bin || !bout || !lock) { free(first); free(star); free(bin); free(bout); free(lock); return -1; }
	for (i = 0; i < x->nv; i++) bin[i] = bout[i] = TM_NONE;
	for (i = 0; i < x->nb; i++) { bout[x->be[2 * i]] = i; bin[x->be[2 * i + 1]] = i; }
	ne = x->ne;
	for (i = 0; i < ne; i++)
	{
		unsigned int t1 = x->eref[i] / 3, k1 = x->eref[i] % 3, a, b, u, w, r[2], n, bi, j, q, cn[64], ncn = 0;
		double mid[2], len;
		int ok = 1;
		if (i > 0 && x->ekey[i] == x->ekey[i - 1]) continue;
		if (x->tri[3 * t1] == TM_NONE) continue;
		a = x->tri[3 * t1 + (k1 + 1) % 3]; b = x->tri[3 * t1 + (k1 + 2) % 3];
		if (rm_key(a, b) != x->ekey[i] || lock[a] || lock[b]) continue;
		mid[0] = 0.5 * (x->p[2 * a] + x->p[2 * b]); mid[1] = 0.5 * (x->p[2 * a + 1] + x->p[2 * b + 1]);
		len = hypot(x->p[2 * b] - x->p[2 * a], x->p[2 * b + 1] - x->p[2 * a + 1]);
		if (!(len < 0.8 * rm_target(x, mid))) continue;
		bi = rm_boundary(x, a, b);
		/* which end goes: never a corner; a boundary vertex only along its
		   boundary edge (into its neighbour on the same curve) */
		{
			int ab = !x->corner[a] && (x->vc[a] == TM_NONE || bi != TM_NONE);
			int ba = !x->corner[b] && (x->vc[b] == TM_NONE || bi != TM_NONE);
			u = TM_NONE; w = TM_NONE;
			if (ab && (x->vc[a] == TM_NONE || !ba)) { u = a; w = b; }
			else if (ba) { u = b; w = a; }
		}
		if (u == TM_NONE) continue;
		n = rm_edge_tris(x, a, b, r);
		if (!n || (bi == TM_NONE && n != 2)) continue;
		/* link condition: the common neighbours of u and w are exactly the
		   vertices opposite the edge */
		for (j = first[u]; j < first[u + 1] && ok; j++)
		{
			unsigned int t = star[j], kk;
			if (x->tri[3 * t] == TM_NONE) continue;
			for (kk = 0; kk < 3; kk++)
			{
				unsigned int v = x->tri[3 * t + kk], z;
				if (v == u || v == w) continue;
				for (z = 0; z < ncn; z++) if (cn[z] == v) break;
				if (z < ncn) continue;
				for (q = first[w]; q < first[w + 1]; q++)
					if (rm_has(x, star[q], v)) break;
				if (q < first[w + 1]) { if (ncn < 64) cn[ncn++] = v; else ok = 0; }
			}
		}
		if (!ok || ncn != n) continue;
		/* the moved triangles stay counter-clockwise and short */
		for (j = first[u]; j < first[u + 1] && ok; j++)
		{
			unsigned int t = star[j], kk, o[3];
			if (x->tri[3 * t] == TM_NONE || rm_has(x, t, w)) continue;
			o[0] = x->tri[3 * t]; o[1] = x->tri[3 * t + 1]; o[2] = x->tri[3 * t + 2];
			for (kk = 0; kk < 3; kk++) if (o[kk] == u) o[kk] = w;
			if (!(tri_cross(x, o[0], o[1], o[2], NULL) > 0)) ok = 0;
			for (kk = 0; kk < 3 && ok; kk++)
				if (o[kk] != w)
				{
					double mm[2];
					mm[0] = 0.5 * (x->p[2 * o[kk]] + x->p[2 * w]); mm[1] = 0.5 * (x->p[2 * o[kk] + 1] + x->p[2 * w + 1]);
					if (hypot(x->p[2 * o[kk]] - x->p[2 * w], x->p[2 * o[kk] + 1] - x->p[2 * w + 1]) > 4.0 / 3.0 * rm_target(x, mm)) ok = 0;
				}
		}
		if (!ok) continue;
		/* lock the neighbourhood for this round, then collapse u into w */
		for (j = first[u]; j < first[u + 1]; j++)
		{
			unsigned int t = star[j], kk;
			if (x->tri[3 * t] == TM_NONE) continue;
			for (kk = 0; kk < 3; kk++) lock[x->tri[3 * t + kk]] = 1;
		}
		for (j = first[w]; j < first[w + 1]; j++)
		{
			unsigned int t = star[j], kk;
			if (x->tri[3 * t] == TM_NONE) continue;
			for (kk = 0; kk < 3; kk++) lock[x->tri[3 * t + kk]] = 1;
		}
		for (j = first[u]; j < first[u + 1]; j++)
		{
			unsigned int t = star[j], kk;
			if (x->tri[3 * t] == TM_NONE) continue;
			if (rm_has(x, t, w)) { x->tri[3 * t] = TM_NONE; continue; }
			for (kk = 0; kk < 3; kk++) if (x->tri[3 * t + kk] == u) x->tri[3 * t + kk] = w;
		}
		/* w takes u's triangles: they join w's star for later checks this round
		   only through the locks (both locked) */
		if (bi != TM_NONE)
		{
			/* the two boundary edges at u merge */
			unsigned int ein = bin[u], eout = bout[u];
			if (ein != TM_NONE && eout != TM_NONE)
			{
				if (x->be[2 * eout + 1] == w)
				{
					/* p -> u -> w : ein becomes p -> w */
					x->be[2 * ein + 1] = w; x->bt[2 * ein + 1] = x->bt[2 * eout + 1];
					bin[w] = ein;
					x->be[2 * eout] = x->be[2 * eout + 1] = TM_NONE;
				}
				else
				{
					/* w -> u -> s : eout becomes w -> s */
					x->be[2 * eout] = w; x->bt[2 * eout] = x->bt[2 * ein];
					bout[w] = eout;
					x->be[2 * ein] = x->be[2 * ein + 1] = TM_NONE;
				}
			}
		}
		done++;
	}
	/* drop dead boundary edges */
	{
		unsigned int nb = 0;
		for (i = 0; i < x->nb; i++)
			if (x->be[2 * i] != TM_NONE)
			{
				x->be[2 * nb] = x->be[2 * i]; x->be[2 * nb + 1] = x->be[2 * i + 1];
				x->bc[nb] = x->bc[i]; x->bt[2 * nb] = x->bt[2 * i]; x->bt[2 * nb + 1] = x->bt[2 * i + 1];
				x->bk[nb] = x->bk[i];
				nb++;
			}
		x->nb = nb;
	}
	free(first); free(star); free(bin); free(bout); free(lock);
	return (int)done;
}

static double rm_incircle(rm_ctx const* x, unsigned int a, unsigned int b, unsigned int c, unsigned int d)
{
	double const* D = x->p + 2 * d;
	double adx = x->p[2 * a] - D[0], ady = x->p[2 * a + 1] - D[1], bdx = x->p[2 * b] - D[0], bdy = x->p[2 * b + 1] - D[1];
	double cdx = x->p[2 * c] - D[0], cdy = x->p[2 * c + 1] - D[1];
	return (adx * adx + ady * ady) * (bdx * cdy - cdx * bdy) + (bdx * bdx + bdy * bdy) * (cdx * ady - adx * cdy) + (cdx * cdx + cdy * cdy) * (adx * bdy - bdx * ady);
}

static int rm_flip(rm_ctx* x)
{
	unsigned int i, ne, flips = 0;
	if (!rm_maps(x)) return -1;
	ne = x->ne;
	memset(x->used, 0, x->capt);
	for (i = 0; i < ne; i++)
	{
		unsigned int r[2], t, k, u, j, a, b, c, d;
		double s, e;
		if (i > 0 && x->ekey[i] == x->ekey[i - 1]) continue;
		t = x->eref[i] / 3; k = x->eref[i] % 3;
		if (x->tri[3 * t] == TM_NONE) continue;
		a = x->tri[3 * t + (k + 1) % 3]; b = x->tri[3 * t + (k + 2) % 3];
		if (rm_key(a, b) != x->ekey[i] || rm_boundary(x, a, b) != TM_NONE) continue;
		if (rm_edge_tris(x, a, b, r) != 2) continue;
		t = r[0] / 3; k = r[0] % 3; u = r[1] / 3; j = r[1] % 3;
		if (x->used[t] || x->used[u]) continue;
		c = x->tri[3 * t + k]; a = x->tri[3 * t + (k + 1) % 3]; b = x->tri[3 * t + (k + 2) % 3];
		d = x->tri[3 * u + j];
		e = hypot(x->p[2 * a] - x->p[2 * b], x->p[2 * a + 1] - x->p[2 * b + 1]);
		s = e * e;
		if (rm_incircle(x, c, a, b, d) <= 1e-9 * s * s) continue;
		if (!(tri_cross(x, c, a, d, NULL) > 0) || !(tri_cross(x, d, b, c, NULL) > 0)) continue;
		x->tri[3 * t] = c; x->tri[3 * t + 1] = a; x->tri[3 * t + 2] = d;
		x->tri[3 * u] = d; x->tri[3 * u + 1] = b; x->tri[3 * u + 2] = c;
		x->used[t] = x->used[u] = 1;
		flips++;
	}
	return (int)flips;
}

static int rm_smooth(rm_ctx* x)
{
	unsigned int *first = NULL, *star = NULL, v, moved = 0;
	unsigned int* bin = (unsigned int*)malloc(sizeof(unsigned int) * (x->nv + 1));
	unsigned int* bout = (unsigned int*)malloc(sizeof(unsigned int) * (x->nv + 1));
	if (!bin || !bout || !rm_stars(x, &first, &star)) { free(bin); free(bout); free(first); free(star); return -1; }
	for (v = 0; v < x->nv; v++) bin[v] = bout[v] = TM_NONE;
	for (v = 0; v < x->nb; v++) { bout[x->be[2 * v]] = v; bin[x->be[2 * v + 1]] = v; }
	for (v = 0; v < x->nv; v++)
	{
		double np[2], nt = 0, acc[2] = { 0, 0 }, cnt = 0;
		unsigned int j;
		int ok = 1;
		if (x->corner[v] || first[v] == first[v + 1]) continue;
		if (x->vc[v] != TM_NONE)
		{
			/* along the curve: halfway to the parameter midpoint of its neighbours */
			unsigned int ein = bin[v], eout = bout[v];
			if (ein == TM_NONE || eout == TM_NONE || x->bc[ein] != x->bc[eout]) continue;
			nt = x->vt[v] + 0.5 * (0.5 * (x->bt[2 * ein] + x->bt[2 * eout + 1]) - x->vt[v]);
			tm_eval(x->m->curves[x->bc[ein]], nt, np);
		}
		else
		{
			/* halfway to the centroid of its neighbours */
			for (j = first[v]; j < first[v + 1]; j++)
			{
				unsigned int t = star[j], k;
				for (k = 0; k < 3; k++)
					if (x->tri[3 * t + k] != v) { acc[0] += x->p[2 * x->tri[3 * t + k]]; acc[1] += x->p[2 * x->tri[3 * t + k] + 1]; cnt++; }
			}
			if (!(cnt > 0)) continue;
			np[0] = x->p[2 * v] + 0.5 * (acc[0] / cnt - x->p[2 * v]);
			np[1] = x->p[2 * v + 1] + 0.5 * (acc[1] / cnt - x->p[2 * v + 1]);
		}
		for (j = first[v]; j < first[v + 1] && ok; j++)
		{
			unsigned int t = star[j], k;
			for (k = 0; k < 3; k++)
				if (x->tri[3 * t + k] == v)
					ok = tri_cross(x, TM_NONE, x->tri[3 * t + (k + 1) % 3], x->tri[3 * t + (k + 2) % 3], np) > 0;
		}
		if (!ok) continue;
		x->p[2 * v] = np[0]; x->p[2 * v + 1] = np[1];
		if (x->vc[v] != TM_NONE)
		{
			x->vt[v] = nt;
			x->bt[2 * bin[v] + 1] = nt;
			x->bt[2 * bout[v]] = nt;
		}
		moved++;
	}
	free(bin); free(bout); free(first); free(star);
	return (int)moved;
}

/* the working mesh back into the output arrays (unused vertices dropped):
   its live triangles, and quads (4 working indices each) if given */
static int rm_store(rm_ctx* x, unsigned int const* quads, unsigned int nq)
{
	qaws_mesh_2d* m = x->m;
	unsigned int* map = (unsigned int*)malloc(sizeof(unsigned int) * (x->nv + 1));
	unsigned int i, k, kind, nt = 0, nv = 0;
	if (!map) return 0;
	for (i = 0; i < x->nv; i++) map[i] = TM_NONE;
	for (i = 0; i < x->nt; i++)
		if (x->tri[3 * i] != TM_NONE)
		{
			nt++;
			for (k = 0; k < 3; k++) if (map[x->tri[3 * i + k]] == TM_NONE) map[x->tri[3 * i + k]] = nv++;
		}
	for (i = 0; i < nq; i++)
		for (k = 0; k < 4; k++) if (map[quads[4 * i + k]] == TM_NONE) map[quads[4 * i + k]] = nv++;
	free(m->ov); free(m->otri); free(m->ovc); free(m->ovt); free(m->obe); free(m->obc); free(m->obt); free(m->obk); free(m->oquad);
	m->ov = (double*)malloc(sizeof(double) * 2 * (nv + 1));
	m->ovc = (unsigned int*)malloc(sizeof(unsigned int) * (nv + 1));
	m->ovt = (double*)malloc(sizeof(double) * (nv + 1));
	m->otri = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (nt + 1));
	m->oquad = (unsigned int*)malloc(sizeof(unsigned int) * 4 * (nq + 1));
	m->obe = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (x->nb + 1));
	m->obc = (unsigned int*)malloc(sizeof(unsigned int) * (x->nb + 1));
	m->obt = (double*)malloc(sizeof(double) * 2 * (x->nb + 1));
	m->obk = (unsigned char*)malloc(x->nb + 1);
	m->onv = m->ont = m->onbe = m->once = m->onq = 0;
	if (!m->ov || !m->ovc || !m->ovt || !m->otri || !m->oquad || !m->obe || !m->obc || !m->obt || !m->obk) { free(map); return 0; }
	for (i = 0; i < x->nv; i++)
		if (map[i] != TM_NONE)
		{
			m->ov[2 * map[i]] = x->p[2 * i]; m->ov[2 * map[i] + 1] = x->p[2 * i + 1];
			m->ovc[map[i]] = x->vc[i]; m->ovt[map[i]] = x->vt[i];
		}
	m->onv = nv;
	for (i = 0; i < x->nt; i++)
		if (x->tri[3 * i] != TM_NONE)
		{
			for (k = 0; k < 3; k++) m->otri[3 * m->ont + k] = map[x->tri[3 * i + k]];
			m->ont++;
		}
	for (i = 0; i < nq; i++)
		for (k = 0; k < 4; k++) m->oquad[4 * i + k] = map[quads[4 * i + k]];
	m->onq = nq;
	for (kind = 0; kind < 2; kind++)
		for (i = 0; i < x->nb; i++)
		{
			unsigned int e = m->onbe + m->once;
			if (x->bk[i] != kind) continue;
			m->obe[2 * e] = map[x->be[2 * i]]; m->obe[2 * e + 1] = map[x->be[2 * i + 1]];
			m->obc[e] = x->bc[i]; m->obt[2 * e] = x->bt[2 * i]; m->obt[2 * e + 1] = x->bt[2 * i + 1];
			m->obk[e] = (unsigned char)kind;
			if (kind) m->once++; else m->onbe++;
		}
	free(map);
	return 1;
}

static void rm_free(rm_ctx* x)
{
	free(x->p); free(x->vc); free(x->vt); free(x->corner); free(x->tri); free(x->be); free(x->bc); free(x->bt); free(x->bk);
	free(x->ekey); free(x->eref); free(x->bkey); free(x->bref); free(x->used);
}

/* the output mesh into a working one; corners: boundary vertices with other
   than two boundary edges, or where two curves meet */
static int rm_load(rm_ctx* x, qaws_mesh_2d* m)
{
	unsigned int i;
	unsigned int *bin, *bout, *deg;
	for (i = 0; i < m->onv; i++)
		if (rm_vert(x, m->ov[2 * i], m->ov[2 * i + 1], m->ovc[i], m->ovt[i]) == TM_NONE) return 0;
	for (i = 0; i < m->ont; i++)
		if (rm_tri(x, m->otri[3 * i], m->otri[3 * i + 1], m->otri[3 * i + 2]) == TM_NONE) return 0;
	for (i = 0; i < m->onbe + m->once; i++)
		if (rm_bedge(x, m->obe[2 * i], m->obe[2 * i + 1], m->obc[i], m->obt[2 * i], m->obt[2 * i + 1], m->obk[i]) == TM_NONE) return 0;
	bin = (unsigned int*)malloc(sizeof(unsigned int) * (x->nv + 1));
	bout = (unsigned int*)malloc(sizeof(unsigned int) * (x->nv + 1));
	deg = (unsigned int*)calloc(x->nv + 1, sizeof(unsigned int));
	if (!bin || !bout || !deg) { free(bin); free(bout); free(deg); return 0; }
	for (i = 0; i < x->nv; i++) bin[i] = bout[i] = TM_NONE;
	for (i = 0; i < x->nb; i++)
	{
		bout[x->be[2 * i]] = i; bin[x->be[2 * i + 1]] = i;
		deg[x->be[2 * i]]++; deg[x->be[2 * i + 1]]++;
	}
	for (i = 0; i < x->nv; i++)
	{
		if (x->vc[i] == TM_NONE) continue;
		/* the seam of a closed curve (its parameter jumps from end to start)
		   counts as a corner too */
		if (deg[i] != 2 || bin[i] == TM_NONE || bout[i] == TM_NONE || x->bc[bin[i]] != x->bc[bout[i]]
			|| x->bt[2 * bin[i] + 1] != x->bt[2 * bout[i]])
			x->corner[i] = 1;
		else
		{
			/* a kink of the curve (a polyline's corner, a C0 knot): where its
			   tangent jumps at a span boundary */
			qaws_curve const* c = x->m->curves[x->vc[i]];
			unsigned int sp;
			double w = c->parameter_range.max_value - c->parameter_range.min_value;
			for (sp = 1; sp < c->span_count; sp++)
				if (fabs(x->vt[i] - c->span_boundaries[sp]) <= 1e-12 * (w + 1))
				{
					/* chords just before and just after (an evaluator may pick either
					   span at the knot itself) */
					double h = 1e-7 * (c->span_boundaries[sp + 1] - c->span_boundaries[sp - 1]), p0[2], p1[2], p2[2], l0, l1;
					tm_eval(c, x->vt[i] - h, p0);
					tm_eval(c, x->vt[i], p1);
					tm_eval(c, x->vt[i] + h, p2);
					l0 = hypot(p1[0] - p0[0], p1[1] - p0[1]); l1 = hypot(p2[0] - p1[0], p2[1] - p1[1]);
					if (!(l0 > 0 && l1 > 0) || ((p1[0] - p0[0]) * (p2[0] - p1[0]) + (p1[1] - p0[1]) * (p2[1] - p1[1])) < 0.99985 * l0 * l1)
						x->corner[i] = 1;
					break;
				}
		}
	}
	free(bin); free(bout); free(deg);
	return 1;
}

qaws_status qaws_mesh_2d_remesh(qaws_mesh_2d* mesh, qaws_scalar target_length, qaws_size_fn size_fn, void* size_user,
	unsigned int iterations)
{
	rm_ctx x;
	unsigned int it;
	int ok = 1;
	if (!mesh || (!(target_length > 0) && !size_fn)) return QAWS_STATUS_INVALID_ARGUMENT;
	if (mesh->onq) return QAWS_STATUS_UNSUPPORTED_OPERATION;
	memset(&x, 0, sizeof(x));
	x.m = mesh;
	x.len = (double)target_length;
	x.size_fn = size_fn;
	x.size_user = size_user;
	if (!rm_load(&x, mesh)) { rm_free(&x); return QAWS_STATUS_ALLOCATION_FAILURE; }
	if (!(x.len > 0))
	{
		/* no target length beside the field: the mean edge to start with */
		unsigned int i, k;
		double s = 0;
		for (i = 0; i < x.nt; i++)
			for (k = 0; k < 3; k++)
				s += hypot(x.p[2 * x.tri[3 * i + k]] - x.p[2 * x.tri[3 * i + (k + 1) % 3]], x.p[2 * x.tri[3 * i + k] + 1] - x.p[2 * x.tri[3 * i + (k + 1) % 3] + 1]);
		x.len = x.nt ? s / (3 * x.nt) : 1;
	}
	if (!iterations) iterations = 5;
	for (it = 0; it < iterations && ok; it++)
	{
		unsigned int r;
		int s;
		for (r = 0; r < 12; r++) { s = rm_split(&x); if (s < 0) { ok = 0; break; } if (!s) break; }
		for (r = 0; r < 4 && ok; r++) { s = rm_collapse(&x); if (s < 0) { ok = 0; break; } if (!s) break; }
		for (r = 0; r < 8 && ok; r++) { s = rm_flip(&x); if (s < 0) { ok = 0; break; } if (!s) break; }
		for (r = 0; r < 3 && ok; r++)
		{
			if (rm_smooth(&x) < 0) { ok = 0; break; }
			if (rm_flip(&x) < 0) ok = 0;
		}
	}
	if (ok) ok = rm_store(&x, NULL, 0);
	rm_free(&x);
	return ok ? QAWS_STATUS_OK : QAWS_STATUS_ALLOCATION_FAILURE;
}
/* ======================================================================== */
/*  Quadrangulation                                                         */
/* ======================================================================== */

/* the corners of face f (4 slots, a triangle's last TM_NONE) */
static unsigned int qd_size(unsigned int const* f)
{
	return f[3] == TM_NONE ? 3u : 4u;
}

/* the worst corner of a polygon off a right angle (quadrilateral) or off 60
   degrees (triangle), in degrees; 1e9 when not strictly convex */
static double qd_shape(rm_ctx const* x, unsigned int const* f, unsigned int n, double const* moved, unsigned int mv)
{
	unsigned int i;
	double worst = 0, ideal = n == 4 ? 90 : 60;
	for (i = 0; i < n; i++)
	{
		unsigned int a = f[(i + n - 1) % n], b = f[i], c = f[(i + 1) % n];
		double const* A = a == mv ? moved : x->p + 2 * a;
		double const* B = b == mv ? moved : x->p + 2 * b;
		double const* C = c == mv ? moved : x->p + 2 * c;
		double ux = A[0] - B[0], uy = A[1] - B[1], wx = C[0] - B[0], wy = C[1] - B[1];
		double cr = wx * uy - wy * ux, ang;
		if (!(cr > 0)) return 1e9;
		ang = atan2(cr, ux * wx + uy * wy) * 180 / TM_PI;
		if (fabs(ang - ideal) > worst) worst = fabs(ang - ideal);
	}
	return worst;
}

typedef struct qd_faces
{
	unsigned int* f;            /* 4 per face */
	unsigned int n, cap;
} qd_faces;

static int qd_add(qd_faces* q, unsigned int a, unsigned int b, unsigned int c, unsigned int d)
{
	if (q->n == q->cap)
	{
		unsigned int nc = q->cap ? 2 * q->cap : 256;
		unsigned int* g = (unsigned int*)realloc(q->f, sizeof(unsigned int) * 4 * nc);
		if (!g) return 0;
		q->f = g;
		q->cap = nc;
	}
	q->f[4 * q->n] = a; q->f[4 * q->n + 1] = b; q->f[4 * q->n + 2] = c; q->f[4 * q->n + 3] = d;
	q->n++;
	return 1;
}

/* every face split through its edge midpoints and its centroid into
   quadrilaterals; feature edges split on their curves */
static int qd_split_all(rm_ctx* x, qd_faces* faces)
{
	uint64_t* key;
	unsigned int* ref;
	unsigned int* mid;
	unsigned int i, k, n = 0, nf = faces->n;
	qd_faces out;
	memset(&out, 0, sizeof(out));
	key = (uint64_t*)malloc(sizeof(uint64_t) * (4 * nf + 1));
	ref = (unsigned int*)malloc(sizeof(unsigned int) * (4 * nf + 1));
	mid = (unsigned int*)malloc(sizeof(unsigned int) * (4 * nf + 1));
	if (!key || !ref || !mid) { free(key); free(ref); free(mid); return 0; }
	for (i = 0; i < nf; i++)
	{
		unsigned int const* f = faces->f + 4 * i;
		unsigned int s = qd_size(f);
		for (k = 0; k < s; k++)
		{
			key[n] = rm_key(f[k], f[(k + 1) % s]);
			ref[n] = 4 * i + k;
			n++;
		}
	}
	rm_sort(key, ref, n);
	/* one midpoint per edge */
	for (i = 0; i < n; i++)
	{
		unsigned int f, k2, a, b, bi, v;
		double p[2], tm = 0;
		if (i > 0 && key[i] == key[i - 1]) { mid[i] = mid[i - 1]; continue; }
		f = ref[i] / 4; k2 = ref[i] % 4;
		a = faces->f[4 * f + k2]; b = faces->f[4 * f + (k2 + 1) % qd_size(faces->f + 4 * f)];
		bi = rm_boundary(x, a, b);
		p[0] = 0.5 * (x->p[2 * a] + x->p[2 * b]); p[1] = 0.5 * (x->p[2 * a + 1] + x->p[2 * b + 1]);
		if (bi != TM_NONE)
		{
			tm = 0.5 * (x->bt[2 * bi] + x->bt[2 * bi + 1]);
			tm_eval(x->m->curves[x->bc[bi]], tm, p);
		}
		v = rm_vert(x, p[0], p[1], bi != TM_NONE ? x->bc[bi] : TM_NONE, tm);
		if (v == TM_NONE) { free(key); free(ref); free(mid); return 0; }
		if (bi != TM_NONE)
		{
			unsigned int bb = x->be[2 * bi + 1], bcv = x->bc[bi];
			double tb = x->bt[2 * bi + 1];
			unsigned char kind = x->bk[bi];
			x->be[2 * bi + 1] = v; x->bt[2 * bi + 1] = tm;
			if (rm_bedge(x, v, bb, bcv, tm, tb, kind) == TM_NONE) { free(key); free(ref); free(mid); return 0; }
		}
		mid[i] = v;
	}
	/* the midpoint of each face edge, back by the face's order */
	{
		unsigned int* fm = (unsigned int*)malloc(sizeof(unsigned int) * (4 * nf + 1));
		if (!fm) { free(key); free(ref); free(mid); return 0; }
		for (i = 0; i < n; i++) fm[ref[i]] = mid[i];
		for (i = 0; i < nf && !x->failed; i++)
		{
			unsigned int const* f = faces->f + 4 * i;
			unsigned int s = qd_size(f), c;
			double cx = 0, cy = 0;
			for (k = 0; k < s; k++) { cx += x->p[2 * f[k]]; cy += x->p[2 * f[k] + 1]; }
			c = rm_vert(x, cx / s, cy / s, TM_NONE, 0);
			if (c == TM_NONE) break;
			for (k = 0; k < s; k++)
				if (!qd_add(&out, f[k], fm[4 * i + k], c, fm[4 * i + (k + s - 1) % s])) { x->failed = 1; break; }
		}
		free(fm);
	}
	free(key); free(ref); free(mid);
	if (x->failed) { free(out.f); return 0; }
	free(faces->f);
	*faces = out;
	return 1;
}

/* smoothing round faces: free vertices halfway to their edge neighbours'
   centroid, feature vertices along their curves; a move that leaves a face
   not convex is not made */
static int qd_smooth(rm_ctx* x, qd_faces const* faces)
{
	unsigned int *first, *star, *fill, *bin, *bout, i, k, v;
	first = (unsigned int*)calloc(x->nv + 2, sizeof(unsigned int));
	star = (unsigned int*)malloc(sizeof(unsigned int) * (4 * faces->n + 1));
	fill = (unsigned int*)malloc(sizeof(unsigned int) * (x->nv + 1));
	bin = (unsigned int*)malloc(sizeof(unsigned int) * (x->nv + 1));
	bout = (unsigned int*)malloc(sizeof(unsigned int) * (x->nv + 1));
	if (!first || !star || !fill || !bin || !bout) { free(first); free(star); free(fill); free(bin); free(bout); return 0; }
	for (i = 0; i < faces->n; i++)
		for (k = 0; k < qd_size(faces->f + 4 * i); k++) first[faces->f[4 * i + k] + 1]++;
	for (i = 0; i < x->nv; i++) first[i + 1] += first[i];
	memcpy(fill, first, sizeof(unsigned int) * (x->nv + 1));
	for (i = 0; i < faces->n; i++)
		for (k = 0; k < qd_size(faces->f + 4 * i); k++) star[fill[faces->f[4 * i + k]]++] = i;
	for (v = 0; v < x->nv; v++) bin[v] = bout[v] = TM_NONE;
	for (i = 0; i < x->nb; i++) { bout[x->be[2 * i]] = i; bin[x->be[2 * i + 1]] = i; }
	for (v = 0; v < x->nv; v++)
	{
		double np[2], nt = 0, acc[2] = { 0, 0 }, cnt = 0, before = 0, after = 0;
		unsigned int j;
		if (x->corner[v] || first[v] == first[v + 1]) continue;
		if (x->vc[v] != TM_NONE)
		{
			unsigned int ein = bin[v], eout = bout[v];
			if (ein == TM_NONE || eout == TM_NONE || x->bc[ein] != x->bc[eout]) continue;
			nt = x->vt[v] + 0.5 * (0.5 * (x->bt[2 * ein] + x->bt[2 * eout + 1]) - x->vt[v]);
			tm_eval(x->m->curves[x->bc[ein]], nt, np);
		}
		else
		{
			/* edge neighbours: the corners next to v in its faces */
			for (j = first[v]; j < first[v + 1]; j++)
			{
				unsigned int const* f = faces->f + 4 * star[j];
				unsigned int s = qd_size(f);
				for (k = 0; k < s; k++)
					if (f[k] == v)
					{
						unsigned int a = f[(k + 1) % s], b = f[(k + s - 1) % s];
						acc[0] += x->p[2 * a] + x->p[2 * b];
						acc[1] += x->p[2 * a + 1] + x->p[2 * b + 1];
						cnt += 2;
					}
			}
			if (!(cnt > 0)) continue;
			np[0] = x->p[2 * v] + 0.5 * (acc[0] / cnt - x->p[2 * v]);
			np[1] = x->p[2 * v + 1] + 0.5 * (acc[1] / cnt - x->p[2 * v + 1]);
		}
		/* kept convex, and no worse at its worst corner */
		for (j = first[v]; j < first[v + 1]; j++)
		{
			unsigned int const* f = faces->f + 4 * star[j];
			double b0 = qd_shape(x, f, qd_size(f), NULL, TM_NONE), a0 = qd_shape(x, f, qd_size(f), np, v);
			if (b0 > before) before = b0;
			if (a0 > after) after = a0;
		}
		if (after >= 1e9 || after > before + 1e-9) continue;
		x->p[2 * v] = np[0]; x->p[2 * v + 1] = np[1];
		if (x->vc[v] != TM_NONE)
		{
			x->vt[v] = nt;
			x->bt[2 * bin[v] + 1] = nt;
			x->bt[2 * bout[v]] = nt;
		}
	}
	free(first); free(star); free(fill); free(bin); free(bout);
	return 1;
}

/* pairing state: each triangle's partner and the candidate joining them */
typedef struct qd_match
{
	unsigned int* mate;
	unsigned int* via;
	unsigned int* first;        /* candidates of triangle t: list[first[t] .. first[t + 1]) */
	unsigned int* list;
	unsigned int const* cq;     /* 6 per candidate: t, u, quad */
	unsigned int* seen;
	unsigned int stamp;
} qd_match;

static void qd_pair(qd_match* q, unsigned int a, unsigned int b, unsigned int c)
{
	q->mate[a] = b; q->mate[b] = a;
	q->via[a] = q->via[b] = c;
}

/* an alternating path from unmatched a, at most depth re-pairings long
   (Blossom-quad's augmenting paths, without blossoms) */
static int qd_augment(qd_match* q, unsigned int a, unsigned int depth)
{
	unsigned int j;
	for (j = q->first[a]; j < q->first[a + 1]; j++)
	{
		unsigned int c = q->list[j], b = q->cq[6 * c] == a ? q->cq[6 * c + 1] : q->cq[6 * c];
		if (q->seen[b] == q->stamp) continue;
		q->seen[b] = q->stamp;
		if (q->mate[b] == TM_NONE) { qd_pair(q, a, b, c); return 1; }
		if (!depth) continue;
		{
			unsigned int m = q->mate[b];
			if (q->seen[m] == q->stamp) continue;
			q->seen[m] = q->stamp;
			q->mate[b] = TM_NONE;     /* b offered to a while m looks further */
			if (qd_augment(q, m, depth - 1)) { qd_pair(q, a, b, c); return 1; }
			q->mate[b] = m;
		}
	}
	return 0;
}

qaws_status qaws_mesh_2d_quadrangulate(qaws_mesh_2d* mesh, int all_quads, unsigned int smoothing_iterations)
{
	rm_ctx x;
	qd_faces faces;
	uint64_t* key = NULL;
	unsigned int* ref = NULL;
	unsigned char* matched = NULL;
	unsigned int i, nc = 0, it, *cq = NULL;
	int ok = 1;
	if (!mesh) return QAWS_STATUS_INVALID_ARGUMENT;
	if (mesh->onq) return QAWS_STATUS_UNSUPPORTED_OPERATION;
	memset(&x, 0, sizeof(x));
	memset(&faces, 0, sizeof(faces));
	x.m = mesh;
	if (!rm_load(&x, mesh) || !rm_maps(&x)) { rm_free(&x); return QAWS_STATUS_ALLOCATION_FAILURE; }
	/* candidate pairs: two triangles over a free edge, convex together */
	key = (uint64_t*)malloc(sizeof(uint64_t) * (x.ne + 1));
	ref = (unsigned int*)malloc(sizeof(unsigned int) * (x.ne + 1));
	cq = (unsigned int*)malloc(sizeof(unsigned int) * 6 * (x.ne + 1));
	matched = (unsigned char*)calloc(x.nt + 1, 1);
	if (!key || !ref || !cq || !matched) ok = 0;
	for (i = 0; ok && i < x.ne; i++)
	{
		unsigned int r[2], t, k, u, j, a, b, c, d, q[4];
		double w;
		if (i > 0 && x.ekey[i] == x.ekey[i - 1]) continue;
		t = x.eref[i] / 3; k = x.eref[i] % 3;
		a = x.tri[3 * t + (k + 1) % 3]; b = x.tri[3 * t + (k + 2) % 3];
		if (rm_boundary(&x, a, b) != TM_NONE || rm_edge_tris(&x, a, b, r) != 2) continue;
		t = r[0] / 3; k = r[0] % 3; u = r[1] / 3; j = r[1] % 3;
		c = x.tri[3 * t + k]; a = x.tri[3 * t + (k + 1) % 3]; b = x.tri[3 * t + (k + 2) % 3]; d = x.tri[3 * u + j];
		q[0] = c; q[1] = a; q[2] = d; q[3] = b;
		w = qd_shape(&x, q, 4, NULL, TM_NONE);
		if (w >= 1e9) continue;
		/* the better shaped first: keys by the worst corner off 90 degrees */
		key[nc] = (uint64_t)(w * 1e6);
		ref[nc] = nc;
		cq[6 * nc] = t; cq[6 * nc + 1] = u;
		cq[6 * nc + 2] = c; cq[6 * nc + 3] = a; cq[6 * nc + 4] = d; cq[6 * nc + 5] = b;
		nc++;
	}
	if (ok) rm_sort(key, ref, nc);
	/* greedy, the best shaped pairs first; then augmenting paths for the
	   triangles left */
	if (ok)
	{
		qd_match q;
		unsigned int t, pass;
		memset(&q, 0, sizeof(q));
		q.mate = (unsigned int*)malloc(sizeof(unsigned int) * (x.nt + 1));
		q.via = (unsigned int*)malloc(sizeof(unsigned int) * (x.nt + 1));
		q.first = (unsigned int*)calloc(x.nt + 2, sizeof(unsigned int));
		q.list = (unsigned int*)malloc(sizeof(unsigned int) * (2 * nc + 1));
		q.seen = (unsigned int*)calloc(x.nt + 1, sizeof(unsigned int));
		q.cq = cq;
		if (!q.mate || !q.via || !q.first || !q.list || !q.seen) ok = 0;
		if (ok)
		{
			unsigned int* fill = (unsigned int*)malloc(sizeof(unsigned int) * (x.nt + 1));
			if (!fill) ok = 0;
			else
			{
				for (t = 0; t < x.nt; t++) q.mate[t] = TM_NONE;
				/* candidates per triangle, best first */
				for (i = 0; i < nc; i++) { q.first[cq[6 * ref[i]] + 1]++; q.first[cq[6 * ref[i] + 1] + 1]++; }
				for (t = 0; t < x.nt; t++) q.first[t + 1] += q.first[t];
				memcpy(fill, q.first, sizeof(unsigned int) * (x.nt + 1));
				for (i = 0; i < nc; i++)
				{
					unsigned int c = ref[i];
					q.list[fill[cq[6 * c]]++] = c;
					q.list[fill[cq[6 * c + 1]]++] = c;
				}
				free(fill);
				for (i = 0; i < nc; i++)
				{
					unsigned int c = ref[i];
					if (q.mate[cq[6 * c]] == TM_NONE && q.mate[cq[6 * c + 1]] == TM_NONE)
						qd_pair(&q, cq[6 * c], cq[6 * c + 1], c);
				}
				for (pass = 0; pass < 6; pass++)
				{
					unsigned int gained = 0;
					for (t = 0; t < x.nt; t++)
					{
						if (x.tri[3 * t] == TM_NONE || q.mate[t] != TM_NONE || q.first[t] == q.first[t + 1]) continue;
						q.stamp++;
						q.seen[t] = q.stamp;
						gained += (unsigned int)qd_augment(&q, t, 8);
					}
					if (!gained) break;
				}
				for (t = 0; ok && t < x.nt; t++)
				{
					if (x.tri[3 * t] == TM_NONE) continue;
					if (q.mate[t] == TM_NONE)
					{
						if (!qd_add(&faces, x.tri[3 * t], x.tri[3 * t + 1], x.tri[3 * t + 2], TM_NONE)) ok = 0;
					}
					else if (t < q.mate[t])
					{
						unsigned int const* c = cq + 6 * q.via[t];
						if (!qd_add(&faces, c[2], c[3], c[4], c[5])) ok = 0;
					}
				}
			}
		}
		free(q.mate); free(q.via); free(q.first); free(q.list); free(q.seen);
	}
	free(key); free(ref); free(cq); free(matched);
	if (ok && all_quads && !qd_split_all(&x, &faces)) ok = 0;
	if (ok) ok = rm_maps(&x);
	for (it = 0; ok && it < smoothing_iterations; it++)
		if (!qd_smooth(&x, &faces)) ok = 0;
	if (ok)
	{
		/* triangles back into x.tri, quads out */
		unsigned int* quads = (unsigned int*)malloc(sizeof(unsigned int) * 4 * (faces.n + 1));
		unsigned int nq = 0;
		if (!quads) ok = 0;
		x.nt = 0;
		for (i = 0; ok && i < faces.n; i++)
		{
			unsigned int const* f = faces.f + 4 * i;
			if (f[3] == TM_NONE)
			{
				if (rm_tri(&x, f[0], f[1], f[2]) == TM_NONE) ok = 0;
			}
			else
			{
				memcpy(quads + 4 * nq, f, sizeof(unsigned int) * 4);
				nq++;
			}
		}
		if (ok) ok = rm_store(&x, quads, nq);
		free(quads);
	}
	free(faces.f);
	rm_free(&x);
	return ok ? QAWS_STATUS_OK : QAWS_STATUS_ALLOCATION_FAILURE;
}

unsigned int qaws_mesh_2d_get_quad_count(qaws_mesh_2d const* m) { return m ? m->onq : 0; }

void qaws_mesh_2d_get_quad(unsigned int out[4], qaws_mesh_2d const* m, unsigned int q)
{
	out[0] = m->oquad[4 * q];
	out[1] = m->oquad[4 * q + 1];
	out[2] = m->oquad[4 * q + 2];
	out[3] = m->oquad[4 * q + 3];
}

unsigned int qaws_mesh_2d_get_constraint_edge_count(qaws_mesh_2d const* m) { return m ? m->once : 0; }

void qaws_mesh_2d_get_constraint_edge(unsigned int out[2], qaws_mesh_2d const* m, unsigned int e)
{
	out[0] = m->obe[2 * (m->onbe + e)];
	out[1] = m->obe[2 * (m->onbe + e) + 1];
}
