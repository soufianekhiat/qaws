/*
 * Quad finishing after INSTANT's extraction (Tarini, Pietroni, Cignoni,
 * Panozzo, Puppo 2010, practical quad mesh simplification; Bommes et al.
 * 2011). Each vertex has an ideal valence: 4 inside; on a crease, the border
 * or a corner, the quads its corner angles hold (plus one on the border), so
 * a box's corners want 3 and a can's rim 4. E, the sum of the squared
 * deviations, is lowered by local operations (edge rotations, diagonal
 * collapses, doublet removals; a rotation and a second operation beside it
 * where none alone helps), each tried first on a copy of its patch, relaxed
 * there, and kept only when its quads stay shapely and face the reference.
 * Corners keep their valence, line vertices never go below theirs, feature
 * edges are never touched. A quad turned inside out (the extraction may
 * leave a few) has its diagonals collapsed whatever E says, kept only when
 * its patch's smallest scaled Jacobian rises.
 *
 * Then relaxation on the reference, vertex by vertex while their
 * neighbourhoods move: a free vertex goes to the better of Zhou and Shimada's
 * angle target (2000: each neighbour's edge turned to bisect the
 * neighbour's two angles beside it, the mean of the ends) and a Newton step
 * down a shape energy (the angles' deviations from the regular polygons',
 * the edges' from the target length, weighted lightly), walked; a line
 * vertex to the middle of its line neighbours along the creases; a move kept
 * only when the smallest scaled Jacobian round the vertex does not drop (and
 * no more of its quads are at or under 0), no edge there gets shorter than
 * half the target length (or than it was) and no quad's edge ratio falls
 * under 0.3 (or than it was), and for free
 * vertices the energy drops. Last, the quads under 0.5 worst first: their
 * free vertices searched in eight directions under the same guards.
 */

#include "qaws_quad_improve.h"
#include "../internal/qaws_internal_memory.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define QI_FS 8                 /* corners of a face at most */
#define QI_FAN 16               /* faces round a vertex at most */
#define QI_LV 96                /* a patch's vertices */
#define QI_LF 64                /* and faces */
#define QI_PI 3.14159265358979323846
#define QI_LW 0.1               /* the lengths' weight against the angles' in the shape energy */
#define QI_RATIO 0.3            /* no move takes a quad's shortest to longest edge under this (under it already: lower) */

static double qi_dot(double const* a, double const* b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void qi_cross(double const* a, double const* b, double* o)
{
	o[0] = a[1] * b[2] - a[2] * b[1];
	o[1] = a[2] * b[0] - a[0] * b[2];
	o[2] = a[0] * b[1] - a[1] * b[0];
}

static void qi_sub(double const* a, double const* b, double* o)
{
	o[0] = a[0] - b[0];
	o[1] = a[1] - b[1];
	o[2] = a[2] - b[2];
}

/* a face's smallest scaled Jacobian about the reference normal nrm ((e1 x e2) . n / |e1| |e2| at each corner; 1 for a
   face that is not a quad) */
static double qi_sj(double const* const* p, unsigned int n, double const* nrm)
{
	double worst = 1;
	unsigned int i;
	if (n != 4)
		return 1;
	for (i = 0; i < 4; i++)
	{
		double e1[3], e2[3], x[3], l, s;
		qi_sub(p[(i + 1) & 3], p[i], e1);
		qi_sub(p[(i + 3) & 3], p[i], e2);
		qi_cross(e1, e2, x);
		l = sqrt(qi_dot(e1, e1) * qi_dot(e2, e2));
		s = l > 0 ? qi_dot(x, nrm) / l : -1;
		if (s < worst)
			worst = s;
	}
	return worst;
}

/* a face's Newell normal along nrm */
static int qi_up(double const* const* p, unsigned int n, double const* nrm)
{
	double s[3] = { 0, 0, 0 };
	unsigned int i;
	for (i = 0; i < n; i++)
	{
		double const *a = p[i], *b = p[(i + 1) % n];
		s[0] += (a[1] - b[1]) * (a[2] + b[2]);
		s[1] += (a[2] - b[2]) * (a[0] + b[0]);
		s[2] += (a[0] - b[0]) * (a[1] + b[1]);
	}
	return qi_dot(s, nrm) > 0;
}

/* the counter-clockwise angle from u to v about n, in [0, 2 pi) */
static double qi_turn(double const* u, double const* v, double const* n)
{
	double x[3], a;
	qi_cross(u, v, x);
	a = atan2(qi_dot(x, n), qi_dot(u, v));
	return a < 0 ? a + 2 * QI_PI : a;
}

typedef struct qi_cand
{
	int de;
	unsigned int type, a, b;        /* 0: the edge a -> b rotated by dir; 1: a quad's diagonal a, b collapsed into a; 2: doublet a */
	int dir;
} qi_cand;

/* an operation as the mesh holds it now: its new faces, the old faces they replace, the change of E */
typedef struct qi_op
{
	unsigned int type, c, f, k, v;
	int dir, de, dn;                /* dn: the change of the number of vertices off their ideal */
	unsigned int nvx[QI_FAN * QI_FS], ns[QI_FAN], nn, old[QI_FAN], nold, merged;
	double mx[3];
	qaws_mesh_location ml;
} qi_op;

/* a patch after an operation: vertices (positions, footpoints, free to move), faces (local vertex indices, the mesh face
   or HM_NONE for a new one, the reference's normal under the centre) */
typedef struct qi_local
{
	unsigned int nv, nf;
	unsigned int v[QI_LV];
	double x[QI_LV][3];
	qaws_mesh_location l[QI_LV];
	unsigned char move[QI_LV];
	unsigned int fs[QI_LF], fv[QI_LF][QI_FS], face[QI_LF];
	double nrm[QI_LF][3];
} qi_local;

/* the faces round a vertex, each from the vertex on counter-clockwise; rho the target length */
typedef struct qi_fan
{
	unsigned int n, fs[QI_FAN], id[QI_FAN][QI_FS];
	double const* p[QI_FAN][QI_FS];
	double const* nrm[QI_FAN];
	double rho;
} qi_fan;

typedef struct qi_work
{
	qaws_hm* m;
	qaws_mesh_reference const* r;
	mr_scratch* s;
	unsigned char *kind, *ideal;
	double rho;
	double* fnrm;                   /* per face: the reference's normal under its centre */
	unsigned int *stamp, *fstamp, gen;
	qi_cand* cand;
	unsigned int ncand, ccand, ops;
	qi_local L;
} qi_work;

/* ---- the reference's normals under faces --------------------------------- */

static void qi_normal_at(qi_work* w, double const* const* p, unsigned int n, qaws_mesh_location hint, double* out)
{
	double cen[3] = { 0, 0, 0 };
	qaws_mesh_location l;
	unsigned int i, k;
	for (i = 0; i < n; i++)
		for (k = 0; k < 3; k++)
			cen[k] += p[i][k] / n;
	mr_project(w->r, cen, hint.patch, w->s, &l);
	if (l.patch == HM_NONE)
		l = hint;
	mr_eval(w->r, &l, NULL, out);
}

static void qi_refresh(qi_work* w, unsigned int f)
{
	qaws_hm const* m = w->m;
	double const* p[QI_FS];
	unsigned int i, n = m->fn[f], hint = m->cv[m->fc[f]];
	if (n > QI_FS)
		return;
	for (i = 0; i < n; i++)
	{
		unsigned int u = m->cv[m->fc[f] + i];
		p[i] = m->x + 3 * (size_t)u;
		if (!w->kind[u] && w->kind[hint])
			hint = u;   /* (a free corner's footpoint is in the face's region) */
	}
	qi_normal_at(w, p, n, m->foot[hint], w->fnrm + 3 * (size_t)f);
}

/* face f's smallest scaled Jacobian; down (may be NULL) counts it when it faces against the reference */
static double qi_face_sj(qi_work const* w, unsigned int f, unsigned int* down)
{
	qaws_hm const* m = w->m;
	double const* p[QI_FS];
	unsigned int i;
	if (m->fn[f] > QI_FS)
		return 1;
	for (i = 0; i < m->fn[f]; i++)
		p[i] = m->x + 3 * (size_t)m->cv[m->fc[f] + i];
	if (down)
		*down += !qi_up(p, m->fn[f], w->fnrm + 3 * (size_t)f);
	return qi_sj(p, m->fn[f], w->fnrm + 3 * (size_t)f);
}

/* ---- fans and their moves ------------------------------------------------ */

static int qi_fan_mesh(qi_work const* w, unsigned int v, qi_fan* fan)
{
	qaws_hm const* m = w->m;
	unsigned int c, j, e;
	fan->n = 0;
	fan->rho = w->rho;
	for (c = m->vc[v]; c != HM_NONE; c = hm_round(m, v, c))
	{
		unsigned int h = m->cf[c];
		if (fan->n == QI_FAN || m->fn[h] > QI_FS)
			return 0;
		fan->fs[fan->n] = m->fn[h];
		for (j = 0, e = c; j < m->fn[h]; j++, e = hm_next(m, e))
		{
			fan->id[fan->n][j] = m->cv[e];
			fan->p[fan->n][j] = m->x + 3 * (size_t)m->cv[e];
		}
		fan->nrm[fan->n] = w->fnrm + 3 * (size_t)h;
		fan->n++;
	}
	return fan->n > 0;
}

static int qi_fan_local(qi_local const* L, unsigned int i, double rho, qi_fan* fan)
{
	unsigned int f, j, k;
	fan->n = 0;
	fan->rho = rho;
	for (f = 0; f < L->nf; f++)
	{
		for (k = 0; k < L->fs[f] && L->fv[f][k] != i; k++) {}
		if (k == L->fs[f])
			continue;
		if (fan->n == QI_FAN)
			return 0;
		fan->fs[fan->n] = L->fs[f];
		for (j = 0; j < L->fs[f]; j++)
		{
			unsigned int u = L->fv[f][(k + j) % L->fs[f]];
			fan->id[fan->n][j] = u;
			fan->p[fan->n][j] = L->x[u];
		}
		fan->nrm[fan->n] = L->nrm[f];
		fan->n++;
	}
	return fan->n > 0;
}

/* the fan with its centre at x: its smallest scaled Jacobian, its shortest edge at the centre (len), its quads' smallest ratio
   of shortest to longest edge (ratio), no face turned against the reference (up), its quads with a corner at SJ <= 0 (neg),
   its shape energy (en, qi_fan_energy) */
typedef struct qi_state
{
	double sj, len, ratio, en;
	int up;
	unsigned int neg;
} qi_state;

/* the shape energy the fan's centre moves, the centre at x: the squared deviations of the corners' angles (its own and its
   two neighbours' in each face) from the regular polygons', and of its edges' lengths from rho, in units of rho (the
   angles alone let the quads squeeze toward triangles) */
static double qi_fan_energy(qi_fan const* fan, double const* x)
{
	double e = 0;
	unsigned int i, j;
	for (i = 0; i < fan->n; i++)
	{
		unsigned int s = fan->fs[i], at[3];
		double ideal = QI_PI * (s - 2) / s, d[3], dl;
		qi_sub(fan->p[i][1], x, d);
		dl = sqrt(qi_dot(d, d)) / fan->rho - 1;
		e += QI_LW * dl * dl;
		qi_sub(fan->p[i][s - 1], x, d);
		dl = sqrt(qi_dot(d, d)) / fan->rho - 1;
		e += QI_LW * dl * dl;
		at[0] = 0;
		at[1] = 1;
		at[2] = s - 1;
		for (j = 0; j < 3; j++)
		{
			unsigned int k = at[j];
			double const *o = k ? fan->p[i][k] : x, *u = (k + 1) % s ? fan->p[i][(k + 1) % s] : x, *v = (k + s - 1) % s ? fan->p[i][(k + s - 1) % s] : x;
			double a[3], b[3], l, c;
			qi_sub(u, o, a);
			qi_sub(v, o, b);
			l = sqrt(qi_dot(a, a) * qi_dot(b, b));
			c = l > 0 ? qi_dot(a, b) / l : 1;
			c = acos(c < -1 ? -1 : c > 1 ? 1 : c) - ideal;
			e += c * c;
		}
	}
	return e;
}

static void qi_fan_eval(qi_fan const* fan, double const* x, qi_state* st)
{
	unsigned int i, k;
	st->sj = 1;
	st->len = HUGE_VAL;
	st->ratio = 1;
	st->up = 1;
	st->neg = 0;
	for (i = 0; i < fan->n; i++)
	{
		double const* q[QI_FS];
		double s, e[3], lo = HUGE_VAL, hi = 0;
		for (k = 0; k < fan->fs[i]; k++)
			q[k] = fan->p[i][k];
		q[0] = x;
		for (k = 0; k < fan->fs[i]; k++)
		{
			qi_sub(q[(k + 1) % fan->fs[i]], q[k], e);
			s = qi_dot(e, e);
			lo = fmin(lo, s);
			hi = fmax(hi, s);
			if (k == 0 || k + 1 == fan->fs[i])
				st->len = fmin(st->len, s);
		}
		if (fan->fs[i] == 4)
			st->ratio = fmin(st->ratio, hi > 0 ? sqrt(lo / hi) : 0);
		s = qi_sj(q, fan->fs[i], fan->nrm[i]);
		st->neg += s <= 0 ? 1u : 0u;
		if (s < st->sj)
			st->sj = s;
		st->up &= qi_up(q, fan->fs[i], fan->nrm[i]);
	}
	st->len = sqrt(st->len);
	st->en = qi_fan_energy(fan, x);
}

/* a move from the state s0 to s1 kept: the smallest scaled Jacobian not lower, the shortest edge at the centre not under
   min(its length before, half the target length), the quads' smallest edge ratio not under min(its value before,
   QI_RATIO) (the scaled Jacobian ignores the lengths, and the quads were squeezed toward triangles), no face turned against
   the reference (one turned may turn back), no more quads at SJ <= 0 (the worst lifted, another may not sink), with en
   the shape energy lower */
static int qi_better(qi_state const* s0, qi_state const* s1, double rho, int en)
{
	return s1->sj >= s0->sj && s1->neg <= s0->neg && s1->len >= fmin(s0->len, 0.5 * rho) && s1->ratio >= fmin(s0->ratio, QI_RATIO)
		&& (s1->up || !s0->up) && (!en || s1->en < s0->en);
}

/* Zhou and Shimada's target for the fan's centre: each neighbour's edge to it turned about n, at the neighbour, by half the
   difference of the neighbour's two angles beside it; the mean of the turned ends. 0 when no neighbour has two faces. */
static int qi_fan_target(qi_fan const* fan, double const* n, double* out)
{
	unsigned int id[2 * QI_FAN], i, j, k, cnt = 0;
	double ang[2 * QI_FAN], sum[3] = { 0, 0, 0 };
	double const* at[2 * QI_FAN];
	double const* x = fan->p[0][0];
	for (i = 0; i < fan->n; i++)
	{
		unsigned int s = fan->fs[i];
		double u[3], e[3];
		/* the next neighbour's angle, from its own next to the centre; the previous one's, from the centre to its previous */
		qi_sub(fan->p[i][2], fan->p[i][1], u);
		qi_sub(x, fan->p[i][1], e);
		id[2 * i] = fan->id[i][1];
		ang[2 * i] = qi_turn(u, e, n);
		at[2 * i] = fan->p[i][1];
		qi_sub(x, fan->p[i][s - 1], e);
		qi_sub(fan->p[i][s - 2], fan->p[i][s - 1], u);
		id[2 * i + 1] = fan->id[i][s - 1];
		ang[2 * i + 1] = qi_turn(e, u, n);
		at[2 * i + 1] = fan->p[i][s - 1];
	}
	for (i = 1; i < 2 * fan->n; i += 2)
		for (j = 0; j < 2 * fan->n; j += 2)
			if (id[j] == id[i])
			{
				double b = 0.5 * (ang[i] - ang[j]), c = cos(b), s = sin(b), d[3], nd[3], dn;
				qi_sub(x, at[i], d);
				qi_cross(n, d, nd);
				dn = qi_dot(n, d);
				for (k = 0; k < 3; k++)
					sum[k] += at[i][k] + d[k] * c + nd[k] * s + n[k] * dn * (1 - c);
				cnt++;
				break;
			}
	if (!cnt)
		return 0;
	for (k = 0; k < 3; k++)
		out[k] = sum[k] / cnt;
	return 1;
}

/* a free vertex's move, walked from l: of the angle target and a Newton step down the shape energy (central differences
   along two tangents, the step halved until it lowers the energy), the one qi_better keeps with the lowest energy (0:
   none) */
static int qi_free_move(qi_work* w, qi_fan const* fan, qaws_mesh_location const* l, double* y, qaws_mesh_location* out)
{
	double n[3], t[2][3], d[3], dn, b[2][3], g[2], h = 1e-3 * w->rho, best = HUGE_VAL, a[3] = { 1, 0, 0 }, z[3], gl, lim;
	double const* x = fan->p[0][0];
	unsigned int k, i, nt, tries;
	qi_state s0, st;
	mr_eval(w->r, l, NULL, n);
	qi_fan_eval(fan, x, &s0);
	nt = (unsigned int)qi_fan_target(fan, n, t[0]);
	if (fabs(n[0]) > 0.9)
	{
		a[0] = 0;
		a[1] = 1;
	}
	qi_cross(n, a, b[0]);
	gl = sqrt(qi_dot(b[0], b[0]));
	for (k = 0; k < 3; k++)
		b[0][k] /= gl;
	qi_cross(n, b[0], b[1]);
	for (i = 0; i < 2; i++)
	{
		double ep, em, c;
		for (k = 0; k < 3; k++)
			z[k] = x[k] + h * b[i][k];
		ep = qi_fan_energy(fan, z);
		for (k = 0; k < 3; k++)
			z[k] = x[k] - h * b[i][k];
		em = qi_fan_energy(fan, z);
		c = (ep + em - 2 * s0.en) / (h * h);
		g[i] = (ep - em) / (2 * h);
		g[i] = c > 0 ? g[i] / c : g[i] * 0.01 * w->rho / (fabs(g[i]) + 1e-300);   /* (not convex that way: a hundredth of the length) */
	}
	/* (at most a quarter of the length) */
	gl = sqrt(g[0] * g[0] + g[1] * g[1]);
	lim = gl > 0.25 * w->rho ? 0.25 * w->rho / gl : 1;
	for (tries = 0; gl > 1e-3 * w->rho && tries < 6; tries++, lim *= 0.5)
	{
		for (k = 0; k < 3; k++)
			z[k] = x[k] - lim * (g[0] * b[0][k] + g[1] * b[1][k]);
		if (qi_fan_energy(fan, z) < s0.en)
		{
			memcpy(t[nt++], z, sizeof z);
			break;
		}
	}
	for (i = 0; i < nt; i++)
	{
		qaws_mesh_location to;
		double yy[3];
		qi_sub(t[i], x, d);
		dn = qi_dot(d, n);
		for (k = 0; k < 3; k++)
			d[k] -= dn * n[k];
		if (mr_walk(w->r, *l, d, 1, &to))
			continue;
		mr_eval(w->r, &to, yy, NULL);
		qi_fan_eval(fan, yy, &st);
		if (qi_better(&s0, &st, w->rho, 1) && st.en < best)
		{
			best = st.en;
			memcpy(y, yy, sizeof yy);
			*out = to;
		}
	}
	return best < HUGE_VAL;
}

/* ---- valences ------------------------------------------------------------ */

static int qi_line(qi_work const* w, unsigned int v)
{
	return w->kind[v] || hm_is_border_vertex(w->m, v);
}

static void qi_ideal(qi_work* w)
{
	qaws_hm const* m = w->m;
	unsigned int v, c;
	for (v = 0; v < m->nv; v++)
	{
		double x = 0;
		int id;
		w->ideal[v] = 4;
		if (m->vc[v] == HM_NONE || !qi_line(w, v))
			continue;
		/* the quads the corner angles hold (a corner of an n-gon as n - 2 of n straight angles) */
		for (c = m->vc[v]; c != HM_NONE; c = hm_round(m, v, c))
		{
			unsigned int n = m->fn[m->cf[c]];
			double a[3], b[3], la, lb;
			qi_sub(m->x + 3 * (size_t)hm_to(m, c), m->x + 3 * (size_t)v, a);
			qi_sub(m->x + 3 * (size_t)m->cv[hm_prev(m, c)], m->x + 3 * (size_t)v, b);
			la = sqrt(qi_dot(a, a) * qi_dot(b, b));
			if (la > 0)
			{
				lb = qi_dot(a, b) / la;
				x += acos(lb < -1 ? -1 : lb > 1 ? 1 : lb) / (QI_PI * (n - 2) / n);
			}
		}
		id = (int)floor(x + 0.5);
		w->ideal[v] = (unsigned char)((id < 1 ? 1 : id > 200 ? 200 : id) + hm_is_border_vertex(m, v));
	}
}

static int qi_dev(qi_work const* w, unsigned int v)
{
	return (int)hm_valence(w->m, v) - (int)w->ideal[v];
}

/* the change of E when u's valence moves by dv; acc[0] cleared when it may not (a corner's never changes, a line
   vertex's never drops below its ideal), acc[1] counting the change of the vertices off their ideal */
static int qi_term(qi_work const* w, unsigned int u, int dv, int* acc)
{
	int val = (int)hm_valence(w->m, u), id = w->ideal[u], nv = val + dv;
	if (dv && w->kind[u] == 2)
		acc[0] = 0;
	if (dv < 0 && nv < id && qi_line(w, u))
		acc[0] = 0;
	acc[1] += (nv != id) - (val != id);
	return (nv - id) * (nv - id) - (val - id) * (val - id);
}

/* ---- the operations as the mesh holds them now --------------------------- */

static int qi_op_rotate(qi_work const* w, unsigned int c, int dir, qi_op* o)
{
	qaws_hm const* m = w->m;
	unsigned int h[6], d, s = dir > 0 ? 1 : 5, i, j;
	int acc[2] = { 1, 0 };
	if (!hm_can_rotate(m, c, dir))
		return 0;
	d = m->co[c];
	h[0] = hm_to(m, c);
	h[1] = m->cv[hm_next(m, hm_next(m, c))];
	h[2] = m->cv[hm_prev(m, c)];
	h[3] = m->cv[c];
	h[4] = m->cv[hm_next(m, hm_next(m, d))];
	h[5] = m->cv[hm_prev(m, d)];
	o->de = qi_term(w, h[0], -1, acc) + qi_term(w, h[3], -1, acc) + qi_term(w, h[s], 1, acc) + qi_term(w, h[(s + 3) % 6], 1, acc);
	if (!acc[0])
		return 0;
	o->type = 0;
	o->dn = acc[1];
	o->c = c;
	o->dir = dir;
	for (i = 0; i < 2; i++)
		for (j = 0; j < 4; j++)
			o->nvx[4 * i + j] = h[(s + 3 * i + j) % 6];
	o->ns[0] = o->ns[1] = 4;
	o->nn = 2;
	o->old[0] = m->cf[c];
	o->old[1] = m->cf[d];
	o->nold = 2;
	o->merged = HM_NONE;
	return 1;
}

/* quad f's diagonal from corner k: the end on a line stays (two on lines never merge: shortcut, two ends of one line merged
   would shorten it, when a crease's vertices should be thinned), else the two meet at their middle in the quad's region
   (placed only when place) */
static int qi_op_collapse(qi_work* w, unsigned int f, unsigned int k, int place, qi_op* o)
{
	qaws_hm const* m = w->m;
	unsigned int c0 = m->fc[f], A, B, C, D, c, i;
	int acc[2] = { 1, 0 };
	if (m->fn[f] != 4)
		return 0;
	if (qi_line(w, m->cv[c0 + (k + 2) % 4]))
		k = (k + 2) % 4;
	A = m->cv[c0 + k];
	B = m->cv[c0 + (k + 1) % 4];
	C = m->cv[c0 + (k + 2) % 4];
	D = m->cv[c0 + (k + 3) % 4];
	if (qi_line(w, C) || !hm_can_collapse_diagonal(m, f, k))
		return 0;
	o->de = qi_term(w, A, (int)hm_valence(m, C) - 2, acc) - qi_dev(w, C) * qi_dev(w, C) + qi_term(w, B, -1, acc) + qi_term(w, D, -1, acc);
	if (!acc[0])
		return 0;
	o->type = 1;
	o->dn = acc[1] - (qi_dev(w, C) != 0);
	o->dir = 0;
	o->f = f;
	o->k = k;
	o->merged = A;
	o->nn = o->nold = 0;
	/* C's faces but f, A for C */
	for (c = m->vc[C], i = 0; c != HM_NONE; c = hm_round(m, C, c))
	{
		unsigned int h = m->cf[c], e;
		if (o->nold == QI_FAN || m->fn[h] > QI_FS)
			return 0;
		o->old[o->nold++] = h;
		if (h == f)
			continue;
		for (e = m->fc[h]; e < m->fc[h] + m->fn[h]; e++)
			o->nvx[i++] = m->cv[e] == C ? A : m->cv[e];
		o->ns[o->nn++] = m->fn[h];
	}
	if (!place)
		return 1;
	if (qi_line(w, A))
	{
		memcpy(o->mx, m->x + 3 * (size_t)A, sizeof o->mx);
		o->ml = m->foot[A];
	}
	else
	{
		double mid[3];
		unsigned int t = m->foot[A].patch;
		for (i = 0; i < 3; i++)
			mid[i] = 0.5 * (m->x[3 * (size_t)A + i] + m->x[3 * (size_t)C + i]);
		mr_project_region(w->r, mid, t, t < w->r->nt ? w->r->treg[t] : HM_NONE, w->s, &o->ml);
		if (o->ml.patch == HM_NONE)
			mr_project(w->r, mid, HM_NONE, w->s, &o->ml);
		if (o->ml.patch == HM_NONE)
			return 0;
		mr_eval(w->r, &o->ml, o->mx, NULL);
	}
	return 1;
}

/* as hm_remove_doublet checks it: f = v a x b and g = v b y a become a x b y */
static int qi_op_doublet(qi_work const* w, unsigned int v, qi_op* o)
{
	qaws_hm const* m = w->m;
	unsigned int c, d, a, b, x, y;
	int acc[2] = { 1, 0 };
	if (m->vc[v] == HM_NONE || (m->vflag[v] & HM_CORNER) || qi_line(w, v) || hm_valence(m, v) != 2 || hm_feature_degree(m, v))
		return 0;
	c = m->vc[v];
	d = hm_ccw(m, c);
	if (d == HM_NONE || m->cf[c] == m->cf[d] || m->fn[m->cf[c]] != 4 || m->fn[m->cf[d]] != 4)
		return 0;
	a = hm_to(m, c);
	x = m->cv[hm_next(m, hm_next(m, c))];
	b = m->cv[hm_prev(m, c)];
	y = m->cv[hm_next(m, hm_next(m, d))];
	if (hm_to(m, d) != b || m->cv[hm_prev(m, d)] != a || a == b || x == y || hm_valence(m, a) < 3 || hm_valence(m, b) < 3)
		return 0;
	o->de = -qi_dev(w, v) * qi_dev(w, v) + qi_term(w, a, -1, acc) + qi_term(w, b, -1, acc);
	if (!acc[0])
		return 0;
	o->type = 2;
	o->dn = acc[1] - (qi_dev(w, v) != 0);
	o->dir = 0;
	o->v = v;
	o->nvx[0] = a;
	o->nvx[1] = x;
	o->nvx[2] = b;
	o->nvx[3] = y;
	o->ns[0] = 4;
	o->nn = 1;
	o->old[0] = m->cf[c];
	o->old[1] = m->cf[d];
	o->nold = 2;
	o->merged = HM_NONE;
	return 1;
}

/* ---- trying an operation on a copy of its patch -------------------------- */

static unsigned int qi_lv(qi_work const* w, qi_local* L, unsigned int u)
{
	unsigned int i;
	for (i = 0; i < L->nv; i++)
		if (L->v[i] == u)
			return i;
	if (L->nv == QI_LV)
		return HM_NONE;
	L->v[i] = u;
	memcpy(L->x[i], w->m->x + 3 * (size_t)u, sizeof L->x[i]);
	L->l[i] = w->m->foot[u];
	L->move[i] = 0;
	L->nv++;
	return i;
}

/* the patch o leaves: its new faces, then every other face round their vertices but the old ones; 0 when it does not fit.
   shortcut: only the new faces' vertices are relaxed in it (their neighbours stay), a ring wider when the patch should
   settle in full before it is judged. */
static int qi_patch(qi_work* w, qi_op const* o, qi_local* L)
{
	qaws_hm const* m = w->m;
	unsigned int i, j, at = 0, n0, c;
	L->nv = L->nf = 0;
	for (i = 0; i < o->nn; at += o->ns[i], i++)
	{
		double const* q[QI_FS];
		if (L->nf == QI_LF)
			return 0;
		for (j = 0; j < o->ns[i]; j++)
		{
			unsigned int k = qi_lv(w, L, o->nvx[at + j]);
			if (k == HM_NONE)
				return 0;
			if (o->nvx[at + j] == o->merged)
			{
				memcpy(L->x[k], o->mx, sizeof L->x[k]);
				L->l[k] = o->ml;
			}
			L->fv[i][j] = k;
			q[j] = L->x[k];
		}
		L->fs[i] = o->ns[i];
		L->face[i] = HM_NONE;
		L->nf++;
		qi_normal_at(w, q, o->ns[i], L->l[L->fv[i][0]], L->nrm[i]);
	}
	n0 = L->nv;
	for (i = 0; i < n0; i++)
	{
		unsigned int u = L->v[i];
		L->move[i] = (unsigned char)!qi_line(w, u);
		for (c = m->vc[u]; c != HM_NONE; c = hm_round(m, u, c))
		{
			unsigned int h = m->cf[c], e;
			for (j = 0; j < o->nold && o->old[j] != h; j++) {}
			if (j < o->nold)
				continue;
			for (j = 0; j < L->nf && L->face[j] != h; j++) {}
			if (j < L->nf)
				continue;
			if (L->nf == QI_LF || m->fn[h] > QI_FS)
				return 0;
			for (j = 0, e = m->fc[h]; e < m->fc[h] + m->fn[h]; e++, j++)
			{
				unsigned int k = qi_lv(w, L, m->cv[e]);
				if (k == HM_NONE)
					return 0;
				L->fv[L->nf][j] = k;
			}
			L->fs[L->nf] = m->fn[h];
			L->face[L->nf] = h;
			memcpy(L->nrm[L->nf], w->fnrm + 3 * (size_t)h, sizeof L->nrm[0]);
			L->nf++;
		}
	}
	return 1;
}

static int qi_commit(qi_work* w, qi_op const* o, qi_local const* L)
{
	qaws_hm* m = w->m;
	unsigned int i, c;
	int done = 0;
	if (o->type == 0)
		done = hm_rotate(m, o->c, o->dir);
	else if (o->type == 1)
	{
		unsigned int A = m->cv[m->fc[o->f] + o->k], C = m->cv[m->fc[o->f] + (o->k + 2) % 4];
		done = hm_collapse_diagonal(m, o->f, o->k, o->mx);
		if (done && w->kind[C] > w->kind[A])
			w->kind[A] = w->kind[C];
	}
	else
		done = hm_remove_doublet(m, o->v);
	if (!done)
		return 0;
	for (i = 0; i < L->nv; i++)
		if (L->move[i] || L->v[i] == o->merged)
		{
			memcpy(m->x + 3 * (size_t)L->v[i], L->x[i], sizeof L->x[i]);
			m->foot[L->v[i]] = L->l[i];
		}
	for (i = 0; i < L->nv; i++)
		for (c = m->vc[L->v[i]]; c != HM_NONE; c = hm_round(m, L->v[i], c))
			qi_refresh(w, m->cf[c]);
	w->ops++;
	return 1;
}

/* the second of a pair: the quads g its first rotation made, and their smallest scaled Jacobian and the count of them
   facing against the reference before it */
typedef struct qi_pre
{
	unsigned int g[2], down;
	double sj;
} qi_pre;

static int qi_in(unsigned int const* a, unsigned int n, unsigned int x)
{
	unsigned int i;
	for (i = 0; i < n && a[i] != x; i++) {}
	return i < n;
}

/* o evaluated on a copy of its patch (its free vertices relaxed there three times), made when E drops with no more
   vertices off their ideal, the patch's smallest scaled Jacobian stays at least max(0.3, its value before - 0.1) (where
   it was under 0.3: no lower) and no more of its faces turn against the reference than before; for the second of a pair
   (pre), before is before the first rotation and its quads are judged after too. On a patch holding a quad turned inside
   out (SJ <= 0 before) E may stay or rise: the shape alone judges it, and its smallest scaled Jacobian must then rise. */
static int qi_try(qi_work* w, qi_op const* o, qi_pre const* pre)
{
	qi_local* L = &w->L;
	double before = 1, after = 1;
	unsigned int i, f, sweep, down0 = 0, down1 = 0;
	int shape;
	if (!qi_patch(w, o, L))
		return 0;
	for (i = 0; i < o->nold; i++)
		if (!pre || !qi_in(pre->g, 2, o->old[i]))
			before = fmin(before, qi_face_sj(w, o->old[i], &down0));
	for (f = 0; f < L->nf; f++)
		if (L->face[f] != HM_NONE && (!pre || !qi_in(pre->g, 2, L->face[f])))
			before = fmin(before, qi_face_sj(w, L->face[f], &down0));
	if (pre)
	{
		before = fmin(before, pre->sj);
		down0 += pre->down;
		for (i = 0; i < 2; i++)
			if (!qi_in(o->old, o->nold, pre->g[i]) && !qi_in(L->face, L->nf, pre->g[i]))
				after = fmin(after, qi_face_sj(w, pre->g[i], &down1));
	}
	shape = o->de >= 0 || o->dn > 0;
	if (shape && before > 0)
		return 0;
	for (sweep = 0; sweep < 3; sweep++)
		for (i = 0; i < L->nv; i++)
		{
			qi_fan fan;
			qaws_mesh_location l;
			double y[3];
			if (!L->move[i] || !qi_fan_local(L, i, w->rho, &fan) || !qi_free_move(w, &fan, &L->l[i], y, &l))
				continue;
			memcpy(L->x[i], y, sizeof y);
			L->l[i] = l;
		}
	for (f = 0; f < L->nf; f++)
	{
		double const* q[QI_FS];
		for (i = 0; i < L->fs[f]; i++)
			q[i] = L->x[L->fv[f][i]];
		after = fmin(after, qi_sj(q, L->fs[f], L->nrm[f]));
		down1 += !qi_up(q, L->fs[f], L->nrm[f]);
	}
	if (down1 > down0 || after < (before < 0.3 ? before : fmax(0.3, before - 0.1)) || (shape && before <= 0 && after <= before))
		return 0;
	return qi_commit(w, o, L);
}

/* ---- passes ------------------------------------------------------------- */

static int qi_push(qi_work* w, qi_op const* o, unsigned int a, unsigned int b)
{
	if (w->ncand == w->ccand)
	{
		unsigned int c = w->ccand ? 2 * w->ccand : 256;
		qi_cand* p = (qi_cand*)qaws_internal_realloc(w->cand, sizeof(qi_cand) * c);
		if (!p)
			return 0;
		w->cand = p;
		w->ccand = c;
	}
	w->cand[w->ncand].de = o->de;
	w->cand[w->ncand].type = o->type;
	w->cand[w->ncand].a = a;
	w->cand[w->ncand].b = b;
	w->cand[w->ncand].dir = o->dir;
	w->ncand++;
	return 1;
}

static int qi_cand_cmp(void const* x, void const* y)
{
	qi_cand const *p = (qi_cand const*)x, *q = (qi_cand const*)y;
	if (p->de != q->de)
		return p->de < q->de ? -1 : 1;
	return p->a < q->a ? -1 : p->a > q->a ? 1 : 0;
}

/* a candidate as the mesh holds it now, placed */
static int qi_rebuild(qi_work* w, qi_cand const* q, qi_op* o)
{
	qaws_hm const* m = w->m;
	unsigned int c;
	if ((m->vflag[q->a] & HM_DELETED) || m->vc[q->a] == HM_NONE)
		return 0;
	if (q->type == 0)
	{
		c = hm_find(m, q->a, q->b);
		return c != HM_NONE && qi_op_rotate(w, c, q->dir, o);
	}
	if (q->type == 2)
		return qi_op_doublet(w, q->a, o);
	for (c = m->vc[q->a]; c != HM_NONE; c = hm_round(m, q->a, c))
	{
		unsigned int h = m->cf[c];
		if (m->fn[h] == 4 && m->cv[hm_next(m, hm_next(m, c))] == q->b)
			return qi_op_collapse(w, h, c - m->fc[h], 1, o);
	}
	return 0;
}

/* the single operations round the vertices off their ideal (and their neighbours) that lower E, best first (ties by
   vertex), each checked again as the mesh holds it when its turn comes; the number made, HM_NONE on allocation failure.
   shortcut: collected once per pass, not again round each commit, so what a commit opens waits for the next pass (at most
   rounds); a queue fed after each commit when single operations commit often. */
static unsigned int qi_pass(qi_work* w)
{
	qaws_hm const* m = w->m;
	unsigned int v, c, i, f, k, done = 0;
	qi_op o;
	w->ncand = 0;
	w->gen++;
	for (v = 0; v < m->nv; v++)
		if (m->vc[v] != HM_NONE && qi_dev(w, v))
		{
			w->stamp[v] = w->gen;
			for (c = m->vc[v]; c != HM_NONE; c = hm_round(m, v, c))
			{
				w->stamp[hm_to(m, c)] = w->gen;
				w->stamp[m->cv[hm_prev(m, c)]] = w->gen;   /* (the border neighbour too) */
			}
		}
	for (v = 0; v < m->nv; v++)
	{
		if (w->stamp[v] != w->gen || m->vc[v] == HM_NONE)
			continue;
		for (c = m->vc[v]; c != HM_NONE; c = hm_round(m, v, c))
		{
			unsigned int u = hm_to(m, c), h = m->cf[c];
			int dir;
			if (m->co[c] != HM_NONE && !(w->stamp[u] == w->gen && u < v))
				for (dir = -1; dir <= 1; dir += 2)
					if (qi_op_rotate(w, c, dir, &o) && o.de < 0 && !qi_push(w, &o, v, u))
						return HM_NONE;
			if (w->fstamp[h] == w->gen)
				continue;
			w->fstamp[h] = w->gen;
			for (k = 0; k < 2; k++)
				if (qi_op_collapse(w, h, k, 0, &o) && o.de < 0
					&& !qi_push(w, &o, m->cv[m->fc[h] + o.k], m->cv[m->fc[h] + (o.k + 2) % 4]))
					return HM_NONE;
		}
		if (qi_op_doublet(w, v, &o) && o.de < 0 && !qi_push(w, &o, v, v))
			return HM_NONE;
	}
	/* the quads turned inside out: their diagonals collapsed whatever E says (qi_try judges those by shape) */
	for (f = 0; f < m->nf; f++)
		if (!(m->fflag[f] & HM_DELETED) && m->fn[f] == 4 && qi_face_sj(w, f, NULL) <= 0)
			for (k = 0; k < 2; k++)
				if (qi_op_collapse(w, f, k, 0, &o) && (o.de >= 0 || w->fstamp[f] != w->gen)
					&& !qi_push(w, &o, m->cv[m->fc[f] + o.k], m->cv[m->fc[f] + (o.k + 2) % 4]))
					return HM_NONE;
	qsort(w->cand, w->ncand, sizeof(qi_cand), qi_cand_cmp);
	for (i = 0; i < w->ncand; i++)
		if (qi_rebuild(w, &w->cand[i], &o))
			done += (unsigned int)qi_try(w, &o, NULL);
	return done;
}

/* after a first rotation (its quads pre->g, its new diagonal p q, its change of E de1): a second operation near it that lowers E
   with it, made: a rotation of an edge of the two quads, a diagonal collapse of a face round their corners, a doublet there */
static int qi_second(qi_work* w, qi_pre const* pre, unsigned int p, unsigned int q, int de1, int dn1)
{
	qaws_hm* m = w->m;
	unsigned int side, e, c, k;
	unsigned int const* g = pre->g;
	qi_op o;
	for (side = 0; side < 2; side++)
		for (e = m->fc[g[side]]; e < m->fc[g[side]] + 4; e++)
		{
			unsigned int u = m->cv[e];
			int d2;
			if ((u != p || hm_to(m, e) != q) && (u != q || hm_to(m, e) != p))
				for (d2 = -1; d2 <= 1; d2 += 2)
					if (qi_op_rotate(w, e, d2, &o) && de1 + o.de < 0 && dn1 + o.dn <= 0)
					{
						o.de += de1;
						o.dn += dn1;
						if (qi_try(w, &o, pre))
							return 1;
					}
			for (c = m->vc[u]; c != HM_NONE; c = hm_round(m, u, c))
				for (k = 0; k < 2; k++)
					if (qi_op_collapse(w, m->cf[c], k, 0, &o) && de1 + o.de < 0 && dn1 + o.dn <= 0 && qi_op_collapse(w, m->cf[c], k, 1, &o))
					{
						o.de += de1;
						o.dn += dn1;
						if (qi_try(w, &o, pre))
							return 1;
					}
			if (qi_op_doublet(w, u, &o) && de1 + o.de < 0 && dn1 + o.dn <= 0)
			{
				o.de += de1;
				o.dn += dn1;
				if (qi_try(w, &o, pre))
					return 1;
			}
		}
	return 0;
}

/* Two operations round each vertex still off its ideal where none alone helps: an edge at it turned, then a second
   operation beside the two quads that made (a singularity stepped along a grid line onto its partner); the first turned
   back when none goes with it. The number of pairs made, HM_NONE when a turn back fails (no room). */
static unsigned int qi_pairs(qi_work* w)
{
	qaws_hm* m = w->m;
	unsigned int v, done = 0;
	for (v = 0; v < m->nv; v++)
	{
		unsigned int ring[QI_FAN], n, i, c;
		int dir, made = 0;
		if (m->vc[v] == HM_NONE || !qi_dev(w, v))
			continue;
		n = hm_ring(m, v, ring, QI_FAN);
		for (i = 0; i < n && i < QI_FAN && !made; i++)
			for (dir = -1; dir <= 1 && !made; dir += 2)
			{
				qi_op o1;
				unsigned int p, q, j;
				qi_pre pre;
				double keep[6];
				c = hm_find(m, v, ring[i]);
				if (c == HM_NONE || !qi_op_rotate(w, c, dir, &o1))
					continue;
				p = o1.nvx[0];
				q = o1.nvx[3];
				/* (turned back only when p and q keep three edges with the one they gain) */
				if (hm_valence(m, p) + 1 < 4u - (unsigned int)hm_is_border_vertex(m, p) || hm_valence(m, q) + 1 < 4u - (unsigned int)hm_is_border_vertex(m, q))
					continue;
				pre.g[0] = o1.old[0];
				pre.g[1] = o1.old[1];
				pre.down = 0;
				pre.sj = fmin(qi_face_sj(w, pre.g[0], &pre.down), qi_face_sj(w, pre.g[1], &pre.down));
				memcpy(keep, w->fnrm + 3 * (size_t)pre.g[0], sizeof(double) * 3);
				memcpy(keep + 3, w->fnrm + 3 * (size_t)pre.g[1], sizeof(double) * 3);
				if (!hm_rotate(m, c, dir))
					continue;
				qi_refresh(w, pre.g[0]);
				qi_refresh(w, pre.g[1]);
				made = qi_second(w, &pre, p, q, o1.de, o1.dn);
				if (made)
				{
					w->ops++;
					done++;
					break;
				}
				/* (the turn back restores the faces at their indices) */
				c = hm_find(m, q, p);
				if (c == HM_NONE || !hm_rotate(m, c, -dir))
					return HM_NONE;
				for (j = 0; j < 2; j++)
					memcpy(w->fnrm + 3 * (size_t)pre.g[j], keep + 3 * j, sizeof(double) * 3);
			}
	}
	return done;
}

/* ---- relaxation --------------------------------------------------------- */

/* the two neighbours of line vertex v along its feature or border edges; 0 unless exactly two */
static int qi_line_neighbours(qaws_hm const* m, unsigned int v, unsigned int* nb)
{
	unsigned int c, nn = 0, j;
	for (c = m->vc[v]; c != HM_NONE && nn <= 2; c = hm_round(m, v, c))
	{
		unsigned int p = hm_prev(m, c), u[2];
		u[0] = (m->cflag[c] & HM_FEATURE) || m->co[c] == HM_NONE ? hm_to(m, c) : HM_NONE;
		u[1] = (m->cflag[p] & HM_FEATURE) || m->co[p] == HM_NONE ? m->cv[p] : HM_NONE;
		for (j = 0; j < 2; j++)
			if (u[j] != HM_NONE && (nn == 0 || nb[0] != u[j]) && (nn < 2 || nb[1] != u[j]))
			{
				if (nn < 2)
					nb[nn] = u[j];
				nn++;
			}
	}
	return nn == 2;
}

/* the best of eight tangent steps from v, halving from a quarter of the length to a thousandth, by the smallest scaled
   Jacobian round it (its shortest edge, edge ratio and quads at SJ <= 0 held as qi_better holds them) */
static void qi_search(qi_work* w, unsigned int v)
{
	qaws_hm* m = w->m;
	qi_fan fan;
	qi_state st;
	double best, step, n[3], t1[3], t2[3], a[3] = { 1, 0, 0 }, l, len, ratio;
	unsigned int k, tries, neg;
	if (!qi_fan_mesh(w, v, &fan))
		return;
	qi_fan_eval(&fan, m->x + 3 * (size_t)v, &st);
	best = st.sj;
	len = fmin(st.len, 0.5 * w->rho);
	ratio = fmin(st.ratio, QI_RATIO);
	neg = st.neg;
	for (step = 0.25 * w->rho; step >= 1e-3 * w->rho; step *= 0.5)
		for (tries = 0; tries < 4; tries++)
		{
			int moved = 0;
			mr_eval(w->r, &m->foot[v], NULL, n);
			if (fabs(n[0]) > 0.9)
			{
				a[0] = 0;
				a[1] = 1;
			}
			qi_cross(n, a, t1);
			l = sqrt(qi_dot(t1, t1));
			for (k = 0; k < 3; k++)
				t1[k] /= l;
			qi_cross(n, t1, t2);
			for (k = 0; k < 8 && !moved; k++)
			{
				double d[3], y[3], ca = cos(k * QI_PI / 4), sa = sin(k * QI_PI / 4);
				qaws_mesh_location to;
				unsigned int i;
				for (i = 0; i < 3; i++)
					d[i] = step * (ca * t1[i] + sa * t2[i]);
				if (mr_walk(w->r, m->foot[v], d, 1, &to))
					continue;
				mr_eval(w->r, &to, y, NULL);
				qi_fan_eval(&fan, y, &st);
				if (st.sj > best + 1e-9 && st.up && st.neg <= neg && st.len >= len && st.ratio >= ratio)
				{
					best = st.sj;
					memcpy(m->x + 3 * (size_t)v, y, sizeof y);
					m->foot[v] = to;
					moved = 1;
				}
			}
			if (!moved)
				break;
		}
}

typedef struct qi_bad
{
	double sj;
	unsigned int f;
} qi_bad;

static int qi_bad_cmp(void const* x, void const* y)
{
	double a = ((qi_bad const*)x)->sj, b = ((qi_bad const*)y)->sj;
	return a < b ? -1 : a > b ? 1 : 0;
}

/* sweeps over the vertices whose neighbourhood moved (at most 30, until no move reaches 0.3% of the length), then the
   worst quads; 0 on allocation failure */
static int qi_relax(qi_work* w)
{
	qaws_hm* m = w->m;
	unsigned int sweep, v, f, nb = 0, i, k;
	unsigned char* busy = (unsigned char*)qaws_internal_malloc((size_t)m->nv + 1);
	qi_bad* bad = (qi_bad*)qaws_internal_malloc(sizeof(qi_bad) * (m->nf + 1));
	if (!busy || !bad)
	{
		qaws_internal_free(busy);
		qaws_internal_free(bad);
		return 0;
	}
	memset(busy, 1, m->nv);
	for (sweep = 0; sweep < 30; sweep++)
	{
		int moved = 0;
		for (v = 0; v < m->nv; v++)
		{
			qi_fan fan;
			qaws_mesh_location l;
			double y[3], d[3];
			if (!busy[v])
				continue;
			busy[v] = 0;
			if (m->vc[v] == HM_NONE || w->kind[v] == 2 || !qi_fan_mesh(w, v, &fan))
				continue;
			if (qi_line(w, v))
			{
				/* along its line, to the middle of its two neighbours there */
				unsigned int ln[2];
				double mid[3];
				qi_state s0, s1;
				if (!qi_line_neighbours(m, v, ln))
					continue;
				for (k = 0; k < 3; k++)
					mid[k] = 0.5 * (m->x[3 * (size_t)ln[0] + k] + m->x[3 * (size_t)ln[1] + k]);
				if (mr_project_crease(w->r, mid, m->foot[v].patch, w->s, &l) == HUGE_VAL || l.patch == HM_NONE)
					continue;
				mr_eval(w->r, &l, y, NULL);
				qi_fan_eval(&fan, m->x + 3 * (size_t)v, &s0);
				qi_fan_eval(&fan, y, &s1);
				if (!qi_better(&s0, &s1, w->rho, 0))
					continue;
			}
			else if (!qi_free_move(w, &fan, &m->foot[v], y, &l))
				continue;
			qi_sub(y, m->x + 3 * (size_t)v, d);
			memcpy(m->x + 3 * (size_t)v, y, sizeof y);
			m->foot[v] = l;
			if (qi_dot(d, d) < 1e-5 * w->rho * w->rho)
				continue;
			/* its neighbourhood again next sweep */
			moved = 1;
			for (i = 0; i < fan.n; i++)
				for (k = 0; k < fan.fs[i]; k++)
					busy[fan.id[i][k]] = 1;
		}
		if (!moved)
			break;
	}
	/* the quads under 0.5, worst first */
	for (f = 0; f < m->nf; f++)
		if (!(m->fflag[f] & HM_DELETED) && m->fn[f] == 4)
		{
			bad[nb].sj = qi_face_sj(w, f, NULL);
			bad[nb].f = f;
			nb += bad[nb].sj < 0.5;
		}
	qsort(bad, nb, sizeof(qi_bad), qi_bad_cmp);
	for (i = 0; i < nb; i++)
		for (k = 0; k < 4; k++)
		{
			v = m->cv[m->fc[bad[i].f] + k];
			if (!qi_line(w, v))
				qi_search(w, v);
		}
	qaws_internal_free(busy);
	qaws_internal_free(bad);
	return 1;
}

unsigned int qi_improve(qaws_hm* m, qaws_mesh_reference const* r, mr_scratch* s, unsigned char* kind, double rho, unsigned int rounds)
{
	qi_work w;
	unsigned int f, v, pass, a, ok = 0;
	memset(&w, 0, sizeof w);
	w.m = m;
	w.r = r;
	w.s = s;
	w.kind = kind;
	w.rho = rho;
	w.ideal = (unsigned char*)qaws_internal_malloc((size_t)m->nv + 1);
	w.stamp = (unsigned int*)qaws_internal_zalloc((size_t)m->nv + 1, sizeof(unsigned int));
	w.fstamp = (unsigned int*)qaws_internal_zalloc((size_t)m->nf + 1, sizeof(unsigned int));
	w.fnrm = (double*)qaws_internal_malloc(sizeof(double) * (3 * (size_t)m->nf + 3));
	if (!w.ideal || !w.stamp || !w.fstamp || !w.fnrm)
		goto done;
	qi_ideal(&w);
	for (f = 0; f < m->nf; f++)
		if (!(m->fflag[f] & HM_DELETED))
			qi_refresh(&w, f);
	for (pass = 0; pass < rounds; pass++)
	{
		unsigned int n = qi_pass(&w), n2;
		if (n == HM_NONE)
			goto done;
		n2 = qi_pairs(&w);
		if (n2 == HM_NONE)
			goto done;
		if (!n && !n2)
			break;
	}
	if (!qi_relax(&w))
		goto done;
	/* the kinds carried through compaction as a vertex attribute */
	a = hm_attr_add(m, NULL, QAWS_MESH_VERTEX, QAWS_ATTRIBUTE_NEAREST, 1);
	if (a == HM_NONE)
		goto done;
	for (v = 0; v < m->nv; v++)
		hm_attr_at(m, a, v)[0] = kind[v];
	hm_compact(m);
	for (v = 0; v < m->nv; v++)
		kind[v] = (unsigned char)hm_attr_at(m, a, v)[0];
	hm_attr_remove(m, a);
	ok = 1;
done:
	qaws_internal_free(w.ideal);
	qaws_internal_free(w.stamp);
	qaws_internal_free(w.fstamp);
	qaws_internal_free(w.fnrm);
	qaws_internal_free(w.cand);
	return ok ? w.ops : HM_NONE;
}
