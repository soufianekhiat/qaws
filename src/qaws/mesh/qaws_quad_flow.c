/*
 * The lattice steps of INSTANT quads made consistent by a min-cost flow
 * (qaws_quad_flow.h).
 *
 * Face f's equation is its closure sum S_f, turned into a combed frame:
 * E_f = R^g_f S_f, the g_f chosen breadth first across the faces so that
 * the step x of an edge between faces F and G enters E_F as +M x and E_G as
 * -M x. Each component of x then enters two rows with opposite signs: a
 * unit of it is a unit of flow between those rows (nodes). Supplies are
 * -E; edges with one equation (the border, a singular face's) lead to a
 * sink. Where the combing cannot agree (round a singular face), a
 * component entering both rows with one sign is kept as it is. Then odd
 * cones are paired by unit paths, each cone moved onto a vertex of its
 * face, and the lattice untangled where the changes turned faces over.
 */

#include "qaws_quad_flow.h"
#include "../internal/qaws_internal_memory.h"
#include "../internal/qaws_internal_flow.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

static double qf_dot(double const* a, double const* b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void qf_cross(double const* a, double const* b, double* o)
{
	o[0] = a[1] * b[2] - a[2] * b[1];
	o[1] = a[2] * b[0] - a[0] * b[2];
	o[2] = a[0] * b[1] - a[1] * b[0];
}

void qf_turn(int const* c, int k, int* o)
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

void qf_steps(qaws_hm const* dm, mf_level const* L, double rho, int* T, unsigned char* K, double* Tc)
{
	unsigned int c, k;
	for (c = 0; c < dm->nc; c++)
	{
		unsigned int a, b, o = dm->co[c];
		double ob[3], xb[3], d[3], g[3], t[3], re;
		if (dm->cf[c] == HM_NONE || (o != HM_NONE && o < c))
			continue;
		a = dm->cv[c];
		b = hm_to(dm, c);
		re = mf_pair_scale(L, a, b, rho);
		K[c] = (unsigned char)mf_rotation(L->q + 3 * (size_t)a, L->nrm + 3 * (size_t)a, L->q + 3 * (size_t)b, L->nrm + 3 * (size_t)b);
		mf_unfold(L->x + 3 * (size_t)a, L->nrm + 3 * (size_t)a, L->x + 3 * (size_t)b, L->nrm + 3 * (size_t)b, L->o + 3 * (size_t)b, ob);
		mf_unfold(L->x + 3 * (size_t)a, L->nrm + 3 * (size_t)a, L->x + 3 * (size_t)b, L->nrm + 3 * (size_t)b, L->x + 3 * (size_t)b, xb);
		for (k = 0; k < 3; k++)
			d[k] = ob[k] - L->o[3 * (size_t)a + k];
		qf_cross(L->nrm + 3 * (size_t)a, L->q + 3 * (size_t)a, t);
		T[2 * (size_t)c] = (int)floor(qf_dot(d, L->q + 3 * (size_t)a) / re + 0.5);
		T[2 * (size_t)c + 1] = (int)floor(qf_dot(d, t) / re + 0.5);
		/* each end's place in its cell, in a's frame (the lattice points' own step is whole but for noise) */
		for (k = 0; k < 3; k++)
		{
			d[k] = L->x[3 * (size_t)a + k] - L->o[3 * (size_t)a + k];
			g[k] = xb[k] - ob[k];
		}
		Tc[4 * (size_t)c] = qf_dot(d, L->q + 3 * (size_t)a) / re;
		Tc[4 * (size_t)c + 1] = qf_dot(d, t) / re;
		Tc[4 * (size_t)c + 2] = qf_dot(g, L->q + 3 * (size_t)a) / re;
		Tc[4 * (size_t)c + 3] = qf_dot(g, t) / re;
		if (o != HM_NONE)
		{
			int m[2];
			m[0] = -T[2 * (size_t)c];
			m[1] = -T[2 * (size_t)c + 1];
			K[o] = (unsigned char)((4 - K[c]) & 3);
			qf_turn(m, K[c], T + 2 * (size_t)o);
		}
	}
}

int qf_closure(qaws_hm const* dm, unsigned int f, int const* T, unsigned char const* K, unsigned char* fr, int* S)
{
	unsigned int c0 = dm->fc[f], i;
	int turn = 0, s[2];
	S[0] = S[1] = 0;
	for (i = 0; i < 3; i++)
	{
		fr[c0 + i] = (unsigned char)(turn & 3);
		qf_turn(T + 2 * (size_t)(c0 + i), turn, s);
		S[0] += s[0];
		S[1] += s[1];
		turn -= K[c0 + i];
	}
	return turn & 3;
}

int qf_axis(qaws_hm const* dm, mf_level const* L, unsigned int c)
{
	unsigned int a = dm->cv[c], b = hm_to(dm, c), k;
	double d[3], t[3];
	for (k = 0; k < 3; k++)
		d[k] = L->x[3 * (size_t)b + k] - L->x[3 * (size_t)a + k];
	qf_cross(L->nrm + 3 * (size_t)a, L->q + 3 * (size_t)a, t);
	return fabs(qf_dot(d, L->q + 3 * (size_t)a)) >= fabs(qf_dot(d, t)) ? 0 : 1;
}

/* the cost of moving a step t by s, its ends at fa and fb in their cells (in the step's frame): how much farther
   from its own lattice point the end that moves at less cost gets (nothing when nearer), plus a little for every
   change; a step of two (a lattice point left without a vertex) and a second unit ten times dearer, a last resort.
   In millionths. */
static long long qf_cost(int t, double fa, double fb, int s, int second)
{
	double d = fabs(fb - s) - fabs(fb), e = fabs(fa + s) - fabs(fa);
	d = d < e ? d : e;
	if (d < 0)
		d = 0;
	if (second || abs(t + s) >= 2)
		d = 10 * (1 + d);
	return (long long)((d + 1e-3) * 1e6 + 0.5);
}

/* The network: nodes 2 f + row for a face where the field does not turn, 2 f for a singular face, the sink Z = 2 nf (the
   border; closing, the cones too, at a price). rec per arc: corner * 4 + component * 2 + (1: lowers the step), HM_NONE
   for a cone's arcs to the sink. */
#define QF_CONE 3000000
#define QF_SLIP 3000000                 /* a dislocation left where the lattice's spacing changes (its extra row ends there) */
#define QF_GRADED 0.02                  /* across a face whose spacings differ by more (relatively) */
typedef struct qf_net
{
	qaws_flow* flow;
	unsigned int* rec;
	unsigned int na, cut;
} qf_net;

static void qf_net_free(qf_net* n)
{
	qaws_flow_destroy(n->flow);
	qaws_internal_free(n->rec);
	memset(n, 0, sizeof *n);
}

/* where component j of an edge's step x enters face f's equation, f's sum holding R^m (s x): the node, and the sign
   (0: a free node, the sink or a singular face's own, which takes either) */
static unsigned int qf_node(unsigned int f, unsigned char const* sing, unsigned int Z, int m, int j, int s, int* sign)
{
	int u[2] = { 0, 0 }, w[2];
	*sign = 0;
	if (f == HM_NONE)
		return Z;
	if (sing[f])
		return 2 * f;
	u[j] = s;
	qf_turn(u, m, w);
	*sign = w[0] ? w[0] : w[1];
	return 2 * f + (w[0] ? 0 : 1);
}

static qaws_status qf_network(qaws_hm const* dm, mf_level const* L, int const* T, unsigned char const* K, double const* Tc, unsigned char const* fr,
	unsigned char const* sing, unsigned char const* g, int parity, qf_net* net)
{
	unsigned int nf = dm->nf, Z = 2 * nf, c, ne = 0, j, open = 0;
	qaws_status st;
	memset(net, 0, sizeof *net);
	for (c = 0; c < dm->nc; c++)
	{
		ne += dm->cf[c] != HM_NONE && (dm->co[c] == HM_NONE || c < dm->co[c]);
		open |= dm->cf[c] != HM_NONE && dm->co[c] == HM_NONE;
	}
	if (ne > 0xFFFFFFFFu / 16 || nf > 0xFFFFFFFFu / 32)   /* 8 ne + 10 nf + 1 arcs count in unsigned int */
		return QAWS_STATUS_ALLOCATION_FAILURE;
	st = qaws_flow_create(&net->flow, 8 * ne + 10 * nf + 1, 2 * nf + 1);
	net->rec = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (8 * (size_t)ne + 10 * (size_t)nf + 1)));
	if (st != QAWS_STATUS_OK || !net->rec)
	{
		qf_net_free(net);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (c = 0; c < dm->nc; c++)
	{
		unsigned int o = dm->co[c], F = dm->cf[c], G;
		int lined, ax;
		if (F == HM_NONE || (o != HM_NONE && o < c))
			continue;
		G = o != HM_NONE ? dm->cf[o] : HM_NONE;
		if (open && ((L->con[dm->cv[c]] | L->con[hm_to(dm, c)] | L->con[dm->cv[hm_prev(dm, c)]] | (o != HM_NONE ? L->con[dm->cv[hm_prev(dm, o)]] : 0))
			& MF_CONE))
			continue;   /* (an edge of a placed cone's star kept, the cone on its vertex's lattice point: the border takes what
			               cannot close; on a closed surface the cones must, at their price below) */
		/* feature and border edges are lattice lines: steps along themselves only */
		lined = o == HM_NONE || ((dm->cflag[c] | dm->cflag[o]) & HM_FEATURE);
		ax = lined ? qf_axis(dm, L, c) : 0;
		if (!lined && (L->con[dm->cv[c]] & L->con[hm_to(dm, c)] & MF_GUIDE))
		{
			/* between two vertices of a guide line (either side of it): steps along the line only */
			double tq[3];
			double const* cq = L->cq + 3 * (size_t)dm->cv[c];
			qf_cross(L->nrm + 3 * (size_t)dm->cv[c], L->q + 3 * (size_t)dm->cv[c], tq);
			lined = 1;
			ax = fabs(qf_dot(cq, L->q + 3 * (size_t)dm->cv[c])) >= fabs(qf_dot(cq, tq)) ? 0 : 1;
		}
		for (j = 0; j < 2; j++)
		{
			unsigned int nF, nG, P, Q, a, k;
			int sF, sG = 0, t = T[2 * (size_t)c + j];
			double fa = Tc[4 * (size_t)c + j], fb = Tc[4 * (size_t)c + 2 + j];
			if (lined && (int)j != ax)
				continue;
			nF = qf_node(F, sing, Z, g[F] + fr[c], (int)j, 1, &sF);
			nG = G == HM_NONE ? Z : qf_node(G, sing, Z, g[G] + fr[o] + K[c], (int)j, -1, &sG);
			if (nF == nG)
				continue;   /* (both the sink) */
			if (sF && sG && sF == sG)
			{
				net->cut++;
				continue;
			}
			/* raising x by a unit adds +1 to P's row and -1 to Q's: a unit of flow P -> Q */
			if (sF > 0 || (!sF && sG < 0))
			{
				P = nF;
				Q = nG;
			}
			else
			{
				P = nG;
				Q = nF;
			}
			for (k = 0; k < 4; k++)
			{
				int s = k < 2 ? 1 : -1;
				long long w = (k & 1) ? qf_cost(t + s, fa + s, fb - s, s, 1) : qf_cost(t, fa, fb, s, 0);
				a = s > 0 ? qaws_flow_arc(net->flow, P, Q, 1, w) : qaws_flow_arc(net->flow, Q, P, 1, w);
				if (a != QAWS_FLOW_NONE)
				{
					net->rec[a] = c * 4 + j * 2 + (s < 0);
					net->na = a + 1;
				}
			}
		}
	}
	/* closing, a cone takes units either way at a price (it moves off its face): a few edges for the first, three times
	   that for the second, ten times for any more */
	for (j = 0; !parity && j < nf; j++)
		for (c = 0; sing[j] && c < 6; c++)
		{
			/* (a placed cone, MF_CONE, a hundred times dearer: it moves only where the lattice cannot close otherwise) */
			long long w = QF_CONE * (c < 2 ? 1 : c < 4 ? 3 : 10) * (((L->con[dm->cv[dm->fc[j]]] | L->con[dm->cv[dm->fc[j] + 1]]
				| L->con[dm->cv[dm->fc[j] + 2]]) & MF_CONE) ? 100 : 1);
			unsigned int a = qaws_flow_arc(net->flow, c & 1 ? Z : 2 * j, c & 1 ? 2 * j : Z, c < 4 ? 1 : (int)nf, w);
			if (a != QAWS_FLOW_NONE)
			{
				net->rec[a] = HM_NONE;
				net->na = a + 1;
			}
		}
	/* a graded lattice needs extra rows (dislocations) where its spacing changes: a face there may keep one, at a price (a
	   few edges' changes: pairs of slips the field left by noise still close) */
	for (j = 0; !parity && L->sc && j < nf; j++)
	{
		double lo = HUGE_VAL, hi = 0;
		if (sing[j])
			continue;
		for (c = dm->fc[j]; c < dm->fc[j] + 3; c++)
		{
			lo = fmin(lo, L->sc[dm->cv[c]]);
			hi = fmax(hi, L->sc[dm->cv[c]]);
		}
		for (c = 0; hi > (1 + QF_GRADED) * lo && c < 4; c++)
		{
			unsigned int a = qaws_flow_arc(net->flow, c & 1 ? Z : 2 * j + (c >> 1), c & 1 ? 2 * j + (c >> 1) : Z, 1, QF_SLIP);
			if (a != QAWS_FLOW_NONE)
			{
				net->rec[a] = HM_NONE;
				net->na = a + 1;
			}
		}
	}
	return QAWS_STATUS_OK;
}

/* the flow written back onto the steps (both corners of each edge); returns the units changed */
static unsigned int qf_apply(qaws_hm const* dm, int* T, unsigned char const* K, qf_net const* net)
{
	unsigned int a, c, n = 0;
	for (a = 0; a < net->na; a++)
	{
		int fl = qaws_flow_get(net->flow, a);
		unsigned int r = net->rec[a];
		if (!fl || r == HM_NONE)
			continue;
		T[2 * (size_t)(r >> 2) + ((r >> 1) & 1)] += (r & 1) ? -fl : fl;
		n += (unsigned int)fl;
	}
	for (c = 0; c < dm->nc; c++)
	{
		unsigned int o = dm->co[c];
		int m[2];
		if (dm->cf[c] == HM_NONE || o == HM_NONE || o < c)
			continue;
		m[0] = -T[2 * (size_t)c];
		m[1] = -T[2 * (size_t)c + 1];
		qf_turn(m, K[c], T + 2 * (size_t)o);
	}
	return n;
}

/* a cone where the frames turn by a quarter is a lattice point when its sum's components add up to an even number
   (det(I - R) = 2) */
static int qf_odd(int turn, int const* S)
{
	return (turn & 1) && ((S[0] + S[1]) & 1);
}

int qf_image(qaws_hm const* dm, unsigned int f, int const* T, unsigned char const* K, int* e)
{
	unsigned int c0 = dm->fc[f], i, j;
	for (i = 0; i < 3; i++)
	{
		int turn = 0, S[2] = { 0, 0 };
		for (j = 0; j < 3; j++)
		{
			unsigned int c = c0 + (i + j) % 3;
			qf_turn(T + 2 * (size_t)c, turn, e + 2 * j);
			S[0] += e[2 * j];
			S[1] += e[2 * j + 1];
			turn -= K[c];
		}
		if (!S[0] && !S[1])
			return (int)i;
		if (!(turn & 3))
			break;   /* (the field does not turn: closed from every corner or from none) */
	}
	return -1;
}

/* the lattice area (twice) of a face whose steps close, 0 for one that does not */
static int qf_area(qaws_hm const* dm, unsigned int f, int const* T, unsigned char const* K)
{
	int e[6];
	return qf_image(dm, f, T, K, e) < 0 ? 0 : e[0] * e[3] - e[1] * e[2];
}

/* vertex v's lattice point moved by d (in v's frame): its steps out less d, in more d (every face round v keeps its sum) */
static void qf_move(qaws_hm const* dm, int* T, unsigned char const* K, unsigned int v, int const* d)
{
	unsigned int c, last = HM_NONE;
	for (c = dm->vc[v]; c != HM_NONE; c = hm_round(dm, v, c))
	{
		unsigned int o = dm->co[c];
		T[2 * (size_t)c] -= d[0];
		T[2 * (size_t)c + 1] -= d[1];
		if (o != HM_NONE)
		{
			int m[2];
			m[0] = -T[2 * (size_t)c];
			m[1] = -T[2 * (size_t)c + 1];
			qf_turn(m, K[c], T + 2 * (size_t)o);
		}
		last = c;
	}
	if (last != HM_NONE && dm->co[hm_prev(dm, last)] == HM_NONE)
	{
		/* (the border edge coming in, stored on its far end only) */
		unsigned int p = hm_prev(dm, last);
		int w[2];
		qf_turn(d, -(int)K[p], w);
		T[2 * (size_t)p] += w[0];
		T[2 * (size_t)p + 1] += w[1];
	}
}

/* round v: its faces turned over, its steps of two, its cones off their vertices (out may be NULL); returns the three
   weighed, a cone moved off its vertex dearest */
static int qf_bad(qaws_hm const* dm, unsigned int v, int const* T, unsigned char const* K, unsigned char* fr, int* out)
{
	unsigned int c, last = HM_NONE;
	int n[3] = { 0, 0, 0 };
	for (c = dm->vc[v]; c != HM_NONE; c = hm_round(dm, v, c))
	{
		int e[6], S[2];
		if (qf_image(dm, dm->cf[c], T, K, e) < 0)
			n[2] += qf_closure(dm, dm->cf[c], T, K, fr, S) != 0;
		else
			n[0] += e[0] * e[3] - e[1] * e[2] < 0;
		n[1] += abs(T[2 * (size_t)c]) > 1 || abs(T[2 * (size_t)c + 1]) > 1;
		last = c;
	}
	if (last != HM_NONE && dm->co[hm_prev(dm, last)] == HM_NONE)
	{
		unsigned int p = hm_prev(dm, last);
		n[1] += abs(T[2 * (size_t)p]) > 1 || abs(T[2 * (size_t)p + 1]) > 1;
	}
	if (out)
		memcpy(out, n, sizeof n);
	return n[0] + 2 * n[1] + 8 * n[2];
}

/* the vertices of v's lattice point (joined to it by zero steps; at most cap, none pinned or on a line, else 0) and the
   quarter turns from v's frame to each one's */
static unsigned int qf_cluster(qaws_hm const* dm, mf_level const* L, int const* T, unsigned char const* K, unsigned int v, unsigned char* mark,
	unsigned int* list, unsigned char* rot, unsigned int cap)
{
	unsigned int n = 1, i, c, ok = 1;
	list[0] = v;
	rot[0] = 0;
	mark[v] = 1;
	for (i = 0; i < n && ok; i++)
	{
		unsigned int u = list[i];
		if (L->con[u] & (MF_PINNED | MF_ON_LINE))
			ok = 0;
		for (c = dm->vc[u]; c != HM_NONE && ok; c = hm_round(dm, u, c))
		{
			unsigned int w = hm_to(dm, c);
			if (mark[w] || T[2 * (size_t)c] || T[2 * (size_t)c + 1])
				continue;
			if (n == cap)
				ok = 0;
			else
			{
				mark[w] = 1;
				list[n] = w;
				rot[n++] = (unsigned char)((rot[i] + K[c]) & 3);
			}
		}
	}
	for (i = 0; i < n; i++)
		mark[list[i]] = 0;
	return ok ? n : 0;
}

/* The cone of a singular face whose sum is even (a lattice point) put on one of its vertices: that vertex's lattice
   point moved onto the cone, (I - R^k)^-1 S away (its sum then 0), alone or with the rest of its lattice point (the
   vertices joined to it by zero steps), whichever adds the least harm and no more than two steps of two and two faces
   turned over (flip removal mends those; a cone left off its vertices is a 3-5 pair for sure: with one step of two
   allowed, test 134's six blobs pure at 0.11 left 6 cones off, with two 2); returns the cones placed.
   shortcut: a cone no vertex can reach that cheaply stays off its vertices (its hole is filled by a fan, a 3-5 pair
   beside it); route units to change its face's sum instead when those pairs matter */
static unsigned int qf_cones(qaws_hm const* dm, mf_level const* L, int* T, unsigned char const* K, unsigned char* fr)
{
	unsigned int f, i, j, q, way, m, bm, placed = 0, list[64], bl[64];
	unsigned char rot[64], brot[64], *mark = (unsigned char*)qaws_internal_zalloc((size_t)dm->nv + 1, 1);
	for (f = 0; f < dm->nf; f++)
	{
		int S[2], t = qf_closure(dm, f, T, K, fr, S), bb = 0, bd[2] = { 0, 0 }, e[6];
		if (!t || qf_image(dm, f, T, K, e) >= 0)
			continue;   /* (the field does not turn, or the cone is on a vertex already) */
		for (i = 0, bm = 0; i < 3; i++)
		{
			unsigned int v = dm->cv[dm->fc[f] + i];
			int s[2] = { 0, 0 }, d[2], turn = 0;
			/* the sum from corner i, in its vertex's frame */
			for (j = 0; j < 3; j++)
			{
				unsigned int c = dm->fc[f] + (i + j) % 3;
				int w[2];
				qf_turn(T + 2 * (size_t)c, turn, w);
				s[0] += w[0];
				s[1] += w[1];
				turn -= K[c];
			}
			if (t == 2 ? ((s[0] | s[1]) & 1) : ((s[0] + s[1]) & 1))
				break;   /* (odd: no lattice point) */
			d[0] = t == 2 ? s[0] / 2 : t == 1 ? (s[0] - s[1]) / 2 : (s[0] + s[1]) / 2;
			d[1] = t == 2 ? s[1] / 2 : t == 1 ? (s[0] + s[1]) / 2 : (s[1] - s[0]) / 2;
			for (way = 0; way < 2; way++)
			{
				int before[3] = { 0, 0, 0 }, after[3] = { 0, 0, 0 }, n3[3], b0 = 0, b1 = 0, dq[2], closed;
				if (way == 0)
				{
					m = 1;
					list[0] = v;
					rot[0] = 0;
					if (L->con[v] & MF_PINNED)
						continue;
					if (L->con[v] & MF_ON_LINE)
					{
						double tq[3];
						double const* cq = L->cq + 3 * (size_t)v;
						qf_cross(L->nrm + 3 * (size_t)v, L->q + 3 * (size_t)v, tq);
						if (d[fabs(qf_dot(cq, L->q + 3 * (size_t)v)) >= fabs(qf_dot(cq, tq)) ? 1 : 0])
							continue;   /* (off its line) */
					}
				}
				else if (!mark || (m = qf_cluster(dm, L, T, K, v, mark, list, rot, 64)) < 2)
					break;
				for (q = 0; q < m; q++)
				{
					b0 += qf_bad(dm, list[q], T, K, fr, n3);
					for (j = 0; j < 3; j++)
						before[j] += n3[j];
				}
				for (q = 0; q < m; q++)
				{
					qf_turn(d, rot[q], dq);
					qf_move(dm, T, K, list[q], dq);
				}
				for (q = 0; q < m; q++)
				{
					b1 += qf_bad(dm, list[q], T, K, fr, n3);
					for (j = 0; j < 3; j++)
						after[j] += n3[j];
				}
				closed = qf_image(dm, f, T, K, e) >= 0;
				for (q = 0; q < m; q++)
				{
					qf_turn(d, rot[q] + 2, dq);
					qf_move(dm, T, K, list[q], dq);
				}
				if (!closed || after[1] > before[1] + 2 || after[0] > before[0] + 2)
					continue;
				if (!bm || b1 - b0 < bb)
				{
					bb = b1 - b0;
					bm = m;
					bd[0] = d[0];
					bd[1] = d[1];
					memcpy(bl, list, sizeof(unsigned int) * m);
					memcpy(brot, rot, m);
				}
			}
		}
		if (bm && i == 3)
		{
			for (q = 0; q < bm; q++)
			{
				int dq[2];
				qf_turn(bd, brot[q], dq);
				qf_move(dm, T, K, bl[q], dq);
			}
			placed++;
		}
	}
	qaws_internal_free(mark);
	return placed;
}

/* QuadriFlow's flip removal: the vertices of faces turned over in the lattice (or with a step of two) moved greedily by
   the step that leaves round them the fewest turned over and steps of two (along their line for those on one, pinned
   ones kept); else with a neighbour, else with the rest of their lattice point. Returns the faces still turned over. */
static unsigned int qf_unflip(qaws_hm const* dm, mf_level const* L, int* T, unsigned char const* K, unsigned char* fr, unsigned char* flip)
{
	static int const step[8][2] = { { 1, 0 }, { -1, 0 }, { 0, 1 }, { 0, -1 }, { 1, 1 }, { -1, -1 }, { 1, -1 }, { -1, 1 } };
	unsigned int f, k, pass, n = 0, moved = 1, list[64], m;
	unsigned char rot[64], *mark = (unsigned char*)qaws_internal_zalloc((size_t)dm->nv + 1, 1);   /* (none: no moves of whole lattice points) */
	for (pass = 0; pass < 16 && moved; pass++)
	{
		moved = 0;
		for (f = 0; f < dm->nf; f++)
		{
			unsigned int c = dm->fc[f];
			flip[f] = qf_area(dm, f, T, K) < 0;
			for (k = 0; k < 3; k++)
				flip[f] |= abs(T[2 * (size_t)(c + k)]) > 1 || abs(T[2 * (size_t)(c + k) + 1]) > 1;
		}
		for (f = 0; f < dm->nf; f++)
			for (k = 0; k < 3 && flip[f]; k++)
			{
				unsigned int v = dm->cv[dm->fc[f] + k], j, best = 8, ax = 2;
				int b0, bb;
				if (L->con[v] & MF_PINNED)
					continue;
				if (L->con[v] & MF_ON_LINE)
				{
					double t[3];
					double const* cq = L->cq + 3 * (size_t)v;
					qf_cross(L->nrm + 3 * (size_t)v, L->q + 3 * (size_t)v, t);
					ax = fabs(qf_dot(cq, L->q + 3 * (size_t)v)) >= fabs(qf_dot(cq, t)) ? 0 : 1;
				}
				b0 = bb = qf_bad(dm, v, T, K, fr, NULL);
				for (j = 0; j < 8; j++)
				{
					int back[2], b;
					if (ax < 2 && (!step[j][ax] || step[j][ax ^ 1]))
						continue;
					qf_move(dm, T, K, v, step[j]);
					b = qf_bad(dm, v, T, K, fr, NULL);
					back[0] = -step[j][0];
					back[1] = -step[j][1];
					qf_move(dm, T, K, v, back);
					if (b < bb)
					{
						bb = b;
						best = j;
					}
				}
				if (best < 8 && bb < b0)
				{
					qf_move(dm, T, K, v, step[best]);
					moved++;
				}
				else if (ax == 2)
				{
					/* no move of v alone helps (a fold wider than a face): v and a neighbour moved together */
					unsigned int c, bw = HM_NONE, bc = 0;
					for (c = dm->vc[v]; c != HM_NONE; c = hm_round(dm, v, c))
					{
						unsigned int w = hm_to(dm, c);
						int p0;
						if (L->con[w] & (MF_PINNED | MF_ON_LINE))
							continue;
						p0 = b0 + qf_bad(dm, w, T, K, fr, NULL);
						for (j = 0; j < 8; j++)
						{
							int dw[2], back[2], b;
							qf_turn(step[j], K[c], dw);
							qf_move(dm, T, K, v, step[j]);
							qf_move(dm, T, K, w, dw);
							b = qf_bad(dm, v, T, K, fr, NULL) + qf_bad(dm, w, T, K, fr, NULL);
							back[0] = -dw[0];
							back[1] = -dw[1];
							qf_move(dm, T, K, w, back);
							back[0] = -step[j][0];
							back[1] = -step[j][1];
							qf_move(dm, T, K, v, back);
							if (b < p0 && (bw == HM_NONE || b < bb))
							{
								bb = b;
								bw = w;
								bc = c;
								best = j;
							}
						}
					}
					if (bw != HM_NONE)
					{
						int dw[2];
						qf_turn(step[best], K[bc], dw);
						qf_move(dm, T, K, v, step[best]);
						qf_move(dm, T, K, bw, dw);
						moved++;
					}
					else if (mark && (m = qf_cluster(dm, L, T, K, v, mark, list, rot, 64)) > 1)
					{
						/* nor with one neighbour: v's whole lattice point (the vertices joined to it by zero steps) moved */
						int c0b = 0, b, d[2];
						unsigned int q;
						for (q = 0; q < m; q++)
							c0b += qf_bad(dm, list[q], T, K, fr, NULL);
						for (j = 0, best = 8, bb = c0b; j < 8; j++)
						{
							for (q = 0; q < m; q++)
							{
								qf_turn(step[j], rot[q], d);
								qf_move(dm, T, K, list[q], d);
							}
							for (q = 0, b = 0; q < m; q++)
								b += qf_bad(dm, list[q], T, K, fr, NULL);
							for (q = 0; q < m; q++)
							{
								qf_turn(step[j], rot[q] + 2, d);
								qf_move(dm, T, K, list[q], d);
							}
							if (b < bb)
							{
								bb = b;
								best = j;
							}
						}
						for (q = 0; best < 8 && q < m; q++)
						{
							qf_turn(step[best], rot[q], d);
							qf_move(dm, T, K, list[q], d);
						}
						moved += best < 8;
					}
				}
			}
	}
	qaws_internal_free(mark);
	for (f = 0; f < dm->nf; f++)
		n += qf_area(dm, f, T, K) < 0;
	return n;
}

unsigned int qf_consistent(qaws_hm const* dm, mf_level const* L, int* T, unsigned char const* K, double const* Tc, unsigned int* out_odd,
	qf_report* out_report)
{
	unsigned int nf = dm->nf, f, i, left = HM_NONE, odd = 0, head, tail;
	unsigned char *fr = (unsigned char*)qaws_internal_scratch_malloc((size_t)dm->nc + 1), *sing = (unsigned char*)qaws_internal_malloc((size_t)nf + 1);
	unsigned char *g = (unsigned char*)qaws_internal_scratch_malloc((size_t)nf + 1), *seen = (unsigned char*)qaws_internal_zalloc((size_t)nf + 1, 1);
	unsigned int *queue = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (nf + 2))), *targets = NULL;
	int* S = (int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(int), (2 * (size_t)nf + 2)));
	qf_report rep;
	qf_net net;
	memset(&rep, 0, sizeof rep);
	memset(&net, 0, sizeof net);
	if (out_odd)
		*out_odd = 0;
	if (!fr || !sing || !g || !seen || !queue || !S)
		goto done;
	for (f = 0; f < nf; f++)
		sing[f] = (unsigned char)qf_closure(dm, f, T, K, fr, S + 2 * (size_t)f);
	/* the frames combed breadth first across edges between faces where the field does not turn */
	for (f = 0; f < nf; f++)
	{
		if (sing[f] || seen[f])
			continue;
		g[f] = 0;
		seen[f] = 1;
		queue[0] = f;
		for (head = 0, tail = 1; head < tail; head++)
		{
			unsigned int F = queue[head], c;
			for (c = dm->fc[F]; c < dm->fc[F] + 3; c++)
			{
				unsigned int o = dm->co[c], G;
				if (o == HM_NONE)
					continue;
				G = dm->cf[o];
				if (sing[G] || seen[G])
					continue;
				g[G] = (unsigned char)((g[F] + fr[c] - fr[o] - K[c]) & 3);
				seen[G] = 1;
				queue[tail++] = G;
			}
		}
	}
	for (f = 0; f < nf; f++)
		if (sing[f])
			g[f] = 0;
	/* 1. every face where the field does not turn closed: supplies -E, the sink the rest */
	if (qf_network(dm, L, T, K, Tc, fr, sing, g, 0, &net) != QAWS_STATUS_OK)
		goto done;
	{
		int sum = 0, un = 0;
		for (f = 0; f < nf; f++)
		{
			int E[2];
			if (sing[f])
				continue;
			qf_turn(S + 2 * (size_t)f, g[f], E);
			qaws_flow_supply(net.flow, 2 * f, -E[0]);
			qaws_flow_supply(net.flow, 2 * f + 1, -E[1]);
			sum += E[0] + E[1];
			rep.supply += (unsigned int)(abs(E[0]) + abs(E[1]));
		}
		qaws_flow_supply(net.flow, 2 * nf, sum);
		qaws_flow_solve(&un, &rep.cost, net.flow);
		rep.unrouted = (unsigned int)un;
	}
	rep.cut = net.cut;
	rep.changed = qf_apply(dm, T, K, &net);
	qf_net_free(&net);
	/* 2. odd cones paired by paths of unit changes (the faces between keep their sums), or sent to the border.
	   shortcut: in index order, each to the nearest odd one left; a min-cost matching when pairs should be cheapest
	   together */
	targets = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (nf + 2)));
	if (!targets)
		goto done;
	for (f = 0, odd = 0; f < nf; f++)
	{
		qf_closure(dm, f, T, K, fr, S + 2 * (size_t)f);
		seen[f] = (unsigned char)qf_odd(sing[f], S + 2 * (size_t)f);
		odd += seen[f];
	}
	if (odd)
	{
		if (qf_network(dm, L, T, K, Tc, fr, sing, g, 1, &net) != QAWS_STATUS_OK)
			goto done;
		for (f = 0; f < nf; f++)
		{
			unsigned int nt = 0, reached, k;
			if (!seen[f])
				continue;
			for (k = 0; k < nf; k++)
				if (seen[k] && k != f)
					targets[nt++] = 2 * k;
			targets[nt++] = 2 * nf;
			qaws_flow_augment(&reached, NULL, net.flow, targets, nt, 2 * f);
			if (reached == QAWS_FLOW_NONE)
				continue;
			rep.paths++;
			seen[f] = 0;
			if (reached != 2 * nf)
				seen[reached / 2] = 0;
		}
		rep.changed += qf_apply(dm, T, K, &net);
		qf_net_free(&net);
	}
	/* 3. the cones on vertices, the lattice untangled where the changes turned faces over */
	rep.cones = qf_cones(dm, L, T, K, fr);
	for (f = 0; f < nf; f++)
		rep.flips += qf_area(dm, f, T, K) < 0;
	rep.flips_left = qf_unflip(dm, L, T, K, fr, seen);
	/* what is left */
	for (f = 0, left = 0, odd = 0; f < nf; f++)
	{
		int* s = S + 2 * (size_t)f;
		qf_closure(dm, f, T, K, fr, s);
		left += !sing[f] && (s[0] || s[1]);
		odd += (unsigned int)qf_odd(sing[f], s);
	}
	for (i = 0; i < dm->nc; i++)
		if (dm->cf[i] != HM_NONE && (dm->co[i] == HM_NONE || i < dm->co[i]))
			rep.twos += abs(T[2 * (size_t)i]) >= 2 || abs(T[2 * (size_t)i + 1]) >= 2;
	if (out_odd)
		*out_odd = odd;
done:
	qf_net_free(&net);
	qaws_internal_free(fr); qaws_internal_free(sing); qaws_internal_free(g); qaws_internal_free(seen); qaws_internal_free(queue); qaws_internal_free(targets); qaws_internal_free(S);
	if (out_report)
		*out_report = rep;
	return left;
}

void qf_refit(qaws_hm const* dm, mf_level* L, double rho, int const* T, unsigned char const* K, unsigned int sweeps)
{
	unsigned int it, a, c, k;
	for (it = 0; it < sweeps; it++)
		for (a = 0; a < dm->nv; a++)
		{
			double const *xa = L->x + 3 * (size_t)a, *na = L->nrm + 3 * (size_t)a, *qa = L->q + 3 * (size_t)a;
			double sum[3] = { 0, 0, 0 }, ta[3], d;
			double* o = L->o + 3 * (size_t)a;
			unsigned int n = 0, last = HM_NONE;
			if (dm->vc[a] == HM_NONE)
				continue;
			if (L->con[a] & MF_PINNED)
			{
				memcpy(o, xa, sizeof(double) * 3);
				continue;
			}
			qf_cross(na, qa, ta);
			/* each neighbour b: its o unfolded into a's plane, less the step a -> b */
			for (c = dm->vc[a]; c != HM_NONE; c = hm_round(dm, a, c))
			{
				unsigned int b = hm_to(dm, c);
				int const* t = T + 2 * (size_t)c;
				double ob[3];
				mf_unfold(xa, na, L->x + 3 * (size_t)b, L->nrm + 3 * (size_t)b, L->o + 3 * (size_t)b, ob);
				for (k = 0; k < 3; k++)
					sum[k] += ob[k] - mf_pair_scale(L, a, b, rho) * (t[0] * qa[k] + t[1] * ta[k]);
				n++;
				last = c;
			}
			if (last != HM_NONE)
			{
				/* on the border, the neighbour across its other border edge */
				unsigned int p = hm_prev(dm, last), b = dm->cv[p];
				if (dm->co[p] == HM_NONE)
				{
					int m[2], t[2];
					double ob[3];
					m[0] = -T[2 * (size_t)p];
					m[1] = -T[2 * (size_t)p + 1];
					qf_turn(m, K[p], t);
					mf_unfold(xa, na, L->x + 3 * (size_t)b, L->nrm + 3 * (size_t)b, L->o + 3 * (size_t)b, ob);
					for (k = 0; k < 3; k++)
						sum[k] += ob[k] - mf_pair_scale(L, a, b, rho) * (t[0] * qa[k] + t[1] * ta[k]);
					n++;
				}
			}
			if (!n)
				continue;
			for (k = 0; k < 3; k++)
				o[k] = sum[k] / n;
			d = (o[0] - xa[0]) * na[0] + (o[1] - xa[1]) * na[1] + (o[2] - xa[2]) * na[2];
			for (k = 0; k < 3; k++)
				o[k] -= d * na[k];
			if (L->con[a] & MF_ON_LINE)
			{
				double const* cq = L->cq + 3 * (size_t)a;
				d = (o[0] - xa[0]) * cq[0] + (o[1] - xa[1]) * cq[1] + (o[2] - xa[2]) * cq[2];
				for (k = 0; k < 3; k++)
					o[k] = xa[k] + d * cq[k];
			}
		}
}
