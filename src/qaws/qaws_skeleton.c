/*
 * Straight skeletons (see qaws_skeleton.h).
 *
 *   1. region: the paths united under the fill rule (qaws_clip): outer loops
 *      counter-clockwise, holes clockwise, the interior on the left
 *   2. planes: every straight edge, n . p = d + w t with n its inward normal
 *   3. corners: the wavefront's vertices in circular lists, each between
 *      the planes of its edges, moving at a constant velocity
 *   4. events: for every new corner, the times it meets its neighbours
 *      (three consecutive planes) and, when reflex, every plane (split);
 *      a heap by time
 *   5. the generalised intersection event at each popped point: the corners
 *      there and the edges through it make chains; chain inner edges end;
 *      chains sorted by angle are joined; corners between antiparallel
 *      edges collapse along them (spikes)
 *   6. faces: each plane's arcs walked from its edge's end back to its start
 */

#include "qaws_skeleton.h"
#include "qaws_curve.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "internal/qaws_internal_types.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define SK_NONE 0xFFFFFFFFu
#if QAWS_SCALAR_IS_FLOAT
#define SK_TOL_REL 1e-6          /* default tolerance, of the extent: the input's float rounding */
#else
#define SK_TOL_REL 1e-9          /* default tolerance, of the extent */
#endif

typedef struct sk_plane
{
	double n[2], d, w;          /* n . p = d + w t, n the inward unit normal */
	unsigned int na, nb;        /* the input edge's nodes, start and end */
	unsigned int path, curve, span;
} sk_plane;

typedef struct sk_corner
{
	double p[2], t0, v[2];      /* the position at time t0, the velocity */
	unsigned int in, out;       /* the planes of the edges before and after */
	unsigned int prev, next;
	unsigned int node;          /* where it started */
	int alive;
} sk_corner;

typedef struct sk_event
{
	double t, p[2];
	unsigned int corner;
	unsigned int seq;
} sk_event;

struct qaws_skeleton
{
	sk_plane* planes;
	unsigned int nplanes;
	unsigned int* loop_first;   /* loop l: planes loop_first[l] .. loop_first[l + 1], in order */
	unsigned int nloops;
	qaws_skeleton_node* nodes;
	unsigned int nnodes, capnodes;
	qaws_skeleton_arc* arcs;
	unsigned int narcs, caparcs;
	unsigned int* face_first;   /* face f: face_nodes[face_first[f] .. face_first[f + 1]] */
	unsigned int* face_nodes;
	double tol, tmax;
};

/* one simulation of the wavefront, to the end (recording the skeleton) or
   to a time */
typedef struct sk_run
{
	qaws_skeleton* s;
	int record;
	double stop, now, eps, eps_t;
	sk_corner* c;
	unsigned int nc, capc;
	unsigned int* refs;         /* per plane: alive corners whose out edge is on it */
	sk_event* heap;
	unsigned int nh, caph, seq;
	unsigned int* tag;          /* per corner: the stamp of the event it is part of */
	unsigned int captag;
	unsigned int stamp;
	unsigned int node_from;     /* nodes of the current time start here */
	qaws_status st;
} sk_run;

static int sk_grow(void** p, unsigned int* cap, unsigned int need, size_t size)
{
	unsigned int n;
	void* q;
	if (need <= *cap)
		return 1;
	n = *cap ? *cap : 16;
	while (n < need) n *= 2;
	q = realloc(*p, size * n);
	if (!q)
		return 0;
	*p = q;
	*cap = n;
	return 1;
}

/* the edge's direction: the interior on its left */
static void sk_dir(sk_plane const* pl, double* d)
{
	d[0] = pl->n[1];
	d[1] = -pl->n[0];
}

static void sk_pos(sk_corner const* c, double t, double* p)
{
	p[0] = c->p[0] + c->v[0] * (t - c->t0);
	p[1] = c->p[1] + c->v[1] * (t - c->t0);
}

/* three planes meeting: n_i . p - w_i t = d_i */
static int sk_solve3(sk_plane const* a, sk_plane const* b, sk_plane const* c, double* p, double* t)
{
	double m[3][4] = {
		{ a->n[0], a->n[1], -a->w, a->d },
		{ b->n[0], b->n[1], -b->w, b->d },
		{ c->n[0], c->n[1], -c->w, c->d } };
	double det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
		+ m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
	double dx, dy, dt;
	if (!(fabs(det) > 1e-14))
		return 0;
	dx = m[0][3] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][3] * m[2][2] - m[1][2] * m[2][3])
		+ m[0][2] * (m[1][3] * m[2][1] - m[1][1] * m[2][3]);
	dy = m[0][0] * (m[1][3] * m[2][2] - m[1][2] * m[2][3]) - m[0][3] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
		+ m[0][2] * (m[1][0] * m[2][3] - m[1][3] * m[2][0]);
	dt = m[0][0] * (m[1][1] * m[2][3] - m[1][3] * m[2][1]) - m[0][1] * (m[1][0] * m[2][3] - m[1][3] * m[2][0])
		+ m[0][3] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
	p[0] = dx / det;
	p[1] = dy / det;
	*t = dt / det;
	return 1;
}

/* the velocity of a corner between planes a and b: 0, or 1 when a and b
   are antiparallel (a spike) */
static int sk_velocity(sk_plane const* a, sk_plane const* b, double* v)
{
	double det = a->n[0] * b->n[1] - a->n[1] * b->n[0], dot = a->n[0] * b->n[0] + a->n[1] * b->n[1];
	if (fabs(det) <= 1e-12)
	{
		if (dot < 0)
		{
			v[0] = v[1] = 0;
			return 1;
		}
		/* collinear: the average speed (exact for equal weights) */
		v[0] = 0.5 * (a->w + b->w) * a->n[0];
		v[1] = 0.5 * (a->w + b->w) * a->n[1];
		return 0;
	}
	v[0] = (a->w * b->n[1] - a->n[1] * b->w) / det;
	v[1] = (a->n[0] * b->w - a->w * b->n[0]) / det;
	return 0;
}

/* ------------------------------------------------------------------------ */
/*  Nodes, arcs, the heap                                                   */
/* ------------------------------------------------------------------------ */

static unsigned int sk_node(sk_run* r, double const* p, double t)
{
	qaws_skeleton* s = r->s;
	unsigned int i;
	if (!r->record)
		return 0;
	for (i = s->nnodes; i-- > r->node_from;)
	{
		qaws_skeleton_node const* n = &s->nodes[i];
		if (fabs((double)n->time - t) <= r->eps_t && hypot((double)n->position.x - p[0], (double)n->position.y - p[1]) <= r->eps)
			return i;
	}
	if (!sk_grow((void**)&s->nodes, &s->capnodes, s->nnodes + 1, sizeof(qaws_skeleton_node)))
	{
		r->st = QAWS_STATUS_ALLOCATION_FAILURE;
		return 0;
	}
	s->nodes[s->nnodes].position.x = (qaws_scalar)p[0];
	s->nodes[s->nnodes].position.y = (qaws_scalar)p[1];
	s->nodes[s->nnodes].time = (qaws_scalar)t;
	if (t > s->tmax) s->tmax = t;
	return s->nnodes++;
}

static void sk_arc(sk_run* r, unsigned int a, unsigned int b, unsigned int fa, unsigned int fb)
{
	qaws_skeleton* s = r->s;
	if (!r->record || a == b)
		return;
	if (!sk_grow((void**)&s->arcs, &s->caparcs, s->narcs + 1, sizeof(qaws_skeleton_arc)))
	{
		r->st = QAWS_STATUS_ALLOCATION_FAILURE;
		return;
	}
	s->arcs[s->narcs].node_a = a;
	s->arcs[s->narcs].node_b = b;
	s->arcs[s->narcs].face_a = fa;
	s->arcs[s->narcs].face_b = fb;
	s->narcs++;
}

static int sk_before(sk_event const* a, sk_event const* b)
{
	return a->t < b->t || (a->t == b->t && a->seq < b->seq);
}

static void sk_push(sk_run* r, double t, double const* p, unsigned int corner)
{
	unsigned int i;
	sk_event e;
	if (!sk_grow((void**)&r->heap, &r->caph, r->nh + 1, sizeof(sk_event)))
	{
		r->st = QAWS_STATUS_ALLOCATION_FAILURE;
		return;
	}
	e.t = t; e.p[0] = p[0]; e.p[1] = p[1]; e.corner = corner; e.seq = r->seq++;
	i = r->nh++;
	while (i > 0 && sk_before(&e, &r->heap[(i - 1) / 2]))
	{
		r->heap[i] = r->heap[(i - 1) / 2];
		i = (i - 1) / 2;
	}
	r->heap[i] = e;
}

static sk_event sk_pop(sk_run* r)
{
	sk_event top = r->heap[0], last = r->heap[--r->nh];
	unsigned int i = 0;
	for (;;)
	{
		unsigned int a = 2 * i + 1, b = a + 1, m = i;
		sk_event const* mv = &last;
		if (a < r->nh && sk_before(&r->heap[a], mv)) { m = a; mv = &r->heap[a]; }
		if (b < r->nh && sk_before(&r->heap[b], mv)) { m = b; mv = &r->heap[b]; }
		if (m == i)
			break;
		r->heap[i] = r->heap[m];
		i = m;
	}
	if (r->nh) r->heap[i] = last;
	return top;
}

/* ------------------------------------------------------------------------ */
/*  Corners                                                                 */
/* ------------------------------------------------------------------------ */

/* a new corner at p, time t, between planes in and out, linked between
   prev and next; SK_NONE on failure */
static unsigned int sk_corner_new(sk_run* r, double const* p, double t, unsigned int in, unsigned int out,
	unsigned int prev, unsigned int next, unsigned int node, int* spike)
{
	sk_corner* c;
	unsigned int i;
	if (!sk_grow((void**)&r->c, &r->capc, r->nc + 1, sizeof(sk_corner)) ||
		!sk_grow((void**)&r->tag, &r->captag, r->nc + 1, sizeof(unsigned int)))
	{
		r->st = QAWS_STATUS_ALLOCATION_FAILURE;
		return SK_NONE;
	}
	i = r->nc++;
	r->tag[i] = 0;
	c = &r->c[i];
	c->p[0] = p[0]; c->p[1] = p[1]; c->t0 = t;
	c->in = in; c->out = out;
	c->prev = prev == SK_NONE ? i : prev;
	c->next = next == SK_NONE ? i : next;
	c->node = node;
	c->alive = 1;
	*spike = sk_velocity(&r->s->planes[in], &r->s->planes[out], c->v);
	if (prev != SK_NONE) r->c[prev].next = i;
	if (next != SK_NONE) r->c[next].prev = i;
	r->refs[out]++;
	return i;
}

/* the corner ends at node */
static void sk_kill(sk_run* r, unsigned int i, unsigned int node)
{
	sk_corner* c = &r->c[i];
	if (!c->alive)
		return;
	c->alive = 0;
	r->refs[c->out]--;
	sk_arc(r, c->node, node, c->in, c->out);
}

/* the events of a new corner: meeting either neighbour, and when reflex
   every plane ahead of it */
static void sk_candidates(sk_run* r, unsigned int i)
{
	sk_plane const* pl = r->s->planes;
	sk_corner const* c = &r->c[i];
	double p[2], t, din[2], dout[2], here[2];
	unsigned int j = c->next, k = c->prev, e;
	if (j != i && r->c[j].out != c->in && sk_solve3(&pl[c->in], &pl[c->out], &pl[r->c[j].out], p, &t) && t >= r->now - r->eps_t)
		sk_push(r, t, p, i);
	if (k != i && r->c[k].in != c->out && sk_solve3(&pl[r->c[k].in], &pl[c->in], &pl[c->out], p, &t) && t >= r->now - r->eps_t)
		sk_push(r, t, p, i);
	sk_dir(&pl[c->in], din);
	sk_dir(&pl[c->out], dout);
	if (din[0] * dout[1] - din[1] * dout[0] >= 0)
		return;
	sk_pos(c, r->now, here);
	for (e = 0; e < r->s->nplanes; e++)
	{
		sk_plane const* q = &pl[e];
		if (!r->refs[e] || e == c->in || e == c->out)
			continue;
		/* ahead of the plane's wavefront, and closing on it */
		if (q->n[0] * here[0] + q->n[1] * here[1] - q->d - q->w * r->now < -r->eps)
			continue;
		if (q->n[0] * c->v[0] + q->n[1] * c->v[1] - q->w >= 0)
			continue;
		if (sk_solve3(&pl[c->in], &pl[c->out], q, p, &t) && t >= r->now - r->eps_t)
			sk_push(r, t, p, i);
	}
}

static void sk_settle(sk_run* r, unsigned int i, int spike, double t);
static void sk_event_at(sk_run* r, double const* p, double t);

/* a corner between antiparallel edges at time t: the edges lie on one line,
   and the corner runs along it to the nearer end; both ends at once when
   they are as near */
static void sk_spike(sk_run* r, unsigned int i, double t)
{
	unsigned int guard = 0;
	while (r->st == QAWS_STATUS_OK && guard++ < 1000000)
	{
		sk_corner* c = &r->c[i];
		unsigned int a = c->prev, b = c->next, m, node, in = c->in, out = c->out;
		double pa[2], pb[2], here[2], sa, sb;
		int sp = 0;
		sk_pos(c, t, here);
		if (a == i)
		{
			sk_kill(r, i, sk_node(r, here, t));
			return;
		}
		sk_pos(&r->c[a], t, pa);
		sk_pos(&r->c[b], t, pb);
		if (a == b)
		{
			/* the loop is these two corners */
			node = sk_node(r, pa, t);
			sk_kill(r, i, node);
			sk_kill(r, a, node);
			return;
		}
		sa = hypot(pa[0] - here[0], pa[1] - here[1]);
		sb = hypot(pb[0] - here[0], pb[1] - here[1]);
		if (fabs(sa - sb) <= r->eps)
		{
			unsigned int aa = r->c[a].prev, bb = r->c[b].next, ain = r->c[a].in, bout = r->c[b].out;
			double q[2];
			q[0] = 0.5 * (pa[0] + pb[0]); q[1] = 0.5 * (pa[1] + pb[1]);
			node = sk_node(r, q, t);
			sk_kill(r, i, node);
			sk_kill(r, a, node);
			sk_kill(r, b, node);
			if (aa == b)
				return;                    /* the loop was these three */
			if (ain == bout)
			{
				r->c[aa].next = bb;
				r->c[bb].prev = aa;
				return;
			}
			m = sk_corner_new(r, q, t, ain, bout, aa, bb, node, &sp);
		}
		else if (sb < sa)
		{
			unsigned int bb = r->c[b].next, bout = r->c[b].out;
			node = sk_node(r, pb, t);
			sk_kill(r, i, node);
			sk_kill(r, b, node);
			if (bout == in)
			{
				r->c[a].next = bb;
				r->c[bb].prev = a;
				return;
			}
			m = sk_corner_new(r, pb, t, in, bout, a, bb, node, &sp);
		}
		else
		{
			unsigned int aa = r->c[a].prev, ain = r->c[a].in;
			node = sk_node(r, pa, t);
			sk_kill(r, i, node);
			sk_kill(r, a, node);
			if (ain == out)
			{
				r->c[aa].next = b;
				r->c[b].prev = aa;
				return;
			}
			m = sk_corner_new(r, pa, t, ain, out, aa, b, node, &sp);
		}
		if (m == SK_NONE)
			return;
		if (!sp)
		{
			sk_settle(r, m, 0, t);
			return;
		}
		i = m;
	}
}

static void sk_settle(sk_run* r, unsigned int i, int spike, double t)
{
	sk_corner* c;
	double q[2];
	if (i == SK_NONE || !r->c[i].alive)
		return;
	c = &r->c[i];
	sk_pos(c, t, q);
	if (c->next == i)
	{
		sk_kill(r, i, sk_node(r, q, t));
		return;
	}
	if (spike)
	{
		sk_spike(r, i, t);
		return;
	}
	if (c->next == c->prev)
	{
		/* two corners left, on the same two planes: one point */
		unsigned int node = sk_node(r, q, t), o = c->next;
		sk_kill(r, i, node);
		sk_kill(r, o, node);
		return;
	}
	/* other corners already here (a spike run to an event's point): one
	   event with them; corners born with this one are not counted */
	{
		unsigned int j;
		for (j = 0; j < r->nc; j++)
		{
			sk_corner const* o = &r->c[j];
			double w[2];
			if (j == i || !o->alive || (o->t0 == c->t0 && o->p[0] == c->p[0] && o->p[1] == c->p[1]))
				continue;
			sk_pos(o, t, w);
			if (hypot(w[0] - q[0], w[1] - q[1]) <= r->eps)
			{
				sk_event_at(r, q, t);
				return;
			}
		}
	}
	sk_candidates(r, i);
}

/* ------------------------------------------------------------------------ */
/*  The generalised intersection event                                      */
/* ------------------------------------------------------------------------ */

typedef struct sk_chain
{
	unsigned int in, out;       /* planes: the first edge, the last */
	unsigned int cin, cout;     /* the corners before the first edge and after the last */
	double ang;                 /* of the ray back along the first edge */
} sk_chain;

static int sk_cmp_chain(void const* a, void const* b)
{
	double x = ((sk_chain const*)a)->ang, y = ((sk_chain const*)b)->ang;
	return x < y ? -1 : (x > y ? 1 : 0);
}

static void sk_event_at(sk_run* r, double const* p, double t)
{
	sk_plane const* pl = r->s->planes;
	unsigned int i, ncons = 0, nsplit = 0, nch = 0, node, stamp = ++r->stamp, k;
	unsigned int *cons = NULL, *split = NULL;
	sk_chain* ch = NULL;
	unsigned int* fresh = NULL;
	int* spikes = NULL;
	cons = (unsigned int*)malloc(sizeof(unsigned int) * (r->nc + 1));
	split = (unsigned int*)malloc(sizeof(unsigned int) * (r->nc + 1));
	if (!cons || !split)
	{
		r->st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	/* the corners at p */
	for (i = 0; i < r->nc; i++)
	{
		double q[2];
		if (!r->c[i].alive)
			continue;
		sk_pos(&r->c[i], t, q);
		if (hypot(q[0] - p[0], q[1] - p[1]) <= r->eps)
		{
			r->tag[i] = stamp;
			cons[ncons++] = i;
		}
	}
	/* the edges through p (between two corners not at p) */
	for (i = 0; i < r->nc; i++)
	{
		sk_corner const* c = &r->c[i];
		sk_plane const* e;
		double a[2], b[2], dir[2], s, l;
		if (!c->alive || r->tag[i] == stamp || r->tag[c->next] == stamp || c->next == i)
			continue;
		e = &pl[c->out];
		if (fabs(e->n[0] * p[0] + e->n[1] * p[1] - e->d - e->w * t) > r->eps)
			continue;
		sk_pos(c, t, a);
		sk_pos(&r->c[c->next], t, b);
		sk_dir(e, dir);
		s = (p[0] - a[0]) * dir[0] + (p[1] - a[1]) * dir[1];
		l = (b[0] - a[0]) * dir[0] + (b[1] - a[1]) * dir[1];
		if (s > r->eps && s < l - r->eps)
			split[nsplit++] = i;
	}
	if (ncons < 2 && !(ncons == 1 && nsplit > 0))
		goto done;
	node = sk_node(r, p, t);
	ch = (sk_chain*)malloc(sizeof(sk_chain) * (ncons + nsplit + 1));
	fresh = (unsigned int*)malloc(sizeof(unsigned int) * (ncons + nsplit + 1));
	spikes = (int*)malloc(sizeof(int) * (ncons + nsplit + 1));
	if (!ch || !fresh || !spikes)
	{
		r->st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	/* runs of consecutive corners at p; a loop all at p just ends */
	for (k = 0; k < ncons; k++)
	{
		unsigned int a = cons[k], b;
		if (r->c[r->c[a].prev].alive && r->tag[r->c[a].prev] == stamp && r->c[a].prev != a)
		{
			/* not a run's first, unless the loop is all here */
			unsigned int w = r->c[a].prev, n = 0;
			while (w != a && r->tag[w] == stamp && n++ < r->nc) w = r->c[w].prev;
			if (w != a)
				continue;
			/* the whole loop: take it once, from its lowest index */
			for (w = r->c[a].next; w != a; w = r->c[w].next)
				if (w < a) break;
			if (w != a)
				continue;
			ch[nch].in = SK_NONE;
			nch++;
			continue;
		}
		b = a;
		while (r->tag[r->c[b].next] == stamp && r->c[b].next != a)
			b = r->c[b].next;
		ch[nch].in = r->c[a].in;
		ch[nch].out = r->c[b].out;
		ch[nch].cin = r->c[a].prev;
		ch[nch].cout = r->c[b].next;
		nch++;
	}
	for (k = 0; k < nsplit; k++)
	{
		unsigned int d = split[k];
		ch[nch].in = ch[nch].out = r->c[d].out;
		ch[nch].cin = d;
		ch[nch].cout = r->c[d].next;
		nch++;
	}
	for (k = 0; k < ncons; k++)
		sk_kill(r, cons[k], node);
	/* drop the closed loops, then the chains by angle round p */
	for (i = k = 0; i < nch; i++)
		if (ch[i].in != SK_NONE)
		{
			double dir[2];
			ch[k] = ch[i];
			sk_dir(&pl[ch[k].in], dir);
			ch[k].ang = atan2(-dir[1], -dir[0]);
			k++;
		}
	nch = k;
	qsort(ch, nch, sizeof(sk_chain), sk_cmp_chain);
	/* the last edge of each chain joined to the first of the next */
	for (i = 0; i < nch; i++)
	{
		sk_chain const* a = &ch[i];
		sk_chain const* b = &ch[(i + 1) % nch];
		fresh[i] = SK_NONE;
		spikes[i] = 0;
		if (b->in == a->out)
		{
			/* one plane on both sides: no corner */
			r->c[b->cin].next = a->cout;
			r->c[a->cout].prev = b->cin;
			continue;
		}
		fresh[i] = sk_corner_new(r, p, t, b->in, a->out, b->cin, a->cout, node, &spikes[i]);
	}
	for (i = 0; i < nch && r->st == QAWS_STATUS_OK; i++)
		if (fresh[i] != SK_NONE && r->c[fresh[i]].alive)
			sk_settle(r, fresh[i], spikes[i], t);
done:
	free(cons); free(split); free(ch); free(fresh); free(spikes);
}

/* ------------------------------------------------------------------------ */
/*  The run                                                                 */
/* ------------------------------------------------------------------------ */

static void sk_run_free(sk_run* r)
{
	free(r->c); free(r->refs); free(r->heap); free(r->tag);
}

static qaws_status sk_simulate(qaws_skeleton* s, int record, double stop, sk_run* r)
{
	unsigned int l, i;
	memset(r, 0, sizeof(*r));
	r->s = s;
	r->record = record;
	r->stop = stop;
	r->eps = s->tol;
	r->eps_t = s->tol;
	r->st = QAWS_STATUS_OK;
	r->refs = (unsigned int*)calloc(s->nplanes + 1, sizeof(unsigned int));
	if (!r->refs)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	{
		double wmax = 0;
		for (i = 0; i < s->nplanes; i++) if (s->planes[i].w > wmax) wmax = s->planes[i].w;
		if (wmax > 0) r->eps_t = s->tol / wmax;
	}
	/* the input's corners, at their nodes */
	for (l = 0; l < s->nloops && r->st == QAWS_STATUS_OK; l++)
	{
		unsigned int f = s->loop_first[l], g = s->loop_first[l + 1], first = r->nc;
		for (i = f; i < g; i++)
		{
			unsigned int in = i == f ? g - 1 : i - 1, ci;
			sk_plane const* a = &s->planes[i];
			int sp;
			double p[2];
			if (record)
			{
				p[0] = (double)s->nodes[a->na].position.x;
				p[1] = (double)s->nodes[a->na].position.y;
			}
			else
			{
				/* the corner from the planes at time 0 */
				sk_plane const* b = &s->planes[in];
				double det = b->n[0] * a->n[1] - b->n[1] * a->n[0];
				if (fabs(det) > 1e-12)
				{
					p[0] = (b->d * a->n[1] - b->n[1] * a->d) / det;
					p[1] = (b->n[0] * a->d - b->d * a->n[0]) / det;
				}
				else
				{
					p[0] = a->n[0] * a->d;
					p[1] = a->n[1] * a->d;
				}
			}
			ci = sk_corner_new(r, p, 0, in, i, i == f ? SK_NONE : r->nc - 1, SK_NONE, a->na, &sp);
			if (ci == SK_NONE)
				break;
			(void)sp;
		}
		if (r->st == QAWS_STATUS_OK && r->nc > first)
		{
			r->c[r->nc - 1].next = first;
			r->c[first].prev = r->nc - 1;
		}
	}
	/* when the input stores corners with the actual node positions, use them
	   (record); else they come from the planes */
	for (i = 0; i < r->nc && r->st == QAWS_STATUS_OK; i++)
		sk_candidates(r, i);
	while (r->nh && r->st == QAWS_STATUS_OK)
	{
		sk_event e = sk_pop(r);
		if (e.t > stop)
			break;
		if (!r->c[e.corner].alive)
			continue;
		if (e.t > r->now)
		{
			r->now = e.t;
			if (record) r->node_from = s->nnodes;
		}
		{
			double q[2];
			sk_pos(&r->c[e.corner], r->now, q);
			sk_event_at(r, q, r->now);
		}
	}
	return r->st;
}

/* ------------------------------------------------------------------------ */
/*  The input                                                               */
/* ------------------------------------------------------------------------ */

/* the knot span of a B-spline holding t (0 for a Bezier) */
static unsigned int sk_span(qaws_curve const* c, double t)
{
	qaws_scalar buf[256], *k = buf;
	unsigned int n = 0, i, span = 0, deg;
	if (c->kind != QAWS_CURVE_KIND_BSPLINE && c->kind != QAWS_CURVE_KIND_NURBS)
		return 0;
	deg = c->degree;
	if (qaws_bspline_get_knots(c, buf, 256, &n) == QAWS_STATUS_BUFFER_TOO_SMALL)
	{
		k = (qaws_scalar*)malloc(sizeof(qaws_scalar) * n);
		if (!k) return 0;
		qaws_bspline_get_knots(c, k, n, &n);
	}
	for (i = deg + 1; i + 1 + deg < n; i++)
		if ((double)k[i] > (double)k[i - 1] && (double)k[i] <= t)
			span++;
	if (k != buf) free(k);
	return span;
}

static int sk_is_line(qaws_curve const* c)
{
	return c->degree == 1 && (c->kind == QAWS_CURVE_KIND_BSPLINE || c->kind == QAWS_CURVE_KIND_BEZIER);
}

void qaws_skeleton_destroy(qaws_skeleton* s)
{
	if (!s)
		return;
	free(s->planes);
	free(s->loop_first);
	free(s->nodes);
	free(s->arcs);
	free(s->face_first);
	free(s->face_nodes);
	free(s);
}

/* the input curves of output curve i of a clip path: the source both its
   vertices lie on */
static void sk_source(qaws_clip_vertex const* v, unsigned int n, unsigned int i, unsigned int* path, unsigned int* curve)
{
	qaws_clip_vertex const* a = &v[i];
	qaws_clip_vertex const* b = &v[(i + 1) % n];
	*path = a->path_a;
	*curve = a->curve_a;
	if ((a->path_a == b->path_a && a->curve_a == b->curve_a) || (a->path_a == b->path_b && a->curve_a == b->curve_b))
		return;
	if (a->curve_b != QAWS_CLIP_NONE_INDEX &&
		((a->path_b == b->path_a && a->curve_b == b->curve_a) || (a->path_b == b->path_b && a->curve_b == b->curve_b)))
	{
		*path = a->path_b;
		*curve = a->curve_b;
	}
}

/* the faces: each plane's arcs walked from its edge's end back to its start */
static qaws_status sk_faces(qaws_skeleton* s)
{
	unsigned int f, i, total = 0, *cnt, *first, *list, *used;
	cnt = (unsigned int*)calloc(s->nplanes + 1, sizeof(unsigned int));
	first = (unsigned int*)calloc(s->nplanes + 1, sizeof(unsigned int));
	list = (unsigned int*)malloc(sizeof(unsigned int) * (2 * s->narcs + 1));
	used = (unsigned int*)calloc(s->narcs + 1, sizeof(unsigned int));
	s->face_first = (unsigned int*)malloc(sizeof(unsigned int) * (s->nplanes + 1));
	s->face_nodes = (unsigned int*)malloc(sizeof(unsigned int) * (2 * s->narcs + 2 * s->nplanes + 1));
	if (!cnt || !first || !list || !used || !s->face_first || !s->face_nodes)
	{
		free(cnt); free(first); free(list); free(used);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < s->narcs; i++)
	{
		cnt[s->arcs[i].face_a]++;
		if (s->arcs[i].face_b != s->arcs[i].face_a) cnt[s->arcs[i].face_b]++;
	}
	for (f = 0; f < s->nplanes; f++) { first[f + 1] = first[f] + cnt[f]; cnt[f] = 0; }
	for (i = 0; i < s->narcs; i++)
	{
		unsigned int a = s->arcs[i].face_a, b = s->arcs[i].face_b;
		list[first[a] + cnt[a]++] = i;
		if (b != a) list[first[b] + cnt[b]++] = i;
	}
	for (f = 0; f < s->nplanes; f++)
	{
		unsigned int cur = s->planes[f].nb, guard = 0;
		s->face_first[f] = total;
		s->face_nodes[total++] = s->planes[f].na;
		s->face_nodes[total++] = cur;
		while (cur != s->planes[f].na && guard++ <= first[f + 1] - first[f])
		{
			unsigned int j, best = SK_NONE;
			for (j = first[f]; j < first[f + 1]; j++)
			{
				qaws_skeleton_arc const* a = &s->arcs[list[j]];
				if (used[list[j]] == f + 1)
					continue;
				if (a->node_a == cur || a->node_b == cur)
				{
					best = list[j];
					break;
				}
			}
			if (best == SK_NONE)
				break;
			used[best] = f + 1;
			cur = s->arcs[best].node_a == cur ? s->arcs[best].node_b : s->arcs[best].node_a;
			if (cur != s->planes[f].na)
				s->face_nodes[total++] = cur;
		}
	}
	s->face_first[s->nplanes] = total;
	free(cnt); free(first); free(list); free(used);
	return QAWS_STATUS_OK;
}

qaws_status qaws_skeleton_create(qaws_skeleton_desc const* desc, qaws_skeleton** out)
{
	qaws_skeleton* s;
	qaws_clip_result* region = NULL;
	qaws_status st;
	unsigned int np, pi, capplanes = 0;
	double lo[2] = { 1e300, 1e300 }, hi[2] = { -1e300, -1e300 };
	sk_run r;
	if (!desc || !out || (desc->path_count && !desc->paths))
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	st = qaws_clip_boolean(QAWS_CLIP_UNION, desc->fill_rule, desc->paths, desc->path_count, NULL, 0, &region);
	if (st != QAWS_STATUS_OK)
		return st;
	s = (qaws_skeleton*)calloc(1, sizeof(qaws_skeleton));
	if (!s)
	{
		qaws_clip_result_destroy(region);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	np = qaws_clip_result_get_path_count(region);
	s->loop_first = (unsigned int*)malloc(sizeof(unsigned int) * (np + 1));
	if (!s->loop_first)
		st = QAWS_STATUS_ALLOCATION_FAILURE;
	/* the planes, loop by loop; a node at every input corner */
	for (pi = 0; pi < np && st == QAWS_STATUS_OK; pi++)
	{
		qaws_path_2d path;
		qaws_clip_vertex const* vx = NULL;
		unsigned int nv = 0, ci, first = s->nplanes, first_node = s->nnodes;
		qaws_clip_result_get_path(region, pi, &path);
		qaws_clip_result_get_vertices(region, 0, pi, &vx, &nv);
		s->loop_first[s->nloops] = first;
		for (ci = 0; ci < path.curve_count && st == QAWS_STATUS_OK; ci++)
		{
			qaws_curve const* c = path.curves[ci];
			qaws_scalar buf[2 * 64], *cp = buf;
			unsigned int n = 0, j, sp = 0, sc = 0;
			if (!sk_is_line(c))
			{
				st = QAWS_STATUS_UNSUPPORTED_OPERATION;
				break;
			}
			if (qaws_curve_get_control_points(c, buf, 64, &n) == QAWS_STATUS_BUFFER_TOO_SMALL)
			{
				cp = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * n);
				if (!cp) { st = QAWS_STATUS_ALLOCATION_FAILURE; break; }
				qaws_curve_get_control_points(c, cp, n, &n);
			}
			if (vx && nv == path.curve_count)
				sk_source(vx, nv, ci, &sp, &sc);
			for (j = 0; j + 1 < n && st == QAWS_STATUS_OK; j++)
			{
				double a[2], b[2], len;
				sk_plane* pl;
				a[0] = (double)cp[2 * j]; a[1] = (double)cp[2 * j + 1];
				b[0] = (double)cp[2 * j + 2]; b[1] = (double)cp[2 * j + 3];
				len = hypot(b[0] - a[0], b[1] - a[1]);
				if (!(len > 0))
					continue;
				if (!sk_grow((void**)&s->planes, &capplanes, s->nplanes + 1, sizeof(sk_plane)))
				{
					st = QAWS_STATUS_ALLOCATION_FAILURE;
					break;
				}
				pl = &s->planes[s->nplanes];
				pl->n[0] = -(b[1] - a[1]) / len;
				pl->n[1] = (b[0] - a[0]) / len;
				pl->d = pl->n[0] * a[0] + pl->n[1] * a[1];
				pl->w = 1;
				pl->path = sp;
				pl->curve = sc;
				pl->span = 0;
				if (sp < desc->path_count && sc < desc->paths[sp].curve_count)
				{
					qaws_curve const* src = desc->paths[sp].curves[sc];
					qaws_vec2 mid;
					qaws_scalar t = 0;
					mid.x = (qaws_scalar)(0.5 * (a[0] + b[0]));
					mid.y = (qaws_scalar)(0.5 * (a[1] + b[1]));
					if (qaws_curve_find_closest_parameter_2d(src, mid, &t) == QAWS_STATUS_OK)
						pl->span = sk_span(src, (double)t);
				}
				if (desc->weight_fn)
				{
					pl->w = (double)desc->weight_fn(desc->weight_user, pl->path, pl->curve, pl->span);
					if (!(pl->w > 0))
						st = QAWS_STATUS_INVALID_ARGUMENT;
				}
				/* the node at the edge's start */
				if (!sk_grow((void**)&s->nodes, &s->capnodes, s->nnodes + 1, sizeof(qaws_skeleton_node)))
				{
					st = QAWS_STATUS_ALLOCATION_FAILURE;
					break;
				}
				s->nodes[s->nnodes].position.x = (qaws_scalar)a[0];
				s->nodes[s->nnodes].position.y = (qaws_scalar)a[1];
				s->nodes[s->nnodes].time = 0;
				pl->na = s->nnodes++;
				s->nplanes++;
				if (a[0] < lo[0]) lo[0] = a[0];
				if (a[1] < lo[1]) lo[1] = a[1];
				if (a[0] > hi[0]) hi[0] = a[0];
				if (a[1] > hi[1]) hi[1] = a[1];
			}
			if (cp != buf) free(cp);
		}
		if (st == QAWS_STATUS_OK && s->nplanes - first < 3)
		{
			s->nplanes = first;               /* nothing to sweep */
			s->nnodes = first_node;
		}
		if (st == QAWS_STATUS_OK && s->nplanes > first)
		{
			unsigned int k;
			for (k = first; k < s->nplanes; k++)
				s->planes[k].nb = s->planes[k + 1 < s->nplanes ? k + 1 : first].na;
			s->nloops++;
		}
	}
	qaws_clip_result_destroy(region);
	if (st == QAWS_STATUS_OK)
	{
		double ext = hypot(hi[0] - lo[0], hi[1] - lo[1]);
		s->loop_first[s->nloops] = s->nplanes;
		s->tol = desc->tolerance > 0 ? (double)desc->tolerance : SK_TOL_REL * (ext > 0 ? ext : 1);
		st = sk_simulate(s, 1, 1e300, &r);
		if (st == QAWS_STATUS_OK)
		{
			unsigned int i;
			for (i = 0; i < r.nc; i++)
				if (r.c[i].alive)
				{
					st = QAWS_STATUS_NUMERICAL_FAILURE;   /* a loop that never closed */
					break;
				}
		}
		sk_run_free(&r);
	}
	if (st == QAWS_STATUS_OK)
		st = sk_faces(s);
	if (st != QAWS_STATUS_OK)
	{
		qaws_skeleton_destroy(s);
		return st;
	}
	*out = s;
	return QAWS_STATUS_OK;
}

unsigned int qaws_skeleton_get_node_count(qaws_skeleton const* s) { return s ? s->nnodes : 0; }
unsigned int qaws_skeleton_get_arc_count(qaws_skeleton const* s) { return s ? s->narcs : 0; }
unsigned int qaws_skeleton_get_face_count(qaws_skeleton const* s) { return s ? s->nplanes : 0; }
qaws_scalar qaws_skeleton_get_max_time(qaws_skeleton const* s) { return s ? (qaws_scalar)s->tmax : 0; }

qaws_status qaws_skeleton_get_node(qaws_skeleton const* s, unsigned int i, qaws_skeleton_node* out)
{
	if (!s || !out || i >= s->nnodes)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = s->nodes[i];
	return QAWS_STATUS_OK;
}

qaws_status qaws_skeleton_get_arc(qaws_skeleton const* s, unsigned int i, qaws_skeleton_arc* out)
{
	if (!s || !out || i >= s->narcs)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = s->arcs[i];
	return QAWS_STATUS_OK;
}

qaws_status qaws_skeleton_get_face(qaws_skeleton const* s, unsigned int i, qaws_skeleton_face* out)
{
	if (!s || !out || i >= s->nplanes)
		return QAWS_STATUS_INVALID_ARGUMENT;
	out->nodes = s->face_nodes + s->face_first[i];
	out->node_count = s->face_first[i + 1] - s->face_first[i];
	out->path = s->planes[i].path;
	out->curve = s->planes[i].curve;
	out->span = s->planes[i].span;
	out->weight = (qaws_scalar)s->planes[i].w;
	return QAWS_STATUS_OK;
}

qaws_status qaws_skeleton_offset(qaws_skeleton const* s, qaws_scalar time, qaws_clip_result** out)
{
	sk_run r;
	qaws_status st;
	unsigned int i, nloops = 0, npts = 0;
	qaws_scalar* pts = NULL;
	unsigned int* lp = NULL;
	qaws_curve** curves = NULL;
	qaws_curve const** cv = NULL;
	qaws_path_2d* paths = NULL;
	if (!s || !out || !(time >= 0))
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	st = sk_simulate((qaws_skeleton*)s, 0, (double)time, &r);
	if (st == QAWS_STATUS_OK)
	{
		/* the alive loops, at the time */
		unsigned int* seen = (unsigned int*)calloc(r.nc + 1, sizeof(unsigned int));
		pts = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * (r.nc + 1));
		lp = (unsigned int*)malloc(sizeof(unsigned int) * (r.nc + 2));
		if (!seen || !pts || !lp)
			st = QAWS_STATUS_ALLOCATION_FAILURE;
		for (i = 0; i < r.nc && st == QAWS_STATUS_OK; i++)
		{
			unsigned int j = i, n = 0;
			if (!r.c[i].alive || seen[i])
				continue;
			lp[nloops] = npts;
			do
			{
				double q[2];
				seen[j] = 1;
				sk_pos(&r.c[j], (double)time, q);
				pts[2 * npts] = (qaws_scalar)q[0];
				pts[2 * npts + 1] = (qaws_scalar)q[1];
				npts++;
				n++;
				j = r.c[j].next;
			} while (j != i && n <= r.nc);
			if (n >= 3) nloops++;
			else npts = lp[nloops];
		}
		lp[nloops] = npts;
		free(seen);
	}
	sk_run_free(&r);
	if (st == QAWS_STATUS_OK)
	{
		curves = (qaws_curve**)calloc(nloops + 1, sizeof(qaws_curve*));
		cv = (qaws_curve const**)malloc(sizeof(qaws_curve*) * (nloops + 1));
		paths = (qaws_path_2d*)malloc(sizeof(qaws_path_2d) * (nloops + 1));
		if (!curves || !cv || !paths)
			st = QAWS_STATUS_ALLOCATION_FAILURE;
		for (i = 0; i < nloops && st == QAWS_STATUS_OK; i++)
		{
			st = qaws_curve_create_polyline_2d(pts + 2 * lp[i], lp[i + 1] - lp[i], 1, &curves[i]);
			cv[i] = curves[i];
			paths[i].curves = &cv[i];
			paths[i].curve_count = 1;
			paths[i].closed = 1;
		}
		if (st == QAWS_STATUS_OK)
			st = qaws_clip_boolean(QAWS_CLIP_UNION, QAWS_FILL_POSITIVE, paths, nloops, NULL, 0, out);
		for (i = 0; curves && i < nloops; i++)
			if (curves[i]) qaws_curve_destroy(curves[i]);
	}
	free(pts); free(lp); free(curves); free(cv); free(paths);
	return st;
}
