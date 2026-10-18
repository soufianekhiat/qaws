/*
 * Min-cost flow by successive shortest paths (qaws_internal_flow.h).
 *
 * Residual arc r of arc a: 2 a runs from -> to with room capacity - flow at
 * cost c, 2 a + 1 runs back with room flow at cost -c. Potentials pi keep
 * every residual arc's reduced cost c + pi(u) - pi(v) non-negative: a
 * search stopped at distance D moves pi by min(d, D) everywhere, which is
 * pi += d - D on the nodes it settled and nothing elsewhere (a shift common
 * to all nodes changes no reduced cost).
 */

#include "qaws_internal_flow.h"
#include <stdlib.h>
#include <string.h>

typedef struct fl_item
{
	long long d;
	unsigned int n;
} fl_item;

struct qaws_flow
{
	unsigned int n, na, cap;
	unsigned int *from, *to;
	int *room, *flow;               /* per arc: its capacity, its flow */
	long long* cost;
	int* excess;                    /* per node: supply not yet sent (negative: demand not yet met) */
	long long *pi, *dist;
	unsigned int *first, *res;      /* the residual arcs leaving each node (CSR), for the first `built` arcs */
	unsigned int built;
	unsigned int *pred, *seen, *done, *settled, gen, nset;
	unsigned char* target;
	fl_item* heap;
};

void qaws_flow_destroy(qaws_flow* f)
{
	if (!f)
		return;
	free(f->from); free(f->to); free(f->room); free(f->flow); free(f->cost); free(f->excess); free(f->pi); free(f->dist);
	free(f->first); free(f->res); free(f->pred); free(f->seen); free(f->done); free(f->settled); free(f->target); free(f->heap);
	free(f);
}

qaws_status qaws_flow_create(unsigned int nodes, unsigned int arc_capacity, qaws_flow** out)
{
	qaws_flow* f;
	size_t n = (size_t)nodes + 1, a = (size_t)arc_capacity + 1;
	if (!out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	f = (qaws_flow*)calloc(1, sizeof *f);
	if (!f)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	f->n = nodes;
	f->cap = arc_capacity;
	f->from = (unsigned int*)malloc(sizeof(unsigned int) * a);
	f->to = (unsigned int*)malloc(sizeof(unsigned int) * a);
	f->room = (int*)malloc(sizeof(int) * a);
	f->flow = (int*)malloc(sizeof(int) * a);
	f->cost = (long long*)malloc(sizeof(long long) * a);
	f->excess = (int*)calloc(n, sizeof(int));
	f->pi = (long long*)calloc(n, sizeof(long long));
	f->dist = (long long*)malloc(sizeof(long long) * n);
	f->first = (unsigned int*)malloc(sizeof(unsigned int) * (n + 1));
	f->res = (unsigned int*)malloc(sizeof(unsigned int) * 2 * a);
	f->pred = (unsigned int*)malloc(sizeof(unsigned int) * n);
	f->seen = (unsigned int*)calloc(n, sizeof(unsigned int));
	f->done = (unsigned int*)calloc(n, sizeof(unsigned int));
	f->settled = (unsigned int*)malloc(sizeof(unsigned int) * n);
	f->target = (unsigned char*)calloc(n, 1);
	f->heap = (fl_item*)malloc(sizeof(fl_item) * (2 * a + 1));   /* (a search relaxes each residual arc at most once) */
	if (!f->from || !f->to || !f->room || !f->flow || !f->cost || !f->excess || !f->pi || !f->dist || !f->first || !f->res || !f->pred
		|| !f->seen || !f->done || !f->settled || !f->target || !f->heap)
	{
		qaws_flow_destroy(f);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	*out = f;
	return QAWS_STATUS_OK;
}

unsigned int qaws_flow_arc(qaws_flow* f, unsigned int from, unsigned int to, int capacity, long long cost)
{
	unsigned int a;
	if (f->na == f->cap || from >= f->n || to >= f->n)
		return QAWS_FLOW_NONE;
	a = f->na++;
	f->from[a] = from;
	f->to[a] = to;
	f->room[a] = capacity > 0 ? capacity : 0;
	f->flow[a] = 0;
	f->cost[a] = cost > 0 ? cost : 0;
	return a;
}

void qaws_flow_supply(qaws_flow* f, unsigned int node, int b)
{
	if (node < f->n)
		f->excess[node] += b;
}

int qaws_flow_get(qaws_flow const* f, unsigned int arc)
{
	return arc < f->na ? f->flow[arc] : 0;
}

static void fl_build(qaws_flow* f)
{
	unsigned int i, r;
	memset(f->first, 0, sizeof(unsigned int) * (f->n + 2));
	for (i = 0; i < f->na; i++)
	{
		f->first[f->from[i] + 2]++;
		f->first[f->to[i] + 2]++;
	}
	for (i = 0; i < f->n; i++)
		f->first[i + 2] += f->first[i + 1];
	for (r = 0; r < 2 * f->na; r++)
		f->res[f->first[(r & 1 ? f->to[r >> 1] : f->from[r >> 1]) + 1]++] = r;
	f->built = f->na;
}

static int fl_less(fl_item const* a, fl_item const* b)
{
	return a->d < b->d || (a->d == b->d && a->n < b->n);
}

static void fl_push(fl_item* h, unsigned int* n, long long d, unsigned int v)
{
	unsigned int i = (*n)++;
	fl_item it;
	it.d = d;
	it.n = v;
	while (i > 0 && fl_less(&it, &h[(i - 1) / 2]))
	{
		h[i] = h[(i - 1) / 2];
		i = (i - 1) / 2;
	}
	h[i] = it;
}

static fl_item fl_pop(fl_item* h, unsigned int* n)
{
	fl_item top = h[0], last = h[--(*n)];
	unsigned int i = 0, j;
	while ((j = 2 * i + 1) < *n)
	{
		if (j + 1 < *n && fl_less(&h[j + 1], &h[j]))
			j++;
		if (!fl_less(&h[j], &last))
			break;
		h[i] = h[j];
		i = j;
	}
	h[i] = last;
	return top;
}

/* Dijkstra on reduced costs over the residual arcs from src. With targets == 0 it stops at the first node with a
   demand; else at the marked node nearest by real cost (the reduced distance less pi(src), plus pi(target); pimin, the
   least pi of the targets, bounds when no nearer one can follow). Returns that node (QAWS_FLOW_NONE: none, nothing
   changed) and moves the potentials. */
static unsigned int fl_search(qaws_flow* f, unsigned int src, int targets, long long pimin)
{
	unsigned int nh = 0, best = QAWS_FLOW_NONE, i, k;
	long long D = 0, bestreal = 0;
	if (f->built != f->na)
		fl_build(f);
	if (++f->gen == 0)
	{
		memset(f->seen, 0, sizeof(unsigned int) * f->n);
		memset(f->done, 0, sizeof(unsigned int) * f->n);
		f->gen = 1;
	}
	f->nset = 0;
	f->seen[src] = f->gen;
	f->dist[src] = 0;
	f->pred[src] = QAWS_FLOW_NONE;
	fl_push(f->heap, &nh, 0, src);
	while (nh)
	{
		fl_item it = fl_pop(f->heap, &nh);
		unsigned int u = it.n;
		if (f->done[u] == f->gen || it.d > f->dist[u])
			continue;
		if (targets && best != QAWS_FLOW_NONE && it.d - f->pi[src] + pimin >= bestreal)
			break;
		f->done[u] = f->gen;
		f->settled[f->nset++] = u;
		if (u != src && (targets ? f->target[u] != 0 : f->excess[u] < 0))
		{
			long long real = it.d - f->pi[src] + f->pi[u];
			if (best == QAWS_FLOW_NONE || real < bestreal)
			{
				best = u;
				bestreal = real;
				D = it.d;
			}
			if (!targets)
				break;
		}
		for (k = f->first[u]; k < f->first[u + 1]; k++)
		{
			unsigned int r = f->res[k], a = r >> 1, v = r & 1 ? f->from[a] : f->to[a];
			int room = r & 1 ? f->flow[a] : f->room[a] - f->flow[a];
			long long nd;
			if (room <= 0 || f->done[v] == f->gen)
				continue;
			nd = it.d + (r & 1 ? -f->cost[a] : f->cost[a]) + f->pi[u] - f->pi[v];
			if (f->seen[v] != f->gen || nd < f->dist[v])
			{
				f->seen[v] = f->gen;
				f->dist[v] = nd;
				f->pred[v] = r;
				fl_push(f->heap, &nh, nd, v);
			}
		}
	}
	if (best == QAWS_FLOW_NONE)
		return best;
	for (i = 0; i < f->nset; i++)
	{
		unsigned int u = f->settled[i];
		if (f->dist[u] <= D)
			f->pi[u] += f->dist[u] - D;
	}
	return best;
}

/* the room along the path the last search found to t */
static int fl_path_room(qaws_flow const* f, unsigned int t)
{
	int m = 0x7FFFFFFF;
	while (f->pred[t] != QAWS_FLOW_NONE)
	{
		unsigned int r = f->pred[t], a = r >> 1;
		int room = r & 1 ? f->flow[a] : f->room[a] - f->flow[a];
		if (room < m)
			m = room;
		t = r & 1 ? f->to[a] : f->from[a];
	}
	return m;
}

/* amount sent along that path; returns its cost */
static long long fl_send(qaws_flow* f, unsigned int t, int amount)
{
	long long c = 0;
	unsigned int v = t;
	while (f->pred[v] != QAWS_FLOW_NONE)
	{
		unsigned int r = f->pred[v], a = r >> 1;
		if (r & 1)
		{
			f->flow[a] -= amount;
			c -= f->cost[a];
			v = f->to[a];
		}
		else
		{
			f->flow[a] += amount;
			c += f->cost[a];
			v = f->from[a];
		}
	}
	f->excess[v] -= amount;
	f->excess[t] += amount;
	return c * amount;
}

qaws_status qaws_flow_solve(qaws_flow* f, int* out_unrouted, long long* out_cost)
{
	unsigned int s, t, a;
	int progress = 1, left = 0;
	long long c = 0;
	if (!f)
		return QAWS_STATUS_INVALID_ARGUMENT;
	/* each excess in turn to its nearest demand (searching from all at once costs more: their balls overlap); a pass
	   that routes nothing ends it, what is left reaching no demand */
	while (progress)
	{
		progress = 0;
		for (s = 0; s < f->n; s++)
			while (f->excess[s] > 0 && (t = fl_search(f, s, 0, 0)) != QAWS_FLOW_NONE)
			{
				int m = fl_path_room(f, t);
				if (m > f->excess[s])
					m = f->excess[s];
				if (m > -f->excess[t])
					m = -f->excess[t];
				fl_send(f, t, m);
				progress = 1;
			}
	}
	for (s = 0; s < f->n; s++)
		left += f->excess[s] > 0 ? f->excess[s] : 0;
	for (a = 0; a < f->na; a++)
		c += f->flow[a] * f->cost[a];
	if (out_unrouted)
		*out_unrouted = left;
	if (out_cost)
		*out_cost = c;
	return QAWS_STATUS_OK;
}

qaws_status qaws_flow_augment(qaws_flow* f, unsigned int from, unsigned int const* targets, unsigned int target_count,
	unsigned int* out_reached, long long* out_cost)
{
	unsigned int i, t, any = 0;
	long long pimin = 0;
	if (out_reached)
		*out_reached = QAWS_FLOW_NONE;
	if (out_cost)
		*out_cost = 0;
	if (!f || from >= f->n || (target_count && !targets))
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < target_count; i++)
		if (targets[i] < f->n && targets[i] != from)
		{
			f->target[targets[i]] = 1;
			pimin = any++ ? (f->pi[targets[i]] < pimin ? f->pi[targets[i]] : pimin) : f->pi[targets[i]];
		}
	t = any ? fl_search(f, from, 1, pimin) : QAWS_FLOW_NONE;
	for (i = 0; i < target_count; i++)
		if (targets[i] < f->n)
			f->target[targets[i]] = 0;
	if (t == QAWS_FLOW_NONE)
		return QAWS_STATUS_OK;
	if (out_cost)
		*out_cost = fl_send(f, t, 1);
	else
		fl_send(f, t, 1);
	if (out_reached)
		*out_reached = t;
	return QAWS_STATUS_OK;
}
