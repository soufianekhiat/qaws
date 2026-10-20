#ifndef QAWS_INTERNAL_FLOW_H
#define QAWS_INTERNAL_FLOW_H

/*
 * Min-cost flow on a directed network by successive shortest paths:
 * Dijkstra on reduced costs (a binary heap), the node potentials updated
 * only where a search went, so each augmentation costs the region it
 * searched. Arcs of integer capacity and non-negative cost, integer
 * supplies at the nodes. Deterministic: ties broken by node and arc index.
 * Generic (its first user makes quads' lattice steps consistent).
 */

#include "../qaws_status.h"

#define QAWS_FLOW_NONE 0xFFFFFFFFu

typedef struct qaws_flow qaws_flow;

/* nodes 0 .. nodes - 1, room for arc_capacity arcs */
qaws_status qaws_flow_create(qaws_flow** out, unsigned int arc_capacity, unsigned int nodes);
/* an arc from -> to of capacity >= 0 and cost >= 0; returns its index (QAWS_FLOW_NONE when full) */
unsigned int qaws_flow_arc(qaws_flow* f, unsigned int from, unsigned int to, int capacity, long long cost);
/* b added to the node's supply (b > 0: an excess to send, b < 0: a demand) */
void qaws_flow_supply(qaws_flow* f, unsigned int node, int b);
/* the cheapest flow meeting the supplies as far as it can: out_unrouted, the excess left (0: all routed); out_cost,
   the flow's total cost; either may be NULL */
qaws_status qaws_flow_solve(int* out_unrouted, long long* out_cost, qaws_flow* f);
/* one more unit from `from` along the cheapest path of the residual network to the nearest of `targets` (the flow
   stays the cheapest for its ends): out_reached, that target (QAWS_FLOW_NONE: none reachable, nothing changed);
   out_cost, the path's cost (may be NULL) */
qaws_status qaws_flow_augment(unsigned int* out_reached, long long* out_cost, qaws_flow* f, unsigned int const* targets,
	unsigned int target_count, unsigned int from);
/* the flow on an arc */
int qaws_flow_get(qaws_flow const* f, unsigned int arc);
void qaws_flow_destroy(qaws_flow* f);

#endif /* QAWS_INTERNAL_FLOW_H */
