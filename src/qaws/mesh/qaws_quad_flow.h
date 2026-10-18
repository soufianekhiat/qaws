#ifndef QAWS_QUAD_FLOW_H
#define QAWS_QUAD_FLOW_H

/*
 * The lattice steps of INSTANT quads made consistent (QuadriFlow: Huang,
 * Zhou, Niessner, Shewchuk, Guibas, SGP 2018). Per corner c of the dense
 * mesh (the edge from its vertex a to the next, b): T, the step from a's
 * lattice point to b's in whole steps of a's frame (q, n x q); K, the
 * quarter turns taking b's frame onto a's (coordinates in b's frame are
 * those in a's turned by K); Tc, where each end sits in its cell. A face
 * whose steps do not close holds a dislocation (an extra half row of the
 * lattice ends there). The steps are mended by a min-cost flow over the
 * faces (qaws_internal_flow.h), the cones put on lattice points and then on
 * vertices, faces the changes turned over untangled, the positions refitted
 * to the mended steps.
 */

#include "qaws_mesh_internal.h"

/* the lattice coordinates c turned by k quarter turns counter-clockwise */
void qf_turn(int const* c, int k, int* o);

/* T and K per corner (2 and 1 each), the two corners of an edge kept opposite; Tc (4 doubles per corner, the edges'
   first corners only): a's place in its cell and b's, in a's frame, in steps (x less its lattice point o) */
void qf_steps(qaws_hm const* dm, mf_level const* L, double rho, int* T, unsigned char* K, double* Tc);

/* face f's steps summed in its first vertex's frame, S (zero: the lattice closes round f); fr per corner: the quarter
   turns from its vertex's frame to the face's. Returns the quarter turns of the frames round f, 0 to 3 (nonzero: a
   singularity of the cross field). */
int qf_closure(qaws_hm const* dm, unsigned int f, int const* T, unsigned char const* K, unsigned char* fr, int* S);

/* face f's lattice image: e (3 x 2), its steps in the frame of the corner from which they close; returns that corner
   (0 to 2), or -1 when they close from none (a dislocation, or a cone off its vertices) */
int qf_image(qaws_hm const* dm, unsigned int f, int const* T, unsigned char const* K, int* e);

/* the lattice axis (0: q, 1: n x q of its vertex) corner c's edge runs along */
int qf_axis(qaws_hm const* dm, mf_level const* L, unsigned int c);

typedef struct qf_report
{
	unsigned int supply;            /* the dislocations' units before the flow */
	unsigned int cut;               /* edge components the combed frames could not balance (kept as they are) */
	unsigned int unrouted;          /* units the flow could not route */
	unsigned int changed;           /* units of change on the steps, closure and parity together */
	unsigned int paths;             /* parity paths (from an odd cone to another or to the border) */
	unsigned int twos;              /* steps of two or more left (extraction drops them) */
	unsigned int flips;             /* faces whose lattice image turned over after the flow */
	unsigned int flips_left;        /* and after vertices were moved to undo them */
	unsigned int cones;             /* cones of even sum put on a vertex of their face */
	long long cost;                 /* the closure flow's cost, in millionths */
} qf_report;

/* The steps T made to close round every face where the field does not turn: a min-cost flow on the faces (two nodes
   per face, the rows of its equation in frames combed across the faces; the border a sink, the singular faces one at
   a price), a change costing how far it takes an end from its own lattice point (Tc). Then each cone of a singular
   face whose frames turn by a quarter put on a lattice point (its sum's parity made even) by a path to another odd
   cone or the border, and on a vertex of its face where moving that vertex's lattice point there harms little; faces
   the changes turned over in the lattice untangled by moving vertices (QuadriFlow's flip removal). Returns the
   dislocations left (HM_NONE on allocation failure); out_odd: the cones left odd; out_report may be NULL. */
unsigned int qf_consistent(qaws_hm const* dm, mf_level const* L, int* T, unsigned char const* K, double const* Tc, unsigned int* out_odd,
	qf_report* out_report);

/* the level-0 positions o re-solved with T fixed (Gauss-Seidel sweeps: each o the mean of its neighbours' o, unfolded,
   less the steps), kept in the tangent plane, on the line for MF_ON_LINE, at x for MF_PINNED; not rounded.
   shortcut: sweeps spread a slip's shear a few edges only; a sparse solve when positions should settle globally */
void qf_refit(qaws_hm const* dm, mf_level* L, double rho, int const* T, unsigned char const* K, unsigned int sweeps);

#endif /* QAWS_QUAD_FLOW_H */
