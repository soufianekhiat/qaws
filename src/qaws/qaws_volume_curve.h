#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_VOLUME_CURVE_H_F32
#define QAWS_VOLUME_CURVE_H_F32
#define QAWS_VOLUME_CURVE_H_BODY
#endif
#else
#ifndef QAWS_VOLUME_CURVE_H_F64
#define QAWS_VOLUME_CURVE_H_F64
#define QAWS_VOLUME_CURVE_H_BODY
#endif
#endif
#ifdef QAWS_VOLUME_CURVE_H_BODY
#undef QAWS_VOLUME_CURVE_H_BODY
#define QAWS_VOLUME_CURVE_H

#include "qaws_volume_mesh.h"
#include "qaws_mesh_reference.h"

/*
 * Curved volume elements: a tetrahedral or hexahedral mesh's quadratic
 * nodes, the border's on a mesh reference, every element's Jacobian kept
 * positive.
 *
 * Elements (VTK's quadratic numbering):
 *   tetrahedron  10 nodes: the corners 0 1 2 3, then the edges' middles
 *                (0 1) (1 2) (2 0) (0 3) (1 3) (2 3);
 *   hexahedron   27 nodes (triquadratic): the corners 0 to 7, the edges'
 *                (0 1) (1 2) (2 3) (3 0) (4 5) (5 6) (6 7) (7 4) (0 4) (1 5)
 *                (2 6) (3 7), the faces' (0 3 7 4) (1 2 6 5) (0 1 5 4)
 *                (3 2 6 7) (0 1 2 3) (4 5 6 7), the centre.
 * The nodes: the mesh's vertices, then one per edge (edges ordered by their
 * smaller vertex index, then the larger), then one per hexahedron's face
 * (faces ordered by their sorted vertices), then one per hexahedron.
 *
 * Placement: every node first at the straight element's point (an edge's
 * middle, a face's or a cell's centre). On the border (facets with no cell
 * across), an edge's node goes to the reference's closest point to it; a
 * face's node to the closest point to its Coons point (the bilinear blend of
 * its curved edges). The inside nodes take the border's displacement,
 * spread through the elements (Jacobi averaging over the elements' nodes,
 * the vertices' displacement zero: it fades inward). Then each element's
 * Jacobian is checked relative to the straight element's: a tetrahedron's
 * bounded from below over the whole element (its determinant, cubic, in the
 * Bernstein basis: the least of its 20 coefficients; Johnen, Remacle,
 * Geuzaine 2013), so a tetrahedron that passes is certified valid; a
 * hexahedron's sampled at its 27 nodes. An element below min_jacobian has
 * its nodes' displacements halved, until it passes (straight at worst,
 * valid when the straight element is).
 * shortcut: a hexahedron's Jacobian is sampled, not bounded; its
 * triquintic determinant's Bernstein coefficients (216) would certify it.
 *
 * The vertices stay where they are: the reference should pass through the
 * border's vertices (a reference made from the volume's border, or a finer
 * mesh of the same surface).
 */

typedef struct qaws_volume_curve_desc
{
	qaws_scalar min_jacobian;           /* the smallest Jacobian an element keeps, relative to its straight self's; 0: 0.1 */
	unsigned int smoothing;             /* Jacobi sweeps spreading the border's displacement inside; 0: 20 */
} qaws_volume_curve_desc;

typedef struct qaws_volume_curve_report
{
	unsigned int border_nodes;          /* nodes placed on the reference */
	unsigned int relaxed;               /* elements whose nodes were pulled back toward straight */
	qaws_scalar min_jacobian;           /* the smallest sampled Jacobian, relative to the straight elements' */
	qaws_scalar max_sagitta;            /* the border nodes' largest distance to the reference, straight */
	qaws_scalar max_distance;           /* and placed (0 unless relaxed) */
} qaws_volume_curve_report;

/* The counts to allocate: nodes, and the elements' node lists' total
   length. UNSUPPORTED_OPERATION for prisms and pyramids. */
qaws_status qaws_volume_mesh_quadratic_node_counts(unsigned int* out_node_count, unsigned int* out_element_node_count, qaws_volume_mesh const* mesh);

/* out_positions: 3 per node; out_element_nodes: per cell in turn its 10 or
   27 nodes. reference NULL: the straight nodes. report and desc may be NULL. */
qaws_status qaws_volume_mesh_quadratic_nodes(qaws_scalar* out_positions, unsigned int* out_element_nodes, qaws_volume_curve_report* report,
	qaws_volume_mesh const* mesh, qaws_mesh_reference const* reference, qaws_volume_curve_desc const* desc);

#endif /* QAWS_VOLUME_CURVE_H */
