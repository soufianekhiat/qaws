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
 * Curved volume elements: a tetrahedral or hexahedral mesh's quadratic or
 * cubic nodes, the border's on a mesh reference, every element's Jacobian
 * kept positive.
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
 * Cubic elements (the corners first, as above):
 *   tetrahedron  20 nodes: the corners, then each edge's two (in the order
 *                above, at 1/3 and 2/3 from its first corner), then each
 *                face's centre, the faces (0 1 3) (1 2 3) (0 2 3) (0 1 2);
 *   hexahedron   64 nodes (tricubic): the corners, each edge's two, then
 *                each face's four (the faces above, at (a, b) / 3 in the
 *                face's frame, a from its first corner toward its second,
 *                b toward its last, a fastest), then the inside's eight
 *                (s fastest, then t, then u).
 * The cubic nodes: the vertices, then two per edge (from its smaller vertex),
 * then each face's (a hexahedron's four in the frame from its smallest
 * vertex toward its smaller neighbour, a tetrahedron's one), then each
 * hexahedron's eight. A cubic face node on the border starts at its Coons
 * point (a tetrahedron's face: its six curved edge nodes' mean times 3 / 2,
 * less half its corners' mean). Their Jacobians are bounded the same way,
 * a tetrahedron's determinant of degree 6 by its 84 Bernstein coefficients,
 * a hexahedron's of degree 8 in each axis by its 729.
 *
 * Placement: every node first at the straight element's point (an edge's
 * middle, a face's or a cell's centre). On the border (facets with no cell
 * across), an edge's nodes go to equal lengths along the curved edge (its
 * segment's samples projected onto the reference, the node at its fraction
 * of their polyline's length, then projected); a face's node to the
 * closest point to its Coons point (the blend of its curved edges). The
 * inside nodes take the border's displacement, spread through the elements
 * (Jacobi averaging over the elements' nodes, the vertices' displacement
 * zero: it fades inward); then the faces inside the mesh take their edges'
 * Coons blend and a hexahedron's own nodes its boundary's (Gordon and Hall's
 * transfinite blend), so they follow the border's curving. An element whose
 * Jacobian bound is then under min_jacobian is untangled first: its free
 * nodes (inside the mesh) moved by pattern search, a step along an axis kept
 * when it raises the least sampled Jacobian of the elements sharing the
 * node, the step halved when none does, until its bound passes; the border
 * left where it is. Then each element's
 * Jacobian is checked relative to the straight element's: a tetrahedron's
 * bounded from below over the whole element (its determinant, cubic, in the
 * Bernstein basis: the least of its 20 coefficients; Johnen, Remacle,
 * Geuzaine 2013); a hexahedron's, of degree 5 in each axis, by its 216
 * tensor Bernstein coefficients, each over the straight hexahedron's in the
 * same basis, the box halved where the bound is under min_jacobian (three
 * times at most: the bound tightens as the boxes shrink). An element that
 * passes is certified valid. An element below
 * min_jacobian has its nodes' displacements halved, until it passes
 * (straight at worst, valid when the straight element is).
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
	qaws_scalar min_jacobian;           /* the smallest Jacobian (a lower bound), relative to the straight elements' */
	qaws_scalar max_sagitta;            /* the border nodes' largest distance to the reference, straight */
	qaws_scalar max_distance;           /* and placed (0 unless relaxed) */
	unsigned int untangled;             /* elements brought over min_jacobian by moving their free nodes alone */
} qaws_volume_curve_report;

/* The counts to allocate: nodes, and the elements' node lists' total
   length. UNSUPPORTED_OPERATION for prisms and pyramids. */
qaws_status qaws_volume_mesh_quadratic_node_counts(unsigned int* out_node_count, unsigned int* out_element_node_count, qaws_volume_mesh const* mesh);

/* out_positions: 3 per node; out_element_nodes: per cell in turn its 10 or
   27 nodes. reference NULL: the straight nodes. report and desc may be NULL. */
qaws_status qaws_volume_mesh_quadratic_nodes(qaws_scalar* out_positions, unsigned int* out_element_nodes, qaws_volume_curve_report* report,
	qaws_volume_mesh const* mesh, qaws_mesh_reference const* reference, qaws_volume_curve_desc const* desc);

/* The cubic counterparts: 20 nodes per tetrahedron, 64 per hexahedron. */
qaws_status qaws_volume_mesh_cubic_node_counts(unsigned int* out_node_count, unsigned int* out_element_node_count, qaws_volume_mesh const* mesh);
qaws_status qaws_volume_mesh_cubic_nodes(qaws_scalar* out_positions, unsigned int* out_element_nodes, qaws_volume_curve_report* report,
	qaws_volume_mesh const* mesh, qaws_mesh_reference const* reference, qaws_volume_curve_desc const* desc);

#endif /* QAWS_VOLUME_CURVE_H */
