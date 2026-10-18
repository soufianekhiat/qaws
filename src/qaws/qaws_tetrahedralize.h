#ifndef QAWS_TETRAHEDRALIZE_H
#define QAWS_TETRAHEDRALIZE_H

#include "qaws_volume_mesh.h"

/*
 * Tetrahedra filling the inside of a closed surface mesh.
 *
 *   DELAUNAY  the Delaunay tetrahedralization of the border's vertices,
 *             computed exactly (the vertices on an integer lattice of 2^30
 *             steps over their box, Bowyer-Watson with exact predicates);
 *             the border recovered in it (an edge missing from it split at
 *             its midpoint, a triangle at its circumcentre, coplanar
 *             neighbours whose shared diagonal the Delaunay choice turned
 *             flipped to match); the cells inside kept; then the points of
 *             a body-centred cubic lattice inside, and the circumcentres of
 *             cells whose circumradius over shortest edge passes the bound,
 *             inserted with the border as walls; then 3-2 and 2-3 flips and
 *             the inside vertices' smoothing, each kept when it raises the
 *             worst dihedral angle. The border's triangles are kept, split
 *             only where the recovery needed it, so the cells by the border
 *             are as long as its edges. The default.
 *   ENVELOPE, STUFFING  reserved (UNSUPPORTED_OPERATION): any triangle soup
 *             within an envelope (fTetWild); isosurface stuffing.
 *
 * The border: a closed qaws_mesh_3d turned outward (polygons fanned into
 * triangles), not self-intersecting. INVALID_TOPOLOGY when it is not
 * closed or its inside cannot be told from its outside; NOT_CONVERGED when
 * the recovery stops short (sharp angles between border triangles);
 * INVALID_ARGUMENT for vertices closer than the lattice step.
 *
 * The report: splits, the border's added vertices; flips, its diagonals
 * flipped; moves, the vertices added inside; iterations, the recovery's
 * passes; quality, the result's border.
 */

typedef enum qaws_tetrahedralize_method
{
	QAWS_TETS_DELAUNAY = 0,
	QAWS_TETS_ENVELOPE,
	QAWS_TETS_STUFFING
} qaws_tetrahedralize_method;

typedef struct qaws_tetrahedralize_desc
{
	qaws_tetrahedralize_method method;
	qaws_scalar target_length;          /* inside; 0: the border's mean edge */
	qaws_scalar radius_edge;            /* refinement bound (circumradius over shortest edge); 0: 2; negative: no refinement */
	int border_only;                    /* 1: no vertex added inside (the border's vertices' cells, flipped) */
} qaws_tetrahedralize_desc;

qaws_status qaws_volume_mesh_tetrahedralize(qaws_mesh_3d const* boundary, qaws_tetrahedralize_desc const* desc, qaws_volume_mesh** out_mesh,
	qaws_mesh_report* out_report);

#endif /* QAWS_TETRAHEDRALIZE_H */
