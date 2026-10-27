#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_TETRAHEDRALIZE_H_F32
#define QAWS_TETRAHEDRALIZE_H_F32
#define QAWS_TETRAHEDRALIZE_H_BODY
#endif
#else
#ifndef QAWS_TETRAHEDRALIZE_H_F64
#define QAWS_TETRAHEDRALIZE_H_F64
#define QAWS_TETRAHEDRALIZE_H_BODY
#endif
#endif
#ifdef QAWS_TETRAHEDRALIZE_H_BODY
#undef QAWS_TETRAHEDRALIZE_H_BODY
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
 *   STUFFING  isosurface stuffing (Labelle and Shewchuk 2007): a
 *             body-centred cubic lattice of spacing target_length over the
 *             border's box, the signed distance to the border at its
 *             vertices (negative where the winding number passes 1/2: a
 *             soup, holes or overlaps do; the border need not be closed),
 *             each lattice edge crossing the surface cut where the
 *             distance vanishes, a vertex closer to a cut than their alphas
 *             (0.24999 of a cube edge, 0.41189 of a diagonal) warped onto
 *             it, then each lattice cell's inside: whole, a corner
 *             tetrahedron or a prism in three (cut through each quad's
 *             smallest vertex, so neighbours conform). The border is the
 *             surface's approximation by the cuts and warped vertices, not
 *             the input's triangles; the inside the lattice's regular
 *             cells. radius_edge, border_only and keep_border unused; the
 *             report's moves the vertices.
 *   ENVELOPE  reserved (UNSUPPORTED_OPERATION): any triangle soup within
 *             an envelope (fTetWild).
 *
 * The border: a closed qaws_mesh_3d turned outward (polygons fanned into
 * triangles), not self-intersecting. INVALID_TOPOLOGY when it is not
 * closed or its inside cannot be told from its outside; NOT_CONVERGED when
 * the recovery stops short (sharp angles between border triangles);
 * INVALID_ARGUMENT for vertices closer than the lattice step.
 *
 * keep_border: the points the recovery added to the border moved inside
 * after (TetGen's Steiner point suppression): each one's border triangles
 * merged (its link re-triangulated in its plane, or in each of two planes
 * along a crease), its star coned from the point inside where the worst
 * cell is best; then the result improved (qaws_volume_mesh_improve). The
 * border is then the input's surface triangulated by its own vertices, but
 * for points that could not be moved (a star no point inside sees whole).
 * Coplanar border triangles may come back with other diagonals.
 *
 * The report: splits, the border's added vertices; collapses, those moved
 * inside again (keep_border); flips, its diagonals flipped; moves, the
 * vertices added inside; iterations, the recovery's passes; quality, the
 * result's border.
 */

#ifndef QAWS_ENUM_QAWS_TETRAHEDRALIZE_METHOD
#define QAWS_ENUM_QAWS_TETRAHEDRALIZE_METHOD
typedef enum qaws_tetrahedralize_method
{
	QAWS_TETS_DELAUNAY = 0,
	QAWS_TETS_ENVELOPE,
	QAWS_TETS_STUFFING
} qaws_tetrahedralize_method;
#endif

typedef struct qaws_tetrahedralize_desc
{
	qaws_tetrahedralize_method method;
	qaws_scalar target_length;          /* inside; 0: the border's mean edge */
	qaws_scalar radius_edge;            /* refinement bound (circumradius over shortest edge); 0: 2; negative: no refinement */
	int border_only;                    /* 1: no vertex added inside (the border's vertices' cells, flipped) */
	int keep_border;                    /* 1: the recovery's border points moved inside where they can be, the border's triangles merged back */
} qaws_tetrahedralize_desc;

qaws_status qaws_volume_mesh_tetrahedralize(qaws_volume_mesh** out_mesh, qaws_mesh_report* out_report, qaws_mesh_3d const* boundary,
	qaws_tetrahedralize_desc const* desc);

#endif /* QAWS_TETRAHEDRALIZE_H */
