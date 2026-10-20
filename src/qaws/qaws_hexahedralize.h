#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_HEXAHEDRALIZE_H_F32
#define QAWS_HEXAHEDRALIZE_H_F32
#define QAWS_HEXAHEDRALIZE_H_BODY
#endif
#else
#ifndef QAWS_HEXAHEDRALIZE_H_F64
#define QAWS_HEXAHEDRALIZE_H_F64
#define QAWS_HEXAHEDRALIZE_H_BODY
#endif
#endif
#ifdef QAWS_HEXAHEDRALIZE_H_BODY
#undef QAWS_HEXAHEDRALIZE_H_BODY
#define QAWS_HEXAHEDRALIZE_H

#include "qaws_volume_mesh.h"

/*
 * Hexahedra filling the inside of a closed surface mesh.
 *
 *   GRID         a grid of cubes over the box (Marechal 2009; Gao, Shen,
 *                Panozzo 2019): the cubes whose corners are all inside kept
 *                (the inside told by the parity of a ray's crossings), the
 *                set made manifold (cubes meeting only at an edge or a
 *                vertex dropped), then a layer of hexahedra padding it to
 *                the surface: each of its border's vertices projected onto
 *                the surface's PN reference (a vertex where the staircase
 *                turns once onto a crease nearby, where it turns twice onto
 *                a corner, so a box keeps its edges), each of its border's
 *                quads extruded to its projection; then every vertex but the
 *                surface's smoothed, each move kept when it raises the worst
 *                scaled Jacobian round it. Never fails on a closed surface
 *                the cubes fit in; the elements by the surface are the
 *                poorest (a sharp concave feature can turn some over:
 *                qaws_volume_mesh_measure_quality tells). The default.
 *   FRAME_FIELD  hexahedra following a frame field (an integer-grid map,
 *                CubeCover's idea): the inside tetrahedralized (DELAUNAY, half
 *                the cell size), its frame field (qaws_volume_mesh_frame_field)
 *                combed into one gauge, then a map u fitted by three Poisson
 *                solves so that its gradient follows the frames at the cell
 *                size, each border patch (facets whose normal follows one
 *                axis) held on a whole-number plane of that axis; the
 *                hexahedra are the unit cubes of the integer lattice in the
 *                map's image, their vertices mapped back. A rotated box comes
 *                out as a grid aligned with it, a bent bar as a bent grid.
 *                Fields with singular edges (a sphere, a round torus) or that
 *                turn round a hole: UNSUPPORTED_OPERATION for now.
 *   POLYCUBE     reserved (UNSUPPORTED_OPERATION): a polycube map.
 *
 * The border: a closed qaws_mesh_3d turned outward. INVALID_TOPOLOGY when it
 * is not closed; INVALID_ARGUMENT when no cube fits inside.
 *
 * The report: iterations, the smoothing sweeps; moves, the vertices moved;
 * quality, the result's border.
 */

#ifndef QAWS_ENUM_QAWS_HEXAHEDRALIZE_METHOD
#define QAWS_ENUM_QAWS_HEXAHEDRALIZE_METHOD
typedef enum qaws_hexahedralize_method
{
	QAWS_HEXES_GRID = 0,
	QAWS_HEXES_POLYCUBE,
	QAWS_HEXES_FRAME_FIELD
} qaws_hexahedralize_method;
#endif

typedef struct qaws_hexahedralize_desc
{
	qaws_hexahedralize_method method;
	qaws_scalar cell_size;              /* the grid's step; 0: the border's mean edge */
} qaws_hexahedralize_desc;

qaws_status qaws_volume_mesh_hexahedralize(qaws_volume_mesh** out_mesh, qaws_mesh_report* out_report, qaws_mesh_3d const* boundary,
	qaws_hexahedralize_desc const* desc);

#endif /* QAWS_HEXAHEDRALIZE_H */
