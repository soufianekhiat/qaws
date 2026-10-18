#ifndef QAWS_HEXAHEDRALIZE_H
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
 *   POLYCUBE, FRAME_FIELD  reserved (UNSUPPORTED_OPERATION): a polycube map;
 *                a frame field's integer-grid map (qaws_volume_mesh_frame_field
 *                computes the field).
 *
 * The border: a closed qaws_mesh_3d turned outward. INVALID_TOPOLOGY when it
 * is not closed; INVALID_ARGUMENT when no cube fits inside.
 *
 * The report: iterations, the smoothing sweeps; moves, the vertices moved;
 * quality, the result's border.
 */

typedef enum qaws_hexahedralize_method
{
	QAWS_HEXES_GRID = 0,
	QAWS_HEXES_POLYCUBE,
	QAWS_HEXES_FRAME_FIELD
} qaws_hexahedralize_method;

typedef struct qaws_hexahedralize_desc
{
	qaws_hexahedralize_method method;
	qaws_scalar cell_size;              /* the grid's step; 0: the border's mean edge */
} qaws_hexahedralize_desc;

qaws_status qaws_volume_mesh_hexahedralize(qaws_mesh_3d const* boundary, qaws_hexahedralize_desc const* desc, qaws_volume_mesh** out_mesh,
	qaws_mesh_report* out_report);

#endif /* QAWS_HEXAHEDRALIZE_H */
