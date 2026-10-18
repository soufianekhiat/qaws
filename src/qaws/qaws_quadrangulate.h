#ifndef QAWS_QUADRANGULATE_H
#define QAWS_QUADRANGULATE_H

#include "qaws_mesh_reference.h"

/*
 * Quadrangulation: a quad mesh of a surface, its edges near a target
 * length, its vertices on the reference.
 *
 *   INSTANT  Instant Meshes (Jakob, Tarini, Panozzo, Sorkine-Hornung 2015)
 *            on the reference: the input remeshed densely first (a third of
 *            the target), then a cross field and a position field on a
 *            lattice of the target spacing (multiresolution, along feature
 *            edges and the border, lattice points pinned at feature
 *            corners), then the mesh read off the fields: vertices whose
 *            lattice points coincide merged, neighbours one lattice step
 *            apart joined, faces traced. Quad-dominant (some triangles and
 *            pentagons where the field turns); `pure` makes all quads by
 *            splitting every face (the fields at twice the target). A few
 *            relaxation steps on the reference follow. The default.
 *   PAIRING, SPLIT, GLOBAL, LAYOUT  reserved (UNSUPPORTED_OPERATION today):
 *            triangles paired by matching, all-quad splitting, integer-grid
 *            maps, patch layouts.
 *
 * As remeshing: footpoints kept, the reference mesh's attributes
 * transferred, the input unchanged. The report counts the singularities
 * (`quality.irregular`) and the faces by size.
 */

typedef enum qaws_quad_method
{
	QAWS_QUADS_INSTANT = 0,
	QAWS_QUADS_PAIRING,
	QAWS_QUADS_SPLIT,
	QAWS_QUADS_GLOBAL,
	QAWS_QUADS_LAYOUT
} qaws_quad_method;

typedef struct qaws_quadrangulate_desc
{
	qaws_quad_method method;
	qaws_mesh_reference const* reference;   /* NULL: the input's own, PN */
	qaws_scalar target_length;              /* the quads' edges; 0: twice the input's mean edge */
	int pure;                               /* all quads */
	unsigned int sweeps;                    /* field smoothing sweeps per level; 0: 8 */
} qaws_quadrangulate_desc;

qaws_status qaws_mesh_quadrangulate(qaws_mesh_3d const* mesh, qaws_quadrangulate_desc const* desc, qaws_mesh_3d** out_mesh,
	qaws_mesh_report* out_report);

#endif /* QAWS_QUADRANGULATE_H */
