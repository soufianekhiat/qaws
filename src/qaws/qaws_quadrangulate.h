#ifndef QAWS_QUADRANGULATE_H
#define QAWS_QUADRANGULATE_H

#include "qaws_mesh_reference.h"

/*
 * Quadrangulation: a quad mesh of a surface, its edges near a target
 * length, its vertices on the reference.
 *
 *   INSTANT  Instant Meshes (Jakob, Tarini, Panozzo, Sorkine-Hornung 2015)
 *            on the reference: the input remeshed densely first (a third of
 *            the target), then a cross field (the fewest singularities of
 *            Instant Meshes' own, free or drawn to the principal curvature
 *            directions, and the globally smoothest field of Knoppel et al.
 *            2013, free or drawn alike, relaxed at the lattice's scale
 *            where lines hold it: a torus comes out as a grid along its
 *            parallels, a can's caps as O-grids) and a position
 *            field on a lattice of the target spacing (multiresolution,
 *            along feature edges and the border, lattice points pinned at
 *            feature corners). The lattice steps between neighbours are made
 *            to agree round every triangle by a min-cost flow (QuadriFlow:
 *            Huang, Zhou, Niessner, Shewchuk, Guibas 2018), each cone of the
 *            field put on a lattice point (and on a vertex where that is
 *            cheap), the folds the changes leave undone, so irregular
 *            vertices stay at the field's singularities; then the mesh read
 *            off them: vertices of one lattice point merged, the triangles'
 *            half cells paired into quads, holes filled, doublets removed.
 *            Quad-dominant; `pure` makes all quads at the target spacing, or
 *            where that leaves another face, from fields at twice the target
 *            with every face split in four. Relaxation on the
 *            reference follows, creases and the border along themselves;
 *            then finishing: valence pairs cancelled by local quad
 *            operations (edge rotations, diagonal collapses, doublets)
 *            against a feature-aware ideal valence (a box's corners keep
 *            3), and the quads relaxed by their angles on the reference,
 *            no move lowering the scaled Jacobian round it.
 *            The default.
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
