#ifndef QAWS_QUADRANGULATE_H
#define QAWS_QUADRANGULATE_H

#include "qaws_mesh_reference.h"
#include "qaws_mesh_field.h"

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
 *            3), and the quads relaxed by their angles and lengths on the
 *            reference, no move lowering the scaled Jacobian round it or
 *            squeezing its quads toward triangles.
 *            The default.
 *   GLOBAL   an integer-grid map (Bommes, Zimmer, Kobbelt 2009; Bommes et
 *            al. 2013) on INSTANT's dense mesh and fields: the lattice
 *            steps INSTANT's flow made consistent are the map's integers
 *            (cones on lattice points, feature lines and the border on
 *            whole isolines, corners on lattice points), each vertex's place
 *            in its cell solved globally (one sparse LDL^T, the energy of
 *            Mixed-Integer Quadrangulation), folds left undone vertex by
 *            vertex; the quads read off the map as QEx does (Ebke et al.
 *            2013): its lattice points the vertices, its unit isolines
 *            traced across the triangles the edges. All quads, irregular
 *            vertices only at the cones; the same finishing as INSTANT.
 *            When the steps leave no map, the map keeps a fold, or its
 *            isolines do not close into quads, INSTANT's mesh is returned:
 *            the report's `iterations` is 1 when the map made the mesh, 0
 *            when it fell back.
 *   PAIRING, SPLIT, LAYOUT  reserved (UNSUPPORTED_OPERATION today):
 *            triangles paired by matching, all-quad splitting, patch
 *            layouts.
 *
 * As remeshing: footpoints kept, the reference mesh's attributes
 * transferred, the input unchanged. The report counts the singularities
 * (`quality.irregular`) and the faces by size.
 *
 * Guidance (INSTANT; all of it indexes the input's vertices):
 *   size        density paint: the quads' edges follow the field (as
 *               remeshing's, over target_length): the dense mesh graded
 *               by it, the lattice's spacing per vertex (Instant Meshes'
 *               position field at a local scale); where the spacing
 *               changes the lattice needs extra rows, which the flow lets
 *               end there (each a pentagon or a 3-5 pair).
 *   guide_directions, guide_weights  as the direction field's guidance
 *               (qaws_mesh_field.h), carried to the dense mesh through its
 *               footpoints (a pull fades across the triangles round its
 *               vertex, a hold covers the dense vertices nearest it); every
 *               candidate field honours it. Needs the reference made from
 *               the input (INVALID_ARGUMENT otherwise).
 *   guide_edges edge loops: pairs of input vertices, segments the quads'
 *               edges run along without their being creases (Instant
 *               Meshes' edge strokes): the dense vertices within half a
 *               dense edge of one held along it, their lattice points on
 *               one lattice line, its quad vertices placed on the segments
 *               and kept there by the relaxation and the finishing.
 *   pinned      input vertices that become quad vertices where they are
 *               (the dense remesh keeps them, as LOCAL does, the lattice
 *               holds a point there); they come out as feature corners.
 *   cone_vertices, cone_indices  singularities placed (INSTANT and
 *               GLOBAL): the dense remesh keeps each such vertex, the
 *               field is the smoothest with its cones there (the
 *               direction field's SMOOTHEST with cones, then relaxed where
 *               lines hold it; Instant Meshes' own fields are not tried),
 *               the lattice points round it set on it. With a border the
 *               flow keeps the steps round the cone (the border takes what
 *               cannot close), so the cone's quad vertex is the input
 *               vertex; on a closed surface the cones must take the
 *               lattice's mismatch and may move a step, and the nearest
 *               quad vertex of valence 4 - index is moved onto the cone.
 *               Either way a quad vertex of valence 4 - index sits exactly
 *               there (kept by the relaxation and the finishing, not
 *               flagged a corner). Poincare-Hopf as the direction
 *               field's: on a closed surface without feature lines the
 *               indices sum to 4 chi (8 on a sphere: a cube-sphere from
 *               eight +1 cones), else INVALID_ARGUMENT; with lines or a
 *               border what they leave over is placed freely.
 *   attractors, attractor_strengths, attractor_radius  the direction
 *               field's attractors on every candidate field (LOCAL's
 *               smoothing and SMOOTHEST's edge weights): free singularities
 *               drawn to or pushed from points.
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
	qaws_size_field const* size;            /* density: the quads' edges, over target_length; on `reference` when sampled on one */
	qaws_scalar const* guide_directions;    /* 3 per input vertex (NULL: none) */
	qaws_scalar const* guide_weights;       /* per input vertex: > 0 a pull, < 0 or infinite held, 0 none; NULL: 1 where given */
	unsigned int const* guide_edges;        /* pairs of input vertices: segments quad edges run along */
	unsigned int guide_edge_count;
	unsigned int const* pinned;             /* input vertices kept as quad vertices */
	unsigned int pinned_count;
	unsigned int const* cone_vertices;      /* input vertices made irregular quad vertices */
	int const* cone_indices;                /* per cone in quarter turns: +1 valence 3, -1 valence 5 (+2 valence 2, ...) */
	unsigned int cone_count;
	qaws_scalar const* attractors;          /* 3 per point: the field's singularities drawn to it or pushed off */
	qaws_scalar const* attractor_strengths; /* per point: > 0 draws, < 0 pushes (the direction field's) */
	unsigned int attractor_count;
	qaws_scalar attractor_radius;           /* 0: a tenth of the input's bounding diagonal */
} qaws_quadrangulate_desc;

qaws_status qaws_mesh_quadrangulate(qaws_mesh_3d** out_mesh, qaws_mesh_report* out_report, qaws_mesh_3d const* mesh,
	qaws_quadrangulate_desc const* desc);

#endif /* QAWS_QUADRANGULATE_H */
