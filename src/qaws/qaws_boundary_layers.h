#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_BOUNDARY_LAYERS_H_F32
#define QAWS_BOUNDARY_LAYERS_H_F32
#define QAWS_BOUNDARY_LAYERS_H_BODY
#endif
#else
#ifndef QAWS_BOUNDARY_LAYERS_H_F64
#define QAWS_BOUNDARY_LAYERS_H_F64
#define QAWS_BOUNDARY_LAYERS_H_BODY
#endif
#endif

#ifdef QAWS_BOUNDARY_LAYERS_H_BODY
#undef QAWS_BOUNDARY_LAYERS_H_BODY
#define QAWS_BOUNDARY_LAYERS_H

#include "qaws_volume_mesh.h"

/*
 * Boundary layers: thin cells grown from a surface, in layers whose
 * thickness grows by a ratio (the first layer_count heights first_height,
 * first_height growth, ...), as flow solvers want them along walls. Each
 * triangle of the surface becomes a stack of prisms, each quad a stack of
 * hexahedra; a vertex moves along one direction through every layer.
 *
 * The direction at a vertex is the mean of its faces' normals weighted by
 * their angles there, then turned until every face round it sees it (a
 * corner of a box gets its diagonal), the height scaled so the layers stay
 * as thick across the faces as asked (at most twice). Where the layers
 * would run into the surface (along the direction, the distance to the
 * nearest face not round the vertex) or into each other (neighbours'
 * directions converging, in a concave region), the vertex's whole stack is
 * thinned, at most 0.4 of the distance, and neighbours' thicknesses kept
 * within 1.25 of each other. A cell turned over (a corner's Jacobian not
 * positive) halves its vertices' stacks until none is.
 *
 * The surface: triangles and quads (INVALID_ARGUMENT for larger faces),
 * oriented (a closed one turned outward grows inward by default, outward
 * with desc->outward); open surfaces grow too. INVALID_ARGUMENT where no
 * direction is seen by every face round a vertex (a spike sharper than a
 * right angle's cone).
 *
 * The result: (layer_count + 1) copies of the surface's vertices, layer k
 * after layer k - 1 (vertex v of layer k is k * vertex_count + v), the
 * cells layer by layer in the surface's face order. Its border is the
 * surface, the last layer's offset copy and, on an open surface, the
 * stacks' sides.
 *
 * The report: iterations, the halving rounds; moves, the vertices thinned
 * below the asked thickness; quality, the result's border.
 */

typedef struct qaws_boundary_layers_desc
{
	unsigned int layer_count;           /* at least 1 */
	qaws_scalar first_height;           /* the first layer's thickness, > 0 */
	qaws_scalar growth;                 /* each layer this times the one before; 0: 1.2; else at least 1 */
	int outward;                        /* 0: against the faces' normals (inside a closed surface turned outward); 1: along them */
} qaws_boundary_layers_desc;

qaws_status qaws_volume_mesh_grow_layers(qaws_volume_mesh** out_mesh, qaws_mesh_report* out_report, qaws_mesh_3d const* surface,
	qaws_boundary_layers_desc const* desc);

#endif /* QAWS_BOUNDARY_LAYERS_H */
