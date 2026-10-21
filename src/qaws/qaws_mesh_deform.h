#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_MESH_DEFORM_H_F32
#define QAWS_MESH_DEFORM_H_F32
#define QAWS_MESH_DEFORM_H_BODY
#endif
#else
#ifndef QAWS_MESH_DEFORM_H_F64
#define QAWS_MESH_DEFORM_H_F64
#define QAWS_MESH_DEFORM_H_BODY
#endif
#endif
#ifdef QAWS_MESH_DEFORM_H_BODY
#undef QAWS_MESH_DEFORM_H_BODY
#define QAWS_MESH_DEFORM_H

#include "qaws_mesh_smooth.h"

/*
 * Shape deformation by handles, in place: the handle vertices go to their
 * targets, the other vertices of the region follow, those outside it stay.
 * The mesh's positions are the rest shape.
 *
 *   ARAP         as rigid as possible (Sorkine, Alexa 2007): each vertex's
 *                ring turns as one, the rotations and the positions found
 *                in turn (local-global); details turn with the shape.
 *                Started from the BIHARMONIC solution; large turns take
 *                tens of rounds.
 *   LAPLACIAN    the rest's Laplacian coordinates kept (Sorkine et al.
 *                2004): linear, details kept but not turned.
 *   BIHARMONIC   the displacement as smooth as possible (Botsch, Kobbelt
 *                2004): linear, the region bends like a thin plate.
 *
 * The Laplacian is the rest shape's cotangent one (or uniform). The region
 * needs a held vertex or a handle in each of its parts.
 */

#ifndef QAWS_ENUM_QAWS_DEFORM_METHOD
#define QAWS_ENUM_QAWS_DEFORM_METHOD
typedef enum qaws_deform_method
{
	QAWS_DEFORM_ARAP = 0,
	QAWS_DEFORM_LAPLACIAN,
	QAWS_DEFORM_BIHARMONIC
} qaws_deform_method;
#endif

typedef struct qaws_deform_desc
{
	qaws_deform_method method;
	qaws_smooth_weights weights;
	unsigned int iterations;            /* ARAP's local-global rounds at most; 0: 200 (they stop once no vertex moves by 1e-6
	                                       of the box's diagonal) */
	unsigned char const* region;        /* per vertex, nonzero: may move; NULL: all of them */
} qaws_deform_desc;

/* handles (handle_count vertices) put at targets (3 each); desc NULL: ARAP. INVALID_ARGUMENT for a handle out of range */
qaws_status qaws_mesh_3d_deform(qaws_mesh_3d* mesh, unsigned int const* handles, qaws_scalar const* targets, unsigned int handle_count,
	qaws_deform_desc const* desc);

#endif /* QAWS_MESH_DEFORM_H */
