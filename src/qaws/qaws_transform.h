#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_TRANSFORM_H_F32
#define QAWS_TRANSFORM_H_F32
#define QAWS_TRANSFORM_H_BODY
#endif
#else
#ifndef QAWS_TRANSFORM_H_F64
#define QAWS_TRANSFORM_H_F64
#define QAWS_TRANSFORM_H_BODY
#endif
#endif
#ifdef QAWS_TRANSFORM_H_BODY
#undef QAWS_TRANSFORM_H_BODY
#define QAWS_TRANSFORM_H

#include "qaws_types.h"
#include "qaws_status.h"

/*
 * A placement of an object built in its own frame: world = translation + R (scale * x), the scale per axis of the
 * object's frame first, then the rotation, then the translation. Kept as its three parts, never as a matrix, so no
 * shear can arise: a scaled box of the object's frame stays a box once the rotation is undone, and queries on a
 * placed object map into its frame without rebuilding what was built there.
 *
 *   rotation     a quaternion (w, x, y, z), normalized where used; all zero: none
 *   scale        per axis; a component 0 counts as 1; negative mirrors (the object's mirror image, its inside kept
 *                inside: a placed solid stays a solid)
 *   translation  added last
 *
 * A zeroed qaws_transform is the identity. In 2D (x, y) the rotation must turn about z (its x and y components 0) and
 * the third components of scale and translation are ignored.
 */
typedef struct qaws_transform
{
	qaws_scalar rotation[4];
	qaws_scalar scale[3];
	qaws_scalar translation[3];
} qaws_transform;

/* count points (3 coordinates each) from the object's frame to the world (out may be points): INVALID_ARGUMENT for a
   NULL transform or points with a count, a value not finite */
qaws_status qaws_transform_points(qaws_scalar* out_points, qaws_transform const* transform, qaws_scalar const* points,
	unsigned int count);

/* the same back, from the world to the object's frame */
qaws_status qaws_transform_points_inverse(qaws_scalar* out_points, qaws_transform const* transform, qaws_scalar const* points,
	unsigned int count);

#endif /* QAWS_TRANSFORM_H */
