#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_MESH_SUBDIVIDE_H_F32
#define QAWS_MESH_SUBDIVIDE_H_F32
#define QAWS_MESH_SUBDIVIDE_H_BODY
#endif
#else
#ifndef QAWS_MESH_SUBDIVIDE_H_F64
#define QAWS_MESH_SUBDIVIDE_H_F64
#define QAWS_MESH_SUBDIVIDE_H_BODY
#endif
#endif
#ifdef QAWS_MESH_SUBDIVIDE_H_BODY
#undef QAWS_MESH_SUBDIVIDE_H_BODY
#define QAWS_MESH_SUBDIVIDE_H

#include "qaws_mesh.h"

/*
 * Subdivision of meshes: Loop (triangles) and Catmull-Clark (any polygons,
 * quads out) as operations on a qaws_mesh_3d, the result a mesh again.
 *
 * Each level splits every face (Loop: a triangle into four; Catmull-Clark: an
 * n-gon into n quads through its centroid) and places the vertices by the
 * scheme's rules (Loop 1987; Catmull and Clark 1978), with the creases of
 * Hoppe et al. 1994 and the semi-sharp creases of DeRose, Kass and Truong
 * 1998:
 *   - the mesh's feature edges are creases: their new points are midpoints,
 *     and a vertex on two of them moves along the crease, as a cubic
 *     B-spline curve (1/8, 6/8, 1/8); on three or more it stays;
 *   - border edges are creases always, and a border vertex of one face (the
 *     corner of a sheet) stays;
 *   - feature corners stay while a crease at them is sharp (but the end of
 *     a lone feature line inside the surface, a dart, is smooth);
 *   - a sharpness per feature edge makes it semi-sharp: sharp for that many
 *     levels, then smooth; a fraction blends the sharp and smooth rules at
 *     the last level (a vertex by the mean sharpness of its creases).
 *
 * The result carries the features (each crease edge's halves, while sharp;
 * the corners) and every attribute by its rule:
 *   per vertex  LINEAR and NORMALIZED by the scheme's own rules, as the
 *               positions (NORMALIZED rescaled to unit length); NEAREST
 *               from the vertex itself, an edge's first end, a face's first
 *               corner;
 *   per corner  linearly within each face (an edge's new corners the mean of
 *               its two corners in the face, the centroid's the mean of all):
 *               UVs keep their seams, follow the faces, not the smoothing;
 *               NEAREST from the corner new ones come from first;
 *   per face    the parent face's value.
 * The result has no footpoints. New vertices follow the old ones: the
 * input's vertices keep their indices, then one per edge (the edges ordered
 * by their smaller vertex index, then the larger), then (Catmull-Clark) one
 * per face. Faces follow their parent's order: Loop's four children are the
 * corners' triangles then the middle one, Catmull-Clark's are the corners'
 * quads in order.
 */

#ifndef QAWS_ENUM_QAWS_MESH_SUBDIVIDE_SCHEME
#define QAWS_ENUM_QAWS_MESH_SUBDIVIDE_SCHEME
typedef enum qaws_mesh_subdivide_scheme
{
	QAWS_MESH_SUBDIVIDE_CATMULL_CLARK = 0,
	QAWS_MESH_SUBDIVIDE_LOOP
} qaws_mesh_subdivide_scheme;
#endif

typedef struct qaws_mesh_subdivide_desc
{
	qaws_mesh_subdivide_scheme scheme;
	unsigned int levels;                /* 0: 1 */
	qaws_scalar const* sharpness;       /* per feature edge (qaws_mesh_3d_get_feature_edges order), >= 0; NULL: infinitely sharp */
} qaws_mesh_subdivide_desc;

/* INVALID_ARGUMENT for a mesh without faces, an unknown scheme, Loop on a face
   that is not a triangle, a negative sharpness; ALLOCATION_FAILURE when the
   result's size does not fit. */
qaws_status qaws_mesh_3d_subdivide(qaws_mesh_3d** out_mesh, qaws_mesh_3d const* mesh, qaws_mesh_subdivide_desc const* desc);

#endif /* QAWS_MESH_SUBDIVIDE_H */
