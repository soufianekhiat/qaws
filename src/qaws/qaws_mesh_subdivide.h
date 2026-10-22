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
 * Subdivision of meshes: Loop (triangles), the modified butterfly
 * (triangles, interpolating) and Catmull-Clark (any polygons, quads out) as
 * operations on a qaws_mesh_3d, the result a mesh again.
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
 * limit moves the last level's vertices onto the limit surface (vertex
 * attributes LINEAR and NORMALIZED with them): Loop's limit mask (the
 * vertex 1 - n chi, each neighbour chi = 1 / (3 / (8 beta) + n)),
 * Catmull-Clark's ((n^2 V + 4 sum E + sum F) / (n (n + 5)), E the edge
 * neighbours, F the diagonal ones), a crease's cubic B-spline limit
 * (1, 4, 1) / 6; corners stay. The vertices then lie on the surface the
 * levels converge to (Halstead, Kass, DeRose 1993).
 *
 * limit_normals adds the limit surface's unit normals as the vertex
 * attribute "limit_normal": at a smooth vertex the cross product of the
 * scheme's two tangent masks over its ring (Loop's cos and sin weights;
 * Catmull-Clark's of Halstead, Kass and DeRose, the edge neighbours by
 * A_n = 1 + cos(2 pi / n) + cos(pi / n) sqrt(2 (9 + cos(2 pi / n))) cos,
 * the diagonal ones by the sum of the two cos beside them), exact; at a
 * vertex on a border, a crease or a corner its faces' normals at the limit
 * positions summed, approximate (the creases' tangent masks are not used).
 *
 * The modified butterfly (Dyn, Levin and Gregory 1990; Zorin, Schroeder and
 * Sweldens 1996) splits triangles as Loop does but keeps every vertex where
 * it is, the surface through them. An edge's point between two regular
 * (valence 6) vertices takes the ten points: 1/2 the ends, 1/8 the two
 * opposite vertices, -1/16 the four wings (a wing across a border
 * reflected); next to an extraordinary vertex of valence n its ring, 3/4
 * the vertex and s_j its j-th neighbour from the edge (n = 3: 5/12, -1/12,
 * -1/12; n = 4: 3/8, 0, -1/8, 0; else (1/4 + cos(2 pi j / n) +
 * 1/2 cos(4 pi j / n)) / n), two such ends averaged. A crease or border
 * edge takes the four points along it (-1/16, 9/16, 9/16, -1/16, an end
 * where the crease stops reflected). Cubic polynomials are reproduced on
 * regular parts. limit leaves the vertices (they are on the surface);
 * limit_normals is refused. Next to a border the stencils' vertices beyond
 * it are virtual, reflected across it (Zorin et al.'s rule): planes are
 * kept up to the border, and the border is the four-point curve.
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
	QAWS_MESH_SUBDIVIDE_LOOP,
	QAWS_MESH_SUBDIVIDE_BUTTERFLY
} qaws_mesh_subdivide_scheme;
#endif

typedef struct qaws_mesh_subdivide_desc
{
	qaws_mesh_subdivide_scheme scheme;
	unsigned int levels;                /* 0: 1 */
	qaws_scalar const* sharpness;       /* per feature edge (qaws_mesh_3d_get_feature_edges order), >= 0; NULL: infinitely sharp */
	int limit;                          /* 1: the last level's vertices moved onto the limit surface */
	int limit_normals;                  /* 1: limit, and a vertex attribute "limit_normal" (NORMALIZED, 3): the limit surface's normals */
} qaws_mesh_subdivide_desc;

/* INVALID_ARGUMENT for a mesh without faces, an unknown scheme, Loop or the
   butterfly on a face that is not a triangle, limit_normals with the
   butterfly, a negative sharpness; ALLOCATION_FAILURE when the
   result's size does not fit. */
qaws_status qaws_mesh_3d_subdivide(qaws_mesh_3d** out_mesh, qaws_mesh_3d const* mesh, qaws_mesh_subdivide_desc const* desc);

#endif /* QAWS_MESH_SUBDIVIDE_H */
