#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_MESH_SMOOTH_H_F32
#define QAWS_MESH_SMOOTH_H_F32
#define QAWS_MESH_SMOOTH_H_BODY
#endif
#else
#ifndef QAWS_MESH_SMOOTH_H_F64
#define QAWS_MESH_SMOOTH_H_F64
#define QAWS_MESH_SMOOTH_H_BODY
#endif
#endif
#ifdef QAWS_MESH_SMOOTH_H_BODY
#undef QAWS_MESH_SMOOTH_H_BODY
#define QAWS_MESH_SMOOTH_H

#include "qaws_mesh.h"

/*
 * Smoothing and fairing on a mesh's vertices, in place (the connectivity
 * and attributes are kept; positions move).
 *
 * Smoothing removes noise:
 *   LAPLACIAN    explicit steps toward the neighbours' weighted mean
 *                (the umbrella operator): simple, but the mesh shrinks;
 *   TAUBIN       a step lambda then a step back mu (Taubin 1995): a
 *                low-pass filter that does not shrink;
 *   IMPLICIT     the mean curvature flow by backward Euler (Desbrun,
 *                Meyer, Schroder, Barr 1999): (M + t L) x' = M x, stable
 *                for any step, one solve smooths at the scale sqrt(t);
 *   BILAPLACIAN  the fourth order flow (M + t L M^-1 L) x' = M x: smooths
 *                while keeping more of the shape's curvature.
 * L is the cotangent Laplacian with lumped masses M (or the uniform graph
 * Laplacian with M = I); the explicit steps clamp negative cotangents.
 *
 * What stays: vertices marked held; the border (unless
 * QAWS_SMOOTH_MOVE_BORDER: it is then smoothed along itself, its corners
 * held); the mesh's feature corners. Vertices on feature lines (the mesh's
 * feature edges, qaws_mesh.h) are smoothed along their line only, from
 * their two line neighbours, so creases stay sharp; vertices on more lines
 * are held. With QAWS_SMOOTH_KEEP_VOLUME a closed mesh with nothing held is
 * scaled about its centroid after each step back to its volume.
 *
 * Fairing makes a region as smooth as possible given the rest: the
 * vertices marked free solve L^k x = 0 (L the Laplacian above, its powers
 * through M^-1), the others held. Order 1 is a membrane (the free region
 * is the harmonic, smallest area spanning of its border), 2 a thin plate
 * (tangent to the surface round it: the default, as hole fills and
 * bridges use), 3 continuous in curvature as well. With cotangent weights
 * (which read the free vertices' positions, however disturbed) the system
 * is solved with uniform weights first, then twice with the cotangents of
 * the last solution. The free region needs
 * k rings of held vertices round it to set its border conditions.
 */

#ifndef QAWS_ENUM_QAWS_SMOOTH_METHOD
#define QAWS_ENUM_QAWS_SMOOTH_METHOD
typedef enum qaws_smooth_method
{
	QAWS_SMOOTH_LAPLACIAN = 0,
	QAWS_SMOOTH_TAUBIN,
	QAWS_SMOOTH_IMPLICIT,
	QAWS_SMOOTH_BILAPLACIAN
} qaws_smooth_method;
#endif

#ifndef QAWS_ENUM_QAWS_SMOOTH_WEIGHTS
#define QAWS_ENUM_QAWS_SMOOTH_WEIGHTS
typedef enum qaws_smooth_weights
{
	QAWS_SMOOTH_WEIGHTS_COTANGENT = 0,  /* the surface's own Laplacian: smooths the shape, not the sampling */
	QAWS_SMOOTH_WEIGHTS_UNIFORM         /* the graph's: also evens out the triangles */
} qaws_smooth_weights;
#endif

#define QAWS_SMOOTH_MOVE_BORDER 1u          /* the border smoothed along itself, not held */
#define QAWS_SMOOTH_KEEP_VOLUME 2u          /* a closed mesh with nothing held keeps its volume */

typedef struct qaws_smooth_desc
{
	qaws_smooth_method method;
	qaws_smooth_weights weights;
	unsigned int iterations;            /* 0: 10 for the explicit methods, 1 for the implicit ones */
	qaws_scalar step;                   /* explicit: the fraction of the way to the mean (0: 0.5); IMPLICIT: t over the mean edge
	                                       length squared, BILAPLACIAN: over its fourth power (0: 1) */
	qaws_scalar mu;                     /* TAUBIN's step back (0: -0.53, with step 0.5) */
	unsigned int flags;                 /* QAWS_SMOOTH_* */
	unsigned char const* held;          /* per vertex, nonzero: held; NULL: none */
} qaws_smooth_desc;

/* desc NULL: 10 Laplacian steps of 0.5, cotangent weights */
qaws_status qaws_mesh_3d_smooth(qaws_mesh_3d* mesh, qaws_smooth_desc const* desc);

typedef struct qaws_fair_desc
{
	unsigned int order;                 /* 1, 2 or 3; 0: 2 */
	qaws_smooth_weights weights;
} qaws_fair_desc;

/* free_vertices: per vertex, nonzero where the vertex moves. INVALID_ARGUMENT for an order past 3; a free part touching no
   held vertex is left as it was */
qaws_status qaws_mesh_3d_fair(qaws_mesh_3d* mesh, unsigned char const* free_vertices, qaws_fair_desc const* desc);

#endif /* QAWS_MESH_SMOOTH_H */
