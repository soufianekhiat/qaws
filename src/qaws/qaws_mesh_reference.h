#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_MESH_REFERENCE_H_F32
#define QAWS_MESH_REFERENCE_H_F32
#define QAWS_MESH_REFERENCE_H_BODY
#endif
#else
#ifndef QAWS_MESH_REFERENCE_H_F64
#define QAWS_MESH_REFERENCE_H_F64
#define QAWS_MESH_REFERENCE_H_BODY
#endif
#endif
#ifdef QAWS_MESH_REFERENCE_H_BODY
#undef QAWS_MESH_REFERENCE_H_BODY
#define QAWS_MESH_REFERENCE_H

#include "qaws_mesh.h"

/*
 * References: the fixed surface the meshing operators work on.
 *
 * A triangle mesh with vertex normals implies a smooth surface; the mesh is
 * a sampling of it. A mesh reference is that surface: each triangle (a
 * polygon is fanned from its first vertex) carries a patch through its
 * corners, built edge curves first, so neighbouring patches share their
 * edge curve and the surface has no cracks:
 *
 *   PN       curved PN triangles (Vlachos, Peters, Boyd, Mitchell 2001): a
 *            cubic patch whose edge curves leave each corner along the
 *            chord projected onto the corner's tangent plane; G1 at the
 *            vertices, C0 across edges; a quadratic normal field. Along a
 *            feature edge (a crease) the curve follows the crease line
 *            instead (its tangent from the neighbouring crease vertices,
 *            straight at feature corners), and each side uses its own
 *            normals: creases stay sharp, without cracks. The border is a
 *            crease too (its curves from the border alone). The default.
 *   LINEAR   the flat triangles, normals interpolated (for shading and
 *            curvature only): for meshes that are the truth.
 *   PHONG    Phong tessellation (Boubekeur and Alexa 2008): each point the
 *            corners' tangent-plane projections of the flat point, blended
 *            by its weights, 3/4 against the flat triangle: a quadratic
 *            patch (cheaper to read than PN, rounder than flat: on a sphere
 *            past it by half the flat sag at an edge's middle); C0 across
 *            edges, normals interpolated linearly; a feature edge straight.
 *   NAGATA   Nagata's patch (2005): each edge the quadratic leaving both
 *            ends in their tangent planes, the patch the quadratic through
 *            them (on a sphere an edge's middle on it to leading order); C0
 *            across edges, normals interpolated linearly; a feature edge
 *            straight.
 *   G1       Walton and Meek's patch (1996): PN's boundary curves, the cross
 *            tangents made the same plane from both sides of an edge
 *            (Chiyokura and Kimura), a quartic Gregory interior; the
 *            normal the surface's own, continuous across edges (G1) but
 *            for creases.
 *
 *   SUBDIVISION  the Loop limit surface of the mesh as its control net (Stam's
 *            box splines where regular, local subdivision elsewhere; creases
 *            and the border by Hoppe et al.'s rules): smooth (C2, C1 at
 *            extraordinary vertices), not through the vertices (it shrinks:
 *            a vertex's location is its limit point). UNSUPPORTED_OPERATION
 *            past about 50 triangles round a vertex.
 *
 *   MLS      the moving least squares surface of the vertices and their
 *            normals (APSS, Guennebaud and Gross 2007: a sphere fitted round each point, exact on spheres), each vertex's
 *            kernel mls_radius times its mean edge length: smooth, near the
 *            vertices but not through them, noise and creases smoothed
 *            (larger radii smooth more); a location's point the flat
 *            triangle's moved along the blended normals onto it.
 *
 * Normals: given per vertex, or computed: the faces' normals round each
 * vertex weighted by their angles there, one normal per run of faces
 * between feature edges (split normals at creases).
 *
 * A location (qaws_mesh_location) is a triangle and the barycentric
 * coordinates of its second and third corners. On it:
 *   evaluate  position, the normal field, principal curvatures and
 *             directions (from the normal field's derivative: the shape
 *             operator; positive where the surface bends away from its
 *             normals, as on a sphere with outward normals);
 *   project   a point's closest location; with hints, the closest one
 *             reachable from the hint's triangle across triangles within
 *             twice the best distance found (a sheet nearby is not jumped
 *             to);
 *   walk      the straightest path from a location along a tangent
 *             vector, on the triangles (unfolded across each edge), its
 *             end carried onto the patch: tangential moves that stay on the
 *             surface without projecting;
 *   attributes  the input's attributes at a location, by their rules.
 *
 * The reference copies its mesh: the input may change or go after creation.
 * Thread-safe for concurrent reads.
 */

#ifndef QAWS_ENUM_QAWS_INTERPOLANT
#define QAWS_ENUM_QAWS_INTERPOLANT
typedef enum qaws_interpolant
{
	QAWS_INTERPOLANT_PN = 0,
	QAWS_INTERPOLANT_LINEAR,
	QAWS_INTERPOLANT_PHONG,
	QAWS_INTERPOLANT_NAGATA,
	QAWS_INTERPOLANT_G1,
	QAWS_INTERPOLANT_SUBDIVISION,
	QAWS_INTERPOLANT_MLS
} qaws_interpolant;
#endif

typedef struct qaws_mesh_reference_desc
{
	qaws_interpolant interpolant;
	qaws_scalar const* normals;          /* 3 per vertex; NULL: computed, split at the mesh's feature edges */
	int seams;                          /* nonzero: the mesh's attribute seams (values jumping across an edge: corner or face
	                                       values, nonzero edge values) are creases of the reference (its regions split there,
	                                       computed normals too) and remembered as seams: a remesh on it keeps them as lines
	                                       that are no features of its result (qaws_mesh_remesh's own reference does this) */
	qaws_scalar mls_radius;             /* MLS: each vertex's kernel radius as a multiple of the mean length of its edges (0: 2.5) */
} qaws_mesh_reference_desc;

typedef struct qaws_mesh_reference qaws_mesh_reference;

/* desc NULL: PN with computed normals */
qaws_status qaws_mesh_reference_create(qaws_mesh_reference** out_reference, qaws_mesh_3d const* mesh, qaws_mesh_reference_desc const* desc);
void qaws_mesh_reference_destroy(qaws_mesh_reference* reference);

unsigned int qaws_mesh_reference_get_triangle_count(qaws_mesh_reference const* reference);
qaws_status qaws_mesh_reference_get_triangle(unsigned int out_vertices[3], qaws_mesh_reference const* reference, unsigned int triangle);
/* the mesh the reference copied */
qaws_mesh_3d const* qaws_mesh_reference_get_mesh(qaws_mesh_reference const* reference);
/* a location of an input vertex */
qaws_status qaws_mesh_reference_get_vertex_location(qaws_mesh_location* out_location, qaws_mesh_reference const* reference,
	unsigned int vertex);

/* per location: a position (3), a unit normal (3), curvatures (8: k1 >= k2,
   then the unit directions d1, d2); each output may be NULL */
qaws_status qaws_mesh_reference_evaluate(qaws_scalar* out_positions, qaws_scalar* out_normals, qaws_scalar* out_curvatures,
	qaws_mesh_reference const* reference, qaws_mesh_location const* locations, unsigned int count);

/* hints and out_distances may be NULL */
qaws_status qaws_mesh_reference_project(qaws_mesh_location* out_locations, qaws_scalar* out_distances, qaws_mesh_reference const* reference,
	qaws_scalar const* points, qaws_mesh_location const* hints, unsigned int count);

/* tangent: a 3D vector, its part in the plane of the triangle walked; its
   length is the distance walked (on the triangles); stops at the border */
qaws_status qaws_mesh_reference_walk(qaws_mesh_location* out_to, qaws_mesh_reference const* reference, qaws_scalar const* tangent,
	qaws_mesh_location from);

/* the input mesh's attribute at the locations (components per location) */
qaws_status qaws_mesh_reference_evaluate_attribute(qaws_scalar* out_values, qaws_mesh_reference const* reference,
	qaws_mesh_location const* locations, unsigned int count, unsigned int attribute);

#endif /* QAWS_MESH_REFERENCE_H */
