#ifndef QAWS_MESH_REFERENCE_H
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
 *
 * PHONG, NAGATA, G1, SUBDIVISION and MLS are reserved for further
 * interpolants (UNSUPPORTED_OPERATION today).
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

typedef struct qaws_mesh_reference_desc
{
	qaws_interpolant interpolant;
	qaws_scalar const* normals;          /* 3 per vertex; NULL: computed, split at the mesh's feature edges */
} qaws_mesh_reference_desc;

typedef struct qaws_mesh_reference qaws_mesh_reference;

/* desc NULL: PN with computed normals */
qaws_status qaws_mesh_reference_create(qaws_mesh_3d const* mesh, qaws_mesh_reference_desc const* desc, qaws_mesh_reference** out_reference);
void qaws_mesh_reference_destroy(qaws_mesh_reference* reference);

unsigned int qaws_mesh_reference_get_triangle_count(qaws_mesh_reference const* reference);
qaws_status qaws_mesh_reference_get_triangle(qaws_mesh_reference const* reference, unsigned int triangle, unsigned int out_vertices[3]);
/* the mesh the reference copied */
qaws_mesh_3d const* qaws_mesh_reference_get_mesh(qaws_mesh_reference const* reference);
/* a location of an input vertex */
qaws_status qaws_mesh_reference_get_vertex_location(qaws_mesh_reference const* reference, unsigned int vertex,
	qaws_mesh_location* out_location);

/* per location: a position (3), a unit normal (3), curvatures (8: k1 >= k2,
   then the unit directions d1, d2); each output may be NULL */
qaws_status qaws_mesh_reference_evaluate(qaws_mesh_reference const* reference, qaws_mesh_location const* locations, unsigned int count,
	qaws_scalar* out_positions, qaws_scalar* out_normals, qaws_scalar* out_curvatures);

/* hints and out_distances may be NULL */
qaws_status qaws_mesh_reference_project(qaws_mesh_reference const* reference, qaws_scalar const* points, unsigned int count,
	qaws_mesh_location const* hints, qaws_mesh_location* out_locations, qaws_scalar* out_distances);

/* tangent: a 3D vector, its part in the plane of the triangle walked; its
   length is the distance walked (on the triangles); stops at the border */
qaws_status qaws_mesh_reference_walk(qaws_mesh_reference const* reference, qaws_mesh_location from, qaws_scalar const* tangent,
	qaws_mesh_location* out_to);

/* the input mesh's attribute at the locations (components per location) */
qaws_status qaws_mesh_reference_evaluate_attribute(qaws_mesh_reference const* reference, unsigned int attribute,
	qaws_mesh_location const* locations, unsigned int count, qaws_scalar* out_values);

#endif /* QAWS_MESH_REFERENCE_H */
