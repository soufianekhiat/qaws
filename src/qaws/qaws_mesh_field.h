#ifndef QAWS_MESH_FIELD_H
#define QAWS_MESH_FIELD_H

#include "qaws_mesh_reference.h"

/*
 * Size fields: the target edge length at each point of a surface, for the
 * meshing operators (adaptive remeshing, generation).
 *
 *   uniform      one length everywhere;
 *   curvature    from a reference's curvature and tolerances: for a largest
 *                curvature k (radius R = 1 / k) an edge of length L leaves
 *                the surface by R - sqrt(R^2 - L^2 / 4) and its normal
 *                turns by 2 asin(L / 2R), so the chord tolerance e gives
 *                L = 2 sqrt(e (2R - e)) and the angle tolerance t gives
 *                L = 2R sin(t / 2); the smaller one, clamped to
 *                [min_length, max_length], sampled at the reference's
 *                vertices, its growth limited (gradation: a size grows by
 *                at most `gradation` times the distance), interpolated
 *                across the triangles;
 *   vertices     sizes given at a reference's mesh vertices, interpolated;
 *   callback     a function of the position;
 *   min          the smaller of two fields.
 *
 * A field sampled on a reference is evaluated at that reference's
 * locations; the others at positions (the location is unused). Fields are
 * immutable; a field built on others or on a reference keeps pointers to
 * them, which must outlive it.
 */

typedef struct qaws_size_field qaws_size_field;

typedef struct qaws_curvature_size_desc
{
	qaws_scalar chord_tolerance;        /* the largest gap between an edge and the surface; 0: unused */
	qaws_scalar angle_tolerance;        /* the largest turn of the normal along an edge, radians; 0: unused */
	qaws_scalar min_length;             /* 0: a thousandth of the reference's size */
	qaws_scalar max_length;             /* 0: a tenth of the reference's size */
	qaws_scalar gradation;              /* the largest growth of the size per unit of distance; 0: unlimited */
} qaws_curvature_size_desc;

qaws_status qaws_size_field_create_uniform(qaws_scalar length, qaws_size_field** out_field);
qaws_status qaws_size_field_create_curvature(qaws_mesh_reference const* reference, qaws_curvature_size_desc const* desc,
	qaws_size_field** out_field);
/* sizes: one per vertex of the reference's mesh (all positive) */
qaws_status qaws_size_field_create_vertices(qaws_mesh_reference const* reference, qaws_scalar const* sizes, qaws_size_field** out_field);
qaws_status qaws_size_field_create_callback(qaws_scalar (*size)(void* user, qaws_scalar const* point), void* user, qaws_size_field** out_field);
qaws_status qaws_size_field_create_min(qaws_size_field const* a, qaws_size_field const* b, qaws_size_field** out_field);
void qaws_size_field_destroy(qaws_size_field* field);

/* the reference a field is sampled on (through min fields: the first found), or NULL */
qaws_mesh_reference const* qaws_size_field_get_reference(qaws_size_field const* field);

/* sizes at points; locations on the field's reference (may be NULL for a
   field on none: then they are unused) */
qaws_status qaws_size_field_evaluate(qaws_size_field const* field, qaws_mesh_location const* locations, qaws_scalar const* points,
	unsigned int count, qaws_scalar* out_sizes);

/*
 * Direction fields: a cross field on a mesh's vertices (4-symmetric: a
 * direction and its quarter turns are one), smoothed by Instant Meshes'
 * multiresolution Gauss-Seidel (pairs of vertices merged level by level,
 * each level smoothed and handed to the finer one), along feature edges and
 * the border. Its index round each face, in quarter turns, marks the
 * singularities: -1 and +1 where a quad mesh following the field has
 * vertices of valence 3 and 5.
 */
typedef struct qaws_direction_field_desc
{
	unsigned int symmetry;              /* 4: crosses (0: 4); 2 and 6 reserved */
	unsigned int sweeps;                /* Gauss-Seidel sweeps per level; 0: 8 */
} qaws_direction_field_desc;

/* per vertex a unit tangent (3); per face the index round it; either output
   may be NULL; desc NULL: the defaults */
qaws_status qaws_mesh_3d_direction_field(qaws_mesh_3d const* mesh, qaws_direction_field_desc const* desc, qaws_scalar* out_directions,
	int* out_face_index);

#endif /* QAWS_MESH_FIELD_H */
