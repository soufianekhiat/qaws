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
 * Direction fields: an n-direction field on a mesh's vertices (a direction
 * and its turns by 1/n of a turn are one; crosses for n = 4), along feature
 * edges and the border.
 *
 *   LOCAL      Instant Meshes' multiresolution Gauss-Seidel (pairs of
 *              vertices merged level by level, each level smoothed and
 *              handed to the finer one); crosses only. The default.
 *   SMOOTHEST  the globally smoothest field (Knoppel, Crane, Pinkall,
 *              Schroder 2013): one global solve (a sparse LDL^T of the
 *              vertices' connection Laplacian), the vertices of feature
 *              lines and the border held along their lines, corners free;
 *              on a closed surface without lines the smallest eigenvector.
 *              With `curvature_aligned`, drawn to the principal curvature
 *              directions where the surface bends one way more than the
 *              other. n = 2, 4 or 6. The relaxed optimum: where symmetry
 *              lets zeros coincide (a disc held by its border) its
 *              singularities sit close together.
 *
 * Its index round each face, in 1/n turns, marks the singularities (they
 * sum to n times the Euler characteristic): for crosses -1 and +1 where a
 * quad mesh following the field has vertices of valence 3 and 5.
 */
typedef enum qaws_direction_field_method
{
	QAWS_DIRECTION_FIELD_LOCAL = 0,
	QAWS_DIRECTION_FIELD_SMOOTHEST
} qaws_direction_field_method;

typedef struct qaws_direction_field_desc
{
	unsigned int symmetry;              /* n: 4 (0: 4); SMOOTHEST also 2 and 6 */
	unsigned int sweeps;                /* LOCAL: Gauss-Seidel sweeps per level; 0: 8 */
	qaws_direction_field_method method; /* LOCAL: Instant Meshes, the default */
	int curvature_aligned;              /* SMOOTHEST: drawn to the principal directions where the surface bends one way more
	                                       than the other (Knoppel's lambda_t = 0) */
} qaws_direction_field_desc;

/* per vertex a unit tangent (3); per face the index round it; either output
   may be NULL; desc NULL: the defaults */
qaws_status qaws_mesh_3d_direction_field(qaws_mesh_3d const* mesh, qaws_direction_field_desc const* desc, qaws_scalar* out_directions,
	int* out_face_index);

/*
 * Frame fields in volumes: per tetrahedron an octahedral frame (three
 * orthonormal axes, one with its 24 turns that map a cube to itself),
 * smooth across facets, the cells on the border holding an axis along its
 * normal (two where they touch two faces of the border meeting square).
 * Each frame is taken as its function on the sphere (band-4 spherical
 * harmonics, Ray, Sokolov, Levy 2016): the smoothest field of such
 * functions comes from one linear solve over all cells, each cell's frame
 * then the one nearest its function, polished by Gauss-Seidel sweeps over
 * the cells' stored adjacency (neighbours matched by the signed permutation
 * of their axes nearest). An edge is singular where the matchings round it
 * compose to a turn: hexahedra following the field meet there three or five
 * round an edge.
 */
typedef struct qaws_frame_field_desc
{
	unsigned int sweeps;                /* the polishing sweeps; 0: 60 */
} qaws_frame_field_desc;

struct qaws_volume_mesh;

/* per cell three unit axes (9); per cell its singular edges (may be NULL):
   bit e for edge e of 0-1, 0-2, 0-3, 1-2, 1-3, 2-3; tetrahedral meshes only
   (UNSUPPORTED_OPERATION otherwise); desc NULL: the defaults */
qaws_status qaws_volume_mesh_frame_field(struct qaws_volume_mesh const* mesh, qaws_frame_field_desc const* desc, qaws_scalar* out_frames,
	unsigned char* out_singular);

#endif /* QAWS_MESH_FIELD_H */
