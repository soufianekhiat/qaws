#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_MESH_H_F32
#define QAWS_MESH_H_F32
#define QAWS_MESH_H_BODY
#endif
#else
#ifndef QAWS_MESH_H_F64
#define QAWS_MESH_H_F64
#define QAWS_MESH_H_BODY
#endif
#endif
#ifdef QAWS_MESH_H_BODY
#undef QAWS_MESH_H_BODY
#define QAWS_MESH_H

#include "qaws_types.h"
#include "qaws_status.h"

/*
 * Polygon meshes in space: the surface type of the meshing module.
 *
 * A qaws_mesh_3d holds triangles, quads or any polygons. It is created from
 * arrays (or by an operator: the stack mesher, remeshing, quadrangulation,
 * ...), read back in bulk, and edited by local operations that check their
 * own validity. Data enters and leaves as arrays: qaws reads and writes no
 * files.
 *
 * Storage is the faces themselves, as geogram's facet corners: a face is a
 * run of corners, each holding its vertex and the opposite corner of the
 * face across its edge (the corner table of Rossignac et al.). A corner is
 * a half-edge (from its vertex to the next corner's), so navigation is a
 * half-edge view of the faces, and faces' adjacency, the dual graph, is
 * stored rather than derived.
 *
 * Faces turn counter-clockwise seen from outside. Indices stay dense: an
 * operation that removes elements moves the last ones of their kind into
 * the freed slots (removing vertex a moves vertex n - 1 to index a), and new
 * elements are appended.
 *
 * Input that is not a manifold is split at creation, never rejected: a face
 * whose edge is already used in the same direction (a third face on an
 * edge, or a neighbour turned the other way) gets its own copies of that
 * edge's vertices, and a vertex shared by several fans of faces gets one
 * copy per fan beyond the first. Copies are appended after the given
 * vertices, so given indices keep their meaning (qaws_mesh_3d_get_vertex_source).
 *
 * Every mesh carries:
 *   attributes  values per vertex, per corner (a vertex as seen from one
 *               face: UVs and split normals) or per face, each with a rule
 *               saying how operations interpolate it;
 *   features    feature edges (creases, boundaries, seams) that operations
 *               keep, and feature corners where they meet or turn sharply;
 *   footpoints  each vertex's location on the reference an operator worked
 *               on (qaws_mesh_reference.h): the correspondence with its input.
 *
 * A mesh is not thread-safe for writes; concurrent reads are safe.
 */

#define QAWS_MESH_NONE 0xFFFFFFFFu

#ifndef QAWS_ENUM_QAWS_MESH_ELEMENT
#define QAWS_ENUM_QAWS_MESH_ELEMENT
typedef enum qaws_mesh_element
{
	QAWS_MESH_VERTEX = 0,
	QAWS_MESH_CORNER,
	QAWS_MESH_FACE,
	QAWS_MESH_EDGE                      /* in qaws_mesh_3d_get_edges order: crease weights, seam ids; an edit keeps an
	                                       edge's value on the edge and its pieces, new edges take zeros */
} qaws_mesh_element;
#endif

#ifndef QAWS_ENUM_QAWS_ATTRIBUTE_RULE
#define QAWS_ENUM_QAWS_ATTRIBUTE_RULE
typedef enum qaws_attribute_rule
{
	QAWS_ATTRIBUTE_LINEAR = 0,          /* interpolated linearly (along an edge, or barycentric at a footpoint) */
	QAWS_ATTRIBUTE_NORMALIZED,          /* linearly, then scaled to unit length (normals, tangents) */
	QAWS_ATTRIBUTE_NEAREST,             /* the nearest value: material ids, labels */
	QAWS_ATTRIBUTE_NORMAL,              /* a normal (3 components): as NORMALIZED, and negated on the faces an operation
	                                       turns over (a Boolean's reversed pieces); rotated with a swept mesh */
	QAWS_ATTRIBUTE_SRGB,                /* a colour stored sRGB-encoded (3 components, or 4 with alpha), blended in linear
	                                       light: decoded, mixed, encoded again; alpha mixed linearly, the colours weighted
	                                       by it (premultiplied) */
	QAWS_ATTRIBUTE_OKLAB,               /* the same colours blended in Oklab (perceptually even, no dark or muddy
	                                       midpoints between hues) */
	QAWS_ATTRIBUTE_TANGENT,             /* a tangent with its handedness (4 components: x y z w): the direction as
	                                       NORMALIZED, w a sign (+1 or -1, the bitangent's side), negated on the faces
	                                       an operation turns over; rotated with a swept mesh */
	QAWS_ATTRIBUTE_VECTOR,              /* a free vector (3 components: a velocity, a displacement, a flow): interpolated
	                                       linearly, its length kept (never normalized) and never negated on turned
	                                       faces; rotated with a swept mesh */
	QAWS_ATTRIBUTE_QUATERNION           /* an orientation (4 components: w x y z, unit): q and -q are one rotation, so
	                                       each value is blended on the first one's side and the mix made unit
	                                       (normalized lerp); composed with a swept mesh's rotation. Subdivision blends
	                                       them on the identity's side (w >= 0): a field near a half turn from the
	                                       identity can flip there */
} qaws_attribute_rule;
#endif

typedef struct qaws_attribute_desc
{
	char const* name;                   /* up to 31 characters; may be NULL */
	qaws_mesh_element element;
	unsigned int components;            /* 1 to 16 */
	qaws_attribute_rule rule;
} qaws_attribute_desc;

/* A point of a reference: on a mesh reference, a triangle and the
   barycentric coordinates of its second and third corners (the first's is
   1 - a - b); on a surface reference, a face and its (u, v). */
typedef struct qaws_mesh_location
{
	unsigned int patch;                 /* QAWS_MESH_NONE: unknown */
	qaws_scalar a, b;
} qaws_mesh_location;

typedef struct qaws_mesh_3d qaws_mesh_3d;

typedef struct qaws_mesh_3d_desc
{
	qaws_scalar const* positions;       /* 3 per vertex */
	unsigned int vertex_count;
	unsigned int const* face_vertices;  /* every face's vertices, counter-clockwise, face after face */
	unsigned int const* face_sizes;     /* vertices per face (3 or more); NULL: all triangles */
	unsigned int face_count;
} qaws_mesh_3d_desc;

/* INVALID_ARGUMENT for a face under 3 vertices, an index out of range or a
   vertex twice in one face. */
qaws_status qaws_mesh_3d_create(qaws_mesh_3d** out_mesh, qaws_mesh_3d_desc const* desc);
qaws_status qaws_mesh_3d_clone(qaws_mesh_3d** out_mesh, qaws_mesh_3d const* mesh);
void qaws_mesh_3d_destroy(qaws_mesh_3d* mesh);

/* ---- reading ---------------------------------------------------------- */

unsigned int qaws_mesh_3d_get_vertex_count(qaws_mesh_3d const* mesh);
unsigned int qaws_mesh_3d_get_edge_count(qaws_mesh_3d const* mesh);
unsigned int qaws_mesh_3d_get_face_count(qaws_mesh_3d const* mesh);
unsigned int qaws_mesh_3d_get_corner_count(qaws_mesh_3d const* mesh);   /* the sum of the face sizes */
qaws_vec3 qaws_mesh_3d_get_vertex(qaws_mesh_3d const* mesh, unsigned int vertex);
qaws_status qaws_mesh_3d_get_positions(qaws_scalar* out_positions, qaws_mesh_3d const* mesh);
qaws_status qaws_mesh_3d_set_positions(qaws_mesh_3d* mesh, qaws_scalar const* positions);
/* out_face_sizes may be NULL for a triangle mesh; corners are listed in this order */
qaws_status qaws_mesh_3d_get_faces(unsigned int* out_face_vertices, unsigned int* out_face_sizes, qaws_mesh_3d const* mesh);
qaws_status qaws_mesh_3d_get_edges(unsigned int* out_edge_vertices, qaws_mesh_3d const* mesh);   /* 2 per edge */

/* The faces fan-triangulated from their first vertex, in face order (for a
   triangle mesh, the faces themselves). */
unsigned int qaws_mesh_3d_get_triangle_count(qaws_mesh_3d const* mesh);
void qaws_mesh_3d_get_triangle(unsigned int out_vertices[3], qaws_mesh_3d const* mesh, unsigned int triangle);
qaws_status qaws_mesh_3d_get_triangles(unsigned int* out_triangles, qaws_mesh_3d const* mesh);   /* 3 per triangle */

/* The given vertex a copy was split from at creation (the vertex itself
   otherwise; QAWS_MESH_NONE for a vertex an operation added). */
unsigned int qaws_mesh_3d_get_vertex_source(qaws_mesh_3d const* mesh, unsigned int vertex);
/* The neighbours of a vertex, counter-clockwise; BUFFER_TOO_SMALL with the
   count when capacity is short. */
qaws_status qaws_mesh_3d_get_vertex_ring(unsigned int* out_vertices, unsigned int* out_count, qaws_mesh_3d const* mesh,
	unsigned int capacity, unsigned int vertex);
int qaws_mesh_3d_is_boundary_vertex(qaws_mesh_3d const* mesh, unsigned int vertex);
/* Unit vertex normals: the face normals around each vertex weighted by
   their angle there (Thurmer and Wuthrich). */
qaws_status qaws_mesh_3d_get_vertex_normals(qaws_scalar* out_normals, qaws_mesh_3d const* mesh);

typedef struct qaws_mesh_topology
{
	unsigned int components;            /* connected by edges */
	unsigned int boundary_loops;
	int euler_characteristic;           /* V - E + F */
	unsigned int genus;                 /* sum over the components of (2 - chi - boundaries) / 2 */
	int closed;                         /* no boundary */
} qaws_mesh_topology;

qaws_status qaws_mesh_3d_get_topology(qaws_mesh_topology* out_topology, qaws_mesh_3d const* mesh);

/* Per corner (qaws_mesh_3d_get_faces order), the face across the corner's
   edge, from its vertex to the next corner's (QAWS_MESH_NONE on the
   border): the dual graph, as stored. */
qaws_status qaws_mesh_3d_get_face_adjacency(unsigned int* out_faces, qaws_mesh_3d const* mesh);
/* The dual mesh: vertex i at face i's centroid, a face round each interior
   vertex through its faces counter-clockwise (a triangle mesh's dual is
   mostly hexagons). Border vertices have none. */
qaws_status qaws_mesh_3d_create_dual(qaws_mesh_3d** out_dual, qaws_mesh_3d const* mesh);
/* The half-edge structure's own invariants: OK, or INTERNAL_ERROR. */
qaws_status qaws_mesh_3d_validate(qaws_mesh_3d const* mesh);

/* ---- editing ------------------------------------------------------------
 *
 * Local operations with their validity checks. A refused operation returns
 * INVALID_TOPOLOGY and leaves the mesh untouched; INVALID_ARGUMENT means
 * there is no such element. Attributes follow their rules: a new vertex
 * interpolates its edge's (at the parameter of its position along it) or
 * its face's; corners on either side of the edge likewise. A vertex added or
 * moved by an edit gets the footpoint given, unknown when NULL.
 */

/* A vertex inserted on edge (a, b) at position (NULL: the midpoint): the
   triangles on either side are split in two, larger faces take the vertex.
   The new vertex is the last. A feature edge's halves are features. */
qaws_status qaws_mesh_3d_split_edge(unsigned int* out_vertex, qaws_mesh_3d* mesh, qaws_scalar const* position, qaws_mesh_location const* footpoint,
	unsigned int a, unsigned int b);

/* a merged into b, which moves to position (NULL: stays): a is removed (the
   last vertex moves to index a) with the triangles on the edge (the last
   faces move into their slots). Refused when the result would not be a
   manifold (the link condition), when a face other than the edge's own
   holds both vertices, when a is a feature corner, or when a is on a
   feature and the edge is not one of its feature edges. */
qaws_status qaws_mesh_3d_collapse_edge(qaws_mesh_3d* mesh, unsigned int a, unsigned int b, qaws_scalar const* position,
	qaws_mesh_location const* footpoint);

/* The edge of two triangles turned to join their other corners. Refused on
   a boundary or feature edge, when either face is not a triangle, or when
   the new edge exists already. */
qaws_status qaws_mesh_3d_flip_edge(qaws_mesh_3d* mesh, unsigned int a, unsigned int b);

/* A vertex inserted in a face at position (NULL: its centroid) and joined to
   each corner: a fan of triangles. The new vertex is the last; the face
   keeps its index for the first triangle, the others are appended. */
qaws_status qaws_mesh_3d_split_face(unsigned int* out_vertex, qaws_mesh_3d* mesh, qaws_scalar const* position,
	qaws_mesh_location const* footpoint, unsigned int face);

int qaws_mesh_3d_can_collapse_edge(qaws_mesh_3d const* mesh, unsigned int a, unsigned int b);
int qaws_mesh_3d_can_flip_edge(qaws_mesh_3d const* mesh, unsigned int a, unsigned int b);

/* ---- features ---------------------------------------------------------- */

/* Replaces the feature edges (2 vertices each; INVALID_ARGUMENT if one is
   not an edge). Corners: vertices on 1 or more than 2 feature edges. */
qaws_status qaws_mesh_3d_set_feature_edges(qaws_mesh_3d* mesh, unsigned int const* edge_vertices, unsigned int edge_count);

#define QAWS_FEATURES_BOUNDARY 1u       /* boundary edges are features */

/* Replaces the features by the edges whose two faces' normals turn by more
   than dihedral_angle (radians; 0: none), and with flags the boundary.
   Corners: vertices on 1 or more than 2 feature edges, and those where the
   feature line turns by more than corner_angle (radians; 0: none). */
qaws_status qaws_mesh_3d_detect_features(qaws_mesh_3d* mesh, qaws_scalar dihedral_angle, qaws_scalar corner_angle,
	unsigned int flags);
unsigned int qaws_mesh_3d_get_feature_edge_count(qaws_mesh_3d const* mesh);
qaws_status qaws_mesh_3d_get_feature_edges(unsigned int* out_edge_vertices, qaws_mesh_3d const* mesh);
int qaws_mesh_3d_is_feature_corner(qaws_mesh_3d const* mesh, unsigned int vertex);
/* The isolines of a vertex attribute (its first component) at the levels
   cut into the mesh as feature edges: each edge a level crosses split where
   the attribute, linear along it, meets the level (the new vertex takes the
   level exactly, its other attributes by their rules), each face cut
   between its two corners on the level and the quads this leaves cut by
   their shorter diagonal; edges already on a level marked. Remeshed
   (qaws_mesh_remesh) the lines stay: as creases of the cut mesh's own
   reference, or as lines across the uncut mesh's reference, the attribute
   on them the level (to the chords' sagitta in the second case). out_edge_count (may
   be NULL): the feature edges made. INVALID_ARGUMENT for an attribute not
   per vertex. */
qaws_status qaws_mesh_3d_cut_isolines(unsigned int* out_edge_count, qaws_mesh_3d* mesh, unsigned int attribute, qaws_scalar const* levels,
	unsigned int level_count);

/* ---- attributes -------------------------------------------------------- */

/* values: components per vertex, per corner (in qaws_mesh_3d_get_faces
   order), per face or per edge (qaws_mesh_3d_get_edges order); NULL: zeros.
   INVALID_ARGUMENT for a NORMAL attribute not of 3 components, an SRGB or
   OKLAB one not of 3 or 4, a TANGENT one not of 4. */
qaws_status qaws_mesh_3d_add_attribute(unsigned int* out_attribute, qaws_mesh_3d* mesh, qaws_scalar const* values,
	qaws_attribute_desc const* desc);
unsigned int qaws_mesh_3d_get_attribute_count(qaws_mesh_3d const* mesh);
qaws_status qaws_mesh_3d_find_attribute(unsigned int* out_attribute, qaws_mesh_3d const* mesh, char const* name);
/* the desc's name points into the mesh, valid while the attribute lives */
qaws_status qaws_mesh_3d_get_attribute_desc(qaws_attribute_desc* out_desc, qaws_mesh_3d const* mesh, unsigned int attribute);
qaws_status qaws_mesh_3d_get_attribute(qaws_scalar* out_values, qaws_mesh_3d const* mesh, unsigned int attribute);
qaws_status qaws_mesh_3d_set_attribute(qaws_mesh_3d* mesh, unsigned int attribute, qaws_scalar const* values);
qaws_status qaws_mesh_3d_remove_attribute(qaws_mesh_3d* mesh, unsigned int attribute);   /* later ones shift down */

/* ---- footpoints -------------------------------------------------------- */

int qaws_mesh_3d_has_footpoints(qaws_mesh_3d const* mesh);
/* INVALID_ARGUMENT when the mesh has none */
qaws_status qaws_mesh_3d_get_footpoints(qaws_mesh_location* out_locations, qaws_mesh_3d const* mesh);
qaws_status qaws_mesh_3d_set_footpoints(qaws_mesh_3d* mesh, qaws_mesh_location const* locations);   /* NULL: none */

/* ---- measures ---------------------------------------------------------- */

typedef struct qaws_mesh_quality
{
	unsigned int triangles, quads, polygons;    /* faces by size */
	qaws_scalar min_angle, max_angle;           /* degrees, over the corners of every face */
	qaws_scalar min_edge, max_edge, mean_edge;
	qaws_scalar max_aspect, mean_aspect;        /* triangles: circumradius over twice the inradius (1 when equilateral) */
	qaws_scalar area;
	unsigned int valence[8];                    /* interior vertices by valence: 3, 4, ..., 9, 10 or more */
	unsigned int irregular;                     /* interior vertices not of valence 6 (4 when quads are the most) */
} qaws_mesh_quality;

qaws_status qaws_mesh_3d_measure_quality(qaws_mesh_quality* out_quality, qaws_mesh_3d const* mesh);

/* What a meshing operator did and reached. A target not reached is not an
   error: the mesh is valid, the report says how far it got. */
typedef struct qaws_mesh_report
{
	unsigned int iterations;
	unsigned int splits, collapses, flips, moves;
	qaws_scalar max_deviation;          /* the result's face centres from the reference */
	qaws_scalar mean_deviation;
	qaws_mesh_quality quality;          /* of the result */
} qaws_mesh_report;

#endif /* QAWS_MESH_H */
