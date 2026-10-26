#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_VOLUME_MESH_H_F32
#define QAWS_VOLUME_MESH_H_F32
#define QAWS_VOLUME_MESH_H_BODY
#endif
#else
#ifndef QAWS_VOLUME_MESH_H_F64
#define QAWS_VOLUME_MESH_H_F64
#define QAWS_VOLUME_MESH_H_BODY
#endif
#endif
#ifdef QAWS_VOLUME_MESH_H_BODY
#undef QAWS_VOLUME_MESH_H_BODY
#define QAWS_VOLUME_MESH_H

#include "qaws_mesh.h"

/*
 * Volume meshes: tetrahedra, hexahedra, prisms and pyramids.
 *
 * Storage follows geogram's cells: a cell is a run of corners, each holding
 * its vertex, and every cell facet holds the facet across it (the facets'
 * adjacency, stored as the corner table stores opposite corners; none on
 * the border). Cells are read back in bulk, data enters and leaves as
 * arrays.
 *
 * Local numbering (VTK's, with every cell positively oriented):
 *   tetrahedron  0 1 2 3, det(v1 - v0, v2 - v0, v3 - v0) > 0;
 *                facets 1 2 3 | 0 3 2 | 0 1 3 | 0 2 1 (facet i opposite vertex i)
 *   hexahedron   0 1 2 3 the base (its right-hand normal toward the top), 4 5 6 7 above them;
 *                facets 0 3 2 1 | 4 5 6 7 | 0 1 5 4 | 1 2 6 5 | 2 3 7 6 | 3 0 4 7
 *   prism        0 1 2 the base (its normal toward the top), 3 4 5 above them;
 *                facets 0 2 1 | 3 4 5 | 0 1 4 3 | 1 2 5 4 | 2 0 3 5
 *   pyramid      0 1 2 3 the base (its normal toward the apex), 4 the apex;
 *                facets 0 3 2 1 | 0 1 4 | 1 2 4 | 2 3 4 | 3 0 4
 * Facets turn counter-clockwise seen from outside their cell. A mesh's
 * facets are its cells' facets in turn (qaws_volume_mesh_get_facet_count).
 *
 * Creation matches facets with the same vertices; a facet shared by more
 * than two cells, or by two cells on the same side, is INVALID_TOPOLOGY.
 */

#ifndef QAWS_ENUM_QAWS_CELL_TYPE
#define QAWS_ENUM_QAWS_CELL_TYPE
typedef enum qaws_cell_type
{
	QAWS_CELL_TETRAHEDRON = 0,
	QAWS_CELL_HEXAHEDRON,
	QAWS_CELL_PRISM,
	QAWS_CELL_PYRAMID
} qaws_cell_type;
#endif

typedef struct qaws_volume_mesh qaws_volume_mesh;

typedef struct qaws_volume_mesh_desc
{
	qaws_scalar const* positions;           /* 3 per vertex */
	unsigned int vertex_count;
	unsigned int const* cell_vertices;      /* each cell's corners in turn */
	qaws_cell_type const* cell_types;       /* per cell; NULL: all tetrahedra */
	unsigned int cell_count;
} qaws_volume_mesh_desc;

qaws_status qaws_volume_mesh_create(qaws_volume_mesh** out_mesh, qaws_volume_mesh_desc const* desc);
qaws_status qaws_volume_mesh_clone(qaws_volume_mesh** out_mesh, qaws_volume_mesh const* mesh);
void qaws_volume_mesh_destroy(qaws_volume_mesh* mesh);

unsigned int qaws_volume_mesh_get_vertex_count(qaws_volume_mesh const* mesh);
unsigned int qaws_volume_mesh_get_cell_count(qaws_volume_mesh const* mesh);
unsigned int qaws_volume_mesh_get_corner_count(qaws_volume_mesh const* mesh);   /* the cells' vertex lists' total length */
unsigned int qaws_volume_mesh_get_facet_count(qaws_volume_mesh const* mesh);    /* 4 per tetrahedron, 6 per hexahedron, 5 per prism or pyramid */
qaws_status qaws_volume_mesh_get_positions(qaws_scalar* out_positions, qaws_volume_mesh const* mesh);
/* each cell's vertices in turn; out_cell_types may be NULL */
qaws_status qaws_volume_mesh_get_cells(unsigned int* out_cell_vertices, qaws_cell_type* out_cell_types, qaws_volume_mesh const* mesh);
/* per facet (cells in turn, each in its type's facet order): the cell across, QAWS_MESH_NONE on the border */
qaws_status qaws_volume_mesh_get_cell_adjacency(unsigned int* out_cells, qaws_volume_mesh const* mesh);

/* The border: its facets as a surface mesh turned outward, with a vertex
   attribute "vertex" (the volume's vertex) and a face attribute "cell".
   The volume's attributes come along by name: vertex ones per vertex, corner
   ones per corner (the facet's corners in its cell), cell and facet ones per
   face (the facet's cell's, the facet's own). */
qaws_status qaws_volume_mesh_get_boundary(qaws_mesh_3d** out_surface, qaws_volume_mesh const* mesh);

/* The structure's invariants: every cell positively oriented, facets
   matched pairwise with the same vertices turned the other way. OK, or
   INVALID_TOPOLOGY. */
qaws_status qaws_volume_mesh_validate(qaws_volume_mesh const* mesh);

/* ---- attributes -------------------------------------------------------- */

/*
 * Values per vertex, per corner (qaws_volume_mesh_get_cells order), per cell
 * (QAWS_MESH_CELL) or per facet (QAWS_MESH_FACET, each side of an inner
 * facet its own), mixed by their rules as on surfaces. Cloning and the border
 * carry them; qaws_volume_mesh_improve keeps them, flipping only inside
 * regions where they are continuous (a facet between cells of different
 * values, with corner values differing at its vertices, or a nonzero facet
 * value stays as a wall, its vertices fixed); qaws_volume_mesh_tetrahedralize
 * and qaws_volume_mesh_hexahedralize take the border surface's (vertex and
 * corner ones per vertex, at each vertex's nearest point of the surface; face
 * ones on the border facets); qaws_volume_mesh_grow_layers maps them up its
 * stacks (face ones per cell); qaws_volume_mesh_connect keeps both volumes'
 * and mixes them across the bridge.
 */

/* values: components per element; NULL: zeros. INVALID_ARGUMENT for an edge or face element, a rule and component count
   that do not go together (as qaws_mesh_3d_add_attribute) */
qaws_status qaws_volume_mesh_add_attribute(unsigned int* out_attribute, qaws_volume_mesh* mesh, qaws_scalar const* values,
	qaws_attribute_desc const* desc);
unsigned int qaws_volume_mesh_get_attribute_count(qaws_volume_mesh const* mesh);
/* OUT_OF_RANGE when there is none of that name */
qaws_status qaws_volume_mesh_find_attribute(unsigned int* out_attribute, qaws_volume_mesh const* mesh, char const* name);
/* the desc's name points into the mesh, valid while the attribute lives */
qaws_status qaws_volume_mesh_get_attribute_desc(qaws_attribute_desc* out_desc, qaws_volume_mesh const* mesh, unsigned int attribute);
qaws_status qaws_volume_mesh_get_attribute(qaws_scalar* out_values, qaws_volume_mesh const* mesh, unsigned int attribute);
qaws_status qaws_volume_mesh_set_attribute(qaws_volume_mesh* mesh, unsigned int attribute, qaws_scalar const* values);
qaws_status qaws_volume_mesh_remove_attribute(qaws_volume_mesh* mesh, unsigned int attribute);   /* later ones shift down */

typedef struct qaws_volume_quality
{
	unsigned int tetrahedra, hexahedra, prisms, pyramids;
	qaws_scalar min_dihedral, max_dihedral;     /* degrees, over the tetrahedra (70.5 regular) */
	qaws_scalar max_radius_edge;                /* tetrahedra: circumradius over shortest edge (0.61 regular) */
	qaws_scalar min_scaled_jacobian;            /* hexahedra, at their corners (1 for a cube) */
	qaws_scalar min_edge, max_edge, mean_edge;  /* over the cells' edges */
	qaws_scalar volume;
} qaws_volume_quality;

qaws_status qaws_volume_mesh_measure_quality(qaws_volume_quality* out_quality, qaws_volume_mesh const* mesh);

#endif /* QAWS_VOLUME_MESH_H */
