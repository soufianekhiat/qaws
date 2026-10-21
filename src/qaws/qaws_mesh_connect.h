#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_MESH_CONNECT_H_F32
#define QAWS_MESH_CONNECT_H_F32
#define QAWS_MESH_CONNECT_H_BODY
#endif
#else
#ifndef QAWS_MESH_CONNECT_H_F64
#define QAWS_MESH_CONNECT_H_F64
#define QAWS_MESH_CONNECT_H_BODY
#endif
#endif
#ifdef QAWS_MESH_CONNECT_H_BODY
#undef QAWS_MESH_CONNECT_H_BODY
#define QAWS_MESH_CONNECT_H

#include "qaws_mesh.h"
#include "qaws_triangulate.h"
#include "qaws_volume_mesh.h"

/*
 * The mesh connector: two borders (one on each mesh, or two on one mesh)
 * joined by a bridge, the result one mesh.
 *
 * Surfaces (qaws_mesh_3d_connect). Each border is named by a vertex on it:
 * its whole loop, or with QAWS_CONNECT_CHAINS the chain from a_vertex to
 * a_end along a's border and the chain from b_vertex to b_end against b's
 * (a_vertex meets b_vertex, a_end meets b_end): a strip instead of a tube.
 * The borders are matched by arc length, b's turned the other way (so the
 * bridge's faces agree with both meshes' orientations), a loop started
 * where the twist is least. Each ring of the bridge is a Hermite blend of
 * the two borders' points at its arc length, leaving each border along its
 * mesh's surface (the direction across the border, in the vertex's tangent
 * plane, away from the mesh), so the bridge meets both meshes tangent to
 * them. The rings' count follows the gap over the borders' mean edge length
 * (or rows), their spacing graded from one side's edge length to the
 * other's, and their vertex counts from one border's to the other's.
 * Between rings of equal counts the faces are quads (with quads set); else
 * triangles zipped by arc length. With QAWS_CONNECT_FAIR the bridge's inner
 * vertices are smoothed (three bi-Laplacian flow steps, qaws_mesh_3d_smooth)
 * with both meshes held: the zipped triangles evened out, the shape kept
 * (a thin plate spanning two borders far apart pinches, as a soap film
 * does).
 *
 * The result's vertices are a's, then b's, then the bridge's; a face
 * attribute "part" says 0 (a), 1 (the bridge) or 2 (b). Both meshes'
 * feature edges are kept. The attributes a and b both have (by name,
 * element and components; all of a's when the borders are on one mesh) are
 * carried: on the bridge a vertex's (or a corner's) value blends the two
 * borders' values at its arc length by its Hermite weight (NEAREST: the
 * nearer border's; NORMALIZED: renormalized), a face's is the nearer
 * border's face's. The report measures the bridge (qaws_quality_report.h's
 * measures): its smallest and largest angles, its smallest scaled Jacobian.
 *
 * Planar meshes (qaws_mesh_2d_connect): the same in their plane (the
 * surfaces at z = 0); attributes and features are not carried.
 *
 * Volumes (qaws_volume_mesh_connect): two patches of border facets (facet
 * indices as qaws_volume_mesh_get_cell_adjacency numbers them) joined by
 * layers of cells, prisms on triangles and hexahedra on quads, the layers
 * Hermite-blended between the patches along their surfaces' normals. The
 * bridge conforms to both meshes, so the patches must be the same
 * triangulation up to placement: each a disk, matched through its Tutte
 * embedding (borders by arc length, the least twist) vertex for vertex,
 * facet for facet; INVALID_ARGUMENT when they are not.
 */

#define QAWS_CONNECT_FAIR 1u                /* the bridge's inner vertices smoothed (bi-Laplacian flow) */
#define QAWS_CONNECT_CHAINS 2u              /* chains of the borders (a_end, b_end), not whole loops: a strip */

typedef struct qaws_connect_desc
{
	unsigned int rows;                  /* faces across the bridge (layers of cells for volumes); 0: the gap over the mean edge */
	int quads;                          /* quads between rings of equal counts */
	unsigned int flags;                 /* QAWS_CONNECT_* */
	unsigned int a_end, b_end;          /* with QAWS_CONNECT_CHAINS: the chains' last vertices */
} qaws_connect_desc;

typedef struct qaws_connect_report
{
	unsigned int rows;
	unsigned int bridge_faces;          /* faces, or cells for volumes */
	unsigned int loop_vertices[2];      /* the two borders' (patches' for volumes) */
	qaws_scalar gap;                    /* the borders' matched points' mean distance */
	qaws_scalar min_angle;              /* the bridge's smallest corner angle (radians; volumes: dihedral) */
	qaws_scalar max_angle;              /* its largest */
	qaws_scalar min_scaled_jacobian;    /* its smallest scaled Jacobian (1 for a square or a cube) */
} qaws_connect_report;

/* b NULL: both borders on a (two different ones). INVALID_ARGUMENT for a vertex not on a border, one loop named twice, a loop
   that does not close, a chain whose end is not on its border */
qaws_status qaws_mesh_3d_connect(qaws_mesh_3d** out_mesh, qaws_connect_report* report, qaws_mesh_3d const* a, unsigned int a_vertex,
	qaws_mesh_3d const* b, unsigned int b_vertex, qaws_connect_desc const* desc);

qaws_status qaws_mesh_2d_connect(qaws_mesh_2d** out_mesh, qaws_connect_report* report, qaws_mesh_2d const* a, unsigned int a_vertex,
	qaws_mesh_2d const* b, unsigned int b_vertex, qaws_connect_desc const* desc);

/* b NULL: both patches on a */
qaws_status qaws_volume_mesh_connect(qaws_volume_mesh** out_mesh, qaws_connect_report* report, qaws_volume_mesh const* a,
	unsigned int const* a_facets, qaws_volume_mesh const* b, unsigned int const* b_facets, unsigned int a_facet_count, unsigned int b_facet_count,
	qaws_connect_desc const* desc);

#endif /* QAWS_MESH_CONNECT_H */
