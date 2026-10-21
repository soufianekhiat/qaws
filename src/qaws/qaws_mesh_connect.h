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

/*
 * The mesh connector: two border loops (one on each mesh, or two on one
 * mesh) joined by a bridge, the result one mesh.
 *
 * Each loop is named by a vertex on it. The loops are matched by arc
 * length, the second turned the other way (so the bridge's faces agree
 * with both meshes' orientations) and started where the twist is least.
 * The bridge is a tube of rings: each point of a ring a Hermite blend of
 * the two loops' points at its arc length, leaving each loop along its
 * mesh's surface (the direction across the border, in the vertex's tangent
 * plane, away from the mesh), so the bridge meets both meshes tangent to
 * them. The rings' count follows the gap over the loops' mean edge length
 * (or rows), their spacing graded from one side's edge length to the
 * other's, and their vertex counts from one loop's to the other's. Between
 * rings of equal counts the faces are quads (with quads set); else
 * triangles zipped by arc length. With QAWS_CONNECT_FAIR the bridge's inner
 * vertices are smoothed (three bi-Laplacian flow steps, qaws_mesh_3d_smooth)
 * with both meshes held: the zipped triangles evened out, the shape kept
 * (a thin plate spanning two loops far apart pinches, as a soap film does).
 * The result's vertices are a's, then b's, then the bridge's; a face
 * attribute "part" says 0 (a), 1 (the bridge) or 2 (b). Other attributes
 * and the features are not carried.
 *
 * Planar meshes (qaws_mesh_2d) and volume bridges (prisms, hexahedra) are
 * reserved.
 */

#define QAWS_CONNECT_FAIR 1u                /* the bridge's inner vertices smoothed (bi-Laplacian flow) */

typedef struct qaws_connect_desc
{
	unsigned int rows;                  /* faces across the bridge; 0: the gap over the loops' mean edge length */
	int quads;                          /* quads between rings of equal counts */
	unsigned int flags;                 /* QAWS_CONNECT_* */
} qaws_connect_desc;

typedef struct qaws_connect_report
{
	unsigned int rows;
	unsigned int bridge_faces;
	unsigned int loop_vertices[2];      /* the two loops' */
	qaws_scalar gap;                    /* the loops' matched points' mean distance */
	qaws_scalar min_angle;              /* the bridge's smallest corner angle (radians) */
} qaws_connect_report;

/* b NULL: both loops on a (two different ones). INVALID_ARGUMENT for a vertex not on a border, one loop named twice, a loop
   that does not close */
qaws_status qaws_mesh_3d_connect(qaws_mesh_3d** out_mesh, qaws_connect_report* report, qaws_mesh_3d const* a, unsigned int a_vertex,
	qaws_mesh_3d const* b, unsigned int b_vertex, qaws_connect_desc const* desc);

#endif /* QAWS_MESH_CONNECT_H */
