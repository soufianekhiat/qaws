#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_MESH_GENERATE_H_F32
#define QAWS_MESH_GENERATE_H_F32
#define QAWS_MESH_GENERATE_H_BODY
#endif
#else
#ifndef QAWS_MESH_GENERATE_H_F64
#define QAWS_MESH_GENERATE_H_F64
#define QAWS_MESH_GENERATE_H_BODY
#endif
#endif
#ifdef QAWS_MESH_GENERATE_H_BODY
#undef QAWS_MESH_GENERATE_H_BODY
#define QAWS_MESH_GENERATE_H

#include "qaws_mesh_reference.h"
#include "qaws_shell_region.h"

/*
 * Meshes of exact geometry: a surface (NURBS, any qaws surface) or a shell
 * (a B-rep: faces, trimmed or not, meeting along shared edges).
 *
 *   REMESH  the shell meshed densely and conformally (qaws_shell_region_
 *           triangulate on the whole shell, half the target edge), the
 *           faces' shared edges and the free sides made features; that mesh's
 *           PN reference remeshed (triangles) or quadrangulated (quads) at
 *           the target size; every vertex then snapped onto its face's exact
 *           surface: its footpoint the face and its (u, v), its position the
 *           surface's there. Elements never straddle a shared edge. The
 *           default.
 *   PARAMETRIC, RESTRICTED_DELAUNAY, STRUCTURED  reserved
 *           (UNSUPPORTED_OPERATION today): meshing in (u, v) with the metric,
 *           with Delaunay guarantees, iso-parametric grids.
 *
 * Sizes: a uniform target length, or adaptive by tolerances (the chord and
 * angle tolerances of qaws_curvature_size_desc, clamped and graded) read
 * from the dense mesh's curvature. The output's footpoints are on the exact
 * geometry: patch = the face (0 for a surface), a = u, b = v. A face
 * attribute "face" holds each element's face. The report's deviation is how
 * far the snapping moved the vertices (the dense reference's error).
 */

#ifndef QAWS_ENUM_QAWS_GENERATE_METHOD
#define QAWS_ENUM_QAWS_GENERATE_METHOD
typedef enum qaws_generate_method
{
	QAWS_GENERATE_REMESH = 0,
	QAWS_GENERATE_PARAMETRIC,
	QAWS_GENERATE_RESTRICTED_DELAUNAY,
	QAWS_GENERATE_STRUCTURED
} qaws_generate_method;
#endif

#ifndef QAWS_ENUM_QAWS_ELEMENTS
#define QAWS_ENUM_QAWS_ELEMENTS
typedef enum qaws_elements
{
	QAWS_ELEMENTS_TRIANGLES = 0,
	QAWS_ELEMENTS_QUADS,                /* all quads */
	QAWS_ELEMENTS_QUAD_DOMINANT
} qaws_elements;
#endif

typedef struct qaws_generate_desc
{
	qaws_generate_method method;
	qaws_elements elements;
	qaws_scalar target_length;          /* uniform; 0: from the tolerances, else a 30th of the square root of the area */
	qaws_scalar chord_tolerance;        /* adaptive (with angle_tolerance): see qaws_curvature_size_desc */
	qaws_scalar angle_tolerance;
	qaws_scalar min_length, max_length;
	qaws_scalar gradation;
} qaws_generate_desc;

qaws_status qaws_mesh_generate_surface(qaws_mesh_3d** out_mesh, qaws_mesh_report* out_report, qaws_surface const* surface,
	qaws_generate_desc const* desc);
qaws_status qaws_mesh_generate_shell(qaws_mesh_3d** out_mesh, qaws_mesh_report* out_report, qaws_shell_atlas const* shell,
	qaws_generate_desc const* desc);

#endif /* QAWS_MESH_GENERATE_H */
