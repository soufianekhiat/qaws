#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_H_F32
#define QAWS_H_F32
#define QAWS_H_BODY
#endif
#else
#ifndef QAWS_H_F64
#define QAWS_H_F64
#define QAWS_H_BODY
#endif
#endif
#ifdef QAWS_H_BODY
#undef QAWS_H_BODY
#define QAWS_H

#define QAWS_VERSION_MAJOR 3
#define QAWS_VERSION_MINOR 0
#define QAWS_VERSION_PATCH 0
#define QAWS_VERSION ((QAWS_VERSION_MAJOR * 10000) + (QAWS_VERSION_MINOR * 100) + QAWS_VERSION_PATCH)
#define QAWS_VERSION_STRING "3.0.0"

#include "qaws_types.h"
#include "qaws_status.h"
#include "qaws_curve.h"
#include "qaws_bezier.h"
#include "qaws_hermite.h"
#include "qaws_catmull_rom.h"
#include "qaws_bspline.h"
#include "qaws_nurbs.h"
#include "qaws_trajectory.h"
#include "qaws_yuksel.h"
#include "qaws_rational_bezier.h"
#include "qaws_composite.h"
#include "qaws_arc.h"
#include "qaws_polynomial.h"
#include "qaws_clothoid.h"
#include "qaws_subdivision.h"
#include "qaws_eval.h"
#include "qaws_sampling.h"
#include "qaws_traversal.h"
#include "qaws_inspect.h"
#include "qaws_curve_batch.h"
#include "qaws_view.h"
#include "qaws_operations.h"
#include "qaws_convert.h"
#include "qaws_export.h"
#include "qaws_import.h"
#include "qaws_surface.h"
#include "qaws_surface_bezier.h"
#include "qaws_surface_bspline.h"
#include "qaws_surface_nurbs.h"
#include "qaws_surface_swept.h"
#include "qaws_surface_ruled.h"
#include "qaws_surface_coons.h"
#include "qaws_surface_extrusion.h"
#include "qaws_surface_revolution.h"
#include "qaws_surface_bilinear.h"
#include "qaws_surface_biquadratic.h"
#include "qaws_surface_pipe.h"
#include "qaws_surface_loft.h"
#include "qaws_surface_gordon.h"
#include "qaws_surface_offset.h"
#include "qaws_surface_trim.h"
#include "qaws_surface_intersect.h"
#include "qaws_surface_tspline.h"
#include "qaws_surface_subdiv.h"
#include "qaws_surface_fillet.h"
#include "qaws_surface_batch.h"
#include "qaws_brep.h"
#include "qaws_boolean_2d.h"
#include "qaws_path.h"
#include "qaws_clip.h"
#include "qaws_clip64.h"
#include "qaws_clip_exact.h"
#include "qaws_offset.h"
#include "qaws_rectclip.h"
#include "qaws_stack.h"
#include "qaws_mesh.h"
#include "qaws_mesh_reference.h"
#include "qaws_mesh_field.h"
#include "qaws_remesh.h"
#include "qaws_quadrangulate.h"
#include "qaws_mesh_generate.h"
#include "qaws_mesh_measure.h"
#include "qaws_mesh_repair.h"
#include "qaws_mesh_bvh.h"
#include "qaws_quality_report.h"
#include "qaws_mesh_watertight.h"
#include "qaws_volume_mesh.h"
#include "qaws_tetrahedralize.h"
#include "qaws_hexahedralize.h"
#include "qaws_boundary_layers.h"
#include "qaws_triangulate.h"
#include "qaws_region_sample.h"
#include "qaws_mean_value.h"
#include "qaws_surface_region.h"
#include "qaws_repulsive.h"
#include "qaws_shell_region.h"
#include "qaws_skeleton.h"
#include "qaws_alloc.h"
#include "qaws_inline.h"
#include "qaws_diff.h"
#include "qaws_diff_geometry.h"
#include "qaws_diff_ops.h"
#include "qaws_diff_map.h"
#include "qaws_diff_functionals.h"
#include "qaws_diff_sampling.h"

#endif /* QAWS_H */
