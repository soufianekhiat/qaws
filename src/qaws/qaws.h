#ifndef QAWS_H
#define QAWS_H

#define QAWS_VERSION_MAJOR 2
#define QAWS_VERSION_MINOR 0
#define QAWS_VERSION_PATCH 0
#define QAWS_VERSION ((QAWS_VERSION_MAJOR * 10000) + (QAWS_VERSION_MINOR * 100) + QAWS_VERSION_PATCH)
#define QAWS_VERSION_STRING "2.0.0"

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
#include "qaws_alloc.h"
#include "qaws_inline.h"

#endif /* QAWS_H */
