#ifndef QAWS_INTERNAL_PATH_H
#define QAWS_INTERNAL_PATH_H

#include "qaws_internal_types.h"

/* Integral of x dy along a 2D curve from t0 to t1 (negative when t0 > t1),
   span by span, to an absolute tolerance. */
qaws_status qaws_internal_curve_area_2d(qaws_curve const* c, double t0, double t1, double tol, double* out);

#endif /* QAWS_INTERNAL_PATH_H */
