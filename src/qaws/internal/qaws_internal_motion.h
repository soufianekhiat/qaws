#ifndef QAWS_INTERNAL_MOTION_H
#define QAWS_INTERNAL_MOTION_H
#include "../qaws_motion.h"

/* bounds over t in [a, b] (qaws_motion.c): a rigid motion's |T''| into *bt and |omega'| + |omega|^2 into *br (so a point p of
   the body moves with |x''| <= *bt + *br |p|), a vertex motion's largest |x''| into *bt (*br 0); 1, or 0 for a motion
   along a curve (no bound known). Used by chabak's sweeps. */
#define mo_second_bound qaws_internal_motion_second_bound
int qaws_internal_motion_second_bound(qaws_motion const* m, double a, double b, double* bt, double* br);

#endif /* QAWS_INTERNAL_MOTION_H */
