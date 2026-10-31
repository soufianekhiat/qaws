#ifndef QAWS_INTERNAL_ARRAY_H
#define QAWS_INTERNAL_ARRAY_H
#include "../qaws_platform.h"
#include "../qaws_view.h"
#include "../qaws_types.h"
#include "../qaws_status.h"
#include "../qaws_path.h"
#include "../qaws_batch_executor.h"
#include "../qaws_motion.h"

/* element i's component c of a host array, read bytewise (qaws_view.c) */
double qaws_internal_array_get(qaws_array_view const* v, unsigned int i, unsigned int c);

#endif /* QAWS_INTERNAL_ARRAY_H */
