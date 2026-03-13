#ifndef QAWS_BOOLEAN_2D_H
#define QAWS_BOOLEAN_2D_H

#include "qaws_types.h"
#include "qaws_status.h"

typedef enum qaws_boolean_op {
	QAWS_BOOLEAN_UNION = 0,
	QAWS_BOOLEAN_INTERSECTION,
	QAWS_BOOLEAN_DIFFERENCE
} qaws_boolean_op;

/* Perform a boolean operation on two planar regions.
   Each region is defined by a closed 2D curve (boundary).
   The result is one or more closed 2D curves representing the boundary
   of the resulting region.
   Uses curve-curve intersection + winding number classification.
   Curves must be 2D and closed (or treated as closed by connecting endpoints). */
qaws_status qaws_boolean_2d(
	qaws_curve const* region_a,
	qaws_curve const* region_b,
	qaws_boolean_op operation,
	qaws_curve** out_boundaries,
	unsigned int boundary_capacity,
	unsigned int* out_count);

#endif /* QAWS_BOOLEAN_2D_H */
