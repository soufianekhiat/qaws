/*
 * Standalone consumer of the amalgamated build.
 *
 * This is the no-build-system path: generate the pair with
 *
 *     python buildsystem/amalgamate.py --output-dir <dir>
 *
 * then compile this file together with the generated qaws.c. Nothing else is
 * needed -- no include paths into src/, no library to link, no CMake:
 *
 *     cc -std=c11 -DQAWS_SCALAR_IS_FLOAT=1 -I<dir> \
 *        examples/amalgam_standalone.c <dir>/qaws.c -lm -o smoke
 *
 * The same two files can simply be dragged into an Xcode or Visual Studio
 * project. CI compiles this on Linux, macOS and Windows to keep that promise
 * honest.
 */

#include "qaws.h"

#include <stdio.h>

int main(void)
{
	qaws_vec3 control_points[4];
	qaws_bezier_desc desc;
	qaws_curve* curve = NULL;
	qaws_eval_result_3d result;
	qaws_status status;

	control_points[0].x = 0; control_points[0].y = 0; control_points[0].z = 0;
	control_points[1].x = 1; control_points[1].y = 2; control_points[1].z = 0;
	control_points[2].x = 2; control_points[2].y = 2; control_points[2].z = 0;
	control_points[3].x = 3; control_points[3].y = 0; control_points[3].z = 0;

	desc.dimension = QAWS_DIMENSION_3D;
	desc.degree = 3;
	desc.control_points = control_points;
	desc.control_point_count = 4;

	status = qaws_curve_create_bezier(&desc, &curve);
	if (status != QAWS_STATUS_OK)
	{
		printf("create failed: %s\n", qaws_status_to_string(status));
		return 1;
	}

	status = qaws_curve_evaluate_3d(
		curve, (qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION, &result);
	if (status != QAWS_STATUS_OK)
	{
		printf("evaluate failed: %s\n", qaws_status_to_string(status));
		qaws_curve_destroy(curve);
		return 1;
	}

	printf("bezier midpoint: %.4f %.4f %.4f\n",
		(double)result.position.x,
		(double)result.position.y,
		(double)result.position.z);

	qaws_curve_destroy(curve);

	/* The midpoint of this control polygon is exactly (1.5, 1.5, 0). */
	if (result.position.x < (qaws_scalar)1.49 || result.position.x > (qaws_scalar)1.51 ||
	    result.position.y < (qaws_scalar)1.49 || result.position.y > (qaws_scalar)1.51)
	{
		printf("unexpected result\n");
		return 1;
	}

	printf("amalgamation OK\n");
	return 0;
}
