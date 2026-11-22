/*
 * Test 181: placements (qaws_transform.h)
 *
 *   - a zeroed transform is the identity; a zero scale component counts as 1;
 *   - the order: the scale per axis of the object's frame, then the rotation, then the translation
 *     (a quarter turn about z of (2, 0, 0) scaled by (3, 1, 1) lands on (0, 6, 0) plus the translation);
 *   - points mapped and mapped back are the points (a turn about a skew axis, a mirror, a non-uniform scale);
 *   - refusals: no transform, a value not finite.
 */
#include "test_common.h"
#include "qaws_transform.h"
#include <math.h>
#include <string.h>

int test_181_transform_main(void)
{
	qaws_transform x;
	qaws_scalar p[6] = { 2, 0, 0, -1, 0.5, 3 }, q[6], r[6];
	double tol = QAWS_SCALAR_IS_FLOAT ? 1e-5 : 1e-13, err = 0;
	unsigned int i;
	char msg[256];
	g_pass = 0;
	g_fail = 0;
	printf("Test 181: placements (qaws_transform)\n");
	memset(&x, 0, sizeof x);
	TEST_ASSERT(qaws_transform_points(q, &x, p, 2) == QAWS_STATUS_OK && memcmp(p, q, sizeof p) == 0, "a zeroed transform is the identity");
	/* a quarter turn about z: w = cos 45, z = sin 45 */
	x.rotation[0] = (qaws_scalar)0.7071067811865476, x.rotation[3] = (qaws_scalar)0.7071067811865476;
	x.scale[0] = 3;
	x.translation[0] = 10, x.translation[1] = 20, x.translation[2] = 30;
	qaws_transform_points(q, &x, p, 1);
	sprintf(msg, "scale, then turn, then move: (2, 0, 0) to (%g, %g, %g), want (10, 26, 30)", (double)q[0], (double)q[1], (double)q[2]);
	TEST_ASSERT(fabs(q[0] - 10) < tol * 30 && fabs(q[1] - 26) < tol * 30 && fabs(q[2] - 30) < tol * 30, msg);
	/* a turn about a skew axis (not normalized), a mirror, a non-uniform scale: there and back */
	x.rotation[0] = 2, x.rotation[1] = 1, x.rotation[2] = -3, x.rotation[3] = (qaws_scalar)0.5;
	x.scale[0] = (qaws_scalar)0.25, x.scale[1] = -4, x.scale[2] = 7;
	qaws_transform_points(q, &x, p, 2);
	qaws_transform_points_inverse(r, &x, q, 2);
	for (i = 0; i < 6; i++)
		err = fmax(err, fabs((double)r[i] - p[i]));
	sprintf(msg, "there and back through a skew turn, a mirror and a non-uniform scale: off by %.2g", err);
	TEST_ASSERT(err < tol * 100, msg);
	TEST_ASSERT(qaws_transform_points(q, NULL, p, 2) == QAWS_STATUS_INVALID_ARGUMENT, "no transform refused");
	x.translation[1] = (qaws_scalar)NAN;
	TEST_ASSERT(qaws_transform_points(q, &x, p, 2) == QAWS_STATUS_INVALID_ARGUMENT, "a value not finite refused");
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
