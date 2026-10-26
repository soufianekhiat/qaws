/*
 * Test 179: motions (qaws_motion.h)
 *
 *   - two poses: each end the given pose (position, rotation, velocity,
 *     angular velocity), the rotation orthonormal throughout, its angular
 *     velocity the finite differences' of R(t);
 *   - quintic: the accelerations and angular accelerations matched too;
 *   - a cubic path given by its ends' values and derivatives reproduced;
 *   - along a curve: its points and derivatives; per vertex: each vertex's
 *     cubic reproduced;
 *   - a half turn refused.
 */

#include "test_common.h"
#include "qaws_motion.h"
#include "qaws_bezier.h"
#include "qaws_eval.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

static qaws_vec3 mt_v(double x, double y, double z)
{
	qaws_vec3 v;
	v.x = (qaws_scalar)x;
	v.y = (qaws_scalar)y;
	v.z = (qaws_scalar)z;
	return v;
}

static void mt_quat(qaws_scalar* q, double ax, double ay, double az, double angle)
{
	double l = sqrt(ax * ax + ay * ay + az * az), s = sin(0.5 * angle) / l;
	q[0] = (qaws_scalar)cos(0.5 * angle);
	q[1] = (qaws_scalar)(ax * s);
	q[2] = (qaws_scalar)(ay * s);
	q[3] = (qaws_scalar)(az * s);
}

static double mt_quat_err(qaws_scalar const* q, qaws_scalar const* R)
{
	double w = q[0], x = q[1], y = q[2], z = q[3], l = sqrt(w * w + x * x + y * y + z * z), M[9], e = 0;
	unsigned int k;
	w /= l; x /= l; y /= l; z /= l;
	M[0] = 1 - 2 * (y * y + z * z); M[1] = 2 * (x * y - w * z);     M[2] = 2 * (x * z + w * y);
	M[3] = 2 * (x * y + w * z);     M[4] = 1 - 2 * (x * x + z * z); M[5] = 2 * (y * z - w * x);
	M[6] = 2 * (x * z - w * y);     M[7] = 2 * (y * z + w * x);     M[8] = 1 - 2 * (x * x + y * y);
	for (k = 0; k < 9; k++) e = fmax(e, fabs(M[k] - (double)R[k]));
	return e;
}

static double mt_d3(qaws_vec3 a, qaws_vec3 b)
{
	return fmax(fabs((double)a.x - b.x), fmax(fabs((double)a.y - b.y), fabs((double)a.z - b.z)));
}

/* the angular velocity from finite differences of R at t: [w]x = R' R^T */
static void mt_fd_omega(qaws_motion const* m, double t, double h, double* w)
{
	qaws_motion_state a, b, c;
	double dR[9], A[9];
	unsigned int i, j, k;
	qaws_motion_evaluate(&a, m, (qaws_scalar)(t - h));
	qaws_motion_evaluate(&b, m, (qaws_scalar)(t + h));
	qaws_motion_evaluate(&c, m, (qaws_scalar)t);
	for (k = 0; k < 9; k++) dR[k] = ((double)b.rotation[k] - (double)a.rotation[k]) / (2 * h);
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
		{
			A[3 * i + j] = 0;
			for (k = 0; k < 3; k++) A[3 * i + j] += dR[3 * i + k] * (double)c.rotation[3 * j + k];
		}
	w[0] = 0.5 * (A[7] - A[5]);
	w[1] = 0.5 * (A[2] - A[6]);
	w[2] = 0.5 * (A[3] - A[1]);
}

static void test_poses(void)
{
	qaws_motion_pose s, e;
	qaws_motion* m = NULL;
	qaws_status st;
	double ends = 0, orth = 0, fdw = 0, fdv = 0, acc = 0;
	unsigned int i, j, k, z;
	char msg[400];
	memset(&s, 0, sizeof s);
	memset(&e, 0, sizeof e);
	s.position = mt_v(0.1, -0.2, 0.3);
	e.position = mt_v(1.5, 0.4, -0.6);
	mt_quat(s.rotation, 1, 2, 0.5, 0.4);
	mt_quat(e.rotation, -0.3, 1, 1, 2.1);
	s.velocity = mt_v(0.5, 1, -0.2);
	e.velocity = mt_v(-0.4, 0.3, 0.9);
	s.angular_velocity = mt_v(0.7, -0.3, 0.2);
	e.angular_velocity = mt_v(-0.2, 0.5, 1.1);
	s.acceleration = mt_v(0.3, -1, 0.4);
	e.acceleration = mt_v(-0.6, 0.2, 0.5);
	s.angular_acceleration = mt_v(0.4, 0.1, -0.5);
	e.angular_acceleration = mt_v(-0.3, 0.6, 0.2);
	for (z = 0; z < 2; z++)
	{
		st = qaws_motion_create_poses(&m, &s, &e, z ? QAWS_MOTION_QUINTIC : 0);
		if (st != QAWS_STATUS_OK)
			break;
		for (k = 0; k < 2; k++)
		{
			qaws_motion_state q;
			qaws_motion_pose const* p = k ? &e : &s;
			qaws_motion_evaluate(&q, m, (qaws_scalar)k);
			ends = fmax(ends, fmax(mt_d3(q.position, p->position), mt_quat_err(p->rotation, q.rotation)));
			ends = fmax(ends, fmax(mt_d3(q.velocity, p->velocity), mt_d3(q.angular_velocity, p->angular_velocity)));
			if (z)
			{
				/* the second derivatives by central differences at the ends (one-sided past them: the polynomials go on) */
				qaws_motion_state a, b;
				double h = 1e-4, wa[3], wb[3];
				qaws_motion_evaluate(&a, m, (qaws_scalar)(k - h));
				qaws_motion_evaluate(&b, m, (qaws_scalar)(k + h));
				acc = fmax(acc, fabs(((double)b.velocity.x - a.velocity.x) / (2 * h) - p->acceleration.x));
				acc = fmax(acc, fabs(((double)b.velocity.y - a.velocity.y) / (2 * h) - p->acceleration.y));
				acc = fmax(acc, fabs(((double)b.velocity.z - a.velocity.z) / (2 * h) - p->acceleration.z));
				mt_fd_omega(m, k - h, 1e-6, wa);
				mt_fd_omega(m, k + h, 1e-6, wb);
				acc = fmax(acc, fabs((wb[0] - wa[0]) / (2 * h) - p->angular_acceleration.x));
				acc = fmax(acc, fabs((wb[1] - wa[1]) / (2 * h) - p->angular_acceleration.y));
				acc = fmax(acc, fabs((wb[2] - wa[2]) / (2 * h) - p->angular_acceleration.z));
			}
		}
		for (k = 1; k < 10; k++)
		{
			qaws_motion_state q, a, b;
			double t = k / 10.0, w[3], h = 1e-6;
			qaws_motion_evaluate(&q, m, (qaws_scalar)t);
			for (i = 0; i < 3; i++)
				for (j = 0; j < 3; j++)
				{
					double d = 0;
					unsigned int r;
					for (r = 0; r < 3; r++) d += (double)q.rotation[3 * r + i] * (double)q.rotation[3 * r + j];
					orth = fmax(orth, fabs(d - (i == j)));
				}
			mt_fd_omega(m, t, h, w);
			fdw = fmax(fdw, fmax(fabs(w[0] - q.angular_velocity.x), fmax(fabs(w[1] - q.angular_velocity.y), fabs(w[2] - q.angular_velocity.z))));
			qaws_motion_evaluate(&a, m, (qaws_scalar)(t - h));
			qaws_motion_evaluate(&b, m, (qaws_scalar)(t + h));
			fdv = fmax(fdv, fabs(((double)b.position.x - a.position.x) / (2 * h) - q.velocity.x));
		}
		qaws_motion_destroy(m);
		m = NULL;
	}
	sprintf(msg, "two poses, cubic and quintic: ends off the poses by %.1e, R(t) off orthonormal %.1e, the angular velocity off R's finite differences %.1e, the velocity %.1e; quintic accelerations off %.1e",
		ends, orth, fdw, fdv, acc);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && ends < 1e-14 && orth < 1e-14 && fdw < 1e-8 && fdv < 1e-8 && acc < 1e-5, msg);
}

/* a cubic path p(t) = (t^3 - t, 2 t^2 - 1, 0.5 t^3 + t) from its ends' values and derivatives: exact at every t; and a
   vertex motion of three points on cubics; a curve's points */
static void test_reproduction(void)
{
	qaws_motion_pose s, e;
	qaws_motion *m = NULL, *v = NULL, *c = NULL;
	qaws_curve* cv = NULL;
	qaws_scalar p0[9], p1[9], v0[9], v1[9], out[9], outv[9];
	double cub = 0, vert = 0, crv = 0;
	unsigned int k, i, j;
	qaws_status st, sv, sc;
	char msg[300];
	memset(&s, 0, sizeof s);
	memset(&e, 0, sizeof e);
	s.position = mt_v(0, -1, 0);           e.position = mt_v(0, 1, 1.5);
	s.velocity = mt_v(-1, 0, 1);           e.velocity = mt_v(2, 4, 2.5);
	s.rotation[0] = e.rotation[0] = 1;
	st = qaws_motion_create_poses(&m, &s, &e, 0);
	for (k = 0; st == QAWS_STATUS_OK && k <= 16; k++)
	{
		qaws_motion_state q;
		double t = k / 16.0;
		qaws_motion_evaluate(&q, m, (qaws_scalar)t);
		cub = fmax(cub, mt_d3(q.position, mt_v(t * t * t - t, 2 * t * t - 1, 0.5 * t * t * t + t)));
	}
	/* three vertices, vertex i on (t^3 + i, t^2 - i t, i t) */
	for (i = 0; i < 3; i++)
	{
		p0[3 * i] = (qaws_scalar)i; p0[3 * i + 1] = 0; p0[3 * i + 2] = 0;
		p1[3 * i] = (qaws_scalar)(1 + i); p1[3 * i + 1] = (qaws_scalar)(1 - (double)i); p1[3 * i + 2] = (qaws_scalar)i;
		v0[3 * i] = 0; v0[3 * i + 1] = -(qaws_scalar)i; v0[3 * i + 2] = (qaws_scalar)i;
		v1[3 * i] = 3; v1[3 * i + 1] = (qaws_scalar)(2 - (double)i); v1[3 * i + 2] = (qaws_scalar)i;
	}
	sv = qaws_motion_create_vertices(&v, p0, p1, v0, v1, NULL, NULL, 3, 0);
	for (k = 0; sv == QAWS_STATUS_OK && k <= 16; k++)
	{
		double t = k / 16.0;
		qaws_motion_transform_points(out, outv, v, (qaws_scalar)t, NULL, 3);
		for (i = 0; i < 3; i++)
		{
			double want[3] = { t * t * t + i, t * t - i * t, (double)i * t };
			for (j = 0; j < 3; j++) vert = fmax(vert, fabs((double)out[3 * i + j] - want[j]));
		}
	}
	{
		/* a cubic Bezier in space, its points and derivatives through the motion */
		qaws_scalar cp[12] = { 0, 0, 0, 1, 2, 0, 2, -1, 1, 3, 0, 2 };
		qaws_bezier_desc bd;
		qaws_scalar rot[4] = { 1, 0, 0, 0 };
		memset(&bd, 0, sizeof bd);
		bd.dimension = QAWS_DIMENSION_3D; bd.degree = 3; bd.control_points = cp; bd.control_point_count = 4;
		qaws_curve_create_bezier(&cv, &bd);
		sc = qaws_motion_create_curve(&c, cv, rot);
		for (k = 0; sc == QAWS_STATUS_OK && k <= 8; k++)
		{
			qaws_motion_state q;
			qaws_eval_result_3d r;
			double t = k / 8.0;
			qaws_motion_evaluate(&q, c, (qaws_scalar)t);
			qaws_curve_evaluate_3d(&r, cv, (qaws_scalar)t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1);
			crv = fmax(crv, fmax(mt_d3(q.position, r.position), mt_d3(q.velocity, r.d1)));
		}
	}
	sprintf(msg, "reproduction: a cubic path from its ends off by %.1e at 17 times; three vertices' cubics off by %.1e; a Bezier curve's points and derivatives off by %.1e",
		cub, vert, crv);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && sv == QAWS_STATUS_OK && sc == QAWS_STATUS_OK && cub < 1e-14 && vert < 1e-14 && crv < 1e-14, msg);
	qaws_motion_destroy(m);
	qaws_motion_destroy(v);
	qaws_motion_destroy(c);
	qaws_curve_destroy(cv);
}

static void test_refusals(void)
{
	qaws_motion_pose s, e;
	qaws_motion* m = NULL;
	qaws_status a, b;
	char msg[200];
	memset(&s, 0, sizeof s);
	memset(&e, 0, sizeof e);
	s.rotation[0] = 1;
	mt_quat(e.rotation, 0, 0, 1, 3.14159265358979323846);
	a = qaws_motion_create_poses(&m, &s, &e, 0);
	e.rotation[0] = e.rotation[1] = e.rotation[2] = e.rotation[3] = 0;
	b = qaws_motion_create_poses(&m, &s, &e, 0);
	sprintf(msg, "refused: a half turn (%d), a zero quaternion (%d)", (int)a, (int)b);
	printf("  %s\n", msg);
	TEST_ASSERT(a == QAWS_STATUS_INVALID_ARGUMENT && b == QAWS_STATUS_INVALID_ARGUMENT && !m, msg);
}

int test_179_motion_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 179: motions\n");
	test_poses();
	test_reproduction();
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
