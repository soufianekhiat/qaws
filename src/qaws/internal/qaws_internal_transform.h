#ifndef QAWS_INTERNAL_TRANSFORM_H
#define QAWS_INTERNAL_TRANSFORM_H

/*
 * A qaws_transform in double, ready to map points and directions both ways: world = t + R (s * x).
 *   to_world   x -> t + R (s * x)
 *   to_frame   y -> (R^T (y - t)) / s   (the object's frame)
 *   to_turned  y -> R^T (y - t)         (the world undone of its rotation and translation: the object's frame scaled by
 *                                        s, where a box of the object's frame scaled per axis is still a box and
 *                                        distances are the world's)
 * Header-only: static inline, included where a query maps through a placement.
 */

#include "qaws_transform.h"
#include <math.h>

typedef struct qaws_tf
{
	double r[9];          /* R by rows */
	double s[3], inv[3];  /* scale, its inverse */
	double t[3];
	double det_sign;      /* the sign of the scale's product: -1 mirrors */
	double smin, smax;    /* the least and largest |s| */
	int identity;         /* nothing to do */
	int rotated;          /* R not the identity */
} qaws_tf;

/* 0 for a value not finite or a rotation of zero length with nonzero parts (never: all zero is none) */
static inline int qaws_tf_init(qaws_tf* f, qaws_transform const* x)
{
	double w = 1, a = 0, b = 0, c = 0, n;
	unsigned int k;
	if (!x)
	{
		for (k = 0; k < 9; k++) f->r[k] = k % 4 == 0;
		for (k = 0; k < 3; k++) f->s[k] = f->inv[k] = 1, f->t[k] = 0;
		f->det_sign = 1, f->smin = f->smax = 1, f->identity = 1, f->rotated = 0;
		return 1;
	}
	for (k = 0; k < 4; k++)
		if (!isfinite((double)x->rotation[k]))
			return 0;
	for (k = 0; k < 3; k++)
		if (!isfinite((double)x->scale[k]) || !isfinite((double)x->translation[k]))
			return 0;
	n = (double)x->rotation[0] * x->rotation[0] + (double)x->rotation[1] * x->rotation[1] + (double)x->rotation[2] * x->rotation[2]
		+ (double)x->rotation[3] * x->rotation[3];
	if (n > 0)
	{
		n = 1 / sqrt(n);
		w = x->rotation[0] * n, a = x->rotation[1] * n, b = x->rotation[2] * n, c = x->rotation[3] * n;
	}
	f->r[0] = 1 - 2 * (b * b + c * c), f->r[1] = 2 * (a * b - w * c), f->r[2] = 2 * (a * c + w * b);
	f->r[3] = 2 * (a * b + w * c), f->r[4] = 1 - 2 * (a * a + c * c), f->r[5] = 2 * (b * c - w * a);
	f->r[6] = 2 * (a * c - w * b), f->r[7] = 2 * (b * c + w * a), f->r[8] = 1 - 2 * (a * a + b * b);
	f->rotated = !(a == 0 && b == 0 && c == 0);
	f->det_sign = 1;
	f->smin = HUGE_VAL, f->smax = 0;
	for (k = 0; k < 3; k++)
	{
		double s = x->scale[k] != 0 ? (double)x->scale[k] : 1;
		f->s[k] = s, f->inv[k] = 1 / s, f->t[k] = x->translation[k];
		if (s < 0) f->det_sign = -f->det_sign;
		if (fabs(s) < f->smin) f->smin = fabs(s);
		if (fabs(s) > f->smax) f->smax = fabs(s);
	}
	f->identity = !f->rotated && f->s[0] == 1 && f->s[1] == 1 && f->s[2] == 1 && f->t[0] == 0 && f->t[1] == 0 && f->t[2] == 0;
	return 1;
}

/* world -> turned (R^T (y - t)) */
static inline void qaws_tf_to_turned(qaws_tf const* f, double const* y, double* o)
{
	double d[3] = { y[0] - f->t[0], y[1] - f->t[1], y[2] - f->t[2] };
	o[0] = f->r[0] * d[0] + f->r[3] * d[1] + f->r[6] * d[2];
	o[1] = f->r[1] * d[0] + f->r[4] * d[1] + f->r[7] * d[2];
	o[2] = f->r[2] * d[0] + f->r[5] * d[1] + f->r[8] * d[2];
}

/* world -> the object's frame */
static inline void qaws_tf_to_frame(qaws_tf const* f, double const* y, double* o)
{
	qaws_tf_to_turned(f, y, o);
	o[0] *= f->inv[0], o[1] *= f->inv[1], o[2] *= f->inv[2];
}

/* a direction, world -> the object's frame (no translation; not normalized: a ray's t is kept) */
static inline void qaws_tf_dir_to_frame(qaws_tf const* f, double const* d, double* o)
{
	o[0] = (f->r[0] * d[0] + f->r[3] * d[1] + f->r[6] * d[2]) * f->inv[0];
	o[1] = (f->r[1] * d[0] + f->r[4] * d[1] + f->r[7] * d[2]) * f->inv[1];
	o[2] = (f->r[2] * d[0] + f->r[5] * d[1] + f->r[8] * d[2]) * f->inv[2];
}

/* the object's frame -> turned (s * x) */
static inline void qaws_tf_frame_to_turned(qaws_tf const* f, double const* x, double* o)
{
	o[0] = f->s[0] * x[0], o[1] = f->s[1] * x[1], o[2] = f->s[2] * x[2];
}

/* turned -> world (t + R z) */
static inline void qaws_tf_turned_to_world(qaws_tf const* f, double const* z, double* o)
{
	o[0] = f->t[0] + f->r[0] * z[0] + f->r[1] * z[1] + f->r[2] * z[2];
	o[1] = f->t[1] + f->r[3] * z[0] + f->r[4] * z[1] + f->r[5] * z[2];
	o[2] = f->t[2] + f->r[6] * z[0] + f->r[7] * z[1] + f->r[8] * z[2];
}

/* the object's frame -> world */
static inline void qaws_tf_to_world(qaws_tf const* f, double const* x, double* o)
{
	double z[3];
	qaws_tf_frame_to_turned(f, x, z);
	qaws_tf_turned_to_world(f, z, o);
}

#endif /* QAWS_INTERNAL_TRANSFORM_H */
