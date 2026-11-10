/*
 * Interpolating splines solved from their points.
 *
 * The C2, Hobby and kappa solvers follow aurimasg/splines:
 *
 * MIT License
 *
 * Copyright (c) 2026 Aurimas Gasiulis
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "qaws_interpolating.h"
#include "qaws_alloc.h"
#include "qaws_bezier.h"
#include "qaws_composite.h"
#include "qaws_curve.h"
#include "internal/qaws_internal_bezier_chain.h"
#include "internal/qaws_internal_memory.h"
#include <float.h>
#include <math.h>
#include <string.h>

#define IS_PI 3.14159265358979323846

qaws_status qaws_internal_bezier_chain(qaws_curve** out_curve, double const* points, unsigned int degree, unsigned int segments,
	qaws_allocator const* allocator)
{
	qaws_curve** seg = (qaws_curve**)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_curve*), (size_t)segments + 1));
	qaws_status st = seg ? QAWS_STATUS_OK : QAWS_STATUS_ALLOCATION_FAILURE;
	unsigned int i, j, made = 0;
	for (i = 0; st == QAWS_STATUS_OK && i < segments; i++)
	{
		qaws_vec2 cp[4];
		qaws_bezier_desc d;
		for (j = 0; j <= degree; j++)
		{
			cp[j].x = (qaws_scalar)points[2 * ((size_t)degree * i + j)];
			cp[j].y = (qaws_scalar)points[2 * ((size_t)degree * i + j) + 1];
		}
		memset(&d, 0, sizeof d);
		d.dimension = QAWS_DIMENSION_2D;
		d.degree = degree;
		d.control_points = cp;
		d.control_point_count = degree + 1;
		st = qaws_curve_create_bezier_ex(&seg[i], &d, allocator);
		made += st == QAWS_STATUS_OK;
	}
	if (st == QAWS_STATUS_OK)
	{
		qaws_composite_desc c;
		memset(&c, 0, sizeof c);
		c.dimension = QAWS_DIMENSION_2D;
		c.segments = seg;
		c.segment_count = segments;
		st = qaws_curve_create_composite_ex(out_curve, &c, allocator);
	}
	if (st != QAWS_STATUS_OK)
		for (i = 0; i < made; i++)
			qaws_curve_destroy(seg[i]);
	qaws_internal_free(seg);
	return st;
}

/* a x[i-1] + b x[i] + c x[i+1] = d; w: 2 n doubles. x may be d. */
static void is_tridiagonal(unsigned int n, double const* a, double const* b, double const* c, double const* d, double* x, double* w)
{
	double* sup = w;
	double* rhs = w + n;
	unsigned int i;
	sup[0] = c[0] / b[0];
	rhs[0] = d[0] / b[0];
	for (i = 1; i < n; i++)
	{
		double den = b[i] - a[i] * sup[i - 1];
		sup[i] = c[i] / den;
		rhs[i] = (d[i] - a[i] * rhs[i - 1]) / den;
	}
	x[n - 1] = rhs[n - 1];
	for (i = n - 1; i > 0; i--)
		x[i - 1] = rhs[i - 1] - sup[i - 1] * x[i];
}

/* The same with a[0] coupling row 0 to x[n-1] and c[n-1] coupling row n-1 to x[0] (n >= 3): the two links are a rank one
   update of a tridiagonal matrix (Sherman-Morrison). w: 4 n doubles. */
static void is_cyclic(unsigned int n, double const* a, double const* b, double const* c, double const* d, double* x, double* w)
{
	double* u = w;
	double* dg = w + n;
	double* work = w + 2 * (size_t)n;
	double factor;
	unsigned int i, last = n - 1;
	memcpy(dg, b, sizeof(double) * n);
	dg[0] -= a[0];
	dg[last] -= c[last];
	memset(u, 0, sizeof(double) * n);
	u[0] = a[0];
	u[last] = c[last];
	is_tridiagonal(n, a, dg, c, u, u, work);
	is_tridiagonal(n, a, dg, c, d, x, work);
	factor = (x[0] + x[last]) / (1.0 + u[0] + u[last]);
	for (i = 0; i < n; i++)
		x[i] -= factor * u[i];
}

/* the C2 cubics through k (n points): 3 m + 1 points in out, m the segment count */
static qaws_status is_cubic_c2(double const* k, unsigned int n, int closed, double* out, unsigned int* out_segments)
{
	unsigned int r = closed ? n : n - 1, last = r - 1, i, axis;
	double* mem = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 9, (size_t)r + 1));
	double *a, *b, *c, *d, *w;
	if (!mem)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	a = mem;
	b = a + r;
	c = b + r;
	d = c + r;
	w = d + r;
	for (i = 0; i < r; i++)
	{
		a[i] = 1.0;
		b[i] = 4.0;
		c[i] = 1.0;
	}
	if (!closed)
	{
		a[0] = 0.0;
		b[0] = 2.0;
		a[last] = 2.0;
		b[last] = 7.0;
		c[last] = 0.0;
	}
	for (axis = 0; axis < 2; axis++)
	{
		for (i = 0; i < r; i++)
			d[i] = 4.0 * k[2 * i + axis] + 2.0 * k[2 * ((i + 1) % n) + axis];
		if (!closed)
		{
			d[0] = k[axis] + 2.0 * k[2 + axis];
			d[last] = 8.0 * k[2 * last + axis] + k[2 * (size_t)r + axis];
		}
		if (closed)
			is_cyclic(r, a, b, c, d, d, w);
		else if (r == 1)
			d[0] = (k[axis] + 2.0 * k[2 + axis]) / 3.0;
		else
			is_tridiagonal(r, a, b, c, d, d, w);
		for (i = 0; i < r; i++)
		{
			size_t at = 3 * (size_t)i;
			unsigned int j = (i + 1) % n;
			out[2 * at + axis] = k[2 * i + axis];
			out[2 * (at + 1) + axis] = d[i];
			if (i < last || closed)
				out[2 * (at + 2) + axis] = 2.0 * k[2 * j + axis] - d[(i + 1) % r];
			else
				out[2 * (at + 2) + axis] = (k[2 * (size_t)r + axis] + d[last]) / 2.0;
			out[2 * (at + 3) + axis] = k[2 * j + axis];
		}
	}
	*out_segments = r;
	qaws_internal_free(mem);
	return QAWS_STATUS_OK;
}

/* Hobby's velocity of the end leaving at angle a toward one arriving at b */
static double is_velocity(double a, double b)
{
	double sa = sin(a), sb = sin(b), ca = cos(a), cb = cos(b), s5 = sqrt(5.0);
	double num = 4.0 + sqrt(8.0) * (sa - sb / 16.0) * (sb - sa / 16.0) * (ca - cb);
	double den = 2.0 + (s5 - 1.0) * ca + (3.0 - s5) * cb;
	return num / den;
}

static void is_rotate_unit(double* ox, double* oy, double x, double y, double angle)
{
	double s = sin(angle), c = cos(angle), rx = x * c - y * s, ry = x * s + y * c, len = sqrt(rx * rx + ry * ry);
	*ox = rx / len;
	*oy = ry / len;
}

static qaws_status is_hobby(double const* k, unsigned int count, int closed, double* out, unsigned int* out_segments)
{
	unsigned int n = closed ? count : count - 1, unknowns = closed ? count : count, i;
	double* mem = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 16, (size_t)count + 2));
	double *D, *dx, *dy, *gamma, *a, *b, *c, *d, *alpha, *beta, *w;
	if (!mem)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	D = mem;
	dx = D + count + 1;
	dy = dx + count + 1;
	gamma = dy + count + 1;
	a = gamma + count + 1;
	b = a + count + 1;
	c = b + count + 1;
	d = c + count + 1;
	alpha = d + count + 1;
	beta = alpha + count + 1;
	w = beta + count + 1;
	for (i = 0; i < n; i++)
	{
		unsigned int j = (i + 1) % count;
		dx[i] = k[2 * j] - k[2 * i];
		dy[i] = k[2 * j + 1] - k[2 * i + 1];
		D[i] = sqrt(dx[i] * dx[i] + dy[i] * dy[i]);
	}
	/* the turn at each point, between the chords on either side of it */
	gamma[0] = 0;
	gamma[n] = 0;
	for (i = closed ? 0 : 1; i < n; i++)
	{
		unsigned int p = closed ? (i + n - 1) % n : i - 1;
		double sn = dy[p] / D[p], cs = dx[p] / D[p];
		gamma[i] = atan2(dy[i] * cs - dx[i] * sn, dx[i] * cs + dy[i] * sn);
	}
	if (closed)
	{
		for (i = 0; i < n; i++)
		{
			unsigned int j = (i + 1) % n, p = (i + n - 1) % n;
			a[i] = 1.0 / D[p];
			b[i] = (2.0 * D[p] + 2.0 * D[i]) / (D[p] * D[i]);
			c[i] = 1.0 / D[i];
			d[i] = -(2.0 * gamma[i] * D[i] + gamma[j] * D[p]) / (D[p] * D[i]);
		}
		is_cyclic(n, a, b, c, d, alpha, w);
		for (i = 0; i < n; i++)
		{
			unsigned int j = (i + 1) % n;
			beta[i] = -gamma[j] - alpha[j];
		}
	}
	else
	{
		/* (the ends' curl at its usual value 1) */
		a[0] = 0.0;
		b[0] = 3.0;
		c[0] = 3.0;
		d[0] = -c[0] * gamma[1];
		for (i = 1; i < n; i++)
		{
			a[i] = 1.0 / D[i - 1];
			b[i] = (2.0 * D[i - 1] + 2.0 * D[i]) / (D[i - 1] * D[i]);
			c[i] = 1.0 / D[i];
			d[i] = -(2.0 * gamma[i] * D[i] + gamma[i + 1] * D[i - 1]) / (D[i - 1] * D[i]);
		}
		a[n] = 3.0;
		b[n] = 3.0;
		c[n] = 0.0;
		d[n] = 0.0;
		is_tridiagonal(unknowns, a, b, c, d, alpha, w);
		for (i = 0; i + 1 < n; i++)
			beta[i] = -gamma[i + 1] - alpha[i + 1];
		beta[n - 1] = -alpha[n];
	}
	for (i = 0; i < n; i++)
	{
		unsigned int j = (i + 1) % count;
		size_t at = 3 * (size_t)i;
		double ra = is_velocity(alpha[i], beta[i]) * D[i] / 3, rb = is_velocity(beta[i], alpha[i]) * D[i] / 3, ux, uy, vx, vy;
		is_rotate_unit(&ux, &uy, dx[i], dy[i], alpha[i]);
		is_rotate_unit(&vx, &vy, dx[i], dy[i], -beta[i]);
		out[2 * at] = k[2 * i];
		out[2 * at + 1] = k[2 * i + 1];
		out[2 * (at + 1)] = k[2 * i] + ux * ra;
		out[2 * (at + 1) + 1] = k[2 * i + 1] + uy * ra;
		out[2 * (at + 2)] = k[2 * j] - vx * rb;
		out[2 * (at + 2) + 1] = k[2 * j + 1] - vy * rb;
		out[2 * (at + 3)] = k[2 * j];
		out[2 * (at + 3) + 1] = k[2 * j + 1];
	}
	*out_segments = n;
	qaws_internal_free(mem);
	return QAWS_STATUS_OK;
}

/* a root of the quadratic in [0, 1] when there is one */
static double is_quadratic_root(double a, double b, double c)
{
	double delta = b * b - 4.0 * a * c;
	if (delta < 0.0)
		return 0.0;
	if (delta > 0.0)
	{
		double s = sqrt(delta), q = -0.5 * (b + (b < 0.0 ? -s : s)), r0 = q / a, r1 = c / q;
		return r0 >= 0.0 && r0 <= 1.0 ? r0 : r1;
	}
	return a != 0.0 ? -0.5 * b / a : 0.5;
}

/* a root of the cubic in [0, 1] when there is one (Numerical Recipes 5.6) */
static double is_cubic_root(double c0, double c1, double c2, double c3)
{
	double a, b, c, Q, R, disc, third, A;
	if (fabs(c0) < DBL_EPSILON)
		return is_quadratic_root(c1, c2, c3);
	a = c1 / c0;
	b = c2 / c0;
	c = c3 / c0;
	Q = (a * a - b * 3.0) / 9.0;
	R = (2.0 * a * a * a - 9.0 * a * b + 27.0 * c) / 54.0;
	disc = R * R - Q * Q * Q;
	third = a / 3.0;
	if (disc < 0.0)
	{
		double ratio = R / sqrt(Q * Q * Q), theta = acos(ratio < -1.0 ? -1.0 : ratio > 1.0 ? 1.0 : ratio), m = -2.0 * sqrt(Q);
		double x0 = m * cos(theta / 3.0) - third, x1 = m * cos((theta + 2.0 * IS_PI) / 3.0) - third, x2 = m * cos((theta - 2.0 * IS_PI) / 3.0) - third;
		if (x0 >= 0.0 && x0 <= 1.0)
			return x0;
		if (x1 >= 0.0 && x1 <= 1.0)
			return x1;
		if (x2 >= 0.0 && x2 <= 1.0)
			return x2;
	}
	A = pow(fabs(R) + sqrt(disc > 0.0 ? disc : 0.0), 1.0 / 3.0);
	if (R > 0.0)
		A = -A;
	if (A != 0.0)
		A += Q / A;
	return A - third;
}

/* the parameter of a quadratic's curvature maximum, held off 0 and 1 where the point would not be reached */
static double is_max_param(double const* c0, double const* cp, double const* c2)
{
	double ex = c2[0] - c0[0], ey = c2[1] - c0[1], gx = c0[0] - cp[0], gy = c0[1] - cp[1];
	double hx = 3.0 * c0[0] - 2.0 * cp[0] - c2[0], hy = 3.0 * c0[1] - 2.0 * cp[1] - c2[1];
	double t = is_cubic_root(ex * ex + ey * ey, 3.0 * (ex * gx + ey * gy), hx * gx + hy * gy, -(gx * gx + gy * gy));
	return !(t >= 1e-6) ? 1e-6 : t > 1.0 - 1e-6 ? 1.0 - 1e-6 : t;
}

static double is_area(double const* a, double const* b, double const* c)
{
	return fabs(a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1])) / 2.0;
}

/* where the join between two quadratics lies along their control points, by the triangles either side of it (halfway
   when both are flat) */
static double is_lambda(double const* c0, double const* c1, double const* c3, double const* c4)
{
	double l = sqrt(is_area(c0, c1, c3)), r = sqrt(is_area(c1, c3, c4));
	return l + r > 0.0 ? l / (l + r) : 0.5;
}

/* quadratics through k: 2 m + 1 points in out. Each quadratic's control point moves until its curvature maximum falls on
   its point and the joins are G2. */
static qaws_status is_kappa(double const* k, unsigned int count, int closed, double* out, unsigned int* out_segments)
{
	enum { ITERATIONS = 60 };
	unsigned int n = closed ? count : count - 2, i, it, axis;
	double const* mid = closed ? k : k + 2;   /* the point each quadratic passes through */
	double* mem = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 18, (size_t)n + 1));
	double *cp, *start, *end, *t, *l, *A, *B, *C, *rhs, *w;
	if (!mem)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	cp = mem;
	start = cp + 2 * (size_t)n;
	end = start + 2 * (size_t)n;
	t = end + 2 * (size_t)n;
	l = t + n;
	A = l + n;
	B = A + n;
	C = B + n;
	rhs = C + n;
	w = rhs + n;
	for (i = 0; i < n; i++)
	{
		unsigned int p = closed ? (i + n - 1) % n : i, q = closed ? (i + 1) % n : i + 2;
		cp[2 * i] = mid[2 * i];
		cp[2 * i + 1] = mid[2 * i + 1];
		for (axis = 0; axis < 2; axis++)
		{
			double before = closed ? k[2 * p + axis] : k[2 * p + axis], after = k[2 * q + axis];
			start[2 * i + axis] = (mid[2 * i + axis] + before) / 2.0;
			end[2 * i + axis] = (mid[2 * i + axis] + after) / 2.0;
		}
	}
	if (!closed)
	{
		start[0] = k[0];
		start[1] = k[1];
		end[2 * ((size_t)n - 1)] = k[2 * ((size_t)n + 1)];
		end[2 * ((size_t)n - 1) + 1] = k[2 * ((size_t)n + 1) + 1];
	}
	for (it = 0; it < ITERATIONS; it++)
	{
		double moved = 0.0;
		unsigned int joins = closed ? n : n - 1;
		for (i = 0; i < n; i++)
			t[i] = is_max_param(start + 2 * (size_t)i, mid + 2 * (size_t)i, end + 2 * (size_t)i);
		for (i = 0; i < joins; i++)
		{
			unsigned int j = (i + 1) % n;
			l[i] = is_lambda(start + 2 * (size_t)i, cp + 2 * (size_t)i, cp + 2 * (size_t)j, end + 2 * (size_t)j);
		}
		for (i = 0; i < n; i++)
		{
			int has_prev = closed || i > 0, has_next = closed || i + 1 < n;
			double lp = has_prev ? l[closed ? (i + n - 1) % n : i - 1] : 0.0, ln = has_next ? l[i] : 0.0;
			double s = 1.0 - t[i];
			A[i] = has_prev ? (1.0 - lp) * s * s : 0.0;
			B[i] = (has_prev ? lp * s * s : 0.0) + 2.0 * s * t[i] + (has_next ? (1.0 - ln) * t[i] * t[i] : 0.0);
			C[i] = has_next ? ln * t[i] * t[i] : 0.0;
		}
		for (axis = 0; axis < 2; axis++)
		{
			for (i = 0; i < n; i++)
				rhs[i] = mid[2 * i + axis];
			if (!closed)
			{
				rhs[0] -= (1.0 - t[0]) * (1.0 - t[0]) * k[axis];
				rhs[n - 1] -= t[n - 1] * t[n - 1] * k[2 * ((size_t)n + 1) + axis];
			}
			if (n == 1)
				rhs[0] /= B[0];
			else if (closed)
				is_cyclic(n, A, B, C, rhs, rhs, w);
			else
				is_tridiagonal(n, A, B, C, rhs, rhs, w);
			for (i = 0; i < n; i++)
				cp[2 * i + axis] = rhs[i];
		}
		for (i = 0; i < joins; i++)
		{
			unsigned int j = (i + 1) % n;
			double x = (1.0 - l[i]) * cp[2 * i] + l[i] * cp[2 * j], y = (1.0 - l[i]) * cp[2 * i + 1] + l[i] * cp[2 * j + 1];
			double ex = x - start[2 * j], ey = y - start[2 * j + 1], m2 = ex * ex + ey * ey;
			if (m2 > moved)
				moved = m2;
			start[2 * j] = end[2 * i] = x;
			start[2 * j + 1] = end[2 * i + 1] = y;
		}
		if (moved <= 1e-16 && it > 0)
			break;
	}
	for (i = 0; i < n; i++)
	{
		size_t at = 2 * (size_t)i;
		out[2 * at] = start[2 * i];
		out[2 * at + 1] = start[2 * i + 1];
		out[2 * (at + 1)] = cp[2 * i];
		out[2 * (at + 1) + 1] = cp[2 * i + 1];
		out[2 * (at + 2)] = end[2 * i];
		out[2 * (at + 2) + 1] = end[2 * i + 1];
	}
	*out_segments = n;
	qaws_internal_free(mem);
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_create_interpolating_spline(
	qaws_curve** out_curve,
	qaws_interpolating_spline_desc const* desc)
{
	double *k = NULL, *unit = NULL, *out = NULL;
	double lox = HUGE_VAL, loy = HUGE_VAL, hix = -HUGE_VAL, hiy = -HUGE_VAL, tx, ty, m;
	unsigned int i, n = 0, segments = 0, degree;
	size_t points;
	qaws_status st;
	if (!out_curve || !desc || !desc->points || desc->point_count < 1)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_curve = NULL;
	if (desc->method != QAWS_INTERPOLATING_CUBIC_C2 && desc->method != QAWS_INTERPOLATING_HOBBY && desc->method != QAWS_INTERPOLATING_KAPPA)
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < desc->point_count; i++)
		if (!(fabs((double)desc->points[i].x) <= DBL_MAX) || !(fabs((double)desc->points[i].y) <= DBL_MAX))
			return QAWS_STATUS_INVALID_ARGUMENT;
	k = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 2, (size_t)desc->point_count + 1));
	if (!k)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < desc->point_count; i++)
	{
		double x = (double)desc->points[i].x, y = (double)desc->points[i].y;
		if (n && k[2 * ((size_t)n - 1)] == x && k[2 * ((size_t)n - 1) + 1] == y)
			continue;
		k[2 * (size_t)n] = x;
		k[2 * (size_t)n + 1] = y;
		n++;
		lox = x < lox ? x : lox;
		hix = x > hix ? x : hix;
		loy = y < loy ? y : loy;
		hiy = y > hiy ? y : hiy;
	}
	if (desc->closed)
		while (n > 1 && k[0] == k[2 * ((size_t)n - 1)] && k[1] == k[2 * ((size_t)n - 1) + 1])
			n--;
	if (n < 2 || (desc->closed && n < 3))
	{
		qaws_internal_free(k);
		return QAWS_STATUS_DEGENERATE_CURVE;
	}
	/* (the solvers' stopping rules are absolute: they work on the points scaled into [-1, 1]) */
	tx = (lox + hix) / 2.0;
	ty = (loy + hiy) / 2.0;
	m = (hix - lox > hiy - loy ? hix - lox : hiy - loy) / 2.0;
	degree = desc->method == QAWS_INTERPOLATING_KAPPA && (n > 2 || desc->closed) ? 2 : 3;
	points = (size_t)degree * ((size_t)n + 1) + 1;
	unit = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 2, (size_t)n + 1));
	out = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 2, points));
	st = unit && out ? QAWS_STATUS_OK : QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; st == QAWS_STATUS_OK && i < n; i++)
	{
		unit[2 * i] = (k[2 * i] - tx) / m;
		unit[2 * i + 1] = (k[2 * i + 1] - ty) / m;
	}
	if (st == QAWS_STATUS_OK && n == 2 && !desc->closed)
	{
		for (i = 0; i < 4; i++)
		{
			double s = (double)i / 3.0;
			out[2 * i] = unit[0] + s * (unit[2] - unit[0]);
			out[2 * i + 1] = unit[1] + s * (unit[3] - unit[1]);
		}
		segments = 1;
	}
	else if (st == QAWS_STATUS_OK && desc->method == QAWS_INTERPOLATING_CUBIC_C2)
		st = is_cubic_c2(unit, n, desc->closed != 0, out, &segments);
	else if (st == QAWS_STATUS_OK && desc->method == QAWS_INTERPOLATING_HOBBY)
		st = is_hobby(unit, n, desc->closed != 0, out, &segments);
	else if (st == QAWS_STATUS_OK)
		st = is_kappa(unit, n, desc->closed != 0, out, &segments);
	if (st == QAWS_STATUS_OK)
	{
		size_t count = 2 * ((size_t)degree * segments + 1), q;
		for (q = 0; q < count; q++)
		{
			if (!(fabs(out[q]) <= DBL_MAX))
				st = QAWS_STATUS_NUMERICAL_FAILURE;
			out[q] = out[q] * m + ((q & 1) ? ty : tx);
		}
		/* (the points themselves back exactly, past the scaling's rounding) */
		if (degree == 3)
			for (i = 0; i <= segments; i++)
			{
				unsigned int j = i % n;
				out[2 * (3 * (size_t)i)] = k[2 * j];
				out[2 * (3 * (size_t)i) + 1] = k[2 * j + 1];
			}
	}
	if (st == QAWS_STATUS_OK)
		st = qaws_internal_bezier_chain(out_curve, out, degree, segments, NULL);
	qaws_internal_free(k);
	qaws_internal_free(unit);
	qaws_internal_free(out);
	return st;
}
