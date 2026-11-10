/*
 * Parallel curves as cubic Bezier segments within a tolerance.
 *
 * Follows aurimasg/cubic-bezier-offsetter:
 *
 * MIT License
 *
 * Copyright (c) 2023 Aurimas Gasiulis
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

#include "qaws_operations.h"
#include "qaws_eval.h"
#include "internal/qaws_internal_bezier_chain.h"
#include "internal/qaws_internal_memory.h"
#include "internal/qaws_internal_types.h"
#include <float.h>
#include <math.h>
#include <string.h>

#define OC_PI 3.14159265358979323846
#define OC_NEAR_CUSP_ITERATIONS 18
#define OC_TINY_PERIMETER2 1e-7          /* a candidate no longer than this adds nothing */
#define OC_MIN_ARC_RADIUS 1e-8
#define OC_MAX_ARC_RADIUS 1e+6
#define OC_CUSP_DERIVATIVE2 1.5e-4       /* a squared derivative under this is a cusp (the curve scaled to about [-1, 1]) */
#define OC_POINT_EPS 1e-14
#define OC_ARC_CENTER_EPS 1e-8
#define OC_NEARLY_STRAIGHT 1e-5
#define OC_STRAIGHT 1e-15
#define OC_MAX_DEPTH 24                  /* halvings of one piece before its end is reached by a straight line */
#define OC_FIT_DEPTH 12                  /* halvings of a span approximated by cubics */
/* a span's derivative is read this far inside its end: at the end itself a family may answer for the next span, across a corner */
#if QAWS_SCALAR_IS_FLOAT
#define OC_SPAN_END 1e-4
#else
#define OC_SPAN_END 1e-9
#endif

typedef struct oc_p
{
	double x, y;
} oc_p;

typedef struct oc_cubic
{
	oc_p p0, p1, p2, p3;
} oc_cubic;

/* The output: a chain of cubics in world coordinates, written from pieces scaled by the current piece's scale and
   translation. */
typedef struct oc_out
{
	double* pts;                       /* x, y; 3 segments + 1 points */
	unsigned int segments, cap;
	int have, failed;
	double scale, tx, ty;
	oc_p prev, prev_world;             /* the last point, in the piece's coordinates and in the world's */
	oc_p cusp;
	int needs_cusp, cusp_clockwise;
} oc_out;

static oc_p oc_xy(double x, double y)
{
	oc_p p;
	p.x = x;
	p.y = y;
	return p;
}

static oc_p oc_add(oc_p a, oc_p b) { return oc_xy(a.x + b.x, a.y + b.y); }
static oc_p oc_sub(oc_p a, oc_p b) { return oc_xy(a.x - b.x, a.y - b.y); }
static oc_p oc_mul(oc_p a, double s) { return oc_xy(a.x * s, a.y * s); }
static double oc_dot(oc_p a, oc_p b) { return a.x * b.x + a.y * b.y; }
static double oc_cross(oc_p a, oc_p b) { return a.x * b.y - a.y * b.x; }
static double oc_len2(oc_p a) { return oc_dot(a, a); }
static double oc_dist(oc_p a, oc_p b) { return sqrt(oc_len2(oc_sub(b, a))); }
static double oc_dist2(oc_p a, oc_p b) { return oc_len2(oc_sub(b, a)); }
static int oc_eq(oc_p a, oc_p b) { return a.x == b.x && a.y == b.y; }
static int oc_near(oc_p a, oc_p b, double e) { return fabs(a.x - b.x) < e && fabs(a.y - b.y) < e; }
static oc_p oc_lerp(oc_p a, oc_p b, double t) { return oc_xy(a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t); }

static oc_p oc_unit(oc_p a)
{
	double m = oc_len2(a);
	if (m != 0.0 && m != 1.0)
	{
		m = sqrt(m);
		return oc_xy(a.x / m, a.y / m);
	}
	return a;
}

static oc_p oc_normal_of(oc_p from, oc_p to) { return oc_xy(to.y - from.y, -(to.x - from.x)); }
static int oc_clockwise(oc_p a, oc_p b, oc_p c) { double t = oc_cross(oc_sub(b, a), oc_sub(c, a)); return fabs(t) >= DBL_EPSILON && t > 0.0; }

static oc_p oc_at(oc_cubic const* c, double t)
{
	oc_p a = oc_lerp(c->p0, c->p1, t), b = oc_lerp(c->p1, c->p2, t), d = oc_lerp(c->p2, c->p3, t);
	return oc_lerp(oc_lerp(a, b, t), oc_lerp(b, d, t), t);
}

static oc_p oc_derived(oc_cubic const* c, double t)
{
	double it = 1.0 - t, d = t * t, a = -it * it, b = 1.0 - 4.0 * t + 3.0 * d, e = 2.0 * t - 3.0 * d;
	return oc_xy(3.0 * (a * c->p0.x + b * c->p1.x + e * c->p2.x + d * c->p3.x), 3.0 * (a * c->p0.y + b * c->p1.y + e * c->p2.y + d * c->p3.y));
}

static oc_p oc_normal(oc_cubic const* c, double t)
{
	oc_p d;
	if (fabs(t) < DBL_EPSILON)
	{
		if (oc_eq(c->p0, c->p1))
			return oc_eq(c->p0, c->p2) ? oc_normal_of(c->p0, c->p3) : oc_normal_of(c->p0, c->p2);
	}
	else if (fabs(t - 1.0) < DBL_EPSILON && oc_eq(c->p2, c->p3))
		return oc_eq(c->p1, c->p3) ? oc_normal_of(c->p0, c->p3) : oc_normal_of(c->p1, c->p3);
	d = oc_derived(c, t);
	return oc_xy(d.y, -d.x);
}

static oc_cubic oc_sub_curve(oc_cubic const* c, double t0, double t1)
{
	oc_cubic r;
	oc_p ab, bc, cd, abc, bcd, abcd, ab1, bc1, cd1, abc1, bcd1;
	double m;
	if (fabs(t0 - t1) < DBL_EPSILON)
	{
		r.p0 = r.p1 = r.p2 = r.p3 = oc_at(c, t0);
		return r;
	}
	if (t0 <= DBL_EPSILON && t1 >= 1.0 - DBL_EPSILON)
		return *c;
	ab = oc_lerp(c->p0, c->p1, t1);
	bc = oc_lerp(c->p1, c->p2, t1);
	cd = oc_lerp(c->p2, c->p3, t1);
	abc = oc_lerp(ab, bc, t1);
	bcd = oc_lerp(bc, cd, t1);
	abcd = oc_lerp(abc, bcd, t1);
	m = t0 / t1;
	ab1 = oc_lerp(c->p0, ab, m);
	bc1 = oc_lerp(ab, abc, m);
	cd1 = oc_lerp(abc, abcd, m);
	abc1 = oc_lerp(ab1, bc1, m);
	bcd1 = oc_lerp(bc1, cd1, m);
	r.p0 = oc_lerp(abc1, bcd1, m);
	r.p1 = bcd1;
	r.p2 = cd1;
	r.p3 = abcd;
	return r;
}

static int oc_accept(double* t, double root)
{
	if (root < -DBL_EPSILON || root > 1.0 + DBL_EPSILON)
		return 0;
	t[0] = root < 0.0 ? 0.0 : root > 1.0 ? 1.0 : root;
	return 1;
}

static int oc_quadratic(double a, double b, double c, double* roots)
{
	double delta = b * b - 4.0 * a * c;
	if (delta < 0.0)
		return 0;
	if (delta > 0.0)
	{
		double d = sqrt(delta), q = -0.5 * (b + (b < 0.0 ? -d : d)), r0 = q / a, r1 = c / q;
		int n;
		if (fabs(r0 - r1) < DBL_EPSILON)
			return oc_accept(roots, r0);
		n = oc_accept(roots, r0 < r1 ? r0 : r1);
		return n + oc_accept(roots + n, r0 < r1 ? r1 : r0);
	}
	return a != 0.0 ? oc_accept(roots, -0.5 * b / a) : 0;
}

/* Numerical Recipes 5.6: the real roots in [0, 1], ascending */
static int oc_cubic_roots(double c0, double c1, double c2, double c3, double* roots)
{
	double a, b, c, Q, R, disc, third, A;
	if (fabs(c0) < DBL_EPSILON)
		return oc_quadratic(c1, c2, c3, roots);
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
		double found[3];
		int n = 0, i, j, out = 0;
		n += oc_accept(found + n, m * cos(theta / 3.0) - third);
		n += oc_accept(found + n, m * cos((theta + 2.0 * OC_PI) / 3.0) - third);
		n += oc_accept(found + n, m * cos((theta - 2.0 * OC_PI) / 3.0) - third);
		for (i = 1; i < n; i++)
			for (j = i; j > 0 && found[j - 1] > found[j]; j--)
			{
				double s = found[j];
				found[j] = found[j - 1];
				found[j - 1] = s;
			}
		for (i = 0; i < n; i++)
			if (!out || fabs(found[i] - roots[out - 1]) >= DBL_EPSILON)
				roots[out++] = found[i];
		return out;
	}
	A = pow(fabs(R) + sqrt(disc), 1.0 / 3.0);
	if (R > 0.0)
		A = -A;
	if (A != 0.0)
		A += Q / A;
	return oc_accept(roots, A - third);
}

static int oc_max_curvature(oc_cubic const* k, double* t)
{
	oc_p a = oc_sub(k->p1, k->p0);
	oc_p b = oc_xy(k->p2.x - 2.0 * k->p1.x + k->p0.x, k->p2.y - 2.0 * k->p1.y + k->p0.y);
	oc_p c = oc_xy(k->p3.x + 3.0 * (k->p1.x - k->p2.x) - k->p0.x, k->p3.y + 3.0 * (k->p1.y - k->p2.y) - k->p0.y);
	return oc_cubic_roots(oc_dot(c, c), 3.0 * oc_dot(b, c), 2.0 * oc_dot(b, b) + oc_dot(c, a), oc_dot(a, b), t);
}

static int oc_inflections(oc_cubic const* k, double* t)
{
	oc_p a = oc_sub(k->p1, k->p0);
	oc_p b = oc_xy(k->p2.x - 2.0 * k->p1.x + k->p0.x, k->p2.y - 2.0 * k->p1.y + k->p0.y);
	oc_p c = oc_xy(k->p3.x + 3.0 * (k->p1.x - k->p2.x) - k->p0.x, k->p3.y + 3.0 * (k->p1.y - k->p2.y) - k->p0.y);
	return oc_quadratic(oc_cross(b, c), oc_cross(a, c), oc_cross(a, b), t);
}

static int oc_ray_hits(oc_cubic const* k, oc_p from, oc_p to, double* t)
{
	oc_p v = oc_sub(to, from);
	double ax = (k->p0.y - from.y) * v.x - (k->p0.x - from.x) * v.y;
	double bx = (k->p1.y - from.y) * v.x - (k->p1.x - from.x) * v.y;
	double cx = (k->p2.y - from.y) * v.x - (k->p2.x - from.x) * v.y;
	double dx = (k->p3.y - from.y) * v.x - (k->p3.x - from.x) * v.y;
	double a = dx, b = cx * 3, c = bx * 3, D = ax;
	return oc_cubic_roots(a - (D - c + b), b + (3 * D - 2 * c), c - 3 * D, D, t);
}

static void oc_append(oc_out* b, oc_p c1, oc_p c2, oc_p to)
{
	if (b->failed)
		return;
	if (b->segments >= b->cap)
	{
		unsigned int cap = b->cap ? 2 * b->cap : 32;
		double* grown = (double*)qaws_internal_realloc(b->pts, qaws_internal_mul3(sizeof(double), 2, 3 * (size_t)cap + 1));
		if (!grown)
		{
			b->failed = 1;
			return;
		}
		b->pts = grown;
		b->cap = cap;
	}
	{
		double* q = b->pts + 2 * (1 + 3 * (size_t)b->segments);
		q[0] = c1.x;
		q[1] = c1.y;
		q[2] = c2.x;
		q[3] = c2.y;
		q[4] = to.x;
		q[5] = to.y;
	}
	b->segments++;
}

static oc_p oc_world(oc_out const* b, oc_p p) { return oc_xy(p.x * b->scale + b->tx, p.y * b->scale + b->ty); }

/* The piece starts at to. A gap from the last piece's end (a corner between spans) is bridged by a line. */
static void oc_move(oc_out* b, oc_p to, double gap)
{
	oc_p w = oc_world(b, to);
	if (!b->have)
	{
		if (!b->pts)
		{
			b->pts = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 2, 3 * (size_t)32 + 1));
			b->cap = b->pts ? 32 : 0;
			b->failed = !b->pts;
		}
		if (b->pts)
		{
			b->pts[0] = w.x;
			b->pts[1] = w.y;
		}
		b->prev_world = w;
		b->have = 1;
	}
	else if (oc_dist(b->prev_world, w) > gap)
	{
		oc_append(b, oc_lerp(b->prev_world, w, 1.0 / 3.0), oc_lerp(b->prev_world, w, 2.0 / 3.0), w);
		b->prev_world = w;
	}
	b->prev = to;
}

static void oc_line(oc_out* b, oc_p to)
{
	if (!oc_eq(b->prev, to))
	{
		oc_p t = oc_world(b, to);
		oc_append(b, oc_lerp(b->prev_world, t, 1.0 / 3.0), oc_lerp(b->prev_world, t, 2.0 / 3.0), t);
		b->prev = to;
		b->prev_world = t;
	}
}

static void oc_curve(oc_out* b, oc_p c1, oc_p c2, oc_p to)
{
	if (!oc_eq(b->prev, c1) || !oc_eq(b->prev, c2) || !oc_eq(b->prev, to))
	{
		oc_p t = oc_world(b, to);
		oc_append(b, oc_world(b, c1), oc_world(b, c2), t);
		b->prev = to;
		b->prev_world = t;
	}
}

/* an arc on the unit circle (Riskus): the cubic through both ends, tangent there */
static void oc_unit_arc(oc_p p0, oc_p p3, oc_p* c1, oc_p* c2)
{
	double q1 = oc_len2(p0), q2 = q1 + oc_dot(p0, p3), k2 = (4.0 / 3.0) * (sqrt(2.0 * q1 * q2) - q2) / oc_cross(p0, p3);
	*c1 = oc_xy(p0.x - k2 * p0.y, p0.y + k2 * p0.x);
	*c2 = oc_xy(p3.x + k2 * p3.y, p3.y - k2 * p3.x);
}

static void oc_arc(oc_out* b, oc_p center, oc_p to, int clockwise)
{
	oc_p from = b->prev, a = oc_sub(from, center), e = oc_sub(to, center);
	double radius = sqrt(oc_len2(a)), start, sweep, step, s, c, cross, den;
	int steps, i;
	if (radius < OC_MIN_ARC_RADIUS)
		return;
	start = atan2(-a.y, a.x);
	if (start < 0)
		start += 2.0 * OC_PI;
	den = sqrt(oc_len2(a) * oc_len2(e));
	if (den < DBL_EPSILON)
		return;
	{
		double co = oc_dot(a, e) / den;
		if (co < -1.0 - FLT_EPSILON || co > 1.0 + FLT_EPSILON)
			return;
		sweep = acos(co < -1.0 ? -1.0 : co > 1.0 ? 1.0 : co);
	}
	if (fabs(sweep) < DBL_EPSILON)
		return;
	cross = oc_cross(a, e);
	if (fabs(cross) >= DBL_EPSILON && (cross > 0.0) != (clockwise != 0))
		sweep = 2.0 * OC_PI - sweep;
	steps = (int)ceil(sweep / (OC_PI / 2.0));
	step = sweep / steps * (clockwise ? -1.0 : 1.0);
	s = -sin(start);
	c = cos(start);
	for (i = 1; i <= steps; i++)
	{
		double a1 = start + step * i, s1 = -sin(a1), c1 = cos(a1);
		oc_p p1, p2;
		oc_unit_arc(oc_xy(c, s), oc_xy(c1, s1), &p1, &p2);
		p1 = oc_add(oc_mul(p1, radius), center);
		p2 = oc_add(oc_mul(p2, radius), center);
		oc_curve(b, p1, p2, i < steps ? oc_add(oc_mul(oc_xy(c1, s1), radius), center) : to);
		s = s1;
		c = c1;
	}
}

static void oc_cusp_arc(oc_out* b, oc_p to)
{
	if (b->needs_cusp)
	{
		b->needs_cusp = 0;
		oc_arc(b, b->cusp, to, b->cusp_clockwise);
		b->cusp = oc_xy(0, 0);
		b->cusp_clockwise = 0;
	}
}

/* whether the candidate lies at the offset from the original, measured along the original's normals */
static int oc_accept_offset(oc_cubic const* o, oc_cubic const* c, double offset, double max_error)
{
	static double const probes[3] = { 0.25, 0.5, 0.75 };
	double hit[3];
	int i;
	/* (with its end tangents almost parallel the fit can come out completely off) */
	if (oc_clockwise(o->p0, o->p1, o->p3) != oc_clockwise(c->p0, c->p1, c->p3)
		|| oc_clockwise(o->p0, o->p2, o->p3) != oc_clockwise(c->p0, c->p2, c->p3))
		return 0;
	for (i = 0; i < 3; i++)
	{
		oc_p op = oc_at(o, probes[i]), n = oc_normal(o, probes[i]);
		if (oc_ray_hits(c, op, oc_add(op, n), hit) != 1 || fabs(oc_dist(op, oc_at(c, hit[0])) - fabs(offset)) > max_error)
			return 0;
	}
	return 1;
}

typedef struct oc_tan
{
	oc_p s0, s1, e0, e1;               /* start tangent, end tangent (from the end backward) */
	double turn1, turn2;
	oc_p sn, en;
} oc_tan;

static double oc_turn(oc_p a, oc_p b, oc_p c) { return oc_cross(oc_unit(oc_sub(b, a)), oc_unit(oc_sub(c, a))); }

static oc_tan oc_tangents(oc_cubic const* c)
{
	oc_tan d;
	d.s0 = c->p0;
	d.s1 = !oc_near(c->p0, c->p1, DBL_EPSILON) ? c->p1 : !oc_near(c->p0, c->p2, DBL_EPSILON) ? c->p2 : c->p3;
	d.e0 = c->p3;
	d.e1 = !oc_near(c->p3, c->p2, DBL_EPSILON) ? c->p2 : !oc_near(c->p3, c->p1, DBL_EPSILON) ? c->p1 : c->p0;
	d.turn1 = oc_turn(d.s0, d.s1, d.e0);
	d.turn2 = oc_turn(d.s0, d.e1, d.e0);
	d.sn = oc_unit(oc_normal_of(d.s0, d.s1));
	d.en = oc_unit(oc_normal_of(d.e0, d.e1));
	return d;
}

static void oc_arc_offset(oc_out* b, double offset, oc_p center, oc_p from, oc_p to, int clockwise)
{
	double push = clockwise ? offset : -offset;
	oc_p a = oc_sub(from, center), e = oc_sub(to, center);
	oc_p p0 = oc_add(from, oc_mul(oc_unit(a), push)), p1 = oc_add(to, oc_mul(oc_unit(e), push));
	oc_cusp_arc(b, p0);
	/* (from the last point written: for a small arc the rounding of the start just computed can flip its turn) */
	oc_arc(b, center, p1, oc_clockwise(center, b->prev, p1));
}

static int oc_try_simple(oc_cubic const* c, oc_tan const* d, oc_out* b, double offset, double max_error)
{
	oc_p d1, d2, p0, p3, p, cc, p1, p2;
	oc_cubic k;
	double div, ca, cb;
	if (!((d->turn1 >= 0 && d->turn2 >= 0) || (d->turn1 <= 0 && d->turn2 <= 0)))
		return 0;
	d1 = oc_sub(c->p1, c->p0);
	d2 = oc_sub(c->p2, c->p3);
	div = oc_cross(d1, d2);
	if (fabs(div) < DBL_EPSILON)
		return 0;
	p0 = oc_add(d->s0, oc_mul(oc_unit(oc_normal_of(d->s0, d->s1)), offset));
	p3 = oc_sub(d->e0, oc_mul(oc_unit(oc_normal_of(d->e0, d->e1)), offset));
	p = oc_add(oc_at(c, 0.5), oc_mul(oc_unit(oc_normal(c, 0.5)), offset));
	cc = oc_mul(oc_sub(p, oc_mul(oc_add(p0, p3), 0.5)), 8.0 / 3.0);
	ca = oc_cross(cc, d2) / div;
	cb = oc_cross(d1, cc) / div;
	p1 = oc_add(p0, oc_mul(d1, ca));
	p2 = oc_add(p3, oc_mul(d2, cb));
	k.p0 = p0;
	k.p1 = p1;
	k.p2 = p2;
	k.p3 = p3;
	if (oc_dist2(k.p0, k.p1) + oc_dist2(k.p1, k.p2) + oc_dist2(k.p2, k.p3) <= OC_TINY_PERIMETER2)
		return 1;
	if (!oc_accept_offset(c, &k, offset, max_error))
		return 0;
	oc_cusp_arc(b, k.p0);
	oc_curve(b, k.p1, k.p2, k.p3);
	return 1;
}

static int oc_circle_hit(oc_p a, oc_p e, oc_p center, double radius)
{
	oc_p d = oc_sub(e, a), g = oc_sub(a, center);
	double A = oc_dot(d, d), B = 2.0 * oc_dot(g, d), C = oc_dot(g, g) - radius * radius, disc = B * B - 4.0 * A * C;
	if (disc > 0)
	{
		double s = sqrt(disc), t1 = (-B - s) / (2.0 * A), t2 = (-B + s) / (2.0 * A);
		return (t1 >= 0.0 && t1 <= 1.0) || (t2 >= 0.0 && t2 <= 1.0);
	}
	return 0;
}

static int oc_good_arc(oc_p center, double radius, oc_cubic const* c, double max_error, double t0, double t1)
{
	static double const probes[4] = { 0.2, 0.4, 0.6, 0.8 };
	double e, me;
	int i;
	if (radius > OC_MAX_ARC_RADIUS)
		return 0;
	e = max_error < radius / 3.0 ? max_error : radius / 3.0;
	me = e * (0.5 + 1e-4);
	for (i = 0; i < 4; i++)
	{
		double t = t0 + (t1 - t0) * probes[i];
		oc_p p = oc_at(c, t), n = oc_unit(oc_normal(c, t));
		if (!oc_circle_hit(oc_add(p, oc_mul(n, me)), oc_sub(p, oc_mul(n, me)), center, radius))
			return 0;
	}
	return 1;
}

/* the point on both lines (a0 a1) and (b0 b1); 0 when parallel */
static int oc_lines(oc_p a0, oc_p a1, oc_p b0, oc_p b1, oc_p* out)
{
	oc_p a = oc_sub(a1, a0), b = oc_sub(b0, b1), c = oc_sub(a0, b0);
	double den = a.y * b.x - a.x * b.y, na;
	if (den == 0)
		return 0;
	na = (b.y * c.x - b.x * c.y) / den;
	*out = oc_add(a0, oc_mul(a, na));
	return 1;
}

static int oc_try_arc(oc_cubic const* c, oc_tan const* d, oc_out* b, double offset, double max_error)
{
	oc_p from, to, V, G, m0, m1, C1, C2, s, e;
	double den, u, v, r0, r3, r03, hit[3], tG, sl, el;
	int clockwise;
	if (!((d->turn1 >= OC_NEARLY_STRAIGHT && d->turn2 >= OC_NEARLY_STRAIGHT) || (d->turn1 <= -OC_NEARLY_STRAIGHT && d->turn2 <= -OC_NEARLY_STRAIGHT)))
		return 0;
	from = oc_unit(oc_sub(d->s1, d->s0));
	to = oc_unit(oc_sub(d->e1, d->e0));
	den = to.x * from.y - to.y * from.x;
	if (fabs(den) < DBL_EPSILON)
		return 0;
	s = d->s0;
	e = d->e0;
	u = ((e.y - s.y) * to.x - (e.x - s.x) * to.y) / den;
	v = ((e.y - s.y) * from.x - (e.x - s.x) * from.y) / den;
	if (u < 0.0 || v < 0.0)
		return 0;
	V = oc_add(s, oc_mul(from, u));
	sl = oc_dist2(d->s0, d->s1);
	el = oc_dist2(d->e0, d->e1);
	if (oc_dist2(c->p0, V) < sl * 0.25 || oc_dist2(c->p3, V) < el * 0.25)
		return 0;
	r3 = oc_dist(c->p3, V);
	r0 = oc_dist(c->p0, V);
	r03 = oc_dist(c->p0, c->p3);
	/* the triangle's incentre: the biarc's split point */
	G = oc_mul(oc_add(oc_add(oc_mul(c->p0, r3), oc_mul(c->p3, r0)), oc_mul(V, r03)), 1.0 / (r3 + r0 + r03));
	m0 = oc_mul(oc_add(c->p0, G), 0.5);
	m1 = oc_mul(oc_add(G, c->p3), 0.5);
	if (!oc_lines(m0, oc_sub(m0, oc_normal_of(c->p0, G)), d->s0, oc_sub(d->s0, oc_normal_of(d->s0, d->s1)), &C1))
		return 0;
	if (oc_ray_hits(c, C1, G, hit) != 1)
		return 0;
	tG = hit[0];
	if (oc_dist(G, oc_at(c, tG)) > max_error)
		return 0;
	if (!oc_lines(m1, oc_sub(m1, oc_normal_of(G, c->p3)), d->e0, oc_add(d->e0, oc_normal_of(d->e0, d->e1)), &C2))
		return 0;
	clockwise = oc_clockwise(c->p0, V, c->p3);
	if (oc_near(C1, C2, OC_ARC_CENTER_EPS))
	{
		if (!oc_good_arc(C1, oc_dist(C1, c->p0), c, max_error, 0, 1))
			return 0;
		oc_arc_offset(b, offset, C1, c->p0, c->p3, clockwise);
		return 1;
	}
	if (!oc_good_arc(C1, oc_dist(C1, c->p0), c, max_error, 0, tG) || !oc_good_arc(C2, oc_dist(C2, c->p3), c, max_error, tG, 1))
		return 0;
	oc_arc_offset(b, offset, C1, c->p0, G, clockwise);
	oc_arc_offset(b, offset, C2, G, c->p3, clockwise);
	return 1;
}

static int oc_nearly_straight(oc_tan const* d)
{
	double minx = d->s0.x < d->e0.x ? d->s0.x : d->e0.x, maxx = d->s0.x > d->e0.x ? d->s0.x : d->e0.x;
	double miny = d->s0.y < d->e0.y ? d->s0.y : d->e0.y, maxy = d->s0.y > d->e0.y ? d->s0.y : d->e0.y;
	return minx <= d->s1.x && miny <= d->s1.y && maxx >= d->s1.x && maxy >= d->s1.y
		&& minx <= d->e1.x && miny <= d->e1.y && maxx >= d->e1.x && maxy >= d->e1.y
		&& fabs(d->turn1) < OC_NEARLY_STRAIGHT && fabs(d->turn2) < OC_NEARLY_STRAIGHT;
}

static int oc_straight(oc_tan const* d) { return fabs(d->turn1) < OC_STRAIGHT && fabs(d->turn2) < OC_STRAIGHT; }

static void oc_approximate(oc_cubic const* c, oc_tan const* d, oc_out* b, double offset, double max_error, int depth)
{
	if (oc_near(c->p0, c->p3, OC_POINT_EPS) && oc_near(c->p0, c->p1, OC_POINT_EPS) && oc_near(c->p0, c->p2, OC_POINT_EPS))
		return;
	if (oc_nearly_straight(d))
	{
		if (oc_straight(d))
		{
			oc_p n = oc_unit(oc_normal_of(c->p0, c->p1));
			oc_cusp_arc(b, oc_add(c->p0, oc_mul(n, offset)));
			oc_line(b, oc_add(c->p1, oc_mul(n, offset)));
		}
		else
		{
			oc_cusp_arc(b, oc_add(d->s0, oc_mul(d->sn, offset)));
			oc_curve(b, oc_add(d->s1, oc_mul(d->sn, offset)), oc_sub(d->e1, oc_mul(d->en, offset)), oc_sub(d->e0, oc_mul(d->en, offset)));
		}
		return;
	}
	if (oc_try_simple(c, d, b, offset, max_error) || oc_try_arc(c, d, b, offset, max_error))
		return;
	if (depth >= OC_MAX_DEPTH)
	{
		/* (the piece is microscopic by now: its end reached directly keeps the chain connected) */
		oc_cusp_arc(b, oc_add(d->s0, oc_mul(d->sn, offset)));
		oc_line(b, oc_sub(d->e0, oc_mul(d->en, offset)));
		return;
	}
	{
		oc_cubic l, r;
		oc_tan dl, dr;
		oc_p m = oc_mul(oc_add(c->p1, c->p2), 0.5);
		l.p0 = c->p0;
		l.p1 = oc_mul(oc_add(c->p0, c->p1), 0.5);
		r.p3 = c->p3;
		r.p2 = oc_mul(oc_add(c->p2, c->p3), 0.5);
		l.p2 = oc_mul(oc_add(l.p1, m), 0.5);
		r.p1 = oc_mul(oc_add(r.p2, m), 0.5);
		l.p3 = r.p0 = oc_mul(oc_add(l.p2, r.p1), 0.5);
		dl = oc_tangents(&l);
		oc_approximate(&l, &dl, b, offset, max_error, depth + 1);
		dr = oc_tangents(&r);
		oc_approximate(&r, &dr, b, offset, max_error, depth + 1);
	}
}

static double oc_search_toward(oc_cubic const* c, double t, double target)
{
	int i;
	for (i = 0; i < OC_NEAR_CUSP_ITERATIONS; i++)
	{
		if (oc_len2(oc_derived(c, t)) < OC_CUSP_DERIVATIVE2 * 2.0)
			return t;
		t = (t + target) / 2.0;
	}
	return t;
}

/* every point on one line: the runs between cusps are offset as lines, a half circle at each cusp */
static void oc_linear(oc_cubic const* c, oc_tan const* d, oc_out* b, double offset, double const* mc, int nmc)
{
	oc_p n = oc_unit(oc_normal_of(d->s0, d->s1)), prev = d->s0, prev_off = oc_add(prev, oc_mul(n, offset));
	int i;
	oc_move(b, prev_off, 0.0);
	for (i = 0; i < nmc; i++)
		if (oc_len2(oc_derived(c, mc[i])) <= 1e-9)
		{
			oc_p at = oc_at(c, mc[i]), nl = oc_unit(oc_normal_of(prev, at)), arc_to = oc_sub(at, oc_mul(nl, offset));
			oc_line(b, oc_add(at, oc_mul(nl, offset)));
			oc_arc(b, at, arc_to, oc_clockwise(prev, prev_off, at));
			prev = at;
			prev_off = arc_to;
		}
	oc_line(b, oc_sub(d->e0, oc_mul(oc_unit(oc_normal_of(d->e0, d->e1)), offset)));
}

static int oc_merge(double* t, int n, double const* s, int count, double eps)
{
	int i, j;
	for (i = 0; i < count; i++)
	{
		int seen = fabs(s[i]) < eps || fabs(s[i] - 1.0) < eps;
		for (j = 0; j < n && !seen; j++)
			seen = fabs(s[i] - t[j]) < eps;
		if (!seen)
			t[n++] = s[i];
	}
	return n;
}

/* split at inflections and curvature maxima, cusps turned into arcs about them */
static void oc_general(oc_cubic const* c, oc_tan const* d, oc_out* b, double offset, double max_error)
{
	double mc[3], inf[2], t[5], prev = 0;
	int nmc = oc_max_curvature(c, mc), ninf, count, i, j;
	oc_cubic k;
	oc_tan nd;
	if (oc_straight(d))
	{
		oc_linear(c, d, b, offset, mc, nmc);
		return;
	}
	ninf = oc_inflections(c, inf);
	count = oc_merge(t, 0, inf, ninf, 1e-7);
	count = oc_merge(t, count, mc, nmc, 1e-5);
	for (i = 1; i < count; i++)
		for (j = i; j > 0 && t[j - 1] > t[j]; j--)
		{
			double s = t[j];
			t[j] = t[j - 1];
			t[j - 1] = s;
		}
	for (i = 0; i < count; i++)
	{
		double T = t[i];
		if (oc_len2(oc_derived(c, T)) < OC_CUSP_DERIVATIVE2)
		{
			double lo = prev + (T - prev) * 0.8, hi = T + ((i == count - 1 ? 1.0 : t[i + 1]) - T) * 0.2, t1, t2;
			t1 = oc_search_toward(c, lo > T - 0.05 ? lo : T - 0.05, T);
			t2 = oc_search_toward(c, hi < T + 0.05 ? hi : T + 0.05, T);
			k = oc_sub_curve(c, prev, t1);
			nd = oc_tangents(&k);
			oc_approximate(&k, &nd, b, offset, max_error, 0);
			b->cusp = oc_at(c, T);
			b->needs_cusp = 1;
			b->cusp_clockwise = oc_clockwise(k.p3, b->cusp, oc_at(c, t2));
			prev = t2;
		}
		else
		{
			k = oc_sub_curve(c, prev, T);
			nd = oc_tangents(&k);
			oc_approximate(&k, &nd, b, offset, max_error, 0);
			prev = T;
		}
	}
	k = oc_sub_curve(c, prev, 1.0);
	nd = oc_tangents(&k);
	oc_approximate(&k, &nd, b, offset, max_error, 0);
}

static void oc_cubic_offset(oc_cubic const* in, oc_out* b, double distance, double tolerance, double gap)
{
	double minx = in->p0.x, maxx = minx, miny = in->p0.y, maxy = miny, m, so, max_error;
	oc_p const* q = &in->p0;
	oc_cubic c;
	oc_tan d;
	int i;
	for (i = 1; i < 4; i++)
	{
		minx = q[i].x < minx ? q[i].x : minx;
		maxx = q[i].x > maxx ? q[i].x : maxx;
		miny = q[i].y < miny ? q[i].y : miny;
		maxy = q[i].y > maxy ? q[i].y : maxy;
	}
	if (maxx - minx < OC_POINT_EPS && maxy - miny < OC_POINT_EPS)
		return;
	m = (maxx - minx > maxy - miny ? maxx - minx : maxy - miny) / 2.0;
	b->scale = m;
	b->tx = (minx + maxx) / 2.0;
	b->ty = (miny + maxy) / 2.0;
	so = distance / m;
	max_error = tolerance / m;
	c.p0 = oc_xy((in->p0.x - b->tx) / m, (in->p0.y - b->ty) / m);
	c.p1 = oc_xy((in->p1.x - b->tx) / m, (in->p1.y - b->ty) / m);
	c.p2 = oc_xy((in->p2.x - b->tx) / m, (in->p2.y - b->ty) / m);
	c.p3 = oc_xy((in->p3.x - b->tx) / m, (in->p3.y - b->ty) / m);
	if (oc_dist2(c.p0, c.p1) < 1e-12)
		c.p1 = c.p0;
	if (oc_dist2(c.p3, c.p2) < 1e-12)
		c.p2 = c.p3;
	d = oc_tangents(&c);
	if (oc_nearly_straight(&d))
	{
		if (oc_straight(&d))
		{
			oc_p n = oc_unit(oc_normal_of(c.p0, c.p3));
			oc_move(b, oc_add(c.p0, oc_mul(n, so)), gap);
			oc_line(b, oc_add(c.p3, oc_mul(n, so)));
		}
		else
		{
			oc_move(b, oc_add(d.s0, oc_mul(d.sn, so)), gap);
			oc_curve(b, oc_add(d.s1, oc_mul(d.sn, so)), oc_sub(d.e1, oc_mul(d.en, so)), oc_sub(d.e0, oc_mul(d.en, so)));
		}
		return;
	}
	oc_move(b, oc_add(d.s0, oc_mul(d.sn, so)), gap);
	/* (the arc first: a curve meant as a circle needs none of the splitting below) */
	if (!oc_try_arc(&c, &d, b, so, max_error))
		oc_general(&c, &d, b, so, max_error);
}

/* a span of the input as cubics within fit (one exactly when the span is a cubic polynomial), each offset */
static qaws_status oc_span(qaws_curve const* curve, unsigned int span, double u0, double u1, int depth, oc_out* b, double distance, double tolerance, double fit, double gap)
{
	qaws_eval_result_2d r0, r1;
	unsigned int const flags = QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1;
	double h = (double)(curve->span_boundaries[span + 1] - curve->span_boundaries[span]) * (u1 - u0);
	qaws_status st = qaws_curve_evaluate_span_2d(&r0, curve, span, (qaws_scalar)u0, flags);
	oc_cubic c;
	int i, fits = 1;
	memset(&r1, 0, sizeof r1);
	if (st == QAWS_STATUS_OK)
		st = qaws_curve_evaluate_span_2d(&r1, curve, span, (qaws_scalar)(u1 < 1.0 ? u1 : 1.0 - OC_SPAN_END), flags);
	if (st == QAWS_STATUS_OK && u1 >= 1.0)
	{
		qaws_eval_result_2d e;
		st = qaws_curve_evaluate_span_2d(&e, curve, span, (qaws_scalar)1, QAWS_EVAL_FLAG_POSITION);
		r1.position = e.position;
	}
	if (st != QAWS_STATUS_OK)
		return st;
	c.p0 = oc_xy((double)r0.position.x, (double)r0.position.y);
	c.p3 = oc_xy((double)r1.position.x, (double)r1.position.y);
	c.p1 = oc_add(c.p0, oc_mul(oc_xy((double)r0.d1.x, (double)r0.d1.y), h / 3.0));
	c.p2 = oc_sub(c.p3, oc_mul(oc_xy((double)r1.d1.x, (double)r1.d1.y), h / 3.0));
	for (i = 1; i < 4 && fits; i++)
	{
		qaws_eval_result_2d m;
		double s = (double)i / 4.0;
		st = qaws_curve_evaluate_span_2d(&m, curve, span, (qaws_scalar)(u0 + (u1 - u0) * s), QAWS_EVAL_FLAG_POSITION);
		if (st != QAWS_STATUS_OK)
			return st;
		fits = oc_dist(oc_at(&c, s), oc_xy((double)m.position.x, (double)m.position.y)) <= fit;
	}
	if (!fits && depth < OC_FIT_DEPTH)
	{
		double mid = (u0 + u1) / 2.0;
		st = oc_span(curve, span, u0, mid, depth + 1, b, distance, tolerance, fit, gap);
		return st == QAWS_STATUS_OK ? oc_span(curve, span, mid, u1, depth + 1, b, distance, tolerance, fit, gap) : st;
	}
	oc_cubic_offset(&c, b, distance, tolerance, gap);
	return b->failed ? QAWS_STATUS_ALLOCATION_FAILURE : QAWS_STATUS_OK;
}

qaws_status qaws_curve_offset_cubic_2d(
	qaws_curve** out_curve,
	qaws_curve const* curve,
	qaws_scalar distance,
	qaws_scalar tolerance)
{
	oc_out b;
	qaws_status st = QAWS_STATUS_OK;
	unsigned int i;
	double gap;
	if (!out_curve || !curve)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_curve = NULL;
	if (curve->dimension != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;
	if (!(fabs((double)distance) <= DBL_MAX) || !((double)tolerance > 0.0) || !((double)tolerance <= DBL_MAX))
		return QAWS_STATUS_INVALID_ARGUMENT;
	memset(&b, 0, sizeof b);
	/* (joints closer than a small part of the tolerance are the same point) */
	gap = (double)tolerance * 1e-3;
	for (i = 0; i < curve->span_count && st == QAWS_STATUS_OK; i++)
		st = oc_span(curve, i, 0.0, 1.0, 0, &b, -(double)distance, (double)tolerance, (double)tolerance * 0.25, gap);
	if (st == QAWS_STATUS_OK && !b.segments)
		st = QAWS_STATUS_DEGENERATE_CURVE;
	if (st == QAWS_STATUS_OK)
		st = qaws_internal_bezier_chain(out_curve, b.pts, 3, b.segments, curve->allocator);
	qaws_internal_free(b.pts);
	return st;
}
