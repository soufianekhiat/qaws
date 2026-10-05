#include "qaws_yuksel.h"
#include "qaws_curve.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_curve.h"
#include "internal/qaws_internal_validation.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef QAWS_PI
#define QAWS_PI ((qaws_scalar)3.14159265358979323846)
#endif

#define QAWS_PI_2 ((qaws_scalar)(QAWS_PI * (qaws_scalar)0.5))

/* ------------------------------------------------------------------ */
/*  Vector helpers (2D/3D, dimension-generic via raw scalars)          */
/* ------------------------------------------------------------------ */

static qaws_scalar vdot(qaws_scalar const* a, qaws_scalar const* b, unsigned int dim)
{
	qaws_scalar r = (qaws_scalar)0;
	unsigned int i;
	for (i = 0; i < dim; i++) r += a[i] * b[i];
	return r;
}

static void vsub(qaws_scalar* out, qaws_scalar const* a, qaws_scalar const* b, unsigned int dim)
{
	unsigned int i;
	for (i = 0; i < dim; i++) out[i] = a[i] - b[i];
}

static void vadd(qaws_scalar* out, qaws_scalar const* a, qaws_scalar const* b, unsigned int dim)
{
	unsigned int i;
	for (i = 0; i < dim; i++) out[i] = a[i] + b[i];
}

static void vscale(qaws_scalar* out, qaws_scalar const* a, qaws_scalar s, unsigned int dim)
{
	unsigned int i;
	for (i = 0; i < dim; i++) out[i] = a[i] * s;
}

static void vcopy(qaws_scalar* out, qaws_scalar const* a, unsigned int dim)
{
	unsigned int i;
	for (i = 0; i < dim; i++) out[i] = a[i];
}

static qaws_scalar vlength(qaws_scalar const* a, unsigned int dim)
{
	return (qaws_scalar)sqrt((double)vdot(a, a, dim));
}

/* ------------------------------------------------------------------ */
/*  Osculating plane helper for 3D circular/elliptical fitting         */
/* ------------------------------------------------------------------ */

/* Build an orthonormal frame from triplet (pj, pi, pk) in 3D.
   Returns 1 on success, 0 if points are collinear/degenerate.
   On success: u, w are unit vectors spanning the osculating plane,
   and normal = u x w. Origin is at pi. */
static int compute_osculating_frame(
	qaws_scalar const* pj, qaws_scalar const* pi, qaws_scalar const* pk,
	qaws_scalar* u_out, qaws_scalar* w_out)
{
	qaws_scalar v1[3], v2[3], normal[3];
	qaws_scalar len_n, len_v1;

	vsub(v1, pj, pi, 3);
	vsub(v2, pk, pi, 3);

	/* normal = v1 x v2 */
	normal[0] = v1[1] * v2[2] - v1[2] * v2[1];
	normal[1] = v1[2] * v2[0] - v1[0] * v2[2];
	normal[2] = v1[0] * v2[1] - v1[1] * v2[0];

	len_n = vlength(normal, 3);
	if (len_n < (qaws_scalar)1e-10)
		return 0; /* collinear / degenerate */

	/* Normalize the normal */
	normal[0] /= len_n;
	normal[1] /= len_n;
	normal[2] /= len_n;

	/* u = normalize(v1) */
	len_v1 = vlength(v1, 3);
	if (len_v1 < (qaws_scalar)1e-10)
		return 0;
	u_out[0] = v1[0] / len_v1;
	u_out[1] = v1[1] / len_v1;
	u_out[2] = v1[2] / len_v1;

	/* w = normal x u (already unit length since normal and u are unit and perpendicular) */
	w_out[0] = normal[1] * u_out[2] - normal[2] * u_out[1];
	w_out[1] = normal[2] * u_out[0] - normal[0] * u_out[2];
	w_out[2] = normal[0] * u_out[1] - normal[1] * u_out[0];

	return 1;
}

/* Project a 3D point onto a 2D frame (origin + u, w axes).
   Returns the 2D coordinates (px, py). */
static void project_to_2d(
	qaws_scalar const* point, qaws_scalar const* origin,
	qaws_scalar const* u, qaws_scalar const* w,
	qaws_scalar* px, qaws_scalar* py)
{
	qaws_scalar rel[3];
	vsub(rel, point, origin, 3);
	*px = vdot(rel, u, 3);
	*py = vdot(rel, w, 3);
}

/* ------------------------------------------------------------------ */
/*  Cubic root via bisection (from Yuksel's reference code)           */
/* ------------------------------------------------------------------ */

static qaws_scalar cubic_root(qaws_scalar d, qaws_scalar c, qaws_scalar b, qaws_scalar a)
{
	qaws_scalar value = (d + (qaws_scalar)3 * c + (qaws_scalar)3 * b + a) / (qaws_scalar)8;
	if (value >= (qaws_scalar)0.000001)
		return cubic_root(d, (d + c) / (qaws_scalar)2,
			(d + (qaws_scalar)2 * c + b) / (qaws_scalar)4, value) / (qaws_scalar)2;
	if (value <= (qaws_scalar)-0.000001)
		return (qaws_scalar)0.5 + cubic_root(value,
			(c + (qaws_scalar)2 * b + a) / (qaws_scalar)4,
			(b + a) / (qaws_scalar)2, a) / (qaws_scalar)2;
	return (qaws_scalar)0.5;
}

/* ------------------------------------------------------------------ */
/*  Find optimal t for quadratic Bezier through 3 points              */
/* ------------------------------------------------------------------ */

static qaws_scalar max_curvature_t(qaws_scalar const* p0, qaws_scalar const* pin,
	qaws_scalar const* p2, unsigned int dim)
{
	qaws_scalar v0[3], v2[3];
	qaws_scalar c;
	vsub(v0, p0, pin, dim);
	vsub(v2, p2, pin, dim);
	c = vdot(v0, v2, dim);
	return cubic_root(
		-vdot(v0, v0, dim),
		-c / (qaws_scalar)3,
		c / (qaws_scalar)3,
		vdot(v2, v2, dim));
}

/* ------------------------------------------------------------------ */
/*  Compute sub-curves for each control point                         */
/* ------------------------------------------------------------------ */

static void compute_bezier_subcurve(
	qaws_scalar const* pts, unsigned int n, unsigned int i,
	int closed, unsigned int dim, qaws_yuksel_subcurve* sc)
{
	unsigned int j, k;
	qaws_scalar const* pj;
	qaws_scalar const* pi;
	qaws_scalar const* pk;
	qaws_scalar t1, omt, denom;
	qaws_scalar tmp1[3], tmp2[3], tmp3[3];

	if (closed) {
		j = (i + n - 1) % n;
		k = (i + 1) % n;
	} else {
		j = (i > 0) ? i - 1 : 0;
		k = (i < n - 1) ? i + 1 : n - 1;
	}

	pj = &pts[j * dim];
	pi = &pts[i * dim];
	pk = &pts[k * dim];

	vcopy(sc->p0, pj, dim);
	vcopy(sc->p2, pk, dim);

	/* Degenerate case: j == i or k == i */
	if (j == i || k == i) {
		vcopy(sc->p1, pi, dim);
		sc->t1 = (qaws_scalar)0.5;
		sc->use_arc = 0;
		return;
	}

	t1 = max_curvature_t(pj, pi, pk, dim);
	/* Clamp to avoid singularities */
	if (t1 < (qaws_scalar)0.01) t1 = (qaws_scalar)0.01;
	if (t1 > (qaws_scalar)0.99) t1 = (qaws_scalar)0.99;

	sc->t1 = t1;
	sc->use_arc = 0;

	/* P1 = (P[i] - (1-t)^2 * P[j] - t^2 * P[k]) / (2*(1-t)*t) */
	omt = (qaws_scalar)1 - t1;
	denom = (qaws_scalar)2 * omt * t1;

	vscale(tmp1, pj, omt * omt, dim);
	vscale(tmp2, pk, t1 * t1, dim);
	vadd(tmp3, tmp1, tmp2, dim);
	vsub(tmp1, pi, tmp3, dim);
	vscale(sc->p1, tmp1, (qaws_scalar)1 / denom, dim);
}

static void compute_circular_subcurve(
	qaws_scalar const* pts, unsigned int n, unsigned int i,
	int closed, unsigned int dim, qaws_yuksel_subcurve* sc)
{
	unsigned int j, k;
	qaws_scalar const* pj;
	qaws_scalar const* pi;
	qaws_scalar const* pk;
	qaws_scalar mid_jk[3], mid_ij[3], dir_jk[3], dir_ij[3];
	qaws_scalar perp_jk[3], perp_ij[3];
	qaws_scalar diff[3], t_param, len_jk, len_ij;
	qaws_scalar center_to_j[3], center_to_k[3];
	qaws_scalar angle_j, angle_k;
	unsigned int d;

	if (closed) {
		j = (i + n - 1) % n;
		k = (i + 1) % n;
	} else {
		j = (i > 0) ? i - 1 : 0;
		k = (i < n - 1) ? i + 1 : n - 1;
	}

	pj = &pts[j * dim];
	pi = &pts[i * dim];
	pk = &pts[k * dim];

	vcopy(sc->p0, pj, dim);
	vcopy(sc->p2, pk, dim);

	if (j == i || k == i) {
		vcopy(sc->p1, pi, dim);
		sc->t1 = (qaws_scalar)0.5;
		sc->use_arc = 0;
		return;
	}

	/* 3D case: project triplet onto osculating plane, fit in 2D, map back */
	if (dim == 3) {
		qaws_scalar frame_u[3], frame_w[3];
		qaws_scalar pj2x, pj2y, pi2x, pi2y, pk2x, pk2y;
		qaws_scalar mid_jk2[2], mid_ij2[2], dir_jk2[2], dir_ij2[2];
		qaws_scalar perp_jk2[2], perp_ij2[2];
		qaws_scalar diff2[2], center2x, center2y;
		qaws_scalar ctj2[2], ctk2[2];
		qaws_scalar denom_cross, len_jk2, len_ij2;

		if (!compute_osculating_frame(pj, pi, pk, frame_u, frame_w)) {
			compute_bezier_subcurve(pts, n, i, closed, dim, sc);
			return;
		}

		/* Project 3D points to 2D in the osculating plane (origin = pi) */
		project_to_2d(pj, pi, frame_u, frame_w, &pj2x, &pj2y);
		pi2x = (qaws_scalar)0;
		pi2y = (qaws_scalar)0;
		project_to_2d(pk, pi, frame_u, frame_w, &pk2x, &pk2y);

		/* Midpoints and directions in 2D */
		mid_jk2[0] = (pj2x + pk2x) / (qaws_scalar)2;
		mid_jk2[1] = (pj2y + pk2y) / (qaws_scalar)2;
		mid_ij2[0] = (pi2x + pj2x) / (qaws_scalar)2;
		mid_ij2[1] = (pi2y + pj2y) / (qaws_scalar)2;
		dir_jk2[0] = pk2x - pj2x;
		dir_jk2[1] = pk2y - pj2y;
		dir_ij2[0] = pj2x - pi2x;
		dir_ij2[1] = pj2y - pi2y;

		len_jk2 = (qaws_scalar)sqrt((double)(dir_jk2[0]*dir_jk2[0] + dir_jk2[1]*dir_jk2[1]));
		len_ij2 = (qaws_scalar)sqrt((double)(dir_ij2[0]*dir_ij2[0] + dir_ij2[1]*dir_ij2[1]));

		if (len_jk2 < (qaws_scalar)1e-10 || len_ij2 < (qaws_scalar)1e-10) {
			compute_bezier_subcurve(pts, n, i, closed, dim, sc);
			return;
		}

		/* Perpendicular bisectors (2D rotate 90) */
		perp_jk2[0] = -dir_jk2[1];
		perp_jk2[1] = dir_jk2[0];
		perp_ij2[0] = -dir_ij2[1];
		perp_ij2[1] = dir_ij2[0];

		/* Find intersection */
		denom_cross = perp_jk2[0] * perp_ij2[1] - perp_jk2[1] * perp_ij2[0];
		if (denom_cross < (qaws_scalar)1e-10 && denom_cross > (qaws_scalar)-1e-10) {
			compute_bezier_subcurve(pts, n, i, closed, dim, sc);
			return;
		}
		diff2[0] = mid_ij2[0] - mid_jk2[0];
		diff2[1] = mid_ij2[1] - mid_jk2[1];
		t_param = (diff2[0] * perp_ij2[1] - diff2[1] * perp_ij2[0]) / denom_cross;

		center2x = mid_jk2[0] + t_param * perp_jk2[0];
		center2y = mid_jk2[1] + t_param * perp_jk2[1];

		ctj2[0] = pj2x - center2x;
		ctj2[1] = pj2y - center2y;
		ctk2[0] = pk2x - center2x;
		ctk2[1] = pk2y - center2y;

		sc->radius = (qaws_scalar)sqrt((double)(ctj2[0]*ctj2[0] + ctj2[1]*ctj2[1]));
		if (sc->radius < (qaws_scalar)1e-10) {
			compute_bezier_subcurve(pts, n, i, closed, dim, sc);
			return;
		}

		/* Map center back to 3D: center_3d = pi + center2x * frame_u + center2y * frame_w */
		for (d = 0; d < 3; d++)
			sc->center[d] = pi[d] + center2x * frame_u[d] + center2y * frame_w[d];

		/* Axis frame in 3D: axis_u points from center toward pj */
		{
			qaws_scalar au_len;
			qaws_scalar au[3];
			vsub(au, pj, sc->center, 3);
			au_len = vlength(au, 3);
			if (au_len < (qaws_scalar)1e-10) {
				compute_bezier_subcurve(pts, n, i, closed, dim, sc);
				return;
			}
			for (d = 0; d < 3; d++)
				sc->axis_u[d] = au[d] / au_len;
		}

		/* axis_v = normal x axis_u (in-plane, perpendicular to axis_u) */
		{
			qaws_scalar normal[3];
			qaws_scalar v1[3], v2[3], len_n;
			vsub(v1, pj, pi, 3);
			vsub(v2, pk, pi, 3);
			normal[0] = v1[1] * v2[2] - v1[2] * v2[1];
			normal[1] = v1[2] * v2[0] - v1[0] * v2[2];
			normal[2] = v1[0] * v2[1] - v1[1] * v2[0];
			len_n = vlength(normal, 3);
			normal[0] /= len_n;
			normal[1] /= len_n;
			normal[2] /= len_n;
			sc->axis_v[0] = normal[1] * sc->axis_u[2] - normal[2] * sc->axis_u[1];
			sc->axis_v[1] = normal[2] * sc->axis_u[0] - normal[0] * sc->axis_u[2];
			sc->axis_v[2] = normal[0] * sc->axis_u[1] - normal[1] * sc->axis_u[0];
		}

		/* Angles: axis_u points to pj, so angle_start = 0 */
		{
			qaws_scalar ctk3[3];
			vsub(ctk3, pk, sc->center, 3);
			angle_k = (qaws_scalar)atan2(
				(double)vdot(ctk3, sc->axis_v, 3),
				(double)vdot(ctk3, sc->axis_u, 3));
		}

		sc->angle_start = (qaws_scalar)0;
		sc->angle_end = angle_k;
		sc->radius_b = sc->radius;
		sc->use_arc = 1;
		sc->t1 = (qaws_scalar)0.5;
		vcopy(sc->p1, pi, dim);
		return;
	}

	/* 2D case: original code */

	/* Midpoints of (j,k) and (i,j) */
	for (d = 0; d < dim; d++) {
		mid_jk[d] = (pj[d] + pk[d]) / (qaws_scalar)2;
		mid_ij[d] = (pi[d] + pj[d]) / (qaws_scalar)2;
		dir_jk[d] = pk[d] - pj[d];
		dir_ij[d] = pj[d] - pi[d];
	}

	len_jk = vlength(dir_jk, dim);
	len_ij = vlength(dir_ij, dim);

	if (len_jk < (qaws_scalar)1e-10 || len_ij < (qaws_scalar)1e-10) {
		compute_bezier_subcurve(pts, n, i, closed, dim, sc);
		return;
	}

	/* Perpendicular bisectors (2D: rotate 90 degrees) */
	perp_jk[0] = -dir_jk[1];
	perp_jk[1] = dir_jk[0];
	perp_ij[0] = -dir_ij[1];
	perp_ij[1] = dir_ij[0];

	/* Find intersection: mid_jk + t * perp_jk = mid_ij + s * perp_ij */
	{
		qaws_scalar denom_cross = perp_jk[0] * perp_ij[1] - perp_jk[1] * perp_ij[0];
		if ((denom_cross < (qaws_scalar)1e-10) && (denom_cross > (qaws_scalar)-1e-10)) {
			/* Collinear: fall back to bezier */
			compute_bezier_subcurve(pts, n, i, closed, dim, sc);
			return;
		}
		vsub(diff, mid_ij, mid_jk, dim);
		t_param = (diff[0] * perp_ij[1] - diff[1] * perp_ij[0]) / denom_cross;
	}

	sc->center[0] = mid_jk[0] + t_param * perp_jk[0];
	sc->center[1] = mid_jk[1] + t_param * perp_jk[1];

	vsub(center_to_j, pj, sc->center, dim);
	vsub(center_to_k, pk, sc->center, dim);

	sc->radius = vlength(center_to_j, dim);
	if (sc->radius < (qaws_scalar)1e-10) {
		compute_bezier_subcurve(pts, n, i, closed, dim, sc);
		return;
	}

	/* Axis frame for the arc */
	sc->axis_u[0] = center_to_j[0] / sc->radius;
	sc->axis_u[1] = center_to_j[1] / sc->radius;
	sc->axis_v[0] = -sc->axis_u[1];
	sc->axis_v[1] = sc->axis_u[0];

	/* Angles from center to p_j and p_k */
	angle_j = (qaws_scalar)0; /* by construction, axis_u points to p_j */
	angle_k = (qaws_scalar)atan2(
		(double)(center_to_k[0] * sc->axis_v[0] + center_to_k[1] * sc->axis_v[1]),
		(double)(center_to_k[0] * sc->axis_u[0] + center_to_k[1] * sc->axis_u[1]));

	sc->angle_start = angle_j;
	sc->angle_end = angle_k;
	sc->radius_b = sc->radius;
	sc->use_arc = 1;

	/* Also compute bezier fallback params for blending derivatives */
	sc->t1 = (qaws_scalar)0.5;
	vcopy(sc->p1, pi, dim); /* Not used directly for arc eval */
}

static void compute_elliptical_subcurve(
	qaws_scalar const* pts, unsigned int n, unsigned int i,
	int closed, unsigned int dim, qaws_yuksel_subcurve* sc)
{
	unsigned int j, k;
	qaws_scalar const* pj;
	qaws_scalar const* pi;
	qaws_scalar const* pk;
	qaws_scalar mid[3], dir[3];
	qaws_scalar len, half_len;
	qaws_scalar pi_rel[3], px, py;
	qaws_scalar a, b;
	int iter;
	unsigned int d;

	if (closed) {
		j = (i + n - 1) % n;
		k = (i + 1) % n;
	} else {
		j = (i > 0) ? i - 1 : 0;
		k = (i < n - 1) ? i + 1 : n - 1;
	}

	pj = &pts[j * dim];
	pi = &pts[i * dim];
	pk = &pts[k * dim];

	vcopy(sc->p0, pj, dim);
	vcopy(sc->p2, pk, dim);

	if (j == i || k == i) {
		vcopy(sc->p1, pi, dim);
		sc->t1 = (qaws_scalar)0.5;
		sc->use_arc = 0;
		return;
	}

	/* 3D case: project triplet onto osculating plane, fit ellipse in 2D, map back */
	if (dim == 3) {
		qaws_scalar frame_u[3], frame_w[3];
		qaws_scalar pj2x, pj2y, pi2x, pi2y, pk2x, pk2y;
		qaws_scalar mid2[2], dir2[2], len2, half_len2;
		qaws_scalar au2[2], av2[2];
		qaws_scalar pirel2[2], px2, py2;
		qaws_scalar a2, b2;
		int iter2;
		qaws_scalar center2x, center2y;

		if (!compute_osculating_frame(pj, pi, pk, frame_u, frame_w)) {
			compute_bezier_subcurve(pts, n, i, closed, dim, sc);
			return;
		}

		/* Project 3D points to 2D (origin = pi) */
		project_to_2d(pj, pi, frame_u, frame_w, &pj2x, &pj2y);
		pi2x = (qaws_scalar)0;
		pi2y = (qaws_scalar)0;
		project_to_2d(pk, pi, frame_u, frame_w, &pk2x, &pk2y);

		/* Center at midpoint of (pj2, pk2) */
		mid2[0] = (pj2x + pk2x) / (qaws_scalar)2;
		mid2[1] = (pj2y + pk2y) / (qaws_scalar)2;
		dir2[0] = pk2x - pj2x;
		dir2[1] = pk2y - pj2y;
		len2 = (qaws_scalar)sqrt((double)(dir2[0]*dir2[0] + dir2[1]*dir2[1]));
		if (len2 < (qaws_scalar)1e-10) {
			compute_bezier_subcurve(pts, n, i, closed, dim, sc);
			return;
		}
		half_len2 = len2 / (qaws_scalar)2;

		/* 2D frame: au along pj->pk, av perpendicular */
		au2[0] = dir2[0] / len2;
		au2[1] = dir2[1] / len2;
		av2[0] = -au2[1];
		av2[1] = au2[0];

		/* Project pi (at origin in 2D) into local ellipse frame */
		pirel2[0] = pi2x - mid2[0];
		pirel2[1] = pi2y - mid2[1];
		px2 = pirel2[0] * au2[0] + pirel2[1] * au2[1];
		py2 = pirel2[0] * av2[0] + pirel2[1] * av2[1];
		if (py2 < (qaws_scalar)0) py2 = -py2;

		a2 = half_len2;
		b2 = py2;
		if (b2 < (qaws_scalar)0.01) b2 = (qaws_scalar)0.01;

		for (iter2 = 0; iter2 < 16; iter2++) {
			qaws_scalar x2a2 = (px2 * px2) / (a2 * a2);
			qaws_scalar denom_val = (qaws_scalar)1 - x2a2;
			if (denom_val < (qaws_scalar)0.001) denom_val = (qaws_scalar)0.001;
			b2 = py2 / (qaws_scalar)sqrt((double)denom_val);
			if (b2 < (qaws_scalar)0.01) b2 = (qaws_scalar)0.01;
		}

		/* Map center back to 3D: first, 2D center is mid2 relative to pi */
		center2x = mid2[0];
		center2y = mid2[1];
		for (d = 0; d < 3; d++)
			sc->center[d] = pi[d] + center2x * frame_u[d] + center2y * frame_w[d];

		/* Map axis_u (along pj->pk direction) back to 3D */
		for (d = 0; d < 3; d++)
			sc->axis_u[d] = au2[0] * frame_u[d] + au2[1] * frame_w[d];

		/* Map axis_v (perpendicular in-plane) back to 3D */
		for (d = 0; d < 3; d++)
			sc->axis_v[d] = av2[0] * frame_u[d] + av2[1] * frame_w[d];

		sc->radius = a2;
		sc->radius_b = b2;

		sc->angle_start = QAWS_PI;
		sc->angle_end = (qaws_scalar)0;

		/* Detect which way pi is relative to the chord in 2D */
		{
			qaws_scalar cross2 = dir2[0] * (pi2y - pj2y) - dir2[1] * (pi2x - pj2x);
			if (cross2 < 0) {
				sc->radius_b = -b2;
			}
		}

		sc->use_arc = 1;
		sc->t1 = (qaws_scalar)0.5;
		vcopy(sc->p1, pi, dim);
		return;
	}

	/* 2D case: original code */

	/* Ellipse fitting: center at midpoint of (pj, pk), major axis along pj->pk */
	for (d = 0; d < dim; d++)
		mid[d] = (pj[d] + pk[d]) / (qaws_scalar)2;

	vsub(dir, pk, pj, dim);
	len = vlength(dir, dim);
	if (len < (qaws_scalar)1e-10) {
		compute_bezier_subcurve(pts, n, i, closed, dim, sc);
		return;
	}

	half_len = len / (qaws_scalar)2;

	/* Normalize */
	sc->axis_u[0] = dir[0] / len;
	sc->axis_u[1] = dir[1] / len;
	sc->axis_v[0] = -sc->axis_u[1];
	sc->axis_v[1] = sc->axis_u[0];

	/* Project pi into local frame */
	vsub(pi_rel, pi, mid, dim);
	px = pi_rel[0] * sc->axis_u[0] + pi_rel[1] * sc->axis_u[1];
	py = pi_rel[0] * sc->axis_v[0] + pi_rel[1] * sc->axis_v[1];

	if (py < (qaws_scalar)0) py = -py; /* ensure positive semi-axis */

	/* Semi-major axis a = half distance between pj and pk */
	a = half_len;

	/* Iterative refinement for semi-minor axis b (16 iterations like the reference) */
	/* Start with b = |py| */
	b = py;
	if (b < (qaws_scalar)0.01) b = (qaws_scalar)0.01;

	for (iter = 0; iter < 16; iter++) {
		/* For the ellipse x^2/a^2 + y^2/b^2 = 1, point pi should lie on it */
		/* Adjust b so that the ellipse passes through pi */
		qaws_scalar x2a2 = (px * px) / (a * a);
		qaws_scalar denom_val = (qaws_scalar)1 - x2a2;
		if (denom_val < (qaws_scalar)0.001) denom_val = (qaws_scalar)0.001;
		b = py / (qaws_scalar)sqrt((double)denom_val);
		if (b < (qaws_scalar)0.01) b = (qaws_scalar)0.01;
	}

	sc->center[0] = mid[0];
	sc->center[1] = mid[1];
	sc->radius = a;
	sc->radius_b = b;

	sc->angle_start = (qaws_scalar)atan2((double)0, (double)(-a)); /* pi angle, i.e. pj */
	sc->angle_end = (qaws_scalar)atan2((double)0, (double)a);       /* 0 angle, i.e. pk */

	/* Determine sign of arc angles: pj corresponds to angle PI, pk to angle 0 */
	sc->angle_start = QAWS_PI;
	sc->angle_end = (qaws_scalar)0;

	/* Detect which way pi is relative to the chord */
	{
		qaws_scalar cross = dir[0] * (pi[1] - pj[1]) - dir[1] * (pi[0] - pj[0]);
		if (cross < 0) {
			/* pi is below the chord, negate b */
			sc->radius_b = -b;
		}
	}

	sc->use_arc = 1;
	sc->t1 = (qaws_scalar)0.5;
	vcopy(sc->p1, pi, dim);
}

static void compute_hybrid_subcurve(
	qaws_scalar const* pts, unsigned int n, unsigned int i,
	int closed, unsigned int dim, qaws_yuksel_subcurve* sc)
{
	/* Try circular first */
	compute_circular_subcurve(pts, n, i, closed, dim, sc);

	if (sc->use_arc) {
		/* Check if arc angle is within [-PI/2, PI/2] */
		qaws_scalar angle_span = sc->angle_end - sc->angle_start;
		if (angle_span < -QAWS_PI_2 || angle_span > QAWS_PI_2) {
			/* Fall back to elliptical */
			compute_elliptical_subcurve(pts, n, i, closed, dim, sc);
		}
	}
}

/* ------------------------------------------------------------------ */
/*  Evaluate a single sub-curve at local parameter t in [0,1]         */
/* ------------------------------------------------------------------ */

static void eval_subcurve_bezier(
	qaws_yuksel_subcurve const* sc, qaws_scalar t, unsigned int dim,
	qaws_scalar* out_pos, qaws_scalar* out_d1, qaws_scalar* out_d2)
{
	qaws_scalar omt = (qaws_scalar)1 - t;
	unsigned int d;

	if (out_pos) {
		for (d = 0; d < dim; d++)
			out_pos[d] = omt * omt * sc->p0[d]
				+ (qaws_scalar)2 * omt * t * sc->p1[d]
				+ t * t * sc->p2[d];
	}

	if (out_d1) {
		for (d = 0; d < dim; d++)
			out_d1[d] = (qaws_scalar)2 * ((t - (qaws_scalar)1) * sc->p0[d]
				+ ((qaws_scalar)1 - (qaws_scalar)2 * t) * sc->p1[d]
				+ t * sc->p2[d]);
	}

	if (out_d2) {
		for (d = 0; d < dim; d++)
			out_d2[d] = (qaws_scalar)2 * (sc->p0[d] - (qaws_scalar)2 * sc->p1[d] + sc->p2[d]);
	}
}

static void eval_subcurve_arc(
	qaws_yuksel_subcurve const* sc, qaws_scalar t, unsigned int dim,
	qaws_scalar* out_pos, qaws_scalar* out_d1, qaws_scalar* out_d2)
{
	qaws_scalar angle, ca, sa, ra, rb;
	qaws_scalar da_dt;
	unsigned int d;

	angle = sc->angle_start + t * (sc->angle_end - sc->angle_start);
	da_dt = sc->angle_end - sc->angle_start;
	ca = (qaws_scalar)cos((double)angle);
	sa = (qaws_scalar)sin((double)angle);
	ra = sc->radius;
	rb = sc->radius_b;

	if (out_pos) {
		for (d = 0; d < dim; d++)
			out_pos[d] = sc->center[d] + ra * ca * sc->axis_u[d] + rb * sa * sc->axis_v[d];
	}

	if (out_d1) {
		for (d = 0; d < dim; d++)
			out_d1[d] = (-ra * sa * sc->axis_u[d] + rb * ca * sc->axis_v[d]) * da_dt;
	}

	if (out_d2) {
		for (d = 0; d < dim; d++)
			out_d2[d] = (-ra * ca * sc->axis_u[d] - rb * sa * sc->axis_v[d]) * da_dt * da_dt;
	}
}

static void eval_subcurve(
	qaws_yuksel_subcurve const* sc, qaws_scalar t, unsigned int dim,
	qaws_scalar* out_pos, qaws_scalar* out_d1, qaws_scalar* out_d2)
{
	if (sc->use_arc)
		eval_subcurve_arc(sc, t, dim, out_pos, out_d1, out_d2);
	else
		eval_subcurve_bezier(sc, t, dim, out_pos, out_d1, out_d2);
}

/* ------------------------------------------------------------------ */
/*  Blended span evaluation                                            */
/*  For span between P[i] and P[i+1]:                                 */
/*    curve1 = sub-curve at P[i],   evaluated at t1 + (1-t1)*t        */
/*    curve2 = sub-curve at P[i+1], evaluated at t2 * t               */
/*  P(t) = cos^2(pi*t/2) * curve1 + sin^2(pi*t/2) * curve2           */
/* ------------------------------------------------------------------ */

static void yuksel_blend_eval(
	qaws_yuksel_impl const* impl,
	unsigned int span_index,
	qaws_scalar t,
	unsigned int dim,
	unsigned int eval_flags,
	qaws_scalar* out_pos,
	qaws_scalar* out_d1,
	qaws_scalar* out_d2,
	qaws_scalar* out_d3)
{
	unsigned int n = impl->control_point_count;
	unsigned int i1, i2;
	qaws_yuksel_subcurve const* sc1;
	qaws_yuksel_subcurve const* sc2;
	qaws_scalar t_sc1, t_sc2;
	qaws_scalar c_val, s_val, c2, s2;
	qaws_scalar p1[3], p2[3];
	qaws_scalar dp1[3], dp2[3];
	qaws_scalar d2p1[3], d2p2[3];
	qaws_scalar dt_sc1, dt_sc2;
	unsigned int d;
	int need_d1, need_d2;
	int is_first_span, is_last_span;

	if (impl->closed) {
		i1 = span_index % n;
		i2 = (span_index + 1) % n;
	} else {
		i1 = span_index;
		i2 = span_index + 1;
	}

	sc1 = &impl->subcurves[i1];
	sc2 = &impl->subcurves[i2];

	/* Map global span t to sub-curve local parameters */
	/* sc1: right half, from t1 to 1.0 -> local = t1 + (1-t1)*t */
	/* sc2: left half, from 0.0 to t2 -> local = t2 * t */
	dt_sc1 = (qaws_scalar)1 - sc1->t1;
	dt_sc2 = sc2->t1;
	t_sc1 = sc1->t1 + dt_sc1 * t;
	t_sc2 = dt_sc2 * t;

	/* Check for boundary spans on open curves */
	is_first_span = (!impl->closed && span_index == 0);
	is_last_span = (!impl->closed && span_index == n - 2);

	need_d1 = (eval_flags & (QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3)) != 0;
	need_d2 = (eval_flags & (QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3)) != 0;

	/* Evaluate both sub-curves */
	memset(p1, 0, sizeof(p1));
	memset(p2, 0, sizeof(p2));
	memset(dp1, 0, sizeof(dp1));
	memset(dp2, 0, sizeof(dp2));
	memset(d2p1, 0, sizeof(d2p1));
	memset(d2p2, 0, sizeof(d2p2));

	eval_subcurve(sc1, t_sc1, dim, p1, need_d1 ? dp1 : NULL, need_d2 ? d2p1 : NULL);
	eval_subcurve(sc2, t_sc2, dim, p2, need_d1 ? dp2 : NULL, need_d2 ? d2p2 : NULL);

	/* Chain rule: dp1/dt = dp1/dt_sc1 * dt_sc1/dt */
	if (need_d1) {
		for (d = 0; d < dim; d++) {
			dp1[d] *= dt_sc1;
			dp2[d] *= dt_sc2;
		}
	}
	if (need_d2) {
		for (d = 0; d < dim; d++) {
			d2p1[d] *= dt_sc1 * dt_sc1;
			d2p2[d] *= dt_sc2 * dt_sc2;
		}
	}

	/* For first/last open-curve spans, use single sub-curve unblended */
	if (is_first_span && !impl->closed) {
		/* Use sc2 only (the curve at P[1] from its left half) */
		if (out_pos) { for (d = 0; d < dim; d++) out_pos[d] = p2[d]; }
		if (out_d1 && need_d1) { for (d = 0; d < dim; d++) out_d1[d] = dp2[d]; }
		if (out_d2 && need_d2) { for (d = 0; d < dim; d++) out_d2[d] = d2p2[d]; }
		if (out_d3) { for (d = 0; d < dim; d++) out_d3[d] = (qaws_scalar)0; }
		return;
	}

	if (is_last_span && !impl->closed) {
		/* Use sc1 only (the curve at P[n-2] from its right half) */
		if (out_pos) { for (d = 0; d < dim; d++) out_pos[d] = p1[d]; }
		if (out_d1 && need_d1) { for (d = 0; d < dim; d++) out_d1[d] = dp1[d]; }
		if (out_d2 && need_d2) { for (d = 0; d < dim; d++) out_d2[d] = d2p1[d]; }
		if (out_d3) { for (d = 0; d < dim; d++) out_d3[d] = (qaws_scalar)0; }
		return;
	}

	/* Blending weights: w1 = cos^2(pi*t/2), w2 = sin^2(pi*t/2) */
	c_val = (qaws_scalar)cos((double)(QAWS_PI_2 * t));
	s_val = (qaws_scalar)sin((double)(QAWS_PI_2 * t));
	c2 = c_val * c_val;
	s2 = s_val * s_val;

	/* Position: P = c2*p1 + s2*p2 */
	if (out_pos) {
		for (d = 0; d < dim; d++)
			out_pos[d] = c2 * p1[d] + s2 * p2[d];
	}

	/* First derivative via product rule:
	   dw1/dt = -pi * cos(pi*t/2) * sin(pi*t/2) = -pi/2 * sin(pi*t)
	   dw2/dt =  pi * cos(pi*t/2) * sin(pi*t/2) =  pi/2 * sin(pi*t)
	   D1 = dw1*p1 + w1*dp1 + dw2*p2 + w2*dp2 */
	if (out_d1 && need_d1) {
		qaws_scalar dw1 = -QAWS_PI * c_val * s_val;
		qaws_scalar dw2 =  QAWS_PI * c_val * s_val;
		for (d = 0; d < dim; d++)
			out_d1[d] = dw1 * p1[d] + c2 * dp1[d]
				+ dw2 * p2[d] + s2 * dp2[d];
	}

	/* Second derivative */
	if (out_d2 && need_d2) {
		qaws_scalar dw1 = -QAWS_PI * c_val * s_val;
		qaws_scalar dw2 =  QAWS_PI * c_val * s_val;
		/* d2w1/dt2 = -pi^2/2 * cos(pi*t) */
		qaws_scalar d2w1 = -QAWS_PI * QAWS_PI / (qaws_scalar)2
			* (qaws_scalar)cos((double)(QAWS_PI * t));
		qaws_scalar d2w2 =  QAWS_PI * QAWS_PI / (qaws_scalar)2
			* (qaws_scalar)cos((double)(QAWS_PI * t));
		for (d = 0; d < dim; d++)
			out_d2[d] = d2w1 * p1[d] + (qaws_scalar)2 * dw1 * dp1[d] + c2 * d2p1[d]
				+ d2w2 * p2[d] + (qaws_scalar)2 * dw2 * dp2[d] + s2 * d2p2[d];
	}

	/* D3: set to zero for now (sub-curves are quadratic, so d3 of sub-curve is 0,
	   but weight d3 is not. For practical use this is sufficient.) */
	if (out_d3) {
		for (d = 0; d < dim; d++)
			out_d3[d] = (qaws_scalar)0;
	}
}

/* ------------------------------------------------------------------ */
/*  Vtable functions                                                   */
/* ------------------------------------------------------------------ */

static qaws_status yuksel_eval_span_2d(
	qaws_curve const* curve,
	unsigned int span_index,
	qaws_scalar local_t,
	unsigned int eval_flags,
	qaws_eval_result_2d* out_result)
{
	qaws_yuksel_impl const* impl = (qaws_yuksel_impl const*)curve->impl;
	qaws_scalar pos[3] = {0}, d1[3] = {0}, d2[3] = {0}, d3[3] = {0};

	memset(out_result, 0, sizeof(*out_result));

	yuksel_blend_eval(impl, span_index, local_t, 2, eval_flags,
		(eval_flags & QAWS_EVAL_FLAG_POSITION) ? pos : NULL,
		(eval_flags & (QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3)) ? d1 : NULL,
		(eval_flags & (QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3)) ? d2 : NULL,
		(eval_flags & QAWS_EVAL_FLAG_D3) ? d3 : NULL);

	if (eval_flags & QAWS_EVAL_FLAG_POSITION) {
		out_result->position.x = pos[0];
		out_result->position.y = pos[1];
		out_result->valid_flags |= QAWS_EVAL_FLAG_POSITION;
	}
	if (eval_flags & QAWS_EVAL_FLAG_D1) {
		out_result->d1.x = d1[0];
		out_result->d1.y = d1[1];
		out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
	}
	if (eval_flags & QAWS_EVAL_FLAG_D2) {
		out_result->d2.x = d2[0];
		out_result->d2.y = d2[1];
		out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
	}
	if (eval_flags & QAWS_EVAL_FLAG_D3) {
		out_result->d3.x = d3[0];
		out_result->d3.y = d3[1];
		out_result->valid_flags |= QAWS_EVAL_FLAG_D3;
	}

	return QAWS_STATUS_OK;
}

static qaws_status yuksel_eval_span_3d(
	qaws_curve const* curve,
	unsigned int span_index,
	qaws_scalar local_t,
	unsigned int eval_flags,
	qaws_eval_result_3d* out_result)
{
	qaws_yuksel_impl const* impl = (qaws_yuksel_impl const*)curve->impl;
	qaws_scalar pos[3] = {0}, d1[3] = {0}, d2[3] = {0}, d3[3] = {0};

	memset(out_result, 0, sizeof(*out_result));

	yuksel_blend_eval(impl, span_index, local_t, 3, eval_flags,
		(eval_flags & QAWS_EVAL_FLAG_POSITION) ? pos : NULL,
		(eval_flags & (QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3)) ? d1 : NULL,
		(eval_flags & (QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3)) ? d2 : NULL,
		(eval_flags & QAWS_EVAL_FLAG_D3) ? d3 : NULL);

	if (eval_flags & QAWS_EVAL_FLAG_POSITION) {
		out_result->position.x = pos[0];
		out_result->position.y = pos[1];
		out_result->position.z = pos[2];
		out_result->valid_flags |= QAWS_EVAL_FLAG_POSITION;
	}
	if (eval_flags & QAWS_EVAL_FLAG_D1) {
		out_result->d1.x = d1[0];
		out_result->d1.y = d1[1];
		out_result->d1.z = d1[2];
		out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
	}
	if (eval_flags & QAWS_EVAL_FLAG_D2) {
		out_result->d2.x = d2[0];
		out_result->d2.y = d2[1];
		out_result->d2.z = d2[2];
		out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
	}
	if (eval_flags & QAWS_EVAL_FLAG_D3) {
		out_result->d3.x = d3[0];
		out_result->d3.y = d3[1];
		out_result->d3.z = d3[2];
		out_result->valid_flags |= QAWS_EVAL_FLAG_D3;
	}

	return QAWS_STATUS_OK;
}

static void yuksel_destroy_impl(void* impl_ptr, qaws_allocator const* allocator)
{
	qaws_yuksel_impl* impl = (qaws_yuksel_impl*)impl_ptr;
	if (impl) {
		qaws_internal_dealloc(allocator, impl->control_points);
		qaws_internal_dealloc(allocator, impl->subcurves);
		qaws_internal_dealloc(allocator, impl);
	}
}

static int yuksel_is_closed(qaws_curve const* curve)
{
	qaws_yuksel_impl const* impl = (qaws_yuksel_impl const*)curve->impl;
	return impl->closed;
}

static int yuksel_is_periodic(qaws_curve const* curve)
{
	qaws_yuksel_impl const* impl = (qaws_yuksel_impl const*)curve->impl;
	return impl->closed;
}

static int yuksel_is_rational(qaws_curve const* curve)
{
	(void)curve;
	return 0;
}

static qaws_continuity yuksel_get_continuity(qaws_curve const* curve)
{
	(void)curve;
	return QAWS_CONTINUITY_C2;
}

/* ------------------------------------------------------------------ */
/*  Vtable                                                             */
/* ------------------------------------------------------------------ */

static qaws_curve_vtable const yuksel_vtable = {
	yuksel_eval_span_2d,
	yuksel_eval_span_3d,
	yuksel_destroy_impl,
	yuksel_is_closed,
	yuksel_is_periodic,
	yuksel_is_rational,
	yuksel_get_continuity,
	NULL /* diff */
};

/* ------------------------------------------------------------------ */
/*  Creation function                                                  */
/* ------------------------------------------------------------------ */

qaws_status qaws_curve_create_yuksel(
	qaws_yuksel_desc const* desc,
	qaws_curve** out_curve)
{
	qaws_status status;
	qaws_curve* curve;
	qaws_yuksel_impl* impl;
	unsigned int dim_count;
	unsigned int span_count;
	unsigned int n;
	size_t pts_size;
	qaws_range range;
	unsigned int i;

	if (!desc || !out_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_curve = NULL;

	status = qaws_internal_validate_dimension(desc->dimension);
	if (status != QAWS_STATUS_OK)
		return status;

	n = desc->control_point_count;

	if (n < 3)
		return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;

	if (!desc->control_points)
		return QAWS_STATUS_INVALID_ARGUMENT;

	dim_count = (desc->dimension == QAWS_DIMENSION_2D) ? 2u : 3u;

	if (desc->closed)
		span_count = n;
	else
		span_count = n - 1;

	range.min_value = (qaws_scalar)0;
	range.max_value = (qaws_scalar)span_count;

	curve = qaws_internal_curve_alloc(
		QAWS_CURVE_KIND_YUKSEL,
		desc->dimension,
		2, /* degree (quadratic sub-curves) */
		span_count,
		range,
		&yuksel_vtable);
	if (!curve)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	for (i = 0; i <= span_count; i++)
		curve->span_boundaries[i] = (qaws_scalar)i;

	impl = (qaws_yuksel_impl*)malloc(sizeof(qaws_yuksel_impl));
	if (!impl) {
		qaws_internal_curve_free(curve);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	pts_size = (size_t)n * (size_t)dim_count * sizeof(qaws_scalar);
	impl->control_points = (qaws_scalar*)malloc(pts_size);
	if (!impl->control_points) {
		free(impl);
		qaws_internal_curve_free(curve);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memcpy(impl->control_points, desc->control_points, pts_size);

	impl->control_point_count = n;
	impl->closed = desc->closed;
	impl->mode = (int)desc->mode;

	impl->subcurves = (qaws_yuksel_subcurve*)calloc(n, sizeof(qaws_yuksel_subcurve));
	if (!impl->subcurves) {
		free(impl->control_points);
		free(impl);
		qaws_internal_curve_free(curve);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	/* Compute sub-curves for each control point */
	for (i = 0; i < n; i++) {
		switch (desc->mode) {
		case QAWS_YUKSEL_MODE_BEZIER:
			compute_bezier_subcurve(impl->control_points, n, i, desc->closed, dim_count, &impl->subcurves[i]);
			break;
		case QAWS_YUKSEL_MODE_CIRCULAR:
			compute_circular_subcurve(impl->control_points, n, i, desc->closed, dim_count, &impl->subcurves[i]);
			break;
		case QAWS_YUKSEL_MODE_ELLIPTICAL:
			compute_elliptical_subcurve(impl->control_points, n, i, desc->closed, dim_count, &impl->subcurves[i]);
			break;
		case QAWS_YUKSEL_MODE_HYBRID:
			compute_hybrid_subcurve(impl->control_points, n, i, desc->closed, dim_count, &impl->subcurves[i]);
			break;
		default:
			compute_bezier_subcurve(impl->control_points, n, i, desc->closed, dim_count, &impl->subcurves[i]);
			break;
		}
	}

	curve->impl = impl;
	*out_curve = curve;

	return QAWS_STATUS_OK;
}
