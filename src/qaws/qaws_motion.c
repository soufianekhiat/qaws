/*
 * Motions (qaws_motion.h): Hermite polynomials for positions, the Cayley
 * map of a Hermite polynomial in the Gibbs vector for rotations.
 */

#include "qaws_motion.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "internal/qaws_internal_memory.h"
#include "internal/qaws_internal_array.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MO_POSES 0
#define MO_CURVE 1
#define MO_VERTICES 2
#define MO_PI 3.14159265358979323846

/* one piece of a pose motion over its own t in [0, 1] */
typedef struct mo_seg
{
	double R0[9];              /* the start's rotation */
	double T[6][3];            /* the translation's p0, v0, a0, p1, v1, a1 */
	double C[6][3];            /* the Gibbs vector's c0 (0), c0', c0'', c1, c1', c1'' */
} mo_seg;

struct qaws_motion
{
	int kind;
	int quintic;
	double R0[9];              /* CURVE: the constant rotation */
	mo_seg* seg;               /* POSES: its pieces, each over an equal share of t (one unless it turns past a half turn) */
	unsigned int nseg;
	qaws_curve const* curve;   /* CURVE: borrowed */
	double u0, u1;             /* CURVE: its parameter range */
	unsigned int n;            /* VERTICES: their count */
	double* V;                 /* VERTICES: 18 each: p0, v0, a0, p1, v1, a1 (NULL: borrowed, P0 to P1 on chords) */
	qaws_array_view P0, P1;    /* VERTICES borrowed: the host's positions at 0 and 1 */
};

/* the Hermite basis at t (cubic: 4 functions on p0, v0, p1, v1, written into slots 0, 1, 3, 4; quintic: 6 on p0, v0, a0,
   p1, v1, a1) and its derivatives */
static void mo_basis(double t, int quintic, double* h, double* dh)
{
	double t2 = t * t, t3 = t2 * t, t4 = t3 * t, t5 = t4 * t;
	if (!quintic)
	{
		h[0] = 2 * t3 - 3 * t2 + 1;  dh[0] = 6 * t2 - 6 * t;
		h[1] = t3 - 2 * t2 + t;      dh[1] = 3 * t2 - 4 * t + 1;
		h[2] = 0;                    dh[2] = 0;
		h[3] = -2 * t3 + 3 * t2;     dh[3] = -6 * t2 + 6 * t;
		h[4] = t3 - t2;              dh[4] = 3 * t2 - 2 * t;
		h[5] = 0;                    dh[5] = 0;
		return;
	}
	h[0] = 1 - 10 * t3 + 15 * t4 - 6 * t5;                    dh[0] = -30 * t2 + 60 * t3 - 30 * t4;
	h[1] = t - 6 * t3 + 8 * t4 - 3 * t5;                      dh[1] = 1 - 18 * t2 + 32 * t3 - 15 * t4;
	h[2] = 0.5 * t2 - 1.5 * t3 + 1.5 * t4 - 0.5 * t5;         dh[2] = t - 4.5 * t2 + 6 * t3 - 2.5 * t4;
	h[3] = 10 * t3 - 15 * t4 + 6 * t5;                        dh[3] = 30 * t2 - 60 * t3 + 30 * t4;
	h[4] = -4 * t3 + 7 * t4 - 3 * t5;                         dh[4] = -12 * t2 + 28 * t3 - 15 * t4;
	h[5] = 0.5 * t3 - t4 + 0.5 * t5;                          dh[5] = 1.5 * t2 - 4 * t3 + 2.5 * t4;
}

/* the Hermite basis' second derivatives at t (slots as mo_basis) */
static void mo_basis2(double t, int quintic, double* d2)
{
	double t2 = t * t, t3 = t2 * t;
	if (!quintic)
	{
		d2[0] = 12 * t - 6;
		d2[1] = 6 * t - 4;
		d2[2] = 0;
		d2[3] = -12 * t + 6;
		d2[4] = 6 * t - 2;
		d2[5] = 0;
		return;
	}
	d2[0] = -60 * t + 180 * t2 - 120 * t3;
	d2[1] = -36 * t + 96 * t2 - 60 * t3;
	d2[2] = 1 - 9 * t + 18 * t2 - 10 * t3;
	d2[3] = 60 * t - 180 * t2 + 120 * t3;
	d2[4] = -24 * t + 84 * t2 - 60 * t3;
	d2[5] = 3 * t - 12 * t2 + 10 * t3;
}

/* the Hermite polynomial of the six data (p0, v0, a0, p1, v1, a1, 3 each) at t, and its derivative */
static void mo_hermite(double const (*D)[3], double t, int quintic, double* x, double* dx)
{
	double h[6], dh[6];
	unsigned int k, j;
	mo_basis(t, quintic, h, dh);
	for (k = 0; k < 3; k++)
	{
		x[k] = 0;
		dx[k] = 0;
		for (j = 0; j < 6; j++)
		{
			x[k] += h[j] * D[j][k];
			dx[k] += dh[j] * D[j][k];
		}
	}
}

static void mo_quat_matrix(double const* q, double* R)
{
	double w = q[0], x = q[1], y = q[2], z = q[3];
	R[0] = 1 - 2 * (y * y + z * z); R[1] = 2 * (x * y - w * z);     R[2] = 2 * (x * z + w * y);
	R[3] = 2 * (x * y + w * z);     R[4] = 1 - 2 * (x * x + z * z); R[5] = 2 * (y * z - w * x);
	R[6] = 2 * (x * z - w * y);     R[7] = 2 * (y * z + w * x);     R[8] = 1 - 2 * (x * x + y * y);
}

/* q normalized into out; 0 for a zero or non-finite quaternion */
static int mo_unit(qaws_scalar const* q, double* out)
{
	double l = 0;
	unsigned int k;
	for (k = 0; k < 4; k++)
	{
		out[k] = (double)q[k];
		if (!isfinite(out[k]))
			return 0;
		l += out[k] * out[k];
	}
	l = sqrt(l);
	if (!(l > 0))
		return 0;
	for (k = 0; k < 4; k++) out[k] /= l;
	return 1;
}

static void mo_cross(double const* a, double const* b, double* o)
{
	o[0] = a[1] * b[2] - a[2] * b[1];
	o[1] = a[2] * b[0] - a[0] * b[2];
	o[2] = a[0] * b[1] - a[1] * b[0];
}

static double mo_dot(double const* a, double const* b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/* the Gibbs vector's derivatives at c from the angular velocity w and acceleration dw (the rotation's own frame of pose 0):
   c' = (w - c x w + (c . w) c) / 2, inverting w = 2 (c' + c x c') / (1 + c . c), and c'' its derivative */
static void mo_gibbs_rates(double const* c, double const* w, double const* dw, double* dc, double* ddc)
{
	double cw[3], dcw[3], cdw[3], cdotw = mo_dot(c, w);
	unsigned int k;
	mo_cross(c, w, cw);
	for (k = 0; k < 3; k++) dc[k] = 0.5 * (w[k] - cw[k] + cdotw * c[k]);
	mo_cross(dc, w, dcw);
	mo_cross(c, dw, cdw);
	for (k = 0; k < 3; k++)
		ddc[k] = 0.5 * (dw[k] - dcw[k] - cdw[k] + (mo_dot(dc, w) + mo_dot(c, dw)) * c[k] + cdotw * dc[k]);
}

static int mo_finite3(qaws_vec3 v)
{
	return isfinite((double)v.x) && isfinite((double)v.y) && isfinite((double)v.z);
}

static void mo_vec(qaws_vec3 v, double* o)
{
	o[0] = (double)v.x;
	o[1] = (double)v.y;
	o[2] = (double)v.z;
}

/* R^T v */
static void mo_tmul(double const* R, double const* v, double* o)
{
	unsigned int k;
	for (k = 0; k < 3; k++) o[k] = R[k] * v[0] + R[3 + k] * v[1] + R[6 + k] * v[2];
}

/* q a b (Hamilton) */
static void mo_qmul(double const* a, double const* b, double* o)
{
	o[0] = a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3];
	o[1] = a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2];
	o[2] = a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1];
	o[3] = a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0];
}

/* one piece from its ends (q unit, the relative rotation under a half turn; T the translation's six data, w and a the
   angular velocities and accelerations, world frame, per unit of the piece's t); 0 for a half turn */
static int mo_piece(mo_seg* g, double const* q0, double const* q1, double const (*T)[3], double const* w0, double const* w1, double const* a0,
	double const* a1, int quintic)
{
	double qc[4], qr[4], wc[3], ac[3], zero[3] = { 0, 0, 0 };
	unsigned int k;
	/* the relative rotation q0^-1 q1, the shorter way round */
	qc[0] = q0[0];
	qc[1] = -q0[1];
	qc[2] = -q0[2];
	qc[3] = -q0[3];
	mo_qmul(qc, q1, qr);
	if (qr[0] < 0)
		for (k = 0; k < 4; k++) qr[k] = -qr[k];
	if (!(qr[0] > 1e-9))
		return 0;   /* (a half turn: no Gibbs vector) */
	mo_quat_matrix(q0, g->R0);
	memcpy(g->T, T, sizeof g->T);
	/* the Gibbs vector's ends: 0 and the relative rotation's; their rates from the angular velocities in pose 0's frame */
	memset(g->C, 0, sizeof g->C);
	for (k = 0; k < 3; k++) g->C[3][k] = qr[1 + k] / qr[0];
	mo_tmul(g->R0, w0, wc);
	mo_tmul(g->R0, a0, ac);
	mo_gibbs_rates(zero, wc, ac, g->C[1], g->C[2]);
	mo_tmul(g->R0, w1, wc);
	mo_tmul(g->R0, a1, ac);
	mo_gibbs_rates(g->C[3], wc, ac, g->C[4], g->C[5]);
	if (!quintic)
		for (k = 0; k < 3; k++) g->C[2][k] = g->C[5][k] = 0;   /* (the cubic uses neither acceleration) */
	return 1;
}

qaws_status qaws_motion_create_poses(qaws_motion** out, qaws_motion_pose const* s, qaws_motion_pose const* e, unsigned int flags)
{
	qaws_motion* m;
	double q0[4], q1[4], qc[4], qr[4], T[6][3], w0[3], w1[3], a0[3] = { 0, 0, 0 }, a1[3] = { 0, 0, 0 }, u[3], wm[3], lu, th, total, phi;
	int quintic = (flags & QAWS_MOTION_QUINTIC) != 0;
	unsigned int k, i, n = 1;
	if (!out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	if (!s || !e || !mo_unit(s->rotation, q0) || !mo_unit(e->rotation, q1) || !mo_finite3(s->position) || !mo_finite3(e->position) ||
		!mo_finite3(s->velocity) || !mo_finite3(e->velocity) || !mo_finite3(s->angular_velocity) || !mo_finite3(e->angular_velocity) ||
		(quintic && (!mo_finite3(s->acceleration) || !mo_finite3(e->acceleration) || !mo_finite3(s->angular_acceleration) ||
		!mo_finite3(e->angular_acceleration))))
		return QAWS_STATUS_INVALID_ARGUMENT;
	memset(T, 0, sizeof T);
	mo_vec(s->position, T[0]);
	mo_vec(s->velocity, T[1]);
	mo_vec(e->position, T[3]);
	mo_vec(e->velocity, T[4]);
	if (quintic)
	{
		mo_vec(s->acceleration, T[2]);
		mo_vec(e->acceleration, T[5]);
		mo_vec(s->angular_acceleration, a0);
		mo_vec(e->angular_acceleration, a1);
	}
	mo_vec(s->angular_velocity, w0);
	mo_vec(e->angular_velocity, w1);
	/* how far it turns: the relative rotation (angle th in [0, pi] about u, world frame) plus whole turns, the total nearest
	   the angular velocities' mean along u (a spin of three turns ends where it started: its velocities say so) */
	qc[0] = q0[0];
	qc[1] = -q0[1];
	qc[2] = -q0[2];
	qc[3] = -q0[3];
	mo_qmul(q1, qc, qr);   /* (q1 q0^-1: the rotation in the world frame) */
	if (qr[0] < 0)
		for (k = 0; k < 4; k++) qr[k] = -qr[k];
	lu = sqrt(qr[1] * qr[1] + qr[2] * qr[2] + qr[3] * qr[3]);
	th = 2 * atan2(lu, qr[0]);
	for (k = 0; k < 3; k++) wm[k] = 0.5 * (w0[k] + w1[k]);
	if (lu > 1e-12)
		for (k = 0; k < 3; k++) u[k] = qr[1 + k] / lu;
	else
	{
		double l = sqrt(mo_dot(wm, wm));
		for (k = 0; k < 3; k++) u[k] = l > 0 ? wm[k] / l : (k == 2);
		th = 0;
	}
	phi = mo_dot(wm, u);
	total = th + 2 * MO_PI * floor((phi - th) / (2 * MO_PI) + 0.5);
	if (fabs(fabs(phi - th) - MO_PI) < 1e-9 * (1 + fabs(phi)) && th > MO_PI - 1e-6)
		return QAWS_STATUS_INVALID_ARGUMENT;   /* (a half turn the velocities do not orient: either way round) */
	if (fabs(total - th) > 1e-12 || !(cos(0.5 * th) > 1e-9))
		n = (unsigned int)ceil(fabs(total) / (0.25 * MO_PI));   /* (past a half turn: pieces of an eighth of a turn at most; under one, as given) */
	if (n > 1u << 20)
		return QAWS_STATUS_INVALID_ARGUMENT;
	m = (qaws_motion*)qaws_internal_zalloc(1, sizeof *m);
	if (!m)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	m->seg = (mo_seg*)qaws_internal_zalloc(n, sizeof(mo_seg));
	if (!m->seg)
	{
		qaws_internal_free(m);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	m->kind = MO_POSES;
	m->quintic = quintic;
	m->nseg = n;
	if (n == 1)
	{
		if (!mo_piece(m->seg, q0, q1, (double const(*)[3])T, w0, w1, a0, a1, quintic))
		{
			qaws_motion_destroy(m);
			return QAWS_STATUS_INVALID_ARGUMENT;
		}
		*out = m;
		return QAWS_STATUS_OK;
	}
	{
		/* pieces: the turn's angle along u a Hermite polynomial (0 to total, its rates the angular velocities' along u), the
		   angular velocities' other parts blended linearly; the translation the same polynomial cut at the joints (exact) */
		double F[6][3], qa[4], wa[3], aa[3], Ta[6][3];
		memset(F, 0, sizeof F);
		F[1][0] = mo_dot(w0, u);
		F[3][0] = total;
		F[4][0] = mo_dot(w1, u);
		F[2][0] = mo_dot(a0, u);
		F[5][0] = mo_dot(a1, u);
		for (i = 0; i <= n; i++)
		{
			double t = (double)i / n, f[3], df[3], ddf, x[3], dx[3], ddx[3], ddh[6], qb[4], qz[4], wb[3], ab[3], Tb[6][3];
			mo_hermite((double const(*)[3])F, t, quintic, f, df);
			mo_hermite((double const(*)[3])T, t, quintic, x, dx);
			/* second derivatives by the basis' (quintic: matched at the joints) */
			mo_basis2(t, quintic, ddh);
			ddf = 0;
			for (k = 0; k < 3; k++) ddx[k] = 0;
			for (k = 0; k < 6; k++)
			{
				unsigned int z;
				ddf += ddh[k] * F[k][0];
				for (z = 0; z < 3; z++) ddx[z] += ddh[k] * T[k][z];
			}
			/* the orientation: q0 turned by f about u */
			qz[0] = cos(0.5 * f[0]);
			for (k = 0; k < 3; k++) qz[1 + k] = sin(0.5 * f[0]) * u[k];
			mo_qmul(qz, q0, qb);
			if (i == n)
				memcpy(qb, q1, sizeof qb);   /* (the end exactly as given) */
			for (k = 0; k < 3; k++)
			{
				double wl = (1 - t) * w0[k] + t * w1[k], al = (1 - t) * a0[k] + t * a1[k];
				wb[k] = (wl - u[k] * ((1 - t) * mo_dot(w0, u) + t * mo_dot(w1, u)) + u[k] * df[0]) / n;
				ab[k] = (al - u[k] * ((1 - t) * mo_dot(a0, u) + t * mo_dot(a1, u)) + u[k] * ddf) / ((double)n * n);
				Tb[0][k] = x[k];
				Tb[1][k] = dx[k] / n;
				Tb[2][k] = ddx[k] / ((double)n * n);
			}
			if (i > 0)
			{
				double TT[6][3];
				for (k = 0; k < 3; k++)
				{
					unsigned int z;
					for (z = 0; z < 3; z++)
					{
						TT[z][k] = Ta[z][k];
						TT[3 + z][k] = Tb[z][k];
					}
				}
				if (!mo_piece(m->seg + i - 1, qa, qb, (double const(*)[3])TT, wa, wb, aa, ab, quintic))
				{
					qaws_motion_destroy(m);
					return QAWS_STATUS_INVALID_ARGUMENT;
				}
			}
			memcpy(qa, qb, sizeof qa);
			memcpy(wa, wb, sizeof wa);
			memcpy(aa, ab, sizeof aa);
			memcpy(Ta, Tb, sizeof Ta);
		}
	}
	*out = m;
	return QAWS_STATUS_OK;
}

qaws_status qaws_motion_create_curve(qaws_motion** out, qaws_curve const* curve, qaws_scalar const* rotation)
{
	qaws_motion* m;
	double q[4] = { 1, 0, 0, 0 };
	qaws_range rg;
	if (!out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	if (!curve || qaws_curve_get_dimension(curve) != QAWS_DIMENSION_3D || (rotation && !mo_unit(rotation, q)))
		return QAWS_STATUS_INVALID_ARGUMENT;
	rg = qaws_curve_get_parameter_range(curve);
	m = (qaws_motion*)qaws_internal_zalloc(1, sizeof *m);
	if (!m)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	m->kind = MO_CURVE;
	m->curve = curve;
	m->u0 = (double)rg.min_value;
	m->u1 = (double)rg.max_value;
	mo_quat_matrix(q, m->R0);
	*out = m;
	return QAWS_STATUS_OK;
}

qaws_status qaws_motion_create_vertices_borrowed(qaws_motion** out, qaws_array_view const* p0, qaws_array_view const* p1)
{
	qaws_motion* m;
	unsigned int i, k;
	if (!out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	if (!p0 || !p1 || qaws_array_view_validate(p0) != QAWS_STATUS_OK || qaws_array_view_validate(p1) != QAWS_STATUS_OK || p0->components != 3 ||
		p1->components != 3 || !p0->count || p0->count != p1->count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < p0->count; i++)
		for (k = 0; k < 3; k++)
			if (!isfinite(qaws_internal_array_get(p0, i, k)) || !isfinite(qaws_internal_array_get(p1, i, k)))
				return QAWS_STATUS_INVALID_ARGUMENT;
	m = (qaws_motion*)qaws_internal_zalloc(1, sizeof *m);
	if (!m)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	m->kind = MO_VERTICES;
	m->n = p0->count;
	m->P0 = *p0;
	m->P1 = *p1;
	*out = m;
	return QAWS_STATUS_OK;
}

qaws_status qaws_motion_create_vertices(qaws_motion** out, qaws_scalar const* p0, qaws_scalar const* p1, qaws_scalar const* v0,
	qaws_scalar const* v1, qaws_scalar const* a0, qaws_scalar const* a1, unsigned int count, unsigned int flags)
{
	qaws_motion* m;
	unsigned int i, k;
	int quintic = (flags & QAWS_MOTION_QUINTIC) != 0;
	if (!out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	if (!count || !p0 || !p1 || !v0 || !v1)
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < 3 * count; i++)
		if (!isfinite((double)p0[i]) || !isfinite((double)p1[i]) || !isfinite((double)v0[i]) || !isfinite((double)v1[i]) ||
			(quintic && ((a0 && !isfinite((double)a0[i])) || (a1 && !isfinite((double)a1[i])))))
			return QAWS_STATUS_INVALID_ARGUMENT;
	m = (qaws_motion*)qaws_internal_zalloc(1, sizeof *m);
	if (!m)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	m->V = (double*)qaws_internal_zalloc(18 * (size_t)count, sizeof(double));
	if (!m->V)
	{
		qaws_internal_free(m);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	m->kind = MO_VERTICES;
	m->quintic = quintic;
	m->n = count;
	for (i = 0; i < count; i++)
		for (k = 0; k < 3; k++)
		{
			double* d = m->V + 18 * (size_t)i;
			d[k] = (double)p0[3 * (size_t)i + k];
			d[3 + k] = (double)v0[3 * (size_t)i + k];
			d[6 + k] = quintic && a0 ? (double)a0[3 * (size_t)i + k] : 0;
			d[9 + k] = (double)p1[3 * (size_t)i + k];
			d[12 + k] = (double)v1[3 * (size_t)i + k];
			d[15 + k] = quintic && a1 ? (double)a1[3 * (size_t)i + k] : 0;
		}
	*out = m;
	return QAWS_STATUS_OK;
}

void qaws_motion_destroy(qaws_motion* m)
{
	if (!m)
		return;
	qaws_internal_free(m->V);
	qaws_internal_free(m->seg);
	qaws_internal_free(m);
}

int qaws_motion_is_per_vertex(qaws_motion const* m)
{
	return m && m->kind == MO_VERTICES;
}

unsigned int qaws_motion_get_vertex_count(qaws_motion const* m)
{
	return m && m->kind == MO_VERTICES ? m->n : 0;
}

/* a rigid motion at t, in double: R (by rows), T, T', the angular velocity */
static qaws_status mo_rigid(qaws_motion const* m, double t, double* R, double* T, double* dT, double* w)
{
	unsigned int k, j;
	if (m->kind == MO_CURVE)
	{
		qaws_eval_result_3d r;
		qaws_status st = qaws_curve_evaluate_3d(&r, m->curve, (qaws_scalar)(m->u0 + t * (m->u1 - m->u0)), QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1);
		if (st != QAWS_STATUS_OK)
			return st;
		memcpy(R, m->R0, sizeof(double) * 9);
		mo_vec(r.position, T);
		mo_vec(r.d1, dT);
		for (k = 0; k < 3; k++)
		{
			dT[k] *= m->u1 - m->u0;
			w[k] = 0;
		}
		return QAWS_STATUS_OK;
	}
	{
		double c[3], dc[3], cc, A[9], K[9], K2[9], cxdc[3], wc[3], sc = m->nseg;
		unsigned int i = t <= 0 ? 0 : t >= 1 ? m->nseg - 1 : (unsigned int)(t * m->nseg);
		mo_seg const* g;
		if (i >= m->nseg)
			i = m->nseg - 1;
		g = m->seg + i;
		t = t * sc - i;   /* (the piece's own t; t outside [0, 1] extends the end pieces) */
		mo_hermite((double const(*)[3])g->T, t, m->quintic, T, dT);
		mo_hermite((double const(*)[3])g->C, t, m->quintic, c, dc);
		for (k = 0; k < 3; k++)
		{
			dT[k] *= sc;
			dc[k] *= sc;
		}
		/* Cay(c) = I + 2 ([c]x + [c]x^2) / (1 + c . c), R = R0 Cay(c) */
		cc = mo_dot(c, c);
		K[0] = 0;     K[1] = -c[2]; K[2] = c[1];
		K[3] = c[2];  K[4] = 0;     K[5] = -c[0];
		K[6] = -c[1]; K[7] = c[0];  K[8] = 0;
		for (k = 0; k < 3; k++)
			for (j = 0; j < 3; j++)
				K2[3 * k + j] = K[3 * k] * K[j] + K[3 * k + 1] * K[3 + j] + K[3 * k + 2] * K[6 + j];
		for (k = 0; k < 9; k++) A[k] = (k % 4 == 0 ? 1 : 0) + 2 * (K[k] + K2[k]) / (1 + cc);
		for (k = 0; k < 3; k++)
			for (j = 0; j < 3; j++)
				R[3 * k + j] = g->R0[3 * k] * A[j] + g->R0[3 * k + 1] * A[3 + j] + g->R0[3 * k + 2] * A[6 + j];
		/* the angular velocity: R0 (2 / (1 + c . c)) (c' + c x c') */
		mo_cross(c, dc, cxdc);
		for (k = 0; k < 3; k++) wc[k] = 2 * (dc[k] + cxdc[k]) / (1 + cc);
		for (k = 0; k < 3; k++) w[k] = g->R0[3 * k] * wc[0] + g->R0[3 * k + 1] * wc[1] + g->R0[3 * k + 2] * wc[2];
	}
	return QAWS_STATUS_OK;
}

qaws_status qaws_motion_evaluate(qaws_motion_state* out, qaws_motion const* m, qaws_scalar t)
{
	double R[9], T[3], dT[3], w[3];
	unsigned int k;
	qaws_status st;
	if (!out || !m || m->kind == MO_VERTICES || !isfinite((double)t))
		return QAWS_STATUS_INVALID_ARGUMENT;
	if ((st = mo_rigid(m, (double)t, R, T, dT, w)) != QAWS_STATUS_OK)
		return st;
	for (k = 0; k < 9; k++) out->rotation[k] = (qaws_scalar)R[k];
	out->position.x = (qaws_scalar)T[0]; out->position.y = (qaws_scalar)T[1]; out->position.z = (qaws_scalar)T[2];
	out->velocity.x = (qaws_scalar)dT[0]; out->velocity.y = (qaws_scalar)dT[1]; out->velocity.z = (qaws_scalar)dT[2];
	out->angular_velocity.x = (qaws_scalar)w[0]; out->angular_velocity.y = (qaws_scalar)w[1]; out->angular_velocity.z = (qaws_scalar)w[2];
	return QAWS_STATUS_OK;
}

qaws_status qaws_motion_transform_points(qaws_scalar* out, qaws_scalar* out_v, qaws_motion const* m, qaws_scalar t, qaws_scalar const* points,
	unsigned int count)
{
	unsigned int i, k;
	if (!out || !m || !isfinite((double)t))
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (m->kind == MO_VERTICES)
	{
		if (count != m->n)
			return QAWS_STATUS_INVALID_ARGUMENT;
		for (i = 0; i < count; i++)
		{
			double x[3], dx[3];
			if (!m->V)
				for (k = 0; k < 3; k++)
				{
					double a = qaws_internal_array_get(&m->P0, i, k), b = qaws_internal_array_get(&m->P1, i, k);
					dx[k] = b - a;
					x[k] = a + (double)t * dx[k];
				}
			else
				mo_hermite((double const(*)[3])(m->V + 18 * (size_t)i), (double)t, m->quintic, x, dx);
			for (k = 0; k < 3; k++)
			{
				out[3 * (size_t)i + k] = (qaws_scalar)x[k];
				if (out_v)
					out_v[3 * (size_t)i + k] = (qaws_scalar)dx[k];
			}
		}
		return QAWS_STATUS_OK;
	}
	if (count && !points)
		return QAWS_STATUS_INVALID_ARGUMENT;
	{
		double R[9], T[3], dT[3], w[3];
		qaws_status st = mo_rigid(m, (double)t, R, T, dT, w);
		if (st != QAWS_STATUS_OK)
			return st;
		for (i = 0; i < count; i++)
		{
			double p[3], r[3], wr[3];
			for (k = 0; k < 3; k++) p[k] = (double)points[3 * (size_t)i + k];
			for (k = 0; k < 3; k++) r[k] = R[3 * k] * p[0] + R[3 * k + 1] * p[1] + R[3 * k + 2] * p[2];
			mo_cross(w, r, wr);
			for (k = 0; k < 3; k++)
			{
				out[3 * (size_t)i + k] = (qaws_scalar)(r[k] + T[k]);
				if (out_v)
					out_v[3 * (size_t)i + k] = (qaws_scalar)(wr[k] + dT[k]);
			}
		}
	}
	return QAWS_STATUS_OK;
}
