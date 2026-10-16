/*
 * Native control-point optimization (qaws_repulsive.h, R5): the working
 * polyline of qaws_repulsive_polyline_from_curves, its samples (curve,
 * parameter) fixed, is a linear function of the curves' control points,
 * x_s = sum_k B_sk P_k, the same weights for every coordinate. B is read
 * from the curves' own tangent evaluation (qaws_diff.h), one control at a
 * time, and checked: one coordinate's tangent moves that coordinate only,
 * and B P reproduces the samples. The energy's derivatives in P are B^T
 * those in x (exact: the map is linear); the descent runs on the core of
 * qaws_repulsive_solve.c with the metric B^T A B per coordinate.
 */

#include "qaws_repulsive.h"
#include "qaws_curve.h"
#include "qaws_diff.h"
#include "internal/qaws_internal_repulsive.h"
#include "internal/qaws_internal_repulsive_solve.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef qaws_internal_rs_core rs_core;

struct qaws_repulsive_controls
{
	unsigned int dim, ncurve, nc, nv, ne;
	qaws_curve const** curves;
	unsigned int* first;            /* each curve's first control, ncurve + 1 */
	unsigned int *edges, *src;
	qaws_scalar* t;
	unsigned int *brow, *bcol;      /* B by rows (vertices), CSR */
	double* bval;
	double* P0;                     /* the controls read at creation, dim each */
	qaws_repulsive_curves* rc;      /* the energy on the polyline */
};

void qaws_repulsive_controls_destroy(qaws_repulsive_controls* ctl)
{
	if (!ctl)
		return;
	qaws_repulsive_curves_destroy(ctl->rc);
	free(ctl->curves); free(ctl->first); free(ctl->edges); free(ctl->src); free(ctl->t);
	free(ctl->brow); free(ctl->bcol); free(ctl->bval); free(ctl->P0);
	free(ctl);
}

/* x = B P (dim each), in double */
static void rc_apply(qaws_repulsive_controls const* ctl, double const* P, double* x)
{
	unsigned int s, a, q, dim = ctl->dim;
	for (s = 0; s < ctl->nv; s++)
	{
		for (q = 0; q < dim; q++) x[dim * s + q] = 0;
		for (a = ctl->brow[s]; a < ctl->brow[s + 1]; a++)
			for (q = 0; q < dim; q++) x[dim * s + q] += ctl->bval[a] * P[dim * ctl->bcol[a] + q];
	}
}

/* P += seed B^T g */
static void rc_apply_t(qaws_repulsive_controls const* ctl, double const* g, double seed, double* P)
{
	unsigned int s, a, q, dim = ctl->dim;
	for (s = 0; s < ctl->nv; s++)
		for (a = ctl->brow[s]; a < ctl->brow[s + 1]; a++)
			for (q = 0; q < dim; q++) P[dim * ctl->bcol[a] + q] += seed * ctl->bval[a] * g[dim * s + q];
}

/* curve ci's B columns from its tangents, one control at a time; the
   triplets appended to (row, col, val) */
static qaws_status rc_columns(qaws_repulsive_controls* ctl, unsigned int ci, unsigned int ncc, unsigned int** tr, unsigned int** tc, double** tv,
	unsigned int* nt, unsigned int* cap)
{
	unsigned int dim = ctl->dim, ns = 0, s, k, i;
	qaws_scalar *ts = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (ctl->nv + 1)), *dir = (qaws_scalar*)calloc(dim * ncc + 1, sizeof(qaws_scalar));
	unsigned int* rows = (unsigned int*)malloc(sizeof(unsigned int) * (ctl->nv + 1));
	qaws_curve_jet_3d* j3 = (qaws_curve_jet_3d*)malloc(sizeof(qaws_curve_jet_3d) * (ctl->nv + 1));
	qaws_curve_jet_2d* j2 = (qaws_curve_jet_2d*)malloc(sizeof(qaws_curve_jet_2d) * (ctl->nv + 1));
	qaws_field_view view;
	qaws_diff_views views;
	qaws_status st = QAWS_STATUS_OK;
	if (!ts || !dir || !rows || !j3 || !j2)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	for (i = 0; i < ctl->nv; i++)
		if (ctl->src[i] == ci) { rows[ns] = i; ts[ns++] = ctl->t[i]; }
	view = qaws_field_view_make(dir, ncc, QAWS_FIELD_CONTROL_POINTS, dim);
	memset(&views, 0, sizeof(views));
	views.fields = &view; views.field_count = 1;
	for (k = 0; k < ncc && st == QAWS_STATUS_OK; k++)
	{
		dir[dim * k] = 1;
		st = dim == 3 ? qaws_curve_eval_batch_tangent_3d(NULL, j3, ctl->curves[ci], ts, NULL, &views, ns, 1, NULL)
			: qaws_curve_eval_batch_tangent_2d(NULL, j2, ctl->curves[ci], ts, NULL, &views, ns, 1, NULL);
		dir[dim * k] = 0;
		for (s = 0; st == QAWS_STATUS_OK && s < ns; s++)
		{
			double w = dim == 3 ? (double)j3[s].d[0].x : (double)j2[s].d[0].x;
			double other = dim == 3 ? fabs((double)j3[s].d[0].y) + fabs((double)j3[s].d[0].z) : fabs((double)j2[s].d[0].y);
			if (other > 1e-9 * (fabs(w) + 1))
			{
				st = QAWS_STATUS_UNSUPPORTED_OPERATION;   /* (the coordinates are not moved alike) */
				break;
			}
			if (w == 0)
				continue;
			if (*nt == *cap)
			{
				unsigned int nc2 = *cap ? 2 * *cap : 1024;
				unsigned int *r2 = (unsigned int*)realloc(*tr, sizeof(unsigned int) * nc2), *c2;
				double* v2;
				if (r2) *tr = r2;
				c2 = (unsigned int*)realloc(*tc, sizeof(unsigned int) * nc2);
				if (c2) *tc = c2;
				v2 = (double*)realloc(*tv, sizeof(double) * nc2);
				if (v2) *tv = v2;
				if (!r2 || !c2 || !v2)
				{
					st = QAWS_STATUS_ALLOCATION_FAILURE;
					break;
				}
				*cap = nc2;
			}
			(*tr)[*nt] = rows[s]; (*tc)[*nt] = ctl->first[ci] + k; (*tv)[*nt] = w; (*nt)++;
		}
	}
done:
	free(ts); free(dir); free(rows); free(j3); free(j2);
	return st;
}

qaws_status qaws_repulsive_controls_create(qaws_repulsive_controls** out, qaws_repulsive_controls_desc const* desc)
{
	qaws_repulsive_controls* ctl;
	qaws_repulsive_curves_desc cd;
	unsigned int i, ci, dim, *tr = NULL, *tc = NULL, nt = 0, cap = 0;
	double* tv = NULL;
	qaws_scalar* pos = NULL;
	qaws_status st;
	if (!desc || !out || !desc->polyline.curves || !desc->polyline.curve_count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	dim = desc->polyline.dimension;
	ctl = (qaws_repulsive_controls*)calloc(1, sizeof(*ctl));
	if (!ctl)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	ctl->dim = dim; ctl->ncurve = desc->polyline.curve_count;
	/* the working polyline, sized then made */
	st = qaws_repulsive_polyline_from_curves(NULL, NULL, NULL, NULL, &ctl->nv, &ctl->ne, 0, 0, &desc->polyline);
	if (st != QAWS_STATUS_OK && st != QAWS_STATUS_BUFFER_TOO_SMALL)
		goto fail;
	pos = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (dim * ctl->nv + 1));
	ctl->edges = (unsigned int*)malloc(sizeof(unsigned int) * (2 * ctl->ne + 1));
	ctl->src = (unsigned int*)malloc(sizeof(unsigned int) * (ctl->nv + 1));
	ctl->t = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (ctl->nv + 1));
	ctl->curves = (qaws_curve const**)malloc(sizeof(qaws_curve const*) * ctl->ncurve);
	ctl->first = (unsigned int*)malloc(sizeof(unsigned int) * (ctl->ncurve + 1));
	if (!pos || !ctl->edges || !ctl->src || !ctl->t || !ctl->curves || !ctl->first)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto fail;
	}
	if ((st = qaws_repulsive_polyline_from_curves(pos, ctl->edges, ctl->src, ctl->t, &ctl->nv, &ctl->ne, ctl->nv, ctl->ne, &desc->polyline)) != QAWS_STATUS_OK)
		goto fail;
	/* the controls: each curve's CONTROL_POINTS field */
	ctl->first[0] = 0;
	for (ci = 0; ci < ctl->ncurve; ci++)
	{
		qaws_field_desc f[16];
		unsigned int nf = 0, ncc = 0;
		ctl->curves[ci] = desc->polyline.curves[ci];
		if (qaws_curve_describe_fields(f, &nf, ctl->curves[ci], 16) != QAWS_STATUS_OK)
		{
			st = QAWS_STATUS_UNSUPPORTED_OPERATION;
			goto fail;
		}
		for (i = 0; i < nf; i++) if (f[i].field == QAWS_FIELD_CONTROL_POINTS) ncc = f[i].count;
		if (!ncc)
		{
			st = QAWS_STATUS_UNSUPPORTED_OPERATION;   /* (no control points) */
			goto fail;
		}
		ctl->first[ci + 1] = ctl->first[ci] + ncc;
	}
	ctl->nc = ctl->first[ctl->ncurve];
	ctl->P0 = (double*)malloc(sizeof(double) * (dim * ctl->nc + 1));
	if (!ctl->P0)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto fail;
	}
	for (ci = 0; ci < ctl->ncurve; ci++)
	{
		unsigned int ncc = ctl->first[ci + 1] - ctl->first[ci], got = 0;
		qaws_scalar* buf = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (dim * ncc + 1));
		if (!buf)
		{
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			goto fail;
		}
		st = qaws_curve_read_field(buf, &got, ctl->curves[ci], dim * ncc, QAWS_FIELD_CONTROL_POINTS);
		for (i = 0; st == QAWS_STATUS_OK && i < dim * ncc; i++) ctl->P0[dim * ctl->first[ci] + i] = (double)buf[i];
		free(buf);
		if (st != QAWS_STATUS_OK || got != dim * ncc)
		{
			st = st != QAWS_STATUS_OK ? st : QAWS_STATUS_UNSUPPORTED_OPERATION;
			goto fail;
		}
		if ((st = rc_columns(ctl, ci, ncc, &tr, &tc, &tv, &nt, &cap)) != QAWS_STATUS_OK)
			goto fail;
	}
	/* B by rows */
	ctl->brow = (unsigned int*)calloc(ctl->nv + 2, sizeof(unsigned int));
	ctl->bcol = (unsigned int*)malloc(sizeof(unsigned int) * (nt + 1));
	ctl->bval = (double*)malloc(sizeof(double) * (nt + 1));
	if (!ctl->brow || !ctl->bcol || !ctl->bval)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto fail;
	}
	for (i = 0; i < nt; i++) ctl->brow[tr[i] + 1]++;
	for (i = 0; i < ctl->nv; i++) ctl->brow[i + 1] += ctl->brow[i];
	{
		unsigned int* fill = (unsigned int*)calloc(ctl->nv + 1, sizeof(unsigned int));
		if (!fill)
		{
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			goto fail;
		}
		for (i = 0; i < nt; i++)
		{
			unsigned int at = ctl->brow[tr[i]] + fill[tr[i]]++;
			ctl->bcol[at] = tc[i]; ctl->bval[at] = tv[i];
		}
		free(fill);
	}
	/* linear in the controls: B P reproduces the samples */
	{
		double* x = (double*)malloc(sizeof(double) * (dim * ctl->nv + 1));
		double err = 0, ext = 0;
		if (!x)
		{
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			goto fail;
		}
		rc_apply(ctl, ctl->P0, x);
		for (i = 0; i < dim * ctl->nv; i++) { err = fmax(err, fabs(x[i] - (double)pos[i])); ext = fmax(ext, fabs((double)pos[i])); }
		free(x);
		if (!(err <= (sizeof(qaws_scalar) == 4 ? 1e-4 : 1e-9) * (ext + 1)))
		{
			st = QAWS_STATUS_UNSUPPORTED_OPERATION;   /* (not linear in its control points) */
			goto fail;
		}
	}
	memset(&cd, 0, sizeof(cd));
	cd.dimension = dim; cd.vertex_count = ctl->nv; cd.edges = ctl->edges; cd.edge_count = ctl->ne;
	cd.alpha = desc->alpha; cd.beta = desc->beta;
	if ((st = qaws_repulsive_curves_create(&ctl->rc, &cd)) != QAWS_STATUS_OK)
		goto fail;
	free(pos); free(tr); free(tc); free(tv);
	*out = ctl;
	return QAWS_STATUS_OK;
fail:
	free(pos); free(tr); free(tc); free(tv);
	qaws_repulsive_controls_destroy(ctl);
	return st;
}

void qaws_repulsive_controls_get_counts(unsigned int* out_controls, unsigned int* out_vertices, unsigned int* out_edges, qaws_repulsive_controls const* ctl)
{
	if (out_controls) *out_controls = ctl ? ctl->nc : 0;
	if (out_vertices) *out_vertices = ctl ? ctl->nv : 0;
	if (out_edges) *out_edges = ctl ? ctl->ne : 0;
}

qaws_repulsive_curves const* qaws_repulsive_controls_get_curves(qaws_repulsive_controls const* ctl)
{
	return ctl ? ctl->rc : NULL;
}

qaws_status qaws_repulsive_controls_read(qaws_scalar* out_controls, qaws_repulsive_controls const* ctl)
{
	unsigned int i;
	if (!ctl || !out_controls)
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < ctl->dim * ctl->nc; i++) out_controls[i] = (qaws_scalar)ctl->P0[i];
	return QAWS_STATUS_OK;
}

/* x = B P with P from qaws_scalar */
static qaws_status rc_positions(qaws_repulsive_controls const* ctl, qaws_scalar const* controls, qaws_scalar const* dcontrols, qaws_scalar* x,
	qaws_scalar* dx)
{
	unsigned int i, n = ctl->dim * ctl->nc, m = ctl->dim * ctl->nv;
	double *P = (double*)malloc(sizeof(double) * (n + 1)), *y = (double*)malloc(sizeof(double) * (m + 1));
	if (!P || !y)
	{
		free(P); free(y);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < n; i++) P[i] = (double)controls[i];
	rc_apply(ctl, P, y);
	for (i = 0; i < m; i++) x[i] = (qaws_scalar)y[i];
	if (dx)
	{
		for (i = 0; i < n; i++) P[i] = (double)dcontrols[i];
		rc_apply(ctl, P, y);
		for (i = 0; i < m; i++) dx[i] = (qaws_scalar)y[i];
	}
	free(P); free(y);
	return QAWS_STATUS_OK;
}

qaws_status qaws_repulsive_controls_positions(qaws_scalar* out_positions, qaws_repulsive_controls const* ctl, qaws_scalar const* controls)
{
	return ctl && controls && out_positions ? rc_positions(ctl, controls, NULL, out_positions, NULL) : QAWS_STATUS_INVALID_ARGUMENT;
}

qaws_status qaws_repulsive_controls_eval(qaws_scalar* out_value, qaws_scalar* out_tangent, qaws_scalar* out_tangent2,
	qaws_repulsive_controls const* ctl, qaws_scalar const* controls, qaws_scalar const* direction)
{
	qaws_scalar *x, *dx;
	qaws_status st;
	if (!ctl || !controls)
		return QAWS_STATUS_INVALID_ARGUMENT;
	x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (ctl->dim * ctl->nv + 1));
	dx = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (ctl->dim * ctl->nv + 1));
	st = !x || !dx ? QAWS_STATUS_ALLOCATION_FAILURE : rc_positions(ctl, controls, direction, x, direction ? dx : NULL);
	if (st == QAWS_STATUS_OK)
		st = qaws_repulsive_curves_eval(out_value, out_tangent, out_tangent2, ctl->rc, x, direction ? dx : NULL);
	free(x); free(dx);
	return st;
}

qaws_status qaws_repulsive_controls_gradient(qaws_scalar* gradient, qaws_scalar* out_value, qaws_repulsive_controls const* ctl, qaws_scalar const* controls,
	qaws_scalar seed)
{
	unsigned int i, m;
	qaws_scalar *x, *g;
	double *gd, *acc;
	qaws_status st;
	if (!ctl || !controls || !gradient)
		return QAWS_STATUS_INVALID_ARGUMENT;
	m = ctl->dim * ctl->nv;
	x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (m + 1));
	g = (qaws_scalar*)calloc(m + 1, sizeof(qaws_scalar));
	gd = (double*)malloc(sizeof(double) * (m + 1));
	acc = (double*)calloc(ctl->dim * ctl->nc + 1, sizeof(double));
	st = !x || !g || !gd || !acc ? QAWS_STATUS_ALLOCATION_FAILURE : rc_positions(ctl, controls, NULL, x, NULL);
	if (st == QAWS_STATUS_OK)
		st = qaws_repulsive_curves_gradient(g, out_value, ctl->rc, x, 1);
	if (st == QAWS_STATUS_OK)
	{
		for (i = 0; i < m; i++) gd[i] = (double)g[i];
		rc_apply_t(ctl, gd, (double)seed, acc);
		for (i = 0; i < ctl->dim * ctl->nc; i++) gradient[i] = (qaws_scalar)((double)gradient[i] + acc[i]);
	}
	free(x); free(g); free(gd); free(acc);
	return st;
}

qaws_status qaws_repulsive_controls_hvp(qaws_scalar* out_hv, qaws_repulsive_controls const* ctl, qaws_scalar const* controls, qaws_scalar const* direction)
{
	unsigned int i, m;
	qaws_scalar *x, *dx, *h;
	double *hd, *acc;
	qaws_status st;
	if (!ctl || !controls || !direction || !out_hv)
		return QAWS_STATUS_INVALID_ARGUMENT;
	m = ctl->dim * ctl->nv;
	x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (m + 1));
	dx = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (m + 1));
	h = (qaws_scalar*)calloc(m + 1, sizeof(qaws_scalar));
	hd = (double*)malloc(sizeof(double) * (m + 1));
	acc = (double*)calloc(ctl->dim * ctl->nc + 1, sizeof(double));
	st = !x || !dx || !h || !hd || !acc ? QAWS_STATUS_ALLOCATION_FAILURE : rc_positions(ctl, controls, direction, x, dx);
	if (st == QAWS_STATUS_OK)
		st = qaws_repulsive_curves_hvp(h, ctl->rc, x, dx);   /* (B^T H B v: the map is linear) */
	if (st == QAWS_STATUS_OK)
	{
		for (i = 0; i < m; i++) hd[i] = (double)h[i];
		rc_apply_t(ctl, hd, 1, acc);
		for (i = 0; i < ctl->dim * ctl->nc; i++) out_hv[i] = (qaws_scalar)((double)out_hv[i] + acc[i]);
	}
	free(x); free(dx); free(h); free(hd); free(acc);
	return st;
}

qaws_status qaws_repulsive_controls_to_curve(qaws_curve** out_curve, qaws_repulsive_controls const* ctl, qaws_scalar const* controls, unsigned int curve)
{
	qaws_field_view view;
	qaws_diff_views views;
	unsigned int ncc;
	if (!ctl || !controls || !out_curve || curve >= ctl->ncurve)
		return QAWS_STATUS_INVALID_ARGUMENT;
	ncc = ctl->first[curve + 1] - ctl->first[curve];
	view = qaws_field_view_make((qaws_scalar*)(controls + ctl->dim * ctl->first[curve]), ncc, QAWS_FIELD_CONTROL_POINTS, ctl->dim);
	memset(&views, 0, sizeof(views));
	views.fields = &view; views.field_count = 1;
	return qaws_curve_clone_with_fields(out_curve, ctl->curves[curve], &views);
}

/* ------------------------------------------------------------------------ */
/*  The solver                                                              */
/* ------------------------------------------------------------------------ */

struct qaws_repulsive_controls_solver
{
	rs_core core;                   /* (first) x: the controls, dim each */
	qaws_repulsive_controls const* ctl;
	qaws_repulsive_curves_solver* inner;   /* the polyline's rows (pins and barycenter aside) */
	unsigned int *pin, npin, kin;
	unsigned char* fixed;
	double *pin_target, *x, *x1, *Cx, *A, *gx, *col, *bsum, bary0[3];
	double theta;                   /* > 0: Barnes-Hut energy and metric columns */
	int bary;
};

static qaws_status rcs_energy(rs_core const* core, double const* P, double* E, double* grad)
{
	qaws_repulsive_controls_solver const* s = (qaws_repulsive_controls_solver const*)core;
	qaws_repulsive_controls const* ctl = s->ctl;
	qaws_status st;
	rc_apply(ctl, P, s->x);
	if (grad) memset(s->gx, 0, sizeof(double) * ctl->dim * ctl->nv);
	st = s->theta > 0 ? qaws_internal_repulsive_curves_bh(ctl->rc, s->x, s->theta, E, grad ? s->gx : NULL)
		: qaws_internal_repulsive_curves_energy(ctl->rc, s->x, E, grad ? s->gx : NULL);
	if (st == QAWS_STATUS_OK && grad)
		rc_apply_t(ctl, s->gx, 1, grad);
	return st;
}

/* B^T A B (nc x nc): from the dense A, or column by column from the
   hierarchical one (nothing of size |V|^2) */
static qaws_status rcs_metric(rs_core const* core, double const* P, double* M)
{
	qaws_repulsive_controls_solver const* s = (qaws_repulsive_controls_solver const*)core;
	qaws_repulsive_controls const* ctl = s->ctl;
	unsigned int nv = ctl->nv, nc = ctl->nc, k, i, a;
	qaws_internal_rbh_op* op = NULL;
	qaws_status st;
	rc_apply(ctl, P, s->x);
	if (s->theta > 0)
		st = qaws_internal_rbh_curves_metric(ctl->rc, s->x, s->theta, &op);
	else
		st = qaws_internal_repulsive_curves_metric(ctl->rc, s->x, s->A);
	if (st != QAWS_STATUS_OK)
		return st;
	for (k = 0; k < nc; k++)
	{
		/* u = B e_k, y = A u, column k = B^T y */
		double* u = s->col;
		double* y = s->col + nv;
		memset(u, 0, sizeof(double) * nv);
		for (i = 0; i < nv; i++)
			for (a = ctl->brow[i]; a < ctl->brow[i + 1]; a++)
				if (ctl->bcol[a] == k) u[i] = ctl->bval[a];
		if (op)
			qaws_internal_rbh_apply(op, u, y);
		else
			for (i = 0; i < nv; i++)
			{
				double sum = 0;
				unsigned int j;
				for (j = 0; j < nv; j++) sum += s->A[(size_t)i * nv + j] * u[j];
				y[i] = sum;
			}
		for (i = 0; i < nc; i++) M[(size_t)i * nc + k] = 0;
		for (i = 0; i < nv; i++)
			for (a = ctl->brow[i]; a < ctl->brow[i + 1]; a++) M[(size_t)ctl->bcol[a] * nc + k] += ctl->bval[a] * y[i];
	}
	qaws_internal_rbh_destroy(op);
	return QAWS_STATUS_OK;
}

/* rows: the barycenter of the samples, the pins (controls), then the
   polyline's rows times B */
static void rcs_constraints(rs_core const* core, double const* P, double* phi, double* C)
{
	qaws_repulsive_controls_solver const* s = (qaws_repulsive_controls_solver const*)core;
	qaws_repulsive_controls const* ctl = s->ctl;
	unsigned int dim = ctl->dim, N = core->N, nv = ctl->nv, r = 0, i, q, k, a, row;
	rs_core const* in = (rs_core const*)s->inner;
	if (C) memset(C, 0, sizeof(double) * core->k * N);
	rc_apply(ctl, P, s->x);
	if (s->bary)
		for (q = 0; q < dim; q++, r++)
		{
			double m = 0;
			for (i = 0; i < nv; i++) m += s->x[dim * i + q] / nv;
			phi[r] = m - s->bary0[q];
			for (k = 0; C && k < ctl->nc; k++) C[(size_t)r * N + dim * k + q] = s->bsum[k] / nv;
		}
	for (i = 0; i < s->npin; i++)
		for (q = 0; q < dim; q++, r++)
		{
			phi[r] = P[dim * s->pin[i] + q] - s->pin_target[dim * i + q];
			if (C) C[(size_t)r * N + dim * s->pin[i] + q] = 1;
		}
	if (!s->kin)
		return;
	in->ops->constraints(in, s->x, phi + r, C ? s->Cx : NULL);
	if (C)
		for (row = 0; row < s->kin; row++)
			for (i = 0; i < nv; i++)
				for (a = ctl->brow[i]; a < ctl->brow[i + 1]; a++)
					for (q = 0; q < dim; q++)
						C[(size_t)(r + row) * N + dim * ctl->bcol[a] + q] += ctl->bval[a] * s->Cx[(size_t)row * dim * nv + dim * i + q];
}

/* the samples move straight when the controls do */
static int rcs_path_ok(rs_core const* core, double const* P0, double const* P1)
{
	qaws_repulsive_controls_solver const* s = (qaws_repulsive_controls_solver const*)core;
	rc_apply(s->ctl, P0, s->x);
	rc_apply(s->ctl, P1, s->x1);
	return qaws_internal_repulsive_curves_path_clear(s->ctl->rc, core->clearance, s->x, s->x1);
}

static qaws_internal_rs_ops const rcs_ops = { rcs_energy, rcs_metric, rcs_constraints, rcs_path_ok, NULL, NULL, NULL, NULL };

void qaws_repulsive_controls_solver_destroy(qaws_repulsive_controls_solver* s)
{
	if (!s)
		return;
	qaws_repulsive_curves_solver_destroy(s->inner);
	free(s->core.x); free(s->pin); free(s->fixed); free(s->pin_target); free(s->x); free(s->x1); free(s->Cx); free(s->A); free(s->gx); free(s->col);
	free(s->bsum);
	free(s);
}

qaws_status qaws_repulsive_controls_solver_create(qaws_repulsive_controls_solver** out_solver, qaws_repulsive_controls const* ctl,
	qaws_scalar const* controls, qaws_repulsive_curves_solve_desc const* desc)
{
	qaws_repulsive_controls_solver* s;
	qaws_repulsive_curves_solve_desc d, di;
	rs_core* k;
	rs_core const* in;
	qaws_scalar* xs = NULL;
	unsigned int i, q, a, dim, nv, nc;
	double E;
	qaws_status st;
	if (!ctl || !controls || !out_solver)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_solver = NULL;
	memset(&d, 0, sizeof(d));
	if (desc) d = *desc;
	dim = ctl->dim; nv = ctl->nv; nc = ctl->nc;
	if ((d.pinned_count && !d.pinned) || !(d.theta >= 0))
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < d.pinned_count; i++) if (d.pinned[i] >= nc) return QAWS_STATUS_INVALID_ARGUMENT;   /* (pins name controls) */
	s = (qaws_repulsive_controls_solver*)calloc(1, sizeof(*s));
	if (!s)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	s->ctl = ctl; s->theta = (double)d.theta;
	k = &s->core;
	k->ops = &rcs_ops; k->dim = dim; k->nv = nc; k->N = dim * nc;
	k->x = (double*)malloc(sizeof(double) * (dim * nc + 1));
	s->pin = (unsigned int*)malloc(sizeof(unsigned int) * (d.pinned_count + 1));
	s->fixed = (unsigned char*)calloc(nc + 1, 1);
	s->pin_target = (double*)malloc(sizeof(double) * (dim * d.pinned_count + 1));
	s->x = (double*)malloc(sizeof(double) * (dim * nv + 1)); s->x1 = (double*)malloc(sizeof(double) * (dim * nv + 1));
	s->A = s->theta > 0 ? NULL : (double*)malloc(sizeof(double) * ((size_t)nv * nv + 1));
	s->gx = (double*)malloc(sizeof(double) * (dim * nv + 1));
	s->col = (double*)malloc(sizeof(double) * (2 * nv + 1));
	s->bsum = (double*)calloc(nc + 1, sizeof(double));
	xs = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (dim * nv + 1));
	if (!k->x || !s->pin || !s->fixed || !s->pin_target || !s->x || !s->x1 || (!s->A && s->theta <= 0) || !s->gx || !s->col || !s->bsum || !xs)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto fail;
	}
	for (i = 0; i < dim * nc; i++) k->x[i] = (double)controls[i];
	rc_apply(ctl, k->x, s->x);
	for (i = 0; i < nv; i++)
		for (a = ctl->brow[i]; a < ctl->brow[i + 1]; a++) s->bsum[ctl->bcol[a]] += ctl->bval[a];
	/* the polyline's rows (targets at the start) from a curve solver on the samples */
	for (i = 0; i < dim * nv; i++) xs[i] = (qaws_scalar)s->x[i];
	di = d;
	di.pinned = NULL; di.pinned_count = 0; di.hold_barycenter = -1; di.theta = 0;
	if ((st = qaws_repulsive_curves_solver_create(&s->inner, ctl->rc, xs, &di)) != QAWS_STATUS_OK)
		goto fail;
	in = (rs_core const*)s->inner;
	s->kin = in->k;
	s->Cx = (double*)malloc(sizeof(double) * ((size_t)s->kin * dim * nv + 1));
	if (!s->Cx)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto fail;
	}
	s->npin = d.pinned_count;
	for (i = 0; i < s->npin; i++)
	{
		s->pin[i] = d.pinned[i];
		s->fixed[d.pinned[i]] = 1;
		for (q = 0; q < dim; q++) s->pin_target[dim * i + q] = k->x[dim * d.pinned[i] + q];
	}
	s->bary = d.hold_barycenter > 0 || (d.hold_barycenter == 0 && d.pinned_count == 0);
	for (q = 0; q < dim; q++)
	{
		s->bary0[q] = 0;
		for (i = 0; i < nv; i++) s->bary0[q] += s->x[dim * i + q] / nv;
	}
	k->fixed = s->fixed;
	k->translates = s->npin == 0 && !s->bary;
	k->k = (s->bary ? dim : 0) + dim * s->npin + s->kin;
	k->checked = in->checked;
	k->lmean = in->lmean; k->clearance = in->clearance; k->tol = in->tol; k->ctol = in->ctol; k->h = in->h;
	k->rep.reason = QAWS_REPULSIVE_RUNNING;
	k->rep.stationarity = 1;
	if ((st = rcs_energy(k, k->x, &E, NULL)) != QAWS_STATUS_OK)
		goto fail;
	k->rep.initial_energy = k->rep.energy = (qaws_scalar)E;
	free(xs);
	*out_solver = s;
	return QAWS_STATUS_OK;
fail:
	free(xs);
	qaws_repulsive_controls_solver_destroy(s);
	return st;
}

qaws_status qaws_repulsive_controls_solver_direction(qaws_scalar* out_direction, qaws_repulsive_controls_solver const* s, qaws_scalar const* controls)
{
	return s ? qaws_internal_rs_direction(&s->core, controls, out_direction) : QAWS_STATUS_INVALID_ARGUMENT;
}

qaws_status qaws_repulsive_controls_solver_step(qaws_repulsive_report* out_report, qaws_scalar* controls, qaws_repulsive_controls_solver* s)
{
	return s ? qaws_internal_rs_step(&s->core, controls, out_report) : QAWS_STATUS_INVALID_ARGUMENT;
}

qaws_status qaws_repulsive_controls_solve(qaws_repulsive_report* out_report, qaws_scalar* controls, qaws_repulsive_controls const* ctl,
	qaws_repulsive_curves_solve_desc const* desc)
{
	qaws_repulsive_controls_solver* s = NULL;
	qaws_repulsive_report rep;
	unsigned int it, max = desc && desc->max_iterations ? desc->max_iterations : 100;
	qaws_status st = qaws_repulsive_controls_solver_create(&s, ctl, controls, desc);
	if (st != QAWS_STATUS_OK)
		return st;
	memset(&rep, 0, sizeof(rep));
	for (it = 0; it < max && st == QAWS_STATUS_OK; it++)
	{
		st = qaws_repulsive_controls_solver_step(&rep, controls, s);
		if (rep.reason != QAWS_REPULSIVE_RUNNING)
			break;
	}
	if (st == QAWS_STATUS_OK && rep.reason == QAWS_REPULSIVE_RUNNING)
		rep.reason = QAWS_REPULSIVE_ITERATION_LIMIT;
	if (out_report) *out_report = rep;
	qaws_repulsive_controls_solver_destroy(s);
	return st;
}
