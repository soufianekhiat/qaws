#include "qaws_shatter.h"
#include "qaws_surface.h"
#include "qaws_surface_trim.h"
#include "qaws_bspline.h"
#include "qaws_curve.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "internal/qaws_internal_memory.h"
#include "internal/qaws_internal_broadphase.h"
#include "internal/qaws_internal_surface.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define SH_NONE 0xFFFFFFFFu

/* ---- the surfaces laid on grids ------------------------------------------------------------------------------------- */

typedef struct sh_surf
{
	qaws_surface const* s;          /* as given */
	qaws_surface const* base;       /* evaluated: a trimmed surface's base */
	double u0, u1, v0, v1;          /* the base's domain */
	unsigned int operand, index;    /* 0: A, 1: B; the index in it */
	unsigned int first_vertex;
} sh_surf;

typedef struct sh_scene
{
	sh_surf* surf;
	unsigned int ns, n;             /* surfaces; grid cells a side */
	double* x;                      /* 3 per vertex */
	double* uv;                     /* 2 per vertex */
	unsigned int* tri;              /* 3 per triangle */
	unsigned int* tsurf;            /* per triangle its surface */
	unsigned int nv, nt;
	double extent, tol;
} sh_scene;

static void sh_scene_free(sh_scene* S)
{
	qaws_internal_free(S->surf);
	qaws_internal_free(S->x);
	qaws_internal_free(S->uv);
	qaws_internal_free(S->tri);
	qaws_internal_free(S->tsurf);
	memset(S, 0, sizeof *S);
}

static int sh_eval(double* p, double* pu, double* pv, qaws_surface const* s, double u, double v)
{
	qaws_surface_eval_result r;
	memset(&r, 0, sizeof r);
	if (qaws_surface_evaluate(&r, s, (qaws_scalar)u, (qaws_scalar)v, pu ? QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV : QAWS_SURFACE_EVAL_POSITION)
		!= QAWS_STATUS_OK)
		return 0;
	p[0] = (double)r.position.x; p[1] = (double)r.position.y; p[2] = (double)r.position.z;
	if (pu)
	{
		pu[0] = (double)r.du.x; pu[1] = (double)r.du.y; pu[2] = (double)r.du.z;
		pv[0] = (double)r.dv.x; pv[1] = (double)r.dv.y; pv[2] = (double)r.dv.z;
	}
	return 1;
}

static qaws_status sh_scene_build(sh_scene* S, qaws_surface_seam_desc const* d)
{
	unsigned int k, i, j, n, t = 0;
	double lo[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL }, hi[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL }, wu, wv;
	memset(S, 0, sizeof *S);
	S->ns = d->a_count + d->b_count;
	S->n = n = d->resolution ? d->resolution : 64u;
	if (n > 4096)
		return QAWS_STATUS_INVALID_ARGUMENT;
	S->nv = S->ns * (n + 1) * (n + 1);
	S->nt = S->ns * 2 * n * n;
	S->surf = (sh_surf*)qaws_internal_zalloc((size_t)S->ns + 1, sizeof(sh_surf));
	S->x = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)S->nv + 1));
	S->uv = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 2, (size_t)S->nv + 1));
	S->tri = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)S->nt + 1));
	S->tsurf = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)S->nt + 1));
	if (!S->surf || !S->x || !S->uv || !S->tri || !S->tsurf)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (k = 0; k < S->ns; k++)
	{
		sh_surf* f = &S->surf[k];
		qaws_range ru, rv;
		f->operand = k < d->a_count ? 0u : 1u;
		f->index = k < d->a_count ? k : k - d->a_count;
		f->s = f->operand ? d->b[f->index] : d->a[f->index];
		f->base = qaws_internal_trim_base(f->s);
		if (!f->base)
			f->base = f->s;
		ru = qaws_surface_get_u_range(f->base);
		rv = qaws_surface_get_v_range(f->base);
		f->u0 = (double)ru.min_value; f->u1 = (double)ru.max_value;
		f->v0 = (double)rv.min_value; f->v1 = (double)rv.max_value;
		f->first_vertex = k * (n + 1) * (n + 1);
		wu = 0.11 * (0.5 + fmod(0.6180339887498949 * (k + 1), 1.0));
		wv = 0.13 * (0.5 + fmod(0.4142135623730950 * (k + 2), 1.0));
		for (j = 0; j <= n; j++)
			for (i = 0; i <= n; i++)
			{
				unsigned int vi = f->first_vertex + j * (n + 1) + i, q;
				double tu = (double)i / n, tv = (double)j / n, u, v;
				/* (each surface's grid warped its own way, its borders kept: two surfaces' grid lines, and a grid row and
				   another surface, do not line up, which would make crossings at vertices) */
				tu += wu * tu * (1 - tu);
				tv += wv * tv * (1 - tv);
				u = f->u0 + (f->u1 - f->u0) * tu;
				v = f->v0 + (f->v1 - f->v0) * tv;
				if (!sh_eval(S->x + 3 * (size_t)vi, NULL, NULL, f->base, u, v))
					return QAWS_STATUS_INVALID_ARGUMENT;
				S->uv[2 * (size_t)vi] = u;
				S->uv[2 * (size_t)vi + 1] = v;
				for (q = 0; q < 3; q++)
				{
					lo[q] = fmin(lo[q], S->x[3 * (size_t)vi + q]);
					hi[q] = fmax(hi[q], S->x[3 * (size_t)vi + q]);
				}
			}
		for (j = 0; j < n; j++)
			for (i = 0; i < n; i++)
			{
				unsigned int a = f->first_vertex + j * (n + 1) + i, b = a + 1, c = a + n + 2, e = a + n + 1;
				S->tri[3 * (size_t)t] = a; S->tri[3 * (size_t)t + 1] = b; S->tri[3 * (size_t)t + 2] = c; S->tsurf[t++] = k;
				S->tri[3 * (size_t)t] = a; S->tri[3 * (size_t)t + 1] = c; S->tri[3 * (size_t)t + 2] = e; S->tsurf[t++] = k;
			}
	}
	S->extent = sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) + (hi[2] - lo[2]) * (hi[2] - lo[2]));
	if (!(S->extent > 0))
		S->extent = 1;
	S->tol = d->tolerance > 0 ? (double)d->tolerance : (QAWS_SCALAR_IS_FLOAT ? 1e-5 : 1e-9) * S->extent;
	return QAWS_STATUS_OK;
}

/* the seam kind of a surface pair */
static unsigned int sh_kind(sh_scene const* S, unsigned int p, unsigned int q)
{
	unsigned int op = S->surf[p].operand, oq = S->surf[q].operand;
	return op != oq ? QAWS_SEAM_A_B : op ? QAWS_SEAM_B_B : QAWS_SEAM_A_A;
}

/* ---- crossings of triangle pairs -------------------------------------------------------------------------------------- */

static double sh_orient(double const* a, double const* b, double const* c, double const* d)
{
	double u[3], v[3], w[3];
	int k;
	for (k = 0; k < 3; k++)
	{
		u[k] = b[k] - a[k];
		v[k] = c[k] - a[k];
		w[k] = d[k] - a[k];
	}
	return u[0] * (v[1] * w[2] - v[2] * w[1]) + u[1] * (v[2] * w[0] - v[0] * w[2]) + u[2] * (v[0] * w[1] - v[1] * w[0]);
}

/* a crossing: edge (e0 < e1, vertices) of one surface's grid through triangle t of another; its point and parameters */
typedef struct sh_event
{
	unsigned long long edge;        /* e0 << 32 | e1 */
	unsigned int tri;
	unsigned int pair;              /* the surface pair's slot (p * ns + q, p < q) */
	double x[3];
	double uv[4];                   /* on the pair's lower surface, then the higher */
	unsigned char border;           /* the edge lies on its grid's border */
} sh_event;

typedef struct sh_seg_job
{
	sh_scene const* S;
	unsigned int kinds;
	sh_event* ev;
	unsigned int nev, cap;
	int oom;
} sh_seg_job;

static int sh_accept(void* user, unsigned int i, unsigned int j)
{
	sh_seg_job const* J = (sh_seg_job const*)user;
	unsigned int p = J->S->tsurf[i], q = J->S->tsurf[j];
	return p != q && (sh_kind(J->S, p, q) & J->kinds);
}

static int sh_on_border(sh_scene const* S, unsigned int a, unsigned int b)
{
	unsigned int n = S->n, k = a / ((n + 1) * (n + 1)), la = a - k * (n + 1) * (n + 1), lb = b - k * (n + 1) * (n + 1);
	unsigned int ia = la % (n + 1), ja = la / (n + 1), ib = lb % (n + 1), jb = lb / (n + 1);
	return (ia == ib && (ia == 0 || ia == n)) || (ja == jb && (ja == 0 || ja == n));
}

/* does edge (a, b) of one grid pierce triangle t of another? (the same answer from either triangle of the edge: the edge's
   vertices in order, a sign of zero read as positive) */
static int sh_pierce(sh_scene const* S, unsigned int a, unsigned int b, unsigned int t, double* x, double* w)
{
	double const *pa = S->x + 3 * (size_t)a, *pb = S->x + 3 * (size_t)b;
	double const *t0 = S->x + 3 * (size_t)S->tri[3 * (size_t)t], *t1 = S->x + 3 * (size_t)S->tri[3 * (size_t)t + 1],
		*t2 = S->x + 3 * (size_t)S->tri[3 * (size_t)t + 2];
	double da = sh_orient(t0, t1, t2, pa), db = sh_orient(t0, t1, t2, pb), w0, w1, w2, s;
	int k;
	if ((da >= 0) == (db >= 0))
		return 0;
	w0 = sh_orient(pa, pb, t1, t2);
	w1 = sh_orient(pa, pb, t2, t0);
	w2 = sh_orient(pa, pb, t0, t1);
	if ((w0 >= 0) != (w1 >= 0) || (w1 >= 0) != (w2 >= 0))
		return 0;
	s = w0 + w1 + w2;
	if (s == 0)
		return 0;
	w[0] = w0 / s; w[1] = w1 / s; w[2] = w2 / s;
	for (k = 0; k < 3; k++)
		x[k] = pa[k] + da / (da - db) * (pb[k] - pa[k]);
	return 1;
}

/* (u, v) on its surface's domain border (to 1e-9 of the domain) */
static int sh_uv_border(sh_surf const* f, double const* uv)
{
	double eu = 1e-9 * (f->u1 - f->u0), ev = 1e-9 * (f->v1 - f->v0);
	return uv[0] <= f->u0 + eu || uv[0] >= f->u1 - eu || uv[1] <= f->v0 + ev || uv[1] >= f->v1 - ev;
}

static int sh_push(sh_seg_job* J, sh_event const* e)
{
	if (J->nev == J->cap)
	{
		unsigned int nc = J->cap ? 2 * J->cap : 1024;
		sh_event* ne = (sh_event*)qaws_internal_realloc(J->ev, qaws_internal_mul2(sizeof(sh_event), nc));
		if (!ne)
			return 0;
		J->ev = ne;
		J->cap = nc;
	}
	J->ev[J->nev++] = *e;
	return 1;
}

/* a triangle pair's crossing: its two events, pushed one after the other (a segment); none when not exactly two */
static qaws_status sh_visit(void* user, unsigned int i, unsigned int j)
{
	sh_seg_job* J = (sh_seg_job*)user;
	sh_scene const* S = J->S;
	sh_event e[6];
	unsigned int ne = 0, side, k, p = S->tsurf[i], q = S->tsurf[j];
	for (side = 0; side < 2; side++)
	{
		unsigned int te = side ? j : i, tt = side ? i : j;
		for (k = 0; k < 3; k++)
		{
			unsigned int a = S->tri[3 * (size_t)te + k], b = S->tri[3 * (size_t)te + (k + 1) % 3], c;
			double x[3], w[3], t2;
			sh_event* E = &e[ne];
			if (a > b)
			{
				c = a; a = b; b = c;
			}
			if (!sh_pierce(S, a, b, tt, x, w))
				continue;
			memcpy(E->x, x, sizeof x);
			E->edge = (unsigned long long)a << 32 | b;
			E->tri = tt;
			E->pair = (p < q ? p : q) * S->ns + (p < q ? q : p);
			E->border = (unsigned char)sh_on_border(S, a, b);
			{
				/* the edge's surface: along the edge; the triangle's: barycentric */
				double const *ua = S->uv + 2 * (size_t)a, *ub = S->uv + 2 * (size_t)b;
				double const *u0 = S->uv + 2 * (size_t)S->tri[3 * (size_t)tt], *u1 = S->uv + 2 * (size_t)S->tri[3 * (size_t)tt + 1],
					*u2 = S->uv + 2 * (size_t)S->tri[3 * (size_t)tt + 2];
				double const* pa = S->x + 3 * (size_t)a, *pb = S->x + 3 * (size_t)b;
				double l2 = 0, d2 = 0, ue[2], ut[2];
				int lower_is_edge = S->tsurf[te] < S->tsurf[tt];
				for (c = 0; c < 3; c++)
				{
					l2 += (pb[c] - pa[c]) * (pb[c] - pa[c]);
					d2 += (x[c] - pa[c]) * (pb[c] - pa[c]);
				}
				t2 = l2 > 0 ? d2 / l2 : 0;
				for (c = 0; c < 2; c++)
				{
					ue[c] = ua[c] + t2 * (ub[c] - ua[c]);
					ut[c] = w[0] * u0[c] + w[1] * u1[c] + w[2] * u2[c];
				}
				memcpy(E->uv, lower_is_edge ? ue : ut, sizeof ue);
				memcpy(E->uv + 2, lower_is_edge ? ut : ue, sizeof ue);
			}
			if (++ne == 6)
				break;
		}
	}
	if (ne != 2)
		return QAWS_STATUS_OK;
	/* (one operand's two surfaces touching along a border of either, as a body's faces do: not a crossing) */
	if (S->surf[p].operand == S->surf[q].operand)
	{
		sh_surf const *lo = &S->surf[p < q ? p : q], *up = &S->surf[p < q ? q : p];
		if ((sh_uv_border(lo, e[0].uv) && sh_uv_border(lo, e[1].uv)) || (sh_uv_border(up, e[0].uv + 2) && sh_uv_border(up, e[1].uv + 2)))
			return QAWS_STATUS_OK;
	}
	if (!sh_push(J, &e[0]) || !sh_push(J, &e[1]))
	{
		J->oom = 1;
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	return QAWS_STATUS_OK;
}

/* ---- Newton onto both surfaces ---------------------------------------------------------------------------------------- */

static double sh_clamp(double x, double lo, double hi)
{
	return x < lo ? lo : x > hi ? hi : x;
}

/* the least change of (u1, v1, u2, v2) making S1 = S2, while the gap shrinks; x: the middle of the two points */
static void sh_newton(sh_surf const* f1, sh_surf const* f2, double* uv, double tol, double* x)
{
	double P[3], Pu[3], Pv[3], Q[3], Qu[3], Qv[3], F[3], gap;
	int it, k;
	if (!sh_eval(P, Pu, Pv, f1->base, uv[0], uv[1]) || !sh_eval(Q, Qu, Qv, f2->base, uv[2], uv[3]))
		return;
	for (k = 0; k < 3; k++)
		F[k] = P[k] - Q[k];
	gap = sqrt(F[0] * F[0] + F[1] * F[1] + F[2] * F[2]);
	for (it = 0; it < 16 && gap > tol; it++)
	{
		/* J = [Pu Pv -Qu -Qv] over the free parameters; delta = -J^T (J J^T)^-1 F. A parameter on its domain's bound that
		   the step would push out is held there and the step taken again over the others (three left: plain Newton) */
		double const* col[4] = { Pu, Pv, Qu, Qv };
		double const lo[4] = { f1->u0, f1->v0, f2->u0, f2->v0 }, hi[4] = { f1->u1, f1->v1, f2->u1, f2->v1 };
		double A[3][3], y[3], det = 0, inv[3][3], step = 1, nu[4], d[4];
		int r, c, h, ok = 0, held[4] = { 0, 0, 0, 0 }, pass, again = 1;
		for (pass = 0; pass < 3 && again; pass++)
		{
			again = 0;
			for (r = 0; r < 3; r++)
				for (c = 0; c < 3; c++)
				{
					A[r][c] = 0;
					for (k = 0; k < 4; k++)
						if (!held[k])
							A[r][c] += col[k][r] * col[k][c];
				}
			inv[0][0] = A[1][1] * A[2][2] - A[1][2] * A[2][1];
			inv[0][1] = A[0][2] * A[2][1] - A[0][1] * A[2][2];
			inv[0][2] = A[0][1] * A[1][2] - A[0][2] * A[1][1];
			inv[1][0] = A[1][2] * A[2][0] - A[1][0] * A[2][2];
			inv[1][1] = A[0][0] * A[2][2] - A[0][2] * A[2][0];
			inv[1][2] = A[0][2] * A[1][0] - A[0][0] * A[1][2];
			inv[2][0] = A[1][0] * A[2][1] - A[1][1] * A[2][0];
			inv[2][1] = A[0][1] * A[2][0] - A[0][0] * A[2][1];
			inv[2][2] = A[0][0] * A[1][1] - A[0][1] * A[1][0];
			det = A[0][0] * inv[0][0] + A[0][1] * inv[1][0] + A[0][2] * inv[2][0];
			if (!(fabs(det) > 0))
				break;
			for (r = 0; r < 3; r++)
				y[r] = (inv[r][0] * F[0] + inv[r][1] * F[1] + inv[r][2] * F[2]) / det;
			for (k = 0; k < 4; k++)
			{
				double span = hi[k] - lo[k];
				d[k] = held[k] ? 0 : -(col[k][0] * y[0] + col[k][1] * y[1] + col[k][2] * y[2]) * (k < 2 ? 1 : -1);
				if (!held[k] && ((uv[k] <= lo[k] + 1e-14 * span && d[k] < 0) || (uv[k] >= hi[k] - 1e-14 * span && d[k] > 0)))
				{
					held[k] = 1;
					again = 1;
				}
			}
		}
		if (!(fabs(det) > 0))
			break;
		for (h = 0; h < 6 && !ok; h++, step *= 0.5)
		{
			double P2[3], Q2[3], g2 = 0;
			for (k = 0; k < 4; k++)
				nu[k] = uv[k] + step * d[k];
			nu[0] = sh_clamp(nu[0], f1->u0, f1->u1); nu[1] = sh_clamp(nu[1], f1->v0, f1->v1);
			nu[2] = sh_clamp(nu[2], f2->u0, f2->u1); nu[3] = sh_clamp(nu[3], f2->v0, f2->v1);
			if (!sh_eval(P2, NULL, NULL, f1->base, nu[0], nu[1]) || !sh_eval(Q2, NULL, NULL, f2->base, nu[2], nu[3]))
				break;
			for (k = 0; k < 3; k++)
				g2 += (P2[k] - Q2[k]) * (P2[k] - Q2[k]);
			if (sqrt(g2) < gap)
				ok = 1;
		}
		if (!ok)
			break;
		memcpy(uv, nu, sizeof nu);
		sh_eval(P, Pu, Pv, f1->base, uv[0], uv[1]);
		sh_eval(Q, Qu, Qv, f2->base, uv[2], uv[3]);
		for (k = 0; k < 3; k++)
			F[k] = P[k] - Q[k];
		gap = sqrt(F[0] * F[0] + F[1] * F[1] + F[2] * F[2]);
	}
	for (k = 0; k < 3; k++)
		x[k] = 0.5 * (P[k] + Q[k]);
}

/* ---- seams ------------------------------------------------------------------------------------------------------------ */

typedef struct sh_seam
{
	unsigned int p, q;              /* surfaces, p < q */
	unsigned int first, count;      /* points */
	int closed;
} sh_seam;

struct qaws_surface_seams
{
	sh_seam* seam;
	unsigned int count;
	qaws_seam_point* pt;
	sh_surf* surf;                  /* (a copy, for the operands and indices) */
	unsigned int ns;
};

static int sh_cmp_key(void const* p, void const* q)
{
	sh_event const *a = (sh_event const*)p, *b = (sh_event const*)q;
	if (a->pair != b->pair)
		return a->pair < b->pair ? -1 : 1;
	if (a->edge != b->edge)
		return a->edge < b->edge ? -1 : 1;
	return a->tri < b->tri ? -1 : a->tri > b->tri;
}

typedef struct sh_ord
{
	sh_event const* e;
	unsigned int i;
} sh_ord;

static int sh_cmp_ord(void const* p, void const* q)
{
	sh_ord const *a = (sh_ord const*)p, *b = (sh_ord const*)q;
	int c = sh_cmp_key(a->e, b->e);
	return c ? c : a->i < b->i ? -1 : a->i > b->i;
}

void qaws_surface_seams_destroy(qaws_surface_seams* seams)
{
	if (!seams)
		return;
	qaws_internal_free(seams->seam);
	qaws_internal_free(seams->pt);
	qaws_internal_free(seams->surf);
	qaws_internal_free(seams);
}

/* the seams of scene S (its kinds), its points moved onto both surfaces and those outside a trim left out */
static qaws_status sh_seams(qaws_surface_seams** out, sh_scene* S, unsigned int kinds)
{
	sh_seg_job J;
	qaws_bp_box* box = NULL;
	sh_ord* ord = NULL;
	unsigned int *id = NULL, *deg = NULL, *start = NULL, *adj = NULL, *nb = NULL, i, k, nid = 0, cap_s = 0, cap_p = 0, np = 0;
	unsigned char* used = NULL;
	qaws_surface_seams* R = NULL;
	double* rx = NULL;
	qaws_status st;
	memset(&J, 0, sizeof J);
	J.S = S;
	J.kinds = kinds;
	box = (qaws_bp_box*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_bp_box), (size_t)S->nt + 1));
	R = (qaws_surface_seams*)qaws_internal_zalloc(1, sizeof *R);
	if (!box || !R)
		goto oom;
	for (i = 0; i < S->nt; i++)
		for (k = 0; k < 3; k++)
		{
			double a = S->x[3 * (size_t)S->tri[3 * (size_t)i] + k], b = S->x[3 * (size_t)S->tri[3 * (size_t)i + 1] + k],
				c = S->x[3 * (size_t)S->tri[3 * (size_t)i + 2] + k], pad = 1e-12 * S->extent;
			box[i].lo[k] = fmin(a, fmin(b, c)) - pad;
			box[i].hi[k] = fmax(a, fmax(b, c)) + pad;
		}
	st = qaws_internal_broadphase(NULL, box, S->nt, 3, sh_accept, sh_visit, &J);
	if (J.oom)
		goto oom;
	if (st != QAWS_STATUS_OK)
		goto done;
	/* one id per event (key: pair, edge, triangle) */
	ord = (sh_ord*)qaws_internal_malloc(qaws_internal_mul2(sizeof(sh_ord), (size_t)J.nev + 1));
	id = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)J.nev + 1));
	if (!ord || !id)
		goto oom;
	for (i = 0; i < J.nev; i++)
	{
		ord[i].e = &J.ev[i];
		ord[i].i = i;
	}
	qsort(ord, J.nev, sizeof *ord, sh_cmp_ord);
	for (i = 0; i < J.nev; i++)
	{
		if (i && sh_cmp_key(ord[i - 1].e, ord[i].e))
			nid++;
		id[ord[i].i] = nid;
	}
	if (J.nev)
		nid++;
	/* the segments (events 2s, 2s + 1) as a graph on the ids; each id's representative event */
	deg = (unsigned int*)qaws_internal_zalloc((size_t)nid + 1, sizeof(unsigned int));
	start = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nid + 2));
	adj = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)J.nev + 1));
	nb = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nid + 1));
	used = (unsigned char*)qaws_internal_zalloc((size_t)J.nev / 2 + 1, 1);
	rx = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 7, (size_t)nid + 1));
	if (!deg || !start || !adj || !nb || !used || !rx)
		goto oom;
	for (i = 0; i < J.nev; i++)
	{
		deg[id[i]]++;
		nb[id[i]] = i;
	}
	for (i = 0, start[0] = 0; i < nid; i++)
		start[i + 1] = start[i] + deg[i];
	memset(deg, 0, sizeof(unsigned int) * nid);
	for (i = 0; i < J.nev; i++)
		adj[start[id[i]] + deg[id[i]]++] = i / 2;
	/* each point moved onto both surfaces (x, then u1 v1 u2 v2) */
	for (i = 0; i < nid; i++)
	{
		sh_event const* e = &J.ev[nb[i]];
		unsigned int p = e->pair / S->ns, q = e->pair % S->ns;
		double uv[4];
		memcpy(uv, e->uv, sizeof uv);
		memcpy(rx + 7 * (size_t)i, e->x, 3 * sizeof(double));
		sh_newton(&S->surf[p], &S->surf[q], uv, S->tol, rx + 7 * (size_t)i);
		memcpy(rx + 7 * (size_t)i + 3, uv, sizeof uv);
	}
	/* chained: from the ends (ids of one segment), then the loops; a point outside a trim breaks the chain */
	{
		unsigned int pass;
		for (pass = 0; pass < 2; pass++)
			for (i = 0; i < J.nev / 2; i++)
			{
				unsigned int a = id[2 * i], b = id[2 * i + 1], cur, seg = i, first_pt = np, closed = 0, broke = 0;
				sh_event const* e0 = &J.ev[2 * i];
				unsigned int p = e0->pair / S->ns, q = e0->pair % S->ns;
				if (used[i] || a == b)
					continue;
				if (pass == 0)
				{
					if (deg[a] == 2 && deg[b] == 2)
						continue;
					cur = deg[a] == 2 ? b : a;
				}
				else
					cur = a;
				for (;;)
				{
					unsigned int nxt, j2, other = SH_NONE;
					int in1 = 1, in2 = 1;
					double const* r = rx + 7 * (size_t)cur;
					/* (the point, unless a trim leaves it out) */
					if (S->surf[p].s != S->surf[p].base)
					{
						qaws_surface_trim_contains(&in1, S->surf[p].s, (qaws_scalar)r[3], (qaws_scalar)r[4]);
					}
					if (S->surf[q].s != S->surf[q].base)
					{
						qaws_surface_trim_contains(&in2, S->surf[q].s, (qaws_scalar)r[5], (qaws_scalar)r[6]);
					}
					if (in1 && in2)
					{
						if (np == cap_p)
						{
							unsigned int nc = cap_p ? 2 * cap_p : 256;
							qaws_seam_point* n2 = (qaws_seam_point*)qaws_internal_realloc(R->pt, qaws_internal_mul2(sizeof(qaws_seam_point), nc));
							if (!n2)
								goto oom;
							R->pt = n2;
							cap_p = nc;
						}
						R->pt[np].position.x = (qaws_scalar)r[0]; R->pt[np].position.y = (qaws_scalar)r[1]; R->pt[np].position.z = (qaws_scalar)r[2];
						R->pt[np].u1 = (qaws_scalar)r[3]; R->pt[np].v1 = (qaws_scalar)r[4];
						R->pt[np].u2 = (qaws_scalar)r[5]; R->pt[np].v2 = (qaws_scalar)r[6];
						np++;
					}
					else if (np > first_pt)
					{
						/* (left out: the seam so far ends here) */
						if (np - first_pt >= 2)
						{
							if (R->count == cap_s)
							{
								unsigned int nc = cap_s ? 2 * cap_s : 64;
								sh_seam* s2 = (sh_seam*)qaws_internal_realloc(R->seam, qaws_internal_mul2(sizeof(sh_seam), nc));
								if (!s2)
									goto oom;
								R->seam = s2;
								cap_s = nc;
							}
							R->seam[R->count].p = p; R->seam[R->count].q = q;
							R->seam[R->count].first = first_pt; R->seam[R->count].count = np - first_pt;
							R->seam[R->count++].closed = 0;
						}
						else
							np = first_pt;
						first_pt = np;
						broke = 1;
					}
					if (seg == SH_NONE)
						break;
					used[seg] = 1;
					nxt = id[2 * seg] == cur ? id[2 * seg + 1] : id[2 * seg];
					cur = nxt;
					seg = SH_NONE;
					if (deg[cur] == 2)
						for (j2 = start[cur]; j2 < start[cur] + 2; j2++)
							if (!used[adj[j2]])
								other = adj[j2];
					seg = other;
					if (seg == SH_NONE && pass == 1)
					{
						closed = !broke && np - first_pt >= 3;
						break;   /* (a loop: back at its first point, not repeated) */
					}
				}
				if (np - first_pt >= 2)
				{
					if (R->count == cap_s)
					{
						unsigned int nc = cap_s ? 2 * cap_s : 64;
						sh_seam* s2 = (sh_seam*)qaws_internal_realloc(R->seam, qaws_internal_mul2(sizeof(sh_seam), nc));
						if (!s2)
							goto oom;
						R->seam = s2;
						cap_s = nc;
					}
					R->seam[R->count].p = p; R->seam[R->count].q = q;
					R->seam[R->count].first = first_pt; R->seam[R->count].count = np - first_pt;
					R->seam[R->count++].closed = (int)closed;
				}
				else
					np = first_pt;
			}
	}
	R->surf = (sh_surf*)qaws_internal_malloc(qaws_internal_mul2(sizeof(sh_surf), (size_t)S->ns + 1));
	if (!R->surf)
		goto oom;
	memcpy(R->surf, S->surf, sizeof(sh_surf) * S->ns);
	R->ns = S->ns;
	*out = R;
	R = NULL;
	st = QAWS_STATUS_OK;
	goto done;
oom:
	st = QAWS_STATUS_ALLOCATION_FAILURE;
done:
	qaws_surface_seams_destroy(R);
	qaws_internal_free(box);
	qaws_internal_free(J.ev);
	qaws_internal_free(ord);
	qaws_internal_free(id);
	qaws_internal_free(deg);
	qaws_internal_free(start);
	qaws_internal_free(adj);
	qaws_internal_free(nb);
	qaws_internal_free(used);
	qaws_internal_free(rx);
	return st;
}

static int sh_desc_ok(qaws_surface_seam_desc const* d)
{
	unsigned int i;
	if (!d || (!d->a_count && !d->b_count) || (d->a_count && !d->a) || (d->b_count && !d->b))
		return 0;
	for (i = 0; i < d->a_count; i++)
		if (!d->a[i])
			return 0;
	for (i = 0; i < d->b_count; i++)
		if (!d->b[i])
			return 0;
	return 1;
}

qaws_status qaws_surface_seams_create(qaws_surface_seams** out_seams, qaws_surface_seam_desc const* desc)
{
	sh_scene S;
	qaws_status st;
	if (!out_seams)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_seams = NULL;
	if (!sh_desc_ok(desc) || !(desc->kinds & (QAWS_SEAM_A_B | QAWS_SEAM_A_A | QAWS_SEAM_B_B)))
		return QAWS_STATUS_INVALID_ARGUMENT;
	st = sh_scene_build(&S, desc);
	if (st == QAWS_STATUS_OK)
		st = sh_seams(out_seams, &S, desc->kinds);
	sh_scene_free(&S);
	return st;
}

unsigned int qaws_surface_seams_get_count(qaws_surface_seams const* seams)
{
	return seams ? seams->count : 0;
}

qaws_status qaws_surface_seams_get(qaws_seam_info* out_info, qaws_seam_point* out_points, qaws_surface_seams const* seams, unsigned int index)
{
	sh_seam const* s;
	if (!seams || index >= seams->count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	s = &seams->seam[index];
	if (out_info)
	{
		sh_surf const *f1 = &seams->surf[s->p], *f2 = &seams->surf[s->q];
		out_info->kind = f1->operand != f2->operand ? QAWS_SEAM_A_B : f1->operand ? QAWS_SEAM_B_B : QAWS_SEAM_A_A;
		out_info->operand1 = f1->operand;
		out_info->surface1 = f1->index;
		out_info->operand2 = f2->operand;
		out_info->surface2 = f2->index;
		out_info->point_count = s->count;
		out_info->closed = s->closed;
	}
	if (out_points)
		memcpy(out_points, seams->pt + s->first, sizeof(qaws_seam_point) * s->count);
	return QAWS_STATUS_OK;
}

/* ---- a planar arrangement of segments: its bounded faces -------------------------------------------------------------- */

typedef struct sh_arr
{
	double* seg;                    /* 4 per segment */
	unsigned int ns, caps;
	/* the result */
	double* vx;                     /* 2 per vertex */
	unsigned int nv;
	unsigned int* loop_first;       /* per cycle its first entry in loop_v, count + 1 */
	unsigned int* loop_v;           /* the cycles' vertices */
	double* loop_area;              /* signed */
	unsigned int nloop;
	unsigned int* hole_of;          /* per cycle: the face it is a hole of (negative cycles), else SH_NONE */
} sh_arr;

static void sh_arr_free(sh_arr* A)
{
	qaws_internal_free(A->seg);
	qaws_internal_free(A->vx);
	qaws_internal_free(A->loop_first);
	qaws_internal_free(A->loop_v);
	qaws_internal_free(A->loop_area);
	qaws_internal_free(A->hole_of);
	memset(A, 0, sizeof *A);
}

static int sh_arr_add(sh_arr* A, double x0, double y0, double x1, double y1)
{
	if (x0 == x1 && y0 == y1)
		return 1;
	if (A->ns == A->caps)
	{
		unsigned int nc = A->caps ? 2 * A->caps : 256;
		double* s2 = (double*)qaws_internal_realloc(A->seg, qaws_internal_mul3(sizeof(double), 4, nc));
		if (!s2)
			return 0;
		A->seg = s2;
		A->caps = nc;
	}
	A->seg[4 * (size_t)A->ns] = x0; A->seg[4 * (size_t)A->ns + 1] = y0;
	A->seg[4 * (size_t)A->ns + 2] = x1; A->seg[4 * (size_t)A->ns + 3] = y1;
	A->ns++;
	return 1;
}

typedef struct sh_cut_job
{
	sh_arr const* A;
	double* cut;                    /* (segment, t) pairs */
	unsigned int nc, cap;
	int oom;
} sh_cut_job;

static int sh_cut_push(sh_cut_job* J, unsigned int s, double t)
{
	if (J->nc == J->cap)
	{
		unsigned int nc = J->cap ? 2 * J->cap : 256;
		double* c2 = (double*)qaws_internal_realloc(J->cut, qaws_internal_mul3(sizeof(double), 2, nc));
		if (!c2)
			return 0;
		J->cut = c2;
		J->cap = nc;
	}
	J->cut[2 * (size_t)J->nc] = (double)s;
	J->cut[2 * (size_t)J->nc + 1] = t;
	J->nc++;
	return 1;
}

static qaws_status sh_cut_visit(void* user, unsigned int i, unsigned int j)
{
	sh_cut_job* J = (sh_cut_job*)user;
	double const *a = J->A->seg + 4 * (size_t)i, *b = J->A->seg + 4 * (size_t)j;
	double rx = a[2] - a[0], ry = a[3] - a[1], sx = b[2] - b[0], sy = b[3] - b[1], den = rx * sy - ry * sx, qx = b[0] - a[0], qy = b[1] - a[1], t, u;
	if (den == 0)
		return QAWS_STATUS_OK;
	t = (qx * sy - qy * sx) / den;
	u = (qx * ry - qy * rx) / den;
	if (t < 0 || t > 1 || u < 0 || u > 1)
		return QAWS_STATUS_OK;
	if (!sh_cut_push(J, i, t) || !sh_cut_push(J, j, u))
	{
		J->oom = 1;
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	return QAWS_STATUS_OK;
}

typedef struct sh_pt
{
	double x, y;
	unsigned int i;
} sh_pt;

static int sh_cmp_pt(void const* p, void const* q)
{
	sh_pt const *a = (sh_pt const*)p, *b = (sh_pt const*)q;
	return a->x < b->x ? -1 : a->x > b->x ? 1 : a->y < b->y ? -1 : a->y > b->y ? 1 : 0;
}

static int sh_cmp_dbl2(void const* p, void const* q)
{
	double const *a = (double const*)p, *b = (double const*)q;
	return a[0] < b[0] ? -1 : a[0] > b[0] ? 1 : a[1] < b[1] ? -1 : a[1] > b[1];
}

static int sh_cmp_u64(void const* p, void const* q)
{
	unsigned long long a = *(unsigned long long const*)p, b = *(unsigned long long const*)q;
	return a < b ? -1 : a > b;
}

static unsigned int sh_find(unsigned int* par, unsigned int i)
{
	while (par[i] != i)
	{
		par[i] = par[par[i]];
		i = par[i];
	}
	return i;
}

typedef struct sh_out_edge
{
	unsigned int from, to;
	double angle;
} sh_out_edge;

static int sh_cmp_out(void const* p, void const* q)
{
	sh_out_edge const *a = (sh_out_edge const*)p, *b = (sh_out_edge const*)q;
	if (a->from != b->from)
		return a->from < b->from ? -1 : 1;
	return a->angle < b->angle ? -1 : a->angle > b->angle;
}

/* winding of polygon (n points, 2 each) around (x, y) */
static int sh_wind(double const* p, unsigned int n, double x, double y)
{
	unsigned int i;
	int w = 0;
	for (i = 0; i < n; i++)
	{
		double x1 = p[2 * i], y1 = p[2 * i + 1], x2 = p[2 * ((i + 1) % n)], y2 = p[2 * ((i + 1) % n) + 1];
		double c = (x2 - x1) * (y - y1) - (x - x1) * (y2 - y1);
		if (y1 <= y)
		{
			if (y2 > y && c > 0)
				w++;
		}
		else if (y2 <= y && c < 0)
			w--;
	}
	return w;
}

/* the segments split where they cross, their ends within eps one vertex, edges with a free end pruned, and the cycles
   traced with their face on the left: positive ones bound faces, each negative one is the hole of the smallest face
   around it (or of none) */
static qaws_status sh_arr_build(sh_arr* A, double eps)
{
	sh_cut_job J;
	qaws_bp_box* box = NULL;
	sh_pt* pt = NULL;
	unsigned int *par = NULL, *pid = NULL, *deg = NULL, *nxt = NULL, *first = NULL, i, j, k, npt, ne = 0, nh;
	unsigned long long* ed = NULL;
	sh_out_edge* he = NULL;
	double* cuts = NULL;
	unsigned char* seen = NULL;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	memset(&J, 0, sizeof J);
	J.A = A;
	box = (qaws_bp_box*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_bp_box), (size_t)A->ns + 1));
	if (!box)
		goto done;
	for (i = 0; i < A->ns; i++)
	{
		double const* s = A->seg + 4 * (size_t)i;
		box[i].lo[0] = fmin(s[0], s[2]) - eps; box[i].hi[0] = fmax(s[0], s[2]) + eps;
		box[i].lo[1] = fmin(s[1], s[3]) - eps; box[i].hi[1] = fmax(s[1], s[3]) + eps;
		box[i].lo[2] = box[i].hi[2] = 0;
	}
	st = qaws_internal_broadphase(NULL, box, A->ns, 2, NULL, sh_cut_visit, &J);
	if (J.oom || st != QAWS_STATUS_OK)
	{
		st = J.oom ? QAWS_STATUS_ALLOCATION_FAILURE : st;
		goto done;
	}
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	/* every segment's points: its ends and its cuts, in order along it */
	for (i = 0; i < A->ns; i++)
	{
		if (!sh_cut_push(&J, i, 0) || !sh_cut_push(&J, i, 1))
			goto done;
	}
	qsort(J.cut, J.nc, 2 * sizeof(double), sh_cmp_dbl2);
	npt = J.nc;
	pt = (sh_pt*)qaws_internal_malloc(qaws_internal_mul2(sizeof(sh_pt), (size_t)npt + 1));
	par = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)npt + 1));
	pid = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)npt + 1));
	if (!pt || !par || !pid)
		goto done;
	for (i = 0; i < npt; i++)
	{
		double const* s = A->seg + 4 * (size_t)J.cut[2 * (size_t)i];
		double t = J.cut[2 * (size_t)i + 1];
		pt[i].x = s[0] + t * (s[2] - s[0]);
		pt[i].y = s[1] + t * (s[3] - s[1]);
		pt[i].i = i;
		par[i] = i;
	}
	/* (points within eps: one vertex) */
	{
		sh_pt* srt = (sh_pt*)qaws_internal_malloc(qaws_internal_mul2(sizeof(sh_pt), (size_t)npt + 1));
		if (!srt)
			goto done;
		memcpy(srt, pt, sizeof(sh_pt) * npt);
		qsort(srt, npt, sizeof *srt, sh_cmp_pt);
		for (i = 0; i < npt; i++)
			for (j = i + 1; j < npt && srt[j].x - srt[i].x <= eps; j++)
				if (fabs(srt[j].y - srt[i].y) <= eps)
				{
					unsigned int a = sh_find(par, srt[i].i), b = sh_find(par, srt[j].i);
					if (a != b)
						par[a < b ? b : a] = a < b ? a : b;
				}
		qaws_internal_free(srt);
	}
	A->nv = 0;
	for (i = 0; i < npt; i++)
		pid[i] = SH_NONE;
	for (i = 0; i < npt; i++)
	{
		unsigned int r = sh_find(par, i);
		if (pid[r] == SH_NONE)
			pid[r] = A->nv++;
		pid[i] = pid[r];
	}
	A->vx = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 2, (size_t)A->nv + 1));
	ed = (unsigned long long*)qaws_internal_malloc(sizeof(unsigned long long) * ((size_t)npt + 1));
	if (!A->vx || !ed)
		goto done;
	for (i = 0; i < npt; i++)
	{
		A->vx[2 * (size_t)pid[i]] = pt[i].x;
		A->vx[2 * (size_t)pid[i] + 1] = pt[i].y;
	}
	/* the edges between consecutive points of a segment, once each */
	for (i = 0; i + 1 < npt; i++)
		if (J.cut[2 * (size_t)i] == J.cut[2 * (size_t)i + 2])
		{
			unsigned long long a = pid[i], b = pid[i + 1];
			if (a != b)
				ed[ne++] = a < b ? (a << 32 | b) : (b << 32 | a);
		}
	qsort(ed, ne, sizeof *ed, sh_cmp_u64);
	for (i = 0, k = 0; i < ne; i++)
		if (!k || ed[k - 1] != ed[i])
			ed[k++] = ed[i];
	ne = k;
	/* (edges with a free end pruned, again and again) */
	deg = (unsigned int*)qaws_internal_zalloc((size_t)A->nv + 1, sizeof(unsigned int));
	seen = (unsigned char*)qaws_internal_zalloc((size_t)ne + 1, 1);
	if (!deg || !seen)
		goto done;
	for (i = 0; i < ne; i++)
	{
		deg[ed[i] >> 32]++;
		deg[ed[i] & 0xFFFFFFFFu]++;
	}
	for (;;)
	{
		unsigned int pruned = 0;
		for (i = 0; i < ne; i++)
			if (!seen[i] && (deg[ed[i] >> 32] < 2 || deg[ed[i] & 0xFFFFFFFFu] < 2))
			{
				seen[i] = 1;
				deg[ed[i] >> 32]--;
				deg[ed[i] & 0xFFFFFFFFu]--;
				pruned++;
			}
		if (!pruned)
			break;
	}
	/* half-edges by vertex, counter-clockwise */
	he = (sh_out_edge*)qaws_internal_malloc(qaws_internal_mul3(sizeof(sh_out_edge), 2, (size_t)ne + 1));
	if (!he)
		goto done;
	for (i = 0, nh = 0; i < ne; i++)
		if (!seen[i])
		{
			unsigned int a = (unsigned int)(ed[i] >> 32), b = (unsigned int)(ed[i] & 0xFFFFFFFFu);
			he[nh].from = a; he[nh].to = b;
			he[nh++].angle = atan2(A->vx[2 * (size_t)b + 1] - A->vx[2 * (size_t)a + 1], A->vx[2 * (size_t)b] - A->vx[2 * (size_t)a]);
			he[nh].from = b; he[nh].to = a;
			he[nh++].angle = atan2(A->vx[2 * (size_t)a + 1] - A->vx[2 * (size_t)b + 1], A->vx[2 * (size_t)a] - A->vx[2 * (size_t)b]);
		}
	qsort(he, nh, sizeof *he, sh_cmp_out);
	first = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)A->nv + 2));
	nxt = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nh + 1));
	qaws_internal_free(seen);
	seen = (unsigned char*)qaws_internal_zalloc((size_t)nh + 1, 1);
	A->loop_first = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nh + 2));
	A->loop_v = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nh + 1));
	A->loop_area = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)nh + 1));
	A->hole_of = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nh + 1));
	if (!first || !nxt || !seen || !A->loop_first || !A->loop_v || !A->loop_area || !A->hole_of)
		goto done;
	for (i = 0, j = 0; i <= A->nv; i++)
	{
		while (j < nh && he[j].from < i)
			j++;
		first[i] = j;
	}
	first[A->nv + 1] = nh;
	/* next of (u -> v): at v, the half-edge before (v -> u) counter-clockwise */
	for (i = 0; i < nh; i++)
	{
		unsigned int v = he[i].to, u = he[i].from, lo = first[v], hi = first[v + 1], q;
		for (q = lo; q < hi && he[q].to != u; q++)
			;
		nxt[i] = q == lo ? hi - 1 : q - 1;
	}
	A->nloop = 0;
	A->loop_first[0] = 0;
	for (i = 0, k = 0; i < nh; i++)
	{
		unsigned int h = i;
		double area = 0;
		if (seen[i])
			continue;
		while (!seen[h])
		{
			unsigned int a = he[h].from, b = he[h].to;
			seen[h] = 1;
			A->loop_v[k++] = a;
			area += A->vx[2 * (size_t)a] * A->vx[2 * (size_t)b + 1] - A->vx[2 * (size_t)b] * A->vx[2 * (size_t)a + 1];
			h = nxt[h];
		}
		A->loop_area[A->nloop] = 0.5 * area;
		A->loop_first[++A->nloop] = k;
	}
	/* holes: each negative cycle in the smallest positive one around a point just left of its leftmost vertex */
	{
		double* poly = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 2, (size_t)k + 1));
		if (!poly)
			goto done;
		for (i = 0; i < A->nloop; i++)
		{
			unsigned int best = SH_NONE, q, lv = A->loop_v[A->loop_first[i]];
			double px, py;
			A->hole_of[i] = SH_NONE;
			if (A->loop_area[i] >= 0)
				continue;
			for (q = A->loop_first[i]; q < A->loop_first[i + 1]; q++)
				if (A->vx[2 * (size_t)A->loop_v[q]] < A->vx[2 * (size_t)lv])
					lv = A->loop_v[q];
			px = A->vx[2 * (size_t)lv] - 4 * eps;
			py = A->vx[2 * (size_t)lv + 1];
			for (j = 0; j < A->nloop; j++)
			{
				unsigned int m = 0;
				if (A->loop_area[j] <= 0 || (best != SH_NONE && A->loop_area[j] >= A->loop_area[best]))
					continue;
				for (q = A->loop_first[j]; q < A->loop_first[j + 1]; q++, m++)
				{
					poly[2 * m] = A->vx[2 * (size_t)A->loop_v[q]];
					poly[2 * m + 1] = A->vx[2 * (size_t)A->loop_v[q] + 1];
				}
				if (sh_wind(poly, m, px, py))
					best = j;
			}
			A->hole_of[i] = best;
		}
		qaws_internal_free(poly);
	}
	st = QAWS_STATUS_OK;
done:
	qaws_internal_free(box);
	qaws_internal_free(J.cut);
	qaws_internal_free(pt);
	qaws_internal_free(par);
	qaws_internal_free(pid);
	qaws_internal_free(deg);
	qaws_internal_free(nxt);
	qaws_internal_free(first);
	qaws_internal_free(ed);
	qaws_internal_free(he);
	qaws_internal_free(cuts);
	qaws_internal_free(seen);
	return st;
}

/* a point inside face f and outside its holes: an edge's middle moved left of it, the first edge that gives one */
static int sh_arr_inner(sh_arr const* A, unsigned int f, double eps, double* out)
{
	unsigned int q, n = A->loop_first[f + 1] - A->loop_first[f], h, m;
	double* poly = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 2, (size_t)A->loop_first[A->nloop] + 1));
	int found = 0;
	if (!poly)
		return 0;
	for (q = 0; q < n && !found; q++)
	{
		unsigned int a = A->loop_v[A->loop_first[f] + q], b = A->loop_v[A->loop_first[f] + (q + 1) % n];
		double ex = A->vx[2 * (size_t)b] - A->vx[2 * (size_t)a], ey = A->vx[2 * (size_t)b + 1] - A->vx[2 * (size_t)a + 1], l = sqrt(ex * ex + ey * ey);
		double px, py;
		int ok;
		if (!(l > 0))
			continue;
		px = 0.5 * (A->vx[2 * (size_t)a] + A->vx[2 * (size_t)b]) - ey / l * fmin(eps, 0.25 * l);
		py = 0.5 * (A->vx[2 * (size_t)a + 1] + A->vx[2 * (size_t)b + 1]) + ex / l * fmin(eps, 0.25 * l);
		for (m = 0; m < n; m++)
		{
			poly[2 * m] = A->vx[2 * (size_t)A->loop_v[A->loop_first[f] + m]];
			poly[2 * m + 1] = A->vx[2 * (size_t)A->loop_v[A->loop_first[f] + m] + 1];
		}
		ok = sh_wind(poly, n, px, py) != 0;
		for (h = 0; ok && h < A->nloop; h++)
			if (A->hole_of[h] == f)
			{
				unsigned int nh2 = A->loop_first[h + 1] - A->loop_first[h];
				for (m = 0; m < nh2; m++)
				{
					poly[2 * m] = A->vx[2 * (size_t)A->loop_v[A->loop_first[h] + m]];
					poly[2 * m + 1] = A->vx[2 * (size_t)A->loop_v[A->loop_first[h] + m] + 1];
				}
				ok = sh_wind(poly, nh2, px, py) == 0;
			}
		if (ok)
		{
			out[0] = px;
			out[1] = py;
			found = 1;
		}
	}
	qaws_internal_free(poly);
	return found;
}

/* ---- shatter ---------------------------------------------------------------------------------------------------------- */

struct qaws_surface_shatter
{
	qaws_surface** piece;
	unsigned int* from;
	int* inside;
	unsigned int count, cap;
	qaws_curve** curve;             /* the pieces' trim polygons, owned */
	unsigned int ncurve, capc;
};

void qaws_surface_shatter_destroy(qaws_surface_shatter* shatter)
{
	unsigned int i;
	if (!shatter)
		return;
	for (i = 0; i < shatter->count; i++)
		qaws_surface_destroy(shatter->piece[i]);
	for (i = 0; i < shatter->ncurve; i++)
		qaws_curve_destroy(shatter->curve[i]);
	qaws_internal_free(shatter->piece);
	qaws_internal_free(shatter->from);
	qaws_internal_free(shatter->inside);
	qaws_internal_free(shatter->curve);
	qaws_internal_free(shatter);
}

/* m edges of cycle loop from its point s (wrapping) as a degree-1 B-spline in (u, v), owned by H */
static qaws_curve* sh_polygon_curve(qaws_surface_shatter* H, sh_arr const* A, unsigned int loop, unsigned int s, unsigned int m)
{
	unsigned int n = A->loop_first[loop + 1] - A->loop_first[loop], i;
	qaws_vec2* cp = (qaws_vec2*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_vec2), (size_t)m + 2));
	qaws_bspline_desc d;
	qaws_curve* c = NULL;
	if (!cp)
		return NULL;
	for (i = 0; i <= m; i++)
	{
		unsigned int v = A->loop_v[A->loop_first[loop] + (s + i) % n];
		cp[i].x = (qaws_scalar)A->vx[2 * (size_t)v];
		cp[i].y = (qaws_scalar)A->vx[2 * (size_t)v + 1];
	}
	memset(&d, 0, sizeof d);
	d.dimension = QAWS_DIMENSION_2D;
	d.degree = 1;
	d.control_points = cp;
	d.control_point_count = m + 1;
	d.is_uniform = 1;
	if (qaws_curve_create_bspline(&c, &d) != QAWS_STATUS_OK)
		c = NULL;
	qaws_internal_free(cp);
	if (c && H->ncurve == H->capc)
	{
		unsigned int nc = H->capc ? 2 * H->capc : 32;
		qaws_curve** c2 = (qaws_curve**)qaws_internal_realloc(H->curve, qaws_internal_mul2(sizeof(qaws_curve*), nc));
		if (!c2)
		{
			qaws_curve_destroy(c);
			return NULL;
		}
		H->curve = c2;
		H->capc = nc;
	}
	if (c)
		H->curve[H->ncurve++] = c;
	return c;
}

/* cycle loop as curves of 50, 25, 20, 10, 5, 4, 2 or 1 edges into cv (as many as edges at most): a trim's containment
   test samples a curve at 100 even parameters, so every corner is one of them; the count, 0 on failure */
static unsigned int sh_polygon_curves(qaws_surface_shatter* H, sh_arr const* A, unsigned int loop, qaws_curve const** cv)
{
	static unsigned int const sizes[8] = { 50, 25, 20, 10, 5, 4, 2, 1 };
	unsigned int n = A->loop_first[loop + 1] - A->loop_first[loop], s = 0, c = 0, k = 0;
	while (s < n)
	{
		while (sizes[k] > n - s)
			k++;
		if (!(cv[c++] = sh_polygon_curve(H, A, loop, s, sizes[k])))
			return 0;
		s += sizes[k];
	}
	return c;
}

/* the winding number of the scene's B triangles (those inside their surface's trim) at p, by solid angles */
static double sh_winding_b(sh_scene const* S, double const* p)
{
	double w = 0;
	unsigned int t;
	for (t = 0; t < S->nt; t++)
	{
		sh_surf const* f = &S->surf[S->tsurf[t]];
		double a[3], b[3], c[3], la, lb, lc, num, den;
		unsigned int k;
		if (!f->operand)
			continue;
		if (f->s != f->base)
		{
			int in = 1;
			double u = 0, v = 0;
			for (k = 0; k < 3; k++)
			{
				u += S->uv[2 * (size_t)S->tri[3 * (size_t)t + k]] / 3;
				v += S->uv[2 * (size_t)S->tri[3 * (size_t)t + k] + 1] / 3;
			}
			qaws_surface_trim_contains(&in, f->s, (qaws_scalar)u, (qaws_scalar)v);
			if (!in)
				continue;
		}
		for (k = 0; k < 3; k++)
		{
			a[k] = S->x[3 * (size_t)S->tri[3 * (size_t)t] + k] - p[k];
			b[k] = S->x[3 * (size_t)S->tri[3 * (size_t)t + 1] + k] - p[k];
			c[k] = S->x[3 * (size_t)S->tri[3 * (size_t)t + 2] + k] - p[k];
		}
		la = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
		lb = sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
		lc = sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
		num = a[0] * (b[1] * c[2] - b[2] * c[1]) + a[1] * (b[2] * c[0] - b[0] * c[2]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
		den = la * lb * lc + (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) * lc + (b[0] * c[0] + b[1] * c[1] + b[2] * c[2]) * la
			+ (c[0] * a[0] + c[1] * a[1] + c[2] * a[2]) * lb;
		w += 2 * atan2(num, den);
	}
	return w / (4 * 3.14159265358979323846);
}

/* surface k of A cut along its seams in R into pieces, added to H */
static qaws_status sh_shatter_surface(qaws_surface_shatter* H, sh_scene const* S, qaws_surface_seams const* R, unsigned int k, int b_is_solid)
{
	sh_surf const* f = &S->surf[k];
	sh_arr A;
	unsigned int i, j, nl = qaws_internal_trim_loop_count(f->s);
	double du = f->u1 - f->u0, dv = f->v1 - f->v0, eps = 1e-9 * (du + dv), ext = 2e-3 * (du + dv);
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	memset(&A, 0, sizeof A);
	/* the domain: the trim loops (sampled) or the rectangle */
	if (nl)
	{
		for (i = 0; i < nl; i++)
		{
			unsigned int nc = 0, c;
			int outer = 0;
			qaws_curve const* const* cv = qaws_internal_trim_loop(&nc, &outer, f->s, i);
			double px = 0, py = 0, fx = 0, fy = 0;
			int have = 0;
			for (c = 0; c < nc; c++)
			{
				qaws_range r = qaws_curve_get_parameter_range(cv[c]);
				unsigned int s2;
				for (s2 = 0; s2 <= 64; s2++)
				{
					qaws_eval_result_2d e;
					double t = (double)r.min_value + ((double)r.max_value - (double)r.min_value) * s2 / 64;
					memset(&e, 0, sizeof e);
					if (qaws_curve_evaluate_2d(&e, cv[c], (qaws_scalar)t, QAWS_EVAL_FLAG_POSITION) != QAWS_STATUS_OK)
						continue;
					if (have && !sh_arr_add(&A, px, py, (double)e.position.x, (double)e.position.y))
						goto done;
					if (!have)
					{
						fx = (double)e.position.x;
						fy = (double)e.position.y;
					}
					px = (double)e.position.x;
					py = (double)e.position.y;
					have = 1;
				}
			}
			if (have && !sh_arr_add(&A, px, py, fx, fy))
				goto done;
		}
	}
	else if (!sh_arr_add(&A, f->u0, f->v0, f->u1, f->v0) || !sh_arr_add(&A, f->u1, f->v0, f->u1, f->v1) || !sh_arr_add(&A, f->u1, f->v1, f->u0, f->v1)
		|| !sh_arr_add(&A, f->u0, f->v1, f->u0, f->v0))
		goto done;
	/* the A-B seams on it, open ones carried a little past their ends (to cross the domain's border they stop on) */
	for (i = 0; i < R->count; i++)
	{
		sh_seam const* s = &R->seam[i];
		unsigned int n = s->count;
		qaws_seam_point const* P = R->pt + s->first;
		if (s->p != k || S->surf[s->q].operand != 1)
			continue;
		for (j = 0; j + 1 < n + (s->closed ? 1u : 0u); j++)
		{
			double x0 = (double)P[j].u1, y0 = (double)P[j].v1, x1 = (double)P[(j + 1) % n].u1, y1 = (double)P[(j + 1) % n].v1;
			if (!s->closed && (j == 0 || j + 2 == n))
			{
				double ex = x1 - x0, ey = y1 - y0, l = sqrt(ex * ex + ey * ey);
				if (l > 0 && j == 0)
				{
					x0 -= ex / l * ext;
					y0 -= ey / l * ext;
				}
				if (l > 0 && j + 2 == n)
				{
					x1 += ex / l * ext;
					y1 += ey / l * ext;
				}
			}
			if (!sh_arr_add(&A, x0, y0, x1, y1))
				goto done;
		}
	}
	st = sh_arr_build(&A, eps);
	if (st != QAWS_STATUS_OK)
		goto done;
	/* each face inside the domain: a piece */
	for (i = 0; i < A.nloop; i++)
	{
		double pin[2], P3[3];
		int in = 1, inside = -1;
		unsigned int nh = 0, h, ncv;
		qaws_trim_loop* loops;
		qaws_curve const** cv;
		qaws_surface_trim_desc td;
		qaws_surface* piece = NULL;
		if (A.loop_area[i] <= 0 || !sh_arr_inner(&A, i, 10 * eps + 1e-6 * (du + dv), pin))
			continue;
		if (nl)
			qaws_surface_trim_contains(&in, f->s, (qaws_scalar)pin[0], (qaws_scalar)pin[1]);
		else
			in = pin[0] > f->u0 && pin[0] < f->u1 && pin[1] > f->v0 && pin[1] < f->v1;
		if (!in)
			continue;
		ncv = A.loop_first[i + 1] - A.loop_first[i];
		for (h = 0; h < A.nloop; h++)
			if (A.hole_of[h] == i)
			{
				nh++;
				ncv += A.loop_first[h + 1] - A.loop_first[h];
			}
		loops = (qaws_trim_loop*)qaws_internal_zalloc((size_t)nh + 1, sizeof(qaws_trim_loop));
		cv = (qaws_curve const**)qaws_internal_zalloc((size_t)ncv + 1, sizeof(qaws_curve const*));
		if (!loops || !cv)
		{
			qaws_internal_free(loops);
			qaws_internal_free(cv);
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			goto done;
		}
		{
			unsigned int at = 0, got;
			got = sh_polygon_curves(H, &A, i, cv);
			loops[0].curves = cv;
			loops[0].curve_count = got;
			loops[0].is_outer = 1;
			if (!got)
				st = QAWS_STATUS_ALLOCATION_FAILURE;
			at = got;
			for (h = 0, j = 1; h < A.nloop; h++)
				if (A.hole_of[h] == i)
				{
					got = sh_polygon_curves(H, &A, h, cv + at);
					loops[j].curves = cv + at;
					loops[j].curve_count = got;
					loops[j++].is_outer = 0;
					if (!got)
						st = QAWS_STATUS_ALLOCATION_FAILURE;
					at += got;
				}
		}
		memset(&td, 0, sizeof td);
		td.base = f->base;
		td.loops = loops;
		td.loop_count = nh + 1;
		if (st != QAWS_STATUS_ALLOCATION_FAILURE)
			st = qaws_surface_create_trimmed(&piece, &td);
		qaws_internal_free(loops);
		qaws_internal_free(cv);
		if (st != QAWS_STATUS_OK)
			goto done;
		if (b_is_solid && sh_eval(P3, NULL, NULL, f->base, pin[0], pin[1]))
			inside = sh_winding_b(S, P3) >= 0.5;
		if (H->count == H->cap)
		{
			unsigned int nc = H->cap ? 2 * H->cap : 16;
			qaws_surface** p2 = (qaws_surface**)qaws_internal_realloc(H->piece, qaws_internal_mul2(sizeof(qaws_surface*), nc));
			unsigned int* f2 = p2 ? (unsigned int*)qaws_internal_realloc(H->from, qaws_internal_mul2(sizeof(unsigned int), nc)) : NULL;
			int* i2 = f2 ? (int*)qaws_internal_realloc(H->inside, qaws_internal_mul2(sizeof(int), nc)) : NULL;
			if (p2)
				H->piece = p2;
			if (f2)
				H->from = f2;
			if (i2)
				H->inside = i2;
			if (!i2)
			{
				qaws_surface_destroy(piece);
				st = QAWS_STATUS_ALLOCATION_FAILURE;
				goto done;
			}
			H->cap = nc;
		}
		H->piece[H->count] = piece;
		H->from[H->count] = f->index;
		H->inside[H->count++] = inside;
	}
	st = QAWS_STATUS_OK;
done:
	sh_arr_free(&A);
	return st;
}

qaws_status qaws_surface_shatter_create(qaws_surface_shatter** out_shatter, qaws_surface_seam_desc const* desc, int b_is_solid)
{
	sh_scene S;
	qaws_surface_seams* R = NULL;
	qaws_surface_shatter* H = NULL;
	qaws_status st;
	unsigned int k;
	if (!out_shatter)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_shatter = NULL;
	if (!sh_desc_ok(desc) || !desc->a_count || !desc->b_count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	st = sh_scene_build(&S, desc);
	if (st == QAWS_STATUS_OK)
		st = sh_seams(&R, &S, QAWS_SEAM_A_B);
	if (st == QAWS_STATUS_OK && !(H = (qaws_surface_shatter*)qaws_internal_zalloc(1, sizeof *H)))
		st = QAWS_STATUS_ALLOCATION_FAILURE;
	for (k = 0; st == QAWS_STATUS_OK && k < desc->a_count; k++)
		st = sh_shatter_surface(H, &S, R, k, b_is_solid);
	if (st == QAWS_STATUS_OK)
	{
		*out_shatter = H;
		H = NULL;
	}
	qaws_surface_shatter_destroy(H);
	qaws_surface_seams_destroy(R);
	sh_scene_free(&S);
	return st;
}

unsigned int qaws_surface_shatter_get_count(qaws_surface_shatter const* shatter)
{
	return shatter ? shatter->count : 0;
}

qaws_status qaws_surface_shatter_get(qaws_surface const** out_piece, unsigned int* out_surface, int* out_inside,
	qaws_surface_shatter const* shatter, unsigned int index)
{
	if (!shatter || index >= shatter->count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (out_piece)
		*out_piece = shatter->piece[index];
	if (out_surface)
		*out_surface = shatter->from[index];
	if (out_inside)
		*out_inside = shatter->inside[index];
	return QAWS_STATUS_OK;
}
