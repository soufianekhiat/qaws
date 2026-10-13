/*
 * qaws_shell_region.c -- curves and regions across the faces of a shell:
 * the atlas of shared sides, straightest geodesics continued over edges,
 * geodesic circles across faces, regions face by face.
 */

#include "qaws_shell_region.h"
#include "qaws_surface.h"
#include "qaws_curve.h"
#include "qaws_bezier.h"
#include "qaws_hermite.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "qaws_operations.h"
#include "qaws_path.h"
#include "internal/qaws_internal_surface_region.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define SH_PI 3.14159265358979323846
#define SH_NONE 0xffffffffu

#if QAWS_SCALAR_IS_FLOAT
#define SH_TOL 1e-4
#define SH_ON_SIDE 1e-5     /* a loop's stretch on a side within rounding */
#else
#define SH_TOL 1e-7
#define SH_ON_SIDE 1e-9
#endif

/* ------------------------------------------------------------------------ */
/*  The atlas                                                               */
/* ------------------------------------------------------------------------ */

/* a face's sides: 0 v = v0 (bottom), 1 u = u1 (right), 2 v = v1 (top),
   3 u = u0 (left); each parametrised by s in [0, 1] with u or v growing */
struct qaws_shell_atlas
{
	unsigned int nf;
	qaws_surface const** s;
	qaws_internal_sg_ctx* c;
	unsigned int (*nb)[4];      /* the neighbour face across each side, or none */
	int (*nbside)[4];
	int (*rev)[4];              /* the neighbour's side runs the other way */
	unsigned int ne;
	double tol;                 /* in space */
};

static void sh_side_uv(qaws_internal_sg_ctx const* c, int side, double s, double* uv)
{
	switch (side)
	{
	case 0: uv[0] = c->u0 + s * (c->u1 - c->u0); uv[1] = c->v0; break;
	case 1: uv[0] = c->u1; uv[1] = c->v0 + s * (c->v1 - c->v0); break;
	case 2: uv[0] = c->u0 + s * (c->u1 - c->u0); uv[1] = c->v1; break;
	default: uv[0] = c->u0; uv[1] = c->v0 + s * (c->v1 - c->v0); break;
	}
}

static qaws_status sh_point(qaws_internal_sg_ctx const* c, double const* uv, unsigned int flags, qaws_surface_eval_result* r)
{
	return qaws_internal_sg_eval(c, uv[0], uv[1], flags, r);
}

static double sh_dist(qaws_vec3 a, qaws_vec3 b)
{
	return sqrt(((double)a.x - b.x) * ((double)a.x - b.x) + ((double)a.y - b.y) * ((double)a.y - b.y) + ((double)a.z - b.z) * ((double)a.z - b.z));
}

/* a side that is an edge: not along a closed direction, not one point */
static int sh_side_real(qaws_internal_sg_ctx const* c, int side, double tol)
{
	qaws_surface_eval_result a, b, m;
	double uv[2];
	if ((side == 1 || side == 3) && c->pu > 0) return 0;
	if ((side == 0 || side == 2) && c->pv > 0) return 0;
	sh_side_uv(c, side, 0, uv); sh_point(c, uv, QAWS_SURFACE_EVAL_POSITION, &a);
	sh_side_uv(c, side, 0.5, uv); sh_point(c, uv, QAWS_SURFACE_EVAL_POSITION, &m);
	sh_side_uv(c, side, 1, uv); sh_point(c, uv, QAWS_SURFACE_EVAL_POSITION, &b);
	return sh_dist(a.position, m.position) > tol || sh_dist(b.position, m.position) > tol;
}

/* the sides coincide in space: +1 the same way, -1 reversed, 0 not */
static int sh_sides_match(qaws_internal_sg_ctx const* a, int sa, qaws_internal_sg_ctx const* b, int sb, double tol)
{
	static double const at[5] = { 0, 0.25, 0.5, 0.75, 1 };
	int same = 1, opp = 1;
	unsigned int k;
	for (k = 0; k < 5; k++)
	{
		qaws_surface_eval_result p, q, r;
		double uv[2];
		sh_side_uv(a, sa, at[k], uv); sh_point(a, uv, QAWS_SURFACE_EVAL_POSITION, &p);
		sh_side_uv(b, sb, at[k], uv); sh_point(b, uv, QAWS_SURFACE_EVAL_POSITION, &q);
		sh_side_uv(b, sb, 1 - at[k], uv); sh_point(b, uv, QAWS_SURFACE_EVAL_POSITION, &r);
		if (sh_dist(p.position, q.position) > tol) same = 0;
		if (sh_dist(p.position, r.position) > tol) opp = 0;
	}
	return same ? 1 : (opp ? -1 : 0);
}

/* the counter-clockwise boundary runs a side with s growing (bottom,
   right) or falling (top, left) */
static int sh_side_ccw(int side)
{
	return side == 0 || side == 1 ? 1 : -1;
}

qaws_status qaws_shell_atlas_create(qaws_surface const* const* surfaces, unsigned int n, qaws_scalar tolerance,
	qaws_shell_atlas** out_atlas)
{
	qaws_shell_atlas* a;
	unsigned int i, j;
	int si, sj;
	double size = 0;
	qaws_status st = QAWS_STATUS_OK;
	if (!surfaces || !n || !out_atlas)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_atlas = NULL;
	a = (qaws_shell_atlas*)calloc(1, sizeof(*a));
	if (!a)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	a->nf = n;
	a->s = (qaws_surface const**)malloc(sizeof(qaws_surface const*) * n);
	a->c = (qaws_internal_sg_ctx*)calloc(n, sizeof(qaws_internal_sg_ctx));
	a->nb = (unsigned int(*)[4])malloc(sizeof(unsigned int[4]) * n);
	a->nbside = (int(*)[4])malloc(sizeof(int[4]) * n);
	a->rev = (int(*)[4])malloc(sizeof(int[4]) * n);
	if (!a->s || !a->c || !a->nb || !a->nbside || !a->rev)
	{
		qaws_shell_atlas_destroy(a);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < n && st == QAWS_STATUS_OK; i++)
	{
		if (!surfaces[i]) { st = QAWS_STATUS_INVALID_ARGUMENT; break; }
		a->s[i] = surfaces[i];
		st = qaws_internal_sg_init(&a->c[i], surfaces[i], 0);
		for (si = 0; si < 4; si++) { a->nb[i][si] = SH_NONE; a->nbside[i][si] = -1; a->rev[i][si] = 0; }
	}
	/* the shell's size: its faces' corners */
	for (i = 0; i < n && st == QAWS_STATUS_OK; i++)
		for (si = 0; si < 4; si++)
		{
			qaws_surface_eval_result p, q;
			double uv[2];
			sh_side_uv(&a->c[i], si, 0, uv); sh_point(&a->c[i], uv, QAWS_SURFACE_EVAL_POSITION, &p);
			sh_side_uv(&a->c[(i + 1) % n], (si + 2) % 4, 0.5, uv); sh_point(&a->c[(i + 1) % n], uv, QAWS_SURFACE_EVAL_POSITION, &q);
			if (sh_dist(p.position, q.position) > size) size = sh_dist(p.position, q.position);
		}
	a->tol = tolerance > 0 ? (double)tolerance : SH_TOL * (size > 0 ? size : 1);
	/* the shared sides */
	for (i = 0; i < n && st == QAWS_STATUS_OK; i++)
		for (si = 0; si < 4 && st == QAWS_STATUS_OK; si++)
		{
			if (a->nb[i][si] != SH_NONE || !sh_side_real(&a->c[i], si, a->tol))
				continue;
			for (j = i; j < n && a->nb[i][si] == SH_NONE; j++)
				for (sj = 0; sj < 4 && a->nb[i][si] == SH_NONE; sj++)
				{
					int m;
					if ((j == i && sj == si) || a->nb[j][sj] != SH_NONE || !sh_side_real(&a->c[j], sj, a->tol))
						continue;
					if (!(m = sh_sides_match(&a->c[i], si, &a->c[j], sj, a->tol)))
						continue;
					/* oriented alike: the edge runs opposite ways round the two faces */
					if (sh_side_ccw(si) * sh_side_ccw(sj) * m != -1)
					{
						st = QAWS_STATUS_INVALID_ARGUMENT;
						break;
					}
					a->nb[i][si] = j; a->nbside[i][si] = sj; a->rev[i][si] = m < 0;
					a->nb[j][sj] = i; a->nbside[j][sj] = si; a->rev[j][sj] = m < 0;
					a->ne++;
				}
		}
	if (st != QAWS_STATUS_OK)
	{
		qaws_shell_atlas_destroy(a);
		return st;
	}
	*out_atlas = a;
	return QAWS_STATUS_OK;
}

qaws_status qaws_shell_atlas_create_from_brep(qaws_brep_shell const* shell, qaws_scalar tolerance, qaws_shell_atlas** out_atlas)
{
	unsigned int n, i;
	qaws_surface const** s;
	qaws_status st = QAWS_STATUS_OK;
	if (!shell || !out_atlas)
		return QAWS_STATUS_INVALID_ARGUMENT;
	n = qaws_brep_get_face_count(shell);
	s = (qaws_surface const**)malloc(sizeof(qaws_surface const*) * (n + 1));
	if (!s)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < n && st == QAWS_STATUS_OK; i++)
		st = qaws_brep_get_face_surface(shell, i, &s[i]);
	if (st == QAWS_STATUS_OK)
		st = qaws_shell_atlas_create(s, n, tolerance, out_atlas);
	free(s);
	return st;
}

void qaws_shell_atlas_destroy(qaws_shell_atlas* a)
{
	if (!a)
		return;
	free(a->s); free(a->c); free(a->nb); free(a->nbside); free(a->rev);
	free(a);
}

unsigned int qaws_shell_atlas_get_face_count(qaws_shell_atlas const* a)
{
	return a ? a->nf : 0;
}

qaws_surface const* qaws_shell_atlas_get_face(qaws_shell_atlas const* a, unsigned int face)
{
	return a && face < a->nf ? a->s[face] : NULL;
}

unsigned int qaws_shell_atlas_get_edge_count(qaws_shell_atlas const* a)
{
	return a ? a->ne : 0;
}

/* ------------------------------------------------------------------------ */
/*  Geodesics over edges                                                    */
/* ------------------------------------------------------------------------ */

/* inside the face's domain (closed directions never leave it) */
static int sh_inside(qaws_internal_sg_ctx const* c, double const* uv)
{
	double eu = 1e-12 * (c->u1 - c->u0), ev = 1e-12 * (c->v1 - c->v0);
	if (c->pu <= 0 && (uv[0] < c->u0 - eu || uv[0] > c->u1 + eu)) return 0;
	if (c->pv <= 0 && (uv[1] < c->v0 - ev || uv[1] > c->v1 + ev)) return 0;
	return 1;
}

/* the side a point just past the domain left by, and its parameter there */
static int sh_exit_side(qaws_internal_sg_ctx const* c, double const* uv, double* s)
{
	double du0 = c->pu > 0 ? 0 : (c->u0 - uv[0]) / (c->u1 - c->u0), du1 = c->pu > 0 ? 0 : (uv[0] - c->u1) / (c->u1 - c->u0);
	double dv0 = c->pv > 0 ? 0 : (c->v0 - uv[1]) / (c->v1 - c->v0), dv1 = c->pv > 0 ? 0 : (uv[1] - c->v1) / (c->v1 - c->v0);
	int side = 0;
	double best = dv0;
	if (du1 > best) { best = du1; side = 1; }
	if (dv1 > best) { best = dv1; side = 2; }
	if (du0 > best) { best = du0; side = 3; }
	if (side == 0 || side == 2) *s = (uv[0] - c->u0) / (c->u1 - c->u0);
	else *s = (uv[1] - c->v0) / (c->v1 - c->v0);
	if (*s < 0) *s = 0;
	if (*s > 1) *s = 1;
	return side;
}

static void sh_add3(double* o, qaws_vec3 a, double ka, qaws_vec3 b, double kb)
{
	o[0] = ka * a.x + kb * b.x; o[1] = ka * a.y + kb * b.y; o[2] = ka * a.z + kb * b.z;
}

static double sh_dot(double const* a, double const* b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void sh_unit(double* a)
{
	double n = sqrt(sh_dot(a, a));
	if (n > 0) { a[0] /= n; a[1] /= n; a[2] /= n; }
}

/* the side's tangent in space (s growing) and the in-plane conormal
   pointing out of the domain there */
static qaws_status sh_side_frame(qaws_internal_sg_ctx const* c, int side, double const* uv, double* E, double* N)
{
	qaws_surface_eval_result r;
	qaws_status st = sh_point(c, uv, QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV, &r);
	double O[3];
	if (st != QAWS_STATUS_OK)
		return st;
	if (side == 0 || side == 2) { sh_add3(E, r.du, 1, r.dv, 0); sh_add3(O, r.dv, side == 0 ? -1 : 1, r.du, 0); }
	else { sh_add3(E, r.dv, 1, r.du, 0); sh_add3(O, r.du, side == 3 ? -1 : 1, r.dv, 0); }
	sh_unit(E);
	{
		double k = sh_dot(O, E);
		N[0] = O[0] - k * E[0]; N[1] = O[1] - k * E[1]; N[2] = O[2] - k * E[2];
		sh_unit(N);
	}
	return QAWS_STATUS_OK;
}

/* the neighbour's side parameter at the point x, from s0, by Gauss-Newton */
static double sh_invert_side(qaws_internal_sg_ctx const* c, int side, double s0, qaws_vec3 x)
{
	double s = s0;
	unsigned int it;
	for (it = 0; it < 30; it++)
	{
		qaws_surface_eval_result r;
		double uv[2], d[3], p[3], dd, step;
		sh_side_uv(c, side, s, uv);
		if (sh_point(c, uv, QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV, &r) != QAWS_STATUS_OK)
			break;
		if (side == 0 || side == 2) sh_add3(d, r.du, c->u1 - c->u0, r.dv, 0);
		else sh_add3(d, r.dv, c->v1 - c->v0, r.du, 0);
		p[0] = r.position.x - x.x; p[1] = r.position.y - x.y; p[2] = r.position.z - x.z;
		dd = sh_dot(d, d);
		if (!(dd > 0))
			break;
		step = sh_dot(p, d) / dd;
		s -= step;
		if (s < 0) s = 0;
		if (s > 1) s = 1;
		if (fabs(step) < 1e-15)
			break;
	}
	return s;
}

/* across the side: y (u, v, u', v') on face f at its side `side`, s, moved to
   the neighbour with the straightest continuation; 0 when there is none */
static int sh_cross(qaws_shell_atlas const* a, unsigned int* f, int side, double s, double* y)
{
	qaws_internal_sg_ctx const* c = &a->c[*f];
	unsigned int g = a->nb[*f][side];
	int sg2 = a->nbside[*f][side];
	qaws_surface_eval_result r, q;
	double uv[2], T[3], EA[3], NA[3], EB[3], NB[3], al, be, Tn[3], g3[3], rhs[2], det, sb;
	if (g == SH_NONE)
		return 0;
	sh_side_uv(c, side, s, uv);
	if (sh_point(c, uv, QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV, &r) != QAWS_STATUS_OK ||
		sh_side_frame(c, side, uv, EA, NA) != QAWS_STATUS_OK)
		return 0;
	sh_add3(T, r.du, y[2], r.dv, y[3]);
	al = sh_dot(T, EA);
	be = sh_dot(T, NA);
	/* the neighbour's point and frame */
	sb = sh_invert_side(&a->c[g], sg2, a->rev[*f][side] ? 1 - s : s, r.position);
	sh_side_uv(&a->c[g], sg2, sb, uv);
	if (sh_side_frame(&a->c[g], sg2, uv, EB, NB) != QAWS_STATUS_OK ||
		sh_point(&a->c[g], uv, QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV, &q) != QAWS_STATUS_OK)
		return 0;
	if (sh_dot(EA, EB) < 0) { EB[0] = -EB[0]; EB[1] = -EB[1]; EB[2] = -EB[2]; }
	/* into the neighbour: against its outward conormal */
	Tn[0] = al * EB[0] - be * NB[0]; Tn[1] = al * EB[1] - be * NB[1]; Tn[2] = al * EB[2] - be * NB[2];
	{
		double su[3], sv[3];
		sh_add3(su, q.du, 1, q.dv, 0);
		sh_add3(sv, q.dv, 1, q.du, 0);
		g3[0] = sh_dot(su, su); g3[1] = sh_dot(su, sv); g3[2] = sh_dot(sv, sv);
		rhs[0] = sh_dot(su, Tn); rhs[1] = sh_dot(sv, Tn);
	}
	det = g3[0] * g3[2] - g3[1] * g3[1];
	if (!(det > 0))
		return 0;
	y[0] = uv[0]; y[1] = uv[1];
	y[2] = (g3[2] * rhs[0] - g3[1] * rhs[1]) / det;
	y[3] = (g3[0] * rhs[1] - g3[1] * rhs[0]) / det;
	*f = g;
	return 1;
}

/* a run of the geodesic on one face: its start state and length */
typedef struct sh_run
{
	unsigned int face;
	double y0[4];
	double len;
} sh_run;

typedef struct sh_runs
{
	sh_run* r;
	unsigned int n, cap;
} sh_runs;

static int sh_runs_add(sh_runs* rs, unsigned int f, double const* y0, double len)
{
	if (!rs)
		return 1;
	if (rs->n + 1 > rs->cap)
	{
		unsigned int cap = rs->cap ? 2 * rs->cap : 8;
		sh_run* q = (sh_run*)realloc(rs->r, sizeof(sh_run) * cap);
		if (!q)
			return 0;
		rs->r = q;
		rs->cap = cap;
	}
	rs->r[rs->n].face = f;
	memcpy(rs->r[rs->n].y0, y0, sizeof(double) * 4);
	rs->r[rs->n].len = len;
	rs->n++;
	return 1;
}

/* the geodesic from state y on face *f over the length L, across edges;
   y and *f become its end */
static qaws_status sh_trace(qaws_shell_atlas const* a, unsigned int* f, double* y, double L, sh_runs* runs, unsigned int* crossings)
{
	double left = L;
	unsigned int hops = 0;
	if (crossings) *crossings = 0;
	while (hops < 256)
	{
		qaws_internal_sg_ctx const* c = &a->c[*f];
		double y0[4], done = 0, speed = fabs(y[2]) + fabs(y[3]), step;
		int crossed = 0;
		memcpy(y0, y, sizeof(y0));
		step = speed > 0 ? 0.05 * c->diag / speed : left;
		if (!(step > 0)) step = left;
		while (done < left)
		{
			double h = step < left - done ? step : left - done, yt[4];
			memcpy(yt, y, sizeof(yt));
			if (qaws_internal_sg_integrate(c, yt, h) == QAWS_STATUS_OK && sh_inside(c, yt))
			{
				memcpy(y, yt, sizeof(yt));
				done += h;
				continue;
			}
			{
				/* it leaves within this step: the exit by bisection on the length */
				double lo = 0, hi = h, ylo[4], s;
				unsigned int it;
				int side;
				memcpy(ylo, y, sizeof(ylo));
				for (it = 0; it < 80 && hi - lo > 1e-15 * (L + 1); it++)
				{
					double mid = 0.5 * (lo + hi), ym[4];
					memcpy(ym, y, sizeof(ym));
					if (qaws_internal_sg_integrate(c, ym, mid) == QAWS_STATUS_OK && sh_inside(c, ym))
					{
						lo = mid;
						memcpy(ylo, ym, sizeof(ylo));
					}
					else
						hi = mid;
				}
				/* the side: from a point just past it */
				{
					double yp[4];
					memcpy(yp, y, sizeof(yp));
					qaws_internal_sg_integrate(c, yp, hi + 1e-9 * c->diag);
					side = sh_exit_side(c, yp, &s);
					/* along a closed direction: wrapped into its period (no corners there) */
					if (side == 0 || side == 2) s = c->pu > 0 ? (ylo[0] - c->u0) / c->pu - floor((ylo[0] - c->u0) / c->pu) : (ylo[0] - c->u0) / (c->u1 - c->u0);
					else s = c->pv > 0 ? (ylo[1] - c->v0) / c->pv - floor((ylo[1] - c->v0) / c->pv) : (ylo[1] - c->v0) / (c->v1 - c->v0);
					if (s < 0) s = 0;
					if (s > 1) s = 1;
				}
				done += lo;
				memcpy(y, ylo, sizeof(ylo));
				if (!sh_runs_add(runs, *f, y0, done))
					return QAWS_STATUS_ALLOCATION_FAILURE;
				/* a corner of the domain: a vertex of the shell */
				if ((s < 1e-9 || s > 1 - 1e-9) && !((side == 0 || side == 2) && c->pu > 0) && !((side == 1 || side == 3) && c->pv > 0))
					return QAWS_STATUS_OUT_OF_RANGE;
				if (!sh_cross(a, f, side, s, y))
					return QAWS_STATUS_OUT_OF_RANGE;
				left -= done;
				hops++;
				crossed = 1;
				if (crossings) (*crossings)++;
				break;
			}
		}
		if (!crossed)
			return sh_runs_add(runs, *f, y0, done) ? QAWS_STATUS_OK : QAWS_STATUS_ALLOCATION_FAILURE;
	}
	return QAWS_STATUS_NUMERICAL_FAILURE;
}

/* the start state: p with its unit (u, v) direction */
static qaws_status sh_start(qaws_shell_atlas const* a, qaws_shell_point p, double const* d, double* y)
{
	qaws_internal_sg_ctx const* c;
	double g[3], n;
	qaws_status st;
	if (!a || p.face >= a->nf)
		return QAWS_STATUS_INVALID_ARGUMENT;
	c = &a->c[p.face];
	if ((st = qaws_internal_sg_first_form(c, p.uv.x, p.uv.y, g)) != QAWS_STATUS_OK)
		return st;
	n = sqrt(g[0] * d[0] * d[0] + 2 * g[1] * d[0] * d[1] + g[2] * d[1] * d[1]);
	if (!(n > 0))
		return QAWS_STATUS_INVALID_ARGUMENT;
	y[0] = p.uv.x; y[1] = p.uv.y; y[2] = d[0] / n; y[3] = d[1] / n;
	return QAWS_STATUS_OK;
}

qaws_status qaws_shell_exp(qaws_shell_atlas const* a, qaws_shell_point p, qaws_vec2 direction, qaws_scalar length,
	qaws_shell_point* out_point, qaws_vec3* out_velocity, unsigned int* out_crossings)
{
	double y[4], d[2];
	unsigned int f;
	qaws_status st;
	d[0] = direction.x; d[1] = direction.y;
	if ((st = sh_start(a, p, d, y)) != QAWS_STATUS_OK)
		return st;
	f = p.face;
	if ((st = sh_trace(a, &f, y, (double)length, NULL, out_crossings)) != QAWS_STATUS_OK)
		return st;
	if (out_point)
	{
		out_point->face = f;
		out_point->uv.x = (qaws_scalar)y[0]; out_point->uv.y = (qaws_scalar)y[1];
	}
	if (out_velocity)
	{
		qaws_surface_eval_result r;
		if ((st = qaws_internal_sg_eval(&a->c[f], y[0], y[1], QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV, &r)) != QAWS_STATUS_OK)
			return st;
		out_velocity->x = (qaws_scalar)(r.du.x * y[2] + r.dv.x * y[3]);
		out_velocity->y = (qaws_scalar)(r.du.y * y[2] + r.dv.y * y[3]);
		out_velocity->z = (qaws_scalar)(r.du.z * y[2] + r.dv.z * y[3]);
	}
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------------ */
/*  Shell curves                                                            */
/* ------------------------------------------------------------------------ */

struct qaws_shell_curve
{
	qaws_shell_atlas const* a;
	unsigned int n, cap;
	unsigned int* face;
	qaws_surface_curve** c;
	qaws_curve** keep;              /* (u, v) curves made here, wrapped by c */
	int closed;
};

static qaws_shell_curve* sh_curve_new(qaws_shell_atlas const* a)
{
	qaws_shell_curve* sc = (qaws_shell_curve*)calloc(1, sizeof(qaws_shell_curve));
	if (sc) sc->a = a;
	return sc;
}

static qaws_status sh_curve_push(qaws_shell_curve* sc, unsigned int face, qaws_surface_curve* c, qaws_curve* keep)
{
	if (sc->n + 1 > sc->cap)
	{
		unsigned int cap = sc->cap ? 2 * sc->cap : 8;
		unsigned int* f = (unsigned int*)realloc(sc->face, sizeof(unsigned int) * cap);
		qaws_surface_curve** cc;
		qaws_curve** kk;
		if (!f) { qaws_surface_curve_destroy(c); qaws_curve_destroy(keep); return QAWS_STATUS_ALLOCATION_FAILURE; }
		sc->face = f;
		cc = (qaws_surface_curve**)realloc(sc->c, sizeof(qaws_surface_curve*) * cap);
		if (!cc) { qaws_surface_curve_destroy(c); qaws_curve_destroy(keep); return QAWS_STATUS_ALLOCATION_FAILURE; }
		sc->c = cc;
		kk = (qaws_curve**)realloc(sc->keep, sizeof(qaws_curve*) * cap);
		if (!kk) { qaws_surface_curve_destroy(c); qaws_curve_destroy(keep); return QAWS_STATUS_ALLOCATION_FAILURE; }
		sc->keep = kk;
		sc->cap = cap;
	}
	sc->face[sc->n] = face;
	sc->c[sc->n] = c;
	sc->keep[sc->n] = keep;
	sc->n++;
	return QAWS_STATUS_OK;
}

void qaws_shell_curve_destroy(qaws_shell_curve* sc)
{
	unsigned int i;
	if (!sc)
		return;
	for (i = 0; i < sc->n; i++) { qaws_surface_curve_destroy(sc->c[i]); qaws_curve_destroy(sc->keep[i]); }
	free(sc->face); free(sc->c); free(sc->keep);
	free(sc);
}

qaws_status qaws_shell_curve_create_geodesic_ray(qaws_shell_atlas const* a, qaws_shell_point p, qaws_vec2 direction,
	qaws_scalar length, qaws_scalar tolerance, qaws_shell_curve** out_curve)
{
	double y[4], d[2];
	unsigned int f, i;
	sh_runs runs;
	qaws_shell_curve* sc;
	qaws_status st;
	if (!out_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_curve = NULL;
	d[0] = direction.x; d[1] = direction.y;
	if ((st = sh_start(a, p, d, y)) != QAWS_STATUS_OK)
		return st;
	memset(&runs, 0, sizeof(runs));
	f = p.face;
	st = sh_trace(a, &f, y, (double)length, &runs, NULL);
	sc = st == QAWS_STATUS_OK ? sh_curve_new(a) : NULL;
	if (st == QAWS_STATUS_OK && !sc)
		st = QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < runs.n && st == QAWS_STATUS_OK; i++)
	{
		qaws_surface_curve* c = NULL;
		qaws_vec2 q, dd;
		if (!(runs.r[i].len > 0))
			continue;
		q.x = (qaws_scalar)runs.r[i].y0[0]; q.y = (qaws_scalar)runs.r[i].y0[1];
		dd.x = (qaws_scalar)runs.r[i].y0[2]; dd.y = (qaws_scalar)runs.r[i].y0[3];
		st = qaws_surface_curve_create_geodesic_ray(a->s[runs.r[i].face], q, dd, (qaws_scalar)runs.r[i].len, tolerance, &c);
		if (st == QAWS_STATUS_OK)
			st = sh_curve_push(sc, runs.r[i].face, c, NULL);
	}
	free(runs.r);
	if (st != QAWS_STATUS_OK)
	{
		qaws_shell_curve_destroy(sc);
		return st;
	}
	*out_curve = sc;
	return QAWS_STATUS_OK;
}

/* the circle's point at angle th: the ray from the centre run r */
typedef struct sh_circ
{
	qaws_shell_atlas const* a;
	qaws_shell_point c;
	double r;
} sh_circ;

static qaws_status sh_circ_at(sh_circ const* k, double th, unsigned int* f, double* uv)
{
	double end[2], vel[2], y[4];
	qaws_status st;
	double cc[2];
	cc[0] = k->c.uv.x; cc[1] = k->c.uv.y;
	if ((st = qaws_internal_sg_shoot(&k->a->c[k->c.face], cc, th, 0, end, vel)) != QAWS_STATUS_OK)
		return st;
	y[0] = end[0]; y[1] = end[1]; y[2] = vel[0]; y[3] = vel[1];
	*f = k->c.face;
	st = sh_trace(k->a, f, y, k->r, NULL, NULL);
	uv[0] = y[0]; uv[1] = y[1];
	return st;
}

/* a run of the circle on face f for angles in [ta, tb], its ends given
   (on the face's sides): a Hermite in (u, v) through points of the
   exponential map, refined until its midpoints agree within the tolerance */
static qaws_status sh_fit_run(sh_circ const* k, unsigned int f, double ta, double tb, double const* ea, double const* eb,
	double tol, qaws_curve** out)
{
	unsigned int n = 16;
	*out = NULL;
	while (n <= 4096)
	{
		double* pts = (double*)malloc(sizeof(double) * 4 * (n + 1));
		double *der, dt = (tb - ta) / n, worst = 0, dl = 1e-6 * (tb - ta);
		unsigned int j, ff;
		qaws_status st = QAWS_STATUS_OK;
		qaws_hermite_desc hd;
		qaws_scalar* buf;
		qaws_curve* cv = NULL;
		if (!pts)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		der = pts + 2 * (n + 1);
		for (j = 0; j <= n && st == QAWS_STATUS_OK; j++)
		{
			double t = ta + dt * j, p[2], qa[2], qb[2], t0, t1;
			if (j == 0) { p[0] = ea[0]; p[1] = ea[1]; }
			else if (j == n) { p[0] = eb[0]; p[1] = eb[1]; }
			else if ((st = sh_circ_at(k, t, &ff, p)) == QAWS_STATUS_OK && ff != f)
				st = QAWS_STATUS_NUMERICAL_FAILURE;
			pts[2 * j] = p[0]; pts[2 * j + 1] = p[1];
			/* d/dj by a difference on the face (one-sided at the ends) */
			t0 = j == 0 ? t : t - dl;
			t1 = j == n ? t : t + dl;
			if (st == QAWS_STATUS_OK && j == 0) { qa[0] = ea[0]; qa[1] = ea[1]; }
			else if (st == QAWS_STATUS_OK && (st = sh_circ_at(k, t0, &ff, qa)) == QAWS_STATUS_OK && ff != f) st = QAWS_STATUS_NUMERICAL_FAILURE;
			if (st == QAWS_STATUS_OK && j == n) { qb[0] = eb[0]; qb[1] = eb[1]; }
			else if (st == QAWS_STATUS_OK && (st = sh_circ_at(k, t1, &ff, qb)) == QAWS_STATUS_OK && ff != f) st = QAWS_STATUS_NUMERICAL_FAILURE;
			if (st == QAWS_STATUS_OK)
			{
				/* the copies across a closed direction: qa, qb near p */
				qaws_internal_sg_ctx const* c = &k->a->c[f];
				if (c->pu > 0) { qa[0] += c->pu * floor((p[0] - qa[0]) / c->pu + 0.5); qb[0] += c->pu * floor((p[0] - qb[0]) / c->pu + 0.5); }
				if (c->pv > 0) { qa[1] += c->pv * floor((p[1] - qa[1]) / c->pv + 0.5); qb[1] += c->pv * floor((p[1] - qb[1]) / c->pv + 0.5); }
				der[2 * j] = (qb[0] - qa[0]) / (t1 - t0) * dt;
				der[2 * j + 1] = (qb[1] - qa[1]) / (t1 - t0) * dt;
			}
		}
		/* the points continuous across a closed direction */
		for (j = 1; j <= n && st == QAWS_STATUS_OK; j++)
		{
			qaws_internal_sg_ctx const* c = &k->a->c[f];
			if (c->pu > 0) pts[2 * j] += c->pu * floor((pts[2 * j - 2] - pts[2 * j]) / c->pu + 0.5);
			if (c->pv > 0) pts[2 * j + 1] += c->pv * floor((pts[2 * j - 1] - pts[2 * j + 1]) / c->pv + 0.5);
		}
		buf = st == QAWS_STATUS_OK ? (qaws_scalar*)malloc(sizeof(qaws_scalar) * 4 * (n + 1)) : NULL;
		if (st == QAWS_STATUS_OK && !buf)
			st = QAWS_STATUS_ALLOCATION_FAILURE;
		if (st == QAWS_STATUS_OK)
		{
			for (j = 0; j < 2 * (n + 1); j++) { buf[j] = (qaws_scalar)pts[j]; buf[2 * (n + 1) + j] = (qaws_scalar)der[j]; }
			memset(&hd, 0, sizeof(hd));
			hd.dimension = QAWS_DIMENSION_2D; hd.degree = 3;
			hd.points = buf; hd.point_count = n + 1;
			hd.derivatives = buf + 2 * (n + 1); hd.derivative_count = n + 1;
			st = qaws_curve_create_hermite(&hd, &cv);
		}
		free(buf);
		for (j = 0; j < n && st == QAWS_STATUS_OK; j++)
		{
			double p[2];
			qaws_eval_result_2d e;
			if ((st = sh_circ_at(k, ta + dt * (j + 0.5), &ff, p)) != QAWS_STATUS_OK)
				break;
			if (ff != f)
			{
				st = QAWS_STATUS_NUMERICAL_FAILURE;
				break;
			}
			if (qaws_curve_evaluate_2d(cv, (qaws_scalar)(j + 0.5), QAWS_EVAL_FLAG_POSITION, &e) == QAWS_STATUS_OK)
			{
				qaws_internal_sg_ctx const* c = &k->a->c[f];
				double du = (double)e.position.x - p[0], dv = (double)e.position.y - p[1];
				if (c->pu > 0) du -= c->pu * floor(du / c->pu + 0.5);
				if (c->pv > 0) dv -= c->pv * floor(dv / c->pv + 0.5);
				if (hypot(du, dv) > worst) worst = hypot(du, dv);
			}
		}
		free(pts);
		if (st != QAWS_STATUS_OK)
		{
			qaws_curve_destroy(cv);
			return st;
		}
		if (worst <= tol || n >= 4096)
		{
			*out = cv;
			return QAWS_STATUS_OK;
		}
		qaws_curve_destroy(cv);
		n *= 2;
	}
	return QAWS_STATUS_NUMERICAL_FAILURE;
}

/* the circle's end on a face's side: the last angle still on face f (by
   bisection from ta on f to tb off it) and that point put on the side */
static qaws_status sh_circ_edge(sh_circ const* k, unsigned int f, double ta, double tb, double* t_out, double* uv_out)
{
	unsigned int it, ff;
	double p[2], lo = ta, hi = tb;
	qaws_internal_sg_ctx const* c = &k->a->c[f];
	qaws_status st = QAWS_STATUS_OK;
	for (it = 0; it < 70; it++)
	{
		double mid = 0.5 * (lo + hi);
		if ((st = sh_circ_at(k, mid, &ff, p)) != QAWS_STATUS_OK)
			return st;
		if (ff == f) lo = mid; else hi = mid;
	}
	if ((st = sh_circ_at(k, lo, &ff, p)) != QAWS_STATUS_OK)
		return st;
	/* onto the nearest side the circle crosses */
	{
		double du0 = fabs(p[0] - c->u0) / (c->u1 - c->u0), du1 = fabs(p[0] - c->u1) / (c->u1 - c->u0);
		double dv0 = fabs(p[1] - c->v0) / (c->v1 - c->v0), dv1 = fabs(p[1] - c->v1) / (c->v1 - c->v0), m = 1e300;
		int side = -1;
		if (c->pv <= 0 && dv0 < m) { m = dv0; side = 0; }
		if (c->pu <= 0 && du1 < m) { m = du1; side = 1; }
		if (c->pv <= 0 && dv1 < m) { m = dv1; side = 2; }
		if (c->pu <= 0 && du0 < m) { m = du0; side = 3; }
		if (side == 0) p[1] = c->v0;
		if (side == 1) p[0] = c->u1;
		if (side == 2) p[1] = c->v1;
		if (side == 3) p[0] = c->u0;
	}
	*t_out = lo;
	uv_out[0] = p[0]; uv_out[1] = p[1];
	return st;
}

qaws_status qaws_shell_curve_create_geodesic_circle(qaws_shell_atlas const* a, qaws_shell_point centre, qaws_scalar radius,
	qaws_scalar tolerance, qaws_shell_curve** out_curve)
{
	sh_circ k;
	unsigned int n = 64, i, *fa = NULL;
	double* th = NULL;
	qaws_shell_curve* sc = NULL;
	qaws_status st = QAWS_STATUS_OK;
	double tol;
	if (!a || !out_curve || centre.face >= a->nf || !(radius > 0))
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_curve = NULL;
	k.a = a; k.c = centre; k.r = (double)radius;
	tol = tolerance > 0 ? (double)tolerance * a->c[centre.face].diag : a->c[centre.face].tol;
	/* the faces the circle visits, at n angles */
	fa = (unsigned int*)malloc(sizeof(unsigned int) * (n + 1));
	th = (double*)malloc(sizeof(double) * (n + 1));
	sc = sh_curve_new(a);
	if (!fa || !th || !sc)
		st = QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < n && st == QAWS_STATUS_OK; i++)
	{
		double p[2];
		th[i] = 2 * SH_PI * i / n;
		st = sh_circ_at(&k, th[i], &fa[i], p);
	}
	if (st == QAWS_STATUS_OK)
	{
		unsigned int first = n;
		for (i = 0; i < n; i++)
			if (fa[i] != fa[(i + n - 1) % n]) { first = i; break; }
		if (first == n)
		{
			/* on one face: the face's own geodesic circle */
			qaws_surface_curve* c = NULL;
			st = qaws_surface_curve_create_geodesic_arc(a->s[centre.face == fa[0] ? centre.face : fa[0]], centre.uv, radius, 0,
				(qaws_scalar)(2 * SH_PI), tolerance, &c);
			if (fa[0] != centre.face && st == QAWS_STATUS_OK)
			{
				qaws_surface_curve_destroy(c);
				st = QAWS_STATUS_UNSUPPORTED_OPERATION;
			}
			else if (st == QAWS_STATUS_OK)
				st = sh_curve_push(sc, centre.face, c, NULL);
		}
		else
		{
			/* the runs, from the first change of face round to it */
			unsigned int start = first, j = 0;
			while (st == QAWS_STATUS_OK && j < n)
			{
				unsigned int f = fa[(start + j) % n], m = 0;
				double ta, tb, ea[2], eb[2], tprev = th[(start + j + n - 1) % n], tcur = th[(start + j) % n];
				qaws_curve* cv = NULL;
				qaws_surface_curve* c = NULL;
				/* the entry: the last angle off f before it, on f's side */
				{
					/* bisection from tcur (on f) back to tprev (off f) */
					if (tprev > tcur) tprev -= 2 * SH_PI;
					st = sh_circ_edge(&k, f, tcur, tprev, &ta, ea);
				}
				while (st == QAWS_STATUS_OK && fa[(start + j + m) % n] == f && m < n) m++;
				if (st != QAWS_STATUS_OK)
					break;
				{
					double tlast = th[(start + j + m - 1) % n], tnext = th[(start + j + m) % n];
					if (tlast < ta) tlast += 2 * SH_PI;
					while (tnext < tlast) tnext += 2 * SH_PI;
					st = sh_circ_edge(&k, f, tlast, tnext, &tb, eb);
					while (tb < ta) tb += 2 * SH_PI;
				}
				if (st == QAWS_STATUS_OK)
				{
					st = sh_fit_run(&k, f, ta, tb, ea, eb, tol, &cv);
				}
				if (st == QAWS_STATUS_OK)
					st = qaws_surface_curve_create_uv(a->s[f], cv, &c);
				if (st == QAWS_STATUS_OK)
					st = sh_curve_push(sc, f, c, cv);
				else
					qaws_curve_destroy(cv);
				j += m;
			}
		}
	}
	free(fa); free(th);
	if (st != QAWS_STATUS_OK)
	{
		qaws_shell_curve_destroy(sc);
		return st;
	}
	sc->closed = 1;
	*out_curve = sc;
	return QAWS_STATUS_OK;
}

unsigned int qaws_shell_curve_get_run_count(qaws_shell_curve const* sc)
{
	return sc ? sc->n : 0;
}

qaws_status qaws_shell_curve_get_run(qaws_shell_curve const* sc, unsigned int index, unsigned int* out_face,
	qaws_surface_curve const** out_curve)
{
	if (!sc || index >= sc->n)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (out_face) *out_face = sc->face[index];
	if (out_curve) *out_curve = sc->c[index];
	return QAWS_STATUS_OK;
}

int qaws_shell_curve_is_closed(qaws_shell_curve const* sc)
{
	return sc ? sc->closed : 0;
}

qaws_status qaws_shell_curve_get_length(qaws_shell_curve const* sc, qaws_scalar* out_length)
{
	unsigned int i;
	double L = 0;
	if (!sc || !out_length)
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < sc->n; i++)
	{
		qaws_scalar l = 0;
		qaws_status st = qaws_surface_curve_get_length(sc->c[i], &l);
		if (st != QAWS_STATUS_OK)
			return st;
		L += (double)l;
	}
	*out_length = (qaws_scalar)L;
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------------ */
/*  Regions                                                                 */
/* ------------------------------------------------------------------------ */

struct qaws_shell_region
{
	qaws_shell_atlas const* a;
	qaws_surface_region** f;        /* per face, NULL: none of it */
};

/* a point on a face's boundary as (circuit, place): an open face has one
   circuit of its four sides (place 0..4, counter-clockwise from (u0, v0));
   a face closed in u has two, its bottom (u growing) and top (u falling),
   each of length 1, and likewise in v */
static void sh_place(qaws_internal_sg_ctx const* c, double const* uv, int* circuit, double* place)
{
	double su = (uv[0] - c->u0) / (c->u1 - c->u0), sv = (uv[1] - c->v0) / (c->v1 - c->v0);
	if (c->pu > 0)
	{
		double w = (uv[0] - c->u0) / c->pu;
		w -= floor(w);
		if (fabs(sv) < fabs(sv - 1)) { *circuit = 0; *place = w; }
		else { *circuit = 1; *place = 1 - w; if (*place >= 1) *place -= 1; }
		return;
	}
	if (c->pv > 0)
	{
		double w = (uv[1] - c->v0) / c->pv;
		w -= floor(w);
		if (fabs(su - 1) < fabs(su)) { *circuit = 0; *place = w; }
		else { *circuit = 1; *place = 1 - w; if (*place >= 1) *place -= 1; }
		return;
	}
	*circuit = 0;
	{
		double d0 = fabs(sv), d1 = fabs(su - 1), d2 = fabs(sv - 1), d3 = fabs(su), m = d0;
		int side = 0;
		if (d1 < m) { m = d1; side = 1; }
		if (d2 < m) { m = d2; side = 2; }
		if (d3 < m) { m = d3; side = 3; }
		if (su < 0) su = 0;
		if (su > 1) su = 1;
		if (sv < 0) sv = 0;
		if (sv > 1) sv = 1;
		*place = side == 0 ? su : side == 1 ? 1 + sv : side == 2 ? 3 - su : 4 - sv;
		if (*place >= 4) *place -= 4;
	}
}

/* the (u, v) of a place on an open face's circuit */
static void sh_place_uv(qaws_internal_sg_ctx const* c, double place, double* uv)
{
	if (place < 1) { uv[0] = c->u0 + place * (c->u1 - c->u0); uv[1] = c->v0; }
	else if (place < 2) { uv[0] = c->u1; uv[1] = c->v0 + (place - 1) * (c->v1 - c->v0); }
	else if (place < 3) { uv[0] = c->u1 - (place - 2) * (c->u1 - c->u0); uv[1] = c->v1; }
	else { uv[0] = c->u0; uv[1] = c->v1 - (place - 3) * (c->v1 - c->v0); }
}

/* the curves made for one face's loops */
typedef struct sh_bag
{
	qaws_curve** c;
	unsigned int n, cap;
	qaws_surface_curve** sc;
	unsigned int ns, caps;
} sh_bag;

static qaws_status sh_bag_line(sh_bag* b, qaws_surface const* s, double const* p, double const* q, qaws_surface_curve const** out)
{
	qaws_bezier_desc bd;
	qaws_scalar cp[4];
	qaws_curve* k = NULL;
	qaws_surface_curve* sc = NULL;
	qaws_status st;
	cp[0] = (qaws_scalar)p[0]; cp[1] = (qaws_scalar)p[1]; cp[2] = (qaws_scalar)q[0]; cp[3] = (qaws_scalar)q[1];
	memset(&bd, 0, sizeof(bd));
	bd.dimension = QAWS_DIMENSION_2D; bd.degree = 1; bd.control_points = cp; bd.control_point_count = 2;
	if ((st = qaws_curve_create_bezier(&bd, &k)) != QAWS_STATUS_OK)
		return st;
	st = qaws_surface_curve_create_uv(s, k, &sc);
	qaws_curve_destroy(k);
	if (st != QAWS_STATUS_OK)
		return st;
	if (b->ns + 1 > b->caps)
	{
		unsigned int cap = b->caps ? 2 * b->caps : 16;
		qaws_surface_curve** q2 = (qaws_surface_curve**)realloc(b->sc, sizeof(qaws_surface_curve*) * cap);
		if (!q2) { qaws_surface_curve_destroy(sc); return QAWS_STATUS_ALLOCATION_FAILURE; }
		b->sc = q2;
		b->caps = cap;
	}
	b->sc[b->ns++] = sc;
	*out = sc;
	return QAWS_STATUS_OK;
}

/* a run moved by whole periods (a new surface curve) */
static qaws_status sh_bag_moved(sh_bag* b, qaws_surface const* s, qaws_surface_curve const* run, double du, double dv,
	qaws_surface_curve const** out)
{
	qaws_scalar t[6];
	qaws_curve* mv = NULL;
	qaws_surface_curve* sc = NULL;
	qaws_status st;
	if (du == 0 && dv == 0)
	{
		*out = run;
		return QAWS_STATUS_OK;
	}
	t[0] = 1; t[1] = 0; t[2] = (qaws_scalar)du; t[3] = 0; t[4] = 1; t[5] = (qaws_scalar)dv;
	if ((st = qaws_curve_transform_2d(qaws_surface_curve_get_uv(run), t, &mv)) != QAWS_STATUS_OK)
		return st;
	st = qaws_surface_curve_create_uv(s, mv, &sc);
	qaws_curve_destroy(mv);
	if (st != QAWS_STATUS_OK)
		return st;
	if (b->ns + 1 > b->caps)
	{
		unsigned int cap = b->caps ? 2 * b->caps : 16;
		qaws_surface_curve** q2 = (qaws_surface_curve**)realloc(b->sc, sizeof(qaws_surface_curve*) * cap);
		if (!q2) { qaws_surface_curve_destroy(sc); return QAWS_STATUS_ALLOCATION_FAILURE; }
		b->sc = q2;
		b->caps = cap;
	}
	b->sc[b->ns++] = sc;
	*out = sc;
	return QAWS_STATUS_OK;
}

static void sh_bag_free(sh_bag* b)
{
	unsigned int i;
	for (i = 0; i < b->ns; i++) qaws_surface_curve_destroy(b->sc[i]);
	for (i = 0; i < b->n; i++) qaws_curve_destroy(b->c[i]);
	free(b->sc); free(b->c);
	memset(b, 0, sizeof(*b));
}

static void sh_run_ends(qaws_surface_curve const* run, double* p0, double* p1)
{
	qaws_curve const* uv = qaws_surface_curve_get_uv(run);
	qaws_range rg = qaws_curve_get_parameter_range(uv);
	qaws_eval_result_2d e;
	qaws_curve_evaluate_2d(uv, rg.min_value, QAWS_EVAL_FLAG_POSITION, &e);
	p0[0] = e.position.x; p0[1] = e.position.y;
	qaws_curve_evaluate_2d(uv, rg.max_value, QAWS_EVAL_FLAG_POSITION, &e);
	p1[0] = e.position.x; p1[1] = e.position.y;
}

/* the whole face as a region */
static qaws_status sh_full_face(qaws_shell_atlas const* a, unsigned int f, qaws_surface_region** out)
{
	qaws_internal_sg_ctx const* c = &a->c[f];
	sh_bag b;
	qaws_surface_curve const* cv[4];
	qaws_surface_path P;
	qaws_status st = QAWS_STATUS_OK;
	unsigned int nc = 0, k;
	memset(&b, 0, sizeof(b));
	if (c->pu > 0 && c->pv > 0)
	{
		/* a torus: below v0 + pv round u, and its constant */
		double p[2], q[2];
		p[0] = c->u0; p[1] = c->v0; q[0] = c->u0 + c->pu; q[1] = c->v0;
		st = sh_bag_line(&b, a->s[f], p, q, &cv[nc++]);
	}
	else if (c->pu > 0)
	{
		double p[2], q[2];
		p[0] = c->u0; p[1] = c->v0; q[0] = c->u0 + c->pu; q[1] = c->v0;
		st = sh_bag_line(&b, a->s[f], p, q, &cv[nc++]);
	}
	else if (c->pv > 0)
	{
		double p[2], q[2];
		p[0] = c->u1; p[1] = c->v0; q[0] = c->u1; q[1] = c->v0 + c->pv;
		st = sh_bag_line(&b, a->s[f], p, q, &cv[nc++]);
	}
	else
		for (k = 0; k < 4 && st == QAWS_STATUS_OK; k++)
		{
			double p[2], q[2];
			sh_place_uv(c, k, p);
			sh_place_uv(c, k + 1 < 4 ? k + 1 : 0, q);
			if (k == 3) { q[0] = c->u0; q[1] = c->v0; }
			st = sh_bag_line(&b, a->s[f], p, q, &cv[nc++]);
		}
	P.curves = cv; P.curve_count = nc; P.closed = 1;
	if (st == QAWS_STATUS_OK)
		st = qaws_surface_clip(a->s[f], QAWS_CLIP_NONE, QAWS_FILL_NON_ZERO, &P, 1, NULL, 0, out);
	sh_bag_free(&b);
	return st;
}

/* one face's region from the runs on it: closed runs as they are, open
   runs (side to side) joined along the face's boundary, counter-clockwise
   from where one leaves to where the next comes in */
static qaws_status sh_face_region(qaws_shell_atlas const* a, unsigned int f, qaws_surface_curve const* const* runs, int const* open,
	unsigned int nr, qaws_surface_region** out)
{
	qaws_internal_sg_ctx const* c = &a->c[f];
	qaws_surface_path* paths = (qaws_surface_path*)malloc(sizeof(qaws_surface_path) * (nr + 1));
	qaws_surface_curve const** all = (qaws_surface_curve const**)malloc(sizeof(qaws_surface_curve const*) * (8 * nr + 8));
	int* used = (int*)calloc(nr + 1, sizeof(int));
	unsigned int np = 0, na = 0, i;
	sh_bag b;
	qaws_status st = QAWS_STATUS_OK;
	memset(&b, 0, sizeof(b));
	if (!paths || !all || !used)
		st = QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < nr && st == QAWS_STATUS_OK; i++)
		if (!open[i])
		{
			all[na] = runs[i];
			paths[np].curves = all + na; paths[np].curve_count = 1; paths[np].closed = 1;
			na++; np++;
			used[i] = 1;
		}
	for (i = 0; i < nr && st == QAWS_STATUS_OK; i++)
	{
		unsigned int cur = i, first = na, guard = 0;
		double shift[2] = { 0, 0 };
		if (used[i])
			continue;
		while (st == QAWS_STATUS_OK && guard++ <= nr)
		{
			double p0[2], p1[2], pe[2];
			int ce;
			double we, best = 1e300;
			unsigned int nxt = SH_NONE, k;
			qaws_surface_curve const* mv;
			used[cur] = 1;
			if ((st = sh_bag_moved(&b, a->s[f], runs[cur], shift[0], shift[1], &mv)) != QAWS_STATUS_OK)
				break;
			all[na++] = mv;
			sh_run_ends(mv, p0, p1);
			pe[0] = p1[0]; pe[1] = p1[1];
			sh_place(c, pe, &ce, &we);
			/* the next entry along the circuit, counter-clockwise */
			for (k = 0; k < nr; k++)
			{
				double q0[2], q1[2], wk, len = (c->pu > 0 || c->pv > 0) ? 1 : 4, d;
				int ck;
				if (!open[k] || (used[k] && k != i))
					continue;
				sh_run_ends(runs[k], q0, q1);
				sh_place(c, q0, &ck, &wk);
				if (ck != ce)
					continue;
				d = wk - we;
				if (d < 0) d += len;
				if (d < best) { best = d; nxt = k; }
			}
			if (nxt == SH_NONE)
			{
				st = QAWS_STATUS_NUMERICAL_FAILURE;
				break;
			}
			/* along the boundary from pe to the next run's start */
			{
				double q0[2], q1[2], qn[2];
				sh_run_ends(runs[nxt], q0, q1);
				if (c->pu > 0 || c->pv > 0)
				{
					/* along a closed direction's side: the next start's copy past pe
					   that way */
					/* circuit 0 runs the closed direction growing, circuit 1 falling */
					int dir = ce == 0 ? 1 : -1;
					double P = c->pu > 0 ? c->pu : c->pv;
					unsigned int ax = c->pu > 0 ? 0 : 1;
					double target;
					if (dir > 0) target = q0[ax] + P * ceil((pe[ax] - q0[ax]) / P - 1e-12);
					else target = q0[ax] + P * floor((pe[ax] - q0[ax]) / P + 1e-12);
					qn[ax] = target;
					qn[1 - ax] = pe[1 - ax];
					if (hypot(qn[0] - pe[0], qn[1] - pe[1]) > 0)
						st = sh_bag_line(&b, a->s[f], pe, qn, &all[na++]);
					shift[0] = ax == 0 ? qn[0] - q0[0] : 0;
					shift[1] = ax == 1 ? qn[1] - q0[1] : 0;
				}
				else
				{
					double wn, at = we, from[2];
					int cn;
					sh_place(c, q0, &cn, &wn);
					if (wn < we) wn += 4;
					from[0] = pe[0]; from[1] = pe[1];
					/* the corners between */
					while (st == QAWS_STATUS_OK && floor(at) + 1 <= wn - 1e-12)
					{
						double cpt[2];
						at = floor(at) + 1;
						sh_place_uv(c, at >= 4 ? at - 4 : at, cpt);
						st = sh_bag_line(&b, a->s[f], from, cpt, &all[na++]);
						from[0] = cpt[0]; from[1] = cpt[1];
					}
					if (st == QAWS_STATUS_OK && hypot(q0[0] - from[0], q0[1] - from[1]) > 0)
						st = sh_bag_line(&b, a->s[f], from, q0, &all[na++]);
				}
			}
			if (nxt == i)
				break;
			cur = nxt;
		}
		if (st == QAWS_STATUS_OK)
		{
			paths[np].curves = all + first; paths[np].curve_count = na - first; paths[np].closed = 1;
			np++;
		}
	}
	if (st == QAWS_STATUS_OK)
		st = qaws_surface_clip(a->s[f], QAWS_CLIP_NONE, QAWS_FILL_NON_ZERO, paths, np, NULL, 0, out);
	sh_bag_free(&b);
	free(paths); free(all); free(used);
	return st;
}

qaws_status qaws_shell_region_create(qaws_shell_atlas const* a, qaws_shell_curve const* const* loops, unsigned int nl,
	qaws_shell_region** out_region)
{
	qaws_shell_region* r;
	unsigned int f, i, j, *state;
	qaws_status st = QAWS_STATUS_OK;
	if (!a || !out_region || (!loops && nl))
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_region = NULL;
	r = (qaws_shell_region*)calloc(1, sizeof(*r));
	state = (unsigned int*)calloc(a->nf + 1, sizeof(unsigned int));   /* 0 unknown, 1 runs, 2 inside, 3 outside */
	if (!r || !state)
	{
		free(r); free(state);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	r->a = a;
	r->f = (qaws_surface_region**)calloc(a->nf + 1, sizeof(qaws_surface_region*));
	if (!r->f)
		st = QAWS_STATUS_ALLOCATION_FAILURE;
	/* the faces with runs */
	for (f = 0; f < a->nf && st == QAWS_STATUS_OK; f++)
	{
		qaws_surface_curve const** runs;
		int* open;
		unsigned int nr = 0;
		for (i = 0; i < nl; i++)
			for (j = 0; j < loops[i]->n; j++)
				if (loops[i]->face[j] == f) nr++;
		if (!nr)
			continue;
		runs = (qaws_surface_curve const**)malloc(sizeof(qaws_surface_curve const*) * nr);
		open = (int*)malloc(sizeof(int) * nr);
		if (!runs || !open)
		{
			free(runs); free(open);
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			break;
		}
		nr = 0;
		for (i = 0; i < nl; i++)
			for (j = 0; j < loops[i]->n; j++)
				if (loops[i]->face[j] == f)
				{
					runs[nr] = loops[i]->c[j];
					open[nr] = loops[i]->n > 1;
					nr++;
				}
		st = sh_face_region(a, f, runs, open, nr, &r->f[f]);
		state[f] = 1;
		free(runs); free(open);
	}
	/* the faces without: inside or outside as a neighbour's side is */
	{
		int changed = 1;
		while (changed && st == QAWS_STATUS_OK)
		{
			changed = 0;
			for (f = 0; f < a->nf && st == QAWS_STATUS_OK; f++)
			{
				int s2;
				if (state[f])
					continue;
				for (s2 = 0; s2 < 4 && !state[f]; s2++)
				{
					unsigned int g = a->nb[f][s2];
					int sg2 = a->nbside[f][s2];
					if (g == SH_NONE || !state[g])
						continue;
					if (state[g] == 2 || state[g] == 3)
						state[f] = state[g];
					else
					{
						/* the neighbour's region just inside its side's middle */
						double uv[2];
						int in = 0;
						qaws_vec2 q;
						qaws_internal_sg_ctx const* cg = &a->c[g];
						sh_side_uv(cg, sg2, 0.5, uv);
						if (sg2 == 0) uv[1] += 1e-6 * (cg->v1 - cg->v0);
						if (sg2 == 1) uv[0] -= 1e-6 * (cg->u1 - cg->u0);
						if (sg2 == 2) uv[1] -= 1e-6 * (cg->v1 - cg->v0);
						if (sg2 == 3) uv[0] += 1e-6 * (cg->u1 - cg->u0);
						q.x = (qaws_scalar)uv[0]; q.y = (qaws_scalar)uv[1];
						qaws_surface_region_contains(r->f[g], q, &in);
						state[f] = in == 1 ? 2 : 3;
					}
					if (state[f] == 2)
						st = sh_full_face(a, f, &r->f[f]);
					changed = 1;
				}
			}
		}
	}
	free(state);
	if (st != QAWS_STATUS_OK)
	{
		qaws_shell_region_destroy(r);
		return st;
	}
	*out_region = r;
	return QAWS_STATUS_OK;
}

void qaws_shell_region_destroy(qaws_shell_region* r)
{
	unsigned int f;
	if (!r)
		return;
	if (r->f)
		for (f = 0; f < r->a->nf; f++) qaws_surface_region_destroy(r->f[f]);
	free(r->f);
	free(r);
}

qaws_surface_region const* qaws_shell_region_get_face(qaws_shell_region const* r, unsigned int face)
{
	return r && face < r->a->nf ? r->f[face] : NULL;
}

qaws_status qaws_shell_region_get_area(qaws_shell_region const* r, qaws_scalar* out_area)
{
	unsigned int f;
	double A = 0;
	if (!r || !out_area)
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (f = 0; f < r->a->nf; f++)
		if (r->f[f])
		{
			qaws_scalar x = 0;
			qaws_status st = qaws_surface_region_get_area(r->f[f], &x);
			if (st != QAWS_STATUS_OK)
				return st;
			A += (double)x;
		}
	*out_area = (qaws_scalar)A;
	return QAWS_STATUS_OK;
}

/* a stretch of a face's loop along a shared side (inside the region on
   both faces, not its boundary): its surface length */
static double sh_shared_length(qaws_shell_atlas const* a, unsigned int f, qaws_curve const* cv)
{
	qaws_internal_sg_ctx const* c = &a->c[f];
	unsigned int n = 0, i, s;
	qaws_scalar probe[2], *cp;
	double L = 0;
	if (qaws_curve_get_degree(cv) != 1)
		return 0;
	qaws_curve_get_control_points(cv, probe, 0, &n);
	cp = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * (n + 1));
	if (!cp)
		return 0;
	qaws_curve_get_control_points(cv, cp, n, &n);
	for (i = 0; i + 1 < n; i++)
		for (s = 0; s < 4; s++)
		{
			double p[2], q[2], eu = SH_ON_SIDE * (c->u1 - c->u0), ev = SH_ON_SIDE * (c->v1 - c->v0);
			int on;
			if (a->nb[f][s] == SH_NONE)
				continue;
			p[0] = cp[2 * i]; p[1] = cp[2 * i + 1]; q[0] = cp[2 * i + 2]; q[1] = cp[2 * i + 3];
			if (s == 0) on = fabs(p[1] - c->v0) < ev && fabs(q[1] - c->v0) < ev;
			else if (s == 1) on = fabs(p[0] - c->u1) < eu && fabs(q[0] - c->u1) < eu;
			else if (s == 2) on = fabs(p[1] - c->v1) < ev && fabs(q[1] - c->v1) < ev;
			else on = fabs(p[0] - c->u0) < eu && fabs(q[0] - c->u0) < eu;
			if (on)
			{
				/* its length by Gauss on the segment */
				static double const x8[4] = { 0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363 };
				static double const w8[4] = { 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763 };
				unsigned int k, sg;
				for (k = 0; k < 4; k++)
					for (sg = 0; sg < 2; sg++)
					{
						double t = 0.5 + (sg ? 0.5 : -0.5) * x8[k], g[3], d[2];
						d[0] = q[0] - p[0]; d[1] = q[1] - p[1];
						if (qaws_internal_sg_first_form(c, p[0] + t * d[0], p[1] + t * d[1], g) == QAWS_STATUS_OK)
							L += 0.5 * w8[k] * sqrt(g[0] * d[0] * d[0] + 2 * g[1] * d[0] * d[1] + g[2] * d[1] * d[1]);
					}
				break;
			}
		}
	free(cp);
	return L;
}

qaws_status qaws_shell_region_get_perimeter(qaws_shell_region const* r, qaws_scalar* out_length)
{
	unsigned int f, i, j;
	double L = 0;
	if (!r || !out_length)
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (f = 0; f < r->a->nf; f++)
		if (r->f[f])
		{
			qaws_scalar x = 0;
			qaws_status st = qaws_surface_region_get_perimeter(r->f[f], &x);
			if (st != QAWS_STATUS_OK)
				return st;
			L += (double)x;
			for (i = 0; i < qaws_surface_region_get_path_count(r->f[f]); i++)
			{
				qaws_path_2d p;
				qaws_surface_region_get_path(r->f[f], i, &p);
				for (j = 0; j < p.curve_count; j++)
					L -= sh_shared_length(r->a, f, p.curves[j]);
			}
		}
	*out_length = (qaws_scalar)L;
	return QAWS_STATUS_OK;
}

qaws_status qaws_shell_region_contains(qaws_shell_region const* r, qaws_shell_point p, int* out_inside)
{
	if (!r || !out_inside || p.face >= r->a->nf)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_inside = 0;
	if (!r->f[p.face])
		return QAWS_STATUS_OK;
	return qaws_surface_region_contains(r->f[p.face], p.uv, out_inside);
}

/* a face region's loops as surface paths (their curves wrapped) */
static qaws_status sh_region_paths(qaws_surface const* s, qaws_surface_region const* r, sh_bag* b, qaws_surface_path** out, unsigned int* n,
	qaws_surface_curve const*** curves)
{
	unsigned int np = qaws_surface_region_get_path_count(r), i, j, total = 0, at = 0;
	qaws_status st = QAWS_STATUS_OK;
	for (i = 0; i < np; i++)
	{
		qaws_path_2d p;
		qaws_surface_region_get_path(r, i, &p);
		total += p.curve_count;
	}
	*out = (qaws_surface_path*)malloc(sizeof(qaws_surface_path) * (np + 1));
	*curves = (qaws_surface_curve const**)malloc(sizeof(qaws_surface_curve const*) * (total + 1));
	*n = np;
	if (!*out || !*curves)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < np && st == QAWS_STATUS_OK; i++)
	{
		qaws_path_2d p;
		qaws_surface_region_get_path(r, i, &p);
		(*out)[i].curves = *curves + at; (*out)[i].curve_count = p.curve_count; (*out)[i].closed = 1;
		for (j = 0; j < p.curve_count && st == QAWS_STATUS_OK; j++)
		{
			qaws_surface_curve* sc = NULL;
			st = qaws_surface_curve_create_uv(s, p.curves[j], &sc);
			if (st == QAWS_STATUS_OK)
			{
				if (b->ns + 1 > b->caps)
				{
					unsigned int cap = b->caps ? 2 * b->caps : 16;
					qaws_surface_curve** q2 = (qaws_surface_curve**)realloc(b->sc, sizeof(qaws_surface_curve*) * cap);
					if (!q2) { qaws_surface_curve_destroy(sc); return QAWS_STATUS_ALLOCATION_FAILURE; }
					b->sc = q2;
					b->caps = cap;
				}
				b->sc[b->ns++] = sc;
				(*curves)[at++] = sc;
			}
		}
	}
	return st;
}

qaws_status qaws_shell_region_clip(qaws_shell_region const* x, qaws_shell_region const* y, qaws_clip_type type,
	qaws_shell_region** out_region)
{
	qaws_shell_region* r;
	unsigned int f;
	qaws_status st = QAWS_STATUS_OK;
	if (!x || !y || !out_region || x->a != y->a)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_region = NULL;
	r = (qaws_shell_region*)calloc(1, sizeof(*r));
	if (!r)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	r->a = x->a;
	r->f = (qaws_surface_region**)calloc(x->a->nf + 1, sizeof(qaws_surface_region*));
	if (!r->f)
		st = QAWS_STATUS_ALLOCATION_FAILURE;
	for (f = 0; f < x->a->nf && st == QAWS_STATUS_OK; f++)
	{
		qaws_surface_region const *p = x->f[f], *q = y->f[f];
		sh_bag b;
		qaws_surface_path *pa = NULL, *pb = NULL;
		qaws_surface_curve const **ca = NULL, **cb = NULL;
		unsigned int na = 0, nb = 0;
		if (!p && !q)
			continue;
		if ((!p || !q) && type == QAWS_CLIP_INTERSECTION)
			continue;
		if (!p && type == QAWS_CLIP_DIFFERENCE)
			continue;
		memset(&b, 0, sizeof(b));
		if (p) st = sh_region_paths(x->a->s[f], p, &b, &pa, &na, &ca);
		if (st == QAWS_STATUS_OK && q && !(type == QAWS_CLIP_DIFFERENCE && !q)) st = sh_region_paths(x->a->s[f], q, &b, &pb, &nb, &cb);
		if (st == QAWS_STATUS_OK)
		{
			if (!p)
				st = qaws_surface_clip(x->a->s[f], QAWS_CLIP_NONE, QAWS_FILL_NON_ZERO, pb, nb, NULL, 0, &r->f[f]);
			else if (!q)
				st = qaws_surface_clip(x->a->s[f], QAWS_CLIP_NONE, QAWS_FILL_NON_ZERO, pa, na, NULL, 0, &r->f[f]);
			else
				st = qaws_surface_clip(x->a->s[f], type, QAWS_FILL_NON_ZERO, pa, na, pb, nb, &r->f[f]);
		}
		/* nothing left on the face: none of it */
		if (st == QAWS_STATUS_OK && r->f[f] && qaws_surface_region_get_path_count(r->f[f]) == 0)
		{
			qaws_surface_region_destroy(r->f[f]);
			r->f[f] = NULL;
		}
		free(pa); free(pb); free(ca); free(cb);
		sh_bag_free(&b);
	}
	if (st != QAWS_STATUS_OK)
	{
		qaws_shell_region_destroy(r);
		return st;
	}
	*out_region = r;
	return QAWS_STATUS_OK;
}
