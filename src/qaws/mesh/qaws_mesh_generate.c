/*
 * Meshes of exact surfaces and shells: a dense conformal mesh of the whole
 * shell (its faces' shared edges and free sides as features), remeshed or
 * quadrangulated on its PN reference, then every vertex snapped onto its
 * face's surface (Newton from the dense mesh's (u, v)). STRUCTURED:
 * iso-parametric grids per face (mg_structured).
 */

#include "qaws_mesh_internal.h"
#include "../internal/qaws_internal_memory.h"
#include "../qaws_mesh_generate.h"
#include "../qaws_remesh.h"
#include "../qaws_quadrangulate.h"
#include "../qaws_surface.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

static double mg_dot(double const* a, double const* b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void mg_v3(qaws_vec3 v, double* out)
{
	out[0] = (double)v.x;
	out[1] = (double)v.y;
	out[2] = (double)v.z;
}

/* the closest point of s to p by Gauss-Newton from (u, v), kept in the domain; fixed[0], fixed[1]: u, v held */
static void mg_newton(qaws_surface const* s, double const* p, double* u, double* v, int const* fixed)
{
	qaws_range ur = qaws_surface_get_u_range(s), vr = qaws_surface_get_v_range(s);
	double lo[2] = { (double)ur.min_value, (double)vr.min_value }, hi[2] = { (double)ur.max_value, (double)vr.max_value };
	int it;
	for (it = 0; it < 30 && !(fixed[0] && fixed[1]); it++)
	{
		qaws_surface_eval_result e;
		double S[3], Su[3], Sv[3], r[3], a, b, c, g0, g1, det, du, dv;
		unsigned int k;
		memset(&e, 0, sizeof e);
		if (qaws_surface_evaluate(&e, s, (qaws_scalar)*u, (qaws_scalar)*v, QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV) != QAWS_STATUS_OK)
			break;
		mg_v3(e.position, S);
		mg_v3(e.du, Su);
		mg_v3(e.dv, Sv);
		for (k = 0; k < 3; k++)
			r[k] = S[k] - p[k];
		a = mg_dot(Su, Su); b = mg_dot(Su, Sv); c = mg_dot(Sv, Sv);
		g0 = mg_dot(Su, r); g1 = mg_dot(Sv, r);
		det = fixed[0] ? c : fixed[1] ? a : a * c - b * b;
		if (!(det > 1e-300))
			break;
		du = fixed[0] ? 0 : fixed[1] ? -g0 / a : -(c * g0 - b * g1) / det;
		dv = fixed[1] ? 0 : fixed[0] ? -g1 / c : -(a * g1 - b * g0) / det;
		*u += du;
		*v += dv;
		*u = *u < lo[0] ? lo[0] : *u > hi[0] ? hi[0] : *u;
		*v = *v < lo[1] ? lo[1] : *v > hi[1] ? hi[1] : *v;
		if (fabs(du) <= 1e-14 * (hi[0] - lo[0] + 1) && fabs(dv) <= 1e-14 * (hi[1] - lo[1] + 1))
			break;
	}
}

/* p onto s from (u, v); `on`: 1 a vertex of the border or a crease, held on the side of the domain nearest
   its free snap (within a hundredth of the domain), 2 a corner, held at the nearest corner of the domain
   shortcut: a trimmed face's vertices on its trim curve are snapped onto the surface, not the curve (off it
   by the dense mesh's error), the trim curve to be followed when trimmed shells need it */
static void mg_snap(qaws_surface const* s, double const* p, double* u, double* v, int on, double* x)
{
	qaws_range ur = qaws_surface_get_u_range(s), vr = qaws_surface_get_v_range(s);
	double lo[2] = { (double)ur.min_value, (double)vr.min_value }, hi[2] = { (double)ur.max_value, (double)vr.max_value }, d[2];
	double* uv[2];
	int fixed[2] = { 0, 0 }, low[2], j;
	uv[0] = u;
	uv[1] = v;
	mg_newton(s, p, u, v, fixed);
	for (j = 0; j < 2; j++)
	{
		double a = (*uv[j] - lo[j]) / (hi[j] - lo[j]), b = (hi[j] - *uv[j]) / (hi[j] - lo[j]);
		d[j] = fmin(a, b);
		low[j] = a <= b;
	}
	for (j = 0; j < 2; j++)
		fixed[j] = on && d[j] <= 0.01 && (on == 2 || d[j] <= d[1 - j]);
	for (j = 0; j < 2; j++)
		if (fixed[j])
			*uv[j] = low[j] ? lo[j] : hi[j];
	if (fixed[0] || fixed[1])
		mg_newton(s, p, u, v, fixed);
	{
		qaws_surface_eval_result e;
		memset(&e, 0, sizeof e);
		qaws_surface_evaluate(&e, s, (qaws_scalar)*u, (qaws_scalar)*v, QAWS_SURFACE_EVAL_POSITION);
		mg_v3(e.position, x);
	}
}

/* a start for the snap onto face F from n corners tv weighted by w (those on F, their u and v
   unwrapped across a seam); none on F: the nearest point of a coarse grid of F's domain */
static void mg_start(qaws_surface const* s, qaws_shell_point const* dv, unsigned int const* tv, double const* w, unsigned int n_corner, unsigned int F,
	double const* p, double* u, double* v)
{
	qaws_range ur = qaws_surface_get_u_range(s), vr = qaws_surface_get_v_range(s);
	double su = (double)ur.max_value - (double)ur.min_value, sv = (double)vr.max_value - (double)vr.min_value, wu = 0, ref_u = 0, ref_v = 0, acc_u = 0, acc_v = 0;
	unsigned int k, n = 0;
	for (k = 0; k < n_corner; k++)
	{
		double cu, cv;
		if (dv[tv[k]].face != F)
			continue;
		cu = (double)dv[tv[k]].uv.x;
		cv = (double)dv[tv[k]].uv.y;
		if (n == 0)
		{
			ref_u = cu;
			ref_v = cv;
		}
		/* across a seam: within half a period of the first */
		if (cu - ref_u > 0.5 * su) cu -= su;
		if (ref_u - cu > 0.5 * su) cu += su;
		if (cv - ref_v > 0.5 * sv) cv -= sv;
		if (ref_v - cv > 0.5 * sv) cv += sv;
		acc_u += w[k] * cu;
		acc_v += w[k] * cv;
		wu += w[k];
		n++;
	}
	if (n && wu > 0)
	{
		*u = acc_u / wu;
		*v = acc_v / wu;
		if (*u > (double)ur.max_value) *u -= su;
		if (*u < (double)ur.min_value) *u += su;
		if (*v > (double)vr.max_value) *v -= sv;
		if (*v < (double)vr.min_value) *v += sv;
		return;
	}
	{
		double best = HUGE_VAL;
		unsigned int i, j;
		for (i = 0; i <= 8; i++)
			for (j = 0; j <= 8; j++)
			{
				qaws_surface_eval_result e;
				double q[3], d, uu = (double)ur.min_value + su * i / 8, vv = (double)vr.min_value + sv * j / 8;
				memset(&e, 0, sizeof e);
				qaws_surface_evaluate(&e, s, (qaws_scalar)uu, (qaws_scalar)vv, QAWS_SURFACE_EVAL_POSITION);
				mg_v3(e.position, q);
				d = (q[0] - p[0]) * (q[0] - p[0]) + (q[1] - p[1]) * (q[1] - p[1]) + (q[2] - p[2]) * (q[2] - p[2]);
				if (d < best)
				{
					best = d;
					*u = uu;
					*v = vv;
				}
			}
	}
}

/* the faces' surfaces: a shell's, or one surface's */
typedef struct mg_src
{
	qaws_shell_atlas const* shell;
	qaws_surface const* one;
} mg_src;

static qaws_surface const* mg_face(mg_src const* s, unsigned int f)
{
	return s->shell ? qaws_shell_atlas_get_face(s->shell, f) : s->one;
}

/* The common part: the dense mesh (positions x, its vertices' points dvp,
   triangles tri on faces tface) as a reference with the faces' boundaries
   as features, the elements on it, every vertex snapped onto its face. */
static qaws_status mg_finish(mg_src const* src, double const* x, qaws_shell_point const* dvp, unsigned int nv, unsigned int const* tri,
	unsigned int const* tface, unsigned int nt, qaws_generate_desc const* desc, double L, qaws_mesh_3d** out_mesh, qaws_mesh_report* rep)
{
	qaws_mesh_3d *dense = NULL, *out = NULL;
	qaws_mesh_reference* ref = NULL;
	qaws_size_field* sf = NULL;
	unsigned int i, k, fa = 0;
	double dmax = 0, dsum = 0;
	qaws_status st = hm_mesh_3d_create(x, nv, tri, NULL, nt, &dense);
	if (st != QAWS_STATUS_OK)
		return st;
	{
		/* each triangle's face as an attribute; the edges between faces as features */
		qaws_attribute_desc ad;
		qaws_scalar* fv = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(qaws_scalar), (nt + 1)));
		qaws_hm* m = &dense->hm;
		unsigned int c;
		if (!fv)
		{
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			goto done;
		}
		for (i = 0; i < nt; i++)
			fv[i] = (qaws_scalar)tface[i];
		memset(&ad, 0, sizeof ad);
		ad.name = "face";
		ad.element = QAWS_MESH_FACE;
		ad.components = 1;
		ad.rule = QAWS_ATTRIBUTE_NEAREST;
		st = qaws_mesh_3d_add_attribute(&fa, dense, fv, &ad);
		qaws_internal_free(fv);
		if (st != QAWS_STATUS_OK)
			goto done;
		for (c = 0; c < m->nc; c++)
		{
			unsigned int o = m->co[c];
			if (m->cf[c] == HM_NONE || o == HM_NONE)
				continue;
			if (hm_attr_at(m, fa, m->cf[c])[0] != hm_attr_at(m, fa, m->cf[o])[0])
			{
				m->cflag[c] |= HM_FEATURE;
				m->cflag[o] |= HM_FEATURE;
			}
		}
		for (i = 0; i < m->nv; i++)
		{
			unsigned int fd = hm_feature_degree(m, i);
			if (fd != 0 && fd != 2)
				m->vflag[i] |= HM_CORNER;
		}
	}
	st = qaws_mesh_reference_create(&ref, dense, NULL);
	if (st != QAWS_STATUS_OK)
		goto done;
	if (desc->target_length == 0 && (desc->chord_tolerance > 0 || desc->angle_tolerance > 0))
	{
		qaws_curvature_size_desc cd;
		memset(&cd, 0, sizeof cd);
		cd.chord_tolerance = desc->chord_tolerance;
		cd.angle_tolerance = desc->angle_tolerance;
		cd.min_length = desc->min_length;
		cd.max_length = desc->max_length > 0 ? desc->max_length : (qaws_scalar)L;
		cd.gradation = desc->gradation;
		st = qaws_size_field_create_curvature(&sf, ref, &cd);
		if (st != QAWS_STATUS_OK)
			goto done;
	}
	if (desc->elements == QAWS_ELEMENTS_TRIANGLES)
	{
		qaws_remesh_desc rd;
		memset(&rd, 0, sizeof rd);
		rd.reference = ref;
		rd.target_length = (qaws_scalar)L;
		rd.size = sf;
		st = qaws_mesh_remesh(&out, rep, dense, &rd);
	}
	else
	{
		qaws_quadrangulate_desc qd;
		memset(&qd, 0, sizeof qd);
		qd.reference = ref;
		qd.target_length = (qaws_scalar)L;
		qd.pure = desc->elements == QAWS_ELEMENTS_QUADS;
		st = qaws_mesh_quadrangulate(&out, rep, dense, &qd);
	}
	if (st != QAWS_STATUS_OK)
		goto done;
	{
		/* every vertex onto its face's surface: the footpoint becomes (face, u, v) */
		qaws_hm* m = &out->hm;
		for (i = 0; i < m->nv; i++)
		{
			qaws_mesh_location l = m->foot[i];
			unsigned int t = l.patch, F, tv[3];
			double w[3], u = 0, v = 0, y[3], d;
			qaws_surface const* s;
			if (t >= ref->nt)
				continue;
			F = (unsigned int)hm_attr_at(&ref->mesh->hm, fa, ref->tf[t])[0];
			s = mg_face(src, F);
			if (!s)
				continue;
			for (k = 0; k < 3; k++)
				tv[k] = ref->tv[3 * (size_t)t + k] < nv ? ref->tv[3 * (size_t)t + k] : 0;
			w[1] = (double)l.a;
			w[2] = (double)l.b;
			w[0] = 1 - w[1] - w[2];
			mg_start(s, dvp, tv, w, 3, F, m->x + 3 * (size_t)i, &u, &v);
			mg_snap(s, m->x + 3 * (size_t)i, &u, &v, (m->vflag[i] & HM_CORNER) ? 2 : hm_feature_degree(m, i) || hm_is_border_vertex(m, i) ? 1 : 0, y);
			d = sqrt((y[0] - m->x[3 * (size_t)i]) * (y[0] - m->x[3 * (size_t)i]) + (y[1] - m->x[3 * (size_t)i + 1]) * (y[1] - m->x[3 * (size_t)i + 1])
				+ (y[2] - m->x[3 * (size_t)i + 2]) * (y[2] - m->x[3 * (size_t)i + 2]));
			dmax = d > dmax ? d : dmax;
			dsum += d;
			memcpy(m->x + 3 * (size_t)i, y, sizeof y);
			m->foot[i].patch = F;
			m->foot[i].a = (qaws_scalar)u;
			m->foot[i].b = (qaws_scalar)v;
		}
		rep->max_deviation = (qaws_scalar)dmax;
		rep->mean_deviation = (qaws_scalar)(m->nv ? dsum / m->nv : 0);
		qaws_mesh_3d_measure_quality(&rep->quality, out);
	}
done:
	qaws_size_field_destroy(sf);
	qaws_mesh_reference_destroy(ref);
	qaws_mesh_3d_destroy(dense);
	if (st != QAWS_STATUS_OK)
	{
		qaws_mesh_3d_destroy(out);
		return st;
	}
	*out_mesh = out;
	return QAWS_STATUS_OK;
}

static int mg_check(qaws_generate_desc const* desc)
{
	return desc && (unsigned int)desc->method <= QAWS_GENERATE_STRUCTURED && (unsigned int)desc->elements <= QAWS_ELEMENTS_QUAD_DOMINANT
		&& desc->target_length >= 0 && desc->chord_tolerance >= 0 && desc->angle_tolerance >= 0 && desc->min_length >= 0 && desc->max_length >= 0;
}

static qaws_status mg_structured(mg_src const* src, unsigned int nf, qaws_generate_desc const* desc, double L, qaws_mesh_3d** out_mesh,
	qaws_mesh_report* rep);
static qaws_status mg_parametric(mg_src const* src, unsigned int nf, qaws_generate_desc const* desc, double L, qaws_mesh_3d** out_mesh,
	qaws_mesh_report* rep);
static qaws_status mg_delaunay(mg_src const* src, unsigned int nf, qaws_generate_desc const* desc, double L, qaws_mesh_3d** out_mesh,
	qaws_mesh_report* rep);

/* the dense mesh's size: half the target, finer for an adaptive one with a small minimum */
static double mg_dense_length(qaws_generate_desc const* desc, double L)
{
	double Ld = 0.5 * L;
	if (desc->target_length == 0 && (desc->chord_tolerance > 0 || desc->angle_tolerance > 0) && desc->min_length > 0)
		Ld = fmax(0.5 * (double)desc->min_length, L / 8);
	return Ld;
}

qaws_status qaws_mesh_generate_shell(qaws_mesh_3d** out_mesh, qaws_mesh_report* out_report, qaws_shell_atlas const* shell,
	qaws_generate_desc const* desc)
{
	qaws_shell_region* whole = NULL;
	qaws_shell_mesh* sm = NULL;
	qaws_surface_mesh_desc smd;
	qaws_mesh_report rep;
	qaws_shell_point* dvp = NULL;
	double *x = NULL, L;
	unsigned int *tri = NULL, *tface = NULL, nv, nt, i;
	qaws_scalar area = 0;
	mg_src src;
	qaws_status st;
	if (!out_mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_mesh = NULL;
	if (!shell || !mg_check(desc))
		return QAWS_STATUS_INVALID_ARGUMENT;
	if ((desc->method == QAWS_GENERATE_PARAMETRIC || desc->method == QAWS_GENERATE_RESTRICTED_DELAUNAY) && desc->elements != QAWS_ELEMENTS_TRIANGLES)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	memset(&rep, 0, sizeof rep);
	st = qaws_shell_region_create_whole(&whole, shell);
	if (st != QAWS_STATUS_OK)
		return st;
	qaws_shell_region_get_area(&area, whole);
	L = desc->target_length > 0 ? (double)desc->target_length : desc->max_length > 0 ? (double)desc->max_length : sqrt((double)area) / 30;
	if (!(L > 0))
	{
		st = QAWS_STATUS_INVALID_ARGUMENT;
		goto done;
	}
	if (desc->method != QAWS_GENERATE_REMESH)
	{
		unsigned int nfc = qaws_shell_atlas_get_face_count(shell);
		src.shell = shell;
		src.one = NULL;
		st = desc->method == QAWS_GENERATE_PARAMETRIC ? mg_parametric(&src, nfc, desc, L, out_mesh, &rep)
			: desc->method == QAWS_GENERATE_RESTRICTED_DELAUNAY ? mg_delaunay(&src, nfc, desc, L, out_mesh, &rep) : mg_structured(&src, nfc, desc, L, out_mesh, &rep);
		goto done;
	}
	/* the dense conformal mesh of the whole shell */
	memset(&smd, 0, sizeof smd);
	smd.edge_length = (qaws_scalar)mg_dense_length(desc, L);
	st = qaws_shell_region_triangulate(&sm, whole, &smd);
	if (st != QAWS_STATUS_OK)
		goto done;
	nv = qaws_shell_mesh_get_vertex_count(sm);
	nt = qaws_shell_mesh_get_triangle_count(sm);
	x = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 3, (nv + 1)));
	dvp = (qaws_shell_point*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_shell_point), (nv + 1)));
	tri = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (nt + 1)));
	tface = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (nt + 1)));
	if (!x || !dvp || !tri || !tface)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	for (i = 0; i < nv; i++)
	{
		qaws_vec3 p;
		qaws_shell_mesh_get_vertex(&dvp[i], &p, sm, i);
		mg_v3(p, x + 3 * (size_t)i);
	}
	for (i = 0; i < nt; i++)
	{
		qaws_shell_mesh_get_triangle(tri + 3 * (size_t)i, sm, i);
		tface[i] = qaws_shell_mesh_get_triangle_face(sm, i);
	}
	src.shell = shell;
	src.one = NULL;
	st = mg_finish(&src, x, dvp, nv, tri, tface, nt, desc, L, out_mesh, &rep);
done:
	qaws_internal_free(x); qaws_internal_free(dvp); qaws_internal_free(tri); qaws_internal_free(tface);
	qaws_shell_mesh_destroy(sm);
	qaws_shell_region_destroy(whole);
	if (st == QAWS_STATUS_OK && out_report)
		*out_report = rep;
	return st;
}

/* the cell key of a coordinate */
static long long mg_cell(double x, double cell)
{
	return (long long)floor(x / cell);
}

/* a surface's dense grid: n per direction by its iso-lines' lengths, the
   samples welded where they meet in space (seams, poles), degenerate
   triangles dropped; the shorter diagonal of each cell */
static qaws_status mg_grid(qaws_surface const* s, double Ld, double** ox, qaws_shell_point** odvp, unsigned int* onv, unsigned int** otri,
	unsigned int* ont, double* out_area)
{
	enum { G = 16 };
	qaws_range ur = qaws_surface_get_u_range(s), vr = qaws_surface_get_v_range(s);
	double u0 = (double)ur.min_value, u1 = (double)ur.max_value, v0 = (double)vr.min_value, v1 = (double)vr.max_value;
	double lu = 0, lv = 0, ext, lo[3] = { 1e300, 1e300, 1e300 }, hi[3] = { -1e300, -1e300, -1e300 }, *x = NULL, area = 0, cell;
	double coarse[3 * (G + 1) * (G + 1)];
	unsigned int nu, nw, i, j, k, n, nt = 0, *map = NULL, *tri = NULL, *slot = NULL, mask, nv = 0;
	qaws_shell_point* dvp = NULL;
	qaws_status st = QAWS_STATUS_OK;
	/* the iso-lines' lengths on a coarse grid */
	for (i = 0; i <= G; i++)
		for (j = 0; j <= G; j++)
		{
			qaws_surface_eval_result e;
			memset(&e, 0, sizeof e);
			qaws_surface_evaluate(&e, s, (qaws_scalar)(u0 + (u1 - u0) * i / G), (qaws_scalar)(v0 + (v1 - v0) * j / G), QAWS_SURFACE_EVAL_POSITION);
			mg_v3(e.position, coarse + 3 * (i * (G + 1) + j));
			for (k = 0; k < 3; k++)
			{
				lo[k] = fmin(lo[k], coarse[3 * (i * (G + 1) + j) + k]);
				hi[k] = fmax(hi[k], coarse[3 * (i * (G + 1) + j) + k]);
			}
		}
	for (j = 0; j <= G; j++)
	{
		double a = 0, b = 0;
		for (i = 0; i < G; i++)
		{
			double const *p = coarse + 3 * (i * (G + 1) + j), *q = coarse + 3 * ((i + 1) * (G + 1) + j);
			double const *r = coarse + 3 * (j * (G + 1) + i), *t = coarse + 3 * (j * (G + 1) + i + 1);
			a += sqrt((p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]) + (p[2] - q[2]) * (p[2] - q[2]));
			b += sqrt((r[0] - t[0]) * (r[0] - t[0]) + (r[1] - t[1]) * (r[1] - t[1]) + (r[2] - t[2]) * (r[2] - t[2]));
		}
		lu = fmax(lu, a);
		lv = fmax(lv, b);
	}
	ext = sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) + (hi[2] - lo[2]) * (hi[2] - lo[2]));
	cell = 1e-6 * (ext > 0 ? ext : 1);
	nu = (unsigned int)fmin(1024, fmax(4, ceil(lu / Ld)));
	nw = (unsigned int)fmin(1024, fmax(4, ceil(lv / Ld)));
	n = (nu + 1) * (nw + 1);
	x = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, n));
	dvp = (qaws_shell_point*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_shell_point), n));
	map = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), n));
	tri = (unsigned int*)qaws_internal_malloc(qaws_internal_mul4(sizeof(unsigned int), 6, nu, nw));
	for (mask = 1; mask < 2 * n; mask <<= 1) {}
	slot = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), mask));
	if (!x || !dvp || !map || !tri || !slot)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	memset(slot, 0xFF, sizeof(unsigned int) * mask);
	mask -= 1;
	/* the samples, welded by position (the cell's neighbours searched too) */
	for (i = 0; i <= nu; i++)
		for (j = 0; j <= nw; j++)
		{
			qaws_surface_eval_result e;
			double p[3], u = u0 + (u1 - u0) * i / nu, v = v0 + (v1 - v0) * j / nw;
			unsigned int q, found = HM_NONE;
			int dx, dy, dz;
			memset(&e, 0, sizeof e);
			qaws_surface_evaluate(&e, s, (qaws_scalar)u, (qaws_scalar)v, QAWS_SURFACE_EVAL_POSITION);
			mg_v3(e.position, p);
			for (dx = -1; dx <= 1 && found == HM_NONE; dx++)
				for (dy = -1; dy <= 1 && found == HM_NONE; dy++)
					for (dz = -1; dz <= 1 && found == HM_NONE; dz++)
					{
						unsigned int hh = (unsigned int)((mg_cell(p[0], cell) + dx) * 73856093LL ^ (mg_cell(p[1], cell) + dy) * 19349663LL
							^ (mg_cell(p[2], cell) + dz) * 83492791LL);
						for (q = hh & mask; slot[q] != HM_NONE; q = (q + 1) & mask)
						{
							double const* o = x + 3 * (size_t)slot[q];
							if (fabs(o[0] - p[0]) <= cell && fabs(o[1] - p[1]) <= cell && fabs(o[2] - p[2]) <= cell)
							{
								found = slot[q];
								break;
							}
						}
					}
			if (found == HM_NONE)
			{
				unsigned int h = (unsigned int)(mg_cell(p[0], cell) * 73856093LL ^ mg_cell(p[1], cell) * 19349663LL ^ mg_cell(p[2], cell) * 83492791LL);
				memcpy(x + 3 * (size_t)nv, p, sizeof p);
				dvp[nv].face = 0;
				dvp[nv].uv.x = (qaws_scalar)u;
				dvp[nv].uv.y = (qaws_scalar)v;
				for (q = h & mask; slot[q] != HM_NONE; q = (q + 1) & mask) {}
				slot[q] = nv;
				found = nv++;
			}
			map[i * (nw + 1) + j] = found;
		}
	/* two triangles a cell, counter-clockwise in (u, v) */
	for (i = 0; i < nu; i++)
		for (j = 0; j < nw; j++)
		{
			unsigned int a = map[i * (nw + 1) + j], b = map[(i + 1) * (nw + 1) + j], c = map[(i + 1) * (nw + 1) + j + 1], d = map[i * (nw + 1) + j + 1];
			unsigned int t6[6], q;
			double const *A = x + 3 * (size_t)a, *B = x + 3 * (size_t)b, *C = x + 3 * (size_t)c, *D = x + 3 * (size_t)d;
			double ac = (A[0] - C[0]) * (A[0] - C[0]) + (A[1] - C[1]) * (A[1] - C[1]) + (A[2] - C[2]) * (A[2] - C[2]);
			double bd = (B[0] - D[0]) * (B[0] - D[0]) + (B[1] - D[1]) * (B[1] - D[1]) + (B[2] - D[2]) * (B[2] - D[2]);
			if (ac <= bd)
			{
				t6[0] = a; t6[1] = b; t6[2] = c; t6[3] = a; t6[4] = c; t6[5] = d;
			}
			else
			{
				t6[0] = a; t6[1] = b; t6[2] = d; t6[3] = b; t6[4] = c; t6[5] = d;
			}
			for (q = 0; q < 2; q++)
			{
				unsigned int const* t3 = t6 + 3 * q;
				double const *P = x + 3 * (size_t)t3[0], *Q = x + 3 * (size_t)t3[1], *R = x + 3 * (size_t)t3[2];
				double e1[3], e2[3], cr[3], ar;
				if (t3[0] == t3[1] || t3[1] == t3[2] || t3[0] == t3[2])
					continue;   /* (at a pole) */
				for (k = 0; k < 3; k++)
				{
					e1[k] = Q[k] - P[k];
					e2[k] = R[k] - P[k];
				}
				cr[0] = e1[1] * e2[2] - e1[2] * e2[1]; cr[1] = e1[2] * e2[0] - e1[0] * e2[2]; cr[2] = e1[0] * e2[1] - e1[1] * e2[0];
				ar = 0.5 * sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]);
				if (!(ar > 1e-14 * ext * ext))
					continue;
				area += ar;
				memcpy(tri + 3 * (size_t)nt, t3, sizeof(unsigned int) * 3);
				nt++;
			}
		}
done:
	qaws_internal_free(map);
	qaws_internal_free(slot);
	if (st != QAWS_STATUS_OK)
	{
		qaws_internal_free(x);
		qaws_internal_free(dvp);
		qaws_internal_free(tri);
		return st;
	}
	*ox = x;
	*odvp = dvp;
	*onv = nv;
	*otri = tri;
	*ont = nt;
	*out_area = area;
	return QAWS_STATUS_OK;
}

/* ---- STRUCTURED: iso-parametric grids ----------------------------------- */

/* side k of the domain dom (u0, u1, v0, v1) at t in [0, 1], by increasing parameter: 0 v = v0, 1 u = u1, 2 v = v1, 3 u = u0 */
static void mg_side_uv(double const* dom, int k, double t, double* uv)
{
	uv[0] = k == 1 ? dom[1] : k == 3 ? dom[0] : dom[0] + t * (dom[1] - dom[0]);
	uv[1] = k == 0 ? dom[2] : k == 2 ? dom[3] : dom[2] + t * (dom[3] - dom[2]);
}

static void mg_at(qaws_surface const* s, double u, double v, double* p)
{
	qaws_surface_eval_result e;
	memset(&e, 0, sizeof e);
	qaws_surface_evaluate(&e, s, (qaws_scalar)u, (qaws_scalar)v, QAWS_SURFACE_EVAL_POSITION);
	mg_v3(e.position, p);
}

/* side k's length, and n + 1 fractions along it at equal lengths (equal in parameter on a side of no length: a pole) */
static double mg_side_params(qaws_surface const* s, double const* dom, int k, unsigned int n, double* out)
{
	enum { K = 128 };
	double cum[K + 1], prev[3], uv[2];
	unsigned int i, j = 0;
	cum[0] = 0;
	for (i = 0; i <= K; i++)
	{
		double p[3];
		mg_side_uv(dom, k, (double)i / K, uv);
		mg_at(s, uv[0], uv[1], p);
		if (i)
			cum[i] = cum[i - 1] + sqrt((p[0] - prev[0]) * (p[0] - prev[0]) + (p[1] - prev[1]) * (p[1] - prev[1]) + (p[2] - prev[2]) * (p[2] - prev[2]));
		memcpy(prev, p, sizeof prev);
	}
	for (i = 0; out && i <= n; i++)
	{
		double target = cum[K] * i / n;
		if (!(cum[K] > 0))
		{
			out[i] = (double)i / n;
			continue;
		}
		while (j < K && cum[j + 1] < target)
			j++;
		out[i] = j >= K ? 1 : (j + (cum[j + 1] > cum[j] ? (target - cum[j]) / (cum[j + 1] - cum[j]) : 0)) / K;
		out[i] = i == 0 ? 0 : i == n ? 1 : out[i];
	}
	return cum[K];
}

static unsigned int mg_root(unsigned int* uf, unsigned int i)
{
	while (uf[i] != i)
		i = uf[i] = uf[uf[i]];
	return i;
}

/* face f's grid: (nu + 1) (nw + 1) nodes' (u, v), node (i, j) at i + (nu + 1) j; each side split at equal lengths, the nodes
   inside blending opposite sides' fractions (u from sides 0 and 2 by j, v from sides 3 and 1 by i) */
static void mg_grid_uv(qaws_surface const* s, double const* dom, unsigned int nu, unsigned int nw, double* a, double* b, double* c, double* d,
	double* uv)
{
	unsigned int i, j;
	mg_side_params(s, dom, 0, nu, a);
	mg_side_params(s, dom, 2, nu, b);
	mg_side_params(s, dom, 3, nw, c);
	mg_side_params(s, dom, 1, nw, d);
	for (j = 0; j <= nw; j++)
		for (i = 0; i <= nu; i++)
		{
			double t = (double)j / nw, r = (double)i / nu, fu = (1 - t) * a[i] + t * b[i], fv = (1 - r) * c[j] + r * d[j];
			uv[2 * (i + (size_t)(nu + 1) * j)] = dom[0] + fu * (dom[1] - dom[0]);
			uv[2 * (i + (size_t)(nu + 1) * j) + 1] = dom[2] + fv * (dom[3] - dom[2]);
		}
}

/* PARAMETRIC: a face's start mesh (its boundary polygon triangulated, flat at z = 0 in (u, v)) remeshed by LOCAL in the
   metric of the surface's first fundamental form over the size squared (sampled on a 64 x 64 grid of the domain D, the
   reference), the start's first npin vertices (the boundary) pinned, so
   faces still meet point to point: the size L, or with tolerances the chord's sqrt(8 chord / k) and the angle's angle / k
   at the larger principal curvature k, within [min_length, L]. Out: the vertices' (u, v) and the triangles,
   counter-clockwise in (u, v). */
static qaws_status mg_param_face(qaws_surface const* s, double const* D, double const* uv, unsigned int nv, unsigned int npin, unsigned int const* tri0, unsigned int nt0,
	double L, qaws_generate_desc const* desc, double** ouv, unsigned int* onv, unsigned int** otri, unsigned int* ont)
{
	enum { GN = 64, GV = (GN + 1) * (GN + 1) };
	unsigned int i, j, *tri = NULL, *pin = NULL;
	double* gx = NULL;
	qaws_mesh_3d* gm = NULL;
	double* x = NULL;
	qaws_scalar *met = NULL, *pos = NULL;
	qaws_mesh_3d *pm = NULL, *om = NULL;
	qaws_mesh_reference* ref = NULL;
	qaws_metric_field* mf = NULL;
	qaws_mesh_reference_desc rd;
	qaws_remesh_desc ld;
	int adaptive = desc->target_length == 0 && (desc->chord_tolerance > 0 || desc->angle_tolerance > 0);
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	x = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)nv + 1));
	pin = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)nv + 1));
	met = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 6, GV));
	gx = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 3, GV));
	tri = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 6, GN * GN));
	if (!x || !pin || !met || !gx || !tri)
		goto done;
	for (j = 0; j <= GN; j++)
		for (i = 0; i <= GN; i++)
		{
			unsigned int v = i + (GN + 1) * j;
			gx[3 * (size_t)v] = D[0] + (D[1] - D[0]) * i / GN;
			gx[3 * (size_t)v + 1] = D[2] + (D[3] - D[2]) * j / GN;
			gx[3 * (size_t)v + 2] = 0;
			if (i < GN && j < GN)
			{
				unsigned int* q = tri + 6 * (size_t)(i + GN * j);
				q[0] = v; q[1] = v + 1; q[2] = v + GN + 2;
				q[3] = v; q[4] = v + GN + 2; q[5] = v + GN + 1;
			}
		}
	if ((st = hm_mesh_3d_create(gx, GV, tri, NULL, 2 * GN * GN, &gm)) != QAWS_STATUS_OK)
		goto done;
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < nv; i++)
	{
		x[3 * (size_t)i] = uv[2 * (size_t)i];
		x[3 * (size_t)i + 1] = uv[2 * (size_t)i + 1];
		x[3 * (size_t)i + 2] = 0;
		pin[i] = i;
	}
	if ((st = hm_mesh_3d_create(x, nv, tri0, NULL, nt0, &pm)) != QAWS_STATUS_OK)
		goto done;
	memset(&rd, 0, sizeof rd);
	rd.interpolant = QAWS_INTERPOLANT_LINEAR;
	if ((st = qaws_mesh_reference_create(&ref, gm, &rd)) != QAWS_STATUS_OK)
		goto done;
	/* each grid node's metric: the first fundamental form over its size squared (z: their mean, the plane's normal) */
	for (i = 0; i < GV; i++)
	{
		qaws_surface_eval_result e;
		double du[3], dv[3], E, F, G, h = L, eps;
		memset(&e, 0, sizeof e);
		qaws_surface_evaluate(&e, s, (qaws_scalar)gx[3 * (size_t)i], (qaws_scalar)gx[3 * (size_t)i + 1],
			QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV | QAWS_SURFACE_EVAL_NORMAL | QAWS_SURFACE_EVAL_DUU | QAWS_SURFACE_EVAL_DUV | QAWS_SURFACE_EVAL_DVV);
		mg_v3(e.du, du);
		mg_v3(e.dv, dv);
		E = mg_dot(du, du);
		F = mg_dot(du, dv);
		G = mg_dot(dv, dv);
		if (adaptive)
		{
			double n[3], uu[3], uw[3], ww[3], l2, m2, n2, det = E * G - F * F, H, K, km;
			mg_v3(e.normal, n);
			mg_v3(e.duu, uu);
			mg_v3(e.duv, uw);
			mg_v3(e.dvv, ww);
			l2 = mg_dot(uu, n);
			m2 = mg_dot(uw, n);
			n2 = mg_dot(ww, n);
			H = det > 0 ? (l2 * G + n2 * E - 2 * m2 * F) / (2 * det) : 0;
			K = det > 0 ? (l2 * n2 - m2 * m2) / det : 0;
			km = fabs(H) + sqrt(fmax(0, H * H - K));
			if (km > 0 && desc->chord_tolerance > 0)
				h = fmin(h, sqrt(8 * (double)desc->chord_tolerance / km));
			if (km > 0 && desc->angle_tolerance > 0)
				h = fmin(h, (double)desc->angle_tolerance / km);
			h = fmax(h, desc->min_length > 0 ? (double)desc->min_length : 1e-3 * L);
		}
		/* (a pole's degenerate form kept invertible) */
		eps = 1e-9 * (E + G) + 1e-300;
		met[6 * (size_t)i] = (qaws_scalar)((E + eps) / (h * h));
		met[6 * (size_t)i + 1] = (qaws_scalar)(F / (h * h));
		met[6 * (size_t)i + 2] = 0;
		met[6 * (size_t)i + 3] = (qaws_scalar)((G + eps) / (h * h));
		met[6 * (size_t)i + 4] = 0;
		met[6 * (size_t)i + 5] = (qaws_scalar)(0.5 * (E + G) / (h * h));
	}
	if ((st = qaws_metric_field_create_vertices(&mf, ref, met)) != QAWS_STATUS_OK)
		goto done;
	memset(&ld, 0, sizeof ld);
	ld.method = QAWS_REMESH_LOCAL;
	ld.reference = ref;
	ld.metric = mf;
	ld.pinned = pin;
	ld.pinned_count = npin;
	if ((st = qaws_mesh_remesh(&om, NULL, pm, &ld)) != QAWS_STATUS_OK)
		goto done;
	*onv = qaws_mesh_3d_get_vertex_count(om);
	*ont = qaws_mesh_3d_get_triangle_count(om);
	pos = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)*onv + 1));
	*ouv = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 2, (size_t)*onv + 1));
	*otri = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)*ont + 1));
	if (!pos || !*ouv || !*otri)
	{
		qaws_internal_free(*ouv);
		qaws_internal_free(*otri);
		*ouv = NULL;
		*otri = NULL;
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	qaws_mesh_3d_get_positions(pos, om);
	for (i = 0; i < *onv; i++)
	{
		(*ouv)[2 * (size_t)i] = (double)pos[3 * (size_t)i];
		(*ouv)[2 * (size_t)i + 1] = (double)pos[3 * (size_t)i + 1];
	}
	qaws_mesh_3d_get_triangles(*otri, om);
	st = QAWS_STATUS_OK;
done:
	qaws_internal_free(x);
	qaws_internal_free(tri);
	qaws_internal_free(pin);
	qaws_internal_free(met);
	qaws_internal_free(pos);
	qaws_metric_field_destroy(mf);
	qaws_mesh_reference_destroy(ref);
	qaws_mesh_3d_destroy(pm);
	qaws_mesh_3d_destroy(om);
	qaws_mesh_3d_destroy(gm);
	qaws_internal_free(gx);
	return st;
}

/* the vertex at p (face f, (u, v)): the one already within cell of it (welded: hashed in slot) or a new one, hashed when
   weld (a point where faces may meet) */
static unsigned int mg_weld_point(double* x, qaws_mesh_location* foot, unsigned int* nv, unsigned int* slot, unsigned int mask, double cell,
	double const* p, unsigned int f, double u, double v, int weld)
{
	unsigned int q, found = HM_NONE, hh;
	int dx, dy, dz;
	for (dx = -1; weld && dx <= 1 && found == HM_NONE; dx++)
		for (dy = -1; dy <= 1 && found == HM_NONE; dy++)
			for (dz = -1; dz <= 1 && found == HM_NONE; dz++)
			{
				hh = (unsigned int)((mg_cell(p[0], cell) + dx) * 73856093LL ^ (mg_cell(p[1], cell) + dy) * 19349663LL ^ (mg_cell(p[2], cell) + dz) * 83492791LL);
				for (q = hh & mask; slot[q] != HM_NONE; q = (q + 1) & mask)
				{
					double const* o = x + 3 * (size_t)slot[q];
					if (fabs(o[0] - p[0]) <= cell && fabs(o[1] - p[1]) <= cell && fabs(o[2] - p[2]) <= cell)
					{
						found = slot[q];
						break;
					}
				}
			}
	if (found != HM_NONE)
		return found;
	memcpy(x + 3 * (size_t)*nv, p, sizeof(double) * 3);
	foot[*nv].patch = f;
	foot[*nv].a = (qaws_scalar)u;
	foot[*nv].b = (qaws_scalar)v;
	if (weld)
	{
		hh = (unsigned int)(mg_cell(p[0], cell) * 73856093LL ^ mg_cell(p[1], cell) * 19349663LL ^ mg_cell(p[2], cell) * 83492791LL);
		for (q = hh & mask; slot[q] != HM_NONE; q = (q + 1) & mask) {}
		slot[q] = *nv;
	}
	return (*nv)++;
}

/* the generated mesh from its arrays: footpoints (face, u, v), the "face" attribute, the edges between faces features and
   their meeting points corners, the quality reported */
static qaws_status mg_assemble(double const* x, qaws_mesh_location const* foot, unsigned int nv, unsigned int const* fv, unsigned int const* fs,
	unsigned int const* fface, unsigned int nface, qaws_mesh_3d** out_mesh, qaws_mesh_report* rep)
{
	qaws_mesh_3d* out = NULL;
	qaws_status st;
	st = hm_mesh_3d_create(x, nv, fv, fs, nface, &out);
	if (st != QAWS_STATUS_OK)
		goto done;
	if (!hm_foot_enable(&out->hm))
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	{
		/* footpoints (face, u, v); the "face" attribute; edges between faces features, their meeting points corners */
		qaws_hm* m = &out->hm;
		qaws_attribute_desc ad;
		qaws_scalar* fa = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(qaws_scalar), (size_t)m->nf + 1));
		unsigned int id, c, v;
		if (!fa)
		{
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			goto done;
		}
		for (v = 0; v < m->nv; v++)
			m->foot[v] = foot[m->vsrc[v] < nv ? m->vsrc[v] : 0];
		for (c = 0; c < m->nf; c++)
			fa[c] = (qaws_scalar)(c < nface ? fface[c] : 0);
		memset(&ad, 0, sizeof ad);
		ad.name = "face";
		ad.element = QAWS_MESH_FACE;
		ad.components = 1;
		ad.rule = QAWS_ATTRIBUTE_NEAREST;
		st = qaws_mesh_3d_add_attribute(&id, out, fa, &ad);
		qaws_internal_free(fa);
		if (st != QAWS_STATUS_OK)
			goto done;
		for (c = 0; c < m->nc; c++)
		{
			unsigned int o = m->co[c];
			if (m->cf[c] == HM_NONE || o == HM_NONE)
				continue;
			if (hm_attr_at(m, id, m->cf[c])[0] != hm_attr_at(m, id, m->cf[o])[0])
			{
				m->cflag[c] |= HM_FEATURE;
				m->cflag[o] |= HM_FEATURE;
			}
		}
		for (v = 0; v < m->nv; v++)
		{
			unsigned int fd = hm_feature_degree(m, v);
			if (fd != 0 && fd != 2)
				m->vflag[v] |= HM_CORNER;
		}
	}
	memset(rep, 0, sizeof *rep);
	qaws_mesh_3d_measure_quality(&rep->quality, out);
	*out_mesh = out;
	return QAWS_STATUS_OK;
done:
	qaws_mesh_3d_destroy(out);
	return st;
}

/* PARAMETRIC: per face its boundary split in space, each side by equal steps of its length over the size (L, or with
   tolerances the side curve's own chord and angle sizes from its curvature, within [min_length, L]: the same curve from
   both faces of an edge, so the same points), the counts made the larger across shared edges; the polygon in (u, v)
   fanned from the domain's centre and remeshed in the surface's metric with the boundary pinned (mg_param_face); every
   vertex lifted, the boundary's welded where faces meet (shared edges, seams, poles). Untrimmed faces. */
static qaws_status mg_parametric(mg_src const* src, unsigned int nf, qaws_generate_desc const* desc, double L, qaws_mesh_3d** out_mesh,
	qaws_mesh_report* rep)
{
	enum { K = 256 };
	unsigned int *n = NULL, f, k, i, nv = 0, nc = 0, nface = 0, cap = 0, ccap = 0, fcap = 0, *fv = NULL, *fs = NULL, *fface = NULL, *slot = NULL, mask, nb[4];
	int nbs[4], adaptive = desc->target_length == 0 && (desc->chord_tolerance > 0 || desc->angle_tolerance > 0);
	double *dom = NULL, *W = NULL, *x = NULL, minstep = HUGE_VAL, cell, hmin = desc->min_length > 0 ? (double)desc->min_length : 1e-3 * L;
	qaws_mesh_location* foot = NULL;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	n = (unsigned int*)qaws_internal_scratch_zalloc(4 * (size_t)nf + 4, sizeof(unsigned int));
	dom = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 4, (size_t)nf + 1));
	W = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 4 * (size_t)nf + 4, K + 1));
	if (!n || !dom || !W)
		goto done;
	/* each side's weight sum_ds / h (the sizes from the side curve's turning between samples) and count */
	for (f = 0; f < nf; f++)
	{
		qaws_surface const* s = mg_face(src, f);
		qaws_range ur = qaws_surface_get_u_range(s), vr = qaws_surface_get_v_range(s);
		if (src->shell && !sh_face_grid_sides(src->shell, f, nb, nbs))
		{
			st = QAWS_STATUS_UNSUPPORTED_OPERATION;   /* (a trimmed face) */
			goto done;
		}
		dom[4 * (size_t)f] = (double)ur.min_value;
		dom[4 * (size_t)f + 1] = (double)ur.max_value;
		dom[4 * (size_t)f + 2] = (double)vr.min_value;
		dom[4 * (size_t)f + 3] = (double)vr.max_value;
		for (k = 0; k < 4; k++)
		{
			double P[K + 1][3], hs[K + 1], len = 0, *w = W + (4 * (size_t)f + k) * (K + 1), uv[2];
			for (i = 0; i <= K; i++)
			{
				mg_side_uv(dom + 4 * (size_t)f, (int)k, (double)i / K, uv);
				mg_at(s, uv[0], uv[1], P[i]);
			}
			for (i = 0; i <= K; i++)
			{
				unsigned int a = i ? i - 1 : 0, b = i < K ? i + 1 : K, q;
				double e1[3], e2[3], l1, l2, c, kap = 0, h = L;
				for (q = 0; q < 3; q++)
				{
					e1[q] = P[i][q] - P[a][q];
					e2[q] = P[b][q] - P[i][q];
				}
				l1 = sqrt(mg_dot(e1, e1));
				l2 = sqrt(mg_dot(e2, e2));
				if (l1 > 0 && l2 > 0)
				{
					c = fmax(-1, fmin(1, mg_dot(e1, e2) / (l1 * l2)));
					kap = acos(c) / (0.5 * (l1 + l2));
				}
				if (adaptive && kap > 0)
				{
					if (desc->chord_tolerance > 0)
						h = fmin(h, sqrt(8 * (double)desc->chord_tolerance / kap));
					if (desc->angle_tolerance > 0)
						h = fmin(h, (double)desc->angle_tolerance / kap);
				}
				hs[i] = fmax(h, hmin);
			}
			w[0] = 0;
			for (i = 1; i <= K; i++)
			{
				double d[3], ds;
				d[0] = P[i][0] - P[i - 1][0];
				d[1] = P[i][1] - P[i - 1][1];
				d[2] = P[i][2] - P[i - 1][2];
				ds = sqrt(mg_dot(d, d));
				len += ds;
				w[i] = w[i - 1] + 0.5 * ds * (1 / hs[i - 1] + 1 / hs[i]);
			}
			n[4 * f + k] = (unsigned int)fmin(65536, fmax(1, ceil(w[K] - 1e-9)));
			if (len > 1e-12 * L)
				minstep = fmin(minstep, len / n[4 * f + k]);
		}
	}
	/* the counts the larger across shared edges (twice: chains of edges settle) */
	for (k = 0; src->shell && k < 2; k++)
		for (f = 0; f < nf; f++)
		{
			unsigned int q;
			sh_face_grid_sides(src->shell, f, nb, nbs);
			for (q = 0; q < 4; q++)
				if (nb[q] < nf && nbs[q] >= 0 && nbs[q] < 4)
				{
					unsigned int a = 4 * f + q, b = 4 * nb[q] + (unsigned int)nbs[q], m = n[a] > n[b] ? n[a] : n[b];
					n[a] = n[b] = m;
				}
		}
	/* the weld: boundary points only, a ten-thousandth of the finest boundary step */
	{
		size_t nbp = 0;
		for (i = 0; i < 4 * nf; i++)
			nbp += n[i] + 1;
		for (mask = 1; mask < 2 * nbp + 2; mask <<= 1) {}
		slot = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), mask));
		if (!slot)
			goto done;
		memset(slot, 0xFF, sizeof(unsigned int) * mask);
		mask -= 1;
		cell = 1e-4 * (minstep < HUGE_VAL ? minstep : L);
	}
	for (f = 0; f < nf; f++)
	{
		qaws_surface const* s = mg_face(src, f);
		double const* D = dom + 4 * (size_t)f;
		double *poly = NULL, *uv0 = NULL, *puv = NULL, eu = 1e-9 * (D[1] - D[0]), ev = 1e-9 * (D[3] - D[2]);
		unsigned int np = 0, *t0 = NULL, n0 = 0, nt0 = 0, pnv = 0, pnt = 0, *pt = NULL, *pmap = NULL, t;
		for (k = 0; k < 4; k++)
			np += n[4 * f + k];
		poly = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 2, (size_t)np + 2));
		if (!poly)
			goto done;
		/* the polygon counter-clockwise: sides 0 and 1 forward, 2 and 3 backward, each without its last point */
		np = 0;
		for (k = 0; k < 4; k++)
		{
			unsigned int m = n[4 * f + k], j, z;
			double const* w = W + (4 * (size_t)f + k) * (K + 1);
			for (j = 0; j < m; j++)
			{
				unsigned int jj = k < 2 ? j : m - j;
				double target = w[K] * jj / m, tt, uv[2];
				for (z = 0; z < K && w[z + 1] < target; z++) {}
				tt = jj == 0 ? 0 : jj == m ? 1 : (z + (w[z + 1] > w[z] ? (target - w[z]) / (w[z + 1] - w[z]) : 0)) / K;
				mg_side_uv(D, (int)k, tt, uv);
				poly[2 * (size_t)np] = uv[0];
				poly[2 * (size_t)np + 1] = uv[1];
				np++;
			}
		}
		/* the start: a fan from the domain's centre (the domain convex, its polygon's points kept: they are pinned) */
		uv0 = poly;
		uv0[2 * (size_t)np] = 0.5 * (D[0] + D[1]);
		uv0[2 * (size_t)np + 1] = 0.5 * (D[2] + D[3]);
		n0 = np + 1;
		nt0 = np;
		t0 = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)np + 1));
		if (!t0)
		{
			qaws_internal_free(poly);
			goto done;
		}
		for (i = 0; i < np; i++)
		{
			t0[3 * (size_t)i] = np;
			t0[3 * (size_t)i + 1] = i;
			t0[3 * (size_t)i + 2] = (i + 1) % np;
		}
		st = mg_param_face(s, D, uv0, n0, np, t0, nt0, L, desc, &puv, &pnv, &pt, &pnt);
		qaws_internal_free(uv0);
		qaws_internal_free(t0);
		if (st != QAWS_STATUS_OK)
			goto done;
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		/* the arrays grown to take the face */
		if (nv + pnv + 1 > cap)
		{
			double* x2;
			qaws_mesh_location* f2;
			cap = 2 * (nv + pnv + 1);
			x2 = (double*)qaws_internal_realloc(x, qaws_internal_mul3(sizeof(double), 3, cap));
			x = x2 ? x2 : x;
			f2 = (qaws_mesh_location*)qaws_internal_realloc(foot, qaws_internal_mul2(sizeof(qaws_mesh_location), cap));
			foot = f2 ? f2 : foot;
			if (!x2 || !f2)
				cap = 0;
		}
		if (nc + 3 * pnt + 1 > ccap)
		{
			unsigned int* v2;
			ccap = 2 * (nc + 3 * pnt + 1);
			v2 = (unsigned int*)qaws_internal_realloc(fv, qaws_internal_mul2(sizeof(unsigned int), ccap));
			fv = v2 ? v2 : fv;
			if (!v2)
				ccap = 0;
		}
		if (nface + pnt + 1 > fcap)
		{
			unsigned int *s2, *g2;
			fcap = 2 * (nface + pnt + 1);
			s2 = (unsigned int*)qaws_internal_realloc(fs, qaws_internal_mul2(sizeof(unsigned int), fcap));
			fs = s2 ? s2 : fs;
			g2 = (unsigned int*)qaws_internal_realloc(fface, qaws_internal_mul2(sizeof(unsigned int), fcap));
			fface = g2 ? g2 : fface;
			if (!s2 || !g2)
				fcap = 0;
		}
		pmap = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)pnv + 1));
		if (!cap || !ccap || !fcap || !pmap)
		{
			qaws_internal_free(puv); qaws_internal_free(pt); qaws_internal_free(pmap);
			goto done;
		}
		for (i = 0; i < pnv; i++)
		{
			double p[3], u = puv[2 * (size_t)i], v = puv[2 * (size_t)i + 1];
			int side = fabs(u - D[0]) <= eu || fabs(u - D[1]) <= eu || fabs(v - D[2]) <= ev || fabs(v - D[3]) <= ev;
			mg_at(s, u, v, p);
			pmap[i] = mg_weld_point(x, foot, &nv, slot, mask, cell, p, f, u, v, side);
		}
		for (t = 0; t < pnt; t++)
		{
			unsigned int a = pmap[pt[3 * (size_t)t]], b = pmap[pt[3 * (size_t)t + 1]], c = pmap[pt[3 * (size_t)t + 2]];
			if (a == b || b == c || a == c)
				continue;   /* (at a pole) */
			fv[nc++] = a;
			fv[nc++] = b;
			fv[nc++] = c;
			fs[nface] = 3;
			fface[nface++] = f;
		}
		qaws_internal_free(puv);
		qaws_internal_free(pt);
		qaws_internal_free(pmap);
	}
	st = mg_assemble(x, foot, nv, fv, fs, fface, nface, out_mesh, rep);
done:
	qaws_internal_free(n); qaws_internal_free(dom); qaws_internal_free(W); qaws_internal_free(x); qaws_internal_free(foot);
	qaws_internal_free(fv); qaws_internal_free(fs); qaws_internal_free(fface); qaws_internal_free(slot);
	return st;
}

/* ---- RESTRICTED_DELAUNAY: Delaunay refinement on the faces ----------------- */

typedef struct mg_rd
{
	mg_src const* src;
	qaws_hm* m;
	unsigned int fa;            /* the "face" attribute */
	double L, chord, cosmax;    /* the size, the distance bound, cos of the angle bound */
	unsigned int* stack;        /* corners whose edges to check */
	unsigned int ns, cap;
} mg_rd;

static unsigned int mg_rd_face(mg_rd const* r, unsigned int f)
{
	return (unsigned int)hm_attr_at(r->m, r->fa, f)[0];
}

static int mg_rd_fixed(qaws_hm const* m, unsigned int c)
{
	return m->co[c] == HM_NONE || (m->cflag[c] & HM_FEATURE);
}

/* the angle at corner c's vertex in its triangle, as its cosine */
static double mg_rd_cos(qaws_hm const* m, unsigned int c)
{
	double const *a = m->x + 3 * (size_t)m->cv[c], *b = m->x + 3 * (size_t)hm_to(m, c), *p = m->x + 3 * (size_t)m->cv[hm_prev(m, c)];
	double u[3], w[3], l;
	int k;
	for (k = 0; k < 3; k++)
	{
		u[k] = b[k] - a[k];
		w[k] = p[k] - a[k];
	}
	l = sqrt(mg_dot(u, u) * mg_dot(w, w));
	return l > 0 ? mg_dot(u, w) / l : 1;
}

/* edge c not Delaunay: the angles facing it sum past pi (its opposite corners' cosines summing below 0) */
static int mg_rd_bad_edge(qaws_hm const* m, unsigned int c)
{
	unsigned int o = m->co[c];
	if (mg_rd_fixed(m, c) || m->fn[m->cf[c]] != 3 || m->fn[m->cf[o]] != 3)
		return 0;
	/* (sin(a + b) < 0 for the angles a, b at the corners facing the edge) */
	{
		double ca = mg_rd_cos(m, hm_prev(m, c)), cb = mg_rd_cos(m, hm_prev(m, o));
		double sa = sqrt(fmax(0, 1 - ca * ca)), sb = sqrt(fmax(0, 1 - cb * cb));
		return sa * cb + ca * sb < -1e-12;
	}
}

static int mg_rd_push(mg_rd* r, unsigned int c)
{
	if (r->ns == r->cap)
	{
		unsigned int nc = 2 * r->cap + 64, *s = (unsigned int*)qaws_internal_realloc(r->stack, qaws_internal_mul2(sizeof(unsigned int), nc));
		if (!s)
			return 0;
		r->stack = s;
		r->cap = nc;
	}
	r->stack[r->ns++] = c;
	return 1;
}

/* Lawson's flips until the stack's edges are Delaunay (each flip's two new outer edges checked in turn) */
static int mg_rd_flips(mg_rd* r)
{
	qaws_hm* m = r->m;
	unsigned int guard = 0;
	while (r->ns && guard++ < 1000000)
	{
		unsigned int c = r->stack[--r->ns], o, f1, f2, q;
		if (c >= m->nc || m->cf[c] == HM_NONE || !mg_rd_bad_edge(m, c) || !hm_can_flip(m, c))
			continue;
		o = m->co[c];
		f1 = m->cf[c];
		f2 = m->cf[o];
		if (!hm_flip(m, c))
			continue;
		/* (the four outer edges of the two triangles) */
		for (q = m->fc[f1]; q < m->fc[f1] + 3; q++)
			if (!mg_rd_push(r, q))
				return 0;
		for (q = m->fc[f2]; q < m->fc[f2] + 3; q++)
			if (!mg_rd_push(r, q))
				return 0;
	}
	return 1;
}

/* the (u, v) of vertex v on face F: its footpoint when there, else snapped from the start (u0, v0) held on F's side */
static void mg_rd_uv(mg_rd const* r, unsigned int v, unsigned int F, double u0, double v0, double* u, double* w)
{
	qaws_hm const* m = r->m;
	double y[3];
	if (m->foot[v].patch == F)
	{
		*u = (double)m->foot[v].a;
		*w = (double)m->foot[v].b;
		return;
	}
	*u = u0;
	*w = v0;
	mg_snap(mg_face(r->src, F), m->x + 3 * (size_t)v, u, w, 1, y);
}

/* a start on face F for a point of triangle f (its corners' footpoints on F, unwrapped across a seam) */
static void mg_rd_start(mg_rd const* r, unsigned int f, unsigned int F, double const* p, double* u, double* v)
{
	qaws_hm const* m = r->m;
	qaws_shell_point dv[3];
	unsigned int tv[3] = { 0, 1, 2 }, k;
	double w[3] = { 1, 1, 1 };
	for (k = 0; k < 3; k++)
	{
		qaws_mesh_location l = m->foot[m->cv[m->fc[f] + k]];
		dv[k].face = l.patch;
		dv[k].uv.x = l.a;
		dv[k].uv.y = l.b;
	}
	mg_start(mg_face(r->src, F), dv, tv, w, 3, F, p, u, v);
}

/* feature or border edge c split at its middle on its curve (the face's side or crease), the new vertex's footpoint
   there; its fans' edges flipped to Delaunay. HM_NONE on failure */
static unsigned int mg_rd_split_edge(mg_rd* r, unsigned int c)
{
	qaws_hm* m = r->m;
	unsigned int f = m->cf[c], F = mg_rd_face(r, f), a = m->cv[c], b = hm_to(m, c), v, q;
	qaws_surface const* s = mg_face(r->src, F);
	qaws_range ur = qaws_surface_get_u_range(s), vr = qaws_surface_get_v_range(s);
	double su = (double)ur.max_value - (double)ur.min_value, sv = (double)vr.max_value - (double)vr.min_value, u0, v0, ua, va, ub, vb, u, w, x[3];
	mg_rd_start(r, f, F, m->x + 3 * (size_t)a, &u0, &v0);
	mg_rd_uv(r, a, F, u0, v0, &ua, &va);
	mg_rd_uv(r, b, F, u0, v0, &ub, &vb);
	/* (across a seam: b within half a period of a) */
	if (ub - ua > 0.5 * su) ub -= su;
	if (ua - ub > 0.5 * su) ub += su;
	if (vb - va > 0.5 * sv) vb -= sv;
	if (va - vb > 0.5 * sv) vb += sv;
	u = 0.5 * (ua + ub);
	w = 0.5 * (va + vb);
	if (u < (double)ur.min_value) u += su;
	if (u > (double)ur.max_value) u -= su;
	if (w < (double)vr.min_value) w += sv;
	if (w > (double)vr.max_value) w -= sv;
	{
		double mid[3];
		int k;
		for (k = 0; k < 3; k++)
			mid[k] = 0.5 * (m->x[3 * (size_t)a + k] + m->x[3 * (size_t)b + k]);
		mg_snap(s, mid, &u, &w, 1, x);
	}
	if (!hm_reserve(m, 1, 12, 4) || (v = hm_split_edge(m, c, 0.5, x, 1)) == HM_NONE)
		return HM_NONE;
	m->foot[v].patch = F;
	m->foot[v].a = (qaws_scalar)u;
	m->foot[v].b = (qaws_scalar)w;
	for (q = m->vc[v]; q != HM_NONE; q = hm_round(m, v, q))
		if (!mg_rd_push(r, hm_next(m, q)))
			return HM_NONE;
	return mg_rd_flips(r) ? v : HM_NONE;
}

/* the triangle's corner whose edge is the first feature or border edge whose diametral ball holds y (within the
   triangle and its neighbours'), or HM_NONE */
static unsigned int mg_rd_encroached(qaws_hm const* m, unsigned int f, double const* y)
{
	unsigned int c, k;
	for (k = 0; k < 4; k++)
	{
		unsigned int g = k == 0 ? f : m->co[m->fc[f] + k - 1] == HM_NONE ? HM_NONE : m->cf[m->co[m->fc[f] + k - 1]];
		if (g == HM_NONE || m->fn[g] != 3)
			continue;
		for (c = m->fc[g]; c < m->fc[g] + 3; c++)
		{
			double const *a = m->x + 3 * (size_t)m->cv[c], *b = m->x + 3 * (size_t)hm_to(m, c);
			double mid[3], d2 = 0, r2 = 0;
			int q;
			if (!mg_rd_fixed(m, c))
				continue;
			for (q = 0; q < 3; q++)
			{
				mid[q] = 0.5 * (a[q] + b[q]);
				d2 += (y[q] - mid[q]) * (y[q] - mid[q]);
				r2 += 0.25 * (a[q] - b[q]) * (a[q] - b[q]);
			}
			if (d2 < r2 * (1 - 1e-9))
				return c;
		}
	}
	return HM_NONE;
}

/* the triangle holding y (projected on each triangle's plane), walked from f; *wall: a feature or border edge the walk
   would cross (HM_NONE when found) */
static unsigned int mg_rd_locate(qaws_hm const* m, unsigned int f, double const* y, unsigned int* wall)
{
	unsigned int steps, prev = HM_NONE;
	*wall = HM_NONE;
	for (steps = 0; steps < 1000; steps++)
	{
		unsigned int c0 = m->fc[f], k, worst = HM_NONE;
		double const* P[3];
		double nrm[3], e1[3], e2[3], bw[3], best = -1e-12, area;
		for (k = 0; k < 3; k++)
			P[k] = m->x + 3 * (size_t)m->cv[c0 + k];
		for (k = 0; k < 3; k++)
		{
			e1[k] = P[1][k] - P[0][k];
			e2[k] = P[2][k] - P[0][k];
		}
		nrm[0] = e1[1] * e2[2] - e1[2] * e2[1];
		nrm[1] = e1[2] * e2[0] - e1[0] * e2[2];
		nrm[2] = e1[0] * e2[1] - e1[1] * e2[0];
		area = mg_dot(nrm, nrm);
		if (!(area > 0))
			return HM_NONE;
		/* (barycentrics: the sub-triangles' signed areas; the edge opposite corner k is corner k + 1's) */
		for (k = 0; k < 3; k++)
		{
			double const *A = P[(k + 1) % 3], *B = P[(k + 2) % 3];
			double a1[3], a2[3], cr[3];
			int q;
			for (q = 0; q < 3; q++)
			{
				a1[q] = B[q] - A[q];
				a2[q] = y[q] - A[q];
			}
			cr[0] = a1[1] * a2[2] - a1[2] * a2[1];
			cr[1] = a1[2] * a2[0] - a1[0] * a2[2];
			cr[2] = a1[0] * a2[1] - a1[1] * a2[0];
			bw[k] = mg_dot(cr, nrm) / area;
			if (bw[k] < best)
			{
				best = bw[k];
				worst = c0 + (k + 1) % 3;
			}
		}
		if (worst == HM_NONE)
			return f;
		if (mg_rd_fixed(m, worst))
		{
			*wall = worst;
			return HM_NONE;
		}
		if (m->cf[m->co[worst]] == prev && steps > 2)
			return f;   /* (oscillating on a bend: this one) */
		prev = f;
		f = m->cf[m->co[worst]];
		if (m->fn[f] != 3)
			return HM_NONE;
	}
	return HM_NONE;
}

/* triangle f's state: 0 good, else bad (an angle under the bound away from two feature edges meeting, its circumradius
   over 0.8 of the size, its centre off the face by more than the chord); its circumcentre */
static int mg_rd_bad(mg_rd const* r, unsigned int f, double* cc)
{
	qaws_hm const* m = r->m;
	unsigned int c0 = m->fc[f], k;
	double const *A = m->x + 3 * (size_t)m->cv[c0], *B = m->x + 3 * (size_t)m->cv[c0 + 1], *C = m->x + 3 * (size_t)m->cv[c0 + 2];
	double u[3], v[3], n[3], uu, vv, nn, t1[3], t2[3], R2;
	int bad = 0;
	for (k = 0; k < 3; k++)
	{
		u[k] = B[k] - A[k];
		v[k] = C[k] - A[k];
	}
	n[0] = u[1] * v[2] - u[2] * v[1]; n[1] = u[2] * v[0] - u[0] * v[2]; n[2] = u[0] * v[1] - u[1] * v[0];
	uu = mg_dot(u, u);
	vv = mg_dot(v, v);
	nn = mg_dot(n, n);
	if (!(nn > 0))
		return 0;
	t1[0] = n[1] * u[2] - n[2] * u[1]; t1[1] = n[2] * u[0] - n[0] * u[2]; t1[2] = n[0] * u[1] - n[1] * u[0];
	t2[0] = v[1] * n[2] - v[2] * n[1]; t2[1] = v[2] * n[0] - v[0] * n[2]; t2[2] = v[0] * n[1] - v[1] * n[0];
	for (k = 0; k < 3; k++)
		cc[k] = A[k] + (vv * t1[k] + uu * t2[k]) / (2 * nn);
	R2 = (cc[0] - A[0]) * (cc[0] - A[0]) + (cc[1] - A[1]) * (cc[1] - A[1]) + (cc[2] - A[2]) * (cc[2] - A[2]);
	if (R2 > 0.64 * r->L * r->L)
		bad = 1;
	for (k = 0; k < 3 && !bad; k++)
	{
		unsigned int c = c0 + k;
		/* (an angle between two feature or border edges is the input's: kept) */
		if (mg_rd_cos(m, c) > r->cosmax && !(mg_rd_fixed(m, c) && mg_rd_fixed(m, hm_prev(m, c))))
			bad = 1;
	}
	if (!bad && r->chord > 0)
	{
		double g[3], uv0, vv0, y[3];
		unsigned int F = mg_rd_face(r, f);
		for (k = 0; k < 3; k++)
			g[k] = (A[k] + B[k] + C[k]) / 3;
		mg_rd_start(r, f, F, g, &uv0, &vv0);
		mg_snap(mg_face(r->src, F), g, &uv0, &vv0, 0, y);
		bad = (y[0] - g[0]) * (y[0] - g[0]) + (y[1] - g[1]) * (y[1] - g[1]) + (y[2] - g[2]) * (y[2] - g[2]) > r->chord * r->chord;
	}
	return bad;
}

/* RESTRICTED_DELAUNAY: the PARAMETRIC mesh made Delaunay (Lawson's flips, never across a feature or border edge) and
   refined as Chew (1993) and Ruppert do: a bad triangle's circumcentre snapped onto its face and inserted (its fans
   flipped back to Delaunay); one encroaching a feature or border edge's diametral ball, or behind one, splits that edge
   at its middle on its curve instead. The angle bound 30 degrees but between two feature edges, the size 0.8 L of
   circumradius, the chord tolerance at the centroid; at most 8 times the start's vertices added. */
static qaws_status mg_delaunay(mg_src const* src, unsigned int nf, qaws_generate_desc const* desc, double L, qaws_mesh_3d** out_mesh,
	qaws_mesh_report* rep)
{
	qaws_mesh_3d* out = NULL;
	mg_rd r;
	unsigned int f, c, added = 0, limit, pass;
	qaws_status st = mg_parametric(src, nf, desc, L, &out, rep);
	if (st != QAWS_STATUS_OK)
		return st;
	memset(&r, 0, sizeof r);
	r.src = src;
	r.m = &out->hm;
	r.L = L;
	r.chord = (double)desc->chord_tolerance;
	r.cosmax = cos(30 * 3.14159265358979323846 / 180);
	if (qaws_mesh_3d_find_attribute(&r.fa, out, "face") != QAWS_STATUS_OK)
	{
		qaws_mesh_3d_destroy(out);
		return QAWS_STATUS_INTERNAL_ERROR;
	}
	limit = 8 * r.m->nv + 64;
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	for (c = 0; c < r.m->nc; c++)
		if (r.m->cf[c] != HM_NONE && !mg_rd_push(&r, c))
			goto done;
	if (!mg_rd_flips(&r))
		goto done;
	for (pass = 0; pass < 64; pass++)
	{
		unsigned int nb = 0;
		for (f = 0; f < r.m->nf && added < limit; f++)
		{
			double cc[3], u, v, y[3];
			unsigned int F, g, wall = HM_NONE, enc, nv;
			if (r.m->fn[f] != 3 || r.m->fc[f] == HM_NONE || !mg_rd_bad(&r, f, cc))
				continue;
			nb++;
			F = mg_rd_face(&r, f);
			mg_rd_start(&r, f, F, cc, &u, &v);
			mg_snap(mg_face(src, F), cc, &u, &v, 0, y);
			enc = mg_rd_encroached(r.m, f, y);
			g = enc == HM_NONE ? mg_rd_locate(r.m, f, y, &wall) : HM_NONE;
			if (enc == HM_NONE && g == HM_NONE)
				enc = wall;
			if (enc != HM_NONE)
			{
				if (mg_rd_split_edge(&r, enc) == HM_NONE)
					goto done;
				added++;
				continue;
			}
			if (g == HM_NONE || mg_rd_face(&r, g) != F)
				continue;
			if (!hm_reserve(r.m, 1, 6, 2) || (nv = hm_split_face(r.m, g, y)) == HM_NONE)
				goto done;
			r.m->foot[nv].patch = F;
			r.m->foot[nv].a = (qaws_scalar)u;
			r.m->foot[nv].b = (qaws_scalar)v;
			for (c = r.m->vc[nv]; c != HM_NONE; c = hm_round(r.m, nv, c))
				if (!mg_rd_push(&r, hm_next(r.m, c)))
					goto done;
			if (!mg_rd_flips(&r))
				goto done;
			added++;
		}
		if (!nb || added >= limit)
			break;
	}
	hm_compact(r.m);
	qaws_mesh_3d_measure_quality(&rep->quality, out);
	*out_mesh = out;
	out = NULL;
	st = QAWS_STATUS_OK;
done:
	qaws_internal_free(r.stack);
	qaws_mesh_3d_destroy(out);
	return st;
}

/* STRUCTURED: per face an iso-parametric grid, its counts made the same across shared sides and opposite sides (classes of
   face directions joined), each class's count its sides' longest length over L; with tolerances (no target) the classes
   of faces whose cells miss them (the surface at a cell's centre off its corners' mean by more than the chord tolerance,
   corner normals turning more than the angle tolerance) doubled until met, not below min_length. Quads, or triangles by
   the shorter diagonal; the nodes welded where faces meet (shared sides, seams, poles), a quad of three nodes a
   triangle. */
static qaws_status mg_structured(mg_src const* src, unsigned int nf, qaws_generate_desc const* desc, double L, qaws_mesh_3d** out_mesh,
	qaws_mesh_report* rep)
{
	unsigned int *uf = NULL, *cnt = NULL, f, k, round, *fv = NULL, *fs = NULL, *map = NULL, *slot = NULL, nv = 0, nface = 0, nc = 0, mask;
	unsigned int *fface = NULL, total = 0, maxn = 0;
	double *dom = NULL, *len = NULL, *x = NULL, *uvbuf = NULL, *par = NULL, ext = 0, cell, lo[3] = { 1e300, 1e300, 1e300 }, hi[3] = { -1e300, -1e300, -1e300 };
	qaws_mesh_location* foot = NULL;
	int adaptive = desc->target_length == 0 && (desc->chord_tolerance > 0 || desc->angle_tolerance > 0);
	qaws_mesh_3d* out = NULL;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	uf = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 2, (size_t)nf + 1));
	cnt = (unsigned int*)qaws_internal_scratch_zalloc(2 * (size_t)nf + 2, sizeof(unsigned int));
	dom = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 4, (size_t)nf + 1));
	len = (double*)qaws_internal_scratch_zalloc(4 * (size_t)nf + 4, sizeof(double));
	if (!uf || !cnt || !dom || !len)
		goto done;
	/* the domains, the sides' lengths, the classes */
	for (f = 0; f < nf; f++)
	{
		qaws_surface const* s = mg_face(src, f);
		qaws_range ur = qaws_surface_get_u_range(s), vr = qaws_surface_get_v_range(s);
		unsigned int nb[4];
		int nbs[4];
		dom[4 * (size_t)f] = (double)ur.min_value;
		dom[4 * (size_t)f + 1] = (double)ur.max_value;
		dom[4 * (size_t)f + 2] = (double)vr.min_value;
		dom[4 * (size_t)f + 3] = (double)vr.max_value;
		for (k = 0; k < 4; k++)
			len[4 * (size_t)f + k] = mg_side_params(s, dom + 4 * (size_t)f, (int)k, 1, NULL);
		uf[2 * f] = 2 * f;
		uf[2 * f + 1] = 2 * f + 1;
		if (src->shell && !sh_face_grid_sides(src->shell, f, nb, nbs))
		{
			st = QAWS_STATUS_UNSUPPORTED_OPERATION;   /* (a trimmed face: no grid) */
			goto done;
		}
		(void)nbs;
	}
	for (f = 0; src->shell && f < nf; f++)
	{
		unsigned int nb[4];
		int nbs[4];
		sh_face_grid_sides(src->shell, f, nb, nbs);
		for (k = 0; k < 4; k++)
			if (nb[k] < nf && nbs[k] >= 0 && nbs[k] < 4)
			{
				unsigned int ra = mg_root(uf, 2 * f + (k & 1)), rb = mg_root(uf, 2 * nb[k] + ((unsigned int)nbs[k] & 1));
				if (ra != rb)
					uf[ra] = rb;
			}
	}
	/* (a class's count: its sides' longest over L; sides 0 and 2 count along u, 1 and 3 along v) */
	for (f = 0; f < nf; f++)
		for (k = 0; k < 4; k++)
		{
			unsigned int r = mg_root(uf, 2 * f + (k & 1)), n = (unsigned int)fmin(4096, fmax(1, ceil(len[4 * (size_t)f + k] / L - 1e-9)));
			cnt[r] = n > cnt[r] ? n : cnt[r];
		}
	/* tolerances: classes of faces missing them doubled */
	for (round = 0; adaptive && round < 8; round++)
	{
		unsigned int grow = 0;
		for (f = 0; f < nf; f++)
		{
			qaws_surface const* s = mg_face(src, f);
			unsigned int ru = mg_root(uf, 2 * f), rw = mg_root(uf, 2 * f + 1), nu = cnt[ru], nw = cnt[rw], i, j, bad = 0;
			double *a, *b, *c, *d, *uv, step;
			step = fmax(len[4 * (size_t)f], len[4 * (size_t)f + 2]) / nu;
			step = fmax(step, fmax(len[4 * (size_t)f + 1], len[4 * (size_t)f + 3]) / nw);
			if (desc->min_length > 0 && step < 2 * (double)desc->min_length)
				continue;
			par = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), 2 * ((size_t)nu + nw + 2)));
			uvbuf = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 2, ((size_t)nu + 1) * ((size_t)nw + 1)));
			if (!par || !uvbuf)
				goto done;
			a = par;
			b = a + nu + 1;
			c = b + nu + 1;
			d = c + nw + 1;
			uv = uvbuf;
			mg_grid_uv(s, dom + 4 * (size_t)f, nu, nw, a, b, c, d, uv);
			for (j = 0; j < nw && bad != 3; j++)
				for (i = 0; i < nu && bad != 3; i++)
				{
					unsigned int q[4] = { i + (nu + 1) * j, i + 1 + (nu + 1) * j, i + 1 + (nu + 1) * (j + 1), i + (nu + 1) * (j + 1) }, z;
					double mid[3] = { 0, 0, 0 }, cu = 0, cv = 0, p[3], nrm[4][3];
					for (z = 0; z < 4; z++)
					{
						qaws_surface_eval_result e;
						memset(&e, 0, sizeof e);
						qaws_surface_evaluate(&e, s, (qaws_scalar)uv[2 * q[z]], (qaws_scalar)uv[2 * q[z] + 1], QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL);
						mg_v3(e.position, p);
						mg_v3(e.normal, nrm[z]);
						for (k = 0; k < 3; k++)
							mid[k] += 0.25 * p[k];
						cu += 0.25 * uv[2 * q[z]];
						cv += 0.25 * uv[2 * q[z] + 1];
					}
					/* each side's middle against its chord (sides 0-1 and 3-2 along u, 1-2 and 0-3 along v), its normals' turn: the
					   direction that misses doubled; the centre missing alone (a twist): both */
					for (z = 0; z < 4; z++)
					{
						unsigned int z1 = (z + 1) % 4, dir = z & 1 ? 2u : 1u;
						double A[3], B[3], m3[3], e2 = 0;
						mg_at(s, uv[2 * q[z]], uv[2 * q[z] + 1], A);
						mg_at(s, uv[2 * q[z1]], uv[2 * q[z1] + 1], B);
						mg_at(s, 0.5 * (uv[2 * q[z]] + uv[2 * q[z1]]), 0.5 * (uv[2 * q[z] + 1] + uv[2 * q[z1] + 1]), m3);
						for (k = 0; k < 3; k++)
							e2 += (m3[k] - 0.5 * (A[k] + B[k])) * (m3[k] - 0.5 * (A[k] + B[k]));
						if (desc->chord_tolerance > 0 && sqrt(e2) > (double)desc->chord_tolerance)
							bad |= dir;
						if (desc->angle_tolerance > 0)
						{
							double dp = mg_dot(nrm[z], nrm[z1]), l2 = mg_dot(nrm[z], nrm[z]) * mg_dot(nrm[z1], nrm[z1]);
							if (l2 > 0 && dp < cos((double)desc->angle_tolerance) * sqrt(l2))
								bad |= dir;
						}
					}
					mg_at(s, cu, cv, p);
					if (!bad && desc->chord_tolerance > 0 && sqrt((p[0] - mid[0]) * (p[0] - mid[0]) + (p[1] - mid[1]) * (p[1] - mid[1]) + (p[2] - mid[2]) * (p[2] - mid[2])) > (double)desc->chord_tolerance)
						bad = 3;
				}
			qaws_internal_free(par);
			qaws_internal_free(uvbuf);
			par = uvbuf = NULL;
			if ((bad & 1) && cnt[ru] < 4096)
				cnt[ru] *= 2;
			if ((bad & 2) && (rw != ru || !(bad & 1)) && cnt[rw] < 4096)
				cnt[rw] *= 2;
			grow += bad != 0;
		}
		if (!grow)
			break;
	}
	/* the nodes, welded by position */
	for (f = 0; f < nf; f++)
	{
		unsigned int nu = cnt[mg_root(uf, 2 * f)], nw = cnt[mg_root(uf, 2 * f + 1)];
		total += (nu + 1) * (nw + 1);
		nc += 6 * nu * nw;   /* (triangles: two a cell) */
		maxn = nu + nw > maxn ? nu + nw : maxn;
	}
	x = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)total + 1));
	foot = (qaws_mesh_location*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_mesh_location), (size_t)total + 1));
	map = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)total + 1));
	fv = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)nc + 1));
	fs = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)nc / 2 + 1));
	fface = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)nc / 2 + 1));
	par = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), 2 * ((size_t)maxn + 2)));
	for (mask = 1; mask < 2 * total + 2; mask <<= 1) {}
	slot = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), mask));
	if (!x || !foot || !map || !fv || !fs || !fface || !par || !slot)
		goto done;
	memset(slot, 0xFF, sizeof(unsigned int) * mask);
	mask -= 1;
	{
		/* (the weld's cell: a millionth of the shell's size, from its corners) */
		for (f = 0; f < nf; f++)
			for (k = 0; k < 4; k++)
			{
				double p[3], uv[2];
				unsigned int q;
				mg_side_uv(dom + 4 * (size_t)f, (int)k, 0, uv);
				mg_at(mg_face(src, f), uv[0], uv[1], p);
				for (q = 0; q < 3; q++)
				{
					lo[q] = fmin(lo[q], p[q]);
					hi[q] = fmax(hi[q], p[q]);
				}
			}
		for (k = 0; k < nf; k++)
			ext = fmax(ext, len[4 * (size_t)k] + len[4 * (size_t)k + 1]);
		ext = fmax(ext, sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) + (hi[2] - lo[2]) * (hi[2] - lo[2])));
		cell = 1e-7 * (ext > 0 ? ext : 1);
	}
	total = 0;
	nc = 0;
	for (f = 0; f < nf; f++)
	{
		qaws_surface const* s = mg_face(src, f);
		unsigned int nu = cnt[mg_root(uf, 2 * f)], nw = cnt[mg_root(uf, 2 * f + 1)], i, j, base = total;
		double *uv = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 2, ((size_t)nu + 1) * ((size_t)nw + 1)));
		if (!uv)
			goto done;
		mg_grid_uv(s, dom + 4 * (size_t)f, nu, nw, par, par + nu + 1, par + 2 * (nu + 1), par + 2 * (nu + 1) + nw + 1, uv);
		for (j = 0; j <= nw; j++)
			for (i = 0; i <= nu; i++)
			{
				double p[3], u = uv[2 * (i + (size_t)(nu + 1) * j)], v = uv[2 * (i + (size_t)(nu + 1) * j) + 1];
				mg_at(s, u, v, p);
				map[base + i + (nu + 1) * j] = mg_weld_point(x, foot, &nv, slot, mask, cell, p, f, u, v, 1);
			}
		total += (nu + 1) * (nw + 1);
		/* the cells, counter-clockwise in (u, v) */
		for (j = 0; j < nw; j++)
			for (i = 0; i < nu; i++)
			{
				unsigned int q4[4], q[4], m = 0, z;
				q4[0] = map[base + i + (nu + 1) * j];
				q4[1] = map[base + i + 1 + (nu + 1) * j];
				q4[2] = map[base + i + 1 + (nu + 1) * (j + 1)];
				q4[3] = map[base + i + (nu + 1) * (j + 1)];
				for (z = 0; z < 4; z++)
					if (q4[z] != q4[(z + 3) % 4])
						q[m++] = q4[z];
				if (m == 4 && q[0] == q[2])
					continue;
				if (m < 3 || (m == 3 && (q[0] == q[2])))
					continue;
				if (m == 4 && desc->elements == QAWS_ELEMENTS_TRIANGLES)
				{
					double const *A = x + 3 * (size_t)q[0], *B = x + 3 * (size_t)q[1], *C = x + 3 * (size_t)q[2], *D = x + 3 * (size_t)q[3];
					double ac = (A[0] - C[0]) * (A[0] - C[0]) + (A[1] - C[1]) * (A[1] - C[1]) + (A[2] - C[2]) * (A[2] - C[2]);
					double bd = (B[0] - D[0]) * (B[0] - D[0]) + (B[1] - D[1]) * (B[1] - D[1]) + (B[2] - D[2]) * (B[2] - D[2]);
					unsigned int t6[6];
					if (ac <= bd)
					{
						t6[0] = q[0]; t6[1] = q[1]; t6[2] = q[2]; t6[3] = q[0]; t6[4] = q[2]; t6[5] = q[3];
					}
					else
					{
						t6[0] = q[0]; t6[1] = q[1]; t6[2] = q[3]; t6[3] = q[1]; t6[4] = q[2]; t6[5] = q[3];
					}
					memcpy(fv + nc, t6, sizeof t6);
					nc += 6;
					fs[nface] = 3;
					fface[nface++] = f;
					fs[nface] = 3;
					fface[nface++] = f;
					continue;
				}
				memcpy(fv + nc, q, sizeof(unsigned int) * m);
				nc += m;
				fs[nface] = m;
				fface[nface++] = f;
			}
		qaws_internal_free(uv);
	}
	st = mg_assemble(x, foot, nv, fv, fs, fface, nface, out_mesh, rep);
done:
	qaws_mesh_3d_destroy(out);
	qaws_internal_free(uf); qaws_internal_free(cnt); qaws_internal_free(dom); qaws_internal_free(len); qaws_internal_free(x); qaws_internal_free(foot);
	qaws_internal_free(map); qaws_internal_free(fv); qaws_internal_free(fs); qaws_internal_free(fface); qaws_internal_free(par); qaws_internal_free(uvbuf);
	qaws_internal_free(slot);
	return st;
}

qaws_status qaws_mesh_generate_surface(qaws_mesh_3d** out_mesh, qaws_mesh_report* out_report, qaws_surface const* surface,
	qaws_generate_desc const* desc)
{
	qaws_mesh_report rep;
	qaws_shell_point* dvp = NULL;
	double *x = NULL, L, area = 0;
	unsigned int *tri = NULL, *tface = NULL, nv = 0, nt = 0;
	mg_src src;
	qaws_status st;
	if (!out_mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_mesh = NULL;
	if (!surface || !mg_check(desc))
		return QAWS_STATUS_INVALID_ARGUMENT;
	if ((desc->method == QAWS_GENERATE_PARAMETRIC || desc->method == QAWS_GENERATE_RESTRICTED_DELAUNAY) && desc->elements != QAWS_ELEMENTS_TRIANGLES)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	memset(&rep, 0, sizeof rep);
	/* no length given: a 30th of the square root of the area, from a coarse pass */
	L = desc->target_length > 0 ? (double)desc->target_length : (double)desc->max_length;
	if (!(L > 0))
	{
		st = mg_grid(surface, 1e300, &x, &dvp, &nv, &tri, &nt, &area);
		qaws_internal_free(x);
		qaws_internal_free(dvp);
		qaws_internal_free(tri);
		x = NULL;
		dvp = NULL;
		tri = NULL;
		if (st != QAWS_STATUS_OK)
			return st;
		L = sqrt(area) / 30;
		if (!(L > 0))
			return QAWS_STATUS_INVALID_ARGUMENT;
	}
	src.shell = NULL;
	src.one = surface;
	if (desc->method != QAWS_GENERATE_REMESH)
	{
		st = desc->method == QAWS_GENERATE_PARAMETRIC ? mg_parametric(&src, 1, desc, L, out_mesh, &rep)
			: desc->method == QAWS_GENERATE_RESTRICTED_DELAUNAY ? mg_delaunay(&src, 1, desc, L, out_mesh, &rep) : mg_structured(&src, 1, desc, L, out_mesh, &rep);
		if (st == QAWS_STATUS_OK && out_report)
			*out_report = rep;
		return st;
	}
	st = mg_grid(surface, mg_dense_length(desc, L), &x, &dvp, &nv, &tri, &nt, &area);
	if (st != QAWS_STATUS_OK)
		return st;
	tface = (unsigned int*)qaws_internal_scratch_zalloc((size_t)nt + 1, sizeof(unsigned int));
	if (!tface)
		st = QAWS_STATUS_ALLOCATION_FAILURE;
	src.shell = NULL;
	src.one = surface;
	if (st == QAWS_STATUS_OK)
		st = mg_finish(&src, x, dvp, nv, tri, tface, nt, desc, L, out_mesh, &rep);
	qaws_internal_free(x);
	qaws_internal_free(dvp);
	qaws_internal_free(tri);
	qaws_internal_free(tface);
	if (st == QAWS_STATUS_OK && out_report)
		*out_report = rep;
	return st;
}

/* ---- curved elements: nodes of order 2 and 3 ------------------------- */

/* nodes inside a face of k corners at order o: (o-1)(o-2)/2 in a triangle, (o-1)^2 in a quad, none in a polygon */
static unsigned int mg_inner(unsigned int k, unsigned int o)
{
	return k == 3 ? (o - 1) * (o - 2) / 2 : k == 4 ? (o - 1) * (o - 1) : 0;
}

qaws_status qaws_mesh_generate_node_counts(unsigned int* out_node_count, unsigned int* out_element_node_count, qaws_mesh_3d const* mesh,
	unsigned int order)
{
	qaws_hm const* m;
	unsigned int f, n, e = 0;
	if (!mesh || order < 1 || order > 3)
		return QAWS_STATUS_INVALID_ARGUMENT;
	m = &mesh->hm;
	n = m->nv + hm_edge_count(m) * (order - 1);
	for (f = 0; f < m->nf; f++)
	{
		n += mg_inner(m->fn[f], order);
		e += m->fn[f] * order + mg_inner(m->fn[f], order);
	}
	if (out_node_count)
		*out_node_count = n;
	if (out_element_node_count)
		*out_element_node_count = e;
	return QAWS_STATUS_OK;
}

/* the node at weights w of vertices vi (n of them): their straight point snapped onto the element's face F,
   started from their (u, v) on F; on: held on the domain's side (a feature or border edge) */
static void mg_node(mg_src const* src, qaws_hm const* m, unsigned int F, unsigned int const* vi, double const* w, unsigned int n, int on,
	qaws_scalar* out_x, qaws_mesh_location* out_loc)
{
	qaws_shell_point cp[4];
	unsigned int tv[4] = { 0, 1, 2, 3 }, k;
	double p[3] = { 0, 0, 0 }, x[3], u = 0, v = 0;
	qaws_surface const* s = mg_face(src, F);
	qaws_mesh_location loc;
	loc.patch = QAWS_MESH_NONE;
	loc.a = loc.b = 0;
	for (k = 0; k < n; k++)
	{
		p[0] += w[k] * m->x[3 * (size_t)vi[k]];
		p[1] += w[k] * m->x[3 * (size_t)vi[k] + 1];
		p[2] += w[k] * m->x[3 * (size_t)vi[k] + 2];
		memset(&cp[k], 0, sizeof cp[k]);
		cp[k].face = m->foot ? m->foot[vi[k]].patch : QAWS_MESH_NONE;
		cp[k].uv.x = m->foot ? m->foot[vi[k]].a : 0;
		cp[k].uv.y = m->foot ? m->foot[vi[k]].b : 0;
	}
	memcpy(x, p, sizeof p);
	if (s)
	{
		mg_start(s, cp, tv, w, n, F, p, &u, &v);
		mg_snap(s, p, &u, &v, on, x);
		loc.patch = F;
		loc.a = (qaws_scalar)u;
		loc.b = (qaws_scalar)v;
	}
	for (k = 0; k < 3; k++)
		out_x[k] = (qaws_scalar)x[k];
	if (out_loc)
		*out_loc = loc;
}

static qaws_status mg_nodes(qaws_scalar* out_positions, qaws_mesh_location* out_locations, unsigned int* out_element_nodes,
	qaws_mesh_3d const* mesh, mg_src const* src, unsigned int o)
{
	qaws_hm const* m;
	unsigned int *base = NULL, fa = HM_NONE, f, c, k, i, j, next, e = 0;
	unsigned char* rev = NULL;
	if (!mesh || !out_positions || !out_element_nodes || o < 1 || o > 3)
		return QAWS_STATUS_INVALID_ARGUMENT;
	m = &mesh->hm;
	if (qaws_mesh_3d_find_attribute(&fa, mesh, "face") != QAWS_STATUS_OK)
	{
		if (src->shell)
			return QAWS_STATUS_INVALID_ARGUMENT;
		fa = HM_NONE;
	}
	/* per corner: its edge's first node, and whether the edge was numbered from its other end */
	base = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)m->nc + 1));
	rev = (unsigned char*)qaws_internal_scratch_zalloc((size_t)m->nc + 1, 1);
	if (!base || !rev)
	{
		qaws_internal_free(base);
		qaws_internal_free(rev);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (c = 0; c < m->nc; c++)
		base[c] = HM_NONE;
	for (i = 0; i < m->nv; i++)
	{
		for (k = 0; k < 3; k++)
			out_positions[3 * (size_t)i + k] = (qaws_scalar)m->x[3 * (size_t)i + k];
		if (out_locations)
		{
			out_locations[i].patch = m->foot ? m->foot[i].patch : QAWS_MESH_NONE;
			out_locations[i].a = m->foot ? m->foot[i].a : 0;
			out_locations[i].b = m->foot ? m->foot[i].b : 0;
		}
	}
	next = m->nv;
	for (f = 0; f < m->nf; f++)
	{
		unsigned int n = m->fn[f], c0 = m->fc[f], F = fa != HM_NONE ? (unsigned int)hm_attr_at(m, fa, f)[0] : 0;
		for (k = 0; k < n; k++)
			out_element_nodes[e++] = m->cv[c0 + k];
		/* edge nodes, from corner k toward corner k + 1 */
		for (k = 0; o > 1 && k < n; k++)
		{
			c = c0 + k;
			if (base[c] == HM_NONE)
			{
				unsigned int ab[2];
				int on = m->co[c] == HM_NONE || (m->cflag[c] & HM_FEATURE);
				ab[0] = m->cv[c];
				ab[1] = hm_to(m, c);
				base[c] = next;
				if (m->co[c] != HM_NONE)
				{
					base[m->co[c]] = next;
					rev[m->co[c]] = 1;
				}
				for (j = 1; j < o; j++, next++)
				{
					double w[2];
					w[1] = (double)j / o;
					w[0] = 1 - w[1];
					mg_node(src, m, F, ab, w, 2, on, out_positions + 3 * (size_t)next, out_locations ? out_locations + next : NULL);
				}
			}
			for (j = 1; j < o; j++)
				out_element_nodes[e++] = base[c] + (rev[c] ? o - 1 - j : j - 1);
		}
		/* inner nodes: a cubic triangle's centroid; a quad's at (i/o, j/o) along its sides 0-1 and 0-3, i fastest */
		if (n == 3 && o == 3)
		{
			double w[3] = { 1.0 / 3, 1.0 / 3, 1.0 / 3 };
			mg_node(src, m, F, m->cv + c0, w, 3, 0, out_positions + 3 * (size_t)next, out_locations ? out_locations + next : NULL);
			out_element_nodes[e++] = next++;
		}
		if (n == 4)
			for (j = 1; j < o; j++)
				for (i = 1; i < o; i++)
				{
					double s = (double)i / o, t = (double)j / o, w[4];
					w[0] = (1 - s) * (1 - t);
					w[1] = s * (1 - t);
					w[2] = s * t;
					w[3] = (1 - s) * t;
					mg_node(src, m, F, m->cv + c0, w, 4, 0, out_positions + 3 * (size_t)next, out_locations ? out_locations + next : NULL);
					out_element_nodes[e++] = next++;
				}
	}
	qaws_internal_free(base);
	qaws_internal_free(rev);
	return QAWS_STATUS_OK;
}

qaws_status qaws_mesh_generate_nodes_surface(qaws_scalar* out_positions, qaws_mesh_location* out_locations, unsigned int* out_element_nodes,
	qaws_mesh_3d const* mesh, qaws_surface const* surface, unsigned int order)
{
	mg_src src;
	if (!surface)
		return QAWS_STATUS_INVALID_ARGUMENT;
	src.shell = NULL;
	src.one = surface;
	return mg_nodes(out_positions, out_locations, out_element_nodes, mesh, &src, order);
}

qaws_status qaws_mesh_generate_nodes_shell(qaws_scalar* out_positions, qaws_mesh_location* out_locations, unsigned int* out_element_nodes,
	qaws_mesh_3d const* mesh, qaws_shell_atlas const* shell, unsigned int order)
{
	mg_src src;
	if (!shell)
		return QAWS_STATUS_INVALID_ARGUMENT;
	src.shell = shell;
	src.one = NULL;
	return mg_nodes(out_positions, out_locations, out_element_nodes, mesh, &src, order);
}
