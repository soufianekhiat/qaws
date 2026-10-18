/*
 * Meshes of exact surfaces and shells: a dense conformal mesh of the whole
 * shell (its faces' shared edges and free sides as features), remeshed or
 * quadrangulated on its PN reference, then every vertex snapped onto its
 * face's surface (Newton from the dense mesh's (u, v)).
 */

#include "qaws_mesh_internal.h"
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

/* the closest point of s to p by Gauss-Newton from (u, v), kept in the domain; a start on a side of the
   domain (the border, or a face's edge) kept on that side
   shortcut: a trimmed face's vertices on its trim curve are snapped onto the surface, not the curve (off it
   by the dense mesh's error), the trim curve to be followed when trimmed shells need it */
static void mg_snap(qaws_surface const* s, double const* p, double* u, double* v, double* x)
{
	qaws_range ur = qaws_surface_get_u_range(s), vr = qaws_surface_get_v_range(s);
	double lo[2] = { (double)ur.min_value, (double)vr.min_value }, hi[2] = { (double)ur.max_value, (double)vr.max_value };
	double* uv[2];
	int it, fixed[2], j;
	uv[0] = u;
	uv[1] = v;
	for (j = 0; j < 2; j++)
	{
		double tol = 1e-9 * (hi[j] - lo[j]);
		fixed[j] = fabs(*uv[j] - lo[j]) <= tol || fabs(*uv[j] - hi[j]) <= tol;
		if (fixed[j])
			*uv[j] = fabs(*uv[j] - lo[j]) <= tol ? lo[j] : hi[j];
	}
	for (it = 0; it < 30 && !(fixed[0] && fixed[1]); it++)
	{
		qaws_surface_eval_result e;
		double S[3], Su[3], Sv[3], r[3], a, b, c, g0, g1, det, du, dv;
		unsigned int k;
		memset(&e, 0, sizeof e);
		if (qaws_surface_evaluate(s, (qaws_scalar)*u, (qaws_scalar)*v, QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV, &e) != QAWS_STATUS_OK)
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
	{
		qaws_surface_eval_result e;
		memset(&e, 0, sizeof e);
		qaws_surface_evaluate(s, (qaws_scalar)*u, (qaws_scalar)*v, QAWS_SURFACE_EVAL_POSITION, &e);
		mg_v3(e.position, x);
	}
}

/* a start for the snap onto face F from reference triangle t's corners (those on F, their u and v
   unwrapped across a seam); none on F: the nearest point of a coarse grid of F's domain */
static void mg_start(qaws_surface const* s, qaws_shell_point const* dv, unsigned int const* tv, double const* w, unsigned int F,
	double const* p, double* u, double* v)
{
	qaws_range ur = qaws_surface_get_u_range(s), vr = qaws_surface_get_v_range(s);
	double su = (double)ur.max_value - (double)ur.min_value, sv = (double)vr.max_value - (double)vr.min_value, wu = 0, ref_u = 0, ref_v = 0, acc_u = 0, acc_v = 0;
	unsigned int k, n = 0;
	for (k = 0; k < 3; k++)
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
				qaws_surface_evaluate(s, (qaws_scalar)uu, (qaws_scalar)vv, QAWS_SURFACE_EVAL_POSITION, &e);
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
		qaws_scalar* fv = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (nt + 1));
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
		st = qaws_mesh_3d_add_attribute(dense, &ad, fv, &fa);
		free(fv);
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
	st = qaws_mesh_reference_create(dense, NULL, &ref);
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
		st = qaws_size_field_create_curvature(ref, &cd, &sf);
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
		st = qaws_mesh_remesh(dense, &rd, &out, rep);
	}
	else
	{
		qaws_quadrangulate_desc qd;
		memset(&qd, 0, sizeof qd);
		qd.reference = ref;
		qd.target_length = (qaws_scalar)L;
		qd.pure = desc->elements == QAWS_ELEMENTS_QUADS;
		st = qaws_mesh_quadrangulate(dense, &qd, &out, rep);
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
			mg_start(s, dvp, tv, w, F, m->x + 3 * (size_t)i, &u, &v);
			mg_snap(s, m->x + 3 * (size_t)i, &u, &v, y);
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
		qaws_mesh_3d_measure_quality(out, &rep->quality);
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

/* the dense mesh's size: half the target, finer for an adaptive one with a small minimum */
static double mg_dense_length(qaws_generate_desc const* desc, double L)
{
	double Ld = 0.5 * L;
	if (desc->target_length == 0 && (desc->chord_tolerance > 0 || desc->angle_tolerance > 0) && desc->min_length > 0)
		Ld = fmax(0.5 * (double)desc->min_length, L / 8);
	return Ld;
}

qaws_status qaws_mesh_generate_shell(qaws_shell_atlas const* shell, qaws_generate_desc const* desc, qaws_mesh_3d** out_mesh,
	qaws_mesh_report* out_report)
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
	if (desc->method != QAWS_GENERATE_REMESH)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	memset(&rep, 0, sizeof rep);
	st = qaws_shell_region_create_whole(shell, &whole);
	if (st != QAWS_STATUS_OK)
		return st;
	qaws_shell_region_get_area(whole, &area);
	L = desc->target_length > 0 ? (double)desc->target_length : desc->max_length > 0 ? (double)desc->max_length : sqrt((double)area) / 30;
	if (!(L > 0))
	{
		st = QAWS_STATUS_INVALID_ARGUMENT;
		goto done;
	}
	/* the dense conformal mesh of the whole shell */
	memset(&smd, 0, sizeof smd);
	smd.edge_length = (qaws_scalar)mg_dense_length(desc, L);
	st = qaws_shell_region_triangulate(whole, &smd, &sm);
	if (st != QAWS_STATUS_OK)
		goto done;
	nv = qaws_shell_mesh_get_vertex_count(sm);
	nt = qaws_shell_mesh_get_triangle_count(sm);
	x = (double*)malloc(sizeof(double) * 3 * (nv + 1));
	dvp = (qaws_shell_point*)malloc(sizeof(qaws_shell_point) * (nv + 1));
	tri = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (nt + 1));
	tface = (unsigned int*)malloc(sizeof(unsigned int) * (nt + 1));
	if (!x || !dvp || !tri || !tface)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	for (i = 0; i < nv; i++)
	{
		qaws_vec3 p;
		qaws_shell_mesh_get_vertex(sm, i, &dvp[i], &p);
		mg_v3(p, x + 3 * (size_t)i);
	}
	for (i = 0; i < nt; i++)
	{
		qaws_shell_mesh_get_triangle(sm, i, tri + 3 * (size_t)i);
		tface[i] = qaws_shell_mesh_get_triangle_face(sm, i);
	}
	src.shell = shell;
	src.one = NULL;
	st = mg_finish(&src, x, dvp, nv, tri, tface, nt, desc, L, out_mesh, &rep);
done:
	free(x); free(dvp); free(tri); free(tface);
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
			qaws_surface_evaluate(s, (qaws_scalar)(u0 + (u1 - u0) * i / G), (qaws_scalar)(v0 + (v1 - v0) * j / G), QAWS_SURFACE_EVAL_POSITION, &e);
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
	x = (double*)malloc(sizeof(double) * 3 * n);
	dvp = (qaws_shell_point*)malloc(sizeof(qaws_shell_point) * n);
	map = (unsigned int*)malloc(sizeof(unsigned int) * n);
	tri = (unsigned int*)malloc(sizeof(unsigned int) * 6 * nu * nw);
	for (mask = 1; mask < 2 * n; mask <<= 1) {}
	slot = (unsigned int*)malloc(sizeof(unsigned int) * mask);
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
			qaws_surface_evaluate(s, (qaws_scalar)u, (qaws_scalar)v, QAWS_SURFACE_EVAL_POSITION, &e);
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
	free(map);
	free(slot);
	if (st != QAWS_STATUS_OK)
	{
		free(x);
		free(dvp);
		free(tri);
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

qaws_status qaws_mesh_generate_surface(qaws_surface const* surface, qaws_generate_desc const* desc, qaws_mesh_3d** out_mesh,
	qaws_mesh_report* out_report)
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
	if (desc->method != QAWS_GENERATE_REMESH)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	memset(&rep, 0, sizeof rep);
	/* no length given: a 30th of the square root of the area, from a coarse pass */
	L = desc->target_length > 0 ? (double)desc->target_length : (double)desc->max_length;
	if (!(L > 0))
	{
		st = mg_grid(surface, 1e300, &x, &dvp, &nv, &tri, &nt, &area);
		free(x);
		free(dvp);
		free(tri);
		x = NULL;
		dvp = NULL;
		tri = NULL;
		if (st != QAWS_STATUS_OK)
			return st;
		L = sqrt(area) / 30;
		if (!(L > 0))
			return QAWS_STATUS_INVALID_ARGUMENT;
	}
	st = mg_grid(surface, mg_dense_length(desc, L), &x, &dvp, &nv, &tri, &nt, &area);
	if (st != QAWS_STATUS_OK)
		return st;
	tface = (unsigned int*)calloc((size_t)nt + 1, sizeof(unsigned int));
	if (!tface)
		st = QAWS_STATUS_ALLOCATION_FAILURE;
	src.shell = NULL;
	src.one = surface;
	if (st == QAWS_STATUS_OK)
		st = mg_finish(&src, x, dvp, nv, tri, tface, nt, desc, L, out_mesh, &rep);
	free(x);
	free(dvp);
	free(tri);
	free(tface);
	if (st == QAWS_STATUS_OK && out_report)
		*out_report = rep;
	return st;
}
