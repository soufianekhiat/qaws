/*
 * Boundary layers (qaws_boundary_layers.h): prisms and hexahedra grown from a
 * surface's triangles and quads along one direction per vertex, the stacks
 * thinned where they would meet the surface or each other.
 *
 * shortcut: the distance to the surface along a vertex's direction is found
 * against every face (each face's box tested first), O(vertices x faces); a
 * bounding volume hierarchy when meshes past 10^5 faces are layered.
 */

#include "qaws_mesh_internal.h"
#include "../internal/qaws_internal_memory.h"
#include "../qaws_boundary_layers.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#define BL_SEE 0.2          /* the least cosine between a vertex's direction and a face round it */
#define BL_REACH 0.4        /* a stack's share of the distance to the surface along its direction */
#define BL_GRADE 1.25       /* neighbours' thicknesses within this ratio */
#define BL_ROUNDS 30        /* halving rounds for cells turned over */

static void bl_sub(double const* a, double const* b, double* o) { o[0] = a[0] - b[0]; o[1] = a[1] - b[1]; o[2] = a[2] - b[2]; }
static double bl_dot(double const* a, double const* b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
static void bl_cross(double const* a, double const* b, double* o)
{
	o[0] = a[1] * b[2] - a[2] * b[1];
	o[1] = a[2] * b[0] - a[0] * b[2];
	o[2] = a[0] * b[1] - a[1] * b[0];
}
static double bl_norm(double* a)
{
	double l = sqrt(bl_dot(a, a));
	if (l > 0) { a[0] /= l; a[1] /= l; a[2] /= l; }
	return l;
}

/* det(b - a, c - a, d - a) */
static double bl_det(double const* a, double const* b, double const* c, double const* d)
{
	double u[3], v[3], w[3], n[3];
	bl_sub(b, a, u); bl_sub(c, a, v); bl_sub(d, a, w);
	bl_cross(u, v, n);
	return bl_dot(n, w);
}

/* every corner's Jacobian positive: a ring of n (3 or 4) below, the same ring above, the base turning toward the top */
static int bl_cell_ok(double const* x, unsigned int const* v, unsigned int n)
{
	unsigned int i;
	for (i = 0; i < n; i++)
	{
		unsigned int j = (i + 1) % n, k = (i + n - 1) % n;
		if (!(bl_det(x + 3 * (size_t)v[i], x + 3 * (size_t)v[j], x + 3 * (size_t)v[k], x + 3 * (size_t)v[n + i]) > 0))
			return 0;
		/* above: the ring seen from below turns the other way, so its next and previous swap */
		if (!(bl_det(x + 3 * (size_t)v[n + i], x + 3 * (size_t)v[n + k], x + 3 * (size_t)v[n + j], x + 3 * (size_t)v[i]) > 0))
			return 0;
	}
	return 1;
}

/* the nearest t > 0 where o + t d meets triangle (a, b, c) (Moller-Trumbore), HUGE_VAL if it does not */
static double bl_ray(double const* o, double const* d, double const* a, double const* b, double const* c)
{
	double e1[3], e2[3], p[3], s[3], q[3], det, u, v, t;
	bl_sub(b, a, e1); bl_sub(c, a, e2);
	bl_cross(d, e2, p);
	det = bl_dot(e1, p);
	if (fabs(det) < 1e-300)
		return HUGE_VAL;
	bl_sub(o, a, s);
	u = bl_dot(s, p) / det;
	if (u < 0 || u > 1)
		return HUGE_VAL;
	bl_cross(s, e1, q);
	v = bl_dot(d, q) / det;
	if (v < 0 || u + v > 1)
		return HUGE_VAL;
	t = bl_dot(e2, q) / det;
	return t > 0 ? t : HUGE_VAL;
}

/* does the segment o + t d, t in [0, tmax], reach the box lo..hi (slabs) */
static int bl_box_hit(double const* o, double const* d, double tmax, double const* lo, double const* hi)
{
	double t0 = 0, t1 = tmax;
	int k;
	for (k = 0; k < 3; k++)
	{
		if (fabs(d[k]) < 1e-300)
		{
			if (o[k] < lo[k] || o[k] > hi[k])
				return 0;
			continue;
		}
		{
			double a = (lo[k] - o[k]) / d[k], b = (hi[k] - o[k]) / d[k];
			if (a > b) { double t = a; a = b; b = t; }
			t0 = a > t0 ? a : t0;
			t1 = b < t1 ? b : t1;
			if (t0 > t1)
				return 0;
		}
	}
	return 1;
}

/* the surface's attributes on the stacks: a vertex's value on each of its copies, a corner's on the cells' corners over it,
   a face's on its stack's cells; edge ones not */
static qaws_status bl_attributes(qaws_volume_mesh* out, qaws_hm const* hm, unsigned int const* fo, unsigned int nl, int outward)
{
	unsigned int a, nv = hm->nv, nf = hm->nf, k, f, j, i;
	for (a = 0; a < hm->na; a++)
	{
		hm_attr const* t = &hm->attr[a];
		qaws_attribute_desc d;
		unsigned int id, n = t->n, at = 0;
		double* o;
		qaws_status st;
		if (t->element == QAWS_MESH_EDGE)
			continue;
		d.name = t->name;
		d.element = t->element == QAWS_MESH_FACE ? QAWS_MESH_CELL : (qaws_mesh_element)t->element;
		d.components = n;
		d.rule = (qaws_attribute_rule)t->rule;
		if ((st = qaws_volume_mesh_add_attribute(&id, out, NULL, &d)) != QAWS_STATUS_OK)
			return st;
		o = out->attr[id].v;
		if (t->element == QAWS_MESH_VERTEX)
			for (k = 0; k <= nl; k++)
				for (i = 0; i < nv; i++)
					memcpy(o + (size_t)n * ((size_t)k * nv + i), hm_attr_at(hm, a, i), sizeof(double) * n);
		else
			for (k = 0; k < nl; k++)
				for (f = 0; f < nf; f++)
				{
					if (t->element == QAWS_MESH_FACE)
					{
						memcpy(o + (size_t)n * ((size_t)k * nf + f), hm_attr_at(hm, a, f), sizeof(double) * n);
						continue;
					}
					for (j = 0; j < fo[f]; j++)
					{
						double const* v = hm_attr_at(hm, a, hm->fc[f] + (outward ? j : (fo[f] - j) % fo[f]));
						memcpy(o + (size_t)n * (at + j), v, sizeof(double) * n);
						memcpy(o + (size_t)n * (at + fo[f] + j), v, sizeof(double) * n);
					}
					at += 2 * fo[f];
				}
	}
	return QAWS_STATUS_OK;
}

qaws_status qaws_volume_mesh_grow_layers(qaws_volume_mesh** out_mesh, qaws_mesh_report* out_report, qaws_mesh_3d const* surface,
	qaws_boundary_layers_desc const* desc)
{
	qaws_hm const* hm;
	qaws_mesh_report rep;
	unsigned int nv, nf, nl, f, i, k, c, round, *fv = NULL, *fo = NULL, *cells = NULL, nt = 0, *tri = NULL, *tface = NULL;
	unsigned char* types = NULL;
	double *x = NULL, *fn = NULL, *dir = NULL, *sc = NULL, *lim = NULL, *y = NULL, *H = NULL, *tbox = NULL, growth, total, sgn;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (!out_mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_mesh = NULL;
	if (!surface || !desc || desc->layer_count < 1 || !(desc->first_height > 0) || (desc->growth != 0 && !(desc->growth >= 1)))
		return QAWS_STATUS_INVALID_ARGUMENT;
	hm = &surface->hm;
	nv = hm->nv;
	nf = hm->nf;
	nl = desc->layer_count;
	if (!nf || qaws_internal_mul2((size_t)nl + 1, nv) > UINT_MAX || qaws_internal_mul2(nl, nf) > UINT_MAX)
		return QAWS_STATUS_INVALID_ARGUMENT;   /* (no faces; more vertices or cells than an unsigned int counts) */
	for (f = 0; f < nf; f++)
		if (hm->fn[f] != 3 && hm->fn[f] != 4)
			return QAWS_STATUS_INVALID_ARGUMENT;
	memset(&rep, 0, sizeof rep);
	growth = desc->growth != 0 ? (double)desc->growth : 1.2;
	sgn = desc->outward ? 1.0 : -1.0;
	x = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)nv + 1));
	fv = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 4, (size_t)nf + 1));
	fo = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)nf + 1));
	fn = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)nf + 1));
	dir = (double*)qaws_internal_scratch_zalloc(3 * (size_t)nv + 3, sizeof(double));
	sc = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), (size_t)nv + 1));
	lim = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), (size_t)nv + 1));
	H = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), (size_t)nl + 2));
	y = (double*)qaws_internal_scratch_malloc(qaws_internal_mul4(sizeof(double), 3, (size_t)nv, (size_t)nl + 1));
	tri = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 6, (size_t)nf + 1));
	tface = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 2, (size_t)nf + 1));
	tbox = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 12, (size_t)nf + 1));
	if (!x || !fv || !fo || !fn || !dir || !sc || !lim || !H || !y || !tri || !tface || !tbox)
		goto done;
	memcpy(x, hm->x, sizeof(double) * 3 * nv);
	/* the faces' vertices (4 slots each), their unit normals (a quad's by its diagonals), the heights of the layers' tops */
	for (f = 0; f < nf; f++)
	{
		unsigned int n = hm->fn[f];
		double a[3], b[3];
		fo[f] = n;
		for (k = 0; k < n; k++)
			fv[4 * f + k] = hm->cv[hm->fc[f] + k];
		if (n == 3)
		{
			bl_sub(x + 3 * (size_t)fv[4 * f + 1], x + 3 * (size_t)fv[4 * f], a);
			bl_sub(x + 3 * (size_t)fv[4 * f + 2], x + 3 * (size_t)fv[4 * f], b);
		}
		else
		{
			bl_sub(x + 3 * (size_t)fv[4 * f + 2], x + 3 * (size_t)fv[4 * f], a);
			bl_sub(x + 3 * (size_t)fv[4 * f + 3], x + 3 * (size_t)fv[4 * f + 1], b);
		}
		bl_cross(a, b, fn + 3 * (size_t)f);
		if (!(bl_norm(fn + 3 * (size_t)f) > 0))
		{
			st = QAWS_STATUS_INVALID_ARGUMENT;   /* (a degenerate face) */
			goto done;
		}
		for (k = 0; k < 3; k++)
			fn[3 * (size_t)f + k] *= sgn;          /* (the way the layers grow) */
	}
	H[0] = 0;
	for (k = 1, total = 0; k <= nl; k++)
	{
		total += (double)desc->first_height * pow(growth, (double)(k - 1));
		H[k] = total;
	}
	/* each vertex's direction: its faces' normals weighted by their angles at it */
	for (f = 0; f < nf; f++)
		for (k = 0; k < fo[f]; k++)
		{
			unsigned int v = fv[4 * f + k], vn = fv[4 * f + (k + 1) % fo[f]], vp = fv[4 * f + (k + fo[f] - 1) % fo[f]];
			double a[3], b[3], ang;
			bl_sub(x + 3 * (size_t)vn, x + 3 * (size_t)v, a);
			bl_sub(x + 3 * (size_t)vp, x + 3 * (size_t)v, b);
			bl_norm(a);
			bl_norm(b);
			ang = acos(fmax(-1.0, fmin(1.0, bl_dot(a, b))));
			for (i = 0; i < 3; i++)
				dir[3 * (size_t)v + i] += ang * fn[3 * (size_t)f + i];
		}
	for (i = 0; i < nv; i++)
		bl_norm(dir + 3 * (size_t)i);
	/* turned until every face round the vertex sees it: pulled toward the faces that see it least, then each vertex's
	   least cosine (which also scales its height, at most twice) */
	for (round = 0; round < 50; round++)
	{
		unsigned int moved = 0;
		double* pull = lim;   /* (lim as scratch: 1 when the vertex moved) */
		for (i = 0; i < nv; i++)
			pull[i] = 0;
		for (f = 0; f < nf; f++)
			for (k = 0; k < fo[f]; k++)
			{
				unsigned int v = fv[4 * f + k];
				double cs = bl_dot(dir + 3 * (size_t)v, fn + 3 * (size_t)f);
				if (cs < BL_SEE + 0.05)
				{
					for (i = 0; i < 3; i++)
						dir[3 * (size_t)v + i] += 0.5 * (BL_SEE + 0.05 - cs) * fn[3 * (size_t)f + i];
					pull[v] = 1;
				}
			}
		for (i = 0; i < nv; i++)
			if (pull[i] != 0)
			{
				bl_norm(dir + 3 * (size_t)i);
				moved++;
			}
		if (!moved)
			break;
	}
	for (i = 0; i < nv; i++)
		sc[i] = 1;
	for (f = 0; f < nf; f++)
		for (k = 0; k < fo[f]; k++)
		{
			unsigned int v = fv[4 * f + k];
			double cs = bl_dot(dir + 3 * (size_t)v, fn + 3 * (size_t)f);
			if (!(cs > 0))
			{
				st = QAWS_STATUS_INVALID_ARGUMENT;   /* (no direction seen by every face round v) */
				goto done;
			}
			sc[v] = fmax(sc[v], 1 / fmax(cs, 0.5));   /* (the stretch, at most 2) */
		}
	for (i = 0; i < nv; i++)
		for (k = 0; k < 3; k++)
			dir[3 * (size_t)i + k] *= sc[i];          /* (the direction now carries the stretch) */
	/* the limits of each stack: the surface ahead along the direction (triangles not round the vertex, their boxes
	   first), and neighbours' directions converging */
	for (f = 0; f < nf; f++)
		for (k = 0; k + 2 < fo[f]; k++)   /* (a fan) */
		{
			unsigned int a = fv[4 * f], b = fv[4 * f + k + 1], cc = fv[4 * f + k + 2], j;
			double* bx = tbox + 6 * (size_t)nt;
			tri[3 * nt] = a; tri[3 * nt + 1] = b; tri[3 * nt + 2] = cc;
			tface[nt] = f;
			for (j = 0; j < 3; j++)
			{
				bx[j] = fmin(x[3 * (size_t)a + j], fmin(x[3 * (size_t)b + j], x[3 * (size_t)cc + j]));
				bx[3 + j] = fmax(x[3 * (size_t)a + j], fmax(x[3 * (size_t)b + j], x[3 * (size_t)cc + j]));
			}
			nt++;
		}
	for (i = 0; i < nv; i++)
	{
		double d[3], best = HUGE_VAL, len;
		unsigned int t;
		memcpy(d, dir + 3 * (size_t)i, sizeof d);
		len = bl_norm(d);
		if (hm->vc[i] == HM_NONE)
		{
			lim[i] = HUGE_VAL;
			continue;
		}
		for (t = 0; t < nt; t++)
		{
			double h;
			if (tri[3 * t] == i || tri[3 * t + 1] == i || tri[3 * t + 2] == i)
				continue;
			if (best < HUGE_VAL ? !bl_box_hit(x + 3 * (size_t)i, d, best, tbox + 6 * (size_t)t, tbox + 6 * (size_t)t + 3)
				: !bl_box_hit(x + 3 * (size_t)i, d, 1e300, tbox + 6 * (size_t)t, tbox + 6 * (size_t)t + 3))
				continue;
			h = bl_ray(x + 3 * (size_t)i, d, x + 3 * (size_t)tri[3 * t], x + 3 * (size_t)tri[3 * t + 1], x + 3 * (size_t)tri[3 * t + 2]);
			best = h < best ? h : best;
		}
		/* (the stack moves len * H[nl] along d) */
		lim[i] = best < HUGE_VAL ? BL_REACH * best / (len * total) : HUGE_VAL;
	}
	for (f = 0; f < nf; f++)
		for (k = 0; k < fo[f]; k++)
		{
			unsigned int a = fv[4 * f + k], b = fv[4 * f + (k + 1) % fo[f]];
			double e[3], dd[3], cv;
			bl_sub(x + 3 * (size_t)b, x + 3 * (size_t)a, e);
			bl_sub(dir + 3 * (size_t)b, dir + 3 * (size_t)a, dd);
			cv = -bl_dot(dd, e);   /* (> 0: converging; the points meet near s = |e|^2 / cv along the unit heights) */
			if (cv > 0)
			{
				double s = 0.5 * bl_dot(e, e) / cv / total;
				lim[a] = fmin(lim[a], s);
				lim[b] = fmin(lim[b], s);
			}
		}
	for (i = 0; i < nv; i++)
		sc[i] = fmin(1.0, lim[i]);
	/* graded, built, and halved where a cell turned over */
	cells = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul4(sizeof(unsigned int), 8, (size_t)nf, (size_t)nl));
	types = (unsigned char*)qaws_internal_scratch_malloc(qaws_internal_add2(qaws_internal_mul2((size_t)nf, (size_t)nl), 1));
	if (!cells || !types)
		goto done;
	for (round = 0; ; round++)
	{
		unsigned int bad = 0, sweep, at = 0;
		for (sweep = 0; sweep < 100; sweep++)
		{
			unsigned int changed = 0;
			for (f = 0; f < nf; f++)
				for (k = 0; k < fo[f]; k++)
				{
					unsigned int a = fv[4 * f + k], b = fv[4 * f + (k + 1) % fo[f]];
					if (sc[a] > BL_GRADE * sc[b]) { sc[a] = BL_GRADE * sc[b]; changed++; }
					if (sc[b] > BL_GRADE * sc[a]) { sc[b] = BL_GRADE * sc[a]; changed++; }
				}
			if (!changed)
				break;
		}
		for (k = 0; k <= nl; k++)
			for (i = 0; i < nv; i++)
				for (c = 0; c < 3; c++)
					y[3 * ((size_t)k * nv + i) + c] = x[3 * (size_t)i + c] + sc[i] * H[k] * dir[3 * (size_t)i + c];
		for (k = 0; k < nl; k++)
			for (f = 0; f < nf; f++)
			{
				unsigned int n = fo[f], j, *cv = cells + at;
				/* the base turning toward the top: the face's order when the layers grow along its normal, else reversed */
				for (j = 0; j < n; j++)
				{
					unsigned int v = fv[4 * f + (desc->outward ? j : (n - j) % n)];
					cv[j] = k * nv + v;
					cv[n + j] = (k + 1) * nv + v;
				}
				if (!bl_cell_ok(y, cv, n))
				{
					bad++;
					for (j = 0; j < n; j++)
						sc[cv[j] % nv] *= 0.5;
				}
				types[(size_t)k * nf + f] = (unsigned char)(n == 3 ? QAWS_CELL_PRISM : QAWS_CELL_HEXAHEDRON);
				at += 2 * n;
			}
		if (!bad)
			break;
		if (round + 1 >= BL_ROUNDS)
		{
			st = QAWS_STATUS_INTERNAL_ERROR;
			goto done;
		}
		rep.iterations++;
	}
	for (i = 0; i < nv; i++)
		rep.moves += hm->vc[i] != HM_NONE && sc[i] < 1;
	st = vm_create(y, (nl + 1) * nv, cells, types, nl * nf, out_mesh);
	if (st == QAWS_STATUS_OK)
		st = bl_attributes(*out_mesh, hm, fo, nl, desc->outward);
	if (st == QAWS_STATUS_OK)
	{
		qaws_mesh_3d* b = NULL;
		if (qaws_volume_mesh_get_boundary(&b, *out_mesh) == QAWS_STATUS_OK)
			qaws_mesh_3d_measure_quality(&rep.quality, b);
		qaws_mesh_3d_destroy(b);
	}
done:
	qaws_internal_free(x); qaws_internal_free(fv); qaws_internal_free(fo); qaws_internal_free(fn); qaws_internal_free(dir);
	qaws_internal_free(sc); qaws_internal_free(lim); qaws_internal_free(H); qaws_internal_free(y); qaws_internal_free(tri);
	qaws_internal_free(tface); qaws_internal_free(tbox); qaws_internal_free(cells); qaws_internal_free(types);
	if (st != QAWS_STATUS_OK)
	{
		qaws_volume_mesh_destroy(*out_mesh);
		*out_mesh = NULL;
		return st;
	}
	if (out_report)
		*out_report = rep;
	return QAWS_STATUS_OK;
}
