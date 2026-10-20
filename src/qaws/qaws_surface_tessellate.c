/*
 * Adaptive surface tessellation (qaws_inspect.h) on an implicit quadtree over the (u, v) domain. Whether a quad splits is a
 * function of the quad alone (its depth and its 3 x 3 test points), so the tree is never stored: both passes, and the walks
 * along neighbours, recompute the same decisions. A leaf's sides carry the vertices of finer neighbours (found by descending
 * the neighbour along the shared side), so the mesh has no T-junction cracks; a leaf with none is two triangles, one with
 * some a fan around its center. Quad corners sit on dyadic parameters, so a vertex shared by leaves is evaluated at the
 * same (u, v) and comes out the same bits in each.
 * shortcut: vertices are written per leaf (shared ones repeated, identical); weld by (u, v) where one index per point is
 * needed.
 * shortcut: no 2:1 grading, so a coarse leaf next to much finer ones fans out to all their vertices (thin triangles);
 * grade when the tree is stored.
 */

#include "qaws_inspect.h"
#include "qaws_surface.h"
#include "internal/qaws_internal_surface.h"
#include <math.h>
#include <string.h>

#define TS_MIN_DEPTH 2
#define TS_MAX_DEPTH_DEFAULT 10
#define TS_MAX_DEPTH 20

typedef struct ts_ctx
{
	qaws_surface const* s;
	double u0, du, v0, dv;
	double chord, angle;
	unsigned int min_depth, max_depth;
	/* output */
	qaws_tessellation_vertex* verts;
	unsigned int* indices;
	unsigned int vcap, icap, vn, in;
	/* the leaf being fanned */
	unsigned int first, prev, count;
	qaws_status status;
} ts_ctx;

static double ts_u(ts_ctx const* c, unsigned int k, unsigned int d)
{
	return c->u0 + c->du * ldexp((double)k, -(int)d);
}

static double ts_v(ts_ctx const* c, unsigned int k, unsigned int d)
{
	return c->v0 + c->dv * ldexp((double)k, -(int)d);
}

static void ts_eval(ts_ctx* c, double u, double v, qaws_surface_eval_result* r)
{
	qaws_status st;
	memset(r, 0, sizeof *r);
	st = qaws_surface_evaluate(r, c->s, (qaws_scalar)u, (qaws_scalar)v, QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL);
	if (st != QAWS_STATUS_OK && c->status == QAWS_STATUS_OK)
		c->status = st;
}

static double ts_angle(qaws_vec3 a, qaws_vec3 b)
{
	double x = (double)a.y * b.z - (double)a.z * b.y, y = (double)a.z * b.x - (double)a.x * b.z, z = (double)a.x * b.y - (double)a.y * b.x;
	double d = (double)a.x * b.x + (double)a.y * b.y + (double)a.z * b.z;
	if ((double)a.x * a.x + (double)a.y * a.y + (double)a.z * a.z == 0 || (double)b.x * b.x + (double)b.y * b.y + (double)b.z * b.z == 0)
		return 0;
	return atan2(sqrt(x * x + y * y + z * z), d);
}

/* the distance from p to the triangle abc (Ericson, Real-Time Collision Detection 5.1.5) */
static double ts_point_triangle(qaws_vec3 p, qaws_vec3 a, qaws_vec3 b, qaws_vec3 c)
{
	double ab[3] = { b.x - a.x, b.y - a.y, b.z - a.z }, ac[3] = { c.x - a.x, c.y - a.y, c.z - a.z };
	double ap[3] = { p.x - a.x, p.y - a.y, p.z - a.z }, bp[3] = { p.x - b.x, p.y - b.y, p.z - b.z }, cp[3] = { p.x - c.x, p.y - c.y, p.z - c.z };
	double d1 = ab[0] * ap[0] + ab[1] * ap[1] + ab[2] * ap[2], d2 = ac[0] * ap[0] + ac[1] * ap[1] + ac[2] * ap[2];
	double d3 = ab[0] * bp[0] + ab[1] * bp[1] + ab[2] * bp[2], d4 = ac[0] * bp[0] + ac[1] * bp[1] + ac[2] * bp[2];
	double d5 = ab[0] * cp[0] + ab[1] * cp[1] + ab[2] * cp[2], d6 = ac[0] * cp[0] + ac[1] * cp[1] + ac[2] * cp[2];
	double va = 0, vb = 0, vc = 0, v, w, q[3];
	unsigned int k;
	if (d1 <= 0 && d2 <= 0)
		v = 0, w = 0;
	else if (d3 >= 0 && d4 <= d3)
		v = 1, w = 0;
	else if (d6 >= 0 && d5 <= d6)
		v = 0, w = 1;
	else if ((vc = d1 * d4 - d3 * d2) <= 0 && d1 >= 0 && d3 <= 0)
		v = d1 / (d1 - d3), w = 0;
	else if ((vb = d5 * d2 - d1 * d6) <= 0 && d2 >= 0 && d6 <= 0)
		v = 0, w = d2 / (d2 - d6);
	else if ((va = d3 * d6 - d5 * d4) <= 0 && d4 - d3 >= 0 && d5 - d6 >= 0)
		w = (d4 - d3) / ((d4 - d3) + (d5 - d6)), v = 1 - w;
	else
	{
		double den = 1 / (va + vb + vc);
		v = vb * den;
		w = vc * den;
	}
	for (k = 0; k < 3; k++)
		q[k] = ap[k] - v * ab[k] - w * ac[k];
	return sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);
}

/* quad (d, i, j) splits: under the minimum depth; or below the maximum, a test point (side midpoints, center, the two
   triangles' centroids in (u, v)) farther than the chord tolerance from the quad's two triangles, or a normal turned from
   the center's by more than half the angle tolerance */
static int ts_split(ts_ctx* c, unsigned int d, unsigned int i, unsigned int j)
{
	qaws_surface_eval_result r[3][3], g[2];
	double u0, u1, v0, v1;
	unsigned int a, b;
	if (d < c->min_depth)
		return 1;
	if (d >= c->max_depth)
		return 0;
	for (a = 0; a < 3; a++)
		for (b = 0; b < 3; b++)
			ts_eval(c, ts_u(c, 2 * i + a, d + 1), ts_v(c, 2 * j + b, d + 1), &r[a][b]);
	if (c->angle > 0)
		for (a = 0; a < 3; a++)
			for (b = 0; b < 3; b++)
				if (ts_angle(r[1][1].normal, r[a][b].normal) > 0.5 * c->angle)
					return 1;
	if (c->chord > 0)
	{
		qaws_vec3 p00 = r[0][0].position, p10 = r[2][0].position, p11 = r[2][2].position, p01 = r[0][2].position;
		u0 = ts_u(c, i, d);
		u1 = ts_u(c, i + 1, d);
		v0 = ts_v(c, j, d);
		v1 = ts_v(c, j + 1, d);
		ts_eval(c, u0 + (u1 - u0) * (2.0 / 3), v0 + (v1 - v0) / 3, &g[0]);
		ts_eval(c, u0 + (u1 - u0) / 3, v0 + (v1 - v0) * (2.0 / 3), &g[1]);
		qaws_vec3 t[7];
		t[0] = r[1][0].position; t[1] = r[2][1].position; t[2] = r[1][2].position; t[3] = r[0][1].position; t[4] = r[1][1].position;
		t[5] = g[0].position; t[6] = g[1].position;
		for (a = 0; a < 7; a++)
		{
			double e0 = ts_point_triangle(t[a], p00, p10, p11), e1 = ts_point_triangle(t[a], p00, p11, p01);
			if ((e0 < e1 ? e0 : e1) > c->chord)
				return 1;
		}
	}
	return 0;
}

/* quad (d, i, j) is a node of the tree: every ancestor splits. Ancestors shared with the leaf (d, li, lj) split already. */
static int ts_exists(ts_ctx* c, unsigned int d, unsigned int i, unsigned int j, unsigned int li, unsigned int lj)
{
	unsigned int k;
	for (k = 1; k <= d; k++)
	{
		if ((i >> k) == (li >> k) && (j >> k) == (lj >> k))
			return 1;
		if (!ts_split(c, d - k, i >> k, j >> k))
			return 0;
	}
	return 1;
}

static void ts_vertex(ts_ctx* c, unsigned int ku, unsigned int kv, unsigned int d)
{
	unsigned int n = c->vn++;
	if (c->verts && n < c->vcap)
	{
		qaws_surface_eval_result r;
		double u = ts_u(c, ku, d), v = ts_v(c, kv, d);
		ts_eval(c, u, v, &r);
		c->verts[n].position = r.position;
		c->verts[n].normal = r.normal;
		c->verts[n].u = (qaws_scalar)u;
		c->verts[n].v = (qaws_scalar)v;
	}
}

static void ts_triangle(ts_ctx* c, unsigned int a, unsigned int b, unsigned int t)
{
	if (c->indices && c->in + 3 <= c->icap)
	{
		c->indices[c->in] = a;
		c->indices[c->in + 1] = b;
		c->indices[c->in + 2] = t;
	}
	c->in += 3;
}

/* the next perimeter vertex of the leaf being fanned (the fan's center is its first vertex) */
static void ts_rim(ts_ctx* c, unsigned int ku, unsigned int kv, unsigned int d)
{
	unsigned int n = c->vn;
	ts_vertex(c, ku, kv, d);
	if (c->count == 0)
		c->first = n;
	else
		ts_triangle(c, c->first - 1, c->prev, n);
	c->prev = n;
	c->count++;
}

/* the vertices inside side `side` (0: v = v0, 1: u = u1, 2: v = v1, 3: u = u0) of node (d, i, j), from its splits; forward
   runs with u or v increasing */
static void ts_side(ts_ctx* c, unsigned int d, unsigned int i, unsigned int j, unsigned int side, int forward)
{
	unsigned int ci[2], cj[2], k;
	if (!ts_split(c, d, i, j))
		return;
	for (k = 0; k < 2; k++)
	{
		ci[k] = side == 1 ? 2 * i + 1 : side == 3 ? 2 * i : 2 * i + k;
		cj[k] = side == 2 ? 2 * j + 1 : side == 0 ? 2 * j : 2 * j + k;
	}
	k = forward ? 0 : 1;
	ts_side(c, d + 1, ci[k], cj[k], side, forward);
	/* the side's midpoint */
	if (side == 0 || side == 2)
		ts_rim(c, 2 * i + 1, side == 0 ? 2 * j : 2 * j + 2, d + 1);
	else
		ts_rim(c, side == 3 ? 2 * i : 2 * i + 2, 2 * j + 1, d + 1);
	ts_side(c, d + 1, ci[1 - k], cj[1 - k], side, forward);
}

/* the neighbour of leaf (d, i, j) across its side `side` contributes its vertices on that side, when it is a node */
static int ts_has_neighbour(ts_ctx* c, unsigned int d, unsigned int i, unsigned int j, unsigned int side, unsigned int* ni, unsigned int* nj)
{
	unsigned int n = 1u << d;
	*ni = i;
	*nj = j;
	if (side == 0) { if (j == 0) return 0; *nj = j - 1; }
	if (side == 1) { if (i + 1 == n) return 0; *ni = i + 1; }
	if (side == 2) { if (j + 1 == n) return 0; *nj = j + 1; }
	if (side == 3) { if (i == 0) return 0; *ni = i - 1; }
	return ts_exists(c, d, *ni, *nj, i, j);
}

static void ts_leaf(ts_ctx* c, unsigned int d, unsigned int i, unsigned int j)
{
	unsigned int ni[4], nj[4], has[4], s, base;
	for (s = 0; s < 4; s++)
		has[s] = ts_has_neighbour(c, d, i, j, s, &ni[s], &nj[s]) && ts_split(c, d, ni[s], nj[s]);
	if (!has[0] && !has[1] && !has[2] && !has[3])
	{
		/* two triangles */
		base = c->vn;
		ts_vertex(c, i, j, d);
		ts_vertex(c, i + 1, j, d);
		ts_vertex(c, i + 1, j + 1, d);
		ts_vertex(c, i, j + 1, d);
		ts_triangle(c, base, base + 1, base + 2);
		ts_triangle(c, base, base + 2, base + 3);
		return;
	}
	/* a fan: the center, then the rim counter-clockwise; each neighbour side seen from the neighbour (its opposite side) */
	ts_vertex(c, 2 * i + 1, 2 * j + 1, d + 1);
	c->count = 0;
	ts_rim(c, i, j, d);
	if (has[0])
		ts_side(c, d, ni[0], nj[0], 2, 1);
	ts_rim(c, i + 1, j, d);
	if (has[1])
		ts_side(c, d, ni[1], nj[1], 3, 1);
	ts_rim(c, i + 1, j + 1, d);
	if (has[2])
		ts_side(c, d, ni[2], nj[2], 0, 0);
	ts_rim(c, i, j + 1, d);
	if (has[3])
		ts_side(c, d, ni[3], nj[3], 1, 0);
	ts_triangle(c, c->first - 1, c->prev, c->first);
}

static void ts_walk(ts_ctx* c, unsigned int d, unsigned int i, unsigned int j)
{
	if (c->status != QAWS_STATUS_OK)
		return;
	if (!ts_split(c, d, i, j))
	{
		ts_leaf(c, d, i, j);
		return;
	}
	ts_walk(c, d + 1, 2 * i, 2 * j);
	ts_walk(c, d + 1, 2 * i + 1, 2 * j);
	ts_walk(c, d + 1, 2 * i, 2 * j + 1);
	ts_walk(c, d + 1, 2 * i + 1, 2 * j + 1);
}

qaws_status qaws_surface_tessellate_adaptive(qaws_surface const* surface, qaws_surface_tessellation_desc const* desc,
	qaws_tessellation_vertex* out_vertices, unsigned int vertex_capacity, unsigned int* out_vertex_count,
	unsigned int* out_indices, unsigned int index_capacity, unsigned int* out_index_count)
{
	ts_ctx c;
	qaws_range ur, vr;
	if (!surface || !desc || !out_vertex_count || !out_index_count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (!(desc->chord_tolerance >= 0) || !(desc->angle_tolerance >= 0) || !(desc->chord_tolerance > 0 || desc->angle_tolerance > 0))
		return QAWS_STATUS_INVALID_ARGUMENT;
	memset(&c, 0, sizeof c);
	c.s = surface;
	ur = qaws_surface_get_u_range(surface);
	vr = qaws_surface_get_v_range(surface);
	c.u0 = ur.min_value;
	c.du = (double)ur.max_value - ur.min_value;
	c.v0 = vr.min_value;
	c.dv = (double)vr.max_value - vr.min_value;
	c.chord = desc->chord_tolerance;
	c.angle = desc->angle_tolerance;
	c.max_depth = desc->max_depth ? desc->max_depth : TS_MAX_DEPTH_DEFAULT;
	c.max_depth = c.max_depth > TS_MAX_DEPTH ? TS_MAX_DEPTH : c.max_depth;
	c.min_depth = desc->min_depth ? desc->min_depth : TS_MIN_DEPTH;
	c.min_depth = c.min_depth > c.max_depth ? c.max_depth : c.min_depth;
	c.verts = out_vertices;
	c.indices = out_indices;
	c.vcap = out_vertices ? vertex_capacity : 0;
	c.icap = out_indices ? index_capacity : 0;
	c.status = QAWS_STATUS_OK;
	ts_walk(&c, 0, 0, 0);
	*out_vertex_count = c.vn;
	*out_index_count = c.in;
	if (c.status != QAWS_STATUS_OK)
		return c.status;
	if ((out_vertices && c.vn > vertex_capacity) || (out_indices && c.in > index_capacity))
		return QAWS_STATUS_BUFFER_TOO_SMALL;
	return QAWS_STATUS_OK;
}
