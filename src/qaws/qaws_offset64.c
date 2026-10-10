/*
 * Offsetting int64 paths: Clipper2's ClipperOffset (Angus Johnson, Boost
 * licence, https://www.angusj.com) ported to C on top of qaws_clip64.
 */

#include "qaws_offset.h"
#include "qaws_platform.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define O64_PI 3.141592653589793238
#define O64_FP_TOL 1e-12
#define O64_ARC_CONST 0.002

typedef struct o64_pt { double x, y; } o64_pt;

typedef struct o64_path
{
	int64_t* p;                 /* x, y pairs */
	unsigned int n, cap;
} o64_path;

typedef struct o64_ctx
{
	o64_path* sol;              /* the raw offset paths */
	unsigned int nsol, capsol;
	o64_path out;               /* the path being built */
	o64_pt* norms;
	unsigned int capnorms;
	double group_delta, step_sin, step_cos, steps_per_rad, temp_lim, arc_tolerance;
	qaws_join_type join_type;
	qaws_end_type end_type;
	int failed;
} o64_ctx;

/* Clipper2's Point64 from doubles: std::round (halves away from zero) */
static int64_t o64_round(double v)
{
	return (int64_t)(v < 0 ? -floor(-v + 0.5) : floor(v + 0.5));
}

static void o64_push(o64_ctx* c, int64_t x, int64_t y)
{
	o64_path* p = &c->out;
	if (p->n == p->cap)
	{
		unsigned int nc = p->cap ? 2 * p->cap : 64;
		int64_t* g = (int64_t*)realloc(p->p, sizeof(int64_t) * 2 * nc);
		if (!g) { c->failed = 1; return; }
		p->p = g;
		p->cap = nc;
	}
	p->p[2 * p->n] = x;
	p->p[2 * p->n + 1] = y;
	p->n++;
}

static void o64_pushd(o64_ctx* c, double x, double y)
{
	o64_push(c, o64_round(x), o64_round(y));
}

/* the path built so far into the solution */
static void o64_emit(o64_ctx* c)
{
	if (c->nsol == c->capsol)
	{
		unsigned int nc = c->capsol ? 2 * c->capsol : 16;
		o64_path* g = (o64_path*)realloc(c->sol, sizeof(o64_path) * nc);
		if (!g) { c->failed = 1; return; }
		c->sol = g;
		c->capsol = nc;
	}
	c->sol[c->nsol++] = c->out;
	memset(&c->out, 0, sizeof(c->out));
}

static o64_pt o64_unit_normal(int64_t const* a, int64_t const* b)
{
	o64_pt r = { 0, 0 };
	double dx, dy, ih;
	if (a[0] == b[0] && a[1] == b[1]) return r;
	dx = (double)(b[0] - a[0]);
	dy = (double)(b[1] - a[1]);
	ih = 1.0 / sqrt(dx * dx + dy * dy);
	dx *= ih;
	dy *= ih;
	r.x = dy;
	r.y = -dx;
	return r;
}

static o64_pt o64_normalize(o64_pt v)
{
	double h = sqrt(v.x * v.x + v.y * v.y);
	o64_pt r = { 0, 0 };
	if (fabs(h) < 0.001) return r;
	r.x = v.x / h;
	r.y = v.y / h;
	return r;
}

/* Clipper2's GetLineIntersectPt for PointD: ip clamped to the first segment */
static int o64_line_intersect(o64_pt a1, o64_pt a2, o64_pt b1, o64_pt b2, o64_pt* ip)
{
	double dx1 = a2.x - a1.x, dy1 = a2.y - a1.y, dx2 = b2.x - b1.x, dy2 = b2.y - b1.y;
	double det = dy1 * dx2 - dy2 * dx1, t;
	if (det == 0.0) return 0;
	t = ((a1.x - b1.x) * dy2 - (a1.y - b1.y) * dx2) / det;
	if (t <= 0.0) *ip = a1;
	else if (t >= 1.0) *ip = a2;
	else { ip->x = a1.x + t * dx1; ip->y = a1.y + t * dy1; }
	return 1;
}

static void o64_build_normals(o64_ctx* c, int64_t const* p, unsigned int n)
{
	unsigned int i;
	if (n > c->capnorms)
	{
		o64_pt* g = (o64_pt*)realloc(c->norms, sizeof(o64_pt) * n);
		if (!g) { c->failed = 1; return; }
		c->norms = g;
		c->capnorms = n;
	}
	for (i = 0; i + 1 < n; i++) c->norms[i] = o64_unit_normal(p + 2 * i, p + 2 * (i + 1));
	if (n) c->norms[n - 1] = o64_unit_normal(p + 2 * (n - 1), p);
}

static void o64_bevel(o64_ctx* c, int64_t const* p, unsigned int j, unsigned int k)
{
	o64_pt const* N = c->norms;
	double px = (double)p[2 * j], py = (double)p[2 * j + 1];
	if (j == k)
	{
		double ad = fabs(c->group_delta);
		o64_pushd(c, px - ad * N[j].x, py - ad * N[j].y);
		o64_pushd(c, px + ad * N[j].x, py + ad * N[j].y);
	}
	else
	{
		o64_pushd(c, px + c->group_delta * N[k].x, py + c->group_delta * N[k].y);
		o64_pushd(c, px + c->group_delta * N[j].x, py + c->group_delta * N[j].y);
	}
}

static void o64_square(o64_ctx* c, int64_t const* p, unsigned int j, unsigned int k)
{
	o64_pt const* N = c->norms;
	o64_pt vec, ptQ, pt1, pt2, pt3, pt4, pt;
	double ad = fabs(c->group_delta), gd = c->group_delta;
	if (j == k)
	{
		vec.x = N[j].y; vec.y = -N[j].x;
	}
	else
	{
		o64_pt s;
		s.x = -N[k].y + N[j].y;
		s.y = N[k].x + -N[j].x;
		vec = o64_normalize(s);
	}
	ptQ.x = (double)p[2 * j] + ad * vec.x;
	ptQ.y = (double)p[2 * j + 1] + ad * vec.y;
	pt1.x = ptQ.x + gd * vec.y; pt1.y = ptQ.y + gd * -vec.x;
	pt2.x = ptQ.x + gd * -vec.y; pt2.y = ptQ.y + gd * vec.x;
	pt3.x = (double)p[2 * k] + N[k].x * gd; pt3.y = (double)p[2 * k + 1] + N[k].y * gd;
	if (j == k)
	{
		pt4.x = pt3.x + vec.x * gd; pt4.y = pt3.y + vec.y * gd;
		pt = ptQ;
		o64_line_intersect(pt1, pt2, pt3, pt4, &pt);
		o64_pushd(c, ptQ.x + (ptQ.x - pt.x), ptQ.y + (ptQ.y - pt.y));
		o64_pushd(c, pt.x, pt.y);
	}
	else
	{
		pt4.x = (double)p[2 * j] + N[k].x * gd; pt4.y = (double)p[2 * j + 1] + N[k].y * gd;
		pt = ptQ;
		o64_line_intersect(pt1, pt2, pt3, pt4, &pt);
		o64_pushd(c, pt.x, pt.y);
		o64_pushd(c, ptQ.x + (ptQ.x - pt.x), ptQ.y + (ptQ.y - pt.y));
	}
}

static void o64_miter(o64_ctx* c, int64_t const* p, unsigned int j, unsigned int k, double cos_a)
{
	o64_pt const* N = c->norms;
	double q = c->group_delta / (cos_a + 1);
	o64_pushd(c, (double)p[2 * j] + (N[k].x + N[j].x) * q, (double)p[2 * j + 1] + (N[k].y + N[j].y) * q);
}

static void o64_round_join(o64_ctx* c, int64_t const* p, unsigned int j, unsigned int k, double angle)
{
	o64_pt const* N = c->norms;
	double px = (double)p[2 * j], py = (double)p[2 * j + 1];
	o64_pt v;
	int steps, i;
	v.x = N[k].x * c->group_delta;
	v.y = N[k].y * c->group_delta;
	if (j == k) { v.x = -v.x; v.y = -v.y; }
	o64_pushd(c, px + v.x, py + v.y);
	steps = (int)ceil(c->steps_per_rad * fabs(angle));
	for (i = 1; i < steps; i++)
	{
		o64_pt w;
		w.x = v.x * c->step_cos - c->step_sin * v.y;
		w.y = v.x * c->step_sin + v.y * c->step_cos;
		v = w;
		o64_pushd(c, px + v.x, py + v.y);
	}
	o64_pushd(c, px + N[j].x * c->group_delta, py + N[j].y * c->group_delta);
}

static void o64_offset_point(o64_ctx* c, int64_t const* p, unsigned int j, unsigned int k)
{
	o64_pt const* N = c->norms;
	double sin_a, cos_a;
	if (p[2 * j] == p[2 * k] && p[2 * j + 1] == p[2 * k + 1]) return;
	/* Clipper2's CrossProduct(vec1, vec2) = vec1.y vec2.x - vec2.y vec1.x */
	sin_a = N[j].y * N[k].x - N[k].y * N[j].x;
	cos_a = N[j].x * N[k].x + N[j].y * N[k].y;
	if (sin_a > 1.0) sin_a = 1.0;
	else if (sin_a < -1.0) sin_a = -1.0;
	if (fabs(c->group_delta) <= O64_FP_TOL)
	{
		o64_push(c, p[2 * j], p[2 * j + 1]);
		return;
	}
	if (cos_a > -0.999 && sin_a * c->group_delta < 0)
	{
		/* concave: through the vertex (the union removes the loops) */
		double px = (double)p[2 * j], py = (double)p[2 * j + 1];
		o64_pushd(c, px + N[k].x * c->group_delta, py + N[k].y * c->group_delta);
		o64_push(c, p[2 * j], p[2 * j + 1]);
		o64_pushd(c, px + N[j].x * c->group_delta, py + N[j].y * c->group_delta);
	}
	else if (cos_a > 0.999 && c->join_type != QAWS_JOIN_ROUND)
		o64_miter(c, p, j, k, cos_a);
	else if (c->join_type == QAWS_JOIN_MITER)
	{
		if (cos_a > c->temp_lim - 1) o64_miter(c, p, j, k, cos_a);
		else o64_square(c, p, j, k);
	}
	else if (c->join_type == QAWS_JOIN_ROUND)
		o64_round_join(c, p, j, k, atan2(sin_a, cos_a));
	else if (c->join_type == QAWS_JOIN_BEVEL)
		o64_bevel(c, p, j, k);
	else
		o64_square(c, p, j, k);
}

static void o64_offset_polygon(o64_ctx* c, int64_t const* p, unsigned int n)
{
	unsigned int j, k;
	c->out.n = 0;
	for (j = 0, k = n - 1; j < n; k = j, j++)
		o64_offset_point(c, p, j, k);
	o64_emit(c);
}

static void o64_offset_open_joined(o64_ctx* c, int64_t const* p, unsigned int n)
{
	int64_t* r = (int64_t*)malloc(sizeof(int64_t) * 2 * n);
	unsigned int i;
	o64_offset_polygon(c, p, n);
	if (!r) { c->failed = 1; return; }
	for (i = 0; i < n; i++) { r[2 * i] = p[2 * (n - 1 - i)]; r[2 * i + 1] = p[2 * (n - 1 - i) + 1]; }
	/* the normals reversed, rotated by one, negated */
	{
		o64_pt* tmp = (o64_pt*)malloc(sizeof(o64_pt) * n);
		if (!tmp) { free(r); c->failed = 1; return; }
		for (i = 0; i < n; i++) tmp[i] = c->norms[n - 1 - i];
		for (i = 0; i < n; i++) { c->norms[i].x = -tmp[(i + 1) % n].x; c->norms[i].y = -tmp[(i + 1) % n].y; }
		free(tmp);
	}
	o64_offset_polygon(c, r, n);
	free(r);
}

static void o64_cap(o64_ctx* c, int64_t const* p, unsigned int i)
{
	if (fabs(c->group_delta) <= O64_FP_TOL)
		o64_push(c, p[2 * i], p[2 * i + 1]);
	else if (c->end_type == QAWS_END_BUTT)
		o64_bevel(c, p, i, i);
	else if (c->end_type == QAWS_END_ROUND)
		o64_round_join(c, p, i, i, O64_PI);
	else
		o64_square(c, p, i, i);
}

static void o64_offset_open(o64_ctx* c, int64_t const* p, unsigned int n)
{
	unsigned int hi = n - 1, i, j, k;
	c->out.n = 0;
	o64_cap(c, p, 0);
	for (j = 1, k = 0; j < hi; k = j, j++)
		o64_offset_point(c, p, j, k);
	for (i = hi; i > 0; i--) { c->norms[i].x = -c->norms[i - 1].x; c->norms[i].y = -c->norms[i - 1].y; }
	c->norms[0] = c->norms[hi];
	o64_cap(c, p, hi);
	for (j = hi - 1, k = hi; j > 0; k = j, j--)
		o64_offset_point(c, p, j, k);
	o64_emit(c);
}

/* a path without consecutive repeats (and for closed ones without a last
   point equal to the first) */
static unsigned int o64_strip(int64_t const* in, unsigned int n, int closed, int64_t* out)
{
	unsigned int i, m = 0;
	for (i = 0; i < n; i++)
	{
		if (m && out[2 * m - 2] == in[2 * i] && out[2 * m - 1] == in[2 * i + 1]) continue;
		out[2 * m] = in[2 * i];
		out[2 * m + 1] = in[2 * i + 1];
		m++;
	}
	if (closed)
		while (m > 1 && out[2 * m - 2] == out[0] && out[2 * m - 1] == out[1]) m--;
	return m;
}

qaws_status qaws_offset64_execute(qaws_offset64_desc const* d, qaws_clip64_result** out)
{
	o64_ctx c;
	unsigned int g, i, total = 0;
	int first_reversed = 0, have_polygon = 0;
	double delta;
	qaws_status st = QAWS_STATUS_OK;
	int64_t** stripped = NULL;
	unsigned int* slen = NULL;
	if (!d || !out || (d->group_count && !d->groups)) return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	memset(&c, 0, sizeof(c));
	delta = d->delta;
	c.temp_lim = (d->miter_limit > 0 ? d->miter_limit : 2.0) <= 1 ? 2.0 : 2.0 / ((d->miter_limit > 0 ? d->miter_limit : 2.0) * (d->miter_limit > 0 ? d->miter_limit : 2.0));
	c.arc_tolerance = d->arc_tolerance;
	for (g = 0; g < d->group_count; g++) total += d->groups[g].path_count;
	stripped = (int64_t**)calloc(total + 1, sizeof(int64_t*));
	slen = (unsigned int*)calloc(total + 1, sizeof(unsigned int));
	if (!stripped || !slen) { free(stripped); free(slen); return QAWS_STATUS_ALLOCATION_FAILURE; }
	total = 0;
	for (g = 0; g < d->group_count && !c.failed; g++)
	{
		qaws_offset64_group const* gr = &d->groups[g];
		int joined = gr->end_type == QAWS_END_POLYGON || gr->end_type == QAWS_END_JOINED, reversed = 0, lowest = -1;
		unsigned int base = total;
		/* strip, and find the lowest closed path (its orientation decides) */
		for (i = 0; i < gr->path_count; i++, total++)
		{
			stripped[total] = (int64_t*)malloc(sizeof(int64_t) * 2 * (gr->paths[i].point_count + 1));
			if (!stripped[total]) { c.failed = 1; break; }
			slen[total] = o64_strip(gr->paths[i].points, gr->paths[i].point_count, joined, stripped[total]);
		}
		if (c.failed) break;
		if (gr->end_type == QAWS_END_POLYGON)
		{
			int64_t bx = INT64_MAX, by = INT64_MIN;
			int neg = 0;
			for (i = 0; i < gr->path_count; i++)
			{
				int64_t const* p = stripped[base + i];
				unsigned int k;
				double a = 1e300;
				for (k = 0; k < slen[base + i]; k++)
				{
					if (p[2 * k + 1] < by || (p[2 * k + 1] == by && p[2 * k] >= bx)) continue;
					if (a == 1e300)
					{
						a = qaws_path64_area2(p, slen[base + i]);
						if (a == 0) break;
						neg = a < 0;
					}
					lowest = (int)i;
					bx = p[2 * k];
					by = p[2 * k + 1];
				}
			}
			reversed = lowest >= 0 && neg;
			if (!have_polygon) { first_reversed = reversed; have_polygon = 1; }
		}
		if (fabs(delta) < 0.5)
		{
			/* an insignificant offset: the paths as they are */
			for (i = 0; i < gr->path_count && !c.failed; i++)
			{
				unsigned int k;
				c.out.n = 0;
				for (k = 0; k < slen[base + i]; k++) o64_push(&c, stripped[base + i][2 * k], stripped[base + i][2 * k + 1]);
				o64_emit(&c);
			}
			continue;
		}
		if (gr->end_type == QAWS_END_POLYGON)
		{
			double dd = lowest < 0 ? fabs(delta) : delta;
			c.group_delta = reversed ? -dd : dd;
		}
		else
			c.group_delta = fabs(delta);
		c.join_type = gr->join_type;
		c.end_type = gr->end_type;
		if (gr->join_type == QAWS_JOIN_ROUND || gr->end_type == QAWS_END_ROUND)
		{
			double ad = fabs(c.group_delta);
			double tol = c.arc_tolerance > O64_FP_TOL ? (c.arc_tolerance < ad ? c.arc_tolerance : ad) : ad * O64_ARC_CONST;
			double steps = O64_PI / acos(1 - tol / ad);
			if (steps > ad * O64_PI) steps = ad * O64_PI;
			c.step_sin = sin(2 * O64_PI / steps);
			c.step_cos = cos(2 * O64_PI / steps);
			if (c.group_delta < 0) c.step_sin = -c.step_sin;
			c.steps_per_rad = steps / (2 * O64_PI);
		}
		for (i = 0; i < gr->path_count && !c.failed; i++)
		{
			int64_t const* p = stripped[base + i];
			unsigned int n = slen[base + i];
			qaws_end_type et = gr->end_type;
			c.out.n = 0;
			if (n == 0) continue;
			if (n == 1)
			{
				double ad = fabs(c.group_delta);
				if (c.group_delta < 1) continue;
				if (gr->join_type == QAWS_JOIN_ROUND)
				{
					/* Clipper2's Ellipse */
					unsigned int steps = c.steps_per_rad > 0 ? (unsigned int)ceil(c.steps_per_rad * 2 * O64_PI) : 0, k;
					double si, co, dx, dy;
					if (steps <= 2) steps = (unsigned int)(O64_PI * sqrt(ad));
					si = sin(2 * O64_PI / steps); co = cos(2 * O64_PI / steps);
					dx = co; dy = si;
					o64_pushd(&c, (double)p[0] + ad, (double)p[1]);
					for (k = 1; k < steps; k++)
					{
						double x;
						o64_pushd(&c, (double)p[0] + ad * dx, (double)p[1] + ad * dy);
						x = dx * co - dy * si;
						dy = dy * co + dx * si;
						dx = x;
					}
				}
				else
				{
					int64_t e = (int64_t)ceil(ad);
					o64_push(&c, p[0] - e, p[1] - e);
					o64_push(&c, p[0] + e, p[1] - e);
					o64_push(&c, p[0] + e, p[1] + e);
					o64_push(&c, p[0] - e, p[1] + e);
				}
				o64_emit(&c);
				continue;
			}
			if (n == 2 && et == QAWS_END_JOINED)
				et = gr->join_type == QAWS_JOIN_ROUND ? QAWS_END_ROUND : QAWS_END_SQUARE;
			c.end_type = et;
			o64_build_normals(&c, p, n);
			if (c.failed) break;
			if (et == QAWS_END_POLYGON) o64_offset_polygon(&c, p, n);
			else if (et == QAWS_END_JOINED) o64_offset_open_joined(&c, p, n);
			else o64_offset_open(&c, p, n);
		}
	}
	if (!c.failed)
	{
		/* the clean-up: one exact union, positive (negative when reversed),
		   keeping the input's orientation */
		qaws_path64* paths = (qaws_path64*)malloc(sizeof(qaws_path64) * (c.nsol + 1));
		qaws_clip64_desc cd;
		if (!paths) st = QAWS_STATUS_ALLOCATION_FAILURE;
		else
		{
			for (i = 0; i < c.nsol; i++) { paths[i].points = c.sol[i].p; paths[i].point_count = c.sol[i].n; }
			memset(&cd, 0, sizeof(cd));
			cd.subjects = paths;
			cd.subject_count = c.nsol;
			cd.clip_type = QAWS_CLIP_UNION;
			cd.fill_rule = first_reversed ? QAWS_FILL_NEGATIVE : QAWS_FILL_POSITIVE;
			cd.flags = (d->flags & QAWS_CLIP_PRESERVE_COLLINEAR)
				| (((d->flags & QAWS_CLIP_REVERSE_SOLUTION) != 0) != (first_reversed != 0) ? QAWS_CLIP_REVERSE_SOLUTION : 0u);
			st = qaws_clip64_execute(&cd, out);
			free(paths);
		}
	}
	else
		st = QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < c.nsol; i++) free(c.sol[i].p);
	free(c.sol);
	free(c.out.p);
	free(c.norms);
	for (i = 0; i < total; i++) free(stripped[i]);
	free(stripped);
	free(slen);
	return st;
}

qaws_status qaws_offset64_paths(qaws_path64 const* paths, unsigned int path_count, double delta,
	qaws_join_type join_type, qaws_end_type end_type, qaws_clip64_result** out)
{
	qaws_offset64_group g;
	qaws_offset64_desc d;
	g.paths = paths;
	g.path_count = path_count;
	g.join_type = join_type;
	g.end_type = end_type;
	memset(&d, 0, sizeof(d));
	d.groups = &g;
	d.group_count = 1;
	d.delta = delta;
	return qaws_offset64_execute(&d, out);
}
