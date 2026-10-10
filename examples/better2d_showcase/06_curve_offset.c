/* Figure 6: offsets of curves past their radii of curvature (qaws_offset,
   after "Fast GPU stroke expansion", Levien and Uguray, HPG 2024). The raw
   offset chain before the union: G1 cubics along the offset, its cusps found
   exactly where 1 + delta k = 0, the evolute (the centres of curvature)
   through each backward run; then the regions after the positive union, and
   the same offset fitted the old way (cubics through offset points) beside
   the new one. */

#include "internal/qaws_internal_offset.h"

#define CO_PI 3.14159265358979323846

/* a curve drawn span by span (chains have many short cubics) */
static void co_curve(svg* s, viewport const* v, qaws_curve const* c, char const* color, double width, double op, int dashed)
{
	unsigned int n = qaws_curve_get_span_count(c), m = 24 * (n ? n : 1), j;
	double* xy = (double*)malloc(sizeof(double) * 2 * (m + 1));
	qaws_range r = qaws_curve_get_parameter_range(c);
	for (j = 0; j <= m; j++)
	{
		qaws_eval_result_2d e;
		qaws_curve_evaluate_2d(c, (qaws_scalar)(r.min_value + (r.max_value - r.min_value) * j / m), QAWS_EVAL_FLAG_POSITION, &e);
		xy[2 * j] = vx(v, e.position.x);
		xy[2 * j + 1] = vy(v, e.position.y);
	}
	svg_polyline(s, xy, (int)m + 1, color, width, op, dashed);
	free(xy);
}

/* the joints of a cubic chain (every third control point) */
static unsigned int co_joints(svg* s, viewport const* v, qaws_curve const* c, char const* fill, double rad)
{
	unsigned int n = 3 * qaws_curve_get_span_count(c) + 1, got = 0, i;
	qaws_scalar* cp = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * (n + 4));
	if (qaws_curve_get_kind(c) == QAWS_CURVE_KIND_BSPLINE && qaws_curve_get_degree(c) == 3)
	{
		qaws_curve_get_control_points(c, cp, n + 4, &got);
		for (i = 0; i < got; i += 3)
			svg_circle(s, vx(v, cp[2 * i]), vy(v, cp[2 * i + 1]), rad, fill, "#ffffff");
	}
	free(cp);
	return got ? (got - 1) / 3 : 0;
}

static void co_result(svg* s, viewport const* v, qaws_clip_result const* r, char const* fill, char const* line, unsigned int* cubics, int joints)
{
	unsigned int i, k;
	if (fill) bo_fill(s, v, r, fill);
	for (i = 0; i < qaws_clip_result_get_path_count(r); i++)
	{
		qaws_path_2d p;
		qaws_clip_result_get_path(r, i, &p);
		for (k = 0; k < p.curve_count; k++)
		{
			co_curve(s, v, p.curves[k], line, 1.8, 0.95, 0);
			if (cubics) *cubics += qaws_curve_get_span_count(p.curves[k]);
			if (joints) co_joints(s, v, p.curves[k], line, 2.4);
		}
	}
}

/* the plain parallel curve C + delta n, swallowtails and all (dashed) and
   the evolute C - n / k where the offset runs backward (dotted purple) */
static void co_parallel(svg* s, viewport const* v, qaws_curve const* c, double delta)
{
	enum { N = 1600 };
	double xy[2 * (N + 1)], ev[2 * (N + 1)];
	qaws_range r = qaws_curve_get_parameter_range(c);
	unsigned int j, ne = 0;
	for (j = 0; j <= N; j++)
	{
		qaws_eval_result_2d e;
		double t = r.min_value + (r.max_value - r.min_value) * j / N, sp, nx, ny, k;
		qaws_curve_evaluate_2d(c, (qaws_scalar)t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2, &e);
		sp = hypot(e.d1.x, e.d1.y);
		nx = e.d1.y / sp; ny = -e.d1.x / sp;
		k = (e.d1.x * e.d2.y - e.d1.y * e.d2.x) / (sp * sp * sp);
		xy[2 * j] = vx(v, e.position.x + delta * nx);
		xy[2 * j + 1] = vy(v, e.position.y + delta * ny);
		if (1 + delta * k < 0)
		{
			ev[2 * ne] = vx(v, e.position.x - nx / k);
			ev[2 * ne + 1] = vy(v, e.position.y - ny / k);
			ne++;
		}
		else if (ne)
		{
			svg_polyline(s, ev, (int)ne, "#8250df", 1.2, 0.8, 1);
			ne = 0;
		}
	}
	if (ne) svg_polyline(s, ev, (int)ne, "#8250df", 1.2, 0.8, 1);
	svg_polyline(s, xy, N + 1, "#8c959f", 1.0, 0.9, 1);
}

static void co_cusps(svg* s, viewport const* v, qaws_curve const* c, double delta, double const* t, unsigned int n)
{
	unsigned int i;
	for (i = 0; i < n; i++)
	{
		qaws_eval_result_2d e;
		double sp;
		qaws_curve_evaluate_2d(c, (qaws_scalar)t[i], QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &e);
		sp = hypot(e.d1.x, e.d1.y);
		svg_circle(s, vx(v, e.position.x + delta * e.d1.y / sp), vy(v, e.position.y - delta * e.d1.x / sp), 4.2, "#24292f", "#ffffff");
	}
}

static qaws_curve* co_bezier(double const* xy, unsigned int n)
{
	qaws_scalar p[8];
	qaws_bezier_desc d;
	qaws_curve* c = NULL;
	unsigned int i;
	for (i = 0; i < 2 * n; i++) p[i] = (qaws_scalar)xy[i];
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D; d.degree = n - 1; d.control_points = p; d.control_point_count = n;
	qaws_curve_create_bezier(&d, &c);
	return c;
}

/* a closed C2 star: the periodic uniform cubic B-spline of m points at
   alternating radii, as Bezier spans in one clamped B-spline */
static qaws_curve* co_star(unsigned int m, double r0, double r1, double cx, double cy)
{
	double Q[32][2];
	qaws_scalar cp[2 * (3 * 32 + 1)], kn[3 * 32 + 5];
	qaws_bspline_desc d;
	qaws_curve* c = NULL;
	unsigned int i, k, n = 0, kc = 0;
	for (i = 0; i < m; i++)
	{
		double a = 2 * CO_PI * i / m + CO_PI / 2, r = (i & 1) ? r1 : r0;
		Q[i][0] = cx + r * cos(a); Q[i][1] = cy + r * sin(a);
	}
	for (i = 0; i < m; i++)
	{
		double const *qm = Q[(i + m - 1) % m], *q0 = Q[i], *q1 = Q[(i + 1) % m], *q2 = Q[(i + 2) % m];
		if (i == 0)
		{
			for (k = 0; k < 2; k++) cp[2 * n + k] = (qaws_scalar)((qm[k] + 4 * q0[k] + q1[k]) / 6);
			n++;
		}
		for (k = 0; k < 2; k++) cp[2 * n + k] = (qaws_scalar)((2 * q0[k] + q1[k]) / 3);
		n++;
		for (k = 0; k < 2; k++) cp[2 * n + k] = (qaws_scalar)((q0[k] + 2 * q1[k]) / 3);
		n++;
		for (k = 0; k < 2; k++) cp[2 * n + k] = (qaws_scalar)((q0[k] + 4 * q1[k] + q2[k]) / 6);
		n++;
	}
	cp[2 * (n - 1)] = cp[0];
	cp[2 * (n - 1) + 1] = cp[1];
	for (i = 0; i < 4; i++) kn[kc++] = 0;
	for (i = 1; i < m; i++) { kn[kc++] = (qaws_scalar)i; kn[kc++] = (qaws_scalar)i; kn[kc++] = (qaws_scalar)i; }
	for (i = 0; i < 4; i++) kn[kc++] = (qaws_scalar)m;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D; d.degree = 3; d.control_points = cp; d.control_point_count = n;
	d.knots = kn; d.knot_count = kc;
	qaws_curve_create_bspline(&d, &c);
	return c;
}

static qaws_scalar co_const(void* user, unsigned int group, unsigned int path, qaws_vec2 point, qaws_vec2 normal)
{
	(void)group; (void)path; (void)point; (void)normal;
	return *(qaws_scalar const*)user;
}

static void co_view(viewport* v, int x0, int y0, int w, int h, double cx, double cy, double span)
{
	v->x0 = x0; v->y0 = y0; v->w = w; v->h = h;
	v->xmin = cx - span / 2; v->xmax = cx + span / 2;
	v->ymin = cy - span / 2 * h / w; v->ymax = cy + span / 2 * h / w;
}

static void demo_curve_offset(void)
{
	int const cols = 3, rows = 2, pw = 420, ph = 360, gap = 16, top = 84;
	int const W = cols * pw + (cols + 1) * gap, H = top + rows * (ph + gap) + 8;
	svg s;
	int k;
	double t0;
	char label[240];
	static double const par[6] = { -1, 1, 0, -1, 1, 1 };
	static double const sxy[8] = { 0, 0, 1, 2, 2, -2, 3, 0 };
	qaws_curve* parabola = co_bezier(par, 3);
	qaws_curve* scurve = co_bezier(sxy, 4);
	qaws_curve* star = co_star(10, 3, 0.6, 0, 0);
	if (!svg_open(&s, "showcase/b2d6_curve_offset.svg", W, H, "Offsets of curves: cusps, evolutes, G1 cubics",
		"After Levien and Uguray: Euler-spiral cuts bracket each cusp (1 + d k = 0), the evolute replaces each backward run, G1 cubics fit the offset. Dashed: plain parallel curve; purple: evolute; dots: cusps"))
		return;
	t0 = (double)clock() / CLOCKS_PER_SEC;
	for (k = 0; k < 6; k++)
	{
		viewport v;
		int x0 = gap + (k % cols) * (pw + gap), y0 = top + (k / cols) * (ph + gap);
		if (k == 0)
		{
			/* the parabola offset 1 toward its focus */
			qaws_curve* ch = NULL;
			double cusps[8];
			unsigned int nc = 0, nk = 0, cub = 0;
			co_view(&v, x0, y0, pw, ph, 0, 0.55, 2.9);
			svg_panel(&s, &v, "parabola, offset 1 toward its focus");
			co_parallel(&s, &v, parabola, -1);
			co_curve(&s, &v, parabola, "#57606a", 2.0, 0.9, 0);
			qaws_internal_offset_curve(parabola, -1, 2e-6, &ch, &nk, cusps, 8, &nc);
			if (ch)
			{
				co_curve(&s, &v, ch, "#cf222e", 2.2, 0.95, 0);
				cub = co_joints(&s, &v, ch, "#cf222e", 2.6);
			}
			co_cusps(&s, &v, parabola, -1, cusps, nc);
			sprintf(label, "cusps t = %.9f, %.9f (Mathematica: same digits); %u cubics",
				nc > 0 ? cusps[0] : 0, nc > 1 ? cusps[1] : 0, cub);
			svg_text(&s, x0 + 10, y0 + ph - 12, 11, "#57606a", "start", label);
			qaws_curve_destroy(ch);
		}
		else if (k == 1)
		{
			/* the S at 0.6 on both sides */
			unsigned int side, cub[2] = { 0, 0 }, ncs[2] = { 0, 0 };
			co_view(&v, x0, y0, pw, ph, 1.5, 0, 4.6);
			svg_panel(&s, &v, "cubic S (radius 0.42 at its bends), offset 0.6 on both sides");
			for (side = 0; side < 2; side++)
			{
				qaws_curve* ch = NULL;
				double cusps[8], d = side ? -0.6 : 0.6;
				unsigned int nk = 0;
				co_parallel(&s, &v, scurve, d);
				qaws_internal_offset_curve(scurve, d, 4e-6, &ch, &nk, cusps, 8, &ncs[side]);
				if (ch)
				{
					co_curve(&s, &v, ch, side ? "#bc4c00" : "#0969da", 2.2, 0.95, 0);
					cub[side] = co_joints(&s, &v, ch, side ? "#bc4c00" : "#0969da", 2.6);
				}
				co_cusps(&s, &v, scurve, d, cusps, ncs[side]);
				qaws_curve_destroy(ch);
			}
			co_curve(&s, &v, scurve, "#57606a", 2.0, 0.9, 0);
			sprintf(label, "right: %u cusps, %u cubics; left: %u cusps, %u cubics (cusps match Mathematica to 1e-10)", ncs[0], cub[0], ncs[1], cub[1]);
			svg_text(&s, x0 + 10, y0 + ph - 12, 11, "#57606a", "start", label);
		}
		else if (k == 2 || k == 3)
		{
			/* the star grown by 1 (valleys of radius 0.85) and shrunk by 0.45 */
			double d = k == 2 ? 1.0 : -0.45;
			qaws_curve const* cc = star;
			qaws_path_2d p;
			qaws_clip_result* r = NULL;
			qaws_curve* ch = NULL;
			double cusps[64];
			unsigned int nc = 0, nk = 0, cub = 0;
			p.curves = &cc; p.curve_count = 1; p.closed = 1;
			co_view(&v, x0, y0, pw, ph, 0, 0, 9.2);
			svg_panel(&s, &v, k == 2 ? "star grown by 1: raw chain, then the union" : "star shrunk by 0.45: raw chain, then the union");
			qaws_internal_offset_curve(star, d, 1e-5, &ch, &nk, cusps, 64, &nc);
			if (qaws_offset_paths(&p, 1, (qaws_scalar)d, QAWS_JOIN_ROUND, QAWS_END_POLYGON, &r) == QAWS_STATUS_OK)
				co_result(&s, &v, r, k == 2 ? "#c8e1ff" : "#dafbe1", k == 2 ? "#0969da" : "#1a7f37", &cub, 0);
			if (ch) co_curve(&s, &v, ch, "#cf222e", 0.9, 0.85, 0);
			co_cusps(&s, &v, star, d, cusps, nc);
			co_curve(&s, &v, star, "#57606a", 1.6, 0.9, 0);
			sprintf(label, "%u cusps; the result: %u cubics after the union", nc, cub);
			svg_text(&s, x0 + 10, y0 + ph - 12, 11, "#57606a", "start", label);
			qaws_clip_result_destroy(r);
			qaws_curve_destroy(ch);
		}
		else if (k == 4)
		{
			/* the S stroked with round ends */
			qaws_curve const* cc = scurve;
			qaws_path_2d p;
			qaws_clip_result* r = NULL;
			unsigned int cub = 0;
			p.curves = &cc; p.curve_count = 1; p.closed = 0;
			co_view(&v, x0, y0, pw, ph, 1.5, 0, 4.6);
			svg_panel(&s, &v, "the S stroked 0.6, round ends: cusps on both sides");
			co_parallel(&s, &v, scurve, 0.6);
			co_parallel(&s, &v, scurve, -0.6);
			if (qaws_offset_paths(&p, 1, (qaws_scalar)0.6, QAWS_JOIN_ROUND, QAWS_END_ROUND, &r) == QAWS_STATUS_OK)
				co_result(&s, &v, r, "#fff1e5", "#bc4c00", &cub, 1);
			co_curve(&s, &v, scurve, "#57606a", 2.0, 0.9, 0);
			sprintf(label, "one region, %u cubics and arcs; each point within 0.6 of the S (grid oracle, test 90)", cub);
			svg_text(&s, x0 + 10, y0 + ph - 12, 11, "#57606a", "start", label);
			qaws_clip_result_destroy(r);
		}
		else
		{
			/* the old fit (cubics through offset points, still used under a
			   delta callback) beside the new one, same offset */
			qaws_curve* a = co_star(8, 2.2, 1.5, -2.45, 0);
			qaws_curve* b = co_star(8, 2.2, 1.5, 2.45, 0);
			qaws_curve const *ca = a, *cb = b;
			qaws_path_2d pa, pb;
			qaws_offset_group g;
			qaws_offset_desc d;
			qaws_clip_result *ra = NULL, *rb = NULL;
			qaws_scalar delta = (qaws_scalar)0.5;
			unsigned int na = 0, nb = 0;
			pa.curves = &ca; pa.curve_count = 1; pa.closed = 1;
			pb.curves = &cb; pb.curve_count = 1; pb.closed = 1;
			co_view(&v, x0, y0, pw, ph, 0, 0, 10.2);
			svg_panel(&s, &v, "the same offset (0.5): old fit | Euler cuts + G1 cubics");
			memset(&g, 0, sizeof(g));
			g.paths = &pa; g.path_count = 1; g.join_type = QAWS_JOIN_ROUND; g.end_type = QAWS_END_POLYGON;
			memset(&d, 0, sizeof(d));
			d.groups = &g; d.group_count = 1; d.delta = delta; d.delta_fn = co_const; d.delta_user = &delta;
			if (qaws_offset_execute(&d, &ra) == QAWS_STATUS_OK)
				co_result(&s, &v, ra, NULL, "#bf8700", &na, 1);
			if (qaws_offset_paths(&pb, 1, delta, QAWS_JOIN_ROUND, QAWS_END_POLYGON, &rb) == QAWS_STATUS_OK)
				co_result(&s, &v, rb, NULL, "#1a7f37", &nb, 1);
			co_curve(&s, &v, a, "#57606a", 1.4, 0.9, 0);
			co_curve(&s, &v, b, "#57606a", 1.4, 0.9, 0);
			sprintf(label, "same tolerance (1e-6 of the extent): %u cubics before, %u after, G1 at every joint", na, nb);
			svg_text(&s, x0 + 10, y0 + ph - 12, 11, "#57606a", "start", label);
			qaws_clip_result_destroy(ra);
			qaws_clip_result_destroy(rb);
			qaws_curve_destroy(a);
			qaws_curve_destroy(b);
		}
	}
	svg_close(&s);
	printf("    curve offset: 6 panels in %.3f s\n", (double)clock() / CLOCKS_PER_SEC - t0);
	qaws_curve_destroy(parabola);
	qaws_curve_destroy(scurve);
	qaws_curve_destroy(star);
}
