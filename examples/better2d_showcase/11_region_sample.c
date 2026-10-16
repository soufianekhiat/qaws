/* Figure 11: sampling regions of curves (qaws_region_sampler): shapes built
   by chains of Boolean operations on paths that mix lines, circle arcs and
   cubic Beziers, sampled uniformly and in proportion to densities (a bump,
   a gradient, stripes, a ring), with every pattern. The sampler is fed the
   Boolean results as they come out of qaws_clip. */

/* a region: a clip result and its paths */
typedef struct rsr
{
	qaws_clip_result* r;
	qaws_path_2d* p;
	unsigned int n;
} rsr;

static void rsr_take(rsr* g, qaws_clip_result* r)
{
	unsigned int i;
	g->r = r;
	g->n = r ? qaws_clip_result_get_path_count(r) : 0;
	g->p = (qaws_path_2d*)malloc(sizeof(qaws_path_2d) * (g->n + 1));
	for (i = 0; i < g->n; i++) qaws_clip_result_get_path(r, i, &g->p[i]);
}

static void rsr_free(rsr* g)
{
	qaws_clip_result_destroy(g->r);
	free(g->p);
	memset(g, 0, sizeof(*g));
}

/* a op b, as a new region (a and b kept) */
static rsr rsr_op(qaws_clip_type op, qaws_path_2d const* a, unsigned int na, qaws_path_2d const* b, unsigned int nb)
{
	rsr g;
	qaws_clip_result* r = NULL;
	qaws_clip_boolean(op, QAWS_FILL_NON_ZERO, a, na, b, nb, &r);
	rsr_take(&g, r);
	return g;
}

/* replace a by a op b */
static void rsr_apply(rsr* a, qaws_clip_type op, qaws_path_2d const* b, unsigned int nb)
{
	rsr g = rsr_op(op, a->p, a->n, b, nb);
	rsr_free(a);
	*a = g;
}

/* curves kept alive for the paths built from them */
static qaws_curve* rs_keep_curves[512];
static qaws_curve const* rs_keep_views[512];
static unsigned int rs_nkeep = 0;

static qaws_curve const* rs_hold(qaws_curve* c)
{
	rs_keep_curves[rs_nkeep] = c;
	rs_keep_views[rs_nkeep] = c;
	return rs_keep_views[rs_nkeep++];
}

static qaws_curve const* rs_line(double x0, double y0, double x1, double y1)
{
	qaws_scalar p[4];
	qaws_bezier_desc d;
	qaws_curve* c = NULL;
	p[0] = (qaws_scalar)x0; p[1] = (qaws_scalar)y0; p[2] = (qaws_scalar)x1; p[3] = (qaws_scalar)y1;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D; d.degree = 1; d.control_points = p; d.control_point_count = 2;
	qaws_curve_create_bezier(&c, &d);
	return rs_hold(c);
}

static qaws_curve const* rs_cubic(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3)
{
	qaws_scalar p[8];
	qaws_bezier_desc d;
	qaws_curve* c = NULL;
	p[0] = (qaws_scalar)x0; p[1] = (qaws_scalar)y0; p[2] = (qaws_scalar)x1; p[3] = (qaws_scalar)y1;
	p[4] = (qaws_scalar)x2; p[5] = (qaws_scalar)y2; p[6] = (qaws_scalar)x3; p[7] = (qaws_scalar)y3;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D; d.degree = 3; d.control_points = p; d.control_point_count = 4;
	qaws_curve_create_bezier(&c, &d);
	return rs_hold(c);
}

/* the arc of circle (cx, cy, r) from angle a0 to a1 (radians, a1 < a0 runs clockwise) */
static qaws_curve const* rs_arc(double cx, double cy, double r, double a0, double a1)
{
	qaws_arc_segment g;
	qaws_arc_desc d;
	qaws_curve* c = NULL;
	memset(&g, 0, sizeof(g));
	g.center[0] = (qaws_scalar)cx; g.center[1] = (qaws_scalar)cy; g.radius = (qaws_scalar)r;
	g.angle_start = (qaws_scalar)a0; g.angle_end = (qaws_scalar)a1;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D; d.segments = &g; d.segment_count = 1;
	qaws_curve_create_arc(&c, &d);
	return rs_hold(c);
}

static qaws_path_2d rs_path(qaws_curve const* const* cs, unsigned int n)
{
	qaws_path_2d p;
	qaws_curve const** keep = &rs_keep_views[rs_nkeep];
	unsigned int i;
	/* the views copied after the held ones, so the path points at stable storage */
	for (i = 0; i < n; i++) rs_keep_views[rs_nkeep + i] = cs[i];
	rs_nkeep += n;
	p.curves = keep; p.curve_count = n; p.closed = 1;
	return p;
}

static qaws_path_2d rs_disc(double cx, double cy, double r)
{
	qaws_curve const* c = rs_arc(cx, cy, r, 0, 2 * PI);
	return rs_path(&c, 1);
}

/* a rectangle with rounded corners: four lines and four quarter arcs */
static qaws_path_2d rs_rounded_rect(double x0, double y0, double x1, double y1, double rr)
{
	qaws_curve const* c[8];
	c[0] = rs_line(x0 + rr, y0, x1 - rr, y0);
	c[1] = rs_arc(x1 - rr, y0 + rr, rr, -PI / 2, 0);
	c[2] = rs_line(x1, y0 + rr, x1, y1 - rr);
	c[3] = rs_arc(x1 - rr, y1 - rr, rr, 0, PI / 2);
	c[4] = rs_line(x1 - rr, y1, x0 + rr, y1);
	c[5] = rs_arc(x0 + rr, y1 - rr, rr, PI / 2, PI);
	c[6] = rs_line(x0, y1 - rr, x0, y0 + rr);
	c[7] = rs_arc(x0 + rr, y0 + rr, rr, PI, 1.5 * PI);
	return rs_path(c, 8);
}

static qaws_path_2d rs_polygon(double const* xy, unsigned int n)
{
	qaws_curve const* c[64];
	unsigned int i;
	for (i = 0; i < n; i++)
		c[i] = rs_line(xy[2 * i], xy[2 * i + 1], xy[2 * ((i + 1) % n)], xy[2 * ((i + 1) % n) + 1]);
	return rs_path(c, n);
}

/* ------------------------------------------------------------------------ */
/*  The shapes                                                              */
/* ------------------------------------------------------------------------ */

/* a gear: a disc with 12 trapezoid teeth (lines), minus a hub, minus three
   bolt holes, minus a keyway */
static rsr rs_gear(void)
{
	double xy[2 * 48];
	qaws_path_2d teeth, body, hub, bolts[3], key;
	unsigned int i, k = 0;
	rsr g;
	for (i = 0; i < 12; i++)
	{
		double a = 2 * PI * i / 12, w = PI / 12;
		xy[k++] = 0.8 * cos(a - w); xy[k++] = 0.8 * sin(a - w);
		xy[k++] = 1.0 * cos(a - 0.45 * w); xy[k++] = 1.0 * sin(a - 0.45 * w);
		xy[k++] = 1.0 * cos(a + 0.45 * w); xy[k++] = 1.0 * sin(a + 0.45 * w);
		xy[k++] = 0.8 * cos(a + w); xy[k++] = 0.8 * sin(a + w);
	}
	teeth = rs_polygon(xy, 48);
	body = rs_disc(0, 0, 0.84);
	hub = rs_disc(0, 0, 0.22);
	for (i = 0; i < 3; i++) bolts[i] = rs_disc(0.5 * cos(2 * PI * i / 3 + 0.5), 0.5 * sin(2 * PI * i / 3 + 0.5), 0.11);
	{
		double kxy[8] = { -0.06, 0.15, 0.06, 0.15, 0.06, 0.3, -0.06, 0.3 };
		key = rs_polygon(kxy, 4);
	}
	g = rsr_op(QAWS_CLIP_UNION, &teeth, 1, &body, 1);
	rsr_apply(&g, QAWS_CLIP_DIFFERENCE, &hub, 1);
	rsr_apply(&g, QAWS_CLIP_DIFFERENCE, bolts, 3);
	rsr_apply(&g, QAWS_CLIP_UNION, &key, 1);
	{
		/* the keyway cut into the hub's rim */
		double kxy[8] = { -0.05, 0.12, 0.05, 0.12, 0.05, 0.26, -0.05, 0.26 };
		qaws_path_2d cut = rs_polygon(kxy, 4);
		rsr_apply(&g, QAWS_CLIP_DIFFERENCE, &cut, 1);
	}
	return g;
}

/* a wrench: a rounded bar united with two discs, minus a slot-shaped jaw
   (lines and arcs) at one end and a hexagon at the other, a Bezier grip cut */
static rsr rs_wrench(void)
{
	qaws_path_2d bar = rs_rounded_rect(-1.0, -0.14, 1.0, 0.14, 0.1);
	qaws_path_2d end0 = rs_disc(-1.05, 0, 0.38), end1 = rs_disc(1.05, 0, 0.34);
	qaws_path_2d jaw, hex, grip;
	rsr g;
	{
		qaws_curve const* c[4];
		c[0] = rs_line(-1.6, -0.12, -1.1, -0.12);
		c[1] = rs_arc(-1.1, 0, 0.12, -PI / 2, PI / 2);
		c[2] = rs_line(-1.1, 0.12, -1.6, 0.12);
		c[3] = rs_line(-1.6, 0.12, -1.6, -0.12);
		jaw = rs_path(c, 4);
	}
	{
		double xy[12];
		unsigned int i;
		for (i = 0; i < 6; i++) { xy[2 * i] = 1.05 + 0.17 * cos(PI / 6 + i * PI / 3); xy[2 * i + 1] = 0.17 * sin(PI / 6 + i * PI / 3); }
		hex = rs_polygon(xy, 6);
	}
	{
		qaws_curve const* c[2];
		c[0] = rs_cubic(-0.4, 0.2, -0.2, -0.05, 0.2, -0.05, 0.4, 0.2);
		c[1] = rs_line(0.4, 0.2, -0.4, 0.2);
		grip = rs_path(c, 2);
	}
	g = rsr_op(QAWS_CLIP_UNION, &bar, 1, &end0, 1);
	rsr_apply(&g, QAWS_CLIP_UNION, &end1, 1);
	rsr_apply(&g, QAWS_CLIP_DIFFERENCE, &jaw, 1);
	rsr_apply(&g, QAWS_CLIP_DIFFERENCE, &hex, 1);
	rsr_apply(&g, QAWS_CLIP_DIFFERENCE, &grip, 1);
	return g;
}

/* a moon and a star: a crescent (disc minus disc) united with a star
   polygon, xor a disc over both, minus a Bezier comet tail */
static rsr rs_moon_star(void)
{
	qaws_path_2d moon = rs_disc(-0.3, 0, 0.9), bite = rs_disc(0.05, 0.15, 0.75), ring = rs_disc(0.55, -0.1, 0.32);
	qaws_path_2d star, tail;
	double xy[20];
	unsigned int i;
	rsr g;
	for (i = 0; i < 10; i++)
	{
		double r = (i % 2) ? 0.2 : 0.5, t = PI / 2 + PI * i / 5;
		xy[2 * i] = 0.55 + r * cos(t); xy[2 * i + 1] = 0.35 + r * sin(t);
	}
	star = rs_polygon(xy, 10);
	{
		qaws_curve const* c[2];
		c[0] = rs_cubic(-1.3, -0.9, -0.8, -0.5, -0.5, -0.75, -0.2, -0.55);
		c[1] = rs_cubic(-0.2, -0.55, -0.6, -0.95, -0.9, -0.75, -1.3, -0.9);
		tail = rs_path(c, 2);
	}
	g = rsr_op(QAWS_CLIP_DIFFERENCE, &moon, 1, &bite, 1);
	rsr_apply(&g, QAWS_CLIP_UNION, &star, 1);
	rsr_apply(&g, QAWS_CLIP_XOR, &ring, 1);
	rsr_apply(&g, QAWS_CLIP_DIFFERENCE, &tail, 1);
	return g;
}

/* a cheese: a rounded square minus seven discs, xor an ellipse breaking out
   of its side, minus a Bezier bite at a corner */
static rsr rs_cheese(void)
{
	qaws_path_2d sq = rs_rounded_rect(-1, -0.8, 1, 0.8, 0.25), holes[7], bite;
	static double const hc[7][3] = { { -0.55, 0.35, 0.18 }, { 0.1, 0.45, 0.12 }, { 0.6, 0.2, 0.2 }, { -0.2, -0.15, 0.22 },
		{ 0.45, -0.45, 0.14 }, { -0.7, -0.5, 0.11 }, { 0.15, 0.05, 0.08 } };
	qaws_curve* e = NULL;
	qaws_path_2d ell;
	unsigned int i;
	rsr g;
	for (i = 0; i < 7; i++) holes[i] = rs_disc(hc[i][0], hc[i][1], hc[i][2]);
	qaws_curve_create_ellipse_2d((qaws_scalar)1.0, (qaws_scalar)-0.1, (qaws_scalar)0.45, (qaws_scalar)0.25, (qaws_scalar)0.4, &e);
	{
		qaws_curve const* ec = rs_hold(e);
		ell = rs_path(&ec, 1);
	}
	{
		qaws_curve const* c[3];
		c[0] = rs_line(-1.1, 0.95, -1.1, 0.35);
		c[1] = rs_cubic(-1.1, 0.35, -0.7, 0.45, -0.6, 0.7, -0.5, 0.95);
		c[2] = rs_line(-0.5, 0.95, -1.1, 0.95);
		bite = rs_path(c, 3);
	}
	g = rsr_op(QAWS_CLIP_DIFFERENCE, &sq, 1, holes, 7);
	rsr_apply(&g, QAWS_CLIP_XOR, &ell, 1);
	rsr_apply(&g, QAWS_CLIP_DIFFERENCE, &bite, 1);
	return g;
}

/* a chain: three annuli (disc minus disc) united, intersected with a
   band bounded by two cubics, then a fourth link xor-ed through */
static rsr rs_chain(void)
{
	rsr g, a;
	qaws_path_2d band, o, in, o4, in4;
	unsigned int i;
	g.r = NULL; g.p = NULL; g.n = 0;
	for (i = 0; i < 3; i++)
	{
		double cx = -0.9 + 0.75 * i;
		o = rs_disc(cx, 0, 0.5);
		in = rs_disc(cx, 0, 0.33);
		a = rsr_op(QAWS_CLIP_DIFFERENCE, &o, 1, &in, 1);
		if (!g.r) g = a;
		else { rsr_apply(&g, QAWS_CLIP_UNION, a.p, a.n); rsr_free(&a); }
	}
	{
		qaws_curve const* c[4];
		c[0] = rs_cubic(-1.6, -0.45, -0.8, -0.2, 0.0, -0.6, 0.8, -0.3);
		c[1] = rs_line(0.8, -0.3, 0.8, 0.45);
		c[2] = rs_cubic(0.8, 0.45, 0.0, 0.7, -0.8, 0.25, -1.6, 0.45);
		c[3] = rs_line(-1.6, 0.45, -1.6, -0.45);
		band = rs_path(c, 4);
	}
	rsr_apply(&g, QAWS_CLIP_INTERSECTION, &band, 1);
	o4 = rs_disc(0.85, 0.05, 0.42);
	in4 = rs_disc(0.85, 0.05, 0.27);
	a = rsr_op(QAWS_CLIP_DIFFERENCE, &o4, 1, &in4, 1);
	rsr_apply(&g, QAWS_CLIP_XOR, a.p, a.n);
	rsr_free(&a);
	return g;
}

/* a leaf: two cubic arcs (a vesica) united with a stem (lines), minus veins
   (thin triangles), minus a round bite (arc) out of its edge */
static rsr rs_leaf(void)
{
	qaws_path_2d leaf, stem, bite, veins[4];
	unsigned int i;
	rsr g;
	{
		qaws_curve const* c[2];
		c[0] = rs_cubic(-1.2, -0.3, -0.6, 0.9, 0.6, 0.9, 1.2, 0.3);
		c[1] = rs_cubic(1.2, 0.3, 0.6, -0.7, -0.6, -1.0, -1.2, -0.3);
		leaf = rs_path(c, 2);
	}
	{
		double xy[8] = { -1.55, -0.55, -1.1, -0.28, -1.12, -0.24, -1.57, -0.51 };
		stem = rs_polygon(xy, 4);
	}
	for (i = 0; i < 4; i++)
	{
		double x = -0.65 + 0.42 * i, xy[6];
		xy[0] = x - 0.05; xy[1] = -0.02 - 0.03 * i; xy[2] = x + 0.32; xy[3] = 0.42 - 0.05 * i; xy[4] = x + 0.02; xy[5] = 0.0 - 0.03 * i;
		veins[i] = rs_polygon(xy, 3);
	}
	bite = rs_disc(0.55, -0.6, 0.25);
	g = rsr_op(QAWS_CLIP_UNION, &leaf, 1, &stem, 1);
	rsr_apply(&g, QAWS_CLIP_DIFFERENCE, veins, 4);
	rsr_apply(&g, QAWS_CLIP_DIFFERENCE, &bite, 1);
	return g;
}

/* ------------------------------------------------------------------------ */
/*  Densities                                                               */
/* ------------------------------------------------------------------------ */

typedef struct rs_dens { double cx, cy, w, k; } rs_dens;

static qaws_scalar rs_bump(void* u, qaws_scalar x, qaws_scalar y)
{
	rs_dens const* d = (rs_dens const*)u;
	return (qaws_scalar)(0.05 + exp(-((x - d->cx) * (x - d->cx) + (y - d->cy) * (y - d->cy)) / (d->w * d->w)));
}

static qaws_scalar rs_gradient(void* u, qaws_scalar x, qaws_scalar y)
{
	rs_dens const* d = (rs_dens const*)u;
	(void)y;
	return (qaws_scalar)exp(d->k * (x - d->cx));
}

static qaws_scalar rs_stripes(void* u, qaws_scalar x, qaws_scalar y)
{
	rs_dens const* d = (rs_dens const*)u;
	double s = sin(d->k * (0.8 * x + 0.6 * y));
	return (qaws_scalar)(0.03 + s * s * s * s);
}

static qaws_scalar rs_ring(void* u, qaws_scalar x, qaws_scalar y)
{
	rs_dens const* d = (rs_dens const*)u;
	double r = sqrt((x - d->cx) * (x - d->cx) + (y - d->cy) * (y - d->cy)) - d->k;
	return (qaws_scalar)(0.04 + exp(-r * r / (d->w * d->w)));
}

/* ------------------------------------------------------------------------ */
/*  Drawing                                                                 */
/* ------------------------------------------------------------------------ */

static void rs_draw_region(svg* s, viewport const* v, rsr const* g)
{
	unsigned int i, k, j;
	fprintf(s->f, "<path fill=\"#ddf4ff\" fill-rule=\"nonzero\" stroke=\"#0969da\" stroke-width=\"1.1\" d=\"");
	for (i = 0; i < g->n; i++)
	{
		for (k = 0; k < g->p[i].curve_count; k++)
		{
			qaws_curve const* c = g->p[i].curves[k];
			qaws_range rg = qaws_curve_get_parameter_range(c);
			unsigned int m = 24 * (qaws_curve_get_span_count(c) ? qaws_curve_get_span_count(c) : 1);
			for (j = (k == 0 ? 0 : 1); j <= m; j++)
			{
				qaws_eval_result_2d e;
				qaws_curve_evaluate_2d(&e, c, (qaws_scalar)(rg.min_value + (rg.max_value - rg.min_value) * j / m), QAWS_EVAL_FLAG_POSITION);
				fprintf(s->f, "%s%.2f,%.2f ", (j == 0 && k == 0) ? "M" : "L", vx(v, e.position.x), vy(v, e.position.y));
			}
		}
		fprintf(s->f, "Z ");
	}
	fprintf(s->f, "\"/>\n");
}

static void rs_draw_points(svg* s, viewport const* v, qaws_vec2 const* p, unsigned int n, char const* color)
{
	unsigned int i;
	fprintf(s->f, "<g fill=\"%s\">", color);
	for (i = 0; i < n; i++)
		fprintf(s->f, "<circle cx=\"%.2f\" cy=\"%.2f\" r=\"1.15\"/>", vx(v, p[i].x), vy(v, p[i].y));
	fprintf(s->f, "</g>\n");
}

static void demo_region_sample(void)
{
	typedef rsr (*maker)(void);
	static maker const makers[6] = { rs_gear, rs_wrench, rs_moon_star, rs_cheese, rs_chain, rs_leaf };
	static char const* const names[6] = {
		"gear: teeth (lines) + disc - hub - 3 bolts + key - keyway",
		"wrench: rounded bar + 2 discs - slot jaw - hexagon - Bezier grip",
		"moon and star: disc - disc + star, xor disc, - Bezier tail",
		"cheese: rounded square - 7 discs, xor ellipse, - Bezier bite",
		"chain: 3 annuli united, &amp; Bezier band, xor a 4th annulus",
		"leaf: 2 cubics + stem - 4 veins - bite" };
	static double const box[6][4] = { { -1.1, 1.1, -1.1, 1.1 }, { -1.65, 1.5, -0.6, 0.6 }, { -1.35, 1.15, -1.0, 1.0 },
		{ -1.15, 1.5, -0.95, 0.95 }, { -1.7, 1.35, -0.65, 0.65 }, { -1.65, 1.3, -1.0, 0.85 } };
	/* the columns: a density (0 uniform, 1 bump, 2 gradient, 3 stripes, 4 ring) and a pattern */
	static unsigned int const dens[5] = { 0, 1, 2, 3, 4 };
	static unsigned int const kinds[5] = { QAWS_REGION_SAMPLE_SOBOL | QAWS_REGION_SAMPLE_BOX, QAWS_REGION_SAMPLE_BLUE_NOISE,
		QAWS_REGION_SAMPLE_SOBOL, QAWS_REGION_SAMPLE_STRATIFIED, QAWS_REGION_SAMPLE_R2 | QAWS_REGION_SAMPLE_BOX };
	static char const* const cols[5] = { "uniform, Sobol over the box", "bump, blue noise", "gradient e^(2.2 x), Sobol",
		"stripes sin^4, stratified", "ring, R2 over the box" };
	unsigned int n = 1800, r, c;
	qaws_vec2* pts = (qaws_vec2*)malloc(sizeof(qaws_vec2) * n);
	svg s;
	double pw = 300, ph = 190;
	clock_t t0 = clock();
	if (!svg_open(&s, "showcase/b2d11_region_sample.svg", (int)(40 + 5 * (pw + 10)), (int)(110 + 6 * (ph + 32)),
		"Sampling regions of curves: uniform and importance sampling, no mesh",
		"qaws_region_sampler on Boolean results (lines, arcs, cubic Beziers), 1800 points per panel"))
	{
		free(pts);
		return;
	}
	for (c = 0; c < 5; c++)
		svg_text(&s, 30 + c * (pw + 10) + pw / 2, 92, 13, "#24292f", "middle", cols[c]);
	for (r = 0; r < 6; r++)
	{
		rsr g = makers[r]();
		qaws_region_sampler_desc d;
		double bx = box[r][1] - box[r][0], by = box[r][3] - box[r][2];
		rs_dens dd;
		dd.cx = box[r][0] + 0.68 * bx; dd.cy = box[r][2] + 0.4 * by; dd.w = 0.28 * (bx < by ? bx : by) + 0.15; dd.k = 0;
		svg_text(&s, 30, 110 + r * (ph + 32) + 14, 13, "#24292f", "start", names[r]);
		for (c = 0; c < 5; c++)
		{
			viewport v;
			qaws_region_sampler* smp = NULL;
			double sc = (pw - 16) / bx < (ph - 12) / by ? (pw - 16) / bx : (ph - 12) / by;
			v.w = bx * sc; v.h = by * sc;
			v.x0 = 30 + c * (pw + 10) + (pw - v.w) / 2;
			v.y0 = 110 + r * (ph + 32) + 22 + (ph - v.h) / 2;
			v.xmin = box[r][0]; v.xmax = box[r][1]; v.ymin = box[r][2]; v.ymax = box[r][3];
			memset(&d, 0, sizeof(d));
			d.paths = g.p; d.path_count = g.n; d.fill_rule = QAWS_FILL_NON_ZERO;
			switch (dens[c])
			{
			case 1: d.density = rs_bump; break;
			case 2: dd.k = 2.2; d.density = rs_gradient; break;
			case 3: dd.k = 7.0; d.density = rs_stripes; break;
			case 4: dd.k = 0.45 * (bx < by ? bx : by); d.density = rs_ring; dd.cx = box[r][0] + 0.5 * bx; dd.cy = box[r][2] + 0.5 * by; dd.w = 0.12 * (bx < by ? bx : by); break;
			default: break;
			}
			d.density_user = &dd;
			rs_draw_region(&s, &v, &g);
			if (qaws_region_sampler_create(&d, &smp) == QAWS_STATUS_OK)
			{
				qaws_region_sampler_generate(smp, n, kinds[c], 7 + r, pts);
				rs_draw_points(&s, &v, pts, n, "#cf222e");
				qaws_region_sampler_destroy(smp);
			}
			/* restore the bump's centre after the ring column moved it */
			dd.cx = box[r][0] + 0.68 * bx; dd.cy = box[r][2] + 0.4 * by; dd.w = 0.28 * (bx < by ? bx : by) + 0.15;
		}
		rsr_free(&g);
	}
	svg_close(&s);
	free(pts);
	{
		unsigned int i;
		for (i = 0; i < rs_nkeep; i++)
			if (rs_keep_curves[i]) qaws_curve_destroy(rs_keep_curves[i]);
	}
	printf("figure 11: region sampling, %.1f s\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
}
