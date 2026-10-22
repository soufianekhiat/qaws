/* Figure 13: curved skeletons (qaws_skeleton on lines and circular arcs).
   Each shape twice: mitered (sigma = 2 pi: reflex corners kept sharp until
   their edges turn tangent) and the medial axis (sigma = pi: the classical
   round offsets). The wavefront at eight times coloured by time, the
   skeleton's conic arcs in black, smooth joints' paths in grey, and the
   arcs that touch the boundary drawn thin (the pure skeleton thick). */

typedef struct csk_seg
{
	int arc;
	double a[5];
} csk_seg;

typedef struct csk_shape
{
	char const* name;
	unsigned int nloops, count[4];
	csk_seg seg[48];
	double box[4];
} csk_shape;

static void csk_line(csk_shape* s, unsigned int* n, double x0, double y0, double x1, double y1)
{
	csk_seg* g = &s->seg[(*n)++];
	g->arc = 0; g->a[0] = x0; g->a[1] = y0; g->a[2] = x1; g->a[3] = y1;
}

static void csk_arc(csk_shape* s, unsigned int* n, double cx, double cy, double r, double a0, double a1)
{
	csk_seg* g = &s->seg[(*n)++];
	g->arc = 1; g->a[0] = cx; g->a[1] = cy; g->a[2] = r; g->a[3] = a0; g->a[4] = a1;
}

static qaws_curve* csk_curve(csk_seg const* g)
{
	qaws_curve* c = NULL;
	if (!g->arc)
	{
		qaws_scalar cp[4];
		qaws_bezier_desc d;
		cp[0] = (qaws_scalar)g->a[0]; cp[1] = (qaws_scalar)g->a[1]; cp[2] = (qaws_scalar)g->a[2]; cp[3] = (qaws_scalar)g->a[3];
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = 1; d.control_points = cp; d.control_point_count = 2;
		qaws_curve_create_bezier(&c, &d);
	}
	else
	{
		qaws_arc_segment s;
		qaws_arc_desc d;
		memset(&s, 0, sizeof(s));
		s.center[0] = (qaws_scalar)g->a[0]; s.center[1] = (qaws_scalar)g->a[1]; s.radius = (qaws_scalar)g->a[2];
		s.angle_start = (qaws_scalar)g->a[3]; s.angle_end = (qaws_scalar)g->a[4];
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.segments = &s; d.segment_count = 1;
		qaws_curve_create_arc(&c, &d);
	}
	return c;
}

static void csk_path_curves(svg* s, viewport const* v, qaws_path_2d const* p, unsigned int np, char const* stroke, double width, char const* fill)
{
	unsigned int i, k, j;
	fprintf(s->f, "<path fill=\"%s\" fill-rule=\"nonzero\" stroke=\"%s\" stroke-width=\"%.2f\" d=\"", fill, stroke, width);
	for (i = 0; i < np; i++)
	{
		for (k = 0; k < p[i].curve_count; k++)
		{
			qaws_curve const* c = p[i].curves[k];
			if (qaws_curve_get_degree(c) == 1 && qaws_curve_get_kind(c) != QAWS_CURVE_KIND_ARC)
			{
				/* a polyline (qaws_clip joins consecutive lines): its points */
				qaws_scalar buf[2 * 256];
				unsigned int m = 0;
				qaws_curve_get_control_points(buf, &m, c, 256);
				for (j = (k == 0 ? 0 : 1); j < m && j < 256; j++)
					fprintf(s->f, "%s%.2f,%.2f ", (j == 0 && k == 0) ? "M" : "L", vx(v, buf[2 * j]), vy(v, buf[2 * j + 1]));
				continue;
			}
			qaws_range rg = qaws_curve_get_parameter_range(c);
			unsigned int m = 48;
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

/* (weight_fn NULL: every weight 1) */
static void csk_draw(svg* s, viewport const* v, csk_shape const* sh, double sigma, qaws_skeleton_weight_fn weight_fn, void* weight_user)
{
	qaws_curve* cs[48];
	qaws_curve const* cv[48];
	qaws_path_2d p[4];
	qaws_skeleton_desc d;
	qaws_skeleton* sk = NULL;
	unsigned int l, i, k, off = 0;
	double tmax;
	char info[128];
	for (l = 0; l < sh->nloops; l++)
	{
		for (i = 0; i < sh->count[l]; i++) { cs[off + i] = csk_curve(&sh->seg[off + i]); cv[off + i] = cs[off + i]; }
		p[l].curves = cv + off; p[l].curve_count = sh->count[l]; p[l].closed = 1;
		off += sh->count[l];
	}
	csk_path_curves(s, v, p, sh->nloops, "none", 0, "#f6f8fa");
	memset(&d, 0, sizeof(d));
	d.paths = p; d.path_count = sh->nloops; d.fill_rule = QAWS_FILL_NON_ZERO; d.splice_angle = (qaws_scalar)sigma;
	d.weight_fn = weight_fn; d.weight_user = weight_user;
	if (qaws_skeleton_create(&sk, &d) != QAWS_STATUS_OK)
	{
		svg_text(s, v->x0 + v->w / 2, v->y0 + v->h / 2, 12, "#cf222e", "middle", "failed");
		for (i = 0; i < off; i++) qaws_curve_destroy(cs[i]);
		return;
	}
	tmax = (double)qaws_skeleton_get_max_time(sk);
	for (k = 1; k <= 8; k++)
	{
		qaws_clip_result* r = NULL;
		char col[32];
		heat((double)k / 9, col);
		if (qaws_skeleton_offset(&r, sk, (qaws_scalar)(tmax * k / 9)) == QAWS_STATUS_OK)
		{
			unsigned int n = qaws_clip_result_get_path_count(r), pi;
			qaws_path_2d* q = (qaws_path_2d*)malloc(sizeof(qaws_path_2d) * (n + 1));
			for (pi = 0; pi < n; pi++) qaws_clip_result_get_path(&q[pi], r, pi);
			csk_path_curves(s, v, q, n, col, 1.0, "none");
			free(q);
		}
		qaws_clip_result_destroy(r);
	}
	for (i = 0; i < qaws_skeleton_get_arc_count(sk); i++)
	{
		qaws_skeleton_arc a;
		qaws_skeleton_get_arc(&a, sk, i);
		fprintf(s->f, "<polyline fill=\"none\" stroke=\"%s\" stroke-width=\"%.1f\"%s points=\"",
			a.smooth ? "#8c959f" : "#1b1f24", a.smooth ? 0.8 : (a.boundary ? 1.0 : 2.0), a.smooth ? " stroke-dasharray=\"3,2\"" : "");
		for (k = 0; k <= 32; k++)
		{
			qaws_vec2 q;
			qaws_scalar t;
			qaws_skeleton_evaluate_arc(&q, &t, sk, i, (qaws_scalar)(k / 32.0));
			fprintf(s->f, "%.2f,%.2f ", vx(v, q.x), vy(v, q.y));
		}
		fprintf(s->f, "\"/>\n");
	}
	csk_path_curves(s, v, p, sh->nloops, "#0969da", 1.6, "none");
	sprintf(info, "%u arcs, top at t = %.3f", qaws_skeleton_get_arc_count(sk), tmax);
	svg_text(s, v->x0 + v->w / 2, v->y0 + v->h + 13, 11, "#57606a", "middle", info);
	qaws_skeleton_destroy(sk);
	for (i = 0; i < off; i++) qaws_curve_destroy(cs[i]);
}

static void demo_curved_skeleton(void)
{
	static csk_shape sh[6];
	unsigned int n, i;
	double const pi = PI;
	svg s;
	double pw = 300, ph = 250;
	clock_t t0 = clock();
	memset(sh, 0, sizeof(sh));
	/* a rounded slot with a bump */
	n = 0; sh[0].name = "square with a round bump";
	csk_line(&sh[0], &n, 0, 0, 2, 0); csk_line(&sh[0], &n, 2, 0, 2, 2); csk_line(&sh[0], &n, 2, 2, 1.6, 2);
	csk_arc(&sh[0], &n, 1, 2, 0.6, 0, pi); csk_line(&sh[0], &n, 0.4, 2, 0, 2); csk_line(&sh[0], &n, 0, 2, 0, 0);
	sh[0].count[0] = n; sh[0].nloops = 1; sh[0].box[0] = -0.1; sh[0].box[1] = 2.1; sh[0].box[2] = -0.1; sh[0].box[3] = 2.7;
	/* an L with a rounded inner corner and a round end */
	n = 0; sh[1].name = "L with round ends and a fillet";
	csk_line(&sh[1], &n, 0, 0, 3, 0); csk_arc(&sh[1], &n, 3, 0.5, 0.5, -pi / 2, pi / 2); csk_line(&sh[1], &n, 3, 1, 1.4, 1);
	csk_arc(&sh[1], &n, 1.4, 1.4, 0.4, -pi / 2, -pi); csk_line(&sh[1], &n, 1, 1.4, 1, 2.6); csk_arc(&sh[1], &n, 0.5, 2.6, 0.5, 0, pi);
	csk_line(&sh[1], &n, 0, 2.6, 0, 0);
	sh[1].count[0] = n; sh[1].nloops = 1; sh[1].box[0] = -0.1; sh[1].box[1] = 3.6; sh[1].box[2] = -0.1; sh[1].box[3] = 3.2;
	/* a disc with an off-centre hole */
	n = 0; sh[2].name = "disc with an off-centre hole";
	csk_arc(&sh[2], &n, 0, 0, 1, 0, 2 * pi); sh[2].count[0] = 1;
	csk_arc(&sh[2], &n, 0.3, 0.1, 0.35, 2 * pi, 0); sh[2].count[1] = 1; sh[2].nloops = 2;
	sh[2].box[0] = -1.05; sh[2].box[1] = 1.05; sh[2].box[2] = -1.05; sh[2].box[3] = 1.05;
	/* a flower: a hexagon whose sides bulge out */
	n = 0; sh[3].name = "hexagon with bulging sides";
	for (i = 0; i < 6; i++)
	{
		double a0 = 2 * pi * i / 6, a1 = 2 * pi * (i + 1) / 6, ax = cos(a0), ay = sin(a0), bx = cos(a1), by = sin(a1);
		double mx = 0.5 * (ax + bx), my = 0.5 * (ay + by), r = 0.62, l = hypot(bx - ax, by - ay), h = sqrt(r * r - 0.25 * l * l);
		double nx = mx / hypot(mx, my), ny = my / hypot(mx, my), cx = mx - h * nx, cy = my - h * ny;
		csk_arc(&sh[3], &n, cx, cy, r, atan2(ay - cy, ax - cx), atan2(by - cy, bx - cx) + (atan2(by - cy, bx - cx) < atan2(ay - cy, ax - cx) ? 2 * pi : 0));
	}
	sh[3].count[0] = n; sh[3].nloops = 1; sh[3].box[0] = -1.25; sh[3].box[1] = 1.25; sh[3].box[2] = -1.25; sh[3].box[3] = 1.25;
	/* a comb of round-tipped teeth */
	n = 0; sh[4].name = "comb with round teeth";
	csk_line(&sh[4], &n, 0, 0, 5, 0); csk_line(&sh[4], &n, 5, 0, 5, 2.4); csk_arc(&sh[4], &n, 4.5, 2.4, 0.5, 0, pi);
	csk_line(&sh[4], &n, 4, 2.4, 4, 1); csk_arc(&sh[4], &n, 3.5, 1, 0.5, 0, -pi); csk_line(&sh[4], &n, 3, 1, 3, 2.4);
	csk_arc(&sh[4], &n, 2.5, 2.4, 0.5, 0, pi); csk_line(&sh[4], &n, 2, 2.4, 2, 1); csk_arc(&sh[4], &n, 1.5, 1, 0.5, 0, -pi);
	csk_line(&sh[4], &n, 1, 1, 1, 2.4); csk_arc(&sh[4], &n, 0.5, 2.4, 0.5, 0, pi); csk_line(&sh[4], &n, 0, 2.4, 0, 0);
	sh[4].count[0] = n; sh[4].nloops = 1; sh[4].box[0] = -0.1; sh[4].box[1] = 5.1; sh[4].box[2] = -0.1; sh[4].box[3] = 3.0;
	/* a stadium with a square hole */
	n = 0; sh[5].name = "stadium with a square hole";
	csk_line(&sh[5], &n, -1, -1, 1, -1); csk_arc(&sh[5], &n, 1, 0, 1, -pi / 2, pi / 2);
	csk_line(&sh[5], &n, 1, 1, -1, 1); csk_arc(&sh[5], &n, -1, 0, 1, pi / 2, 3 * pi / 2); sh[5].count[0] = n;
	csk_line(&sh[5], &n, -0.3, -0.3, -0.3, 0.3); csk_line(&sh[5], &n, -0.3, 0.3, 0.3, 0.3); csk_line(&sh[5], &n, 0.3, 0.3, 0.3, -0.3);
	csk_line(&sh[5], &n, 0.3, -0.3, -0.3, -0.3); sh[5].count[1] = 4; sh[5].nloops = 2;
	sh[5].box[0] = -2.05; sh[5].box[1] = 2.05; sh[5].box[2] = -1.05; sh[5].box[3] = 1.05;
	if (!svg_open(&s, "showcase/b2d13_curved_skeleton.svg", (int)(40 + 4 * (pw + 12)), (int)(110 + 3 * (ph + 44)),
		"Curved skeletons: lines and circular arcs, mitered and medial axis",
		"qaws_skeleton: the wavefront at eight times coloured by time; skeleton black (thick: the pure part, clear of the boundary), smooth joints grey"))
		return;
	for (i = 0; i < 6; i++)
	{
		unsigned int c;
		for (c = 0; c < 2; c++)
		{
			viewport v;
			double bx = sh[i].box[1] - sh[i].box[0], by = sh[i].box[3] - sh[i].box[2], sc;
			unsigned int col = 2 * (i % 2) + c, row = i / 2;
			char label[160];
			sc = (pw - 20) / bx < (ph - 30) / by ? (pw - 20) / bx : (ph - 30) / by;
			v.w = bx * sc; v.h = by * sc;
			v.x0 = 24 + col * (pw + 12) + (pw - v.w) / 2;
			v.y0 = 104 + row * (ph + 44) + 10 + (ph - 30 - v.h) / 2;
			v.xmin = sh[i].box[0]; v.xmax = sh[i].box[1]; v.ymin = sh[i].box[2]; v.ymax = sh[i].box[3];
			sprintf(label, "%s, %s", sh[i].name, c ? "medial axis (sigma = pi)" : "mitered (sigma = 2 pi)");
			svg_text(&s, 24 + col * (pw + 12) + pw / 2, 104 + row * (ph + 44), 12, "#24292f", "middle", label);
			csk_draw(&s, &v, &sh[i], c ? pi : 0, NULL, NULL);
		}
	}
	svg_close(&s);
	printf("figure 13: curved skeletons, %.1f s\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
}
