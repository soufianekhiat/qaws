/* Figure 14: skeletons of free-form regions (qaws_skeleton on Beziers and
   NURBS through G1 arc splines). Each region mitered (sigma = 2 pi) and as
   its medial axis (sigma = pi); the wavefront at eight times, the skeleton
   (thick: the pure part), smooth joints' paths in grey. */

static qaws_curve* ffs_cubic(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3)
{
	qaws_scalar cp[8];
	qaws_bezier_desc d;
	qaws_curve* c = NULL;
	cp[0] = (qaws_scalar)x0; cp[1] = (qaws_scalar)y0; cp[2] = (qaws_scalar)x1; cp[3] = (qaws_scalar)y1;
	cp[4] = (qaws_scalar)x2; cp[5] = (qaws_scalar)y2; cp[6] = (qaws_scalar)x3; cp[7] = (qaws_scalar)y3;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D; d.degree = 3; d.control_points = cp; d.control_point_count = 4;
	qaws_curve_create_bezier(&c, &d);
	return c;
}

/* a closed Catmull-Rom loop through the points, as cubic Beziers */
static unsigned int ffs_loop(double const* xy, unsigned int n, qaws_curve** out)
{
	unsigned int i;
	for (i = 0; i < n; i++)
	{
		double const* a = xy + 2 * ((i + n - 1) % n);
		double const* b = xy + 2 * i;
		double const* c = xy + 2 * ((i + 1) % n);
		double const* d = xy + 2 * ((i + 2) % n);
		out[i] = ffs_cubic(b[0], b[1], b[0] + (c[0] - a[0]) / 6, b[1] + (c[1] - a[1]) / 6,
			c[0] - (d[0] - b[0]) / 6, c[1] - (d[1] - b[1]) / 6, c[0], c[1]);
	}
	return n;
}

typedef struct ffs_shape
{
	char const* name;
	qaws_curve* cs[32];
	unsigned int count[4], nloops;
	double box[4];
} ffs_shape;

static void ffs_draw(svg* s, viewport const* v, ffs_shape const* sh, double sigma)
{
	qaws_curve const* cv[32];
	qaws_path_2d p[4];
	qaws_skeleton_desc d;
	qaws_skeleton* sk = NULL;
	unsigned int l, i, k, off = 0;
	double tmax;
	char info[128];
	for (l = 0; l < sh->nloops; l++)
	{
		for (i = 0; i < sh->count[l]; i++) cv[off + i] = sh->cs[off + i];
		p[l].curves = cv + off; p[l].curve_count = sh->count[l]; p[l].closed = 1;
		off += sh->count[l];
	}
	csk_path_curves(s, v, p, sh->nloops, "none", 0, "#f6f8fa");
	memset(&d, 0, sizeof(d));
	d.paths = p; d.path_count = sh->nloops; d.fill_rule = QAWS_FILL_NON_ZERO; d.splice_angle = (qaws_scalar)sigma;
	if (qaws_skeleton_create(&d, &sk) != QAWS_STATUS_OK)
	{
		svg_text(s, v->x0 + v->w / 2, v->y0 + v->h / 2, 12, "#cf222e", "middle", "failed");
		return;
	}
	tmax = (double)qaws_skeleton_get_max_time(sk);
	for (k = 1; k <= 8; k++)
	{
		qaws_clip_result* r = NULL;
		char col[32];
		heat((double)k / 9, col);
		if (qaws_skeleton_offset(sk, (qaws_scalar)(tmax * k / 9), &r) == QAWS_STATUS_OK)
		{
			unsigned int n = qaws_clip_result_get_path_count(r), pi;
			qaws_path_2d* q = (qaws_path_2d*)malloc(sizeof(qaws_path_2d) * (n + 1));
			for (pi = 0; pi < n; pi++) qaws_clip_result_get_path(r, pi, &q[pi]);
			csk_path_curves(s, v, q, n, col, 1.0, "none");
			free(q);
		}
		qaws_clip_result_destroy(r);
	}
	for (i = 0; i < qaws_skeleton_get_arc_count(sk); i++)
	{
		qaws_skeleton_arc a;
		qaws_skeleton_get_arc(sk, i, &a);
		if (a.smooth)
			continue;
		fprintf(s->f, "<polyline fill=\"none\" stroke=\"#1b1f24\" stroke-width=\"%.1f\" points=\"", a.boundary ? 1.0 : 2.0);
		for (k = 0; k <= 32; k++)
		{
			qaws_vec2 q;
			qaws_scalar t;
			qaws_skeleton_evaluate_arc(sk, i, (qaws_scalar)(k / 32.0), &q, &t);
			fprintf(s->f, "%.2f,%.2f ", vx(v, q.x), vy(v, q.y));
		}
		fprintf(s->f, "\"/>\n");
	}
	csk_path_curves(s, v, p, sh->nloops, "#0969da", 1.6, "none");
	sprintf(info, "%u faces (arc spline pieces), top at t = %.3f", qaws_skeleton_get_face_count(sk), tmax);
	svg_text(s, v->x0 + v->w / 2, v->y0 + v->h + 13, 11, "#57606a", "middle", info);
	qaws_skeleton_destroy(sk);
}

static void demo_freeform_skeleton(void)
{
	static ffs_shape sh[5];
	unsigned int i, n;
	double const pi = PI;
	svg s;
	double pw = 300, ph = 250;
	clock_t t0 = clock();
	memset(sh, 0, sizeof(sh));
	sh[0].name = "ellipse (rational NURBS)";
	qaws_curve_create_ellipse_2d(0, 0, 1.5, 0.8, 0.35, &sh[0].cs[0]);
	sh[0].count[0] = 1; sh[0].nloops = 1;
	sh[0].box[0] = -1.55; sh[0].box[1] = 1.55; sh[0].box[2] = -1.1; sh[0].box[3] = 1.1;
	sh[1].name = "blob of four cubics";
	sh[1].cs[0] = ffs_cubic(0, 0, 1, -0.9, 2.4, -0.5, 2.2, 0.6);
	sh[1].cs[1] = ffs_cubic(2.2, 0.6, 2, 1.6, 1.4, 1.1, 1, 1);
	sh[1].cs[2] = ffs_cubic(1, 1, 0.6, 0.9, 0.2, 1.8, -0.5, 1.2);
	sh[1].cs[3] = ffs_cubic(-0.5, 1.2, -1.2, 0.6, -0.6, 0.2, 0, 0);
	sh[1].count[0] = 4; sh[1].nloops = 1;
	sh[1].box[0] = -1.0; sh[1].box[1] = 2.4; sh[1].box[2] = -0.75; sh[1].box[3] = 1.6;
	{
		double xy[2 * 20];
		for (i = 0; i < 20; i++)
		{
			double a = 2 * pi * i / 20, r = 1 + 0.32 * cos(5 * a);
			xy[2 * i] = r * cos(a); xy[2 * i + 1] = r * sin(a);
		}
		sh[2].name = "flower: closed Catmull-Rom of cubics";
		n = ffs_loop(xy, 20, sh[2].cs);
		sh[2].count[0] = n; sh[2].nloops = 1;
		sh[2].box[0] = -1.35; sh[2].box[1] = 1.35; sh[2].box[2] = -1.35; sh[2].box[3] = 1.35;
	}
	{
		qaws_arc_segment g;
		qaws_arc_desc d;
		sh[3].name = "ellipse with two round holes";
		qaws_curve_create_ellipse_2d(0, 0, 1.6, 1.0, 0, &sh[3].cs[0]);
		sh[3].count[0] = 1;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.segments = &g; d.segment_count = 1;
		memset(&g, 0, sizeof(g));
		g.center[0] = (qaws_scalar)-0.6; g.center[1] = (qaws_scalar)0.1; g.radius = (qaws_scalar)0.3; g.angle_start = (qaws_scalar)(2 * pi); g.angle_end = 0;
		qaws_curve_create_arc(&sh[3].cs[1], &d);
		g.center[0] = (qaws_scalar)0.7; g.center[1] = (qaws_scalar)-0.2; g.radius = (qaws_scalar)0.22;
		qaws_curve_create_arc(&sh[3].cs[2], &d);
		sh[3].count[1] = 1; sh[3].count[2] = 1; sh[3].nloops = 3;
		sh[3].box[0] = -1.65; sh[3].box[1] = 1.65; sh[3].box[2] = -1.05; sh[3].box[3] = 1.05;
	}
	{
		/* a square with a wavy cubic top */
		qaws_scalar cp[4];
		qaws_bezier_desc bd;
		sh[4].name = "square with a wavy cubic edge";
		memset(&bd, 0, sizeof(bd));
		bd.dimension = QAWS_DIMENSION_2D; bd.degree = 1; bd.control_points = cp; bd.control_point_count = 2;
		cp[0] = 0; cp[1] = 0; cp[2] = 2; cp[3] = 0; qaws_curve_create_bezier(&sh[4].cs[0], &bd);
		cp[0] = 2; cp[1] = 0; cp[2] = 2; cp[3] = 1.6; qaws_curve_create_bezier(&sh[4].cs[1], &bd);
		sh[4].cs[2] = ffs_cubic(2, 1.6, 1.4, 2.6, 0.6, 0.6, 0, 1.6);
		cp[0] = 0; cp[1] = 1.6; cp[2] = 0; cp[3] = 0; qaws_curve_create_bezier(&sh[4].cs[3], &bd);
		sh[4].count[0] = 4; sh[4].nloops = 1;
		sh[4].box[0] = -0.05; sh[4].box[1] = 2.05; sh[4].box[2] = -0.05; sh[4].box[3] = 2.15;
	}
	if (!svg_open(&s, "showcase/b2d14_freeform_skeleton.svg", (int)(40 + 4 * (pw + 12)), (int)(110 + 3 * (ph + 44)),
		"Skeletons of free-form regions: Beziers and NURBS through G1 arc splines",
		"qaws_skeleton: each curve as equal-tangent biarcs within 1e-4 of the extent; mitered and medial axis; wavefronts by time"))
		return;
	for (i = 0; i < 5; i++)
	{
		unsigned int c;
		for (c = 0; c < 2; c++)
		{
			viewport v;
			double bx = sh[i].box[1] - sh[i].box[0], by = sh[i].box[3] - sh[i].box[2], sc;
			unsigned int idx = 2 * i + c, col = idx % 4, row = idx / 4;
			char label[160];
			sc = (pw - 20) / bx < (ph - 30) / by ? (pw - 20) / bx : (ph - 30) / by;
			v.w = bx * sc; v.h = by * sc;
			v.x0 = 24 + col * (pw + 12) + (pw - v.w) / 2;
			v.y0 = 104 + row * (ph + 44) + 10 + (ph - 30 - v.h) / 2;
			v.xmin = sh[i].box[0]; v.xmax = sh[i].box[1]; v.ymin = sh[i].box[2]; v.ymax = sh[i].box[3];
			sprintf(label, "%s, %s", sh[i].name, c ? "medial axis" : "mitered");
			svg_text(&s, 24 + col * (pw + 12) + pw / 2, 104 + row * (ph + 44), 12, "#24292f", "middle", label);
			ffs_draw(&s, &v, &sh[i], c ? pi : 0);
		}
	}
	svg_close(&s);
	for (i = 0; i < 5; i++)
	{
		unsigned int k, tot = sh[i].count[0] + sh[i].count[1] + sh[i].count[2] + sh[i].count[3];
		for (k = 0; k < tot; k++) if (sh[i].cs[k]) qaws_curve_destroy(sh[i].cs[k]);
	}
	printf("figure 14: free-form skeletons, %.1f s\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
}
