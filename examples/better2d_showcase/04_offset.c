/* Figure 4: offsetting paths of curves (qaws_offset). Polygons with each
   join type, open polylines with each end type, and curved regions offset
   at a series of distances in both directions, cusps and overlaps removed by
   the final union. */

static void of_draw_result(svg* s, viewport const* v, qaws_clip_result const* r, char const* color, double width)
{
	unsigned int i, k;
	for (i = 0; i < qaws_clip_result_get_path_count(r); i++)
	{
		qaws_path_2d p;
		qaws_clip_result_get_path(r, i, &p);
		for (k = 0; k < p.curve_count; k++)
			bo_curve(s, v, p.curves[k], color, width, 0);
	}
}

static void of_series(svg* s, viewport const* v, qaws_path_2d const* p, unsigned int n, double const* deltas, unsigned int nd,
	qaws_join_type jt, qaws_end_type et, unsigned int* total_paths)
{
	unsigned int i;
	for (i = 0; i < nd; i++)
	{
		qaws_clip_result* r = NULL;
		char col[16];
		heat(0.15 + 0.7 * i / (nd > 1 ? nd - 1 : 1), col);
		if (qaws_offset_paths(p, n, (qaws_scalar)deltas[i], jt, et, &r) == QAWS_STATUS_OK)
		{
			if (i == 0 || deltas[i] * deltas[0] < 0)
				bo_fill(s, v, r, "#c8e1ff");
			of_draw_result(s, v, r, col, 1.6);
			*total_paths += qaws_clip_result_get_path_count(r);
		}
		qaws_clip_result_destroy(r);
	}
}

static void demo_offset(void)
{
	int const cols = 4, rows = 3, pw = 300, ph = 260, gap = 16, top = 84;
	int const W = cols * pw + (cols + 1) * gap, H = top + rows * (ph + gap) + 8;
	static char const* jn[4] = { "square", "bevel", "round", "miter" };
	static char const* en[5] = { "polygon", "joined", "butt", "square", "round" };
	svg s;
	int k;
	double t0;
	unsigned int paths = 0;
	char label[200];
	qaws_curve *zig = NULL, *open = NULL, *blob, *ell = NULL, *spiral = NULL;
	if (!svg_open(&s, "showcase/b2d4_offset.svg", W, H, "Offsetting paths of curves",
		"qaws_offset: lines translated and arcs made concentric exactly, other curves fitted to the true offset; joins and caps as Clipper2; one positive union cleans every loop"))
		return;
	{
		qaws_scalar z[16] = { 0, 0, 4, 0, 4, 3, 2, 1, 1, 4, 3, (qaws_scalar)5.5, 0, 5, (qaws_scalar)-1.5, 2 };
		qaws_curve_create_polyline_2d(z, 8, 1, &zig);
	}
	{
		qaws_scalar o[12] = { 0, 0, 2, 2, 4, 0, 6, 2, 8, 0, 9, (qaws_scalar)1.5 };
		qaws_curve_create_polyline_2d(o, 6, 0, &open);
	}
	blob = bo_blob(0, 0, 2.0, 0.35, 3, 0.2);
	qaws_curve_create_ellipse_2d(0, 0, 3, (qaws_scalar)1.4, (qaws_scalar)0.3, &ell);
	{
		qaws_clothoid_desc cd;
		memset(&cd, 0, sizeof(cd));
		cd.start_curvature = 0; cd.end_curvature = (qaws_scalar)1.6; cd.length = 9;
		qaws_curve_create_clothoid(&cd, &spiral);
	}
	t0 = (double)clock() / CLOCKS_PER_SEC;
	for (k = 0; k < 12; k++)
	{
		viewport v;
		qaws_curve const* c;
		qaws_path_2d p;
		v.x0 = gap + (k % cols) * (pw + gap); v.y0 = top + (k / cols) * (ph + gap); v.w = pw; v.h = ph;
		if (k < 4)
		{
			/* a concave polygon grown with each join type */
			static double const d[3] = { 0.25, 0.5, 0.8 };
			c = zig;
			p.curves = &c; p.curve_count = 1; p.closed = 1;
			v.xmin = -3.0; v.xmax = 6.2; v.ymin = -1.4; v.ymax = v.ymin + (v.xmax - v.xmin) * ph / pw;
			sprintf(label, "polygon grown, %s joins", jn[k]);
			svg_panel(&s, &v, label);
			of_series(&s, &v, &p, 1, d, 3, (qaws_join_type)k, QAWS_END_POLYGON, &paths);
			bo_curve(&s, &v, zig, "#24292f", 1.2, 1);
		}
		else if (k < 8)
		{
			/* an open polyline with each end type (round joins) */
			static double const d[2] = { 0.35, 0.7 };
			qaws_end_type et = (qaws_end_type)(k - 4 + 1);
			c = open;
			p.curves = &c; p.curve_count = 1; p.closed = 0;
			v.xmin = -1.3; v.xmax = 10.3; v.ymin = -2.6; v.ymax = v.ymin + (v.xmax - v.xmin) * ph / pw;
			sprintf(label, "open polyline, %s ends", en[et]);
			svg_panel(&s, &v, label);
			of_series(&s, &v, &p, 1, d, 2, k == 7 ? QAWS_JOIN_MITER : QAWS_JOIN_ROUND, et, &paths);
			bo_curve(&s, &v, open, "#24292f", 1.2, 1);
		}
		else
		{
			static double const din[6] = { -0.9, -0.6, -0.3, 0.3, 0.6, 0.9 };
			static double const dcl[3] = { 0.3, 0.6, 0.9 };
			if (k == 8 || k == 9)
			{
				c = k == 8 ? blob : ell;
				p.curves = &c; p.curve_count = 1; p.closed = 1;
				v.xmin = -4.2; v.xmax = 4.2; v.ymin = -3.6; v.ymax = v.ymin + (v.xmax - v.xmin) * ph / pw;
				sprintf(label, "%s, offsets -0.9 .. +0.9", k == 8 ? "cubic B-spline blob" : "NURBS ellipse");
				svg_panel(&s, &v, label);
				of_series(&s, &v, &p, 1, din, 6, QAWS_JOIN_ROUND, QAWS_END_POLYGON, &paths);
				bo_curve(&s, &v, c, "#24292f", 1.4, 1);
			}
			else if (k == 10)
			{
				c = spiral;
				p.curves = &c; p.curve_count = 1; p.closed = 0;
				v.xmin = -1.5; v.xmax = 7.0; v.ymin = -1.6; v.ymax = v.ymin + (v.xmax - v.xmin) * ph / pw;
				sprintf(label, "clothoid, round ends, 0.3 .. 0.9");
				svg_panel(&s, &v, label);
				of_series(&s, &v, &p, 1, dcl, 3, QAWS_JOIN_ROUND, QAWS_END_ROUND, &paths);
				bo_curve(&s, &v, spiral, "#24292f", 1.4, 1);
			}
			else
			{
				/* a ring: outer blob and inner ellipse as a hole, grown */
				qaws_curve* hole = NULL;
				qaws_curve* outer = bo_blob(0, 0, 3.4, 0.2, 4, 0.0);
				qaws_curve* inner_ccw = NULL;
				qaws_curve const* cs[2];
				qaws_path_2d pp[2];
				static double const dr[3] = { 0.25, 0.5, -0.35 };
				qaws_curve_create_ellipse_2d(0, 0, (qaws_scalar)1.6, (qaws_scalar)1.0, (qaws_scalar)0.4, &inner_ccw);
				qaws_curve_reverse(inner_ccw, &hole);
				cs[0] = outer; cs[1] = hole;
				pp[0].curves = &cs[0]; pp[0].curve_count = 1; pp[0].closed = 1;
				pp[1].curves = &cs[1]; pp[1].curve_count = 1; pp[1].closed = 1;
				v.xmin = -4.6; v.xmax = 4.6; v.ymin = -4.0; v.ymax = v.ymin + (v.xmax - v.xmin) * ph / pw;
				sprintf(label, "region with a hole, +0.25 +0.5 -0.35");
				svg_panel(&s, &v, label);
				of_series(&s, &v, pp, 2, dr, 3, QAWS_JOIN_ROUND, QAWS_END_POLYGON, &paths);
				bo_curve(&s, &v, outer, "#24292f", 1.2, 1);
				bo_curve(&s, &v, hole, "#24292f", 1.2, 1);
				qaws_curve_destroy(outer); qaws_curve_destroy(hole); qaws_curve_destroy(inner_ccw);
			}
		}
	}
	svg_close(&s);
	printf("    offset: 12 panels, %u result paths in %.3f s\n", paths, (double)clock() / CLOCKS_PER_SEC - t0);
	qaws_curve_destroy(zig); qaws_curve_destroy(open); qaws_curve_destroy(blob); qaws_curve_destroy(ell); qaws_curve_destroy(spiral);
}
