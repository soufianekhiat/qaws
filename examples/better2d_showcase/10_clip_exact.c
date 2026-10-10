/* Figure 10: exact Boolean operations on curves (qaws_clip_exact) on the
   degenerate cases a tolerance gets wrong: shared edges, touching corners,
   tangencies, three edges through one rational point, nesting, a duplicate
   curve, a sliver an integer engine rounds away, and the one kind of case it
   refuses (a tangency at an irrational point). Results filled with their
   pieces in alternating colours, inputs dashed. */

typedef struct ce_scene
{
	qaws_curve* curves[16];
	qaws_curve const* views[16];
	qaws_path_2d paths[16];
	unsigned int n, ns;
} ce_scene;

static void ce_add(ce_scene* s, qaws_curve* c)
{
	s->curves[s->n] = c;
	s->views[s->n] = c;
	s->paths[s->n].curves = &s->views[s->n];
	s->paths[s->n].curve_count = 1;
	s->paths[s->n].closed = 1;
	s->n++;
}

static qaws_curve* ce_rect(double x0, double y0, double x1, double y1)
{
	double xy[8];
	xy[0] = x0; xy[1] = y0; xy[2] = x1; xy[3] = y0; xy[4] = x1; xy[5] = y1; xy[6] = x0; xy[7] = y1;
	return bo_polyline(xy, 4);
}

static qaws_curve* ce_bezier(unsigned int degree, double const* p)
{
	qaws_scalar cp[12];
	qaws_bezier_desc d;
	qaws_curve* c = NULL;
	unsigned int i;
	for (i = 0; i < 2 * (degree + 1); i++) cp[i] = (qaws_scalar)p[i];
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D; d.degree = degree; d.control_points = cp; d.control_point_count = degree + 1;
	qaws_curve_create_bezier(&d, &c);
	return c;
}

/* a closed path of one curve and the line back to its start */
static void ce_add_closed(ce_scene* s, qaws_curve* c)
{
	qaws_eval_result_2d a, b;
	qaws_range r = qaws_curve_get_parameter_range(c);
	double xy[4];
	qaws_curve* line = NULL;
	qaws_scalar cp[4];
	qaws_bezier_desc d;
	qaws_curve_evaluate_2d(c, r.max_value, QAWS_EVAL_FLAG_POSITION, &b);
	qaws_curve_evaluate_2d(c, r.min_value, QAWS_EVAL_FLAG_POSITION, &a);
	xy[0] = b.position.x; xy[1] = b.position.y; xy[2] = a.position.x; xy[3] = a.position.y;
	cp[0] = (qaws_scalar)xy[0]; cp[1] = (qaws_scalar)xy[1]; cp[2] = (qaws_scalar)xy[2]; cp[3] = (qaws_scalar)xy[3];
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D; d.degree = 1; d.control_points = cp; d.control_point_count = 2;
	qaws_curve_create_bezier(&d, &line);
	s->curves[s->n] = c;
	s->views[s->n] = c;
	s->curves[s->n + 1] = line;
	s->views[s->n + 1] = line;
	s->paths[s->n].curves = &s->views[s->n];
	s->paths[s->n].curve_count = 2;
	s->paths[s->n].closed = 1;
	s->n += 2;
}

static void ce_free(ce_scene* s)
{
	unsigned int i;
	for (i = 0; i < s->n; i++) qaws_curve_destroy(s->curves[i]);
}

static void demo_clip_exact(void)
{
	int const cols = 3, rows = 3, pw = 400, ph = 340, gap = 16, top = 84;
	int const W = cols * pw + (cols + 1) * gap, H = top + rows * (ph + gap) + 8;
	static char const* status_names[] = { "ok" };
	svg s;
	int k;
	double t0;
	if (!svg_open(&s, "showcase/b2d10_clip_exact.svg", W, H, "Exact Boolean operations on curves",
		"qaws_clip_exact: certified hits (tangencies with multiplicity, overlaps with exact ends), vertices proven equal, order round vertices from exact sub-curves, faces by rays; no tolerance anywhere"))
		return;
	(void)status_names;
	t0 = (double)clock() / CLOCKS_PER_SEC;
	for (k = 0; k < 9; k++)
	{
		ce_scene sc;
		qaws_clip_exact_desc d;
		qaws_clip_exact_stats st;
		qaws_clip_result* r = NULL;
		qaws_status status;
		viewport v;
		char const* title = "";
		char label[256];
		unsigned int i, first_clip = 0;
		double lo[2] = { 1e300, 1e300 }, hi[2] = { -1e300, -1e300 }, mx, my;
		memset(&sc, 0, sizeof(sc));
		memset(&d, 0, sizeof(d));
		d.fill_rule = QAWS_FILL_NON_ZERO;
		switch (k)
		{
		case 0:
			ce_add(&sc, ce_rect(0, 0, 1, 1)); ce_add(&sc, ce_rect(1, 0, 2, 1)); ce_add(&sc, ce_rect(0.5, 1, 1.5, 1.75));
			first_clip = 3; d.clip_type = QAWS_CLIP_UNION;
			title = "union of squares sharing edges (overlaps merged)";
			break;
		case 1:
			ce_add(&sc, ce_rect(0, 0, 1, 1)); ce_add(&sc, ce_rect(1, 1, 2, 2));
			first_clip = 1; d.clip_type = QAWS_CLIP_XOR;
			title = "xor of squares touching at a corner";
			break;
		case 2:
		{
			static double const par[6] = { 0, 1, 1, -1, 2, 1 };
			ce_add(&sc, ce_rect(0, 0, 2, 2));
			ce_add_closed(&sc, ce_bezier(2, par));
			first_clip = 1; d.clip_type = QAWS_CLIP_DIFFERENCE;
			title = "square - parabola touching its side (tangency)";
			break;
		}
		case 3:
		{
			static double const t1[6] = { 0, 0, 1, 1, -1, 1 }, t2[6] = { 1, 0, 0, 0.5, 1, 1 }, t3[6] = { 0, 1, 0.5, 0, 0, 0 };
			ce_add(&sc, bo_polyline(t1, 3)); ce_add(&sc, bo_polyline(t2, 3)); ce_add(&sc, bo_polyline(t3, 3));
			first_clip = 3; d.clip_type = QAWS_CLIP_UNION; d.fill_rule = QAWS_FILL_EVEN_ODD;
			title = "three edges through (1/3, 1/3), even-odd";
			break;
		}
		case 4:
		{
			double hole[8] = { 1, 1, 1, 5, 5, 5, 5, 1 };
			ce_add(&sc, ce_rect(0, 0, 6, 6)); ce_add(&sc, bo_polyline(hole, 4)); ce_add(&sc, ce_rect(2, 2, 4, 4));
			ce_add(&sc, bo_blob(3, 3, 1.6, 0.25, 5, 0.3));
			first_clip = 3; d.clip_type = QAWS_CLIP_XOR;
			title = "a square with a hole and an island, xor a blob";
			break;
		}
		case 5:
		{
			double xy[10];
			unsigned int q;
			for (q = 0; q < 5; q++)
			{
				double a = PI / 2 + q * 4 * PI / 5;
				xy[2 * q] = floor(cos(a) * 1024) / 256;
				xy[2 * q + 1] = floor(sin(a) * 1024) / 256;
			}
			ce_add(&sc, bo_polyline(xy, 5));
			first_clip = 1; d.clip_type = QAWS_CLIP_UNION; d.fill_rule = QAWS_FILL_EVEN_ODD;
			title = "a pentagram, even-odd";
			break;
		}
		case 6:
		{
			unsigned int q;
			for (q = 0; q < 4; q++)
			{
				qaws_curve* e = NULL;
				qaws_curve_create_ellipse_2d((qaws_scalar)(0.75 * q), (qaws_scalar)(0.25 * q), 1.5f, 1.0f, (qaws_scalar)(0.3 * q), &e);
				ce_add(&sc, e);
			}
			{
				qaws_curve* e = NULL;
				qaws_curve_create_ellipse_2d(0.75f, 0.25f, 1.5f, 1.0f, 0.3f, &e);
				ce_add(&sc, e);
			}
			first_clip = 5; d.clip_type = QAWS_CLIP_UNION;
			title = "four ellipses and a duplicate of one, union";
			break;
		}
		case 7:
		{
			double sj[6] = { 14, 15, 16, 12, 10, 12 }, cl[12] = { 15, 14, 11, 14, 13, 16, 17, 10, 10, 17, 18, 13 };
			ce_add(&sc, bo_polyline(sj, 3)); ce_add(&sc, bo_polyline(cl, 6));
			first_clip = 1; d.clip_type = QAWS_CLIP_INTERSECTION; d.fill_rule = QAWS_FILL_EVEN_ODD;
			title = "Clipper2 test 22: the sliver integers round away";
			break;
		}
		default:
		{
			static double const quartic[10] = { -1, 0.75, -0.5, -2.25, 0, 4.75, 0.5, -2.25, 1, 0.75 };
			ce_add_closed(&sc, ce_bezier(4, quartic));
			ce_add(&sc, ce_rect(-1.5, -1, 1.5, 0));
			first_clip = 1; d.clip_type = QAWS_CLIP_UNION;
			title = "touching at x = +-1/sqrt(2): refused, not guessed";
			break;
		}
		}
		/* subjects: the paths before first_clip (paths are at the curve index of their first curve) */
		{
			qaws_path_2d subj[16], clip[16];
			unsigned int nsj = 0, ncl = 0, ci = 0, pi;
			for (pi = 0; pi < sc.n; pi += sc.paths[pi].curve_count ? sc.paths[pi].curve_count : 1, ci++)
			{
				if (ci < first_clip) subj[nsj++] = sc.paths[pi];
				else clip[ncl++] = sc.paths[pi];
			}
			d.subjects = subj; d.subject_count = nsj;
			d.clips = clip; d.clip_count = ncl;
			status = qaws_clip_exact_execute(&d, &r, &st);
		}
		/* the view */
		for (i = 0; i < sc.n; i++)
		{
			qaws_range rr = qaws_curve_get_parameter_range(sc.curves[i]);
			unsigned int j;
			for (j = 0; j <= 64; j++)
			{
				qaws_eval_result_2d e;
				qaws_curve_evaluate_2d(sc.curves[i], (qaws_scalar)(rr.min_value + (rr.max_value - rr.min_value) * j / 64), QAWS_EVAL_FLAG_POSITION, &e);
				if (e.position.x < lo[0]) lo[0] = e.position.x;
				if (e.position.y < lo[1]) lo[1] = e.position.y;
				if (e.position.x > hi[0]) hi[0] = e.position.x;
				if (e.position.y > hi[1]) hi[1] = e.position.y;
			}
		}
		mx = 0.5 * (lo[0] + hi[0]);
		my = 0.5 * (lo[1] + hi[1]);
		{
			double half = 0.56 * fmax((hi[0] - lo[0]) / (pw - 20.0) * (ph - 60.0), hi[1] - lo[1]);
			double halfx = half * (pw - 20.0) / (ph - 60.0);
			v.x0 = gap + (k % cols) * (pw + gap) + 10; v.y0 = top + (k / cols) * (ph + gap) + 28; v.w = pw - 20; v.h = ph - 60;
			v.xmin = mx - halfx; v.xmax = mx + halfx; v.ymin = my - half; v.ymax = my + half;
		}
		{
			viewport pv = v;
			pv.x0 = gap + (k % cols) * (pw + gap); pv.y0 = top + (k / cols) * (ph + gap); pv.w = pw; pv.h = ph;
			svg_panel(&s, &pv, title);
		}
		if (status == QAWS_STATUS_OK)
			bo_fill(&s, &v, r, k % 2 ? "#54aeff" : "#4ac26b");
		for (i = 0; i < sc.n; i++)
			bo_curve(&s, &v, sc.curves[i], "#57606a", 1.1, 1);
		if (status == QAWS_STATUS_OK)
		{
			bo_boundaries(&s, &v, r);
			sprintf(label, "%u paths; %u vertices, %u edges, %u faces; %u overlaps, %u tangencies, %u identity proofs",
				qaws_clip_result_get_path_count(r), st.vertex_count, st.edge_count, st.face_count, st.overlap_count, st.tangency_count, st.identity_tests);
		}
		else
			sprintf(label, "CERTIFICATION_FAILED: tangent at an irrational vertex (%u tangencies)", st.tangency_count);
		svg_text(&s, gap + (k % cols) * (pw + gap) + 10, top + (k / cols) * (ph + gap) + ph - 10, 11,
			status == QAWS_STATUS_OK ? "#57606a" : "#cf222e", "start", label);
		qaws_clip_result_destroy(r);
		ce_free(&sc);
	}
	svg_close(&s);
	printf("    clip exact: 9 panels in %.3f s\n", (double)clock() / CLOCKS_PER_SEC - t0);
}
