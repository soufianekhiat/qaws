/* Figure 2: Boolean operations on regions of curves (qaws_clip). Every
   result path filled (holes cut out by their reversed orientation), its
   boundary pieces drawn in alternating colours to show where the new curves
   start and end, and the inputs dashed. */

#define BO_SAMPLES 160

static void bo_path_points(qaws_path_2d const* p, viewport const* v, double** xy, unsigned int* n, unsigned int* cap)
{
	unsigned int i, j;
	for (i = 0; i < p->curve_count; i++)
	{
		qaws_range r = qaws_curve_get_parameter_range(p->curves[i]);
		for (j = 0; j < BO_SAMPLES; j++)
		{
			qaws_eval_result_2d e;
			if (*n == *cap)
			{
				*cap = *cap ? 2 * *cap : 1024;
				*xy = (double*)realloc(*xy, sizeof(double) * 2 * *cap);
			}
			qaws_curve_evaluate_2d(p->curves[i], (qaws_scalar)(r.min_value + (r.max_value - r.min_value) * j / BO_SAMPLES),
				QAWS_EVAL_FLAG_POSITION, &e);
			(*xy)[2 * *n] = vx(v, e.position.x);
			(*xy)[2 * *n + 1] = vy(v, e.position.y);
			(*n)++;
		}
	}
}

/* every closed path of the result as one SVG path with the non-zero rule:
   holes run the other way, so they come out empty */
static void bo_fill(svg* s, viewport const* v, qaws_clip_result const* r, char const* fill)
{
	unsigned int i, k, n, cap = 0;
	double* xy = NULL;
	fprintf(s->f, "<path fill=\"%s\" fill-opacity=\"0.55\" fill-rule=\"nonzero\" stroke=\"none\" d=\"", fill);
	for (i = 0; i < qaws_clip_result_get_path_count(r); i++)
	{
		qaws_path_2d p;
		qaws_clip_result_get_path(r, i, &p);
		n = 0;
		bo_path_points(&p, v, &xy, &n, &cap);
		for (k = 0; k < n; k++)
			fprintf(s->f, "%s%.2f,%.2f ", k ? "L" : "M", xy[2 * k], xy[2 * k + 1]);
		fprintf(s->f, "Z ");
	}
	fprintf(s->f, "\"/>\n");
	free(xy);
}

static void bo_curve(svg* s, viewport const* v, qaws_curve const* c, char const* color, double width, int dashed)
{
	double xy[2 * (BO_SAMPLES + 1)];
	qaws_range r = qaws_curve_get_parameter_range(c);
	int j;
	for (j = 0; j <= BO_SAMPLES; j++)
	{
		qaws_eval_result_2d e;
		qaws_curve_evaluate_2d(c, (qaws_scalar)(r.min_value + (r.max_value - r.min_value) * j / BO_SAMPLES), QAWS_EVAL_FLAG_POSITION, &e);
		xy[2 * j] = vx(v, e.position.x);
		xy[2 * j + 1] = vy(v, e.position.y);
	}
	svg_polyline(s, xy, BO_SAMPLES + 1, color, width, 0.95, dashed);
}

static void bo_boundaries(svg* s, viewport const* v, qaws_clip_result const* r)
{
	static char const* cols[2] = { "#0969da", "#bf3989" };
	unsigned int i, k, c = 0;
	for (i = 0; i < qaws_clip_result_get_path_count(r); i++)
	{
		qaws_path_2d p;
		qaws_clip_vertex const* vt = NULL;
		unsigned int nv = 0;
		qaws_clip_result_get_path(r, i, &p);
		for (k = 0; k < p.curve_count; k++)
			bo_curve(s, v, p.curves[k], cols[c++ & 1], 2.4, 0);
		qaws_clip_result_get_vertices(r, 0, i, &vt, &nv);
		for (k = 0; k < nv; k++)
			if (vt[k].operand_b != QAWS_CLIP_NONE_INDEX)
				svg_circle(s, vx(v, vt[k].position.x), vy(v, vt[k].position.y), 3.2, "#1a7f37", "#ffffff");
	}
	for (i = 0; i < qaws_clip_result_get_open_path_count(r); i++)
	{
		qaws_path_2d p;
		qaws_clip_result_get_open_path(r, i, &p);
		for (k = 0; k < p.curve_count; k++)
			bo_curve(s, v, p.curves[k], "#cf222e", 3.2, 0);
	}
}

static qaws_curve* bo_polyline(double const* xy, unsigned int n)
{
	qaws_scalar p[128];
	qaws_curve* c = NULL;
	unsigned int i;
	for (i = 0; i < 2 * n; i++) p[i] = (qaws_scalar)xy[i];
	qaws_curve_create_polyline_2d(p, n, 1, &c);
	return c;
}

/* a closed blob: a periodic-looking cubic B-spline through a ring of points */
static qaws_curve* bo_blob(double cx, double cy, double r, double wobble, unsigned int lobes, double phase)
{
	enum { N = 12 };
	qaws_scalar cp[2 * (N + 3)], kn[N + 3 + 4];
	qaws_bspline_desc d;
	qaws_curve* c = NULL;
	unsigned int i;
	for (i = 0; i < N + 3; i++)
	{
		double a = 2 * PI * (i % N) / N + phase, rr = r * (1 + wobble * cos(lobes * a));
		cp[2 * i] = (qaws_scalar)(cx + rr * cos(a));
		cp[2 * i + 1] = (qaws_scalar)(cy + rr * sin(a));
	}
	for (i = 0; i < N + 3 + 4; i++)
		kn[i] = (qaws_scalar)i;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D; d.degree = 3; d.control_points = cp; d.control_point_count = N + 3;
	d.knots = kn; d.knot_count = N + 3 + 4;
	qaws_curve_create_bspline(&d, &c);
	return c;
}

static void demo_boolean(void)
{
	int const cols = 4, rows = 3, pw = 300, ph = 260, gap = 16, top = 84;
	int const W = cols * pw + (cols + 1) * gap, H = top + rows * (ph + gap) + 8;
	svg s;
	unsigned int panel = 0;
	qaws_curve *sq_a, *sq_b, *circ_a, *circ_b, *blob_a, *blob_b, *star, *ring_o, *ring_h, *ring_i, *line;
	double t0;
	char label[200];
	if (!svg_open(&s, "showcase/b2d2_boolean.svg", W, H, "Boolean operations on regions of curves",
		"qaws_clip: one batched intersection, a planar arrangement, Clipper2's truth tables per face; the result is new curves on the inputs (alternating colours show the pieces), green dots where inputs cross"))
		return;
	{
		double a[8] = { 0, 0, 2, 0, 2, 2, 0, 2 }, b[8] = { 1, 1, 3, 1, 3, 3, 1, 3 };
		sq_a = bo_polyline(a, 4);
		sq_b = bo_polyline(b, 4);
	}
	qaws_curve_create_ellipse_2d(0, 0, (qaws_scalar)1.2, (qaws_scalar)1.2, 0, &circ_a);
	qaws_curve_create_ellipse_2d((qaws_scalar)1.0, (qaws_scalar)0.3, (qaws_scalar)1.4, (qaws_scalar)0.8, (qaws_scalar)0.5, &circ_b);
	blob_a = bo_blob(0, 0, 1.2, 0.25, 3, 0.0);
	blob_b = bo_blob(0.9, 0.2, 1.0, 0.3, 4, 0.3);
	{
		double st[10];
		unsigned int i;
		for (i = 0; i < 5; i++)
		{
			double a = PI / 2 + 2 * PI * (2 * i % 5) / 5;
			st[2 * i] = cos(a);
			st[2 * i + 1] = sin(a);
		}
		star = bo_polyline(st, 5);
	}
	qaws_curve_create_ellipse_2d(0, 0, 2, 2, 0, &ring_o);
	qaws_curve_create_ellipse_2d(0, 0, (qaws_scalar)1.3, (qaws_scalar)1.3, 0, &ring_h);
	{
		double sq[8] = { -0.6, -0.6, 0.6, -0.6, 0.6, 0.6, -0.6, 0.6 };
		ring_i = bo_polyline(sq, 4);
	}
	{
		qaws_scalar lp[8] = { -2.5f, -1.0f, -0.5f, 1.5f, 0.8f, -1.6f, 2.5f, 0.9f };
		qaws_curve_create_polyline_2d(lp, 4, 0, &line);
	}
	t0 = (double)clock() / CLOCKS_PER_SEC;
	/* row 1: squares, every clip type; row 2: the same on curves; row 3: fill rules, nesting, open paths */
	{
		struct { qaws_curve const* a; qaws_curve const* b; qaws_clip_type ct; qaws_fill_rule fr; char const* name; double box[4]; } jobs[12];
		static char const* opn[5] = { "none", "intersection", "union", "difference", "xor" };
		int k;
		for (k = 0; k < 4; k++)
		{
			jobs[k].a = sq_a; jobs[k].b = sq_b; jobs[k].ct = (qaws_clip_type)(k + 1); jobs[k].fr = QAWS_FILL_NON_ZERO;
			jobs[k].box[0] = -0.5; jobs[k].box[1] = 3.5; jobs[k].box[2] = -0.5; jobs[k].box[3] = 3.5;
			jobs[k].name = opn[k + 1];
		}
		for (k = 0; k < 4; k++)
		{
			jobs[4 + k].a = k < 2 ? circ_a : blob_a; jobs[4 + k].b = k < 2 ? circ_b : blob_b;
			jobs[4 + k].ct = (qaws_clip_type)(k == 0 ? QAWS_CLIP_UNION : k == 1 ? QAWS_CLIP_XOR : k == 2 ? QAWS_CLIP_INTERSECTION : QAWS_CLIP_DIFFERENCE);
			jobs[4 + k].fr = QAWS_FILL_NON_ZERO;
			jobs[4 + k].box[0] = -1.8; jobs[4 + k].box[1] = 2.6; jobs[4 + k].box[2] = -1.7; jobs[4 + k].box[3] = 1.9;
			jobs[4 + k].name = opn[jobs[4 + k].ct];
		}
		for (k = 0; k < 12; k++)
		{
			viewport v;
			qaws_clip_result* r = NULL;
			int col = k % cols, row = k / cols;
			qaws_path_2d sp[3], cp;
			qaws_curve const* sc[3];
			qaws_curve const* cc;
			qaws_clip_desc d;
			double area = 0;
			unsigned int i;
			v.x0 = gap + col * (pw + gap); v.y0 = top + row * (ph + gap); v.w = pw; v.h = ph;
			memset(&d, 0, sizeof(d));
			if (k < 8)
			{
				sc[0] = jobs[k].a; cc = jobs[k].b;
				sp[0].curves = &sc[0]; sp[0].curve_count = 1; sp[0].closed = 1;
				cp.curves = &cc; cp.curve_count = 1; cp.closed = 1;
				d.subjects = sp; d.subject_count = 1; d.clips = &cp; d.clip_count = 1;
				d.clip_type = jobs[k].ct; d.fill_rule = jobs[k].fr;
				v.xmin = jobs[k].box[0]; v.xmax = jobs[k].box[1];
				v.ymin = jobs[k].box[2] - (jobs[k].box[3] - jobs[k].box[2]) * (ph * 1.0 / pw * (v.xmax - v.xmin) / (jobs[k].box[3] - jobs[k].box[2]) - 1) / 2;
				v.ymax = v.ymin + (v.xmax - v.xmin) * ph / pw;
				sprintf(label, "%s, %s", k < 4 ? "polygons" : (k < 6 ? "NURBS ellipses" : "cubic B-spline blobs"), jobs[k].name);
			}
			else if (k < 10)
			{
				sc[0] = star;
				sp[0].curves = &sc[0]; sp[0].curve_count = 1; sp[0].closed = 1;
				d.subjects = sp; d.subject_count = 1;
				d.clip_type = QAWS_CLIP_UNION; d.fill_rule = k == 8 ? QAWS_FILL_NON_ZERO : QAWS_FILL_EVEN_ODD;
				v.xmin = -1.3; v.xmax = 1.3; v.ymin = -1.3 * ph / pw - 0.1; v.ymax = 1.3 * ph / pw - 0.1;
				sprintf(label, "self-crossing pentagram, %s", k == 8 ? "non-zero" : "even-odd");
			}
			else if (k == 10)
			{
				sc[0] = ring_o; sc[1] = ring_h; sc[2] = ring_i;
				for (i = 0; i < 3; i++) { sp[i].curves = &sc[i]; sp[i].curve_count = 1; sp[i].closed = 1; }
				cc = sq_b;
				cp.curves = &cc; cp.curve_count = 1; cp.closed = 1;
				d.subjects = sp; d.subject_count = 3; d.clips = &cp; d.clip_count = 1;
				d.clip_type = QAWS_CLIP_DIFFERENCE; d.fill_rule = QAWS_FILL_EVEN_ODD;
				v.xmin = -2.4; v.xmax = 3.4; v.ymin = -2.4 * ph / pw - 0.2; v.ymax = 3.4 * ph / pw + 0.2;
				v.ymin = -2.6; v.ymax = v.ymin + (v.xmax - v.xmin) * ph / pw;
				sprintf(label, "ring + island (even-odd) minus a square");
			}
			else
			{
				qaws_path_2d op;
				qaws_curve const* lc = line;
				sc[0] = circ_a;
				sp[0].curves = &sc[0]; sp[0].curve_count = 1; sp[0].closed = 1;
				op.curves = &lc; op.curve_count = 1; op.closed = 0;
				cc = sq_b;
				cp.curves = &sc[0]; cp.curve_count = 1; cp.closed = 1;
				d.open_subjects = &op; d.open_subject_count = 1;
				d.clips = &cp; d.clip_count = 1;
				d.clip_type = QAWS_CLIP_INTERSECTION; d.fill_rule = QAWS_FILL_NON_ZERO;
				v.xmin = -2.7; v.xmax = 2.7; v.ymin = -2.7 * ph / pw; v.ymax = 2.7 * ph / pw;
				sprintf(label, "open polyline clipped by an ellipse");
				svg_panel(&s, &v, label);
				r = NULL;
				qaws_clip_execute(&d, &r);
				bo_curve(&s, &v, circ_a, "#57606a", 1.2, 1);
				bo_curve(&s, &v, line, "#8c959f", 1.2, 1);
				if (r) bo_boundaries(&s, &v, r);
				sprintf(label, "%u open pieces inside", r ? qaws_clip_result_get_open_path_count(r) : 0);
				svg_text(&s, v.x0 + 10, v.y0 + ph - 12, 12, "#57606a", "start", label);
				qaws_clip_result_destroy(r);
				panel++;
				continue;
			}
			svg_panel(&s, &v, label);
			qaws_clip_execute(&d, &r);
			if (r)
			{
				qaws_scalar a = 0;
				unsigned int holes = 0, curves = 0;
				for (i = 0; i < qaws_clip_result_get_path_count(r); i++)
				{
					qaws_path_2d p;
					qaws_clip_result_get_path(r, i, &p);
					qaws_path_compute_area_2d(&p, &a);
					area += a;
					holes += qaws_clip_result_is_hole(r, i);
					curves += p.curve_count;
				}
				bo_fill(&s, &v, r, "#54aeff");
				for (i = 0; i < d.subject_count; i++)
					bo_curve(&s, &v, d.subjects[i].curves[0], "#57606a", 1.2, 1);
				for (i = 0; i < d.clip_count; i++)
					bo_curve(&s, &v, d.clips[i].curves[0], "#57606a", 1.2, 1);
				bo_boundaries(&s, &v, r);
				sprintf(label, "%u path%s (%u hole%s), %u curves, area %.6f", qaws_clip_result_get_path_count(r),
					qaws_clip_result_get_path_count(r) == 1 ? "" : "s", holes, holes == 1 ? "" : "s", curves, area);
				svg_text(&s, v.x0 + 10, v.y0 + ph - 12, 12, "#57606a", "start", label);
				qaws_clip_result_destroy(r);
			}
			panel++;
		}
	}
	svg_close(&s);
	printf("    boolean: %u panels in %.3f s\n", panel, (double)clock() / CLOCKS_PER_SEC - t0);
	qaws_curve_destroy(sq_a); qaws_curve_destroy(sq_b); qaws_curve_destroy(circ_a); qaws_curve_destroy(circ_b);
	qaws_curve_destroy(blob_a); qaws_curve_destroy(blob_b); qaws_curve_destroy(star);
	qaws_curve_destroy(ring_o); qaws_curve_destroy(ring_h); qaws_curve_destroy(ring_i); qaws_curve_destroy(line);
}
