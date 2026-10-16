/* Figure 18: regions across seams and round poles (qaws_surface_clip on
   closed surfaces). On a sphere: a geodesic square across the seam less a
   circle on it, a circle round the north pole with one across the seam,
   a band between two latitude lines xor a circle across the seam, and a
   south cap less a circle centred on the pole. Below each, the parameter
   plane, where the regions come back as pieces of one period. Last row: a
   band round a cylinder less a circle across its seam, and two circles on
   a torus across both of its seams. Uses figures 16 and 17's helpers. */

/* a surface curve in 3D (faded behind) and in the parameter plane, its u
   brought into one period and broken where it wraps */
static void ssf_curve(svg* s, qaws_surface_curve const* c, qaws_surface const* surf, sgf_view const* view, viewport const* vuv,
	double pu, double pv, char const* color)
{
	double uv[2 * 600], xy[2 * 600];
	qaws_range ru = qaws_surface_get_u_range(surf), rv = qaws_surface_get_v_range(surf);
	unsigned int m = sgf_sample(c, uv, 600), k, n = 0;
	if (pv > 0)
		for (k = 0; k < m; k++) { double w = fmod(uv[2 * k + 1] - rv.min_value, pv); if (w < 0) w += pv; uv[2 * k + 1] = w + rv.min_value; }
	sgf_curve3d(s, surf, view, uv, m, pu, color, 2.0, 0);
	if (!vuv)
		return;
	for (k = 0; k < m; k++)
	{
		double u = uv[2 * k];
		if (pu > 0) { u = fmod(u - ru.min_value, pu); if (u < 0) u += pu; u += ru.min_value; }
		if (n && (fabs(vx(vuv, u) - xy[2 * n - 2]) > 0.5 * vuv->w || fabs(vy(vuv, uv[2 * k + 1]) - xy[2 * n - 1]) > 0.5 * vuv->h))
		{
			svg_polyline(s, xy, (int)n, color, 1.8, 1.0, 0);
			n = 0;
		}
		xy[2 * n] = vx(vuv, u); xy[2 * n + 1] = vy(vuv, uv[2 * k + 1]); n++;
	}
	if (n > 1) svg_polyline(s, xy, (int)n, color, 1.8, 1.0, 0);
}

/* a region on the surface: its (u, v) loops sampled, scanned in thin rows
   of v (even-odd crossings), each row's runs drawn as surface quads where
   they face the viewer. (A loop's own outline projected would cover the
   back of a closed surface too.) */
static void ssf_fill3d(svg* s, qaws_surface_region const* r, qaws_surface const* surf, sgf_view const* view, char const* fill)
{
	unsigned int np = qaws_surface_region_get_path_count(r), pi, ci, k, total = 0, rows = 180, j, cap = 0;
	qaws_range ru = qaws_surface_get_u_range(surf), rv = qaws_surface_get_v_range(surf);
	double *pts = NULL, *xs = NULL, du = ((double)ru.max_value - ru.min_value) / 160;
	unsigned int* start = (unsigned int*)malloc(sizeof(unsigned int) * (np + 1));
	for (pi = 0; pi < np; pi++)
	{
		qaws_path_2d p;
		qaws_surface_region_get_path(r, pi, &p);
		start[pi] = total;
		for (ci = 0; ci < p.curve_count; ci++)
		{
			qaws_range rg = qaws_curve_get_parameter_range(p.curves[ci]);
			unsigned int spans = qaws_curve_get_span_count(p.curves[ci]), n = spans * 4 < 24 ? 24 : (spans * 4 > 400 ? 400 : spans * 4);
			for (k = 0; k < n; k++)
			{
				qaws_eval_result_2d e;
				if (total + 1 >= cap) { cap = cap ? 2 * cap : 4096; pts = (double*)realloc(pts, sizeof(double) * 2 * cap); }
				qaws_curve_evaluate_2d(&e, p.curves[ci], (qaws_scalar)(rg.min_value + (rg.max_value - rg.min_value) * k / n), QAWS_EVAL_FLAG_POSITION);
				pts[2 * total] = e.position.x; pts[2 * total + 1] = e.position.y; total++;
			}
		}
	}
	start[np] = total;
	xs = (double*)malloc(sizeof(double) * (total + 1));
	fprintf(s->f, "<g fill=\"%s\" stroke=\"%s\" stroke-width=\"0.6\" stroke-linejoin=\"round\" opacity=\"0.6\">\n", fill, fill);
	for (j = 0; j < rows; j++)
	{
		double va = rv.min_value + ((double)rv.max_value - rv.min_value) * j / rows, vb = rv.min_value + ((double)rv.max_value - rv.min_value) * (j + 1) / rows;
		double vm = 0.5 * (va + vb);
		unsigned int nx = 0, a, b;
		for (pi = 0; pi < np; pi++)
			for (k = start[pi]; k < start[pi + 1]; k++)
			{
				unsigned int l = k + 1 < start[pi + 1] ? k + 1 : start[pi];
				double y0 = pts[2 * k + 1], y1 = pts[2 * l + 1];
				if ((y0 <= vm) != (y1 <= vm))
					xs[nx++] = pts[2 * k] + (vm - y0) / (y1 - y0) * (pts[2 * l] - pts[2 * k]);
			}
		for (a = 1; a < nx; a++)
			for (b = a; b > 0 && xs[b - 1] > xs[b]; b--) { double t = xs[b]; xs[b] = xs[b - 1]; xs[b - 1] = t; }
		for (a = 0; a + 1 < nx; a += 2)
		{
			double u = xs[a];
			while (u < xs[a + 1])
			{
				double ue = u + du < xs[a + 1] ? u + du : xs[a + 1], q[8];
				double cu[4], cv[4];
				unsigned int m;
				if (sgf_facing(surf, view, 0.5 * (u + ue), vm, 0))
				{
					cu[0] = u; cv[0] = va; cu[1] = ue; cv[1] = va; cu[2] = ue; cv[2] = vb; cu[3] = u; cv[3] = vb;
					for (m = 0; m < 4; m++)
					{
						qaws_surface_eval_result e;
						qaws_surface_evaluate(surf, (qaws_scalar)cu[m], (qaws_scalar)cv[m], QAWS_SURFACE_EVAL_POSITION, &e);
						sgf_project(view, e.position, &q[2 * m], &q[2 * m + 1], NULL);
					}
					fprintf(s->f, "<path d=\"M%.2f,%.2f L%.2f,%.2f L%.2f,%.2f L%.2f,%.2f Z\"/>\n", q[0], q[1], q[2], q[3], q[4], q[5], q[6], q[7]);
				}
				u = ue;
			}
		}
	}
	fprintf(s->f, "</g>\n");
	free(pts); free(xs); free(start);
}

static void ssf_box(svg* s, viewport const* vuv, qaws_surface const* surf)
{
	qaws_range ru = qaws_surface_get_u_range(surf), rv = qaws_surface_get_v_range(surf);
	double box[10];
	box[0] = vx(vuv, ru.min_value); box[1] = vy(vuv, rv.min_value);
	box[2] = vx(vuv, ru.max_value); box[3] = vy(vuv, rv.min_value);
	box[4] = vx(vuv, ru.max_value); box[5] = vy(vuv, rv.max_value);
	box[6] = vx(vuv, ru.min_value); box[7] = vy(vuv, rv.max_value);
	box[8] = box[0]; box[9] = box[1];
	svg_polyline(s, box, 5, "#57606a", 1.0, 1.0, 0);
}

static qaws_surface_curve* ssf_line(qaws_surface const* s, double u0, double v0, double u1, double v1, qaws_curve** keep)
{
	qaws_bezier_desc bd;
	qaws_scalar cp[4];
	qaws_surface_curve* c = NULL;
	cp[0] = (qaws_scalar)u0; cp[1] = (qaws_scalar)v0; cp[2] = (qaws_scalar)u1; cp[3] = (qaws_scalar)v1;
	memset(&bd, 0, sizeof(bd));
	bd.dimension = QAWS_DIMENSION_2D; bd.degree = 1; bd.control_points = cp; bd.control_point_count = 2;
	qaws_curve_create_bezier(keep, &bd);
	qaws_surface_curve_create_uv(s, *keep, &c);
	return c;
}

static void ssf_path(qaws_surface_path* p, qaws_surface_curve const** cv, unsigned int n)
{
	p->curves = cv; p->curve_count = n; p->closed = 1;
}

static void demo_surface_seams(void)
{
	clock_t t0 = clock();
	svg s;
	double const pw = 300, ph = 260, pi = 3.14159265358979323846;
	unsigned int i, k;
	char label[220];
	qaws_status st;
	if (!svg_open(&s, "showcase/b2d18_surface_seams.svg", (int)(48 + 4 * (pw + 12)), (int)(110 + 3 * (ph + 44)),
		"Regions across seams and round poles",
		"qaws_surface_clip on closed surfaces: loops cut into period cells and united there; caps round poles; areas by the metric, exact topology"))
		return;
	/* the sphere: u round the axis (period 1), v from the south pole to the north */
	{
		qaws_curve* prof = sgf_arc(0, 0, 1, -pi / 2, pi / 2);
		qaws_surface* sph = sgf_revolve(prof);
		qaws_surface_curve *sq[4] = { 0 }, *cs = NULL, *cn = NULL, *cs2 = NULL, *cb = NULL, *csp = NULL, *lat[3] = { 0 };
		qaws_curve* keep[3] = { 0 };
		qaws_vec2 corner[4];
		corner[0] = sgf_uv(0.92, 0.38); corner[1] = sgf_uv(1.08, 0.38); corner[2] = sgf_uv(1.08, 0.62); corner[3] = sgf_uv(0.92, 0.62);
		for (i = 0; i < 4; i++) qaws_surface_curve_create_geodesic(sph, corner[i], corner[(i + 1) % 4], 0, &sq[i]);
		qaws_surface_curve_create_geodesic_arc(sph, sgf_uv(1.02, 0.56), (qaws_scalar)0.3, 0, (qaws_scalar)(2 * pi), 0, &cs);
		qaws_surface_curve_create_geodesic_arc(sph, sgf_uv(0.3, 1), (qaws_scalar)0.55, 0, (qaws_scalar)(2 * pi), 0, &cn);
		qaws_surface_curve_create_geodesic_arc(sph, sgf_uv(0.0, 0.72), (qaws_scalar)0.4, 0, (qaws_scalar)(2 * pi), 0, &cs2);
		qaws_surface_curve_create_geodesic_arc(sph, sgf_uv(0.03, 0.5), (qaws_scalar)0.42, 0, (qaws_scalar)(2 * pi), 0, &cb);
		qaws_surface_curve_create_geodesic_arc(sph, sgf_uv(0.0, 0), (qaws_scalar)0.35, 0, (qaws_scalar)(2 * pi), 0, &csp);
		/* latitude lines: going +u the inside is north */
		lat[0] = ssf_line(sph, 0, 0.42, 1, 0.42, &keep[0]);
		lat[1] = ssf_line(sph, 1, 0.6, 0, 0.6, &keep[1]);
		lat[2] = ssf_line(sph, 1, 0.3, 0, 0.3, &keep[2]);
		for (k = 0; k < 4; k++)
		{
			qaws_surface_curve const *pa[4], *pb[2];
			qaws_surface_path A, B;
			qaws_surface_region* r = NULL;
			qaws_scalar area = 0;
			viewport v3, vuv;
			sgf_view view;
			char const* name = "";
			unsigned int na = 0, nb = 0;
			qaws_clip_type op = QAWS_CLIP_UNION;
			v3.x0 = 24 + k * (pw + 12); v3.y0 = 104; v3.w = pw; v3.h = ph;
			vuv.x0 = v3.x0; vuv.y0 = 104 + ph + 44; vuv.w = pw; vuv.h = ph;
			vuv.xmin = -0.1; vuv.xmax = 1.1; vuv.ymin = 0.5 - 0.6 * ph / pw; vuv.ymax = 0.5 + 0.6 * ph / pw;
			/* the seam u = 0 faces the viewer at yaw -pi/2 */
			sgf_view_init(&view, -pi / 2 + 0.25, k == 1 ? 0.75 : (k == 3 ? -0.75 : 0.2), 108, v3.x0 + pw / 2, v3.y0 + ph / 2 + 8);
			if (k == 0)
			{
				for (i = 0; i < 4; i++) pa[i] = sq[i];
				na = 4; pb[0] = cs; nb = 1; op = QAWS_CLIP_DIFFERENCE;
				name = "a geodesic square on the seam less a circle";
			}
			else if (k == 1)
			{
				pa[0] = cn; na = 1; pb[0] = cs2; nb = 1; op = QAWS_CLIP_UNION;
				name = "a north pole circle with a circle on the seam";
			}
			else if (k == 2)
			{
				pa[0] = lat[0]; pa[1] = lat[1]; na = 2; pb[0] = cb; nb = 1; op = QAWS_CLIP_XOR;
				name = "a band of latitude lines xor a circle";
			}
			else
			{
				pa[0] = lat[2]; na = 1; pb[0] = csp; nb = 1; op = QAWS_CLIP_DIFFERENCE;
				name = "a south cap less a pole circle";
			}
			if (k == 2)
			{
				/* two loops, one operand: each its own path */
				qaws_surface_path AA[2];
				ssf_path(&AA[0], pa, 1); ssf_path(&AA[1], pa + 1, 1);
				ssf_path(&B, pb, nb);
				st = qaws_surface_clip(sph, op, QAWS_FILL_NON_ZERO, AA, 2, &B, 1, &r);
			}
			else
			{
				ssf_path(&A, pa, na); ssf_path(&B, pb, nb);
				st = qaws_surface_clip(sph, op, QAWS_FILL_NON_ZERO, &A, 1, &B, 1, &r);
			}
			if (st == QAWS_STATUS_OK) qaws_surface_region_get_area(r, &area);
			svg_panel(&s, &v3, name);
			sprintf(label, "(u, v): area %.4f", (double)area);
			svg_panel(&s, &vuv, label);
			sgf_mesh(&s, sph, &view, 12, 1);
			ssf_box(&s, &vuv, sph);
			if (r)
			{
				ssf_fill3d(&s, r, sph, &view, "#54aeff");
				sbf_fill(&s, r, sph, NULL, &vuv, 1, "#54aeff", "none");
			}
			for (i = 0; i < na; i++) ssf_curve(&s, pa[i], sph, &view, &vuv, 1, 0, "#0969da");
			for (i = 0; i < nb; i++) ssf_curve(&s, pb[i], sph, &view, &vuv, 1, 0, "#cf222e");
			if (st != QAWS_STATUS_OK) printf("figure 18: sphere panel %u: status %d\n", k, (int)st);
			qaws_surface_region_destroy(r);
		}
		for (i = 0; i < 4; i++) qaws_surface_curve_destroy(sq[i]);
		qaws_surface_curve_destroy(cs); qaws_surface_curve_destroy(cn); qaws_surface_curve_destroy(cs2);
		qaws_surface_curve_destroy(cb); qaws_surface_curve_destroy(csp);
		for (i = 0; i < 3; i++) { qaws_surface_curve_destroy(lat[i]); qaws_curve_destroy(keep[i]); }
		qaws_surface_destroy(sph);
		qaws_curve_destroy(prof);
	}
	/* a cylinder: a band round it less a circle across its seam */
	{
		double Rc = 0.8, H = 2.0;
		qaws_curve* circle = sgf_arc(0, 0, Rc, 0, 2 * pi);
		qaws_surface_extrusion_desc ed;
		qaws_surface* cy = NULL;
		qaws_curve *k1 = NULL, *k2 = NULL;
		qaws_surface_curve *lo, *hi, *c = NULL;
		qaws_surface_curve const *pa[2], *pb[1];
		qaws_surface_path A[2], B;
		qaws_surface_region* r = NULL;
		qaws_scalar area = 0;
		qaws_range ru;
		viewport v3, vuv;
		sgf_view view;
		double P;
		memset(&ed, 0, sizeof(ed));
		ed.profile = circle; ed.direction.z = 1; ed.length = (qaws_scalar)H;
		qaws_surface_create_extrusion(&ed, &cy);
		ru = qaws_surface_get_u_range(cy);
		P = (double)ru.max_value - ru.min_value;
		lo = ssf_line(cy, ru.min_value, 0.3, ru.max_value, 0.3, &k1);
		hi = ssf_line(cy, ru.max_value, 0.7, ru.min_value, 0.7, &k2);
		qaws_surface_curve_create_geodesic_arc(cy, sgf_uv(ru.min_value + 0.01 * P, 0.62), (qaws_scalar)0.45, 0, (qaws_scalar)(2 * pi), 0, &c);
		pa[0] = lo; pa[1] = hi; pb[0] = c;
		ssf_path(&A[0], pa, 1); ssf_path(&A[1], pa + 1, 1); ssf_path(&B, pb, 1);
		st = qaws_surface_clip(cy, QAWS_CLIP_DIFFERENCE, QAWS_FILL_NON_ZERO, A, 2, &B, 1, &r);
		if (st == QAWS_STATUS_OK) qaws_surface_region_get_area(r, &area);
		v3.x0 = 24; v3.y0 = 104 + 2 * (ph + 44); v3.w = pw; v3.h = ph;
		vuv.x0 = v3.x0 + pw + 12; vuv.y0 = v3.y0; vuv.w = pw; vuv.h = ph;
		vuv.xmin = ru.min_value - 0.1 * P; vuv.xmax = ru.max_value + 0.1 * P;
		vuv.ymin = -0.1; vuv.ymax = 1.1;
		/* (u, v) scaled to the panel: the circle is an ellipse there */
		sgf_view_init(&view, -pi / 2 + 0.3, 0.35, 95, v3.x0 + pw / 2, v3.y0 + ph / 2 + 8);
		view.o[2] = H / 2;
		svg_panel(&s, &v3, "a cylinder: a band less a circle on the seam");
		{
			/* unrolled: the band less the disc, but for its segment above the band */
			double rr = 0.45, d = (0.7 - 0.62) * H, seg = rr * rr * acos(d / rr) - d * sqrt(rr * rr - d * d);
			sprintf(label, "area %.4f (unrolled, closed form %.4f)", (double)area, 2 * pi * Rc * 0.4 * H - (pi * rr * rr - seg));
		}
		svg_panel(&s, &vuv, label);
		sgf_mesh(&s, cy, &view, 12, P);
		ssf_box(&s, &vuv, cy);
		if (r)
		{
			ssf_fill3d(&s, r, cy, &view, "#54aeff");
			sbf_fill(&s, r, cy, NULL, &vuv, P, "#54aeff", "none");
		}
		ssf_curve(&s, lo, cy, &view, &vuv, P, 0, "#0969da");
		ssf_curve(&s, hi, cy, &view, &vuv, P, 0, "#0969da");
		ssf_curve(&s, c, cy, &view, &vuv, P, 0, "#cf222e");
		if (st != QAWS_STATUS_OK) printf("figure 18: cylinder: status %d\n", (int)st);
		qaws_surface_region_destroy(r);
		qaws_surface_curve_destroy(lo); qaws_surface_curve_destroy(hi); qaws_surface_curve_destroy(c);
		qaws_curve_destroy(k1); qaws_curve_destroy(k2);
		qaws_surface_destroy(cy);
		qaws_curve_destroy(circle);
	}
	/* a torus: two circles across both seams at once (the corner of its domain) */
	{
		qaws_curve* prof = sgf_arc(1.0, 0, 0.42, 0, 2 * pi);
		qaws_surface* tor = sgf_revolve(prof);
		qaws_surface_curve *c1 = NULL, *c2 = NULL;
		qaws_surface_curve const *pa[1], *pb[1];
		qaws_surface_path A, B;
		qaws_surface_region* r = NULL;
		qaws_scalar area = 0, a1 = 0, a2 = 0;
		viewport v3, vuv;
		sgf_view view;
		qaws_surface_curve_create_geodesic_arc(tor, sgf_uv(0.0, 0.0), (qaws_scalar)0.32, 0, (qaws_scalar)(2 * pi), 0, &c1);
		qaws_surface_curve_create_geodesic_arc(tor, sgf_uv(0.05, 0.93), (qaws_scalar)0.28, 0, (qaws_scalar)(2 * pi), 0, &c2);
		pa[0] = c1; pb[0] = c2;
		ssf_path(&A, pa, 1); ssf_path(&B, pb, 1);
		{
			qaws_surface_region* q = NULL;
			if (qaws_surface_clip(tor, QAWS_CLIP_NONE, QAWS_FILL_NON_ZERO, &A, 1, NULL, 0, &q) == QAWS_STATUS_OK) qaws_surface_region_get_area(q, &a1);
			qaws_surface_region_destroy(q); q = NULL;
			if (qaws_surface_clip(tor, QAWS_CLIP_NONE, QAWS_FILL_NON_ZERO, &B, 1, NULL, 0, &q) == QAWS_STATUS_OK) qaws_surface_region_get_area(q, &a2);
			qaws_surface_region_destroy(q);
		}
		st = qaws_surface_clip(tor, QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO, &A, 1, &B, 1, &r);
		if (st == QAWS_STATUS_OK) qaws_surface_region_get_area(r, &area);
		v3.x0 = 24 + 2 * (pw + 12); v3.y0 = 104 + 2 * (ph + 44); v3.w = pw; v3.h = ph;
		vuv.x0 = v3.x0 + pw + 12; vuv.y0 = v3.y0; vuv.w = pw; vuv.h = ph;
		vuv.xmin = -0.1; vuv.xmax = 1.1; vuv.ymin = 0.5 - 0.6 * ph / pw; vuv.ymax = 0.5 + 0.6 * ph / pw;
		sgf_view_init(&view, -pi / 2 + 0.35, 0.55, 105, v3.x0 + pw / 2, v3.y0 + ph / 2 + 4);
		svg_panel(&s, &v3, "a torus: circles across both seams, united");
		sprintf(label, "area %.4f (the two: %.4f, %.4f)", (double)area, (double)a1, (double)a2);
		svg_panel(&s, &vuv, label);
		sgf_mesh(&s, tor, &view, 16, 1);
		ssf_box(&s, &vuv, tor);
		if (r)
		{
			ssf_fill3d(&s, r, tor, &view, "#54aeff");
			sbf_fill(&s, r, tor, NULL, &vuv, 1, "#54aeff", "none");
		}
		ssf_curve(&s, c1, tor, &view, &vuv, 1, 1, "#0969da");
		ssf_curve(&s, c2, tor, &view, &vuv, 1, 1, "#cf222e");
		if (st != QAWS_STATUS_OK) printf("figure 18: torus: status %d\n", (int)st);
		qaws_surface_region_destroy(r);
		qaws_surface_curve_destroy(c1); qaws_surface_curve_destroy(c2);
		qaws_surface_destroy(tor);
		qaws_curve_destroy(prof);
	}
	svg_close(&s);
	printf("figure 18: regions across seams and poles, %.1f s\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
}
