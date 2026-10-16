/* Figure 19: geodesic offsets on surfaces (qaws_surface_offset). On a
   sphere: a geodesic triangle grown (strips swept by conormal geodesics,
   geodesic circles round the corners), a square with a circle shrunk, the
   tube of an open geodesic polyline (round ends), and a band across the
   seam grown. Below, their parameter planes. Last row: rings of offsets of
   a region on a wavy bicubic patch, and a tube along a geodesic polyline
   across both seams of a torus. Uses figures 16 to 18's helpers. */

/* a region's boundary in 3D and (u, v), the cuts along seams left out */
static void sof_outline(svg* s, qaws_surface_region const* r, qaws_surface const* surf, sgf_view const* view, viewport const* vuv,
	double pu, double pv, char const* color, double width)
{
	unsigned int np = qaws_surface_region_get_path_count(r), pi, ci, k;
	qaws_range ru = qaws_surface_get_u_range(surf), rv = qaws_surface_get_v_range(surf);
	for (pi = 0; pi < np; pi++)
	{
		qaws_path_2d p;
		qaws_surface_region_get_path(r, pi, &p);
		for (ci = 0; ci < p.curve_count; ci++)
		{
			qaws_range rg = qaws_curve_get_parameter_range(p.curves[ci]);
			if (qaws_curve_get_degree(p.curves[ci]) == 1)
			{
				/* a polyline (lines merged by the clip, cuts along seams among
				   them): segment by segment, the cuts left out */
				unsigned int nc = 0, seg;
				qaws_scalar probe[2], *cp;
				qaws_curve_get_control_points(probe, &nc, p.curves[ci], 0);
				cp = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * (nc + 1));
				qaws_curve_get_control_points(cp, &nc, p.curves[ci], nc);
				for (seg = 0; seg + 1 < nc; seg++)
				{
					double a0 = cp[2 * seg], a1 = cp[2 * seg + 1], b0 = cp[2 * seg + 2], b1 = cp[2 * seg + 3];
					double eu = 1e-9 * (ru.max_value - ru.min_value), ev = 1e-9 * (rv.max_value - rv.min_value);
					double uv[2 * 49], xy[2 * 49];
					int cut = 0;
					if (pu > 0 && fabs(a0 - b0) < eu && (fabs(a0 - ru.min_value) < eu || fabs(a0 - ru.min_value - pu) < eu)) cut = 1;
					if (pv > 0 && fabs(a1 - b1) < ev && (fabs(a1 - rv.min_value) < ev || fabs(a1 - rv.min_value - pv) < ev)) cut = 1;
					if (fabs(a1 - b1) < ev && (fabs(a1 - rv.min_value) < ev || fabs(a1 - rv.max_value) < ev)) cut = 1;
					if (cut)
						continue;
					for (k = 0; k <= 48; k++) { uv[2 * k] = a0 + (b0 - a0) * k / 48; uv[2 * k + 1] = a1 + (b1 - a1) * k / 48; }
					sgf_curve3d(s, surf, view, uv, 49, pu, color, width, 0);
					if (vuv)
					{
						for (k = 0; k <= 48; k++) { xy[2 * k] = vx(vuv, uv[2 * k]); xy[2 * k + 1] = vy(vuv, uv[2 * k + 1]); }
						svg_polyline(s, xy, 49, color, 0.8 * width, 1.0, 0);
					}
				}
				free(cp);
				continue;
			}
			unsigned int spans = qaws_curve_get_span_count(p.curves[ci]), n = spans * 4 < 32 ? 32 : (spans * 4 > 600 ? 600 : spans * 4);
			double* uv = (double*)malloc(sizeof(double) * 2 * (n + 1));
			double* xy = (double*)malloc(sizeof(double) * 2 * (n + 1));
			for (k = 0; k <= n; k++)
			{
				qaws_eval_result_2d e;
				qaws_curve_evaluate_2d(&e, p.curves[ci], (qaws_scalar)(rg.min_value + (rg.max_value - rg.min_value) * k / n), QAWS_EVAL_FLAG_POSITION);
				uv[2 * k] = e.position.x; uv[2 * k + 1] = e.position.y;
			}
			{
				/* on a seam (or the domain's edge) at both ends and the middle: a cut, not the boundary */
				double const* a = uv, *b = uv + 2 * n, *m = uv + 2 * (n / 2);
				double eu = 1e-9 * (ru.max_value - ru.min_value), ev = 1e-9 * (rv.max_value - rv.min_value);
				int seam = 0;
				if (pu > 0)
				{
					double e0 = ru.min_value, e1 = ru.min_value + pu;
					if ((fabs(a[0] - e0) < eu && fabs(b[0] - e0) < eu && fabs(m[0] - e0) < eu) || (fabs(a[0] - e1) < eu && fabs(b[0] - e1) < eu && fabs(m[0] - e1) < eu)) seam = 1;
				}
				if (pv > 0)
				{
					double e0 = rv.min_value, e1 = rv.min_value + pv;
					if ((fabs(a[1] - e0) < ev && fabs(b[1] - e0) < ev && fabs(m[1] - e0) < ev) || (fabs(a[1] - e1) < ev && fabs(b[1] - e1) < ev && fabs(m[1] - e1) < ev)) seam = 1;
				}
				/* a pole: a whole edge of a sphere's domain */
				if ((fabs(a[1] - rv.min_value) < ev && fabs(b[1] - rv.min_value) < ev) || (fabs(a[1] - rv.max_value) < ev && fabs(b[1] - rv.max_value) < ev)) seam = 1;
				if (!seam)
				{
					sgf_curve3d(s, surf, view, uv, n + 1, pu, color, width, 0);
					if (vuv)
					{
						for (k = 0; k <= n; k++) { xy[2 * k] = vx(vuv, uv[2 * k]); xy[2 * k + 1] = vy(vuv, uv[2 * k + 1]); }
						svg_polyline(s, xy, (int)n + 1, color, 0.8 * width, 1.0, 0);
					}
				}
			}
			free(uv); free(xy);
		}
	}
}

static qaws_surface_region* sof_offset(qaws_surface const* s, qaws_surface_curve* const* cv, unsigned int n, int closed, double d)
{
	qaws_surface_curve const* pc[8];
	qaws_surface_path P;
	qaws_surface_region* r = NULL;
	unsigned int i;
	qaws_status st;
	for (i = 0; i < n; i++) pc[i] = cv[i];
	P.curves = pc; P.curve_count = n; P.closed = closed;
	st = qaws_surface_offset(s, QAWS_FILL_NON_ZERO, &P, 1, (qaws_scalar)d, 0, &r);
	if (st != QAWS_STATUS_OK) printf("figure 19: offset by %.3f: status %d\n", d, (int)st);
	return r;
}

static void demo_surface_offset(void)
{
	clock_t t0 = clock();
	svg s;
	double const pw = 300, ph = 260, pi = 3.14159265358979323846;
	unsigned int i, k;
	char label[220];
	if (!svg_open(&s, "showcase/b2d19_surface_offset.svg", (int)(48 + 4 * (pw + 12)), (int)(110 + 3 * (ph + 44)),
		"Geodesic offsets on surfaces",
		"qaws_surface_offset: the boundary carried along its conormal geodesics, geodesic circles round corners and ends; grown, shrunk, tubes; across seams"))
		return;
	/* the sphere */
	{
		qaws_curve* prof = sgf_arc(0, 0, 1, -pi / 2, pi / 2);
		qaws_surface* sph = sgf_revolve(prof);
		qaws_surface_curve *tri[3] = { 0 }, *sq[4] = { 0 }, *circ = NULL, *pl[3] = { 0 }, *lat[2] = { 0 }, *sc[2] = { 0 };
		qaws_curve* keep[2] = { 0 };
		qaws_vec2 a[4];
		a[0] = sgf_uv(0.92, 0.36); a[1] = sgf_uv(1.12, 0.40); a[2] = sgf_uv(1.0, 0.66);
		for (i = 0; i < 3; i++) qaws_surface_curve_create_geodesic(sph, a[i], a[(i + 1) % 3], 0, &tri[i]);
		a[0] = sgf_uv(0.9, 0.36); a[1] = sgf_uv(1.06, 0.36); a[2] = sgf_uv(1.06, 0.6); a[3] = sgf_uv(0.9, 0.6);
		for (i = 0; i < 4; i++) qaws_surface_curve_create_geodesic(sph, a[i], a[(i + 1) % 4], 0, &sq[i]);
		qaws_surface_curve_create_geodesic_arc(sph, sgf_uv(1.08, 0.6), (qaws_scalar)0.3, 0, (qaws_scalar)(2 * pi), 0, &circ);
		a[0] = sgf_uv(0.88, 0.3); a[1] = sgf_uv(0.98, 0.62); a[2] = sgf_uv(1.1, 0.42);
		for (i = 0; i < 2; i++) qaws_surface_curve_create_geodesic(sph, a[i], a[i + 1], 0, &pl[i]);
		lat[0] = ssf_line(sph, 0, 0.45, 1, 0.45, &keep[0]);
		lat[1] = ssf_line(sph, 1, 0.55, 0, 0.55, &keep[1]);
		for (k = 0; k < 4; k++)
		{
			viewport v3, vuv;
			sgf_view view;
			qaws_surface_region* r = NULL;
			qaws_surface_region* base = NULL;
			qaws_scalar area = 0, area0 = 0;
			char const* name = "";
			v3.x0 = 24 + k * (pw + 12); v3.y0 = 104; v3.w = pw; v3.h = ph;
			vuv.x0 = v3.x0; vuv.y0 = 104 + ph + 44; vuv.w = pw; vuv.h = ph;
			vuv.xmin = -0.1; vuv.xmax = 1.1; vuv.ymin = 0.5 - 0.6 * ph / pw; vuv.ymax = 0.5 + 0.6 * ph / pw;
			sgf_view_init(&view, -pi / 2 + 0.2, 0.25, 108, v3.x0 + pw / 2, v3.y0 + ph / 2 + 8);
			if (k == 0)
			{
				r = sof_offset(sph, tri, 3, 1, 0.15);
				{
					qaws_surface_curve const* pc[3];
					qaws_surface_path P;
					for (i = 0; i < 3; i++) pc[i] = tri[i];
					ssf_path(&P, pc, 3);
					qaws_surface_clip(sph, QAWS_CLIP_NONE, QAWS_FILL_NON_ZERO, &P, 1, NULL, 0, &base);
				}
				name = "a geodesic triangle grown by 0.15";
			}
			else if (k == 1)
			{
				/* the square with the circle, shrunk */
				qaws_surface_curve const *pa[4], *pb[1];
				qaws_surface_path A[2];
				for (i = 0; i < 4; i++) pa[i] = sq[i];
				pb[0] = circ;
				ssf_path(&A[0], pa, 4); ssf_path(&A[1], pb, 1);
				qaws_surface_clip(sph, QAWS_CLIP_NONE, QAWS_FILL_NON_ZERO, A, 2, NULL, 0, &base);
				if (qaws_surface_offset(sph, QAWS_FILL_NON_ZERO, A, 2, (qaws_scalar)-0.08, 0, &r) != QAWS_STATUS_OK) printf("figure 19: shrink failed\n");
				name = "a square with a circle, shrunk by 0.08";
			}
			else if (k == 2)
			{
				r = sof_offset(sph, pl, 2, 0, 0.1);
				name = "an open geodesic polyline's tube, 0.1";
			}
			else
			{
				qaws_surface_curve const *p0[1], *p1[1];
				qaws_surface_path B[2];
				p0[0] = lat[0]; p1[0] = lat[1];
				ssf_path(&B[0], p0, 1); ssf_path(&B[1], p1, 1);
				qaws_surface_clip(sph, QAWS_CLIP_NONE, QAWS_FILL_NON_ZERO, B, 2, NULL, 0, &base);
				if (qaws_surface_offset(sph, QAWS_FILL_NON_ZERO, B, 2, (qaws_scalar)0.25, 0, &r) != QAWS_STATUS_OK) printf("figure 19: band failed\n");
				view.pitch = 0.45;
				sgf_view_init(&view, -pi / 2 + 0.2, 0.45, 108, v3.x0 + pw / 2, v3.y0 + ph / 2 + 8);
				name = "a band round the sphere grown by 0.25";
			}
			if (r) qaws_surface_region_get_area(r, &area);
			if (base) qaws_surface_region_get_area(base, &area0);
			svg_panel(&s, &v3, name);
			if (base) sprintf(label, "(u, v): area %.4f from %.4f", (double)area, (double)area0);
			else sprintf(label, "(u, v): area %.4f", (double)area);
			svg_panel(&s, &vuv, label);
			sgf_mesh(&s, sph, &view, 12, 1);
			ssf_box(&s, &vuv, sph);
			if (r)
			{
				ssf_fill3d(&s, r, sph, &view, "#54aeff");
				sbf_fill(&s, r, sph, NULL, &vuv, 1, "#54aeff", "none");
				sof_outline(&s, r, sph, &view, &vuv, 1, 0, "#0969da", 1.8);
			}
			if (base) sof_outline(&s, base, sph, &view, &vuv, 1, 0, "#cf222e", 1.6);
			if (k == 2) for (i = 0; i < 2; i++) ssf_curve(&s, pl[i], sph, &view, &vuv, 1, 0, "#cf222e");
			qaws_surface_region_destroy(r);
			qaws_surface_region_destroy(base);
		}
		for (i = 0; i < 3; i++) qaws_surface_curve_destroy(tri[i]);
		for (i = 0; i < 4; i++) qaws_surface_curve_destroy(sq[i]);
		for (i = 0; i < 2; i++) { qaws_surface_curve_destroy(pl[i]); qaws_surface_curve_destroy(lat[i]); qaws_curve_destroy(keep[i]); qaws_surface_curve_destroy(sc[i]); }
		qaws_surface_curve_destroy(circ);
		qaws_surface_destroy(sph);
		qaws_curve_destroy(prof);
	}
	/* a wavy patch: rings of offsets of a geodesic square with a circle */
	{
		static double const z[16] = { 0, 0.3, -0.2, 0, 0.35, 0.8, 0.5, 0.1, -0.3, 0.6, 0.9, 0.35, 0, -0.15, 0.3, 0 };
		static double const dist[5] = { -0.16, -0.08, 0, 0.08, 0.16 };
		static char const* const col[5] = { "#8250df", "#bf3989", "#cf222e", "#bc4c00", "#9a6700" };
		qaws_vec3 cp[16];
		qaws_surface_bezier_desc bd;
		qaws_surface* patch = NULL;
		qaws_surface_curve *sq[4] = { 0 }, *circ = NULL;
		qaws_surface_curve const *pa[4], *pb[1];
		qaws_surface_path A[2];
		qaws_vec2 corner[4];
		viewport v3, vuv;
		sgf_view view;
		for (i = 0; i < 16; i++) { cp[i].x = (qaws_scalar)((i / 4) / 3.0 * 2); cp[i].y = (qaws_scalar)((i % 4) / 3.0 * 2); cp[i].z = (qaws_scalar)z[i]; }
		memset(&bd, 0, sizeof(bd));
		bd.u_degree = 3; bd.v_degree = 3; bd.control_points = cp; bd.u_point_count = 4; bd.v_point_count = 4;
		qaws_surface_create_bezier(&patch, &bd);
		corner[0] = sgf_uv(0.25, 0.25); corner[1] = sgf_uv(0.6, 0.22); corner[2] = sgf_uv(0.62, 0.58); corner[3] = sgf_uv(0.28, 0.6);
		for (i = 0; i < 4; i++) { qaws_surface_curve_create_geodesic(patch, corner[i], corner[(i + 1) % 4], 0, &sq[i]); pa[i] = sq[i]; }
		qaws_surface_curve_create_geodesic_arc(patch, sgf_uv(0.63, 0.62), (qaws_scalar)0.36, 0, (qaws_scalar)(2 * pi), 0, &circ);
		pb[0] = circ;
		ssf_path(&A[0], pa, 4); ssf_path(&A[1], pb, 1);
		v3.x0 = 24; v3.y0 = 104 + 2 * (ph + 44); v3.w = pw; v3.h = ph;
		vuv.x0 = v3.x0 + pw + 12; vuv.y0 = v3.y0; vuv.w = pw; vuv.h = ph;
		vuv.xmin = -0.1; vuv.xmax = 1.1; vuv.ymin = 0.5 - 0.6 * ph / pw; vuv.ymax = 0.5 + 0.6 * ph / pw;
		sgf_view_init(&view, 0.7, 0.75, 100, v3.x0 + pw / 2, v3.y0 + ph / 2 + 10);
		view.o[0] = 1; view.o[1] = 1; view.o[2] = 0.3;
		svg_panel(&s, &v3, "a wavy patch: offsets -0.16 to 0.16 of a region");
		svg_panel(&s, &vuv, "(u, v): the same rings, bent by the metric");
		sgf_mesh(&s, patch, &view, 10, 0);
		ssf_box(&s, &vuv, patch);
		for (k = 0; k < 5; k++)
		{
			qaws_surface_region* r = NULL;
			if (dist[k] == 0)
				qaws_surface_clip(patch, QAWS_CLIP_NONE, QAWS_FILL_NON_ZERO, A, 2, NULL, 0, &r);
			else if (qaws_surface_offset(patch, QAWS_FILL_NON_ZERO, A, 2, (qaws_scalar)dist[k], 0, &r) != QAWS_STATUS_OK)
				printf("figure 19: patch offset %.2f failed\n", dist[k]);
			if (r)
			{
				if (dist[k] == 0)
				{
					ssf_fill3d(&s, r, patch, &view, "#54aeff");
					sbf_fill(&s, r, patch, NULL, &vuv, 0, "#54aeff", "none");
				}
				sof_outline(&s, r, patch, &view, &vuv, 0, 0, col[k], dist[k] == 0 ? 2.0 : 1.5);
			}
			qaws_surface_region_destroy(r);
		}
		for (i = 0; i < 4; i++) qaws_surface_curve_destroy(sq[i]);
		qaws_surface_curve_destroy(circ);
		qaws_surface_destroy(patch);
	}
	/* a torus: the tube of a geodesic polyline across both seams */
	{
		qaws_curve* prof = sgf_arc(1.0, 0, 0.42, 0, 2 * pi);
		qaws_surface* tor = sgf_revolve(prof);
		qaws_surface_curve* pl[3] = { 0 };
		qaws_surface_region* r;
		qaws_vec2 a[4];
		qaws_scalar area = 0;
		viewport v3, vuv;
		sgf_view view;
		a[0] = sgf_uv(0.86, 0.2); a[1] = sgf_uv(0.97, 0.9); a[2] = sgf_uv(1.08, 1.12); a[3] = sgf_uv(1.16, 0.95);
		for (i = 0; i < 3; i++) qaws_surface_curve_create_geodesic(tor, a[i], a[i + 1], 0, &pl[i]);
		r = sof_offset(tor, pl, 3, 0, 0.09);
		if (r) qaws_surface_region_get_area(r, &area);
		v3.x0 = 24 + 2 * (pw + 12); v3.y0 = 104 + 2 * (ph + 44); v3.w = pw; v3.h = ph;
		vuv.x0 = v3.x0 + pw + 12; vuv.y0 = v3.y0; vuv.w = pw; vuv.h = ph;
		vuv.xmin = -0.1; vuv.xmax = 1.1; vuv.ymin = 0.5 - 0.6 * ph / pw; vuv.ymax = 0.5 + 0.6 * ph / pw;
		sgf_view_init(&view, -pi / 2 + 0.35, 0.55, 105, v3.x0 + pw / 2, v3.y0 + ph / 2 + 4);
		svg_panel(&s, &v3, "a torus: a polyline's tube across both seams");
		sprintf(label, "(u, v): area %.4f", (double)area);
		svg_panel(&s, &vuv, label);
		sgf_mesh(&s, tor, &view, 16, 1);
		ssf_box(&s, &vuv, tor);
		if (r)
		{
			ssf_fill3d(&s, r, tor, &view, "#54aeff");
			sbf_fill(&s, r, tor, NULL, &vuv, 1, "#54aeff", "none");
			sof_outline(&s, r, tor, &view, &vuv, 1, 1, "#0969da", 1.6);
		}
		for (i = 0; i < 3; i++) ssf_curve(&s, pl[i], tor, &view, &vuv, 1, 1, "#cf222e");
		qaws_surface_region_destroy(r);
		for (i = 0; i < 3; i++) qaws_surface_curve_destroy(pl[i]);
		qaws_surface_destroy(tor);
		qaws_curve_destroy(prof);
	}
	svg_close(&s);
	printf("figure 19: geodesic offsets on surfaces, %.1f s\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
}
