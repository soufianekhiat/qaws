/* Figure 17: Booleans of regions on surfaces (qaws_surface_clip). A
   geodesic square and a geodesic circle on a wavy bicubic patch: union,
   intersection, difference and xor, shaded on the surface and in the
   parameter plane, each with its surface area. Below, on a sphere: a
   lens of two geodesic circles, a geodesic triangle less a circle, and three
   circles' union. Uses figure 16's drawing helpers. */

/* the region's loops as one SVG path (even-odd), in 3D or in (u, v) */
static void sbf_fill(svg* s, qaws_surface_region const* r, qaws_surface const* surf, sgf_view const* view, viewport const* vuv,
	double pu, char const* fill, char const* stroke)
{
	unsigned int np = qaws_surface_region_get_path_count(r), pi, ci, k;
	qaws_range ru = qaws_surface_get_u_range(surf);
	fprintf(s->f, "<path fill=\"%s\" fill-opacity=\"0.55\" fill-rule=\"evenodd\" stroke=\"%s\" stroke-width=\"1.6\" d=\"", fill, stroke);
	for (pi = 0; pi < np; pi++)
	{
		qaws_path_2d p;
		int first = 1;
		qaws_surface_region_get_path(r, pi, &p);
		for (ci = 0; ci < p.curve_count; ci++)
		{
			qaws_range rg = qaws_curve_get_parameter_range(p.curves[ci]);
			unsigned int spans = qaws_curve_get_span_count(p.curves[ci]), n = spans * 4 < 24 ? 24 : (spans * 4 > 400 ? 400 : spans * 4);
			for (k = 0; k <= n; k++)
			{
				qaws_eval_result_2d e;
				double sx, sy;
				qaws_curve_evaluate_2d(p.curves[ci], (qaws_scalar)(rg.min_value + (rg.max_value - rg.min_value) * k / n), QAWS_EVAL_FLAG_POSITION, &e);
				if (view)
				{
					qaws_surface_eval_result q;
					double uu = e.position.x;
					if (pu > 0) { uu = fmod(uu - ru.min_value, pu); if (uu < 0) uu += pu; uu += ru.min_value; }
					qaws_surface_evaluate(surf, (qaws_scalar)uu, e.position.y, QAWS_SURFACE_EVAL_POSITION, &q);
					sgf_project(view, q.position, &sx, &sy, NULL);
				}
				else
				{
					sx = vx(vuv, e.position.x);
					sy = vy(vuv, e.position.y);
				}
				fprintf(s->f, "%s%.2f,%.2f ", first ? "M" : "L", sx, sy);
				first = 0;
			}
		}
		fprintf(s->f, "Z ");
	}
	fprintf(s->f, "\"/>\n");
}

static void demo_surface_boolean(void)
{
	clock_t t0 = clock();
	svg s;
	double const pw = 300, ph = 260, pi = 3.14159265358979323846;
	static char const* const opname[4] = { "union", "intersection", "difference", "xor" };
	static qaws_clip_type const ops[4] = { QAWS_CLIP_UNION, QAWS_CLIP_INTERSECTION, QAWS_CLIP_DIFFERENCE, QAWS_CLIP_XOR };
	static double const z[16] = { 0, 0.3, -0.2, 0, 0.35, 0.8, 0.5, 0.1, -0.3, 0.6, 0.9, 0.35, 0, -0.15, 0.3, 0 };
	qaws_vec3 cp[16];
	qaws_surface_bezier_desc bd;
	qaws_surface* patch = NULL;
	qaws_surface_curve *sq[4] = { 0 }, *circle = NULL;
	qaws_surface_curve const *psq[4], *pci[1];
	qaws_surface_path SQ, CI;
	qaws_vec2 corner[4];
	unsigned int i, k;
	char label[200];
	if (!svg_open(&s, "showcase/b2d17_surface_boolean.svg", (int)(48 + 4 * (pw + 12)), (int)(110 + 3 * (ph + 44)),
		"Booleans of regions on surfaces, with their surface areas",
		"qaws_surface_clip: geodesic edges and circles combined through their parameter-plane images (exact topology); areas by Green's theorem with the metric"))
		return;
	for (i = 0; i < 16; i++) { cp[i].x = (qaws_scalar)((i / 4) / 3.0 * 2); cp[i].y = (qaws_scalar)((i % 4) / 3.0 * 2); cp[i].z = (qaws_scalar)z[i]; }
	memset(&bd, 0, sizeof(bd));
	bd.u_degree = 3; bd.v_degree = 3; bd.control_points = cp; bd.u_point_count = 4; bd.v_point_count = 4;
	qaws_surface_create_bezier(&bd, &patch);
	corner[0] = sgf_uv(0.15, 0.2); corner[1] = sgf_uv(0.62, 0.15); corner[2] = sgf_uv(0.67, 0.62); corner[3] = sgf_uv(0.2, 0.67);
	for (i = 0; i < 4; i++) { qaws_surface_curve_create_geodesic(patch, corner[i], corner[(i + 1) % 4], 0, &sq[i]); psq[i] = sq[i]; }
	qaws_surface_curve_create_geodesic_arc(patch, sgf_uv(0.62, 0.62), (qaws_scalar)0.55, 0, (qaws_scalar)(2 * pi), 0, &circle);
	pci[0] = circle;
	SQ.curves = psq; SQ.curve_count = 4; SQ.closed = 1;
	CI.curves = pci; CI.curve_count = 1; CI.closed = 1;
	for (k = 0; k < 4; k++)
	{
		qaws_surface_region* r = NULL;
		qaws_scalar area = 0;
		viewport v3, vuv;
		sgf_view view;
		v3.x0 = 24 + k * (pw + 12); v3.y0 = 104; v3.w = pw; v3.h = ph;
		vuv.x0 = v3.x0; vuv.y0 = 104 + ph + 44; vuv.w = pw; vuv.h = ph;
		vuv.xmin = -0.05; vuv.xmax = 1.05; vuv.ymin = -0.05 - 0.55 * (ph / pw - 1); vuv.ymax = 1.05;
		vuv.xmin = -0.1; vuv.xmax = 1.1; vuv.ymin = 0.5 - 0.6 * ph / pw; vuv.ymax = 0.5 + 0.6 * ph / pw;
		sgf_view_init(&view, 0.7, 0.75, 100, v3.x0 + pw / 2, v3.y0 + ph / 2 + 10);
		view.o[0] = 1; view.o[1] = 1; view.o[2] = 0.3;
		qaws_surface_clip(patch, ops[k], QAWS_FILL_NON_ZERO, &SQ, 1, &CI, 1, &r);
		qaws_surface_region_get_area(r, &area);
		sprintf(label, "%s: area %.4f", opname[k], (double)area);
		svg_panel(&s, &v3, label);
		svg_panel(&s, &vuv, k == 0 ? "the same in the parameter plane" : NULL);
		sgf_mesh(&s, patch, &view, 10, 0);
		sbf_fill(&s, r, patch, &view, NULL, 0, "#54aeff", "#0969da");
		{
			double box[10];
			box[0] = vx(&vuv, 0); box[1] = vy(&vuv, 0); box[2] = vx(&vuv, 1); box[3] = vy(&vuv, 0);
			box[4] = vx(&vuv, 1); box[5] = vy(&vuv, 1); box[6] = vx(&vuv, 0); box[7] = vy(&vuv, 1); box[8] = box[0]; box[9] = box[1];
			svg_polyline(&s, box, 5, "#57606a", 1.0, 1.0, 0);
		}
		sbf_fill(&s, r, patch, NULL, &vuv, 0, "#54aeff", "#0969da");
		qaws_surface_region_destroy(r);
	}
	for (i = 0; i < 4; i++) qaws_surface_curve_destroy(sq[i]);
	qaws_surface_curve_destroy(circle);
	qaws_surface_destroy(patch);
	/* the sphere */
	{
		qaws_curve* prof = sgf_arc(0, 0, 1, -pi / 2, pi / 2);
		qaws_surface* sph = sgf_revolve(prof);
		qaws_surface_curve *c1 = NULL, *c2 = NULL, *c3 = NULL, *t[3] = { 0 };
		qaws_surface_curve const *p1[1], *p2[1], *p3[1], *pt[3];
		qaws_surface_path P1, P2, P3, PT, P12[2];
		qaws_vec2 tri[3];
		static char const* const names[4] = { "a lens of two geodesic circles", "a geodesic triangle less a circle",
			"three geodesic circles' union", "their xor (even-odd)" };
		qaws_surface_curve_create_geodesic_arc(sph, sgf_uv(0.36, 0.52), (qaws_scalar)0.55, 0, (qaws_scalar)(2 * pi), 0, &c1);
		qaws_surface_curve_create_geodesic_arc(sph, sgf_uv(0.44, 0.48), (qaws_scalar)0.5, 0, (qaws_scalar)(2 * pi), 0, &c2);
		qaws_surface_curve_create_geodesic_arc(sph, sgf_uv(0.40, 0.62), (qaws_scalar)0.45, 0, (qaws_scalar)(2 * pi), 0, &c3);
		tri[0] = sgf_uv(0.26, 0.4); tri[1] = sgf_uv(0.52, 0.36); tri[2] = sgf_uv(0.4, 0.75);
		for (i = 0; i < 3; i++) { qaws_surface_curve_create_geodesic(sph, tri[i], tri[(i + 1) % 3], 0, &t[i]); pt[i] = t[i]; }
		p1[0] = c1; p2[0] = c2; p3[0] = c3;
		P1.curves = p1; P1.curve_count = 1; P1.closed = 1;
		P2.curves = p2; P2.curve_count = 1; P2.closed = 1;
		P3.curves = p3; P3.curve_count = 1; P3.closed = 1;
		PT.curves = pt; PT.curve_count = 3; PT.closed = 1;
		P12[0] = P1; P12[1] = P2;
		for (k = 0; k < 4; k++)
		{
			qaws_surface_region* r = NULL;
			qaws_scalar area = 0;
			viewport v3;
			sgf_view view;
			v3.x0 = 24 + k * (pw + 12); v3.y0 = 104 + 2 * (ph + 44); v3.w = pw; v3.h = ph;
			sgf_view_init(&view, 3.0 - 2 * pi * 0.1, 0.25, 105, v3.x0 + pw / 2, v3.y0 + ph / 2 + 10);
			if (k == 0) qaws_surface_clip(sph, QAWS_CLIP_INTERSECTION, QAWS_FILL_NON_ZERO, &P1, 1, &P2, 1, &r);
			else if (k == 1) qaws_surface_clip(sph, QAWS_CLIP_DIFFERENCE, QAWS_FILL_NON_ZERO, &PT, 1, &P3, 1, &r);
			else if (k == 2) qaws_surface_clip(sph, QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO, P12, 2, &P3, 1, &r);
			else
			{
				qaws_surface_path all[3];
				all[0] = P1; all[1] = P2; all[2] = P3;
				qaws_surface_clip(sph, QAWS_CLIP_NONE, QAWS_FILL_EVEN_ODD, all, 3, NULL, 0, &r);
			}
			qaws_surface_region_get_area(r, &area);
			sprintf(label, "%s: area %.4f", names[k], (double)area);
			svg_panel(&s, &v3, label);
			sgf_mesh(&s, sph, &view, 12, 1);
			sbf_fill(&s, r, sph, &view, NULL, 1, "#ff8182", "#cf222e");
			qaws_surface_region_destroy(r);
		}
		qaws_surface_curve_destroy(c1); qaws_surface_curve_destroy(c2); qaws_surface_curve_destroy(c3);
		for (i = 0; i < 3; i++) qaws_surface_curve_destroy(t[i]);
		qaws_surface_destroy(sph);
		qaws_curve_destroy(prof);
	}
	svg_close(&s);
	printf("figure 17: Booleans on surfaces, %.1f s\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
}
