/* Figure 22: medial axes of regions on surfaces
   (qaws_surface_region_medial_axis). The points at the same geodesic
   distance from two nearest boundary points, with the largest inscribed
   geodesic circle: a geodesic triangle and a band on a sphere, a geodesic
   quadrilateral on a wavy patch, and a square less a circle across a
   torus's seam. Below, the same in the parameter planes. Uses figures 16
   to 19's helpers. */

/* the axis's segments in 3D and in (u, v), junctions as dots */
static void smf_axis(svg* s, qaws_surface_medial_axis const* m, qaws_surface const* surf, sgf_view const* view, viewport const* vuv,
	double pu, char const* color)
{
	unsigned int i, n = qaws_surface_medial_axis_get_segment_count(m);
	for (i = 0; i < n; i++)
	{
		qaws_vec2 p[2];
		qaws_scalar r[2];
		double uv[2 * 9], xy[4];
		unsigned int k;
		qaws_surface_medial_axis_get_segment(m, i, p, r);
		for (k = 0; k <= 8; k++)
		{
			uv[2 * k] = p[0].x + (p[1].x - p[0].x) * k / 8.0;
			uv[2 * k + 1] = p[0].y + (p[1].y - p[0].y) * k / 8.0;
		}
		sgf_curve3d(s, surf, view, uv, 9, pu, color, 1.6, 0);
		xy[0] = vx(vuv, p[0].x); xy[1] = vy(vuv, p[0].y); xy[2] = vx(vuv, p[1].x); xy[3] = vy(vuv, p[1].y);
		svg_polyline(s, xy, 2, color, 1.3, 1.0, 0);
	}
	for (i = 0; i < qaws_surface_medial_axis_get_junction_count(m); i++)
	{
		qaws_vec2 p;
		qaws_scalar r;
		qaws_surface_eval_result e;
		double sx, sy;
		qaws_surface_medial_axis_get_junction(m, i, &p, &r);
		qaws_surface_evaluate(&e, surf, p.x, p.y, QAWS_SURFACE_EVAL_POSITION);
		if (sgf_facing(surf, view, p.x, p.y, pu))
		{
			sgf_project(view, e.position, &sx, &sy, NULL);
			fprintf(s->f, "<circle cx=\"%.2f\" cy=\"%.2f\" r=\"3\" fill=\"%s\"/>\n", sx, sy, color);
		}
		fprintf(s->f, "<circle cx=\"%.2f\" cy=\"%.2f\" r=\"2.5\" fill=\"%s\"/>\n", vx(vuv, p.x), vy(vuv, p.y), color);
	}
}

/* the largest inscribed geodesic circle */
static void smf_inscribed(svg* s, qaws_surface_medial_axis const* m, qaws_surface const* surf, sgf_view const* view, viewport const* vuv, double pu)
{
	qaws_vec2 c;
	qaws_scalar r;
	qaws_surface_curve* circ = NULL;
	if (qaws_surface_medial_axis_get_max_inscribed(m, &c, &r) != QAWS_STATUS_OK || !(r > 0))
		return;
	if (qaws_surface_curve_create_geodesic_arc(surf, c, (qaws_scalar)(0.999 * r), 0, (qaws_scalar)(2 * 3.14159265358979323846), 0, &circ) == QAWS_STATUS_OK)
		ssf_curve(s, circ, surf, view, vuv, pu, 0, "#1a7f37");
	qaws_surface_curve_destroy(circ);
}

static void demo_surface_medial(void)
{
	clock_t t0 = clock();
	svg s;
	double const pw = 300, ph = 260, pi = 3.14159265358979323846;
	unsigned int i, k;
	char label[220];
	if (!svg_open(&s, "showcase/b2d22_surface_medial.svg", (int)(48 + 4 * (pw + 12)), (int)(110 + 2 * (ph + 44)),
		"Medial axes of regions on surfaces",
		"qaws_surface_region_medial_axis: points equally far (geodesically) from two nearest boundary points; junctions; the largest inscribed geodesic circle (green)"))
		return;
	for (k = 0; k < 4; k++)
	{
		viewport v3, vuv;
		sgf_view view;
		qaws_surface* surf = NULL;
		qaws_curve* prof = NULL;
		qaws_surface_region* r = NULL;
		qaws_surface_medial_axis* m = NULL;
		qaws_surface_medial_desc md;
		qaws_surface_curve* cv[8] = { 0 };
		qaws_curve* keep[8] = { 0 };
		unsigned int ncv = 0, mesh_lines = 12;
		double pu = 0;
		char const* name = "";
		v3.x0 = 24 + k * (pw + 12); v3.y0 = 104; v3.w = pw; v3.h = ph;
		vuv = v3; vuv.y0 = 104 + ph + 44;
		vuv.xmin = -0.1; vuv.xmax = 1.1; vuv.ymin = 0.5 - 0.6 * ph / pw; vuv.ymax = 0.5 + 0.6 * ph / pw;
		memset(&md, 0, sizeof(md));
		md.resolution = (qaws_scalar)0.04;
		if (k == 0 || k == 1)
		{
			prof = sgf_arc(0, 0, 1, -pi / 2, pi / 2);
			surf = sgf_revolve(prof);
			pu = 1;
			sgf_view_init(&view, -pi / 2 - 2 * pi * 0.38 + 0.2, k == 1 ? 0.45 : 0.25, 108, v3.x0 + pw / 2, v3.y0 + ph / 2 + 8);
			if (k == 0)
			{
				qaws_vec2 v[3];
				qaws_surface_curve const* pc[3];
				qaws_surface_path T;
				v[0] = sgf_uv(0.28, 0.38); v[1] = sgf_uv(0.5, 0.42); v[2] = sgf_uv(0.37, 0.72);
				for (i = 0; i < 3; i++) { qaws_surface_curve_create_geodesic(surf, v[i], v[(i + 1) % 3], 0, &cv[i]); pc[i] = cv[i]; }
				ncv = 3;
				ssf_path(&T, pc, 3);
				qaws_surface_clip(surf, QAWS_CLIP_NONE, QAWS_FILL_NON_ZERO, &T, 1, NULL, 0, &r);
				name = "a geodesic triangle on a sphere";
			}
			else
			{
				/* a band with a bite: latitudes less a circle */
				qaws_surface_curve const *pb[1], *pc2[1], *pcirc[1];
				qaws_surface_path B[2], C;
				cv[0] = ssf_line(surf, 0, 0.38, 1, 0.38, &keep[0]);
				cv[1] = ssf_line(surf, 1, 0.62, 0, 0.62, &keep[1]);
				qaws_surface_curve_create_geodesic_arc(surf, sgf_uv(0.38, 0.62), (qaws_scalar)0.3, 0, (qaws_scalar)(2 * pi), 0, &cv[2]);
				ncv = 3;
				pb[0] = cv[0]; pc2[0] = cv[1]; pcirc[0] = cv[2];
				ssf_path(&B[0], pb, 1); ssf_path(&B[1], pc2, 1); ssf_path(&C, pcirc, 1);
				qaws_surface_clip(surf, QAWS_CLIP_DIFFERENCE, QAWS_FILL_NON_ZERO, B, 2, &C, 1, &r);
				name = "a band round a sphere, less a circle";
			}
		}
		else if (k == 2)
		{
			static double const z[16] = { 0, 0.3, -0.2, 0, 0.35, 0.8, 0.5, 0.1, -0.3, 0.6, 0.9, 0.35, 0, -0.15, 0.3, 0 };
			qaws_vec3 cp[16];
			qaws_surface_bezier_desc bd;
			qaws_vec2 c[4];
			qaws_surface_curve const* pc[4];
			qaws_surface_path Q;
			for (i = 0; i < 16; i++) { cp[i].x = (qaws_scalar)((i / 4) / 3.0 * 2); cp[i].y = (qaws_scalar)((i % 4) / 3.0 * 2); cp[i].z = (qaws_scalar)z[i]; }
			memset(&bd, 0, sizeof(bd));
			bd.u_degree = 3; bd.v_degree = 3; bd.control_points = cp; bd.u_point_count = 4; bd.v_point_count = 4;
			qaws_surface_create_bezier(&surf, &bd);
			c[0] = sgf_uv(0.15, 0.2); c[1] = sgf_uv(0.85, 0.25); c[2] = sgf_uv(0.75, 0.7); c[3] = sgf_uv(0.2, 0.8);
			for (i = 0; i < 4; i++) { qaws_surface_curve_create_geodesic(surf, c[i], c[(i + 1) % 4], 0, &cv[i]); pc[i] = cv[i]; }
			ncv = 4;
			ssf_path(&Q, pc, 4);
			qaws_surface_clip(surf, QAWS_CLIP_NONE, QAWS_FILL_NON_ZERO, &Q, 1, NULL, 0, &r);
			sgf_view_init(&view, 0.7, 0.75, 100, v3.x0 + pw / 2, v3.y0 + ph / 2 + 10);
			view.o[0] = 1; view.o[1] = 1; view.o[2] = 0.3;
			mesh_lines = 10;
			name = "a geodesic quadrilateral on a wavy patch";
		}
		else
		{
			/* a torus: a square drawn in (u, v) across the u seam less a circle */
			qaws_surface_curve const *pq[4], *pc[1];
			qaws_surface_path Q, C;
			prof = sgf_arc(1.0, 0, 0.42, 0, 2 * pi);
			surf = sgf_revolve(prof);
			pu = 1;
			cv[0] = ssf_line(surf, -0.1, 0.05, 0.12, 0.05, &keep[0]);
			cv[1] = ssf_line(surf, 0.12, 0.05, 0.12, 0.3, &keep[1]);
			cv[2] = ssf_line(surf, 0.12, 0.3, -0.1, 0.3, &keep[2]);
			cv[3] = ssf_line(surf, -0.1, 0.3, -0.1, 0.05, &keep[3]);
			for (i = 0; i < 4; i++) pq[i] = cv[i];
			qaws_surface_curve_create_geodesic_arc(surf, sgf_uv(0.02, 0.18), (qaws_scalar)0.12, 0, (qaws_scalar)(2 * pi), 0, &cv[4]);
			ncv = 5;
			pc[0] = cv[4];
			ssf_path(&Q, pq, 4); ssf_path(&C, pc, 1);
			qaws_surface_clip(surf, QAWS_CLIP_DIFFERENCE, QAWS_FILL_NON_ZERO, &Q, 1, &C, 1, &r);
			sgf_view_init(&view, -pi / 2 + 0.35, 0.75, 105, v3.x0 + pw / 2, v3.y0 + ph / 2 + 4);
			mesh_lines = 16;
			name = "a torus: a square less a circle, across the seam";
		}
		svg_panel(&s, &v3, name);
		sgf_mesh(&s, surf, &view, mesh_lines, pu);
		ssf_box(&s, &vuv, surf);
		if (r && qaws_surface_region_medial_axis(r, &md, &m) == QAWS_STATUS_OK)
		{
			qaws_vec2 c;
			qaws_scalar rr = 0;
			qaws_surface_medial_axis_get_max_inscribed(m, &c, &rr);
			sprintf(label, "(u, v): %u segments, %u junctions, largest circle %.4f", qaws_surface_medial_axis_get_segment_count(m),
				qaws_surface_medial_axis_get_junction_count(m), (double)rr);
			svg_panel(&s, &vuv, label);
			ssf_fill3d(&s, r, surf, &view, "#54aeff");
			sbf_fill(&s, r, surf, NULL, &vuv, pu, "#54aeff", "none");
			sof_outline(&s, r, surf, &view, &vuv, pu, 0, "#0969da", 1.4);
			smf_inscribed(&s, m, surf, &view, &vuv, pu);
			smf_axis(&s, m, surf, &view, &vuv, pu, "#cf222e");
		}
		qaws_surface_medial_axis_destroy(m);
		qaws_surface_region_destroy(r);
		for (i = 0; i < ncv; i++) qaws_surface_curve_destroy(cv[i]);
		for (i = 0; i < 8; i++) qaws_curve_destroy(keep[i]);
		qaws_surface_destroy(surf);
		qaws_curve_destroy(prof);
	}
	svg_close(&s);
	printf("figure 22: medial axes on surfaces, %.1f s\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
}
