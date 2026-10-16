/* Figure 30: mitered skeletons (qaws_surface_region_mitered_skeleton)
   against medial axes, an L on figure 22's wavy patch and a dart of
   shortest paths on figure 27's bump; geodesic triangles certified
   (qaws_surface_geodesics_certify) with one's area bounded
   (qaws_surface_geodesic_polygon_area); shortest paths on a torus
   seeded by the distance field (qaws_surface_curve_create_shortest_path
   with no obstacles) against the geodesics the chords lead to. Below, the
   same in the parameter planes. Uses figures 16 to 22 and 27's helpers. */

/* the geodesic segment along a shortest path from a to b (its curve kept) */
static qaws_geodesic_segment s30_seg(qaws_surface* f, qaws_vec2 a, qaws_vec2 b, qaws_surface_curve** out)
{
	qaws_geodesic_segment g;
	qaws_eval_result_2d e;
	qaws_range r;
	memset(&g, 0, sizeof(g));
	*out = NULL;
	if (qaws_surface_curve_create_shortest_path(f, NULL, a, b, 0, out) != QAWS_STATUS_OK)
		return g;
	r = qaws_curve_get_parameter_range(qaws_surface_curve_get_uv(*out));
	qaws_curve_evaluate_2d(&e, qaws_surface_curve_get_uv(*out), r.min_value, QAWS_EVAL_FLAG_D1);
	g.start = a; g.direction = e.d1;
	qaws_surface_curve_get_length(*out, &g.length);
	return g;
}

static void demo_mitered_certified(void)
{
	clock_t t0 = clock();
	svg s;
	double const pw = 300, ph = 260, pi = 3.14159265358979323846;
	unsigned int i, k;
	char label[260];
	if (!svg_open(&s, "showcase/b2d30_mitered_certified.svg", (int)(48 + 4 * (pw + 12)), (int)(110 + 2 * (ph + 44)),
		"Mitered skeletons, certified arrangements, shortest paths on a torus",
		"red: mitered skeletons (edges moved as their geodesics), purple: medial axes; certified crossings and area bounds; seeded shortest paths (red) against the chords' geodesics (grey)"))
		return;
	for (k = 0; k < 4; k++)
	{
		viewport v3, vuv;
		sgf_view view;
		qaws_surface* surf = NULL;
		qaws_curve* prof = NULL;
		qaws_surface_region* r = NULL;
		qaws_surface_medial_axis *mm = NULL, *ma = NULL;
		qaws_surface_curve* cv[12] = { 0 };
		qaws_curve* keep[8] = { 0 };
		unsigned int ncv = 0, mesh_lines = 12;
		double pu = 0;
		char const* name = "";
		v3.x0 = 24 + k * (pw + 12); v3.y0 = 104; v3.w = pw; v3.h = ph;
		vuv = v3; vuv.y0 = 104 + ph + 44;
		vuv.xmin = -0.1; vuv.xmax = 1.1; vuv.ymin = 0.5 - 0.6 * ph / pw; vuv.ymax = 0.5 + 0.6 * ph / pw;
		label[0] = 0;
		if (k == 0)
		{
			/* an L drawn in (u, v) on the wavy patch */
			static double const z[16] = { 0, 0.3, -0.2, 0, 0.35, 0.8, 0.5, 0.1, -0.3, 0.6, 0.9, 0.35, 0, -0.15, 0.3, 0 };
			static double const L[6][2] = { { 0.15, 0.15 }, { 0.85, 0.15 }, { 0.85, 0.42 }, { 0.42, 0.42 }, { 0.42, 0.85 }, { 0.15, 0.85 } };
			qaws_vec3 cp[16];
			qaws_surface_bezier_desc bd;
			qaws_surface_curve const* pc[6];
			qaws_surface_path Q;
			for (i = 0; i < 16; i++) { cp[i].x = (qaws_scalar)((i / 4) / 3.0 * 2); cp[i].y = (qaws_scalar)((i % 4) / 3.0 * 2); cp[i].z = (qaws_scalar)z[i]; }
			memset(&bd, 0, sizeof(bd));
			bd.u_degree = 3; bd.v_degree = 3; bd.control_points = cp; bd.u_point_count = 4; bd.v_point_count = 4;
			qaws_surface_create_bezier(&surf, &bd);
			for (i = 0; i < 6; i++) { cv[i] = ssf_line(surf, L[i][0], L[i][1], L[(i + 1) % 6][0], L[(i + 1) % 6][1], &keep[i]); pc[i] = cv[i]; }
			ncv = 6;
			ssf_path(&Q, pc, 6);
			qaws_surface_clip(surf, QAWS_CLIP_NONE, QAWS_FILL_NON_ZERO, &Q, 1, NULL, 0, &r);
			sgf_view_init(&view, 0.7, 0.75, 100, v3.x0 + pw / 2, v3.y0 + ph / 2 + 10);
			view.o[0] = 1; view.o[1] = 1; view.o[2] = 0.3;
			mesh_lines = 10;
			name = "an L on a wavy patch: mitered against medial";
		}
		else if (k == 1)
		{
			/* a dart of shortest paths on the bump: a reflex corner */
			qaws_vec2 v[4];
			qaws_surface_curve const* pc[4];
			qaws_surface_path Q;
			qaws_geodesic_segment g;
			surf = s27_bump();
			v[0] = sgf_uv(0.12, 0.15); v[1] = sgf_uv(0.88, 0.3); v[2] = sgf_uv(0.55, 0.5); v[3] = sgf_uv(0.35, 0.88);
			for (i = 0; i < 4; i++) { g = s30_seg(surf, v[i], v[(i + 1) % 4], &cv[i]); pc[i] = cv[i]; }
			(void)g;
			ncv = 4;
			ssf_path(&Q, pc, 4);
			qaws_surface_clip(surf, QAWS_CLIP_NONE, QAWS_FILL_NON_ZERO, &Q, 1, NULL, 0, &r);
			sgf_view_init(&view, -0.7, 0.75, 110, v3.x0 + pw / 2, v3.y0 + ph / 2 + 10);
			name = "a dart of shortest paths on the bump";
		}
		else if (k == 2)
		{
			/* two geodesic triangles, their crossings certified; one's area bounded */
			static double const T[2][3][2] = { { { 0.15, 0.2 }, { 0.75, 0.3 }, { 0.4, 0.8 } }, { { 0.3, 0.45 }, { 0.85, 0.55 }, { 0.6, 0.15 } } };
			qaws_geodesic_segment g[6];
			static unsigned int const corner[6] = { 1, 2, 0, 4, 5, 3 };   /* each triangle closed: an edge ends where the next starts */
			qaws_surface_curve const* pc[3];
			qaws_surface_path Q;
			qaws_scalar lo = 0, hi = 0, A = 0;
			unsigned int n = 0, j;
			int sc, sa;
			surf = s27_bump();
			for (i = 0; i < 6; i++)
				g[i] = s30_seg(surf, sgf_uv(T[i / 3][i % 3][0], T[i / 3][i % 3][1]), sgf_uv(T[i / 3][(i + 1) % 3][0], T[i / 3][(i + 1) % 3][1]), &cv[i]);
			ncv = 6;
			for (i = 0; i < 3; i++) pc[i] = cv[i];
			ssf_path(&Q, pc, 3);
			qaws_surface_clip(surf, QAWS_CLIP_NONE, QAWS_FILL_NON_ZERO, &Q, 1, NULL, 0, &r);
			if (r) qaws_surface_region_get_area(r, &A);
			sc = (int)qaws_surface_geodesics_certify(surf, g, 6, corner, 0, &n);
			sa = (int)qaws_surface_geodesic_polygon_area(surf, g, 3, 0, 512, &lo, &hi);
			sgf_view_init(&view, -0.7, 0.75, 110, v3.x0 + pw / 2, v3.y0 + ph / 2 + 10);
			name = "geodesic triangles certified on the bump";
			sprintf(label, "%u crossings certified%s; area [%.3f, %.3f]%s", n, sc ? " (failed)" : "", (double)lo, (double)hi, sa ? " (failed)" : "");
			(void)j;
		}
		else
		{
			/* a torus: shortest paths seeded by the distance field */
			static double const P[4][4] = { { 0.2, 0, 0.8, 0 }, { 0.2, 0, 0.8, 0.05 }, { 0.1, 0.3, 0.6, 0.7 }, { 0.05, 0.5, 0.55, 0.5 } };
			prof = sgf_arc(1.0, 0, 0.35, 0, 2 * pi);
			surf = sgf_revolve(prof);
			pu = 1;
			for (i = 0; i < 4; i++)
			{
				qaws_surface_curve_create_shortest_path(surf, NULL, sgf_uv(P[i][0], P[i][1]), sgf_uv(P[i][2], P[i][3]), 0, &cv[2 * i]);
				qaws_surface_curve_create_geodesic(surf, sgf_uv(P[i][0], P[i][1]), sgf_uv(P[i][2], P[i][3]), 0, &cv[2 * i + 1]);
			}
			ncv = 8;
			sgf_view_init(&view, -pi / 2 + 0.35, 0.75, 105, v3.x0 + pw / 2, v3.y0 + ph / 2 + 4);
			mesh_lines = 16;
			name = "a torus: shortest paths seeded";
			for (i = 0; i < 4; i++)
			{
				qaws_scalar a = 0, b = 0;
				char one[40];
				if (cv[2 * i]) qaws_surface_curve_get_length(cv[2 * i], &a);
				if (cv[2 * i + 1]) qaws_surface_curve_get_length(cv[2 * i + 1], &b);
				sprintf(one, "%s%.2f/%.2f", i ? ", " : "seeded/chord ", (double)a, (double)b);
				strcat(label, one);
			}
		}
		if (k < 2 && r)
		{
			qaws_surface_medial_desc md;
			memset(&md, 0, sizeof(md));
			md.resolution = (qaws_scalar)0.04;
			qaws_surface_region_medial_axis(r, &md, &ma);
			if (qaws_surface_region_mitered_skeleton(r, &md, &mm) == QAWS_STATUS_OK)
				sprintf(label, "mitered: %u segments, %u junctions (medial %u)", qaws_surface_medial_axis_get_segment_count(mm),
					qaws_surface_medial_axis_get_junction_count(mm), qaws_surface_medial_axis_get_segment_count(ma));
		}
		svg_panel(&s, &v3, name);
		svg_panel(&s, &vuv, label);
		sgf_mesh(&s, surf, &view, mesh_lines, pu);
		ssf_box(&s, &vuv, surf);
		if (k < 2 && r)
		{
			ssf_fill3d(&s, r, surf, &view, "#54aeff");
			sbf_fill(&s, r, surf, NULL, &vuv, pu, "#54aeff", "none");
			sof_outline(&s, r, surf, &view, &vuv, pu, 0, "#0969da", 1.4);
			if (ma) smf_axis(&s, ma, surf, &view, &vuv, pu, "#8250df");
			if (mm) smf_axis(&s, mm, surf, &view, &vuv, pu, "#cf222e");
		}
		else if (k == 2)
		{
			if (r)
			{
				ssf_fill3d(&s, r, surf, &view, "#54aeff");
				sbf_fill(&s, r, surf, NULL, &vuv, pu, "#54aeff", "none");
			}
			for (i = 0; i < 6; i++)
				if (cv[i]) ssf_curve(&s, cv[i], surf, &view, &vuv, pu, 0, i < 3 ? "#0969da" : "#cf222e");
			for (i = 0; i < 3; i++)
			{
				unsigned int j, nx = 0, q;
				qaws_surface_crossing x[4];
				for (j = 3; j < 6; j++)
					if (cv[i] && cv[j] && qaws_surface_curve_intersect(cv[i], cv[j], x, 4, &nx) == QAWS_STATUS_OK)
						for (q = 0; q < nx && q < 4; q++)
						{
							qaws_surface_eval_result e;
							double sx, sy;
							qaws_surface_evaluate(&e, surf, x[q].point.x, x[q].point.y, QAWS_SURFACE_EVAL_POSITION);
							sgf_project(&view, e.position, &sx, &sy, NULL);
							fprintf(s.f, "<circle cx=\"%.2f\" cy=\"%.2f\" r=\"3.5\" fill=\"#1b1f24\"/>\n", sx, sy);
							fprintf(s.f, "<circle cx=\"%.2f\" cy=\"%.2f\" r=\"3\" fill=\"#1b1f24\"/>\n", vx(&vuv, x[q].point.x), vy(&vuv, x[q].point.y));
						}
			}
		}
		else
			for (i = 0; i < 8; i++)
				if (cv[i]) ssf_curve(&s, cv[i], surf, &view, &vuv, pu, 1, i % 2 ? "#8c959f" : "#cf222e");
		qaws_surface_medial_axis_destroy(mm);
		qaws_surface_medial_axis_destroy(ma);
		qaws_surface_region_destroy(r);
		for (i = 0; i < ncv; i++) qaws_surface_curve_destroy(cv[i]);
		for (i = 0; i < 8; i++) qaws_curve_destroy(keep[i]);
		qaws_surface_destroy(surf);
		qaws_curve_destroy(prof);
		printf("  panel %u: %s\n", k, label);
	}
	svg_close(&s);
	printf("figure 30: mitered skeletons, certified arrangements, torus paths, %.1f s\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
}
