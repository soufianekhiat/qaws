/* Figure 20: sampling and meshing regions on surfaces in their metric
   (qaws_surface_region_sample, qaws_surface_region_triangulate). On a
   sphere, a geodesic circle across the seam: Sobol points uniform in
   surface area, blue noise spaced on the surface, and a band sampled with
   the density 1 + z/R; the circle's parameter plane, where the uniform
   points crowd toward the pole. Below: the circle meshed to a surface edge
   length and its mesh in (u, v), larger where the metric shrinks; a wavy
   patch's square less a circle, and a torus's circle across both seams.
   Uses figures 16 to 19's helpers. */

static void ssm_points(svg* s, qaws_surface const* surf, sgf_view const* view, qaws_vec2 const* uv, qaws_vec3 const* p, unsigned int n,
	viewport const* vuv, char const* color, double radius)
{
	unsigned int i;
	for (i = 0; i < n; i++)
	{
		double sx, sy;
		if (view)
		{
			if (!sgf_facing(surf, view, uv[i].x, uv[i].y, 0))
				continue;
			sgf_project(view, p[i], &sx, &sy, NULL);
		}
		else
		{
			sx = vx(vuv, uv[i].x);
			sy = vy(vuv, uv[i].y);
		}
		fprintf(s->f, "<circle cx=\"%.2f\" cy=\"%.2f\" r=\"%.2f\" fill=\"%s\"/>\n", sx, sy, radius, color);
	}
}

/* the mesh's edges in space (front-facing triangles) or in (u, v) */
static void ssm_mesh(svg* s, qaws_surface const* surf, sgf_view const* view, viewport const* vuv, qaws_mesh_2d const* m, char const* color)
{
	unsigned int t, k, nt = qaws_mesh_2d_get_triangle_count(m);
	qaws_range ru = qaws_surface_get_u_range(surf), rv = qaws_surface_get_v_range(surf);
	fprintf(s->f, "<path fill=\"%s\" fill-opacity=\"0.18\" stroke=\"%s\" stroke-width=\"0.45\" stroke-linejoin=\"round\" d=\"", color, color);
	for (t = 0; t < nt; t++)
	{
		unsigned int v[3];
		double xy[6], cu = 0, cv = 0;
		qaws_mesh_2d_get_triangle(m, t, v);
		for (k = 0; k < 3; k++)
		{
			qaws_vec2 w = qaws_mesh_2d_get_vertex(m, v[k]);
			double uu = w.x, vv = w.y;
			cu += uu / 3; cv += vv / 3;
			if (view)
			{
				qaws_surface_eval_result e;
				if (uu < ru.min_value) uu = ru.min_value;
				if (uu > ru.max_value) uu = ru.max_value;
				if (vv < rv.min_value) vv = rv.min_value;
				if (vv > rv.max_value) vv = rv.max_value;
				qaws_surface_evaluate(&e, surf, (qaws_scalar)uu, (qaws_scalar)vv, QAWS_SURFACE_EVAL_POSITION);
				sgf_project(view, e.position, &xy[2 * k], &xy[2 * k + 1], NULL);
			}
			else
			{
				xy[2 * k] = vx(vuv, uu);
				xy[2 * k + 1] = vy(vuv, vv);
			}
		}
		if (view && !sgf_facing(surf, view, cu, cv, 0))
			continue;
		fprintf(s->f, "M%.2f,%.2f L%.2f,%.2f L%.2f,%.2f Z ", xy[0], xy[1], xy[2], xy[3], xy[4], xy[5]);
	}
	fprintf(s->f, "\"/>\n");
}

static qaws_scalar ssm_density(void* user, qaws_scalar u, qaws_scalar v)
{
	(void)user; (void)u;
	return (qaws_scalar)(1 + sin(3.14159265358979323846 * (double)v - 3.14159265358979323846 / 2));
}

static void demo_surface_sample(void)
{
	clock_t t0 = clock();
	svg s;
	double const pw = 300, ph = 260, pi = 3.14159265358979323846;
	unsigned int i, k;
	char label[220];
	if (!svg_open(&s, "showcase/b2d20_surface_sample.svg", (int)(48 + 4 * (pw + 12)), (int)(110 + 2 * (ph + 44)),
		"Sampling and meshing regions on surfaces",
		"qaws_surface_region_sample and _triangulate: points uniform in surface area (or a density), blue noise spaced on the surface, meshes sized by the metric"))
		return;
	/* the sphere */
	{
		qaws_curve* prof = sgf_arc(0, 0, 1, -pi / 2, pi / 2);
		qaws_surface* sph = sgf_revolve(prof);
		qaws_surface_curve *c = NULL, *lat[2] = { 0 };
		qaws_curve* keep[2] = { 0 };
		qaws_surface_curve const *pc[1], *pl[2];
		qaws_surface_path P, B[2];
		qaws_surface_region *rc = NULL, *rb = NULL;
		qaws_surface_sample_desc sd;
		unsigned int n = 1500;
		qaws_vec2* w = (qaws_vec2*)malloc(sizeof(qaws_vec2) * n);
		qaws_vec3* p = (qaws_vec3*)malloc(sizeof(qaws_vec3) * n);
		qaws_surface_curve_create_geodesic_arc(sph, sgf_uv(0.02, 0.62), (qaws_scalar)0.62, 0, (qaws_scalar)(2 * pi), 0, &c);
		pc[0] = c;
		ssf_path(&P, pc, 1);
		qaws_surface_clip(sph, QAWS_CLIP_NONE, QAWS_FILL_NON_ZERO, &P, 1, NULL, 0, &rc);
		lat[0] = ssf_line(sph, 0, 0.3, 1, 0.3, &keep[0]);
		lat[1] = ssf_line(sph, 1, 0.75, 0, 0.75, &keep[1]);
		pl[0] = lat[0]; pl[1] = lat[1];
		ssf_path(&B[0], pl, 1); ssf_path(&B[1], pl + 1, 1);
		qaws_surface_clip(sph, QAWS_CLIP_NONE, QAWS_FILL_NON_ZERO, B, 2, NULL, 0, &rb);
		for (k = 0; k < 4; k++)
		{
			viewport v3, vuv;
			sgf_view view;
			unsigned int m = k == 1 ? 600 : n;
			v3.x0 = 24 + k * (pw + 12); v3.y0 = 104; v3.w = pw; v3.h = ph;
			vuv = v3;
			vuv.xmin = -0.1; vuv.xmax = 1.1; vuv.ymin = 0.5 - 0.6 * ph / pw; vuv.ymax = 0.5 + 0.6 * ph / pw;
			sgf_view_init(&view, -pi / 2 + 0.25, 0.45, 108, v3.x0 + pw / 2, v3.y0 + ph / 2 + 8);
			memset(&sd, 0, sizeof(sd));
			sd.seed = 3;
			if (k == 0) { sd.kind = QAWS_REGION_SAMPLE_SOBOL; svg_panel(&s, &v3, "Sobol points, uniform in surface area"); }
			else if (k == 1) { sd.kind = QAWS_REGION_SAMPLE_BLUE_NOISE; svg_panel(&s, &v3, "blue noise, spaced on the surface"); }
			else if (k == 2) { sd.kind = QAWS_REGION_SAMPLE_SOBOL; sd.density = ssm_density; svg_panel(&s, &v3, "a band with the density 1 + z/R"); }
			else { sd.kind = QAWS_REGION_SAMPLE_SOBOL; svg_panel(&s, &v3, "(u, v): the uniform points crowd to the pole"); }
			if (k < 3)
			{
				sgf_mesh(&s, sph, &view, 12, 1);
				if (qaws_surface_region_sample(k == 2 ? rb : rc, &sd, m, w, p) == QAWS_STATUS_OK)
					ssm_points(&s, sph, &view, w, p, m, NULL, k == 1 ? "#cf222e" : "#0969da", 1.2);
				sof_outline(&s, k == 2 ? rb : rc, sph, &view, NULL, 1, 0, "#57606a", 1.2);
			}
			else
			{
				ssf_box(&s, &vuv, sph);
				if (qaws_surface_region_sample(rc, &sd, m, w, p) == QAWS_STATUS_OK)
					ssm_points(&s, sph, NULL, w, p, m, &vuv, "#0969da", 1.0);
				sbf_fill(&s, rc, sph, NULL, &vuv, 1, "none", "#57606a");
			}
		}
		/* meshes: the circle to 0.06, in space and in (u, v) */
		{
			qaws_surface_mesh_desc md;
			qaws_mesh_2d* mesh = NULL;
			viewport v3, vuv;
			sgf_view view;
			v3.x0 = 24; v3.y0 = 104 + ph + 44; v3.w = pw; v3.h = ph;
			vuv = v3; vuv.x0 = v3.x0 + pw + 12;
			vuv.xmin = -0.1; vuv.xmax = 1.1; vuv.ymin = 0.5 - 0.6 * ph / pw; vuv.ymax = 0.5 + 0.6 * ph / pw;
			sgf_view_init(&view, -pi / 2 + 0.25, 0.45, 108, v3.x0 + pw / 2, v3.y0 + ph / 2 + 8);
			memset(&md, 0, sizeof(md));
			md.edge_length = (qaws_scalar)0.06;
			md.remesh_iterations = 3;
			if (qaws_surface_region_triangulate(rc, &md, &mesh) == QAWS_STATUS_OK)
			{
				sprintf(label, "meshed to edges of 0.06: %u triangles", qaws_mesh_2d_get_triangle_count(mesh));
				svg_panel(&s, &v3, label);
				svg_panel(&s, &vuv, "(u, v): larger where the metric shrinks");
				sgf_mesh(&s, sph, &view, 12, 1);
				ssm_mesh(&s, sph, &view, NULL, mesh, "#0969da");
				ssf_box(&s, &vuv, sph);
				ssm_mesh(&s, sph, NULL, &vuv, mesh, "#0969da");
			}
			qaws_mesh_2d_destroy(mesh);
		}
		free(w); free(p);
		qaws_surface_region_destroy(rc); qaws_surface_region_destroy(rb);
		qaws_surface_curve_destroy(c);
		for (i = 0; i < 2; i++) { qaws_surface_curve_destroy(lat[i]); qaws_curve_destroy(keep[i]); }
		qaws_surface_destroy(sph);
		qaws_curve_destroy(prof);
	}
	/* a wavy patch: a square less a circle */
	{
		static double const z[16] = { 0, 0.3, -0.2, 0, 0.35, 0.8, 0.5, 0.1, -0.3, 0.6, 0.9, 0.35, 0, -0.15, 0.3, 0 };
		qaws_vec3 cp[16];
		qaws_surface_bezier_desc bd;
		qaws_surface* patch = NULL;
		qaws_surface_curve *sq[4] = { 0 }, *circ = NULL;
		qaws_surface_curve const *pa[4], *pb[1];
		qaws_surface_path A, Bp;
		qaws_surface_region* r = NULL;
		qaws_vec2 corner[4];
		qaws_surface_mesh_desc md;
		qaws_mesh_2d* mesh = NULL;
		viewport v3;
		sgf_view view;
		for (i = 0; i < 16; i++) { cp[i].x = (qaws_scalar)((i / 4) / 3.0 * 2); cp[i].y = (qaws_scalar)((i % 4) / 3.0 * 2); cp[i].z = (qaws_scalar)z[i]; }
		memset(&bd, 0, sizeof(bd));
		bd.u_degree = 3; bd.v_degree = 3; bd.control_points = cp; bd.u_point_count = 4; bd.v_point_count = 4;
		qaws_surface_create_bezier(&patch, &bd);
		corner[0] = sgf_uv(0.15, 0.18); corner[1] = sgf_uv(0.85, 0.15); corner[2] = sgf_uv(0.82, 0.85); corner[3] = sgf_uv(0.18, 0.82);
		for (i = 0; i < 4; i++) { qaws_surface_curve_create_geodesic(patch, corner[i], corner[(i + 1) % 4], 0, &sq[i]); pa[i] = sq[i]; }
		qaws_surface_curve_create_geodesic_arc(patch, sgf_uv(0.5, 0.52), (qaws_scalar)0.4, 0, (qaws_scalar)(2 * pi), 0, &circ);
		pb[0] = circ;
		ssf_path(&A, pa, 4); ssf_path(&Bp, pb, 1);
		qaws_surface_clip(patch, QAWS_CLIP_DIFFERENCE, QAWS_FILL_NON_ZERO, &A, 1, &Bp, 1, &r);
		v3.x0 = 24 + 2 * (pw + 12); v3.y0 = 104 + ph + 44; v3.w = pw; v3.h = ph;
		sgf_view_init(&view, 0.7, 0.75, 100, v3.x0 + pw / 2, v3.y0 + ph / 2 + 10);
		view.o[0] = 1; view.o[1] = 1; view.o[2] = 0.3;
		memset(&md, 0, sizeof(md));
		md.edge_length = (qaws_scalar)0.08;
		md.remesh_iterations = 3;
		if (r && qaws_surface_region_triangulate(r, &md, &mesh) == QAWS_STATUS_OK)
		{
			sprintf(label, "a wavy patch: a square less a circle, %u triangles", qaws_mesh_2d_get_triangle_count(mesh));
			svg_panel(&s, &v3, label);
			sgf_mesh(&s, patch, &view, 10, 0);
			ssm_mesh(&s, patch, &view, NULL, mesh, "#8250df");
		}
		qaws_mesh_2d_destroy(mesh);
		qaws_surface_region_destroy(r);
		for (i = 0; i < 4; i++) qaws_surface_curve_destroy(sq[i]);
		qaws_surface_curve_destroy(circ);
		qaws_surface_destroy(patch);
	}
	/* a torus: a circle across both seams */
	{
		qaws_curve* prof = sgf_arc(1.0, 0, 0.42, 0, 2 * pi);
		qaws_surface* tor = sgf_revolve(prof);
		qaws_surface_curve* c = NULL;
		qaws_surface_curve const* pc[1];
		qaws_surface_path P;
		qaws_surface_region* r = NULL;
		qaws_surface_mesh_desc md;
		qaws_mesh_2d* mesh = NULL;
		viewport v3;
		sgf_view view;
		qaws_surface_curve_create_geodesic_arc(tor, sgf_uv(0.0, 0.0), (qaws_scalar)0.4, 0, (qaws_scalar)(2 * pi), 0, &c);
		pc[0] = c;
		ssf_path(&P, pc, 1);
		qaws_surface_clip(tor, QAWS_CLIP_NONE, QAWS_FILL_NON_ZERO, &P, 1, NULL, 0, &r);
		v3.x0 = 24 + 3 * (pw + 12); v3.y0 = 104 + ph + 44; v3.w = pw; v3.h = ph;
		sgf_view_init(&view, -pi / 2 + 0.35, 0.45, 105, v3.x0 + pw / 2, v3.y0 + ph / 2 + 4);
		memset(&md, 0, sizeof(md));
		md.edge_length = (qaws_scalar)0.045;
		md.remesh_iterations = 3;
		if (r && qaws_surface_region_triangulate(r, &md, &mesh) == QAWS_STATUS_OK)
		{
			sprintf(label, "a torus: a circle across both seams, %u triangles", qaws_mesh_2d_get_triangle_count(mesh));
			svg_panel(&s, &v3, label);
			sgf_mesh(&s, tor, &view, 16, 1);
			ssm_mesh(&s, tor, &view, NULL, mesh, "#bc4c00");
		}
		qaws_mesh_2d_destroy(mesh);
		qaws_surface_region_destroy(r);
		qaws_surface_curve_destroy(c);
		qaws_surface_destroy(tor);
		qaws_curve_destroy(prof);
	}
	svg_close(&s);
	printf("figure 20: sampling and meshing on surfaces, %.1f s\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
}
