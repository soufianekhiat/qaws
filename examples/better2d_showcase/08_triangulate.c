/* Figure 8: meshing regions of curves (qaws_triangulate): constrained
   Delaunay triangulations with the boundary sampled on the true curves,
   holes, Ruppert refinement, a size field, a self-crossing input read with
   the fill rule, isotropic remeshing, constraint curves inside the region,
   and quadrangulation (quad-dominant, then all quads, constrained). */

static void tr_draw(svg* s, viewport const* v, qaws_mesh_2d const* m, char const* fill, char const* line)
{
	unsigned int i, k, nt = qaws_mesh_2d_get_triangle_count(m), nq = qaws_mesh_2d_get_quad_count(m);
	unsigned int nb = qaws_mesh_2d_get_boundary_edge_count(m), nc = qaws_mesh_2d_get_constraint_edge_count(m);
	fprintf(s->f, "<path fill=\"%s\" fill-opacity=\"0.35\" stroke=\"%s\" stroke-width=\"0.5\" stroke-opacity=\"0.85\" stroke-linejoin=\"round\" d=\"", fill, line);
	for (i = 0; i < nq; i++)
	{
		unsigned int q[4];
		qaws_mesh_2d_get_quad(m, i, q);
		for (k = 0; k < 4; k++)
		{
			qaws_vec2 p = qaws_mesh_2d_get_vertex(m, q[k]);
			fprintf(s->f, "%s%.2f,%.2f ", k ? "L" : "M", vx(v, p.x), vy(v, p.y));
		}
		fprintf(s->f, "Z ");
	}
	fprintf(s->f, "\"/>\n");
	/* triangles (left over beside quads in a darker tone) */
	fprintf(s->f, "<path fill=\"%s\" fill-opacity=\"%s\" stroke=\"%s\" stroke-width=\"0.5\" stroke-opacity=\"0.85\" stroke-linejoin=\"round\" d=\"",
		nq ? "#cf222e" : fill, nq ? "0.30" : "0.35", line);
	for (i = 0; i < nt; i++)
	{
		unsigned int t[3];
		qaws_mesh_2d_get_triangle(m, i, t);
		for (k = 0; k < 3; k++)
		{
			qaws_vec2 p = qaws_mesh_2d_get_vertex(m, t[k]);
			fprintf(s->f, "%s%.2f,%.2f ", k ? "L" : "M", vx(v, p.x), vy(v, p.y));
		}
		fprintf(s->f, "Z ");
	}
	fprintf(s->f, "\"/>\n");
	for (i = 0; i < nb + nc; i++)
	{
		unsigned int e[2];
		qaws_vec2 a, b;
		if (i < nb) qaws_mesh_2d_get_boundary_edge(m, i, e);
		else qaws_mesh_2d_get_constraint_edge(m, i - nb, e);
		a = qaws_mesh_2d_get_vertex(m, e[0]);
		b = qaws_mesh_2d_get_vertex(m, e[1]);
		svg_line(s, vx(v, a.x), vy(v, a.y), vx(v, b.x), vy(v, b.y), i < nb ? "#24292f" : "#cf222e", i < nb ? 1.4 : 1.8, 0.95);
	}
}

static double tr_min_angle(qaws_mesh_2d const* m)
{
	unsigned int i, k, nt = qaws_mesh_2d_get_triangle_count(m);
	double best = 180;
	for (i = 0; i < nt; i++)
	{
		unsigned int t[3];
		qaws_vec2 p[3];
		double l[3];
		qaws_mesh_2d_get_triangle(m, i, t);
		for (k = 0; k < 3; k++) p[k] = qaws_mesh_2d_get_vertex(m, t[k]);
		for (k = 0; k < 3; k++) l[k] = hypot((double)p[(k + 2) % 3].x - p[(k + 1) % 3].x, (double)p[(k + 2) % 3].y - p[(k + 1) % 3].y);
		for (k = 0; k < 3; k++)
		{
			double b = l[(k + 1) % 3], c = l[(k + 2) % 3], cv = (b * b + c * c - l[k] * l[k]) / (2 * b * c);
			double a = acos(cv < -1 ? -1 : (cv > 1 ? 1 : cv)) * 180 / CO_PI;
			if (a < best) best = a;
		}
	}
	return best;
}

static qaws_scalar tr_size(void* user, qaws_vec2 p)
{
	(void)user;
	/* fine round the point (1.4, 1.2), coarse away from it */
	return (qaws_scalar)(0.06 + 0.22 * hypot((double)p.x - 1.4, (double)p.y - 1.2));
}

static void demo_triangulate(void)
{
	int const cols = 3, rows = 3, pw = 420, ph = 360, gap = 16, top = 84;
	int const W = cols * pw + (cols + 1) * gap, H = top + rows * (ph + gap) + 8;
	svg s;
	int k;
	double t0;
	char label[240];
	qaws_curve* star = co_star(10, 3, 1.6, 0, 0);
	qaws_curve *h1 = NULL, *h2 = NULL, *h1r = NULL, *h2r = NULL, *disc = NULL, *wave = NULL, *ring = NULL, *rbox = NULL;
	if (!svg_open(&s, "showcase/b2d8_triangulate.svg", W, H, "Meshing regions of curves",
		"qaws_triangulate: constrained Delaunay with exact predicates, boundary vertices on the true curves, Ruppert refinement, remeshing, constraint curves (red), quadrangulation"))
		return;
	qaws_curve_create_ellipse_2d((qaws_scalar)-0.9, (qaws_scalar)0.4, (qaws_scalar)0.55, (qaws_scalar)0.35, 0, &h1);
	qaws_curve_create_ellipse_2d((qaws_scalar)0.9, (qaws_scalar)-0.5, (qaws_scalar)0.4, (qaws_scalar)0.4, 0, &h2);
	qaws_curve_reverse(&h1r, h1);
	qaws_curve_reverse(&h2r, h2);
	qaws_curve_create_ellipse_2d(0, 0, (qaws_scalar)3.2, (qaws_scalar)2.4, 0, &disc);
	{
		/* a rounded box: a closed cubic B-spline through the corners' arcs */
		static double const w[8] = { -3.2, -0.7, -1.1, 1.5, 1.1, -1.9, 3.2, 0.3 };
		wave = co_bezier(w, 4);
	}
	qaws_curve_create_ellipse_2d((qaws_scalar)-1.1, (qaws_scalar)0.75, (qaws_scalar)0.8, (qaws_scalar)0.8, 0, &ring);
	rbox = co_star(8, 2.95, 2.7, 0, -0.15);
	t0 = (double)clock() / CLOCKS_PER_SEC;
	for (k = 0; k < 9; k++)
	{
		viewport v;
		int x0 = gap + (k % cols) * (pw + gap), y0 = top + (k / cols) * (ph + gap);
		qaws_curve const* cs[3];
		qaws_curve const* cn[2];
		qaws_path_2d p[3], cp[2];
		qaws_triangulate_desc d;
		qaws_mesh_2d* m = NULL;
		qaws_curve* poly = NULL;
		unsigned int i, n = 1;
		char const* title = "";
		char const *fill = "#c8e1ff", *line = "#0969da";
		double remesh = 0;
		int quads = -1;
		memset(&d, 0, sizeof(d));
		d.fill_rule = QAWS_FILL_NON_ZERO;
		co_view(&v, x0, y0, pw, ph, 0, -0.15, 7.6);
		cs[0] = star; cs[1] = h1r; cs[2] = h2r; n = 3;
		d.tolerance = (qaws_scalar)0.01;
		switch (k)
		{
		case 0:
			cs[0] = disc; n = 1;
			title = "an ellipse: constrained Delaunay, chord sag 0.01";
			break;
		case 1:
			title = "a star with two holes: the triangles inside only";
			fill = "#dafbe1"; line = "#1a7f37";
			break;
		case 2:
			d.min_angle = 30;
			title = "refined to 30 degrees (Ruppert)";
			fill = "#fbefff"; line = "#8250df";
			break;
		case 3:
			d.tolerance = (qaws_scalar)0.005; d.min_angle = 25; d.size_fn = tr_size;
			title = "a size field: fine round one point";
			fill = "#fff8c5"; line = "#9a6700";
			break;
		case 4:
		{
			static qaws_scalar const z[14] = { -3, -2, 3, 2, 3, -2, -3, 2, -1.5f, -2.6f, 0, 2.6f, 1.5f, -2.6f };
			qaws_curve_create_polyline_2d(z, 7, 1, &poly);
			cs[0] = poly; n = 1;
			d.fill_rule = QAWS_FILL_EVEN_ODD; d.min_angle = 25; d.tolerance = 0;
			title = "a self-crossing polygon, even-odd, 25 degrees";
			fill = "#ffd8b5"; line = "#bc4c00";
			break;
		}
		case 5:
			remesh = 0.22;
			title = "remeshed to edges of 0.22 (split, collapse, flip, smooth)";
			fill = "#dafbe1"; line = "#1a7f37";
			break;
		case 6:
			cs[0] = rbox; n = 1;
			cn[0] = wave; cn[1] = ring;
			cp[0].curves = &cn[0]; cp[0].curve_count = 1; cp[0].closed = 0;
			cp[1].curves = &cn[1]; cp[1].curve_count = 1; cp[1].closed = 1;
			d.constraints = cp; d.constraint_count = 2; d.min_angle = 28;
			title = "constraint curves (red) followed, refined to 28 degrees";
			break;
		case 7:
			remesh = 0.24; quads = 0;
			title = "quad-dominant: pairs by shape, augmenting paths";
			fill = "#fff8c5"; line = "#9a6700";
			break;
		default:
			cs[0] = rbox; n = 1;
			cn[0] = wave; cn[1] = ring;
			cp[0].curves = &cn[0]; cp[0].curve_count = 1; cp[0].closed = 0;
			cp[1].curves = &cn[1]; cp[1].curve_count = 1; cp[1].closed = 1;
			d.constraints = cp; d.constraint_count = 2;
			remesh = 0.3; quads = 1;
			title = "all quads, constrained: every face a convex quad";
			fill = "#fbefff"; line = "#8250df";
			break;
		}
		for (i = 0; i < n; i++) { p[i].curves = &cs[i]; p[i].curve_count = 1; p[i].closed = 1; }
		d.paths = p; d.path_count = n;
		svg_panel(&s, &v, title);
		if (qaws_triangulate_region_2d(&d, &m) == QAWS_STATUS_OK)
		{
			if (remesh > 0) qaws_mesh_2d_remesh(m, (qaws_scalar)remesh, NULL, NULL, 6);
			if (quads >= 0) qaws_mesh_2d_quadrangulate(m, quads, 4);
			tr_draw(&s, &v, m, fill, line);
		}
		if (m && quads >= 0)
			sprintf(label, "%u quads, %u triangles, %u vertices", qaws_mesh_2d_get_quad_count(m), qaws_mesh_2d_get_triangle_count(m), qaws_mesh_2d_get_vertex_count(m));
		else
			sprintf(label, "%u vertices, %u triangles, smallest angle %.1f", m ? qaws_mesh_2d_get_vertex_count(m) : 0, m ? qaws_mesh_2d_get_triangle_count(m) : 0, m ? tr_min_angle(m) : 0);
		svg_text(&s, x0 + 10, y0 + ph - 12, 11, "#57606a", "start", label);
		qaws_mesh_2d_destroy(m);
		qaws_curve_destroy(poly);
	}
	svg_close(&s);
	printf("    triangulate: 9 panels in %.3f s\n", (double)clock() / CLOCKS_PER_SEC - t0);
	qaws_curve_destroy(star); qaws_curve_destroy(h1); qaws_curve_destroy(h2);
	qaws_curve_destroy(h1r); qaws_curve_destroy(h2r); qaws_curve_destroy(disc);
	qaws_curve_destroy(wave); qaws_curve_destroy(ring); qaws_curve_destroy(rbox);
}
