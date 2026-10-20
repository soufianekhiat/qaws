/*
 * Figure 5: meshes of exact geometry (qaws_mesh_generate_surface/_shell).
 *   - a wavy bicubic B-spline patch: uniform triangles, all quads, and
 *     triangles adaptive by a chord tolerance;
 *   - a rational NURBS torus (two NURBS circles' tensor product): quads,
 *     every vertex on the torus;
 *   - a box of six bilinear faces as a shell: quads, the box's edges kept,
 *     the faces coloured.
 * Every vertex lies on its face's surface (snapped: its footpoint the face
 * and its (u, v)).
 */

static qaws_surface* ms_wavy(void)
{
	qaws_vec3 cp[49];
	qaws_surface_bspline_desc d;
	qaws_surface* s = NULL;
	unsigned int i, j;
	for (i = 0; i < 7; i++)
		for (j = 0; j < 7; j++)
		{
			cp[i * 7 + j].x = (qaws_scalar)(i / 3.0 - 1);
			cp[i * 7 + j].y = (qaws_scalar)(j / 3.0 - 1);
			cp[i * 7 + j].z = (qaws_scalar)(0.55 * sin(1.7 * i) * cos(1.3 * j) + (i == 3 && j == 3 ? 0.6 : 0));
		}
	memset(&d, 0, sizeof d);
	d.u_degree = 3;
	d.v_degree = 3;
	d.control_points = cp;
	d.u_point_count = 7;
	d.v_point_count = 7;
	qaws_surface_create_bspline(&s, &d);
	return s;
}

/* the torus of radii R, r as a rational NURBS surface: the unit circle's 9 points times the profile circle's */
static qaws_surface* ms_nurbs_torus(double R, double r)
{
	static double const cx[9] = { 1, 1, 0, -1, -1, -1, 0, 1, 1 }, cy[9] = { 0, 1, 1, 1, 0, -1, -1, -1, 0 };
	static qaws_scalar const kn[12] = { 0, 0, 0, 0.25f, 0.25f, 0.5f, 0.5f, 0.75f, 0.75f, 1, 1, 1 };
	qaws_vec3 cp[81];
	qaws_scalar w[81];
	qaws_surface_nurbs_desc d;
	qaws_surface* s = NULL;
	unsigned int i, j;
	double h = sqrt(0.5);
	for (i = 0; i < 9; i++)
		for (j = 0; j < 9; j++)
		{
			double rad = R + r * cx[j];
			cp[i * 9 + j].x = (qaws_scalar)(rad * cx[i]);
			cp[i * 9 + j].y = (qaws_scalar)(rad * cy[i]);
			cp[i * 9 + j].z = (qaws_scalar)(r * cy[j]);
			w[i * 9 + j] = (qaws_scalar)((i & 1 ? h : 1) * (j & 1 ? h : 1));
		}
	memset(&d, 0, sizeof d);
	d.u_degree = 2;
	d.v_degree = 2;
	d.control_points = cp;
	d.u_point_count = 9;
	d.v_point_count = 9;
	d.weights = w;
	d.u_knots = kn;
	d.u_knot_count = 12;
	d.v_knots = kn;
	d.v_knot_count = 12;
	qaws_surface_create_nurbs(&s, &d);
	return s;
}

/* a box's six faces (bilinear), outward */
static void ms_box_faces(qaws_surface** f, double a, double b, double c)
{
	static int const q[6][4][3] = {
		{ { -1, -1, -1 }, { -1, 1, -1 }, { 1, -1, -1 }, { 1, 1, -1 } }, { { -1, -1, 1 }, { 1, -1, 1 }, { -1, 1, 1 }, { 1, 1, 1 } },
		{ { -1, -1, -1 }, { 1, -1, -1 }, { -1, -1, 1 }, { 1, -1, 1 } }, { { -1, 1, -1 }, { -1, 1, 1 }, { 1, 1, -1 }, { 1, 1, 1 } },
		{ { -1, -1, -1 }, { -1, -1, 1 }, { -1, 1, -1 }, { -1, 1, 1 } }, { { 1, -1, -1 }, { 1, 1, -1 }, { 1, -1, 1 }, { 1, 1, 1 } } };
	unsigned int i;
	for (i = 0; i < 6; i++)
	{
		qaws_surface_bilinear_desc d;
		qaws_vec3* p[4] = { &d.p00, &d.p10, &d.p01, &d.p11 };
		unsigned int k;
		for (k = 0; k < 4; k++)
		{
			p[k]->x = (qaws_scalar)(a * q[i][k][0]);
			p[k]->y = (qaws_scalar)(b * q[i][k][1]);
			p[k]->z = (qaws_scalar)(c * q[i][k][2]);
		}
		f[i] = NULL;
		qaws_surface_create_bilinear(&f[i], &d);
	}
}

static qaws_mesh_3d* ms_generate(qaws_surface const* s, qaws_shell_atlas const* a, qaws_elements el, double L, double chord, qaws_mesh_report* rep)
{
	qaws_generate_desc d;
	qaws_mesh_3d* out = NULL;
	qaws_status st;
	memset(&d, 0, sizeof d);
	d.elements = el;
	d.target_length = (qaws_scalar)L;
	if (chord > 0)
	{
		d.chord_tolerance = (qaws_scalar)chord;
		d.min_length = (qaws_scalar)0.03;
		d.max_length = (qaws_scalar)0.3;
		d.gradation = (qaws_scalar)0.4;
	}
	st = s ? qaws_mesh_generate_surface(s, &d, &out, rep) : qaws_mesh_generate_shell(a, &d, &out, rep);
	if (st != QAWS_STATUS_OK)
		printf("  generate failed: %s\n", qaws_status_to_string(st));
	return out;
}

static void demo_generate(void)
{
	static double const wire[3] = { 40, 45, 60 }, tint[3] = { 178, 192, 212 }, gold[3] = { 232, 196, 120 }, green[3] = { 170, 205, 160 };
	static double const pal[6][3] = { { 230, 97, 1 }, { 253, 184, 99 }, { 178, 171, 210 }, { 94, 60, 153 }, { 120, 198, 121 }, { 49, 130, 189 } };
	er_canvas* c = er_create(2750, 640, 3, 252, 252, 250);
	qaws_surface *wavy = ms_wavy(), *torus = ms_nurbs_torus(2, 0.7), *box[6];
	qaws_shell_atlas* atlas = NULL;
	qaws_mesh_3d *g1 = NULL, *g2 = NULL, *g3 = NULL, *g4 = NULL, *g5 = NULL;
	qaws_mesh_report rep;
	char sub[200];
	if (!c || !wavy || !torus)
		return;
	ms_box_faces(box, 1, 0.7, 0.5);
	qaws_shell_atlas_create(&atlas, (qaws_surface const* const*)box, 6, 0);
	g1 = ms_generate(wavy, NULL, QAWS_ELEMENTS_TRIANGLES, 0.1, 0, &rep);
	if (g1)
	{
		sprintf(sub, "%u triangles, min angle %.0f", qaws_mesh_3d_get_face_count(g1), (double)rep.quality.min_angle);
		ms_panel(c, 0, 0, 550, 640, g1, 0.6, 0.65, "a B-spline patch", sub);
		ms_draw(c, g1, tint, NULL, 1, wire, 0.9);
		printf("  wavy triangles: %u, angles %.1f..%.1f, snapped max %.2e\n", qaws_mesh_3d_get_face_count(g1), (double)rep.quality.min_angle,
			(double)rep.quality.max_angle, (double)rep.max_deviation);
	}
	g2 = ms_generate(wavy, NULL, QAWS_ELEMENTS_QUADS, 0.1, 0, &rep);
	if (g2)
	{
		sprintf(sub, "%u quads, %u singular", rep.quality.quads, rep.quality.irregular);
		ms_panel(c, 550, 0, 550, 640, g2, 0.6, 0.65, "in quads", sub);
		ms_draw(c, g2, tint, NULL, 1, wire, 0.9);
		printf("  wavy quads: %u quads, %u others, snapped max %.2e\n", rep.quality.quads, rep.quality.triangles + rep.quality.polygons, (double)rep.max_deviation);
	}
	g3 = ms_generate(wavy, NULL, QAWS_ELEMENTS_TRIANGLES, 0, 0.002, &rep);
	if (g3)
	{
		sprintf(sub, "chord 0.002: %u triangles", qaws_mesh_3d_get_face_count(g3));
		ms_panel(c, 1100, 0, 550, 640, g3, 0.6, 0.65, "adaptive", sub);
		ms_draw(c, g3, gold, NULL, 1, wire, 0.9);
		printf("  wavy adaptive: %u triangles, edges %.3f..%.3f, snapped max %.2e\n", qaws_mesh_3d_get_face_count(g3), (double)rep.quality.min_edge,
			(double)rep.quality.max_edge, (double)rep.max_deviation);
	}
	g4 = ms_generate(torus, NULL, QAWS_ELEMENTS_QUADS, 0.25, 0, &rep);
	if (g4)
	{
		double* p = ms_positions(g4), off = 0;
		unsigned int i, nv = qaws_mesh_3d_get_vertex_count(g4);
		for (i = 0; i < nv; i++)
		{
			double rho = sqrt(p[3 * i] * p[3 * i] + p[3 * i + 1] * p[3 * i + 1]);
			off = fmax(off, fabs(sqrt((rho - 2) * (rho - 2) + p[3 * i + 2] * p[3 * i + 2]) - 0.7));
		}
		sprintf(sub, "%u quads, on the torus to %.0e", rep.quality.quads, off);
		ms_panel(c, 1650, 0, 550, 640, g4, 0.5, 0.75, "a NURBS torus", sub);
		ms_draw(c, g4, green, NULL, 1, wire, 0.9);
		printf("  nurbs torus: %u quads, %u others, off the torus %.2e\n", rep.quality.quads, rep.quality.triangles + rep.quality.polygons, off);
		free(p);
	}
	g5 = atlas ? ms_generate(NULL, atlas, QAWS_ELEMENTS_QUADS, 0.15, 0, &rep) : NULL;
	if (g5)
	{
		unsigned int nf = qaws_mesh_3d_get_face_count(g5), nt = qaws_mesh_3d_get_triangle_count(g5), *t, f, a, q = 0, ti = 0;
		qaws_scalar* fid = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nf);
		unsigned int* fs = (unsigned int*)malloc(sizeof(unsigned int) * nf);
		unsigned int* fv = (unsigned int*)malloc(sizeof(unsigned int) * qaws_mesh_3d_get_corner_count(g5));
		double* p = ms_positions(g5);
		t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nt);
		qaws_mesh_3d_find_attribute(g5, "face", &a);
		qaws_mesh_3d_get_attribute(g5, a, fid);
		qaws_mesh_3d_get_faces(g5, fv, fs);
		sprintf(sub, "a shell of 6 faces: %u quads", rep.quality.quads);
		ms_panel(c, 2200, 0, 550, 640, g5, 0.62, 0.55, "a box", sub);
		for (f = 0; f < nf; f++)
		{
			unsigned int k, id = (unsigned int)fid[f] % 6;
			for (k = 1; k + 1 < fs[f]; k++)
				er_triangle(c, p + 3 * fv[q], p + 3 * fv[q + k], p + 3 * fv[q + k + 1], NULL, NULL, NULL, pal[id], pal[id], pal[id]);
			q += fs[f];
		}
		(void)ti;
		ms_draw(c, g5, NULL, NULL, 0, wire, 1.0);
		printf("  box shell: %u quads, %u others, snapped max %.2e\n", rep.quality.quads, rep.quality.triangles + rep.quality.polygons, (double)rep.max_deviation);
		free(fid); free(fs); free(fv); free(p); free(t);
	}
	er_save(c, "showcase/mesh5_generate.png");
	printf("mesh5_generate.png\n");
	qaws_mesh_3d_destroy(g1);
	qaws_mesh_3d_destroy(g2);
	qaws_mesh_3d_destroy(g3);
	qaws_mesh_3d_destroy(g4);
	qaws_mesh_3d_destroy(g5);
	qaws_shell_atlas_destroy(atlas);
	qaws_surface_destroy(wavy);
	qaws_surface_destroy(torus);
	{
		unsigned int i;
		for (i = 0; i < 6; i++)
			qaws_surface_destroy(box[i]);
	}
	er_destroy(c);
}
