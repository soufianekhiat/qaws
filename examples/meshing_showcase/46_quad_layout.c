/*
 * Figure 46: quad layouts (QAWS_QUADS_LAYOUT): GLOBAL's quads parted by
 * their motorcycle graph (separatrices from the irregular vertices, straight
 * on until the border, a crease or an earlier one), each patch a grid of
 * quads in its own colour.
 *   - figure 2's blob;
 *   - the can (its rims creases);
 *   - the torus: no irregular vertex, one patch round it;
 *   - the wavy patch of figure 5, its border's corners.
 */

/* the quads coloured by their face attribute "patch", outlined */
static void ms_ly_draw(er_canvas* c, qaws_mesh_3d const* m, double const* wire)
{
	unsigned int nf = qaws_mesh_3d_get_face_count(m), nc = qaws_mesh_3d_get_corner_count(m), f, q = 0, k, id;
	unsigned int *fv = (unsigned int*)malloc(sizeof(unsigned int) * nc), *fs = (unsigned int*)malloc(sizeof(unsigned int) * nf);
	qaws_scalar* pa = (qaws_scalar*)calloc(nf + 1, sizeof(qaws_scalar));
	double* x = ms_positions(m);
	qaws_mesh_3d_get_faces(fv, fs, m);
	if (qaws_mesh_3d_find_attribute(&id, m, "patch") == QAWS_STATUS_OK)
		qaws_mesh_3d_get_attribute(pa, m, id);
	for (f = 0; f < nf; q += fs[f++])
	{
		double col[3], t = fmod(0.6180339887 * (double)pa[f] + 0.13, 1.0);
		ms_colormap(t, col);
		for (k = 1; k + 1 < fs[f]; k++)
			er_triangle(c, x + 3 * fv[q], x + 3 * fv[q + k], x + 3 * fv[q + k + 1], NULL, NULL, NULL, col, col, col);
		for (k = 0; k < fs[f]; k++)
			er_line(c, x + 3 * fv[q + k], x + 3 * fv[q + (k + 1) % fs[f]], wire, 0.6);
	}
	free(fv);
	free(fs);
	free(pa);
	free(x);
}

static void ms_ly_panel(er_canvas* c, int px, qaws_mesh_3d const* m, double L, double yaw, double pitch, char const* title, double const* wire)
{
	qaws_quadrangulate_desc d;
	qaws_mesh_3d* q = NULL;
	qaws_mesh_report r;
	char sub[200];
	memset(&d, 0, sizeof d);
	d.method = QAWS_QUADS_LAYOUT;
	d.target_length = (qaws_scalar)L;
	if (qaws_mesh_quadrangulate(&q, &r, m, &d) != QAWS_STATUS_OK)
	{
		printf("  %s: failed\n", title);
		return;
	}
	sprintf(sub, "%u quads in %u patch%s%s", r.quality.quads, r.splits, r.splits == 1 ? "" : "es", r.iterations ? "" : " (INSTANT's)");
	ms_panel(c, px, 0, 450, 640, q, yaw, pitch, title, sub);
	ms_ly_draw(c, q, wire);
	printf("  %s: %s, %u irregular\n", title, sub, r.quality.irregular);
	qaws_mesh_3d_destroy(q);
}

static void demo_quad_layout(void)
{
	static double const wire[3] = { 40, 45, 60 };
	er_canvas* c = er_create(1800, 640, 3, 252, 252, 250);
	qaws_mesh_3d *blob = ms_blob(1), *can = ms_can(8, 1, 1.6), *torus = ms_torus(24, 12, 1, 0.4), *wavy = NULL;
	if (!c || !blob || !can || !torus)
		return;
	qaws_mesh_3d_detect_features(can, (qaws_scalar)(50 * MS_PI / 180), 0, 0);
	ms_ly_panel(c, 0, blob, 0.16, 0.5, 0.45, "a blob", wire);
	ms_ly_panel(c, 450, can, 0.25, 0.6, 0.6, "a can", wire);
	ms_ly_panel(c, 900, torus, 0.2, 0.4, 0.8, "a torus", wire);
	{
		/* figure 5's wavy patch: a bicubic B-spline sampled on a grid */
		qaws_vec3 cp[49], s[24 * 24];
		qaws_surface_bspline_desc bd;
		qaws_surface* srf = NULL;
		unsigned int i, j, k = 0, got = 0, n = 24, f[6 * 23 * 23];
		double x[3 * 24 * 24];
		for (i = 0; i < 7; i++)
			for (j = 0; j < 7; j++)
			{
				cp[i * 7 + j].x = (qaws_scalar)(i / 3.0 - 1);
				cp[i * 7 + j].y = (qaws_scalar)(j / 3.0 - 1);
				cp[i * 7 + j].z = (qaws_scalar)(0.55 * sin(1.7 * i) * cos(1.3 * j) + (i == 3 && j == 3 ? 0.6 : 0));
			}
		memset(&bd, 0, sizeof bd);
		bd.u_degree = 3;
		bd.v_degree = 3;
		bd.control_points = cp;
		bd.u_point_count = 7;
		bd.v_point_count = 7;
		if (qaws_surface_create_bspline(&srf, &bd) == QAWS_STATUS_OK && qaws_surface_sample(s, &got, srf, n, n, n * n) == QAWS_STATUS_OK && got == n * n)
		{
			for (i = 0; i < n * n; i++)
			{
				x[3 * i] = (double)s[i].x;
				x[3 * i + 1] = (double)s[i].y;
				x[3 * i + 2] = (double)s[i].z;
			}
			for (i = 0; i + 1 < n; i++)
				for (j = 0; j + 1 < n; j++)
				{
					unsigned int a = i * n + j, b = (i + 1) * n + j, cc = (i + 1) * n + j + 1, e = i * n + j + 1;
					f[k++] = a; f[k++] = b; f[k++] = cc;
					f[k++] = a; f[k++] = cc; f[k++] = e;
				}
			wavy = ms_make(x, n * n, f, NULL, k / 3);
		}
		qaws_surface_destroy(srf);
	}
	if (wavy)
		ms_ly_panel(c, 1350, wavy, 0.18, 0.5, 0.7, "a wavy patch", wire);
	qaws_mesh_3d_destroy(blob);
	qaws_mesh_3d_destroy(can);
	qaws_mesh_3d_destroy(torus);
	qaws_mesh_3d_destroy(wavy);
	er_save(c, "showcase/mesh46_quad_layout.png");
	printf("mesh46_quad_layout.png\n");
	er_destroy(c);
}
