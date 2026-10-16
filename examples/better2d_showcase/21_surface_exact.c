/* Figure 21: Booleans on surfaces with exact topology
   (qaws_surface_clip_exact). On a sphere, shapes drawn in (u, v) whose
   edges coincide: a ring of sixteen squares sharing their sides (across the
   seam too) united into one band, a checkerboard's xor, squares touching
   at corners, and two geodesic circles' intersection. Shared and touching
   edges are decided exactly, not within a tolerance. Uses figures 16 to
   19's helpers. */

static void sxf_square(qaws_surface const* s, double u0, double v0, double u1, double v1, qaws_curve** k, qaws_surface_curve** c)
{
	double x[5] = { u0, u1, u1, u0, u0 }, y[5] = { v0, v0, v1, v1, v0 };
	unsigned int i;
	for (i = 0; i < 4; i++)
	{
		qaws_bezier_desc bz;
		qaws_scalar cp[4];
		cp[0] = (qaws_scalar)x[i]; cp[1] = (qaws_scalar)y[i]; cp[2] = (qaws_scalar)x[i + 1]; cp[3] = (qaws_scalar)y[i + 1];
		memset(&bz, 0, sizeof(bz));
		bz.dimension = QAWS_DIMENSION_2D; bz.degree = 1; bz.control_points = cp; bz.control_point_count = 2;
		k[i] = NULL; c[i] = NULL;
		qaws_curve_create_bezier(&k[i], &bz);
		qaws_surface_curve_create_uv(s, k[i], &c[i]);
	}
}

static void demo_surface_exact(void)
{
	clock_t t0 = clock();
	svg s;
	double const pw = 300, ph = 260, pi = 3.14159265358979323846;
	unsigned int i, j, k;
	char label[220];
	qaws_curve* prof = sgf_arc(0, 0, 1, -pi / 2, pi / 2);
	qaws_surface* sph = sgf_revolve(prof);
	/* sixteen squares round the sphere, two rows of eight */
	qaws_curve* kk[16][4];
	qaws_surface_curve* cc[16][4];
	qaws_surface_curve const* pc[16][4];
	qaws_surface_path P[16];
	if (!svg_open(&s, "showcase/b2d21_surface_exact.svg", (int)(48 + 4 * (pw + 12)), (int)(110 + 2 * (ph + 44)),
		"Booleans on surfaces with exact topology",
		"qaws_surface_clip_exact: shared and touching edges decided exactly (qaws_clip_exact on the (u, v) images, the seam's cuts included)"))
		return;
	for (i = 0; i < 16; i++)
	{
		double u0 = (i % 8) / 8.0, v0 = 0.375 + 0.125 * (i / 8);
		sxf_square(sph, u0, v0, u0 + 0.125, v0 + 0.125, kk[i], cc[i]);
		for (j = 0; j < 4; j++) pc[i][j] = cc[i][j];
		P[i].curves = pc[i]; P[i].curve_count = 4; P[i].closed = 1;
	}
	for (k = 0; k < 4; k++)
	{
		viewport v3, vuv;
		sgf_view view;
		qaws_surface_region* r = NULL;
		qaws_scalar area = 0;
		qaws_status st = QAWS_STATUS_OK;
		char const* name = "";
		double expect = -1;
		v3.x0 = 24 + k * (pw + 12); v3.y0 = 104; v3.w = pw; v3.h = ph;
		vuv = v3; vuv.y0 = 104 + ph + 44;
		vuv.xmin = -0.1; vuv.xmax = 1.1; vuv.ymin = 0.5 - 0.6 * ph / pw; vuv.ymax = 0.5 + 0.6 * ph / pw;
		sgf_view_init(&view, -pi / 2 + 0.2, 0.3, 108, v3.x0 + pw / 2, v3.y0 + ph / 2 + 8);
		if (k == 0)
		{
			/* the ring united: the band between latitudes -pi/8 and pi/8 */
			st = qaws_surface_clip_exact(sph, QAWS_CLIP_NONE, QAWS_FILL_NON_ZERO, P, 16, NULL, 0, NULL, &r);
			name = "16 squares sharing sides, united";
			expect = 2 * pi * (sin(pi / 8) - sin(-pi / 8));
		}
		else if (k == 1)
		{
			/* the checkerboard: every other square, xor with the ring */
			qaws_surface_path even[8];
			unsigned int n = 0;
			for (i = 0; i < 16; i++) if (((i % 8) + i / 8) % 2 == 0) even[n++] = P[i];
			st = qaws_surface_clip_exact(sph, QAWS_CLIP_XOR, QAWS_FILL_NON_ZERO, P, 16, even, n, NULL, &r);
			name = "the ring xor a checkerboard of it";
			expect = 0.5 * 2 * pi * (sin(pi / 8) - sin(-pi / 8));
		}
		else if (k == 2)
		{
			/* squares meeting only at corners: a checkerboard's one colour */
			qaws_surface_path even[8];
			unsigned int n = 0;
			for (i = 0; i < 16; i++) if (((i % 8) + i / 8) % 2 == 1) even[n++] = P[i];
			st = qaws_surface_clip_exact(sph, QAWS_CLIP_NONE, QAWS_FILL_NON_ZERO, even, n, NULL, 0, NULL, &r);
			name = "squares touching at corners only";
			expect = 0.5 * 2 * pi * (sin(pi / 8) - sin(-pi / 8));
		}
		else
		{
			qaws_surface_curve *c1 = NULL, *c2 = NULL;
			qaws_surface_curve const *p1[1], *p2[1];
			qaws_surface_path A, B;
			qaws_surface_curve_create_geodesic_arc(sph, sgf_uv(0.24, 0.52), (qaws_scalar)0.5, 0, (qaws_scalar)(2 * pi), 0, &c1);
			qaws_surface_curve_create_geodesic_arc(sph, sgf_uv(0.32, 0.47), (qaws_scalar)0.45, 0, (qaws_scalar)(2 * pi), 0, &c2);
			p1[0] = c1; p2[0] = c2;
			ssf_path(&A, p1, 1); ssf_path(&B, p2, 1);
			st = qaws_surface_clip_exact(sph, QAWS_CLIP_INTERSECTION, QAWS_FILL_NON_ZERO, &A, 1, &B, 1, NULL, &r);
			sgf_view_init(&view, -pi / 2 - 2 * pi * 0.28 + 0.2, 0.3, 108, v3.x0 + pw / 2, v3.y0 + ph / 2 + 8);
			name = "two geodesic circles' lens";
			qaws_surface_curve_destroy(c1); qaws_surface_curve_destroy(c2);
		}
		if (st == QAWS_STATUS_OK) qaws_surface_region_get_area(r, &area);
		svg_panel(&s, &v3, name);
		if (expect >= 0)
			sprintf(label, "%u loops, area %.10f, exact %.10f", r ? qaws_surface_region_get_path_count(r) : 0, (double)area, expect);
		else
			sprintf(label, "%u loop, area %.6f (status %d)", r ? qaws_surface_region_get_path_count(r) : 0, (double)area, (int)st);
		svg_panel(&s, &vuv, label);
		sgf_mesh(&s, sph, &view, 12, 1);
		ssf_box(&s, &vuv, sph);
		if (r)
		{
			ssf_fill3d(&s, r, sph, &view, "#54aeff");
			sbf_fill(&s, r, sph, NULL, &vuv, 1, "#54aeff", "none");
			sof_outline(&s, r, sph, &view, &vuv, 1, 0, "#0969da", 1.6);
		}
		qaws_surface_region_destroy(r);
	}
	for (i = 0; i < 16; i++)
		for (j = 0; j < 4; j++) { qaws_surface_curve_destroy(cc[i][j]); qaws_curve_destroy(kk[i][j]); }
	qaws_surface_destroy(sph);
	qaws_curve_destroy(prof);
	svg_close(&s);
	printf("figure 21: exact Booleans on surfaces, %.1f s\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
}
