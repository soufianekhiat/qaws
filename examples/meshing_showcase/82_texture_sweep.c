/*
 * Figure 82: a sweep's attributes per texel (qaws_mesh_3d_sweep_texture).
 * A unit cube moves 1 along x slowing to a stop; its swept solid's shape,
 * a 2 x 1 x 1 box of 12 triangles, coloured by the time each place was
 * first touched (blue 0 to yellow 1):
 *   - from its eight vertices' exact times, interpolated: the top a linear
 *     ramp, the slowing lost;
 *   - from a 128 x 128 texture of the time evaluated per texel (the box's
 *     charts packed, a gutter of 4): the ramp bunched toward the stop, as
 *     the motion is.
 * Printed: the texels written and the top's colours against the exact time.
 */

/* the first touch at x on the box's top: 1 - (2 - x)^(1/3) past the cube at rest */
static double ms_ts_time(double x)
{
	return x <= 1 ? 0 : 1 - cbrt(2 - x);
}

static void demo_texture_sweep(void)
{
	enum { R = 128, K = 16 };
	static double const wire[3] = { 60, 66, 80 };
	er_canvas* c = er_create(1200, 640, 3, 252, 252, 250);
	qaws_mesh_3d *b = ms_sw_box(0, 1, 0, 1, 0, 1), *o = ms_sw_box(0, 2, 0, 1, 0, 1);
	qaws_motion* m = NULL;
	qaws_motion_pose s, e;
	qaws_sweep_desc d;
	qaws_texture t;
	qaws_scalar* tex = (qaws_scalar*)calloc(R * R, sizeof(qaws_scalar));
	unsigned char* cov = (unsigned char*)malloc(R * R);
	unsigned int nt = o ? qaws_mesh_3d_get_triangle_count(o) : 0, id, i, written = 0, panel;
	qaws_scalar* uv = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 6 * (nt + 1));
	unsigned int* tri = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (nt + 1));
	double* x = o ? ms_positions(o) : NULL;
	qaws_param_desc pd;
	qaws_attribute_desc ad;
	qaws_texture_desc td;
	char sub[200];
	if (!c || !b || !o || !tex || !cov || !uv || !tri || !x)
		return;
	memset(&s, 0, sizeof s);
	memset(&e, 0, sizeof e);
	s.rotation[0] = e.rotation[0] = 1;
	e.position.x = 1;
	s.velocity.x = 3;
	qaws_motion_create_poses(&m, &s, &e, 0);
	memset(&pd, 0, sizeof pd);
	pd.auto_seams = 1;
	pd.pack = 1;
	qaws_mesh_3d_parameterize(uv, NULL, o, NULL, 0, &pd);
	memset(&ad, 0, sizeof ad);
	ad.name = "uv"; ad.element = QAWS_MESH_CORNER; ad.components = 2; ad.rule = QAWS_ATTRIBUTE_LINEAR;
	qaws_mesh_3d_add_attribute(&id, o, uv, &ad);
	memset(&t, 0, sizeof t);
	t.texels = tex; t.width = R; t.height = R; t.components = 1;
	memset(&d, 0, sizeof d);
	d.tolerance = (qaws_scalar)0.002;
	memset(&td, 0, sizeof td);
	td.gutter = 4;   /* (bilinear reads by a chart's border take the gutter, not texels of no chart) */
	if (m && qaws_mesh_3d_sweep_texture(cov, &t, o, b, m, "time", &d, &td) == QAWS_STATUS_OK)
		for (i = 0; i < R * R; i++)
			written += cov[i] == 1;
	qaws_mesh_3d_get_triangles(tri, o);
	for (panel = 0; panel < 2; panel++)
	{
		double worst = 0;
		ms_panel(c, 600 * (int)panel, 0, 600, 640, o, -0.6, 0.6, panel ? "time from 128 x 128 texels" : "time from 8 vertices", NULL);
		/* (each triangle cut in K x K for the drawing, its pieces' corners coloured by the time read there) */
		for (i = 0; i < nt; i++)
		{
			unsigned int a2, b2, h, z, k, q;
			double P[3][3], U[3][2], T[3];
			for (q = 0; q < 3; q++)
			{
				for (k = 0; k < 3; k++)
					P[q][k] = x[3 * tri[3 * i + q] + k];
				U[q][0] = (double)uv[6 * i + 2 * q];
				U[q][1] = (double)uv[6 * i + 2 * q + 1];
				T[q] = ms_ts_time(P[q][0]);
			}
			for (b2 = 0; b2 < K; b2++)
				for (a2 = 0; a2 + b2 < K; a2++)
					for (h = 0; h < (a2 + b2 + 1 < K ? 2u : 1u); h++)
					{
						unsigned int g[3][2];
						double Q[3][3], C[3][3];
						if (!h)
						{
							g[0][0] = a2; g[0][1] = b2; g[1][0] = a2 + 1; g[1][1] = b2; g[2][0] = a2; g[2][1] = b2 + 1;
						}
						else
						{
							g[0][0] = a2 + 1; g[0][1] = b2; g[1][0] = a2 + 1; g[1][1] = b2 + 1; g[2][0] = a2; g[2][1] = b2 + 1;
						}
						for (z = 0; z < 3; z++)
						{
							double w1 = (double)g[z][0] / K, w2 = (double)g[z][1] / K, w0 = 1 - w1 - w2, tt;
							for (k = 0; k < 3; k++)
								Q[z][k] = w0 * P[0][k] + w1 * P[1][k] + w2 * P[2][k];
							if (panel)
							{
								qaws_scalar suv[2], y[1];
								suv[0] = (qaws_scalar)(w0 * U[0][0] + w1 * U[1][0] + w2 * U[2][0]);
								suv[1] = (qaws_scalar)(w0 * U[0][1] + w1 * U[1][1] + w2 * U[2][1]);
								qaws_texture_sample(y, &t, suv, 1);
								tt = (double)y[0];
							}
							else
								tt = w0 * T[0] + w1 * T[1] + w2 * T[2];
							/* (the top, away from its rims: against the exact time) */
							if (fabs(Q[z][2] - 1) < 1e-12 && Q[z][0] > 1.05 && Q[z][0] < 1.95 && Q[z][1] > 0.05 && Q[z][1] < 0.95)
								worst = fmax(worst, fabs(tt - ms_ts_time(Q[z][0])));
							ms_sg_ramp(tt, C[z]);
						}
						er_triangle(c, Q[0], Q[1], Q[2], NULL, NULL, NULL, C[0], C[1], C[2]);
					}
		}
		ms_draw(c, o, NULL, NULL, 0, wire, 0.6);
		if (panel)
			sprintf(sub, "%u texels; the top off the time by %.4f", written, worst);
		else
			sprintf(sub, "12 triangles; the top off the time by %.4f", worst);
		er_text(c, 600 * (int)panel + 14, 46, sub, 2, 95, 95, 95);
		printf("  %s: %s\n", panel ? "texels" : "vertices", sub);
	}
	er_save(c, "showcase/mesh82_texture_sweep.png");
	printf("mesh82_texture_sweep.png\n");
	qaws_motion_destroy(m);
	qaws_mesh_3d_destroy(b); qaws_mesh_3d_destroy(o);
	free(tex); free(cov); free(uv); free(tri); free(x);
	er_destroy(c);
}
