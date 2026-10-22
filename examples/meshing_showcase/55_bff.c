/*
 * Figure 55: boundary first flattening (QAWS_PARAM_BFF), figure 30's hill:
 *   - the hill with BFF's (u, v) lines, its border mapped to a circle;
 *   - the layouts in the plane: LSCM (free border), BFF with the border's
 *     lengths kept, BFF to a disc, BFF to a square whose corners are the
 *     midpoints of the hill's sides (bff_turning: pi / 2 there, 0 elsewhere;
 *     the 2 turned: the old corners' triangles, flat on a side).
 * Printed: the angle and area distortions.
 */

static void demo_bff(void)
{
	static double const warm[3] = { 225, 190, 150 }, wire[3] = { 90, 80, 70 }, red[3] = { 200, 60, 50 }, blue[3] = { 50, 90, 200 };
	static char const* const titles[4] = { "LSCM", "BFF, lengths kept", "BFF, a disc", "BFF, a square" };
	enum { N = 32 };
	er_canvas* c = er_create(2250, 640, 3, 252, 252, 250);
	double* x = (double*)malloc(sizeof(double) * 3 * (N + 1) * (N + 1));
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 6 * N * N), i, j, kk = 0, nt, k;
	qaws_mesh_3d* m;
	char sub[200];
	qaws_scalar turn[(N + 1) * (N + 1)];
	if (!c)
		return;
	for (j = 0; j <= N; j++)
		for (i = 0; i <= N; i++)
		{
			double u = 2.0 * i / N - 1, v = 2.0 * j / N - 1;
			x[3 * (j * (N + 1) + i)] = u;
			x[3 * (j * (N + 1) + i) + 1] = v;
			x[3 * (j * (N + 1) + i) + 2] = 0.7 * exp(-(u * u + v * v) / 0.18);
		}
	for (j = 0; j < N; j++)
		for (i = 0; i < N; i++)
		{
			unsigned int a = j * (N + 1) + i;
			t[kk++] = a; t[kk++] = a + 1; t[kk++] = a + N + 2;
			t[kk++] = a; t[kk++] = a + N + 2; t[kk++] = a + N + 1;
		}
	memset(turn, 0, sizeof turn);
	turn[N / 2] = turn[(N / 2) * (N + 1)] = turn[(N / 2) * (N + 1) + N] = turn[N * (N + 1) + N / 2] = (qaws_scalar)1.5707963267948966;
	m = ms_make(x, (N + 1) * (N + 1), t, NULL, 2 * N * N);
	nt = qaws_mesh_3d_get_triangle_count(m);
	for (k = 0; k < 4; k++)
	{
		qaws_scalar* uv = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 6 * nt);
		qaws_param_desc d;
		qaws_param_report r;
		memset(&d, 0, sizeof d);
		d.method = k ? QAWS_PARAM_BFF : QAWS_PARAM_LSCM;
		d.bff_disk = k == 2;
		d.bff_turning = k == 3 ? turn : NULL;
		if (qaws_mesh_3d_parameterize(uv, &r, m, NULL, 0, &d) == QAWS_STATUS_OK)
		{
			double* fx = (double*)malloc(sizeof(double) * 9 * nt);
			unsigned int* ft = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nt);
			qaws_mesh_3d* flat;
			if (k == 2)
			{
				double lo = HUGE_VAL, hi = -HUGE_VAL;
				for (i = 0; i < 6 * nt; i++)
				{
					lo = fmin(lo, (double)uv[i]);
					hi = fmax(hi, (double)uv[i]);
				}
				ms_panel(c, 0, 0, 450, 640, m, -0.6, 0.7, "a hill, BFF to a disc", "its (u, v) lines");
				ms_draw(c, m, warm, NULL, 1, NULL, 0);
				ms_df_uvlines(c, m, uv, 0, (hi - lo) / 16, red);
				ms_df_uvlines(c, m, uv, 1, (hi - lo) / 16, blue);
			}
			for (i = 0; i < 3 * nt; i++)
			{
				fx[3 * i] = uv[2 * i];
				fx[3 * i + 1] = uv[2 * i + 1];
				fx[3 * i + 2] = 0;
				ft[i] = i;
			}
			flat = ms_make(fx, 3 * nt, ft, NULL, nt);
			sprintf(sub, "angles %.3f, areas %.3f, %u turned", (double)r.angle_distortion, (double)r.area_distortion, r.flipped);
			ms_panel(c, 450 + 450 * (int)k, 0, 450, 640, flat, 0.0, 1.5707, titles[k], sub);
			ms_draw(c, flat, warm, NULL, 0, wire, 0.3);
			printf("  %s: %s\n", titles[k], sub);
			qaws_mesh_3d_destroy(flat);
			free(fx);
			free(ft);
		}
		free(uv);
	}
	qaws_mesh_3d_destroy(m);
	free(x);
	free(t);
	er_save(c, "showcase/mesh55_bff.png");
	printf("mesh55_bff.png\n");
	er_destroy(c);
}
