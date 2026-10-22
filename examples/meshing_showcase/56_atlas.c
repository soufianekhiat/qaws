/*
 * Figure 56: an atlas (auto_seams, pack): figure 2's blob, closed, cut by
 * itself into charts (each triangle grouped by the axis direction nearest
 * its normal) and packed into the unit square:
 *   - the blob with the atlas's (u, v) lines;
 *   - the atlas, each chart coloured by its group.
 * Printed: the charts, the fill and the angle distortion.
 */

static void demo_atlas(void)
{
	static double const warm[3] = { 225, 190, 150 }, wire[3] = { 90, 80, 70 }, red[3] = { 200, 60, 50 }, blue[3] = { 50, 90, 200 };
	static double const group[6][3] = { { 230, 120, 110 }, { 240, 180, 100 }, { 130, 200, 120 }, { 110, 190, 200 }, { 140, 150, 230 }, { 210, 140, 210 } };
	er_canvas* c = er_create(1280, 640, 3, 252, 252, 250);
	qaws_mesh_3d* m = ms_blob(3);
	unsigned int nt, i, k;
	qaws_scalar* uv;
	qaws_param_desc d;
	qaws_param_report r;
	char sub[200];
	if (!c || !m)
		return;
	nt = qaws_mesh_3d_get_triangle_count(m);
	uv = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 6 * nt);
	memset(&d, 0, sizeof d);
	d.auto_seams = 1;
	d.pack = 1;
	if (qaws_mesh_3d_parameterize(uv, &r, m, NULL, 0, &d) == QAWS_STATUS_OK)
	{
		double* fx = (double*)malloc(sizeof(double) * 9 * nt);
		double* rgb = (double*)malloc(sizeof(double) * 9 * nt);
		unsigned int* ft = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nt);
		qaws_mesh_3d* flat;
		sprintf(sub, "%u charts, fill %.2f, angles %.3f, %u turned", r.charts, (double)r.fill, (double)r.angle_distortion, r.flipped);
		ms_panel(c, 0, 0, 640, 640, m, -0.6, 0.5, "a closed blob, cut by itself", sub);
		ms_draw(c, m, warm, NULL, 1, NULL, 0);
		ms_df_uvlines(c, m, uv, 0, 1.0 / 48, red);
		ms_df_uvlines(c, m, uv, 1, 1.0 / 48, blue);
		for (i = 0; i < nt; i++)
		{
			unsigned int t[3], g = 0;
			qaws_vec3 a, b, e;
			double u[3], w[3], nr[3], best = -1;
			qaws_mesh_3d_get_triangle(t, m, i);
			a = qaws_mesh_3d_get_vertex(m, t[0]);
			b = qaws_mesh_3d_get_vertex(m, t[1]);
			e = qaws_mesh_3d_get_vertex(m, t[2]);
			u[0] = b.x - a.x; u[1] = b.y - a.y; u[2] = b.z - a.z;
			w[0] = e.x - a.x; w[1] = e.y - a.y; w[2] = e.z - a.z;
			nr[0] = u[1] * w[2] - u[2] * w[1];
			nr[1] = u[2] * w[0] - u[0] * w[2];
			nr[2] = u[0] * w[1] - u[1] * w[0];
			for (k = 0; k < 3; k++)
				if (fabs(nr[k]) > best)
				{
					best = fabs(nr[k]);
					g = 2 * k + (nr[k] < 0);
				}
			for (k = 0; k < 3; k++)
			{
				fx[9 * i + 3 * k] = uv[6 * i + 2 * k];
				fx[9 * i + 3 * k + 1] = uv[6 * i + 2 * k + 1];
				fx[9 * i + 3 * k + 2] = 0;
				memcpy(rgb + 9 * i + 3 * k, group[g], sizeof group[g]);
				ft[3 * i + k] = 3 * i + k;
			}
		}
		flat = ms_make(fx, 3 * nt, ft, NULL, nt);
		ms_panel(c, 640, 0, 640, 640, flat, 0.0, 1.5707, "the atlas in the unit square", "the charts on shelves, one scale");
		ms_draw(c, flat, NULL, rgb, 0, wire, 0.3);
		printf("  %s\n", sub);
		qaws_mesh_3d_destroy(flat);
		free(fx);
		free(rgb);
		free(ft);
	}
	free(uv);
	qaws_mesh_3d_destroy(m);
	er_save(c, "showcase/mesh56_atlas.png");
	printf("mesh56_atlas.png\n");
	er_destroy(c);
}
