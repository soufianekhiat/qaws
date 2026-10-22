/*
 * Figure 56: an atlas (auto_seams, pack, max_stretch): figure 2's blob,
 * closed, cut by itself into charts (each triangle grouped by the axis
 * direction nearest its normal, the groups' borders smoothed) and packed
 * into the unit square:
 *   - the blob with the split atlas's (u, v) lines;
 *   - the atlas;
 *   - its charts split until their area scales stay within 1.3.
 * The atlases' triangles are coloured by their area scale (plane over mesh)
 * over the atlas's mean: blue shrunk, red grown (a factor 2 at either end).
 * Printed: the charts, the fill and the distortions.
 */

static void demo_atlas(void)
{
	static double const warm[3] = { 225, 190, 150 }, wire[3] = { 90, 80, 70 }, red[3] = { 200, 60, 50 }, blue[3] = { 50, 90, 200 };
	static char const* const titles[2] = { "the atlas in the unit square", "charts split, stretch within 1.3" };
	er_canvas* c = er_create(1920, 640, 3, 252, 252, 250);
	qaws_mesh_3d* m = ms_blob(3);
	unsigned int nt, i, k, z;
	qaws_scalar* uv;
	qaws_param_desc d;
	qaws_param_report r;
	char sub[200];
	if (!c || !m)
		return;
	nt = qaws_mesh_3d_get_triangle_count(m);
	uv = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 6 * nt);
	for (k = 0; k < 2; k++)
	{
		double* fx = (double*)malloc(sizeof(double) * 9 * nt);
		double* rgb = (double*)malloc(sizeof(double) * 9 * nt);
		double* s = (double*)malloc(sizeof(double) * nt);
		unsigned int* ft = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nt);
		double a2 = 0, a3 = 0;
		qaws_mesh_3d* flat;
		memset(&d, 0, sizeof d);
		d.auto_seams = 1;
		d.pack = 1;
		d.max_stretch = k ? 1.3 : 0;
		if (qaws_mesh_3d_parameterize(uv, &r, m, NULL, 0, &d) != QAWS_STATUS_OK)
		{
			free(fx);
			free(rgb);
			free(s);
			free(ft);
			continue;
		}
		sprintf(sub, "%u charts, fill %.2f, angles %.3f, areas %.3f", r.charts, (double)r.fill, (double)r.angle_distortion,
			(double)r.area_distortion);
		if (k)
		{
			ms_panel(c, 0, 0, 640, 640, m, -0.6, 0.5, "a closed blob, cut by itself", "the split atlas's (u, v) lines");
			ms_draw(c, m, warm, NULL, 1, NULL, 0);
			ms_df_uvlines(c, m, uv, 0, 1.0 / 48, red);
			ms_df_uvlines(c, m, uv, 1, 1.0 / 48, blue);
		}
		/* each triangle's area scale; the colour: its log over the mean's, +-log 2 at the ends */
		for (i = 0; i < nt; i++)
		{
			unsigned int t[3];
			qaws_vec3 a, b, e;
			double u[3], w[3], nr[3], am;
			qaws_scalar const* q = uv + 6 * i;
			qaws_mesh_3d_get_triangle(t, m, i);
			a = qaws_mesh_3d_get_vertex(m, t[0]);
			b = qaws_mesh_3d_get_vertex(m, t[1]);
			e = qaws_mesh_3d_get_vertex(m, t[2]);
			u[0] = b.x - a.x; u[1] = b.y - a.y; u[2] = b.z - a.z;
			w[0] = e.x - a.x; w[1] = e.y - a.y; w[2] = e.z - a.z;
			nr[0] = u[1] * w[2] - u[2] * w[1];
			nr[1] = u[2] * w[0] - u[0] * w[2];
			nr[2] = u[0] * w[1] - u[1] * w[0];
			am = 0.5 * sqrt(nr[0] * nr[0] + nr[1] * nr[1] + nr[2] * nr[2]);
			s[i] = 0.5 * fabs(((double)q[2] - q[0]) * ((double)q[5] - q[1]) - ((double)q[4] - q[0]) * ((double)q[3] - q[1]));
			a2 += s[i];
			a3 += am;
			s[i] = am > 0 ? s[i] / am : 1;
		}
		for (i = 0; i < nt; i++)
		{
			double col[3];
			ms_colormap(0.5 + 0.5 * log(s[i] * a3 / a2) / log(2.0), col);
			for (z = 0; z < 3; z++)
			{
				fx[9 * i + 3 * z] = uv[6 * i + 2 * z];
				fx[9 * i + 3 * z + 1] = uv[6 * i + 2 * z + 1];
				fx[9 * i + 3 * z + 2] = 0;
				memcpy(rgb + 9 * i + 3 * z, col, sizeof col);
				ft[3 * i + z] = 3 * i + z;
			}
		}
		flat = ms_make(fx, 3 * nt, ft, NULL, nt);
		ms_panel(c, 640 + 640 * (int)k, 0, 640, 640, flat, 0.0, 1.5707, titles[k], sub);
		ms_draw(c, flat, NULL, rgb, 0, wire, 0.3);
		printf("  %s\n", sub);
		qaws_mesh_3d_destroy(flat);
		free(fx);
		free(rgb);
		free(s);
		free(ft);
	}
	free(uv);
	qaws_mesh_3d_destroy(m);
	er_save(c, "showcase/mesh56_atlas.png");
	printf("mesh56_atlas.png\n");
	er_destroy(c);
}
