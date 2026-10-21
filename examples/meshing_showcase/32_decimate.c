/*
 * Figure 32: decimation by quadric error (qaws_mesh_3d_decimate).
 *   - a rippled torus decimated to 20%, 5% and 1.5% of its triangles: the
 *     triangles stay where it bends, the certified deviation from the input
 *     printed under each;
 *   - a cube in fine grids with its creases as features: down to a few
 *     triangles a side, its edges and corners kept.
 */

static qaws_mesh_3d* ms_dc_torus(void)
{
	enum { NU = 160, NV = 64 };
	double* x = (double*)malloc(sizeof(double) * 3 * NU * NV);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 6 * NU * NV), i, j;
	qaws_mesh_3d* m;
	ms_torus_grid(NU, NV, x, t);
	for (i = 0; i < NU; i++)
		for (j = 0; j < NV; j++)
		{
			/* ripples round the ring, on the outer half only */
			double u = 2 * MS_PI * i / NU, v = 2 * MS_PI * j / NV, a = 0.06 * sin(9 * u) * (cos(v) > 0 ? cos(v) : 0);
			double* p = x + 3 * (i * NV + j);
			p[0] += a * cos(v) * cos(u);
			p[1] += a * cos(v) * sin(u);
			p[2] += a * sin(v);
		}
	m = ms_make(x, NU * NV, t, NULL, 2 * NU * NV);
	free(x);
	free(t);
	return m;
}

static void demo_decimate(void)
{
	static double const col[3] = { 180, 190, 214 }, warm[3] = { 226, 196, 150 }, wire[3] = { 60, 66, 80 };
	static double const share[4] = { 1, 0.2, 0.05, 0.015 };
	er_canvas* c = er_create(2500, 720, 3, 252, 252, 250);
	qaws_mesh_3d* frame = ms_dc_torus();
	char sub[200], title[80];
	unsigned int k;
	if (!c)
		return;
	for (k = 0; k < 4; k++)
	{
		qaws_mesh_3d* m = ms_dc_torus();
		qaws_decimate_desc d;
		qaws_decimate_report r;
		memset(&d, 0, sizeof d);
		memset(&r, 0, sizeof r);
		d.target_faces = (unsigned int)(share[k] * qaws_mesh_3d_get_face_count(m));
		d.flags = QAWS_DECIMATE_MEASURE;
		if (k)
		{
			qaws_mesh_3d_decimate(&r, m, &d);
			sprintf(sub, "%u triangles, off by %.4f, certified", r.faces_after, (double)r.deviation.hausdorff[1]);
		}
		else
			sprintf(sub, "%u triangles", qaws_mesh_3d_get_face_count(m));
		sprintf(title, k ? "%.1f%% of it" : "a rippled torus", 100 * share[k]);
		ms_panel(c, 500 * (int)k, 0, 500, 720, frame, -0.3, 0.9, title, sub);
		ms_draw(c, m, col, NULL, 0, wire, k ? 0.6 : 0.2);
		printf("  torus, %s: %s\n", title, sub);
		qaws_mesh_3d_destroy(m);
	}
	{
		/* the cube: 8 x 8 grids a side, creases detected */
		unsigned int n = 8, s = n + 1, nv = 0, nt = 0, f, i, j, q;
		unsigned int* id = (unsigned int*)malloc(sizeof(unsigned int) * s * s * s);
		double* x = (double*)malloc(sizeof(double) * 3 * s * s * s);
		unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 36 * n * n);
		qaws_mesh_3d* m;
		qaws_decimate_desc d;
		qaws_decimate_report r;
		for (i = 0; i < s * s * s; i++)
			id[i] = 0xFFFFFFFFu;
		for (f = 0; f < 6; f++)
		{
			unsigned int ax = f / 2, u = (ax + 1) % 3, w = (ax + 2) % 3, hi = f % 2;
			for (j = 0; j < n; j++)
				for (i = 0; i < n; i++)
				{
					unsigned int qv[4], cc;
					for (cc = 0; cc < 4; cc++)
					{
						unsigned int g[3], key;
						g[ax] = hi ? n : 0;
						g[u] = i + (cc == 1 || cc == 2);
						g[w] = j + (cc >= 2);
						key = (g[0] * s + g[1]) * s + g[2];
						if (id[key] == 0xFFFFFFFFu)
						{
							for (q = 0; q < 3; q++)
								x[3 * nv + q] = (double)g[q] / n - 0.5;
							id[key] = nv++;
						}
						qv[cc] = id[key];
					}
					t[3 * nt] = qv[0]; t[3 * nt + 1] = hi ? qv[1] : qv[2]; t[3 * nt + 2] = hi ? qv[2] : qv[1]; nt++;
					t[3 * nt] = qv[0]; t[3 * nt + 1] = hi ? qv[2] : qv[3]; t[3 * nt + 2] = hi ? qv[3] : qv[2]; nt++;
				}
		}
		m = ms_make(x, nv, t, NULL, nt);
		qaws_mesh_3d_detect_features(m, (qaws_scalar)0.5, (qaws_scalar)0.5, 0);
		memset(&d, 0, sizeof d);
		d.target_faces = 12;
		qaws_mesh_3d_decimate(&r, m, &d);
		sprintf(sub, "%u -> %u triangles, creases kept", r.faces_before, r.faces_after);
		ms_panel(c, 2000, 0, 500, 720, m, 0.6, 0.55, "a fine cube, decimated", sub);
		ms_draw(c, m, warm, NULL, 0, wire, 0.8);
		printf("  cube: %s\n", sub);
		qaws_mesh_3d_destroy(m);
		free(id); free(x); free(t);
	}
	qaws_mesh_3d_destroy(frame);
	er_save(c, "showcase/mesh32_decimate.png");
	printf("mesh32_decimate.png\n");
	er_destroy(c);
}
