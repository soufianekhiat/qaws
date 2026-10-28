/*
 * Figure 67: features through sweeps, edge values through the connector.
 *   - a featured cube swept along a cubic Hermite curve while turning a
 *     quarter: its features are the solid's creases, the curves its corners
 *     trace included (red);
 *   - a block carved by a featured box tool along a curve: the slot's
 *     edges and the tool's traced corners are features;
 *   - two tubes with a value per edge (one hue per edge) joined by the
 *     connector: each tube's edges keep their values, the bridge's are 0
 *     (thin);
 *   - a 16-sided cylinder (its rims features) swept along a curve: the
 *     lines its rims trace where they run along the motion (their envelope,
 *     traced by rim vertices that are no corners) are features too.
 * Printed: feature edges, faces; for the cylinder the edges a feature
 * without turning by more than 45 degrees or the reverse.
 */

/* m's edges a feature without its faces turning by more than 45 degrees, or the reverse */
static unsigned int ms_sf_off(qaws_mesh_3d const* m)
{
	unsigned int ne = qaws_mesh_3d_get_edge_count(m), nt = qaws_mesh_3d_get_triangle_count(m), nfe = qaws_mesh_3d_get_feature_edge_count(m), k, f, j, off = 0;
	unsigned int *e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (ne + 1)), *fe = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (nfe + 1));
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (nt + 1));
	double* x = ms_positions(m);
	double* n = (double*)malloc(sizeof(double) * 3 * (nt + 1));
	qaws_mesh_3d_get_edges(e, m);
	qaws_mesh_3d_get_feature_edges(fe, m);
	qaws_mesh_3d_get_triangles(t, m);
	for (f = 0; f < nt; f++)
	{
		double const *a = x + 3 * t[3 * f], *b = x + 3 * t[3 * f + 1], *c = x + 3 * t[3 * f + 2];
		double u[3] = { b[0] - a[0], b[1] - a[1], b[2] - a[2] }, w[3] = { c[0] - a[0], c[1] - a[1], c[2] - a[2] }, l;
		n[3 * f] = u[1] * w[2] - u[2] * w[1]; n[3 * f + 1] = u[2] * w[0] - u[0] * w[2]; n[3 * f + 2] = u[0] * w[1] - u[1] * w[0];
		l = sqrt(n[3 * f] * n[3 * f] + n[3 * f + 1] * n[3 * f + 1] + n[3 * f + 2] * n[3 * f + 2]);
		for (j = 0; j < 3; j++) n[3 * f + j] /= l > 0 ? l : 1;
	}
	for (k = 0; k < ne; k++)
	{
		unsigned int g[2] = { 0, 0 }, side = 0;
		int feat = 0, crease;
		for (j = 0; j < nfe; j++)
			feat |= (fe[2 * j] == e[2 * k] && fe[2 * j + 1] == e[2 * k + 1]) || (fe[2 * j] == e[2 * k + 1] && fe[2 * j + 1] == e[2 * k]);
		for (f = 0; f < nt && side < 2; f++)
			for (j = 0; j < 3; j++)
				if ((t[3 * f + j] == e[2 * k] && t[3 * f + (j + 1) % 3] == e[2 * k + 1]) || (t[3 * f + j] == e[2 * k + 1] && t[3 * f + (j + 1) % 3] == e[2 * k]))
					g[side++] = f;
		crease = side == 2 && n[3 * g[0]] * n[3 * g[1]] + n[3 * g[0] + 1] * n[3 * g[1] + 1] + n[3 * g[0] + 2] * n[3 * g[1] + 2] < cos(MS_PI / 4);
		off += (unsigned int)(feat != crease);
	}
	free(e); free(fe); free(t); free(x); free(n);
	return off;
}

static void demo_sweep_features(void)
{
	static double const face[3] = { 226, 222, 212 }, red[3] = { 210, 40, 40 }, wire[3] = { 150, 150, 150 };
	er_canvas* c = er_create(2400, 640, 3, 252, 252, 250);
	char sub[200];
	if (!c)
		return;
	{
		/* a featured cube along a curve, turning */
		qaws_mesh_3d *cube = ms_sw_box(-0.3, 0.3, -0.3, 0.3, -0.3, 0.3), *o = NULL;
		qaws_motion* m = ms_at3_poses(1);
		qaws_sweep_desc d;
		qaws_mesh_3d_detect_features(cube, (qaws_scalar)0.5, 0, 0);
		memset(&d, 0, sizeof d);
		d.tolerance = (qaws_scalar)0.004;
		d.simplify = 1;   /* (the flat regions and straight creases without the slabs' extra vertices) */
		if (m && qaws_mesh_3d_sweep(&o, NULL, cube, m, &d) == QAWS_STATUS_OK && o)
		{
			sprintf(sub, "%u feature edges, %u faces", qaws_mesh_3d_get_feature_edge_count(o), qaws_mesh_3d_get_face_count(o));
			ms_panel(c, 0, 0, 600, 640, o, -0.35, 0.95, "a cube's swept creases", sub);
			ms_draw(c, o, face, NULL, 0, NULL, 0);
			ms_at2_features(c, o, red, 2.5);
			printf("  cube: %s\n", sub);
		}
		qaws_motion_destroy(m);
		qaws_mesh_3d_destroy(cube);
		qaws_mesh_3d_destroy(o);
	}
	{
		/* a block carved by a featured box tool */
		qaws_mesh_3d *tool = ms_sw_box(-0.3, 0.3, -0.3, 0.3, -0.3, 0.3), *block = ms_sw_box(0, 4, 0, 2, 0, 1), *o = NULL;
		qaws_motion* m = ms_at3_poses(2);
		qaws_motion const* mv[1];
		qaws_sweep_desc d;
		qaws_mesh_3d_detect_features(tool, (qaws_scalar)0.5, 0, 0);
		qaws_mesh_3d_detect_features(block, (qaws_scalar)0.5, 0, 0);
		mv[0] = m;
		memset(&d, 0, sizeof d);
		d.tolerance = (qaws_scalar)0.004;
		d.simplify = 1;   /* (the flat regions and straight creases without the slabs' extra vertices) */
		if (m && qaws_mesh_sweep_boolean(&o, NULL, tool, mv, block, 1, QAWS_BOOLEAN_DIFFERENCE, &d) == QAWS_STATUS_OK && o)
		{
			sprintf(sub, "%u feature edges, %u faces", qaws_mesh_3d_get_feature_edge_count(o), qaws_mesh_3d_get_face_count(o));
			ms_panel(c, 600, 0, 600, 640, o, -0.35, 0.95, "carving: the slot's creases", sub);
			ms_draw(c, o, face, NULL, 0, NULL, 0);
			ms_at2_features(c, o, red, 2.5);
			printf("  carving: %s\n", sub);
		}
		qaws_motion_destroy(m);
		qaws_mesh_3d_destroy(tool); qaws_mesh_3d_destroy(block); qaws_mesh_3d_destroy(o);
	}
	{
		/* two tubes with edge values joined */
		qaws_mesh_3d* t[2] = { NULL, NULL };
		qaws_mesh_3d* o = NULL;
		qaws_connect_desc cd;
		unsigned int s;
		for (s = 0; s < 2; s++)
		{
			enum { N = 16 };
			double x[3 * 2 * N], w[3 * 2 * N + 8];
			unsigned int f[6 * N], i, k = 0, ne;
			for (i = 0; i < 2 * N; i++)
			{
				double a = 2 * MS_PI * (i % N) / N;
				x[3 * i] = 0.6 * cos(a); x[3 * i + 1] = 0.6 * sin(a); x[3 * i + 2] = (s ? 2.2 : 0) + (i < N ? 0 : 1);
			}
			for (i = 0; i < N; i++)
			{
				unsigned int a = i, b = (i + 1) % N;
				f[k++] = a; f[k++] = b; f[k++] = b + N;
				f[k++] = a; f[k++] = b + N; f[k++] = a + N;
			}
			t[s] = ms_make(x, 2 * N, f, NULL, 2 * N);
			ne = qaws_mesh_3d_get_edge_count(t[s]);
			for (i = 0; i < ne; i++) w[i] = 1 + (i + 7 * s) % 12;
			ms_at_add(t[s], "w", QAWS_MESH_EDGE, 1, QAWS_ATTRIBUTE_NEAREST, w, ne);
		}
		memset(&cd, 0, sizeof cd);
		cd.quads = 1;
		cd.rows = 4;
		if (qaws_mesh_3d_connect(&o, NULL, t[0], 16, t[1], 0, &cd) == QAWS_STATUS_OK && o)
		{
			unsigned int valued = 0;
			ms_panel(c, 1200, 0, 600, 640, o, -0.5, 0.35, "edge values, joined", NULL);
			ms_sr_edges(c, o, 12, &valued);
			sprintf(sub, "%u edges, %u keep their tube's value", qaws_mesh_3d_get_edge_count(o), valued);
			er_text(c, 1200 + 14, 46, sub, 2, 95, 95, 95);
			printf("  joined: %s\n", sub);
		}
		(void)wire;
		qaws_mesh_3d_destroy(t[0]); qaws_mesh_3d_destroy(t[1]); qaws_mesh_3d_destroy(o);
	}
	{
		/* a 16-sided cylinder along a curve, not turning */
		enum { N = 16 };
		double x[3 * 2 * N];
		unsigned int f[4 * N + 2 * N], fs[N + 2], i, k = 0;
		qaws_mesh_3d *cyl, *o = NULL;
		qaws_motion_pose s, e;
		qaws_motion* m = NULL;
		qaws_sweep_desc d;
		for (i = 0; i < 2 * N; i++)
		{
			double a = 2 * MS_PI * (i % N) / N;
			x[3 * i] = 0.4 * cos(a); x[3 * i + 1] = 0.4 * sin(a); x[3 * i + 2] = i < N ? -0.3 : 0.3;
		}
		for (i = 0; i < N; i++)
		{
			unsigned int a = i, b = (i + 1) % N;
			f[k++] = a; f[k++] = b; f[k++] = b + N; f[k++] = a + N;
			fs[i] = 4;
		}
		for (i = 0; i < N; i++)
			f[k++] = N - 1 - i;
		for (i = 0; i < N; i++)
			f[k++] = N + i;
		fs[N] = fs[N + 1] = N;
		cyl = ms_make(x, 2 * N, f, fs, N + 2);
		qaws_mesh_3d_detect_features(cyl, (qaws_scalar)0.5, 0, 0);
		memset(&s, 0, sizeof s);
		memset(&e, 0, sizeof e);
		s.rotation[0] = e.rotation[0] = 1;
		s.position.x = -0.5; s.position.y = 1; s.position.z = 1.05;
		e.position.x = 4.5; e.position.y = 1; e.position.z = 1.05;
		s.velocity.x = e.velocity.x = 5; s.velocity.y = 2.4; e.velocity.y = -2.4;
		qaws_motion_create_poses(&m, &s, &e, 0);
		memset(&d, 0, sizeof d);
		d.tolerance = (qaws_scalar)0.004;
		d.simplify = 1;
		if (cyl && m && qaws_mesh_3d_sweep(&o, NULL, cyl, m, &d) == QAWS_STATUS_OK && o)
		{
			unsigned int off = ms_sf_off(o);
			sprintf(sub, "%u feature edges, %u faces, %u off", qaws_mesh_3d_get_feature_edge_count(o), qaws_mesh_3d_get_face_count(o), off);
			ms_panel(c, 1800, 0, 600, 640, o, -0.35, 0.95, "a cylinder's rims' envelope", sub);
			ms_draw(c, o, face, NULL, 0, NULL, 0);
			ms_at2_features(c, o, red, 2.5);
			printf("  cylinder: %s\n", sub);
		}
		qaws_motion_destroy(m);
		qaws_mesh_3d_destroy(cyl);
		qaws_mesh_3d_destroy(o);
	}
	er_save(c, "showcase/mesh67_sweep_features.png");
	printf("mesh67_sweep_features.png\n");
	er_destroy(c);
}
