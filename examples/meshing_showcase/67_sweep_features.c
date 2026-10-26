/*
 * Figure 67: features through sweeps, edge values through the connector.
 *   - a featured cube swept along a cubic Hermite curve while turning a
 *     quarter: its features are the solid's creases, the curves its corners
 *     trace included (red);
 *   - a block carved by a featured box tool along a curve: the slot's
 *     edges and the tool's traced corners are features;
 *   - two tubes with a value per edge (one hue per edge) joined by the
 *     connector: each tube's edges keep their values, the bridge's are 0
 *     (thin).
 * Printed: feature edges, faces.
 */

static void demo_sweep_features(void)
{
	static double const face[3] = { 226, 222, 212 }, red[3] = { 210, 40, 40 }, wire[3] = { 150, 150, 150 };
	er_canvas* c = er_create(1800, 640, 3, 252, 252, 250);
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
	er_save(c, "showcase/mesh67_sweep_features.png");
	printf("mesh67_sweep_features.png\n");
	er_destroy(c);
}
