/*
 * Figure 48: watertight surfaces from oriented point clouds
 * (qaws_mesh_make_watertight_points): each point a dipole of its share of
 * the surface, the winding number's 1/2 level dual contoured.
 *   - the points: figure 40's ring smoothed (its vertices and normals);
 *   - the surface from them;
 *   - the same points shaken and their normals tilted: still closed;
 *   - a cube sampled on its faces: its edges and corners kept sharp;
 *   - the shaken points without normals: estimated from their neighbours
 *     and oriented along a spanning tree.
 */

static void demo_points(void)
{
	static double const wire[3] = { 90, 95, 110 }, rose[3] = { 225, 180, 180 }, green[3] = { 150, 200, 160 }, blue[3] = { 140, 170, 220 };
	er_canvas* c = er_create(2250, 640, 3, 252, 252, 250);
	qaws_mesh_3d *cage = ms_sd_cage(0), *ring = NULL, *r1 = NULL, *r2 = NULL, *cb = NULL, *r3 = NULL;
	qaws_mesh_subdivide_desc sd;
	qaws_watertight_desc wd;
	qaws_watertight_report rep;
	qaws_scalar *x = NULL, *nr = NULL, *x2 = NULL, *n2 = NULL;
	unsigned int n = 0, i, seed = 21;
	char sub[200];
	double const yaw = 2.6, pitch = 0.85;
	if (!c || !cage)
		return;
	memset(&sd, 0, sizeof sd);
	sd.scheme = QAWS_MESH_SUBDIVIDE_CATMULL_CLARK;
	sd.levels = 3;
	sd.limit = 1;
	qaws_mesh_3d_subdivide(&ring, cage, &sd);
	if (!ring)
		return;
	n = qaws_mesh_3d_get_vertex_count(ring);
	x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n);
	nr = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n);
	x2 = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n);
	n2 = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n);
	qaws_mesh_3d_get_positions(x, ring);
	qaws_mesh_3d_get_vertex_normals(nr, ring);
	{
		/* the points drawn as dots, coloured by their normals */
		double* p = (double*)malloc(sizeof(double) * 3 * n);
		for (i = 0; i < 3 * n; i++) p[i] = (double)x[i];
		er_panel(c, 0, 72, 450, 568, p, n, yaw, pitch, 0.06);
		er_text(c, 14, 12, "points", 3, 30, 30, 30);
		sprintf(sub, "%u, with normals", n);
		er_text(c, 14, 46, sub, 2, 95, 95, 95);
		for (i = 0; i < n; i++)
		{
			double col[3] = { 127 + 110 * (double)nr[3 * i], 127 + 110 * (double)nr[3 * i + 1], 127 + 110 * (double)nr[3 * i + 2] };
			er_point(c, p + 3 * i, col, 1.6);
		}
		free(p);
	}
	memset(&wd, 0, sizeof wd);
	wd.resolution = 64;
	if (qaws_mesh_make_watertight_points(&r1, &rep, x, nr, n, &wd) == QAWS_STATUS_OK && r1)
	{
		sprintf(sub, "%u triangles, closed %d", qaws_mesh_3d_get_face_count(r1), rep.closed);
		ms_panel(c, 450, 0, 450, 640, r1, yaw, pitch, "the surface", sub);
		ms_draw(c, r1, rose, NULL, 1, NULL, 0);
		printf("  points: %s\n", sub);
	}
	for (i = 0; i < 3 * n; i++)
	{
		x2[i] = x[i] + (qaws_scalar)(0.03 * (ms_uniform(&seed) - 0.5));
		n2[i] = nr[i] + (qaws_scalar)(0.5 * (ms_uniform(&seed) - 0.5));
	}
	if (qaws_mesh_make_watertight_points(&r2, &rep, x2, n2, n, &wd) == QAWS_STATUS_OK && r2)
	{
		sprintf(sub, "shaken, normals tilted: closed %d", rep.closed);
		ms_panel(c, 900, 0, 450, 640, r2, yaw, pitch, "noisy points", sub);
		ms_draw(c, r2, green, NULL, 1, NULL, 0);
		printf("  noisy: %s\n", sub);
	}
	wd.resolution = 64;
	if (qaws_mesh_make_watertight_points(&r3, &rep, x2, NULL, n, &wd) == QAWS_STATUS_OK && r3)
	{
		sprintf(sub, "shaken, normals estimated: closed %d", rep.closed);
		ms_panel(c, 1800, 0, 450, 640, r3, yaw, pitch, "no normals", sub);
		ms_draw(c, r3, rose, NULL, 1, NULL, 0);
		printf("  no normals: %s\n", sub);
	}
	{
		/* a cube sampled on its faces */
		unsigned int k, m = 0, j, ii;
		qaws_scalar *cx = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * 6 * 900), *cn = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * 6 * 900);
		for (k = 0; k < 6; k++)
		{
			unsigned int a = k / 2, u = (a + 1) % 3, w = (a + 2) % 3;
			double sg = k % 2 ? 0.5 : -0.5;
			for (ii = 0; ii < 30; ii++)
				for (j = 0; j < 30; j++)
				{
					cx[3 * m + a] = (qaws_scalar)sg;
					cx[3 * m + u] = (qaws_scalar)((ii + ms_uniform(&seed)) / 30 - 0.5);
					cx[3 * m + w] = (qaws_scalar)((j + ms_uniform(&seed)) / 30 - 0.5);
					cn[3 * m] = cn[3 * m + 1] = cn[3 * m + 2] = 0;
					cn[3 * m + a] = (qaws_scalar)(sg > 0 ? 1 : -1);
					m++;
				}
		}
		wd.resolution = 40;
		if (qaws_mesh_make_watertight_points(&cb, &rep, cx, cn, m, &wd) == QAWS_STATUS_OK && cb)
		{
			sprintf(sub, "%u points: %u sharp vertices", m, rep.sharp);
			ms_panel(c, 1350, 0, 450, 640, cb, 0.6, 0.5, "a cube's points", sub);
			ms_draw(c, cb, blue, NULL, 0, wire, 0.3);
			printf("  cube: %s\n", sub);
		}
		free(cx);
		free(cn);
	}
	free(x); free(nr); free(x2); free(n2);
	qaws_mesh_3d_destroy(cage);
	qaws_mesh_3d_destroy(ring);
	qaws_mesh_3d_destroy(r1);
	qaws_mesh_3d_destroy(r2);
	qaws_mesh_3d_destroy(cb);
	qaws_mesh_3d_destroy(r3);
	er_save(c, "showcase/mesh48_points.png");
	printf("mesh48_points.png\n");
	er_destroy(c);
}
