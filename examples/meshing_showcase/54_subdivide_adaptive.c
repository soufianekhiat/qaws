/*
 * Figure 54: adaptive subdivision (adaptive_angle): a cube of 4 x 4 squares
 * a side, three levels.
 *   - the cube;
 *   - Loop, uniform: every face split each level;
 *   - Loop, adaptive at 0.3 radians: only where the surface turns (its edges
 *     and corners), the flat middles left coarse, closed without
 *     T-junctions;
 *   - Catmull-Clark on the cube's quads, uniform and adaptive at 0.3: the
 *     split faces' quads, transition polygons (a face keeping its split
 *     edges' points) between them and the coarse middles.
 */

/* the cube of n x n quads a side: ms_gx_cube's triangles paired back */
static qaws_mesh_3d* ms_quad_cube(unsigned int n)
{
	qaws_mesh_3d* t = ms_gx_cube(n);
	qaws_mesh_3d* q = NULL;
	unsigned int nv, nt, i, *tv, *fv, *fs;
	qaws_scalar* xs;
	double* x;
	if (!t)
		return NULL;
	nv = qaws_mesh_3d_get_vertex_count(t);
	nt = qaws_mesh_3d_get_triangle_count(t);
	tv = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nt);
	fv = (unsigned int*)malloc(sizeof(unsigned int) * 2 * nt);
	fs = (unsigned int*)malloc(sizeof(unsigned int) * nt / 2);
	xs = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	x = (double*)malloc(sizeof(double) * 3 * nv);
	qaws_mesh_3d_get_triangles(tv, t);
	qaws_mesh_3d_get_positions(xs, t);
	for (i = 0; i < 3 * nv; i++) x[i] = (double)xs[i];
	for (i = 0; i < nt / 2; i++)
	{
		/* (triangles q0 q1 q2 and q0 q2 q3) */
		fv[4 * i] = tv[6 * i]; fv[4 * i + 1] = tv[6 * i + 1]; fv[4 * i + 2] = tv[6 * i + 2]; fv[4 * i + 3] = tv[6 * i + 5];
		fs[i] = 4;
	}
	q = ms_make(x, nv, fv, fs, nt / 2);
	qaws_mesh_3d_destroy(t);
	free(tv); free(fv); free(fs); free(xs); free(x);
	return q;
}

static void demo_subdivide_adaptive(void)
{
	static double const tone[5][3] = { { 200, 205, 215 }, { 120, 150, 210 }, { 130, 190, 140 }, { 210, 160, 110 }, { 200, 140, 170 } }, wire[3] = { 60, 66, 80 };
	static char const* const titles[5] = { "a cube", "Loop, uniform", "Loop, adaptive", "Catmull-Clark, uniform", "Catmull-Clark, adaptive" };
	er_canvas* c = er_create(2250, 640, 3, 252, 252, 250);
	unsigned int k;
	char sub[200];
	if (!c)
		return;
	for (k = 0; k < 5; k++)
	{
		qaws_mesh_3d *cube = k >= 3 ? ms_quad_cube(4) : ms_gx_cube(4), *o = NULL;
		qaws_mesh_subdivide_desc d;
		memset(&d, 0, sizeof d);
		d.scheme = k >= 3 ? QAWS_MESH_SUBDIVIDE_CATMULL_CLARK : QAWS_MESH_SUBDIVIDE_LOOP;
		d.levels = 3;
		d.adaptive_angle = (qaws_scalar)(k == 2 || k == 4 ? 0.3 : 0);
		if (k)
			qaws_mesh_3d_subdivide(&o, cube, &d);
		else
			o = cube, cube = NULL;
		if (o)
		{
			sprintf(sub, "%u %s", qaws_mesh_3d_get_face_count(o), k >= 3 ? "faces" : "triangles");
			ms_panel(c, 450 * (int)k, 0, 450, 640, o, 0.6, 0.55, titles[k], sub);
			ms_draw(c, o, tone[k], NULL, k > 0, wire, 0.5);
			printf("  %s: %s\n", titles[k], sub);
		}
		qaws_mesh_3d_destroy(o);
		qaws_mesh_3d_destroy(cube);
	}
	er_save(c, "showcase/mesh54_subdivide_adaptive.png");
	printf("mesh54_subdivide_adaptive.png\n");
	er_destroy(c);
}
