/*
 * Figure 54: adaptive subdivision (adaptive_angle): a cube of 4 x 4 squares
 * a side, three levels of Loop.
 *   - the cube;
 *   - uniform: every face split each level;
 *   - adaptive at 0.3 radians: only where the surface turns (its edges and
 *     corners), the flat middles left coarse, closed without T-junctions.
 */

static void demo_subdivide_adaptive(void)
{
	static double const tone[3][3] = { { 200, 205, 215 }, { 120, 150, 210 }, { 130, 190, 140 } }, wire[3] = { 60, 66, 80 };
	static char const* const titles[3] = { "a cube", "Loop, uniform", "Loop, adaptive" };
	er_canvas* c = er_create(1350, 640, 3, 252, 252, 250);
	unsigned int k;
	char sub[200];
	if (!c)
		return;
	for (k = 0; k < 3; k++)
	{
		qaws_mesh_3d *cube = ms_gx_cube(4), *o = NULL;
		qaws_mesh_subdivide_desc d;
		memset(&d, 0, sizeof d);
		d.scheme = QAWS_MESH_SUBDIVIDE_LOOP;
		d.levels = 3;
		d.adaptive_angle = (qaws_scalar)(k == 2 ? 0.3 : 0);
		if (k)
			qaws_mesh_3d_subdivide(&o, cube, &d);
		else
			o = cube, cube = NULL;
		if (o)
		{
			sprintf(sub, "%u triangles", qaws_mesh_3d_get_face_count(o));
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
