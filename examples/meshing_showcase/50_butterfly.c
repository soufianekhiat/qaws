/*
 * Figure 50: the modified butterfly and sqrt(3) beside Loop
 * (qaws_mesh_3d_subdivide).
 *   - the cage cut into triangles;
 *   - Loop, three levels: approximating, the surface inside the cage, its
 *     vertices (red) off it;
 *   - the modified butterfly, three levels: interpolating, through every
 *     cage vertex;
 *   - sqrt(3), four levels (81 triangles for one, against Loop's 64 in
 *     three): approximating, refining by three each level.
 * Printed: how far each surface's vertex is from its cage vertex, at most.
 */

static void demo_butterfly(void)
{
	static double const cage[3] = { 200, 205, 215 }, loop_rgb[3] = { 120, 150, 210 }, bfly_rgb[3] = { 220, 170, 110 }, s3_rgb[3] = { 140, 190, 150 }, wire[3] = { 60, 66, 80 },
		red[3] = { 200, 50, 40 };
	static qaws_mesh_subdivide_scheme const schemes[3] = { QAWS_MESH_SUBDIVIDE_LOOP, QAWS_MESH_SUBDIVIDE_BUTTERFLY, QAWS_MESH_SUBDIVIDE_SQRT3 };
	static char const* const titles[3] = { "Loop", "the modified butterfly", "sqrt(3)" };
	er_canvas* c = er_create(1800, 640, 3, 252, 252, 250);
	qaws_mesh_3d* t = ms_sd_cage(1);
	double const yaw = 2.6, pitch = 0.85;
	char sub[200];
	unsigned int s, nv, i;
	if (!c || !t)
		return;
	nv = qaws_mesh_3d_get_vertex_count(t);
	sprintf(sub, "%u triangles", qaws_mesh_3d_get_face_count(t));
	ms_panel(c, 0, 0, 450, 640, t, yaw, pitch, "the cage, triangles", sub);
	ms_draw(c, t, cage, NULL, 0, wire, 0.9);
	for (s = 0; s < 3; s++)
	{
		qaws_mesh_subdivide_desc d;
		qaws_mesh_3d* o = NULL;
		double off = 0;
		memset(&d, 0, sizeof d);
		d.scheme = schemes[s];
		d.levels = s == 2 ? 4 : 3;
		qaws_mesh_3d_subdivide(&o, t, &d);
		if (!o)
			continue;
		for (i = 0; i < nv; i++)
		{
			/* (the cage's vertices keep their indices) */
			qaws_vec3 p = qaws_mesh_3d_get_vertex(t, i), q = qaws_mesh_3d_get_vertex(o, i);
			off = fmax(off, sqrt(((double)p.x - q.x) * ((double)p.x - q.x) + ((double)p.y - q.y) * ((double)p.y - q.y) + ((double)p.z - q.z) * ((double)p.z - q.z)));
		}
		sprintf(sub, "the cage's vertices %.3f off", off);
		ms_panel(c, 450 + 450 * (int)s, 0, 450, 640, o, yaw, pitch, titles[s], sub);
		ms_draw(c, o, s == 2 ? s3_rgb : s ? bfly_rgb : loop_rgb, NULL, 1, NULL, 0);
		for (i = 0; i < nv; i++)
		{
			qaws_vec3 p = qaws_mesh_3d_get_vertex(t, i);
			double a[3];
			a[0] = p.x;
			a[1] = p.y;
			a[2] = p.z;
			er_point(c, a, red, 3.5);
		}
		printf("  %s: %s\n", titles[s], sub);
		qaws_mesh_3d_destroy(o);
	}
	qaws_mesh_3d_destroy(t);
	er_save(c, "showcase/mesh50_butterfly.png");
	printf("mesh50_butterfly.png\n");
	er_destroy(c);
}
