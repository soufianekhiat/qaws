/*
 * Figure 52: rounded offsets (QAWS_MESH_OFFSET_ROUND), the U of figure 39.
 *   - moved out by 0.12, mitred: the faces' planes moved, corners sharp;
 *   - rounded out by 0.12: the points within 0.12, cylinders about its
 *     convex edges and spheres at its corners, its concave edges sharp;
 *   - rounded in by 0.08: the points 0.08 deep, its concave edges and
 *     corners now rounded, its convex ones sharp.
 * Printed: the distances of each result's vertices to the U.
 */

static void demo_offset_round(void)
{
	static double const tone[3][3] = { { 130, 190, 140 }, { 120, 150, 210 }, { 230, 160, 100 } }, wire[3] = { 60, 66, 80 };
	static qaws_mesh_offset_method const how[3] = { QAWS_MESH_OFFSET_MOVE, QAWS_MESH_OFFSET_ROUND, QAWS_MESH_OFFSET_ROUND };
	static double const dist[3] = { 0.12, 0.12, -0.08 };
	static char const* const titles[3] = { "moved, mitred", "rounded out", "rounded in" };
	er_canvas* c = er_create(1350, 640, 3, 252, 252, 250);
	qaws_mesh_3d* u = ms_of_u();
	unsigned int k;
	char sub[200];
	if (!c || !u)
		return;
	for (k = 0; k < 3; k++)
	{
		qaws_mesh_offset_desc d;
		qaws_mesh_offset_report r;
		qaws_mesh_3d* o = NULL;
		clock_t t0 = clock();
		memset(&d, 0, sizeof d);
		d.kind = QAWS_MESH_OFFSET_SURFACE;
		d.method = how[k];
		d.distance = (qaws_scalar)dist[k];
		qaws_mesh_3d_offset(&o, &r, u, &d);
		if (!o)
			continue;
		sprintf(sub, "%+.2f: %.3f to %.3f away, %.0f ms", dist[k], (double)r.min_distance, (double)r.max_distance,
			1000.0 * (double)(clock() - t0) / CLOCKS_PER_SEC);
		ms_panel(c, 450 * (int)k, 0, 450, 640, k == 2 ? u : o, 0.6, 0.55, titles[k], sub);
		ms_draw(c, o, tone[k], NULL, k > 0, k ? NULL : wire, 0.6);
		if (k == 2)
			ms_draw(c, u, NULL, NULL, 0, wire, 0.5);
		printf("  %s: %s\n", titles[k], sub);
		qaws_mesh_3d_destroy(o);
	}
	qaws_mesh_3d_destroy(u);
	er_save(c, "showcase/mesh52_offset_round.png");
	printf("mesh52_offset_round.png\n");
	er_destroy(c);
}
