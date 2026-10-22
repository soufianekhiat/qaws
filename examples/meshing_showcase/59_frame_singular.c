/*
 * Figure 59: hexahedra following a frame field through its singular edges
 * (FRAME_FIELD): the field relaxed by MBO so its singular curves come apart
 * and straight, an integer-grid map with transitions across a cut
 * (CubeCover), its lattice read off as HexEx does, the hexahedra relaxed.
 * The cut faces coloured by their cell's smallest scaled Jacobian.
 *   - a can (radius 0.5, height 1.5) at 0.3: its four singular lines along
 *     its side, a grid turning round it; half;
 *   - a ball at 0.4: its singular graph a cube's; whole, then half;
 *   - a torus (radii 2 and 0.7, its NURBS surface meshed at 0.6) at 0.6: the
 *     cut keeps its disc across the hole.
 */

static void ms_ff_panel(er_canvas* c, int px, qaws_mesh_3d const* s, double h, double yaw, double pitch, char const* title, int (*keep)(double const*),
	double const* rgb, double const* wire)
{
	char sub[200];
	qaws_volume_mesh* v = ms_hexes_by(s, QAWS_HEXES_FRAME_FIELD, h, sub, title);
	if (!v)
		return;
	{
		qaws_volume_quality q;
		qaws_volume_mesh_measure_quality(&q, v);
		sprintf(sub, "%u cells, Jacobian >= %.2f", q.hexahedra, (double)q.min_scaled_jacobian);
	}
	ms_vpanel(c, px, 0, 450, 700, v, yaw, pitch, title, sub);
	ms_hex_cutaway(c, v, keep, rgb, wire);
	qaws_volume_mesh_destroy(v);
}

static void demo_frame_singular(void)
{
	static double const wire[3] = { 40, 45, 60 }, tint[3] = { 178, 192, 212 }, tan[3] = { 232, 196, 120 }, rose[3] = { 225, 180, 180 };
	er_canvas* c = er_create(1800, 700, 3, 252, 252, 250);
	qaws_mesh_3d *can = ms_can(48, 0.5, 1.5), *ball = ms_rm_sphere(32, 16), *torus = NULL;
	qaws_surface* ts = ms_nurbs_torus(2, 0.7);
	if (!c || !can || !ball)
		return;
	if (ts)
		torus = ms_generate(ts, NULL, QAWS_ELEMENTS_TRIANGLES, 0.6, 0, NULL);
	ms_ff_panel(c, 0, can, 0.3, 0.5, 0.6, "a can, half", ms_pc_keep, tan, wire);
	ms_ff_panel(c, 450, ball, 0.4, 0.6, 0.5, "a ball, whole", ms_keep_all, tint, wire);
	ms_ff_panel(c, 900, ball, 0.4, 0.6, 0.5, "a ball, half", ms_pc_keep, tint, wire);
	if (torus)
		ms_ff_panel(c, 1350, torus, 0.6, 0.4, 0.9, "a torus", ms_keep_all, rose, wire);
	qaws_mesh_3d_destroy(can);
	qaws_mesh_3d_destroy(ball);
	qaws_mesh_3d_destroy(torus);
	qaws_surface_destroy(ts);
	er_save(c, "showcase/mesh59_frame_singular.png");
	printf("mesh59_frame_singular.png\n");
	er_destroy(c);
}
