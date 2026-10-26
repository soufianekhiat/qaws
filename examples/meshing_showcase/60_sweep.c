/*
 * Figure 60: swept solids (qaws_mesh_3d_sweep) from motions (qaws_motion):
 * time cut where the trajectories leave their chords, each slab's front
 * faces extruded, all united by the Boolean.
 *   - a cube moving along a curved Hermite path while turning half way round
 *     its vertical axis (the rotation the Cayley map of a Hermite Gibbs
 *     vector);
 *   - a ball along a Hermite curve: a tube;
 *   - a thin box turned a quarter about an axis beside it: an annular sector;
 *   - carving (qaws_mesh_sweep_boolean): a block less a ball tool's sweep along
 *     a curved path dipping into it, then a straight slot (two moves).
 * Printed: slabs, faces, closed, the chords' deviation.
 */

static qaws_mesh_3d* ms_sw_box(double x0, double x1, double y0, double y1, double z0, double z1)
{
	static unsigned int const f[36] = { 0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 2, 3, 7, 2, 7, 6, 1, 2, 6, 1, 6, 5, 0, 4, 7, 0, 7, 3 };
	double x[24];
	unsigned int i;
	for (i = 0; i < 8; i++)
	{
		x[3 * i] = (i == 1 || i == 2 || i == 5 || i == 6) ? x1 : x0;
		x[3 * i + 1] = (i == 2 || i == 3 || i == 6 || i == 7) ? y1 : y0;
		x[3 * i + 2] = i >= 4 ? z1 : z0;
	}
	return ms_make(x, 8, f, NULL, 12);
}

static void ms_sw_panel(er_canvas* c, int px, qaws_mesh_3d const* body, qaws_motion const* m, double tol, double yaw, double pitch, char const* title,
	double const* rgb)
{
	static double const wire[3] = { 60, 66, 80 };
	qaws_mesh_3d* o = NULL;
	qaws_sweep_desc d;
	qaws_sweep_report r;
	char sub[200];
	memset(&d, 0, sizeof d);
	d.tolerance = (qaws_scalar)tol;
	if (qaws_mesh_3d_sweep(&o, &r, body, m, &d) != QAWS_STATUS_OK || !o)
		return;
	sprintf(sub, "%u slabs, %u faces, closed %d", r.slabs, qaws_mesh_3d_get_face_count(o), r.closed);
	ms_panel(c, px, 0, 600, 640, o, yaw, pitch, title, sub);
	ms_draw(c, o, rgb, NULL, 0, wire, 0.4);
	printf("  %s: %s\n", title, sub);
	qaws_mesh_3d_destroy(o);
}

static void demo_sweep(void)
{
	static double const tan[3] = { 232, 196, 120 }, blue[3] = { 140, 170, 220 }, green[3] = { 150, 200, 150 };
	er_canvas* c = er_create(2400, 640, 3, 252, 252, 250);
	qaws_motion_pose s, e;
	qaws_motion* m = NULL;
	if (!c)
		return;
	{
		/* a cube on a curved path, turning half way round z */
		qaws_mesh_3d* cube = ms_sw_box(-0.25, 0.25, -0.25, 0.25, -0.25, 0.25);
		memset(&s, 0, sizeof s);
		memset(&e, 0, sizeof e);
		s.rotation[0] = 1;
		e.rotation[0] = (qaws_scalar)cos(0.4 * 3.14159265358979323846);
		e.rotation[3] = (qaws_scalar)sin(0.4 * 3.14159265358979323846);
		e.position.x = 3; e.position.y = 0.5;
		s.velocity.x = 2; s.velocity.y = 3;
		e.velocity.x = 2; e.velocity.y = -2;
		s.angular_velocity.z = e.angular_velocity.z = (qaws_scalar)(0.8 * 3.14159265358979323846);
		qaws_motion_create_poses(&m, &s, &e, 0);
		if (m)
			ms_sw_panel(c, 0, cube, m, 0.005, -0.3, 0.9, "a cube along a curve, turning", tan);
		qaws_motion_destroy(m);
		m = NULL;
		qaws_mesh_3d_destroy(cube);
	}
	{
		/* a ball along a Hermite curve rising */
		qaws_mesh_3d* ball = ms_rm_sphere(24, 12);
		unsigned int nv = ball ? qaws_mesh_3d_get_vertex_count(ball) : 0, i;
		if (ball)
		{
			qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
			qaws_mesh_3d_get_positions(x, ball);
			for (i = 0; i < 3 * nv; i++) x[i] = (qaws_scalar)(0.25 * (double)x[i]);
			qaws_mesh_3d_set_positions(ball, x);
			free(x);
		}
		memset(&s, 0, sizeof s);
		memset(&e, 0, sizeof e);
		s.rotation[0] = e.rotation[0] = 1;
		e.position.x = 3; e.position.z = 1;
		s.velocity.x = 1; s.velocity.y = 4;
		e.velocity.x = 1; e.velocity.y = -4; e.velocity.z = 2;
		qaws_motion_create_poses(&m, &s, &e, 0);
		if (m && ball)
			ms_sw_panel(c, 600, ball, m, 0.005, -0.3, 0.8, "a ball along a Hermite curve", blue);
		qaws_motion_destroy(m);
		m = NULL;
		qaws_mesh_3d_destroy(ball);
	}
	{
		/* a thin box turned a quarter about the z axis beside it */
		qaws_mesh_3d* box = ms_sw_box(1, 2, -0.1, 0.1, 0, 0.6);
		memset(&s, 0, sizeof s);
		memset(&e, 0, sizeof e);
		s.rotation[0] = 1;
		e.rotation[0] = (qaws_scalar)cos(3.14159265358979323846 / 4);
		e.rotation[3] = (qaws_scalar)sin(3.14159265358979323846 / 4);
		s.angular_velocity.z = e.angular_velocity.z = (qaws_scalar)(3.14159265358979323846 / 2);
		qaws_motion_create_poses(&m, &s, &e, 0);
		if (m)
			ms_sw_panel(c, 1200, box, m, 0.002, -0.5, 0.9, "a thin box turned a quarter", green);
		qaws_motion_destroy(m);
		qaws_mesh_3d_destroy(box);
	}
	{
		/* carving: a block less a ball tool along a curve dipping into it, then a straight slot */
		static double const wood[3] = { 214, 178, 130 }, wire[3] = { 60, 66, 80 };
		qaws_mesh_3d *block = ms_sw_box(0, 4, 0, 2, 0, 1), *ball = ms_rm_sphere(24, 12), *o = NULL;
		qaws_motion* mv[2] = { NULL, NULL };
		qaws_motion const* cm[2];
		qaws_sweep_desc d;
		qaws_sweep_boolean_report br;
		unsigned int nv = ball ? qaws_mesh_3d_get_vertex_count(ball) : 0, i;
		char sub[200];
		if (ball)
		{
			qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
			qaws_mesh_3d_get_positions(x, ball);
			for (i = 0; i < 3 * nv; i++) x[i] = (qaws_scalar)(0.3 * (double)x[i]);
			qaws_mesh_3d_set_positions(ball, x);
			free(x);
		}
		memset(&s, 0, sizeof s);
		memset(&e, 0, sizeof e);
		s.rotation[0] = e.rotation[0] = 1;
		s.position.x = -0.5; s.position.y = 0.6; s.position.z = 1.2;
		e.position.x = 4.5; e.position.y = 0.6; e.position.z = 1.2;
		s.velocity.x = 5; s.velocity.y = 3; s.velocity.z = -2;
		e.velocity.x = 5; e.velocity.y = -3; e.velocity.z = 2;
		qaws_motion_create_poses(&mv[0], &s, &e, 0);
		memset(&s, 0, sizeof s);
		memset(&e, 0, sizeof e);
		s.rotation[0] = e.rotation[0] = 1;
		s.position.x = 3; s.position.y = -0.5; s.position.z = 1;
		e.position.x = 3; e.position.y = 2.5; e.position.z = 1;
		s.velocity.y = e.velocity.y = 3;
		qaws_motion_create_poses(&mv[1], &s, &e, 0);
		cm[0] = mv[0];
		cm[1] = mv[1];
		memset(&d, 0, sizeof d);
		d.tolerance = (qaws_scalar)0.005;
		if (block && ball && mv[0] && mv[1] && qaws_mesh_sweep_boolean(&o, &br, ball, cm, block, 2, QAWS_BOOLEAN_DIFFERENCE, &d) == QAWS_STATUS_OK && o)
		{
			sprintf(sub, "%u moves, %u slabs, removed %.4f of %.1f, closed %d", br.moves, br.slabs, (double)br.changed_volume, (double)br.stock_volume, br.closed);
			ms_panel(c, 1800, 0, 600, 640, o, -0.4, 0.8, "carving: two moves", sub);
			ms_draw(c, o, wood, NULL, 0, wire, 0.3);
			printf("  carving: %s\n", sub);
		}
		qaws_mesh_3d_destroy(o);
		qaws_mesh_3d_destroy(block);
		qaws_mesh_3d_destroy(ball);
		qaws_motion_destroy(mv[0]);
		qaws_motion_destroy(mv[1]);
	}
	er_save(c, "showcase/mesh60_sweep.png");
	printf("mesh60_sweep.png\n");
	er_destroy(c);
}
