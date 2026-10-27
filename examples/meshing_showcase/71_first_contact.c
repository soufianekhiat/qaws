/*
 * Figure 71: continuous collision detection (qaws_mesh_3d_first_contact).
 *   - a ball falling on a slab: its positions every tenth of the motion
 *     (faint), the pose of first contact (solid) and its point (red);
 *   - a cube shot through a thin wall: the frames at 0, 0.15 and 0.3 never
 *     meet the wall (a check by frames passes through it); the first contact,
 *     at 0.19, is found;
 *   - a rod turning a quarter toward a wall: its contact pose, the corner on
 *     the wall.
 * Printed: the contact times.
 */

/* a copy of m moved by the motion to t */
static qaws_mesh_3d* ms_fc_at(qaws_mesh_3d const* m, qaws_motion const* mo, double t)
{
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m);
	qaws_scalar *x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv), *y = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	qaws_mesh_3d* c = NULL;
	qaws_mesh_3d_get_positions(x, m);
	qaws_motion_transform_points(y, NULL, mo, (qaws_scalar)t, x, nv);
	qaws_mesh_3d_clone(&c, m);
	if (c)
		qaws_mesh_3d_set_positions(c, y);
	free(x); free(y);
	return c;
}

/* a panel: the obstacle, the moving mesh at the given times (faint), at the contact (solid), the contact point */
static void ms_fc_panel(er_canvas* c, int px, qaws_mesh_3d const* m, qaws_motion const* mo, qaws_mesh_3d const* obstacle, double const* ts, unsigned int nts,
	char const* title, double yaw, double pitch)
{
	static double const grey[3] = { 200, 196, 188 }, ghost[3] = { 170, 190, 230 }, solid[3] = { 70, 120, 220 }, red[3] = { 220, 40, 40 }, wire[3] = { 110, 110, 120 };
	qaws_contact_report r;
	qaws_sweep_desc d;
	unsigned int i, n = 0, nv = qaws_mesh_3d_get_vertex_count(m), no = qaws_mesh_3d_get_vertex_count(obstacle);
	double* fit = (double*)malloc(sizeof(double) * 3 * ((nts + 1) * nv + no + 1));
	char sub[200];
	memset(&d, 0, sizeof d);
	d.tolerance = (qaws_scalar)1e-4;
	qaws_mesh_3d_first_contact(&r, m, mo, obstacle, &d);
	{
		double* x = ms_positions(obstacle);
		memcpy(fit, x, sizeof(double) * 3 * no);
		n = no;
		free(x);
	}
	for (i = 0; i <= nts; i++)
	{
		qaws_mesh_3d* g = ms_fc_at(m, mo, i < nts ? ts[i] : (double)r.time);
		double* x = g ? ms_positions(g) : NULL;
		if (x)
		{
			memcpy(fit + 3 * n, x, sizeof(double) * 3 * nv);
			n += nv;
		}
		free(x);
		qaws_mesh_3d_destroy(g);
	}
	if (r.contact)
		sprintf(sub, "first contact in [%.5f, %.5f]", (double)r.time, (double)r.time_upper);
	else
		sprintf(sub, "no contact");
	er_panel(c, px, 72, 600, 640 - 72, fit, n, yaw, pitch, 0.06);
	er_text(c, px + 14, 12, title, 3, 30, 30, 30);
	er_text(c, px + 14, 46, sub, 2, 95, 95, 95);
	ms_draw(c, obstacle, grey, NULL, 0, wire, 0.4);
	for (i = 0; i < nts; i++)
	{
		qaws_mesh_3d* g = ms_fc_at(m, mo, ts[i]);
		if (g)
			ms_draw(c, g, ghost, NULL, 0, wire, 0.2);
		qaws_mesh_3d_destroy(g);
	}
	if (r.contact)
	{
		qaws_mesh_3d* g = ms_fc_at(m, mo, (double)r.time);
		double pt[3] = { r.point.x, r.point.y, r.point.z };
		double bias = c->depth_bias;
		if (g)
			ms_draw(c, g, solid, NULL, 0, wire, 0.4);
		qaws_mesh_3d_destroy(g);
		c->depth_bias = 1e30;
		er_point(c, pt, red, 7);
		c->depth_bias = bias;
	}
	printf("  %s: %s\n", title, sub);
	free(fit);
}

static void demo_first_contact(void)
{
	er_canvas* c = er_create(1800, 640, 3, 252, 252, 250);
	qaws_motion_pose s, e;
	qaws_motion* m = NULL;
	if (!c)
		return;
	{
		qaws_mesh_3d *b = ms_at_ball(0.5, 0, 0, 0), *slab = ms_sw_box(-1.5, 1.5, -1, 1, -0.4, 0);
		double ts[6] = { 0, 0.1, 0.2, 0.3, 0.4, 0.45 };
		memset(&s, 0, sizeof s);
		memset(&e, 0, sizeof e);
		s.rotation[0] = e.rotation[0] = 1;
		s.position.z = 2;
		e.position.z = -1;
		s.velocity.z = e.velocity.z = -3;
		qaws_motion_create_poses(&m, &s, &e, 0);
		ms_fc_panel(c, 0, b, m, slab, ts, 6, "a ball falling on a slab", -0.6, 0.35);
		qaws_motion_destroy(m);
		m = NULL;
		qaws_mesh_3d_destroy(b); qaws_mesh_3d_destroy(slab);
	}
	{
		qaws_mesh_3d *cube = ms_sw_box(-0.1, 0.1, -0.1, 0.1, -0.1, 0.1), *wall = ms_sw_box(2, 2.02, -0.6, 0.6, -0.6, 0.6);
		double ts[3] = { 0, 0.15, 0.3 };
		memset(&s, 0, sizeof s);
		memset(&e, 0, sizeof e);
		s.rotation[0] = e.rotation[0] = 1;
		e.position.x = 10;
		s.velocity.x = e.velocity.x = 10;
		qaws_motion_create_poses(&m, &s, &e, 0);
		ms_fc_panel(c, 600, cube, m, wall, ts, 3, "frames miss a thin wall", 0.7, 0.4);
		qaws_motion_destroy(m);
		m = NULL;
		qaws_mesh_3d_destroy(cube); qaws_mesh_3d_destroy(wall);
	}
	{
		qaws_mesh_3d *rod = ms_sw_box(0, 1, -0.03, 0.03, -0.03, 0.03), *wall = ms_sw_box(-1.2, 1.2, 0.5, 0.6, -0.3, 0.3);
		double ts[4] = { 0, 0.1, 0.2, 0.3 };
		memset(&s, 0, sizeof s);
		memset(&e, 0, sizeof e);
		s.rotation[0] = 1;
		e.rotation[0] = (qaws_scalar)cos(MS_PI / 4);
		e.rotation[3] = (qaws_scalar)sin(MS_PI / 4);
		s.angular_velocity.z = e.angular_velocity.z = (qaws_scalar)(MS_PI / 2);
		qaws_motion_create_poses(&m, &s, &e, 0);
		ms_fc_panel(c, 1200, rod, m, wall, ts, 4, "a rod turning into a wall", 0, 1.3);
		qaws_motion_destroy(m);
		qaws_mesh_3d_destroy(rod); qaws_mesh_3d_destroy(wall);
	}
	er_save(c, "showcase/mesh71_first_contact.png");
	printf("mesh71_first_contact.png\n");
	er_destroy(c);
}
