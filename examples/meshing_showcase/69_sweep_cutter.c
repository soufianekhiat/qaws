/*
 * Figure 69: a one-flute cutter. A blade on one side of its axis (not
 * symmetric about it) spinning 20 turns while fed 0.05 a turn into a block,
 * the tolerance 0.002: the drift of a turn passes the tolerance, so by slabs
 * the sweep's hulls pile up turn after turn (it did not finish in minutes);
 * with qaws_sweep_desc.revolve_spins the blade is revolved (a disc) and swept
 * along the feed, the result past the true sweep by at most a turn's drift.
 *   - the cutter over the block, coloured by distance from its axis;
 *   - the slot carved (the floor in rings: each place the blade's point that
 *     first reached it);
 *   - the slot from above, the blade's true positions over it every half
 *     turn (black): the slot holds them all, its end within a turn's drift.
 * Printed: turns, a turn's drift, faces, time.
 */

static void demo_sweep_cutter(void)
{
	static double const wire[3] = { 110, 110, 120 }, ink[3] = { 20, 20, 20 };
	er_canvas* c = er_create(1800, 640, 3, 252, 252, 250);
	double blade[6] = { 0, -0.03, 0, 0.3, 0.03, 0.6 }, col[24], wood[24];
	qaws_mesh_3d *tool = ms_sw_box(blade[0], blade[3], blade[1], blade[4], blade[2], blade[5]), *block = ms_sw_box(-0.6, 1.6, -0.6, 0.6, -0.6, 0.25), *o = NULL;
	qaws_motion_pose p0, p1;
	qaws_motion* m = NULL;
	qaws_motion const* mv[1];
	qaws_sweep_desc d;
	qaws_sweep_boolean_report r;
	unsigned int i, el;
	char sub[200];
	clock_t t0;
	if (!c || !tool || !block)
		return;
	for (i = 0; i < 8; i++)
	{
		qaws_vec3 v = qaws_mesh_3d_get_vertex(tool, i);
		double h = 6 * (double)v.x;   /* (the hue by distance from the axis) */
		col[3 * i] = 0.55 + 0.4 * cos(h);
		col[3 * i + 1] = 0.55 + 0.4 * cos(h - 2.0944);
		col[3 * i + 2] = 0.55 + 0.4 * cos(h + 2.0944);
		wood[3 * i] = 0.84; wood[3 * i + 1] = 0.76; wood[3 * i + 2] = 0.62;
	}
	qaws_mesh_3d_detect_features(tool, (qaws_scalar)0.5, 0, 0);
	qaws_mesh_3d_detect_features(block, (qaws_scalar)0.5, 0, 0);
	ms_at_add(tool, "color", QAWS_MESH_VERTEX, 3, QAWS_ATTRIBUTE_OKLAB, col, 8);
	ms_at_add(block, "color", QAWS_MESH_VERTEX, 3, QAWS_ATTRIBUTE_OKLAB, wood, 8);
	/* 20 turns about z, fed 1 along x */
	memset(&p0, 0, sizeof p0);
	memset(&p1, 0, sizeof p1);
	p0.rotation[0] = p1.rotation[0] = 1;
	p0.angular_velocity.z = p1.angular_velocity.z = (qaws_scalar)(40 * MS_PI);
	p1.position.x = 1;
	p0.velocity.x = p1.velocity.x = 1;
	qaws_motion_create_poses(&m, &p0, &p1, 0);
	mv[0] = m;
	{
		/* the cutter over the block */
		double* fit = (double*)malloc(sizeof(double) * 48);
		double* x = ms_positions(block);
		double* cc;
		memcpy(fit, x, sizeof(double) * 24);
		free(x);
		x = ms_positions(tool);
		memcpy(fit + 24, x, sizeof(double) * 24);
		free(x);
		er_panel(c, 0, 72, 600, 640 - 72, fit, 16, -0.35, 0.95, 0.06);
		er_text(c, 14, 12, "a one-flute cutter", 3, 30, 30, 30);
		er_text(c, 14, 46, "20 turns, fed 0.05 a turn (25 tolerances)", 2, 95, 95, 95);
		cc = ms_spin_colours(block);
		if (cc)
			ms_at_draw(c, block, cc, NULL, NULL, wire, 0.3);
		free(cc);
		cc = ms_spin_colours(tool);
		if (cc)
			ms_at_draw(c, tool, cc, NULL, NULL, wire, 0.3);
		free(cc); free(fit);
	}
	memset(&d, 0, sizeof d);
	d.tolerance = (qaws_scalar)0.002;
	d.revolve_spins = 1;
	d.simplify = 1;
	t0 = clock();
	if (m && qaws_mesh_sweep_boolean(&o, &r, tool, mv, block, 1, QAWS_BOOLEAN_DIFFERENCE, &d) == QAWS_STATUS_OK && o)
	{
		double* cc = ms_at_colours(o, &el);
		double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
		sprintf(sub, "revolve_spins: %u faces, closed %d, %.2f s", qaws_mesh_3d_get_face_count(o), r.closed, secs);
		ms_panel(c, 600, 0, 600, 640, o, -0.35, 0.95, "the slot", sub);
		if (cc)
			ms_at_draw(c, o, cc, NULL, NULL, wire, 0.2);
		printf("  cutter: %s\n", sub);
		{
			/* from above: the slot and the blade's true positions every 0.3 of a turn */
			double* x = ms_positions(o);
			unsigned int nv = qaws_mesh_3d_get_vertex_count(o), k, n = 0;
			double* fit = (double*)malloc(sizeof(double) * 3 * (nv + 1));
			for (k = 0; k < nv; k++)
				if (x[3 * k + 2] > 0.24 && fabs(x[3 * k + 1]) < 0.45 && x[3 * k] > -0.45 && x[3 * k] < 1.45)
				{
					fit[3 * n] = x[3 * k]; fit[3 * n + 1] = x[3 * k + 1]; fit[3 * n + 2] = x[3 * k + 2];
					n++;
				}
			er_panel(c, 1200, 72, 600, 640 - 72, fit, n, 0, 1.5, 0.06);
			er_text(c, 1214, 12, "the blade every 0.3 turn", 3, 30, 30, 30);
			er_text(c, 1214, 46, "each inside the slot (its end within 0.05)", 2, 95, 95, 95);
			if (cc)
				ms_at_draw(c, o, cc, NULL, NULL, wire, 0.2);
			c->depth_bias = 1e30;
			for (k = 0; k <= 67; k++)
			{
				qaws_scalar tx[24], ox[24];
				double a[3], b[3];
				unsigned int q;
				qaws_mesh_3d_get_positions(tx, tool);
				qaws_motion_transform_points(ox, NULL, m, (qaws_scalar)(k / 67.0), tx, 8);   /* (every 0.3 of a turn) */
				/* (the blade's top: corners 4 5 6 7) */
				for (q = 0; q < 4; q++)
				{
					unsigned int i0 = 4 + q, i1 = 4 + (q + 1) % 4;
					a[0] = ox[3 * i0]; a[1] = ox[3 * i0 + 1]; a[2] = 0.26;
					b[0] = ox[3 * i1]; b[1] = ox[3 * i1 + 1]; b[2] = 0.26;
					er_line(c, a, b, ink, 1.2);
				}
			}
			free(x); free(fit);
		}
		free(cc);
	}
	er_save(c, "showcase/mesh69_sweep_cutter.png");
	printf("mesh69_sweep_cutter.png\n");
	qaws_motion_destroy(m);
	qaws_mesh_3d_destroy(tool); qaws_mesh_3d_destroy(block); qaws_mesh_3d_destroy(o);
	er_destroy(c);
}
