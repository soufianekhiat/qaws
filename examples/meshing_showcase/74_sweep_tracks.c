/*
 * Figure 74: attributes animated over a sweep (qaws_sweep_desc.tracks).
 *   - a cube along a curve while turning, its colour a track: white hot at
 *     the start, orange, red, then dark grey (Oklab, blended between keys):
 *     the sweep shows it cooling along the path (by slabs);
 *   - the same swept implicitly: each place the colour at its first touch;
 *   - a box beside z turning three turns while drifting (implicitly), its
 *     colour held over each turn (blue, yellow, red): each place shows the
 *     pass that first reached it.
 * Printed: faces.
 */

/* a vertex colour track of n keys over the mesh's vertices, each key one colour */
static void ms_tk_keys(qaws_scalar* values, unsigned int nv, double const (*col)[3], unsigned int n)
{
	unsigned int k, i, j;
	for (k = 0; k < n; k++)
		for (i = 0; i < nv; i++)
			for (j = 0; j < 3; j++)
				values[((size_t)k * nv + i) * 3 + j] = (qaws_scalar)col[k][j];
}

static void demo_sweep_tracks(void)
{
	static double const glow[4][3] = { { 1.0, 0.97, 0.85 }, { 1.0, 0.6, 0.15 }, { 0.75, 0.12, 0.08 }, { 0.22, 0.22, 0.24 } };
	er_canvas* c = er_create(1800, 640, 3, 252, 252, 250);
	qaws_motion_pose s, e;
	qaws_motion* m = NULL;
	qaws_sweep_desc d;
	qaws_attribute_track tr;
	char sub[200];
	if (!c)
		return;
	{
		qaws_mesh_3d *cube = ms_sw_box(-0.2, 0.2, -0.2, 0.2, -0.2, 0.2), *o = NULL, *oi = NULL;
		double col[24];
		qaws_scalar times[4] = { 0, (qaws_scalar)0.3, (qaws_scalar)0.65, 1 }, keys[4 * 24];
		unsigned int i, el;
		for (i = 0; i < 24; i++)
			col[i] = glow[0][i % 3];
		ms_at_add(cube, "color", QAWS_MESH_VERTEX, 3, QAWS_ATTRIBUTE_OKLAB, col, 8);
		ms_tk_keys(keys, 8, glow, 4);
		memset(&s, 0, sizeof s);
		memset(&e, 0, sizeof e);
		s.rotation[0] = 1;
		e.rotation[0] = (qaws_scalar)cos(MS_PI / 4);
		e.rotation[3] = (qaws_scalar)sin(MS_PI / 4);
		s.angular_velocity.z = e.angular_velocity.z = (qaws_scalar)(MS_PI / 2);
		e.position.x = 2;
		s.velocity.x = 1; s.velocity.y = 2;
		e.velocity.x = 1; e.velocity.y = -2;
		qaws_motion_create_poses(&m, &s, &e, 0);
		memset(&d, 0, sizeof d);
		d.tolerance = (qaws_scalar)0.004;
		d.resolution = 64;
		tr.name = "color"; tr.times = times; tr.values = keys; tr.key_count = 4;
		d.tracks = &tr;
		d.track_count = 1;
		if (qaws_mesh_3d_sweep(&o, NULL, cube, m, &d) == QAWS_STATUS_OK && o)
		{
			double* cc = ms_at_colours(o, &el);
			sprintf(sub, "by slabs: %u faces", qaws_mesh_3d_get_face_count(o));
			ms_panel(c, 0, 0, 600, 640, o, -0.5, 0.75, "a tool cooling", sub);
			if (cc)
				ms_at_draw(c, o, cc, NULL, NULL, NULL, 0);
			printf("  cooling %s\n", sub);
			free(cc);
		}
		if (qaws_mesh_3d_sweep_implicit(&oi, NULL, cube, m, &d) == QAWS_STATUS_OK && oi)
		{
			double* cc = ms_at_colours(oi, &el);
			sprintf(sub, "implicit, grid 64: %u faces", qaws_mesh_3d_get_face_count(oi));
			ms_panel(c, 600, 0, 600, 640, oi, -0.5, 0.75, "the same, implicitly", sub);
			if (cc)
				ms_at_draw(c, oi, cc, NULL, NULL, NULL, 0);
			printf("  cooling %s\n", sub);
			free(cc);
		}
		qaws_motion_destroy(m);
		m = NULL;
		qaws_mesh_3d_destroy(cube); qaws_mesh_3d_destroy(o); qaws_mesh_3d_destroy(oi);
	}
	{
		static double const three[3][3] = { { 0.2, 0.45, 0.95 }, { 0.98, 0.78, 0.2 }, { 0.85, 0.2, 0.25 } };
		qaws_mesh_3d *b = ms_sw_box(0.4, 1, -0.15, 0.15, 0, 0.3), *o = NULL;
		double col[24], held[6][3];
		qaws_scalar times[6], keys[6 * 24];
		unsigned int i, el;
		for (i = 0; i < 24; i++)
			col[i] = three[0][i % 3];
		ms_at_add(b, "color", QAWS_MESH_VERTEX, 3, QAWS_ATTRIBUTE_OKLAB, col, 8);
		/* (a colour held over each turn: two keys a hair apart where a turn ends) */
		for (i = 0; i < 6; i++)
		{
			times[i] = (qaws_scalar)(i % 2 ? (i / 2 + 1) / 3.0 - (i == 5 ? 0 : 1e-6) : i / 2 / 3.0);
			memcpy(held[i], three[i / 2], sizeof held[i]);
		}
		ms_tk_keys(keys, 8, (double const(*)[3])held, 6);
		memset(&s, 0, sizeof s);
		memset(&e, 0, sizeof e);
		s.rotation[0] = e.rotation[0] = 1;
		s.angular_velocity.z = e.angular_velocity.z = (qaws_scalar)(6 * MS_PI);
		e.position.x = (qaws_scalar)0.9;
		s.velocity.x = e.velocity.x = (qaws_scalar)0.9;
		qaws_motion_create_poses(&m, &s, &e, 0);
		memset(&d, 0, sizeof d);
		d.tolerance = (qaws_scalar)0.004;
		d.resolution = 64;
		tr.name = "color"; tr.times = times; tr.values = keys; tr.key_count = 6;
		d.tracks = &tr;
		d.track_count = 1;
		if (qaws_mesh_3d_sweep(&o, NULL, b, m, &d) == QAWS_STATUS_OK && o)
		{
			double* cc = ms_at_colours(o, &el);
			sprintf(sub, "blue, yellow, red: turns 1 to 3, %u faces", qaws_mesh_3d_get_face_count(o));
			ms_panel(c, 1200, 0, 600, 640, o, -0.35, 1.0, "colour by pass", sub);
			if (cc)
				ms_at_draw(c, o, cc, NULL, NULL, NULL, 0);
			printf("  passes %s\n", sub);
			free(cc);
		}
		qaws_motion_destroy(m);
		qaws_mesh_3d_destroy(b); qaws_mesh_3d_destroy(o);
	}
	er_save(c, "showcase/mesh74_sweep_tracks.png");
	printf("mesh74_sweep_tracks.png\n");
	er_destroy(c);
}
