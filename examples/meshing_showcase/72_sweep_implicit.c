/*
 * Figure 72: implicit sweeps (qaws_mesh_3d_sweep_implicit) against the
 * explicit ones (qaws_mesh_3d_sweep).
 *   - a box turned a quarter, swept by slabs (hulls and their union);
 *   - the same, implicit on a grid of 48: coloured by the distance to the
 *     explicit sweep (white 0, red half a cell);
 *   - the one-flute cutter of figure 69 (20 turns, fed 0.05 a turn), implicit:
 *     the true sweep within the grid, no revolution assumed;
 *   - the same coloured by the distance to its revolved sweep (revolve_spins,
 *     red at 0.03; a turn's drift is 0.05): the revolution passes the true sweep where
 *     the blade's next turn has not reached yet;
 *   - attributes: a cube with a colour per side and a gradient over its top
 *     (corner values, sharp at its edges) along a curve while turning a
 *     quarter, implicit: each corner the colour of the face that made it,
 *     the top the gradient's colours where the top face passed, the walls
 *     its edges sweep split between their two faces' colours; the same swept
 *     by slabs to compare (its top the same streaks: the side faces' hulls'
 *     top walls lie in the top's plane and take its values).
 * Printed: volumes, distances, samples, times.
 */

/* per vertex of a: white at 0 from b's surface to red at scale; *far the farthest */
static double* ms_si_heat(qaws_mesh_3d const* a, qaws_mesh_3d const* b, double scale, double* far)
{
	qaws_mesh_bvh* t = NULL;
	unsigned int n = qaws_mesh_3d_get_vertex_count(a), i;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (n + 1));
	qaws_mesh_hit* h = (qaws_mesh_hit*)malloc(sizeof(qaws_mesh_hit) * (n + 1));
	double* rgb = (double*)malloc(sizeof(double) * 3 * (n + 1));
	*far = 0;
	qaws_mesh_bvh_create_mesh_3d(&t, b, NULL);
	qaws_mesh_3d_get_positions(x, a);
	qaws_mesh_bvh_closest(h, t, x, n, 0, NULL);
	for (i = 0; i < n; i++)
	{
		double f = fmin(1, (double)h[i].distance / scale);
		*far = fmax(*far, (double)h[i].distance);
		rgb[3 * i] = 245 - 25 * f;
		rgb[3 * i + 1] = 245 - 205 * f;
		rgb[3 * i + 2] = 245 - 205 * f;
	}
	qaws_mesh_bvh_destroy(t);
	free(x); free(h);
	return rgb;
}

static double ms_si_volume(qaws_mesh_3d const* m)
{
	double* x = ms_positions(m);
	unsigned int nt = qaws_mesh_3d_get_triangle_count(m), i, *t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (nt + 1));
	double v = 0;
	qaws_mesh_3d_get_triangles(t, m);
	for (i = 0; i < nt; i++)
	{
		double const *a = x + 3 * t[3 * i], *b = x + 3 * t[3 * i + 1], *e = x + 3 * t[3 * i + 2];
		v += (a[0] * (b[1] * e[2] - b[2] * e[1]) + a[1] * (b[2] * e[0] - b[0] * e[2]) + a[2] * (b[0] * e[1] - b[1] * e[0])) / 6;
	}
	free(x); free(t);
	return v;
}

static void demo_sweep_implicit(void)
{
	static double const blue[3] = { 120, 160, 225 }, orange[3] = { 230, 170, 100 }, wire[3] = { 110, 110, 120 };
	er_canvas* c = er_create(3600, 640, 3, 252, 252, 250);
	qaws_motion_pose s, e;
	qaws_motion* m = NULL;
	qaws_sweep_desc d;
	qaws_sweep_report re, ri;
	char sub[200];
	clock_t t0;
	double secs;
	if (!c)
		return;
	{
		qaws_mesh_3d *b = ms_sw_box(1, 2, -0.1, 0.1, 0, 0.6), *ex = NULL, *im = NULL;
		memset(&s, 0, sizeof s);
		memset(&e, 0, sizeof e);
		s.rotation[0] = 1;
		e.rotation[0] = (qaws_scalar)cos(MS_PI / 4);
		e.rotation[3] = (qaws_scalar)sin(MS_PI / 4);
		s.angular_velocity.z = e.angular_velocity.z = (qaws_scalar)(MS_PI / 2);
		qaws_motion_create_poses(&m, &s, &e, 0);
		memset(&d, 0, sizeof d);
		d.tolerance = (qaws_scalar)0.002;
		d.simplify = 1;
		d.resolution = 48;
		t0 = clock();
		qaws_mesh_3d_sweep(&ex, &re, b, m, &d);
		secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
		t0 = clock();
		qaws_mesh_3d_sweep_implicit(&im, &ri, b, m, &d);
		if (ex && im)
		{
			double far, back, cell = 2.0 / 48, *heat;
			sprintf(sub, "by slabs: volume %.4f, %.2f s", ms_si_volume(ex), secs);
			ms_panel(c, 0, 0, 600, 640, ex, -0.4, 0.8, "a box turned a quarter", sub);
			ms_draw(c, ex, blue, NULL, 0, wire, 0.3);
			printf("  %s\n", sub);
			secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
			heat = ms_si_heat(im, ex, 0.5 * cell, &far);
			{
				double* h2 = ms_si_heat(ex, im, 0.5 * cell, &back);
				free(h2);
			}
			sprintf(sub, "grid 48: volume %.4f, apart %.3f / %.3f", ms_si_volume(im), far, back);
			ms_panel(c, 600, 0, 600, 640, im, -0.4, 0.8, "implicit, by distance", sub);
			ms_draw(c, im, NULL, heat, 1, NULL, 0);
			printf("  implicit %s (%u samples, closed %d, %.2f s)\n", sub, ri.slabs, ri.closed, secs);
			free(heat);
		}
		qaws_motion_destroy(m);
		m = NULL;
		qaws_mesh_3d_destroy(b); qaws_mesh_3d_destroy(ex); qaws_mesh_3d_destroy(im);
	}
	{
		qaws_mesh_3d *tool = ms_sw_box(0, 0.3, -0.03, 0.03, 0, 0.6), *rv = NULL, *im = NULL;
		memset(&s, 0, sizeof s);
		memset(&e, 0, sizeof e);
		s.rotation[0] = e.rotation[0] = 1;
		s.angular_velocity.z = e.angular_velocity.z = (qaws_scalar)(40 * MS_PI);
		e.position.x = 1;
		s.velocity.x = e.velocity.x = 1;
		qaws_motion_create_poses(&m, &s, &e, 0);
		memset(&d, 0, sizeof d);
		d.tolerance = (qaws_scalar)0.002;
		d.resolution = 64;
		t0 = clock();
		qaws_mesh_3d_sweep_implicit(&im, &ri, tool, m, &d);
		secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
		d.revolve_spins = 1;
		qaws_mesh_3d_sweep(&rv, &re, tool, m, &d);
		if (im && rv)
		{
			double far, back, *heat;
			sprintf(sub, "implicit, grid 64: %u times, %.1f s", ri.slabs, secs);
			ms_panel(c, 1200, 0, 600, 640, im, -0.35, 0.95, "a one-flute cutter", sub);
			ms_draw(c, im, orange, NULL, 1, NULL, 0);
			printf("  cutter %s, volume %.4f, closed %d\n", sub, ms_si_volume(im), ri.closed);
			heat = ms_si_heat(im, rv, 0.03, &far);
			{
				double* h2 = ms_si_heat(rv, im, 0.03, &back);
				free(h2);
			}
			sprintf(sub, "revolved: volume %.4f, past by %.3f", ms_si_volume(rv), back);
			ms_panel(c, 1800, 0, 600, 640, im, -0.35, 0.95, "against its revolution", sub);
			ms_draw(c, im, NULL, heat, 1, NULL, 0);
			printf("  %s (implicit from it %.3f, drift a turn %.3f)\n", sub, far, (double)re.spin_drift);
			free(heat);
		}
		qaws_motion_destroy(m);
		qaws_mesh_3d_destroy(tool); qaws_mesh_3d_destroy(rv); qaws_mesh_3d_destroy(im);
	}
	{
		/* a colour per face of a cube, along a curve while turning a quarter */
		static double const pal[6][3] = { { 0.85, 0.3, 0.3 }, { 0.95, 0.8, 0.3 }, { 0.3, 0.7, 0.4 }, { 0.3, 0.5, 0.9 }, { 0.7, 0.4, 0.85 }, { 0.95, 0.6, 0.25 } };
		qaws_mesh_3d *cube = ms_sw_box(-0.2, 0.2, -0.2, 0.2, -0.2, 0.2), *im = NULL;
		double col[108];
		unsigned int i, el;
		static double const top[8][3] = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0.95, 0.25, 0.25 }, { 0.25, 0.4, 0.95 }, { 0.3, 0.85, 0.4 }, { 0.98, 0.98, 0.95 } };
		unsigned int tv[36];
		qaws_mesh_3d_get_triangles(tv, cube);
		for (i = 0; i < 108; i++)
			col[i] = pal[i / 18][i % 3];
		/* (the top face: a colour per corner, blended across it) */
		for (i = 6; i < 12; i++)
			for (el = 0; el < 3; el++)
				col[3 * i + el] = top[tv[i]][el];
		ms_at_add(cube, "color", QAWS_MESH_CORNER, 3, QAWS_ATTRIBUTE_OKLAB, col, 36);
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
		d.resolution = 64;
		t0 = clock();
		qaws_mesh_3d_sweep_implicit(&im, &ri, cube, m, &d);
		secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
		if (im)
		{
			double* cc = ms_at_colours(im, &el);
			sprintf(sub, "implicit, grid 64, colours carried: %.1f s", secs);
			ms_panel(c, 2400, 0, 600, 640, im, -0.5, 0.75, "a cube's faces", sub);
			if (cc)
				ms_at_draw(c, im, cc, NULL, NULL, NULL, 0);
			{
				qaws_mesh_3d* ex = NULL;
				d.tolerance = (qaws_scalar)0.002;
				if (qaws_mesh_3d_sweep(&ex, &re, cube, m, &d) == QAWS_STATUS_OK && ex)
				{
					double* ce = ms_at_colours(ex, &el);
					ms_panel(c, 3000, 0, 600, 640, ex, -0.5, 0.75, "by slabs, to compare", "the explicit sweep's colours");
					if (ce)
						ms_at_draw(c, ex, ce, NULL, NULL, NULL, 0);
					free(ce);
				}
				qaws_mesh_3d_destroy(ex);
			}
			printf("  attributes: %s, closed %d\n", sub, ri.closed);
			free(cc);
		}
		qaws_motion_destroy(m);
		qaws_mesh_3d_destroy(cube); qaws_mesh_3d_destroy(im);
	}
	er_save(c, "showcase/mesh72_sweep_implicit.png");
	printf("mesh72_sweep_implicit.png\n");
	er_destroy(c);
}
