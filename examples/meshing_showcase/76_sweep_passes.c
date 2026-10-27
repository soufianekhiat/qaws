/*
 * Figure 76: passes over one place (qaws_sweep_desc.combines, the
 * generated "passes" and "dwell"), swept implicitly.
 *   - a box beside z turning three turns while drifting along x: each place
 *     coloured by how many turns passed over it;
 *   - the same by the time it stayed within reach (dwell);
 *   - a ball out along x and back the way it came: "time" combined FIRST,
 *     the way out (0 to 0.5);
 *   - the same combined LAST: the way back (1 to 0.5), the far end once.
 * Printed: the ranges.
 */

/* per corner colours of a 1-component corner attribute on the ramp, its value from lo to hi */
static double* ms_ps_colours(qaws_mesh_3d const* m, char const* name, double lo, double hi, double* vlo, double* vhi)
{
	unsigned int nt = qaws_mesh_3d_get_triangle_count(m), i, el;
	double *v = ms_at_get(m, name, &el), *o;
	if (!v)
		return NULL;
	o = (double*)malloc(sizeof(double) * 9 * (nt + 1));
	*vlo = 1e30;
	*vhi = -1e30;
	for (i = 0; i < 3 * nt; i++)
	{
		*vlo = fmin(*vlo, v[i]);
		*vhi = fmax(*vhi, v[i]);
		ms_sg_ramp((v[i] - lo) / (hi - lo), o + 3 * i);
	}
	free(v);
	return o;
}

static void demo_sweep_passes(void)
{
	er_canvas* c = er_create(2400, 640, 3, 252, 252, 250);
	qaws_motion_pose s, e;
	qaws_motion* m = NULL;
	qaws_sweep_desc d;
	qaws_attribute_combine cb;
	char sub[200];
	double lo, hi, *cc;
	if (!c)
		return;
	{
		qaws_mesh_3d *b = ms_sw_box(0.4, 1, -0.15, 0.15, 0, 0.3), *o = NULL;
		memset(&s, 0, sizeof s);
		memset(&e, 0, sizeof e);
		s.rotation[0] = e.rotation[0] = 1;
		s.angular_velocity.z = e.angular_velocity.z = (qaws_scalar)(6 * MS_PI);
		e.position.x = (qaws_scalar)0.9;
		s.velocity.x = e.velocity.x = (qaws_scalar)0.9;
		qaws_motion_create_poses(&m, &s, &e, 0);
		memset(&d, 0, sizeof d);
		d.resolution = 64;
		d.generate = QAWS_SWEEP_GENERATE_PASSES | QAWS_SWEEP_GENERATE_DWELL;
		if (qaws_mesh_3d_sweep_implicit(&o, NULL, b, m, &d) == QAWS_STATUS_OK && o)
		{
			if ((cc = ms_ps_colours(o, "passes", 1, 3, &lo, &hi)) != NULL)
			{
				sprintf(sub, "\"passes\": %.0f to %.0f", lo, hi);
				ms_panel(c, 0, 0, 600, 640, o, -0.35, 1.0, "turns over each place", sub);
				ms_at_draw(c, o, cc, NULL, NULL, NULL, 0);
				printf("  passes %s\n", sub);
				free(cc);
			}
			if ((cc = ms_ps_colours(o, "dwell", 0, 0.12, &lo, &hi)) != NULL)
			{
				sprintf(sub, "\"dwell\": %.3f to %.3f of the motion", lo, hi);
				ms_panel(c, 600, 0, 600, 640, o, -0.35, 1.0, "time within reach", sub);
				ms_at_draw(c, o, cc, NULL, NULL, NULL, 0);
				printf("  dwell %s\n", sub);
				free(cc);
			}
		}
		qaws_motion_destroy(m);
		m = NULL;
		qaws_mesh_3d_destroy(b); qaws_mesh_3d_destroy(o);
	}
	{
		qaws_mesh_3d *ball = ms_at_ball(0.2, 0, 0, 0), *o = NULL;
		unsigned int pass;
		memset(&s, 0, sizeof s);
		memset(&e, 0, sizeof e);
		s.rotation[0] = e.rotation[0] = 1;
		s.velocity.x = 3;
		e.velocity.x = -3;   /* (x = 3 t (1 - t): back the way it came) */
		qaws_motion_create_poses(&m, &s, &e, 0);
		for (pass = 0; pass < 2; pass++)
		{
			memset(&d, 0, sizeof d);
			d.resolution = 64;
			d.generate = QAWS_SWEEP_GENERATE_TIME;
			cb.name = "time";
			cb.combine = pass ? QAWS_SWEEP_COMBINE_LAST : QAWS_SWEEP_COMBINE_FIRST;
			d.combines = &cb;
			d.combine_count = 1;
			o = NULL;
			if (qaws_mesh_3d_sweep_implicit(&o, NULL, ball, m, &d) == QAWS_STATUS_OK && o && (cc = ms_ps_colours(o, "time", 0, 1, &lo, &hi)) != NULL)
			{
				sprintf(sub, "\"time\" %s: %.2f to %.2f", pass ? "LAST" : "FIRST", lo, hi);
				ms_panel(c, 1200 + 600 * (int)pass, 0, 600, 640, o, -0.3, 1.1, pass ? "out and back: the last pass" : "out and back: the first pass", sub);
				ms_at_draw(c, o, cc, NULL, NULL, NULL, 0);
				printf("  %s\n", sub);
				free(cc);
			}
			qaws_mesh_3d_destroy(o);
		}
		qaws_motion_destroy(m);
		qaws_mesh_3d_destroy(ball);
	}
	er_save(c, "showcase/mesh76_sweep_passes.png");
	printf("mesh76_sweep_passes.png\n");
	er_destroy(c);
}
