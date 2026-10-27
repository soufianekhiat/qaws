/*
 * Figure 73: attributes a sweep generates (qaws_sweep_desc.generate).
 *   - a cube along a curve while turning a quarter (by slabs), coloured by
 *     "time": when each place was swept;
 *   - the same by "angle", the turn made by then (a hue per quarter turn);
 *   - a box beside z turning three turns while drifting along x faster than
 *     the tolerance a turn (swept implicitly), coloured by "index": the turn
 *     that first reached each place;
 *   - the cube's sweep textured through "source_face" and "source_uv": a
 *     checker on each face of the cube, its own hue per face, read where the
 *     point that swept each place lies on the cube (the sweep remeshed finely,
 *     its attributes carried, so the checker is read densely).
 * Printed: faces, the ranges.
 */

/* a ramp from dark blue through green to yellow, x in [0, 1] (0 to 255) */
static void ms_sg_ramp(double x, double* rgb)
{
	x = fmin(1, fmax(0, x));
	rgb[0] = 255 * fmin(1, fmax(0, 1.6 * x - 0.4));
	rgb[1] = 255 * (0.15 + 0.75 * sin(3.14159265358979 * 0.5 * fmin(1, 1.25 * x)));
	rgb[2] = 255 * fmax(0, 0.75 - 0.9 * x);
}

static void ms_sg_hue(double h, double* rgb)
{
	rgb[0] = 255 * (0.55 + 0.4 * cos(h));
	rgb[1] = 255 * (0.55 + 0.4 * cos(h - 2.0944));
	rgb[2] = 255 * (0.55 + 0.4 * cos(h + 2.0944));
}

/* per corner colours of m from a generated attribute: mode 0 time, 1 angle, 2 index, 3 the checker */
static double* ms_sg_colours(qaws_mesh_3d const* m, int mode, double* lo, double* hi)
{
	unsigned int nt = qaws_mesh_3d_get_face_count(m), i, el;
	double *v = ms_at_get(m, mode == 0 ? "time" : mode == 1 ? "angle" : mode == 2 ? "index" : "source_uv", &el), *f = mode == 3 ? ms_at_get(m, "source_face", &el) : NULL;
	double* o;
	if (!v || (mode == 3 && !f))
	{
		free(v); free(f);
		return NULL;
	}
	o = (double*)malloc(sizeof(double) * 9 * (nt + 1));
	*lo = 1e30;
	*hi = -1e30;
	for (i = 0; i < 3 * nt; i++)
	{
		double* c = o + 3 * i;
		if (mode == 3)
		{
			/* (a 4 x 4 checker on the face's triangle, the face's hue) */
			int cu = (int)floor(4 * v[2 * i]), cv = (int)floor(4 * v[2 * i + 1]);
			ms_sg_hue(0.9 * floor(f[i] / 2 + 0.5), c);
			if ((cu + cv) & 1)
			{
				c[0] *= 0.55; c[1] *= 0.55; c[2] *= 0.55;
			}
			continue;
		}
		*lo = fmin(*lo, v[i]);
		*hi = fmax(*hi, v[i]);
		if (mode == 0)
			ms_sg_ramp(v[i], c);
		else if (mode == 1)
			ms_sg_hue(4 * v[i], c);
		else
			ms_sg_ramp(v[i] / 2, c);
	}
	free(v); free(f);
	return o;
}

static void demo_sweep_generate(void)
{
	er_canvas* c = er_create(2400, 640, 3, 252, 252, 250);
	qaws_motion_pose s, e;
	qaws_motion* m = NULL;
	qaws_sweep_desc d;
	qaws_sweep_report r;
	char sub[200];
	double lo, hi, *cc;
	if (!c)
		return;
	{
		qaws_mesh_3d *cube = ms_sw_box(-0.2, 0.2, -0.2, 0.2, -0.2, 0.2), *o = NULL;
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
		d.generate = QAWS_SWEEP_GENERATE_TIME | QAWS_SWEEP_GENERATE_ANGLE | QAWS_SWEEP_GENERATE_SOURCE;
		if (qaws_mesh_3d_sweep(&o, &r, cube, m, &d) == QAWS_STATUS_OK && o)
		{
			if ((cc = ms_sg_colours(o, 0, &lo, &hi)) != NULL)
			{
				sprintf(sub, "\"time\": %.2f to %.2f, %u faces", lo, hi, qaws_mesh_3d_get_face_count(o));
				ms_panel(c, 0, 0, 600, 640, o, -0.5, 0.75, "time", sub);
				ms_at_draw(c, o, cc, NULL, NULL, NULL, 0);
				printf("  time: %s\n", sub);
				free(cc);
			}
			if ((cc = ms_sg_colours(o, 1, &lo, &hi)) != NULL)
			{
				sprintf(sub, "\"angle\": %.3f to %.3f rad", lo, hi);
				ms_panel(c, 600, 0, 600, 640, o, -0.5, 0.75, "angle turned", sub);
				ms_at_draw(c, o, cc, NULL, NULL, NULL, 0);
				printf("  angle: %s\n", sub);
				free(cc);
			}
			{
				/* (remeshed finely, its attributes carried and seams kept: the checker read densely) */
				qaws_remesh_desc rd;
				qaws_mesh_3d* fine = NULL;
				memset(&rd, 0, sizeof rd);
				rd.target_length = (qaws_scalar)0.012;
				if (qaws_mesh_remesh(&fine, NULL, o, &rd) == QAWS_STATUS_OK && fine && (cc = ms_sg_colours(fine, 3, &lo, &hi)) != NULL)
				{
					sprintf(sub, "checker by source_uv, remeshed: %u faces", qaws_mesh_3d_get_face_count(fine));
					ms_panel(c, 1800, 0, 600, 640, fine, -0.5, 0.75, "the source point", sub);
					ms_at_draw(c, fine, cc, NULL, NULL, NULL, 0);
					free(cc);
				}
				qaws_mesh_3d_destroy(fine);
			}
		}
		qaws_motion_destroy(m);
		m = NULL;
		qaws_mesh_3d_destroy(cube); qaws_mesh_3d_destroy(o);
	}
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
		d.tolerance = (qaws_scalar)0.004;
		d.resolution = 64;
		d.generate = QAWS_SWEEP_GENERATE_INDEX;
		if (qaws_mesh_3d_sweep(&o, &r, b, m, &d) == QAWS_STATUS_OK && o && (cc = ms_sg_colours(o, 2, &lo, &hi)) != NULL)
		{
			sprintf(sub, "\"index\" %.0f to %.0f (%.0f turns, implicit %d)", lo, hi, (double)r.turns, r.implicit);
			ms_panel(c, 1200, 0, 600, 640, o, -0.35, 1.0, "the turn first there", sub);
			ms_at_draw(c, o, cc, NULL, NULL, NULL, 0);
			printf("  index: %s\n", sub);
			free(cc);
		}
		qaws_motion_destroy(m);
		qaws_mesh_3d_destroy(b); qaws_mesh_3d_destroy(o);
	}
	er_save(c, "showcase/mesh73_sweep_generate.png");
	printf("mesh73_sweep_generate.png\n");
	er_destroy(c);
}
