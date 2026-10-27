/*
 * Figure 68: spins. A wheel turning about its own axis passes over the same
 * places again and again; its sweep is its solid of revolution swept along
 * the drift (qaws_mesh_3d_sweep finds the spin from the motion).
 *   - a torus 100 mm across, its tube 10 mm, coloured by the angle round its
 *     ring, spinning at 10 turns a second and moving 1 mm a second, its rim
 *     3 mm into a 100 mm cube: 3 turns (0.3 s) as one motion;
 *   - the cube carved: the groove takes the colours of the pass that reaches
 *     each point first (passes coincide: the first stands);
 *   - the groove from above;
 *   - an octagonal wheel the same way: each point of its round groove takes
 *     the colour of the corner the turn brings there first.
 * Printed: the turns found, faces, the removed volume.
 */

/* a ring (torus about y through the origin, radius R, tube r) or a prism of nu sides (r its half width), coloured by the
   angle round the axis (stored sRGB, blended in Oklab) */
static qaws_mesh_3d* ms_spin_tool(unsigned int nu, unsigned int nv, double R, double r, int prism)
{
	unsigned int n = prism ? 2 * nu + 2 : nu * nv, nf = prism ? 4 * nu : 2 * nu * nv, i, j, k = 0;
	double* x = (double*)malloc(sizeof(double) * 3 * n);
	double* col = (double*)malloc(sizeof(double) * 3 * n);
	unsigned int* f = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nf);
	qaws_mesh_3d* m;
	double vol = 0;
	for (j = 0; j < (prism ? 2u : nv); j++)
		for (i = 0; i < nu; i++)
		{
			double u = 2 * MS_PI * i / nu, v = 2 * MS_PI * (j + 0.5 * (i & 1)) / nv, rr = prism ? R + r : R + r * cos(v), h = 3 * u;
			unsigned int q = j * nu + i;
			x[3 * q] = rr * cos(u);
			x[3 * q + 1] = prism ? (j ? r : -r) : r * sin(v);
			x[3 * q + 2] = rr * sin(u);
			col[3 * q] = 0.55 + 0.4 * cos(h);
			col[3 * q + 1] = 0.55 + 0.4 * cos(h - 2.0944);
			col[3 * q + 2] = 0.55 + 0.4 * cos(h + 2.0944);
		}
	if (prism)
		for (j = 0; j < 2; j++)
		{
			x[3 * (2 * nu + j)] = 0; x[3 * (2 * nu + j) + 1] = j ? r : -r; x[3 * (2 * nu + j) + 2] = 0;
			col[3 * (2 * nu + j)] = col[3 * (2 * nu + j) + 1] = col[3 * (2 * nu + j) + 2] = 0.6;
		}
	for (i = 0; i < nu; i++)
		if (prism)
		{
			unsigned int a = i, b = (i + 1) % nu, c = nu + (i + 1) % nu, e = nu + i;
			f[k++] = a; f[k++] = b; f[k++] = c;
			f[k++] = a; f[k++] = c; f[k++] = e;
			f[k++] = 2 * nu; f[k++] = b; f[k++] = a;
			f[k++] = 2 * nu + 1; f[k++] = e; f[k++] = c;
		}
		else
			for (j = 0; j < nv; j++)
			{
				unsigned int a = j * nu + i, b = j * nu + (i + 1) % nu, c = ((j + 1) % nv) * nu + (i + 1) % nu, e = ((j + 1) % nv) * nu + i;
				f[k++] = a; f[k++] = c; f[k++] = b;
				f[k++] = a; f[k++] = e; f[k++] = c;
			}
	for (i = 0; i < nf; i++)
	{
		double const *a = x + 3 * f[3 * i], *b = x + 3 * f[3 * i + 1], *c = x + 3 * f[3 * i + 2];
		vol += a[0] * (b[1] * c[2] - b[2] * c[1]) + a[1] * (b[2] * c[0] - b[0] * c[2]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
	}
	if (vol < 0)
		for (i = 0; i < nf; i++)
		{
			unsigned int z = f[3 * i + 1];
			f[3 * i + 1] = f[3 * i + 2];
			f[3 * i + 2] = z;
		}
	m = ms_make(x, n, f, NULL, nf);
	if (m)
	{
		qaws_mesh_3d_detect_features(m, (qaws_scalar)0.5, 0, 0);
		ms_at_add(m, "color", QAWS_MESH_VERTEX, 3, QAWS_ATTRIBUTE_OKLAB, col, n);
	}
	free(x); free(col); free(f);
	return m;
}

/* the colours per corner (0 to 255), a vertex attribute spread to its corners */
static double* ms_spin_colours(qaws_mesh_3d const* m)
{
	unsigned int el, nt = qaws_mesh_3d_get_triangle_count(m), i, k, *t;
	double *v = ms_at_get(m, "color", &el), *o;
	if (!v || el != QAWS_MESH_VERTEX)
	{
		free(v);
		return ms_at_colours(m, &el);
	}
	t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (nt + 1));
	o = (double*)malloc(sizeof(double) * 9 * (nt + 1));
	qaws_mesh_3d_get_triangles(t, m);
	for (i = 0; i < 3 * nt; i++)
		for (k = 0; k < 3; k++)
			o[3 * i + k] = 255 * fmin(1, fmax(0, v[3 * t[i] + k]));
	free(v); free(t);
	return o;
}

/* 10 turns a second about y, 1 mm a second along x, over T seconds, the wheel's centre at height cz */
static qaws_motion* ms_spin_motion(double T, double cz)
{
	qaws_motion_pose p0, p1;
	qaws_motion* m = NULL;
	double ab = 20 * MS_PI * T;
	memset(&p0, 0, sizeof p0);
	memset(&p1, 0, sizeof p1);
	p0.position.z = p1.position.z = (qaws_scalar)cz;
	p1.position.x = (qaws_scalar)T;
	p0.rotation[0] = 1;
	p1.rotation[0] = (qaws_scalar)cos(ab / 2);
	p1.rotation[2] = (qaws_scalar)sin(ab / 2);
	p0.velocity.x = p1.velocity.x = (qaws_scalar)T;
	p0.angular_velocity.y = p1.angular_velocity.y = (qaws_scalar)(20 * MS_PI * T);
	qaws_motion_create_poses(&m, &p0, &p1, 0);
	return m;
}

/* the cube carved by the spinning tool, drawn in a panel fitted to the given points (NULL: the result) */
static qaws_mesh_3d* ms_spin_carve(qaws_mesh_3d const* tool, qaws_mesh_3d const* cube, double cz, qaws_sweep_boolean_report* r)
{
	qaws_motion* m = ms_spin_motion(0.3, cz);
	qaws_motion const* mv[1];
	qaws_sweep_desc d;
	qaws_mesh_3d* o = NULL;
	mv[0] = m;
	memset(&d, 0, sizeof d);
	d.tolerance = (qaws_scalar)0.15;   /* (above the drift of a turn, 0.1 mm: a spin) */
	d.simplify = 1;
	if (m)
		qaws_mesh_sweep_boolean(&o, r, tool, mv, cube, 1, QAWS_BOOLEAN_DIFFERENCE, &d);
	qaws_motion_destroy(m);
	return o;
}

/* the result's vertices in and round the groove, for a close view */
static double* ms_spin_groove(qaws_mesh_3d const* o, unsigned int* n)
{
	unsigned int nv = qaws_mesh_3d_get_vertex_count(o), i, k = 0;
	double* x = ms_positions(o);
	double* g = (double*)malloc(sizeof(double) * 3 * (nv + 1));
	for (i = 0; i < nv; i++)
		if (x[3 * i + 2] < 50 + 1e-6 && x[3 * i + 2] > 40 && fabs(x[3 * i]) < 40 && fabs(x[3 * i + 1]) < 12)   /* (the groove and its rim on the top) */
		{
			g[3 * k] = x[3 * i]; g[3 * k + 1] = x[3 * i + 1]; g[3 * k + 2] = x[3 * i + 2];
			k++;
		}
	free(x);
	*n = k;
	return g;
}

static void demo_sweep_spin(void)
{
	static double const wire[3] = { 110, 110, 120 };
	er_canvas* c = er_create(2400, 640, 3, 252, 252, 250);
	qaws_mesh_3d *torus = ms_spin_tool(96, 16, 45, 5, 0), *wheel = ms_spin_tool(8, 2, 45, 5, 1), *cube = ms_sw_box(-50, 50, -50, 50, -50, 50), *o = NULL, *ow = NULL;
	double cz = 50 + 45 + 5 - 3, wood[24];
	unsigned int i, el, ng;
	char sub[200];
	qaws_sweep_boolean_report r, rw;
	if (!c || !torus || !wheel || !cube)
		return;
	for (i = 0; i < 8; i++)
	{
		wood[3 * i] = 0.84; wood[3 * i + 1] = 0.76; wood[3 * i + 2] = 0.62;
	}
	qaws_mesh_3d_detect_features(cube, (qaws_scalar)0.5, 0, 0);
	ms_at_add(cube, "color", QAWS_MESH_VERTEX, 3, QAWS_ATTRIBUTE_OKLAB, wood, 8);
	{
		/* the tool where it starts, over the cube */
		qaws_mesh_3d* placed = NULL;
		unsigned int nv = qaws_mesh_3d_get_vertex_count(torus);
		double* x = ms_positions(torus);
		qaws_scalar* xs = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
		double *cc, *fit;
		for (i = 0; i < nv; i++)
		{
			xs[3 * i] = (qaws_scalar)x[3 * i];
			xs[3 * i + 1] = (qaws_scalar)x[3 * i + 1];
			xs[3 * i + 2] = (qaws_scalar)(x[3 * i + 2] + cz);
		}
		qaws_mesh_3d_clone(&placed, torus);
		if (placed)
			qaws_mesh_3d_set_positions(placed, xs);
		fit = (double*)malloc(sizeof(double) * 3 * (nv + 8));
		for (i = 0; i < 3 * nv; i++) fit[i] = (double)xs[i];
		{
			double* cx = ms_positions(cube);
			memcpy(fit + 3 * nv, cx, sizeof(double) * 24);
			free(cx);
		}
		er_panel(c, 0, 72, 600, 640 - 72, fit, nv + 8, -0.35, 0.95, 0.04);
		er_text(c, 14, 12, "a spinning wheel", 3, 30, 30, 30);
		er_text(c, 14, 46, "torus 100/10 mm, 10 turns/s, 1 mm/s", 2, 95, 95, 95);
		cc = ms_spin_colours(cube);
		if (cc)
			ms_at_draw(c, cube, cc, NULL, NULL, wire, 0.3);
		free(cc);
		cc = placed ? ms_spin_colours(placed) : NULL;
		if (cc)
			ms_at_draw(c, placed, cc, NULL, NULL, wire, 0.15);
		free(cc); free(fit); free(x); free(xs);
		qaws_mesh_3d_destroy(placed);
	}
	o = ms_spin_carve(torus, cube, cz, &r);
	if (o)
	{
		double* cc = ms_at_colours(o, &el);
		double* g = ms_spin_groove(o, &ng);
		sprintf(sub, "removed %.1f mm3, %u faces, closed %d", (double)r.changed_volume, qaws_mesh_3d_get_face_count(o), r.closed);
		ms_panel(c, 600, 0, 600, 640, o, -0.35, 0.95, "3 turns as one motion", sub);
		if (cc)
			ms_at_draw(c, o, cc, NULL, NULL, wire, 0.15);
		printf("  torus: %s\n", sub);
		er_panel(c, 1200, 72, 600, 640 - 72, g, ng, 0, 1.5, 0.08);
		er_text(c, 1214, 12, "the groove from above", 3, 30, 30, 30);
		er_text(c, 1214, 46, "each point: the first pass there", 2, 95, 95, 95);
		if (cc)
			ms_at_draw(c, o, cc, NULL, NULL, wire, 0.15);
		free(cc); free(g);
	}
	ow = ms_spin_carve(wheel, cube, cz, &rw);
	if (ow)
	{
		double* cc = ms_at_colours(ow, &el);
		double* g = ms_spin_groove(ow, &ng);
		sprintf(sub, "removed %.1f mm3, %u faces, closed %d", (double)rw.changed_volume, qaws_mesh_3d_get_face_count(ow), rw.closed);
		er_panel(c, 1800, 72, 600, 640 - 72, g, ng, 0, 1.5, 0.08);
		er_text(c, 1814, 12, "an octagon wheel, from above", 3, 30, 30, 30);
		er_text(c, 1814, 46, sub, 2, 95, 95, 95);
		if (cc)
			ms_at_draw(c, ow, cc, NULL, NULL, wire, 0.15);
		printf("  octagon: %s\n", sub);
		free(cc); free(g);
	}
	er_save(c, "showcase/mesh68_sweep_spin.png");
	printf("mesh68_sweep_spin.png\n");
	qaws_mesh_3d_destroy(torus); qaws_mesh_3d_destroy(wheel); qaws_mesh_3d_destroy(cube); qaws_mesh_3d_destroy(o); qaws_mesh_3d_destroy(ow);
	er_destroy(c);
}
