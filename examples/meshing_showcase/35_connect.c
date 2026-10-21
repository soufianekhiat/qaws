/*
 * Figure 35: the mesh connector (qaws_mesh_3d_connect).
 *   - two tubes at an angle, their ends joined by a bridge leaving each
 *     along its surface: in quads, then smoothed;
 *   - a wide tube of 32 sides and a narrow one of 12 offset from it: a
 *     funnel of triangles, its vertex counts graded from one loop to the
 *     other.
 * The bridges are tinted.
 */

/* an open tube of radius r along the axis through c in direction d (unit), length len, nu sides and nv rows */
static qaws_mesh_3d* ms_cn_tube(double const* c, double const* d, double r, double len, unsigned int nu, unsigned int nv)
{
	double* x = (double*)malloc(sizeof(double) * 3 * nu * (nv + 1));
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 6 * nu * nv), i, j, k = 0;
	double u[3], w[3], a0[3] = { 0, 0, 0 }, l;
	qaws_mesh_3d* m;
	a0[fabs(d[0]) < 0.9 ? 0 : 1] = 1;
	u[0] = d[1] * a0[2] - d[2] * a0[1];
	u[1] = d[2] * a0[0] - d[0] * a0[2];
	u[2] = d[0] * a0[1] - d[1] * a0[0];
	l = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
	u[0] /= l; u[1] /= l; u[2] /= l;
	w[0] = d[1] * u[2] - d[2] * u[1];
	w[1] = d[2] * u[0] - d[0] * u[2];
	w[2] = d[0] * u[1] - d[1] * u[0];
	for (j = 0; j <= nv; j++)
		for (i = 0; i < nu; i++)
		{
			double a = 2 * MS_PI * i / nu, s = len * j / nv;
			unsigned int q;
			for (q = 0; q < 3; q++)
				x[3 * (j * nu + i) + q] = c[q] + s * d[q] + r * (cos(a) * u[q] + sin(a) * w[q]);
		}
	for (j = 0; j < nv; j++)
		for (i = 0; i < nu; i++)
		{
			unsigned int a = j * nu + i, b = j * nu + (i + 1) % nu;
			t[k++] = a; t[k++] = b; t[k++] = b + nu;
			t[k++] = a; t[k++] = b + nu; t[k++] = a + nu;
		}
	m = ms_make(x, nu * (nv + 1), t, NULL, 2 * nu * nv);
	free(x);
	free(t);
	return m;
}

/* the result drawn: the bridge (part 1) tinted */
static void ms_cn_draw(er_canvas* c, qaws_mesh_3d const* m, double const* col, double const* tint, double const* wire)
{
	unsigned int nf = qaws_mesh_3d_get_face_count(m), nv = qaws_mesh_3d_get_vertex_count(m), at = 0, f, k, off = 0;
	qaws_scalar* part = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (nf + 1));
	unsigned int* fs = (unsigned int*)malloc(sizeof(unsigned int) * (nf + 1));
	unsigned int* fv = (unsigned int*)malloc(sizeof(unsigned int) * (qaws_mesh_3d_get_corner_count(m) + 1));
	double* rgb = (double*)malloc(sizeof(double) * 3 * nv);
	qaws_mesh_3d_find_attribute(&at, m, "part");
	qaws_mesh_3d_get_attribute(part, m, at);
	qaws_mesh_3d_get_faces(fv, fs, m);
	for (k = 0; k < nv; k++)
		memcpy(rgb + 3 * k, col, sizeof(double) * 3);
	for (f = 0; f < nf; off += fs[f], f++)
		if (part[f] == 1)
			for (k = 0; k < fs[f]; k++)
				memcpy(rgb + 3 * fv[off + k], tint, sizeof(double) * 3);
	ms_draw(c, m, NULL, rgb, 1, wire, 0.4);
	free(part); free(fs); free(fv); free(rgb);
}

static void demo_connect(void)
{
	static double const col[3] = { 180, 190, 214 }, tint[3] = { 236, 160, 90 }, wire[3] = { 60, 66, 80 };
	er_canvas* c = er_create(1800, 720, 3, 252, 252, 250);
	char sub[200];
	unsigned int k;
	if (!c)
		return;
	for (k = 0; k < 2; k++)
	{
		/* a tube up the z axis, another coming in from the side at 60 degrees */
		double const ca[3] = { 0, 0, -1.2 }, da[3] = { 0, 0, 1 }, cb[3] = { 1.6, 0, 1.4 };
		double db[3] = { 0.8660254037844386, 0, 0.5 };
		qaws_mesh_3d *a = ms_cn_tube(ca, da, 0.3, 1.2, 24, 8), *b = ms_cn_tube(cb, db, 0.3, 1.2, 24, 8), *m = NULL;
		qaws_connect_desc d;
		qaws_connect_report r;
		memset(&d, 0, sizeof d);
		d.quads = 1;
		d.flags = k ? QAWS_CONNECT_FAIR : 0;
		/* (a's far end ring at the top, b's first ring) */
		qaws_mesh_3d_connect(&m, &r, a, 24 * 8, b, 0, &d);
		if (m)
		{
			sprintf(sub, "%u rows, gap %.2f, min angle %.0f deg", r.rows, (double)r.gap, (double)r.min_angle * 180 / MS_PI);
			ms_panel(c, 600 * (int)k, 0, 600, 720, m, 0.3, 0.2, k ? "smoothed (bi-Laplacian flow)" : "two tubes joined", sub);
			ms_cn_draw(c, m, col, tint, wire);
			printf("  tubes%s: %s\n", k ? ", faired" : "", sub);
		}
		qaws_mesh_3d_destroy(m);
		qaws_mesh_3d_destroy(a);
		qaws_mesh_3d_destroy(b);
	}
	{
		double const ca[3] = { 0, 0, -1.4 }, da[3] = { 0, 0, 1 }, cb[3] = { 0.25, 0, 0.6 };
		qaws_mesh_3d *a = ms_cn_tube(ca, da, 0.55, 1.2, 32, 6), *b = ms_cn_tube(cb, da, 0.18, 1.0, 12, 6), *m = NULL;
		qaws_connect_desc d;
		qaws_connect_report r;
		memset(&d, 0, sizeof d);
		d.flags = QAWS_CONNECT_FAIR;
		qaws_mesh_3d_connect(&m, &r, a, 32 * 6, b, 0, &d);
		if (m)
		{
			sprintf(sub, "%u to %u sides, %u rows", r.loop_vertices[0], r.loop_vertices[1], r.rows);
			ms_panel(c, 1200, 0, 600, 720, m, 0.4, 0.3, "a funnel", sub);
			ms_cn_draw(c, m, col, tint, wire);
			printf("  funnel: %s\n", sub);
		}
		qaws_mesh_3d_destroy(m);
		qaws_mesh_3d_destroy(a);
		qaws_mesh_3d_destroy(b);
	}
	er_save(c, "showcase/mesh35_connect.png");
	printf("mesh35_connect.png\n");
	er_destroy(c);
}
