/*
 * Figure 2: the surface between the vertices.
 *   - a coarse blob of 80 triangles, flat;
 *   - its PN reference, sampled finely: smooth through the same vertices,
 *     the patches' edge curves drawn;
 *   - a coarse can (8 sides) whose rims are not features: PN rounds it all;
 *   - the same can with its creases detected: rims sharp, sides round,
 *     caps flat, no cracks;
 *   - mean curvature read from the reference of a finer blob.
 */

/* the reference sampled n x n per triangle: smooth shaded by its normal
   field, coloured by a value per sample (curv: mean curvature mapped from
   [lo, hi]) or one colour; its triangles' edge curves drawn when wire */
static void ms_draw_reference(er_canvas* c, qaws_mesh_reference const* r, unsigned int n, double const* rgb, int curv, double lo, double hi,
	double const* wire, double wire_width)
{
	unsigned int nt = qaws_mesh_reference_get_triangle_count(r), t, i, j, k, m = (n + 1) * (n + 2) / 2;
	qaws_mesh_location* loc = (qaws_mesh_location*)malloc(sizeof(qaws_mesh_location) * m);
	qaws_scalar *x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * m), *nr = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * m),
		*kv = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 8 * m);
	double *X = (double*)malloc(sizeof(double) * 3 * m), *N = (double*)malloc(sizeof(double) * 3 * m), *C = (double*)malloc(sizeof(double) * 3 * m);
	for (t = 0; t < nt; t++)
	{
		unsigned int q = 0;
		for (j = 0; j <= n; j++)
			for (i = 0; i + j <= n; i++, q++)
			{
				loc[q].patch = t;
				loc[q].a = (qaws_scalar)((double)i / n);
				loc[q].b = (qaws_scalar)((double)j / n);
			}
		qaws_mesh_reference_evaluate(r, loc, m, x, nr, curv ? kv : NULL);
		for (q = 0; q < m; q++)
		{
			for (k = 0; k < 3; k++)
			{
				X[3 * q + k] = (double)x[3 * q + k];
				N[3 * q + k] = (double)nr[3 * q + k];
				C[3 * q + k] = rgb ? rgb[k] : 200;
			}
			if (curv)
				ms_colormap((0.5 * (kv[8 * q] + kv[8 * q + 1]) - lo) / (hi - lo), C + 3 * q);
		}
		/* the sample (i, j) of row j is at j (2n + 3 - j) / 2 + i */
		for (j = 0; j < n; j++)
			for (i = 0; i + j < n; i++)
			{
				unsigned int a = j * (2 * n + 3 - j) / 2 + i, b = a + 1, d = (j + 1) * (2 * n + 2 - j) / 2 + i;
				er_triangle(c, X + 3 * a, X + 3 * b, X + 3 * d, N + 3 * a, N + 3 * b, N + 3 * d, C + 3 * a, C + 3 * b, C + 3 * d);
				if (i + j + 1 < n)
				{
					unsigned int e = d + 1;
					er_triangle(c, X + 3 * b, X + 3 * e, X + 3 * d, N + 3 * b, N + 3 * e, N + 3 * d, C + 3 * b, C + 3 * e, C + 3 * d);
				}
			}
		if (wire)
		{
			/* the three edge curves */
			for (k = 0; k < 3; k++)
				for (i = 0; i < n; i++)
				{
					qaws_mesh_location e2[2];
					qaws_scalar y[6];
					double y0[3], y1[3];
					unsigned int s;
					for (s = 0; s < 2; s++)
					{
						double u = (double)(i + s) / n, bb[3] = { 0, 0, 0 };
						bb[k] = 1 - u;
						bb[(k + 1) % 3] = u;
						e2[s].patch = t;
						e2[s].a = (qaws_scalar)bb[1];
						e2[s].b = (qaws_scalar)bb[2];
					}
					qaws_mesh_reference_evaluate(r, e2, 2, y, NULL, NULL);
					for (s = 0; s < 3; s++)
					{
						y0[s] = (double)y[s];
						y1[s] = (double)y[3 + s];
					}
					er_line(c, y0, y1, wire, wire_width);
				}
		}
	}
	free(loc); free(x); free(nr); free(kv); free(X); free(N); free(C);
}

/* an icosahedron split once, its points pushed in and out smoothly */
static qaws_mesh_3d* ms_blob(unsigned int levels)
{
	static unsigned int const f[60] = { 0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11, 1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8,
		3, 9, 4, 3, 4, 2, 3, 2, 6, 3, 6, 8, 3, 8, 9, 4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7, 9, 8, 1 };
	double g = (1 + sqrt(5.0)) / 2, x[36];
	double const v[36] = { -1, g, 0, 1, g, 0, -1, -g, 0, 1, -g, 0, 0, -1, g, 0, 1, g, 0, -1, -g, 0, 1, -g, g, 0, -1, g, 0, 1, -g, 0, -1, -g, 0, 1 };
	qaws_mesh_3d* m;
	unsigned int i, n, edits = 0;
	double* p;
	for (i = 0; i < 36; i++)
		x[i] = v[i] / sqrt(1 + g * g);
	m = ms_make(x, 12, f, NULL, 20);
	for (i = 0; i < levels; i++)
		ms_subdivide(m, &edits);
	n = qaws_mesh_3d_get_vertex_count(m);
	p = ms_positions(m);
	for (i = 0; i < n; i++)
	{
		double* q = p + 3 * i, l = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]), s;
		q[0] /= l; q[1] /= l; q[2] /= l;
		s = 1 + 0.32 * sin(2.2 * q[0] + 0.4) * cos(1.7 * q[1]) + 0.2 * q[2] * q[2];
		q[0] *= 1.25 * s; q[1] *= s; q[2] *= 0.9 * s;
	}
	{
		qaws_scalar* ps = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n);
		for (i = 0; i < 3 * n; i++)
			ps[i] = (qaws_scalar)p[i];
		qaws_mesh_3d_set_positions(m, ps);
		free(ps);
	}
	free(p);
	return m;
}

/* a can: n sides, a ring at the bottom and the top, each cap a fan round a centre */
static qaws_mesh_3d* ms_can(unsigned int n, double r, double h)
{
	double* x = (double*)malloc(sizeof(double) * 3 * (2 * n + 2));
	unsigned int* f = (unsigned int*)malloc(sizeof(unsigned int) * 3 * 4 * n), i, k = 0;
	qaws_mesh_3d* m;
	for (i = 0; i < n; i++)
	{
		double a = 2 * MS_PI * i / n;
		x[3 * i] = r * cos(a); x[3 * i + 1] = r * sin(a); x[3 * i + 2] = -h / 2;
		x[3 * (n + i)] = r * cos(a); x[3 * (n + i) + 1] = r * sin(a); x[3 * (n + i) + 2] = h / 2;
	}
	x[6 * n] = 0; x[6 * n + 1] = 0; x[6 * n + 2] = -h / 2;
	x[6 * n + 3] = 0; x[6 * n + 4] = 0; x[6 * n + 5] = h / 2;
	for (i = 0; i < n; i++)
	{
		unsigned int a = i, b = (i + 1) % n;
		f[k++] = a; f[k++] = b; f[k++] = n + b;
		f[k++] = a; f[k++] = n + b; f[k++] = n + a;
		f[k++] = 2 * n; f[k++] = b; f[k++] = a;             /* bottom, facing down */
		f[k++] = 2 * n + 1; f[k++] = n + a; f[k++] = n + b; /* top */
	}
	m = ms_make(x, 2 * n + 2, f, NULL, 4 * n);
	free(x);
	free(f);
	return m;
}

static void demo_reference(void)
{
	static double const wire[3] = { 40, 45, 60 }, tint[3] = { 178, 192, 212 }, warm[3] = { 226, 200, 160 }, curve[3] = { 70, 60, 50 };
	er_canvas* c = er_create(2750, 640, 3, 252, 252, 250);
	qaws_mesh_3d *blob = ms_blob(1), *fine = ms_blob(3), *can = ms_can(8, 1, 1.6), *can2 = NULL;
	qaws_mesh_reference *rb = NULL, *rc = NULL, *rc2 = NULL, *rl = NULL, *rf = NULL;
	qaws_mesh_reference_desc d;
	char sub[200];
	double klo = 1e300, khi = -1e300;
	if (!c || !blob || !can)
		return;
	memset(&d, 0, sizeof d);
	d.interpolant = QAWS_INTERPOLANT_PN;
	qaws_mesh_reference_create(blob, &d, &rb);
	d.interpolant = QAWS_INTERPOLANT_LINEAR;
	qaws_mesh_reference_create(blob, &d, &rl);
	d.interpolant = QAWS_INTERPOLANT_PN;
	qaws_mesh_reference_create(fine, &d, &rf);
	qaws_mesh_reference_create(can, &d, &rc);
	qaws_mesh_3d_clone(can, &can2);
	qaws_mesh_3d_detect_features(can2, (qaws_scalar)(50 * MS_PI / 180), 0, 0);
	qaws_mesh_reference_create(can2, &d, &rc2);
	if (!rb || !rl || !rc || !rc2 || !rf)
		return;
	sprintf(sub, "%u triangles, flat", qaws_mesh_3d_get_face_count(blob));
	ms_panel(c, 0, 0, 550, 640, blob, 0.6, 0.5, "a coarse blob", sub);
	ms_draw(c, blob, tint, NULL, 0, wire, 1.6);
	ms_panel(c, 550, 0, 550, 640, blob, 0.6, 0.5, "its PN reference", "smooth through the same vertices");
	ms_draw_reference(c, rb, 12, tint, 0, 0, 0, curve, 1.4);
	sprintf(sub, "%u triangles, no features: all rounded", qaws_mesh_3d_get_face_count(can));
	ms_panel(c, 1100, 0, 550, 640, can, 0.35, 0.5, "a coarse can", sub);
	ms_draw_reference(c, rc, 12, warm, 0, 0, 0, curve, 1.4);
	sprintf(sub, "%u crease edges: rims sharp, caps flat", qaws_mesh_3d_get_feature_edge_count(can2));
	ms_panel(c, 1650, 0, 550, 640, can2, 0.35, 0.5, "with its creases", sub);
	ms_draw_reference(c, rc2, 12, warm, 0, 0, 0, curve, 1.4);
	{
		/* the blob's mean curvature range, then its colours */
		unsigned int nt = qaws_mesh_reference_get_triangle_count(rf), t, i;
		for (t = 0; t < nt; t++)
			for (i = 0; i < 6; i++)
			{
				qaws_mesh_location l;
				qaws_scalar k[8];
				l.patch = t;
				l.a = (qaws_scalar)(0.15 + 0.12 * (i % 3));
				l.b = (qaws_scalar)(0.15 + 0.12 * (i / 3));
				qaws_mesh_reference_evaluate(rf, &l, 1, NULL, NULL, k);
				klo = fmin(klo, 0.5 * (k[0] + k[1]));
				khi = fmax(khi, 0.5 * (k[0] + k[1]));
			}
		sprintf(sub, "%u triangles: %.2f to %.2f", qaws_mesh_3d_get_face_count(fine), klo, khi);
		ms_panel(c, 2200, 0, 550, 640, fine, 0.6, 0.5, "mean curvature", sub);
		ms_draw_reference(c, rf, 4, NULL, 1, klo, khi, NULL, 0);
	}
	er_save(c, "showcase/mesh2_reference.png");
	printf("mesh2_reference.png: blob %u triangles, can %u triangles, %u creases; mean curvature %.3f to %.3f\n", qaws_mesh_3d_get_face_count(blob),
		qaws_mesh_3d_get_face_count(can), qaws_mesh_3d_get_feature_edge_count(can2), klo, khi);
	qaws_mesh_reference_destroy(rb);
	qaws_mesh_reference_destroy(rl);
	qaws_mesh_reference_destroy(rf);
	qaws_mesh_3d_destroy(fine);
	qaws_mesh_reference_destroy(rc);
	qaws_mesh_reference_destroy(rc2);
	qaws_mesh_3d_destroy(blob);
	qaws_mesh_3d_destroy(can);
	qaws_mesh_3d_destroy(can2);
	er_destroy(c);
}
