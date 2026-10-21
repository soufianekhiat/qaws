/*
 * Figure 21: anisotropic remeshing to a metric (qaws_remesh_desc.metric).
 *   - a solution with a steep layer, f = tanh(12 (x - 1/2)), on the fine
 *     grid it was computed on;
 *   - that grid remeshed to f's Hessian metric (for f's linear
 *     interpolation within 0.01): long thin triangles along the layer,
 *     large ones where f is flat;
 *   - a cylinder remeshed to a constant metric, 0.08 around and 0.4 along.
 */

/* an n x n grid of the unit square in z = 0, each cell two triangles */
static qaws_mesh_3d* ms_square_grid(unsigned int n)
{
	unsigned int nv = (n + 1) * (n + 1), i, j, t = 0;
	double* x = (double*)malloc(sizeof(double) * 3 * nv);
	unsigned int* f = (unsigned int*)malloc(sizeof(unsigned int) * 6 * n * n);
	qaws_mesh_3d* m;
	for (j = 0; j <= n; j++)
		for (i = 0; i <= n; i++)
		{
			x[3 * (j * (n + 1) + i)] = (double)i / n;
			x[3 * (j * (n + 1) + i) + 1] = (double)j / n;
			x[3 * (j * (n + 1) + i) + 2] = 0;
		}
	for (j = 0; j < n; j++)
		for (i = 0; i < n; i++)
		{
			unsigned int a = j * (n + 1) + i, b = a + 1, c = a + n + 2, e = a + n + 1;
			f[t++] = a; f[t++] = b; f[t++] = c;
			f[t++] = a; f[t++] = c; f[t++] = e;
		}
	m = ms_make(x, nv, f, NULL, 2 * n * n);
	free(x);
	free(f);
	return m;
}

/* an open cylinder of radius 1, height 2, na around and nz up; its radial normals into nrm */
static qaws_mesh_3d* ms_cylinder_grid(unsigned int na, unsigned int nz, qaws_scalar* nrm)
{
	unsigned int nv = na * (nz + 1), i, j, t = 0;
	double* x = (double*)malloc(sizeof(double) * 3 * nv);
	unsigned int* f = (unsigned int*)malloc(sizeof(unsigned int) * 6 * na * nz);
	qaws_mesh_3d* m;
	for (j = 0; j <= nz; j++)
		for (i = 0; i < na; i++)
		{
			double a = 2 * 3.14159265358979323846 * i / na;
			unsigned int v = j * na + i;
			x[3 * v] = cos(a);
			x[3 * v + 1] = sin(a);
			x[3 * v + 2] = 2.0 * j / nz;
			nrm[3 * v] = (qaws_scalar)cos(a);
			nrm[3 * v + 1] = (qaws_scalar)sin(a);
			nrm[3 * v + 2] = 0;
		}
	for (j = 0; j < nz; j++)
		for (i = 0; i < na; i++)
		{
			unsigned int a = j * na + i, b = j * na + (i + 1) % na, c = b + na, e = a + na;
			f[t++] = a; f[t++] = b; f[t++] = c;
			f[t++] = a; f[t++] = c; f[t++] = e;
		}
	m = ms_make(x, nv, f, NULL, 2 * na * nz);
	free(x);
	free(f);
	return m;
}

static double ms_layer(double x)
{
	return tanh(12 * (x - 0.5));
}

/* each vertex coloured by the layer function of its x */
static double* ms_layer_colours(qaws_mesh_3d const* m)
{
	unsigned int n = qaws_mesh_3d_get_vertex_count(m), i;
	double *x = ms_positions(m), *rgb = (double*)malloc(sizeof(double) * 3 * (n ? n : 1));
	for (i = 0; i < n; i++)
		ms_colormap(0.5 + 0.5 * ms_layer(x[3 * i]), rgb + 3 * i);
	free(x);
	return rgb;
}

static void demo_metric(void)
{
	static double const wire[3] = { 40, 45, 60 }, tint[3] = { 178, 192, 212 };
	er_canvas* c = er_create(2200, 700, 3, 252, 252, 250);
	qaws_scalar* nrm = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * 64 * 17);
	qaws_mesh_3d *sq = ms_square_grid(64), *cyl = nrm ? ms_cylinder_grid(64, 16, nrm) : NULL, *a = NULL, *b = NULL;
	qaws_mesh_reference *rs = NULL, *rc = NULL;
	qaws_metric_field *ms = NULL, *mc = NULL;
	qaws_mesh_reference_desc rd;
	qaws_hessian_metric_desc hd;
	qaws_remesh_desc d;
	char sub[200];
	unsigned int i;
	if (!c || !sq || !cyl)
		return;
	/* the layer's Hessian metric on the square */
	qaws_mesh_reference_create(&rs, sq, NULL);
	{
		unsigned int nv = qaws_mesh_3d_get_vertex_count(sq);
		double* x = ms_positions(sq);
		qaws_scalar* f = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv);
		for (i = 0; i < nv; i++)
			f[i] = (qaws_scalar)ms_layer(x[3 * i]);
		memset(&hd, 0, sizeof hd);
		hd.error = 0.01f;
		hd.min_length = 0.01f;
		hd.max_length = 0.2f;
		if (rs)
			qaws_metric_field_create_hessian(&ms, rs, f, &hd);
		free(f);
		free(x);
	}
	memset(&d, 0, sizeof d);
	d.reference = rs;
	d.metric = ms;
	if (ms)
		qaws_mesh_remesh(&a, NULL, sq, &d);
	/* the cylinder's constant metric: 0.08 around, 0.4 along */
	memset(&rd, 0, sizeof rd);
	rd.normals = nrm;
	qaws_mesh_reference_create(&rc, cyl, &rd);
	{
		unsigned int nv = qaws_mesh_3d_get_vertex_count(cyl);
		qaws_scalar* M = (qaws_scalar*)calloc(6 * (size_t)nv, sizeof(qaws_scalar));
		for (i = 0; M && i < nv; i++)
		{
			M[6 * i] = M[6 * i + 3] = (qaws_scalar)(1 / (0.08 * 0.08));
			M[6 * i + 5] = (qaws_scalar)(1 / (0.4 * 0.4));
		}
		if (rc && M)
			qaws_metric_field_create_vertices(&mc, rc, M);
		free(M);
	}
	d.reference = rc;
	d.metric = mc;
	if (mc)
		qaws_mesh_remesh(&b, NULL, cyl, &d);
	{
		double* rgb = ms_layer_colours(sq);
		sprintf(sub, "%u triangles", qaws_mesh_3d_get_triangle_count(sq));
		ms_panel(c, 0, 0, 733, 700, sq, 0, 1.5707963, "a layer, tanh 12 (x - 1/2)", sub);
		ms_draw(c, sq, NULL, rgb, 0, wire, 0.5);
		printf("  the layer: %s\n", sub);
		free(rgb);
	}
	if (a)
	{
		double* rgb = ms_layer_colours(a);
		sprintf(sub, "%u triangles, for an interpolation error of 0.01", qaws_mesh_3d_get_triangle_count(a));
		ms_panel(c, 733, 0, 733, 700, a, 0, 1.5707963, "remeshed to its Hessian", sub);
		ms_draw(c, a, NULL, rgb, 0, wire, 1.2);
		printf("  its Hessian: %s\n", sub);
		free(rgb);
	}
	if (b)
	{
		sprintf(sub, "%u triangles, 0.08 around and 0.4 along", qaws_mesh_3d_get_triangle_count(b));
		ms_panel(c, 1466, 0, 733, 700, b, 0.6, 0.45, "a cylinder, stretched", sub);
		ms_draw(c, b, tint, NULL, 1, wire, 1.2);
		printf("  the cylinder: %s\n", sub);
	}
	er_save(c, "showcase/mesh21_metric.png");
	printf("mesh21_metric.png\n");
	qaws_mesh_3d_destroy(a);
	qaws_mesh_3d_destroy(b);
	qaws_metric_field_destroy(ms);
	qaws_metric_field_destroy(mc);
	qaws_mesh_reference_destroy(rs);
	qaws_mesh_reference_destroy(rc);
	qaws_mesh_3d_destroy(sq);
	qaws_mesh_3d_destroy(cyl);
	free(nrm);
	er_destroy(c);
}
