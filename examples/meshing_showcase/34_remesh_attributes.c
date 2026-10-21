/*
 * Figure 34: remeshing from attributes (qaws_mesh_3d_cut_isolines,
 * qaws_metric_field_create_intersection, qaws_metric_field_create_graded).
 *   - a sheet with a hill, its attribute x y cut along five isolines and
 *     remeshed (the lines creases of its own reference): they stay (red);
 *   - a flat square remeshed to the Hessian metric of an attribute with a
 *     steep ridge across x; then to the intersection with a second
 *     attribute's (a ridge across y): both ridges get their stretched
 *     triangles; then the intersection graded by 1.3: the sizes change
 *     gently away from the ridges.
 */

/* an n x n grid over [-1, 1]^2, height from z(x, y) */
static qaws_mesh_3d* ms_ra_sheet(unsigned int n, double (*z)(double, double))
{
	double* x = (double*)malloc(sizeof(double) * 3 * (n + 1) * (n + 1));
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 6 * n * n), i, j, k = 0;
	qaws_mesh_3d* m;
	for (j = 0; j <= n; j++)
		for (i = 0; i <= n; i++)
		{
			double u = 2.0 * i / n - 1, v = 2.0 * j / n - 1;
			x[3 * (j * (n + 1) + i)] = u;
			x[3 * (j * (n + 1) + i) + 1] = v;
			x[3 * (j * (n + 1) + i) + 2] = z(u, v);
		}
	for (j = 0; j < n; j++)
		for (i = 0; i < n; i++)
		{
			unsigned int a = j * (n + 1) + i;
			t[k++] = a; t[k++] = a + 1; t[k++] = a + n + 2;
			t[k++] = a; t[k++] = a + n + 2; t[k++] = a + n + 1;
		}
	m = ms_make(x, (n + 1) * (n + 1), t, NULL, 2 * n * n);
	free(x);
	free(t);
	return m;
}

static double ms_ra_hill(double u, double v)
{
	return 0.5 * exp(-(u * u + v * v) / 0.3);
}

static double ms_ra_flat(double u, double v)
{
	(void)u;
	(void)v;
	return 0;
}

static void ms_ra_features(er_canvas* c, qaws_mesh_3d const* m, double const* col)
{
	unsigned int n = qaws_mesh_3d_get_feature_edge_count(m), i;
	unsigned int* e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (n + 1));
	qaws_mesh_3d_get_feature_edges(e, m);
	for (i = 0; i < n; i++)
	{
		qaws_vec3 a = qaws_mesh_3d_get_vertex(m, e[2 * i]), b = qaws_mesh_3d_get_vertex(m, e[2 * i + 1]);
		double A[3] = { a.x, a.y, a.z }, B[3] = { b.x, b.y, b.z };
		er_line(c, A, B, col, 3);
	}
	free(e);
}

static void demo_remesh_attributes(void)
{
	static double const col[3] = { 180, 190, 214 }, warm[3] = { 226, 196, 150 }, wire[3] = { 60, 66, 80 }, red[3] = { 210, 40, 40 };
	static char const* const title[3] = { "the metric of a ridge across x", "and of one across y, intersected", "intersected, graded by 1.3" };
	er_canvas* c = er_create(2400, 720, 3, 252, 252, 250);
	char sub[200];
	unsigned int k;
	if (!c)
		return;
	{
		/* isolines of x y on the hill, kept through remeshing */
		qaws_mesh_3d *sheet = ms_ra_sheet(40, ms_ra_hill), *cut = NULL, *out = NULL;
		qaws_mesh_reference* ref = NULL;
		unsigned int nv = qaws_mesh_3d_get_vertex_count(sheet), i, a = 0, ne = 0;
		qaws_scalar* f = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv);
		qaws_scalar const levels[5] = { (qaws_scalar)-0.5, (qaws_scalar)-0.2, 0, (qaws_scalar)0.2, (qaws_scalar)0.5 };
		qaws_attribute_desc ad;
		qaws_remesh_desc d;
		for (i = 0; i < nv; i++)
		{
			qaws_vec3 p = qaws_mesh_3d_get_vertex(sheet, i);
			f[i] = (qaws_scalar)(p.x * p.y);
		}
		memset(&ad, 0, sizeof ad);
		ad.name = "f";
		ad.element = QAWS_MESH_VERTEX;
		ad.components = 1;
		qaws_mesh_3d_add_attribute(&a, sheet, f, &ad);
		qaws_mesh_reference_create(&ref, sheet, NULL);
		qaws_mesh_3d_clone(&cut, sheet);
		qaws_mesh_3d_cut_isolines(&ne, cut, a, levels, 5);
		memset(&d, 0, sizeof d);
		d.target_length = (qaws_scalar)0.07;
		qaws_mesh_remesh(&out, NULL, cut, &d);
		if (out)
		{
			sprintf(sub, "isolines of x y, remeshed: %u faces", qaws_mesh_3d_get_face_count(out));
			ms_panel(c, 0, 0, 600, 720, out, -0.5, 0.75, "isolines kept as lines", sub);
			ms_draw(c, out, warm, NULL, 1, wire, 0.4);
			ms_ra_features(c, out, red);
			printf("  isolines: %s\n", sub);
		}
		qaws_mesh_3d_destroy(out);
		qaws_mesh_3d_destroy(cut);
		qaws_mesh_reference_destroy(ref);
		qaws_mesh_3d_destroy(sheet);
		free(f);
	}
	for (k = 0; k < 3; k++)
	{
		/* a flat square: the Hessian metrics of two ridges, alone, intersected, graded */
		qaws_mesh_3d *sq = ms_ra_sheet(60, ms_ra_flat), *out = NULL;
		qaws_mesh_reference* ref = NULL;
		qaws_metric_field *m1 = NULL, *m2 = NULL, *mi = NULL, *mg = NULL;
		qaws_hessian_metric_desc hd;
		qaws_remesh_desc d;
		unsigned int nv = qaws_mesh_3d_get_vertex_count(sq), i;
		qaws_scalar* f1 = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv);
		qaws_scalar* f2 = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv);
		qaws_mesh_reference_desc rd;
		for (i = 0; i < nv; i++)
		{
			qaws_vec3 p = qaws_mesh_3d_get_vertex(sq, i);
			f1[i] = (qaws_scalar)tanh(12 * (p.x - 0.3 * p.y - 0.2));
			f2[i] = (qaws_scalar)tanh(12 * (p.y + 0.3));
		}
		memset(&rd, 0, sizeof rd);
		rd.interpolant = QAWS_INTERPOLANT_LINEAR;
		qaws_mesh_reference_create(&ref, sq, &rd);
		memset(&hd, 0, sizeof hd);
		hd.error = (qaws_scalar)0.02;
		hd.min_length = (qaws_scalar)0.01;
		hd.max_length = (qaws_scalar)0.25;
		qaws_metric_field_create_hessian(&m1, ref, f1, &hd);
		qaws_metric_field_create_hessian(&m2, ref, f2, &hd);
		{
			qaws_metric_field const* p[2] = { m1, m2 };
			qaws_metric_field_create_intersection(&mi, p, 2);
		}
		qaws_metric_field_create_graded(&mg, mi, (qaws_scalar)1.3);
		memset(&d, 0, sizeof d);
		d.reference = ref;
		d.metric = k == 0 ? m1 : k == 1 ? mi : mg;
		qaws_mesh_remesh(&out, NULL, sq, &d);
		if (out)
		{
			sprintf(sub, "%u faces", qaws_mesh_3d_get_face_count(out));
			ms_panel(c, 600 + 600 * (int)k, 0, 600, 720, sq, 0.0, 1.5707, title[k], sub);
			ms_draw(c, out, col, NULL, 0, wire, 0.5);
			printf("  %s: %s\n", title[k], sub);
		}
		qaws_mesh_3d_destroy(out);
		qaws_metric_field_destroy(m1); qaws_metric_field_destroy(m2); qaws_metric_field_destroy(mi); qaws_metric_field_destroy(mg);
		qaws_mesh_reference_destroy(ref);
		qaws_mesh_3d_destroy(sq);
		free(f1);
		free(f2);
	}
	er_save(c, "showcase/mesh34_remesh_attributes.png");
	printf("mesh34_remesh_attributes.png\n");
	er_destroy(c);
}
