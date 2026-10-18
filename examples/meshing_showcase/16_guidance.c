/*
 * Figure 16: user guidance of INSTANT quads (as test 137).
 *   - density paint: a sphere painted fine (0.07) on its upper half and
 *     coarse (0.14) below, all quads, the faces tinted by side;
 *   - a stroke of directions round the equator at 30 degrees to it (the
 *     guide's directions drawn at its input vertices): the quads without
 *     it, then pulled by it;
 *   - an edge loop: a parallel of the input given as guide edges (red), the
 *     quads' edges along it;
 *   - pinned vertices (red balls) kept as quad vertices where they are.
 */

/* a latitude-longitude sphere: the poles, then nlat - 1 parallels of nlon vertices (vertex 1 + (j - 1) nlon + i) */
static qaws_mesh_3d* ms_uv_sphere(unsigned int nlat, unsigned int nlon, double** out_x, unsigned int* out_nv)
{
	unsigned int nv = 2 + (nlat - 1) * nlon, nf = 2 * nlon * (nlat - 1), i, j, k = 0;
	double* x = (double*)malloc(sizeof(double) * 3 * nv);
	unsigned int* f = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nf);
	qaws_mesh_3d* m;
	x[0] = 0; x[1] = 0; x[2] = 1;
	x[3 * (nv - 1)] = 0; x[3 * (nv - 1) + 1] = 0; x[3 * (nv - 1) + 2] = -1;
	for (j = 1; j < nlat; j++)
		for (i = 0; i < nlon; i++)
		{
			double th = MS_PI * j / nlat, ph = 2 * MS_PI * (i + 0.5 * (j & 1)) / nlon;
			unsigned int v = 1 + (j - 1) * nlon + i;
			x[3 * v] = sin(th) * cos(ph);
			x[3 * v + 1] = sin(th) * sin(ph);
			x[3 * v + 2] = cos(th);
		}
	for (i = 0; i < nlon; i++)
	{
		unsigned int a = 1 + i, b = 1 + (i + 1) % nlon, c = 1 + (nlat - 2) * nlon + i, e = 1 + (nlat - 2) * nlon + (i + 1) % nlon;
		f[k++] = 0; f[k++] = a; f[k++] = b;
		f[k++] = nv - 1; f[k++] = e; f[k++] = c;
	}
	for (j = 1; j + 1 < nlat; j++)
		for (i = 0; i < nlon; i++)
		{
			unsigned int a = 1 + (j - 1) * nlon + i, b = 1 + (j - 1) * nlon + (i + 1) % nlon, c = 1 + j * nlon + i, e = 1 + j * nlon + (i + 1) % nlon;
			if (j & 1)
			{
				f[k++] = a; f[k++] = c; f[k++] = e;
				f[k++] = a; f[k++] = e; f[k++] = b;
			}
			else
			{
				f[k++] = a; f[k++] = c; f[k++] = b;
				f[k++] = b; f[k++] = c; f[k++] = e;
			}
		}
	m = ms_make(x, nv, f, NULL, nf);
	free(f);
	*out_x = x;
	*out_nv = nv;
	return m;
}

/* the stroke's direction at a point of the unit sphere: 30 degrees from the parallel toward the north */
static void ms_stroke_dir(double const* p, double* d)
{
	double r = sqrt(p[0] * p[0] + p[1] * p[1]);
	unsigned int k;
	double e[3], nth[3];
	e[0] = -p[1] / r; e[1] = p[0] / r; e[2] = 0;
	nth[0] = -p[0] * p[2] / r; nth[1] = -p[1] * p[2] / r; nth[2] = r;
	for (k = 0; k < 3; k++)
		d[k] = cos(MS_PI / 6) * e[k] + sin(MS_PI / 6) * nth[k];
}

/* the mean angle (degrees, 0 to 45) of the quads' edges within the stroke to its crosses */
static double ms_stroke_angle(qaws_mesh_3d const* m)
{
	unsigned int ne = qaws_mesh_3d_get_edge_count(m), *e = (unsigned int*)malloc(sizeof(unsigned int) * (2 * ne + 2)), i, k, n = 0;
	double* x = ms_positions(m), sum = 0;
	qaws_mesh_3d_get_edges(m, e);
	for (i = 0; i < ne; i++)
	{
		double const *p = x + 3 * e[2 * i], *u = x + 3 * e[2 * i + 1];
		double mid[3], t[3], d[3], s[3], l, dn, a, g;
		for (k = 0; k < 3; k++)
		{
			mid[k] = 0.5 * (p[k] + u[k]);
			t[k] = u[k] - p[k];
		}
		l = sqrt(mid[0] * mid[0] + mid[1] * mid[1] + mid[2] * mid[2]);
		for (k = 0; k < 3; k++)
			mid[k] /= l;
		if (fabs(mid[2]) > 0.15)
			continue;
		dn = t[0] * mid[0] + t[1] * mid[1] + t[2] * mid[2];
		for (k = 0; k < 3; k++)
			t[k] -= dn * mid[k];
		ms_stroke_dir(mid, d);
		s[0] = mid[1] * d[2] - mid[2] * d[1];
		s[1] = mid[2] * d[0] - mid[0] * d[2];
		s[2] = mid[0] * d[1] - mid[1] * d[0];
		a = atan2(t[0] * s[0] + t[1] * s[1] + t[2] * s[2], t[0] * d[0] + t[1] * d[1] + t[2] * d[2]);
		g = fmod(fabs(a) * 180 / MS_PI, 90);
		sum += g > 45 ? 90 - g : g;
		n++;
	}
	free(e);
	free(x);
	return n ? sum / n : 0;
}

static void demo_guidance(void)
{
	static double const wire[3] = { 40, 45, 60 }, tint[3] = { 178, 192, 212 }, warm[3] = { 232, 196, 150 }, cool[3] = { 160, 190, 222 };
	static double const red[3] = { 214, 39, 40 }, dark[3] = { 30, 30, 40 }, green[3] = { 170, 205, 160 };
	er_canvas* c = er_create(2750, 640, 3, 252, 252, 250);
	double* x0;
	unsigned int nv0, i, k, nlon = 48, ring = 8, pins[6] = { 101, 377, 610, 845, 290, 700 };
	qaws_mesh_3d* s = ms_uv_sphere(24, nlon, &x0, &nv0);
	qaws_mesh_3d* q[5] = { NULL, NULL, NULL, NULL, NULL };
	qaws_mesh_reference* r = NULL;
	qaws_size_field* sf = NULL;
	qaws_scalar *sizes = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv0), *g = (qaws_scalar*)calloc(3 * nv0, sizeof(qaws_scalar));
	qaws_scalar* w = (qaws_scalar*)calloc(nv0, sizeof(qaws_scalar));
	unsigned int* ge = (unsigned int*)malloc(sizeof(unsigned int) * 2 * nlon);
	qaws_quadrangulate_desc d;
	qaws_mesh_report rep[5];
	char sub[200];
	if (!c || !s)
		return;
	for (i = 0; i < nv0; i++)
	{
		sizes[i] = (qaws_scalar)(x0[3 * i + 2] > 0 ? 0.07 : 0.14);
		if (fabs(x0[3 * i + 2]) < 0.2)
		{
			double dd[3];
			ms_stroke_dir(x0 + 3 * i, dd);
			for (k = 0; k < 3; k++)
				g[3 * i + k] = (qaws_scalar)dd[k];
			w[i] = 1;
		}
	}
	for (i = 0; i < nlon; i++)
	{
		ge[2 * i] = 1 + (ring - 1) * nlon + i;
		ge[2 * i + 1] = 1 + (ring - 1) * nlon + (i + 1) % nlon;
	}
	qaws_mesh_reference_create(s, NULL, &r);
	qaws_size_field_create_vertices(r, sizes, &sf);
	for (i = 0; i < 5; i++)
	{
		clock_t t0 = clock();
		qaws_status st;
		memset(&d, 0, sizeof d);
		d.target_length = (qaws_scalar)0.12;
		d.pure = 1;
		if (i == 0)
		{
			d.reference = r;
			d.size = sf;
		}
		if (i == 2)
		{
			d.guide_directions = g;
			d.guide_weights = w;
		}
		if (i == 3)
		{
			d.guide_edges = ge;
			d.guide_edge_count = nlon;
		}
		if (i == 4)
		{
			d.pinned = pins;
			d.pinned_count = 6;
		}
		st = qaws_mesh_quadrangulate(s, &d, &q[i], &rep[i]);
		printf("  panel %u: %s, %u quads, %u others, %u singular, %.1f s\n", i, qaws_status_to_string(st), rep[i].quality.quads,
			rep[i].quality.triangles + rep[i].quality.polygons, rep[i].quality.irregular, (double)(clock() - t0) / CLOCKS_PER_SEC);
	}
	/* 1. density */
	if (q[0])
	{
		unsigned int nf = qaws_mesh_3d_get_face_count(q[0]), ne = qaws_mesh_3d_get_edge_count(q[0]), f, c0;
		unsigned int *fv = (unsigned int*)malloc(sizeof(unsigned int) * (qaws_mesh_3d_get_corner_count(q[0]) + 1));
		unsigned int *fs = (unsigned int*)malloc(sizeof(unsigned int) * (nf + 1)), *e = (unsigned int*)malloc(sizeof(unsigned int) * (2 * ne + 2));
		double* p = ms_positions(q[0]), sum[2] = { 0, 0 };
		unsigned int cnt[2] = { 0, 0 };
		qaws_mesh_3d_get_faces(q[0], fv, fs);
		qaws_mesh_3d_get_edges(q[0], e);
		for (i = 0; i < ne; i++)
		{
			double const *a = p + 3 * e[2 * i], *b = p + 3 * e[2 * i + 1];
			double z = 0.5 * (a[2] + b[2]);
			if (fabs(z) < 0.3)
				continue;
			sum[z < 0] += sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
			cnt[z < 0]++;
		}
		sprintf(sub, "painted 0.07 / 0.14: edges %.3f / %.3f", cnt[0] ? sum[0] / cnt[0] : 0, cnt[1] ? sum[1] / cnt[1] : 0);
		ms_panel(c, 0, 0, 550, 640, q[0], 0.6, 0.15, "density paint", sub);
		for (f = 0, c0 = 0; f < nf; c0 += fs[f], f++)
		{
			double z = 0;
			for (k = 0; k < fs[f]; k++)
				z += p[3 * fv[c0 + k] + 2] / fs[f];
			for (k = 1; k + 1 < fs[f]; k++)
				er_triangle(c, p + 3 * fv[c0], p + 3 * fv[c0 + k], p + 3 * fv[c0 + k + 1], p + 3 * fv[c0], p + 3 * fv[c0 + k], p + 3 * fv[c0 + k + 1],
					z > 0 ? warm : cool, z > 0 ? warm : cool, z > 0 ? warm : cool);
		}
		for (i = 0; i < ne; i++)
			er_line(c, p + 3 * e[2 * i], p + 3 * e[2 * i + 1], wire, 1.0);
		free(fv); free(fs); free(e); free(p);
	}
	/* 2, 3. the stroke, without and with */
	for (i = 1; i < 3; i++)
		if (q[i])
		{
			sprintf(sub, "edges %.1f deg off the stroke, %u irregular", ms_stroke_angle(q[i]), rep[i].quality.irregular);
			ms_panel(c, 550 * (int)i, 0, 550, 640, q[i], 0.6, 0.15, i == 1 ? "without guidance" : "a stroke at 30 degrees", sub);
			ms_draw(c, q[i], tint, NULL, 1, wire, 1.0);
			for (k = 0; k < nv0; k++)
				if (w[k] != 0)
				{
					double a[3], b[3], h = 0.045, dd[3];
					unsigned int j;
					ms_stroke_dir(x0 + 3 * k, dd);
					for (j = 0; j < 3; j++)
					{
						a[j] = 1.01 * x0[3 * k + j] - h * dd[j];
						b[j] = 1.01 * x0[3 * k + j] + h * dd[j];
					}
					er_line(c, a, b, i == 2 ? red : dark, 2.2);
				}
		}
	/* 4. the edge loop */
	if (q[3])
	{
		double lift[6];
		sprintf(sub, "a parallel as guide edges, %u irregular", rep[3].quality.irregular);
		ms_panel(c, 1650, 0, 550, 640, q[3], 0.6, 0.55, "an edge loop", sub);
		ms_draw(c, q[3], green, NULL, 1, wire, 1.0);
		for (i = 0; i < nlon; i++)
		{
			for (k = 0; k < 3; k++)
			{
				lift[k] = 1.004 * x0[3 * ge[2 * i] + k];
				lift[3 + k] = 1.004 * x0[3 * ge[2 * i + 1] + k];
			}
			er_line(c, lift, lift + 3, red, 2.6);
		}
	}
	/* 5. pinned */
	if (q[4])
	{
		sprintf(sub, "6 input vertices kept, %u irregular", rep[4].quality.irregular);
		ms_panel(c, 2200, 0, 550, 640, q[4], 0.6, 0.15, "pinned vertices", sub);
		ms_draw(c, q[4], tint, NULL, 1, wire, 1.0);
		for (i = 0; i < 6; i++)
			er_point(c, x0 + 3 * pins[i], red, 7);
	}
	er_save(c, "showcase/mesh16_guidance.png");
	printf("mesh16_guidance.png\n");
	for (i = 0; i < 5; i++)
		qaws_mesh_3d_destroy(q[i]);
	qaws_size_field_destroy(sf);
	qaws_mesh_reference_destroy(r);
	qaws_mesh_3d_destroy(s);
	free(sizes); free(g); free(w); free(ge); free(x0);
	er_destroy(c);
}
