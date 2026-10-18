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
 * Second row, singularities (irregular vertices: blue valence 3, orange 5):
 *   - eight +1 cones placed at a cube's corner directions (red balls), by
 *     INSTANT and GLOBAL: a cube-sphere, its valence-3 vertices on them;
 *   - a +1 (red) and -1 (purple) pair placed on a square patch;
 *   - an attractor (strength 1, reach 0.1) drawing a free cone to a point:
 *     the quads without it, then with it.
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


/* a square [-1, 1]^2 at z = 0, an (n + 1)^2 grid split along alternating diagonals (vertex j (n + 1) + i) */
static qaws_mesh_3d* ms_square(unsigned int n, double** out_x, unsigned int* out_nv)
{
	unsigned int nv = (n + 1) * (n + 1), nf = 2 * n * n, i, j, k = 0;
	double* x = (double*)malloc(sizeof(double) * 3 * nv);
	unsigned int* f = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nf);
	qaws_mesh_3d* m;
	for (j = 0; j <= n; j++)
		for (i = 0; i <= n; i++)
		{
			x[3 * (j * (n + 1) + i)] = -1 + 2.0 * i / n;
			x[3 * (j * (n + 1) + i) + 1] = -1 + 2.0 * j / n;
			x[3 * (j * (n + 1) + i) + 2] = 0;
		}
	for (j = 0; j < n; j++)
		for (i = 0; i < n; i++)
		{
			unsigned int a = j * (n + 1) + i, b = a + 1, c = a + n + 1, e = c + 1;
			if ((i + j) & 1)
			{
				f[k++] = a; f[k++] = b; f[k++] = e;
				f[k++] = a; f[k++] = e; f[k++] = c;
			}
			else
			{
				f[k++] = a; f[k++] = b; f[k++] = c;
				f[k++] = b; f[k++] = e; f[k++] = c;
			}
		}
	m = ms_make(x, nv, f, NULL, nf);
	free(f);
	*out_x = x;
	*out_nv = nv;
	return m;
}

/* the vertex of x (nv of them) nearest p */
static unsigned int ms_nearest(double const* x, unsigned int nv, double const* p)
{
	unsigned int i, best = 0;
	double bd = HUGE_VAL;
	for (i = 0; i < nv; i++)
	{
		double a = x[3 * i] - p[0], b = x[3 * i + 1] - p[1], d = x[3 * i + 2] - p[2];
		if (a * a + b * b + d * d < bd)
		{
			bd = a * a + b * b + d * d;
			best = i;
		}
	}
	return best;
}

/* the interior vertices of valence 3 (blue) and 5 or more (orange) as balls, lifted by lift along their position (0: none);
   returns the distance from to (may be NULL) to the nearest of them */
static double ms_mark_irregular(er_canvas* c, qaws_mesh_3d const* m, double lift, double const* to)
{
	static double const three[3] = { 40, 110, 210 }, five[3] = { 240, 150, 30 };
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), nf = qaws_mesh_3d_get_face_count(m), ne = qaws_mesh_3d_get_edge_count(m), i, f, c0;
	unsigned int *fv = (unsigned int*)malloc(sizeof(unsigned int) * (qaws_mesh_3d_get_corner_count(m) + 1)), *fs = (unsigned int*)malloc(sizeof(unsigned int) * (nf + 1));
	unsigned int *e = (unsigned int*)malloc(sizeof(unsigned int) * (2 * ne + 2)), *deg = (unsigned int*)calloc(nv + 1, sizeof(unsigned int));
	unsigned int* ed = (unsigned int*)calloc(nv + 1, sizeof(unsigned int));
	double* x = ms_positions(m), near = HUGE_VAL;
	qaws_mesh_3d_get_faces(m, fv, fs);
	qaws_mesh_3d_get_edges(m, e);
	for (f = 0, c0 = 0; f < nf; c0 += fs[f], f++)
		for (i = 0; i < fs[f]; i++)
			deg[fv[c0 + i]]++;
	for (i = 0; i < ne; i++)
	{
		ed[e[2 * i]]++;
		ed[e[2 * i + 1]]++;
	}
	for (i = 0; i < nv; i++)
		if (ed[i] == deg[i] && deg[i] != 4)
		{
			double p[3];
			unsigned int k;
			for (k = 0; k < 3; k++)
				p[k] = x[3 * i + k] * (1 + lift);
			if (to)
				near = fmin(near, sqrt((x[3 * i] - to[0]) * (x[3 * i] - to[0]) + (x[3 * i + 1] - to[1]) * (x[3 * i + 1] - to[1]) + (x[3 * i + 2] - to[2]) * (x[3 * i + 2] - to[2])));
			if (lift == 0)
				p[2] += 0.01;
			er_point(c, p, deg[i] < 4 ? three : five, 6);
		}
	free(fv); free(fs); free(e); free(deg); free(ed); free(x);
	return near;
}

/* the second row: placed cones (a cube-sphere by INSTANT and GLOBAL, a +1 / -1 pair on a square) and an attractor */
static void demo_cones(er_canvas* c, int y)
{
	static double const wire[3] = { 40, 45, 60 }, tint[3] = { 178, 192, 212 }, red[3] = { 214, 39, 40 }, sand[3] = { 226, 214, 186 };
	static double const purple[3] = { 150, 60, 170 };
	double *x0, *xs, ap[3];
	unsigned int nv0, nvs, i, k, cv[8], pv[2];
	int ci[8] = { 1, 1, 1, 1, 1, 1, 1, 1 }, pi[2] = { 1, -1 };
	qaws_mesh_3d *s = ms_uv_sphere(24, 48, &x0, &nv0), *sq = ms_square(24, &xs, &nvs), *q[5] = { NULL, NULL, NULL, NULL, NULL };
	qaws_scalar apt[3], str = 1;
	qaws_quadrangulate_desc d;
	qaws_mesh_report rep[5];
	char sub[200];
	for (i = 0; i < 8; i++)
	{
		double dd[3];
		dd[0] = (i & 1 ? 1 : -1) / sqrt(3.0);
		dd[1] = (i & 2 ? 1 : -1) / sqrt(3.0);
		dd[2] = (i & 4 ? 1 : -1) / sqrt(3.0);
		cv[i] = ms_nearest(x0, nv0, dd);
	}
	{
		double da[3] = { -0.5, 0.1, 0 }, db[3] = { 0.5, -0.1, 0 };
		pv[0] = ms_nearest(xs, nvs, da);
		pv[1] = ms_nearest(xs, nvs, db);
	}
	ap[0] = 0.62; ap[1] = -0.45; ap[2] = 0.64;
	{
		double l = sqrt(ap[0] * ap[0] + ap[1] * ap[1] + ap[2] * ap[2]);
		for (k = 0; k < 3; k++)
		{
			ap[k] /= l;
			apt[k] = (qaws_scalar)ap[k];
		}
	}
	for (i = 0; i < 5; i++)
	{
		clock_t t0 = clock();
		qaws_status st;
		memset(&d, 0, sizeof d);
		d.target_length = (qaws_scalar)(i == 2 ? 0.12 : 0.15);
		d.pure = 1;
		d.method = i == 1 ? QAWS_QUADS_GLOBAL : QAWS_QUADS_INSTANT;
		if (i < 2)
		{
			d.cone_vertices = cv;
			d.cone_indices = ci;
			d.cone_count = 8;
		}
		if (i == 2)
		{
			d.cone_vertices = pv;
			d.cone_indices = pi;
			d.cone_count = 2;
		}
		if (i == 4)
		{
			d.attractors = apt;
			d.attractor_strengths = &str;
			d.attractor_count = 1;
			d.attractor_radius = (qaws_scalar)0.1;
		}
		st = qaws_mesh_quadrangulate(i == 2 ? sq : s, &d, &q[i], &rep[i]);
		printf("  cones panel %u: %s, %u quads, %u others, %u singular, %.1f s\n", i, qaws_status_to_string(st), rep[i].quality.quads,
			rep[i].quality.triangles + rep[i].quality.polygons, rep[i].quality.irregular, (double)(clock() - t0) / CLOCKS_PER_SEC);
	}
	for (i = 0; i < 2; i++)
		if (q[i])
		{
			sprintf(sub, "+1 at a cube's corners: %u irregular", rep[i].quality.irregular);
			ms_panel(c, 550 * (int)i, y, 550, 640, q[i], 0.5, 0.45, i ? "placed cones, GLOBAL" : "placed cones, INSTANT", sub);
			ms_draw(c, q[i], tint, NULL, 1, wire, 1.0);
			for (k = 0; k < 8; k++)
				er_point(c, x0 + 3 * cv[k], red, 9);
			ms_mark_irregular(c, q[i], 0.01, NULL);
		}
	if (q[2])
	{
		sprintf(sub, "+1 (red), -1 (purple): %u irregular", rep[2].quality.irregular);
		ms_panel(c, 1100, y, 550, 640, q[2], 0.0, 1.25, "a cone pair on a patch", sub);
		ms_draw(c, q[2], sand, NULL, 0, wire, 1.0);
		er_point(c, xs + 3 * pv[0], red, 9);
		er_point(c, xs + 3 * pv[1], purple, 9);
		ms_mark_irregular(c, q[2], 0, NULL);
	}
	for (i = 3; i < 5; i++)
		if (q[i])
		{
			double lifted[3], near;
			ms_panel(c, 550 * (int)i, y, 550, 640, q[i], 0.5, 0.45, i == 3 ? "no attractor" : "an attractor (red, 1, 0.1)", NULL);
			ms_draw(c, q[i], tint, NULL, 1, wire, 1.0);
			for (k = 0; k < 3; k++)
				lifted[k] = 1.02 * ap[k];
			er_point(c, lifted, i == 4 ? red : wire, 7);
			near = ms_mark_irregular(c, q[i], 0.01, ap);
			sprintf(sub, "%u irregular, the nearest %.3f from it", rep[i].quality.irregular, near);
			er_text(c, 550 * (int)i + 14, y + 46, sub, 2, 95, 95, 95);
			printf("  %s\n", sub);
		}
	for (i = 0; i < 5; i++)
		qaws_mesh_3d_destroy(q[i]);
	qaws_mesh_3d_destroy(s);
	qaws_mesh_3d_destroy(sq);
	free(x0);
	free(xs);
}

static void demo_guidance(void)
{
	static double const wire[3] = { 40, 45, 60 }, tint[3] = { 178, 192, 212 }, warm[3] = { 232, 196, 150 }, cool[3] = { 160, 190, 222 };
	static double const red[3] = { 214, 39, 40 }, dark[3] = { 30, 30, 40 }, green[3] = { 170, 205, 160 };
	er_canvas* c = er_create(2750, 1280, 3, 252, 252, 250);
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
	demo_cones(c, 640);
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
