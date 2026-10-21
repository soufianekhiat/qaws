/*
 * Figure 25: watertight surfaces (qaws_mesh_make_watertight).
 *   - a sphere's triangles with three holes punched in it;
 *   - made watertight: the winding number's inside contoured on a grid,
 *     the holes closed by smooth caps;
 *   - two cubes crossing (a self-intersecting soup): their union, its
 *     edges and corners sharp (the crossings' planes from the input's own
 *     triangles);
 *   - an open wavy sheet wrapped at an offset: a closed shell round it.
 */

static qaws_mesh_3d* ms_wt_make(qaws_scalar const* x, unsigned int nv, unsigned int const* t, unsigned int nt, unsigned int res, double off,
	qaws_watertight_report* r)
{
	qaws_watertight_desc d;
	qaws_mesh_3d* m = NULL;
	memset(&d, 0, sizeof d);
	d.resolution = res;
	d.offset = (qaws_scalar)off;
	qaws_mesh_make_watertight(&m, r, x, t, NULL, nv, nt, &d);
	return m;
}

static void demo_watertight(void)
{
	static double const tint[3] = { 178, 192, 212 }, warm[3] = { 226, 196, 150 }, wire[3] = { 50, 55, 70 };
	static double const holes[3][3] = { { 0.8, 0.6, 0 }, { 0, 0.7071, 0.7071 }, { -0.6, 0.6, 0.53 } };
	er_canvas* c = er_create(2200, 640, 3, 252, 252, 250);
	unsigned int na = 40, nb = 20, nv = 2 + na * (nb - 1), i, j, t = 0, kept = 0;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	unsigned int* f = (unsigned int*)malloc(sizeof(unsigned int) * 6 * na * nb);
	unsigned int* k2 = (unsigned int*)malloc(sizeof(unsigned int) * 6 * na * nb);
	qaws_watertight_report r;
	qaws_mesh_3d *in = NULL, *out = NULL;
	char sub[200];
	if (!c || !x || !f || !k2)
		return;
	/* the sphere, its holes */
	x[0] = 0; x[1] = 0; x[2] = 1;
	x[3] = 0; x[4] = 0; x[5] = -1;
	for (j = 1; j < nb; j++)
		for (i = 0; i < na; i++)
		{
			double th = MS_PI * j / nb, ph = 2 * MS_PI * i / na;
			unsigned int v = 2 + (j - 1) * na + i;
			x[3 * v] = (qaws_scalar)(sin(th) * cos(ph));
			x[3 * v + 1] = (qaws_scalar)(sin(th) * sin(ph));
			x[3 * v + 2] = (qaws_scalar)cos(th);
		}
	for (i = 0; i < na; i++)
	{
		f[t++] = 0; f[t++] = 2 + i; f[t++] = 2 + (i + 1) % na;
		f[t++] = 1; f[t++] = 2 + (nb - 2) * na + (i + 1) % na; f[t++] = 2 + (nb - 2) * na + i;
	}
	for (j = 1; j + 1 < nb; j++)
		for (i = 0; i < na; i++)
		{
			unsigned int a = 2 + (j - 1) * na + i, b = 2 + (j - 1) * na + (i + 1) % na, cc = b + na, d = a + na;
			f[t++] = a; f[t++] = d; f[t++] = cc;
			f[t++] = a; f[t++] = cc; f[t++] = b;
		}
	for (i = 0; i < t / 3; i++)
	{
		double g[3] = { 0, 0, 0 };
		unsigned int k, h, drop = 0;
		for (k = 0; k < 3; k++)
			for (j = 0; j < 3; j++)
				g[j] += x[3 * f[3 * i + k] + j] / 3;
		for (h = 0; h < 3; h++)
			drop |= (g[0] - holes[h][0]) * (g[0] - holes[h][0]) + (g[1] - holes[h][1]) * (g[1] - holes[h][1]) + (g[2] - holes[h][2]) * (g[2] - holes[h][2])
				< 0.35 * 0.35;
		if (!drop)
		{
			memcpy(k2 + 3 * kept, f + 3 * i, sizeof(unsigned int) * 3);
			kept++;
		}
	}
	{
		double* xd = (double*)malloc(sizeof(double) * 3 * nv);
		for (i = 0; i < 3 * nv; i++)
			xd[i] = x[i];
		in = ms_make(xd, nv, k2, NULL, kept);
		free(xd);
	}
	out = ms_wt_make(x, nv, k2, kept, 48, 0, &r);
	if (in)
	{
		sprintf(sub, "%u triangles, three holes", kept);
		ms_panel(c, 0, 0, 550, 640, in, 3.3, 0.35, "a sphere with holes", sub);
		ms_draw(c, in, tint, NULL, 0, wire, 0.6);
	}
	if (out)
	{
		sprintf(sub, "grid %u, closed %d, %u faces", r.grid[0], r.closed, qaws_mesh_3d_get_face_count(out));
		ms_panel(c, 550, 0, 550, 640, out, 3.3, 0.35, "made watertight", sub);
		ms_draw(c, out, warm, NULL, 0, wire, 0.4);
		printf("  sphere: %s, %u evaluated\n", sub, r.evaluated);
	}
	qaws_mesh_3d_destroy(in);
	qaws_mesh_3d_destroy(out);
	/* two cubes crossing */
	{
		static unsigned int const bt[36] = { 0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 2, 3, 7, 2, 7, 6, 0, 4, 7, 0, 7, 3, 1, 2, 6, 1, 6, 5 };
		static int const bx[24] = { 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1 };
		static double const lo[2][3] = { { 0, 0, 0 }, { 0.4, 0.3, 0.2 } };
		qaws_scalar bxs[48];
		unsigned int bts[72], b, k;
		for (b = 0; b < 2; b++)
		{
			for (i = 0; i < 8; i++)
				for (k = 0; k < 3; k++)
					bxs[3 * (8 * b + i) + k] = (qaws_scalar)(lo[b][k] + bx[3 * i + k]);
			for (i = 0; i < 36; i++)
				bts[36 * b + i] = 8 * b + bt[i];
		}
		out = ms_wt_make(bxs, 16, bts, 24, 40, 0, &r);
		if (out)
		{
			sprintf(sub, "their union, %u sharp cells", r.sharp);
			ms_panel(c, 1100, 0, 550, 640, out, 0.6, 0.45, "two cubes crossing", sub);
			ms_draw(c, out, warm, NULL, 0, wire, 0.4);
			printf("  cubes: %s, closed %d\n", sub, r.closed);
		}
		qaws_mesh_3d_destroy(out);
	}
	/* an open wavy sheet, wrapped at 0.06 */
	{
		unsigned int n = 16, sv = (n + 1) * (n + 1), st = 0;
		qaws_scalar* sx = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * sv);
		unsigned int* sf = (unsigned int*)malloc(sizeof(unsigned int) * 6 * n * n);
		for (j = 0; j <= n; j++)
			for (i = 0; i <= n; i++)
			{
				double u = (double)i / n, v = (double)j / n;
				sx[3 * (j * (n + 1) + i)] = (qaws_scalar)u;
				sx[3 * (j * (n + 1) + i) + 1] = (qaws_scalar)v;
				sx[3 * (j * (n + 1) + i) + 2] = (qaws_scalar)(0.15 * sin(2 * MS_PI * u) * cos(MS_PI * v));
			}
		for (j = 0; j < n; j++)
			for (i = 0; i < n; i++)
			{
				unsigned int a = j * (n + 1) + i;
				sf[st++] = a; sf[st++] = a + 1; sf[st++] = a + n + 2;
				sf[st++] = a; sf[st++] = a + n + 2; sf[st++] = a + n + 1;
			}
		out = ms_wt_make(sx, sv, sf, st / 3, 56, 0.06, &r);
		if (out)
		{
			sprintf(sub, "offset 0.06: a closed shell, %u faces", qaws_mesh_3d_get_face_count(out));
			ms_panel(c, 1650, 0, 550, 640, out, 0.5, 0.6, "an open sheet, wrapped", sub);
			ms_draw(c, out, warm, NULL, 0, wire, 0.4);
			printf("  sheet: %s, closed %d\n", sub, r.closed);
		}
		qaws_mesh_3d_destroy(out);
		free(sx);
		free(sf);
	}
	er_save(c, "showcase/mesh25_watertight.png");
	printf("mesh25_watertight.png\n");
	free(x);
	free(f);
	free(k2);
	er_destroy(c);
}
