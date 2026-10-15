/* Figure 33: the surfaces' tangent-point energy minimized
   (qaws_repulsive_mesh_solve, R3). A bumpy sphere with its area and volume
   held, before and after; a sheet pressed under a pinned one (an
   obstacle), collision-checked. Flat-shaded triangles drawn far to near. */

/* the mesh's faces, shaded by their normals against a light, far first */
static void s33_mesh(svg* s, projection const* pr, qaws_scalar const* x, unsigned int const* tris, unsigned int nf, int const* rgb)
{
	unsigned int* order = (unsigned int*)malloc(sizeof(unsigned int) * nf);
	double* depth = (double*)malloc(sizeof(double) * nf);
	unsigned int i, k;
	static double const light[3] = { 0.35, -0.45, 0.82 };
	for (i = 0; i < nf; i++)
	{
		unsigned int const* t = tris + 3 * i;
		order[i] = i;
		depth[i] = (double)x[3 * t[0]] + x[3 * t[1]] + x[3 * t[2]] + x[3 * t[0] + 1] + x[3 * t[1] + 1] + x[3 * t[2] + 1] - 0.3 * ((double)x[3 * t[0] + 2] + x[3 * t[1] + 2] + x[3 * t[2] + 2]);
	}
	for (i = 1; i < nf; i++)
	{
		unsigned int key = order[i];
		k = i;
		while (k > 0 && depth[order[k - 1]] > depth[key]) { order[k] = order[k - 1]; k--; }
		order[k] = key;
	}
	for (i = 0; i < nf; i++)
	{
		unsigned int const* t = tris + 3 * order[i];
		double p[3][2], u[3], v[3], n[3], l, sh;
		unsigned int c, q;
		for (c = 0; c < 3; c++) project(pr, x[3 * t[c]], x[3 * t[c] + 1], x[3 * t[c] + 2], &p[c][0], &p[c][1]);
		for (q = 0; q < 3; q++) { u[q] = (double)x[3 * t[1] + q] - x[3 * t[0] + q]; v[q] = (double)x[3 * t[2] + q] - x[3 * t[0] + q]; }
		n[0] = u[1] * v[2] - u[2] * v[1]; n[1] = u[2] * v[0] - u[0] * v[2]; n[2] = u[0] * v[1] - u[1] * v[0];
		l = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
		sh = l > 0 ? 0.35 + 0.65 * fabs(n[0] * light[0] + n[1] * light[1] + n[2] * light[2]) / l : 0.35;
		fprintf(s->f, "<polygon points=\"%.2f,%.2f %.2f,%.2f %.2f,%.2f\" fill=\"rgb(%d,%d,%d)\" stroke=\"rgb(%d,%d,%d)\" stroke-width=\"0.4\"/>\n",
			p[0][0], p[0][1], p[1][0], p[1][1], p[2][0], p[2][1],
			(int)(rgb[0] * sh), (int)(rgb[1] * sh), (int)(rgb[2] * sh), (int)(rgb[0] * sh * 0.8), (int)(rgb[1] * sh * 0.8), (int)(rgb[2] * sh * 0.8));
	}
	free(order); free(depth);
}

/* an octahedron subdivided `level` times on the unit sphere */
static void s33_sphere(unsigned int level, double* x, unsigned int* tris, unsigned int* nv, unsigned int* nf)
{
	static double const ov[6][3] = { { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 }, { 0, 0, 1 }, { 0, 0, -1 } };
	static unsigned int const of[8][3] = { { 0, 2, 4 }, { 2, 1, 4 }, { 1, 3, 4 }, { 3, 0, 4 }, { 2, 0, 5 }, { 1, 2, 5 }, { 3, 1, 5 }, { 0, 3, 5 } };
	unsigned int l, i, n = 6, f = 8, q;
	unsigned int* t2 = (unsigned int*)malloc(sizeof(unsigned int) * 3 * 4 * 8 * (1u << (2 * level)));
	memcpy(x, ov, sizeof(ov));
	memcpy(tris, of, sizeof(of));
	for (l = 0; l < level; l++)
	{
		unsigned int nf2 = 0;
		for (i = 0; i < f; i++)
		{
			unsigned int c[3], mid[3], k;
			for (k = 0; k < 3; k++) c[k] = tris[3 * i + k];
			for (k = 0; k < 3; k++)
			{
				unsigned int a = c[k], b = c[(k + 1) % 3], j;
				double m[3], len;
				for (q = 0; q < 3; q++) m[q] = 0.5 * (x[3 * a + q] + x[3 * b + q]);
				len = sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]);
				for (q = 0; q < 3; q++) m[q] /= len;
				for (j = 0; j < n && !(fabs(x[3 * j] - m[0]) < 1e-12 && fabs(x[3 * j + 1] - m[1]) < 1e-12 && fabs(x[3 * j + 2] - m[2]) < 1e-12); j++)
					;
				if (j == n) { for (q = 0; q < 3; q++) x[3 * n + q] = m[q]; n++; }
				mid[k] = j;
			}
			t2[3 * nf2] = c[0]; t2[3 * nf2 + 1] = mid[0]; t2[3 * nf2 + 2] = mid[2]; nf2++;
			t2[3 * nf2] = mid[0]; t2[3 * nf2 + 1] = c[1]; t2[3 * nf2 + 2] = mid[1]; nf2++;
			t2[3 * nf2] = mid[2]; t2[3 * nf2 + 1] = mid[1]; t2[3 * nf2 + 2] = c[2]; nf2++;
			t2[3 * nf2] = mid[0]; t2[3 * nf2 + 1] = mid[1]; t2[3 * nf2 + 2] = mid[2]; nf2++;
		}
		memcpy(tris, t2, sizeof(unsigned int) * 3 * nf2);
		f = nf2;
	}
	free(t2);
	*nv = n; *nf = f;
}

static void demo_repulsive_mesh(void)
{
	clock_t t0 = clock();
	svg s;
	double const pw = 330, ph = 330;
	static int const blue[3] = { 84, 174, 255 }, grey[3] = { 200, 205, 210 }, green[3] = { 120, 200, 140 };
	char label[200];
	if (!svg_open(&s, "showcase/b2d33_repulsive_mesh.svg", (int)(48 + 3 * (pw + 12)), (int)(110 + ph + 30),
		"Repulsive surfaces",
		"qaws_repulsive_mesh_solve: the surfaces' tangent-point energy (p = 6) in the fractional Sobolev metric"))
		return;
	/* a bumpy sphere, before and after */
	{
		double xd[3 * 260];
		qaws_scalar x[3 * 260], x0[3 * 260];
		unsigned int tris[3 * 520], nv, nf, i, k;
		qaws_repulsive_mesh_desc d;
		qaws_repulsive_mesh* m = NULL;
		qaws_repulsive_mesh_solve_desc sd;
		qaws_repulsive_report r;
		s33_sphere(3, xd, tris, &nv, &nf);
		for (i = 0; i < nv; i++)
		{
			double rr = 1 + 0.22 * sin(4 * xd[3 * i]) * sin(3 * xd[3 * i + 2]) + 0.1 * cos(5 * xd[3 * i + 1]);
			x[3 * i] = (qaws_scalar)(xd[3 * i] * rr); x[3 * i + 1] = (qaws_scalar)(xd[3 * i + 1] * rr); x[3 * i + 2] = (qaws_scalar)(xd[3 * i + 2] * rr * 0.85);
		}
		memcpy(x0, x, sizeof(x));
		memset(&d, 0, sizeof(d));
		d.vertex_count = nv; d.triangles = tris; d.triangle_count = nf;
		memset(&sd, 0, sizeof(sd));
		sd.hold_area = 1; sd.hold_volume = 1; sd.max_iterations = 40;
		memset(&r, 0, sizeof(r));
		if (qaws_repulsive_mesh_create(&d, &m) == QAWS_STATUS_OK)
			qaws_repulsive_mesh_solve(m, &sd, x, &r);
		for (k = 0; k < 2; k++)
		{
			viewport v;
			projection pr;
			v.x0 = 24 + k * (pw + 12); v.y0 = 84; v.w = pw; v.h = ph; v.xmin = v.ymin = -1; v.xmax = v.ymax = 1;
			pr.cx = v.x0 + pw / 2; pr.cy = v.y0 + ph / 2 + 5; pr.scale = 95; pr.zscale = 1;
			sprintf(label, k ? "after %u steps: E %.2f (area, volume held)" : "a bumpy sphere (%u faces): E %.2f", k ? r.iterations : nf, k ? (double)r.energy : (double)r.initial_energy);
			svg_panel(&s, &v, label);
			s33_mesh(&s, &pr, k ? x : x0, tris, nf, k ? blue : grey);
		}
		qaws_repulsive_mesh_destroy(m);
	}
	/* a sheet under a pinned one */
	{
		enum { n = 9 };
		qaws_scalar x[3 * 2 * n * n], x0[3 * 2 * n * n];
		unsigned int tris[3 * 4 * (n - 1) * (n - 1)], pins[n * n + 4 * n], i, j, nf = 0, np = 0, sh;
		qaws_repulsive_mesh_desc d;
		qaws_repulsive_mesh* m = NULL;
		qaws_repulsive_mesh_solve_desc sd;
		qaws_repulsive_report r;
		viewport v;
		projection pr;
		for (sh = 0; sh < 2; sh++)
		{
			unsigned int base = sh * n * n;
			for (i = 0; i < n; i++)
				for (j = 0; j < n; j++)
				{
					double u = (double)i / (n - 1) - 0.5, w = (double)j / (n - 1) - 0.5;
					x[3 * (base + n * i + j)] = (qaws_scalar)u; x[3 * (base + n * i + j) + 1] = (qaws_scalar)w;
					x[3 * (base + n * i + j) + 2] = (qaws_scalar)(sh ? 0.22 : 0.15 * cos(3.1 * u) * cos(3.1 * w));
				}
			for (i = 0; i + 1 < n; i++)
				for (j = 0; j + 1 < n; j++)
				{
					unsigned int a = base + n * i + j;
					tris[3 * nf] = a; tris[3 * nf + 1] = a + n; tris[3 * nf + 2] = a + n + 1; nf++;
					tris[3 * nf] = a; tris[3 * nf + 1] = a + n + 1; tris[3 * nf + 2] = a + 1; nf++;
				}
		}
		for (i = 0; i < n * n; i++) pins[np++] = n * n + i;
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
				if (i == 0 || j == 0 || i + 1 == n || j + 1 == n) pins[np++] = n * i + j;   /* the membrane's rim */
		memcpy(x0, x, sizeof(x));
		memset(&d, 0, sizeof(d));
		d.vertex_count = 2 * n * n; d.triangles = tris; d.triangle_count = nf;
		memset(&sd, 0, sizeof(sd));
		sd.pinned = pins; sd.pinned_count = np; sd.hold_area = 1; sd.collision_checked = 1; sd.max_iterations = 30;
		memset(&r, 0, sizeof(r));
		if (qaws_repulsive_mesh_create(&d, &m) == QAWS_STATUS_OK)
			qaws_repulsive_mesh_solve(m, &sd, x, &r);
		v.x0 = 24 + 2 * (pw + 12); v.y0 = 84; v.w = pw; v.h = ph; v.xmin = v.ymin = -1; v.xmax = v.ymax = 1;
		pr.cx = v.x0 + pw / 2; pr.cy = v.y0 + ph / 2 - 10; pr.scale = 150; pr.zscale = 1.6;
		sprintf(label, "a membrane under a fixed sheet: E %.0f -> %.2f", (double)r.initial_energy, (double)r.energy);
		svg_panel(&s, &v, label);
		s33_mesh(&s, &pr, x, tris, nf / 2, blue);
		s33_mesh(&s, &pr, x, tris + 3 * (nf / 2), nf / 2, green);
		svg_text(&s, v.x0 + 10, v.y0 + ph - 10, 12, "#57606a", "start", "collision-checked; green fixed, blue rim pinned, area held");
		qaws_repulsive_mesh_destroy(m);
	}
	svg_close(&s);
	printf("figure 33: repulsive surfaces, %.1f s\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
}
