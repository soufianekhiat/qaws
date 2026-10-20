/*
 * Figure 18: repair for farq.
 *   - a torus broken as exported meshes come: a soup (each triangle its own
 *     three points, off by 1e-7), three holes, three faces in ten turned
 *     over (drawn dark: their back shows), duplicates (some turned the
 *     other way), a quarter shifted 0.004 off the rest (a gap); the holes'
 *     borders red, the gap's orange;
 *   - qaws_mesh_repair (gap tolerance 0.01, the other policies by
 *     default): closed and outward, the fills tinted; seen from above and
 *     from below. Printed: the report.
 */

/* the torus (R 1, r 0.4) on a grid of nu by nv, triangles counter-clockwise from outside */
static void ms_torus_grid(unsigned int nu, unsigned int nv, double* x, unsigned int* t)
{
	unsigned int i, j;
	for (i = 0; i < nu; i++)
		for (j = 0; j < nv; j++)
		{
			double u = 2 * MS_PI * i / nu, v = 2 * MS_PI * j / nv;
			double* p = x + 3 * (i * nv + j);
			p[0] = (1 + 0.4 * cos(v)) * cos(u);
			p[1] = (1 + 0.4 * cos(v)) * sin(u);
			p[2] = 0.4 * sin(v);
		}
	for (i = 0; i < nu; i++)
		for (j = 0; j < nv; j++)
		{
			unsigned int a = i * nv + j, b = ((i + 1) % nu) * nv + j, c = ((i + 1) % nu) * nv + (j + 1) % nv, d = i * nv + (j + 1) % nv;
			unsigned int* q = t + 6 * (i * nv + j);
			q[0] = a; q[1] = b; q[2] = c;
			q[3] = a; q[4] = c; q[5] = d;
		}
}

static void demo_repair(void)
{
	enum { NU = 64, NV = 24, NT = 2 * NU * NV };
	static double const hole[3][3] = { { 0.3, 1.2, 0.26 }, { 2.0, 0.0, 0.16 }, { 4.2, 2.0, 0.34 } };   /* u, v, radius */
	static double const shift[3] = { 0.004, 0.0025, 0.003 }, wire[3] = { 60, 66, 80 };
	static double const base[3] = { 176, 186, 206 }, fillc[3] = { 236, 140, 60 }, red[3] = { 210, 30, 30 }, orange[3] = { 245, 150, 0 };
	double *x = (double*)malloc(sizeof(double) * 3 * NU * NV), *sx = NULL;
	unsigned int *t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * NT), *keep = (unsigned int*)malloc(sizeof(unsigned int) * NT);
	unsigned int *cnt = (unsigned int*)calloc(NU * NV, sizeof(unsigned int)), *st = NULL, *flip = NULL;
	unsigned int f, k, nk = 0, ns = 0, seed = 2024, i;
	qaws_scalar* sp = NULL;
	qaws_mesh_3d* m = NULL;
	qaws_mesh_repair_desc d;
	qaws_mesh_repair_report r;
	qaws_status status;
	er_canvas* c = er_create(2200, 700, 3, 252, 252, 250);
	char sub[220];
	clock_t t0;
	double secs;
	if (!c || !x || !t || !keep || !cnt)
		goto done;
	ms_torus_grid(NU, NV, x, t);
	/* the holes: faces whose centre is near a point of the torus */
	for (f = 0; f < NT; f++)
	{
		double cc[3] = { 0, 0, 0 };
		int out = 0;
		for (k = 0; k < 3; k++)
		{
			cc[0] += x[3 * t[3 * f + k]] / 3;
			cc[1] += x[3 * t[3 * f + k] + 1] / 3;
			cc[2] += x[3 * t[3 * f + k] + 2] / 3;
		}
		for (k = 0; k < 3; k++)
		{
			double p[3] = { (1 + 0.4 * cos(hole[k][1])) * cos(hole[k][0]), (1 + 0.4 * cos(hole[k][1])) * sin(hole[k][0]), 0.4 * sin(hole[k][1]) };
			out |= (cc[0] - p[0]) * (cc[0] - p[0]) + (cc[1] - p[1]) * (cc[1] - p[1]) + (cc[2] - p[2]) * (cc[2] - p[2]) < hole[k][2] * hole[k][2];
		}
		if (!out)
			keep[nk++] = f;
	}
	/* the soup: every kept face its own points, 30% turned, a quarter shifted, 40 copies and 6 copies turned over */
	sx = (double*)malloc(sizeof(double) * 9 * (nk + 46));
	st = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (nk + 46));
	flip = (unsigned int*)calloc(nk + 46, sizeof(unsigned int));
	if (!sx || !st || !flip)
		goto done;
	for (i = 0; i < nk + 46; i++)
	{
		unsigned int g = i < nk ? keep[i] : keep[(i - nk) * 37 % nk], q = t[3 * g] / NV;
		int moved = q >= NU / 4 && q < NU / 2;
		ms_rand(&seed);
		flip[i] = i < nk ? (ms_rand(&seed) % 10) < 3 : i >= nk + 40;
		for (k = 0; k < 3; k++)
		{
			unsigned int v = t[3 * g + (flip[i] && k ? 3 - k : k)], j;
			for (j = 0; j < 3; j++)
				sx[9 * i + 3 * k + j] = x[3 * v + j] + (moved ? shift[j] : 0) + 1e-7 * (ms_uniform(&seed) - 0.5);
			st[3 * i + k] = 3 * i + k;
			if (i < nk && !moved)
				cnt[v] |= 1;
			if (i < nk && moved)
				cnt[v] |= 2;
		}
	}
	ns = nk + 46;
	sp = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 9 * ns);
	if (!sp)
		goto done;
	for (i = 0; i < 9 * ns; i++)
		sp[i] = (qaws_scalar)sx[i];
	memset(&d, 0, sizeof d);
	d.gap_tolerance = (qaws_scalar)0.01;
	t0 = clock();
	status = qaws_mesh_repair(&m, &r, sp, st, NULL, 3 * ns, ns, &d);
	secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
	for (i = 0, k = 0; i < nk; i++)
		k += flip[i];
	printf("  a broken torus: %u faces given (%u turned, 46 copies, 3 holes, a quarter shifted 0.004)\n", ns, k);
	if (status != QAWS_STATUS_OK || !m)
	{
		printf("  repair failed (%d)\n", (int)status);
		goto done;
	}
	printf("  repaired in %.3f s: %u welded, %u duplicates, %u dropped, %u stitched, %u holes filled (%u faces), %u turned back, closed %d, "
		"%u boundary, %u non-manifold edges; %u vertices, %u faces\n", secs, r.welded_vertices, r.duplicate_faces, r.dropped_faces,
		r.stitched_vertices, r.holes_filled, r.fill_faces, r.flipped_faces, r.closed, r.boundary_edges, r.non_manifold_edges,
		qaws_mesh_3d_get_vertex_count(m), qaws_mesh_3d_get_face_count(m));
	/* the soup: flipped faces show their darker back; hole borders red, the gap's seam orange */
	er_panel(c, 0, 72, 733, 628, sx, 3 * ns, 0.5, 0.95, 0.04);
	er_text(c, 14, 12, "a broken soup", 3, 30, 30, 30);
	sprintf(sub, "%u faces, %u points, 3 in 10 turned, 46 copies", ns, 3 * ns);
	er_text(c, 14, 46, sub, 2, 95, 95, 95);
	for (i = 0; i < ns; i++)
		er_triangle(c, sx + 9 * i, sx + 9 * i + 3, sx + 9 * i + 6, NULL, NULL, NULL, base, base, base);
	for (f = 0; f < NT; f++)
		for (k = 0; k < 3; k++)
		{
			/* an edge of the grid on one kept face only (a hole's border), or between the shifted quarter and the rest */
			unsigned int a = t[3 * f + k], b = t[3 * f + (k + 1) % 3], j, n = 0;
			for (j = 0; j < nk && n < 3; j++)
			{
				unsigned int const* q = t + 3 * keep[j];
				n += (q[0] == a || q[1] == a || q[2] == a) && (q[0] == b || q[1] == b || q[2] == b);
			}
			if (n == 1 && a < b)
				er_line(c, x + 3 * a, x + 3 * b, red, 2.5);
			else if (n == 2 && a < b && cnt[a] == 3 && cnt[b] == 3)
				er_line(c, x + 3 * a, x + 3 * b, orange, 2.5);
		}
	/* the repaired mesh, fills tinted, from above and below */
	for (k = 0; k < 2; k++)
	{
		unsigned int nt = qaws_mesh_3d_get_face_count(m), nv = qaws_mesh_3d_get_vertex_count(m), a = 0;
		unsigned int* tt = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nt);
		qaws_scalar *src = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nt), *ns3 = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
		double *p = ms_positions(m), *n = (double*)malloc(sizeof(double) * 3 * nv);
		if (!tt || !src || !ns3 || !n)
		{
			free(tt); free(src); free(ns3); free(p); free(n);
			break;
		}
		qaws_mesh_3d_get_triangles(tt, m);
		qaws_mesh_3d_find_attribute(&a, m, "source_face");
		qaws_mesh_3d_get_attribute(src, m, a);
		qaws_mesh_3d_get_vertex_normals(ns3, m);
		for (i = 0; i < 3 * nv; i++)
			n[i] = (double)ns3[i];
		if (k == 0)
			sprintf(sub, "%u holes filled (%u faces), %u turned back", r.holes_filled, r.fill_faces, r.flipped_faces);
		else
			sprintf(sub, "%u stitched, %u duplicates removed, %.2f s", r.stitched_vertices, r.duplicate_faces, secs);
		ms_panel(c, 733 + 733 * (int)k, 0, 733, 700, m, 0.5, k ? -0.95 : 0.95, k ? "repaired, from below" : "repaired: closed, outward", sub);
		for (i = 0; i < nt; i++)
		{
			double const* col = (double)src[i] < 0 ? fillc : base;
			er_triangle(c, p + 3 * tt[3 * i], p + 3 * tt[3 * i + 1], p + 3 * tt[3 * i + 2], n + 3 * tt[3 * i], n + 3 * tt[3 * i + 1],
				n + 3 * tt[3 * i + 2], col, col, col);
		}
		ms_draw(c, m, NULL, NULL, 0, wire, 0.5);
		free(tt); free(src); free(ns3); free(p); free(n);
	}
	er_text(c, 14, 670, "hole borders red, the gap orange; fills tinted", 2, 95, 95, 95);
	er_save(c, "showcase/mesh18_repair.png");
done:
	qaws_mesh_3d_destroy(m);
	er_destroy(c);
	free(x); free(t); free(keep); free(cnt); free(sx); free(st); free(flip); free(sp);
}
