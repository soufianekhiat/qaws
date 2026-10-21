/*
 * Figure 33: what decimation keeps (qaws_mesh_3d_decimate's options).
 *   - a sphere rippled round its axis, seen down the axis, to 250 faces:
 *     plainly and keeping the silhouette seen along the axis (the input's
 *     outline drawn over each in red); with 8 x 8 octahedral views from
 *     above, seen side on: the faces stay where the views see;
 *   - a can (a cylinder with flat lids) to 60 faces: plainly the rims are
 *     cut; with creases found by angle they stay;
 *   - a rippled torus to 600 faces with importance 30 on its right half:
 *     the detail kept where it matters.
 */

/* the sphere: an icosphere (level 4) rippled 12 times round z */
static qaws_mesh_3d* ms_dv_rippled(void)
{
	static int const ic[60] = { 0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11, 1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8, 3, 9, 4, 3, 4, 2,
		3, 2, 6, 3, 6, 8, 3, 8, 9, 4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7, 9, 8, 1 };
	double p = (1 + sqrt(5.0)) / 2;
	double const v0[36] = { -1, p, 0, 1, p, 0, -1, -p, 0, 1, -p, 0, 0, -1, p, 0, 1, p, 0, -1, -p, 0, 1, -p, p, 0, -1, p, 0, 1, -p, 0, -1, -p, 0, 1 };
	unsigned int nt = 20, nv = 12, i, j, lv, nm;
	double* x = (double*)malloc(sizeof(double) * 3 * 2600);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * 5120);
	unsigned int* t2 = (unsigned int*)malloc(sizeof(unsigned int) * 3 * 5120);
	unsigned int(*mid)[3] = (unsigned int(*)[3])malloc(sizeof(unsigned int) * 3 * 8000);
	qaws_mesh_3d* m;
	for (i = 0; i < 12; i++)
	{
		double l = sqrt(v0[3 * i] * v0[3 * i] + v0[3 * i + 1] * v0[3 * i + 1] + v0[3 * i + 2] * v0[3 * i + 2]);
		for (j = 0; j < 3; j++)
			x[3 * i + j] = v0[3 * i + j] / l;
	}
	for (i = 0; i < 60; i++)
		t[i] = (unsigned int)ic[i];
	for (lv = 0; lv < 4; lv++)
	{
		nm = 0;
		for (i = 0; i < nt; i++)
		{
			unsigned int c[3], e;
			for (e = 0; e < 3; e++)
			{
				unsigned int a = t[3 * i + e], b = t[3 * i + (e + 1) % 3], lo = a < b ? a : b, hi = a < b ? b : a, s;
				double l;
				for (s = 0; s < nm && !(mid[s][0] == lo && mid[s][1] == hi); s++)
					;
				if (s == nm)
				{
					mid[nm][0] = lo;
					mid[nm][1] = hi;
					mid[nm++][2] = nv;
					for (j = 0; j < 3; j++)
						x[3 * nv + j] = 0.5 * (x[3 * a + j] + x[3 * b + j]);
					l = sqrt(x[3 * nv] * x[3 * nv] + x[3 * nv + 1] * x[3 * nv + 1] + x[3 * nv + 2] * x[3 * nv + 2]);
					for (j = 0; j < 3; j++)
						x[3 * nv + j] /= l;
					nv++;
				}
				c[e] = mid[s][2];
			}
			{
				unsigned int const q[12] = { t[3 * i], c[0], c[2], c[0], t[3 * i + 1], c[1], c[2], c[1], t[3 * i + 2], c[0], c[1], c[2] };
				memcpy(t2 + 12 * i, q, sizeof q);
			}
		}
		nt *= 4;
		memcpy(t, t2, sizeof(unsigned int) * 3 * nt);
	}
	for (i = 0; i < nv; i++)
	{
		double r = 1 + 0.07 * sin(12 * atan2(x[3 * i + 1], x[3 * i]));
		for (j = 0; j < 3; j++)
			x[3 * i + j] *= r;
	}
	m = ms_make(x, nv, t, NULL, nt);
	free(x); free(t); free(t2); free(mid);
	return m;
}

/* the input's outline seen along z: its silhouette edges, drawn above everything */
static void ms_dv_outline(er_canvas* c, qaws_mesh_3d const* m, double const* col)
{
	unsigned int nt = qaws_mesh_3d_get_triangle_count(m), f, k;
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nt);
	double* x = ms_positions(m);
	double* nz = (double*)malloc(sizeof(double) * nt);
	qaws_mesh_3d_get_triangles(t, m);
	for (f = 0; f < nt; f++)
	{
		double const *a = x + 3 * t[3 * f], *b = x + 3 * t[3 * f + 1], *e = x + 3 * t[3 * f + 2];
		nz[f] = (b[0] - a[0]) * (e[1] - a[1]) - (b[1] - a[1]) * (e[0] - a[0]);
	}
	/* (a pair of faces across an edge turned opposite ways along z: found by brute force, the mesh is small) */
	for (f = 0; f < nt; f++)
		for (k = 0; k < 3; k++)
		{
			unsigned int a = t[3 * f + k], b = t[3 * f + (k + 1) % 3], g, q;
			if (nz[f] <= 0)
				continue;
			for (g = 0; g < nt; g++)
				for (q = 0; q < 3; q++)
					if (t[3 * g + q] == b && t[3 * g + (q + 1) % 3] == a && nz[g] <= 0)
					{
						double A[3] = { x[3 * a], x[3 * a + 1], 5 }, B[3] = { x[3 * b], x[3 * b + 1], 5 };
						er_line(c, A, B, col, 2.4);
					}
		}
	free(t);
	free(x);
	free(nz);
}

/* a can: a cylinder of radius 0.5, height 1, nu round and nv along, with lids fanned to their centres */
static qaws_mesh_3d* ms_dv_can(unsigned int nu, unsigned int nv)
{
	unsigned int n = nu * (nv + 1) + 2, i, j, k = 0;
	double* x = (double*)malloc(sizeof(double) * 3 * n);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (2 * nu * nv + 2 * nu));
	qaws_mesh_3d* m;
	for (j = 0; j <= nv; j++)
		for (i = 0; i < nu; i++)
		{
			double a = 2 * MS_PI * i / nu;
			x[3 * (j * nu + i)] = 0.5 * cos(a);
			x[3 * (j * nu + i) + 1] = 0.5 * sin(a);
			x[3 * (j * nu + i) + 2] = (double)j / nv - 0.5;
		}
	x[3 * (n - 2)] = x[3 * (n - 2) + 1] = 0;
	x[3 * (n - 2) + 2] = -0.5;
	x[3 * (n - 1)] = x[3 * (n - 1) + 1] = 0;
	x[3 * (n - 1) + 2] = 0.5;
	for (j = 0; j < nv; j++)
		for (i = 0; i < nu; i++)
		{
			unsigned int a = j * nu + i, b = j * nu + (i + 1) % nu;
			t[k++] = a; t[k++] = b; t[k++] = b + nu;
			t[k++] = a; t[k++] = b + nu; t[k++] = a + nu;
		}
	for (i = 0; i < nu; i++)
	{
		t[k++] = n - 2; t[k++] = (i + 1) % nu; t[k++] = i;
		t[k++] = n - 1; t[k++] = nv * nu + i; t[k++] = nv * nu + (i + 1) % nu;
	}
	m = ms_make(x, n, t, NULL, k / 3);
	free(x);
	free(t);
	return m;
}

static void demo_decimate_views(void)
{
	static double const col[3] = { 180, 190, 214 }, warm[3] = { 226, 196, 150 }, wire[3] = { 60, 66, 80 }, red[3] = { 210, 40, 40 };
	static double const keep[3] = { 236, 150, 90 };
	static char const* const sname[4] = { "a rippled sphere, from above", "250 faces, plainly", "its outline from above kept",
		"8 x 8 views from above" };
	er_canvas* c = er_create(2400, 1440, 3, 252, 252, 250);
	qaws_mesh_3d* in = ms_dv_rippled();
	char sub[200];
	unsigned int k;
	if (!c)
		return;
	for (k = 0; k < 4; k++)
	{
		qaws_mesh_3d* m = ms_dv_rippled();
		qaws_decimate_desc d;
		qaws_decimate_report r;
		qaws_scalar const z[3] = { 0, 0, 1 };
		memset(&d, 0, sizeof d);
		memset(&r, 0, sizeof r);
		d.target_faces = 250;
		if (k == 2)
		{
			d.views = z;
			d.view_count = 1;
		}
		if (k == 3)
		{
			d.octahedral_views = 8;
			d.flags = QAWS_DECIMATE_HEMISPHERE;
		}
		if (k)
			qaws_mesh_3d_decimate(&r, m, &d);
		sprintf(sub, k == 3 ? "%u faces, seen side on: the bottom goes" : k ? "%u faces, %u silhouette edges kept in" : "%u faces", k ? r.faces_after : qaws_mesh_3d_get_face_count(m), r.silhouette_edges);
		ms_panel(c, 600 * (int)k, 0, 600, 720, in, 0.0, k == 3 ? 0.15 : 1.5707, sname[k], sub);
		ms_draw(c, m, col, NULL, 0, wire, k ? 0.6 : 0.15);
		if (k && k < 3)
			ms_dv_outline(c, in, red);
		printf("  rippled sphere, %s: %s\n", sname[k], sub);
		qaws_mesh_3d_destroy(m);
	}
	qaws_mesh_3d_destroy(in);
	for (k = 0; k < 2; k++)
	{
		qaws_mesh_3d* m = ms_dv_can(48, 12);
		qaws_decimate_desc d;
		qaws_decimate_report r;
		memset(&d, 0, sizeof d);
		d.target_faces = 60;
		d.crease_angle = (qaws_scalar)(k ? 0.6 : 0);
		qaws_mesh_3d_decimate(&r, m, &d);
		sprintf(sub, k ? "%u faces, %u creases found and kept" : "%u faces: the rims are cut", r.faces_after, r.creases);
		ms_panel(c, 600 * (int)k, 720, 600, 720, m, 0.5, 0.45, k ? "creases found by angle" : "a can to 60 faces, plainly", sub);
		ms_draw(c, m, warm, NULL, 0, wire, 0.8);
		printf("  can, %s\n", sub);
		qaws_mesh_3d_destroy(m);
	}
	{
		/* a rippled torus, its right half 30 times as important */
		enum { NU = 120, NV = 48 };
		double* x = (double*)malloc(sizeof(double) * 3 * NU * NV);
		unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 6 * NU * NV), i, j;
		qaws_scalar* imp = (qaws_scalar*)malloc(sizeof(qaws_scalar) * NU * NV);
		qaws_mesh_3d* m;
		qaws_decimate_desc d;
		qaws_decimate_report r;
		double* rgb = (double*)malloc(sizeof(double) * 3 * NU * NV);
		ms_torus_grid(NU, NV, x, t);
		for (i = 0; i < NU; i++)
			for (j = 0; j < NV; j++)
			{
				double u = 2 * MS_PI * i / NU, v = 2 * MS_PI * j / NV, a = 0.05 * sin(9 * u) * (cos(v) > 0 ? cos(v) : 0);
				double* p = x + 3 * (i * NV + j);
				p[0] += a * cos(v) * cos(u);
				p[1] += a * cos(v) * sin(u);
				p[2] += a * sin(v);
				imp[i * NV + j] = (qaws_scalar)(p[0] > 0 ? 30 : 1);
			}
		m = ms_make(x, NU * NV, t, NULL, 2 * NU * NV);
		memset(&d, 0, sizeof d);
		d.target_faces = 600;
		d.importance = imp;
		qaws_mesh_3d_decimate(&r, m, &d);
		for (i = 0; i < qaws_mesh_3d_get_vertex_count(m); i++)
			memcpy(rgb + 3 * i, qaws_mesh_3d_get_vertex(m, i).x > 0 ? keep : col, sizeof col);
		sprintf(sub, "%u faces, the right half weighs 30 times more", r.faces_after);
		ms_panel(c, 1200, 720, 1200, 720, m, -0.3, 0.9, "importance: detail kept where it matters", sub);
		ms_draw(c, m, NULL, rgb, 0, wire, 0.6);
		printf("  torus with importance: %s\n", sub);
		qaws_mesh_3d_destroy(m);
		free(x); free(t); free(imp); free(rgb);
	}
	er_save(c, "showcase/mesh33_decimate_views.png");
	printf("mesh33_decimate_views.png\n");
	er_destroy(c);
}
