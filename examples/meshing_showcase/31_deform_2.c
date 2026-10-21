/*
 * Figure 31: the ARAP paper's own tests (Sorkine, Alexa 2007, figures 2
 * and 6) on qaws_mesh_3d_deform.
 *   - a plate with 4 x 4 bumps, one side held and the opposite one lifted
 *     and drawn in: as rigidly as possible the bumps stay upright on the
 *     bent plate; by Laplacian editing they keep their rest direction and
 *     lean; the biharmonic bend stretches them;
 *   - a bar of square section held at its foot, its top turned a quarter turn:
 *     ARAP twists it keeping its section; the linear methods pinch it.
 * Printed: the edge lengths' change.
 */

enum { MS_D2_N = 48 };

/* the plate: [0, 1] x [0, 1] in N x N cells, 4 x 4 bumps */
static qaws_mesh_3d* ms_d2_plate(void)
{
	enum { N = MS_D2_N };
	double* x = (double*)malloc(sizeof(double) * 3 * (N + 1) * (N + 1));
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 6 * N * N), i, j, k = 0, a, b;
	qaws_mesh_3d* m;
	for (j = 0; j <= N; j++)
		for (i = 0; i <= N; i++)
		{
			double u = (double)i / N, v = (double)j / N, z = 0;
			for (a = 0; a < 4; a++)
				for (b = 0; b < 4; b++)
				{
					double cu = 0.2 + 0.2 * a, cv = 0.2 + 0.2 * b, r2 = (u - cu) * (u - cu) + (v - cv) * (v - cv);
					z += 0.09 * exp(-r2 / 0.0025);
				}
			x[3 * (j * (N + 1) + i)] = u;
			x[3 * (j * (N + 1) + i) + 1] = v;
			x[3 * (j * (N + 1) + i) + 2] = z;
		}
	for (j = 0; j < N; j++)
		for (i = 0; i < N; i++)
		{
			unsigned int c0 = j * (N + 1) + i;
			t[k++] = c0; t[k++] = c0 + 1; t[k++] = c0 + N + 2;
			t[k++] = c0; t[k++] = c0 + N + 2; t[k++] = c0 + N + 1;
		}
	m = ms_make(x, (N + 1) * (N + 1), t, NULL, 2 * N * N);
	free(x);
	free(t);
	return m;
}

/* a bar of square section (side 0.3, nu round, nv along, height 1.6), open at its ends */
static qaws_mesh_3d* ms_d2_bar(unsigned int nu, unsigned int nv)
{
	double* x = (double*)malloc(sizeof(double) * 3 * nu * (nv + 1));
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 6 * nu * nv), i, j, k = 0;
	qaws_mesh_3d* m;
	for (j = 0; j <= nv; j++)
		for (i = 0; i < nu; i++)
		{
			/* round the square's perimeter, a side per quarter */
			double s = 4.0 * i / nu, f = s - floor(s), px, py;
			switch ((int)s)
			{
			case 0: px = -0.15 + 0.3 * f; py = -0.15; break;
			case 1: px = 0.15; py = -0.15 + 0.3 * f; break;
			case 2: px = 0.15 - 0.3 * f; py = 0.15; break;
			default: px = -0.15; py = 0.15 - 0.3 * f; break;
			}
			x[3 * (j * nu + i)] = px;
			x[3 * (j * nu + i) + 1] = py;
			x[3 * (j * nu + i) + 2] = 1.6 * j / nv;
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

static double ms_d2_lengths(qaws_mesh_3d const* m, qaws_mesh_3d const* rest)
{
	unsigned int ne = qaws_mesh_3d_get_edge_count(m), i;
	unsigned int* e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne);
	double s = 0;
	qaws_mesh_3d_get_edges(e, m);
	for (i = 0; i < ne; i++)
	{
		qaws_vec3 a = qaws_mesh_3d_get_vertex(m, e[2 * i]), b = qaws_mesh_3d_get_vertex(m, e[2 * i + 1]);
		qaws_vec3 a0 = qaws_mesh_3d_get_vertex(rest, e[2 * i]), b0 = qaws_mesh_3d_get_vertex(rest, e[2 * i + 1]);
		double l = sqrt((double)(a.x - b.x) * (a.x - b.x) + (double)(a.y - b.y) * (a.y - b.y) + (double)(a.z - b.z) * (a.z - b.z));
		double l0 = sqrt((double)(a0.x - b0.x) * (a0.x - b0.x) + (double)(a0.y - b0.y) * (a0.y - b0.y) + (double)(a0.z - b0.z) * (a0.z - b0.z));
		s += fabs(l / l0 - 1);
	}
	free(e);
	return s / ne;
}

static void demo_deform_2(void)
{
	static double const col[3] = { 180, 190, 214 }, wire[3] = { 70, 75, 90 }, hand[3] = { 220, 60, 40 };
	static char const* const pname[4] = { "a plate with 4 x 4 bumps", "ARAP: one side lifted", "Laplacian editing", "biharmonic" };
	static char const* const bname[4] = { "a bar of square section", "ARAP: twisted", "Laplacian editing", "biharmonic" };
	er_canvas* c = er_create(2400, 1440, 3, 252, 252, 250);
	char sub[200];
	unsigned int k;
	if (!c)
		return;
	{
		/* the plate: its left two columns held, its right two lifted and drawn in */
		enum { N = MS_D2_N };
		qaws_mesh_3d* rest = ms_d2_plate();
		unsigned int h[4 * (N + 1)], nh = 0, i, j;
		qaws_scalar tg[12 * (N + 1)];
		for (j = 0; j <= N; j++)
			for (i = 0; i <= N; i++)
				if (i <= 1 || i >= N - 1)
				{
					qaws_vec3 p = qaws_mesh_3d_get_vertex(rest, j * (N + 1) + i);
					h[nh] = j * (N + 1) + i;
					tg[3 * nh] = (qaws_scalar)(i <= 1 ? p.x : p.x - 0.12);
					tg[3 * nh + 1] = p.y;
					tg[3 * nh + 2] = (qaws_scalar)(i <= 1 ? p.z : p.z + 0.5);
					nh++;
				}
		for (k = 0; k < 4; k++)
		{
			qaws_mesh_3d* m = ms_d2_plate();
			qaws_deform_desc d;
			double* rgb = (double*)malloc(sizeof(double) * 3 * (N + 1) * (N + 1));
			memset(&d, 0, sizeof d);
			d.method = k == 1 ? QAWS_DEFORM_ARAP : k == 2 ? QAWS_DEFORM_LAPLACIAN : QAWS_DEFORM_BIHARMONIC;
			if (k)
				qaws_mesh_3d_deform(m, h, tg, nh, &d);
			for (i = 0; i < (N + 1) * (N + 1); i++)
				memcpy(rgb + 3 * i, (i % (N + 1)) <= 1 || (i % (N + 1)) >= N - 1 ? hand : col, sizeof col);
			if (k)
				sprintf(sub, "edge lengths changed by %.1f%% on average", 100 * ms_d2_lengths(m, rest));
			else
				sprintf(sub, "red: held (left), lifted (right)");
			ms_panel(c, 600 * (int)k, 0, 600, 720, rest, -0.25, 0.3, pname[k], sub);
			ms_draw(c, m, NULL, rgb, 1, wire, 0.2);
			printf("  plate, %s: %s\n", pname[k], sub);
			qaws_mesh_3d_destroy(m);
			free(rgb);
		}
		qaws_mesh_3d_destroy(rest);
	}
	for (k = 0; k < 4; k++)
	{
		/* the bar: its foot held, its top turned half a turn about its axis */
		enum { NU = 32, NV = 48 };
		qaws_mesh_3d* m = ms_d2_bar(NU, NV);
		qaws_mesh_3d* rest = ms_d2_bar(NU, NV);
		unsigned int h[2 * NU], i;
		qaws_scalar tg[6 * NU];
		qaws_deform_desc d;
		for (i = 0; i < NU; i++)
		{
			qaws_vec3 p = qaws_mesh_3d_get_vertex(m, i), q = qaws_mesh_3d_get_vertex(m, NV * NU + i);
			h[i] = i;
			tg[3 * i] = p.x;
			tg[3 * i + 1] = p.y;
			tg[3 * i + 2] = p.z;
			h[NU + i] = NV * NU + i;
			tg[3 * (NU + i)] = -q.y;
			tg[3 * (NU + i) + 1] = q.x;
			tg[3 * (NU + i) + 2] = q.z;
		}
		memset(&d, 0, sizeof d);
		d.method = k == 1 ? QAWS_DEFORM_ARAP : k == 2 ? QAWS_DEFORM_LAPLACIAN : QAWS_DEFORM_BIHARMONIC;
		if (k)
			qaws_mesh_3d_deform(m, h, tg, 2 * NU, &d);
		if (k)
			sprintf(sub, "edge lengths changed by %.1f%%", 100 * ms_d2_lengths(m, rest));
		else
			sprintf(sub, "its foot held, its top to turn a quarter turn");
		ms_panel(c, 600 * (int)k, 720, 600, 720, rest, 0.35, 0.25, bname[k], sub);
		ms_draw(c, m, col, NULL, 0, wire, 0.4);
		printf("  bar, %s: %s\n", bname[k], sub);
		qaws_mesh_3d_destroy(m);
		qaws_mesh_3d_destroy(rest);
	}
	er_save(c, "showcase/mesh31_deform_2.png");
	printf("mesh31_deform_2.png\n");
	er_destroy(c);
}
