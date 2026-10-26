/*
 * Figure 64: subdivision by the type of data (qaws_mesh_3d_subdivide).
 *   - a ring's phase (each vertex its angle round the centre) as LINEAR
 *     values and as an ANGLE, three Catmull-Clark levels, coloured by hue:
 *     LINEAR values run through every hue where the phase wraps from pi to
 *     -pi, the ANGLE goes the short way round;
 *   - edge values (one hue per edge of a cube) through Catmull-Clark, Loop
 *     (on the cube's triangles, the diagonals 0) and sqrt(3) (the cube's
 *     creases features: its other edges flip away): each edge's pieces keep
 *     its value, the new edges 0 (thin).
 * Printed: vertices, edges, valued edges.
 */

static void ms_sr_hue(double h, double* rgb)
{
	unsigned int k;
	for (k = 0; k < 3; k++)
	{
		double x = fmod(h * 3 / MS_PI + 12 - 2.0 * k, 6.0), w = x < 1 ? x : x < 3 ? 1 : x < 4 ? 4 - x : 0;
		rgb[k] = 255 * (0.25 + 0.7 * w);
	}
}

/* a mesh's faces flat in one colour, its edges by an edge attribute "w" (0: thin grey, k: hue k of n) */
static void ms_sr_edges(er_canvas* c, qaws_mesh_3d const* m, unsigned int n, unsigned int* valued)
{
	static double const face[3] = { 236, 234, 228 }, thin[3] = { 150, 150, 150 };
	unsigned int ne = qaws_mesh_3d_get_edge_count(m), i, el;
	unsigned int* e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (ne + 1));
	double *x = ms_positions(m), *w = ms_at_get(m, "w", &el);
	ms_draw(c, m, face, NULL, 0, NULL, 0);
	qaws_mesh_3d_get_edges(e, m);
	*valued = 0;
	for (i = 0; w && i < ne; i++)
		if (w[i] == 0)
			er_line(c, x + 3 * e[2 * i], x + 3 * e[2 * i + 1], thin, 1);
	for (i = 0; w && i < ne; i++)
		if (w[i] != 0)
		{
			double rgb[3];
			ms_sr_hue(2 * MS_PI * (w[i] - 1) / n, rgb);
			er_line(c, x + 3 * e[2 * i], x + 3 * e[2 * i + 1], rgb, 4);
			(*valued)++;
		}
	free(e); free(x); free(w);
}

/* a cube (quads, or triangles), an edge value 1 + i on its box edges in order, 0 on the diagonals */
static qaws_mesh_3d* ms_sr_cube(int tris)
{
	static unsigned int const q[24] = { 0, 3, 2, 1, 4, 5, 6, 7, 0, 1, 5, 4, 2, 3, 7, 6, 1, 2, 6, 5, 0, 4, 7, 3 };
	static unsigned int const s[6] = { 4, 4, 4, 4, 4, 4 };
	double x[24];
	unsigned int t[36], i, ne, k = 0, *e;
	double* w;
	qaws_mesh_3d* m;
	for (i = 0; i < 8; i++)
	{
		x[3 * i] = (i == 1 || i == 2 || i == 5 || i == 6) ? 1 : -1;
		x[3 * i + 1] = (i == 2 || i == 3 || i == 6 || i == 7) ? 1 : -1;
		x[3 * i + 2] = i >= 4 ? 1 : -1;
	}
	for (i = 0; i < 6; i++)
	{
		t[6 * i] = q[4 * i]; t[6 * i + 1] = q[4 * i + 1]; t[6 * i + 2] = q[4 * i + 2];
		t[6 * i + 3] = q[4 * i]; t[6 * i + 4] = q[4 * i + 2]; t[6 * i + 5] = q[4 * i + 3];
	}
	m = tris ? ms_make(x, 8, t, NULL, 12) : ms_make(x, 8, q, s, 6);
	ne = qaws_mesh_3d_get_edge_count(m);
	e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne);
	w = (double*)malloc(sizeof(double) * ne);
	qaws_mesh_3d_get_edges(e, m);
	for (i = 0; i < ne; i++)
	{
		double const *a = x + 3 * e[2 * i], *b = x + 3 * e[2 * i + 1];
		int same = (a[0] == b[0]) + (a[1] == b[1]) + (a[2] == b[2]);
		w[i] = same == 2 ? ++k : 0;
	}
	ms_at_add(m, "w", QAWS_MESH_EDGE, 1, QAWS_ATTRIBUTE_LINEAR, w, ne);
	free(e); free(w);
	return m;
}

static void ms_sr_cube_panel(er_canvas* c, int px, qaws_mesh_subdivide_scheme scheme, unsigned int levels, int tris, int features, char const* title)
{
	qaws_mesh_3d *m = ms_sr_cube(tris), *o = NULL;
	qaws_mesh_subdivide_desc sd;
	unsigned int valued = 0;
	char sub[200];
	if (features)
		qaws_mesh_3d_detect_features(m, (qaws_scalar)0.5, 0, 0);
	memset(&sd, 0, sizeof sd);
	sd.scheme = scheme;
	sd.levels = levels;
	if (m && qaws_mesh_3d_subdivide(&o, m, &sd) == QAWS_STATUS_OK && o)
	{
		ms_panel(c, px, 0, 600, 640, o, -0.6, 0.6, title, NULL);
		ms_sr_edges(c, o, 12, &valued);
		sprintf(sub, "%u vertices, %u edges, %u valued", qaws_mesh_3d_get_vertex_count(o), qaws_mesh_3d_get_edge_count(o), valued);
		er_text(c, px + 14, 46, sub, 2, 95, 95, 95);
		printf("  %s: %s\n", title, sub);
	}
	qaws_mesh_3d_destroy(m);
	qaws_mesh_3d_destroy(o);
}

static void demo_subdivide_rules(void)
{
	er_canvas* c = er_create(2400, 640, 3, 252, 252, 250);
	if (!c)
		return;
	{
		/* a ring of 12 quads, its phase as LINEAR values (left) and as an ANGLE (right) */
		enum { N = 12 };
		double x[4 * N * 3 + 3], ph[2 * N];
		unsigned int f[4 * N], s[N], i, r, side;
		qaws_mesh_3d* o[2] = { NULL, NULL };
		double* pts = NULL;
		unsigned int np = 0;
		qaws_mesh_subdivide_desc sd;
		memset(&sd, 0, sizeof sd);
		sd.scheme = QAWS_MESH_SUBDIVIDE_CATMULL_CLARK;
		sd.levels = 3;
		for (side = 0; side < 2; side++)
		{
			qaws_mesh_3d* ring;
			for (i = 0; i < N; i++)
				for (r = 0; r < 2; r++)
				{
					double a = 2 * MS_PI * i / N, rad = r ? 1.0 : 0.55;
					x[3 * (2 * i + r)] = rad * cos(a) + (side ? 1.25 : -1.25);
					x[3 * (2 * i + r) + 1] = rad * sin(a);
					x[3 * (2 * i + r) + 2] = 0;
					ph[2 * i + r] = atan2(sin(a), cos(a));
				}
			for (i = 0; i < N; i++)
			{
				unsigned int j = (i + 1) % N;
				f[4 * i] = 2 * i; f[4 * i + 1] = 2 * j; f[4 * i + 2] = 2 * j + 1; f[4 * i + 3] = 2 * i + 1;
				s[i] = 4;
			}
			ring = ms_make(x, 2 * N, f, s, N);
			ms_at_add(ring, "phase", QAWS_MESH_VERTEX, 1, side ? QAWS_ATTRIBUTE_ANGLE : QAWS_ATTRIBUTE_LINEAR, ph, 2 * N);
			qaws_mesh_3d_subdivide(&o[side], ring, &sd);
			qaws_mesh_3d_destroy(ring);
		}
		if (o[0] && o[1])
		{
			double *p0 = ms_positions(o[0]), *p1 = ms_positions(o[1]);
			unsigned int n0 = qaws_mesh_3d_get_vertex_count(o[0]), n1 = qaws_mesh_3d_get_vertex_count(o[1]);
			pts = (double*)malloc(sizeof(double) * 3 * (n0 + n1));
			memcpy(pts, p0, sizeof(double) * 3 * n0);
			memcpy(pts + 3 * n0, p1, sizeof(double) * 3 * n1);
			np = n0 + n1;
			free(p0); free(p1);
			er_panel(c, 0, 72, 1200, 640 - 72, pts, np, 0, 1.45, 0.08);
			er_text(c, 14, 12, "a ring's phase subdivided: LINEAR, ANGLE", 3, 30, 30, 30);
			er_text(c, 14, 46, "left: through every hue where the phase wraps; right: the short way round", 2, 95, 95, 95);
			for (side = 0; side < 2; side++)
			{
				unsigned int nv = qaws_mesh_3d_get_vertex_count(o[side]), el;
				double *v = ms_at_get(o[side], "phase", &el), *col = (double*)malloc(sizeof(double) * 3 * nv);
				for (i = 0; v && i < nv; i++) ms_sr_hue(v[i], col + 3 * i);
				if (v)
					ms_draw(c, o[side], NULL, col, 0, NULL, 0);
				free(v); free(col);
			}
			printf("  ring: %u vertices per side\n", qaws_mesh_3d_get_vertex_count(o[0]));
		}
		free(pts);
		qaws_mesh_3d_destroy(o[0]); qaws_mesh_3d_destroy(o[1]);
	}
	ms_sr_cube_panel(c, 1200, QAWS_MESH_SUBDIVIDE_CATMULL_CLARK, 2, 0, 0, "edge values, Catmull-Clark");
	ms_sr_cube_panel(c, 1800, QAWS_MESH_SUBDIVIDE_LOOP, 2, 1, 0, "edge values, Loop");
	er_save(c, "showcase/mesh64_subdivide_rules.png");
	printf("mesh64_subdivide_rules.png\n");
	er_destroy(c);
	{
		/* sqrt(3) printed only (its creases kept as features) */
		qaws_mesh_3d *m = ms_sr_cube(1), *o = NULL;
		qaws_mesh_subdivide_desc sd;
		memset(&sd, 0, sizeof sd);
		sd.scheme = QAWS_MESH_SUBDIVIDE_SQRT3;
		sd.levels = 2;
		qaws_mesh_3d_detect_features(m, (qaws_scalar)0.5, 0, 0);
		if (qaws_mesh_3d_subdivide(&o, m, &sd) == QAWS_STATUS_OK && o)
		{
			unsigned int ne = qaws_mesh_3d_get_edge_count(o), i, el, n = 0;
			double* w = ms_at_get(o, "w", &el);
			for (i = 0; w && i < ne; i++) n += w[i] != 0;
			printf("  sqrt(3) x2, creases as features: %u edges, %u valued\n", ne, n);
			free(w);
		}
		qaws_mesh_3d_destroy(m); qaws_mesh_3d_destroy(o);
	}
}
