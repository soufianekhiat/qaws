/* Figure 24: geodesic disks through saddle vertices (qaws_shell_region.h).
   Faces round a vertex with more than 2 pi of angle: five squares (5 pi/2)
   and six (the monkey saddle, 3 pi), flat, then as bicubic patches bent by
   a smooth map (curved faces and curved edges). Rays from a centre part at
   the vertex; a disk through it is the plane disk and, behind the vertex,
   the disk about it of the radius left. On the bent monkey saddle, two
   disks through the vertex combined by every Boolean. The surfaces are
   shaded and drawn back to front (each face tinted), regions painted on
   them. Uses figure 23's helpers. */

/* faces round the origin, face i spanned by axis directions d[i] and
   d[i + 1] (scaled by size): flat bilinear squares (bend 0), or bicubic
   patches whose control grid is the flat one bent (shared edges share
   their control points, so the faces still meet along them) */
static void ssd_bend(double const* p, double k, qaws_vec3* q)
{
	q->x = (qaws_scalar)(p[0] + k * 0.35 * sin(2.2 * p[1] + 0.7 * p[2]));
	q->y = (qaws_scalar)(p[1] + k * 0.35 * sin(2.2 * p[2] + 0.7 * p[0]));
	q->z = (qaws_scalar)(p[2] + k * 0.45 * sin(2.0 * p[0] + 1.3 * p[1]) * cos(1.1 * p[2]));
}

static void ssd_fan(qaws_surface** f, double const (*d)[3], unsigned int n, double size, double bend)
{
	static double const O[3] = { 0, 0, 0 };
	unsigned int i, j, q, m;
	for (i = 0; i < n; i++)
	{
		double const* a = d[i];
		double const* b = d[(i + 1) % n];
		if (bend == 0)
		{
			double pa[3], pb[3], pc[3];
			for (m = 0; m < 3; m++) { pa[m] = size * a[m]; pb[m] = size * b[m]; pc[m] = pa[m] + pb[m]; }
			f[i] = shf_face(O, pa, pb, pc);
		}
		else
		{
			qaws_vec3 cp[16];
			qaws_surface_bezier_desc bd;
			for (j = 0; j < 4; j++)
				for (q = 0; q < 4; q++)
				{
					double p[3];
					for (m = 0; m < 3; m++) p[m] = size * (a[m] * j / 3.0 + b[m] * q / 3.0);
					ssd_bend(p, bend, &cp[4 * j + q]);
				}
			memset(&bd, 0, sizeof(bd));
			bd.u_degree = 3; bd.v_degree = 3; bd.control_points = cp; bd.u_point_count = 4; bd.v_point_count = 4;
			f[i] = NULL;
			qaws_surface_create_bezier(&bd, &f[i]);
		}
	}
}

/* ------------------------------------------------------------------ */
/*  A painter's scene: shaded quads and segments, drawn far to near     */
/* ------------------------------------------------------------------ */

typedef struct ssd_item
{
	double depth;
	double xy[8];
	unsigned int n;             /* 4: a quad, 2: a segment, 1: a dot */
	int rgb[3];
	double width;
} ssd_item;

typedef struct ssd_scene
{
	ssd_item* it;
	unsigned int n, cap;
	sgf_view const* view;
	/* optional, per face: its trimmed part (quads outside left out), the
	   (u, v) box to tessellate, the quads along each side (0: 34) */
	qaws_surface_region const* const* trims;
	double const (*box)[4];
	unsigned int res;
	int const (*tints)[3];      /* optional, per face: its colour (else a pale tint by its number) */
	int const* inner;           /* optional, per face: +1 (its normal points out) or -1 (in), its inner side darker */
	int (*keep)(unsigned int face, double u, double v);   /* optional: the quads drawn (by their middle) */
} ssd_scene;

static ssd_item* ssd_new(ssd_scene* sc)
{
	if (sc->n + 1 > sc->cap)
	{
		unsigned int cap = sc->cap ? 2 * sc->cap : 4096;
		ssd_item* q = (ssd_item*)realloc(sc->it, sizeof(ssd_item) * cap);
		if (!q)
			return NULL;
		sc->it = q;
		sc->cap = cap;
	}
	return &sc->it[sc->n++];
}

static int ssd_far_first(void const* a, void const* b)
{
	double da = ((ssd_item const*)a)->depth, db = ((ssd_item const*)b)->depth;
	return da < db ? 1 : (da > db ? -1 : 0);
}

static void ssd_draw(svg* s, ssd_scene* sc)
{
	unsigned int i;
	qsort(sc->it, sc->n, sizeof(ssd_item), ssd_far_first);
	for (i = 0; i < sc->n; i++)
	{
		ssd_item const* t = &sc->it[i];
		if (t->n == 4)
			fprintf(s->f, "<polygon points=\"%.1f,%.1f %.1f,%.1f %.1f,%.1f %.1f,%.1f\" fill=\"rgb(%d,%d,%d)\" stroke=\"rgb(%d,%d,%d)\" stroke-width=\"0.6\"/>\n",
				t->xy[0], t->xy[1], t->xy[2], t->xy[3], t->xy[4], t->xy[5], t->xy[6], t->xy[7],
				t->rgb[0], t->rgb[1], t->rgb[2], t->rgb[0], t->rgb[1], t->rgb[2]);
		else if (t->n == 2)
			fprintf(s->f, "<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"rgb(%d,%d,%d)\" stroke-width=\"%.2f\" stroke-linecap=\"round\"/>\n",
				t->xy[0], t->xy[1], t->xy[2], t->xy[3], t->rgb[0], t->rgb[1], t->rgb[2], t->width);
		else
			fprintf(s->f, "<circle cx=\"%.1f\" cy=\"%.1f\" r=\"%.1f\" fill=\"rgb(%d,%d,%d)\" stroke=\"#ffffff\" stroke-width=\"1.2\"/>\n",
				t->xy[0], t->xy[1], t->width, t->rgb[0], t->rgb[1], t->rgb[2]);
	}
	free(sc->it);
	memset(sc, 0, sizeof(*sc));
}

/* the surface at (u, v), brought into its parameter ranges by whole ranges
   when past them (a closed surface's loops may run past its seam) */
static void ssd_eval(qaws_surface const* f, double u, double v, unsigned int flags, qaws_surface_eval_result* r)
{
	qaws_range ru = qaws_surface_get_u_range(f), rv = qaws_surface_get_v_range(f);
	double lu = ru.max_value - ru.min_value, lv = rv.max_value - rv.min_value;
	if (lu > 0 && (u < ru.min_value - 1e-9 * lu || u > ru.max_value + 1e-9 * lu)) u -= lu * floor((u - ru.min_value) / lu);
	if (lv > 0 && (v < rv.min_value - 1e-9 * lv || v > rv.max_value + 1e-9 * lv)) v -= lv * floor((v - rv.min_value) / lv);
	qaws_surface_evaluate(f, (qaws_scalar)u, (qaws_scalar)v, flags, r);
}

/* a region painted on a face: its colour where the face is inside it */
typedef struct ssd_paint
{
	qaws_shell_region const* r;
	int rgb[3];
} ssd_paint;

/* the faces, tessellated and shaded on both sides (a light from upper
   left, towards the viewer); each quad tinted per face and painted by the
   regions holding its middle */
static void ssd_surfaces(ssd_scene* sc, qaws_surface** f, unsigned int nf, ssd_paint const* paint, unsigned int np)
{
	static int const tint[6][3] = { { 236, 229, 214 }, { 214, 226, 236 }, { 222, 236, 214 }, { 236, 216, 226 }, { 222, 220, 238 }, { 238, 230, 206 } };
	unsigned int const N = sc->res ? sc->res : 34;
	double L[3] = { -0.35, 0.55, -0.75 }, ll = sqrt(0.35 * 0.35 + 0.55 * 0.55 + 0.75 * 0.75);
	unsigned int k, i, j, m;
	sgf_view const* v = sc->view;
	L[0] /= ll; L[1] /= ll; L[2] /= ll;
	for (k = 0; k < nf; k++)
	{
		qaws_range ru = qaws_surface_get_u_range(f[k]), rv = qaws_surface_get_v_range(f[k]);
		if (sc->box)
		{
			ru.min_value = (qaws_scalar)sc->box[k][0]; ru.max_value = (qaws_scalar)sc->box[k][1];
			rv.min_value = (qaws_scalar)sc->box[k][2]; rv.max_value = (qaws_scalar)sc->box[k][3];
		}
		for (i = 0; i < N; i++)
			for (j = 0; j < N; j++)
			{
				ssd_item* t;
				if (sc->trims && sc->trims[k])
				{
					int in = 0;
					qaws_vec2 q0;
					q0.x = (qaws_scalar)(ru.min_value + (ru.max_value - ru.min_value) * (i + 0.5) / N);
					q0.y = (qaws_scalar)(rv.min_value + (rv.max_value - rv.min_value) * (j + 0.5) / N);
					qaws_surface_region_contains(sc->trims[k], q0, &in);
					if (in != 1)
						continue;
				}
				if (sc->keep && !sc->keep(k, ru.min_value + (ru.max_value - ru.min_value) * (i + 0.5) / N, rv.min_value + (rv.max_value - rv.min_value) * (j + 0.5) / N))
					continue;
				t = ssd_new(sc);
				double uu[4], vv[4], dsum = 0, cu, cv, sh, col[3];
				qaws_surface_eval_result r;
				if (!t)
					return;
				uu[0] = i; vv[0] = j; uu[1] = i + 1; vv[1] = j; uu[2] = i + 1; vv[2] = j + 1; uu[3] = i; vv[3] = j + 1;
				for (m = 0; m < 4; m++)
				{
					double u = ru.min_value + (ru.max_value - ru.min_value) * uu[m] / N, w = rv.min_value + (rv.max_value - rv.min_value) * vv[m] / N, dep;
					ssd_eval(f[k], u, w, QAWS_SURFACE_EVAL_POSITION, &r);
					sgf_project(v, r.position, &t->xy[2 * m], &t->xy[2 * m + 1], &dep);
					dsum += dep;
				}
				t->n = 4;
				t->depth = dsum / 4;
				cu = ru.min_value + (ru.max_value - ru.min_value) * (i + 0.5) / N;
				cv = rv.min_value + (rv.max_value - rv.min_value) * (j + 0.5) / N;
				ssd_eval(f[k], cu, cv, QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV, &r);
				{
					double nx = (double)r.du.y * r.dv.z - (double)r.du.z * r.dv.y, ny = (double)r.du.z * r.dv.x - (double)r.du.x * r.dv.z;
					double nz = (double)r.du.x * r.dv.y - (double)r.du.y * r.dv.x, nn = sqrt(nx * nx + ny * ny + nz * nz);
					double a = v->rot[0] * nx + v->rot[1] * ny + v->rot[2] * nz, b = v->rot[3] * nx + v->rot[4] * ny + v->rot[5] * nz;
					double c = v->rot[6] * nx + v->rot[7] * ny + v->rot[8] * nz;
					sh = nn > 0 ? 0.55 + 0.45 * fabs(a * L[0] + b * L[1] + c * L[2]) / nn : 1;
					/* (a face's inner side, seen through a hole, darker) */
					if (sc->inner && sc->inner[k] * c > 0) sh *= 0.5;
				}
				for (m = 0; m < 3; m++) col[m] = sc->tints ? sc->tints[k][m] : tint[k % 6][m];
				for (m = 0; m < np; m++)
				{
					int in = 0;
					qaws_vec2 q;
					qaws_surface_region const* g = paint[m].r ? qaws_shell_region_get_face(paint[m].r, k) : NULL;
					q.x = (qaws_scalar)cu; q.y = (qaws_scalar)cv;
					if (g) qaws_surface_region_contains(g, q, &in);
					if (in == 1)
					{
						unsigned int c;
						for (c = 0; c < 3; c++) col[c] = 0.3 * col[c] + 0.7 * paint[m].rgb[c];
					}
				}
				for (m = 0; m < 3; m++)
				{
					double x = col[m] * sh;
					t->rgb[m] = (int)(x > 255 ? 255 : x);
				}
			}
	}
}

/* a (u, v) polyline on a face as segments, lifted a hair towards the viewer
   so it shows over the face it lies on */
static void ssd_uv_line(ssd_scene* sc, qaws_surface const* f, double const* uv, unsigned int n, int const* rgb, double width)
{
	unsigned int i;
	double px = 0, py = 0, pd = 0;
	for (i = 0; i < n; i++)
	{
		qaws_surface_eval_result r;
		double x, y, d;
		ssd_eval(f, uv[2 * i], uv[2 * i + 1], QAWS_SURFACE_EVAL_POSITION, &r);
		sgf_project(sc->view, r.position, &x, &y, &d);
		if (i > 0)
		{
			ssd_item* t = ssd_new(sc);
			if (!t)
				return;
			t->n = 2;
			t->xy[0] = px; t->xy[1] = py; t->xy[2] = x; t->xy[3] = y;
			t->depth = 0.5 * (pd + d) - 0.03;
			t->rgb[0] = rgb[0]; t->rgb[1] = rgb[1]; t->rgb[2] = rgb[2];
			t->width = width;
		}
		px = x; py = y; pd = d;
	}
}

static void ssd_curve_line(ssd_scene* sc, qaws_surface const* f, qaws_curve const* cv, int const* rgb, double width)
{
	double uv[2 * 97];
	qaws_range rg = qaws_curve_get_parameter_range(cv);
	unsigned int k;
	for (k = 0; k <= 96; k++)
	{
		qaws_eval_result_2d e;
		qaws_curve_evaluate_2d(&e, cv, (qaws_scalar)(rg.min_value + (rg.max_value - rg.min_value) * k / 96.0), QAWS_EVAL_FLAG_POSITION);
		uv[2 * k] = e.position.x; uv[2 * k + 1] = e.position.y;
	}
	ssd_uv_line(sc, f, uv, 97, rgb, width);
}

/* the faces' sides, dark */
static void ssd_edges(ssd_scene* sc, qaws_surface** f, unsigned int nf)
{
	static int const ink[3] = { 87, 96, 106 };
	unsigned int k, s, q;
	for (k = 0; k < nf; k++)
	{
		qaws_range ru = qaws_surface_get_u_range(f[k]), rv = qaws_surface_get_v_range(f[k]);
		for (s = 0; s < 4; s++)
		{
			double uv[2 * 41];
			for (q = 0; q <= 40; q++)
			{
				double t = q / 40.0;
				uv[2 * q] = s == 0 || s == 2 ? ru.min_value + t * (ru.max_value - ru.min_value) : (s == 1 ? ru.max_value : ru.min_value);
				uv[2 * q + 1] = s == 1 || s == 3 ? rv.min_value + t * (rv.max_value - rv.min_value) : (s == 2 ? rv.max_value : rv.min_value);
			}
			ssd_uv_line(sc, f[k], uv, 41, ink, 1.1);
		}
	}
}

/* a shell region's boundary on every face */
/* a region's cut along a seam: the curve on a side of the domain (u or v at
   0 or 1 at its ends and middle) whose opposite side is the same place */
static int ssd_on_seam(qaws_surface const* f, qaws_curve const* c)
{
	qaws_range rg = qaws_curve_get_parameter_range(c);
	qaws_eval_result_2d e[3];
	qaws_surface_eval_result a, b;
	unsigned int ax, i;
	for (i = 0; i < 3; i++)
		qaws_curve_evaluate_2d(&e[i], c, rg.min_value + (rg.max_value - rg.min_value) * (qaws_scalar)(0.5 * i), QAWS_EVAL_FLAG_POSITION);
	for (ax = 0; ax < 2; ax++)
	{
		double w = ax ? e[1].position.y : e[1].position.x, q[2];
		int side = 1;
		for (i = 0; i < 3; i++)
			side &= fabs((ax ? e[i].position.y : e[i].position.x) - w) < 1e-6 && (fabs(w) < 1e-6 || fabs(w - 1) < 1e-6);
		if (!side)
			continue;
		q[0] = e[1].position.x; q[1] = e[1].position.y;
		q[ax] = 1 - q[ax];
		ssd_eval(f, e[1].position.x, e[1].position.y, QAWS_SURFACE_EVAL_POSITION, &a);
		ssd_eval(f, q[0], q[1], QAWS_SURFACE_EVAL_POSITION, &b);
		if (fabs((double)(a.position.x - b.position.x)) + fabs((double)(a.position.y - b.position.y)) + fabs((double)(a.position.z - b.position.z)) < 1e-5)
			return 1;
	}
	return 0;
}

static void ssd_outline(ssd_scene* sc, qaws_shell_region const* g, qaws_surface** f, unsigned int nf, int const* rgb)
{
	unsigned int k, p, c;
	if (!g)
		return;
	for (k = 0; k < nf; k++)
	{
		qaws_surface_region const* r = qaws_shell_region_get_face(g, k);
		if (!r)
			continue;
		for (p = 0; p < qaws_surface_region_get_path_count(r); p++)
		{
			qaws_path_2d path;
			qaws_surface_region_get_path(r, p, &path);
			for (c = 0; c < path.curve_count; c++)
				if (!ssd_on_seam(f[k], path.curves[c]))
					ssd_curve_line(sc, f[k], path.curves[c], rgb, 2.0);
		}
	}
}

/* a point of the shell as a dot (the vertex, a centre) */
static void ssd_dot(ssd_scene* sc, qaws_surface const* f, double u, double v, int const* rgb, double radius)
{
	qaws_surface_eval_result r;
	ssd_item* t = ssd_new(sc);
	double d;
	if (!t)
		return;
	ssd_eval(f, u, v, QAWS_SURFACE_EVAL_POSITION, &r);
	sgf_project(sc->view, r.position, &t->xy[0], &t->xy[1], &d);
	t->n = 1;
	t->depth = d - 0.1;
	t->rgb[0] = rgb[0]; t->rgb[1] = rgb[1]; t->rgb[2] = rgb[2];
	t->width = radius;
}

/* straightest rays from c over len, at angles about the direction th0:
   left of it red, right of it purple */
static void ssd_rays(ssd_scene* sc, qaws_shell_atlas const* a, qaws_surface** f, qaws_shell_point c, double th0, double len)
{
	static int const red[3] = { 207, 34, 46 }, purple[3] = { 130, 80, 223 };
	double const pi = 3.14159265358979323846;
	unsigned int i, w;
	for (i = 0; i < 24; i++)
	{
		qaws_shell_curve* g = NULL;
		qaws_vec2 d;
		double th = th0 + 0.5 * pi * ((i + 0.5) / 24.0 - 0.5);
		d.x = (qaws_scalar)cos(th); d.y = (qaws_scalar)sin(th);
		if (qaws_shell_curve_create_geodesic_ray(a, c, d, (qaws_scalar)len, 0, &g) == QAWS_STATUS_OK)
			for (w = 0; w < qaws_shell_curve_get_run_count(g); w++)
			{
				unsigned int fw;
				qaws_surface_curve const* run;
				qaws_shell_curve_get_run(g, w, &fw, &run);
				ssd_curve_line(sc, f[fw], qaws_surface_curve_get_uv(run), i < 12 ? purple : red, 1.3);
			}
		qaws_shell_curve_destroy(g);
	}
}

static qaws_shell_region* ssd_disk(qaws_shell_atlas const* a, qaws_shell_point c, double r, double* area)
{
	qaws_shell_region* g = NULL;
	qaws_scalar A = 0;
	qaws_shell_region_create_geodesic_disk(a, c, (qaws_scalar)r, 0, &g);
	if (g) qaws_shell_region_get_area(g, &A);
	if (area) *area = (double)A;
	return g;
}

/* a panel: its title and a line saying what to look at */
static void ssd_panel(svg* s, viewport const* v, char const* title, char const* note)
{
	svg_panel(s, v, title);
	if (note)
		fprintf(s->f, "<text x=\"%.1f\" y=\"%.1f\" font-size=\"12\" fill=\"#57606a\">%s</text>\n", v->x0 + 10, v->y0 + 38, note);
}

static void demo_shell_saddle(void)
{
	clock_t t0 = clock();
	svg s;
	double const pw = 400, ph = 400, pi = 3.14159265358979323846, rho0 = 0.3 * sqrt(2.0);
	static double const D5[5][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { -1, 0, 0 }, { 0, -1, 0 }, { 0, 0, -1 } };
	static double const D6[6][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { -1, 0, 0 }, { 0, -1, 0 }, { 0, 0, -1 } };
	static int const blue[3] = { 9, 105, 218 }, red[3] = { 207, 34, 46 }, black[3] = { 27, 31, 36 };
	static int const fblue[3] = { 84, 174, 255 }, fred[3] = { 255, 129, 130 }, fpurple[3] = { 170, 120, 250 };
	unsigned int k, i;
	char label[220];
	if (!svg_open(&s, "showcase/b2d24_shell_saddle.svg", (int)(48 + 4 * (pw + 12)), (int)(110 + 3 * (ph + 44)),
		"Geodesic disks through saddle vertices, flat and curved",
		"Faces round one vertex with more than 2 pi of angle (black dot): rays from a centre part there, and a disk through it fills in behind. Rows: flat, bent, Booleans on the bent one."))
		return;
	for (k = 0; k < 12; k++)
	{
		/* panels 0-3 flat, 4-7 bent, 8-11 Booleans on the bent monkey saddle */
		int bent = k >= 4, six = !(k == 0 || k == 4);
		unsigned int n = six ? 6 : 5;
		qaws_surface* f[6];
		qaws_shell_atlas* a = NULL;
		qaws_shell_point c = bent ? shf_pt(0, 0.2, 0.2) : shf_pt(0, 0.3, 0.3), cb = shf_pt(2, 0.22, 0.19);
		viewport v3;
		sgf_view view;
		ssd_scene sc;
		ssd_paint paint[2];
		unsigned int np = 0;
		qaws_shell_region *g = NULL, *h = NULL, *x = NULL;
		double A = 0, B = 0;
		qaws_scalar X = 0;
		memset(&sc, 0, sizeof(sc));
		ssd_fan(f, six ? D6 : D5, n, bent ? 1.6 : 1, bent ? 1 : 0);
		qaws_shell_atlas_create((qaws_surface const* const*)f, n, 0, &a);
		v3.x0 = 24 + (k % 4) * (pw + 12); v3.y0 = 104 + (k / 4) * (ph + 44); v3.w = pw; v3.h = ph;
		if (six)
			sgf_view_init(&view, -0.25 * pi, 0.615, bent ? 92 : 125, v3.x0 + pw / 2, v3.y0 + ph / 2 + 18);
		else
			sgf_view_init(&view, 0.25 * pi, 0.5, bent ? 95 : 135, v3.x0 + pw / 2, v3.y0 + ph / 2 + 5);
		sc.view = &view;
		if (k == 0 || k == 3)
		{
			/* growing disks: short of the vertex, through it, past it */
			static double const R5[3] = { 0.75, 0.55, 0.35 }, R6[3] = { 0.65, 0.5, 0.35 };
			double const* R = k == 0 ? R5 : R6;
			static int const F[3][3] = { { 182, 227, 255 }, { 128, 204, 255 }, { 54, 150, 240 } };
			qaws_shell_region* gr[3];
			ssd_paint pp[3];
			if (k == 0)
				ssd_panel(&s, &v3, "five squares round a vertex (5 pi/2)", "disks r 0.35, 0.55, 0.75: past the vertex, one reaches down the walls");
			else
				ssd_panel(&s, &v3, "monkey saddle: disks r 0.35, 0.5, 0.65", "short of the vertex a plain disk; past it, reaching every face");
			for (i = 0; i < 3; i++) { gr[i] = ssd_disk(a, c, R[i], NULL); pp[i].r = gr[i]; memcpy(pp[i].rgb, F[i], sizeof(F[i])); }
			ssd_surfaces(&sc, f, n, pp, 3);
			for (i = 0; i < 3; i++) { ssd_outline(&sc, gr[i], f, n, blue); qaws_shell_region_destroy(gr[i]); }
		}
		else if (k == 1 || k == 5)
		{
			ssd_panel(&s, &v3, bent ? "bent monkey saddle: rays part at the vertex" : "monkey saddle (six squares, 3 pi): rays part",
				"rays just left (red) and right (purple) of the vertex go different ways");
			ssd_surfaces(&sc, f, n, NULL, 0);
			ssd_rays(&sc, a, f, c, 1.25 * pi, bent ? 1.5 : 1.3);
		}
		else if (k == 2 || k == 4 || k == 6)
		{
			double r = 0.6;
			g = ssd_disk(a, c, r, &A);
			if (k == 2)
				sprintf(label, "monkey saddle: disk r %.2f, area %.6f", r, A);
			else
				sprintf(label, "bent %s: disk r 0.6, area %.6f", six ? "monkey saddle" : "five-face saddle", A);
			if (k == 2)
			{
				char note[160];
				sprintf(note, "closed form pi r^2 + (pi/2)(r - rho)^2 = %.6f", pi * r * r + 0.5 * pi * (r - rho0) * (r - rho0));
				ssd_panel(&s, &v3, label, note);
			}
			else
				ssd_panel(&s, &v3, label, "curved faces and curved edges: no closed form; f32 agrees with f64");
			paint[np].r = g; memcpy(paint[np].rgb, fblue, sizeof(fblue)); np++;
			ssd_surfaces(&sc, f, n, paint, np);
			ssd_outline(&sc, g, f, n, blue);
		}
		else
		{
			/* two disks through the vertex of the bent monkey saddle */
			static qaws_clip_type const T[4] = { QAWS_CLIP_UNION, QAWS_CLIP_INTERSECTION, QAWS_CLIP_DIFFERENCE, QAWS_CLIP_XOR };
			static char const* const TN[4] = { "A union B", "A intersect B", "A minus B", "A xor B" };
			g = ssd_disk(a, c, 0.6, &A);
			h = ssd_disk(a, cb, 0.55, &B);
			if (k == 7)
			{
				sprintf(label, "bent monkey saddle: disks A %.4f and B %.4f", A, B);
				ssd_panel(&s, &v3, label, "A (blue) and B (red) both reach round the vertex");
				paint[0].r = g; memcpy(paint[0].rgb, fblue, sizeof(fblue));
				paint[1].r = h; memcpy(paint[1].rgb, fred, sizeof(fred));
				np = 2;
			}
			else
			{
				if (g && h) qaws_shell_region_clip(g, h, T[k - 8], &x);
				if (x) qaws_shell_region_get_area(x, &X);
				sprintf(label, "%s: area %.6f", TN[k - 8], (double)X);
				ssd_panel(&s, &v3, label, "the result filled; A and B outlined");
				paint[0].r = x; memcpy(paint[0].rgb, fpurple, sizeof(fpurple));
				np = 1;
			}
			ssd_surfaces(&sc, f, n, paint, np);
			ssd_outline(&sc, g, f, n, blue);
			ssd_outline(&sc, h, f, n, red);
			ssd_dot(&sc, f[cb.face], cb.uv.x, cb.uv.y, red, 3.5);
		}
		ssd_edges(&sc, f, n);
		ssd_dot(&sc, f[0], 0, 0, black, 4.5);
		ssd_dot(&sc, f[0], c.uv.x, c.uv.y, k >= 8 || k == 7 ? blue : black, 3.5);
		ssd_draw(&s, &sc);
		qaws_shell_region_destroy(g); qaws_shell_region_destroy(h); qaws_shell_region_destroy(x);
		qaws_shell_atlas_destroy(a);
		for (i = 0; i < n; i++) qaws_surface_destroy(f[i]);
	}
	svg_close(&s);
	printf("figure 24: geodesic disks through saddle vertices, %.1f s\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
}
