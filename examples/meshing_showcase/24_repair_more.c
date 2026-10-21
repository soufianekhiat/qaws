/*
 * Figure 24: repair, completed (qaws_mesh_repair).
 *   - a unit box of 20 x 20 grids with a hole cut round a corner;
 *   - filled along its creases: the three crease lines entering the hole
 *     carried to the corner they meet at, each sector filled and faired
 *     on its own (the fills tinted): the corner back, sharp;
 *   - two cubes crossing, with QAWS_MESH_REPAIR_FIND_INTERSECTIONS: the
 *     faces in a crossing pair (exact predicates) in red.
 */

static double const ms_rb_o[6][3] = { { 0, 0, 0 }, { 0, 0, 1 }, { 0, 0, 0 }, { 0, 1, 0 }, { 0, 0, 0 }, { 1, 0, 0 } };
static double const ms_rb_u[6][3] = { { 0, 1, 0 }, { 1, 0, 0 }, { 1, 0, 0 }, { 0, 0, 1 }, { 0, 0, 1 }, { 0, 1, 0 } };
static double const ms_rb_v[6][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };

/* a box [lo, lo + s]^3 of n x n grids per face, each face its own points (a soup to weld); appended */
static void ms_rb_box(double const* lo, double s, unsigned int n, double* x, unsigned int* nv, unsigned int* t, unsigned int* nt)
{
	unsigned int f, i, j, k;
	for (f = 0; f < 6; f++)
	{
		unsigned int b = *nv;
		for (i = 0; i <= n; i++)
			for (j = 0; j <= n; j++)
			{
				for (k = 0; k < 3; k++)
					x[3 * *nv + k] = lo[k] + s * (ms_rb_o[f][k] + (double)i / n * ms_rb_u[f][k] + (double)j / n * ms_rb_v[f][k]);
				(*nv)++;
			}
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
			{
				unsigned int p00 = b + i * (n + 1) + j, p10 = p00 + n + 1;
				t[3 * *nt] = p00; t[3 * *nt + 1] = p10; t[3 * *nt + 2] = p10 + 1;
				t[3 * *nt + 3] = p00; t[3 * *nt + 4] = p10 + 1; t[3 * *nt + 5] = p00 + 1;
				*nt += 2;
			}
	}
}

static qaws_mesh_3d* ms_rb_repair(double const* x, unsigned int nv, unsigned int const* t, unsigned int nt, unsigned int flags,
	qaws_mesh_repair_report* r)
{
	qaws_scalar* p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	qaws_mesh_repair_desc d;
	qaws_mesh_3d* m = NULL;
	unsigned int i;
	for (i = 0; i < 3 * nv; i++)
		p[i] = (qaws_scalar)x[i];
	memset(&d, 0, sizeof d);
	d.flags = flags;
	qaws_mesh_repair(&m, r, p, t, NULL, nv, nt, &d);
	free(p);
	return m;
}

/* faces coloured by a face attribute: base where it is (or is not) the key, accent otherwise */
static void ms_rb_by_attr(er_canvas* c, qaws_mesh_3d const* m, char const* name, double key, int equal, double const* base, double const* accent)
{
	static double const wire[3] = { 50, 55, 70 };
	unsigned int nf = qaws_mesh_3d_get_face_count(m), a = QAWS_MESH_NONE, f;
	qaws_scalar* v = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (nf + 1));
	double* rgb = (double*)malloc(sizeof(double) * 3 * (nf + 1));
	qaws_mesh_3d_find_attribute(&a, m, name);
	if (a != QAWS_MESH_NONE)
		qaws_mesh_3d_get_attribute(v, m, a);
	for (f = 0; f < nf; f++)
	{
		int hit = a != QAWS_MESH_NONE && ((v[f] == key) == equal);
		memcpy(rgb + 3 * f, hit ? accent : base, sizeof(double) * 3);
	}
	ms_draw_faces(c, m, rgb);
	ms_draw(c, m, NULL, NULL, 0, wire, 0.6);
	free(v);
	free(rgb);
}

static void demo_repair_more(void)
{
	static double const tint[3] = { 200, 208, 222 }, fill[3] = { 244, 178, 92 }, red[3] = { 214, 52, 48 }, wire[3] = { 50, 55, 70 };
	static double const lo0[3] = { 0, 0, 0 }, lo1[3] = { 0.45, 0.4, 0.35 };
	er_canvas* c = er_create(2200, 700, 3, 252, 252, 250);
	double* x = (double*)malloc(sizeof(double) * 3 * 12 * 441);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * 12 * 800);
	unsigned int* kt = (unsigned int*)malloc(sizeof(unsigned int) * 3 * 12 * 800);
	unsigned int nv = 0, nt = 0, kn = 0, f, i, k;
	qaws_mesh_repair_report r;
	qaws_mesh_3d *holed, *filled, *crossed;
	char sub[200];
	if (!c || !x || !t || !kt)
		return;
	ms_rb_box(lo0, 1, 20, x, &nv, t, &nt);
	for (f = 0; f < nt; f++)
	{
		double g[3] = { 0, 0, 0 }, dd = 0;
		for (k = 0; k < 3; k++)
			for (i = 0; i < 3; i++)
				g[i] += x[3 * t[3 * f + k] + i] / 3;
		for (i = 0; i < 3; i++)
			dd += (g[i] - 1) * (g[i] - 1);
		if (dd > 0.35 * 0.35)
		{
			memcpy(kt + 3 * kn, t + 3 * f, sizeof(unsigned int) * 3);
			kn++;
		}
	}
	/* the holed box as it is (welded only: no fills) */
	{
		qaws_scalar* p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
		qaws_mesh_repair_desc d;
		for (i = 0; i < 3 * nv; i++)
			p[i] = (qaws_scalar)x[i];
		memset(&d, 0, sizeof d);
		d.hole_perimeter = -1;
		holed = NULL;
		qaws_mesh_repair(&holed, &r, p, kt, NULL, nv, kn, &d);
		free(p);
	}
	filled = ms_rb_repair(x, nv, kt, kn, 0, &r);
	if (holed)
	{
		sprintf(sub, "%u triangles, a hole round the corner", qaws_mesh_3d_get_face_count(holed));
		ms_panel(c, 0, 0, 733, 700, holed, 3.9, 0.55, "a box with a hole", sub);
		ms_draw(c, holed, tint, NULL, 0, wire, 0.6);
	}
	if (filled)
	{
		sprintf(sub, "%u fill faces, three creases carried to the corner", r.fill_faces);
		ms_panel(c, 733, 0, 733, 700, filled, 3.9, 0.55, "filled along its creases", sub);
		ms_rb_by_attr(c, filled, "source_face", -1, 1, tint, fill);
		printf("  filled: %s\n", sub);
	}
	/* two cubes crossing */
	nv = nt = 0;
	ms_rb_box(lo0, 1, 6, x, &nv, t, &nt);
	ms_rb_box(lo1, 1, 6, x, &nv, t, &nt);
	crossed = ms_rb_repair(x, nv, t, nt, QAWS_MESH_REPAIR_FIND_INTERSECTIONS, &r);
	if (crossed)
	{
		sprintf(sub, "%u crossing pairs, %u faces (exact predicates)", r.self_intersections, r.intersecting_faces);
		ms_panel(c, 1466, 0, 733, 700, crossed, 0.75, 0.55, "two cubes crossing", sub);
		ms_rb_by_attr(c, crossed, "intersecting", 0, 0, tint, red);
		printf("  crossing: %s\n", sub);
	}
	er_save(c, "showcase/mesh24_repair_more.png");
	printf("mesh24_repair_more.png\n");
	qaws_mesh_3d_destroy(holed);
	qaws_mesh_3d_destroy(filled);
	qaws_mesh_3d_destroy(crossed);
	free(x);
	free(t);
	free(kt);
	er_destroy(c);
}
