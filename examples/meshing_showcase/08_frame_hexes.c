/*
 * Figure 8: hexahedra following a frame field (FRAME_FIELD, an integer-grid
 * map) against the grid's (GRID):
 *   - a box turned in space: on the grid, a staircase padded to its faces;
 *     along its frame field, a grid aligned with it, every cell a cube;
 *   - a bent bar and a twisted bar along their frame fields: the cells
 *     follow the bend and the twist.
 * Half of each mesh is shown, the cut faces coloured by their cell's
 * smallest scaled Jacobian (blue at 0.2 to red at 1).
 */

typedef void (*ms_map_fn)(double const* p, double* x);

static void ms_map_turned_box(double const* p, double* x)
{
	double a = p[0] * 2 - 1, b = (p[1] * 2 - 1) * 0.7, c = (p[2] * 2 - 1) * 0.5, y0 = a * cos(0.5) - b * sin(0.5), y1 = a * sin(0.5) + b * cos(0.5);
	x[0] = y0;
	x[1] = y1 * cos(0.3) - c * sin(0.3);
	x[2] = y1 * sin(0.3) + c * cos(0.3);
}

static void ms_map_bent_bar(double const* p, double* x)
{
	double r = 1 + p[0], t = MS_PI / 2 * p[1];
	x[0] = r * cos(t);
	x[1] = r * sin(t);
	x[2] = p[2];
}

static void ms_map_twisted_bar(double const* p, double* x)
{
	double a = (p[0] * 2 - 1) * 0.5, b = (p[1] * 2 - 1) * 0.5, t = MS_PI / 3 * p[2];
	x[0] = a * cos(t) - b * sin(t);
	x[1] = a * sin(t) + b * cos(t);
	x[2] = 3 * p[2];
}

/* the border of the lattice block n[0] x n[1] x n[2] as quads, mapped */
static qaws_mesh_3d* ms_lattice_surface(unsigned int const* n, ms_map_fn f)
{
	static unsigned char const F[6][4] = { { 0, 3, 2, 1 }, { 4, 5, 6, 7 }, { 0, 1, 5, 4 }, { 1, 2, 6, 5 }, { 2, 3, 7, 6 }, { 3, 0, 4, 7 } };
	static int const FN[6][3] = { { 0, 0, -1 }, { 0, 0, 1 }, { 0, -1, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { -1, 0, 0 } };
	static int const CO[8][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { 1, 0, 1 }, { 1, 1, 1 }, { 0, 1, 1 } };
	unsigned int nv = (n[0] + 1) * (n[1] + 1) * (n[2] + 1), *map = (unsigned int*)malloc(sizeof(unsigned int) * nv), *fv = NULL, nf = 0, cap = 0, i, j, k, q, used = 0;
	double* x = (double*)malloc(sizeof(double) * 3 * nv);
	unsigned int* fs;
	qaws_mesh_3d* m;
	for (i = 0; i < nv; i++)
		map[i] = QAWS_MESH_NONE;
	for (k = 0; k < n[2]; k++)
		for (j = 0; j < n[1]; j++)
			for (i = 0; i < n[0]; i++)
				for (q = 0; q < 6; q++)
				{
					long a = (long)i + FN[q][0], b = (long)j + FN[q][1], c = (long)k + FN[q][2];
					unsigned int e;
					if (a >= 0 && b >= 0 && c >= 0 && a < (long)n[0] && b < (long)n[1] && c < (long)n[2])
						continue;
					if (nf + 1 > cap)
					{
						cap = cap ? 2 * cap : 256;
						fv = (unsigned int*)realloc(fv, sizeof(unsigned int) * 4 * cap);
					}
					for (e = 0; e < 4; e++)
					{
						unsigned int ci = i + CO[F[q][e]][0], cj = j + CO[F[q][e]][1], ck = k + CO[F[q][e]][2];
						unsigned int v = ci + (n[0] + 1) * (cj + (n[1] + 1) * ck);
						if (map[v] == QAWS_MESH_NONE)
						{
							double p[3];
							p[0] = (double)ci / n[0];
							p[1] = (double)cj / n[1];
							p[2] = (double)ck / n[2];
							f(p, x + 3 * used);
							map[v] = used++;
						}
						fv[4 * nf + e] = map[v];
					}
					nf++;
				}
	fs = (unsigned int*)malloc(sizeof(unsigned int) * nf);
	for (i = 0; i < nf; i++)
		fs[i] = 4;
	m = ms_make(x, used, fv, fs, nf);
	free(map); free(x); free(fv); free(fs);
	return m;
}

static int ms_keep_all(double const* p)
{
	(void)p;
	return 1;
}

/* the turned box's own half: its local y (the map turned back) positive */
static int ms_keep_box_half(double const* p)
{
	double y1 = p[1] * cos(0.3) + p[2] * sin(0.3);
	return -p[0] * sin(0.5) + y1 * cos(0.5) > 0;
}

static qaws_volume_mesh* ms_hexes_by(qaws_mesh_3d const* b, qaws_hexahedralize_method method, double h, char* sub, char const* what)
{
	qaws_hexahedralize_desc d;
	qaws_volume_mesh* v = NULL;
	qaws_volume_quality q;
	qaws_status st;
	memset(&d, 0, sizeof d);
	d.method = method;
	d.cell_size = (qaws_scalar)h;
	st = qaws_volume_mesh_hexahedralize(&v, NULL, b, &d);
	if (st != QAWS_STATUS_OK)
	{
		printf("  %s: hexahedralize failed: %s\n", what, qaws_status_to_string(st));
		return NULL;
	}
	qaws_volume_mesh_measure_quality(&q, v);
	sprintf(sub, "%u hexahedra, scaled Jacobian from %.2f", q.hexahedra, (double)q.min_scaled_jacobian);
	printf("  %s: %u hexahedra, scaled Jacobian %.3f, volume %.5f, edges %.3f..%.3f, %s\n", what, q.hexahedra, (double)q.min_scaled_jacobian,
		(double)q.volume, (double)q.min_edge, (double)q.max_edge, qaws_status_to_string(qaws_volume_mesh_validate(v)));
	return v;
}

static void demo_frame_hexes(void)
{
	static double const wire[3] = { 40, 45, 60 }, tan[3] = { 232, 196, 120 }, green[3] = { 170, 205, 160 }, tint[3] = { 178, 192, 212 };
	static unsigned int const nbox[3] = { 16, 12, 8 }, nbent[3] = { 8, 20, 8 }, ntwist[3] = { 8, 8, 24 };
	er_canvas* c = er_create(2750, 700, 3, 252, 252, 250);
	qaws_mesh_3d *box = ms_lattice_surface(nbox, ms_map_turned_box), *bent = ms_lattice_surface(nbent, ms_map_bent_bar);
	qaws_mesh_3d* twist = ms_lattice_surface(ntwist, ms_map_twisted_bar);
	qaws_volume_mesh *g1 = NULL, *f1 = NULL, *f2 = NULL, *f3 = NULL;
	char sub[200];
	if (!c || !box || !bent || !twist)
		return;
	g1 = ms_hexes_by(box, QAWS_HEXES_GRID, 0.15, sub, "turned box, grid");
	if (g1)
	{
		ms_vpanel(c, 0, 0, 687, 700, g1, 0.62, 0.55, "a turned box, on a grid", sub);
		ms_hex_cutaway(c, g1, ms_keep_box_half, tan, wire);
	}
	f1 = ms_hexes_by(box, QAWS_HEXES_FRAME_FIELD, 0.15, sub, "turned box, frame field");
	if (f1)
	{
		ms_vpanel(c, 687, 0, 687, 700, f1, 0.62, 0.55, "along its frame field", sub);
		ms_hex_cutaway(c, f1, ms_keep_box_half, tan, wire);
	}
	f2 = ms_hexes_by(bent, QAWS_HEXES_FRAME_FIELD, 0.15, sub, "bent bar");
	if (f2)
	{
		ms_vpanel(c, 1374, 0, 687, 700, f2, 0.9, 0.6, "a bent bar", sub);
		ms_hex_cutaway(c, f2, ms_keep_all, green, wire);
	}
	f3 = ms_hexes_by(twist, QAWS_HEXES_FRAME_FIELD, 0.15, sub, "twisted bar");
	if (f3)
	{
		ms_vpanel(c, 2061, 0, 687, 700, f3, 0.5, 1.1, "a twisted bar", sub);
		ms_hex_cutaway(c, f3, ms_keep_all, tint, wire);
	}
	er_save(c, "showcase/mesh8_frame_hexes.png");
	printf("mesh8_frame_hexes.png\n");
	qaws_volume_mesh_destroy(g1);
	qaws_volume_mesh_destroy(f1);
	qaws_volume_mesh_destroy(f2);
	qaws_volume_mesh_destroy(f3);
	qaws_mesh_3d_destroy(box);
	qaws_mesh_3d_destroy(bent);
	qaws_mesh_3d_destroy(twist);
	er_destroy(c);
}
