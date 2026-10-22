/*
 * Figure 57: polycube labels by a graph cut: a 2 x 1 x 1 box with a bump on
 * its top (a Gaussian 0.2 wide), as polycube hexahedra (half shown, the cut
 * faces coloured by their cell's smallest scaled Jacobian). The facets are
 * labelled by an alpha-expansion graph cut (their normals against the axes,
 * the labels' borders short): the bump's steep sides stay in the top's
 * patch, where nearest axes cut them off as small patches (scaled Jacobians
 * 0.36 and 0.29 for the bumps of 0.3 and 0.4).
 *   - the box with a bump of 0.3, its surface;
 *   - its hexahedra; then a bump of 0.4 and one of 0.45.
 */

static double g_ms_bump = 0.3;

static void ms_map_bumped_box(double const* p, double* x)
{
	x[0] = 2 * p[0];
	x[1] = p[1];
	x[2] = p[2];
	if (p[2] > 0.999)
		x[2] += g_ms_bump * exp(-((x[0] - 0.7) * (x[0] - 0.7) + (x[1] - 0.5) * (x[1] - 0.5)) / 0.04);
}

static void demo_polycube_cut(void)
{
	static double const wire[3] = { 40, 45, 60 }, tint[3] = { 178, 192, 212 }, warm[3] = { 225, 190, 150 };
	static unsigned int const n[3] = { 24, 12, 6 };
	static double const bumps[3] = { 0.3, 0.4, 0.45 };
	er_canvas* c = er_create(1800, 700, 3, 252, 252, 250);
	unsigned int k;
	char title[100];
	if (!c)
		return;
	for (k = 0; k < 3; k++)
	{
		qaws_mesh_3d* s;
		g_ms_bump = bumps[k];
		s = ms_lattice_surface(n, ms_map_bumped_box);
		if (!s)
			continue;
		if (k == 0)
		{
			ms_panel(c, 0, 0, 450, 700, s, 0.6, 0.45, "a box with a bump", "0.3 high, 0.2 wide");
			ms_draw(c, s, warm, NULL, 0, wire, 0.4);
		}
		sprintf(title, "a bump of %.2f", bumps[k]);
		ms_pc_panel(c, 450 + 450 * (int)k, s, 0.25, 0.6, 0.45, title, ms_keep_box_half, tint, wire);
		qaws_mesh_3d_destroy(s);
	}
	er_save(c, "showcase/mesh57_polycube_cut.png");
	printf("mesh57_polycube_cut.png\n");
	er_destroy(c);
}
