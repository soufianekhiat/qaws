/*
 * Figure 44: hexahedra of a polycube map (QAWS_HEXES_POLYCUBE): the inside
 * deformed until every border facet faces an axis, the integer-grid map's
 * lattice read off with the polycube's faces as its patches. Half of each
 * mesh is shown, the cut faces coloured by their cell's smallest scaled
 * Jacobian (blue at 0.2 to red at 1).
 *   - a sphere: its polycube a cube (FRAME_FIELD refuses it: its field has
 *     singular edges);
 *   - a box rounded by Catmull-Clark (its edges semi-sharp): a box again;
 *   - figure 40's ring smoothed by Catmull-Clark: a polycube of genus one;
 *   - figure 8's bent bar: labelled by nearest axes, its polycube turns a
 *     corner at 45 degrees.
 */

static int ms_pc_keep(double const* p)
{
	return p[1] > 0;
}

/* a panel of a polycube's hexahedra at cell size h, half shown when keep says */
static void ms_pc_panel(er_canvas* c, int px, qaws_mesh_3d const* s, double h, double yaw, double pitch, char const* title, int (*keep)(double const*),
	double const* rgb, double const* wire)
{
	char sub[200];
	qaws_volume_mesh* v = ms_hexes_by(s, QAWS_HEXES_POLYCUBE, h, sub, title);
	qaws_volume_quality q;
	if (!v)
		return;
	qaws_volume_mesh_measure_quality(&q, v);
	sprintf(sub, "%u cells, Jacobian >= %.2f", q.hexahedra, (double)q.min_scaled_jacobian);
	ms_vpanel(c, px, 0, 450, 700, v, yaw, pitch, title, sub);
	ms_hex_cutaway(c, v, keep, rgb, wire);
	qaws_volume_mesh_destroy(v);
}

static void demo_polycube(void)
{
	static double const wire[3] = { 40, 45, 60 }, tan[3] = { 232, 196, 120 }, green[3] = { 170, 205, 160 }, tint[3] = { 178, 192, 212 },
		rose[3] = { 225, 180, 180 };
	static unsigned int const nbent[3] = { 8, 20, 8 }, nbox[3] = { 2, 1, 1 };
	er_canvas* c = er_create(1800, 700, 3, 252, 252, 250);
	qaws_mesh_3d *s = ms_rm_sphere(32, 16), *cube = ms_lattice_surface(nbox, ms_map_turned_box), *cage = ms_sd_cage(0), *ring = NULL, *rounded = NULL;
	qaws_mesh_3d* bent = ms_lattice_surface(nbent, ms_map_bent_bar);
	qaws_mesh_subdivide_desc sd;
	qaws_scalar sh[64];
	unsigned int i;
	if (!c || !s || !cube || !cage || !bent)
		return;
	memset(&sd, 0, sizeof sd);
	sd.scheme = QAWS_MESH_SUBDIVIDE_CATMULL_CLARK;
	sd.levels = 2;
	qaws_mesh_3d_subdivide(&ring, cage, &sd);
	/* the turned box's 12 edges semi-sharp: rounded, not melted */
	qaws_mesh_3d_detect_features(cube, (qaws_scalar)0.5, 0, 0);
	for (i = 0; i < 64; i++) sh[i] = (qaws_scalar)1.5;
	sd.levels = 3;
	sd.sharpness = qaws_mesh_3d_get_feature_edge_count(cube) <= 64 ? sh : NULL;
	qaws_mesh_3d_subdivide(&rounded, cube, &sd);
	ms_pc_panel(c, 0, s, 0.2, 0.6, 0.5, "a sphere: a cube", ms_pc_keep, tint, wire);
	if (rounded)
		ms_pc_panel(c, 450, rounded, 0.15, 0.62, 0.55, "a rounded box", ms_keep_box_half, green, wire);
	if (ring)
		ms_pc_panel(c, 900, ring, 0.25, 2.6, 0.85, "a ring: genus one", ms_keep_all, rose, wire);
	ms_pc_panel(c, 1350, bent, 0.15, 0.5, 0.75, "a bent bar: a corner", ms_keep_all, tan, wire);
	qaws_mesh_3d_destroy(s);
	qaws_mesh_3d_destroy(cube);
	qaws_mesh_3d_destroy(rounded);
	qaws_mesh_3d_destroy(cage);
	qaws_mesh_3d_destroy(ring);
	qaws_mesh_3d_destroy(bent);
	er_save(c, "showcase/mesh44_polycube.png");
	printf("mesh44_polycube.png\n");
	er_destroy(c);
}
