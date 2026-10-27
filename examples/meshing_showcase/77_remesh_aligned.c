/*
 * Figure 77: remeshing by FIELD_ALIGNED (Instant Meshes' triangles: a
 * 6-symmetric field and its hexagonal lattice), the vertices coloured by
 * valence (5 blue, 6 white, 7 red).
 *   - the blob by LOCAL and by FIELD_ALIGNED at the same target: the
 *     lattice's rows run along the field, the defects only at its
 *     singularities;
 *   - the cube with its creases: rows along the creases, corners kept;
 *   - a torus (no singularity needed: all valence 6 where the field
 *     closes).
 */

static void demo_remesh_aligned(void)
{
	er_canvas* c = er_create(2200, 640, 3, 252, 252, 250);
	qaws_mesh_3d *blob = ms_blob(1), *cube = ms_gx_cube(2), *torus = ms_torus(30, 14, 1, 0.42), *o;
	qaws_mesh_reference *rb = NULL, *rc = NULL, *rt = NULL;
	qaws_mesh_report rep;
	double sec;
	if (!c || !blob || !cube || !torus)
		return;
	qaws_mesh_reference_create(&rb, blob, NULL);
	if ((o = ms_remesh_by(blob, rb, QAWS_REMESH_LOCAL, 0.11, &rep, &sec)) != NULL)
		ms_valence_panel(c, 0, o, &rep, sec, "the blob by LOCAL");
	qaws_mesh_3d_destroy(o);
	if ((o = ms_remesh_by(blob, rb, QAWS_REMESH_FIELD_ALIGNED, 0.11, &rep, &sec)) != NULL)
		ms_valence_panel(c, 550, o, &rep, sec, "the blob by FIELD_ALIGNED");
	qaws_mesh_3d_destroy(o);
	qaws_mesh_3d_detect_features(cube, (qaws_scalar)(MS_PI / 6), 0, 0);
	qaws_mesh_reference_create(&rc, cube, NULL);
	if ((o = rc ? ms_remesh_by(cube, rc, QAWS_REMESH_FIELD_ALIGNED, 0.25, &rep, &sec) : NULL) != NULL)
		ms_valence_panel(c, 1100, o, &rep, sec, "the cube with creases");
	qaws_mesh_3d_destroy(o);
	qaws_mesh_reference_create(&rt, torus, NULL);
	if ((o = rt ? ms_remesh_by(torus, rt, QAWS_REMESH_FIELD_ALIGNED, 0.12, &rep, &sec) : NULL) != NULL)
		ms_valence_panel(c, 1650, o, &rep, sec, "a torus");
	qaws_mesh_3d_destroy(o);
	er_save(c, "showcase/mesh77_remesh_aligned.png");
	printf("mesh77_remesh_aligned.png\n");
	qaws_mesh_reference_destroy(rb);
	qaws_mesh_reference_destroy(rc);
	qaws_mesh_reference_destroy(rt);
	qaws_mesh_3d_destroy(blob);
	qaws_mesh_3d_destroy(cube);
	qaws_mesh_3d_destroy(torus);
	er_destroy(c);
}
