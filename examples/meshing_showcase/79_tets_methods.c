/*
 * Figure 79: the tetrahedralization methods beside DELAUNAY, cut open (the
 * cut's faces coloured by their cell's smallest dihedral angle, blue below
 * 10 degrees to red at 70).
 *   - the NURBS torus's surface filled by DELAUNAY (its border the input's
 *     triangles);
 *   - the same by STUFFING (isosurface stuffing: a body-centred cubic
 *     lattice cut by the surface, vertices near it warped onto it): the
 *     lattice's regular cells inside, the border the surface's cuts;
 *   - a blob's surface with every 12th triangle taken out, by STUFFING: the
 *     winding number tells the inside of the soup.
 */

static qaws_volume_mesh* ms_tm_fill(qaws_mesh_3d const* s, qaws_tetrahedralize_method method, double L, char* sub, char const* name)
{
	qaws_tetrahedralize_desc d;
	qaws_volume_mesh* v = NULL;
	qaws_volume_quality q;
	qaws_status st;
	memset(&d, 0, sizeof d);
	d.method = method;
	d.target_length = (qaws_scalar)L;
	st = qaws_volume_mesh_tetrahedralize(&v, NULL, s, &d);
	if (st != QAWS_STATUS_OK || !v)
	{
		printf("  %s: %s\n", name, qaws_status_to_string(st));
		return NULL;
	}
	qaws_volume_mesh_measure_quality(&q, v);
	sprintf(sub, "%u tetrahedra, dihedral %.1f..%.1f", q.tetrahedra, (double)q.min_dihedral, (double)q.max_dihedral);
	printf("  %s: %u tetrahedra, %u vertices, dihedral %.2f..%.2f, radius-edge %.2f, volume %.4f\n", name, q.tetrahedra,
		qaws_volume_mesh_get_vertex_count(v), (double)q.min_dihedral, (double)q.max_dihedral, (double)q.max_radius_edge, (double)q.volume);
	return v;
}

static void demo_tets_methods(void)
{
	static double const wire[3] = { 40, 45, 60 }, green[3] = { 170, 205, 160 }, tint[3] = { 178, 192, 212 };
	er_canvas* c = er_create(2200, 700, 3, 252, 252, 250);
	qaws_surface* torus = ms_nurbs_torus(2, 0.7);
	qaws_mesh_3d *ts = NULL, *blob = ms_blob(2), *holed = NULL;
	qaws_volume_mesh* v;
	char sub[200];
	if (!c || !torus || !blob)
		return;
	ts = ms_generate(torus, NULL, QAWS_ELEMENTS_TRIANGLES, 0.25, 0, NULL);
	if (ts && (v = ms_tm_fill(ts, QAWS_TETS_DELAUNAY, 0.25, sub, "torus by DELAUNAY")) != NULL)
	{
		ms_vpanel(c, 0, 0, 733, 700, v, 0.5, 0.75, "DELAUNAY", sub);
		ms_cutaway(c, v, ms_keep_half, green, wire);
		qaws_volume_mesh_destroy(v);
	}
	if (ts && (v = ms_tm_fill(ts, QAWS_TETS_STUFFING, 0.25, sub, "torus by STUFFING")) != NULL)
	{
		ms_vpanel(c, 733, 0, 733, 700, v, 0.5, 0.75, "STUFFING", sub);
		ms_cutaway(c, v, ms_keep_half, green, wire);
		qaws_volume_mesh_destroy(v);
	}
	{
		/* (every 12th triangle out) */
		unsigned int nf = qaws_mesh_3d_get_face_count(blob), nv = qaws_mesh_3d_get_vertex_count(blob), nc = qaws_mesh_3d_get_corner_count(blob), f, k, o = 0, o2 = 0, kept = 0;
		unsigned int *fv = (unsigned int*)malloc(sizeof(unsigned int) * (nc + 1)), *fs = (unsigned int*)malloc(sizeof(unsigned int) * (nf + 1));
		unsigned int *fv2 = (unsigned int*)malloc(sizeof(unsigned int) * (nc + 1)), *fs2 = (unsigned int*)malloc(sizeof(unsigned int) * (nf + 1));
		qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (nv + 1));
		qaws_mesh_3d_desc md;
		qaws_mesh_3d_get_faces(fv, fs, blob);
		qaws_mesh_3d_get_positions(x, blob);
		for (f = 0; f < nf; f++)
		{
			if (f % 12 != 5)
			{
				for (k = 0; k < fs[f]; k++)
					fv2[o2 + k] = fv[o + k];
				o2 += fs[f];
				fs2[kept++] = fs[f];
			}
			o += fs[f];
		}
		memset(&md, 0, sizeof md);
		md.positions = x;
		md.vertex_count = nv;
		md.face_vertices = fv2;
		md.face_sizes = fs2;
		md.face_count = kept;
		qaws_mesh_3d_create(&holed, &md);
		free(fv); free(fs); free(fv2); free(fs2); free(x);
	}
	if (holed && (v = ms_tm_fill(holed, QAWS_TETS_STUFFING, 0.12, sub, "holed blob by STUFFING")) != NULL)
	{
		ms_vpanel(c, 1466, 0, 733, 700, v, 0.4, 0.6, "a soup with holes, STUFFING", sub);
		ms_cutaway(c, v, ms_keep_half, tint, wire);
		qaws_volume_mesh_destroy(v);
	}
	er_save(c, "showcase/mesh79_tets_methods.png");
	printf("mesh79_tets_methods.png\n");
	qaws_mesh_3d_destroy(ts);
	qaws_mesh_3d_destroy(blob);
	qaws_mesh_3d_destroy(holed);
	qaws_surface_destroy(torus);
	er_destroy(c);
}
