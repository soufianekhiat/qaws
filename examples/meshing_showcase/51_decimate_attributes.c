/*
 * Figure 51: attributes in the decimation's quadric (attribute_weight,
 * Hoppe 1999).
 *   - a torus painted in two colours along slanted bands, a colour per
 *     vertex;
 *   - decimated to a twelfth by the shape alone: the bands' edges smeared;
 *   - the same with the colour in the quadric: the faces gather on the
 *     bands' edges and they stay sharp.
 * Printed: the face counts and each result's colour error against the
 * original vertices' (root mean square, through the nearest vertex left).
 */

static void ms_da_paint(qaws_mesh_3d* m)
{
	static double const band[2][3] = { { 0.95, 0.6, 0.2 }, { 0.2, 0.35, 0.8 } };
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), i, k, id;
	qaws_scalar* col = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	qaws_attribute_desc ad;
	for (i = 0; i < nv; i++)
	{
		qaws_vec3 p = qaws_mesh_3d_get_vertex(m, i);
		double u = atan2((double)p.y, (double)p.x), v = atan2((double)p.z, sqrt((double)p.x * p.x + (double)p.y * p.y) - 1);
		int b = sin(5 * u + 2 * v) > 0;
		for (k = 0; k < 3; k++) col[3 * i + k] = (qaws_scalar)band[b][k];
	}
	memset(&ad, 0, sizeof ad);
	ad.name = "colour";
	ad.element = QAWS_MESH_VERTEX;
	ad.components = 3;
	ad.rule = QAWS_ATTRIBUTE_LINEAR;
	qaws_mesh_3d_add_attribute(&id, m, col, &ad);
	free(col);
}

/* the colours per vertex, 0 to 255 */
static double* ms_da_colours(qaws_mesh_3d const* m)
{
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), i, id;
	qaws_scalar* col = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	double* rgb = (double*)malloc(sizeof(double) * 3 * nv);
	qaws_mesh_3d_find_attribute(&id, m, "colour");
	qaws_mesh_3d_get_attribute(col, m, id);
	for (i = 0; i < 3 * nv; i++) rgb[i] = 255 * fmin(1, fmax(0, (double)col[i]));
	free(col);
	return rgb;
}

static double const ms_da_wire[3] = { 40, 40, 50 };

static void demo_decimate_attributes(void)
{
	enum { NU = 120, NV = 48 };
	er_canvas* c = er_create(2250, 720, 3, 252, 252, 250);
	double* x = (double*)malloc(sizeof(double) * 3 * NU * NV);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 6 * NU * NV), r;
	qaws_mesh_3d* m0;
	double* rgb0;
	char sub[200];
	if (!c)
		return;
	ms_torus_grid(NU, NV, x, t);
	m0 = ms_make(x, NU * NV, t, NULL, 2 * NU * NV);
	ms_da_paint(m0);
	rgb0 = ms_da_colours(m0);
	sprintf(sub, "%u faces, a colour per vertex", qaws_mesh_3d_get_face_count(m0));
	ms_panel(c, 0, 0, 750, 720, m0, 0.5, 0.9, "a torus painted in bands", sub);
	ms_draw(c, m0, NULL, rgb0, 1, NULL, 0);
	for (r = 0; r < 2; r++)
	{
		qaws_mesh_3d* m = NULL;
		qaws_decimate_desc d;
		qaws_decimate_report rep;
		double* rgb;
		double err = 0;
		unsigned int nv0 = qaws_mesh_3d_get_vertex_count(m0), i, k;
		qaws_mesh_3d_clone(&m, m0);
		memset(&d, 0, sizeof d);
		d.target_faces = 2 * NU * NV / 12;
		d.attribute_weight = (qaws_scalar)(r ? 0.03 : 0);
		qaws_mesh_3d_decimate(&rep, m, &d);
		rgb = ms_da_colours(m);
		for (i = 0; i < nv0; i++)
		{
			/* (the original vertex against the nearest vertex left) */
			qaws_vec3 p = qaws_mesh_3d_get_vertex(m0, i);
			unsigned int nv = qaws_mesh_3d_get_vertex_count(m), j, best = 0;
			double bd = HUGE_VAL;
			for (j = 0; j < nv; j++)
			{
				qaws_vec3 q = qaws_mesh_3d_get_vertex(m, j);
				double dd = ((double)p.x - q.x) * ((double)p.x - q.x) + ((double)p.y - q.y) * ((double)p.y - q.y) + ((double)p.z - q.z) * ((double)p.z - q.z);
				if (dd < bd)
				{
					bd = dd;
					best = j;
				}
			}
			for (k = 0; k < 3; k++) err += (rgb[3 * best + k] - rgb0[3 * i + k]) * (rgb[3 * best + k] - rgb0[3 * i + k]) / (255.0 * 255.0);
		}
		err = sqrt(err / (3.0 * nv0));
		sprintf(sub, "%u faces; colour error %.3f", rep.faces_after, err);
		ms_panel(c, 750 + 750 * (int)r, 0, 750, 720, m, 0.5, 0.9, r ? "with the colour in the quadric" : "by the shape alone", sub);
		ms_draw(c, m, NULL, rgb, 0, NULL, 0);
		ms_draw(c, m, NULL, NULL, 0, ms_da_wire, 0.6);
		printf("  %s: %s\n", r ? "with colour" : "shape alone", sub);
		free(rgb);
		qaws_mesh_3d_destroy(m);
	}
	er_save(c, "showcase/mesh51_decimate_attributes.png");
	printf("mesh51_decimate_attributes.png\n");
	free(rgb0);
	free(x);
	free(t);
	qaws_mesh_3d_destroy(m0);
	er_destroy(c);
}
