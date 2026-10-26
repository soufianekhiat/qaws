/*
 * Figure 61: attributes carried through Booleans and sweeps
 * (qaws_mesh_3d_boolean, qaws_mesh_3d_sweep, qaws_mesh_sweep_boolean).
 *   - a ball with sRGB vertex colours (blended in linear light) and radial
 *     NORMAL normals less a box of its own colour: the cut stays sharp, each
 *     side its own colours (the box's carved faces' normals turned with them);
 *   - UVs per corner: a ball's latitude-longitude chart less a smaller
 *     ball with its own, a checker drawn from the carried UVs;
 *   - a coloured ball swept along a curve: the tube carries the ball's
 *     colours and normals (turned with the motion);
 *   - carving: the slot's faces carry the tool's material.
 * Printed: faces, the attributes found and their elements.
 */

/* a mesh drawn with per-corner colours (or per face: fc) and per-corner normals (NULL: flat) */
static void ms_at_draw(er_canvas* c, qaws_mesh_3d const* m, double const* cc, double const* fc, double const* cn, double const* wire, double wire_width)
{
	unsigned int nt = qaws_mesh_3d_get_face_count(m), ne = qaws_mesh_3d_get_edge_count(m), i;
	double* x = ms_positions(m);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (nt ? nt : 1));
	unsigned int* e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (ne ? ne : 1));
	qaws_mesh_3d_get_triangles(t, m);
	for (i = 0; i < nt; i++)
		er_triangle(c, x + 3 * t[3 * i], x + 3 * t[3 * i + 1], x + 3 * t[3 * i + 2], cn ? cn + 9 * i : NULL, cn ? cn + 9 * i + 3 : NULL,
			cn ? cn + 9 * i + 6 : NULL, fc ? fc + 3 * i : cc + 9 * i, fc ? fc + 3 * i : cc + 9 * i + 3, fc ? fc + 3 * i : cc + 9 * i + 6);
	if (wire)
	{
		qaws_mesh_3d_get_edges(e, m);
		for (i = 0; i < ne; i++)
			er_line(c, x + 3 * e[2 * i], x + 3 * e[2 * i + 1], wire, wire_width);
	}
	free(x); free(t); free(e);
}

/* an attribute's values as doubles (NULL: absent) */
static double* ms_at_get(qaws_mesh_3d const* m, char const* name, unsigned int* element)
{
	unsigned int id, count, i;
	qaws_attribute_desc d;
	qaws_scalar* v;
	double* o;
	if (qaws_mesh_3d_find_attribute(&id, m, name) != QAWS_STATUS_OK)
		return NULL;
	qaws_mesh_3d_get_attribute_desc(&d, m, id);
	count = d.components * (d.element == QAWS_MESH_VERTEX ? qaws_mesh_3d_get_vertex_count(m) : d.element == QAWS_MESH_CORNER ? qaws_mesh_3d_get_corner_count(m) : d.element == QAWS_MESH_EDGE ? qaws_mesh_3d_get_edge_count(m) : qaws_mesh_3d_get_face_count(m));
	v = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (count + 1));
	o = (double*)malloc(sizeof(double) * (count + 1));
	qaws_mesh_3d_get_attribute(v, m, id);
	for (i = 0; i < count; i++) o[i] = (double)v[i];
	free(v);
	*element = (unsigned int)d.element;
	return o;
}

static void ms_at_add(qaws_mesh_3d* m, char const* name, qaws_mesh_element e, unsigned int n, qaws_attribute_rule rule, double const* v, unsigned int count)
{
	qaws_attribute_desc d;
	unsigned int id, i;
	qaws_scalar* s = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (n * count + 1));
	for (i = 0; i < n * count; i++) s[i] = (qaws_scalar)v[i];
	d.name = name; d.element = e; d.components = n; d.rule = rule;
	qaws_mesh_3d_add_attribute(&id, m, s, &d);
	free(s);
}

/* a ball of radius r at o: sRGB colours from a hue wheel round z, radial normals */
static qaws_mesh_3d* ms_at_ball(double r, double ox, double oy, double oz)
{
	qaws_mesh_3d* b = ms_rm_sphere(32, 16);
	unsigned int nv = b ? qaws_mesh_3d_get_vertex_count(b) : 0, i;
	double *col = (double*)malloc(sizeof(double) * 3 * (nv + 1)), *nrm = (double*)malloc(sizeof(double) * 3 * (nv + 1));
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (nv + 1));
	if (!b)
		return NULL;
	qaws_mesh_3d_get_positions(x, b);
	for (i = 0; i < nv; i++)
	{
		double p[3] = { (double)x[3 * i], (double)x[3 * i + 1], (double)x[3 * i + 2] }, h = atan2(p[1], p[0]), k;
		unsigned int j;
		for (j = 0; j < 3; j++)
		{
			k = 0.5 + 0.5 * cos(h - 2.0943951 * j);
			col[3 * i + j] = 0.15 + 0.8 * k * (0.6 + 0.4 * p[2]);
			nrm[3 * i + j] = p[j];
		}
		x[3 * i] = (qaws_scalar)(r * p[0] + ox);
		x[3 * i + 1] = (qaws_scalar)(r * p[1] + oy);
		x[3 * i + 2] = (qaws_scalar)(r * p[2] + oz);
	}
	qaws_mesh_3d_set_positions(b, x);
	ms_at_add(b, "color", QAWS_MESH_VERTEX, 3, QAWS_ATTRIBUTE_SRGB, col, nv);
	ms_at_add(b, "normal", QAWS_MESH_VERTEX, 3, QAWS_ATTRIBUTE_NORMAL, nrm, nv);
	free(col); free(nrm); free(x);
	return b;
}

static void ms_at_face_value(qaws_mesh_3d* m, char const* name, unsigned int n, double const* v)
{
	unsigned int nf = qaws_mesh_3d_get_face_count(m), i, k;
	double* a = (double*)malloc(sizeof(double) * (n * nf + 1));
	for (i = 0; i < nf; i++)
		for (k = 0; k < n; k++) a[n * i + k] = v[k];
	ms_at_add(m, name, QAWS_MESH_FACE, n, n == 1 ? QAWS_ATTRIBUTE_NEAREST : QAWS_ATTRIBUTE_SRGB, a, nf);
	free(a);
}

/* per-corner colours (0..255) of a carried "color" (sRGB, per corner or face) */
static double* ms_at_colours(qaws_mesh_3d const* m, unsigned int* el)
{
	unsigned int nt = qaws_mesh_3d_get_face_count(m), i;
	double *v = ms_at_get(m, "color", el), *o;
	if (!v)
		return NULL;
	o = (double*)malloc(sizeof(double) * 9 * (nt + 1));
	for (i = 0; i < 9 * nt; i++)
		o[i] = 255 * fmin(1, fmax(0, *el == QAWS_MESH_FACE ? v[3 * (i / 9) + i % 3] : v[i]));
	free(v);
	return o;
}

static void demo_attributes(void)
{
	static double const wire[3] = { 60, 66, 80 };
	er_canvas* c = er_create(2400, 640, 3, 252, 252, 250);
	char sub[200];
	if (!c)
		return;
	{
		/* a coloured ball less a box of its own colour */
		static double const orange[3] = { 0.95, 0.55, 0.2 };
		qaws_mesh_3d *ball = ms_at_ball(1, 0, 0, 0), *box = ms_sw_box(0, 1.4, -0.2, 1.4, 0, 1.4), *o = NULL;
		unsigned int ec = 0, en = 0;
		if (box)
		{
			/* (per vertex, as the ball's: one name, one element) */
			unsigned int nv = qaws_mesh_3d_get_vertex_count(box), i;
			double* v = (double*)malloc(sizeof(double) * 3 * nv);
			for (i = 0; i < 3 * nv; i++) v[i] = orange[i % 3];
			ms_at_add(box, "color", QAWS_MESH_VERTEX, 3, QAWS_ATTRIBUTE_SRGB, v, nv);
			free(v);
		}
		if (ball && box && qaws_mesh_3d_boolean(&o, NULL, ball, box, QAWS_BOOLEAN_DIFFERENCE, NULL) == QAWS_STATUS_OK && o)
		{
			double *cc = ms_at_colours(o, &ec), *cn = ms_at_get(o, "normal", &en);
			sprintf(sub, "%u faces; colour, normals per %s", qaws_mesh_3d_get_face_count(o), en == QAWS_MESH_CORNER ? "corner" : "?");
			ms_panel(c, 0, 0, 600, 640, o, -0.6, 0.75, "colours and normals carried", sub);
			ms_at_draw(c, o, cc, NULL, cn, NULL, 0);
			printf("  colours and normals: %s\n", sub);
			free(cc); free(cn);
		}
		qaws_mesh_3d_destroy(ball); qaws_mesh_3d_destroy(box); qaws_mesh_3d_destroy(o);
	}
	{
		/* UVs: a ball's latitude-longitude chart less a smaller ball's, a checker from them */
		qaws_mesh_3d *ball = ms_at_ball(1, 0, 0, 0), *box = ms_at_ball(0.75, 0.55, 0.45, 0.6), *o = NULL;
		unsigned int i, j, el = 0;
		{
			unsigned int nf = qaws_mesh_3d_get_face_count(ball), nc = qaws_mesh_3d_get_corner_count(ball);
			unsigned int* fv = (unsigned int*)malloc(sizeof(unsigned int) * nc);
			double* uv = (double*)malloc(sizeof(double) * 2 * nc);
			qaws_mesh_3d_get_faces(fv, NULL, ball);
			for (i = 0; i < nf; i++)
			{
				double lo = 1, hi = 0;
				for (j = 0; j < 3; j++)
				{
					qaws_vec3 p = qaws_mesh_3d_get_vertex(ball, fv[3 * i + j]);
					double l = sqrt((double)p.x * p.x + (double)p.y * p.y + (double)p.z * p.z);
					uv[6 * i + 2 * j] = 0.5 + atan2((double)p.y, (double)p.x) / (2 * 3.14159265358979323846);
					uv[6 * i + 2 * j + 1] = acos((double)p.z / l) / 3.14159265358979323846;
					lo = fmin(lo, uv[6 * i + 2 * j]);
					hi = fmax(hi, uv[6 * i + 2 * j]);
				}
				for (j = 0; hi - lo > 0.5 && j < 3; j++)
					if (uv[6 * i + 2 * j] < 0.5) uv[6 * i + 2 * j] += 1;   /* (the seam: one side of the chart) */
			}
			ms_at_add(ball, "uv", QAWS_MESH_CORNER, 2, QAWS_ATTRIBUTE_LINEAR, uv, nc);
			free(uv);
			nc = qaws_mesh_3d_get_corner_count(box);
			nf = qaws_mesh_3d_get_face_count(box);
			fv = (unsigned int*)realloc(fv, sizeof(unsigned int) * nc);
			uv = (double*)malloc(sizeof(double) * 2 * nc);
			qaws_mesh_3d_get_faces(fv, NULL, box);
			for (i = 0; i < nf; i++)
			{
				double lo = 1, hi = 0;
				for (j = 0; j < 3; j++)
				{
					qaws_vec3 p = qaws_mesh_3d_get_vertex(box, fv[3 * i + j]);
					double q[3] = { (double)p.x - 0.55, (double)p.y - 0.45, (double)p.z - 0.6 }, l = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);
					uv[6 * i + 2 * j] = 0.5 + atan2(q[1], q[0]) / (2 * 3.14159265358979323846);
					uv[6 * i + 2 * j + 1] = acos(q[2] / l) / 3.14159265358979323846;
					lo = fmin(lo, uv[6 * i + 2 * j]);
					hi = fmax(hi, uv[6 * i + 2 * j]);
				}
				for (j = 0; hi - lo > 0.5 && j < 3; j++)
					if (uv[6 * i + 2 * j] < 0.5) uv[6 * i + 2 * j] += 1;
			}
			ms_at_add(box, "uv", QAWS_MESH_CORNER, 2, QAWS_ATTRIBUTE_LINEAR, uv, nc);
			free(uv);
			free(fv);
		}
		if (ball && box && qaws_mesh_3d_boolean(&o, NULL, ball, box, QAWS_BOOLEAN_DIFFERENCE, NULL) == QAWS_STATUS_OK && o)
		{
			unsigned int nt = qaws_mesh_3d_get_face_count(o), es = 0;
			double *uv = ms_at_get(o, "uv", &el), *src = ms_at_get(o, "source", &es), *fc = (double*)malloc(sizeof(double) * 3 * (nt + 1));
			unsigned int nball = qaws_mesh_3d_get_face_count(ball);
			for (i = 0; uv && src && i < nt; i++)
			{
				/* the checker at the face's centre in its own chart: the ball's 16 x 8, the bite's 24 x 12 */
				double u = (uv[6 * i] + uv[6 * i + 2] + uv[6 * i + 4]) / 3, v = (uv[6 * i + 1] + uv[6 * i + 3] + uv[6 * i + 5]) / 3;
				int on_ball = (unsigned int)src[i] < nball, k = ((int)floor(u * (on_ball ? 16 : 24)) + (int)floor(v * (on_ball ? 8 : 12))) & 1;
				fc[3 * i] = on_ball ? (k ? 235 : 70) : (k ? 245 : 200);
				fc[3 * i + 1] = on_ball ? (k ? 235 : 90) : (k ? 200 : 110);
				fc[3 * i + 2] = on_ball ? (k ? 240 : 150) : (k ? 150 : 50);
			}
			sprintf(sub, "%u faces; uv per %s, a checker from it", nt, el == QAWS_MESH_CORNER ? "corner" : "?");
			ms_panel(c, 600, 0, 600, 640, o, -0.6, 0.75, "UVs carried", sub);
			ms_at_draw(c, o, NULL, fc, NULL, NULL, 0);
			printf("  UVs: %s\n", sub);
			free(uv); free(src); free(fc);
		}
		qaws_mesh_3d_destroy(ball); qaws_mesh_3d_destroy(box); qaws_mesh_3d_destroy(o);
	}
	{
		/* a coloured ball swept along a Hermite curve */
		qaws_mesh_3d *ball = ms_at_ball(0.3, 0, 0, 0), *o = NULL;
		qaws_motion_pose s, e;
		qaws_motion* m = NULL;
		qaws_sweep_desc d;
		qaws_sweep_report r;
		unsigned int ec = 0, en = 0;
		memset(&s, 0, sizeof s);
		memset(&e, 0, sizeof e);
		s.rotation[0] = 1;
		e.rotation[0] = (qaws_scalar)cos(0.4 * 3.14159265358979323846);
		e.rotation[3] = (qaws_scalar)sin(0.4 * 3.14159265358979323846);
		e.position.x = 3; e.position.z = 0.8;
		s.velocity.x = 1; s.velocity.y = 4;
		e.velocity.x = 1; e.velocity.y = -4; e.velocity.z = 2;
		s.angular_velocity.z = e.angular_velocity.z = (qaws_scalar)(0.8 * 3.14159265358979323846);
		qaws_motion_create_poses(&m, &s, &e, 0);
		memset(&d, 0, sizeof d);
		d.tolerance = (qaws_scalar)0.005;
		if (ball && m && qaws_mesh_3d_sweep(&o, &r, ball, m, &d) == QAWS_STATUS_OK && o)
		{
			double *cc = ms_at_colours(o, &ec), *cn = ms_at_get(o, "normal", &en);
			sprintf(sub, "%u slabs, %u faces; colour, normals per %s", r.slabs, qaws_mesh_3d_get_face_count(o), en == QAWS_MESH_CORNER ? "corner" : "?");
			ms_panel(c, 1200, 0, 600, 640, o, -0.3, 0.8, "a sweep carries them", sub);
			ms_at_draw(c, o, cc, NULL, cn, NULL, 0);
			printf("  sweep: %s\n", sub);
			free(cc); free(cn);
		}
		qaws_motion_destroy(m);
		qaws_mesh_3d_destroy(ball); qaws_mesh_3d_destroy(o);
	}
	{
		/* carving: the slot takes the tool's material (its colour) */
		static double const wood[3] = { 0.84, 0.70, 0.51 }, red[3] = { 0.85, 0.25, 0.22 };
		qaws_mesh_3d *block = ms_sw_box(0, 4, 0, 2, 0, 1), *tool = ms_sw_box(-0.3, 0.3, -0.3, 0.3, -0.3, 0.3), *o = NULL;
		qaws_motion_pose s, e;
		qaws_motion *mv[2] = { NULL, NULL };
		qaws_motion const* cm[2];
		qaws_sweep_boolean_report br;
		unsigned int ec = 0;
		ms_at_face_value(block, "color", 3, wood);
		ms_at_face_value(tool, "color", 3, red);
		memset(&s, 0, sizeof s);
		memset(&e, 0, sizeof e);
		s.rotation[0] = e.rotation[0] = 1;
		s.position.x = -0.5; s.position.y = 0.6; s.position.z = 1.1;
		e.position.x = 4.5; e.position.y = 0.6; e.position.z = 1.1;
		s.velocity.x = e.velocity.x = 5;
		qaws_motion_create_poses(&mv[0], &s, &e, 0);
		s.position.x = 3; s.position.y = -0.5; s.position.z = 1;
		e.position.x = 3; e.position.y = 2.5; e.position.z = 1;
		s.velocity.x = e.velocity.x = 0;
		s.velocity.y = e.velocity.y = 3;
		qaws_motion_create_poses(&mv[1], &s, &e, 0);
		cm[0] = mv[0];
		cm[1] = mv[1];
		if (block && tool && mv[0] && mv[1] && qaws_mesh_sweep_boolean(&o, &br, tool, cm, block, 2, QAWS_BOOLEAN_DIFFERENCE, NULL) == QAWS_STATUS_OK && o)
		{
			double* cc = ms_at_colours(o, &ec);
			sprintf(sub, "%u moves, %u faces; colour per %s", br.moves, qaws_mesh_3d_get_face_count(o), ec == QAWS_MESH_FACE ? "face" : "?");
			ms_panel(c, 1800, 0, 600, 640, o, -0.4, 0.8, "carving: the tool's colour", sub);
			ms_at_draw(c, o, cc, NULL, NULL, wire, 0.3);
			printf("  carving: %s\n", sub);
			free(cc);
		}
		qaws_motion_destroy(mv[0]); qaws_motion_destroy(mv[1]);
		qaws_mesh_3d_destroy(block); qaws_mesh_3d_destroy(tool); qaws_mesh_3d_destroy(o);
	}
	er_save(c, "showcase/mesh61_attributes.png");
	printf("mesh61_attributes.png\n");
	er_destroy(c);
}
