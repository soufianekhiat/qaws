/*
 * Figure 63: sweep schemes and attribute rules by the type of data.
 *   Top row, one coloured cube tool carving one block (qaws_mesh_sweep_boolean)
 *   by four motion schemes, the slot taking the tool's colours (Oklab):
 *   - linear: two poses, the velocity the displacement, no turn;
 *   - cubic Hermite: the same ends with side velocities and a quarter turn;
 *   - quintic Hermite: the same ends and velocities, accelerations added;
 *   - along a Bezier curve: a constant rotation.
 *   The tool's path drawn over each, the tool at t = 0, 1/2 and 1.
 *   Bottom row, the same data under different rules:
 *   - a box less a ball, the ball's radial directions given twice: as a
 *     NORMAL (blue, negated on the bite: it points out of the result) and as
 *     a VECTOR (orange, kept: it still points away from the ball's centre);
 *   - a box swept through a quarter turn: a VECTOR (orange) turns with the
 *     motion, the same values LINEAR (grey) stay as given;
 *   - one quad subdivided (Catmull-Clark, 4 levels), its left side the
 *     identity, its right side minus a third of a turn about z: as a
 *     QUATERNION the frames stay rigid and turn the short way, as LINEAR
 *     values they shrink through the middle and turn the long way;
 *   - a QUATERNION through the cubic sweep: each frame the tool's orientation
 *     where that piece of the solid was swept.
 * Printed: slabs, removed volumes, glyph counts.
 */

/* the faces' centroids and unit normals (3 each, caller frees) */
static void ms_at3_frames(qaws_mesh_3d const* m, double** cen, double** nrm)
{
	unsigned int nt = qaws_mesh_3d_get_face_count(m), i, k;
	double* x = ms_positions(m);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (nt + 1));
	*cen = (double*)malloc(sizeof(double) * 3 * (nt + 1));
	*nrm = (double*)malloc(sizeof(double) * 3 * (nt + 1));
	qaws_mesh_3d_get_triangles(t, m);
	for (i = 0; i < nt; i++)
	{
		double const *a = x + 3 * t[3 * i], *b = x + 3 * t[3 * i + 1], *d = x + 3 * t[3 * i + 2];
		double u[3], w[3], n[3], l;
		for (k = 0; k < 3; k++)
		{
			(*cen)[3 * i + k] = (a[k] + b[k] + d[k]) / 3;
			u[k] = b[k] - a[k];
			w[k] = d[k] - a[k];
		}
		n[0] = u[1] * w[2] - u[2] * w[1]; n[1] = u[2] * w[0] - u[0] * w[2]; n[2] = u[0] * w[1] - u[1] * w[0];
		l = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
		for (k = 0; k < 3; k++) (*nrm)[3 * i + k] = l > 0 ? n[k] / l : 0;
	}
	free(x); free(t);
}

/* faces kept (keep[f]) one per cell of a grid of the given size: their indices, the count returned */
static unsigned int ms_at3_pick(unsigned int* out, double const* cen, unsigned char const* keep, unsigned int nt, double cell)
{
	int* used = (int*)malloc(sizeof(int) * 3 * (nt + 1));
	unsigned int n = 0, i, j;
	for (i = 0; i < nt; i++)
	{
		int k[3], dup = 0;
		if (!keep[i])
			continue;
		for (j = 0; j < 3; j++) k[j] = (int)floor(cen[3 * i + j] / cell);
		for (j = 0; j < n && !dup; j++)
			dup = used[3 * j] == k[0] && used[3 * j + 1] == k[1] && used[3 * j + 2] == k[2];
		if (dup)
			continue;
		memcpy(used + 3 * n, k, sizeof k);
		out[n++] = i;
	}
	free(used);
	return n;
}

/* an arrow from p along v (scaled by len) */
static void ms_at3_arrow(er_canvas* c, double const* p, double const* v, double len, double const* col)
{
	double q[3];
	unsigned int k;
	for (k = 0; k < 3; k++) q[k] = p[k] + len * v[k];
	er_line(c, p, q, col, 3);
	er_point(c, q, col, 4.5);
}

/* a frame from a quaternion (w x y z), not made unit: its axes scaled by |q|^2 (x red, y green, z blue) */
static void ms_at3_triad(er_canvas* c, double const* p, double const* q, double len)
{
	static double const col[3][3] = { { 215, 50, 50 }, { 40, 160, 70 }, { 50, 90, 220 } };
	double w = q[0], x = q[1], y = q[2], z = q[3], R[9];
	unsigned int a, k;
	R[0] = w * w + x * x - y * y - z * z; R[1] = 2 * (x * y - w * z); R[2] = 2 * (x * z + w * y);
	R[3] = 2 * (x * y + w * z); R[4] = w * w - x * x + y * y - z * z; R[5] = 2 * (y * z - w * x);
	R[6] = 2 * (x * z - w * y); R[7] = 2 * (y * z + w * x); R[8] = w * w - x * x - y * y + z * z;
	for (a = 0; a < 3; a++)
	{
		double e[3];
		for (k = 0; k < 3; k++) e[k] = p[k] + len * R[3 * k + a];
		er_line(c, p, e, col[a], 3);
	}
	er_point(c, p, col[0], 2.5);
}

/* the mean of a per-corner attribute over triangle f's corners */
static void ms_at3_mean(double* out, double const* v, unsigned int n, unsigned int f)
{
	unsigned int k;
	for (k = 0; k < n; k++) out[k] = (v[n * (3 * f) + k] + v[n * (3 * f + 1) + k] + v[n * (3 * f + 2) + k]) / 3;
}

static void ms_at3_box_edges(er_canvas* c, double const* x, double const* col)
{
	static unsigned int const e[24] = { 0, 1, 1, 2, 2, 3, 3, 0, 4, 5, 5, 6, 6, 7, 7, 4, 0, 4, 1, 5, 2, 6, 3, 7 };
	unsigned int i;
	for (i = 0; i < 12; i++) er_line(c, x + 3 * e[2 * i], x + 3 * e[2 * i + 1], col, 2);
}

/* the cube tool of half size h: eight Oklab-blended colours (stored sRGB), the identity quaternion */
static qaws_mesh_3d* ms_at3_tool(double h)
{
	static double const col[24] = { 0.95, 0.25, 0.30, 0.98, 0.70, 0.15, 0.30, 0.75, 0.35, 0.15, 0.55, 0.95,
		0.75, 0.30, 0.85, 0.20, 0.80, 0.85, 0.98, 0.95, 0.30, 0.95, 0.45, 0.65 };
	double q[32];
	unsigned int i;
	qaws_mesh_3d* m = ms_sw_box(-h, h, -h, h, -h, h);
	for (i = 0; i < 8; i++)
	{
		q[4 * i] = 1; q[4 * i + 1] = q[4 * i + 2] = q[4 * i + 3] = 0;
	}
	ms_at_add(m, "color", QAWS_MESH_VERTEX, 3, QAWS_ATTRIBUTE_OKLAB, col, 8);
	ms_at_add(m, "q", QAWS_MESH_VERTEX, 4, QAWS_ATTRIBUTE_QUATERNION, q, 8);
	return m;
}

/* the poses of the top row's schemes: 0 linear, 1 cubic, 2 quintic */
static qaws_motion* ms_at3_poses(int scheme)
{
	qaws_motion_pose s, e;
	qaws_motion* m = NULL;
	memset(&s, 0, sizeof s);
	memset(&e, 0, sizeof e);
	s.rotation[0] = 1;
	e.rotation[0] = 1;
	s.position.x = -0.5; s.position.y = 1; s.position.z = 1.05;
	e.position.x = 4.5; e.position.y = 1; e.position.z = 1.05;
	s.velocity.x = e.velocity.x = 5;
	if (scheme > 0)
	{
		s.velocity.y = 2.4; e.velocity.y = -2.4;
		e.rotation[0] = (qaws_scalar)cos(MS_PI / 4); e.rotation[3] = (qaws_scalar)sin(MS_PI / 4);
		s.angular_velocity.z = e.angular_velocity.z = (qaws_scalar)(MS_PI / 2);
	}
	if (scheme == 2)
		s.acceleration.y = e.acceleration.y = -40;
	qaws_motion_create_poses(&m, &s, &e, scheme == 2 ? QAWS_MOTION_QUINTIC : 0);
	return m;
}

/* the motion borrows the curve: *cv kept until the motion is destroyed */
static qaws_motion* ms_at3_curve(qaws_curve** cv)
{
	qaws_scalar cp[12] = { -0.5f, 0.35f, 1.05f, 1.6f, 2.6f, 1.05f, 2.4f, -0.6f, 1.05f, 4.5f, 1.65f, 1.05f };
	qaws_scalar rot[4];
	qaws_bezier_desc bd;
	qaws_motion* m = NULL;
	rot[0] = (qaws_scalar)cos(MS_PI / 8); rot[1] = rot[2] = 0; rot[3] = (qaws_scalar)sin(MS_PI / 8);
	memset(&bd, 0, sizeof bd);
	bd.dimension = QAWS_DIMENSION_3D; bd.degree = 3; bd.control_points = cp; bd.control_point_count = 4;
	*cv = NULL;
	if (qaws_curve_create_bezier(cv, &bd) == QAWS_STATUS_OK)
		qaws_motion_create_curve(&m, *cv, rot);
	return m;
}

static void ms_at3_carve(er_canvas* c, int px, qaws_mesh_3d const* tool, qaws_mesh_3d const* block, qaws_motion const* m, char const* title)
{
	static double const ink[3] = { 40, 40, 48 }, ghost[3] = { 70, 70, 90 }, wire[3] = { 90, 80, 70 };
	qaws_mesh_3d* o = NULL;
	qaws_motion const* mv[1];
	qaws_sweep_desc d;
	qaws_sweep_boolean_report r;
	qaws_scalar tx[24], ox[24];
	unsigned int el, i, k;
	double* cc;
	char sub[200];
	mv[0] = m;
	memset(&d, 0, sizeof d);
	d.tolerance = (qaws_scalar)0.004;
	{
		qaws_status st = m ? qaws_mesh_sweep_boolean(&o, &r, tool, mv, block, 1, QAWS_BOOLEAN_DIFFERENCE, &d) : QAWS_STATUS_INVALID_ARGUMENT;
		if (st != QAWS_STATUS_OK || !o)
		{
			printf("  %s: failed (%d)\n", title, (int)st);
			return;
		}
	}
	sprintf(sub, "%u slabs, removed %.3f, %u faces", r.slabs, (double)r.changed_volume, qaws_mesh_3d_get_face_count(o));
	{
		/* the view fitted to the result and the tool at its ends */
		unsigned int nv = qaws_mesh_3d_get_vertex_count(o);
		double* x = ms_positions(o);
		double* fit = (double*)malloc(sizeof(double) * (3 * (size_t)nv + 48));
		memcpy(fit, x, sizeof(double) * 3 * nv);
		qaws_mesh_3d_get_positions(tx, tool);
		for (i = 0; i < 2; i++)
		{
			qaws_motion_transform_points(ox, NULL, m, (qaws_scalar)i, tx, 8);
			for (k = 0; k < 24; k++) fit[3 * nv + 24 * i + k] = (double)ox[k];
		}
		er_panel(c, px, 72, 600, 640 - 72, fit, nv + 16, -0.35, 0.95, 0.04);
		er_text(c, px + 14, 12, title, 3, 30, 30, 30);
		er_text(c, px + 14, 46, sub, 2, 95, 95, 95);
		free(x); free(fit);
	}
	cc = ms_at_colours(o, &el);
	if (cc)
		ms_at_draw(c, o, cc, NULL, NULL, wire, 0.25);
	free(cc);
	/* the path (the tool's centre) lifted above the block, the tool at 0, 1/2 and 1 */
	{
		qaws_scalar zero[3] = { 0, 0, 0 }, p[3];
		double a[3], b[3];
		for (i = 0; i <= 64; i++)
		{
			qaws_motion_transform_points(p, NULL, m, (qaws_scalar)(i / 64.0), zero, 1);
			for (k = 0; k < 3; k++) b[k] = (double)p[k] + (k == 2 ? 0.45 : 0);
			if (i)
				er_line(c, a, b, ink, 2.5);
			memcpy(a, b, sizeof a);
		}
	}
	qaws_mesh_3d_get_positions(tx, tool);
	for (i = 0; i < 3; i++)
	{
		double x[24];
		qaws_motion_transform_points(ox, NULL, m, (qaws_scalar)(i / 2.0), tx, 8);
		for (k = 0; k < 24; k++) x[k] = (double)ox[k];
		ms_at3_box_edges(c, x, ghost);
	}
	printf("  %s: %s\n", title, sub);
	qaws_mesh_3d_destroy(o);
}

static void demo_sweep_schemes(void)
{
	static double const orange[3] = { 235, 120, 20 }, blue[3] = { 40, 90, 225 }, grey[3] = { 120, 120, 120 }, wire[3] = { 110, 110, 120 };
	er_canvas* c = er_create(2400, 1280, 3, 252, 252, 250);
	qaws_mesh_3d *tool = ms_at3_tool(0.3), *block = ms_sw_box(0, 4, 0, 2, 0, 1);
	qaws_motion* m[4];
	qaws_curve* path = NULL;
	unsigned int i, k;
	if (!c || !tool || !block)
		return;
	{
		double wood[24];
		for (i = 0; i < 8; i++)
		{
			wood[3 * i] = 0.84; wood[3 * i + 1] = 0.72; wood[3 * i + 2] = 0.55;
		}
		ms_at_add(block, "color", QAWS_MESH_VERTEX, 3, QAWS_ATTRIBUTE_OKLAB, wood, 8);
	}
	m[0] = ms_at3_poses(0);
	m[1] = ms_at3_poses(1);
	m[2] = ms_at3_poses(2);
	m[3] = ms_at3_curve(&path);
	ms_at3_carve(c, 0, tool, block, m[0], "linear: two poses");
	ms_at3_carve(c, 600, tool, block, m[1], "cubic Hermite, a quarter turn");
	ms_at3_carve(c, 1200, tool, block, m[2], "quintic: accelerations too");
	ms_at3_carve(c, 1800, tool, block, m[3], "along a Bezier curve");

	{
		/* NORMAL against VECTOR through a bite */
		qaws_mesh_3d *a = ms_sw_box(0, 2, 0, 2, 0, 1), *b = ms_at_ball(0.8, 1, 1, 1.15), *o = NULL;
		unsigned int nv = b ? qaws_mesh_3d_get_vertex_count(b) : 0, el;
		double grey8[24];
		for (i = 0; i < 24; i++) grey8[i] = 0.78;
		ms_at_add(a, "color", QAWS_MESH_VERTEX, 3, QAWS_ATTRIBUTE_SRGB, grey8, 8);
		if (b)
		{
			double* nb = ms_at_get(b, "normal", &el);
			ms_at_add(b, "vector", QAWS_MESH_VERTEX, 3, QAWS_ATTRIBUTE_VECTOR, nb, nv);
			free(nb);
		}
		if (b && qaws_mesh_3d_boolean(&o, NULL, a, b, QAWS_BOOLEAN_DIFFERENCE, NULL) == QAWS_STATUS_OK && o)
		{
			unsigned int nt = qaws_mesh_3d_get_face_count(o), np, *pick = (unsigned int*)malloc(sizeof(unsigned int) * (nt + 1));
			double *cc = ms_at_colours(o, &el), *nn = ms_at_get(o, "normal", &el), *vv = ms_at_get(o, "vector", &el), *src = ms_at_get(o, "source", &el), *cen, *fn;
			unsigned char* keep = (unsigned char*)malloc(nt + 1);
			char sub[200];
			double bias;
			ms_at3_frames(o, &cen, &fn);
			for (i = 0; i < nt; i++) keep[i] = src && src[i] >= 12 && fn[3 * i + 2] > 0.2;
			np = ms_at3_pick(pick, cen, keep, nt, 0.28);
			sprintf(sub, "blue NORMAL negated, orange VECTOR kept");
			ms_panel(c, 0, 640, 600, 640, o, -0.5, 1.05, "NORMAL and VECTOR in a bite", sub);
			ms_at_draw(c, o, cc, NULL, NULL, wire, 0.2);
			bias = c->depth_bias;
			c->depth_bias = 1e30;   /* (drawn over: the VECTORs point into the solid) */
			for (i = 0; nn && vv && i < np; i++)
			{
				double p[3], n3[3], v3[3];
				for (k = 0; k < 3; k++) p[k] = cen[3 * pick[i] + k] + 0.01 * fn[3 * pick[i] + k];
				ms_at3_mean(n3, nn, 3, pick[i]);
				ms_at3_mean(v3, vv, 3, pick[i]);
				ms_at3_arrow(c, p, n3, 0.22, blue);
				ms_at3_arrow(c, p, v3, 0.22, orange);
			}
			c->depth_bias = bias;
			printf("  bite: %s\n", sub);
			free(pick); free(cc); free(nn); free(vv); free(src); free(cen); free(fn); free(keep);
		}
		qaws_mesh_3d_destroy(a); qaws_mesh_3d_destroy(b); qaws_mesh_3d_destroy(o);
	}
	{
		/* VECTOR against LINEAR through a quarter turn */
		static double const cols[24] = { 0.25, 0.55, 0.95, 0.95, 0.45, 0.30, 0.95, 0.45, 0.30, 0.25, 0.55, 0.95,
			0.25, 0.55, 0.95, 0.95, 0.45, 0.30, 0.95, 0.45, 0.30, 0.25, 0.55, 0.95 };
		qaws_mesh_3d *box = ms_sw_box(1, 2.4, -0.15, 0.15, 0, 0.5), *o = NULL;
		qaws_motion_pose s, e;
		qaws_motion* mq = NULL;
		double ex[24];
		for (i = 0; i < 8; i++)
		{
			ex[3 * i] = 1; ex[3 * i + 1] = 0; ex[3 * i + 2] = 0;
		}
		ms_at_add(box, "color", QAWS_MESH_VERTEX, 3, QAWS_ATTRIBUTE_OKLAB, cols, 8);
		ms_at_add(box, "vector", QAWS_MESH_VERTEX, 3, QAWS_ATTRIBUTE_VECTOR, ex, 8);
		ms_at_add(box, "linear", QAWS_MESH_VERTEX, 3, QAWS_ATTRIBUTE_LINEAR, ex, 8);
		memset(&s, 0, sizeof s);
		memset(&e, 0, sizeof e);
		s.rotation[0] = 1;
		e.rotation[0] = (qaws_scalar)cos(MS_PI / 4); e.rotation[3] = (qaws_scalar)sin(MS_PI / 4);
		s.angular_velocity.z = e.angular_velocity.z = (qaws_scalar)(MS_PI / 2);
		qaws_motion_create_poses(&mq, &s, &e, 0);
		if (mq && qaws_mesh_3d_sweep(&o, NULL, box, mq, NULL) == QAWS_STATUS_OK && o)
		{
			unsigned int nt = qaws_mesh_3d_get_face_count(o), np, el, *pick = (unsigned int*)malloc(sizeof(unsigned int) * (nt + 1));
			double *cc = ms_at_colours(o, &el), *vv = ms_at_get(o, "vector", &el), *ll = ms_at_get(o, "linear", &el), *cen, *fn;
			unsigned char* keep = (unsigned char*)malloc(nt + 1);
			char sub[200];
			ms_at3_frames(o, &cen, &fn);
			for (i = 0; i < nt; i++) keep[i] = fn[3 * i + 2] > 0.99;
			np = ms_at3_pick(pick, cen, keep, nt, 0.3);
			sprintf(sub, "orange VECTOR turned, grey LINEAR as given");
			ms_panel(c, 600, 640, 600, 640, o, -0.6, 0.75, "VECTOR and LINEAR swept", sub);
			ms_at_draw(c, o, cc, NULL, NULL, wire, 0.2);
			for (i = 0; vv && ll && i < np; i++)
			{
				double p[3], v3[3], l3[3];
				for (k = 0; k < 3; k++) p[k] = cen[3 * pick[i] + k] + (k == 2 ? 0.01 : 0);
				ms_at3_mean(v3, vv, 3, pick[i]);
				ms_at3_mean(l3, ll, 3, pick[i]);
				ms_at3_arrow(c, p, l3, 0.24, grey);
				ms_at3_arrow(c, p, v3, 0.24, orange);
			}
			printf("  quarter turn: %s\n", sub);
			free(pick); free(cc); free(vv); free(ll); free(cen); free(fn); free(keep);
		}
		qaws_motion_destroy(mq);
		qaws_mesh_3d_destroy(box); qaws_mesh_3d_destroy(o);
	}
	{
		/* QUATERNION against LINEAR through subdivision: the identity on the left, minus a third of a turn on the right */
		static unsigned int const f4[4] = { 0, 1, 2, 3 }, s4[1] = { 4 };
		double const cs = cos(MS_PI / 3), sn = sin(MS_PI / 3);
		double q[16] = { 1, 0, 0, 0, -cs, 0, 0, -sn, -cs, 0, 0, -sn, 1, 0, 0, 0 }, pts[6 * 300];
		static double const strip_col[2][3] = { { 200, 214, 236 }, { 226, 214, 200 } };
		qaws_mesh_3d* sub_m[2] = { NULL, NULL };
		qaws_mesh_subdivide_desc sd;
		unsigned int npts = 0, ng[2] = { 0, 0 }, s;
		memset(&sd, 0, sizeof sd);
		sd.scheme = QAWS_MESH_SUBDIVIDE_CATMULL_CLARK;
		sd.levels = 4;
		for (s = 0; s < 2; s++)
		{
			double x[12] = { 0, 0, 0, 4, 0, 0, 4, 0.5, 0, 0, 0.5, 0 };
			qaws_mesh_3d* quad;
			for (i = 0; i < 4; i++) x[3 * i + 1] += s ? 0 : 1.3;
			quad = ms_make(x, 4, f4, s4, 1);
			ms_at_add(quad, "q", QAWS_MESH_VERTEX, 4, s ? QAWS_ATTRIBUTE_LINEAR : QAWS_ATTRIBUTE_QUATERNION, q, 4);
			qaws_mesh_3d_subdivide(&sub_m[s], quad, &sd);
			qaws_mesh_3d_destroy(quad);
			if (sub_m[s])
			{
				double* p = ms_positions(sub_m[s]);
				unsigned int nv = qaws_mesh_3d_get_vertex_count(sub_m[s]);
				for (i = 0; i < nv && npts < 300; i++)
					memcpy(pts + 3 * npts++, p + 3 * i, sizeof(double) * 3);
				free(p);
			}
		}
		pts[3 * npts] = 2; pts[3 * npts + 1] = 0.9; pts[3 * npts + 2] = 0.8;   /* (room above for the frames) */
		npts++;
		er_panel(c, 1200, 640 + 72, 600, 640 - 72, pts, npts, -0.25, 0.9, 0.06);
		er_text(c, 1200 + 14, 640 + 12, "QUATERNION, LINEAR subdivided", 3, 30, 30, 30);
		for (s = 0; s < 2; s++)
		{
			unsigned int el, nv;
			double *p, *qq;
			if (!sub_m[s])
				continue;
			ms_draw(c, sub_m[s], strip_col[s], NULL, 0, wire, 0.2);
			p = ms_positions(sub_m[s]);
			qq = ms_at_get(sub_m[s], "q", &el);
			nv = qaws_mesh_3d_get_vertex_count(sub_m[s]);
			for (i = 0; qq && i < nv; i++)
				if (fabs(p[3 * i + 1] - (s ? 0.25 : 1.55)) < 0.02)
				{
					ms_at3_triad(c, p + 3 * i, qq + 4 * i, 0.3);
					ng[s]++;
				}
			free(p); free(qq);
		}
		{
			char sub[200];
			sprintf(sub, "top QUATERNION, bottom LINEAR (%u, %u)", ng[0], ng[1]);
			er_text(c, 1200 + 14, 640 + 46, sub, 2, 95, 95, 95);
			printf("  subdivided: %s\n", sub);
		}
		qaws_mesh_3d_destroy(sub_m[0]); qaws_mesh_3d_destroy(sub_m[1]);
	}
	{
		/* frames through the cubic sweep */
		qaws_mesh_3d* o = NULL;
		if (m[1] && qaws_mesh_3d_sweep(&o, NULL, tool, m[1], NULL) == QAWS_STATUS_OK && o)
		{
			unsigned int nt = qaws_mesh_3d_get_face_count(o), np, el, *pick = (unsigned int*)malloc(sizeof(unsigned int) * (nt + 1));
			double *cc = ms_at_colours(o, &el), *qq = ms_at_get(o, "q", &el), *cen, *fn;
			unsigned char* keep = (unsigned char*)malloc(nt + 1);
			char sub[200];
			ms_at3_frames(o, &cen, &fn);
			for (i = 0; i < nt; i++) keep[i] = fn[3 * i + 2] > 0.99;
			np = ms_at3_pick(pick, cen, keep, nt, 0.45);
			sprintf(sub, "the tool's orientation (%u frames)", np);
			ms_panel(c, 1800, 640, 600, 640, o, -0.35, 0.95, "QUATERNION through a sweep", sub);
			ms_at_draw(c, o, cc, NULL, NULL, wire, 0.2);
			for (i = 0; qq && i < np; i++)
			{
				double p[3], q4[4];
				for (k = 0; k < 3; k++) p[k] = cen[3 * pick[i] + k] + (k == 2 ? 0.01 : 0);
				ms_at3_mean(q4, qq, 4, pick[i]);
				ms_at3_triad(c, p, q4, 0.3);
			}
			printf("  frames: %s\n", sub);
			free(pick); free(cc); free(qq); free(cen); free(fn); free(keep);
		}
		qaws_mesh_3d_destroy(o);
	}
	for (i = 0; i < 4; i++) qaws_motion_destroy(m[i]);
	qaws_curve_destroy(path);
	qaws_mesh_3d_destroy(tool);
	qaws_mesh_3d_destroy(block);
	er_save(c, "showcase/mesh63_sweep_schemes.png");
	printf("mesh63_sweep_schemes.png\n");
	er_destroy(c);
}
