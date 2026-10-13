/* Figure 16: curves on surfaces traced by the surfaces' own geometry
   (qaws_surface_region.h). On a sphere, a wavy Bezier patch and a torus:
   geodesic polygons (each edge solved on the surface) and geodesic circles
   (every point at one geodesic distance from the centre). Each surface is
   drawn in 3D (back-facing parts faded) beside its parameter plane, where
   the same curves bend: they are traced intrinsically, not projected. */

typedef struct sgf_view
{
	double yaw, pitch, scale, cx, cy;
	double o[3];                /* the object's centre */
	double rot[9];
} sgf_view;

static void sgf_view_init(sgf_view* v, double yaw, double pitch, double scale, double cx, double cy)
{
	double cy_ = cos(yaw), sy = sin(yaw), cp = cos(pitch), sp = sin(pitch);
	v->yaw = yaw; v->pitch = pitch; v->scale = scale; v->cx = cx; v->cy = cy;
	v->o[0] = v->o[1] = v->o[2] = 0;
	/* yaw about z, then pitch about x: screen x, screen y (up), depth */
	v->rot[0] = cy_; v->rot[1] = -sy; v->rot[2] = 0;
	v->rot[3] = sy * sp; v->rot[4] = cy_ * sp; v->rot[5] = cp;
	v->rot[6] = sy * cp; v->rot[7] = cy_ * cp; v->rot[8] = -sp;
}

static void sgf_project(sgf_view const* v, qaws_vec3 p, double* sx, double* sy, double* depth)
{
	double x = p.x - v->o[0], y = p.y - v->o[1], z = p.z - v->o[2];
	*sx = v->cx + v->scale * (v->rot[0] * x + v->rot[1] * y + v->rot[2] * z);
	*sy = v->cy - v->scale * (v->rot[3] * x + v->rot[4] * y + v->rot[5] * z);
	if (depth) *depth = v->rot[6] * x + v->rot[7] * y + v->rot[8] * z;
}

/* facing the viewer at (u, v): the normal's component towards the eye */
static int sgf_facing(qaws_surface const* s, sgf_view const* v, double u, double w, double pu)
{
	qaws_surface_eval_result r;
	qaws_range ru = qaws_surface_get_u_range(s);
	double uu = u, nz;
	if (pu > 0) { uu = fmod(u - ru.min_value, pu); if (uu < 0) uu += pu; uu += ru.min_value; }
	if (qaws_surface_evaluate(s, (qaws_scalar)uu, (qaws_scalar)w, QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV | QAWS_SURFACE_EVAL_POSITION, &r) != QAWS_STATUS_OK)
		return 1;
	{
		double nx = (double)r.du.y * r.dv.z - (double)r.du.z * r.dv.y, ny = (double)r.du.z * r.dv.x - (double)r.du.x * r.dv.z;
		double nn = (double)r.du.x * r.dv.y - (double)r.du.y * r.dv.x;
		nz = -(v->rot[6] * nx + v->rot[7] * ny + v->rot[8] * nn);
	}
	return nz >= 0;
}

/* a polyline of surface points, split where it turns away (faded there) */
static void sgf_curve3d(svg* s, qaws_surface const* surf, sgf_view const* v, double const* uv, unsigned int n, double pu,
	char const* color, double width, int two_sided)
{
	double* xy = (double*)malloc(sizeof(double) * 2 * (n + 1));
	unsigned int k, m = 0;
	int prev = -1;
	qaws_range ru = qaws_surface_get_u_range(surf);
	for (k = 0; k < n; k++)
	{
		qaws_surface_eval_result r;
		double uu = uv[2 * k], sx, sy;
		int f = two_sided ? 1 : sgf_facing(surf, v, uv[2 * k], uv[2 * k + 1], pu);
		if (pu > 0) { uu = fmod(uu - ru.min_value, pu); if (uu < 0) uu += pu; uu += ru.min_value; }
		qaws_surface_evaluate(surf, (qaws_scalar)uu, (qaws_scalar)uv[2 * k + 1], QAWS_SURFACE_EVAL_POSITION, &r);
		sgf_project(v, r.position, &sx, &sy, NULL);
		if (prev >= 0 && f != prev)
		{
			xy[2 * m] = sx; xy[2 * m + 1] = sy; m++;
			if (m > 1) svg_polyline(s, xy, (int)m, color, prev ? width : 0.7 * width, prev ? 1.0 : 0.25, !prev);
			m = 0;
		}
		xy[2 * m] = sx; xy[2 * m + 1] = sy; m++;
		prev = f;
	}
	if (m > 1) svg_polyline(s, xy, (int)m, color, prev ? width : 0.7 * width, prev ? 1.0 : 0.25, !prev);
	free(xy);
}

/* the surface's parameter lines */
static void sgf_mesh(svg* s, qaws_surface const* surf, sgf_view const* v, unsigned int lines, double pu)
{
	qaws_range ru = qaws_surface_get_u_range(surf), rv = qaws_surface_get_v_range(surf);
	unsigned int i, k;
	double uv[2 * 161];
	for (i = 0; i <= lines; i++)
	{
		for (k = 0; k <= 160; k++)
		{
			uv[2 * k] = ru.min_value + (ru.max_value - ru.min_value) * i / lines;
			uv[2 * k + 1] = rv.min_value + (rv.max_value - rv.min_value) * k / 160;
		}
		sgf_curve3d(s, surf, v, uv, 161, pu, "#afb8c1", 0.7, 0);
		for (k = 0; k <= 160; k++)
		{
			uv[2 * k] = ru.min_value + (ru.max_value - ru.min_value) * k / 160;
			uv[2 * k + 1] = rv.min_value + (rv.max_value - rv.min_value) * i / lines;
		}
		sgf_curve3d(s, surf, v, uv, 161, pu, "#afb8c1", 0.7, 0);
	}
}

/* a surface curve sampled in (u, v) */
static unsigned int sgf_sample(qaws_surface_curve const* c, double* uv, unsigned int n)
{
	qaws_range rg = qaws_curve_get_parameter_range(qaws_surface_curve_get_uv(c));
	unsigned int k;
	for (k = 0; k < n; k++)
	{
		qaws_vec2 p;
		qaws_surface_curve_evaluate(c, (qaws_scalar)(rg.min_value + (rg.max_value - rg.min_value) * k / (n - 1)), &p, NULL);
		uv[2 * k] = p.x; uv[2 * k + 1] = p.y;
	}
	return n;
}

typedef struct sgf_item
{
	qaws_surface_curve* c;
	char const* color;
} sgf_item;

static void sgf_panel(svg* s, viewport const* v3, viewport const* vuv, qaws_surface const* surf, sgf_view const* view,
	sgf_item const* items, unsigned int n, double pu, char const* label3, char const* labeluv)
{
	unsigned int i;
	double uv[2 * 400];
	qaws_range ru = qaws_surface_get_u_range(surf), rv = qaws_surface_get_v_range(surf);
	svg_panel(s, v3, label3);
	svg_panel(s, vuv, labeluv);
	sgf_mesh(s, surf, view, 12, pu);
	/* the parameter plane's domain */
	{
		double box[10];
		box[0] = vx(vuv, ru.min_value); box[1] = vy(vuv, rv.min_value);
		box[2] = vx(vuv, ru.max_value); box[3] = vy(vuv, rv.min_value);
		box[4] = vx(vuv, ru.max_value); box[5] = vy(vuv, rv.max_value);
		box[6] = vx(vuv, ru.min_value); box[7] = vy(vuv, rv.max_value);
		box[8] = box[0]; box[9] = box[1];
		svg_polyline(s, box, 5, "#57606a", 1.0, 1.0, 0);
	}
	for (i = 0; i < n; i++)
	{
		unsigned int m = sgf_sample(items[i].c, uv, 400), k;
		double xy[2 * 400];
		sgf_curve3d(s, surf, view, uv, m, pu, items[i].color, 2.0, 0);
		for (k = 0; k < m; k++) { xy[2 * k] = vx(vuv, uv[2 * k]); xy[2 * k + 1] = vy(vuv, uv[2 * k + 1]); }
		svg_polyline(s, xy, (int)m, items[i].color, 1.8, 1.0, 0);
	}
}

static qaws_surface* sgf_revolve(qaws_curve* profile)
{
	qaws_surface_revolution_desc rd;
	qaws_surface* s = NULL;
	memset(&rd, 0, sizeof(rd));
	rd.profile = profile; rd.axis_direction.z = 1;
	qaws_surface_create_revolution(&rd, &s);
	return s;
}

static qaws_curve* sgf_arc(double cx, double cy, double r, double a0, double a1)
{
	qaws_arc_segment g;
	qaws_arc_desc ad;
	qaws_curve* c = NULL;
	memset(&g, 0, sizeof(g));
	g.center[0] = (qaws_scalar)cx; g.center[1] = (qaws_scalar)cy; g.radius = (qaws_scalar)r;
	g.angle_start = (qaws_scalar)a0; g.angle_end = (qaws_scalar)a1;
	memset(&ad, 0, sizeof(ad));
	ad.dimension = QAWS_DIMENSION_2D; ad.segments = &g; ad.segment_count = 1;
	qaws_curve_create_arc(&ad, &c);
	return c;
}

static qaws_vec2 sgf_uv(double u, double v)
{
	qaws_vec2 p;
	p.x = (qaws_scalar)u; p.y = (qaws_scalar)v;
	return p;
}

static void demo_surface_geodesic(void)
{
	clock_t t0 = clock();
	svg s;
	double const pw = 420, ph = 330, pi = 3.14159265358979323846;
	static char const* const ring[5] = { "#1f6feb", "#2da44e", "#bf8700", "#cf222e", "#8250df" };
	unsigned int row, i;
	if (!svg_open(&s, "showcase/b2d16_surface_geodesic.svg", (int)(48 + 2 * (pw + 12) + pw + 12), (int)(110 + 3 * (ph + 44)),
		"Curves on surfaces, traced by the surfaces' own geometry",
		"qaws_surface_region: geodesic polygons (edges solved on the surface) and geodesic circles (one geodesic distance from the centre); right: their parameter-plane images"))
		return;
	for (row = 0; row < 3; row++)
	{
		qaws_curve* prof = NULL;
		qaws_surface* surf = NULL;
		qaws_vec3 cp[16];
		sgf_item items[16];
		unsigned int n = 0;
		double pu = 0;
		sgf_view view;
		viewport v3, vuv;
		char const* name = "";
		qaws_range ru, rv;
		if (row == 0)
		{
			prof = sgf_arc(0, 0, 1, -pi / 2, pi / 2);
			surf = sgf_revolve(prof);
			pu = 1;
			name = "a sphere: a geodesic triangle and circles of radius 0.1 to 0.5 about one centre";
		}
		else if (row == 1)
		{
			static double const z[16] = { 0, 0.3, -0.2, 0, 0.35, 0.8, 0.5, 0.1, -0.3, 0.6, 0.9, 0.35, 0, -0.15, 0.3, 0 };
			qaws_surface_bezier_desc bd;
			for (i = 0; i < 16; i++) { cp[i].x = (qaws_scalar)((i / 4) / 3.0 * 2); cp[i].y = (qaws_scalar)((i % 4) / 3.0 * 2); cp[i].z = (qaws_scalar)z[i]; }
			memset(&bd, 0, sizeof(bd));
			bd.u_degree = 3; bd.v_degree = 3; bd.control_points = cp; bd.u_point_count = 4; bd.v_point_count = 4;
			qaws_surface_create_bezier(&bd, &surf);
			name = "a wavy bicubic patch: a geodesic square and circles about its centre";
		}
		else
		{
			prof = sgf_arc(1.0, 0, 0.42, 0, 2 * pi);
			surf = sgf_revolve(prof);
			pu = 1;
			name = "a torus: a geodesic quadrilateral round its outer equator and circles on top of the tube";
		}
		ru = qaws_surface_get_u_range(surf); rv = qaws_surface_get_v_range(surf);
		if (row == 0)
		{
			qaws_vec2 a = sgf_uv(0.14, 0.68), b = sgf_uv(0.3, 0.74), c = sgf_uv(0.24, 0.5), ctr = sgf_uv(0.4, 0.42);
			qaws_surface_curve_create_geodesic(surf, a, b, 0, &items[n].c); items[n++].color = "#24292f";
			qaws_surface_curve_create_geodesic(surf, b, c, 0, &items[n].c); items[n++].color = "#24292f";
			qaws_surface_curve_create_geodesic(surf, c, a, 0, &items[n].c); items[n++].color = "#24292f";
			for (i = 0; i < 5; i++)
				if (qaws_surface_curve_create_geodesic_arc(surf, ctr, (qaws_scalar)(0.1 * (i + 1)), 0, (qaws_scalar)(2 * pi), 0, &items[n].c) == QAWS_STATUS_OK)
					items[n++].color = ring[i];
			sgf_view_init(&view, 3.0, 0.3, 125, 0, 0);
		}
		else if (row == 1)
		{
			qaws_vec2 q[4];
			q[0] = sgf_uv(0.15, 0.2); q[1] = sgf_uv(0.85, 0.15); q[2] = sgf_uv(0.8, 0.85); q[3] = sgf_uv(0.2, 0.8);
			for (i = 0; i < 4; i++)
				if (qaws_surface_curve_create_geodesic(surf, q[i], q[(i + 1) % 4], 0, &items[n].c) == QAWS_STATUS_OK)
					items[n++].color = "#24292f";
			for (i = 0; i < 4; i++)
				if (qaws_surface_curve_create_geodesic_arc(surf, sgf_uv(0.5, 0.5), (qaws_scalar)(0.15 + 0.15 * i), 0, (qaws_scalar)(2 * pi), 0, &items[n].c) == QAWS_STATUS_OK)
					items[n++].color = ring[i];
			sgf_view_init(&view, 0.7, 0.75, 150, 0, 0);
			view.o[0] = 1; view.o[1] = 1; view.o[2] = 0.3;
		}
		else
		{
			qaws_vec2 q[4];
			q[0] = sgf_uv(0.0, 0.0); q[1] = sgf_uv(0.25, 0.08); q[2] = sgf_uv(0.5, 0.0); q[3] = sgf_uv(0.75, -0.08);
			for (i = 0; i < 4; i++)
			{
				qaws_vec2 b = q[(i + 1) % 4];
				if (i == 3) b.x += 1;      /* back round to the start, in the covering plane */
				if (qaws_surface_curve_create_geodesic(surf, q[i], b, 0, &items[n].c) == QAWS_STATUS_OK)
					items[n++].color = "#24292f";
			}
			for (i = 0; i < 3; i++)
				if (qaws_surface_curve_create_geodesic_arc(surf, sgf_uv(0.62, 0.25), (qaws_scalar)(0.12 + 0.12 * i), 0, (qaws_scalar)(2 * pi), 0, &items[n].c) == QAWS_STATUS_OK)
					items[n++].color = ring[i];
			sgf_view_init(&view, 0.81, 0.85, 105, 0, 0);
		}
		v3.x0 = 24; v3.y0 = 104 + row * (ph + 44); v3.w = 2 * pw + 12; v3.h = ph;
		v3.xmin = 0; v3.xmax = 1; v3.ymin = 0; v3.ymax = 1;
		view.cx = v3.x0 + v3.w / 2;
		view.cy = v3.y0 + v3.h / 2 + 8;
		vuv.x0 = 24 + 2 * (pw + 12); vuv.y0 = v3.y0; vuv.w = pw; vuv.h = ph;
		{
			double padu = 0.06 * (ru.max_value - ru.min_value), padv = 0.08 * (rv.max_value - rv.min_value);
			vuv.xmin = ru.min_value - padu; vuv.xmax = ru.max_value + padu;
			vuv.ymin = rv.min_value - padv; vuv.ymax = rv.max_value + padv;
			if (row == 2) { vuv.ymin = rv.min_value - padv - 0.15; vuv.ymax = rv.max_value + padv; }
		}
		sgf_panel(&s, &v3, &vuv, surf, &view, items, n, pu, name, "the parameter plane (u, v)");
		for (i = 0; i < n; i++) qaws_surface_curve_destroy(items[i].c);
		qaws_surface_destroy(surf);
		qaws_curve_destroy(prof);
	}
	svg_close(&s);
	printf("figure 16: curves on surfaces, %.1f s\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
}
