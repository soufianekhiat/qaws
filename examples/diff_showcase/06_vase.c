/* ================================================================== */
/*  6. Vase from scanned points: fitting a revolution profile         */
/* ================================================================== */

#define VASE_CP 7
#define VASE_RINGS 22
#define VASE_SPOKES 14
#define VASE_N (VASE_RINGS * VASE_SPOKES)
#define VASE_ITERS 500

static double vase_radius(double z)
{
	return 0.85 + 0.38 * sin(2.1 * z + 0.3) - 0.08 * z;
}

static qaws_surface* vase_surface(qaws_curve* profile)
{
	qaws_surface_revolution_desc d;
	qaws_surface* s = NULL;
	memset(&d, 0, sizeof(d));
	d.profile = profile;
	d.axis_origin = v3(0, 0, 0);
	d.axis_direction = v3(0, 0, 1);
	d.angle = 0; /* full turn */
	qaws_surface_create_revolution(&d, &s);
	return s;
}

static double vase_height_value(void const* user, qaws_scalar u, qaws_scalar v)
{
	(void)user;
	(void)u;
	return v;
}

static void demo_vase(void)
{
	qaws_scalar cps[VASE_CP * 2], grad[VASE_CP * 2];
	qaws_scalar us[VASE_N], vs[VASE_N];
	qaws_vec3 targets[VASE_N];
	double loss[VASE_ITERS + 1];
	qaws_scalar initial[VASE_CP * 2];
	adam opt;
	int i, j, it;
	svg s;
	char buf[200];

	memset(&opt, 0, sizeof(opt));
	for (i = 0; i < VASE_CP; i++)
	{
		cps[2 * i] = 1.0f;                                  /* radius */
		cps[2 * i + 1] = (qaws_scalar)(3.0 * i / (VASE_CP - 1)); /* height */
	}
	memcpy(initial, cps, sizeof(cps));
	for (i = 0; i < VASE_RINGS; i++)
		for (j = 0; j < VASE_SPOKES; j++)
		{
			int k = i * VASE_SPOKES + j;
			double v = (i + 0.5) / VASE_RINGS, u = (j + 0.25 * (i % 2)) / VASE_SPOKES;
			double z = 3.0 * v, r = vase_radius(z) * (1 + 0.015 * sin(13.0 * k));
			us[k] = (qaws_scalar)u;
			vs[k] = (qaws_scalar)v;
			targets[k] = v3((qaws_scalar)(r * cos(2 * PI * u)), (qaws_scalar)(r * sin(2 * PI * u)), (qaws_scalar)z);
		}

	for (it = 0; it <= VASE_ITERS; it++)
	{
		qaws_curve* profile = bspline_2d(cps, VASE_CP);
		qaws_surface* vase = vase_surface(profile);
		static qaws_surface_jet primal[VASE_N], tangent[VASE_N], ybar[VASE_N];
		qaws_field_view fv;
		qaws_diff_views child, views;
		double l = 0;

		qaws_surface_eval_batch_tangent(NULL, vase, us, vs, NULL, NULL, VASE_N, QAWS_SJET_P, NULL, primal, tangent);
		memset(ybar, 0, sizeof(ybar));
		for (i = 0; i < VASE_N; i++)
		{
			qaws_vec3 d = v3(primal[i].d[0].x - targets[i].x, primal[i].d[0].y - targets[i].y, primal[i].d[0].z - targets[i].z);
			l += (d.x * d.x + d.y * d.y + d.z * d.z) / VASE_N;
			ybar[i].d[0] = v3((qaws_scalar)(2 * d.x / VASE_N), (qaws_scalar)(2 * d.y / VASE_N), (qaws_scalar)(2 * d.z / VASE_N));
		}
		loss[it] = l;

		/* The revolution has no own fields in this fit; its profile is
		   child 0, so the profile control point adjoints live there. */
		memset(grad, 0, sizeof(grad));
		child = one_field(&fv, QAWS_FIELD_CONTROL_POINTS, grad, VASE_CP, 2);
		views.fields = NULL;
		views.field_count = 0;
		views.children = &child;
		views.child_count = 1;
		qaws_surface_eval_batch_adjoint(NULL, vase, us, vs, VASE_N, QAWS_SJET_P, ybar, &views, NULL, NULL);
		qaws_surface_destroy(vase);
		qaws_curve_destroy(profile);
		if (it < VASE_ITERS)
			adam_step(&opt, cps, grad, VASE_CP * 2, 0.01);
	}

	svg_open(&s, "showcase/6_vase.svg", 1200, 600, "Vase from scanned points: fitting a surface of revolution",
		"308 noisy points; the gradient flows surface -> revolution rule -> child[0] (profile B-spline, 7 control points). Color = height parameter.");
	{
		viewport a = { 30, 80, 370, 500, 0, 0, 0, 0 };
		viewport b = { 415, 80, 370, 500, 0, 0, 0, 0 };
		viewport c = { 800, 80, 370, 280, -0.2, 1.8, -0.2, 3.2 };
		viewport lv = { 800, 380, 370, 200, 0, 0, 0, 0 };
		projection pa = { 215, 470, 95, 1.0 }, pb = { 600, 470, 95, 1.0 };
		qaws_curve* p0 = bspline_2d(initial, VASE_CP);
		qaws_curve* p1 = bspline_2d(cps, VASE_CP);
		qaws_surface* s0 = vase_surface(p0);
		qaws_surface* s1 = vase_surface(p1);
		double xy[2 * 200];

		svg_panel(&s, &a, "initial profile (cylinder) + scan points");
		draw_quads(&s, &pa, s0, 28, 16, vase_height_value, NULL, 0, 1);
		for (i = 0; i < VASE_N; i++)
		{
			double x, y;
			project(&pa, targets[i].x, targets[i].y, targets[i].z, &x, &y);
			svg_circle(&s, x, y, 1.5, "#24292f", "none");
		}
		svg_panel(&s, &b, "fitted vase");
		draw_quads(&s, &pb, s1, 28, 16, vase_height_value, NULL, 0, 1);

		/* Profile plot: (radius, height). */
		svg_panel(&s, &c, "profile r(z): dashed initial, blue fitted, green truth");
		for (i = 0; i < 200; i++)
		{
			double z = 3.0 * i / 199.0;
			xy[2 * i] = vx(&c, vase_radius(z));
			xy[2 * i + 1] = vy(&c, z);
		}
		svg_polyline(&s, xy, 200, "#2da44e", 5, 0.45, 0);
		curve_polyline(p0, &c, xy, 200);
		svg_polyline(&s, xy, 200, "#8c959f", 1.6, 1, 1);
		curve_polyline(p1, &c, xy, 200);
		svg_polyline(&s, xy, 200, "#0969da", 2.4, 1, 0);
		control_polygon(&s, &c, cps, VASE_CP, "#0969da");
		svg_loss_plot(&s, &lv, loss, VASE_ITERS + 1, "#cf222e", "point-to-surface MSE (log)");

		qaws_surface_destroy(s0);
		qaws_surface_destroy(s1);
		qaws_curve_destroy(p0);
		qaws_curve_destroy(p1);
	}
	svg_close(&s);
	sprintf(buf, "%.4e -> %.4e", loss[0], loss[VASE_ITERS]);
	printf("6_vase: MSE %s\n", buf);
}

