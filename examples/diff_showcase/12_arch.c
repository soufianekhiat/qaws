/* ================================================================== */
/*  12. Arch design from implicit features: intersections with a deck */
/*      and the apex (extremum) placed by adjoints of their solves     */
/* ================================================================== */

#define AR_N 8
#define AR_ITERS 160
#define AR_DECK 1.0
#define AR_APEX_X 3.0
#define AR_APEX_Y 2.4
#define AR_X0 1.1
#define AR_X1 4.9

static qaws_curve* arch_curve(qaws_scalar const* cps)
{
	return bspline_2d(cps, AR_N);
}

static qaws_curve* arch_deck(void)
{
	static qaws_scalar const d[8] = { -0.5, AR_DECK, 2.0, AR_DECK, 4.0, AR_DECK, 6.5, AR_DECK };
	return bspline_2d(d, 4);
}

typedef struct arch_features
{
	int ok;
	qaws_scalar t_apex, t_cross[2], tb_cross[2];
	qaws_vec3 apex, cross[2];
} arch_features;

/* Finds the features with the discrete finders (the active set), then
   accumulates the gradient of the design loss through the adjoints of
   the implicit solves. Returns the loss. */
static double arch_loss(qaws_curve const* arch, qaws_curve const* deck, qaws_scalar* grad, arch_features* f)
{
	qaws_intersection_2d hits[4];
	qaws_scalar ext[8];
	unsigned int nh = 0, ne = 0, i;
	double loss = 0;
	qaws_field_view fv;
	qaws_diff_views gv = one_field(&fv, QAWS_FIELD_CONTROL_POINTS, grad, AR_N, 2);
	qaws_vec3 up = v3(0, 1, 0);

	memset(f, 0, sizeof(*f));
	if (grad)
		memset(grad, 0, sizeof(qaws_scalar) * AR_N * 2);
	qaws_curve_find_intersections_2d(arch, deck, hits, 4, &nh);
	qaws_curve_find_extrema(arch, 1, ext, 8, &ne);
	if (nh < 2 || ne < 1)
		return 1e9;
	/* apex: the highest extremum; crossings: left-most and right-most */
	{
		qaws_scalar best = -1e30f;
		for (i = 0; i < ne; i++)
		{
			qaws_eval_result_2d r;
			qaws_curve_evaluate_2d(arch, ext[i], QAWS_EVAL_FLAG_POSITION, &r);
			if (r.position.y > best)
			{
				best = r.position.y;
				f->t_apex = ext[i];
			}
		}
		f->t_cross[0] = f->t_cross[1] = hits[0].parameter_a;
		f->tb_cross[0] = f->tb_cross[1] = hits[0].parameter_b;
		for (i = 1; i < nh; i++)
		{
			if (hits[i].parameter_a < f->t_cross[0]) { f->t_cross[0] = hits[i].parameter_a; f->tb_cross[0] = hits[i].parameter_b; }
			if (hits[i].parameter_a > f->t_cross[1]) { f->t_cross[1] = hits[i].parameter_a; f->tb_cross[1] = hits[i].parameter_b; }
		}
	}
	{
		qaws_curve_extremum e, eb;
		qaws_curve_extremum_tangent(NULL, arch, up, f->t_apex, NULL, NULL, &e, NULL);
		f->apex = e.position;
		loss += (e.position.x - AR_APEX_X) * (e.position.x - AR_APEX_X) + (e.value - AR_APEX_Y) * (e.value - AR_APEX_Y);
		if (grad)
		{
			memset(&eb, 0, sizeof(eb));
			eb.position.x = (qaws_scalar)(2 * (e.position.x - AR_APEX_X));
			eb.value = (qaws_scalar)(2 * (e.value - AR_APEX_Y));
			qaws_curve_extremum_adjoint(NULL, arch, up, f->t_apex, &eb, &gv, NULL);
		}
	}
	for (i = 0; i < 2; i++)
	{
		qaws_curve_pair_point p, pb;
		double target = i == 0 ? AR_X0 : AR_X1;
		qaws_curve_pair_point_tangent(NULL, arch, deck, f->t_cross[i], f->tb_cross[i], NULL, NULL, &p, NULL);
		f->cross[i] = p.position_a;
		loss += (p.position_a.x - target) * (p.position_a.x - target);
		if (grad)
		{
			memset(&pb, 0, sizeof(pb));
			pb.position_a.x = (qaws_scalar)(2 * (p.position_a.x - target));
			qaws_curve_pair_point_adjoint(NULL, arch, deck, f->t_cross[i], f->tb_cross[i], &pb, &gv, NULL);
		}
	}
	/* a little fairness */
	{
		qaws_scalar bg[AR_N * 2], bend = 0;
		qaws_field_view bf;
		qaws_diff_views bv = one_field(&bf, QAWS_FIELD_CONTROL_POINTS, bg, AR_N, 2);
		memset(bg, 0, sizeof(bg));
		qaws_curve_functional_gradient(NULL, arch, QAWS_FUNCTIONAL_BENDING, 0, &bv, &bend);
		loss += 1e-3 * bend;
		if (grad)
			for (i = 0; i < AR_N * 2; i++)
				grad[i] += (qaws_scalar)(1e-3 * bg[i]);
	}
	f->ok = 1;
	return loss;
}

static void demo_arch(void)
{
	qaws_scalar cps[AR_N * 2], g[AR_N * 2], init[AR_N * 2];
	qaws_scalar ghosts[12][AR_N * 2];
	arch_features f0, f1;
	double loss[AR_ITERS + 1];
	qaws_curve* deck = arch_deck();
	adam opt;
	int it, i, ng = 0;
	svg s;
	char buf[200];
	static double const start[AR_N][2] = {
		{ 0.0, 0.0 }, { 0.3, 0.9 }, { 0.9, 1.9 }, { 2.4, 1.6 }, { 3.6, 3.0 }, { 4.8, 2.2 }, { 5.6, 1.0 }, { 6.0, 0.0 } };

	for (i = 0; i < AR_N; i++)
	{
		cps[2 * i] = (qaws_scalar)start[i][0];
		cps[2 * i + 1] = (qaws_scalar)start[i][1];
	}
	memcpy(init, cps, sizeof(cps));
	memset(&opt, 0, sizeof(opt));
	for (it = 0; it <= AR_ITERS; it++)
	{
		qaws_curve* c = arch_curve(cps);
		arch_features f;
		loss[it] = arch_loss(c, deck, g, &f);
		qaws_curve_destroy(c);
		if (it % 16 == 0 && ng < 12)
			memcpy(ghosts[ng++], cps, sizeof(cps));
		if (it == AR_ITERS)
			break;
		/* feet stay on the ground */
		g[0] = g[1] = g[2 * AR_N - 2] = g[2 * AR_N - 1] = 0;
		adam_step(&opt, cps, g, AR_N * 2, 0.04 * (1.0 - 0.8 * it / (double)AR_ITERS));
	}

	svg_open(&s, "showcase/12_arch_features.svg", 1200, 600, "Arch design from implicit features",
		"Deck crossings (curve pair solves) and the apex (extremum solve) are steered to targets by adjoints of the implicit "
		"function theorem; the finders only pick which solutions exist.");
	{
		viewport v = { 24, 80, 760, 500, -0.6, 6.6, -0.5, 3.6 };
		viewport lv = { 800, 80, 376, 300, 0, 0, 0, 0 };
		double xy[2 * 300];
		qaws_curve* c;

		svg_panel(&s, &v, "initial (grey), iterations (light blue), final (blue)");
		/* ground and deck */
		svg_line(&s, vx(&v, -0.6), vy(&v, 0), vx(&v, 6.6), vy(&v, 0), "#8c959f", 2, 1);
		svg_line(&s, vx(&v, -0.5), vy(&v, AR_DECK), vx(&v, 6.5), vy(&v, AR_DECK), "#57606a", 5, 0.85);
		svg_text(&s, vx(&v, 6.45), vy(&v, AR_DECK) - 8, 12, "#57606a", "end", "deck");
		for (i = 1; i < ng - 1; i++)
		{
			c = arch_curve(ghosts[i]);
			curve_polyline(c, &v, xy, 300);
			svg_polyline(&s, xy, 300, "#54aeff", 1.4, 0.45, 0);
			qaws_curve_destroy(c);
		}
		c = arch_curve(init);
		arch_loss(c, deck, NULL, &f0);
		curve_polyline(c, &v, xy, 300);
		svg_polyline(&s, xy, 300, "#8c959f", 2, 1, 1);
		qaws_curve_destroy(c);
		c = arch_curve(cps);
		arch_loss(c, deck, NULL, &f1);
		curve_polyline(c, &v, xy, 300);
		svg_polyline(&s, xy, 300, "#0969da", 3.2, 1, 0);
		control_polygon(&s, &v, cps, AR_N, "#0969da");
		qaws_curve_destroy(c);

		/* targets */
		svg_line(&s, vx(&v, AR_X0), vy(&v, AR_DECK) - 18, vx(&v, AR_X0), vy(&v, AR_DECK) + 18, "#cf222e", 2, 0.9);
		svg_line(&s, vx(&v, AR_X1), vy(&v, AR_DECK) - 18, vx(&v, AR_X1), vy(&v, AR_DECK) + 18, "#cf222e", 2, 0.9);
		svg_line(&s, vx(&v, AR_APEX_X) - 12, vy(&v, AR_APEX_Y), vx(&v, AR_APEX_X) + 12, vy(&v, AR_APEX_Y), "#cf222e", 2, 0.9);
		svg_line(&s, vx(&v, AR_APEX_X), vy(&v, AR_APEX_Y) - 12, vx(&v, AR_APEX_X), vy(&v, AR_APEX_Y) + 12, "#cf222e", 2, 0.9);
		svg_text(&s, vx(&v, AR_APEX_X) + 16, vy(&v, AR_APEX_Y) - 8, 12, "#cf222e", "start", "apex target");
		svg_text(&s, vx(&v, AR_X0), vy(&v, AR_DECK) + 34, 12, "#cf222e", "middle", "crossing target");
		svg_text(&s, vx(&v, AR_X1), vy(&v, AR_DECK) + 34, 12, "#cf222e", "middle", "crossing target");
		/* features: initial (grey) and final (blue) */
		svg_circle(&s, vx(&v, f0.apex.x), vy(&v, f0.apex.y), 5, "#ffffff", "#57606a");
		svg_circle(&s, vx(&v, f0.cross[0].x), vy(&v, f0.cross[0].y), 5, "#ffffff", "#57606a");
		svg_circle(&s, vx(&v, f0.cross[1].x), vy(&v, f0.cross[1].y), 5, "#ffffff", "#57606a");
		svg_circle(&s, vx(&v, f1.apex.x), vy(&v, f1.apex.y), 6, "#0969da", "#ffffff");
		svg_circle(&s, vx(&v, f1.cross[0].x), vy(&v, f1.cross[0].y), 6, "#0969da", "#ffffff");
		svg_circle(&s, vx(&v, f1.cross[1].x), vy(&v, f1.cross[1].y), 6, "#0969da", "#ffffff");

		svg_loss_plot(&s, &lv, loss, AR_ITERS + 1, "#0969da", "design loss (log scale)");
		sprintf(buf, "apex (%.3f, %.3f)   target (%.1f, %.1f)", f1.apex.x, f1.apex.y, AR_APEX_X, AR_APEX_Y);
		svg_text(&s, lv.x0, lv.y0 + lv.h + 26, 13, "#24292f", "start", buf);
		sprintf(buf, "crossings x = %.3f, %.3f   targets %.1f, %.1f", f1.cross[0].x, f1.cross[1].x, AR_X0, AR_X1);
		svg_text(&s, lv.x0, lv.y0 + lv.h + 48, 13, "#24292f", "start", buf);
		svg_text(&s, lv.x0, lv.y0 + lv.h + 76, 12, "#57606a", "start", "qaws_curve_pair_point_adjoint (deck crossings)");
		svg_text(&s, lv.x0, lv.y0 + lv.h + 94, 12, "#57606a", "start", "qaws_curve_extremum_adjoint (apex height, position)");
		svg_text(&s, lv.x0, lv.y0 + lv.h + 112, 12, "#57606a", "start", "+ 1e-3 bending (functional gradient), feet masked");
	}
	svg_close(&s);
	qaws_curve_destroy(deck);
	printf("12_arch_features: loss %.4e -> %.4e, apex (%.4f, %.4f), crossings %.4f %.4f\n", loss[0], loss[AR_ITERS],
		f1.apex.x, f1.apex.y, f1.cross[0].x, f1.cross[1].x);
}


