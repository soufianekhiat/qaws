/* ================================================================== */
/*  9. Projection-based fitting and the validity map of closest point */
/* ================================================================== */

#define PRJ_CP 10
#define PRJ_N 70
#define PRJ_ITERS 400

static void prj_target(int i, double* x, double* y)
{
	double s = i / (double)(PRJ_N - 1);
	*x = 0.7 + 4.4 * s;
	*y = 2.1 + 1.25 * sin(2 * PI * 1.25 * s) + 0.5 * s;
}

static void demo_projection(void)
{
	qaws_scalar cps[PRJ_CP * 2], grad[PRJ_CP * 2];
	qaws_vec3 targets[PRJ_N];
	double loss[PRJ_ITERS + 1];
	adam opt;
	int i, it;
	svg s;
	viewport a = { 30, 80, 560, 500, 0.2, 5.8, -0.3, 4.6 };
	viewport b = { 610, 80, 560, 500, 0.2, 5.8, -0.3, 4.6 };
	viewport lv = { 360, 470, 220, 100, 0, 0, 0, 0 };
	char buf[160];

	memset(&opt, 0, sizeof(opt));
	for (i = 0; i < PRJ_N; i++)
	{
		double x, y;
		prj_target(i, &x, &y);
		targets[i] = v3((qaws_scalar)x, (qaws_scalar)y, 0);
	}
	for (i = 0; i < PRJ_CP; i++)
	{
		cps[2 * i] = (qaws_scalar)(1.0 + 3.6 * i / (PRJ_CP - 1));
		cps[2 * i + 1] = (qaws_scalar)2.1;
	}

	for (it = 0; it <= PRJ_ITERS; it++)
	{
		qaws_curve* c = bspline_2d(cps, PRJ_CP);
		qaws_field_view fv;
		qaws_diff_views views;
		double l = 0;
		memset(grad, 0, sizeof(grad));
		views = one_field(&fv, QAWS_FIELD_CONTROL_POINTS, grad, PRJ_CP, 2);
		for (i = 0; i < PRJ_N; i++)
		{
			/* d^2 to the curve: no correspondences, the foot point is
			   re-solved and differentiated through its optimality. */
			qaws_curve_closest_point val, bar;
			qaws_curve_closest_point_tangent(NULL, c, targets[i], NULL, NULL, &val, NULL);
			l += val.distance * val.distance / PRJ_N;
			memset(&bar, 0, sizeof(bar));
			bar.distance = (qaws_scalar)(2 * val.distance / PRJ_N);
			qaws_curve_closest_point_adjoint(NULL, c, targets[i], &bar, &views, NULL);
		}
		loss[it] = l;
		qaws_curve_destroy(c);
		if (it < PRJ_ITERS)
			adam_step(&opt, cps, grad, PRJ_CP * 2, 0.02);
	}

	svg_open(&s, "showcase/9_projection.svg", 1200, 600, "Closest point through the implicit function theorem",
		"Left: fit without correspondences, minimizing squared point-to-curve distance (closest-point adjoint). Right: validity reported for a grid of queries.");

	svg_panel(&s, &a, "fit by distance to the curve");
	{
		qaws_curve* c = bspline_2d(cps, PRJ_CP);
		double xy[2 * 300];
		curve_polyline(c, &a, xy, 300);
		for (i = 0; i < PRJ_N; i++)
		{
			qaws_curve_closest_point val;
			qaws_curve_closest_point_tangent(NULL, c, targets[i], NULL, NULL, &val, NULL);
			svg_line(&s, vx(&a, targets[i].x), vy(&a, targets[i].y), vx(&a, val.position.x), vy(&a, val.position.y), "#cf222e", 1, 0.7);
			svg_circle(&s, vx(&a, targets[i].x), vy(&a, targets[i].y), 2.6, "#d0d7de", "#57606a");
		}
		svg_polyline(&s, xy, 300, "#0969da", 2.6, 1, 0);
		control_polygon(&s, &a, cps, PRJ_CP, "#0969da");
		{
			qaws_scalar init[PRJ_CP * 2];
			qaws_curve* c0;
			for (i = 0; i < PRJ_CP; i++)
			{
				init[2 * i] = (qaws_scalar)(1.0 + 3.6 * i / (PRJ_CP - 1));
				init[2 * i + 1] = (qaws_scalar)2.1;
			}
			c0 = bspline_2d(init, PRJ_CP);
			curve_polyline(c0, &a, xy, 300);
			svg_polyline(&s, xy, 300, "#8c959f", 1.6, 1, 1);
			qaws_curve_destroy(c0);
		}
		svg_loss_plot(&s, &lv, loss, PRJ_ITERS + 1, "#cf222e", "mean d^2");
		qaws_curve_destroy(c);
	}

	/* Validity map: the medial axis shows up as the ambiguous region. */
	svg_panel(&s, &b, "validity of d(foot point)/d(query) over the plane");
	{
		qaws_curve* c = bspline_2d(cps, PRJ_CP);
		double xy[2 * 300];
		int gx, gy, nx = 112, ny = 100, counts[3] = { 0, 0, 0 };
		double cw = b.w / nx, ch = b.h / ny;
		for (gx = 0; gx < nx; gx++)
			for (gy = 0; gy < ny; gy++)
			{
				qaws_diff_context ctx;
				qaws_diff_report report;
				qaws_curve_closest_point val;
				double wx = b.xmin + (b.xmax - b.xmin) * (gx + 0.5) / nx;
				double wy = b.ymax - (b.ymax - b.ymin) * (gy + 0.5) / ny;
				char col[40];
				double rel;
				qaws_diff_context_init(&ctx);
				qaws_diff_report_reset(&report);
				ctx.report = &report;
				qaws_curve_closest_point_tangent(&ctx, c, v3((qaws_scalar)wx, (qaws_scalar)wy, 0), NULL, NULL, &val, NULL);
				/* Shade by the reported gap to the best competing foot point,
				   relative to the distance: the medial axis is where it vanishes. */
				rel = report.branch_gap / (val.distance + 1e-9);
				if (rel > 0.3) rel = 0.3;
				rel /= 0.3;
				if (report.validity == QAWS_DIFF_AMBIGUOUS || report.validity == QAWS_DIFF_ILL_CONDITIONED)
				{
					sprintf(col, "#cf222e");
					counts[1]++;
				}
				else if (report.validity == QAWS_DIFF_VALID_LOCALLY)
				{
					sprintf(col, "rgb(%d,%d,%d)", (int)(150 + 70 * rel), (int)(200 + 40 * rel), 255);
					counts[2]++;
				}
				else
				{
					sprintf(col, "rgb(%d,%d,%d)", (int)(235 - 17 * rel), (int)(120 + 131 * rel), (int)(110 + 115 * rel));
					counts[0]++;
				}
				fprintf(s.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.2f\" height=\"%.2f\" fill=\"%s\"/>\n",
					b.x0 + gx * cw, b.y0 + gy * ch, cw + 0.3, ch + 0.3, col);
			}
		curve_polyline(c, &b, xy, 300);
		svg_polyline(&s, xy, 300, "#0969da", 2.6, 1, 0);
		fprintf(s.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"330\" height=\"64\" fill=\"#ffffff\" stroke=\"#d0d7de\"/>\n",
			b.x0 + 10, b.y0 + b.h - 74);
		fprintf(s.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"12\" height=\"12\" fill=\"#dafbe1\" stroke=\"#8c959f\"/>\n", b.x0 + 20, b.y0 + b.h - 66);
		svg_text(&s, b.x0 + 38, b.y0 + b.h - 56, 12, "#24292f", "start", "valid; pinker = smaller gap to a competing foot point");
		fprintf(s.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"12\" height=\"12\" fill=\"#b6e3ff\" stroke=\"#8c959f\"/>\n", b.x0 + 20, b.y0 + b.h - 48);
		svg_text(&s, b.x0 + 38, b.y0 + b.h - 38, 12, "#24292f", "start", "valid locally: endpoint (active set frozen)");
		fprintf(s.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"12\" height=\"12\" fill=\"#cf222e\"/>\n", b.x0 + 20, b.y0 + b.h - 30);
		svg_text(&s, b.x0 + 38, b.y0 + b.h - 20, 12, "#24292f", "start", "ambiguous: gap &lt; 0.5% of the distance (medial axis)");
		qaws_curve_destroy(c);
		sprintf(buf, "%d valid, %d ambiguous, %d endpoint", counts[0], counts[1], counts[2]);
		printf("9_projection: fit %.4e -> %.4e, map: %s\n", loss[0], loss[PRJ_ITERS], buf);
	}
	svg_close(&s);
}

