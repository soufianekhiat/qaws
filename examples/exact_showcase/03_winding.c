/* ================================================================== */
/*  3. Certified point-in-region on a loop of rational conics         */
/* ================================================================== */

#define WD_N 160
#define WD_R 1048576.0   /* 2^20 lattice units */

/* Sampled winding: angle sum over a polyline of `per` samples per piece
   (the approach of qaws_curve_compute_winding_number_2d, 256 in total). */
static int wd_sampled(qaws_curve* const* pieces, int count, int per, double const* p)
{
	double sum = 0, px = 0, py = 0;
	int k, i, first = 1;
	for (k = 0; k < count; k++)
		for (i = (k == 0 ? 0 : 1); i <= per; i++)
		{
			qaws_eval_result_2d r;
			double dx, dy;
			qaws_curve_evaluate_2d(pieces[k], (qaws_scalar)((double)i / per), QAWS_EVAL_FLAG_POSITION, &r);
			dx = r.position.x - p[0];
			dy = r.position.y - p[1];
			if (!first)
				sum += atan2(px * dy - py * dx, px * dx + py * dy);
			px = dx;
			py = dy;
			first = 0;
		}
	return (int)floor(sum / (2 * 3.14159265358979) + 0.5);
}

static void demo_winding(void)
{
	static signed char sampled[WD_N * WD_N], certified[WD_N * WD_N], wrong[WD_N * WD_N];
	static char const* const in_colors[3] = { "#cf222e", "#f6f8fa", "#0969da" };   /* undecided, outside, inside */
	static char const* const wrong_colors[3] = { "#f6f8fa", "#f6f8fa", "#cf222e" };
	qaws_curve* pieces[4];
	qaws_exact_curve* exact[4];
	qaws_exact_desc desc;
	double w1 = nearbyint(WD_R / sqrt(2.0)), cx = 1000, cy = -2000, ang = 0.6, bx, by, half = 200;
	int k, i, j, n_wrong = 0, n_failed = 0;
	char buf[300];
	svg s;
	viewport vf = { 20, 90, 300, 300, -1.25 * WD_R, 1.25 * WD_R, -1.25 * WD_R, 1.25 * WD_R };
	qaws_exact_desc_default(&desc);
	desc.space_exp2 = 0;
	for (k = 0; k < 4; k++)
	{
		double a = k * 3.14159265358979 / 2, b = (k + 1) * 3.14159265358979 / 2;
		qaws_scalar cps[6], ws[3];
		qaws_rational_bezier_desc d;
		cps[0] = (qaws_scalar)(cx + nearbyint(WD_R * cos(a)));
		cps[1] = (qaws_scalar)(cy + nearbyint(WD_R * sin(a)));
		cps[2] = (qaws_scalar)(cx + nearbyint(WD_R * (cos(a) - sin(a))));
		cps[3] = (qaws_scalar)(cy + nearbyint(WD_R * (sin(a) + cos(a))));
		cps[4] = (qaws_scalar)(cx + nearbyint(WD_R * cos(b)));
		cps[5] = (qaws_scalar)(cy + nearbyint(WD_R * sin(b)));
		ws[0] = (qaws_scalar)WD_R;
		ws[1] = (qaws_scalar)w1;
		ws[2] = (qaws_scalar)WD_R;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D;
		d.degree = 2;
		d.control_points = cps;
		d.control_point_count = 3;
		d.weights = ws;
		d.weight_count = 3;
		qaws_curve_create_rational_bezier(&d, &pieces[k]);
		qaws_exact_curve_prepare(&desc, pieces[k], &exact[k], NULL);
	}
	/* zoom window centred on the boundary at angle `ang` (the conic is
	   within a fraction of a lattice unit of the circle there) */
	bx = cx + WD_R * cos(ang);
	by = cy + WD_R * sin(ang);
	for (j = 0; j < WD_N; j++)
		for (i = 0; i < WD_N; i++)
		{
			double p[2];
			int w = 0, kk = j * WD_N + i;
			qaws_status st;
			p[0] = bx - half + 2 * half * (i + 0.5) / WD_N;
			p[1] = by - half + 2 * half * (j + 0.5) / WD_N;
			sampled[kk] = (signed char)(wd_sampled(pieces, 4, 64, p) != 0);
			st = qaws_exact_winding_2d((qaws_exact_curve const* const*)exact, 4, p, &w);
			certified[kk] = (signed char)(st == QAWS_STATUS_OK ? (w != 0) : -1);
			n_failed += st != QAWS_STATUS_OK;
			wrong[kk] = (signed char)(certified[kk] >= 0 && sampled[kk] != certified[kk]);
			n_wrong += wrong[kk];
		}
	svg_open(&s, "showcase/exact3_winding.svg", 1280, 470, "Certified point-in-region on a loop of rational conics",
		"Four rational quadratics (radius 2^20 lattice units). Zoom: 400 x 400 lattice units across the boundary, 160 x 160 points. Blue: inside, light: outside.");
	svg_panel(&s, &vf, "the loop and the zoom window");
	{
		double xy[2 * 401];
		int q = 0;
		for (k = 0; k < 4; k++)
			for (i = 0; i <= 100; i++)
			{
				double e[2];
				if (k > 0 && i == 0)
					continue;
				qaws_exact_curve_evaluate(exact[k], i / 100.0, 0, e, NULL);
				xy[2 * q] = vx(&vf, e[0] - cx);
				xy[2 * q + 1] = vy(&vf, e[1] - cy);
				q++;
			}
		svg_polyline(&s, xy, q, "#24292f", 2, 1, 0);
		fprintf(s.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"8\" height=\"8\" fill=\"none\" stroke=\"#cf222e\" stroke-width=\"2\"/>\n",
			vx(&vf, bx - cx) - 4, vy(&vf, by - cy) - 4);
	}
	or_map(&s, 340, 120, 1.8, WD_N, sampled, in_colors, "sampled winding (256-point polyline)");
	or_map(&s, 650, 120, 1.8, WD_N, certified, in_colors, "qaws_exact_winding_2d (certified)");
	or_map(&s, 960, 120, 1.8, WD_N, wrong, wrong_colors, "points the sampled test gets wrong");
	sprintf(buf, "%d of %d points misclassified by the sampled test (the chord cuts up to ~80 lattice units inside the conic); %d undecided by the certified test (on the curve)",
		n_wrong, WD_N * WD_N, n_failed);
	svg_text(&s, 20, 445, 13, "#57606a", "start", buf);
	svg_close(&s);
	printf("exact3_winding: sampled wrong on %d of %d points, certified undecided %d\n", n_wrong, WD_N * WD_N, n_failed);
	for (k = 0; k < 4; k++)
	{
		qaws_exact_curve_destroy(exact[k]);
		qaws_curve_destroy(pieces[k]);
	}
}
