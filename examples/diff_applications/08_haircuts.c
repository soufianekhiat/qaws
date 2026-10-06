/* ================================================================== */
/*  8. Real haircuts: vector hair strands from portrait photos        */
/* ================================================================== */

typedef struct hair_spec
{
	char const* name;
	char const* label;
	int hair[8][2], nhair;       /* hair color seeds (pixels of the 360-wide photo) */
	int non[8][2], nnon;         /* skin / clothes / background seeds */
	int box[4];                  /* x0, y0, x1, y1: where hair may be */
	int face[4];                 /* face ellipse: cx, cy, rx, ry (never hair) */
} hair_spec;

static hair_spec const g_hair_specs[6] = {
	{ "hair_straight", "straight, shoulder length",
		{ { 120, 90 }, { 100, 180 }, { 240, 200 }, { 170, 60 }, { 110, 240 }, { 250, 240 } }, 6,
		{ { 175, 150 }, { 175, 200 }, { 40, 80 }, { 330, 100 }, { 180, 330 }, { 60, 300 }, { 300, 300 } }, 7,
		{ 60, 20, 300, 252 }, { 175, 150, 42, 55 } },
	{ "hair_short_curly", "short, curly",
		{ { 120, 40 }, { 80, 90 }, { 150, 30 }, { 80, 140 }, { 190, 40 } }, 5,
		{ { 130, 110 }, { 150, 150 }, { 60, 250 }, { 200, 200 }, { 330, 60 }, { 20, 100 }, { 260, 120 } }, 7,
		{ 40, 0, 240, 175 }, { 135, 118, 38, 52 } },
	{ "hair_coily", "voluminous, coily",
		{ { 100, 120 }, { 300, 120 }, { 200, 70 }, { 80, 190 }, { 300, 190 }, { 260, 60 }, { 130, 60 } }, 7,
		{ { 190, 150 }, { 190, 200 }, { 160, 170 }, { 220, 170 }, { 190, 115 }, { 40, 60 }, { 340, 250 }, { 200, 330 } }, 8,
		{ 50, 30, 345, 265 }, { 190, 165, 48, 62 } },
	{ "hair_long", "long, straight (profile)",
		{ { 260, 80 }, { 300, 250 }, { 320, 400 }, { 250, 160 }, { 200, 40 }, { 240, 300 } }, 6,
		{ { 200, 180 }, { 100, 40 }, { 60, 140 }, { 100, 330 }, { 220, 380 }, { 40, 420 } }, 6,
		{ 150, 0, 360, 456 }, { 195, 185, 45, 72 } },
	{ "hair_twists", "long twists",
		{ { 140, 140 }, { 110, 250 }, { 280, 180 }, { 270, 260 }, { 210, 80 }, { 150, 90 }, { 100, 330 } }, 7,
		{ { 210, 180 }, { 220, 220 }, { 180, 190 }, { 240, 190 }, { 210, 135 }, { 60, 100 }, { 330, 150 }, { 240, 400 } }, 8,
		{ 60, 60, 330, 440 }, { 213, 185, 45, 58 } },
	{ "hair_wavy", "long, wavy, grey (with beard)",
		{ { 100, 60 }, { 200, 30 }, { 60, 250 }, { 320, 240 }, { 300, 120 }, { 50, 150 }, { 180, 330 } }, 7,
		{ { 180, 180 }, { 150, 170 }, { 80, 400 }, { 340, 420 }, { 350, 30 } }, 5,
		{ 15, 0, 360, 420 }, { 185, 175, 52, 58 } }
};

/* Seeded hair mask: patch-averaged color, nearest hair seed vs nearest
   other seed, restricted to the box, kept connected to the hair seeds,
   then softened. */
static float* hair_mask(image const* img, hair_spec const* sp)
{
	int w = img->w, h = img->h, n = w * h, i, x, y, k;
	float* avg = (float*)malloc(sizeof(float) * 3 * (size_t)n);
	float* m = (float*)calloc((size_t)n, sizeof(float));
	unsigned char* lab = (unsigned char*)calloc((size_t)n, 1);
	int* stack = (int*)malloc(sizeof(int) * (size_t)n);
	double hc[8][3], nc[8][3];
	int top = 0;
	for (k = 0; k < 3; k++)
	{
		float* ch = (float*)malloc(sizeof(float) * (size_t)n);
		for (i = 0; i < n; i++)
			ch[i] = img->rgb[3 * i + k];
		blur(ch, w, h, 2.0);
		for (i = 0; i < n; i++)
			avg[3 * i + k] = ch[i];
		free(ch);
	}
	for (k = 0; k < sp->nhair; k++)
		for (i = 0; i < 3; i++)
			hc[k][i] = avg[3 * (sp->hair[k][1] * w + sp->hair[k][0]) + i];
	for (k = 0; k < sp->nnon; k++)
		for (i = 0; i < 3; i++)
			nc[k][i] = avg[3 * (sp->non[k][1] * w + sp->non[k][0]) + i];
	for (y = 0; y < h; y++)
		for (x = 0; x < w; x++)
		{
			double dh = 1e9, dn = 1e9;
			float const* c = &avg[3 * (y * w + x)];
			if (x < sp->box[0] || x >= sp->box[2] || y < sp->box[1] || y >= sp->box[3])
				continue;
			{
				double ex = (x - sp->face[0]) / (double)sp->face[2], ey = (y - sp->face[1]) / (double)sp->face[3];
				if (ex * ex + ey * ey < 1)
					continue;
			}
			for (k = 0; k < sp->nhair; k++)
			{
				double d = (c[0] - hc[k][0]) * (c[0] - hc[k][0]) + (c[1] - hc[k][1]) * (c[1] - hc[k][1]) + (c[2] - hc[k][2]) * (c[2] - hc[k][2]);
				if (d < dh) dh = d;
			}
			for (k = 0; k < sp->nnon; k++)
			{
				double d = (c[0] - nc[k][0]) * (c[0] - nc[k][0]) + (c[1] - nc[k][1]) * (c[1] - nc[k][1]) + (c[2] - nc[k][2]) * (c[2] - nc[k][2]);
				if (d < dn) dn = d;
			}
			lab[y * w + x] = dh < dn ? 1 : 0;
		}
	/* connected to the hair seeds */
	for (k = 0; k < sp->nhair; k++)
	{
		int id = sp->hair[k][1] * w + sp->hair[k][0];
		if (lab[id] == 1)
		{
			lab[id] = 2;
			stack[top++] = id;
		}
	}
	while (top > 0)
	{
		int id = stack[--top], px = id % w, py = id / w, d;
		static int const nb[4][2] = { { 1, 0 }, { -1, 0 }, { 0, 1 }, { 0, -1 } };
		m[id] = 1;
		for (d = 0; d < 4; d++)
		{
			int qx = px + nb[d][0], qy = py + nb[d][1], q;
			if (qx < 0 || qy < 0 || qx >= w || qy >= h)
				continue;
			q = qy * w + qx;
			if (lab[q] == 1)
			{
				lab[q] = 2;
				stack[top++] = q;
			}
		}
	}
	blur(m, w, h, 1.5);
	free(avg);
	free(lab);
	free(stack);
	return m;
}

/* Fits and optimizes strands (as in application 1); accumulates stats. */
static void hair_fit_strands(tensor_field const* tf, strand* strands, int ns, double* a0, double* a1, int* points, unsigned int* cps_total)
{
	int i;
	*a0 = *a1 = 0;
	*points = 0;
	*cps_total = 0;
	for (i = 0; i < ns; i++)
	{
		strand* st = &strands[i];
		unsigned int ncp = strand_cp_count(st), got = 0;
		qaws_bspline_fit_desc d;
		qaws_scalar cps[48], g[48];
		adam opt;
		double ea, aa;
		int it;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D;
		d.data_points = st->xy;
		d.data_point_count = (unsigned int)st->n;
		d.degree = 3;
		d.control_point_count = ncp;
		st->fit = st->opt = NULL;
		if (qaws_curve_fit_bspline(&d, &st->fit) != QAWS_STATUS_OK)
		{
			st->fit = NULL;
			continue;
		}
		qaws_curve_read_field(st->fit, QAWS_FIELD_CONTROL_POINTS, cps, ncp * 2, &got);
		qaws_curve_read_field(st->fit, QAWS_FIELD_KNOTS, st->knots, ncp + 4, &got);
		qaws_curve_destroy(st->fit);
		st->fit = strand_curve(cps, ncp, st->knots);
		strand_alignment(tf, st->fit, &ea, &aa);
		*a0 += aa;
		memset(&opt, 0, sizeof(opt));
		for (it = 0; it < HS_ITERS; it++)
		{
			memset(g, 0, sizeof(g));
			strand_energy(tf, st, cps, ncp, g);
			adam_step(&opt, cps, g, (int)(2 * ncp), 0.25 * (1.0 - 0.8 * it / (double)HS_ITERS));
		}
		st->opt = strand_curve(cps, ncp, st->knots);
		strand_alignment(tf, st->opt, &ea, &aa);
		*a1 += aa;
		*points += st->n;
		*cps_total += ncp;
	}
	if (ns)
	{
		*a0 /= ns;
		*a1 /= ns;
	}
}

/* Strand drawn in segments colored by direction (hue = angle mod 180). */
static void strand_svg_hue(svg* s, imgmap const* m, qaws_curve const* c, double width)
{
	qaws_range r = qaws_curve_get_parameter_range(c);
	double px = 0, py = 0;
	int k, n = 40;
	for (k = 0; k <= n; k++)
	{
		qaws_eval_result_2d e;
		double sx, sy;
		qaws_curve_evaluate_2d(c, (qaws_scalar)(r.min_value + (r.max_value - r.min_value) * k / n),
			QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &e);
		map_pt(m, e.position.x, e.position.y, &sx, &sy);
		if (k > 0)
		{
			char col[48];
			double ang = atan2(e.d1.y, e.d1.x) * 180.0 / PI;
			if (ang < 0) ang += 180;
			sprintf(col, "hsl(%d,85%%,60%%)", (int)(2 * ang) % 360);
			svg_line(s, px, py, sx, sy, col, width, 0.95);
		}
		px = sx;
		py = sy;
	}
}

static void app_haircuts(void)
{
	static strand strands[HS_MAX_STRANDS];
	svg s;
	char buf[256], path[512];
	int p;

	svg_open(&s, "showcase/app8_haircuts.svg", 1240, 795, "Real haircuts: vector hair from portrait photos",
		"Seeded hair mask, orientation field, evenly spaced streamlines in the mask, B-spline strands optimized to "
		"follow the hair through unit-tangent adjoints (right: strands colored by direction).");
	for (p = 0; p < 6; p++)
	{
		hair_spec const* sp = &g_hair_specs[p];
		image img, masked;
		tensor_field tf;
		float* mask;
		int ns, i, pts;
		unsigned int ncps;
		double a0, a1, sc, cw = 400, x0 = 20 + (p % 3) * 405, y0 = 90 + (p / 3) * 342;
		imgmap mp, mv;

		sprintf(path, "%s/%s.ppm", g_photos, sp->name);
		if (!image_load_ppm(path, &img))
		{
			printf("8_haircuts: %s missing\n", path);
			continue;
		}
		mask = hair_mask(&img, sp);
		tensor_build(&img, 1.0, 2.5, &tf);
		/* tracing reads the mask through the luminance channel */
		masked = img;
		masked.lum = mask;
		ns = trace_strands(&masked, &tf, 2.0, 0.2, 0.5, strands, HS_MAX_STRANDS);
		hair_fit_strands(&tf, strands, ns, &a0, &a1, &pts, &ncps);

		sc = (cw / 2 - 4) / img.w;
		mp.x0 = x0; mp.y0 = y0 + 22; mp.ox = 0; mp.oy = 0; mp.scale = sc;
		mv = mp;
		mv.x0 = x0 + cw / 2 + 2;
		sprintf(buf, "%s", sp->label);
		svg_text(&s, x0, y0 + 14, 13, "#24292f", "start", buf);
		fprintf(s.f, "<image href=\"../%s/%s.png\" x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\"/>\n",
			g_photos, sp->name, mp.x0, mp.y0, img.w * sc, img.h * sc);
		/* mask outline: dim the non-hair pixels on the right panel */
		fprintf(s.f, "<image href=\"../%s/%s.png\" x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" opacity=\"0.22\"/>\n",
			g_photos, sp->name, mv.x0, mv.y0, img.w * sc, img.h * sc);
		fprintf(s.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"#0d1117\" opacity=\"0.85\"/>\n",
			mv.x0, mv.y0, img.w * sc, img.h * sc);
		for (i = 0; i < ns; i++)
			if (strands[i].opt)
				strand_svg_hue(&s, &mv, strands[i].opt, 1.1);
		sprintf(buf, "%d strands, %d points -> %u control points; angle %.1f -> %.1f deg", ns, pts, ncps, a0, a1);
		svg_text(&s, x0, y0 + 22 + img.h * sc + 16, 11, "#57606a", "start", buf);
		printf("8_haircuts: %-18s %4d strands, %6d points -> %4u CPs, angle %.2f -> %.2f deg\n", sp->name, ns, pts, ncps, a0, a1);

		for (i = 0; i < ns; i++)
		{
			free(strands[i].xy);
			if (strands[i].fit) qaws_curve_destroy(strands[i].fit);
			if (strands[i].opt) qaws_curve_destroy(strands[i].opt);
		}
		free(mask);
		free(tf.xx); free(tf.xy); free(tf.yy);
		image_free(&img);
	}
	svg_text(&s, 20, 782, 11, "#57606a", "start",
		"Photos: Wikimedia Commons (CC0, public domain, CC BY 2.0, CC BY-SA 2.0/4.0); authors and licenses in photos/CREDITS.txt.");
	svg_close(&s);
}

