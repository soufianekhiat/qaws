/* ================================================================== */
/*  2. Photo vectorization with interpolating splines                 */
/* ================================================================== */

#define VZ_LEVELS 9
#define VZ_MAX_CURVES 400
#define VZ_MAX_PTS 160
#define VZ_ITERS 80

typedef struct contour
{
	int n;                 /* polyline points */
	double* xy;
	int closed;
	double level;
	int np;                /* interpolation points */
	qaws_scalar pts[VZ_MAX_PTS * 2];
	qaws_scalar init[VZ_MAX_PTS * 2];
	int family;            /* 0: centripetal Catmull-Rom, 1: Yuksel C2 */
	float color[3];
} contour;

/* Edge crossing point of an isoline on the pixel lattice. Horizontal edge
   (x,y)-(x+1,y) has id y*w+x, vertical (x,y)-(x,y+1) has id w*h+y*w+x. */
static void edge_point(float const* f, int w, int h, int id, double level, double* px, double* py)
{
	int vert = id >= w * h, k = vert ? id - w * h : id, x = k % w, y = k / w;
	double a = f[y * w + x], b = vert ? f[(y + 1) * w + x] : f[y * w + x + 1];
	double t = (level - a) / (b - a);
	(void)h;
	*px = x + (vert ? 0 : t);
	*py = y + (vert ? t : 0);
}

/* Marching squares at one level, chained into polylines. */
static int extract_contours(float const* f, int w, int h, double level, int min_points, contour* out, int cap)
{
	int nedges = 2 * w * h, x, y, i, count = 0, nseg = 0;
	int* link = (int*)malloc(sizeof(int) * (size_t)nedges * 2);   /* two neighbours per edge */
	unsigned char* used = (unsigned char*)calloc((size_t)nedges, 1);
	int* chain = (int*)malloc(sizeof(int) * (size_t)nedges);
	for (i = 0; i < nedges * 2; i++)
		link[i] = -1;
	for (y = 0; y + 1 < h; y++)
		for (x = 0; x + 1 < w; x++)
		{
			int e[4], c = 0, seg[4], ns = 0, s;
			double v00 = f[y * w + x], v10 = f[y * w + x + 1], v11 = f[(y + 1) * w + x + 1], v01 = f[(y + 1) * w + x];
			e[0] = y * w + x;                 /* bottom */
			e[1] = w * h + y * w + x + 1;     /* right */
			e[2] = (y + 1) * w + x;           /* top */
			e[3] = w * h + y * w + x;         /* left */
			if (v00 > level) c |= 1;
			if (v10 > level) c |= 2;
			if (v11 > level) c |= 4;
			if (v01 > level) c |= 8;
			switch (c)
			{
			case 1: case 14: seg[ns++] = e[3]; seg[ns++] = e[0]; break;
			case 2: case 13: seg[ns++] = e[0]; seg[ns++] = e[1]; break;
			case 3: case 12: seg[ns++] = e[3]; seg[ns++] = e[1]; break;
			case 4: case 11: seg[ns++] = e[1]; seg[ns++] = e[2]; break;
			case 6: case 9: seg[ns++] = e[0]; seg[ns++] = e[2]; break;
			case 7: case 8: seg[ns++] = e[3]; seg[ns++] = e[2]; break;
			case 5: case 10:
			{
				int center_high = (v00 + v10 + v11 + v01) * 0.25 > level;
				if ((c == 5) == center_high)
				{
					seg[ns++] = e[3]; seg[ns++] = e[2]; seg[ns++] = e[0]; seg[ns++] = e[1];
				}
				else
				{
					seg[ns++] = e[3]; seg[ns++] = e[0]; seg[ns++] = e[1]; seg[ns++] = e[2];
				}
				break;
			}
			default: break;
			}
			for (s = 0; s < ns; s += 2)
			{
				int a = seg[s], b = seg[s + 1];
				if (link[2 * a] < 0) link[2 * a] = b; else link[2 * a + 1] = b;
				if (link[2 * b] < 0) link[2 * b] = a; else link[2 * b + 1] = a;
				nseg++;
			}
		}
	for (i = 0; i < nedges && count < cap; i++)
	{
		int start = i, n = 0, cur, prev, closed = 0, k;
		if (used[i] || link[2 * i] < 0)
			continue;
		/* walk to one end of an open chain first */
		cur = i;
		prev = -1;
		for (;;)
		{
			int next = link[2 * cur] != prev ? link[2 * cur] : link[2 * cur + 1];
			if (next < 0 || next == i)
				break;
			prev = cur;
			cur = next;
			if (++n > nedges)
				break;
		}
		start = link[2 * cur] >= 0 && (link[2 * cur] == i || link[2 * cur + 1] == i) && n > 0 && link[2 * cur + 1] >= 0 ? i : cur;
		/* collect */
		n = 0;
		cur = start;
		prev = -1;
		while (cur >= 0 && !used[cur])
		{
			int next;
			used[cur] = 1;
			chain[n++] = cur;
			next = link[2 * cur] != prev ? link[2 * cur] : link[2 * cur + 1];
			if (next == prev)
				next = -1;
			prev = cur;
			if (next == start)
			{
				closed = 1;
				break;
			}
			cur = next;
		}
		if (n < min_points)
			continue;
		{
			contour* c = &out[count++];
			c->n = n;
			c->closed = closed;
			c->level = level;
			c->xy = (double*)malloc(sizeof(double) * 2 * (size_t)n);
			for (k = 0; k < n; k++)
				edge_point(f, w, h, chain[k], level, &c->xy[2 * k], &c->xy[2 * k + 1]);
		}
	}
	(void)nseg;
	free(link);
	free(used);
	free(chain);
	return count;
}

/* Points every `spacing` pixels of arc length along the polyline. */
static void contour_resample(contour* c, double spacing)
{
	double total = 0, acc = 0, step;
	int k, np, j = 0;
	for (k = 1; k < c->n; k++)
		total += hypot(c->xy[2 * k] - c->xy[2 * k - 2], c->xy[2 * k + 1] - c->xy[2 * k - 1]);
	if (c->closed)
		total += hypot(c->xy[0] - c->xy[2 * c->n - 2], c->xy[1] - c->xy[2 * c->n - 1]);
	np = (int)(total / spacing) + (c->closed ? 0 : 1);
	if (np < (c->closed ? 4 : 4)) np = 4;
	if (np > VZ_MAX_PTS - 2) np = VZ_MAX_PTS - 2;
	step = total / (c->closed ? np : np - 1);
	for (k = 0; k < np; k++)
	{
		double target = k * step;
		while (j + 1 < c->n)
		{
			double seg = hypot(c->xy[2 * j + 2] - c->xy[2 * j], c->xy[2 * j + 3] - c->xy[2 * j + 1]);
			if (acc + seg >= target)
				break;
			acc += seg;
			j++;
		}
		if (j + 1 < c->n)
		{
			double seg = hypot(c->xy[2 * j + 2] - c->xy[2 * j], c->xy[2 * j + 3] - c->xy[2 * j + 1]);
			double t = seg > 0 ? (target - acc) / seg : 0;
			c->pts[2 * k] = (qaws_scalar)(c->xy[2 * j] + t * (c->xy[2 * j + 2] - c->xy[2 * j]));
			c->pts[2 * k + 1] = (qaws_scalar)(c->xy[2 * j + 1] + t * (c->xy[2 * j + 3] - c->xy[2 * j + 1]));
		}
		else
		{
			c->pts[2 * k] = (qaws_scalar)c->xy[2 * c->n - 2];
			c->pts[2 * k + 1] = (qaws_scalar)c->xy[2 * c->n - 1];
		}
	}
	c->np = np;
	if (!c->closed)
	{
		/* phantom end points so the open spline reaches both ends */
		memmove(c->pts + 2, c->pts, sizeof(qaws_scalar) * 2 * (size_t)np);
		c->pts[0] = 2 * c->pts[2] - c->pts[4];
		c->pts[1] = 2 * c->pts[3] - c->pts[5];
		c->pts[2 * (np + 1)] = 2 * c->pts[2 * np] - c->pts[2 * np - 2];
		c->pts[2 * (np + 1) + 1] = 2 * c->pts[2 * np + 1] - c->pts[2 * np - 1];
		c->np = np + 2;
	}
	memcpy(c->init, c->pts, sizeof(qaws_scalar) * 2 * (size_t)c->np);
}

static qaws_curve* contour_curve(contour const* c, qaws_scalar const* pts)
{
	qaws_catmull_rom_desc d;
	qaws_curve* cr = NULL;
	if (c->family == 1)
	{
		qaws_yuksel_desc y;
		memset(&y, 0, sizeof(y));
		y.dimension = QAWS_DIMENSION_2D;
		y.control_points = pts;
		y.control_point_count = (unsigned int)c->np;
		y.mode = QAWS_YUKSEL_MODE_BEZIER;
		y.closed = c->closed;
		qaws_curve_create_yuksel(&y, &cr);
		return cr;
	}
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D;
	d.control_points = pts;
	d.control_point_count = (unsigned int)c->np;
	d.parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL;
	d.closed = c->closed;
	qaws_curve_create_catmull_rom(&d, &cr);
	return cr;
}

/* Mean (I(C) - level)^2 and mean distance estimate |I - level| / |grad I|. */
static double contour_energy(float const* f, int w, int h, contour const* c, qaws_scalar const* pts,
	qaws_scalar* grad, double* out_dist)
{
	enum { MAXS = VZ_MAX_PTS * 4 };
	static qaws_scalar ts[MAXS];
	static qaws_curve_jet_2d prim[MAXS], tan[MAXS], bar[MAXS];
	qaws_curve* cr = contour_curve(c, pts);
	qaws_range r = qaws_curve_get_parameter_range(cr);
	unsigned int m = (unsigned int)(r.max_value * 4), k;
	double e = 0, dist = 0, spring = 2e-6;
	if (m > MAXS) m = MAXS;
	for (k = 0; k < m; k++)
		ts[k] = (qaws_scalar)(r.max_value * (k + 0.5) / m);
	qaws_curve_eval_batch_tangent_2d(NULL, cr, ts, NULL, m, QAWS_EVAL_FLAG_POSITION, NULL, prim, tan);
	for (k = 0; k < m; k++)
	{
		double gx, gy, v = sample(f, w, h, prim[k].d[0].x, prim[k].d[0].y, &gx, &gy) - c->level;
		e += v * v / m;
		dist += fabs(v) / (sqrt(gx * gx + gy * gy) + 1e-6) / m;
		memset(&bar[k], 0, sizeof(bar[k]));
		bar[k].d[0].x = (qaws_scalar)(2 * v * gx / m);
		bar[k].d[0].y = (qaws_scalar)(2 * v * gy / m);
		bar[k].channels = QAWS_EVAL_FLAG_POSITION;
	}
	if (grad)
	{
		qaws_field_view fv;
		qaws_diff_views views = one_field(&fv, QAWS_FIELD_POINTS, grad, (unsigned int)c->np, 2);
		int i;
		memset(grad, 0, sizeof(qaws_scalar) * 2 * (size_t)c->np);
		qaws_curve_eval_batch_adjoint_2d(NULL, cr, ts, m, QAWS_EVAL_FLAG_POSITION, bar, &views, NULL);
		/* keep the points from sliding: spring to their initial positions */
		for (i = 0; i < 2 * c->np; i++)
		{
			double d = pts[i] - c->init[i];
			e += spring * d * d;
			grad[i] += (qaws_scalar)(2 * spring * d);
		}
	}
	if (out_dist)
		*out_dist = dist;
	qaws_curve_destroy(cr);
	return e;
}

static void contour_svg(svg* s, imgmap const* m, qaws_curve const* c, char const* color, double width)
{
	qaws_range r = qaws_curve_get_parameter_range(c);
	int n = (int)(r.max_value * 8) + 2, k;
	double* xy = (double*)malloc(sizeof(double) * 2 * (size_t)n);
	for (k = 0; k < n; k++)
	{
		qaws_eval_result_2d e;
		qaws_curve_evaluate_2d(c, (qaws_scalar)(r.max_value * k / (n - 1)), QAWS_EVAL_FLAG_POSITION, &e);
		map_pt(m, e.position.x, e.position.y, &xy[2 * k], &xy[2 * k + 1]);
	}
	svg_polyline(s, xy, n, color, width, 1.0, 0);
	free(xy);
}

/* Signed curvature of a 2D curve at t. */
static double curve_kappa(qaws_curve const* c, double t)
{
	qaws_eval_result_2d e;
	double s;
	qaws_curve_evaluate_2d(c, (qaws_scalar)t, QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2, &e);
	s = sqrt(e.d1.x * e.d1.x + e.d1.y * e.d1.y);
	return s > 1e-12 ? (e.d1.x * e.d2.y - e.d1.y * e.d2.x) / (s * s * s) : 0;
}

/* Mean |curvature jump| across the interpolation points (1/px). */
static double curvature_jumps(qaws_curve const* c, int* count)
{
	qaws_range r = qaws_curve_get_parameter_range(c);
	int k, n = (int)(r.max_value + 0.5);
	double acc = 0;
	*count = 0;
	for (k = 1; k < n; k++)
	{
		acc += fabs(curve_kappa(c, k - 1e-4) - curve_kappa(c, k + 1e-4));
		(*count)++;
	}
	return acc;
}

/* Curve plus its curvature comb (normals scaled by curvature). */
static void comb_svg(svg* s, imgmap const* m, qaws_curve const* c, double comb_scale, char const* color, char const* comb_color)
{
	qaws_range r = qaws_curve_get_parameter_range(c);
	int n = (int)(r.max_value * 12) + 1, k;
	double* tips = (double*)malloc(sizeof(double) * 2 * (size_t)n);
	for (k = 0; k < n; k++)
	{
		qaws_eval_result_2d e;
		double t = r.max_value * k / (n - 1), sp, nx, ny, kap, x0, y0;
		qaws_curve_evaluate_2d(c, (qaws_scalar)t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2, &e);
		sp = sqrt(e.d1.x * e.d1.x + e.d1.y * e.d1.y) + 1e-12;
		nx = -e.d1.y / sp;
		ny = e.d1.x / sp;
		kap = (e.d1.x * e.d2.y - e.d1.y * e.d2.x) / (sp * sp * sp);
		/* clamp the comb length so tight turns do not swamp the figure */
		if (comb_scale * kap > 9) kap = 9 / comb_scale;
		if (comb_scale * kap < -9) kap = -9 / comb_scale;
		map_pt(m, e.position.x, e.position.y, &x0, &y0);
		map_pt(m, e.position.x - comb_scale * kap * nx, e.position.y - comb_scale * kap * ny, &tips[2 * k], &tips[2 * k + 1]);
		svg_line(s, x0, y0, tips[2 * k], tips[2 * k + 1], comb_color, 0.7, 0.55);
	}
	svg_polyline(s, tips, n, comb_color, 1.0, 0.9, 0);
	contour_svg(s, m, c, color, 2.0);
	free(tips);
}

/* Same contours, same points, same energy: centripetal Catmull-Rom (C1)
   against Yuksel C2 interpolating splines. */
static void yuksel_compare(contour* cs, int nc, float const* f, image const* img)
{
	static contour ys[VZ_MAX_CURVES];
	double dist_cr = 0, dist_yk = 0, jump_cr = 0, jump_yk = 0;
	int ncr = 0, nyk = 0, i, it;
	svg s;
	char buf[256];

	for (i = 0; i < nc; i++)
	{
		contour* c = &ys[i];
		qaws_scalar g[VZ_MAX_PTS * 2];
		adam opt;
		double dd;
		int cnt;
		*c = cs[i];
		c->family = 1;
		if (!c->closed)
		{
			/* interpolating splines through the real points: drop the phantoms */
			c->np -= 2;
			memmove(c->pts, cs[i].init + 2, sizeof(qaws_scalar) * 2 * (size_t)c->np);
			memcpy(c->init, c->pts, sizeof(qaws_scalar) * 2 * (size_t)c->np);
		}
		else
		{
			memcpy(c->pts, cs[i].init, sizeof(qaws_scalar) * 2 * (size_t)c->np);
		}
		memset(&opt, 0, sizeof(opt));
		for (it = 0; it < VZ_ITERS; it++)
		{
			contour_energy(f, img->w, img->h, c, c->pts, g, NULL);
			adam_step(&opt, c->pts, g, 2 * c->np, 0.25 * (1.0 - 0.8 * it / (double)VZ_ITERS));
		}
		contour_energy(f, img->w, img->h, c, c->pts, NULL, &dd);
		dist_yk += dd;
		contour_energy(f, img->w, img->h, &cs[i], cs[i].pts, NULL, &dd);
		dist_cr += dd;
		{
			qaws_curve* a = contour_curve(&cs[i], cs[i].pts);
			qaws_curve* b = contour_curve(c, c->pts);
			jump_cr += curvature_jumps(a, &cnt);
			ncr += cnt;
			jump_yk += curvature_jumps(b, &cnt);
			nyk += cnt;
			qaws_curve_destroy(a);
			qaws_curve_destroy(b);
		}
	}

	svg_open(&s, "showcase/app2b_yuksel_vs_catmull_rom.svg", 1240, 700, "Interpolating splines: Catmull-Rom vs Yuksel C2",
		"Same isocontours, same interpolation points, same optimization; curvature combs show C1 kinks of centripetal "
		"Catmull-Rom against the C2 Yuksel splines (curvature peaks at the points).");
	{
		double zx0 = 170, zy0 = 105, zs = 5.2, zw = 112, zh = 84;
		imgmap ma = { 20, 100, zx0, zy0, zs }, mb = { 640, 100, zx0, zy0, zs };
		double zpw = zw * zs, zph = zh * zs;
		int p;
		for (p = 0; p < 2; p++)
		{
			imgmap const* m = p ? &mb : &ma;
			fprintf(s.f, "<text x=\"%.1f\" y=\"92\" font-size=\"13\" font-weight=\"600\" fill=\"%s\">%s</text>\n", m->x0,
				p ? "#0969da" : "#9a6700", p ? "Yuksel C2 (Bezier mode), optimized through its new adjoints"
				: "centripetal Catmull-Rom (C1), optimized");
			fprintf(s.f, "<clipPath id=\"ykclip%d\"><rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\"/></clipPath>\n",
				p, m->x0, m->y0, zpw, zph);
			fprintf(s.f, "<g clip-path=\"url(#ykclip%d)\">\n", p);
			fprintf(s.f, "<image href=\"../%s/folds.png\" x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" opacity=\"0.35\"/>\n",
				g_photos, m->x0 - zx0 * zs, m->y0 - zy0 * zs, img->w * zs, img->h * zs);
			for (i = 0; i < nc; i++)
			{
				contour const* c = p ? &ys[i] : &cs[i];
				qaws_curve* cv;
				int k, inside = 0, rank = 0, j;
				/* only the 5 contours with the most points in the crop */
				for (k = 0; k < cs[i].np; k++)
					inside += cs[i].init[2 * k] > zx0 && cs[i].init[2 * k] < zx0 + zw &&
						cs[i].init[2 * k + 1] > zy0 && cs[i].init[2 * k + 1] < zy0 + zh;
				for (j = 0; j < nc; j++)
				{
					int in2 = 0;
					for (k = 0; k < cs[j].np; k++)
						in2 += cs[j].init[2 * k] > zx0 && cs[j].init[2 * k] < zx0 + zw &&
							cs[j].init[2 * k + 1] > zy0 && cs[j].init[2 * k + 1] < zy0 + zh;
					rank += in2 > inside || (in2 == inside && j < i);
				}
				if (rank >= 5 || inside == 0)
					continue;
				cv = contour_curve(c, c->pts);
				comb_svg(&s, m, cv, 260.0, p ? "#0969da" : "#9a6700", p ? "#54aeff" : "#d4a72c");
				for (k = 0; k < c->np; k++)
				{
					double sx, sy;
					map_pt(m, c->pts[2 * k], c->pts[2 * k + 1], &sx, &sy);
					svg_circle(&s, sx, sy, 2.4, "#ffffff", "#24292f");
				}
				qaws_curve_destroy(cv);
			}
			fprintf(s.f, "</g>\n");
			fprintf(s.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"none\" stroke=\"#d0d7de\"/>\n", m->x0, m->y0, zpw, zph);
		}
		sprintf(buf, "mean curvature jump at the points: Catmull-Rom %.4f / px, Yuksel %.2e / px", jump_cr / ncr, jump_yk / nyk);
		svg_text(&s, 20, 100 + zph + 34, 15, "#24292f", "start", buf);
		sprintf(buf, "mean distance to the isocontour after optimization: Catmull-Rom %.3f px, Yuksel %.3f px", dist_cr / nc, dist_yk / nc);
		svg_text(&s, 20, 100 + zph + 60, 15, "#0969da", "start", buf);
		svg_text(&s, 20, 100 + zph + 90, 13, "#57606a", "start",
			"Yuksel adjoints: implicit derivative of the max-curvature parameter (cubic root), dual-number P1, sub-curve "
			"reparameterization and trigonometric blend");
	}
	svg_close(&s);
	printf("2b_yuksel: curvature jump CR %.4f vs Yuksel %.2e /px, distance CR %.3f vs Yuksel %.3f px\n",
		jump_cr / ncr, jump_yk / nyk, dist_cr / nc, dist_yk / nc);
}

static void app_vectorize(void)
{
	static contour cs[VZ_MAX_CURVES];
	char path[512], buf[256];
	image img;
	float* f;
	int nc = 0, l, i, total_poly = 0, total_pts = 0;
	double d0 = 0, d1 = 0;
	svg s;

	sprintf(path, "%s/folds.ppm", g_photos);
	if (!image_load_ppm(path, &img))
	{
		printf("2_vectorize: %s not found (see examples/photo_to_ppm.ps1)\n", path);
		return;
	}
	f = (float*)malloc(sizeof(float) * (size_t)img.w * img.h);
	memcpy(f, img.lum, sizeof(float) * (size_t)img.w * img.h);
	blur(f, img.w, img.h, 1.6);
	{
		/* levels at luminance quantiles */
		int n = img.w * img.h, hist[256] = { 0 }, acc = 0, q = 1, b;
		double levels[VZ_LEVELS];
		for (i = 0; i < n; i++)
			hist[(int)(f[i] * 255.0f) > 255 ? 255 : (int)(f[i] * 255.0f)]++;
		for (b = 0; b < 256 && q <= VZ_LEVELS; b++)
		{
			acc += hist[b];
			while (q <= VZ_LEVELS && acc >= n * q / (VZ_LEVELS + 1))
				levels[q++ - 1] = (b + 0.5) / 255.0;
		}
		for (l = 0; l < VZ_LEVELS; l++)
			nc += extract_contours(f, img.w, img.h, levels[l], 50, cs + nc, VZ_MAX_CURVES - nc);
	}
	for (i = 0; i < nc; i++)
	{
		contour* c = &cs[i];
		qaws_scalar g[VZ_MAX_PTS * 2];
		adam opt;
		double dd;
		int it, k;
		contour_resample(c, 20.0);
		contour_energy(f, img.w, img.h, c, c->pts, NULL, &dd);
		d0 += dd;
		memset(&opt, 0, sizeof(opt));
		for (it = 0; it < VZ_ITERS; it++)
		{
			contour_energy(f, img.w, img.h, c, c->pts, g, NULL);
			adam_step(&opt, c->pts, g, 2 * c->np, 0.25 * (1.0 - 0.8 * it / (double)VZ_ITERS));
		}
		contour_energy(f, img.w, img.h, c, c->pts, NULL, &dd);
		d1 += dd;
		total_poly += c->n;
		total_pts += c->np;
		{
			/* mean photo color along the contour */
			double rgb[3] = { 0, 0, 0 };
			for (k = 0; k < c->n; k++)
			{
				int x = (int)c->xy[2 * k], y = (int)c->xy[2 * k + 1], ch;
				for (ch = 0; ch < 3; ch++)
					rgb[ch] += img.rgb[(y * img.w + x) * 3 + ch] / c->n;
			}
			for (k = 0; k < 3; k++)
				c->color[k] = (float)rgb[k];
		}
	}

	svg_open(&s, "showcase/app2_vectorize.svg", 1240, 840, "Photo vectorization with interpolating splines",
		"Isocontours of the photo become centripetal Catmull-Rom splines; their interpolation points are optimized "
		"onto the contours through Catmull-Rom adjoints and the image gradient.");
	{
		double pw = 590, sc = pw / img.w, ph = img.h * sc;
		imgmap ma = { 20, 100, 0, 0, sc }, mb = { 630, 100, 0, 0, sc };
		double zx0 = 150, zy0 = 95, zs = 3.0, zw = 160, zh = 95;
		imgmap mz = { 20, 100 + ph + 50, zx0, zy0, zs };
		double zpw = zw * zs, zph = zh * zs;
		(void)ma;
		fprintf(s.f, "<text x=\"20\" y=\"92\" font-size=\"13\" font-weight=\"600\" fill=\"#24292f\">photo</text>\n");
		fprintf(s.f, "<text x=\"630\" y=\"92\" font-size=\"13\" font-weight=\"600\" fill=\"#24292f\">%d splines through %d points (vector)</text>\n", nc, total_pts);
		fprintf(s.f, "<image href=\"../%s/folds.png\" x=\"20\" y=\"100\" width=\"%.1f\" height=\"%.1f\"/>\n", g_photos, pw, ph);
		fprintf(s.f, "<rect x=\"630\" y=\"100\" width=\"%.1f\" height=\"%.1f\" fill=\"#ffffff\" stroke=\"#d0d7de\"/>\n", pw, ph);
		for (i = 0; i < nc; i++)
		{
			char col[40];
			qaws_curve* cr = contour_curve(&cs[i], cs[i].pts);
			double shade = 0.75;
			sprintf(col, "rgb(%d,%d,%d)", (int)(255 * cs[i].color[0] * shade), (int)(255 * cs[i].color[1] * shade),
				(int)(255 * cs[i].color[2] * shade));
			contour_svg(&s, &mb, cr, col, 1.4);
			qaws_curve_destroy(cr);
		}
		fprintf(s.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"none\" stroke=\"#cf222e\" stroke-width=\"1.5\"/>\n",
			20 + zx0 * sc, 100 + zy0 * sc, zw * sc, zh * sc);

		fprintf(s.f, "<text x=\"20\" y=\"%.1f\" font-size=\"13\" font-weight=\"600\" fill=\"#24292f\">zoom: isocontour (grey), "
			"initial spline (orange), optimized spline and points (blue)</text>\n", mz.y0 - 8);
		fprintf(s.f, "<clipPath id=\"vzclip\"><rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\"/></clipPath>\n", mz.x0, mz.y0, zpw, zph);
		fprintf(s.f, "<g clip-path=\"url(#vzclip)\">\n");
		fprintf(s.f, "<image href=\"../%s/folds.png\" x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" opacity=\"0.45\"/>\n",
			g_photos, mz.x0 - zx0 * zs, mz.y0 - zy0 * zs, img.w * zs, img.h * zs);
		for (i = 0; i < nc; i++)
		{
			contour* c = &cs[i];
			double* xy = (double*)malloc(sizeof(double) * 2 * (size_t)c->n);
			qaws_curve *c0 = contour_curve(c, c->init), *c1 = contour_curve(c, c->pts);
			int k;
			for (k = 0; k < c->n; k++)
				map_pt(&mz, c->xy[2 * k], c->xy[2 * k + 1], &xy[2 * k], &xy[2 * k + 1]);
			svg_polyline(&s, xy, c->n, "#57606a", 1.0, 0.9, 1);
			contour_svg(&s, &mz, c0, "#ff8c42", 2.0);
			contour_svg(&s, &mz, c1, "#0969da", 2.2);
			for (k = 0; k < c->np; k++)
			{
				double sx, sy;
				map_pt(&mz, c->pts[2 * k], c->pts[2 * k + 1], &sx, &sy);
				svg_circle(&s, sx, sy, 2.6, "#ffffff", "#0969da");
			}
			qaws_curve_destroy(c0);
			qaws_curve_destroy(c1);
			free(xy);
		}
		fprintf(s.f, "</g>\n");
		fprintf(s.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"none\" stroke=\"#d0d7de\"/>\n", mz.x0, mz.y0, zpw, zph);
		{
			double tx = mz.x0 + zpw + 30, ty = mz.y0 + 20;
			sprintf(buf, "%d contour polylines (%d vertices) -> %d spline points", nc, total_poly, total_pts);
			svg_text(&s, tx, ty, 15, "#24292f", "start", buf);
			sprintf(buf, "mean distance to the isocontour: %.3f px -> %.3f px", d0 / nc, d1 / nc);
			svg_text(&s, tx, ty + 26, 15, "#0969da", "start", buf);
			svg_text(&s, tx, ty + 62, 13, "#57606a", "start", "pipeline");
			svg_text(&s, tx, ty + 82, 13, "#24292f", "start", "1. blur, 9 luminance levels at quantiles");
			svg_text(&s, tx, ty + 102, 13, "#24292f", "start", "2. marching squares, chained into closed / open polylines");
			svg_text(&s, tx, ty + 122, 13, "#24292f", "start", "3. one interpolation point every 20 px (centripetal Catmull-Rom)");
			svg_text(&s, tx, ty + 142, 13, "#24292f", "start", "4. minimize mean (I(C(t)) - level)^2 over the points:");
			svg_text(&s, tx, ty + 162, 13, "#24292f", "start", "   2 (I - level) grad I -> jet adjoint -> qaws_curve_eval_batch_adjoint_2d");
			svg_text(&s, tx, ty + 182, 13, "#24292f", "start", "   -> QAWS_FIELD_POINTS (non-linear centripetal rule)");
		}
	}
	svg_close(&s);
	yuksel_compare(cs, nc, f, &img);
	printf("2_vectorize: %d curves, %d polyline vertices -> %d points, distance %.3f -> %.3f px\n",
		nc, total_poly, total_pts, d0 / nc, d1 / nc);
	for (i = 0; i < nc; i++)
		free(cs[i].xy);
	free(f);
	image_free(&img);
}

