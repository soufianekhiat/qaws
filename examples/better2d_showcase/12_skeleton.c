/* Figure 12: straight skeletons (qaws_skeleton). Each panel shows the faces
   lit as the roof they make (each face rises at its edge's weight), the
   skeleton over them, and the wavefront at eight times as contours coloured
   by time. Plans, a plus, a comb, a star with a hole, random polygons, unions
   of grid rectangles, and weights. */

typedef struct sk_fig
{
	char const* name;
	unsigned int nloops;
	unsigned int count[16];
	double pts[2 * 400];
	double w[64];               /* per span of the first loop; 0 = 1 */
} sk_fig;

static void skf_add(sk_fig* f, double const* xy, unsigned int n)
{
	unsigned int i, off = 0;
	for (i = 0; i < f->nloops; i++) off += f->count[i];
	memcpy(f->pts + 2 * off, xy, sizeof(double) * 2 * n);
	f->count[f->nloops++] = n;
}

static qaws_scalar skf_weight(void* user, unsigned int path, unsigned int curve, unsigned int span)
{
	sk_fig const* f = (sk_fig const*)user;
	(void)curve;
	if (path == 0 && span < 64 && f->w[span] > 0)
		return (qaws_scalar)f->w[span];
	return 1;
}

static void skf_draw(svg* s, viewport const* v, sk_fig* f)
{
	qaws_curve* cs[16];
	qaws_curve const* cv[16];
	qaws_path_2d p[16];
	qaws_skeleton_desc d;
	qaws_skeleton* sk = NULL;
	unsigned int l, i, k, off = 0;
	double tmax;
	for (l = 0; l < f->nloops; l++)
	{
		qaws_scalar buf[2 * 400];
		for (i = 0; i < 2 * f->count[l]; i++) buf[i] = (qaws_scalar)f->pts[2 * off + i];
		qaws_curve_create_polyline_2d(&cs[l], buf, f->count[l], 1);
		cv[l] = cs[l];
		p[l].curves = &cv[l]; p[l].curve_count = 1; p[l].closed = 1;
		off += f->count[l];
	}
	memset(&d, 0, sizeof(d));
	d.paths = p; d.path_count = f->nloops; d.fill_rule = QAWS_FILL_NON_ZERO;
	d.weight_fn = skf_weight; d.weight_user = f;
	if (qaws_skeleton_create(&d, &sk) != QAWS_STATUS_OK)
	{
		svg_text(s, v->x0 + v->w / 2, v->y0 + v->h / 2, 12, "#cf222e", "middle", "failed");
		for (l = 0; l < f->nloops; l++) qaws_curve_destroy(cs[l]);
		return;
	}
	tmax = (double)qaws_skeleton_get_max_time(sk);
	/* the faces, lit as a roof: the face rises along its inward normal at
	   slope w, light from the upper left */
	for (i = 0; i < qaws_skeleton_get_face_count(sk); i++)
	{
		qaws_skeleton_face fc;
		qaws_skeleton_node a, b;
		double nx, ny, len, w, nz, lit;
		int g;
		qaws_skeleton_get_face(sk, i, &fc);
		qaws_skeleton_get_node(sk, fc.nodes[0], &a);
		qaws_skeleton_get_node(sk, fc.nodes[1], &b);
		len = hypot((double)b.position.x - a.position.x, (double)b.position.y - a.position.y);
		nx = -((double)b.position.y - a.position.y) / len;
		ny = ((double)b.position.x - a.position.x) / len;
		w = (double)fc.weight;
		/* the roof's normal: (-w n, 1) normalised */
		nz = 1 / sqrt(1 + w * w);
		lit = (-w * nx * -0.55 + -w * ny * 0.55 + 1 * 0.63) * nz;
		if (lit < 0) lit = 0;
		g = (int)(150 + 100 * lit);
		if (g > 250) g = 250;
		fprintf(s->f, "<path fill=\"rgb(%d,%d,%d)\" stroke=\"#ffffff\" stroke-width=\"0.4\" d=\"", g, g, (int)(g * 0.92 + 18));
		for (k = 0; k < fc.node_count; k++)
		{
			qaws_skeleton_node nd;
			qaws_skeleton_get_node(sk, fc.nodes[k], &nd);
			fprintf(s->f, "%s%.2f,%.2f ", k ? "L" : "M", vx(v, nd.position.x), vy(v, nd.position.y));
		}
		fprintf(s->f, "Z\"/>\n");
	}
	/* the wavefront at eight times */
	for (k = 1; k <= 8; k++)
	{
		qaws_clip_result* r = NULL;
		char col[32];
		double t = tmax * k / 9;
		heat((double)k / 9, col);
		if (qaws_skeleton_offset(sk, (qaws_scalar)t, &r) == QAWS_STATUS_OK)
		{
			unsigned int pi, n = qaws_clip_result_get_path_count(r);
			fprintf(s->f, "<path fill=\"none\" stroke=\"%s\" stroke-width=\"1\" d=\"", col);
			for (pi = 0; pi < n; pi++)
			{
				qaws_path_2d q;
				unsigned int ci;
				qaws_clip_result_get_path(&q, r, pi);
				for (ci = 0; ci < q.curve_count; ci++)
				{
					qaws_scalar buf[2 * 512];
					unsigned int m = 0, j;
					qaws_curve_get_control_points(buf, &m, q.curves[ci], 512);
					for (j = 0; j < m && j < 512; j++)
						fprintf(s->f, "%s%.2f,%.2f ", (ci == 0 && j == 0) ? "M" : "L", vx(v, buf[2 * j]), vy(v, buf[2 * j + 1]));
				}
				fprintf(s->f, "Z ");
			}
			fprintf(s->f, "\"/>\n");
		}
		qaws_clip_result_destroy(r);
	}
	/* the skeleton */
	fprintf(s->f, "<g stroke=\"#1b1f24\" stroke-width=\"1.3\" stroke-linecap=\"round\">");
	for (i = 0; i < qaws_skeleton_get_arc_count(sk); i++)
	{
		qaws_skeleton_arc a;
		qaws_skeleton_node x, y;
		qaws_skeleton_get_arc(sk, i, &a);
		qaws_skeleton_get_node(sk, a.node_a, &x);
		qaws_skeleton_get_node(sk, a.node_b, &y);
		fprintf(s->f, "<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\"/>", vx(v, x.position.x), vy(v, x.position.y),
			vx(v, y.position.x), vy(v, y.position.y));
	}
	fprintf(s->f, "</g>\n");
	/* the region's outline: the faces' edges */
	fprintf(s->f, "<g stroke=\"#0969da\" stroke-width=\"1.6\" stroke-linecap=\"round\">");
	for (i = 0; i < qaws_skeleton_get_face_count(sk); i++)
	{
		qaws_skeleton_face fc;
		qaws_skeleton_node a, b;
		qaws_skeleton_get_face(sk, i, &fc);
		qaws_skeleton_get_node(sk, fc.nodes[0], &a);
		qaws_skeleton_get_node(sk, fc.nodes[1], &b);
		fprintf(s->f, "<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\"/>", vx(v, a.position.x), vy(v, a.position.y),
			vx(v, b.position.x), vy(v, b.position.y));
	}
	fprintf(s->f, "</g>\n");
	{
		char info[96];
		sprintf(info, "%u nodes, %u arcs, top at t = %.3f", qaws_skeleton_get_node_count(sk), qaws_skeleton_get_arc_count(sk), tmax);
		svg_text(s, v->x0 + v->w / 2, v->y0 + v->h + 14, 11, "#57606a", "middle", info);
	}
	qaws_skeleton_destroy(sk);
	for (l = 0; l < f->nloops; l++) qaws_curve_destroy(cs[l]);
}

static void skf_bounds(sk_fig const* f, double* box)
{
	unsigned int i, n = 0;
	for (i = 0; i < f->nloops; i++) n += f->count[i];
	box[0] = box[2] = 1e300; box[1] = box[3] = -1e300;
	for (i = 0; i < n; i++)
	{
		if (f->pts[2 * i] < box[0]) box[0] = f->pts[2 * i];
		if (f->pts[2 * i] > box[1]) box[1] = f->pts[2 * i];
		if (f->pts[2 * i + 1] < box[2]) box[2] = f->pts[2 * i + 1];
		if (f->pts[2 * i + 1] > box[3]) box[3] = f->pts[2 * i + 1];
	}
}

static void demo_skeleton(void)
{
	static sk_fig figs[8];
	static double const plan[2 * 16] = { 0, 0, 6, 0, 6, 2, 9, 2, 9, 6, 7, 6, 7, 8, 3, 8, 3, 5, 1, 5, 1, 4, 0, 4, 0, 2.5, -1.5, 2.5, -1.5, 1, 0, 1 };
	static double const plus[24] = { 1, 0, 2, 0, 2, 1, 3, 1, 3, 2, 2, 2, 2, 3, 1, 3, 1, 2, 0, 2, 0, 1, 1, 1 };
	static double const comb[2 * 20] = { 0, 0, 9, 0, 9, 4, 8, 4, 8, 1.2, 6.8, 1.2, 6.8, 4, 5.6, 4, 5.6, 1.6, 3.4, 1.6, 3.4, 4, 2.2, 4, 2.2, 1.2, 1, 1.2, 1, 4, 0, 4 };
	static double const spiral[2 * 14] = { 0, 0, 8, 0, 8, 7, 2, 7, 2, 3, 5, 3, 5, 4.2, 3.2, 4.2, 3.2, 5.8, 6.8, 5.8, 6.8, 1.2, 1.2, 1.2, 1.2, 7, 0, 7 };
	static double const house[2 * 6] = { 0, 0, 10, 0, 10, 4, 6, 4, 6, 7, 0, 7 };
	unsigned int i, k;
	unsigned long long st = 2024;
	svg s;
	double pw = 360, ph = 300;
	clock_t t0 = clock();
	memset(figs, 0, sizeof(figs));
	figs[0].name = "floor plan (orthogonal, 16 edges)"; skf_add(&figs[0], plan, 16);
	figs[1].name = "plus: every event degenerate"; skf_add(&figs[1], plus, 12);
	figs[2].name = "comb: split events"; skf_add(&figs[2], comb, 16);
	figs[3].name = "orthogonal spiral: corridors of equal width"; skf_add(&figs[3], spiral, 14);
	{
		double xy[2 * 14], hole[2 * 5];
		for (i = 0; i < 14; i++)
		{
			double a = 2 * PI * i / 14 + PI / 2, r = (i & 1) ? 0.5 : 1.0;
			xy[2 * i] = r * cos(a); xy[2 * i + 1] = r * sin(a);
		}
		for (i = 0; i < 5; i++)
		{
			double a = -2 * PI * i / 5 + PI / 2;
			hole[2 * i] = 0.16 * cos(a); hole[2 * i + 1] = 0.16 * sin(a) - 0.02;
		}
		figs[4].name = "7-point star with a pentagonal hole"; skf_add(&figs[4], xy, 14); skf_add(&figs[4], hole, 5);
	}
	{
		double xy[2 * 64];
		for (i = 0; i < 64; i++)
		{
			double a = 2 * PI * i / 64, r;
			st = st * 6364136223846793005ull + 1442695040888963407ull;
			r = 1.0 + 0.55 * ((st >> 11) * (1.0 / 9007199254740992.0));
			xy[2 * i] = r * cos(a); xy[2 * i + 1] = r * sin(a);
		}
		figs[5].name = "random star-shaped polygon, 64 edges"; skf_add(&figs[5], xy, 64);
	}
	figs[6].name = "union of 9 grid rectangles";
	for (k = 0; k < 9; k++)
	{
		double x0, y0, w, h, xy[8];
		st = st * 6364136223846793005ull + 1442695040888963407ull; x0 = (double)((st >> 33) % 9);
		st = st * 6364136223846793005ull + 1442695040888963407ull; y0 = (double)((st >> 33) % 7);
		st = st * 6364136223846793005ull + 1442695040888963407ull; w = 1 + (double)((st >> 33) % 5);
		st = st * 6364136223846793005ull + 1442695040888963407ull; h = 1 + (double)((st >> 33) % 4);
		xy[0] = x0; xy[1] = y0; xy[2] = x0 + w; xy[3] = y0; xy[4] = x0 + w; xy[5] = y0 + h; xy[6] = x0; xy[7] = y0 + h;
		skf_add(&figs[6], xy, 4);
	}
	figs[7].name = "weights: end walls at speed 0.25 (near gables)";
	skf_add(&figs[7], house, 6);
	figs[7].w[1] = 0.25; figs[7].w[5] = 0.25; figs[7].w[3] = 0.25;
	if (!svg_open(&s, "showcase/b2d12_skeleton.svg", (int)(40 + 4 * (pw + 12)), (int)(100 + 2 * (ph + 44)),
		"Straight skeletons: roofs, wavefronts and degenerate events",
		"qaws_skeleton: faces lit as roofs, the skeleton (black), the wavefront at eight times coloured by time"))
		return;
	for (i = 0; i < 8; i++)
	{
		viewport v;
		double box[4], bx, by, sc;
		unsigned int col = i % 4, row = i / 4;
		skf_bounds(&figs[i], box);
		bx = box[1] - box[0]; by = box[3] - box[2];
		sc = (pw - 20) / bx < (ph - 30) / by ? (pw - 20) / bx : (ph - 30) / by;
		v.w = bx * sc; v.h = by * sc;
		v.x0 = 24 + col * (pw + 12) + (pw - v.w) / 2;
		v.y0 = 96 + row * (ph + 44) + 10 + (ph - 30 - v.h) / 2;
		v.xmin = box[0]; v.xmax = box[1]; v.ymin = box[2]; v.ymax = box[3];
		svg_text(&s, 24 + col * (pw + 12) + pw / 2, 96 + row * (ph + 44), 13, "#24292f", "middle", figs[i].name);
		skf_draw(&s, &v, &figs[i]);
	}
	svg_close(&s);
	printf("figure 12: straight skeletons, %.1f s\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
}
