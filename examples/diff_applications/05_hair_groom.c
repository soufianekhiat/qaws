/* ================================================================== */
/*  5. Hair grooming on a head: 3D strands from scalp roots           */
/* ================================================================== */

#define HG_STRANDS 1200
#define HG_CP 7
#define HG_SAMPLES 24
#define HG_ITERS 160

static double const g_head[3] = { 1.0, 1.15, 1.2 };   /* ellipsoid semi-axes */

static double g_groom_align = 1.0;   /* weight of the combing guide */
static double g_groom_len = 1.0;     /* strand length scale */
static double g_groom_bend = 0.02;   /* bending weight */
static double g_photo_mask_weight = 2.0; /* hair-mask attraction in the photo term */
static double g_groom_grav = 0.6;    /* gravity weight */
static double g_groom_root_tilt = 0; /* 0: roots along the normal, 1: along the guide */

static void hg_add(qaws_vec3* v, int c, double x)
{
	if (c == 0) v->x += (qaws_scalar)x; else if (c == 1) v->y += (qaws_scalar)x; else v->z += (qaws_scalar)x;
}

typedef struct groom_strand
{
	double root[3], normal[3];
	double rest_len;
	qaws_scalar cps[HG_CP * 3];
	qaws_scalar init[HG_CP * 3];
	double color[3];
} groom_strand;

/* Ellipsoid "radius" sqrt(phi) and its gradient. */
static double head_phi(double const* x, double* g)
{
	double q = 0, r;
	int a;
	for (a = 0; a < 3; a++)
		q += x[a] * x[a] / (g_head[a] * g_head[a]);
	r = sqrt(q) + 1e-12;
	if (g)
		for (a = 0; a < 3; a++)
			g[a] = x[a] / (g_head[a] * g_head[a]) / r;
	return r;
}

/* Combing guide: down and back, parted at the center line. */
static void groom_guide(double const* x, double* g)
{
	double side = x[0] > 0 ? 1 : -1, len;
	g[0] = 0.25 * side * (0.35 + fabs(x[0]));
	g[1] = -0.55;
	g[2] = -1.0;
	len = sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
	g[0] /= len; g[1] /= len; g[2] /= len;
}

static qaws_curve* groom_curve(qaws_scalar const* cps)
{
	qaws_bspline_desc d;
	qaws_curve* c = NULL;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_3D;
	d.degree = 3;
	d.control_points = cps;
	d.control_point_count = HG_CP;
	d.is_uniform = 1;
	qaws_curve_create_bspline(&d, &c);
	return c;
}

/* Energy of one strand and its gradient on the control points:
   (L - L0)^2 + bending + gravity + head collision + guide alignment of
   the unit tangent. */
static double groom_energy(groom_strand const* s, qaws_scalar const* cps, qaws_scalar* grad, double* parts)
{
	qaws_curve* c = groom_curve(cps);
	qaws_range r = qaws_curve_get_parameter_range(c);
	qaws_scalar ts[HG_SAMPLES];
	qaws_curve_jet_3d prim[HG_SAMPLES], tan[HG_SAMPLES], bar[HG_SAMPLES];
	qaws_field_view fv;
	qaws_diff_views views = one_field(&fv, QAWS_FIELD_CONTROL_POINTS, grad, HG_CP, 3);
	double e_len, e_bend, e_grav = 0, e_col = 0, e_align = 0;
	double w_len = 20, w_bend = g_groom_bend, w_grav = g_groom_grav, w_col = 400, w_align = g_groom_align;
	int k, a;

	memset(grad, 0, sizeof(qaws_scalar) * HG_CP * 3);
	for (k = 0; k < HG_SAMPLES; k++)
		ts[k] = (qaws_scalar)(r.min_value + (r.max_value - r.min_value) * (k + 0.5) / HG_SAMPLES);
	qaws_curve_eval_batch_tangent_3d(NULL, c, ts, NULL, HG_SAMPLES, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2,
		NULL, prim, tan);
	memset(bar, 0, sizeof(bar));
	for (k = 0; k < HG_SAMPLES; k++)
	{
		double x[3] = { prim[k].d[0].x, prim[k].d[0].y, prim[k].d[0].z }, gphi[3], g[3], rr, pen, dot;
		qaws_curve_geometry_3d geo, gb;
		qaws_diff_validity validity;
		/* gravity: mean height */
		e_grav += x[2] / HG_SAMPLES;
		bar[k].d[0].z += (qaws_scalar)(w_grav / HG_SAMPLES);
		/* collision: keep outside 1.05 times the head */
		rr = head_phi(x, gphi);
		pen = 1.05 - rr;
		if (pen > 0)
		{
			e_col += pen * pen / HG_SAMPLES;
			for (a = 0; a < 3; a++)
				hg_add(&bar[k].d[0], a, -(qaws_scalar)(w_col * 2 * pen * gphi[a] / HG_SAMPLES));
		}
		/* alignment of the unit tangent with the guide (outer half of the strand) */
		qaws_curve_geometry_eval_3d(&prim[k], NULL, NULL, &geo, NULL, NULL, &validity);
		groom_guide(x, g);
		{
			double wk = (k + 0.5) / HG_SAMPLES;
			dot = geo.tangent.x * g[0] + geo.tangent.y * g[1] + geo.tangent.z * g[2];
			e_align += wk * (1 - dot) / HG_SAMPLES;
			memset(&gb, 0, sizeof(gb));
			gb.tangent = v3((qaws_scalar)(-w_align * wk * g[0] / HG_SAMPLES), (qaws_scalar)(-w_align * wk * g[1] / HG_SAMPLES),
				(qaws_scalar)(-w_align * wk * g[2] / HG_SAMPLES));
			qaws_curve_geometry_adjoint_3d(&prim[k], &gb, &bar[k], &validity);
		}
		bar[k].channels |= QAWS_EVAL_FLAG_POSITION;
	}
	qaws_curve_eval_batch_adjoint_3d(NULL, c, ts, HG_SAMPLES, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2,
		bar, &views, NULL);
	{
		/* length and bending from the integral functionals */
		qaws_scalar gl[HG_CP * 3], gbnd[HG_CP * 3], len = 0, bend = 0;
		qaws_field_view f1, f2;
		qaws_diff_views v1 = one_field(&f1, QAWS_FIELD_CONTROL_POINTS, gl, HG_CP, 3);
		qaws_diff_views v2 = one_field(&f2, QAWS_FIELD_CONTROL_POINTS, gbnd, HG_CP, 3);
		memset(gl, 0, sizeof(gl));
		memset(gbnd, 0, sizeof(gbnd));
		qaws_curve_functional_gradient(NULL, c, QAWS_FUNCTIONAL_LENGTH, 0, &v1, &len);
		qaws_curve_functional_gradient(NULL, c, QAWS_FUNCTIONAL_BENDING, 0, &v2, &bend);
		e_len = (len - s->rest_len) * (len - s->rest_len);
		e_bend = bend;
		for (k = 0; k < HG_CP * 3; k++)
			grad[k] += (qaws_scalar)(w_len * 2 * (len - s->rest_len) * gl[k] + w_bend * gbnd[k]);
	}
	/* roots are fixed: first two control points */
	for (k = 0; k < 6; k++)
		grad[k] = 0;
	if (parts)
	{
		parts[0] = e_len; parts[1] = e_bend; parts[2] = e_grav; parts[3] = e_col; parts[4] = e_align;
	}
	qaws_curve_destroy(c);
	return w_len * e_len + w_bend * e_bend + w_grav * e_grav + w_col * e_col + w_align * e_align;
}

static void groom_setup(groom_strand* st, int* count)
{
	int i, n = 0;
	unsigned int rng = 2024u;
	for (i = 0; i < 4 * HG_STRANDS && n < HG_STRANDS; i++)
	{
		/* Fibonacci directions, scalp region only */
		double z = 1 - 2 * (i + 0.5) / (4.0 * HG_STRANDS), rad = sqrt(1 - z * z), phi = i * 2.399963229728653;
		double d[3] = { rad * cos(phi), rad * sin(phi), z }, p[3], gp[3], len;
		int a, k;
		groom_strand* s;
		if (z < 0.05 || (d[1] > 0.45 && z < 0.55))
			continue;
		for (a = 0; a < 3; a++)
			p[a] = d[a] * g_head[a];
		head_phi(p, gp);
		len = sqrt(gp[0] * gp[0] + gp[1] * gp[1] + gp[2] * gp[2]);
		s = &st[n++];
		for (a = 0; a < 3; a++)
		{
			s->root[a] = p[a] * 1.01;
			s->normal[a] = gp[a] / len;
		}
		rng = rng * 1664525u + 1013904223u;
		s->rest_len = g_groom_len * (0.75 + 0.35 * ((rng >> 8) / 16777216.0));
		/* initial strand: straight along the root direction. Hair leaves the
		   scalp at a shallow angle: mostly along the combing guide projected
		   on the tangent plane, slightly outwards (g_groom_root_tilt = 0
		   keeps the plain normal). */
		{
			double dir[3], g[3], gn, tl;
			groom_guide(s->root, g);
			gn = g[0] * s->normal[0] + g[1] * s->normal[1] + g[2] * s->normal[2];
			for (a = 0; a < 3; a++)
				g[a] -= gn * s->normal[a];
			tl = sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
			for (a = 0; a < 3; a++)
				dir[a] = s->normal[a] * (1 - g_groom_root_tilt) + (tl > 1e-6 ? g[a] / tl : 0) * g_groom_root_tilt;
			tl = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
			for (k = 0; k < HG_CP; k++)
				for (a = 0; a < 3; a++)
					s->cps[3 * k + a] = (qaws_scalar)(s->root[a] + dir[a] / tl * s->rest_len * k / (HG_CP - 1));
		}
		memcpy(s->init, s->cps, sizeof(s->cps));
		rng = rng * 1664525u + 1013904223u;
		{
			double t = (rng >> 8) / 16777216.0;
			s->color[0] = 0.42 + 0.25 * t;
			s->color[1] = 0.26 + 0.17 * t;
			s->color[2] = 0.12 + 0.08 * t;
		}
	}
	*count = n;
}

/* Depth-sorted primitives: head quads and shaded hair segments. */
typedef struct prim3
{
	double xy[8];
	int n;              /* 2 = segment, 4 = quad */
	double depth;
	char color[32];
	double width;
} prim3;

static int prim3_cmp(void const* a, void const* b)
{
	double da = ((prim3 const*)a)->depth, db = ((prim3 const*)b)->depth;
	return da < db ? 1 : (da > db ? -1 : 0);
}

static void groom_render(svg* s, view3 const* v, groom_strand const* st, int ns, int use_init)
{
	int nh = 48, nq = nh * nh, cap = nq + ns * 40, n = 0, i, j, k;
	prim3* pr = (prim3*)malloc(sizeof(prim3) * (size_t)cap);
	double light[3] = { 0.35, 0.75, 0.55 }, eye[3], ln;
	double cyw = cos(v->yaw), syw = sin(v->yaw), cp = cos(v->pitch), sp = sin(v->pitch);
	ln = sqrt(light[0] * light[0] + light[1] * light[1] + light[2] * light[2]);
	for (k = 0; k < 3; k++) light[k] /= ln;
	/* view direction (towards the viewer) in world coordinates */
	eye[0] = -syw * cp; eye[1] = -cyw * cp; eye[2] = sp;
	/* head */
	for (i = 0; i < nh; i++)
		for (j = 0; j < nh; j++)
		{
			double p[4][3], dsum = 0, nrm[3], d, cen[3] = { 0, 0, 0 };
			int c;
			double th[2] = { PI * i / nh, PI * (i + 1) / nh }, ph[2] = { 2 * PI * j / nh, 2 * PI * (j + 1) / nh };
			int idx[4][2] = { { 0, 0 }, { 1, 0 }, { 1, 1 }, { 0, 1 } };
			for (c = 0; c < 4; c++)
			{
				double t = th[idx[c][0]], f = ph[idx[c][1]];
				p[c][0] = g_head[0] * sin(t) * cos(f);
				p[c][1] = g_head[1] * sin(t) * sin(f);
				p[c][2] = g_head[2] * cos(t);
				view_xform(v, p[c], &pr[n].xy[2 * c], &pr[n].xy[2 * c + 1], &d);
				dsum += d;
				for (k = 0; k < 3; k++) cen[k] += p[c][k] / 4;
			}
			head_phi(cen, nrm);
			ln = sqrt(nrm[0] * nrm[0] + nrm[1] * nrm[1] + nrm[2] * nrm[2]);
			for (k = 0; k < 3; k++) nrm[k] /= ln;
			if (nrm[0] * eye[0] + nrm[1] * eye[1] + nrm[2] * eye[2] < -0.05)
				continue;
			{
				double dif = nrm[0] * light[0] + nrm[1] * light[1] + nrm[2] * light[2], kk = 0.35 + 0.65 * (dif > 0 ? dif : 0);
				sprintf(pr[n].color, "rgb(%d,%d,%d)", (int)(240 * kk), (int)(200 * kk), (int)(178 * kk));
			}
			pr[n].n = 4;
			pr[n].depth = dsum / 4;
			pr[n].width = 0.4;
			n++;
		}
	/* strands: Kajiya-Kay style shading from the tangent */
	for (i = 0; i < ns; i++)
	{
		qaws_curve* c = groom_curve(use_init ? st[i].init : st[i].cps);
		qaws_range r = qaws_curve_get_parameter_range(c);
		double prev[3], ps[2], pd = 0;
		int m = 28;
		for (k = 0; k <= m && n < cap; k++)
		{
			qaws_eval_result_3d e;
			double x[3], sx, sy, d;
			qaws_curve_evaluate_3d(c, (qaws_scalar)(r.min_value + (r.max_value - r.min_value) * k / m), QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &e);
			x[0] = e.position.x; x[1] = e.position.y; x[2] = e.position.z;
			view_xform(v, x, &sx, &sy, &d);
			if (k > 0)
			{
				double t[3] = { e.d1.x, e.d1.y, e.d1.z }, tl = sqrt(t[0] * t[0] + t[1] * t[1] + t[2] * t[2]) + 1e-12;
				double tdl, h[3], tdh, dif, spec, kk;
				for (j = 0; j < 3; j++) { t[j] /= tl; h[j] = light[j] + eye[j]; }
				ln = sqrt(h[0] * h[0] + h[1] * h[1] + h[2] * h[2]);
				tdl = t[0] * light[0] + t[1] * light[1] + t[2] * light[2];
				tdh = (t[0] * h[0] + t[1] * h[1] + t[2] * h[2]) / ln;
				dif = sqrt(1 - tdl * tdl);
				spec = pow(sqrt(fabs(1 - tdh * tdh)), 60.0);
				kk = 0.30 + 0.70 * dif;
				sprintf(pr[n].color, "rgb(%d,%d,%d)",
					(int)fmin(255, 255 * (st[i].color[0] * kk + 0.45 * spec)),
					(int)fmin(255, 255 * (st[i].color[1] * kk + 0.38 * spec)),
					(int)fmin(255, 255 * (st[i].color[2] * kk + 0.28 * spec)));
				pr[n].n = 2;
				pr[n].xy[0] = ps[0]; pr[n].xy[1] = ps[1]; pr[n].xy[2] = sx; pr[n].xy[3] = sy;
				pr[n].depth = 0.5 * (pd + d) - 0.01;
				pr[n].width = 1.3;
				n++;
			}
			ps[0] = sx; ps[1] = sy; pd = d;
			prev[0] = x[0]; prev[1] = x[1]; prev[2] = x[2];
		}
		(void)prev;
		qaws_curve_destroy(c);
	}
	qsort(pr, (size_t)n, sizeof(prim3), prim3_cmp);
	for (i = 0; i < n; i++)
	{
		if (pr[i].n == 4)
			fprintf(s->f, "<polygon points=\"%.1f,%.1f %.1f,%.1f %.1f,%.1f %.1f,%.1f\" fill=\"%s\" stroke=\"%s\" stroke-width=\"0.4\"/>\n",
				pr[i].xy[0], pr[i].xy[1], pr[i].xy[2], pr[i].xy[3], pr[i].xy[4], pr[i].xy[5], pr[i].xy[6], pr[i].xy[7], pr[i].color, pr[i].color);
		else
			fprintf(s->f, "<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"%s\" stroke-width=\"%.2f\" stroke-linecap=\"round\"/>\n",
				pr[i].xy[0], pr[i].xy[1], pr[i].xy[2], pr[i].xy[3], pr[i].color, pr[i].width);
	}
	free(pr);
}

static void app_hair_groom(void)
{
	static groom_strand st[HG_STRANDS];
	int ns, i, it;
	double e0[5] = { 0 }, e1[5] = { 0 }, total0 = 0, total1 = 0;
	svg s;
	char buf[256];

	groom_setup(st, &ns);
	for (i = 0; i < ns; i++)
	{
		qaws_scalar g[HG_CP * 3];
		double parts[5];
		adam opt;
		int k;
		total0 += groom_energy(&st[i], st[i].cps, g, parts);
		for (k = 0; k < 5; k++) e0[k] += parts[k] / ns;
		memset(&opt, 0, sizeof(opt));
		for (it = 0; it < HG_ITERS; it++)
		{
			groom_energy(&st[i], st[i].cps, g, NULL);
			adam_step(&opt, st[i].cps, g, HG_CP * 3, 0.03 * (1.0 - 0.8 * it / (double)HG_ITERS));
		}
		total1 += groom_energy(&st[i], st[i].cps, g, parts);
		for (k = 0; k < 5; k++) e1[k] += parts[k] / ns;
	}

	svg_open(&s, "showcase/app5_hair_groom.svg", 1240, 620, "Hair grooming: 3D strands from scalp roots",
		"Cubic B-spline strands with fixed roots; energy = length + bending (integral functionals) + gravity + head "
		"collision + guide alignment of the unit tangent (3D geometry adjoint).");
	{
		view3 va = { 220, 330, 100, PI + 0.55, 0.22 }, vb = { 620, 330, 100, PI + 0.55, 0.22 }, vc = { 1020, 330, 100, 0.45, 0.22 };
		fprintf(s.f, "<rect x=\"20\" y=\"80\" width=\"1200\" height=\"520\" fill=\"#f3f1ee\"/>\n");
		fprintf(s.f, "<text x=\"30\" y=\"100\" font-size=\"13\" font-weight=\"600\" fill=\"#24292f\">initial: %d strands along the scalp normals</text>\n", ns);
		groom_render(&s, &va, st, ns, 1);
		fprintf(s.f, "<text x=\"440\" y=\"100\" font-size=\"13\" font-weight=\"600\" fill=\"#24292f\">groomed (front three-quarter)</text>\n");
		groom_render(&s, &vb, st, ns, 0);
		fprintf(s.f, "<text x=\"850\" y=\"100\" font-size=\"13\" font-weight=\"600\" fill=\"#24292f\">groomed (back)</text>\n");
		groom_render(&s, &vc, st, ns, 0);
		sprintf(buf, "per strand: length error %.4f -> %.4f, collision %.4f -> %.5f, alignment %.3f -> %.3f, energy %.2f -> %.2f",
			e0[0], e1[0], e0[3], e1[3], e0[4], e1[4], total0 / ns, total1 / ns);
		svg_text(&s, 30, 550, 14, "#24292f", "start", buf);
		svg_text(&s, 30, 576, 13, "#57606a", "start",
			"gradients: qaws_curve_functional_gradient (LENGTH, BENDING), unit tangent -> qaws_curve_geometry_adjoint_3d, "
			"positions -> jet adjoint -> qaws_curve_eval_batch_adjoint_3d");
		svg_text(&s, 30, 596, 13, "#57606a", "start",
			"shading: Kajiya-Kay style diffuse sqrt(1 - (T.L)^2) and highlight from the strand tangent; roots fixed by zeroing two control points");
	}
	svg_close(&s);
	printf("5_hair_groom: %d strands, energy %.3f -> %.3f, collision %.4f -> %.5f, alignment %.3f -> %.3f\n",
		ns, total0 / ns, total1 / ns, e0[3], e1[3], e0[4], e1[4]);
}

