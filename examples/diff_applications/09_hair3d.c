/* ================================================================== */
/*  9. 3D hair from frontal portraits                                 */
/* ================================================================== */

#define H3_ITERS 120

typedef struct hair3d_case
{
	int spec;            /* index in g_hair_specs */
	double length;       /* strand length scale (head height ~ 2.4) */
} hair3d_case;

static void app_hair3d(void)
{
	static hair3d_case const cases[4] = { { 0, 1.7 }, { 2, 1.1 }, { 4, 2.4 }, { 5, 2.2 } };
	static groom_strand st[HG_STRANDS];
	static hair_photo hp;
	svg s;
	char buf[256], path[512];
	int c;

	svg_open(&s, "showcase/app9_hair3d.svg", 1240, 1010, "3D hair from frontal portraits",
		"Head placed from the face ellipse; 3D strands keep their priors (length, bending, gravity, collision) while their "
		"projections follow each photo's orientation field and hair mask.");
	for (c = 0; c < 4; c++)
	{
		hair_spec const* sp = &g_hair_specs[cases[c].spec];
		double x0 = 20 + c * 305, sc, a0 = 0, a1 = 0, hair_mean[3];
		int ns, i, it, n0 = 0, n1 = 0;

		sprintf(path, "%s/%s.ppm", g_photos, sp->name);
		if (!image_load_ppm(path, &hp.img))
		{
			printf("9_hair3d: %s missing\n", path);
			continue;
		}
		tensor_build(&hp.img, 1.0, 2.5, &hp.tf);
		hp.mask = hair_mask(&hp.img, sp);
		blur(hp.mask, hp.img.w, hp.img.h, 2.5);
		{
			int n = hp.img.w * hp.img.h, k;
			float mx = 0;
			double hc[3] = { 0, 0, 0 }, cnt = 0;
			hp.pull = (float*)malloc(sizeof(float) * (size_t)n);
			memcpy(hp.pull, hp.mask, sizeof(float) * (size_t)n);
			blur(hp.pull, hp.img.w, hp.img.h, 10.0);
			for (k = 0; k < n; k++)
				if (hp.pull[k] > mx) mx = hp.pull[k];
			for (k = 0; k < n; k++)
			{
				hp.pull[k] = hp.pull[k] / (mx > 0 ? mx : 1);
				if (hp.mask[k] > 0.5f)
				{
					hc[0] += hp.img.rgb[3 * k]; hc[1] += hp.img.rgb[3 * k + 1]; hc[2] += hp.img.rgb[3 * k + 2];
					cnt++;
				}
			}
			for (k = 0; k < 3; k++)
				hair_mean[k] = cnt > 0 ? 0.08 + hc[k] / cnt : 0.4;
		}
		{
			/* frontal camera from the face: face half-width ~ 0.78 head units,
			   face center ~ 0.15 below the head center */
			double s_px = sp->face[2] / 0.78;
			photo_cam_make(&hp.cam, sp->face[0], sp->face[1] - 0.15 * s_px, s_px, PI, 0.0);
		}
		g_groom_len = cases[c].length;
		g_groom_align = 0.15;
		g_groom_bend = 0.1;
		g_groom_grav = 1.5;
		g_photo_mask_weight = 8.0;
		g_groom_root_tilt = 0.8;
		groom_setup(st, &ns);
		for (i = 0; i < ns; i++)
		{
			qaws_scalar g[HG_CP * 3];
			adam opt;
			double ang;
			int n_in;
			memset(&opt, 0, sizeof(opt));
			for (it = 0; it < H3_ITERS; it++)
			{
				groom_energy(&st[i], st[i].cps, g, NULL);
				adam_step(&opt, st[i].cps, g, HG_CP * 3, 0.03 * (1.0 - 0.8 * it / (double)H3_ITERS));
			}
			photo_energy(&hp, st[i].cps, NULL, &ang, &n_in);
			a0 += ang * n_in;
			n0 += n_in;
			memset(&opt, 0, sizeof(opt));
			for (it = 0; it < H3_ITERS; it++)
			{
				groom_energy(&st[i], st[i].cps, g, NULL);
				photo_energy(&hp, st[i].cps, g, NULL, NULL);
				adam_step(&opt, st[i].cps, g, HG_CP * 3, 0.02 * (1.0 - 0.8 * it / (double)H3_ITERS));
			}
			photo_energy(&hp, st[i].cps, NULL, &ang, &n_in);
			a1 += ang * n_in;
			n1 += n_in;
			{
				/* strand color: the photo under the strand's projected middle */
				qaws_curve* cv = groom_curve(st[i].cps);
				qaws_range r = qaws_curve_get_parameter_range(cv);
				qaws_eval_result_3d e;
				double x[3], px, py, d;
				int ix, iy, ch;
				for (ch = 0; ch < 3; ch++)
					st[i].color[ch] = hair_mean[ch];
				qaws_curve_evaluate_3d(cv, (qaws_scalar)(0.5 * (r.min_value + r.max_value)), QAWS_EVAL_FLAG_POSITION, &e);
				x[0] = e.position.x; x[1] = e.position.y; x[2] = e.position.z;
				view_xform(&hp.cam.v, x, &px, &py, &d);
				ix = (int)px; iy = (int)py;
				if (ix >= 0 && iy >= 0 && ix < hp.img.w && iy < hp.img.h && hp.mask[iy * hp.img.w + ix] > 0.5)
					for (ch = 0; ch < 3; ch++)
						st[i].color[ch] = 0.08 + hp.img.rgb[3 * (iy * hp.img.w + ix) + ch];
				qaws_curve_destroy(cv);
			}
		}
		g_groom_len = 1.0;
		g_groom_align = 1.0;
		g_groom_bend = 0.02;
		g_groom_grav = 0.6;
		g_photo_mask_weight = 2.0;
		g_groom_root_tilt = 0;

		/* column: photo, camera-view projection, 3D side view */
		sc = 290.0 / hp.img.w;
		svg_text(&s, x0, 92, 13, "#24292f", "start", sp->label);
		fprintf(s.f, "<image href=\"../%s/%s.png\" x=\"%.1f\" y=\"100\" width=\"290\" height=\"%.1f\"/>\n", g_photos, sp->name, x0,
			hp.img.h * sc);
		{
			double y1 = 100 + 290 * 4.0 / 3 + 12;
			view3 vf = { x0 + 145, y1 + 140, 64, PI, 0.12 }, vs = { x0 + 145, y1 + 140 + 280, 64, PI / 2 + 0.3, 0.12 };
			fprintf(s.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"290\" height=\"560\" fill=\"#f3f1ee\"/>\n", x0, y1);
			fprintf(s.f, "<clipPath id=\"h3c%d\"><rect x=\"%.1f\" y=\"%.1f\" width=\"290\" height=\"560\"/></clipPath><g clip-path=\"url(#h3c%d)\">\n", c, x0, y1, c);
			groom_render(&s, &vf, st, ns, 0);
			groom_render(&s, &vs, st, ns, 0);
			fprintf(s.f, "</g>\n");
			svg_text(&s, x0 + 6, y1 + 16, 11, "#57606a", "start", "3D strands, front");
			svg_text(&s, x0 + 6, y1 + 300, 11, "#57606a", "start", "3D strands, side");
			sprintf(buf, "%d strands, angle %.1f -> %.1f deg", ns, n0 ? a0 / n0 : 0, n1 ? a1 / n1 : 0);
			svg_text(&s, x0, y1 + 580, 12, "#0969da", "start", buf);
		}
		printf("9_hair3d: %-14s %d strands, angle to photo flow %.1f -> %.1f deg\n", sp->name, ns, n0 ? a0 / n0 : 0, n1 ? a1 / n1 : 0);
		free(hp.mask);
		free(hp.pull);
		hp.pull = NULL;
		free(hp.tf.xx); free(hp.tf.xy); free(hp.tf.yy);
		image_free(&hp.img);
	}
	svg_text(&s, 20, 1000, 11, "#57606a", "start",
		"Photos: Wikimedia Commons (see photos/CREDITS.txt). Gradients: projected unit tangent vs structure tensor and hair mask -> "
		"camera projection -> qaws_curve_geometry_adjoint_3d -> batch adjoint; priors through curve functionals.");
	svg_close(&s);
}

