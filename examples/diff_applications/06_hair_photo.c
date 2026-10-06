/* ================================================================== */
/*  6. Hair from a photo: 3D strands whose projections follow the     */
/*     photo's orientation field (single-view hair modeling)          */
/* ================================================================== */

#define HP_ITERS 150

typedef struct photo_cam
{
	view3 v;            /* maps world to image pixels */
	double P[2][3];     /* linear part of the projection */
} photo_cam;

static void photo_cam_make(photo_cam* c, double cx, double cy, double scale, double yaw, double pitch)
{
	double cyw = cos(yaw), syw = sin(yaw), cp = cos(pitch), sp = sin(pitch);
	c->v.cx = cx; c->v.cy = cy; c->v.scale = scale; c->v.yaw = yaw; c->v.pitch = pitch;
	/* view_xform: sx = cx + s (cyw x - syw y), sy = cy - s (sp (syw x + cyw y) + cp z) */
	c->P[0][0] = scale * cyw;       c->P[0][1] = -scale * syw;      c->P[0][2] = 0;
	c->P[1][0] = -scale * sp * syw; c->P[1][1] = -scale * sp * cyw; c->P[1][2] = -scale * cp;
}

typedef struct hair_photo
{
	image img;
	tensor_field tf;
	float* mask;        /* blurred hair mask, 0..1 */
	float* pull;        /* wide-range attraction to the mask (NULL: the mask) */
	photo_cam cam;
} hair_photo;

/* Image terms of one strand: orientation of the projected unit tangent
   against the structure tensor, and staying inside the hair mask.
   Adds to grad (control points); returns the energy and, optionally, the
   mean angle (degrees) over samples inside the mask. */
static double photo_energy(hair_photo const* hp, qaws_scalar const* cps, qaws_scalar* grad, double* angle, int* inside)
{
	enum { M = 24 };
	qaws_curve* c = groom_curve(cps);
	qaws_range r = qaws_curve_get_parameter_range(c);
	qaws_scalar ts[M];
	qaws_curve_jet_3d prim[M], tan[M], bar[M];
	qaws_field_view fv;
	qaws_scalar g[HG_CP * 3];
	qaws_diff_views views = one_field(&fv, QAWS_FIELD_CONTROL_POINTS, g, HG_CP, 3);
	double e = 0, w_img = 3.0, w_mask = g_photo_mask_weight, ang = 0;
	int k, a, n_in = 0;
	int W = hp->img.w, H = hp->img.h;

	for (k = 0; k < M; k++)
		ts[k] = (qaws_scalar)(r.min_value + (r.max_value - r.min_value) * (k + 0.5) / M);
	qaws_curve_eval_batch_tangent_3d(NULL, c, ts, NULL, M, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2,
		NULL, prim, tan);
	memset(bar, 0, sizeof(bar));
	for (k = 0; k < M; k++)
	{
		double x[3] = { prim[k].d[0].x, prim[k].d[0].y, prim[k].d[0].z }, px, py, depth, mdx, mdy, m;
		qaws_curve_geometry_3d geo, gb;
		qaws_diff_validity validity;
		double t2[2], len2, u[2], dex, dey, detx, dety, ea, dt2[2], dpos[2];
		view_xform(&hp->cam.v, x, &px, &py, &depth);
		if (px < 1 || py < 1 || px > W - 2 || py > H - 2)
			continue;
		/* the photo only sees the front: points behind the head's center
		   plane follow the 3D priors alone */
		if (depth > 0)
			continue;
		/* stay on the hair: (1 - pull)^2 with a wide-range pull field (the
		   sharp mask is flat inside the face and would not push strands out) */
		{
			double pl = sample(hp->pull ? hp->pull : hp->mask, W, H, px, py, &mdx, &mdy);
			e += w_mask * (1 - pl) * (1 - pl) / M;
			dpos[0] = -w_mask * 2 * (1 - pl) * mdx / M;
			dpos[1] = -w_mask * 2 * (1 - pl) * mdy / M;
		}
		m = sample(hp->mask, W, H, px, py, NULL, NULL);
		/* projected tangent against the orientation field, where there is hair */
		qaws_curve_geometry_eval_3d(&prim[k], NULL, NULL, &geo, NULL, NULL, &validity);
		t2[0] = hp->cam.P[0][0] * geo.tangent.x + hp->cam.P[0][1] * geo.tangent.y + hp->cam.P[0][2] * geo.tangent.z;
		t2[1] = hp->cam.P[1][0] * geo.tangent.x + hp->cam.P[1][1] * geo.tangent.y + hp->cam.P[1][2] * geo.tangent.z;
		len2 = sqrt(t2[0] * t2[0] + t2[1] * t2[1]);
		memset(&gb, 0, sizeof(gb));
		if (m > 0.5 && len2 > 1e-6)
		{
			double wk = m, du;
			u[0] = t2[0] / len2;
			u[1] = t2[1] / len2;
			ea = align_energy(&hp->tf, px, py, u[0], u[1], &dex, &dey, &detx, &dety);
			e += w_img * wk * ea / M;
			dpos[0] += w_img * wk * dex / M;
			dpos[1] += w_img * wk * dey / M;
			/* d/dt2 of e(t2 / |t2|) = (I - u u^T) de/du / |t2| */
			du = u[0] * detx + u[1] * dety;
			dt2[0] = w_img * wk * (detx - du * u[0]) / len2 / M;
			dt2[1] = w_img * wk * (dety - du * u[1]) / len2 / M;
			gb.tangent = v3((qaws_scalar)(hp->cam.P[0][0] * dt2[0] + hp->cam.P[1][0] * dt2[1]),
				(qaws_scalar)(hp->cam.P[0][1] * dt2[0] + hp->cam.P[1][1] * dt2[1]),
				(qaws_scalar)(hp->cam.P[0][2] * dt2[0] + hp->cam.P[1][2] * dt2[1]));
			qaws_curve_geometry_adjoint_3d(&prim[k], &gb, &bar[k], &validity);
			if (angle)
			{
				double ox, oy, d;
				tensor_direction(&hp->tf, px, py, &ox, &oy);
				d = fabs(u[0] * ox + u[1] * oy);
				ang += acos(d > 1 ? 1 : d) * 180.0 / PI;
				n_in++;
			}
		}
		/* image-plane position gradient back to 3D: P^T dpos */
		for (a = 0; a < 3; a++)
			hg_add(&bar[k].d[0], a, hp->cam.P[0][a] * dpos[0] + hp->cam.P[1][a] * dpos[1]);
		bar[k].channels |= QAWS_EVAL_FLAG_POSITION;
	}
	memset(g, 0, sizeof(g));
	qaws_curve_eval_batch_adjoint_3d(NULL, c, ts, M, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2,
		bar, &views, NULL);
	if (grad)
		for (k = 6; k < HG_CP * 3; k++)   /* roots fixed */
			grad[k] += g[k];
	if (angle)
		*angle = n_in ? ang / n_in : 0;
	if (inside)
		*inside = n_in;
	qaws_curve_destroy(c);
	return e;
}

/* Strands drawn in the camera view, over the photo. */
static void photo_overlay(svg* s, hair_photo const* hp, groom_strand const* st, int ns, int use_init, double x0, double y0, double sc,
	char const* color, double opacity)
{
	int i, k;
	for (i = 0; i < ns; i++)
	{
		qaws_curve* c = groom_curve(use_init ? st[i].init : st[i].cps);
		qaws_range r = qaws_curve_get_parameter_range(c);
		double xy[2 * 25];
		int n = 0;
		for (k = 0; k <= 24; k++)
		{
			qaws_eval_result_3d e;
			double x[3], px, py, d;
			qaws_curve_evaluate_3d(c, (qaws_scalar)(r.min_value + (r.max_value - r.min_value) * k / 24), QAWS_EVAL_FLAG_POSITION, &e);
			x[0] = e.position.x; x[1] = e.position.y; x[2] = e.position.z;
			view_xform(&hp->cam.v, x, &px, &py, &d);
			xy[2 * n] = x0 + px * sc;
			xy[2 * n + 1] = y0 + py * sc;
			n++;
		}
		svg_polyline(s, xy, n, color, 0.7, opacity, 0);
		qaws_curve_destroy(c);
	}
}

static void app_hair_photo(void)
{
	static groom_strand st[HG_STRANDS];
	static hair_photo hp;
	char path[512], buf[256], name[64];
	int ns, i, it, n_in0 = 0, n_in1 = 0;
	double a0 = 0, a1 = 0;
	svg s;

	/* a hair photo when provided, otherwise the ribbons as the flow source */
	strcpy(name, "hair");
	sprintf(path, "%s/hair.ppm", g_photos);
	if (!image_load_ppm(path, &hp.img))
	{
		strcpy(name, "ribbons");
		sprintf(path, "%s/ribbons.ppm", g_photos);
		if (!image_load_ppm(path, &hp.img))
		{
			printf("6_hair_photo: no photo (see examples/photo_to_ppm.ps1)\n");
			return;
		}
	}
	tensor_build(&hp.img, 1.0, 3.0, &hp.tf);
	{
		int n = hp.img.w * hp.img.h;
		hp.mask = (float*)malloc(sizeof(float) * (size_t)n);
		for (i = 0; i < n; i++)
			hp.mask[i] = hp.img.lum[i] > 0.12f ? 1.0f : 0.0f;
		blur(hp.mask, hp.img.w, hp.img.h, 4.0);
	}
	/* camera: the photo is the front three-quarter view of the head */
	photo_cam_make(&hp.cam, hp.img.w * 0.52, hp.img.h * 0.80, hp.img.h * 0.48, PI + 0.55, 0.22);

	/* start from the generic groom, then follow the photo */
	groom_setup(st, &ns);
	for (i = 0; i < ns; i++)
	{
		qaws_scalar g[HG_CP * 3];
		adam opt;
		double ang;
		int n_in;
		memset(&opt, 0, sizeof(opt));
		for (it = 0; it < HG_ITERS; it++)
		{
			groom_energy(&st[i], st[i].cps, g, NULL);
			adam_step(&opt, st[i].cps, g, HG_CP * 3, 0.03 * (1.0 - 0.8 * it / (double)HG_ITERS));
		}
		memcpy(st[i].init, st[i].cps, sizeof(st[i].cps));
		photo_energy(&hp, st[i].cps, NULL, &ang, &n_in);
		a0 += ang * n_in;
		n_in0 += n_in;
		memset(&opt, 0, sizeof(opt));
		for (it = 0; it < HP_ITERS; it++)
		{
			groom_energy(&st[i], st[i].cps, g, NULL);
			photo_energy(&hp, st[i].cps, g, NULL, NULL);
			adam_step(&opt, st[i].cps, g, HG_CP * 3, 0.02 * (1.0 - 0.8 * it / (double)HP_ITERS));
		}
		photo_energy(&hp, st[i].cps, NULL, &ang, &n_in);
		a1 += ang * n_in;
		n_in1 += n_in;
	}

	svg_open(&s, "showcase/app6_hair_photo.svg", 1240, 640, "Hair from a photo: projected strands follow the image",
		"Single-view hair modeling: 3D strands keep their 3D priors while their projections align with the photo's "
		"orientation field and stay on the hair mask.");
	{
		double sc = 400.0 / hp.img.w, ph = hp.img.h * sc;
		view3 vside = { 1030, 330, 100, PI / 2 + 0.25, 0.18 };
		fprintf(s.f, "<text x=\"20\" y=\"92\" font-size=\"13\" font-weight=\"600\" fill=\"#24292f\">photo (%s) and orientation field</text>\n", name);
		fprintf(s.f, "<image href=\"../%s/%s.png\" x=\"20\" y=\"100\" width=\"400\" height=\"%.1f\"/>\n", g_photos, name, ph);
		{
			int gx, gy;
			for (gy = 6; gy < hp.img.h; gy += 12)
				for (gx = 6; gx < hp.img.w; gx += 12)
				{
					double ox, oy, sx = 20 + gx * sc, sy = 100 + gy * sc;
					char col[32];
					if (hp.mask[gy * hp.img.w + gx] < 0.5)
						continue;
					heat(tensor_direction(&hp.tf, gx, gy, &ox, &oy), col);
					svg_line(&s, sx - 3.5 * ox, sy - 3.5 * oy, sx + 3.5 * ox, sy + 3.5 * oy, col, 1.2, 0.9);
				}
		}
		fprintf(s.f, "<text x=\"430\" y=\"92\" font-size=\"13\" font-weight=\"600\" fill=\"#24292f\">camera view: generic groom (orange) and photo-driven (cyan)</text>\n");
		fprintf(s.f, "<image href=\"../%s/%s.png\" x=\"430\" y=\"100\" width=\"400\" height=\"%.1f\" opacity=\"0.55\"/>\n", g_photos, name, ph);
		photo_overlay(&s, &hp, st, ns, 1, 430, 100, sc, "#ff8c42", 0.35);
		photo_overlay(&s, &hp, st, ns, 0, 430, 100, sc, "#39d0ff", 0.8);
		fprintf(s.f, "<text x=\"850\" y=\"92\" font-size=\"13\" font-weight=\"600\" fill=\"#24292f\">3D result, side view</text>\n");
		fprintf(s.f, "<rect x=\"845\" y=\"100\" width=\"380\" height=\"420\" fill=\"#f3f1ee\"/>\n");
		groom_render(&s, &vside, st, ns, 0);
		sprintf(buf, "mean angle between projected strands and the photo flow: %.1f deg (generic groom) -> %.1f deg",
			n_in0 ? a0 / n_in0 : 0, n_in1 ? a1 / n_in1 : 0);
		svg_text(&s, 20, 100 + ph + 40, 15, "#0969da", "start", buf);
		sprintf(buf, "%d strands; samples on the hair mask: %d -> %d", ns, n_in0, n_in1);
		svg_text(&s, 20, 100 + ph + 64, 13, "#24292f", "start", buf);
		svg_text(&s, 20, 100 + ph + 90, 13, "#57606a", "start",
			"image term: e(P T / |P T|) with e = u^T J u / tr J (structure tensor J), projection P of the camera; gradient through");
		svg_text(&s, 20, 100 + ph + 108, 13, "#57606a", "start",
			"(I - u u^T) / |P T|, P^T, qaws_curve_geometry_adjoint_3d and the batch adjoint; mask term (1 - m)^2 pulls strands onto the hair.");
		svg_text(&s, 20, 100 + ph + 126, 13, "#57606a", "start",
			"3D priors kept: length, bending, gravity, head collision. Put any hair photo at photos/hair.ppm (photo_to_ppm.ps1).");
	}
	svg_close(&s);
	printf("6_hair_photo: %d strands, angle %.1f -> %.1f deg, samples on mask %d -> %d\n", ns,
		n_in0 ? a0 / n_in0 : 0, n_in1 ? a1 / n_in1 : 0, n_in0, n_in1);
	free(hp.mask);
	free(hp.tf.xx); free(hp.tf.xy); free(hp.tf.yy);
	image_free(&hp.img);
}

