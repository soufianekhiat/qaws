/*
 * Figure 22: spatial queries (qaws_mesh_bvh), on a torus soup with three
 * holes punched in it.
 *   - ray traced through the BVH: a camera ray per pixel (watertight, both
 *     sides), a shadow ray toward the light (any hit); the outside blue, the
 *     inside seen through the holes orange;
 *   - its generalized winding number on the plane z = 0 (through two of the
 *     holes): 1 in the tube, 0 outside, the 1/2 level (dark) closing the
 *     holes smoothly;
 *   - its signed distance on the same plane (negative where the winding
 *     number passes 1/2), isolines every 0.1.
 */

#define MS_BV_R 1.0
#define MS_BV_r 0.42

/* the torus as arrays, faces near the three hole centres dropped (u, v of each centre) */
static unsigned int ms_bvh_torus(qaws_scalar** ox, unsigned int* onv, unsigned int** of)
{
	static double const hole[3][2] = { { 0.0, 0.0 }, { 2.1, MS_PI }, { 4.2, 0.5 * MS_PI } };
	unsigned int nu = 120, nv = 48, i, j, k = 0, h;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nu * nv);
	unsigned int* f = (unsigned int*)malloc(sizeof(unsigned int) * 6 * nu * nv);
	for (j = 0; j < nv; j++)
		for (i = 0; i < nu; i++)
		{
			double u = 2 * MS_PI * i / nu, v = 2 * MS_PI * j / nv;
			x[3 * (j * nu + i)] = (qaws_scalar)((MS_BV_R + MS_BV_r * cos(v)) * cos(u));
			x[3 * (j * nu + i) + 1] = (qaws_scalar)((MS_BV_R + MS_BV_r * cos(v)) * sin(u));
			x[3 * (j * nu + i) + 2] = (qaws_scalar)(MS_BV_r * sin(v));
		}
	for (j = 0; j < nv; j++)
		for (i = 0; i < nu; i++)
		{
			unsigned int a = j * nu + i, b = j * nu + (i + 1) % nu, c = ((j + 1) % nv) * nu + (i + 1) % nu, d = ((j + 1) % nv) * nu + i, drop = 0;
			double u = 2 * MS_PI * (i + 0.5) / nu, v = 2 * MS_PI * (j + 0.5) / nv;
			for (h = 0; h < 3; h++)
			{
				double du = fabs(u - hole[h][0]), dv = fabs(v - hole[h][1]);
				du = fmin(du, 2 * MS_PI - du) * MS_BV_R;
				dv = fmin(dv, 2 * MS_PI - dv) * MS_BV_r;
				if (du * du + dv * dv < 0.3 * 0.3)
					drop = 1;
			}
			if (drop)
				continue;
			f[k++] = a; f[k++] = b; f[k++] = c;
			f[k++] = a; f[k++] = c; f[k++] = d;
		}
	*ox = x;
	*onv = nu * nv;
	*of = f;
	return k / 3;
}

/* the field over [-L, L]^2 at z = 0 into a panel, with its level lines (step; 0: the level 1/2 only) */
static void ms_bvh_slice(er_canvas* c, int px, int py, int pw, int ph, double const* field, int gn, double L, int winding, double step)
{
	int X, Y, side = (pw < ph ? pw : ph) * c->ss, x0 = px * c->ss + (pw * c->ss - side) / 2, y0 = py * c->ss + (ph * c->ss - side) / 2;
	for (Y = 0; Y < side; Y++)
		for (X = 0; X < side; X++)
		{
			int gx = X * gn / side, gy = (side - 1 - Y) * gn / side, line = 0;
			double v = field[gy * gn + gx], rgb[3];
			unsigned char* o = c->rgb + 3 * ((size_t)(y0 + Y) * c->W + x0 + X);
			if (gx + 1 < gn && gy + 1 < gn)
			{
				double a = v, b = field[gy * gn + gx + 1], d = field[(gy + 1) * gn + gx];
				if (winding)
					line = (a - 0.5) * (b - 0.5) < 0 || (a - 0.5) * (d - 0.5) < 0;
				else
					line = floor(a / step) != floor(b / step) || floor(a / step) != floor(d / step);
				if (!winding && ((a < 0) != (b < 0) || (a < 0) != (d < 0)))
					line = 2;
			}
			if (winding)
				ms_colormap(v, rgb);
			else
				ms_colormap(0.5 - v / 0.7, rgb);
			if (line)
			{
				double k = line == 2 || winding ? 0.15 : 0.55;
				rgb[0] *= k;
				rgb[1] *= k;
				rgb[2] *= k;
			}
			o[0] = (unsigned char)rgb[0];
			o[1] = (unsigned char)rgb[1];
			o[2] = (unsigned char)rgb[2];
		}
	(void)L;
}

static void demo_bvh(void)
{
	er_canvas* c = er_create(2200, 700, 2, 252, 252, 250);
	qaws_scalar* x;
	unsigned int nv, nt, *f, i, X, Y, gn = 360;
	double L = 1.6, t0, t1;
	qaws_mesh_bvh_desc d;
	qaws_mesh_bvh* b = NULL;
	char sub[200];
	if (!c)
		return;
	nt = ms_bvh_torus(&x, &nv, &f);
	memset(&d, 0, sizeof d);
	d.positions = x;
	d.vertex_count = nv;
	d.dimension = 3;
	d.simplices = f;
	d.simplex_size = 3;
	d.simplex_count = nt;
	t0 = (double)clock();
	qaws_mesh_bvh_create(&b, &d);
	t1 = (double)clock();
	if (!b)
		return;
	printf("  %u triangles, the BVH built in %.3f s\n", nt, (t1 - t0) / CLOCKS_PER_SEC);
	/* the ray traced panel: the canvas's own view, one ray per internal pixel */
	{
		double box[24], key[3];
		int pw, ph, ox, oy;
		unsigned int n, k;
		qaws_scalar *o, *dir, *so, *sd;
		qaws_mesh_hit *h, *sh;
		for (k = 0; k < 8; k++)
		{
			box[3 * k] = (k & 1) ? 1.5 : -1.5;
			box[3 * k + 1] = (k & 2) ? 1.5 : -1.5;
			box[3 * k + 2] = (k & 4) ? 0.45 : -0.45;
		}
		er_panel(c, 0, 72, 733, 628, box, 8, 0.5, 0.75, 0.02);
		pw = c->pw;
		ph = c->ph;
		ox = c->px;
		oy = c->py;
		n = (unsigned int)(pw * ph);
		o = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n);
		dir = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n);
		so = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n);
		sd = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n);
		h = (qaws_mesh_hit*)malloc(sizeof(qaws_mesh_hit) * n);
		sh = (qaws_mesh_hit*)malloc(sizeof(qaws_mesh_hit) * n);
		/* the key light's direction in the world: the view's key turned back */
		for (k = 0; k < 3; k++)
			key[k] = c->rot[k] * c->key[0] + c->rot[3 + k] * c->key[1] + c->rot[6 + k] * c->key[2];
		for (Y = 0; Y < (unsigned int)ph; Y++)
			for (X = 0; X < (unsigned int)pw; X++)
			{
				double q[3] = { c->cx + (X + 0.5 - pw * 0.5) / c->scale, c->cy - (Y + 0.5 - ph * 0.5) / c->scale, 10 };
				unsigned int p = Y * pw + X;
				for (k = 0; k < 3; k++)
				{
					o[3 * p + k] = (qaws_scalar)(c->rot[k] * q[0] + c->rot[3 + k] * q[1] + c->rot[6 + k] * q[2]);
					dir[3 * p + k] = (qaws_scalar)(-c->rot[6 + k]);
				}
			}
		t0 = (double)clock();
		qaws_mesh_bvh_ray(h, b, o, dir, n, 0, 0, NULL);
		for (i = 0; i < n; i++)
			for (k = 0; k < 3; k++)
			{
				so[3 * i + k] = (qaws_scalar)(h[i].point[k] + 1e-4 * key[k]);
				sd[3 * i + k] = (qaws_scalar)key[k];
			}
		qaws_mesh_bvh_ray(sh, b, so, sd, n, 0, 1, NULL);
		t1 = (double)clock();
		for (i = 0; i < n; i++)
		{
			unsigned int const* tv = h[i].vertices;
			double e1[3], e2[3], nw[3], nview[3], len, col[3], out[3], facing, q[3], rxy;
			unsigned char* px;
			if (h[i].element == QAWS_MESH_NONE)
				continue;
			for (k = 0; k < 3; k++)
			{
				e1[k] = x[3 * tv[1] + k] - x[3 * tv[0] + k];
				e2[k] = x[3 * tv[2] + k] - x[3 * tv[0] + k];
				q[k] = h[i].point[k];
			}
			/* the side seen from the face's own normal, the shading from the torus's (smooth) */
			nw[0] = e1[1] * e2[2] - e1[2] * e2[1];
			nw[1] = e1[2] * e2[0] - e1[0] * e2[2];
			nw[2] = e1[0] * e2[1] - e1[1] * e2[0];
			er_to_view(c, nw, nview);
			facing = nview[2];
			rxy = sqrt(q[0] * q[0] + q[1] * q[1]);
			nw[0] = q[0] - MS_BV_R * q[0] / rxy;
			nw[1] = q[1] - MS_BV_R * q[1] / rxy;
			nw[2] = q[2];
			len = sqrt(nw[0] * nw[0] + nw[1] * nw[1] + nw[2] * nw[2]);
			er_to_view(c, nw, nview);
			for (k = 0; k < 3; k++)
				nview[k] = (facing < 0 ? -nview[k] : nview[k]) / len;
			if (facing >= 0)
			{
				col[0] = 120; col[1] = 150; col[2] = 205;
			}
			else
			{
				col[0] = 235; col[1] = 150; col[2] = 70;
			}
			er_shade(c, nview, col, out);
			if (sh[i].element != QAWS_MESH_NONE)
				for (k = 0; k < 3; k++)
					out[k] *= 0.55;
			px = c->rgb + 3 * ((size_t)(oy + (int)(i / pw)) * c->W + ox + (int)(i % pw));
			px[0] = (unsigned char)out[0];
			px[1] = (unsigned char)out[1];
			px[2] = (unsigned char)out[2];
		}
		sprintf(sub, "%u triangles; %.1f M camera, %.1f M shadow rays: %.2f s", nt, n * 1e-6, n * 1e-6, (t1 - t0) / CLOCKS_PER_SEC);
		er_text(c, 14, 12, "ray traced through the BVH", 3, 30, 30, 30);
		er_text(c, 14, 46, sub, 2, 95, 95, 95);
		printf("  rays: %s\n", sub);
		free(o); free(dir); free(so); free(sd); free(h); free(sh);
	}
	/* the slices: the winding number, then the signed distance */
	{
		qaws_scalar* p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * gn * gn);
		qaws_scalar* v = (qaws_scalar*)malloc(sizeof(qaws_scalar) * gn * gn);
		double* field = (double*)malloc(sizeof(double) * gn * gn);
		for (Y = 0; Y < gn; Y++)
			for (X = 0; X < gn; X++)
			{
				p[3 * (Y * gn + X)] = (qaws_scalar)(-L + 2 * L * (X + 0.5) / gn);
				p[3 * (Y * gn + X) + 1] = (qaws_scalar)(-L + 2 * L * (Y + 0.5) / gn);
				p[3 * (Y * gn + X) + 2] = 0;
			}
		t0 = (double)clock();
		qaws_mesh_bvh_winding_number(v, b, p, gn * gn, NULL);
		t1 = (double)clock();
		for (i = 0; i < gn * gn; i++)
			field[i] = v[i];
		ms_bvh_slice(c, 733, 72, 733, 628, field, (int)gn, L, 1, 0);
		sprintf(sub, "%u points in %.2f s, the 1/2 level dark", gn * gn, (t1 - t0) / CLOCKS_PER_SEC);
		er_text(c, 747, 12, "winding number at z = 0", 3, 30, 30, 30);
		er_text(c, 747, 46, sub, 2, 95, 95, 95);
		printf("  winding: %s\n", sub);
		t0 = (double)clock();
		qaws_mesh_bvh_signed_distance(v, NULL, b, p, gn * gn, NULL);
		t1 = (double)clock();
		for (i = 0; i < gn * gn; i++)
			field[i] = v[i];
		ms_bvh_slice(c, 1466, 72, 733, 628, field, (int)gn, L, 0, 0.1);
		sprintf(sub, "%u points in %.2f s, isolines every 0.1, zero dark", gn * gn, (t1 - t0) / CLOCKS_PER_SEC);
		er_text(c, 1480, 12, "signed distance at z = 0", 3, 30, 30, 30);
		er_text(c, 1480, 46, sub, 2, 95, 95, 95);
		printf("  signed distance: %s\n", sub);
		free(p); free(v); free(field);
	}
	er_save(c, "showcase/mesh22_bvh.png");
	printf("mesh22_bvh.png\n");
	qaws_mesh_bvh_destroy(b);
	free(x);
	free(f);
	er_destroy(c);
}
