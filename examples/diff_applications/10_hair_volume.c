/* ================================================================== */
/*  10. Single-view hair modeling: strands lifted onto a hair volume  */
/*      and fitted in 3D with projection adjoints                      */
/* ================================================================== */

#define HV_MAXCP 20
#define HV_MAX (3 * HS_MAX_STRANDS)
#define HV_ITERS 40

/* --- strand rasterizer --------------------------------------------- */

/* PNG writer with stored (uncompressed) deflate blocks: no zlib needed. */
static unsigned int png_crc(unsigned int crc, unsigned char const* p, size_t n)
{
	static unsigned int table[256];
	static int ready = 0;
	size_t i;
	if (!ready)
	{
		unsigned int c, k, b;
		for (k = 0; k < 256; k++)
		{
			c = k;
			for (b = 0; b < 8; b++)
				c = c & 1 ? 0xedb88320u ^ (c >> 1) : c >> 1;
			table[k] = c;
		}
		ready = 1;
	}
	crc = ~crc;
	for (i = 0; i < n; i++)
		crc = table[(crc ^ p[i]) & 0xff] ^ (crc >> 8);
	return ~crc;
}

static void png_be32(unsigned char* p, unsigned int v)
{
	p[0] = (unsigned char)(v >> 24); p[1] = (unsigned char)(v >> 16); p[2] = (unsigned char)(v >> 8); p[3] = (unsigned char)v;
}

static void png_chunk(FILE* f, char const* type, unsigned char const* data, size_t n)
{
	unsigned char hdr[8];
	unsigned int crc;
	png_be32(hdr, (unsigned int)n);
	memcpy(hdr + 4, type, 4);
	fwrite(hdr, 1, 8, f);
	if (n)
		fwrite(data, 1, n, f);
	crc = png_crc(0, hdr + 4, 4);
	crc = png_crc(crc, data, n);
	png_be32(hdr, crc);
	fwrite(hdr, 1, 4, f);
}

static int png_write(char const* path, int w, int h, unsigned char const* rgb)
{
	static unsigned char const sig[8] = { 137, 80, 78, 71, 13, 10, 26, 10 };
	size_t raw_n = (size_t)h * (3 * (size_t)w + 1), blocks = (raw_n + 65534) / 65535, n = 2 + raw_n + 5 * blocks + 4, pos = 0, i;
	unsigned char* raw = (unsigned char*)malloc(raw_n);
	unsigned char* z = (unsigned char*)malloc(n);
	unsigned char ihdr[13];
	unsigned int a = 1, b = 0;
	FILE* f = fopen(path, "wb");
	int y;
	if (!f || !raw || !z)
	{
		if (f) fclose(f);
		free(raw); free(z);
		return 0;
	}
	for (y = 0; y < h; y++)
	{
		raw[(size_t)y * (3 * w + 1)] = 0;
		memcpy(raw + (size_t)y * (3 * w + 1) + 1, rgb + (size_t)y * 3 * w, 3 * (size_t)w);
	}
	z[pos++] = 0x78;
	z[pos++] = 0x01;
	for (i = 0; i < raw_n; i += 65535)
	{
		size_t len = raw_n - i < 65535 ? raw_n - i : 65535;
		z[pos++] = (unsigned char)(i + len == raw_n ? 1 : 0);
		z[pos++] = (unsigned char)(len & 0xff);
		z[pos++] = (unsigned char)(len >> 8);
		z[pos++] = (unsigned char)(~len & 0xff);
		z[pos++] = (unsigned char)((~len >> 8) & 0xff);
		memcpy(z + pos, raw + i, len);
		pos += len;
	}
	for (i = 0; i < raw_n; i++)
	{
		a = (a + raw[i]) % 65521;
		b = (b + a) % 65521;
	}
	png_be32(z + pos, (b << 16) | a);
	pos += 4;
	fwrite(sig, 1, 8, f);
	png_be32(ihdr, (unsigned int)w);
	png_be32(ihdr + 4, (unsigned int)h);
	ihdr[8] = 8; ihdr[9] = 2; ihdr[10] = 0; ihdr[11] = 0; ihdr[12] = 0;
	png_chunk(f, "IHDR", ihdr, 13);
	png_chunk(f, "IDAT", z, pos);
	png_chunk(f, "IEND", NULL, 0);
	fclose(f);
	free(raw);
	free(z);
	return 1;
}

typedef struct hv_canvas
{
	int w, h;          /* supersampled size */
	float* rgb;
	float* depth;      /* head / neck depth, 1e9 where empty */
} hv_canvas;

static void hv_canvas_init(hv_canvas* c, int w, int h, double const* bg)
{
	int i;
	c->w = w;
	c->h = h;
	c->rgb = (float*)malloc(sizeof(float) * 3 * (size_t)w * h);
	c->depth = (float*)malloc(sizeof(float) * (size_t)w * h);
	for (i = 0; i < w * h; i++)
	{
		c->rgb[3 * i] = (float)bg[0]; c->rgb[3 * i + 1] = (float)bg[1]; c->rgb[3 * i + 2] = (float)bg[2];
		c->depth[i] = 1e9f;
	}
}

/* Box-filters the ss x ss supersamples and writes a PNG. */
static void hv_canvas_save(hv_canvas* c, int ss, char const* path)
{
	int w = c->w / ss, h = c->h / ss, x, y, i, j, k;
	unsigned char* out = (unsigned char*)malloc(3 * (size_t)w * h);
	for (y = 0; y < h; y++)
		for (x = 0; x < w; x++)
			for (k = 0; k < 3; k++)
			{
				double s = 0;
				for (j = 0; j < ss; j++)
					for (i = 0; i < ss; i++)
						s += c->rgb[3 * ((y * ss + j) * c->w + x * ss + i) + k];
				s = s / (ss * ss);
				out[3 * (y * w + x) + k] = (unsigned char)(255 * (s < 0 ? 0 : (s > 1 ? 1 : pow(s, 1 / 1.1))));
			}
	png_write(path, w, h, out);
	free(out);
	free(c->rgb);
	free(c->depth);
}

/* World point of the pixel ray at view depth y2 (inverse of view_xform). */
static void hv_ray(view3 const* v, double sx, double sy, double y2, double* p)
{
	double cy = cos(v->yaw), sy_ = sin(v->yaw), cp = cos(v->pitch), sp = sin(v->pitch);
	double x = (sx - v->cx) / v->scale, z2 = (v->cy - sy) / v->scale;
	double y = cp * y2 + sp * z2, z = -sp * y2 + cp * z2;
	p[0] = cy * x + sy_ * y;
	p[1] = -sy_ * x + cy * y;
	p[2] = z;
}

/* Head ellipsoid and neck cylinder ray-cast per pixel: depth buffer and
   shading. scalp_color, when given, tints the head points scalp_test accepts. */
static void hv_canvas_head(hv_canvas* c, view3 const* v, double const* light, double const* skin,
	int (*scalp_test)(void const*, double const*), void const* scalp_ctx, double const* scalp_color)
{
	int x, y, a;
	for (y = 0; y < c->h; y++)
		for (x = 0; x < c->w; x++)
		{
			double p0[3], p1[3], d[3], best = 1e9, nrm[3], hit[3], A, B, C, disc, t;
			int which = 0;
			hv_ray(v, x + 0.5, y + 0.5, 0, p0);
			hv_ray(v, x + 0.5, y + 0.5, 1, p1);
			for (a = 0; a < 3; a++)
				d[a] = p1[a] - p0[a];
			/* ellipsoid */
			A = B = 0;
			C = -1;
			for (a = 0; a < 3; a++)
			{
				double h2 = g_head[a] * g_head[a];
				A += d[a] * d[a] / h2;
				B += 2 * p0[a] * d[a] / h2;
				C += p0[a] * p0[a] / h2;
			}
			disc = B * B - 4 * A * C;
			if (disc >= 0)
			{
				t = (-B - sqrt(disc)) / (2 * A);
				best = t;
				which = 1;
			}
			/* neck: (x / 0.42)^2 + ((y + 0.12) / 0.45)^2 = 1, -2.62 < z < -0.7 */
			A = d[0] * d[0] / (0.42 * 0.42) + d[1] * d[1] / (0.45 * 0.45);
			B = 2 * p0[0] * d[0] / (0.42 * 0.42) + 2 * (p0[1] + 0.12) * d[1] / (0.45 * 0.45);
			C = p0[0] * p0[0] / (0.42 * 0.42) + (p0[1] + 0.12) * (p0[1] + 0.12) / (0.45 * 0.45) - 1;
			disc = B * B - 4 * A * C;
			if (A > 1e-12 && disc >= 0)
			{
				t = (-B - sqrt(disc)) / (2 * A);
				if (t < best)
				{
					double z = p0[2] + t * d[2];
					if (z < -0.7 && z > -2.62)
					{
						best = t;
						which = 2;
					}
				}
			}
			if (!which)
				continue;
			for (a = 0; a < 3; a++)
				hit[a] = p0[a] + best * d[a];
			if (which == 1)
				for (a = 0; a < 3; a++)
					nrm[a] = hit[a] / (g_head[a] * g_head[a]);
			else
			{
				nrm[0] = hit[0] / (0.42 * 0.42);
				nrm[1] = (hit[1] + 0.12) / (0.45 * 0.45);
				nrm[2] = 0;
			}
			{
				double ln = sqrt(nrm[0] * nrm[0] + nrm[1] * nrm[1] + nrm[2] * nrm[2]) + 1e-12, dif, kk;
				double const* col = skin;
				dif = (nrm[0] * light[0] + nrm[1] * light[1] + nrm[2] * light[2]) / ln;
				kk = 0.3 + 0.7 * (dif > 0 ? dif : 0);
				if (which == 1 && scalp_test && scalp_test(scalp_ctx, hit))
				{
					col = scalp_color;
					kk *= 0.6;
				}
				for (a = 0; a < 3; a++)
					c->rgb[3 * (y * c->w + x) + a] = (float)(col[a] * kk);
				c->depth[y * c->w + x] = (float)best;
			}
		}
}

typedef struct hv_seg
{
	float x0, y0, x1, y1, d0, d1;
	float rgb[3];
	float depth;
} hv_seg;

static int hv_seg_cmp(void const* a, void const* b)
{
	float da = ((hv_seg const*)a)->depth, db = ((hv_seg const*)b)->depth;
	return da < db ? 1 : (da > db ? -1 : 0);
}

/* Far-to-near alpha blending of anti-aliased segments, depth tested
   against the head. */
static void hv_canvas_segments(hv_canvas* c, hv_seg* segs, int n, double width, double alpha)
{
	int i, x, y, k;
	double hw = (width > 1 ? width : 1) * 0.5, wa = (width < 1 ? width : 1) * alpha;
	qsort(segs, (size_t)n, sizeof(hv_seg), hv_seg_cmp);
	for (i = 0; i < n; i++)
	{
		hv_seg const* s = &segs[i];
		double dx = s->x1 - s->x0, dy = s->y1 - s->y0, l2 = dx * dx + dy * dy + 1e-12;
		int xa = (int)floor(fmin(s->x0, s->x1) - hw - 1), xb = (int)ceil(fmax(s->x0, s->x1) + hw + 1);
		int ya = (int)floor(fmin(s->y0, s->y1) - hw - 1), yb = (int)ceil(fmax(s->y0, s->y1) + hw + 1);
		if (xa < 0) xa = 0;
		if (ya < 0) ya = 0;
		if (xb > c->w - 1) xb = c->w - 1;
		if (yb > c->h - 1) yb = c->h - 1;
		for (y = ya; y <= yb; y++)
			for (x = xa; x <= xb; x++)
			{
				double px = x + 0.5 - s->x0, py = y + 0.5 - s->y0, t = (px * dx + py * dy) / l2, ex, ey, dist, cov, dd;
				if (t < 0) t = 0;
				if (t > 1) t = 1;
				ex = px - t * dx;
				ey = py - t * dy;
				dist = sqrt(ex * ex + ey * ey);
				cov = hw + 0.5 - dist;
				if (cov <= 0)
					continue;
				if (cov > 1) cov = 1;
				dd = s->d0 + t * (s->d1 - s->d0);
				if (dd > c->depth[y * c->w + x] + 0.02)
					continue;
				cov *= wa;
				for (k = 0; k < 3; k++)
				{
					float* p = &c->rgb[3 * (y * c->w + x) + k];
					*p = (float)(*p * (1 - cov) + s->rgb[k] * cov);
				}
			}
	}
}


typedef struct hv_cam
{
	double cx, cy, s, yaw;      /* orthographic camera, pitch 0 */
} hv_cam;

/* world -> pixel and depth toward the camera */
static void hv_project(hv_cam const* c, double const* x, double* px, double* py, double* t)
{
	double cw = cos(c->yaw), sw = sin(c->yaw);
	*px = c->cx + c->s * (cw * x[0] - sw * x[1]);
	*py = c->cy - c->s * x[2];
	*t = -(sw * x[0] + cw * x[1]);
}

static void hv_unproject(hv_cam const* c, double px, double py, double t, double* x)
{
	double cw = cos(c->yaw), sw = sin(c->yaw), u = (px - c->cx) / c->s;
	x[0] = cw * u - sw * t;
	x[1] = -sw * u - cw * t;
	x[2] = (c->cy - py) / c->s;
}

/* Depths (toward the camera) where the pixel ray meets the head; 0 if missed. */
static int hv_head_hit(hv_cam const* c, double px, double py, double* t_front, double* t_back)
{
	double p0[3], p1[3], A = 0, B = 0, C = -1, disc;
	int a;
	hv_unproject(c, px, py, 0, p0);
	hv_unproject(c, px, py, 1, p1);
	for (a = 0; a < 3; a++)
	{
		double d = p1[a] - p0[a], h2 = g_head[a] * g_head[a];
		A += d * d / h2;
		B += 2 * p0[a] * d / h2;
		C += p0[a] * p0[a] / h2;
	}
	disc = B * B - 4 * A * C;
	if (disc < 0 || A <= 0)
		return 0;
	*t_front = (-B + sqrt(disc)) / (2 * A);
	*t_back = (-B - sqrt(disc)) / (2 * A);
	return 1;
}

/* Two-pass chamfer distance (pixels) to the outside of the mask. */
static float* hv_distance(float const* mask, int w, int h)
{
	float* d = (float*)malloc(sizeof(float) * (size_t)w * h);
	int x, y;
	for (y = 0; y < h; y++)
		for (x = 0; x < w; x++)
			d[y * w + x] = mask[y * w + x] > 0.5f ? 1e6f : 0.0f;
	for (y = 0; y < h; y++)
		for (x = 0; x < w; x++)
		{
			float v = d[y * w + x];
			if (x > 0 && d[y * w + x - 1] + 1 < v) v = d[y * w + x - 1] + 1;
			if (y > 0 && d[(y - 1) * w + x] + 1 < v) v = d[(y - 1) * w + x] + 1;
			if (x > 0 && y > 0 && d[(y - 1) * w + x - 1] + 1.4142f < v) v = d[(y - 1) * w + x - 1] + 1.4142f;
			if (x + 1 < w && y > 0 && d[(y - 1) * w + x + 1] + 1.4142f < v) v = d[(y - 1) * w + x + 1] + 1.4142f;
			d[y * w + x] = v;
		}
	for (y = h - 1; y >= 0; y--)
		for (x = w - 1; x >= 0; x--)
		{
			float v = d[y * w + x];
			if (x + 1 < w && d[y * w + x + 1] + 1 < v) v = d[y * w + x + 1] + 1;
			if (y + 1 < h && d[(y + 1) * w + x] + 1 < v) v = d[(y + 1) * w + x] + 1;
			if (x + 1 < w && y + 1 < h && d[(y + 1) * w + x + 1] + 1.4142f < v) v = d[(y + 1) * w + x + 1] + 1.4142f;
			if (x > 0 && y + 1 < h && d[(y + 1) * w + x - 1] + 1.4142f < v) v = d[(y + 1) * w + x - 1] + 1.4142f;
			d[y * w + x] = v;
		}
	return d;
}

typedef struct hv_strand
{
	int ncp, layer;             /* 0 front, 1 back, 2 middle */
	qaws_scalar cps[3 * HV_MAXCP];
	qaws_scalar knots[HV_MAXCP + 4];
	double color[3];
	strand const* src;          /* traced 2D strand */
} hv_strand;

static qaws_curve* hv_curve(hv_strand const* s, qaws_scalar const* cps)
{
	qaws_bspline_desc d;
	qaws_curve* c = NULL;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_3D;
	d.degree = 3;
	d.control_points = cps;
	d.control_point_count = (unsigned int)s->ncp;
	d.knots = s->knots;
	d.knot_count = (unsigned int)s->ncp + 4;
	qaws_curve_create_bspline(&d, &c);
	return c;
}

typedef struct hv_scene
{
	image img;
	tensor_field tf;
	float* mask;
	float* surf[3];      /* depth toward the camera of the front, back and middle hair surfaces */
	double hair_col[3];  /* mean hair color */
	float* scalp;        /* the hair mask, blurred: where the scalp takes the hair color */
	int frontal;
	float* back_mask;    /* frontal photos: where hair hangs behind the head (NULL otherwise) */
	hv_cam cam;
} hv_scene;

/* E = projection fit to the traced strand + distance to the strand's hair
   surface + head collision + bending. Exact gradient through the 3D batch
   adjoint and the bending functional. */
static double hv_energy(hv_scene const* sc, hv_strand const* s, qaws_scalar const* cps, qaws_scalar* grad)
{
	enum { MAXS = HV_MAXCP * 4 + 8 };
	qaws_curve* c = hv_curve(s, cps);
	qaws_range r = qaws_curve_get_parameter_range(c);
	qaws_scalar ts[MAXS];
	qaws_curve_jet_3d prim[MAXS], tan[MAXS], bar[MAXS];
	qaws_field_view fv;
	qaws_diff_views views = one_field(&fv, QAWS_FIELD_CONTROL_POINTS, grad, (unsigned int)s->ncp, 3);
	int m = s->ncp * 4 + 8, k, a, W = sc->img.w, H = sc->img.h;
	double e = 0, cw = cos(sc->cam.yaw), sw = sin(sc->cam.yaw), S = sc->cam.s;
	double w_proj = 1.0, w_depth = 2.0, w_col = 200, w_bend = 0.02;
	float const* surf = sc->surf[s->layer];

	for (k = 0; k < m; k++)
		ts[k] = (qaws_scalar)(r.min_value + (r.max_value - r.min_value) * (k + 0.5) / m);
	qaws_curve_eval_batch_tangent_3d(NULL, c, ts, NULL, (unsigned int)m, QAWS_EVAL_FLAG_POSITION, NULL, prim, tan);
	memset(bar, 0, sizeof(qaws_curve_jet_3d) * (size_t)m);
	memset(grad, 0, sizeof(qaws_scalar) * 3 * (size_t)s->ncp);
	for (k = 0; k < m; k++)
	{
		double x[3] = { prim[k].d[0].x, prim[k].d[0].y, prim[k].d[0].z }, px, py, t, gphi[3], rr;
		double f = (k + 0.5) / m * (s->src->n - 1), wq, qx, qy, dpx[3], dpy[3], dt[3], g[3] = { 0, 0, 0 };
		int i0 = (int)f, i1;
		i1 = i0 + 1 < s->src->n ? i0 + 1 : i0;
		wq = f - i0;
		qx = s->src->xy[2 * i0] * (1 - wq) + s->src->xy[2 * i1] * wq;
		qy = s->src->xy[2 * i0 + 1] * (1 - wq) + s->src->xy[2 * i1 + 1] * wq;
		hv_project(&sc->cam, x, &px, &py, &t);
		dpx[0] = S * cw; dpx[1] = -S * sw; dpx[2] = 0;
		dpy[0] = 0; dpy[1] = 0; dpy[2] = -S;
		dt[0] = -sw; dt[1] = -cw; dt[2] = 0;
		/* projection: |(p - q) / S|^2 */
		{
			double ex = (px - qx) / S, ey = (py - qy) / S;
			e += w_proj * (ex * ex + ey * ey) / m;
			for (a = 0; a < 3; a++)
				g[a] += w_proj * 2 * (ex * dpx[a] + ey * dpy[a]) / S / m;
		}
		/* surface: (t - D(p))^2 */
		if (px > 1 && py > 1 && px < W - 2 && py < H - 2)
		{
			double gx, gy, D = sample(surf, W, H, px, py, &gx, &gy), ed = t - D;
			e += w_depth * ed * ed / m;
			for (a = 0; a < 3; a++)
				g[a] += w_depth * 2 * ed * (dt[a] - gx * dpx[a] - gy * dpy[a]) / m;
		}
		/* head collision */
		rr = head_phi(x, gphi);
		if (rr < 1.02)
		{
			double pen = 1.02 - rr;
			e += w_col * pen * pen / m;
			for (a = 0; a < 3; a++)
				g[a] -= w_col * 2 * pen * gphi[a] / m;
		}
		bar[k].d[0] = v3((qaws_scalar)g[0], (qaws_scalar)g[1], (qaws_scalar)g[2]);
		bar[k].channels = QAWS_EVAL_FLAG_POSITION;
	}
	qaws_curve_eval_batch_adjoint_3d(NULL, c, ts, (unsigned int)m, QAWS_EVAL_FLAG_POSITION, bar, &views, NULL);
	{
		/* bending, scaled by the strand length so short strands may curl */
		qaws_scalar bg[3 * HV_MAXCP], bend = 0;
		qaws_field_view bf;
		qaws_diff_views bv = one_field(&bf, QAWS_FIELD_CONTROL_POINTS, bg, (unsigned int)s->ncp, 3);
		double L = s->src->n * 1.0 / S + 1e-3, scale = w_bend / (L * L * L);
		memset(bg, 0, sizeof(bg));
		qaws_curve_functional_gradient(NULL, c, QAWS_FUNCTIONAL_BENDING, 0, &bv, &bend);
		e += scale * bend;
		for (k = 0; k < 3 * s->ncp; k++)
			grad[k] += (qaws_scalar)(scale * bg[k]);
	}
	qaws_curve_destroy(c);
	return e;
}

/* --- Gabor orientation (Chai et al. 2012) --------------------------- */

#define HV_ANGLES 32

/* Orientation and confidence of the stripes of f by a bank of HV_ANGLES
   even Gabor kernels (sigma_u 1.8, sigma_v 2.4, lambda 4); theta is the
   stripe direction, conf the spread of the responses around the maximum.
   Only pixels where roi > 0.2 are filtered. */
static void hv_gabor_pass(float const* f, float const* roi, int w, int h, float* theta, float* conf)
{
	enum { R = 8, K = 2 * R + 1 };
	static float ker[HV_ANGLES][K * K];
	static int ready = 0;
	int a, x, y, i, j;
	if (!ready)
	{
		for (a = 0; a < HV_ANGLES; a++)
		{
			double th = PI * a / HV_ANGLES, c = cos(th), s = sin(th), mean = 0;
			for (j = 0; j < K; j++)
				for (i = 0; i < K; i++)
				{
					double u = (i - R) * c + (j - R) * s, v = -(i - R) * s + (j - R) * c;
					ker[a][j * K + i] = (float)(exp(-0.5 * (u * u / (1.8 * 1.8) + v * v / (2.4 * 2.4))) * cos(2 * PI * u / 4.0));
					mean += ker[a][j * K + i];
				}
			mean /= K * K;
			for (i = 0; i < K * K; i++)
				ker[a][i] -= (float)mean;
		}
		ready = 1;
	}
	for (y = 0; y < h; y++)
		for (x = 0; x < w; x++)
		{
			double F[HV_ANGLES], best = -1, wsum = 0;
			int ab = 0;
			theta[y * w + x] = 0;
			conf[y * w + x] = 0;
			if (roi[y * w + x] < 0.2f || x < R || y < R || x >= w - R || y >= h - R)
				continue;
			for (a = 0; a < HV_ANGLES; a++)
			{
				double r = 0;
				float const* k = ker[a];
				for (j = 0; j < K; j++)
				{
					float const* row = f + (y + j - R) * w + x - R;
					for (i = 0; i < K; i++)
						r += k[j * K + i] * row[i];
				}
				F[a] = fabs(r);
				if (F[a] > best)
				{
					best = F[a];
					ab = a;
				}
			}
			for (a = 0; a < HV_ANGLES; a++)
			{
				int da = abs(a - ab);
				double d = PI * (da < HV_ANGLES - da ? da : HV_ANGLES - da) / HV_ANGLES;
				wsum += d * (F[a] - best) * (F[a] - best);
			}
			/* the kernel oscillates across the stripes: they run 90 degrees off */
			theta[y * w + x] = (float)(PI * ab / HV_ANGLES + 0.5 * PI);
			conf[y * w + x] = (float)sqrt(wsum);
		}
}

/* Builds the tracer's tensor field from the Gabor orientation: one pass on
   the luminance, one refinement pass on the confidence map (which removes
   the noise responses), confidences normalized to their 95th percentile. */
static void hv_gabor_field(image const* img, float const* mask, tensor_field* t)
{
	int w = img->w, h = img->h, n = w * h, i;
	float* roi = (float*)malloc(sizeof(float) * (size_t)n);
	float* th = (float*)malloc(sizeof(float) * (size_t)n);
	float* cf = (float*)malloc(sizeof(float) * (size_t)n);
	float* th2 = (float*)malloc(sizeof(float) * (size_t)n);
	float* cf2 = (float*)malloc(sizeof(float) * (size_t)n);
	float* tmp = (float*)malloc(sizeof(float) * (size_t)n);
	double scale;
	memcpy(roi, mask, sizeof(float) * (size_t)n);
	blur(roi, w, h, 4.0);
	hv_gabor_pass(img->lum, roi, w, h, th, cf);
	hv_gabor_pass(cf, roi, w, h, th2, cf2);
	{
		/* 95th percentile of the confidences inside the region */
		int cnt = 0, hist[1001], acc = 0, b;
		double mx = 0;
		memset(hist, 0, sizeof(hist));
		for (i = 0; i < n; i++)
			if (roi[i] > 0.2f && cf2[i] > mx)
				mx = cf2[i];
		for (i = 0; i < n; i++)
			if (roi[i] > 0.2f)
			{
				hist[(int)(1000 * cf2[i] / (mx + 1e-12))]++;
				cnt++;
			}
		for (b = 0; b <= 1000 && acc < 0.95 * cnt; b++)
			acc += hist[b];
		scale = 1.0 / (b / 1000.0 * mx + 1e-12);
	}
	t->w = w;
	t->h = h;
	t->xx = (float*)calloc((size_t)n, sizeof(float));
	t->xy = (float*)calloc((size_t)n, sizeof(float));
	t->yy = (float*)calloc((size_t)n, sizeof(float));
	for (i = 0; i < n; i++)
	{
		/* tensor_direction reads a gradient tensor: the normal of the stripes */
		double c = fmin(1.0, cf2[i] * scale), nx = -sin(th2[i]), ny = cos(th2[i]);
		t->xx[i] = (float)(c * nx * nx);
		t->xy[i] = (float)(c * nx * ny);
		t->yy[i] = (float)(c * ny * ny);
	}
	blur(t->xx, w, h, 1.2);
	blur(t->xy, w, h, 1.2);
	blur(t->yy, w, h, 1.2);
	(void)tmp;
	free(tmp); free(roi); free(th); free(cf); free(th2); free(cf2);
}

/* --- 3D orientation field and strands grown from the scalp ---------- */

#define HV_GN 100          /* grid cells per axis (x, y); z has HV_GZ */
#define HV_GZ 110
#define HV_GH 0.05         /* cell size, head units */
#define HV_GROWN_MAX 28000
#define HV_GROWN_PTS 160

typedef struct hv_volume
{
	float* occ;            /* 1 inside the hair volume */
	float* T;              /* 6 components of the orientation tensor per cell */
	unsigned char* fixed;  /* constrained cells */
	float* dens;           /* strand points per cell, for self-shadowing */
} hv_volume;

static double const g_hv_lo[3] = { -0.5 * HV_GN * HV_GH, -0.5 * HV_GN * HV_GH, -3.6 };

static int hv_cell(double const* p, int* i, int* j, int* k)
{
	*i = (int)floor((p[0] - g_hv_lo[0]) / HV_GH);
	*j = (int)floor((p[1] - g_hv_lo[1]) / HV_GH);
	*k = (int)floor((p[2] - g_hv_lo[2]) / HV_GH);
	return *i >= 0 && *j >= 0 && *k >= 0 && *i < HV_GN && *j < HV_GN && *k < HV_GZ;
}

#define HV_IDX(i, j, k) (((size_t)(k) * HV_GN + (size_t)(j)) * HV_GN + (size_t)(i))

/* Trilinear read of a scalar grid (or of component c of a 6-vector grid). */
static double hv_trilinear(float const* g, int stride, int c, double const* p)
{
	double f[3], r = 0;
	int b[3], a, di, dj, dk;
	for (a = 0; a < 3; a++)
	{
		double u = (p[a] - g_hv_lo[a]) / HV_GH - 0.5;
		b[a] = (int)floor(u);
		f[a] = u - b[a];
	}
	for (dk = 0; dk < 2; dk++)
		for (dj = 0; dj < 2; dj++)
			for (di = 0; di < 2; di++)
			{
				int i = b[0] + di, j = b[1] + dj, k = b[2] + dk;
				double wgt = (di ? f[0] : 1 - f[0]) * (dj ? f[1] : 1 - f[1]) * (dk ? f[2] : 1 - f[2]);
				if (i < 0 || j < 0 || k < 0 || i >= HV_GN || j >= HV_GN || k >= HV_GZ)
					continue;
				r += wgt * g[stride * HV_IDX(i, j, k) + c];
			}
	return r;
}

/* Dominant direction of the field at p, signed to agree with prev. */
static double hv_field_dir(hv_volume const* v, double const* p, double const* prev, double* out)
{
	double T[6], q[3], r[3], ln;
	int c, it;
	for (c = 0; c < 6; c++)
		T[c] = hv_trilinear(v->T, 6, c, p);
	q[0] = prev[0]; q[1] = prev[1]; q[2] = prev[2];
	for (it = 0; it < 3; it++)
	{
		r[0] = T[0] * q[0] + T[1] * q[1] + T[2] * q[2];
		r[1] = T[1] * q[0] + T[3] * q[1] + T[4] * q[2];
		r[2] = T[2] * q[0] + T[4] * q[1] + T[5] * q[2];
		ln = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
		if (ln < 1e-9)
		{
			out[0] = prev[0]; out[1] = prev[1]; out[2] = prev[2];
			return 0;
		}
		q[0] = r[0] / ln; q[1] = r[1] / ln; q[2] = r[2] / ln;
	}
	if (q[0] * prev[0] + q[1] * prev[1] + q[2] * prev[2] < 0)
	{
		q[0] = -q[0]; q[1] = -q[1]; q[2] = -q[2];
	}
	out[0] = q[0]; out[1] = q[1]; out[2] = q[2];
	return T[0] + T[3] + T[5];
}

/* Adds the unit tangent t as t t^T to the constraint accumulators of the cell at p. */
static void hv_splat(hv_volume* v, float* acc, float* wacc, double const* p, double const* t)
{
	int i, j, k;
	size_t id;
	if (!hv_cell(p, &i, &j, &k))
		return;
	id = HV_IDX(i, j, k);
	acc[6 * id + 0] += (float)(t[0] * t[0]);
	acc[6 * id + 1] += (float)(t[0] * t[1]);
	acc[6 * id + 2] += (float)(t[0] * t[2]);
	acc[6 * id + 3] += (float)(t[1] * t[1]);
	acc[6 * id + 4] += (float)(t[1] * t[2]);
	acc[6 * id + 5] += (float)(t[2] * t[2]);
	wacc[id] += 1;
	(void)v;
}

typedef struct hv_grown
{
	int n;
	float p[3 * HV_GROWN_PTS];
	float rgb[3 * HV_GROWN_PTS];
} hv_grown;

/* Under the hair, the scalp takes the hair color: the head point projects
   into the photo's hair mask, or (frontal photos) lies on the hidden back
   half above the nape. */
static int hv_scalp(hv_scene const* sc, double const* x)
{
	double px, py, t;
	int ix, iy;
	hv_project(&sc->cam, x, &px, &py, &t);
	if (sc->frontal && t < 0 && x[2] > -0.45)
		return 1;
	ix = (int)px;
	iy = (int)py;
	if (ix < 0 || iy < 0 || ix >= sc->img.w || iy >= sc->img.h)
		return 0;
	return sc->scalp[iy * sc->img.w + ix] > 0.25f;
}

static int g_hv_flat_dark = 0;   /* reject flat dark pixels (a black jacket next to black hair) */

/* Hair mask for a plain background: the background is flooded from the
   border pixels near the border's median color (small steps, bounded drift
   from the starting color); the remaining pixels are hair when their color
   is nearer a hair seed than a skin / clothes seed, inside the box, off the
   face, connected to the hair seeds; a closing fills the holes. */
static float* hv_mask(image const* img, hair_spec const* sp)
{
	int w = img->w, h = img->h, n = w * h, i, k, q, top = 0;
	float* avg = (float*)malloc(sizeof(float) * 3 * (size_t)n);
	float* m = (float*)calloc((size_t)n, sizeof(float));
	unsigned char* lab = (unsigned char*)calloc((size_t)n, 1);   /* 1 background, 2 hair, 3 other, 4 hair candidate */
	int* stack = (int*)malloc(sizeof(int) * (size_t)n);
	int* origin = (int*)malloc(sizeof(int) * (size_t)n);
	double hc[8][3], nc[8][3], bgc[3];
	int hs[8][2];
	float *lm = NULL, *lv = NULL;
	for (k = 0; k < 3; k++)
	{
		float* ch = (float*)malloc(sizeof(float) * (size_t)n);
		for (i = 0; i < n; i++)
			ch[i] = img->rgb[3 * i + k];
		blur(ch, w, h, 1.2);
		for (i = 0; i < n; i++)
			avg[3 * i + k] = ch[i];
		free(ch);
	}
#define HV_D2(a, b) ((avg[3 * (a)] - avg[3 * (b)]) * (avg[3 * (a)] - avg[3 * (b)]) + \
	(avg[3 * (a) + 1] - avg[3 * (b) + 1]) * (avg[3 * (a) + 1] - avg[3 * (b) + 1]) + \
	(avg[3 * (a) + 2] - avg[3 * (b) + 2]) * (avg[3 * (a) + 2] - avg[3 * (b) + 2]))
	{
		/* the background color: per-channel median of the border; hair and
		   clothes touching the border must not seed the flood */
		int nb = 2 * (w + h), cnt[3][256], c;
		memset(cnt, 0, sizeof(cnt));
		for (i = 0; i < n; i++)
		{
			int x = i % w, y = i / w;
			if (x == 0 || y == 0 || x == w - 1 || y == h - 1)
				for (q = 0; q < 3; q++)
					cnt[q][(int)(255 * fmin(1, fmax(0, avg[3 * i + q])))]++;
		}
		for (q = 0; q < 3; q++)
		{
			int acc = 0;
			for (c = 0; c < 256 && acc < nb / 2; c++)
				acc += cnt[q][c];
			bgc[q] = c / 255.0;
		}
	}
	for (i = 0; i < n; i++)
	{
		int x = i % w, y = i / w;
		double db = 0;
		for (q = 0; q < 3; q++)
			db += (avg[3 * i + q] - bgc[q]) * (avg[3 * i + q] - bgc[q]);
		if ((x == 0 || y == 0 || x == w - 1 || y == h - 1) && db < 0.02)
		{
			lab[i] = 1;
			origin[i] = i;
			stack[top++] = i;
		}
	}
	while (top > 0)
	{
		int p = stack[--top], x = p % w, y = p / w;
		int nb[4] = { x > 0 ? p - 1 : -1, x < w - 1 ? p + 1 : -1, y > 0 ? p - w : -1, y < h - 1 ? p + w : -1 };
		for (q = 0; q < 4; q++)
		{
			int r = nb[q];
			if (r < 0 || lab[r])
				continue;
			if (HV_D2(r, p) < 0.0012 && HV_D2(r, origin[p]) < 0.03)
			{
				lab[r] = 1;
				origin[r] = origin[p];
				stack[top++] = r;
			}
		}
	}
#undef HV_D2
	/* hair seeds that landed on the background move to the nearby pixel
	   farthest from the background color */
	for (k = 0; k < sp->nhair; k++)
	{
		int bx = sp->hair[k][0], by = sp->hair[k][1], dx, dy;
		double best = -1;
		hs[k][0] = bx;
		hs[k][1] = by;
		if (lab[by * w + bx] != 1)
			continue;
		for (dy = -10; dy <= 10; dy++)
			for (dx = -10; dx <= 10; dx++)
			{
				int x = bx + dx, y = by + dy;
				double db = 0;
				if (x < 0 || y < 0 || x >= w || y >= h || lab[y * w + x] == 1)
					continue;
				for (q = 0; q < 3; q++)
					db += (avg[3 * (y * w + x) + q] - bgc[q]) * (avg[3 * (y * w + x) + q] - bgc[q]);
				if (db > best)
				{
					best = db;
					hs[k][0] = x;
					hs[k][1] = y;
				}
			}
	}
	for (k = 0; k < sp->nhair; k++)
		for (q = 0; q < 3; q++)
			hc[k][q] = avg[3 * (hs[k][1] * w + hs[k][0]) + q];
	for (k = 0; k < sp->nnon; k++)
		for (q = 0; q < 3; q++)
			nc[k][q] = avg[3 * (sp->non[k][1] * w + sp->non[k][0]) + q];
	if (g_hv_flat_dark)
	{
		/* local luminance mean and variance: flat dark cloth is not hair */
		lm = (float*)malloc(sizeof(float) * (size_t)n);
		lv = (float*)malloc(sizeof(float) * (size_t)n);
		for (i = 0; i < n; i++)
		{
			lm[i] = (float)((avg[3 * i] + avg[3 * i + 1] + avg[3 * i + 2]) / 3);
			lv[i] = lm[i] * lm[i];
		}
		blur(lm, w, h, 2.5);
		blur(lv, w, h, 2.5);
	}
	for (i = 0; i < n; i++)
	{
		int x = i % w, y = i / w;
		double dh = 1e9, dn = 1e9, fx, fy;
		if (lab[i] == 1)
			continue;
		lab[i] = 3;
		if (x < sp->box[0] || x > sp->box[2] || y < sp->box[1] || y > sp->box[3])
			continue;
		fx = (x - sp->face[0]) / (double)sp->face[2];
		fy = (y - sp->face[1]) / (double)sp->face[3];
		if (fx * fx + fy * fy < 1)
			continue;
		for (k = 0; k < sp->nhair; k++)
		{
			double d = 0;
			for (q = 0; q < 3; q++)
				d += (avg[3 * i + q] - hc[k][q]) * (avg[3 * i + q] - hc[k][q]);
			if (d < dh) dh = d;
		}
		for (k = 0; k < sp->nnon; k++)
		{
			double d = 0;
			for (q = 0; q < 3; q++)
				d += (avg[3 * i + q] - nc[k][q]) * (avg[3 * i + q] - nc[k][q]);
			if (d < dn) dn = d;
		}
		if (dh < dn && !(lm && lm[i] < 0.18 && lv[i] - lm[i] * lm[i] < 0.0005))
			lab[i] = 4;
	}
	/* keep the candidates connected to a hair seed */
	top = 0;
	for (k = 0; k < sp->nhair; k++)
	{
		int p = hs[k][1] * w + hs[k][0];
		if (lab[p] == 4)
		{
			lab[p] = 2;
			stack[top++] = p;
		}
	}
	while (top > 0)
	{
		int p = stack[--top], x = p % w, y = p / w;
		int nb[4] = { x > 0 ? p - 1 : -1, x < w - 1 ? p + 1 : -1, y > 0 ? p - w : -1, y < h - 1 ? p + w : -1 };
		for (q = 0; q < 4; q++)
			if (nb[q] >= 0 && lab[nb[q]] == 4)
			{
				lab[nb[q]] = 2;
				stack[top++] = nb[q];
			}
	}
	/* closing: fill the holes of the hair region (highlights, dark gaps),
	   never into the background */
	for (i = 0; i < n; i++)
		m[i] = lab[i] == 2 ? 1.0f : 0.0f;
	blur(m, w, h, 5.0);
	for (i = 0; i < n; i++)
		m[i] = (m[i] > 0.3f && lab[i] != 1) || lab[i] == 2 ? 1.0f : 0.0f;
	blur(m, w, h, 1.0);
	free(avg);
	free(lab);
	free(stack);
	free(origin);
	free(lm);
	free(lv);
	return m;
}

typedef struct hv_case
{
	hair_spec spec;
	double yaw;          /* camera yaw: PI = frontal, PI / 2 = left profile */
	double cx, cy, s;    /* head center (pixels) and pixels per head unit */
	int flat_dark;       /* black clothes touch the hair */
	double tie[3];       /* ponytail tie: pixel x, y and depth toward the camera (x = 0: none) */
	double curl[2];      /* curl period and radius, head units (0: straight); set per photo */
} hv_case;

/* Builds the hair volume and the 3D strands of one portrait. */
static int hv_reconstruct(hv_scene* sc, hv_case const* hc, strand* traced, int* ntraced, hv_strand* hs, double* proj_rms)
{
	hair_spec const* sp = &hc->spec;
	int W = sc->img.w, H = sc->img.h, x, y, i, n = 0, it, nt, nreal, l;
	float* dist;
	double s_px = hc->s, R = 0.6, err = 0, hair_col[3] = { 0, 0, 0 }, wsum = 0;
	int cnt = 0, frontal = fabs(hc->yaw - PI) < 0.6;
	sc->cam.cx = hc->cx;
	sc->cam.cy = hc->cy;
	sc->cam.s = s_px;
	sc->cam.yaw = hc->yaw;
	g_hv_flat_dark = hc->flat_dark;
	sc->mask = hv_mask(&sc->img, sp);
	hv_gabor_field(&sc->img, sc->mask, &sc->tf);
	sc->back_mask = NULL;
	dist = hv_distance(sc->mask, W, H);
	for (l = 0; l < 3; l++)
		sc->surf[l] = (float*)malloc(sizeof(float) * (size_t)W * H);
	for (y = 0; y < H; y++)
		for (x = 0; x < W; x++)
		{
			/* front: the head (plus a hair layer) or the silhouette inflated with
			   a circular profile of radius R. back: behind the head; below it,
			   in a frontal view, a curtain hanging down the back, closing in on
			   the neck as it falls. Hair over the scalp is 0.1 thick. */
			double d = dist[y * W + x] / s_px, infl, tf, tb, f, b, u, z;
			if (d > R) d = R;
			infl = sqrt(fmax(0, 2 * R * d - d * d));
			u = (x - hc->cx) / s_px;
			z = (hc->cy - y) / s_px;
			f = infl;
			b = -infl;
			if (frontal && z < 0.3)
			{
				double depth = 0.9 * g_head[1] * sqrt(fmax(0, 1 - u * u / 2.6));
				if (z < -0.6)
					depth = 0.7 + (depth - 0.7) * exp(1.5 * (z + 0.6));
				b = fmin(b, -depth);
			}
			if (hv_head_hit(&sc->cam, x, y, &tf, &tb))
			{
				f = fmax(tf + 0.1, f);
				b = fmin(tb - 0.1, b);
			}
			sc->surf[0][y * W + x] = (float)f;
			sc->surf[1][y * W + x] = (float)b;
			sc->surf[2][y * W + x] = (float)(0.5 * (f + b));
			if (sc->mask[y * W + x] > 0.5)
			{
				for (l = 0; l < 3; l++)
					hair_col[l] += sc->img.rgb[3 * (y * W + x) + l];
				wsum++;
			}
		}
	for (l = 0; l < 3; l++)
	{
		blur(sc->surf[l], W, H, 2.0);
		hair_col[l] /= wsum + 1e-9;
		sc->hair_col[l] = hair_col[l];
	}
	sc->frontal = frontal;
	sc->scalp = (float*)malloc(sizeof(float) * (size_t)W * H);
	memcpy(sc->scalp, sc->mask, sizeof(float) * (size_t)W * H);
	blur(sc->scalp, W, H, 6.0);
	free(dist);
	{
		image masked = sc->img;
		masked.lum = sc->mask;
		nt = trace_strands(&masked, &sc->tf, 1.0, 0.12, 0.5, traced, HS_MAX_STRANDS);
	}
	nreal = nt;
	/* The face hides the back of the head: in a frontal view, comb synthetic
	   strands down the back, from the top of the hair in each column to the
	   hair's lower end on the left and right of the face. */
	if (frontal)
	{
		int bl = 0, br = 0, x0 = (int)(hc->cx - 0.95 * s_px), x1 = (int)(hc->cx + 0.95 * s_px);
		unsigned int rng = 77u;
		sc->back_mask = (float*)calloc((size_t)W * H, sizeof(float));
		for (y = 0; y < H; y++)
			for (x = 0; x < W; x++)
				if (sc->mask[y * W + x] > 0.5)
				{
					if (x < hc->cx - 0.6 * s_px && y > bl) bl = y;
					if (x > hc->cx + 0.6 * s_px && y > br) br = y;
				}
		for (x = x0 < 1 ? 1 : x0; x <= x1 && x < W - 1 && nt < HS_MAX_STRANDS; x += 2)
		{
			int top = -1, bot, k, np;
			double f = (x - x0) / (double)(x1 - x0 + 1), amp, freq, phase;
			for (y = 0; y < H && top < 0; y++)
				if (sc->mask[y * W + x] > 0.5)
					top = y;
			if (top < 0)
				continue;
			bot = (int)(bl * (1 - f) + br * f);
			if (bot > H - 2) bot = H - 2;
			for (y = top; y <= bot; y++)
				for (k = x - 1; k <= x + 1; k++)
					if (k >= 0 && k < W)
						sc->back_mask[y * W + k] = 1;
			/* a random sway so the strands do not run parallel */
			rng = rng * 1664525u + 1013904223u;
			amp = 1.0 + 3.0 * ((rng >> 8) / 16777216.0);
			rng = rng * 1664525u + 1013904223u;
			freq = 0.02 + 0.05 * ((rng >> 8) / 16777216.0);
			phase = (rng & 1023) * 0.01;
			top += (int)((rng >> 12) % 6);
			np = bot - top;
			if (np < 8)
				continue;
			traced[nt].n = np;
			traced[nt].xy = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * (size_t)np);
			traced[nt].fit = traced[nt].opt = NULL;
			for (k = 0; k < np; k++)
			{
				traced[nt].xy[2 * k] = (qaws_scalar)(x + amp * sin(freq * k + phase));
				traced[nt].xy[2 * k + 1] = (qaws_scalar)(top + k);
			}
			nt++;
		}
	}
	*ntraced = nt;
	/* lift, fit and optimize: front, back and middle copies of the traced
	   strands (middle only off the head), back copies of the synthetic ones */
	for (i = 0; i < nt; i++)
	{
		strand const* st = &traced[i];
		int layer;
		for (layer = 0; layer < 3 && n < HV_MAX; layer++)
		{
			hv_strand* h = &hs[n];
			qaws_scalar* pts;
			qaws_bspline_fit_desc d;
			qaws_curve* c = NULL;
			unsigned int got = 0;
			int k, mid = st->n / 2;
			double tf, tb;
			adam opt;
			if (i >= nreal && layer != 1)
				continue;
			if (layer == 2 && hv_head_hit(&sc->cam, st->xy[2 * mid], st->xy[2 * mid + 1], &tf, &tb))
				continue;
			pts = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (size_t)st->n);
			for (k = 0; k < st->n; k++)
			{
				double p[3], px = st->xy[2 * k], py = st->xy[2 * k + 1];
				double t = sample(sc->surf[layer], W, H, px, py, NULL, NULL);
				hv_unproject(&sc->cam, px, py, t, p);
				pts[3 * k] = (qaws_scalar)p[0]; pts[3 * k + 1] = (qaws_scalar)p[1]; pts[3 * k + 2] = (qaws_scalar)p[2];
			}
			h->ncp = st->n / 14 + 3;
			if (h->ncp < 4) h->ncp = 4;
			if (h->ncp > HV_MAXCP) h->ncp = HV_MAXCP;
			h->layer = layer;
			h->src = st;
			memset(&d, 0, sizeof(d));
			d.dimension = QAWS_DIMENSION_3D;
			d.data_points = pts;
			d.data_point_count = (unsigned int)st->n;
			d.degree = 3;
			d.control_point_count = (unsigned int)h->ncp;
			if (qaws_curve_fit_bspline(&d, &c) != QAWS_STATUS_OK)
			{
				free(pts);
				continue;
			}
			qaws_curve_read_field(c, QAWS_FIELD_CONTROL_POINTS, h->cps, 3 * h->ncp, &got);
			qaws_curve_read_field(c, QAWS_FIELD_KNOTS, h->knots, h->ncp + 4, &got);
			qaws_curve_destroy(c);
			free(pts);
			memset(&opt, 0, sizeof(opt));
			for (it = 0; it < HV_ITERS; it++)
			{
				qaws_scalar g[3 * HV_MAXCP];
				hv_energy(sc, h, h->cps, g);
				adam_step(&opt, h->cps, g, 3 * h->ncp, 0.01 * (1.0 - 0.7 * it / (double)HV_ITERS));
			}
			{
				/* color: the photo under the strand's hair pixels; residual projection error */
				double col[3] = { 0, 0, 0 }, cw = 0, shade = layer == 0 ? 1.0 : (layer == 1 ? 0.7 : 0.8);
				int q;
				for (k = 0; k < st->n; k++)
				{
					int ix = (int)st->xy[2 * k], iy = (int)st->xy[2 * k + 1], j = iy * W + ix;
					if (ix < 0 || iy < 0 || ix >= W || iy >= H || sc->mask[j] < 0.8)
						continue;
					for (q = 0; q < 3; q++)
						col[q] += sc->img.rgb[3 * j + q];
					cw++;
				}
				for (q = 0; q < 3; q++)
					h->color[q] = shade * (cw > 3 ? col[q] / cw : hair_col[q]);
				if (layer == 0)
				{
					qaws_curve* cc = hv_curve(h, h->cps);
					qaws_range rr = qaws_curve_get_parameter_range(cc);
					for (k = 0; k < 16; k++)
					{
						qaws_eval_result_3d e;
						double p[3], px, py, t, f = (k + 0.5) / 16 * (st->n - 1);
						int i0 = (int)f;
						qaws_curve_evaluate_3d(cc, (qaws_scalar)(rr.min_value + (rr.max_value - rr.min_value) * (k + 0.5) / 16), QAWS_EVAL_FLAG_POSITION, &e);
						p[0] = e.position.x; p[1] = e.position.y; p[2] = e.position.z;
						hv_project(&sc->cam, p, &px, &py, &t);
						err += (px - st->xy[2 * i0]) * (px - st->xy[2 * i0]) + (py - st->xy[2 * i0 + 1]) * (py - st->xy[2 * i0 + 1]);
						cnt++;
					}
					qaws_curve_destroy(cc);
				}
			}
			n++;
		}
	}
	if (getenv("QAWS_HV_DEBUG"))
	{
		/* the hair mask in red over the photo */
		char name[256];
		FILE* f;
		sprintf(name, "showcase/hv_mask_%s.ppm", sp->name);
		f = fopen(name, "wb");
		if (f)
		{
			fprintf(f, "P6\n%d %d\n255\n", W, H);
			for (i = 0; i < W * H; i++)
				for (l = 0; l < 3; l++)
					fputc((int)(255 * (l == 0 ? 0.5 * sc->img.rgb[3 * i] + 0.5 * sc->mask[i] : sc->img.rgb[3 * i + l] * (1 - 0.5 * sc->mask[i]))), f);
			fclose(f);
		}
	}
	*proj_rms = cnt ? sqrt(err / cnt) : 0;
	return n;
}

static void hv_free(hv_scene* sc, strand* traced, int nt)
{
	int i;
	for (i = 0; i < nt; i++)
		free(traced[i].xy);
	free(sc->back_mask);
	free(sc->mask); free(sc->scalp); free(sc->surf[0]); free(sc->surf[1]); free(sc->surf[2]);
	free(sc->tf.xx); free(sc->tf.xy); free(sc->tf.yy);
	image_free(&sc->img);
}

static hv_case const g_hv_cases[5] = {
	{ { "plain_wavy", "long, wavy",
		{ { 85, 120 }, { 90, 250 }, { 260, 200 }, { 270, 280 }, { 150, 60 }, { 200, 55 }, { 70, 330 }, { 290, 340 } }, 8,
		{ { 175, 160 }, { 175, 205 }, { 150, 330 }, { 180, 400 }, { 330, 400 }, { 40, 300 } }, 6,
		{ 20, 20, 345, 435 }, { 175, 162, 46, 64 } }, PI, 177, 130, 71, 0, { 0, 0, 0 }, { 0.8, 0.045 } },
	{ { "plain_curly", "curly, voluminous",
		{ { 80, 150 }, { 60, 250 }, { 300, 250 }, { 310, 330 }, { 170, 40 }, { 240, 60 }, { 90, 380 }, { 280, 400 } }, 8,
		{ { 190, 200 }, { 190, 250 }, { 160, 160 }, { 220, 160 }, { 210, 460 }, { 190, 330 } }, 6,
		{ 20, 10, 358, 430 }, { 190, 200, 56, 78 } }, PI, 190, 172, 94, 1, { 0, 0, 0 }, { 0.24, 0.06 } },
	{ { "plain_long", "long, straight",
		{ { 150, 40 }, { 220, 40 }, { 100, 160 }, { 75, 240 }, { 250, 200 }, { 55, 300 }, { 110, 120 }, { 240, 110 } }, 8,
		{ { 180, 150 }, { 180, 200 }, { 250, 235 }, { 200, 280 }, { 300, 280 }, { 170, 300 }, { 140, 230 }, { 330, 300 } }, 8,
		{ 20, 0, 300, 341 }, { 180, 140, 42, 58 } }, PI, 180, 107, 69, 0, { 0, 0, 0 }, { 0, 0 } },
	{ { "plain_ponytail", "high ponytail, held up",
		{ { 50, 135 }, { 80, 137 }, { 110, 128 }, { 140, 112 }, { 170, 86 }, { 200, 52 }, { 222, 25 }, { 213, 70 } }, 8,
		{ { 255, 100 }, { 245, 132 }, { 290, 180 }, { 25, 150 }, { 270, 70 } }, 5,
		{ 15, 0, 300, 152 }, { 255, 100, 32, 45 } }, PI + 0.35, 245, 82, 48, 0, { 212, 22, -0.25 }, { 0, 0 } },
	{ { "plain_bob_profile", "bob, profile view",
		{ { 150, 60 }, { 250, 80 }, { 300, 200 }, { 280, 300 }, { 100, 120 }, { 320, 330 }, { 200, 120 } }, 7,
		{ { 130, 250 }, { 100, 330 }, { 340, 340 }, { 250, 450 }, { 215, 230 } }, 5,
		{ 40, 15, 360, 380 }, { 125, 260, 70, 90 } }, PI / 2, 205, 190, 118, 0, { 0, 0, 0 }, { 0, 0 } }
};

/* The hair volume (cells between the back and front depth surfaces over the
   hair mask, or behind the head for frontal photos, outside the head), the
   orientation tensors of the fitted strands as constraints, and a
   Gauss-Seidel diffusion of the tensors through the rest of the volume. */
static void hv_volume_build(hv_scene const* sc, hv_strand const* hs, int n, hv_volume* v)
{
	size_t N = (size_t)HV_GN * HV_GN * HV_GZ, id;
	int W = sc->img.w, H = sc->img.h, i, j, k, s, it, c;
	float* acc = (float*)calloc(6 * N, sizeof(float));
	float* wacc = (float*)calloc(N, sizeof(float));
	v->occ = (float*)calloc(N, sizeof(float));
	v->T = (float*)calloc(6 * N, sizeof(float));
	v->fixed = (unsigned char*)calloc(N, 1);
	v->dens = (float*)calloc(N, sizeof(float));
	for (k = 0; k < HV_GZ; k++)
		for (j = 0; j < HV_GN; j++)
			for (i = 0; i < HV_GN; i++)
			{
				double p[3], px, py, t, m, f, b;
				p[0] = g_hv_lo[0] + (i + 0.5) * HV_GH;
				p[1] = g_hv_lo[1] + (j + 0.5) * HV_GH;
				p[2] = g_hv_lo[2] + (k + 0.5) * HV_GH;
				hv_project(&sc->cam, p, &px, &py, &t);
				if (px < 0 || py < 0 || px > W - 1 || py > H - 1)
					continue;
				m = sample(sc->mask, W, H, px, py, NULL, NULL);
				if (sc->back_mask && t < 0)
					m = fmax(m, sample(sc->back_mask, W, H, px, py, NULL, NULL));
				if (m < 0.5)
					continue;
				f = sample(sc->surf[0], W, H, px, py, NULL, NULL);
				b = sample(sc->surf[1], W, H, px, py, NULL, NULL);
				if (t > f + 0.04 || t < b - 0.04 || head_phi(p, NULL) < 1.0)
					continue;
				v->occ[HV_IDX(i, j, k)] = 1;
			}
	/* constraints: tangents of the fitted strands */
	for (s = 0; s < n; s++)
	{
		qaws_curve* cv = hv_curve(&hs[s], hs[s].cps);
		qaws_range r = qaws_curve_get_parameter_range(cv);
		int m = hs[s].ncp * 10, q;
		for (q = 0; q <= m; q++)
		{
			qaws_eval_result_3d e;
			double p[3], t[3], tl;
			qaws_curve_evaluate_3d(cv, (qaws_scalar)(r.min_value + (r.max_value - r.min_value) * q / m),
				QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &e);
			p[0] = e.position.x; p[1] = e.position.y; p[2] = e.position.z;
			tl = sqrt(e.d1.x * e.d1.x + e.d1.y * e.d1.y + e.d1.z * e.d1.z) + 1e-12;
			t[0] = e.d1.x / tl; t[1] = e.d1.y / tl; t[2] = e.d1.z / tl;
			hv_splat(v, acc, wacc, p, t);
		}
		qaws_curve_destroy(cv);
	}
	for (id = 0; id < N; id++)
	{
		if (!v->occ[id])
			continue;
		if (wacc[id] > 0)
		{
			for (c = 0; c < 6; c++)
				v->T[6 * id + c] = acc[6 * id + c] / wacc[id];
			v->fixed[id] = 1;
		}
		else
			v->T[6 * id + 5] = 1;   /* start hanging down */
	}
	for (it = 0; it < 150; it++)
		for (k = 1; k < HV_GZ - 1; k++)
			for (j = 1; j < HV_GN - 1; j++)
				for (i = 1; i < HV_GN - 1; i++)
				{
					size_t nb[6];
					double sum[6] = { 0, 0, 0, 0, 0, 0 };
					int cnt = 0, q;
					id = HV_IDX(i, j, k);
					if (!v->occ[id] || v->fixed[id])
						continue;
					nb[0] = id - 1; nb[1] = id + 1; nb[2] = id - HV_GN; nb[3] = id + HV_GN;
					nb[4] = id - (size_t)HV_GN * HV_GN; nb[5] = id + (size_t)HV_GN * HV_GN;
					for (q = 0; q < 6; q++)
						if (v->occ[nb[q]])
						{
							for (c = 0; c < 6; c++)
								sum[c] += v->T[6 * nb[q] + c];
							cnt++;
						}
					if (cnt)
						for (c = 0; c < 6; c++)
							v->T[6 * id + c] = (float)(sum[c] / cnt);
				}
	free(acc);
	free(wacc);
}

static void hv_volume_free(hv_volume* v)
{
	free(v->occ); free(v->T); free(v->fixed); free(v->dens);
}

/* Traces one strand through the field from pos along dir (RK2, half-cell
   steps): it stays out of the head, gathers toward the tie of a ponytail
   until it reaches it (then leaves it away from the head), and stops when it
   leaves the volume, turns sharply or fills max points (kept every 2 steps). */
static int hv_trace(hv_volume const* v, double const* tie, double* pos, double* dir, float* pts, int max)
{
	double h = 0.5 * HV_GH, ln;
	int step, outside = 0, passed = tie ? 0 : 1, free_turn = 0, n = 0, a;
	for (step = 0; step < 2 * max - 1; step++)
	{
		double d1[3], d2[3], mid[3], gp[3], r;
		int pass;
		for (pass = 0; pass < 2; pass++)
		{
			double const* at = pass ? mid : pos;
			double* dd = pass ? d2 : d1;
			hv_field_dir(v, at, pass ? d1 : dir, dd);
			if (!passed)
			{
				/* gather toward the tie */
				double tv[3] = { tie[0] - at[0], tie[1] - at[1], tie[2] - at[2] };
				double tl = sqrt(tv[0] * tv[0] + tv[1] * tv[1] + tv[2] * tv[2]) + 1e-12;
				for (a = 0; a < 3; a++)
					dd[a] += 1.5 * tv[a] / tl;
				ln = sqrt(dd[0] * dd[0] + dd[1] * dd[1] + dd[2] * dd[2]) + 1e-12;
				for (a = 0; a < 3; a++)
					dd[a] /= ln;
			}
			if (!pass)
				for (a = 0; a < 3; a++)
					mid[a] = pos[a] + 0.5 * h * d1[a];
		}
		if (free_turn > 0)
			free_turn--;
		else if (d2[0] * dir[0] + d2[1] * dir[1] + d2[2] * dir[2] < 0.6)
			break;
		for (a = 0; a < 3; a++)
			pos[a] += h * d2[a];
		r = head_phi(pos, gp);
		if (r < 1.02)
		{
			ln = sqrt(gp[0] * gp[0] + gp[1] * gp[1] + gp[2] * gp[2]) + 1e-12;
			for (a = 0; a < 3; a++)
				pos[a] += (1.02 - r) * gp[a] / ln * 1.1;
		}
		for (a = 0; a < 3; a++)
			dir[a] = d2[a];
		if (!passed)
		{
			double dx = pos[0] - tie[0], dy = pos[1] - tie[1], dz = pos[2] - tie[2];
			if (dx * dx + dy * dy + dz * dz < 0.15 * 0.15)
			{
				/* through the tie: follow the tail, away from the head */
				double e[3];
				passed = 1;
				free_turn = 12;
				hv_field_dir(v, pos, pos, e);
				for (a = 0; a < 3; a++)
					dir[a] = e[a];
			}
		}
		if (hv_trilinear(v->occ, 1, 0, pos) < 0.5)
		{
			if (++outside > 4)
				break;
		}
		else
			outside = 0;
		if (step % 2 == 0 && n < max)
		{
			for (a = 0; a < 3; a++)
				pts[3 * n + a] = (float)pos[a];
			n++;
		}
	}
	/* trim the samples that left the volume */
	while (n > 0)
	{
		double q[3] = { pts[3 * (n - 1)], pts[3 * (n - 1) + 1], pts[3 * (n - 1) + 2] };
		if (hv_trilinear(v->occ, 1, 0, q) >= 0.5)
			break;
		n--;
	}
	return n;
}

/* Vertex colors (the blurred photo where the vertex is on the visible
   front, the strand's mean elsewhere) and the density of the shadow grid. */
static void hv_strand_finish(hv_scene const* sc, float const* rgb_blur, hv_volume* v, hv_grown* s)
{
	int W = sc->img.w, H = sc->img.h, k, a;
	double col[3] = { 0, 0, 0 }, cw = 0;
	for (k = 0; k < s->n; k++)
	{
		double q[3] = { s->p[3 * k], s->p[3 * k + 1], s->p[3 * k + 2] }, px, py, t;
		int ix, iy, vis = 0;
		hv_project(&sc->cam, q, &px, &py, &t);
		ix = (int)px;
		iy = (int)py;
		if (ix >= 0 && iy >= 0 && ix < W && iy < H && sc->mask[iy * W + ix] > 0.9f &&
			t > sample(sc->surf[0], W, H, px, py, NULL, NULL) - 0.25)
			vis = 1;
		for (a = 0; a < 3; a++)
			s->rgb[3 * k + a] = vis ? rgb_blur[3 * (iy * W + ix) + a] : -1.0f;
		if (vis)
		{
			for (a = 0; a < 3; a++)
				col[a] += s->rgb[3 * k + a];
			cw++;
		}
	}
	for (a = 0; a < 3; a++)
		col[a] = cw > 0 ? col[a] / cw : 0.85 * sc->hair_col[a];
	for (k = 0; k < s->n; k++)
	{
		int ci, cj, ck;
		double q[3] = { s->p[3 * k], s->p[3 * k + 1], s->p[3 * k + 2] };
		if (s->rgb[3 * k] < 0)
			for (a = 0; a < 3; a++)
				s->rgb[3 * k + a] = (float)(0.9 * col[a]);
		if (hv_cell(q, &ci, &cj, &ck))
			v->dens[HV_IDX(ci, cj, ck)] += 1;
	}
}

/* Curl and fuzz (HairNet, Choe and Ko): the strand is resampled finely
   enough for its curl, then offset on a parallel-transported frame by a
   helix of the given period and radius, ramped in from the root, with a
   random phase and a +-15% period jitter; straight hair only gets a slow
   random wobble of 0.012. */
static void hv_curl(hv_grown* s, double const* curl, unsigned int* rng)
{
	static float q[3 * HV_GROWN_PTS];
	double len = 0, ds, period, radius, phase, nrm[3], acc = 0, wob[2];
	int k, a, m, j = 0;
	for (k = 1; k < s->n; k++)
	{
		double dx = s->p[3 * k] - s->p[3 * k - 3], dy = s->p[3 * k + 1] - s->p[3 * k - 2], dz = s->p[3 * k + 2] - s->p[3 * k - 1];
		len += sqrt(dx * dx + dy * dy + dz * dz);
	}
	*rng = *rng * 1664525u + 1013904223u;
	phase = 2 * PI * ((*rng >> 8) / 16777216.0);
	*rng = *rng * 1664525u + 1013904223u;
	period = curl[0] > 0 ? curl[0] * (0.85 + 0.3 * ((*rng >> 8) / 16777216.0)) : 0.6 + 0.6 * ((*rng >> 8) / 16777216.0);
	radius = curl[0] > 0 ? curl[1] : 0.012;
	*rng = *rng * 1664525u + 1013904223u;
	wob[0] = (*rng >> 8) / 16777216.0;
	wob[1] = 1 - wob[0];
	ds = fmin(0.05, period / 7);
	m = (int)(len / ds) + 1;
	if (m > HV_GROWN_PTS)
		m = HV_GROWN_PTS;
	if (m < 2 || len < 1e-6)
		return;
	/* uniform resampling */
	for (k = 0; k < m; k++)
	{
		double target = k * ds, seg;
		while (j + 2 < s->n)
		{
			double dx = s->p[3 * j + 3] - s->p[3 * j], dy = s->p[3 * j + 4] - s->p[3 * j + 1], dz = s->p[3 * j + 5] - s->p[3 * j + 2];
			seg = sqrt(dx * dx + dy * dy + dz * dz);
			if (acc + seg >= target)
				break;
			acc += seg;
			j++;
		}
		{
			double dx = s->p[3 * j + 3] - s->p[3 * j], dy = s->p[3 * j + 4] - s->p[3 * j + 1], dz = s->p[3 * j + 5] - s->p[3 * j + 2];
			double f;
			seg = sqrt(dx * dx + dy * dy + dz * dz) + 1e-12;
			f = (target - acc) / seg;
			if (f > 1) f = 1;
			if (f < 0) f = 0;
			for (a = 0; a < 3; a++)
				q[3 * k + a] = (float)(s->p[3 * j + a] + f * (s->p[3 * j + 3 + a] - s->p[3 * j + a]));
		}
	}
	/* parallel-transported frame and helix offsets */
	nrm[0] = 1; nrm[1] = 0; nrm[2] = 0;
	for (k = 0; k < m; k++)
	{
		double t[3], b[3], tl, dn, ramp = fmin(1.0, k * ds / 0.3), ang = 2 * PI * k * ds / period + phase, r, c1, c2;
		int k0 = k > 0 ? k - 1 : 0, k1 = k + 1 < m ? k + 1 : m - 1;
		for (a = 0; a < 3; a++)
			t[a] = q[3 * k1 + a] - q[3 * k0 + a];
		tl = sqrt(t[0] * t[0] + t[1] * t[1] + t[2] * t[2]) + 1e-12;
		for (a = 0; a < 3; a++)
			t[a] /= tl;
		dn = nrm[0] * t[0] + nrm[1] * t[1] + nrm[2] * t[2];
		for (a = 0; a < 3; a++)
			nrm[a] -= dn * t[a];
		tl = sqrt(nrm[0] * nrm[0] + nrm[1] * nrm[1] + nrm[2] * nrm[2]);
		if (tl < 1e-6)
		{
			nrm[0] = t[1]; nrm[1] = -t[0]; nrm[2] = 0;
			tl = sqrt(nrm[0] * nrm[0] + nrm[1] * nrm[1]) + 1e-12;
		}
		for (a = 0; a < 3; a++)
			nrm[a] /= tl;
		b[0] = t[1] * nrm[2] - t[2] * nrm[1];
		b[1] = t[2] * nrm[0] - t[0] * nrm[2];
		b[2] = t[0] * nrm[1] - t[1] * nrm[0];
		r = radius * ramp;
		c1 = curl[0] > 0 ? cos(ang) : wob[0] * sin(ang);
		c2 = curl[0] > 0 ? sin(ang) : wob[1] * sin(0.7 * ang + 1.3);
		for (a = 0; a < 3; a++)
			s->p[3 * k + a] = (float)(q[3 * k + a] + r * (c1 * nrm[a] + c2 * b[a]));
	}
	s->n = m;
}

/* Strands grown from scalp roots (signed by the combing guide, or toward the
   tie of a ponytail), then fill strands (Chai et al. 2012) seeded in the
   empty cells of the volume and traced both ways, the higher end as root. */
static int hv_grow(hv_scene const* sc, double const* tie, double const* curl, hv_volume* v, hv_grown* out, int scalp_cap, int cap,
	int* scalp_count)
{
	int W = sc->img.w, H = sc->img.h, cnt = 0, i, a, cand = 4 * scalp_cap;
	unsigned int crng = 4242u;
	size_t N = (size_t)HV_GN * HV_GN * HV_GZ, id;
	float* rgb_blur = (float*)malloc(sizeof(float) * 3 * (size_t)W * H);
	float* ch = (float*)malloc(sizeof(float) * (size_t)W * H);
	for (a = 0; a < 3; a++)
	{
		for (i = 0; i < W * H; i++)
			ch[i] = sc->img.rgb[3 * i + a];
		blur(ch, W, H, 1.0);
		for (i = 0; i < W * H; i++)
			rgb_blur[3 * i + a] = ch[i];
	}
	free(ch);
	for (i = 0; i < cand && cnt < scalp_cap; i++)
	{
		double z = 1 - 2 * (i + 0.5) / cand, rad = sqrt(1 - z * z), phi = i * 2.399963229728653;
		double d[3] = { rad * cos(phi), rad * sin(phi), z }, pos[3], nrm[3], g[3], dir[3], ln, gn;
		hv_grown* s = &out[cnt];
		if (z < -0.2 || (d[1] > 0.3 && z < 0.62))
			continue;
		for (a = 0; a < 3; a++)
			pos[a] = d[a] * g_head[a] * 1.03;
		if (!hv_scalp(sc, pos) || hv_trilinear(v->occ, 1, 0, pos) < 0.25)
			continue;
		head_phi(pos, nrm);
		ln = sqrt(nrm[0] * nrm[0] + nrm[1] * nrm[1] + nrm[2] * nrm[2]);
		for (a = 0; a < 3; a++)
			nrm[a] /= ln;
		if (tie)
			for (a = 0; a < 3; a++)
				g[a] = tie[a] - pos[a];
		else
			groom_guide(pos, g);
		gn = g[0] * nrm[0] + g[1] * nrm[1] + g[2] * nrm[2];
		for (a = 0; a < 3; a++)
			g[a] -= gn * nrm[a];
		ln = sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]) + 1e-12;
		for (a = 0; a < 3; a++)
			g[a] /= ln;
		hv_field_dir(v, pos, g, dir);
		gn = dir[0] * nrm[0] + dir[1] * nrm[1] + dir[2] * nrm[2];
		for (a = 0; a < 3; a++)
			dir[a] = dir[a] - gn * nrm[a] + 0.3 * nrm[a] + 0.5 * g[a];
		ln = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]) + 1e-12;
		for (a = 0; a < 3; a++)
			dir[a] /= ln;
		s->n = hv_trace(v, tie, pos, dir, s->p, HV_GROWN_PTS);
		if (s->n < 4)
			continue;
		hv_curl(s, curl, &crng);
		hv_strand_finish(sc, rgb_blur, v, s);
		cnt++;
	}
	*scalp_count = cnt;
	{
		/* fill: empty cells in a random order */
		unsigned int rng = 99u;
		size_t* cells = (size_t*)malloc(sizeof(size_t) * N);
		size_t nc = 0, q;
		static float back[3 * HV_GROWN_PTS];
		for (id = 0; id < N; id++)
			if (v->occ[id])
				cells[nc++] = id;
		for (q = nc; q > 1; q--)
		{
			size_t j, tmp;
			rng = rng * 1664525u + 1013904223u;
			j = (size_t)(rng % (unsigned int)q);
			tmp = cells[q - 1]; cells[q - 1] = cells[j]; cells[j] = tmp;
		}
		for (q = 0; q < nc && cnt < cap; q++)
		{
			double pos[3], dir[3], down[3] = { 0, 0, -1 }, rev[3];
			int ci, cj, ck, nf, nb, k;
			hv_grown* s = &out[cnt];
			id = cells[q];
			if (v->dens[id] > 0)
				continue;
			ci = (int)(id % HV_GN);
			cj = (int)((id / HV_GN) % HV_GN);
			ck = (int)(id / ((size_t)HV_GN * HV_GN));
			pos[0] = g_hv_lo[0] + (ci + 0.5) * HV_GH;
			pos[1] = g_hv_lo[1] + (cj + 0.5) * HV_GH;
			pos[2] = g_hv_lo[2] + (ck + 0.5) * HV_GH;
			if (hv_field_dir(v, pos, down, dir) < 0.2)
				continue;
			{
				double p2[3] = { pos[0], pos[1], pos[2] };
				for (a = 0; a < 3; a++)
					rev[a] = -dir[a];
				nb = hv_trace(v, NULL, p2, rev, back, HV_GROWN_PTS / 2);
			}
			nf = hv_trace(v, NULL, pos, dir, s->p + 3 * nb, HV_GROWN_PTS - nb);
			/* backward part reversed in front */
			for (k = 0; k < nb; k++)
				for (a = 0; a < 3; a++)
					s->p[3 * k + a] = back[3 * (nb - 1 - k) + a];
			s->n = nb + nf;
			if (s->n < 8)
				continue;
			if (s->p[2] < s->p[3 * (s->n - 1) + 2])
				for (k = 0; k < s->n / 2; k++)
					for (a = 0; a < 3; a++)
					{
						float tmp = s->p[3 * k + a];
						s->p[3 * k + a] = s->p[3 * (s->n - 1 - k) + a];
						s->p[3 * (s->n - 1 - k) + a] = tmp;
					}
			hv_curl(s, curl, &crng);
			hv_strand_finish(sc, rgb_blur, v, s);
			cnt++;
		}
		free(cells);
	}
	free(rgb_blur);
	return cnt;
}

static int hv_scalp_cb(void const* ctx, double const* x)
{
	return hv_scalp((hv_scene const*)ctx, x);
}

/* One view of the grown strands as a PNG: ray-cast head, depth-sorted
   anti-aliased segments with Kajiya-Kay shading, self-shadowed by the strand
   density along the light (deep-opacity style), roots slightly darker. */
static void hv_render_png(char const* path, int w, int h, view3 const* view, hv_scene const* sc, hv_volume const* v,
	hv_grown const* gs, int n, double const* bg)
{
	enum { SS = 2 };
	hv_canvas cv;
	view3 vs = *view;
	double light[3], eye[3], ln, hh[3], skin[3] = { 0.91, 0.77, 0.67 };
	double cyw = cos(vs.yaw), syw = sin(vs.yaw), cp = cos(vs.pitch), sp = sin(vs.pitch);
	int i, k, a, ns = 0, cap = 0;
	hv_seg* segs;
	vs.cx *= SS; vs.cy *= SS; vs.scale *= SS;
	eye[0] = -syw * cp; eye[1] = -cyw * cp; eye[2] = sp;
	light[0] = 0.7 * eye[0] + 0.35 * cyw;
	light[1] = 0.7 * eye[1] - 0.35 * syw;
	light[2] = 0.7 * eye[2] + 0.55;
	ln = sqrt(light[0] * light[0] + light[1] * light[1] + light[2] * light[2]);
	for (a = 0; a < 3; a++) light[a] /= ln;
	for (a = 0; a < 3; a++) hh[a] = light[a] + eye[a];
	ln = sqrt(hh[0] * hh[0] + hh[1] * hh[1] + hh[2] * hh[2]);
	for (a = 0; a < 3; a++) hh[a] /= ln;
	hv_canvas_init(&cv, w * SS, h * SS, bg);
	hv_canvas_head(&cv, &vs, light, skin, hv_scalp_cb, sc, sc->hair_col);
	for (i = 0; i < n; i++)
		cap += gs[i].n;
	segs = (hv_seg*)malloc(sizeof(hv_seg) * (size_t)cap);
	for (i = 0; i < n; i++)
	{
		hv_grown const* s = &gs[i];
		for (k = 0; k + 1 < s->n; k++)
		{
			double p0[3], p1[3], t[3], tl, tdl, tdh, dif, spec, sh = 0, m[3], tr, x0, y0, d0, x1, y1, d1, fade;
			int q;
			hv_seg* sg = &segs[ns];
			for (a = 0; a < 3; a++)
			{
				p0[a] = s->p[3 * k + a];
				p1[a] = s->p[3 * (k + 1) + a];
				t[a] = p1[a] - p0[a];
				m[a] = 0.5 * (p0[a] + p1[a]);
			}
			tl = sqrt(t[0] * t[0] + t[1] * t[1] + t[2] * t[2]) + 1e-12;
			for (a = 0; a < 3; a++) t[a] /= tl;
			tdl = t[0] * light[0] + t[1] * light[1] + t[2] * light[2];
			tdh = t[0] * hh[0] + t[1] * hh[1] + t[2] * hh[2];
			dif = sqrt(fabs(1 - tdl * tdl));
			spec = pow(sqrt(fabs(1 - tdh * tdh)), 60.0);
			for (q = 1; q <= 24; q++)
			{
				double sp3[3];
				int ci, cj, ck;
				for (a = 0; a < 3; a++)
					sp3[a] = m[a] + light[a] * HV_GH * q;
				if (!hv_cell(sp3, &ci, &cj, &ck))
					break;
				sh += v->dens[HV_IDX(ci, cj, ck)];
			}
			tr = 0.3 + 0.7 * exp(-0.012 * sh);
			fade = 0.8 + 0.2 * k / (double)(s->n - 1);
			for (a = 0; a < 3; a++)
			{
				double base = s->rgb[3 * k + a];
				sg->rgb[a] = (float)fmin(1.0, (base * (0.45 + 0.55 * dif) * fade + 0.22 * spec) * tr);
			}
			view_xform(&vs, p0, &x0, &y0, &d0);
			view_xform(&vs, p1, &x1, &y1, &d1);
			sg->x0 = (float)x0; sg->y0 = (float)y0; sg->d0 = (float)d0;
			sg->x1 = (float)x1; sg->y1 = (float)y1; sg->d1 = (float)d1;
			sg->depth = (float)(0.5 * (d0 + d1));
			ns++;
		}
	}
	hv_canvas_segments(&cv, segs, ns, 1.1, 0.8);
	free(segs);
	hv_canvas_save(&cv, SS, path);
}

/* Mean angle (degrees) between the projected visible strand tangents and
   the photo's orientation field, and the share of hair pixels they cover. */
static void hv_grown_stats(hv_scene const* sc, hv_grown const* gs, int n, double* angle, double* coverage)
{
	int W = sc->img.w, H = sc->img.h, i, k, cnt = 0, hair = 0, cov = 0;
	unsigned char* hit = (unsigned char*)calloc((size_t)W * H, 1);
	double sum = 0;
	for (i = 0; i < n; i++)
		for (k = 0; k + 1 < gs[i].n; k++)
		{
			double p0[3] = { gs[i].p[3 * k], gs[i].p[3 * k + 1], gs[i].p[3 * k + 2] };
			double p1[3] = { gs[i].p[3 * k + 3], gs[i].p[3 * k + 4], gs[i].p[3 * k + 5] };
			double x0, y0, t0, x1, y1, t1, dx, dy, dl, ox, oy, c, dot;
			hv_project(&sc->cam, p0, &x0, &y0, &t0);
			hv_project(&sc->cam, p1, &x1, &y1, &t1);
			if (x0 < 2 || y0 < 2 || x0 > W - 3 || y0 > H - 3 || sc->mask[(int)y0 * W + (int)x0] < 0.5f ||
				t0 < sample(sc->surf[0], W, H, x0, y0, NULL, NULL) - 0.25)
				continue;
			hit[(int)y0 * W + (int)x0] = 1;
			dx = x1 - x0;
			dy = y1 - y0;
			dl = sqrt(dx * dx + dy * dy);
			if (dl < 0.3)
				continue;
			c = tensor_direction(&sc->tf, x0, y0, &ox, &oy);
			if (c < 0.3)
				continue;
			dot = fabs(dx * ox + dy * oy) / dl;
			sum += acos(dot > 1 ? 1 : dot) * 180.0 / PI;
			cnt++;
		}
	for (i = 0; i < W * H; i++)
		if (sc->mask[i] > 0.5f)
		{
			int x = i % W, y = i / W, dx, dy, any = 0;
			hair++;
			for (dy = -1; dy <= 1 && !any; dy++)
				for (dx = -1; dx <= 1 && !any; dx++)
					if (x + dx >= 0 && y + dy >= 0 && x + dx < W && y + dy < H && hit[(y + dy) * W + x + dx])
						any = 1;
			cov += any;
		}
	*angle = cnt ? sum / cnt : 0;
	*coverage = hair ? (double)cov / hair : 0;
	free(hit);
}

static void app_hair_volume(void)
{
	static strand traced[HS_MAX_STRANDS];
	static hv_strand hs[HV_MAX];
	static hv_grown grown[HV_GROWN_MAX];
	static hv_scene sc;
	static hv_volume vol;
	double const bg_cam[3] = { 1, 1, 1 }, bg_new[3] = { 0.953, 0.945, 0.933 };
	int c;
	double y0 = 80;
	svg s;
	char buf[256], path[512], png[512];

	MAKE_DIR("showcase/app10");
	svg_open(&s, "showcase/app10_hair_volume.svg", 1230, 1560, "Single-view hair modeling",
		"Plain-background portraits. Gabor orientation, image strands lifted onto a hair volume and fitted as 3D B-splines "
		"(projection adjoints); their tangents constrain a diffused 3D orientation field; strands grown from the scalp.");
	for (c = 0; c < 5; c++)
	{
		hv_case const* hc = &g_hv_cases[c];
		hair_spec const* sp = &hc->spec;
		int nt = 0, n, ng, nscalp, v;
		double rms, k = 240.0 / 360, ih, angle, coverage, tie[3];
		sprintf(path, "%s/%s.ppm", g_photos, sp->name);
		if (!image_load_ppm(path, &sc.img))
		{
			printf("10_hair_volume: %s missing\n", path);
			continue;
		}
		n = hv_reconstruct(&sc, hc, traced, &nt, hs, &rms);
		hv_volume_build(&sc, hs, n, &vol);
		if (hc->tie[0] > 0)
			hv_unproject(&sc.cam, hc->tie[0], hc->tie[1], hc->tie[2], tie);
		ng = hv_grow(&sc, hc->tie[0] > 0 ? tie : NULL, hc->curl, &vol, grown, 14000, HV_GROWN_MAX, &nscalp);
		hv_grown_stats(&sc, grown, ng, &angle, &coverage);
		ih = sc.img.h * k;
		if (ih > 315) ih = 315;
		if (ih < 200) ih = 200;
		/* 1. photo */
		fprintf(s.f, "<clipPath id=\"hvp%d\"><rect x=\"20\" y=\"%.1f\" width=\"240\" height=\"%.1f\"/></clipPath>\n", c, y0, ih);
		fprintf(s.f, "<image href=\"../%s/%s.png\" x=\"20\" y=\"%.1f\" width=\"240\" height=\"%.1f\" clip-path=\"url(#hvp%d)\"/>\n",
			g_photos, sp->name, y0, sc.img.h * k, c);
		/* 2. the grown strands from the photo camera; 3. new views */
		for (v = 0; v < 4; v++)
		{
			view3 vw;
			int pw = v == 0 ? 240 : 220;
			double x0 = v == 0 ? 270 : 520 + (v - 1) * 230;
			if (v == 0)
			{
				vw.cx = sc.cam.cx * k; vw.cy = sc.cam.cy * k; vw.scale = sc.cam.s * k; vw.yaw = sc.cam.yaw; vw.pitch = 0;
			}
			else
			{
				double yaws[3] = { hc->yaw - 0.9, hc->yaw + PI / 2, hc->yaw + PI };
				vw.cx = 110; vw.cy = ih * 0.36; vw.scale = ih / 6.2; vw.yaw = yaws[v - 1]; vw.pitch = 0.12;
			}
			sprintf(png, "showcase/app10/%s_%d.png", sp->name, v);
			hv_render_png(png, pw, (int)ih, &vw, &sc, &vol, grown, ng, v == 0 ? bg_cam : bg_new);
			fprintf(s.f, "<image href=\"app10/%s_%d.png\" x=\"%.1f\" y=\"%.1f\" width=\"%d\" height=\"%d\"/>\n",
				sp->name, v, x0, y0, pw, (int)ih);
		}
		sprintf(buf, "%s: %d image strands -> %d fitted -> %d strands (%d from the scalp); angle to the photo %.1f deg, hair covered %.0f%%",
			sp->label, nt, n, ng, nscalp, angle, 100 * coverage);
		svg_text(&s, 20, y0 + ih + 16, 12, "#0969da", "start", buf);
		printf("10_hair_volume: %-18s %d image strands, %d fitted (RMS %.2f px), %d strands (%d scalp), angle %.1f deg, coverage %.0f%%\n",
			sp->name, nt, n, rms, ng, nscalp, angle, 100 * coverage);
		hv_volume_free(&vol);
		hv_free(&sc, traced, nt);
		y0 += ih + 32;
	}
	svg_text(&s, 20, y0 + 8, 11, "#57606a", "start",
		"Columns: photo; grown strands from the photo camera; three new views. Photos: Wikimedia Commons (see photos/CREDITS.txt).");
	svg_close(&s);
}

