/*
 * Applications of the differentiation API on real data.
 *
 *   1. Hair strands from a photo: orientation field (structure tensor),
 *      evenly spaced streamlines, B-spline strands optimized to follow the
 *      field through unit-tangent adjoints (qaws_diff_geometry.h).
 *
 * Photos are read as binary PPM; examples/photo_to_ppm.ps1 converts any
 * image (and writes a PNG used as figure background). Figures are written
 * to showcase/ in the working directory.
 *
 *   qaws_diff_applications [photos directory]   (default: photos)
 */

#include "example_svg.h"

#ifdef _WIN32
#include <direct.h>
#define MAKE_DIR(p) _mkdir(p)
#else
#include <sys/stat.h>
#define MAKE_DIR(p) mkdir(p, 0755)
#endif

static char const* g_photos = "photos";

/* ================================================================== */
/*  Images                                                            */
/* ================================================================== */

typedef struct image
{
	int w, h;
	float* rgb;   /* w * h * 3, 0..1 */
	float* lum;   /* w * h */
} image;

static int read_token(FILE* f, char* buf, int cap)
{
	int c, n = 0;
	do
	{
		c = fgetc(f);
		if (c == '#')
			while (c != '\n' && c != EOF)
				c = fgetc(f);
	} while (c == ' ' || c == '\n' || c == '\r' || c == '\t');
	while (c != EOF && c != ' ' && c != '\n' && c != '\r' && c != '\t' && n < cap - 1)
	{
		buf[n++] = (char)c;
		c = fgetc(f);
	}
	buf[n] = 0;
	return n;
}

static int image_load_ppm(char const* path, image* img)
{
	FILE* f = fopen(path, "rb");
	char tok[32];
	int maxv, i, n;
	unsigned char* raw;
	memset(img, 0, sizeof(*img));
	if (!f)
		return 0;
	read_token(f, tok, 32);
	if (strcmp(tok, "P6") != 0)
	{
		fclose(f);
		return 0;
	}
	read_token(f, tok, 32); img->w = atoi(tok);
	read_token(f, tok, 32); img->h = atoi(tok);
	read_token(f, tok, 32); maxv = atoi(tok);
	n = img->w * img->h;
	raw = (unsigned char*)malloc((size_t)n * 3);
	img->rgb = (float*)malloc(sizeof(float) * (size_t)n * 3);
	img->lum = (float*)malloc(sizeof(float) * (size_t)n);
	if (fread(raw, 1, (size_t)n * 3, f) != (size_t)n * 3)
	{
		fclose(f);
		free(raw);
		return 0;
	}
	fclose(f);
	for (i = 0; i < n * 3; i++)
		img->rgb[i] = raw[i] / (float)maxv;
	for (i = 0; i < n; i++)
		img->lum[i] = 0.299f * img->rgb[3 * i] + 0.587f * img->rgb[3 * i + 1] + 0.114f * img->rgb[3 * i + 2];
	free(raw);
	return 1;
}

static void image_free(image* img)
{
	free(img->rgb);
	free(img->lum);
}

/* Separable Gaussian blur of a w x h field. */
static void blur(float* f, int w, int h, double sigma)
{
	int r = (int)ceil(3 * sigma), x, y, k;
	float* tmp = (float*)malloc(sizeof(float) * (size_t)w * h);
	float* ker = (float*)malloc(sizeof(float) * (size_t)(2 * r + 1));
	double sum = 0;
	for (k = -r; k <= r; k++)
		sum += ker[k + r] = (float)exp(-k * k / (2 * sigma * sigma));
	for (k = 0; k <= 2 * r; k++)
		ker[k] = (float)(ker[k] / sum);
	for (y = 0; y < h; y++)
		for (x = 0; x < w; x++)
		{
			double a = 0;
			for (k = -r; k <= r; k++)
			{
				int xx = x + k < 0 ? 0 : (x + k >= w ? w - 1 : x + k);
				a += ker[k + r] * f[y * w + xx];
			}
			tmp[y * w + x] = (float)a;
		}
	for (y = 0; y < h; y++)
		for (x = 0; x < w; x++)
		{
			double a = 0;
			for (k = -r; k <= r; k++)
			{
				int yy = y + k < 0 ? 0 : (y + k >= h ? h - 1 : y + k);
				a += ker[k + r] * tmp[yy * w + x];
			}
			f[y * w + x] = (float)a;
		}
	free(tmp);
	free(ker);
}

/* Bilinear sample with its spatial gradient (pixel units, clamped). */
static double sample(float const* f, int w, int h, double x, double y, double* dx, double* dy)
{
	int x0, y0;
	double fx, fy, a, b, c, d;
	if (x < 0) x = 0;
	if (y < 0) y = 0;
	if (x > w - 1.001) x = w - 1.001;
	if (y > h - 1.001) y = h - 1.001;
	x0 = (int)x;
	y0 = (int)y;
	fx = x - x0;
	fy = y - y0;
	a = f[y0 * w + x0];
	b = f[y0 * w + x0 + 1];
	c = f[(y0 + 1) * w + x0];
	d = f[(y0 + 1) * w + x0 + 1];
	if (dx) *dx = (1 - fy) * (b - a) + fy * (d - c);
	if (dy) *dy = (1 - fx) * (c - a) + fx * (d - b);
	return (1 - fy) * ((1 - fx) * a + fx * b) + fy * ((1 - fx) * c + fx * d);
}

/* ================================================================== */
/*  1. Hair strands from a photo                                      */
/* ================================================================== */

#define HS_MAX_STRANDS 2000
#define HS_MAX_POINTS 900
#define HS_ITERS 120
#define HS_SAMPLES_PER_CP 4

typedef struct tensor_field
{
	int w, h;
	float *xx, *xy, *yy;
} tensor_field;

static void tensor_build(image const* img, double sigma_d, double sigma_t, tensor_field* t)
{
	int n = img->w * img->h, x, y;
	float* l = (float*)malloc(sizeof(float) * (size_t)n);
	memcpy(l, img->lum, sizeof(float) * (size_t)n);
	blur(l, img->w, img->h, sigma_d);
	t->w = img->w;
	t->h = img->h;
	t->xx = (float*)calloc((size_t)n, sizeof(float));
	t->xy = (float*)calloc((size_t)n, sizeof(float));
	t->yy = (float*)calloc((size_t)n, sizeof(float));
	for (y = 1; y < img->h - 1; y++)
		for (x = 1; x < img->w - 1; x++)
		{
			float const* p = l + y * img->w + x;
			double gx = (p[1 - img->w] + 2 * p[1] + p[1 + img->w] - p[-1 - img->w] - 2 * p[-1] - p[-1 + img->w]) / 8;
			double gy = (p[img->w - 1] + 2 * p[img->w] + p[img->w + 1] - p[-img->w - 1] - 2 * p[-img->w] - p[-img->w + 1]) / 8;
			t->xx[y * img->w + x] = (float)(gx * gx);
			t->xy[y * img->w + x] = (float)(gx * gy);
			t->yy[y * img->w + x] = (float)(gy * gy);
		}
	blur(t->xx, img->w, img->h, sigma_t);
	blur(t->xy, img->w, img->h, sigma_t);
	blur(t->yy, img->w, img->h, sigma_t);
	free(l);
}

/* Unit direction along the stripes and the coherence at (x, y). */
static double tensor_direction(tensor_field const* t, double x, double y, double* ox, double* oy)
{
	double a = sample(t->xx, t->w, t->h, x, y, NULL, NULL);
	double b = sample(t->xy, t->w, t->h, x, y, NULL, NULL);
	double c = sample(t->yy, t->w, t->h, x, y, NULL, NULL);
	double theta = 0.5 * atan2(2 * b, a - c) + 0.5 * PI;   /* gradient direction + 90 degrees */
	*ox = cos(theta);
	*oy = sin(theta);
	return sqrt((a - c) * (a - c) + 4 * b * b) / (a + c + 1e-12);
}

/* Alignment energy e = T^T J T / tr(J) of a unit tangent T at x, with
   its gradients (zero when T follows the stripes). */
static double align_energy(tensor_field const* t, double x, double y, double tx, double ty,
	double* de_dx, double* de_dy, double* de_dtx, double* de_dty)
{
	double ax, ay, bx, by, cx, cy;
	double a = sample(t->xx, t->w, t->h, x, y, &ax, &ay);
	double b = sample(t->xy, t->w, t->h, x, y, &bx, &by);
	double c = sample(t->yy, t->w, t->h, x, y, &cx, &cy);
	double tr = a + c + 1e-12;
	double q = a * tx * tx + 2 * b * tx * ty + c * ty * ty;
	double e = q / tr;
	/* d(q / tr) = (dq tr - q dtr) / tr^2 */
	double dqx = ax * tx * tx + 2 * bx * tx * ty + cx * ty * ty;
	double dqy = ay * tx * tx + 2 * by * tx * ty + cy * ty * ty;
	*de_dx = (dqx - e * (ax + cx)) / tr;
	*de_dy = (dqy - e * (ay + cy)) / tr;
	*de_dtx = 2 * (a * tx + b * ty) / tr;
	*de_dty = 2 * (b * tx + c * ty) / tr;
	return e;
}

typedef struct strand
{
	int n;          /* traced points */
	qaws_scalar* xy;  /* traced polyline */
	qaws_curve* fit;
	qaws_curve* opt;
	float color[3];
	qaws_scalar knots[32];   /* knot vector of the fit, kept while optimizing */
} strand;

/* Evenly spaced streamlines (Jobard-Lefer style) along the stripes. */
static int trace_strands(image const* img, tensor_field const* t, double dsep, double min_coherence,
	double min_lum, strand* out, int cap)
{
	int gw = (int)(img->w / dsep) + 1, gh = (int)(img->h / dsep) + 1;
	unsigned char* occ = (unsigned char*)calloc((size_t)gw * gh, 1);
	float* buf = (float*)malloc(sizeof(float) * 2 * 2 * HS_MAX_POINTS);
	int count = 0, sx, sy, k;
	unsigned int rng = 12345u;
	int* order = (int*)malloc(sizeof(int) * (size_t)gw * gh);

	for (k = 0; k < gw * gh; k++)
		order[k] = k;
	for (k = gw * gh - 1; k > 0; k--)
	{
		int j, tmp;
		rng = rng * 1664525u + 1013904223u;
		j = (int)(rng % (unsigned int)(k + 1));
		tmp = order[k]; order[k] = order[j]; order[j] = tmp;
	}
	for (k = 0; k < gw * gh && count < cap; k++)
	{
		double x0, y0;
		int dir, nb = 0, nf = 0, i;
		float* back = buf;
		float* fwd = buf + 2 * HS_MAX_POINTS;
		sx = order[k] % gw;
		sy = order[k] / gw;
		if (occ[sy * gw + sx])
			continue;
		x0 = (sx + 0.5) * dsep;
		y0 = (sy + 0.5) * dsep;
		if (x0 >= img->w - 2 || y0 >= img->h - 2)
			continue;
		{
			double ox, oy;
			if (tensor_direction(t, x0, y0, &ox, &oy) < min_coherence ||
			    sample(img->lum, img->w, img->h, x0, y0, NULL, NULL) < min_lum)
				continue;
		}
		for (dir = 0; dir < 2; dir++)
		{
			double x = x0, y = y0, px = 0, py = 0;
			float* line = dir ? back : fwd;
			int* nline = dir ? &nb : &nf;
			int step;
			{
				double ox, oy;
				tensor_direction(t, x, y, &ox, &oy);
				px = dir ? -ox : ox;
				py = dir ? -oy : oy;
			}
			for (step = 0; step < HS_MAX_POINTS; step++)
			{
				double ox, oy, mx, my, c;
				int gx, gy;
				/* midpoint step, orientation kept consistent with the previous step */
				tensor_direction(t, x, y, &ox, &oy);
				if (ox * px + oy * py < 0) { ox = -ox; oy = -oy; }
				mx = x + 0.75 * ox;
				my = y + 0.75 * oy;
				c = tensor_direction(t, mx, my, &ox, &oy);
				if (ox * px + oy * py < 0) { ox = -ox; oy = -oy; }
				x += 1.5 * ox;
				y += 1.5 * oy;
				px = ox;
				py = oy;
				if (x < 2 || y < 2 || x > img->w - 3 || y > img->h - 3 || c < min_coherence ||
				    sample(img->lum, img->w, img->h, x, y, NULL, NULL) < min_lum)
					break;
				gx = (int)(x / dsep);
				gy = (int)(y / dsep);
				if (occ[gy * gw + gx] && step > 2)
					break;
				line[2 * *nline] = (float)x;
				line[2 * *nline + 1] = (float)y;
				(*nline)++;
			}
		}
		if (nb + nf + 1 < 20)
			continue;
		{
			strand* s = &out[count++];
			s->n = nb + nf + 1;
			s->xy = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * (size_t)s->n);
			for (i = 0; i < nb; i++)
			{
				s->xy[2 * i] = (qaws_scalar)back[2 * (nb - 1 - i)];
				s->xy[2 * i + 1] = back[2 * (nb - 1 - i) + 1];
			}
			s->xy[2 * nb] = (qaws_scalar)x0;
			s->xy[2 * nb + 1] = (qaws_scalar)y0;
			for (i = 0; i < nf; i++)
			{
				s->xy[2 * (nb + 1 + i)] = fwd[2 * i];
				s->xy[2 * (nb + 1 + i) + 1] = fwd[2 * i + 1];
			}
			for (i = 0; i < s->n; i++)
				occ[(int)(s->xy[2 * i + 1] / dsep) * gw + (int)(s->xy[2 * i] / dsep)] = 1;
			{
				int mid = s->n / 2, c;
				int px = (int)s->xy[2 * mid], py = (int)s->xy[2 * mid + 1];
				for (c = 0; c < 3; c++)
					s->color[c] = img->rgb[(py * img->w + px) * 3 + c];
			}
		}
	}
	free(occ);
	free(buf);
	free(order);
	return count;
}

static unsigned int strand_cp_count(strand const* s)
{
	unsigned int n = (unsigned int)(s->n / 40) + 1;
	return n < 4 ? 4 : (n > 24 ? 24 : n);
}

static qaws_curve* strand_curve(qaws_scalar const* cps, unsigned int n, qaws_scalar const* knots)
{
	qaws_bspline_desc d;
	qaws_curve* c = NULL;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D;
	d.degree = 3;
	d.control_points = cps;
	d.control_point_count = n;
	d.knots = knots;
	d.knot_count = n + 4;
	qaws_curve_create_bspline(&d, &c);
	return c;
}

/* Mean alignment energy and angle (degrees) of a strand curve. */
static void strand_alignment(tensor_field const* t, qaws_curve const* c, double* energy, double* angle)
{
	qaws_range r = qaws_curve_get_parameter_range(c);
	int k, m = 64;
	*energy = 0;
	*angle = 0;
	for (k = 0; k < m; k++)
	{
		qaws_eval_result_2d e;
		double len, tx, ty, ox, oy, a, b, cx, cy, d;
		qaws_curve_evaluate_2d(c, (qaws_scalar)(r.min_value + (r.max_value - r.min_value) * (k + 0.5) / m),
			QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &e);
		len = sqrt(e.d1.x * e.d1.x + e.d1.y * e.d1.y) + 1e-12;
		tx = e.d1.x / len;
		ty = e.d1.y / len;
		*energy += align_energy(t, e.position.x, e.position.y, tx, ty, &a, &b, &cx, &cy) / m;
		tensor_direction(t, e.position.x, e.position.y, &ox, &oy);
		d = fabs(tx * ox + ty * oy);
		*angle += acos(d > 1 ? 1 : d) * 180.0 / PI / m;
	}
}

/* Strand energy: alignment of the unit tangent with the stripes (through
   qaws_curve_geometry_adjoint_2d), a spring to the traced polyline and
   bending. Returns the energy and accumulates the control point gradient. */
static double strand_energy(tensor_field const* t, strand const* s, qaws_scalar const* cps, unsigned int ncp,
	qaws_scalar* grad)
{
	enum { MAXS = 24 * HS_SAMPLES_PER_CP + 8 };
	qaws_curve* c = strand_curve(cps, ncp, s->knots);
	qaws_range r = qaws_curve_get_parameter_range(c);
	unsigned int m = ncp * HS_SAMPLES_PER_CP + 8, k;
	qaws_scalar ts[MAXS];
	qaws_curve_jet_2d prim[MAXS], tan[MAXS], bar[MAXS];
	qaws_field_view fv;
	qaws_diff_views views = one_field(&fv, QAWS_FIELD_CONTROL_POINTS, grad, ncp, 2);
	double e = 0, spring = 0.0005, bend_w = 0.3;

	for (k = 0; k < m; k++)
		ts[k] = (qaws_scalar)(r.min_value + (r.max_value - r.min_value) * (k + 0.5) / m);
	qaws_curve_eval_batch_tangent_2d(NULL, c, ts, NULL, m, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2,
		NULL, prim, tan);
	memset(bar, 0, sizeof(qaws_curve_jet_2d) * m);
	for (k = 0; k < m; k++)
	{
		qaws_curve_geometry_2d g, gb;
		qaws_diff_validity validity;
		double dex, dey, detx, dety, ea;
		/* spring target: the traced point at the same relative position */
		double f = (k + 0.5) / m * (s->n - 1);
		int i0 = (int)f;
		double w = f - i0;
		double qx = s->xy[2 * i0] * (1 - w) + s->xy[2 * (i0 + 1 < s->n ? i0 + 1 : i0)] * w;
		double qy = s->xy[2 * i0 + 1] * (1 - w) + s->xy[2 * (i0 + 1 < s->n ? i0 + 1 : i0) + 1] * w;
		double rx = prim[k].d[0].x - qx, ry = prim[k].d[0].y - qy;

		qaws_curve_geometry_eval_2d(&prim[k], NULL, NULL, &g, NULL, NULL, &validity);
		ea = align_energy(t, prim[k].d[0].x, prim[k].d[0].y, g.tangent.x, g.tangent.y, &dex, &dey, &detx, &dety);
		e += (ea + spring * (rx * rx + ry * ry)) / m;
		bar[k].d[0].x = (qaws_scalar)((dex + 2 * spring * rx) / m);
		bar[k].d[0].y = (qaws_scalar)((dey + 2 * spring * ry) / m);
		bar[k].channels = QAWS_EVAL_FLAG_POSITION;
		memset(&gb, 0, sizeof(gb));
		gb.tangent.x = (qaws_scalar)(detx / m);
		gb.tangent.y = (qaws_scalar)(dety / m);
		qaws_curve_geometry_adjoint_2d(&prim[k], &gb, &bar[k], &validity);
	}
	qaws_curve_eval_batch_adjoint_2d(NULL, c, ts, m, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2,
		bar, &views, NULL);
	{
		/* bending per unit length cubed keeps the weight scale free */
		qaws_scalar bg[48], bend = 0;
		qaws_field_view bf;
		qaws_diff_views bv = one_field(&bf, QAWS_FIELD_CONTROL_POINTS, bg, ncp, 2);
		double len = s->n * 1.5, scale = bend_w / (len * len * len);
		memset(bg, 0, sizeof(bg));
		qaws_curve_functional_gradient(NULL, c, QAWS_FUNCTIONAL_BENDING, 0, &bv, &bend);
		e += scale * bend;
		for (k = 0; k < 2 * ncp; k++)
			grad[k] += (qaws_scalar)(scale * bg[k]);
	}
	qaws_curve_destroy(c);
	return e;
}

/* Image pixel -> screen mapping of a panel (with an optional crop). */
typedef struct imgmap
{
	double x0, y0, ox, oy, scale;
} imgmap;

static void map_pt(imgmap const* m, double x, double y, double* sx, double* sy)
{
	*sx = m->x0 + (x - m->ox) * m->scale;
	*sy = m->y0 + (y - m->oy) * m->scale;
}

static void strand_svg(svg* s, imgmap const* m, qaws_curve const* c, char const* color, double width, double opacity)
{
	double xy[2 * 160];
	qaws_range r = qaws_curve_get_parameter_range(c);
	int k;
	for (k = 0; k < 160; k++)
	{
		qaws_eval_result_2d e;
		qaws_curve_evaluate_2d(c, (qaws_scalar)(r.min_value + (r.max_value - r.min_value) * k / 159.0), QAWS_EVAL_FLAG_POSITION, &e);
		map_pt(m, e.position.x, e.position.y, &xy[2 * k], &xy[2 * k + 1]);
	}
	svg_polyline(s, xy, 160, color, width, opacity, 0);
}

/* Strand drawn in short segments taking the photo color under them. */
static void strand_svg_photo(svg* s, imgmap const* m, qaws_curve const* c, image const* img, double width)
{
	qaws_range r = qaws_curve_get_parameter_range(c);
	double px = 0, py = 0;
	int k, n = 48;
	for (k = 0; k <= n; k++)
	{
		qaws_eval_result_2d e;
		double sx, sy;
		qaws_curve_evaluate_2d(c, (qaws_scalar)(r.min_value + (r.max_value - r.min_value) * k / n), QAWS_EVAL_FLAG_POSITION, &e);
		map_pt(m, e.position.x, e.position.y, &sx, &sy);
		if (k > 0)
		{
			char col[40];
			int ix = (int)(e.position.x + 0.5), iy = (int)(e.position.y + 0.5);
			float const* p;
			ix = ix < 0 ? 0 : (ix >= img->w ? img->w - 1 : ix);
			iy = iy < 0 ? 0 : (iy >= img->h ? img->h - 1 : iy);
			p = img->rgb + 3 * (iy * img->w + ix);
			sprintf(col, "rgb(%d,%d,%d)", (int)(255 * p[0]), (int)(255 * p[1]), (int)(255 * p[2]));
			svg_line(s, px, py, sx, sy, col, width, 0.95);
		}
		px = sx;
		py = sy;
	}
}

static void app_hair(void)
{
	static strand strands[HS_MAX_STRANDS];
	char path[512], buf[256];
	image img;
	tensor_field tf;
	int ns, i, total_points = 0;
	unsigned int total_cp = 0;
	double e0 = 0, a0 = 0, e1 = 0, a1 = 0;
	svg s;

	sprintf(path, "%s/ribbons.ppm", g_photos);
	if (!image_load_ppm(path, &img))
	{
		printf("1_hair: %s not found (see examples/photo_to_ppm.ps1)\n", path);
		return;
	}
	tensor_build(&img, 1.0, 3.0, &tf);
	ns = trace_strands(&img, &tf, 2.2, 0.35, 0.10, strands, HS_MAX_STRANDS);

	for (i = 0; i < ns; i++)
	{
		strand* st = &strands[i];
		unsigned int ncp = strand_cp_count(st), k;
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
		if (qaws_curve_fit_bspline(&d, &st->fit) != QAWS_STATUS_OK)
		{
			st->fit = NULL;
			continue;
		}
		/* the fit as a curve we can rebuild from control points */
		{
			unsigned int got = 0;
			qaws_curve_read_field(st->fit, QAWS_FIELD_CONTROL_POINTS, cps, ncp * 2, &got);
			qaws_curve_read_field(st->fit, QAWS_FIELD_KNOTS, st->knots, ncp + 4, &got);
		}
		qaws_curve_destroy(st->fit);
		st->fit = strand_curve(cps, ncp, st->knots);
		strand_alignment(&tf, st->fit, &ea, &aa);
		e0 += ea; a0 += aa;
		memset(&opt, 0, sizeof(opt));
		for (it = 0; it < HS_ITERS; it++)
		{
			memset(g, 0, sizeof(g));
			strand_energy(&tf, st, cps, ncp, g);
			adam_step(&opt, cps, g, (int)(2 * ncp), 0.35 * (1.0 - 0.8 * it / (double)HS_ITERS));
		}
		st->opt = strand_curve(cps, ncp, st->knots);
		strand_alignment(&tf, st->opt, &ea, &aa);
		e1 += ea; a1 += aa;
		total_points += st->n;
		total_cp += ncp;
		(void)k;
	}

	svg_open(&s, "showcase/app1_hair_strands.svg", 1240, 830, "Hair strands from a photo",
		"Orientation field from the structure tensor, evenly spaced streamlines, B-spline strands optimized "
		"to follow it through unit-tangent adjoints.");
	{
		double pw = 590, sc = pw / img.w, ph = img.h * sc;
		imgmap ma = { 20, 100, 0, 0, sc }, mb = { 630, 100, 0, 0, sc };
		/* zoom: a crop of the photo at 3.4x */
		double zx0 = 196, zy0 = 120, zs = 3.4, zw = 150, zh = 88;
		imgmap mz = { 20, 100 + ph + 50, zx0, zy0, zs };
		double zpw = zw * zs, zph = zh * zs;

		fprintf(s.f, "<text x=\"20\" y=\"92\" font-size=\"13\" font-weight=\"600\" fill=\"#24292f\">photo</text>\n");
		fprintf(s.f, "<text x=\"630\" y=\"92\" font-size=\"13\" font-weight=\"600\" fill=\"#24292f\">%d optimized B-spline strands, colored by the photo</text>\n", ns);
		fprintf(s.f, "<image href=\"../%s/ribbons.png\" x=\"20\" y=\"100\" width=\"%.1f\" height=\"%.1f\"/>\n", g_photos, pw, ph);
		fprintf(s.f, "<rect x=\"630\" y=\"100\" width=\"%.1f\" height=\"%.1f\" fill=\"#05070a\"/>\n", pw, ph);
		for (i = 0; i < ns; i++)
			if (strands[i].opt)
				strand_svg_photo(&s, &mb, strands[i].opt, &img, 1.25);
		/* zoom frame on the photo */
		fprintf(s.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"none\" stroke=\"#ffffff\" stroke-width=\"1.5\"/>\n",
			20 + zx0 * sc, 100 + zy0 * sc, zw * sc, zh * sc);

		/* zoom panel: orientation glyphs, streamlines, initial fits, optimized strands */
		fprintf(s.f, "<text x=\"20\" y=\"%.1f\" font-size=\"13\" font-weight=\"600\" fill=\"#24292f\">zoom: orientation field, "
			"streamline (grey), initial fit (orange), optimized (cyan)</text>\n", mz.y0 - 8);
		fprintf(s.f, "<clipPath id=\"zoomclip\"><rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\"/></clipPath>\n",
			mz.x0, mz.y0, zpw, zph);
		fprintf(s.f, "<g clip-path=\"url(#zoomclip)\">\n");
		fprintf(s.f, "<image href=\"../%s/ribbons.png\" x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" opacity=\"0.35\"/>\n",
			g_photos, mz.x0 - zx0 * zs, mz.y0 - zy0 * zs, img.w * zs, img.h * zs);
		fprintf(s.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"#05070a\" opacity=\"0.55\"/>\n", mz.x0, mz.y0, zpw, zph);
		{
			double gx, gy;
			for (gy = zy0 + 2; gy < zy0 + zh; gy += 5)
				for (gx = zx0 + 2; gx < zx0 + zw; gx += 5)
				{
					double ox, oy, sx, sy, co = tensor_direction(&tf, gx, gy, &ox, &oy);
					char col[32];
					heat(co, col);
					map_pt(&mz, gx, gy, &sx, &sy);
					svg_line(&s, sx - 6 * ox, sy - 6 * oy, sx + 6 * ox, sy + 6 * oy, col, 1.1, 0.55);
				}
		}
		for (i = 0; i < ns; i++)
		{
			strand* st = &strands[i];
			double xy[2 * HS_MAX_POINTS * 2];
			int k, m = 0;
			if (!st->opt)
				continue;
			for (k = 0; k < st->n; k++, m++)
				map_pt(&mz, st->xy[2 * k], st->xy[2 * k + 1], &xy[2 * m], &xy[2 * m + 1]);
			svg_polyline(&s, xy, m, "#c9d1d9", 1.0, 0.7, 1);
			strand_svg(&s, &mz, st->fit, "#ff8c42", 2.0, 0.95);
			strand_svg(&s, &mz, st->opt, "#39d0ff", 2.2, 0.95);
		}
		fprintf(s.f, "</g>\n");
		fprintf(s.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"none\" stroke=\"#d0d7de\"/>\n", mz.x0, mz.y0, zpw, zph);

		{
			double tx = mz.x0 + zpw + 30, ty = mz.y0 + 20;
			sprintf(buf, "%d strands: %d traced points -> %u control points", ns, total_points, total_cp);
			svg_text(&s, tx, ty, 15, "#24292f", "start", buf);
			sprintf(buf, "mean angle to the stripes: initial fits %.2f deg -> optimized %.2f deg", a0 / ns, a1 / ns);
			svg_text(&s, tx, ty + 26, 15, "#0969da", "start", buf);
			sprintf(buf, "alignment energy T^T J T / tr J: %.4f -> %.4f", e0 / ns, e1 / ns);
			svg_text(&s, tx, ty + 50, 14, "#57606a", "start", buf);
			svg_text(&s, tx, ty + 86, 13, "#57606a", "start", "pipeline");
			svg_text(&s, tx, ty + 106, 13, "#24292f", "start", "1. structure tensor J = G * (grad I grad I^T), stripes = minor eigenvector");
			svg_text(&s, tx, ty + 126, 13, "#24292f", "start", "2. evenly spaced streamlines (separation 2.2 px)");
			svg_text(&s, tx, ty + 146, 13, "#24292f", "start", "3. least-squares cubic B-spline per streamline (1 control point / 60 px)");
			svg_text(&s, tx, ty + 166, 13, "#24292f", "start", "4. minimize T^T J T / tr J + spring + bending with Adam:");
			svg_text(&s, tx, ty + 186, 13, "#24292f", "start", "   unit tangent T -> qaws_curve_geometry_adjoint_2d -> jet adjoint");
			svg_text(&s, tx, ty + 206, 13, "#24292f", "start", "   -> qaws_curve_eval_batch_adjoint_2d -> control points");
		}
	}
	svg_close(&s);
	printf("1_hair: %d strands, %d points -> %u CPs, angle %.2f -> %.2f deg, energy %.4f -> %.4f\n",
		ns, total_points, total_cp, a0 / ns, a1 / ns, e0 / ns, e1 / ns);
	for (i = 0; i < ns; i++)
	{
		free(strands[i].xy);
		if (strands[i].fit) qaws_curve_destroy(strands[i].fit);
		if (strands[i].opt) qaws_curve_destroy(strands[i].opt);
	}
	free(tf.xx); free(tf.xy); free(tf.yy);
	image_free(&img);
}

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
	printf("2_vectorize: %d curves, %d polyline vertices -> %d points, distance %.3f -> %.3f px\n",
		nc, total_poly, total_pts, d0 / nc, d1 / nc);
	for (i = 0; i < nc; i++)
		free(cs[i].xy);
	free(f);
	image_free(&img);
}

int main(int argc, char** argv)
{
	setvbuf(stdout, NULL, _IONBF, 0);
	if (argc > 1)
		g_photos = argv[1];
	MAKE_DIR("showcase");
	app_hair();
	app_vectorize();
	return 0;
}
