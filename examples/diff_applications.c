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

/* ================================================================== */
/*  3. Triangulated mesh to NURBS patches with ADMM                   */
/* ================================================================== */

#define MP_PATCHES 6
#define MP_N 7                       /* control points per direction */
#define MP_CP (MP_N * MP_N)
#define MP_GRID 26                   /* mesh vertices per face edge */
#define MP_OUTER 40                  /* ADMM iterations */
#define MP_INNER 12                  /* Adam steps per x-update */

typedef struct mesh
{
	int nv, nt;
	double* v;      /* nv * 3 */
	int* t;         /* nt * 3 */
} mesh;

/* A bumpy closed shape (radius as a function of direction). */
static double blob_radius(double x, double y, double z)
{
	double phi = atan2(y, x), th = acos(z / sqrt(x * x + y * y + z * z));
	return 1.0 + 0.16 * sin(3 * phi) * sin(th) * sin(th) + 0.12 * cos(2 * th) + 0.08 * sin(5 * phi + 2 * th) * sin(th);
}

/* Cube face f: axis a = f / 2, sign s, tangent axes b, c. */
static void cube_point(int f, double S, double T, double* p)
{
	int a = f / 2, b = (a + 1) % 3, c = (a + 2) % 3;
	double s = (f & 1) ? -1.0 : 1.0;
	p[a] = s;
	p[b] = (f & 1) ? T : S;   /* flip on negative faces so every face is outward-oriented */
	p[c] = (f & 1) ? S : T;
}

static void mesh_blob(mesh* m)
{
	int f, i, j, k = 0, t = 0;
	m->nv = MP_PATCHES * MP_GRID * MP_GRID;
	m->nt = MP_PATCHES * (MP_GRID - 1) * (MP_GRID - 1) * 2;
	m->v = (double*)malloc(sizeof(double) * 3 * (size_t)m->nv);
	m->t = (int*)malloc(sizeof(int) * 3 * (size_t)m->nt);
	for (f = 0; f < MP_PATCHES; f++)
	{
		for (i = 0; i < MP_GRID; i++)
			for (j = 0; j < MP_GRID; j++)
			{
				double p[3], len, r;
				/* equal-angle cube map for even sampling */
				cube_point(f, tan((i / (double)(MP_GRID - 1) - 0.5) * PI / 2), tan((j / (double)(MP_GRID - 1) - 0.5) * PI / 2), p);
				len = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
				r = blob_radius(p[0], p[1], p[2]);
				m->v[3 * k] = p[0] / len * r;
				m->v[3 * k + 1] = p[1] / len * r;
				m->v[3 * k + 2] = p[2] / len * r;
				k++;
			}
		for (i = 0; i + 1 < MP_GRID; i++)
			for (j = 0; j + 1 < MP_GRID; j++)
			{
				int a = f * MP_GRID * MP_GRID + i * MP_GRID + j;
				m->t[3 * t] = a; m->t[3 * t + 1] = a + MP_GRID; m->t[3 * t + 2] = a + MP_GRID + 1; t++;
				m->t[3 * t] = a; m->t[3 * t + 1] = a + MP_GRID + 1; m->t[3 * t + 2] = a + 1; t++;
			}
	}
}

/* Optional OBJ input (v / f lines, triangles or polygons fanned). */
static int mesh_load_obj(char const* path, mesh* m)
{
	FILE* f = fopen(path, "r");
	char line[512];
	int cv = 0, ct = 0;
	if (!f)
		return 0;
	m->nv = m->nt = 0;
	while (fgets(line, sizeof(line), f))
	{
		if (line[0] == 'v' && line[1] == ' ') m->nv++;
		else if (line[0] == 'f' && line[1] == ' ')
		{
			int n = 0;
			char* p = line + 1;
			while (*p)
			{
				while (*p == ' ') p++;
				if (*p && *p != '\n' && *p != '\r') { n++; while (*p && *p != ' ') p++; }
				else break;
			}
			m->nt += n - 2;
		}
	}
	m->v = (double*)malloc(sizeof(double) * 3 * (size_t)m->nv);
	m->t = (int*)malloc(sizeof(int) * 3 * (size_t)m->nt);
	rewind(f);
	while (fgets(line, sizeof(line), f))
	{
		if (line[0] == 'v' && line[1] == ' ')
		{
			sscanf(line + 2, "%lf %lf %lf", &m->v[3 * cv], &m->v[3 * cv + 1], &m->v[3 * cv + 2]);
			cv++;
		}
		else if (line[0] == 'f' && line[1] == ' ')
		{
			int idx[64], n = 0, k;
			char* p = line + 1;
			while (*p && n < 64)
			{
				while (*p == ' ') p++;
				if (!*p || *p == '\n' || *p == '\r') break;
				idx[n++] = atoi(p) - 1;
				while (*p && *p != ' ') p++;
			}
			for (k = 1; k + 1 < n; k++)
			{
				m->t[3 * ct] = idx[0]; m->t[3 * ct + 1] = idx[k]; m->t[3 * ct + 2] = idx[k + 1];
				ct++;
			}
		}
	}
	fclose(f);
	m->nt = ct;
	{
		/* center and scale to radius ~1 */
		double c[3] = { 0, 0, 0 }, rmax = 0;
		int i, a;
		for (i = 0; i < m->nv; i++)
			for (a = 0; a < 3; a++)
				c[a] += m->v[3 * i + a] / m->nv;
		for (i = 0; i < m->nv; i++)
		{
			double r = 0;
			for (a = 0; a < 3; a++)
			{
				m->v[3 * i + a] -= c[a];
				r += m->v[3 * i + a] * m->v[3 * i + a];
			}
			if (r > rmax) rmax = r;
		}
		rmax = sqrt(rmax);
		for (i = 0; i < 3 * m->nv; i++)
			m->v[i] /= rmax;
	}
	return 1;
}

typedef struct patch_fit
{
	qaws_scalar cps[MP_CP * 3];
	qaws_scalar dual[MP_CP * 3];     /* scaled ADMM duals (boundary nodes) */
	int node[MP_CP];                 /* global node, -1 when interior */
	int n;                           /* assigned mesh points */
	int* pts;
	qaws_scalar* uv;                 /* foot points */
	adam opt;
} patch_fit;

static qaws_scalar const g_mp_knots[MP_N + 4] = { 0, 0, 0, 0, 0.25f, 0.5f, 0.75f, 1, 1, 1, 1 };
static qaws_scalar g_mp_weights[MP_CP];

static qaws_surface* mp_surface(qaws_scalar const* cps)
{
	/* bicubic B-spline: a NURBS patch with unit weights, linear in its
	   control points (exact x-updates, direct thin-plate HVP) */
	qaws_surface_bspline_desc d;
	qaws_surface* s = NULL;
	memset(&d, 0, sizeof(d));
	d.u_degree = 3;
	d.v_degree = 3;
	d.control_points = (qaws_vec3 const*)cps;
	d.u_point_count = MP_N;
	d.v_point_count = MP_N;
	d.u_knots = g_mp_knots;
	d.u_knot_count = MP_N + 4;
	d.v_knots = g_mp_knots;
	d.v_knot_count = MP_N + 4;
	qaws_surface_create_bspline(&d, &s);
	return s;
}

/* Newton steps on the foot point of every assigned point (clamped to the
   patch domain); returns the patch's mean squared distance. */
static double mp_project(patch_fit* p, mesh const* m, qaws_surface const* s, int steps)
{
	double e = 0;
	int i, it;
	for (i = 0; i < p->n; i++)
	{
		double const* x = &m->v[3 * p->pts[i]];
		qaws_scalar* uv = &p->uv[2 * i];
		qaws_surface_jet j;
		double rx, ry, rz;
		for (it = 0; it <= steps; it++)
		{
			double g0, g1, a, b, c, det, du, dv;
			qaws_surface_eval_jet(s, uv[0], uv[1], QAWS_SJET_ORDER2, &j);
			rx = j.d[0].x - x[0]; ry = j.d[0].y - x[1]; rz = j.d[0].z - x[2];
			if (it == steps)
				break;
			g0 = rx * j.d[1].x + ry * j.d[1].y + rz * j.d[1].z;
			g1 = rx * j.d[2].x + ry * j.d[2].y + rz * j.d[2].z;
			a = j.d[1].x * j.d[1].x + j.d[1].y * j.d[1].y + j.d[1].z * j.d[1].z + rx * j.d[3].x + ry * j.d[3].y + rz * j.d[3].z;
			b = j.d[1].x * j.d[2].x + j.d[1].y * j.d[2].y + j.d[1].z * j.d[2].z + rx * j.d[4].x + ry * j.d[4].y + rz * j.d[4].z;
			c = j.d[2].x * j.d[2].x + j.d[2].y * j.d[2].y + j.d[2].z * j.d[2].z + rx * j.d[5].x + ry * j.d[5].y + rz * j.d[5].z;
			det = a * c - b * b;
			if (!(det > 1e-12) || !(a > 0))
				break;
			du = -(c * g0 - b * g1) / det;
			dv = -(a * g1 - b * g0) / det;
			uv[0] = (qaws_scalar)(uv[0] + du < 0 ? 0 : (uv[0] + du > 1 ? 1 : uv[0] + du));
			uv[1] = (qaws_scalar)(uv[1] + dv < 0 ? 0 : (uv[1] + dv > 1 ? 1 : uv[1] + dv));
		}
		e += rx * rx + ry * ry + rz * rz;
	}
	return p->n ? e / p->n : 0;
}
/* Builds the patches, the point assignment and the seam nodes. */
static int mp_setup(mesh const* m, patch_fit* pf, int* node_count)
{
	double keys[MP_PATCHES * MP_CP][3];
	int nkeys = 0, f, i, j, k;
	for (i = 0; i < MP_CP; i++)
		g_mp_weights[i] = 1;
	for (f = 0; f < MP_PATCHES; f++)
	{
		patch_fit* p = &pf[f];
		memset(p, 0, sizeof(*p));
		for (i = 0; i < MP_N; i++)
			for (j = 0; j < MP_N; j++)
			{
				double c[3], len;
				int b = i == 0 || j == 0 || i == MP_N - 1 || j == MP_N - 1, idx = -1;
				cube_point(f, -1 + 2.0 * i / (MP_N - 1), -1 + 2.0 * j / (MP_N - 1), c);
				len = sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
				p->cps[(i * MP_N + j) * 3 + 0] = (qaws_scalar)(c[0] / len);
				p->cps[(i * MP_N + j) * 3 + 1] = (qaws_scalar)(c[1] / len);
				p->cps[(i * MP_N + j) * 3 + 2] = (qaws_scalar)(c[2] / len);
				if (b)
				{
					for (k = 0; k < nkeys; k++)
						if (fabs(keys[k][0] - c[0]) + fabs(keys[k][1] - c[1]) + fabs(keys[k][2] - c[2]) < 1e-9)
							idx = k;
					if (idx < 0)
					{
						idx = nkeys++;
						keys[idx][0] = c[0]; keys[idx][1] = c[1]; keys[idx][2] = c[2];
					}
				}
				p->node[i * MP_N + j] = idx;
			}
	}
	*node_count = nkeys;
	/* assignment by dominant direction, initial foot points from the cube map */
	for (f = 0; f < MP_PATCHES; f++)
	{
		pf[f].pts = (int*)malloc(sizeof(int) * (size_t)m->nv);
		pf[f].uv = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * (size_t)m->nv);
	}
	for (i = 0; i < m->nv; i++)
	{
		double const* x = &m->v[3 * i];
		int a = 0, b, c, ff;
		double s, S, T;
		for (k = 1; k < 3; k++)
			if (fabs(x[k]) > fabs(x[a])) a = k;
		s = x[a] > 0 ? 1 : -1;
		ff = 2 * a + (s < 0);
		b = (a + 1) % 3;
		c = (a + 2) % 3;
		S = x[b] / fabs(x[a]);
		T = x[c] / fabs(x[a]);
		if (ff & 1) { double tmp = S; S = T; T = tmp; }
		pf[ff].pts[pf[ff].n] = i;
		pf[ff].uv[2 * pf[ff].n] = (qaws_scalar)((S + 1) / 2);
		pf[ff].uv[2 * pf[ff].n + 1] = (qaws_scalar)((T + 1) / 2);
		pf[ff].n++;
	}
	return 1;
}

/* Largest distance between copies of a seam node. */
static double mp_seam_gap(patch_fit const* pf, int nodes)
{
	double gap = 0;
	int g, f, k, f2, k2;
	for (g = 0; g < nodes; g++)
		for (f = 0; f < MP_PATCHES; f++)
			for (k = 0; k < MP_CP; k++)
				if (pf[f].node[k] == g)
					for (f2 = f + 1; f2 < MP_PATCHES; f2++)
						for (k2 = 0; k2 < MP_CP; k2++)
							if (pf[f2].node[k2] == g)
							{
								double dx = pf[f].cps[3 * k] - pf[f2].cps[3 * k2];
								double dy = pf[f].cps[3 * k + 1] - pf[f2].cps[3 * k2 + 1];
								double dz = pf[f].cps[3 * k + 2] - pf[f2].cps[3 * k2 + 2];
								double d = sqrt(dx * dx + dy * dy + dz * dz);
								if (d > gap) gap = d;
							}
	return gap;
}

/* Thin-plate Hessian of a patch (one coordinate; the energy decouples
   over x, y, z), assembled column by column from direct HVPs. */
static void mp_thin_plate_hessian(double* H)
{
	qaws_scalar cps[MP_CP * 3], dir[MP_CP * 3], out[MP_CP * 3];
	qaws_surface* s;
	int k, j;
	memset(cps, 0, sizeof(cps));
	s = mp_surface(cps);
	for (k = 0; k < MP_CP; k++)
	{
		qaws_field_view fd, fo;
		qaws_diff_views vd, vo;
		memset(dir, 0, sizeof(dir));
		memset(out, 0, sizeof(out));
		dir[3 * k] = 1;
		vd = one_field(&fd, QAWS_FIELD_CONTROL_POINTS, dir, MP_CP, 3);
		vo = one_field(&fo, QAWS_FIELD_CONTROL_POINTS, out, MP_CP, 3);
		qaws_surface_functional_hvp(NULL, s, QAWS_FUNCTIONAL_THIN_PLATE, 4, &vd, &vo);
		for (j = 0; j < MP_CP; j++)
			H[j * MP_CP + k] = out[3 * j];
	}
	qaws_surface_destroy(s);
}

/* In-place Cholesky solve of an SPD n x n system for three right-hand sides. */
static void mp_cholesky_solve(double* A, double b[3][MP_CP], int n)
{
	int i, j, k, c;
	for (i = 0; i < n; i++)
		for (j = 0; j <= i; j++)
		{
			double v = A[i * n + j];
			for (k = 0; k < j; k++)
				v -= A[i * n + k] * A[j * n + k];
			A[i * n + j] = (i == j) ? sqrt(v > 1e-18 ? v : 1e-18) : v / A[j * n + j];
		}
	for (c = 0; c < 3; c++)
	{
		for (i = 0; i < n; i++)
		{
			double v = b[c][i];
			for (k = 0; k < i; k++)
				v -= A[i * n + k] * b[c][k];
			b[c][i] = v / A[i * n + i];
		}
		for (i = n; i-- > 0;)
		{
			double v = b[c][i];
			for (k = i + 1; k < n; k++)
				v -= A[k * n + i] * b[c][k];
			b[c][i] = v / A[i * n + i];
		}
	}
}

/* Exact x-update at fixed foot points:
   (2 A + tp H + rho D) P = 2 b + rho D (z - u), with A, b from the exact
   basis weights of every foot point (qaws_surface_local_support). */
static void mp_xupdate(patch_fit* p, mesh const* m, double const* H, double tp, double rho, qaws_scalar const* z)
{
	static double A[MP_CP * MP_CP];
	double b[3][MP_CP];
	qaws_surface* s = mp_surface(p->cps);
	int i, a, c, k;
	memset(A, 0, sizeof(A));
	memset(b, 0, sizeof(b));
	for (i = 0; i < p->n; i++)
	{
		qaws_surface_support sup;
		double const* x = &m->v[3 * p->pts[i]];
		int idx[16], n = 0;
		double w[16];
		qaws_surface_local_support(s, p->uv[2 * i], p->uv[2 * i + 1], 0, &sup);
		for (a = 0; a < (int)sup.u_count; a++)
			for (c = 0; c < (int)sup.v_count; c++)
			{
				idx[n] = (int)((sup.u_first + a) * sup.u_stride + (sup.v_first + c) * sup.v_stride);
				w[n++] = sup.u_weights[0][a] * sup.v_weights[0][c];
			}
		for (a = 0; a < n; a++)
		{
			for (c = 0; c < n; c++)
				A[idx[a] * MP_CP + idx[c]] += 2 * w[a] * w[c] / p->n;
			for (k = 0; k < 3; k++)
				b[k][idx[a]] += 2 * w[a] * x[k] / p->n;
		}
	}
	for (k = 0; k < MP_CP * MP_CP; k++)
		A[k] += tp * H[k];
	if (rho > 0)
		for (k = 0; k < MP_CP; k++)
			if (p->node[k] >= 0)
			{
				A[k * MP_CP + k] += rho;
				for (c = 0; c < 3; c++)
					b[c][k] += rho * (z[3 * p->node[k] + c] - p->dual[3 * k + c]);
			}
	mp_cholesky_solve(A, b, MP_CP);
	for (k = 0; k < MP_CP; k++)
		for (c = 0; c < 3; c++)
			p->cps[3 * k + c] = (qaws_scalar)b[c][k];
	qaws_surface_destroy(s);
}

/* ADMM on the seam copies: foot points, x-updates (independent per patch),
   z = mean of (x + u) over copies, u += x - z. rho = 0 gives the
   independent fits. */
static void mp_solve(mesh const* m, patch_fit* pf, int nodes, double rho, double* rms, double* gap)
{
	static double H[MP_CP * MP_CP];
	double tp = 2e-6;
	qaws_scalar* z = (qaws_scalar*)calloc((size_t)nodes * 3, sizeof(qaws_scalar));
	int* cnt = (int*)calloc((size_t)nodes, sizeof(int));
	int it, f, k, c;

	mp_thin_plate_hessian(H);
	for (f = 0; f < MP_PATCHES; f++)
		for (k = 0; k < MP_CP; k++)
			if (pf[f].node[k] >= 0)
				for (c = 0; c < 3; c++)
					z[3 * pf[f].node[k] + c] += pf[f].cps[3 * k + c];
	for (f = 0; f < MP_PATCHES; f++)
		for (k = 0; k < MP_CP; k++)
			if (pf[f].node[k] >= 0)
				cnt[pf[f].node[k]]++;
	for (k = 0; k < nodes; k++)
		for (c = 0; c < 3; c++)
			z[3 * k + c] /= cnt[k];
	for (it = 0; it < MP_OUTER; it++)
	{
		double e = 0;
		int npts = 0;
		for (f = 0; f < MP_PATCHES; f++)
		{
			qaws_surface* s = mp_surface(pf[f].cps);
			e += mp_project(&pf[f], m, s, 3) * pf[f].n;
			npts += pf[f].n;
			qaws_surface_destroy(s);
		}
		rms[it] = sqrt(e / npts);
		for (f = 0; f < MP_PATCHES; f++)
			mp_xupdate(&pf[f], m, H, tp, rho, z);
		memset(z, 0, sizeof(qaws_scalar) * (size_t)nodes * 3);
		for (f = 0; f < MP_PATCHES; f++)
			for (k = 0; k < MP_CP; k++)
				if (pf[f].node[k] >= 0)
					for (c = 0; c < 3; c++)
						z[3 * pf[f].node[k] + c] += (pf[f].cps[3 * k + c] + pf[f].dual[3 * k + c]) / cnt[pf[f].node[k]];
		if (rho > 0)
			for (f = 0; f < MP_PATCHES; f++)
				for (k = 0; k < MP_CP; k++)
					if (pf[f].node[k] >= 0)
						for (c = 0; c < 3; c++)
							pf[f].dual[3 * k + c] += pf[f].cps[3 * k + c] - z[3 * pf[f].node[k] + c];
		gap[it] = mp_seam_gap(pf, nodes);
	}
	free(z);
	free(cnt);
}

/* --- rendering ---------------------------------------------------- */

typedef struct view3
{
	double cx, cy, scale, yaw, pitch;
} view3;

static void view_xform(view3 const* v, double const* p, double* sx, double* sy, double* depth)
{
	double cy = cos(v->yaw), sy_ = sin(v->yaw), cp = cos(v->pitch), sp = sin(v->pitch);
	double x = cy * p[0] - sy_ * p[1];
	double y = sy_ * p[0] + cy * p[1];
	double z = p[2];
	double y2 = cp * y - sp * z, z2 = sp * y + cp * z;
	*sx = v->cx + x * v->scale;
	*sy = v->cy - z2 * v->scale;
	*depth = y2;   /* larger = farther */
}

typedef struct poly3
{
	double xy[8];
	int n;
	double depth;
	char color[40];
} poly3;

static int poly_cmp(void const* a, void const* b)
{
	double da = ((poly3 const*)a)->depth, db = ((poly3 const*)b)->depth;
	return da < db ? 1 : (da > db ? -1 : 0);
}

static void shade(double const* n, double const* base, char* out)
{
	double l[3] = { -0.45, -0.55, 0.70 }, d = n[0] * l[0] + n[1] * l[1] + n[2] * l[2], k;
	double ln = sqrt(l[0] * l[0] + l[1] * l[1] + l[2] * l[2]);
	d /= ln;
	k = 0.30 + 0.70 * (d > 0 ? d : 0);
	sprintf(out, "rgb(%d,%d,%d)", (int)(255 * base[0] * k), (int)(255 * base[1] * k), (int)(255 * base[2] * k));
}

static void draw_polys(svg* s, poly3* polys, int n, double stroke)
{
	int i;
	qsort(polys, (size_t)n, sizeof(poly3), poly_cmp);
	for (i = 0; i < n; i++)
	{
		poly3 const* p = &polys[i];
		if (p->n == 3)
			fprintf(s->f, "<polygon points=\"%.1f,%.1f %.1f,%.1f %.1f,%.1f\" fill=\"%s\" stroke=\"%s\" stroke-width=\"%.2f\"/>\n",
				p->xy[0], p->xy[1], p->xy[2], p->xy[3], p->xy[4], p->xy[5], p->color, p->color, stroke);
		else
			fprintf(s->f, "<polygon points=\"%.1f,%.1f %.1f,%.1f %.1f,%.1f %.1f,%.1f\" fill=\"%s\" stroke=\"%s\" stroke-width=\"%.2f\"/>\n",
				p->xy[0], p->xy[1], p->xy[2], p->xy[3], p->xy[4], p->xy[5], p->xy[6], p->xy[7], p->color, p->color, stroke);
	}
}

static void render_mesh(svg* s, view3 const* v, mesh const* m)
{
	poly3* polys = (poly3*)malloc(sizeof(poly3) * (size_t)m->nt);
	int t, n = 0;
	double base[3] = { 0.78, 0.80, 0.84 };
	for (t = 0; t < m->nt; t++)
	{
		double const* a = &m->v[3 * m->t[3 * t]];
		double const* b = &m->v[3 * m->t[3 * t + 1]];
		double const* c = &m->v[3 * m->t[3 * t + 2]];
		double e1[3], e2[3], nn[3], len, d0, d1, d2, cen[3];
		int k;
		for (k = 0; k < 3; k++) { e1[k] = b[k] - a[k]; e2[k] = c[k] - a[k]; cen[k] = (a[k] + b[k] + c[k]) / 3; }
		nn[0] = e1[1] * e2[2] - e1[2] * e2[1];
		nn[1] = e1[2] * e2[0] - e1[0] * e2[2];
		nn[2] = e1[0] * e2[1] - e1[1] * e2[0];
		len = sqrt(nn[0] * nn[0] + nn[1] * nn[1] + nn[2] * nn[2]) + 1e-12;
		if (nn[0] * cen[0] + nn[1] * cen[1] + nn[2] * cen[2] < 0) len = -len;
		for (k = 0; k < 3; k++) nn[k] /= len;
		view_xform(v, a, &polys[n].xy[0], &polys[n].xy[1], &d0);
		view_xform(v, b, &polys[n].xy[2], &polys[n].xy[3], &d1);
		view_xform(v, c, &polys[n].xy[4], &polys[n].xy[5], &d2);
		{
			/* back faces away from the viewer are skipped */
			double vd[3] = { 0, 0, 0 }, cy = cos(v->yaw), sy_ = sin(v->yaw), cp = cos(v->pitch), sp = sin(v->pitch);
			double y = sy_ * nn[0] + cy * nn[1];
			if (cp * y - sp * nn[2] > 0.02)
				continue;
			(void)vd;
		}
		polys[n].n = 3;
		polys[n].depth = (d0 + d1 + d2) / 3;
		shade(nn, base, polys[n].color);
		n++;
	}
	draw_polys(s, polys, n, 0.4);
	free(polys);
}

static void render_patches(svg* s, view3 const* v, patch_fit const* pf, int grid, int heat_error, mesh const* m)
{
	static double const colors[MP_PATCHES][3] = {
		{ 0.35, 0.55, 0.95 }, { 0.95, 0.55, 0.30 }, { 0.40, 0.80, 0.45 },
		{ 0.85, 0.40, 0.70 }, { 0.95, 0.80, 0.30 }, { 0.40, 0.80, 0.85 } };
	poly3* polys = (poly3*)malloc(sizeof(poly3) * MP_PATCHES * (size_t)grid * grid);
	int f, i, j, n = 0;
	(void)heat_error;
	(void)m;
	for (f = 0; f < MP_PATCHES; f++)
	{
		qaws_surface* srf = mp_surface(pf[f].cps);
		for (i = 0; i < grid; i++)
			for (j = 0; j < grid; j++)
			{
				qaws_surface_eval_result r[4], mid;
				double p[3], d, dsum = 0, nn[3], cy = cos(v->yaw), sy_ = sin(v->yaw), cp = cos(v->pitch), sp = sin(v->pitch);
				int c;
				qaws_scalar uu[4] = { (qaws_scalar)(i / (double)grid), (qaws_scalar)((i + 1) / (double)grid),
					(qaws_scalar)((i + 1) / (double)grid), (qaws_scalar)(i / (double)grid) };
				qaws_scalar vv[4] = { (qaws_scalar)(j / (double)grid), (qaws_scalar)(j / (double)grid),
					(qaws_scalar)((j + 1) / (double)grid), (qaws_scalar)((j + 1) / (double)grid) };
				for (c = 0; c < 4; c++)
				{
					qaws_surface_evaluate(srf, uu[c], vv[c], QAWS_SURFACE_EVAL_POSITION, &r[c]);
					p[0] = r[c].position.x; p[1] = r[c].position.y; p[2] = r[c].position.z;
					view_xform(v, p, &polys[n].xy[2 * c], &polys[n].xy[2 * c + 1], &d);
					dsum += d;
				}
				qaws_surface_evaluate(srf, (qaws_scalar)((i + 0.5) / grid), (qaws_scalar)((j + 0.5) / grid),
					QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL, &mid);
				nn[0] = mid.normal.x; nn[1] = mid.normal.y; nn[2] = mid.normal.z;
				if (nn[0] * mid.position.x + nn[1] * mid.position.y + nn[2] * mid.position.z < 0)
				{
					nn[0] = -nn[0]; nn[1] = -nn[1]; nn[2] = -nn[2];
				}
				if (cp * (sy_ * nn[0] + cy * nn[1]) - sp * nn[2] > 0.02)
					continue;
				polys[n].n = 4;
				polys[n].depth = dsum / 4;
				shade(nn, colors[f], polys[n].color);
				n++;
			}
		qaws_surface_destroy(srf);
	}
	draw_polys(s, polys, n, 0.5);
	free(polys);
}

static void app_mesh_patches(char const* obj_path)
{
	static patch_fit indep[MP_PATCHES], admm[MP_PATCHES];
	double rms0[MP_OUTER], gap0[MP_OUTER], rms1[MP_OUTER], gap1[MP_OUTER];
	mesh m;
	int nodes, f;
	svg s;
	char buf[256];

	if (!obj_path || !mesh_load_obj(obj_path, &m))
		mesh_blob(&m);
	mp_setup(&m, indep, &nodes);
	mp_setup(&m, admm, &nodes);
	mp_solve(&m, indep, nodes, 0.0, rms0, gap0);
	mp_solve(&m, admm, nodes, 0.03, rms1, gap1);

	svg_open(&s, "showcase/app3_mesh_patches.svg", 1240, 800, "Triangulated mesh to 6 bicubic patches with ADMM",
		"Each patch fits its points independently (exact least squares at foot points + thin plate); shared seam "
		"control points are ADMM consensus variables, so the patch set closes up without a global solve.");
	{
		view3 va = { 215, 320, 135, 0.65, 0.42 }, vb = va, vc = va;
		viewport lv = { 40, 560, 520, 210, 0, 0, 0, 0 };
		vb.cx = 620;
		vc.cx = 1025;
		fprintf(s.f, "<text x=\"20\" y=\"92\" font-size=\"13\" font-weight=\"600\" fill=\"#24292f\">input mesh: %d vertices, %d triangles</text>\n", m.nv, m.nt);
		render_mesh(&s, &va, &m);
		fprintf(s.f, "<text x=\"425\" y=\"92\" font-size=\"13\" font-weight=\"600\" fill=\"#cf222e\">independent patch fits: open, overlapping seams</text>\n");
		render_patches(&s, &vb, indep, 22, 0, &m);
		fprintf(s.f, "<text x=\"830\" y=\"92\" font-size=\"13\" font-weight=\"600\" fill=\"#0969da\">ADMM: %d control points, closed seams</text>\n", MP_PATCHES * MP_CP);
		render_patches(&s, &vc, admm, 22, 0, &m);

		svg_panel(&s, &lv, "log10 seam gap (red: independent, blue: ADMM) and RMS (dashed)");
		{
			double xy0[2 * MP_OUTER], xy1[2 * MP_OUTER], xr0[2 * MP_OUTER], xr1[2 * MP_OUTER], lo = 1e300, hi = -1e300;
			int it;
			for (it = 0; it < MP_OUTER; it++)
			{
				double v[4] = { gap0[it], gap1[it], rms0[it], rms1[it] }, l;
				int q;
				for (q = 0; q < 4; q++)
				{
					l = log10(v[q] > 1e-12 ? v[q] : 1e-12);
					if (l < lo) lo = l;
					if (l > hi) hi = l;
				}
			}
			for (it = 0; it < MP_OUTER; it++)
			{
				double x = lv.x0 + 14 + (lv.w - 28) * it / (double)(MP_OUTER - 1);
				xy0[2 * it] = xy1[2 * it] = xr0[2 * it] = xr1[2 * it] = x;
				xy0[2 * it + 1] = lv.y0 + 34 + (lv.h - 54) * (hi - log10(gap0[it] > 1e-12 ? gap0[it] : 1e-12)) / (hi - lo);
				xy1[2 * it + 1] = lv.y0 + 34 + (lv.h - 54) * (hi - log10(gap1[it] > 1e-12 ? gap1[it] : 1e-12)) / (hi - lo);
				xr0[2 * it + 1] = lv.y0 + 34 + (lv.h - 54) * (hi - log10(rms0[it])) / (hi - lo);
				xr1[2 * it + 1] = lv.y0 + 34 + (lv.h - 54) * (hi - log10(rms1[it])) / (hi - lo);
			}
			svg_polyline(&s, xr0, MP_OUTER, "#cf222e", 1.6, 0.8, 1);
			svg_polyline(&s, xr1, MP_OUTER, "#0969da", 1.6, 0.8, 1);
			svg_polyline(&s, xy0, MP_OUTER, "#cf222e", 2.2, 1, 0);
			svg_polyline(&s, xy1, MP_OUTER, "#0969da", 2.6, 1, 0);
			svg_text(&s, lv.x0 + lv.w - 14, lv.y0 + lv.h - 8, 11, "#57606a", "end", "ADMM iteration");
		}
		{
			double tx = 600, ty = 585;
			sprintf(buf, "RMS distance to the mesh: independent %.4f, ADMM %.4f (shape radius ~1)", rms0[MP_OUTER - 1], rms1[MP_OUTER - 1]);
			svg_text(&s, tx, ty, 14, "#24292f", "start", buf);
			sprintf(buf, "max seam gap: independent %.4f, ADMM %.1e", gap0[MP_OUTER - 1], gap1[MP_OUTER - 1]);
			svg_text(&s, tx, ty + 24, 14, "#0969da", "start", buf);
			svg_text(&s, tx, ty + 58, 13, "#57606a", "start", "per ADMM iteration, every patch in parallel:");
			svg_text(&s, tx, ty + 78, 13, "#24292f", "start", "  foot points: Newton on |S(u,v) - X|^2 (surface jets)");
			svg_text(&s, tx, ty + 98, 13, "#24292f", "start", "  x-update: (2A + t H + rho D) P = 2b + rho D (z - u)");
			svg_text(&s, tx, ty + 118, 13, "#24292f", "start", "    A, b from exact basis weights (qaws_surface_local_support)");
			svg_text(&s, tx, ty + 138, 13, "#24292f", "start", "    H: thin-plate Hessian from direct HVPs (qaws_surface_functional_hvp)");
			svg_text(&s, tx, ty + 158, 13, "#57606a", "start", "then z = mean of (x + u) over seam copies, u += x - z");
		}
	}
	svg_close(&s);
	printf("3_mesh_patches: %d points, rms indep %.4f admm %.4f, seam gap indep %.4f admm %.2e\n",
		m.nv, rms0[MP_OUTER - 1], rms1[MP_OUTER - 1], gap0[MP_OUTER - 1], gap1[MP_OUTER - 1]);
	for (f = 0; f < MP_PATCHES; f++)
	{
		free(indep[f].pts); free(indep[f].uv);
		free(admm[f].pts); free(admm[f].uv);
	}
	free(m.v);
	free(m.t);
}

int main(int argc, char** argv)
{
	setvbuf(stdout, NULL, _IONBF, 0);
	if (argc > 1)
		g_photos = argv[1];
	MAKE_DIR("showcase");
	app_hair();
	app_vectorize();
	app_mesh_patches(argc > 2 ? argv[2] : NULL);
	return 0;
}
