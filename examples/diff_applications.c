/*
 * Applications of the differentiation API on real data.
 *
 *   1. Hair strands from a photo: orientation field (structure tensor),
 *      evenly spaced streamlines, B-spline strands optimized to follow the
 *      field through unit-tangent adjoints (qaws_diff_geometry.h).
 *   2. Photo vectorization: isocontours become interpolating splines
 *      (centripetal Catmull-Rom and Yuksel C2) whose points are optimized
 *      onto the contours; curvature combs compare the two families.
 *   3. Triangulated mesh to six bicubic patches with ADMM: exact per-patch
 *      least squares at foot points, seam control points in consensus.
 *   4. Non-rigid registration of a template curve to a scan: CMA-ES on the
 *      pose with gradient refinements inside the fitness (basin hopping).
 *   5. Hair grooming on a head: 3D B-spline strands with fixed roots under
 *      length, bending, gravity, collision and guide-alignment energies.
 *
 *   6. Hair from a photo: strands keep their 3D priors while their
 *      projections follow the photo orientation field (photos/hair.ppm,
 *      or the ribbons photo as the flow source).
 *
 * QAWS_APP=n in the environment runs only application n.
 *
 * Photos are read as binary PPM; examples/photo_to_ppm.ps1 converts any
 * image (and writes a PNG used as figure background): ribbons.ppm and
 * folds.ppm in the photos directory. Figures are written to showcase/ in
 * the working directory.
 *
 *   qaws_diff_applications [photos directory] [mesh.obj]
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

/* ================================================================== */
/*  4. Non-rigid registration: CMA-ES basin hopping + gradients       */
/* ================================================================== */

#define RG_N 16                /* template interpolation points */
#define RG_SCAN 120            /* scan points */
#define RG_LAMBDA 10           /* CMA-ES population */
#define RG_GEN 18
#define RG_LOCAL 10            /* gradient steps inside the fitness */
#define RG_FINAL 250
#define RG_TRIALS 16

typedef struct reg_problem
{
	qaws_scalar tmpl[RG_N * 2];        /* template points */
	double scan[RG_SCAN * 2];
	int nscan;
	double cx, cy;                     /* scan centroid */
} reg_problem;

/* pose: theta, log scale, tx, ty; plus per-point offsets (non-rigid) */
typedef struct reg_state
{
	double pose[4];
	double off[RG_N * 2];
	double foot[RG_SCAN];      /* foot parameters of the scan points */
} reg_state;

static double rg_template_radius(double phi)
{
	return 1.0 + 0.30 * cos(2 * phi + 0.4) + 0.20 * sin(3 * phi) + 0.12 * cos(5 * phi + 1.0);
}

static void rg_points(reg_problem const* pb, reg_state const* st, qaws_scalar* out)
{
	double c = cos(st->pose[0]), s = sin(st->pose[0]), sc = exp(st->pose[1]);
	int i;
	for (i = 0; i < RG_N; i++)
	{
		double x = pb->tmpl[2 * i] + st->off[2 * i], y = pb->tmpl[2 * i + 1] + st->off[2 * i + 1];
		out[2 * i] = (qaws_scalar)(sc * (c * x - s * y) + st->pose[2]);
		out[2 * i + 1] = (qaws_scalar)(sc * (s * x + c * y) + st->pose[3]);
	}
}

static qaws_curve* rg_curve(qaws_scalar const* pts)
{
	qaws_catmull_rom_desc d;
	qaws_curve* c = NULL;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D;
	d.control_points = pts;
	d.control_point_count = RG_N;
	d.parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL;
	d.closed = 1;
	qaws_curve_create_catmull_rom(&d, &c);
	return c;
}

/* Foot points: nearest of 8 samples per span (global), refined later by
   Newton. Called when the pose jumps (CMA-ES candidates). */
static void rg_seed(reg_problem const* pb, reg_state* st)
{
	enum { M = RG_N * 8 };
	qaws_scalar pts[RG_N * 2], ts[M];
	qaws_curve_jet_2d prim[M], tan[M];
	qaws_curve* c;
	int k, j;
	rg_points(pb, st, pts);
	c = rg_curve(pts);
	for (j = 0; j < M; j++)
		ts[j] = (qaws_scalar)(RG_N * (j + 0.5) / M);
	qaws_curve_eval_batch_tangent_2d(NULL, c, ts, NULL, M, QAWS_EVAL_FLAG_POSITION, NULL, prim, tan);
	for (k = 0; k < pb->nscan; k++)
	{
		double best = 1e300;
		for (j = 0; j < M; j++)
		{
			double dx = prim[j].d[0].x - pb->scan[2 * k], dy = prim[j].d[0].y - pb->scan[2 * k + 1];
			if (dx * dx + dy * dy < best)
			{
				best = dx * dx + dy * dy;
				st->foot[k] = ts[j];
			}
		}
	}
	qaws_curve_destroy(c);
}

/* E = mean squared distance of the scan to the deformed template at the
   foot points + smoothness of the offsets. The foot points take two Newton
   steps first; at a foot point the distance gradient is the pullback of
   2 (C - q) alone, through the Catmull-Rom adjoint onto the points, then
   onto pose and offsets. */
static double rg_energy(reg_problem const* pb, reg_state* st, reg_state* g, int nonrigid)
{
	qaws_scalar pts[RG_N * 2], pbar[RG_N * 2], ts[RG_SCAN];
	qaws_curve_jet_2d prim[RG_SCAN], tan[RG_SCAN], bar[RG_SCAN];
	qaws_curve* c;
	qaws_field_view fv;
	qaws_diff_views views = one_field(&fv, QAWS_FIELD_POINTS, pbar, RG_N, 2);
	double e = 0, lam = 0.4;
	int k, i, it;
	rg_points(pb, st, pts);
	c = rg_curve(pts);
	for (it = 0; it < 3; it++)
	{
		for (k = 0; k < pb->nscan; k++)
			ts[k] = (qaws_scalar)st->foot[k];
		qaws_curve_eval_batch_tangent_2d(NULL, c, ts, NULL, (unsigned int)pb->nscan,
			QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2, NULL, prim, tan);
		if (it == 2)
			break;
		for (k = 0; k < pb->nscan; k++)
		{
			double rx = prim[k].d[0].x - pb->scan[2 * k], ry = prim[k].d[0].y - pb->scan[2 * k + 1];
			double gv = rx * prim[k].d[1].x + ry * prim[k].d[1].y;
			double h = prim[k].d[1].x * prim[k].d[1].x + prim[k].d[1].y * prim[k].d[1].y + rx * prim[k].d[2].x + ry * prim[k].d[2].y;
			double t = st->foot[k] - (h > 1e-9 ? gv / h : 0);
			if (t - st->foot[k] > 0.5) t = st->foot[k] + 0.5;
			if (t - st->foot[k] < -0.5) t = st->foot[k] - 0.5;
			t = fmod(t, (double)RG_N);
			if (t < 0) t += RG_N;
			st->foot[k] = t;
		}
	}
	for (k = 0; k < pb->nscan; k++)
	{
		double rx = prim[k].d[0].x - pb->scan[2 * k], ry = prim[k].d[0].y - pb->scan[2 * k + 1];
		e += (rx * rx + ry * ry) / pb->nscan;
		memset(&bar[k], 0, sizeof(bar[k]));
		bar[k].d[0].x = (qaws_scalar)(2 * rx / pb->nscan);
		bar[k].d[0].y = (qaws_scalar)(2 * ry / pb->nscan);
		bar[k].channels = QAWS_EVAL_FLAG_POSITION;
	}
	for (i = 0; i < RG_N; i++)
	{
		int j = (i + 1) % RG_N;
		double dx = st->off[2 * i] - st->off[2 * j], dy = st->off[2 * i + 1] - st->off[2 * j + 1];
		e += lam * (dx * dx + dy * dy) / RG_N;
	}
	if (g)
	{
		double cth = cos(st->pose[0]), sth = sin(st->pose[0]), sc = exp(st->pose[1]);
		memset(pbar, 0, sizeof(pbar));
		qaws_curve_eval_batch_adjoint_2d(NULL, c, ts, (unsigned int)pb->nscan, QAWS_EVAL_FLAG_POSITION, bar, &views, NULL);
		memset(g->pose, 0, sizeof(g->pose));
		memset(g->off, 0, sizeof(g->off));
		for (i = 0; i < RG_N; i++)
		{
			double x = pb->tmpl[2 * i] + st->off[2 * i], y = pb->tmpl[2 * i + 1] + st->off[2 * i + 1];
			double gx = pbar[2 * i], gy = pbar[2 * i + 1];
			/* P' = s R p + t */
			g->pose[0] += sc * (gx * (-sth * x - cth * y) + gy * (cth * x - sth * y));
			g->pose[1] += sc * (gx * (cth * x - sth * y) + gy * (sth * x + cth * y));
			g->pose[2] += gx;
			g->pose[3] += gy;
			if (nonrigid)
			{
				int j = (i + 1) % RG_N, h = (i + RG_N - 1) % RG_N;
				g->off[2 * i] += sc * (cth * gx + sth * gy);
				g->off[2 * i + 1] += sc * (-sth * gx + cth * gy);
				g->off[2 * i] += lam * 2 * (2 * st->off[2 * i] - st->off[2 * j] - st->off[2 * h]) / RG_N;
				g->off[2 * i + 1] += lam * 2 * (2 * st->off[2 * i + 1] - st->off[2 * j + 1] - st->off[2 * h + 1]) / RG_N;
			}
		}
	}
	qaws_curve_destroy(c);
	return e;
}

/* Local refinement with Adam; returns the final energy. */
static double rg_refine(reg_problem const* pb, reg_state* st, int steps, int nonrigid, double lr)
{
	adam a;
	int it, k;
	double e = 0;
	memset(&a, 0, sizeof(a));
	for (it = 0; it <= steps; it++)
	{
		reg_state g;
		qaws_scalar x[4 + RG_N * 2], gg[4 + RG_N * 2];
		e = rg_energy(pb, st, it < steps ? &g : NULL, nonrigid);
		if (it == steps)
			break;
		for (k = 0; k < 4; k++) { x[k] = (qaws_scalar)st->pose[k]; gg[k] = (qaws_scalar)g.pose[k]; }
		for (k = 0; k < RG_N * 2; k++) { x[4 + k] = (qaws_scalar)st->off[k]; gg[4 + k] = (qaws_scalar)g.off[k]; }
		adam_step(&a, x, gg, nonrigid ? 4 + RG_N * 2 : 4, lr * (1.0 - 0.8 * it / (double)steps));
		for (k = 0; k < 4; k++) st->pose[k] = x[k];
		if (nonrigid)
			for (k = 0; k < RG_N * 2; k++) st->off[k] = x[4 + k];
	}
	return e;
}

/* --- CMA-ES (4D pose, rank-one and rank-mu updates) ---------------- */

static double rg_gauss(unsigned int* rng)
{
	double u1, u2;
	*rng = *rng * 1664525u + 1013904223u;
	u1 = ((*rng >> 8) + 0.5) / 16777216.0;
	*rng = *rng * 1664525u + 1013904223u;
	u2 = ((*rng >> 8) + 0.5) / 16777216.0;
	return sqrt(-2 * log(u1)) * cos(2 * PI * u2);
}

/* Jacobi eigen decomposition of a symmetric 4x4: C = B diag(d) B^T. */
static void rg_eigen(double C[4][4], double B[4][4], double* d)
{
	double A[4][4];
	int i, j, k, sweep;
	memcpy(A, C, sizeof(A));
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			B[i][j] = i == j;
	for (sweep = 0; sweep < 30; sweep++)
		for (i = 0; i < 4; i++)
			for (j = i + 1; j < 4; j++)
			{
				double th, c, s, t;
				if (fabs(A[i][j]) < 1e-15)
					continue;
				th = (A[j][j] - A[i][i]) / (2 * A[i][j]);
				t = (th >= 0 ? 1 : -1) / (fabs(th) + sqrt(th * th + 1));
				c = 1 / sqrt(t * t + 1);
				s = t * c;
				for (k = 0; k < 4; k++)
				{
					double aki = A[k][i], akj = A[k][j];
					A[k][i] = c * aki - s * akj;
					A[k][j] = s * aki + c * akj;
				}
				for (k = 0; k < 4; k++)
				{
					double aik = A[i][k], ajk = A[j][k];
					A[i][k] = c * aik - s * ajk;
					A[j][k] = s * aik + c * ajk;
				}
				for (k = 0; k < 4; k++)
				{
					double bki = B[k][i], bkj = B[k][j];
					B[k][i] = c * bki - s * bkj;
					B[k][j] = s * bki + c * bkj;
				}
			}
	for (i = 0; i < 4; i++)
		d[i] = A[i][i] > 1e-20 ? A[i][i] : 1e-20;
}

typedef struct cma_trace
{
	double best[RG_GEN];
	double cand[RG_GEN][RG_LAMBDA][4];
} cma_trace;

/* Minimizes the refined energy over the pose; scaled coordinates
   y = (pose - pose0) / scale. Returns the best pose. */
static void rg_cmaes(reg_problem const* pb, double const* pose0, double const* scale, unsigned int seed, double* best_pose,
	cma_trace* trace)
{
	enum { N = 4, MU = RG_LAMBDA / 2 };
	double m[N] = { 0, 0, 0, 0 }, sigma = 0.6, C[N][N], B[N][N], D[N], pc[N] = { 0 }, ps[N] = { 0 };
	double w[MU], wsum = 0, mueff, cc, cs, c1, cmu, damps, chi = sqrt((double)N) * (1 - 1.0 / (4 * N) + 1.0 / (21 * N * N));
	double best_f = 1e300;
	unsigned int rng = seed;
	int i, j, k, gen;
	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
			C[i][j] = i == j;
	for (i = 0; i < MU; i++)
	{
		w[i] = log(MU + 0.5) - log(i + 1.0);
		wsum += w[i];
	}
	mueff = 0;
	for (i = 0; i < MU; i++)
	{
		w[i] /= wsum;
		mueff += w[i] * w[i];
	}
	mueff = 1 / mueff;
	cc = (4 + mueff / N) / (N + 4 + 2 * mueff / N);
	cs = (mueff + 2) / (N + mueff + 5);
	c1 = 2 / ((N + 1.3) * (N + 1.3) + mueff);
	cmu = 2 * (mueff - 2 + 1 / mueff) / ((N + 2) * (N + 2) + mueff);
	if (cmu > 1 - c1) cmu = 1 - c1;
	damps = 1 + 2 * (mueff > N + 1 ? sqrt((mueff - 1) / (N + 1)) - 1 : 0) + cs;
	memcpy(best_pose, pose0, sizeof(double) * N);

	for (gen = 0; gen < RG_GEN; gen++)
	{
		double y[RG_LAMBDA][N], z[RG_LAMBDA][N], f[RG_LAMBDA], mold[N];
		int order[RG_LAMBDA];
		rg_eigen(C, B, D);
		for (k = 0; k < RG_LAMBDA; k++)
		{
			reg_state st;
			for (i = 0; i < N; i++)
				z[k][i] = rg_gauss(&rng);
			for (i = 0; i < N; i++)
			{
				double v = 0;
				for (j = 0; j < N; j++)
					v += B[i][j] * sqrt(D[j]) * z[k][j];
				y[k][i] = m[i] + sigma * v;
			}
			memset(&st, 0, sizeof(st));
			for (i = 0; i < N; i++)
				st.pose[i] = pose0[i] + scale[i] * y[k][i];
			rg_seed(pb, &st);
			/* fitness: energy after a short rigid refinement (basin hopping) */
			f[k] = rg_refine(pb, &st, RG_LOCAL, 0, 0.05);
			for (i = 0; i < N; i++)
				trace->cand[gen][k][i] = st.pose[i];
			if (f[k] < best_f)
			{
				best_f = f[k];
				memcpy(best_pose, st.pose, sizeof(double) * N);
			}
			order[k] = k;
		}
		trace->best[gen] = best_f;
		for (i = 1; i < RG_LAMBDA; i++)
		{
			int key = order[i];
			j = i - 1;
			while (j >= 0 && f[order[j]] > f[key]) { order[j + 1] = order[j]; j--; }
			order[j + 1] = key;
		}
		memcpy(mold, m, sizeof(m));
		for (i = 0; i < N; i++)
		{
			m[i] = 0;
			for (k = 0; k < MU; k++)
				m[i] += w[k] * y[order[k]][i];
		}
		{
			/* evolution paths */
			double zm[N] = { 0 }, bz[N], hs, norm = 0;
			for (k = 0; k < MU; k++)
				for (i = 0; i < N; i++)
					zm[i] += w[k] * z[order[k]][i];
			for (i = 0; i < N; i++)
			{
				bz[i] = 0;
				for (j = 0; j < N; j++)
					bz[i] += B[i][j] * zm[j];
			}
			for (i = 0; i < N; i++)
			{
				ps[i] = (1 - cs) * ps[i] + sqrt(cs * (2 - cs) * mueff) * bz[i];
				norm += ps[i] * ps[i];
			}
			norm = sqrt(norm);
			hs = norm / sqrt(1 - pow(1 - cs, 2.0 * (gen + 1))) / chi < 1.4 + 2.0 / (N + 1) ? 1 : 0;
			for (i = 0; i < N; i++)
				pc[i] = (1 - cc) * pc[i] + hs * sqrt(cc * (2 - cc) * mueff) * (m[i] - mold[i]) / sigma;
			for (i = 0; i < N; i++)
				for (j = 0; j < N; j++)
				{
					double rmu = 0;
					for (k = 0; k < MU; k++)
						rmu += w[k] * (y[order[k]][i] - mold[i]) * (y[order[k]][j] - mold[j]) / (sigma * sigma);
					C[i][j] = (1 - c1 - cmu) * C[i][j] + c1 * (pc[i] * pc[j] + (1 - hs) * cc * (2 - cc) * C[i][j]) + cmu * rmu;
				}
			sigma *= exp((cs / damps) * (norm / chi - 1));
		}
	}
}

static void rg_make_problem(reg_problem* pb, double theta, unsigned int seed)
{
	unsigned int rng = seed;
	int i, k = 0;
	for (i = 0; i < RG_N; i++)
	{
		double phi = 2 * PI * i / RG_N, r = rg_template_radius(phi);
		pb->tmpl[2 * i] = (qaws_scalar)(r * cos(phi));
		pb->tmpl[2 * i + 1] = (qaws_scalar)(r * sin(phi));
	}
	/* scan: rotated, scaled, bent, noisy, with a missing arc */
	for (i = 0; i < RG_SCAN; i++)
	{
		double phi = 2 * PI * i / RG_SCAN, r = rg_template_radius(phi), x, y, bx, by;
		if (phi > 4.3 && phi < 5.3)
			continue;
		x = r * cos(phi);
		y = r * sin(phi);
		bx = x + 0.12 * sin(1.3 * y);
		by = y + 0.08 * cos(1.1 * x);
		pb->scan[2 * k] = 1.15 * (cos(theta) * bx - sin(theta) * by) + 3.0 + 0.012 * rg_gauss(&rng);
		pb->scan[2 * k + 1] = 1.15 * (sin(theta) * bx + cos(theta) * by) + 0.4 + 0.012 * rg_gauss(&rng);
		k++;
	}
	pb->nscan = k;
	pb->cx = pb->cy = 0;
	for (i = 0; i < k; i++)
	{
		pb->cx += pb->scan[2 * i] / k;
		pb->cy += pb->scan[2 * i + 1] / k;
	}
}

static void rg_curve_svg(svg* s, viewport const* v, reg_problem const* pb, reg_state const* st, char const* color, double width,
	double opacity, int dashed)
{
	qaws_scalar pts[RG_N * 2];
	qaws_curve* c;
	double xy[2 * 200];
	int k;
	rg_points(pb, st, pts);
	c = rg_curve(pts);
	for (k = 0; k < 200; k++)
	{
		qaws_eval_result_2d e;
		qaws_curve_evaluate_2d(c, (qaws_scalar)(RG_N * k / 199.0), QAWS_EVAL_FLAG_POSITION, &e);
		xy[2 * k] = vx(v, e.position.x);
		xy[2 * k + 1] = vy(v, e.position.y);
	}
	svg_polyline(s, xy, 200, color, width, opacity, dashed);
	qaws_curve_destroy(c);
}

static void app_registration(void)
{
	reg_problem pb;
	reg_state gd, cm, init;
	cma_trace trace;
	double scale[4] = { 2.0, 0.3, 0.6, 0.6 }, pose0[4], best_pose[4], e_gd, e_cm;
	int success_gd = 0, success_cm = 0, t;
	double trials_theta[RG_TRIALS], trials_gd[RG_TRIALS], trials_cm[RG_TRIALS];
	svg s;
	char buf[256];

	/* the featured case: a large rotation */
	rg_make_problem(&pb, 2.6, 7u);
	memset(&init, 0, sizeof(init));
	init.pose[2] = pb.cx;
	init.pose[3] = pb.cy;
	rg_seed(&pb, &init);
	gd = init;
	e_gd = rg_refine(&pb, &gd, RG_FINAL, 1, 0.03);
	memcpy(pose0, init.pose, sizeof(pose0));
	rg_cmaes(&pb, pose0, scale, 99u, best_pose, &trace);
	memset(&cm, 0, sizeof(cm));
	memcpy(cm.pose, best_pose, sizeof(best_pose));
	rg_seed(&pb, &cm);
	e_cm = rg_refine(&pb, &cm, RG_FINAL, 1, 0.03);

	/* statistics over random rotations */
	for (t = 0; t < RG_TRIALS; t++)
	{
		reg_problem q;
		reg_state a, b;
		double th = 2 * PI * (t + 0.5) / RG_TRIALS - PI, bp[4];
		cma_trace tr;
		rg_make_problem(&q, th, 100u + (unsigned int)t);
		memset(&a, 0, sizeof(a));
		a.pose[2] = q.cx;
		a.pose[3] = q.cy;
		rg_seed(&q, &a);
		b = a;
		trials_theta[t] = th;
		trials_gd[t] = sqrt(rg_refine(&q, &a, RG_FINAL, 1, 0.03));
		rg_cmaes(&q, b.pose, scale, 500u + (unsigned int)t, bp, &tr);
		memcpy(b.pose, bp, sizeof(bp));
		rg_seed(&q, &b);
		trials_cm[t] = sqrt(rg_refine(&q, &b, RG_FINAL, 1, 0.03));
		success_gd += trials_gd[t] < 0.03;
		success_cm += trials_cm[t] < 0.03;
	}

	svg_open(&s, "showcase/app4_registration.svg", 1240, 760, "Non-rigid registration: CMA-ES basin hopping + adjoints",
		"A template (closed centripetal Catmull-Rom) is registered to a partial, noisy, bent and rotated scan. Energy: "
		"squared distances at tracked foot points (exact batch adjoints) + offset smoothness.");
	{
		viewport a = { 20, 90, 390, 390, 0.6, 5.4, -2.0, 2.8 }, b = a;
		viewport lv = { 840, 90, 380, 220, 0, 0, 0, 0 }, sv = { 840, 330, 380, 230, 0, 0, 0, 0 };
		int k, gen;
		b.x0 = 425;
		svg_panel(&s, &a, "gradient descent from the identity pose");
		svg_panel(&s, &b, "CMA-ES on the pose, gradient inside the fitness");
		for (k = 0; k < pb.nscan; k++)
		{
			svg_circle(&s, vx(&a, pb.scan[2 * k]), vy(&a, pb.scan[2 * k + 1]), 2.2, "#24292f", "#24292f");
			svg_circle(&s, vx(&b, pb.scan[2 * k]), vy(&b, pb.scan[2 * k + 1]), 2.2, "#24292f", "#24292f");
		}
		rg_curve_svg(&s, &a, &pb, &init, "#8c959f", 1.6, 1, 1);
		rg_curve_svg(&s, &a, &pb, &gd, "#cf222e", 2.6, 1, 0);
		for (gen = 0; gen < RG_GEN; gen += 3)
			for (k = 0; k < RG_LAMBDA; k++)
			{
				reg_state st;
				memset(&st, 0, sizeof(st));
				memcpy(st.pose, trace.cand[gen][k], sizeof(st.pose));
				rg_curve_svg(&s, &b, &pb, &st, "#54aeff", 0.8, 0.35, 0);
			}
		rg_curve_svg(&s, &b, &pb, &cm, "#0969da", 2.8, 1, 0);
		sprintf(buf, "RMS %.4f (wrong basin)", sqrt(e_gd));
		svg_text(&s, a.x0 + a.w - 10, a.y0 + a.h - 10, 13, "#cf222e", "end", buf);
		sprintf(buf, "RMS %.4f", sqrt(e_cm));
		svg_text(&s, b.x0 + b.w - 10, b.y0 + b.h - 10, 13, "#0969da", "end", buf);
		svg_text(&s, a.x0 + 10, a.y0 + a.h - 10, 12, "#57606a", "start", "dots: scan, dashed: template");
		svg_text(&s, b.x0 + 10, b.y0 + a.h - 10, 12, "#57606a", "start", "light: candidates every 3 generations");

		{
			double lo = 1e300, hi = -1e300, xy[2 * RG_GEN];
			for (gen = 0; gen < RG_GEN; gen++)
			{
				double l = log10(sqrt(trace.best[gen]));
				if (l < lo) lo = l;
				if (l > hi) hi = l;
			}
			if (hi - lo < 0.1) hi = lo + 0.1;
			svg_panel(&s, &lv, "best refined RMS per CMA-ES generation (log)");
			for (gen = 0; gen < RG_GEN; gen++)
			{
				xy[2 * gen] = lv.x0 + 14 + (lv.w - 28) * gen / (double)(RG_GEN - 1);
				xy[2 * gen + 1] = lv.y0 + 34 + (lv.h - 54) * (hi - log10(sqrt(trace.best[gen]))) / (hi - lo);
			}
			svg_polyline(&s, xy, RG_GEN, "#0969da", 2.4, 1, 0);
		}
		svg_panel(&s, &sv, "final RMS over 16 scan rotations (-180 to 180 deg)");
		{
			double bw = (sv.w - 40) / RG_TRIALS;
			for (t = 0; t < RG_TRIALS; t++)
			{
				double x = sv.x0 + 20 + t * bw, base = sv.y0 + sv.h - 30, hs = sv.h - 70;
				double h0 = hs * (trials_gd[t] > 0.25 ? 1 : trials_gd[t] / 0.25), h1 = hs * (trials_cm[t] > 0.25 ? 1 : trials_cm[t] / 0.25);
				fprintf(s.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"#cf222e\" opacity=\"0.8\"/>\n", x + 1, base - h0, bw * 0.42, h0 + 0.5);
				fprintf(s.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"#0969da\"/>\n", x + bw * 0.46, base - h1, bw * 0.42, h1 + 0.5);
			}
			svg_line(&s, sv.x0 + 20, sv.y0 + sv.h - 30 - (sv.h - 70) * 0.03 / 0.25, sv.x0 + sv.w - 20, sv.y0 + sv.h - 30 - (sv.h - 70) * 0.03 / 0.25, "#57606a", 1, 0.8);
			svg_text(&s, sv.x0 + 20, sv.y0 + sv.h - 12, 11, "#57606a", "start", "red: gradient descent, blue: CMA-ES + gradient; line: success threshold 0.03");
		}
		{
			double tx = 840, ty = 612;
			sprintf(buf, "registered (RMS &lt; 0.03): gradient descent %d / %d", success_gd, RG_TRIALS);
			svg_text(&s, tx, ty - 22, 13, "#cf222e", "start", buf);
			sprintf(buf, "registered (RMS &lt; 0.03): CMA-ES + gradient %d / %d", success_cm, RG_TRIALS);
			svg_text(&s, tx, ty, 13, "#24292f", "start", buf);
			svg_text(&s, tx, ty + 24, 12, "#57606a", "start", "fitness: pose energy after 10 rigid Adam steps (Baldwinian)");
			svg_text(&s, tx, ty + 42, 12, "#57606a", "start", "final: pose + 32 offsets, 250 Adam steps");
			svg_text(&s, tx, ty + 60, 12, "#57606a", "start", "gradient: 2 (C - q) -> Catmull-Rom adjoint -> points -> pose, offsets");
		}
	}
	svg_close(&s);
	printf("4_registration: featured RMS gd %.4f cma %.4f; success gd %d/%d cma %d/%d\n", sqrt(e_gd), sqrt(e_cm),
		success_gd, RG_TRIALS, success_cm, RG_TRIALS);
	(void)trials_theta;
}

/* ================================================================== */
/*  5. Hair grooming on a head: 3D strands from scalp roots           */
/* ================================================================== */

#define HG_STRANDS 1200
#define HG_CP 7
#define HG_SAMPLES 24
#define HG_ITERS 160

static double const g_head[3] = { 1.0, 1.15, 1.2 };   /* ellipsoid semi-axes */

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
	double w_len = 20, w_bend = 0.02, w_grav = 0.6, w_col = 400, w_align = 1.0;
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
		s->rest_len = 0.75 + 0.35 * ((rng >> 8) / 16777216.0);
		/* initial strand: straight out along the normal ("hedgehog") */
		for (k = 0; k < HG_CP; k++)
			for (a = 0; a < 3; a++)
				s->cps[3 * k + a] = (qaws_scalar)(s->root[a] + s->normal[a] * s->rest_len * k / (HG_CP - 1));
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
	double e = 0, w_img = 3.0, w_mask = 2.0, ang = 0;
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
		/* stay on the hair: (1 - m)^2 */
		m = sample(hp->mask, W, H, px, py, &mdx, &mdy);
		e += w_mask * (1 - m) * (1 - m) / M;
		dpos[0] = -w_mask * 2 * (1 - m) * mdx / M;
		dpos[1] = -w_mask * 2 * (1 - m) * mdy / M;
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

int main(int argc, char** argv)
{
	setvbuf(stdout, NULL, _IONBF, 0);
	if (argc > 1)
		g_photos = argv[1];
	MAKE_DIR("showcase");
	{
		/* QAWS_APP=n runs only application n */
		char const* only = getenv("QAWS_APP");
		int pick = only ? atoi(only) : 0;
		if (!pick || pick == 1) app_hair();
		if (!pick || pick == 2) app_vectorize();
		if (!pick || pick == 3) app_mesh_patches(argc > 2 ? argv[2] : NULL);
		if (!pick || pick == 4) app_registration();
		if (!pick || pick == 5) app_hair_groom();
		if (!pick || pick == 6) app_hair_photo();
	}
	return 0;
}
