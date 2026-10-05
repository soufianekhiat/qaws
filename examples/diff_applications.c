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
 *   7. Rational patches: NURBS weights freed after the ADMM fit (blob and
 *      sphere meshes), seams kept closed.
 *
 *   8. Real haircuts: vector hair strands from six portrait photos with a
 *      seeded hair mask (photos/hair_*.ppm, see photos/CREDITS.txt).
 *
 *   9. 3D hair from four frontal portraits: head placed from the face,
 *      visible strands follow the photo, hidden ones the 3D priors.
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

static double g_blob_amp = 1.0;   /* 0: exact sphere, 1: bumpy blob */

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
				r = 1.0 + g_blob_amp * (blob_radius(p[0], p[1], p[2]) - 1.0);
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

/* --- rational refinement: weights become free ----------------------- */

static qaws_surface* mp_surface_nurbs(qaws_scalar const* cps, qaws_scalar const* w)
{
	qaws_surface_nurbs_desc d;
	qaws_surface* s = NULL;
	memset(&d, 0, sizeof(d));
	d.u_degree = 3;
	d.v_degree = 3;
	d.control_points = (qaws_vec3 const*)cps;
	d.u_point_count = MP_N;
	d.v_point_count = MP_N;
	d.weights = w;
	d.u_knots = g_mp_knots;
	d.u_knot_count = MP_N + 4;
	d.v_knots = g_mp_knots;
	d.v_knot_count = MP_N + 4;
	qaws_surface_create_nurbs(&d, &s);
	return s;
}

/* Mean squared distance at the foot points of a NURBS patch with its
   gradient on control points and weights (exact rational adjoints). */
static double mp_rational_energy(patch_fit const* p, mesh const* m, qaws_scalar const* cps, qaws_scalar const* w,
	qaws_scalar* g_cp, qaws_scalar* g_w)
{
	qaws_surface* s = mp_surface_nurbs(cps, w);
	qaws_surface_jet* bars = (qaws_surface_jet*)malloc(sizeof(qaws_surface_jet) * (size_t)(p->n + 1));
	qaws_scalar* us = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * (size_t)(p->n + 1));
	qaws_scalar* vs = us + p->n;
	qaws_field_view fv[2];
	qaws_diff_views views;
	double e = 0;
	int i;
	for (i = 0; i < p->n; i++)
	{
		double const* x = &m->v[3 * p->pts[i]];
		qaws_surface_jet j;
		double rx, ry, rz;
		us[i] = p->uv[2 * i];
		vs[i] = p->uv[2 * i + 1];
		qaws_surface_eval_jet(s, us[i], vs[i], QAWS_SJET_P, &j);
		rx = j.d[0].x - x[0]; ry = j.d[0].y - x[1]; rz = j.d[0].z - x[2];
		e += (rx * rx + ry * ry + rz * rz) / p->n;
		memset(&bars[i], 0, sizeof(bars[i]));
		bars[i].d[0] = v3((qaws_scalar)(2 * rx / p->n), (qaws_scalar)(2 * ry / p->n), (qaws_scalar)(2 * rz / p->n));
		bars[i].channels = QAWS_SJET_P;
	}
	memset(g_cp, 0, sizeof(qaws_scalar) * MP_CP * 3);
	memset(g_w, 0, sizeof(qaws_scalar) * MP_CP);
	fv[0] = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, g_cp, MP_CP, 3);
	fv[1] = qaws_field_view_make(QAWS_FIELD_WEIGHTS, g_w, MP_CP, 1);
	views.fields = fv;
	views.field_count = 2;
	views.children = NULL;
	views.child_count = 0;
	if (p->n)
		qaws_surface_eval_batch_adjoint(NULL, s, us, vs, (unsigned int)p->n, QAWS_SJET_P, bars, &views, NULL, NULL);
	free(bars);
	free(us);
	qaws_surface_destroy(s);
	return e;
}

/* Joint Adam on control points and log-weights of all patches; after each
   step the copies of every seam node (positions and weights) are replaced
   by their mean, so the patch set stays closed. Foot points are
   re-projected every 10 steps. Returns the final RMS. */
static double mp_refine_rational(mesh const* m, patch_fit* pf, int nodes, qaws_scalar (*weights)[MP_CP], int iters, int free_weights,
	double* rms_hist)
{
	static adam opt[MP_PATCHES];
	qaws_scalar* zsum = (qaws_scalar*)calloc((size_t)nodes * 4, sizeof(qaws_scalar));
	int* cnt = (int*)calloc((size_t)nodes, sizeof(int));
	int it, f, k, c;
	double rms = 0;
	for (f = 0; f < MP_PATCHES; f++)
	{
		memset(&opt[f], 0, sizeof(adam));
		for (k = 0; k < MP_CP; k++)
		{
			weights[f][k] = 1;
			if (pf[f].node[k] >= 0)
				cnt[pf[f].node[k]]++;
		}
	}
	for (it = 0; it <= iters; it++)
	{
		double e = 0;
		int npts = 0;
		if (it % 10 == 0 || it == iters)
			for (f = 0; f < MP_PATCHES; f++)
			{
				qaws_surface* s = mp_surface_nurbs(pf[f].cps, weights[f]);
				double ef = mp_project(&pf[f], m, s, 3);
				qaws_surface_destroy(s);
				e += ef * pf[f].n;
				npts += pf[f].n;
			}
		if (it % 10 == 0 || it == iters)
		{
			rms = sqrt(e / npts);
			if (rms_hist)
				rms_hist[it / 10] = rms;
		}
		if (it == iters)
			break;
		for (f = 0; f < MP_PATCHES; f++)
		{
			qaws_scalar g_cp[MP_CP * 3], g_w[MP_CP], x[MP_CP * 4], g[MP_CP * 4];
			mp_rational_energy(&pf[f], m, pf[f].cps, weights[f], g_cp, g_w);
			/* parameters: positions and log-weights (d/dlog w = w d/dw) */
			for (k = 0; k < MP_CP * 3; k++) { x[k] = pf[f].cps[k]; g[k] = g_cp[k]; }
			for (k = 0; k < MP_CP; k++)
			{
				x[MP_CP * 3 + k] = (qaws_scalar)log(weights[f][k]);
				g[MP_CP * 3 + k] = free_weights ? weights[f][k] * g_w[k] : 0;
			}
			adam_step(&opt[f], x, g, MP_CP * 4, 0.002 * (1.0 - 0.7 * it / (double)iters));
			for (k = 0; k < MP_CP * 3; k++) pf[f].cps[k] = x[k];
			for (k = 0; k < MP_CP; k++) weights[f][k] = (qaws_scalar)exp(x[MP_CP * 3 + k]);
		}
		/* seam consensus by projection */
		memset(zsum, 0, sizeof(qaws_scalar) * (size_t)nodes * 4);
		for (f = 0; f < MP_PATCHES; f++)
			for (k = 0; k < MP_CP; k++)
				if (pf[f].node[k] >= 0)
				{
					for (c = 0; c < 3; c++)
						zsum[4 * pf[f].node[k] + c] += pf[f].cps[3 * k + c] / cnt[pf[f].node[k]];
					zsum[4 * pf[f].node[k] + 3] += weights[f][k] / cnt[pf[f].node[k]];
				}
		for (f = 0; f < MP_PATCHES; f++)
			for (k = 0; k < MP_CP; k++)
				if (pf[f].node[k] >= 0)
				{
					for (c = 0; c < 3; c++)
						pf[f].cps[3 * k + c] = zsum[4 * pf[f].node[k] + c];
					weights[f][k] = zsum[4 * pf[f].node[k] + 3];
				}
	}
	free(zsum);
	free(cnt);
	return rms;
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

/* ================================================================== */
/*  7. Rational patches: freeing the NURBS weights after ADMM         */
/* ================================================================== */

#define RP_ITERS 300

typedef struct rational_run
{
	double rms_admm, rms_ctl, rms_rat, gap_rat, wmin, wmax;
	double hist_ctl[RP_ITERS / 10 + 1], hist_rat[RP_ITERS / 10 + 1];
} rational_run;

static void rational_experiment(double amp, rational_run* out, patch_fit* rat_out, qaws_scalar (*w_out)[MP_CP], mesh* m_out)
{
	static patch_fit admm[MP_PATCHES], ctl[MP_PATCHES];
	static qaws_scalar wc[MP_PATCHES][MP_CP];
	double rms[MP_OUTER], gap[MP_OUTER];
	int nodes, f, k;
	g_blob_amp = amp;
	mesh_blob(m_out);
	mp_setup(m_out, admm, &nodes);
	mp_solve(m_out, admm, nodes, 0.03, rms, gap);
	out->rms_admm = rms[MP_OUTER - 1];
	for (f = 0; f < MP_PATCHES; f++)
	{
		size_t n = (size_t)admm[f].n;
		ctl[f] = admm[f];
		rat_out[f] = admm[f];
		/* own copies of the foot points */
		ctl[f].uv = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * n);
		rat_out[f].uv = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * n);
		memcpy(ctl[f].uv, admm[f].uv, sizeof(qaws_scalar) * 2 * n);
		memcpy(rat_out[f].uv, admm[f].uv, sizeof(qaws_scalar) * 2 * n);
	}
	out->rms_ctl = mp_refine_rational(m_out, ctl, nodes, wc, RP_ITERS, 0, out->hist_ctl);
	out->rms_rat = mp_refine_rational(m_out, rat_out, nodes, w_out, RP_ITERS, 1, out->hist_rat);
	out->gap_rat = mp_seam_gap(rat_out, nodes);
	out->wmin = 1e30;
	out->wmax = 0;
	for (f = 0; f < MP_PATCHES; f++)
		for (k = 0; k < MP_CP; k++)
		{
			if (w_out[f][k] < out->wmin) out->wmin = w_out[f][k];
			if (w_out[f][k] > out->wmax) out->wmax = w_out[f][k];
		}
	for (f = 0; f < MP_PATCHES; f++)
	{
		free(admm[f].pts); free(admm[f].uv);
		free(ctl[f].uv);
	}
}

static void app_rational_patches(void)
{
	static patch_fit rat_blob[MP_PATCHES], rat_sphere[MP_PATCHES];
	static qaws_scalar w_blob[MP_PATCHES][MP_CP], w_sphere[MP_PATCHES][MP_CP];
	rational_run rb, rs;
	mesh mb, ms;
	svg s;
	char buf[256];
	int q;
	rational_experiment(1.0, &rb, rat_blob, w_blob, &mb);
	rational_experiment(0.0, &rs, rat_sphere, w_sphere, &ms);

	svg_open(&s, "showcase/app7_rational_patches.svg", 1240, 700, "Rational patches: freeing the NURBS weights",
		"After ADMM: refine control points only (weights 1) or control points + log-weights through the rational "
		"surface adjoints; seam copies stay averaged, so the patch set stays closed.");
	for (q = 0; q < 2; q++)
	{
		rational_run const* r = q ? &rs : &rb;
		qaws_scalar (*w)[MP_CP] = q ? w_sphere : w_blob;
		viewport lv = { 20 + q * 610, 90, 590, 250, 0, 0, 0, 0 };
		double lo = 1e300, hi = -1e300, x0[2 * (RP_ITERS / 10 + 1)], x1[2 * (RP_ITERS / 10 + 1)];
		int n = RP_ITERS / 10 + 1, k, f, i, j;
		for (k = 0; k < n; k++)
		{
			double a = log10(r->hist_ctl[k]), b = log10(r->hist_rat[k]);
			if (a < lo) lo = a; if (b < lo) lo = b;
			if (a > hi) hi = a; if (b > hi) hi = b;
		}
		sprintf(buf, "%s: log10 RMS distance per step", q ? "sphere mesh" : "blob mesh");
		svg_panel(&s, &lv, buf);
		for (k = 0; k < n; k++)
		{
			x0[2 * k] = x1[2 * k] = lv.x0 + 14 + (lv.w - 28) * k / (double)(n - 1);
			x0[2 * k + 1] = lv.y0 + 34 + (lv.h - 54) * (hi - log10(r->hist_ctl[k])) / (hi - lo + 1e-12);
			x1[2 * k + 1] = lv.y0 + 34 + (lv.h - 54) * (hi - log10(r->hist_rat[k])) / (hi - lo + 1e-12);
		}
		svg_polyline(&s, x0, n, "#8c959f", 2.2, 1, 0);
		svg_polyline(&s, x1, n, "#0969da", 2.6, 1, 0);
		svg_text(&s, lv.x0 + lv.w - 14, x0[2 * n - 1] - 8, 11, "#57606a", "end", "control points only");
		svg_text(&s, lv.x0 + lv.w - 14, x1[2 * n - 1] + 16, 11, "#0969da", "end", "control points + weights");
		sprintf(buf, "RMS: ADMM %.2e, control points %.2e, + weights %.2e (%.0f%% lower)", r->rms_admm, r->rms_ctl, r->rms_rat,
			100.0 * (1 - r->rms_rat / r->rms_ctl));
		svg_text(&s, lv.x0, lv.y0 + lv.h + 24, 13, "#24292f", "start", buf);
		sprintf(buf, "weights %.3f .. %.3f, max seam gap %.1e", r->wmin, r->wmax, r->gap_rat);
		svg_text(&s, lv.x0, lv.y0 + lv.h + 44, 13, "#57606a", "start", buf);
		/* weight grids per patch */
		for (f = 0; f < MP_PATCHES; f++)
		{
			double gx = lv.x0 + 8 + f * 97, gy = lv.y0 + lv.h + 70, cs = 12;
			double wl = q ? 0.96 : 0.75, wh = q ? 1.04 : 1.25;
			for (i = 0; i < MP_N; i++)
				for (j = 0; j < MP_N; j++)
				{
					char col[32];
					heat((w[f][i * MP_N + j] - wl) / (wh - wl), col);
					fprintf(s.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"%s\"/>\n",
						gx + j * cs, gy + i * cs, cs - 1, cs - 1, col);
				}
			sprintf(buf, "patch %d", f + 1);
			svg_text(&s, gx, gy + MP_N * cs + 14, 11, "#57606a", "start", buf);
		}
		sprintf(buf, "%.2f", q ? 0.96 : 0.75);
		{
			char hb[16];
			sprintf(hb, "%.2f", q ? 1.04 : 1.25);
			svg_colorbar(&s, lv.x0 + 8, lv.y0 + lv.h + 175, 220, 8, buf, hb);
		}
	}
	svg_text(&s, 20, 660, 13, "#57606a", "start",
		"gradient: 2 (S - X) at the foot points -> qaws_surface_eval_batch_adjoint with CONTROL_POINTS and WEIGHTS views "
		"(rational quotient rule), weights parameterized by their log to stay positive");
	svg_close(&s);
	printf("7_rational: blob admm %.5f ctl %.5f rat %.5f (w %.3f..%.3f, gap %.1e); sphere admm %.5f ctl %.5f rat %.5f (w %.3f..%.3f)\n",
		rb.rms_admm, rb.rms_ctl, rb.rms_rat, rb.wmin, rb.wmax, rb.gap_rat, rs.rms_admm, rs.rms_ctl, rs.rms_rat, rs.wmin, rs.wmax);
}

/* ================================================================== */
/*  8. Real haircuts: vector hair strands from portrait photos        */
/* ================================================================== */

typedef struct hair_spec
{
	char const* name;
	char const* label;
	int hair[8][2], nhair;       /* hair color seeds (pixels of the 360-wide photo) */
	int non[8][2], nnon;         /* skin / clothes / background seeds */
	int box[4];                  /* x0, y0, x1, y1: where hair may be */
	int face[4];                 /* face ellipse: cx, cy, rx, ry (never hair) */
} hair_spec;

static hair_spec const g_hair_specs[6] = {
	{ "hair_straight", "straight, shoulder length",
		{ { 120, 90 }, { 100, 180 }, { 240, 200 }, { 170, 60 }, { 110, 240 }, { 250, 240 } }, 6,
		{ { 175, 150 }, { 175, 200 }, { 40, 80 }, { 330, 100 }, { 180, 330 }, { 60, 300 }, { 300, 300 } }, 7,
		{ 60, 20, 300, 252 }, { 175, 150, 42, 55 } },
	{ "hair_short_curly", "short, curly",
		{ { 120, 40 }, { 80, 90 }, { 150, 30 }, { 80, 140 }, { 190, 40 } }, 5,
		{ { 130, 110 }, { 150, 150 }, { 60, 250 }, { 200, 200 }, { 330, 60 }, { 20, 100 }, { 260, 120 } }, 7,
		{ 40, 0, 240, 175 }, { 135, 118, 38, 52 } },
	{ "hair_coily", "voluminous, coily",
		{ { 100, 120 }, { 300, 120 }, { 200, 70 }, { 80, 190 }, { 300, 190 }, { 260, 60 }, { 130, 60 } }, 7,
		{ { 190, 150 }, { 190, 200 }, { 160, 170 }, { 220, 170 }, { 190, 115 }, { 40, 60 }, { 340, 250 }, { 200, 330 } }, 8,
		{ 50, 30, 345, 265 }, { 190, 165, 48, 62 } },
	{ "hair_long", "long, straight (profile)",
		{ { 260, 80 }, { 300, 250 }, { 320, 400 }, { 250, 160 }, { 200, 40 }, { 240, 300 } }, 6,
		{ { 200, 180 }, { 100, 40 }, { 60, 140 }, { 100, 330 }, { 220, 380 }, { 40, 420 } }, 6,
		{ 150, 0, 360, 456 }, { 195, 185, 45, 72 } },
	{ "hair_twists", "long twists",
		{ { 140, 140 }, { 110, 250 }, { 280, 180 }, { 270, 260 }, { 210, 80 }, { 150, 90 }, { 100, 330 } }, 7,
		{ { 210, 180 }, { 220, 220 }, { 180, 190 }, { 240, 190 }, { 210, 135 }, { 60, 100 }, { 330, 150 }, { 240, 400 } }, 8,
		{ 60, 60, 330, 440 }, { 213, 185, 45, 58 } },
	{ "hair_wavy", "long, wavy, grey (with beard)",
		{ { 100, 60 }, { 200, 30 }, { 60, 250 }, { 320, 240 }, { 300, 120 }, { 50, 150 }, { 180, 330 } }, 7,
		{ { 180, 180 }, { 150, 170 }, { 80, 400 }, { 340, 420 }, { 350, 30 } }, 5,
		{ 15, 0, 360, 420 }, { 185, 175, 52, 58 } }
};

/* Seeded hair mask: patch-averaged color, nearest hair seed vs nearest
   other seed, restricted to the box, kept connected to the hair seeds,
   then softened. */
static float* hair_mask(image const* img, hair_spec const* sp)
{
	int w = img->w, h = img->h, n = w * h, i, x, y, k;
	float* avg = (float*)malloc(sizeof(float) * 3 * (size_t)n);
	float* m = (float*)calloc((size_t)n, sizeof(float));
	unsigned char* lab = (unsigned char*)calloc((size_t)n, 1);
	int* stack = (int*)malloc(sizeof(int) * (size_t)n);
	double hc[8][3], nc[8][3];
	int top = 0;
	for (k = 0; k < 3; k++)
	{
		float* ch = (float*)malloc(sizeof(float) * (size_t)n);
		for (i = 0; i < n; i++)
			ch[i] = img->rgb[3 * i + k];
		blur(ch, w, h, 2.0);
		for (i = 0; i < n; i++)
			avg[3 * i + k] = ch[i];
		free(ch);
	}
	for (k = 0; k < sp->nhair; k++)
		for (i = 0; i < 3; i++)
			hc[k][i] = avg[3 * (sp->hair[k][1] * w + sp->hair[k][0]) + i];
	for (k = 0; k < sp->nnon; k++)
		for (i = 0; i < 3; i++)
			nc[k][i] = avg[3 * (sp->non[k][1] * w + sp->non[k][0]) + i];
	for (y = 0; y < h; y++)
		for (x = 0; x < w; x++)
		{
			double dh = 1e9, dn = 1e9;
			float const* c = &avg[3 * (y * w + x)];
			if (x < sp->box[0] || x >= sp->box[2] || y < sp->box[1] || y >= sp->box[3])
				continue;
			{
				double ex = (x - sp->face[0]) / (double)sp->face[2], ey = (y - sp->face[1]) / (double)sp->face[3];
				if (ex * ex + ey * ey < 1)
					continue;
			}
			for (k = 0; k < sp->nhair; k++)
			{
				double d = (c[0] - hc[k][0]) * (c[0] - hc[k][0]) + (c[1] - hc[k][1]) * (c[1] - hc[k][1]) + (c[2] - hc[k][2]) * (c[2] - hc[k][2]);
				if (d < dh) dh = d;
			}
			for (k = 0; k < sp->nnon; k++)
			{
				double d = (c[0] - nc[k][0]) * (c[0] - nc[k][0]) + (c[1] - nc[k][1]) * (c[1] - nc[k][1]) + (c[2] - nc[k][2]) * (c[2] - nc[k][2]);
				if (d < dn) dn = d;
			}
			lab[y * w + x] = dh < dn ? 1 : 0;
		}
	/* connected to the hair seeds */
	for (k = 0; k < sp->nhair; k++)
	{
		int id = sp->hair[k][1] * w + sp->hair[k][0];
		if (lab[id] == 1)
		{
			lab[id] = 2;
			stack[top++] = id;
		}
	}
	while (top > 0)
	{
		int id = stack[--top], px = id % w, py = id / w, d;
		static int const nb[4][2] = { { 1, 0 }, { -1, 0 }, { 0, 1 }, { 0, -1 } };
		m[id] = 1;
		for (d = 0; d < 4; d++)
		{
			int qx = px + nb[d][0], qy = py + nb[d][1], q;
			if (qx < 0 || qy < 0 || qx >= w || qy >= h)
				continue;
			q = qy * w + qx;
			if (lab[q] == 1)
			{
				lab[q] = 2;
				stack[top++] = q;
			}
		}
	}
	blur(m, w, h, 1.5);
	free(avg);
	free(lab);
	free(stack);
	return m;
}

/* Fits and optimizes strands (as in application 1); accumulates stats. */
static void hair_fit_strands(tensor_field const* tf, strand* strands, int ns, double* a0, double* a1, int* points, unsigned int* cps_total)
{
	int i;
	*a0 = *a1 = 0;
	*points = 0;
	*cps_total = 0;
	for (i = 0; i < ns; i++)
	{
		strand* st = &strands[i];
		unsigned int ncp = strand_cp_count(st), got = 0;
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
		st->fit = st->opt = NULL;
		if (qaws_curve_fit_bspline(&d, &st->fit) != QAWS_STATUS_OK)
		{
			st->fit = NULL;
			continue;
		}
		qaws_curve_read_field(st->fit, QAWS_FIELD_CONTROL_POINTS, cps, ncp * 2, &got);
		qaws_curve_read_field(st->fit, QAWS_FIELD_KNOTS, st->knots, ncp + 4, &got);
		qaws_curve_destroy(st->fit);
		st->fit = strand_curve(cps, ncp, st->knots);
		strand_alignment(tf, st->fit, &ea, &aa);
		*a0 += aa;
		memset(&opt, 0, sizeof(opt));
		for (it = 0; it < HS_ITERS; it++)
		{
			memset(g, 0, sizeof(g));
			strand_energy(tf, st, cps, ncp, g);
			adam_step(&opt, cps, g, (int)(2 * ncp), 0.25 * (1.0 - 0.8 * it / (double)HS_ITERS));
		}
		st->opt = strand_curve(cps, ncp, st->knots);
		strand_alignment(tf, st->opt, &ea, &aa);
		*a1 += aa;
		*points += st->n;
		*cps_total += ncp;
	}
	if (ns)
	{
		*a0 /= ns;
		*a1 /= ns;
	}
}

/* Strand drawn in segments colored by direction (hue = angle mod 180). */
static void strand_svg_hue(svg* s, imgmap const* m, qaws_curve const* c, double width)
{
	qaws_range r = qaws_curve_get_parameter_range(c);
	double px = 0, py = 0;
	int k, n = 40;
	for (k = 0; k <= n; k++)
	{
		qaws_eval_result_2d e;
		double sx, sy;
		qaws_curve_evaluate_2d(c, (qaws_scalar)(r.min_value + (r.max_value - r.min_value) * k / n),
			QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &e);
		map_pt(m, e.position.x, e.position.y, &sx, &sy);
		if (k > 0)
		{
			char col[48];
			double ang = atan2(e.d1.y, e.d1.x) * 180.0 / PI;
			if (ang < 0) ang += 180;
			sprintf(col, "hsl(%d,85%%,60%%)", (int)(2 * ang) % 360);
			svg_line(s, px, py, sx, sy, col, width, 0.95);
		}
		px = sx;
		py = sy;
	}
}

static void app_haircuts(void)
{
	static strand strands[HS_MAX_STRANDS];
	svg s;
	char buf[256], path[512];
	int p;

	svg_open(&s, "showcase/app8_haircuts.svg", 1240, 795, "Real haircuts: vector hair from portrait photos",
		"Seeded hair mask, orientation field, evenly spaced streamlines in the mask, B-spline strands optimized to "
		"follow the hair through unit-tangent adjoints (right: strands colored by direction).");
	for (p = 0; p < 6; p++)
	{
		hair_spec const* sp = &g_hair_specs[p];
		image img, masked;
		tensor_field tf;
		float* mask;
		int ns, i, pts;
		unsigned int ncps;
		double a0, a1, sc, cw = 400, x0 = 20 + (p % 3) * 405, y0 = 90 + (p / 3) * 342;
		imgmap mp, mv;

		sprintf(path, "%s/%s.ppm", g_photos, sp->name);
		if (!image_load_ppm(path, &img))
		{
			printf("8_haircuts: %s missing\n", path);
			continue;
		}
		mask = hair_mask(&img, sp);
		tensor_build(&img, 1.0, 2.5, &tf);
		/* tracing reads the mask through the luminance channel */
		masked = img;
		masked.lum = mask;
		ns = trace_strands(&masked, &tf, 2.0, 0.2, 0.5, strands, HS_MAX_STRANDS);
		hair_fit_strands(&tf, strands, ns, &a0, &a1, &pts, &ncps);

		sc = (cw / 2 - 4) / img.w;
		mp.x0 = x0; mp.y0 = y0 + 22; mp.ox = 0; mp.oy = 0; mp.scale = sc;
		mv = mp;
		mv.x0 = x0 + cw / 2 + 2;
		sprintf(buf, "%s", sp->label);
		svg_text(&s, x0, y0 + 14, 13, "#24292f", "start", buf);
		fprintf(s.f, "<image href=\"../%s/%s.png\" x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\"/>\n",
			g_photos, sp->name, mp.x0, mp.y0, img.w * sc, img.h * sc);
		/* mask outline: dim the non-hair pixels on the right panel */
		fprintf(s.f, "<image href=\"../%s/%s.png\" x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" opacity=\"0.22\"/>\n",
			g_photos, sp->name, mv.x0, mv.y0, img.w * sc, img.h * sc);
		fprintf(s.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"#0d1117\" opacity=\"0.85\"/>\n",
			mv.x0, mv.y0, img.w * sc, img.h * sc);
		for (i = 0; i < ns; i++)
			if (strands[i].opt)
				strand_svg_hue(&s, &mv, strands[i].opt, 1.1);
		sprintf(buf, "%d strands, %d points -> %u control points; angle %.1f -> %.1f deg", ns, pts, ncps, a0, a1);
		svg_text(&s, x0, y0 + 22 + img.h * sc + 16, 11, "#57606a", "start", buf);
		printf("8_haircuts: %-18s %4d strands, %6d points -> %4u CPs, angle %.2f -> %.2f deg\n", sp->name, ns, pts, ncps, a0, a1);

		for (i = 0; i < ns; i++)
		{
			free(strands[i].xy);
			if (strands[i].fit) qaws_curve_destroy(strands[i].fit);
			if (strands[i].opt) qaws_curve_destroy(strands[i].opt);
		}
		free(mask);
		free(tf.xx); free(tf.xy); free(tf.yy);
		image_free(&img);
	}
	svg_text(&s, 20, 782, 11, "#57606a", "start",
		"Photos: Wikimedia Commons (CC0, public domain, CC BY 2.0, CC BY-SA 2.0/4.0); authors and licenses in photos/CREDITS.txt.");
	svg_close(&s);
}

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
#define HV_GROWN_MAX 16000
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
			   in a frontal view, a curtain hanging down the back. */
			double d = dist[y * W + x] / s_px, infl, tf, tb, f, b, u, z;
			if (d > R) d = R;
			infl = sqrt(fmax(0, 2 * R * d - d * d));
			u = (x - hc->cx) / s_px;
			z = (hc->cy - y) / s_px;
			f = infl;
			b = -infl;
			if (frontal && z < 0.3)
				b = fmin(b, -0.9 * g_head[1] * sqrt(fmax(0, 1 - u * u / 2.6)));
			if (hv_head_hit(&sc->cam, x, y, &tf, &tb))
			{
				f = fmax(tf + 0.05, f);
				b = fmin(tb - 0.05, b);
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
		{ 20, 20, 345, 435 }, { 175, 162, 46, 64 } }, PI, 177, 130, 71, 0, { 0, 0, 0 } },
	{ { "plain_curly", "curly, voluminous",
		{ { 80, 150 }, { 60, 250 }, { 300, 250 }, { 310, 330 }, { 170, 40 }, { 240, 60 }, { 90, 380 }, { 280, 400 } }, 8,
		{ { 190, 200 }, { 190, 250 }, { 160, 160 }, { 220, 160 }, { 210, 460 }, { 190, 330 } }, 6,
		{ 20, 10, 358, 430 }, { 190, 200, 56, 78 } }, PI, 190, 172, 94, 1, { 0, 0, 0 } },
	{ { "plain_long", "long, straight",
		{ { 150, 40 }, { 220, 40 }, { 100, 160 }, { 75, 240 }, { 250, 200 }, { 55, 300 }, { 110, 120 }, { 240, 110 } }, 8,
		{ { 180, 150 }, { 180, 200 }, { 250, 235 }, { 200, 280 }, { 300, 280 }, { 170, 300 }, { 140, 230 }, { 330, 300 } }, 8,
		{ 20, 0, 300, 341 }, { 180, 140, 42, 58 } }, PI, 180, 107, 69, 0, { 0, 0, 0 } },
	{ { "plain_ponytail", "high ponytail, held up",
		{ { 50, 135 }, { 80, 137 }, { 110, 128 }, { 140, 112 }, { 170, 86 }, { 200, 52 }, { 222, 25 }, { 213, 70 } }, 8,
		{ { 255, 100 }, { 245, 132 }, { 290, 180 }, { 25, 150 }, { 270, 70 } }, 5,
		{ 15, 0, 300, 152 }, { 255, 100, 32, 45 } }, PI + 0.35, 245, 82, 48, 0, { 212, 22, -0.25 } },
	{ { "plain_bob_profile", "bob, profile view",
		{ { 150, 60 }, { 250, 80 }, { 300, 200 }, { 280, 300 }, { 100, 120 }, { 320, 330 }, { 200, 120 } }, 7,
		{ { 130, 250 }, { 100, 330 }, { 340, 340 }, { 250, 450 }, { 215, 230 } }, 5,
		{ 40, 15, 360, 380 }, { 125, 260, 70, 90 } }, PI / 2, 195, 185, 130, 0, { 0, 0, 0 } }
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

/* Strands grown from scalp roots through the field (RK2, half-cell steps):
   they leave the scalp along the field (signed by the combing guide, or
   toward the tie of a ponytail), stay out of the head, and stop when they
   leave the volume, turn sharply or reach the point budget. Vertex colors
   come from the photo where the strand is visible, its mean elsewhere. */
static int hv_grow(hv_scene const* sc, double const* tie, hv_volume* v, hv_grown* out, int cap)
{
	int W = sc->img.w, H = sc->img.h, cnt = 0, i, a, cand = 4 * cap;
	double h = 0.5 * HV_GH;
	for (i = 0; i < cand && cnt < cap; i++)
	{
		double z = 1 - 2 * (i + 0.5) / cand, rad = sqrt(1 - z * z), phi = i * 2.399963229728653;
		double d[3] = { rad * cos(phi), rad * sin(phi), z }, pos[3], nrm[3], g[3], dir[3], ln, gn, col[3] = { 0, 0, 0 }, cw = 0;
		hv_grown* s = &out[cnt];
		int step, outside = 0, passed = 0;
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
		s->n = 0;
		for (step = 0; step < 2 * HV_GROWN_PTS - 1; step++)
		{
			double d1[3], d2[3], mid[3], gp[3], r;
			int pass;
			for (pass = 0; pass < 2; pass++)
			{
				double const* at = pass ? mid : pos;
				double* dd = pass ? d2 : d1;
				hv_field_dir(v, at, pass ? d1 : dir, dd);
				if (tie && !passed)
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
			if (d2[0] * dir[0] + d2[1] * dir[1] + d2[2] * dir[2] < 0.6)
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
			if (tie && !passed)
			{
				double dx = pos[0] - tie[0], dy = pos[1] - tie[1], dz = pos[2] - tie[2];
				if (dx * dx + dy * dy + dz * dz < 0.15 * 0.15)
					passed = 1;
			}
			if (hv_trilinear(v->occ, 1, 0, pos) < 0.5)
			{
				if (++outside > 4)
					break;
			}
			else
				outside = 0;
			if (step % 2 == 0 && s->n < HV_GROWN_PTS)
			{
				for (a = 0; a < 3; a++)
					s->p[3 * s->n + a] = (float)pos[a];
				s->n++;
			}
		}
		/* trim the samples that left the volume */
		while (s->n > 0)
		{
			double q[3] = { s->p[3 * (s->n - 1)], s->p[3 * (s->n - 1) + 1], s->p[3 * (s->n - 1) + 2] };
			if (hv_trilinear(v->occ, 1, 0, q) >= 0.5)
				break;
			s->n--;
		}
		if (s->n < 4)
			continue;
		/* colors: the photo where the vertex is on the visible front */
		for (step = 0; step < s->n; step++)
		{
			double q[3] = { s->p[3 * step], s->p[3 * step + 1], s->p[3 * step + 2] }, px, py, t;
			int ix, iy, vis = 0;
			hv_project(&sc->cam, q, &px, &py, &t);
			ix = (int)px;
			iy = (int)py;
			if (ix >= 0 && iy >= 0 && ix < W && iy < H && sc->mask[iy * W + ix] > 0.5f &&
				t > sample(sc->surf[0], W, H, px, py, NULL, NULL) - 0.25)
				vis = 1;
			for (a = 0; a < 3; a++)
				s->rgb[3 * step + a] = vis ? sc->img.rgb[3 * (iy * W + ix) + a] : -1.0f;
			if (vis)
			{
				for (a = 0; a < 3; a++)
					col[a] += s->rgb[3 * step + a];
				cw++;
			}
		}
		for (a = 0; a < 3; a++)
			col[a] = cw > 0 ? col[a] / cw : 0.85 * sc->hair_col[a];
		for (step = 0; step < s->n; step++)
			if (s->rgb[3 * step] < 0)
				for (a = 0; a < 3; a++)
					s->rgb[3 * step + a] = (float)(0.9 * col[a]);
		for (step = 0; step < s->n; step++)
		{
			int ci, cj, ck;
			double q[3] = { s->p[3 * step], s->p[3 * step + 1], s->p[3 * step + 2] };
			if (hv_cell(q, &ci, &cj, &ck))
				v->dens[HV_IDX(ci, cj, ck)] += 1;
		}
		cnt++;
	}
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
		int nt = 0, n, ng, v;
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
		ng = hv_grow(&sc, hc->tie[0] > 0 ? tie : NULL, &vol, grown, HV_GROWN_MAX);
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
		sprintf(buf, "%s: %d image strands -> %d fitted -> %d grown strands; projected angle to the photo %.1f deg, hair covered %.0f%%",
			sp->label, nt, n, ng, angle, 100 * coverage);
		svg_text(&s, 20, y0 + ih + 16, 12, "#0969da", "start", buf);
		printf("10_hair_volume: %-18s %d image strands, %d fitted (RMS %.2f px), %d grown, angle %.1f deg, coverage %.0f%%\n",
			sp->name, nt, n, rms, ng, angle, 100 * coverage);
		hv_volume_free(&vol);
		hv_free(&sc, traced, nt);
		y0 += ih + 32;
	}
	svg_text(&s, 20, y0 + 8, 11, "#57606a", "start",
		"Columns: photo; grown strands from the photo camera; three new views. Photos: Wikimedia Commons (see photos/CREDITS.txt).");
	svg_close(&s);
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
		if (!pick || pick == 7) app_rational_patches();
		if (!pick || pick == 8) app_haircuts();
		if (!pick || pick == 9) app_hair3d();
		if (!pick || pick == 10) app_hair_volume();
	}
	return 0;
}
