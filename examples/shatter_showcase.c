/*
 * shatter_showcase.c - seams and shatter of curved surfaces (qaws_shatter.h)
 *
 * Writes PNG figures into ./shots/:
 *   shatter1_seams.png    seams: two bicubic patches crossing, a saddle and a
 *                         plane (two hyperbola branches), three patches whose
 *                         seams meet at a point
 *   shatter2_pieces.png   shatter: a wavy patch by a box (the piece inside
 *                         lifted out), a saddle by a wavy patch, a wavy patch
 *                         by two crossing patches (four pieces pulled apart)
 */

#include "qaws.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
#include <direct.h>
#define MAKE_DIR(p) _mkdir(p)
#else
#include <sys/stat.h>
#define MAKE_DIR(p) mkdir(p, 0755)
#endif

#include "example_raster.h"

static qaws_vec3 ss_v(double x, double y, double z)
{
	qaws_vec3 v;
	v.x = (qaws_scalar)x;
	v.y = (qaws_scalar)y;
	v.z = (qaws_scalar)z;
	return v;
}

static qaws_surface* ss_bilinear(qaws_vec3 p00, qaws_vec3 p10, qaws_vec3 p01, qaws_vec3 p11)
{
	qaws_surface_bilinear_desc d;
	qaws_surface* s = NULL;
	d.p00 = p00; d.p10 = p10; d.p01 = p01; d.p11 = p11;
	qaws_surface_create_bilinear(&s, &d);
	return s;
}

/* a bicubic patch over (a, b) in [-1, 1]^2, the third coordinate h[i][j] placed on axis (2: z, 0: x) */
static qaws_surface* ss_bicubic(double const h[4][4], int axis, double sa, double sb)
{
	qaws_vec3 cp[16];
	qaws_surface_bezier_desc d;
	qaws_surface* s = NULL;
	unsigned int i, j;
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
		{
			double a = sa * (-1 + 2.0 * i / 3), b = sb * (-1 + 2.0 * j / 3);
			cp[4 * i + j] = axis == 2 ? ss_v(a, b, h[i][j]) : axis == 0 ? ss_v(h[i][j], a, b) : ss_v(a, h[i][j], b);
		}
	memset(&d, 0, sizeof d);
	d.u_degree = 3;
	d.v_degree = 3;
	d.control_points = cp;
	d.u_point_count = 4;
	d.v_point_count = 4;
	qaws_surface_create_bezier(&s, &d);
	return s;
}

/* the six faces of a box, normals outward */
static void ss_box(double const* lo, double const* hi, qaws_surface** f)
{
	double x0 = lo[0], y0 = lo[1], z0 = lo[2], x1 = hi[0], y1 = hi[1], z1 = hi[2];
	f[0] = ss_bilinear(ss_v(x0, y0, z0), ss_v(x0, y1, z0), ss_v(x1, y0, z0), ss_v(x1, y1, z0));
	f[1] = ss_bilinear(ss_v(x0, y0, z1), ss_v(x1, y0, z1), ss_v(x0, y1, z1), ss_v(x1, y1, z1));
	f[2] = ss_bilinear(ss_v(x0, y0, z0), ss_v(x1, y0, z0), ss_v(x0, y0, z1), ss_v(x1, y0, z1));
	f[3] = ss_bilinear(ss_v(x0, y1, z0), ss_v(x0, y1, z1), ss_v(x1, y1, z0), ss_v(x1, y1, z1));
	f[4] = ss_bilinear(ss_v(x0, y0, z0), ss_v(x0, y0, z1), ss_v(x0, y1, z0), ss_v(x0, y1, z1));
	f[5] = ss_bilinear(ss_v(x1, y0, z0), ss_v(x1, y1, z0), ss_v(x1, y0, z1), ss_v(x1, y1, z1));
}

/* a surface (or a trimmed piece over base, its trim kept) on an n x n grid: positions moved by off, the triangles inside the
   trim; out: 9 doubles a triangle, the count */
static double* ss_triangles(qaws_surface const* s, qaws_surface const* base, unsigned int n, double const* off, unsigned int* count)
{
	qaws_range ru = qaws_surface_get_u_range(base), rv = qaws_surface_get_v_range(base);
	double* g = (double*)malloc(sizeof(double) * 3 * (n + 1) * (n + 1));
	double* t = (double*)malloc(sizeof(double) * 18 * n * n);
	unsigned int i, j, k, m = 0;
	for (j = 0; j <= n; j++)
		for (i = 0; i <= n; i++)
		{
			qaws_surface_eval_result r;
			double u = ru.min_value + (ru.max_value - ru.min_value) * i / n, v = rv.min_value + (rv.max_value - rv.min_value) * j / n;
			qaws_surface_evaluate(&r, base, (qaws_scalar)u, (qaws_scalar)v, QAWS_SURFACE_EVAL_POSITION);
			g[3 * (j * (n + 1) + i)] = r.position.x + off[0];
			g[3 * (j * (n + 1) + i) + 1] = r.position.y + off[1];
			g[3 * (j * (n + 1) + i) + 2] = r.position.z + off[2];
		}
	for (j = 0; j < n; j++)
		for (i = 0; i < n; i++)
		{
			unsigned int q[2][3] = { { j * (n + 1) + i, j * (n + 1) + i + 1, (j + 1) * (n + 1) + i + 1 }, { j * (n + 1) + i, (j + 1) * (n + 1) + i + 1, (j + 1) * (n + 1) + i } };
			double cu[2] = { (i + 2.0 / 3) / n, (i + 1.0 / 3) / n }, cv[2] = { (j + 1.0 / 3) / n, (j + 2.0 / 3) / n };
			unsigned int h;
			for (h = 0; h < 2; h++)
			{
				int in = 1;
				if (s != base)
					qaws_surface_trim_contains(&in, s, (qaws_scalar)(ru.min_value + (ru.max_value - ru.min_value) * cu[h]),
						(qaws_scalar)(rv.min_value + (rv.max_value - rv.min_value) * cv[h]));
				if (!in)
					continue;
				for (k = 0; k < 3; k++)
					memcpy(t + 9 * m + 3 * k, g + 3 * q[h][k], 3 * sizeof(double));
				m++;
			}
		}
	free(g);
	*count = m;
	return t;
}

typedef struct ss_draw
{
	double* tri[24];
	unsigned int n[24];
	double col[24][3];
	unsigned int count;
} ss_draw;

static void ss_add(ss_draw* D, qaws_surface const* s, qaws_surface const* base, unsigned int n, double const* off, double const* col)
{
	if (D->count == 24)
		return;
	D->tri[D->count] = ss_triangles(s, base, n, off, &D->n[D->count]);
	memcpy(D->col[D->count], col, sizeof D->col[0]);
	D->count++;
}

/* the panel fitted to all the triangles (and the seams' points), drawn, titled */
static void ss_panel(er_canvas* c, int x, int w, ss_draw* D, double const* const* lines, unsigned int const* line_n, unsigned int nl,
	char const* title, char const* sub)
{
	static double const seam[3] = { 200, 30, 45 }, wire[3] = { 70, 74, 86 };
	unsigned int i, k, tot = 0, at = 0;
	double* pts;
	for (i = 0; i < D->count; i++)
		tot += 3 * D->n[i];
	pts = (double*)malloc(sizeof(double) * 3 * (tot + 1));
	for (i = 0; i < D->count; i++)
	{
		memcpy(pts + 3 * at, D->tri[i], sizeof(double) * 9 * D->n[i]);
		at += 3 * D->n[i];
	}
	er_panel(c, x, 72, w, 628, pts, tot, -0.85, 0.5, 0.06);
	free(pts);
	er_text(c, x + 14, 12, title, 3, 30, 30, 30);
	er_text(c, x + 14, 46, sub, 2, 95, 95, 95);
	for (i = 0; i < D->count; i++)
		for (k = 0; k < D->n[i]; k++)
		{
			double const* t = D->tri[i] + 9 * k;
			er_triangle(c, t, t + 3, t + 6, NULL, NULL, NULL, D->col[i], D->col[i], D->col[i]);
		}
	for (i = 0; i < D->count; i++)
		for (k = 0; k < D->n[i]; k += 7)
		{
			double const* t = D->tri[i] + 9 * k;
			er_line(c, t, t + 3, wire, 0.2);
		}
	for (i = 0; i < nl; i++)
		for (k = 0; k + 1 < line_n[i]; k++)
			er_line(c, lines[i] + 3 * k, lines[i] + 3 * k + 3, seam, 4.5);
	for (i = 0; i < D->count; i++)
		free(D->tri[i]);
	D->count = 0;
}

/* the seams of A and B as point arrays (closed ones closed again), their count */
static unsigned int ss_seams(qaws_surface_seam_desc const* d, double** lines, unsigned int* line_n, unsigned int cap, unsigned int* points)
{
	qaws_surface_seams* S = NULL;
	unsigned int n, i, j;
	*points = 0;
	if (qaws_surface_seams_create(&S, d) != QAWS_STATUS_OK)
		return 0;
	n = qaws_surface_seams_get_count(S);
	if (n > cap)
		n = cap;
	for (i = 0; i < n; i++)
	{
		qaws_seam_info in;
		qaws_seam_point* p;
		qaws_surface_seams_get(&in, NULL, S, i);
		p = (qaws_seam_point*)malloc(sizeof(qaws_seam_point) * (in.point_count + 1));
		qaws_surface_seams_get(NULL, p, S, i);
		lines[i] = (double*)malloc(sizeof(double) * 3 * (in.point_count + 2));
		for (j = 0; j < in.point_count + (in.closed ? 1u : 0u); j++)
		{
			qaws_seam_point const* q = &p[j % in.point_count];
			lines[i][3 * j] = q->position.x;
			lines[i][3 * j + 1] = q->position.y;
			lines[i][3 * j + 2] = q->position.z;
		}
		line_n[i] = in.point_count + (in.closed ? 1u : 0u);
		*points += in.point_count;
		free(p);
	}
	qaws_surface_seams_destroy(S);
	return n;
}

static double const ss_wavy[4][4] = { { 0.05, 0.25, -0.15, 0.0 }, { 0.2, -0.25, 0.3, -0.1 }, { -0.1, 0.3, -0.25, 0.15 }, { 0.0, -0.15, 0.25, 0.05 } };
static double const ss_bend[4][4] = { { 0.1, -0.2, 0.2, -0.05 }, { -0.15, 0.25, -0.2, 0.1 }, { 0.2, -0.1, 0.25, -0.2 }, { -0.05, 0.15, -0.1, 0.1 } };

int main(void)
{
	static double const grey[3] = { 200, 202, 208 }, warm[3] = { 236, 170, 104 }, blue[3] = { 140, 172, 222 }, green[3] = { 140, 200, 140 };
	static double const pal[6][3] = { { 140, 172, 222 }, { 236, 170, 104 }, { 140, 200, 140 }, { 214, 140, 180 }, { 230, 214, 120 }, { 150, 200, 210 } };
	static double const zero[3] = { 0, 0, 0 };
	qaws_surface *W = ss_bicubic(ss_wavy, 2, 1, 1), *V = ss_bicubic(ss_bend, 0, 1, 0.8), *U = ss_bicubic(ss_bend, 1, 1, 0.8);
	qaws_surface *saddle = ss_bilinear(ss_v(-1, -1, 0.5), ss_v(1, -1, -0.5), ss_v(-1, 1, -0.5), ss_v(1, 1, 0.5));
	qaws_surface* plane = ss_bilinear(ss_v(-1.3, -1.3, 0.2), ss_v(1.3, -1.3, 0.2), ss_v(-1.3, 1.3, 0.2), ss_v(1.3, 1.3, 0.2));
	qaws_surface* box[6];
	double const lo[3] = { -0.45, -0.4, -0.3 }, hi[3] = { 0.45, 0.4, 0.3 };
	ss_draw D;
	double* lines[64];
	unsigned int line_n[64], nl, pts, i;
	char sub[160];
	er_canvas* c;
	qaws_surface_seam_desc d;
	memset(&D, 0, sizeof D);
	MAKE_DIR("shots");
	ss_box(lo, hi, box);

	/* ---- seams ---- */
	c = er_create(2100, 700, 3, 252, 252, 250);
	memset(&d, 0, sizeof d);
	d.kinds = QAWS_SEAM_A_B;
	d.a = (qaws_surface const* const*)&W; d.a_count = 1;
	d.b = (qaws_surface const* const*)&V; d.b_count = 1;
	nl = ss_seams(&d, lines, line_n, 64, &pts);
	sprintf(sub, "%u seam%s, %u points on both", nl, nl == 1 ? "" : "s", pts);
	ss_add(&D, W, W, 90, zero, grey);
	ss_add(&D, V, V, 90, zero, warm);
	ss_panel(c, 0, 700, &D, (double const* const*)lines, line_n, nl, "two bicubic patches", sub);
	printf("  seams, two bicubic patches: %s\n", sub);
	for (i = 0; i < nl; i++)
		free(lines[i]);
	d.a = (qaws_surface const* const*)&saddle;
	d.b = (qaws_surface const* const*)&plane;
	nl = ss_seams(&d, lines, line_n, 64, &pts);
	sprintf(sub, "%u seams (hyperbola branches), %u points", nl, pts);
	ss_add(&D, saddle, saddle, 90, zero, green);
	ss_add(&D, plane, plane, 30, zero, grey);
	ss_panel(c, 700, 700, &D, (double const* const*)lines, line_n, nl, "a saddle and a plane", sub);
	printf("  seams, a saddle and a plane: %s\n", sub);
	for (i = 0; i < nl; i++)
		free(lines[i]);
	{
		qaws_surface* two[2];
		two[0] = V;
		two[1] = U;
		d.a = (qaws_surface const* const*)&W;
		d.b = (qaws_surface const* const*)two;
		d.b_count = 2;
		d.kinds = QAWS_SEAM_A_B | QAWS_SEAM_B_B;
		nl = ss_seams(&d, lines, line_n, 64, &pts);
		sprintf(sub, "%u seams meeting, %u points", nl, pts);
		ss_add(&D, W, W, 90, zero, grey);
		ss_add(&D, V, V, 90, zero, warm);
		ss_add(&D, U, U, 90, zero, blue);
		ss_panel(c, 1400, 700, &D, (double const* const*)lines, line_n, nl, "three patches", sub);
		printf("  seams, three patches: %s\n", sub);
		for (i = 0; i < nl; i++)
			free(lines[i]);
	}
	er_save(c, "shots/shatter1_seams.png");
	printf("shatter1_seams.png\n");
	er_destroy(c);

	/* ---- shatter ---- */
	c = er_create(2100, 700, 3, 252, 252, 250);
	{
		qaws_surface_shatter* H = NULL;
		unsigned int n;
		memset(&d, 0, sizeof d);
		d.a = (qaws_surface const* const*)&W; d.a_count = 1;
		d.b = (qaws_surface const* const*)box; d.b_count = 6;
		qaws_surface_shatter_create(&H, &d, 1);
		n = qaws_surface_shatter_get_count(H);
		for (i = 0; i < n; i++)
		{
			qaws_surface const* p;
			int inside;
			double const up[3] = { 0, 0, 0.45 };
			qaws_surface_shatter_get(&p, NULL, &inside, H, i);
			ss_add(&D, p, W, 140, inside == 1 ? up : zero, inside == 1 ? warm : blue);
		}
		sprintf(sub, "%u pieces, the one inside the box lifted", n);
		ss_panel(c, 0, 700, &D, NULL, NULL, 0, "a wavy patch by a box", sub);
		{
			/* (the box's edges) */
			static double const edge_col[3] = { 90, 90, 100 };
			unsigned int e;
			for (e = 0; e < 12; e++)
			{
				double a[3], b[3];
				unsigned int ax = e / 4, k1 = e % 4;
				a[ax] = lo[ax];
				b[ax] = hi[ax];
				a[(ax + 1) % 3] = b[(ax + 1) % 3] = (k1 & 1) ? hi[(ax + 1) % 3] : lo[(ax + 1) % 3];
				a[(ax + 2) % 3] = b[(ax + 2) % 3] = (k1 & 2) ? hi[(ax + 2) % 3] : lo[(ax + 2) % 3];
				er_line(c, a, b, edge_col, 1.6);
			}
		}
		printf("  shatter, a wavy patch by a box: %s\n", sub);
		qaws_surface_shatter_destroy(H);
		H = NULL;
		d.a = (qaws_surface const* const*)&saddle;
		d.b = (qaws_surface const* const*)&W; d.b_count = 1;
		qaws_surface_shatter_create(&H, &d, 0);
		n = qaws_surface_shatter_get_count(H);
		for (i = 0; i < n; i++)
		{
			qaws_surface const* p;
			double off[3] = { 0, 0, 0.18 * (double)i };
			qaws_surface_shatter_get(&p, NULL, NULL, H, i);
			ss_add(&D, p, saddle, 140, off, pal[i % 6]);
		}
		sprintf(sub, "%u pieces, pulled apart", n);
		ss_panel(c, 700, 700, &D, NULL, NULL, 0, "a saddle by a wavy patch", sub);
		printf("  shatter, a saddle by a wavy patch: %s\n", sub);
		qaws_surface_shatter_destroy(H);
		H = NULL;
		{
			qaws_surface* two[2];
			two[0] = V;
			two[1] = U;
			d.a = (qaws_surface const* const*)&W;
			d.b = (qaws_surface const* const*)two;
			d.b_count = 2;
			qaws_surface_shatter_create(&H, &d, 0);
			n = qaws_surface_shatter_get_count(H);
			for (i = 0; i < n; i++)
			{
				qaws_surface const* p;
				qaws_surface_eval_result r;
				double off[3];
				qaws_surface_shatter_get(&p, NULL, NULL, H, i);
				/* (each piece moved out from the middle, by where a point of it lies) */
				{
					unsigned int m;
					double* t = ss_triangles(p, W, 40, zero, &m);
					double gx = 0, gy = 0;
					unsigned int k;
					for (k = 0; k < m; k++)
					{
						gx += (t[9 * k] + t[9 * k + 3] + t[9 * k + 6]) / 3;
						gy += (t[9 * k + 1] + t[9 * k + 4] + t[9 * k + 7]) / 3;
					}
					free(t);
					off[0] = m ? 0.35 * (gx / m > 0.15 ? 1 : -1) : 0;
					off[1] = m ? 0.35 * (gy / m > 0 ? 1 : -1) : 0;
					off[2] = 0;
				}
				(void)r;
				ss_add(&D, p, W, 140, off, pal[i % 6]);
			}
			sprintf(sub, "%u pieces, pulled apart", n);
			ss_panel(c, 1400, 700, &D, NULL, NULL, 0, "a wavy patch by two patches", sub);
			printf("  shatter, a wavy patch by two patches: %s\n", sub);
			qaws_surface_shatter_destroy(H);
		}
	}
	er_save(c, "shots/shatter2_pieces.png");
	printf("shatter2_pieces.png\n");
	er_destroy(c);
	qaws_surface_destroy(W);
	qaws_surface_destroy(V);
	qaws_surface_destroy(U);
	qaws_surface_destroy(saddle);
	qaws_surface_destroy(plane);
	for (i = 0; i < 6; i++)
		qaws_surface_destroy(box[i]);
	return 0;
}
