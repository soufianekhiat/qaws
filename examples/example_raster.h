#ifndef EXAMPLE_RASTER_H
#define EXAMPLE_RASTER_H

/*
 * A small z-buffered rasterizer for the meshing figures: shaded triangles
 * (smooth or flat, colour per face or per vertex), hidden-line wireframes
 * (edges depth-tested against the faces), points, labels in a 5x7 bitmap
 * font, panels in one canvas, supersampled and written as an RGBA PNG
 * (stored deflate blocks: no zlib). Orthographic views, z up.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct er_canvas
{
	int w, h, ss, W, H;             /* output size, supersampling, internal size */
	unsigned char* rgb;             /* W x H */
	float* z;
	/* the current panel's view */
	int px, py, pw, ph;
	double rot[9], scale, cx, cy, depth_bias;
	double key[3], fill[3], half[3];
} er_canvas;

static er_canvas* er_create(int w, int h, int ss, int r, int g, int b)
{
	er_canvas* c = (er_canvas*)calloc(1, sizeof(er_canvas));
	int i;
	double l;
	if (!c)
		return NULL;
	c->w = w; c->h = h; c->ss = ss; c->W = w * ss; c->H = h * ss;
	c->rgb = (unsigned char*)malloc(3 * (size_t)c->W * c->H);
	c->z = (float*)malloc(sizeof(float) * (size_t)c->W * c->H);
	if (!c->rgb || !c->z)
	{
		free(c->rgb); free(c->z); free(c);
		return NULL;
	}
	for (i = 0; i < c->W * c->H; i++)
	{
		c->rgb[3 * (size_t)i] = (unsigned char)r; c->rgb[3 * (size_t)i + 1] = (unsigned char)g; c->rgb[3 * (size_t)i + 2] = (unsigned char)b;
		c->z[i] = 1e30f;
	}
	/* lights in view space (x right, y up, z toward the eye) */
	c->key[0] = -0.45; c->key[1] = 0.6; c->key[2] = 0.66;
	c->fill[0] = 0.6; c->fill[1] = -0.1; c->fill[2] = 0.8;
	l = sqrt(c->key[0] * c->key[0] + c->key[1] * c->key[1] + c->key[2] * c->key[2]);
	for (i = 0; i < 3; i++) c->key[i] /= l;
	l = sqrt(c->fill[0] * c->fill[0] + c->fill[1] * c->fill[1] + c->fill[2] * c->fill[2]);
	for (i = 0; i < 3; i++) c->fill[i] /= l;
	c->half[0] = c->key[0]; c->half[1] = c->key[1]; c->half[2] = c->key[2] + 1;
	l = sqrt(c->half[0] * c->half[0] + c->half[1] * c->half[1] + c->half[2] * c->half[2]);
	for (i = 0; i < 3; i++) c->half[i] /= l;
	return c;
}

static void er_destroy(er_canvas* c)
{
	if (!c)
		return;
	free(c->rgb); free(c->z); free(c);
}

/* world to view: turn by yaw round z, tilt by pitch round x; view x right, y up, z toward the eye */
static void er_view_rotation(double yaw, double pitch, double* r)
{
	double cy = cos(yaw), sy = sin(yaw), cp = cos(pitch), sp = sin(pitch);
	/* after the yaw: x' = cy x - sy y, y' = sy x + cy y (depth), z' = z */
	r[0] = cy; r[1] = -sy; r[2] = 0;
	r[3] = sp * sy; r[4] = sp * cy; r[5] = cp;       /* screen up: z tilted toward the viewer */
	r[6] = -cp * sy; r[7] = -cp * cy; r[8] = sp;     /* toward the eye */
}

static void er_to_view(er_canvas const* c, double const* p, double* q)
{
	q[0] = c->rot[0] * p[0] + c->rot[1] * p[1] + c->rot[2] * p[2];
	q[1] = c->rot[3] * p[0] + c->rot[4] * p[1] + c->rot[5] * p[2];
	q[2] = c->rot[6] * p[0] + c->rot[7] * p[1] + c->rot[8] * p[2];
}

/* view to internal pixels: x right, y down, depth (smaller is nearer) */
static void er_project(er_canvas const* c, double const* p, double* s)
{
	double q[3];
	er_to_view(c, p, q);
	s[0] = c->px + c->pw * 0.5 + (q[0] - c->cx) * c->scale;
	s[1] = c->py + c->ph * 0.5 - (q[1] - c->cy) * c->scale;
	s[2] = -q[2] * c->scale;
}

/* a panel at output pixels (x, y, w, h) showing points x[0..n) (3 each) fitted with a margin */
static void er_panel(er_canvas* c, int x, int y, int w, int h, double const* pts, unsigned int n, double yaw, double pitch, double margin)
{
	double lo[2] = { 1e300, 1e300 }, hi[2] = { -1e300, -1e300 }, q[3], ext;
	unsigned int i;
	c->px = x * c->ss; c->py = y * c->ss; c->pw = w * c->ss; c->ph = h * c->ss;
	er_view_rotation(yaw, pitch, c->rot);
	for (i = 0; i < n; i++)
	{
		er_to_view(c, pts + 3 * (size_t)i, q);
		if (q[0] < lo[0]) lo[0] = q[0];
		if (q[1] < lo[1]) lo[1] = q[1];
		if (q[0] > hi[0]) hi[0] = q[0];
		if (q[1] > hi[1]) hi[1] = q[1];
	}
	c->cx = 0.5 * (lo[0] + hi[0]);
	c->cy = 0.5 * (lo[1] + hi[1]);
	ext = fmax(hi[0] - lo[0], (hi[1] - lo[1]) * w / (double)h);
	c->scale = ext > 0 ? (1 - 2 * margin) * c->pw / ext : 1;
	c->depth_bias = 0.004 * c->scale * (ext > 0 ? ext : 1);
}

/* Blinn-Phong with a key and a fill light; n the unit view-space normal facing the eye */
static void er_shade(er_canvas const* c, double const* n, double const* col, double* out)
{
	double dk = fmax(0, n[0] * c->key[0] + n[1] * c->key[1] + n[2] * c->key[2]);
	double df = fmax(0, n[0] * c->fill[0] + n[1] * c->fill[1] + n[2] * c->fill[2]);
	double sp = pow(fmax(0, n[0] * c->half[0] + n[1] * c->half[1] + n[2] * c->half[2]), 40);
	int k;
	for (k = 0; k < 3; k++)
		out[k] = fmin(255, col[k] * (0.30 + 0.62 * dk + 0.22 * df) + 255 * 0.28 * sp);
}

/* a triangle; normals (world, may be NULL: flat) and colours per vertex (cols[3 * k]) */
static void er_triangle(er_canvas* c, double const* p0, double const* p1, double const* p2, double const* n0, double const* n1,
	double const* n2, double const* c0, double const* c1, double const* c2)
{
	double A[3], B[3], C[3], area, fn[3], u[3], v[3], vn[3][3], gn[3];
	double const* P[3] = { p0, p1, p2 };
	double const* N[3] = { n0, n1, n2 };
	int xa, xb, ya, yb, x, y, k, j;
	er_project(c, p0, A); er_project(c, p1, B); er_project(c, p2, C);
	area = (B[0] - A[0]) * (C[1] - A[1]) - (B[1] - A[1]) * (C[0] - A[0]);
	if (fabs(area) < 1e-12)
		return;
	for (k = 0; k < 3; k++) { u[k] = P[1][k] - P[0][k]; v[k] = P[2][k] - P[0][k]; }
	fn[0] = u[1] * v[2] - u[2] * v[1]; fn[1] = u[2] * v[0] - u[0] * v[2]; fn[2] = u[0] * v[1] - u[1] * v[0];
	{
		/* the geometric normal in view space: which side faces the eye */
		double l;
		er_to_view(c, fn, gn);
		l = sqrt(gn[0] * gn[0] + gn[1] * gn[1] + gn[2] * gn[2]);
		if (l > 0)
			for (k = 0; k < 3; k++) gn[k] /= l;
	}
	for (j = 0; j < 3; j++)
	{
		double const* src = N[j] ? N[j] : fn;
		double w[3], l;
		er_to_view(c, src, w);
		l = sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
		for (k = 0; k < 3; k++) vn[j][k] = l > 0 ? w[k] / l : 0;
	}
	xa = (int)fmax(c->px, floor(fmin(A[0], fmin(B[0], C[0])))); xb = (int)fmin(c->px + c->pw - 1, ceil(fmax(A[0], fmax(B[0], C[0]))));
	ya = (int)fmax(c->py, floor(fmin(A[1], fmin(B[1], C[1])))); yb = (int)fmin(c->py + c->ph - 1, ceil(fmax(A[1], fmax(B[1], C[1]))));
	for (y = ya; y <= yb; y++)
		for (x = xa; x <= xb; x++)
		{
			double sx = x + 0.5, sy = y + 0.5, l0, l1, l2, d, n[3], nl, col[3], out[3];
			size_t i;
			l0 = ((B[0] - sx) * (C[1] - sy) - (B[1] - sy) * (C[0] - sx)) / area;
			l1 = ((C[0] - sx) * (A[1] - sy) - (C[1] - sy) * (A[0] - sx)) / area;
			l2 = 1 - l0 - l1;
			if (l0 < 0 || l1 < 0 || l2 < 0)
				continue;
			i = (size_t)y * c->W + x;
			d = l0 * A[2] + l1 * B[2] + l2 * C[2];
			if (d >= c->z[i])
				continue;
			c->z[i] = (float)d;
			for (k = 0; k < 3; k++)
			{
				n[k] = l0 * vn[0][k] + l1 * vn[1][k] + l2 * vn[2][k];
				col[k] = l0 * c0[k] + l1 * c1[k] + l2 * c2[k];
			}
			nl = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
			if (nl == 0)
				continue;
			for (k = 0; k < 3; k++) n[k] /= nl;
			if (gn[2] < 0)
			{
				/* the back side: darker */
				for (k = 0; k < 3; k++) { n[k] = -n[k]; col[k] *= 0.6; }
			}
			if (n[2] < 0)
			{
				/* a shading normal turned away on a side facing the eye: the face's own */
				for (k = 0; k < 3; k++) n[k] = gn[2] < 0 ? -gn[k] : gn[k];
			}
			er_shade(c, n, col, out);
			for (k = 0; k < 3; k++) c->rgb[3 * i + k] = (unsigned char)out[k];
		}
}

/* a line of the given width (output pixels), hidden where faces are in front */
static void er_line(er_canvas* c, double const* a, double const* b, double const* col, double width)
{
	double A[3], B[3], len, r = 0.5 * width * c->ss;
	int steps, s, x, y, k;
	er_project(c, a, A);
	er_project(c, b, B);
	len = sqrt((B[0] - A[0]) * (B[0] - A[0]) + (B[1] - A[1]) * (B[1] - A[1]));
	steps = (int)ceil(len / fmax(0.5, r * 0.5)) + 1;
	for (s = 0; s <= steps; s++)
	{
		double t = (double)s / steps, X = A[0] + t * (B[0] - A[0]), Y = A[1] + t * (B[1] - A[1]), D = A[2] + t * (B[2] - A[2]);
		int xa = (int)floor(X - r), xb = (int)ceil(X + r), ya = (int)floor(Y - r), yb = (int)ceil(Y + r);
		for (y = ya; y <= yb; y++)
			for (x = xa; x <= xb; x++)
			{
				size_t i;
				if (x < c->px || y < c->py || x >= c->px + c->pw || y >= c->py + c->ph)
					continue;
				if ((x + 0.5 - X) * (x + 0.5 - X) + (y + 0.5 - Y) * (y + 0.5 - Y) > r * r)
					continue;
				i = (size_t)y * c->W + x;
				if (D > c->z[i] + c->depth_bias)
					continue;
				for (k = 0; k < 3; k++) c->rgb[3 * i + k] = (unsigned char)col[k];
			}
	}
}

/* a shaded ball of radius r (output pixels) */
static void er_point(er_canvas* c, double const* p, double const* col, double r)
{
	double P[3], R = r * c->ss;
	int x, y, k;
	er_project(c, p, P);
	for (y = (int)floor(P[1] - R); y <= (int)ceil(P[1] + R); y++)
		for (x = (int)floor(P[0] - R); x <= (int)ceil(P[0] + R); x++)
		{
			double u = (x + 0.5 - P[0]) / R, v = (y + 0.5 - P[1]) / R, s2 = u * u + v * v, n[3], out[3], d;
			size_t i;
			if (s2 >= 1 || x < c->px || y < c->py || x >= c->px + c->pw || y >= c->py + c->ph)
				continue;
			n[0] = u; n[1] = -v; n[2] = sqrt(1 - s2);
			d = P[2] - n[2] * R - c->depth_bias;
			i = (size_t)y * c->W + x;
			if (d >= c->z[i])
				continue;
			c->z[i] = (float)d;
			er_shade(c, n, col, out);
			for (k = 0; k < 3; k++) c->rgb[3 * i + k] = (unsigned char)out[k];
		}
}

/* the classic 5x7 font, columns, bit 0 at the top: ' ' to '~' */
static unsigned char const er_font[95][5] = {
	{0,0,0,0,0},{0,0,0x5F,0,0},{0,7,0,7,0},{0x14,0x7F,0x14,0x7F,0x14},{0x24,0x2A,0x7F,0x2A,0x12},{0x23,0x13,8,0x64,0x62},
	{0x36,0x49,0x55,0x22,0x50},{0,5,3,0,0},{0,0x1C,0x22,0x41,0},{0,0x41,0x22,0x1C,0},{8,0x2A,0x1C,0x2A,8},{8,8,0x3E,8,8},
	{0,0x50,0x30,0,0},{8,8,8,8,8},{0,0x60,0x60,0,0},{0x20,0x10,8,4,2},{0x3E,0x51,0x49,0x45,0x3E},{0,0x42,0x7F,0x40,0},
	{0x42,0x61,0x51,0x49,0x46},{0x21,0x41,0x45,0x4B,0x31},{0x18,0x14,0x12,0x7F,0x10},{0x27,0x45,0x45,0x45,0x39},
	{0x3C,0x4A,0x49,0x49,0x30},{1,0x71,9,5,3},{0x36,0x49,0x49,0x49,0x36},{6,0x49,0x49,0x29,0x1E},{0,0x36,0x36,0,0},
	{0,0x56,0x36,0,0},{0,8,0x14,0x22,0x41},{0x14,0x14,0x14,0x14,0x14},{0x41,0x22,0x14,8,0},{2,1,0x51,9,6},
	{0x32,0x49,0x79,0x41,0x3E},{0x7E,0x11,0x11,0x11,0x7E},{0x7F,0x49,0x49,0x49,0x36},{0x3E,0x41,0x41,0x41,0x22},
	{0x7F,0x41,0x41,0x22,0x1C},{0x7F,0x49,0x49,0x49,0x41},{0x7F,9,9,1,1},{0x3E,0x41,0x41,0x51,0x32},{0x7F,8,8,8,0x7F},
	{0,0x41,0x7F,0x41,0},{0x20,0x40,0x41,0x3F,1},{0x7F,8,0x14,0x22,0x41},{0x7F,0x40,0x40,0x40,0x40},{0x7F,2,4,2,0x7F},
	{0x7F,4,8,0x10,0x7F},{0x3E,0x41,0x41,0x41,0x3E},{0x7F,9,9,9,6},{0x3E,0x41,0x51,0x21,0x5E},{0x7F,9,0x19,0x29,0x46},
	{0x46,0x49,0x49,0x49,0x31},{1,1,0x7F,1,1},{0x3F,0x40,0x40,0x40,0x3F},{0x1F,0x20,0x40,0x20,0x1F},{0x7F,0x20,0x18,0x20,0x7F},
	{0x63,0x14,8,0x14,0x63},{3,4,0x78,4,3},{0x61,0x51,0x49,0x45,0x43},{0,0,0x7F,0x41,0x41},{2,4,8,0x10,0x20},
	{0x41,0x41,0x7F,0,0},{4,2,1,2,4},{0x40,0x40,0x40,0x40,0x40},{0,1,2,4,0},{0x20,0x54,0x54,0x54,0x78},
	{0x7F,0x48,0x44,0x44,0x38},{0x38,0x44,0x44,0x44,0x20},{0x38,0x44,0x44,0x48,0x7F},{0x38,0x54,0x54,0x54,0x18},
	{8,0x7E,9,1,2},{8,0x14,0x54,0x54,0x3C},{0x7F,8,4,4,0x78},{0,0x44,0x7D,0x40,0},{0x20,0x40,0x44,0x3D,0},
	{0,0x7F,0x10,0x28,0x44},{0,0x41,0x7F,0x40,0},{0x7C,4,0x18,4,0x78},{0x7C,8,4,4,0x78},{0x38,0x44,0x44,0x44,0x38},
	{0x7C,0x14,0x14,0x14,8},{8,0x14,0x14,0x18,0x7C},{0x7C,8,4,4,8},{0x48,0x54,0x54,0x54,0x20},{4,0x3F,0x44,0x40,0x20},
	{0x3C,0x40,0x40,0x20,0x7C},{0x1C,0x20,0x40,0x20,0x1C},{0x3C,0x40,0x30,0x40,0x3C},{0x44,0x28,0x10,0x28,0x44},
	{0x0C,0x50,0x50,0x50,0x3C},{0x44,0x64,0x54,0x4C,0x44},{0,8,0x36,0x41,0},{0,0,0x7F,0,0},{0,0x41,0x36,8,0},
	{8,4,8,0x10,8} };

/* text at output pixel (x, y) (its top left), each font pixel `scale` output pixels */
static void er_text(er_canvas* c, int x, int y, char const* s, int scale, int r, int g, int b)
{
	int k = 0, col, row, u, v, X, Y;
	for (; *s; s++, k++)
	{
		unsigned char ch = (unsigned char)*s;
		if (ch < 32 || ch > 126)
			continue;
		for (col = 0; col < 5; col++)
			for (row = 0; row < 7; row++)
			{
				if (!(er_font[ch - 32][col] >> row & 1))
					continue;
				for (v = 0; v < scale * c->ss; v++)
					for (u = 0; u < scale * c->ss; u++)
					{
						size_t i;
						X = (x + (k * 6 + col) * scale) * c->ss + u;
						Y = (y + row * scale) * c->ss + v;
						if (X < 0 || Y < 0 || X >= c->W || Y >= c->H)
							continue;
						i = (size_t)Y * c->W + X;
						c->rgb[3 * i] = (unsigned char)r; c->rgb[3 * i + 1] = (unsigned char)g; c->rgb[3 * i + 2] = (unsigned char)b;
					}
			}
	}
}

/* ---- PNG, stored deflate blocks ---------------------------------------- */

static unsigned int er_crc(unsigned int crc, unsigned char const* p, size_t n)
{
	static unsigned int table[256];
	static int ready = 0;
	size_t i;
	if (!ready)
	{
		unsigned int c, k, b;
		for (k = 0; k < 256; k++)
		{
			for (c = k, b = 0; b < 8; b++) c = c & 1 ? 0xedb88320u ^ (c >> 1) : c >> 1;
			table[k] = c;
		}
		ready = 1;
	}
	crc = ~crc;
	for (i = 0; i < n; i++) crc = table[(crc ^ p[i]) & 0xff] ^ (crc >> 8);
	return ~crc;
}

static void er_be32(unsigned char* p, unsigned int v)
{
	p[0] = (unsigned char)(v >> 24); p[1] = (unsigned char)(v >> 16); p[2] = (unsigned char)(v >> 8); p[3] = (unsigned char)v;
}

static void er_chunk(FILE* f, char const* type, unsigned char const* data, size_t n)
{
	unsigned char hdr[8];
	er_be32(hdr, (unsigned int)n);
	memcpy(hdr + 4, type, 4);
	fwrite(hdr, 1, 8, f);
	if (n)
		fwrite(data, 1, n, f);
	er_be32(hdr, er_crc(er_crc(0, hdr + 4, 4), data, n));
	fwrite(hdr, 1, 4, f);
}

/* the canvas box-filtered to its output size and written as an RGB PNG */
static int er_save(er_canvas const* c, char const* path)
{
	static unsigned char const sig[8] = { 137, 80, 78, 71, 13, 10, 26, 10 };
	size_t row = 3 * (size_t)c->w + 1, raw_n = (size_t)c->h * row, blocks = (raw_n + 65534) / 65535, pos = 0, i;
	unsigned char* raw = (unsigned char*)malloc(raw_n);
	unsigned char* z = (unsigned char*)malloc(2 + raw_n + 5 * blocks + 4);
	unsigned char ihdr[13];
	unsigned int a = 1, b = 0;
	FILE* f = fopen(path, "wb");
	int x, y, u, v, k;
	if (!f || !raw || !z)
	{
		if (f) fclose(f);
		free(raw); free(z);
		return 0;
	}
	for (y = 0; y < c->h; y++)
	{
		raw[(size_t)y * row] = 0;
		for (x = 0; x < c->w; x++)
			for (k = 0; k < 3; k++)
			{
				unsigned int s = 0;
				for (v = 0; v < c->ss; v++)
					for (u = 0; u < c->ss; u++)
						s += c->rgb[3 * ((size_t)(y * c->ss + v) * c->W + x * c->ss + u) + k];
				raw[(size_t)y * row + 1 + 3 * (size_t)x + k] = (unsigned char)((s + c->ss * c->ss / 2) / (c->ss * c->ss));
			}
	}
	z[pos++] = 0x78; z[pos++] = 0x01;
	for (i = 0; i < raw_n; i += 65535)
	{
		size_t len = raw_n - i < 65535 ? raw_n - i : 65535;
		z[pos++] = (unsigned char)(i + len == raw_n);
		z[pos++] = (unsigned char)(len & 0xff); z[pos++] = (unsigned char)(len >> 8);
		z[pos++] = (unsigned char)(~len & 0xff); z[pos++] = (unsigned char)((~len >> 8) & 0xff);
		memcpy(z + pos, raw + i, len);
		pos += len;
	}
	for (i = 0; i < raw_n; i++) { a = (a + raw[i]) % 65521; b = (b + a) % 65521; }
	er_be32(z + pos, (b << 16) | a);
	pos += 4;
	fwrite(sig, 1, 8, f);
	er_be32(ihdr, (unsigned int)c->w); er_be32(ihdr + 4, (unsigned int)c->h);
	ihdr[8] = 8; ihdr[9] = 2; ihdr[10] = 0; ihdr[11] = 0; ihdr[12] = 0;   /* 8-bit RGB */
	er_chunk(f, "IHDR", ihdr, 13);
	er_chunk(f, "IDAT", z, pos);
	er_chunk(f, "IEND", NULL, 0);
	fclose(f);
	free(raw); free(z);
	return 1;
}

#endif /* EXAMPLE_RASTER_H */
