/*
 * Test 187: texture attributes (qaws_mesh_texture.h), small textures
 *
 *   - a linear field on a unit square of two triangles (uv its x and y) baked into 16 x 16 texels: every texel the field
 *     at its centre, all covered; read back between the centres exactly (bilinear of a linear field), and sampled onto
 *     a finer mesh's corners exactly there;
 *   - a label per face (NEAREST) baked and read: each texel its triangle's label, never a mix;
 *   - a unit normal (NORMAL) from two corners' directions: every texel unit;
 *   - a chart on the lower left quarter of the square, a gutter of 2: two rings grown past its border (coverage 2), their
 *     values within the chart's range, the rest untouched;
 *   - repeated along u: a triangle across u = 1 (its corners at u 0.9 and 0.1 + 1) written on both sides of the edge;
 *   - only the faces of one texture id baked;
 *   - refusals: an attribute of other components, no uv, an edge attribute;
 *   - from a detailed mesh (qaws_mesh_3d_bake_texture_from): a wave on a 32 x 32 grid onto two triangles' 32 x 32 texture,
 *     each texel the fine mesh's interpolation at its point exactly (the two triangles' vertices could not hold it);
 *     lifted 0.01 off, max_distance 0.005 writes none, 0.02 all;
 *   - texture ids through a Boolean: two boxes (ids 0 and 1, a linear and a constant texture) united, each texture read onto
 *     its own faces: A's corners its linear field exactly (cut corners too), B's its constant;
 *   - a surface patch's texture: a generated cylinder, a periodic field over its domain repeated round it, read through the
 *     footpoints within the bilinear error, seam vertices too; vertices without footpoints skip their faces, or given the
 *     surface are read at their closest points.
 */

#include "test_common.h"
#include "qaws_mesh_texture.h"
#include "qaws_mesh_boolean.h"
#include "qaws_mesh_generate.h"
#include "qaws_bezier.h"
#include "qaws_surface.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

/* an n x n grid of the unit square in the xy plane (2 n^2 triangles), "uv" per corner its x and y scaled by s */
static qaws_mesh_3d* tx_grid(unsigned int n, double s)
{
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (n + 1) * (n + 1));
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 6 * n * n);
	qaws_scalar* uv;
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	qaws_attribute_desc ad;
	unsigned int i, j, k = 0, id, nc;
	for (j = 0; j <= n; j++)
		for (i = 0; i <= n; i++)
		{
			x[3 * (j * (n + 1) + i)] = (qaws_scalar)((double)i / n);
			x[3 * (j * (n + 1) + i) + 1] = (qaws_scalar)((double)j / n);
			x[3 * (j * (n + 1) + i) + 2] = 0;
		}
	for (j = 0; j < n; j++)
		for (i = 0; i < n; i++)
		{
			unsigned int a = j * (n + 1) + i;
			t[k++] = a; t[k++] = a + 1; t[k++] = a + n + 2;
			t[k++] = a; t[k++] = a + n + 2; t[k++] = a + n + 1;
		}
	memset(&d, 0, sizeof d);
	d.positions = x; d.vertex_count = (n + 1) * (n + 1); d.face_vertices = t; d.face_count = 2 * n * n;
	qaws_mesh_3d_create(&m, &d);
	nc = qaws_mesh_3d_get_corner_count(m);
	uv = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * nc);
	for (i = 0; i < nc; i++)
	{
		uv[2 * i] = (qaws_scalar)(s * (double)x[3 * t[i]]);
		uv[2 * i + 1] = (qaws_scalar)(s * (double)x[3 * t[i] + 1]);
	}
	memset(&ad, 0, sizeof ad);
	ad.name = "uv"; ad.element = QAWS_MESH_CORNER; ad.components = 2; ad.rule = QAWS_ATTRIBUTE_LINEAR;
	qaws_mesh_3d_add_attribute(&id, m, uv, &ad);
	free(x); free(t); free(uv);
	return m;
}

static void tx_add(qaws_mesh_3d* m, char const* name, qaws_mesh_element el, unsigned int n, qaws_attribute_rule rule, qaws_scalar const* v)
{
	qaws_attribute_desc ad;
	unsigned int id;
	memset(&ad, 0, sizeof ad);
	ad.name = name; ad.element = el; ad.components = n; ad.rule = rule;
	qaws_mesh_3d_add_attribute(&id, m, v, &ad);
}

static double tx_field(double x, double y)
{
	return 3 * x + 2 * y + 1;
}

static void test_linear(void)
{
	qaws_mesh_3d *m = tx_grid(1, 1), *fine = tx_grid(4, 1);
	qaws_scalar f[4], tex[256], uv[2], y[1];
	unsigned char cov[256];
	qaws_texture t;
	unsigned int i, j, covered = 0, id, nc;
	double e_tex = 0, e_read = 0, e_fine = 0;
	qaws_status st;
	char msg[300];
	for (i = 0; i < 4; i++)
	{
		qaws_vec3 p = qaws_mesh_3d_get_vertex(m, i);
		f[i] = (qaws_scalar)tx_field(p.x, p.y);
	}
	tx_add(m, "f", QAWS_MESH_VERTEX, 1, QAWS_ATTRIBUTE_LINEAR, f);
	memset(&t, 0, sizeof t);
	t.texels = tex; t.width = 16; t.height = 16; t.components = 1;
	st = qaws_mesh_3d_bake_texture(cov, &t, m, "f", NULL);
	for (j = 0; j < 16; j++)
		for (i = 0; i < 16; i++)
		{
			covered += cov[j * 16 + i] == 1;
			e_tex = fmax(e_tex, fabs((double)tex[j * 16 + i] - tx_field((i + 0.5) / 16, (j + 0.5) / 16)));
		}
	for (i = 0; i < 50; i++)
	{
		uv[0] = (qaws_scalar)(1.0 / 32 + fmod(i * 0.6180339887, 1.0) * 30.0 / 32);
		uv[1] = (qaws_scalar)(1.0 / 32 + fmod(i * 0.7548776662, 1.0) * 30.0 / 32);
		qaws_texture_sample(y, &t, uv, 1);
		e_read = fmax(e_read, fabs((double)y[0] - tx_field(uv[0], uv[1])));
	}
	/* sampled onto a finer mesh's corners (those inside the texel centres' square) */
	if (qaws_mesh_3d_sample_texture(fine, &t, "f", NULL) == QAWS_STATUS_OK && qaws_mesh_3d_find_attribute(&id, fine, "f") == QAWS_STATUS_OK)
	{
		qaws_scalar* fv = (qaws_scalar*)malloc(sizeof(qaws_scalar) * qaws_mesh_3d_get_corner_count(fine));
		qaws_scalar* fu = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * qaws_mesh_3d_get_corner_count(fine));
		unsigned int uid;
		nc = qaws_mesh_3d_get_corner_count(fine);
		qaws_mesh_3d_get_attribute(fv, fine, id);
		qaws_mesh_3d_find_attribute(&uid, fine, "uv");
		qaws_mesh_3d_get_attribute(fu, fine, uid);
		for (i = 0; i < nc; i++)
			if (fu[2 * i] > 0.03 && fu[2 * i] < 0.97 && fu[2 * i + 1] > 0.03 && fu[2 * i + 1] < 0.97)
				e_fine = fmax(e_fine, fabs((double)fv[i] - tx_field(fu[2 * i], fu[2 * i + 1])));
		free(fv); free(fu);
	}
	else
		e_fine = 1;
	sprintf(msg, "a linear field baked into 16 x 16: %u texels covered, off at the centres by %.1e, read between them %.1e, sampled onto a finer mesh %.1e",
		covered, e_tex, e_read, e_fine);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && covered == 256 && e_tex < 1e-12 && e_read < 1e-12 && e_fine < 1e-12, msg);
	qaws_mesh_3d_destroy(m);
	qaws_mesh_3d_destroy(fine);
}

static void test_rules(void)
{
	qaws_mesh_3d* m = tx_grid(1, 1);
	qaws_scalar lab[2] = { 3, 7 }, nrm[18], tl[256], tn[3 * 256], uv[2], y[3];
	qaws_texture t;
	unsigned int i, j, bad = 0, mixed = 0;
	double unit = 0;
	char msg[300];
	tx_add(m, "label", QAWS_MESH_FACE, 1, QAWS_ATTRIBUTE_NEAREST, lab);
	/* (corner normals: the first triangle's toward +x, the second's toward +y, both tilted up) */
	for (i = 0; i < 6; i++)
	{
		double a = i < 3 ? 0 : 1.5707963267948966, l = sqrt(1 + 0.25);
		nrm[3 * i] = (qaws_scalar)(cos(a) / l); nrm[3 * i + 1] = (qaws_scalar)(sin(a) / l); nrm[3 * i + 2] = (qaws_scalar)(0.5 / l);
	}
	nrm[0] = 0; nrm[1] = 0; nrm[2] = 1;   /* (one corner straight up: values differ across each triangle) */
	tx_add(m, "n", QAWS_MESH_CORNER, 3, QAWS_ATTRIBUTE_NORMAL, nrm);
	memset(&t, 0, sizeof t);
	t.texels = tl; t.width = 16; t.height = 16; t.components = 1;
	qaws_mesh_3d_bake_texture(NULL, &t, m, "label", NULL);
	for (j = 0; j < 16; j++)
		for (i = 0; i < 16; i++)
		{
			/* (the first triangle below the diagonal: x > y) */
			double want = (i + 0.5) > (j + 0.5) ? 3 : (i + 0.5) < (j + 0.5) ? 7 : -1;
			if (want >= 0 && (double)tl[j * 16 + i] != want)
				bad++;
		}
	for (i = 0; i < 40; i++)
	{
		uv[0] = (qaws_scalar)fmod(i * 0.6180339887, 1.0);
		uv[1] = (qaws_scalar)fmod(i * 0.7548776662, 1.0);
		qaws_texture_sample(y, &t, uv, 1);
		mixed += (double)y[0] != 3 && (double)y[0] != 7;
	}
	t.texels = tn; t.components = 3;
	qaws_mesh_3d_bake_texture(NULL, &t, m, "n", NULL);
	for (i = 0; i < 256; i++)
		unit = fmax(unit, fabs(sqrt((double)tn[3 * i] * tn[3 * i] + (double)tn[3 * i + 1] * tn[3 * i + 1] + (double)tn[3 * i + 2] * tn[3 * i + 2]) - 1));
	for (i = 0; i < 40; i++)
	{
		uv[0] = (qaws_scalar)fmod(i * 0.6180339887, 1.0);
		uv[1] = (qaws_scalar)fmod(i * 0.7548776662, 1.0);
		qaws_texture_sample(y, &t, uv, 1);
		unit = fmax(unit, fabs(sqrt((double)y[0] * y[0] + (double)y[1] * y[1] + (double)y[2] * y[2]) - 1));
	}
	sprintf(msg, "rules: %u texels off their triangle's label, %u reads a mix of labels; normals off unit length by %.1e (rule %d)", bad, mixed, unit, (int)t.rule);
	printf("  %s\n", msg);
	TEST_ASSERT(bad == 0 && mixed == 0 && unit < 1e-12 && t.rule == QAWS_ATTRIBUTE_NORMAL, msg);
	qaws_mesh_3d_destroy(m);
}

static void test_gutter_repeat_ids(void)
{
	qaws_mesh_3d* m = tx_grid(1, 0.5);   /* (uv the lower left quarter: u and v in [0, 0.5]) */
	qaws_scalar f[4], tex[256], sentinel = -100, ids[2] = { 0, 1 };
	unsigned char cov[256];
	qaws_texture t;
	unsigned int i, j, in = 0, gut = 0, left = 0, untouched = 0, out_of_range = 0, wrap_left = 0, wrap_right = 0, idn = 0;
	char msg[400];
	qaws_status s1, s2, s3;
	for (i = 0; i < 4; i++)
	{
		qaws_vec3 p = qaws_mesh_3d_get_vertex(m, i);
		f[i] = (qaws_scalar)tx_field(p.x, p.y);
	}
	tx_add(m, "f", QAWS_MESH_VERTEX, 1, QAWS_ATTRIBUTE_LINEAR, f);
	memset(&t, 0, sizeof t);
	t.texels = tex; t.width = 16; t.height = 16; t.components = 1;
	{
		qaws_texture_desc d;
		memset(&d, 0, sizeof d);
		d.gutter = 2;
		for (i = 0; i < 256; i++)
			tex[i] = sentinel;
		qaws_mesh_3d_bake_texture(cov, &t, m, "f", &d);
	}
	for (j = 0; j < 16; j++)
		for (i = 0; i < 16; i++)
		{
			unsigned char c = cov[j * 16 + i];
			double v = (double)tex[j * 16 + i];
			in += c == 1 && i < 8;
			gut += c == 2;
			left += c == 2 && (i == 8 || i == 9 || j == 8 || j == 9);
			untouched += c == 0 && v == sentinel;
			out_of_range += c != 0 && (v < 1 || v > 6);
		}
	{
		/* repeated along u: one triangle from u 0.9 across the edge to 0.1 (given as 1.1) */
		qaws_mesh_3d* w = tx_grid(1, 1);
		qaws_scalar uvw[6] = { 0.9f, 0.2f, 0.1f, 0.2f, 0.1f, 0.8f }, two[6 * 2], one[4] = { 1, 1, 1, 1 };
		unsigned int uid;
		unsigned char cw[256];
		qaws_mesh_3d_find_attribute(&uid, w, "uv");
		qaws_mesh_3d_get_attribute(two, w, uid);
		memcpy(two, uvw, sizeof uvw);
		two[6] = two[8] = two[10] = 0.5f; two[7] = two[9] = two[11] = 0.5f;   /* (the other triangle a point: nothing) */
		qaws_mesh_3d_set_attribute(w, uid, two);
		tx_add(w, "one", QAWS_MESH_VERTEX, 1, QAWS_ATTRIBUTE_LINEAR, one);
		t.repeat_u = 1;
		s1 = qaws_mesh_3d_bake_texture(cw, &t, w, "one", NULL);
		for (j = 0; j < 16; j++)
			for (i = 0; i < 16; i++)
				if (cw[j * 16 + i])
				{
					wrap_left += i < 2;
					wrap_right += i > 13;
				}
		t.repeat_u = 0;
		/* (only the faces of id 1: the second triangle, above the diagonal) */
		{
			qaws_mesh_3d* g = tx_grid(1, 1);
			qaws_texture_desc d;
			unsigned char cg[256];
			tx_add(g, "one", QAWS_MESH_VERTEX, 1, QAWS_ATTRIBUTE_LINEAR, one);
			tx_add(g, "texture", QAWS_MESH_FACE, 1, QAWS_ATTRIBUTE_NEAREST, ids);
			memset(&d, 0, sizeof d);
			d.texture_id = "texture";
			d.id = 1;
			qaws_mesh_3d_bake_texture(cg, &t, g, "one", &d);
			for (j = 0; j < 16; j++)
				for (i = 0; i < 16; i++)
					idn += cg[j * 16 + i] && i > j;   /* (a texel of the other face) */
			qaws_mesh_3d_destroy(g);
		}
		/* refusals: three components into one, no uv, an edge attribute */
		{
			qaws_scalar e3[12] = { 0 }, ed[5] = { 0 };
			tx_add(w, "v3", QAWS_MESH_VERTEX, 3, QAWS_ATTRIBUTE_LINEAR, e3);
			tx_add(w, "ed", QAWS_MESH_EDGE, 1, QAWS_ATTRIBUTE_LINEAR, ed);
			s2 = qaws_mesh_3d_bake_texture(NULL, &t, w, "v3", NULL);
			s3 = qaws_mesh_3d_bake_texture(NULL, &t, w, "ed", NULL);
		}
		qaws_mesh_3d_destroy(w);
	}
	{
		qaws_texture_desc d;
		memset(&d, 0, sizeof d);
		d.uv = "none";
		s3 = s3 == QAWS_STATUS_INVALID_ARGUMENT ? qaws_mesh_3d_bake_texture(NULL, &t, m, "f", &d) : QAWS_STATUS_OK;
	}
	sprintf(msg, "a chart on the lower left quarter: %u texels in it, %u in a gutter of 2 (%u in the two rows and columns past it), %u beyond untouched, %u off the "
		"chart's range; repeated: %u texels left of u 1/8 and %u right of 7/8; one id: %u texels of the other face; refusals %s %s",
		in, gut, left, untouched, out_of_range, wrap_left, wrap_right, idn, qaws_status_to_string(s2), qaws_status_to_string(s3));
	printf("  %s\n", msg);
	TEST_ASSERT(in == 64 && left == 36 && gut == 36 && untouched == 156 && out_of_range == 0 && s1 == QAWS_STATUS_OK && wrap_left > 0 && wrap_right > 0
		&& idn == 0 && s2 == QAWS_STATUS_INVALID_ARGUMENT && s3 == QAWS_STATUS_INVALID_ARGUMENT, msg);
	qaws_mesh_3d_destroy(m);
}

/* the fine grid's (n x n, tx_grid's triangles) linear interpolation of per-vertex values f at (x, y) */
static double tx_interp(qaws_scalar const* f, unsigned int n, double x, double y)
{
	unsigned int i = (unsigned int)fmin(n - 1, floor(x * n)), j = (unsigned int)fmin(n - 1, floor(y * n)), a = j * (n + 1) + i;
	double fx = x * n - i, fy = y * n - j;
	/* (below the cell's diagonal: a, a + 1, a + n + 2; above: a, a + n + 2, a + n + 1) */
	if (fx >= fy)
		return (double)f[a] + fx * ((double)f[a + 1] - f[a]) + fy * ((double)f[a + n + 2] - f[a + 1]);
	return (double)f[a] + fy * ((double)f[a + n + 1] - f[a]) + fx * ((double)f[a + n + 2] - f[a + n + 1]);
}

static double tx_wave(double x, double y)
{
	return sin(10 * x) * cos(8 * y);
}

/* a detailed mesh's values onto a coarse one's texture: a 32 x 32 grid with a wave per vertex, the unit square of two
   triangles (its uv x and y) baked from it into 32 x 32 texels: each texel the fine mesh's own interpolation at its
   centre; the coarse mesh's four vertex values off the wave by far more; lifted 0.01 off the source, a max_distance of
   0.005 leaves every texel untouched, 0.02 writes them all */
static void test_transfer(void)
{
	enum { N = 32, R = 32 };
	qaws_mesh_3d *fine = tx_grid(N, 1), *coarse = tx_grid(1, 1);
	qaws_scalar* f = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (N + 1) * (N + 1));
	qaws_scalar tex[R * R];
	unsigned char cov[R * R];
	qaws_texture t;
	qaws_texture_desc d;
	unsigned int i, j, written = 0, near0 = 0, near1 = 0;
	double e_tex = 0, e_vert = 0;
	qaws_status st;
	char msg[400];
	for (i = 0; i < (N + 1) * (N + 1); i++)
	{
		qaws_vec3 p = qaws_mesh_3d_get_vertex(fine, i);
		f[i] = (qaws_scalar)tx_wave(p.x, p.y);
	}
	tx_add(fine, "f", QAWS_MESH_VERTEX, 1, QAWS_ATTRIBUTE_LINEAR, f);
	memset(&t, 0, sizeof t);
	t.texels = tex; t.width = R; t.height = R; t.components = 1;
	st = qaws_mesh_3d_bake_texture_from(cov, &t, coarse, fine, "f", NULL);
	for (j = 0; j < R; j++)
		for (i = 0; i < R; i++)
		{
			double x = (i + 0.5) / R, y = (j + 0.5) / R;
			written += cov[j * R + i] == 1;
			e_tex = fmax(e_tex, fabs((double)tex[j * R + i] - tx_interp(f, N, x, y)));
			/* (what the coarse mesh's vertices can hold: the wave at its four corners, interpolated) */
			{
				double a = tx_wave(0, 0), b = tx_wave(1, 0), c = tx_wave(1, 1), e = tx_wave(0, 1);
				double lin = x >= y ? a + x * (b - a) + y * (c - b) : a + y * (e - a) + x * (c - e);
				e_vert = fmax(e_vert, fabs(lin - tx_interp(f, N, x, y)));
			}
		}
	{
		/* lifted 0.01 above the source */
		unsigned int nv = qaws_mesh_3d_get_vertex_count(coarse);
		qaws_scalar x[12];
		qaws_mesh_3d_get_positions(x, coarse);
		for (i = 0; i < nv; i++)
			x[3 * i + 2] = (qaws_scalar)0.01;
		qaws_mesh_3d_set_positions(coarse, x);
		memset(&d, 0, sizeof d);
		d.max_distance = (qaws_scalar)0.005;
		qaws_mesh_3d_bake_texture_from(cov, &t, coarse, fine, "f", &d);
		for (i = 0; i < R * R; i++)
			near0 += cov[i] != 0;
		d.max_distance = (qaws_scalar)0.02;
		qaws_mesh_3d_bake_texture_from(cov, &t, coarse, fine, "f", &d);
		for (i = 0; i < R * R; i++)
			near1 += cov[i] == 1;
	}
	sprintf(msg, "from a 32 x 32 grid onto two triangles' 32 x 32 texture: %u texels, off the fine mesh by %.1e (its four vertices could only hold it "
		"within %.2f); lifted 0.01: %u texels within 0.005, %u within 0.02", written, e_tex, e_vert, near0, near1);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && written == R * R && e_tex < 1e-12 && e_vert > 0.5 && near0 == 0 && near1 == R * R, msg);
	free(f);
	qaws_mesh_3d_destroy(fine);
	qaws_mesh_3d_destroy(coarse);
}

/* a box [x0, x1] x [y0, y1] x [z0, z1] of 12 triangles, "uv" per corner (0.8 x + 0.1, 0.8 y + 0.1), "texture" per face id */
static qaws_mesh_3d* tx_box(double x0, double x1, double y0, double y1, double z0, double z1, double id)
{
	static unsigned int const f[36] = { 0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 2, 3, 7, 2, 7, 6, 1, 2, 6, 1, 6, 5, 0, 4, 7, 0, 7, 3 };
	qaws_scalar x[24], uv[72], ids[12];
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	unsigned int i;
	for (i = 0; i < 8; i++)
	{
		x[3 * i] = (qaws_scalar)((i == 1 || i == 2 || i == 5 || i == 6) ? x1 : x0);
		x[3 * i + 1] = (qaws_scalar)((i == 2 || i == 3 || i == 6 || i == 7) ? y1 : y0);
		x[3 * i + 2] = (qaws_scalar)(i >= 4 ? z1 : z0);
	}
	memset(&d, 0, sizeof d);
	d.positions = x; d.vertex_count = 8; d.face_vertices = f; d.face_count = 12;
	qaws_mesh_3d_create(&m, &d);
	for (i = 0; i < 36; i++)
	{
		uv[2 * i] = (qaws_scalar)(0.8 * (double)x[3 * f[i]] + 0.1);
		uv[2 * i + 1] = (qaws_scalar)(0.8 * (double)x[3 * f[i] + 1] + 0.1);
	}
	for (i = 0; i < 12; i++)
		ids[i] = (qaws_scalar)id;
	tx_add(m, "uv", QAWS_MESH_CORNER, 2, QAWS_ATTRIBUTE_LINEAR, uv);
	tx_add(m, "texture", QAWS_MESH_FACE, 1, QAWS_ATTRIBUTE_NEAREST, ids);
	return m;
}

/* texture ids through a Boolean: box A (id 0) with a texture linear in u, box B (id 1) with a constant one, united
   (qaws_mesh_3d_boolean carries "uv" per corner and "texture" per face); each texture read onto its own faces of the
   result: A's corners the linear field at their uv (0.8 x + 0.1: exact, cut corners too), B's the constant */
static void test_ids_boolean(void)
{
	qaws_mesh_3d *a = tx_box(0, 1, 0, 1, 0, 1, 0), *b = tx_box(0.5, 1.5, 0.25, 0.75, 0.25, 0.75, 1), *o = NULL;
	qaws_scalar ta[16 * 16], tb[16] = { 0 };
	qaws_texture xa, xb;
	qaws_texture_desc d;
	unsigned int i, j, nc = 0, from_a = 0, from_b = 0;
	double e_a = 0, e_b = 0;
	qaws_status st;
	char msg[300];
	for (j = 0; j < 16; j++)
		for (i = 0; i < 16; i++)
			ta[j * 16 + i] = (qaws_scalar)((i + 0.5) / 16);
	for (i = 0; i < 16; i++)
		tb[i] = 2;
	memset(&xa, 0, sizeof xa);
	xa.texels = ta; xa.width = 16; xa.height = 16; xa.components = 1;
	xb = xa;
	xb.texels = tb; xb.width = 4; xb.height = 4;
	st = qaws_mesh_3d_boolean(&o, NULL, a, b, QAWS_BOOLEAN_UNION, NULL);
	memset(&d, 0, sizeof d);
	d.texture_id = "texture";
	if (st == QAWS_STATUS_OK && o)
	{
		d.id = 0;
		st = qaws_mesh_3d_sample_texture(o, &xa, "value", &d);
		d.id = 1;
		if (st == QAWS_STATUS_OK)
			st = qaws_mesh_3d_sample_texture(o, &xb, "value", &d);
	}
	if (st == QAWS_STATUS_OK && o)
	{
		unsigned int nf = qaws_mesh_3d_get_face_count(o), vid, iid, uid, f, c = 0;
		qaws_scalar *val, *ids, *uv;
		unsigned int *fv, *fs;
		nc = qaws_mesh_3d_get_corner_count(o);
		val = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nc);
		uv = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * nc);
		ids = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nf);
		fv = (unsigned int*)malloc(sizeof(unsigned int) * nc);
		fs = (unsigned int*)malloc(sizeof(unsigned int) * nf);
		qaws_mesh_3d_find_attribute(&vid, o, "value");
		qaws_mesh_3d_find_attribute(&iid, o, "texture");
		qaws_mesh_3d_find_attribute(&uid, o, "uv");
		qaws_mesh_3d_get_attribute(val, o, vid);
		qaws_mesh_3d_get_attribute(ids, o, iid);
		qaws_mesh_3d_get_attribute(uv, o, uid);
		qaws_mesh_3d_get_faces(fv, fs, o);
		for (f = 0; f < nf; f++)
			for (j = 0; j < fs[f]; j++, c++)
			{
				qaws_vec3 p = qaws_mesh_3d_get_vertex(o, fv[c]);
				if ((double)ids[f] == 0)
				{
					from_a++;
					e_a = fmax(e_a, fabs((double)val[c] - (0.8 * (double)p.x + 0.1)) + fabs((double)uv[2 * c] - (0.8 * (double)p.x + 0.1)));
				}
				else
				{
					from_b++;
					e_b = fmax(e_b, fabs((double)val[c] - 2));
				}
			}
		free(val); free(uv); free(ids); free(fv); free(fs);
	}
	sprintf(msg, "texture ids through a union: %u corners of A off its linear texture by %.1e, %u of B off its constant by %.1e", from_a, e_a, from_b, e_b);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && from_a > 0 && from_b > 0 && e_a < 1e-12 && e_b < 1e-12, msg);
	qaws_mesh_3d_destroy(a); qaws_mesh_3d_destroy(b); qaws_mesh_3d_destroy(o);
}

static qaws_curve* tx_line(double x0, double y0, double x1, double y1)
{
	qaws_bezier_desc bd;
	qaws_scalar cp[4];
	qaws_curve* c = NULL;
	cp[0] = (qaws_scalar)x0; cp[1] = (qaws_scalar)y0; cp[2] = (qaws_scalar)x1; cp[3] = (qaws_scalar)y1;
	memset(&bd, 0, sizeof bd);
	bd.dimension = QAWS_DIMENSION_2D; bd.degree = 1; bd.control_points = cp; bd.control_point_count = 2;
	qaws_curve_create_bezier(&c, &bd);
	return c;
}

static double tx_periodic(double s, double t)
{
	return cos(2 * 3.14159265358979323846 * s) + 3 * t;
}

/* a texture over a surface patch: a cylinder (a revolution, u round it) meshed by PARAMETRIC, a 64 x 64 texture of
   cos(2 pi s) + 3 t over its domain repeated along u, read onto the corners through the footpoints: each the field at its
   vertex's footpoint within the bilinear reading's error (seam vertices too); with a few vertices' footpoints taken
   away, their faces skipped without the surface and, given it, read at their closest points as well */
static void test_patch(void)
{
	enum { R = 64 };
	qaws_curve* prof = tx_line(0.5, 0, 0.5, -3.14159265358979323846);
	qaws_surface* s = NULL;
	qaws_surface_revolution_desc rd;
	qaws_mesh_3d* m = NULL;
	qaws_generate_desc gd;
	qaws_scalar tex[R * R];
	qaws_texture x;
	qaws_texture_desc d;
	qaws_range ur, vr;
	unsigned int i, j, nv = 0, nc = 0, cut = 0, skipped = 0;
	double worst = 0, worst_proj = 0;
	qaws_status st;
	char msg[400];
	memset(&rd, 0, sizeof rd);
	rd.profile = prof;
	rd.axis_direction.z = 1;
	qaws_surface_create_revolution(&s, &rd);
	memset(&gd, 0, sizeof gd);
	gd.method = QAWS_GENERATE_PARAMETRIC;
	gd.target_length = (qaws_scalar)0.3;
	st = s ? qaws_mesh_generate_surface(&m, NULL, s, &gd) : QAWS_STATUS_INVALID_ARGUMENT;
	for (j = 0; j < R; j++)
		for (i = 0; i < R; i++)
			tex[j * R + i] = (qaws_scalar)tx_periodic((i + 0.5) / R, (j + 0.5) / R);
	memset(&x, 0, sizeof x);
	x.texels = tex; x.width = R; x.height = R; x.components = 1; x.repeat_u = 1;
	memset(&d, 0, sizeof d);
	d.footpoints = 1;
	d.patch = 0;
	if (s)
	{
		ur = qaws_surface_get_u_range(s);
		vr = qaws_surface_get_v_range(s);
		d.domain[0] = ur.min_value; d.domain[1] = ur.max_value; d.domain[2] = vr.min_value; d.domain[3] = vr.max_value;
	}
	if (st == QAWS_STATUS_OK && m)
	{
		qaws_mesh_location* loc;
		qaws_scalar* val;
		unsigned int *fv, *fs, nf = qaws_mesh_3d_get_face_count(m), id, f, c;
		nv = qaws_mesh_3d_get_vertex_count(m);
		nc = qaws_mesh_3d_get_corner_count(m);
		loc = (qaws_mesh_location*)malloc(sizeof(qaws_mesh_location) * nv);
		val = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nc);
		fv = (unsigned int*)malloc(sizeof(unsigned int) * nc);
		fs = (unsigned int*)malloc(sizeof(unsigned int) * nf);
		qaws_mesh_3d_get_footpoints(loc, m);
		qaws_mesh_3d_get_faces(fv, fs, m);
		st = qaws_mesh_3d_sample_texture(m, &x, "value", &d);
		qaws_mesh_3d_find_attribute(&id, m, "value");
		qaws_mesh_3d_get_attribute(val, m, id);
		for (c = 0; c < nc; c++)
		{
			qaws_mesh_location const* l = loc + fv[c];
			double sv = ((double)l->a - (double)d.domain[0]) / ((double)d.domain[1] - (double)d.domain[0]);
			double tv = ((double)l->b - (double)d.domain[2]) / ((double)d.domain[3] - (double)d.domain[2]);
			/* (away from the clamped rows at v's ends) */
			if (tv > 1.0 / R && tv < 1 - 1.0 / R)
				worst = fmax(worst, fabs((double)val[c] - tx_periodic(sv, tv)));
		}
		/* every seventh vertex's footpoint taken away: without the surface its faces kept as they were (the attribute reset to
		   -5 first), given the surface read at its closest point */
		{
			qaws_mesh_location* l2 = (qaws_mesh_location*)malloc(sizeof(qaws_mesh_location) * nv);
			qaws_scalar* minus = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nc);
			memcpy(l2, loc, sizeof(qaws_mesh_location) * nv);
			for (i = 0; i < nv; i += 7, cut++)
				l2[i].patch = QAWS_MESH_NONE;
			qaws_mesh_3d_set_footpoints(m, l2);
			for (c = 0; c < nc; c++)
				minus[c] = -5;
			qaws_mesh_3d_set_attribute(m, id, minus);
			if (st == QAWS_STATUS_OK)
				st = qaws_mesh_3d_sample_texture(m, &x, "value", &d);
			qaws_mesh_3d_get_attribute(val, m, id);
			for (f = 0, c = 0; f < nf; c += fs[f++])
				skipped += (double)val[c] == -5;
			d.surface = s;
			if (st == QAWS_STATUS_OK)
				st = qaws_mesh_3d_sample_texture(m, &x, "value", &d);
			qaws_mesh_3d_get_attribute(val, m, id);
			for (c = 0; c < nc; c++)
			{
				qaws_mesh_location const* l = loc + fv[c];
				double sv = ((double)l->a - (double)d.domain[0]) / ((double)d.domain[1] - (double)d.domain[0]);
				double tv = ((double)l->b - (double)d.domain[2]) / ((double)d.domain[3] - (double)d.domain[2]);
				if ((double)val[c] == -5)
					worst_proj = 1e300;
				else if (tv > 1.0 / R && tv < 1 - 1.0 / R)
					worst_proj = fmax(worst_proj, fabs((double)val[c] - tx_periodic(sv, tv)));
			}
			free(l2); free(minus);
		}
		free(loc); free(val); free(fv); free(fs);
	}
	sprintf(msg, "a cylinder's patch texture (64 x 64, repeated round it) read through %u vertices' footpoints: off cos(2 pi s) + 3 t by %.1e; %u "
		"footpoints taken away: %u faces skipped, given the surface all read, off by %.1e", nv, worst, cut, skipped, worst_proj);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && nv > 0 && worst < 3e-3 && skipped > 0 && worst_proj < 3e-3, msg);
	qaws_mesh_3d_destroy(m);
	qaws_surface_destroy(s);
	qaws_curve_destroy(prof);
}

int test_187_mesh_texture_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 187: texture attributes\n");
	test_linear();
	test_rules();
	test_gutter_repeat_ids();
	test_transfer();
	test_ids_boolean();
	test_patch();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
