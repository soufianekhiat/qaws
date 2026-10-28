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
 *     lifted 0.01 off, max_distance 0.005 writes none, 0.02 all.
 */

#include "test_common.h"
#include "qaws_mesh_texture.h"
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

int test_187_mesh_texture_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 187: texture attributes\n");
	test_linear();
	test_rules();
	test_gutter_repeat_ids();
	test_transfer();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
