/*
 * Figure 47: feature-preserving denoising (QAWS_SMOOTH_BILATERAL).
 *   - a box of triangles, its vertices shaken (its edges not marked);
 *   - Taubin smoothing: the noise goes, the edges round off;
 *   - bilateral: the noise goes, the edges stay sharp;
 *   - the same on a noisy blob.
 */

/* a box's faces as n x n grids of triangle pairs, every vertex shaken along a random direction by amp */
static qaws_mesh_3d* ms_dn_box(unsigned int n, double amp)
{
	unsigned int N = n + 1, *id = (unsigned int*)malloc(sizeof(unsigned int) * N * N * N), nv = 0, nf = 0, a, s, u, v, k, seed = 99;
	double* x = (double*)malloc(sizeof(double) * 3 * N * N * N);
	unsigned int* f = (unsigned int*)malloc(sizeof(unsigned int) * 36 * n * n);
	static double const size[3] = { 1.4, 1, 0.8 };
	qaws_mesh_3d* m;
	for (k = 0; k < N * N * N; k++)
	{
		unsigned int p[3];
		p[0] = k % N;
		p[1] = (k / N) % N;
		p[2] = k / (N * N);
		id[k] = 0xFFFFFFFFu;
		if (p[0] % n && p[1] % n && p[2] % n)
			continue;
		id[k] = nv;
		for (a = 0; a < 3; a++)
			x[3 * nv + a] = size[a] * ((double)p[a] / n - 0.5) + amp * (ms_uniform(&seed) - 0.5);
		nv++;
	}
	for (a = 0; a < 3; a++)
		for (s = 0; s < 2; s++)
			for (u = 0; u < n; u++)
				for (v = 0; v < n; v++)
				{
					unsigned int q[4], c, b = (a + 1) % 3, d = (a + 2) % 3;
					for (c = 0; c < 4; c++)
					{
						unsigned int p[3];
						p[a] = s * n;
						p[b] = u + (c == 1 || c == 2);
						p[d] = v + (c >= 2);
						q[c] = id[p[0] + N * (p[1] + N * p[2])];
					}
					if (!s)
					{
						c = q[1];
						q[1] = q[3];
						q[3] = c;
					}
					f[3 * nf] = q[0]; f[3 * nf + 1] = q[1]; f[3 * nf + 2] = q[2]; nf++;
					f[3 * nf] = q[0]; f[3 * nf + 1] = q[2]; f[3 * nf + 2] = q[3]; nf++;
				}
	m = ms_make(x, nv, f, NULL, nf);
	free(id);
	free(x);
	free(f);
	return m;
}

static void demo_denoise(void)
{
	static double const skin[3] = { 200, 205, 215 }, blue[3] = { 140, 170, 220 }, green[3] = { 150, 200, 160 }, wire[3] = { 90, 95, 110 };
	er_canvas* c = er_create(1800, 640, 3, 252, 252, 250);
	qaws_mesh_3d *noisy = ms_dn_box(16, 0.03), *taubin = ms_dn_box(16, 0.03), *bil = ms_dn_box(16, 0.03), *blob = ms_blob(3);
	qaws_smooth_desc d;
	double const yaw = 0.6, pitch = 0.5;
	if (!c || !noisy || !taubin || !bil || !blob)
		return;
	memset(&d, 0, sizeof d);
	d.method = QAWS_SMOOTH_TAUBIN;
	d.iterations = 30;
	qaws_mesh_3d_smooth(taubin, &d);
	memset(&d, 0, sizeof d);
	d.method = QAWS_SMOOTH_BILATERAL;
	qaws_mesh_3d_smooth(bil, &d);
	ms_panel(c, 0, 0, 450, 640, noisy, yaw, pitch, "noisy", "a box, its edges not marked");
	ms_draw(c, noisy, skin, NULL, 0, wire, 0.3);
	ms_panel(c, 450, 0, 450, 640, taubin, yaw, pitch, "Taubin", "30 steps: the edges round");
	ms_draw(c, taubin, blue, NULL, 0, wire, 0.3);
	ms_panel(c, 900, 0, 450, 640, bil, yaw, pitch, "bilateral", "the edges stay sharp");
	ms_draw(c, bil, green, NULL, 0, wire, 0.3);
	{
		/* the blob shaken, then denoised */
		unsigned int n = qaws_mesh_3d_get_vertex_count(blob), i, seed = 5;
		qaws_scalar* p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n);
		qaws_mesh_3d_get_positions(p, blob);
		for (i = 0; i < 3 * n; i++) p[i] += (qaws_scalar)(0.02 * (ms_uniform(&seed) - 0.5));
		qaws_mesh_3d_set_positions(blob, p);
		free(p);
		qaws_mesh_3d_smooth(blob, &d);
		ms_panel(c, 1350, 0, 450, 640, blob, 0.5, 0.45, "a blob, bilateral", "shaken by 0.02, denoised");
		ms_draw(c, blob, green, NULL, 1, NULL, 0);
	}
	qaws_mesh_3d_destroy(noisy);
	qaws_mesh_3d_destroy(taubin);
	qaws_mesh_3d_destroy(bil);
	qaws_mesh_3d_destroy(blob);
	er_save(c, "showcase/mesh47_denoise.png");
	printf("mesh47_denoise.png\n");
	er_destroy(c);
}
