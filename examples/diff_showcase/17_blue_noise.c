/* ================================================================== */
/*  17. Blue noise on a patch through the xi adjoints of the warp    */
/* ================================================================== */

#define BN_N 220
#define BN_ITERS 70

/* Gaussian repulsion of the warped samples and its gradient on them. */
static double bn_energy(qaws_surface_cdf_sample const* s, double sigma, qaws_surface_cdf_sample* adj)
{
	double e = 0, inv = 1.0 / (2 * sigma * sigma);
	int i, j;
	for (i = 0; i < BN_N; i++)
	{
		adj[i].u = 0;
		adj[i].v = 0;
		adj[i].position = v3(0, 0, 0);
	}
	for (i = 0; i < BN_N; i++)
		for (j = i + 1; j < BN_N; j++)
		{
			double dx = s[i].position.x - s[j].position.x, dy = s[i].position.y - s[j].position.y, dz = s[i].position.z - s[j].position.z;
			double k = exp(-(dx * dx + dy * dy + dz * dz) * inv), c = -2 * inv * k;
			e += k;
			adj[i].position.x += (qaws_scalar)(c * dx); adj[i].position.y += (qaws_scalar)(c * dy); adj[i].position.z += (qaws_scalar)(c * dz);
			adj[j].position.x -= (qaws_scalar)(c * dx); adj[j].position.y -= (qaws_scalar)(c * dy); adj[j].position.z -= (qaws_scalar)(c * dz);
		}
	return e;
}

/* Mean and minimum nearest-neighbor distance on the surface. */
static void bn_spacing(qaws_surface_cdf_sample const* s, double* mean, double* minimum)
{
	int i, j;
	*mean = 0;
	*minimum = 1e30;
	for (i = 0; i < BN_N; i++)
	{
		double best = 1e30;
		for (j = 0; j < BN_N; j++)
			if (j != i)
			{
				double dx = s[i].position.x - s[j].position.x, dy = s[i].position.y - s[j].position.y, dz = s[i].position.z - s[j].position.z;
				double d = sqrt(dx * dx + dy * dy + dz * dz);
				if (d < best) best = d;
			}
		*mean += best / BN_N;
		if (best < *minimum) *minimum = best;
	}
}

static void demo_blue_noise(void)
{
	static double const colx[4] = { 0.0, 0.15, 0.6, 3.0 };
	qaws_scalar cps[48], xi[2 * BN_N], g[2 * BN_N];
	static qaws_surface_cdf_sample s0[BN_N], s[BN_N], adj[BN_N];
	double loss[BN_ITERS + 1], m0, n0, m1, n1, sigma, total;
	qaws_surface_bezier_desc d;
	qaws_surface* surf = NULL;
	qaws_scalar tot = 0;
	unsigned int rng = 99u;
	char buf[200];
	int i, j, it;
	svg sv;
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
		{
			cps[3 * (i * 4 + j)] = (qaws_scalar)colx[j];
			cps[3 * (i * 4 + j) + 1] = (qaws_scalar)(i * 0.7);
			cps[3 * (i * 4 + j) + 2] = (qaws_scalar)(0.5 * sin(1.1 * i + 0.9 * j));
		}
	d.u_degree = 3;
	d.v_degree = 3;
	d.control_points = (qaws_vec3 const*)cps;
	d.u_point_count = 4;
	d.v_point_count = 4;
	qaws_surface_create_bezier(&d, &surf);
	for (i = 0; i < 2 * BN_N; i++)
	{
		rng = rng * 1664525u + 1013904223u;
		xi[i] = (qaws_scalar)((rng >> 8) / 16777216.0);
	}
	qaws_surface_cdf_sample_tangent(NULL, surf, NULL, xi, NULL, BN_N, 4, 4, NULL, s0, NULL, NULL, &tot);
	total = tot;
	sigma = 0.55 * sqrt(total / BN_N);
	for (it = 0; it <= BN_ITERS; it++)
	{
		qaws_surface_cdf_sample_tangent(NULL, surf, NULL, xi, NULL, BN_N, 4, 4, NULL, s, NULL, NULL, NULL);
		loss[it] = bn_energy(s, sigma, adj);
		if (it == BN_ITERS)
			break;
		memset(g, 0, sizeof(g));
		qaws_surface_cdf_sample_adjoint(NULL, surf, NULL, xi, BN_N, 4, 4, adj, NULL, g);
		/* sign-like normalized steps on the points, shrinking over the run */
		{
			int k;
			for (k = 0; k < 2 * BN_N; k++)
			{
				/* soft walls at the square's sides keep points off the clamped edges */
				double x, tau = 0.25 / sqrt((double)BN_N), beta = 0.2 * loss[0] / BN_N;
				g[k] += (qaws_scalar)(beta * (exp(-(1 - xi[k]) / tau) - exp(-xi[k] / tau)) / tau);
				x = xi[k] - 0.012 * (1.0 - 0.8 * it / (double)BN_ITERS) * g[k] / (fabs(g[k]) + 1e-3 * loss[0] / BN_N);
				xi[k] = (qaws_scalar)(x < 0 ? 0 : (x > 1 ? 1 : x));
			}
		}
	}
	bn_spacing(s0, &m0, &n0);
	bn_spacing(s, &m1, &n1);
	svg_open(&sv, "showcase/17_blue_noise.svg", 1200, 560, "Blue noise on a patch through the warp adjoints",
		"220 random points of the unit square warped by area; their Gaussian repulsion on the surface pulled back to the points (xi adjoints) and descended.");
	for (i = 0; i < 2; i++)
	{
		viewport v = { 30 + i * 390, 80, 370, 440, 0, 1, 0, 1 };
		projection pr;
		qaws_surface_cdf_sample const* ss = i ? s : s0;
		pr.cx = v.x0 + 165;
		pr.cy = v.y0 + 150;
		pr.scale = 82;
		pr.zscale = 1.0;
		svg_panel(&sv, &v, i ? "after 70 steps: blue noise on the surface" : "white noise warped by area");
		sw_draw_patch(&sv, &pr, surf);
		for (j = 0; j < BN_N; j++)
		{
			double sx, sy;
			project(&pr, ss[j].position.x, ss[j].position.y, ss[j].position.z, &sx, &sy);
			svg_circle(&sv, sx, sy, 2.4, i ? "#0969da" : "#cf222e", i ? "#0969da" : "#cf222e");
		}
		sprintf(buf, "nearest neighbor: mean %.3f, min %.3f", i ? m1 : m0, i ? n1 : n0);
		svg_text(&sv, v.x0 + 12, v.y0 + v.h - 12, 12, "#57606a", "start", buf);
	}
	{
		viewport lp = { 810, 80, 360, 440, 0, 1, 0, 1 };
		svg_loss_plot(&sv, &lp, loss, BN_ITERS + 1, "#0969da", "repulsion energy (log10) per step");
	}
	svg_close(&sv);
	qaws_surface_destroy(surf);
	printf("17_blue_noise: energy %.4g -> %.4g, nearest neighbor mean %.3f -> %.3f, min %.3f -> %.3f\n", loss[0], loss[BN_ITERS], m0,
		m1, n0, n1);
}

