/* ================================================================== */
/*  16. Surfaces: stratified points warped by area, density, curvature */
/* ================================================================== */

#define SW_GRID 18
#define SW_N (SW_GRID * SW_GRID)

/* rho = 0.2 + 4 exp(-|x - c|^2 / (2 s^2)) around a point of the patch */
static qaws_scalar sw_spot(qaws_vec3 p, void* user, qaws_vec3* g, qaws_scalar* H)
{
	double cx = 2.1, cy = 1.0, s2 = 0.45 * 0.45, dx = p.x - cx, dy = p.y - cy;
	double e = 4.0 * exp(-(dx * dx + dy * dy) / (2 * s2));
	(void)user;
	g->x = (qaws_scalar)(-e * dx / s2);
	g->y = (qaws_scalar)(-e * dy / s2);
	g->z = 0;
	H[0] = (qaws_scalar)(e * (dx * dx / s2 - 1) / s2);
	H[1] = (qaws_scalar)(e * dx * dy / (s2 * s2));
	H[2] = 0;
	H[3] = (qaws_scalar)(e * (dy * dy / s2 - 1) / s2);
	H[4] = 0;
	H[5] = 0;
	return (qaws_scalar)(0.2 + e);
}

static void sw_draw_patch(svg* s, projection const* pr, qaws_surface const* surf)
{
	int i, k;
	for (i = 0; i <= 10; i++)
	{
		double xy[2 * 41];
		int dirn;
		for (dirn = 0; dirn < 2; dirn++)
		{
			for (k = 0; k <= 40; k++)
			{
				qaws_surface_jet j;
				qaws_scalar u = (qaws_scalar)(dirn ? i / 10.0 : k / 40.0), v = (qaws_scalar)(dirn ? k / 40.0 : i / 10.0);
				qaws_surface_eval_jet(surf, u, v, QAWS_SJET_P, &j);
				project(pr, j.d[0].x, j.d[0].y, j.d[0].z, &xy[2 * k], &xy[2 * k + 1]);
			}
			svg_polyline(s, xy, 41, "#8c959f", 0.8, 0.8, 0);
		}
	}
}

static void demo_surface_cdf(void)
{
	/* columns crowded toward x = 0 in parameter space: uniform (u, v) bunches up there */
	static double const colx[4] = { 0.0, 0.15, 0.6, 3.0 };
	qaws_scalar cps[48], xi[2 * SW_N];
	qaws_sample_measure_desc dens = { QAWS_MEASURE_DENSITY, 0, sw_spot, NULL };
	qaws_sample_measure_desc curv = { QAWS_MEASURE_CURVATURE, (qaws_scalar)0.3, NULL, NULL };
	qaws_sample_measure_desc const* ms[4] = { NULL, NULL, &dens, &curv };
	char const* cols[4] = { "#cf222e", "#1a7f37", "#8250df", "#bc4c00" };
	qaws_surface_bezier_desc d;
	qaws_surface* surf = NULL;
	qaws_surface_cdf_sample out[SW_N];
	unsigned int rng = 2026u;
	char buf[160];
	int i, j, panel;
	svg s;
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
	/* jittered stratified points */
	for (i = 0; i < SW_GRID; i++)
		for (j = 0; j < SW_GRID; j++)
		{
			rng = rng * 1664525u + 1013904223u;
			xi[2 * (i * SW_GRID + j)] = (qaws_scalar)((i + (rng >> 8) / 16777216.0) / SW_GRID);
			rng = rng * 1664525u + 1013904223u;
			xi[2 * (i * SW_GRID + j) + 1] = (qaws_scalar)((j + (rng >> 8) / 16777216.0) / SW_GRID);
		}
	svg_open(&s, "showcase/16_surface_cdf.svg", 1590, 520, "Surfaces: stratified points warped by inverse CDFs",
		"324 jittered points of the unit square: used as (u, v) directly, then warped by the area, a density field and the curvature sqrt(0.3^2 + k1^2 + k2^2) (marginal and conditional CDFs).");
	for (panel = 0; panel < 4; panel++)
	{
		viewport v = { 30 + panel * 390, 80, 370, 420, 0, 1, 0, 1 };
		projection pr;
		pr.cx = v.x0 + 165;
		pr.cy = v.y0 + 150;
		pr.scale = 82;
		pr.zscale = 1.0;
		svg_panel(&s, &v, panel == 0 ? "(u, v) = xi: crowded where the parameterization is" :
			(panel == 1 ? "area measure: even on the surface" : (panel == 2 ? "density measure: a bright spot" : "curvature measure: gathers where the patch bends")));
		sw_draw_patch(&s, &pr, surf);
		if (panel == 0)
		{
			for (i = 0; i < SW_N; i++)
			{
				qaws_surface_jet jt;
				double sx, sy;
				qaws_surface_eval_jet(surf, xi[2 * i], xi[2 * i + 1], QAWS_SJET_P, &jt);
				project(&pr, jt.d[0].x, jt.d[0].y, jt.d[0].z, &sx, &sy);
				svg_circle(&s, sx, sy, 2.1, "#cf222e", "#cf222e");
			}
		}
		else
		{
			qaws_scalar total = 0;
			qaws_surface_cdf_sample_tangent(NULL, surf, ms[panel], xi, NULL, SW_N, 0, 0, NULL, out, NULL, NULL, &total);
			for (i = 0; i < SW_N; i++)
			{
				double sx, sy;
				project(&pr, out[i].position.x, out[i].position.y, out[i].position.z, &sx, &sy);
				svg_circle(&s, sx, sy, 2.1, cols[panel], cols[panel]);
			}
			sprintf(buf, "total measure %.3f", (double)total);
			svg_text(&s, v.x0 + 12, v.y0 + v.h - 12, 12, "#57606a", "start", buf);
		}
	}
	svg_close(&s);
	qaws_surface_destroy(surf);
	printf("16_surface_cdf: written\n");
}

