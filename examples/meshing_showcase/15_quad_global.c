/*
 * Figure 15: quads by integer-grid maps (GLOBAL) against INSTANT.
 *   - the blob (0.11) and the can from below (0.175, where the flow's
 *     integers give no map and GLOBAL's own quantization does), all quads,
 *     by GLOBAL and by INSTANT side by side: each quad coloured by its smallest
 *     scaled Jacobian (light at 1, amber toward 0.5, red under), irregular
 *     vertices dotted (valence 3 blue, 5 red);
 *   - the blob's map: its dense triangles with the lattice's isolines
 *     (u whole red, v whole blue), the quads' edges before they are placed
 *     and relaxed; the cones where four isolines do not meet.
 * GLOBAL's map through the internal entry (qg_parametrization), the rest
 * through qaws_mesh_quadrangulate.
 */

#include "mesh/qaws_quad_global.h"   /* (the map only) */

static qaws_mesh_3d* ms_quads_by(qaws_mesh_3d const* m, double L, qaws_quad_method method, qaws_mesh_report* rep)
{
	qaws_quadrangulate_desc d;
	qaws_mesh_3d* out = NULL;
	memset(&d, 0, sizeof d);
	d.method = method;
	d.target_length = (qaws_scalar)L;
	d.pure = 1;
	if (qaws_mesh_quadrangulate(m, &d, &out, rep) != QAWS_STATUS_OK)
		printf("  quadrangulate failed\n");
	return out;
}

/* the map's triangles shaded, the segments where u or v is whole */
static void ms_isolines(er_canvas* c, qg_map const* map)
{
	static double const shade[3] = { 236, 236, 232 }, red[3] = { 200, 40, 40 }, blue[3] = { 40, 70, 200 };
	unsigned int t, a, k;
	for (t = 0; t < map->nt; t++)
	{
		double const *x = map->x + 9 * (size_t)t, *u = map->uv + 6 * (size_t)t;
		double e1[3], e2[3], n[3], l;
		for (k = 0; k < 3; k++)
		{
			e1[k] = x[3 + k] - x[k];
			e2[k] = x[6 + k] - x[k];
		}
		n[0] = e1[1] * e2[2] - e1[2] * e2[1];
		n[1] = e1[2] * e2[0] - e1[0] * e2[2];
		n[2] = e1[0] * e2[1] - e1[1] * e2[0];
		l = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
		for (k = 0; k < 3 && l > 0; k++)
			n[k] /= l;
		er_triangle(c, x, x + 3, x + 6, n, n, n, shade, shade, shade);
		for (a = 0; a < 2; a++)
		{
			double lo = fmin(u[a], fmin(u[2 + a], u[4 + a])), hi = fmax(u[a], fmax(u[2 + a], u[4 + a]));
			int i;
			for (i = (int)ceil(lo); i <= (int)floor(hi); i++)
			{
				double p[2][3];
				unsigned int np = 0, e;
				for (e = 0; e < 3 && np < 2; e++)
				{
					double fa = u[2 * e + a] - i, fb = u[2 * ((e + 1) % 3) + a] - i, s;
					if ((fa < 0) == (fb < 0) || fa == fb)
						continue;
					s = fa / (fa - fb);
					for (k = 0; k < 3; k++)
						p[np][k] = x[3 * e + k] + s * (x[3 * ((e + 1) % 3) + k] - x[3 * e + k]);
					np++;
				}
				if (np == 2)
					er_line(c, p[0], p[1], a ? blue : red, 1.4);
			}
		}
	}
}

static void demo_quad_global(void)
{
	er_canvas* c = er_create(2750, 640, 3, 252, 252, 250);
	qaws_mesh_3d *can = ms_can(8, 1, 1.6), *blob = ms_blob(1), *q[4] = { NULL, NULL, NULL, NULL };
	static char const* const title[4] = { "the blob, GLOBAL", "the blob, INSTANT", "the can from below, GLOBAL", "the can, INSTANT" };
	qaws_mesh_report rep;
	ms_qi_stats st;
	char sub[200];
	unsigned int i;
	if (!c || !can || !blob)
		return;
	qaws_mesh_3d_detect_features(can, (qaws_scalar)(50 * MS_PI / 180), 0, 0);
	for (i = 0; i < 4; i++)
	{
		memset(&rep, 0, sizeof rep);
		q[i] = ms_quads_by(i < 2 ? blob : can, i < 2 ? 0.11 : 0.175, (i & 1) ? QAWS_QUADS_INSTANT : QAWS_QUADS_GLOBAL, &rep);
		if (!q[i])
			continue;
		ms_panel(c, 550 * (int)i, 0, 550, 640, q[i], i < 2 ? 0.5 : 0.35, i < 2 ? 0.5 : -0.55, title[i], NULL);
		ms_sj_draw(c, q[i], &st);
		sprintf(sub, "%u quads, %u irr, SJ %.2f, AD %.1f, ED %.2f%s", rep.quality.quads, st.irr, st.smin, st.ad, st.ed,
			(i & 1) ? "" : rep.iterations ? "" : " (fell back)");
		er_text(c, 550 * (int)i + 14, 46, sub, 2, 95, 95, 95);
		printf("  %s: %u quads, irregular %u, SJ min %.3f mean %.3f, AD %.2f deg, ED %.3f%s\n", title[i] + 4, rep.quality.quads, st.irr, st.smin,
			st.smean, st.ad, st.ed, (i & 1) ? "" : rep.iterations ? ", by the map" : ", INSTANT's (fell back)");
	}
	{
		qaws_quadrangulate_desc d;
		qg_map map;
		memset(&d, 0, sizeof d);
		d.method = QAWS_QUADS_GLOBAL;
		d.target_length = (qaws_scalar)0.11;
		d.pure = 1;
		if (q[0] && qg_parametrization(blob, &d, &map) == QAWS_STATUS_OK)
		{
			ms_panel(c, 2200, 0, 550, 640, q[0], 0.5, 0.5, "the blob's map", "isolines: u whole red, v whole blue");
			ms_isolines(c, &map);
			printf("  the blob's map: %u triangles\n", map.nt);
			qg_map_free(&map);
		}
	}
	er_save(c, "showcase/mesh15_quad_global.png");
	printf("mesh15_quad_global.png\n");
	for (i = 0; i < 4; i++)
		qaws_mesh_3d_destroy(q[i]);
	qaws_mesh_3d_destroy(can);
	qaws_mesh_3d_destroy(blob);
	er_destroy(c);
}
