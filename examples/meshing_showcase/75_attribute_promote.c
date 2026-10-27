/*
 * Figure 75: attribute promotion (qaws_mesh_promote.h) on a ball.
 *   - a colour per face: six regions, by the axis its centre leans to;
 *   - promoted to the vertices by AVERAGE (Oklab): the borders blend;
 *   - by MODE: each vertex the most common colour about it, the borders kept;
 *   - pieces: each face's x + 0.6 y + 0.3 z averaged over its region (a
 *     piece attribute on the faces) and shown on a ramp: one value a region.
 * Printed: faces, the regions' heights.
 */

/* per corner colours (0 to 255) of a colour attribute of any element but edges */
static double* ms_pm_colours(qaws_mesh_3d const* m, char const* name)
{
	unsigned int nt = qaws_mesh_3d_get_triangle_count(m), i, k, el, t[3];
	double *v = ms_at_get(m, name, &el), *o;
	if (!v)
		return NULL;
	o = (double*)malloc(sizeof(double) * 9 * (nt + 1));
	for (i = 0; i < nt; i++)
	{
		qaws_mesh_3d_get_triangle(t, m, i);
		for (k = 0; k < 9; k++)
		{
			unsigned int c = k / 3, e = el == QAWS_MESH_FACE ? i : el == QAWS_MESH_VERTEX ? t[c] : 3 * i + c;
			o[9 * i + k] = 255 * fmin(1, fmax(0, v[3 * e + k % 3]));
		}
	}
	free(v);
	return o;
}

static void demo_attribute_promote(void)
{
	static double const hue[6][3] = { { 0.9, 0.3, 0.3 }, { 0.3, 0.75, 0.4 }, { 0.3, 0.45, 0.95 }, { 0.95, 0.8, 0.25 }, { 0.7, 0.35, 0.85 }, { 0.3, 0.8, 0.85 } };
	er_canvas* c = er_create(2400, 640, 3, 252, 252, 250);
	qaws_mesh_3d* ball = ms_at_ball(1, 0, 0, 0);
	unsigned int nf, i, k, t[3];
	qaws_scalar *col, *reg, *hgt;
	qaws_attribute_desc ad;
	qaws_attribute_promote_desc pd;
	double* cc;
	char sub[200];
	if (!c || !ball)
		return;
	nf = qaws_mesh_3d_get_triangle_count(ball);
	col = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nf);
	reg = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nf);
	hgt = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nf);
	for (i = 0; i < nf; i++)
	{
		double g[3] = { 0, 0, 0 };
		unsigned int ax = 0, r;
		qaws_mesh_3d_get_triangle(t, ball, i);
		for (k = 0; k < 3; k++)
		{
			qaws_vec3 p = qaws_mesh_3d_get_vertex(ball, t[k]);
			g[0] += p.x / 3; g[1] += p.y / 3; g[2] += p.z / 3;
		}
		for (k = 1; k < 3; k++)
			if (fabs(g[k]) > fabs(g[ax]))
				ax = k;
		r = 2 * ax + (g[ax] < 0);
		reg[i] = (qaws_scalar)r;
		hgt[i] = (qaws_scalar)(g[0] + 0.6 * g[1] + 0.3 * g[2]);
		for (k = 0; k < 3; k++)
			col[3 * i + k] = (qaws_scalar)hue[r][k];
	}
	memset(&ad, 0, sizeof ad);
	ad.name = "paint"; ad.element = QAWS_MESH_FACE; ad.components = 3; ad.rule = QAWS_ATTRIBUTE_OKLAB;
	qaws_mesh_3d_add_attribute(&k, ball, col, &ad);
	ad.name = "region"; ad.components = 1; ad.rule = QAWS_ATTRIBUTE_NEAREST;
	qaws_mesh_3d_add_attribute(&k, ball, reg, &ad);
	ad.name = "height"; ad.rule = QAWS_ATTRIBUTE_LINEAR;
	qaws_mesh_3d_add_attribute(&k, ball, hgt, &ad);
	if ((cc = ms_pm_colours(ball, "paint")) != NULL)
	{
		sprintf(sub, "a colour per face: %u faces, six regions", nf);
		ms_panel(c, 0, 0, 600, 640, ball, -0.6, 0.45, "per face", sub);
		ms_at_draw(c, ball, cc, NULL, NULL, NULL, 0);
		free(cc);
	}
	memset(&pd, 0, sizeof pd);
	pd.name = "paint"; pd.element = QAWS_MESH_VERTEX; pd.method = QAWS_PROMOTE_AVERAGE; pd.new_name = "avg";
	if (qaws_mesh_3d_promote_attribute(ball, &pd) == QAWS_STATUS_OK && (cc = ms_pm_colours(ball, "avg")) != NULL)
	{
		ms_panel(c, 600, 0, 600, 640, ball, -0.6, 0.45, "to vertices: AVERAGE", "each vertex its faces' mean (Oklab)");
		ms_at_draw(c, ball, cc, NULL, NULL, NULL, 0);
		free(cc);
	}
	pd.method = QAWS_PROMOTE_MODE; pd.new_name = "mode";
	if (qaws_mesh_3d_promote_attribute(ball, &pd) == QAWS_STATUS_OK && (cc = ms_pm_colours(ball, "mode")) != NULL)
	{
		ms_panel(c, 1200, 0, 600, 640, ball, -0.6, 0.45, "to vertices: MODE", "each vertex its faces' commonest");
		ms_at_draw(c, ball, cc, NULL, NULL, NULL, 0);
		free(cc);
	}
	memset(&pd, 0, sizeof pd);
	pd.name = "height"; pd.element = QAWS_MESH_FACE; pd.method = QAWS_PROMOTE_AVERAGE; pd.new_name = "region_height"; pd.piece = "region";
	if (qaws_mesh_3d_promote_attribute(ball, &pd) == QAWS_STATUS_OK)
	{
		unsigned int el;
		double* rh = ms_at_get(ball, "region_height", &el);
		double* rc = (double*)malloc(sizeof(double) * 9 * nf);
		for (i = 0; rh && i < nf; i++)
			for (k = 0; k < 3; k++)
				ms_sg_ramp(0.5 + 0.45 * rh[i], rc + 9 * i + 3 * k);
		ms_panel(c, 1800, 0, 600, 640, ball, -0.6, 0.45, "pieces: a region's mean", "x + 0.6 y + 0.3 z averaged per region");
		if (rh)
			ms_at_draw(c, ball, rc, NULL, NULL, NULL, 0);
		for (k = 0; rh && k < 6; k++)
			for (i = 0; i < nf; i++)
				if ((unsigned int)reg[i] == k)
				{
					printf("  region %u: height %.4f\n", k, rh[i]);
					break;
				}
		free(rh); free(rc);
	}
	er_save(c, "showcase/mesh75_attribute_promote.png");
	printf("mesh75_attribute_promote.png\n");
	free(col); free(reg); free(hgt);
	qaws_mesh_3d_destroy(ball);
	er_destroy(c);
}
