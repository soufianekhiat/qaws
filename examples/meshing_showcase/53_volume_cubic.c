/*
 * Figure 53: cubic volume elements (qaws_volume_mesh_cubic_nodes) beside
 * quadratic ones, a coarse ball of hexahedra (an O-grid of 32 cells) cut in
 * half, each border face drawn finely from its nodes:
 *   - straight;
 *   - triquadratic (27 nodes a cell), the border on the sphere;
 *   - tricubic (64 nodes a cell), the border closer still (the corner cells,
 *     at 0.085 of their straight Jacobians, untangled over 0.1 by moving
 *     their inside nodes, the border kept).
 * Printed: how far the drawn border strays from the unit sphere.
 */

/* the hexahedron's nodes of degree p as the library numbers them: (i, j, k) in [0, p]^3 per slot */
static unsigned int ms_vq_lattice(unsigned int p, unsigned char (*lat)[3])
{
	unsigned char c[8][3];
	unsigned int q, m, a, b, i, j, k, n = 0, d;
	for (q = 0; q < 8; q++)
	{
		c[q][0] = (unsigned char)((q == 1 || q == 2 || q == 5 || q == 6) ? p : 0);
		c[q][1] = (unsigned char)((q == 2 || q == 3 || q == 6 || q == 7) ? p : 0);
		c[q][2] = (unsigned char)(q >= 4 ? p : 0);
		memcpy(lat[n++], c[q], 3);
	}
	for (q = 0; q < 12; q++)
		for (m = 1; m < p; m++, n++)
			for (d = 0; d < 3; d++) lat[n][d] = (unsigned char)((c[ms_vc_he[q][0]][d] * (p - m) + c[ms_vc_he[q][1]][d] * m) / p);
	for (q = 0; q < 6; q++)
		for (b = 1; b < p; b++)
			for (a = 1; a < p; a++, n++)
				for (d = 0; d < 3; d++)
					lat[n][d] = (unsigned char)(((int)p * c[ms_vc_hvf[q][0]][d] + (c[ms_vc_hvf[q][1]][d] - c[ms_vc_hvf[q][0]][d]) * (int)a +
						(c[ms_vc_hvf[q][3]][d] - c[ms_vc_hvf[q][0]][d]) * (int)b) / (int)p);
	for (k = 1; k < p; k++)
		for (j = 1; j < p; j++)
			for (i = 1; i < p; i++, n++)
			{
				lat[n][0] = (unsigned char)i;
				lat[n][1] = (unsigned char)j;
				lat[n][2] = (unsigned char)k;
			}
	return n;
}

/* 1D Lagrange on p + 1 equal steps of [0, 1] at u */
static void ms_vq_lag(unsigned int p, double u, double* v)
{
	unsigned int m, r;
	for (m = 0; m <= p; m++)
	{
		v[m] = 1;
		for (r = 0; r <= p; r++)
			if (r != m)
				v[m] *= (u - (double)r / p) / (((double)m - (double)r) / p);
	}
}

/* the cells with their centroid at y > 0: their faces not shared with another kept cell, drawn n x n from the face's
   (p + 1)^2 nodes; the border's largest distance from the unit sphere into *off */
static void ms_vq_draw(er_canvas* c, qaws_volume_mesh const* v, qaws_scalar const* pos, unsigned int const* en, unsigned int p, double const* rgb,
	double const* cutc, double const* wire, double* off)
{
	unsigned char lat[64][3];
	unsigned int nc = qaws_volume_mesh_get_cell_count(v), nf = qaws_volume_mesh_get_facet_count(v), nn = ms_vq_lattice(p, lat), cc, f, n = 8;
	unsigned int* adj = (unsigned int*)malloc(sizeof(unsigned int) * nf);
	unsigned char* kept = (unsigned char*)malloc(nc);
	qaws_volume_mesh_get_cell_adjacency(adj, v);
	for (cc = 0; cc < nc; cc++)
	{
		double y = 0;
		unsigned int k;
		for (k = 0; k < 8; k++) y += (double)pos[3 * en[nn * cc + k] + 1];
		kept[cc] = (unsigned char)(y > 0);
	}
	for (cc = 0; cc < nc; cc++)
	{
		unsigned int const* e = en + nn * cc;
		if (!kept[cc])
			continue;
		for (f = 0; f < 6; f++)
		{
			/* (facet f has corners ms_vc_hf[f]; the library's face of the same corners gives the frame) */
			unsigned int o = adj[6 * cc + f], s, hit = 0, a, b, i, j, k, z, w;
			double g[9][9][3], va[4], vb[4], nodes[4][4][3];
			if (o != QAWS_MESH_NONE && kept[o])
				continue;
			for (s = 0; s < 6; s++)
			{
				unsigned int all = 1;
				for (z = 0; z < 4; z++)
				{
					int found = 0;
					for (w = 0; w < 4; w++) found |= ms_vc_hvf[s][w] == ms_vc_hf[f][z];
					all &= (unsigned int)found;
				}
				if (all)
					hit = s;
			}
			for (b = 0; b <= p; b++)
				for (a = 0; a <= p; a++)
				{
					unsigned char want[3];
					unsigned int d, q;
					for (d = 0; d < 3; d++)
						want[d] = (unsigned char)(((int)p * lat[ms_vc_hvf[hit][0]][d] + ((int)lat[ms_vc_hvf[hit][1]][d] - lat[ms_vc_hvf[hit][0]][d]) * (int)a +
							((int)lat[ms_vc_hvf[hit][3]][d] - lat[ms_vc_hvf[hit][0]][d]) * (int)b) / (int)p);
					for (q = 0; q < nn; q++)
						if (!memcmp(lat[q], want, 3))
							for (d = 0; d < 3; d++) nodes[b][a][d] = (double)pos[3 * e[q] + d];
				}
			for (j = 0; j <= n; j++)
				for (i = 0; i <= n; i++)
				{
					ms_vq_lag(p, (double)i / n, va);
					ms_vq_lag(p, (double)j / n, vb);
					for (k = 0; k < 3; k++)
					{
						g[j][i][k] = 0;
						for (b = 0; b <= p; b++)
							for (a = 0; a <= p; a++) g[j][i][k] += va[a] * vb[b] * nodes[b][a][k];
					}
					if (o == QAWS_MESH_NONE)
						*off = fmax(*off, fabs(sqrt(g[j][i][0] * g[j][i][0] + g[j][i][1] * g[j][i][1] + g[j][i][2] * g[j][i][2]) - 1));
				}
			for (j = 0; j < n; j++)
				for (i = 0; i < n; i++)
				{
					double const* col = o == QAWS_MESH_NONE ? rgb : cutc;
					er_triangle(c, g[j][i], g[j][i + 1], g[j + 1][i], NULL, NULL, NULL, col, col, col);
					er_triangle(c, g[j][i + 1], g[j + 1][i + 1], g[j + 1][i], NULL, NULL, NULL, col, col, col);
				}
			for (i = 0; i < n; i++)
			{
				er_line(c, g[0][i], g[0][i + 1], wire, 0.7);
				er_line(c, g[n][i], g[n][i + 1], wire, 0.7);
				er_line(c, g[i][0], g[i + 1][0], wire, 0.7);
				er_line(c, g[i][n], g[i + 1][n], wire, 0.7);
			}
		}
	}
	free(adj);
	free(kept);
}

static void demo_volume_cubic(void)
{
	static double const wire[3] = { 60, 66, 80 }, skin[3][3] = { { 200, 205, 215 }, { 130, 190, 140 }, { 120, 150, 210 } }, cutc[3] = { 232, 196, 120 };
	static char const* const titles[3] = { "straight", "triquadratic", "tricubic" };
	er_canvas* c = er_create(1350, 700, 3, 252, 252, 250);
	qaws_mesh_3d* fine = ms_rm_sphere(96, 48);
	qaws_mesh_reference* ref = NULL;
	qaws_volume_mesh* b = ms_vc_ball(2, 1, 0.5);
	unsigned int k;
	if (!c || !fine || !b)
		return;
	qaws_mesh_reference_create(&ref, fine, NULL);
	for (k = 0; k < 3 && ref; k++)
	{
		unsigned int p = k == 2 ? 3 : 2, nn = 0, nen = 0, i;
		qaws_scalar* pos;
		unsigned int* en;
		double* pd;
		double off = 0;
		qaws_volume_curve_report r;
		char sub[200];
		if ((p == 3 ? qaws_volume_mesh_cubic_node_counts(&nn, &nen, b) : qaws_volume_mesh_quadratic_node_counts(&nn, &nen, b)) != QAWS_STATUS_OK)
			continue;
		pos = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nn);
		en = (unsigned int*)malloc(sizeof(unsigned int) * nen);
		pd = (double*)malloc(sizeof(double) * 3 * nn);
		memset(&r, 0, sizeof r);
		if (p == 3)
			qaws_volume_mesh_cubic_nodes(pos, en, &r, b, ref, NULL);
		else
			qaws_volume_mesh_quadratic_nodes(pos, en, &r, b, k ? ref : NULL, NULL);
		for (i = 0; i < 3 * nn; i++) pd[i] = (double)pos[i];
		er_panel(c, 450 * (int)k, 72, 450, 628, pd, nn, 3.8, 0.5, 0.06);
		ms_vq_draw(c, b, pos, en, p, skin[k], cutc, wire, &off);
		sprintf(sub, "off the sphere %.4f", off);
		er_text(c, 450 * (int)k + 14, 12, titles[k], 3, 30, 30, 30);
		er_text(c, 450 * (int)k + 14, 46, sub, 2, 95, 95, 95);
		printf("  %s: %s, Jacobian >= %.2f, %u untangled, %u relaxed\n", titles[k], sub, (double)r.min_jacobian, r.untangled, r.relaxed);
		free(pos);
		free(en);
		free(pd);
	}
	qaws_volume_mesh_destroy(b);
	qaws_mesh_reference_destroy(ref);
	qaws_mesh_3d_destroy(fine);
	er_save(c, "showcase/mesh53_volume_cubic.png");
	printf("mesh53_volume_cubic.png\n");
	er_destroy(c);
}
