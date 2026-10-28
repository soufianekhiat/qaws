/*
 * Figure 81: texture attributes (qaws_mesh_texture.h): a colour held on
 * texels rather than vertices.
 *   - a detailed ball (8064 faces), a checker of colours per vertex (sRGB);
 *   - decimated to 400 faces, the colours carried per vertex: the checker
 *     smeared, its vertices too few for it;
 *   - the decimated ball's charts (auto seams, packed) as a 256 x 256
 *     texture baked from the detailed ball (qaws_mesh_3d_bake_texture_from,
 *     a gutter of 4);
 *   - the decimated ball drawn with that texture (each triangle cut in 12 x
 *     12 for the drawing, each piece's corners read from the texture): the
 *     checker back on 400 faces.
 * Printed: faces, charts, texels written, and the share of the decimated
 * ball's points whose colour is off the checker there (by more than 0.25),
 * from the texture and from the vertices.
 */

/* the checker: blue and orange cells, 12 round z by 8 down it */
static void ms_tx_checker(double const* p, double* rgb)
{
	double ph = atan2(p[1], p[0]), th = acos(fmax(-1, fmin(1, p[2] / sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]))));
	int a = (int)floor(ph / (2 * MS_PI) * 12 + 12), b = (int)floor(th / MS_PI * 8);
	static double const c0[3] = { 0.20, 0.35, 0.80 }, c1[3] = { 0.95, 0.60, 0.15 };
	double const* c = (a + b) % 2 ? c1 : c0;
	rgb[0] = c[0]; rgb[1] = c[1]; rgb[2] = c[2];
}

/* per corner colours (9 per triangle, 0..255) of a vertex colour attribute */
static double* ms_tx_corner_rgb(qaws_mesh_3d const* m, char const* name)
{
	unsigned int el, nt = qaws_mesh_3d_get_triangle_count(m), i, k;
	double* v = ms_at_get(m, name, &el);
	double* o = (double*)malloc(sizeof(double) * 9 * (nt + 1));
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (nt + 1));
	qaws_mesh_3d_get_triangles(t, m);
	for (i = 0; v && i < 3 * nt; i++)
		for (k = 0; k < 3; k++)
			o[3 * i + k] = 255 * v[3 * (el == QAWS_MESH_VERTEX ? t[i] : i) + k];
	free(v); free(t);
	return o;
}

static void demo_textures(void)
{
	static double const wire[3] = { 60, 66, 80 };
	er_canvas* c = er_create(2400, 640, 3, 252, 252, 250);
	qaws_mesh_3d *hi = ms_rm_sphere(96, 43), *lo = NULL;
	unsigned int nv = hi ? qaws_mesh_3d_get_vertex_count(hi) : 0, i, j, k;
	enum { R = 256 };
	qaws_scalar* tex = (qaws_scalar*)calloc(3 * R * R, sizeof(qaws_scalar));
	unsigned char* cov = (unsigned char*)malloc(R * R);
	qaws_texture t;
	char sub[200];
	if (!c || !hi || !tex || !cov)
		return;
	{
		/* the detailed ball's checker */
		double* x = ms_positions(hi);
		double* col = (double*)malloc(sizeof(double) * 3 * nv);
		double* cc;
		for (i = 0; i < nv; i++)
			ms_tx_checker(x + 3 * i, col + 3 * i);
		ms_at_add(hi, "color", QAWS_MESH_VERTEX, 3, QAWS_ATTRIBUTE_SRGB, col, nv);
		cc = ms_tx_corner_rgb(hi, "color");
		sprintf(sub, "%u faces, a colour per vertex", qaws_mesh_3d_get_face_count(hi));
		ms_panel(c, 0, 0, 600, 640, hi, -0.5, 0.35, "a detailed ball", sub);
		ms_at_draw(c, hi, cc, NULL, NULL, NULL, 0);
		printf("  detailed: %s\n", sub);
		free(x); free(col); free(cc);
	}
	{
		/* decimated, the colours carried per vertex */
		qaws_decimate_desc dd;
		qaws_mesh_3d_clone(&lo, hi);
		memset(&dd, 0, sizeof dd);
		dd.target_faces = 400;
		if (lo && qaws_mesh_3d_decimate(NULL, lo, &dd) == QAWS_STATUS_OK)
		{
			double* cc = ms_tx_corner_rgb(lo, "color");
			sprintf(sub, "%u faces, colours per vertex", qaws_mesh_3d_get_face_count(lo));
			ms_panel(c, 600, 0, 600, 640, lo, -0.5, 0.35, "decimated", sub);
			ms_at_draw(c, lo, cc, NULL, NULL, wire, 0.3);
			printf("  decimated: %s\n", sub);
			free(cc);
		}
	}
	if (lo)
	{
		/* the decimated ball's charts, the texture baked from the detailed ball */
		unsigned int nt = qaws_mesh_3d_get_triangle_count(lo), written = 0;
		qaws_scalar* uv = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 6 * (nt + 1));
		qaws_param_desc pd;
		qaws_param_report pr;
		qaws_texture_desc td;
		qaws_attribute_desc ad;
		unsigned int id;
		memset(&pd, 0, sizeof pd);
		pd.method = QAWS_PARAM_LSCM;
		pd.auto_seams = 1;
		pd.pack = 1;
		memset(&pr, 0, sizeof pr);
		if (qaws_mesh_3d_parameterize(uv, &pr, lo, NULL, 0, &pd) == QAWS_STATUS_OK)
		{
			memset(&ad, 0, sizeof ad);
			ad.name = "uv"; ad.element = QAWS_MESH_CORNER; ad.components = 2; ad.rule = QAWS_ATTRIBUTE_LINEAR;
			qaws_mesh_3d_add_attribute(&id, lo, uv, &ad);
			memset(&t, 0, sizeof t);
			t.texels = tex; t.width = R; t.height = R; t.components = 3;
			for (i = 0; i < 3 * R * R; i++)
				tex[i] = 1;   /* (white where no chart) */
			memset(&td, 0, sizeof td);
			td.gutter = 4;
			if (qaws_mesh_3d_bake_texture_from(cov, &t, lo, hi, "color", &td) == QAWS_STATUS_OK)
			{
				/* the texture itself: a quad per texel on a unit square */
				qaws_mesh_3d* sq = NULL;
				qaws_scalar sx[12] = { 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0 };
				unsigned int st[6] = { 0, 1, 2, 0, 2, 3 };
				qaws_mesh_3d_desc md;
				for (i = 0; i < R * R; i++)
					written += cov[i] == 1;
				memset(&md, 0, sizeof md);
				md.positions = sx; md.vertex_count = 4; md.face_vertices = st; md.face_count = 2;
				qaws_mesh_3d_create(&sq, &md);
				sprintf(sub, "%u x %u, %u charts, %u texels written", R, R, pr.charts, written);
				ms_panel(c, 1200, 0, 600, 640, sq, 0, 1.5707963, "its texture, baked", sub);
				for (j = 0; j < R; j++)
					for (i = 0; i < R; i++)
					{
						double p00[3] = { (double)i / R, (double)j / R, 0 }, p10[3] = { (double)(i + 1) / R, (double)j / R, 0 };
						double p11[3] = { (double)(i + 1) / R, (double)(j + 1) / R, 0 }, p01[3] = { (double)i / R, (double)(j + 1) / R, 0 }, col[3];
						for (k = 0; k < 3; k++)
							col[k] = 255 * (double)tex[3 * (j * R + i) + k];
						er_triangle(c, p00, p10, p11, NULL, NULL, NULL, col, col, col);
						er_triangle(c, p00, p11, p01, NULL, NULL, NULL, col, col, col);
					}
				printf("  texture: %s\n", sub);
				qaws_mesh_3d_destroy(sq);
				{
					/* the decimated ball with it: each triangle cut for the drawing, its pieces' corners read from the texture;
					   the colours at the pieces' points against the detailed ball's checker there */
					enum { K = 12 };
					double* x = ms_positions(lo);
					unsigned int* tri = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (nt + 1));
					unsigned int far = 0, farv = 0, pts = 0, el;
					double* vc = ms_at_get(lo, "color", &el);
					qaws_mesh_3d_get_triangles(tri, lo);
					ms_panel(c, 1800, 0, 600, 640, lo, -0.5, 0.35, "decimated, textured", NULL);
					for (i = 0; i < nt; i++)
					{
						unsigned int a, b, q;
						double P[3][3], U[3][2];
						for (q = 0; q < 3; q++)
						{
							for (k = 0; k < 3; k++)
								P[q][k] = x[3 * tri[3 * i + q] + k];
							U[q][0] = (double)uv[6 * i + 2 * q];
							U[q][1] = (double)uv[6 * i + 2 * q + 1];
						}
						for (b = 0; b < K; b++)
							for (a = 0; a + b < K; a++)
							{
								/* the piece's corners (a, b), (a + 1, b), (a, b + 1) and, past the first, the other half */
								unsigned int half;
								for (half = 0; half < (a + b + 1 < K ? 2u : 1u); half++)
								{
									unsigned int g[3][2], z;
									double Q[3][3], C[3][3];
									if (!half)
									{
										g[0][0] = a; g[0][1] = b; g[1][0] = a + 1; g[1][1] = b; g[2][0] = a; g[2][1] = b + 1;
									}
									else
									{
										g[0][0] = a + 1; g[0][1] = b; g[1][0] = a + 1; g[1][1] = b + 1; g[2][0] = a; g[2][1] = b + 1;
									}
									for (z = 0; z < 3; z++)
									{
										double w1 = (double)g[z][0] / K, w2 = (double)g[z][1] / K, w0 = 1 - w1 - w2, ruv[2];
										qaws_scalar suv[2], sc[3];
										for (k = 0; k < 3; k++)
											Q[z][k] = w0 * P[0][k] + w1 * P[1][k] + w2 * P[2][k];
										ruv[0] = w0 * U[0][0] + w1 * U[1][0] + w2 * U[2][0];
										ruv[1] = w0 * U[0][1] + w1 * U[1][1] + w2 * U[2][1];
										suv[0] = (qaws_scalar)ruv[0]; suv[1] = (qaws_scalar)ruv[1];
										qaws_texture_sample(sc, &t, suv, 1);
										for (k = 0; k < 3; k++)
											C[z][k] = 255 * (double)sc[k];
										if (!half)
										{
											double want[3], d2 = 0, dv2 = 0;
											ms_tx_checker(Q[z], want);
											for (k = 0; k < 3; k++)
											{
												/* (the per-vertex colours there, interpolated) */
												double pv = vc ? w0 * vc[3 * tri[3 * i] + k] + w1 * vc[3 * tri[3 * i + 1] + k] + w2 * vc[3 * tri[3 * i + 2] + k] : 0;
												d2 += (C[z][k] / 255 - want[k]) * (C[z][k] / 255 - want[k]);
												dv2 += (pv - want[k]) * (pv - want[k]);
											}
											far += sqrt(d2) > 0.25;
											farv += sqrt(dv2) > 0.25;
											pts++;
										}
									}
									er_triangle(c, Q[0], Q[1], Q[2], NULL, NULL, NULL, C[0], C[1], C[2]);
								}
							}
					}
					sprintf(sub, "off the checker %.1f%% (vertices %.1f%%)", 100.0 * far / pts, 100.0 * farv / pts);
					er_text(c, 1800 + 14, 46, sub, 2, 95, 95, 95);
					printf("  textured: %s\n", sub);
					free(x); free(tri); free(vc);
				}
			}
		}
		free(uv);
	}
	er_save(c, "showcase/mesh81_textures.png");
	printf("mesh81_textures.png\n");
	qaws_mesh_3d_destroy(hi);
	qaws_mesh_3d_destroy(lo);
	free(tex); free(cov);
	er_destroy(c);
}
