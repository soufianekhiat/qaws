/*
 * Figure 14: quad finishing (valence pairs cancelled by local quad
 * operations, relaxation by angles guarded by the scaled Jacobian).
 *   - the can (from below: its rims' slivers stay) and the blob, all quads,
 *     and the wavy patch generated in quads, as INSTANT leaves them now:
 *     each quad coloured by its smallest scaled Jacobian (light at 1, amber
 *     toward 0.5, red under), irregular vertices dotted (valence 3 blue, 5
 *     red, others purple); the numbers before the finishing (the plain
 *     relaxation, measured once in test 135) in the subtitles;
 *   - a planar 8 x 8 grid with one edge turned (two 3-5 pairs), and the
 *     same after the finishing: the grid again.
 * The grid goes through the internal entry (qi_improve), the rest through
 * qaws_mesh_quadrangulate and qaws_mesh_generate_surface.
 */

#include "mesh/qaws_quad_improve.h"   /* (the synthetic grid only) */
#include <time.h>

/* a quad's smallest scaled Jacobian about its Newell normal */
static double ms_quad_sj(double const* p, unsigned int const* v)
{
	double n[3] = { 0, 0, 0 }, ln, worst = 1;
	unsigned int k, i;
	for (k = 0; k < 4; k++)
	{
		double const *a = p + 3 * v[k], *b = p + 3 * v[(k + 1) & 3];
		n[0] += (a[1] - b[1]) * (a[2] + b[2]);
		n[1] += (a[2] - b[2]) * (a[0] + b[0]);
		n[2] += (a[0] - b[0]) * (a[1] + b[1]);
	}
	ln = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
	for (k = 0; k < 4; k++)
	{
		double const *o = p + 3 * v[k], *a = p + 3 * v[(k + 1) & 3], *b = p + 3 * v[(k + 3) & 3];
		double e1[3], e2[3], l1, l2, s;
		for (i = 0; i < 3; i++)
		{
			e1[i] = a[i] - o[i];
			e2[i] = b[i] - o[i];
		}
		l1 = sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
		l2 = sqrt(e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2]);
		s = l1 > 0 && l2 > 0 && ln > 0 ? ((e1[1] * e2[2] - e1[2] * e2[1]) * n[0] + (e1[2] * e2[0] - e1[0] * e2[2]) * n[1] 
			+ (e1[0] * e2[1] - e1[1] * e2[0]) * n[2]) / (ln * l1 * l2) : -1;
		worst = fmin(worst, s);
	}
	return worst;
}

/* the mesh drawn face by face coloured by its smallest scaled Jacobian, its edges, its irregular vertices dotted; the
   numbers: SJ min and mean, AD (mean |angle - 90| at the quads' corners), irregular count */
static void ms_sj_draw(er_canvas* c, qaws_mesh_3d const* m, double* out_min, double* out_mean, double* out_ad, unsigned int* out_irr)
{
	static double const wire[3] = { 45, 50, 62 }, good[3] = { 196, 210, 228 }, warn[3] = { 240, 196, 120 }, red[3] = { 214, 39, 40 },
		blue[3] = { 31, 90, 200 }, purple[3] = { 140, 60, 170 };
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), nf = qaws_mesh_3d_get_face_count(m), ne = qaws_mesh_3d_get_edge_count(m);
	unsigned int *fv = (unsigned int*)malloc(sizeof(unsigned int) * (qaws_mesh_3d_get_corner_count(m) + 1));
	unsigned int *fs = (unsigned int*)malloc(sizeof(unsigned int) * (nf + 1));
	unsigned int *e = (unsigned int*)malloc(sizeof(unsigned int) * (2 * ne + 2)), *deg = (unsigned int*)calloc(nv + 1, sizeof(unsigned int));
	unsigned char* border = (unsigned char*)calloc(nv + 1, 1);
	unsigned int *across = (unsigned int*)malloc(sizeof(unsigned int) * (qaws_mesh_3d_get_corner_count(m) + 1)), f, c0, k, nq = 0, nc = 0;
	double* p = ms_positions(m);
	double* n = (double*)malloc(sizeof(double) * 3 * (nv + 1));
	qaws_scalar* ns = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (nv + 1));
	double smin = 1, ssum = 0, ad = 0;
	qaws_mesh_3d_get_faces(m, fv, fs);
	qaws_mesh_3d_get_edges(m, e);
	qaws_mesh_3d_get_face_adjacency(m, across);
	qaws_mesh_3d_get_vertex_normals(m, ns);
	for (k = 0; k < 3 * nv; k++)
		n[k] = (double)ns[k];
	for (f = 0, c0 = 0; f < nf; c0 += fs[f], f++)
	{
		double col[3], s = fs[f] == 4 ? ms_quad_sj(p, fv + c0) : 1, t;
		for (k = 0; k < fs[f]; k++)
		{
			double const *o = p + 3 * fv[c0 + k], *a = p + 3 * fv[c0 + (k + 1) % fs[f]], *b = p + 3 * fv[c0 + (k + fs[f] - 1) % fs[f]];
			double d = ((a[0] - o[0]) * (b[0] - o[0]) + (a[1] - o[1]) * (b[1] - o[1]) + (a[2] - o[2]) * (b[2] - o[2]))
				/ sqrt(((a[0] - o[0]) * (a[0] - o[0]) + (a[1] - o[1]) * (a[1] - o[1]) + (a[2] - o[2]) * (a[2] - o[2]))
					* ((b[0] - o[0]) * (b[0] - o[0]) + (b[1] - o[1]) * (b[1] - o[1]) + (b[2] - o[2]) * (b[2] - o[2])));
			ad += fabs(acos(fmax(-1, fmin(1, d))) * 180 / MS_PI - 180.0 * (fs[f] - 2) / fs[f]);
			nc++;
			if (across[c0 + k] == QAWS_MESH_NONE)
				border[fv[c0 + k]] = border[fv[c0 + (k + 1) % fs[f]]] = 1;
		}
		if (fs[f] == 4)
		{
			smin = fmin(smin, s);
			ssum += s;
			nq++;
		}
		/* light at 1, amber toward 0.5, red under */
		t = (1 - s) / 0.5;
		for (k = 0; k < 3; k++)
			col[k] = s < 0.5 ? red[k] : good[k] + fmin(1, t) * (warn[k] - good[k]);
		for (k = 1; k + 1 < fs[f]; k++)
			er_triangle(c, p + 3 * fv[c0], p + 3 * fv[c0 + k], p + 3 * fv[c0 + k + 1], n + 3 * fv[c0], n + 3 * fv[c0 + k], n + 3 * fv[c0 + k + 1],
				col, col, col);
	}
	for (k = 0; k < ne; k++)
	{
		er_line(c, p + 3 * e[2 * k], p + 3 * e[2 * k + 1], wire, 1.0);
		deg[e[2 * k]]++;
		deg[e[2 * k + 1]]++;
	}
	*out_irr = 0;
	for (k = 0; k < nv; k++)
		if (!border[k] && deg[k] != 4)
		{
			er_point(c, p + 3 * k, deg[k] == 3 ? blue : deg[k] == 5 ? red : purple, 6);
			(*out_irr)++;
		}
	*out_min = smin;
	*out_mean = nq ? ssum / nq : 0;
	*out_ad = nc ? ad / nc : 0;
	free(fv); free(fs); free(e); free(deg); free(border); free(across); free(p); free(n); free(ns);
}

/* an 8 x 8 grid of unit quads in z = 0 */
static qaws_mesh_3d* ms_grid8(void)
{
	double x[3 * 81];
	unsigned int f[4 * 64], fs[64], i, j, k = 0;
	for (j = 0; j <= 8; j++)
		for (i = 0; i <= 8; i++)
		{
			x[3 * (9 * j + i)] = i;
			x[3 * (9 * j + i) + 1] = j;
			x[3 * (9 * j + i) + 2] = 0;
		}
	for (j = 0; j < 8; j++)
		for (i = 0; i < 8; i++)
		{
			f[k++] = 9 * j + i; f[k++] = 9 * j + i + 1; f[k++] = 9 * j + i + 10; f[k++] = 9 * j + i + 9;
			fs[8 * j + i] = 4;
		}
	return ms_make(x, 81, f, fs, 64);
}

static void demo_quad_improve(void)
{
	/* before the finishing (qd_relax alone, measured once in test 135) */
	static char const* const was[3] = { "before: 22 irregular, SJ min -0.60, AD 5.8", "before: 12 irregular, SJ min 0.72, AD 5.7",
		"before: 4 irregular, SJ min 0.80, AD 4.6" };
	er_canvas* c = er_create(2750, 640, 3, 252, 252, 250);
	qaws_mesh_3d *can = ms_can(8, 1, 1.6), *blob = ms_blob(1), *q[3] = { NULL, NULL, NULL }, *g = NULL;
	qaws_surface* wavy = ms_wavy();
	qaws_mesh_report rep;
	char sub[200];
	double smin, smean, ad, ms[3] = { 0, 0, 0 };
	unsigned int irr, i;
	static char const* const title[3] = { "the can, from below", "the blob", "the wavy patch" };
	if (!c || !can || !blob || !wavy)
		return;
	qaws_mesh_3d_detect_features(can, (qaws_scalar)(50 * MS_PI / 180), 0, 0);
	for (i = 0; i < 3; i++)
	{
		clock_t t0 = clock();
		q[i] = i == 0 ? ms_quads(can, 0.14, 1, &rep) : i == 1 ? ms_quads(blob, 0.11, 1, &rep) : ms_generate(wavy, NULL, QAWS_ELEMENTS_QUADS, 0.1, 0, &rep);
		ms[i] = 1000.0 * (double)(clock() - t0) / CLOCKS_PER_SEC;
		if (!q[i])
			continue;
		ms_panel(c, 550 * (int)i, 0, 550, 640, q[i], i == 0 ? 0.35 : 0.6, i == 0 ? -0.55 : i == 2 ? 0.65 : 0.5, title[i], was[i]);
		ms_sj_draw(c, q[i], &smin, &smean, &ad, &irr);
		sprintf(sub, "after:  %u irregular, SJ min %.2f, AD %.1f", irr, smin, ad);
		er_text(c, 550 * (int)i + 14, 76, sub, 2, 95, 95, 95);
		printf("  %s: %u faces, irregular %u, SJ min %.3f mean %.3f, AD %.2f deg, edges %.3f..%.3f, %.0f ms\n", title[i] + 4, rep.quality.quads, irr, smin,
			smean, ad, (double)rep.quality.min_edge, (double)rep.quality.max_edge, ms[i]);
	}
	/* the grid with an edge turned, then finished */
	g = ms_grid8();
	if (g)
	{
		qaws_mesh_reference* r = NULL;
		qaws_hm* m = &g->hm;
		unsigned char kind[82];
		unsigned int v, ops = 0;
		qaws_mesh_reference_create(g, NULL, &r);
		hm_foot_enable(m);
		for (v = 0; v < m->nv; v++)
		{
			qaws_mesh_reference_get_vertex_location(r, v, &m->foot[v]);
			kind[v] = (unsigned char)(hm_is_border_vertex(m, v) ? hm_valence(m, v) == 2 ? 2 : 1 : 0);
			if (kind[v] == 2)
				m->vflag[v] |= HM_CORNER;
		}
		hm_rotate(m, hm_find(m, 39, 40), 1);
		for (i = 0; i < 2; i++)
		{
			ms_panel(c, 1650 + 550 * (int)i, 0, 550, 640, g, 0.0, 1.2, i ? "the grid, finished" : "a grid, one edge turned", NULL);
			ms_sj_draw(c, g, &smin, &smean, &ad, &irr);
			if (i)
				sprintf(sub, "%u operation(s): %u irregular, SJ min %.2f", ops, irr, smin);
			else
				sprintf(sub, "two 3-5 pairs: %u irregular, SJ min %.2f", irr, smin);
			er_text(c, 1650 + 550 * (int)i + 14, 46, sub, 2, 95, 95, 95);
			printf("  grid %s: irregular %u, SJ min %.3f\n", i ? "finished" : "turned", irr, smin);
			if (!i && r)
			{
				mr_scratch s;
				if (mr_scratch_init(&s, r->nt))
				{
					ops = qi_improve(m, r, &s, kind, 1, 4);
					mr_scratch_free(&s);
				}
				printf("  grid: %u operation(s)\n", ops);
			}
		}
		qaws_mesh_reference_destroy(r);
	}
	er_save(c, "showcase/mesh14_quad_improve.png");
	printf("mesh14_quad_improve.png\n");
	for (i = 0; i < 3; i++)
		qaws_mesh_3d_destroy(q[i]);
	qaws_mesh_3d_destroy(can);
	qaws_mesh_3d_destroy(blob);
	qaws_mesh_3d_destroy(g);
	qaws_surface_destroy(wavy);
	er_destroy(c);
}
