/*
 * Figure 13: quads by min-cost flow (QuadriFlow's consistency).
 *   - the blob's dense mesh: the faces holding a dislocation once the steps
 *     are rounded (orange), the edges the flow changed (by cost, blue to
 *     red), the cones (singular faces, black);
 *   - the blob quad-dominant: faces other than quads in red;
 *   - the blob all quads at the target spacing, the irregular vertices
 *     marked (valence 3 blue, 5 red);
 *   - the can's cap, all quads, close up.
 */

#include "mesh/qaws_quad_flow.h"
#include <time.h>

/* a dense remesh at a third of rho, its fields (free cross field), its steps before (T0) and after consistency (T) */
typedef struct mq_dense
{
	qaws_mesh_3d* dense;
	mf_hier h;
	int *T0, *T;
	unsigned char* K;
	double* Tc;
	qf_report rep;
	unsigned int left, odd;
	double seconds;
} mq_dense;

static int mq_run(qaws_mesh_3d const* m, double rho, mq_dense* d)
{
	qaws_remesh_desc rd;
	unsigned int nc;
	clock_t t0;
	memset(d, 0, sizeof *d);
	memset(&rd, 0, sizeof rd);
	rd.target_length = (qaws_scalar)(rho / 3);
	rd.iterations = 6;
	if (qaws_mesh_remesh(&d->dense, NULL, m, &rd) != QAWS_STATUS_OK || mf_build(&d->h, &d->dense->hm, NULL) != QAWS_STATUS_OK)
		return 0;
	mf_orientations(&d->h, 8, NULL, 0);
	mf_positions(&d->h, rho, 8);
	nc = d->dense->hm.nc;
	d->T0 = (int*)malloc(sizeof(int) * 2 * nc);
	d->T = (int*)malloc(sizeof(int) * 2 * nc);
	d->K = (unsigned char*)malloc(nc);
	d->Tc = (double*)malloc(sizeof(double) * 4 * nc);
	qf_steps(&d->dense->hm, &d->h.lv[0], rho, d->T, d->K, d->Tc);
	memcpy(d->T0, d->T, sizeof(int) * 2 * nc);
	t0 = clock();
	d->left = qf_consistent(&d->dense->hm, &d->h.lv[0], d->T, d->K, d->Tc, &d->odd, &d->rep);
	d->seconds = (double)(clock() - t0) / CLOCKS_PER_SEC;
	return 1;
}

static void mq_free(mq_dense* d)
{
	qaws_mesh_3d_destroy(d->dense);
	mf_free(&d->h);
	free(d->T0); free(d->T); free(d->K); free(d->Tc);
}

/* the dense mesh, its dislocated faces before the flow, the changed edges by cost, the cones */
static void mq_draw_dense(er_canvas* c, mq_dense const* d)
{
	static double const tint[3] = { 214, 218, 226 }, orange[3] = { 240, 140, 30 }, black[3] = { 20, 20, 20 };
	qaws_hm const* m = &d->dense->hm;
	unsigned char* fr = (unsigned char*)malloc(m->nc);
	unsigned int f, k, cc;
	ms_draw(c, d->dense, tint, NULL, 1, NULL, 0);
	for (f = 0; f < m->nf; f++)
	{
		int S[2], t = qf_closure(m, f, d->T0, d->K, fr, S);
		double cen[3] = { 0, 0, 0 };
		for (k = 0; k < 3; k++)
			for (cc = 0; cc < 3; cc++)
				cen[cc] += m->x[3 * (size_t)m->cv[m->fc[f] + k] + cc] / 3;
		if (t)
			er_point(c, cen, black, 6);
		else if (S[0] || S[1])
			er_point(c, cen, orange, 5);
	}
	for (cc = 0; cc < m->nc; cc++)
	{
		double cost = 0, rgb[3];
		unsigned int j;
		if (m->cf[cc] == HM_NONE || (m->co[cc] != HM_NONE && m->co[cc] < cc))
			continue;
		if (d->T[2 * cc] == d->T0[2 * cc] && d->T[2 * cc + 1] == d->T0[2 * cc + 1])
			continue;
		for (j = 0; j < 2; j++)
		{
			/* (as the flow priced it: how far the cheaper end moves off its lattice point) */
			int s = d->T[2 * cc + j] - d->T0[2 * cc + j];
			double fa = d->Tc[4 * cc + j], fb = d->Tc[4 * cc + 2 + j], u = fabs(fb - s) - fabs(fb), w = fabs(fa + s) - fabs(fa);
			cost += s ? (u < w ? u : w) : 0;
		}
		ms_colormap(cost, rgb);
		er_line(c, m->x + 3 * (size_t)m->cv[cc], m->x + 3 * (size_t)hm_to(m, cc), rgb, 2.6);
	}
	free(fr);
}

/* the irregular vertices of a quad mesh off its border: valence 3 blue, 5 or more red */
static void mq_mark_irregular(er_canvas* c, qaws_mesh_3d const* q)
{
	static double const red[3] = { 214, 39, 40 }, blue[3] = { 31, 90, 200 };
	qaws_hm const* m = &q->hm;
	unsigned int v;
	for (v = 0; v < m->nv; v++)
	{
		unsigned int val;
		if (m->vc[v] == HM_NONE || hm_is_border_vertex(m, v))
			continue;
		val = hm_valence(m, v);
		if (val != 4)
			er_point(c, m->x + 3 * (size_t)v, val < 4 ? blue : red, 8);
	}
}

/* the faces other than quads, filled red */
static void mq_mark_others(er_canvas* c, qaws_mesh_3d const* q)
{
	static double const red[3] = { 214, 39, 40 };
	qaws_hm const* m = &q->hm;
	unsigned int f, k;
	for (f = 0; f < m->nf; f++)
	{
		double const* a;
		if (m->fn[f] == 4)
			continue;
		a = m->x + 3 * (size_t)m->cv[m->fc[f]];
		for (k = 1; k + 1 < m->fn[f]; k++)
			er_triangle(c, a, m->x + 3 * (size_t)m->cv[m->fc[f] + k], m->x + 3 * (size_t)m->cv[m->fc[f] + k + 1], NULL, NULL, NULL, red, red, red);
	}
}

static void demo_quad_flow(void)
{
	static double const wire[3] = { 40, 45, 60 }, tint[3] = { 178, 192, 212 }, warm[3] = { 226, 200, 160 };
	er_canvas* c = er_create(2200, 640, 3, 252, 252, 250);
	qaws_mesh_3d *blob = ms_blob(1), *can = ms_can(8, 1, 1.6), *q1 = NULL, *q2 = NULL, *q3 = NULL;
	qaws_mesh_report rep;
	mq_dense d;
	char sub[200];
	clock_t t0;
	if (!c || !blob || !can)
		return;
	/* 1. the flow on the blob's dense mesh (the spacing of panel 2, a free cross field) */
	if (mq_run(blob, 0.16, &d))
	{
		sprintf(sub, "%u dislocations, %u unit changes", d.rep.supply, d.rep.changed);
		ms_panel(c, 0, 0, 550, 640, d.dense, 0.6, 0.5, "the flow's changes", sub);
		mq_draw_dense(c, &d);
		printf("  blob dense (%u faces) at 0.16: supplies %u, cut %u, unrouted %u, cost %.1f, changed %u, parity paths %u, cones placed %u, "
			"flips %u -> %u, dislocations left %u, odd cones %u, steps of two %u, %.3f s\n", d.dense->hm.nf, d.rep.supply, d.rep.cut, d.rep.unrouted,
			d.rep.cost * 1e-6, d.rep.changed, d.rep.paths, d.rep.cones, d.rep.flips, d.rep.flips_left, d.left, d.odd, d.rep.twos, d.seconds);
		mq_free(&d);
	}
	/* 2. quad-dominant */
	t0 = clock();
	q1 = ms_quads(blob, 0.16, 0, &rep);
	if (q1)
	{
		/* (the quad count on stdout: the line fits the panel without it) */
		sprintf(sub, "%u others (were 6), %u irregular (were 8)", rep.quality.triangles + rep.quality.polygons, rep.quality.irregular);
		ms_panel(c, 550, 0, 550, 640, q1, 0.6, 0.5, "quad-dominant", sub);
		ms_draw(c, q1, tint, NULL, 1, wire, 1.1);
		mq_mark_others(c, q1);
		printf("  blob quad-dominant 0.16: %u quads, %u triangles, %u polygons, irregular %u, mean edge %.3f, %.2f s\n", rep.quality.quads,
			rep.quality.triangles, rep.quality.polygons, rep.quality.irregular, (double)rep.quality.mean_edge, (double)(clock() - t0) / CLOCKS_PER_SEC);
	}
	/* 3. all quads at the target */
	t0 = clock();
	q2 = ms_quads(blob, 0.11, 1, &rep);
	if (q2)
	{
		sprintf(sub, "%u quads, %u irregular (were 12)", rep.quality.quads, rep.quality.irregular);
		ms_panel(c, 1100, 0, 550, 640, q2, 0.6, 0.5, "all quads at 0.11", sub);
		ms_draw(c, q2, tint, NULL, 1, wire, 1.0);
		mq_mark_irregular(c, q2);
		printf("  blob pure 0.11: %u quads, irregular %u, mean edge %.3f, %.2f s\n", rep.quality.quads, rep.quality.irregular,
			(double)rep.quality.mean_edge, (double)(clock() - t0) / CLOCKS_PER_SEC);
	}
	/* 4. the can's cap, close up */
	qaws_mesh_3d_detect_features(can, (qaws_scalar)(50 * MS_PI / 180), 0, 0);
	t0 = clock();
	q3 = ms_quads(can, 0.14, 1, &rep);
	if (q3)
	{
		unsigned int n = qaws_mesh_3d_get_vertex_count(q3), i, k = 0;
		double *p = ms_positions(q3), *cap = (double*)malloc(sizeof(double) * 3 * (n + 1));
		for (i = 0; i < n; i++)
			if (p[3 * i + 2] > 0.7)
			{
				memcpy(cap + 3 * k, p + 3 * i, sizeof(double) * 3);
				k++;
			}
		er_panel(c, 1650, 72, 550, 640 - 72, cap, k, 0.3, 1.2, 0.08);
		sprintf(sub, "%u quads, %u irregular (were 22)", rep.quality.quads, rep.quality.irregular);
		er_text(c, 1650 + 14, 12, "the can's cap", 3, 30, 30, 30);
		er_text(c, 1650 + 14, 46, sub, 2, 95, 95, 95);
		ms_draw(c, q3, warm, NULL, 1, wire, 1.0);
		mq_mark_irregular(c, q3);
		printf("  can pure 0.14: %u quads, irregular %u, mean edge %.3f, %.2f s\n", rep.quality.quads, rep.quality.irregular,
			(double)rep.quality.mean_edge, (double)(clock() - t0) / CLOCKS_PER_SEC);
		free(p);
		free(cap);
	}
	er_save(c, "showcase/mesh13_quad_flow.png");
	printf("mesh13_quad_flow.png\n");
	qaws_mesh_3d_destroy(blob);
	qaws_mesh_3d_destroy(can);
	qaws_mesh_3d_destroy(q1);
	qaws_mesh_3d_destroy(q2);
	qaws_mesh_3d_destroy(q3);
	er_destroy(c);
}
