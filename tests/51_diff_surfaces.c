/*
 * Test 51: Differentiable surface evaluation
 *
 * Tensor-product families (Bezier, B-spline, bilinear, biquadratic, NURBS):
 *   - primal jets match qaws_surface_evaluate, third order matches
 *     finite differences of second order
 *   - parameter and coordinate tangents, second tangent
 *   - batch adjoint identity and accumulation equivalence
 *   - knot-line report, support index, normal composition
 */

#include "test_diff.h"

#define SURF_MAX_CP 64
#define SURF_SAMPLES 7

/* Scalars in the combined parameter vector: xyz per point, plus one weight when rational. */
#define PARAMS(f) ((f)->cp_count * ((f)->rational ? 4u : 3u))

typedef struct surface_family
{
	char const* name;
	unsigned int cp_count;
	unsigned int u_count, v_count, u_degree, v_degree;
	qaws_scalar cps[SURF_MAX_CP * 4];   /* control points, then weights when rational */
	qaws_scalar u_knots[16], v_knots[16];
	unsigned int u_knot_count, v_knot_count;
	int kind; /* 0 bezier, 1 bspline, 2 bilinear, 3 biquadratic, 4 nurbs */
	int rational;
} surface_family;

static qaws_surface* family_surface(surface_family const* f, qaws_scalar const* cps)
{
	qaws_surface* s = NULL;
	switch (f->kind)
	{
	case 0:
	{
		qaws_surface_bezier_desc d;
		d.u_degree = f->u_degree;
		d.v_degree = f->v_degree;
		d.control_points = (qaws_vec3 const*)cps;
		d.u_point_count = f->u_count;
		d.v_point_count = f->v_count;
		qaws_surface_create_bezier(&d, &s);
		break;
	}
	case 1:
	{
		qaws_surface_bspline_desc d;
		d.u_degree = f->u_degree;
		d.v_degree = f->v_degree;
		d.control_points = (qaws_vec3 const*)cps;
		d.u_point_count = f->u_count;
		d.v_point_count = f->v_count;
		d.u_knots = f->u_knots;
		d.u_knot_count = f->u_knot_count;
		d.v_knots = f->v_knots;
		d.v_knot_count = f->v_knot_count;
		qaws_surface_create_bspline(&d, &s);
		break;
	}
	case 2:
	{
		qaws_surface_bilinear_desc d;
		memcpy(&d.p00, cps + 0, sizeof(qaws_vec3));
		memcpy(&d.p10, cps + 3, sizeof(qaws_vec3));
		memcpy(&d.p01, cps + 6, sizeof(qaws_vec3));
		memcpy(&d.p11, cps + 9, sizeof(qaws_vec3));
		qaws_surface_create_bilinear(&d, &s);
		break;
	}
	case 4:
	{
		qaws_surface_nurbs_desc d;
		d.u_degree = f->u_degree;
		d.v_degree = f->v_degree;
		d.control_points = (qaws_vec3 const*)cps;
		d.u_point_count = f->u_count;
		d.v_point_count = f->v_count;
		d.weights = cps + f->cp_count * 3;
		d.u_knots = f->u_knots;
		d.u_knot_count = f->u_knot_count;
		d.v_knots = f->v_knots;
		d.v_knot_count = f->v_knot_count;
		qaws_surface_create_nurbs(&d, &s);
		break;
	}
	default:
	{
		qaws_surface_biquadratic_desc d;
		memcpy(d.control_points, cps, sizeof(qaws_vec3) * 9);
		qaws_surface_create_biquadratic(&d, &s);
		break;
	}
	}
	return s;
}

static qaws_surface* family_surface_shifted(surface_family const* f, qaws_scalar const* dir, double h)
{
	qaws_scalar p[SURF_MAX_CP * 4];
	unsigned int i;
	for (i = 0; i < PARAMS(f); i++)
		p[i] = (qaws_scalar)(f->cps[i] + h * dir[i]);
	return family_surface(f, p);
}

static void make_surface_families(surface_family* fams)
{
	unsigned int i;
	memset(fams, 0, sizeof(surface_family) * 5);
	diff_seed(31337);

	fams[0].name = "bezier";
	fams[0].kind = 0;
	fams[0].u_degree = 2; fams[0].v_degree = 3;
	fams[0].u_count = 3; fams[0].v_count = 4;

	fams[1].name = "bspline";
	fams[1].kind = 1;
	fams[1].u_degree = 3; fams[1].v_degree = 2;
	fams[1].u_count = 6; fams[1].v_count = 5;
	{
		static qaws_scalar const uk[10] = { 0, 0, 0, 0, 0.3f, 0.7f, 1, 1, 1, 1 };
		static qaws_scalar const vk[8] = { 0, 0, 0, 0.25f, 0.55f, 1, 1, 1 };
		for (i = 0; i < 10; i++) fams[1].u_knots[i] = uk[i];
		for (i = 0; i < 8; i++) fams[1].v_knots[i] = vk[i];
		fams[1].u_knot_count = 10;
		fams[1].v_knot_count = 8;
	}

	fams[2].name = "bilinear";
	fams[2].kind = 2;
	fams[2].u_degree = 1; fams[2].v_degree = 1;
	fams[2].u_count = 2; fams[2].v_count = 2;

	fams[3].name = "biquadratic";
	fams[3].kind = 3;
	fams[3].u_degree = 2; fams[3].v_degree = 2;
	fams[3].u_count = 3; fams[3].v_count = 3;

	fams[4] = fams[1];
	fams[4].name = "nurbs";
	fams[4].kind = 4;
	fams[4].rational = 1;

	for (i = 0; i < 5; i++)
	{
		unsigned int a, b;
		fams[i].cp_count = fams[i].u_count * fams[i].v_count;
		/* A gently curved sheet with random relief keeps normals well defined. */
		for (a = 0; a < fams[i].u_count; a++)
			for (b = 0; b < fams[i].v_count; b++)
			{
				unsigned int e = (fams[i].kind == 2 || fams[i].kind == 3) ? b * fams[i].u_count + a : a * fams[i].v_count + b;
				fams[i].cps[e * 3 + 0] = (qaws_scalar)a + diff_rand() * (qaws_scalar)0.2;
				fams[i].cps[e * 3 + 1] = (qaws_scalar)b + diff_rand() * (qaws_scalar)0.2;
				fams[i].cps[e * 3 + 2] = diff_rand();
			}
		if (fams[i].rational)
			for (a = 0; a < fams[i].cp_count; a++)
				fams[i].cps[fams[i].cp_count * 3 + a] = (qaws_scalar)1.25 + (qaws_scalar)0.75 * diff_rand();
	}
}

static void sample_uv(unsigned int i, qaws_scalar* u, qaws_scalar* v)
{
	*u = ((qaws_scalar)i + (qaws_scalar)0.41) / (qaws_scalar)SURF_SAMPLES;
	*v = ((qaws_scalar)(SURF_SAMPLES - 1 - i) + (qaws_scalar)0.23) / (qaws_scalar)SURF_SAMPLES;
}

/* Componentwise, relative to the larger vector magnitude. */
static int vec_close(qaws_vec3 a, qaws_vec3 b, double tol)
{
	double scale = 1.0;
	double ma = fabs(a.x) > fabs(a.y) ? fabs(a.x) : fabs(a.y);
	double mb = fabs(b.x) > fabs(b.y) ? fabs(b.x) : fabs(b.y);
	if (fabs(a.z) > ma) ma = fabs(a.z);
	if (fabs(b.z) > mb) mb = fabs(b.z);
	if (ma > scale) scale = ma;
	if (mb > scale) scale = mb;
	return fabs(a.x - b.x) <= tol * scale && fabs(a.y - b.y) <= tol * scale && fabs(a.z - b.z) <= tol * scale;
}

static double sjet_dot(qaws_surface_jet const* a, qaws_surface_jet const* b, unsigned int channels)
{
	double s = 0;
	unsigned int i;
	for (i = 0; i < QAWS_SURFACE_JET_COUNT; i++)
		if (channels & (1u << i))
			s += (double)a->d[i].x * b->d[i].x + (double)a->d[i].y * b->d[i].y + (double)a->d[i].z * b->d[i].z;
	return s;
}

static void rand_sjet(qaws_surface_jet* j)
{
	unsigned int i;
	for (i = 0; i < QAWS_SURFACE_JET_COUNT; i++)
		j->d[i] = diff_rand_vec3();
	j->channels = QAWS_SJET_ORDER3;
}

static void cp_views(surface_family const* f, qaws_scalar* data, qaws_field_view* v, qaws_diff_views* views)
{
	v[0] = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, data, f->cp_count, 3);
	v[1] = qaws_field_view_make(QAWS_FIELD_WEIGHTS, data + f->cp_count * 3, f->cp_count, 1);
	views->fields = v;
	views->field_count = f->rational ? 2u : 1u;
	views->children = NULL;
	views->child_count = 0;
}

/* ------------------------------------------------------------------ */

static void check_primal(surface_family const* f, qaws_surface const* s)
{
	unsigned int i;
	int ok = 1, ok3 = 1;
	double h = DIFF_FD_STEP;
	for (i = 0; i < SURF_SAMPLES; i++)
	{
		qaws_scalar u, v;
		qaws_surface_jet j, jp, jm;
		qaws_surface_eval_result r;
		sample_uv(i, &u, &v);
		TEST_ASSERT_STATUS(qaws_surface_eval_jet(s, u, v, QAWS_SJET_ORDER3, &j));
		qaws_surface_evaluate(s, u, v, QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV |
			QAWS_SURFACE_EVAL_DUU | QAWS_SURFACE_EVAL_DUV | QAWS_SURFACE_EVAL_DVV, &r);
		if (!vec_close(j.d[0], r.position, 1e-4) || !vec_close(j.d[1], r.du, 1e-4) || !vec_close(j.d[2], r.dv, 1e-4) ||
		    !vec_close(j.d[3], r.duu, 1e-4) || !vec_close(j.d[4], r.duv, 1e-4) || !vec_close(j.d[5], r.dvv, 1e-4))
			ok = 0;

		/* Third order: finite differences of the analytic second order. */
		qaws_surface_eval_jet(s, (qaws_scalar)(u + h), v, QAWS_SJET_ORDER2, &jp);
		qaws_surface_eval_jet(s, (qaws_scalar)(u - h), v, QAWS_SJET_ORDER2, &jm);
		if (!diff_close((jp.d[3].x - jm.d[3].x) / (2 * h), j.d[6].x, DIFF_TOL * 100) ||
		    !diff_close((jp.d[4].y - jm.d[4].y) / (2 * h), j.d[7].y, DIFF_TOL * 100) ||
		    !diff_close((jp.d[5].z - jm.d[5].z) / (2 * h), j.d[8].z, DIFF_TOL * 100))
			ok3 = 0;
		qaws_surface_eval_jet(s, u, (qaws_scalar)(v + h), QAWS_SJET_ORDER2, &jp);
		qaws_surface_eval_jet(s, u, (qaws_scalar)(v - h), QAWS_SJET_ORDER2, &jm);
		if (!diff_close((jp.d[5].x - jm.d[5].x) / (2 * h), j.d[9].x, DIFF_TOL * 100))
			ok3 = 0;
	}
	printf("    %s: jets match evaluate: %s, third order consistent: %s\n", f->name, ok ? "yes" : "NO", ok3 ? "yes" : "NO");
	TEST_ASSERT(ok, "surface jet matches qaws_surface_evaluate");
	TEST_ASSERT(ok3, "third-order partials consistent with second order");
}

static void check_tangents(surface_family const* f, qaws_surface const* s)
{
	qaws_scalar dir[SURF_MAX_CP * 4];
	qaws_field_view fv[2];
	qaws_diff_views views;
	qaws_surface *sp, *sm;
	unsigned int i, ch;
	int ok_p = 1, ok_c = 1;
	double h = DIFF_FD_STEP;

	diff_rand_fill(dir, PARAMS(f));
	cp_views(f, dir, fv, &views);
	sp = family_surface_shifted(f, dir, h);
	sm = family_surface_shifted(f, dir, -h);

	for (i = 0; i < SURF_SAMPLES; i++)
	{
		qaws_scalar u, v;
		qaws_surface_jet p, t, jp, jm, tc;
		sample_uv(i, &u, &v);
		qaws_surface_eval_tangent(NULL, s, u, v, 0, 0, QAWS_SJET_ORDER3, &views, &p, &t);
		qaws_surface_eval_jet(sp, u, v, QAWS_SJET_ORDER3, &jp);
		qaws_surface_eval_jet(sm, u, v, QAWS_SJET_ORDER3, &jm);
		for (ch = 0; ch < QAWS_SURFACE_JET_COUNT; ch++)
		{
			qaws_vec3 fd;
			fd.x = (qaws_scalar)((jp.d[ch].x - jm.d[ch].x) / (2 * h));
			fd.y = (qaws_scalar)((jp.d[ch].y - jm.d[ch].y) / (2 * h));
			fd.z = (qaws_scalar)((jp.d[ch].z - jm.d[ch].z) / (2 * h));
			if (!vec_close(fd, t.d[ch], DIFF_TOL * 50))
				ok_p = 0;
		}

		/* u' = 1: S -> Su, Su -> Suu, Sv -> Suv, Suv -> Suuv, Svv -> Suvv */
		qaws_surface_eval_tangent(NULL, s, u, v, 1, 0, QAWS_SJET_ORDER2, NULL, &p, &tc);
		{
			qaws_surface_jet full;
			qaws_surface_eval_jet(s, u, v, QAWS_SJET_ORDER3, &full);
			if (!vec_close(tc.d[0], full.d[1], 1e-4) || !vec_close(tc.d[1], full.d[3], 1e-4) ||
			    !vec_close(tc.d[2], full.d[4], 1e-4) || !vec_close(tc.d[4], full.d[7], 1e-4) ||
			    !vec_close(tc.d[5], full.d[8], 1e-4))
				ok_c = 0;
		}
	}
	qaws_surface_destroy(sp);
	qaws_surface_destroy(sm);
	printf("    %s: parameter tangent matches finite differences: %s\n", f->name, ok_p ? "yes" : "NO");
	TEST_ASSERT(ok_p, "surface parameter tangent matches finite differences");
	TEST_ASSERT(ok_c, "surface coordinate tangent equals next partials");
}

static qaws_status run_surface_adjoint(surface_family const* f, qaws_surface const* s,
	qaws_diff_accumulation acc, qaws_scalar const* us, qaws_scalar const* vs,
	qaws_surface_jet const* ybar, qaws_scalar* pbar, qaws_scalar* ubar, qaws_scalar* vbar)
{
	qaws_diff_context ctx;
	qaws_field_view fv[2];
	qaws_diff_views views;
	qaws_diff_context_init(&ctx);
	ctx.accumulation = acc;
	ctx.tile_size = 3;
	memset(pbar, 0, sizeof(qaws_scalar) * PARAMS(f));
	memset(ubar, 0, sizeof(qaws_scalar) * SURF_SAMPLES);
	memset(vbar, 0, sizeof(qaws_scalar) * SURF_SAMPLES);
	cp_views(f, pbar, fv, &views);
	return qaws_surface_eval_batch_adjoint(&ctx, s, us, vs, SURF_SAMPLES, QAWS_SJET_ORDER3,
		ybar, &views, ubar, vbar);
}

static void check_adjoint(surface_family const* f, qaws_surface const* s)
{
	qaws_scalar us[SURF_SAMPLES], vs[SURF_SAMPLES], udot[SURF_SAMPLES], vdot[SURF_SAMPLES];
	qaws_scalar ubar[SURF_SAMPLES], vbar[SURF_SAMPLES];
	qaws_scalar dir[SURF_MAX_CP * 4], pbar[SURF_MAX_CP * 4];
	qaws_surface_jet ybar[SURF_SAMPLES], tan[SURF_SAMPLES];
	qaws_field_view fv[2];
	qaws_diff_views views;
	unsigned int i;
	double lhs = 0, rhs;

	for (i = 0; i < SURF_SAMPLES; i++)
	{
		sample_uv(i, &us[i], &vs[i]);
		udot[i] = diff_rand();
		vdot[i] = diff_rand();
		rand_sjet(&ybar[i]);
	}
	diff_rand_fill(dir, PARAMS(f));
	cp_views(f, dir, fv, &views);
	TEST_ASSERT_STATUS(qaws_surface_eval_batch_tangent(NULL, s, us, vs, udot, vdot, SURF_SAMPLES,
		QAWS_SJET_ORDER3, &views, NULL, tan));
	for (i = 0; i < SURF_SAMPLES; i++)
		lhs += sjet_dot(&ybar[i], &tan[i], QAWS_SJET_ORDER3);

	TEST_ASSERT_STATUS(run_surface_adjoint(f, s, QAWS_ACCUMULATE_SCATTER, us, vs, ybar, pbar, ubar, vbar));
	rhs = diff_dot(pbar, dir, PARAMS(f)) + diff_dot(ubar, udot, SURF_SAMPLES) + diff_dot(vbar, vdot, SURF_SAMPLES);
	printf("    %s: <ybar, J xdot> = %.9g  <J^T ybar, xdot> = %.9g\n", f->name, lhs, rhs);
	TEST_ASSERT(diff_close(lhs, rhs, DIFF_TOL), "surface adjoint identity over a batch");

	{
		qaws_scalar pbar_t[SURF_MAX_CP * 4], pbar_g[SURF_MAX_CP * 4];
		qaws_scalar ubar2[SURF_SAMPLES], vbar2[SURF_SAMPLES];
		int ok = 1;
		unsigned int n;
		TEST_ASSERT_STATUS(run_surface_adjoint(f, s, QAWS_ACCUMULATE_TILED, us, vs, ybar, pbar_t, ubar2, vbar2));
		TEST_ASSERT_STATUS(run_surface_adjoint(f, s, QAWS_ACCUMULATE_GATHER, us, vs, ybar, pbar_g, ubar2, vbar2));
		for (n = 0; n < PARAMS(f); n++)
			if (!diff_close(pbar[n], pbar_t[n], DIFF_TOL) || !diff_close(pbar[n], pbar_g[n], DIFF_TOL))
				ok = 0;
		for (n = 0; n < SURF_SAMPLES; n++)
			if (!diff_close(ubar[n], ubar2[n], DIFF_TOL) || !diff_close(vbar[n], vbar2[n], DIFF_TOL))
				ok = 0;
		TEST_ASSERT(ok, "surface scatter, tiled and gather accumulation agree");
	}
}

static void check_tangent2(surface_family const* f, qaws_surface const* s)
{
	qaws_scalar dir[SURF_MAX_CP * 4];
	qaws_field_view fv[2];
	qaws_diff_views views;
	qaws_surface *sp, *sm;
	unsigned int i, ch;
	int ok = 1;
	double h = DIFF_FD_STEP;
	qaws_scalar udot = (qaws_scalar)0.6, vdot = (qaws_scalar)-0.35;

	diff_rand_fill(dir, PARAMS(f));
	cp_views(f, dir, fv, &views);
	sp = family_surface_shifted(f, dir, h);
	sm = family_surface_shifted(f, dir, -h);
	for (i = 0; i < SURF_SAMPLES; i++)
	{
		qaws_scalar u, v;
		qaws_surface_jet p, t, tt, tp, tm;
		sample_uv(i, &u, &v);
		TEST_ASSERT_STATUS(qaws_surface_eval_batch_tangent2(NULL, s, &u, &v, &udot, &vdot, 1,
			QAWS_SJET_ORDER3, &views, &p, &t, &tt));
		qaws_surface_eval_tangent(NULL, sp, (qaws_scalar)(u + h * udot), (qaws_scalar)(v + h * vdot),
			udot, vdot, QAWS_SJET_ORDER3, &views, &p, &tp);
		qaws_surface_eval_tangent(NULL, sm, (qaws_scalar)(u - h * udot), (qaws_scalar)(v - h * vdot),
			udot, vdot, QAWS_SJET_ORDER3, &views, &p, &tm);
		for (ch = 0; ch < QAWS_SURFACE_JET_COUNT; ch++)
		{
			qaws_vec3 fd;
			fd.x = (qaws_scalar)((tp.d[ch].x - tm.d[ch].x) / (2 * h));
			fd.y = (qaws_scalar)((tp.d[ch].y - tm.d[ch].y) / (2 * h));
			fd.z = (qaws_scalar)((tp.d[ch].z - tm.d[ch].z) / (2 * h));
			if (!vec_close(fd, tt.d[ch], DIFF_TOL * 100))
			{
				if (ok)
					printf("    %s: sample %u channel %u fd (%.9g %.9g %.9g) tangent2 (%.9g %.9g %.9g)\n", f->name, i, ch,
						fd.x, fd.y, fd.z, tt.d[ch].x, tt.d[ch].y, tt.d[ch].z);
				ok = 0;
			}
		}
	}
	qaws_surface_destroy(sp);
	qaws_surface_destroy(sm);
	printf("    %s: second tangent matches finite differences: %s\n", f->name, ok ? "yes" : "NO");
	TEST_ASSERT(ok, "surface second tangent matches finite differences");
}

/* N = normalize(Su x Sv): surface rule composed with dual kernels. */
static void check_normal_composition(surface_family const* f, qaws_surface const* s)
{
	qaws_scalar dir[SURF_MAX_CP * 4], pbar[SURF_MAX_CP * 4];
	qaws_field_view fv[2];
	qaws_diff_views views;
	unsigned int i;
	int ok = 1;

	for (i = 0; i < SURF_SAMPLES; i++)
	{
		qaws_scalar u, v, udot = diff_rand(), vdot = diff_rand(), ubar = 0, vbar = 0;
		qaws_surface_jet p, t, ybar;
		qaws_dual3 su, sv, n;
		qaws_vec3 nbar = diff_rand_vec3();
		qaws_vec3_pair g;
		double lhs, rhs;

		sample_uv(i, &u, &v);
		diff_rand_fill(dir, PARAMS(f));
		cp_views(f, dir, fv, &views);
		qaws_surface_eval_tangent(NULL, s, u, v, udot, vdot, QAWS_SJET_U | QAWS_SJET_V, &views, &p, &t);
		su = qaws_dual3_make(p.d[1], t.d[1], qaws_v3_zero());
		sv = qaws_dual3_make(p.d[2], t.d[2], qaws_v3_zero());
		n = qaws_dual3_normalize(qaws_dual3_cross(su, sv));
		lhs = qaws_v3_dot(nbar, n.t);

		memset(&ybar, 0, sizeof(ybar));
		g = qaws_cross_adjoint(p.d[1], p.d[2], qaws_normalize_adjoint(qaws_v3_cross(p.d[1], p.d[2]), nbar));
		ybar.d[1] = g.a;
		ybar.d[2] = g.b;
		memset(pbar, 0, sizeof(pbar));
		cp_views(f, pbar, fv, &views);
		qaws_surface_eval_adjoint(NULL, s, u, v, QAWS_SJET_U | QAWS_SJET_V, &ybar, &views, &ubar, &vbar);
		rhs = diff_dot(pbar, dir, PARAMS(f)) + (double)ubar * udot + (double)vbar * vdot;
		if (!diff_close(lhs, rhs, DIFF_TOL))
			ok = 0;
	}
	TEST_ASSERT(ok, "unit normal composition satisfies the adjoint identity");
}

static void test_surface_families(void)
{
	surface_family fams[5];
	unsigned int i;
	make_surface_families(fams);
	for (i = 0; i < 5; i++)
	{
		qaws_surface* s = family_surface(&fams[i], fams[i].cps);
		TEST_ASSERT(s != NULL, "surface created");
		if (!s)
			continue;
		TEST_ASSERT(((qaws_surface_get_diff_capabilities(s) & QAWS_CAP_LINEAR) != 0) == !fams[i].rational,
			"linear flag only for polynomial surfaces");
		TEST_ASSERT((qaws_surface_get_diff_capabilities(s) & QAWS_CAP_TANGENT2) != 0, "second tangent capability");
		check_primal(&fams[i], s);
		check_tangents(&fams[i], s);
		check_adjoint(&fams[i], s);
		check_tangent2(&fams[i], s);
		check_normal_composition(&fams[i], s);
		qaws_surface_destroy(s);
	}
}

static void test_bspline_surface_details(void)
{
	surface_family fams[5];
	qaws_surface* s;
	qaws_diff_context ctx;
	qaws_diff_report report;
	qaws_surface_jet p, t;
	qaws_field_desc fields[4];
	unsigned int n = 0;

	make_surface_families(fams);
	s = family_surface(&fams[1], fams[1].cps);

	qaws_diff_context_init(&ctx);
	qaws_diff_report_reset(&report);
	ctx.report = &report;
	qaws_surface_eval_tangent(&ctx, s, (qaws_scalar)0.5, (qaws_scalar)0.4, 0, 0, QAWS_SJET_ORDER1, NULL, &p, &t);
	TEST_ASSERT(report.validity == QAWS_DIFF_VALID, "interior sample is valid");
	qaws_surface_eval_tangent(&ctx, s, fams[1].u_knots[4], (qaws_scalar)0.4, 1, 0, QAWS_SJET_ORDER1, NULL, &p, &t);
	TEST_ASSERT(report.validity == QAWS_DIFF_AT_BOUNDARY, "coordinate tangent on a knot line is one-sided");
	TEST_ASSERT((report.frozen_used & QAWS_FREEZE_SPAN) != 0, "knot span recorded as frozen");

	TEST_ASSERT_STATUS(qaws_surface_describe_fields(s, fields, 4, &n));
	TEST_ASSERT(n == 3 && fields[1].field == QAWS_FIELD_U_KNOTS && fields[2].field == QAWS_FIELD_V_KNOTS,
		"B-spline surface schema lists both knot vectors");

	{
		qaws_scalar us[4] = { 0.1f, 0.5f, 0.7f, 0.95f };
		qaws_scalar vs[4] = { 0.2f, 0.1f, 0.8f, 0.5f };
		unsigned int offsets[31], samples[64], entries = 0;
		TEST_ASSERT_STATUS(qaws_surface_build_support_index(s, us, vs, 4, QAWS_FIELD_CONTROL_POINTS,
			offsets, 31, samples, 64, &entries));
		TEST_ASSERT(entries == 4 * 12, "support index has (p+1)(q+1) entries per sample");
		TEST_ASSERT(offsets[30] == entries, "support index offsets are cumulative");
	}
	qaws_surface_destroy(s);
}

int test_51_diff_surfaces_main(void)
{
	g_pass = 0;
	g_fail = 0;

	printf("Test 51: Differentiable surface evaluation\n");
	test_surface_families();
	test_bspline_surface_details();

	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
