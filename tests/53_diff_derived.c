/*
 * Test 53: Differentiable derived surfaces
 *
 * Derived surfaces chain their rules into the curves they reference.
 * Tangent and adjoint storage mirror that structure: the surface's own
 * fields in views->fields, child i in views->children[i].
 *
 * For each derived kind:
 *   - jets match qaws_surface_evaluate
 *   - tangents (own fields, children, u/v) match finite differences of
 *     rebuilt surfaces
 *   - adjoint identity across own fields, children and coordinates
 *   - second tangent matches finite differences of the tangent
 * Plus a composition through the surface normal and the child API.
 */

#include "test_diff.h"

#define DER_SAMPLES 5
#define DER_MAX 32
#define DER_FIXTURES 5

/* A derived fixture: up to four child curves plus up to two own fields. */
typedef struct derived_fixture
{
	char const* name;
	int kind; /* 0 extrusion (2D profile), 1 extrusion (3D profile), 2 ruled, 3 revolution, 4 coons */
	unsigned int child_count;
	unsigned int child_dim[4];
	unsigned int child_cp[4];
	qaws_scalar child_params[4][DER_MAX * 3];
	unsigned int own_count;      /* total scalars of the own fields */
	unsigned int own_nf;
	qaws_diff_field own_field[2];
	unsigned int own_comp[2];
	qaws_scalar own[4];
} derived_fixture;

typedef struct derived_instance
{
	qaws_curve* children[4];
	qaws_surface* surface;
} derived_instance;

static qaws_scalar const g_knots_6[10] = { 0, 0, 0, 0, 0.4f, 0.75f, 1.5f, 1.5f, 1.5f, 1.5f };

static qaws_curve* make_child(unsigned int dim, unsigned int cp_count, qaws_scalar const* p)
{
	qaws_curve* c = NULL;
	if (dim == 2)
	{
		qaws_bezier_desc d;
		d.dimension = QAWS_DIMENSION_2D;
		d.degree = cp_count - 1;
		d.control_points = p;
		d.control_point_count = cp_count;
		qaws_curve_create_bezier(&d, &c);
	}
	else
	{
		qaws_bspline_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_3D;
		d.degree = 3;
		d.control_points = p;
		d.control_point_count = cp_count;
		d.knots = g_knots_6;
		d.knot_count = cp_count + 4;
		qaws_curve_create_bspline(&d, &c);
	}
	return c;
}

static int instance_build(derived_fixture const* f, qaws_scalar const (*child_params)[DER_MAX * 3],
	qaws_scalar const* own, derived_instance* out)
{
	unsigned int i;
	memset(out, 0, sizeof(*out));
	for (i = 0; i < f->child_count; i++)
	{
		out->children[i] = make_child(f->child_dim[i], f->child_cp[i], child_params[i]);
		if (!out->children[i])
			return 0;
	}
	if (f->kind <= 1)
	{
		qaws_surface_extrusion_desc d;
		d.profile = out->children[0];
		d.direction = qaws_v3(own[0], own[1], own[2]);
		d.length = 0;
		qaws_surface_create_extrusion(&d, &out->surface);
	}
	else if (f->kind == 2)
	{
		qaws_surface_ruled_desc d;
		memset(&d, 0, sizeof(d));
		d.curve_a = out->children[0];
		d.curve_b = out->children[1];
		qaws_surface_create_ruled(&d, &out->surface);
	}
	else if (f->kind == 3)
	{
		qaws_surface_revolution_desc d;
		memset(&d, 0, sizeof(d));
		d.profile = out->children[0];
		d.axis_origin = qaws_v3(own[0], own[1], own[2]);
		d.axis_direction = qaws_v3(0, 0, 1);
		d.angle = own[3];
		qaws_surface_create_revolution(&d, &out->surface);
	}
	else
	{
		qaws_surface_coons_desc d;
		memset(&d, 0, sizeof(d));
		d.c0 = out->children[0];
		d.c1 = out->children[1];
		d.d0 = out->children[2];
		d.d1 = out->children[3];
		qaws_surface_create_coons(&d, &out->surface);
	}
	return out->surface != NULL;
}

static void instance_destroy(derived_instance* in)
{
	unsigned int i;
	if (in->surface)
		qaws_surface_destroy(in->surface);
	for (i = 0; i < 4; i++)
		if (in->children[i])
			qaws_curve_destroy(in->children[i]);
}

/* Direction tangent / adjoint storage plus child storage. */
typedef struct derived_storage
{
	qaws_scalar child[4][DER_MAX * 3];
	qaws_scalar own[4];
	qaws_field_view own_view[2];
	qaws_field_view child_view[4];
	qaws_diff_views child_views[4];
	qaws_diff_views views;
} derived_storage;

static void storage_bind(derived_fixture const* f, derived_storage* s)
{
	unsigned int i;
	s->own_view[0] = qaws_field_view_make(f->own_field[0], s->own, 1, f->own_comp[0]);
	s->own_view[1] = qaws_field_view_make(f->own_field[1], s->own + f->own_comp[0], 1, f->own_comp[1]);
	for (i = 0; i < f->child_count; i++)
	{
		s->child_view[i] = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, s->child[i], f->child_cp[i], f->child_dim[i]);
		s->child_views[i].fields = &s->child_view[i];
		s->child_views[i].field_count = 1;
		s->child_views[i].children = NULL;
		s->child_views[i].child_count = 0;
	}
	s->views.fields = s->own_view;
	s->views.field_count = f->own_nf;
	s->views.children = s->child_views;
	s->views.child_count = f->child_count;
}

static void storage_random(derived_fixture const* f, derived_storage* s)
{
	unsigned int i;
	memset(s, 0, sizeof(*s));
	for (i = 0; i < f->child_count; i++)
		diff_rand_fill(s->child[i], f->child_cp[i] * f->child_dim[i]);
	diff_rand_fill(s->own, 4);
	storage_bind(f, s);
}

static double storage_dot(derived_fixture const* f, derived_storage const* a, derived_storage const* b)
{
	double r = 0;
	unsigned int i;
	for (i = 0; i < f->child_count; i++)
		r += diff_dot(a->child[i], b->child[i], f->child_cp[i] * f->child_dim[i]);
	r += diff_dot(a->own, b->own, f->own_count);
	return r;
}

static int shifted_build(derived_fixture const* f, derived_storage const* dir, double h, derived_instance* out)
{
	qaws_scalar cp[4][DER_MAX * 3], own[4];
	unsigned int i, n;
	for (i = 0; i < f->child_count; i++)
		for (n = 0; n < f->child_cp[i] * f->child_dim[i]; n++)
			cp[i][n] = (qaws_scalar)(f->child_params[i][n] + h * dir->child[i][n]);
	for (n = 0; n < 4; n++)
		own[n] = (qaws_scalar)(f->own[n] + (n < f->own_count ? h * dir->own[n] : 0));
	return instance_build(f, (qaws_scalar const (*)[DER_MAX * 3])cp, own, out);
}

static int vclose(qaws_vec3 a, qaws_vec3 b, double tol)
{
	double ma = fabs(a.x), mb = fabs(b.x), scale = 1;
	if (fabs(a.y) > ma) ma = fabs(a.y);
	if (fabs(a.z) > ma) ma = fabs(a.z);
	if (fabs(b.y) > mb) mb = fabs(b.y);
	if (fabs(b.z) > mb) mb = fabs(b.z);
	if (ma > scale) scale = ma;
	if (mb > scale) scale = mb;
	return fabs(a.x - b.x) <= tol * scale && fabs(a.y - b.y) <= tol * scale && fabs(a.z - b.z) <= tol * scale;
}

static void make_fixtures(derived_fixture* fx)
{
	unsigned int i;
	memset(fx, 0, sizeof(derived_fixture) * DER_FIXTURES);
	diff_seed(909);

	fx[0].name = "extrusion_2d";
	fx[0].kind = 0;
	fx[0].child_count = 1;
	fx[0].child_dim[0] = 2;
	fx[0].child_cp[0] = 5;

	fx[1].name = "extrusion_3d";
	fx[1].kind = 1;
	fx[1].child_count = 1;
	fx[1].child_dim[0] = 3;
	fx[1].child_cp[0] = 6;

	fx[2].name = "ruled";
	fx[2].kind = 2;
	fx[2].child_count = 2;
	fx[2].child_dim[0] = 3;
	fx[2].child_dim[1] = 3;
	fx[2].child_cp[0] = 6;
	fx[2].child_cp[1] = 6;

	for (i = 0; i < 3; i++)
	{
		unsigned int c, n;
		for (c = 0; c < fx[i].child_count; c++)
			for (n = 0; n < fx[i].child_cp[c]; n++)
			{
				qaws_scalar* p = &fx[i].child_params[c][n * fx[i].child_dim[c]];
				p[0] = (qaws_scalar)n + diff_rand() * (qaws_scalar)0.3;
				p[1] = (qaws_scalar)c * 2 + diff_rand() * (qaws_scalar)0.6;
				if (fx[i].child_dim[c] == 3)
					p[2] = diff_rand();
			}
		if (fx[i].kind <= 1)
		{
			fx[i].own_nf = 1;
			fx[i].own_field[0] = QAWS_FIELD_DIRECTION;
			fx[i].own_comp[0] = 3;
			fx[i].own_count = 3;
		}
		fx[i].own[0] = (qaws_scalar)0.2;
		fx[i].own[1] = (qaws_scalar)0.3;
		fx[i].own[2] = (qaws_scalar)1.7;
	}

	/* Revolution: 2D profile (radius, height), own center and angle. */
	fx[3].name = "revolution";
	fx[3].kind = 3;
	fx[3].child_count = 1;
	fx[3].child_dim[0] = 2;
	fx[3].child_cp[0] = 5;
	for (i = 0; i < 5; i++)
	{
		fx[3].child_params[0][2 * i] = (qaws_scalar)1.0 + (qaws_scalar)0.4 * diff_rand();
		fx[3].child_params[0][2 * i + 1] = (qaws_scalar)i * (qaws_scalar)0.6;
	}
	fx[3].own_nf = 2;
	fx[3].own_field[0] = QAWS_FIELD_CENTER;
	fx[3].own_comp[0] = 3;
	fx[3].own_field[1] = QAWS_FIELD_ANGLE_END;
	fx[3].own_comp[1] = 1;
	fx[3].own_count = 4;
	fx[3].own[0] = (qaws_scalar)0.1;
	fx[3].own[1] = (qaws_scalar)-0.2;
	fx[3].own[2] = (qaws_scalar)0.3;
	fx[3].own[3] = (qaws_scalar)4.0;

	/* Coons: four boundary curves sharing corners (0,0) (3,0) (0,3) (3,3). */
	fx[4].name = "coons";
	fx[4].kind = 4;
	fx[4].child_count = 4;
	for (i = 0; i < 4; i++)
	{
		unsigned int n;
		fx[4].child_dim[i] = 3;
		fx[4].child_cp[i] = 6;
		for (n = 0; n < 6; n++)
		{
			qaws_scalar s = (qaws_scalar)n * (qaws_scalar)0.6;
			qaws_scalar* p = &fx[4].child_params[i][3 * n];
			qaws_scalar bump = (n == 0 || n == 5) ? QAWS_ZERO : (qaws_scalar)0.4 * diff_rand();
			if (i == 0) { p[0] = s; p[1] = 0; }       /* c0: v = 0 */
			else if (i == 1) { p[0] = s; p[1] = 3; }  /* c1: v = 1 */
			else if (i == 2) { p[0] = 0; p[1] = s; }  /* d0: u = 0 */
			else { p[0] = 3; p[1] = s; }              /* d1: u = 1 */
			p[2] = bump + (i == 1 ? (qaws_scalar)0.5 * (qaws_scalar)n / 5 : QAWS_ZERO)
			     + (i == 3 ? (qaws_scalar)0.5 * (qaws_scalar)n / 5 : QAWS_ZERO);
		}
	}
}

static void sample(unsigned int i, qaws_scalar* u, qaws_scalar* v)
{
	*u = (qaws_scalar)0.12 + (qaws_scalar)0.17 * (qaws_scalar)i;
	*v = (qaws_scalar)0.81 - (qaws_scalar)0.14 * (qaws_scalar)i;
}

static void check_fixture(derived_fixture const* f)
{
	derived_instance base;
	unsigned int i, ch;
	int ok_jet = 1, ok_fd = 1, ok_adj = 1, ok_t2 = 1;
	double h = DIFF_FD_STEP;

	if (!instance_build(f, (qaws_scalar const (*)[DER_MAX * 3])f->child_params, f->own, &base))
	{
		TEST_ASSERT(0, "derived surface built");
		return;
	}

	for (i = 0; i < DER_SAMPLES; i++)
	{
		qaws_scalar u, v, udot = diff_rand(), vdot = diff_rand(), ubar = 0, vbar = 0;
		derived_storage dir, bar;
		derived_instance sp, sm;
		qaws_surface_jet p, t, tt, jp, jm, ybar;
		qaws_surface_eval_result r;
		double lhs = 0, rhs;

		sample(i, &u, &v);
		qaws_surface_eval_jet(base.surface, u, v, QAWS_SJET_ORDER2, &p);
		qaws_surface_evaluate(base.surface, u, v, QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV |
			QAWS_SURFACE_EVAL_DUU | QAWS_SURFACE_EVAL_DUV | QAWS_SURFACE_EVAL_DVV, &r);
		if (!vclose(p.d[0], r.position, 1e-4) || !vclose(p.d[1], r.du, 1e-4) || !vclose(p.d[2], r.dv, 1e-4) ||
		    !vclose(p.d[3], r.duu, 1e-4) || !vclose(p.d[4], r.duv, 1e-4) || !vclose(p.d[5], r.dvv, 1e-4))
			ok_jet = 0;

		storage_random(f, &dir);
		TEST_ASSERT_STATUS(qaws_surface_eval_batch_tangent2(NULL, base.surface, &u, &v, &udot, &vdot, 1,
			QAWS_SJET_ORDER3, &dir.views, &p, &t, &tt));

		shifted_build(f, &dir, h, &sp);
		shifted_build(f, &dir, -h, &sm);
		qaws_surface_eval_jet(sp.surface, (qaws_scalar)(u + h * udot), (qaws_scalar)(v + h * vdot), QAWS_SJET_ORDER3, &jp);
		qaws_surface_eval_jet(sm.surface, (qaws_scalar)(u - h * udot), (qaws_scalar)(v - h * vdot), QAWS_SJET_ORDER3, &jm);
		for (ch = 0; ch < QAWS_SURFACE_JET_COUNT; ch++)
		{
			qaws_vec3 fd = qaws_v3_scale(qaws_v3_sub(jp.d[ch], jm.d[ch]), (qaws_scalar)(0.5 / h));
			if (!vclose(fd, t.d[ch], DIFF_TOL * 100))
				ok_fd = 0;
		}

		{
			qaws_surface_jet pp, tp, pm, tm;
			qaws_surface_eval_tangent(NULL, sp.surface, (qaws_scalar)(u + h * udot), (qaws_scalar)(v + h * vdot),
				udot, vdot, QAWS_SJET_ORDER3, &dir.views, &pp, &tp);
			qaws_surface_eval_tangent(NULL, sm.surface, (qaws_scalar)(u - h * udot), (qaws_scalar)(v - h * vdot),
				udot, vdot, QAWS_SJET_ORDER3, &dir.views, &pm, &tm);
			for (ch = 0; ch < QAWS_SURFACE_JET_COUNT; ch++)
			{
				qaws_vec3 fd = qaws_v3_scale(qaws_v3_sub(tp.d[ch], tm.d[ch]), (qaws_scalar)(0.5 / h));
				if (!vclose(fd, tt.d[ch], DIFF_TOL * 300))
					ok_t2 = 0;
			}
		}
		instance_destroy(&sp);
		instance_destroy(&sm);

		for (ch = 0; ch < QAWS_SURFACE_JET_COUNT; ch++)
			ybar.d[ch] = diff_rand_vec3();
		ybar.channels = QAWS_SJET_ORDER3;
		for (ch = 0; ch < QAWS_SURFACE_JET_COUNT; ch++)
			lhs += qaws_v3_dot(ybar.d[ch], t.d[ch]);

		memset(&bar, 0, sizeof(bar));
		storage_bind(f, &bar);
		TEST_ASSERT_STATUS(qaws_surface_eval_adjoint(NULL, base.surface, u, v, QAWS_SJET_ORDER3, &ybar,
			&bar.views, &ubar, &vbar));
		rhs = storage_dot(f, &bar, &dir) + (double)ubar * udot + (double)vbar * vdot;
		if (!diff_close(lhs, rhs, DIFF_TOL * 10))
		{
			printf("    %s: adjoint mismatch %.9g vs %.9g\n", f->name, lhs, rhs);
			ok_adj = 0;
		}
	}

	printf("    %s: jets %s, tangent %s, second tangent %s, adjoint %s\n", f->name,
		ok_jet ? "ok" : "NO", ok_fd ? "ok" : "NO", ok_t2 ? "ok" : "NO", ok_adj ? "ok" : "NO");
	TEST_ASSERT(ok_jet, "derived jets match evaluate");
	TEST_ASSERT(ok_fd, "derived tangent matches finite differences");
	TEST_ASSERT(ok_t2, "derived second tangent matches finite differences");
	TEST_ASSERT(ok_adj, "derived adjoint identity across children and own fields");
	instance_destroy(&base);
}

static void test_children_api(void)
{
	derived_fixture fx[DER_FIXTURES];
	derived_instance in;
	qaws_diff_child children[4];
	unsigned int n = 0;
	qaws_field_desc fields[2];
	unsigned int nf = 0;

	make_fixtures(fx);
	instance_build(&fx[2], (qaws_scalar const (*)[DER_MAX * 3])fx[2].child_params, fx[2].own, &in);
	TEST_ASSERT_STATUS(qaws_surface_diff_children(in.surface, children, 4, &n));
	TEST_ASSERT(n == 2 && children[0].curve == in.children[0] && children[1].curve == in.children[1] &&
		!children[0].surface, "ruled surface exposes both curves as children");
	TEST_ASSERT_STATUS(qaws_surface_describe_fields(in.surface, fields, 2, &nf));
	TEST_ASSERT(nf == 0, "ruled surface has no own fields");
	instance_destroy(&in);

	instance_build(&fx[1], (qaws_scalar const (*)[DER_MAX * 3])fx[1].child_params, fx[1].own, &in);
	TEST_ASSERT_STATUS(qaws_surface_describe_fields(in.surface, fields, 2, &nf));
	TEST_ASSERT(nf == 1 && fields[0].field == QAWS_FIELD_DIRECTION && fields[0].domain == QAWS_DOMAIN_DIRECTION,
		"extrusion exposes its direction vector");
	{
		char buf[64];
		qaws_param_key k = qaws_param_key_make(QAWS_FIELD_CONTROL_POINTS, 2, 0, 3);
		qaws_param_key_prepend_child(&k, 0);
		qaws_param_key_to_string(&k, buf, sizeof(buf));
		TEST_ASSERT(strcmp(buf, "child[0]/control_points/2/x") == 0, "child parameter path");
	}
	instance_destroy(&in);
}

/* Profile control points -> extrusion -> unit normal: one chain, one identity. */
static void test_normal_through_extrusion(void)
{
	derived_fixture fx[DER_FIXTURES];
	derived_instance in;
	unsigned int i;
	int ok = 1;

	make_fixtures(fx);
	instance_build(&fx[1], (qaws_scalar const (*)[DER_MAX * 3])fx[1].child_params, fx[1].own, &in);
	for (i = 0; i < DER_SAMPLES; i++)
	{
		qaws_scalar u, v, ubar = 0, vbar = 0;
		derived_storage dir, bar;
		qaws_surface_jet p, t, ybar;
		qaws_surface_geometry g, gt, gbar;
		double lhs, rhs;

		sample(i, &u, &v);
		storage_random(&fx[1], &dir);
		qaws_surface_eval_tangent(NULL, in.surface, u, v, 0, 0, QAWS_SJET_ORDER2, &dir.views, &p, &t);
		qaws_surface_geometry_eval(&p, &t, NULL, &g, &gt, NULL, NULL);
		memset(&gbar, 0, sizeof(gbar));
		gbar.normal = diff_rand_vec3();
		gbar.mean = diff_rand();
		lhs = (double)qaws_v3_dot(gbar.normal, gt.normal) + (double)gbar.mean * gt.mean;

		memset(&ybar, 0, sizeof(ybar));
		qaws_surface_geometry_adjoint(&p, &gbar, &ybar, NULL);
		memset(&bar, 0, sizeof(bar));
		storage_bind(&fx[1], &bar);
		qaws_surface_eval_adjoint(NULL, in.surface, u, v, ybar.channels, &ybar, &bar.views, &ubar, &vbar);
		rhs = storage_dot(&fx[1], &bar, &dir);
		if (!diff_close(lhs, rhs, DIFF_TOL * 10))
			ok = 0;
	}
	TEST_ASSERT(ok, "profile -> extrusion -> normal and mean curvature adjoint identity");
	instance_destroy(&in);
}

/* ------------------------------------------------------------------ */
/*  Offset surface: own distance + a NURBS base surface as child 0    */
/* ------------------------------------------------------------------ */

static qaws_scalar const g_off_knots[8] = { 0, 0, 0, 0, 1, 1, 1, 1 };

typedef struct offset_params
{
	qaws_scalar cps[48];
	qaws_scalar w[16];
	qaws_scalar distance;
} offset_params;

typedef struct offset_instance
{
	qaws_surface* base;
	qaws_surface* offset;
} offset_instance;

static int offset_build(offset_params const* p, offset_instance* out)
{
	qaws_surface_nurbs_desc d;
	qaws_surface_offset_desc od;
	out->base = NULL;
	out->offset = NULL;
	d.u_degree = 3;
	d.v_degree = 3;
	d.control_points = (qaws_vec3 const*)p->cps;
	d.u_point_count = 4;
	d.v_point_count = 4;
	d.weights = p->w;
	d.u_knots = g_off_knots;
	d.u_knot_count = 8;
	d.v_knots = g_off_knots;
	d.v_knot_count = 8;
	if (qaws_surface_create_nurbs(&d, &out->base) != QAWS_STATUS_OK)
		return 0;
	od.base = out->base;
	od.distance = p->distance;
	return qaws_surface_create_offset(&od, &out->offset) == QAWS_STATUS_OK;
}

static void offset_destroy(offset_instance* in)
{
	if (in->offset) qaws_surface_destroy(in->offset);
	if (in->base) qaws_surface_destroy(in->base);
}

typedef struct offset_storage
{
	offset_params data;
	qaws_field_view own_view;
	qaws_field_view base_views[2];
	qaws_diff_views child;
	qaws_diff_views views;
} offset_storage;

static void offset_bind(offset_storage* s)
{
	s->own_view = qaws_field_view_make(QAWS_FIELD_OFFSET_DISTANCE, &s->data.distance, 1, 1);
	s->base_views[0] = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, s->data.cps, 16, 3);
	s->base_views[1] = qaws_field_view_make(QAWS_FIELD_WEIGHTS, s->data.w, 16, 1);
	s->child.fields = s->base_views;
	s->child.field_count = 2;
	s->child.children = NULL;
	s->child.child_count = 0;
	s->views.fields = &s->own_view;
	s->views.field_count = 1;
	s->views.children = &s->child;
	s->views.child_count = 1;
}

static double offset_dot(offset_params const* a, offset_params const* b)
{
	return diff_dot(a->cps, b->cps, 48) + diff_dot(a->w, b->w, 16) + (double)a->distance * b->distance;
}

static void offset_shift(offset_params const* p, offset_params const* dir, double h, offset_params* out)
{
	unsigned int n;
	for (n = 0; n < 48; n++) out->cps[n] = (qaws_scalar)(p->cps[n] + h * dir->cps[n]);
	for (n = 0; n < 16; n++) out->w[n] = (qaws_scalar)(p->w[n] + h * dir->w[n]);
	out->distance = (qaws_scalar)(p->distance + h * dir->distance);
}

static void test_offset(void)
{
	offset_params base;
	offset_instance in;
	unsigned int i, a, b, ch;
	int ok_pos = 1, ok_fd = 1, ok_t2 = 1, ok_adj = 1;
	double h = DIFF_FD_STEP;

	diff_seed(515);
	for (a = 0; a < 4; a++)
		for (b = 0; b < 4; b++)
		{
			qaws_scalar* p = &base.cps[(a * 4 + b) * 3];
			p[0] = (qaws_scalar)a + diff_rand() * (qaws_scalar)0.1;
			p[1] = (qaws_scalar)b + diff_rand() * (qaws_scalar)0.1;
			p[2] = (qaws_scalar)0.3 * ((qaws_scalar)a - (qaws_scalar)1.5) * ((qaws_scalar)b - (qaws_scalar)1.5) + diff_rand() * (qaws_scalar)0.1;
			base.w[a * 4 + b] = (qaws_scalar)1 + (qaws_scalar)0.25 * diff_rand();
		}
	base.distance = (qaws_scalar)0.35;
	offset_build(&base, &in);

	for (i = 0; i < DER_SAMPLES; i++)
	{
		qaws_scalar u, v, udot = diff_rand(), vdot = diff_rand(), ubar = 0, vbar = 0;
		offset_storage dir, bar;
		offset_params pp, pm;
		offset_instance sp, sm;
		qaws_surface_jet p, t, tt, jp, jm, ybar;
		double lhs = 0, rhs;

		sample(i, &u, &v);

		/* Position = base + d * base normal; partials consistent with each other. */
		{
			qaws_surface_eval_result br;
			qaws_surface_jet j0, ju1, ju0;
			qaws_surface_evaluate(in.base, u, v, QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL, &br);
			qaws_surface_eval_jet(in.offset, u, v, QAWS_SJET_ORDER2, &j0);
			if (!vclose(j0.d[0], qaws_v3_axpy(br.position, br.normal, base.distance), 1e-6))
				ok_pos = 0;
			qaws_surface_eval_jet(in.offset, (qaws_scalar)(u + h), v, QAWS_SJET_ORDER2, &ju1);
			qaws_surface_eval_jet(in.offset, (qaws_scalar)(u - h), v, QAWS_SJET_ORDER2, &ju0);
			if (!vclose(qaws_v3_scale(qaws_v3_sub(ju1.d[0], ju0.d[0]), (qaws_scalar)(0.5 / h)), j0.d[1], DIFF_TOL * 100) ||
			    !vclose(qaws_v3_scale(qaws_v3_sub(ju1.d[1], ju0.d[1]), (qaws_scalar)(0.5 / h)), j0.d[3], DIFF_TOL * 100) ||
			    !vclose(qaws_v3_scale(qaws_v3_sub(ju1.d[2], ju0.d[2]), (qaws_scalar)(0.5 / h)), j0.d[4], DIFF_TOL * 100))
				ok_pos = 0;
		}

		memset(&dir, 0, sizeof(dir));
		diff_rand_fill(dir.data.cps, 48);
		diff_rand_fill(dir.data.w, 16);
		dir.data.distance = diff_rand();
		offset_bind(&dir);
		TEST_ASSERT_STATUS(qaws_surface_eval_batch_tangent2(NULL, in.offset, &u, &v, &udot, &vdot, 1,
			QAWS_SJET_ORDER2, &dir.views, &p, &t, &tt));

		offset_shift(&base, &dir.data, h, &pp);
		offset_shift(&base, &dir.data, -h, &pm);
		offset_build(&pp, &sp);
		offset_build(&pm, &sm);
		qaws_surface_eval_jet(sp.offset, (qaws_scalar)(u + h * udot), (qaws_scalar)(v + h * vdot), QAWS_SJET_ORDER2, &jp);
		qaws_surface_eval_jet(sm.offset, (qaws_scalar)(u - h * udot), (qaws_scalar)(v - h * vdot), QAWS_SJET_ORDER2, &jm);
		for (ch = 0; ch < 6; ch++)
			if (!vclose(qaws_v3_scale(qaws_v3_sub(jp.d[ch], jm.d[ch]), (qaws_scalar)(0.5 / h)), t.d[ch], DIFF_TOL * 100))
				ok_fd = 0;
		{
			qaws_surface_jet p1, t1, p2, t2;
			qaws_surface_eval_tangent(NULL, sp.offset, (qaws_scalar)(u + h * udot), (qaws_scalar)(v + h * vdot),
				udot, vdot, QAWS_SJET_ORDER2, &dir.views, &p1, &t1);
			qaws_surface_eval_tangent(NULL, sm.offset, (qaws_scalar)(u - h * udot), (qaws_scalar)(v - h * vdot),
				udot, vdot, QAWS_SJET_ORDER2, &dir.views, &p2, &t2);
			for (ch = 0; ch < 6; ch++)
				if (!vclose(qaws_v3_scale(qaws_v3_sub(t1.d[ch], t2.d[ch]), (qaws_scalar)(0.5 / h)), tt.d[ch], DIFF_TOL * 300))
					ok_t2 = 0;
		}
		offset_destroy(&sp);
		offset_destroy(&sm);

		memset(&ybar, 0, sizeof(ybar));
		for (ch = 0; ch < 6; ch++)
		{
			ybar.d[ch] = diff_rand_vec3();
			lhs += qaws_v3_dot(ybar.d[ch], t.d[ch]);
		}
		ybar.channels = QAWS_SJET_ORDER2;
		memset(&bar, 0, sizeof(bar));
		offset_bind(&bar);
		TEST_ASSERT_STATUS(qaws_surface_eval_adjoint(NULL, in.offset, u, v, QAWS_SJET_ORDER2, &ybar, &bar.views, &ubar, &vbar));
		rhs = offset_dot(&bar.data, &dir.data) + (double)ubar * udot + (double)vbar * vdot;
		if (!diff_close(lhs, rhs, DIFF_TOL * 10))
		{
			printf("    offset: adjoint mismatch %.9g vs %.9g\n", lhs, rhs);
			ok_adj = 0;
		}
	}

	{
		qaws_surface_jet p, t;
		TEST_ASSERT(qaws_surface_eval_tangent(NULL, in.offset, (qaws_scalar)0.5, (qaws_scalar)0.5, 0, 0,
			QAWS_SJET_ORDER3, NULL, &p, &t) == QAWS_STATUS_UNSUPPORTED_OPERATION,
			"offset refuses third-order jets instead of approximating");
	}

	printf("    offset: analytic jets %s, tangent %s, second tangent %s, adjoint %s\n",
		ok_pos ? "ok" : "NO", ok_fd ? "ok" : "NO", ok_t2 ? "ok" : "NO", ok_adj ? "ok" : "NO");
	TEST_ASSERT(ok_pos, "offset analytic jets: position and partial consistency");
	TEST_ASSERT(ok_fd, "offset tangent matches finite differences (distance, base points and weights)");
	TEST_ASSERT(ok_t2, "offset second tangent matches finite differences");
	TEST_ASSERT(ok_adj, "offset adjoint identity through the base surface");
	offset_destroy(&in);
}

int test_53_diff_derived_main(void)
{
	derived_fixture fx[DER_FIXTURES];
	unsigned int i;

	g_pass = 0;
	g_fail = 0;

	printf("Test 53: Differentiable derived surfaces\n");
	make_fixtures(fx);
	for (i = 0; i < DER_FIXTURES; i++)
		check_fixture(&fx[i]);
	test_children_api();
	test_normal_through_extrusion();
	test_offset();

	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
