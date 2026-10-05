/*
 * Test 55: Differential maps of geometry-building operations
 *
 * For each operation (split, join, conversions, degree change) and input
 * family (Bezier, Hermite, B-spline, NURBS):
 *   - the map tangent matches finite differences of the operation applied
 *     to rebuilt inputs (qaws_curve_clone_with_fields)
 *   - adjoint identity <ybar, J xdot> = <J^T ybar, xdot>
 * plus the Bezier split parameter and the refusal for other families.
 */

#include "test_diff.h"

#define MAP_FIELDS 4
#define MAP_OBJECTS 3
#define MAP_MAX_SCALARS 256

/* Storage for every differentiable field of one object. */
typedef struct object_storage
{
	unsigned int field_count;
	qaws_diff_field fields[MAP_FIELDS];
	unsigned int counts[MAP_FIELDS], comps[MAP_FIELDS];
	qaws_scalar data[MAP_FIELDS][MAP_MAX_SCALARS];
	qaws_field_view views[MAP_FIELDS];
	qaws_diff_views dv;
} object_storage;

static void storage_layout(object_storage* s, qaws_curve const* c)
{
	qaws_field_desc fd[8];
	unsigned int n = 0, i;
	memset(s, 0, sizeof(*s));
	qaws_curve_describe_fields(c, fd, 8, &n);
	for (i = 0; i < n && s->field_count < MAP_FIELDS; i++)
	{
		if (!fd[i].capabilities)
			continue;
		s->fields[s->field_count] = fd[i].field;
		s->counts[s->field_count] = fd[i].count;
		s->comps[s->field_count] = (unsigned int)fd[i].value_type;
		s->field_count++;
	}
}

static void storage_bind(object_storage* s)
{
	unsigned int f;
	for (f = 0; f < s->field_count; f++)
		s->views[f] = qaws_field_view_make(s->fields[f], s->data[f], s->counts[f], s->comps[f]);
	s->dv.fields = s->views;
	s->dv.field_count = s->field_count;
	s->dv.children = NULL;
	s->dv.child_count = 0;
}

static unsigned int storage_size(object_storage const* s, unsigned int f)
{
	return s->counts[f] * s->comps[f];
}

static double storage_dot(object_storage const* a, object_storage const* b)
{
	double r = 0;
	unsigned int f;
	for (f = 0; f < a->field_count; f++)
		r += diff_dot(a->data[f], b->data[f], storage_size(a, f));
	return r;
}

/* Primal values of an object's fields into storage. */
static void storage_read(object_storage* s, qaws_curve const* c)
{
	unsigned int f, n;
	for (f = 0; f < s->field_count; f++)
		qaws_curve_read_field(c, s->fields[f], s->data[f], MAP_MAX_SCALARS, &n);
}

/* Curve with values + h * dir. */
static qaws_curve* shifted(qaws_curve const* c, object_storage const* dir, double h)
{
	object_storage v;
	unsigned int f, k;
	qaws_curve* out = NULL;
	storage_layout(&v, c);
	storage_read(&v, c);
	for (f = 0; f < v.field_count; f++)
		for (k = 0; k < storage_size(&v, f); k++)
			v.data[f][k] = (qaws_scalar)(v.data[f][k] + h * dir->data[f][k]);
	storage_bind(&v);
	qaws_curve_clone_with_fields(c, &v.dv, &out);
	return out;
}

typedef qaws_status (*op_fn)(qaws_curve const* const* in, qaws_scalar param, qaws_curve** out, qaws_diff_map** map);

static qaws_status op_split(qaws_curve const* const* in, qaws_scalar p, qaws_curve** out, qaws_diff_map** map)
{
	return qaws_curve_split_diff(in[0], p, &out[0], &out[1], map);
}
static qaws_status op_join(qaws_curve const* const* in, qaws_scalar p, qaws_curve** out, qaws_diff_map** map)
{
	(void)p;
	return qaws_curve_join_diff(in[0], in[1], &out[0], map);
}
static qaws_status op_hermite_to_bezier(qaws_curve const* const* in, qaws_scalar p, qaws_curve** out, qaws_diff_map** map)
{
	(void)p;
	return qaws_curve_convert_hermite_to_bezier_diff(in[0], 1, &out[0], map);
}
static qaws_status op_bezier_to_bspline(qaws_curve const* const* in, qaws_scalar p, qaws_curve** out, qaws_diff_map** map)
{
	(void)p;
	return qaws_curve_convert_bezier_to_bspline_diff(in[0], &out[0], map);
}
static qaws_status op_bspline_to_nurbs(qaws_curve const* const* in, qaws_scalar p, qaws_curve** out, qaws_diff_map** map)
{
	(void)p;
	return qaws_curve_convert_bspline_to_nurbs_diff(in[0], &out[0], map);
}
static qaws_status op_elevate(qaws_curve const* const* in, qaws_scalar p, qaws_curve** out, qaws_diff_map** map)
{
	(void)p;
	return qaws_curve_elevate_degree_diff(in[0], &out[0], map);
}
static qaws_status op_reduce(qaws_curve const* const* in, qaws_scalar p, qaws_curve** out, qaws_diff_map** map)
{
	(void)p;
	return qaws_curve_reduce_degree_diff(in[0], &out[0], map);
}

/* Checks one operation: tangent against finite differences and the
   adjoint identity. param_tangent != 0 also moves the operation parameter
   (input object n_in). */
static void check_op(char const* name, op_fn op, qaws_curve const* const* in, unsigned int n_in,
	unsigned int n_out, qaws_scalar param, qaws_scalar param_tangent)
{
	qaws_curve* out[MAP_OBJECTS] = { NULL, NULL, NULL };
	qaws_curve* outp[MAP_OBJECTS] = { NULL, NULL, NULL };
	qaws_curve* outm[MAP_OBJECTS] = { NULL, NULL, NULL };
	qaws_curve const* inp[MAP_OBJECTS];
	qaws_curve const* inm[MAP_OBJECTS];
	qaws_curve* shifted_p[MAP_OBJECTS];
	qaws_curve* shifted_m[MAP_OBJECTS];
	qaws_diff_map* map = NULL;
	object_storage dir[MAP_OBJECTS + 1], tan[MAP_OBJECTS], ybar[MAP_OBJECTS], bar[MAP_OBJECTS + 1];
	object_storage vp, vm;
	qaws_diff_views const* in_views[MAP_OBJECTS + 1];
	qaws_diff_views* out_views[MAP_OBJECTS];
	qaws_diff_views const* out_adj[MAP_OBJECTS];
	qaws_diff_views* in_adj[MAP_OBJECTS + 1];
	qaws_scalar pdot_storage = param_tangent, pbar_storage = 0;
	qaws_field_view pview, pbar_view;
	qaws_diff_views pviews, pbar_views;
	unsigned int i, f, k, n_in_views = n_in + (param_tangent != 0 ? 1u : 0u);
	double h = DIFF_FD_STEP, lhs = 0, rhs = 0;
	int ok_fd = 1;
	qaws_status st;

	st = op(in, param, out, &map);
	TEST_ASSERT_STATUS(st);
	if (st != QAWS_STATUS_OK)
		return;

	/* Random input tangents. */
	for (i = 0; i < n_in; i++)
	{
		storage_layout(&dir[i], in[i]);
		for (f = 0; f < dir[i].field_count; f++)
			diff_rand_fill(dir[i].data[f], storage_size(&dir[i], f));
		storage_bind(&dir[i]);
		in_views[i] = &dir[i].dv;
	}
	pview = qaws_field_view_make(QAWS_FIELD_PARAMETER, &pdot_storage, 1, 1);
	pviews.fields = &pview;
	pviews.field_count = 1;
	pviews.children = NULL;
	pviews.child_count = 0;
	if (param_tangent != 0)
		in_views[n_in] = &pviews;

	for (i = 0; i < n_out; i++)
	{
		storage_layout(&tan[i], out[i]);
		storage_bind(&tan[i]);
		out_views[i] = &tan[i].dv;
	}
	TEST_ASSERT_STATUS(qaws_diff_map_tangent(map, NULL, in_views, n_in_views, out_views, n_out));

	/* Finite differences of the operation on rebuilt inputs. */
	for (i = 0; i < n_in; i++)
	{
		shifted_p[i] = shifted(in[i], &dir[i], h);
		shifted_m[i] = shifted(in[i], &dir[i], -h);
		inp[i] = shifted_p[i];
		inm[i] = shifted_m[i];
	}
	op(inp, (qaws_scalar)(param + h * param_tangent), outp, NULL);
	op(inm, (qaws_scalar)(param - h * param_tangent), outm, NULL);
	for (i = 0; i < n_out; i++)
	{
		storage_layout(&vp, out[i]);
		storage_layout(&vm, out[i]);
		storage_read(&vp, outp[i]);
		storage_read(&vm, outm[i]);
		for (f = 0; f < vp.field_count; f++)
			for (k = 0; k < storage_size(&vp, f); k++)
			{
				double fd = (vp.data[f][k] - vm.data[f][k]) / (2 * h);
				if (!diff_close(fd, tan[i].data[f][k], DIFF_TOL * 300))
				{
					if (ok_fd)
						printf("    %s: out %u field %s[%u] fd %.9g map %.9g\n", name, i,
							qaws_diff_field_name(vp.fields[f]), k, fd, (double)tan[i].data[f][k]);
					ok_fd = 0;
				}
			}
	}
	for (i = 0; i < n_in; i++)
	{
		qaws_curve_destroy(shifted_p[i]);
		qaws_curve_destroy(shifted_m[i]);
	}
	for (i = 0; i < n_out; i++)
	{
		qaws_curve_destroy(outp[i]);
		qaws_curve_destroy(outm[i]);
	}

	/* Adjoint identity. */
	for (i = 0; i < n_out; i++)
	{
		storage_layout(&ybar[i], out[i]);
		for (f = 0; f < ybar[i].field_count; f++)
			diff_rand_fill(ybar[i].data[f], storage_size(&ybar[i], f));
		storage_bind(&ybar[i]);
		out_adj[i] = &ybar[i].dv;
		lhs += storage_dot(&ybar[i], &tan[i]);
	}
	for (i = 0; i < n_in; i++)
	{
		storage_layout(&bar[i], in[i]);
		storage_bind(&bar[i]);
		in_adj[i] = &bar[i].dv;
	}
	pbar_view = qaws_field_view_make(QAWS_FIELD_PARAMETER, &pbar_storage, 1, 1);
	pbar_views.fields = &pbar_view;
	pbar_views.field_count = 1;
	pbar_views.children = NULL;
	pbar_views.child_count = 0;
	in_adj[n_in] = &pbar_views;
	TEST_ASSERT_STATUS(qaws_diff_map_adjoint(map, NULL, out_adj, n_out, in_adj, n_in + 1));
	for (i = 0; i < n_in; i++)
		rhs += storage_dot(&bar[i], &dir[i]);
	rhs += (double)pbar_storage * param_tangent;

	{
		unsigned int entries = 0;
		qaws_diff_map_get_entries(map, NULL, 0, &entries);
		printf("    %s: %u sparse entries, tangent %s, adjoint %.9g vs %.9g\n", name, entries,
			ok_fd ? "ok" : "NO", lhs, rhs);
	}
	TEST_ASSERT(ok_fd, "map tangent matches finite differences of the operation");
	TEST_ASSERT(diff_close(lhs, rhs, DIFF_TOL * 10), "map adjoint identity");

	qaws_diff_map_destroy(map);
	for (i = 0; i < n_out; i++)
		qaws_curve_destroy(out[i]);
}

/* ------------------------------------------------------------------ */

static qaws_scalar const g_map_knots[11] = { 0, 0, 0, 0, 0.6f, 1.3f, 2, 3, 3, 3, 3 };

static qaws_curve* make_bezier(unsigned int degree)
{
	qaws_scalar cps[3 * 8];
	qaws_bezier_desc d;
	qaws_curve* c = NULL;
	diff_rand_fill(cps, 3 * (degree + 1));
	d.dimension = QAWS_DIMENSION_3D;
	d.degree = degree;
	d.control_points = cps;
	d.control_point_count = degree + 1;
	qaws_curve_create_bezier(&d, &c);
	return c;
}

static qaws_curve* make_bspline(void)
{
	qaws_scalar cps[3 * 7];
	qaws_bspline_desc d;
	qaws_curve* c = NULL;
	diff_rand_fill(cps, 21);
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_3D;
	d.degree = 3;
	d.control_points = cps;
	d.control_point_count = 7;
	d.knots = g_map_knots;
	d.knot_count = 11;
	qaws_curve_create_bspline(&d, &c);
	return c;
}

static qaws_curve* make_nurbs(void)
{
	qaws_scalar cps[3 * 7], w[7];
	qaws_nurbs_desc d;
	qaws_curve* c = NULL;
	unsigned int i;
	diff_rand_fill(cps, 21);
	for (i = 0; i < 7; i++)
		w[i] = (qaws_scalar)1.2 + (qaws_scalar)0.5 * diff_rand();
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_3D;
	d.degree = 3;
	d.control_points = cps;
	d.control_point_count = 7;
	d.knots = g_map_knots;
	d.knot_count = 11;
	d.weights = w;
	d.weight_count = 7;
	qaws_curve_create_nurbs(&d, &c);
	return c;
}

static qaws_curve* make_hermite(void)
{
	qaws_scalar pts[3 * 4], ders[3 * 4];
	qaws_hermite_desc d;
	qaws_curve* c = NULL;
	diff_rand_fill(pts, 12);
	diff_rand_fill(ders, 12);
	d.dimension = QAWS_DIMENSION_3D;
	d.degree = 3;
	d.points = pts;
	d.derivatives = ders;
	d.point_count = 4;
	d.derivative_count = 4;
	qaws_curve_create_hermite(&d, &c);
	return c;
}

static void test_clone(void)
{
	qaws_curve* c = make_nurbs();
	qaws_curve* clone = NULL;
	qaws_scalar w[7] = { 1, 2, 1, 2, 1, 2, 1 };
	qaws_scalar got[7];
	unsigned int n = 0;
	qaws_field_view v = qaws_field_view_make(QAWS_FIELD_WEIGHTS, w, 7, 1);
	qaws_diff_views vs;
	vs.fields = &v;
	vs.field_count = 1;
	vs.children = NULL;
	vs.child_count = 0;
	TEST_ASSERT_STATUS(qaws_curve_clone_with_fields(c, &vs, &clone));
	qaws_curve_read_field(clone, QAWS_FIELD_WEIGHTS, got, 7, &n);
	TEST_ASSERT(n == 7 && got[1] == 2 && got[2] == 1, "clone replaces the given field");
	{
		qaws_scalar a[21], b[21];
		qaws_curve_read_field(c, QAWS_FIELD_CONTROL_POINTS, a, 21, &n);
		qaws_curve_read_field(clone, QAWS_FIELD_CONTROL_POINTS, b, 21, &n);
		TEST_ASSERT(memcmp(a, b, sizeof(a)) == 0, "clone keeps absent fields");
	}
	qaws_curve_destroy(clone);
	qaws_curve_destroy(c);
}

int test_55_diff_maps_main(void)
{
	qaws_curve* c[2];

	g_pass = 0;
	g_fail = 0;
	printf("Test 55: Differential maps of geometry-building operations\n");
	diff_seed(5555);
	test_clone();

	c[0] = make_bezier(5);
	check_op("split bezier (with parameter)", op_split, (qaws_curve const* const*)c, 1, 2, (qaws_scalar)0.37, (qaws_scalar)0.8);
	check_op("elevate bezier", op_elevate, (qaws_curve const* const*)c, 1, 1, 0, 0);
	check_op("reduce bezier", op_reduce, (qaws_curve const* const*)c, 1, 1, 0, 0);
	check_op("bezier -> bspline", op_bezier_to_bspline, (qaws_curve const* const*)c, 1, 1, 0, 0);
	qaws_curve_destroy(c[0]);

	c[0] = make_bspline();
	check_op("split bspline", op_split, (qaws_curve const* const*)c, 1, 2, (qaws_scalar)1.7, 0);
	check_op("bspline -> nurbs", op_bspline_to_nurbs, (qaws_curve const* const*)c, 1, 1, 0, 0);
	c[1] = make_bspline();
	check_op("join bspline", op_join, (qaws_curve const* const*)c, 2, 1, 0, 0);
	qaws_curve_destroy(c[1]);
	{
		/* a non-Bezier split refuses a parameter tangent */
		qaws_curve *l = NULL, *r = NULL;
		qaws_diff_map* map = NULL;
		qaws_scalar pdot = 1;
		qaws_field_view pv = qaws_field_view_make(QAWS_FIELD_PARAMETER, &pdot, 1, 1);
		qaws_diff_views pvs, empty;
		qaws_diff_views const* ins[2];
		qaws_diff_views* outs[2];
		pvs.fields = &pv; pvs.field_count = 1; pvs.children = NULL; pvs.child_count = 0;
		empty.fields = NULL; empty.field_count = 0; empty.children = NULL; empty.child_count = 0;
		qaws_curve_split_diff(c[0], (qaws_scalar)1.7, &l, &r, &map);
		ins[0] = NULL; ins[1] = &pvs;
		outs[0] = &empty; outs[1] = &empty;
		TEST_ASSERT(qaws_diff_map_tangent(map, NULL, ins, 2, outs, 2) == QAWS_STATUS_UNSUPPORTED_OPERATION,
			"B-spline split parameter tangent is refused");
		qaws_diff_map_destroy(map);
		qaws_curve_destroy(l);
		qaws_curve_destroy(r);
	}
	qaws_curve_destroy(c[0]);

	c[0] = make_nurbs();
	printf("    nurbs block\n");
	check_op("split nurbs (homogeneous)", op_split, (qaws_curve const* const*)c, 1, 2, (qaws_scalar)1.1, 0);
	c[1] = make_nurbs();
	check_op("join nurbs", op_join, (qaws_curve const* const*)c, 2, 1, 0, 0);
	qaws_curve_destroy(c[1]);
	qaws_curve_destroy(c[0]);

	c[0] = make_hermite();
	check_op("split hermite", op_split, (qaws_curve const* const*)c, 1, 2, (qaws_scalar)1.4, 0);
	check_op("hermite -> bezier", op_hermite_to_bezier, (qaws_curve const* const*)c, 1, 1, 0, 0);
	qaws_curve_destroy(c[0]);

	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
