/*
 * Test 92: the surface of a stack as a 3D triangle mesh (qaws_stack_mesh_3d)
 *
 *   - closed: every directed edge has its reverse exactly once
 *   - outward: the volume by the divergence theorem positive, and the
 *     closed form's: a cylinder 2 pi, a cone pi / 3, a tube with a hole
 *     3 pi, a sphere from 17 circles with point poles 4 pi / 3
 *   - a split (an ellipse into two circles) stays closed
 *   - the same mesh on an executor
 */

#include "test_common.h"
#include <math.h>
#include <string.h>

#define SMT_PI 3.14159265358979323846

typedef struct smt_stack
{
	qaws_curve* curves[64];
	qaws_curve const* views[64];
	qaws_path_2d paths[64];
	qaws_stack_level levels[32];
	unsigned int nc, nl;
	qaws_stack_2d st;
} smt_stack;

static void smt_level(smt_stack* s, double z, qaws_curve* a, qaws_curve* b)
{
	qaws_stack_level* l = &s->levels[s->nl++];
	unsigned int first = s->nc;
	s->curves[s->nc] = a; s->views[s->nc] = a;
	s->paths[s->nc].curves = &s->views[s->nc]; s->paths[s->nc].curve_count = 1; s->paths[s->nc].closed = 1;
	s->nc++;
	if (b)
	{
		s->curves[s->nc] = b; s->views[s->nc] = b;
		s->paths[s->nc].curves = &s->views[s->nc]; s->paths[s->nc].curve_count = 1; s->paths[s->nc].closed = 1;
		s->nc++;
	}
	l->z = (qaws_scalar)z;
	l->paths = &s->paths[first];
	l->path_count = s->nc - first;
}

static void smt_done(smt_stack* s, qaws_stack_interp in, qaws_stack_blend bl)
{
	memset(&s->st, 0, sizeof(s->st));
	s->st.levels = s->levels;
	s->st.level_count = s->nl;
	s->st.fill_rule = QAWS_FILL_NON_ZERO;
	s->st.interp = in;
	s->st.blend = bl;
}

static void smt_free(smt_stack* s)
{
	unsigned int i;
	for (i = 0; i < s->nc; i++) qaws_curve_destroy(s->curves[i]);
}

static qaws_curve* smt_circle(double cx, double cy, double r)
{
	qaws_curve* c = NULL;
	if (r < 1e-12)
	{
		qaws_scalar p[4] = { (qaws_scalar)cx, (qaws_scalar)cy, (qaws_scalar)cx, (qaws_scalar)cy };
		qaws_curve_create_polyline_2d(p, 2, 1, &c);
	}
	else
		qaws_curve_create_ellipse_2d((qaws_scalar)cx, (qaws_scalar)cy, (qaws_scalar)r, (qaws_scalar)r, 0, &c);
	return c;
}

static qaws_curve* smt_hole(double r)
{
	qaws_curve *c = smt_circle(0, 0, r), *rv = NULL;
	qaws_curve_reverse(c, &rv);
	qaws_curve_destroy(c);
	return rv;
}

/* directed edges without their reverse (or doubled), and the volume */
static void smt_check(qaws_mesh_3d const* m, unsigned int* open, double* volume)
{
	unsigned int nt = qaws_mesh_3d_get_triangle_count(m), i, k, n = 0;
	uint64_t* e = (uint64_t*)malloc(sizeof(uint64_t) * 3 * (nt + 1));
	double v = 0;
	for (i = 0; i < nt; i++)
	{
		unsigned int t[3];
		qaws_vec3 a, b, c;
		qaws_mesh_3d_get_triangle(m, i, t);
		a = qaws_mesh_3d_get_vertex(m, t[0]); b = qaws_mesh_3d_get_vertex(m, t[1]); c = qaws_mesh_3d_get_vertex(m, t[2]);
		v += ((double)a.x * ((double)b.y * c.z - (double)b.z * c.y) - (double)a.y * ((double)b.x * c.z - (double)b.z * c.x)
			+ (double)a.z * ((double)b.x * c.y - (double)b.y * c.x)) / 6;
		for (k = 0; k < 3; k++) e[n++] = ((uint64_t)t[k] << 32) | t[(k + 1) % 3];
	}
	/* sort, then every edge must find its reverse once */
	{
		unsigned int a, b;
		for (a = 1; a < n; a++)
		{
			uint64_t key = e[a];
			for (b = a; b > 0 && e[b - 1] > key; b--) e[b] = e[b - 1];
			e[b] = key;
		}
	}
	*open = 0;
	for (i = 0; i < n; i++)
	{
		uint64_t rev = ((e[i] & 0xFFFFFFFFu) << 32) | (e[i] >> 32);
		unsigned int lo = 0, hi = n, cnt = 0;
		while (lo < hi) { unsigned int mid = (lo + hi) / 2; if (e[mid] < rev) lo = mid + 1; else hi = mid; }
		while (lo < n && e[lo] == rev) { cnt++; lo++; }
		if (cnt != 1 || (i > 0 && e[i] == e[i - 1])) (*open)++;
	}
	*volume = v;
	free(e);
}

static qaws_mesh_3d* smt_mesh(qaws_stack_2d const* st, unsigned int sections, qaws_batch_executor const* ex)
{
	qaws_stack_mesh_desc d;
	qaws_mesh_3d* m = NULL;
	memset(&d, 0, sizeof(d));
	d.stack = st; d.sections = sections; d.caps = 1; d.executor = ex;
	qaws_stack_mesh_3d(&d, &m);
	return m;
}

static void smt_backwards(void* user, unsigned int count, qaws_batch_task_fn task, void* ctx)
{
	unsigned int i;
	(void)user;
	for (i = count; i-- > 0;) task(ctx, i, i + 1);
}

static void smt_case(char const* name, qaws_stack_2d const* st, unsigned int sections, double want, double rel)
{
	qaws_mesh_3d* m = smt_mesh(st, sections, NULL);
	unsigned int open = 0;
	double vol = 0;
	char msg[300];
	if (m) smt_check(m, &open, &vol);
	sprintf(msg, "%s: %u triangles, %u open edges, volume %.6f (%.6f)", name, m ? qaws_mesh_3d_get_triangle_count(m) : 0, open, vol, want);
	TEST_ASSERT(m && qaws_mesh_3d_get_triangle_count(m) > 0 && open == 0 && fabs(vol - want) < rel * want, msg);
	qaws_mesh_3d_destroy(m);
}

int test_92_stack_mesh_main(void)
{
	char msg[300];
	g_pass = 0;
	g_fail = 0;
	printf("Test 92: the surface of a stack as a 3D mesh\n");
	{
		smt_stack s;
		memset(&s, 0, sizeof(s));
		smt_level(&s, 0, smt_circle(0, 0, 1), NULL);
		smt_level(&s, 2, smt_circle(0, 0, 1), NULL);
		smt_done(&s, QAWS_STACK_LINEAR, QAWS_STACK_BLEND_AUTO);
		smt_case("cylinder r 1, height 2, capped", &s.st, 0, 2 * SMT_PI, 1e-3);
		smt_free(&s);
	}
	{
		smt_stack s;
		memset(&s, 0, sizeof(s));
		smt_level(&s, 0, smt_circle(0, 0, 1), NULL);
		smt_level(&s, 1, smt_circle(0, 0, 0), NULL);
		smt_done(&s, QAWS_STACK_LINEAR, QAWS_STACK_BLEND_AUTO);
		smt_case("cone r 1, height 1 (a circle and a point)", &s.st, 0, SMT_PI / 3, 2e-3);
		smt_free(&s);
	}
	{
		smt_stack s;
		memset(&s, 0, sizeof(s));
		smt_level(&s, 0, smt_circle(0, 0, 2), smt_hole(1));
		smt_level(&s, 1, smt_circle(0, 0, 2), smt_hole(1));
		smt_done(&s, QAWS_STACK_LINEAR, QAWS_STACK_BLEND_AUTO);
		smt_case("tube r 2 with a hole r 1, height 1", &s.st, 0, 3 * SMT_PI, 1e-3);
		smt_free(&s);
	}
	{
		smt_stack s;
		unsigned int i;
		memset(&s, 0, sizeof(s));
		for (i = 0; i <= 16; i++)
		{
			double th = SMT_PI * (16 - i) / 16;
			smt_level(&s, cos(th), smt_circle(0, 0, sin(th)), NULL);
		}
		smt_done(&s, QAWS_STACK_CUBIC, QAWS_STACK_BLEND_AUTO);
		smt_case("sphere from 17 circles, cubic, 4 sections between levels", &s.st, 4, 4 * SMT_PI / 3, 6e-3);
		/* the same on an executor running backwards */
		{
			qaws_batch_executor ex;
			qaws_mesh_3d *a = smt_mesh(&s.st, 2, NULL), *b;
			unsigned int same = 1, k;
			memset(&ex, 0, sizeof(ex));
			ex.parallel_for = smt_backwards;
			b = smt_mesh(&s.st, 2, &ex);
			same = a && b && qaws_mesh_3d_get_vertex_count(a) == qaws_mesh_3d_get_vertex_count(b) && qaws_mesh_3d_get_triangle_count(a) == qaws_mesh_3d_get_triangle_count(b);
			for (k = 0; same && k < qaws_mesh_3d_get_vertex_count(a); k++)
			{
				qaws_vec3 p = qaws_mesh_3d_get_vertex(a, k), q = qaws_mesh_3d_get_vertex(b, k);
				same = p.x == q.x && p.y == q.y && p.z == q.z;
			}
			sprintf(msg, "sphere: the same mesh on an executor (%d)", same);
			TEST_ASSERT(same, msg);
			qaws_mesh_3d_destroy(a); qaws_mesh_3d_destroy(b);
		}
		smt_free(&s);
	}
	{
		/* an ellipse at z 0 splitting into two circles at z 1 */
		smt_stack s;
		qaws_curve* e = NULL;
		qaws_mesh_3d* m;
		unsigned int open = 0;
		double vol = 0;
		memset(&s, 0, sizeof(s));
		qaws_curve_create_ellipse_2d(0, 0, 3, (qaws_scalar)1.2, 0, &e);
		smt_level(&s, 0, e, NULL);
		smt_level(&s, 1, smt_circle(-2, 0, 1), smt_circle(2, 0, 1));
		smt_done(&s, QAWS_STACK_LINEAR, QAWS_STACK_BLEND_DISTANCE);
		m = smt_mesh(&s.st, 8, NULL);
		if (m) smt_check(m, &open, &vol);
		sprintf(msg, "ellipse splitting into two circles, 8 sections: %u triangles, %u open edges, volume %.4f (between %.4f and %.4f)",
			m ? qaws_mesh_3d_get_triangle_count(m) : 0, open, vol, 2 * SMT_PI, 3.6 * SMT_PI);
		TEST_ASSERT(m && open == 0 && vol > 2 * SMT_PI && vol < 3.6 * SMT_PI, msg);
		qaws_mesh_3d_destroy(m);
		smt_free(&s);
	}
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
