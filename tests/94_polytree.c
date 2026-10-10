/*
 * Test 94: nesting parity with Clipper2 (its TestPolytreeHoles.cpp in C):
 * PolytreeHoleOwner.txt and PolytreeHoleOwner2.txt, and the polytree
 * regressions #618, #942, #957, #973, on both engines (qaws_clip64 exact,
 * qaws_clip on polylines). Clipper2's tests read children by position;
 * here the same shapes are checked without depending on the order.
 */

#include "test_common.h"
#include <math.h>
#include <string.h>

/* a result as paths of points with parents and hole flags */
typedef struct p94_tree
{
	unsigned int n;
	double** xy;
	unsigned int* np;
	unsigned int* parent;
	int* hole;
} p94_tree;

static void p94_free(p94_tree* t)
{
	unsigned int i;
	for (i = 0; i < t->n; i++) free(t->xy[i]);
	free(t->xy); free(t->np); free(t->parent); free(t->hole);
	memset(t, 0, sizeof(*t));
}

static void p94_alloc(p94_tree* t, unsigned int n)
{
	t->n = n;
	t->xy = (double**)calloc(n + 1, sizeof(double*));
	t->np = (unsigned int*)calloc(n + 1, sizeof(unsigned int));
	t->parent = (unsigned int*)calloc(n + 1, sizeof(unsigned int));
	t->hole = (int*)calloc(n + 1, sizeof(int));
}

typedef struct p94_input
{
	int64_t* sub[8192]; unsigned int nsub[8192]; unsigned int subjects;
	int64_t* clp[64]; unsigned int nclp[64]; unsigned int clips;
	qaws_clip_type ct;
	qaws_fill_rule fr;
} p94_input;

static void p94_add(int64_t** arr, unsigned int* counts, unsigned int* n, long long const* xy, unsigned int pairs)
{
	unsigned int i;
	arr[*n] = (int64_t*)malloc(sizeof(int64_t) * 2 * (pairs + 1));
	for (i = 0; i < 2 * pairs; i++) arr[*n][i] = (int64_t)xy[i];
	counts[*n] = pairs;
	(*n)++;
}

static void p94_input_free(p94_input* in)
{
	unsigned int i;
	for (i = 0; i < in->subjects; i++) free(in->sub[i]);
	for (i = 0; i < in->clips; i++) free(in->clp[i]);
	in->subjects = in->clips = 0;
}

/* exact engine */
static int p94_run64(p94_input const* in, p94_tree* t)
{
	qaws_path64* s = (qaws_path64*)malloc(sizeof(qaws_path64) * (in->subjects + 1));
	qaws_path64 c[64];
	qaws_clip64_desc d;
	qaws_clip64_result* r = NULL;
	unsigned int i;
	int ok;
	for (i = 0; i < in->subjects; i++) { s[i].points = in->sub[i]; s[i].point_count = in->nsub[i]; }
	for (i = 0; i < in->clips; i++) { c[i].points = in->clp[i]; c[i].point_count = in->nclp[i]; }
	memset(&d, 0, sizeof(d));
	d.subjects = s; d.subject_count = in->subjects; d.clips = c; d.clip_count = in->clips;
	d.clip_type = in->ct; d.fill_rule = in->fr;
	ok = qaws_clip64_execute(&d, &r) == QAWS_STATUS_OK && r;
	free(s);
	if (!ok) return 0;
	p94_alloc(t, qaws_clip64_result_get_path_count(r));
	for (i = 0; i < t->n; i++)
	{
		int64_t const* p;
		unsigned int n = 0, k;
		qaws_clip64_result_get_path(r, i, &p, &n);
		t->xy[i] = (double*)malloc(sizeof(double) * 2 * (n + 1));
		for (k = 0; k < 2 * n; k++) t->xy[i][k] = (double)p[k];
		t->np[i] = n;
		t->parent[i] = qaws_clip64_result_get_parent(r, i);
		t->hole[i] = qaws_clip64_result_is_hole(r, i);
	}
	qaws_clip64_result_destroy(r);
	return 1;
}

/* float engine on polylines */
static int p94_runf(p94_input const* in, p94_tree* t)
{
	qaws_curve** cs = (qaws_curve**)malloc(sizeof(qaws_curve*) * (in->subjects + in->clips + 1));
	qaws_curve const** cv = (qaws_curve const**)malloc(sizeof(qaws_curve*) * (in->subjects + in->clips + 1));
	qaws_path_2d* s = (qaws_path_2d*)malloc(sizeof(qaws_path_2d) * (in->subjects + 1));
	qaws_path_2d c[64];
	qaws_clip_desc d;
	qaws_clip_result* r = NULL;
	unsigned int i, k, nc = 0;
	for (i = 0; i < in->subjects + in->clips; i++)
	{
		int64_t const* p = i < in->subjects ? in->sub[i] : in->clp[i - in->subjects];
		unsigned int n = i < in->subjects ? in->nsub[i] : in->nclp[i - in->subjects];
		qaws_scalar* xy = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * (n + 1));
		qaws_path_2d* pp = i < in->subjects ? &s[i] : &c[i - in->subjects];
		for (k = 0; k < 2 * n; k++) xy[k] = (qaws_scalar)p[k];
		cs[nc] = NULL;
		qaws_curve_create_polyline_2d(xy, n, 1, &cs[nc]);
		free(xy);
		cv[nc] = cs[nc];
		pp->curves = &cv[nc]; pp->curve_count = cs[nc] ? 1 : 0; pp->closed = 1;
		nc++;
	}
	memset(&d, 0, sizeof(d));
	d.subjects = s; d.subject_count = in->subjects; d.clips = c; d.clip_count = in->clips;
	d.clip_type = in->ct; d.fill_rule = in->fr;
	if (qaws_clip_execute(&d, &r) != QAWS_STATUS_OK || !r) { for (i = 0; i < nc; i++) qaws_curve_destroy(cs[i]); free(cs); free((void*)cv); free(s); return 0; }
	p94_alloc(t, qaws_clip_result_get_path_count(r));
	for (i = 0; i < t->n; i++)
	{
		qaws_clip_vertex const* vt = NULL;
		unsigned int nv = 0;
		qaws_clip_result_get_vertices(r, 0, i, &vt, &nv);
		t->xy[i] = (double*)malloc(sizeof(double) * 2 * (nv + 1));
		for (k = 0; k < nv; k++) { t->xy[i][2 * k] = (double)vt[k].position.x; t->xy[i][2 * k + 1] = (double)vt[k].position.y; }
		t->np[i] = nv;
		t->parent[i] = qaws_clip_result_get_parent(r, i);
		t->hole[i] = qaws_clip_result_is_hole(r, i);
	}
	qaws_clip_result_destroy(r);
	for (i = 0; i < nc; i++) qaws_curve_destroy(cs[i]);
	free(cs); free((void*)cv); free(s);
	return 1;
}

/* point in polygon: 1 inside, 0 outside, -1 on (Clipper2's three results) */
static int p94_pip(double const* p, unsigned int n, double x, double y)
{
	unsigned int i;
	int in = 0;
	for (i = 0; i < n; i++)
	{
		unsigned int j = (i + n - 1) % n;
		double ax = p[2 * j], ay = p[2 * j + 1], bx = p[2 * i], by = p[2 * i + 1];
		double cr = (bx - ax) * (y - ay) - (by - ay) * (x - ax);
		if (cr == 0 && x >= fmin(ax, bx) && x <= fmax(ax, bx) && y >= fmin(ay, by) && y <= fmax(ay, by)) return -1;
		if ((ay > y) != (by > y) && x < ax + (y - ay) * (bx - ax) / (by - ay)) in = !in;
	}
	return in;
}

/* Clipper2's CheckPolytreeFullyContainsChildren: no child with more than
   one vertex (net) outside its parent */
static int p94_contains(p94_tree const* t)
{
	unsigned int i, k;
	for (i = 0; i < t->n; i++)
	{
		int out = 0;
		unsigned int p = t->parent[i];
		if (p == QAWS_CLIP_NONE_INDEX) continue;
		for (k = 0; k < t->np[i]; k++)
		{
			int r = p94_pip(t->xy[p], t->np[p], t->xy[i][2 * k], t->xy[i][2 * k + 1]);
			if (r == 1) out--;
			else if (r == 0) out++;
			if (out > 1) return 0;
			if (out < -1) break;
		}
	}
	return 1;
}

static unsigned int p94_children(p94_tree const* t, unsigned int p)
{
	unsigned int i, c = 0;
	for (i = 0; i < t->n; i++) if (t->parent[i] == p) c++;
	return c;
}

/* Clipper2's PolytreeContainsPoint: outers +1, holes -1 */
static int p94_tree_contains(p94_tree const* t, double x, double y)
{
	unsigned int i;
	int c = 0;
	for (i = 0; i < t->n; i++)
		if (p94_pip(t->xy[i], t->np[i], x, y) != 0) c += t->hole[i] ? -1 : 1;
	return c != 0;
}

static double p94_area(double const* p, unsigned int n)
{
	unsigned int i;
	double a = 0;
	for (i = 0; i < n; i++)
	{
		unsigned int j = (i + n - 1) % n;
		a += (p[2 * j + 1] + p[2 * i + 1]) * (p[2 * j] - p[2 * i]);
	}
	return 0.5 * a;
}

/* ------------------------------------------------------------- loading */

static int p94_load(char const* file, unsigned int num, p94_input* in)
{
	FILE* f = fopen(file, "r");
	char* line = (char*)malloc(1 << 20);
	int section = 0, found = 0;
	unsigned int seen = 0;
	if (!f || !line) { if (f) fclose(f); free(line); return 0; }
	in->ct = QAWS_CLIP_UNION; in->fr = QAWS_FILL_EVEN_ODD;
	while (fgets(line, 1 << 20, f))
	{
		/* tests by their order in the file, as Clipper2's loader counts them */
		if (!strncmp(line, "CAPTION:", 8)) { if (found) break; found = ++seen == num; section = 0; continue; }
		if (!found) continue;
		if (!strncmp(line, "CLIPTYPE:", 9))
		{
			if (strstr(line, "INTERSECTION")) in->ct = QAWS_CLIP_INTERSECTION;
			else if (strstr(line, "UNION")) in->ct = QAWS_CLIP_UNION;
			else if (strstr(line, "DIFFERENCE")) in->ct = QAWS_CLIP_DIFFERENCE;
			else if (strstr(line, "XOR")) in->ct = QAWS_CLIP_XOR;
			continue;
		}
		if (!strncmp(line, "FILLRULE:", 9))
		{
			if (strstr(line, "EVENODD")) in->fr = QAWS_FILL_EVEN_ODD;
			else if (strstr(line, "NONZERO")) in->fr = QAWS_FILL_NON_ZERO;
			else if (strstr(line, "POSITIVE")) in->fr = QAWS_FILL_POSITIVE;
			else if (strstr(line, "NEGATIVE")) in->fr = QAWS_FILL_NEGATIVE;
			continue;
		}
		if (!strncmp(line, "SUBJECTS_OPEN", 13)) { section = 3; continue; }
		if (!strncmp(line, "SUBJECTS", 8)) { section = 1; continue; }
		if (!strncmp(line, "CLIPS", 5)) { section = 2; continue; }
		if (section == 1 || section == 2)
		{
			long long* tmp = (long long*)malloc(sizeof(long long) * (1 << 17));
			unsigned int k = 0;
			char* s = line;
			long long x, y;
			int used;
			while (tmp && sscanf(s, "%lld,%lld%n", &x, &y, &used) == 2 && k + 2 < (1u << 17))
			{
				tmp[k++] = x; tmp[k++] = y;
				s += used;
				while (*s == ' ' || *s == ',') s++;
			}
			if (k)
			{
				if (section == 1 && in->subjects < 8192) p94_add(in->sub, in->nsub, &in->subjects, tmp, k / 2);
				if (section == 2 && in->clips < 64) p94_add(in->clp, in->nclp, &in->clips, tmp, k / 2);
			}
			free(tmp);
		}
	}
	fclose(f);
	free(line);
	return found;
}

/* --------------------------------------------------------------- tests */

static char const* const p94_engine[2] = { "exact", "float" };

static int p94_run(p94_input const* in, unsigned int e, p94_tree* t)
{
	return e ? p94_runf(in, t) : p94_run64(in, t);
}

static void test_holes1(void)
{
	char file[512], msg[200];
	unsigned int e;
	p94_input in;
	memset(&in, 0, sizeof(in));
	sprintf(file, "%s/clipper2/PolytreeHoleOwner.txt", QAWS_TEST_DATA_DIR);
	if (!p94_load(file, 1, &in)) { TEST_ASSERT(0, "PolytreeHoleOwner.txt loaded"); return; }
	for (e = 0; e < 2; e++)
	{
		p94_tree t;
		memset(&t, 0, sizeof(t));
		sprintf(msg, "TestPolytreeHoles1 (%s): every child inside its parent", p94_engine[e]);
		TEST_ASSERT(p94_run(&in, e, &t) && t.n > 0 && p94_contains(&t), msg);
		p94_free(&t);
	}
	p94_input_free(&in);
}

static void test_holes2(void)
{
	static double const outside[4][2] = { { 21887, 10420 }, { 21726, 10825 }, { 21662, 10845 }, { 21617, 10890 } };
	static double const inside[4][2] = { { 21887, 10430 }, { 21843, 10520 }, { 21810, 10686 }, { 21900, 10461 } };
	char file[512], msg[300];
	unsigned int e, i;
	double sub_area = 0;
	p94_input in;
	memset(&in, 0, sizeof(in));
	sprintf(file, "%s/clipper2/PolytreeHoleOwner2.txt", QAWS_TEST_DATA_DIR);
	if (!p94_load(file, 1, &in)) { TEST_ASSERT(0, "PolytreeHoleOwner2.txt loaded"); return; }
	in.fr = QAWS_FILL_NEGATIVE;
	for (i = 0; i < in.subjects; i++)
	{
		double* d = (double*)malloc(sizeof(double) * 2 * in.nsub[i]);
		unsigned int k;
		for (k = 0; k < 2 * in.nsub[i]; k++) d[k] = (double)in.sub[i][k];
		sub_area -= p94_area(d, in.nsub[i]);
		free(d);
	}
	for (e = 0; e < 2; e++)
	{
		p94_tree t;
		double area = 0;
		int out_ok = 1, in_ok = 1;
		memset(&t, 0, sizeof(t));
		if (!p94_run(&in, e, &t)) { TEST_ASSERT(0, "PolytreeHoleOwner2 clipped"); continue; }
		for (i = 0; i < t.n; i++) area += p94_area(t.xy[i], t.np[i]);
		/* negative fill: the result runs clockwise; Clipper2 compares magnitudes */
		area = fabs(area);
		for (i = 0; i < 4; i++)
		{
			out_ok = out_ok && !p94_tree_contains(&t, outside[i][0], outside[i][1]);
			in_ok = in_ok && p94_tree_contains(&t, inside[i][0], inside[i][1]);
		}
		sprintf(msg, "TestPolytreeHoles2 (%s): area %.0f between 0.92 and 1 times the subjects' %.0f, children inside parents, points of interest in and out (%d, %d)",
			p94_engine[e], area, sub_area, out_ok, in_ok);
		TEST_ASSERT(area < sub_area && area > 0.92 * sub_area && p94_contains(&t) && out_ok && in_ok, msg);
		p94_free(&t);
	}
	p94_input_free(&in);
}

/* the shape of the tree, order apart */
typedef struct p94_expect
{
	char const* name;
	unsigned int roots;                 /* top-level paths */
	unsigned int root_children[4];      /* children counts of the roots, as a multiset (sorted descending) */
	unsigned int grandchild;            /* some child of a root has this many children (SIZE_MAX: not checked) */
} p94_expect;

static int p94_cmp_desc(void const* a, void const* b)
{
	unsigned int x = *(unsigned int const*)a, y = *(unsigned int const*)b;
	return x > y ? -1 : (x < y ? 1 : 0);
}

static void p94_check(p94_input const* in, p94_expect const* x)
{
	unsigned int e;
	for (e = 0; e < 2; e++)
	{
		p94_tree t;
		unsigned int i, roots = 0, rc[64], grand_ok = x->grandchild == 0xFFFFFFFFu;
		int ok;
		char msg[300];
		memset(&t, 0, sizeof(t));
		if (!p94_run(in, e, &t)) { TEST_ASSERT(0, x->name); continue; }
		for (i = 0; i < t.n; i++)
			if (t.parent[i] == QAWS_CLIP_NONE_INDEX && roots < 64) rc[roots++] = p94_children(&t, i);
		qsort(rc, roots, sizeof(unsigned int), p94_cmp_desc);
		ok = roots == x->roots;
		for (i = 0; ok && i < roots && i < 4; i++) ok = rc[i] == x->root_children[i];
		for (i = 0; !grand_ok && i < t.n; i++)
			if (t.parent[i] != QAWS_CLIP_NONE_INDEX && t.parent[t.parent[i]] == QAWS_CLIP_NONE_INDEX && p94_children(&t, i) == x->grandchild)
				grand_ok = 1;
		sprintf(msg, "%s (%s): %u top-level paths, children %u %u %u, grandchildren as Clipper2's, children inside parents",
			x->name, p94_engine[e], roots, roots > 0 ? rc[0] : 0, roots > 1 ? rc[1] : 0, roots > 2 ? rc[2] : 0);
		TEST_ASSERT(ok && grand_ok && p94_contains(&t), msg);
		p94_free(&t);
	}
}

static void test_holes3_to_10(void)
{
	p94_input in;
	memset(&in, 0, sizeof(in));
	{
		static long long const s[] = { 1072,501, 1072,501, 1072,539, 1072,539, 1072,539, 870,539, 870,539, 870,539, 870,520, 894,520, 898,524, 911,524, 915,520, 915,520, 936,520,
			940,524, 953,524, 957,520, 957,520, 978,520, 983,524, 995,524, 1000,520, 1021,520, 1025,524, 1038,524, 1042,520, 1038,516, 1025,516, 1021,520, 1000,520, 995,516,
			983,516, 978,520, 957,520, 953,516, 940,516, 936,520, 915,520, 911,516, 898,516, 894,520, 870,520, 870,516, 870,501, 870,501, 870,501, 1072,501 };
		static long long const c[] = { 870,501, 971,501, 971,539, 870,539 };
		static p94_expect const x = { "TestPolytreeHoles3", 1, { 2 }, 0xFFFFFFFFu };
		p94_add(in.sub, in.nsub, &in.subjects, s, sizeof(s) / sizeof(s[0]) / 2);
		p94_add(in.clp, in.nclp, &in.clips, c, 4);
		in.ct = QAWS_CLIP_INTERSECTION; in.fr = QAWS_FILL_NON_ZERO;
		p94_check(&in, &x);
		p94_input_free(&in);
	}
	{
		static long long const a[] = { 50,500, 50,300, 100,300, 100,350, 150,350, 150,250, 200,250, 200,450, 350,450, 350,200, 400,200, 400,225, 450,225,
			450,175, 400,175, 400,200, 350,200, 350,175, 200,175, 200,250, 150,250, 150,200, 100,200, 100,300, 50,300, 50,125, 500,125, 500,500 };
		static long long const b[] = { 250,425, 250,375, 300,375, 300,425 };
		static p94_expect const x = { "TestPolytreeHoles4 (#618)", 1, { 3 }, 1 };
		p94_add(in.sub, in.nsub, &in.subjects, a, sizeof(a) / sizeof(a[0]) / 2);
		p94_add(in.sub, in.nsub, &in.subjects, b, 4);
		in.ct = QAWS_CLIP_UNION; in.fr = QAWS_FILL_NON_ZERO;
		p94_check(&in, &x);
		p94_input_free(&in);
	}
	{
		static long long const s[] = { 0,30, 400,30, 400,100, 0,100 };
		static long long const c1[] = { 20,30, 30,30, 30,150, 20,150 };
		static long long const c2[] = { 200,0, 300,0, 300,30, 280,30, 280,20, 220,20, 220,30, 200,30 };
		static long long const c3[] = { 200,50, 300,50, 300,80, 200,80 };
		static p94_expect const x = { "TestPolytreeHoles5", 3, { 2, 0, 0 }, 0xFFFFFFFFu };
		p94_add(in.sub, in.nsub, &in.subjects, s, 4);
		p94_add(in.clp, in.nclp, &in.clips, c1, 4);
		p94_add(in.clp, in.nclp, &in.clips, c2, 8);
		p94_add(in.clp, in.nclp, &in.clips, c3, 4);
		in.ct = QAWS_CLIP_XOR; in.fr = QAWS_FILL_NON_ZERO;
		p94_check(&in, &x);
		p94_input_free(&in);
	}
	{
		static long long const s1[] = { 150,50, 200,50, 200,100, 150,100 };
		static long long const s2[] = { 125,100, 150,100, 150,150, 125,150 };
		static long long const s3[] = { 225,50, 300,50, 300,80, 225,80 };
		static long long const s4[] = { 225,100, 300,100, 300,150, 275,150, 275,175, 260,175, 260,250, 235,250, 235,300, 275,300, 275,275, 300,275, 300,350, 225,350 };
		static long long const s5[] = { 300,150, 350,150, 350,175, 300,175 };
		static long long const c1[] = { 0,0, 400,0, 400,50, 0,50 };
		static long long const c2[] = { 0,100, 400,100, 400,150, 0,150 };
		static long long const c3[] = { 260,175, 325,175, 325,275, 260,275 };
		static p94_expect const x = { "TestPolytreeHoles6 (#618)", 3, { 1, 0, 0 }, 0xFFFFFFFFu };
		p94_add(in.sub, in.nsub, &in.subjects, s1, 4); p94_add(in.sub, in.nsub, &in.subjects, s2, 4);
		p94_add(in.sub, in.nsub, &in.subjects, s3, 4); p94_add(in.sub, in.nsub, &in.subjects, s4, 14);
		p94_add(in.sub, in.nsub, &in.subjects, s5, 4);
		p94_add(in.clp, in.nclp, &in.clips, c1, 4); p94_add(in.clp, in.nclp, &in.clips, c2, 4); p94_add(in.clp, in.nclp, &in.clips, c3, 4);
		in.ct = QAWS_CLIP_XOR; in.fr = QAWS_FILL_NON_ZERO;
		p94_check(&in, &x);
		p94_input_free(&in);
	}
	{
		static long long const a[] = { 0,0, 100000,0, 100000,100000, 200000,100000, 200000,0, 300000,0, 300000,200000, 0,200000 };
		static long long const b[] = { 0,0, 0,-100000, 250000,-100000, 250000,0 };
		static p94_expect const x = { "TestPolytreeHoles7 (#618)", 1, { 1 }, 0xFFFFFFFFu };
		p94_add(in.sub, in.nsub, &in.subjects, a, 8);
		p94_add(in.sub, in.nsub, &in.subjects, b, 4);
		in.ct = QAWS_CLIP_UNION; in.fr = QAWS_FILL_NON_ZERO;
		p94_check(&in, &x);
		p94_input_free(&in);
	}
	{
		static long long const a[] = { 1588700,-8717600, 1616200,-8474800, 1588700,-8474800 };
		static long long const b[] = { 13583800,-15601600, 13582800,-15508500, 13555300,-15508500, 13555500,-15182200, 13010900,-15185400 };
		static long long const c[] = { 956700,-3092300, 1152600,3147400, 25600,3151700 };
		static long long const d[] = { 22575900,-16604000, 31286800,-12171900, 31110200,4882800, 30996200,4826300, 30414400,5447400, 30260000,5391500,
			29662200,5805400, 28844500,5337900, 28435000,5789300, 27721400,5026400, 22876300,5034300, 21977700,4414900, 21148000,4654700, 20917600,4653400,
			19334300,12411000, -2591700,12177200, 53200,3151100, -2564300,12149800, 7819400,4692400, 10116000,5228600, 6975500,3120100, 7379700,3124700,
			11037900,596200, 12257000,2587800, 12257000,596200, 15227300,2352700, 18444400,1112100, 19961100,5549400, 20173200,5078600, 20330000,5079300,
			20970200,4544300, 20989600,4563700, 19465500,1112100, 21611600,4182100, 22925100,1112200, 22952700,1637200, 23059000,1112200, 24908100,4181200,
			27070100,3800600, 27238000,3800700, 28582200,520300, 29367800,1050100, 29291400,179400, 29133700,360700, 29056700,312600, 29121900,332500,
			29269900,162300, 28941400,213100, 27491300,-3041500, 27588700,-2997800, 22104900,-16142800, 13010900,-15603000, 13555500,-15182200,
			13555300,-15508500, 13582800,-15508500, 13583100,-15154700, 1588700,-8822800, 1588700,-8379900, 1588700,-8474800, 1616200,-8474800,
			1003900,-630100, 1253300,-12284500, 12983400,-16239900 };
		static long long const e2[] = { 198200,12149800, 1010600,12149800, 1011500,11859600 };
		static long long const f[] = { 21996700,-7432000, 22096700,-7432000, 22096700,-7332000 };
		static p94_expect const x = { "TestPolytreeHoles8 (#942)", 1, { 2 }, 1 };
		p94_add(in.sub, in.nsub, &in.subjects, a, 3); p94_add(in.sub, in.nsub, &in.subjects, b, 5);
		p94_add(in.sub, in.nsub, &in.subjects, c, 3); p94_add(in.sub, in.nsub, &in.subjects, d, sizeof(d) / sizeof(d[0]) / 2);
		p94_add(in.sub, in.nsub, &in.subjects, e2, 3); p94_add(in.sub, in.nsub, &in.subjects, f, 3);
		in.ct = QAWS_CLIP_UNION; in.fr = QAWS_FILL_NON_ZERO;
		p94_check(&in, &x);
		p94_input_free(&in);
	}
	{
		static long long const r[17][12] = {
			{ 77910,46865, 78720,46865, 78720,48000, 77910,48000, 77910,46865 },
			{ 82780,53015, 93600,53015, 93600,54335, 82780,54335, 82780,53015 },
			{ 82780,48975, 84080,48975, 84080,53015, 82780,53015, 82780,48975 },
			{ 77910,48000, 84080,48000, 84080,48975, 77910,48975, 77910,48000 },
			{ 89880,40615, 90700,40615, 90700,46865, 89880,46865, 89880,40615 },
			{ 92700,54335, 93600,54335, 93600,61420, 92700,61420, 92700,54335 },
			{ 78950,47425, 84080,47425, 84080,47770, 78950,47770, 78950,47425 },
			{ 82780,61420, 93600,61420, 93600,62435, 82780,62435, 82780,61420 },
			{ 0 },
			{ 76195,39880, 89880,39880, 89880,41045, 76195,41045, 76195,39880 },
			{ 85490,56145, 90520,56145, 90520,59235, 85490,59235, 85490,56145 },
			{ 89880,39880, 101680,39880, 101680,40615, 89880,40615, 89880,39880 },
			{ 89880,46865, 100680,46865, 100680,47770, 89880,47770, 89880,46865 },
			{ 82780,54335, 83280,54335, 83280,61420, 82780,61420, 82780,54335 },
			{ 76195,41045, 76855,41045, 76855,62665, 76195,62665, 76195,41045 },
			{ 76195,62665, 100675,62665, 100675,63085, 76195,63085, 76195,62665 },
			{ 82780,41045, 84080,41045, 84080,47425, 82780,47425, 82780,41045 } };
		static long long const big[] = { 101680,63085, 100675,63085, 100675,47770, 100680,47770, 100680,40615, 101680,40615, 101680,63085 };
		static p94_expect const x = { "TestPolytreeHoles9 (#957)", 1, { 2 }, 1 };
		unsigned int i;
		for (i = 0; i < 17; i++)
			if (i == 8) p94_add(in.sub, in.nsub, &in.subjects, big, 7);
			else p94_add(in.sub, in.nsub, &in.subjects, r[i], 5);
		in.ct = QAWS_CLIP_UNION; in.fr = QAWS_FILL_NON_ZERO;
		p94_check(&in, &x);
		p94_input_free(&in);
	}
	{
		static long long const a[] = { 0,0, 79530,0, 79530,940, 0,940, 0,0 };
		static long long const b[] = { 0,33360, 79530,33360, 79530,34300, 0,34300, 0,33360 };
		static long long const c[] = { 78470,940, 79530,940, 79530,33360, 78470,33360, 78470,940 };
		static long long const d[] = { 0,940, 940,940, 940,33360, 0,33360, 0,940 };
		static long long const e2[] = { 29290,940, 30350,940, 30350,33360, 29290,33360, 29290,940 };
		static p94_expect const x = { "TestPolytreeHoles10 (#973)", 1, { 2 }, 0xFFFFFFFFu };
		p94_add(in.sub, in.nsub, &in.subjects, a, 5); p94_add(in.sub, in.nsub, &in.subjects, b, 5);
		p94_add(in.sub, in.nsub, &in.subjects, c, 5); p94_add(in.sub, in.nsub, &in.subjects, d, 5);
		p94_add(in.sub, in.nsub, &in.subjects, e2, 5);
		in.ct = QAWS_CLIP_UNION; in.fr = QAWS_FILL_NON_ZERO;
		p94_check(&in, &x);
		p94_input_free(&in);
	}
}

int test_94_polytree_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 94: nesting parity with Clipper2 (PolyTree)\n");
	test_holes1();
	test_holes2();
	test_holes3_to_10();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
