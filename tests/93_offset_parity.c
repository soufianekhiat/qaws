/*
 * Test 93: offset parity with Clipper2 (qaws_offset64, and the corpus on
 * qaws_offset): Clipper2's TestOffsets.cpp and TestOffsetOrientation.cpp
 * re-expressed in C, on its Offsets.txt corpus and the paths of its
 * regression issues (#405, #424, #448, #456, #482, #593, #715, #724, #733,
 * #873, #965), with Clipper2's own expectations.
 */

#include "test_common.h"
#include "data/clipper2/offsets_paths.h"
#include <math.h>
#include <string.h>

#define O93_PI 3.141592653589793238

typedef struct o93_paths
{
	int64_t* xy[64];
	unsigned int n[64];
	unsigned int count;
} o93_paths;

static void o93_add(o93_paths* p, long long const* xy, unsigned int pairs)
{
	unsigned int i;
	p->xy[p->count] = (int64_t*)malloc(sizeof(int64_t) * 2 * (pairs + 1));
	for (i = 0; i < 2 * pairs; i++) p->xy[p->count][i] = (int64_t)xy[i];
	p->n[p->count] = pairs;
	p->count++;
}

static void o93_reverse(o93_paths* p, unsigned int k)
{
	unsigned int i, n = p->n[k];
	for (i = 0; i < n / 2; i++)
	{
		int64_t x = p->xy[k][2 * i], y = p->xy[k][2 * i + 1];
		p->xy[k][2 * i] = p->xy[k][2 * (n - 1 - i)]; p->xy[k][2 * i + 1] = p->xy[k][2 * (n - 1 - i) + 1];
		p->xy[k][2 * (n - 1 - i)] = x; p->xy[k][2 * (n - 1 - i) + 1] = y;
	}
}

static void o93_free(o93_paths* p)
{
	unsigned int i;
	for (i = 0; i < p->count; i++) free(p->xy[i]);
	p->count = 0;
}

static qaws_clip64_result* o93_inflate(o93_paths const* p, double delta, qaws_join_type jt, qaws_end_type et, double miter, double arc_tol, unsigned int flags)
{
	qaws_path64 paths[64];
	qaws_offset64_group g;
	qaws_offset64_desc d;
	qaws_clip64_result* r = NULL;
	unsigned int i;
	for (i = 0; i < p->count; i++) { paths[i].points = p->xy[i]; paths[i].point_count = p->n[i]; }
	g.paths = paths; g.path_count = p->count; g.join_type = jt; g.end_type = et;
	memset(&d, 0, sizeof(d));
	d.groups = &g; d.group_count = 1; d.delta = delta; d.miter_limit = miter; d.arc_tolerance = arc_tol; d.flags = flags;
	qaws_offset64_execute(&d, &r);
	return r;
}

static unsigned int o93_count(qaws_clip64_result const* r) { return r ? qaws_clip64_result_get_path_count(r) : 0; }

static unsigned int o93_size(qaws_clip64_result const* r, unsigned int i)
{
	int64_t const* p;
	unsigned int n = 0;
	if (!r || i >= o93_count(r)) return 0;
	qaws_clip64_result_get_path(r, i, &p, &n);
	return n;
}

/* Clipper2's Area (half twice the signed area) */
static double o93_area(qaws_clip64_result const* r, unsigned int i)
{
	int64_t const* p;
	unsigned int n = 0;
	qaws_clip64_result_get_path(r, i, &p, &n);
	return 0.5 * qaws_path64_area2(p, n);
}

static double o93_total_area(qaws_clip64_result const* r)
{
	unsigned int i;
	double a = 0;
	for (i = 0; i < o93_count(r); i++) a += o93_area(r, i);
	return a;
}

/* ---------------------------------------------------------------- corpus */

/* subjects of test `num` of a Clipper2 test file */
static int o93_load(char const* file, unsigned int num, o93_paths* out)
{
	FILE* f = fopen(file, "r");
	char line[1 << 16];
	int in = 0, found = 0;
	if (!f) return 0;
	while (fgets(line, sizeof(line), f))
	{
		unsigned int c;
		if (sscanf(line, "CAPTION: %u", &c) == 1) { if (found) break; found = c == num; in = 0; continue; }
		if (!found) continue;
		if (!strncmp(line, "SUBJECTS", 8)) { in = 1; continue; }
		if (!strncmp(line, "CLIPS", 5) || !strncmp(line, "SUBJECTS_OPEN", 13)) { in = 0; continue; }
		if (in)
		{
			long long tmp[1 << 14];
			unsigned int k = 0;
			char* s = line;
			long long x, y;
			int used;
			while (sscanf(s, "%lld,%lld%n", &x, &y, &used) == 2 && k + 2 < (1u << 14))
			{
				tmp[k++] = x; tmp[k++] = y;
				s += used;
				while (*s == ' ' || *s == ',') s++;
			}
			if (k && out->count < 64) o93_add(out, tmp, k / 2);
		}
	}
	fclose(f);
	return found;
}

static void test_corpus(void)
{
	char file[512], msg[300];
	unsigned int t;
	sprintf(file, "%s/clipper2/Offsets.txt", QAWS_TEST_DATA_DIR);
	for (t = 1; t <= 2; t++)
	{
		o93_paths p;
		qaws_clip64_result* r;
		unsigned int i, pos = 0, neg = 0;
		double total;
		memset(&p, 0, sizeof(p));
		if (!o93_load(file, t, &p)) { TEST_ASSERT(0, "Offsets.txt loaded"); continue; }
		r = o93_inflate(&p, 1, QAWS_JOIN_ROUND, QAWS_END_POLYGON, 0, 0, 0);
		total = r ? o93_total_area(r) : 0;
		for (i = 0; i < o93_count(r); i++) { if (o93_area(r, i) > 0) pos++; else neg++; }
		sprintf(msg, "Offsets.txt %u (round, 1): exactly one outer path (%u positive, %u negative, area %s)", t, pos, neg, total > 0 ? "positive" : "negative");
		TEST_ASSERT(r && (total > 0 ? pos == 1 : neg == 1), msg);
		qaws_clip64_result_destroy(r);
		/* the same on curves (qaws_offset): one outer path */
		{
			qaws_curve* cs[64];
			qaws_curve const* cv[64];
			qaws_path_2d pp[64];
			qaws_clip_result* fr = NULL;
			unsigned int outer = 0;
			for (i = 0; i < p.count; i++)
			{
				qaws_scalar* xy = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * p.n[i]);
				unsigned int k;
				for (k = 0; k < 2 * p.n[i]; k++) xy[k] = (qaws_scalar)p.xy[i][k];
				cs[i] = NULL;
				qaws_curve_create_polyline_2d(xy, p.n[i], 1, &cs[i]);
				free(xy);
				cv[i] = cs[i];
				pp[i].curves = &cv[i]; pp[i].curve_count = 1; pp[i].closed = 1;
			}
			qaws_offset_paths(pp, p.count, 1, QAWS_JOIN_ROUND, QAWS_END_POLYGON, &fr);
			for (i = 0; fr && i < qaws_clip_result_get_path_count(fr); i++)
				if (qaws_clip_result_get_parent(fr, i) == QAWS_CLIP_NONE_INDEX) outer++;
			sprintf(msg, "Offsets.txt %u on curves (qaws_offset, round, 1): one outer path (%u)", t, outer);
			TEST_ASSERT(fr && outer == 1, msg);
			qaws_clip_result_destroy(fr);
			for (i = 0; i < p.count; i++) qaws_curve_destroy(cs[i]);
		}
		o93_free(&p);
	}
}

/* ----------------------------------------------------- Clipper2's tests */

static void test_offsets2(void) /* #448, #456 */
{
	static long long const s[] = { 500,500, 1000,500, 1000,1500, 500,1500, 0,1000 };
	double delta = 100, arc_tol = 2.5, min_dist = delta * 2;
	o93_paths p;
	qaws_clip64_result* r;
	char msg[200];
	memset(&p, 0, sizeof(p));
	o93_add(&p, s, 5);
	r = o93_inflate(&p, delta, QAWS_JOIN_ROUND, QAWS_END_POLYGON, 0, arc_tol, 0);
	if (r && o93_count(r))
	{
		int64_t const* q;
		unsigned int n = 0, i, k;
		qaws_clip64_result_get_path(r, 0, &q, &n);
		for (k = 0; k < 5; k++)
			for (i = 0; i < n; i++)
			{
				unsigned int j = (i + n - 1) % n;
				double mx = (double)((q[2 * j] + q[2 * i]) / 2), my = (double)((q[2 * j + 1] + q[2 * i + 1]) / 2);
				double dd = hypot(mx - s[2 * k], my - s[2 * k + 1]);
				if (dd < delta * 2 && dd < min_dist) min_dist = dd;
			}
	}
	sprintf(msg, "TestOffsets2 (#448, #456): nearest midpoint %.2f >= delta - arc tolerance - 1 (%.2f), %u vertices <= 21", min_dist, delta - arc_tol - 1, o93_size(r, 0));
	TEST_ASSERT(r && min_dist + 1 >= delta - arc_tol && o93_size(r, 0) <= 21, msg);
	qaws_clip64_result_destroy(r);
	o93_free(&p);
}

static void test_offsets3(void) /* #424 */
{
	static long long const s[] = {
		1525311078, 1352369439, 1526632284, 1366692987, 1519397110, 1367437476, 1520246456, 1380177674, 1520613458, 1385913385,
		1517383844, 1386238444, 1517771817, 1392099983, 1518233190, 1398758441, 1518421934, 1401883197, 1518694564, 1406612275,
		1520267428, 1430289121, 1520770744, 1438027612, 1521148232, 1443438264, 1521441833, 1448964260, 1521683005, 1452518932,
		1521819320, 1454374912, 1527943004, 1454154711, 1527649403, 1448523858, 1535901696, 1447989084, 1535524209, 1442788147,
		1538953052, 1442463089, 1541553521, 1442242888, 1541459149, 1438855987, 1538764308, 1439076188, 1538575565, 1436832236,
		1538764308, 1436832236, 1536509870, 1405374956, 1550497874, 1404347351, 1550214758, 1402428457, 1543818445, 1402868859,
		1543734559, 1402124370, 1540672717, 1402344571, 1540473487, 1399995761, 1524996506, 1400981422, 1524807762, 1398223667,
		1530092585, 1397898609, 1531675935, 1397783265, 1531392819, 1394920653, 1529809469, 1395025510, 1529348096, 1388880855,
		1531099218, 1388660654, 1530826588, 1385158410, 1532955197, 1384938209, 1532661596, 1379003269, 1532472852, 1376235028,
		1531277476, 1376350372, 1530050642, 1361806623, 1599487345, 1352704983, 1602758902, 1378489467, 1618990858, 1376350372,
		1615058698, 1344085688, 1603230761, 1345700495, 1598648484, 1346329641, 1598931599, 1348667965, 1596698132, 1348993024,
		1595775386, 1342722540 };
	o93_paths p;
	qaws_clip64_result* r;
	char msg[200];
	memset(&p, 0, sizeof(p));
	o93_add(&p, s, (unsigned int)(sizeof(s) / sizeof(s[0]) / 2));
	r = o93_inflate(&p, -209715, QAWS_JOIN_MITER, QAWS_END_POLYGON, 0, 0, 0);
	sprintf(msg, "TestOffsets3 (#424): %u vertices against the subject's %u (at most 1 more)", o93_size(r, 0), p.n[0]);
	TEST_ASSERT(r && o93_count(r) >= 1 && (int)o93_size(r, 0) - (int)p.n[0] <= 1, msg);
	qaws_clip64_result_destroy(r);
	o93_free(&p);
}

static void test_offsets4(void) /* #482 */
{
	static long long const a[] = { 0,0, 20000,200, 40000,0, 40000,50000, 0,50000, 0,0 };
	static long long const b[] = { 0,0, 20000,400, 40000,0, 40000,50000, 0,50000, 0,0 };
	static long long const c[] = { 0,0, 20000,1500, 40000,0, 40000,50000, 0,50000, 0,0 };
	o93_paths p;
	qaws_clip64_result *r1, *r2, *r3, *r4;
	char msg[300];
	memset(&p, 0, sizeof(p)); o93_add(&p, a, 6);
	r1 = o93_inflate(&p, -5000, QAWS_JOIN_SQUARE, QAWS_END_POLYGON, 0, 0, 0); o93_free(&p);
	o93_add(&p, b, 6);
	r2 = o93_inflate(&p, -5000, QAWS_JOIN_SQUARE, QAWS_END_POLYGON, 0, 0, 0);
	r3 = o93_inflate(&p, -5000, QAWS_JOIN_ROUND, QAWS_END_POLYGON, 2, 100, 0); o93_free(&p);
	o93_add(&p, c, 6);
	r4 = o93_inflate(&p, -5000, QAWS_JOIN_ROUND, QAWS_END_POLYGON, 2, 100, 0); o93_free(&p);
	sprintf(msg, "TestOffsets4 (#482): square 5 and 5 vertices (%u, %u), round more than 5 (%u, %u)", o93_size(r1, 0), o93_size(r2, 0), o93_size(r3, 0), o93_size(r4, 0));
	TEST_ASSERT(o93_size(r1, 0) == 5 && o93_size(r2, 0) == 5 && o93_size(r3, 0) > 5 && o93_size(r4, 0) > 5, msg);
	qaws_clip64_result_destroy(r1); qaws_clip64_result_destroy(r2); qaws_clip64_result_destroy(r3); qaws_clip64_result_destroy(r4);
}

static void test_offsets5(void) /* #593: clean-up */
{
	static long long const box[] = { -47877,-47877, 84788,-47877, 84788,81432, -47877,81432 };
	o93_paths p;
	qaws_clip64_result* r;
	char msg[200];
	memset(&p, 0, sizeof(p));
	o93_add(&p, off5_a, (unsigned int)(sizeof(off5_a) / sizeof(off5_a[0]) / 2));
	o93_add(&p, box, 4);
	r = o93_inflate(&p, -10000, QAWS_JOIN_ROUND, QAWS_END_POLYGON, 0, 0, 0);
	sprintf(msg, "TestOffsets5 (#593 clean-up): 2 paths (%u)", o93_count(r));
	TEST_ASSERT(o93_count(r) == 2, msg);
	qaws_clip64_result_destroy(r);
	o93_free(&p);
}

static void test_offsets6(void) /* #593: rounded ends */
{
	static long long const sq[] = { 620,620, -620,620, -620,-620, 620,-620 };
	static long long const arc[] = { 20,-277, 42,-275, 59,-272, 80,-266, 97,-261, 114,-254, 135,-243, 149,-235, 167,-222, 182,-211, 197,-197,
		212,-181, 223,-167, 234,-150, 244,-133, 253,-116, 260,-99, 267,-78, 272,-61, 275,-40, 278,-18, 276,-39,
		272,-61, 267,-79, 260,-99, 253,-116, 245,-133, 235,-150, 223,-167, 212,-181, 197,-197, 182,-211, 168,-222, 152,-233,
		135,-243, 114,-254, 97,-261, 80,-267, 59,-272, 42,-275, 20,-278 };
	o93_paths p;
	qaws_clip64_result* r;
	char msg[200];
	double a = 0;
	unsigned int i;
	memset(&p, 0, sizeof(p));
	o93_add(&p, sq, 4);
	o93_add(&p, arc, (unsigned int)(sizeof(arc) / sizeof(arc[0]) / 2));
	r = o93_inflate(&p, -50, QAWS_JOIN_ROUND, QAWS_END_POLYGON, 0, 0, 0);
	/* the hole (negative area): Clipper2 lists it second */
	for (i = 0; i < o93_count(r); i++) if (o93_area(r, i) < a) a = o93_area(r, i);
	sprintf(msg, "TestOffsets6 (#593 rounded ends): 2 paths (%u), the hole's area below -47500 (%.0f)", o93_count(r), a);
	TEST_ASSERT(o93_count(r) == 2 && a < -47500, msg);
	qaws_clip64_result_destroy(r);
	o93_free(&p);
}

static void test_offsets7(void) /* #593, #715 */
{
	static long long const sq[] = { 0,0, 100,0, 100,100, 0,100 };
	static long long const hole[] = { 40,60, 60,60, 60,40, 40,40 };
	o93_paths p;
	qaws_clip64_result *r1, *r2, *r3, *r4;
	char msg[200];
	memset(&p, 0, sizeof(p));
	o93_add(&p, sq, 4);
	r1 = o93_inflate(&p, -50, QAWS_JOIN_MITER, QAWS_END_POLYGON, 0, 0, 0);
	o93_add(&p, hole, 4);
	r2 = o93_inflate(&p, 10, QAWS_JOIN_MITER, QAWS_END_POLYGON, 0, 0, 0);
	o93_reverse(&p, 0); o93_reverse(&p, 1);
	r3 = o93_inflate(&p, 10, QAWS_JOIN_MITER, QAWS_END_POLYGON, 0, 0, 0);
	p.count = 1;
	r4 = o93_inflate(&p, -50, QAWS_JOIN_MITER, QAWS_END_POLYGON, 0, 0, 0);
	p.count = 2;
	sprintf(msg, "TestOffsets7 (#593, #715): sizes 0, 1, 1, 0 (%u, %u, %u, %u)", o93_count(r1), o93_count(r2), o93_count(r3), o93_count(r4));
	TEST_ASSERT(o93_count(r1) == 0 && o93_count(r2) == 1 && o93_count(r3) == 1 && o93_count(r4) == 0, msg);
	qaws_clip64_result_destroy(r1); qaws_clip64_result_destroy(r2); qaws_clip64_result_destroy(r3); qaws_clip64_result_destroy(r4);
	o93_free(&p);
}

/* Clipper2's GetOffsetQuality: for each solution vertex and edge midpoint,
   its distance to the nearest subject segment; the smallest and largest */
static void o93_quality(int64_t const* sub, unsigned int ns, int64_t const* sol, unsigned int nsol, double* smallest, double* largest)
{
	unsigned int i, j, h;
	*smallest = 1e300;
	*largest = 0;
	for (i = 0; i < nsol; i++)
		for (h = 0; h < 2; h++)
		{
			unsigned int pi = (i + nsol - 1) % nsol;
			double px = h ? 0.5 * ((double)sol[2 * pi] + sol[2 * i]) : (double)sol[2 * i];
			double py = h ? 0.5 * ((double)sol[2 * pi + 1] + sol[2 * i + 1]) : (double)sol[2 * i + 1];
			double best = 1e300;
			for (j = 0; j < ns; j++)
			{
				unsigned int k = (j + 1) % ns;
				double ax = (double)sub[2 * j], ay = (double)sub[2 * j + 1], dx = (double)sub[2 * k] - ax, dy = (double)sub[2 * k + 1] - ay;
				double q = dx * dx + dy * dy > 0 ? ((px - ax) * dx + (py - ay) * dy) / (dx * dx + dy * dy) : 0, d2;
				if (q < 0) q = 0;
				if (q > 1) q = 1;
				d2 = (px - ax - q * dx) * (px - ax - q * dx) + (py - ay - q * dy) * (py - ay - q * dy);
				if (d2 < best) best = d2;
			}
			best = sqrt(best);
			if (best < *smallest) *smallest = best;
			if (best > *largest) *largest = best;
		}
}

static void test_offsets8(void) /* #724 */
{
	o93_paths p;
	qaws_clip64_result* r;
	double offset = -50329979.277800001, arc_tol = 5000, smallest = 0, largest = 0;
	char msg[300];
	memset(&p, 0, sizeof(p));
	o93_add(&p, off8, (unsigned int)(sizeof(off8) / sizeof(off8[0]) / 2));
	r = o93_inflate(&p, offset, QAWS_JOIN_ROUND, QAWS_END_POLYGON, 2, arc_tol, 0);
	if (r && o93_count(r))
	{
		int64_t const* q;
		unsigned int n = 0;
		qaws_clip64_result_get_path(r, 0, &q, &n);
		o93_quality(p.xy[0], p.n[0], q, n, &smallest, &largest);
	}
	offset = fabs(offset);
	sprintf(msg, "TestOffsets8 (#724): distances to the subject from %.1f to %.1f, within the arc tolerance %.0f (+1) of %.1f", smallest, largest, arc_tol, offset);
	TEST_ASSERT(r && o93_count(r) && offset - smallest - 1 <= arc_tol && largest - offset - 1 <= arc_tol, msg);
	qaws_clip64_result_destroy(r);
	o93_free(&p);
}

static void test_offsets9(void) /* #733: orientation */
{
	static long long const s[] = { 100,100, 200,100, 200,400, 100,400 };
	static long long const h[] = { 130,130, 170,130, 170,370, 130,370 };
	o93_paths p;
	qaws_clip64_result* r;
	double sub_area;
	int ok = 1;
	char msg[200];
	memset(&p, 0, sizeof(p));
	o93_add(&p, s, 4);
	sub_area = fabs(0.5 * qaws_path64_area2(p.xy[0], 4));
	r = o93_inflate(&p, 50, QAWS_JOIN_MITER, QAWS_END_POLYGON, 0, 0, 0);
	ok = ok && o93_count(r) == 1 && o93_area(r, 0) > 0;
	qaws_clip64_result_destroy(r);
	o93_reverse(&p, 0);
	r = o93_inflate(&p, 50, QAWS_JOIN_MITER, QAWS_END_POLYGON, 0, 0, 0);
	ok = ok && o93_count(r) == 1 && fabs(o93_area(r, 0)) > sub_area && o93_area(r, 0) < 0;
	qaws_clip64_result_destroy(r);
	r = o93_inflate(&p, 50, QAWS_JOIN_MITER, QAWS_END_POLYGON, 2, 0, QAWS_CLIP_REVERSE_SOLUTION);
	ok = ok && o93_count(r) == 1 && fabs(o93_area(r, 0)) > sub_area && o93_area(r, 0) > 0;
	qaws_clip64_result_destroy(r);
	o93_add(&p, h, 4);
	r = o93_inflate(&p, 30, QAWS_JOIN_MITER, QAWS_END_POLYGON, 0, 0, 0);
	ok = ok && o93_count(r) == 1 && o93_area(r, 0) < 0;
	qaws_clip64_result_destroy(r);
	r = o93_inflate(&p, 30, QAWS_JOIN_MITER, QAWS_END_POLYGON, 2, 0, QAWS_CLIP_REVERSE_SOLUTION);
	ok = ok && o93_count(r) == 1 && fabs(o93_area(r, 0)) > sub_area && o93_area(r, 0) > 0;
	qaws_clip64_result_destroy(r);
	r = o93_inflate(&p, -15, QAWS_JOIN_MITER, QAWS_END_POLYGON, 0, 0, 0);
	ok = ok && o93_count(r) == 0;
	qaws_clip64_result_destroy(r);
	sprintf(msg, "TestOffsets9 (#733): orientation kept, reversed on request, over-shrunk to nothing (%d)", ok);
	TEST_ASSERT(ok, msg);
	o93_free(&p);
}

static void test_offsets10(void) /* #715 */
{
	o93_paths p;
	qaws_clip64_result* r;
	char msg[200];
	memset(&p, 0, sizeof(p));
	o93_add(&p, off10_a, (unsigned int)(sizeof(off10_a) / sizeof(off10_a[0]) / 2));
	o93_add(&p, off10_b, (unsigned int)(sizeof(off10_b) / sizeof(off10_b[0]) / 2));
	r = o93_inflate(&p, -2212495.6382562499, QAWS_JOIN_ROUND, QAWS_END_POLYGON, 2, 104857.61318750000, 0);
	sprintf(msg, "TestOffsets10 (#715): 2 paths (%u)", o93_count(r));
	TEST_ASSERT(o93_count(r) == 2, msg);
	qaws_clip64_result_destroy(r);
	o93_free(&p);
}

static void test_offsets11_12_13(void) /* #405, #873, #965 */
{
	static long long const a[] = { -1,-1, -1,11, 11,11, 11,-1 };
	static long long const b[] = { 667680768,-36382704, 737202688,-87034880, 742581888,-86055680, 747603968,-84684800 };
	static long long const c[] = { 0,0, 0,10, 10,0 };
	static long long const pt[] = { 0,20 };
	o93_paths p;
	qaws_clip64_result *r1, *r2, *r3, *r4;
	char msg[300];
	memset(&p, 0, sizeof(p)); o93_add(&p, a, 4);
	r1 = o93_inflate(&p, -50, QAWS_JOIN_MITER, QAWS_END_POLYGON, 0, 0, 0); o93_free(&p);
	o93_add(&p, b, 4);
	r2 = o93_inflate(&p, -249561088, QAWS_JOIN_MITER, QAWS_END_POLYGON, 0, 0, 0); o93_free(&p);
	o93_add(&p, c, 3);
	r3 = o93_inflate(&p, 2, QAWS_JOIN_MITER, QAWS_END_POLYGON, 0, 0, 0);
	o93_add(&p, pt, 1);
	r4 = o93_inflate(&p, 2, QAWS_JOIN_MITER, QAWS_END_POLYGON, 0, 0, 0); o93_free(&p);
	sprintf(msg, "TestOffsets11/12 (#405, #873): over-shrunk to nothing (%u, %u); TestOffsets13 (#965): area 122 with and without a lone point (%.0f, %.0f)",
		o93_count(r1), o93_count(r2), r3 ? fabs(o93_total_area(r3)) : 0, r4 ? fabs(o93_total_area(r4)) : 0);
	TEST_ASSERT(o93_count(r1) == 0 && o93_count(r2) == 0 && r3 && fabs(o93_total_area(r3)) == 122 && r4 && fabs(o93_total_area(r4)) == 122, msg);
	qaws_clip64_result_destroy(r1); qaws_clip64_result_destroy(r2); qaws_clip64_result_destroy(r3); qaws_clip64_result_destroy(r4);
}

static void test_orientation(void)
{
	static long long const a[] = { 0,0, 0,5, 5,5, 5,0 };
	static long long const b0[] = { 20,220, 280,220, 280,280, 20,280 };
	static long long const b1[] = { 0,200, 0,300, 300,300, 300,200 };
	o93_paths p;
	qaws_clip64_result *r1, *r2;
	int ok;
	char msg[200];
	memset(&p, 0, sizeof(p)); o93_add(&p, a, 4);
	r1 = o93_inflate(&p, 1, QAWS_JOIN_ROUND, QAWS_END_POLYGON, 0, 0, 0);
	ok = o93_count(r1) == 1 && ((qaws_path64_area2(p.xy[0], 4) > 0) == (o93_area(r1, 0) > 0));
	o93_free(&p);
	o93_add(&p, b0, 4); o93_add(&p, b1, 4);
	r2 = o93_inflate(&p, 5, QAWS_JOIN_ROUND, QAWS_END_POLYGON, 0, 0, QAWS_CLIP_REVERSE_SOLUTION);
	/* the outer path (the larger) runs against subject[1] */
	{
		unsigned int i, outer = 0;
		for (i = 1; i < o93_count(r2); i++) if (fabs(o93_area(r2, i)) > fabs(o93_area(r2, outer))) outer = i;
		ok = ok && o93_count(r2) == 2 && ((qaws_path64_area2(p.xy[1], 4) > 0) != (o93_area(r2, outer) > 0));
	}
	sprintf(msg, "TestOffsettingOrientation1/2: the solution keeps the input's orientation, reversed on request (%d)", ok);
	TEST_ASSERT(ok, msg);
	qaws_clip64_result_destroy(r1); qaws_clip64_result_destroy(r2);
	o93_free(&p);
}

int test_93_offset_parity_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 93: offset parity with Clipper2\n");
	test_corpus();
	test_offsets2();
	test_offsets3();
	test_offsets4();
	test_offsets5();
	test_offsets6();
	test_offsets7();
	test_offsets8();
	test_offsets9();
	test_offsets10();
	test_offsets11_12_13();
	test_orientation();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
