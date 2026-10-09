/*
 * Test 81: Batched curve / surface intersection
 *
 *   - paraboloids z = x^2 + y^2 + c on [-1, 1]^2 (biquadratic Bezier patches)
 *     against random 3D lines: every root of the line's quadratic inside
 *     the patch found once, at the closed-form point
 *   - random 3D cubics: every hit of the pairwise
 *     qaws_surface_find_curve_intersections also found by the batch
 *   - timing: the batch against every curve / surface pair
 *   - the certified batch: dyadic lines against the same paraboloids, each
 *     closed-form root in exactly one enclosure, as the pairwise exact call
 */

#include "test_common.h"
#include "qaws_surface_batch.h"
#include "qaws_exact.h"
#include <math.h>
#include <string.h>
#include <time.h>

#define SBT_SURF 3

static double sbt_lift(unsigned int k) { return 0.75 * k; }

static qaws_surface* sbt_paraboloid(double c)
{
	static double const a[3] = { 1, -1, 1 };
	qaws_vec3 cp[9];
	qaws_surface_bezier_desc d;
	qaws_surface* s = NULL;
	unsigned int i, j;
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
		{
			cp[i * 3 + j].x = (qaws_scalar)(-1.0 + i);
			cp[i * 3 + j].y = (qaws_scalar)(-1.0 + j);
			cp[i * 3 + j].z = (qaws_scalar)(a[i] + a[j] + c);
		}
	memset(&d, 0, sizeof(d));
	d.u_degree = 2;
	d.v_degree = 2;
	d.control_points = cp;
	d.u_point_count = 3;
	d.v_point_count = 3;
	qaws_surface_create_bezier(&d, &s);
	return s;
}

static unsigned long long g_sbt_state = 0xD1B54A32D192ED03ull;
static double sbt_rand(void)
{
	g_sbt_state ^= g_sbt_state << 13; g_sbt_state ^= g_sbt_state >> 7; g_sbt_state ^= g_sbt_state << 17;
	return (double)(g_sbt_state >> 11) / 9007199254740992.0;
}

static qaws_curve* sbt_segment(double const* p, double const* q)
{
	qaws_scalar cp[6];
	qaws_bezier_desc d;
	qaws_curve* c = NULL;
	unsigned int k;
	for (k = 0; k < 3; k++)
	{
		cp[k] = (qaws_scalar)p[k];
		cp[3 + k] = (qaws_scalar)q[k];
	}
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_3D;
	d.degree = 1;
	d.control_points = cp;
	d.control_point_count = 2;
	qaws_curve_create_bezier(&d, &c);
	return c;
}

/* roots s in [0, 1] of |xy(s)|^2 + c - z(s) = 0 inside the patch; 0 when a
   root is too close to a border or a tangency (ambiguous for a test) */
static int sbt_roots(double const* p, double const* q, double c, double* s_out, unsigned int* n)
{
	double d[3], A, B, C, disc, r[2];
	unsigned int k, m = 0;
	for (k = 0; k < 3; k++)
		d[k] = q[k] - p[k];
	A = d[0] * d[0] + d[1] * d[1];
	B = 2 * (p[0] * d[0] + p[1] * d[1]) - d[2];
	C = p[0] * p[0] + p[1] * p[1] + c - p[2];
	disc = B * B - 4 * A * C;
	*n = 0;
	if (fabs(disc) < 1e-6)
		return 0;
	if (disc < 0)
		return 1;
	r[0] = (-B - sqrt(disc)) / (2 * A);
	r[1] = (-B + sqrt(disc)) / (2 * A);
	for (k = 0; k < 2; k++)
	{
		double x = p[0] + r[k] * d[0], y = p[1] + r[k] * d[1];
		if (fabs(r[k]) < 1e-4 || fabs(r[k] - 1) < 1e-4 || fabs(fabs(x) - 1) < 1e-4 || fabs(fabs(y) - 1) < 1e-4)
			return 0;
		if (r[k] > 0 && r[k] < 1 && fabs(x) < 1 && fabs(y) < 1)
			s_out[m++] = r[k];
	}
	*n = m;
	return 1;
}

#define SBT_LINES 300

static double sbt_now(void)
{
	return (double)clock() / CLOCKS_PER_SEC;
}

static void test_lines(void)
{
	qaws_surface* sf[SBT_SURF];
	qaws_curve* cs[SBT_LINES];
	static double P[SBT_LINES][3], Q[SBT_LINES][3];
	static qaws_curve_surface_batch_hit h[4096];
	qaws_curve_surface_batch_desc d;
	qaws_surface_batch_stats st;
	unsigned int i, k, j, n = 0, expected = 0, matched = 0, np = 0;
	double worst = 0, t0, tb, tp;
	qaws_status s;
	char msg[200];
	for (k = 0; k < SBT_SURF; k++)
		sf[k] = sbt_paraboloid(sbt_lift(k));
	for (i = 0; i < SBT_LINES; i++)
	{
		for (;;)
		{
			unsigned int c, ok = 1;
			double tmp[2];
			for (c = 0; c < 3; c++)
			{
				P[i][c] = c < 2 ? -1.3 + 2.6 * sbt_rand() : -0.5 + 4.5 * sbt_rand();
				Q[i][c] = c < 2 ? -1.3 + 2.6 * sbt_rand() : -0.5 + 4.5 * sbt_rand();
			}
			for (k = 0; k < SBT_SURF && ok; k++)
			{
				unsigned int m;
				ok = sbt_roots(P[i], Q[i], sbt_lift(k), tmp, &m);
			}
			if (ok)
				break;
		}
		cs[i] = sbt_segment(P[i], Q[i]);
	}
	memset(&d, 0, sizeof(d));
	d.curves = (qaws_curve const* const*)cs;
	d.curve_count = SBT_LINES;
	d.surfaces = (qaws_surface const* const*)sf;
	d.surface_count = SBT_SURF;
	t0 = sbt_now();
	s = qaws_curve_surface_batch_find_intersections(&d, h, 4096, &n, &st);
	tb = sbt_now() - t0;
	for (i = 0; i < SBT_LINES; i++)
		for (k = 0; k < SBT_SURF; k++)
		{
			double r[2];
			unsigned int m = 0, a;
			sbt_roots(P[i], Q[i], sbt_lift(k), r, &m);
			for (a = 0; a < m; a++)
			{
				unsigned int found = 0;
				double best = 1e30;
				expected++;
				for (j = 0; j < n; j++)
					if (h[j].curve == i && h[j].surface == k)
					{
						double e = fabs(h[j].t - r[a]);
						if (e < 1e-4)
							found++;
						if (e < best)
							best = e;
					}
				matched += found == 1;
				if (best > worst)
					worst = best;
			}
		}
	/* the pairwise call over every curve / surface pair */
	t0 = sbt_now();
	for (i = 0; i < SBT_LINES; i++)
		for (k = 0; k < SBT_SURF; k++)
		{
			qaws_surface_curve_intersection buf[16];
			unsigned int m = 0;
			qaws_surface_find_curve_intersections(sf[k], cs[i], buf, 16, &m);
			np += m;
		}
	tp = sbt_now() - t0;
	printf("    %u lines x %u paraboloids: %u segments, %u patches, %u cells, %u candidates, %u refined; %u hits in %.4f s (worst t error %.1e); pairwise %u hits in %.3f s\n",
		SBT_LINES, SBT_SURF, st.segment_count, st.patch_count, st.cell_count, st.candidate_count, st.newton_count, n, tb, worst, np, tp);
	sprintf(msg, "every one of the %u closed-form line / paraboloid hits found once", expected);
	TEST_ASSERT(s == QAWS_STATUS_OK && matched == expected && n == expected, msg);
	TEST_ASSERT(worst < (QAWS_SCALAR_IS_FLOAT ? 1e-4 : 1e-12), "hits agree with the quadratic's roots");
	for (i = 0; i < SBT_LINES; i++)
		qaws_curve_destroy(cs[i]);
	for (k = 0; k < SBT_SURF; k++)
		qaws_surface_destroy(sf[k]);
}

#define SBT_CUBICS 40

static void test_cubics(void)
{
	qaws_surface* sf[SBT_SURF];
	qaws_curve* cs[SBT_CUBICS];
	static qaws_curve_surface_batch_hit h[4096];
	qaws_curve_surface_batch_desc d;
	unsigned int i, k, j, n = 0, total = 0, found = 0;
	char msg[200];
	for (k = 0; k < SBT_SURF; k++)
		sf[k] = sbt_paraboloid(sbt_lift(k));
	for (i = 0; i < SBT_CUBICS; i++)
	{
		qaws_scalar cp[12];
		qaws_bezier_desc bd;
		unsigned int c;
		for (c = 0; c < 12; c++)
			cp[c] = (qaws_scalar)(c % 3 < 2 ? -1.2 + 2.4 * sbt_rand() : -0.5 + 4.5 * sbt_rand());
		memset(&bd, 0, sizeof(bd));
		bd.dimension = QAWS_DIMENSION_3D;
		bd.degree = 3;
		bd.control_points = cp;
		bd.control_point_count = 4;
		cs[i] = NULL;
		qaws_curve_create_bezier(&bd, &cs[i]);
	}
	memset(&d, 0, sizeof(d));
	d.curves = (qaws_curve const* const*)cs;
	d.curve_count = SBT_CUBICS;
	d.surfaces = (qaws_surface const* const*)sf;
	d.surface_count = SBT_SURF;
	qaws_curve_surface_batch_find_intersections(&d, h, 4096, &n, NULL);
	for (i = 0; i < SBT_CUBICS; i++)
		for (k = 0; k < SBT_SURF; k++)
		{
			qaws_surface_curve_intersection buf[16];
			unsigned int m = 0, a;
			qaws_surface_find_curve_intersections(sf[k], cs[i], buf, 16, &m);
			for (a = 0; a < m && a < 16; a++)
			{
				total++;
				for (j = 0; j < n; j++)
					if (h[j].curve == i && h[j].surface == k && fabs(h[j].position.x - buf[a].position.x) + fabs(h[j].position.y - buf[a].position.y)
						+ fabs(h[j].position.z - buf[a].position.z) < 1e-3)
					{
						found++;
						break;
					}
			}
		}
	printf("    %u random cubics x %u paraboloids: batch %u hits, pairwise %u, pairwise hits found by the batch %u\n", SBT_CUBICS, SBT_SURF, n, total, found);
	sprintf(msg, "batch finds every one of the %u pairwise curve / surface hits", total);
	TEST_ASSERT(found == total && n >= total, msg);
	for (i = 0; i < SBT_CUBICS; i++)
		qaws_curve_destroy(cs[i]);
	for (k = 0; k < SBT_SURF; k++)
		qaws_surface_destroy(sf[k]);
}

#define SBT_EXACT_LINES 60

/* the certified batch: dyadic lines (exact) against the dyadic paraboloids */
static void test_exact(void)
{
	qaws_surface* sf[SBT_SURF];
	qaws_exact_surface* es[SBT_SURF];
	qaws_curve* cs[SBT_EXACT_LINES];
	qaws_exact_curve* ec[SBT_EXACT_LINES];
	static double P[SBT_EXACT_LINES][3], Q[SBT_EXACT_LINES][3];
	static qaws_exact_curve_surface_batch_hit h[1024];
	qaws_exact_curve_surface_hit hp[16];
	qaws_exact_surface_batch_desc d;
	qaws_exact_batch_stats st;
	qaws_exact_desc ed;
	unsigned int i, k, j, n = 0, expected = 0, matched = 0, np = 0;
	double t0, tb, tp;
	qaws_status s;
	char msg[200];
	qaws_exact_desc_default(&ed);
	ed.space_exp2 = -10;
	for (k = 0; k < SBT_SURF; k++)
	{
		sf[k] = sbt_paraboloid(sbt_lift(k));
		qaws_exact_surface_prepare(&ed, sf[k], &es[k], NULL);
	}
	for (i = 0; i < SBT_EXACT_LINES; i++)
	{
		for (;;)
		{
			unsigned int c, ok = 1;
			double tmp[2];
			for (c = 0; c < 3; c++)
			{
				/* on the 2^-10 lattice */
				P[i][c] = ldexp(nearbyint(ldexp(c < 2 ? -1.3 + 2.6 * sbt_rand() : -0.5 + 4.5 * sbt_rand(), 10)), -10);
				Q[i][c] = ldexp(nearbyint(ldexp(c < 2 ? -1.3 + 2.6 * sbt_rand() : -0.5 + 4.5 * sbt_rand(), 10)), -10);
			}
			for (k = 0; k < SBT_SURF && ok; k++)
			{
				unsigned int m;
				ok = sbt_roots(P[i], Q[i], sbt_lift(k), tmp, &m);
			}
			if (ok)
				break;
		}
		cs[i] = sbt_segment(P[i], Q[i]);
		qaws_exact_curve_prepare(&ed, cs[i], &ec[i], NULL);
	}
	memset(&d, 0, sizeof(d));
	d.curves = (qaws_exact_curve const* const*)ec;
	d.curve_count = SBT_EXACT_LINES;
	d.surfaces = (qaws_exact_surface const* const*)es;
	d.surface_count = SBT_SURF;
	t0 = sbt_now();
	s = qaws_exact_curve_surface_batch_hits(&d, h, 1024, &n, &st);
	tb = sbt_now() - t0;
	for (i = 0; i < SBT_EXACT_LINES; i++)
		for (k = 0; k < SBT_SURF; k++)
		{
			double r[2];
			unsigned int m = 0, a;
			sbt_roots(P[i], Q[i], sbt_lift(k), r, &m);
			for (a = 0; a < m; a++)
			{
				unsigned int found = 0;
				expected++;
				for (j = 0; j < n; j++)
					if (h[j].curve == i && h[j].surface == k && h[j].hit.t_lo - 1e-12 <= r[a] && r[a] <= h[j].hit.t_hi + 1e-12)
						found++;
				matched += found == 1;
			}
		}
	t0 = sbt_now();
	for (i = 0; i < SBT_EXACT_LINES; i++)
		for (k = 0; k < SBT_SURF; k++)
		{
			unsigned int m = 0;
			qaws_exact_curve_surface_hits(ec[i], es[k], hp, 16, &m);
			np += m;
		}
	tp = sbt_now() - t0;
	printf("    exact: %u lines x %u paraboloids: %u spans, %u patches, %u candidates; %u certified hits in %.3f s; pairwise exact %u in %.3f s\n",
		SBT_EXACT_LINES, SBT_SURF, st.span_count, st.patch_count, st.candidate_count, n, tb, np, tp);
	sprintf(msg, "certified batch: each of the %u closed-form roots in exactly one enclosure", expected);
	TEST_ASSERT(s == QAWS_STATUS_OK && matched == expected && n == expected && np == n, msg);
	for (i = 0; i < SBT_EXACT_LINES; i++)
	{
		qaws_exact_curve_destroy(ec[i]);
		qaws_curve_destroy(cs[i]);
	}
	for (k = 0; k < SBT_SURF; k++)
	{
		qaws_exact_surface_destroy(es[k]);
		qaws_surface_destroy(sf[k]);
	}
}

static qaws_surface* sbt_plane(double h, double s)
{
	/* z = h + s x over [-1.5, 1.5]^2, a bilinear Bezier patch */
	qaws_vec3 cp[4];
	qaws_surface_bezier_desc d;
	qaws_surface* out = NULL;
	unsigned int i, j;
	for (i = 0; i < 2; i++)
		for (j = 0; j < 2; j++)
		{
			double x = i ? 1.5 : -1.5, y = j ? 1.5 : -1.5;
			cp[i * 2 + j].x = (qaws_scalar)x;
			cp[i * 2 + j].y = (qaws_scalar)y;
			cp[i * 2 + j].z = (qaws_scalar)(h + s * x);
		}
	memset(&d, 0, sizeof(d));
	d.u_degree = 1;
	d.v_degree = 1;
	d.control_points = cp;
	d.u_point_count = 2;
	d.v_point_count = 2;
	qaws_surface_create_bezier(&d, &out);
	return out;
}

#define SBT_PLANES 6

/* paraboloids z = x^2 + y^2 + c (family 0) against planes z = h + s x
   (family 1): circles of centre (s / 2, 0), radius^2 = h - c + s^2 / 4 */
static void test_ssi(void)
{
	static double const ph[SBT_PLANES] = { 0.30, 0.55, 0.80, 1.10, 1.35, 1.44 };
	static double const ps[SBT_PLANES] = { 0.0, 0.25, -0.3, 0.0, 0.2, 0.0 };
	qaws_surface* sf[SBT_SURF + SBT_PLANES];
	unsigned int fam[SBT_SURF + SBT_PLANES], i, k, nc = 0, np = 0, closed_expected = 0, closed_ok = 0, open_expected = 0, open_found = 0, bad_curves = 0;
	static qaws_surface_batch_curve cv[256];
	static qaws_ssi_point pt[1 << 16];
	qaws_surface_batch_desc d;
	qaws_surface_batch_stats st;
	double worst = 0, worst_gap = 0, t0, tb, tp;
	qaws_status s;
	char msg[240];
	for (k = 0; k < SBT_SURF; k++)
	{
		sf[k] = sbt_paraboloid(sbt_lift(k));
		fam[k] = 0;
	}
	for (k = 0; k < SBT_PLANES; k++)
	{
		sf[SBT_SURF + k] = sbt_plane(ph[k], ps[k]);
		fam[SBT_SURF + k] = 1;
	}
	memset(&d, 0, sizeof(d));
	d.surfaces = (qaws_surface const* const*)sf;
	d.surface_count = SBT_SURF + SBT_PLANES;
	d.families = fam;
	t0 = sbt_now();
	s = qaws_surface_batch_find_intersections(&d, cv, 256, &nc, pt, 1 << 16, &np, &st);
	tb = sbt_now() - t0;
	for (i = 0; i < SBT_SURF; i++)
		for (k = 0; k < SBT_PLANES; k++)
		{
			double cx = ps[k] / 2, r2 = ph[k] - sbt_lift(i) + ps[k] * ps[k] / 4, r = r2 > 0 ? sqrt(r2) : 0;
			unsigned int c, ncurve = 0, nopen = 0;
			for (c = 0; c < nc && c < 256; c++)
			{
				unsigned int a;
				double ang[4096], gap = 0;
				unsigned int na = 0;
				if (cv[c].surface_a != i || cv[c].surface_b != SBT_SURF + k)
					continue;
				ncurve++;
				nopen += !cv[c].closed;
				for (a = 0; a < cv[c].count; a++)
				{
					qaws_ssi_point const* p = &pt[cv[c].first + a];
					double e = fabs(hypot(p->position.x - cx, p->position.y) - r);
					if (e > worst)
						worst = e;
					if (na < 4096)
						ang[na++] = atan2(p->position.y, p->position.x - cx);
				}
				if (cv[c].closed)
				{
					/* largest angular step along the closed polyline */
					for (a = 0; a < na; a++)
					{
						double g = fabs(ang[(a + 1) % na] - ang[a]);
						if (g > 3.14159265358979)
							g = 2 * 3.14159265358979 - g;
						if (g > gap)
							gap = g;
					}
					if (gap > worst_gap)
						worst_gap = gap;
				}
			}
			/* inside the paraboloid patch [-1, 1]^2: one closed curve */
			if (r > 0 && fabs(cx) + r < 0.95)
			{
				closed_expected++;
				closed_ok += ncurve == 1 && nopen == 0;
			}
			else if (r > 1.02 && r < 1.38 && fabs(cx) < 1e-12)
			{
				/* the circle leaves the patch [-1, 1]^2 through its four sides: four corner arcs */
				open_expected += 4;
				open_found += ncurve == 4 && nopen == 4 ? 4 : 0;
			}
			else if (r == 0)
				bad_curves += ncurve;
		}
	t0 = sbt_now();
	{
		unsigned int pairs = 0;
		static qaws_ssi_curve pc[64];
		static qaws_ssi_point pp[1 << 14];
		for (i = 0; i < SBT_SURF; i++)
			for (k = 0; k < SBT_PLANES; k++)
			{
				qaws_ssi_desc sd;
				unsigned int n = 0;
				memset(&sd, 0, sizeof(sd));
				sd.surface_a = sf[i];
				sd.surface_b = sf[SBT_SURF + k];
				qaws_surface_intersect(&sd, pc, 64, &n, pp, 1 << 14);
				pairs += n;
			}
		tp = sbt_now() - t0;
		printf("    %u paraboloids x %u planes: %u patches, %u candidates, %u refined; %u curves, %u points in %.4f s (worst radius error %.1e, widest angular step %.3f); pairwise marching %u curves in %.3f s\n",
			SBT_SURF, SBT_PLANES, st.patch_count, st.candidate_count, st.newton_count, nc, np, tb, worst, worst_gap, pairs, tp);
	}
	sprintf(msg, "every interior circle is one closed curve (%u of %u), the boundary-crossing circles four arcs (%u of %u)", closed_ok, closed_expected, open_found, open_expected);
	TEST_ASSERT(s == QAWS_STATUS_OK && closed_ok == closed_expected && open_found == open_expected && bad_curves == 0, msg);
	TEST_ASSERT(worst < (QAWS_SCALAR_IS_FLOAT ? 1e-4 : 1e-10) && worst_gap < 0.35, "points on the circles to rounding, no gap along them");
	for (k = 0; k < SBT_SURF + SBT_PLANES; k++)
		qaws_surface_destroy(sf[k]);
}

#define SBT_EPLANES 7

/* the certified surface batch: dyadic paraboloids x dyadic horizontal planes */
static void test_exact_ssi(void)
{
	static double const ph[SBT_EPLANES] = { 0.25, 0.5, 1.0625, 1.25, 1.5, 1.6875, 1.9375 };
	qaws_surface* sf[SBT_SURF + SBT_EPLANES];
	qaws_exact_surface* es[SBT_SURF + SBT_EPLANES];
	unsigned int fam[SBT_SURF + SBT_EPLANES], i, k, np = 0, nb = 0, closed_expected = 0, closed_ok = 0, outside = 0, same = 1, pair_points = 0;
	static qaws_exact_ssi_point pts[1 << 14], pp[1 << 13];
	static qaws_exact_ssi_batch_branch br[512];
	static qaws_exact_ssi_branch pb[64];
	qaws_exact_ssi_batch_desc d;
	qaws_exact_batch_stats st;
	qaws_exact_desc ed;
	double t0, tb, tp;
	qaws_status s;
	char msg[240];
	qaws_exact_desc_default(&ed);
	ed.space_exp2 = -10;
	for (k = 0; k < SBT_SURF + SBT_EPLANES; k++)
	{
		sf[k] = k < SBT_SURF ? sbt_paraboloid(sbt_lift(k)) : sbt_plane(ph[k - SBT_SURF], 0);
		fam[k] = k >= SBT_SURF;
		qaws_exact_surface_prepare(&ed, sf[k], &es[k], NULL);
	}
	memset(&d, 0, sizeof(d));
	d.surfaces = (qaws_exact_surface const* const*)es;
	d.surface_count = SBT_SURF + SBT_EPLANES;
	d.families = fam;
	t0 = sbt_now();
	s = qaws_exact_surface_batch_hits(&d, pts, 1 << 14, &np, br, 512, &nb, &st);
	tb = sbt_now() - t0;
	for (i = 0; i < SBT_SURF; i++)
		for (k = 0; k < SBT_EPLANES; k++)
		{
			double r2 = ph[k] - sbt_lift(i);
			unsigned int c, nbr = 0, nclosed = 0, a;
			for (c = 0; c < nb; c++)
			{
				if (br[c].surface_a != i || br[c].surface_b != SBT_SURF + k)
					continue;
				nbr++;
				nclosed += br[c].branch.closed != 0;
				for (a = 0; a < br[c].branch.count; a++)
				{
					/* x = 2 u - 1, y = 2 v - 1 on the paraboloid: x^2 + y^2 must reach r^2 */
					qaws_exact_ssi_point const* p = &pts[br[c].branch.first + a];
					double x0 = 2 * p->u1_lo - 1, x1 = 2 * p->u1_hi - 1, y0 = 2 * p->v1_lo - 1, y1 = 2 * p->v1_hi - 1;
					double xl = x0 <= 0 && x1 >= 0 ? 0 : fmin(x0 * x0, x1 * x1), xh = fmax(x0 * x0, x1 * x1);
					double yl = y0 <= 0 && y1 >= 0 ? 0 : fmin(y0 * y0, y1 * y1), yh = fmax(y0 * y0, y1 * y1);
					if (xl + yl > r2 + 1e-9 || xh + yh < r2 - 1e-9)
						outside++;
				}
			}
			if (r2 > 0 && r2 < 0.95)
			{
				closed_expected++;
				closed_ok += nbr == 1 && nclosed == 1;
			}
		}
	/* pair by pair, every paraboloid / plane pair */
	t0 = sbt_now();
	for (i = 0; i < SBT_SURF; i++)
		for (k = 0; k < SBT_EPLANES; k++)
		{
			unsigned int n = 0, m = 0, c, bp = 0, bb = 0;
			qaws_status s2 = qaws_exact_surface_surface_hits(es[i], es[SBT_SURF + k], 0, pp, 1 << 13, &n, pb, 64, &m);
			if (s2 != QAWS_STATUS_OK)
				continue;
			pair_points += n;
			for (c = 0; c < nb; c++)
				if (br[c].surface_a == i && br[c].surface_b == SBT_SURF + k)
				{
					bb++;
					bp += br[c].branch.count;
				}
			same &= bb == m && bp == n;
		}
	tp = sbt_now() - t0;
	printf("    exact: %u paraboloids x %u planes: %u patches, %u candidate patch pairs; %u branches, %u certified points in %.3f s (status %d, %u uncertified pair); pairwise exact %u points in %.3f s\n",
		SBT_SURF, SBT_EPLANES, st.patch_count, st.candidate_count, nb, np, tb, (int)s, st.uncertified_count, pair_points, tp);
	sprintf(msg, "certified surface batch: %u of %u interior circles one closed branch, every point box on its circle, the tangent pair alone uncertified", closed_ok, closed_expected);
	TEST_ASSERT(s == QAWS_STATUS_CERTIFICATION_FAILED && st.uncertified_count == 1 && closed_ok == closed_expected && outside == 0, msg);
	TEST_ASSERT(same && pair_points == np, "the same branches and points as the pairwise certified call");
	for (k = 0; k < SBT_SURF + SBT_EPLANES; k++)
	{
		qaws_exact_surface_destroy(es[k]);
		qaws_surface_destroy(sf[k]);
	}
}

#define SBT_FRAMES 4

/* paraboloids prepared once, against planes that change and lines */
static void test_surface_sets(void)
{
	qaws_surface* par[SBT_SURF];
	qaws_surface* all[SBT_SURF + SBT_PLANES];
	unsigned int fam[SBT_SURF + SBT_PLANES], i, k, f, same = 1, same_cs, total = 0;
	static qaws_surface_batch_curve cs1[256], cs2[256];
	static qaws_ssi_point ps1[1 << 15], ps2[1 << 15];
	qaws_surface_set* pset = NULL;
	qaws_surface_batch_desc d;
	char msg[200];
	for (k = 0; k < SBT_SURF; k++)
		par[k] = sbt_paraboloid(sbt_lift(k));
	memset(&d, 0, sizeof(d));
	d.surfaces = (qaws_surface const* const*)par;
	d.surface_count = SBT_SURF;
	d.flatness = (qaws_scalar)0.01;
	qaws_surface_set_create(&d, &pset);
	for (f = 0; f < SBT_FRAMES; f++)
	{
		qaws_surface* pl[SBT_PLANES];
		qaws_surface_set* qset = NULL;
		qaws_surface_batch_desc qd;
		unsigned int n1 = 0, p1 = 0, n2 = 0, p2 = 0, c;
		for (k = 0; k < SBT_PLANES; k++)
			pl[k] = sbt_plane(0.3 + 0.25 * k + 0.07 * f, 0.1 * f - 0.15);
		memset(&qd, 0, sizeof(qd));
		qd.surfaces = (qaws_surface const* const*)pl;
		qd.surface_count = SBT_PLANES;
		qd.flatness = (qaws_scalar)0.01;
		qaws_surface_set_create(&qd, &qset);
		qaws_surface_set_find_intersections(pset, qset, cs1, 256, &n1, ps1, 1 << 15, &p1, NULL);
		/* the one-shot call on everything, two families */
		for (i = 0; i < SBT_SURF + SBT_PLANES; i++)
		{
			all[i] = i < SBT_SURF ? par[i] : pl[i - SBT_SURF];
			fam[i] = i >= SBT_SURF;
		}
		memset(&qd, 0, sizeof(qd));
		qd.surfaces = (qaws_surface const* const*)all;
		qd.surface_count = SBT_SURF + SBT_PLANES;
		qd.families = fam;
		qd.flatness = (qaws_scalar)0.01;
		qaws_surface_batch_find_intersections(&qd, cs2, 256, &n2, ps2, 1 << 15, &p2, NULL);
		same &= n1 == n2 && p1 == p2;
		for (c = 0; c < n1 && c < n2 && same; c++)
			same &= cs1[c].surface_a == cs2[c].surface_a && cs1[c].surface_b + SBT_SURF == cs2[c].surface_b && cs1[c].count == cs2[c].count
				&& cs1[c].closed == cs2[c].closed;
		for (c = 0; c < p1 && c < p2 && same; c++)
			same &= fabs(ps1[c].position.x - ps2[c].position.x) < 1e-12 && fabs(ps1[c].position.z - ps2[c].position.z) < 1e-12;
		total += n1;
		qaws_surface_set_destroy(qset);
		for (k = 0; k < SBT_PLANES; k++)
			qaws_surface_destroy(pl[k]);
	}
	/* prepared lines against the prepared paraboloids */
	{
		qaws_curve* ln[40];
		qaws_curve_set* cset = NULL;
		qaws_curve_batch_desc cd;
		qaws_curve_surface_batch_desc od;
		static qaws_curve_surface_batch_hit h1[512], h2[512];
		unsigned int n1 = 0, n2 = 0, c;
		for (i = 0; i < 40; i++)
		{
			double p[3], q[3];
			unsigned int k2;
			for (k2 = 0; k2 < 3; k2++)
			{
				p[k2] = k2 < 2 ? -1.3 + 2.6 * sbt_rand() : -0.5 + 4.5 * sbt_rand();
				q[k2] = k2 < 2 ? -1.3 + 2.6 * sbt_rand() : -0.5 + 4.5 * sbt_rand();
			}
			ln[i] = sbt_segment(p, q);
		}
		memset(&cd, 0, sizeof(cd));
		cd.curves = (qaws_curve const* const*)ln;
		cd.curve_count = 40;
		cd.flatness = (qaws_scalar)0.01;
		qaws_curve_set_create(&cd, &cset);
		qaws_curve_set_find_surface_intersections(cset, pset, h1, 512, &n1, NULL);
		memset(&od, 0, sizeof(od));
		od.curves = (qaws_curve const* const*)ln;
		od.curve_count = 40;
		od.surfaces = (qaws_surface const* const*)par;
		od.surface_count = SBT_SURF;
		od.flatness = (qaws_scalar)0.01;
		qaws_curve_surface_batch_find_intersections(&od, h2, 512, &n2, NULL);
		same_cs = n1 == n2;
		for (c = 0; c < n1 && c < n2; c++)
			same_cs &= h1[c].curve == h2[c].curve && h1[c].surface == h2[c].surface && fabs(h1[c].t - h2[c].t) < 1e-12;
		printf("    prepared paraboloids (%u patches): %u frames of %u planes, %u curves as the one-shot calls: %s; prepared lines x paraboloids %u hits (one-shot %u)\n",
			qaws_surface_set_get_patch_count(pset), SBT_FRAMES, SBT_PLANES, total, same ? "yes" : "no", n1, n2);
		qaws_curve_set_destroy(cset);
		for (i = 0; i < 40; i++)
			qaws_curve_destroy(ln[i]);
	}
	sprintf(msg, "prepared surface sets equal the one-shot surface batch on all %u frames", SBT_FRAMES);
	TEST_ASSERT(same && total > 0, msg);
	TEST_ASSERT(same_cs, "a prepared curve set against a prepared surface set equals the one-shot curve / surface batch");
	qaws_surface_set_destroy(pset);
	for (k = 0; k < SBT_SURF; k++)
		qaws_surface_destroy(par[k]);
}

#define SBT_CQ 400
#define SBT_CS (SBT_SURF + 2)

/* closest points on paraboloids and tilted planes */
static void test_surface_closest(void)
{
	static double const ph[2] = { 2.5, -0.3 }, ps[2] = { 0.2, -0.1 };
	qaws_surface* sf[SBT_CS];
	static qaws_scalar pts[3 * SBT_CQ];
	static qaws_surface_batch_closest out[SBT_CQ], outs[SBT_CQ];
	qaws_surface_batch_closest_desc d;
	qaws_surface_batch_stats st;
	qaws_surface_set* set = NULL;
	qaws_surface_batch_desc sd;
	unsigned int i, k, a, b, never_worse = 0, close = 0, same = 1, pair_worse = 0;
	double t0, tb, tp, worst_gain = 0;
	char msg[200];
	for (k = 0; k < SBT_SURF; k++)
		sf[k] = sbt_paraboloid(sbt_lift(k));
	for (k = 0; k < 2; k++)
		sf[SBT_SURF + k] = sbt_plane(ph[k], ps[k]);
	for (i = 0; i < SBT_CQ; i++)
	{
		pts[3 * i] = (qaws_scalar)(-1.1 + 2.2 * sbt_rand());
		pts[3 * i + 1] = (qaws_scalar)(-1.1 + 2.2 * sbt_rand());
		pts[3 * i + 2] = (qaws_scalar)(-0.4 + 3.6 * sbt_rand());
	}
	memset(&d, 0, sizeof(d));
	d.surfaces = (qaws_surface const* const*)sf;
	d.surface_count = SBT_CS;
	d.points = pts;
	d.point_count = SBT_CQ;
	t0 = sbt_now();
	qaws_surface_batch_find_closest(&d, out, &st);
	tb = sbt_now() - t0;
	for (i = 0; i < SBT_CQ; i++)
	{
		double x = pts[3 * i], y = pts[3 * i + 1], z = pts[3 * i + 2], ref = 1e30;
		/* planes: the orthogonal distance (the projection stays inside the patch) */
		for (k = 0; k < 2; k++)
		{
			double dd = fabs(z - ph[k] - ps[k] * x) / sqrt(1 + ps[k] * ps[k]);
			if (dd < ref)
				ref = dd;
		}
		/* paraboloids: dense sampling of [-1, 1]^2 */
		for (k = 0; k < SBT_SURF; k++)
			for (a = 0; a <= 160; a++)
				for (b = 0; b <= 160; b++)
				{
					double px = -1 + 2.0 * a / 160, py = -1 + 2.0 * b / 160, pz = px * px + py * py + sbt_lift(k);
					double dd = sqrt((px - x) * (px - x) + (py - y) * (py - y) + (pz - z) * (pz - z));
					if (dd < ref)
						ref = dd;
				}
		never_worse += out[i].surface != QAWS_CURVE_BATCH_NONE && out[i].distance <= ref + (QAWS_SCALAR_IS_FLOAT ? 2e-5 : 1e-9);
		close += ref - out[i].distance < 2e-2;   /* the sampling step on the steep flanks */
		if (ref - out[i].distance > worst_gain)
			worst_gain = ref - out[i].distance;
	}
	/* the prepared set gives the same answers */
	memset(&sd, 0, sizeof(sd));
	sd.surfaces = (qaws_surface const* const*)sf;
	sd.surface_count = SBT_CS;
	qaws_surface_set_create(&sd, &set);
	qaws_surface_set_find_closest(set, pts, SBT_CQ, 0, outs, NULL);
	for (i = 0; i < SBT_CQ; i++)
		same &= outs[i].surface == out[i].surface && fabs(outs[i].distance - out[i].distance) < 1e-12;
	qaws_surface_set_destroy(set);
	/* every point against every surface with the pairwise call */
	t0 = sbt_now();
	for (i = 0; i < SBT_CQ; i++)
	{
		qaws_vec3 q, c;
		double best = 1e30;
		q.x = pts[3 * i]; q.y = pts[3 * i + 1]; q.z = pts[3 * i + 2];
		for (k = 0; k < SBT_CS; k++)
		{
			qaws_scalar u = 0, v = 0;
			if (qaws_surface_find_closest_point(sf[k], q, &u, &v, &c) == QAWS_STATUS_OK)
			{
				double dd = sqrt((c.x - q.x) * (c.x - q.x) + (c.y - q.y) * (c.y - q.y) + (c.z - q.z) * (c.z - q.z));
				if (dd < best)
					best = dd;
			}
		}
		pair_worse += best > out[i].distance + 1e-6;
	}
	tp = sbt_now() - t0;
	printf("    closest on surfaces: %u points x %u surfaces in %.4f s (%u patches, %u refined); never worse than dense sampling %u, within its resolution %u (largest gain %.1e); pairwise %.3f s, farther on %u points\n",
		SBT_CQ, SBT_CS, tb, st.patch_count, st.candidate_count, never_worse, close, worst_gain, tp, pair_worse);
	sprintf(msg, "closest points on surfaces: never farther than a dense sampling (%u of %u) and within its resolution", never_worse, SBT_CQ);
	TEST_ASSERT(never_worse == SBT_CQ && close == SBT_CQ, msg);
	TEST_ASSERT(same, "the prepared surface set gives the same closest points");
	for (k = 0; k < SBT_CS; k++)
		qaws_surface_destroy(sf[k]);
}

#define SBT_ECQ 120

/* certified closest points on the dyadic paraboloids and two planes, against the float batch */
static void test_exact_surface_closest(void)
{
	qaws_surface* sf[SBT_SURF + 2];
	qaws_exact_surface* es[SBT_SURF + 2];
	static double pts[3 * SBT_ECQ];
	static qaws_scalar fpts[3 * SBT_ECQ];
	static qaws_exact_surface_closest_point eo[SBT_ECQ];
	static qaws_surface_batch_closest fo[SBT_ECQ];
	qaws_exact_ssi_batch_desc d;
	qaws_exact_batch_stats st;
	qaws_surface_batch_closest_desc fd;
	qaws_exact_desc ed;
	unsigned int i, k, n = SBT_SURF + 2, agree = 0, certified = 0;
	double widest = 0, t0, te, tol = QAWS_SCALAR_IS_FLOAT ? 2e-5 : 1e-9;
	qaws_status s;
	char msg[200];
	for (k = 0; k < SBT_SURF; k++)
		sf[k] = sbt_paraboloid(sbt_lift(k));
	sf[SBT_SURF] = sbt_plane(-0.5, 0);
	sf[SBT_SURF + 1] = sbt_plane(3.75, 0);
	qaws_exact_desc_default(&ed);
	ed.space_exp2 = -16;
	for (k = 0; k < n; k++)
		qaws_exact_surface_prepare(&ed, sf[k], &es[k], NULL);
	for (i = 0; i < SBT_ECQ; i++)
	{
		pts[3 * i] = ldexp(nearbyint(ldexp(-1.1 + 2.2 * sbt_rand(), 10)), -10);
		pts[3 * i + 1] = ldexp(nearbyint(ldexp(-1.1 + 2.2 * sbt_rand(), 10)), -10);
		pts[3 * i + 2] = ldexp(nearbyint(ldexp(-0.4 + 3.6 * sbt_rand(), 10)), -10);
		for (k = 0; k < 3; k++)
			fpts[3 * i + k] = (qaws_scalar)pts[3 * i + k];
	}
	memset(&d, 0, sizeof(d));
	d.surfaces = (qaws_exact_surface const* const*)es;
	d.surface_count = n;
	t0 = sbt_now();
	s = qaws_exact_surface_batch_closest(&d, pts, SBT_ECQ, eo, &st);
	te = sbt_now() - t0;
	memset(&fd, 0, sizeof(fd));
	fd.surfaces = (qaws_surface const* const*)sf;
	fd.surface_count = n;
	fd.points = fpts;
	fd.point_count = SBT_ECQ;
	qaws_surface_batch_find_closest(&fd, fo, NULL);
	for (i = 0; i < SBT_ECQ; i++)
	{
		certified += eo[i].certified != 0;
		if (eo[i].distance_hi - eo[i].distance_lo > widest)
			widest = eo[i].distance_hi - eo[i].distance_lo;
		if (eo[i].surface == fo[i].surface && fo[i].distance >= eo[i].distance_lo - tol && fo[i].distance <= eo[i].distance_hi + tol)
			agree++;
	}
	printf("    certified closest on surfaces: %u points x %u surfaces, %u patches solved, %u certified, widest enclosure %.1e, in %.3f s (status %d); float batch inside %u\n",
		SBT_ECQ, n, st.candidate_count, certified, widest, te, (int)s, agree);
	sprintf(msg, "certified closest points on surfaces: every float answer inside its enclosure (%u of %u), all certified", agree, SBT_ECQ);
	TEST_ASSERT(s == QAWS_STATUS_OK && agree == SBT_ECQ && certified == SBT_ECQ && widest < 1e-9, msg);
	/* midway between the two planes: a tie */
	{
		double p[3] = { 0.25, -0.125, 1.625 };
		qaws_exact_surface const* planes[2];
		qaws_exact_surface_closest_point o;
		planes[0] = es[SBT_SURF];
		planes[1] = es[SBT_SURF + 1];
		d.surfaces = planes;
		d.surface_count = 2;
		s = qaws_exact_surface_batch_closest(&d, p, 1, &o, NULL);
		TEST_ASSERT(s == QAWS_STATUS_OK && !o.certified && o.distance_lo <= 2.125 && o.distance_hi >= 2.125,
			"a point midway between two planes: a tie, reported uncertified at the exact distance");
	}
	for (k = 0; k < n; k++)
	{
		qaws_exact_surface_destroy(es[k]);
		qaws_surface_destroy(sf[k]);
	}
}

#define SBT_RAYS 1000

/* first hit of a ray: the paraboloids z = x^2 + y^2 + c on [-1, 1]^2 and the
   plane z = -0.25 on [-1.5, 1.5]^2; 0 when a hit is too close to a border,
   a tangency or another surface's hit (ambiguous for a test) */
static int sbt_first_hit(double const* o, double const* d, double* t_out, int* s_out)
{
	double best = 1e30;
	int bs = -1, k;
	double second = 1e30;
	for (k = 0; k < SBT_SURF; k++)
	{
		double c = sbt_lift((unsigned int)k), A = d[0] * d[0] + d[1] * d[1], B = 2 * (o[0] * d[0] + o[1] * d[1]) - d[2];
		double C = o[0] * o[0] + o[1] * o[1] + c - o[2], disc = B * B - 4 * A * C, r[2];
		int j;
		if (fabs(disc) < 1e-6)
			return 0;
		if (disc < 0)
			continue;
		r[0] = (-B - sqrt(disc)) / (2 * A);
		r[1] = (-B + sqrt(disc)) / (2 * A);
		for (j = 0; j < 2; j++)
		{
			double x = o[0] + r[j] * d[0], y = o[1] + r[j] * d[1];
			if (r[j] < 0)
				continue;
			if (fabs(fabs(x) - 1) < 1e-4 || fabs(fabs(y) - 1) < 1e-4)
				return 0;
			if (fabs(x) < 1 && fabs(y) < 1)
			{
				if (r[j] < best) { second = best; best = r[j]; bs = k; }
				else if (r[j] < second) second = r[j];
			}
		}
	}
	{
		double t = (-0.25 - o[2]) / d[2], x = o[0] + t * d[0], y = o[1] + t * d[1];
		if (t >= 0 && fabs(x) < 1.5 && fabs(y) < 1.5)
		{
			if (fabs(fabs(x) - 1.5) < 1e-4 || fabs(fabs(y) - 1.5) < 1e-4)
				return 0;
			if (t < best) { second = best; best = t; bs = SBT_SURF; }
			else if (t < second) second = t;
		}
	}
	if (second - best < 1e-4)
		return 0;
	*t_out = best;
	*s_out = bs;
	return 1;
}

static void test_raycast(void)
{
	qaws_surface* sf[SBT_SURF + 1];
	static qaws_scalar org[3 * SBT_RAYS], dir[3 * SBT_RAYS];
	static double t_ref[SBT_RAYS];
	static int s_ref[SBT_RAYS];
	static qaws_surface_ray_hit h[SBT_RAYS], hs[SBT_RAYS];
	qaws_surface_ray_desc d;
	qaws_surface_batch_stats st;
	qaws_surface_set* set = NULL;
	qaws_surface_batch_desc sd;
	unsigned int i, k, ok = 0, same = 1, pair_ok = 0;
	double worst = 0, t0, tb, tp, tol = QAWS_SCALAR_IS_FLOAT ? 1e-4 : 1e-9;
	char msg[200];
	for (k = 0; k < SBT_SURF; k++)
		sf[k] = sbt_paraboloid(sbt_lift(k));
	sf[SBT_SURF] = sbt_plane(-0.25, 0);
	for (i = 0; i < SBT_RAYS; i++)
		for (;;)
		{
			double o[3], dv[3];
			o[0] = -1.4 + 2.8 * sbt_rand(); o[1] = -1.4 + 2.8 * sbt_rand(); o[2] = 4.5;
			dv[0] = 0.6 * (sbt_rand() - 0.5); dv[1] = 0.6 * (sbt_rand() - 0.5); dv[2] = -1;
			if (!sbt_first_hit(o, dv, &t_ref[i], &s_ref[i]))
				continue;
			for (k = 0; k < 3; k++)
			{
				org[3 * i + k] = (qaws_scalar)o[k];
				dir[3 * i + k] = (qaws_scalar)dv[k];
			}
			break;
		}
	memset(&d, 0, sizeof(d));
	d.surfaces = (qaws_surface const* const*)sf;
	d.surface_count = SBT_SURF + 1;
	d.origins = org;
	d.directions = dir;
	d.ray_count = SBT_RAYS;
	t0 = sbt_now();
	qaws_surface_batch_raycast(&d, h, &st);
	tb = sbt_now() - t0;
	for (i = 0; i < SBT_RAYS; i++)
	{
		int good = s_ref[i] < 0 ? h[i].surface == QAWS_CURVE_BATCH_NONE : (h[i].surface == (unsigned int)s_ref[i] && fabs(h[i].t - t_ref[i]) < tol);
		ok += good;
		if (s_ref[i] >= 0 && h[i].surface == (unsigned int)s_ref[i] && fabs(h[i].t - t_ref[i]) > worst)
			worst = fabs(h[i].t - t_ref[i]);
	}
	/* the prepared set gives the same hits */
	memset(&sd, 0, sizeof(sd));
	sd.surfaces = (qaws_surface const* const*)sf;
	sd.surface_count = SBT_SURF + 1;
	qaws_surface_set_create(&sd, &set);
	qaws_surface_set_raycast(set, org, dir, SBT_RAYS, 0, hs, NULL);
	for (i = 0; i < SBT_RAYS; i++)
		same &= hs[i].surface == h[i].surface && (h[i].surface == QAWS_CURVE_BATCH_NONE || fabs(hs[i].t - h[i].t) < 1e-12);
	qaws_surface_set_destroy(set);
	/* the pairwise curve / surface call on segment curves, a tenth of the rays */
	t0 = sbt_now();
	for (i = 0; i < SBT_RAYS / 10; i++)
	{
		double p[3], q[3], best = 1e30;
		qaws_curve* seg;
		for (k = 0; k < 3; k++)
		{
			p[k] = org[3 * i + k];
			q[k] = org[3 * i + k] + 10 * dir[3 * i + k];
		}
		seg = sbt_segment(p, q);
		for (k = 0; k < SBT_SURF + 1; k++)
		{
			qaws_surface_curve_intersection buf[8];
			unsigned int m = 0, j;
			qaws_surface_find_curve_intersections(sf[k], seg, buf, 8, &m);
			for (j = 0; j < m && j < 8; j++)
				if (10 * buf[j].t < best)
					best = 10 * buf[j].t;
		}
		pair_ok += s_ref[i] < 0 ? best > 1e29 : fabs(best - t_ref[i]) < 1e-6;
		qaws_curve_destroy(seg);
	}
	tp = (sbt_now() - t0) * 10;
	printf("    ray casting: %u rays x %u surfaces in %.4f s (%u patches, %u met, %u refined), first hit right on %u (worst t error %.1e); pairwise (extrapolated) %.3f s, right on %u of %u\n",
		SBT_RAYS, SBT_SURF + 1, tb, st.patch_count, st.candidate_count, st.newton_count, ok, worst, tp, pair_ok, SBT_RAYS / 10);
	sprintf(msg, "ray casting: the first hit of all %u rays (surface and closed-form t)", SBT_RAYS);
	TEST_ASSERT(ok == SBT_RAYS, msg);
	TEST_ASSERT(same, "the prepared surface set casts the same rays");
	for (k = 0; k <= SBT_SURF; k++)
		qaws_surface_destroy(sf[k]);
}

#define SBT_ERAYS 300

/* certified ray casting: dyadic rays against the exact paraboloids and plane */
static void test_exact_raycast(void)
{
	qaws_surface* sf[SBT_SURF + 1];
	qaws_exact_surface* es[SBT_SURF + 1];
	static double p0[3 * SBT_ERAYS], p1[3 * SBT_ERAYS], t_ref[SBT_ERAYS];
	static int s_ref[SBT_ERAYS];
	static qaws_exact_ray_hit h[SBT_ERAYS];
	qaws_exact_surface_hit sh[16];
	qaws_exact_ssi_batch_desc d;
	qaws_exact_batch_stats st;
	qaws_exact_desc ed;
	unsigned int i, k, ok = 0, certified = 0, misses = 0, pair_ok = 0;
	double t0, tb, tp;
	qaws_status s;
	char msg[200];
	for (k = 0; k < SBT_SURF; k++)
		sf[k] = sbt_paraboloid(sbt_lift(k));
	sf[SBT_SURF] = sbt_plane(-0.25, 0);
	qaws_exact_desc_default(&ed);
	ed.space_exp2 = -16;
	for (k = 0; k <= SBT_SURF; k++)
		qaws_exact_surface_prepare(&ed, sf[k], &es[k], NULL);
	for (i = 0; i < SBT_ERAYS; i++)
		for (;;)
		{
			double o[3], dv[3];
			o[0] = ldexp(nearbyint(ldexp(-1.4 + 2.8 * sbt_rand(), 10)), -10);
			o[1] = ldexp(nearbyint(ldexp(-1.4 + 2.8 * sbt_rand(), 10)), -10);
			o[2] = 4.5;
			dv[0] = ldexp(nearbyint(ldexp(0.8 * (sbt_rand() - 0.5), 10)), -10);
			dv[1] = ldexp(nearbyint(ldexp(0.8 * (sbt_rand() - 0.5), 10)), -10);
			dv[2] = -1;
			if (!sbt_first_hit(o, dv, &t_ref[i], &s_ref[i]))
				continue;
			for (k = 0; k < 3; k++)
			{
				p0[3 * i + k] = o[k];
				p1[3 * i + k] = o[k] + dv[k];   /* exact: dyadic */
			}
			break;
		}
	memset(&d, 0, sizeof(d));
	d.surfaces = (qaws_exact_surface const* const*)es;
	d.surface_count = SBT_SURF + 1;
	t0 = sbt_now();
	s = qaws_exact_surface_batch_raycast(&d, p0, p1, SBT_ERAYS, 0, h, &st);
	tb = sbt_now() - t0;
	for (i = 0; i < SBT_ERAYS; i++)
	{
		certified += h[i].certified != 0;
		if (s_ref[i] < 0)
		{
			misses++;
			ok += h[i].surface == QAWS_EXACT_CLOSEST_NONE;
		}
		else
			ok += h[i].surface == (unsigned int)s_ref[i] && h[i].hit.t_lo - 1e-12 <= t_ref[i] && t_ref[i] <= h[i].hit.t_hi + 1e-12;
	}
	/* every surface with the single-ray call, a third of the rays */
	t0 = sbt_now();
	for (i = 0; i < SBT_ERAYS / 3; i++)
	{
		double best = 1e30;
		for (k = 0; k <= SBT_SURF; k++)
		{
			unsigned int n = 0, j;
			qaws_exact_surface_line_hits(es[k], p0 + 3 * i, p1 + 3 * i, sh, 16, &n);
			for (j = 0; j < n && j < 16; j++)
				if (sh[j].t_hi >= 0 && sh[j].t_lo < best)
					best = sh[j].t_lo;
		}
		pair_ok += s_ref[i] < 0 ? best > 1e29 : fabs(best - t_ref[i]) < 1e-9;
	}
	tp = (sbt_now() - t0) * 3;
	printf("    certified ray casting: %u rays (%u misses) x %u surfaces in %.3f s, %u patches solved; first hit right on %u, certified %u; single-ray calls (extrapolated) %.3f s\n",
		SBT_ERAYS, misses, SBT_SURF + 1, tb, st.candidate_count, ok, certified, tp);
	sprintf(msg, "certified ray casting: every first hit right and certified (%u, %u of %u); the single-ray calls agree on %u of %u", ok, certified, SBT_ERAYS,
		pair_ok, SBT_ERAYS / 3);
	TEST_ASSERT(s == QAWS_STATUS_OK && ok == SBT_ERAYS && certified == SBT_ERAYS && pair_ok == SBT_ERAYS / 3, msg);
	/* many patches: a bicubic height surface of 9 x 9 patches, the batch against the single-ray call */
	{
		static qaws_vec3 cp[12 * 12];
		qaws_surface_bspline_desc bd;
		qaws_surface* hs = NULL;
		qaws_exact_surface* ehs = NULL;
		qaws_exact_surface const* one[1];
		static double r0[300], r1[300];
		static qaws_exact_ray_hit hr[100];
		unsigned int a, b, nsame = 0, ncert = 0;
		double tb2, ts2;
		for (a = 0; a < 12; a++)
			for (b = 0; b < 12; b++)
			{
				cp[a * 12 + b].x = (qaws_scalar)(-1.5 + 3.0 * a / 11);
				cp[a * 12 + b].y = (qaws_scalar)(-1.5 + 3.0 * b / 11);
				cp[a * 12 + b].z = (qaws_scalar)ldexp(nearbyint(ldexp(0.4 * sin(1.7 * a) * cos(1.3 * b), 8)), -8);
			}
		memset(&bd, 0, sizeof(bd));
		bd.u_degree = 3;
		bd.v_degree = 3;
		bd.control_points = cp;
		bd.u_point_count = 12;
		bd.v_point_count = 12;
		qaws_surface_create_bspline(&bd, &hs);
		qaws_exact_surface_prepare(&ed, hs, &ehs, NULL);
		for (i = 0; i < 100; i++)
		{
			r0[3 * i] = ldexp(nearbyint(ldexp(-1.2 + 2.4 * sbt_rand(), 10)), -10);
			r0[3 * i + 1] = ldexp(nearbyint(ldexp(-1.2 + 2.4 * sbt_rand(), 10)), -10);
			r0[3 * i + 2] = 2;
			r1[3 * i] = r0[3 * i] + ldexp(nearbyint(ldexp(0.4 * (sbt_rand() - 0.5), 10)), -10);
			r1[3 * i + 1] = r0[3 * i + 1] + ldexp(nearbyint(ldexp(0.4 * (sbt_rand() - 0.5), 10)), -10);
			r1[3 * i + 2] = 1;
		}
		one[0] = ehs;
		d.surfaces = one;
		d.surface_count = 1;
		t0 = sbt_now();
		qaws_exact_surface_batch_raycast(&d, r0, r1, 100, 0, hr, &st);
		tb2 = sbt_now() - t0;
		t0 = sbt_now();
		for (i = 0; i < 100; i++)
		{
			unsigned int n = 0, j;
			double best = 1e30;
			qaws_exact_surface_line_hits(ehs, r0 + 3 * i, r1 + 3 * i, sh, 16, &n);
			for (j = 0; j < n && j < 16; j++)
				if (sh[j].t_hi >= 0 && sh[j].t_lo < best)
					best = sh[j].t_lo;
			/* the same first hit, or both miss */
			nsame += best > 1e29 ? hr[i].surface == QAWS_EXACT_CLOSEST_NONE : (hr[i].surface == 0 && hr[i].hit.t_lo == best);
			ncert += hr[i].certified != 0;
		}
		ts2 = sbt_now() - t0;
		printf("    certified rays on a %u-patch bicubic surface: 100 rays in %.3f s (%u patches solved), single-ray calls %.3f s; same first hit %u, certified %u\n",
			st.patch_count, tb2, st.candidate_count, ts2, nsame, ncert);
		TEST_ASSERT(nsame == 100 && ncert == 100, "many patches: the batch finds the single-ray call's first hit, certified, solving only the patches met");
		qaws_exact_surface_destroy(ehs);
		qaws_surface_destroy(hs);
		d.surfaces = (qaws_exact_surface const* const*)es;
		d.surface_count = SBT_SURF + 1;
	}
	/* rays upward from above everything: certified misses */
	{
		double q0[6] = { 0.25, -0.5, 4.5, -1.25, 0.75, 4.5 }, q1[6] = { 0.375, -0.5, 5.5, -1.25, 1.0, 5.5 };
		qaws_exact_ray_hit m[2];
		s = qaws_exact_surface_batch_raycast(&d, q0, q1, 2, 0, m, NULL);
		TEST_ASSERT(s == QAWS_STATUS_OK && m[0].surface == QAWS_EXACT_CLOSEST_NONE && m[1].surface == QAWS_EXACT_CLOSEST_NONE && m[0].certified && m[1].certified,
			"rays leaving the surfaces: certified misses");
	}
	for (k = 0; k <= SBT_SURF; k++)
	{
		qaws_exact_surface_destroy(es[k]);
		qaws_surface_destroy(sf[k]);
	}
}

/* an executor running the chunks one by one in a scrambled order */
static void sbt_shuffle_for(void* user, unsigned int count, qaws_batch_task_fn task, void* ctx)
{
	unsigned int k, step = 7;
	(void)user;
	while (count > 1 && count % step == 0)
		step += 2;
	for (k = 0; k < count; k++)
	{
		unsigned int c = (k * step + 3) % count;
		task(ctx, c, c + 1);
	}
}

#define SBT_THREADS 4

typedef struct sbt_thr
{
	qaws_batch_task_fn task;
	void* ctx;
	unsigned int begin, end;
} sbt_thr;

#ifdef _WIN32
#include <windows.h>
#include <process.h>
static unsigned __stdcall sbt_thr_main(void* p)
{
	sbt_thr* a = (sbt_thr*)p;
	a->task(a->ctx, a->begin, a->end);
	return 0;
}
#else
#include <pthread.h>
static void* sbt_thr_main(void* p)
{
	sbt_thr* a = (sbt_thr*)p;
	a->task(a->ctx, a->begin, a->end);
	return NULL;
}
#endif

/* a real executor: SBT_THREADS threads, contiguous ranges */
static void sbt_threads_for(void* user, unsigned int count, qaws_batch_task_fn task, void* ctx)
{
	sbt_thr arg[SBT_THREADS];
	unsigned int t;
#ifdef _WIN32
	HANDLE h[SBT_THREADS];
#else
	pthread_t h[SBT_THREADS];
#endif
	(void)user;
	for (t = 0; t < SBT_THREADS; t++)
	{
		arg[t].task = task;
		arg[t].ctx = ctx;
		arg[t].begin = count * t / SBT_THREADS;
		arg[t].end = count * (t + 1) / SBT_THREADS;
	}
	for (t = 1; t < SBT_THREADS; t++)
#ifdef _WIN32
		h[t] = (HANDLE)_beginthreadex(NULL, 0, sbt_thr_main, &arg[t], 0, NULL);
#else
		pthread_create(&h[t], NULL, sbt_thr_main, &arg[t]);
#endif
	task(ctx, arg[0].begin, arg[0].end);
	for (t = 1; t < SBT_THREADS; t++)
	{
#ifdef _WIN32
		WaitForSingleObject(h[t], INFINITE);
		CloseHandle(h[t]);
#else
		pthread_join(h[t], NULL);
#endif
	}
}

static qaws_scalar sbt_field(void* user, qaws_scalar const* p, qaws_scalar* g)
{
	(void)user;
	if (g)
	{
		g[0] = 2 * p[0];
		g[1] = 4 * p[1];
	}
	return p[0] * p[0] + 2 * p[1] * p[1];
}

/* every executor gives the serial results; threads give them faster */
static void test_executors(void)
{
	qaws_surface* sf[SBT_SURF + 1];
	qaws_curve* ring[24];
	unsigned int n = 20000, i, k, e, same[2] = { 1, 1 };
	qaws_scalar* pts = (qaws_scalar*)malloc(3 * n * sizeof(qaws_scalar));
	qaws_scalar* org = (qaws_scalar*)malloc(3 * n * sizeof(qaws_scalar));
	qaws_scalar* dir = (qaws_scalar*)malloc(3 * n * sizeof(qaws_scalar));
	qaws_scalar* pts2 = (qaws_scalar*)malloc(2 * n * sizeof(qaws_scalar));
	qaws_scalar lv[16];
	qaws_surface_batch_closest* sc[3];
	qaws_surface_ray_hit* rh[3];
	qaws_closest_point* cc[3];
	qaws_level_crossing* lc[3];
	unsigned int nlc[3];
	double tser = 0, tthr = 0;
	qaws_batch_executor ex[2];
	char msg[200];
	ex[0].parallel_for = sbt_shuffle_for; ex[0].user = NULL;
	ex[1].parallel_for = sbt_threads_for; ex[1].user = NULL;
	for (k = 0; k < SBT_SURF; k++)
		sf[k] = sbt_paraboloid(sbt_lift(k));
	sf[SBT_SURF] = sbt_plane(-0.25, 0);
	for (k = 0; k < 24; k++)
	{
		/* parabolas y = x^2 + c, c from -1.15 up */
		qaws_scalar co[6] = { 0, 0, 1, 0, 0, 0 };
		qaws_polynomial_desc pd;
		co[5] = (qaws_scalar)(-1.15 + 0.1 * k);
		memset(&pd, 0, sizeof(pd));
		pd.dimension = QAWS_DIMENSION_2D;
		pd.degree = 2;
		pd.coefficients = co;
		pd.coefficient_count = 3;
		pd.t_min = -1.5;
		pd.t_max = 1.5;
		ring[k] = NULL;
		qaws_curve_create_polynomial(&pd, &ring[k]);
	}
	for (i = 0; i < 16; i++)
		lv[i] = (qaws_scalar)(0.1 + 0.2 * i);
	for (i = 0; i < n; i++)
	{
		pts[3 * i] = (qaws_scalar)(-1.1 + 2.2 * sbt_rand());
		pts[3 * i + 1] = (qaws_scalar)(-1.1 + 2.2 * sbt_rand());
		pts[3 * i + 2] = (qaws_scalar)(-0.4 + 3.6 * sbt_rand());
		org[3 * i] = (qaws_scalar)(-1.3 + 2.6 * sbt_rand());
		org[3 * i + 1] = (qaws_scalar)(-1.3 + 2.6 * sbt_rand());
		org[3 * i + 2] = (qaws_scalar)4.5;
		dir[3 * i] = (qaws_scalar)(0.6 * (sbt_rand() - 0.5));
		dir[3 * i + 1] = (qaws_scalar)(0.6 * (sbt_rand() - 0.5));
		dir[3 * i + 2] = -1;
		pts2[2 * i] = pts[3 * i];
		pts2[2 * i + 1] = pts[3 * i + 1];
	}
	for (e = 0; e < 3; e++)
	{
		qaws_surface_batch_closest_desc cd;
		qaws_surface_ray_desc rd;
		qaws_closest_desc kd;
		qaws_level_crossing_desc ld;
		qaws_batch_executor const* x = e == 0 ? NULL : &ex[e - 1];
		double t0 = sbt_now();
		sc[e] = (qaws_surface_batch_closest*)malloc(n * sizeof(qaws_surface_batch_closest));
		rh[e] = (qaws_surface_ray_hit*)malloc(n * sizeof(qaws_surface_ray_hit));
		cc[e] = (qaws_closest_point*)malloc(n * sizeof(qaws_closest_point));
		lc[e] = (qaws_level_crossing*)malloc(4096 * sizeof(qaws_level_crossing));
		memset(&cd, 0, sizeof(cd));
		cd.surfaces = (qaws_surface const* const*)sf;
		cd.surface_count = SBT_SURF + 1;
		cd.points = pts;
		cd.point_count = n;
		cd.executor = x;
		qaws_surface_batch_find_closest(&cd, sc[e], NULL);
		memset(&rd, 0, sizeof(rd));
		rd.surfaces = (qaws_surface const* const*)sf;
		rd.surface_count = SBT_SURF + 1;
		rd.origins = org;
		rd.directions = dir;
		rd.ray_count = n;
		rd.executor = x;
		qaws_surface_batch_raycast(&rd, rh[e], NULL);
		memset(&kd, 0, sizeof(kd));
		kd.curves = (qaws_curve const* const*)ring;
		kd.curve_count = 24;
		kd.points = pts2;
		kd.point_count = n;
		kd.executor = x;
		qaws_curve_batch_find_closest(&kd, cc[e], NULL);
		memset(&ld, 0, sizeof(ld));
		ld.curves = (qaws_curve const* const*)ring;
		ld.curve_count = 24;
		ld.field = sbt_field;
		ld.levels = lv;
		ld.level_count = 16;
		ld.executor = x;
		qaws_curve_batch_find_level_crossings(&ld, lc[e], 4096, &nlc[e]);
		if (e == 0) tser = sbt_now() - t0;
		if (e == 2) tthr = sbt_now() - t0;
		if (e > 0)
		{
			int ok = nlc[e] == nlc[0];
			for (i = 0; i < n && ok; i++)
				ok = sc[e][i].surface == sc[0][i].surface && sc[e][i].distance == sc[0][i].distance && rh[e][i].surface == rh[0][i].surface
					&& rh[e][i].t == rh[0][i].t && cc[e][i].curve == cc[0][i].curve && cc[e][i].distance == cc[0][i].distance;
			for (i = 0; i < nlc[e] && i < 4096 && ok; i++)
				ok = lc[e][i].curve == lc[0][i].curve && lc[e][i].level == lc[0][i].level && lc[e][i].parameter == lc[0][i].parameter;
			same[e - 1] = ok;
		}
	}
	printf("    executors: %u points + %u rays + %u curve points + %u level crossings; serial %.3f s, %u threads %.3f s (x%.1f); shuffled order identical %d, threads identical %d\n",
		n, n, n, nlc[0], tser, SBT_THREADS, tthr, tser / (tthr > 1e-6 ? tthr : 1e-6), same[0], same[1]);
	sprintf(msg, "a scrambled executor and %u threads give the serial results bit for bit", SBT_THREADS);
	TEST_ASSERT(same[0] && same[1], msg);
	for (e = 0; e < 3; e++)
	{
		free(sc[e]);
		free(rh[e]);
		free(cc[e]);
		free(lc[e]);
	}
	for (k = 0; k <= SBT_SURF; k++)
		qaws_surface_destroy(sf[k]);
	for (k = 0; k < 24; k++)
		qaws_curve_destroy(ring[k]);
	free(pts);
	free(org);
	free(dir);
	free(pts2);
}

int test_81_surface_batch_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 81: Batched curve / surface intersection\n");
	test_lines();
	test_cubics();
	test_exact();
	test_ssi();
	test_exact_ssi();
	test_surface_sets();
	test_surface_closest();
	test_exact_surface_closest();
	test_raycast();
	test_exact_raycast();
	test_executors();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
