/*
 * Test 104: sampling and meshing regions on surfaces in their metric
 * (qaws_surface_region_sample, qaws_surface_region_triangulate)
 *
 *   - a sphere: Sobol points over a geodesic circle across the seam give
 *     the mean of cos(angle from the centre), (1 + cos(r/R)) / 2; over a
 *     band between latitudes, z is uniform (Archimedes): its mean and mean
 *     square; with a density 1 + z / R the mean of z moves to its closed
 *     form; blue noise keeps its points apart (in space) by more than the
 *     packing radius for that many points over the area
 *   - meshes: a geodesic circle on the sphere (the triangles' areas in
 *     space against the region's area, the longest edge in space against
 *     the target), a parallelogram on a flat sheared patch (its area
 *     exactly), and a circle across both seams of a torus
 */

#include "test_common.h"
#include "qaws_surface_region.h"
#include "qaws_region_sample.h"
#include "qaws_triangulate.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define SST_PI 3.14159265358979323846
#define SST_F32 (sizeof(qaws_scalar) == 4)

static qaws_vec2 sst_uv(double u, double v)
{
	qaws_vec2 p;
	p.x = (qaws_scalar)u; p.y = (qaws_scalar)v;
	return p;
}

static qaws_surface* sst_revolve(double cx, double r, double a0, double a1, qaws_curve** prof)
{
	qaws_arc_segment g;
	qaws_arc_desc ad;
	qaws_surface_revolution_desc rd;
	qaws_surface* s = NULL;
	memset(&g, 0, sizeof(g));
	g.center[0] = (qaws_scalar)cx; g.radius = (qaws_scalar)r; g.angle_start = (qaws_scalar)a0; g.angle_end = (qaws_scalar)a1;
	memset(&ad, 0, sizeof(ad));
	ad.dimension = QAWS_DIMENSION_2D; ad.segments = &g; ad.segment_count = 1;
	qaws_curve_create_arc(prof, &ad);
	memset(&rd, 0, sizeof(rd));
	rd.profile = *prof; rd.axis_direction.z = 1;
	qaws_surface_create_revolution(&s, &rd);
	return s;
}

static qaws_surface_region* sst_region(qaws_surface const* s, qaws_surface_curve* const* cv, unsigned int n, unsigned int paths)
{
	qaws_surface_curve const* pc[8];
	qaws_surface_path P[8];
	qaws_surface_region* r = NULL;
	unsigned int i;
	for (i = 0; i < n; i++) pc[i] = cv[i];
	if (paths == 1)
	{
		P[0].curves = pc; P[0].curve_count = n; P[0].closed = 1;
	}
	else
		for (i = 0; i < paths; i++) { P[i].curves = pc + i; P[i].curve_count = 1; P[i].closed = 1; }
	qaws_surface_clip(&r, s, P, NULL, paths, 0, QAWS_CLIP_NONE, QAWS_FILL_NON_ZERO);
	return r;
}

/* the mesh in space: its area and longest edge */
static void sst_mesh_measure(qaws_surface const* s, qaws_mesh_2d const* m, double* area, double* longest)
{
	unsigned int t, k;
	*area = 0; *longest = 0;
	for (t = 0; t < qaws_mesh_2d_get_triangle_count(m); t++)
	{
		unsigned int v[3];
		qaws_vec3 q[3];
		double ax, ay, az, bx, by, bz, nx, ny, nz;
		qaws_mesh_2d_get_triangle(v, m, t);
		for (k = 0; k < 3; k++)
		{
			qaws_vec2 w = qaws_mesh_2d_get_vertex(m, v[k]);
			qaws_range ru = qaws_surface_get_u_range(s), rv = qaws_surface_get_v_range(s);
			qaws_surface_eval_result e;
			double uu = w.x, vv = w.y;
			/* (a closed surface's mesh vertices may sit a rounding past its domain) */
			if (uu < ru.min_value) uu = ru.min_value;
			if (uu > ru.max_value) uu = ru.max_value;
			if (vv < rv.min_value) vv = rv.min_value;
			if (vv > rv.max_value) vv = rv.max_value;
			qaws_surface_evaluate(&e, s, (qaws_scalar)uu, (qaws_scalar)vv, QAWS_SURFACE_EVAL_POSITION);
			q[k] = e.position;
		}
		ax = q[1].x - q[0].x; ay = q[1].y - q[0].y; az = q[1].z - q[0].z;
		bx = q[2].x - q[0].x; by = q[2].y - q[0].y; bz = q[2].z - q[0].z;
		nx = ay * bz - az * by; ny = az * bx - ax * bz; nz = ax * by - ay * bx;
		*area += 0.5 * sqrt(nx * nx + ny * ny + nz * nz);
		for (k = 0; k < 3; k++)
		{
			qaws_vec3 a = q[k], b = q[(k + 1) % 3];
			double d = sqrt((double)(a.x - b.x) * (a.x - b.x) + (double)(a.y - b.y) * (a.y - b.y) + (double)(a.z - b.z) * (a.z - b.z));
			if (d > *longest) *longest = d;
		}
	}
}

static qaws_scalar sst_density(void* user, qaws_scalar u, qaws_scalar v)
{
	(void)user; (void)u;
	/* 1 + z / R on the sphere: z / R = sin(latitude), latitude = pi v - pi / 2 */
	return (qaws_scalar)(1 + sin(SST_PI * (double)v - SST_PI / 2));
}

static void test_sampling(void)
{
	double R = 1.3, rho = 0.6;
	qaws_curve* prof = NULL;
	qaws_surface* s = sst_revolve(0, R, -SST_PI / 2, SST_PI / 2, &prof);
	qaws_surface_curve *c = NULL, *lat[2] = { 0 };
	qaws_surface_region *rc, *rb;
	qaws_surface_sample_desc sd;
	unsigned int n = 4096, i;
	qaws_vec3* p = (qaws_vec3*)malloc(sizeof(qaws_vec3) * n);
	qaws_vec2* w = (qaws_vec2*)malloc(sizeof(qaws_vec2) * n);
	qaws_curve *k1 = NULL, *k2 = NULL;
	double cx[3], tol = 2e-3;
	qaws_status s1, s2, s3, s4;
	char msg[300];
	qaws_surface_curve_create_geodesic_arc(&c, s, sst_uv(0.02, 0.55), (qaws_scalar)rho, 0, (qaws_scalar)(2 * SST_PI), 0);
	rc = sst_region(s, &c, 1, 1);
	{
		double la = SST_PI * 0.55 - SST_PI / 2, lo = 2 * SST_PI * 0.02;
		cx[0] = cos(la) * cos(lo); cx[1] = cos(la) * sin(lo); cx[2] = sin(la);
	}
	/* the band between latitude lines drawn in (u, v) */
	{
		qaws_bezier_desc bz;
		qaws_scalar cp[4];
		memset(&bz, 0, sizeof(bz));
		bz.dimension = QAWS_DIMENSION_2D; bz.degree = 1; bz.control_points = cp; bz.control_point_count = 2;
		cp[0] = 0; cp[1] = (qaws_scalar)0.35; cp[2] = 1; cp[3] = (qaws_scalar)0.35;
		qaws_curve_create_bezier(&k1, &bz);
		cp[0] = 1; cp[1] = (qaws_scalar)0.7; cp[2] = 0; cp[3] = (qaws_scalar)0.7;
		qaws_curve_create_bezier(&k2, &bz);
		qaws_surface_curve_create_uv(&lat[0], s, k1);
		qaws_surface_curve_create_uv(&lat[1], s, k2);
		rb = sst_region(s, lat, 2, 2);
	}
	memset(&sd, 0, sizeof(sd));
	sd.kind = QAWS_REGION_SAMPLE_SOBOL; sd.seed = 7;
	s1 = qaws_surface_region_sample(NULL, p, rc, n, &sd);
	{
		double mean = 0;
		for (i = 0; i < n; i++) mean += ((double)p[i].x * cx[0] + (double)p[i].y * cx[1] + (double)p[i].z * cx[2]) / R;
		mean /= n;
		s2 = qaws_surface_region_sample(w, p, rb, n, &sd);
		{
			double z1 = R * sin(SST_PI * 0.35 - SST_PI / 2), z2 = R * sin(SST_PI * 0.7 - SST_PI / 2), mz = 0, mz2 = 0;
			int inside = 1;
			for (i = 0; i < n; i++)
			{
				mz += p[i].z; mz2 += (double)p[i].z * p[i].z;
				if (!(w[i].y >= 0.35 - 1e-6 && w[i].y <= 0.7 + 1e-6)) inside = 0;
			}
			mz /= n; mz2 /= n;
			sprintf(msg, "sphere, Sobol: a circle across the seam, mean cos %.2g off; a band, z uniform (Archimedes): mean %.2g, mean square %.2g off",
				fabs(mean - (1 + cos(rho / R)) / 2), fabs(mz - 0.5 * (z1 + z2)) / R, fabs(mz2 - (z1 * z1 + z1 * z2 + z2 * z2) / 3) / (R * R));
			TEST_ASSERT(!s1 && !s2 && inside && fabs(mean - (1 + cos(rho / R)) / 2) <= tol && fabs(mz - 0.5 * (z1 + z2)) / R <= tol &&
				fabs(mz2 - (z1 * z1 + z1 * z2 + z2 * z2) / 3) / (R * R) <= tol, msg);
			/* the density 1 + z / R: z's mean over the band becomes
			   (int z (1 + z/R) dz) / (int (1 + z/R) dz) */
			sd.density = sst_density;
			s3 = qaws_surface_region_sample(NULL, p, rb, n, &sd);
			sd.density = NULL;
			{
				double num = 0.5 * (z2 * z2 - z1 * z1) + (z2 * z2 * z2 - z1 * z1 * z1) / (3 * R), den = (z2 - z1) + 0.5 * (z2 * z2 - z1 * z1) / R, md = 0;
				for (i = 0; i < n; i++) md += p[i].z;
				md /= n;
				sprintf(msg, "sphere, Sobol with the density 1 + z/R over the band: mean z %.2g off its closed form", fabs(md - num / den) / R);
				TEST_ASSERT(!s3 && fabs(md - num / den) / R <= tol, msg);
			}
		}
	}
	/* blue noise: apart by more than the packing radius sqrt(A / (2 sqrt 3 n)) */
	{
		unsigned int nb = 500, j;
		qaws_scalar a = 0;
		double dmin = 1e30, rmax;
		int inside = 1;
		sd.kind = QAWS_REGION_SAMPLE_BLUE_NOISE;
		s4 = qaws_surface_region_sample(w, p, rc, nb, &sd);
		qaws_surface_region_get_area(&a, rc);
		rmax = sqrt((double)a / (2 * sqrt(3.0) * nb));
		for (i = 0; i < nb; i++)
		{
			int in = 0;
			qaws_surface_region_contains(&in, rc, w[i]);
			if (in == 0) inside = 0;
			for (j = i + 1; j < nb; j++)
			{
				double d = sqrt(((double)p[i].x - p[j].x) * ((double)p[i].x - p[j].x) + ((double)p[i].y - p[j].y) * ((double)p[i].y - p[j].y) +
					((double)p[i].z - p[j].z) * ((double)p[i].z - p[j].z));
				if (d < dmin) dmin = d;
			}
		}
		sprintf(msg, "sphere, blue noise over a circle across the seam: 500 points apart by %.4f, %.2f of the packing radius", dmin, dmin / rmax);
		TEST_ASSERT(!s4 && inside && dmin >= 0.75 * rmax, msg);
	}
	free(p); free(w);
	qaws_surface_region_destroy(rc); qaws_surface_region_destroy(rb);
	qaws_surface_curve_destroy(c); qaws_surface_curve_destroy(lat[0]); qaws_surface_curve_destroy(lat[1]);
	qaws_curve_destroy(k1); qaws_curve_destroy(k2);
	qaws_surface_destroy(s);
	qaws_curve_destroy(prof);
}

static void test_meshing(void)
{
	char msg[300];
	/* the sphere: a geodesic circle across the seam */
	{
		double R = 1.3, L = 0.05, am = 0, em = 0;
		qaws_curve* prof = NULL;
		qaws_surface* s = sst_revolve(0, R, -SST_PI / 2, SST_PI / 2, &prof);
		qaws_surface_curve* c = NULL;
		qaws_surface_region* r;
		qaws_surface_mesh_desc md;
		qaws_mesh_2d* m = NULL;
		qaws_scalar a = 0;
		qaws_status st;
		qaws_surface_curve_create_geodesic_arc(&c, s, sst_uv(0.02, 0.55), (qaws_scalar)0.6, 0, (qaws_scalar)(2 * SST_PI), 0);
		r = sst_region(s, &c, 1, 1);
		qaws_surface_region_get_area(&a, r);
		memset(&md, 0, sizeof(md));
		md.edge_length = (qaws_scalar)L;
		st = qaws_surface_region_triangulate(&m, r, &md);
		if (st == QAWS_STATUS_OK) sst_mesh_measure(s, m, &am, &em);
		sprintf(msg, "sphere: a circle across the seam meshed to 0.05: %u triangles, their area %.2g off (relative), longest edge %.4f",
			m ? qaws_mesh_2d_get_triangle_count(m) : 0, st ? 1.0 : fabs(am - a) / a, st ? 0.0 : em);
		TEST_ASSERT(!st && fabs(am - a) / a <= 5e-4 && em <= 1.05 * L, msg);
		qaws_mesh_2d_destroy(m);
		qaws_surface_region_destroy(r);
		qaws_surface_curve_destroy(c);
		qaws_surface_destroy(s);
		qaws_curve_destroy(prof);
	}
	/* a flat sheared patch: a geodesic parallelogram's area exactly */
	{
		qaws_surface_bilinear_desc bd;
		qaws_surface* s = NULL;
		qaws_surface_curve* sq[4] = { 0 };
		qaws_vec2 corner[4];
		qaws_surface_region* r;
		qaws_surface_mesh_desc md;
		qaws_mesh_2d* m = NULL;
		double am = 0, em = 0;
		unsigned int i;
		qaws_status st;
		memset(&bd, 0, sizeof(bd));
		bd.p10.x = 2; bd.p01.x = (qaws_scalar)0.6; bd.p01.y = 3; bd.p11.x = (qaws_scalar)2.6; bd.p11.y = 3;
		qaws_surface_create_bilinear(&s, &bd);
		corner[0] = sst_uv(0.3, 0.3); corner[1] = sst_uv(0.62, 0.3); corner[2] = sst_uv(0.62, 0.58); corner[3] = sst_uv(0.3, 0.58);
		for (i = 0; i < 4; i++) qaws_surface_curve_create_geodesic(&sq[i], s, corner[i], corner[(i + 1) % 4], 0);
		r = sst_region(s, sq, 4, 1);
		memset(&md, 0, sizeof(md));
		md.edge_length = (qaws_scalar)0.05;
		st = qaws_surface_region_triangulate(&m, r, &md);
		if (st == QAWS_STATUS_OK) sst_mesh_measure(s, m, &am, &em);
		sprintf(msg, "flat sheared patch: a parallelogram's mesh, its area against 0.64 x 0.84 %.2g, longest edge %.4f", fabs(am - 0.64 * 0.84), em);
		TEST_ASSERT(!st && fabs(am - 0.64 * 0.84) <= (SST_F32 ? 1e-5 : 1e-9) && em <= 0.0525, msg);
		qaws_mesh_2d_destroy(m);
		qaws_surface_region_destroy(r);
		for (i = 0; i < 4; i++) qaws_surface_curve_destroy(sq[i]);
		qaws_surface_destroy(s);
	}
	/* a torus: a circle across both seams */
	{
		qaws_curve* prof = NULL;
		qaws_surface* s = sst_revolve(1.0, 0.42, 0, 2 * SST_PI, &prof);
		qaws_surface_curve* c = NULL;
		qaws_surface_region* r;
		qaws_surface_mesh_desc md;
		qaws_mesh_2d* m = NULL;
		qaws_scalar a = 0;
		double am = 0, em = 0;
		qaws_status st;
		qaws_surface_curve_create_geodesic_arc(&c, s, sst_uv(0.0, 0.0), (qaws_scalar)0.3, 0, (qaws_scalar)(2 * SST_PI), 0);
		r = sst_region(s, &c, 1, 1);
		qaws_surface_region_get_area(&a, r);
		memset(&md, 0, sizeof(md));
		md.edge_length = (qaws_scalar)0.03;
		st = qaws_surface_region_triangulate(&m, r, &md);
		if (st == QAWS_STATUS_OK) sst_mesh_measure(s, m, &am, &em);
		sprintf(msg, "torus: a circle across both seams meshed to 0.03 (four pieces in (u, v)): area %.2g off (relative), longest edge %.4f",
			st ? 1.0 : fabs(am - a) / a, em);
		TEST_ASSERT(!st && fabs(am - a) / a <= 5e-4 && em <= 1.05 * 0.03, msg);
		qaws_mesh_2d_destroy(m);
		qaws_surface_region_destroy(r);
		qaws_surface_curve_destroy(c);
		qaws_surface_destroy(s);
		qaws_curve_destroy(prof);
	}
}

int test_104_surface_sample_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 104: sampling and meshing regions on surfaces\n");
	test_sampling();
	test_meshing();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
