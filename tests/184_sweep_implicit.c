/*
 * Test 184: implicit sweeps (qaws_mesh_3d_sweep_implicit), small meshes, coarse grids
 *
 *   - a unit cube translated by (1, 1, 0): closed, its volume within 2% of the Minkowski sum's 3, every vertex of either
 *     sweep within a cell of the other (the explicit one exact here);
 *   - a box turned a quarter about z: the implicit and the explicit sweeps within a cell of each other both ways, their
 *     volumes within 2%;
 *   - a vertex motion (a box's top rising from 1 to 2): closed, volume 2 within 3%;
 *   - refusals: no mesh, a vertex motion of another count;
 *   - an attribute track: a vertex "heat" animated from 0 to 1 equals the generated time at every corner;
 *   - passes: a ball out along x and back (x = 3 t (1 - t), even about t = 0.5): most corners passed twice, at t and 1 - t
 *     (a ball's sweep is touched at an instant, unlike a face sliding along itself, touched over a span), so the heat track
 *     combined by SUM is 1 there (within the pass times' refinement), the others once; the same swept by slabs (tolerance
 *     0.01: the slabs know only the first touch, the passes found from the motion within twice the tolerance);
 *   - generated attributes (qaws_sweep_desc.generate): the cube moved by (1, 1, 0) at a constant rate: velocity (1, 1, 0)
 *     at every corner, distance sqrt 2 times the time, orientation the identity, and each corner its source point (its
 *     source_face and source_uv on the cube) moved by (time, time, 0), within a cell;
 *   - attributes: the cube translated with its positions as a vertex attribute (each corner's value the source point: its
 *     z and x - y those of the corner, the translation along (1, 1, 0), within half a cell), each face's corners its side's
 *     number (a corner attribute: no value blended across the cube's edges, nine corners in ten their triangle's own side,
 *     the rest on the walls the cube's vertical edges sweep, either side's), a face attribute the same; a box turned a
 *     quarter with a constant VECTOR (1, 0, 0): unit, flat, turned between 0 and a quarter;
 *   - per texel (qaws_mesh_3d_sweep_texture): a unit cube slowing to a stop along x, its swept solid's 12-triangle box
 *     flattened and its "time" on 128 x 128 texels: on the top face within 0.01 of 1 - (2 - x)^(1/3), where the
 *     vertices' exact times interpolated are off by over 0.5.
 */

#include "test_common.h"
#include "qaws_mesh_sweep.h"
#include "qaws_mesh_bvh.h"
#include "qaws_mesh_param.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

static double si_volume(qaws_mesh_3d const* m)
{
	unsigned int nt = m ? qaws_mesh_3d_get_triangle_count(m) : 0, i;
	double v = 0;
	for (i = 0; i < nt; i++)
	{
		unsigned int c[3];
		qaws_vec3 a, b, e;
		qaws_mesh_3d_get_triangle(c, m, i);
		a = qaws_mesh_3d_get_vertex(m, c[0]);
		b = qaws_mesh_3d_get_vertex(m, c[1]);
		e = qaws_mesh_3d_get_vertex(m, c[2]);
		v += ((double)a.x * ((double)b.y * e.z - (double)b.z * e.y) + (double)a.y * ((double)b.z * e.x - (double)b.x * e.z) +
			(double)a.z * ((double)b.x * e.y - (double)b.y * e.x)) / 6;
	}
	return v;
}

/* the farthest vertex of a from b's surface */
static double si_far(qaws_mesh_3d const* a, qaws_mesh_3d const* b)
{
	qaws_mesh_bvh* t = NULL;
	unsigned int n = a && b ? qaws_mesh_3d_get_vertex_count(a) : 0, i;
	double mx = 0;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (3 * (size_t)n + 1));
	qaws_mesh_hit* h = (qaws_mesh_hit*)malloc(sizeof(qaws_mesh_hit) * ((size_t)n + 1));
	if (!n || qaws_mesh_bvh_create_mesh_3d(&t, b, NULL) != QAWS_STATUS_OK)
		mx = 1e30;
	else
	{
		qaws_mesh_3d_get_positions(x, a);
		qaws_mesh_bvh_closest(h, t, x, n, 0, NULL);
		for (i = 0; i < n; i++)
			mx = fmax(mx, (double)h[i].distance);
	}
	qaws_mesh_bvh_destroy(t);
	free(x); free(h);
	return mx;
}

static qaws_mesh_3d* si_box(double x0, double x1, double y0, double y1, double z0, double z1)
{
	static unsigned int const f[36] = { 0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 2, 3, 7, 2, 7, 6, 1, 2, 6, 1, 6, 5, 0, 4, 7, 0, 7, 3 };
	qaws_scalar x[24];
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	unsigned int i;
	for (i = 0; i < 8; i++)
	{
		x[3 * i] = (qaws_scalar)((i == 1 || i == 2 || i == 5 || i == 6) ? x1 : x0);
		x[3 * i + 1] = (qaws_scalar)((i == 2 || i == 3 || i == 6 || i == 7) ? y1 : y0);
		x[3 * i + 2] = (qaws_scalar)(i >= 4 ? z1 : z0);
	}
	memset(&d, 0, sizeof d);
	d.positions = x; d.vertex_count = 8; d.face_vertices = f; d.face_count = 12;
	qaws_mesh_3d_create(&m, &d);
	return m;
}

/* the explicit and the implicit sweeps of b by m: volumes, each one's farthest vertex from the other, the cell */
static void si_both(qaws_mesh_3d const* b, qaws_motion const* m, double tol, unsigned int res, double* ve, double* vi, double* ie, double* ei, double* cell,
	int* closed)
{
	qaws_mesh_3d *e = NULL, *i = NULL;
	qaws_sweep_desc d;
	qaws_sweep_report r;
	memset(&d, 0, sizeof d);
	d.tolerance = (qaws_scalar)tol;
	d.resolution = res;
	qaws_mesh_3d_sweep(&e, NULL, b, m, &d);
	qaws_mesh_3d_sweep_implicit(&i, &r, b, m, &d);
	*ve = si_volume(e);
	*vi = si_volume(i);
	*ie = si_far(i, e);
	*ei = si_far(e, i);
	*closed = r.closed;
	*cell = 1;
	if (e)
	{
		/* (the explicit sweep's box: about the grid's) */
		unsigned int n = qaws_mesh_3d_get_vertex_count(e), j, k;
		double lo[3] = { 1e30, 1e30, 1e30 }, hi[3] = { -1e30, -1e30, -1e30 }, ext = 0;
		for (j = 0; j < n; j++)
		{
			qaws_vec3 p = qaws_mesh_3d_get_vertex(e, j);
			double q[3] = { p.x, p.y, p.z };
			for (k = 0; k < 3; k++)
			{
				lo[k] = fmin(lo[k], q[k]);
				hi[k] = fmax(hi[k], q[k]);
			}
		}
		for (k = 0; k < 3; k++)
			ext = fmax(ext, hi[k] - lo[k]);
		*cell = ext / res;
	}
	qaws_mesh_3d_destroy(e);
	qaws_mesh_3d_destroy(i);
}

static void test_rigid(void)
{
	qaws_mesh_3d *c = si_box(0, 1, 0, 1, 0, 1), *b = si_box(1, 2, -0.1, 0.1, 0, 0.6);
	qaws_motion *mt = NULL, *mr = NULL;
	qaws_motion_pose s, e;
	double ve, vi, ie, ei, cell, ve2, vi2, ie2, ei2, cell2;
	int cl, cl2;
	char msg[400];
	memset(&s, 0, sizeof s);
	memset(&e, 0, sizeof e);
	s.rotation[0] = e.rotation[0] = 1;
	e.position.x = e.position.y = 1;
	s.velocity = e.velocity = e.position;
	qaws_motion_create_poses(&mt, &s, &e, 0);
	si_both(c, mt, 1e-3, 32, &ve, &vi, &ie, &ei, &cell, &cl);
	memset(&s, 0, sizeof s);
	memset(&e, 0, sizeof e);
	s.rotation[0] = 1;
	e.rotation[0] = (qaws_scalar)cos(3.14159265358979323846 / 4);
	e.rotation[3] = (qaws_scalar)sin(3.14159265358979323846 / 4);
	s.angular_velocity.z = e.angular_velocity.z = (qaws_scalar)(3.14159265358979323846 / 2);
	qaws_motion_create_poses(&mr, &s, &e, 0);
	si_both(b, mr, 1e-2, 32, &ve2, &vi2, &ie2, &ei2, &cell2, &cl2);
	sprintf(msg, "a cube translated: volume %.4f (3), apart %.4f and %.4f (cell %.4f), closed %d; a box turned a quarter: volumes %.4f and %.4f, apart %.4f and %.4f (cell %.4f), closed %d",
		vi, ie, ei, cell, cl, ve2, vi2, ie2, ei2, cell2, cl2);
	printf("  %s\n", msg);
	TEST_ASSERT(cl && fabs(vi - 3) < 0.06 && ie < cell && ei < cell && cl2 && fabs(vi2 - ve2) < 0.02 * ve2 && ie2 < cell2 && ei2 < cell2, msg);
	qaws_motion_destroy(mt); qaws_motion_destroy(mr);
	qaws_mesh_3d_destroy(c); qaws_mesh_3d_destroy(b);
}

static void test_vertex_motion(void)
{
	qaws_mesh_3d *b = si_box(0, 1, 0, 1, 0, 1), *o = NULL;
	qaws_motion *mv = NULL, *m7 = NULL;
	qaws_scalar p0[24], p1[24], v[24];
	qaws_sweep_desc d;
	qaws_sweep_report r;
	qaws_status s1, s2;
	unsigned int i;
	double vo;
	char msg[300];
	qaws_mesh_3d_get_positions(p0, b);
	memcpy(p1, p0, sizeof p1);
	for (i = 4; i < 8; i++)
		p1[3 * i + 2] = 2;
	for (i = 0; i < 24; i++)
		v[i] = p1[i] - p0[i];
	qaws_motion_create_vertices(&mv, p0, p1, v, v, NULL, NULL, 8, 0);
	qaws_motion_create_vertices(&m7, p0, p1, v, v, NULL, NULL, 7, 0);
	memset(&d, 0, sizeof d);
	d.resolution = 24;
	qaws_mesh_3d_sweep_implicit(&o, &r, b, mv, &d);
	vo = si_volume(o);
	qaws_mesh_3d_destroy(o);
	o = NULL;
	s1 = qaws_mesh_3d_sweep_implicit(&o, NULL, NULL, mv, &d);
	s2 = qaws_mesh_3d_sweep_implicit(&o, NULL, b, m7, &d);
	sprintf(msg, "a box's top rising: volume %.4f (2), %u samples, closed %d; refusals %d %d", vo, r.slabs, r.closed, (int)s1, (int)s2);
	printf("  %s\n", msg);
	TEST_ASSERT(r.closed && fabs(vo - 2) < 0.06 && s1 == QAWS_STATUS_INVALID_ARGUMENT && s2 == QAWS_STATUS_INVALID_ARGUMENT && !o, msg);
	qaws_motion_destroy(mv); qaws_motion_destroy(m7);
	qaws_mesh_3d_destroy(b);
}

static qaws_scalar* si_values(qaws_mesh_3d const* m, char const* name, unsigned int n)
{
	unsigned int id;
	qaws_scalar* v = (qaws_scalar*)malloc(sizeof(qaws_scalar) * n * (3 * (size_t)qaws_mesh_3d_get_triangle_count(m) + 1));
	if (v && qaws_mesh_3d_find_attribute(&id, m, name) == QAWS_STATUS_OK)
		qaws_mesh_3d_get_attribute(v, m, id);
	else
	{
		free(v);
		v = NULL;
	}
	return v;
}

static void test_attributes(void)
{
	qaws_mesh_3d *c = si_box(0, 1, 0, 1, 0, 1), *b = si_box(1, 2, -0.1, 0.1, 0, 0.6), *o = NULL, *ot = NULL;
	qaws_motion *mt = NULL, *mr = NULL;
	qaws_motion_pose s, e;
	qaws_sweep_desc d;
	qaws_attribute_desc ad;
	qaws_scalar pos[24], side[36], fside[12], dir[24], *pv = NULL, *cv = NULL, *fv = NULL, *dv = NULL;
	unsigned int i, k, nt = 0, ntt = 0, id, blended = 0, own = 0;
	double perr = 0, cell = 3.0 / 32, worst = 0, amin = 9, amax = -9;
	char msg[400];
	qaws_mesh_3d_get_positions(pos, c);
	for (i = 0; i < 36; i++)
		side[i] = (qaws_scalar)(i / 6);
	for (i = 0; i < 12; i++)
		fside[i] = (qaws_scalar)(i / 2);
	for (i = 0; i < 24; i++)
		dir[i] = (qaws_scalar)(i % 3 == 0);
	ad.name = "pos"; ad.element = QAWS_MESH_VERTEX; ad.components = 3; ad.rule = QAWS_ATTRIBUTE_LINEAR;
	qaws_mesh_3d_add_attribute(&id, c, pos, &ad);
	ad.name = "side"; ad.element = QAWS_MESH_CORNER; ad.components = 1;
	qaws_mesh_3d_add_attribute(&id, c, side, &ad);
	ad.name = "fside"; ad.element = QAWS_MESH_FACE;
	qaws_mesh_3d_add_attribute(&id, c, fside, &ad);
	ad.name = "dir"; ad.element = QAWS_MESH_VERTEX; ad.components = 3; ad.rule = QAWS_ATTRIBUTE_VECTOR;
	qaws_mesh_3d_add_attribute(&id, b, dir, &ad);
	memset(&s, 0, sizeof s);
	memset(&e, 0, sizeof e);
	s.rotation[0] = e.rotation[0] = 1;
	e.position.x = e.position.y = 1;
	s.velocity = e.velocity = e.position;
	qaws_motion_create_poses(&mt, &s, &e, 0);
	memset(&s, 0, sizeof s);
	memset(&e, 0, sizeof e);
	s.rotation[0] = 1;
	e.rotation[0] = (qaws_scalar)cos(3.14159265358979323846 / 4);
	e.rotation[3] = (qaws_scalar)sin(3.14159265358979323846 / 4);
	s.angular_velocity.z = e.angular_velocity.z = (qaws_scalar)(3.14159265358979323846 / 2);
	qaws_motion_create_poses(&mr, &s, &e, 0);
	memset(&d, 0, sizeof d);
	d.resolution = 32;
	qaws_mesh_3d_sweep_implicit(&o, NULL, c, mt, &d);
	d.resolution = 24;
	qaws_mesh_3d_sweep_implicit(&ot, NULL, b, mr, &d);
	if (o)
	{
		nt = qaws_mesh_3d_get_triangle_count(o);
		pv = si_values(o, "pos", 3);
		cv = si_values(o, "side", 1);
		fv = si_values(o, "fside", 1);
	}
	for (i = 0; pv && cv && fv && i < nt; i++)
	{
		unsigned int t[3];
		qaws_mesh_3d_get_triangle(t, o, i);
		for (k = 0; k < 3; k++)
		{
			qaws_vec3 x = qaws_mesh_3d_get_vertex(o, t[k]);
			double z = (double)cv[3 * i + k];
			perr = fmax(perr, fabs((double)pv[3 * (3 * i + k) + 2] - x.z));
			perr = fmax(perr, fabs(((double)pv[3 * (3 * i + k)] - pv[3 * (3 * i + k) + 1]) - ((double)x.x - x.y)));
			if (fabs(z - floor(z + 0.5)) > 1e-6)
				blended++;
			else if (fabs(z - (double)fv[i]) < 1e-6)
				own++;
		}
	}
	if (ot)
	{
		ntt = qaws_mesh_3d_get_triangle_count(ot);
		dv = si_values(ot, "dir", 3);
	}
	for (i = 0; dv && i < 3 * ntt; i++)
	{
		double a = atan2((double)dv[3 * i + 1], (double)dv[3 * i]);
		worst = fmax(worst, fabs(sqrt((double)dv[3 * i] * dv[3 * i] + (double)dv[3 * i + 1] * dv[3 * i + 1] + (double)dv[3 * i + 2] * dv[3 * i + 2]) - 1) + fabs((double)dv[3 * i + 2]));
		amin = fmin(amin, a);
		amax = fmax(amax, a);
	}
	sprintf(msg, "attributes: positions off by %.4f (cell %.4f), corners blended %u, their triangle's side %u of %u; a turned VECTOR unit and flat within %.1e, turned %.4f to %.4f",
		perr, cell, blended, own, 3 * nt, worst, amin, amax);
	printf("  %s\n", msg);
	TEST_ASSERT(pv && cv && fv && dv && perr < 0.5 * cell && !blended && own > 0.85 * 3 * nt && worst < 1e-9 && amin > -1e-6 && amax < 1.5708 + 1e-6 && amax > 1.5, msg);
	free(pv); free(cv); free(fv); free(dv);
	qaws_motion_destroy(mt); qaws_motion_destroy(mr);
	qaws_mesh_3d_destroy(c); qaws_mesh_3d_destroy(b); qaws_mesh_3d_destroy(o); qaws_mesh_3d_destroy(ot);
}

/* a result corner's source point on the body (rest coordinates): face source_face's triangle, weights source_uv */
static void si_source_point(qaws_mesh_3d const* body, unsigned int face, double const* uv, double* p)
{
	unsigned int t[3], k;
	qaws_vec3 a, b, c;
	qaws_mesh_3d_get_triangle(t, body, face);
	a = qaws_mesh_3d_get_vertex(body, t[0]);
	b = qaws_mesh_3d_get_vertex(body, t[1]);
	c = qaws_mesh_3d_get_vertex(body, t[2]);
	for (k = 0; k < 3; k++)
	{
		double pa = k == 0 ? a.x : k == 1 ? a.y : a.z, pb = k == 0 ? b.x : k == 1 ? b.y : b.z, pc = k == 0 ? c.x : k == 1 ? c.y : c.z;
		p[k] = (1 - uv[0] - uv[1]) * pa + uv[0] * pb + uv[1] * pc;
	}
}

/* attribute name's values per corner of o (n each), or NULL */
static qaws_scalar* si_corner_values(qaws_mesh_3d const* o, char const* name, unsigned int n)
{
	unsigned int id;
	qaws_scalar* v;
	if (qaws_mesh_3d_find_attribute(&id, o, name) != QAWS_STATUS_OK)
		return NULL;
	v = (qaws_scalar*)malloc(sizeof(qaws_scalar) * n * (qaws_mesh_3d_get_corner_count(o) + 1));
	qaws_mesh_3d_get_attribute(v, o, id);
	return v;
}

static void test_generate(void)
{
	qaws_mesh_3d *c = si_box(0, 1, 0, 1, 0, 1), *o = NULL;
	qaws_motion* m = NULL;
	qaws_motion_pose s, e;
	qaws_sweep_desc d;
	qaws_scalar *tv = NULL, *vv = NULL, *dv = NULL, *qv = NULL, *fv = NULL, *uv = NULL;
	unsigned int i, k, nt = 0;
	double wv = 0, wd = 0, wq = 0, wp = 0, cell = 3.0 / 24;
	char msg[300];
	memset(&s, 0, sizeof s);
	memset(&e, 0, sizeof e);
	s.rotation[0] = e.rotation[0] = 1;
	e.position.x = e.position.y = 1;
	s.velocity = e.velocity = e.position;
	qaws_motion_create_poses(&m, &s, &e, 0);
	memset(&d, 0, sizeof d);
	d.resolution = 24;
	d.generate = QAWS_SWEEP_GENERATE_TIME | QAWS_SWEEP_GENERATE_VELOCITY | QAWS_SWEEP_GENERATE_DISTANCE | QAWS_SWEEP_GENERATE_ORIENTATION | QAWS_SWEEP_GENERATE_SOURCE;
	qaws_mesh_3d_sweep_implicit(&o, NULL, c, m, &d);
	if (o)
	{
		nt = qaws_mesh_3d_get_triangle_count(o);
		tv = si_corner_values(o, "time", 1); vv = si_corner_values(o, "velocity", 3); dv = si_corner_values(o, "distance", 1);
		qv = si_corner_values(o, "orientation", 4); fv = si_corner_values(o, "source_face", 1); uv = si_corner_values(o, "source_uv", 2);
	}
	for (i = 0; tv && vv && dv && qv && fv && uv && i < 3 * nt; i++)
	{
		unsigned int t3[3];
		double sp[3], u2[2] = { uv[2 * i], uv[2 * i + 1] }, t = (double)tv[i];
		qaws_vec3 p;
		qaws_mesh_3d_get_triangle(t3, o, i / 3);
		p = qaws_mesh_3d_get_vertex(o, t3[i % 3]);
		wv = fmax(wv, fabs((double)vv[3 * i] - 1) + fabs((double)vv[3 * i + 1] - 1) + fabs((double)vv[3 * i + 2]));
		wd = fmax(wd, fabs((double)dv[i] - sqrt(2.0) * t));
		wq = fmax(wq, fabs((double)qv[4 * i] - 1) + fabs((double)qv[4 * i + 1]) + fabs((double)qv[4 * i + 2]) + fabs((double)qv[4 * i + 3]));
		si_source_point(c, (unsigned int)lround((double)fv[i]), u2, sp);
		for (k = 0; k < 3; k++)
			sp[k] += k < 2 ? t : 0;
		wp = fmax(wp, sqrt(((double)p.x - sp[0]) * ((double)p.x - sp[0]) + ((double)p.y - sp[1]) * ((double)p.y - sp[1]) + ((double)p.z - sp[2]) * ((double)p.z - sp[2])));
	}
	sprintf(msg, "generated: velocity off (1, 1, 0) by %.1e, distance off sqrt 2 t by %.1e, orientation off the identity by %.1e, corners off their moved source by %.4f (cell %.4f)",
		wv, wd, wq, wp, cell);
	printf("  %s\n", msg);
	TEST_ASSERT(o && tv && vv && dv && qv && fv && uv && wv < 1e-9 && wd < 1e-3 && wq < 1e-9 && wp < cell, msg);
	free(tv); free(vv); free(dv); free(qv); free(fv); free(uv);
	qaws_motion_destroy(m);
	qaws_mesh_3d_destroy(c); qaws_mesh_3d_destroy(o);
}

/* a sweep of b by m with "heat" (a vertex attribute) animated from 0 at t = 0 to 1 at t = 1, the generated time beside it:
   the largest gap between them over the corners (-1: missing); st the call's status */
static double si_heat_gap(qaws_mesh_3d* b, qaws_motion const* m, int implicit, double tol, qaws_status* st)
{
	unsigned int nv = qaws_mesh_3d_get_vertex_count(b), i, id;
	qaws_scalar* keys = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * nv);
	qaws_scalar times[2] = { 0, 1 }, *zero = (qaws_scalar*)calloc(nv, sizeof(qaws_scalar)), *hv, *tv;
	qaws_attribute_desc ad;
	qaws_attribute_track tr;
	qaws_sweep_desc d;
	qaws_mesh_3d* o = NULL;
	double gap = -1;
	memset(&ad, 0, sizeof ad);
	ad.name = "heat"; ad.element = QAWS_MESH_VERTEX; ad.components = 1; ad.rule = QAWS_ATTRIBUTE_LINEAR;
	if (qaws_mesh_3d_find_attribute(&id, b, "heat") != QAWS_STATUS_OK)
		qaws_mesh_3d_add_attribute(&id, b, zero, &ad);
	for (i = 0; i < nv; i++)
	{
		keys[i] = 0;
		keys[nv + i] = 1;
	}
	tr.name = "heat"; tr.times = times; tr.values = keys; tr.key_count = 2;
	memset(&d, 0, sizeof d);
	d.tolerance = (qaws_scalar)tol;
	d.resolution = 24;
	d.generate = QAWS_SWEEP_GENERATE_TIME;
	d.tracks = &tr;
	d.track_count = 1;
	*st = implicit ? qaws_mesh_3d_sweep_implicit(&o, NULL, b, m, &d) : qaws_mesh_3d_sweep(&o, NULL, b, m, &d);
	if (o && (hv = si_corner_values(o, "heat", 1)) != NULL)
	{
		if ((tv = si_corner_values(o, "time", 1)) != NULL)
		{
			gap = 0;
			for (i = 0; i < qaws_mesh_3d_get_corner_count(o); i++)
				gap = fmax(gap, fabs((double)hv[i] - (double)tv[i]));
			free(tv);
		}
		free(hv);
	}
	qaws_mesh_3d_destroy(o);
	free(keys); free(zero);
	return gap;
}

static void test_track(void)
{
	qaws_mesh_3d* c = si_box(0, 1, 0, 1, 0, 1);
	qaws_motion* m = NULL;
	qaws_motion_pose s, e;
	qaws_status st;
	double g;
	char msg[200];
	memset(&s, 0, sizeof s);
	memset(&e, 0, sizeof e);
	s.rotation[0] = e.rotation[0] = 1;
	e.position.x = e.position.y = 1;
	s.velocity = e.velocity = e.position;
	qaws_motion_create_poses(&m, &s, &e, 0);
	g = si_heat_gap(c, m, 1, 0, &st);
	sprintf(msg, "a track implicitly: heat off the time by %.1e", g);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && g >= 0 && g < 1e-9, msg);
	qaws_motion_destroy(m);
	qaws_mesh_3d_destroy(c);
}

/* a UV ball of radius r at the origin, poles on z */
static qaws_mesh_3d* si_ball(double r, unsigned int nu, unsigned int nv)
{
	unsigned int n = nu * (nv - 1) + 2, i, j, k = 0;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n);
	unsigned int* f = (unsigned int*)malloc(sizeof(unsigned int) * 6 * nu * nv);
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	for (j = 1; j < nv; j++)
		for (i = 0; i < nu; i++)
		{
			double th = 3.14159265358979323846 * j / nv, ph = 2 * 3.14159265358979323846 * i / nu;
			unsigned int q = (j - 1) * nu + i;
			x[3 * q] = (qaws_scalar)(r * sin(th) * cos(ph));
			x[3 * q + 1] = (qaws_scalar)(r * sin(th) * sin(ph));
			x[3 * q + 2] = (qaws_scalar)(r * cos(th));
		}
	x[3 * (n - 2)] = x[3 * (n - 2) + 1] = 0;
	x[3 * (n - 2) + 2] = (qaws_scalar)r;
	x[3 * (n - 1)] = x[3 * (n - 1) + 1] = 0;
	x[3 * (n - 1) + 2] = (qaws_scalar)-r;
	for (i = 0; i < nu; i++)
	{
		unsigned int i1 = (i + 1) % nu;
		f[k++] = n - 2; f[k++] = i; f[k++] = i1;
		for (j = 1; j + 1 < nv; j++)
		{
			unsigned int a = (j - 1) * nu + i, b = (j - 1) * nu + i1, c = j * nu + i1, e = j * nu + i;
			f[k++] = a; f[k++] = e; f[k++] = c;
			f[k++] = a; f[k++] = c; f[k++] = b;
		}
		f[k++] = n - 1; f[k++] = (nv - 2) * nu + i1; f[k++] = (nv - 2) * nu + i;
	}
	memset(&d, 0, sizeof d);
	d.positions = x; d.vertex_count = n; d.face_vertices = f; d.face_count = k / 3;
	qaws_mesh_3d_create(&m, &d);
	free(x); free(f);
	return m;
}

static void test_passes(int slabs)
{
	qaws_mesh_3d* c = si_ball(0.2, 24, 12), *o = NULL;
	qaws_motion* m = NULL;
	qaws_motion_pose s, e;
	qaws_sweep_desc d;
	qaws_attribute_track tr;
	qaws_attribute_combine cb;
	qaws_attribute_desc ad;
	unsigned int nb = qaws_mesh_3d_get_vertex_count(c), i, id, nc = 0, two = 0, one = 0, more = 0;
	qaws_scalar *zero = (qaws_scalar*)calloc(nb + 1, sizeof(qaws_scalar)), *keys = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * (nb + 1)), times[2] = { 0, 1 }, *pv = NULL, *hv = NULL;
	double wh = 0;
	char msg[300];
	for (i = 0; i < 2 * nb; i++)
		keys[i] = (qaws_scalar)(i < nb ? 0 : 1);
	memset(&ad, 0, sizeof ad);
	ad.name = "heat"; ad.element = QAWS_MESH_VERTEX; ad.components = 1; ad.rule = QAWS_ATTRIBUTE_LINEAR;
	qaws_mesh_3d_add_attribute(&id, c, zero, &ad);
	memset(&s, 0, sizeof s);
	memset(&e, 0, sizeof e);
	s.rotation[0] = e.rotation[0] = 1;
	s.velocity.x = 3;
	e.velocity.x = -3;
	qaws_motion_create_poses(&m, &s, &e, 0);
	tr.name = "heat"; tr.times = times; tr.values = keys; tr.key_count = 2;
	cb.name = "heat"; cb.combine = QAWS_SWEEP_COMBINE_SUM;
	memset(&d, 0, sizeof d);
	d.resolution = 24;
	d.tracks = &tr; d.track_count = 1;
	d.combines = &cb; d.combine_count = 1;
	d.generate = QAWS_SWEEP_GENERATE_PASSES;
	if (slabs)
	{
		d.tolerance = (qaws_scalar)0.01;
		qaws_mesh_3d_sweep(&o, NULL, c, m, &d);
	}
	else
		qaws_mesh_3d_sweep_implicit(&o, NULL, c, m, &d);
	if (o)
	{
		nc = qaws_mesh_3d_get_corner_count(o);
		pv = si_corner_values(o, "passes", 1); hv = si_corner_values(o, "heat", 1);
	}
	for (i = 0; pv && hv && i < nc; i++)
	{
		if (pv[i] == 2)
		{
			two++;
			wh = fmax(wh, fabs((double)hv[i] - 1));
		}
		else if (pv[i] == 1)
			one++;
		else
			more++;
	}
	sprintf(msg, "passes there and back%s: %u corners passed twice (heat summed off 1 by %.3f), %u once, %u other", slabs ? " (by slabs)" : "", two, wh, one, more);
	printf("  %s\n", msg);
	TEST_ASSERT(o && pv && hv && two > nc / 2 && wh < 0.05 && !more, msg);
	free(pv); free(hv); free(zero); free(keys);
	qaws_motion_destroy(m);
	qaws_mesh_3d_destroy(c); qaws_mesh_3d_destroy(o);
}

/* the first touch of a point at x on the top face: the cube's front reaches x at 1 - (2 - x)^(1/3) (0 inside it at rest) */
static double si_first_touch(double x)
{
	return x <= 1 ? 0 : 1 - cbrt(2 - x);
}

/* a sweep's time per texel: a unit cube moved by 1 along x, slowing to a stop (x = 1 - (1 - t)^3); its swept solid's
   shape, a 2 x 1 x 1 box of 12 triangles, flattened into charts and its "time" evaluated on 128 x 128 texels
   (qaws_mesh_3d_sweep_texture): on the top face a point at x is first touched at 1 - (2 - x)^(1/3), which two
   triangles cannot hold by their vertices; read from the texture at a lattice of points inside the top triangles
   (past x 1.05 and short of 1.95) it is within 0.01 of it, the vertices' exact times interpolated far from it */
static void test_sweep_texture(void)
{
	enum { R = 128 };
	qaws_mesh_3d *b = si_box(0, 1, 0, 1, 0, 1), *o = si_box(0, 2, 0, 1, 0, 1);
	qaws_motion* m = NULL;
	qaws_motion_pose s, e;
	qaws_sweep_desc d;
	qaws_texture t;
	qaws_scalar* tex = (qaws_scalar*)calloc(R * R, sizeof(qaws_scalar));
	unsigned char* cov = (unsigned char*)malloc(R * R);
	unsigned int nt = qaws_mesh_3d_get_triangle_count(o), id, written = 0, pts = 0, i, k, q;
	qaws_scalar* uv = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 6 * (nt + 1));
	unsigned int* tri = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (nt + 1));
	double e_tex = 0, e_vert = 0;
	qaws_param_desc pd;
	qaws_attribute_desc ad;
	qaws_status st;
	char msg[400];
	memset(&s, 0, sizeof s);
	memset(&e, 0, sizeof e);
	s.rotation[0] = e.rotation[0] = 1;
	e.position.x = 1;
	s.velocity.x = 3;
	qaws_motion_create_poses(&m, &s, &e, 0);
	memset(&d, 0, sizeof d);
	d.tolerance = (qaws_scalar)0.002;
	memset(&pd, 0, sizeof pd);
	pd.auto_seams = 1;
	pd.pack = 1;
	st = qaws_mesh_3d_parameterize(uv, NULL, o, NULL, 0, &pd);
	memset(&ad, 0, sizeof ad);
	ad.name = "uv"; ad.element = QAWS_MESH_CORNER; ad.components = 2; ad.rule = QAWS_ATTRIBUTE_LINEAR;
	if (st == QAWS_STATUS_OK)
		st = qaws_mesh_3d_add_attribute(&id, o, uv, &ad);
	memset(&t, 0, sizeof t);
	t.texels = tex; t.width = R; t.height = R; t.components = 1;
	if (st == QAWS_STATUS_OK)
		st = qaws_mesh_3d_sweep_texture(cov, &t, o, b, m, "time", &d, NULL);
	for (i = 0; st == QAWS_STATUS_OK && i < R * R; i++)
		written += cov[i] == 1;
	qaws_mesh_3d_get_triangles(tri, o);
	for (i = 0; st == QAWS_STATUS_OK && i < nt; i++)
	{
		qaws_vec3 pv[3];
		for (k = 0; k < 3; k++)
			pv[k] = qaws_mesh_3d_get_vertex(o, tri[3 * i + k]);
		if (fabs((double)pv[0].z - 1) > 1e-12 || fabs((double)pv[1].z - 1) > 1e-12 || fabs((double)pv[2].z - 1) > 1e-12)
			continue;
		for (q = 0; q < 64; q++)
		{
			/* (a lattice of points inside the triangle: weights a / 10, b / 10 and the rest, none under 0.1) */
			double bw[3] = { (q % 8 + 1) / 10.0, (q / 8 + 1) / 10.0, 0 }, x = 0, want, vert = 0;
			qaws_scalar suv[2] = { 0, 0 }, y[1];
			bw[2] = 1 - bw[0] - bw[1];
			if (bw[2] < 0.1 - 1e-12)
				continue;
			for (k = 0; k < 3; k++)
			{
				x += bw[k] * (double)pv[k].x;
				suv[0] += (qaws_scalar)(bw[k] * (double)uv[6 * i + 2 * k]);
				suv[1] += (qaws_scalar)(bw[k] * (double)uv[6 * i + 2 * k + 1]);
				vert += bw[k] * si_first_touch((double)pv[k].x);
			}
			if (x < 1.05 || x > 1.95)
				continue;
			want = si_first_touch(x);
			qaws_texture_sample(y, &t, suv, 1);
			e_tex = fmax(e_tex, fabs((double)y[0] - want));
			e_vert = fmax(e_vert, fabs(vert - want));
			pts++;
		}
	}
	sprintf(msg, "a sweep's time on a 12-triangle box's 128 x 128 texels: %u written; at %u top points off 1 - (2 - x)^(1/3) by %.4f (its vertices' "
		"exact times interpolated by %.4f)", written, pts, e_tex, e_vert);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && written > 0 && pts > 0 && e_tex < 0.01 && e_vert > 5 * e_tex, msg);
	free(tex); free(cov); free(uv); free(tri);
	qaws_motion_destroy(m);
	qaws_mesh_3d_destroy(b); qaws_mesh_3d_destroy(o);
}

int test_184_sweep_implicit_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 184: implicit sweeps\n");
	test_rigid();
	test_vertex_motion();
	test_attributes();
	test_generate();
	test_track();
	test_passes(0);
	test_passes(1);
	test_sweep_texture();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
