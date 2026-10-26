/*
 * Swept solids (qaws_mesh_sweep.h): time cut into slabs by the chords'
 * deviation, each slab's front faces extruded, all united by the Boolean.
 */

#include "../internal/qaws_internal_memory.h"
#include "../qaws_mesh_sweep.h"
#include "qaws_mesh_internal.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct swp_ctx
{
	qaws_motion const* m;
	unsigned int nv;
	qaws_scalar const* body;   /* a rigid motion's points (the mesh's positions) */
	qaws_scalar* buf;          /* scratch: nv points */
} swp_ctx;

/* the vertices at t, in double */
static qaws_status swp_at(swp_ctx* c, double t, double* X)
{
	unsigned int i;
	qaws_status st = qaws_motion_transform_points(c->buf, NULL, c->m, (qaws_scalar)t, qaws_motion_is_per_vertex(c->m) ? NULL : c->body, c->nv);
	for (i = 0; st == QAWS_STATUS_OK && i < 3 * c->nv; i++) X[i] = (double)c->buf[i];
	return st;
}

static int swp_cmp_u64(void const* a, void const* b)
{
	unsigned long long x = *(unsigned long long const*)a, y = *(unsigned long long const*)b;
	return x < y ? -1 : x > y;
}

static int swp_cmp_double(void const* a, void const* b)
{
	double x = *(double const*)a, y = *(double const*)b;
	return x < y ? -1 : x > y;
}

qaws_status qaws_mesh_3d_sweep(qaws_mesh_3d** out, qaws_sweep_report* report, qaws_mesh_3d const* mesh, qaws_motion const* motion,
	qaws_sweep_desc const* desc)
{
	qaws_sweep_report local;
	swp_ctx c;
	unsigned int nv, nt, i, k, *t = NULL, *twin = NULL, ntimes = 0, captimes = 0, maxs, ns;
	unsigned long long* key = NULL;
	double *Xa = NULL, *Xb = NULL, *Xm = NULL, *times = NULL, tol, lo[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL }, hi[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL }, diag = 0;
	qaws_scalar *px = NULL, *sx = NULL;
	unsigned int *sf = NULL, snv = 0, snf = 0, capv = 0, capf = 0;
	unsigned char *front = NULL, *fr = NULL;
	qaws_intersect_desc id;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (!report)
		report = &local;
	memset(report, 0, sizeof *report);
	if (!out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	if (!mesh || !motion)
		return QAWS_STATUS_INVALID_ARGUMENT;
	nv = qaws_mesh_3d_get_vertex_count(mesh);
	nt = qaws_mesh_3d_get_triangle_count(mesh);
	if (!nt || (qaws_motion_is_per_vertex(motion) && qaws_motion_get_vertex_count(motion) != nv))
		return QAWS_STATUS_INVALID_ARGUMENT;
	memset(&c, 0, sizeof c);
	c.m = motion;
	c.nv = nv;
	px = (qaws_scalar*)qaws_internal_scratch_malloc(sizeof(qaws_scalar) * (3 * (size_t)nv + 1));
	c.buf = (qaws_scalar*)qaws_internal_scratch_malloc(sizeof(qaws_scalar) * (3 * (size_t)nv + 1));
	t = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * (3 * (size_t)nt + 1));
	twin = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * (3 * (size_t)nt + 1));
	key = (unsigned long long*)qaws_internal_scratch_malloc(sizeof(unsigned long long) * (6 * (size_t)nt + 1));
	Xa = (double*)qaws_internal_scratch_malloc(sizeof(double) * (3 * (size_t)nv + 1));
	Xb = (double*)qaws_internal_scratch_malloc(sizeof(double) * (3 * (size_t)nv + 1));
	Xm = (double*)qaws_internal_scratch_malloc(sizeof(double) * (3 * (size_t)nv + 1));
	front = (unsigned char*)qaws_internal_scratch_malloc((size_t)nt + 1);
	if (!px || !c.buf || !t || !twin || !key || !Xa || !Xb || !Xm || !front)
		goto done;
	qaws_mesh_3d_get_positions(px, mesh);
	qaws_mesh_3d_get_triangles(t, mesh);
	c.body = px;
	maxs = desc && desc->max_slabs ? desc->max_slabs : 256;
	/* the half-edges' twins: each directed edge (u, v) keyed, (v, u) found */
	for (i = 0; i < 3 * nt; i++)
	{
		unsigned int u = t[i], v = t[3 * (i / 3) + (i % 3 + 1) % 3];
		key[2 * (size_t)i] = ((unsigned long long)u << 32) | v;
		key[2 * (size_t)i + 1] = i;
	}
	qsort(key, 3 * (size_t)nt, 2 * sizeof(unsigned long long), swp_cmp_u64);
	for (i = 0; i < 3 * nt; i++)
	{
		unsigned int u = t[i], v = t[3 * (i / 3) + (i % 3 + 1) % 3];
		unsigned long long want = ((unsigned long long)v << 32) | u;
		size_t a = 0, b = 3 * (size_t)nt;
		while (a < b)
		{
			size_t mid = a + (b - a) / 2;
			if (key[2 * mid] < want) a = mid + 1; else b = mid;
		}
		twin[i] = a < 3 * (size_t)nt && key[2 * a] == want ? (unsigned int)(key[2 * a + 1] / 3) : 0xFFFFFFFFu;
	}
	/* the tolerance from the solid's box at 0 */
	if ((st = swp_at(&c, 0, Xa)) != QAWS_STATUS_OK)
		goto done;
	for (i = 0; i < nv; i++)
		for (k = 0; k < 3; k++)
		{
			lo[k] = fmin(lo[k], Xa[3 * (size_t)i + k]);
			hi[k] = fmax(hi[k], Xa[3 * (size_t)i + k]);
		}
	for (k = 0; k < 3; k++) diag += (hi[k] - lo[k]) * (hi[k] - lo[k]);
	diag = sqrt(diag);
	tol = desc && desc->tolerance > 0 ? (double)desc->tolerance : 1e-3 * diag;
	/* the slabs: [0, 1] halved until the chords' deviation (at the middle and quarters) is within tol */
	{
		double* stack = (double*)qaws_internal_scratch_malloc(sizeof(double) * 2 * (2 * (size_t)maxs + 2));
		unsigned int nsk = 0;
		if (!stack)
			goto oom;
		stack[nsk++] = 0;
		stack[nsk++] = 1;
		ns = 1;
		while (nsk)
		{
			double b = stack[--nsk], a = stack[--nsk], dev = 0;
			unsigned int q;
			if ((st = swp_at(&c, a, Xa)) != QAWS_STATUS_OK || (st = swp_at(&c, b, Xb)) != QAWS_STATUS_OK)
				break;
			for (q = 1; q <= 3 && st == QAWS_STATUS_OK; q++)
			{
				double f = 0.25 * q;
				if ((st = swp_at(&c, a + f * (b - a), Xm)) != QAWS_STATUS_OK)
					break;
				for (i = 0; i < nv; i++)
				{
					double d2 = 0;
					for (k = 0; k < 3; k++)
					{
						double e = Xm[3 * (size_t)i + k] - ((1 - f) * Xa[3 * (size_t)i + k] + f * Xb[3 * (size_t)i + k]);
						d2 += e * e;
					}
					dev = fmax(dev, sqrt(d2));
				}
			}
			if (st != QAWS_STATUS_OK)
				break;
			if (dev > tol && ns < maxs)
			{
				stack[nsk++] = 0.5 * (a + b);
				stack[nsk++] = b;
				stack[nsk++] = a;
				stack[nsk++] = 0.5 * (a + b);
				ns++;
				continue;
			}
			report->deviation = (qaws_scalar)fmax((double)report->deviation, dev);
			if (ntimes + 2 > captimes)
			{
				unsigned int cap = captimes ? 2 * captimes : 64;
				double* g = (double*)qaws_internal_realloc(times, sizeof(double) * cap);
				if (!g)
				{
					qaws_internal_free(stack);
					goto oom;
				}
				times = g;
				captimes = cap;
			}
			times[ntimes++] = a;
			times[ntimes++] = b;
		}
		qaws_internal_free(stack);
		if (st != QAWS_STATUS_OK)
			goto done;
	}
	/* the slabs in time order (pairs sorted by their start) */
	qsort(times, ntimes / 2, 2 * sizeof(double), swp_cmp_double);
	report->slabs = ntimes / 2;
	/* the soup on layers of vertices, one per slab end (the solid at 0 on layer 0): each slab's front faces extruded between
	   its two layers; a face front in two slabs in a row is not closed between them (the end of one and the start of the
	   next would coincide, turned), so consecutive extrusions join */
	{
		size_t S = report->slabs, need_v = (size_t)nv * (S + 1), need_f = (size_t)nt * (1 + 8 * S);   /* (per slab a front face's two caps and 6 wall triangles at most) */
		if (need_v > 0xFFFFFFF0u || need_f > 0xFFFFFFF0u)
			goto oom;
		capv = (unsigned int)need_v;
		capf = (unsigned int)need_f;
		sx = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)capv + 1));
		sf = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)capf + 1));
		fr = (unsigned char*)qaws_internal_scratch_malloc(qaws_internal_mul2((size_t)nt, S) + 1);
		if (!sx || !sf || !fr)
			goto oom;
	}
	for (ns = 0; ns <= report->slabs; ns++)
	{
		/* layer ns: the vertices at the slab ends */
		if ((st = swp_at(&c, ns < report->slabs ? times[2 * (size_t)ns] : times[2 * (size_t)ns - 1], Xa)) != QAWS_STATUS_OK)
			goto done;
		for (i = 0; i < 3 * nv; i++) sx[3 * (size_t)ns * nv + i] = (qaws_scalar)Xa[i];
	}
	snv = nv * (report->slabs + 1);
	memcpy(sf, t, sizeof(unsigned int) * 3 * (size_t)nt);
	snf = nt;
	/* front faces per slab: the normal at its start against the vertices' mean displacement over it */
	for (ns = 0; ns < report->slabs; ns++)
	{
		unsigned int f;
		qaws_scalar const *A = sx + 3 * (size_t)ns * nv, *Bv = sx + 3 * ((size_t)ns + 1) * nv;
		unsigned char* fs = fr + (size_t)ns * nt;
		for (f = 0; f < nt; f++)
		{
			double u[3], w[3], n[3], d[3] = { 0, 0, 0 }, dn;
			unsigned int a = t[3 * f], b = t[3 * f + 1], e = t[3 * f + 2];
			for (k = 0; k < 3; k++)
			{
				u[k] = (double)A[3 * (size_t)b + k] - (double)A[3 * (size_t)a + k];
				w[k] = (double)A[3 * (size_t)e + k] - (double)A[3 * (size_t)a + k];
			}
			n[0] = u[1] * w[2] - u[2] * w[1];
			n[1] = u[2] * w[0] - u[0] * w[2];
			n[2] = u[0] * w[1] - u[1] * w[0];
			for (i = 0; i < 3; i++)
				for (k = 0; k < 3; k++) d[k] += ((double)Bv[3 * (size_t)t[3 * f + i] + k] - (double)A[3 * (size_t)t[3 * f + i] + k]) / 3;
			dn = n[0] * d[0] + n[1] * d[1] + n[2] * d[2];
			fs[f] = dn > 1e-14 * sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]) * sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
			report->front_faces += fs[f];
		}
	}
	for (ns = 0; ns < report->slabs; ns++)
	{
		unsigned int f, oa = ns * nv, ob = (ns + 1) * nv;
		unsigned char const *fs = fr + (size_t)ns * nt, *fp = ns ? fs - nt : NULL, *fn = ns + 1 < report->slabs ? fs + nt : NULL;
		for (f = 0; f < nt; f++)
		{
			unsigned int j;
			if (!fs[f])
				continue;
			/* the face at the start turned (unless the slab before extruded it), at the end (unless the next one does) */
			if (!fp || !fp[f])
			{
				sf[3 * (size_t)snf] = oa + t[3 * f]; sf[3 * (size_t)snf + 1] = oa + t[3 * f + 2]; sf[3 * (size_t)snf + 2] = oa + t[3 * f + 1]; snf++;
			}
			if (!fn || !fn[f])
			{
				sf[3 * (size_t)snf] = ob + t[3 * f]; sf[3 * (size_t)snf + 1] = ob + t[3 * f + 1]; sf[3 * (size_t)snf + 2] = ob + t[3 * f + 2]; snf++;
			}
			/* walls on the border of the front region: edge u -> v, the quad (u, v, v', u') */
			for (j = 0; j < 3; j++)
			{
				unsigned int g = twin[3 * f + j], u = t[3 * f + j], v = t[3 * f + (j + 1) % 3];
				if (g != 0xFFFFFFFFu && fs[g])
					continue;
				sf[3 * (size_t)snf] = oa + u; sf[3 * (size_t)snf + 1] = oa + v; sf[3 * (size_t)snf + 2] = ob + v; snf++;
				sf[3 * (size_t)snf] = oa + u; sf[3 * (size_t)snf + 1] = ob + v; sf[3 * (size_t)snf + 2] = ob + u; snf++;
			}
		}
	}
	{
		qaws_mesh_3d_desc md;
		memset(&md, 0, sizeof md);
		md.positions = sx;
		md.vertex_count = snv;
		md.face_vertices = sf;
		md.face_count = snf;
		memset(&id, 0, sizeof id);
		id.executor = desc ? desc->executor : NULL;
		st = qaws_mesh_boolean(out, &report->boolean, &md, NULL, QAWS_BOOLEAN_UNION, &id);
		report->closed = report->boolean.closed;
	}
	goto done;
oom:
	st = QAWS_STATUS_ALLOCATION_FAILURE;
done:
	qaws_internal_free(px); qaws_internal_free(c.buf); qaws_internal_free(t); qaws_internal_free(twin); qaws_internal_free(key);
	qaws_internal_free(Xa); qaws_internal_free(Xb); qaws_internal_free(Xm); qaws_internal_free(front); qaws_internal_free(times);
	qaws_internal_free(sx); qaws_internal_free(sf); qaws_internal_free(fr);
	return st;
}

/* a closed mesh's signed volume (its triangles fanned from the origin) */
static double swp_volume(qaws_mesh_3d const* m)
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

qaws_status qaws_mesh_sweep_boolean(qaws_mesh_3d** out, qaws_sweep_boolean_report* report, qaws_mesh_3d const* tool, qaws_motion const* const* moves,
	unsigned int move_count, qaws_mesh_3d const* stock, qaws_boolean_op op, qaws_sweep_desc const* desc)
{
	qaws_sweep_boolean_report local;
	qaws_scalar *bx = NULL, *ax = NULL;
	unsigned int *bt = NULL, *at = NULL, nbv = 0, nbt = 0, i, k;
	qaws_mesh_3d_desc a, b;
	qaws_intersect_desc id;
	qaws_status st = QAWS_STATUS_OK;
	if (!report)
		report = &local;
	memset(report, 0, sizeof *report);
	if (!out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	if (!tool || !moves || !move_count || !stock ||
		(op != QAWS_BOOLEAN_DIFFERENCE && op != QAWS_BOOLEAN_UNION && op != QAWS_BOOLEAN_INTERSECTION))
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < move_count; i++)
		if (!moves[i])
			return QAWS_STATUS_INVALID_ARGUMENT;
	/* each move's sweep, into one operand */
	for (i = 0; i < move_count && st == QAWS_STATUS_OK; i++)
	{
		qaws_mesh_3d* s = NULL;
		qaws_sweep_report r;
		unsigned int nv, nt;
		st = qaws_mesh_3d_sweep(&s, &r, tool, moves[i], desc);
		report->slabs += r.slabs;
		if (st != QAWS_STATUS_OK || !s)
		{
			qaws_mesh_3d_destroy(s);
			continue;
		}
		nv = qaws_mesh_3d_get_vertex_count(s);
		nt = qaws_mesh_3d_get_triangle_count(s);
		{
			qaws_scalar* g = (qaws_scalar*)qaws_internal_realloc(bx, qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)nbv + nv + 1));
			unsigned int* h = g ? (unsigned int*)qaws_internal_realloc(bt, qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)nbt + nt + 1)) : NULL;
			if (g) bx = g;
			if (!g || !h)
			{
				qaws_mesh_3d_destroy(s);
				st = QAWS_STATUS_ALLOCATION_FAILURE;
				break;
			}
			bt = h;
		}
		qaws_mesh_3d_get_positions(bx + 3 * (size_t)nbv, s);
		qaws_mesh_3d_get_triangles(bt + 3 * (size_t)nbt, s);
		for (k = 0; k < 3 * nt; k++) bt[3 * (size_t)nbt + k] += nbv;
		nbv += nv;
		nbt += nt;
		report->moves++;
		qaws_mesh_3d_destroy(s);
	}
	if (st != QAWS_STATUS_OK)
		goto done;
	report->stock_volume = (qaws_scalar)swp_volume(stock);
	{
		unsigned int nv = qaws_mesh_3d_get_vertex_count(stock), nt = qaws_mesh_3d_get_triangle_count(stock);
		ax = (qaws_scalar*)qaws_internal_scratch_malloc(sizeof(qaws_scalar) * (3 * (size_t)nv + 1));
		at = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * (3 * (size_t)nt + 1));
		if (!ax || !at)
		{
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			goto done;
		}
		qaws_mesh_3d_get_positions(ax, stock);
		qaws_mesh_3d_get_triangles(at, stock);
		memset(&a, 0, sizeof a);
		a.positions = ax;
		a.vertex_count = nv;
		a.face_vertices = at;
		a.face_count = nt;
	}
	memset(&id, 0, sizeof id);
	id.executor = desc ? desc->executor : NULL;
	if (!nbt)
	{
		/* (the tool swept nothing: the stock, or nothing, as the op says) */
		if (op == QAWS_BOOLEAN_INTERSECTION)
			goto done;
		st = qaws_mesh_boolean(out, &report->boolean, &a, NULL, QAWS_BOOLEAN_UNION, &id);
	}
	else
	{
		memset(&b, 0, sizeof b);
		b.positions = bx;
		b.vertex_count = nbv;
		b.face_vertices = bt;
		b.face_count = nbt;
		/* (sweep op stock: the stock first, taken from for the difference) */
		st = qaws_mesh_boolean(out, &report->boolean, &a, &b, op, &id);
	}
	report->closed = report->boolean.closed;
	report->result_volume = (qaws_scalar)swp_volume(*out);
	report->changed_volume = (qaws_scalar)(op == QAWS_BOOLEAN_DIFFERENCE ? (double)report->stock_volume - (double)report->result_volume :
		op == QAWS_BOOLEAN_UNION ? (double)report->result_volume - (double)report->stock_volume : (double)report->result_volume);
done:
	qaws_internal_free(bx); qaws_internal_free(bt); qaws_internal_free(ax); qaws_internal_free(at);
	return st;
}
