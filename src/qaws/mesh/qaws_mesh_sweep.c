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

/* the soup as a mesh carrying the swept mesh's attributes: per soup face its triangle (sof; bit 31: turned, 30: a wall), per
   soup corner the triangle's corner it repeats (soc), on the layer of its vertex; vertex attributes become corner ones,
   directions (NORMAL, and NORMALIZED of 3 components) turn with a rigid motion, NORMAL ones are negated on turned faces;
   a wall's NORMAL face value is its own normal (its face's is not: a wall stands across it); a soup edge repeating an
   edge of the mesh at one layer (a cap's, a wall's top or bottom) takes its edge values and feature flag, the others
   zeros */
static qaws_status swp_soup(qaws_mesh_3d** out, qaws_mesh_3d const* mesh, qaws_motion const* motion, qaws_mesh_3d_desc const* md,
	unsigned int const* sof, unsigned int const* soc, unsigned int nv, double const* layer_t, unsigned int layers)
{
	unsigned int na = qaws_mesh_3d_get_attribute_count(mesh), nf = qaws_mesh_3d_get_face_count(mesh), nc = qaws_mesh_3d_get_corner_count(mesh);
	unsigned int nt = qaws_mesh_3d_get_triangle_count(mesh), a, i, j, k, *fv = NULL, *fs = NULL, *tf = NULL, *tc = NULL;
	unsigned int ine = qaws_mesh_3d_get_edge_count(mesh), sne, *ie = NULL, *se = NULL, *smap = NULL, nfe = 0;
	unsigned long long* ik = NULL;
	unsigned char* ifeat = NULL;
	qaws_scalar *val = NULL, *sv = NULL;
	double* rot = NULL;
	qaws_status st;
	int rigid = !qaws_motion_is_per_vertex(motion);
	if ((st = qaws_mesh_3d_create(out, md)) != QAWS_STATUS_OK)
		return st;
	fv = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * ((size_t)nc + 1));
	fs = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * ((size_t)nf + 1));
	tf = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * ((size_t)nt + 1));
	tc = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)nt + 1));
	val = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 16, (size_t)(nv > nc ? (nv > nf ? nv : nf) : (nc > nf ? nc : nf)) + 1));
	sv = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 16, 3 * (size_t)md->face_count + 1));
	rot = (double*)qaws_internal_scratch_malloc(sizeof(double) * (9 * (size_t)layers + 1));
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (!fv || !fs || !tf || !tc || !val || !sv || !rot)
		goto done;
	/* the soup's edges on the mesh's (one layer, an edge of the mesh), their features */
	sne = qaws_mesh_3d_get_edge_count(*out);
	{
		unsigned int nfeat = qaws_mesh_3d_get_feature_edge_count(mesh), m = ine > nfeat ? ine : nfeat;
		ie = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 2, (size_t)m + 1));
		se = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 2, (size_t)sne + 1));
		smap = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * ((size_t)sne + 1));
		ifeat = (unsigned char*)qaws_internal_zalloc((size_t)ine + 1, 1);
		if (!ie || !se || !smap || !ifeat)
			goto done;
		qaws_mesh_3d_get_edges(ie, mesh);
		if (!(ik = cy_keys(ie, ine)))
			goto done;
		qaws_mesh_3d_get_feature_edges(ie, mesh);
		for (i = 0; i < nfeat; i++)
		{
			unsigned int x = cy_find(ik, ine, ie[2 * (size_t)i], ie[2 * (size_t)i + 1]);
			if (x != QAWS_MESH_NONE)
				ifeat[x] = 1;
		}
		qaws_mesh_3d_get_edges(se, *out);
		for (i = 0; i < sne; i++)
		{
			unsigned int x = se[2 * (size_t)i], y = se[2 * (size_t)i + 1];
			smap[i] = x / nv == y / nv ? cy_find(ik, ine, x % nv, y % nv) : QAWS_MESH_NONE;
			if (smap[i] != QAWS_MESH_NONE && ifeat[smap[i]])
			{
				se[2 * (size_t)nfe] = x;
				se[2 * (size_t)nfe + 1] = y;
				nfe++;
			}
		}
		if (nfe && qaws_mesh_3d_set_feature_edges(*out, se, nfe) != QAWS_STATUS_OK)
			goto done;
	}
	/* the triangles (fanned from each face's first corner) back to their faces' corners */
	qaws_mesh_3d_get_faces(fv, fs, mesh);
	for (i = 0, j = 0, k = 0; i < nf; k += fs[i++])
		for (a = 0; a + 2 < fs[i]; a++, j++)
		{
			tf[j] = i;
			tc[3 * (size_t)j] = k;
			tc[3 * (size_t)j + 1] = k + a + 1;
			tc[3 * (size_t)j + 2] = k + a + 2;
		}
	for (i = 0; i < layers; i++)
	{
		qaws_motion_state s;
		if (rigid && qaws_motion_evaluate(&s, motion, (qaws_scalar)layer_t[i]) == QAWS_STATUS_OK)
			for (k = 0; k < 9; k++)
				rot[9 * (size_t)i + k] = (double)s.rotation[k];
		else
			for (k = 0; k < 9; k++)
				rot[9 * (size_t)i + k] = k % 4 == 0;
	}
	st = QAWS_STATUS_OK;
	for (a = 0; a < na && st == QAWS_STATUS_OK; a++)
	{
		qaws_attribute_desc d;
		unsigned int n, id, cnt;
		int dir;
		qaws_mesh_3d_get_attribute_desc(&d, mesh, a);
		if (d.name && !strcmp(d.name, "source"))
			continue;
		n = d.components;
		dir = d.rule == QAWS_ATTRIBUTE_NORMAL || d.rule == QAWS_ATTRIBUTE_TANGENT || (d.rule == QAWS_ATTRIBUTE_NORMALIZED && n == 3);
		qaws_mesh_3d_get_attribute(val, mesh, a);
		if (d.element == QAWS_MESH_EDGE)
		{
			for (i = 0; i < sne; i++)
				for (k = 0; k < n; k++)
					sv[(size_t)n * i + k] = smap[i] != QAWS_MESH_NONE ? val[(size_t)n * smap[i] + k] : 0;
			st = qaws_mesh_3d_add_attribute(&id, *out, sv, &d);
			continue;
		}
		cnt = d.element == QAWS_MESH_FACE ? md->face_count : 3 * md->face_count;
		for (i = 0; i < cnt; i++)
		{
			unsigned int g = d.element == QAWS_MESH_FACE ? sof[i] : sof[i / 3], f = g & 0x3FFFFFFFu, turned = g >> 31;
			unsigned int c = d.element == QAWS_MESH_FACE ? 0 : tc[soc[i]], e = d.element == QAWS_MESH_FACE ? tf[f] : d.element == QAWS_MESH_VERTEX ? fv[c] : c;
			unsigned int layer = d.element == QAWS_MESH_FACE ? 0 : md->face_vertices[i] / nv;
			qaws_scalar const* x = val + (size_t)n * e;
			qaws_scalar* o = sv + (size_t)n * i;
			for (k = 0; k < n; k++)
				o[k] = x[k];
			if (dir && d.element != QAWS_MESH_FACE)
			{
				double const* R = rot + 9 * (size_t)layer;
				for (k = 0; k < 3; k++)
					o[k] = (qaws_scalar)(R[3 * k] * (double)x[0] + R[3 * k + 1] * (double)x[1] + R[3 * k + 2] * (double)x[2]);
			}
			if (d.rule == QAWS_ATTRIBUTE_NORMAL && turned)
				for (k = 0; k < 3; k++)
					o[k] = -o[k];
			if (d.rule == QAWS_ATTRIBUTE_TANGENT && turned)
				o[3] = -o[3];
		}
		for (i = 0; d.rule == QAWS_ATTRIBUTE_NORMAL && d.element == QAWS_MESH_FACE && i < cnt; i++)
			if (sof[i] >> 30 & 1)
			{
				qaws_scalar const *p0 = md->positions + 3 * (size_t)md->face_vertices[3 * (size_t)i], *p1 = md->positions + 3 * (size_t)md->face_vertices[3 * (size_t)i + 1],
					*p2 = md->positions + 3 * (size_t)md->face_vertices[3 * (size_t)i + 2];
				double u[3], w[3], nn[3], l;
				for (k = 0; k < 3; k++)
				{
					u[k] = (double)p1[k] - (double)p0[k];
					w[k] = (double)p2[k] - (double)p0[k];
				}
				nn[0] = u[1] * w[2] - u[2] * w[1];
				nn[1] = u[2] * w[0] - u[0] * w[2];
				nn[2] = u[0] * w[1] - u[1] * w[0];
				l = sqrt(nn[0] * nn[0] + nn[1] * nn[1] + nn[2] * nn[2]);
				for (k = 0; k < 3; k++)
					sv[3 * (size_t)i + k] = (qaws_scalar)(l > 0 ? nn[k] / l : 0);
			}
		/* (a face attribute's directions: the motion at 0, as the solid at 0's) */
		if (dir && d.element == QAWS_MESH_FACE)
			for (i = 0; i < cnt; i++)
			{
				qaws_scalar* o = sv + (size_t)n * i;
				double y[3];
				for (k = 0; k < 3; k++)
					y[k] = rot[3 * k] * (double)o[0] + rot[3 * k + 1] * (double)o[1] + rot[3 * k + 2] * (double)o[2];
				for (k = 0; k < 3; k++)
					o[k] = (qaws_scalar)y[k];
			}
		if (d.element == QAWS_MESH_VERTEX)
			d.element = QAWS_MESH_CORNER;
		st = qaws_mesh_3d_add_attribute(&id, *out, sv, &d);
	}
done:
	if (st != QAWS_STATUS_OK)
	{
		qaws_mesh_3d_destroy(*out);
		*out = NULL;
	}
	qaws_internal_free(fv); qaws_internal_free(fs); qaws_internal_free(tf); qaws_internal_free(tc); qaws_internal_free(val);
	qaws_internal_free(sv); qaws_internal_free(rot); qaws_internal_free(ie); qaws_internal_free(se); qaws_internal_free(smap);
	qaws_internal_free(ik); qaws_internal_free(ifeat);
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
	unsigned int *sf = NULL, *sof = NULL, *soc = NULL, snv = 0, snf = 0, capv = 0, capf = 0;
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
		sof = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)capf + 1));
		soc = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)capf + 1));
		if (!sx || !sf || !fr || !sof || !soc)
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
	for (i = 0; i < 3 * nt; i++)
	{
		sof[i / 3] = i / 3;
		soc[i] = i;
	}
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
				sf[3 * (size_t)snf] = oa + t[3 * f]; sf[3 * (size_t)snf + 1] = oa + t[3 * f + 2]; sf[3 * (size_t)snf + 2] = oa + t[3 * f + 1];
				sof[snf] = f | 0x80000000u;
				soc[3 * (size_t)snf] = 3 * f; soc[3 * (size_t)snf + 1] = 3 * f + 2; soc[3 * (size_t)snf + 2] = 3 * f + 1; snf++;
			}
			if (!fn || !fn[f])
			{
				sf[3 * (size_t)snf] = ob + t[3 * f]; sf[3 * (size_t)snf + 1] = ob + t[3 * f + 1]; sf[3 * (size_t)snf + 2] = ob + t[3 * f + 2];
				sof[snf] = f;
				soc[3 * (size_t)snf] = 3 * f; soc[3 * (size_t)snf + 1] = 3 * f + 1; soc[3 * (size_t)snf + 2] = 3 * f + 2; snf++;
			}
			/* walls on the border of the front region: edge u -> v, the quad (u, v, v', u') */
			for (j = 0; j < 3; j++)
			{
				unsigned int g = twin[3 * f + j], u = t[3 * f + j], v = t[3 * f + (j + 1) % 3];
				if (g != 0xFFFFFFFFu && fs[g])
					continue;
				sf[3 * (size_t)snf] = oa + u; sf[3 * (size_t)snf + 1] = oa + v; sf[3 * (size_t)snf + 2] = ob + v;
				sof[snf] = f | 0x40000000u;
				soc[3 * (size_t)snf] = 3 * f + j; soc[3 * (size_t)snf + 1] = 3 * f + (j + 1) % 3; soc[3 * (size_t)snf + 2] = 3 * f + (j + 1) % 3; snf++;
				sf[3 * (size_t)snf] = oa + u; sf[3 * (size_t)snf + 1] = ob + v; sf[3 * (size_t)snf + 2] = ob + u;
				sof[snf] = f | 0x40000000u;
				soc[3 * (size_t)snf] = 3 * f + j; soc[3 * (size_t)snf + 1] = 3 * f + (j + 1) % 3; soc[3 * (size_t)snf + 2] = 3 * f + j; snf++;
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
		{
			/* the soup as a mesh carrying the attributes and features, its self-union carrying them on (the cuts features too) */
			qaws_mesh_3d* soup = NULL;
			for (ns = 0; ns <= report->slabs; ns++)
				times[ns] = ns < report->slabs ? times[2 * (size_t)ns] : times[2 * (size_t)ns - 1];   /* (each layer's time, in place) */
			st = swp_soup(&soup, mesh, motion, &md, sof, soc, nv, times, report->slabs + 1);
			if (st == QAWS_STATUS_OK)
				st = cy_boolean(out, &report->boolean, (qaws_mesh_3d const* const*)&soup, 1, NULL, 0, QAWS_BOOLEAN_UNION, &id);
			qaws_mesh_3d_destroy(soup);
		}
		report->closed = report->boolean.closed;
	}
	goto done;
oom:
	st = QAWS_STATUS_ALLOCATION_FAILURE;
done:
	qaws_internal_free(px); qaws_internal_free(c.buf); qaws_internal_free(t); qaws_internal_free(twin); qaws_internal_free(key);
	qaws_internal_free(Xa); qaws_internal_free(Xb); qaws_internal_free(Xm); qaws_internal_free(front); qaws_internal_free(times);
	qaws_internal_free(sx); qaws_internal_free(sf); qaws_internal_free(fr); qaws_internal_free(sof); qaws_internal_free(soc);
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
	qaws_mesh_3d const* stock, unsigned int move_count, qaws_boolean_op op, qaws_sweep_desc const* desc)
{
	qaws_sweep_boolean_report local;
	qaws_mesh_3d** sw = NULL;
	unsigned int i, ns = 0;
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
	sw = (qaws_mesh_3d**)qaws_internal_zalloc(move_count, sizeof *sw);
	if (!sw)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	/* each move's sweep (its attributes the tool's), all of them one operand */
	for (i = 0; i < move_count && st == QAWS_STATUS_OK; i++)
	{
		qaws_sweep_report r;
		st = qaws_mesh_3d_sweep(&sw[ns], &r, tool, moves[i], desc);
		report->slabs += r.slabs;
		if (st == QAWS_STATUS_OK && sw[ns])
		{
			ns++;
			report->moves++;
		}
	}
	if (st != QAWS_STATUS_OK)
		goto done;
	report->stock_volume = (qaws_scalar)swp_volume(stock);
	memset(&id, 0, sizeof id);
	id.executor = desc ? desc->executor : NULL;
	if (!ns)
	{
		/* (the tool swept nothing: the stock, or nothing, as the op says) */
		if (op == QAWS_BOOLEAN_INTERSECTION)
			goto done;
		st = cy_boolean(out, &report->boolean, &stock, 1, NULL, 0, QAWS_BOOLEAN_UNION, &id);
	}
	else
		st = cy_boolean(out, &report->boolean, &stock, 1, (qaws_mesh_3d const* const*)sw, ns, op, &id);   /* (sweep op stock: the stock first, taken from for the difference) */
	report->closed = report->boolean.closed;
	report->result_volume = (qaws_scalar)swp_volume(*out);
	report->changed_volume = (qaws_scalar)(op == QAWS_BOOLEAN_DIFFERENCE ? (double)report->stock_volume - (double)report->result_volume :
		op == QAWS_BOOLEAN_UNION ? (double)report->result_volume - (double)report->stock_volume : (double)report->result_volume);
done:
	for (i = 0; i < ns; i++)
		qaws_mesh_3d_destroy(sw[i]);
	qaws_internal_free(sw);
	return st;
}
