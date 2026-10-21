/*
 * The mesh connector (qaws_mesh_connect.h): border loops walked on the
 * corner store, matched by arc length, rings Hermite-blended between them
 * along the meshes' surfaces, zipped into faces, the result built anew.
 */

#include "../internal/qaws_internal_memory.h"
#include "../qaws_mesh_connect.h"
#include "../qaws_mesh_smooth.h"
#include "qaws_mesh_internal.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* a border loop: its vertices (in the border's direction), positions, the directions leaving the mesh across it, and its arc
   length's share at each vertex */
typedef struct cn_loop
{
	unsigned int n;
	unsigned int* v;
	double *x, *o, *s, h;           /* h: the mean edge length */
} cn_loop;

static void cn_loop_free(cn_loop* L)
{
	qaws_internal_free(L->v);
	qaws_internal_free(L->x);
	qaws_internal_free(L->o);
	qaws_internal_free(L->s);
	memset(L, 0, sizeof *L);
}

static void cn_cross(double const* a, double const* b, double* o)
{
	o[0] = a[1] * b[2] - a[2] * b[1];
	o[1] = a[2] * b[0] - a[0] * b[2];
	o[2] = a[0] * b[1] - a[1] * b[0];
}

static double cn_norm(double* a)
{
	double l = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
	if (l > 0)
	{
		a[0] /= l;
		a[1] /= l;
		a[2] /= l;
	}
	return l;
}

/* the loop through v0 (along the border's own corners), its geometry; INVALID_ARGUMENT when v0 is not on a border or the walk
   does not close */
static qaws_status cn_walk(cn_loop* L, qaws_mesh_3d const* mesh, unsigned int v0)
{
	qaws_hm const* m = &mesh->hm;
	unsigned int v, i, k, cap;
	qaws_scalar* nrm = NULL;
	double len = 0;
	memset(L, 0, sizeof *L);
	if (v0 >= m->nv || !hm_is_border_vertex(m, v0))
		return QAWS_STATUS_INVALID_ARGUMENT;
	cap = m->nv;
	L->v = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)cap + 1));
	if (!L->v)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	v = v0;
	do
	{
		unsigned int c = m->vc[v];
		if (L->n >= cap || c == HM_NONE || m->co[c] != HM_NONE)
		{
			cn_loop_free(L);
			return QAWS_STATUS_INVALID_ARGUMENT;
		}
		L->v[L->n++] = v;
		v = hm_to(m, c);
	} while (v != v0);
	if (L->n < 3)
	{
		cn_loop_free(L);
		return QAWS_STATUS_INVALID_ARGUMENT;
	}
	L->x = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, L->n));
	L->o = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, L->n));
	L->s = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)L->n + 1));
	nrm = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)m->nv + 1));
	if (!L->x || !L->o || !L->s || !nrm || qaws_mesh_3d_get_vertex_normals(nrm, mesh) != QAWS_STATUS_OK)
	{
		qaws_internal_free(nrm);
		cn_loop_free(L);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < L->n; i++)
		memcpy(L->x + 3 * (size_t)i, m->x + 3 * (size_t)L->v[i], 3 * sizeof(double));
	L->s[0] = 0;
	for (i = 0; i < L->n; i++)
	{
		double const *p = L->x + 3 * (size_t)i, *q = L->x + 3 * (size_t)((i + 1) % L->n);
		len += sqrt((q[0] - p[0]) * (q[0] - p[0]) + (q[1] - p[1]) * (q[1] - p[1]) + (q[2] - p[2]) * (q[2] - p[2]));
		L->s[i + 1] = len;
	}
	for (i = 0; i <= L->n; i++)
		L->s[i] /= len;
	L->h = len / L->n;
	/* leaving the mesh: across the border in the tangent plane, away from the face at the vertex */
	for (i = 0; i < L->n; i++)
	{
		double const *pp = L->x + 3 * (size_t)((i + L->n - 1) % L->n), *pn = L->x + 3 * (size_t)((i + 1) % L->n), *x = L->x + 3 * (size_t)i;
		double t[3] = { pn[0] - pp[0], pn[1] - pp[1], pn[2] - pp[2] }, n[3], cen[3] = { 0, 0, 0 }, *o = L->o + 3 * (size_t)i;
		unsigned int f = m->cf[m->vc[L->v[i]]], c;
		for (k = 0; k < 3; k++)
			n[k] = (double)nrm[3 * (size_t)L->v[i] + k];
		cn_norm(t);
		cn_cross(t, n, o);
		cn_norm(o);
		for (c = m->fc[f]; c < m->fc[f] + m->fn[f]; c++)
			for (k = 0; k < 3; k++)
				cen[k] += m->x[3 * (size_t)m->cv[c] + k] / m->fn[f];
		if (o[0] * (cen[0] - x[0]) + o[1] * (cen[1] - x[1]) + o[2] * (cen[2] - x[2]) > 0)
			for (k = 0; k < 3; k++)
				o[k] = -o[k];
	}
	qaws_internal_free(nrm);
	return QAWS_STATUS_OK;
}

/* the loop reversed and started at its vertex j0 (the shares recomputed from there) */
static qaws_status cn_turn(cn_loop* L, unsigned int j0)
{
	cn_loop R;
	unsigned int i, n = L->n;
	double len = 0;
	memset(&R, 0, sizeof R);
	R.n = n;
	R.h = L->h;
	R.v = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * n);
	R.x = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, n));
	R.o = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, n));
	R.s = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)n + 1));
	if (!R.v || !R.x || !R.o || !R.s)
	{
		cn_loop_free(&R);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < n; i++)
	{
		unsigned int j = (j0 + n - i) % n;
		R.v[i] = L->v[j];
		memcpy(R.x + 3 * (size_t)i, L->x + 3 * (size_t)j, 3 * sizeof(double));
		memcpy(R.o + 3 * (size_t)i, L->o + 3 * (size_t)j, 3 * sizeof(double));
	}
	R.s[0] = 0;
	for (i = 0; i < n; i++)
	{
		double const *p = R.x + 3 * (size_t)i, *q = R.x + 3 * (size_t)((i + 1) % n);
		len += sqrt((q[0] - p[0]) * (q[0] - p[0]) + (q[1] - p[1]) * (q[1] - p[1]) + (q[2] - p[2]) * (q[2] - p[2]));
		R.s[i + 1] = len;
	}
	for (i = 0; i <= n; i++)
		R.s[i] /= len;
	cn_loop_free(L);
	*L = R;
	return QAWS_STATUS_OK;
}

/* the loop's point and leaving direction at share s of its length (linear between its vertices) */
static void cn_at(cn_loop const* L, double s, double* x, double* o)
{
	unsigned int lo = 0, hi = L->n, k;
	double t;
	s -= floor(s);
	while (hi - lo > 1)
	{
		unsigned int mid = (lo + hi) / 2;
		if (L->s[mid] <= s)
			lo = mid;
		else
			hi = mid;
	}
	t = L->s[lo + 1] > L->s[lo] ? (s - L->s[lo]) / (L->s[lo + 1] - L->s[lo]) : 0;
	for (k = 0; k < 3; k++)
	{
		double const *p = L->x + 3 * (size_t)lo, *q = L->x + 3 * (size_t)((lo + 1) % L->n);
		double const *a = L->o + 3 * (size_t)lo, *b = L->o + 3 * (size_t)((lo + 1) % L->n);
		x[k] = p[k] + t * (q[k] - p[k]);
		if (o)
			o[k] = a[k] + t * (b[k] - a[k]);
	}
	if (o)
		cn_norm(o);
}

static double cn_angle(double const* o, double const* p, double const* q)
{
	double u[3] = { p[0] - o[0], p[1] - o[1], p[2] - o[2] }, w[3] = { q[0] - o[0], q[1] - o[1], q[2] - o[2] }, c[3];
	cn_cross(u, w, c);
	return atan2(sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]), u[0] * w[0] + u[1] * w[1] + u[2] * w[2]);
}

qaws_status qaws_mesh_3d_connect(qaws_mesh_3d** out, qaws_connect_report* report, qaws_mesh_3d const* a, unsigned int av,
	qaws_mesh_3d const* b, unsigned int bv, qaws_connect_desc const* desc)
{
	qaws_connect_report local;
	cn_loop A, B;
	unsigned int nva, nvb, R, r, i, j, k, nfa = 0, nfb = 0, ca = 0, cb = 0, nbv = 0, nbf = 0, nbc = 0, *ring_start = NULL,
		*ring_n = NULL, *fv = NULL, *fs = NULL;
	double *bx = NULL, gap = 0, *t = NULL;
	qaws_scalar* pos = NULL;
	unsigned int *afs = NULL, *bfs = NULL;
	int same = !b || b == a;
	qaws_status st;
	memset(&A, 0, sizeof A);
	memset(&B, 0, sizeof B);
	if (!report)
		report = &local;
	memset(report, 0, sizeof *report);
	if (!out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	if (!a)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (same)
		b = a;
	st = cn_walk(&A, a, av);
	if (st == QAWS_STATUS_OK)
		st = cn_walk(&B, b, bv);
	if (st != QAWS_STATUS_OK)
		goto done;
	if (same)
		for (i = 0; i < A.n; i++)
			if (A.v[i] == bv)
			{
				st = QAWS_STATUS_INVALID_ARGUMENT;   /* (one loop named twice) */
				goto done;
			}
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	{
		/* b's loop turned the other way, started where the matched points are nearest */
		unsigned int best = 0;
		double bc = HUGE_VAL;
		for (j = 0; j < B.n; j++)
		{
			cn_loop T;
			double c = 0;
			memset(&T, 0, sizeof T);
			T.n = B.n;
			T.v = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * B.n);
			T.x = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, B.n));
			T.o = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, B.n));
			T.s = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)B.n + 1));
			if (!T.v || !T.x || !T.o || !T.s)
			{
				cn_loop_free(&T);
				goto done;
			}
			memcpy(T.v, B.v, sizeof(unsigned int) * B.n);
			memcpy(T.x, B.x, sizeof(double) * 3 * B.n);
			memcpy(T.o, B.o, sizeof(double) * 3 * B.n);
			memcpy(T.s, B.s, sizeof(double) * (B.n + 1));
			if (cn_turn(&T, j) != QAWS_STATUS_OK)
			{
				cn_loop_free(&T);
				goto done;
			}
			for (i = 0; i < A.n; i++)
			{
				double y[3];
				cn_at(&T, A.s[i], y, NULL);
				c += (y[0] - A.x[3 * i]) * (y[0] - A.x[3 * i]) + (y[1] - A.x[3 * i + 1]) * (y[1] - A.x[3 * i + 1])
					+ (y[2] - A.x[3 * i + 2]) * (y[2] - A.x[3 * i + 2]);
			}
			cn_loop_free(&T);
			if (c < bc)
			{
				bc = c;
				best = j;
			}
		}
		if (cn_turn(&B, best) != QAWS_STATUS_OK)
			goto done;
	}
	for (i = 0; i < A.n; i++)
	{
		double y[3];
		cn_at(&B, A.s[i], y, NULL);
		gap += sqrt((y[0] - A.x[3 * i]) * (y[0] - A.x[3 * i]) + (y[1] - A.x[3 * i + 1]) * (y[1] - A.x[3 * i + 1])
			+ (y[2] - A.x[3 * i + 2]) * (y[2] - A.x[3 * i + 2]));
	}
	gap /= A.n;
	R = desc && desc->rows ? desc->rows : (unsigned int)fmax(1, floor(gap / (0.5 * (A.h + B.h)) + 0.5));
	report->rows = R;
	report->loop_vertices[0] = A.n;
	report->loop_vertices[1] = B.n;
	report->gap = (qaws_scalar)gap;
	/* the rings: 0 is a's loop, R is b's; between, counts from one to the other, spacing graded from a's edges to b's */
	ring_start = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)R + 2));
	ring_n = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)R + 2));
	t = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)R + 2));
	if (!ring_start || !ring_n || !t)
		goto done;
	{
		double tot = 0;
		t[0] = 0;
		for (r = 0; r < R; r++)
		{
			tot += A.h + (B.h - A.h) * (r + 0.5) / R;
			t[r + 1] = tot;
		}
		for (r = 0; r <= R; r++)
			t[r] /= tot;
	}
	nva = a->hm.nv;
	nvb = same ? 0 : b->hm.nv;
	for (r = 0; r <= R; r++)
	{
		ring_n[r] = r == 0 ? A.n : r == R ? B.n : (unsigned int)floor(A.n + ((double)B.n - A.n) * r / R + 0.5);
		if (r > 0 && r < R)
		{
			ring_start[r] = nva + nvb + nbv;
			nbv += ring_n[r];
		}
	}
	bx = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)nbv + 1));
	if (!bx)
		goto done;
	for (r = 1, k = 0; r < R; r++)
		for (i = 0; i < ring_n[r]; i++, k++)
		{
			/* the Hermite blend at t[r] between the loops' points at share s, leaving a and arriving into b along their surfaces */
			double s = (double)i / ring_n[r], pa[3], oa[3], pb[3], ob[3], L, tt = t[r], h00, h10, h01, h11;
			unsigned int q;
			cn_at(&A, s, pa, oa);
			cn_at(&B, s, pb, ob);
			L = sqrt((pb[0] - pa[0]) * (pb[0] - pa[0]) + (pb[1] - pa[1]) * (pb[1] - pa[1]) + (pb[2] - pa[2]) * (pb[2] - pa[2]));
			h00 = 2 * tt * tt * tt - 3 * tt * tt + 1;
			h10 = tt * tt * tt - 2 * tt * tt + tt;
			h01 = -2 * tt * tt * tt + 3 * tt * tt;
			h11 = tt * tt * tt - tt * tt;
			for (q = 0; q < 3; q++)
				bx[3 * (size_t)k + q] = h00 * pa[q] + h10 * L * oa[q] + h01 * pb[q] - h11 * L * ob[q];
		}
	/* faces: a's, b's, then the bridge's rows */
	nfa = qaws_mesh_3d_get_face_count(a);
	ca = qaws_mesh_3d_get_corner_count(a);
	if (!same)
	{
		nfb = qaws_mesh_3d_get_face_count(b);
		cb = qaws_mesh_3d_get_corner_count(b);
	}
	for (r = 0; r < R; r++)
	{
		int quads = desc && desc->quads && ring_n[r] == ring_n[r + 1];
		nbf += quads ? ring_n[r] : ring_n[r] + ring_n[r + 1];
		nbc += quads ? 4 * ring_n[r] : 3 * (ring_n[r] + ring_n[r + 1]);
	}
	pos = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)nva + nvb + nbv + 1));
	fv = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)ca + cb + nbc + 1));
	fs = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nfa + nfb + nbf + 1));
	afs = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nfa + 1));
	if (!pos || !fv || !fs || !afs)
		goto done;
	if (qaws_mesh_3d_get_positions(pos, a) != QAWS_STATUS_OK || qaws_mesh_3d_get_faces(fv, afs, a) != QAWS_STATUS_OK)
		goto done;
	memcpy(fs, afs, sizeof(unsigned int) * nfa);
	if (!same)
	{
		bfs = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nfb + 1));
		if (!bfs || qaws_mesh_3d_get_positions(pos + 3 * (size_t)nva, b) != QAWS_STATUS_OK || qaws_mesh_3d_get_faces(fv + ca, bfs, b) != QAWS_STATUS_OK)
			goto done;
		for (i = 0; i < cb; i++)
			fv[ca + i] += nva;
		memcpy(fs + nfa, bfs, sizeof(unsigned int) * nfb);
	}
	for (i = 0; i < 3 * nbv; i++)
		pos[3 * ((size_t)nva + nvb) + i] = (qaws_scalar)bx[i];
	{
		unsigned int nf = nfa + nfb, nc = ca + cb;
		for (r = 0; r < R; r++)
		{
			/* ring r's vertices (lower, nearer a) and ring r + 1's (upper), with their shares; lower edges run as a's border does */
			unsigned int nl = ring_n[r], nu = ring_n[r + 1], il = 0, iu = 0;
			unsigned int *lv = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * nl), *uv = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * nu);
			double *ls = (double*)qaws_internal_malloc(sizeof(double) * (nl + 1)), *us = (double*)qaws_internal_malloc(sizeof(double) * (nu + 1));
			int quads = desc && desc->quads && nl == nu;
			if (!lv || !uv || !ls || !us)
			{
				qaws_internal_free(lv); qaws_internal_free(uv); qaws_internal_free(ls); qaws_internal_free(us);
				goto done;
			}
			for (i = 0; i < nl; i++)
			{
				lv[i] = r == 0 ? A.v[i] : ring_start[r] + i;
				ls[i] = r == 0 ? A.s[i] : (double)i / nl;
			}
			for (i = 0; i < nu; i++)
			{
				uv[i] = r + 1 == R ? (same ? B.v[i] : nva + B.v[i]) : ring_start[r + 1] + i;
				us[i] = r + 1 == R ? B.s[i] : (double)i / nu;
			}
			ls[nl] = 1;
			us[nu] = 1;
			if (quads)
				for (i = 0; i < nl; i++)
				{
					fv[nc++] = lv[(i + 1) % nl];
					fv[nc++] = lv[i];
					fv[nc++] = uv[i];
					fv[nc++] = uv[(i + 1) % nu];
					fs[nf++] = 4;
				}
			else
				while (il < nl || iu < nu)
				{
					/* (the ring whose next vertex comes first by arc length advances) */
					double sl = il < nl ? ls[il + 1] : HUGE_VAL, su = iu < nu ? us[iu + 1] : HUGE_VAL;
					if (sl <= su)
					{
						fv[nc++] = lv[(il + 1) % nl];
						fv[nc++] = lv[il % nl];
						fv[nc++] = uv[iu % nu];
						il++;
					}
					else
					{
						fv[nc++] = uv[iu % nu];
						fv[nc++] = uv[(iu + 1) % nu];
						fv[nc++] = lv[il % nl];
						iu++;
					}
					fs[nf++] = 3;
				}
			qaws_internal_free(lv); qaws_internal_free(uv); qaws_internal_free(ls); qaws_internal_free(us);
		}
		report->bridge_faces = nf - nfa - nfb;
		{
			qaws_mesh_3d_desc d;
			memset(&d, 0, sizeof d);
			d.positions = pos;
			d.vertex_count = nva + nvb + nbv;
			d.face_vertices = fv;
			d.face_sizes = fs;
			d.face_count = nf;
			st = qaws_mesh_3d_create(out, &d);
			if (st != QAWS_STATUS_OK)
				goto done;
		}
		{
			/* the faces' parts; the bridge's smallest angle */
			qaws_scalar* part = (qaws_scalar*)qaws_internal_malloc(sizeof(qaws_scalar) * ((size_t)nf + 1));
			qaws_attribute_desc ad;
			unsigned int at = 0, f;
			double mina = HUGE_VAL;
			if (!part)
			{
				st = QAWS_STATUS_ALLOCATION_FAILURE;
				goto done;
			}
			for (f = 0; f < nf; f++)
				part[f] = (qaws_scalar)(f < nfa ? 0 : f < nfa + nfb ? 2 : 1);
			memset(&ad, 0, sizeof ad);
			ad.name = "part";
			ad.element = QAWS_MESH_FACE;
			ad.components = 1;
			ad.rule = QAWS_ATTRIBUTE_NEAREST;
			st = qaws_mesh_3d_add_attribute(&at, *out, part, &ad);
			qaws_internal_free(part);
			if (st != QAWS_STATUS_OK)
				goto done;
			if (desc && (desc->flags & QAWS_CONNECT_FAIR) && nbv)
			{
				/* (a smoothing flow, not a thin plate: spanning two loops far apart, a thin plate pinches like a soap film; the flow
				   evens the bridge out and leaves its shape) */
				unsigned char* held = (unsigned char*)qaws_internal_malloc((size_t)nva + nvb + nbv);
				qaws_smooth_desc sd;
				if (!held)
				{
					st = QAWS_STATUS_ALLOCATION_FAILURE;
					goto done;
				}
				memset(held, 1, (size_t)nva + nvb);
				memset(held + nva + nvb, 0, nbv);
				memset(&sd, 0, sizeof sd);
				sd.method = QAWS_SMOOTH_BILAPLACIAN;
				sd.iterations = 3;
				sd.held = held;
				st = qaws_mesh_3d_smooth(*out, &sd);
				qaws_internal_free(held);
				if (st != QAWS_STATUS_OK)
					goto done;
			}
			for (f = nfa + nfb, k = ca + cb; f < nf; k += fs[f], f++)
				for (j = 0; j < fs[f]; j++)
				{
					qaws_vec3 p = qaws_mesh_3d_get_vertex(*out, fv[k + j]), q = qaws_mesh_3d_get_vertex(*out, fv[k + (j + 1) % fs[f]]);
					qaws_vec3 w = qaws_mesh_3d_get_vertex(*out, fv[k + (j + fs[f] - 1) % fs[f]]);
					double P[3] = { p.x, p.y, p.z }, Q[3] = { q.x, q.y, q.z }, W[3] = { w.x, w.y, w.z }, ang = cn_angle(P, Q, W);
					mina = ang < mina ? ang : mina;
				}
			report->min_angle = (qaws_scalar)mina;
		}
	}
	st = QAWS_STATUS_OK;
done:
	if (st != QAWS_STATUS_OK && out && *out)
	{
		qaws_mesh_3d_destroy(*out);
		*out = NULL;
	}
	cn_loop_free(&A);
	cn_loop_free(&B);
	qaws_internal_free(ring_start); qaws_internal_free(ring_n); qaws_internal_free(t); qaws_internal_free(bx);
	qaws_internal_free(pos); qaws_internal_free(fv); qaws_internal_free(fs); qaws_internal_free(afs);
	qaws_internal_free(bfs);
	return st;
}
