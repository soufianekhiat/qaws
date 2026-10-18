#include "qaws_cross_field.h"
#include "../internal/qaws_internal_ldlt.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define CF_PI 3.14159265358979323846

struct cf_system
{
	qaws_hm const* m;
	double const* nrm;
	unsigned int n, nfree;
	int held, loose;                /* some vertex held to its line; some unknown in a component without one */
	unsigned int ncp, *cp;          /* the connected components, per vertex its own */
	unsigned char* ch;              /* per component: some vertex held */
	double *th, *as, *ar;          /* per corner slot: its edge's rescaled angle at its vertex; its rescaled and real angle */
	double *rr, *we;                /* per corner slot: its edge's transport (cos, sin of n rho); its edge's weight */
	unsigned int *e0, *enb;         /* per vertex its edges (CSR): the neighbour, */
	double *ew, *ez;                /* the weight, the transport from the neighbour (2) */
	double *mass, *u, area;         /* per vertex: the lumped mass; the field (2), held vertices their values */
	unsigned int* id;               /* per vertex: its unknown, HM_NONE when held */
	double *b, *x;                  /* 2 per unknown */
	double nnz_a;                   /* the matrix's lower triangle */
	qaws_internal_ldlt f, g;        /* A_FF + 1e-8 M_FF; the relaxation's A_FF + M_FF / geps^2, factored once for its geps */
	double geps;
	double *pw, *pg;               /* per vertex (NULL: none): a guide's penalty weight, its value (2) */
};

static double cf_dot(double const* a, double const* b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void cf_cross(double const* a, double const* b, double* o)
{
	o[0] = a[1] * b[2] - a[2] * b[1];
	o[1] = a[2] * b[0] - a[0] * b[2];
	o[2] = a[0] * b[1] - a[1] * b[0];
}

static void cf_sub(double const* a, double const* b, double* o)
{
	o[0] = a[0] - b[0];
	o[1] = a[1] - b[1];
	o[2] = a[2] - b[2];
}

static double cf_unit(double* v)
{
	double l = sqrt(cf_dot(v, v));
	if (l > 0)
	{
		v[0] /= l; v[1] /= l; v[2] /= l;
	}
	return l;
}

/* v made a unit tangent of the plane of unit normal n (any one when it vanishes) */
static void cf_tangent(double* v, double const* n)
{
	double d = cf_dot(v, n);
	v[0] -= d * n[0]; v[1] -= d * n[1]; v[2] -= d * n[2];
	if (cf_unit(v) == 0)
	{
		double a[3] = { 1, 0, 0 };
		if (fabs(n[0]) > 0.9)
		{
			a[0] = 0;
			a[1] = 1;
		}
		cf_cross(n, a, v);
		cf_unit(v);
	}
}

void cf_free(cf_system* s)
{
	if (!s)
		return;
	free(s->th); free(s->as); free(s->ar); free(s->rr); free(s->we); free(s->e0); free(s->enb); free(s->ew); free(s->ez);
	free(s->mass); free(s->u); free(s->id); free(s->b); free(s->x); free(s->cp); free(s->ch); free(s->pw); free(s->pg);
	qaws_internal_ldlt_free(&s->f);
	qaws_internal_ldlt_free(&s->g);
	free(s);
}

/* A_FF + sigma M_FF in real 2 x 2 blocks ([a -b; b a] for a + ib), each unknown's two rows at its vertex for the ordering;
   factored into f */
static qaws_status cf_factor(cf_system* s, double sigma, qaws_internal_ldlt* f)
{
	qaws_hm const* m = s->m;
	unsigned int nrow = 2 * s->nfree, v, k, at;
	unsigned int *rp = (unsigned int*)malloc(sizeof(unsigned int) * (nrow + 1)), *cj = NULL;
	double *cx = NULL, *xyz = (double*)malloc(sizeof(double) * 3 * (nrow + 1));
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	memset(f, 0, sizeof *f);
	if (!rp || !xyz)
		goto done;
	for (v = 0, at = 0; v < m->nv; v++)
	{
		unsigned int nf = 0, a = s->id[v];
		if (a == HM_NONE)
			continue;
		for (k = s->e0[v]; k < s->e0[v + 1]; k++)
			nf += s->id[s->enb[k]] != HM_NONE;
		rp[2 * a] = at;
		rp[2 * a + 1] = at + 1 + 2 * nf;
		at += 2 + 4 * nf;
		memcpy(xyz + 6 * (size_t)a, m->x + 3 * (size_t)v, sizeof(double) * 3);
		memcpy(xyz + 6 * (size_t)a + 3, m->x + 3 * (size_t)v, sizeof(double) * 3);
	}
	rp[nrow] = at;
	cj = (unsigned int*)malloc(sizeof(unsigned int) * (at + 1));
	cx = (double*)malloc(sizeof(double) * (at + 1));
	if (!cj || !cx)
		goto done;
	for (v = 0; v < m->nv; v++)
	{
		unsigned int a = s->id[v], p0, p1;
		double d = sigma * s->mass[v];
		if (a == HM_NONE)
			continue;
		for (k = s->e0[v]; k < s->e0[v + 1]; k++)
			d += s->ew[k];
		if (s->pw)
			d += s->pw[v];
		p0 = rp[2 * a];
		p1 = rp[2 * a + 1];
		cj[p0] = 2 * a;
		cx[p0++] = d;
		cj[p1] = 2 * a + 1;
		cx[p1++] = d;
		for (k = s->e0[v]; k < s->e0[v + 1]; k++)
		{
			unsigned int bb = s->id[s->enb[k]];
			double zr = -s->ew[k] * s->ez[2 * (size_t)k], zi = -s->ew[k] * s->ez[2 * (size_t)k + 1];
			if (bb == HM_NONE)
				continue;
			cj[p0] = 2 * bb; cx[p0++] = zr;
			cj[p0] = 2 * bb + 1; cx[p0++] = -zi;
			cj[p1] = 2 * bb; cx[p1++] = zi;
			cj[p1] = 2 * bb + 1; cx[p1++] = zr;
		}
	}
	for (k = 0; f == &s->f && k < nrow; k++)
		for (at = rp[k]; at < rp[k + 1]; at++)
			s->nnz_a += cj[at] <= k;
	st = nrow ? qaws_internal_ldlt_factor(f, nrow, rp, cj, cx, xyz, 0) : QAWS_STATUS_OK;
done:
	free(rp); free(cj); free(cx); free(xyz);
	return st;
}

static double cf_angle(cf_system const* s, unsigned int v, double const* g);

qaws_status cf_create(qaws_hm const* m, double const* nrm, mf_level const* L0, unsigned int n, cf_system** out)
{
	cf_system* s = (cf_system*)calloc(1, sizeof(cf_system));
	unsigned int v, f, c, k, nc = m->nc, nv = m->nv, ne = 0, at;
	double mean = 0;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	*out = NULL;
	if (!s)
		return st;
	s->m = m;
	s->nrm = nrm;
	s->n = n;
	s->th = (double*)calloc((size_t)nc + 1, sizeof(double));
	s->as = (double*)calloc((size_t)nc + 1, sizeof(double));
	s->ar = (double*)calloc((size_t)nc + 1, sizeof(double));
	s->we = (double*)calloc((size_t)nc + 1, sizeof(double));
	s->rr = (double*)calloc(2 * (size_t)nc + 2, sizeof(double));
	s->e0 = (unsigned int*)malloc(sizeof(unsigned int) * (nv + 1));
	s->mass = (double*)calloc((size_t)nv + 1, sizeof(double));
	s->u = (double*)malloc(sizeof(double) * 2 * (nv + 1));
	s->id = (unsigned int*)malloc(sizeof(unsigned int) * (nv + 1));
	if (!s->th || !s->as || !s->ar || !s->we || !s->rr || !s->e0 || !s->mass || !s->u || !s->id)
		goto done;
	/* each vertex's corner angles, rescaled to sum to 2 pi (pi on the border) from its first edge on */
	for (v = 0; v < nv; v++)
	{
		double sum = 0, sc, t = 0;
		if (m->vc[v] == HM_NONE)
			continue;
		for (c = m->vc[v]; c != HM_NONE; c = hm_round(m, v, c))
		{
			double a[3], b[3], x[3];
			cf_sub(m->x + 3 * (size_t)hm_to(m, c), m->x + 3 * (size_t)v, a);
			cf_sub(m->x + 3 * (size_t)m->cv[hm_prev(m, c)], m->x + 3 * (size_t)v, b);
			cf_cross(a, b, x);
			s->ar[c] = atan2(sqrt(cf_dot(x, x)), cf_dot(a, b));
			sum += s->ar[c];
		}
		sc = sum > 0 ? (hm_is_border_vertex(m, v) ? CF_PI : 2 * CF_PI) / sum : 1;
		for (c = m->vc[v]; c != HM_NONE; c = hm_round(m, v, c))
		{
			s->th[c] = t;
			s->as[c] = sc * s->ar[c];
			t += s->as[c];
		}
	}
	/* weights: half the cotangent opposite the edge in each triangle (in a polygon its centroid's distance to the edge's
	   midpoint over the edge's length, the same for regular ones); masses: a share of each face's area per corner.
	   shortcut: polygons by that dual length, not a polygon Laplacian (Alexa, Wardetzky 2011); upgrade when SMOOTHEST
	   serves quad meshes that are far from regular */
	for (f = 0; f < m->nf; f++)
	{
		unsigned int c0 = m->fc[f], fn = m->fn[f];
		double nn[3], cen[3] = { 0, 0, 0 }, A;
		hm_face_normal(m, f, nn);
		A = sqrt(cf_dot(nn, nn));
		s->area += A;
		for (c = c0; c < c0 + fn; c++)
		{
			s->mass[m->cv[c]] += A / fn;
			for (k = 0; k < 3; k++)
				cen[k] += m->x[3 * (size_t)m->cv[c] + k] / fn;
		}
		for (c = c0; c < c0 + fn; c++)
		{
			double const *p = m->x + 3 * (size_t)m->cv[c], *q = m->x + 3 * (size_t)hm_to(m, c);
			double a[3], b[3], x[3], l;
			if (fn == 3)
			{
				double const* o = m->x + 3 * (size_t)m->cv[hm_prev(m, c)];
				cf_sub(p, o, a);
				cf_sub(q, o, b);
				cf_cross(a, b, x);
				l = sqrt(cf_dot(x, x));
				s->we[c] += l > 0 ? 0.5 * cf_dot(a, b) / l : 0;
			}
			else
			{
				for (k = 0; k < 3; k++)
				{
					a[k] = cen[k] - 0.5 * (p[k] + q[k]);
					b[k] = q[k] - p[k];
				}
				l = sqrt(cf_dot(b, b));
				s->we[c] += l > 0 ? sqrt(cf_dot(a, a)) / l : 0;
			}
		}
	}
	/* both sides of an edge summed, the weights kept above a millionth of their mean; the transports */
	for (f = 0; f < m->nf; f++)
		for (c = m->fc[f]; c < m->fc[f] + m->fn[f]; c++)
		{
			unsigned int o = m->co[c], nx = hm_next(m, c);
			double rho = s->th[nx] + s->as[nx] - s->th[c] + CF_PI;
			if (o != HM_NONE && o < c)
				continue;
			s->rr[2 * (size_t)c] = cos(n * rho);
			s->rr[2 * (size_t)c + 1] = sin(n * rho);
			if (o != HM_NONE)
			{
				s->we[c] = s->we[o] = s->we[c] + s->we[o];
				s->rr[2 * (size_t)o] = s->rr[2 * (size_t)c];
				s->rr[2 * (size_t)o + 1] = -s->rr[2 * (size_t)c + 1];
			}
			mean += fabs(s->we[c]);
			ne++;
		}
	mean = ne ? 1e-6 * mean / ne : 0;
	for (f = 0; f < m->nf; f++)
		for (c = m->fc[f]; c < m->fc[f] + m->fn[f]; c++)
			s->we[c] = s->we[c] > mean ? s->we[c] : mean;
	/* each vertex's edges: out along its corners (the transport back the conjugate), and on the border in along the last */
	for (v = 0, at = 0; v < nv; v++)
		at += m->vc[v] == HM_NONE ? 0 : hm_valence(m, v);
	s->enb = (unsigned int*)malloc(sizeof(unsigned int) * (at + 1));
	s->ew = (double*)malloc(sizeof(double) * (at + 1));
	s->ez = (double*)malloc(sizeof(double) * 2 * (at + 1));
	if (!s->enb || !s->ew || !s->ez)
		goto done;
	for (v = 0, at = 0; v < nv; v++)
	{
		unsigned int last = HM_NONE;
		s->e0[v] = at;
		if (m->vc[v] == HM_NONE)
			continue;
		for (c = m->vc[v]; c != HM_NONE; c = hm_round(m, v, c))
		{
			s->enb[at] = hm_to(m, c);
			s->ew[at] = s->we[c];
			s->ez[2 * (size_t)at] = s->rr[2 * (size_t)c];
			s->ez[2 * (size_t)at + 1] = -s->rr[2 * (size_t)c + 1];
			at++;
			last = c;
		}
		if (hm_is_border_vertex(m, v))
		{
			unsigned int p = hm_prev(m, last);
			s->enb[at] = m->cv[p];
			s->ew[at] = s->we[p];
			s->ez[2 * (size_t)at] = s->rr[2 * (size_t)p];
			s->ez[2 * (size_t)at + 1] = s->rr[2 * (size_t)p + 1];
			at++;
		}
	}
	s->e0[nv] = at;
	/* the connected components (id the stack, set below): each solved by its lines when it has some, else for its
	   eigenvector */
	s->cp = (unsigned int*)malloc(sizeof(unsigned int) * (nv + 1));
	if (!s->cp)
		goto done;
	for (v = 0; v < nv; v++)
		s->cp[v] = HM_NONE;
	for (v = 0; v < nv; v++)
	{
		unsigned int top = 0;
		if (s->cp[v] != HM_NONE)
			continue;
		s->cp[v] = s->ncp;
		s->id[top++] = v;
		while (top)
			for (c = s->id[--top], k = s->e0[c]; k < s->e0[c + 1]; k++)
				if (s->cp[s->enb[k]] == HM_NONE)
				{
					s->cp[s->enb[k]] = s->ncp;
					s->id[top++] = s->enb[k];
				}
		s->ncp++;
	}
	s->ch = (unsigned char*)calloc((size_t)s->ncp + 1, 1);
	if (!s->ch)
		goto done;
	/* held: off the corners, along the vertex's two feature (or border) edges, the mean of their directions' n-th powers
	   (alike on a straight line: n pi is a whole turn for even n); left free where the two disagree by more than 45 degrees */
	for (v = 0; v < nv; v++)
	{
		s->u[2 * (size_t)v] = 1;
		s->u[2 * (size_t)v + 1] = 0;
		if (m->vc[v] == HM_NONE)
		{
			s->id[v] = HM_NONE;
			continue;
		}
		if ((L0->con[v] & MF_ALONG) && !(L0->con[v] & MF_PINNED))
		{
			double z[2][2];
			unsigned int cnt = 0, last = HM_NONE;
			for (c = m->vc[v]; c != HM_NONE; c = hm_round(m, v, c))
			{
				if (((m->cflag[c] & HM_FEATURE) || m->co[c] == HM_NONE) && cnt++ < 2)
				{
					z[cnt - 1][0] = cos(n * s->th[c]);
					z[cnt - 1][1] = sin(n * s->th[c]);
				}
				last = c;
			}
			if (hm_is_border_vertex(m, v) && cnt++ < 2)
			{
				z[cnt - 1][0] = cos(n * (s->th[last] + s->as[last]));
				z[cnt - 1][1] = sin(n * (s->th[last] + s->as[last]));
			}
			if (cnt == 2 && z[0][0] * z[1][0] + z[0][1] * z[1][1] >= cos(CF_PI / 4))
			{
				double r = hypot(z[0][0] + z[1][0], z[0][1] + z[1][1]);
				s->u[2 * (size_t)v] = (z[0][0] + z[1][0]) / r;
				s->u[2 * (size_t)v + 1] = (z[0][1] + z[1][1]) / r;
				s->id[v] = HM_NONE;
				s->held = s->ch[s->cp[v]] = 1;
				continue;
			}
			if (cnt == 0)
			{
				/* held off the lines: a guide's direction */
				double phi = cf_angle(s, v, L0->cq + 3 * (size_t)v);
				s->u[2 * (size_t)v] = cos(n * phi);
				s->u[2 * (size_t)v + 1] = sin(n * phi);
				s->id[v] = HM_NONE;
				s->held = s->ch[s->cp[v]] = 1;
				continue;
			}
		}
		s->id[v] = s->nfree++;
	}
	for (v = 0; v < nv; v++)
		s->loose |= s->id[v] != HM_NONE && !s->ch[s->cp[v]];
	s->b = (double*)malloc(sizeof(double) * (2 * (size_t)s->nfree + 1));
	s->x = (double*)malloc(sizeof(double) * (2 * (size_t)s->nfree + 1));
	if (s->b && s->x)
		st = cf_factor(s, 1e-8, &s->f);
done:
	if (st != QAWS_STATUS_OK)
	{
		cf_free(s);
		return st;
	}
	*out = s;
	return QAWS_STATUS_OK;
}

void cf_fill(cf_system const* s, double* out_factor, double* out_matrix)
{
	*out_factor = s->f.nnz + s->f.n;
	*out_matrix = s->nnz_a;
}

/* the free unknowns from M q - A_FC u_C (q: 2 per vertex, NULL: none) by f (A_FF + sigma M_FF); held_only: only the
   components held by lines take theirs (the matrix is block diagonal across components) */
static void cf_solve(cf_system* s, qaws_internal_ldlt const* f, double const* q, int held_only)
{
	qaws_hm const* m = s->m;
	unsigned int v, k;
	for (v = 0; v < m->nv; v++)
	{
		unsigned int a = s->id[v];
		double br = 0, bi = 0;
		if (a == HM_NONE)
			continue;
		if (q)
		{
			br = s->mass[v] * q[2 * (size_t)v];
			bi = s->mass[v] * q[2 * (size_t)v + 1];
		}
		if (s->pw)
		{
			br += s->pw[v] * s->pg[2 * (size_t)v];
			bi += s->pw[v] * s->pg[2 * (size_t)v + 1];
		}
		for (k = s->e0[v]; k < s->e0[v + 1]; k++)
		{
			unsigned int j = s->enb[k];
			double zr = s->ew[k] * s->ez[2 * (size_t)k], zi = s->ew[k] * s->ez[2 * (size_t)k + 1];
			if (s->id[j] != HM_NONE)
				continue;
			br += zr * s->u[2 * (size_t)j] - zi * s->u[2 * (size_t)j + 1];
			bi += zr * s->u[2 * (size_t)j + 1] + zi * s->u[2 * (size_t)j];
		}
		s->b[2 * (size_t)a] = br;
		s->b[2 * (size_t)a + 1] = bi;
	}
	qaws_internal_ldlt_solve(f, s->b, s->x);
	for (v = 0; v < m->nv; v++)
		if (s->id[v] != HM_NONE && (!held_only || s->ch[s->cp[v]]))
		{
			s->u[2 * (size_t)v] = s->x[2 * (size_t)s->id[v]];
			s->u[2 * (size_t)v + 1] = s->x[2 * (size_t)s->id[v] + 1];
		}
}

/* back to 3D: the direction psi = arg(u) / n found in the wedge between consecutive edges that holds it, its share of
   the wedge's rescaled angle turned from the wedge's first edge within the wedge's face by that share of its real angle,
   then turned with the wedge's plane into the tangent plane */
static void cf_directions(cf_system const* s, double* out_q)
{
	qaws_hm const* m = s->m;
	unsigned int v, c, k;
	for (v = 0; v < m->nv; v++)
	{
		double const *nv = s->nrm + 3 * (size_t)v, *p = m->x + 3 * (size_t)v;
		double* d = out_q + 3 * (size_t)v;
		double psi = atan2(s->u[2 * (size_t)v + 1], s->u[2 * (size_t)v]) / s->n, t = 1, e1[3], e2[3], dd;
		unsigned int w = HM_NONE;
		d[0] = 1; d[1] = 0; d[2] = 0;
		if (m->vc[v] == HM_NONE)
		{
			cf_tangent(d, nv);
			continue;
		}
		if (psi < 0)
			psi += 2 * CF_PI / s->n;
		for (c = m->vc[v]; c != HM_NONE; c = hm_round(m, v, c))
		{
			w = c;
			if (psi < s->th[c] + s->as[c])
			{
				t = s->as[c] > 0 ? (psi - s->th[c]) / s->as[c] : 0;
				break;
			}
		}
		t = t < 0 ? 0 : t > 1 ? 1 : t;
		cf_sub(m->x + 3 * (size_t)hm_to(m, w), p, e1);
		cf_sub(m->x + 3 * (size_t)m->cv[hm_prev(m, w)], p, e2);
		cf_unit(e1);
		dd = cf_dot(e2, e1);
		for (k = 0; k < 3; k++)
			e2[k] -= dd * e1[k];
		cf_unit(e2);
		for (k = 0; k < 3; k++)
			d[k] = cos(t * s->ar[w]) * e1[k] + sin(t * s->ar[w]) * e2[k];
		{
			/* turned from the wedge's plane into the tangent plane, not projected: on a crease the normal is one side's,
			   and a direction across the crease on the other side is that normal (projected, it would be noise) */
			double f[3], a[3], ad[3], cs, da;
			cf_cross(e1, e2, f);
			cf_cross(f, nv, a);
			cs = cf_dot(f, nv);
			if (cs > -0.5)
			{
				cf_cross(a, d, ad);
				da = cf_dot(a, d) / (1 + cs);
				for (k = 0; k < 3; k++)
					d[k] = cs * d[k] + ad[k] + da * a[k];
			}
		}
		cf_tangent(d, nv);
	}
}

qaws_status cf_smoothest(cf_system* s, double* out_q)
{
	qaws_hm const* m = s->m;
	unsigned int v, c, it;
	double *q = (double*)calloc(2 * (size_t)m->nv + 2, sizeof(double)), *nn = (double*)malloc(sizeof(double) * ((size_t)s->ncp + 1));
	if (!q || !nn)
	{
		free(q);
		free(nn);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	/* a component held by lines: one solve, A_FF u_F = -A_FC u_C (q 0 there); one without: inverse power iteration toward
	   its smallest eigenvector, u = A^-1 M u normalized in M, 20 steps or until a step moves u by less than 1e-6 in M. The
	   matrix is block diagonal across components: each solve serves both */
	for (it = 0; s->nfree && it < 20; it++)
	{
		double moved = 0;
		for (v = 0; v < m->nv; v++)
			if (s->id[v] != HM_NONE && !s->ch[s->cp[v]])
			{
				q[2 * (size_t)v] = it ? s->u[2 * (size_t)v] : 1;
				q[2 * (size_t)v + 1] = it ? s->u[2 * (size_t)v + 1] : 0;
			}
		cf_solve(s, &s->f, q, 0);
		if (!s->loose)
			break;
		memset(nn, 0, sizeof(double) * s->ncp);
		for (v = 0; v < m->nv; v++)
			if (s->id[v] != HM_NONE && !s->ch[s->cp[v]])
				nn[s->cp[v]] += s->mass[v] * (s->u[2 * (size_t)v] * s->u[2 * (size_t)v] + s->u[2 * (size_t)v + 1] * s->u[2 * (size_t)v + 1]);
		for (c = 0; c < s->ncp; c++)
			nn[c] = nn[c] > 0 ? 1 / sqrt(nn[c]) : 1;
		for (v = 0; v < m->nv; v++)
			if (s->id[v] != HM_NONE && !s->ch[s->cp[v]])
			{
				double a, b;
				s->u[2 * (size_t)v] *= nn[s->cp[v]];
				s->u[2 * (size_t)v + 1] *= nn[s->cp[v]];
				a = s->u[2 * (size_t)v] - q[2 * (size_t)v];
				b = s->u[2 * (size_t)v + 1] - q[2 * (size_t)v + 1];
				moved += s->mass[v] * (a * a + b * b);
			}
		if (it && moved < 1e-12)
			break;
	}
	free(q);
	free(nn);
	cf_directions(s, out_q);
	return QAWS_STATUS_OK;
}

/* the rescaled angle of tangent g at v: its angle in the tangent plane placed in the wedge of projected edges that holds it */
static double cf_angle(cf_system const* s, unsigned int v, double const* g)
{
	qaws_hm const* m = s->m;
	double const *nv = s->nrm + 3 * (size_t)v, *p = m->x + 3 * (size_t)v;
	double e[3], f[3], gam;
	unsigned int c;
	cf_sub(m->x + 3 * (size_t)hm_to(m, m->vc[v]), p, e);
	cf_tangent(e, nv);
	cf_cross(nv, e, f);
	gam = atan2(cf_dot(g, f), cf_dot(g, e));
	for (c = m->vc[v]; c != HM_NONE; c = hm_round(m, v, c))
	{
		double a[3], b[3], ba, bb, rel, del;
		cf_sub(m->x + 3 * (size_t)hm_to(m, c), p, a);
		cf_sub(m->x + 3 * (size_t)m->cv[hm_prev(m, c)], p, b);
		ba = atan2(cf_dot(a, f), cf_dot(a, e));
		bb = atan2(cf_dot(b, f), cf_dot(b, e));
		del = fmod(bb - ba + 4 * CF_PI, 2 * CF_PI);
		rel = fmod(gam - ba + 4 * CF_PI, 2 * CF_PI);
		if (rel <= del && del > 0)
			return s->th[c] + s->as[c] * rel / del;
	}
	return gam;
}

qaws_status cf_aligned(cf_system* s, double const* guide, double* out_q)
{
	qaws_hm const* m = s->m;
	double *q = (double*)malloc(sizeof(double) * 2 * (m->nv + 1)), *gc = (double*)calloc((size_t)s->ncp + 1, sizeof(double)), nn = 0;
	unsigned int v;
	int bare = 0;
	qaws_status st;
	if (!q || !gc)
	{
		free(q);
		free(gc);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	/* q = |g|^2 exp(i n phi), normalized to 1 / sqrt(area) in M (the pull then weighs as the held lines at any scale) */
	for (v = 0; v < m->nv; v++)
	{
		double const* g = guide + 3 * (size_t)v;
		double g2 = cf_dot(g, g), phi;
		q[2 * (size_t)v] = q[2 * (size_t)v + 1] = 0;
		if (m->vc[v] == HM_NONE || g2 == 0)
			continue;
		phi = cf_angle(s, v, g);
		q[2 * (size_t)v] = g2 * cos(s->n * phi);
		q[2 * (size_t)v + 1] = g2 * sin(s->n * phi);
		nn += s->mass[v] * g2 * g2;
		gc[s->cp[v]] += g2;
	}
	/* a component with neither lines nor guide would come out 0: drawn to its smoothest field instead */
	for (v = 0; v < m->nv; v++)
		bare |= s->id[v] != HM_NONE && !s->ch[s->cp[v]] && !(gc[s->cp[v]] > 0);
	if (!(nn > 0) || !(s->area > 0) || bare)
	{
		st = cf_smoothest(s, out_q);
		for (v = 0; st == QAWS_STATUS_OK && v < m->nv; v++)
			if (!s->ch[s->cp[v]] && !(gc[s->cp[v]] > 0))
			{
				q[2 * (size_t)v] = s->u[2 * (size_t)v];
				q[2 * (size_t)v + 1] = s->u[2 * (size_t)v + 1];
			}
		if (st != QAWS_STATUS_OK || !(nn > 0) || !(s->area > 0))
		{
			free(q);
			free(gc);
			return st;
		}
	}
	free(gc);
	nn = 1 / sqrt(nn * s->area);
	for (v = 0; v < 2 * m->nv; v++)
		q[v] *= nn;
	if (s->nfree)
		cf_solve(s, &s->f, q, 0);
	free(q);
	cf_directions(s, out_q);
	return QAWS_STATUS_OK;
}

qaws_status cf_guide(cf_system* s, double const* g)
{
	qaws_hm const* m = s->m;
	unsigned int v, k, any = 0;
	if (!s->pw)
	{
		s->pw = (double*)calloc((size_t)m->nv + 1, sizeof(double));
		s->pg = (double*)calloc(2 * (size_t)m->nv + 2, sizeof(double));
		if (!s->pw || !s->pg)
			return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	/* w (the vertex's edge weights) |u - exp(i n phi)|^2: as strong as w times the neighbours together, as LOCAL's pull */
	for (v = 0; v < m->nv; v++)
	{
		double w = sqrt(cf_dot(g + 3 * (size_t)v, g + 3 * (size_t)v)), d = 0, phi;
		s->pw[v] = 0;
		if (s->id[v] == HM_NONE || !(w > 0))
			continue;
		for (k = s->e0[v]; k < s->e0[v + 1]; k++)
			d += s->ew[k];
		phi = cf_angle(s, v, g + 3 * (size_t)v);
		s->pw[v] = w * d;
		s->pg[2 * (size_t)v] = cos(s->n * phi);
		s->pg[2 * (size_t)v + 1] = sin(s->n * phi);
		s->held = s->ch[s->cp[v]] = 1;   /* (its component then solved once, not for its eigenvector) */
		any = 1;
	}
	if (!any)
		return QAWS_STATUS_OK;
	s->loose = 0;
	for (v = 0; v < m->nv; v++)
		s->loose |= s->id[v] != HM_NONE && !s->ch[s->cp[v]];
	qaws_internal_ldlt_free(&s->f);
	qaws_internal_ldlt_free(&s->g);
	s->geps = 0;
	s->nnz_a = 0;
	return s->nfree ? cf_factor(s, 1e-8, &s->f) : QAWS_STATUS_OK;
}

/* Ginzburg-Landau by MBO steps (Merriman, Bence, Osher 1994; Viertel, Osting 2019 for cross fields): the last field
   normalized and diffused for a time eps^2 ((A + M / eps^2) v = M u / eps^2, the held values kept), at most steps
   times, until no normalized value moves by 1e-3. The relaxed field's coinciding zeros (a disc held by its border: one
   zero of order n at its centre, z^n) part into simple ones some eps apart. Components without held vertices are left
   as they are (the eigenvector's zeros are simple). The factor is kept for the next call at the same eps */
qaws_status cf_relax(cf_system* s, double eps, unsigned int steps, double* out_q)
{
	qaws_hm const* m = s->m;
	double *q, sigma = 1 / (eps * eps);
	unsigned int v, it;
	qaws_status st = QAWS_STATUS_OK;
	if (s->nfree && s->held && eps > 0 && steps)
	{
		q = (double*)calloc(2 * (size_t)m->nv + 2, sizeof(double));
		if (!q)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		if (s->geps != eps)
		{
			qaws_internal_ldlt_free(&s->g);
			s->geps = 0;
			st = cf_factor(s, sigma, &s->g);
			s->geps = st == QAWS_STATUS_OK ? eps : 0;
		}
		for (it = 0; st == QAWS_STATUS_OK && it < steps; it++)
		{
			double moved = 0;
			for (v = 0; v < m->nv; v++)
			{
				double l = hypot(s->u[2 * (size_t)v], s->u[2 * (size_t)v + 1]), a = q[2 * (size_t)v], b = q[2 * (size_t)v + 1];
				int on = s->id[v] != HM_NONE && s->ch[s->cp[v]] && l > 0;
				q[2 * (size_t)v] = on ? sigma * s->u[2 * (size_t)v] / l : 0;
				q[2 * (size_t)v + 1] = on ? sigma * s->u[2 * (size_t)v + 1] / l : 0;
				moved = fmax(moved, hypot(q[2 * (size_t)v] - a, q[2 * (size_t)v + 1] - b) / sigma);
			}
			if (moved < 1e-3)
				break;
			cf_solve(s, &s->g, q, 1);
		}
		free(q);
	}
	if (st == QAWS_STATUS_OK)
		cf_directions(s, out_q);
	return st;
}

/* per face: the field's turns along its edges against the transport, with the transport's own turn round it, n times
   the face's rescaled angle excess, taken whole (wrapped into a half turn it loses whole turns on coarse faces; whole,
   the faces' excesses sum to 2 pi chi) */
void cf_index(cf_system const* s, int* out)
{
	qaws_hm const* m = s->m;
	unsigned int f, c;
	for (f = 0; f < m->nf; f++)
	{
		double sum = -(double)s->n * (m->fn[f] - 2.0) * CF_PI;
		for (c = m->fc[f]; c < m->fc[f] + m->fn[f]; c++)
		{
			double const *ui = s->u + 2 * (size_t)m->cv[c], *uj = s->u + 2 * (size_t)hm_to(m, c), *r = s->rr + 2 * (size_t)c;
			double ar = r[0] * ui[0] - r[1] * ui[1], ai = r[0] * ui[1] + r[1] * ui[0];
			sum += atan2(uj[1] * ar - uj[0] * ai, uj[0] * ar + uj[1] * ai) + s->n * s->as[c];
		}
		out[f] = (int)floor(sum / (2 * CF_PI) + 0.5);
	}
}

void cf_restrict(mf_hier* h)
{
	unsigned int l, i, c, k;
	for (l = 1; l < h->nl; l++)
	{
		mf_level *L = &h->lv[l], *F = &h->lv[l - 1];
		for (i = 0; i < L->n; i++)
		{
			double *q = L->q + 3 * (size_t)i, w;
			double const* n = L->nrm + 3 * (size_t)i;
			if (L->con[i] & MF_ALONG)
			{
				memcpy(q, L->cq + 3 * (size_t)i, sizeof(double) * 3);
				continue;
			}
			c = L->child0[i];
			memcpy(q, F->q + 3 * (size_t)L->child[c], sizeof(double) * 3);
			cf_tangent(q, n);
			w = F->area[L->child[c]];
			for (c++; c < L->child0[i + 1]; c++)
			{
				unsigned int ch = L->child[c];
				double a[3], b[3], wc = F->area[ch];
				mf_compat_orientation(q, n, F->q + 3 * (size_t)ch, F->nrm + 3 * (size_t)ch, a, b);
				for (k = 0; k < 3; k++)
					q[k] = a[k] * w + b[k] * wc;
				w += wc;
				cf_tangent(q, n);
			}
		}
	}
}

qaws_status cf_direction_field(qaws_mesh_3d const* mesh, mf_hier const* h, unsigned int n, int aligned, qaws_scalar* out_directions,
	int* out_face_index)
{
	qaws_hm const* m = &mesh->hm;
	cf_system* sys = NULL;
	qaws_mesh_reference* r = NULL;
	double *q = NULL, *g = NULL;
	unsigned int v, k;
	qaws_status st = cf_create(m, h->lv[0].nrm, &h->lv[0], n, &sys);
	if (st != QAWS_STATUS_OK)
		return st;
	q = (double*)malloc(sizeof(double) * 3 * (m->nv + 1));
	if (!q)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	for (v = 0; v < m->nv && !(h->lv[0].con[v] & MF_PREFER); v++) {}
	if (v < m->nv)
	{
		/* the guidance's pulls (mf_guide on level 0) as penalties */
		for (v = 0; v < 3 * m->nv; v++)
			q[v] = (h->lv[0].con[v / 3] & MF_PREFER) ? h->lv[0].cq[v] : 0;
		st = cf_guide(sys, q);
		if (st != QAWS_STATUS_OK)
			goto done;
	}
	if (aligned)
	{
		/* the guide: the direction of largest curvature on the PN reference, as strong as the surface bends one way more */
		g = (double*)calloc(3 * (size_t)m->nv + 3, sizeof(double));
		if (!g)
		{
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			goto done;
		}
		st = qaws_mesh_reference_create(mesh, NULL, &r);
		if (st != QAWS_STATUS_OK)
			goto done;
		for (v = 0; v < m->nv; v++)
		{
			unsigned int l = r->vloc[v];
			double xa[3], xb[3], nn[3], na[3], nb[3], k1, k2, d1[3], d2[3], an;
			if (l == HM_NONE)
				continue;
			mi_eval(&r->patch[l / 3], l % 3 == 1, l % 3 == 2, NULL, xa, xb, NULL, NULL, NULL, nn, na, nb);
			mi_curvature(xa, xb, nn, na, nb, &k1, &k2, d1, d2);
			an = (fabs(k1) + fabs(k2)) * r->extent > 1e-3 ? (k1 - k2) / (fabs(k1) + fabs(k2)) : 0;
			for (k = 0; k < 3; k++)
				g[3 * (size_t)v + k] = d1[k] * an;
		}
		st = cf_aligned(sys, g, q);
	}
	else
		st = cf_smoothest(sys, q);
	if (st != QAWS_STATUS_OK)
		goto done;
	for (v = 0; out_directions && v < 3 * m->nv; v++)
		out_directions[v] = (qaws_scalar)q[v];
	if (out_face_index)
		cf_index(sys, out_face_index);
done:
	free(q);
	free(g);
	qaws_mesh_reference_destroy(r);
	cf_free(sys);
	return st;
}
