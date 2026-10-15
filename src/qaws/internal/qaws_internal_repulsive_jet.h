#ifndef QAWS_INTERNAL_REPULSIVE_JET_H
#define QAWS_INTERNAL_REPULSIVE_JET_H

/* Second-order forward jets for the repulsive energies (qaws_repulsive.c)
   and their hierarchical evaluation (qaws_repulsive_bh.c). */

/* (static inline: unused functions draw no warning) */
#define RP_INLINE static inline
#include <math.h>

#define RP_MAXN 20      /* a term's variables: at most 6 vertices of 3 coordinates, or a far term's 19 */

/* ------------------------------------------------------------------------ */
/*  Second-order jets                                                       */
/* ------------------------------------------------------------------------ */

/* f and its derivatives: g along each variable, t and t2 along the
   direction, h the derivative of g along the direction */
typedef struct rp_jet
{
	double v, t, t2;
	double g[RP_MAXN], h[RP_MAXN];
} rp_jet;

RP_INLINE void rp_const(rp_jet* r, double c, unsigned int n)
{
	unsigned int k;
	r->v = c; r->t = r->t2 = 0;
	for (k = 0; k < n; k++) r->g[k] = r->h[k] = 0;
}

/* variable number slot (n = 0: no gradient wanted), moving at rate dx */
RP_INLINE void rp_var(rp_jet* r, double x, double dx, unsigned int slot, unsigned int n)
{
	rp_const(r, x, n);
	r->t = dx;
	if (slot < n) r->g[slot] = 1;
}

RP_INLINE void rp_add(rp_jet* r, rp_jet const* a, rp_jet const* b, unsigned int n)
{
	unsigned int k;
	r->v = a->v + b->v; r->t = a->t + b->t; r->t2 = a->t2 + b->t2;
	for (k = 0; k < n; k++) { r->g[k] = a->g[k] + b->g[k]; r->h[k] = a->h[k] + b->h[k]; }
}

RP_INLINE void rp_sub(rp_jet* r, rp_jet const* a, rp_jet const* b, unsigned int n)
{
	unsigned int k;
	r->v = a->v - b->v; r->t = a->t - b->t; r->t2 = a->t2 - b->t2;
	for (k = 0; k < n; k++) { r->g[k] = a->g[k] - b->g[k]; r->h[k] = a->h[k] - b->h[k]; }
}

RP_INLINE void rp_scale(rp_jet* r, rp_jet const* a, double c, unsigned int n)
{
	unsigned int k;
	r->v = c * a->v; r->t = c * a->t; r->t2 = c * a->t2;
	for (k = 0; k < n; k++) { r->g[k] = c * a->g[k]; r->h[k] = c * a->h[k]; }
}

RP_INLINE void rp_mul(rp_jet* r, rp_jet const* a, rp_jet const* b, unsigned int n)
{
	rp_jet q;
	unsigned int k;
	q.v = a->v * b->v;
	q.t = a->t * b->v + a->v * b->t;
	q.t2 = a->t2 * b->v + 2 * a->t * b->t + a->v * b->t2;
	for (k = 0; k < n; k++)
	{
		q.g[k] = a->g[k] * b->v + a->v * b->g[k];
		q.h[k] = a->h[k] * b->v + a->g[k] * b->t + a->t * b->g[k] + a->v * b->h[k];
	}
	*r = q;
}

/* u^s for u >= 0. At u = 0 (u a square, so its derivatives vanish there
   too) the limits of the smooth formulas: u itself for s = 1, all zero for
   s > 1 (|w|^(2s) has zero value, gradient and Hessian at w = 0). 0 when
   undefined (u < 0, or u = 0 with s < 1). */
RP_INLINE int rp_pow(rp_jet* r, rp_jet const* u, double s, unsigned int n)
{
	double f0, f1, f2;
	unsigned int k;
	rp_jet q;
	if (u->v <= 0)
	{
		if (u->v < 0 || s < 1)
		{
			rp_const(r, 0, n);   /* (written all the same) */
			return 0;
		}
		if (s == 1)
			*r = *u;
		else
			rp_const(r, 0, n);
		return 1;
	}
	f0 = pow(u->v, s);
	f1 = s * f0 / u->v;
	f2 = (s - 1) * f1 / u->v;
	q.v = f0;
	q.t = f1 * u->t;
	q.t2 = f1 * u->t2 + f2 * u->t * u->t;
	for (k = 0; k < n; k++)
	{
		q.g[k] = f1 * u->g[k];
		q.h[k] = f1 * u->h[k] + f2 * u->g[k] * u->t;
	}
	*r = q;
	return 1;
}

RP_INLINE void rp_dot(rp_jet* r, rp_jet const* a, rp_jet const* b, unsigned int dim, unsigned int n)
{
	rp_jet p;
	unsigned int c;
	rp_mul(r, &a[0], &b[0], n);
	for (c = 1; c < dim; c++)
	{
		rp_mul(&p, &a[c], &b[c], n);
		rp_add(r, r, &p, n);
	}
}

RP_INLINE void rp_cross3(rp_jet* r, rp_jet const* a, rp_jet const* b, unsigned int n)
{
	rp_jet p, q;
	unsigned int c;
	for (c = 0; c < 3; c++)
	{
		rp_mul(&p, &a[(c + 1) % 3], &b[(c + 2) % 3], n);
		rp_mul(&q, &a[(c + 2) % 3], &b[(c + 1) % 3], n);
		rp_sub(&r[c], &p, &q, n);
	}
}

/* |a x b|^2 (2D: the scalar cross product squared) */
RP_INLINE void rp_cross2(rp_jet* r, rp_jet const* a, rp_jet const* b, unsigned int dim, unsigned int n)
{
	rp_jet c[3];
	if (dim == 3)
	{
		rp_cross3(c, a, b, n);
		rp_dot(r, c, c, 3, n);
	}
	else
	{
		rp_jet p, q;
		rp_mul(&p, &a[0], &b[1], n);
		rp_mul(&q, &a[1], &b[0], n);
		rp_sub(&c[0], &p, &q, n);
		rp_mul(r, &c[0], &c[0], n);
	}
}

#endif /* QAWS_INTERNAL_REPULSIVE_JET_H */
