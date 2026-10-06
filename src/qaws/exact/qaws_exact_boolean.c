#include "qaws_exact_curve.h"
#include "../qaws_boolean_2d.h"
#include "../internal/qaws_internal_types.h"
#include "../internal/qaws_internal_curve.h"
#include <math.h>
#include <string.h>

#define TRY(x) do { st = (x); if (st != QAWS_STATUS_OK) return st; } while (0)
#define BO_MAX_CROSSINGS 256

/* Is the curve closed (start and end the same point, exactly)? */
static int is_closed(qaws_exact_curve const* c)
{
	qaws_exact_span const* f = &c->spans[0];
	qaws_exact_span const* l = &c->spans[c->span_count - 1];
	qaws_exact_int x, y;
	unsigned int i;
	for (i = 0; i < 2; i++)
	{
		if (qaws_exact_int_mul(&x, &f->h[i], &l->h[l->degree * 3 + 2]) != QAWS_STATUS_OK ||
		    qaws_exact_int_mul(&y, &l->h[l->degree * 3 + i], &f->h[2]) != QAWS_STATUS_OK || qaws_exact_int_cmp(&x, &y) != 0)
			return 0;
	}
	return 1;
}

/*
 * The homogeneous point of the curve at the double parameter t (any dyadic,
 * not only the parameter lattice): integer De Casteljau at the rational
 * local s = (t 2^shift - a) / (b - a), times a positive factor.
 */
static qaws_status point_at(qaws_exact_curve const* c, double t, qaws_exact_int* P)
{
	qaws_exact_span const* sp = &c->spans[c->span_count - 1];
	qaws_exact_int num, den, d[(QAWS_EXACT_MAX_DEGREE + 1) * 3], omn, t1, t2;
	unsigned int s, n, r, i, k;
	int64_t m;
	int e, sh;
	qaws_status st;
	for (s = 0; s < c->span_count; s++)
		if (t <= ldexp((double)c->spans[s].b, -c->param_shift))
		{
			sp = &c->spans[s];
			break;
		}
	/* t 2^shift = m 2^(e + shift) */
	qaws_exact_split_double(t, &m, &e);
	sh = e + c->param_shift;
	qaws_exact_int_from_i64(&num, m);
	qaws_exact_int_from_i64(&den, sp->b - sp->a);
	qaws_exact_int_from_i64(&t1, sp->a);
	if (m == 0)
		qaws_exact_int_zero(&num);
	if (sh >= 0)
	{
		if (m != 0)
			TRY(qaws_exact_int_shl(&num, &num, (unsigned int)sh));
	}
	else
	{
		TRY(qaws_exact_int_shl(&t1, &t1, (unsigned int)-sh));
		TRY(qaws_exact_int_shl(&den, &den, (unsigned int)-sh));
	}
	TRY(qaws_exact_int_sub(&num, &num, &t1));   /* s = num / den */
	TRY(qaws_exact_int_sub(&omn, &den, &num));
	n = sp->degree;
	for (i = 0; i < (n + 1) * 3; i++)
		d[i] = sp->h[i];
	for (r = 1; r <= n; r++)
		for (i = 0; i + r <= n; i++)
			for (k = 0; k < 3; k++)
			{
				TRY(qaws_exact_int_mul(&t1, &d[i * 3 + k], &omn));
				TRY(qaws_exact_int_mul(&t2, &d[(i + 1) * 3 + k], &num));
				TRY(qaws_exact_int_add(&d[i * 3 + k], &t1, &t2));
			}
	for (k = 0; k < 3; k++)
		P[k] = d[k];
	/* reduce by the common gcd (positive) */
	qaws_exact_int_gcd(&t1, &P[0], &P[1]);
	qaws_exact_int_gcd(&t1, &t1, &P[2]);
	if (qaws_exact_int_bits(&t1) > 1)
		for (k = 0; k < 3; k++)
			TRY(qaws_exact_int_divmod(&P[k], NULL, &P[k], &t1));
	return QAWS_STATUS_OK;
}

typedef struct bo_side
{
	qaws_exact_curve const* c;
	unsigned int n;                        /* crossings */
	unsigned int order[BO_MAX_CROSSINGS];  /* crossing indices sorted along this curve */
	unsigned int pos[BO_MAX_CROSSINGS];    /* crossing -> its rank along this curve */
	int inside[BO_MAX_CROSSINGS];          /* piece j (from order[j] to order[j + 1]) inside the other region */
	int keep[BO_MAX_CROSSINGS];
	int used[BO_MAX_CROSSINGS];
} bo_side;

typedef struct bo_cross
{
	double lo[2], hi[2];   /* parameter enclosures on a (0) and b (1) */
} bo_cross;

/* An interior parameter of piece j of a side, strictly between its end crossings. */
static qaws_status piece_param(bo_side const* sd, bo_cross const* x, unsigned int side, unsigned int j, double* t)
{
	double t_start = ldexp((double)sd->c->spans[0].a, -sd->c->param_shift);
	double t_end = ldexp((double)sd->c->spans[sd->c->span_count - 1].b, -sd->c->param_shift);
	double lo, hi;
	if (sd->n == 0)
	{
		*t = 0.5 * (t_start + t_end);
		return QAWS_STATUS_OK;
	}
	lo = x[sd->order[j]].hi[side];
	hi = x[sd->order[(j + 1) % sd->n]].lo[side];
	if (j + 1 == sd->n)
	{
		/* the piece wraps through the closing point */
		if (lo < t_end)
			hi = t_end;
		else
		{
			lo = t_start;
		}
	}
	*t = 0.5 * (lo + hi);
	if (!(lo < *t && *t < hi))
		return QAWS_STATUS_CERTIFICATION_FAILED;   /* crossings closer than the doubles separate */
	return QAWS_STATUS_OK;
}

static int keep_piece(unsigned int op, unsigned int side, int inside)
{
	if (op == QAWS_BOOLEAN_UNION)
		return !inside;
	if (op == QAWS_BOOLEAN_INTERSECTION)
		return inside;
	return side == 0 ? !inside : inside;   /* a minus b: a outside b, b inside a */
}

static qaws_status push_piece(qaws_exact_piece* out, unsigned int capacity, unsigned int* count, bo_side const* sd, bo_cross const* x, unsigned int side,
	unsigned int j, int reversed)
{
	qaws_exact_piece p;
	double s_lo, s_hi, e_lo, e_hi;
	if (sd->n == 0)
	{
		s_lo = s_hi = ldexp((double)sd->c->spans[0].a, -sd->c->param_shift);
		e_lo = e_hi = ldexp((double)sd->c->spans[sd->c->span_count - 1].b, -sd->c->param_shift);
	}
	else
	{
		bo_cross const* a = &x[sd->order[j]];
		bo_cross const* b = &x[sd->order[(j + 1) % sd->n]];
		s_lo = a->lo[side];
		s_hi = a->hi[side];
		e_lo = b->lo[side];
		e_hi = b->hi[side];
	}
	if (*count >= capacity)
		return QAWS_STATUS_BUFFER_TOO_SMALL;
	p.region = side;
	p.reversed = reversed;
	p.t0_lo = reversed ? e_lo : s_lo;
	p.t0_hi = reversed ? e_hi : s_hi;
	p.t1_lo = reversed ? s_lo : e_lo;
	p.t1_hi = reversed ? s_hi : e_hi;
	out[(*count)++] = p;
	return QAWS_STATUS_OK;
}

qaws_status qaws_exact_boolean_2d(qaws_exact_curve const* a, qaws_exact_curve const* b, unsigned int operation, qaws_exact_piece* out_pieces,
	unsigned int piece_capacity, unsigned int* out_piece_count, qaws_exact_loop* out_loops, unsigned int loop_capacity, unsigned int* out_loop_count)
{
	qaws_exact_pair* pairs;
	bo_cross* x;
	bo_side* sd;
	unsigned int np = 0, i, j, s, npieces = 0, nloops = 0;
	qaws_status st = QAWS_STATUS_OK;
	if (!a || !b || !out_piece_count || !out_loop_count || (!out_pieces && piece_capacity) || (!out_loops && loop_capacity) || a->dimension != 2 ||
	    b->dimension != 2 || operation > QAWS_BOOLEAN_DIFFERENCE)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_piece_count = 0;
	*out_loop_count = 0;
	if (!is_closed(a) || !is_closed(b))
		return QAWS_STATUS_INVALID_ARGUMENT;
	pairs = (qaws_exact_pair*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_exact_pair) * BO_MAX_CROSSINGS));
	x = (bo_cross*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(bo_cross) * BO_MAX_CROSSINGS));
	sd = (bo_side*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(bo_side) * 2));
	if (!pairs || !x || !sd)
		st = QAWS_STATUS_ALLOCATION_FAILURE;
	/* both boundaries simple */
	for (s = 0; s < 2 && st == QAWS_STATUS_OK; s++)
	{
		st = qaws_exact_curve_self_hits(s == 0 ? a : b, pairs, BO_MAX_CROSSINGS, &np);
		if (st == QAWS_STATUS_OK && np > 0)
			st = QAWS_STATUS_CERTIFICATION_FAILED;
	}
	if (st == QAWS_STATUS_OK)
		st = qaws_exact_curve_curve_hits(a, b, pairs, BO_MAX_CROSSINGS, &np);
	if (st == QAWS_STATUS_OK)
	{
		double t_end[2], t_start[2];
		t_start[0] = ldexp((double)a->spans[0].a, -a->param_shift);
		t_end[0] = ldexp((double)a->spans[a->span_count - 1].b, -a->param_shift);
		t_start[1] = ldexp((double)b->spans[0].a, -b->param_shift);
		t_end[1] = ldexp((double)b->spans[b->span_count - 1].b, -b->param_shift);
		for (i = 0; i < np; i++)
		{
			x[i].lo[0] = pairs[i].a_lo;
			x[i].hi[0] = pairs[i].a_hi;
			x[i].lo[1] = pairs[i].b_lo;
			x[i].hi[1] = pairs[i].b_hi;
			/* the closing point is one point: keep it at the start */
			for (s = 0; s < 2; s++)
				if (x[i].lo[s] == t_end[s] && x[i].hi[s] == t_end[s])
					x[i].lo[s] = x[i].hi[s] = t_start[s];
		}
	}
	/* order the crossings along each boundary */
	for (s = 0; s < 2 && st == QAWS_STATUS_OK; s++)
	{
		bo_side* d = &sd[s];
		d->c = s == 0 ? a : b;
		d->n = np;
		for (i = 0; i < np; i++)
			d->order[i] = i;
		for (i = 1; i < np; i++)
		{
			unsigned int key = d->order[i];
			int k = (int)i - 1;
			while (k >= 0 && x[d->order[k]].lo[s] > x[key].lo[s])
			{
				d->order[k + 1] = d->order[k];
				k--;
			}
			d->order[k + 1] = key;
		}
		for (i = 0; i < np; i++)
			d->pos[d->order[i]] = i;
		/* classify each piece by the other boundary's winding number at an exact interior point */
		for (j = 0; j < (np ? np : 1) && st == QAWS_STATUS_OK; j++)
		{
			double t;
			qaws_exact_int P[3];
			qaws_exact_curve const* other = s == 0 ? b : a;
			int w = 0;
			st = piece_param(d, x, s, j, &t);
			if (st == QAWS_STATUS_OK) st = point_at(d->c, t, P);
			if (st == QAWS_STATUS_OK) st = qaws_exact_winding_2d_hom(&other, 1, P, &w);
			d->inside[j] = w != 0;
			d->keep[j] = keep_piece(operation, s, d->inside[j]);
			d->used[j] = 0;
		}
	}
	if (st == QAWS_STATUS_OK && np == 0)
	{
		/* no crossing: whole boundaries */
		for (s = 0; s < 2 && st == QAWS_STATUS_OK; s++)
			if (sd[s].keep[0])
			{
				if (nloops >= loop_capacity)
				{
					st = QAWS_STATUS_BUFFER_TOO_SMALL;
					break;
				}
				out_loops[nloops].first = npieces;
				out_loops[nloops].count = 1;
				nloops++;
				st = push_piece(out_pieces, piece_capacity, &npieces, &sd[s], x, s, 0, 0);
			}
	}
	else if (st == QAWS_STATUS_OK)
	{
		/* link: at each crossing two kept pieces meet; at a touch the same curve goes on */
		for (s = 0; s < 2 && st == QAWS_STATUS_OK; s++)
			for (j = 0; j < np && st == QAWS_STATUS_OK; j++)
			{
				unsigned int cs = s, cj = j, guard = 0, first = npieces;
				int rev = 0;
				if (!sd[s].keep[j] || sd[s].used[j])
					continue;
				while (st == QAWS_STATUS_OK && !sd[cs].used[cj])
				{
					unsigned int node, o = 1 - cs, nxt_same, in_o, out_o;
					int touch;
					sd[cs].used[cj] = 1;
					st = push_piece(out_pieces, piece_capacity, &npieces, &sd[cs], x, cs, cj, rev);
					if (st != QAWS_STATUS_OK || ++guard > 4 * np + 4)
					{
						if (st == QAWS_STATUS_OK)
							st = QAWS_STATUS_INTERNAL_ERROR;
						break;
					}
					/* the crossing where this piece ends */
					node = rev ? sd[cs].order[cj] : sd[cs].order[(cj + 1) % np];
					nxt_same = rev ? (cj + np - 1) % np : (cj + 1) % np;
					out_o = sd[o].pos[node];                /* the other's piece starting here */
					in_o = (out_o + np - 1) % np;           /* and the one ending here */
					touch = sd[cs].inside[cj] == sd[cs].inside[nxt_same] && sd[o].inside[in_o] == sd[o].inside[out_o];
					if (touch && sd[cs].keep[nxt_same])
						cj = nxt_same;   /* continue on the same curve, same direction */
					else if (sd[o].keep[out_o])
					{
						cs = o;
						cj = out_o;
						rev = 0;
					}
					else if (sd[o].keep[in_o])
					{
						cs = o;
						cj = in_o;
						rev = 1;
					}
					else if (sd[cs].keep[nxt_same])
						cj = nxt_same;
					else
					{
						st = QAWS_STATUS_INTERNAL_ERROR;
						break;
					}
				}
				if (st == QAWS_STATUS_OK)
				{
					if (nloops >= loop_capacity)
						st = QAWS_STATUS_BUFFER_TOO_SMALL;
					else
					{
						out_loops[nloops].first = first;
						out_loops[nloops].count = npieces - first;
						nloops++;
					}
				}
			}
	}
	qaws_internal_dealloc(NULL, pairs);
	qaws_internal_dealloc(NULL, x);
	qaws_internal_dealloc(NULL, sd);
	*out_piece_count = npieces;
	*out_loop_count = nloops;
	return st;
}
