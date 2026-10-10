#include "qaws_exact_surface.h"
#include "qaws_exact_solve.h"
#include "../internal/qaws_internal_types.h"
#include "../internal/qaws_internal_curve.h"
#include <math.h>
#include <string.h>

#define TRY(x) do { st = (x); if (st != QAWS_STATUS_OK) return st; } while (0)
#define SSI_MAX_SIZE 625          /* (p1 + 1)(q1 + 1)(p2 + 1)(q2 + 1) */
#define SSI_MAX_SPLITS 40         /* total 4D subdivision depth */
#define SSI_MAX_BOXES 200000
#define SSI_FACE_BOXES 20000
#define SSI_MAX_POINTS 4096

typedef struct ssi_ctx
{
	unsigned int n[4], size, stride[4];
} ssi_ctx;

typedef struct ssi_box
{
	uint64_t lo[4], hi[4];
	int dep[4];
	qaws_exact_int* F;   /* 3 * size */
} ssi_box;

/* A point: local 4D box of the patch pair; the global enclosure. */
typedef struct ssi_pt
{
	qaws_exact_ssi_point g;
} ssi_pt;

static void* ss_alloc(size_t b)
{
	return qaws_internal_alloc(NULL, (unsigned long)b);
}

static void ss_free(void* p)
{
	qaws_internal_dealloc(NULL, p);
}

/*
 * Divides a tensor by the common power of two of its coefficients (a
 * positive factor: signs and roots are kept). Subdivision only ever
 * introduces powers of two (2^n halving, 16^n the 7/16 cut), so this keeps
 * the integers as small as a full gcd would, at the cost of a shift.
 */
static qaws_status ssi_normalize(qaws_exact_int* T, unsigned int size)
{
	unsigned int i, k = ~0u;
	for (i = 0; i < size; i++)
		if (!qaws_exact_int_is_zero(&T[i]))
		{
			unsigned int z = qaws_exact_int_ctz(&T[i]);
			if (z < k)
				k = z;
			if (k == 0)
				return QAWS_STATUS_OK;
		}
	if (k == ~0u)
		return QAWS_STATUS_OK;
	for (i = 0; i < size; i++)
		qaws_exact_int_shr(&T[i], &T[i], k);
	return QAWS_STATUS_OK;
}

/*
 * The cut: at SSI_CUT_A / SSI_CUT_B of a box, not at its middle, so that
 * symmetric or exactly dyadic positions of the curve (1/2, a seam crossing
 * a mid-plane) do not fall on box edges, where no face point certifies.
 */
#define SSI_CUT_A 7
#define SSI_CUT_B 16
#define SSI_CUT_BITS 4

/* Splits T along d at a / b (integer De Casteljau, weights (b - a, a)): both parts times b^n. */
static qaws_status ssi_split_dir(ssi_ctx const* cx, qaws_exact_int const* T, unsigned int d, qaws_exact_int* L, qaws_exact_int* R)
{
	unsigned int n = cx->n[d], sd = cx->stride[d], i, r, o, k;
	qaws_exact_int line[17], t1, t2;
	qaws_status st;
	for (o = 0; o < cx->size; o++)
	{
		if ((o / sd) % (n + 1) != 0)
			continue;
		for (i = 0; i <= n; i++)
			line[i] = T[o + i * sd];
		L[o] = line[0];
		R[o + n * sd] = line[n];
		for (r = 1; r <= n; r++)
		{
			for (i = 0; i + r <= n; i++)
			{
				TRY(qaws_exact_int_mul_i64(&t1, &line[i], SSI_CUT_B - SSI_CUT_A));
				TRY(qaws_exact_int_mul_i64(&t2, &line[i + 1], SSI_CUT_A));
				TRY(qaws_exact_int_add(&line[i], &t1, &t2));
			}
			L[o + r * sd] = line[0];
			R[o + (n - r) * sd] = line[n - r];
		}
		/* point r (left) and n - r (right) came from level r: bring all to b^n */
		for (r = 0; r <= n; r++)
			for (k = r; k < n; k++)
			{
				TRY(qaws_exact_int_mul_i64(&L[o + r * sd], &L[o + r * sd], SSI_CUT_B));
				TRY(qaws_exact_int_mul_i64(&R[o + (n - r) * sd], &R[o + (n - r) * sd], SSI_CUT_B));
			}
	}
	return QAWS_STATUS_OK;
}

static int ssi_excluded(ssi_ctx const* cx, qaws_exact_int const* F)
{
	unsigned int c, i;
	for (c = 0; c < 3; c++)
	{
		int s = qaws_exact_int_sign(&F[c * cx->size]);
		if (s == 0)
			continue;
		for (i = 1; i < cx->size; i++)
			if (qaws_exact_int_sign(&F[c * cx->size + i]) != s)
				break;
		if (i == cx->size)
			return 1;
	}
	return 0;
}

/* ------------------------------------------------------------------ */
/*  Exact interval arithmetic                                          */
/* ------------------------------------------------------------------ */

typedef struct ival
{
	qaws_exact_int lo, hi;
} ival;

static qaws_status iv_mul(ival* r, ival const* a, ival const* b)
{
	qaws_exact_int p[4];
	unsigned int i;
	qaws_status st;
	TRY(qaws_exact_int_mul(&p[0], &a->lo, &b->lo));
	TRY(qaws_exact_int_mul(&p[1], &a->lo, &b->hi));
	TRY(qaws_exact_int_mul(&p[2], &a->hi, &b->lo));
	TRY(qaws_exact_int_mul(&p[3], &a->hi, &b->hi));
	r->lo = p[0];
	r->hi = p[0];
	for (i = 1; i < 4; i++)
	{
		if (qaws_exact_int_cmp(&p[i], &r->lo) < 0) r->lo = p[i];
		if (qaws_exact_int_cmp(&p[i], &r->hi) > 0) r->hi = p[i];
	}
	return QAWS_STATUS_OK;
}

/* r += sign a */
static qaws_status iv_acc(ival* r, ival const* a, int sign)
{
	qaws_status st;
	if (sign > 0)
	{
		TRY(qaws_exact_int_add(&r->lo, &r->lo, &a->lo));
		return qaws_exact_int_add(&r->hi, &r->hi, &a->hi);
	}
	TRY(qaws_exact_int_sub(&r->lo, &r->lo, &a->hi));
	return qaws_exact_int_sub(&r->hi, &r->hi, &a->lo);
}

/*
 * Regular box: some 3 x 3 minor of the interval Jacobian (rows: the
 * equations, columns: the four local variables, bounded by Bernstein
 * derivative coefficients) has an interval determinant not containing 0.
 * Returns the free variable (the ssi_excluded column), or -1.
 */
static qaws_status regular_var(ssi_ctx const* cx, qaws_exact_int const* F, int* out)
{
	ival J[3][4];
	unsigned int c, j, o, free_var;
	qaws_status st;
	*out = -1;
	for (c = 0; c < 3; c++)
		for (j = 0; j < 4; j++)
		{
			int first = 1;
			for (o = 0; o < cx->size; o++)
			{
				qaws_exact_int dlt;
				if ((o / cx->stride[j]) % (cx->n[j] + 1) == cx->n[j])
					continue;
				TRY(qaws_exact_int_sub(&dlt, &F[c * cx->size + o + cx->stride[j]], &F[c * cx->size + o]));
				if (first || qaws_exact_int_cmp(&dlt, &J[c][j].lo) < 0) J[c][j].lo = dlt;
				if (first || qaws_exact_int_cmp(&dlt, &J[c][j].hi) > 0) J[c][j].hi = dlt;
				first = 0;
			}
			/* times the degree: a positive factor per column, the sign test does not need it */
		}
	for (free_var = 0; free_var < 4; free_var++)
	{
		unsigned int col[3], k = 0, p;
		static unsigned int const perm[6][3] = { { 0, 1, 2 }, { 1, 2, 0 }, { 2, 0, 1 }, { 0, 2, 1 }, { 2, 1, 0 }, { 1, 0, 2 } };
		ival det;
		for (j = 0; j < 4; j++)
			if (j != free_var)
				col[k++] = j;
		qaws_exact_int_zero(&det.lo);
		qaws_exact_int_zero(&det.hi);
		for (p = 0; p < 6; p++)
		{
			ival t;
			TRY(iv_mul(&t, &J[0][col[perm[p][0]]], &J[1][col[perm[p][1]]]));
			TRY(iv_mul(&t, &t, &J[2][col[perm[p][2]]]));
			TRY(iv_acc(&det, &t, p < 3 ? 1 : -1));
		}
		if (qaws_exact_int_sign(&det.lo) > 0 || qaws_exact_int_sign(&det.hi) < 0)
		{
			*out = (int)free_var;
			return QAWS_STATUS_OK;
		}
	}
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Points                                                             */
/* ------------------------------------------------------------------ */

static qaws_status to_param(int64_t a, int64_t b, int shift, uint64_t ilo, uint64_t ihi, int dep, double* lo, double* hi)
{
	qaws_exact_span sp;
	int e1, e2;
	qaws_status st;
	sp.degree = 1;
	sp.a = a;
	sp.b = b;
	sp.h = NULL;
	TRY(qaws_exact_span_param_to_double(&sp, shift, ilo, dep, lo, &e1));
	TRY(qaws_exact_span_param_to_double(&sp, shift, ihi, dep, hi, &e2));
	if (!e1) *lo = nextafter(*lo, -HUGE_VAL);
	if (!e2) *hi = nextafter(*hi, HUGE_VAL);
	return QAWS_STATUS_OK;
}

static int ov(double a0, double a1, double b0, double b1)
{
	return a0 <= b1 && b0 <= a1;
}

static int ssi_same_point(qaws_exact_ssi_point const* a, qaws_exact_ssi_point const* b)
{
	return ov(a->u1_lo, a->u1_hi, b->u1_lo, b->u1_hi) && ov(a->v1_lo, a->v1_hi, b->v1_lo, b->v1_hi) && ov(a->u2_lo, a->u2_hi, b->u2_lo, b->u2_hi) &&
	       ov(a->v2_lo, a->v2_hi, b->v2_lo, b->v2_hi);
}

/* The certified points of one face of an unsplit patch-pair box. Adjacent
   patches join exactly, so the face is the same three equations over the
   same region for both boxes that share it: solved once. */
typedef struct ssi_face
{
	unsigned int key[5];      /* direction, boundary index along it, the other three patch indices; key[0] = ~0u: empty */
	qaws_status status;
	unsigned int npts;
	qaws_exact_ssi_point pts[8];
} ssi_face;

typedef struct ssi_out
{
	qaws_exact_ssi_point* pts;
	unsigned int npts;
	ssi_face* faces;           /* open addressing, face_cap a power of two */
	unsigned int face_cap, face_count;
	unsigned int (*arcs)[2];
	unsigned int narcs, arc_cap;
	int seam[4];               /* the parameter (u1, v1, u2, v2) wraps around: its two ends are one place */
	double start[4], end[4];
} ssi_out;

/* Same point, a seam's two ends identified (points exactly on them). */
static int same_point_seam(ssi_out const* o, qaws_exact_ssi_point const* a, qaws_exact_ssi_point const* b)
{
	double const* pa = &a->u1_lo;
	double const* pb = &b->u1_lo;
	unsigned int k;
	for (k = 0; k < 4; k++)
	{
		double alo = pa[2 * k], ahi = pa[2 * k + 1], blo = pb[2 * k], bhi = pb[2 * k + 1];
		if (ov(alo, ahi, blo, bhi))
			continue;
		if (o->seam[k] && ((alo == o->start[k] && ahi == o->start[k] && blo == o->end[k] && bhi == o->end[k]) ||
		                   (alo == o->end[k] && ahi == o->end[k] && blo == o->start[k] && bhi == o->start[k])))
			continue;
		return 0;
	}
	return 1;
}

/* Proportional homogeneous points (exact). */
static int same_hpoint(qaws_exact_int const* P, qaws_exact_int const* Q)
{
	unsigned int c;
	for (c = 0; c < 3; c++)
	{
		qaws_exact_int x, y;
		if (qaws_exact_int_mul(&x, &P[c], &Q[3]) != QAWS_STATUS_OK || qaws_exact_int_mul(&y, &Q[c], &P[3]) != QAWS_STATUS_OK || qaws_exact_int_cmp(&x, &y) != 0)
			return 0;
	}
	return 1;
}

/* Is the surface closed in u (dir 0) or v (dir 1): its first and last patch edges the same, exactly? */
static int surface_closed(qaws_exact_surface const* s, unsigned int dir)
{
	unsigned int p = s->p, q = s->q, i, k;
	if (dir == 0)
	{
		for (i = 0; i < s->nv; i++)
		{
			qaws_exact_int const* f = s->patch[0 * s->nv + i];
			qaws_exact_int const* l = s->patch[(s->nu - 1) * s->nv + i];
			for (k = 0; k <= q; k++)
				if (!same_hpoint(&f[(0 * (q + 1) + k) * 4], &l[(p * (q + 1) + k) * 4]))
					return 0;
		}
		return 1;
	}
	for (i = 0; i < s->nu; i++)
	{
		qaws_exact_int const* f = s->patch[i * s->nv + 0];
		qaws_exact_int const* l = s->patch[i * s->nv + s->nv - 1];
		for (k = 0; k <= p; k++)
			if (!same_hpoint(&f[(k * (q + 1) + 0) * 4], &l[(k * (q + 1) + q) * 4]))
				return 0;
	}
	return 1;
}

/* The index of p among the points, appended when new. */
static qaws_status point_index(ssi_out* o, qaws_exact_ssi_point const* p, unsigned int* idx)
{
	unsigned int i;
	for (i = 0; i < o->npts; i++)
		if (same_point_seam(o, &o->pts[i], p))
		{
			*idx = i;
			return QAWS_STATUS_OK;
		}
	if (o->npts >= SSI_MAX_POINTS)
		return QAWS_STATUS_BUFFER_TOO_SMALL;
	o->pts[o->npts] = *p;
	*idx = o->npts++;
	return QAWS_STATUS_OK;
}

typedef struct ssi_pair
{
	qaws_exact_surface const* a;
	qaws_exact_surface const* b;
	unsigned int iu1, iv1, iu2, iv2;
} ssi_pair;

static unsigned int face_hash(unsigned int const* key)
{
	unsigned int h = 2166136261u, i;
	for (i = 0; i < 5; i++)
		h = (h ^ key[i]) * 16777619u;
	return h;
}

/* the cached face, or the empty slot where it goes (NULL when out of memory) */
static ssi_face* face_slot(ssi_out* o, unsigned int const* key)
{
	unsigned int i, mask;
	if (2 * (o->face_count + 1) > o->face_cap)
	{
		/* grow and rehash */
		unsigned int ncap = o->face_cap ? 2 * o->face_cap : 256, j;
		ssi_face* nf = (ssi_face*)ss_alloc(sizeof(ssi_face) * ncap);
		if (!nf)
			return NULL;
		for (j = 0; j < ncap; j++)
			nf[j].key[0] = ~0u;
		for (j = 0; j < o->face_cap; j++)
			if (o->faces[j].key[0] != ~0u)
			{
				unsigned int m = face_hash(o->faces[j].key) & (ncap - 1);
				while (nf[m].key[0] != ~0u)
					m = (m + 1) & (ncap - 1);
				nf[m] = o->faces[j];
			}
		ss_free(o->faces);
		o->faces = nf;
		o->face_cap = ncap;
	}
	mask = o->face_cap - 1;
	for (i = face_hash(key) & mask;; i = (i + 1) & mask)
	{
		ssi_face* f = &o->faces[i];
		if (f->key[0] == ~0u || memcmp(f->key, key, sizeof(f->key)) == 0)
			return f;
	}
}

/*
 * The certified points where solution arcs cross the faces of a regular
 * box: on each face, three equations in the other three variables.
 */
static qaws_status face_points(ssi_ctx const* cx, ssi_box const* B, ssi_pair const* pr, ssi_out* out, qaws_exact_ssi_point* pts, unsigned int* npts,
	qaws_exact_int* face)
{
	qaws_exact_box3 fb[32];
	unsigned int k, side, i, idx[4];
	int unsplit = B->dep[0] == 0 && B->dep[1] == 0 && B->dep[2] == 0 && B->dep[3] == 0;
	qaws_status st = QAWS_STATUS_OK;
	*npts = 0;
	idx[0] = pr->iu1;
	idx[1] = pr->iv1;
	idx[2] = pr->iu2;
	idx[3] = pr->iv2;
	for (k = 0; k < 4 && st == QAWS_STATUS_OK; k++)
		for (side = 0; side < 2 && st == QAWS_STATUS_OK; side++)
		{
			unsigned int n3[3], m = 0, o, c, fsize = 1, nf = 0, q, key[5], nfp = 0;
			qaws_exact_ssi_point fp[8];
			ssi_face* cached = NULL;
			for (i = 0; i < 4; i++)
				if (i != k)
				{
					n3[m++] = cx->n[i];
					fsize *= cx->n[i] + 1;
				}
			if (unsplit)
			{
				key[0] = k;
				key[1] = idx[k] + side;
				for (i = 0, m = 2; i < 4; i++)
					if (i != k)
						key[m++] = idx[i];
				cached = face_slot(out, key);
				if (!cached)
					return QAWS_STATUS_ALLOCATION_FAILURE;
				if (cached->key[0] != ~0u)
				{
					/* solved by the neighbouring box */
					st = cached->status;
					nfp = cached->npts;
					memcpy(fp, cached->pts, sizeof(qaws_exact_ssi_point) * nfp);
					goto add;
				}
			}
			/* slice: index_k fixed at 0 or n_k, the rest in order */
			for (c = 0; c < 3; c++)
			{
				unsigned int w = 0;
				for (o = 0; o < cx->size; o++)
					if ((o / cx->stride[k]) % (cx->n[k] + 1) == (side ? cx->n[k] : 0))
						face[c * fsize + w++] = B->F[c * cx->size + o];
			}
			st = qaws_exact_solve3(n3, face, fb, 32, &nf, SSI_FACE_BOXES);
			for (q = 0; q < nf && st == QAWS_STATUS_OK; q++)
			{
				/* back to the patch pair's local cube: x = (lo + (hi - lo) y) 2^-dep, y the face box */
				uint64_t lo4[4], hi4[4];
				int dep4[4];
				qaws_exact_ssi_point p;
				for (i = 0; i < 4; i++)
				{
					if (i == k)
					{
						uint64_t v = side ? B->hi[i] : B->lo[i];
						lo4[i] = hi4[i] = v;
						dep4[i] = B->dep[i];
					}
					else
					{
						/* x = lo + (hi - lo) y at the combined depth; coarsen the
						   face box outward when the combined depth would pass 60 bits */
						unsigned int fi = i < k ? i : i - 1;
						int fd = fb[q].dep[fi], ex = B->dep[i] + fd - 60;
						uint64_t flo = fb[q].lo[fi], fhi = fb[q].hi[fi];
						if (ex > 0)
						{
							flo >>= ex;
							fhi = (fhi + (((uint64_t)1 << ex) - 1)) >> ex;
							fd -= ex;
						}
						lo4[i] = (B->lo[i] << fd) + (B->hi[i] - B->lo[i]) * flo;
						hi4[i] = (B->lo[i] << fd) + (B->hi[i] - B->lo[i]) * fhi;
						dep4[i] = B->dep[i] + fd;
					}
				}
				st = to_param(pr->a->ub[pr->iu1], pr->a->ub[pr->iu1 + 1], pr->a->u_shift, lo4[0], hi4[0], dep4[0], &p.u1_lo, &p.u1_hi);
				if (st == QAWS_STATUS_OK) st = to_param(pr->a->vb[pr->iv1], pr->a->vb[pr->iv1 + 1], pr->a->v_shift, lo4[1], hi4[1], dep4[1], &p.v1_lo, &p.v1_hi);
				if (st == QAWS_STATUS_OK) st = to_param(pr->b->ub[pr->iu2], pr->b->ub[pr->iu2 + 1], pr->b->u_shift, lo4[2], hi4[2], dep4[2], &p.u2_lo, &p.u2_hi);
				if (st == QAWS_STATUS_OK) st = to_param(pr->b->vb[pr->iv2], pr->b->vb[pr->iv2 + 1], pr->b->v_shift, lo4[3], hi4[3], dep4[3], &p.v2_lo, &p.v2_hi);
				if (st != QAWS_STATUS_OK)
					break;
				if (nfp >= 8)
					return QAWS_STATUS_OK;   /* too many: the caller subdivides */
				fp[nfp++] = p;
			}
			if (cached && st != QAWS_STATUS_ALLOCATION_FAILURE)
			{
				memcpy(cached->key, key, sizeof(key));
				cached->status = st;
				cached->npts = nfp;
				memcpy(cached->pts, fp, sizeof(qaws_exact_ssi_point) * nfp);
				out->face_count++;
			}
		add:
			/* a point on an edge or corner of the box belongs to several faces: once */
			for (q = 0; q < nfp && st == QAWS_STATUS_OK; q++)
			{
				unsigned int r, dup = 0;
				for (r = 0; r < *npts && !dup; r++)
					dup = ssi_same_point(&pts[r], &fp[q]);
				if (!dup)
				{
					if (*npts >= 8)
						return QAWS_STATUS_OK;   /* too many: the caller subdivides */
					pts[(*npts)++] = fp[q];
				}
			}
		}
	return st;
}

static unsigned int ssi_shallowest(ssi_box const* B)
{
	unsigned int d = 0, k;
	for (k = 1; k < 4; k++)
		if (B->dep[k] < B->dep[d])
			d = k;
	return d;
}

/* One patch pair: subdivide the 4D box; regular boxes give arcs between their face points. */
static qaws_status patch_pair(ssi_pair const* pr, unsigned int min_depth, ssi_out* out, unsigned int* boxes)
{
	ssi_ctx cx;
	qaws_exact_int const* ha = pr->a->patch[pr->iu1 * pr->a->nv + pr->iv1];
	qaws_exact_int const* hb = pr->b->patch[pr->iu2 * pr->b->nv + pr->iv2];
	unsigned int nstack = 0, a1, b1, a2, b2, c, cap = 16;
	ssi_box* stack = NULL;
	qaws_exact_int* store = NULL;
	qaws_exact_int* face = NULL;
	qaws_status st = QAWS_STATUS_OK;
	cx.n[0] = pr->a->p;
	cx.n[1] = pr->a->q;
	cx.n[2] = pr->b->p;
	cx.n[3] = pr->b->q;
	cx.size = (cx.n[0] + 1) * (cx.n[1] + 1) * (cx.n[2] + 1) * (cx.n[3] + 1);
	cx.stride[3] = 1;
	cx.stride[2] = cx.n[3] + 1;
	cx.stride[1] = (cx.n[2] + 1) * (cx.n[3] + 1);
	cx.stride[0] = (cx.n[1] + 1) * cx.stride[1];
	if (cx.size > SSI_MAX_SIZE)
		return QAWS_STATUS_EXACT_UNSUPPORTED;
	stack = (ssi_box*)ss_alloc(sizeof(ssi_box) * cap);
	store = (qaws_exact_int*)ss_alloc(sizeof(qaws_exact_int) * 3 * cx.size * cap);
	face = (qaws_exact_int*)ss_alloc(sizeof(qaws_exact_int) * 3 * cx.size);
	if (!stack || !store || !face)
		st = QAWS_STATUS_ALLOCATION_FAILURE;
	/* F_c = X^a_c W^b - X^b_c W^a on the product Bernstein basis */
	for (c = 0; c < 3 && st == QAWS_STATUS_OK; c++)
	{
		for (a1 = 0; a1 <= cx.n[0] && st == QAWS_STATUS_OK; a1++)
			for (b1 = 0; b1 <= cx.n[1] && st == QAWS_STATUS_OK; b1++)
				for (a2 = 0; a2 <= cx.n[2] && st == QAWS_STATUS_OK; a2++)
					for (b2 = 0; b2 <= cx.n[3] && st == QAWS_STATUS_OK; b2++)
					{
						qaws_exact_int t1, t2;
						unsigned int ia = (a1 * (cx.n[1] + 1) + b1) * 4, ib = (a2 * (cx.n[3] + 1) + b2) * 4;
						unsigned int o = a1 * cx.stride[0] + b1 * cx.stride[1] + a2 * cx.stride[2] + b2;
						st = qaws_exact_int_mul(&t1, &ha[ia + c], &hb[ib + 3]);
						if (st == QAWS_STATUS_OK) st = qaws_exact_int_mul(&t2, &hb[ib + c], &ha[ia + 3]);
						if (st == QAWS_STATUS_OK) st = qaws_exact_int_sub(&store[c * cx.size + o], &t1, &t2);
					}
		if (st == QAWS_STATUS_OK)
			st = ssi_normalize(&store[c * cx.size], cx.size);
	}
	if (st == QAWS_STATUS_OK)
	{
		unsigned int k;
		nstack = 1;
		stack[0].F = store;
		for (k = 0; k < 4; k++)
		{
			stack[0].lo[k] = 0;
			stack[0].hi[k] = 1;
			stack[0].dep[k] = 0;
		}
	}
	while (st == QAWS_STATUS_OK && nstack > 0)
	{
		ssi_box B = stack[--nstack];
		int fv = -1;
		if (++*boxes > SSI_MAX_BOXES)
		{
			st = QAWS_STATUS_CERTIFICATION_FAILED;
			break;
		}
		if (ssi_excluded(&cx, B.F))
			continue;
		st = regular_var(&cx, B.F, &fv);
		if (st != QAWS_STATUS_OK)
			break;
		if (fv >= 0 && (unsigned int)B.dep[0] >= SSI_CUT_BITS * min_depth && (unsigned int)B.dep[1] >= SSI_CUT_BITS * min_depth && (unsigned int)B.dep[2] >= SSI_CUT_BITS * min_depth &&
		    (unsigned int)B.dep[3] >= SSI_CUT_BITS * min_depth)
		{
			qaws_exact_ssi_point pts[8];
			unsigned int np = 0;
			qaws_status fs = face_points(&cx, &B, pr, out, pts, &np, face);
			if (fs == QAWS_STATUS_OK && np <= 2)   /* one point: the curve only touches the box (an arc needs two distinct ends) */
			{
				if (np == 2)
				{
					unsigned int i0 = 0, i1 = 0;
					st = point_index(out, &pts[0], &i0);
					if (st == QAWS_STATUS_OK) st = point_index(out, &pts[1], &i1);
					if (st == QAWS_STATUS_OK)
					{
						if (out->narcs >= out->arc_cap)
							st = QAWS_STATUS_BUFFER_TOO_SMALL;
						else
						{
							out->arcs[out->narcs][0] = i0;
							out->arcs[out->narcs][1] = i1;
							out->narcs++;
						}
					}
				}
				continue;
			}
			if (fs != QAWS_STATUS_OK && fs != QAWS_STATUS_CERTIFICATION_FAILED)
			{
				st = fs;
				break;
			}
			/* an odd count, several arcs, or a face point that does not certify: split */
		}
		if ((unsigned int)(B.dep[0] + B.dep[1] + B.dep[2] + B.dep[3]) >= SSI_CUT_BITS * (SSI_MAX_SPLITS + 4 * min_depth > 56 ? 56 : SSI_MAX_SPLITS + 4 * min_depth))
		{
			st = QAWS_STATUS_CERTIFICATION_FAILED;   /* tangential contact, overlap or a singular point */
			break;
		}
		{
			unsigned int d = ssi_shallowest(&B), slot = nstack, k;
			ssi_box L, R;
			if (slot + 3 > cap)
			{
				/* grow the stack and its tensor slots (box i keeps slot i); B sits in slot nstack */
				unsigned int ncap = cap * 2, q;
				ssi_box* nstk = (ssi_box*)ss_alloc(sizeof(ssi_box) * ncap);
				qaws_exact_int* nsto = (qaws_exact_int*)ss_alloc(sizeof(qaws_exact_int) * 3 * cx.size * ncap);
				if (!nstk || !nsto)
				{
					ss_free(nstk);
					ss_free(nsto);
					st = QAWS_STATUS_ALLOCATION_FAILURE;
					break;
				}
				memcpy(nsto, store, sizeof(qaws_exact_int) * 3 * cx.size * (slot + 1));
				memcpy(nstk, stack, sizeof(ssi_box) * nstack);
				for (q = 0; q < nstack; q++)
					nstk[q].F = nsto + 3 * cx.size * q;
				B.F = nsto + 3 * cx.size * slot;
				ss_free(stack);
				ss_free(store);
				stack = nstk;
				store = nsto;
				cap = ncap;
			}
			L.F = store + 3 * cx.size * (slot + 1);
			R.F = store + 3 * cx.size * (slot + 2);
			for (c = 0; c < 3 && st == QAWS_STATUS_OK; c++)
			{
				st = ssi_split_dir(&cx, &B.F[c * cx.size], d, &L.F[c * cx.size], &R.F[c * cx.size]);
				if (st == QAWS_STATUS_OK) st = ssi_normalize(&L.F[c * cx.size], cx.size);
				if (st == QAWS_STATUS_OK) st = ssi_normalize(&R.F[c * cx.size], cx.size);
			}
			if (st != QAWS_STATUS_OK)
				break;
			for (k = 0; k < 4; k++)
			{
				L.lo[k] = R.lo[k] = B.lo[k];
				L.hi[k] = R.hi[k] = B.hi[k];
				L.dep[k] = R.dep[k] = B.dep[k];
			}
			L.lo[d] = B.lo[d] << SSI_CUT_BITS;
			L.hi[d] = (B.lo[d] << SSI_CUT_BITS) + (B.hi[d] - B.lo[d]) * SSI_CUT_A;
			R.lo[d] = L.hi[d];
			R.hi[d] = B.hi[d] << SSI_CUT_BITS;
			L.dep[d] = R.dep[d] = B.dep[d] + SSI_CUT_BITS;
			memmove(store + 3 * cx.size * slot, R.F, sizeof(qaws_exact_int) * 3 * cx.size);
			R.F = store + 3 * cx.size * slot;
			stack[nstack++] = R;
			stack[nstack++] = L;
		}
	}
	ss_free(stack);
	ss_free(store);
	ss_free(face);
	return st;
}

qaws_status qaws_exact_ssi_solve(qaws_exact_surface const* a, qaws_exact_surface const* b, unsigned int min_depth, unsigned int const* pairs, unsigned int pair_count,
	qaws_exact_ssi_point* out_points,
	unsigned int point_capacity, unsigned int* out_point_count, qaws_exact_ssi_branch* out_branches, unsigned int branch_capacity,
	unsigned int* out_branch_count)
{
	ssi_out o;
	ssi_pair pr;
	unsigned int boxes = 0, i, j, np = 0, nb = 0;
	unsigned int* deg = NULL;
	int* used = NULL;
	qaws_status st = QAWS_STATUS_OK;
	if (!a || !b || !out_point_count || !out_branch_count || (!out_points && point_capacity) || (!out_branches && branch_capacity))
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_point_count = 0;
	*out_branch_count = 0;
	if (a->space_exp2 != b->space_exp2)
		return QAWS_STATUS_EXACT_INCOMPATIBLE_SPACE;
	if (min_depth > 10)
		min_depth = 10;
	memset(&o, 0, sizeof(o));
	/* closed surfaces: a point on the seam has two parameter values */
	o.seam[0] = surface_closed(a, 0);
	o.seam[1] = surface_closed(a, 1);
	o.seam[2] = surface_closed(b, 0);
	o.seam[3] = surface_closed(b, 1);
	o.start[0] = ldexp((double)a->ub[0], -a->u_shift);
	o.end[0] = ldexp((double)a->ub[a->nu], -a->u_shift);
	o.start[1] = ldexp((double)a->vb[0], -a->v_shift);
	o.end[1] = ldexp((double)a->vb[a->nv], -a->v_shift);
	o.start[2] = ldexp((double)b->ub[0], -b->u_shift);
	o.end[2] = ldexp((double)b->ub[b->nu], -b->u_shift);
	o.start[3] = ldexp((double)b->vb[0], -b->v_shift);
	o.end[3] = ldexp((double)b->vb[b->nv], -b->v_shift);
	o.pts = (qaws_exact_ssi_point*)ss_alloc(sizeof(qaws_exact_ssi_point) * SSI_MAX_POINTS);
	o.arc_cap = 2 * SSI_MAX_POINTS;
	o.arcs = (unsigned int (*)[2])ss_alloc(sizeof(unsigned int) * 2 * o.arc_cap);
	if (!o.pts || !o.arcs)
		st = QAWS_STATUS_ALLOCATION_FAILURE;
	pr.a = a;
	pr.b = b;
	/* only the listed patch pairs: the others are proven apart */
	for (i = 0; i < pair_count && st == QAWS_STATUS_OK; i++)
	{
		pr.iu1 = pairs[4 * i];
		pr.iv1 = pairs[4 * i + 1];
		pr.iu2 = pairs[4 * i + 2];
		pr.iv2 = pairs[4 * i + 3];
		st = patch_pair(&pr, min_depth, &o, &boxes);
	}
	/* chain the arcs: points of degree 1 end branches, degree 2 continue them */
	if (st == QAWS_STATUS_OK)
	{
		deg = (unsigned int*)ss_alloc(sizeof(unsigned int) * (o.npts + 1));
		used = (int*)ss_alloc(sizeof(int) * (o.narcs + 1));
		if (!deg || !used)
			st = QAWS_STATUS_ALLOCATION_FAILURE;
	}
	if (st == QAWS_STATUS_OK)
	{
		memset(deg, 0, sizeof(unsigned int) * (o.npts + 1));
		memset(used, 0, sizeof(int) * (o.narcs + 1));
		for (i = 0; i < o.narcs; i++)
		{
			if (o.arcs[i][0] == o.arcs[i][1])
			{
				st = QAWS_STATUS_CERTIFICATION_FAILED;   /* two face points that cannot be told apart */
			}
			deg[o.arcs[i][0]]++;
			deg[o.arcs[i][1]]++;
		}
		for (i = 0; i < o.npts; i++)
			if (deg[i] > 2)
			{
				st = QAWS_STATUS_CERTIFICATION_FAILED;   /* a branching: a singular intersection point */
			}
	}
	/* open branches first (from a degree-1 point), then the closed ones */
	for (j = 0; j < 2 && st == QAWS_STATUS_OK; j++)
		for (i = 0; i < o.npts && st == QAWS_STATUS_OK; i++)
		{
			unsigned int cur = i, first = np, k, closed = 0;
			if (j == 0 ? deg[i] != 1 : deg[i] != 2)
				continue;
			/* a point already placed: skip */
			for (k = 0; k < o.narcs; k++)
				if (!used[k] && (o.arcs[k][0] == i || o.arcs[k][1] == i))
					break;
			if (k == o.narcs)
				continue;
			if (nb >= branch_capacity)
			{
				st = QAWS_STATUS_BUFFER_TOO_SMALL;
				break;
			}
			for (;;)
			{
				if (np >= point_capacity)
				{
					st = QAWS_STATUS_BUFFER_TOO_SMALL;
					break;
				}
				out_points[np++] = o.pts[cur];
				for (k = 0; k < o.narcs; k++)
					if (!used[k] && (o.arcs[k][0] == cur || o.arcs[k][1] == cur))
						break;
				if (k == o.narcs)
					break;
				used[k] = 1;
				cur = o.arcs[k][0] == cur ? o.arcs[k][1] : o.arcs[k][0];
				if (cur == i)
				{
					closed = 1;
					break;
				}
			}
			out_branches[nb].first = first;
			out_branches[nb].count = np - first;
			out_branches[nb].closed = (int)closed;
			nb++;
		}
	ss_free(o.pts);
	ss_free(o.faces);
	ss_free(o.arcs);
	ss_free(deg);
	ss_free(used);
	*out_point_count = np;
	*out_branch_count = nb;
	return st;
}

qaws_status qaws_exact_surface_surface_hits(qaws_exact_surface const* a, qaws_exact_surface const* b, unsigned int min_depth, qaws_exact_ssi_point* out_points,
	unsigned int point_capacity, unsigned int* out_point_count, qaws_exact_ssi_branch* out_branches, unsigned int branch_capacity,
	unsigned int* out_branch_count)
{
	unsigned int na, nb, i, j, k, n = 0;
	unsigned int* pairs;
	double (*box)[2][3];
	qaws_status st;
	if (!a || !b || !out_point_count || !out_branch_count || (!out_points && point_capacity) || (!out_branches && branch_capacity))
		return QAWS_STATUS_INVALID_ARGUMENT;
	na = a->nu * a->nv;
	nb = b->nu * b->nv;
	/* patch pairs whose control boxes overlap */
	pairs = (unsigned int*)ss_alloc(sizeof(unsigned int) * 4 * na * nb + 4);
	box = (double (*)[2][3])ss_alloc(sizeof(double) * 6 * (na + nb));
	if (!pairs || !box)
	{
		ss_free(pairs);
		ss_free(box);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < na; i++)
		qaws_exact_patch_box(a, i / a->nv, i % a->nv, box[i][0], box[i][1]);
	for (j = 0; j < nb; j++)
		qaws_exact_patch_box(b, j / b->nv, j % b->nv, box[na + j][0], box[na + j][1]);
	for (i = 0; i < na; i++)
		for (j = 0; j < nb; j++)
		{
			double const* P = box[i][0];
			double const* Q = box[na + j][0];
			int apart = 0;
			for (k = 0; k < 3; k++)
				if (P[k] <= box[i][1][k] && Q[k] <= box[na + j][1][k] && (box[i][1][k] < Q[k] || box[na + j][1][k] < P[k]))
					apart = 1;
			if (apart)
				continue;
			pairs[4 * n] = i / a->nv;
			pairs[4 * n + 1] = i % a->nv;
			pairs[4 * n + 2] = j / b->nv;
			pairs[4 * n + 3] = j % b->nv;
			n++;
		}
	st = qaws_exact_ssi_solve(a, b, min_depth, pairs, n, out_points, point_capacity, out_point_count, out_branches, branch_capacity, out_branch_count);
	ss_free(pairs);
	ss_free(box);
	return st;
}
