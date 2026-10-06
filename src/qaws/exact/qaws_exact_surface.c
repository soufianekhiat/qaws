#include "qaws_exact_surface.h"
#include "../qaws_surface.h"
#include "../qaws_diff.h"
#include "../internal/qaws_internal_types.h"
#include "../internal/qaws_internal_curve.h"
#include <math.h>
#include <string.h>

#define TRY(x) do { st = (x); if (st != QAWS_STATUS_OK) return st; } while (0)
#define QAWS_EXACT_MAX_NET 64

static void exact_surface_free(qaws_exact_surface* es)
{
	unsigned int i;
	if (!es)
		return;
	for (i = 0; es->patch && i < es->nu * es->nv; i++)
		qaws_internal_dealloc(NULL, es->patch[i]);
	qaws_internal_dealloc(NULL, es->patch);
	qaws_internal_dealloc(NULL, es->ub);
	qaws_internal_dealloc(NULL, es->vb);
	qaws_internal_dealloc(NULL, es);
}

void qaws_exact_surface_destroy(qaws_exact_surface* surface)
{
	exact_surface_free(surface);
}

/* Knots onto the lattice of their largest magnitude; non-empty spans counted. */
static qaws_status lattice_knots(double const* kn, unsigned int nk, unsigned int deg, unsigned int ncp, unsigned int param_bits,
	int64_t* K, int* shift, unsigned int* spans, double* max_err)
{
	double tmax = 0, err;
	unsigned int i;
	qaws_status st;
	for (i = 0; i < nk; i++)
		if (fabs(kn[i]) > tmax)
			tmax = fabs(kn[i]);
	*shift = qaws_exact_param_shift_for(tmax > 0 ? tmax : 1, param_bits);
	for (i = 0; i < nk; i++)
	{
		TRY(qaws_exact_quantize(kn[i], -*shift, param_bits + 1, &K[i], &err));
		if (err > *max_err) *max_err = err;
	}
	*spans = 0;
	for (i = deg; i < ncp; i++)
		if (K[i] < K[i + 1])
			(*spans)++;
	return *spans ? QAWS_STATUS_OK : QAWS_STATUS_EXACT_UNSUPPORTED;
}

/*
 * Patch (su, sv) by two blossom passes: every control column b gives the
 * u-Bezier points of span su (rows i = 0..p), then each row i gives the
 * v-Bezier points of span sv. One lcm clears the whole net.
 */
static qaws_status extract_patch(qaws_exact_hfrac const* net, unsigned int ncv, int64_t const* Ku, int64_t const* Kv, unsigned int su,
	unsigned int sv, unsigned int p, unsigned int q, qaws_exact_int* h)
{
	qaws_exact_hfrac col[QAWS_EXACT_MAX_DEGREE + 1], row[QAWS_EXACT_MAX_DEGREE + 1];
	qaws_exact_hfrac* rowsrc = NULL;
	qaws_exact_hfrac* out = NULL;
	unsigned int a, b, i, j;
	qaws_status st = QAWS_STATUS_OK;
	rowsrc = (qaws_exact_hfrac*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_exact_hfrac) * (p + 1) * (q + 1)));
	out = (qaws_exact_hfrac*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_exact_hfrac) * (p + 1) * (q + 1)));
	if (!rowsrc || !out)
		st = QAWS_STATUS_ALLOCATION_FAILURE;
	/* u pass on the q + 1 columns the v span touches */
	for (b = 0; b <= q && st == QAWS_STATUS_OK; b++)
	{
		for (a = 0; a <= p; a++)
			col[a] = net[(su - p + a) * ncv + (sv - q + b)];
		for (i = 0; i <= p && st == QAWS_STATUS_OK; i++)
			st = qaws_exact_blossom(col, Ku, su, p, i, 4, &rowsrc[i * (q + 1) + b]);
	}
	/* v pass on each row */
	for (i = 0; i <= p && st == QAWS_STATUS_OK; i++)
	{
		for (b = 0; b <= q; b++)
			row[b] = rowsrc[i * (q + 1) + b];
		for (j = 0; j <= q && st == QAWS_STATUS_OK; j++)
			st = qaws_exact_blossom(row, Kv, sv, q, j, 4, &out[i * (q + 1) + j]);
	}
	if (st == QAWS_STATUS_OK)
		st = qaws_exact_clear_denominators(out, (p + 1) * (q + 1) - 1, 4, h);
	qaws_internal_dealloc(NULL, rowsrc);
	qaws_internal_dealloc(NULL, out);
	return st;
}

qaws_status qaws_exact_surface_prepare(qaws_exact_desc const* desc, qaws_surface const* surface, qaws_exact_surface** out_surface,
	qaws_exact_report* out_report)
{
	qaws_exact_desc d;
	qaws_surface_kind kind;
	qaws_field_desc fields[8];
	unsigned int nf = 0, f, ncp = 0, nku = 0, nkv = 0, got = 0, p, q, ncu, ncv, i, c, su, sv, k, spu = 0, spv = 0, storage = 0;
	int rational, wexp = 0;
	double max_pos = 0, max_w = 0, max_k = 0, wmax = 0;
	double kud[QAWS_EXACT_MAX_KNOTS], kvd[QAWS_EXACT_MAX_KNOTS];
	int64_t Ku[QAWS_EXACT_MAX_KNOTS], Kv[QAWS_EXACT_MAX_KNOTS];
	qaws_scalar* cps = NULL;
	qaws_scalar* ws = NULL;
	qaws_scalar kn[QAWS_EXACT_MAX_KNOTS];
	qaws_exact_hfrac* net = NULL;
	qaws_exact_surface* es = NULL;
	qaws_status st = QAWS_STATUS_OK;
	if (!surface || !out_surface)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_surface = NULL;
	if (desc)
		d = *desc;
	else
		qaws_exact_desc_default(&d);
	if (d.param_bits == 0 || d.param_bits > 56 || d.coord_bits == 0 || d.coord_bits > 32 || d.weight_bits == 0 || d.weight_bits > 30)
		return QAWS_STATUS_INVALID_ARGUMENT;
	kind = qaws_surface_get_kind(surface);
	if (kind != QAWS_SURFACE_KIND_BEZIER && kind != QAWS_SURFACE_KIND_BSPLINE && kind != QAWS_SURFACE_KIND_NURBS)
		return QAWS_STATUS_EXACT_UNSUPPORTED;
	rational = kind == QAWS_SURFACE_KIND_NURBS;
	p = qaws_surface_get_u_degree(surface);
	q = qaws_surface_get_v_degree(surface);
	if (qaws_surface_describe_fields(surface, fields, 8, &nf) != QAWS_STATUS_OK)
		return QAWS_STATUS_EXACT_UNSUPPORTED;
	for (f = 0; f < nf && f < 8; f++)
	{
		if (fields[f].field == QAWS_FIELD_CONTROL_POINTS) ncp = fields[f].count;
		if (fields[f].field == QAWS_FIELD_U_KNOTS) nku = fields[f].count;
		if (fields[f].field == QAWS_FIELD_V_KNOTS) nkv = fields[f].count;
	}
	if (kind == QAWS_SURFACE_KIND_BEZIER)
	{
		ncu = p + 1;
		ncv = q + 1;
		nku = 2 * (p + 1);
		nkv = 2 * (q + 1);
	}
	else
	{
		ncu = nku - p - 1;
		ncv = nkv - q - 1;
	}
	if (p < 1 || q < 1 || p > QAWS_EXACT_MAX_DEGREE || q > QAWS_EXACT_MAX_DEGREE || nku > QAWS_EXACT_MAX_KNOTS || nkv > QAWS_EXACT_MAX_KNOTS ||
	    ncu > QAWS_EXACT_MAX_NET || ncv > QAWS_EXACT_MAX_NET || ncu * ncv != ncp)
		return QAWS_STATUS_EXACT_RANGE_EXCEEDED;

	/* knots (a Bezier patch is the clamped single span on [0, 1]) */
	if (kind == QAWS_SURFACE_KIND_BEZIER)
	{
		for (i = 0; i < nku; i++) kud[i] = i <= p ? 0 : 1;
		for (i = 0; i < nkv; i++) kvd[i] = i <= q ? 0 : 1;
	}
	else
	{
		if (qaws_surface_read_field(surface, QAWS_FIELD_U_KNOTS, kn, nku, &got) != QAWS_STATUS_OK)
			return QAWS_STATUS_EXACT_UNSUPPORTED;
		for (i = 0; i < nku; i++) kud[i] = (double)kn[i];
		if (qaws_surface_read_field(surface, QAWS_FIELD_V_KNOTS, kn, nkv, &got) != QAWS_STATUS_OK)
			return QAWS_STATUS_EXACT_UNSUPPORTED;
		for (i = 0; i < nkv; i++) kvd[i] = (double)kn[i];
	}
	es = (qaws_exact_surface*)qaws_internal_alloc(NULL, (unsigned long)sizeof(qaws_exact_surface));
	if (!es)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	memset(es, 0, sizeof(*es));
	es->p = p;
	es->q = q;
	es->space_exp2 = d.space_exp2;
	st = lattice_knots(kud, nku, p, ncu, d.param_bits, Ku, &es->u_shift, &spu, &max_k);
	if (st == QAWS_STATUS_OK)
		st = lattice_knots(kvd, nkv, q, ncv, d.param_bits, Kv, &es->v_shift, &spv, &max_k);

	/* homogeneous integer net */
	cps = (qaws_scalar*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_scalar) * ncp * 3));
	ws = (qaws_scalar*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_scalar) * ncp));
	net = (qaws_exact_hfrac*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_exact_hfrac) * ncp));
	if (st == QAWS_STATUS_OK && (!cps || !ws || !net))
		st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (st == QAWS_STATUS_OK && qaws_surface_read_field(surface, QAWS_FIELD_CONTROL_POINTS, cps, ncp * 3, &got) != QAWS_STATUS_OK)
		st = QAWS_STATUS_EXACT_UNSUPPORTED;
	for (i = 0; st == QAWS_STATUS_OK && i < ncp; i++)
		ws[i] = 1;
	if (st == QAWS_STATUS_OK && rational && qaws_surface_read_field(surface, QAWS_FIELD_WEIGHTS, ws, ncp, &got) != QAWS_STATUS_OK)
		st = QAWS_STATUS_EXACT_UNSUPPORTED;
	for (i = 0; st == QAWS_STATUS_OK && i < ncp; i++)
	{
		if (!(ws[i] > 0))
			st = QAWS_STATUS_EXACT_UNSUPPORTED;
		else if (ws[i] > wmax)
			wmax = ws[i];
	}
	if (st == QAWS_STATUS_OK)
	{
		frexp(wmax, &wexp);
		wexp = (int)d.weight_bits - wexp;
	}
	for (i = 0; st == QAWS_STATUS_OK && i < ncp; i++)
	{
		int64_t W = 1, X;
		double err;
		if (rational)
		{
			double r = nearbyint(ldexp((double)ws[i], wexp));
			if (!(r >= 1))
				st = QAWS_STATUS_EXACT_RANGE_EXCEEDED;
			W = (int64_t)r;
			err = fabs(ldexp(r, -wexp) - (double)ws[i]) / (double)ws[i];
			if (err > max_w) max_w = err;
		}
		for (c = 0; c < 3 && st == QAWS_STATUS_OK; c++)
		{
			st = qaws_exact_quantize((double)cps[i * 3 + c], d.space_exp2, d.coord_bits, &X, &err);
			if (err > max_pos) max_pos = err;
			qaws_exact_int_from_i64(&net[i].num[c], W * X);
		}
		qaws_exact_int_from_i64(&net[i].num[3], W);
		qaws_exact_int_from_i64(&net[i].den, 1);
	}

	/* patches */
	if (st == QAWS_STATUS_OK)
	{
		es->ub = (int64_t*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(int64_t) * (spu + 1)));
		es->vb = (int64_t*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(int64_t) * (spv + 1)));
		es->patch = (qaws_exact_int**)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_exact_int*) * spu * spv));
		if (!es->ub || !es->vb || !es->patch)
			st = QAWS_STATUS_ALLOCATION_FAILURE;
		else
		{
			memset(es->patch, 0, sizeof(qaws_exact_int*) * spu * spv);
			es->nu = spu;
			es->nv = spv;
		}
	}
	if (st == QAWS_STATUS_OK)
	{
		unsigned int iu = 0, iv;
		for (su = p; su < ncu && st == QAWS_STATUS_OK; su++)
		{
			if (!(Ku[su] < Ku[su + 1]))
				continue;
			es->ub[iu] = Ku[su];
			es->ub[iu + 1] = Ku[su + 1];
			iv = 0;
			for (sv = q; sv < ncv && st == QAWS_STATUS_OK; sv++)
			{
				qaws_exact_int* h;
				if (!(Kv[sv] < Kv[sv + 1]))
					continue;
				es->vb[iv] = Kv[sv];
				es->vb[iv + 1] = Kv[sv + 1];
				h = (qaws_exact_int*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_exact_int) * (p + 1) * (q + 1) * 4));
				es->patch[iu * spv + iv] = h;
				if (!h)
				{
					st = QAWS_STATUS_ALLOCATION_FAILURE;
					break;
				}
				st = extract_patch(net, ncv, Ku, Kv, su, sv, p, q, h);
				for (k = 0; k < (p + 1) * (q + 1) * 4 && st == QAWS_STATUS_OK; k++)
					if (qaws_exact_int_bits(&h[k]) > storage)
						storage = qaws_exact_int_bits(&h[k]);
				iv++;
			}
			iu++;
		}
	}
	qaws_internal_dealloc(NULL, cps);
	qaws_internal_dealloc(NULL, ws);
	qaws_internal_dealloc(NULL, net);
	if (st != QAWS_STATUS_OK)
	{
		exact_surface_free(es);
		return st;
	}
	if (out_report)
	{
		memset(out_report, 0, sizeof(*out_report));
		out_report->quality = QAWS_NUMERIC_EXACT_RATIONAL;
		out_report->flags = (max_pos > 0 || max_w > 0 || max_k > 0) ? QAWS_EXACT_FLAG_INPUT_QUANTIZED : QAWS_EXACT_FLAG_NONE;
		out_report->storage_bits = storage;
		out_report->max_position_quantization_error = max_pos;
		out_report->max_weight_quantization_error = max_w;
		out_report->parameter_quantization_error = max_k;
	}
	*out_surface = es;
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Evaluation                                                         */
/* ------------------------------------------------------------------ */

/* Index of the break interval holding X (clamped): lower end inclusive, the last one at the end. */
static unsigned int find_interval(int64_t const* br, unsigned int n, int64_t* X)
{
	unsigned int i;
	if (*X < br[0]) *X = br[0];
	if (*X > br[n]) *X = br[n];
	for (i = 0; i + 1 < n; i++)
		if (*X < br[i + 1])
			return i;
	return n - 1;
}

/*
 * The six homogeneous partials H^(a,b), a + b <= 2, at the local point,
 * each scaled by Lu^p Lv^q (a common factor: it cancels in S and its
 * derivatives). Order: 00, 10, 01, 20, 11, 02.
 */
static qaws_status homogeneous_partials(qaws_exact_surface const* s, int64_t U, int64_t V, qaws_exact_int H[6][4], unsigned int* iu_out,
	unsigned int* iv_out)
{
	static unsigned int const da[6] = { 0, 1, 0, 2, 1, 0 }, db[6] = { 0, 0, 1, 0, 1, 2 };
	qaws_exact_int col[(QAWS_EXACT_MAX_DEGREE + 1) * 4];
	unsigned int iu = find_interval(s->ub, s->nu, &U), iv = find_interval(s->vb, s->nv, &V), k, i;
	qaws_exact_int const* h = s->patch[iu * s->nv + iv];
	qaws_exact_span row, cs;
	qaws_status st;
	row.degree = s->q;
	row.a = s->vb[iv];
	row.b = s->vb[iv + 1];
	cs.degree = s->p;
	cs.a = s->ub[iu];
	cs.b = s->ub[iu + 1];
	cs.h = col;
	for (k = 0; k < 6; k++)
	{
		for (i = 0; i <= s->p; i++)
		{
			row.h = (qaws_exact_int*)&h[i * (s->q + 1) * 4];
			TRY(qaws_exact_homogeneous_derivative(&row, 4, V - row.a, db[k], &col[i * 4]));
		}
		TRY(qaws_exact_homogeneous_derivative(&cs, 4, U - cs.a, da[k], H[k]));
	}
	*iu_out = iu;
	*iv_out = iv;
	return QAWS_STATUS_OK;
}

/* x <- x * y - z * w (scratch t) */
static qaws_status mul_sub(qaws_exact_int* x, qaws_exact_int const* y, qaws_exact_int const* z, qaws_exact_int const* w)
{
	qaws_exact_int t1, t2;
	qaws_status st;
	TRY(qaws_exact_int_mul(&t1, x, y));
	TRY(qaws_exact_int_mul(&t2, z, w));
	return qaws_exact_int_sub(x, &t1, &t2);
}

/*
 * Numerators of the partials over powers of W (S = X / W):
 *   S    = X / W
 *   Su   = Nu / W^2,  Nu = Xu W - X Wu      (Sv alike)
 *   Suu  = (Xuu W^2 - Wuu X W - 2 Wu Nu) / W^3
 *   Suv  = (Xuv W^2 - Wuv X W - Wu Nv - Wv Nu) / W^3
 *   Svv  = (Xvv W^2 - Wvv X W - 2 Wv Nv) / W^3
 * then d/du = (2^u_shift / Lu) d/ds and d/dv alike.
 */
static qaws_status partial_rational(qaws_exact_surface const* s, int64_t U, int64_t V, unsigned int a, unsigned int b,
	qaws_exact_int* num, qaws_exact_int* den, qaws_exact_int* nu_out, qaws_exact_int* nv_out)
{
	qaws_exact_int H[6][4], W2, W3, t1, t2;
	unsigned int iu, iv, c, k;
	qaws_status st;
	if (a + b > 2)
		return QAWS_STATUS_INVALID_ARGUMENT;
	TRY(homogeneous_partials(s, U, V, H, &iu, &iv));
	TRY(qaws_exact_int_mul(&W2, &H[0][3], &H[0][3]));
	TRY(qaws_exact_int_mul(&W3, &W2, &H[0][3]));
	for (c = 0; c < 3; c++)
	{
		qaws_exact_int Nu = H[1][c], Nv = H[2][c];
		TRY(mul_sub(&Nu, &H[0][3], &H[0][c], &H[1][3]));
		TRY(mul_sub(&Nv, &H[0][3], &H[0][c], &H[2][3]));
		if (nu_out) nu_out[c] = Nu;
		if (nv_out) nv_out[c] = Nv;
		if (a + b == 0)
			num[c] = H[0][c];
		else if (a + b == 1)
			num[c] = a ? Nu : Nv;
		else
		{
			/* k: index of the second partial (3: uu, 4: uv, 5: vv) */
			k = a == 2 ? 3 : (a == 1 ? 4 : 5);
			TRY(qaws_exact_int_mul(&num[c], &H[k][c], &W2));
			TRY(qaws_exact_int_mul(&t1, &H[k][3], &H[0][c]));
			TRY(qaws_exact_int_mul(&t1, &t1, &H[0][3]));
			TRY(qaws_exact_int_sub(&num[c], &num[c], &t1));
			if (k == 4)
			{
				TRY(qaws_exact_int_mul(&t1, &H[1][3], &Nv));
				TRY(qaws_exact_int_mul(&t2, &H[2][3], &Nu));
				TRY(qaws_exact_int_add(&t1, &t1, &t2));
			}
			else
			{
				TRY(qaws_exact_int_mul(&t1, &H[k == 3 ? 1 : 2][3], k == 3 ? &Nu : &Nv));
				TRY(qaws_exact_int_mul_i64(&t1, &t1, 2));
			}
			TRY(qaws_exact_int_sub(&num[c], &num[c], &t1));
		}
	}
	*den = a + b == 0 ? H[0][3] : (a + b == 1 ? W2 : W3);
	for (k = 0; k < a; k++)
	{
		for (c = 0; c < 3; c++)
			TRY(qaws_exact_int_shl(&num[c], &num[c], (unsigned int)s->u_shift));
		TRY(qaws_exact_int_mul_i64(den, den, s->ub[iu + 1] - s->ub[iu]));
	}
	for (k = 0; k < b; k++)
	{
		for (c = 0; c < 3; c++)
			TRY(qaws_exact_int_shl(&num[c], &num[c], (unsigned int)s->v_shift));
		TRY(qaws_exact_int_mul_i64(den, den, s->vb[iv + 1] - s->vb[iv]));
	}
	return QAWS_STATUS_OK;
}

qaws_status qaws_exact_surface_eval_rational(qaws_exact_surface const* s, int64_t U, int64_t V, unsigned int a, unsigned int b,
	qaws_exact_int* num, qaws_exact_int* den)
{
	if (!s || !num || !den)
		return QAWS_STATUS_INVALID_ARGUMENT;
	return partial_rational(s, U, V, a, b, num, den, NULL, NULL);
}

qaws_status qaws_exact_surface_normal_rational(qaws_exact_surface const* s, int64_t U, int64_t V, qaws_exact_int* num, qaws_exact_int* den)
{
	qaws_exact_int pos[3], Nu[3], Nv[3], W4, t1, t2;
	unsigned int c, iu, iv;
	qaws_status st;
	if (!s || !num || !den)
		return QAWS_STATUS_INVALID_ARGUMENT;
	TRY(partial_rational(s, U, V, 0, 0, pos, &W4, Nu, Nv));
	/* (Nu / W^2) x (Nv / W^2) */
	for (c = 0; c < 3; c++)
	{
		unsigned int c1 = (c + 1) % 3, c2 = (c + 2) % 3;
		TRY(qaws_exact_int_mul(&t1, &Nu[c1], &Nv[c2]));
		TRY(qaws_exact_int_mul(&t2, &Nu[c2], &Nv[c1]));
		TRY(qaws_exact_int_sub(&num[c], &t1, &t2));
		TRY(qaws_exact_int_shl(&num[c], &num[c], (unsigned int)(s->u_shift + s->v_shift)));
	}
	TRY(qaws_exact_int_mul(&W4, &W4, &W4));
	TRY(qaws_exact_int_mul(&W4, &W4, &W4));
	iu = find_interval(s->ub, s->nu, &U);
	iv = find_interval(s->vb, s->nv, &V);
	TRY(qaws_exact_int_mul_i64(&W4, &W4, s->ub[iu + 1] - s->ub[iu]));
	TRY(qaws_exact_int_mul_i64(den, &W4, s->vb[iv + 1] - s->vb[iv]));
	return QAWS_STATUS_OK;
}

static int64_t lattice_param(double x, int shift)
{
	double r = nearbyint(ldexp(x, shift));
	if (r < -9.0e18) r = -9.0e18;
	if (r > 9.0e18) r = 9.0e18;
	return (int64_t)r;
}

qaws_status qaws_exact_surface_evaluate(qaws_exact_surface const* surface, double u, double v, unsigned int order, double* out,
	double* out_normal)
{
	static unsigned int const da[6] = { 0, 1, 0, 2, 1, 0 }, db[6] = { 0, 0, 1, 0, 1, 2 };
	unsigned int count, k, c;
	int64_t U, V;
	double scale;
	qaws_status st;
	if (!surface || !out || order > 2 || !(u - u == 0.0) || !(v - v == 0.0))
		return QAWS_STATUS_INVALID_ARGUMENT;
	U = lattice_param(u, surface->u_shift);
	V = lattice_param(v, surface->v_shift);
	scale = ldexp(1.0, surface->space_exp2);
	count = order == 0 ? 1 : (order == 1 ? 3 : 6);
	for (k = 0; k < count; k++)
	{
		qaws_exact_int num[3], den;
		TRY(qaws_exact_surface_eval_rational(surface, U, V, da[k], db[k], num, &den));
		for (c = 0; c < 3; c++)
			out[k * 3 + c] = qaws_exact_ratio_to_double(&num[c], &den) * scale;
	}
	if (out_normal)
	{
		qaws_exact_int num[3], den;
		TRY(qaws_exact_surface_normal_rational(surface, U, V, num, &den));
		for (c = 0; c < 3; c++)
			out_normal[c] = qaws_exact_ratio_to_double(&num[c], &den) * scale * scale;
	}
	return QAWS_STATUS_OK;
}

void qaws_exact_surface_patch_count(qaws_exact_surface const* surface, unsigned int* out_u_count, unsigned int* out_v_count)
{
	if (out_u_count) *out_u_count = surface ? surface->nu : 0;
	if (out_v_count) *out_v_count = surface ? surface->nv : 0;
}

qaws_status qaws_exact_surface_patch_bezier(qaws_exact_surface const* surface, unsigned int iu, unsigned int iv, unsigned int* out_p,
	unsigned int* out_q, double out_rect[4], double* out_points, double* out_weights)
{
	qaws_exact_int const* h;
	qaws_exact_int wmax;
	unsigned int n, i, c;
	double scale;
	if (!surface || iu >= surface->nu || iv >= surface->nv)
		return QAWS_STATUS_INVALID_ARGUMENT;
	h = surface->patch[iu * surface->nv + iv];
	n = (surface->p + 1) * (surface->q + 1);
	scale = ldexp(1.0, surface->space_exp2);
	if (out_p) *out_p = surface->p;
	if (out_q) *out_q = surface->q;
	if (out_rect)
	{
		out_rect[0] = ldexp((double)surface->ub[iu], -surface->u_shift);
		out_rect[1] = ldexp((double)surface->ub[iu + 1], -surface->u_shift);
		out_rect[2] = ldexp((double)surface->vb[iv], -surface->v_shift);
		out_rect[3] = ldexp((double)surface->vb[iv + 1], -surface->v_shift);
	}
	wmax = h[3];
	for (i = 1; i < n; i++)
		if (qaws_exact_int_cmp(&h[i * 4 + 3], &wmax) > 0)
			wmax = h[i * 4 + 3];
	for (i = 0; i < n; i++)
	{
		for (c = 0; c < 3 && out_points; c++)
			out_points[i * 3 + c] = qaws_exact_ratio_to_double(&h[i * 4 + c], &h[i * 4 + 3]) * scale;
		if (out_weights)
			out_weights[i] = qaws_exact_ratio_to_double(&h[i * 4 + 3], &wmax);
	}
	return QAWS_STATUS_OK;
}
