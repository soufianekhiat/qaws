/*
 * Spatial queries on meshes (qaws_mesh_bvh.h): a binary BVH over simplices,
 * built by binned SAH, with dipole moments per node for the winding number.
 * Children are allocated when their parent is split, so a parent's index is
 * below its children's: refit runs over the nodes backwards.
 */

#include "../internal/qaws_internal_memory.h"
#include "../qaws_mesh_bvh.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <stddef.h>
#include "../simd/qaws_simd.h"
#include "../internal/qaws_internal_cpu.h"
#if defined(_MSC_VER)
#include <intrin.h>
#endif

#if defined(__SSE2__) || defined(_M_X64) || defined(_M_AMD64) || defined(__AVX2__)
#define BV_SSE 1
#else
#define BV_SSE 0
#endif

/* AVX2 kernels, compiled per function (the library keeps its baseline flags) and chosen at creation by CPUID */
#if BV_SSE && (defined(_MSC_VER) || defined(__GNUC__) || defined(__clang__))
#include <immintrin.h>
#define BV_AVX2 1
#if defined(_MSC_VER) && !defined(__clang__)
#define BV_AVX2_FN
#else
#define BV_AVX2_FN __attribute__((target("avx2,fma")))
#endif
#else
#define BV_AVX2 0
#endif

#define BV_BINS 8
#define BV_CHUNK 64
#define BV_PI 3.14159265358979323846
#define BV_WMAX 8

typedef struct bv_node
{
	double lo[3], hi[3];
	unsigned int a, b;   /* interior: the children; leaf: the first entry of ord and the count */
	unsigned int leaf, first, cnt, lf;   /* (first, cnt: its entries, contiguous in ord; lf: a leaf in the wide tree) */
} bv_node;

/* Queries run on w-wide nodes collapsed from the binary SAH tree (4-wide for every query; on the AVX2 path an 8-wide copy
   with the same leaves for the rays): per slot a float box
   relative to the scene's centre, rounded outward and widened (so the floats' rounding of rays and points never loses a
   hit), tested w at a time. A node is 8 w floats: lo per axis (w each), hi per axis, then w children and w counts
   (inner slot: a node and 0; leaf slot: its first entry of ord and its entries; empty: QAWS_MESH_NONE and 0).
   Empty slots have lo +inf and hi -inf. */

typedef struct bv_wide
{
	float* nf;                  /* the nodes, nst floats each, in mem aligned to 64 bytes */
	void* mem;
	unsigned int n4, w, nst, depth;
	double e;                   /* the expected node visits of a ray (the nodes' areas over the root's) */
} bv_wide;

typedef struct bv_stk
{
	unsigned int node;
	float t;                    /* the node's entry distance (rays) or squared distance (points) */
} bv_stk;

/* per wide node: centre (3), area vector sum (3), first moment C_ij = sum a_i d_j (9), radius, area, second moment T_ijk = sum a_i d_j d_k (27) */
#define BV_WM 44
/* per leaf slot, in floats: centre (relative to the scene's), N, C, radius, T symmetric (18) */
#define BV_LM 34

struct qaws_mesh_bvh
{
	unsigned int dim, k, nv, ns, nn, depth, leaf;
	unsigned int simd;    /* the path: 0 scalar, 1 SSE2, 2 AVX2 */
	double* x;            /* 3 per vertex (z 0 in 2D) */
	unsigned int* s;      /* k per simplex */
	unsigned int* el;     /* per simplex */
	unsigned int* ord;    /* the simplices in leaf order (entries) */
	unsigned int* ee;     /* per entry its element */
	unsigned int* ev;     /* per entry its k vertices */
	double* px;           /* per entry its k vertices' positions, contiguous */
	float* scr;           /* triangles in 3D, for the rays' screen: per entry (ns + 8) 12 floats, v0 - c, e1, e2 each padded by its L1 norm */
	bv_node* nodes;       /* the binary tree (build only) */
	bv_wide t;            /* 4-wide: closest points, winding numbers, rays when r is empty */
	bv_wide r;            /* 8-wide, the same leaves: rays on the AVX2 path */
	double c[3], pad;     /* the scene's centre, the boxes' widening */
	double beta;
	double* wm;           /* BV_WM per node of t when the winding number applies */
	float* lm;            /* BV_LM per slot of t (leaf slots only filled) */
};

/* node i: its slots' lower bounds along axis a, upper bounds, children, counts */
static float* bv_nlo(bv_wide const* t, unsigned int i, unsigned int a)
{
	return t->nf + (size_t)i * t->nst + a * t->w;
}

static float* bv_nhi(bv_wide const* t, unsigned int i, unsigned int a)
{
	return t->nf + (size_t)i * t->nst + (3 + a) * t->w;
}

static unsigned int* bv_nch(bv_wide const* t, unsigned int i)
{
	return (unsigned int*)(t->nf + (size_t)i * t->nst + 6 * t->w);
}

static unsigned int* bv_ncnt(bv_wide const* t, unsigned int i)
{
	return (unsigned int*)(t->nf + (size_t)i * t->nst + 7 * t->w);
}

/* the lowest set bit of a nonzero mask */
static unsigned int bv_ctz(unsigned int m)
{
#if defined(_MSC_VER)
	unsigned long r;
	_BitScanForward(&r, m);
	return (unsigned int)r;
#elif defined(__GNUC__) || defined(__clang__)
	return (unsigned int)__builtin_ctz(m);
#else
	unsigned int r = 0;
	while (!(m >> r & 1))
		r++;
	return r;
#endif
}

static int bv_winds(qaws_mesh_bvh const* b)
{
	return (b->k == 3 && b->dim == 3) || (b->k == 2 && b->dim == 2);
}

static double bv_dot(double const* a, double const* b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void bv_sub(double const* a, double const* b, double* o)
{
	o[0] = a[0] - b[0];
	o[1] = a[1] - b[1];
	o[2] = a[2] - b[2];
}

static void bv_cross(double const* a, double const* b, double* o)
{
	o[0] = a[1] * b[2] - a[2] * b[1];
	o[1] = a[2] * b[0] - a[0] * b[2];
	o[2] = a[0] * b[1] - a[1] * b[0];
}

static double bv_area(double const* lo, double const* hi, unsigned int dim)
{
	double dx = hi[0] - lo[0], dy = hi[1] - lo[1], dz = hi[2] - lo[2];
	if (dx < 0)
		return 0;
	return dim == 2 ? dx + dy : dx * dy + dy * dz + dz * dx;
}

static void bv_grow(double* lo, double* hi, double const* plo, double const* phi)
{
	unsigned int i;
	for (i = 0; i < 3; i++)
	{
		if (plo[i] < lo[i])
			lo[i] = plo[i];
		if (phi[i] > hi[i])
			hi[i] = phi[i];
	}
}

static void bv_empty(double* lo, double* hi)
{
	lo[0] = lo[1] = lo[2] = HUGE_VAL;
	hi[0] = hi[1] = hi[2] = -HUGE_VAL;
}

static void bv_prim_box(qaws_mesh_bvh const* b, unsigned int p, double* lo, double* hi)
{
	unsigned int j;
	bv_empty(lo, hi);
	for (j = 0; j < b->k; j++)
	{
		double const* v = b->x + 3 * (size_t)b->s[(size_t)b->k * p + j];
		bv_grow(lo, hi, v, v);
	}
}

/* ---- build ---- */

/* the float next to f toward +inf (up) or -inf, by its bits */
static float bv_step(float f, int up)
{
	unsigned int u;
	if (f == 0)
		return up ? 1.4e-45f : -1.4e-45f;
	memcpy(&u, &f, sizeof u);
	u = (f > 0) == (up != 0) ? u + 1 : u - 1;
	memcpy(&f, &u, sizeof u);
	return f;
}

/* doubles to floats rounded outward */
static float bv_down(double v)
{
	float f = (float)v;
	return (double)f > v ? bv_step(f, 0) : f;
}

static float bv_up(double v)
{
	float f = (float)v;
	return (double)f < v ? bv_step(f, 1) : f;
}

/* a simplex while building: its box in floats (rounded outward; the tree only groups, refit takes the exact boxes) and its index,
   moved with the partitions so every pass reads in order; lo and hi each 4 floats for SSE (lane 3 unused) */
typedef struct bv_ref
{
	float lo[3];
	unsigned int id;
	float hi[3];
	unsigned int pad;
} bv_ref;

/* a node to split: its refs and its bounds (box, and centroids doubled: lo + hi) */
typedef struct bv_task
{
	unsigned int node, begin, end, depth;
	float lo[4], hi[4], clo[4], chi[4];
} bv_task;

static void bv_ref_swap(bv_ref* r, unsigned int i, unsigned int j)
{
	bv_ref t = r[i];
	r[i] = r[j];
	r[j] = t;
}

static int bv_bin(bv_ref const* r, unsigned int axis, float clo, float scale)
{
	int bin = (int)((r->lo[axis] + r->hi[axis] - clo) * scale);
	return bin < 0 ? 0 : bin >= BV_BINS ? BV_BINS - 1 : bin;
}

/* refs [begin, end) whose bin along axis is at most `last` first; returns the first of the rest */
static unsigned int bv_partition(bv_ref* r, unsigned int axis, float clo, float scale, int last, unsigned int begin, unsigned int end)
{
	unsigned int i = begin, j = end;
	while (i < j)
	{
		if (bv_bin(r + i, axis, clo, scale) <= last)
			i++;
		else
			bv_ref_swap(r, i, --j);
	}
	return i;
}

/* refs [begin, end) arranged so the middle one has the median centroid along axis (quickselect) */
static void bv_median(bv_ref* r, unsigned int axis, unsigned int begin, unsigned int end)
{
	unsigned int lo = begin, hi = end - 1, mid = begin + (end - begin) / 2;
	while (lo < hi)
	{
		float pv = r[lo + (hi - lo) / 2].lo[axis] + r[lo + (hi - lo) / 2].hi[axis];
		unsigned int i = lo, j = hi;
		while (i <= j)
		{
			while (r[i].lo[axis] + r[i].hi[axis] < pv)
				i++;
			while (r[j].lo[axis] + r[j].hi[axis] > pv)
				j--;
			if (i <= j)
			{
				bv_ref_swap(r, i, j);
				i++;
				if (j == 0)
					break;
				j--;
			}
		}
		if (mid <= j)
			hi = j;
		else if (mid >= i)
			lo = i;
		else
			break;
	}
}

static float bv_farea(float const* lo, float const* hi, unsigned int dim)
{
	float dx = hi[0] - lo[0], dy = hi[1] - lo[1], dz = hi[2] - lo[2];
	if (dx < 0)
		return 0;
	return dim == 2 ? dx + dy : dx * dy + dy * dz + dz * dx;
}

/* a 4-float min / max accumulation; the same arithmetic with or without SSE */
static void bv_acc(float* lo, float* hi, float const* plo, float const* phi)
{
#if BV_SSE
	_mm_storeu_ps(lo, _mm_min_ps(_mm_loadu_ps(lo), _mm_loadu_ps(plo)));
	_mm_storeu_ps(hi, _mm_max_ps(_mm_loadu_ps(hi), _mm_loadu_ps(phi)));
#else
	unsigned int i;
	for (i = 0; i < 4; i++)
	{
		lo[i] = lo[i] < plo[i] ? lo[i] : plo[i];
		hi[i] = hi[i] > phi[i] ? hi[i] : phi[i];
	}
#endif
}

static void bv_fempty4(float* lo, float* hi)
{
	lo[0] = lo[1] = lo[2] = lo[3] = HUGE_VALF;
	hi[0] = hi[1] = hi[2] = hi[3] = -HUGE_VALF;
}

/* bounds of refs [begin, end) into t (the median's children) */
static void bv_task_bounds(bv_ref const* r, bv_task* t)
{
	unsigned int i;
	bv_fempty4(t->lo, t->hi);
	bv_fempty4(t->clo, t->chi);
	for (i = t->begin; i < t->end; i++)
	{
		float c[4] = { r[i].lo[0] + r[i].hi[0], r[i].lo[1] + r[i].hi[1], r[i].lo[2] + r[i].hi[2], 0 };
		bv_acc(t->lo, t->hi, r[i].lo, r[i].hi);
		bv_acc(t->clo, t->chi, c, c);
	}
}

#if BV_AVX2
/* refs [begin, end) into the bins of every axis with scale > 0, one 8-float max per ref and axis on (-lo, id, hi, pad):
   the same bins as bv_acc's (min and max are exact), lane 3 of each bound left 0 */
BV_AVX2_FN static void bv_bin_avx2(bv_ref const* r, unsigned int begin, unsigned int end, unsigned int dim, float const* clo, float const* scale,
	float (*blo)[BV_BINS][4], float (*bhi)[BV_BINS][4], unsigned int (*bcount)[BV_BINS])
{
	__m256 acc[3][BV_BINS], neg = _mm256_setr_ps(-0.0f, -0.0f, -0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f), none = _mm256_set1_ps(-HUGE_VALF);
	unsigned int i, a, j, c;
	for (a = 0; a < 3; a++)
		for (j = 0; j < BV_BINS; j++)
			acc[a][j] = none;
	/* the bins of the 3 axes at once, as bv_bin: (lo + hi - clo) scale truncated, clamped to [0, BV_BINS - 1] */
	__m128 c4 = _mm_setr_ps(clo[0], clo[1], clo[2], 0), s4 = _mm_setr_ps(scale[0], scale[1], scale[2], 0);
	__m128i top = _mm_set1_epi32(BV_BINS - 1), zero = _mm_setzero_si128();
	for (i = begin; i < end; i++)
	{
		__m256 v = _mm256_loadu_ps((float const*)(r + i));
		__m128 lo = _mm256_castps256_ps128(v), hi = _mm256_extractf128_ps(v, 1);
		__m128i bi = _mm_cvttps_epi32(_mm_mul_ps(_mm_sub_ps(_mm_add_ps(lo, hi), c4), s4));
		int bins[4];
		bi = _mm_min_epi32(_mm_max_epi32(bi, zero), top);
		_mm_storeu_si128((__m128i*)bins, bi);
		v = _mm256_xor_ps(v, neg);
		for (a = 0; a < dim; a++)
			if (scale[a] > 0)
			{
				bcount[a][bins[a]]++;
				acc[a][bins[a]] = _mm256_max_ps(acc[a][bins[a]], v);
			}
	}
	for (a = 0; a < 3; a++)
		for (j = 0; j < BV_BINS; j++)
		{
			float o[8];
			_mm256_storeu_ps(o, acc[a][j]);
			for (c = 0; c < 3; c++)
			{
				blo[a][j][c] = -o[c];
				bhi[a][j][c] = o[4 + c];
			}
			blo[a][j][3] = bhi[a][j][3] = 0;
		}
	_mm256_zeroupper();
}
#endif

static qaws_status bv_build(qaws_mesh_bvh* b, unsigned int leaf_size)
{
	bv_ref* r = (bv_ref*)qaws_internal_scratch_malloc(qaws_internal_mul2((size_t)b->ns + 1, sizeof(bv_ref)));
	unsigned int top = 0, scap = 64, i, cap = b->ns ? 2 * b->ns - 1 : 0;
	bv_task* stack = (bv_task*)qaws_internal_malloc(scap * sizeof(bv_task));
	b->ord = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2((size_t)b->ns + 1, sizeof(unsigned int)));
	b->nodes = (bv_node*)qaws_internal_malloc(qaws_internal_mul2((size_t)cap + 1, sizeof(bv_node)));
	if (!r || !stack || !b->ord || !b->nodes)
	{
		qaws_internal_free(r);
		qaws_internal_free(stack);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < b->ns; i++)
	{
		double lo[3], hi[3];
		unsigned int a;
		bv_prim_box(b, i, lo, hi);
		for (a = 0; a < 3; a++)
		{
			r[i].lo[a] = bv_down(lo[a]);
			r[i].hi[a] = bv_up(hi[a]);
		}
		r[i].id = i;
		r[i].pad = 0;
	}
	b->nn = 0;
	b->depth = 0;
	if (b->ns)
	{
		b->nn = 1;
		stack[0].node = 0;
		stack[0].begin = 0;
		stack[0].end = b->ns;
		stack[0].depth = 1;
		bv_task_bounds(r, stack);
		top = 1;
	}
	while (top)
	{
		bv_task t = stack[--top];
		bv_node* nd = b->nodes + t.node;
		unsigned int n = t.end - t.begin, axis, best_axis = 3, mid, j, k;
		int best_bin = 0;
		float best = HUGE_VALF, scale[3];
		float blo[3][BV_BINS][4], bhi[3][BV_BINS][4];
		unsigned int bcount[3][BV_BINS];
		bv_task kid[2];
		nd->first = t.begin;
		nd->cnt = n;
		if (t.depth > b->depth)
			b->depth = t.depth;
		for (axis = 0; axis < 3; axis++)
		{
			nd->lo[axis] = t.lo[axis];
			nd->hi[axis] = t.hi[axis];
		}
		if (n <= leaf_size)
		{
			nd->leaf = 1;
			nd->a = t.begin;
			nd->b = n;
			continue;
		}
		/* binned SAH, every axis the centroids spread along binned in one pass, each bin with its box and centroids' box */
		for (axis = 0; axis < 3; axis++)
		{
			scale[axis] = n > 4 && axis < b->dim && t.chi[axis] > t.clo[axis] ? BV_BINS / (t.chi[axis] - t.clo[axis]) : 0;
			for (j = 0; j < BV_BINS; j++)
			{
				bv_fempty4(blo[axis][j], bhi[axis][j]);
				bcount[axis][j] = 0;
			}
		}
		if (n > 4)
		{
#if BV_AVX2
			if (b->simd == 2)
				bv_bin_avx2(r, t.begin, t.end, b->dim, t.clo, scale, blo, bhi, bcount);
			else
#endif
				for (i = t.begin; i < t.end; i++)
					for (axis = 0; axis < b->dim; axis++)
						if (scale[axis] > 0)
						{
							int bin = bv_bin(r + i, axis, t.clo[axis], scale[axis]);
							bcount[axis][bin]++;
							bv_acc(blo[axis][bin], bhi[axis][bin], r[i].lo, r[i].hi);
						}
		}
		for (axis = 0; axis < b->dim; axis++)
		{
			float rlo[4], rhi[4], llo[4], lhi[4], rarea[BV_BINS];
			unsigned int rcount[BV_BINS], left = 0;
			if (!(scale[axis] > 0))
				continue;
			/* the right side of the plane after bin j - 1: bins j.. */
			bv_fempty4(rlo, rhi);
			for (j = BV_BINS - 1; j > 0; j--)
			{
				bv_acc(rlo, rhi, blo[axis][j], bhi[axis][j]);
				rcount[j - 1] = (j + 1 < BV_BINS ? rcount[j] : 0) + bcount[axis][j];
				rarea[j - 1] = bv_farea(rlo, rhi, b->dim);
			}
			bv_fempty4(llo, lhi);
			for (j = 0; j + 1 < BV_BINS; j++)
			{
				float cost;
				left += bcount[axis][j];
				bv_acc(llo, lhi, blo[axis][j], bhi[axis][j]);
				if (!left || !rcount[j])
					continue;
				cost = bv_farea(llo, lhi, b->dim) * left + rarea[j] * rcount[j];
				if (cost < best)
				{
					best = cost;
					best_axis = axis;
					best_bin = (int)j;
				}
			}
		}
		for (k = 0; k < 2; k++)
			kid[k].depth = t.depth + 1;
		mid = t.begin;
		if (best_axis < 3)
		{
			/* the children's bounds from the bins on each side */
			mid = bv_partition(r, best_axis, t.clo[best_axis], scale[best_axis], best_bin, t.begin, t.end);
			for (k = 0; k < 2; k++)
			{
				bv_fempty4(kid[k].lo, kid[k].hi);
				bv_fempty4(kid[k].clo, kid[k].chi);
			}
			for (j = 0; j < BV_BINS; j++)
			{
				k = (int)j <= best_bin ? 0 : 1;
				bv_acc(kid[k].lo, kid[k].hi, blo[best_axis][j], bhi[best_axis][j]);
			}
			/* the centroids (doubled) within the box, doubled */
			for (k = 0; k < 2; k++)
				for (j = 0; j < 4; j++)
				{
					kid[k].clo[j] = 2 * kid[k].lo[j];
					kid[k].chi[j] = 2 * kid[k].hi[j];
				}
		}
		kid[0].begin = t.begin;
		kid[1].end = t.end;
		if (mid == t.begin || mid == t.end)
		{
			/* no split by the bins (centroids together): the median along the widest spread */
			float w = -1;
			for (axis = 0, best_axis = 0; axis < b->dim; axis++)
				if (t.chi[axis] - t.clo[axis] > w)
				{
					w = t.chi[axis] - t.clo[axis];
					best_axis = axis;
				}
			bv_median(r, best_axis, t.begin, t.end);
			mid = t.begin + n / 2;
			kid[0].end = kid[1].begin = mid;
			bv_task_bounds(r, kid);
			bv_task_bounds(r, kid + 1);
		}
		kid[0].end = kid[1].begin = mid;
		nd->leaf = 0;
		nd->a = b->nn;
		nd->b = b->nn + 1;
		kid[0].node = nd->a;
		kid[1].node = nd->b;
		b->nn += 2;
		if (top + 2 > scap)
		{
			bv_task* grown = (bv_task*)qaws_internal_realloc(stack, 2 * (size_t)scap * sizeof(bv_task));
			if (!grown)
			{
				qaws_internal_free(r);
				qaws_internal_free(stack);
				return QAWS_STATUS_ALLOCATION_FAILURE;
			}
			stack = grown;
			scap *= 2;
		}
		stack[top++] = kid[1];
		stack[top++] = kid[0];
	}
	for (i = 0; i < b->ns; i++)
		b->ord[i] = r[i].id;
	qaws_internal_free(r);
	qaws_internal_free(stack);
	return QAWS_STATUS_OK;
}

/* entry e's (a position in ord) vertices, contiguous */
static double const* bv_pv(qaws_mesh_bvh const* b, unsigned int e, unsigned int j)
{
	return b->px + 3 * ((size_t)b->k * e + j);
}

/* entry e's oriented measure vector (area times unit normal; in 2D the length times the right-hand normal) and centroid */
static void bv_entry_moment(qaws_mesh_bvh const* b, unsigned int e, double* av, double* c)
{
	double const *A = bv_pv(b, e, 0), *B = bv_pv(b, e, 1);
	if (b->k == 3)
	{
		double const* C = bv_pv(b, e, 2);
		double e1[3], e2[3];
		bv_sub(B, A, e1);
		bv_sub(C, A, e2);
		bv_cross(e1, e2, av);
		av[0] *= 0.5;
		av[1] *= 0.5;
		av[2] *= 0.5;
		c[0] = (A[0] + B[0] + C[0]) / 3;
		c[1] = (A[1] + B[1] + C[1]) / 3;
		c[2] = (A[2] + B[2] + C[2]) / 3;
	}
	else
	{
		av[0] = B[1] - A[1];
		av[1] = -(B[0] - A[0]);
		av[2] = 0;
		c[0] = 0.5 * (A[0] + B[0]);
		c[1] = 0.5 * (A[1] + B[1]);
		c[2] = 0;
	}
}

/* the moments of entries [first, first + count) into m (BV_WM); T's lower half (k < j) left to bv_mirror */
static void bv_leaf_moment(qaws_mesh_bvh const* b, unsigned int first, unsigned int count, double* m)
{
	double area = 0, cs[3] = { 0, 0, 0 }, mean[3] = { 0, 0, 0 }, avs[8][3], cen[8][3], sw = b->k == 3 ? 1.0 / 12 : 1.0 / 6;
	unsigned int j, q, r, c2, c3;
	memset(m, 0, sizeof(double) * BV_WM);
	for (j = 0; j < count; j++)
	{
		double av[3], c[3], ar;
		bv_entry_moment(b, first + j, av, c);
		if (j < 8)
		{
			memcpy(avs[j], av, sizeof av);
			memcpy(cen[j], c, sizeof c);
		}
		ar = sqrt(bv_dot(av, av));
		area += ar;
		for (q = 0; q < 3; q++)
		{
			cs[q] += ar * c[q];
			mean[q] += c[q] / count;
			m[3 + q] += av[q];
		}
	}
	for (q = 0; q < 3; q++)
		m[q] = area > 0 ? cs[q] / area : mean[q];
	for (j = 0; j < count; j++)
	{
		double av[3], c[3], d[3], S[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		if (j < 8)
		{
			memcpy(av, avs[j], sizeof av);
			memcpy(c, cen[j], sizeof c);
		}
		else
			bv_entry_moment(b, first + j, av, c);
		bv_sub(c, m, d);
		for (r = 0; r < 3; r++)
			for (c2 = 0; c2 < 3; c2++)
				m[6 + 3 * r + c2] += av[r] * d[c2];
		for (r = 0; r < b->k; r++)
		{
			double e[3], v[3], l;
			bv_sub(bv_pv(b, first + j, r), m, e);
			l = bv_dot(e, e);
			if (l > m[15] * m[15])
				m[15] = sqrt(l);
			/* the simplex's own spread: its second moment about its centroid per unit measure */
			bv_sub(bv_pv(b, first + j, r), c, v);
			for (c2 = 0; c2 < 3; c2++)
				for (c3 = c2; c3 < 3; c3++)
					S[3 * c2 + c3] += sw * v[c2] * v[c3];
		}
		for (c2 = 0; c2 < 3; c2++)
			for (c3 = c2; c3 < 3; c3++)
			{
				double dd = d[c2] * d[c3] + S[3 * c2 + c3];
				for (r = 0; r < 3; r++)
					m[17 + 9 * r + 3 * c2 + c3] += av[r] * dd;
			}
	}
	m[16] = area;
}

/* T's lower half (k < j) from its upper one */
static void bv_mirror(double* m)
{
	unsigned int r, j, k;
	for (r = 0; r < 3; r++)
		for (j = 0; j < 3; j++)
			for (k = 0; k < j; k++)
				m[17 + 9 * r + 3 * j + k] = m[17 + 9 * r + 3 * k + j];
}

/* the moments of up to 4 parts merged into m: centre by area, C and T carried to it (T + C_ij s_k + C_ik s_j + N_i s_j s_k, s the shift);
   reads the parts' T upper halves, writes m's whole */
static void bv_merge_moments(double const* const* part, unsigned int n, double* m)
{
	double area = 0;
	unsigned int i, q, r, c2, c3;
	memset(m, 0, sizeof(double) * BV_WM);
	for (i = 0; i < n; i++)
		area += part[i][16];
	for (q = 0; q < 3; q++)
		for (i = 0; i < n; i++)
			m[q] += area > 0 ? part[i][16] * part[i][q] / area : part[i][q] / n;
	for (i = 0; i < n; i++)
	{
		double const* p = part[i];
		double s[3], rad;
		bv_sub(p, m, s);
		for (q = 0; q < 3; q++)
			m[3 + q] += p[3 + q];
		for (r = 0; r < 3; r++)
		{
			for (c2 = 0; c2 < 3; c2++)
				m[6 + 3 * r + c2] += p[6 + 3 * r + c2] + p[3 + r] * s[c2];
			for (c2 = 0; c2 < 3; c2++)
				for (c3 = c2; c3 < 3; c3++)
					m[17 + 9 * r + 3 * c2 + c3] += p[17 + 9 * r + 3 * c2 + c3] + p[6 + 3 * r + c2] * s[c3] + p[6 + 3 * r + c3] * s[c2]
						+ p[3 + r] * s[c2] * s[c3];
		}
		rad = p[15] + sqrt(bv_dot(s, s));
		if (rad > m[15])
			m[15] = rad;
	}
	m[16] = area;
	bv_mirror(m);
}

static void bv_slot_set(qaws_mesh_bvh const* b, bv_wide* t, unsigned int i, unsigned int s, double const* lo, double const* hi)
{
	unsigned int a;
	for (a = 0; a < 3; a++)
	{
		bv_nlo(t, i, a)[s] = bv_down(lo[a] - b->c[a] - b->pad);
		bv_nhi(t, i, a)[s] = bv_up(hi[a] - b->c[a] + b->pad);
	}
}

static void bv_slot_clear(bv_wide* t, unsigned int i, unsigned int s)
{
	unsigned int a;
	for (a = 0; a < 3; a++)
	{
		/* (lo +inf, hi -inf: every slab test, whatever the signs, and every box distance miss it) */
		bv_nlo(t, i, a)[s] = HUGE_VALF;
		bv_nhi(t, i, a)[s] = -HUGE_VALF;
	}
	bv_nch(t, i)[s] = QAWS_MESH_NONE;
	bv_ncnt(t, i)[s] = 0;
}

/* the binary tree (its leaves marked) collapsed into t->w-wide nodes, the largest inner child opened first, parents before children */
static qaws_status bv_emit(qaws_mesh_bvh const* b, bv_wide* t)
{
	unsigned int *stk, top;
	double root = bv_area(b->nodes[0].lo, b->nodes[0].hi, 3);
	t->nst = 8 * t->w;
	t->e = 0;
	if (!(root > 0))
		root = 1;
	t->n4 = 0;
	t->depth = 0;
	/* the nodes on cache lines; a wide node per interior binary node at most (each opens one), and the root */
	t->mem =qaws_internal_malloc(qaws_internal_add2(qaws_internal_mul3((size_t)b->nn / 2 + 2, t->nst, sizeof(float)), 64));
	t->nf = t->mem ? (float*)(((uintptr_t)t->mem + 63) & ~(uintptr_t)63) : NULL;
	stk = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3((size_t)b->nn + 1, 3, sizeof(unsigned int)));
	if (!t->nf || !stk)
	{
		qaws_internal_free(stk);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	t->n4 = 1;
	stk[0] = 0;
	stk[1] = 0;
	stk[2] = 1;
	top = 3;
	while (top)
	{
		unsigned int depth = stk[--top], n4 = stk[--top], bin = stk[--top], list[BV_WMAX], cnt = 0, s;
		t->e += bv_area(b->nodes[bin].lo, b->nodes[bin].hi, 3) / root;
		if (depth > t->depth)
			t->depth = depth;
		if (b->nodes[bin].lf)
			list[cnt++] = bin;
		else
		{
			list[cnt++] = b->nodes[bin].a;
			list[cnt++] = b->nodes[bin].b;
		}
		while (cnt < t->w)
		{
			unsigned int best = BV_WMAX;
			double ba = -1;
			for (s = 0; s < cnt; s++)
				if (!b->nodes[list[s]].lf && bv_area(b->nodes[list[s]].lo, b->nodes[list[s]].hi, 3) > ba)
				{
					ba = bv_area(b->nodes[list[s]].lo, b->nodes[list[s]].hi, 3);
					best = s;
				}
			if (best == BV_WMAX)
				break;
			s = list[best];
			list[best] = b->nodes[s].a;
			list[cnt++] = b->nodes[s].b;
		}
		for (s = 0; s < t->w; s++)
		{
			bv_node const* c;
			if (s >= cnt)
			{
				bv_slot_clear(t, n4, s);
				continue;
			}
			c = b->nodes + list[s];
			bv_slot_set(b, t, n4, s, c->lo, c->hi);
			if (c->lf)
			{
				bv_nch(t, n4)[s] = c->first;
				bv_ncnt(t, n4)[s] = c->cnt;
			}
			else
			{
				bv_nch(t, n4)[s] = t->n4++;
				bv_ncnt(t, n4)[s] = 0;
				stk[top++] = list[s];
				stk[top++] = bv_nch(t, n4)[s];
				stk[top++] = depth + 1;
			}
		}
	}
	qaws_internal_free(stk);
	return QAWS_STATUS_OK;
}

/* the leaves chosen on the binary tree, the entry arrays, then the wide trees: 4-wide for every query (8-wide too for rays on
   the AVX2 path); the binary nodes freed */
static qaws_status bv_collapse(qaws_mesh_bvh* b)
{
	unsigned int i;
	qaws_status st;
	double lo[3], hi[3];
	b->t.w = 4;
	if (!b->nn)
		return QAWS_STATUS_OK;
	/* the leaves: a subtree of at most leaf simplices becomes one where its surface area costs no more than its split
	   (a leaf screened 4 simplices at a time: its area per 4; a split: half its area for the node, plus its children) */
	{
		double* cost = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(b->nn, sizeof(double)));
		if (!cost)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		for (i = b->nn; i-- > 0;)
		{
			bv_node* nd = b->nodes + i;
			double A = bv_area(nd->lo, nd->hi, b->dim), as_leaf = A * ((nd->cnt + 3) / 4);
			if (nd->leaf)
			{
				nd->lf = 1;
				cost[i] = as_leaf;
			}
			else
			{
				double split = 0.5 * A + cost[nd->a] + cost[nd->b];
				nd->lf = nd->cnt <= b->leaf && as_leaf <= split;
				cost[i] = nd->lf ? as_leaf : split;
			}
		}
		qaws_internal_free(cost);
	}
	b->px = (double*)qaws_internal_malloc(qaws_internal_mul3((size_t)b->ns + 1, 3 * (size_t)b->k, sizeof(double)));
	if (b->k == 3 && b->dim == 3)
		b->scr = (float*)qaws_internal_zalloc(qaws_internal_mul2((size_t)b->ns + 8, 12), sizeof(float));
	b->ee = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2((size_t)b->ns + 1, sizeof(unsigned int)));
	b->ev = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3((size_t)b->ns + 1, b->k, sizeof(unsigned int)));
	if (!b->px || !b->ee || !b->ev || (b->k == 3 && b->dim == 3 && !b->scr))
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < b->ns; i++)
	{
		b->ee[i] = b->el[b->ord[i]];
		memcpy(b->ev + (size_t)b->k * i, b->s + (size_t)b->k * b->ord[i], b->k * sizeof(unsigned int));
	}
	/* the scene's centre and the boxes' widening (the floats' rounding of rays and points within it) */
	bv_empty(lo, hi);
	bv_grow(lo, hi, b->nodes[0].lo, b->nodes[0].hi);
	for (i = 0; i < 3; i++)
		b->c[i] = 0.5 * (lo[i] + hi[i]);
	b->pad = 1e-6 * (fabs(hi[0] - lo[0]) + fabs(hi[1] - lo[1]) + fabs(hi[2] - lo[2])) + 1e-30;
	st = bv_emit(b, &b->t);
	if (st == QAWS_STATUS_OK && b->simd == 2 && bv_winds(b) && b->k == 3)
	{
		b->r.w = 8;
		st = bv_emit(b, &b->r);
		/* kept where it saves rays enough node visits to pay for its wider tests (its expected visits under 0.69 of the
		   4-wide tree's: about 0.66 on smooth meshes, 0.72 on a soup where it measured slower); shortcut: a threshold from 3
		   scenes, upgrade to a cost model if a mesh family is found on the wrong side */
		if (st == QAWS_STATUS_OK && !(b->r.e < 0.69 * b->t.e))
		{
			qaws_internal_free(b->r.mem);
			memset(&b->r, 0, sizeof b->r);
		}
	}
	qaws_internal_free(b->nodes);
	b->nodes = NULL;
	return st;
}

/* a leaf slot's moments packed in floats (BV_LM): centre relative to the scene's, N, C, radius, T symmetric in its last two indices */
static void bv_pack_moment(qaws_mesh_bvh const* b, double const* m, float* p)
{
	unsigned int i, j, k, n = 16;
	for (i = 0; i < 3; i++)
		p[i] = (float)(m[i] - b->c[i]);
	for (i = 3; i < 16; i++)
		p[i] = (float)m[i];
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			for (k = j; k < 3; k++)
				p[n++] = (float)m[17 + 9 * i + 3 * j + k];
}

/* a wide tree's boxes (and, when w, its moments) from the leaves up */
static void bv_refit_tree(qaws_mesh_bvh* b, bv_wide* t, unsigned int w)
{
	unsigned int i, j, s, a;
	double leafm[BV_WMAX][BV_WM];
	for (i = t->n4; i-- > 0;)
	{
		unsigned int* ch = bv_nch(t, i), *cn = bv_ncnt(t, i);
		double const* part[BV_WMAX];
		unsigned int np = 0;
		for (s = 0; s < t->w; s++)
		{
			if (ch[s] == QAWS_MESH_NONE)
				continue;
			if (cn[s])
			{
				double lo[3], hi[3];
				bv_empty(lo, hi);
				for (j = 0; j < cn[s] * b->k; j++)
					bv_grow(lo, hi, b->px + 3 * ((size_t)b->k * ch[s] + j), b->px + 3 * ((size_t)b->k * ch[s] + j));
				bv_slot_set(b, t, i, s, lo, hi);
				if (w)
				{
					bv_leaf_moment(b, ch[s], cn[s], leafm[s]);
					bv_pack_moment(b, leafm[s], b->lm + BV_LM * ((size_t)t->w * i + s));
					part[np++] = leafm[s];
				}
			}
			else
			{
				/* the child's slots' union (its floats already rounded outward) */
				unsigned int u, c = ch[s];
				for (a = 0; a < 3; a++)
				{
					float l = HUGE_VALF, h = -HUGE_VALF;
					for (u = 0; u < t->w; u++)
						if (bv_nch(t, c)[u] != QAWS_MESH_NONE)
						{
							l = bv_nlo(t, c, a)[u] < l ? bv_nlo(t, c, a)[u] : l;
							h = bv_nhi(t, c, a)[u] > h ? bv_nhi(t, c, a)[u] : h;
						}
					bv_nlo(t, i, a)[s] = l;
					bv_nhi(t, i, a)[s] = h;
				}
				if (w)
					part[np++] = b->wm + BV_WM * (size_t)c;
			}
		}
		if (w)
		{
			/* the radius also bounded by the node's box: its farthest corner from the centre */
			double* m = b->wm + BV_WM * (size_t)i, r2 = 0;
			bv_merge_moments(part, np, m);
			for (a = 0; a < 3; a++)
			{
				double l = HUGE_VAL, h = -HUGE_VAL, e;
				for (s = 0; s < t->w; s++)
					if (ch[s] != QAWS_MESH_NONE)
					{
						l = fmin(l, (double)bv_nlo(t, i, a)[s] + b->c[a]);
						h = fmax(h, (double)bv_nhi(t, i, a)[s] + b->c[a]);
					}
				e = fmax(fabs(m[a] - l), fabs(h - m[a]));
				r2 += e * e;
			}
			if (sqrt(r2) < m[15])
				m[15] = sqrt(r2);
		}
	}
}

/* the leaf-ordered positions and the rays' screen, then the trees: at creation the 8-wide one keeps the binary tree's boxes
   it was emitted with (rounded outward and widened, so conservative), after moved vertices it is refit too */
static void bv_refit_nodes(qaws_mesh_bvh* b, int moved)
{
	unsigned int i, j, a;
	for (i = 0; i < b->ns; i++)
		for (j = 0; j < b->k; j++)
			memcpy(b->px + 3 * ((size_t)b->k * i + j), b->x + 3 * (size_t)b->s[(size_t)b->k * b->ord[i] + j], 3 * sizeof(double));
	if (b->scr)
		for (i = 0; i < b->ns; i++)
		{
			double const* v = b->px + 9 * (size_t)i;
			float e[12];
			for (a = 0; a < 3; a++)
			{
				e[a] = (float)(v[a] - b->c[a]);
				e[4 + a] = (float)(v[3 + a] - v[a]);
				e[8 + a] = (float)(v[6 + a] - v[a]);
			}
			/* the three vectors' L1 norms (the screen's bounds) */
			for (a = 0; a < 3; a++)
				e[4 * a + 3] = (fabsf(e[4 * a]) + fabsf(e[4 * a + 1])) + fabsf(e[4 * a + 2]);
			memcpy(b->scr + 12 * (size_t)i, e, sizeof e);
		}
	bv_refit_tree(b, &b->t, bv_winds(b) && b->wm);
	if (b->r.n4 && moved)
		bv_refit_tree(b, &b->r, 0);
}

/* ---- creation ---- */

/* takes x, s, el (freed on failure) */
static qaws_status bv_make(qaws_mesh_bvh** out, double* x, unsigned int nv, unsigned int dim, unsigned int* s, unsigned int k, unsigned int ns,
	unsigned int* el, qaws_mesh_bvh_desc const* desc)
{
	qaws_mesh_bvh* b = (qaws_mesh_bvh*)qaws_internal_zalloc(1, sizeof *b);
	qaws_status st;
	unsigned int leaf;
	double beta = desc && desc->winding_accuracy > 0 ? (double)desc->winding_accuracy : 2.25;
	if (!b)
	{
		qaws_internal_free(x);
		qaws_internal_free(s);
		qaws_internal_free(el);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	b->x = x;
	b->s = s;
	b->el = el;
	b->nv = nv;
	b->dim = dim;
	b->k = k;
	b->ns = ns;
	b->beta = beta;
	/* the path, once: AVX2 with 8-wide nodes, else SSE2 (or scalar) with 4 */
	b->simd = BV_SSE ? 1 : 0;
#if BV_AVX2
	if ((qaws_internal_cpu_features() & QAWS_CPU_AVX2) == QAWS_CPU_AVX2)
		b->simd = 2;
#endif
	leaf = desc && desc->leaf_size ? desc->leaf_size : 4;
	b->leaf = leaf;
	st = bv_build(b, 1);   /* (to single simplices: the collapse forms the leaves by cost) */
	if (st == QAWS_STATUS_OK)
		st = bv_collapse(b);
	if (st == QAWS_STATUS_OK && bv_winds(b) && b->t.n4)
	{
		b->wm = (double*)qaws_internal_malloc(qaws_internal_mul3((size_t)b->t.n4, BV_WM, sizeof(double)));
		b->lm = (float*)qaws_internal_malloc(qaws_internal_mul3((size_t)b->t.n4, 4 * BV_LM, sizeof(float)));
		if (!b->wm || !b->lm)
			st = QAWS_STATUS_ALLOCATION_FAILURE;
	}
	if (st != QAWS_STATUS_OK)
	{
		qaws_mesh_bvh_destroy(b);
		return st;
	}
	bv_refit_nodes(b, 0);
	*out = b;
	return QAWS_STATUS_OK;
}

/* positions (dim per vertex) into 3 doubles per vertex; INVALID_ARGUMENT when one is not finite */
static qaws_status bv_positions(double** out, qaws_scalar const* p, unsigned int nv, unsigned int dim)
{
	unsigned int i, j;
	double* x = (double*)qaws_internal_malloc(qaws_internal_mul3((size_t)nv + 1, 3, sizeof(double)));
	if (!x)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < nv; i++)
		for (j = 0; j < 3; j++)
		{
			x[3 * (size_t)i + j] = j < dim ? (double)p[(size_t)dim * i + j] : 0.0;
			if (!isfinite(x[3 * (size_t)i + j]))
			{
				qaws_internal_free(x);
				return QAWS_STATUS_INVALID_ARGUMENT;
			}
		}
	*out = x;
	return QAWS_STATUS_OK;
}

qaws_status qaws_mesh_bvh_create(qaws_mesh_bvh** out, qaws_mesh_bvh_desc const* d)
{
	double* x = NULL;
	unsigned int *s, *el, i;
	size_t n;
	qaws_status st;
	if (!out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	if (!d || (d->dimension != 2 && d->dimension != 3) || d->simplex_size < 2 || d->simplex_size > 4 || (d->simplex_size == 4 && d->dimension != 3)
		|| (d->vertex_count && !d->positions) || (d->simplex_count && !d->simplices))
		return QAWS_STATUS_INVALID_ARGUMENT;
	n = (size_t)d->simplex_size * d->simplex_count;
	for (i = 0; i < n; i++)
		if (d->simplices[i] >= d->vertex_count)
			return QAWS_STATUS_INVALID_ARGUMENT;
	st = bv_positions(&x, d->positions, d->vertex_count, d->dimension);
	if (st != QAWS_STATUS_OK)
		return st;
	s = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(n + 1, sizeof(unsigned int)));
	el = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2((size_t)d->simplex_count + 1, sizeof(unsigned int)));
	if (!s || !el)
	{
		qaws_internal_free(x);
		qaws_internal_free(s);
		qaws_internal_free(el);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	if (n)
		memcpy(s, d->simplices, n * sizeof(unsigned int));
	for (i = 0; i < d->simplex_count; i++)
		el[i] = d->elements ? d->elements[i] : i;
	return bv_make(out, x, d->vertex_count, d->dimension, s, d->simplex_size, d->simplex_count, el, d);
}

qaws_status qaws_mesh_bvh_create_mesh_3d(qaws_mesh_bvh** out, qaws_mesh_3d const* mesh, qaws_mesh_bvh_desc const* desc)
{
	unsigned int nv, nf, nc, ns = 0, f, c = 0, j, *fv, *fs, *s, *el;
	qaws_scalar* p;
	double* x = NULL;
	qaws_status st;
	if (!out || !mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	nv = qaws_mesh_3d_get_vertex_count(mesh);
	nf = qaws_mesh_3d_get_face_count(mesh);
	nc = qaws_mesh_3d_get_corner_count(mesh);
	p = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3((size_t)nv + 1, 3, sizeof(qaws_scalar)));
	fv = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2((size_t)nc + 1, sizeof(unsigned int)));
	fs = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2((size_t)nf + 1, sizeof(unsigned int)));
	s = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3((size_t)nc + 1, 3, sizeof(unsigned int)));
	el = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2((size_t)nc + 1, sizeof(unsigned int)));
	st = !p || !fv || !fs || !s || !el ? QAWS_STATUS_ALLOCATION_FAILURE : QAWS_STATUS_OK;
	if (st == QAWS_STATUS_OK)
		st = qaws_mesh_3d_get_positions(p, mesh);
	if (st == QAWS_STATUS_OK)
		st = qaws_mesh_3d_get_faces(fv, fs, mesh);
	if (st == QAWS_STATUS_OK)
		st = bv_positions(&x, p, nv, 3);
	if (st == QAWS_STATUS_OK)
		for (f = 0; f < nf; c += fs[f++])
			for (j = 1; j + 1 < fs[f]; j++, ns++)
			{
				s[3 * (size_t)ns] = fv[c];
				s[3 * (size_t)ns + 1] = fv[c + j];
				s[3 * (size_t)ns + 2] = fv[c + j + 1];
				el[ns] = f;
			}
	qaws_internal_free(p);
	qaws_internal_free(fv);
	qaws_internal_free(fs);
	if (st != QAWS_STATUS_OK)
	{
		qaws_internal_free(s);
		qaws_internal_free(el);
		return st;
	}
	return bv_make(out, x, nv, 3, s, 3, ns, el, desc);
}

qaws_status qaws_mesh_bvh_create_mesh_2d(qaws_mesh_bvh** out, qaws_mesh_2d const* mesh, int border, qaws_mesh_bvh_desc const* desc)
{
	unsigned int nv, nt, nq, ne, n, i, k, *s, *el;
	double* x;
	if (!out || !mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	nv = qaws_mesh_2d_get_vertex_count(mesh);
	nt = qaws_mesh_2d_get_triangle_count(mesh);
	nq = qaws_mesh_2d_get_quad_count(mesh);
	ne = qaws_mesh_2d_get_boundary_edge_count(mesh);
	k = border ? 2 : 3;
	n = border ? ne : nt + 2 * nq;
	x = (double*)qaws_internal_malloc(qaws_internal_mul3((size_t)nv + 1, 3, sizeof(double)));
	s = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3((size_t)n + 1, k, sizeof(unsigned int)));
	el = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2((size_t)n + 1, sizeof(unsigned int)));
	if (!x || !s || !el)
	{
		qaws_internal_free(x);
		qaws_internal_free(s);
		qaws_internal_free(el);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < nv; i++)
	{
		qaws_vec2 v = qaws_mesh_2d_get_vertex(mesh, i);
		x[3 * (size_t)i] = v.x;
		x[3 * (size_t)i + 1] = v.y;
		x[3 * (size_t)i + 2] = 0;
	}
	if (border)
		for (i = 0; i < ne; i++)
		{
			qaws_mesh_2d_get_boundary_edge(s + 2 * (size_t)i, mesh, i);
			el[i] = i;
		}
	else
	{
		for (i = 0; i < nt; i++)
		{
			qaws_mesh_2d_get_triangle(s + 3 * (size_t)i, mesh, i);
			el[i] = i;
		}
		for (i = 0; i < nq; i++)
		{
			unsigned int q[4], *t = s + 3 * ((size_t)nt + 2 * (size_t)i);
			qaws_mesh_2d_get_quad(q, mesh, i);
			t[0] = q[0];
			t[1] = q[1];
			t[2] = q[2];
			t[3] = q[0];
			t[4] = q[2];
			t[5] = q[3];
			el[nt + 2 * i] = el[nt + 2 * i + 1] = nt + i;
		}
	}
	return bv_make(out, x, nv, 2, s, k, n, el, desc);
}

/* each cell type's corners, then its facets not through corner 0 (as fans from corner 0 they tile the cell) */
static unsigned int const bv_cell_size[4] = { 4, 8, 6, 5 };
static unsigned int const bv_cell_tets[4][6][3] = {
	{ { 1, 2, 3 } },
	{ { 4, 5, 6 }, { 4, 6, 7 }, { 1, 2, 6 }, { 1, 6, 5 }, { 2, 3, 7 }, { 2, 7, 6 } },
	{ { 3, 4, 5 }, { 1, 2, 5 }, { 1, 5, 4 } },
	{ { 1, 2, 4 }, { 2, 3, 4 } }
};
static unsigned int const bv_cell_tet_count[4] = { 1, 6, 3, 2 };

qaws_status qaws_mesh_bvh_create_volume(qaws_mesh_bvh** out, qaws_volume_mesh const* mesh, qaws_mesh_bvh_desc const* desc)
{
	unsigned int nv, nc, ncorner, ns = 0, i, j, c = 0, *cv, *s, *el;
	qaws_cell_type* ct;
	qaws_scalar* p;
	double* x = NULL;
	qaws_status st;
	if (!out || !mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	nv = qaws_volume_mesh_get_vertex_count(mesh);
	nc = qaws_volume_mesh_get_cell_count(mesh);
	ncorner = qaws_volume_mesh_get_corner_count(mesh);
	p = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3((size_t)nv + 1, 3, sizeof(qaws_scalar)));
	cv = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2((size_t)ncorner + 1, sizeof(unsigned int)));
	ct = (qaws_cell_type*)qaws_internal_scratch_malloc(qaws_internal_mul2((size_t)nc + 1, sizeof(qaws_cell_type)));
	s = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3((size_t)nc + 1, 24, sizeof(unsigned int)));
	el = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3((size_t)nc + 1, 6, sizeof(unsigned int)));
	st = !p || !cv || !ct || !s || !el ? QAWS_STATUS_ALLOCATION_FAILURE : QAWS_STATUS_OK;
	if (st == QAWS_STATUS_OK)
		st = qaws_volume_mesh_get_positions(p, mesh);
	if (st == QAWS_STATUS_OK)
		st = qaws_volume_mesh_get_cells(cv, ct, mesh);
	if (st == QAWS_STATUS_OK)
		st = bv_positions(&x, p, nv, 3);
	if (st == QAWS_STATUS_OK)
		for (i = 0; i < nc; c += bv_cell_size[ct[i++]])
			for (j = 0; j < bv_cell_tet_count[ct[i]]; j++, ns++)
			{
				unsigned int const* f = bv_cell_tets[ct[i]][j];
				s[4 * (size_t)ns] = cv[c];
				s[4 * (size_t)ns + 1] = cv[c + f[0]];
				s[4 * (size_t)ns + 2] = cv[c + f[1]];
				s[4 * (size_t)ns + 3] = cv[c + f[2]];
				el[ns] = i;
			}
	qaws_internal_free(p);
	qaws_internal_free(cv);
	qaws_internal_free(ct);
	if (st != QAWS_STATUS_OK)
	{
		qaws_internal_free(s);
		qaws_internal_free(el);
		return st;
	}
	return bv_make(out, x, nv, 3, s, 4, ns, el, desc);
}

void qaws_mesh_bvh_destroy(qaws_mesh_bvh* b)
{
	if (!b)
		return;
	qaws_internal_free(b->x);
	qaws_internal_free(b->s);
	qaws_internal_free(b->el);
	qaws_internal_free(b->ord);
	qaws_internal_free(b->nodes);
	qaws_internal_free(b->t.mem);
	qaws_internal_free(b->r.mem);
	qaws_internal_free(b->px);
	qaws_internal_free(b->scr);
	qaws_internal_free(b->ee);
	qaws_internal_free(b->ev);
	qaws_internal_free(b->wm);
	qaws_internal_free(b->lm);
	qaws_internal_free(b);
}

qaws_status qaws_mesh_bvh_refit(qaws_mesh_bvh* b, qaws_scalar const* positions)
{
	double* x = NULL;
	qaws_status st;
	if (!b || (b->nv && !positions))
		return QAWS_STATUS_INVALID_ARGUMENT;
	st = bv_positions(&x, positions, b->nv, b->dim);
	if (st != QAWS_STATUS_OK)
		return st;
	qaws_internal_free(b->x);
	b->x = x;
	bv_refit_nodes(b, 1);
	return QAWS_STATUS_OK;
}

/* ---- point to simplex ---- */

static double bv_seg_closest(double const* p, double const* a, double const* b, double* t_out)
{
	double e[3], w[3], l2, t = 0, q[3];
	bv_sub(b, a, e);
	bv_sub(p, a, w);
	l2 = bv_dot(e, e);
	if (l2 > 0)
	{
		t = bv_dot(w, e) / l2;
		t = t < 0 ? 0 : t > 1 ? 1 : t;
	}
	q[0] = a[0] + t * e[0] - p[0];
	q[1] = a[1] + t * e[1] - p[1];
	q[2] = a[2] + t * e[2] - p[2];
	*t_out = t;
	return bv_dot(q, q);
}

/* squared distance from p to triangle abc, the closest point's weights on a, b, c (Ericson's regions; edges for a degenerate one) */
static double bv_tri_closest(double const* p, double const* a, double const* b, double const* c, double* w)
{
	double ab[3], ac[3], ap[3], bp[3], cp[3], d1, d2, d3, d4, d5, d6, va, vb, vc, den, v, u, q[3], n[3];
	bv_sub(b, a, ab);
	bv_sub(c, a, ac);
	bv_cross(ab, ac, n);
	if (!(bv_dot(n, n) > 1e-30 * bv_dot(ab, ab) * bv_dot(ac, ac)) || !(bv_dot(n, n) > 0))
	{
		/* degenerate: the nearest of its three edges */
		double t, best, dd;
		best = bv_seg_closest(p, a, b, &t);
		w[0] = 1 - t;
		w[1] = t;
		w[2] = 0;
		dd = bv_seg_closest(p, b, c, &t);
		if (dd < best)
		{
			best = dd;
			w[0] = 0;
			w[1] = 1 - t;
			w[2] = t;
		}
		dd = bv_seg_closest(p, c, a, &t);
		if (dd < best)
		{
			best = dd;
			w[0] = t;
			w[1] = 0;
			w[2] = 1 - t;
		}
		return best;
	}
	bv_sub(p, a, ap);
	d1 = bv_dot(ab, ap);
	d2 = bv_dot(ac, ap);
	if (d1 <= 0 && d2 <= 0)
	{
		w[0] = 1; w[1] = 0; w[2] = 0;
		return bv_dot(ap, ap);
	}
	bv_sub(p, b, bp);
	d3 = bv_dot(ab, bp);
	d4 = bv_dot(ac, bp);
	if (d3 >= 0 && d4 <= d3)
	{
		w[0] = 0; w[1] = 1; w[2] = 0;
		return bv_dot(bp, bp);
	}
	vc = d1 * d4 - d3 * d2;
	if (vc <= 0 && d1 >= 0 && d3 <= 0)
	{
		v = d1 / (d1 - d3);
		w[0] = 1 - v; w[1] = v; w[2] = 0;
		goto done;
	}
	bv_sub(p, c, cp);
	d5 = bv_dot(ab, cp);
	d6 = bv_dot(ac, cp);
	if (d6 >= 0 && d5 <= d6)
	{
		w[0] = 0; w[1] = 0; w[2] = 1;
		return bv_dot(cp, cp);
	}
	vb = d5 * d2 - d1 * d6;
	if (vb <= 0 && d2 >= 0 && d6 <= 0)
	{
		v = d2 / (d2 - d6);
		w[0] = 1 - v; w[1] = 0; w[2] = v;
		goto done;
	}
	va = d3 * d6 - d5 * d4;
	if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0)
	{
		v = (d4 - d3) / ((d4 - d3) + (d5 - d6));
		w[0] = 0; w[1] = 1 - v; w[2] = v;
		goto done;
	}
	den = 1 / (va + vb + vc);
	v = vb * den;
	u = vc * den;
	w[0] = 1 - v - u; w[1] = v; w[2] = u;
done:
	q[0] = w[0] * a[0] + w[1] * b[0] + w[2] * c[0] - p[0];
	q[1] = w[0] * a[1] + w[1] * b[1] + w[2] * c[1] - p[1];
	q[2] = w[0] * a[2] + w[1] * b[2] + w[2] * c[2] - p[2];
	return bv_dot(q, q);
}

static double bv_det3(double const* a, double const* b, double const* c)
{
	double n[3];
	bv_cross(b, c, n);
	return bv_dot(a, n);
}

/* squared distance from p to tetrahedron v[0..3] (0 inside), the weights on its corners */
static double bv_tet_closest(double const* p, double const* const* v, double* w)
{
	static unsigned int const face[4][3] = { { 1, 2, 3 }, { 0, 3, 2 }, { 0, 1, 3 }, { 0, 2, 1 } };
	double e1[3], e2[3], e3[3], q[3], vol, best = HUGE_VAL;
	unsigned int f, j;
	bv_sub(v[1], v[0], e1);
	bv_sub(v[2], v[0], e2);
	bv_sub(v[3], v[0], e3);
	bv_sub(p, v[0], q);
	vol = bv_det3(e1, e2, e3);
	if (vol != 0)
	{
		double l1 = bv_det3(q, e2, e3) / vol, l2 = bv_det3(e1, q, e3) / vol, l3 = bv_det3(e1, e2, q) / vol, l0 = 1 - l1 - l2 - l3;
		if (l0 >= 0 && l1 >= 0 && l2 >= 0 && l3 >= 0)
		{
			w[0] = l0; w[1] = l1; w[2] = l2; w[3] = l3;
			return 0;
		}
	}
	for (f = 0; f < 4; f++)
	{
		double fw[3], d = bv_tri_closest(p, v[face[f][0]], v[face[f][1]], v[face[f][2]], fw);
		if (d < best)
		{
			best = d;
			for (j = 0; j < 4; j++)
				w[j] = 0;
			for (j = 0; j < 3; j++)
				w[face[f][j]] = fw[j];
		}
	}
	return best;
}

static void bv_hit_none(qaws_mesh_hit* h)
{
	unsigned int j;
	h->element = h->simplex = QAWS_MESH_NONE;
	for (j = 0; j < 4; j++)
	{
		h->vertices[j] = QAWS_MESH_NONE;
		h->weights[j] = 0;
	}
	h->point[0] = h->point[1] = h->point[2] = 0;
	h->distance = (qaws_scalar)HUGE_VAL;
}

/* a hit on entry e: from the entry-ordered arrays (coherent with the traversal); the point from the weights unless given */
static void bv_hit_fill(qaws_mesh_bvh const* b, qaws_mesh_hit* h, unsigned int e, double const* w, double dist, double const* point)
{
	unsigned int j, q;
	double pt[3] = { 0, 0, 0 };
	h->element = b->ee[e];
	h->simplex = b->ord[e];
	for (j = 0; j < 4; j++)
	{
		h->vertices[j] = j < b->k ? b->ev[(size_t)b->k * e + j] : QAWS_MESH_NONE;
		h->weights[j] = j < b->k ? (qaws_scalar)w[j] : 0;
		if (j < b->k && !point)
			for (q = 0; q < 3; q++)
				pt[q] += w[j] * bv_pv(b, e, j)[q];
	}
	for (q = 0; q < 3; q++)
		h->point[q] = (qaws_scalar)(point ? point[q] : pt[q]);
	h->distance = (qaws_scalar)dist;
}

static double bv_entry_closest(qaws_mesh_bvh const* b, unsigned int e, double const* q, double* w)
{
	double const* v[4];
	unsigned int j;
	for (j = 0; j < b->k; j++)
		v[j] = bv_pv(b, e, j);
	if (b->k == 2)
	{
		double t, d = bv_seg_closest(q, v[0], v[1], &t);
		w[0] = 1 - t;
		w[1] = t;
		return d;
	}
	if (b->k == 3)
		return bv_tri_closest(q, v[0], v[1], v[2], w);
	return bv_tet_closest(q, v, w);
}

/* A node's slots against a point (squared distances) or a ray (entry distances and a mask of the hit slots), on the
   path's width: SSE2 for 4, AVX2 for 8, scalar for either. Every path does the same float operations in the same order,
   so they agree lane for lane. The AVX2 kernels clear the upper halves before returning (no AVX to SSE transitions). */

static void bv_slots_d2_any(float const* n, unsigned int w, float const* q, float* out)
{
	unsigned int a, s;
	for (s = 0; s < w; s++)
	{
		float d2 = 0;
		for (a = 0; a < 3; a++)
		{
			float dl = n[a * w + s] - q[a], dh = q[a] - n[(3 + a) * w + s], m = dl > dh ? dl : dh;
			m = m > 0 ? m : 0;
			d2 = d2 + m * m;
		}
		out[s] = d2;
	}
}

#if BV_SSE
static void bv_slots_d2_sse(float const* n, float const* q, float* out)
{
	__m128 d2 = _mm_setzero_ps(), z = _mm_setzero_ps();
	unsigned int a;
	for (a = 0; a < 3; a++)
	{
		__m128 qa = _mm_set1_ps(q[a]), dl = _mm_sub_ps(_mm_load_ps(n + 4 * a), qa), dh = _mm_sub_ps(qa, _mm_load_ps(n + 12 + 4 * a));
		__m128 m = _mm_max_ps(_mm_max_ps(dl, dh), z);
		d2 = _mm_add_ps(d2, _mm_mul_ps(m, m));
	}
	_mm_storeu_ps(out, d2);
}
#endif

static void bv_slots_d2(qaws_mesh_bvh const* b, float const* n, float const* q, float* out)
{
#if BV_SSE
	if (b->simd)
	{
		bv_slots_d2_sse(n, q, out);
		return;
	}
#endif
	bv_slots_d2_any(n, 4, q, out);
}

/* a ray as the node tests read it: origin (relative to the centre) and inverse direction, and the float offsets of each
   axis' near and far rows in a node (picked by the direction's signs, so no min or max per axis) */
typedef struct bv_rq
{
	float o[3], id[3];
	unsigned int nr[3], fr[3];
} bv_rq;

static unsigned int bv_slots_ray_any(float const* n, unsigned int w, bv_rq const* q, float tmax, float* tn, unsigned int* lf)
{
	unsigned int a, s, m = 0;
	unsigned int const* cn = (unsigned int const*)(n + 7 * w);
	*lf = 0;
	for (s = 0; s < w; s++)
	{
		*lf |= (unsigned int)(cn[s] != 0) << s;
		float t0 = 0, t1 = tmax;
		for (a = 0; a < 3; a++)
		{
			float tl = (n[q->nr[a] + s] - q->o[a]) * q->id[a], th = (n[q->fr[a] + s] - q->o[a]) * q->id[a];
			t0 = t0 > tl ? t0 : tl;
			t1 = t1 < th ? t1 : th;
		}
		tn[s] = t0;
		m |= (unsigned int)(t0 <= t1) << s;
	}
	return m;
}

#if BV_SSE
static unsigned int bv_slots_ray_sse(float const* n, bv_rq const* q, float tmax, float* tn, unsigned int* lf)
{
	__m128 t0 = _mm_setzero_ps(), t1 = _mm_set1_ps(tmax);
	unsigned int a;
	for (a = 0; a < 3; a++)
	{
		__m128 oa = _mm_set1_ps(q->o[a]), ia = _mm_set1_ps(q->id[a]);
		t0 = _mm_max_ps(t0, _mm_mul_ps(_mm_sub_ps(_mm_load_ps(n + q->nr[a]), oa), ia));
		t1 = _mm_min_ps(t1, _mm_mul_ps(_mm_sub_ps(_mm_load_ps(n + q->fr[a]), oa), ia));
	}
	_mm_storeu_ps(tn, t0);
	*lf = (unsigned int)_mm_movemask_ps(_mm_castsi128_ps(_mm_cmpgt_epi32(_mm_load_si128((__m128i const*)(n + 28)), _mm_setzero_si128())));
	return (unsigned int)_mm_movemask_ps(_mm_cmple_ps(t0, t1));
}
#endif

#if BV_AVX2
BV_AVX2_FN static unsigned int bv_slots_ray_avx2(float const* n, bv_rq const* q, float tmax, float* tn, unsigned int* lf)
{
	__m256 t0 = _mm256_setzero_ps(), t1 = _mm256_set1_ps(tmax);
	unsigned int a, m;
	for (a = 0; a < 3; a++)
	{
		__m256 oa = _mm256_set1_ps(q->o[a]), ia = _mm256_set1_ps(q->id[a]);
		t0 = _mm256_max_ps(t0, _mm256_mul_ps(_mm256_sub_ps(_mm256_load_ps(n + q->nr[a]), oa), ia));
		t1 = _mm256_min_ps(t1, _mm256_mul_ps(_mm256_sub_ps(_mm256_load_ps(n + q->fr[a]), oa), ia));
	}
	_mm256_storeu_ps(tn, t0);
	m = (unsigned int)_mm256_movemask_ps(_mm256_cmp_ps(t0, t1, _CMP_LE_OQ));
	*lf = (unsigned int)_mm256_movemask_ps(_mm256_castsi256_ps(_mm256_cmpgt_epi32(_mm256_load_si256((__m256i const*)(n + 56)), _mm256_setzero_si256())));
	_mm256_zeroupper();
	return m;
}
#endif

static unsigned int bv_slots_ray(qaws_mesh_bvh const* b, bv_wide const* t, float const* n, bv_rq const* q, float tmax, float* tn, unsigned int* lf)
{
#if BV_AVX2
	if (t->w == 8)
		return bv_slots_ray_avx2(n, q, tmax, tn, lf);
#endif
#if BV_SSE
	if (b->simd && t->w == 4)
		return bv_slots_ray_sse(n, q, tmax, tn, lf);
#endif
	return bv_slots_ray_any(n, t->w, q, tmax, tn, lf);
}

/* a node's leaf slots as a mask */
/* the w slots with a finite key, nearest first */
static unsigned int bv_order(float const* key, unsigned int w, unsigned int* slot)
{
	unsigned int n = 0, s, i;
	for (s = 0; s < w; s++)
		if (key[s] < HUGE_VALF)
		{
			for (i = n; i > 0 && key[slot[i - 1]] > key[s]; i--)
				slot[i] = slot[i - 1];
			slot[i] = s;
			n++;
		}
	return n;
}

/* a float bound just above a double (for pruning against floats) */
static float bv_fbound(double v)
{
	float f;
	if (!(v < 3e38))
		return HUGE_VALF;
	f = (float)(v * (1 + 1e-6));
	return (double)f >= v ? f : bv_up(v);
}

/* the closest entry to q with squared distance within *best2 (updated); its weights into bw. Returns the entry. */
static unsigned int bv_closest_one(qaws_mesh_bvh const* b, double const* q, double* best2, double* bw, bv_stk* stk)
{
	unsigned int top = 0, found = QAWS_MESH_NONE, i, s;
	float qf[3] = { (float)(q[0] - b->c[0]), (float)(q[1] - b->c[1]), (float)(q[2] - b->c[2]) }, bf = bv_fbound(*best2);
	if (!b->t.n4)
		return found;
	stk[top].node = 0;
	stk[top++].t = 0;
	while (top)
	{
		bv_stk e = stk[--top];
		float const* n;
		unsigned int const *ch, *cn;
		float key[BV_WMAX];
		unsigned int slot[BV_WMAX], ns;
		if (e.t > bf)
			continue;
		n = b->t.nf + (size_t)e.node * b->t.nst;
		ch = (unsigned int const*)(n + 24);
		cn = ch + 4;
		bv_slots_d2(b, n, qf, key);
		for (s = 0; s < 4; s++)
			if (key[s] > bf)
				key[s] = HUGE_VALF;
		ns = bv_order(key, 4, slot);
		for (i = 0; i < ns; i++)
		{
			unsigned int sl = slot[i], c;
			if (!cn[sl] || key[sl] > bf)
				continue;
			for (c = 0; c < cn[sl]; c++)
			{
				double w[4], d = bv_entry_closest(b, ch[sl] + c, q, w);
				if (d <= *best2 && (found == QAWS_MESH_NONE || d < *best2))
				{
					*best2 = d;
					found = ch[sl] + c;
					memcpy(bw, w, sizeof w);
					bf = bv_fbound(d);
				}
			}
			if (*best2 == 0)
				return found;
		}
		for (i = ns; i-- > 0;)
			if (!cn[slot[i]] && key[slot[i]] <= bf)
			{
				stk[top].node = ch[slot[i]];
				stk[top++].t = key[slot[i]];
			}
	}
	return found;
}

/* The leaf triangles' float screen: Moller-Trumbore on 4 entries at a time (v0 relative to the centre, e1, e2), each
   test widened by 1e-5 of a bound on its terms (30 times the floats' rounding, inputs included), so it can only let
   through more than the exact test; the exact one decides. Bit l of the result: entry first + l is a candidate. */
static unsigned int bv_screen(qaws_mesh_bvh const* b, unsigned int first, unsigned int n, float const* o, float const* d, float tmax)
{
	float const* tf = b->scr;
#if BV_SSE
	__m128 sign = _mm_set1_ps(-0.0f), tol = _mm_set1_ps(1e-5f), z = _mm_setzero_ps();
	__m128 dx = _mm_set1_ps(d[0]), dy = _mm_set1_ps(d[1]), dz = _mm_set1_ps(d[2]);
	__m128 ox = _mm_set1_ps(o[0]), oy = _mm_set1_ps(o[1]), oz = _mm_set1_ps(o[2]);
	__m128 da = _mm_set1_ps(fabsf(d[0]) + fabsf(d[1]) + fabsf(d[2])), oa = _mm_set1_ps(fabsf(o[0]) + fabsf(o[1]) + fabsf(o[2]));
	/* 4 entries of 12 floats (v0, e1, e2, each padded to 4), transposed to lanes */
	float const* p0 = tf + 12 * (size_t)first;
	__m128 vx = _mm_loadu_ps(p0), vy = _mm_loadu_ps(p0 + 12), vz = _mm_loadu_ps(p0 + 24), vw = _mm_loadu_ps(p0 + 36);
	__m128 ax = _mm_loadu_ps(p0 + 4), ay = _mm_loadu_ps(p0 + 16), az = _mm_loadu_ps(p0 + 28), aw = _mm_loadu_ps(p0 + 40);
	__m128 bx = _mm_loadu_ps(p0 + 8), by = _mm_loadu_ps(p0 + 20), bz = _mm_loadu_ps(p0 + 32), bw = _mm_loadu_ps(p0 + 44);
	_MM_TRANSPOSE4_PS(vx, vy, vz, vw);
	_MM_TRANSPOSE4_PS(ax, ay, az, aw);
	_MM_TRANSPOSE4_PS(bx, by, bz, bw);
	/* pvec = d x e2, det = e1 . pvec, tvec = o - v0, u = tvec . pvec, qvec = tvec x e1, v = d . qvec, t = e2 . qvec */
	__m128 px = _mm_sub_ps(_mm_mul_ps(dy, bz), _mm_mul_ps(dz, by)), py = _mm_sub_ps(_mm_mul_ps(dz, bx), _mm_mul_ps(dx, bz)),
		pz = _mm_sub_ps(_mm_mul_ps(dx, by), _mm_mul_ps(dy, bx));
	__m128 det = _mm_add_ps(_mm_add_ps(_mm_mul_ps(ax, px), _mm_mul_ps(ay, py)), _mm_mul_ps(az, pz));
	__m128 tx = _mm_sub_ps(ox, vx), ty = _mm_sub_ps(oy, vy), tz = _mm_sub_ps(oz, vz);
	__m128 u = _mm_add_ps(_mm_add_ps(_mm_mul_ps(tx, px), _mm_mul_ps(ty, py)), _mm_mul_ps(tz, pz));
	__m128 qx = _mm_sub_ps(_mm_mul_ps(ty, az), _mm_mul_ps(tz, ay)), qy = _mm_sub_ps(_mm_mul_ps(tz, ax), _mm_mul_ps(tx, az)),
		qz = _mm_sub_ps(_mm_mul_ps(tx, ay), _mm_mul_ps(ty, ax));
	__m128 v = _mm_add_ps(_mm_add_ps(_mm_mul_ps(dx, qx), _mm_mul_ps(dy, qy)), _mm_mul_ps(dz, qz));
	__m128 t = _mm_add_ps(_mm_add_ps(_mm_mul_ps(bx, qx), _mm_mul_ps(by, qy)), _mm_mul_ps(bz, qz));
	/* the bounds: |tvec| with the rounding of o and v0, |e1|, |e2| (L1) */
	__m128 tn = _mm_add_ps(_mm_add_ps(_mm_andnot_ps(sign, tx), _mm_andnot_ps(sign, ty)), _mm_add_ps(_mm_andnot_ps(sign, tz), oa));
	__m128 e1n = aw, e2n = bw, vn = vw;   /* (the norms, stored in the pads) */
	__m128 tu, tt, s, m;
	tn = _mm_add_ps(tn, vn);
	tu = _mm_mul_ps(tol, _mm_mul_ps(_mm_add_ps(tn, e1n), _mm_mul_ps(da, _mm_add_ps(e1n, e2n))));
	tt = _mm_mul_ps(tol, _mm_mul_ps(tn, _mm_mul_ps(e1n, e2n)));
	/* to det's sign */
	s = _mm_and_ps(det, sign);
	det = _mm_xor_ps(det, s);
	u = _mm_xor_ps(u, s);
	v = _mm_xor_ps(v, s);
	t = _mm_xor_ps(t, s);
	m = _mm_and_ps(_mm_cmpge_ps(_mm_add_ps(u, tu), z), _mm_cmpge_ps(_mm_add_ps(v, tu), z));
	m = _mm_and_ps(m, _mm_cmple_ps(_mm_add_ps(u, v), _mm_add_ps(det, _mm_add_ps(tu, tu))));
	m = _mm_and_ps(m, _mm_cmpge_ps(_mm_add_ps(t, tt), z));
	m = _mm_and_ps(m, _mm_cmple_ps(t, _mm_add_ps(_mm_mul_ps(_mm_set1_ps(tmax), _mm_add_ps(det, tu)), tt)));
	return (unsigned int)_mm_movemask_ps(m) & ((1u << n) - 1);
#else
	unsigned int l, mask = 0;
	float da = fabsf(d[0]) + fabsf(d[1]) + fabsf(d[2]), oa = fabsf(o[0]) + fabsf(o[1]) + fabsf(o[2]);
	for (l = 0; l < n; l++)
	{
		float const* p0 = tf + 12 * ((size_t)first + l);
		float vx = p0[0], vy = p0[1], vz = p0[2], ax = p0[4], ay = p0[5], az = p0[6], bx = p0[8], by = p0[9], bz = p0[10];
		float px = d[1] * bz - d[2] * by, py = d[2] * bx - d[0] * bz, pz = d[0] * by - d[1] * bx;
		float det = ax * px + ay * py + az * pz, tx = o[0] - vx, ty = o[1] - vy, tz = o[2] - vz;
		float u = tx * px + ty * py + tz * pz;
		float qx = ty * az - tz * ay, qy = tz * ax - tx * az, qz = tx * ay - ty * ax;
		float v = d[0] * qx + d[1] * qy + d[2] * qz, t = bx * qx + by * qy + bz * qz;
		float tn = (fabsf(tx) + fabsf(ty)) + (fabsf(tz) + oa) + p0[3], e1n = p0[7], e2n = p0[11];
		float tu = 1e-5f * ((tn + e1n) * (da * (e1n + e2n))), tt = 1e-5f * (tn * (e1n * e2n));
		if (det < 0)
		{
			det = -det;
			u = -u;
			v = -v;
			t = -t;
		}
		if (u + tu >= 0 && v + tu >= 0 && u + v <= det + (tu + tu) && t + tt >= 0 && t <= tmax * (det + tu) + tt)
			mask |= 1u << l;
	}
	return mask;
#endif
}

/* Woop, Benthin and Wald's watertight test; t in [0, tmax] and the weights on a, b, c */
static int bv_ray_tri(double const* o, unsigned int const* ax, double const* sh, double const* a, double const* b, double const* c,
	double tmax, double* t_out, double* w)
{
	unsigned int kx = ax[0], ky = ax[1], kz = ax[2];
	double A[3], B[3], C[3], Ax, Ay, Bx, By, Cx, Cy, U, V, W, det, T;
	bv_sub(a, o, A);
	bv_sub(b, o, B);
	bv_sub(c, o, C);
	Ax = A[kx] - sh[0] * A[kz];
	Ay = A[ky] - sh[1] * A[kz];
	Bx = B[kx] - sh[0] * B[kz];
	By = B[ky] - sh[1] * B[kz];
	Cx = C[kx] - sh[0] * C[kz];
	Cy = C[ky] - sh[1] * C[kz];
	U = Cx * By - Cy * Bx;
	V = Ax * Cy - Ay * Cx;
	W = Bx * Ay - By * Ax;
	if ((U < 0 || V < 0 || W < 0) && (U > 0 || V > 0 || W > 0))
		return 0;
	det = U + V + W;
	if (det == 0)
		return 0;
	T = U * sh[2] * A[kz] + V * sh[2] * B[kz] + W * sh[2] * C[kz];
	*t_out = T / det;
	if (!(*t_out >= 0 && *t_out <= tmax))
		return 0;
	w[0] = U / det;
	w[1] = V / det;
	w[2] = W / det;
	return 1;
}

/* a 2D ray against segment ab: sides by orientation with the same rule at a shared vertex, so a ray through it hits one side */
static int bv_ray_seg(double const* o, double const* d, double const* a, double const* b, double tmax, double* t_out, double* w)
{
	double sa = d[0] * (a[1] - o[1]) - d[1] * (a[0] - o[0]), sb = d[0] * (b[1] - o[1]) - d[1] * (b[0] - o[0]), e[2], den, t, s;
	if ((sa >= 0) == (sb >= 0))
		return 0;
	e[0] = b[0] - a[0];
	e[1] = b[1] - a[1];
	den = d[0] * e[1] - d[1] * e[0];
	if (den == 0)
		return 0;
	t = ((a[0] - o[0]) * e[1] - (a[1] - o[1]) * e[0]) / den;
	if (!(t >= 0 && t <= tmax))
		return 0;
	s = sa / (sa - sb);
	*t_out = t;
	w[0] = 1 - s;
	w[1] = s;
	return 1;
}

static unsigned int bv_ray_one(qaws_mesh_bvh const* b, double const* o, double const* d, double tmax, int any, double* t_best, double* bw,
	bv_stk* stk)
{
	double sh[3] = { 0, 0, 0 };
	bv_wide const* tr = b->r.n4 ? &b->r : &b->t;
	unsigned int ax[3] = { 0, 1, 2 }, top = 0, found = QAWS_MESH_NONE, i, cand = 0, node = 0, w = tr->w;
	float of[3], df[3], tf;
	bv_rq q;
	for (i = 0; i < 3; i++)
	{
		of[i] = q.o[i] = (float)(o[i] - b->c[i]);
		df[i] = (float)d[i];
		q.id[i] = d[i] != 0 ? (float)(1 / d[i]) : 1e38f;
		q.nr[i] = (q.id[i] < 0 ? 3 + i : i) * w;
		q.fr[i] = (q.id[i] < 0 ? i : 3 + i) * w;
	}
	if (b->k == 3)
	{
		unsigned int kz = fabs(d[0]) > fabs(d[1]) ? (fabs(d[0]) > fabs(d[2]) ? 0 : 2) : (fabs(d[1]) > fabs(d[2]) ? 1 : 2);
		unsigned int kx = kz == 2 ? 0 : kz + 1, ky = kx == 2 ? 0 : kx + 1;
		if (d[kz] < 0)
		{
			unsigned int t = kx;
			kx = ky;
			ky = t;
		}
		ax[0] = kx;
		ax[1] = ky;
		ax[2] = kz;
		sh[0] = d[kx] / d[kz];
		sh[1] = d[ky] / d[kz];
		sh[2] = 1 / d[kz];
	}
	*t_best = tmax;
	tf = bv_fbound(tmax);
	if (!tr->n4)
		return found;
	for (;;)
	{
		float const* n = tr->nf + (size_t)node * tr->nst;
		unsigned int const* ch = (unsigned int const*)(n + 6 * w), *cn = ch + w;
		float tn[BV_WMAX], near = HUGE_VALF;
		unsigned int lf, m = bv_slots_ray(b, tr, n, &q, tf, tn, &lf), s, next = QAWS_MESH_NONE;
		/* the leaves now, nearest first (any hit: in slot order) */
		unsigned int lm = m & lf, im = m & ~lm;
		while (lm)
		{
			unsigned int sl = bv_ctz(lm), c;
			if (!any && (lm & (lm - 1)))
				for (s = lm & (lm - 1); s; s &= s - 1)
					if (tn[bv_ctz(s)] < tn[sl])
						sl = bv_ctz(s);
			lm &= ~(1u << sl);
			if (tn[sl] > tf)
				continue;
			for (c = 0; c < cn[sl]; c++)
			{
				unsigned int en = ch[sl] + c;
				double t, wt[4] = { 0, 0, 0, 0 };
				int hit;
				if (b->scr && !(c & 3))
					cand = bv_screen(b, en, cn[sl] - c < 4 ? cn[sl] - c : 4, of, df, tf);
				if (b->scr && !(cand >> (c & 3) & 1))
					continue;
				hit = b->k == 3 ? bv_ray_tri(o, ax, sh, bv_pv(b, en, 0), bv_pv(b, en, 1), bv_pv(b, en, 2), *t_best, &t, wt)
				                : bv_ray_seg(o, d, bv_pv(b, en, 0), bv_pv(b, en, 1), *t_best, &t, wt);
				if (hit && (found == QAWS_MESH_NONE || t < *t_best))
				{
					*t_best = t;
					found = en;
					memcpy(bw, wt, sizeof wt);
					tf = bv_fbound(t);
					if (any)
						return found;
				}
			}
		}
		/* the inner slots: on into the nearest, the others pushed */
		for (; im; im &= im - 1)
		{
			s = bv_ctz(im);
			if (tn[s] > tf)
				continue;
#if BV_SSE
			/* (the child on its way while this node's others are sorted out) */
			{
				char const* cp = (char const*)(tr->nf + (size_t)ch[s] * tr->nst);
				_mm_prefetch(cp, _MM_HINT_T0);
				_mm_prefetch(cp + 64, _MM_HINT_T0);
				if (w == 8)
				{
					_mm_prefetch(cp + 128, _MM_HINT_T0);
					_mm_prefetch(cp + 192, _MM_HINT_T0);
				}
			}
#endif
			if (tn[s] < near)
			{
				if (next != QAWS_MESH_NONE)
				{
					stk[top].node = next;
					stk[top++].t = near;
				}
				next = ch[s];
				near = tn[s];
			}
			else
			{
				stk[top].node = ch[s];
				stk[top++].t = tn[s];
			}
		}
		if (next == QAWS_MESH_NONE)
		{
			do
			{
				if (!top)
					return found;
				--top;
			} while (stk[top].t > tf);
			next = stk[top].node;
		}
		node = next;
	}
}

/* ---- winding number ---- */

static double bv_entry_winding(qaws_mesh_bvh const* b, unsigned int e, double const* q)
{
	double A[3], B[3];
	bv_sub(bv_pv(b, e, 0), q, A);
	bv_sub(bv_pv(b, e, 1), q, B);
	if (b->k == 3)
	{
		/* Van Oosterom and Strackee's solid angle */
		double C[3], la, lb, lc, num, den;
		bv_sub(bv_pv(b, e, 2), q, C);
		la = sqrt(bv_dot(A, A));
		lb = sqrt(bv_dot(B, B));
		lc = sqrt(bv_dot(C, C));
		num = bv_det3(A, B, C);
		den = la * lb * lc + bv_dot(A, B) * lc + bv_dot(B, C) * la + bv_dot(C, A) * lb;
		return 2 * atan2(num, den) / (4 * BV_PI);
	}
	return atan2(A[0] * B[1] - A[1] * B[0], A[0] * B[0] + A[1] * B[1]) / (2 * BV_PI);
}

/* the expansion to second order at r = centre - q (R2 = |r|^2): g . N + J : C + H : T / 2, g(r) = r / (4 pi R^3) (2D: r / (2 pi R^2)),
   J_ij = a1 d_ij + a2 r_i r_j and H_ijk = a2 (d_ij r_k + d_ik r_j + d_jk r_i) + a3 r_i r_j r_k contracted in closed form;
   Ts: T's 6 entries (j <= k) per i */
static double bv_expand(unsigned int dim, double const* rv, double R2, double const* N, double const* C, double const* Ts)
{
	double R = sqrt(R2), g, a1, a2, a3, rCr = 0, trC = C[0] + C[4] + C[8], diag = 0, rtr = 0, rTr = 0;
	unsigned int i, j;
	if (dim == 3)
	{
		double R3 = R2 * R, R5 = R3 * R2, R7 = R5 * R2;
		g = 1 / (4 * BV_PI * R3);
		a1 = 1 / R3 / (4 * BV_PI);
		a2 = -3 / R5 / (4 * BV_PI);
		a3 = 15 / R7 / (4 * BV_PI);
	}
	else
	{
		double R4 = R2 * R2, R6 = R4 * R2;
		g = 1 / (2 * BV_PI * R2);
		a1 = 1 / R2 / (2 * BV_PI);
		a2 = -2 / R4 / (2 * BV_PI);
		a3 = 8 / R6 / (2 * BV_PI);
	}
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			rCr += rv[i] * C[3 * i + j] * rv[j];
	for (i = 0; i < 3; i++)
	{
		double const* t = Ts + 6 * i;   /* (00, 01, 02, 11, 12, 22) */
		double u0 = t[0] * rv[0] + t[1] * rv[1] + t[2] * rv[2], u1 = t[1] * rv[0] + t[3] * rv[1] + t[4] * rv[2],
			u2 = t[2] * rv[0] + t[4] * rv[1] + t[5] * rv[2], u[3] = { u0, u1, u2 };
		diag += u[i];
		rtr += rv[i] * (t[0] + t[3] + t[5]);
		rTr += rv[i] * (rv[0] * u0 + rv[1] * u1 + rv[2] * u2);
	}
	return g * bv_dot(rv, N) + (a1 * trC + a2 * rCr) + 0.5 * (a2 * (2 * diag + rtr) + a3 * rTr);
}

/* node moments m (BV_WM, T whole) at q: when q is far (beyond beta times the radius), adds the expansion to *w and returns 1 */
static int bv_far(qaws_mesh_bvh const* b, double const* m, double const* q, double* w)
{
	double rv[3], R2, Ts[18];
	unsigned int i, j, k, n = 0;
	bv_sub(m, q, rv);
	R2 = bv_dot(rv, rv);
	if (!(R2 > b->beta * b->beta * m[15] * m[15]))
		return 0;
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			for (k = j; k < 3; k++)
				Ts[n++] = m[17 + 9 * i + 3 * j + k];
	*w += bv_expand(b->dim, rv, R2, m + 3, m + 6, Ts);
	return 1;
}

/* the same on a leaf slot's packed moments (BV_LM floats) */
static int bv_far_leaf(qaws_mesh_bvh const* b, float const* p, double const* q, double* w)
{
	double rv[3], R2, N[3], C[9], Ts[18];
	unsigned int i;
	for (i = 0; i < 3; i++)
		rv[i] = ((double)p[i] + b->c[i]) - q[i];
	R2 = bv_dot(rv, rv);
	if (!(R2 > b->beta * b->beta * (double)p[15] * (double)p[15]))
		return 0;
	for (i = 0; i < 3; i++)
		N[i] = p[3 + i];
	for (i = 0; i < 9; i++)
		C[i] = p[6 + i];
	for (i = 0; i < 18; i++)
		Ts[i] = p[16 + i];
	*w += bv_expand(b->dim, rv, R2, N, C, Ts);
	return 1;
}

static double bv_winding_one(qaws_mesh_bvh const* b, double const* q, bv_stk* stk)
{
	unsigned int top = 0, i, s;
	double w = 0;
	if (!b->t.n4)
		return 0;
	stk[top++].node = 0;
	while (top)
	{
		unsigned int id = stk[--top].node, *ch = bv_nch(&b->t, id), *cn = bv_ncnt(&b->t, id);
		if (bv_far(b, b->wm + BV_WM * (size_t)id, q, &w))
			continue;
		for (s = 0; s < 4; s++)
		{
			if (ch[s] == QAWS_MESH_NONE)
				continue;
			if (cn[s])
			{
				if (!bv_far_leaf(b, b->lm + BV_LM * (4 * (size_t)id + s), q, &w))
					for (i = 0; i < cn[s]; i++)
						w += bv_entry_winding(b, ch[s] + i, q);
			}
			else
				stk[top++].node = ch[s];
		}
	}
	return w;
}

/* ---- batches ---- */

typedef struct bv_job
{
	qaws_mesh_bvh const* b;
	qaws_scalar const* p;     /* points or origins */
	qaws_scalar const* d;     /* directions */
	unsigned int count;
	double tmax2;             /* closest: squared max distance; ray: tmax */
	int kind, any;            /* 0 closest, 1 ray, 2 winding, 3 signed */
	qaws_mesh_hit* hits;
	qaws_scalar* values;
	bv_stk* stacks;           /* bv_stack_len per chunk */
} bv_job;

static void bv_read(qaws_mesh_bvh const* b, qaws_scalar const* p, unsigned int i, double* q)
{
	q[0] = (double)p[(size_t)b->dim * i];
	q[1] = (double)p[(size_t)b->dim * i + 1];
	q[2] = b->dim == 3 ? (double)p[(size_t)b->dim * i + 2] : 0.0;
}

/* a query's stack: (w - 1) depth + 2 entries for the larger of the trees */
static size_t bv_stack_len(qaws_mesh_bvh const* b)
{
	size_t a = 3 * (size_t)b->t.depth + 2, r = 7 * (size_t)b->r.depth + 2;
	return a > r ? a : r;
}

static void bv_task_run(void* ctx, unsigned int begin, unsigned int end)
{
	bv_job const* j = (bv_job const*)ctx;
	qaws_mesh_bvh const* b = j->b;
	unsigned int ch, i;
	for (ch = begin; ch < end; ch++)
	{
		bv_stk* stack = j->stacks + (size_t)ch * bv_stack_len(b);
		unsigned int last = (ch + 1) * BV_CHUNK < j->count ? (ch + 1) * BV_CHUNK : j->count;
		for (i = ch * BV_CHUNK; i < last; i++)
		{
			double q[3], w[4] = { 0, 0, 0, 0 };
			if (j->kind == 0 || j->kind == 3)
			{
				double best2 = j->kind == 0 ? j->tmax2 : HUGE_VAL;
				unsigned int p;
				qaws_mesh_hit h;
				bv_read(b, j->p, i, q);
				p = bv_closest_one(b, q, &best2, w, stack);
				bv_hit_none(&h);
				if (p != QAWS_MESH_NONE)
					bv_hit_fill(b, &h, p, w, sqrt(best2), NULL);
				if (j->kind == 3)
				{
					double wn = bv_winding_one(b, q, stack);
					j->values[i] = (qaws_scalar)(wn >= 0.5 ? -sqrt(best2) : sqrt(best2));
					if (p == QAWS_MESH_NONE)
						j->values[i] = (qaws_scalar)HUGE_VAL;
				}
				if (j->hits)
					j->hits[i] = h;
			}
			else if (j->kind == 1)
			{
				double d[3], t = 0;
				unsigned int p;
				bv_read(b, j->p, i, q);
				bv_read(b, j->d, i, d);
				p = (d[0] != 0 || d[1] != 0 || d[2] != 0) ? bv_ray_one(b, q, d, j->tmax2, j->any, &t, w, stack) : QAWS_MESH_NONE;
				bv_hit_none(&j->hits[i]);
				if (p != QAWS_MESH_NONE)
				{
					double pt[3] = { q[0] + t * d[0], q[1] + t * d[1], q[2] + t * d[2] };
					bv_hit_fill(b, &j->hits[i], p, w, t, pt);
				}
			}
			else
			{
				bv_read(b, j->p, i, q);
				j->values[i] = (qaws_scalar)bv_winding_one(b, q, stack);
			}
		}
	}
}

static qaws_status bv_run(bv_job* j, qaws_batch_executor const* executor)
{
	unsigned int chunks = (j->count + BV_CHUNK - 1) / BV_CHUNK;
	if (!j->count)
		return QAWS_STATUS_OK;
	j->stacks = (bv_stk*)qaws_internal_scratch_malloc(qaws_internal_mul3(chunks, bv_stack_len(j->b), sizeof(bv_stk)));
	if (!j->stacks)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	if (executor && executor->parallel_for)
		executor->parallel_for(executor->user, chunks, bv_task_run, j);
	else
		bv_task_run(j, 0, chunks);
	qaws_internal_free(j->stacks);
	return QAWS_STATUS_OK;
}

qaws_status qaws_mesh_bvh_closest(qaws_mesh_hit* out, qaws_mesh_bvh const* b, qaws_scalar const* points, unsigned int count, qaws_scalar max_distance,
	qaws_batch_executor const* executor)
{
	bv_job j;
	if (!b || (count && (!out || !points)) || !(max_distance >= 0))
		return QAWS_STATUS_INVALID_ARGUMENT;
	memset(&j, 0, sizeof j);
	j.b = b;
	j.p = points;
	j.count = count;
	j.tmax2 = max_distance > 0 ? (double)max_distance * max_distance : HUGE_VAL;
	j.kind = 0;
	j.hits = out;
	return bv_run(&j, executor);
}

qaws_status qaws_mesh_bvh_ray(qaws_mesh_hit* out, qaws_mesh_bvh const* b, qaws_scalar const* origins, qaws_scalar const* directions, unsigned int count,
	qaws_scalar tmax, int any_hit, qaws_batch_executor const* executor)
{
	bv_job j;
	if (!b || (count && (!out || !origins || !directions)) || !(tmax >= 0))
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (!((b->k == 3 && b->dim == 3) || (b->k == 2 && b->dim == 2)))
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	memset(&j, 0, sizeof j);
	j.b = b;
	j.p = origins;
	j.d = directions;
	j.count = count;
	j.tmax2 = tmax > 0 ? (double)tmax : HUGE_VAL;
	j.kind = 1;
	j.any = any_hit;
	j.hits = out;
	return bv_run(&j, executor);
}

qaws_status qaws_mesh_bvh_winding_number(qaws_scalar* out, qaws_mesh_bvh const* b, qaws_scalar const* points, unsigned int count,
	qaws_batch_executor const* executor)
{
	bv_job j;
	if (!b || (count && (!out || !points)))
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (!bv_winds(b))
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	memset(&j, 0, sizeof j);
	j.b = b;
	j.p = points;
	j.count = count;
	j.kind = 2;
	j.values = out;
	return bv_run(&j, executor);
}

qaws_status qaws_mesh_bvh_signed_distance(qaws_scalar* out, qaws_mesh_hit* out_hits, qaws_mesh_bvh const* b, qaws_scalar const* points,
	unsigned int count, qaws_batch_executor const* executor)
{
	bv_job j;
	if (!b || (count && (!out || !points)))
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (!bv_winds(b))
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	memset(&j, 0, sizeof j);
	j.b = b;
	j.p = points;
	j.count = count;
	j.kind = 3;
	j.values = out;
	j.hits = out_hits;
	return bv_run(&j, executor);
}
