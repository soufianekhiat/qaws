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
#include "../simd/qaws_simd.h"

#if defined(__SSE2__) || defined(_M_X64) || defined(_M_AMD64) || defined(__AVX2__)
#define BV_SSE 1
#else
#define BV_SSE 0
#endif

#define BV_BINS 8
#define BV_CHUNK 64
#define BV_PI 3.14159265358979323846

typedef struct bv_node
{
	double lo[3], hi[3];
	unsigned int a, b;   /* interior: the children; leaf: the first entry of ord and the count */
	unsigned int leaf, first, cnt;   /* (first, cnt: its entries, contiguous in ord) */
} bv_node;

/* Queries run on 4-wide nodes collapsed from the binary SAH tree: per slot a float box relative to the scene's centre,
   rounded outward and widened (so the floats' rounding of rays and points never loses a hit), tested 4 at a time.
   Empty slots have lo +inf and hi -inf. */
typedef struct bv_node4
{
	float lo[3][4], hi[3][4];   /* per axis, per slot */
	unsigned int child[4];      /* inner slot: a node; leaf slot: its first entry of ord; QAWS_MESH_NONE: empty */
	unsigned int count[4];      /* 0: inner slot; else the leaf's entries */
} bv_node4;

typedef struct bv_stk
{
	unsigned int node;
	float t;                    /* the node's entry distance (rays) or squared distance (points) */
} bv_stk;

/* per 4-wide node: centre (3), area vector sum (3), first moment C_ij = sum a_i d_j (9), radius, area, second moment T_ijk = sum a_i d_j d_k (27) */
#define BV_WM 44
/* per leaf slot, in floats: centre (relative to the scene's), N, C, radius, T symmetric (18) */
#define BV_LM 34

struct qaws_mesh_bvh
{
	unsigned int dim, k, nv, ns, nn, n4, depth, leaf;
	double* x;            /* 3 per vertex (z 0 in 2D) */
	unsigned int* s;      /* k per simplex */
	unsigned int* el;     /* per simplex */
	unsigned int* ord;    /* the simplices in leaf order (entries) */
	unsigned int* ee;     /* per entry its element */
	unsigned int* ev;     /* per entry its k vertices */
	double* px;           /* per entry its k vertices' positions, contiguous */
	float* scr;           /* triangles in 3D: 9 arrays of ns + 4 floats (v0 - c, e1, e2 per axis) for the rays' screen */
	bv_node* nodes;       /* the binary tree (build only) */
	bv_node4* nodes4;
	double c[3], pad;     /* the scene's centre, the boxes' widening */
	double beta;
	double* wm;           /* BV_WM per 4-wide node when the winding number applies */
	float* lm;            /* BV_LM per slot (4 per node; leaf slots only filled) */
};

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

typedef struct bv_task
{
	unsigned int node, begin, end, depth;
} bv_task;

/* the entries of ord[begin, end) with centroid along axis below `split` first; returns the first of the rest */
static unsigned int bv_partition(unsigned int* ord, double const* cen, unsigned int axis, double split, unsigned int begin, unsigned int end)
{
	unsigned int i = begin, j = end;
	while (i < j)
	{
		if (cen[3 * (size_t)ord[i] + axis] < split)
			i++;
		else
		{
			unsigned int t = ord[--j];
			ord[j] = ord[i];
			ord[i] = t;
		}
	}
	return i;
}

/* ord[begin, end) arranged so its middle entry has the median centroid along axis (quickselect) */
static void bv_median(unsigned int* ord, double const* cen, unsigned int axis, unsigned int begin, unsigned int end)
{
	unsigned int lo = begin, hi = end - 1, mid = begin + (end - begin) / 2;
	while (lo < hi)
	{
		double pv = cen[3 * (size_t)ord[lo + (hi - lo) / 2] + axis];
		unsigned int i = lo, j = hi;
		while (i <= j)
		{
			while (cen[3 * (size_t)ord[i] + axis] < pv)
				i++;
			while (cen[3 * (size_t)ord[j] + axis] > pv)
				j--;
			if (i <= j)
			{
				unsigned int t = ord[i];
				ord[i] = ord[j];
				ord[j] = t;
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

static qaws_status bv_build(qaws_mesh_bvh* b, unsigned int leaf_size)
{
	double* box = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2((size_t)b->ns, 6 * sizeof(double)));
	double* cen = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2((size_t)b->ns, 3 * sizeof(double)));
	bv_task* stack = (bv_task*)qaws_internal_scratch_malloc(qaws_internal_mul2((size_t)b->ns + 1, sizeof(bv_task)));
	unsigned int top = 0, i, cap = b->ns ? 2 * b->ns - 1 : 0;
	b->ord = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2((size_t)b->ns + 1, sizeof(unsigned int)));
	b->nodes = (bv_node*)qaws_internal_malloc(qaws_internal_mul2((size_t)cap + 1, sizeof(bv_node)));
	if (!box || !cen || !stack || !b->ord || !b->nodes)
	{
		qaws_internal_free(box);
		qaws_internal_free(cen);
		qaws_internal_free(stack);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < b->ns; i++)
	{
		bv_prim_box(b, i, box + 6 * (size_t)i, box + 6 * (size_t)i + 3);
		cen[3 * (size_t)i] = 0.5 * (box[6 * (size_t)i] + box[6 * (size_t)i + 3]);
		cen[3 * (size_t)i + 1] = 0.5 * (box[6 * (size_t)i + 1] + box[6 * (size_t)i + 4]);
		cen[3 * (size_t)i + 2] = 0.5 * (box[6 * (size_t)i + 2] + box[6 * (size_t)i + 5]);
		b->ord[i] = i;
	}
	b->nn = 0;
	b->depth = 0;
	if (b->ns)
	{
		b->nn = 1;
		stack[top].node = 0;
		stack[top].begin = 0;
		stack[top].end = b->ns;
		stack[top++].depth = 1;
	}
	while (top)
	{
		bv_task t = stack[--top];
		bv_node* nd = b->nodes + t.node;
		nd->first = t.begin;
		nd->cnt = t.end - t.begin;
		unsigned int n = t.end - t.begin, axis, best_axis = 3, best_bin = 0, mid;
		double clo[3], chi[3], best = HUGE_VAL;
		if (t.depth > b->depth)
			b->depth = t.depth;
		bv_empty(nd->lo, nd->hi);
		bv_empty(clo, chi);
		for (i = t.begin; i < t.end; i++)
		{
			unsigned int p = b->ord[i];
			bv_grow(nd->lo, nd->hi, box + 6 * (size_t)p, box + 6 * (size_t)p + 3);
			bv_grow(clo, chi, cen + 3 * (size_t)p, cen + 3 * (size_t)p);
		}
		if (n <= leaf_size)
		{
			nd->leaf = 1;
			nd->a = t.begin;
			nd->b = n;
			continue;
		}
		/* binned SAH over each axis the centroids spread along */
		for (axis = 0; axis < b->dim; axis++)
		{
			double blo[BV_BINS][3], bhi[BV_BINS][3], ext = chi[axis] - clo[axis], scale, rlo[3], rhi[3], rarea[BV_BINS];
			unsigned int bcount[BV_BINS], rcount[BV_BINS], j, left = 0;
			double llo[3], lhi[3];
			if (!(ext > 0))
				continue;
			scale = BV_BINS / ext;
			for (j = 0; j < BV_BINS; j++)
			{
				bv_empty(blo[j], bhi[j]);
				bcount[j] = 0;
			}
			for (i = t.begin; i < t.end; i++)
			{
				unsigned int p = b->ord[i], bin = (unsigned int)((cen[3 * (size_t)p + axis] - clo[axis]) * scale);
				if (bin >= BV_BINS)
					bin = BV_BINS - 1;
				bcount[bin]++;
				bv_grow(blo[bin], bhi[bin], box + 6 * (size_t)p, box + 6 * (size_t)p + 3);
			}
			/* the right side of the plane after bin j - 1: bins j.. */
			bv_empty(rlo, rhi);
			for (j = BV_BINS - 1; j > 0; j--)
			{
				bv_grow(rlo, rhi, blo[j], bhi[j]);
				rcount[j - 1] = (j + 1 < BV_BINS ? rcount[j] : 0) + bcount[j];
				rarea[j - 1] = bv_area(rlo, rhi, b->dim);
			}
			bv_empty(llo, lhi);
			for (j = 0; j + 1 < BV_BINS; j++)
			{
				double cost;
				left += bcount[j];
				bv_grow(llo, lhi, blo[j], bhi[j]);
				if (!left || !rcount[j])
					continue;
				cost = bv_area(llo, lhi, b->dim) * left + rarea[j] * rcount[j];
				if (cost < best)
				{
					best = cost;
					best_axis = axis;
					best_bin = j;
				}
			}
		}
		mid = t.begin;
		if (best_axis < 3)
		{
			double split = clo[best_axis] + (chi[best_axis] - clo[best_axis]) * (best_bin + 1) / BV_BINS;
			mid = bv_partition(b->ord, cen, best_axis, split, t.begin, t.end);
		}
		if (mid == t.begin || mid == t.end)
		{
			/* no split by the bins (centroids together): the median along the widest spread, or a leaf when they coincide */
			double w = -1;
			for (axis = 0, best_axis = 0; axis < b->dim; axis++)
				if (chi[axis] - clo[axis] > w)
				{
					w = chi[axis] - clo[axis];
					best_axis = axis;
				}
			bv_median(b->ord, cen, best_axis, t.begin, t.end);
			mid = t.begin + n / 2;
		}
		nd->leaf = 0;
		nd->a = b->nn;
		nd->b = b->nn + 1;
		b->nn += 2;
		stack[top].node = nd->b;
		stack[top].begin = mid;
		stack[top].end = t.end;
		stack[top++].depth = t.depth + 1;
		stack[top].node = nd->a;
		stack[top].begin = t.begin;
		stack[top].end = mid;
		stack[top++].depth = t.depth + 1;
	}
	qaws_internal_free(box);
	qaws_internal_free(cen);
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

/* the moments of entries [first, first + count) into m (BV_WM) */
static void bv_leaf_moment(qaws_mesh_bvh const* b, unsigned int first, unsigned int count, double* m)
{
	double area = 0, cs[3] = { 0, 0, 0 }, mean[3] = { 0, 0, 0 };
	unsigned int j, q, r, c2, c3;
	memset(m, 0, sizeof(double) * BV_WM);
	for (j = 0; j < count; j++)
	{
		double av[3], c[3], ar;
		bv_entry_moment(b, first + j, av, c);
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
		double av[3], c[3], d[3], S[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 }, sw = b->k == 3 ? 1.0 / 12 : 1.0 / 6;
		bv_entry_moment(b, first + j, av, c);
		bv_sub(c, m, d);
		for (r = 0; r < 3; r++)
			for (c2 = 0; c2 < 3; c2++)
				m[6 + 3 * r + c2] += av[r] * d[c2];
		for (r = 0; r < b->k; r++)
		{
			double e[3], v[3];
			bv_sub(bv_pv(b, first + j, r), m, e);
			if (sqrt(bv_dot(e, e)) > m[15])
				m[15] = sqrt(bv_dot(e, e));
			/* the simplex's own spread: its second moment about its centroid per unit measure */
			bv_sub(bv_pv(b, first + j, r), c, v);
			for (c2 = 0; c2 < 3; c2++)
				for (c3 = 0; c3 < 3; c3++)
					S[3 * c2 + c3] += sw * v[c2] * v[c3];
		}
		for (r = 0; r < 3; r++)
			for (c2 = 0; c2 < 3; c2++)
				for (c3 = 0; c3 < 3; c3++)
					m[17 + 9 * r + 3 * c2 + c3] += av[r] * (d[c2] * d[c3] + S[3 * c2 + c3]);
	}
	m[16] = area;
}

/* the moments of up to 4 parts merged into m: centre by area, C and T carried to it (T + C_ij s_k + C_ik s_j + N_i s_j s_k, s the shift) */
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
			for (c2 = 0; c2 < 3; c2++)
			{
				m[6 + 3 * r + c2] += p[6 + 3 * r + c2] + p[3 + r] * s[c2];
				for (c3 = 0; c3 < 3; c3++)
					m[17 + 9 * r + 3 * c2 + c3] += p[17 + 9 * r + 3 * c2 + c3] + p[6 + 3 * r + c2] * s[c3] + p[6 + 3 * r + c3] * s[c2]
						+ p[3 + r] * s[c2] * s[c3];
			}
		rad = p[15] + sqrt(bv_dot(s, s));
		if (rad > m[15])
			m[15] = rad;
	}
	m[16] = area;
}

/* a double box to floats relative to the centre, rounded outward and widened by pad */
static float bv_down(double v)
{
	float f = (float)v;
	return (double)f > v ? nextafterf(f, -HUGE_VALF) : f;
}

static float bv_up(double v)
{
	float f = (float)v;
	return (double)f < v ? nextafterf(f, HUGE_VALF) : f;
}

static void bv_slot_set(qaws_mesh_bvh* b, bv_node4* n, unsigned int s, double const* lo, double const* hi)
{
	unsigned int a;
	for (a = 0; a < 3; a++)
	{
		n->lo[a][s] = bv_down(lo[a] - b->c[a] - b->pad);
		n->hi[a][s] = bv_up(hi[a] - b->c[a] + b->pad);
	}
}

static void bv_slot_clear(bv_node4* n, unsigned int s)
{
	unsigned int a;
	for (a = 0; a < 3; a++)
	{
		/* (lo +inf, hi -inf: every slab test, whatever the signs, and every box distance miss it) */
		n->lo[a][s] = HUGE_VALF;
		n->hi[a][s] = -HUGE_VALF;
	}
	n->child[s] = QAWS_MESH_NONE;
	n->count[s] = 0;
}

/* the binary tree collapsed into 4-wide nodes (the largest inner child opened first), parents before children; the binary nodes freed */
static qaws_status bv_collapse(qaws_mesh_bvh* b)
{
	unsigned int *stk = NULL, top = 0, i;
	double lo[3], hi[3];
	b->n4 = 0;
	b->depth = 0;
	if (!b->nn)
		return QAWS_STATUS_OK;
	b->nodes4 = (bv_node4*)qaws_internal_malloc(qaws_internal_mul2((size_t)b->nn + 1, sizeof(bv_node4)));
	stk = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3((size_t)b->nn + 1, 3, sizeof(unsigned int)));
	b->px = (double*)qaws_internal_malloc(qaws_internal_mul3((size_t)b->ns + 1, 3 * (size_t)b->k, sizeof(double)));
	if (b->k == 3 && b->dim == 3)
		b->scr = (float*)qaws_internal_zalloc(qaws_internal_mul2((size_t)b->ns + 4, 9), sizeof(float));
	b->ee = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2((size_t)b->ns + 1, sizeof(unsigned int)));
	b->ev = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3((size_t)b->ns + 1, b->k, sizeof(unsigned int)));
	if (!b->nodes4 || !stk || !b->px || !b->ee || !b->ev || (b->k == 3 && b->dim == 3 && !b->scr))
	{
		qaws_internal_free(stk);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
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
	b->n4 = 1;
	stk[0] = 0;
	stk[1] = 0;
	stk[2] = 1;
	top = 3;
	while (top)
	{
		unsigned int depth = stk[--top], n4 = stk[--top], bin = stk[--top], list[4], cnt = 0, s;
		bv_node4* n = b->nodes4 + n4;
		if (depth > b->depth)
			b->depth = depth;
		if (b->nodes[bin].cnt <= b->leaf)
			list[cnt++] = bin;
		else
		{
			list[cnt++] = b->nodes[bin].a;
			list[cnt++] = b->nodes[bin].b;
		}
		while (cnt < 4)
		{
			unsigned int best = 4;
			double ba = -1;
			for (s = 0; s < cnt; s++)
				if (b->nodes[list[s]].cnt > b->leaf && bv_area(b->nodes[list[s]].lo, b->nodes[list[s]].hi, 3) > ba)
				{
					ba = bv_area(b->nodes[list[s]].lo, b->nodes[list[s]].hi, 3);
					best = s;
				}
			if (best == 4)
				break;
			s = list[best];
			list[best] = b->nodes[s].a;
			list[cnt++] = b->nodes[s].b;
		}
		for (s = 0; s < 4; s++)
		{
			bv_node const* c;
			if (s >= cnt)
			{
				bv_slot_clear(n, s);
				continue;
			}
			c = b->nodes + list[s];
			bv_slot_set(b, n, s, c->lo, c->hi);
			if (c->cnt <= b->leaf)
			{
				n->child[s] = c->first;
				n->count[s] = c->cnt;
			}
			else
			{
				n->child[s] = b->n4++;
				n->count[s] = 0;
				stk[top++] = list[s];
				stk[top++] = n->child[s];
				stk[top++] = depth + 1;
			}
		}
	}
	qaws_internal_free(stk);
	qaws_internal_free(b->nodes);
	b->nodes = NULL;
	return QAWS_STATUS_OK;
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

static void bv_unpack_moment(qaws_mesh_bvh const* b, float const* p, double* m)
{
	unsigned int i, j, k, n = 16;
	for (i = 0; i < 3; i++)
		m[i] = (double)p[i] + b->c[i];
	for (i = 3; i < 16; i++)
		m[i] = p[i];
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			for (k = j; k < 3; k++)
				m[17 + 9 * i + 3 * j + k] = m[17 + 9 * i + 3 * k + j] = p[n++];
}

/* the leaf-ordered positions, then the boxes (and moments) from the leaves up */
static void bv_refit_nodes(qaws_mesh_bvh* b)
{
	unsigned int i, j, s, a, w = bv_winds(b) && b->wm;
	double leafm[4][BV_WM];
	for (i = 0; i < b->ns; i++)
		for (j = 0; j < b->k; j++)
			memcpy(b->px + 3 * ((size_t)b->k * i + j), b->x + 3 * (size_t)b->s[(size_t)b->k * b->ord[i] + j], 3 * sizeof(double));
	if (b->scr)
		for (i = 0; i < b->ns; i++)
		{
			size_t st = (size_t)b->ns + 4;
			double const* v = b->px + 9 * (size_t)i;
			for (a = 0; a < 3; a++)
			{
				b->scr[a * st + i] = (float)(v[a] - b->c[a]);
				b->scr[(3 + a) * st + i] = (float)(v[3 + a] - v[a]);
				b->scr[(6 + a) * st + i] = (float)(v[6 + a] - v[a]);
			}
		}
	for (i = b->n4; i-- > 0;)
	{
		bv_node4* n = b->nodes4 + i;
		double const* part[4];
		unsigned int np = 0;
		for (s = 0; s < 4; s++)
		{
			if (n->child[s] == QAWS_MESH_NONE)
				continue;
			if (n->count[s])
			{
				double lo[3], hi[3];
				bv_empty(lo, hi);
				for (j = 0; j < n->count[s] * b->k; j++)
					bv_grow(lo, hi, b->px + 3 * ((size_t)b->k * n->child[s] + j), b->px + 3 * ((size_t)b->k * n->child[s] + j));
				bv_slot_set(b, n, s, lo, hi);
				if (w)
				{
					bv_leaf_moment(b, n->child[s], n->count[s], leafm[s]);
					bv_pack_moment(b, leafm[s], b->lm + BV_LM * (4 * (size_t)i + s));
					part[np++] = leafm[s];
				}
			}
			else
			{
				/* the child's slots' union (its floats already rounded outward) */
				bv_node4 const* c = b->nodes4 + n->child[s];
				unsigned int t;
				for (a = 0; a < 3; a++)
				{
					float l = HUGE_VALF, h = -HUGE_VALF;
					for (t = 0; t < 4; t++)
						if (c->child[t] != QAWS_MESH_NONE)
						{
							l = c->lo[a][t] < l ? c->lo[a][t] : l;
							h = c->hi[a][t] > h ? c->hi[a][t] : h;
						}
					n->lo[a][s] = l;
					n->hi[a][s] = h;
				}
				if (w)
					part[np++] = b->wm + BV_WM * (size_t)n->child[s];
			}
		}
		if (w)
			bv_merge_moments(part, np, b->wm + BV_WM * (size_t)i);
	}
}

/* ---- creation ---- */

/* takes x, s, el (freed on failure) */
static qaws_status bv_make(qaws_mesh_bvh** out, double* x, unsigned int nv, unsigned int dim, unsigned int* s, unsigned int k, unsigned int ns,
	unsigned int* el, qaws_mesh_bvh_desc const* desc)
{
	qaws_mesh_bvh* b = (qaws_mesh_bvh*)qaws_internal_zalloc(1, sizeof *b);
	qaws_status st;
	unsigned int leaf = desc && desc->leaf_size ? desc->leaf_size : 4;
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
	b->leaf = leaf;
	st = bv_build(b, leaf);
	if (st == QAWS_STATUS_OK)
		st = bv_collapse(b);
	if (st == QAWS_STATUS_OK && bv_winds(b) && b->n4)
	{
		b->wm = (double*)qaws_internal_malloc(qaws_internal_mul3((size_t)b->n4, BV_WM, sizeof(double)));
		b->lm = (float*)qaws_internal_malloc(qaws_internal_mul3((size_t)b->n4, 4 * BV_LM, sizeof(float)));
		if (!b->wm || !b->lm)
			st = QAWS_STATUS_ALLOCATION_FAILURE;
	}
	if (st != QAWS_STATUS_OK)
	{
		qaws_mesh_bvh_destroy(b);
		return st;
	}
	bv_refit_nodes(b);
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
	qaws_internal_free(b->nodes4);
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
	bv_refit_nodes(b);
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

/* a node's 4 slots against a point (squared distances) or a ray (entry distances, misses at +inf); same arithmetic with or without SSE */
static void bv_slots_d2(bv_node4 const* n, float const* q, float* out)
{
#if BV_SSE
	__m128 d2 = _mm_setzero_ps(), z = _mm_setzero_ps();
	unsigned int a;
	for (a = 0; a < 3; a++)
	{
		__m128 qa = _mm_set1_ps(q[a]), dl = _mm_sub_ps(_mm_loadu_ps(n->lo[a]), qa), dh = _mm_sub_ps(qa, _mm_loadu_ps(n->hi[a]));
		__m128 m = _mm_max_ps(_mm_max_ps(dl, dh), z);
		d2 = _mm_add_ps(d2, _mm_mul_ps(m, m));
	}
	_mm_storeu_ps(out, d2);
#else
	unsigned int a, s;
	for (s = 0; s < 4; s++)
	{
		float d2 = 0;
		for (a = 0; a < 3; a++)
		{
			float dl = n->lo[a][s] - q[a], dh = q[a] - n->hi[a][s], m = dl > dh ? dl : dh;
			m = m > 0 ? m : 0;
			d2 = d2 + m * m;
		}
		out[s] = d2;
	}
#endif
}

/* a node's 4 slots against a ray: the entry distances into tn and the hit slots as a mask; near and far planes picked by the
   direction's signs (nf: per axis 0 when it is positive), so no min or max per axis. Same arithmetic with or without SSE. */
static unsigned int bv_slots_ray(bv_node4 const* n, float const* o, float const* id, unsigned int const* nf, float tmax, float* tn)
{
#if BV_SSE
	__m128 t0 = _mm_setzero_ps(), t1 = _mm_set1_ps(tmax);
	unsigned int a;
	for (a = 0; a < 3; a++)
	{
		float const* nr = nf[a] ? n->hi[a] : n->lo[a];
		float const* fr = nf[a] ? n->lo[a] : n->hi[a];
		__m128 oa = _mm_set1_ps(o[a]), ia = _mm_set1_ps(id[a]);
		t0 = _mm_max_ps(t0, _mm_mul_ps(_mm_sub_ps(_mm_loadu_ps(nr), oa), ia));
		t1 = _mm_min_ps(t1, _mm_mul_ps(_mm_sub_ps(_mm_loadu_ps(fr), oa), ia));
	}
	_mm_storeu_ps(tn, t0);
	return (unsigned int)_mm_movemask_ps(_mm_cmple_ps(t0, t1));
#else
	unsigned int a, s, m = 0;
	for (s = 0; s < 4; s++)
	{
		float t0 = 0, t1 = tmax;
		for (a = 0; a < 3; a++)
		{
			float tl = ((nf[a] ? n->hi[a][s] : n->lo[a][s]) - o[a]) * id[a], th = ((nf[a] ? n->lo[a][s] : n->hi[a][s]) - o[a]) * id[a];
			t0 = t0 > tl ? t0 : tl;
			t1 = t1 < th ? t1 : th;
		}
		tn[s] = t0;
		m |= (unsigned int)(t0 <= t1) << s;
	}
	return m;
#endif
}

/* the slots with a finite key, nearest first */
static unsigned int bv_order(float const* key, unsigned int* slot)
{
	unsigned int n = 0, s, i;
	for (s = 0; s < 4; s++)
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

/* the closest entry to q with squared distance within *best2 (updated); its weights into bw. Returns the simplex. */
static unsigned int bv_closest_one(qaws_mesh_bvh const* b, double const* q, double* best2, double* bw, bv_stk* stk)
{
	unsigned int top = 0, found = QAWS_MESH_NONE, i, s;
	float qf[3] = { (float)(q[0] - b->c[0]), (float)(q[1] - b->c[1]), (float)(q[2] - b->c[2]) }, bf = bv_fbound(*best2);
	if (!b->n4)
		return found;
	stk[top].node = 0;
	stk[top++].t = 0;
	while (top)
	{
		bv_stk e = stk[--top];
		bv_node4 const* n;
		float key[4];
		unsigned int slot[4], ns;
		if (e.t > bf)
			continue;
		n = b->nodes4 + e.node;
		bv_slots_d2(n, qf, key);
		for (s = 0; s < 4; s++)
			if (key[s] > bf)
				key[s] = HUGE_VALF;
		ns = bv_order(key, slot);
		for (i = 0; i < ns; i++)
		{
			unsigned int sl = slot[i], c;
			if (!n->count[sl] || key[sl] > bf)
				continue;
			for (c = 0; c < n->count[sl]; c++)
			{
				double w[4], d = bv_entry_closest(b, n->child[sl] + c, q, w);
				if (d <= *best2 && (found == QAWS_MESH_NONE || d < *best2))
				{
					*best2 = d;
					found = n->child[sl] + c;
					memcpy(bw, w, sizeof w);
					bf = bv_fbound(d);
				}
			}
			if (*best2 == 0)
				return found;
		}
		for (i = ns; i-- > 0;)
			if (!n->count[slot[i]] && key[slot[i]] <= bf)
			{
				stk[top].node = n->child[slot[i]];
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
	size_t st = (size_t)b->ns + 4;
#if BV_SSE
	__m128 sign = _mm_set1_ps(-0.0f), tol = _mm_set1_ps(1e-5f), z = _mm_setzero_ps();
	__m128 dx = _mm_set1_ps(d[0]), dy = _mm_set1_ps(d[1]), dz = _mm_set1_ps(d[2]);
	__m128 ox = _mm_set1_ps(o[0]), oy = _mm_set1_ps(o[1]), oz = _mm_set1_ps(o[2]);
	__m128 da = _mm_set1_ps(fabsf(d[0]) + fabsf(d[1]) + fabsf(d[2])), oa = _mm_set1_ps(fabsf(o[0]) + fabsf(o[1]) + fabsf(o[2]));
	__m128 vx = _mm_loadu_ps(tf + first), vy = _mm_loadu_ps(tf + st + first), vz = _mm_loadu_ps(tf + 2 * st + first);
	__m128 ax = _mm_loadu_ps(tf + 3 * st + first), ay = _mm_loadu_ps(tf + 4 * st + first), az = _mm_loadu_ps(tf + 5 * st + first);
	__m128 bx = _mm_loadu_ps(tf + 6 * st + first), by = _mm_loadu_ps(tf + 7 * st + first), bz = _mm_loadu_ps(tf + 8 * st + first);
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
	__m128 e1n = _mm_add_ps(_mm_add_ps(_mm_andnot_ps(sign, ax), _mm_andnot_ps(sign, ay)), _mm_andnot_ps(sign, az));
	__m128 e2n = _mm_add_ps(_mm_add_ps(_mm_andnot_ps(sign, bx), _mm_andnot_ps(sign, by)), _mm_andnot_ps(sign, bz));
	__m128 vn = _mm_add_ps(_mm_add_ps(_mm_andnot_ps(sign, vx), _mm_andnot_ps(sign, vy)), _mm_andnot_ps(sign, vz));
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
		size_t e = first + l;
		float vx = tf[e], vy = tf[st + e], vz = tf[2 * st + e], ax = tf[3 * st + e], ay = tf[4 * st + e], az = tf[5 * st + e];
		float bx = tf[6 * st + e], by = tf[7 * st + e], bz = tf[8 * st + e];
		float px = d[1] * bz - d[2] * by, py = d[2] * bx - d[0] * bz, pz = d[0] * by - d[1] * bx;
		float det = ax * px + ay * py + az * pz, tx = o[0] - vx, ty = o[1] - vy, tz = o[2] - vz;
		float u = tx * px + ty * py + tz * pz;
		float qx = ty * az - tz * ay, qy = tz * ax - tx * az, qz = tx * ay - ty * ax;
		float v = d[0] * qx + d[1] * qy + d[2] * qz, t = bx * qx + by * qy + bz * qz;
		float tn = fabsf(tx) + fabsf(ty) + fabsf(tz) + oa + (fabsf(vx) + fabsf(vy) + fabsf(vz));
		float e1n = fabsf(ax) + fabsf(ay) + fabsf(az), e2n = fabsf(bx) + fabsf(by) + fabsf(bz);
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
	unsigned int ax[3] = { 0, 1, 2 }, nf[3], top = 0, found = QAWS_MESH_NONE, i, cand = 0, node = 0;
	float of[3], id[3], df[3], tf;
	for (i = 0; i < 3; i++)
	{
		of[i] = (float)(o[i] - b->c[i]);
		df[i] = (float)d[i];
		id[i] = d[i] != 0 ? (float)(1 / d[i]) : 1e38f;
		nf[i] = id[i] < 0;
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
	if (!b->n4)
		return found;
	for (;;)
	{
		bv_node4 const* n = b->nodes4 + node;
		float tn[4], near = HUGE_VALF;
		unsigned int m = bv_slots_ray(n, of, id, nf, tf, tn), s, next = QAWS_MESH_NONE;
		/* the leaves now, nearest first */
		while (m)
		{
			unsigned int sl = 4, c;
			float best = HUGE_VALF;
			for (s = 0; s < 4; s++)
				if ((m >> s & 1) && n->count[s] && tn[s] <= best)
				{
					best = tn[s];
					sl = s;
				}
			if (sl == 4)
				break;
			m &= ~(1u << sl);
			if (tn[sl] > tf)
				continue;
			for (c = 0; c < n->count[sl]; c++)
			{
				unsigned int en = n->child[sl] + c;
				double t, w[4] = { 0, 0, 0, 0 };
				int hit;
				if (b->scr && !(c & 3))
					cand = bv_screen(b, en, n->count[sl] - c < 4 ? n->count[sl] - c : 4, of, df, tf);
				if (b->scr && !(cand >> (c & 3) & 1))
					continue;
				hit = b->k == 3 ? bv_ray_tri(o, ax, sh, bv_pv(b, en, 0), bv_pv(b, en, 1), bv_pv(b, en, 2), *t_best, &t, w)
				                : bv_ray_seg(o, d, bv_pv(b, en, 0), bv_pv(b, en, 1), *t_best, &t, w);
				if (hit && (found == QAWS_MESH_NONE || t < *t_best))
				{
					*t_best = t;
					found = en;
					memcpy(bw, w, sizeof w);
					tf = bv_fbound(t);
					if (any)
						return found;
				}
			}
		}
		/* the inner slots: on into the nearest, the others pushed */
		for (s = 0; s < 4; s++)
			if ((m >> s & 1) && tn[s] <= tf)
			{
				if (tn[s] < near)
				{
					if (next != QAWS_MESH_NONE)
					{
						stk[top].node = next;
						stk[top++].t = near;
					}
					next = n->child[s];
					near = tn[s];
				}
				else
				{
					stk[top].node = n->child[s];
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

/* the expansion of moments m at q when q is far from them (beyond beta times their radius): adds it to *w and returns 1 */
static int bv_far(qaws_mesh_bvh const* b, double const* m, double const* q, double* w)
{
	double rv[3], R2, R, g, a1, a2, a3, t = 0, h = 0;
	unsigned int r, c, k3;
	bv_sub(m, q, rv);
	R2 = bv_dot(rv, rv);
	R = sqrt(R2);
	if (!(R > b->beta * m[15]))
		return 0;
	/* the expansion to second order: g . N + J : C + H : T / 2, g(r) = r / (4 pi R^3) (2D: r / (2 pi R^2)),
	   J and H its first and second derivatives */
	if (b->dim == 3)
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
	/* J_ij = a1 d_ij + a2 r_i r_j; H_ijk = a2 (d_ij r_k + d_ik r_j + d_jk r_i) + a3 r_i r_j r_k */
	for (r = 0; r < 3; r++)
		for (c = 0; c < 3; c++)
		{
			t += ((r == c ? a1 : 0) + a2 * rv[r] * rv[c]) * m[6 + 3 * r + c];
			for (k3 = 0; k3 < 3; k3++)
				h += (a2 * ((r == c ? rv[k3] : 0) + (r == k3 ? rv[c] : 0) + (c == k3 ? rv[r] : 0)) + a3 * rv[r] * rv[c] * rv[k3])
					* m[17 + 9 * r + 3 * c + k3];
		}
	*w += g * bv_dot(rv, m + 3) + t + 0.5 * h;
	return 1;
}

static double bv_winding_one(qaws_mesh_bvh const* b, double const* q, bv_stk* stk)
{
	unsigned int top = 0, i, s;
	double w = 0;
	if (!b->n4)
		return 0;
	stk[top++].node = 0;
	while (top)
	{
		unsigned int id = stk[--top].node;
		bv_node4 const* n = b->nodes4 + id;
		if (bv_far(b, b->wm + BV_WM * (size_t)id, q, &w))
			continue;
		for (s = 0; s < 4; s++)
		{
			if (n->child[s] == QAWS_MESH_NONE)
				continue;
			if (n->count[s])
			{
				double m[BV_WM];
				bv_unpack_moment(b, b->lm + BV_LM * (4 * (size_t)id + s), m);
				if (!bv_far(b, m, q, &w))
					for (i = 0; i < n->count[s]; i++)
						w += bv_entry_winding(b, n->child[s] + i, q);
			}
			else
				stk[top++].node = n->child[s];
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
	bv_stk* stacks;           /* 3 depth + 2 per chunk */
} bv_job;

static void bv_read(qaws_mesh_bvh const* b, qaws_scalar const* p, unsigned int i, double* q)
{
	q[0] = (double)p[(size_t)b->dim * i];
	q[1] = (double)p[(size_t)b->dim * i + 1];
	q[2] = b->dim == 3 ? (double)p[(size_t)b->dim * i + 2] : 0.0;
}

static void bv_task_run(void* ctx, unsigned int begin, unsigned int end)
{
	bv_job const* j = (bv_job const*)ctx;
	qaws_mesh_bvh const* b = j->b;
	unsigned int ch, i;
	for (ch = begin; ch < end; ch++)
	{
		bv_stk* stack = j->stacks + (size_t)ch * (3 * (size_t)b->depth + 2);
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
	j->stacks = (bv_stk*)qaws_internal_scratch_malloc(qaws_internal_mul3(chunks, 3 * (size_t)j->b->depth + 2, sizeof(bv_stk)));
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
