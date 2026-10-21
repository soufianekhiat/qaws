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

#define BV_BINS 8
#define BV_CHUNK 64
#define BV_PI 3.14159265358979323846

typedef struct bv_node
{
	double lo[3], hi[3];
	unsigned int a, b;   /* interior: the children; leaf: the first entry of ord and the count */
	unsigned int leaf;
} bv_node;

/* per node: centre (3), area vector sum (3), first moment C_ij = sum a_i d_j (9), radius, area, second moment T_ijk = sum a_i d_j d_k (27) */
#define BV_WM 44

struct qaws_mesh_bvh
{
	unsigned int dim, k, nv, ns, nn, depth;
	double* x;            /* 3 per vertex (z 0 in 2D) */
	unsigned int* s;      /* k per simplex */
	unsigned int* el;     /* per simplex */
	unsigned int* ord;    /* the simplices in leaf order */
	bv_node* nodes;
	double beta;
	double* wm;           /* BV_WM per node when the winding number applies */
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

/* a simplex's oriented measure vector (area times unit normal; in 2D the length times the right-hand normal) and centroid */
static void bv_prim_moment(qaws_mesh_bvh const* b, unsigned int p, double* av, double* c)
{
	unsigned int const* s = b->s + (size_t)b->k * p;
	double const *A = b->x + 3 * (size_t)s[0], *B = b->x + 3 * (size_t)s[1];
	if (b->k == 3)
	{
		double const* C = b->x + 3 * (size_t)s[2];
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

/* boxes (and moments) from the leaves up */
static void bv_refit_nodes(qaws_mesh_bvh* b)
{
	unsigned int i, j, q, w = bv_winds(b) && b->wm;
	for (i = b->nn; i-- > 0;)
	{
		bv_node* nd = b->nodes + i;
		double* m = w ? b->wm + BV_WM * (size_t)i : NULL;
		bv_empty(nd->lo, nd->hi);
		if (nd->leaf)
		{
			double area = 0, cs[3] = { 0, 0, 0 }, mean[3] = { 0, 0, 0 };
			for (j = 0; j < nd->b; j++)
			{
				double lo[3], hi[3];
				bv_prim_box(b, b->ord[nd->a + j], lo, hi);
				bv_grow(nd->lo, nd->hi, lo, hi);
			}
			if (!m)
				continue;
			memset(m, 0, sizeof(double) * BV_WM);
			for (j = 0; j < nd->b; j++)
			{
				double av[3], c[3], ar;
				bv_prim_moment(b, b->ord[nd->a + j], av, c);
				ar = sqrt(bv_dot(av, av));
				area += ar;
				for (q = 0; q < 3; q++)
				{
					cs[q] += ar * c[q];
					mean[q] += c[q] / nd->b;
					m[3 + q] += av[q];
				}
			}
			for (q = 0; q < 3; q++)
				m[q] = area > 0 ? cs[q] / area : mean[q];
			for (j = 0; j < nd->b; j++)
			{
				unsigned int p = b->ord[nd->a + j], r, c2;
				double av[3], c[3], d[3];
				bv_prim_moment(b, p, av, c);
				double S[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 }, sw = b->k == 3 ? 1.0 / 12 : 1.0 / 6;
				unsigned int c3;
				bv_sub(c, m, d);
				for (r = 0; r < 3; r++)
					for (c2 = 0; c2 < 3; c2++)
						m[6 + 3 * r + c2] += av[r] * d[c2];
				for (r = 0; r < b->k; r++)
				{
					double e[3], v[3];
					bv_sub(b->x + 3 * (size_t)b->s[(size_t)b->k * p + r], m, e);
					if (sqrt(bv_dot(e, e)) > m[15])
						m[15] = sqrt(bv_dot(e, e));
					/* the simplex's own spread: its second moment about its centroid per unit measure */
					bv_sub(b->x + 3 * (size_t)b->s[(size_t)b->k * p + r], c, v);
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
		else
		{
			bv_node const *l = b->nodes + nd->a, *r = b->nodes + nd->b;
			bv_grow(nd->lo, nd->hi, l->lo, l->hi);
			bv_grow(nd->lo, nd->hi, r->lo, r->hi);
			if (!m)
				continue;
			{
				double const *ml = b->wm + BV_WM * (size_t)nd->a, *mr = b->wm + BV_WM * (size_t)nd->b;
				double area = ml[16] + mr[16], dl[3], dr[3];
				unsigned int r2, c2, c3, side;
				for (q = 0; q < 3; q++)
				{
					m[q] = area > 0 ? (ml[16] * ml[q] + mr[16] * mr[q]) / area : 0.5 * (ml[q] + mr[q]);
					m[3 + q] = ml[3 + q] + mr[3 + q];
				}
				bv_sub(ml, m, dl);
				bv_sub(mr, m, dr);
				for (r2 = 0; r2 < 3; r2++)
					for (c2 = 0; c2 < 3; c2++)
						m[6 + 3 * r2 + c2] = ml[6 + 3 * r2 + c2] + mr[6 + 3 * r2 + c2] + ml[3 + r2] * dl[c2] + mr[3 + r2] * dr[c2];
				/* T about the new centre: T + C_ij s_k + C_ik s_j + N_i s_j s_k per child, s its shift */
				for (r2 = 0; r2 < 27; r2++)
					m[17 + r2] = 0;
				for (side = 0; side < 2; side++)
				{
					double const* mc = side ? mr : ml;
					double const* s = side ? dr : dl;
					for (r2 = 0; r2 < 3; r2++)
						for (c2 = 0; c2 < 3; c2++)
							for (c3 = 0; c3 < 3; c3++)
								m[17 + 9 * r2 + 3 * c2 + c3] += mc[17 + 9 * r2 + 3 * c2 + c3] + mc[6 + 3 * r2 + c2] * s[c3] + mc[6 + 3 * r2 + c3] * s[c2]
									+ mc[3 + r2] * s[c2] * s[c3];
				}
				m[15] = ml[15] + sqrt(bv_dot(dl, dl));
				if (mr[15] + sqrt(bv_dot(dr, dr)) > m[15])
					m[15] = mr[15] + sqrt(bv_dot(dr, dr));
				m[16] = area;
			}
		}
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
	double beta = desc && desc->winding_accuracy > 0 ? (double)desc->winding_accuracy : 2.0;
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
	st = bv_build(b, leaf);
	if (st == QAWS_STATUS_OK && bv_winds(b) && b->nn)
	{
		b->wm = (double*)qaws_internal_malloc(qaws_internal_mul3((size_t)b->nn, BV_WM, sizeof(double)));
		if (!b->wm)
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
	qaws_internal_free(b->wm);
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

static double bv_box_d2(bv_node const* n, double const* p)
{
	double d = 0, t;
	unsigned int i;
	for (i = 0; i < 3; i++)
	{
		t = p[i] < n->lo[i] ? n->lo[i] - p[i] : p[i] > n->hi[i] ? p[i] - n->hi[i] : 0;
		d += t * t;
	}
	return d;
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

static void bv_hit_fill(qaws_mesh_bvh const* b, qaws_mesh_hit* h, unsigned int p, double const* w, double dist)
{
	unsigned int j, q;
	double pt[3] = { 0, 0, 0 };
	h->element = b->el[p];
	h->simplex = p;
	for (j = 0; j < 4; j++)
	{
		h->vertices[j] = j < b->k ? b->s[(size_t)b->k * p + j] : QAWS_MESH_NONE;
		h->weights[j] = j < b->k ? (qaws_scalar)w[j] : 0;
		if (j < b->k)
			for (q = 0; q < 3; q++)
				pt[q] += w[j] * b->x[3 * (size_t)h->vertices[j] + q];
	}
	for (q = 0; q < 3; q++)
		h->point[q] = (qaws_scalar)pt[q];
	h->distance = (qaws_scalar)dist;
}

static double bv_prim_closest(qaws_mesh_bvh const* b, unsigned int p, double const* q, double* w)
{
	unsigned int const* s = b->s + (size_t)b->k * p;
	double const* v[4];
	unsigned int j;
	for (j = 0; j < b->k; j++)
		v[j] = b->x + 3 * (size_t)s[j];
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

/* the closest simplex to q with squared distance under *best2 (the best so far; updated); stack: depth + 1 entries */
static unsigned int bv_closest_one(qaws_mesh_bvh const* b, double const* q, double* best2, double* bw, unsigned int* stack)
{
	unsigned int top = 0, found = QAWS_MESH_NONE, i;
	if (!b->nn || bv_box_d2(b->nodes, q) > *best2)
		return found;
	stack[top++] = 0;
	while (top)
	{
		bv_node const* n = b->nodes + stack[--top];
		if (bv_box_d2(n, q) > *best2)
			continue;
		if (n->leaf)
		{
			for (i = 0; i < n->b; i++)
			{
				double w[4], d = bv_prim_closest(b, b->ord[n->a + i], q, w);
				if (d <= *best2 && (found == QAWS_MESH_NONE || d < *best2))
				{
					*best2 = d;
					found = b->ord[n->a + i];
					memcpy(bw, w, sizeof w);
				}
			}
			if (*best2 == 0)
				break;
		}
		else
		{
			double dl = bv_box_d2(b->nodes + n->a, q), dr = bv_box_d2(b->nodes + n->b, q);
			unsigned int nearer = dl <= dr ? n->a : n->b, farther = dl <= dr ? n->b : n->a;
			if ((dl <= dr ? dr : dl) <= *best2)
				stack[top++] = farther;
			if ((dl <= dr ? dl : dr) <= *best2)
				stack[top++] = nearer;
		}
	}
	return found;
}

/* ---- rays ---- */

/* the ray's entry into the box, or HUGE_VAL when it misses it within [0, tmax] */
static double bv_box_ray(bv_node const* n, double const* o, double const* inv, double tmax)
{
	double t0 = 0, t1 = tmax;
	unsigned int i;
	for (i = 0; i < 3; i++)
	{
		double a, c;
		if (inv[i] == HUGE_VAL || inv[i] == -HUGE_VAL)
		{
			if (o[i] < n->lo[i] || o[i] > n->hi[i])
				return HUGE_VAL;
			continue;
		}
		a = (n->lo[i] - o[i]) * inv[i];
		c = (n->hi[i] - o[i]) * inv[i];
		if (a > c)
		{
			double t = a;
			a = c;
			c = t;
		}
		if (a > t0)
			t0 = a;
		if (c < t1)
			t1 = c;
		if (t0 > t1)
			return HUGE_VAL;
	}
	return t0;
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
	unsigned int* stack)
{
	double inv[3], sh[3] = { 0, 0, 0 };
	unsigned int ax[3] = { 0, 1, 2 }, top = 0, found = QAWS_MESH_NONE, i;
	for (i = 0; i < 3; i++)
		inv[i] = d[i] != 0 ? 1 / d[i] : HUGE_VAL;
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
	if (!b->nn || bv_box_ray(b->nodes, o, inv, tmax) == HUGE_VAL)
		return found;
	stack[top++] = 0;
	while (top)
	{
		bv_node const* n = b->nodes + stack[--top];
		if (bv_box_ray(n, o, inv, *t_best) == HUGE_VAL)
			continue;
		if (n->leaf)
		{
			for (i = 0; i < n->b; i++)
			{
				unsigned int p = b->ord[n->a + i];
				unsigned int const* s = b->s + (size_t)b->k * p;
				double t, w[4] = { 0, 0, 0, 0 };
				int hit = b->k == 3 ? bv_ray_tri(o, ax, sh, b->x + 3 * (size_t)s[0], b->x + 3 * (size_t)s[1], b->x + 3 * (size_t)s[2], *t_best, &t, w)
				                    : bv_ray_seg(o, d, b->x + 3 * (size_t)s[0], b->x + 3 * (size_t)s[1], *t_best, &t, w);
				if (hit && (found == QAWS_MESH_NONE || t < *t_best))
				{
					*t_best = t;
					found = p;
					memcpy(bw, w, sizeof w);
					if (any)
						return found;
				}
			}
		}
		else
		{
			double tl = bv_box_ray(b->nodes + n->a, o, inv, *t_best), tr = bv_box_ray(b->nodes + n->b, o, inv, *t_best);
			unsigned int nearer = tl <= tr ? n->a : n->b, farther = tl <= tr ? n->b : n->a;
			if ((tl <= tr ? tr : tl) != HUGE_VAL)
				stack[top++] = farther;
			if ((tl <= tr ? tl : tr) != HUGE_VAL)
				stack[top++] = nearer;
		}
	}
	return found;
}

/* ---- winding number ---- */

static double bv_prim_winding(qaws_mesh_bvh const* b, unsigned int p, double const* q)
{
	unsigned int const* s = b->s + (size_t)b->k * p;
	double A[3], B[3];
	bv_sub(b->x + 3 * (size_t)s[0], q, A);
	bv_sub(b->x + 3 * (size_t)s[1], q, B);
	if (b->k == 3)
	{
		/* Van Oosterom and Strackee's solid angle */
		double C[3], la, lb, lc, num, den;
		bv_sub(b->x + 3 * (size_t)s[2], q, C);
		la = sqrt(bv_dot(A, A));
		lb = sqrt(bv_dot(B, B));
		lc = sqrt(bv_dot(C, C));
		num = bv_det3(A, B, C);
		den = la * lb * lc + bv_dot(A, B) * lc + bv_dot(B, C) * la + bv_dot(C, A) * lb;
		return 2 * atan2(num, den) / (4 * BV_PI);
	}
	return atan2(A[0] * B[1] - A[1] * B[0], A[0] * B[0] + A[1] * B[1]) / (2 * BV_PI);
}

static double bv_winding_one(qaws_mesh_bvh const* b, double const* q, unsigned int* stack)
{
	unsigned int top = 0, i, r, c;
	double w = 0;
	if (!b->nn)
		return 0;
	stack[top++] = 0;
	while (top)
	{
		unsigned int id = stack[--top];
		bv_node const* n = b->nodes + id;
		double const* m = b->wm + BV_WM * (size_t)id;
		double rv[3], R2, R;
		bv_sub(m, q, rv);
		R2 = bv_dot(rv, rv);
		R = sqrt(R2);
		if (R > b->beta * m[15])
		{
			/* the expansion to second order: g . N + J : C + H : T / 2, g(r) = r / (4 pi R^3) (2D: r / (2 pi R^2)),
			   J and H its first and second derivatives */
			double g, a1, a2, a3, t = 0, h = 0;
			unsigned int k3;
			if (b->dim == 3)
			{
				double R3 = R2 * R, R5 = R3 * R2, R7 = R5 * R2;
				g = 1 / (4 * BV_PI * R3);
				a1 = 1 / R3;
				a2 = -3 / R5;
				a3 = 15 / R7;
				a1 /= 4 * BV_PI;
				a2 /= 4 * BV_PI;
				a3 /= 4 * BV_PI;
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
			w += g * bv_dot(rv, m + 3) + t + 0.5 * h;
			continue;
		}
		if (n->leaf)
			for (i = 0; i < n->b; i++)
				w += bv_prim_winding(b, b->ord[n->a + i], q);
		else
		{
			stack[top++] = n->b;
			stack[top++] = n->a;
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
	unsigned int* stacks;     /* depth + 1 per chunk */
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
		unsigned int* stack = j->stacks + (size_t)ch * (b->depth + 1);
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
					bv_hit_fill(b, &h, p, w, sqrt(best2));
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
					bv_hit_fill(b, &j->hits[i], p, w, t);
					j->hits[i].point[0] = (qaws_scalar)(q[0] + t * d[0]);
					j->hits[i].point[1] = (qaws_scalar)(q[1] + t * d[1]);
					j->hits[i].point[2] = (qaws_scalar)(q[2] + t * d[2]);
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
	j->stacks = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(chunks, (size_t)j->b->depth + 1, sizeof(unsigned int)));
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
