/*
 * Curved volume elements (qaws_volume_curve.h): quadratic nodes numbered
 * from the cells' edges and faces, the border's projected onto a mesh
 * reference, the displacement spread inside, elements relaxed toward
 * straight until their sampled Jacobians pass.
 */

#include "../internal/qaws_internal_memory.h"
#include "../qaws_volume_curve.h"
#include "qaws_mesh_internal.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define VCU_ROUNDS 12u   /* halvings before an element is made straight */

static unsigned char const vcu_tet_edge[6][2] = { { 0, 1 }, { 1, 2 }, { 2, 0 }, { 0, 3 }, { 1, 3 }, { 2, 3 } };
static unsigned char const vcu_hex_edge[12][2] = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 }, { 4, 5 }, { 5, 6 }, { 6, 7 }, { 7, 4 }, { 0, 4 }, { 1, 5 },
	{ 2, 6 }, { 3, 7 } };
static unsigned char const vcu_hex_face[6][4] = { { 0, 3, 7, 4 }, { 1, 2, 6, 5 }, { 0, 1, 5, 4 }, { 3, 2, 6, 7 }, { 0, 1, 2, 3 }, { 4, 5, 6, 7 } };
/* the triquadratic hexahedron's nodes in its parameters (-1, 0, 1) */
static signed char const vcu_hex_at[27][3] = { { -1, -1, -1 }, { 1, -1, -1 }, { 1, 1, -1 }, { -1, 1, -1 }, { -1, -1, 1 }, { 1, -1, 1 }, { 1, 1, 1 }, { -1, 1, 1 },
	{ 0, -1, -1 }, { 1, 0, -1 }, { 0, 1, -1 }, { -1, 0, -1 }, { 0, -1, 1 }, { 1, 0, 1 }, { 0, 1, 1 }, { -1, 0, 1 }, { -1, -1, 0 }, { 1, -1, 0 }, { 1, 1, 0 },
	{ -1, 1, 0 }, { -1, 0, 0 }, { 1, 0, 0 }, { 0, -1, 0 }, { 0, 1, 0 }, { 0, 0, -1 }, { 0, 0, 1 }, { 0, 0, 0 } };

typedef struct vcu_key
{
	unsigned int k[4], slot;   /* sorted vertices (an edge's last two HM_NONE); the element node slot */
} vcu_key;

static int vcu_key_cmp(void const* p, void const* q)
{
	vcu_key const *a = (vcu_key const*)p, *b = (vcu_key const*)q;
	unsigned int i;
	for (i = 0; i < 4; i++)
		if (a->k[i] != b->k[i])
			return a->k[i] < b->k[i] ? -1 : 1;
	return 0;
}

static void vcu_sort4(unsigned int* s, unsigned int n)
{
	unsigned int i, j;
	for (i = 1; i < n; i++)
		for (j = i; j > 0 && s[j - 1] > s[j]; j--)
		{
			unsigned int t = s[j];
			s[j] = s[j - 1];
			s[j - 1] = t;
		}
}

static double vcu_det(double const* a, double const* b, double const* c)
{
	return a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
}

/* a quadratic tetrahedron's Jacobian determinant at barycentric l (x: its 10 nodes) */
static double vcu_tet_det(double const (*x)[3], double const* l)
{
	double g[10][4], col[3][3];
	unsigned int n, k, d;
	memset(g, 0, sizeof g);
	for (k = 0; k < 4; k++) g[k][k] = 4 * l[k] - 1;
	for (n = 0; n < 6; n++)
	{
		unsigned int i = vcu_tet_edge[n][0], j = vcu_tet_edge[n][1];
		g[4 + n][i] = 4 * l[j];
		g[4 + n][j] = 4 * l[i];
	}
	for (k = 0; k < 3; k++)
		for (d = 0; d < 3; d++)
		{
			double s = 0;
			for (n = 0; n < 10; n++) s += x[n][d] * (g[n][k + 1] - g[n][0]);
			col[k][d] = s;
		}
	return vcu_det(col[0], col[1], col[2]);
}

static void vcu_lag(double s, double* v, double* dv)
{
	v[0] = 0.5 * s * (s - 1); v[1] = 1 - s * s; v[2] = 0.5 * s * (s + 1);
	dv[0] = s - 0.5; dv[1] = -2 * s; dv[2] = s + 0.5;
}

/* a triquadratic hexahedron's Jacobian determinant at (s, t, u) */
static double vcu_hex_det(double const (*x)[3], double s, double t, double u)
{
	double a[3], da[3], b[3], db[3], c[3], dc[3], col[3][3];
	unsigned int n, d;
	vcu_lag(s, a, da);
	vcu_lag(t, b, db);
	vcu_lag(u, c, dc);
	memset(col, 0, sizeof col);
	for (n = 0; n < 27; n++)
	{
		unsigned int i = (unsigned int)(vcu_hex_at[n][0] + 1), j = (unsigned int)(vcu_hex_at[n][1] + 1), k = (unsigned int)(vcu_hex_at[n][2] + 1);
		double w0 = da[i] * b[j] * c[k], w1 = a[i] * db[j] * c[k], w2 = a[i] * b[j] * dc[k];
		for (d = 0; d < 3; d++)
		{
			col[0][d] += w0 * x[n][d];
			col[1][d] += w1 * x[n][d];
			col[2][d] += w2 * x[n][d];
		}
	}
	return vcu_det(col[0], col[1], col[2]);
}

/* an element's smallest sampled Jacobian over its straight self's (nodes en, positions X straight + D) */
static double vcu_ratio(unsigned int const* en, unsigned int nn, double const* X, double const* D)
{
	double xc[27][3], xs[27][3], r = HUGE_VAL;
	unsigned int n, d;
	for (n = 0; n < nn; n++)
		for (d = 0; d < 3; d++)
		{
			xs[n][d] = X[3 * (size_t)en[n] + d];
			xc[n][d] = xs[n][d] + D[3 * (size_t)en[n] + d];
		}
	if (nn == 10)
	{
		static double const samp[15][4] = { { 1, 0, 0, 0 }, { 0, 1, 0, 0 }, { 0, 0, 1, 0 }, { 0, 0, 0, 1 }, { 0.5, 0.5, 0, 0 }, { 0, 0.5, 0.5, 0 },
			{ 0.5, 0, 0.5, 0 }, { 0.5, 0, 0, 0.5 }, { 0, 0.5, 0, 0.5 }, { 0, 0, 0.5, 0.5 }, { 0, 1.0 / 3, 1.0 / 3, 1.0 / 3 }, { 1.0 / 3, 0, 1.0 / 3, 1.0 / 3 },
			{ 1.0 / 3, 1.0 / 3, 0, 1.0 / 3 }, { 1.0 / 3, 1.0 / 3, 1.0 / 3, 0 }, { 0.25, 0.25, 0.25, 0.25 } };
		for (n = 0; n < 15; n++)
		{
			double ds = vcu_tet_det((double const(*)[3])xs, samp[n]), dc = vcu_tet_det((double const(*)[3])xc, samp[n]);
			r = fmin(r, ds > 0 ? dc / ds : (dc > 0 ? 1 : -1));
		}
	}
	else
		for (n = 0; n < 27; n++)
		{
			double s = vcu_hex_at[n][0], t = vcu_hex_at[n][1], u = vcu_hex_at[n][2];
			double ds = vcu_hex_det((double const(*)[3])xs, s, t, u), dc = vcu_hex_det((double const(*)[3])xc, s, t, u);
			r = fmin(r, ds > 0 ? dc / ds : (dc > 0 ? 1 : -1));
		}
	return r;
}

static int vcu_kinds(qaws_volume_mesh const* mesh, unsigned int* nt, unsigned int* nh)
{
	unsigned int c;
	*nt = *nh = 0;
	for (c = 0; c < mesh->nc; c++)
	{
		if (mesh->type[c] == QAWS_CELL_TETRAHEDRON) (*nt)++;
		else if (mesh->type[c] == QAWS_CELL_HEXAHEDRON) (*nh)++;
		else return 0;
	}
	return 1;
}

/* the local edge slot of corners i and j in a cell of its type (the element node index) */
static unsigned int vcu_edge_slot(unsigned int type, unsigned int i, unsigned int j)
{
	unsigned int e;
	if (type == QAWS_CELL_TETRAHEDRON)
	{
		for (e = 0; e < 6; e++)
			if ((vcu_tet_edge[e][0] == i && vcu_tet_edge[e][1] == j) || (vcu_tet_edge[e][0] == j && vcu_tet_edge[e][1] == i))
				return 4 + e;
	}
	else
		for (e = 0; e < 12; e++)
			if ((vcu_hex_edge[e][0] == i && vcu_hex_edge[e][1] == j) || (vcu_hex_edge[e][0] == j && vcu_hex_edge[e][1] == i))
				return 8 + e;
	return HM_NONE;
}

qaws_status qaws_volume_mesh_quadratic_node_counts(unsigned int* out_node_count, unsigned int* out_element_node_count, qaws_volume_mesh const* mesh)
{
	unsigned int nt, nh, c, k, ne = 0, nf = 0, nk = 0;
	vcu_key* key;
	if (!out_node_count || !out_element_node_count || !mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (!vcu_kinds(mesh, &nt, &nh))
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	key = (vcu_key*)qaws_internal_scratch_malloc(sizeof(vcu_key) * (6 * (size_t)nt + 18 * (size_t)nh + 1));
	if (!key)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (c = 0; c < mesh->nc; c++)
	{
		unsigned int const* v = mesh->cv + mesh->c0[c];
		int tet = mesh->type[c] == QAWS_CELL_TETRAHEDRON;
		for (k = 0; k < (tet ? 6u : 12u); k++)
		{
			key[nk].k[0] = v[tet ? vcu_tet_edge[k][0] : vcu_hex_edge[k][0]];
			key[nk].k[1] = v[tet ? vcu_tet_edge[k][1] : vcu_hex_edge[k][1]];
			key[nk].k[2] = key[nk].k[3] = HM_NONE;
			vcu_sort4(key[nk].k, 2);
			nk++;
		}
		for (k = 0; !tet && k < 6; k++)
		{
			unsigned int i;
			for (i = 0; i < 4; i++) key[nk].k[i] = v[vcu_hex_face[k][i]];
			vcu_sort4(key[nk].k, 4);
			nk++;
		}
	}
	qsort(key, nk, sizeof(vcu_key), vcu_key_cmp);
	for (k = 0; k < nk; k++)
		if (k == 0 || vcu_key_cmp(&key[k], &key[k - 1]))
		{
			if (key[k].k[2] == HM_NONE) ne++;
			else nf++;
		}
	qaws_internal_free(key);
	*out_node_count = mesh->nv + ne + nf + nh;
	*out_element_node_count = 10 * nt + 27 * nh;
	return QAWS_STATUS_OK;
}

qaws_status qaws_volume_mesh_quadratic_nodes(qaws_scalar* out_positions, unsigned int* out_element_nodes, qaws_volume_curve_report* report,
	qaws_volume_mesh const* mesh, qaws_mesh_reference const* reference, qaws_volume_curve_desc const* desc)
{
	unsigned int nt, nh, c, k, i, d, nk = 0, nn, nen, ne = 0, nf = 0, nb = 0, *eo = NULL, *en = NULL;
	unsigned int *ea = NULL, *fc = NULL, *bord = NULL;   /* per edge node its ends; per face node its 4 corners and 4 edge nodes; border nodes */
	unsigned char *isb = NULL, *flag = NULL;
	double *X = NULL, *D = NULL, *acc = NULL, minj = desc && desc->min_jacobian > 0 ? (double)desc->min_jacobian : 0.1;
	unsigned int sweeps = desc && desc->smoothing ? desc->smoothing : 20, *cnt = NULL;
	vcu_key* key = NULL;
	qaws_volume_curve_report rep;
	qaws_status st;
	memset(&rep, 0, sizeof rep);
	if (!out_positions || !out_element_nodes || !mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	st = qaws_volume_mesh_quadratic_node_counts(&nn, &nen, mesh);
	if (st != QAWS_STATUS_OK)
		return st;
	vcu_kinds(mesh, &nt, &nh);
	key = (vcu_key*)qaws_internal_scratch_malloc(sizeof(vcu_key) * (6 * (size_t)nt + 18 * (size_t)nh + 1));
	eo = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * ((size_t)mesh->nc + 1));
	en = out_element_nodes;
	ea = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * (2 * (size_t)nn + 1));
	fc = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * (8 * (size_t)nn + 1));
	isb = (unsigned char*)qaws_internal_scratch_zalloc((size_t)nn + 1, 1);
	flag = (unsigned char*)qaws_internal_scratch_zalloc((size_t)mesh->nc + 1, 1);
	X = (double*)qaws_internal_scratch_zalloc(3 * (size_t)nn + 1, sizeof(double));
	D = (double*)qaws_internal_scratch_zalloc(3 * (size_t)nn + 1, sizeof(double));
	acc = (double*)qaws_internal_scratch_zalloc(3 * (size_t)nn + 1, sizeof(double));
	cnt = (unsigned int*)qaws_internal_scratch_zalloc((size_t)nn + 1, sizeof(unsigned int));
	bord = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * ((size_t)nn + 1));
	if (!key || !eo || !ea || !fc || !isb || !flag || !X || !D || !acc || !cnt || !bord)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}

	/* the numbering: corners, then edges and faces by their sorted vertices, then the hexahedra's centres */
	eo[0] = 0;
	for (c = 0; c < mesh->nc; c++)
	{
		unsigned int const* v = mesh->cv + mesh->c0[c];
		int tet = mesh->type[c] == QAWS_CELL_TETRAHEDRON;
		eo[c + 1] = eo[c] + (tet ? 10u : 27u);
		for (k = 0; k < (tet ? 4u : 8u); k++) en[eo[c] + k] = v[k];
		for (k = 0; k < (tet ? 6u : 12u); k++)
		{
			key[nk].k[0] = v[tet ? vcu_tet_edge[k][0] : vcu_hex_edge[k][0]];
			key[nk].k[1] = v[tet ? vcu_tet_edge[k][1] : vcu_hex_edge[k][1]];
			key[nk].k[2] = key[nk].k[3] = HM_NONE;
			vcu_sort4(key[nk].k, 2);
			key[nk++].slot = eo[c] + (tet ? 4u : 8u) + k;
		}
		for (k = 0; !tet && k < 6; k++)
		{
			for (i = 0; i < 4; i++) key[nk].k[i] = v[vcu_hex_face[k][i]];
			vcu_sort4(key[nk].k, 4);
			key[nk++].slot = eo[c] + 20 + k;
		}
	}
	qsort(key, nk, sizeof(vcu_key), vcu_key_cmp);
	{
		unsigned int id = HM_NONE;
		for (k = 0; k < nk; k++)
		{
			if (k == 0 || vcu_key_cmp(&key[k], &key[k - 1]))
				id = key[k].k[2] == HM_NONE ? mesh->nv + ne++ : HM_NONE;
			if (key[k].k[2] == HM_NONE)
			{
				en[key[k].slot] = id;
				ea[2 * (size_t)id] = key[k].k[0];
				ea[2 * (size_t)id + 1] = key[k].k[1];
			}
		}
		for (k = 0; k < nk; k++)
		{
			if (key[k].k[2] == HM_NONE)
				continue;
			if (k == 0 || vcu_key_cmp(&key[k], &key[k - 1]))
				id = mesh->nv + ne + nf++;
			en[key[k].slot] = id;
		}
	}
	for (c = 0, i = 0; c < mesh->nc; c++)
		if (mesh->type[c] == QAWS_CELL_HEXAHEDRON)
		{
			en[eo[c] + 26] = mesh->nv + ne + nf + i++;
			/* each face node's corners and edge nodes (from the first hexahedron numbering it) */
			for (k = 0; k < 6; k++)
			{
				unsigned int fid = en[eo[c] + 20 + k], *f = fc + 8 * (size_t)fid;
				for (d = 0; d < 4; d++)
				{
					f[d] = en[eo[c] + vcu_hex_face[k][d]];
					f[4 + d] = en[eo[c] + vcu_edge_slot(QAWS_CELL_HEXAHEDRON, vcu_hex_face[k][d], vcu_hex_face[k][(d + 1) % 4])];
				}
			}
		}

	/* the straight points */
	for (i = 0; i < mesh->nv; i++)
		for (d = 0; d < 3; d++) X[3 * (size_t)i + d] = mesh->x[3 * (size_t)i + d];
	for (i = 0; i < ne; i++)
		for (d = 0; d < 3; d++)
			X[3 * ((size_t)mesh->nv + i) + d] = 0.5 * (mesh->x[3 * (size_t)ea[2 * ((size_t)mesh->nv + i)] + d] + mesh->x[3 * (size_t)ea[2 * ((size_t)mesh->nv + i) + 1] + d]);
	for (i = 0; i < nf; i++)
	{
		unsigned int const* f = fc + 8 * ((size_t)mesh->nv + ne + i);
		for (d = 0; d < 3; d++)
			X[3 * ((size_t)mesh->nv + ne + i) + d] = 0.25 * (mesh->x[3 * (size_t)f[0] + d] + mesh->x[3 * (size_t)f[1] + d] + mesh->x[3 * (size_t)f[2] + d] + mesh->x[3 * (size_t)f[3] + d]);
	}
	for (c = 0; c < mesh->nc; c++)
		if (mesh->type[c] == QAWS_CELL_HEXAHEDRON)
			for (d = 0; d < 3; d++)
			{
				double s = 0;
				for (k = 0; k < 8; k++) s += mesh->x[3 * (size_t)en[eo[c] + k] + d];
				X[3 * (size_t)en[eo[c] + 26] + d] = s / 8;
			}

	/* the border's nodes: on facets with no cell across */
	for (c = 0; c < mesh->nc; c++)
	{
		unsigned int nfa = mesh->f0[c + 1] - mesh->f0[c], f;
		for (f = 0; f < nfa; f++)
		{
			unsigned int n, j;
			unsigned char const* loc;
			if (mesh->fo[mesh->f0[c] + f] != HM_NONE)
				continue;
			loc = vm_facet(mesh->type[c], f, &n);
			for (j = 0; j < n; j++)
			{
				unsigned int s = en[eo[c] + vcu_edge_slot(mesh->type[c], loc[j], loc[(j + 1) % n])];
				if (!isb[s])
				{
					isb[s] = 1;
					bord[nb++] = s;
				}
			}
			if (n == 4)
			{
				unsigned int sv[4], s = HM_NONE;
				for (j = 0; j < 4; j++) sv[j] = loc[j];
				vcu_sort4(sv, 4);
				for (j = 0; j < 6 && s == HM_NONE; j++)
				{
					unsigned int fv[4];
					for (d = 0; d < 4; d++) fv[d] = vcu_hex_face[j][d];
					vcu_sort4(fv, 4);
					if (!memcmp(fv, sv, sizeof fv))
						s = en[eo[c] + 20 + j];
				}
				if (s != HM_NONE && !isb[s])
				{
					isb[s] = 2;
					bord[nb++] = s;
				}
			}
		}
	}
	rep.border_nodes = reference ? nb : 0;

	/* edge nodes, then face nodes (from their Coons points) onto the reference */
	if (reference && nb)
	{
		qaws_scalar* pts = (qaws_scalar*)qaws_internal_scratch_malloc(sizeof(qaws_scalar) * 4 * (size_t)nb);
		qaws_mesh_location* loc = (qaws_mesh_location*)qaws_internal_scratch_malloc(sizeof(qaws_mesh_location) * (size_t)nb);
		unsigned int pass;
		if (!pts || !loc)
		{
			qaws_internal_free(pts);
			qaws_internal_free(loc);
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			goto done;
		}
		for (pass = 1; pass <= 2 && st == QAWS_STATUS_OK; pass++)
		{
			unsigned int m = 0, j;
			for (j = 0; j < nb; j++)
			{
				unsigned int s = bord[j];
				double p[3];
				if (isb[s] != pass)
					continue;
				for (d = 0; d < 3; d++) p[d] = X[3 * (size_t)s + d];
				if (pass == 2)
				{
					/* Coons: the curved edges' middles' mean twice, less the corners' mean */
					unsigned int const* f = fc + 8 * (size_t)s;
					for (d = 0; d < 3; d++)
					{
						double e = 0, v = 0;
						for (k = 0; k < 4; k++)
						{
							e += X[3 * (size_t)f[4 + k] + d] + D[3 * (size_t)f[4 + k] + d];
							v += X[3 * (size_t)f[k] + d];
						}
						p[d] = e / 2 - v / 4;
					}
				}
				for (d = 0; d < 3; d++) pts[3 * (size_t)m + d] = (qaws_scalar)p[d];
				m++;
			}
			if (!m)
				continue;
			st = qaws_mesh_reference_project(loc, NULL, reference, pts, NULL, m);
			if (st == QAWS_STATUS_OK)
				st = qaws_mesh_reference_evaluate(pts, NULL, NULL, reference, loc, m);
			for (j = 0, m = 0; st == QAWS_STATUS_OK && j < nb; j++)
			{
				unsigned int s = bord[j];
				if (isb[s] != pass)
					continue;
				for (d = 0; d < 3; d++) D[3 * (size_t)s + d] = (double)pts[3 * (size_t)m + d] - X[3 * (size_t)s + d];
				m++;
			}
		}
		/* the straight nodes' distances: the sagitta */
		if (st == QAWS_STATUS_OK)
		{
			unsigned int j;
			qaws_scalar* dist = pts + 3 * (size_t)nb;
			for (j = 0; j < nb; j++)
				for (d = 0; d < 3; d++) pts[3 * (size_t)j + d] = (qaws_scalar)X[3 * (size_t)bord[j] + d];
			st = qaws_mesh_reference_project(loc, dist, reference, pts, NULL, nb);
			for (j = 0; st == QAWS_STATUS_OK && j < nb; j++)
				rep.max_sagitta = (qaws_scalar)fmax((double)rep.max_sagitta, (double)dist[j]);
		}
		qaws_internal_free(pts);
		qaws_internal_free(loc);
		if (st != QAWS_STATUS_OK)
			goto done;

		/* the displacement spread to the inside nodes: Jacobi averaging over the elements' other nodes (the vertices still: it fades inward) */
		for (k = 0; k < sweeps; k++)
		{
			memset(acc, 0, sizeof(double) * 3 * (size_t)nn);
			memset(cnt, 0, sizeof(unsigned int) * (size_t)nn);
			for (c = 0; c < mesh->nc; c++)
			{
				unsigned int n0 = mesh->type[c] == QAWS_CELL_TETRAHEDRON ? 4u : 8u, n1 = eo[c + 1] - eo[c], a, b;
				for (a = n0; a < n1; a++)
				{
					unsigned int s = en[eo[c] + a];
					if (isb[s])
						continue;
					for (b = 0; b < n1; b++)
						if (b != a)
						{
							unsigned int t = en[eo[c] + b];
							for (d = 0; d < 3; d++) acc[3 * (size_t)s + d] += D[3 * (size_t)t + d];
							cnt[s]++;
						}
				}
			}
			for (i = mesh->nv; i < nn; i++)
				if (!isb[i] && cnt[i])
					for (d = 0; d < 3; d++) D[3 * (size_t)i + d] = acc[3 * (size_t)i + d] / cnt[i];
		}

		/* the Jacobians: failing elements' displacements halved, then straight */
		for (k = 0; k <= VCU_ROUNDS; k++)
		{
			unsigned int bad = 0;
			for (c = 0; c < mesh->nc; c++)
			{
				unsigned int a;
				if (vcu_ratio(en + eo[c], eo[c + 1] - eo[c], X, D) >= minj)
					continue;
				bad++;
				flag[c] = 1;
				for (a = mesh->type[c] == QAWS_CELL_TETRAHEDRON ? 4u : 8u; a < eo[c + 1] - eo[c]; a++)
					for (d = 0; d < 3; d++) D[3 * (size_t)en[eo[c] + a] + d] *= k < VCU_ROUNDS ? 0.5 : 0;
			}
			if (!bad)
				break;
		}
		for (c = 0; c < mesh->nc; c++) rep.relaxed += flag[c];
		{
			unsigned int j;
			qaws_scalar* pts2 = (qaws_scalar*)qaws_internal_scratch_malloc(sizeof(qaws_scalar) * 4 * (size_t)nb);
			qaws_mesh_location* loc2 = (qaws_mesh_location*)qaws_internal_scratch_malloc(sizeof(qaws_mesh_location) * (size_t)nb);
			if (!pts2 || !loc2)
			{
				qaws_internal_free(pts2);
				qaws_internal_free(loc2);
				st = QAWS_STATUS_ALLOCATION_FAILURE;
				goto done;
			}
			for (j = 0; j < nb; j++)
				for (d = 0; d < 3; d++) pts2[3 * (size_t)j + d] = (qaws_scalar)(X[3 * (size_t)bord[j] + d] + D[3 * (size_t)bord[j] + d]);
			st = qaws_mesh_reference_project(loc2, pts2 + 3 * (size_t)nb, reference, pts2, NULL, nb);
			for (j = 0; st == QAWS_STATUS_OK && j < nb; j++)
				rep.max_distance = (qaws_scalar)fmax((double)rep.max_distance, (double)pts2[3 * (size_t)nb + j]);
			qaws_internal_free(pts2);
			qaws_internal_free(loc2);
			if (st != QAWS_STATUS_OK)
				goto done;
		}
	}
	rep.min_jacobian = (qaws_scalar)HUGE_VAL;
	for (c = 0; c < mesh->nc; c++)
		rep.min_jacobian = (qaws_scalar)fmin((double)rep.min_jacobian, vcu_ratio(en + eo[c], eo[c + 1] - eo[c], X, D));
	for (i = 0; i < 3 * nn; i++) out_positions[i] = (qaws_scalar)(X[i] + D[i]);
done:
	if (report)
		*report = rep;
	qaws_internal_free(key);
	qaws_internal_free(eo);
	qaws_internal_free(ea);
	qaws_internal_free(fc);
	qaws_internal_free(isb);
	qaws_internal_free(flag);
	qaws_internal_free(X);
	qaws_internal_free(D);
	qaws_internal_free(acc);
	qaws_internal_free(cnt);
	qaws_internal_free(bord);
	return st;
}
