/*
 * Octahedral frame fields on tetrahedral meshes (qaws_mesh_field.h): per
 * cell a rotation (its columns the frame's axes), the border's normals held
 * in the cells on it, Gauss-Seidel over the facets' adjacency with each
 * neighbour matched by the signed permutation of its axes nearest the
 * cell's, the average projected back onto the rotations (polar
 * decomposition by Newton's iteration); singular edges where the matchings
 * round an edge compose to a turn.
 *
 * shortcut: frames smoothed by local averaging (Huang et al. 2011 and Ray et
 * al. 2016 optimize the spherical-harmonic representation globally, with
 * fewer and better placed singularities), the integer-grid map that turns
 * a field into hexahedra (FRAME_FIELD) to come.
 */

#include "qaws_mesh_internal.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* matrices are 3 x 3, column-major: m[3 j + i] is row i of column j (column j the axis j) */

static double ff_dot(double const* a, double const* b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void ff_cross(double const* a, double const* b, double* o)
{
	o[0] = a[1] * b[2] - a[2] * b[1];
	o[1] = a[2] * b[0] - a[0] * b[2];
	o[2] = a[0] * b[1] - a[1] * b[0];
}

static double ff_unit(double* v)
{
	double l = sqrt(ff_dot(v, v));
	if (l > 0)
	{
		v[0] /= l; v[1] /= l; v[2] /= l;
	}
	return l;
}

/* the signed permutation P (P[i] = j + 1 or -(j + 1)) taking b's axes onto a's: axis i of a ~ sign * axis j of b, a
   rotation (determinant 1) nearest b^T a */
static void ff_match(double const* a, double const* b, int* P)
{
	static int const perm[6][3] = { { 0, 1, 2 }, { 0, 2, 1 }, { 1, 0, 2 }, { 1, 2, 0 }, { 2, 0, 1 }, { 2, 1, 0 } };
	double D[3][3], best = -1;
	int p, i, bp = 0;
	for (i = 0; i < 3; i++)
	{
		int j;
		for (j = 0; j < 3; j++)
			D[i][j] = ff_dot(a + 3 * i, b + 3 * j);
	}
	for (p = 0; p < 6; p++)
	{
		double s = fabs(D[0][perm[p][0]]) + fabs(D[1][perm[p][1]]) + fabs(D[2][perm[p][2]]);
		if (s > best)
		{
			best = s;
			bp = p;
		}
	}
	for (i = 0; i < 3; i++)
		P[i] = D[i][perm[bp][i]] < 0 ? -(perm[bp][i] + 1) : perm[bp][i] + 1;
	/* a proper rotation: an odd permutation with even signs (or the reverse) flips the weakest axis */
	{
		int odd = bp == 1 || bp == 2 || bp == 5, neg = (P[0] < 0) + (P[1] < 0) + (P[2] < 0), w = 0;
		if ((odd + neg) & 1)
		{
			for (i = 1; i < 3; i++)
				if (fabs(D[i][perm[bp][i]]) < fabs(D[w][perm[bp][w]]))
					w = i;
			P[w] = -P[w];
		}
	}
}

/* b's axes matched to a's, added to m */
static void ff_add_matched(double const* a, double const* b, double* m)
{
	int P[3], i, k;
	ff_match(a, b, P);
	for (i = 0; i < 3; i++)
	{
		int j = (P[i] < 0 ? -P[i] : P[i]) - 1;
		double s = P[i] < 0 ? -1 : 1;
		for (k = 0; k < 3; k++)
			m[3 * i + k] += s * b[3 * j + k];
	}
}

/* the rotation nearest m (polar decomposition, Newton: R <- (R + R^-T) / 2) */
static void ff_polar(double* m)
{
	int it;
	for (it = 0; it < 12; it++)
	{
		double c[9], det, inv[9];
		int k;
		/* the cofactors: R^-T = cof(R) / det(R) */
		ff_cross(m + 3, m + 6, c);
		ff_cross(m + 6, m, c + 3);
		ff_cross(m, m + 3, c + 6);
		det = ff_dot(m, c);
		if (!(fabs(det) > 1e-300))
			break;
		for (k = 0; k < 9; k++)
		{
			inv[k] = c[k] / det;
			m[k] = 0.5 * (m[k] + inv[k]);
		}
	}
	ff_unit(m);
	ff_unit(m + 3);
	ff_cross(m, m + 3, m + 6);
	ff_unit(m + 6);
	ff_cross(m + 6, m, m + 3);
}

/* a frame with axis 0 along n (and axis 1 along n2 made orthogonal, when given) */
static void ff_from_normals(double const* n, double const* n2, double* m)
{
	double t[3];
	memcpy(m, n, sizeof(double) * 3);
	ff_unit(m);
	if (n2)
		memcpy(t, n2, sizeof t);
	else
	{
		t[0] = fabs(m[0]) < 0.6 ? 1 : 0;
		t[1] = fabs(m[0]) < 0.6 ? 0 : 1;
		t[2] = 0;
	}
	{
		double d = ff_dot(t, m);
		t[0] -= d * m[0]; t[1] -= d * m[1]; t[2] -= d * m[2];
	}
	ff_unit(t);
	memcpy(m + 3, t, sizeof t);
	ff_cross(m, m + 3, m + 6);
}

/* the averaged frame bent back to the cell's constraints: the axis nearest each normal set onto it */
static void ff_constrain(double* m, double const* n, unsigned int nn)
{
	double keep[9];
	int i0 = 0, i1 = -1, i;
	if (!nn)
		return;
	memcpy(keep, m, sizeof keep);
	for (i = 1; i < 3; i++)
		if (fabs(ff_dot(m + 3 * i, n)) > fabs(ff_dot(m + 3 * i0, n)))
			i0 = i;
	if (nn > 1)
		for (i = 0; i < 3; i++)
			if (i != i0 && (i1 < 0 || fabs(ff_dot(m + 3 * i, n + 3)) > fabs(ff_dot(m + 3 * i1, n + 3))))
				i1 = i;
	{
		double a0[3], a1[3], a2[3], s = ff_dot(m + 3 * i0, n) < 0 ? -1 : 1;
		int k, i2;
		for (k = 0; k < 3; k++)
			a0[k] = s * n[k];
		if (nn > 1)
		{
			double s1 = ff_dot(m + 3 * i1, n + 3) < 0 ? -1 : 1, d;
			for (k = 0; k < 3; k++)
				a1[k] = s1 * n[3 + k];
			d = ff_dot(a1, a0);
			for (k = 0; k < 3; k++)
				a1[k] -= d * a0[k];
		}
		else
		{
			i1 = (i0 + 1) % 3;
			{
				double d = ff_dot(keep + 3 * i1, a0);
				for (k = 0; k < 3; k++)
					a1[k] = keep[3 * i1 + k] - d * a0[k];
			}
		}
		if (ff_unit(a1) == 0)
		{
			ff_from_normals(a0, NULL, m);
			return;
		}
		i2 = 3 - i0 - i1;
		ff_cross(a0, a1, a2);
		/* the third axis keeps the handedness: (i0, i1, i2) an even order of (0, 1, 2) or not */
		if (!((i1 - i0 + 3) % 3 == 1))
		{
			a2[0] = -a2[0]; a2[1] = -a2[1]; a2[2] = -a2[2];
		}
		memcpy(m + 3 * i0, a0, sizeof a0);
		memcpy(m + 3 * i1, a1, sizeof a1);
		memcpy(m + 3 * i2, a2, sizeof a2);
	}
}

qaws_status qaws_volume_mesh_frame_field(struct qaws_volume_mesh const* mesh, qaws_frame_field_desc const* desc, qaws_scalar* out_frames,
	unsigned char* out_singular)
{
	qaws_volume_mesh const* m = mesh;
	double *R = NULL, *nrm = NULL;
	unsigned char* ncon = NULL;
	unsigned int c, s, i, k, sweeps = desc && desc->sweeps ? desc->sweeps : 60;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (!m || (m->nc && !out_frames))
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (c = 0; c < m->nc; c++)
		if (m->type[c] != QAWS_CELL_TETRAHEDRON)
			return QAWS_STATUS_UNSUPPORTED_OPERATION;
	R = (double*)malloc(sizeof(double) * (9 * (size_t)m->nc + 1));
	nrm = (double*)calloc(6 * (size_t)m->nc + 1, sizeof(double));
	ncon = (unsigned char*)calloc((size_t)m->nc + 1, 1);
	if (!R || !nrm || !ncon)
		goto done;
	/* the border's normals in its cells: one, or two meeting near square */
	for (c = 0; c < m->nc; c++)
	{
		unsigned int f;
		for (f = 0; f < 4; f++)
		{
			unsigned int n;
			unsigned char const* loc;
			double u[3], w[3], q[3];
			double const *a, *b, *d;
			if (m->fo[m->f0[c] + f] != HM_NONE)
				continue;
			loc = vm_facet(QAWS_CELL_TETRAHEDRON, f, &n);
			a = m->x + 3 * (size_t)m->cv[m->c0[c] + loc[0]];
			b = m->x + 3 * (size_t)m->cv[m->c0[c] + loc[1]];
			d = m->x + 3 * (size_t)m->cv[m->c0[c] + loc[2]];
			for (k = 0; k < 3; k++)
			{
				u[k] = b[k] - a[k];
				w[k] = d[k] - a[k];
			}
			ff_cross(u, w, q);
			if (ff_unit(q) == 0)
				continue;
			if (ncon[c] == 0)
			{
				memcpy(nrm + 6 * (size_t)c, q, sizeof q);
				ncon[c] = 1;
			}
			else if (ncon[c] == 1 && fabs(ff_dot(q, nrm + 6 * (size_t)c)) < 0.5)
			{
				memcpy(nrm + 6 * (size_t)c + 3, q, sizeof q);
				ncon[c] = 2;
			}
		}
	}
	/* the start: the border's cells from their normals (the tangent axes from the world's), the rest the world's axes */
	for (c = 0; c < m->nc; c++)
		if (ncon[c])
			ff_from_normals(nrm + 6 * (size_t)c, ncon[c] > 1 ? nrm + 6 * (size_t)c + 3 : NULL, R + 9 * (size_t)c);
		else
		{
			memset(R + 9 * (size_t)c, 0, sizeof(double) * 9);
			R[9 * (size_t)c] = R[9 * (size_t)c + 4] = R[9 * (size_t)c + 8] = 1;
		}
	/* smoothed: each cell its neighbours' frames matched to its own, averaged, made a rotation, bent to its normals */
	for (s = 0; s < sweeps; s++)
		for (i = 0; i < m->nc; i++)
		{
			unsigned int cc = s & 1 ? m->nc - 1 - i : i, f, n = 0;
			double acc[9];
			memset(acc, 0, sizeof acc);
			for (f = 0; f < 4; f++)
			{
				unsigned int o = m->fo[m->f0[cc] + f];
				if (o == HM_NONE)
					continue;
				ff_add_matched(R + 9 * (size_t)cc, R + 9 * (size_t)(o / 4), acc);   /* (tetrahedra: facet o of cell o / 4) */
				n++;
			}
			if (!n)
				continue;
			ff_polar(acc);
			ff_constrain(acc, nrm + 6 * (size_t)cc, ncon[cc]);
			memcpy(R + 9 * (size_t)cc, acc, sizeof acc);
		}
	for (c = 0; c < 9 * m->nc; c++)
		out_frames[c] = (qaws_scalar)R[c];
	/* singular edges: round each inner edge, the matchings composed */
	if (out_singular)
	{
		memset(out_singular, 0, m->nc);
		for (c = 0; c < m->nc; c++)
		{
			unsigned int const* v = m->cv + m->c0[c];
			unsigned int a, b;
			for (a = 0; a < 4; a++)
				for (b = a + 1; b < 4; b++)
				{
					/* the cells round edge {v[a], v[b]}: across the facet opposite one of the other two vertices, and on */
					unsigned int ea = v[a], eb = v[b], cur = c, prev = HM_NONE, guard, lowest = c, ring[64], nr = 0;
					int T[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 }, inner = 1;
					for (guard = 0; guard < 64; guard++)
					{
						unsigned int f, nx = HM_NONE;
						ring[nr++] = cur;
						for (f = 0; f < 4; f++)
						{
							unsigned int w = m->cv[m->c0[cur] + f], o;
							if (w == ea || w == eb)
								continue;
							o = m->fo[m->f0[cur] + f];
							if (o == HM_NONE)
							{
								inner = 0;
								break;
							}
							if (o / 4 != prev && nx == HM_NONE)
								nx = o / 4;
						}
						if (!inner || nx == HM_NONE)
						{
							inner = 0;
							break;
						}
						/* composed: R_cur ~ R_nx Q, so round the cycle R_c ~ R_c (Q_last ... Q_first); T = Q T */
						{
							int P[3], U[9], r, q, z;
							ff_match(R + 9 * (size_t)cur, R + 9 * (size_t)nx, P);
							for (r = 0; r < 3; r++)
								for (q = 0; q < 3; q++)
								{
									int sum = 0;
									/* Q's column z holds its sign at row |P[z]| - 1 */
									for (z = 0; z < 3; z++)
										if ((P[z] < 0 ? -P[z] : P[z]) - 1 == r)
											sum += (P[z] < 0 ? -1 : 1) * T[3 * q + z];
									U[3 * q + r] = sum;
								}
							memcpy(T, U, sizeof T);
						}
						lowest = nx < lowest ? nx : lowest;
						prev = cur;
						cur = nx;
						if (cur == c)
							break;
					}
					if (!inner || cur != c || lowest != c)
						continue;   /* (on the border, or seen from its lowest cell) */
					if (!(T[0] == 1 && T[4] == 1 && T[8] == 1))
						for (k = 0; k < nr; k++)
						{
							/* the bit of the edge in each cell round it: edges 0-1 0-2 0-3 1-2 1-3 2-3 */
							unsigned int const* w = m->cv + m->c0[ring[k]];
							unsigned int x0, x1, e = 0;
							for (x0 = 0; x0 < 4; x0++)
								for (x1 = x0 + 1; x1 < 4; x1++, e++)
									if ((w[x0] == ea && w[x1] == eb) || (w[x0] == eb && w[x1] == ea))
										out_singular[ring[k]] |= (unsigned char)(1u << e);
						}
				}
		}
	}
	st = QAWS_STATUS_OK;
done:
	free(R);
	free(nrm);
	free(ncon);
	return st;
}
