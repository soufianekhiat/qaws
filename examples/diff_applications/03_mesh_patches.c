/* ================================================================== */
/*  3. Triangulated mesh to NURBS patches with ADMM                   */
/* ================================================================== */

#define MP_PATCHES 6
#define MP_N 7                       /* control points per direction */
#define MP_CP (MP_N * MP_N)
#define MP_GRID 26                   /* mesh vertices per face edge */
#define MP_OUTER 40                  /* ADMM iterations */
#define MP_INNER 12                  /* Adam steps per x-update */

typedef struct mesh
{
	int nv, nt;
	double* v;      /* nv * 3 */
	int* t;         /* nt * 3 */
} mesh;

/* A bumpy closed shape (radius as a function of direction). */
static double blob_radius(double x, double y, double z)
{
	double phi = atan2(y, x), th = acos(z / sqrt(x * x + y * y + z * z));
	return 1.0 + 0.16 * sin(3 * phi) * sin(th) * sin(th) + 0.12 * cos(2 * th) + 0.08 * sin(5 * phi + 2 * th) * sin(th);
}

/* Cube face f: axis a = f / 2, sign s, tangent axes b, c. */
static void cube_point(int f, double S, double T, double* p)
{
	int a = f / 2, b = (a + 1) % 3, c = (a + 2) % 3;
	double s = (f & 1) ? -1.0 : 1.0;
	p[a] = s;
	p[b] = (f & 1) ? T : S;   /* flip on negative faces so every face is outward-oriented */
	p[c] = (f & 1) ? S : T;
}

static double g_blob_amp = 1.0;   /* 0: exact sphere, 1: bumpy blob */

static void mesh_blob(mesh* m)
{
	int f, i, j, k = 0, t = 0;
	m->nv = MP_PATCHES * MP_GRID * MP_GRID;
	m->nt = MP_PATCHES * (MP_GRID - 1) * (MP_GRID - 1) * 2;
	m->v = (double*)malloc(sizeof(double) * 3 * (size_t)m->nv);
	m->t = (int*)malloc(sizeof(int) * 3 * (size_t)m->nt);
	for (f = 0; f < MP_PATCHES; f++)
	{
		for (i = 0; i < MP_GRID; i++)
			for (j = 0; j < MP_GRID; j++)
			{
				double p[3], len, r;
				/* equal-angle cube map for even sampling */
				cube_point(f, tan((i / (double)(MP_GRID - 1) - 0.5) * PI / 2), tan((j / (double)(MP_GRID - 1) - 0.5) * PI / 2), p);
				len = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
				r = 1.0 + g_blob_amp * (blob_radius(p[0], p[1], p[2]) - 1.0);
				m->v[3 * k] = p[0] / len * r;
				m->v[3 * k + 1] = p[1] / len * r;
				m->v[3 * k + 2] = p[2] / len * r;
				k++;
			}
		for (i = 0; i + 1 < MP_GRID; i++)
			for (j = 0; j + 1 < MP_GRID; j++)
			{
				int a = f * MP_GRID * MP_GRID + i * MP_GRID + j;
				m->t[3 * t] = a; m->t[3 * t + 1] = a + MP_GRID; m->t[3 * t + 2] = a + MP_GRID + 1; t++;
				m->t[3 * t] = a; m->t[3 * t + 1] = a + MP_GRID + 1; m->t[3 * t + 2] = a + 1; t++;
			}
	}
}

/* Optional OBJ input (v / f lines, triangles or polygons fanned). */
static int mesh_load_obj(char const* path, mesh* m)
{
	FILE* f = fopen(path, "r");
	char line[512];
	int cv = 0, ct = 0;
	if (!f)
		return 0;
	m->nv = m->nt = 0;
	while (fgets(line, sizeof(line), f))
	{
		if (line[0] == 'v' && line[1] == ' ') m->nv++;
		else if (line[0] == 'f' && line[1] == ' ')
		{
			int n = 0;
			char* p = line + 1;
			while (*p)
			{
				while (*p == ' ') p++;
				if (*p && *p != '\n' && *p != '\r') { n++; while (*p && *p != ' ') p++; }
				else break;
			}
			m->nt += n - 2;
		}
	}
	m->v = (double*)malloc(sizeof(double) * 3 * (size_t)m->nv);
	m->t = (int*)malloc(sizeof(int) * 3 * (size_t)m->nt);
	rewind(f);
	while (fgets(line, sizeof(line), f))
	{
		if (line[0] == 'v' && line[1] == ' ')
		{
			sscanf(line + 2, "%lf %lf %lf", &m->v[3 * cv], &m->v[3 * cv + 1], &m->v[3 * cv + 2]);
			cv++;
		}
		else if (line[0] == 'f' && line[1] == ' ')
		{
			int idx[64], n = 0, k;
			char* p = line + 1;
			while (*p && n < 64)
			{
				while (*p == ' ') p++;
				if (!*p || *p == '\n' || *p == '\r') break;
				idx[n++] = atoi(p) - 1;
				while (*p && *p != ' ') p++;
			}
			for (k = 1; k + 1 < n; k++)
			{
				m->t[3 * ct] = idx[0]; m->t[3 * ct + 1] = idx[k]; m->t[3 * ct + 2] = idx[k + 1];
				ct++;
			}
		}
	}
	fclose(f);
	m->nt = ct;
	{
		/* center and scale to radius ~1 */
		double c[3] = { 0, 0, 0 }, rmax = 0;
		int i, a;
		for (i = 0; i < m->nv; i++)
			for (a = 0; a < 3; a++)
				c[a] += m->v[3 * i + a] / m->nv;
		for (i = 0; i < m->nv; i++)
		{
			double r = 0;
			for (a = 0; a < 3; a++)
			{
				m->v[3 * i + a] -= c[a];
				r += m->v[3 * i + a] * m->v[3 * i + a];
			}
			if (r > rmax) rmax = r;
		}
		rmax = sqrt(rmax);
		for (i = 0; i < 3 * m->nv; i++)
			m->v[i] /= rmax;
	}
	return 1;
}

typedef struct patch_fit
{
	qaws_scalar cps[MP_CP * 3];
	qaws_scalar dual[MP_CP * 3];     /* scaled ADMM duals (boundary nodes) */
	int node[MP_CP];                 /* global node, -1 when interior */
	int n;                           /* assigned mesh points */
	int* pts;
	qaws_scalar* uv;                 /* foot points */
	adam opt;
} patch_fit;

static qaws_scalar const g_mp_knots[MP_N + 4] = { 0, 0, 0, 0, 0.25f, 0.5f, 0.75f, 1, 1, 1, 1 };
static qaws_scalar g_mp_weights[MP_CP];

static qaws_surface* mp_surface(qaws_scalar const* cps)
{
	/* bicubic B-spline: a NURBS patch with unit weights, linear in its
	   control points (exact x-updates, direct thin-plate HVP) */
	qaws_surface_bspline_desc d;
	qaws_surface* s = NULL;
	memset(&d, 0, sizeof(d));
	d.u_degree = 3;
	d.v_degree = 3;
	d.control_points = (qaws_vec3 const*)cps;
	d.u_point_count = MP_N;
	d.v_point_count = MP_N;
	d.u_knots = g_mp_knots;
	d.u_knot_count = MP_N + 4;
	d.v_knots = g_mp_knots;
	d.v_knot_count = MP_N + 4;
	qaws_surface_create_bspline(&d, &s);
	return s;
}

/* Newton steps on the foot point of every assigned point (clamped to the
   patch domain); returns the patch's mean squared distance. */
static double mp_project(patch_fit* p, mesh const* m, qaws_surface const* s, int steps)
{
	double e = 0;
	int i, it;
	for (i = 0; i < p->n; i++)
	{
		double const* x = &m->v[3 * p->pts[i]];
		qaws_scalar* uv = &p->uv[2 * i];
		qaws_surface_jet j;
		double rx, ry, rz;
		for (it = 0; it <= steps; it++)
		{
			double g0, g1, a, b, c, det, du, dv;
			qaws_surface_eval_jet(s, uv[0], uv[1], QAWS_SJET_ORDER2, &j);
			rx = j.d[0].x - x[0]; ry = j.d[0].y - x[1]; rz = j.d[0].z - x[2];
			if (it == steps)
				break;
			g0 = rx * j.d[1].x + ry * j.d[1].y + rz * j.d[1].z;
			g1 = rx * j.d[2].x + ry * j.d[2].y + rz * j.d[2].z;
			a = j.d[1].x * j.d[1].x + j.d[1].y * j.d[1].y + j.d[1].z * j.d[1].z + rx * j.d[3].x + ry * j.d[3].y + rz * j.d[3].z;
			b = j.d[1].x * j.d[2].x + j.d[1].y * j.d[2].y + j.d[1].z * j.d[2].z + rx * j.d[4].x + ry * j.d[4].y + rz * j.d[4].z;
			c = j.d[2].x * j.d[2].x + j.d[2].y * j.d[2].y + j.d[2].z * j.d[2].z + rx * j.d[5].x + ry * j.d[5].y + rz * j.d[5].z;
			det = a * c - b * b;
			if (!(det > 1e-12) || !(a > 0))
				break;
			du = -(c * g0 - b * g1) / det;
			dv = -(a * g1 - b * g0) / det;
			uv[0] = (qaws_scalar)(uv[0] + du < 0 ? 0 : (uv[0] + du > 1 ? 1 : uv[0] + du));
			uv[1] = (qaws_scalar)(uv[1] + dv < 0 ? 0 : (uv[1] + dv > 1 ? 1 : uv[1] + dv));
		}
		e += rx * rx + ry * ry + rz * rz;
	}
	return p->n ? e / p->n : 0;
}
/* Builds the patches, the point assignment and the seam nodes. */
static int mp_setup(mesh const* m, patch_fit* pf, int* node_count)
{
	double keys[MP_PATCHES * MP_CP][3];
	int nkeys = 0, f, i, j, k;
	for (i = 0; i < MP_CP; i++)
		g_mp_weights[i] = 1;
	for (f = 0; f < MP_PATCHES; f++)
	{
		patch_fit* p = &pf[f];
		memset(p, 0, sizeof(*p));
		for (i = 0; i < MP_N; i++)
			for (j = 0; j < MP_N; j++)
			{
				double c[3], len;
				int b = i == 0 || j == 0 || i == MP_N - 1 || j == MP_N - 1, idx = -1;
				cube_point(f, -1 + 2.0 * i / (MP_N - 1), -1 + 2.0 * j / (MP_N - 1), c);
				len = sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
				p->cps[(i * MP_N + j) * 3 + 0] = (qaws_scalar)(c[0] / len);
				p->cps[(i * MP_N + j) * 3 + 1] = (qaws_scalar)(c[1] / len);
				p->cps[(i * MP_N + j) * 3 + 2] = (qaws_scalar)(c[2] / len);
				if (b)
				{
					for (k = 0; k < nkeys; k++)
						if (fabs(keys[k][0] - c[0]) + fabs(keys[k][1] - c[1]) + fabs(keys[k][2] - c[2]) < 1e-9)
							idx = k;
					if (idx < 0)
					{
						idx = nkeys++;
						keys[idx][0] = c[0]; keys[idx][1] = c[1]; keys[idx][2] = c[2];
					}
				}
				p->node[i * MP_N + j] = idx;
			}
	}
	*node_count = nkeys;
	/* assignment by dominant direction, initial foot points from the cube map */
	for (f = 0; f < MP_PATCHES; f++)
	{
		pf[f].pts = (int*)malloc(sizeof(int) * (size_t)m->nv);
		pf[f].uv = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * (size_t)m->nv);
	}
	for (i = 0; i < m->nv; i++)
	{
		double const* x = &m->v[3 * i];
		int a = 0, b, c, ff;
		double s, S, T;
		for (k = 1; k < 3; k++)
			if (fabs(x[k]) > fabs(x[a])) a = k;
		s = x[a] > 0 ? 1 : -1;
		ff = 2 * a + (s < 0);
		b = (a + 1) % 3;
		c = (a + 2) % 3;
		S = x[b] / fabs(x[a]);
		T = x[c] / fabs(x[a]);
		if (ff & 1) { double tmp = S; S = T; T = tmp; }
		pf[ff].pts[pf[ff].n] = i;
		pf[ff].uv[2 * pf[ff].n] = (qaws_scalar)((S + 1) / 2);
		pf[ff].uv[2 * pf[ff].n + 1] = (qaws_scalar)((T + 1) / 2);
		pf[ff].n++;
	}
	return 1;
}

/* Largest distance between copies of a seam node. */
static double mp_seam_gap(patch_fit const* pf, int nodes)
{
	double gap = 0;
	int g, f, k, f2, k2;
	for (g = 0; g < nodes; g++)
		for (f = 0; f < MP_PATCHES; f++)
			for (k = 0; k < MP_CP; k++)
				if (pf[f].node[k] == g)
					for (f2 = f + 1; f2 < MP_PATCHES; f2++)
						for (k2 = 0; k2 < MP_CP; k2++)
							if (pf[f2].node[k2] == g)
							{
								double dx = pf[f].cps[3 * k] - pf[f2].cps[3 * k2];
								double dy = pf[f].cps[3 * k + 1] - pf[f2].cps[3 * k2 + 1];
								double dz = pf[f].cps[3 * k + 2] - pf[f2].cps[3 * k2 + 2];
								double d = sqrt(dx * dx + dy * dy + dz * dz);
								if (d > gap) gap = d;
							}
	return gap;
}

/* Thin-plate Hessian of a patch (one coordinate; the energy decouples
   over x, y, z), assembled column by column from direct HVPs. */
static void mp_thin_plate_hessian(double* H)
{
	qaws_scalar cps[MP_CP * 3], dir[MP_CP * 3], out[MP_CP * 3];
	qaws_surface* s;
	int k, j;
	memset(cps, 0, sizeof(cps));
	s = mp_surface(cps);
	for (k = 0; k < MP_CP; k++)
	{
		qaws_field_view fd, fo;
		qaws_diff_views vd, vo;
		memset(dir, 0, sizeof(dir));
		memset(out, 0, sizeof(out));
		dir[3 * k] = 1;
		vd = one_field(&fd, QAWS_FIELD_CONTROL_POINTS, dir, MP_CP, 3);
		vo = one_field(&fo, QAWS_FIELD_CONTROL_POINTS, out, MP_CP, 3);
		qaws_surface_functional_hvp(NULL, s, QAWS_FUNCTIONAL_THIN_PLATE, 4, &vd, &vo);
		for (j = 0; j < MP_CP; j++)
			H[j * MP_CP + k] = out[3 * j];
	}
	qaws_surface_destroy(s);
}

/* In-place Cholesky solve of an SPD n x n system for three right-hand sides. */
static void mp_cholesky_solve(double* A, double b[3][MP_CP], int n)
{
	int i, j, k, c;
	for (i = 0; i < n; i++)
		for (j = 0; j <= i; j++)
		{
			double v = A[i * n + j];
			for (k = 0; k < j; k++)
				v -= A[i * n + k] * A[j * n + k];
			A[i * n + j] = (i == j) ? sqrt(v > 1e-18 ? v : 1e-18) : v / A[j * n + j];
		}
	for (c = 0; c < 3; c++)
	{
		for (i = 0; i < n; i++)
		{
			double v = b[c][i];
			for (k = 0; k < i; k++)
				v -= A[i * n + k] * b[c][k];
			b[c][i] = v / A[i * n + i];
		}
		for (i = n; i-- > 0;)
		{
			double v = b[c][i];
			for (k = i + 1; k < n; k++)
				v -= A[k * n + i] * b[c][k];
			b[c][i] = v / A[i * n + i];
		}
	}
}

/* Exact x-update at fixed foot points:
   (2 A + tp H + rho D) P = 2 b + rho D (z - u), with A, b from the exact
   basis weights of every foot point (qaws_surface_local_support). */
static void mp_xupdate(patch_fit* p, mesh const* m, double const* H, double tp, double rho, qaws_scalar const* z)
{
	static double A[MP_CP * MP_CP];
	double b[3][MP_CP];
	qaws_surface* s = mp_surface(p->cps);
	int i, a, c, k;
	memset(A, 0, sizeof(A));
	memset(b, 0, sizeof(b));
	for (i = 0; i < p->n; i++)
	{
		qaws_surface_support sup;
		double const* x = &m->v[3 * p->pts[i]];
		int idx[16], n = 0;
		double w[16];
		qaws_surface_local_support(s, p->uv[2 * i], p->uv[2 * i + 1], 0, &sup);
		for (a = 0; a < (int)sup.u_count; a++)
			for (c = 0; c < (int)sup.v_count; c++)
			{
				idx[n] = (int)((sup.u_first + a) * sup.u_stride + (sup.v_first + c) * sup.v_stride);
				w[n++] = sup.u_weights[0][a] * sup.v_weights[0][c];
			}
		for (a = 0; a < n; a++)
		{
			for (c = 0; c < n; c++)
				A[idx[a] * MP_CP + idx[c]] += 2 * w[a] * w[c] / p->n;
			for (k = 0; k < 3; k++)
				b[k][idx[a]] += 2 * w[a] * x[k] / p->n;
		}
	}
	for (k = 0; k < MP_CP * MP_CP; k++)
		A[k] += tp * H[k];
	if (rho > 0)
		for (k = 0; k < MP_CP; k++)
			if (p->node[k] >= 0)
			{
				A[k * MP_CP + k] += rho;
				for (c = 0; c < 3; c++)
					b[c][k] += rho * (z[3 * p->node[k] + c] - p->dual[3 * k + c]);
			}
	mp_cholesky_solve(A, b, MP_CP);
	for (k = 0; k < MP_CP; k++)
		for (c = 0; c < 3; c++)
			p->cps[3 * k + c] = (qaws_scalar)b[c][k];
	qaws_surface_destroy(s);
}

/* ADMM on the seam copies: foot points, x-updates (independent per patch),
   z = mean of (x + u) over copies, u += x - z. rho = 0 gives the
   independent fits. */
static void mp_solve(mesh const* m, patch_fit* pf, int nodes, double rho, double* rms, double* gap)
{
	static double H[MP_CP * MP_CP];
	double tp = 2e-6;
	qaws_scalar* z = (qaws_scalar*)calloc((size_t)nodes * 3, sizeof(qaws_scalar));
	int* cnt = (int*)calloc((size_t)nodes, sizeof(int));
	int it, f, k, c;

	mp_thin_plate_hessian(H);
	for (f = 0; f < MP_PATCHES; f++)
		for (k = 0; k < MP_CP; k++)
			if (pf[f].node[k] >= 0)
				for (c = 0; c < 3; c++)
					z[3 * pf[f].node[k] + c] += pf[f].cps[3 * k + c];
	for (f = 0; f < MP_PATCHES; f++)
		for (k = 0; k < MP_CP; k++)
			if (pf[f].node[k] >= 0)
				cnt[pf[f].node[k]]++;
	for (k = 0; k < nodes; k++)
		for (c = 0; c < 3; c++)
			z[3 * k + c] /= cnt[k];
	for (it = 0; it < MP_OUTER; it++)
	{
		double e = 0;
		int npts = 0;
		for (f = 0; f < MP_PATCHES; f++)
		{
			qaws_surface* s = mp_surface(pf[f].cps);
			e += mp_project(&pf[f], m, s, 3) * pf[f].n;
			npts += pf[f].n;
			qaws_surface_destroy(s);
		}
		rms[it] = sqrt(e / npts);
		for (f = 0; f < MP_PATCHES; f++)
			mp_xupdate(&pf[f], m, H, tp, rho, z);
		memset(z, 0, sizeof(qaws_scalar) * (size_t)nodes * 3);
		for (f = 0; f < MP_PATCHES; f++)
			for (k = 0; k < MP_CP; k++)
				if (pf[f].node[k] >= 0)
					for (c = 0; c < 3; c++)
						z[3 * pf[f].node[k] + c] += (pf[f].cps[3 * k + c] + pf[f].dual[3 * k + c]) / cnt[pf[f].node[k]];
		if (rho > 0)
			for (f = 0; f < MP_PATCHES; f++)
				for (k = 0; k < MP_CP; k++)
					if (pf[f].node[k] >= 0)
						for (c = 0; c < 3; c++)
							pf[f].dual[3 * k + c] += pf[f].cps[3 * k + c] - z[3 * pf[f].node[k] + c];
		gap[it] = mp_seam_gap(pf, nodes);
	}
	free(z);
	free(cnt);
}

/* --- rational refinement: weights become free ----------------------- */

static qaws_surface* mp_surface_nurbs(qaws_scalar const* cps, qaws_scalar const* w)
{
	qaws_surface_nurbs_desc d;
	qaws_surface* s = NULL;
	memset(&d, 0, sizeof(d));
	d.u_degree = 3;
	d.v_degree = 3;
	d.control_points = (qaws_vec3 const*)cps;
	d.u_point_count = MP_N;
	d.v_point_count = MP_N;
	d.weights = w;
	d.u_knots = g_mp_knots;
	d.u_knot_count = MP_N + 4;
	d.v_knots = g_mp_knots;
	d.v_knot_count = MP_N + 4;
	qaws_surface_create_nurbs(&d, &s);
	return s;
}

/* Mean squared distance at the foot points of a NURBS patch with its
   gradient on control points and weights (exact rational adjoints). */
static double mp_rational_energy(patch_fit const* p, mesh const* m, qaws_scalar const* cps, qaws_scalar const* w,
	qaws_scalar* g_cp, qaws_scalar* g_w)
{
	qaws_surface* s = mp_surface_nurbs(cps, w);
	qaws_surface_jet* bars = (qaws_surface_jet*)malloc(sizeof(qaws_surface_jet) * (size_t)(p->n + 1));
	qaws_scalar* us = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * (size_t)(p->n + 1));
	qaws_scalar* vs = us + p->n;
	qaws_field_view fv[2];
	qaws_diff_views views;
	double e = 0;
	int i;
	for (i = 0; i < p->n; i++)
	{
		double const* x = &m->v[3 * p->pts[i]];
		qaws_surface_jet j;
		double rx, ry, rz;
		us[i] = p->uv[2 * i];
		vs[i] = p->uv[2 * i + 1];
		qaws_surface_eval_jet(s, us[i], vs[i], QAWS_SJET_P, &j);
		rx = j.d[0].x - x[0]; ry = j.d[0].y - x[1]; rz = j.d[0].z - x[2];
		e += (rx * rx + ry * ry + rz * rz) / p->n;
		memset(&bars[i], 0, sizeof(bars[i]));
		bars[i].d[0] = v3((qaws_scalar)(2 * rx / p->n), (qaws_scalar)(2 * ry / p->n), (qaws_scalar)(2 * rz / p->n));
		bars[i].channels = QAWS_SJET_P;
	}
	memset(g_cp, 0, sizeof(qaws_scalar) * MP_CP * 3);
	memset(g_w, 0, sizeof(qaws_scalar) * MP_CP);
	fv[0] = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, g_cp, MP_CP, 3);
	fv[1] = qaws_field_view_make(QAWS_FIELD_WEIGHTS, g_w, MP_CP, 1);
	views.fields = fv;
	views.field_count = 2;
	views.children = NULL;
	views.child_count = 0;
	if (p->n)
		qaws_surface_eval_batch_adjoint(NULL, s, us, vs, (unsigned int)p->n, QAWS_SJET_P, bars, &views, NULL, NULL);
	free(bars);
	free(us);
	qaws_surface_destroy(s);
	return e;
}

/* Joint Adam on control points and log-weights of all patches; after each
   step the copies of every seam node (positions and weights) are replaced
   by their mean, so the patch set stays closed. Foot points are
   re-projected every 10 steps. Returns the final RMS. */
static double mp_refine_rational(mesh const* m, patch_fit* pf, int nodes, qaws_scalar (*weights)[MP_CP], int iters, int free_weights,
	double* rms_hist)
{
	static adam opt[MP_PATCHES];
	qaws_scalar* zsum = (qaws_scalar*)calloc((size_t)nodes * 4, sizeof(qaws_scalar));
	int* cnt = (int*)calloc((size_t)nodes, sizeof(int));
	int it, f, k, c;
	double rms = 0;
	for (f = 0; f < MP_PATCHES; f++)
	{
		memset(&opt[f], 0, sizeof(adam));
		for (k = 0; k < MP_CP; k++)
		{
			weights[f][k] = 1;
			if (pf[f].node[k] >= 0)
				cnt[pf[f].node[k]]++;
		}
	}
	for (it = 0; it <= iters; it++)
	{
		double e = 0;
		int npts = 0;
		if (it % 10 == 0 || it == iters)
			for (f = 0; f < MP_PATCHES; f++)
			{
				qaws_surface* s = mp_surface_nurbs(pf[f].cps, weights[f]);
				double ef = mp_project(&pf[f], m, s, 3);
				qaws_surface_destroy(s);
				e += ef * pf[f].n;
				npts += pf[f].n;
			}
		if (it % 10 == 0 || it == iters)
		{
			rms = sqrt(e / npts);
			if (rms_hist)
				rms_hist[it / 10] = rms;
		}
		if (it == iters)
			break;
		for (f = 0; f < MP_PATCHES; f++)
		{
			qaws_scalar g_cp[MP_CP * 3], g_w[MP_CP], x[MP_CP * 4], g[MP_CP * 4];
			mp_rational_energy(&pf[f], m, pf[f].cps, weights[f], g_cp, g_w);
			/* parameters: positions and log-weights (d/dlog w = w d/dw) */
			for (k = 0; k < MP_CP * 3; k++) { x[k] = pf[f].cps[k]; g[k] = g_cp[k]; }
			for (k = 0; k < MP_CP; k++)
			{
				x[MP_CP * 3 + k] = (qaws_scalar)log(weights[f][k]);
				g[MP_CP * 3 + k] = free_weights ? weights[f][k] * g_w[k] : 0;
			}
			adam_step(&opt[f], x, g, MP_CP * 4, 0.002 * (1.0 - 0.7 * it / (double)iters));
			for (k = 0; k < MP_CP * 3; k++) pf[f].cps[k] = x[k];
			for (k = 0; k < MP_CP; k++) weights[f][k] = (qaws_scalar)exp(x[MP_CP * 3 + k]);
		}
		/* seam consensus by projection */
		memset(zsum, 0, sizeof(qaws_scalar) * (size_t)nodes * 4);
		for (f = 0; f < MP_PATCHES; f++)
			for (k = 0; k < MP_CP; k++)
				if (pf[f].node[k] >= 0)
				{
					for (c = 0; c < 3; c++)
						zsum[4 * pf[f].node[k] + c] += pf[f].cps[3 * k + c] / cnt[pf[f].node[k]];
					zsum[4 * pf[f].node[k] + 3] += weights[f][k] / cnt[pf[f].node[k]];
				}
		for (f = 0; f < MP_PATCHES; f++)
			for (k = 0; k < MP_CP; k++)
				if (pf[f].node[k] >= 0)
				{
					for (c = 0; c < 3; c++)
						pf[f].cps[3 * k + c] = zsum[4 * pf[f].node[k] + c];
					weights[f][k] = zsum[4 * pf[f].node[k] + 3];
				}
	}
	free(zsum);
	free(cnt);
	return rms;
}

/* --- rendering ---------------------------------------------------- */

typedef struct view3
{
	double cx, cy, scale, yaw, pitch;
} view3;

static void view_xform(view3 const* v, double const* p, double* sx, double* sy, double* depth)
{
	double cy = cos(v->yaw), sy_ = sin(v->yaw), cp = cos(v->pitch), sp = sin(v->pitch);
	double x = cy * p[0] - sy_ * p[1];
	double y = sy_ * p[0] + cy * p[1];
	double z = p[2];
	double y2 = cp * y - sp * z, z2 = sp * y + cp * z;
	*sx = v->cx + x * v->scale;
	*sy = v->cy - z2 * v->scale;
	*depth = y2;   /* larger = farther */
}

typedef struct poly3
{
	double xy[8];
	int n;
	double depth;
	char color[40];
} poly3;

static int poly_cmp(void const* a, void const* b)
{
	double da = ((poly3 const*)a)->depth, db = ((poly3 const*)b)->depth;
	return da < db ? 1 : (da > db ? -1 : 0);
}

static void shade(double const* n, double const* base, char* out)
{
	double l[3] = { -0.45, -0.55, 0.70 }, d = n[0] * l[0] + n[1] * l[1] + n[2] * l[2], k;
	double ln = sqrt(l[0] * l[0] + l[1] * l[1] + l[2] * l[2]);
	d /= ln;
	k = 0.30 + 0.70 * (d > 0 ? d : 0);
	sprintf(out, "rgb(%d,%d,%d)", (int)(255 * base[0] * k), (int)(255 * base[1] * k), (int)(255 * base[2] * k));
}

static void draw_polys(svg* s, poly3* polys, int n, double stroke)
{
	int i;
	qsort(polys, (size_t)n, sizeof(poly3), poly_cmp);
	for (i = 0; i < n; i++)
	{
		poly3 const* p = &polys[i];
		if (p->n == 3)
			fprintf(s->f, "<polygon points=\"%.1f,%.1f %.1f,%.1f %.1f,%.1f\" fill=\"%s\" stroke=\"%s\" stroke-width=\"%.2f\"/>\n",
				p->xy[0], p->xy[1], p->xy[2], p->xy[3], p->xy[4], p->xy[5], p->color, p->color, stroke);
		else
			fprintf(s->f, "<polygon points=\"%.1f,%.1f %.1f,%.1f %.1f,%.1f %.1f,%.1f\" fill=\"%s\" stroke=\"%s\" stroke-width=\"%.2f\"/>\n",
				p->xy[0], p->xy[1], p->xy[2], p->xy[3], p->xy[4], p->xy[5], p->xy[6], p->xy[7], p->color, p->color, stroke);
	}
}

static void render_mesh(svg* s, view3 const* v, mesh const* m)
{
	poly3* polys = (poly3*)malloc(sizeof(poly3) * (size_t)m->nt);
	int t, n = 0;
	double base[3] = { 0.78, 0.80, 0.84 };
	for (t = 0; t < m->nt; t++)
	{
		double const* a = &m->v[3 * m->t[3 * t]];
		double const* b = &m->v[3 * m->t[3 * t + 1]];
		double const* c = &m->v[3 * m->t[3 * t + 2]];
		double e1[3], e2[3], nn[3], len, d0, d1, d2, cen[3];
		int k;
		for (k = 0; k < 3; k++) { e1[k] = b[k] - a[k]; e2[k] = c[k] - a[k]; cen[k] = (a[k] + b[k] + c[k]) / 3; }
		nn[0] = e1[1] * e2[2] - e1[2] * e2[1];
		nn[1] = e1[2] * e2[0] - e1[0] * e2[2];
		nn[2] = e1[0] * e2[1] - e1[1] * e2[0];
		len = sqrt(nn[0] * nn[0] + nn[1] * nn[1] + nn[2] * nn[2]) + 1e-12;
		if (nn[0] * cen[0] + nn[1] * cen[1] + nn[2] * cen[2] < 0) len = -len;
		for (k = 0; k < 3; k++) nn[k] /= len;
		view_xform(v, a, &polys[n].xy[0], &polys[n].xy[1], &d0);
		view_xform(v, b, &polys[n].xy[2], &polys[n].xy[3], &d1);
		view_xform(v, c, &polys[n].xy[4], &polys[n].xy[5], &d2);
		{
			/* back faces away from the viewer are skipped */
			double vd[3] = { 0, 0, 0 }, cy = cos(v->yaw), sy_ = sin(v->yaw), cp = cos(v->pitch), sp = sin(v->pitch);
			double y = sy_ * nn[0] + cy * nn[1];
			if (cp * y - sp * nn[2] > 0.02)
				continue;
			(void)vd;
		}
		polys[n].n = 3;
		polys[n].depth = (d0 + d1 + d2) / 3;
		shade(nn, base, polys[n].color);
		n++;
	}
	draw_polys(s, polys, n, 0.4);
	free(polys);
}

static void render_patches(svg* s, view3 const* v, patch_fit const* pf, int grid, int heat_error, mesh const* m)
{
	static double const colors[MP_PATCHES][3] = {
		{ 0.35, 0.55, 0.95 }, { 0.95, 0.55, 0.30 }, { 0.40, 0.80, 0.45 },
		{ 0.85, 0.40, 0.70 }, { 0.95, 0.80, 0.30 }, { 0.40, 0.80, 0.85 } };
	poly3* polys = (poly3*)malloc(sizeof(poly3) * MP_PATCHES * (size_t)grid * grid);
	int f, i, j, n = 0;
	(void)heat_error;
	(void)m;
	for (f = 0; f < MP_PATCHES; f++)
	{
		qaws_surface* srf = mp_surface(pf[f].cps);
		for (i = 0; i < grid; i++)
			for (j = 0; j < grid; j++)
			{
				qaws_surface_eval_result r[4], mid;
				double p[3], d, dsum = 0, nn[3], cy = cos(v->yaw), sy_ = sin(v->yaw), cp = cos(v->pitch), sp = sin(v->pitch);
				int c;
				qaws_scalar uu[4] = { (qaws_scalar)(i / (double)grid), (qaws_scalar)((i + 1) / (double)grid),
					(qaws_scalar)((i + 1) / (double)grid), (qaws_scalar)(i / (double)grid) };
				qaws_scalar vv[4] = { (qaws_scalar)(j / (double)grid), (qaws_scalar)(j / (double)grid),
					(qaws_scalar)((j + 1) / (double)grid), (qaws_scalar)((j + 1) / (double)grid) };
				for (c = 0; c < 4; c++)
				{
					qaws_surface_evaluate(srf, uu[c], vv[c], QAWS_SURFACE_EVAL_POSITION, &r[c]);
					p[0] = r[c].position.x; p[1] = r[c].position.y; p[2] = r[c].position.z;
					view_xform(v, p, &polys[n].xy[2 * c], &polys[n].xy[2 * c + 1], &d);
					dsum += d;
				}
				qaws_surface_evaluate(srf, (qaws_scalar)((i + 0.5) / grid), (qaws_scalar)((j + 0.5) / grid),
					QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL, &mid);
				nn[0] = mid.normal.x; nn[1] = mid.normal.y; nn[2] = mid.normal.z;
				if (nn[0] * mid.position.x + nn[1] * mid.position.y + nn[2] * mid.position.z < 0)
				{
					nn[0] = -nn[0]; nn[1] = -nn[1]; nn[2] = -nn[2];
				}
				if (cp * (sy_ * nn[0] + cy * nn[1]) - sp * nn[2] > 0.02)
					continue;
				polys[n].n = 4;
				polys[n].depth = dsum / 4;
				shade(nn, colors[f], polys[n].color);
				n++;
			}
		qaws_surface_destroy(srf);
	}
	draw_polys(s, polys, n, 0.5);
	free(polys);
}

static void app_mesh_patches(char const* obj_path)
{
	static patch_fit indep[MP_PATCHES], admm[MP_PATCHES];
	double rms0[MP_OUTER], gap0[MP_OUTER], rms1[MP_OUTER], gap1[MP_OUTER];
	mesh m;
	int nodes, f;
	svg s;
	char buf[256];

	if (!obj_path || !mesh_load_obj(obj_path, &m))
		mesh_blob(&m);
	mp_setup(&m, indep, &nodes);
	mp_setup(&m, admm, &nodes);
	mp_solve(&m, indep, nodes, 0.0, rms0, gap0);
	mp_solve(&m, admm, nodes, 0.03, rms1, gap1);

	svg_open(&s, "showcase/app3_mesh_patches.svg", 1240, 800, "Triangulated mesh to 6 bicubic patches with ADMM",
		"Each patch fits its points independently (exact least squares at foot points + thin plate); shared seam "
		"control points are ADMM consensus variables, so the patch set closes up without a global solve.");
	{
		view3 va = { 215, 320, 135, 0.65, 0.42 }, vb = va, vc = va;
		viewport lv = { 40, 560, 520, 210, 0, 0, 0, 0 };
		vb.cx = 620;
		vc.cx = 1025;
		fprintf(s.f, "<text x=\"20\" y=\"92\" font-size=\"13\" font-weight=\"600\" fill=\"#24292f\">input mesh: %d vertices, %d triangles</text>\n", m.nv, m.nt);
		render_mesh(&s, &va, &m);
		fprintf(s.f, "<text x=\"425\" y=\"92\" font-size=\"13\" font-weight=\"600\" fill=\"#cf222e\">independent patch fits: open, overlapping seams</text>\n");
		render_patches(&s, &vb, indep, 22, 0, &m);
		fprintf(s.f, "<text x=\"830\" y=\"92\" font-size=\"13\" font-weight=\"600\" fill=\"#0969da\">ADMM: %d control points, closed seams</text>\n", MP_PATCHES * MP_CP);
		render_patches(&s, &vc, admm, 22, 0, &m);

		svg_panel(&s, &lv, "log10 seam gap (red: independent, blue: ADMM) and RMS (dashed)");
		{
			double xy0[2 * MP_OUTER], xy1[2 * MP_OUTER], xr0[2 * MP_OUTER], xr1[2 * MP_OUTER], lo = 1e300, hi = -1e300;
			int it;
			for (it = 0; it < MP_OUTER; it++)
			{
				double v[4] = { gap0[it], gap1[it], rms0[it], rms1[it] }, l;
				int q;
				for (q = 0; q < 4; q++)
				{
					l = log10(v[q] > 1e-12 ? v[q] : 1e-12);
					if (l < lo) lo = l;
					if (l > hi) hi = l;
				}
			}
			for (it = 0; it < MP_OUTER; it++)
			{
				double x = lv.x0 + 14 + (lv.w - 28) * it / (double)(MP_OUTER - 1);
				xy0[2 * it] = xy1[2 * it] = xr0[2 * it] = xr1[2 * it] = x;
				xy0[2 * it + 1] = lv.y0 + 34 + (lv.h - 54) * (hi - log10(gap0[it] > 1e-12 ? gap0[it] : 1e-12)) / (hi - lo);
				xy1[2 * it + 1] = lv.y0 + 34 + (lv.h - 54) * (hi - log10(gap1[it] > 1e-12 ? gap1[it] : 1e-12)) / (hi - lo);
				xr0[2 * it + 1] = lv.y0 + 34 + (lv.h - 54) * (hi - log10(rms0[it])) / (hi - lo);
				xr1[2 * it + 1] = lv.y0 + 34 + (lv.h - 54) * (hi - log10(rms1[it])) / (hi - lo);
			}
			svg_polyline(&s, xr0, MP_OUTER, "#cf222e", 1.6, 0.8, 1);
			svg_polyline(&s, xr1, MP_OUTER, "#0969da", 1.6, 0.8, 1);
			svg_polyline(&s, xy0, MP_OUTER, "#cf222e", 2.2, 1, 0);
			svg_polyline(&s, xy1, MP_OUTER, "#0969da", 2.6, 1, 0);
			svg_text(&s, lv.x0 + lv.w - 14, lv.y0 + lv.h - 8, 11, "#57606a", "end", "ADMM iteration");
		}
		{
			double tx = 600, ty = 585;
			sprintf(buf, "RMS distance to the mesh: independent %.4f, ADMM %.4f (shape radius ~1)", rms0[MP_OUTER - 1], rms1[MP_OUTER - 1]);
			svg_text(&s, tx, ty, 14, "#24292f", "start", buf);
			sprintf(buf, "max seam gap: independent %.4f, ADMM %.1e", gap0[MP_OUTER - 1], gap1[MP_OUTER - 1]);
			svg_text(&s, tx, ty + 24, 14, "#0969da", "start", buf);
			svg_text(&s, tx, ty + 58, 13, "#57606a", "start", "per ADMM iteration, every patch in parallel:");
			svg_text(&s, tx, ty + 78, 13, "#24292f", "start", "  foot points: Newton on |S(u,v) - X|^2 (surface jets)");
			svg_text(&s, tx, ty + 98, 13, "#24292f", "start", "  x-update: (2A + t H + rho D) P = 2b + rho D (z - u)");
			svg_text(&s, tx, ty + 118, 13, "#24292f", "start", "    A, b from exact basis weights (qaws_surface_local_support)");
			svg_text(&s, tx, ty + 138, 13, "#24292f", "start", "    H: thin-plate Hessian from direct HVPs (qaws_surface_functional_hvp)");
			svg_text(&s, tx, ty + 158, 13, "#57606a", "start", "then z = mean of (x + u) over seam copies, u += x - z");
		}
	}
	svg_close(&s);
	printf("3_mesh_patches: %d points, rms indep %.4f admm %.4f, seam gap indep %.4f admm %.2e\n",
		m.nv, rms0[MP_OUTER - 1], rms1[MP_OUTER - 1], gap0[MP_OUTER - 1], gap1[MP_OUTER - 1]);
	for (f = 0; f < MP_PATCHES; f++)
	{
		free(indep[f].pts); free(indep[f].uv);
		free(admm[f].pts); free(admm[f].uv);
	}
	free(m.v);
	free(m.t);
}

