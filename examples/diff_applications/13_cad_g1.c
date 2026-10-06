/* ================================================================== */
/* 13. CAD reverse engineering: scan to B-spline patches with G1 seams */
/* ================================================================== */

/*
 * The ADMM fit of application 3 closes the six patches of the blob scan
 * (C0 seams) but their tangent planes still break along the seams. This
 * refinement treats the patch set as one model: shared seam control points
 * are single variables (C0 stays exact), and the energy
 *   mean |S(u_i, v_i) - X_i|^2 + t thin plate + g sum_seams w(s) (1 - (n_a . n_b)^2)
 * pulls the unit normals of both patches together at seam samples (G1).
 * The data term pulls back through the batch adjoint, the thin plate
 * through qaws_surface_functional_gradient, the G1 term through the
 * surface adjoint on the S_u, S_v channels:
 *   N = S_u x S_v, n = N / |N|, c = n_a . n_b,
 *   N_bar = (I - n n^T)(-2 c n_other) / |N|,
 *   S_u_bar = S_v x N_bar,  S_v_bar = N_bar x S_u.
 * Zebra stripes (isophotes of a fixed light) show the result: they kink at
 * a G0 seam and run through a G1 one.
 */

#define CG_SEAMS 24
#define CG_SAMPLES 24
static double const g_cg_gamma = 1.0;
#define CG_ITERS 1500

typedef struct cg_seam
{
	int fa, fb;          /* patches */
	int ea, eb;          /* edges: 0 u=0, 1 u=1, 2 v=0, 3 v=1 */
	int flip;            /* b runs the other way */
} cg_seam;

static void cg_edge_uv(int e, double s, double* u, double* v)
{
	*u = e == 0 ? 0 : (e == 1 ? 1 : s);
	*v = e == 2 ? 0 : (e == 3 ? 1 : s);
}

static void cg_eval(qaws_surface const* srf, double u, double v, double* p)
{
	qaws_surface_eval_result r;
	qaws_surface_evaluate(srf, (qaws_scalar)u, (qaws_scalar)v, QAWS_SURFACE_EVAL_POSITION, &r);
	p[0] = r.position.x;
	p[1] = r.position.y;
	p[2] = r.position.z;
}

/* Pairs of patch edges that coincide (by their end and mid points). */
static int cg_find_seams(patch_fit const* pf, cg_seam* seams)
{
	qaws_surface* s[MP_PATCHES];
	int n = 0, fa, fb, ea, eb, f;
	for (f = 0; f < MP_PATCHES; f++)
		s[f] = mp_surface(pf[f].cps);
	for (fa = 0; fa < MP_PATCHES; fa++)
		for (ea = 0; ea < 4; ea++)
			for (fb = fa + 1; fb < MP_PATCHES; fb++)
				for (eb = 0; eb < 4; eb++)
				{
					double pa[3][3], pb[3][3], u, v, d0 = 0, d1 = 0;
					int k, c;
					for (k = 0; k < 3; k++)
					{
						cg_edge_uv(ea, 0.5 * k, &u, &v);
						cg_eval(s[fa], u, v, pa[k]);
						cg_edge_uv(eb, 0.5 * k, &u, &v);
						cg_eval(s[fb], u, v, pb[k]);
					}
					for (k = 0; k < 3; k++)
						for (c = 0; c < 3; c++)
						{
							d0 += fabs(pa[k][c] - pb[k][c]);
							d1 += fabs(pa[k][c] - pb[2 - k][c]);
						}
					if ((d0 < 1e-3 || d1 < 1e-3) && n < CG_SEAMS)
					{
						seams[n].fa = fa;
						seams[n].fb = fb;
						seams[n].ea = ea;
						seams[n].eb = eb;
						seams[n].flip = d1 < d0;
						n++;
					}
				}
	for (f = 0; f < MP_PATCHES; f++)
		qaws_surface_destroy(s[f]);
	return n;
}

static void cg_normal(qaws_surface_jet const* j, double* N, double* len)
{
	qaws_vec3 a = j->d[1], b = j->d[2];
	N[0] = (double)a.y * b.z - (double)a.z * b.y;
	N[1] = (double)a.z * b.x - (double)a.x * b.z;
	N[2] = (double)a.x * b.y - (double)a.y * b.x;
	*len = sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
}

/* Crease angle (degrees) between the patches at seam parameter s. */
static double cg_crease(qaws_surface const* sa, qaws_surface const* sb, cg_seam const* sm, double s)
{
	qaws_surface_jet ja, jb;
	double ua, va, ub, vb, Na[3], Nb[3], la, lb, c;
	cg_edge_uv(sm->ea, s, &ua, &va);
	cg_edge_uv(sm->eb, sm->flip ? 1 - s : s, &ub, &vb);
	qaws_surface_eval_jet(sa, (qaws_scalar)ua, (qaws_scalar)va, QAWS_SJET_ORDER1, &ja);
	qaws_surface_eval_jet(sb, (qaws_scalar)ub, (qaws_scalar)vb, QAWS_SJET_ORDER1, &jb);
	cg_normal(&ja, Na, &la);
	cg_normal(&jb, Nb, &lb);
	c = fabs(Na[0] * Nb[0] + Na[1] * Nb[1] + Na[2] * Nb[2]) / (la * lb);
	return acos(c > 1 ? 1 : c) * 180 / 3.14159265358979;
}

/* Largest and mean crease over the seams, corners excluded (s in [0.1, 0.9]:
   three patches meet at a cube corner, where bicubic G1 is not possible). */
static void cg_creases(patch_fit const* pf, cg_seam const* seams, int ns, double* mx, double* mean)
{
	qaws_surface* s[MP_PATCHES];
	int f, k, i, cnt = 0;
	*mx = 0;
	*mean = 0;
	for (f = 0; f < MP_PATCHES; f++)
		s[f] = mp_surface(pf[f].cps);
	for (k = 0; k < ns; k++)
		for (i = 0; i <= 40; i++)
		{
			double a = cg_crease(s[seams[k].fa], s[seams[k].fb], &seams[k], 0.1 + 0.8 * i / 40.0);
			if (a > *mx) *mx = a;
			*mean += a;
			cnt++;
		}
	*mean /= cnt;
	for (f = 0; f < MP_PATCHES; f++)
		qaws_surface_destroy(s[f]);
}

/* G1 term and its control point gradient (added to grad[f]). */
static double cg_g1_energy(qaws_surface* const* s, cg_seam const* seams, int ns, double gamma, qaws_scalar (*grad)[MP_CP * 3])
{
	double e = 0;
	int k, i;
	for (k = 0; k < ns; k++)
		for (i = 0; i < CG_SAMPLES; i++)
		{
			cg_seam const* sm = &seams[k];
			double sp = (i + 0.5) / CG_SAMPLES, w = sin(3.14159265358979 * sp), uv[2][2], N[2][3], len[2], n[2][3], c;
			qaws_surface_jet j[2];
			int side, a;
			w = gamma * w * w / (ns * CG_SAMPLES);
			cg_edge_uv(sm->ea, sp, &uv[0][0], &uv[0][1]);
			cg_edge_uv(sm->eb, sm->flip ? 1 - sp : sp, &uv[1][0], &uv[1][1]);
			for (side = 0; side < 2; side++)
			{
				qaws_surface_eval_jet(s[side ? sm->fb : sm->fa], (qaws_scalar)uv[side][0], (qaws_scalar)uv[side][1], QAWS_SJET_ORDER1, &j[side]);
				cg_normal(&j[side], N[side], &len[side]);
				for (a = 0; a < 3; a++)
					n[side][a] = N[side][a] / len[side];
			}
			c = n[0][0] * n[1][0] + n[0][1] * n[1][1] + n[0][2] * n[1][2];
			e += w * (1 - c * c);
			for (side = 0; side < 2; side++)
			{
				double const* no = n[1 - side];
				double nb[3], pr, Nbar[3];
				qaws_surface_jet bar;
				qaws_field_view fv;
				qaws_diff_views views;
				qaws_vec3 su = j[side].d[1], sv = j[side].d[2];
				for (a = 0; a < 3; a++)
					nb[a] = -2 * c * no[a] * w;
				pr = nb[0] * n[side][0] + nb[1] * n[side][1] + nb[2] * n[side][2];
				for (a = 0; a < 3; a++)
					Nbar[a] = (nb[a] - pr * n[side][a]) / len[side];
				memset(&bar, 0, sizeof(bar));
				/* S_u_bar = S_v x N_bar, S_v_bar = N_bar x S_u */
				bar.d[1].x = (qaws_scalar)(sv.y * Nbar[2] - sv.z * Nbar[1]);
				bar.d[1].y = (qaws_scalar)(sv.z * Nbar[0] - sv.x * Nbar[2]);
				bar.d[1].z = (qaws_scalar)(sv.x * Nbar[1] - sv.y * Nbar[0]);
				bar.d[2].x = (qaws_scalar)(Nbar[1] * su.z - Nbar[2] * su.y);
				bar.d[2].y = (qaws_scalar)(Nbar[2] * su.x - Nbar[0] * su.z);
				bar.d[2].z = (qaws_scalar)(Nbar[0] * su.y - Nbar[1] * su.x);
				bar.channels = QAWS_SJET_U | QAWS_SJET_V;
				views = one_field(&fv, QAWS_FIELD_CONTROL_POINTS, grad[side ? sm->fb : sm->fa], MP_CP, 3);
				qaws_surface_eval_adjoint(NULL, s[side ? sm->fb : sm->fa], (qaws_scalar)uv[side][0], (qaws_scalar)uv[side][1],
					QAWS_SJET_U | QAWS_SJET_V, &bar, &views, NULL, NULL);
			}
		}
	return e;
}

/* Full energy and gradient on the patch control points. */
static double cg_energy(patch_fit* pf, mesh const* m, cg_seam const* seams, int ns, double gamma, qaws_scalar (*grad)[MP_CP * 3],
	double* out_rms)
{
	static qaws_surface_jet bars[MP_GRID * MP_GRID * 2];
	static qaws_scalar us[MP_GRID * MP_GRID * 2], vs[MP_GRID * MP_GRID * 2];
	qaws_surface* s[MP_PATCHES];
	double e = 0, d2 = 0, tp = 2e-6;
	int f, i, npts = 0;
	for (f = 0; f < MP_PATCHES; f++)
	{
		s[f] = mp_surface(pf[f].cps);
		npts += pf[f].n;
	}
	for (f = 0; f < MP_PATCHES; f++)
	{
		qaws_field_view fv;
		qaws_diff_views views = one_field(&fv, QAWS_FIELD_CONTROL_POINTS, grad[f], MP_CP, 3);
		qaws_scalar value = 0;
		memset(grad[f], 0, sizeof(grad[f]));
		for (i = 0; i < pf[f].n && i < MP_GRID * MP_GRID * 2; i++)
		{
			qaws_surface_jet j;
			double const* x = &m->v[3 * pf[f].pts[i]];
			double r[3];
			us[i] = pf[f].uv[2 * i];
			vs[i] = pf[f].uv[2 * i + 1];
			qaws_surface_eval_jet(s[f], us[i], vs[i], QAWS_SJET_P, &j);
			r[0] = j.d[0].x - x[0];
			r[1] = j.d[0].y - x[1];
			r[2] = j.d[0].z - x[2];
			d2 += r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
			memset(&bars[i], 0, sizeof(bars[i]));
			bars[i].d[0].x = (qaws_scalar)(2 * r[0] / npts);
			bars[i].d[0].y = (qaws_scalar)(2 * r[1] / npts);
			bars[i].d[0].z = (qaws_scalar)(2 * r[2] / npts);
			bars[i].channels = QAWS_SJET_P;
		}
		qaws_surface_eval_batch_adjoint(NULL, s[f], us, vs, (unsigned int)i, QAWS_SJET_P, bars, &views, NULL, NULL);
		{
			/* thin plate, scaled into the gradient */
			qaws_scalar tg[MP_CP * 3];
			qaws_field_view ft;
			qaws_diff_views vt = one_field(&ft, QAWS_FIELD_CONTROL_POINTS, tg, MP_CP, 3);
			int k;
			memset(tg, 0, sizeof(tg));
			qaws_surface_functional_gradient(NULL, s[f], QAWS_FUNCTIONAL_THIN_PLATE, 4, &vt, &value);
			for (k = 0; k < MP_CP * 3; k++)
				grad[f][k] += (qaws_scalar)(tp * tg[k]);
			e += tp * value;
		}
	}
	e += d2 / npts;
	e += cg_g1_energy(s, seams, ns, gamma, grad);
	for (f = 0; f < MP_PATCHES; f++)
		qaws_surface_destroy(s[f]);
	if (out_rms)
		*out_rms = sqrt(d2 / npts);
	return e;
}

/* Seam nodes made exactly shared (the mean of their copies). */
static void cg_share_nodes(patch_fit* pf, int nodes)
{
	double* acc = (double*)calloc((size_t)nodes * 4, sizeof(double));
	int f, k, c;
	for (f = 0; f < MP_PATCHES; f++)
		for (k = 0; k < MP_CP; k++)
			if (pf[f].node[k] >= 0)
			{
				for (c = 0; c < 3; c++)
					acc[4 * pf[f].node[k] + c] += pf[f].cps[3 * k + c];
				acc[4 * pf[f].node[k] + 3] += 1;
			}
	for (f = 0; f < MP_PATCHES; f++)
		for (k = 0; k < MP_CP; k++)
			if (pf[f].node[k] >= 0)
				for (c = 0; c < 3; c++)
					pf[f].cps[3 * k + c] = (qaws_scalar)(acc[4 * pf[f].node[k] + c] / acc[4 * pf[f].node[k] + 3]);
	free(acc);
}

/* Adam on the shared variables: interior control points per patch, seam
   nodes once (their gradients summed over copies). */
static void cg_refine(patch_fit* pf, mesh const* m, int nodes, cg_seam const* seams, int ns, double gamma, double* hist_e,
	double* hist_crease)
{
	static qaws_scalar grad[MP_PATCHES][MP_CP * 3];
	int nv = MP_PATCHES * MP_CP + nodes, f, k, c, it;
	qaws_scalar* x = (qaws_scalar*)calloc((size_t)nv * 3, sizeof(qaws_scalar));
	qaws_scalar* g = (qaws_scalar*)calloc((size_t)nv * 3, sizeof(qaws_scalar));
	double* mo = (double*)calloc((size_t)nv * 3, sizeof(double));
	double* vo = (double*)calloc((size_t)nv * 3, sizeof(double));
	for (it = 0; it < CG_ITERS; it++)
	{
		double e, mx, mean;
		if (it % 25 == 0)
			for (f = 0; f < MP_PATCHES; f++)
			{
				qaws_surface* s = mp_surface(pf[f].cps);
				mp_project(&pf[f], m, s, 3);
				qaws_surface_destroy(s);
			}
		e = cg_energy(pf, m, seams, ns, gamma, grad, NULL);
		if (hist_e) hist_e[it] = e;
		if (hist_crease && it % 10 == 0)
		{
			cg_creases(pf, seams, ns, &mx, &mean);
			hist_crease[it / 10] = mx;
		}
		/* gather into the shared variables */
		memset(g, 0, sizeof(qaws_scalar) * (size_t)nv * 3);
		for (f = 0; f < MP_PATCHES; f++)
			for (k = 0; k < MP_CP; k++)
			{
				int var = pf[f].node[k] >= 0 ? MP_PATCHES * MP_CP + pf[f].node[k] : f * MP_CP + k;
				for (c = 0; c < 3; c++)
				{
					g[3 * var + c] += grad[f][3 * k + c];
					x[3 * var + c] = pf[f].cps[3 * k + c];
				}
			}
		/* Adam (flat vectors larger than the helper's) */
		{
			double lr = it < 800 ? 2e-3 : 4e-4, b1t = 1 - pow(0.9, it + 1), b2t = 1 - pow(0.999, it + 1);
			for (k = 0; k < 3 * nv; k++)
			{
				mo[k] = 0.9 * mo[k] + 0.1 * g[k];
				vo[k] = 0.999 * vo[k] + 0.001 * g[k] * g[k];
				x[k] = (qaws_scalar)(x[k] - lr * (mo[k] / b1t) / (sqrt(vo[k] / b2t) + 1e-12));
			}
		}
		for (f = 0; f < MP_PATCHES; f++)
			for (k = 0; k < MP_CP; k++)
			{
				int var = pf[f].node[k] >= 0 ? MP_PATCHES * MP_CP + pf[f].node[k] : f * MP_CP + k;
				for (c = 0; c < 3; c++)
					pf[f].cps[3 * k + c] = x[3 * var + c];
			}
	}
	free(x);
	free(g);
	free(mo);
	free(vo);
}

/* Patches shaded with zebra stripes: isophotes of a light direction. */
static void cg_render_zebra(svg* s, view3 const* v, patch_fit const* pf, int grid)
{
	poly3* polys = (poly3*)malloc(sizeof(poly3) * MP_PATCHES * (size_t)grid * grid);
	double L[3] = { 0.30, -0.80, 0.52 };
	int f, i, j, n = 0;
	for (f = 0; f < MP_PATCHES; f++)
	{
		qaws_surface* srf = mp_surface(pf[f].cps);
		for (i = 0; i < grid; i++)
			for (j = 0; j < grid; j++)
			{
				qaws_surface_eval_result r[4], mid;
				double p[3], d, dsum = 0, nn[3], cy = cos(v->yaw), sy_ = sin(v->yaw), cp = cos(v->pitch), sp = sin(v->pitch), ang, k;
				int c, stripe;
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
				ang = acos((nn[0] * L[0] + nn[1] * L[1] + nn[2] * L[2]) / sqrt(L[0] * L[0] + L[1] * L[1] + L[2] * L[2]));
				stripe = (int)floor(ang * 14 / 3.14159265358979);
				k = 0.55 + 0.45 * (-(cp * (sy_ * nn[0] + cy * nn[1]) - sp * nn[2]));
				if (stripe & 1)
					sprintf(polys[n].color, "rgb(%d,%d,%d)", (int)(30 * k), (int)(32 * k), (int)(40 * k));
				else
					sprintf(polys[n].color, "rgb(%d,%d,%d)", (int)(250 * k), (int)(250 * k), (int)(252 * k));
				polys[n].n = 4;
				polys[n].depth = dsum / 4;
				n++;
			}
		qaws_surface_destroy(srf);
	}
	draw_polys(s, polys, n, 0.35);
	free(polys);
}

static void app_cad_g1(void)
{
	static patch_fit pf[MP_PATCHES], before[MP_PATCHES];
	static double hist_e[CG_ITERS], hist_cr[CG_ITERS / 10 + 1];
	double rms0[MP_OUTER], gap0[MP_OUTER], mx0, mean0, mx1, mean1, rms_a, rms_b;
	static qaws_scalar grad[MP_PATCHES][MP_CP * 3];
	cg_seam seams[CG_SEAMS];
	mesh m;
	int nodes, ns, f;
	svg s;
	char buf[256];
	g_blob_amp = 1.0;
	mesh_blob(&m);
	mp_setup(&m, pf, &nodes);
	mp_solve(&m, pf, nodes, 0.03, rms0, gap0);
	cg_share_nodes(pf, nodes);
	ns = cg_find_seams(pf, seams);
	for (f = 0; f < MP_PATCHES; f++)
		before[f] = pf[f];
	cg_creases(pf, seams, ns, &mx0, &mean0);
	cg_energy(pf, &m, seams, ns, 0, grad, &rms_b);
	cg_refine(pf, &m, nodes, seams, ns, g_cg_gamma, hist_e, hist_cr);
	cg_creases(pf, seams, ns, &mx1, &mean1);
	cg_energy(pf, &m, seams, ns, 0, grad, &rms_a);
	printf("13_cad_g1: %d seams; crease max %.2f deg -> %.2f deg, mean %.2f -> %.3f deg; rms %.4f -> %.4f; seam gap %.1e\n",
		ns, mx0, mx1, mean0, mean1, rms_b, rms_a, mp_seam_gap(pf, nodes));

	svg_open(&s, "showcase/app13_cad_g1.svg", 1240, 720, "CAD reverse engineering: scan to six B-spline patches with G1 seams",
		"Zebra stripes kink where tangent planes break. ADMM closes the seams (C0); a joint refinement adds 1 - (n_a . n_b)^2 at seam samples (surface adjoint).");
	{
		view3 va = { 300, 360, 190, 0.65, 0.42 }, vb = va;
		vb.cx = 900;
		sprintf(buf, "ADMM fit (C0): max crease %.2f deg, mean %.2f deg", mx0, mean0);
		fprintf(s.f, "<text x=\"60\" y=\"100\" font-size=\"14\" font-weight=\"600\" fill=\"#cf222e\">%s</text>\n", buf);
		cg_render_zebra(&s, &va, before, 64);
		sprintf(buf, "G1 refinement: max crease %.2f deg, mean %.3f deg", mx1, mean1);
		fprintf(s.f, "<text x=\"660\" y=\"100\" font-size=\"14\" font-weight=\"600\" fill=\"#1a7f37\">%s</text>\n", buf);
		cg_render_zebra(&s, &vb, pf, 64);
		sprintf(buf, "RMS distance to the scan: %.4f -> %.4f (shape radius ~1); seams stay closed (shared control points); cube corners excluded (valence 3)",
			rms_b, rms_a);
		svg_text(&s, 60, 690, 13, "#57606a", "start", buf);
	}
	svg_close(&s);
	for (f = 0; f < MP_PATCHES; f++)
	{
		free(pf[f].pts);
		free(pf[f].uv);
	}
	free(m.v);
	free(m.t);
}
