/* ================================================================== */
/* 12. Shape optimization: 2D vesicles by Newton-CG on exact HVPs      */
/* ================================================================== */

/*
 * A closed curve of fixed length L0 = 2 pi and enclosed area A0 = nu pi
 * (reduced area nu = 4 pi A / L^2) minimizes its bending energy
 * integral kappa^2 ds: the 2D vesicle problem. Ellipses at nu near 1 give
 * way to peanuts and biconcave, red-blood-cell-like sections as nu drops.
 *
 * The curve is a periodic cubic B-spline (uniform knots, first three
 * control points repeated). The objective is an augmented Lagrangian
 *   K + alpha B + lambda_L c_L + mu/2 c_L^2 + lambda_A c_A + mu/2 c_A^2,
 * with K = integral kappa^2 ds, B = integral |C''|^2 dt (keeps the
 * parameterization even), c_L = L - L0, c_A = A - A0. K, B and L and their
 * gradients and Hessian-vector products come from qaws_curve_functional_*;
 * the area A = 1/2 integral (x y' - y x') dt is bilinear in the jets, so its
 * gradient and HVP are one batch adjoint of the jets and of their tangents.
 * Newton-CG: truncated conjugate gradients on H p = -g with exact HVPs,
 * then a backtracking line search; lambda updates between solves.
 */

#define VS_N 36
#define VS_Q 8
#define VS_RULE 4
#define VS_CG_MAX (6 * VS_N)

typedef struct vs_state
{
	double L0, A0, alpha, beta, s2, mu, lam_L, lam_A;
	/* constraint values and gradients at the current Newton point, for
	   the Hessian-vector products of the CG solve */
	double cl, ca, gL[2 * VS_N], gA[2 * VS_N];
} vs_state;

static qaws_scalar g_vs_knots[VS_N + 7];

static qaws_curve* vs_curve(qaws_scalar const* p)
{
	qaws_scalar ext[(VS_N + 3) * 2];
	qaws_bspline_desc d;
	qaws_curve* c = NULL;
	int i;
	for (i = 0; i < VS_N + 7; i++)
		g_vs_knots[i] = (qaws_scalar)i;
	for (i = 0; i < VS_N + 3; i++)
	{
		ext[2 * i] = p[2 * (i % VS_N)];
		ext[2 * i + 1] = p[2 * (i % VS_N) + 1];
	}
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D;
	d.degree = 3;
	d.control_points = ext;
	d.control_point_count = VS_N + 3;
	d.knots = g_vs_knots;
	d.knot_count = VS_N + 7;
	qaws_curve_create_bspline(&d, &c);
	return c;
}

static void vs_unfold(qaws_scalar const* d, qaws_scalar* ext)
{
	int i;
	for (i = 0; i < VS_N + 3; i++)
	{
		ext[2 * i] = d[2 * (i % VS_N)];
		ext[2 * i + 1] = d[2 * (i % VS_N) + 1];
	}
}

static void vs_fold_add(qaws_scalar const* ext, double s, double* out)
{
	int i;
	for (i = 0; i < VS_N + 3; i++)
	{
		out[2 * (i % VS_N)] += s * ext[2 * i];
		out[2 * (i % VS_N) + 1] += s * ext[2 * i + 1];
	}
}

/* Value, gradient and (if d) HVP of one library functional, on the unique points. */
static double vs_functional(qaws_curve const* c, qaws_curve_functional f, double scale, double const* d, double* grad, double* hv,
	double* raw_grad)
{
	qaws_scalar ext[(VS_N + 3) * 2], dext[(VS_N + 3) * 2], dd[VS_N * 2];
	qaws_field_view fv, fd;
	qaws_diff_views vg, vd;
	qaws_scalar value = 0;
	int i;
	memset(ext, 0, sizeof(ext));
	vg = one_field(&fv, QAWS_FIELD_CONTROL_POINTS, ext, VS_N + 3, 2);
	if (grad || raw_grad)
		qaws_curve_functional_gradient(NULL, c, f, VS_RULE, &vg, &value);
	else if (!(d && hv))
		qaws_curve_functional_eval(NULL, c, f, VS_RULE, NULL, &value, NULL, NULL);
	if (grad)
		vs_fold_add(ext, scale, grad);
	if (raw_grad)
	{
		memset(raw_grad, 0, sizeof(double) * 2 * VS_N);
		vs_fold_add(ext, 1, raw_grad);
	}
	if (d && hv)
	{
		for (i = 0; i < 2 * VS_N; i++)
			dd[i] = (qaws_scalar)d[i];
		vs_unfold(dd, dext);
		vd = one_field(&fd, QAWS_FIELD_CONTROL_POINTS, dext, VS_N + 3, 2);
		memset(ext, 0, sizeof(ext));
		qaws_curve_functional_hvp(NULL, c, f, VS_RULE, &vd, &vg);
		vs_fold_add(ext, scale, hv);
	}
	return value;
}

/* Enclosed area 1/2 integral (x y' - y x') dt, its gradient and (if d) H d. */
static double vs_area(qaws_curve const* c, double const* d, double* grad, double* hv)
{
	static qaws_scalar ts[VS_N * VS_Q];
	static qaws_curve_jet_2d prim[VS_N * VS_Q], tan[VS_N * VS_Q], bar[VS_N * VS_Q];
	static double ws[VS_N * VS_Q];
	qaws_scalar ext[(VS_N + 3) * 2], dext[(VS_N + 3) * 2], dd[VS_N * 2];
	qaws_field_view fv, fd;
	qaws_diff_views vg, vd;
	unsigned int m = VS_N * VS_Q, k;
	double a = 0;
	for (k = 0; k < m; k++)
	{
		static double const gx[4] = { 0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363 };
		static double const gw[4] = { 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763 };
		int q = (int)(k % VS_Q), j = q < 4 ? 3 - q : q - 4;
		ts[k] = (qaws_scalar)(3 + k / VS_Q + 0.5 * (1 + (q < 4 ? -gx[j] : gx[j])));
		ws[k] = 0.5 * gw[j];
	}
	vd = one_field(&fd, QAWS_FIELD_CONTROL_POINTS, dext, VS_N + 3, 2);
	if (d)
	{
		for (k = 0; k < 2 * VS_N; k++)
			dd[k] = (qaws_scalar)d[k];
		vs_unfold(dd, dext);
	}
	qaws_curve_eval_batch_tangent_2d(NULL, c, ts, NULL, m, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, d ? &vd : NULL, prim, tan);
	for (k = 0; k < m; k++)
	{
		qaws_vec2 p = prim[k].d[0], q = prim[k].d[1];
		a += 0.5 * ws[k] * (p.x * q.y - p.y * q.x);
		memset(&bar[k], 0, sizeof(bar[k]));
		bar[k].d[0].x = (qaws_scalar)(0.5 * ws[k] * q.y);
		bar[k].d[0].y = (qaws_scalar)(-0.5 * ws[k] * q.x);
		bar[k].d[1].x = (qaws_scalar)(-0.5 * ws[k] * p.y);
		bar[k].d[1].y = (qaws_scalar)(0.5 * ws[k] * p.x);
		bar[k].channels = QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1;
	}
	vg = one_field(&fv, QAWS_FIELD_CONTROL_POINTS, ext, VS_N + 3, 2);
	if (grad)
	{
		memset(ext, 0, sizeof(ext));
		qaws_curve_eval_batch_adjoint_2d(NULL, c, ts, m, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, bar, &vg, NULL);
		memset(grad, 0, sizeof(double) * 2 * VS_N);
		vs_fold_add(ext, 1, grad);
	}
	if (d && hv)
	{
		/* the area is bilinear in the jets: the gradient bars of the jet tangents */
		for (k = 0; k < m; k++)
		{
			qaws_vec2 p = tan[k].d[0], q = tan[k].d[1];
			bar[k].d[0].x = (qaws_scalar)(0.5 * ws[k] * q.y);
			bar[k].d[0].y = (qaws_scalar)(-0.5 * ws[k] * q.x);
			bar[k].d[1].x = (qaws_scalar)(-0.5 * ws[k] * p.y);
			bar[k].d[1].y = (qaws_scalar)(0.5 * ws[k] * p.x);
		}
		memset(ext, 0, sizeof(ext));
		qaws_curve_eval_batch_adjoint_2d(NULL, c, ts, m, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, bar, &vg, NULL);
		memset(hv, 0, sizeof(double) * 2 * VS_N);
		vs_fold_add(ext, 1, hv);
	}
	return a;
}

/*
 * Gauge term R = integral (|C'|^2 - s2)^2 dt (s2 the squared speed of an
 * arc-length parameterization, held fixed): it vanishes when the curve is
 * parameterized by arc length, so it leaves the shape alone, but pins the
 * tangential sliding to which the bending energy is blind. Polynomial in
 * C': with r = |C'|^2 - s2, dR/dC' = 4 w r C', and along a tangent
 * d(dR/dC') = 4 w (2 (C' . C'_dot) C' + r C'_dot).
 */
static double vs_speed(qaws_curve const* c, double s2, double const* d, double* grad, double* hv)
{
	static qaws_scalar ts[VS_N * VS_Q];
	static qaws_curve_jet_2d prim[VS_N * VS_Q], tan[VS_N * VS_Q], bar[VS_N * VS_Q];
	static double ws[VS_N * VS_Q];
	qaws_scalar ext[(VS_N + 3) * 2], dext[(VS_N + 3) * 2], dd[VS_N * 2];
	qaws_field_view fv, fd;
	qaws_diff_views vg, vd;
	unsigned int m = VS_N * VS_Q, k;
	double R = 0;
	for (k = 0; k < m; k++)
	{
		static double const gx[4] = { 0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363 };
		static double const gw[4] = { 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763 };
		int q = (int)(k % VS_Q), j = q < 4 ? 3 - q : q - 4;
		ts[k] = (qaws_scalar)(3 + k / VS_Q + 0.5 * (1 + (q < 4 ? -gx[j] : gx[j])));
		ws[k] = 0.5 * gw[j];
	}
	vd = one_field(&fd, QAWS_FIELD_CONTROL_POINTS, dext, VS_N + 3, 2);
	if (d)
	{
		for (k = 0; k < 2 * VS_N; k++)
			dd[k] = (qaws_scalar)d[k];
		vs_unfold(dd, dext);
	}
	qaws_curve_eval_batch_tangent_2d(NULL, c, ts, NULL, m, QAWS_EVAL_FLAG_D1, d ? &vd : NULL, prim, tan);
	vg = one_field(&fv, QAWS_FIELD_CONTROL_POINTS, ext, VS_N + 3, 2);
	for (k = 0; k < m; k++)
	{
		qaws_vec2 q = prim[k].d[1];
		double r = (double)q.x * q.x + (double)q.y * q.y - s2;
		R += ws[k] * r * r;
		memset(&bar[k], 0, sizeof(bar[k]));
		bar[k].channels = QAWS_EVAL_FLAG_D1;
		if (d && hv)
		{
			qaws_vec2 qd = tan[k].d[1];
			double rd = 2 * ((double)q.x * qd.x + (double)q.y * qd.y);
			bar[k].d[1].x = (qaws_scalar)(4 * ws[k] * (rd * q.x + r * qd.x));
			bar[k].d[1].y = (qaws_scalar)(4 * ws[k] * (rd * q.y + r * qd.y));
		}
		else
		{
			bar[k].d[1].x = (qaws_scalar)(4 * ws[k] * r * q.x);
			bar[k].d[1].y = (qaws_scalar)(4 * ws[k] * r * q.y);
		}
	}
	if (grad || (d && hv))
	{
		memset(ext, 0, sizeof(ext));
		qaws_curve_eval_batch_adjoint_2d(NULL, c, ts, m, QAWS_EVAL_FLAG_D1, bar, &vg, NULL);
		vs_fold_add(ext, 1, (d && hv) ? hv : grad);
	}
	return R;
}

/* Objective and (if grad) its gradient; the gradient pass also caches the
   constraint values and gradients in s for vs_hvp. */
static double vs_objective(vs_state* s, qaws_scalar const* p, double* grad, double* out_cl, double* out_ca, double* out_k)
{
	qaws_curve* c = vs_curve(p);
	double gL[2 * VS_N], gA[2 * VS_N], gR[2 * VS_N], K, B, L, A, R, cl, ca;
	int i;
	if (grad) memset(grad, 0, sizeof(double) * 2 * VS_N);
	K = vs_functional(c, QAWS_FUNCTIONAL_CURVATURE_SQUARED, 1, NULL, grad, NULL, NULL);
	B = vs_functional(c, QAWS_FUNCTIONAL_BENDING, s->alpha, NULL, grad, NULL, NULL);
	L = vs_functional(c, QAWS_FUNCTIONAL_LENGTH, 0, NULL, NULL, NULL, grad ? gL : NULL);
	A = vs_area(c, NULL, grad ? gA : NULL, NULL);
	memset(gR, 0, sizeof(gR));
	R = vs_speed(c, s->s2, NULL, grad ? gR : NULL, NULL);
	cl = L - s->L0;
	ca = A - s->A0;
	if (grad)
	{
		for (i = 0; i < 2 * VS_N; i++)
			grad[i] += s->beta * gR[i] + (s->lam_L + s->mu * cl) * gL[i] + (s->lam_A + s->mu * ca) * gA[i];
		s->cl = cl;
		s->ca = ca;
		memcpy(s->gL, gL, sizeof(gL));
		memcpy(s->gA, gA, sizeof(gA));
	}
	qaws_curve_destroy(c);
	if (out_cl) *out_cl = cl;
	if (out_ca) *out_ca = ca;
	if (out_k) *out_k = K;
	return K + s->alpha * B + s->beta * R + s->lam_L * cl + 0.5 * s->mu * cl * cl + s->lam_A * ca + 0.5 * s->mu * ca * ca;
}

/* Hessian of the objective times d at the point of the last gradient pass:
   H_K d + alpha H_B d + mu (gL gL^T + gA gA^T) d + (lambda + mu c) (H_L, H_A) d. */
static void vs_hvp(vs_state const* s, qaws_scalar const* p, double const* d, double* hv)
{
	qaws_curve* c = vs_curve(p);
	double hL[2 * VS_N], hA[2 * VS_N], hR[2 * VS_N], dl = 0, da = 0;
	int i;
	memset(hv, 0, sizeof(double) * 2 * VS_N);
	memset(hL, 0, sizeof(hL));
	vs_functional(c, QAWS_FUNCTIONAL_CURVATURE_SQUARED, 1, d, NULL, hv, NULL);
	vs_functional(c, QAWS_FUNCTIONAL_BENDING, s->alpha, d, NULL, hv, NULL);
	vs_functional(c, QAWS_FUNCTIONAL_LENGTH, 1, d, NULL, hL, NULL);
	vs_area(c, d, NULL, hA);
	memset(hR, 0, sizeof(hR));
	vs_speed(c, s->s2, d, NULL, hR);
	for (i = 0; i < 2 * VS_N; i++)
	{
		dl += s->gL[i] * d[i];
		da += s->gA[i] * d[i];
	}
	for (i = 0; i < 2 * VS_N; i++)
		hv[i] += s->beta * hR[i] + s->mu * (dl * s->gL[i] + da * s->gA[i]) + (s->lam_L + s->mu * s->cl) * hL[i] + (s->lam_A + s->mu * s->ca) * hA[i];
	qaws_curve_destroy(c);
}

static double vs_dot(double const* a, double const* b)
{
	double r = 0;
	int i;
	for (i = 0; i < 2 * VS_N; i++)
		r += a[i] * b[i];
	return r;
}

/* Truncated CG on H p = -g (stops at negative curvature). */
static int vs_newton_direction(vs_state const* s, qaws_scalar const* p, double const* g, double* step)
{
	double r[2 * VS_N], dir[2 * VS_N], hd[2 * VS_N], gn = sqrt(vs_dot(g, g)), tol, rr;
	int i, it;
	tol = gn * (gn < 0.5 ? sqrt(gn) : 0.5);
	memset(step, 0, sizeof(double) * 2 * VS_N);
	for (i = 0; i < 2 * VS_N; i++)
	{
		r[i] = -g[i];
		dir[i] = r[i];
	}
	rr = vs_dot(r, r);
	for (it = 0; it < VS_CG_MAX; it++)
	{
		double curv, al, rr2;
		vs_hvp(s, p, dir, hd);
		curv = vs_dot(dir, hd);
		if (curv <= 1e-12 * vs_dot(dir, dir))
		{
			if (it == 0)
				for (i = 0; i < 2 * VS_N; i++)
					step[i] = -g[i];
			break;
		}
		al = rr / curv;
		for (i = 0; i < 2 * VS_N; i++)
		{
			step[i] += al * dir[i];
			r[i] -= al * hd[i];
		}
		rr2 = vs_dot(r, r);
		if (sqrt(rr2) <= tol)
			return it + 1;
		for (i = 0; i < 2 * VS_N; i++)
			dir[i] = r[i] + rr2 / rr * dir[i];
		rr = rr2;
	}
	return it;
}

/* Minimizes from p; records |gradient| per Newton step. */
static int vs_solve(vs_state* s, qaws_scalar* p, double* gnorms, int cap, int* cg_total)
{
	double g[2 * VS_N], step[2 * VS_N];
	qaws_scalar trial[2 * VS_N];
	int outer, iters = 0, i;
	*cg_total = 0;
	for (outer = 0; outer < 6; outer++)
	{
		int it;
		double cl, ca;
		for (it = 0; it < 40 && iters < cap; it++)
		{
			double f = vs_objective(s, p, g, NULL, NULL, NULL), gn = sqrt(vs_dot(g, g)), t = 1, slope;
			gnorms[iters++] = gn;
			if (gn < 1e-8)
				break;
			*cg_total += vs_newton_direction(s, p, g, step);
			slope = vs_dot(g, step);
			if (slope >= 0)
			{
				for (i = 0; i < 2 * VS_N; i++)
					step[i] = -g[i];
				slope = -gn * gn;
			}
			for (;;)
			{
				for (i = 0; i < 2 * VS_N; i++)
					trial[i] = (qaws_scalar)(p[i] + t * step[i]);
				double ft = vs_objective(s, trial, NULL, NULL, NULL, NULL);
				/* Armijo, or a full Newton step whose decrease is below the
				   round-off of f (near the solution) */
				if (ft <= f + 1e-4 * t * slope || (t == 1 && ft <= f + 1e-14 * fabs(f)) || t < 1e-8)
					break;
				t *= 0.5;
			}
			memcpy(p, trial, sizeof(trial));
		}
		vs_objective(s, p, NULL, &cl, &ca, NULL);
		s->lam_L += s->mu * cl;
		s->lam_A += s->mu * ca;
		if (fabs(cl) < 1e-9 && fabs(ca) < 1e-9)
			break;
	}
	return iters;
}

static void app_vesicles(void)
{
	static double const nus[4] = { 0.9, 0.75, 0.6, 0.5 };
	static char const* const colors[4] = { "#0969da", "#1a7f37", "#bc4c00", "#8250df" };
	double gn[4][256];
	int iters[4];
	svg s;
	int v;
	svg_open(&s, "showcase/app12_vesicles.svg", 1500, 640, "Shape optimization: 2D vesicles by Newton-CG on exact Hessian-vector products",
		"Bending energy of a closed cubic B-spline (36 control points) at length 2 pi and reduced area nu = 4 pi A / L^2; exact HVPs of the functionals and the area.");
	for (v = 0; v < 4; v++)
	{
		vs_state st;
		qaws_scalar p[2 * VS_N];
		viewport vp = { 20 + v * 300, 70, 290, 330, -1.9, 1.9, -1.9, 1.9 };
		double cl, ca, K, xy[2 * 401];
		char buf[160];
		int i, cg = 0;
		qaws_curve* c;
		qaws_range r;
		/* start: an ellipse of the target area, slightly lopsided */
		double a = 1.0 + 0.9 * (1 - nus[v]), b = nus[v] / a;
		for (i = 0; i < VS_N; i++)
		{
			double th = 2 * 3.14159265358979 * (i + 0.5) / VS_N;
			p[2 * i] = (qaws_scalar)(a * cos(th) * (1 + 0.04 * sin(3 * th)));
			p[2 * i + 1] = (qaws_scalar)(b * sin(th));
		}
		st.L0 = 2 * 3.14159265358979;
		st.A0 = nus[v] * 3.14159265358979;
		st.alpha = 1e-6;
		st.beta = 10;
		st.s2 = (2 * 3.14159265358979 / VS_N) * (2 * 3.14159265358979 / VS_N);
		st.mu = 50;
		st.lam_L = 0;
		st.lam_A = 0;
		{
			double d[2 * VS_N], g0[2 * VS_N], gp[2 * VS_N], gm[2 * VS_N], hv[2 * VS_N], h = 1e-5, err = 0, nrm = 0;
			qaws_scalar pp[2 * VS_N], pm[2 * VS_N];
			for (i = 0; i < 2 * VS_N; i++)
				d[i] = sin(1.7 * i + 0.3);
			for (i = 0; i < 2 * VS_N; i++)
			{
				pp[i] = (qaws_scalar)(p[i] + h * d[i]);
				pm[i] = (qaws_scalar)(p[i] - h * d[i]);
			}
			vs_objective(&st, pp, gp, NULL, NULL, NULL);
			vs_objective(&st, pm, gm, NULL, NULL, NULL);
			vs_objective(&st, p, g0, NULL, NULL, NULL);
			vs_hvp(&st, p, d, hv);
			for (i = 0; i < 2 * VS_N; i++)
			{
				double fd = (gp[i] - gm[i]) / (2 * h);
				err += (hv[i] - fd) * (hv[i] - fd);
				nrm += fd * fd;
			}
			/* the exact HVP against central differences of the gradient */
			printf("12_vesicles: nu %.2f: HVP vs finite differences of the gradient, relative error %.2e\n", nus[v], sqrt(err / nrm));
		}
		iters[v] = vs_solve(&st, p, gn[v], 256, &cg);
		vs_objective(&st, p, NULL, &cl, &ca, &K);
		sprintf(buf, "nu = %.2f", nus[v]);
		svg_panel(&s, &vp, buf);
		c = vs_curve(p);
		r = qaws_curve_get_parameter_range(c);
		fprintf(s.f, "<path fill=\"%s\" fill-opacity=\"0.18\" stroke=\"%s\" stroke-width=\"2.2\" d=\"", colors[v], colors[v]);
		for (i = 0; i <= 400; i++)
		{
			qaws_eval_result_2d e;
			qaws_curve_evaluate_2d(c, (qaws_scalar)(r.min_value + (r.max_value - r.min_value) * i / 400.0), QAWS_EVAL_FLAG_POSITION, &e);
			xy[2 * i] = vx(&vp, e.position.x);
			xy[2 * i + 1] = vy(&vp, e.position.y);
			fprintf(s.f, "%s%.2f %.2f", i ? "L" : "M", xy[2 * i], xy[2 * i + 1]);
		}
		fprintf(s.f, "Z\"/>\n");
		for (i = 0; i < VS_N; i++)
			svg_circle(&s, vx(&vp, p[2 * i]), vy(&vp, p[2 * i + 1]), 2.0, "#ffffff", colors[v]);
		qaws_curve_destroy(c);
		sprintf(buf, "energy %.4f, |L - 2pi| %.1e, |A - A0| %.1e", K, fabs(cl), fabs(ca));
		svg_text(&s, vp.x0 + 10, vp.y0 + vp.h - 26, 11, "#57606a", "start", buf);
		sprintf(buf, "%d Newton steps, %d HVP-CG products", iters[v], cg);
		svg_text(&s, vp.x0 + 10, vp.y0 + vp.h - 10, 11, "#57606a", "start", buf);
		printf("12_vesicles: nu %.2f: energy %.6f (circle 2 pi = 6.283185), |cL| %.1e, |cA| %.1e, %d Newton steps, %d CG products\n",
			nus[v], K, fabs(cl), fabs(ca), iters[v], cg);
	}
	{
		/* convergence: |gradient| per Newton step, log scale */
		viewport lp = { 1220, 70, 260, 330, 0, 1, 0, 1 };
		double lo = -12, hi = 3;
		int i;
		svg_panel(&s, &lp, "|gradient| per Newton step (log10)");
		for (v = 0; v < 4; v++)
		{
			double xy[2 * 256];
			for (i = 0; i < iters[v]; i++)
			{
				double l = log10(gn[v][i] > 1e-12 ? gn[v][i] : 1e-12);
				xy[2 * i] = lp.x0 + 14 + (lp.w - 28) * i / (double)(iters[v] > 1 ? iters[v] - 1 : 1);
				xy[2 * i + 1] = lp.y0 + 34 + (lp.h - 60) * (hi - l) / (hi - lo);
			}
			svg_polyline(&s, xy, iters[v], colors[v], 1.8, 0.9, 0);
		}
		svg_text(&s, lp.x0 + 12, lp.y0 + lp.h - 10, 11, "#57606a", "start", "1e3 .. 1e-12; jumps: lambda updates");
	}
	svg_close(&s);
}
