/* ================================================================== */
/*  4. Non-rigid registration: CMA-ES basin hopping + gradients       */
/* ================================================================== */

#define RG_N 16                /* template interpolation points */
#define RG_SCAN 120            /* scan points */
#define RG_LAMBDA 10           /* CMA-ES population */
#define RG_GEN 18
#define RG_LOCAL 10            /* gradient steps inside the fitness */
#define RG_FINAL 250
#define RG_TRIALS 16

typedef struct reg_problem
{
	qaws_scalar tmpl[RG_N * 2];        /* template points */
	double scan[RG_SCAN * 2];
	int nscan;
	double cx, cy;                     /* scan centroid */
} reg_problem;

/* pose: theta, log scale, tx, ty; plus per-point offsets (non-rigid) */
typedef struct reg_state
{
	double pose[4];
	double off[RG_N * 2];
	double foot[RG_SCAN];      /* foot parameters of the scan points */
} reg_state;

static double rg_template_radius(double phi)
{
	return 1.0 + 0.30 * cos(2 * phi + 0.4) + 0.20 * sin(3 * phi) + 0.12 * cos(5 * phi + 1.0);
}

static void rg_points(reg_problem const* pb, reg_state const* st, qaws_scalar* out)
{
	double c = cos(st->pose[0]), s = sin(st->pose[0]), sc = exp(st->pose[1]);
	int i;
	for (i = 0; i < RG_N; i++)
	{
		double x = pb->tmpl[2 * i] + st->off[2 * i], y = pb->tmpl[2 * i + 1] + st->off[2 * i + 1];
		out[2 * i] = (qaws_scalar)(sc * (c * x - s * y) + st->pose[2]);
		out[2 * i + 1] = (qaws_scalar)(sc * (s * x + c * y) + st->pose[3]);
	}
}

static qaws_curve* rg_curve(qaws_scalar const* pts)
{
	qaws_catmull_rom_desc d;
	qaws_curve* c = NULL;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D;
	d.control_points = pts;
	d.control_point_count = RG_N;
	d.parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL;
	d.closed = 1;
	qaws_curve_create_catmull_rom(&d, &c);
	return c;
}

/* Foot points: nearest of 8 samples per span (global), refined later by
   Newton. Called when the pose jumps (CMA-ES candidates). */
static void rg_seed(reg_problem const* pb, reg_state* st)
{
	enum { M = RG_N * 8 };
	qaws_scalar pts[RG_N * 2], ts[M];
	qaws_curve_jet_2d prim[M], tan[M];
	qaws_curve* c;
	int k, j;
	rg_points(pb, st, pts);
	c = rg_curve(pts);
	for (j = 0; j < M; j++)
		ts[j] = (qaws_scalar)(RG_N * (j + 0.5) / M);
	qaws_curve_eval_batch_tangent_2d(NULL, c, ts, NULL, M, QAWS_EVAL_FLAG_POSITION, NULL, prim, tan);
	for (k = 0; k < pb->nscan; k++)
	{
		double best = 1e300;
		for (j = 0; j < M; j++)
		{
			double dx = prim[j].d[0].x - pb->scan[2 * k], dy = prim[j].d[0].y - pb->scan[2 * k + 1];
			if (dx * dx + dy * dy < best)
			{
				best = dx * dx + dy * dy;
				st->foot[k] = ts[j];
			}
		}
	}
	qaws_curve_destroy(c);
}

/* E = mean squared distance of the scan to the deformed template at the
   foot points + smoothness of the offsets. The foot points take two Newton
   steps first; at a foot point the distance gradient is the pullback of
   2 (C - q) alone, through the Catmull-Rom adjoint onto the points, then
   onto pose and offsets. */
static double rg_energy(reg_problem const* pb, reg_state* st, reg_state* g, int nonrigid)
{
	qaws_scalar pts[RG_N * 2], pbar[RG_N * 2], ts[RG_SCAN];
	qaws_curve_jet_2d prim[RG_SCAN], tan[RG_SCAN], bar[RG_SCAN];
	qaws_curve* c;
	qaws_field_view fv;
	qaws_diff_views views = one_field(&fv, QAWS_FIELD_POINTS, pbar, RG_N, 2);
	double e = 0, lam = 0.4;
	int k, i, it;
	rg_points(pb, st, pts);
	c = rg_curve(pts);
	for (it = 0; it < 3; it++)
	{
		for (k = 0; k < pb->nscan; k++)
			ts[k] = (qaws_scalar)st->foot[k];
		qaws_curve_eval_batch_tangent_2d(NULL, c, ts, NULL, (unsigned int)pb->nscan,
			QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2, NULL, prim, tan);
		if (it == 2)
			break;
		for (k = 0; k < pb->nscan; k++)
		{
			double rx = prim[k].d[0].x - pb->scan[2 * k], ry = prim[k].d[0].y - pb->scan[2 * k + 1];
			double gv = rx * prim[k].d[1].x + ry * prim[k].d[1].y;
			double h = prim[k].d[1].x * prim[k].d[1].x + prim[k].d[1].y * prim[k].d[1].y + rx * prim[k].d[2].x + ry * prim[k].d[2].y;
			double t = st->foot[k] - (h > 1e-9 ? gv / h : 0);
			if (t - st->foot[k] > 0.5) t = st->foot[k] + 0.5;
			if (t - st->foot[k] < -0.5) t = st->foot[k] - 0.5;
			t = fmod(t, (double)RG_N);
			if (t < 0) t += RG_N;
			st->foot[k] = t;
		}
	}
	for (k = 0; k < pb->nscan; k++)
	{
		double rx = prim[k].d[0].x - pb->scan[2 * k], ry = prim[k].d[0].y - pb->scan[2 * k + 1];
		e += (rx * rx + ry * ry) / pb->nscan;
		memset(&bar[k], 0, sizeof(bar[k]));
		bar[k].d[0].x = (qaws_scalar)(2 * rx / pb->nscan);
		bar[k].d[0].y = (qaws_scalar)(2 * ry / pb->nscan);
		bar[k].channels = QAWS_EVAL_FLAG_POSITION;
	}
	for (i = 0; i < RG_N; i++)
	{
		int j = (i + 1) % RG_N;
		double dx = st->off[2 * i] - st->off[2 * j], dy = st->off[2 * i + 1] - st->off[2 * j + 1];
		e += lam * (dx * dx + dy * dy) / RG_N;
	}
	if (g)
	{
		double cth = cos(st->pose[0]), sth = sin(st->pose[0]), sc = exp(st->pose[1]);
		memset(pbar, 0, sizeof(pbar));
		qaws_curve_eval_batch_adjoint_2d(NULL, c, ts, (unsigned int)pb->nscan, QAWS_EVAL_FLAG_POSITION, bar, &views, NULL);
		memset(g->pose, 0, sizeof(g->pose));
		memset(g->off, 0, sizeof(g->off));
		for (i = 0; i < RG_N; i++)
		{
			double x = pb->tmpl[2 * i] + st->off[2 * i], y = pb->tmpl[2 * i + 1] + st->off[2 * i + 1];
			double gx = pbar[2 * i], gy = pbar[2 * i + 1];
			/* P' = s R p + t */
			g->pose[0] += sc * (gx * (-sth * x - cth * y) + gy * (cth * x - sth * y));
			g->pose[1] += sc * (gx * (cth * x - sth * y) + gy * (sth * x + cth * y));
			g->pose[2] += gx;
			g->pose[3] += gy;
			if (nonrigid)
			{
				int j = (i + 1) % RG_N, h = (i + RG_N - 1) % RG_N;
				g->off[2 * i] += sc * (cth * gx + sth * gy);
				g->off[2 * i + 1] += sc * (-sth * gx + cth * gy);
				g->off[2 * i] += lam * 2 * (2 * st->off[2 * i] - st->off[2 * j] - st->off[2 * h]) / RG_N;
				g->off[2 * i + 1] += lam * 2 * (2 * st->off[2 * i + 1] - st->off[2 * j + 1] - st->off[2 * h + 1]) / RG_N;
			}
		}
	}
	qaws_curve_destroy(c);
	return e;
}

/* Local refinement with Adam; returns the final energy. */
static double rg_refine(reg_problem const* pb, reg_state* st, int steps, int nonrigid, double lr)
{
	adam a;
	int it, k;
	double e = 0;
	memset(&a, 0, sizeof(a));
	for (it = 0; it <= steps; it++)
	{
		reg_state g;
		qaws_scalar x[4 + RG_N * 2], gg[4 + RG_N * 2];
		e = rg_energy(pb, st, it < steps ? &g : NULL, nonrigid);
		if (it == steps)
			break;
		for (k = 0; k < 4; k++) { x[k] = (qaws_scalar)st->pose[k]; gg[k] = (qaws_scalar)g.pose[k]; }
		for (k = 0; k < RG_N * 2; k++) { x[4 + k] = (qaws_scalar)st->off[k]; gg[4 + k] = (qaws_scalar)g.off[k]; }
		adam_step(&a, x, gg, nonrigid ? 4 + RG_N * 2 : 4, lr * (1.0 - 0.8 * it / (double)steps));
		for (k = 0; k < 4; k++) st->pose[k] = x[k];
		if (nonrigid)
			for (k = 0; k < RG_N * 2; k++) st->off[k] = x[4 + k];
	}
	return e;
}

/* --- CMA-ES (4D pose, rank-one and rank-mu updates) ---------------- */

static double rg_gauss(unsigned int* rng)
{
	double u1, u2;
	*rng = *rng * 1664525u + 1013904223u;
	u1 = ((*rng >> 8) + 0.5) / 16777216.0;
	*rng = *rng * 1664525u + 1013904223u;
	u2 = ((*rng >> 8) + 0.5) / 16777216.0;
	return sqrt(-2 * log(u1)) * cos(2 * PI * u2);
}

/* Jacobi eigen decomposition of a symmetric 4x4: C = B diag(d) B^T. */
static void rg_eigen(double C[4][4], double B[4][4], double* d)
{
	double A[4][4];
	int i, j, k, sweep;
	memcpy(A, C, sizeof(A));
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			B[i][j] = i == j;
	for (sweep = 0; sweep < 30; sweep++)
		for (i = 0; i < 4; i++)
			for (j = i + 1; j < 4; j++)
			{
				double th, c, s, t;
				if (fabs(A[i][j]) < 1e-15)
					continue;
				th = (A[j][j] - A[i][i]) / (2 * A[i][j]);
				t = (th >= 0 ? 1 : -1) / (fabs(th) + sqrt(th * th + 1));
				c = 1 / sqrt(t * t + 1);
				s = t * c;
				for (k = 0; k < 4; k++)
				{
					double aki = A[k][i], akj = A[k][j];
					A[k][i] = c * aki - s * akj;
					A[k][j] = s * aki + c * akj;
				}
				for (k = 0; k < 4; k++)
				{
					double aik = A[i][k], ajk = A[j][k];
					A[i][k] = c * aik - s * ajk;
					A[j][k] = s * aik + c * ajk;
				}
				for (k = 0; k < 4; k++)
				{
					double bki = B[k][i], bkj = B[k][j];
					B[k][i] = c * bki - s * bkj;
					B[k][j] = s * bki + c * bkj;
				}
			}
	for (i = 0; i < 4; i++)
		d[i] = A[i][i] > 1e-20 ? A[i][i] : 1e-20;
}

typedef struct cma_trace
{
	double best[RG_GEN];
	double cand[RG_GEN][RG_LAMBDA][4];
} cma_trace;

/* Minimizes the refined energy over the pose; scaled coordinates
   y = (pose - pose0) / scale. Returns the best pose. */
static void rg_cmaes(reg_problem const* pb, double const* pose0, double const* scale, unsigned int seed, double* best_pose,
	cma_trace* trace)
{
	enum { N = 4, MU = RG_LAMBDA / 2 };
	double m[N] = { 0, 0, 0, 0 }, sigma = 0.6, C[N][N], B[N][N], D[N], pc[N] = { 0 }, ps[N] = { 0 };
	double w[MU], wsum = 0, mueff, cc, cs, c1, cmu, damps, chi = sqrt((double)N) * (1 - 1.0 / (4 * N) + 1.0 / (21 * N * N));
	double best_f = 1e300;
	unsigned int rng = seed;
	int i, j, k, gen;
	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
			C[i][j] = i == j;
	for (i = 0; i < MU; i++)
	{
		w[i] = log(MU + 0.5) - log(i + 1.0);
		wsum += w[i];
	}
	mueff = 0;
	for (i = 0; i < MU; i++)
	{
		w[i] /= wsum;
		mueff += w[i] * w[i];
	}
	mueff = 1 / mueff;
	cc = (4 + mueff / N) / (N + 4 + 2 * mueff / N);
	cs = (mueff + 2) / (N + mueff + 5);
	c1 = 2 / ((N + 1.3) * (N + 1.3) + mueff);
	cmu = 2 * (mueff - 2 + 1 / mueff) / ((N + 2) * (N + 2) + mueff);
	if (cmu > 1 - c1) cmu = 1 - c1;
	damps = 1 + 2 * (mueff > N + 1 ? sqrt((mueff - 1) / (N + 1)) - 1 : 0) + cs;
	memcpy(best_pose, pose0, sizeof(double) * N);

	for (gen = 0; gen < RG_GEN; gen++)
	{
		double y[RG_LAMBDA][N], z[RG_LAMBDA][N], f[RG_LAMBDA], mold[N];
		int order[RG_LAMBDA];
		rg_eigen(C, B, D);
		for (k = 0; k < RG_LAMBDA; k++)
		{
			reg_state st;
			for (i = 0; i < N; i++)
				z[k][i] = rg_gauss(&rng);
			for (i = 0; i < N; i++)
			{
				double v = 0;
				for (j = 0; j < N; j++)
					v += B[i][j] * sqrt(D[j]) * z[k][j];
				y[k][i] = m[i] + sigma * v;
			}
			memset(&st, 0, sizeof(st));
			for (i = 0; i < N; i++)
				st.pose[i] = pose0[i] + scale[i] * y[k][i];
			rg_seed(pb, &st);
			/* fitness: energy after a short rigid refinement (basin hopping) */
			f[k] = rg_refine(pb, &st, RG_LOCAL, 0, 0.05);
			for (i = 0; i < N; i++)
				trace->cand[gen][k][i] = st.pose[i];
			if (f[k] < best_f)
			{
				best_f = f[k];
				memcpy(best_pose, st.pose, sizeof(double) * N);
			}
			order[k] = k;
		}
		trace->best[gen] = best_f;
		for (i = 1; i < RG_LAMBDA; i++)
		{
			int key = order[i];
			j = i - 1;
			while (j >= 0 && f[order[j]] > f[key]) { order[j + 1] = order[j]; j--; }
			order[j + 1] = key;
		}
		memcpy(mold, m, sizeof(m));
		for (i = 0; i < N; i++)
		{
			m[i] = 0;
			for (k = 0; k < MU; k++)
				m[i] += w[k] * y[order[k]][i];
		}
		{
			/* evolution paths */
			double zm[N] = { 0 }, bz[N], hs, norm = 0;
			for (k = 0; k < MU; k++)
				for (i = 0; i < N; i++)
					zm[i] += w[k] * z[order[k]][i];
			for (i = 0; i < N; i++)
			{
				bz[i] = 0;
				for (j = 0; j < N; j++)
					bz[i] += B[i][j] * zm[j];
			}
			for (i = 0; i < N; i++)
			{
				ps[i] = (1 - cs) * ps[i] + sqrt(cs * (2 - cs) * mueff) * bz[i];
				norm += ps[i] * ps[i];
			}
			norm = sqrt(norm);
			hs = norm / sqrt(1 - pow(1 - cs, 2.0 * (gen + 1))) / chi < 1.4 + 2.0 / (N + 1) ? 1 : 0;
			for (i = 0; i < N; i++)
				pc[i] = (1 - cc) * pc[i] + hs * sqrt(cc * (2 - cc) * mueff) * (m[i] - mold[i]) / sigma;
			for (i = 0; i < N; i++)
				for (j = 0; j < N; j++)
				{
					double rmu = 0;
					for (k = 0; k < MU; k++)
						rmu += w[k] * (y[order[k]][i] - mold[i]) * (y[order[k]][j] - mold[j]) / (sigma * sigma);
					C[i][j] = (1 - c1 - cmu) * C[i][j] + c1 * (pc[i] * pc[j] + (1 - hs) * cc * (2 - cc) * C[i][j]) + cmu * rmu;
				}
			sigma *= exp((cs / damps) * (norm / chi - 1));
		}
	}
}

static void rg_make_problem(reg_problem* pb, double theta, unsigned int seed)
{
	unsigned int rng = seed;
	int i, k = 0;
	for (i = 0; i < RG_N; i++)
	{
		double phi = 2 * PI * i / RG_N, r = rg_template_radius(phi);
		pb->tmpl[2 * i] = (qaws_scalar)(r * cos(phi));
		pb->tmpl[2 * i + 1] = (qaws_scalar)(r * sin(phi));
	}
	/* scan: rotated, scaled, bent, noisy, with a missing arc */
	for (i = 0; i < RG_SCAN; i++)
	{
		double phi = 2 * PI * i / RG_SCAN, r = rg_template_radius(phi), x, y, bx, by;
		if (phi > 4.3 && phi < 5.3)
			continue;
		x = r * cos(phi);
		y = r * sin(phi);
		bx = x + 0.12 * sin(1.3 * y);
		by = y + 0.08 * cos(1.1 * x);
		pb->scan[2 * k] = 1.15 * (cos(theta) * bx - sin(theta) * by) + 3.0 + 0.012 * rg_gauss(&rng);
		pb->scan[2 * k + 1] = 1.15 * (sin(theta) * bx + cos(theta) * by) + 0.4 + 0.012 * rg_gauss(&rng);
		k++;
	}
	pb->nscan = k;
	pb->cx = pb->cy = 0;
	for (i = 0; i < k; i++)
	{
		pb->cx += pb->scan[2 * i] / k;
		pb->cy += pb->scan[2 * i + 1] / k;
	}
}

static void rg_curve_svg(svg* s, viewport const* v, reg_problem const* pb, reg_state const* st, char const* color, double width,
	double opacity, int dashed)
{
	qaws_scalar pts[RG_N * 2];
	qaws_curve* c;
	double xy[2 * 200];
	int k;
	rg_points(pb, st, pts);
	c = rg_curve(pts);
	for (k = 0; k < 200; k++)
	{
		qaws_eval_result_2d e;
		qaws_curve_evaluate_2d(c, (qaws_scalar)(RG_N * k / 199.0), QAWS_EVAL_FLAG_POSITION, &e);
		xy[2 * k] = vx(v, e.position.x);
		xy[2 * k + 1] = vy(v, e.position.y);
	}
	svg_polyline(s, xy, 200, color, width, opacity, dashed);
	qaws_curve_destroy(c);
}

static void app_registration(void)
{
	reg_problem pb;
	reg_state gd, cm, init;
	cma_trace trace;
	double scale[4] = { 2.0, 0.3, 0.6, 0.6 }, pose0[4], best_pose[4], e_gd, e_cm;
	int success_gd = 0, success_cm = 0, t;
	double trials_theta[RG_TRIALS], trials_gd[RG_TRIALS], trials_cm[RG_TRIALS];
	svg s;
	char buf[256];

	/* the featured case: a large rotation */
	rg_make_problem(&pb, 2.6, 7u);
	memset(&init, 0, sizeof(init));
	init.pose[2] = pb.cx;
	init.pose[3] = pb.cy;
	rg_seed(&pb, &init);
	gd = init;
	e_gd = rg_refine(&pb, &gd, RG_FINAL, 1, 0.03);
	memcpy(pose0, init.pose, sizeof(pose0));
	rg_cmaes(&pb, pose0, scale, 99u, best_pose, &trace);
	memset(&cm, 0, sizeof(cm));
	memcpy(cm.pose, best_pose, sizeof(best_pose));
	rg_seed(&pb, &cm);
	e_cm = rg_refine(&pb, &cm, RG_FINAL, 1, 0.03);

	/* statistics over random rotations */
	for (t = 0; t < RG_TRIALS; t++)
	{
		reg_problem q;
		reg_state a, b;
		double th = 2 * PI * (t + 0.5) / RG_TRIALS - PI, bp[4];
		cma_trace tr;
		rg_make_problem(&q, th, 100u + (unsigned int)t);
		memset(&a, 0, sizeof(a));
		a.pose[2] = q.cx;
		a.pose[3] = q.cy;
		rg_seed(&q, &a);
		b = a;
		trials_theta[t] = th;
		trials_gd[t] = sqrt(rg_refine(&q, &a, RG_FINAL, 1, 0.03));
		rg_cmaes(&q, b.pose, scale, 500u + (unsigned int)t, bp, &tr);
		memcpy(b.pose, bp, sizeof(bp));
		rg_seed(&q, &b);
		trials_cm[t] = sqrt(rg_refine(&q, &b, RG_FINAL, 1, 0.03));
		success_gd += trials_gd[t] < 0.03;
		success_cm += trials_cm[t] < 0.03;
	}

	svg_open(&s, "showcase/app4_registration.svg", 1240, 760, "Non-rigid registration: CMA-ES basin hopping + adjoints",
		"A template (closed centripetal Catmull-Rom) is registered to a partial, noisy, bent and rotated scan. Energy: "
		"squared distances at tracked foot points (exact batch adjoints) + offset smoothness.");
	{
		viewport a = { 20, 90, 390, 390, 0.6, 5.4, -2.0, 2.8 }, b = a;
		viewport lv = { 840, 90, 380, 220, 0, 0, 0, 0 }, sv = { 840, 330, 380, 230, 0, 0, 0, 0 };
		int k, gen;
		b.x0 = 425;
		svg_panel(&s, &a, "gradient descent from the identity pose");
		svg_panel(&s, &b, "CMA-ES on the pose, gradient inside the fitness");
		for (k = 0; k < pb.nscan; k++)
		{
			svg_circle(&s, vx(&a, pb.scan[2 * k]), vy(&a, pb.scan[2 * k + 1]), 2.2, "#24292f", "#24292f");
			svg_circle(&s, vx(&b, pb.scan[2 * k]), vy(&b, pb.scan[2 * k + 1]), 2.2, "#24292f", "#24292f");
		}
		rg_curve_svg(&s, &a, &pb, &init, "#8c959f", 1.6, 1, 1);
		rg_curve_svg(&s, &a, &pb, &gd, "#cf222e", 2.6, 1, 0);
		for (gen = 0; gen < RG_GEN; gen += 3)
			for (k = 0; k < RG_LAMBDA; k++)
			{
				reg_state st;
				memset(&st, 0, sizeof(st));
				memcpy(st.pose, trace.cand[gen][k], sizeof(st.pose));
				rg_curve_svg(&s, &b, &pb, &st, "#54aeff", 0.8, 0.35, 0);
			}
		rg_curve_svg(&s, &b, &pb, &cm, "#0969da", 2.8, 1, 0);
		sprintf(buf, "RMS %.4f (wrong basin)", sqrt(e_gd));
		svg_text(&s, a.x0 + a.w - 10, a.y0 + a.h - 10, 13, "#cf222e", "end", buf);
		sprintf(buf, "RMS %.4f", sqrt(e_cm));
		svg_text(&s, b.x0 + b.w - 10, b.y0 + b.h - 10, 13, "#0969da", "end", buf);
		svg_text(&s, a.x0 + 10, a.y0 + a.h - 10, 12, "#57606a", "start", "dots: scan, dashed: template");
		svg_text(&s, b.x0 + 10, b.y0 + a.h - 10, 12, "#57606a", "start", "light: candidates every 3 generations");

		{
			double lo = 1e300, hi = -1e300, xy[2 * RG_GEN];
			for (gen = 0; gen < RG_GEN; gen++)
			{
				double l = log10(sqrt(trace.best[gen]));
				if (l < lo) lo = l;
				if (l > hi) hi = l;
			}
			if (hi - lo < 0.1) hi = lo + 0.1;
			svg_panel(&s, &lv, "best refined RMS per CMA-ES generation (log)");
			for (gen = 0; gen < RG_GEN; gen++)
			{
				xy[2 * gen] = lv.x0 + 14 + (lv.w - 28) * gen / (double)(RG_GEN - 1);
				xy[2 * gen + 1] = lv.y0 + 34 + (lv.h - 54) * (hi - log10(sqrt(trace.best[gen]))) / (hi - lo);
			}
			svg_polyline(&s, xy, RG_GEN, "#0969da", 2.4, 1, 0);
		}
		svg_panel(&s, &sv, "final RMS over 16 scan rotations (-180 to 180 deg)");
		{
			double bw = (sv.w - 40) / RG_TRIALS;
			for (t = 0; t < RG_TRIALS; t++)
			{
				double x = sv.x0 + 20 + t * bw, base = sv.y0 + sv.h - 30, hs = sv.h - 70;
				double h0 = hs * (trials_gd[t] > 0.25 ? 1 : trials_gd[t] / 0.25), h1 = hs * (trials_cm[t] > 0.25 ? 1 : trials_cm[t] / 0.25);
				fprintf(s.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"#cf222e\" opacity=\"0.8\"/>\n", x + 1, base - h0, bw * 0.42, h0 + 0.5);
				fprintf(s.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"#0969da\"/>\n", x + bw * 0.46, base - h1, bw * 0.42, h1 + 0.5);
			}
			svg_line(&s, sv.x0 + 20, sv.y0 + sv.h - 30 - (sv.h - 70) * 0.03 / 0.25, sv.x0 + sv.w - 20, sv.y0 + sv.h - 30 - (sv.h - 70) * 0.03 / 0.25, "#57606a", 1, 0.8);
			svg_text(&s, sv.x0 + 20, sv.y0 + sv.h - 12, 11, "#57606a", "start", "red: gradient descent, blue: CMA-ES + gradient; line: success threshold 0.03");
		}
		{
			double tx = 840, ty = 612;
			sprintf(buf, "registered (RMS &lt; 0.03): gradient descent %d / %d", success_gd, RG_TRIALS);
			svg_text(&s, tx, ty - 22, 13, "#cf222e", "start", buf);
			sprintf(buf, "registered (RMS &lt; 0.03): CMA-ES + gradient %d / %d", success_cm, RG_TRIALS);
			svg_text(&s, tx, ty, 13, "#24292f", "start", buf);
			svg_text(&s, tx, ty + 24, 12, "#57606a", "start", "fitness: pose energy after 10 rigid Adam steps (Baldwinian)");
			svg_text(&s, tx, ty + 42, 12, "#57606a", "start", "final: pose + 32 offsets, 250 Adam steps");
			svg_text(&s, tx, ty + 60, 12, "#57606a", "start", "gradient: 2 (C - q) -> Catmull-Rom adjoint -> points -> pose, offsets");
		}
	}
	svg_close(&s);
	printf("4_registration: featured RMS gd %.4f cma %.4f; success gd %d/%d cma %d/%d\n", sqrt(e_gd), sqrt(e_cm),
		success_gd, RG_TRIALS, success_cm, RG_TRIALS);
	(void)trials_theta;
}

