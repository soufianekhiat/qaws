/* Figure 8: every batch on the caller's threads. One small executor (a few
   threads taking chunks from a shared counter) handed to each batch through
   its desc; the speedup over the serial run for 1 to 16 threads, and the
   results compared bit for bit with the serial ones. */

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#include <process.h>
#else
#include <pthread.h>
#include <unistd.h>
#endif

#define TH_MAX 16
#define TH_OPS 9
#define TH_POINTS 20000
#define TH_RAYS (240 * 150)
#define TH_EXACT 128

/* the executor: `threads` workers (the caller is one) take the chunks in
   turn from a shared counter */
typedef struct th_run
{
	qaws_batch_task_fn task;
	void* ctx;
	unsigned int count;
	volatile long next;
} th_run;

static long th_take(th_run* r)
{
#ifdef _WIN32
	return InterlockedIncrement(&r->next) - 1;
#else
	return __atomic_fetch_add(&r->next, 1, __ATOMIC_RELAXED);
#endif
}

static void th_work(th_run* r)
{
	long k;
	while ((k = th_take(r)) < (long)r->count)
		r->task(r->ctx, (unsigned int)k, (unsigned int)k + 1);
}

#ifdef _WIN32
static unsigned __stdcall th_main(void* p) { th_work((th_run*)p); return 0; }
#else
static void* th_main(void* p) { th_work((th_run*)p); return NULL; }
#endif

static void th_parallel_for(void* user, unsigned int count, qaws_batch_task_fn task, void* ctx)
{
	unsigned int threads = *(unsigned int const*)user, k;
	th_run r;
#ifdef _WIN32
	HANDLE h[TH_MAX];
#else
	pthread_t h[TH_MAX];
#endif
	r.task = task;
	r.ctx = ctx;
	r.count = count;
	r.next = 0;
	if (threads > count)
		threads = count;
	for (k = 1; k < threads; k++)
#ifdef _WIN32
		h[k] = (HANDLE)_beginthreadex(NULL, 0, th_main, &r, 0, NULL);
#else
		pthread_create(&h[k], NULL, th_main, &r);
#endif
	th_work(&r);
	for (k = 1; k < threads; k++)
	{
#ifdef _WIN32
		WaitForSingleObject(h[k], INFINITE);
		CloseHandle(h[k]);
#else
		pthread_join(h[k], NULL);
#endif
	}
}

static unsigned int th_cores(void)
{
#ifdef _WIN32
	SYSTEM_INFO si;
	GetSystemInfo(&si);
	return (unsigned int)si.dwNumberOfProcessors;
#else
	long n = sysconf(_SC_NPROCESSORS_ONLN);
	return n > 0 ? (unsigned int)n : 1;
#endif
}

/* wall clock (clock() is processor time on some systems) */
static double th_now(void)
{
	struct timespec ts;
	timespec_get(&ts, TIME_UTC);
	return (double)ts.tv_sec + 1e-9 * ts.tv_nsec;
}

/* the workloads, all built once */
typedef struct th_scene
{
	hf_scene hf;
	qaws_curve const** contours;
	qaws_curve const** gradients;
	unsigned int ncontour, ngradient;
	qaws_scalar levels[96];
	qaws_surface* sf[1 + TR_LEVELS];
	qaws_exact_surface* es[1 + TR_LEVELS];
	unsigned int fam[1 + TR_LEVELS];
	qaws_scalar* pts2;              /* TH_POINTS points in the plane */
	qaws_scalar* pts3;              /* TH_POINTS points over the terrain */
	qaws_scalar* org;               /* TH_RAYS camera rays */
	qaws_scalar* dir;
	double epts[3 * TH_EXACT];           /* certified queries */
	double ep0[3 * TH_EXACT], ep1[3 * TH_EXACT];
} th_scene;


/* one run of operation op with the executor x: a checksum of the results
   (equal checksums for equal results) */
static unsigned long long th_hash(unsigned long long h, void const* p, size_t n)
{
	unsigned char const* b = (unsigned char const*)p;
	size_t i;
	for (i = 0; i < n; i++)
		h = (h ^ b[i]) * 1099511628211ull;
	return h;
}

static unsigned long long th_op(th_scene* sc, unsigned int op, qaws_batch_executor const* x)
{
	static qaws_curve_batch_hit_2d ch[1 << 16];
	static qaws_level_crossing lc[1 << 16];
	static qaws_closest_point cp[TH_POINTS];
	static qaws_surface_batch_closest sp[TH_POINTS];
	static qaws_surface_ray_hit rh[TH_RAYS];
	static qaws_surface_batch_curve cv[1024];
	static qaws_ssi_point pt[1 << 17];
	static qaws_exact_ssi_point ept[1 << 16];
	static qaws_exact_ssi_batch_branch ebr[1024];
	static qaws_exact_surface_closest_point ecp[TH_EXACT];
	static qaws_exact_ray_hit erh[TH_EXACT];
	unsigned long long h = 14695981039346656037ull;
	unsigned int n = 0, m = 0, i;
	switch (op)
	{
	case 0:
	{
		qaws_curve_batch_desc d;
		memset(&d, 0, sizeof(d));
		d.curves = (qaws_curve const* const*)sc->hf.curves;
		d.curve_count = sc->hf.count;
		d.families = sc->hf.family;
		d.executor = x;
		qaws_curve_batch_find_intersections_2d(&d, ch, 1 << 16, &n, NULL);
		return th_hash(h, ch, n * sizeof(ch[0]));
	}
	case 1:
	{
		qaws_level_crossing_desc d;
		memset(&d, 0, sizeof(d));
		d.curves = sc->gradients;
		d.curve_count = sc->ngradient;
		d.field = hf_field;
		d.levels = sc->levels;
		d.level_count = 96;
		d.executor = x;
		qaws_curve_batch_find_level_crossings(&d, lc, 1 << 16, &n);
		for (i = 0; i < n; i++)
		{
			h = th_hash(h, &lc[i].curve, 2 * sizeof(unsigned int));
			h = th_hash(h, &lc[i].parameter, sizeof(qaws_scalar));
		}
		return h;
	}
	case 2:
	{
		qaws_closest_desc d;
		memset(&d, 0, sizeof(d));
		d.curves = sc->contours;
		d.curve_count = sc->ncontour;
		d.points = sc->pts2;
		d.point_count = TH_POINTS;
		d.executor = x;
		qaws_curve_batch_find_closest(&d, cp, NULL);
		for (i = 0; i < TH_POINTS; i++)
		{
			h = th_hash(h, &cp[i].curve, sizeof(unsigned int));
			h = th_hash(h, &cp[i].parameter, 2 * sizeof(qaws_scalar));
		}
		return h;
	}
	case 3:
	{
		qaws_surface_batch_closest_desc d;
		memset(&d, 0, sizeof(d));
		d.surfaces = (qaws_surface const* const*)sc->sf;
		d.surface_count = 1;
		d.points = sc->pts3;
		d.point_count = TH_POINTS;
		d.executor = x;
		qaws_surface_batch_find_closest(&d, sp, NULL);
		for (i = 0; i < TH_POINTS; i++)
		{
			h = th_hash(h, &sp[i].surface, sizeof(unsigned int));
			h = th_hash(h, &sp[i].distance, sizeof(qaws_scalar));
		}
		return h;
	}
	case 4:
	{
		qaws_surface_ray_desc d;
		memset(&d, 0, sizeof(d));
		d.surfaces = (qaws_surface const* const*)sc->sf;
		d.surface_count = 1;
		d.origins = sc->org;
		d.directions = sc->dir;
		d.ray_count = TH_RAYS;
		d.executor = x;
		qaws_surface_batch_raycast(&d, rh, NULL);
		for (i = 0; i < TH_RAYS; i++)
		{
			h = th_hash(h, &rh[i].surface, sizeof(unsigned int));
			h = th_hash(h, &rh[i].t, sizeof(qaws_scalar));
		}
		return h;
	}
	case 5:
	{
		qaws_surface_batch_desc d;
		memset(&d, 0, sizeof(d));
		d.surfaces = (qaws_surface const* const*)sc->sf;
		d.surface_count = 1 + TR_LEVELS;
		d.families = sc->fam;
		d.executor = x;
		qaws_surface_batch_find_intersections(&d, cv, 1024, &n, pt, 1 << 17, &m, NULL);
		return th_hash(th_hash(h, cv, n * sizeof(cv[0])), pt, m * sizeof(pt[0]));
	}
	case 6:
	{
		qaws_exact_ssi_batch_desc d;
		memset(&d, 0, sizeof(d));
		d.surfaces = (qaws_exact_surface const* const*)sc->es;
		d.surface_count = 1 + TR_LEVELS;
		d.families = sc->fam;
		d.executor = x;
		qaws_exact_surface_batch_hits(&d, ept, 1 << 16, &m, ebr, 1024, &n, NULL);
		return th_hash(th_hash(h, ebr, n * sizeof(ebr[0])), ept, m * sizeof(ept[0]));
	}
	case 7:
	{
		qaws_exact_ssi_batch_desc d;
		memset(&d, 0, sizeof(d));
		d.surfaces = (qaws_exact_surface const* const*)sc->es;
		d.surface_count = 1;
		d.executor = x;
		qaws_exact_surface_batch_closest(&d, sc->epts, TH_EXACT, ecp, NULL);
		return th_hash(h, ecp, sizeof(ecp));
	}
	default:
	{
		qaws_exact_ssi_batch_desc d;
		memset(&d, 0, sizeof(d));
		d.surfaces = (qaws_exact_surface const* const*)sc->es;
		d.surface_count = 1;
		d.executor = x;
		qaws_exact_surface_batch_raycast(&d, sc->ep0, sc->ep1, TH_EXACT, 0, erh, NULL);
		for (i = 0; i < TH_EXACT; i++)
		{
			h = th_hash(h, &erh[i].surface, sizeof(unsigned int));
			h = th_hash(h, &erh[i].hit.t_lo, 2 * sizeof(double));
			h = th_hash(h, &erh[i].certified, sizeof(int));
		}
		return h;
	}
	}
}

static void th_build(th_scene* sc)
{
	static double f[(HF_GRID + 1) * (HF_GRID + 1)];
	unsigned long long state = 0xBB67AE8584CAA73Bull;
	double zlo = 1e30, zhi = -1e30, eye[3] = { 2.3, -2.7, 2.1 };
	qaws_exact_desc xd;
	unsigned int i, j;
	for (j = 0; j <= HF_GRID; j++)
		for (i = 0; i <= HF_GRID; i++)
			f[j * (HF_GRID + 1) + i] = hf_height(-1 + 2.0 * i / HF_GRID, -1 + 2.0 * j / HF_GRID, NULL, NULL);
	hf_build(&sc->hf, 96, f);
	sc->contours = (qaws_curve const**)malloc(sc->hf.count * sizeof(qaws_curve*));
	sc->gradients = (qaws_curve const**)malloc(sc->hf.count * sizeof(qaws_curve*));
	sc->ncontour = sc->ngradient = 0;
	for (i = 0; i < sc->hf.count; i++)
	{
		if (sc->hf.family[i] == 0)
			sc->contours[sc->ncontour++] = sc->hf.curves[i];
		else
			sc->gradients[sc->ngradient++] = sc->hf.curves[i];
	}
	for (i = 0; i < 96; i++)
		sc->levels[i] = (qaws_scalar)(-0.65 + 1.55 * (i + 0.5) / 96);
	sc->sf[0] = tr_terrain();
	sc->fam[0] = 0;
	for (i = 0; i <= 40; i++)
		for (j = 0; j <= 40; j++)
		{
			double z = tr_value(sc->sf[0], (qaws_scalar)(i / 40.0), (qaws_scalar)(j / 40.0));
			if (z < zlo) zlo = z;
			if (z > zhi) zhi = z;
		}
	for (i = 0; i < TR_LEVELS; i++)
	{
		sc->sf[1 + i] = tr_plane(ldexp(nearbyint(ldexp(zlo + (zhi - zlo) * (i + 0.5) / TR_LEVELS, 12)), -12));
		sc->fam[1 + i] = 1;
	}
	qaws_exact_desc_default(&xd);
	xd.space_exp2 = -20;
	for (i = 0; i <= TR_LEVELS; i++)
		qaws_exact_surface_prepare(&xd, sc->sf[i], &sc->es[i], NULL);
	sc->pts2 = (qaws_scalar*)malloc(2 * TH_POINTS * sizeof(qaws_scalar));
	sc->pts3 = (qaws_scalar*)malloc(3 * TH_POINTS * sizeof(qaws_scalar));
	sc->org = (qaws_scalar*)malloc(3 * TH_RAYS * sizeof(qaws_scalar));
	sc->dir = (qaws_scalar*)malloc(3 * TH_RAYS * sizeof(qaws_scalar));
	for (i = 0; i < TH_POINTS; i++)
	{
		double r[3];
		unsigned int k;
		for (k = 0; k < 3; k++)
		{
			state ^= state << 13; state ^= state >> 7; state ^= state << 17;
			r[k] = (double)(state >> 11) / 9007199254740992.0;
		}
		sc->pts2[2 * i] = (qaws_scalar)(-0.98 + 1.96 * r[0]);
		sc->pts2[2 * i + 1] = (qaws_scalar)(-0.98 + 1.96 * r[1]);
		sc->pts3[3 * i] = sc->pts2[2 * i];
		sc->pts3[3 * i + 1] = sc->pts2[2 * i + 1];
		sc->pts3[3 * i + 2] = (qaws_scalar)(-0.8 + 1.6 * r[2]);
		if (i < TH_EXACT)
			for (k = 0; k < 3; k++)
			{
				/* dyadic: exact as given */
				sc->epts[3 * i + k] = ldexp(nearbyint(ldexp((double)(k < 2 ? sc->pts2[2 * i + k] : sc->pts3[3 * i + 2]), 12)), -12);
				sc->ep0[3 * i + k] = eye[k];
				sc->ep1[3 * i + k] = ldexp(nearbyint(ldexp(k < 2 ? 0.9 * (2 * r[k] - 1) : 0.0, 12)), -12);
			}
	}
	for (j = 0; j < 150; j++)
		for (i = 0; i < 240; i++)
		{
			unsigned int p = j * 240 + i, k;
			double t[3];
			t[0] = -1 + 2.0 * (i + 0.5) / 240;
			t[1] = -1 + 2.0 * (j + 0.5) / 150;
			t[2] = -0.1;
			for (k = 0; k < 3; k++)
			{
				sc->org[3 * p + k] = (qaws_scalar)eye[k];
				sc->dir[3 * p + k] = (qaws_scalar)(t[k] - eye[k]);
			}
		}
}

static void demo_threads(void)
{
	static char const* const name[TH_OPS] = { "curve intersections", "level crossings", "closest points on curves", "closest points on a surface",
		"ray casting", "surface intersections", "certified surface intersections", "certified closest points", "certified rays" };
	static char const* const what[TH_OPS] = { "contours x gradient lines", "gradient lines x 96 levels", "20000 points, contours", "20000 points, terrain",
		"36000 rays, terrain", "terrain x 16 planes", "terrain x 16 planes", "128 points, terrain", "128 rays, terrain" };
	static unsigned int const counts[5] = { 1, 2, 4, 8, 16 };
	th_scene* sc = (th_scene*)calloc(1, sizeof(th_scene));
	double t[TH_OPS][5];
	int same[TH_OPS];
	unsigned int op, c, i;
	svg s;
	char sub[400], buf[200];
	th_build(sc);
	for (op = 0; op < TH_OPS; op++)
	{
		unsigned long long ref = th_op(sc, op, NULL);
		same[op] = 1;
		for (c = 0; c < 5; c++)
		{
			qaws_batch_executor ex;
			unsigned int threads = counts[c], reps = 0;
			double t0, best = 1e30;
			ex.parallel_for = th_parallel_for;
			ex.user = &threads;
			/* the best of the runs over 0.3 s (at least one) */
			t0 = th_now();
			do
			{
				double t1 = th_now(), dt;
				same[op] &= th_op(sc, op, &ex) == ref;
				dt = th_now() - t1;
				if (dt < best)
					best = dt;
				reps++;
			} while (th_now() - t0 < 0.3 && reps < 50);
			t[op][c] = best;
		}
		printf("    %-32s 1 thread %8.1f ms, x%.1f on 2, x%.1f on 4, x%.1f on 8, x%.1f on 16; same results %d\n", name[op], t[op][0] * 1000,
			t[op][0] / t[op][1], t[op][0] / t[op][2], t[op][0] / t[op][3], t[op][0] / t[op][4], same[op]);
	}

	sprintf(sub, "one small executor handed to every batch through its desc: chunks of points, rays, grid cells or surface pairs on 1 to 16 threads (%u cores); results identical to the serial run, bit for bit",
		th_cores());
	if (!svg_open(&s, "showcase/batch8_threads.svg", 1400, 760, "Batches on the caller's threads", sub))
		return;
	{
		static char const* const tc[5] = { "#9ecae1", "#6baed6", "#3182bd", "#08519c", "#08306b" };
		double x0 = 420, x1 = 1330, y0 = 100, row = 66, smax = 16;
		/* the speedup axis */
		for (i = 0; i <= 16; i += 2)
		{
			double x = x0 + (x1 - x0) * i / smax;
			svg_line(&s, x, y0 - 6, x, y0 + row * TH_OPS, "#d0d7de", 1, 1);
			sprintf(buf, "x%u", i);
			svg_text(&s, x, y0 + row * TH_OPS + 18, 12, "#57606a", "middle", buf);
		}
		svg_text(&s, (x0 + x1) / 2, y0 + row * TH_OPS + 38, 13, "#24292f", "middle", "speedup over one thread");
		for (op = 0; op < TH_OPS; op++)
		{
			double y = y0 + row * op;
			svg_text(&s, 30, y + 22, 15, "#1b1f24", "start", name[op]);
			sprintf(buf, "%s; serial %.1f ms%s", what[op], t[op][0] * 1000, same[op] ? "" : " (RESULTS DIFFER)");
			svg_text(&s, 30, y + 42, 12, same[op] ? "#57606a" : "#cf222e", "start", buf);
			for (c = 0; c < 5; c++)
			{
				double sp = t[op][0] / t[op][c], w = (x1 - x0) * (sp > smax ? smax : sp) / smax;
				fprintf(s.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"10\" fill=\"%s\"/>\n", x0, y + 4 + 11 * c, w, tc[c]);
				sprintf(buf, "x%.1f", sp);
				svg_text(&s, x0 + w + 6, y + 13 + 11 * c, 10, "#24292f", "start", buf);
			}
		}
		for (c = 0; c < 5; c++)
		{
			fprintf(s.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"14\" height=\"10\" fill=\"%s\"/>\n", 30.0 + 100 * c, 735.0, tc[c]);
			sprintf(buf, "%u thread%s", counts[c], c ? "s" : "");
			svg_text(&s, 50 + 100 * c, 744, 11, "#24292f", "start", buf);
		}
	}
	svg_close(&s);
	printf("  -> showcase/batch8_threads.svg\n");
	for (i = 0; i <= TR_LEVELS; i++)
	{
		qaws_exact_surface_destroy(sc->es[i]);
		qaws_surface_destroy(sc->sf[i]);
	}
	hf_free(&sc->hf);
	free((void*)sc->contours);
	free((void*)sc->gradients);
	free(sc->pts2);
	free(sc->pts3);
	free(sc->org);
	free(sc->dir);
	free(sc);
}
