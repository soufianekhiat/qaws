/*
 * repulsive_paper.c - figures of Repulsive Curves (Yu, Schumacher, Crane,
 * ACM TOG 2021) reproduced with qaws_repulsive, from the paper's own scenes
 * (examples/data/repulsive_curves, MIT license, see LICENSE there).
 *
 * Each scene is the paper's: its curve, its (alpha, beta), every edge's
 * length and the barycenter held, solved by the dense fractional Sobolev
 * descent on the caller's threads. Writes into ./showcase/:
 *   paper1_freedman.svg   Freedman's unknot untangled to the circle (Fig. 2)
 *   paper2_k5.svg         K5 embedded in space (Fig. 27)
 *   paper3_k33.svg        K3,3 embedded in space (Fig. 27)
 *   paper4_earbuds.svg    a tangled earbud cable untangled (Fig. 33)
 *   paper5_planar.svg     a tree drawn in the plane (Fig. 26)
 *
 * Usage: qaws_repulsive_paper [data directory] (default
 * ../examples/data/repulsive_curves, from working_dir). Minutes, not
 * seconds: the scenes have 500 to 800 vertices and take hundreds of steps
 * (about 20 minutes on 20 threads).
 */

#include "qaws.h"
#include "qaws_repulsive.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#include <process.h>
#include <direct.h>
#define MAKE_DIR(p) _mkdir(p)
#else
#include <pthread.h>
#include <unistd.h>
#include <sys/stat.h>
#define MAKE_DIR(p) mkdir(p, 0755)
#endif

#include "example_svg.h"

/* ------------------------------------------------------------------ */
/*  The executor: workers taking chunks from a shared counter         */
/* ------------------------------------------------------------------ */

#define RP_THREADS_MAX 64

typedef struct rp_run
{
	qaws_batch_task_fn task;
	void* ctx;
	unsigned int count;
	volatile long next;
} rp_run;

static long rp_take(rp_run* r)
{
#ifdef _WIN32
	return InterlockedIncrement(&r->next) - 1;
#else
	return __atomic_fetch_add(&r->next, 1, __ATOMIC_RELAXED);
#endif
}

static void rp_work(rp_run* r)
{
	long k;
	while ((k = rp_take(r)) < (long)r->count)
		r->task(r->ctx, (unsigned int)k, (unsigned int)k + 1);
}

#ifdef _WIN32
static unsigned __stdcall rp_main(void* p) { rp_work((rp_run*)p); return 0; }
#else
static void* rp_main(void* p) { rp_work((rp_run*)p); return NULL; }
#endif

static void rp_parallel_for(void* user, unsigned int count, qaws_batch_task_fn task, void* ctx)
{
	unsigned int threads = *(unsigned int const*)user, k;
	rp_run r;
#ifdef _WIN32
	HANDLE h[RP_THREADS_MAX];
#else
	pthread_t h[RP_THREADS_MAX];
#endif
	r.task = task; r.ctx = ctx; r.count = count; r.next = 0;
	if (threads > count)
		threads = count;
	for (k = 1; k < threads; k++)
#ifdef _WIN32
		h[k] = (HANDLE)_beginthreadex(NULL, 0, rp_main, &r, 0, NULL);
#else
		pthread_create(&h[k], NULL, rp_main, &r);
#endif
	rp_work(&r);
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

static unsigned int rp_thread_count(void)
{
	long n;
#ifdef _WIN32
	SYSTEM_INFO si;
	GetSystemInfo(&si);
	n = (long)si.dwNumberOfProcessors;
#else
	n = sysconf(_SC_NPROCESSORS_ONLN);
#endif
	return n < 1 ? 1 : n > RP_THREADS_MAX ? RP_THREADS_MAX : (unsigned int)n;
}

/* ------------------------------------------------------------------ */
/*  Scenes                                                            */
/* ------------------------------------------------------------------ */

/* a snapshot's step count meaning: until the solver stops */
#define RP_UNTIL 5000u

typedef struct rp_curve
{
	double* x;             /* 3 per vertex, as read */
	unsigned int* e;
	unsigned int nv, ne;
} rp_curve;

/* an OBJ curve: v and l lines (1-based) */
static int rp_load(char const* dir, char const* name, rp_curve* c)
{
	char path[1024], line[256];
	FILE* f;
	unsigned int cv = 0, ce = 0;
	memset(c, 0, sizeof(*c));
	sprintf(path, "%s/%s", dir, name);
	if (!(f = fopen(path, "r")))
	{
		printf("cannot read %s\n", path);
		return 0;
	}
	while (fgets(line, sizeof(line), f))
	{
		if (line[0] == 'v' && line[1] == ' ') c->nv++;
		else if (line[0] == 'l' && line[1] == ' ') c->ne++;
	}
	c->x = (double*)malloc(sizeof(double) * 3 * (c->nv + 1));
	c->e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (c->ne + 1));
	rewind(f);
	while (c->x && c->e && fgets(line, sizeof(line), f))
	{
		double a, b, d;
		unsigned int i, j;
		if (line[0] == 'v' && line[1] == ' ' && sscanf(line + 2, "%lf %lf %lf", &a, &b, &d) == 3 && cv < c->nv)
		{
			c->x[3 * cv] = a; c->x[3 * cv + 1] = b; c->x[3 * cv + 2] = d; cv++;
		}
		else if (line[0] == 'l' && sscanf(line + 2, "%u %u", &i, &j) == 2 && ce < c->ne && i >= 1 && j >= 1 && i <= c->nv && j <= c->nv)
		{
			c->e[2 * ce] = i - 1; c->e[2 * ce + 1] = j - 1; ce++;
		}
	}
	fclose(f);
	c->ne = ce;
	return c->x && c->e && cv == c->nv && ce > 0;
}

static void rp_free(rp_curve* c) { free(c->x); free(c->e); }

/* solve in dim 3, or 2 on (x, z) (the paper's y-plane scene); snapshots
   after each of steps[0..nsnap) iterations (positions, 3 per vertex) */
static void rp_solve(rp_curve const* c, unsigned int dim, double alpha, double beta, unsigned int const* steps, unsigned int nsnap,
	double** snaps, qaws_repulsive_report* reps)
{
	unsigned int threads = rp_thread_count(), i, k, it = 0;
	qaws_batch_executor ex;
	qaws_repulsive_curves_desc cd;
	qaws_repulsive_curves* rc = NULL;
	qaws_repulsive_curves_solve_desc sd;
	qaws_repulsive_curves_solver* s = NULL;
	qaws_repulsive_report r;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * c->nv);
	memset(&r, 0, sizeof(r));
	for (i = 0; i < c->nv; i++)
		for (k = 0; k < dim; k++) x[dim * i + k] = (qaws_scalar)c->x[3 * i + (dim == 3 ? k : 2 * k)];
	ex.parallel_for = rp_parallel_for; ex.user = &threads;
	memset(&cd, 0, sizeof(cd));
	cd.dimension = dim; cd.vertex_count = c->nv; cd.edges = c->e; cd.edge_count = c->ne;
	cd.alpha = (qaws_scalar)alpha; cd.beta = (qaws_scalar)beta; cd.executor = &ex;
	memset(&sd, 0, sizeof(sd));
	sd.hold_edge_lengths = 1;   /* (the barycenter too: by default, nothing pinned) */
	if (qaws_repulsive_curves_create(&cd, &rc) != QAWS_STATUS_OK || qaws_repulsive_curves_solver_create(rc, &sd, x, &s) != QAWS_STATUS_OK)
		printf("  the scene was refused\n");
	for (k = 0; k < nsnap; k++)
	{
		while (s && it < steps[k] && r.reason == QAWS_REPULSIVE_RUNNING)
		{
			qaws_repulsive_curves_solver_step(s, x, &r);
			it = r.iterations;
		}
		reps[k] = r;
		snaps[k] = (double*)malloc(sizeof(double) * 3 * c->nv);
		for (i = 0; i < c->nv; i++)
		{
			snaps[k][3 * i] = (double)x[dim * i];
			snaps[k][3 * i + 1] = dim == 3 ? (double)x[dim * i + 1] : 0;
			snaps[k][3 * i + 2] = (double)x[dim * i + dim - 1];
		}
		printf("  %u steps: E %.6g, stationarity %.3g, reason %d\n", it, (double)r.energy, (double)r.stationarity, (int)r.reason);
	}
	qaws_repulsive_curves_solver_destroy(s);
	qaws_repulsive_curves_destroy(rc);
	free(x);
}

/* ------------------------------------------------------------------ */
/*  Drawing: a knot diagram, segments far to front, each over a halo   */
/* ------------------------------------------------------------------ */

typedef struct rp_view
{
	double yaw, pitch;     /* radians */
} rp_view;

static void rp_rotate(rp_view const* r, double const* p, double* q)
{
	double cy = cos(r->yaw), sy = sin(r->yaw), cp = cos(r->pitch), sp = sin(r->pitch);
	double x = cy * p[0] + sy * p[1], y = -sy * p[0] + cy * p[1], z = p[2];
	q[0] = x;
	q[1] = cp * z - sp * y;    /* up */
	q[2] = sp * z + cp * y;    /* depth: larger is farther */
}

typedef struct rp_seg
{
	double a[2], b[2], depth;
	unsigned int edge;
} rp_seg;

static int rp_seg_cmp(void const* p, void const* q)
{
	double a = ((rp_seg const*)p)->depth, b = ((rp_seg const*)q)->depth;
	return a < b ? 1 : a > b ? -1 : 0;
}

/* the curve x (3 per vertex) in panel v, fitted, coloured along t (one per
   edge, in [0, 1]; NULL: one colour) */
static void rp_draw(svg* s, viewport const* v, rp_view const* view, double const* x, rp_curve const* c, double const* t, double width)
{
	rp_seg* seg = (rp_seg*)malloc(sizeof(rp_seg) * c->ne);
	double lo[2] = { 1e300, 1e300 }, hi[2] = { -1e300, -1e300 }, dlo = 1e300, dhi = -1e300, scale, cx, cy;
	unsigned int i, k;
	for (i = 0; i < c->ne; i++)
	{
		double p[3], q[3];
		rp_rotate(view, x + 3 * c->e[2 * i], p);
		rp_rotate(view, x + 3 * c->e[2 * i + 1], q);
		seg[i].a[0] = p[0]; seg[i].a[1] = p[1]; seg[i].b[0] = q[0]; seg[i].b[1] = q[1];
		seg[i].depth = 0.5 * (p[2] + q[2]);
		seg[i].edge = i;
		for (k = 0; k < 2; k++)
		{
			lo[k] = fmin(lo[k], fmin(p[k], q[k])); hi[k] = fmax(hi[k], fmax(p[k], q[k]));
		}
		dlo = fmin(dlo, seg[i].depth); dhi = fmax(dhi, seg[i].depth);
	}
	scale = 0.88 * fmin(v->w / fmax(hi[0] - lo[0], 1e-12), v->h / fmax(hi[1] - lo[1], 1e-12));
	cx = v->x0 + 0.5 * v->w - scale * 0.5 * (lo[0] + hi[0]);
	cy = v->y0 + 0.5 * v->h + scale * 0.5 * (lo[1] + hi[1]);
	qsort(seg, c->ne, sizeof(rp_seg), rp_seg_cmp);
	for (i = 0; i < c->ne; i++)
	{
		double ax = cx + scale * seg[i].a[0], ay = cy - scale * seg[i].a[1], bx = cx + scale * seg[i].b[0], by = cy - scale * seg[i].b[1];
		double front = dhi > dlo ? (dhi - seg[i].depth) / (dhi - dlo) : 1, w = width * (0.75 + 0.5 * front);
		char col[32];
		if (t)
			heat(0.1 + 0.8 * t[seg[i].edge], col);
		else
			sprintf(col, "rgb(%d,%d,%d)", (int)(20 + 30 * (1 - front)), (int)(90 + 60 * (1 - front)), (int)(200 + 30 * (1 - front)));
		fprintf(s->f, "<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"#ffffff\" stroke-width=\"%.2f\" stroke-linecap=\"butt\"/>\n",
			ax, ay, bx, by, w + 2.6);
		fprintf(s->f, "<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"%s\" stroke-width=\"%.2f\" stroke-linecap=\"round\"/>\n",
			ax, ay, bx, by, col, w);
	}
	/* a graph's nodes: the vertices not on exactly two edges */
	{
		unsigned char* deg = (unsigned char*)calloc(c->nv, 1);
		for (i = 0; deg && i < 2 * c->ne; i++)
			if (deg[c->e[i]] < 255) deg[c->e[i]]++;
		for (i = 0; deg && i < c->nv; i++)
			if (deg[i] != 2)
			{
				double p[3];
				rp_rotate(view, x + 3 * i, p);
				svg_circle(s, cx + scale * p[0], cy - scale * p[1], 0.9 * width + 1.5, "#24292f", "#ffffff");
			}
		free(deg);
	}
	free(seg);
}

/* each edge's position in the file, in [0, 1], for colouring: the scenes
   list their edges along their chains */
static double* rp_chain_param(rp_curve const* c)
{
	double* t = (double*)calloc(c->ne, sizeof(double));
	unsigned int i;
	for (i = 0; t && i < c->ne; i++) t[i] = c->ne > 1 ? (double)i / (c->ne - 1) : 0;
	return t;
}

/* a figure: nsnap panels of one scene */
static void rp_figure(char const* dir, char const* file, char const* svgname, char const* title, char const* subtitle, unsigned int dim,
	double alpha, double beta, unsigned int const* steps, unsigned int nsnap, rp_view const* view, int colour_along)
{
	rp_curve c;
	double* snaps[8];
	qaws_repulsive_report reps[8];
	svg s;
	double const pw = 300, ph = 300;
	time_t t0 = time(NULL);
	unsigned int k;
	char label[200];
	double* t = NULL;
	if (!rp_load(dir, file, &c))
		return;
	printf("%s: %u vertices, %u edges\n", file, c.nv, c.ne);
	rp_solve(&c, dim, alpha, beta, steps, nsnap, snaps, reps);
	if (colour_along)
		t = rp_chain_param(&c);
	if (svg_open(&s, svgname, (int)(48 + (nsnap + 1) * (pw + 12)), (int)(110 + ph + 30), title, subtitle))
	{
		viewport v;
		v.y0 = 84; v.w = pw; v.h = ph;
		v.x0 = 24; v.xmin = v.ymin = 0; v.xmax = v.ymax = 1;
		sprintf(label, "input: %u vertices, E %.4g", c.nv, (double)reps[0].initial_energy);
		svg_panel(&s, &v, label);
		rp_draw(&s, &v, view, c.x, &c, t, 2.2);
		for (k = 0; k < nsnap; k++)
		{
			v.x0 = 24 + (k + 1) * (pw + 12);
			sprintf(label, "%s%u steps: E %.4g", reps[k].reason == QAWS_REPULSIVE_CONVERGED ? "converged, " : "", reps[k].iterations, (double)reps[k].energy);
			svg_panel(&s, &v, label);
			rp_draw(&s, &v, view, snaps[k], &c, t, 2.2);
		}
		svg_close(&s);
	}
	for (k = 0; k < nsnap; k++) free(snaps[k]);
	free(t);
	rp_free(&c);
	printf("  %s, %.0f s\n", svgname, difftime(time(NULL), t0));
}

int main(int argc, char** argv)
{
	char const* dir = argc > 1 ? argv[1] : "../examples/data/repulsive_curves";
	setvbuf(stdout, NULL, _IONBF, 0);
	MAKE_DIR("showcase");
	{
		static unsigned int const steps[3] = { 300, 400, RP_UNTIL };
		rp_view view = { 0.0, 1.15 };
		rp_figure(dir, "freedmancurve600.obj", "showcase/paper1_freedman.svg", "Repulsive Curves, Fig. 2: Freedman's unknot to the circle",
			"alpha 2, beta 4.5; every edge's length and the barycenter held", 3, 2, 4.5, steps, 3, &view, 1);
	}
	{
		static unsigned int const k5[2] = { 50, RP_UNTIL }, k33[2] = { 200, RP_UNTIL };
		rp_view view = { 0.3, 0.5 };
		rp_figure(dir, "k5.obj", "showcase/paper2_k5.svg", "Repulsive Curves, Fig. 27: K5 embedded in space",
			"a jittered planar drawing; alpha 3, beta 6; every edge's length (its weight) and the barycenter held", 3, 3, 6, k5, 2, &view, 1);
		rp_figure(dir, "k33.obj", "showcase/paper3_k33.svg", "Repulsive Curves, Fig. 27: K3,3 embedded in space",
			"a jittered planar drawing; alpha 3, beta 6; every edge's length (its weight) and the barycenter held", 3, 3, 6, k33, 2, &view, 1);
	}
	{
		static unsigned int const steps[3] = { 300, 450, RP_UNTIL };
		rp_view view = { 0.3, 0.5 };
		rp_figure(dir, "earbuds.obj", "showcase/paper4_earbuds.svg", "Repulsive Curves, Fig. 33: untangling a pair of earbuds",
			"alpha 3, beta 6; every edge's length and the barycenter held", 3, 3, 6, steps, 3, &view, 1);
	}
	{
		static unsigned int const steps[2] = { 100, RP_UNTIL };
		rp_view view = { 0.0, 0.0 };
		rp_figure(dir, "graph-sampled.obj", "showcase/paper5_planar.svg", "Repulsive Curves, Fig. 26: a tree drawn in the plane",
			"its edges as repulsive curves: alpha 3, beta 6; every edge's length and the barycenter held; 2D", 2, 3, 6, steps, 2, &view, 0);
	}
	return 0;
}
