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
 *   paper6_paths.svg      four agents' paths planned in space-time (Fig. 31)
 * and from Repulsive Surfaces (Yu, Brakensiek, Schumacher, Crane, ACM TOG 2021):
 *   paper7_cylinder.svg   a cylinder with Dirichlet ends (Fig. 9, left)
 *   paper8_hopf.svg       a Hopf link to its canonical shape (Fig. 11)
 *   paper9_borromean.svg  Borromean rings to their canonical shape (Fig. 11)
 *   paper10_neck.svg      area minimization regularized: no pinch-off (Fig. 9, centre)
 *
 * Usage: qaws_repulsive_paper [data directory [figure 1 to 10 [most steps]]] (default
 * ../examples/data/repulsive_curves, from working_dir; every figure). Minutes, not
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

/* at most this many steps in every solve (a preview; argument 3) */
static unsigned int rp_cap = RP_UNTIL;


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

/* solve in dim 3, or 2 on (x, z) (the paper's y-plane scene), with sd0's
   settings (NULL: every edge's length held); snapshots after each of
   steps[0..nsnap) iterations (positions, 3 per vertex) */
static void rp_solve(rp_curve const* c, unsigned int dim, double alpha, double beta, unsigned int const* steps, unsigned int nsnap,
	double** snaps, qaws_repulsive_report* reps, qaws_repulsive_curves_solve_desc const* sd0)
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
	if (sd0)
		sd = *sd0;
	else
		sd.hold_edge_lengths = 1;   /* (the barycenter too: by default, nothing pinned) */
	if (qaws_repulsive_curves_create(&cd, &rc) != QAWS_STATUS_OK || qaws_repulsive_curves_solver_create(rc, &sd, x, &s) != QAWS_STATUS_OK)
		printf("  the scene was refused\n");
	for (k = 0; k < nsnap; k++)
	{
		while (s && it < steps[k] && it < rp_cap && r.reason == QAWS_REPULSIVE_RUNNING)
		{
			qaws_repulsive_curves_solver_step(s, x, &r);
			it = r.iterations;
			if (it % 50 == 0 && r.reason == QAWS_REPULSIVE_RUNNING)
				printf("    %u steps: E %.6g, stationarity %.3g\n", it, (double)r.energy, (double)r.stationarity);
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

/* a panel's map from the view's plane: points (3 each) fitted into v */
typedef struct rp_fit
{
	double scale, cx, cy;
} rp_fit;

static void rp_fit_points(viewport const* v, rp_view const* view, double const* pts, unsigned int n, rp_fit* f)
{
	double lo[2] = { 1e300, 1e300 }, hi[2] = { -1e300, -1e300 };
	unsigned int i, k;
	for (i = 0; i < n; i++)
	{
		double p[3];
		rp_rotate(view, pts + 3 * i, p);
		for (k = 0; k < 2; k++) { lo[k] = fmin(lo[k], p[k]); hi[k] = fmax(hi[k], p[k]); }
	}
	f->scale = 0.88 * fmin(v->w / fmax(hi[0] - lo[0], 1e-12), v->h / fmax(hi[1] - lo[1], 1e-12));
	f->cx = v->x0 + 0.5 * v->w - f->scale * 0.5 * (lo[0] + hi[0]);
	f->cy = v->y0 + 0.5 * v->h + f->scale * 0.5 * (lo[1] + hi[1]);
}

/* the curve x (3 per vertex) through the fit, coloured along t (one per
   edge, in [0, 1]; NULL: one colour shaded by depth) */
static void rp_draw_fit(svg* s, rp_fit const* fit, rp_view const* view, double const* x, rp_curve const* c, double const* t, double width)
{
	rp_seg* seg = (rp_seg*)malloc(sizeof(rp_seg) * c->ne);
	double dlo = 1e300, dhi = -1e300, scale = fit->scale, cx = fit->cx, cy = fit->cy;
	unsigned int i;
	for (i = 0; i < c->ne; i++)
	{
		double p[3], q[3];
		rp_rotate(view, x + 3 * c->e[2 * i], p);
		rp_rotate(view, x + 3 * c->e[2 * i + 1], q);
		seg[i].a[0] = p[0]; seg[i].a[1] = p[1]; seg[i].b[0] = q[0]; seg[i].b[1] = q[1];
		seg[i].depth = 0.5 * (p[2] + q[2]);
		seg[i].edge = i;
		dlo = fmin(dlo, seg[i].depth); dhi = fmax(dhi, seg[i].depth);
	}
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

/* the curve fitted into panel v */
static void rp_draw(svg* s, viewport const* v, rp_view const* view, double const* x, rp_curve const* c, double const* t, double width)
{
	rp_fit fit;
	rp_fit_points(v, view, x, c->nv, &fit);
	rp_draw_fit(s, &fit, view, x, c, t, width);
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
	rp_solve(&c, dim, alpha, beta, steps, nsnap, snaps, reps, NULL);
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

/* ------------------------------------------------------------------ */
/*  Meshes (Repulsive Surfaces)                                       */
/* ------------------------------------------------------------------ */

typedef struct rp_mesh
{
	double* x;             /* 3 per vertex */
	unsigned int* t;       /* 3 per face */
	unsigned int nv, nf;
} rp_mesh;

/* an OBJ mesh: v and f lines (1-based, v/vt/vn accepted), polygons fanned */
static int rp_load_mesh(char const* dir, char const* name, rp_mesh* m)
{
	char path[1024], line[1024];
	FILE* f;
	unsigned int cv = 0, cf = 0, pass;
	memset(m, 0, sizeof(*m));
	sprintf(path, "%s/%s", dir, name);
	if (!(f = fopen(path, "r")))
	{
		printf("cannot read %s\n", path);
		return 0;
	}
	for (pass = 0; pass < 2; pass++)
	{
		cv = cf = 0;
		rewind(f);
		while (fgets(line, sizeof(line), f))
		{
			if (line[0] == 'v' && line[1] == ' ')
			{
				double a, b, d;
				if (pass && sscanf(line + 2, "%lf %lf %lf", &a, &b, &d) == 3) { m->x[3 * cv] = a; m->x[3 * cv + 1] = b; m->x[3 * cv + 2] = d; }
				cv++;
			}
			else if (line[0] == 'f' && line[1] == ' ')
			{
				unsigned int idx[32], n = 0, k;
				char* p = line + 2;
				while (n < 32)
				{
					while (*p == ' ' || *p == '\t') p++;
					if (*p < '0' || *p > '9') break;
					idx[n++] = (unsigned int)strtoul(p, &p, 10) - 1;
					while (*p && *p != ' ' && *p != '\t') p++;
				}
				for (k = 1; k + 1 < n; k++, cf++)
					if (pass) { m->t[3 * cf] = idx[0]; m->t[3 * cf + 1] = idx[k]; m->t[3 * cf + 2] = idx[k + 1]; }
			}
		}
		if (!pass)
		{
			m->nv = cv; m->nf = cf;
			m->x = (double*)malloc(sizeof(double) * 3 * (cv + 1));
			m->t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (cf + 1));
			if (!m->x || !m->t)
				break;
		}
	}
	fclose(f);
	for (cf = 0; m->x && m->t && cf < 3 * m->nf; cf++)
		if (m->t[cf] >= m->nv)
			return 0;
	return m->x && m->t && m->nf > 0;
}

static void rp_free_mesh(rp_mesh* m) { free(m->x); free(m->t); }

typedef struct rp_face
{
	double depth;
	unsigned int f;
} rp_face;

static int rp_face_cmp(void const* p, void const* q)
{
	double a = ((rp_face const*)p)->depth, b = ((rp_face const*)q)->depth;
	return a < b ? 1 : a > b ? -1 : 0;
}

/* faces nf (3 vertex indices each) of x through the fit, far to near, flat
   shaded (two-sided) in colour rgb; a face's colour from rgbs[face] if given */
static void rp_draw_mesh(svg* s, rp_fit const* fit, rp_view const* view, double const* x, unsigned int const* tri, unsigned int nf,
	int const* rgb, int const* rgbs)
{
	static double const light[3] = { -0.35, 0.55, -0.76 };   /* (in the view: from the upper left, toward the eye) */
	rp_face* order = (rp_face*)malloc(sizeof(rp_face) * nf);
	unsigned int i, k;
	for (i = 0; order && i < nf; i++)
	{
		double d = 0, q[3];
		for (k = 0; k < 3; k++) { rp_rotate(view, x + 3 * tri[3 * i + k], q); d += q[2]; }
		order[i].depth = d; order[i].f = i;
	}
	if (order)
		qsort(order, nf, sizeof(rp_face), rp_face_cmp);
	for (i = 0; order && i < nf; i++)
	{
		unsigned int const* t = tri + 3 * order[i].f;
		int const* c = rgbs ? rgbs + 3 * order[i].f : rgb;
		double p[3][3], u[3], w[3], n[3], l, sh;
		for (k = 0; k < 3; k++) rp_rotate(view, x + 3 * t[k], p[k]);
		for (k = 0; k < 3; k++) { u[k] = p[1][k] - p[0][k]; w[k] = p[2][k] - p[0][k]; }
		n[0] = u[1] * w[2] - u[2] * w[1]; n[1] = u[2] * w[0] - u[0] * w[2]; n[2] = u[0] * w[1] - u[1] * w[0];
		l = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
		sh = l > 0 ? 0.3 + 0.7 * fabs(n[0] * light[0] + n[1] * light[1] + n[2] * light[2]) / l : 0.3;
		fprintf(s->f, "<polygon points=\"%.2f,%.2f %.2f,%.2f %.2f,%.2f\" fill=\"rgb(%d,%d,%d)\" stroke=\"rgb(%d,%d,%d)\" stroke-width=\"0.35\" stroke-linejoin=\"round\"/>\n",
			fit->cx + fit->scale * p[0][0], fit->cy - fit->scale * p[0][1], fit->cx + fit->scale * p[1][0], fit->cy - fit->scale * p[1][1],
			fit->cx + fit->scale * p[2][0], fit->cy - fit->scale * p[2][1],
			(int)(c[0] * sh), (int)(c[1] * sh), (int)(c[2] * sh), (int)(c[0] * sh * 0.82), (int)(c[1] * sh * 0.82), (int)(c[2] * sh * 0.82));
	}
	free(order);
}

/* an open cylinder: nu around, nh rings from z = -h/2 to h/2, radius r */
static int rp_cylinder(unsigned int nu, unsigned int nh, double r, double h, rp_mesh* m)
{
	unsigned int i, j, f = 0;
	m->nv = nu * nh; m->nf = 2 * nu * (nh - 1);
	m->x = (double*)malloc(sizeof(double) * 3 * m->nv);
	m->t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * m->nf);
	if (!m->x || !m->t)
		return 0;
	for (j = 0; j < nh; j++)
		for (i = 0; i < nu; i++)
		{
			double a = 6.283185307179586 * (i + 0.5 * (j % 2)) / nu;   /* (alternate rings shifted: near-equilateral faces) */
			unsigned int v = j * nu + i;
			m->x[3 * v] = r * cos(a); m->x[3 * v + 1] = r * sin(a); m->x[3 * v + 2] = -0.5 * h + h * j / (nh - 1);
		}
	for (j = 0; j + 1 < nh; j++)
		for (i = 0; i < nu; i++)
		{
			unsigned int a = j * nu + i, b = j * nu + (i + 1) % nu, c = (j + 1) * nu + i, d = (j + 1) * nu + (i + 1) % nu;
			if (j % 2 == 0)
			{
				m->t[3 * f] = a; m->t[3 * f + 1] = b; m->t[3 * f + 2] = c; f++;
				m->t[3 * f] = b; m->t[3 * f + 1] = d; m->t[3 * f + 2] = c; f++;
			}
			else
			{
				m->t[3 * f] = a; m->t[3 * f + 1] = d; m->t[3 * f + 2] = c; f++;
				m->t[3 * f] = a; m->t[3 * f + 1] = b; m->t[3 * f + 2] = d; f++;
			}
		}
	return 1;
}

/* the mesh's vertices on a boundary (an edge of one face) */
static unsigned int* rp_boundary(rp_mesh const* m, unsigned int* out_count)
{
	unsigned int* key = (unsigned int*)malloc(sizeof(unsigned int) * 6 * m->nf);
	unsigned int* pins = (unsigned int*)malloc(sizeof(unsigned int) * (m->nv + 1));
	unsigned char* b = (unsigned char*)calloc(m->nv + 1, 1);
	unsigned int f, k, n = 0, i;
	*out_count = 0;
	if (!key || !pins || !b)
	{
		free(key); free(b);
		return pins;
	}
	for (f = 0; f < m->nf; f++)
		for (k = 0; k < 3; k++)
		{
			unsigned int u = m->t[3 * f + k], v = m->t[3 * f + (k + 1) % 3];
			key[2 * n] = u < v ? u : v; key[2 * n + 1] = u < v ? v : u; n++;
		}
	/* (an edge counted once: quadratic, but run once per figure) */
	for (i = 0; i < n; i++)
	{
		unsigned int c = 0, j;
		for (j = 0; j < n && c < 2; j++) c += key[2 * j] == key[2 * i] && key[2 * j + 1] == key[2 * i + 1];
		if (c == 1) b[key[2 * i]] = b[key[2 * i + 1]] = 1;
	}
	for (i = 0; i < m->nv; i++) if (b[i]) pins[(*out_count)++] = i;
	free(key); free(b);
	return pins;
}

/* each face's colour by its component */
static int* rp_face_colours(rp_mesh const* m)
{
	static int const pal[4][3] = { { 110, 165, 235 }, { 240, 150, 80 }, { 120, 200, 130 }, { 200, 130, 210 } };
	unsigned int* lab = (unsigned int*)malloc(sizeof(unsigned int) * m->nv);
	int* rgb = (int*)malloc(sizeof(int) * 3 * m->nf);
	unsigned int i, k, changed = 1;
	if (!lab || !rgb)
	{
		free(lab); free(rgb);
		return NULL;
	}
	for (i = 0; i < m->nv; i++) lab[i] = i;
	while (changed)
		for (changed = 0, i = 0; i < m->nf; i++)
		{
			unsigned int a = lab[m->t[3 * i]], b = lab[m->t[3 * i + 1]], c = lab[m->t[3 * i + 2]], lo = a < b ? (a < c ? a : c) : (b < c ? b : c);
			for (k = 0; k < 3; k++)
				if (lab[m->t[3 * i + k]] != lo) { lab[m->t[3 * i + k]] = lo; changed = 1; }
		}
	for (i = 0; i < m->nf; i++)
	{
		unsigned int r = 0, v;
		for (v = 0; v < lab[m->t[3 * i]]; v++) r += lab[v] == v;   /* (the component's rank) */
		for (k = 0; k < 3; k++) rgb[3 * i + k] = pal[r % 4][k];
	}
	free(lab);
	return rgb;
}

/* solve the mesh (p, sd's settings) on the
   caller's threads; snapshots after each of steps[0..nsnap) */
static void rp_solve_mesh(rp_mesh const* mesh, double p, qaws_repulsive_mesh_solve_desc const* sd0, unsigned int const* steps, unsigned int nsnap,
	double** snaps, qaws_repulsive_report* reps)
{
	unsigned int threads = rp_thread_count(), i, k, it = 0;
	qaws_batch_executor ex;
	qaws_repulsive_mesh_desc md;
	qaws_repulsive_mesh* rm = NULL;
	qaws_repulsive_mesh_solve_desc sd = *sd0;
	qaws_repulsive_mesh_solver* s = NULL;
	qaws_repulsive_report r;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * mesh->nv);
	memset(&r, 0, sizeof(r));
	for (i = 0; i < 3 * mesh->nv; i++) x[i] = (qaws_scalar)mesh->x[i];
	ex.parallel_for = rp_parallel_for; ex.user = &threads;
	memset(&md, 0, sizeof(md));
	md.vertex_count = mesh->nv; md.triangles = mesh->t; md.triangle_count = mesh->nf; md.p = (qaws_scalar)p; md.executor = &ex;
	if (qaws_repulsive_mesh_create(&md, &rm) != QAWS_STATUS_OK || qaws_repulsive_mesh_solver_create(rm, &sd, x, &s) != QAWS_STATUS_OK)
		printf("  the scene was refused\n");
	for (k = 0; k < nsnap; k++)
	{
		while (s && it < steps[k] && it < rp_cap && r.reason == QAWS_REPULSIVE_RUNNING)
		{
			qaws_repulsive_mesh_solver_step(s, x, &r);
			it = r.iterations;
			if (it % 25 == 0 && r.reason == QAWS_REPULSIVE_RUNNING)
				printf("    %u steps: E %.6g, stationarity %.3g\n", it, (double)r.energy, (double)r.stationarity);
		}
		reps[k] = r;
		snaps[k] = (double*)malloc(sizeof(double) * 3 * mesh->nv);
		for (i = 0; i < 3 * mesh->nv; i++) snaps[k][i] = (double)x[i];
		printf("  %u steps: E %.6g, stationarity %.3g, reason %d\n", it, (double)r.energy, (double)r.stationarity, (int)r.reason);
	}
	qaws_repulsive_mesh_solver_destroy(s);
	qaws_repulsive_mesh_destroy(rm);
	free(x);
}

/* a surfaces figure: the input and nsnap snapshots, every panel fitted to
   the input's view */
/* panels of the mesh at positions xs[0..n), labelled, one scale for all */
static void rp_mesh_panels(char const* svgname, char const* title, char const* subtitle, rp_mesh const* m, double const* const* xs,
	char const (*labels)[200], unsigned int n, rp_view const* view, int const* rgbs)
{
	static int const rgb[3] = { 110, 165, 235 };
	double const pw = 300, ph = 300;
	svg s;
	unsigned int k;
	if (svg_open(&s, svgname, (int)(48 + n * (pw + 12)), (int)(110 + ph + 30), title, subtitle))
	{
		viewport v;
		rp_fit fit;
		double scale = 1e300;
		v.y0 = 84; v.w = pw; v.h = ph; v.xmin = v.ymin = 0; v.xmax = v.ymax = 1; v.x0 = 24;
		for (k = 0; k < n; k++)
		{
			rp_fit_points(&v, view, xs[k], m->nv, &fit);
			scale = fmin(scale, fit.scale);   /* (one scale for every panel: sizes compare) */
		}
		for (k = 0; k < n; k++)
		{
			v.x0 = 24 + k * (pw + 12);
			svg_panel(&s, &v, labels[k]);
			rp_fit_points(&v, view, xs[k], m->nv, &fit);
			fit.cx = v.x0 + 0.5 * v.w + (fit.cx - v.x0 - 0.5 * v.w) * scale / fit.scale;
			fit.cy = v.y0 + 0.5 * v.h + (fit.cy - v.y0 - 0.5 * v.h) * scale / fit.scale;
			fit.scale = scale;
			rp_draw_mesh(&s, &fit, view, xs[k], m->t, m->nf, rgb, rgbs);
		}
		svg_close(&s);
	}
}

static void rp_surface_figure(rp_mesh const* m, char const* svgname, char const* title, char const* subtitle, double p,
	qaws_repulsive_mesh_solve_desc const* sd, unsigned int const* steps, unsigned int nsnap, rp_view const* view)
{
	double* snaps[8];
	qaws_repulsive_report reps[8];
	time_t t0 = time(NULL);
	unsigned int k;
	int* rgbs = rp_face_colours(m);
	double const* xs[9];
	char labels[9][200];
	printf("%s: %u vertices, %u faces\n", svgname, m->nv, m->nf);
	rp_solve_mesh(m, p, sd, steps, nsnap, snaps, reps);
	xs[0] = m->x;
	sprintf(labels[0], "input: %u faces, E %.4g", m->nf, (double)reps[0].initial_energy);
	for (k = 0; k < nsnap; k++)
	{
		xs[k + 1] = snaps[k];
		sprintf(labels[k + 1], "%s%u steps: E %.4g", reps[k].reason == QAWS_REPULSIVE_CONVERGED ? "converged, " : "", reps[k].iterations, (double)reps[k].energy);
	}
	rp_mesh_panels(svgname, title, subtitle, m, xs, labels, nsnap + 1, view, rgbs);
	for (k = 0; k < nsnap; k++) free(snaps[k]);
	free(rgbs);
	printf("  %s, %.0f s\n", svgname, difftime(time(NULL), t0));
}

/* ------------------------------------------------------------------ */
/*  Path planning in space-time (Fig. 31)                             */
/* ------------------------------------------------------------------ */

/* the room of the paper's scene: an H in (x, z), extruded along time y */
static qaws_scalar const rp_h_room[24] = { -3, -3, -1, -3, -1, -1, 1, -1, 1, -3, 3, -3, 3, 3, 1, 3, 1, 1, -1, 1, -1, 3, -3, 3 };

static void rp_line3(svg* s, rp_fit const* fit, rp_view const* view, double const* p, double const* q, char const* color, double width)
{
	double a[3], b[3];
	rp_rotate(view, p, a);
	rp_rotate(view, q, b);
	fprintf(s->f, "<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"%s\" stroke-width=\"%.2f\" stroke-linecap=\"round\"/>\n",
		fit->cx + fit->scale * a[0], fit->cy - fit->scale * a[1], fit->cx + fit->scale * b[0], fit->cy - fit->scale * b[1], color, width);
}

/* the room's prism from t0 to t1 (its two H outlines and its vertical edges) */
static void rp_room3(svg* s, rp_fit const* fit, rp_view const* view, double t0, double t1, int caps_only)
{
	unsigned int j;
	for (j = 0; j < 12; j++)
	{
		qaws_scalar const *a = rp_h_room + 2 * j, *b = rp_h_room + 2 * ((j + 1) % 12);
		double p0[3] = { a[0], t0, a[1] }, q0[3] = { b[0], t0, b[1] }, p1[3] = { a[0], t1, a[1] }, q1[3] = { b[0], t1, b[1] };
		rp_line3(s, fit, view, p0, q0, "#8c959f", 1.2);
		if (!caps_only)
		{
			rp_line3(s, fit, view, p1, q1, "#8c959f", 1.2);
			rp_line3(s, fit, view, p0, p1, "#d0d7de", 1);
		}
	}
}

/* each edge's component (an agent), over the count, for colouring */
static double* rp_components(rp_curve const* c, unsigned int* out_count)
{
	unsigned int* lab = (unsigned int*)malloc(sizeof(unsigned int) * c->nv);
	double* t = (double*)malloc(sizeof(double) * c->ne);
	unsigned int i, n = 0, changed = 1;
	for (i = 0; lab && i < c->nv; i++) lab[i] = i;
	while (lab && changed)
		for (changed = 0, i = 0; i < c->ne; i++)
		{
			unsigned int a = c->e[2 * i], b = c->e[2 * i + 1], m = lab[a] < lab[b] ? lab[a] : lab[b];
			if (lab[a] != m || lab[b] != m) { lab[a] = lab[b] = m; changed = 1; }
		}
	for (i = 0; lab && i < c->nv; i++) n += lab[i] == i;
	for (i = 0; lab && t && i < c->ne; i++)
	{
		unsigned int k = 0, v;
		for (v = 0; v < lab[c->e[2 * i]]; v++) k += lab[v] == v;   /* (the label's rank) */
		t[i] = n > 1 ? (double)k / (n - 1) : 0;
	}
	*out_count = n;
	free(lab);
	return t;
}

/* the smallest distance in (x, z) between two agents at one time: from each
   vertex to the other agents' edges spanning its time y, interpolated */
static double rp_slice_closest(double const* x, rp_curve const* c, double const* comp)
{
	double* col = (double*)malloc(sizeof(double) * c->nv), m = HUGE_VAL;
	unsigned int i, e;
	for (i = 0; col && i < c->ne; i++) col[c->e[2 * i]] = col[c->e[2 * i + 1]] = comp[i];
	for (i = 0; col && i < c->nv; i++)
		for (e = 0; e < c->ne; e++)
		{
			double const *p = x + 3 * c->e[2 * e], *q = x + 3 * c->e[2 * e + 1];
			double t = x[3 * i + 1], s;
			if (comp[e] == col[i] || (t - p[1]) * (t - q[1]) > 0 || p[1] == q[1])
				continue;
			s = (t - p[1]) / (q[1] - p[1]);
			m = fmin(m, hypot(x[3 * i] - (p[0] + s * (q[0] - p[0])), x[3 * i + 2] - (p[2] + s * (q[2] - p[2]))));
		}
	free(col);
	return m;
}

static void rp_paths_figure(char const* dir)
{
	static unsigned int const steps[1] = { 400 };   /* (by then set; the tail to stationarity 1e-4 is slow: 3e-4 at 750) */
	rp_curve c;
	double* snaps[1];
	qaws_repulsive_report reps[1];
	qaws_repulsive_curves_solve_desc sd;
	unsigned int* pins, npin = 0, i, nagents;
	unsigned char* deg;
	double *comp, before, after, prism[3 * 24];
	time_t t0 = time(NULL);
	svg s;
	double const pw = 300, ph = 330;
	char label[200];
	if (!rp_load(dir, "paths-center.obj", &c))
		return;
	printf("paths-center.obj: %u vertices, %u edges\n", c.nv, c.ne);
	/* the paper's scene: the ends pinned (fix_special_vertices), time held,
	   the length weighted (optimize_length 1), the H room's walls repelling
	   (repel_surface h-shape.obj 0.05: twice that here, the time integral of
	   r^-3 being 2 / rho^2) */
	deg = (unsigned char*)calloc(c.nv, 1);
	pins = (unsigned int*)malloc(sizeof(unsigned int) * c.nv);
	for (i = 0; i < 2 * c.ne; i++) deg[c.e[i]]++;
	for (i = 0; i < c.nv; i++) if (deg[i] != 2) pins[npin++] = i;
	memset(&sd, 0, sizeof(sd));
	sd.pinned = pins; sd.pinned_count = npin;
	sd.hold_coordinate = 2; sd.length_weight = 1;
	sd.room = rp_h_room; sd.room_count = 12; sd.room_weight = (qaws_scalar)0.1;
	rp_solve(&c, 3, 3, 6, steps, 1, snaps, reps, &sd);
	comp = rp_components(&c, &nagents);
	before = rp_slice_closest(c.x, &c, comp);
	after = rp_slice_closest(snaps[0], &c, comp);
	printf("  %u agents; the closest in a time slice %.3f -> %.3f\n", nagents, before, after);
	for (i = 0; i < 12; i++)
	{
		prism[3 * i] = prism[3 * (i + 12)] = rp_h_room[2 * i];
		prism[3 * i + 1] = 0; prism[3 * (i + 12) + 1] = 8;
		prism[3 * i + 2] = prism[3 * (i + 12) + 2] = rp_h_room[2 * i + 1];
	}
	if (svg_open(&s, "showcase/paper6_paths.svg", (int)(48 + 3 * (pw + 12)), (int)(110 + ph + 30),
		"Repulsive Curves, Fig. 31: four agents' paths planned in space-time",
		"paths as curves in (x, time, z): time held, the ends pinned, the length weighted, the H room's walls repelling; alpha 3, beta 6"))
	{
		viewport v;
		rp_fit fit;
		rp_view side = { 0.55, -1.2 }, top = { 0, 0 };
		unsigned int k, e;
		v.y0 = 84; v.w = pw; v.h = ph; v.xmin = v.ymin = 0; v.xmax = v.ymax = 1;
		for (k = 0; k < 2; k++)
		{
			v.x0 = 24 + k * (pw + 12);
			if (k) sprintf(label, "%s%u steps: E %.4g; apart %.2f", reps[0].reason == QAWS_REPULSIVE_CONVERGED ? "converged, " : "",
				reps[0].iterations, (double)reps[0].energy, after);
			else sprintf(label, "input: E %.4g; the agents %.2f apart", (double)reps[0].initial_energy, before);
			svg_panel(&s, &v, label);
			rp_fit_points(&v, &side, prism, 24, &fit);
			rp_room3(&s, &fit, &side, 0, 8, 0);
			rp_draw_fit(&s, &fit, &side, k ? snaps[0] : c.x, &c, comp, 2.0);
		}
		/* from above: the paths over the room, the input's in grey */
		v.x0 = 24 + 2 * (pw + 12);
		svg_panel(&s, &v, "from above (along time); the input in grey");
		rp_fit_points(&v, &top, prism, 24, &fit);
		rp_room3(&s, &fit, &top, 0, 8, 1);
		for (e = 0; e < c.ne; e++) rp_line3(&s, &fit, &top, c.x + 3 * c.e[2 * e], c.x + 3 * c.e[2 * e + 1], "#afb8c1", 1.2);
		rp_draw_fit(&s, &fit, &top, snaps[0], &c, comp, 2.0);
		svg_close(&s);
	}
	free(snaps[0]); free(comp); free(pins); free(deg);
	rp_free(&c);
	printf("  showcase/paper6_paths.svg, %.0f s\n", difftime(time(NULL), t0));
}

int main(int argc, char** argv)
{
	char const* dir = argc > 1 ? argv[1] : "../examples/data/repulsive_curves";
	int only = argc > 2 ? atoi(argv[2]) : 0;   /* (one figure: 1 to 10) */
	if (argc > 3)
		rp_cap = (unsigned int)atoi(argv[3]);
	setvbuf(stdout, NULL, _IONBF, 0);
	MAKE_DIR("showcase");
	if (!only || only == 1)
	{
		static unsigned int const steps[3] = { 300, 400, RP_UNTIL };
		rp_view view = { 0.0, 1.15 };
		rp_figure(dir, "freedmancurve600.obj", "showcase/paper1_freedman.svg", "Repulsive Curves, Fig. 2: Freedman's unknot to the circle",
			"alpha 2, beta 4.5; every edge's length and the barycenter held", 3, 2, 4.5, steps, 3, &view, 1);
	}
	if (!only || only == 2 || only == 3)
	{
		static unsigned int const k5[2] = { 50, RP_UNTIL }, k33[2] = { 200, RP_UNTIL };
		rp_view view = { 0.3, 0.5 };
		if (!only || only == 2)
			rp_figure(dir, "k5.obj", "showcase/paper2_k5.svg", "Repulsive Curves, Fig. 27: K5 embedded in space",
				"a jittered planar drawing; alpha 3, beta 6; every edge's length (its weight) and the barycenter held", 3, 3, 6, k5, 2, &view, 1);
		if (!only || only == 3)
			rp_figure(dir, "k33.obj", "showcase/paper3_k33.svg", "Repulsive Curves, Fig. 27: K3,3 embedded in space",
				"a jittered planar drawing; alpha 3, beta 6; every edge's length (its weight) and the barycenter held", 3, 3, 6, k33, 2, &view, 1);
	}
	if (!only || only == 4)
	{
		static unsigned int const steps[3] = { 300, 450, RP_UNTIL };
		rp_view view = { 0.3, 0.5 };
		rp_figure(dir, "earbuds.obj", "showcase/paper4_earbuds.svg", "Repulsive Curves, Fig. 33: untangling a pair of earbuds",
			"alpha 3, beta 6; every edge's length and the barycenter held", 3, 3, 6, steps, 3, &view, 1);
	}
	if (!only || only == 5)
	{
		static unsigned int const steps[2] = { 100, RP_UNTIL };
		rp_view view = { 0.0, 0.0 };
		rp_figure(dir, "graph-sampled.obj", "showcase/paper5_planar.svg", "Repulsive Curves, Fig. 26: a tree drawn in the plane",
			"its edges as repulsive curves: alpha 3, beta 6; every edge's length and the barycenter held; 2D", 2, 3, 6, steps, 2, &view, 0);
	}
	if (!only || only == 6)
		rp_paths_figure(dir);
	if (!only || only == 7)
	{
		/* Repulsive Surfaces, Fig. 9 (left): a cylinder, both rims pinned */
		static unsigned int const steps[2] = { 10, RP_UNTIL };
		rp_view view = { 0.5, 0.45 };
		rp_mesh m;
		qaws_repulsive_mesh_solve_desc sd;
		unsigned int npin;
		unsigned int* pins;
		if (rp_cylinder(48, 19, 1, 2, &m) && (pins = rp_boundary(&m, &npin)) != NULL)
		{
			memset(&sd, 0, sizeof(sd));
			sd.pinned = pins; sd.pinned_count = npin; sd.theta = (qaws_scalar)0.25;
			rp_surface_figure(&m, "showcase/paper7_cylinder.svg", "Repulsive Surfaces, Fig. 9 (left): a cylinder with Dirichlet ends",
				"both rims pinned; p 6; the hierarchical fractional Sobolev descent (theta 0.25)", 6, &sd, steps, 2, &view);
			free(pins);
		}
		rp_free_mesh(&m);
	}
	if (!only || only == 10)
	{
		/* Repulsive Surfaces, Fig. 9 (centre): area minimization between two
		   rings past the catenoid's limit (1.6 apart for radius 1: a minimal
		   surface pinches off), regularized by the tangent-point energy */
		static unsigned int const steps[1] = { RP_UNTIL };
		static double const weights[3] = { 1, 10, 100 };
		rp_view view = { 0.5, 0.45 };
		rp_mesh m;
		qaws_repulsive_mesh_solve_desc sd;
		unsigned int npin, k;
		unsigned int* pins;
		time_t t0 = time(NULL);
		if (rp_cylinder(48, 16, 1, 1.6, &m) && (pins = rp_boundary(&m, &npin)) != NULL)
		{
			double* snaps[3];
			qaws_repulsive_report reps[3];
			double const* xs[4];
			char labels[4][200];
			memset(&sd, 0, sizeof(sd));
			sd.pinned = pins; sd.pinned_count = npin; sd.theta = (qaws_scalar)0.25;
			xs[0] = m.x;
			sprintf(labels[0], "input: %u faces", m.nf);
			for (k = 0; k < 3; k++)
			{
				sd.area_weight = (qaws_scalar)weights[k];
				printf("paper10_neck: the area weighted by %g\n", weights[k]);
				rp_solve_mesh(&m, 6, &sd, steps, 1, snaps + k, reps + k);
				xs[k + 1] = snaps[k];
				sprintf(labels[k + 1], "area weight %g: %u steps, E %.4g", weights[k], reps[k].iterations, (double)reps[k].energy);
			}
			rp_mesh_panels("showcase/paper10_neck.svg", "Repulsive Surfaces, Fig. 9 (centre): area minimization without pinching off",
				"two rings 1.6 apart (radius 1: past the catenoid's limit), pinned; the area weighted, the tangent-point energy (p 6) keeping a neck",
				&m, xs, labels, 4, &view, NULL);
			for (k = 0; k < 3; k++) free(snaps[k]);
			free(pins);
		}
		rp_free_mesh(&m);
		printf("  showcase/paper10_neck.svg, %.0f s\n", difftime(time(NULL), t0));
	}
	/* Repulsive Surfaces, Fig. 11 (the gallery): linked surfaces flowing to
	   their canonical embeddings; the paper's meshes, the area held (the
	   paper's flow holds nothing and is stopped by eye) */
	{
		static char const* const files[2] = { "hopf.obj", "borromean.obj" };
		static char const* const svgs[2] = { "showcase/paper8_hopf.svg", "showcase/paper9_borromean.svg" };
		static char const* const titles[2] = { "Repulsive Surfaces, Fig. 11: a Hopf link", "Repulsive Surfaces, Fig. 11: Borromean rings" };
		/* shortcut: the rings stopped at 400 steps (9 s a step, every step backtracked to 1/32), far from the stationarity floor */
		static unsigned int const steps[2][2] = { { 25, RP_UNTIL }, { 200, 400 } };
		char sdir[1024];
		unsigned int k;
		sprintf(sdir, "%s/../repulsive_surfaces", dir);
		for (k = 0; k < 2; k++)
		{
			rp_view view = { 0.6, 0.5 };
			rp_mesh m;
			qaws_repulsive_mesh_solve_desc sd;
			if (only && only != 8 + (int)k)
				continue;
			if (!rp_load_mesh(sdir, files[k], &m))
				continue;
			memset(&sd, 0, sizeof(sd));
			sd.hold_area = 1; sd.theta = (qaws_scalar)0.25;
			rp_surface_figure(&m, svgs[k], titles[k], "p 6; the area and the barycenter held; the hierarchical fractional Sobolev descent (theta 0.25)", 6, &sd,
				steps[k], 2, &view);
			rp_free_mesh(&m);
		}
	}
	return 0;
}
