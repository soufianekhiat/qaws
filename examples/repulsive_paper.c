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
 *   paper8_hopf.svg       a Hopf link to its canonical shape (a scene of the paper's code)
 *   paper9_borromean.svg  Borromean rings to their canonical shape (a scene of the paper's code)
 *   paper10_neck.svg      area minimization regularized: no pinch-off (Fig. 9, centre)
 *   paper11_handcuffs.svg linked handcuffs to the unlinked double torus (Fig. 1)
 *   paper12_torus_ring.svg a double torus and a ring (a scene of the paper's code)
 *   paper13_figure8.svg   a broken figure-eight tube (a scene of the paper's code)
 *   paper14_walnut.svg    growing a walnut inside its shell (Fig. 22, top)
 *   paper15_shrinkwrap.svg shrink wrapping a point cloud (Fig. 20)
 *   paper16_envelopes.svg nested envelopes around it (Fig. 21; after figure 15)
 *   paper17_cobbles.svg   cobblestones grown in a slab (Fig. 22, bottom)
 *   paper15_shrinkwrap_views.svg the wrap see-through and cut open over the points (with 15; alone: 18)
 *   paper19_vines.svg     vines climbing a vase: curves growing on a surface (an application)
 *
 * Usage: qaws_repulsive_paper [data directory [figure 1 to 19 [most steps]]] (default
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
	if (qaws_repulsive_curves_create(&rc, &cd) != QAWS_STATUS_OK || qaws_repulsive_curves_solver_create(&s, rc, x, &sd) != QAWS_STATUS_OK)
		printf("  the scene was refused\n");
	for (k = 0; k < nsnap; k++)
	{
		while (s && it < steps[k] && it < rp_cap && r.reason == QAWS_REPULSIVE_RUNNING)
		{
			qaws_repulsive_curves_solver_step(&r, x, s);
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
	return m->x && m->t && m->nv > 0;   /* (no faces: a point cloud) */
}

static void rp_free_mesh(rp_mesh* m) { free(m->x); free(m->t); }

/* ---- rendered panels: a z-buffered rasterizer to RGBA PNG ------------ */

/* PNG with stored (uncompressed) deflate blocks: no zlib needed */
static unsigned int rp_crc(unsigned int crc, unsigned char const* p, size_t n)
{
	static unsigned int table[256];
	static int ready = 0;
	size_t i;
	if (!ready)
	{
		unsigned int c, k, b;
		for (k = 0; k < 256; k++)
		{
			for (c = k, b = 0; b < 8; b++) c = c & 1 ? 0xedb88320u ^ (c >> 1) : c >> 1;
			table[k] = c;
		}
		ready = 1;
	}
	crc = ~crc;
	for (i = 0; i < n; i++) crc = table[(crc ^ p[i]) & 0xff] ^ (crc >> 8);
	return ~crc;
}

static void rp_be32(unsigned char* p, unsigned int v)
{
	p[0] = (unsigned char)(v >> 24); p[1] = (unsigned char)(v >> 16); p[2] = (unsigned char)(v >> 8); p[3] = (unsigned char)v;
}

static void rp_png_chunk(FILE* f, char const* type, unsigned char const* data, size_t n)
{
	unsigned char hdr[8];
	rp_be32(hdr, (unsigned int)n);
	memcpy(hdr + 4, type, 4);
	fwrite(hdr, 1, 8, f);
	if (n)
		fwrite(data, 1, n, f);
	rp_be32(hdr, rp_crc(rp_crc(0, hdr + 4, 4), data, n));
	fwrite(hdr, 1, 4, f);
}

static int rp_png_write(char const* path, int w, int h, unsigned char const* rgba)
{
	static unsigned char const sig[8] = { 137, 80, 78, 71, 13, 10, 26, 10 };
	size_t row = 4 * (size_t)w + 1, raw_n = (size_t)h * row, blocks = (raw_n + 65534) / 65535, pos = 0, i;
	unsigned char* raw = (unsigned char*)malloc(raw_n);
	unsigned char* z = (unsigned char*)malloc(2 + raw_n + 5 * blocks + 4);
	unsigned char ihdr[13];
	unsigned int a = 1, b = 0;
	FILE* f = fopen(path, "wb");
	int y;
	if (!f || !raw || !z)
	{
		if (f) fclose(f);
		free(raw); free(z);
		return 0;
	}
	for (y = 0; y < h; y++)
	{
		raw[(size_t)y * row] = 0;
		memcpy(raw + (size_t)y * row + 1, rgba + (size_t)y * 4 * w, 4 * (size_t)w);
	}
	z[pos++] = 0x78; z[pos++] = 0x01;
	for (i = 0; i < raw_n; i += 65535)
	{
		size_t len = raw_n - i < 65535 ? raw_n - i : 65535;
		z[pos++] = (unsigned char)(i + len == raw_n);
		z[pos++] = (unsigned char)(len & 0xff); z[pos++] = (unsigned char)(len >> 8);
		z[pos++] = (unsigned char)(~len & 0xff); z[pos++] = (unsigned char)((~len >> 8) & 0xff);
		memcpy(z + pos, raw + i, len);
		pos += len;
	}
	for (i = 0; i < raw_n; i++) { a = (a + raw[i]) % 65521; b = (b + a) % 65521; }
	rp_be32(z + pos, (b << 16) | a);
	pos += 4;
	fwrite(sig, 1, 8, f);
	rp_be32(ihdr, (unsigned int)w); rp_be32(ihdr + 4, (unsigned int)h);
	ihdr[8] = 8; ihdr[9] = 6; ihdr[10] = 0; ihdr[11] = 0; ihdr[12] = 0;   /* 8-bit RGBA */
	rp_png_chunk(f, "IHDR", ihdr, 13);
	rp_png_chunk(f, "IDAT", z, pos);
	rp_png_chunk(f, "IEND", NULL, 0);
	fclose(f);
	free(raw); free(z);
	return 1;
}

#define RP_SS 3   /* supersamples per output pixel, each way */

/* what a panel adds to its mesh (rp_render): points drawn as small shaded
   balls (a point cloud), the mesh drawn see-through over them, and its
   faces on one side of a plane cut away */
typedef struct rp_extra
{
	double const* points;
	unsigned int npoints;
	double radius;                  /* the balls', in output pixels */
	int prgb[3];
	double alpha;                   /* the mesh's opacity; 0 or 1: opaque */
	double cut[4];                  /* faces whose centroid c has cut . (c, 1) > 0 dropped; all 0: none */
	int const* face_rgb;            /* each face's colour (3 each), over the components' and the tint */
} rp_extra;

/* Blinn-Phong with a key and a fill light, n the unit normal facing the eye */
static void rp_shade(double const* n, int const* c, double shade, double const* key, double const* fill, double const* hk, double* out)
{
	double dk = fmax(0, n[0] * key[0] + n[1] * key[1] + n[2] * key[2]);
	double df = fmax(0, n[0] * fill[0] + n[1] * fill[1] + n[2] * fill[2]);
	double sp = pow(fmax(0, n[0] * hk[0] + n[1] * hk[1] + n[2] * hk[2]), 36);
	unsigned int k;
	for (k = 0; k < 3; k++) out[k] = fmin(255, fmax(0, shade * (c[k] * (0.34 + 0.10 * n[1] + 0.62 * dk + 0.30 * df) + 255 * 0.35 * sp)));
}

/* the faces of x rendered into a w x h RGBA image, transparent where empty:
   pixel (px, py) is the panel point (x0 + px / res, y0 + py / res) under the
   fit. Smooth shaded by the area-weighted vertex normals (Blinn-Phong, a key
   and a fill light); both sides lit, the side facing away darker; a face's
   colour from rgbs[face] if given; ex (may be NULL) adds points, see-through
   and a cut. See-through: the points drawn first, then each pixel's mesh
   layers in front of them averaged and laid over with 1 - (1 - alpha)^layers */
static int rp_render(char const* path, int w, int h, double x0, double y0, double res, rp_fit const* fit, rp_view const* view,
	double const* x, unsigned int nv, unsigned int const* tri, unsigned int nf, int const* rgb, int const* rgbs, rp_extra const* ex)
{
	static double const key[3] = { -0.42, 0.62, -0.66 }, fill[3] = { 0.65, 0.05, -0.76 };   /* (view space: y up, the eye toward -z) */
	int W = w * RP_SS, H = h * RP_SS, px, py;
	int glass = ex && ex->alpha > 0 && ex->alpha < 1, cut = ex && (ex->cut[0] != 0 || ex->cut[1] != 0 || ex->cut[2] != 0);
	double r = res * RP_SS, hk[3], hl;
	double* q = (double*)malloc(sizeof(double) * 3 * nv);
	double* vn = (double*)calloc(3 * (size_t)nv, sizeof(double));
	float* zb = (float*)malloc(sizeof(float) * (size_t)W * H);
	unsigned char* hi = (unsigned char*)calloc(4 * (size_t)W * H, 1);
	unsigned char* out = (unsigned char*)malloc(4 * (size_t)w * h);
	float* acc = glass ? (float*)calloc(3 * (size_t)W * H, sizeof(float)) : NULL;   /* (see-through: the layers' colours summed, and counted) */
	unsigned short* cnt = glass ? (unsigned short*)calloc((size_t)W * H, sizeof(unsigned short)) : NULL;
	unsigned int i, k;
	int ok = 0;
	if (!q || !vn || !zb || !hi || !out || (glass && (!acc || !cnt)))
		goto done;
	hk[0] = key[0]; hk[1] = key[1]; hk[2] = key[2] - 1;   /* (the half vector toward the eye) */
	hl = sqrt(hk[0] * hk[0] + hk[1] * hk[1] + hk[2] * hk[2]);
	for (k = 0; k < 3; k++) hk[k] /= hl;
	for (i = 0; i < nv; i++)
	{
		rp_rotate(view, x + 3 * i, q + 3 * i);
		q[3 * i] = (fit->cx + fit->scale * q[3 * i] - x0) * r;
		q[3 * i + 1] = (fit->cy - fit->scale * q[3 * i + 1] - y0) * r;   /* (pixels: y down) */
		q[3 * i + 2] *= fit->scale * r;
	}
	for (i = 0; i < nf; i++)
	{
		/* area-weighted normals, in view space (y back up) */
		unsigned int const* t = tri + 3 * i;
		double a[3], b[3], n[3];
		for (k = 0; k < 3; k++)
		{
			a[k] = q[3 * t[1] + k] - q[3 * t[0] + k];
			b[k] = q[3 * t[2] + k] - q[3 * t[0] + k];
		}
		a[1] = -a[1]; b[1] = -b[1];
		n[0] = a[1] * b[2] - a[2] * b[1]; n[1] = a[2] * b[0] - a[0] * b[2]; n[2] = a[0] * b[1] - a[1] * b[0];
		for (k = 0; k < 3; k++) { vn[3 * t[0] + k] += n[k]; vn[3 * t[1] + k] += n[k]; vn[3 * t[2] + k] += n[k]; }
	}
	for (i = 0; i < (unsigned int)(W * H); i++) zb[i] = 1e30f;
	for (i = 0; ex && i < ex->npoints; i++)
	{
		/* a point: a ball of the given radius, its depth and normal from its disc */
		double p[3], R = ex->radius * RP_SS, cx, cy;
		int xa, xb, ya, yb;
		rp_rotate(view, ex->points + 3 * i, p);
		cx = (fit->cx + fit->scale * p[0] - x0) * r; cy = (fit->cy - fit->scale * p[1] - y0) * r; p[2] *= fit->scale * r;
		xa = (int)fmax(0, floor(cx - R)); xb = (int)fmin(W - 1, ceil(cx + R)); ya = (int)fmax(0, floor(cy - R)); yb = (int)fmin(H - 1, ceil(cy + R));
		for (py = ya; py <= yb; py++)
			for (px = xa; px <= xb; px++)
			{
				double u = (px + 0.5 - cx) / R, v = (py + 0.5 - cy) / R, s2 = u * u + v * v, n[3], col[3], d;
				unsigned char* o;
				if (s2 >= 1)
					continue;
				n[0] = u; n[1] = -v; n[2] = -sqrt(1 - s2);
				d = p[2] + n[2] * R;
				if (d >= zb[(size_t)py * W + px])
					continue;
				zb[(size_t)py * W + px] = (float)d;
				rp_shade(n, ex->prgb, 1, key, fill, hk, col);
				o = hi + 4 * ((size_t)py * W + px);
				for (k = 0; k < 3; k++) o[k] = (unsigned char)col[k];
				o[3] = 255;
			}
	}
	for (i = 0; i < nf; i++)
	{
		unsigned int const* t = tri + 3 * i;
		double const *A = q + 3 * t[0], *B = q + 3 * t[1], *C = q + 3 * t[2];
		double area = (B[0] - A[0]) * (C[1] - A[1]) - (B[1] - A[1]) * (C[0] - A[0]);
		int const* c = rgbs ? rgbs + 3 * i : rgb;
		int xa, xb, ya, yb;
		if (fabs(area) < 1e-12)
			continue;
		if (cut)
		{
			double s = ex->cut[3];
			for (k = 0; k < 3; k++) s += ex->cut[k] * (x[3 * t[0] + k] + x[3 * t[1] + k] + x[3 * t[2] + k]) / 3;
			if (s > 0)
				continue;
		}
		xa = (int)fmax(0, floor(fmin(A[0], fmin(B[0], C[0])))); xb = (int)fmin(W - 1, ceil(fmax(A[0], fmax(B[0], C[0]))));
		ya = (int)fmax(0, floor(fmin(A[1], fmin(B[1], C[1])))); yb = (int)fmin(H - 1, ceil(fmax(A[1], fmax(B[1], C[1]))));
		for (py = ya; py <= yb; py++)
			for (px = xa; px <= xb; px++)
			{
				double sx = px + 0.5, sy = py + 0.5, l0, l1, l2, d, n[3], nl, col[3], shade = 1;
				unsigned char* o;
				l0 = ((B[0] - sx) * (C[1] - sy) - (B[1] - sy) * (C[0] - sx)) / area;
				l1 = ((C[0] - sx) * (A[1] - sy) - (C[1] - sy) * (A[0] - sx)) / area;
				l2 = 1 - l0 - l1;
				if (l0 < 0 || l1 < 0 || l2 < 0)
					continue;
				d = l0 * A[2] + l1 * B[2] + l2 * C[2];
				if (d >= zb[(size_t)py * W + px])
					continue;
				if (!glass)
					zb[(size_t)py * W + px] = (float)d;
				for (k = 0; k < 3; k++) n[k] = l0 * vn[3 * t[0] + k] + l1 * vn[3 * t[1] + k] + l2 * vn[3 * t[2] + k];
				nl = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
				if (nl == 0)
					continue;
				for (k = 0; k < 3; k++) n[k] /= nl;
				if (n[2] > 0)
				{
					for (k = 0; k < 3; k++) n[k] = -n[k];   /* (the side facing away, as seen) */
					shade = 0.72;
				}
				rp_shade(n, c, shade, key, fill, hk, col);
				if (glass)
				{
					float* a = acc + 3 * ((size_t)py * W + px);
					for (k = 0; k < 3; k++) a[k] += (float)col[k];
					if (cnt[(size_t)py * W + px] < 65535) cnt[(size_t)py * W + px]++;
					continue;
				}
				o = hi + 4 * ((size_t)py * W + px);
				for (k = 0; k < 3; k++) o[k] = (unsigned char)col[k];
				o[3] = 255;
			}
	}
	for (i = 0; glass && i < (unsigned int)(W * H); i++)
	{
		/* the see-through layers over what lies behind them */
		unsigned char* o = hi + 4 * (size_t)i;
		double T, m;
		if (!cnt[i])
			continue;
		T = pow(1 - ex->alpha, cnt[i]);
		for (k = 0; k < 3; k++)
		{
			m = acc[3 * (size_t)i + k] / cnt[i];
			o[k] = (unsigned char)fmin(255, o[3] ? o[k] * T + m * (1 - T) : m);
		}
		o[3] = (unsigned char)(o[3] ? 255 : 255 * (1 - T));
	}
	/* box-filter the supersamples, colour weighted by coverage */
	for (py = 0; py < h; py++)
		for (px = 0; px < w; px++)
		{
			double s[4] = { 0, 0, 0, 0 };
			int u, v;
			for (v = 0; v < RP_SS; v++)
				for (u = 0; u < RP_SS; u++)
				{
					unsigned char const* p = hi + 4 * ((size_t)(py * RP_SS + v) * W + px * RP_SS + u);
					for (k = 0; k < 3; k++) s[k] += p[k] * (p[3] / 255.0);
					s[3] += p[3];
				}
			for (k = 0; k < 3; k++) out[4 * ((size_t)py * w + px) + k] = (unsigned char)(s[3] > 0 ? fmin(255, s[k] * 255 / s[3] + 0.5) : 0);
			out[4 * ((size_t)py * w + px) + 3] = (unsigned char)(s[3] / (RP_SS * RP_SS) + 0.5);
		}
	ok = rp_png_write(path, w, h, out);
done:
	free(q); free(vn); free(zb); free(hi); free(out); free(acc); free(cnt);
	return ok;
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
	static int const pal[4][3] = { { 72, 132, 222 }, { 236, 140, 64 }, { 120, 200, 130 }, { 200, 130, 210 } };
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

/* ---- dynamic remeshing, as the paper's program does after every step -- */

typedef struct rp_edge_key { unsigned int lo, hi, h; } rp_edge_key;

static int rp_edge_cmp(void const* a, void const* b)
{
	rp_edge_key const *p = (rp_edge_key const*)a, *q = (rp_edge_key const*)b;
	return p->lo != q->lo ? (p->lo < q->lo ? -1 : 1) : p->hi != q->hi ? (p->hi < q->hi ? -1 : 1) : 0;
}

/* half-edge h = 3 f + c runs from t[h] to the face's next corner; twin[h]
   is the opposite one (a closed manifold mesh); keys sorted by edge */
static void rp_twins(rp_mesh const* m, rp_edge_key* keys, unsigned int* twin)
{
	unsigned int h, n = 3 * m->nf;
	for (h = 0; h < n; h++)
	{
		unsigned int a = m->t[h], b = m->t[3 * (h / 3) + (h + 1) % 3];
		keys[h].lo = a < b ? a : b; keys[h].hi = a < b ? b : a; keys[h].h = h;
	}
	qsort(keys, n, sizeof(rp_edge_key), rp_edge_cmp);
	for (h = 0; h + 1 < n; h += 2) { twin[keys[h].h] = keys[h + 1].h; twin[keys[h + 1].h] = keys[h].h; }
}

static double rp_corner_angle(double const* a, double const* b, double const* c)   /* at a */
{
	double u[3], v[3], cr[3];
	unsigned int k;
	for (k = 0; k < 3; k++) { u[k] = b[k] - a[k]; v[k] = c[k] - a[k]; }
	cr[0] = u[1] * v[2] - u[2] * v[1]; cr[1] = u[2] * v[0] - u[0] * v[2]; cr[2] = u[0] * v[1] - u[1] * v[0];
	return atan2(sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]), u[0] * v[0] + u[1] * v[1] + u[2] * v[2]);
}

/* flips every edge whose opposite angles sum past pi, in passes that touch a
   vertex once each (so the twins stay valid), until none is left */
static void rp_delaunay(rp_mesh* m, rp_edge_key* keys, unsigned int* twin, unsigned char* touched)
{
	unsigned int pass, flips = 1;
	for (pass = 0; pass < 100 && flips; pass++)
	{
		unsigned int h;
		rp_twins(m, keys, twin);
		memset(touched, 0, m->nv);
		for (flips = 0, h = 0; h < 3 * m->nf; h++)
		{
			unsigned int g = twin[h], f = h / 3, e = g / 3;
			unsigned int u = m->t[h], v = m->t[3 * f + (h + 1) % 3], w = m->t[3 * f + (h + 2) % 3], z = m->t[3 * e + (g + 2) % 3];
			rp_edge_key key;
			if (g < h || touched[u] || touched[v] || touched[w] || touched[z] || w == z)
				continue;
			if (rp_corner_angle(m->x + 3 * w, m->x + 3 * u, m->x + 3 * v) + rp_corner_angle(m->x + 3 * z, m->x + 3 * u, m->x + 3 * v) <= PI)
				continue;
			key.lo = w < z ? w : z; key.hi = w < z ? z : w;
			if (bsearch(&key, keys, 3 * m->nf, sizeof(rp_edge_key), rp_edge_cmp))
				continue;   /* (w z already an edge: the flip would make the mesh non-manifold) */
			m->t[3 * f] = w; m->t[3 * f + 1] = u; m->t[3 * f + 2] = z;
			m->t[3 * e] = z; m->t[3 * e + 1] = v; m->t[3 * e + 2] = w;
			touched[u] = touched[v] = touched[w] = touched[z] = 1;
			flips++;
		}
	}
}

static void rp_face_normal(rp_mesh const* m, double const* x, unsigned int f, double* n)
{
	unsigned int const* t = m->t + 3 * f;
	double const *p = x + 3 * t[0], *q = x + 3 * t[1], *r = x + 3 * t[2];
	double u[3] = { q[0] - p[0], q[1] - p[1], q[2] - p[2] }, w[3] = { r[0] - p[0], r[1] - p[1], r[2] - p[2] };
	n[0] = u[1] * w[2] - u[2] * w[1]; n[1] = u[2] * w[0] - u[0] * w[2]; n[2] = u[0] * w[1] - u[1] * w[0];
}

/* each vertex half way toward the area-weighted circumcentres of its faces,
   in its tangent plane (the reference's smoothByCircumcenter); a face that
   would turn by more than 90 degrees keeps its vertices where they were
   (an obtuse face's circumcentre lies outside it: unguarded, the smoothing
   can fold a face through its neighbours) */
static void rp_smooth(rp_mesh* m, double* acc)
{
	unsigned int f, v, k, round, reverted = 1;
	double* nrm = acc + 3 * m->nv;
	double* area = nrm + 3 * m->nv;
	double* old = area + m->nv;
	memcpy(old, m->x, sizeof(double) * 3 * m->nv);
	memset(acc, 0, sizeof(double) * 7 * m->nv);
	for (f = 0; f < m->nf; f++)
	{
		unsigned int const* t = m->t + 3 * f;
		double const *p1 = m->x + 3 * t[0], *p2 = m->x + 3 * t[1], *p3 = m->x + 3 * t[2];
		double u[3], w[3], n[3], a2 = 0, b2 = 0, c2 = 0, o[3], A, s;
		for (k = 0; k < 3; k++)
		{
			u[k] = p2[k] - p1[k]; w[k] = p3[k] - p1[k];
			a2 += (p3[k] - p2[k]) * (p3[k] - p2[k]); b2 += w[k] * w[k]; c2 += u[k] * u[k];
		}
		n[0] = u[1] * w[2] - u[2] * w[1]; n[1] = u[2] * w[0] - u[0] * w[2]; n[2] = u[0] * w[1] - u[1] * w[0];
		A = 0.5 * sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
		o[0] = a2 * (b2 + c2 - a2); o[1] = b2 * (c2 + a2 - b2); o[2] = c2 * (a2 + b2 - c2);
		s = o[0] + o[1] + o[2];
		if (!(A > 0) || s == 0)
			continue;
		for (v = 0; v < 3; v++)
			for (k = 0; k < 3; k++)
			{
				double cc = (o[0] * p1[k] + o[1] * p2[k] + o[2] * p3[k]) / s;
				acc[3 * t[v] + k] += A * (cc - m->x[3 * t[v] + k]);
				nrm[3 * t[v] + k] += n[k];
			}
		for (v = 0; v < 3; v++) area[t[v]] += A;
	}
	for (v = 0; v < m->nv; v++)
	{
		double* d = acc + 3 * v;
		double const* n = nrm + 3 * v;
		double nn = n[0] * n[0] + n[1] * n[1] + n[2] * n[2], dn;
		if (!(area[v] > 0) || !(nn > 0))
			continue;
		dn = (d[0] * n[0] + d[1] * n[1] + d[2] * n[2]) / nn;
		for (k = 0; k < 3; k++) m->x[3 * v + k] += 0.5 * (d[k] - dn * n[k]) / area[v];
	}
	for (round = 0; round < 20 && reverted; round++)
		for (reverted = 0, f = 0; f < m->nf; f++)
		{
			double a[3], b[3];
			rp_face_normal(m, old, f, a);
			rp_face_normal(m, m->x, f, b);
			if (a[0] * b[0] + a[1] * b[1] + a[2] * b[2] > 0)
				continue;
			for (v = 0; v < 3; v++)
			{
				unsigned int i = m->t[3 * f + v];
				reverted += memcmp(m->x + 3 * i, old + 3 * i, sizeof(double) * 3) != 0;
				memcpy(m->x + 3 * i, old + 3 * i, sizeof(double) * 3);
			}
		}
}

static double rp_mesh_area(rp_mesh const* m)
{
	double a = 0;
	unsigned int f;
	for (f = 0; f < m->nf; f++)
	{
		double n[3];
		rp_face_normal(m, m->x, f, n);
		a += 0.5 * sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
	}
	return a;
}

static double rp_dist(double const* a, double const* b)
{
	return sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
}

/* edges longer than 1.5 L split at their midpoints, each face once a pass
   (so the twins stay valid), passes until none is left; 0: no memory */
static int rp_split(rp_mesh* m, double L)
{
	unsigned int pass, splits = 1;
	for (pass = 0; pass < 8 && splits; pass++)
	{
		unsigned int nf = m->nf, h;
		rp_edge_key* keys = (rp_edge_key*)malloc(sizeof(rp_edge_key) * 3 * nf);
		unsigned int* twin = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nf);
		unsigned char* done = (unsigned char*)calloc(nf, 1);
		/* (at most nf / 2 splits a pass, each 1 vertex and 2 faces) */
		double* x = (double*)realloc(m->x, sizeof(double) * 3 * (m->nv + nf / 2 + 1));
		unsigned int* t = x ? (unsigned int*)realloc(m->t, sizeof(unsigned int) * 3 * (2 * nf + 1)) : NULL;
		if (x) m->x = x;
		if (t) m->t = t;
		if (!keys || !twin || !done || !x || !t)
		{
			free(keys); free(twin); free(done);
			return 0;
		}
		rp_twins(m, keys, twin);
		for (splits = 0, h = 0; h < 3 * nf; h++)
		{
			unsigned int g = twin[h], f = h / 3, e = g / 3, u, v, w, z, c = m->nv, q;
			if (g < h || done[f] || done[e])
				continue;
			u = t[h]; v = t[3 * f + (h + 1) % 3]; w = t[3 * f + (h + 2) % 3]; z = t[3 * e + (g + 2) % 3];
			if (!(rp_dist(x + 3 * u, x + 3 * v) > 1.5 * L))
				continue;
			for (q = 0; q < 3; q++) x[3 * c + q] = 0.5 * (x[3 * u + q] + x[3 * v + q]);
			m->nv++;
			t[3 * f] = u; t[3 * f + 1] = c; t[3 * f + 2] = w;
			t[3 * m->nf] = c; t[3 * m->nf + 1] = v; t[3 * m->nf + 2] = w; m->nf++;
			t[3 * e] = v; t[3 * e + 1] = c; t[3 * e + 2] = z;
			t[3 * m->nf] = c; t[3 * m->nf + 1] = u; t[3 * m->nf + 2] = z; m->nf++;
			done[f] = done[e] = 1;
			splits++;
		}
		free(keys); free(twin); free(done);
	}
	return 1;
}

/* the vertices around v (its faces' other corners), up to cap; their count */
static unsigned int rp_ring(rp_mesh const* m, unsigned int const* vf, unsigned int const* vo, unsigned int v, unsigned int* out, unsigned int cap)
{
	unsigned int n = 0, i, c, j;
	for (i = vo[v]; i < vo[v + 1]; i++)
		for (c = 0; c < 3; c++)
		{
			unsigned int a = m->t[3 * vf[i] + c];
			for (j = 0; j < n && out[j] != a; j++) {}
			if (a != v && j == n && n < cap) out[n++] = a;
		}
	return n;
}

/* edges shorter than L / 2 collapsed to their midpoints, passes until none
   is left: an edge u v only when its ends' rings share just its two
   opposite corners (the link condition: the mesh stays manifold), those
   keep more than 3 faces, and no face around turns over; a collapse's
   neighbourhood untouched for the rest of its pass. 0: no memory */
static int rp_collapse(rp_mesh* m, double L)
{
	enum { RING = 64 };
	unsigned int pass, collapses = 1;
	for (pass = 0; pass < 8 && collapses; pass++)
	{
		unsigned int nv = m->nv, nf = m->nf, h, i, c, nn, nk, a;
		rp_edge_key* keys = (rp_edge_key*)malloc(sizeof(rp_edge_key) * 3 * nf);
		unsigned int* twin = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nf);
		unsigned int* vo = (unsigned int*)calloc(nv + 1, sizeof(unsigned int));
		unsigned int* vf = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nf);
		unsigned int* remap = (unsigned int*)malloc(sizeof(unsigned int) * (nv + 1));
		unsigned char* busy = (unsigned char*)calloc(nv, 1);
		unsigned char* dead = (unsigned char*)calloc(nf, 1);
		unsigned char* gone = (unsigned char*)calloc(nv, 1);
		if (!keys || !twin || !vo || !vf || !remap || !busy || !dead || !gone)
		{
			free(keys); free(twin); free(vo); free(vf); free(remap); free(busy); free(dead); free(gone);
			return 0;
		}
		rp_twins(m, keys, twin);
		for (i = 0; i < 3 * nf; i++) vo[m->t[i] + 1]++;
		for (i = 0; i < nv; i++) vo[i + 1] += vo[i];
		for (i = 0; i < nf; i++)
			for (c = 0; c < 3; c++) vf[vo[m->t[3 * i + c]]++] = i;
		for (i = nv; i > 0; i--) vo[i] = vo[i - 1];   /* (the offsets back to the starts) */
		vo[0] = 0;
		for (collapses = 0, h = 0; h < 3 * nf; h++)
		{
			unsigned int g = twin[h], f = h / 3, e = g / 3, u, v, w, z, ru[RING], rv[RING], a, b, shared = 0, ok = 1;
			double mid[3];
			if (g < h || dead[f] || dead[e])
				continue;
			u = m->t[h]; v = m->t[3 * f + (h + 1) % 3]; w = m->t[3 * f + (h + 2) % 3]; z = m->t[3 * e + (g + 2) % 3];
			if (busy[u] || busy[v] || !(rp_dist(m->x + 3 * u, m->x + 3 * v) < 0.5 * L) || vo[w + 1] - vo[w] <= 3 || vo[z + 1] - vo[z] <= 3)
				continue;
			nn = rp_ring(m, vf, vo, u, ru, RING);
			nk = rp_ring(m, vf, vo, v, rv, RING);
			for (a = 0; a < nn && ok; a++)
			{
				ok = !busy[ru[a]];
				for (b = 0; b < nk; b++) shared += ru[a] == rv[b];
			}
			for (b = 0; b < nk && ok; b++) ok = !busy[rv[b]];
			if (!ok || shared != 2 || nn >= RING || nk >= RING)
				continue;
			/* no face around u or v (but the two going) turns over */
			for (a = 0; a < 3; a++) mid[a] = 0.5 * (m->x[3 * u + a] + m->x[3 * v + a]);
			for (a = 0; a < 2 && ok; a++)
				for (i = vo[a ? v : u]; i < vo[(a ? v : u) + 1] && ok; i++)
				{
					unsigned int fi = vf[i], k;
					double n0[3], n1[3], y[9];
					if (fi == f || fi == e)
						continue;
					rp_face_normal(m, m->x, fi, n0);
					for (k = 0; k < 3; k++)
					{
						unsigned int p = m->t[3 * fi + k];
						for (c = 0; c < 3; c++) y[3 * k + c] = p == u || p == v ? mid[c] : m->x[3 * p + c];
					}
					{
						double e1[3] = { y[3] - y[0], y[4] - y[1], y[5] - y[2] }, e2[3] = { y[6] - y[0], y[7] - y[1], y[8] - y[2] };
						n1[0] = e1[1] * e2[2] - e1[2] * e2[1]; n1[1] = e1[2] * e2[0] - e1[0] * e2[2]; n1[2] = e1[0] * e2[1] - e1[1] * e2[0];
					}
					ok = n0[0] * n1[0] + n0[1] * n1[1] + n0[2] * n1[2] > 0;
				}
			if (!ok)
				continue;
			for (c = 0; c < 3; c++) m->x[3 * u + c] = mid[c];
			for (i = vo[v]; i < vo[v + 1]; i++)
				for (c = 0; c < 3; c++) if (m->t[3 * vf[i] + c] == v) m->t[3 * vf[i] + c] = u;
			dead[f] = dead[e] = 1; gone[v] = 1;
			busy[u] = busy[v] = 1;
			for (a = 0; a < nn; a++) busy[ru[a]] = 1;
			for (b = 0; b < nk; b++) busy[rv[b]] = 1;
			collapses++;
		}
		/* the dead faces and gone vertices dropped */
		for (i = 0, c = 0; i < nv; i++)
			if (!gone[i])
			{
				remap[i] = c;
				for (a = 0; a < 3; a++) m->x[3 * c + a] = m->x[3 * i + a];
				c++;
			}
		m->nv = c;
		for (i = 0, c = 0; i < nf; i++)
			if (!dead[i])
			{
				for (a = 0; a < 3; a++) m->t[3 * c + a] = remap[m->t[3 * i + a]];
				c++;
			}
		m->nf = c;
		free(keys); free(twin); free(vo); free(vf); free(remap); free(busy); free(dead); free(gone);
	}
	return 1;
}

/* the reference's Remesh(5, true) (its interface's default after every
   step): edges split above 1.5 L and collapsed below L / 2 (L the input's
   mean edge), then five rounds of smoothing and Delaunay flips; then the
   area scaled back to what it was before (the smoothing moves it slightly)
   about the vertex mean */
static int rp_remesh(rp_mesh* m, double L)
{
	rp_edge_key* keys;
	unsigned int* twin;
	unsigned char* touched;
	double* acc;
	double c[3] = { 0, 0, 0 }, s, area0 = rp_mesh_area(m);
	unsigned int i, k;
	int ok = rp_split(m, L) && rp_collapse(m, L);
	keys = (rp_edge_key*)malloc(sizeof(rp_edge_key) * 3 * m->nf);
	twin = (unsigned int*)malloc(sizeof(unsigned int) * 3 * m->nf);
	touched = (unsigned char*)malloc(m->nv);
	acc = (double*)malloc(sizeof(double) * 10 * m->nv);
	ok = ok && keys && twin && touched && acc;
	for (i = 0; ok && i < 5; i++)
	{
		rp_smooth(m, acc);
		rp_delaunay(m, keys, twin, touched);
	}
	free(keys); free(twin); free(touched); free(acc);
	if (!ok)
		return 0;
	for (i = 0; i < m->nv; i++) for (k = 0; k < 3; k++) c[k] += m->x[3 * i + k] / m->nv;
	s = sqrt(area0 / rp_mesh_area(m));
	for (i = 0; i < m->nv; i++) for (k = 0; k < 3; k++) m->x[3 * i + k] = c[k] + s * (m->x[3 * i + k] - c[k]);
	return 1;
}

/* a small deterministic generator: [0, 1) */
static double rp_rand(unsigned int* state)
{
	*state = *state * 1664525u + 1013904223u;
	return (*state >> 8) / 16777216.0;
}

/* the cobblestones' scene: n stones scattered without overlap (the seed
   fixes them) in the slab |x|, |y| <= w / 2, |z| <= h / 2, as a wall built
   by hand: each an ellipsoid (an octahedron split twice, onto the sphere)
   of radius 0.18 to 0.36 times w / 4, stretched by up to 12 % across, its
   height 0.4 to 0.65 of that, turned about z; and the slab, its six sides
   as open grids of cells about w / 40 wide (an obstacle: small faces near
   the stones, where one point per face must do) */
static int rp_cobbles(unsigned int n, double w, double h, unsigned int seed, rp_mesh* m, rp_mesh* box)
{
	static double const ov[6][3] = { { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 }, { 0, 0, 1 }, { 0, 0, -1 } };
	static unsigned int const of[8][3] = { { 0, 2, 4 }, { 2, 1, 4 }, { 1, 3, 4 }, { 3, 0, 4 }, { 2, 0, 5 }, { 1, 2, 5 }, { 3, 1, 5 }, { 0, 3, 5 } };
	double sx[66 * 3], half = 0.5 * w, s = 0.25 * w, (*stone)[7];
	unsigned int st[128 * 3], t2[128 * 3], nv = 6, nf = 8, l, i, j, k, q, side, f = 0, v = 0, placed = 0, tries;
	unsigned int cells[3], gn[6][2];
	memcpy(sx, ov, sizeof(ov)); memcpy(st, of, sizeof(of));
	for (l = 0; l < 2; l++)
	{
		unsigned int nf2 = 0;
		for (i = 0; i < nf; i++)
		{
			unsigned int mid[3];
			for (k = 0; k < 3; k++)
			{
				unsigned int a = st[3 * i + k], b = st[3 * i + (k + 1) % 3];
				double p[3], len;
				for (q = 0; q < 3; q++) p[q] = 0.5 * (sx[3 * a + q] + sx[3 * b + q]);
				len = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
				for (q = 0; q < 3; q++) p[q] /= len;
				for (j = 0; j < nv && rp_dist(sx + 3 * j, p) > 1e-12; j++) {}
				if (j == nv) { for (q = 0; q < 3; q++) sx[3 * nv + q] = p[q]; nv++; }
				mid[k] = j;
			}
			{
				unsigned int const c0 = st[3 * i], c1 = st[3 * i + 1], c2 = st[3 * i + 2];
				unsigned int nt[12] = { c0, mid[0], mid[2], mid[0], c1, mid[1], mid[2], mid[1], c2, mid[0], mid[1], mid[2] };
				memcpy(t2 + 3 * nf2, nt, sizeof(nt));
				nf2 += 4;
			}
		}
		nf = nf2;
		memcpy(st, t2, sizeof(unsigned int) * 3 * nf);
	}
	/* the stones: centre x, y, z, a radius across r, the other over r, the height over r, the turn; placed where none overlaps */
	stone = (double(*)[7])malloc(sizeof(double[7]) * (n + 1));
	m->nv = n * nv; m->nf = n * nf;
	m->x = (double*)malloc(sizeof(double) * 3 * m->nv);
	m->t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * m->nf);
	if (!stone || !m->x || !m->t)
	{
		free(stone);
		return 0;
	}
	for (tries = 0; placed < n && tries < 100000; tries++)
	{
		double r = s * (0.18 + 0.18 * rp_rand(&seed)), c[2], ok = 1;
		c[0] = -half + r + (w - 2 * r) * rp_rand(&seed);
		c[1] = -half + r + (w - 2 * r) * rp_rand(&seed);
		for (j = 0; j < placed && ok; j++)
			ok = hypot(c[0] - stone[j][0], c[1] - stone[j][1]) > r + stone[j][3] + 0.04 * s;
		if (!ok)
			continue;
		stone[placed][0] = c[0]; stone[placed][1] = c[1]; stone[placed][2] = 0.04 * s * (rp_rand(&seed) - 0.5);
		stone[placed][3] = r;
		stone[placed][4] = 0.88 + 0.12 * rp_rand(&seed);                 /* (the other radius across, over r) */
		stone[placed][5] = 0.4 + 0.25 * rp_rand(&seed);                   /* (the height over r) */
		stone[placed][6] = 6.283185307179586 * rp_rand(&seed);             /* (the turn about z) */
		placed++;
	}
	if (placed < n)
	{
		free(stone);
		return 0;
	}
	for (i = 0; i < n; i++)
	{
		double r = stone[i][3], ra = r, rb = r * stone[i][4], hz = r * stone[i][5], ang = stone[i][6];
		double ca = cos(ang), sa = sin(ang);
		for (j = 0; j < nv; j++)
		{
			double px = ra * sx[3 * j], py = rb * sx[3 * j + 1];
			double* p = m->x + 3 * (i * nv + j);
			p[0] = stone[i][0] + ca * px - sa * py; p[1] = stone[i][1] + sa * px + ca * py; p[2] = stone[i][2] + hz * sx[3 * j + 2];
		}
		for (j = 0; j < 3 * nf; j++) m->t[3 * i * nf + j] = st[j] + i * nv;
	}
	free(stone);
	/* the box: per side, its two in-plane axes' cell counts */
	cells[0] = cells[1] = (unsigned int)ceil(2 * half / (0.1 * s)); cells[2] = (unsigned int)ceil(h / (0.1 * s));
	for (side = 0; side < 6; side++) { gn[side][0] = cells[side < 2 ? 0 : side < 4 ? 1 : 0]; gn[side][1] = cells[side < 2 ? 1 : 2]; }
	box->nv = box->nf = 0;
	for (side = 0; side < 6; side++) { box->nv += (gn[side][0] + 1) * (gn[side][1] + 1); box->nf += 2 * gn[side][0] * gn[side][1]; }
	box->x = (double*)malloc(sizeof(double) * 3 * box->nv);
	box->t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * box->nf);
	if (!box->x || !box->t)
		return 0;
	for (side = 0; side < 6; side++)
	{
		/* sides 0, 1: z = -+h/2 over x, y; 2, 3: x = -+half over y, z; 4, 5: y = -+half over x, z */
		unsigned int a, b, base = v, na = gn[side][0], nb = gn[side][1];
		double sg = side % 2 ? 1 : -1;
		for (b = 0; b <= nb; b++)
			for (a = 0; a <= na; a++, v++)
			{
				double u = (double)a / na, w = (double)b / nb, * p = box->x + 3 * v;
				if (side < 2) { p[0] = -half + 2 * half * u; p[1] = -half + 2 * half * w; p[2] = sg * 0.5 * h; }
				else if (side < 4) { p[0] = sg * half; p[1] = -half + 2 * half * u; p[2] = -0.5 * h + h * w; }
				else { p[0] = -half + 2 * half * u; p[1] = sg * half; p[2] = -0.5 * h + h * w; }
			}
		for (b = 0; b < nb; b++)
			for (a = 0; a < na; a++)
			{
				unsigned int c0 = base + b * (na + 1) + a, c1 = c0 + 1, c2 = c0 + na + 1, c3 = c2 + 1;
				box->t[3 * f] = c0; box->t[3 * f + 1] = c1; box->t[3 * f + 2] = c3; f++;
				box->t[3 * f] = c0; box->t[3 * f + 1] = c3; box->t[3 * f + 2] = c2; f++;
			}
	}
	return 1;
}

static int rp_copy_mesh(rp_mesh const* a, rp_mesh* b)
{
	b->nv = a->nv; b->nf = a->nf;
	b->x = (double*)malloc(sizeof(double) * 3 * a->nv);
	b->t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * a->nf);
	if (!b->x || !b->t)
	{
		rp_free_mesh(b);
		return 0;
	}
	memcpy(b->x, a->x, sizeof(double) * 3 * a->nv);
	memcpy(b->t, a->t, sizeof(unsigned int) * 3 * a->nf);
	return 1;
}

/* how a surface runs beyond its solve settings */
typedef struct rp_opts
{
	int remesh;                     /* remeshed after every step and the solver made again (nothing pinned) */
	double area_mult, volume_mult;  /* the held area and volume raised to these multiples of the input's (0: kept) */
	unsigned int ramp;              /* a little each step over this many (the paper's program's growth schedules) */
	double edge;                    /* the remeshing's target edge (0: the input's mean edge) */
} rp_opts;

static double rp_mesh_volume(rp_mesh const* m)
{
	double V = 0;
	unsigned int f;
	for (f = 0; f < m->nf; f++)
	{
		double const *a = m->x + 3 * m->t[3 * f], *b = m->x + 3 * m->t[3 * f + 1], *c = m->x + 3 * m->t[3 * f + 2];
		V += (a[0] * (b[1] * c[2] - b[2] * c[1]) + a[1] * (b[2] * c[0] - b[0] * c[2]) + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6;
	}
	return V;
}

/* solve the mesh (p, sd's settings) on the caller's threads; snapshots
   (their own copies of the mesh) after each of steps[0..nsnap), snapped (if
   given) called after each; o (may be NULL) its remeshing and schedules */
static void rp_solve_mesh(rp_mesh const* mesh, double p, qaws_repulsive_mesh_solve_desc const* sd0, unsigned int const* steps, unsigned int nsnap,
	rp_mesh* snaps, qaws_repulsive_report* reps, void (*snapped)(void*, unsigned int), void* user, rp_opts const* o)
{
	unsigned int threads = rp_thread_count(), i, k, it = 0;
	double A0 = rp_mesh_area(mesh), V0 = rp_mesh_volume(mesh), L = 0;
	int remesh = o && o->remesh;
	qaws_batch_executor ex;
	qaws_repulsive_mesh_desc md;
	qaws_repulsive_mesh* rm = NULL;
	qaws_repulsive_mesh_solve_desc sd = *sd0;
	qaws_repulsive_mesh_solver* s = NULL;
	qaws_repulsive_report r;
	qaws_scalar e0 = 0;
	rp_mesh cur;
	qaws_scalar* x = NULL;
	memset(&r, 0, sizeof(r));
	if (!rp_copy_mesh(mesh, &cur))
		return;
	for (i = 0; i < 3 * mesh->nf; i++) L += rp_dist(mesh->x + 3 * mesh->t[i], mesh->x + 3 * mesh->t[3 * (i / 3) + (i + 1) % 3]) / (3 * mesh->nf);
	if (o && o->edge > 0)
		L = o->edge;
	ex.parallel_for = rp_parallel_for; ex.user = &threads;
	memset(&md, 0, sizeof(md));
	md.p = (qaws_scalar)p; md.executor = &ex;
	for (k = 0; k < nsnap; k++)
	{
		while (it < steps[k] && it < rp_cap && r.reason == QAWS_REPULSIVE_RUNNING)
		{
			if (!s)
			{
				/* (made again after each remeshing: the triangles, and their count, changed) */
				qaws_scalar* nx = (qaws_scalar*)realloc(x, sizeof(qaws_scalar) * 3 * cur.nv);
				if (!nx)
				{
					r.reason = QAWS_REPULSIVE_NUMERICAL_FAILURE;
					break;
				}
				x = nx;
				for (i = 0; i < 3 * cur.nv; i++) x[i] = (qaws_scalar)cur.x[i];
				md.vertex_count = cur.nv; md.triangles = cur.t; md.triangle_count = cur.nf;
				if (it)
					sd.initial_step = 2 * r.last_step;   /* (where the last solver would have tried next) */
				if (qaws_repulsive_mesh_create(&rm, &md) != QAWS_STATUS_OK || qaws_repulsive_mesh_solver_create(&s, rm, x, &sd) != QAWS_STATUS_OK)
				{
					printf("  the scene was refused\n");
					r.reason = QAWS_REPULSIVE_NUMERICAL_FAILURE;
					break;
				}
			}
			if (o && o->ramp)
			{
				/* the growth schedules: this step's targets */
				double f = (double)(it + 1 < o->ramp ? it + 1 : o->ramp) / o->ramp;
				qaws_scalar a = (qaws_scalar)(A0 * (1 + (o->area_mult - 1) * f)), v = (qaws_scalar)(V0 * (1 + (o->volume_mult - 1) * f));
				qaws_repulsive_mesh_solver_set_targets(s, sd.hold_area && o->area_mult > 0 ? &a : NULL, sd.hold_volume && o->volume_mult > 0 ? &v : NULL);
			}
			qaws_repulsive_mesh_solver_step(&r, x, s);
			if (it == 0)
				e0 = r.initial_energy;
			r.initial_energy = e0;
			r.iterations = ++it;
			for (i = 0; i < 3 * cur.nv; i++) cur.x[i] = (double)x[i];
			if (remesh && r.reason == QAWS_REPULSIVE_RUNNING)
			{
				qaws_repulsive_mesh_solver_destroy(s); s = NULL;
				qaws_repulsive_mesh_destroy(rm); rm = NULL;
				rp_remesh(&cur, L);
			}
			if (it % 25 == 0 && r.reason == QAWS_REPULSIVE_RUNNING)
				printf("    %u steps: E %.6g, stationarity %.3g\n", it, (double)r.energy, (double)r.stationarity);
		}
		reps[k] = r;
		if (!rp_copy_mesh(&cur, snaps + k))
			memset(snaps + k, 0, sizeof(rp_mesh));
		printf("  %u steps: E %.6g, stationarity %.3g, reason %d\n", it, (double)r.energy, (double)r.stationarity, (int)r.reason);
		fflush(stdout);
		if (snapped)
			snapped(user, k);
	}
	qaws_repulsive_mesh_solver_destroy(s);
	qaws_repulsive_mesh_destroy(rm);
	rp_free_mesh(&cur);
	free(x);
}

/* panels of the meshes ms[0..n), labelled, one scale for all; each face
   coloured by its component, or tint for all; ex (NULL, or one per panel)
   their points, see-through and cuts */
static void rp_mesh_panels(char const* svgname, char const* title, char const* subtitle, rp_mesh const* ms,
	char const (*labels)[200], unsigned int n, rp_view const* view, int const* tint, rp_extra const* ex)
{
	static int const rgb[3] = { 72, 132, 222 };   /* (the paper's blue) */
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
			rp_fit_points(&v, view, ms[k].x, ms[k].nv, &fit);
			scale = fmin(scale, fit.scale);   /* (one scale for every panel: sizes compare) */
		}
		for (k = 0; k < n; k++)
		{
			char png[600];
			char const* base;
			size_t stem = strlen(svgname) > 4 ? strlen(svgname) - 4 : strlen(svgname);   /* (without .svg) */
			rp_mesh const* m = ms + k;
			int* rgbs = rp_face_colours(m);
			unsigned int i;
			for (i = 0; tint && rgbs && i < 3 * m->nf; i++) rgbs[i] = tint[i % 3];
			for (i = 0; ex && ex[k].face_rgb && rgbs && i < 3 * m->nf; i++) rgbs[i] = ex[k].face_rgb[i];
			v.x0 = 24 + k * (pw + 12);
			svg_panel(&s, &v, NULL);
			rp_fit_points(&v, view, m->x, m->nv, &fit);
			fit.cx = v.x0 + 0.5 * v.w + (fit.cx - v.x0 - 0.5 * v.w) * scale / fit.scale;
			fit.cy = v.y0 + 0.5 * v.h + (fit.cy - v.y0 - 0.5 * v.h) * scale / fit.scale;
			fit.scale = scale;
			/* the panel rendered at twice its size, beside the SVG */
			sprintf(png, "%.*s_%u.png", (int)stem, svgname, k);
			base = strrchr(png, '/') ? strrchr(png, '/') + 1 : png;
			if (rp_render(png, (int)(2 * pw), (int)(2 * ph), v.x0, v.y0, 2, &fit, view, m->x, m->nv, m->t, m->nf, rgb, rgbs, ex ? ex + k : NULL))
				fprintf(s.f, "<image href=\"%s\" x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\"/>\n", base, v.x0, v.y0, pw, ph);
			fprintf(s.f, "<text x=\"%.1f\" y=\"%.1f\" font-size=\"13\" font-weight=\"600\" fill=\"#24292f\">%s</text>\n", v.x0 + 10, v.y0 + 20, labels[k]);
			free(rgbs);
		}
		svg_close(&s);
	}
}

/* the shrink wrap seen with what it wraps (from figure 15's snapshots
   showcase/paper15_shrinkwrap_2.obj and _4.obj): the bunny's points alone,
   the wrap see-through over them at 100 steps and at its end, and there cut open
   by the plane through the points' centre facing the eye */
static void rp_wrap_views(char const* sdir)
{
	static int const clay[3] = { 205, 196, 186 };
	rp_view view = { 0.3, -1.4 };   /* (the bunny side on: its data has z up) */
	rp_mesh ms[4], pts;
	rp_extra ex[4];
	char labels[4][200];
	double c[3] = { 0, 0, 0 }, eye[3] = { 0, 0, 0 };
	unsigned int i, k;
	memset(ms, 0, sizeof(ms)); memset(ex, 0, sizeof(ex));
	if (!rp_load_mesh(sdir, "bunny_points.obj", &pts) || !rp_load_mesh("showcase", "paper15_shrinkwrap_2.obj", ms + 1) ||
		!rp_load_mesh("showcase", "paper15_shrinkwrap_4.obj", ms + 2) || !rp_copy_mesh(ms + 2, ms + 3))
	{
		printf("paper15_shrinkwrap_views: run figure 15 first (its snapshots)\n");
		rp_free_mesh(&pts);
		for (k = 0; k < 4; k++) rp_free_mesh(ms + k);
		return;
	}
	ms[0] = pts; ms[0].nf = 0;   /* (the points alone: a mesh of no faces, drawn by its points) */
	for (i = 0; i < pts.nv; i++) for (k = 0; k < 3; k++) c[k] += pts.x[3 * i + k] / pts.nv;
	for (k = 0; k < 3; k++)
	{
		/* toward the eye in the mesh's frame: the rotation's transpose applied to (0, 0, -1) */
		double e[3] = { 0, 0, 0 }, re[3];
		e[k] = 1;
		rp_rotate(&view, e, re);
		eye[k] = -re[2];
	}
	for (k = 0; k < 4; k++)
	{
		ex[k].points = pts.x; ex[k].npoints = pts.nv; ex[k].radius = 1.1;
		ex[k].prgb[0] = 176; ex[k].prgb[1] = 64; ex[k].prgb[2] = 52;
	}
	ex[1].alpha = ex[2].alpha = 0.3;
	for (k = 0; k < 3; k++) ex[3].cut[k] = eye[k];
	ex[3].cut[3] = -(eye[0] * c[0] + eye[1] * c[1] + eye[2] * c[2]);
	sprintf(labels[0], "the bunny's %u points", pts.nv);
	sprintf(labels[1], "100 steps, see-through");
	sprintf(labels[2], "final, see-through");
	sprintf(labels[3], "final, cut open");
	rp_mesh_panels("showcase/paper15_shrinkwrap_views.svg", "Repulsive Surfaces, Fig. 20: the shrink wrap around the points",
		"the wrap drawn see-through (opacity 0.3 a layer), then cut by the plane through the points' centre facing the eye", ms, labels, 4, &view, clay, ex);
	for (k = 1; k < 4; k++) rp_free_mesh(ms + k);
	rp_free_mesh(&pts);
}

typedef struct rp_surface_run
{
	rp_mesh const* m;
	char const *svgname, *title, *subtitle;
	rp_view const* view;
	int const* tint;
	rp_mesh ms[9];   /* the input, then the snapshots */
	qaws_repulsive_report reps[8];
} rp_surface_run;

/* the figure so far, rewritten at each snapshot: a long run shows its
   progress, and a run stopped early keeps the panels it reached */
static void rp_surface_snapped(void* user, unsigned int last)
{
	rp_surface_run const* f = (rp_surface_run const*)user;
	char labels[9][200];
	unsigned int k;
	sprintf(labels[0], "input: %u faces, E %.4g", f->m->nf, (double)f->reps[0].initial_energy);
	for (k = 0; k <= last; k++)
		sprintf(labels[k + 1], "%s%u steps: E %.4g", f->reps[k].reason == QAWS_REPULSIVE_CONVERGED ? "converged, " : "", f->reps[k].iterations, (double)f->reps[k].energy);
	rp_mesh_panels(f->svgname, f->title, f->subtitle, f->ms, labels, last + 2, f->view, f->tint, NULL);
	{
		/* the snapshot as an OBJ beside the figure (to score it elsewhere, or to start from it) */
		char path[600];
		FILE* o;
		rp_mesh const* m = f->ms + last + 1;
		sprintf(path, "%.*s_%u.obj", (int)strlen(f->svgname) - 4, f->svgname, last + 1);
		if ((o = fopen(path, "w")) != NULL)
		{
			for (k = 0; k < m->nv; k++) fprintf(o, "v %.9g %.9g %.9g\n", m->x[3 * k], m->x[3 * k + 1], m->x[3 * k + 2]);
			for (k = 0; k < m->nf; k++) fprintf(o, "f %u %u %u\n", m->t[3 * k] + 1, m->t[3 * k + 1] + 1, m->t[3 * k + 2] + 1);
			fclose(o);
		}
	}
}

static void rp_surface_figure(rp_mesh const* m, char const* svgname, char const* title, char const* subtitle, double p,
	qaws_repulsive_mesh_solve_desc const* sd, unsigned int const* steps, unsigned int nsnap, rp_view const* view, int const* tint, rp_opts const* o)
{
	rp_surface_run f;
	time_t t0 = time(NULL);
	unsigned int k;
	f.m = m; f.svgname = svgname; f.title = title; f.subtitle = subtitle; f.view = view; f.tint = tint;
	f.ms[0] = *m;
	printf("%s: %u vertices, %u faces\n", svgname, m->nv, m->nf);
	rp_solve_mesh(m, p, sd, steps, nsnap, f.ms + 1, f.reps, rp_surface_snapped, &f, o);
	for (k = 0; k < nsnap; k++) rp_free_mesh(f.ms + 1 + k);
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

/* ------------------------------------------------------------------ */
/*  Vines: repulsive curves growing over a surface                    */
/* ------------------------------------------------------------------ */

enum { RP_VMAX = 900 };   /* the most vine vertices */

typedef struct rp_vines
{
	qaws_repulsive_surface_desc sd;
	qaws_scalar uv[2 * RP_VMAX];        /* each vertex's surface parameters */
	unsigned int e[2 * RP_VMAX];        /* edges, each from the vertex it grew from */
	int parent[RP_VMAX];                /* -1: a root */
	unsigned int age[RP_VMAX];          /* the stage it grew in */
	unsigned char grow[RP_VMAX];        /* a tip still growing */
	unsigned int nv, ne;
} rp_vines;

/* the surface point and its unit normal (outward from the z axis) at uv */
static void rp_vine_frame(rp_vines const* g, qaws_scalar const* uv, double* x, double* n)
{
	qaws_range rv = qaws_surface_get_v_range(g->sd.surface);
	double h = 1e-4 * ((double)rv.max_value - rv.min_value), dv = (double)uv[1] + h > (double)rv.max_value ? -h : h;
	qaws_scalar q[6] = { uv[0], uv[1], (qaws_scalar)((double)uv[0] + h), uv[1], uv[0], (qaws_scalar)((double)uv[1] + dv) }, p[9];
	double a[3], b[3], l;
	unsigned int k;
	qaws_repulsive_surface_positions(p, q, 3, &g->sd);
	for (k = 0; k < 3; k++) { x[k] = (double)p[k]; a[k] = (double)p[3 + k] - p[k]; b[k] = ((double)p[6 + k] - p[k]) * (dv > 0 ? 1 : -1); }
	n[0] = a[1] * b[2] - a[2] * b[1]; n[1] = a[2] * b[0] - a[0] * b[2]; n[2] = a[0] * b[1] - a[1] * b[0];
	l = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
	if (n[0] * x[0] + n[1] * x[1] < 0) l = -l;
	for (k = 0; k < 3; k++) n[k] /= l;
}

/* a new vertex grown from vertex t along the parameter direction d, one
   edge length L away in 3D (the direction rescaled by its image's length);
   -1 when it would pass the top of the surface's v range (the rim), 0 the
   bottom, or when there is no room */
static int rp_vine_add(rp_vines* g, unsigned int t, double du, double dv, double L, unsigned int stage)
{
	qaws_range rv = qaws_surface_get_v_range(g->sd.surface);
	qaws_scalar q[4];
	double x0[3], x1[3], n[3], s, len;
	unsigned int n2 = g->nv;
	if (n2 >= RP_VMAX)
		return 0;
	q[0] = g->uv[2 * t]; q[1] = g->uv[2 * t + 1];
	q[2] = (qaws_scalar)((double)q[0] + du); q[3] = (qaws_scalar)((double)q[1] + dv);
	rp_vine_frame(g, q, x0, n);
	rp_vine_frame(g, q + 2, x1, n);
	len = rp_dist(x0, x1);
	if (!(len > 0))
		return 0;
	s = L / len;
	q[2] = (qaws_scalar)((double)q[0] + s * du); q[3] = (qaws_scalar)((double)q[1] + s * dv);
	if ((double)q[3] > (double)rv.max_value - 0.02 * ((double)rv.max_value - rv.min_value))
		return -1;
	if ((double)q[3] < (double)rv.min_value)
		return 0;
	{
		/* no room: it would land on another vine (growth is not collision-checked; the descent then starts from a near contact) */
		unsigned int j;
		rp_vine_frame(g, q + 2, x1, n);
		for (j = 0; j < g->nv; j++)
		{
			double y[3];
			if (j == t)
				continue;
			rp_vine_frame(g, g->uv + 2 * j, y, n);
			if (rp_dist(x1, y) < 0.6 * L)
				return 0;
		}
	}
	g->uv[2 * n2] = q[2]; g->uv[2 * n2 + 1] = q[3];
	g->parent[n2] = (int)t; g->age[n2] = stage; g->grow[n2] = 1;
	g->e[2 * g->ne] = t; g->e[2 * g->ne + 1] = n2; g->ne++;
	g->nv++;
	return 1;
}

/* the vase and the vines (6-sided tubes, older segments thicker, sitting on
   the surface) and leaves (every third vertex, alternating sides) as one
   mesh, each face's colour in rgb */
static int rp_vine_mesh(rp_vines const* g, unsigned int stage, rp_mesh* m, int** rgb)
{
	enum { NU = 96, NV = 48, SIDES = 6 };
	static int const vase[3] = { 196, 122, 86 }, bark[3] = { 104, 76, 48 };
	qaws_range ru = qaws_surface_get_u_range(g->sd.surface), rv = qaws_surface_get_v_range(g->sd.surface);
	unsigned int maxv = NU * (NV + 1) + g->ne * 2 * SIDES + g->nv * 4, maxf = 2 * NU * NV + g->ne * 2 * SIDES + g->nv * 2;
	double* X = (double*)malloc(sizeof(double) * 3 * g->nv);
	double* N = (double*)malloc(sizeof(double) * 3 * g->nv);
	double* R = (double*)malloc(sizeof(double) * g->nv);
	unsigned int i, j, k, q, v = 0, f = 0, seed = 11;
	m->x = (double*)malloc(sizeof(double) * 3 * maxv);
	m->t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * maxf);
	*rgb = (int*)malloc(sizeof(int) * 3 * maxf);
	if (!X || !N || !R || !m->x || !m->t || !*rgb)
	{
		free(X); free(N); free(R);
		return 0;
	}
	/* the vase */
	for (j = 0; j <= NV; j++)
		for (i = 0; i < NU; i++)
		{
			qaws_scalar q2[2];
			double n[3];
			q2[0] = (qaws_scalar)((double)ru.min_value + ((double)ru.max_value - ru.min_value) * i / NU);
			q2[1] = (qaws_scalar)((double)rv.min_value + ((double)rv.max_value - rv.min_value) * j / NV);
			rp_vine_frame(g, q2, m->x + 3 * v, n);
			v++;
		}
	for (j = 0; j < NV; j++)
		for (i = 0; i < NU; i++)
		{
			unsigned int a = j * NU + i, b = j * NU + (i + 1) % NU, c = a + NU, d = b + NU;
			m->t[3 * f] = a; m->t[3 * f + 1] = b; m->t[3 * f + 2] = d; for (q = 0; q < 3; q++) (*rgb)[3 * f + q] = vase[q]; f++;
			m->t[3 * f] = a; m->t[3 * f + 1] = d; m->t[3 * f + 2] = c; for (q = 0; q < 3; q++) (*rgb)[3 * f + q] = vase[q]; f++;
		}
	/* the vines: radius by age, the tube's axis lifted off the surface by it */
	for (i = 0; i < g->nv; i++)
	{
		double frac = stage > 0 ? (double)(stage - (g->age[i] < stage ? g->age[i] : stage)) / stage : 0;
		R[i] = 0.004 + 0.011 * frac;
		rp_vine_frame(g, g->uv + 2 * i, X + 3 * i, N + 3 * i);
		for (k = 0; k < 3; k++) X[3 * i + k] += R[i] * N[3 * i + k];
	}
	for (i = 0; i < g->ne; i++)
	{
		unsigned int a = g->e[2 * i], b = g->e[2 * i + 1], base = v;
		double t[3], w[3], u2[3], l;
		for (k = 0; k < 3; k++) t[k] = X[3 * b + k] - X[3 * a + k];
		l = sqrt(t[0] * t[0] + t[1] * t[1] + t[2] * t[2]);
		if (!(l > 0))
			continue;
		for (k = 0; k < 3; k++) t[k] /= l;
		w[0] = t[1] * N[3 * a + 2] - t[2] * N[3 * a + 1]; w[1] = t[2] * N[3 * a] - t[0] * N[3 * a + 2]; w[2] = t[0] * N[3 * a + 1] - t[1] * N[3 * a];
		l = sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
		for (k = 0; k < 3; k++) w[k] /= l;
		u2[0] = t[1] * w[2] - t[2] * w[1]; u2[1] = t[2] * w[0] - t[0] * w[2]; u2[2] = t[0] * w[1] - t[1] * w[0];
		for (j = 0; j < 2 * SIDES; j++)
		{
			unsigned int end = j / SIDES;
			double ang = 6.283185307179586 * (j % SIDES) / SIDES, r = R[end ? b : a];
			for (k = 0; k < 3; k++) m->x[3 * v + k] = X[3 * (end ? b : a) + k] + r * (cos(ang) * w[k] + sin(ang) * u2[k]);
			v++;
		}
		for (j = 0; j < SIDES; j++)
		{
			unsigned int p0 = base + j, p1 = base + (j + 1) % SIDES, p2 = p0 + SIDES, p3 = p1 + SIDES;
			m->t[3 * f] = p0; m->t[3 * f + 1] = p1; m->t[3 * f + 2] = p3; for (q = 0; q < 3; q++) (*rgb)[3 * f + q] = bark[q]; f++;
			m->t[3 * f] = p0; m->t[3 * f + 1] = p3; m->t[3 * f + 2] = p2; for (q = 0; q < 3; q++) (*rgb)[3 * f + q] = bark[q]; f++;
		}
	}
	/* the leaves: a diamond in the surface's tangent plane, turned 65 degrees off the vine, alternating sides */
	for (i = 0; i < g->nv; i++)
	{
		double t[3], d[3], s[3], c[3], len = 0.07, wid = 0.035, sg = (i / 3) % 2 ? 1 : -1, l;
		int const* n = NULL;
		int green[3];
		if (g->parent[i] < 0 || i % 3 || g->age[i] + 2 > stage)
			continue;
		for (k = 0; k < 3; k++) t[k] = X[3 * i + k] - X[3 * g->parent[i] + k];
		l = sqrt(t[0] * t[0] + t[1] * t[1] + t[2] * t[2]);
		if (!(l > 0))
			continue;
		for (k = 0; k < 3; k++) t[k] /= l;
		s[0] = N[3 * i + 1] * t[2] - N[3 * i + 2] * t[1]; s[1] = N[3 * i + 2] * t[0] - N[3 * i] * t[2]; s[2] = N[3 * i] * t[1] - N[3 * i + 1] * t[0];
		for (k = 0; k < 3; k++) d[k] = cos(1.13) * t[k] + sg * sin(1.13) * s[k];
		s[0] = N[3 * i + 1] * d[2] - N[3 * i + 2] * d[1]; s[1] = N[3 * i + 2] * d[0] - N[3 * i] * d[2]; s[2] = N[3 * i] * d[1] - N[3 * i + 1] * d[0];
		for (k = 0; k < 3; k++) c[k] = X[3 * i + k] + 0.6 * R[i] * N[3 * i + k];
		for (k = 0; k < 3; k++)
		{
			m->x[3 * v + k] = c[k];
			m->x[3 * (v + 1) + k] = c[k] + 0.45 * len * d[k] + 0.5 * wid * s[k] + 0.012 * N[3 * i + k];
			m->x[3 * (v + 2) + k] = c[k] + len * d[k] + 0.004 * N[3 * i + k];
			m->x[3 * (v + 3) + k] = c[k] + 0.45 * len * d[k] - 0.5 * wid * s[k] + 0.012 * N[3 * i + k];
		}
		green[0] = 58 + (int)(30 * rp_rand(&seed)); green[1] = 118 + (int)(40 * rp_rand(&seed)); green[2] = 44 + (int)(20 * rp_rand(&seed));
		n = green;
		m->t[3 * f] = v; m->t[3 * f + 1] = v + 1; m->t[3 * f + 2] = v + 2; for (q = 0; q < 3; q++) (*rgb)[3 * f + q] = n[q]; f++;
		m->t[3 * f] = v; m->t[3 * f + 1] = v + 2; m->t[3 * f + 2] = v + 3; for (q = 0; q < 3; q++) (*rgb)[3 * f + q] = n[q]; f++;
		v += 4;
	}
	m->nv = v; m->nf = f;
	free(X); free(N); free(R);
	return 1;
}

/* vines climbing a vase: roots pinned around its foot, every stage a few
   descent steps (the total length held, collision-checked, the curves
   on the surface exactly) then each tip grown by an edge, and now and then
   a branch from three vertices back; the curves' repulsion spreads them */
static void rp_vines_figure(void)
{
	enum { ROOTS = 7, STAGES = 150, SNAPS = 4 };
	static unsigned int const snap_at[SNAPS] = { 20, 60, 100, 150 };
	static qaws_scalar const prof[14] = { (qaws_scalar)0.5, 0, (qaws_scalar)0.78, (qaws_scalar)0.2, (qaws_scalar)0.95, (qaws_scalar)0.62,
		(qaws_scalar)0.62, (qaws_scalar)1.15, (qaws_scalar)0.38, (qaws_scalar)1.5, (qaws_scalar)0.5, (qaws_scalar)1.82, (qaws_scalar)0.68, 2 };
	static qaws_scalar const knots[11] = { 0, 0, 0, 0, (qaws_scalar)0.25, (qaws_scalar)0.5, (qaws_scalar)0.75, 1, 1, 1, 1 };
	double const L = 0.045;
	qaws_curve* profile = NULL;
	qaws_surface* vase = NULL;
	qaws_bspline_desc bd;
	qaws_surface_revolution_desc rd;
	qaws_range ru, rv;
	rp_vines* g = (rp_vines*)calloc(1, sizeof(rp_vines));
	rp_mesh ms[SNAPS];
	rp_extra ex[SNAPS];
	int* rgbs[SNAPS];
	char labels[SNAPS][200];
	unsigned int pins[RP_VMAX], np = 0, s, k, i, ns = 0, seed = 5;   /* (the roots, then the vertices that reached the foot or the rim) */
	unsigned char held[RP_VMAX];
	rp_view view = { 0.35, -0.3 };
	time_t t0 = time(NULL);
	memset(ms, 0, sizeof(ms)); memset(ex, 0, sizeof(ex)); memset(rgbs, 0, sizeof(rgbs)); memset(held, 0, sizeof(held));
	memset(&bd, 0, sizeof(bd));
	bd.dimension = QAWS_DIMENSION_2D; bd.degree = 3; bd.control_points = prof; bd.control_point_count = 7; bd.knots = knots; bd.knot_count = 11;
	memset(&rd, 0, sizeof(rd));
	if (!g || qaws_curve_create_bspline(&profile, &bd) != QAWS_STATUS_OK)
	{
		free(g);
		return;
	}
	rd.profile = profile; rd.axis_direction.z = 1;
	if (qaws_surface_create_revolution(&vase, &rd) != QAWS_STATUS_OK)
	{
		qaws_curve_destroy(profile);
		free(g);
		return;
	}
	g->sd.surface = vase; g->sd.periodic_u = 1;
	ru = qaws_surface_get_u_range(vase); rv = qaws_surface_get_v_range(vase);
	printf("showcase/paper19_vines.svg: %u roots, %u stages\n", ROOTS, STAGES);
	/* the roots around the foot, each with a first shoot upward */
	for (k = 0; k < ROOTS; k++)
	{
		g->uv[2 * g->nv] = (qaws_scalar)((double)ru.min_value + ((double)ru.max_value - ru.min_value) * (k + 0.7 * rp_rand(&seed)) / ROOTS);
		g->uv[2 * g->nv + 1] = (qaws_scalar)((double)rv.min_value + 0.01 * ((double)rv.max_value - rv.min_value));
		g->parent[g->nv] = -1; g->age[g->nv] = 0; g->grow[g->nv] = 0;
		held[g->nv] = 1;
		pins[np++] = g->nv++;
		rp_vine_add(g, pins[np - 1], 0.02 * (rp_rand(&seed) - 0.5), 0.01, L, 0);
	}
	for (s = 1; s <= STAGES && ns < SNAPS; s++)
	{
		qaws_repulsive_curves_desc cd;
		qaws_repulsive_curves* c = NULL;
		qaws_repulsive_curves_solve_desc d;
		qaws_repulsive_report r;
		unsigned int n0 = g->nv;
		qaws_status st;
		/* a vine at the foot or the rim holds on there: a vertex at the parameter range's edge could not move outward at all
		   (the surface solver refuses any step leaving it), so every step would fail */
		for (i = 0; i < g->nv; i++)
		{
			double v = (double)g->uv[2 * i + 1], span = (double)rv.max_value - rv.min_value;
			if (!held[i] && (v < (double)rv.min_value + 0.015 * span || v > (double)rv.max_value - 0.015 * span))
			{
				held[i] = 1;
				pins[np++] = i;
			}
		}
		memset(&cd, 0, sizeof(cd));
		cd.dimension = 3; cd.vertex_count = g->nv; cd.edges = g->e; cd.edge_count = g->ne;
		memset(&d, 0, sizeof(d));
		d.hold_total_length = 1; d.collision_checked = 1; d.max_iterations = 3; d.pinned = pins; d.pinned_count = np;
		/* (the total length, one row: every edge's, hundreds of rows pulled back through the surface, stalled the projection once the vines filled out) */
		memset(&r, 0, sizeof(r));
		st = qaws_repulsive_curves_create(&c, &cd);
		st = st == QAWS_STATUS_OK ? qaws_repulsive_surface_curves_solve(&r, g->uv, c, &g->sd, &d) : st;
		qaws_repulsive_curves_destroy(c);
		if (st != QAWS_STATUS_OK)
		{
			printf("  stage %u: the solve failed (%d)\n", s, (int)st);
			break;
		}
		/* each tip grown one edge on: its direction kept, drawn a little upward, turned a little */
		for (i = 0; i < n0; i++)
		{
			double du, dv, a = 0.5 * (rp_rand(&seed) - 0.5), ca = cos(a), sa = sin(a), l;
			if (!g->grow[i])
				continue;
			g->grow[i] = 0;
			du = (double)g->uv[2 * i] - g->uv[2 * g->parent[i]];
			dv = (double)g->uv[2 * i + 1] - g->uv[2 * g->parent[i] + 1];
			l = sqrt(du * du + dv * dv);
			dv += 0.3 * l;
			if (rp_vine_add(g, i, ca * du - sa * dv, sa * du + ca * dv, L, s) < 0)
				held[i] = 1, pins[np++] = i;   /* (at the rim: it holds on there) */
		}
		/* now and then a branch from three vertices behind a new tip, turned 60 degrees off the vine there */
		if (s % 12 == 6)
		{
			unsigned int n1 = g->nv, t;
			for (t = n0; t < n1; t++)
			{
				int b = (int)t, up;
				unsigned int deg = 0, j;
				double du, dv, a = (rp_rand(&seed) < 0.5 ? -1 : 1) * 1.05;
				if (rp_rand(&seed) > 0.25)
					continue;
				for (k = 0; k < 3 && g->parent[b] >= 0; k++) b = g->parent[b];
				up = g->parent[b];
				for (j = 0; j < g->ne; j++) deg += g->e[2 * j] == (unsigned int)b || g->e[2 * j + 1] == (unsigned int)b;
				if (up < 0 || deg >= 3)
					continue;
				du = (double)g->uv[2 * b] - g->uv[2 * up]; dv = (double)g->uv[2 * b + 1] - g->uv[2 * up + 1];
				rp_vine_add(g, (unsigned int)b, cos(a) * du - sin(a) * dv, sin(a) * du + cos(a) * dv, L, s);
			}
		}
		if (s % 10 == 0)
			printf("    stage %u: %u vertices, E %.5g, reason %d, %u steps (refused: energy %u, collision %u, projection %u)\n", s, g->nv, (double)r.energy,
				(int)r.reason, r.iterations, r.rejected_energy, r.rejected_collision, r.rejected_projection);
		if (s == snap_at[ns])
		{
			if (rp_vine_mesh(g, s, ms + ns, rgbs + ns))
			{
				ex[ns].face_rgb = rgbs[ns];
				sprintf(labels[ns], "stage %u: %u vertices, %u edges", s, g->nv, g->ne);
				ns++;
				rp_mesh_panels("showcase/paper19_vines.svg", "Vines climbing a vase: repulsive curves growing on a surface",
					"7 roots pinned; each stage 3 steps (the total length held, collision-checked, on the surface exactly), then every tip grown and now and then a branch",
					ms, labels, ns, &view, NULL, ex);
			}
		}
	}
	for (k = 0; k < ns; k++) { rp_free_mesh(ms + k); free(rgbs[k]); }
	qaws_surface_destroy(vase);
	qaws_curve_destroy(profile);
	free(g);
	printf("  showcase/paper19_vines.svg, %.0f s\n", difftime(time(NULL), t0));
}

int main(int argc, char** argv)
{
	char const* dir = argc > 1 ? argv[1] : "../examples/data/repulsive_curves";
	int only = argc > 2 ? atoi(argv[2]) : 0;   /* (one figure: 1 to 19) */
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
		rp_view view = { 0.5, -0.32 };   /* (from slightly above, as the paper) */
		rp_mesh m;
		qaws_repulsive_mesh_solve_desc sd;
		unsigned int npin;
		unsigned int* pins;
		if (rp_cylinder(48, 36, 1, 4, &m) && (pins = rp_boundary(&m, &npin)) != NULL)   /* (twice as tall as wide, as the paper draws it) */
		{
			memset(&sd, 0, sizeof(sd));
			sd.pinned = pins; sd.pinned_count = npin; sd.theta = (qaws_scalar)0.25;
			rp_surface_figure(&m, "showcase/paper7_cylinder.svg", "Repulsive Surfaces, Fig. 9 (left): a cylinder with Dirichlet ends",
				"both rims pinned; p 6; the hierarchical fractional Sobolev descent (theta 0.25)", 6, &sd, steps, 2, &view, NULL, NULL);
			free(pins);
		}
		rp_free_mesh(&m);
	}
	if (!only || only == 10)
	{
		/* Repulsive Surfaces, Fig. 9 (centre): area minimization between two
		   rings past the catenoid's limit (3 apart for radius 1: a minimal
		   surface pinches off), regularized by the tangent-point energy */
		static unsigned int const steps[1] = { RP_UNTIL };
		static double const weights[3] = { 1, 10, 100 };
		rp_view view = { 0.5, -0.32 };   /* (from slightly above, as the paper) */
		rp_mesh m;
		qaws_repulsive_mesh_solve_desc sd;
		unsigned int npin, k;
		unsigned int* pins;
		time_t t0 = time(NULL);
		if (rp_cylinder(48, 28, 1, 3, &m) && (pins = rp_boundary(&m, &npin)) != NULL)
		{
			rp_mesh ms[4];   /* the input, then each weight's result */
			qaws_repulsive_report reps[3];
			char labels[4][200];
			memset(&sd, 0, sizeof(sd));
			sd.pinned = pins; sd.pinned_count = npin; sd.theta = (qaws_scalar)0.25;
			ms[0] = m;
			sprintf(labels[0], "input: %u faces", m.nf);
			for (k = 0; k < 3; k++)
			{
				sd.area_weight = (qaws_scalar)weights[k];
				printf("paper10_neck: the area weighted by %g\n", weights[k]);
				rp_solve_mesh(&m, 6, &sd, steps, 1, ms + 1 + k, reps + k, NULL, NULL, NULL);
				sprintf(labels[k + 1], "area weight %g: %u steps, E %.4g", weights[k], reps[k].iterations, (double)reps[k].energy);
			}
			rp_mesh_panels("showcase/paper10_neck.svg", "Repulsive Surfaces, Fig. 9 (centre): area minimization without pinching off",
				"two rings 3 apart (radius 1: far past the catenoid's limit of 1.33), pinned; the area weighted, the tangent-point energy (p 6) keeping a neck",
				ms, labels, 4, &view, NULL, NULL);
			for (k = 0; k < 3; k++) rp_free_mesh(ms + 1 + k);
			free(pins);
		}
		rp_free_mesh(&m);
		printf("  showcase/paper10_neck.svg, %.0f s\n", difftime(time(NULL), t0));
	}
	/* Repulsive Surfaces, scenes of the paper's code (not figures): linked surfaces flowing to
	   their canonical embeddings; the paper's meshes, the area and the
	   barycenter held as in the paper's scene files */
	{
		static char const* const files[2] = { "hopf.obj", "borromean.obj" };
		static char const* const svgs[2] = { "showcase/paper8_hopf.svg", "showcase/paper9_borromean.svg" };
		static char const* const titles[2] = { "Repulsive Surfaces, a scene of its code: a Hopf link", "Repulsive Surfaces, a scene of its code: Borromean rings" };
		/* shortcut: the rings stopped at the scene's 350 steps (9 s a step, every step backtracked to 1/32), far from the stationarity floor */
		static unsigned int const steps[2][2] = { { 25, RP_UNTIL }, { 175, 350 } };
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
				steps[k], 2, &view, NULL, NULL);
			rp_free_mesh(&m);
		}
	}
	/* Repulsive Surfaces: genus-2 scenes of the paper's code, its meshes and
	   settings (the area and the barycenter held, its iteration limits).
	   Fig. 1: linked handcuffs flowing to the unlinked double torus, in its
	   copper; a double torus with a ring through it; a figure-eight tube */
	{
		static char const* const files[3] = { "doubletorus10k.obj", "doubletorusring.obj", "figure8tube.obj" };
		static char const* const svgs[3] = { "showcase/paper11_handcuffs.svg", "showcase/paper12_torus_ring.svg", "showcase/paper13_figure8.svg" };
		static char const* const titles[3] = { "Repulsive Surfaces, Fig. 1: linked handcuffs to the unlinked double torus",
			"Repulsive Surfaces, a scene of its code: a double torus and a ring", "Repulsive Surfaces, a scene of its code: a broken figure-eight tube" };
		/* (the ring: past its scene's 200 iterations, where the paper's code
		   too is still threaded through both holes; one handle comes free
		   between 500 and 750 steps there) */
		static unsigned int const steps[3][4] = { { 100, 200, 350, 700 }, { 100, 300, 600, 1000 }, { 50, 100, 250, 500 } };
		static int const copper[3] = { 214, 128, 96 };
		static rp_opts const remeshed = { 1, 0, 0, 0, 0 };
		char sdir[1024];
		unsigned int k;
		sprintf(sdir, "%s/../repulsive_surfaces", dir);
		for (k = 0; k < 3; k++)
		{
			rp_view view = { 0.3, k == 0 ? 1.4 : -0.3 };   /* (the double torus from above: both holes seen) */
			rp_mesh m;
			qaws_repulsive_mesh_solve_desc sd;
			if (only && only != 11 + (int)k)
				continue;
			if (!rp_load_mesh(sdir, files[k], &m))
				continue;
			memset(&sd, 0, sizeof(sd));
			sd.hold_area = 1; sd.theta = (qaws_scalar)0.25;
			rp_surface_figure(&m, svgs[k], titles[k], "the paper's mesh; p 6; the area and the barycenter held; the hierarchical fractional Sobolev descent (theta 0.25)",
				6, &sd, steps[k], 4, &view, k == 0 ? copper : NULL, &remeshed);
			rp_free_mesh(&m);
		}
	}
	/* Repulsive Surfaces, Fig. 22 (top) and Fig. 20: growing a surface inside
	   an obstacle, and shrink wrapping a point cloud; the paper's scenes
	   (Walnut; Reconstruction with its bunny, the sphere at a quarter of the
	   faces), remeshed after every step as its program does */
	{
		char sdir[1024];
		rp_mesh m, ob;
		qaws_repulsive_mesh_solve_desc sd;
		qaws_repulsive_obstacle od;
		qaws_scalar* ox = NULL;
		unsigned int i;
		sprintf(sdir, "%s/../repulsive_surfaces", dir);
		memset(&ob, 0, sizeof(ob));
		memset(&od, 0, sizeof(od));
		if ((!only || only == 14) && rp_load_mesh(sdir, "walnut_inside.obj", &m))
		{
			/* the area grown 8x over 200 steps inside a sphere (weight 0.5), the volume kept softly (10) */
			static unsigned int const steps[4] = { 25, 75, 150, 220 };
			static rp_opts const grow = { 1, 8, 0, 200, 0 };
			static int const kernel[3] = { 200, 160, 110 };
			rp_view view = { 0.4, -0.35 };
			if (rp_load_mesh(sdir, "walnut_shell.obj", &ob) && (ox = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * ob.nv)) != NULL)
			{
				for (i = 0; i < 3 * ob.nv; i++) ox[i] = (qaws_scalar)ob.x[i];
				od.positions = ox; od.vertex_count = ob.nv; od.triangles = ob.t; od.triangle_count = ob.nf; od.weight = (qaws_scalar)0.5;
				memset(&sd, 0, sizeof(sd));
				sd.hold_area = 1; sd.theta = (qaws_scalar)0.25;
				sd.volume_weight = 10; sd.volume_target = (qaws_scalar)rp_mesh_volume(&m);   /* (absolute: the solver is made again each step) */
				sd.obstacles = &od; sd.obstacle_count = 1;
				rp_surface_figure(&m, "showcase/paper14_walnut.svg", "Repulsive Surfaces, Fig. 22 (top): growing a walnut inside its shell",
					"the paper's scene; p 6; the area grown 8x over 200 steps, the volume kept softly (10), a sphere for the shell (obstacle, 0.5); remeshed",
					6, &sd, steps, 4, &view, kernel, &grow);
			}
			free(ox); ox = NULL;
			rp_free_mesh(&ob); rp_free_mesh(&m);
		}
		if ((!only || only == 15) && rp_load_mesh(sdir, "sphere20k.obj", &m))
		{
			/* a sphere around the bunny's points, its volume brought down to 0.11x over 200 steps, its barycenter free */
			static unsigned int const steps[4] = { 50, 100, 150, 220 };
			static rp_opts const shrink = { 1, 0, 0.11, 200, 0 };
			static int const clay[3] = { 205, 196, 186 };
			rp_view view = { 0.3, -1.4 };   /* (the bunny side on: its data has z up) */
			for (i = 0; i < 3 * m.nv; i++) m.x[i] *= 0.6;   /* (the scene's sphere has radius 1.32) */
			if (rp_load_mesh(sdir, "bunny_points.obj", &ob) && (ox = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * ob.nv)) != NULL)
			{
				for (i = 0; i < 3 * ob.nv; i++) ox[i] = (qaws_scalar)ob.x[i];
				od.positions = ox; od.vertex_count = ob.nv; od.triangles = NULL; od.triangle_count = 0; od.weight = 1;
				memset(&sd, 0, sizeof(sd));
				sd.hold_volume = 1; sd.translations = 1; sd.theta = (qaws_scalar)0.25;
				sd.max_step = (qaws_scalar)0.4;   /* (as far as the paper's program moves: longer steps passed between the points) */
				sd.collision_checked = 1;   /* (no step through the points: their energy alone let some pass) */
				sd.obstacles = &od; sd.obstacle_count = 1;
				rp_surface_figure(&m, "showcase/paper15_shrinkwrap.svg", "Repulsive Surfaces, Fig. 20: shrink wrapping a point cloud",
					"the paper's scene: a sphere around the bunny's 14290 points (obstacle, 1), its volume to 0.11x over 200 steps, its barycenter free; p 6; remeshed",
					6, &sd, steps, 4, &view, clay, &shrink);
				rp_wrap_views(sdir);
			}
			free(ox); ox = NULL;
			rp_free_mesh(&ob); rp_free_mesh(&m);
		}
		if (!only || only == 19)
			rp_vines_figure();
		if (only == 18)
			rp_wrap_views(sdir);   /* (figure 15's views alone, from its snapshots) */
		if (!only || only == 16)
		{
			/* Fig. 21: nested envelopes by recursive shrink wrapping, around
			   the shrink-wrapped bunny (figure 15's last snapshot): each level
			   the scene's sphere shrunk onto the one inside it (an obstacle,
			   1), its volume brought down over 150 steps to 1.15x that one's,
			   remeshed coarser a level (2x, 2.5x, 3.1x the sphere's edge) */
			enum { LEVELS = 3 };
			static unsigned int const steps[1] = { 170 };
			rp_view view = { 0.3, -1.4 };
			rp_mesh lv[LEVELS + 2], sph;
			rp_extra ex[LEVELS + 2];
			char labels[LEVELS + 2][200];
			unsigned int nl, k, base;
			double L = 0, Vsph;
			memset(lv, 0, sizeof(lv)); memset(ex, 0, sizeof(ex)); memset(&sph, 0, sizeof(sph));
			if (rp_load_mesh("showcase", "paper15_shrinkwrap_4.obj", lv) && rp_load_mesh(sdir, "sphere20k.obj", &sph))
			{
				time_t t0 = time(NULL);
				for (i = 0; i < 3 * sph.nv; i++) sph.x[i] *= 0.6;   /* (as figure 15's) */
				for (i = 0; i < 3 * sph.nf; i++) L += rp_dist(sph.x + 3 * sph.t[i], sph.x + 3 * sph.t[3 * (i / 3) + (i + 1) % 3]) / (3 * sph.nf);
				Vsph = rp_mesh_volume(&sph);
				sprintf(labels[0], "level 0, the shrink wrap: %u faces", lv[0].nf);
				printf("showcase/paper16_envelopes.svg: around %u faces\n", lv[0].nf);
				for (nl = 1; nl <= LEVELS; nl++)
				{
					rp_mesh const* in = lv + nl - 1;
					qaws_repulsive_report rep;
					rp_opts o;
					memset(&o, 0, sizeof(o));
					o.remesh = 1; o.ramp = 150; o.edge = 2 * L * pow(1.25, nl - 1);
					o.volume_mult = 1.15 * rp_mesh_volume(in) / Vsph;
					ox = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * in->nv);
					if (!ox)
						break;
					for (i = 0; i < 3 * in->nv; i++) ox[i] = (qaws_scalar)in->x[i];
					od.positions = ox; od.vertex_count = in->nv; od.triangles = in->t; od.triangle_count = in->nf; od.weight = 1;
					memset(&sd, 0, sizeof(sd));
					sd.hold_volume = 1; sd.translations = 1; sd.theta = (qaws_scalar)0.25; sd.max_step = (qaws_scalar)0.4;
					sd.collision_checked = 1;   /* (strictly nested: no step through the level inside) */
					sd.obstacles = &od; sd.obstacle_count = 1;
					printf("  level %u: the volume to %.4g\n", nl, 1.15 * rp_mesh_volume(in));
					rp_solve_mesh(&sph, 6, &sd, steps, 1, lv + nl, &rep, NULL, NULL, &o);
					free(ox); ox = NULL;
					if (!lv[nl].nv)
						break;
					sprintf(labels[nl], "level %u: %.3gx the level inside, %u faces", nl, rp_mesh_volume(lv + nl) / rp_mesh_volume(in), lv[nl].nf);
				}
				/* all levels at once, cut open by the plane through their centre facing the eye (one colour each) */
				for (k = 0, base = 0; k < nl; k++) { lv[nl].nv += lv[k].nv; lv[nl].nf += lv[k].nf; }
				lv[nl].x = (double*)malloc(sizeof(double) * 3 * (lv[nl].nv + 1));
				lv[nl].t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (lv[nl].nf + 1));
				if (lv[nl].x && lv[nl].t)
				{
					double c[3] = { 0, 0, 0 }, eye[3];
					unsigned int f = 0;
					for (k = 0; k < nl; k++)
					{
						memcpy(lv[nl].x + 3 * base, lv[k].x, sizeof(double) * 3 * lv[k].nv);
						for (i = 0; i < 3 * lv[k].nf; i++) lv[nl].t[3 * f + i] = lv[k].t[i] + base;
						base += lv[k].nv; f += lv[k].nf;
					}
					for (i = 0; i < lv[0].nv; i++) for (k = 0; k < 3; k++) c[k] += lv[0].x[3 * i + k] / lv[0].nv;
					for (k = 0; k < 3; k++)
					{
						double e[3] = { 0, 0, 0 }, re[3];
						e[k] = 1;
						rp_rotate(&view, e, re);
						eye[k] = -re[2];
					}
					for (k = 0; k < 3; k++) ex[nl].cut[k] = eye[k];
					ex[nl].cut[3] = -(eye[0] * c[0] + eye[1] * c[1] + eye[2] * c[2]);
					sprintf(labels[nl], "all levels, cut open");
					rp_mesh_panels("showcase/paper16_envelopes.svg", "Repulsive Surfaces, Fig. 21: nested envelopes by recursive shrink wrapping",
						"each level the sphere shrunk onto the one inside it (obstacle, 1) to 1.15x its volume, remeshed coarser; p 6; the last panel cut open",
						lv, labels, nl + 1, &view, NULL, ex);
				}
				for (k = 0; k <= nl; k++) rp_free_mesh(lv + k);
				printf("  showcase/paper16_envelopes.svg, %.0f s\n", difftime(time(NULL), t0));
			}
			else
				printf("paper16_envelopes: run figure 15 first (the envelopes start from its shrink wrap)\n");
			rp_free_mesh(&sph);
		}
		if (!only || only == 17)
		{
			/* Fig. 22 (bottom): many small stones grown inside a slab, the
			   walnut's technique (its scene is not among the published ones):
			   24 flattened ellipsoids of different sizes, shapes and turns
			   scattered as by hand (a sphere already encloses the most for its
			   area, so its volume could not outgrow its area), their total
			   volume raised to 60 % of the slab's and their total area 1.15x
			   what that volume's growth would take on its own, over 200 steps:
			   each stone fills the gaps its neighbours leave; the slab an
			   obstacle (0.1, as the program's Hamburger) */
			static unsigned int const steps[4] = { 50, 100, 200, 250 };
			static int const stone[3] = { 168, 162, 152 };
			rp_view view = { 0.35, 1.15 };   /* (from above, tilted) */
			rp_opts grow;
			memset(&m, 0, sizeof(m));
			memset(&grow, 0, sizeof(grow));
			if (rp_cobbles(24, 4, 0.7, 7, &m, &ob) && (ox = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * ob.nv)) != NULL)
			{
				grow.remesh = 1; grow.ramp = 200;
				grow.volume_mult = 0.6 * 4 * 4 * 0.7 / rp_mesh_volume(&m);
				grow.area_mult = 1.15 * pow(grow.volume_mult, 2.0 / 3);
				for (i = 0; i < 3 * ob.nv; i++) ox[i] = (qaws_scalar)ob.x[i];
				od.positions = ox; od.vertex_count = ob.nv; od.triangles = ob.t; od.triangle_count = ob.nf; od.weight = (qaws_scalar)0.1;
				memset(&sd, 0, sizeof(sd));
				sd.hold_area = 1; sd.hold_volume = 1; sd.theta = (qaws_scalar)0.25;
				sd.obstacles = &od; sd.obstacle_count = 1;
				rp_surface_figure(&m, "showcase/paper17_cobbles.svg", "Repulsive Surfaces, Fig. 22 (bottom): cobblestones grown in a slab",
					"24 flattened stones scattered in a 4 x 4 x 0.7 slab (obstacle, 0.1); the volume raised to 60 % of the slab's over 200 steps; p 6; remeshed",
					6, &sd, steps, 4, &view, stone, &grow);
			}
			free(ox); ox = NULL;
			rp_free_mesh(&ob); rp_free_mesh(&m);
		}
	}
	return 0;
}
