/* Figure 3: Clipper2's own test records (tests/data/clipper2/Polygons.txt)
   through qaws_clip: subjects dashed blue, clips dashed red, the result
   filled, with the path count and area next to the ones Clipper2 stores. */

#ifndef QAWS_EXAMPLE_DATA_DIR
#define QAWS_EXAMPLE_DATA_DIR "../tests/data"
#endif

typedef struct c2s_rec
{
	int number, ct, fr;
	double area;
	long long count;
	unsigned int ns, nc;
	unsigned int sstart[512], cstart[512];
	double* pts;
	unsigned int npts, cap;
} c2s_rec;

static int c2s_push(c2s_rec* r, double x, double y)
{
	if (r->npts == r->cap)
	{
		unsigned int nc = r->cap ? 2 * r->cap : 1024;
		double* g = (double*)realloc(r->pts, sizeof(double) * 2 * nc);
		if (!g) return 0;
		r->pts = g;
		r->cap = nc;
	}
	r->pts[2 * r->npts] = x;
	r->pts[2 * r->npts + 1] = y;
	r->npts++;
	return 1;
}

/* the record with the given number */
static int c2s_load(char const* path, int number, c2s_rec* r)
{
	FILE* f = fopen(path, "r");
	static char line[1 << 16];
	int in = 0, section = 0;
	memset(r, 0, sizeof(*r));
	if (!f) return 0;
	while (fgets(line, sizeof(line), f))
	{
		if (!strncmp(line, "CAPTION:", 8))
		{
			if (in) break;
			in = atoi(line + 8) == number;
			r->number = number;
			continue;
		}
		if (!in) continue;
		if (!strncmp(line, "CLIPTYPE:", 9))
			r->ct = strstr(line, "INTERSECTION") ? QAWS_CLIP_INTERSECTION : strstr(line, "UNION") ? QAWS_CLIP_UNION :
				strstr(line, "DIFFERENCE") ? QAWS_CLIP_DIFFERENCE : QAWS_CLIP_XOR;
		else if (!strncmp(line, "FILLRULE:", 9))
			r->fr = strstr(line, "EVENODD") ? QAWS_FILL_EVEN_ODD : strstr(line, "NONZERO") ? QAWS_FILL_NON_ZERO :
				strstr(line, "POSITIVE") ? QAWS_FILL_POSITIVE : QAWS_FILL_NEGATIVE;
		else if (!strncmp(line, "SOL_AREA:", 9)) r->area = atof(line + 9);
		else if (!strncmp(line, "SOL_COUNT:", 10)) r->count = atoll(line + 10);
		else if (!strncmp(line, "SUBJECTS_OPEN", 13)) section = 0;
		else if (!strncmp(line, "SUBJECTS", 8)) section = 1;
		else if (!strncmp(line, "CLIPS", 5)) section = 2;
		else if (section && strpbrk(line, "0123456789"))
		{
			char* s = line;
			double v[2];
			int k = 0;
			if (section == 1 && r->ns < 511) r->sstart[r->ns++] = r->npts;
			if (section == 2 && r->nc < 511) r->cstart[r->nc++] = r->npts;
			while (*s)
			{
				char* e;
				double d;
				while (*s == ' ' || *s == ',' || *s == '\t' || *s == '\r' || *s == '\n') s++;
				if (!*s) break;
				d = strtod(s, &e);
				if (e == s) { s++; continue; }
				s = e;
				v[k++] = d;
				if (k == 2) { c2s_push(r, v[0], v[1]); k = 0; }
			}
		}
	}
	fclose(f);
	r->sstart[r->ns] = r->cstart[0];
	if (r->nc) r->cstart[r->nc] = r->npts;
	else r->sstart[r->ns] = r->npts;
	return r->number == number && (r->ns || r->nc);
}

static void c2s_poly(svg* s, viewport const* v, double const* p, unsigned int n, char const* color)
{
	double xy[2 * 2049];
	unsigned int i;
	if (n > 2048) n = 2048;
	for (i = 0; i < n; i++)
	{
		xy[2 * i] = vx(v, p[2 * i]);
		xy[2 * i + 1] = vy(v, p[2 * i + 1]);
	}
	xy[2 * n] = xy[0];
	xy[2 * n + 1] = xy[1];
	svg_polyline(s, xy, (int)n + 1, color, 1.0, 0.7, 1);
}

static void demo_clipper2(void)
{
	static int const picks[8] = { 2, 14, 18, 62, 101, 120, 160, 193 };
	static char const* ctn[5] = { "none", "intersection", "union", "difference", "xor" };
	static char const* frn[4] = { "even-odd", "non-zero", "positive", "negative" };
	int const cols = 4, rows = 2, pw = 300, ph = 300, gap = 16, top = 84;
	int const W = cols * pw + (cols + 1) * gap, H = top + rows * (ph + gap) + 8;
	char const* path = QAWS_EXAMPLE_DATA_DIR "/clipper2/Polygons.txt";
	svg s;
	int k;
	double total_time = 0.0;
	if (!svg_open(&s, "showcase/b2d3_clipper2.svg", W, H, "Clipper2's test records through qaws_clip",
		"subjects dashed blue, clips dashed red, the result filled; our path count and area against the ones Clipper2 stores (it rounds every vertex to integers)"))
		return;
	for (k = 0; k < 8; k++)
	{
		c2s_rec r;
		viewport v;
		qaws_curve** cs;
		qaws_curve const** views;
		qaws_path_2d* paths;
		qaws_clip_desc d;
		qaws_clip_result* res = NULL;
		unsigned int i, n;
		double lo[2] = { 1e300, 1e300 }, hi[2] = { -1e300, -1e300 }, cx, cy, half, area = 0, t0;
		char label[256];
		if (!c2s_load(path, picks[k], &r))
			continue;
		for (i = 0; i < r.npts; i++)
		{
			if (r.pts[2 * i] < lo[0]) lo[0] = r.pts[2 * i];
			if (r.pts[2 * i] > hi[0]) hi[0] = r.pts[2 * i];
			if (r.pts[2 * i + 1] < lo[1]) lo[1] = r.pts[2 * i + 1];
			if (r.pts[2 * i + 1] > hi[1]) hi[1] = r.pts[2 * i + 1];
		}
		cx = 0.5 * (lo[0] + hi[0]); cy = 0.5 * (lo[1] + hi[1]);
		half = 0.56 * (hi[0] - lo[0] > hi[1] - lo[1] ? hi[0] - lo[0] : hi[1] - lo[1]);
		v.x0 = gap + (k % cols) * (pw + gap); v.y0 = top + (k / cols) * (ph + gap); v.w = pw; v.h = ph;
		v.xmin = cx - half; v.xmax = cx + half; v.ymin = cy - half * 1.08; v.ymax = cy + half * 0.92;
		n = r.ns + r.nc;
		cs = (qaws_curve**)calloc(n + 1, sizeof(qaws_curve*));
		views = (qaws_curve const**)calloc(n + 1, sizeof(qaws_curve*));
		paths = (qaws_path_2d*)calloc(n + 1, sizeof(qaws_path_2d));
		for (i = 0; i < n; i++)
		{
			unsigned int a = i < r.ns ? r.sstart[i] : r.cstart[i - r.ns];
			unsigned int b = i < r.ns ? r.sstart[i + 1] : r.cstart[i - r.ns + 1];
			qaws_scalar* p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * (b - a + 1));
			unsigned int j;
			for (j = 0; j < 2 * (b - a); j++) p[j] = (qaws_scalar)r.pts[2 * a + j];
			if (b - a >= 3) qaws_curve_create_polyline_2d(p, b - a, 1, &cs[i]);
			free(p);
			views[i] = cs[i];
			paths[i].curves = &views[i];
			paths[i].curve_count = cs[i] ? 1 : 0;
			paths[i].closed = 1;
		}
		memset(&d, 0, sizeof(d));
		d.subjects = paths; d.subject_count = r.ns;
		d.clips = paths + r.ns; d.clip_count = r.nc;
		d.clip_type = (qaws_clip_type)r.ct; d.fill_rule = (qaws_fill_rule)r.fr;
		t0 = (double)clock() / CLOCKS_PER_SEC;
		qaws_clip_execute(&d, &res);
		total_time += (double)clock() / CLOCKS_PER_SEC - t0;
		sprintf(label, "#%d %s, %s", r.number, ctn[r.ct], frn[r.fr]);
		svg_panel(&s, &v, label);
		if (res)
		{
			bo_fill(&s, &v, res, "#54aeff");
			for (i = 0; i < qaws_clip_result_get_path_count(res); i++)
			{
				qaws_path_2d p;
				qaws_scalar a = 0;
				unsigned int c;
				qaws_clip_result_get_path(res, i, &p);
				qaws_path_compute_area_2d(&p, &a);
				area += a;
				for (c = 0; c < p.curve_count; c++)
					bo_curve(&s, &v, p.curves[c], "#0969da", 1.4, 0);
			}
		}
		for (i = 0; i < n; i++)
		{
			unsigned int a = i < r.ns ? r.sstart[i] : r.cstart[i - r.ns];
			unsigned int b = i < r.ns ? r.sstart[i + 1] : r.cstart[i - r.ns + 1];
			c2s_poly(&s, &v, &r.pts[2 * a], b - a, i < r.ns ? "#0550ae" : "#cf222e");
		}
		sprintf(label, "%u paths, area %.0f (Clipper2: %lld, %.0f)", res ? qaws_clip_result_get_path_count(res) : 0, area, r.count, r.area);
		svg_text(&s, v.x0 + 10, v.y0 + ph - 12, 12, "#57606a", "start", label);
		qaws_clip_result_destroy(res);
		for (i = 0; i < n; i++) qaws_curve_destroy(cs[i]);
		free(cs); free((void*)views); free(paths); free(r.pts);
	}
	svg_close(&s);
	printf("    clipper2 records: 8 in %.3f s\n", total_time);
}
