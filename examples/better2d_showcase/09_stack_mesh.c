/* Figure 9: stacks of contours as 3D surface meshes (qaws_stack_mesh_3d),
   drawn with flat shading: a cone from a circle and a point, a sphere from
   9 circles (cubic, sections between), a tube whose hole drifts, an
   ellipse splitting in two, cone + sphere meshed level by level, and a
   twisting star. */

typedef struct sm3_view
{
	viewport v;
	double az, el, scale, cx, cy, cz;
} sm3_view;

static void sm3_rot(sm3_view const* w, double const* p, double* q)
{
	double x = p[0] - w->cx, y = p[1] - w->cy, z = p[2] - w->cz;
	double u = x * cos(w->az) - y * sin(w->az), d = x * sin(w->az) + y * cos(w->az);
	q[0] = u;
	q[1] = z * cos(w->el) + d * sin(w->el);     /* up on screen */
	q[2] = d * cos(w->el) - z * sin(w->el);     /* away from the viewer */
}

/* a z-buffer rendering of the mesh (2x supersampled, flat shading), into
   the panel as an embedded image */
static void sm3_b64(FILE* f, unsigned char const* d, size_t n)
{
	static char const tab[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
	size_t i;
	for (i = 0; i + 2 < n; i += 3)
	{
		unsigned int v = ((unsigned int)d[i] << 16) | ((unsigned int)d[i + 1] << 8) | d[i + 2];
		fputc(tab[(v >> 18) & 63], f); fputc(tab[(v >> 12) & 63], f); fputc(tab[(v >> 6) & 63], f); fputc(tab[v & 63], f);
	}
	if (n - i == 1)
	{
		unsigned int v = (unsigned int)d[i] << 16;
		fputc(tab[(v >> 18) & 63], f); fputc(tab[(v >> 12) & 63], f); fputs("==", f);
	}
	else if (n - i == 2)
	{
		unsigned int v = ((unsigned int)d[i] << 16) | ((unsigned int)d[i + 1] << 8);
		fputc(tab[(v >> 18) & 63], f); fputc(tab[(v >> 12) & 63], f); fputc(tab[(v >> 6) & 63], f); fputc('=', f);
	}
}

static void sm3_draw(svg* s, sm3_view const* w, qaws_mesh_3d const* m, double r0, double g0, double b0)
{
	int const SS = 2, W = (int)w->v.w - 4, H = (int)w->v.h - 44, BW = W * SS, BH = H * SS;
	unsigned int nt = qaws_mesh_3d_get_triangle_count(m), i, k;
	float* zb = (float*)malloc(sizeof(float) * (size_t)BW * BH);
	unsigned char* rgb = (unsigned char*)malloc((size_t)BW * BH * 3);
	unsigned char* bmp;
	double L[3] = { -0.45, -0.55, 0.7 }, ll = sqrt(L[0] * L[0] + L[1] * L[1] + L[2] * L[2]);
	int x, y, row;
	size_t size;
	if (!zb || !rgb) { free(zb); free(rgb); return; }
	L[0] /= ll; L[1] /= ll; L[2] /= ll;
	for (i = 0; i < (unsigned int)(BW * BH); i++) { zb[i] = 1e30f; rgb[3 * i] = 0xf6; rgb[3 * i + 1] = 0xf8; rgb[3 * i + 2] = 0xfa; }
	for (i = 0; i < nt; i++)
	{
		unsigned int t[3];
		double q[3][3], sx[3], sy[3], nrm[3], nl, sh, area;
		unsigned char c[3];
		int x0, x1, y0, y1;
		qaws_mesh_3d_get_triangle(m, i, t);
		for (k = 0; k < 3; k++)
		{
			qaws_vec3 v = qaws_mesh_3d_get_vertex(m, t[k]);
			double p[3] = { v.x, v.y, v.z };
			sm3_rot(w, p, q[k]);
			sx[k] = SS * (W / 2.0 + w->scale * q[k][0]);
			sy[k] = SS * (H / 2.0 + 6 - w->scale * q[k][1]);
		}
		/* counter-clockwise on the screen (y down here): facing the viewer */
		area = (sx[1] - sx[0]) * (sy[2] - sy[0]) - (sy[1] - sy[0]) * (sx[2] - sx[0]);
		if (!(area < 0)) continue;
		{
			qaws_vec3 a = qaws_mesh_3d_get_vertex(m, t[0]), b = qaws_mesh_3d_get_vertex(m, t[1]), cc = qaws_mesh_3d_get_vertex(m, t[2]);
			double u[3] = { (double)b.x - a.x, (double)b.y - a.y, (double)b.z - a.z }, v[3] = { (double)cc.x - a.x, (double)cc.y - a.y, (double)cc.z - a.z };
			nrm[0] = u[1] * v[2] - u[2] * v[1];
			nrm[1] = u[2] * v[0] - u[0] * v[2];
			nrm[2] = u[0] * v[1] - u[1] * v[0];
		}
		nl = sqrt(nrm[0] * nrm[0] + nrm[1] * nrm[1] + nrm[2] * nrm[2]);
		if (!(nl > 0)) continue;
		sh = 0.32 + 0.68 * fmax(0, (nrm[0] * L[0] + nrm[1] * L[1] + nrm[2] * L[2]) / nl);
		c[0] = (unsigned char)(255 * fmin(1, r0 * sh)); c[1] = (unsigned char)(255 * fmin(1, g0 * sh)); c[2] = (unsigned char)(255 * fmin(1, b0 * sh));
		x0 = (int)floor(fmin(sx[0], fmin(sx[1], sx[2]))); x1 = (int)ceil(fmax(sx[0], fmax(sx[1], sx[2])));
		y0 = (int)floor(fmin(sy[0], fmin(sy[1], sy[2]))); y1 = (int)ceil(fmax(sy[0], fmax(sy[1], sy[2])));
		if (x0 < 0) x0 = 0;
		if (y0 < 0) y0 = 0;
		if (x1 > BW - 1) x1 = BW - 1;
		if (y1 > BH - 1) y1 = BH - 1;
		for (y = y0; y <= y1; y++)
			for (x = x0; x <= x1; x++)
			{
				double px = x + 0.5, py = y + 0.5;
				double w0 = ((sx[1] - px) * (sy[2] - py) - (sy[1] - py) * (sx[2] - px)) / area;
				double w1 = ((sx[2] - px) * (sy[0] - py) - (sy[2] - py) * (sx[0] - px)) / area;
				double w2 = 1 - w0 - w1;
				float z;
				if (w0 < -1e-9 || w1 < -1e-9 || w2 < -1e-9) continue;
				z = (float)(w0 * q[0][2] + w1 * q[1][2] + w2 * q[2][2]);
				if (z < zb[y * BW + x])
				{
					zb[y * BW + x] = z;
					rgb[3 * (y * BW + x)] = c[0]; rgb[3 * (y * BW + x) + 1] = c[1]; rgb[3 * (y * BW + x) + 2] = c[2];
				}
			}
	}
	/* down to the panel's pixels, as a 24-bit BMP (rows bottom up, padded to 4) */
	row = (3 * W + 3) & ~3;
	size = 54 + (size_t)row * H;
	bmp = (unsigned char*)calloc(size, 1);
	if (bmp)
	{
		unsigned int fs = (unsigned int)size;
		bmp[0] = 'B'; bmp[1] = 'M';
		bmp[2] = (unsigned char)fs; bmp[3] = (unsigned char)(fs >> 8); bmp[4] = (unsigned char)(fs >> 16); bmp[5] = (unsigned char)(fs >> 24);
		bmp[10] = 54; bmp[14] = 40;
		bmp[18] = (unsigned char)W; bmp[19] = (unsigned char)(W >> 8);
		bmp[22] = (unsigned char)H; bmp[23] = (unsigned char)(H >> 8);
		bmp[26] = 1; bmp[28] = 24;
		for (y = 0; y < H; y++)
			for (x = 0; x < W; x++)
			{
				int acc[3] = { 0, 0, 0 }, a, b;
				unsigned char* o = bmp + 54 + (size_t)(H - 1 - y) * row + 3 * x;
				for (b = 0; b < SS; b++)
					for (a = 0; a < SS; a++)
					{
						unsigned char const* p = rgb + 3 * ((size_t)(y * SS + b) * BW + x * SS + a);
						acc[0] += p[0]; acc[1] += p[1]; acc[2] += p[2];
					}
				o[0] = (unsigned char)(acc[2] / (SS * SS));
				o[1] = (unsigned char)(acc[1] / (SS * SS));
				o[2] = (unsigned char)(acc[0] / (SS * SS));
			}
		fprintf(s->f, "<image x=\"%.1f\" y=\"%.1f\" width=\"%d\" height=\"%d\" href=\"data:image/bmp;base64,", w->v.x0 + 2, w->v.y0 + 26, W, H);
		sm3_b64(s->f, bmp, size);
		fprintf(s->f, "\"/>\n");
		free(bmp);
	}
	free(zb);
	free(rgb);
}

typedef struct sm3_stack
{
	qaws_curve* c[160];
	qaws_curve const* v[160];
	qaws_path_2d p[160];
	qaws_stack_level l[64];
	unsigned int nc, nl;
} sm3_stack;

static void sm3_level(sm3_stack* s, double z, qaws_curve* const* cs, unsigned int n)
{
	unsigned int i, first = s->nc;
	for (i = 0; i < n; i++)
	{
		s->c[s->nc] = cs[i]; s->v[s->nc] = cs[i];
		s->p[s->nc].curves = &s->v[s->nc]; s->p[s->nc].curve_count = 1; s->p[s->nc].closed = 1;
		s->nc++;
	}
	s->l[s->nl].z = (qaws_scalar)z;
	s->l[s->nl].paths = &s->p[first];
	s->l[s->nl].path_count = n;
	s->nl++;
}

static void sm3_free(sm3_stack* s)
{
	unsigned int i;
	for (i = 0; i < s->nc; i++) qaws_curve_destroy(s->c[i]);
}

static qaws_curve* sm3_circle(double x, double y, double r, int hole)
{
	qaws_curve *c = sk_circ(x, y, r), *rv = NULL;
	if (!hole) return c;
	qaws_curve_reverse(c, &rv);
	qaws_curve_destroy(c);
	return rv;
}

static void demo_stack_mesh(void)
{
	int const cols = 3, rows = 2, pw = 420, ph = 380, gap = 16, top = 84;
	int const W = cols * pw + (cols + 1) * gap, H = top + rows * (ph + gap) + 8;
	svg s;
	int k;
	double t0;
	char label[200];
	if (!svg_open(&s, "showcase/b2d9_stack_mesh.svg", W, H, "Stacks of contours as surface meshes",
		"qaws_stack_mesh_3d: sections sampled on their curves, contours paired up the stack and stitched; fans at tips, caps at the ends and where contours split; closed meshes"))
		return;
	t0 = (double)clock() / CLOCKS_PER_SEC;
	for (k = 0; k < 6; k++)
	{
		sm3_view w;
		sm3_stack st;
		qaws_stack_2d sd;
		qaws_stack_mesh_desc d;
		qaws_mesh_3d* m = NULL;
		qaws_stack_result* br = NULL;
		qaws_curve* cs[3];
		unsigned int i, sections = 0;
		double col[3] = { 0.45, 0.65, 0.95 };
		char const* title = "";
		memset(&st, 0, sizeof(st));
		memset(&sd, 0, sizeof(sd));
		memset(&w, 0, sizeof(w));
		w.v.x0 = gap + (k % cols) * (pw + gap); w.v.y0 = top + (k / cols) * (ph + gap); w.v.w = pw; w.v.h = ph;
		w.az = 0.55; w.el = 0.42; w.scale = 62; w.cz = 1.1;
		sd.fill_rule = QAWS_FILL_NON_ZERO;
		if (k == 0)
		{
			cs[0] = sm3_circle(0, 0, 1.7, 0); sm3_level(&st, 0, cs, 1);
			cs[0] = sm3_circle(0, 0, 0, 0); sm3_level(&st, 2.6, cs, 1);
			title = "a cone: a circle and a point";
			w.cz = 1.3;
		}
		else if (k == 1)
		{
			for (i = 0; i <= 8; i++)
			{
				double th = CO_PI * (8 - i) / 8;
				cs[0] = sm3_circle(0, 0, 1.6 * sin(th), 0);
				sm3_level(&st, 1.6 + 1.6 * cos(th), cs, 1);
			}
			sd.interp = QAWS_STACK_CUBIC;
			sections = 4;
			title = "a sphere from 9 circles, cubic, 4 sections between";
			w.cz = 1.6; col[0] = 0.4; col[1] = 0.8; col[2] = 0.5;
		}
		else if (k == 2)
		{
			for (i = 0; i <= 6; i++)
			{
				double z = 0.45 * i, a = 0.9 * i;
				cs[0] = sm3_circle(0, 0, 1.8 - 0.08 * i, 0);
				cs[1] = sm3_circle(0.7 * cos(a), 0.7 * sin(a), 0.6, 1);
				sm3_level(&st, z, cs, 2);
			}
			sd.interp = QAWS_STACK_CUBIC;
			sections = 3;
			title = "a tube whose hole drifts round";
			w.cz = 1.35; col[0] = 0.95; col[1] = 0.7; col[2] = 0.35;
		}
		else if (k == 3)
		{
			qaws_curve* e = NULL;
			qaws_curve_create_ellipse_2d((qaws_scalar)0.2, 0, (qaws_scalar)2.2, (qaws_scalar)1.0, 0, &e);
			cs[0] = e; sm3_level(&st, 0, cs, 1);
			cs[0] = sm3_circle(-1.2, 0, 0.7, 0); cs[1] = sm3_circle(1.6, 0, 0.9, 0); sm3_level(&st, 2.2, cs, 2);
			sd.blend = QAWS_STACK_BLEND_DISTANCE; sd.grid = 120;
			sections = 10;
			title = "an ellipse splitting into two circles";
			w.cz = 1.1; col[0] = 0.75; col[1] = 0.55; col[2] = 0.95;
		}
		else if (k == 4)
		{
			/* cone + sphere, the Boolean level by level, then meshed */
			sm3_stack a, b;
			qaws_stack_2d sa, sb;
			memset(&a, 0, sizeof(a)); memset(&b, 0, sizeof(b));
			memset(&sa, 0, sizeof(sa)); memset(&sb, 0, sizeof(sb));
			cs[0] = sm3_circle(0, 0, 1.5, 0); sm3_level(&a, 0, cs, 1);
			cs[0] = sm3_circle(0, 0, 0, 0); sm3_level(&a, 2.6, cs, 1);
			for (i = 0; i <= 12; i++)
			{
				double th = CO_PI * (12 - i) / 12;
				cs[0] = sm3_circle(0.9, 0, 1.0 * sin(th), 0);
				sm3_level(&b, 1.3 + 1.0 * cos(th), cs, 1);
			}
			sa.levels = a.l; sa.level_count = a.nl; sa.fill_rule = QAWS_FILL_NON_ZERO;
			sb.levels = b.l; sb.level_count = b.nl; sb.fill_rule = QAWS_FILL_NON_ZERO; sb.interp = QAWS_STACK_CUBIC;
			if (qaws_stack_boolean_2d(QAWS_CLIP_UNION, &sa, &sb, NULL, &br) == QAWS_STATUS_OK)
			{
				/* the result's levels as a stack (matched by distance) */
				static qaws_path_2d rp[64][8];
				for (i = 0; i < qaws_stack_result_get_level_count(br) && i < 64; i++)
				{
					qaws_clip_result const* lr = qaws_stack_result_get_level(br, i);
					unsigned int j, np = qaws_clip_result_get_path_count(lr);
					for (j = 0; j < np && j < 8; j++) qaws_clip_result_get_path(lr, j, &rp[i][j]);
					st.l[i].z = qaws_stack_result_get_z(br, i);
					st.l[i].paths = rp[i];
					st.l[i].path_count = np < 8 ? np : 8;
				}
				st.nl = i;
			}
			sm3_free(&a); sm3_free(&b);
			sd.blend = QAWS_STACK_BLEND_DISTANCE; sd.grid = 120;
			sections = 2;
			title = "cone + sphere: the Boolean level by level, meshed";
			w.cz = 1.3; col[0] = 0.95; col[1] = 0.5; col[2] = 0.45;
		}
		else
		{
			/* a star twisting as it rises: matched contours of one structure */
			for (i = 0; i <= 6; i++)
			{
				double z = 0.42 * i, sc = 1.0 - 0.05 * i, a = 0.18 * i;
				qaws_curve* star = co_star(10, 1.9 * sc, 1.0 * sc, 0, 0);
				qaws_curve* rot = NULL;
				qaws_scalar mtx[6];
				mtx[0] = (qaws_scalar)cos(a); mtx[1] = (qaws_scalar)-sin(a); mtx[2] = 0;
				mtx[3] = (qaws_scalar)sin(a); mtx[4] = (qaws_scalar)cos(a); mtx[5] = 0;
				qaws_curve_transform_2d(star, mtx, &rot);
				qaws_curve_destroy(star);
				cs[0] = rot;
				sm3_level(&st, z, cs, 1);
			}
			sd.interp = QAWS_STACK_CUBIC;
			sections = 3;
			title = "a star twisting as it rises (matched contours)";
			w.cz = 1.25; col[0] = 0.45; col[1] = 0.85; col[2] = 0.9;
		}
		sd.levels = st.l; sd.level_count = st.nl;
		memset(&d, 0, sizeof(d));
		d.stack = &sd; d.sections = sections; d.caps = 1;
		svg_panel(&s, &w.v, title);
		if (qaws_stack_mesh_3d(&d, &m) == QAWS_STATUS_OK)
			sm3_draw(&s, &w, m, col[0], col[1], col[2]);
		sprintf(label, "%u vertices, %u triangles, %u levels%s", m ? qaws_mesh_3d_get_vertex_count(m) : 0, m ? qaws_mesh_3d_get_triangle_count(m) : 0,
			st.nl, sections ? ", sections between" : "");
		svg_text(&s, w.v.x0 + 10, w.v.y0 + ph - 12, 11, "#57606a", "start", label);
		qaws_mesh_3d_destroy(m);
		qaws_stack_result_destroy(br);
		sm3_free(&st);
	}
	svg_close(&s);
	printf("    stack mesh: 6 panels in %.3f s\n", (double)clock() / CLOCKS_PER_SEC - t0);
}
