/* Figure 5: 2.5D stacks of parallel contours (qaws_stack), drawn in an
   oblique view: the input levels dark, sections interpolated between them
   light. A cone from a circle and a point, a sphere from 9 circles blended
   linearly and cubically, a blob splitting in two (distance blend), and
   cone x sphere solved level by level. */

typedef struct sk_view
{
	viewport v;
	double sx, sy;    /* screen offset per unit of z */
} sk_view;

static void sk_draw_region(svg* s, sk_view const* sv, qaws_clip_result const* r, double z, char const* color, double width, double op)
{
	unsigned int i, k, j;
	for (i = 0; r && i < qaws_clip_result_get_path_count(r); i++)
	{
		qaws_path_2d p;
		qaws_clip_result_get_path(r, i, &p);
		for (k = 0; k < p.curve_count; k++)
		{
			double xy[2 * 97];
			qaws_range rg = qaws_curve_get_parameter_range(p.curves[k]);
			for (j = 0; j <= 96; j++)
			{
				qaws_eval_result_2d e;
				qaws_curve_evaluate_2d(p.curves[k], (qaws_scalar)(rg.min_value + (rg.max_value - rg.min_value) * j / 96), QAWS_EVAL_FLAG_POSITION, &e);
				xy[2 * j] = vx(&sv->v, e.position.x) + sv->sx * z;
				xy[2 * j + 1] = vy(&sv->v, e.position.y * 0.45) - sv->sy * z;
			}
			svg_polyline(s, xy, 97, color, width, op, 0);
		}
	}
}

static void sk_draw_stack(svg* s, sk_view const* sv, qaws_stack_2d const* st, unsigned int sections, char const* scol)
{
	unsigned int i;
	double z0 = st->levels[0].z, z1 = st->levels[st->level_count - 1].z;
	for (i = 1; i < sections; i++)
	{
		double z = z0 + (z1 - z0) * i / sections;
		qaws_clip_result* r = NULL;
		if (qaws_stack_section_2d(st, (qaws_scalar)z, &r) == QAWS_STATUS_OK)
			sk_draw_region(s, sv, r, z, scol, 1.0, 0.75);
		qaws_clip_result_destroy(r);
	}
	for (i = 0; i < st->level_count; i++)
	{
		qaws_clip_result* r = NULL;
		if (qaws_stack_section_2d(st, st->levels[i].z, &r) == QAWS_STATUS_OK)
			sk_draw_region(s, sv, r, st->levels[i].z, "#24292f", 1.6, 0.95);
		qaws_clip_result_destroy(r);
	}
}

typedef struct sk_store
{
	qaws_curve* c[128];
	qaws_curve const* v[128];
	qaws_path_2d p[128];
	qaws_stack_level l[128];
	unsigned int n, nl;
} sk_store;

static void sk_level1(sk_store* s, double z, qaws_curve* c)
{
	unsigned int i = s->n++;
	s->c[i] = c; s->v[i] = c;
	s->p[i].curves = &s->v[i]; s->p[i].curve_count = 1; s->p[i].closed = 1;
	s->l[s->nl].z = (qaws_scalar)z; s->l[s->nl].paths = &s->p[i]; s->l[s->nl].path_count = 1;
	s->nl++;
}

static qaws_curve* sk_circ(double x, double y, double r)
{
	qaws_curve* c = NULL;
	if (r < 1e-9)
	{
		qaws_scalar p[4] = { (qaws_scalar)x, (qaws_scalar)y, (qaws_scalar)x, (qaws_scalar)y };
		qaws_curve_create_polyline_2d(p, 2, 1, &c);
	}
	else
		qaws_curve_create_ellipse_2d((qaws_scalar)x, (qaws_scalar)y, (qaws_scalar)r, (qaws_scalar)r, 0, &c);
	return c;
}

static void sk_free(sk_store* s)
{
	unsigned int i;
	for (i = 0; i < s->n; i++) qaws_curve_destroy(s->c[i]);
}

static qaws_stack_2d sk_make(sk_store* s, qaws_stack_interp in)
{
	qaws_stack_2d st;
	memset(&st, 0, sizeof(st));
	st.levels = s->l; st.level_count = s->nl; st.fill_rule = QAWS_FILL_NON_ZERO; st.interp = in;
	return st;
}

static void demo_stack(void)
{
	int const cols = 3, rows = 2, pw = 400, ph = 330, gap = 16, top = 84;
	int const W = cols * pw + (cols + 1) * gap, H = top + rows * (ph + gap) + 8;
	svg s;
	sk_store cone, sph, two;
	qaws_stack_2d st_cone, st_sph_l, st_sph_c;
	unsigned int i, k;
	double t0;
	char label[200];
	if (!svg_open(&s, "showcase/b2d5_stack.svg", W, H, "2.5D: stacks of parallel contours",
		"qaws_stack: levels at uneven heights (dark), sections interpolated between them (light); matched contours blend exactly, distance fields handle splits; Booleans solved level by level"))
		return;
	memset(&cone, 0, sizeof(cone));
	memset(&sph, 0, sizeof(sph));
	memset(&two, 0, sizeof(two));
	sk_level1(&cone, 0, sk_circ(0, 0, 1.5));
	sk_level1(&cone, 2.2, sk_circ(0, 0, 0));
	st_cone = sk_make(&cone, QAWS_STACK_LINEAR);
	for (i = 0; i <= 8; i++)
	{
		double th = PI * (8 - i) / 8;
		sk_level1(&sph, 1.1 + 1.1 * cos(th), sk_circ(0.6, 0, 1.1 * sin(th)));
	}
	st_sph_l = sk_make(&sph, QAWS_STACK_LINEAR);
	st_sph_c = sk_make(&sph, QAWS_STACK_CUBIC);
	t0 = (double)clock() / CLOCKS_PER_SEC;
	for (k = 0; k < 6; k++)
	{
		sk_view sv;
		sv.v.x0 = gap + (k % cols) * (pw + gap); sv.v.y0 = top + (k / cols) * (ph + gap); sv.v.w = pw; sv.v.h = ph;
		sv.v.xmin = -2.6; sv.v.xmax = 3.4; sv.v.ymin = -2.2 + 0.4; sv.v.ymax = sv.v.ymin + (sv.v.xmax - sv.v.xmin) * ph / pw;
		sv.sx = 18; sv.sy = 62;
		if (k == 0)
		{
			svg_panel(&s, &sv.v, "cone: a circle and a point, 2 levels");
			sk_draw_stack(&s, &sv, &st_cone, 12, "#0969da");
		}
		else if (k == 1)
		{
			sprintf(label, "sphere from 9 circles, linear blend");
			svg_panel(&s, &sv.v, label);
			sk_draw_stack(&s, &sv, &st_sph_l, 30, "#1a7f37");
		}
		else if (k == 2)
		{
			sprintf(label, "the same 9 circles, cubic blend");
			svg_panel(&s, &sv.v, label);
			sk_draw_stack(&s, &sv, &st_sph_c, 30, "#1a7f37");
		}
		else if (k == 3)
		{
			/* a long ellipse at z 0 splitting into two circles at z 2 */
			qaws_curve* e = NULL;
			qaws_curve *c1 = sk_circ(-1.2, 0, 0.7), *c2 = sk_circ(1.6, 0, 0.9);
			qaws_curve const* v[3];
			qaws_path_2d p[3];
			qaws_stack_level lv[2];
			qaws_stack_2d st;
			qaws_curve_create_ellipse_2d((qaws_scalar)0.2, 0, (qaws_scalar)2.2, (qaws_scalar)1.0, 0, &e);
			v[0] = e; v[1] = c1; v[2] = c2;
			for (i = 0; i < 3; i++) { p[i].curves = &v[i]; p[i].curve_count = 1; p[i].closed = 1; }
			lv[0].z = 0; lv[0].paths = &p[0]; lv[0].path_count = 1;
			lv[1].z = 2; lv[1].paths = &p[1]; lv[1].path_count = 2;
			memset(&st, 0, sizeof(st));
			st.levels = lv; st.level_count = 2; st.fill_rule = QAWS_FILL_NON_ZERO; st.grid = 160;
			svg_panel(&s, &sv.v, "an ellipse splitting into two circles (distance blend)");
			sk_draw_stack(&s, &sv, &st, 12, "#8250df");
			qaws_curve_destroy(e); qaws_curve_destroy(c1); qaws_curve_destroy(c2);
		}
		else
		{
			/* cone x sphere and cone + sphere, level by level */
			qaws_stack_result* r = NULL;
			qaws_clip_type ct = k == 4 ? QAWS_CLIP_INTERSECTION : QAWS_CLIP_UNION;
			sprintf(label, "cone %s sphere, solved at every level of either", k == 4 ? "x" : "+");
			svg_panel(&s, &sv.v, label);
			if (qaws_stack_boolean_2d(ct, &st_cone, &st_sph_c, NULL, &r) == QAWS_STATUS_OK)
			{
				for (i = 0; i < qaws_stack_result_get_level_count(r); i++)
					sk_draw_region(&s, &sv, qaws_stack_result_get_level(r, i), qaws_stack_result_get_z(r, i),
						k == 4 ? "#cf222e" : "#bc4c00", 1.8, 0.95);
				sprintf(label, "%u levels", qaws_stack_result_get_level_count(r));
				svg_text(&s, sv.v.x0 + 10, sv.v.y0 + ph - 12, 12, "#57606a", "start", label);
			}
			qaws_stack_result_destroy(r);
		}
	}
	svg_close(&s);
	printf("    stack: 6 panels in %.3f s\n", (double)clock() / CLOCKS_PER_SEC - t0);
	sk_free(&cone);
	sk_free(&sph);
	sk_free(&two);
}
