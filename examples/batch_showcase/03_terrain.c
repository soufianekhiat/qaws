/* Figure 3: the heightfield of figure 1 as a B-spline terrain. Its contour
   lines as the intersection curves with a stack of horizontal planes, all in
   one call; and ballistic arcs dropped on it, every impact in one call. */

#define TR_GRID 24
#define TR_LEVELS 16
#define TR_ARCS 150

static qaws_surface* tr_terrain(void)
{
	static qaws_vec3 cp[TR_GRID * TR_GRID];
	qaws_surface_bspline_desc d;
	qaws_surface* s = NULL;
	unsigned int i, j;
	for (i = 0; i < TR_GRID; i++)
		for (j = 0; j < TR_GRID; j++)
		{
			double x = -1 + 2.0 * i / (TR_GRID - 1), y = -1 + 2.0 * j / (TR_GRID - 1);
			cp[i * TR_GRID + j].x = (qaws_scalar)x;
			cp[i * TR_GRID + j].y = (qaws_scalar)y;
			cp[i * TR_GRID + j].z = (qaws_scalar)(0.6 * hf_height(x, y, NULL, NULL));
		}
	memset(&d, 0, sizeof(d));
	d.u_degree = 3;
	d.v_degree = 3;
	d.control_points = cp;
	d.u_point_count = TR_GRID;
	d.v_point_count = TR_GRID;
	qaws_surface_create_bspline(&d, &s);
	return s;
}

static qaws_surface* tr_plane(double z)
{
	qaws_vec3 cp[4];
	qaws_surface_bezier_desc d;
	qaws_surface* s = NULL;
	unsigned int i, j;
	for (i = 0; i < 2; i++)
		for (j = 0; j < 2; j++)
		{
			cp[i * 2 + j].x = (qaws_scalar)(i ? 1.1 : -1.1);
			cp[i * 2 + j].y = (qaws_scalar)(j ? 1.1 : -1.1);
			cp[i * 2 + j].z = (qaws_scalar)z;
		}
	memset(&d, 0, sizeof(d));
	d.u_degree = 1;
	d.v_degree = 1;
	d.control_points = cp;
	d.u_point_count = 2;
	d.v_point_count = 2;
	qaws_surface_create_bezier(&d, &s);
	return s;
}

static double tr_value(void const* user, qaws_scalar u, qaws_scalar v)
{
	qaws_surface_eval_result r;
	qaws_surface_evaluate((qaws_surface const*)user, u, v, QAWS_SURFACE_EVAL_POSITION, &r);
	return r.position.z;
}

static void demo_terrain(void)
{
	qaws_surface* sf[1 + TR_LEVELS];
	unsigned int fam[1 + TR_LEVELS], i, j, nc = 0, np = 0, nh = 0, pair_curves = 0, pair_hits = 0;
	static qaws_surface_batch_curve cv[1024];
	static qaws_ssi_point pt[1 << 17];
	static qaws_curve_surface_batch_hit hits[4096];
	qaws_curve* arcs[TR_ARCS];
	qaws_surface_batch_desc sd;
	qaws_curve_surface_batch_desc cd;
	qaws_surface_batch_stats sst, cst;
	double t0, t_ssi, t_ssi_pair, t_cs, t_cs_pair, zlo = 1e30, zhi = -1e30;
	svg s;
	projection pr;
	char sub[400], buf[200];
	unsigned long long state = 0x2545F4914F6CDD1Dull;

	sf[0] = tr_terrain();
	fam[0] = 0;
	for (i = 0; i <= 40; i++)
		for (j = 0; j <= 40; j++)
		{
			double z = tr_value(sf[0], (qaws_scalar)(i / 40.0), (qaws_scalar)(j / 40.0));
			if (z < zlo) zlo = z;
			if (z > zhi) zhi = z;
		}
	for (i = 0; i < TR_LEVELS; i++)
	{
		sf[1 + i] = tr_plane(zlo + (zhi - zlo) * (i + 0.5) / TR_LEVELS);
		fam[1 + i] = 1;
	}

	/* contours: terrain x planes in one call */
	memset(&sd, 0, sizeof(sd));
	sd.surfaces = (qaws_surface const* const*)sf;
	sd.surface_count = 1 + TR_LEVELS;
	sd.families = fam;
	t0 = hf_now();
	qaws_surface_batch_find_intersections(&sd, cv, 1024, &nc, pt, 1 << 17, &np, &sst);
	t_ssi = hf_now() - t0;
	t0 = hf_now();
	for (i = 0; i < TR_LEVELS; i++)
	{
		static qaws_ssi_curve pc[256];
		static qaws_ssi_point pp[1 << 15];
		qaws_ssi_desc d;
		unsigned int n = 0;
		memset(&d, 0, sizeof(d));
		d.surface_a = sf[0];
		d.surface_b = sf[1 + i];
		qaws_surface_intersect(&d, pc, 256, &n, pp, 1 << 15);
		pair_curves += n;
	}
	t_ssi_pair = hf_now() - t0;

	/* ballistic arcs: quadratic Beziers falling across the terrain */
	for (i = 0; i < TR_ARCS; i++)
	{
		qaws_scalar cp[9];
		qaws_bezier_desc bd;
		double r[4];
		unsigned int k;
		for (k = 0; k < 4; k++)
		{
			state ^= state << 13; state ^= state >> 7; state ^= state << 17;
			r[k] = (double)(state >> 11) / 9007199254740992.0;
		}
		cp[0] = (qaws_scalar)(-0.95 + 1.9 * r[0]); cp[1] = (qaws_scalar)(-0.95 + 1.9 * r[1]); cp[2] = (qaws_scalar)(zhi + 0.2);
		cp[3] = (qaws_scalar)(cp[0] + 0.5 * (r[2] - 0.5)); cp[4] = (qaws_scalar)(cp[1] + 0.5 * (r[3] - 0.5)); cp[5] = (qaws_scalar)(zhi + 0.4);
		cp[6] = (qaws_scalar)(cp[0] + 1.0 * (r[2] - 0.5)); cp[7] = (qaws_scalar)(cp[1] + 1.0 * (r[3] - 0.5)); cp[8] = (qaws_scalar)(zlo - 0.3);
		memset(&bd, 0, sizeof(bd));
		bd.dimension = QAWS_DIMENSION_3D;
		bd.degree = 2;
		bd.control_points = cp;
		bd.control_point_count = 3;
		arcs[i] = NULL;
		qaws_curve_create_bezier(&bd, &arcs[i]);
	}
	memset(&cd, 0, sizeof(cd));
	cd.curves = (qaws_curve const* const*)arcs;
	cd.curve_count = TR_ARCS;
	cd.surfaces = (qaws_surface const* const*)sf;
	cd.surface_count = 1;
	t0 = hf_now();
	qaws_curve_surface_batch_find_intersections(&cd, hits, 4096, &nh, &cst);
	t_cs = hf_now() - t0;
	t0 = hf_now();
	for (i = 0; i < TR_ARCS; i++)
	{
		qaws_surface_curve_intersection b[16];
		unsigned int n = 0;
		qaws_surface_find_curve_intersections(sf[0], arcs[i], b, 16, &n);
		pair_hits += n;
	}
	t_cs_pair = hf_now() - t0;
	printf("    terrain: %u contour curves (%u points) from %u planes in %.3f s (%u patches, %u candidates); pairwise marching %u curves in %.3f s\n",
		nc, np, TR_LEVELS, t_ssi, sst.patch_count, sst.candidate_count, pair_curves, t_ssi_pair);
	printf("    terrain: %u impacts of %u arcs in %.4f s (%u segments, %u patches, %u candidates); pairwise %u in %.3f s\n",
		nh, TR_ARCS, t_cs, cst.segment_count, cst.patch_count, cst.candidate_count, pair_hits, t_cs_pair);

	sprintf(sub, "left: %u contour curves = terrain x %u planes, one call (%.0f ms; pair by pair %.0f ms, %u curves).   right: %u arcs, %u impacts in one call (%.1f ms; pair by pair %.0f ms)",
		nc, TR_LEVELS, t_ssi * 1000, t_ssi_pair * 1000, pair_curves, TR_ARCS, nh, t_cs * 1000, t_cs_pair * 1000);
	if (!svg_open(&s, "showcase/batch3_terrain.svg", 1400, 760, "Batched surface intersections on a B-spline terrain", sub))
		return;
	/* left: contours */
	pr.cx = 350; pr.cy = 360; pr.scale = 205; pr.zscale = 1.0;
	draw_quads(&s, &pr, sf[0], 44, 44, tr_value, sf[0], zlo, zhi);
	for (i = 0; i < nc && i < 1024; i++)
	{
		double* xy = (double*)malloc(2 * (cv[i].count + 1) * sizeof(double));
		unsigned int k;
		for (k = 0; k < cv[i].count; k++)
		{
			qaws_ssi_point const* p = &pt[cv[i].first + k];
			project(&pr, p->position.x, p->position.y, p->position.z + 0.004, &xy[2 * k], &xy[2 * k + 1]);
		}
		if (cv[i].closed)
		{
			xy[2 * k] = xy[0];
			xy[2 * k + 1] = xy[1];
			k++;
		}
		svg_polyline(&s, xy, (int)k, "#1b1f24", 1.3, 0.9, 0);
		free(xy);
	}
	svg_text(&s, 350, 720, 14, "#24292f", "middle", "contour lines: terrain x 16 horizontal planes (families: terrain / planes)");
	/* right: impacts */
	pr.cx = 1050; pr.cy = 360;
	draw_quads(&s, &pr, sf[0], 44, 44, tr_value, sf[0], zlo, zhi);
	for (i = 0; i < TR_ARCS; i++)
	{
		double xy[2 * 64];
		unsigned int k, hit_end = 64;
		qaws_range rg = qaws_curve_get_parameter_range(arcs[i]);
		double tend = rg.max_value;
		/* draw each arc up to its first impact */
		for (k = 0; k < nh && k < 4096; k++)
			if (hits[k].curve == i)
			{
				tend = hits[k].t;
				break;
			}
		for (k = 0; k < hit_end; k++)
		{
			qaws_eval_result_3d e;
			qaws_curve_evaluate_3d(arcs[i], (qaws_scalar)(rg.min_value + (tend - rg.min_value) * k / (hit_end - 1)), QAWS_EVAL_FLAG_POSITION, &e);
			project(&pr, e.position.x, e.position.y, e.position.z, &xy[2 * k], &xy[2 * k + 1]);
		}
		svg_polyline(&s, xy, (int)hit_end, "#0969da", 0.9, 0.55, 0);
	}
	for (i = 0; i < nh && i < 4096; i++)
	{
		double x, y;
		project(&pr, hits[i].position.x, hits[i].position.y, hits[i].position.z, &x, &y);
		svg_circle(&s, x, y, 3.0, "#cf222e", "#ffffff");
	}
	sprintf(buf, "%u ballistic arcs x terrain: %u impacts (red), arcs drawn to their first impact", TR_ARCS, nh);
	svg_text(&s, 1050, 720, 14, "#24292f", "middle", buf);
	svg_close(&s);
	printf("  -> showcase/batch3_terrain.svg\n");
	for (i = 0; i < TR_ARCS; i++)
		qaws_curve_destroy(arcs[i]);
	for (i = 0; i <= TR_LEVELS; i++)
		qaws_surface_destroy(sf[i]);
}
