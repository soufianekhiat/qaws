/*
 * Figure 58: a torus's periods quantized together (QAWS_QUADS_GLOBAL): tori
 * of radius 1 round the hole and 0.3, 0.5 and 0.7 round the tube, all
 * quads. The flow rounds each edge's steps alone and its grid could close
 * round the hole only after a shift along the tube (a spiral); the two
 * periods are now chosen together on the map's energy. One straight edge
 * loop round the hole (red) and one round the tube (blue) drawn on each.
 * Printed: the quads, the edge loops and those that spiral.
 */

/* the straight edge loop from corner s (through valence 4 vertices) into path (its vertices, n), its turns round the hole
   (*it) and the tube (*ip); seen marks its corners both ways; 0 when it does not close */
static int ms_tp_loop(qaws_mesh_3d const* q, unsigned int const* fv, unsigned int const* c0, unsigned int const* pos, unsigned int nc, unsigned int s,
	unsigned char* seen, unsigned int* path, unsigned int* n, int* it, int* ip)
{
	unsigned int e = s, i;
	double wt = 0, wp = 0;
	*n = 0;
	do
	{
		unsigned int k = pos[e] & 255, m = pos[e] >> 8, a = fv[e], b = fv[c0[e] + (k + 1) % m], c = fv[c0[e] + (k + 2) % m], g, d;
		qaws_vec3 A = qaws_mesh_3d_get_vertex(q, a), B = qaws_mesh_3d_get_vertex(q, b);
		double dt = atan2(B.y, B.x) - atan2(A.y, A.x), dp = atan2(B.z, sqrt(B.x * B.x + B.y * B.y) - 1) - atan2(A.z, sqrt(A.x * A.x + A.y * A.y) - 1);
		wt += dt - 2 * MS_PI * floor(dt / (2 * MS_PI) + 0.5);
		wp += dp - 2 * MS_PI * floor(dp / (2 * MS_PI) + 0.5);
		seen[e] = 1;
		path[(*n)++] = a;
		for (g = 0; g < nc && !(fv[g] == c && fv[c0[g] + ((pos[g] & 255) + 1) % (pos[g] >> 8)] == b); g++)
			;
		for (i = 0; i < nc; i++)
			if (fv[i] == b && fv[c0[i] + ((pos[i] & 255) + 1) % (pos[i] >> 8)] == a)
				seen[i] = 1;
		if (g == nc)
			return 0;
		d = fv[c0[g] + ((pos[g] & 255) + 2) % (pos[g] >> 8)];
		for (i = 0; i < nc && !(fv[i] == b && fv[c0[i] + ((pos[i] & 255) + 1) % (pos[i] >> 8)] == d); i++)
			;
		e = i;
	} while (e < nc && e != s && *n < nc);
	*it = (int)floor(wt / (2 * MS_PI) + 0.5);
	*ip = (int)floor(wp / (2 * MS_PI) + 0.5);
	return e == s;
}

static void demo_torus_periods(void)
{
	static double const tube[3] = { 0.3, 0.5, 0.7 }, len[3] = { 0.17, 0.3, 0.22 };
	static double const tint[3] = { 178, 192, 212 }, wire[3] = { 50, 55, 70 }, red[3] = { 210, 40, 40 }, blue[3] = { 40, 80, 210 };
	er_canvas* c = er_create(1800, 640, 3, 252, 252, 250);
	unsigned int k;
	if (!c)
		return;
	for (k = 0; k < 3; k++)
	{
		qaws_mesh_3d* m = ms_torus(96, 48, 1, tube[k]);
		qaws_quadrangulate_desc d;
		qaws_mesh_report rep;
		qaws_mesh_3d* q = NULL;
		char title[100], sub[200];
		memset(&d, 0, sizeof d);
		d.method = QAWS_QUADS_GLOBAL;
		d.target_length = (qaws_scalar)len[k];
		d.pure = 1;
		if (m && qaws_mesh_quadrangulate(&q, &rep, m, &d) == QAWS_STATUS_OK)
		{
			unsigned int nf = qaws_mesh_3d_get_face_count(q), nc = qaws_mesh_3d_get_corner_count(q), f, i, j, s, loops = 0, spirals = 0, n;
			unsigned int *fv = (unsigned int*)malloc(sizeof(unsigned int) * (nc + 1)), *fs = (unsigned int*)malloc(sizeof(unsigned int) * (nf + 1));
			unsigned int *c0 = (unsigned int*)malloc(sizeof(unsigned int) * (nc + 1)), *pos = (unsigned int*)malloc(sizeof(unsigned int) * (nc + 1));
			unsigned int* path = (unsigned int*)malloc(sizeof(unsigned int) * (nc + 1));
			unsigned char* seen = (unsigned char*)calloc(nc + 1, 1);
			int drawn[2] = { 0, 0 }, it, ip;
			double* x = ms_positions(q);
			qaws_mesh_3d_get_faces(fv, fs, q);
			for (f = 0, i = 0; f < nf; f++)
				for (j = 0; j < fs[f]; j++, i++)
				{
					c0[i] = i - j;
					pos[i] = j | (fs[f] << 8);
				}
			sprintf(title, "radii 1 and %.1f", tube[k]);
			ms_panel(c, 600 * (int)k, 0, 600, 640, q, 0.3, 0.9, title, NULL);
			ms_draw(c, q, tint, NULL, 0, wire, 0.4);
			for (s = 0; s < nc; s++)
			{
				int closed;
				if (seen[s])
					continue;
				closed = ms_tp_loop(q, fv, c0, pos, nc, s, seen, path, &n, &it, &ip);
				loops++;
				if (!closed || !((abs(it) == 1 && ip == 0) || (it == 0 && abs(ip) == 1)))
				{
					spirals++;
					continue;
				}
				if (!drawn[it == 0])
				{
					drawn[it == 0] = 1;
					for (i = 0; i < n; i++)
						er_line(c, x + 3 * (size_t)path[i], x + 3 * (size_t)path[(i + 1) % n], it == 0 ? blue : red, 2.5);
				}
			}
			sprintf(sub, "%u quads, %u edge loops, %u spiral", qaws_mesh_3d_get_face_count(q), loops, spirals);
			er_text(c, 600 * (int)k + 14, 46, sub, 2, 95, 95, 95);
			printf("  %s: %s\n", title, sub);
			free(fv); free(fs); free(c0); free(pos); free(path); free(seen); free(x);
		}
		qaws_mesh_3d_destroy(q);
		qaws_mesh_3d_destroy(m);
	}
	er_save(c, "showcase/mesh58_torus_periods.png");
	printf("mesh58_torus_periods.png\n");
	er_destroy(c);
}
