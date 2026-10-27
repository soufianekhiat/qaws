/*
 * Attributes a sweep generates (qaws_sweep_desc.generate): per result corner,
 * from the time its place was swept and the body's point that swept it.
 */

#include "../internal/qaws_internal_memory.h"
#include "../qaws_mesh_sweep.h"
#include "qaws_mesh_internal.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define SG_STEPS 2048u   /* the motion's turned angle and travelled distance tabulated over [0, 1] */

/* the turned angle and the travelled distance of a rigid motion at SG_STEPS + 1 times (trapezoids of the speeds) */
static qaws_status sg_table(qaws_motion const* m, double* ang, double* dist)
{
	unsigned int i;
	double wp = 0, vp = 0;
	ang[0] = dist[0] = 0;
	for (i = 0; i <= SG_STEPS; i++)
	{
		qaws_motion_state s;
		double w, v;
		qaws_status st = qaws_motion_evaluate(&s, m, (qaws_scalar)((double)i / SG_STEPS));
		if (st != QAWS_STATUS_OK)
			return st;
		w = sqrt((double)s.angular_velocity.x * s.angular_velocity.x + (double)s.angular_velocity.y * s.angular_velocity.y + (double)s.angular_velocity.z * s.angular_velocity.z);
		v = sqrt((double)s.velocity.x * s.velocity.x + (double)s.velocity.y * s.velocity.y + (double)s.velocity.z * s.velocity.z);
		if (i)
		{
			ang[i] = ang[i - 1] + 0.5 * (w + wp) / SG_STEPS;
			dist[i] = dist[i - 1] + 0.5 * (v + vp) / SG_STEPS;
		}
		wp = w;
		vp = v;
	}
	return QAWS_STATUS_OK;
}

static double sg_at(double const* tab, double t)
{
	double x = fmin(fmax(t, 0), 1) * SG_STEPS;
	unsigned int i = (unsigned int)x;
	if (i >= SG_STEPS)
		return tab[SG_STEPS];
	return tab[i] + (x - i) * (tab[i + 1] - tab[i]);
}

unsigned int swp_combine_of(qaws_sweep_desc const* desc, char const* name)
{
	unsigned int i;
	if (!desc || !name)
		return QAWS_SWEEP_COMBINE_FIRST;
	for (i = 0; i < desc->combine_count; i++)
		if (desc->combines[i].name && !strcmp(desc->combines[i].name, name))
			return (unsigned int)desc->combines[i].combine;
	return QAWS_SWEEP_COMBINE_FIRST;
}

int swp_passes_wanted(qaws_sweep_desc const* desc)
{
	unsigned int i;
	if (!desc)
		return 0;
	if (desc->generate & (QAWS_SWEEP_GENERATE_PASSES | QAWS_SWEEP_GENERATE_DWELL))
		return 1;
	for (i = 0; desc->combines && i < desc->combine_count; i++)
		if (desc->combines[i].combine != QAWS_SWEEP_COMBINE_FIRST)
			return 1;
	return 0;
}

void swp_combine(double* out, double const* vals, unsigned int np, unsigned int n, unsigned int rule, unsigned int combine)
{
	unsigned int p, k;
	if (!np)
	{
		for (k = 0; k < n; k++)
			out[k] = 0;
		return;
	}
	if (combine == QAWS_SWEEP_COMBINE_AVERAGE && np > 1)
	{
		unsigned int pos[256];
		double w[256];
		unsigned int m = np < 256 ? np : 256;
		for (p = 0; p < m; p++)
		{
			pos[p] = p;
			w[p] = 1.0 / m;
		}
		hm_blend(out, vals, n, pos, w, m, n, rule);
		return;
	}
	for (k = 0; k < n; k++)
	{
		double v = vals[k];
		if (combine == QAWS_SWEEP_COMBINE_LAST)
			v = vals[(size_t)n * (np - 1) + k];
		else if (combine == QAWS_SWEEP_COMBINE_MINIMUM || combine == QAWS_SWEEP_COMBINE_MAXIMUM || combine == QAWS_SWEEP_COMBINE_SUM)
			for (p = 1; p < np; p++)
			{
				double x = vals[(size_t)n * p + k];
				v = combine == QAWS_SWEEP_COMBINE_SUM ? v + x : combine == QAWS_SWEEP_COMBINE_MINIMUM ? fmin(v, x) : fmax(v, x);
			}
		out[k] = v;
	}
}

double* swp_angle_table(qaws_motion const* m)
{
	double* tab = (double*)qaws_internal_scratch_malloc(sizeof(double) * 2 * (SG_STEPS + 1));
	if (tab && sg_table(m, tab, tab + SG_STEPS + 1) != QAWS_STATUS_OK)
	{
		qaws_internal_free(tab);
		tab = NULL;
	}
	return tab;
}

double swp_table_total(double const* tab)
{
	return tab[SG_STEPS];
}

double swp_table_time(double const* ang, double angle)
{
	unsigned int lo = 0, hi = SG_STEPS;
	if (!(ang[SG_STEPS] > 0) || angle <= 0)
		return 0;
	if (angle >= ang[SG_STEPS])
		return 1;
	while (hi - lo > 1)
	{
		unsigned int mid = (lo + hi) / 2;
		if (ang[mid] <= angle)
			lo = mid;
		else
			hi = mid;
	}
	return (lo + (ang[hi] > ang[lo] ? (angle - ang[lo]) / (ang[hi] - ang[lo]) : 0)) / SG_STEPS;
}

qaws_attribute_track const* swp_track(qaws_sweep_desc const* desc, char const* name)
{
	unsigned int i;
	if (!desc || !name)
		return NULL;
	for (i = 0; i < desc->track_count; i++)
		if (desc->tracks[i].name && !strcmp(desc->tracks[i].name, name))
			return &desc->tracks[i];
	return NULL;
}

qaws_status swp_tracks_check(qaws_sweep_desc const* desc, qaws_mesh_3d const* mesh)
{
	unsigned int i, k, id;
	if (!desc || !desc->track_count)
		return QAWS_STATUS_OK;
	if (!desc->tracks)
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < desc->track_count; i++)
	{
		qaws_attribute_track const* tr = &desc->tracks[i];
		qaws_attribute_desc d;
		if (!tr->name || !tr->key_count || !tr->times || !tr->values || qaws_mesh_3d_find_attribute(&id, mesh, tr->name) != QAWS_STATUS_OK)
			return QAWS_STATUS_INVALID_ARGUMENT;
		qaws_mesh_3d_get_attribute_desc(&d, mesh, id);
		if (d.element == QAWS_MESH_EDGE)
			return QAWS_STATUS_INVALID_ARGUMENT;
		for (k = 0; k < tr->key_count; k++)
			if (!isfinite((double)tr->times[k]) || (k && !((double)tr->times[k] > (double)tr->times[k - 1])))
				return QAWS_STATUS_INVALID_ARGUMENT;
	}
	return QAWS_STATUS_OK;
}

void swp_track_value(qaws_attribute_track const* tr, unsigned int count, unsigned int n, unsigned int rule, unsigned int e, double t, double* out)
{
	unsigned int k = 0, lo = 0, hi = tr->key_count - 1, j;
	double a[16], b[16], s = 0;
	if (tr->key_count == 1 || t <= (double)tr->times[0])
		k = 0;
	else if (t >= (double)tr->times[hi])
		k = hi;
	else
	{
		while (hi - lo > 1)
		{
			unsigned int mid = (lo + hi) / 2;
			if ((double)tr->times[mid] <= t)
				lo = mid;
			else
				hi = mid;
		}
		k = lo;
		s = (t - (double)tr->times[lo]) / ((double)tr->times[hi] - (double)tr->times[lo]);
	}
	for (j = 0; j < n; j++)
	{
		a[j] = (double)tr->values[((size_t)k * count + e) * n + j];
		b[j] = s > 0 ? (double)tr->values[((size_t)(k + 1) * count + e) * n + j] : a[j];
	}
	if (s > 0)
		hm_lerp(out, a, b, s, n, rule);
	else
		memcpy(out, a, sizeof(double) * n);
}

double swp_time_at_angle(qaws_motion const* m, double angle)
{
	double* ang = (double*)qaws_internal_scratch_malloc(sizeof(double) * 2 * (SG_STEPS + 1));
	double t = 0;
	unsigned int lo = 0, hi = SG_STEPS;
	if (!ang || sg_table(m, ang, ang + SG_STEPS + 1) != QAWS_STATUS_OK || !(ang[SG_STEPS] > 0))
	{
		qaws_internal_free(ang);
		return 0;
	}
	if (angle >= ang[SG_STEPS])
		t = 1;
	else if (angle > 0)
	{
		while (hi - lo > 1)
		{
			unsigned int mid = (lo + hi) / 2;
			if (ang[mid] <= angle)
				lo = mid;
			else
				hi = mid;
		}
		t = (lo + (ang[hi] > ang[lo] ? (angle - ang[lo]) / (ang[hi] - ang[lo]) : 0)) / SG_STEPS;
	}
	qaws_internal_free(ang);
	return t;
}

/* a rotation's quaternion (w x y z, w not negative) */
static void sg_quat(double const* R, double* q)
{
	double tr = R[0] + R[4] + R[8], s;
	unsigned int k;
	if (tr > 0)
	{
		s = 2 * sqrt(tr + 1);
		q[0] = s / 4; q[1] = (R[7] - R[5]) / s; q[2] = (R[2] - R[6]) / s; q[3] = (R[3] - R[1]) / s;
	}
	else if (R[0] > R[4] && R[0] > R[8])
	{
		s = 2 * sqrt(1 + R[0] - R[4] - R[8]);
		q[0] = (R[7] - R[5]) / s; q[1] = s / 4; q[2] = (R[1] + R[3]) / s; q[3] = (R[2] + R[6]) / s;
	}
	else if (R[4] > R[8])
	{
		s = 2 * sqrt(1 + R[4] - R[0] - R[8]);
		q[0] = (R[2] - R[6]) / s; q[1] = (R[1] + R[3]) / s; q[2] = s / 4; q[3] = (R[5] + R[7]) / s;
	}
	else
	{
		s = 2 * sqrt(1 + R[8] - R[0] - R[4]);
		q[0] = (R[3] - R[1]) / s; q[1] = (R[2] + R[6]) / s; q[2] = (R[5] + R[7]) / s; q[3] = s / 4;
	}
	if (q[0] < 0)
		for (k = 0; k < 4; k++)
			q[k] = -q[k];
}

typedef struct sg_attr
{
	unsigned int flag;
	char const* name;
	unsigned int n;
	qaws_attribute_rule rule;
} sg_attr;

static sg_attr const sg_list[] = {
	{ QAWS_SWEEP_GENERATE_TIME, "time", 1, QAWS_ATTRIBUTE_LINEAR },
	{ QAWS_SWEEP_GENERATE_INDEX, "index", 1, QAWS_ATTRIBUTE_NEAREST },
	{ QAWS_SWEEP_GENERATE_ANGLE, "angle", 1, QAWS_ATTRIBUTE_LINEAR },
	{ QAWS_SWEEP_GENERATE_ORIENTATION, "orientation", 4, QAWS_ATTRIBUTE_QUATERNION },
	{ QAWS_SWEEP_GENERATE_POSITION, "pose_position", 3, QAWS_ATTRIBUTE_LINEAR },
	{ QAWS_SWEEP_GENERATE_VELOCITY, "velocity", 3, QAWS_ATTRIBUTE_VECTOR },
	{ QAWS_SWEEP_GENERATE_DISTANCE, "distance", 1, QAWS_ATTRIBUTE_LINEAR },
	{ QAWS_SWEEP_GENERATE_SOURCE, "source_face", 1, QAWS_ATTRIBUTE_NEAREST },
	{ QAWS_SWEEP_GENERATE_SOURCE, "source_uv", 2, QAWS_ATTRIBUTE_LINEAR }
};

static int sg_cmp_time(void const* a, void const* b, void const* ctx)
{
	double const* ct = (double const*)ctx;
	double x = ct[*(unsigned int const*)a], y = ct[*(unsigned int const*)b];
	return x < y ? -1 : x > y;
}

/* (qsort without a context: the times through a file-scope pointer would not be thread-safe, so an insertion-merge sort here) */
static void sg_sort(unsigned int* idx, unsigned int* tmp, unsigned int n, double const* ct)
{
	unsigned int w, i;
	for (i = 0; i < n; i++)
		idx[i] = i;
	for (w = 1; w < n; w *= 2)
	{
		for (i = 0; i < n; i += 2 * w)
		{
			unsigned int a = i, b = i + w < n ? i + w : n, ea = b, eb = i + 2 * w < n ? i + 2 * w : n, k = i;
			while (a < ea || b < eb)
				tmp[k++] = b >= eb || (a < ea && sg_cmp_time(&idx[a], &idx[b], ct) <= 0) ? idx[a++] : idx[b++];
		}
		memcpy(idx, tmp, sizeof(unsigned int) * n);
	}
}

#define SG_REC 18u   /* per pass: time, index, angle, orientation (4), pose position (3), velocity (3), distance, face, uv (2) */

qaws_status swp_generate(qaws_mesh_3d* out, qaws_mesh_3d const* mesh, qaws_motion const* motion, unsigned int flags, unsigned int nt,
	unsigned int const* cpo, unsigned int const* csrc, double const* cw, double const* ct, double const* cdw, qaws_sweep_desc const* desc)
{
	unsigned int nto = qaws_mesh_3d_get_triangle_count(out), nc = 3 * nto, np = cpo ? cpo[nc] : nc, nv = qaws_mesh_3d_get_vertex_count(mesh),
		nf = qaws_mesh_3d_get_face_count(mesh), a, i, k;
	int pv = qaws_motion_is_per_vertex(motion);
	double *tab = NULL, *val = NULL, *pvv = NULL;
	unsigned int *tri = NULL, *tf = NULL, *fs = NULL, *fv = NULL, *idx = NULL, *tmp = NULL, maxp = 1;
	qaws_scalar *x = NULL, *vel = NULL, *pos = NULL, *sv = NULL;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	double last = -1;
	if (!flags)
		return QAWS_STATUS_OK;
	for (i = 0; cpo && i < nc; i++)
		maxp = cpo[i + 1] - cpo[i] > maxp ? cpo[i + 1] - cpo[i] : maxp;
	tab = (double*)qaws_internal_scratch_malloc(sizeof(double) * 2 * (SG_STEPS + 1));
	val = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), SG_REC, (size_t)np + 1));
	pvv = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 4, (size_t)maxp + 1));
	tri = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)nt + 1));
	tf = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)nt + 1));
	fs = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)nf + 1));
	fv = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)qaws_mesh_3d_get_corner_count(mesh) + 1));
	idx = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 2, (size_t)np + 1));
	x = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)nv + 1));
	vel = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)nv + 1));
	pos = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)nv + 1));
	sv = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 4, (size_t)nc + 1));
	if (!tab || !val || !pvv || !tri || !tf || !fs || !fv || !idx || !x || !vel || !pos || !sv)
		goto done;
	tmp = idx + np;
	qaws_mesh_3d_get_triangles(tri, mesh);
	qaws_mesh_3d_get_positions(x, mesh);
	qaws_mesh_3d_get_faces(fv, fs, mesh);
	/* (each triangle's face: faces fanned from their first corner) */
	for (i = 0, k = 0; i < nf; i++)
	{
		unsigned int j;
		for (j = 0; j + 2 < fs[i] && k < nt; j++)
			tf[k++] = i;
	}
	if (!pv && (st = sg_table(motion, tab, tab + SG_STEPS + 1)) != QAWS_STATUS_OK)
		goto done;
	/* every pass's values, the passes in time order (a vertex motion's velocities evaluated once per time) */
	sg_sort(idx, tmp, np, ct);
	for (i = 0; i < np; i++)
	{
		unsigned int c = idx[i], f = csrc[c];
		double t = ct[c], p[3] = { 0, 0, 0 }, R[9], T[3] = { 0, 0, 0 }, v[3] = { 0, 0, 0 }, q[4] = { 1, 0, 0, 0 }, angle = 0, d = 0;
		double* o = val + SG_REC * (size_t)c;
		for (k = 0; k < 9; k++)
			R[k] = k % 4 == 0;
		for (k = 0; k < 3; k++)
		{
			unsigned int j;
			for (j = 0; j < 3; j++)
				p[j] += cw[3 * (size_t)c + k] * (double)x[3 * (size_t)tri[3 * (size_t)f + k] + j];
		}
		if (pv)
		{
			if (t != last && (st = qaws_motion_transform_points(pos, vel, motion, (qaws_scalar)t, NULL, nv)) != QAWS_STATUS_OK)
				goto done;
			last = t;
			for (k = 0; k < 3; k++)
			{
				unsigned int j;
				for (j = 0; j < 3; j++)
					v[j] += cw[3 * (size_t)c + k] * (double)vel[3 * (size_t)tri[3 * (size_t)f + k] + j];
			}
		}
		else
		{
			qaws_motion_state s;
			qaws_scalar pp[3] = { (qaws_scalar)p[0], (qaws_scalar)p[1], (qaws_scalar)p[2] }, wp[3], wv[3];
			if ((st = qaws_motion_evaluate(&s, motion, (qaws_scalar)t)) != QAWS_STATUS_OK ||
				(st = qaws_motion_transform_points(wp, wv, motion, (qaws_scalar)t, pp, 1)) != QAWS_STATUS_OK)
				goto done;
			for (k = 0; k < 9; k++)
				R[k] = (double)s.rotation[k];
			T[0] = s.position.x; T[1] = s.position.y; T[2] = s.position.z;
			for (k = 0; k < 3; k++)
				v[k] = (double)wv[k];
			sg_quat(R, q);
			angle = sg_at(tab, t);
			d = sg_at(tab + SG_STEPS + 1, t);
		}
		o[0] = t;
		o[1] = floor(angle / (2 * 3.14159265358979323846) + 1e-9);
		o[2] = angle;
		for (k = 0; k < 4; k++)
			o[3 + k] = q[k];
		for (k = 0; k < 3; k++)
		{
			o[7 + k] = T[k];
			o[10 + k] = v[k];
		}
		o[13] = d;
		o[14] = tf[f];
		o[15] = cw[3 * (size_t)c + 1];
		o[16] = cw[3 * (size_t)c + 2];
		o[17] = 0;
	}
	st = QAWS_STATUS_OK;
	for (a = 0; st == QAWS_STATUS_OK && a < sizeof sg_list / sizeof *sg_list + 2; a++)
	{
		static unsigned int const at[9] = { 0, 1, 2, 3, 7, 10, 13, 14, 15 };
		static char const* const extra[2] = { "passes", "dwell" };
		unsigned int flag = a < 9 ? sg_list[a].flag : a == 9 ? (unsigned int)QAWS_SWEEP_GENERATE_PASSES : (unsigned int)QAWS_SWEEP_GENERATE_DWELL;
		char const* name = a < 9 ? sg_list[a].name : extra[a - 9];
		unsigned int n = a < 9 ? sg_list[a].n : 1, rule = a < 9 ? (unsigned int)sg_list[a].rule : a == 9 ? (unsigned int)QAWS_ATTRIBUTE_NEAREST : (unsigned int)QAWS_ATTRIBUTE_LINEAR;
		unsigned int comb = swp_combine_of(desc, name);
		qaws_attribute_desc d;
		unsigned int id, j;
		if (!(flags & flag))
			continue;
		for (i = 0; i < nc; i++)
		{
			unsigned int p0 = cpo ? cpo[i] : i, p1 = cpo ? cpo[i + 1] : i + 1, pp;
			double y[4];
			if (a >= 9)
			{
				/* (the corner's passes: their count, their summed time within reach) */
				double s2 = 0;
				for (pp = p0; pp < p1 && a == 10; pp++)
					s2 += cdw ? cdw[pp] : 0;
				sv[i] = (qaws_scalar)(a == 9 ? (double)(p1 - p0) : s2);
				continue;
			}
			if (comb == QAWS_SWEEP_COMBINE_FIRST)
				p1 = p0 + 1;
			else if (comb == QAWS_SWEEP_COMBINE_LAST)
				p0 = p1 - 1;
			for (pp = p0; pp < p1; pp++)
				for (j = 0; j < n; j++)
					pvv[(size_t)n * (pp - p0) + j] = val[SG_REC * (size_t)pp + at[a] + j];
			swp_combine(y, pvv, p1 - p0, n, rule, comb);
			for (j = 0; j < n; j++)
				sv[(size_t)n * i + j] = (qaws_scalar)y[j];
		}
		memset(&d, 0, sizeof d);
		d.name = name;
		d.element = QAWS_MESH_CORNER;
		d.components = n;
		d.rule = (qaws_attribute_rule)rule;
		if (qaws_mesh_3d_find_attribute(&id, out, d.name) == QAWS_STATUS_OK)
			qaws_mesh_3d_remove_attribute(out, id);   /* (a generated one given back: the new values) */
		st = qaws_mesh_3d_add_attribute(&id, out, sv, &d);
	}
done:
	qaws_internal_free(tab); qaws_internal_free(val); qaws_internal_free(pvv); qaws_internal_free(tri); qaws_internal_free(tf); qaws_internal_free(fs);
	qaws_internal_free(fv); qaws_internal_free(idx); qaws_internal_free(x); qaws_internal_free(vel); qaws_internal_free(pos); qaws_internal_free(sv);
	return st;
}
