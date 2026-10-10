/*
 * Rectangle clipping and Minkowski sums on the two Boolean engines.
 */

#include "qaws_rectclip.h"
#include "qaws_curve.h"
#include "internal/qaws_internal_clip.h"
#include <stdlib.h>
#include <string.h>

/* ======================================================================== */
/*  Rectangles                                                              */
/* ======================================================================== */

static qaws_status rc_clip(qaws_scalar const rect[4], qaws_path_2d const* paths, unsigned int n, int open,
	qaws_clip_result** out)
{
	qaws_scalar box[8];
	qaws_curve* r = NULL;
	qaws_curve const* rv;
	qaws_path_2d rp;
	qaws_clip_result* all = NULL;
	unsigned int i;
	qaws_status s;
	if (!rect || (!paths && n) || !out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	box[0] = rect[0]; box[1] = rect[1]; box[2] = rect[2]; box[3] = rect[1];
	box[4] = rect[2]; box[5] = rect[3]; box[6] = rect[0]; box[7] = rect[3];
	if (!(rect[2] > rect[0]) || !(rect[3] > rect[1]))
		return qaws_internal_clip_result_empty(out);
	s = qaws_curve_create_polyline_2d(box, 4, 1, &r);
	if (s != QAWS_STATUS_OK)
		return s;
	rv = r;
	rp.curves = &rv; rp.curve_count = 1; rp.closed = 1;
	s = qaws_internal_clip_result_empty(&all);
	/* every path on its own, as Clipper2's RectClip */
	for (i = 0; i < n && s == QAWS_STATUS_OK; i++)
	{
		qaws_clip_desc d;
		qaws_clip_result* one = NULL;
		memset(&d, 0, sizeof(d));
		if (open)
		{
			d.open_subjects = &paths[i];
			d.open_subject_count = 1;
		}
		else
		{
			d.subjects = &paths[i];
			d.subject_count = 1;
		}
		d.clips = &rp;
		d.clip_count = 1;
		d.clip_type = QAWS_CLIP_INTERSECTION;
		d.fill_rule = QAWS_FILL_NON_ZERO;
		s = qaws_clip_execute(&d, &one);
		if (s == QAWS_STATUS_OK)
			s = qaws_internal_clip_result_append(all, one);
	}
	qaws_curve_destroy(r);
	if (s != QAWS_STATUS_OK)
	{
		qaws_clip_result_destroy(all);
		return s;
	}
	*out = all;
	return QAWS_STATUS_OK;
}

qaws_status qaws_rect_clip_2d(qaws_scalar const rect[4], qaws_path_2d const* paths, unsigned int n, qaws_clip_result** out)
{
	return rc_clip(rect, paths, n, 0, out);
}

qaws_status qaws_rect_clip_lines_2d(qaws_scalar const rect[4], qaws_path_2d const* paths, unsigned int n, qaws_clip_result** out)
{
	return rc_clip(rect, paths, n, 1, out);
}

static qaws_status rc_clip64(int64_t const rect[4], qaws_path64 const* paths, unsigned int n, int open, qaws_clip64_result** out)
{
	int64_t box[8];
	qaws_path64 rp;
	qaws_clip64_result* all = NULL;
	unsigned int i;
	qaws_status s;
	if (!rect || (!paths && n) || !out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	if (!(rect[2] > rect[0]) || !(rect[3] > rect[1]))
		return qaws_internal_clip64_result_empty(out);
	box[0] = rect[0]; box[1] = rect[1]; box[2] = rect[2]; box[3] = rect[1];
	box[4] = rect[2]; box[5] = rect[3]; box[6] = rect[0]; box[7] = rect[3];
	rp.points = box;
	rp.point_count = 4;
	s = qaws_internal_clip64_result_empty(&all);
	for (i = 0; i < n && s == QAWS_STATUS_OK; i++)
	{
		qaws_clip64_desc d;
		qaws_clip64_result* one = NULL;
		memset(&d, 0, sizeof(d));
		if (open) { d.open_subjects = &paths[i]; d.open_subject_count = 1; }
		else { d.subjects = &paths[i]; d.subject_count = 1; }
		d.clips = &rp;
		d.clip_count = 1;
		d.clip_type = QAWS_CLIP_INTERSECTION;
		d.fill_rule = QAWS_FILL_NON_ZERO;
		s = qaws_clip64_execute(&d, &one);
		if (s == QAWS_STATUS_OK)
			s = qaws_internal_clip64_result_append(all, one);
	}
	if (s != QAWS_STATUS_OK)
	{
		qaws_clip64_result_destroy(all);
		return s;
	}
	*out = all;
	return QAWS_STATUS_OK;
}

qaws_status qaws_rect_clip64(int64_t const rect[4], qaws_path64 const* paths, unsigned int n, qaws_clip64_result** out)
{
	return rc_clip64(rect, paths, n, 0, out);
}

qaws_status qaws_rect_clip_lines64(int64_t const rect[4], qaws_path64 const* paths, unsigned int n, qaws_clip64_result** out)
{
	return rc_clip64(rect, paths, n, 1, out);
}

/* ======================================================================== */
/*  Minkowski                                                               */
/* ======================================================================== */

/* the quads of path edge i (q_i -> q_i+1) and pattern edge j (p_j -> p_j+1),
   each made counter-clockwise: 4 points per quad into out */
static unsigned int mk_quads(double const* pat, unsigned int np, double const* path, unsigned int nq, int closed, double sign,
	double* out)
{
	unsigned int i, j, m = 0, edges = closed ? nq : nq - 1;
	for (i = 0; i < edges; i++)
	{
		unsigned int i2 = (i + 1) % nq;
		for (j = 0; j < np; j++)
		{
			unsigned int j2 = (j + 1) % np, k;
			double q[8], area = 0.0;
			q[0] = path[2 * i] + sign * pat[2 * j];    q[1] = path[2 * i + 1] + sign * pat[2 * j + 1];
			q[2] = path[2 * i2] + sign * pat[2 * j];   q[3] = path[2 * i2 + 1] + sign * pat[2 * j + 1];
			q[4] = path[2 * i2] + sign * pat[2 * j2];  q[5] = path[2 * i2 + 1] + sign * pat[2 * j2 + 1];
			q[6] = path[2 * i] + sign * pat[2 * j2];   q[7] = path[2 * i + 1] + sign * pat[2 * j2 + 1];
			for (k = 0; k < 4; k++)
				area += q[2 * k] * q[2 * ((k + 1) % 4) + 1] - q[2 * ((k + 1) % 4)] * q[2 * k + 1];
			if (area < 0)
			{
				/* reverse to counter-clockwise */
				double t0 = q[2], t1 = q[3];
				q[2] = q[6]; q[3] = q[7];
				q[6] = t0; q[7] = t1;
			}
			memcpy(&out[8 * m], q, sizeof(q));
			m++;
		}
	}
	return m;
}

static qaws_status mk_float(qaws_scalar const* pattern, unsigned int np, qaws_scalar const* path, unsigned int nq, int closed,
	double sign, qaws_clip_result** out)
{
	double *pat, *pth, *quads;
	unsigned int nquad, i, k;
	qaws_curve** cs;
	qaws_curve const** views;
	qaws_path_2d* paths;
	qaws_status s = QAWS_STATUS_OK;
	if (!pattern || !path || !out || np < 1 || nq < 1)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	pat = (double*)malloc(sizeof(double) * 2 * np);
	pth = (double*)malloc(sizeof(double) * 2 * nq);
	quads = (double*)malloc(sizeof(double) * 8 * (nq * np + 1));
	if (!pat || !pth || !quads)
	{
		free(pat); free(pth); free(quads);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < 2 * np; i++) pat[i] = pattern[i];
	for (i = 0; i < 2 * nq; i++) pth[i] = path[i];
	nquad = nq >= 2 ? mk_quads(pat, np, pth, nq, closed, sign, quads) : 0;
	cs = (qaws_curve**)calloc(nquad + 1, sizeof(qaws_curve*));
	views = (qaws_curve const**)calloc(nquad + 1, sizeof(qaws_curve*));
	paths = (qaws_path_2d*)calloc(nquad + 1, sizeof(qaws_path_2d));
	if (!cs || !views || !paths)
		s = QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < nquad && s == QAWS_STATUS_OK; i++)
	{
		qaws_scalar q[8];
		for (k = 0; k < 8; k++) q[k] = (qaws_scalar)quads[8 * i + k];
		s = qaws_curve_create_polyline_2d(q, 4, 1, &cs[i]);
		views[i] = cs[i];
		paths[i].curves = &views[i];
		paths[i].curve_count = 1;
		paths[i].closed = 1;
	}
	if (s == QAWS_STATUS_OK)
		s = qaws_clip_boolean(QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO, paths, nquad, NULL, 0, out);
	for (i = 0; cs && i < nquad; i++) qaws_curve_destroy(cs[i]);
	free(cs); free((void*)views); free(paths);
	free(pat); free(pth); free(quads);
	return s;
}

static qaws_status mk_64(int64_t const* pattern, unsigned int np, int64_t const* path, unsigned int nq, int closed,
	double sign, qaws_clip64_result** out)
{
	unsigned int nquad, i;
	int64_t* q64;
	qaws_path64* paths;
	qaws_clip64_desc d;
	qaws_status s;
	if (!pattern || !path || !out || np < 1 || nq < 1)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	nquad = nq >= 2 ? (closed ? nq : nq - 1) * np : 0;
	q64 = (int64_t*)malloc(sizeof(int64_t) * 8 * (nquad + 1));
	paths = (qaws_path64*)malloc(sizeof(qaws_path64) * (nquad + 1));
	if (!q64 || !paths)
	{
		free(q64); free(paths);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	{
		unsigned int edges = closed ? nq : nq - 1, m = 0, j;
		int64_t sg = sign < 0 ? -1 : 1;
		for (i = 0; nq >= 2 && i < edges; i++)
		{
			unsigned int i2 = (i + 1) % nq;
			for (j = 0; j < np; j++)
			{
				unsigned int j2 = (j + 1) % np;
				int64_t* q = &q64[8 * m];
				double area = 0.0;
				q[0] = path[2 * i] + sg * pattern[2 * j];    q[1] = path[2 * i + 1] + sg * pattern[2 * j + 1];
				q[2] = path[2 * i2] + sg * pattern[2 * j];   q[3] = path[2 * i2 + 1] + sg * pattern[2 * j + 1];
				q[4] = path[2 * i2] + sg * pattern[2 * j2];  q[5] = path[2 * i2 + 1] + sg * pattern[2 * j2 + 1];
				q[6] = path[2 * i] + sg * pattern[2 * j2];   q[7] = path[2 * i + 1] + sg * pattern[2 * j2 + 1];
				area = qaws_path64_area2(q, 4);
				if (area < 0)
				{
					int64_t t0 = q[2], t1 = q[3];
					q[2] = q[6]; q[3] = q[7];
					q[6] = t0; q[7] = t1;
				}
				paths[m].points = q;
				paths[m].point_count = 4;
				m++;
			}
		}
		nquad = m;
	}

	memset(&d, 0, sizeof(d));
	d.subjects = paths;
	d.subject_count = nquad;
	d.clip_type = QAWS_CLIP_UNION;
	d.fill_rule = QAWS_FILL_NON_ZERO;
	s = qaws_clip64_execute(&d, out);
	free(q64); free(paths);
	return s;
}

qaws_status qaws_minkowski_sum_2d(qaws_scalar const* pattern, unsigned int np, qaws_scalar const* path, unsigned int nq, int closed,
	qaws_clip_result** out)
{
	return mk_float(pattern, np, path, nq, closed, 1.0, out);
}

qaws_status qaws_minkowski_diff_2d(qaws_scalar const* pattern, unsigned int np, qaws_scalar const* path, unsigned int nq, int closed,
	qaws_clip_result** out)
{
	return mk_float(pattern, np, path, nq, closed, -1.0, out);
}

qaws_status qaws_minkowski_sum64(int64_t const* pattern, unsigned int np, int64_t const* path, unsigned int nq, int closed,
	qaws_clip64_result** out)
{
	return mk_64(pattern, np, path, nq, closed, 1.0, out);
}

qaws_status qaws_minkowski_diff64(int64_t const* pattern, unsigned int np, int64_t const* path, unsigned int nq, int closed,
	qaws_clip64_result** out)
{
	return mk_64(pattern, np, path, nq, closed, -1.0, out);
}
