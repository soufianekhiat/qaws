#include "test_common.h"

static void test_svg_rational_bezier(void)
{
	printf("test_svg_rational_bezier\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/28_rational_bezier.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-0.5, (qaws_scalar)2.5, (qaws_scalar)2.5,
		(qaws_scalar)500, (qaws_scalar)500))
		return;

	/* Draw several rational Bezier quarter-circle arcs with different weights */
	{
		static const qaws_scalar weights[] = {
			(qaws_scalar)0.3, (qaws_scalar)0.70710678118, (qaws_scalar)1.0, (qaws_scalar)2.0
		};
		static const char* colors[] = { "#6272a4", "#e94560", "#50fa7b", "#f5a623" };
		unsigned int wi;

		for (wi = 0; wi < 4; ++wi)
		{
			qaws_scalar cp[] = { 1, 0,   1, 1,   0, 1 };
			qaws_scalar w[] = { 1, 0, 1 };
			qaws_rational_bezier_desc desc;
			qaws_curve *crv = NULL;

			w[1] = weights[wi];
			memset(&desc, 0, sizeof(desc));
			desc.dimension = QAWS_DIMENSION_2D;
			desc.degree = 2;
			desc.control_points = cp;
			desc.control_point_count = 3;
			desc.weights = w;
			desc.weight_count = 3;

			if (qaws_curve_create_rational_bezier(&desc, &crv) == QAWS_STATUS_OK)
			{
				qaws_vec2 buf[SVG_SAMPLES];
				unsigned int n = svg_sample_curve(crv, buf, SVG_SAMPLES);
				svg_polyline(&svg, buf, n, colors[wi], (qaws_scalar)2.0);
				qaws_curve_destroy(crv);
			}
		}
	}

	/* Control polygon */
	{
		qaws_vec2 cp_vis[] = { {1,0}, {1,1}, {0,1} };
		svg_polyline(&svg, cp_vis, 3, "#ffffff30", (qaws_scalar)1.0);
		svg_dots(&svg, cp_vis, 3, "#ffffff60", (qaws_scalar)3);
	}

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)2.2,
		"Rational Bezier: w=0.3 0.707 1.0 2.0", "#8888aa");
	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/28_rational_bezier.svg\n");
}

static void test_svg_composite(void)
{
	printf("test_svg_composite\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/28_composite.svg",
		(qaws_scalar)-1.0, (qaws_scalar)-1.5, (qaws_scalar)7.0, (qaws_scalar)5.0,
		(qaws_scalar)600, (qaws_scalar)430))
		return;

	/* Composite: Bezier S-curve + Hermite loop + linear segment */
	{
		qaws_vec2 cp1[] = { {0,0}, {0,2}, {2,2}, {2,0} };
		qaws_bezier_desc bd;
		qaws_curve *seg1 = NULL;
		memset(&bd, 0, sizeof(bd));
		bd.dimension = QAWS_DIMENSION_2D;
		bd.degree = 3;
		bd.control_points = cp1;
		bd.control_point_count = 4;
		qaws_curve_create_bezier(&bd, &seg1);

		qaws_vec2 hp[] = { {2,0}, {4,2} };
		qaws_vec2 ht[] = { {2,2}, {2,-2} };
		qaws_hermite_desc hd;
		qaws_curve *seg2 = NULL;
		memset(&hd, 0, sizeof(hd));
		hd.dimension = QAWS_DIMENSION_2D;
		hd.degree = 3;
		hd.points = hp;
		hd.derivatives = ht;
		hd.point_count = 2;
		hd.derivative_count = 2;
		qaws_curve_create_hermite(&hd, &seg2);

		qaws_vec2 cp3[] = { {4,2}, {5,0} };
		qaws_bezier_desc bd3;
		qaws_curve *seg3 = NULL;
		memset(&bd3, 0, sizeof(bd3));
		bd3.dimension = QAWS_DIMENSION_2D;
		bd3.degree = 1;
		bd3.control_points = cp3;
		bd3.control_point_count = 2;
		qaws_curve_create_bezier(&bd3, &seg3);

		if (seg1 && seg2 && seg3)
		{
			qaws_curve *segs[3];
			qaws_composite_desc cd;
			qaws_curve *comp = NULL;

			segs[0] = seg1; segs[1] = seg2; segs[2] = seg3;
			memset(&cd, 0, sizeof(cd));
			cd.dimension = QAWS_DIMENSION_2D;
			cd.segments = segs;
			cd.segment_count = 3;

			if (qaws_curve_create_composite(&cd, &comp) == QAWS_STATUS_OK)
			{
				qaws_vec2 buf[SVG_SAMPLES];
				unsigned int n = svg_sample_curve(comp, buf, SVG_SAMPLES);
				svg_polyline(&svg, buf, n, "#e94560", (qaws_scalar)2.5);

				/* Mark segment boundaries */
				{
					qaws_eval_result_2d r;
					qaws_curve_evaluate_2d(comp, 0, QAWS_EVAL_FLAG_POSITION, &r);
					svg_dot(&svg, r.position.x, r.position.y, "#50fa7b", (qaws_scalar)5);
					qaws_curve_evaluate_2d(comp, 1, QAWS_EVAL_FLAG_POSITION, &r);
					svg_dot(&svg, r.position.x, r.position.y, "#f5a623", (qaws_scalar)5);
					qaws_curve_evaluate_2d(comp, 2, QAWS_EVAL_FLAG_POSITION, &r);
					svg_dot(&svg, r.position.x, r.position.y, "#f5a623", (qaws_scalar)5);
					qaws_curve_evaluate_2d(comp, 3, QAWS_EVAL_FLAG_POSITION, &r);
					svg_dot(&svg, r.position.x, r.position.y, "#50fa7b", (qaws_scalar)5);
				}
				qaws_curve_destroy(comp);
			}
			else
			{
				qaws_curve_destroy(seg1);
				qaws_curve_destroy(seg2);
				qaws_curve_destroy(seg3);
			}
		}
	}

	svg_label(&svg, (qaws_scalar)-0.5, (qaws_scalar)3.0,
		"Composite: Bezier + Hermite + Linear", "#8888aa");
	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/28_composite.svg\n");
}

static void test_svg_arc(void)
{
	printf("test_svg_arc\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/28_arc.svg",
		(qaws_scalar)-4.0, (qaws_scalar)-4.0, (qaws_scalar)12.0, (qaws_scalar)8.0,
		(qaws_scalar)750, (qaws_scalar)500))
		return;

	{
		qaws_scalar pi = (qaws_scalar)3.14159265358979;

		/* 1. Concentric quarter arcs at different radii */
		{
			static const qaws_scalar radii[] = { 0.5, 1.0, 1.5, 2.0, 2.5, 3.0 };
			static const char* colors[] = {
				"#e94560", "#f5a623", "#50fa7b", "#78dce8", "#ab9df2", "#ff6188"
			};
			unsigned int ri;
			for (ri = 0; ri < 6; ++ri) {
				qaws_arc_segment seg;
				qaws_arc_desc desc;
				qaws_curve *crv = NULL;
				memset(&seg, 0, sizeof(seg));
				seg.center[0] = -1; seg.center[1] = 0;
				seg.radius = radii[ri];
				seg.angle_start = (qaws_scalar)0;
				seg.angle_end = pi * (qaws_scalar)0.5;
				memset(&desc, 0, sizeof(desc));
				desc.dimension = QAWS_DIMENSION_2D;
				desc.segments = &seg;
				desc.segment_count = 1;
				if (qaws_curve_create_arc(&desc, &crv) == QAWS_STATUS_OK) {
					qaws_vec2 buf[SVG_SAMPLES];
					unsigned int n = svg_sample_curve(crv, buf, SVG_SAMPLES);
					svg_polyline(&svg, buf, n, colors[ri], (qaws_scalar)2.0);
					qaws_curve_destroy(crv);
				}
			}
			svg_dot(&svg, (qaws_scalar)-1, (qaws_scalar)0, "#ffffff", (qaws_scalar)3);
			svg_label(&svg, (qaws_scalar)-3.5, (qaws_scalar)3.5,
				"Concentric quarter arcs", "#8888aa");
		}

		/* 2. S-shaped path from alternating arcs */
		{
			qaws_arc_segment segs[4];
			qaws_arc_desc desc;
			qaws_curve *crv = NULL;
			unsigned int i;

			/* Arc 1: right turn, center at (4, 0) */
			memset(&segs[0], 0, sizeof(segs[0]));
			segs[0].center[0] = 4; segs[0].center[1] = 0;
			segs[0].radius = (qaws_scalar)1.5;
			segs[0].angle_start = pi;
			segs[0].angle_end = pi * (qaws_scalar)0.5;

			/* Arc 2: left turn, center at (4, 3) */
			memset(&segs[1], 0, sizeof(segs[1]));
			segs[1].center[0] = 4; segs[1].center[1] = 3;
			segs[1].radius = (qaws_scalar)1.5;
			segs[1].angle_start = pi * (qaws_scalar)1.5;
			segs[1].angle_end = pi * (qaws_scalar)2.0;

			/* Arc 3: right turn, center at (7, 3) */
			memset(&segs[2], 0, sizeof(segs[2]));
			segs[2].center[0] = 7; segs[2].center[1] = 3;
			segs[2].radius = (qaws_scalar)1.5;
			segs[2].angle_start = pi;
			segs[2].angle_end = pi * (qaws_scalar)0.5;

			/* Arc 4: left turn */
			memset(&segs[3], 0, sizeof(segs[3]));
			segs[3].center[0] = 7; segs[3].center[1] = 0;
			segs[3].radius = (qaws_scalar)1.5;
			segs[3].angle_start = pi * (qaws_scalar)0.5;
			segs[3].angle_end = 0;

			memset(&desc, 0, sizeof(desc));
			desc.dimension = QAWS_DIMENSION_2D;
			desc.segments = segs;
			desc.segment_count = 4;

			if (qaws_curve_create_arc(&desc, &crv) == QAWS_STATUS_OK) {
				qaws_vec2 buf[SVG_SAMPLES];
				unsigned int n = svg_sample_curve(crv, buf, SVG_SAMPLES);
				svg_polyline(&svg, buf, n, "#78dce8", (qaws_scalar)2.5);
				qaws_curve_destroy(crv);
			}

			/* Mark centers */
			for (i = 0; i < 4; ++i)
				svg_ring(&svg, segs[i].center[0], segs[i].center[1],
					"#78dce8", (qaws_scalar)3, (qaws_scalar)1);

			svg_label(&svg, (qaws_scalar)3.0, (qaws_scalar)3.5,
				"S-path (4 arcs)", "#8888aa");
		}

		/* 3. Full circle from varying-size arcs */
		{
			qaws_arc_segment segs[6];
			static const char* colors[] = {
				"#e94560", "#f5a623", "#50fa7b", "#78dce8", "#ab9df2", "#ff6188"
			};
			qaws_scalar a = 0;
			unsigned int i;

			for (i = 0; i < 6; ++i) {
				qaws_scalar span = pi * (qaws_scalar)2 / (qaws_scalar)6;
				memset(&segs[i], 0, sizeof(segs[i]));
				segs[i].center[0] = 4; segs[i].center[1] = -2;
				segs[i].radius = (qaws_scalar)1.5;
				segs[i].angle_start = a;
				segs[i].angle_end = a + span;
				a += span;
			}

			/* Draw each segment individually for different colors */
			for (i = 0; i < 6; ++i) {
				qaws_arc_desc d1;
				qaws_curve *seg_crv = NULL;
				memset(&d1, 0, sizeof(d1));
				d1.dimension = QAWS_DIMENSION_2D;
				d1.segments = &segs[i];
				d1.segment_count = 1;
				if (qaws_curve_create_arc(&d1, &seg_crv) == QAWS_STATUS_OK) {
					qaws_vec2 buf[SVG_SAMPLES];
					unsigned int n = svg_sample_curve(seg_crv, buf, SVG_SAMPLES);
					svg_polyline(&svg, buf, n, colors[i], (qaws_scalar)3.0);
					qaws_curve_destroy(seg_crv);
				}
			}
			svg_label(&svg, (qaws_scalar)2.5, (qaws_scalar)-3.5,
				"6-segment circle", "#8888aa");
		}
	}

	svg_label(&svg, (qaws_scalar)-3.5, (qaws_scalar)4.2,
		"Piecewise Circular Arcs", "#8888aa");
	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/28_arc.svg\n");
}

static void test_svg_polynomial(void)
{
	printf("test_svg_polynomial\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/28_polynomial.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-1.5, (qaws_scalar)7.0, (qaws_scalar)4.0,
		(qaws_scalar)600, (qaws_scalar)350))
		return;

	/* Lissajous-like: x(t) = 3t, y(t) = sin-like cubic approx */
	{
		/* P(t) = (3t, 4t(1-t)(2t-1)) for t in [0,2] */
		/* Coefficients: c0 + c1*t + c2*t^2 + c3*t^3 */
		/* x: 0 + 3t => c0x=0, c1x=3, c2x=0, c3x=0 */
		/* y: 4t(1-t)(2t-1) = 4t(2t-1-2t^2+t) = 4t(3t-1-2t^2) = 12t^2-4t-8t^3 */
		/* => c0y=0, c1y=-4, c2y=12, c3y=-8 */
		qaws_scalar coeffs[] = {
			0, 0,     /* c0: (0, 0) */
			3, -4,    /* c1: (3, -4) */
			0, 12,    /* c2: (0, 12) */
			0, -8     /* c3: (0, -8) */
		};
		qaws_polynomial_desc desc;
		qaws_curve *crv = NULL;

		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.coefficients = coeffs;
		desc.coefficient_count = 4;
		desc.t_min = 0;
		desc.t_max = 2;

		if (qaws_curve_create_polynomial(&desc, &crv) == QAWS_STATUS_OK)
		{
			qaws_vec2 buf[SVG_SAMPLES];
			unsigned int n = svg_sample_curve(crv, buf, SVG_SAMPLES);
			svg_polyline(&svg, buf, n, "#bd93f9", (qaws_scalar)2.5);
			qaws_curve_destroy(crv);
		}
	}

	/* Parabola: (t, t^2) */
	{
		qaws_scalar coeffs[] = { 1, 0,   2, 0,   0, 1 };
		qaws_polynomial_desc desc;
		qaws_curve *crv = NULL;

		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 2;
		desc.coefficients = coeffs;
		desc.coefficient_count = 3;
		desc.t_min = -1;
		desc.t_max = 1;

		if (qaws_curve_create_polynomial(&desc, &crv) == QAWS_STATUS_OK)
		{
			qaws_vec2 buf[SVG_SAMPLES];
			unsigned int n = svg_sample_curve(crv, buf, SVG_SAMPLES);
			svg_polyline(&svg, buf, n, "#f5a623", (qaws_scalar)2.0);
			qaws_curve_destroy(crv);
		}
	}

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)2.2,
		"Polynomial: cubic wave + parabola", "#8888aa");
	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/28_polynomial.svg\n");
}

static void test_svg_clothoid(void)
{
	printf("test_svg_clothoid\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/28_clothoid.svg",
		(qaws_scalar)-1.5, (qaws_scalar)-1.0, (qaws_scalar)4.0, (qaws_scalar)3.5,
		(qaws_scalar)600, (qaws_scalar)525))
		return;

	/* Classic Euler spiral: kappa goes from 0 to a positive value */
	{
		qaws_clothoid_desc desc;
		qaws_curve *crv = NULL;

		memset(&desc, 0, sizeof(desc));
		desc.origin_x = 0;
		desc.origin_y = 0;
		desc.start_angle = 0;
		desc.start_curvature = 0;
		desc.end_curvature = 3;
		desc.length = 4;

		if (qaws_curve_create_clothoid(&desc, &crv) == QAWS_STATUS_OK)
		{
			qaws_vec2 buf[SVG_SAMPLES];
			unsigned int n = svg_sample_curve(crv, buf, SVG_SAMPLES);
			svg_polyline(&svg, buf, n, "#e94560", (qaws_scalar)2.5);
			svg_dot(&svg, buf[0].x, buf[0].y, "#50fa7b", (qaws_scalar)4);
			svg_dot(&svg, buf[n-1].x, buf[n-1].y, "#f5a623", (qaws_scalar)4);
			qaws_curve_destroy(crv);
		}
	}

	/* Reverse spiral: kappa from positive to 0 */
	{
		qaws_clothoid_desc desc;
		qaws_curve *crv = NULL;

		memset(&desc, 0, sizeof(desc));
		desc.origin_x = 0;
		desc.origin_y = 0;
		desc.start_angle = 0;
		desc.start_curvature = 3;
		desc.end_curvature = 0;
		desc.length = 4;

		if (qaws_curve_create_clothoid(&desc, &crv) == QAWS_STATUS_OK)
		{
			qaws_vec2 buf[SVG_SAMPLES];
			unsigned int n = svg_sample_curve(crv, buf, SVG_SAMPLES);
			svg_polyline(&svg, buf, n, "#6272a4", (qaws_scalar)2.0);
			qaws_curve_destroy(crv);
		}
	}

	svg_label(&svg, (qaws_scalar)-1.2, (qaws_scalar)2.2,
		"Clothoid / Euler spiral", "#8888aa");
	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/28_clothoid.svg\n");
}

static void test_svg_subdivision(void)
{
	printf("test_svg_subdivision\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/28_subdivision.svg",
		(qaws_scalar)-1.0, (qaws_scalar)-1.0, (qaws_scalar)6.0, (qaws_scalar)5.0,
		(qaws_scalar)600, (qaws_scalar)500))
		return;

	/* Control polygon */
	qaws_scalar cp[] = { 0,0,  1,3,  2,1,  3.5,3.5,  5,0 };
	qaws_vec2 cp_vis[] = { {0,0}, {1,3}, {2,1}, {(qaws_scalar)3.5,(qaws_scalar)3.5}, {5,0} };
	svg_polyline(&svg, cp_vis, 5, "#ffffff20", (qaws_scalar)1.0);
	svg_dots(&svg, cp_vis, 5, "#ffffff50", (qaws_scalar)3);

	/* Chaikin */
	{
		qaws_subdivision_desc desc;
		qaws_curve *crv = NULL;
		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.scheme = QAWS_SUBDIVISION_CHAIKIN;
		desc.control_points = cp;
		desc.control_point_count = 5;
		desc.closed = 0;
		desc.refinement_levels = 6;

		if (qaws_curve_create_subdivision(&desc, &crv) == QAWS_STATUS_OK)
		{
			qaws_vec2 buf[SVG_SAMPLES];
			unsigned int n = svg_sample_curve(crv, buf, SVG_SAMPLES);
			svg_polyline(&svg, buf, n, "#e94560", (qaws_scalar)2.5);
			qaws_curve_destroy(crv);
		}
	}

	/* Lane-Riesenfeld order 3 */
	{
		qaws_subdivision_desc desc;
		qaws_curve *crv = NULL;
		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.scheme = QAWS_SUBDIVISION_LANE_RIESENFELD_3;
		desc.control_points = cp;
		desc.control_point_count = 5;
		desc.closed = 0;
		desc.refinement_levels = 6;

		if (qaws_curve_create_subdivision(&desc, &crv) == QAWS_STATUS_OK)
		{
			qaws_vec2 buf[SVG_SAMPLES];
			unsigned int n = svg_sample_curve(crv, buf, SVG_SAMPLES);
			svg_polyline(&svg, buf, n, "#50fa7b", (qaws_scalar)2.0);
			qaws_curve_destroy(crv);
		}
	}

	/* Lane-Riesenfeld order 4 */
	{
		qaws_subdivision_desc desc;
		qaws_curve *crv = NULL;
		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.scheme = QAWS_SUBDIVISION_LANE_RIESENFELD_4;
		desc.control_points = cp;
		desc.control_point_count = 5;
		desc.closed = 0;
		desc.refinement_levels = 6;

		if (qaws_curve_create_subdivision(&desc, &crv) == QAWS_STATUS_OK)
		{
			qaws_vec2 buf[SVG_SAMPLES];
			unsigned int n = svg_sample_curve(crv, buf, SVG_SAMPLES);
			svg_polyline(&svg, buf, n, "#6272a4", (qaws_scalar)2.0);
			qaws_curve_destroy(crv);
		}
	}

	svg_label(&svg, (qaws_scalar)-0.5, (qaws_scalar)4.0,
		"Subdivision: Chaikin(red) LR3(green) LR4(blue)", "#8888aa");
	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/28_subdivision.svg\n");
}

static void test_svg_all_new_families(void)
{
	printf("test_svg_all_new_families\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/28_all_new_families.svg",
		(qaws_scalar)-2.5, (qaws_scalar)-2.0, (qaws_scalar)10.0, (qaws_scalar)6.0,
		(qaws_scalar)800, (qaws_scalar)480))
		return;

	/* 1. Rational Bezier quarter circle */
	{
		qaws_scalar cp[] = { -2, -1,   -2, 1,   0, 1 };
		qaws_scalar w[] = { 1, (qaws_scalar)0.70710678118, 1 };
		qaws_rational_bezier_desc desc;
		qaws_curve *crv = NULL;
		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 2;
		desc.control_points = cp;
		desc.control_point_count = 3;
		desc.weights = w;
		desc.weight_count = 3;
		if (qaws_curve_create_rational_bezier(&desc, &crv) == QAWS_STATUS_OK)
		{
			qaws_vec2 buf[SVG_SAMPLES];
			unsigned int n = svg_sample_curve(crv, buf, SVG_SAMPLES);
			svg_polyline(&svg, buf, n, "#e94560", (qaws_scalar)2.0);
			qaws_curve_destroy(crv);
		}
	}
	svg_label(&svg, (qaws_scalar)-2.3, (qaws_scalar)1.5, "RatBez", "#e94560");

	/* 2. Circular arc (half circle) */
	{
		qaws_scalar pi = (qaws_scalar)3.14159265358979;
		qaws_arc_segment seg;
		qaws_arc_desc desc;
		qaws_curve *crv = NULL;
		memset(&seg, 0, sizeof(seg));
		seg.center[0] = 1; seg.center[1] = 0;
		seg.radius = (qaws_scalar)1.5;
		seg.angle_start = (qaws_scalar)0.0;
		seg.angle_end = pi;
		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.segments = &seg;
		desc.segment_count = 1;
		if (qaws_curve_create_arc(&desc, &crv) == QAWS_STATUS_OK)
		{
			qaws_vec2 buf[SVG_SAMPLES];
			unsigned int n = svg_sample_curve(crv, buf, SVG_SAMPLES);
			svg_polyline(&svg, buf, n, "#50fa7b", (qaws_scalar)2.0);
			qaws_curve_destroy(crv);
		}
	}
	svg_label(&svg, (qaws_scalar)0.3, (qaws_scalar)2.0, "Arc", "#50fa7b");

	/* 3. Clothoid */
	{
		qaws_clothoid_desc desc;
		qaws_curve *crv = NULL;
		memset(&desc, 0, sizeof(desc));
		desc.origin_x = 3;
		desc.origin_y = -1;
		desc.start_angle = (qaws_scalar)0.3;
		desc.start_curvature = 0;
		desc.end_curvature = 2;
		desc.length = 3;
		if (qaws_curve_create_clothoid(&desc, &crv) == QAWS_STATUS_OK)
		{
			qaws_vec2 buf[SVG_SAMPLES];
			unsigned int n = svg_sample_curve(crv, buf, SVG_SAMPLES);
			svg_polyline(&svg, buf, n, "#f5a623", (qaws_scalar)2.0);
			qaws_curve_destroy(crv);
		}
	}
	svg_label(&svg, (qaws_scalar)3.0, (qaws_scalar)2.5, "Clothoid", "#f5a623");

	/* 4. Polynomial parabola */
	{
		qaws_scalar coeffs[] = { 5, -1,   1, 0,   0, 1 };
		qaws_polynomial_desc desc;
		qaws_curve *crv = NULL;
		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 2;
		desc.coefficients = coeffs;
		desc.coefficient_count = 3;
		desc.t_min = (qaws_scalar)-1.5;
		desc.t_max = (qaws_scalar)1.5;
		if (qaws_curve_create_polynomial(&desc, &crv) == QAWS_STATUS_OK)
		{
			qaws_vec2 buf[SVG_SAMPLES];
			unsigned int n = svg_sample_curve(crv, buf, SVG_SAMPLES);
			svg_polyline(&svg, buf, n, "#bd93f9", (qaws_scalar)2.0);
			qaws_curve_destroy(crv);
		}
	}
	svg_label(&svg, (qaws_scalar)5.2, (qaws_scalar)1.5, "Poly", "#bd93f9");

	/* 5. Subdivision Chaikin */
	{
		qaws_scalar cp[] = { 5.5,3,  6,1,  6.5,3,  7,0.5 };
		qaws_subdivision_desc desc;
		qaws_curve *crv = NULL;
		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.scheme = QAWS_SUBDIVISION_CHAIKIN;
		desc.control_points = cp;
		desc.control_point_count = 4;
		desc.closed = 0;
		desc.refinement_levels = 6;
		if (qaws_curve_create_subdivision(&desc, &crv) == QAWS_STATUS_OK)
		{
			qaws_vec2 buf[SVG_SAMPLES];
			unsigned int n = svg_sample_curve(crv, buf, SVG_SAMPLES);
			svg_polyline(&svg, buf, n, "#ff79c6", (qaws_scalar)2.0);
			qaws_curve_destroy(crv);
		}
	}
	svg_label(&svg, (qaws_scalar)5.5, (qaws_scalar)3.5, "Subdiv", "#ff79c6");

	svg_label(&svg, (qaws_scalar)-2.0, (qaws_scalar)-1.5,
		"All new curve families", "#8888aa");
	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/28_all_new_families.svg\n");
}

int test_28_svg_families_main(void)
{
	svg_ensure_output_dir();
	test_svg_rational_bezier();
	test_svg_composite();
	test_svg_arc();
	test_svg_polynomial();
	test_svg_clothoid();
	test_svg_subdivision();
	test_svg_all_new_families();
	return 0;
}
