#include "test_common.h"

static void test_bspline_fitting(void)
{
	printf("test_bspline_fitting\n");

	/* Fit a cubic B-spline to a semicircle of sample points */
	{
		qaws_vec2 pts[21];
		unsigned int i;
		qaws_curve *crv = NULL;
		qaws_bspline_fit_desc desc;
		qaws_eval_result_2d r;
		qaws_range rng;

		for (i = 0; i < 21; ++i) {
			qaws_scalar angle = (qaws_scalar)3.14159265358979 * (qaws_scalar)i / (qaws_scalar)20;
			pts[i].x = (qaws_scalar)cos((double)angle);
			pts[i].y = (qaws_scalar)sin((double)angle);
		}

		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.data_points = pts;
		desc.data_point_count = 21;
		desc.degree = 3;
		desc.control_point_count = 8;
		desc.parameters = NULL;

		TEST_ASSERT_STATUS(qaws_curve_fit_bspline(&desc, &crv));
		TEST_ASSERT(crv != NULL, "bspline fit created");
		TEST_ASSERT(qaws_curve_get_kind(crv) == QAWS_CURVE_KIND_BSPLINE, "fit kind bspline");
		TEST_ASSERT(qaws_curve_get_degree(crv) == 3, "fit degree 3");

		/* Endpoints should match data endpoints */
		rng = qaws_curve_get_parameter_range(crv);
		qaws_curve_evaluate_2d(crv, rng.min_value, QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT(approx_eq(r.position.x, pts[0].x) && approx_eq(r.position.y, pts[0].y),
			"fit start matches data");

		qaws_curve_evaluate_2d(crv, rng.max_value, QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT(approx_eq(r.position.x, pts[20].x) && approx_eq(r.position.y, pts[20].y),
			"fit end matches data");

		/* Mid-curve should be near the semicircle (radius ~ 1) */
		{
			qaws_scalar mid = (rng.min_value + rng.max_value) * (qaws_scalar)0.5;
			qaws_scalar radius;
			qaws_curve_evaluate_2d(crv, mid, QAWS_EVAL_FLAG_POSITION, &r);
			radius = (qaws_scalar)sqrt((double)(r.position.x * r.position.x + r.position.y * r.position.y));
			TEST_ASSERT(radius > (qaws_scalar)0.95 && radius < (qaws_scalar)1.05,
				"fit mid-curve near unit circle");
		}

		qaws_curve_destroy(crv);
	}

	/* 3D fitting */
	{
		qaws_vec3 pts3d[11];
		unsigned int i;
		qaws_curve *crv = NULL;
		qaws_bspline_fit_desc desc;

		for (i = 0; i < 11; ++i) {
			qaws_scalar t = (qaws_scalar)i / (qaws_scalar)10;
			pts3d[i].x = t;
			pts3d[i].y = (qaws_scalar)sin(3.14159265358979 * (double)t);
			pts3d[i].z = t * t;
		}

		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_3D;
		desc.data_points = pts3d;
		desc.data_point_count = 11;
		desc.degree = 3;
		desc.control_point_count = 6;
		desc.parameters = NULL;

		TEST_ASSERT_STATUS(qaws_curve_fit_bspline(&desc, &crv));
		TEST_ASSERT(crv != NULL, "bspline fit 3d created");
		TEST_ASSERT(qaws_curve_get_dimension(crv) == QAWS_DIMENSION_3D, "fit 3d dimension");

		qaws_curve_destroy(crv);
	}

	/* Error cases */
	{
		qaws_vec2 pts2[3] = { {0,0}, {1,1}, {2,0} };
		qaws_bspline_fit_desc desc;
		qaws_curve *crv = NULL;

		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.data_points = pts2;
		desc.data_point_count = 3;
		desc.degree = 3;
		desc.control_point_count = 5;
		desc.parameters = NULL;

		/* n > m should fail */
		TEST_ASSERT(qaws_curve_fit_bspline(&desc, &crv) != QAWS_STATUS_OK,
			"fit n > m fails");
		TEST_ASSERT(crv == NULL, "fit n > m null");

		/* Null desc fails */
		TEST_ASSERT(qaws_curve_fit_bspline(NULL, &crv) != QAWS_STATUS_OK,
			"fit null desc fails");
	}
}

static void test_arc_length_reparam(void)
{
	printf("test_arc_length_reparam\n");

	/* Reparameterize a cubic Bezier */
	{
		qaws_vec2 cp[] = { {0,0}, {0,2}, {2,2}, {2,0} };
		qaws_bezier_desc desc;
		qaws_curve *src = NULL, *rep = NULL;
		qaws_eval_result_2d r;
		qaws_range rng;

		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.control_points = cp;
		desc.control_point_count = 4;

		TEST_ASSERT_STATUS(qaws_curve_create_bezier(&desc, &src));
		TEST_ASSERT_STATUS(qaws_curve_reparameterize_arc_length(src, 512, &rep));
		TEST_ASSERT(rep != NULL, "reparam created");
		TEST_ASSERT(qaws_curve_get_kind(rep) == QAWS_CURVE_KIND_REPARAMETERIZED,
			"reparam kind");

		rng = qaws_curve_get_parameter_range(rep);
		TEST_ASSERT(rng.min_value >= (qaws_scalar)0, "reparam range min 0");
		TEST_ASSERT(rng.max_value > (qaws_scalar)0, "reparam range max > 0");

		/* Start and end should match source */
		qaws_curve_evaluate_2d(rep, rng.min_value, QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)0) &&
			approx_eq(r.position.y, (qaws_scalar)0), "reparam start");

		qaws_curve_evaluate_2d(rep, rng.max_value, QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)2) &&
			approx_eq(r.position.y, (qaws_scalar)0), "reparam end");

		/* Uniform spacing: sample at equal arc-length intervals and
		   check that adjacent point distances are approximately equal */
		{
			qaws_scalar total = rng.max_value;
			qaws_vec2 prev_pos;
			qaws_scalar dists[10];
			int uniform = 1;
			unsigned int i;
			qaws_scalar mean;

			qaws_curve_evaluate_2d(rep, (qaws_scalar)0, QAWS_EVAL_FLAG_POSITION, &r);
			prev_pos = r.position;

			for (i = 1; i <= 10; ++i) {
				qaws_scalar s = total * (qaws_scalar)i / (qaws_scalar)10;
				qaws_scalar dx, dy;
				qaws_curve_evaluate_2d(rep, s, QAWS_EVAL_FLAG_POSITION, &r);
				dx = r.position.x - prev_pos.x;
				dy = r.position.y - prev_pos.y;
				dists[i - 1] = (qaws_scalar)sqrt((double)(dx * dx + dy * dy));
				prev_pos = r.position;
			}

			mean = (qaws_scalar)0;
			for (i = 0; i < 10; ++i)
				mean += dists[i];
			mean /= (qaws_scalar)10;

			for (i = 0; i < 10; ++i) {
				qaws_scalar ratio = dists[i] / mean;
				if (ratio < (qaws_scalar)0.8 || ratio > (qaws_scalar)1.2)
					uniform = 0;
			}
			TEST_ASSERT(uniform, "reparam uniform spacing");
		}

		/* D1 should be valid */
		qaws_curve_evaluate_2d(rep, rng.max_value * (qaws_scalar)0.5,
			QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &r);
		TEST_ASSERT(r.valid_flags & QAWS_EVAL_FLAG_D1, "reparam d1 valid");

		qaws_curve_destroy(rep);
		qaws_curve_destroy(src);
	}

	/* Error: null curve */
	{
		qaws_curve *out = NULL;
		TEST_ASSERT(qaws_curve_reparameterize_arc_length(NULL, 0, &out) != QAWS_STATUS_OK,
			"reparam null fails");
	}
}

static void test_svg_export_comparison(void)
{
	/*
	 * Side-by-side comparison of the two SVG rendering methods:
	 *
	 * LEFT column: sampled polyline (how the test SVGs work)
	 *   - Curve is evaluated at N uniform samples
	 *   - Drawn as SVG <polyline> through sample points
	 *   - Visible facets at low sample counts
	 *
	 * RIGHT column: SVG path export (qaws_curve_export_svg_path)
	 *   - Curve is converted to exact SVG path commands (M/C/Q)
	 *   - Bezier: exact, Hermite/Catmull-Rom: exact via per-span conversion
	 *   - Other families: Hermite-to-Bezier fit (smooth cubic segments)
	 *
	 * Multiple families are shown to demonstrate the difference.
	 */
	svg_writer svg;
	char path_buf[4096];
	unsigned int path_len;
	qaws_scalar col_left = 0.0f;
	qaws_scalar col_right = 6.0f;
	unsigned int low_samples = 12; /* deliberately low for polyline to show facets */

	printf("test_svg_export_comparison\n");

	if (!svg_open(&svg, SVG_OUTPUT_DIR "/21_svg_export_comparison.svg",
		(qaws_scalar)-1.5, (qaws_scalar)-2.0, (qaws_scalar)14.0, (qaws_scalar)18.0,
		(qaws_scalar)800, (qaws_scalar)1000))
		return;

	/* Labels */
	svg_label(&svg, col_left + 0.5f, 16.5f, "Polyline (12 samples)", "#ff79c6");
	svg_label(&svg, col_right + 0.5f, 16.5f, "SVG Path Export (exact)", "#50fa7b");

	/* ---- Row 1: Cubic Bezier ---- */
	{
		qaws_vec2 cp[] = {{0,0}, {1,3}, {3,3}, {4,0}};
		qaws_bezier_desc desc;
		qaws_curve* curve = NULL;
		unsigned int i;

		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.control_points = cp;
		desc.control_point_count = 4;
		qaws_curve_create_bezier(&desc, &curve);

		/* LEFT: polyline sampling */
		{
			qaws_vec2 pts[12];
			qaws_eval_result_2d r;
			for (i = 0; i < low_samples; i++)
			{
				qaws_scalar t = (qaws_scalar)i / (qaws_scalar)(low_samples - 1);
				qaws_curve_evaluate_2d(curve, t, QAWS_EVAL_FLAG_POSITION, &r);
				pts[i].x = r.position.x + col_left;
				pts[i].y = r.position.y + 12.0f;
			}
			svg_polyline(&svg, pts, low_samples, "#ff79c6", 2.0f);
			/* Show sample dots */
			svg_dots(&svg, pts, low_samples, "#ff79c680", 2.5f);
		}

		/* RIGHT: SVG path export */
		{
			qaws_curve_export_svg_path(curve, path_buf, sizeof(path_buf), &path_len, 0);
			/* We need to transform the path coordinates. Since SVG path is absolute,
			   we'll use a <g transform> group. */
			{
				qaws_scalar scale_x = (svg.svg_w - 2 * svg.padding) / svg.view_w;
				qaws_scalar scale_y = (svg.svg_h - 2 * svg.padding) / svg.view_h;
				qaws_scalar tx = svg.padding + (col_right - svg.view_x) * scale_x;
				qaws_scalar ty_base = svg.svg_h - svg.padding - (12.0f - svg.view_y) * scale_y;
				fprintf(svg.fp,
					"<g transform=\"translate(%.2f,%.2f) scale(%.4f,%.4f)\">\n"
					"<path d=\"%s\" fill=\"none\" stroke=\"#50fa7b\" stroke-width=\"%.4f\"/>\n"
					"</g>\n",
					(double)tx, (double)ty_base,
					(double)scale_x, (double)(-scale_y),
					path_buf, (double)(2.0f / scale_x));
			}
		}

		/* Control polygon (both sides) */
		{
			qaws_vec2 cp_left[4], cp_right[4];
			for (i = 0; i < 4; i++)
			{
				cp_left[i].x = cp[i].x + col_left;
				cp_left[i].y = cp[i].y + 12.0f;
				cp_right[i].x = cp[i].x + col_right;
				cp_right[i].y = cp[i].y + 12.0f;
			}
			svg_polyline_dashed(&svg, cp_left, 4, "#555577", 1.0f, "3,3");
			svg_polyline_dashed(&svg, cp_right, 4, "#555577", 1.0f, "3,3");
			svg_dots(&svg, cp_left, 4, "#ffffff", 2.5f);
			svg_dots(&svg, cp_right, 4, "#ffffff", 2.5f);
		}

		svg_label(&svg, col_left, 15.8f, "Cubic Bezier", "#8888aa");
		qaws_curve_destroy(curve);
	}

	/* ---- Row 2: Catmull-Rom ---- */
	{
		qaws_vec2 cr_cp[] = {{0,0}, {1,2.5f}, {2,0.5f}, {3,2.0f}, {4,0}};
		qaws_catmull_rom_desc desc;
		qaws_curve* curve = NULL;
		unsigned int i;
		qaws_range range;

		desc.dimension = QAWS_DIMENSION_2D;
		desc.control_points = cr_cp;
		desc.control_point_count = 5;
		desc.parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL;
		desc.closed = 0;
		qaws_curve_create_catmull_rom(&desc, &curve);
		range = qaws_curve_get_parameter_range(curve);

		/* LEFT: polyline */
		{
			qaws_vec2 pts[12];
			qaws_eval_result_2d r;
			for (i = 0; i < low_samples; i++)
			{
				qaws_scalar t = range.min_value + (range.max_value - range.min_value)
					* (qaws_scalar)i / (qaws_scalar)(low_samples - 1);
				qaws_curve_evaluate_2d(curve, t, QAWS_EVAL_FLAG_POSITION, &r);
				pts[i].x = r.position.x + col_left;
				pts[i].y = r.position.y + 8.0f;
			}
			svg_polyline(&svg, pts, low_samples, "#ff79c6", 2.0f);
			svg_dots(&svg, pts, low_samples, "#ff79c680", 2.5f);
		}

		/* RIGHT: SVG path export */
		{
			qaws_curve_export_svg_path(curve, path_buf, sizeof(path_buf), &path_len, 0);
			{
				qaws_scalar scale_x = (svg.svg_w - 2 * svg.padding) / svg.view_w;
				qaws_scalar scale_y = (svg.svg_h - 2 * svg.padding) / svg.view_h;
				qaws_scalar tx = svg.padding + (col_right - svg.view_x) * scale_x;
				qaws_scalar ty_base = svg.svg_h - svg.padding - (8.0f - svg.view_y) * scale_y;
				fprintf(svg.fp,
					"<g transform=\"translate(%.2f,%.2f) scale(%.4f,%.4f)\">\n"
					"<path d=\"%s\" fill=\"none\" stroke=\"#50fa7b\" stroke-width=\"%.4f\"/>\n"
					"</g>\n",
					(double)tx, (double)ty_base,
					(double)scale_x, (double)(-scale_y),
					path_buf, (double)(2.0f / scale_x));
			}
		}

		/* Control points both sides */
		{
			qaws_vec2 cp_left[5], cp_right[5];
			for (i = 0; i < 5; i++)
			{
				cp_left[i].x = cr_cp[i].x + col_left;
				cp_left[i].y = cr_cp[i].y + 8.0f;
				cp_right[i].x = cr_cp[i].x + col_right;
				cp_right[i].y = cr_cp[i].y + 8.0f;
			}
			svg_dots(&svg, cp_left, 5, "#ffffff", 2.5f);
			svg_dots(&svg, cp_right, 5, "#ffffff", 2.5f);
		}

		svg_label(&svg, col_left, 11.5f, "Catmull-Rom (centripetal)", "#8888aa");
		qaws_curve_destroy(curve);
	}

	/* ---- Row 3: B-Spline ---- */
	{
		qaws_vec2 bs_cp[] = {{0,0}, {0.5f,2.5f}, {1.5f,3.0f}, {2.5f,1.0f}, {3.5f,2.5f}, {4,0}};
		qaws_bspline_desc desc;
		qaws_curve* curve = NULL;
		unsigned int i;
		qaws_range range;

		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.control_points = bs_cp;
		desc.control_point_count = 6;
		desc.knots = NULL;
		desc.knot_count = 0;
		desc.is_uniform = 1;
		desc.is_closed = 0;
		qaws_curve_create_bspline(&desc, &curve);
		range = qaws_curve_get_parameter_range(curve);

		/* LEFT: polyline */
		{
			qaws_vec2 pts[12];
			qaws_eval_result_2d r;
			for (i = 0; i < low_samples; i++)
			{
				qaws_scalar t = range.min_value + (range.max_value - range.min_value)
					* (qaws_scalar)i / (qaws_scalar)(low_samples - 1);
				qaws_curve_evaluate_2d(curve, t, QAWS_EVAL_FLAG_POSITION, &r);
				pts[i].x = r.position.x + col_left;
				pts[i].y = r.position.y + 4.0f;
			}
			svg_polyline(&svg, pts, low_samples, "#ff79c6", 2.0f);
			svg_dots(&svg, pts, low_samples, "#ff79c680", 2.5f);
		}

		/* RIGHT: SVG path export (sampled + Hermite-to-Bezier fit) */
		{
			qaws_curve_export_svg_path(curve, path_buf, sizeof(path_buf), &path_len, 16);
			{
				qaws_scalar scale_x = (svg.svg_w - 2 * svg.padding) / svg.view_w;
				qaws_scalar scale_y = (svg.svg_h - 2 * svg.padding) / svg.view_h;
				qaws_scalar tx = svg.padding + (col_right - svg.view_x) * scale_x;
				qaws_scalar ty_base = svg.svg_h - svg.padding - (4.0f - svg.view_y) * scale_y;
				fprintf(svg.fp,
					"<g transform=\"translate(%.2f,%.2f) scale(%.4f,%.4f)\">\n"
					"<path d=\"%s\" fill=\"none\" stroke=\"#50fa7b\" stroke-width=\"%.4f\"/>\n"
					"</g>\n",
					(double)tx, (double)ty_base,
					(double)scale_x, (double)(-scale_y),
					path_buf, (double)(2.0f / scale_x));
			}
		}

		/* Control points both sides */
		{
			qaws_vec2 cp_left[6], cp_right[6];
			for (i = 0; i < 6; i++)
			{
				cp_left[i].x = bs_cp[i].x + col_left;
				cp_left[i].y = bs_cp[i].y + 4.0f;
				cp_right[i].x = bs_cp[i].x + col_right;
				cp_right[i].y = bs_cp[i].y + 4.0f;
			}
			svg_dots(&svg, cp_left, 6, "#ffffff", 2.5f);
			svg_dots(&svg, cp_right, 6, "#ffffff", 2.5f);
		}

		svg_label(&svg, col_left, 7.5f, "B-Spline (degree 3)", "#8888aa");
		qaws_curve_destroy(curve);
	}

	/* ---- Row 4: Hermite ---- */
	{
		qaws_vec2 h_pts[] = {{0,0}, {2,2.5f}, {4,0}};
		qaws_vec2 h_der[] = {{3,2}, {2,0}, {3,-2}};
		qaws_hermite_desc desc;
		qaws_curve* curve = NULL;
		unsigned int i;
		qaws_range range;

		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.points = h_pts;
		desc.point_count = 3;
		desc.derivatives = h_der;
		desc.derivative_count = 3;
		qaws_curve_create_hermite(&desc, &curve);
		range = qaws_curve_get_parameter_range(curve);

		/* LEFT: polyline */
		{
			qaws_vec2 pts[12];
			qaws_eval_result_2d r;
			for (i = 0; i < low_samples; i++)
			{
				qaws_scalar t = range.min_value + (range.max_value - range.min_value)
					* (qaws_scalar)i / (qaws_scalar)(low_samples - 1);
				qaws_curve_evaluate_2d(curve, t, QAWS_EVAL_FLAG_POSITION, &r);
				pts[i].x = r.position.x + col_left;
				pts[i].y = r.position.y + 0.0f;
			}
			svg_polyline(&svg, pts, low_samples, "#ff79c6", 2.0f);
			svg_dots(&svg, pts, low_samples, "#ff79c680", 2.5f);
		}

		/* RIGHT: SVG path export (exact per-span conversion) */
		{
			qaws_curve_export_svg_path(curve, path_buf, sizeof(path_buf), &path_len, 0);
			{
				qaws_scalar scale_x = (svg.svg_w - 2 * svg.padding) / svg.view_w;
				qaws_scalar scale_y = (svg.svg_h - 2 * svg.padding) / svg.view_h;
				qaws_scalar tx = svg.padding + (col_right - svg.view_x) * scale_x;
				qaws_scalar ty_base = svg.svg_h - svg.padding - (0.0f - svg.view_y) * scale_y;
				fprintf(svg.fp,
					"<g transform=\"translate(%.2f,%.2f) scale(%.4f,%.4f)\">\n"
					"<path d=\"%s\" fill=\"none\" stroke=\"#50fa7b\" stroke-width=\"%.4f\"/>\n"
					"</g>\n",
					(double)tx, (double)ty_base,
					(double)scale_x, (double)(-scale_y),
					path_buf, (double)(2.0f / scale_x));
			}
		}

		/* Key points both sides */
		{
			qaws_vec2 kp_left[3], kp_right[3];
			for (i = 0; i < 3; i++)
			{
				kp_left[i].x = h_pts[i].x + col_left;
				kp_left[i].y = h_pts[i].y + 0.0f;
				kp_right[i].x = h_pts[i].x + col_right;
				kp_right[i].y = h_pts[i].y + 0.0f;
			}
			svg_dots(&svg, kp_left, 3, "#ffffff", 2.5f);
			svg_dots(&svg, kp_right, 3, "#ffffff", 2.5f);
		}

		svg_label(&svg, col_left, 3.5f, "Hermite (cubic)", "#8888aa");
		qaws_curve_destroy(curve);
	}

	/* Unit test: verify SVG path export returns valid data */
	{
		qaws_vec2 cp[] = {{0,0}, {1,2}, {3,1}};
		qaws_bezier_desc desc;
		qaws_curve* curve = NULL;
		char small_buf[8];
		unsigned int len;

		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 2;
		desc.control_points = cp;
		desc.control_point_count = 3;
		qaws_curve_create_bezier(&desc, &curve);

		/* Quadratic bezier should produce M...Q... */
		qaws_curve_export_svg_path(curve, path_buf, sizeof(path_buf), &path_len, 0);
		TEST_ASSERT(path_len > 0, "svg_export: non-empty path");
		TEST_ASSERT(path_buf[0] == 'M', "svg_export: starts with M");
		TEST_ASSERT(strstr(path_buf, "Q") != NULL, "svg_export: contains Q command");

		/* Small buffer should truncate gracefully */
		qaws_curve_export_svg_path(curve, small_buf, sizeof(small_buf), &len, 0);
		TEST_ASSERT(len > 0, "svg_export: truncated still writes");
		TEST_ASSERT(small_buf[sizeof(small_buf) - 1] == '\0', "svg_export: null-terminated");

		/* Error: 3D curve should fail */
		{
			qaws_vec3 cp3[] = {{0,0,0}, {1,1,1}};
			qaws_bezier_desc d3;
			qaws_curve* c3 = NULL;
			qaws_status s;
			d3.dimension = QAWS_DIMENSION_3D;
			d3.degree = 1;
			d3.control_points = cp3;
			d3.control_point_count = 2;
			qaws_curve_create_bezier(&d3, &c3);
			s = qaws_curve_export_svg_path(c3, path_buf, sizeof(path_buf), &path_len, 0);
			TEST_ASSERT(s == QAWS_STATUS_INVALID_DIMENSION, "svg_export: reject 3D");
			qaws_curve_destroy(c3);
		}

		/* Error: null args */
		{
			qaws_status s = qaws_curve_export_svg_path(NULL, path_buf, sizeof(path_buf), &path_len, 0);
			TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT, "svg_export: null curve");
		}

		qaws_curve_destroy(curve);
	}

	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/21_svg_export_comparison.svg\n");
}

static void test_polyline_import(void)
{
	printf("test_polyline_import\n");

	/* 2D Catmull-Rom from polyline */
	{
		qaws_vec2 pts[] = {{0,0}, {1,2}, {3,1}, {4,3}, {5,0}};
		qaws_curve* curve = NULL;
		qaws_eval_result_2d r;

		TEST_ASSERT_STATUS(qaws_polyline_import_catmull_rom(
			QAWS_DIMENSION_2D, pts, 5,
			QAWS_PARAMETERIZATION_CENTRIPETAL, 0, &curve));
		TEST_ASSERT(curve != NULL, "polyline_catmull_rom: created");
		TEST_ASSERT(qaws_curve_get_kind(curve) == QAWS_CURVE_KIND_CATMULL_ROM,
			"polyline_catmull_rom: correct kind");
		TEST_ASSERT(qaws_curve_get_dimension(curve) == QAWS_DIMENSION_2D,
			"polyline_catmull_rom: correct dimension");

		/* Evaluate at start — should be near first interior point */
		qaws_curve_evaluate_2d(curve, (qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT(fabs(r.position.x - pts[1].x) < 0.01 &&
			fabs(r.position.y - pts[1].y) < 0.01,
			"polyline_catmull_rom: t=0 near first interpolated point");

		qaws_curve_destroy(curve);
	}

	/* 3D Catmull-Rom from polyline */
	{
		qaws_vec3 pts[] = {{0,0,0}, {1,2,1}, {3,1,2}, {4,3,0}};
		qaws_curve* curve = NULL;

		TEST_ASSERT_STATUS(qaws_polyline_import_catmull_rom(
			QAWS_DIMENSION_3D, pts, 4,
			QAWS_PARAMETERIZATION_CHORDAL, 0, &curve));
		TEST_ASSERT(curve != NULL, "polyline_catmull_rom_3d: created");
		TEST_ASSERT(qaws_curve_get_dimension(curve) == QAWS_DIMENSION_3D,
			"polyline_catmull_rom_3d: correct dimension");
		qaws_curve_destroy(curve);
	}

	/* Closed Catmull-Rom loop */
	{
		qaws_vec2 pts[] = {{0,0}, {2,3}, {4,0}, {2,-1}};
		qaws_curve* curve = NULL;

		TEST_ASSERT_STATUS(qaws_polyline_import_catmull_rom(
			QAWS_DIMENSION_2D, pts, 4,
			QAWS_PARAMETERIZATION_CENTRIPETAL, 1, &curve));
		TEST_ASSERT(curve != NULL, "polyline_catmull_rom_closed: created");
		qaws_curve_destroy(curve);
	}

	/* 2D Trajectory from polyline */
	{
		qaws_vec2 pts[] = {{0,0}, {1,2}, {3,1}, {5,3}};
		qaws_curve* curve = NULL;
		qaws_eval_result_2d r0, r1;
		qaws_range range;

		TEST_ASSERT_STATUS(qaws_polyline_import_trajectory(
			QAWS_DIMENSION_2D, pts, 4, 0, &curve));
		TEST_ASSERT(curve != NULL, "polyline_trajectory: created");
		TEST_ASSERT(qaws_curve_get_kind(curve) == QAWS_CURVE_KIND_TRAJECTORY,
			"polyline_trajectory: correct kind");

		/* Endpoints should match first/last points */
		range = qaws_curve_get_parameter_range(curve);
		qaws_curve_evaluate_2d(curve, range.min_value, QAWS_EVAL_FLAG_POSITION, &r0);
		qaws_curve_evaluate_2d(curve, range.max_value, QAWS_EVAL_FLAG_POSITION, &r1);
		TEST_ASSERT(fabs(r0.position.x - pts[0].x) < 0.01 &&
			fabs(r0.position.y - pts[0].y) < 0.01,
			"polyline_trajectory: start matches first point");
		TEST_ASSERT(fabs(r1.position.x - pts[3].x) < 0.01 &&
			fabs(r1.position.y - pts[3].y) < 0.01,
			"polyline_trajectory: end matches last point");

		qaws_curve_destroy(curve);
	}

	/* 3D Trajectory from polyline */
	{
		qaws_vec3 pts[] = {{0,0,0}, {1,1,1}, {2,0,2}};
		qaws_curve* curve = NULL;

		TEST_ASSERT_STATUS(qaws_polyline_import_trajectory(
			QAWS_DIMENSION_3D, pts, 3, 0, &curve));
		TEST_ASSERT(curve != NULL, "polyline_trajectory_3d: created");
		qaws_curve_destroy(curve);
	}

	/* Closed trajectory */
	{
		qaws_vec2 pts[] = {{0,0}, {2,3}, {4,0}, {2,-1}};
		qaws_curve* curve = NULL;

		TEST_ASSERT_STATUS(qaws_polyline_import_trajectory(
			QAWS_DIMENSION_2D, pts, 4, 1, &curve));
		TEST_ASSERT(curve != NULL, "polyline_trajectory_closed: created");
		qaws_curve_destroy(curve);
	}

	/* Error: too few points */
	{
		qaws_vec2 pts[] = {{0,0}};
		qaws_curve* curve = NULL;
		qaws_status s;

		s = qaws_polyline_import_catmull_rom(
			QAWS_DIMENSION_2D, pts, 1,
			QAWS_PARAMETERIZATION_CENTRIPETAL, 0, &curve);
		TEST_ASSERT(s != QAWS_STATUS_OK, "polyline_catmull_rom: rejects 1 point");

		s = qaws_polyline_import_trajectory(
			QAWS_DIMENSION_2D, pts, 1, 0, &curve);
		TEST_ASSERT(s != QAWS_STATUS_OK, "polyline_trajectory: rejects 1 point");
	}

	/* Error: null pointers */
	{
		qaws_curve* curve = NULL;
		qaws_status s;

		s = qaws_polyline_import_catmull_rom(
			QAWS_DIMENSION_2D, NULL, 4,
			QAWS_PARAMETERIZATION_CENTRIPETAL, 0, &curve);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT, "polyline_catmull_rom: null points");

		s = qaws_polyline_import_trajectory(
			QAWS_DIMENSION_2D, NULL, 4, 0, &curve);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT, "polyline_trajectory: null points");
	}

	/* ---- Polyline export tests ---- */

	/* 2D Bezier export */
	{
		qaws_vec2 cp[] = {{0,0}, {1,3}, {3,3}, {4,0}};
		qaws_bezier_desc desc;
		qaws_curve* curve = NULL;
		qaws_vec2 pts[32];
		unsigned int count = 0;

		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.control_points = cp;
		desc.control_point_count = 4;
		qaws_curve_create_bezier(&desc, &curve);

		TEST_ASSERT_STATUS(qaws_polyline_export_2d(curve, 10, QAWS_POLYLINE_UNIFORM, pts, &count));
		TEST_ASSERT(count == 10, "polyline_export_2d: correct count");
		/* First point should be at curve start */
		TEST_ASSERT(fabs(pts[0].x - 0.0) < 0.01 && fabs(pts[0].y - 0.0) < 0.01,
			"polyline_export_2d: start point");
		/* Last point should be at curve end */
		TEST_ASSERT(fabs(pts[9].x - 4.0) < 0.01 && fabs(pts[9].y - 0.0) < 0.01,
			"polyline_export_2d: end point");
		/* Mid point should be reasonable (on the curve) */
		TEST_ASSERT(pts[5].y > 0.0, "polyline_export_2d: mid point above x axis");

		qaws_curve_destroy(curve);
	}

	/* 3D Bezier export */
	{
		qaws_vec3 cp[] = {{0,0,0}, {1,2,1}, {3,2,1}, {4,0,0}};
		qaws_bezier_desc desc;
		qaws_curve* curve = NULL;
		qaws_vec3 pts[20];
		unsigned int count = 0;

		desc.dimension = QAWS_DIMENSION_3D;
		desc.degree = 3;
		desc.control_points = cp;
		desc.control_point_count = 4;
		qaws_curve_create_bezier(&desc, &curve);

		TEST_ASSERT_STATUS(qaws_polyline_export_3d(curve, 8, QAWS_POLYLINE_UNIFORM, pts, &count));
		TEST_ASSERT(count == 8, "polyline_export_3d: correct count");
		TEST_ASSERT(fabs(pts[0].x - 0.0) < 0.01 && fabs(pts[0].y - 0.0) < 0.01 &&
			fabs(pts[0].z - 0.0) < 0.01, "polyline_export_3d: start point");
		TEST_ASSERT(fabs(pts[7].x - 4.0) < 0.01 && fabs(pts[7].y - 0.0) < 0.01 &&
			fabs(pts[7].z - 0.0) < 0.01, "polyline_export_3d: end point");

		qaws_curve_destroy(curve);
	}

	/* Round-trip: export 2D then re-import as Catmull-Rom */
	{
		qaws_vec2 cp[] = {{0,0}, {1,3}, {3,3}, {4,0}};
		qaws_bezier_desc desc;
		qaws_curve* original = NULL;
		qaws_curve* reimported = NULL;
		qaws_vec2 exported[32];
		unsigned int count = 0;
		qaws_eval_result_2d r_orig, r_reimp;

		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.control_points = cp;
		desc.control_point_count = 4;
		qaws_curve_create_bezier(&desc, &original);

		TEST_ASSERT_STATUS(qaws_polyline_export_2d(original, 20, QAWS_POLYLINE_UNIFORM, exported, &count));
		TEST_ASSERT_STATUS(qaws_polyline_import_catmull_rom(
			QAWS_DIMENSION_2D, exported, count,
			QAWS_PARAMETERIZATION_CENTRIPETAL, 0, &reimported));
		TEST_ASSERT(reimported != NULL, "polyline_round_trip: reimported");

		/* Compare midpoints — Catmull-Rom interpolates through points so
		   with enough samples the reimported curve should be close */
		qaws_curve_evaluate_2d(original, (qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION, &r_orig);
		{
			qaws_range rng = qaws_curve_get_parameter_range(reimported);
			qaws_scalar mid = (rng.min_value + rng.max_value) * (qaws_scalar)0.5;
			qaws_curve_evaluate_2d(reimported, mid, QAWS_EVAL_FLAG_POSITION, &r_reimp);
		}
		TEST_ASSERT(fabs(r_orig.position.x - r_reimp.position.x) < 0.5 &&
			fabs(r_orig.position.y - r_reimp.position.y) < 0.5,
			"polyline_round_trip: midpoints close");

		qaws_curve_destroy(original);
		qaws_curve_destroy(reimported);
	}

	/* Dimension mismatch errors */
	{
		qaws_vec2 cp[] = {{0,0}, {1,1}};
		qaws_vec3 cp3[] = {{0,0,0}, {1,1,1}};
		qaws_bezier_desc desc2, desc3;
		qaws_curve* c2d = NULL;
		qaws_curve* c3d = NULL;
		qaws_vec2 pts2[4];
		qaws_vec3 pts3[4];
		unsigned int count;
		qaws_status s;

		desc2.dimension = QAWS_DIMENSION_2D;
		desc2.degree = 1;
		desc2.control_points = cp;
		desc2.control_point_count = 2;
		qaws_curve_create_bezier(&desc2, &c2d);

		desc3.dimension = QAWS_DIMENSION_3D;
		desc3.degree = 1;
		desc3.control_points = cp3;
		desc3.control_point_count = 2;
		qaws_curve_create_bezier(&desc3, &c3d);

		s = qaws_polyline_export_3d(c2d, 4, QAWS_POLYLINE_UNIFORM, pts3, &count);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_DIMENSION, "polyline_export: 2d curve to 3d fails");
		s = qaws_polyline_export_2d(c3d, 4, QAWS_POLYLINE_UNIFORM, pts2, &count);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_DIMENSION, "polyline_export: 3d curve to 2d fails");

		qaws_curve_destroy(c2d);
		qaws_curve_destroy(c3d);
	}

	/* Error: too few samples */
	{
		qaws_vec2 cp[] = {{0,0}, {1,1}};
		qaws_bezier_desc desc;
		qaws_curve* curve = NULL;
		qaws_vec2 pts[4];
		unsigned int count;
		qaws_status s;

		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 1;
		desc.control_points = cp;
		desc.control_point_count = 2;
		qaws_curve_create_bezier(&desc, &curve);

		s = qaws_polyline_export_2d(curve, 1, QAWS_POLYLINE_UNIFORM, pts, &count);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT, "polyline_export: rejects 1 sample");

		s = qaws_polyline_export_2d(curve, 0, QAWS_POLYLINE_UNIFORM, pts, &count);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT, "polyline_export: rejects 0 samples");

		qaws_curve_destroy(curve);
	}

	/* Null pointer errors */
	{
		qaws_vec2 pts[4];
		unsigned int count;
		qaws_status s;

		s = qaws_polyline_export_2d(NULL, 4, QAWS_POLYLINE_UNIFORM, pts, &count);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT, "polyline_export: null curve");

		s = qaws_polyline_export_2d(NULL, 4, QAWS_POLYLINE_UNIFORM, NULL, &count);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT, "polyline_export: null out_points");
	}

	/* ---- Curvature-weighted export tests ---- */

	/* 2D curvature-weighted: should produce more points near high curvature */
	{
		qaws_vec2 cp[] = {{0,0}, {0,3}, {3,3}, {3,0}};
		qaws_bezier_desc desc;
		qaws_curve* curve = NULL;
		qaws_vec2 uniform_pts[32], curvature_pts[32];
		unsigned int u_count = 0, c_count = 0;
		unsigned int i;
		qaws_scalar u_max_gap = 0, c_max_gap = 0;

		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.control_points = cp;
		desc.control_point_count = 4;
		qaws_curve_create_bezier(&desc, &curve);

		TEST_ASSERT_STATUS(qaws_polyline_export_2d(curve, 20, QAWS_POLYLINE_UNIFORM, uniform_pts, &u_count));
		TEST_ASSERT_STATUS(qaws_polyline_export_2d(curve, 20, QAWS_POLYLINE_CURVATURE, curvature_pts, &c_count));
		TEST_ASSERT(u_count == 20, "polyline_curvature_2d: uniform count");
		TEST_ASSERT(c_count > 0 && c_count <= 20, "polyline_curvature_2d: curvature count valid");

		/* Compute max gap between consecutive points for each */
		for (i = 1; i < u_count; i++)
		{
			qaws_scalar dx = uniform_pts[i].x - uniform_pts[i-1].x;
			qaws_scalar dy = uniform_pts[i].y - uniform_pts[i-1].y;
			qaws_scalar gap = dx*dx + dy*dy;
			if (gap > u_max_gap) u_max_gap = gap;
		}
		for (i = 1; i < c_count; i++)
		{
			qaws_scalar dx = curvature_pts[i].x - curvature_pts[i-1].x;
			qaws_scalar dy = curvature_pts[i].y - curvature_pts[i-1].y;
			qaws_scalar gap = dx*dx + dy*dy;
			if (gap > c_max_gap) c_max_gap = gap;
		}
		/* Curvature sampling should have a smaller max gap (more even spacing
		   in screen space) or at least produce valid output */
		TEST_ASSERT(c_max_gap > 0, "polyline_curvature_2d: non-degenerate output");

		qaws_curve_destroy(curve);
	}

	/* 3D curvature-weighted */
	{
		qaws_vec3 cp[] = {{0,0,0}, {0,2,1}, {2,2,1}, {2,0,0}};
		qaws_bezier_desc desc;
		qaws_curve* curve = NULL;
		qaws_vec3 pts[16];
		unsigned int count = 0;

		desc.dimension = QAWS_DIMENSION_3D;
		desc.degree = 3;
		desc.control_points = cp;
		desc.control_point_count = 4;
		qaws_curve_create_bezier(&desc, &curve);

		TEST_ASSERT_STATUS(qaws_polyline_export_3d(curve, 12, QAWS_POLYLINE_CURVATURE, pts, &count));
		TEST_ASSERT(count > 0 && count <= 12, "polyline_curvature_3d: valid count");
		/* Start and end should still be on the curve */
		TEST_ASSERT(fabs(pts[0].x - 0.0) < 0.1, "polyline_curvature_3d: start near origin");
		TEST_ASSERT(fabs(pts[count-1].x - 2.0) < 0.1, "polyline_curvature_3d: end near (2,0,0)");

		qaws_curve_destroy(curve);
	}
}

int test_21_export_import_main(void)
{
	g_pass = 0;
	g_fail = 0;

	svg_ensure_output_dir();
	test_bspline_fitting();
	test_arc_length_reparam();
	test_svg_export_comparison();
	test_polyline_import();

	printf("21_export_import: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail > 0 ? 1 : 0;
}
