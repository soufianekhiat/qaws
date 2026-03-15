#include "test_common.h"

static void test_obj_surface_bezier(void)
{
	obj_writer w;
	qaws_surface* surf = NULL;
	qaws_surface_bezier_desc desc;
	/* Bicubic Bezier patch (4x4 CPs) shaped like a saddle */
	qaws_vec3 cp[16];
	unsigned int i, j;

	printf("test_obj_surface_bezier\n");

	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++)
		{
			unsigned int idx = i * 4 + j;
			qaws_scalar x = (qaws_scalar)i - (qaws_scalar)1.5;
			qaws_scalar y = (qaws_scalar)j - (qaws_scalar)1.5;
			cp[idx].x = x;
			cp[idx].y = y;
			cp[idx].z = (qaws_scalar)(0.5 * (x * x - y * y)); /* saddle */
		}
	}

	desc.u_degree = 3;
	desc.v_degree = 3;
	desc.control_points = cp;
	desc.u_point_count = 4;
	desc.v_point_count = 4;

	if (qaws_surface_create_bezier(&desc, &surf) != QAWS_STATUS_OK) return;

	if (!obj_open(&w,
		OBJ_OUTPUT_DIR "/30_surface_bezier.obj",
		OBJ_OUTPUT_DIR "/30_surface_bezier.mtl"))
	{
		qaws_surface_destroy(surf);
		return;
	}

	obj_material(&w, "patch", 0.2, 0.6, 1.0);
	obj_material(&w, "cp_grid", 1.0, 0.3, 0.1);

	obj_group(&w, "bezier_patch");
	obj_use_material(&w, "patch");
	obj_surface_mesh(&w, surf, 32, 32);

	obj_group(&w, "control_points");
	obj_use_material(&w, "cp_grid");
	obj_surface_cp_grid(&w, cp, 4, 4, (qaws_scalar)0.06);

	obj_close(&w);
	qaws_surface_destroy(surf);
	printf("  -> " OBJ_OUTPUT_DIR "/30_surface_bezier.obj\n");
}

static void test_obj_surface_bspline(void)
{
	obj_writer w;
	qaws_surface* surf = NULL;
	qaws_surface_bspline_desc desc;
	/* 6x6 control points, degree (3,3) - rolling terrain */
	qaws_vec3 cp[36];
	unsigned int i, j;

	printf("test_obj_surface_bspline\n");

	for (i = 0; i < 6; i++)
	{
		for (j = 0; j < 6; j++)
		{
			unsigned int idx = i * 6 + j;
			qaws_scalar x = (qaws_scalar)i;
			qaws_scalar y = (qaws_scalar)j;
			cp[idx].x = x;
			cp[idx].y = y;
			cp[idx].z = (qaws_scalar)(sin(x * 1.2) * cos(y * 0.8) * 0.8);
		}
	}

	memset(&desc, 0, sizeof(desc));
	desc.u_degree = 3;
	desc.v_degree = 3;
	desc.control_points = cp;
	desc.u_point_count = 6;
	desc.v_point_count = 6;
	desc.u_knot_count = 0;
	desc.v_knot_count = 0;

	if (qaws_surface_create_bspline(&desc, &surf) != QAWS_STATUS_OK) return;

	if (!obj_open(&w,
		OBJ_OUTPUT_DIR "/30_surface_bspline.obj",
		OBJ_OUTPUT_DIR "/30_surface_bspline.mtl"))
	{
		qaws_surface_destroy(surf);
		return;
	}

	obj_material(&w, "terrain", 0.3, 0.8, 0.3);
	obj_material(&w, "cp_grid", 0.8, 0.4, 0.1);

	obj_group(&w, "bspline_surface");
	obj_use_material(&w, "terrain");
	obj_surface_mesh(&w, surf, 40, 40);

	obj_group(&w, "control_points");
	obj_use_material(&w, "cp_grid");
	obj_surface_cp_grid(&w, cp, 6, 6, (qaws_scalar)0.05);

	obj_close(&w);
	qaws_surface_destroy(surf);
	printf("  -> " OBJ_OUTPUT_DIR "/30_surface_bspline.obj\n");
}

static void test_obj_surface_nurbs(void)
{
	obj_writer w;
	qaws_surface* surf = NULL;
	qaws_surface_nurbs_desc desc;
	qaws_scalar wm = (qaws_scalar)0.707106781; /* 1/sqrt(2) */
	qaws_scalar R = 2;
	unsigned int row, col;

	/* Full hemisphere: quarter-circle profile (3 CPs in u) revolved 360 degrees
	   (9 CPs in v). Degree (2,2). */
	qaws_vec3 cp[27]; /* 3 x 9 */
	qaws_scalar weights[27];
	qaws_scalar u_knots[6] = {0, 0, 0, 1, 1, 1};
	/* 9 CPs, degree 2 -> 12 knots, double knots at each 90-degree break */
	qaws_scalar v_knots[12] = {
		0, 0, 0,
		0.25f, 0.25f,
		0.5f, 0.5f,
		0.75f, 0.75f,
		1, 1, 1
	};

	printf("test_obj_surface_nurbs\n");

	/* Profile: quarter circle in XZ plane from (R,0,0) up to (0,0,R).
	   3 rows: P0=(R,0,0) w=1, P1=(R,0,R) w=1/sqrt2, P2=(0,0,R) w=1 */
	{
		qaws_scalar profile_r[3] = {R, R, 0}; /* radial distance from Z axis */
		qaws_scalar profile_z[3] = {0, R, R};
		qaws_scalar profile_w[3] = {1, wm, 1};

		/* 9 revolution CPs at 0, 45, 90, 135, 180, 225, 270, 315, 360 degrees */
		double angles[9] = {0, 45, 90, 135, 180, 225, 270, 315, 360};
		qaws_scalar rev_w[9] = {1, wm, 1, wm, 1, wm, 1, wm, 1};

		for (row = 0; row < 3; row++)
		{
			qaws_scalar r_val = profile_r[row];
			qaws_scalar z_val = profile_z[row];
			qaws_scalar pw = profile_w[row];
			for (col = 0; col < 9; col++)
			{
				double a = angles[col] * 3.14159265358979323846 / 180.0;
				unsigned int idx = row * 9 + col;
				cp[idx].x = (qaws_scalar)(r_val * cos(a));
				cp[idx].y = (qaws_scalar)(r_val * sin(a));
				cp[idx].z = z_val;
				weights[idx] = pw * rev_w[col];
			}
		}
	}

	memset(&desc, 0, sizeof(desc));
	desc.u_degree = 2;
	desc.v_degree = 2;
	desc.control_points = cp;
	desc.u_point_count = 3;
	desc.v_point_count = 9;
	desc.weights = weights;
	desc.u_knots = u_knots;
	desc.u_knot_count = 6;
	desc.v_knots = v_knots;
	desc.v_knot_count = 12;

	if (qaws_surface_create_nurbs(&desc, &surf) != QAWS_STATUS_OK)
	{
		printf("  SKIP: nurbs surface create failed\n");
		return;
	}

	if (!obj_open(&w,
		OBJ_OUTPUT_DIR "/30_surface_nurbs.obj",
		OBJ_OUTPUT_DIR "/30_surface_nurbs.mtl"))
	{
		qaws_surface_destroy(surf);
		return;
	}

	obj_material(&w, "shell", 0.8, 0.2, 0.6);
	obj_material(&w, "cp_grid", 0.3, 0.3, 0.9);

	obj_group(&w, "nurbs_surface");
	obj_use_material(&w, "shell");
	obj_surface_mesh(&w, surf, 32, 32);

	obj_group(&w, "control_points");
	obj_use_material(&w, "cp_grid");
	obj_surface_cp_grid(&w, cp, 3, 9, (qaws_scalar)0.08);

	obj_close(&w);
	qaws_surface_destroy(surf);
	printf("  -> " OBJ_OUTPUT_DIR "/30_surface_nurbs.obj\n");
}

static void test_obj_surface_gallery(void)
{
	/* All three surface types in one OBJ file for comparison */
	obj_writer w;
	qaws_surface* bezier_surf = NULL;
	qaws_surface* bspline_surf = NULL;
	qaws_surface* nurbs_surf = NULL;
	unsigned int i, j;

	printf("test_obj_surface_gallery\n");

	if (!obj_open(&w,
		OBJ_OUTPUT_DIR "/30_surface_gallery.obj",
		OBJ_OUTPUT_DIR "/30_surface_gallery.mtl"))
		return;

	obj_material(&w, "bezier", 0.2, 0.5, 1.0);
	obj_material(&w, "bspline", 0.2, 0.9, 0.3);
	obj_material(&w, "nurbs", 1.0, 0.3, 0.5);
	obj_material(&w, "cp", 0.6, 0.6, 0.6);

	/* Bezier: wavy patch, offset to x=[-3, 0] */
	{
		qaws_vec3 cp[16];
		qaws_surface_bezier_desc desc;
		for (i = 0; i < 4; i++)
			for (j = 0; j < 4; j++)
			{
				unsigned int idx = i * 4 + j;
				cp[idx].x = (qaws_scalar)i - 5;
				cp[idx].y = (qaws_scalar)j;
				cp[idx].z = (qaws_scalar)(sin((double)i * 1.2) * sin((double)j * 1.2));
			}
		desc.u_degree = 3;
		desc.v_degree = 3;
		desc.control_points = cp;
		desc.u_point_count = 4;
		desc.v_point_count = 4;
		if (qaws_surface_create_bezier(&desc, &bezier_surf) == QAWS_STATUS_OK)
		{
			obj_group(&w, "bezier_patch");
			obj_use_material(&w, "bezier");
			obj_surface_mesh(&w, bezier_surf, 24, 24);
			obj_group(&w, "bezier_cp");
			obj_use_material(&w, "cp");
			obj_surface_cp_grid(&w, cp, 4, 4, (qaws_scalar)0.05);
		}
	}

	/* B-spline: terrain, centered at x=[0, 5] */
	{
		qaws_vec3 cp[36];
		qaws_surface_bspline_desc desc;
		for (i = 0; i < 6; i++)
			for (j = 0; j < 6; j++)
			{
				unsigned int idx = i * 6 + j;
				cp[idx].x = (qaws_scalar)i;
				cp[idx].y = (qaws_scalar)j;
				cp[idx].z = (qaws_scalar)(0.6 * cos((double)i * 0.9) * sin((double)j * 0.7));
			}
		memset(&desc, 0, sizeof(desc));
		desc.u_degree = 3;
		desc.v_degree = 3;
		desc.control_points = cp;
		desc.u_point_count = 6;
		desc.v_point_count = 6;
		if (qaws_surface_create_bspline(&desc, &bspline_surf) == QAWS_STATUS_OK)
		{
			obj_group(&w, "bspline_surface");
			obj_use_material(&w, "bspline");
			obj_surface_mesh(&w, bspline_surf, 32, 32);
			obj_group(&w, "bspline_cp");
			obj_use_material(&w, "cp");
			obj_surface_cp_grid(&w, cp, 6, 6, (qaws_scalar)0.05);
		}
	}

	/* NURBS: weighted dome, offset to x=[7, 10] */
	{
		qaws_vec3 cp[9];
		qaws_scalar wt[9];
		qaws_surface_nurbs_desc desc;
		for (i = 0; i < 3; i++)
			for (j = 0; j < 3; j++)
			{
				unsigned int idx = i * 3 + j;
				cp[idx].x = (qaws_scalar)i + 7;
				cp[idx].y = (qaws_scalar)j;
				cp[idx].z = (i == 1 && j == 1) ? (qaws_scalar)2.0 : 0;
				wt[idx] = (i == 1 && j == 1) ? (qaws_scalar)2.0 : (qaws_scalar)1.0;
			}
		memset(&desc, 0, sizeof(desc));
		desc.u_degree = 2;
		desc.v_degree = 2;
		desc.control_points = cp;
		desc.u_point_count = 3;
		desc.v_point_count = 3;
		desc.weights = wt;
		if (qaws_surface_create_nurbs(&desc, &nurbs_surf) == QAWS_STATUS_OK)
		{
			obj_group(&w, "nurbs_surface");
			obj_use_material(&w, "nurbs");
			obj_surface_mesh(&w, nurbs_surf, 24, 24);
			obj_group(&w, "nurbs_cp");
			obj_use_material(&w, "cp");
			obj_surface_cp_grid(&w, cp, 3, 3, (qaws_scalar)0.06);
		}
	}

	obj_close(&w);
	if (bezier_surf) qaws_surface_destroy(bezier_surf);
	if (bspline_surf) qaws_surface_destroy(bspline_surf);
	if (nurbs_surf) qaws_surface_destroy(nurbs_surf);
	printf("  -> " OBJ_OUTPUT_DIR "/30_surface_gallery.obj\n");
}

static void test_obj_surface_ruled(void)
{
	obj_writer w;
	qaws_curve* ca = NULL;
	qaws_curve* cb = NULL;
	qaws_surface* surf = NULL;
	qaws_surface_ruled_desc rdesc;

	printf("test_obj_surface_ruled\n");

	/* Curve A: wavy line at z=0 */
	{
		qaws_scalar pts[] = {0,0,0, 1,2,0, 3,-2,0, 5,1,0, 7,0,0};
		qaws_bezier_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_3D;
		d.degree = 4;
		d.control_points = pts;
		d.control_point_count = 5;
		qaws_curve_create_bezier(&d, &ca);
	}

	/* Curve B: wavy line at z=3, different shape */
	{
		qaws_scalar pts[] = {0,1,3, 2,-1,4, 4,2,3, 6,-1,3.5f, 7,0,3};
		qaws_bezier_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_3D;
		d.degree = 4;
		d.control_points = pts;
		d.control_point_count = 5;
		qaws_curve_create_bezier(&d, &cb);
	}

	rdesc.curve_a = ca;
	rdesc.curve_b = cb;
	if (qaws_surface_create_ruled(&rdesc, &surf) != QAWS_STATUS_OK)
	{
		qaws_curve_destroy(ca);
		qaws_curve_destroy(cb);
		return;
	}

	if (!obj_open(&w,
		OBJ_OUTPUT_DIR "/30_surface_ruled.obj",
		OBJ_OUTPUT_DIR "/30_surface_ruled.mtl"))
	{
		qaws_surface_destroy(surf);
		qaws_curve_destroy(ca);
		qaws_curve_destroy(cb);
		return;
	}

	obj_material(&w, "ruled", 0.9, 0.7, 0.2);
	obj_material(&w, "edge_a", 1.0, 0.2, 0.2);
	obj_material(&w, "edge_b", 0.2, 0.2, 1.0);

	obj_group(&w, "ruled_surface");
	obj_use_material(&w, "ruled");
	obj_surface_mesh(&w, surf, 40, 20);

	/* Draw boundary curves */
	obj_group(&w, "curve_a");
	obj_use_material(&w, "edge_a");
	obj_curve_polyline(&w, ca, 80);

	obj_group(&w, "curve_b");
	obj_use_material(&w, "edge_b");
	obj_curve_polyline(&w, cb, 80);

	obj_close(&w);
	qaws_surface_destroy(surf);
	qaws_curve_destroy(ca);
	qaws_curve_destroy(cb);
	printf("  -> " OBJ_OUTPUT_DIR "/30_surface_ruled.obj\n");
}

static void test_obj_surface_swept(void)
{
	obj_writer w;
	qaws_curve* path = NULL;
	qaws_curve* profile = NULL;
	qaws_surface* surf = NULL;
	qaws_surface_swept_desc sdesc;

	printf("test_obj_surface_swept\n");

	/* Path: helix-like 3D curve */
	{
		unsigned int i;
		unsigned int n = 25;
		qaws_scalar pts[75]; /* 25 * 3 */
		qaws_scalar der[75];
		qaws_hermite_desc d;
		double dt_scale = 4.0 * M_PI / (double)(n - 1);
		for (i = 0; i < n; i++)
		{
			double t = 4.0 * M_PI * (double)i / (double)(n - 1);
			pts[i * 3 + 0] = (qaws_scalar)(2.0 * cos(t));
			pts[i * 3 + 1] = (qaws_scalar)(2.0 * sin(t));
			pts[i * 3 + 2] = (qaws_scalar)(t * 0.3);
			der[i * 3 + 0] = (qaws_scalar)(-2.0 * sin(t) * dt_scale);
			der[i * 3 + 1] = (qaws_scalar)(2.0 * cos(t) * dt_scale);
			der[i * 3 + 2] = (qaws_scalar)(0.3 * dt_scale);
		}
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_3D;
		d.degree = 3;
		d.points = pts;
		d.derivatives = der;
		d.point_count = n;
		d.derivative_count = n;
		qaws_curve_create_hermite(&d, &path);
	}

	/* Profile: 5-pointed star cross-section */
	{
		unsigned int i;
		unsigned int n = 11; /* 5 outer + 5 inner + close */
		qaws_scalar pts[22]; /* 11 * 2 */
		qaws_scalar der[22];
		qaws_hermite_desc d;
		for (i = 0; i < n; i++)
		{
			double a = 2.0 * M_PI * (double)i / (double)(n - 1);
			/* Alternate between outer radius (1.0) and inner radius (0.4) */
			double r_val = (i % 2 == 0) ? 1.0 : 0.4;
			double r_next = ((i + 1) % 2 == 0) ? 1.0 : 0.4;
			double a_next = 2.0 * M_PI * (double)(i + 1) / (double)(n - 1);
			pts[i * 2 + 0] = (qaws_scalar)(r_val * cos(a));
			pts[i * 2 + 1] = (qaws_scalar)(r_val * sin(a));
			/* Derivative: finite difference direction toward next point, scaled */
			der[i * 2 + 0] = (qaws_scalar)((r_next * cos(a_next) - r_val * cos(a)) * 0.5);
			der[i * 2 + 1] = (qaws_scalar)((r_next * sin(a_next) - r_val * sin(a)) * 0.5);
		}
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D;
		d.degree = 3;
		d.points = pts;
		d.derivatives = der;
		d.point_count = n;
		d.derivative_count = n;
		qaws_curve_create_hermite(&d, &profile);
	}

	memset(&sdesc, 0, sizeof(sdesc));
	sdesc.path = path;
	sdesc.profile = profile;
	sdesc.scale = (qaws_scalar)0.3;

	if (qaws_surface_create_swept(&sdesc, &surf) != QAWS_STATUS_OK)
	{
		if (path) qaws_curve_destroy(path);
		if (profile) qaws_curve_destroy(profile);
		return;
	}

	if (!obj_open(&w,
		OBJ_OUTPUT_DIR "/30_surface_swept.obj",
		OBJ_OUTPUT_DIR "/30_surface_swept.mtl"))
	{
		qaws_surface_destroy(surf);
		qaws_curve_destroy(path);
		qaws_curve_destroy(profile);
		return;
	}

	obj_material(&w, "swept", 0.3, 0.7, 0.9);
	obj_material(&w, "path_line", 1.0, 0.3, 0.1);

	obj_group(&w, "swept_surface");
	obj_use_material(&w, "swept");
	obj_surface_mesh(&w, surf, 120, 24);

	/* Draw path curve as wireframe */
	obj_group(&w, "path");
	obj_use_material(&w, "path_line");
	obj_curve_polyline(&w, path, 200);

	obj_close(&w);
	qaws_surface_destroy(surf);
	qaws_curve_destroy(path);
	qaws_curve_destroy(profile);
	printf("  -> " OBJ_OUTPUT_DIR "/30_surface_swept.obj\n");
}

int test_30_obj_surfaces_main(void)
{
	svg_ensure_output_dir();
	test_obj_surface_bezier();
	test_obj_surface_bspline();
	test_obj_surface_nurbs();
	test_obj_surface_gallery();
	test_obj_surface_ruled();
	test_obj_surface_swept();
	return 0;
}
