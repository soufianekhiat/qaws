/*
 * Qaws - Dependency-free C11 curve evaluation and traversal library
 * Copyright (c) 2026 Soufiane KHIAT
 * SPDX-License-Identifier: MIT
 *
 * Main test runner - runs all unit test suites consecutively
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Forward declarations of test main functions */
extern int test_00_status_main(void);
extern int test_01_bezier_main(void);
extern int test_02_hermite_main(void);
extern int test_03_catmull_rom_main(void);
extern int test_04_bspline_main(void);
extern int test_05_nurbs_main(void);
extern int test_06_trajectory_main(void);
extern int test_07_yuksel_main(void);
extern int test_08_rational_bezier_main(void);
extern int test_09_arc_main(void);
extern int test_10_polynomial_main(void);
extern int test_11_clothoid_main(void);
extern int test_12_subdivision_main(void);
extern int test_13_composite_main(void);
extern int test_14_evaluation_main(void);
extern int test_15_sampling_main(void);
extern int test_16_traversal_main(void);
extern int test_17_inspection_main(void);
extern int test_18_operations_main(void);
extern int test_19_conversion_main(void);
extern int test_20_intersection_main(void);
extern int test_21_export_import_main(void);
extern int test_22_surfaces_main(void);
extern int test_23_performance_main(void);
extern int test_24_error_handling_main(void);
extern int test_25_svg_curves_main(void);
extern int test_26_svg_operations_main(void);
extern int test_27_svg_analysis_main(void);
extern int test_28_svg_families_main(void);
extern int test_29_obj_curves_main(void);
extern int test_30_obj_surfaces_main(void);
extern int test_31_surface_analysis_main(void);
extern int test_32_surface_modeling_main(void);
extern int test_33_advanced_inspection_main(void);
extern int test_34_surface_patches_main(void);
extern int test_35_surface_pipe_main(void);
extern int test_36_surface_loft_main(void);
extern int test_37_surface_gordon_main(void);
extern int test_38_surface_offset_main(void);
extern int test_39_surface_trim_main(void);
extern int test_40_surface_intersect_main(void);
extern int test_41_curve_projection_main(void);
extern int test_42_curve_operations_main(void);
extern int test_43_boolean_2d_main(void);
extern int test_44_geodesic_main(void);
extern int test_45_tspline_main(void);
extern int test_46_subdiv_main(void);
extern int test_47_brep_main(void);
extern int test_48_fillet_main(void);
extern int test_49_diff_model_main(void);
extern int test_50_diff_curves_main(void);
extern int test_51_diff_surfaces_main(void);
extern int test_52_diff_geometry_main(void);
extern int test_53_diff_derived_main(void);
extern int test_54_diff_implicit_main(void);
extern int test_55_diff_maps_main(void);
extern int test_56_diff_functionals_main(void);
extern int test_57_diff_core_main(void);
extern int test_58_diff_sampling_main(void);
extern int test_59_diff_surface_sampling_main(void);
extern int test_60_diff_sampling_core_main(void);
extern int test_61_diff_surface_functional_core_main(void);
extern int test_62_diff_surface_sampling_core_main(void);
extern int test_63_exact_int_main(void);
extern int test_64_exact_predicates_main(void);
extern int test_65_exact_bezier_main(void);
extern int test_66_exact_winding_main(void);
extern int test_67_exact_spline_main(void);
extern int test_68_exact_families_main(void);
extern int test_69_exact_surface_main(void);
extern int test_70_exact_hits_main(void);
extern int test_71_exact_curve_hits_main(void);
extern int test_72_exact_curve_hits_3d_main(void);
extern int test_73_exact_self_hits_main(void);
extern int test_74_exact_surface_hits_main(void);
extern int test_75_exact_boolean_main(void);

/* Test registry */
typedef struct {
	char const *name;
	int (*test_fn)(void);
} test_suite;

static test_suite const g_test_suites[] = {
	{"00_status", test_00_status_main},
	{"01_bezier", test_01_bezier_main},
	{"02_hermite", test_02_hermite_main},
	{"03_catmull_rom", test_03_catmull_rom_main},
	{"04_bspline", test_04_bspline_main},
	{"05_nurbs", test_05_nurbs_main},
	{"06_trajectory", test_06_trajectory_main},
	{"07_yuksel", test_07_yuksel_main},
	{"08_rational_bezier", test_08_rational_bezier_main},
	{"09_arc", test_09_arc_main},
	{"10_polynomial", test_10_polynomial_main},
	{"11_clothoid", test_11_clothoid_main},
	{"12_subdivision", test_12_subdivision_main},
	{"13_composite", test_13_composite_main},
	{"14_evaluation", test_14_evaluation_main},
	{"15_sampling", test_15_sampling_main},
	{"16_traversal", test_16_traversal_main},
	{"17_inspection", test_17_inspection_main},
	{"18_operations", test_18_operations_main},
	{"19_conversion", test_19_conversion_main},
	{"20_intersection", test_20_intersection_main},
	{"21_export_import", test_21_export_import_main},
	{"22_surfaces", test_22_surfaces_main},
	{"23_performance", test_23_performance_main},
	{"24_error_handling", test_24_error_handling_main},
	{"25_svg_curves", test_25_svg_curves_main},
	{"26_svg_operations", test_26_svg_operations_main},
	{"27_svg_analysis", test_27_svg_analysis_main},
	{"28_svg_families", test_28_svg_families_main},
	{"29_obj_curves", test_29_obj_curves_main},
	{"30_obj_surfaces", test_30_obj_surfaces_main},
	{"31_surface_analysis", test_31_surface_analysis_main},
	{"32_surface_modeling", test_32_surface_modeling_main},
	{"33_advanced_inspection", test_33_advanced_inspection_main},
	{"34_surface_patches", test_34_surface_patches_main},
	{"35_surface_pipe", test_35_surface_pipe_main},
	{"36_surface_loft", test_36_surface_loft_main},
	{"37_surface_gordon", test_37_surface_gordon_main},
	{"38_surface_offset", test_38_surface_offset_main},
	{"39_surface_trim", test_39_surface_trim_main},
	{"40_surface_intersect", test_40_surface_intersect_main},
	{"41_curve_projection", test_41_curve_projection_main},
	{"42_curve_operations", test_42_curve_operations_main},
	{"43_boolean_2d", test_43_boolean_2d_main},
	{"44_geodesic", test_44_geodesic_main},
	{"45_tspline", test_45_tspline_main},
	{"46_subdiv", test_46_subdiv_main},
	{"47_brep", test_47_brep_main},
	{"48_fillet", test_48_fillet_main},
	{"49_diff_model", test_49_diff_model_main},
	{"50_diff_curves", test_50_diff_curves_main},
	{"51_diff_surfaces", test_51_diff_surfaces_main},
	{"52_diff_geometry", test_52_diff_geometry_main},
	{"53_diff_derived", test_53_diff_derived_main},
	{"54_diff_implicit", test_54_diff_implicit_main},
	{"55_diff_maps", test_55_diff_maps_main},
	{"56_diff_functionals", test_56_diff_functionals_main},
	{"57_diff_core", test_57_diff_core_main},
	{"58_diff_sampling", test_58_diff_sampling_main},
	{"59_diff_surface_sampling", test_59_diff_surface_sampling_main},
	{"60_diff_sampling_core", test_60_diff_sampling_core_main},
	{"61_diff_surface_functional_core", test_61_diff_surface_functional_core_main},
	{"62_diff_surface_sampling_core", test_62_diff_surface_sampling_core_main},
	{"63_exact_int", test_63_exact_int_main},
	{"64_exact_predicates", test_64_exact_predicates_main},
	{"65_exact_bezier", test_65_exact_bezier_main},
	{"66_exact_winding", test_66_exact_winding_main},
	{"67_exact_spline", test_67_exact_spline_main},
	{"68_exact_families", test_68_exact_families_main},
	{"69_exact_surface", test_69_exact_surface_main},
	{"70_exact_hits", test_70_exact_hits_main},
	{"71_exact_curve_hits", test_71_exact_curve_hits_main},
	{"72_exact_curve_hits_3d", test_72_exact_curve_hits_3d_main},
	{"73_exact_self_hits", test_73_exact_self_hits_main},
	{"74_exact_surface_hits", test_74_exact_surface_hits_main},
	{"75_exact_boolean", test_75_exact_boolean_main},
};

int main(void) {
	/* unbuffered, so a crash shows where it happened; QAWS_TEST=substring
	   runs only the matching suites */
	char const* only = getenv("QAWS_TEST");
	setvbuf(stdout, NULL, _IONBF, 0);
	printf("========================================\n");
	printf("Qaws Unit Test Runner\n");
	printf("========================================\n\n");

	int total_failed = 0;
	int total_passed = 0;
	size_t const num_suites = sizeof(g_test_suites) / sizeof(g_test_suites[0]);

	for (size_t i = 0; i < num_suites; i++) {
		if (only && !strstr(g_test_suites[i].name, only))
			continue;
		printf("\n");
		printf("========================================\n");
		printf("Running test suite: %s\n", g_test_suites[i].name);
		printf("========================================\n");

		int result = g_test_suites[i].test_fn();

		if (result == 0) {
			total_passed++;
			printf("[PASS] Test suite '%s' passed\n", g_test_suites[i].name);
		} else {
			total_failed++;
			printf("[FAIL] Test suite '%s' failed with code %d\n",
				g_test_suites[i].name, result);
		}
	}

	printf("\n");
	printf("========================================\n");
	printf("Overall Results\n");
	printf("========================================\n");
	printf("Test suites passed: %d\n", total_passed);
	printf("Test suites failed: %d\n", total_failed);
	printf("Total test suites:  %zu\n", num_suites);
	printf("========================================\n");

	if (total_failed > 0) {
		printf("\nSome tests FAILED!\n");
		return 1;
	} else {
		printf("\nAll tests PASSED!\n");
		return 0;
	}
}
