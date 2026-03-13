/*
 * Qaws - Dependency-free C11 curve evaluation and traversal library
 * Copyright (c) 2026 Soufiane KHIAT
 * SPDX-License-Identifier: MIT
 *
 * Main test runner - runs all unit test suites consecutively
 */

#include <stdio.h>
#include <stdlib.h>

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
};

int main(void) {
	printf("========================================\n");
	printf("Qaws Unit Test Runner\n");
	printf("========================================\n\n");

	int total_failed = 0;
	int total_passed = 0;
	size_t const num_suites = sizeof(g_test_suites) / sizeof(g_test_suites[0]);

	for (size_t i = 0; i < num_suites; i++) {
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
