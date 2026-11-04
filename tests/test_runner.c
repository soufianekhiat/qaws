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
#include <time.h>
#include "qaws_alloc.h"

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
extern int test_21_export_import_main(void);
extern int test_22_surfaces_main(void);
extern int test_23_performance_main(void);
extern int test_24_error_handling_main(void);
extern int test_25_svg_curves_main(void);
extern int test_26_svg_operations_main(void);
extern int test_28_svg_families_main(void);
extern int test_30_obj_surfaces_main(void);
extern int test_31_surface_analysis_main(void);
extern int test_32_surface_modeling_main(void);
extern int test_34_surface_patches_main(void);
extern int test_35_surface_pipe_main(void);
extern int test_36_surface_loft_main(void);
extern int test_37_surface_gordon_main(void);
extern int test_38_surface_offset_main(void);
extern int test_39_surface_trim_main(void);
extern int test_42_curve_operations_main(void);
extern int test_45_tspline_main(void);
extern int test_46_subdiv_main(void);
extern int test_47_brep_main(void);
extern int test_49_diff_model_main(void);
extern int test_50_diff_curves_main(void);
extern int test_51_diff_surfaces_main(void);
extern int test_52_diff_geometry_main(void);
extern int test_53_diff_derived_main(void);
extern int test_55_diff_maps_main(void);
extern int test_57_diff_core_main(void);
extern int test_59_diff_surface_sampling_main(void);
extern int test_60_diff_sampling_core_main(void);
extern int test_62_diff_surface_sampling_core_main(void);
extern int test_63_exact_int_main(void);
extern int test_64_exact_predicates_main(void);
extern int test_65_exact_bezier_main(void);
extern int test_66_exact_winding_main(void);
extern int test_67_exact_spline_main(void);
extern int test_68_exact_families_main(void);
extern int test_69_exact_surface_main(void);
extern int test_76_exact_composite_main(void);
extern int test_82_curve_extract_main(void);
extern int test_83_path_main(void);
extern int test_140_eval_limits_main(void);
extern int test_141_eval_threads_main(void);
extern int test_142_curve_view_main(void);
extern int test_143_surface_view_main(void);
extern int test_144_eval_consistency_main(void);
extern int test_147_frames_main(void);
extern int test_148_arc_length_main(void);
extern int test_149_closest_point_main(void);
extern int test_150_curve_tessellation_main(void);
extern int test_151_surface_tessellation_main(void);
extern int test_152_sweep_derivatives_main(void);
extern int test_153_parameterization_main(void);
extern int test_155_precision_both_main(void);
extern int test_156_surface_enclose_main(void);
extern int test_179_motion_main(void);

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
	{"21_export_import", test_21_export_import_main},
	{"22_surfaces", test_22_surfaces_main},
	{"23_performance", test_23_performance_main},
	{"24_error_handling", test_24_error_handling_main},
	{"25_svg_curves", test_25_svg_curves_main},
	{"26_svg_operations", test_26_svg_operations_main},
	{"28_svg_families", test_28_svg_families_main},
	{"30_obj_surfaces", test_30_obj_surfaces_main},
	{"31_surface_analysis", test_31_surface_analysis_main},
	{"32_surface_modeling", test_32_surface_modeling_main},
	{"34_surface_patches", test_34_surface_patches_main},
	{"35_surface_pipe", test_35_surface_pipe_main},
	{"36_surface_loft", test_36_surface_loft_main},
	{"37_surface_gordon", test_37_surface_gordon_main},
	{"38_surface_offset", test_38_surface_offset_main},
	{"39_surface_trim", test_39_surface_trim_main},
	{"42_curve_operations", test_42_curve_operations_main},
	{"45_tspline", test_45_tspline_main},
	{"46_subdiv", test_46_subdiv_main},
	{"47_brep", test_47_brep_main},
	{"49_diff_model", test_49_diff_model_main},
	{"50_diff_curves", test_50_diff_curves_main},
	{"51_diff_surfaces", test_51_diff_surfaces_main},
	{"52_diff_geometry", test_52_diff_geometry_main},
	{"53_diff_derived", test_53_diff_derived_main},
	{"55_diff_maps", test_55_diff_maps_main},
	{"57_diff_core", test_57_diff_core_main},
	{"59_diff_surface_sampling", test_59_diff_surface_sampling_main},
	{"60_diff_sampling_core", test_60_diff_sampling_core_main},
	{"62_diff_surface_sampling_core", test_62_diff_surface_sampling_core_main},
	{"63_exact_int", test_63_exact_int_main},
	{"64_exact_predicates", test_64_exact_predicates_main},
	{"65_exact_bezier", test_65_exact_bezier_main},
	{"66_exact_winding", test_66_exact_winding_main},
	{"67_exact_spline", test_67_exact_spline_main},
	{"68_exact_families", test_68_exact_families_main},
	{"69_exact_surface", test_69_exact_surface_main},
	{"76_exact_composite", test_76_exact_composite_main},
	{"82_curve_extract", test_82_curve_extract_main},
	{"83_path", test_83_path_main},
{"140_eval_limits", test_140_eval_limits_main},
{"141_eval_threads", test_141_eval_threads_main},
{"142_curve_view", test_142_curve_view_main},
{"143_surface_view", test_143_surface_view_main},
{"144_eval_consistency", test_144_eval_consistency_main},
{"147_frames", test_147_frames_main},
{"148_arc_length", test_148_arc_length_main},
{"149_closest_point", test_149_closest_point_main},
{"150_curve_tessellation", test_150_curve_tessellation_main},
{"151_surface_tessellation", test_151_surface_tessellation_main},
{"152_sweep_derivatives", test_152_sweep_derivatives_main},
{"153_parameterization", test_153_parameterization_main},
{"155_precision_both", test_155_precision_both_main},
{"156_surface_enclose", test_156_surface_enclose_main},
	{"179_motion", test_179_motion_main},
};

int main(void) {
	/* unbuffered, so a crash shows where it happened; QAWS_TEST=substring
	   runs only the matching suites, QAWS_TEST==name the one so named (as
	   CTest runs each suite) */
	/* a copy: a suite may change the variable (154 starts child runners), which frees the string getenv returned */
	static char only_copy[256];
	char const* only = getenv("QAWS_TEST");
	if (only)
	{
		strncpy(only_copy, only, sizeof only_copy - 1);
		only = only_copy;
	}
	setvbuf(stdout, NULL, _IONBF, 0);
	/* QAWS_TEST_SCRATCH=1: the main thread's scratch from an arena that must hold nothing between suites (a temporary
	   kept by a result that outlives its suite shows as a block held) */
	static qaws_arena scratch;
	unsigned char* scratch_memory = NULL;
	if (getenv("QAWS_TEST_SCRATCH") && atoi(getenv("QAWS_TEST_SCRATCH")))
	{
		scratch_memory = (unsigned char*)malloc((size_t)1 << 30);
		if (scratch_memory)
		{
			qaws_arena_init(&scratch, scratch_memory, (size_t)1 << 30);
			qaws_exchange_thread_scratch_allocator(NULL, &scratch.allocator);
		}
	}
	printf("========================================\n");
	printf("Qaws Unit Test Runner\n");
	printf("========================================\n\n");

	int total_failed = 0;
	int total_passed = 0;
	size_t const num_suites = sizeof(g_test_suites) / sizeof(g_test_suites[0]);

	for (size_t i = 0; i < num_suites; i++) {
		if (only && (only[0] == '=' ? strcmp(g_test_suites[i].name, only + 1) != 0 : !strstr(g_test_suites[i].name, only)))
			continue;
		printf("\n");
		printf("========================================\n");
		printf("Running test suite: %s\n", g_test_suites[i].name);
		printf("========================================\n");

		struct timespec t0, t1;
		timespec_get(&t0, TIME_UTC);
		int result = g_test_suites[i].test_fn();
		timespec_get(&t1, TIME_UTC);
		double const seconds = (double)(t1.tv_sec - t0.tv_sec) + 1e-9 * (double)(t1.tv_nsec - t0.tv_nsec);
		if (scratch_memory)
		{
			printf("  scratch arena: peak %lu bytes, %lu refused, %lu invalid frees, %lu held\n", (unsigned long)scratch.peak,
				(unsigned long)scratch.failures, (unsigned long)scratch.invalid_frees, (unsigned long)scratch.live);
			if (scratch.live || scratch.invalid_frees || qaws_arena_reset(&scratch) != QAWS_STATUS_OK)
				result += 1;
			qaws_arena_init(&scratch, scratch_memory, (size_t)1 << 30);
		}

		if (result == 0) {
			total_passed++;
			printf("[PASS] Test suite '%s' passed (%.2f s)\n", g_test_suites[i].name, seconds);
		} else {
			total_failed++;
			printf("[FAIL] Test suite '%s' failed with code %d (%.2f s)\n",
				g_test_suites[i].name, result, seconds);
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

	if (only && total_passed + total_failed == 0) {
		printf("\nNo suite matches QAWS_TEST=%s\n", only);
		return 1;
	}
	if (total_failed > 0) {
		printf("\nSome tests FAILED!\n");
		return 1;
	} else {
		printf("\nAll tests PASSED!\n");
		return 0;
	}
}
