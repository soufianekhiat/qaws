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
extern int test_76_exact_composite_main(void);
extern int test_77_exact_curve_surface_main(void);
extern int test_78_exact_ssi_main(void);
extern int test_79_curve_batch_main(void);
extern int test_80_exact_batch_main(void);
extern int test_81_surface_batch_main(void);
extern int test_82_curve_extract_main(void);
extern int test_83_path_main(void);
extern int test_84_clip_main(void);
extern int test_85_clipper2_parity_main(void);
extern int test_86_clip64_main(void);
extern int test_87_offset_main(void);
extern int test_88_rectclip_main(void);
extern int test_89_stack_main(void);
extern int test_90_offset_curves_main(void);
extern int test_91_triangulate_main(void);
extern int test_92_stack_mesh_main(void);
extern int test_93_offset_parity_main(void);
extern int test_94_polytree_main(void);
extern int test_95_jobs_main(void);
extern int test_96_exact_contacts_main(void);
extern int test_97_clip_exact_main(void);
extern int test_98_region_sample_main(void);
extern int test_99_skeleton_main(void);
extern int test_100_mean_value_main(void);
extern int test_101_surface_geodesic_main(void);
extern int test_102_surface_region_main(void);
extern int test_103_surface_offset_main(void);
extern int test_104_surface_sample_main(void);
extern int test_105_surface_exact_main(void);
extern int test_106_surface_medial_main(void);
extern int test_107_shell_region_main(void);
extern int test_108_shell_trimmed_main(void);
extern int test_109_geodesic_enclosure_main(void);
extern int test_110_geodesic_query_main(void);
extern int test_111_shell_mesh_main(void);
extern int test_112_repulsive_main(void);
extern int test_113_repulsive_solve_main(void);
extern int test_114_repulsive_mesh_solve_main(void);
extern int test_115_repulsive_bh_main(void);
extern int test_116_repulsive_hsolve_main(void);
extern int test_117_repulsive_surface_main(void);
extern int test_118_repulsive_controls_main(void);
extern int test_119_repulsive_core_main(void);
extern int test_120_repulsive_spacetime_main(void);
extern int test_121_mesh_main(void);
extern int test_122_mesh_reference_main(void);
extern int test_123_size_field_main(void);
extern int test_124_remesh_main(void);
extern int test_125_quadrangulate_main(void);
extern int test_126_mesh_generate_main(void);
extern int test_127_tetrahedralize_main(void);
extern int test_128_hexahedralize_main(void);
extern int test_129_mesh_2d_main(void);
extern int test_130_mesh_measure_main(void);
extern int test_131_remesh_quality_main(void);
extern int test_132_remesh_cvt_main(void);
extern int test_133_cross_field_main(void);
extern int test_134_quad_flow_main(void);
extern int test_135_quad_improve_main(void);
extern int test_136_quad_global_main(void);
extern int test_137_mesh_guidance_main(void);
extern int test_138_feature_size_main(void);
extern int test_139_mesh_repair_main(void);
extern int test_140_eval_limits_main(void);
extern int test_141_eval_threads_main(void);
extern int test_142_curve_view_main(void);
extern int test_143_surface_view_main(void);
extern int test_144_eval_consistency_main(void);
extern int test_145_repulsive_control_main(void);
extern int test_146_allocator_main(void);
extern int test_147_frames_main(void);
extern int test_148_arc_length_main(void);
extern int test_149_closest_point_main(void);
extern int test_150_curve_tessellation_main(void);
extern int test_151_surface_tessellation_main(void);
extern int test_152_sweep_derivatives_main(void);
extern int test_153_parameterization_main(void);
extern int test_154_memory_main(void);
extern int test_155_precision_both_main(void);
extern int test_156_boundary_layers_main(void);
extern int test_157_curved_elements_main(void);
extern int test_158_metric_remesh_main(void);
extern int test_159_mesh_bvh_main(void);
extern int test_160_quality_report_main(void);
extern int test_161_watertight_main(void);
extern int test_162_intersect_main(void);
extern int test_163_boolean_main(void);
extern int test_164_mesh_geodesic_main(void);
extern int test_165_mesh_smooth_main(void);
extern int test_166_mesh_deform_main(void);
extern int test_167_mesh_decimate_main(void);
extern int test_168_remesh_attributes_main(void);
extern int test_169_mesh_connect_main(void);
extern int test_170_repulsive_mesh_curves_main(void);
extern int test_171_mesh_offset_main(void);
extern int test_172_mesh_subdivide_main(void);
extern int test_173_volume_improve_main(void);
extern int test_174_volume_curve_main(void);
extern int test_175_tetrahedra_border_main(void);
extern int test_176_hexahedra_polycube_main(void);
extern int test_177_quads_pairing_main(void);
extern int test_178_quad_layout_main(void);

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
	{"76_exact_composite", test_76_exact_composite_main},
	{"77_exact_curve_surface", test_77_exact_curve_surface_main},
	{"78_exact_ssi", test_78_exact_ssi_main},
	{"79_curve_batch", test_79_curve_batch_main},
	{"80_exact_batch", test_80_exact_batch_main},
	{"81_surface_batch", test_81_surface_batch_main},
	{"82_curve_extract", test_82_curve_extract_main},
	{"83_path", test_83_path_main},
	{"84_clip", test_84_clip_main},
	{"85_clipper2_parity", test_85_clipper2_parity_main},
	{"86_clip64", test_86_clip64_main},
	{"87_offset", test_87_offset_main},
	{"88_rectclip", test_88_rectclip_main},
	{"89_stack", test_89_stack_main},
	{"90_offset_curves", test_90_offset_curves_main},
	{"91_triangulate", test_91_triangulate_main},
	{"92_stack_mesh", test_92_stack_mesh_main},
	{"93_offset_parity", test_93_offset_parity_main},
	{"94_polytree", test_94_polytree_main},
	{"95_jobs", test_95_jobs_main},
	{"96_exact_contacts", test_96_exact_contacts_main},
	{"97_clip_exact", test_97_clip_exact_main},
	{"98_region_sample", test_98_region_sample_main},
	{"99_skeleton", test_99_skeleton_main},
	{"100_mean_value", test_100_mean_value_main},
	{"101_surface_geodesic", test_101_surface_geodesic_main},
	{"102_surface_region", test_102_surface_region_main},
	{"103_surface_offset", test_103_surface_offset_main},
	{"104_surface_sample", test_104_surface_sample_main},
	{"105_surface_exact", test_105_surface_exact_main},
	{"106_surface_medial", test_106_surface_medial_main},
	{"107_shell_region", test_107_shell_region_main},
	{"108_shell_trimmed", test_108_shell_trimmed_main},
	{"109_geodesic_enclosure", test_109_geodesic_enclosure_main},
	{"110_geodesic_query", test_110_geodesic_query_main},
	{"111_shell_mesh", test_111_shell_mesh_main},
	{"112_repulsive", test_112_repulsive_main},
	{"113_repulsive_solve", test_113_repulsive_solve_main},
	{"114_repulsive_mesh_solve", test_114_repulsive_mesh_solve_main},
	{"115_repulsive_bh", test_115_repulsive_bh_main},
	{"116_repulsive_hsolve", test_116_repulsive_hsolve_main},
	{"117_repulsive_surface", test_117_repulsive_surface_main},
	{"118_repulsive_controls", test_118_repulsive_controls_main},
	{"119_repulsive_core", test_119_repulsive_core_main},
	{"120_repulsive_spacetime", test_120_repulsive_spacetime_main},
	{"121_mesh", test_121_mesh_main},
	{"122_mesh_reference", test_122_mesh_reference_main},
	{"123_size_field", test_123_size_field_main},
	{"124_remesh", test_124_remesh_main},
	{"125_quadrangulate", test_125_quadrangulate_main},
	{"126_mesh_generate", test_126_mesh_generate_main},
	{"127_tetrahedralize", test_127_tetrahedralize_main},
	{"128_hexahedralize", test_128_hexahedralize_main},
	{"129_mesh_2d", test_129_mesh_2d_main},
{"130_mesh_measure", test_130_mesh_measure_main},
{"131_remesh_quality", test_131_remesh_quality_main},
{"132_remesh_cvt", test_132_remesh_cvt_main},
{"133_cross_field", test_133_cross_field_main},
{"134_quad_flow", test_134_quad_flow_main},
{"135_quad_improve", test_135_quad_improve_main},
{"136_quad_global", test_136_quad_global_main},
{"137_mesh_guidance", test_137_mesh_guidance_main},
{"138_feature_size", test_138_feature_size_main},
{"139_mesh_repair", test_139_mesh_repair_main},
{"140_eval_limits", test_140_eval_limits_main},
{"141_eval_threads", test_141_eval_threads_main},
{"142_curve_view", test_142_curve_view_main},
{"143_surface_view", test_143_surface_view_main},
{"144_eval_consistency", test_144_eval_consistency_main},
{"145_repulsive_control", test_145_repulsive_control_main},
{"146_allocator", test_146_allocator_main},
{"147_frames", test_147_frames_main},
{"148_arc_length", test_148_arc_length_main},
{"149_closest_point", test_149_closest_point_main},
{"150_curve_tessellation", test_150_curve_tessellation_main},
{"151_surface_tessellation", test_151_surface_tessellation_main},
{"152_sweep_derivatives", test_152_sweep_derivatives_main},
{"153_parameterization", test_153_parameterization_main},
{"154_memory", test_154_memory_main},
{"155_precision_both", test_155_precision_both_main},
{"156_boundary_layers", test_156_boundary_layers_main},
{"157_curved_elements", test_157_curved_elements_main},
{"158_metric_remesh", test_158_metric_remesh_main},
{"159_mesh_bvh", test_159_mesh_bvh_main},
{"160_quality_report", test_160_quality_report_main},
{"161_watertight", test_161_watertight_main},
{"162_intersect", test_162_intersect_main},
{"163_boolean", test_163_boolean_main},
	{"164_mesh_geodesic", test_164_mesh_geodesic_main},
	{"165_mesh_smooth", test_165_mesh_smooth_main},
	{"166_mesh_deform", test_166_mesh_deform_main},
	{"167_mesh_decimate", test_167_mesh_decimate_main},
	{"168_remesh_attributes", test_168_remesh_attributes_main},
	{"169_mesh_connect", test_169_mesh_connect_main},
	{"170_repulsive_mesh_curves", test_170_repulsive_mesh_curves_main},
	{"171_mesh_offset", test_171_mesh_offset_main},
	{"172_mesh_subdivide", test_172_mesh_subdivide_main},
	{"173_volume_improve", test_173_volume_improve_main},
	{"174_volume_curve", test_174_volume_curve_main},
	{"175_tetrahedra_border", test_175_tetrahedra_border_main},
	{"176_hexahedra_polycube", test_176_hexahedra_polycube_main},
	{"177_quads_pairing", test_177_quads_pairing_main},
	{"178_quad_layout", test_178_quad_layout_main},
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
