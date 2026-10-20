add_executable(
  qaws_tests
  test_runner.c
  00_status.c
  01_bezier.c
  02_hermite.c
  03_catmull_rom.c
  04_bspline.c
  05_nurbs.c
  06_trajectory.c
  07_yuksel.c
  08_rational_bezier.c
  09_arc.c
  10_polynomial.c
  11_clothoid.c
  12_subdivision.c
  13_composite.c
  14_evaluation.c
  15_sampling.c
  16_traversal.c
  17_inspection.c
  18_operations.c
  19_conversion.c
  20_intersection.c
  21_export_import.c
  22_surfaces.c
  23_performance.c
  24_error_handling.c
  25_svg_curves.c
  26_svg_operations.c
  27_svg_analysis.c
  28_svg_families.c
  29_obj_curves.c
  30_obj_surfaces.c
  31_surface_analysis.c
  32_surface_modeling.c
  33_advanced_inspection.c
  34_surface_patches.c
  35_surface_pipe.c
  36_surface_loft.c
  37_surface_gordon.c
  38_surface_offset.c
  39_surface_trim.c
  40_surface_intersect.c
  41_curve_projection.c
  42_curve_operations.c
  43_boolean_2d.c
  44_geodesic.c
  45_tspline.c
  46_subdiv.c
  47_brep.c
  48_fillet.c
  49_diff_model.c
  50_diff_curves.c
  51_diff_surfaces.c
  52_diff_geometry.c
  53_diff_derived.c
  54_diff_implicit.c
  55_diff_maps.c
  56_diff_functionals.c
  57_diff_core.c
  58_diff_sampling.c
  59_diff_surface_sampling.c
  60_diff_sampling_core.c
  61_diff_surface_functional_core.c
  62_diff_surface_sampling_core.c
  63_exact_int.c
  64_exact_predicates.c
  65_exact_bezier.c
  66_exact_winding.c
  67_exact_spline.c
  68_exact_families.c
  69_exact_surface.c
  70_exact_hits.c
  71_exact_curve_hits.c
  72_exact_curve_hits_3d.c
  73_exact_self_hits.c
  74_exact_surface_hits.c
  75_exact_boolean.c
  76_exact_composite.c
  77_exact_curve_surface.c
  78_exact_ssi.c
  79_curve_batch.c
  80_exact_batch.c
  81_surface_batch.c
  82_curve_extract.c
  83_path.c
  84_clip.c
  85_clipper2_parity.c
  86_clip64.c
  87_offset.c
  88_rectclip.c
  89_stack.c
  90_offset_curves.c
  91_triangulate.c
  92_stack_mesh.c
  93_offset_parity.c
  94_polytree.c
  95_jobs.c
  96_exact_contacts.c
  97_clip_exact.c
  98_region_sample.c
  99_skeleton.c
  100_mean_value.c
  101_surface_geodesic.c
  102_surface_region.c
  103_surface_offset.c
  104_surface_sample.c
  105_surface_exact.c
  106_surface_medial.c
  107_shell_region.c
  108_shell_trimmed.c
  109_geodesic_enclosure.c
  110_geodesic_query.c
  111_shell_mesh.c
  112_repulsive.c
  113_repulsive_solve.c
  114_repulsive_mesh_solve.c
  115_repulsive_bh.c
  116_repulsive_hsolve.c
  117_repulsive_surface.c
  118_repulsive_controls.c
  119_repulsive_core.c
  120_repulsive_spacetime.c
  121_mesh.c
  122_mesh_reference.c
  123_size_field.c
  124_remesh.c
  125_quadrangulate.c
  126_mesh_generate.c
  127_tetrahedralize.c
  128_hexahedralize.c
  129_mesh_2d.c
130_mesh_measure.c
131_remesh_quality.c
132_remesh_cvt.c
133_cross_field.c
134_quad_flow.c
135_quad_improve.c
136_quad_global.c
137_mesh_guidance.c
138_feature_size.c
139_mesh_repair.c
140_eval_limits.c
141_eval_threads.c
142_curve_view.c
143_surface_view.c
144_eval_consistency.c
145_repulsive_control.c
146_allocator.c
147_frames.c
148_arc_length.c
149_closest_point.c
150_curve_tessellation.c
151_surface_tessellation.c
152_sweep_derivatives.c
153_parameterization.c
154_memory.c
155_precision_both.c
156_boundary_layers.c
)

# The amalgamated library exports exactly the same symbols, so the tests can
# link it instead and validate the generated qaws.c against the full suite.
# They keep including the per-file headers for declarations either way.
if(QAWS_TEST_AMALGAMATION)
  target_link_libraries(qaws_tests PRIVATE Qaws::amalgam)
else()
  target_link_libraries(qaws_tests PRIVATE Qaws::qaws)
endif()

target_include_directories(qaws_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/src/qaws")

# test 81 runs the batches on a few threads (an example executor)
if(NOT WIN32)
  find_package(Threads REQUIRED)
  target_link_libraries(qaws_tests PRIVATE Threads::Threads)
endif()
set_target_properties(qaws_tests PROPERTIES C_STANDARD 11 C_STANDARD_REQUIRED YES C_EXTENSIONS NO)

if(MSVC)
  target_compile_definitions(qaws_tests PRIVATE _CRT_SECURE_NO_WARNINGS NOMINMAX)
endif()
target_compile_definitions(qaws_tests PRIVATE QAWS_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")

set(QAWS_WORKING_DIR "${PROJECT_SOURCE_DIR}/working_dir")
file(MAKE_DIRECTORY "${QAWS_WORKING_DIR}")

set_target_properties(qaws_tests PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${QAWS_WORKING_DIR}"
  RUNTIME_OUTPUT_DIRECTORY_DEBUG "${QAWS_WORKING_DIR}"
  RUNTIME_OUTPUT_DIRECTORY_RELEASE "${QAWS_WORKING_DIR}"
  VS_DEBUGGER_WORKING_DIRECTORY "${QAWS_WORKING_DIR}"
)

# one CTest test per suite, named as its file: `ctest -R "^13"` picks some, `ctest -j 4` runs a few at once (a laptop CPU gains little past its performance cores)
# (QAWS_TEST_HEAVY=1 in the environment adds the heavy checks); qaws_tests run by hand still runs them all in turn
get_target_property(_qaws_test_sources qaws_tests SOURCES)
foreach(_qaws_test_source IN LISTS _qaws_test_sources)
  if(_qaws_test_source MATCHES "^([0-9]+_[A-Za-z0-9_]+)\\.c$")
    add_test(NAME ${CMAKE_MATCH_1} COMMAND qaws_tests)
    set_tests_properties(${CMAKE_MATCH_1} PROPERTIES WORKING_DIRECTORY "${QAWS_WORKING_DIR}" ENVIRONMENT "QAWS_TEST==${CMAKE_MATCH_1}")
  endif()
endforeach()

# The API checks (docs/parameter_convention_plan.md): every declaration in
# the parameter convention, less the shrinking allow-list; qaws_names.h the
# map of the headers' precision-bearing names. The once-only headers
# (precision-free: the status, the executor, the allocator and the memory
# module) are declared once and their names left unsuffixed. The mesh
# module's internal functions have their own prefixes (hm_, cf_, ...): -any.
add_executable(qaws_api_tool "${PROJECT_SOURCE_DIR}/tools/qaws_api_tool.c")
if(MSVC)
  target_compile_definitions(qaws_api_tool PRIVATE _CRT_SECURE_NO_WARNINGS)
endif()
set(_QAWS_SRC "${PROJECT_SOURCE_DIR}/src/qaws")
file(GLOB _QAWS_PUBLIC_HEADERS "${_QAWS_SRC}/*.h")
list(REMOVE_ITEM _QAWS_PUBLIC_HEADERS "${_QAWS_SRC}/qaws_names.h")
file(GLOB _QAWS_INTERNAL_HEADERS "${_QAWS_SRC}/internal/*.h")
file(GLOB _QAWS_MESH_HEADERS "${_QAWS_SRC}/mesh/*.h")
file(GLOB _QAWS_CORE_HEADERS "${_QAWS_SRC}/core/*.h")
file(GLOB _QAWS_OTHER_HEADERS "${_QAWS_SRC}/exact/*.h" "${_QAWS_SRC}/simd/*.h")
set(_QAWS_ONCE_HEADERS "${_QAWS_SRC}/qaws_status.h" "${_QAWS_SRC}/qaws_batch_executor.h" "${_QAWS_SRC}/qaws_allocator.h"
  "${_QAWS_SRC}/internal/qaws_internal_memory.h")
set(_QAWS_NAMED_HEADERS ${_QAWS_PUBLIC_HEADERS})
list(REMOVE_ITEM _QAWS_NAMED_HEADERS ${_QAWS_ONCE_HEADERS})
set(_QAWS_NAMED_INTERNAL_HEADERS ${_QAWS_INTERNAL_HEADERS})
list(REMOVE_ITEM _QAWS_NAMED_INTERNAL_HEADERS ${_QAWS_ONCE_HEADERS})
add_test(NAME qaws_api_convention
  COMMAND qaws_api_tool survey -allow "${PROJECT_SOURCE_DIR}/tools/qaws_api_allow.txt"
    ${_QAWS_PUBLIC_HEADERS} ${_QAWS_INTERNAL_HEADERS} ${_QAWS_MESH_HEADERS} ${_QAWS_CORE_HEADERS})
add_test(NAME qaws_api_names
  COMMAND qaws_api_tool names -check "${_QAWS_SRC}/qaws_names.h"
    -once ${_QAWS_ONCE_HEADERS}
    -- ${_QAWS_NAMED_HEADERS} ${_QAWS_NAMED_INTERNAL_HEADERS} ${_QAWS_OTHER_HEADERS} ${_QAWS_CORE_HEADERS}
    -any ${_QAWS_MESH_HEADERS})
