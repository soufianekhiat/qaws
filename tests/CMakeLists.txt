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
set_target_properties(qaws_tests PROPERTIES C_STANDARD 11 C_STANDARD_REQUIRED YES C_EXTENSIONS NO)

if(MSVC)
  target_compile_definitions(qaws_tests PRIVATE _CRT_SECURE_NO_WARNINGS NOMINMAX)
endif()

set(QAWS_WORKING_DIR "${PROJECT_SOURCE_DIR}/working_dir")
file(MAKE_DIRECTORY "${QAWS_WORKING_DIR}")

set_target_properties(qaws_tests PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${QAWS_WORKING_DIR}"
  RUNTIME_OUTPUT_DIRECTORY_DEBUG "${QAWS_WORKING_DIR}"
  RUNTIME_OUTPUT_DIRECTORY_RELEASE "${QAWS_WORKING_DIR}"
  VS_DEBUGGER_WORKING_DIRECTORY "${QAWS_WORKING_DIR}"
)

add_test(NAME qaws_tests COMMAND qaws_tests)
set_tests_properties(qaws_tests PROPERTIES WORKING_DIRECTORY "${QAWS_WORKING_DIR}")
