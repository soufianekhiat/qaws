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
  21_export_import.c
  22_surfaces.c
  23_performance.c
  24_error_handling.c
  25_svg_curves.c
  26_svg_operations.c
  28_svg_families.c
  30_obj_surfaces.c
  31_surface_analysis.c
  32_surface_modeling.c
  34_surface_patches.c
  35_surface_pipe.c
  36_surface_loft.c
  37_surface_gordon.c
  38_surface_offset.c
  39_surface_trim.c
  42_curve_operations.c
  45_tspline.c
  46_subdiv.c
  47_brep.c
  49_diff_model.c
  50_diff_curves.c
  51_diff_surfaces.c
  52_diff_geometry.c
  53_diff_derived.c
  55_diff_maps.c
  57_diff_core.c
  59_diff_surface_sampling.c
  60_diff_sampling_core.c
  62_diff_surface_sampling_core.c
  63_exact_int.c
  64_exact_predicates.c
  65_exact_bezier.c
  66_exact_winding.c
  67_exact_spline.c
  68_exact_families.c
  69_exact_surface.c
  76_exact_composite.c
  82_curve_extract.c
  83_path.c
  84_interpolating_offset.c
140_eval_limits.c
141_eval_threads.c
142_curve_view.c
143_surface_view.c
144_eval_consistency.c
147_frames.c
148_arc_length.c
149_closest_point.c
150_curve_tessellation.c
151_surface_tessellation.c
152_sweep_derivatives.c
153_parameterization.c
155_precision_both.c
156_surface_enclose.c
157_curve_enclose.c
179_motion.c
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
  "${_QAWS_SRC}/internal/qaws_internal_memory.h" "${_QAWS_SRC}/internal/qaws_internal_cpu.h")
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
