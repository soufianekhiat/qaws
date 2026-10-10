cmake_minimum_required(VERSION 3.20)

project(
  Qaws
  VERSION 2.0.0
  DESCRIPTION "Dependency-free C11 curve evaluation and traversal library"
  LANGUAGES C
)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(CTest) # provides BUILD_TESTING
include(CheckCCompilerFlag)

option(QAWS_SCALAR_IS_FLOAT "Use float (f32) instead of double (f64)" ON)
option(QAWS_ENABLE_WARNINGS "Enable compiler warnings" ON)

# f64 builds are strict by default, f32 builds are permissive by default.
set(_QAWS_WARNINGS_AS_ERRORS_DEFAULT ON)
if(QAWS_SCALAR_IS_FLOAT)
  set(_QAWS_WARNINGS_AS_ERRORS_DEFAULT OFF)
endif()
option(
  QAWS_WARNINGS_AS_ERRORS
  "Treat warnings as errors (default: ON for f64, OFF for f32)"
  ${_QAWS_WARNINGS_AS_ERRORS_DEFAULT}
)

option(QAWS_DISABLE_STRICT_ALIASING "Compile with -fno-strict-aliasing on GCC/Clang" ON)
option(QAWS_ENABLE_SIMD "Enable SIMD batch evaluation (auto-detect ISA)" OFF)
option(QAWS_INSTALL "Enable install rules" ON)
option(QAWS_INSTALL_PKGCONFIG "Install pkg-config .pc file" ON)

# Single translation unit builds.
option(QAWS_UNITY_BUILD "Compile the library as one translation unit" OFF)
option(QAWS_AMALGAMATION "Generate and build the qaws.h + qaws.c amalgamation" OFF)
option(QAWS_TEST_AMALGAMATION "Run the test suite against the amalgamated library" OFF)

if(QAWS_TEST_AMALGAMATION AND NOT QAWS_AMALGAMATION)
  set(QAWS_AMALGAMATION ON CACHE BOOL "" FORCE)
  message(STATUS "QAWS_TEST_AMALGAMATION implies QAWS_AMALGAMATION")
endif()

set(QAWS_DEBUG_POSTFIX "d" CACHE STRING "Debug postfix for output names (e.g. 'd' -> qawsd)")
set(QAWS_DEBUG_PREFIX "" CACHE STRING "Debug prefix for output names (e.g. 'dbg_' -> dbg_qaws)")

set(_QAWS_SCALAR_IS_FLOAT_VALUE 0)
if(QAWS_SCALAR_IS_FLOAT)
  set(_QAWS_SCALAR_IS_FLOAT_VALUE 1)
endif()

set(QAWS_SOURCES
  src/qaws/qaws_status.c
  src/qaws/qaws_curve.c
  src/qaws/qaws_eval.c
  src/qaws/qaws_sampling.c
  src/qaws/qaws_traversal.c
  src/qaws/qaws_inspect.c
  src/qaws/qaws_curve_batch.c
  src/qaws/qaws_surface_batch.c
  src/qaws/qaws_bezier.c
  src/qaws/qaws_hermite.c
  src/qaws/qaws_catmull_rom.c
  src/qaws/qaws_bspline.c
  src/qaws/qaws_nurbs.c
  src/qaws/qaws_trajectory.c
  src/qaws/qaws_yuksel.c
  src/qaws/qaws_operations.c
  src/qaws/qaws_extract.c
  src/qaws/qaws_path.c
  src/qaws/qaws_clip.c
  src/qaws/qaws_clip64.c
  src/qaws/qaws_offset.c
  src/qaws/qaws_rectclip.c
  src/qaws/qaws_stack.c
  src/qaws/qaws_convert.c
  src/qaws/qaws_export.c
  src/qaws/qaws_import.c
  src/qaws/qaws_rational_bezier.c
  src/qaws/qaws_composite.c
  src/qaws/qaws_arc.c
  src/qaws/qaws_polynomial.c
  src/qaws/qaws_clothoid.c
  src/qaws/qaws_subdivision.c
  src/qaws/qaws_surface.c
  src/qaws/qaws_surface_bezier.c
  src/qaws/qaws_surface_bspline.c
  src/qaws/qaws_surface_nurbs.c
  src/qaws/qaws_surface_swept.c
  src/qaws/qaws_surface_ruled.c
  src/qaws/qaws_surface_coons.c
  src/qaws/qaws_surface_extrusion.c
  src/qaws/qaws_surface_revolution.c
  src/qaws/qaws_surface_bilinear.c
  src/qaws/qaws_surface_biquadratic.c
  src/qaws/qaws_surface_pipe.c
  src/qaws/qaws_surface_loft.c
  src/qaws/qaws_surface_gordon.c
  src/qaws/qaws_surface_offset.c
  src/qaws/qaws_surface_trim.c
  src/qaws/qaws_surface_intersect.c
  src/qaws/qaws_surface_tspline.c
  src/qaws/qaws_surface_subdiv.c
  src/qaws/qaws_surface_fillet.c
  src/qaws/qaws_brep.c
  src/qaws/qaws_boolean_2d.c
  src/qaws/qaws_inline.c
  src/qaws/qaws_prepare.c
  src/qaws/qaws_diff.c
  src/qaws/qaws_diff_surface.c
  src/qaws/qaws_diff_geometry.c
  src/qaws/qaws_diff_ops.c
  src/qaws/qaws_diff_map.c
  src/qaws/qaws_diff_functionals.c
  src/qaws/qaws_diff_sampling.c
  src/qaws/internal/qaws_internal_curve.c
  src/qaws/internal/qaws_internal_basis.c
  src/qaws/internal/qaws_internal_span.c
  src/qaws/internal/qaws_internal_arc_length.c
  src/qaws/internal/qaws_internal_solver.c
  src/qaws/internal/qaws_internal_validation.c
  src/qaws/internal/qaws_internal_surface.c
  src/qaws/internal/qaws_internal_fit.c
  src/qaws/internal/qaws_internal_mesh.c
  src/qaws/internal/qaws_internal_sparse.c
  src/qaws/internal/qaws_internal_diff.c
  src/qaws/internal/qaws_internal_broadphase.c
  src/qaws/internal/qaws_internal_flatten.c
  src/qaws/internal/qaws_internal_parallel.c
)

if(QAWS_ENABLE_SIMD)
  list(APPEND QAWS_SOURCES
    src/qaws/batch/qaws_batch_bezier.c
    src/qaws/batch/qaws_batch_bspline.c
    src/qaws/batch/qaws_batch_nurbs.c
    src/qaws/batch/qaws_batch_surface.c
  )
endif()

# ===================================================================
# Toolchain probing. Done once, then applied to every library target.
# ===================================================================

if(QAWS_ENABLE_SIMD)
  include(CheckCSourceCompiles)
  check_c_source_compiles("
    #include <immintrin.h>
    int main() { __m256 v = _mm256_setzero_ps(); (void)v; return 0; }" QAWS_HAS_AVX2)
  check_c_source_compiles("
    #include <emmintrin.h>
    int main() { __m128 v = _mm_setzero_ps(); (void)v; return 0; }" QAWS_HAS_SSE2)
  check_c_source_compiles("
    #include <arm_neon.h>
    int main() { float32x4_t v = vdupq_n_f32(0); (void)v; return 0; }" QAWS_HAS_NEON)

  set(QAWS_SIMD_DEFINE "")
  set(QAWS_SIMD_OPTIONS "")
  if(QAWS_HAS_AVX2)
    set(QAWS_SIMD_DEFINE QAWS_SIMD_AVX2=1)
    if(MSVC)
      set(QAWS_SIMD_OPTIONS /arch:AVX2)
    else()
      set(QAWS_SIMD_OPTIONS -mavx2 -mfma)
    endif()
  elseif(QAWS_HAS_SSE2)
    set(QAWS_SIMD_DEFINE QAWS_SIMD_SSE2=1)
  elseif(QAWS_HAS_NEON)
    set(QAWS_SIMD_DEFINE QAWS_SIMD_NEON=1)
  endif()
endif()

if(QAWS_DISABLE_STRICT_ALIASING AND NOT MSVC)
  check_c_compiler_flag(-fno-strict-aliasing QAWS_HAS_NO_STRICT_ALIASING)
endif()

# Everything that must hold for both the per-file library and the
# amalgamated one. Keeping it in a function stops the two drifting apart.
function(qaws_configure_target target)
  set_target_properties(
    ${target}
    PROPERTIES
      C_STANDARD 11
      C_STANDARD_REQUIRED YES
      C_EXTENSIONS NO
      POSITION_INDEPENDENT_CODE ON
  )

  target_compile_definitions(
    ${target}
    PUBLIC
      QAWS_SCALAR_IS_FLOAT=${_QAWS_SCALAR_IS_FLOAT_VALUE}
  )

  if(MSVC)
    target_compile_definitions(${target} PRIVATE _CRT_SECURE_NO_WARNINGS NOMINMAX)
  endif()

  if(NOT WIN32)
    target_link_libraries(${target} PUBLIC m)
  endif()

  if(WIN32 AND BUILD_SHARED_LIBS)
    set_target_properties(${target} PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
  endif()

  if(QAWS_ENABLE_WARNINGS)
    if(MSVC)
      target_compile_options(${target} PRIVATE /W4)
      if(QAWS_WARNINGS_AS_ERRORS)
        target_compile_options(${target} PRIVATE /WX)
      endif()
    else()
      target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
      if(QAWS_WARNINGS_AS_ERRORS)
        target_compile_options(${target} PRIVATE -Werror)
      endif()
    endif()
  endif()

  if(QAWS_ENABLE_SIMD AND QAWS_SIMD_DEFINE)
    target_compile_definitions(${target} PRIVATE ${QAWS_SIMD_DEFINE})
    if(QAWS_SIMD_OPTIONS)
      target_compile_options(${target} PRIVATE ${QAWS_SIMD_OPTIONS})
    endif()
  endif()

  if(QAWS_HAS_NO_STRICT_ALIASING)
    target_compile_options(${target} PRIVATE -fno-strict-aliasing)
  endif()
endfunction()

# Exact/certified CPU geometry (C only; independent of QAWS_SCALAR_IS_FLOAT).
option(QAWS_ENABLE_EXACT "Enable exact/certified CPU geometry" ON)
if(QAWS_ENABLE_EXACT)
  list(APPEND QAWS_SOURCES
    src/qaws/exact/qaws_exact_int.c
    src/qaws/exact/qaws_exact_predicates.c
    src/qaws/exact/qaws_exact_curve.c
    src/qaws/exact/qaws_exact_winding.c
    src/qaws/exact/qaws_exact_surface.c
    src/qaws/exact/qaws_exact_roots.c
    src/qaws/exact/qaws_exact_intersect.c
    src/qaws/exact/qaws_exact_curve_curve.c
    src/qaws/exact/qaws_exact_poly.c
    src/qaws/exact/qaws_exact_solve.c
    src/qaws/exact/qaws_exact_surface_hits.c
    src/qaws/exact/qaws_exact_boolean.c
    src/qaws/exact/qaws_exact_curve_surface.c
    src/qaws/exact/qaws_exact_ssi.c
    src/qaws/exact/qaws_exact_batch.c
    src/qaws/exact/qaws_exact_closest.c
  )
endif()

add_library(qaws ${QAWS_SOURCES})
add_library(Qaws::qaws ALIAS qaws)
if(QAWS_ENABLE_EXACT)
  target_compile_definitions(qaws PUBLIC QAWS_ENABLE_EXACT=1)
endif()

qaws_configure_target(qaws)

set_target_properties(
  qaws
  PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)

if(QAWS_UNITY_BUILD)
  # BATCH_SIZE 0 folds every source into a single translation unit.
  set_target_properties(
    qaws
    PROPERTIES
      UNITY_BUILD ON
      UNITY_BUILD_MODE BATCH
      UNITY_BUILD_BATCH_SIZE 0
  )
endif()

if(QAWS_DEBUG_POSTFIX)
  set_target_properties(qaws PROPERTIES DEBUG_POSTFIX "${QAWS_DEBUG_POSTFIX}")
endif()
if(QAWS_DEBUG_PREFIX)
  set_target_properties(qaws PROPERTIES OUTPUT_NAME_DEBUG "${QAWS_DEBUG_PREFIX}qaws")
endif()

target_include_directories(
  qaws
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/qaws>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# ===================================================================
# Amalgamation: qaws.h + qaws.c, for consumers with no build system.
# ===================================================================

set(QAWS_AMALGAMATION_DIR "${CMAKE_CURRENT_BINARY_DIR}/amalgam"
    CACHE PATH "Directory the generated qaws.h/qaws.c are written to")

set(_QAWS_AMALGAMATE_SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/buildsystem/amalgamate.cmake")
set(_QAWS_AMALGAMATE_ARGS
  "-DQAWS_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/src/qaws"
  "-DQAWS_OUTPUT_DIR=${QAWS_AMALGAMATION_DIR}"
)
if(QAWS_ENABLE_SIMD)
  list(APPEND _QAWS_AMALGAMATE_ARGS "-DQAWS_SIMD=ON")
endif()

set(QAWS_AMALGAMATION_HEADER "${QAWS_AMALGAMATION_DIR}/qaws.h")
set(QAWS_AMALGAMATION_SOURCE "${QAWS_AMALGAMATION_DIR}/qaws.c")

# The generator inlines headers as well as sources, so both must retrigger it.
file(GLOB_RECURSE _QAWS_AMALGAMATE_INPUTS CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/src/qaws/*.c"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/qaws/*.h"
)

add_custom_command(
  OUTPUT ${QAWS_AMALGAMATION_HEADER} ${QAWS_AMALGAMATION_SOURCE}
  COMMAND ${CMAKE_COMMAND} ${_QAWS_AMALGAMATE_ARGS} -P "${_QAWS_AMALGAMATE_SCRIPT}"
  DEPENDS ${_QAWS_AMALGAMATE_SCRIPT} ${_QAWS_AMALGAMATE_INPUTS}
  COMMENT "Generating amalgamation in ${QAWS_AMALGAMATION_DIR}"
  VERBATIM
)

# Convenience target: `cmake --build . --target amalgamate` and nothing else.
add_custom_target(
  amalgamate
  DEPENDS ${QAWS_AMALGAMATION_HEADER} ${QAWS_AMALGAMATION_SOURCE}
)

if(QAWS_AMALGAMATION)
  add_library(qaws_amalgam ${QAWS_AMALGAMATION_SOURCE})
  add_library(Qaws::amalgam ALIAS qaws_amalgam)

  qaws_configure_target(qaws_amalgam)

  # The generated qaws.c includes "qaws.h" from its own directory.
  target_include_directories(
    qaws_amalgam
    PUBLIC $<BUILD_INTERFACE:${QAWS_AMALGAMATION_DIR}>
  )
endif()

option(QAWS_BUILD_EXAMPLES "Build the C example programs" OFF)
if(QAWS_BUILD_EXAMPLES)
  add_executable(qaws_diff_showcase examples/diff_showcase.c)
  target_link_libraries(qaws_diff_showcase PRIVATE Qaws::qaws)
  target_include_directories(qaws_diff_showcase PRIVATE "${PROJECT_SOURCE_DIR}/src/qaws")
  if(MSVC)
    target_compile_definitions(qaws_diff_showcase PRIVATE _CRT_SECURE_NO_WARNINGS)
  endif()
  set_target_properties(qaws_diff_showcase PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/working_dir"
    RUNTIME_OUTPUT_DIRECTORY_RELEASE "${PROJECT_SOURCE_DIR}/working_dir"
    RUNTIME_OUTPUT_DIRECTORY_DEBUG "${PROJECT_SOURCE_DIR}/working_dir")
  add_executable(qaws_diff_applications examples/diff_applications.c)
  target_link_libraries(qaws_diff_applications PRIVATE Qaws::qaws)
  target_include_directories(qaws_diff_applications PRIVATE "${PROJECT_SOURCE_DIR}/src/qaws")
  if(MSVC)
    target_compile_definitions(qaws_diff_applications PRIVATE _CRT_SECURE_NO_WARNINGS)
  endif()
  set_target_properties(qaws_diff_applications PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/working_dir"
    RUNTIME_OUTPUT_DIRECTORY_RELEASE "${PROJECT_SOURCE_DIR}/working_dir"
    RUNTIME_OUTPUT_DIRECTORY_DEBUG "${PROJECT_SOURCE_DIR}/working_dir")
  add_executable(qaws_batch_showcase examples/batch_showcase.c)
  target_link_libraries(qaws_batch_showcase PRIVATE Qaws::qaws)
  target_include_directories(qaws_batch_showcase PRIVATE "${PROJECT_SOURCE_DIR}/src/qaws")
  if(NOT WIN32)
    find_package(Threads REQUIRED)
    target_link_libraries(qaws_batch_showcase PRIVATE Threads::Threads)
  endif()
  if(MSVC)
    target_compile_definitions(qaws_batch_showcase PRIVATE _CRT_SECURE_NO_WARNINGS)
  endif()
  set_target_properties(qaws_batch_showcase PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/working_dir"
    RUNTIME_OUTPUT_DIRECTORY_RELEASE "${PROJECT_SOURCE_DIR}/working_dir"
    RUNTIME_OUTPUT_DIRECTORY_DEBUG "${PROJECT_SOURCE_DIR}/working_dir")
  add_executable(qaws_better2d_showcase examples/better2d_showcase.c)
  target_link_libraries(qaws_better2d_showcase PRIVATE Qaws::qaws)
  target_include_directories(qaws_better2d_showcase PRIVATE "${PROJECT_SOURCE_DIR}/src/qaws")
  if(MSVC)
    target_compile_definitions(qaws_better2d_showcase PRIVATE _CRT_SECURE_NO_WARNINGS)
  endif()
  set_target_properties(qaws_better2d_showcase PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/working_dir"
    RUNTIME_OUTPUT_DIRECTORY_RELEASE "${PROJECT_SOURCE_DIR}/working_dir"
    RUNTIME_OUTPUT_DIRECTORY_DEBUG "${PROJECT_SOURCE_DIR}/working_dir")
  if(QAWS_ENABLE_EXACT)
    add_executable(qaws_exact_showcase examples/exact_showcase.c)
    target_link_libraries(qaws_exact_showcase PRIVATE Qaws::qaws)
    target_include_directories(qaws_exact_showcase PRIVATE "${PROJECT_SOURCE_DIR}/src/qaws")
    if(MSVC)
      target_compile_definitions(qaws_exact_showcase PRIVATE _CRT_SECURE_NO_WARNINGS)
    endif()
    set_target_properties(qaws_exact_showcase PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/working_dir"
      RUNTIME_OUTPUT_DIRECTORY_RELEASE "${PROJECT_SOURCE_DIR}/working_dir"
      RUNTIME_OUTPUT_DIRECTORY_DEBUG "${PROJECT_SOURCE_DIR}/working_dir")
  endif()
endif()


# Clipper2 comparison (optional; Clipper2 is compiled into the benchmark only).
option(QAWS_BUILD_CLIPPER2_BENCH "Build the Clipper2 comparison benchmark" OFF)
set(QAWS_CLIPPER2_DIR "" CACHE PATH "A Clipper2 checkout (the directory holding CPP/Clipper2Lib)")
if(QAWS_BUILD_CLIPPER2_BENCH)
  if(NOT EXISTS "${QAWS_CLIPPER2_DIR}/CPP/Clipper2Lib/src/clipper.engine.cpp")
    message(FATAL_ERROR "QAWS_CLIPPER2_DIR must point at a Clipper2 checkout")
  endif()
  enable_language(CXX)
  add_executable(qaws_clipper2_compare
    bench/clipper2_compare.cpp
    "${QAWS_CLIPPER2_DIR}/CPP/Clipper2Lib/src/clipper.engine.cpp"
    "${QAWS_CLIPPER2_DIR}/CPP/Clipper2Lib/src/clipper.offset.cpp"
    "${QAWS_CLIPPER2_DIR}/CPP/Clipper2Lib/src/clipper.rectclip.cpp")
  target_include_directories(qaws_clipper2_compare PRIVATE
    "${QAWS_CLIPPER2_DIR}/CPP/Clipper2Lib/include" "${PROJECT_SOURCE_DIR}/src/qaws")
  target_link_libraries(qaws_clipper2_compare PRIVATE Qaws::qaws)
  set_target_properties(qaws_clipper2_compare PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/working_dir"
    RUNTIME_OUTPUT_DIRECTORY_RELEASE "${PROJECT_SOURCE_DIR}/working_dir"
    RUNTIME_OUTPUT_DIRECTORY_DEBUG "${PROJECT_SOURCE_DIR}/working_dir")
endif()
if(BUILD_TESTING)
  add_subdirectory(tests)
  add_custom_target(
    check
    COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
    DEPENDS qaws_tests
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
endif()

if(QAWS_INSTALL)
  install(
    TARGETS qaws
    EXPORT QawsTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )

  install(
    FILES
      src/qaws/qaws.h
      src/qaws/qaws_types.h
      src/qaws/qaws_status.h
      src/qaws/qaws_curve.h
      src/qaws/qaws_bezier.h
      src/qaws/qaws_hermite.h
      src/qaws/qaws_catmull_rom.h
      src/qaws/qaws_bspline.h
      src/qaws/qaws_nurbs.h
      src/qaws/qaws_trajectory.h
      src/qaws/qaws_yuksel.h
      src/qaws/qaws_rational_bezier.h
      src/qaws/qaws_composite.h
      src/qaws/qaws_arc.h
      src/qaws/qaws_polynomial.h
      src/qaws/qaws_clothoid.h
      src/qaws/qaws_subdivision.h
      src/qaws/qaws_eval.h
      src/qaws/qaws_sampling.h
      src/qaws/qaws_traversal.h
      src/qaws/qaws_inspect.h
      src/qaws/qaws_curve_batch.h
      src/qaws/qaws_batch_executor.h
      src/qaws/qaws_surface_batch.h
      src/qaws/qaws_operations.h
      src/qaws/qaws_convert.h
      src/qaws/qaws_export.h
      src/qaws/qaws_import.h
      src/qaws/qaws_surface.h
      src/qaws/qaws_surface_types.h
      src/qaws/qaws_surface_bezier.h
      src/qaws/qaws_surface_bspline.h
      src/qaws/qaws_surface_nurbs.h
      src/qaws/qaws_surface_swept.h
      src/qaws/qaws_surface_ruled.h
      src/qaws/qaws_surface_coons.h
      src/qaws/qaws_surface_extrusion.h
      src/qaws/qaws_surface_revolution.h
      src/qaws/qaws_surface_bilinear.h
      src/qaws/qaws_surface_biquadratic.h
      src/qaws/qaws_surface_pipe.h
      src/qaws/qaws_surface_loft.h
      src/qaws/qaws_surface_gordon.h
      src/qaws/qaws_surface_offset.h
      src/qaws/qaws_surface_trim.h
      src/qaws/qaws_surface_intersect.h
      src/qaws/qaws_surface_tspline.h
      src/qaws/qaws_surface_subdiv.h
      src/qaws/qaws_surface_fillet.h
      src/qaws/qaws_brep.h
      src/qaws/qaws_boolean_2d.h
      src/qaws/qaws_path.h
      src/qaws/qaws_clip.h
      src/qaws/qaws_clip64.h
      src/qaws/qaws_offset.h
      src/qaws/qaws_rectclip.h
      src/qaws/qaws_stack.h
      src/qaws/qaws_alloc.h
      src/qaws/qaws_inline.h
      src/qaws/qaws_diff.h
      src/qaws/qaws_diff_types.h
      src/qaws/qaws_diff_geometry.h
      src/qaws/qaws_diff_ops.h
      src/qaws/qaws_diff_map.h
      src/qaws/qaws_diff_functionals.h
      src/qaws/qaws_diff_sampling.h
      src/qaws/qaws_platform.h
      src/qaws/qaws_core_types.h
      src/qaws/qaws_prepare.h
      src/qaws/qaws_hlsl.h
      src/qaws/qaws_glsl.h
      src/qaws/qaws_halide.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )

  install(
    FILES
      src/qaws/core/qaws_decasteljau_core.h
      src/qaws/core/qaws_cubic_poly_core.h
      src/qaws/core/qaws_bspline_basis_core.h
      src/qaws/core/qaws_bspline_eval_core.h
      src/qaws/core/qaws_nurbs_eval_core.h
      src/qaws/core/qaws_rational_bezier_core.h
      src/qaws/core/qaws_horner_core.h
      src/qaws/core/qaws_arc_core.h
      src/qaws/core/qaws_surface_core.h
      src/qaws/core/qaws_dual_core.h
      src/qaws/core/qaws_bspline_diff_core.h
      src/qaws/core/qaws_bspline_sampling_core.h
      src/qaws/core/qaws_bspline_surface_functional_core.h
      src/qaws/core/qaws_bspline_surface_sampling_core.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/core
  )

  set(QAWS_CMAKE_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/Qaws")

  install(
    EXPORT QawsTargets
    FILE QawsTargets.cmake
    NAMESPACE Qaws::
    DESTINATION ${QAWS_CMAKE_INSTALL_DIR}
  )

  configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/buildsystem/cmake/QawsConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/QawsConfig.cmake"
    INSTALL_DESTINATION ${QAWS_CMAKE_INSTALL_DIR}
  )

  write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/QawsConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
  )

  install(
    FILES
      "${CMAKE_CURRENT_BINARY_DIR}/QawsConfig.cmake"
      "${CMAKE_CURRENT_BINARY_DIR}/QawsConfigVersion.cmake"
    DESTINATION ${QAWS_CMAKE_INSTALL_DIR}
  )

  if(QAWS_INSTALL_PKGCONFIG)
    set(_QAWS_PKGCONFIG_LIBS_PRIVATE "")
    if(NOT WIN32)
      set(_QAWS_PKGCONFIG_LIBS_PRIVATE "-lm")
    endif()

    configure_file(
      "${CMAKE_CURRENT_SOURCE_DIR}/buildsystem/cmake/qaws.pc.in"
      "${CMAKE_CURRENT_BINARY_DIR}/qaws.pc"
      @ONLY
    )

    install(
      FILES "${CMAKE_CURRENT_BINARY_DIR}/qaws.pc"
      DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig"
    )
  endif()
endif()
