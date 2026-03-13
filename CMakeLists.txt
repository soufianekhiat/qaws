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
  src/qaws/qaws_bezier.c
  src/qaws/qaws_hermite.c
  src/qaws/qaws_catmull_rom.c
  src/qaws/qaws_bspline.c
  src/qaws/qaws_nurbs.c
  src/qaws/qaws_trajectory.c
  src/qaws/qaws_yuksel.c
  src/qaws/qaws_operations.c
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
  src/qaws/qaws_inline.c
  src/qaws/qaws_prepare.c
  src/qaws/internal/qaws_internal_curve.c
  src/qaws/internal/qaws_internal_basis.c
  src/qaws/internal/qaws_internal_span.c
  src/qaws/internal/qaws_internal_arc_length.c
  src/qaws/internal/qaws_internal_solver.c
  src/qaws/internal/qaws_internal_validation.c
  src/qaws/internal/qaws_internal_surface.c
)

if(QAWS_ENABLE_SIMD)
  list(APPEND QAWS_SOURCES
    src/qaws/batch/qaws_batch_bezier.c
    src/qaws/batch/qaws_batch_bspline.c
    src/qaws/batch/qaws_batch_nurbs.c
    src/qaws/batch/qaws_batch_surface.c
  )
endif()

add_library(qaws ${QAWS_SOURCES})
add_library(Qaws::qaws ALIAS qaws)

set_target_properties(
  qaws
  PROPERTIES
    C_STANDARD 11
    C_STANDARD_REQUIRED YES
    C_EXTENSIONS NO
    POSITION_INDEPENDENT_CODE ON
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)

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

target_compile_definitions(
  qaws
  PUBLIC
    QAWS_SCALAR_IS_FLOAT=${_QAWS_SCALAR_IS_FLOAT_VALUE}
)

if(MSVC)
  target_compile_definitions(qaws PRIVATE _CRT_SECURE_NO_WARNINGS NOMINMAX)
endif()

if(NOT WIN32)
  target_link_libraries(qaws PUBLIC m)
endif()

if(WIN32 AND BUILD_SHARED_LIBS)
  set_target_properties(qaws PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

if(QAWS_ENABLE_WARNINGS)
  if(MSVC)
    target_compile_options(qaws PRIVATE /W4)
    if(QAWS_WARNINGS_AS_ERRORS)
      target_compile_options(qaws PRIVATE /WX)
    endif()
  else()
    target_compile_options(qaws PRIVATE -Wall -Wextra -Wpedantic)
    if(QAWS_WARNINGS_AS_ERRORS)
      target_compile_options(qaws PRIVATE -Werror)
    endif()
  endif()
endif()

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
  if(QAWS_HAS_AVX2)
    target_compile_definitions(qaws PRIVATE QAWS_SIMD_AVX2=1)
    if(MSVC)
      target_compile_options(qaws PRIVATE /arch:AVX2)
    else()
      target_compile_options(qaws PRIVATE -mavx2 -mfma)
    endif()
  elseif(QAWS_HAS_SSE2)
    target_compile_definitions(qaws PRIVATE QAWS_SIMD_SSE2=1)
  elseif(QAWS_HAS_NEON)
    target_compile_definitions(qaws PRIVATE QAWS_SIMD_NEON=1)
  endif()
endif()

if(QAWS_DISABLE_STRICT_ALIASING AND NOT MSVC)
  check_c_compiler_flag(-fno-strict-aliasing QAWS_HAS_NO_STRICT_ALIASING)
  if(QAWS_HAS_NO_STRICT_ALIASING)
    target_compile_options(qaws PRIVATE -fno-strict-aliasing)
  endif()
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
      src/qaws/qaws_alloc.h
      src/qaws/qaws_inline.h
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
