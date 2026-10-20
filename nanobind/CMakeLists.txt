# Copyright 2026 Nick Thompson
# Distributed under the Boost Software License, Version 1.0.
# https://www.boost.org/LICENSE_1_0.txt
#
# Dev-only nanobind layer, not built by CI. Usage (venv with nanobind numpy mpmath pytest):
#   cmake -S nanobind -B build-nb -G Ninja -DPython_EXECUTABLE=venv/bin/python && cmake --build build-nb
#   PYTHONPATH=build-nb venv/bin/python -m pytest nanobind/tests

cmake_minimum_required(VERSION 3.18)
project(boost_math_nanobind LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
# Same optimisation level as CI.
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

find_package(Python 3.9 REQUIRED COMPONENTS Interpreter Development.Module)
if(NOT nanobind_DIR)
  # -I: don't let this directory (also called nanobind) shadow the package.
  execute_process(
    COMMAND "${Python_EXECUTABLE}" -I -m nanobind --cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_DIR
    COMMAND_ERROR_IS_FATAL ANY)
endif()
find_package(nanobind CONFIG REQUIRED)

set(BOOST_MATH_NB_SOURCES
  src/module.cpp
  src/gamma.cpp
  src/beta.cpp
  src/erf_expint_zeta.cpp
  src/elementary.cpp
  src/bessel.cpp
  src/elliptic.cpp
  src/jacobi_theta.cpp
  src/polynomials.cpp
  src/hypergeometric.cpp
  src/distributions_continuous.cpp
  src/distributions_heavy.cpp
  src/distributions_discrete.cpp
  src/distributions_noncentral.cpp
  src/fft.cpp
  src/quadrature.cpp
  src/optimization.cpp
  src/polynomial_roots.cpp
)

nanobind_add_module(_ext NOMINSIZE ${BOOST_MATH_NB_SOURCES})
target_include_directories(_ext PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../include" "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_compile_definitions(_ext PRIVATE BOOST_MATH_STANDALONE)
# Optional: cma_es needs Eigen.
find_package(Eigen3 CONFIG QUIET)
if(Eigen3_FOUND)
  target_link_libraries(_ext PRIVATE Eigen3::Eigen)
endif()
# Optional: boost.math.ulps needs Boost.Multiprecision (standalone, which needs Boost.Config) and zlib. Only the
# sibling checkouts (or the cache variables) are searched, never a system Boost of another release.
find_path(BOOST_MATH_NB_MULTIPRECISION_INCLUDE boost/multiprecision/cpp_bin_float.hpp
  HINTS "${CMAKE_CURRENT_SOURCE_DIR}/../../multiprecision/include" "$ENV{HOME}/boost/libs/multiprecision/include"
  DOC "Boost.Multiprecision's include directory, for boost.math.ulps" NO_DEFAULT_PATH)
find_path(BOOST_MATH_NB_CONFIG_INCLUDE boost/config.hpp
  HINTS "${CMAKE_CURRENT_SOURCE_DIR}/../../config/include" "$ENV{HOME}/boost/libs/config/include"
  DOC "Boost.Config's include directory, for boost.math.ulps" NO_DEFAULT_PATH)
find_package(ZLIB)
if(BOOST_MATH_NB_MULTIPRECISION_INCLUDE AND BOOST_MATH_NB_CONFIG_INCLUDE AND ZLIB_FOUND)
  set(ulps_sources src/ulps.cpp src/ulps_gamma.cpp src/ulps_special.cpp src/ulps_elliptic.cpp)
  target_sources(_ext PRIVATE ${ulps_sources})
  # Only these see Multiprecision, so that its presence changes nothing else:
  set_source_files_properties(${ulps_sources} PROPERTIES
    INCLUDE_DIRECTORIES "${BOOST_MATH_NB_MULTIPRECISION_INCLUDE};${BOOST_MATH_NB_CONFIG_INCLUDE}"
    COMPILE_DEFINITIONS BOOST_MP_STANDALONE)
  set_source_files_properties(src/module.cpp PROPERTIES COMPILE_DEFINITIONS BMPY_HAVE_ULPS)
  target_link_libraries(_ext PRIVATE ZLIB::ZLIB)
else()
  message(STATUS "boost.math.ulps disabled: needs Boost.Multiprecision, Boost.Config and zlib "
                 "(set BOOST_MATH_NB_MULTIPRECISION_INCLUDE and BOOST_MATH_NB_CONFIG_INCLUDE)")
endif()
if(NOT MSVC)
  # Match CI and test/Makefile: contraction allowed everywhere.
  target_compile_options(_ext PRIVATE -ffp-contract=fast)
endif()

set(pkg_dir "${CMAKE_CURRENT_BINARY_DIR}/boost/math")
set_target_properties(_ext PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${pkg_dir}")
# The pure-Python files, copied on every build so that edits to them need no relink:
add_custom_target(boost_math_py ALL
  COMMAND "${CMAKE_COMMAND}" -E make_directory "${pkg_dir}"
  COMMAND "${CMAKE_COMMAND}" -E copy_if_different
          "${CMAKE_CURRENT_SOURCE_DIR}/boost/math/__init__.py" "${CMAKE_CURRENT_SOURCE_DIR}/boost/math/ulps.py" "${pkg_dir}")
