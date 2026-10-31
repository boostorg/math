# Copyright 2026 Nick Thompson
# Distributed under the Boost Software License, Version 1.0.
# https://www.boost.org/LICENSE_1_0.txt
#
# The nanobind layer; .github/workflows/wheels.yml builds its wheels. Usage (venv with nanobind numpy mpmath pytest):
#   cmake -S nanobind -B build-nb -G Ninja -DPython_EXECUTABLE=venv/bin/python && cmake --build build-nb
#   PYTHONPATH=build-nb venv/bin/python -m pytest nanobind/tests
# A wheel, from a checkout (an sdist would lack ../include): pip wheel --no-deps -w dist ./nanobind

cmake_minimum_required(VERSION 3.18)
project(boost_math_nanobind LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# -std=c++20, not gnu++20: in GNU mode Boost.Config enables __float128, which would need libquadmath.
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
# Same optimisation level as CI.
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

find_package(Python 3.12 REQUIRED COMPONENTS Interpreter Development.Module OPTIONAL_COMPONENTS Development.SABIModule)
# scikit-build-core puts nanobind on the prefix path itself.
if(NOT nanobind_DIR AND NOT SKBUILD)
  # -I: don't let this directory (also called nanobind) shadow the package.
  execute_process(
    COMMAND "${Python_EXECUTABLE}" -I -m nanobind --cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_DIR
    COMMAND_ERROR_IS_FATAL ANY)
endif()
find_package(nanobind CONFIG REQUIRED)

set(BOOST_MATH_NB_SOURCES
  src/module.cpp
  src/errors.cpp
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
  src/roots.cpp
  src/statistics.cpp
  src/interpolators.cpp
  src/special_extras.cpp
  src/daubechies.cpp
  src/constants.cpp
  src/differentiation.cpp
  src/expansions.cpp
  src/naive_monte_carlo.cpp
)

nanobind_add_module(_ext NOMINSIZE STABLE_ABI ${BOOST_MATH_NB_SOURCES})
target_include_directories(_ext PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../include" "${CMAKE_CURRENT_SOURCE_DIR}/src")
# Errors go to the handlers in src/errors.cpp, which serve boost.math.errstate. Every source, the ulps ones too:
# a policy<> that differed between sources would be an ODR violation.
target_compile_definitions(_ext PRIVATE BOOST_MATH_STANDALONE
  BOOST_MATH_DOMAIN_ERROR_POLICY=user_error BOOST_MATH_POLE_ERROR_POLICY=user_error
  BOOST_MATH_OVERFLOW_ERROR_POLICY=user_error BOOST_MATH_EVALUATION_ERROR_POLICY=user_error
  BOOST_MATH_ROUNDING_ERROR_POLICY=user_error BOOST_MATH_INDETERMINATE_RESULT_ERROR_POLICY=user_error)
# boost.math.ulps needs Boost.Multiprecision (standalone, which needs Boost.Config): optional, except in a wheel. Only
# the sibling checkouts (or the cache variables) are searched, never a system Boost of another release.
find_path(BOOST_MATH_NB_MULTIPRECISION_INCLUDE boost/multiprecision/cpp_bin_float.hpp
  HINTS "${CMAKE_CURRENT_SOURCE_DIR}/../../multiprecision/include" "$ENV{HOME}/boost/libs/multiprecision/include"
  DOC "Boost.Multiprecision's include directory, for boost.math.ulps" NO_DEFAULT_PATH)
find_path(BOOST_MATH_NB_CONFIG_INCLUDE boost/config.hpp
  HINTS "${CMAKE_CURRENT_SOURCE_DIR}/../../config/include" "$ENV{HOME}/boost/libs/config/include"
  DOC "Boost.Config's include directory, for boost.math.ulps" NO_DEFAULT_PATH)
if(BOOST_MATH_NB_MULTIPRECISION_INCLUDE AND BOOST_MATH_NB_CONFIG_INCLUDE)
  set(ulps_sources src/ulps.cpp src/ulps_gamma.cpp src/ulps_special.cpp src/ulps_elliptic.cpp)
  target_sources(_ext PRIVATE ${ulps_sources})
  # Only these see Multiprecision, so that its presence changes nothing else (errors.cpp instantiates the error
  # handlers for the reference type):
  set_source_files_properties(${ulps_sources} src/errors.cpp PROPERTIES
    INCLUDE_DIRECTORIES "${BOOST_MATH_NB_MULTIPRECISION_INCLUDE};${BOOST_MATH_NB_CONFIG_INCLUDE}"
    COMPILE_DEFINITIONS BOOST_MP_STANDALONE)
  set_source_files_properties(src/module.cpp PROPERTIES COMPILE_DEFINITIONS BMPY_HAVE_ULPS)
  set_property(SOURCE src/errors.cpp APPEND PROPERTY COMPILE_DEFINITIONS BMPY_HAVE_ULPS)
  install(FILES boost/math/ulps.py DESTINATION boost/math)
else()
  # A wheel without it would be published as complete, so it is an error there:
  if(SKBUILD)
    set(ulps_missing FATAL_ERROR)
  else()
    set(ulps_missing STATUS)
  endif()
  message(${ulps_missing} "boost.math.ulps disabled: needs Boost.Multiprecision and Boost.Config "
                 "(set BOOST_MATH_NB_MULTIPRECISION_INCLUDE and BOOST_MATH_NB_CONFIG_INCLUDE)")
endif()
if(SKBUILD)
  install(FILES ../LICENSE DESTINATION "${SKBUILD_METADATA_DIR}/licenses")
endif()
# MSVC needs only /utf-8, for the docstrings: nanobind adds /bigobj, and CMake /EHsc.
if(MSVC)
  target_compile_options(_ext PRIVATE /utf-8)
else()
  # Match CI and test/Makefile: contraction allowed everywhere.
  target_compile_options(_ext PRIVATE -ffp-contract=fast)
endif()

set(pkg_dir "${CMAKE_CURRENT_BINARY_DIR}/boost/math")
# $<1:...> keeps multi-config generators (Visual Studio) from adding a Release/ subdirectory.
set_target_properties(_ext PROPERTIES LIBRARY_OUTPUT_DIRECTORY "$<1:${pkg_dir}>")
# The pure-Python files, copied on every build so that edits to them need no relink:
add_custom_target(boost_math_py ALL
  COMMAND "${CMAKE_COMMAND}" -E make_directory "${pkg_dir}"
  COMMAND "${CMAKE_COMMAND}" -E copy_if_different
          "${CMAKE_CURRENT_SOURCE_DIR}/boost/math/__init__.py" "${CMAKE_CURRENT_SOURCE_DIR}/boost/math/ulps.py" "${pkg_dir}")
install(TARGETS _ext LIBRARY DESTINATION boost/math)
install(FILES boost/math/__init__.py DESTINATION boost/math)
# Type stubs, boost/math/<module>.pyi, and py.typed, generated from the module just built. The submodules exist only in
# sys.modules, so each has its own stubgen run.
set(marker MARKER_FILE "${pkg_dir}/py.typed")
set(stubs "${pkg_dir}/py.typed")
foreach(module _ext distributions quadrature optimize fft roots statistics interpolators constants differentiation expansions)
  nanobind_add_stub(boost_math_stub_${module} MODULE boost.math.${module} OUTPUT "${pkg_dir}/${module}.pyi" ${marker}
    PYTHON_PATH "${CMAKE_CURRENT_BINARY_DIR}" PATTERN_FILE "${CMAKE_CURRENT_SOURCE_DIR}/stub_patterns.txt"
    DEPENDS _ext boost_math_py)
  list(APPEND stubs "${pkg_dir}/${module}.pyi")
  set(marker "")
endforeach()
install(FILES ${stubs} DESTINATION boost/math)
