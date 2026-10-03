# Copyright 2026 Matt Borland
# Distributed under the Boost Software License, Version 1.0.
# https://www.boost.org/LICENSE_1_0.txt
#
# Validates the configurable top level namespace (https://github.com/boostorg/math/issues/769)
# by relocating the whole library into ::my_lib::math.

cmake_minimum_required(VERSION 3.8...3.20)

project(cmake_namespace_test LANGUAGES CXX)

option(BOOST_MATH_NAMESPACE_TEST_HEADERS "Compile every public header on its own with the relocated namespace" ON)

# Math is header-only and standalone by default (no deps needed)
add_subdirectory(../.. boostorg/math)

set(namespace_definition "BOOST_MATH_NAMESPACE_COMPONENTS=my_lib,math")

enable_testing()

# Behavior through ::my_lib::math
add_executable(namespace_test namespace_test.cpp)
target_link_libraries(namespace_test Boost::math)
target_compile_definitions(namespace_test PRIVATE ${namespace_definition})
add_test(NAME namespace_test COMMAND namespace_test)

# BOOST_MATH_NAMESPACE::name and my_lib::math::name are the same entities
add_executable(qualified_name_test qualified_name_test.cpp)
target_link_libraries(qualified_name_test Boost::math)
target_compile_definitions(qualified_name_test PRIVATE ${namespace_definition})
add_test(NAME qualified_name_test COMMAND qualified_name_test)

# No symbol of the test binary may be mangled into ::boost
if(CMAKE_NM AND NOT MSVC)
  add_test(NAME namespace_symbols COMMAND ${CMAKE_COMMAND} -DNM=${CMAKE_NM} -DFILE=$<TARGET_FILE:namespace_test> -P ${CMAKE_CURRENT_SOURCE_DIR}/check_symbols.cmake)
endif()

if(BOOST_MATH_NAMESPACE_TEST_HEADERS)

  set(include_root "${CMAKE_CURRENT_SOURCE_DIR}/../../include")
  file(GLOB_RECURSE headers RELATIVE "${include_root}" "${include_root}/boost/math/*.hpp")
  list(FILTER headers EXCLUDE REGEX "/detail/")

  # Headers that need third party libraries, need non-standalone Boost, or are fragments of another header
  list(FILTER headers EXCLUDE REGEX "^boost/math/bindings/")
  list(REMOVE_ITEM headers
    boost/math/common_factor.hpp
    boost/math/common_factor_ct.hpp
    boost/math/common_factor_rt.hpp
    boost/math/concepts/distributions.hpp
    boost/math/constants/calculate_constants.hpp
    boost/math/cstdfloat/cstdfloat_complex_std.hpp
    boost/math/differentiation/autodiff_cpp11.hpp
    boost/math/interpolators/cardinal_trigonometric.hpp
    boost/math/optimization/cma_es.hpp
    boost/math/special_functions/chebyshev_transform.hpp
    boost/math/tools/signal_statistics.hpp
  )

  set(header_sources "")
  foreach(header IN LISTS headers)
    string(MAKE_C_IDENTIFIER "${header}" name)
    set(source "${CMAKE_CURRENT_BINARY_DIR}/headers/${name}.cpp")
    configure_file(header_test.cpp.in "${source}" @ONLY)
    list(APPEND header_sources "${source}")
  endforeach()

  add_library(namespace_header_tests OBJECT ${header_sources})
  target_link_libraries(namespace_header_tests PRIVATE Boost::math)
  target_compile_definitions(namespace_header_tests PRIVATE ${namespace_definition})
  target_compile_features(namespace_header_tests PRIVATE cxx_std_17)

endif()

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -C $<CONFIG>)
