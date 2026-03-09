# Copyright 2025 Matt Borland
# Distributed under the Boost Software License, Version 1.0.
# https://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 3.8...3.20)

project(cmake_subdir_test_legacy LANGUAGES CXX)

# config must be added before math (compiled targets depend on Boost::config)
add_subdirectory(../../../config boostorg/config)

# Enable legacy compiled targets before adding math
set(BOOST_MATH_BUILD_WITH_LEGACY_FUNCTIONS ON CACHE BOOL "" FORCE)

add_subdirectory(../.. boostorg/math)

add_executable(quick_tr1 ../quick_tr1.cpp)
target_link_libraries(quick_tr1 Boost::math Boost::math_c99 Boost::math_tr1)

enable_testing()
add_test(quick_tr1 quick_tr1)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -C $<CONFIG>)
