# Copyright 2025 Matt Borland
# Distributed under the Boost Software License, Version 1.0.
# https://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 3.8...3.16)

project(cmake_install_test LANGUAGES CXX)

find_package(boost_math REQUIRED)

add_executable(quick ../quick.cpp)
target_link_libraries(quick Boost::math)

enable_testing()
add_test(quick quick)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -C $<CONFIG>)
