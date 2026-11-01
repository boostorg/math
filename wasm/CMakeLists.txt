# Copyright 2026 Nick Thompson
# Distributed under the Boost Software License, Version 1.0.
# https://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 3.16)
project(boost_math_wasm LANGUAGES CXX)

if(NOT EMSCRIPTEN)
  message(FATAL_ERROR "Configure with emcmake")
endif()

option(BOOST_MATH_WASM_RELAXED_SIMD "Use relaxed SIMD (not supported by Safari)" OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXECUTABLE_SUFFIX .mjs)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_executable(boost_math bindings.cpp)
target_include_directories(boost_math PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)
target_compile_definitions(boost_math PRIVATE BOOST_MATH_STANDALONE)
# compensated_arithmetic.hpp's #pragma float_control is unknown to wasm.
target_compile_options(boost_math PRIVATE -Wno-ignored-pragmas)
target_link_options(boost_math PRIVATE -lembind -sALLOW_MEMORY_GROWTH)
if(BOOST_MATH_WASM_RELAXED_SIMD)
  # The SLP vectorizer's 128-bit code runs slower in V8 than scalar code.
  target_compile_options(boost_math PRIVATE -mrelaxed-simd -fno-slp-vectorize)
endif()

