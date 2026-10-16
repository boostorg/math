# Copyright 2021 Matt Borland
# Distributed under the Boost Software License, Version 1.0.
# https://www.boost.org/LICENSE_1_0.txt

include_directories(../../include_private)
file(GLOB SRC_FILES CONFIGURE_DEPENDS "*.cpp")
# Users compile these headers with warnings as errors, so we do too. -Wno-deprecated-declarations
# because some of these tests exercise deprecated functions on purpose.
# gcc's plain -Wshadow also flags header locals that shadow the single-letter globals of
# test_compile_result.hpp, so gcc gets -Wshadow=local (locals and parameters only).
set(BOOST_MATH_COMPILE_TEST_WARNINGS "")
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(BOOST_MATH_COMPILE_TEST_WARNINGS -Wall -Wextra -Wshadow=local -Werror -Wno-deprecated-declarations)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(BOOST_MATH_COMPILE_TEST_WARNINGS -Wall -Wextra -Wshadow -Werror -Wno-deprecated-declarations)
endif()
boost_test(TYPE "compile" SOURCES ${SRC_FILES} COMPILE_DEFINITIONS BOOST_MATH_STANDALONE COMPILE_OPTIONS ${BOOST_MATH_COMPILE_TEST_WARNINGS} COMPILE_FEATURES cxx_std_17 LINK_LIBRARIES Boost::math Boost::multiprecision Boost::numeric_ublas Boost::unit_test_framework )
