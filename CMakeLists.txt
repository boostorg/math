# Copyright 2020 Peter Dimov
# Copyright 2021 Matt Borland
# Distributed under the Boost Software License, Version 1.0.
# https://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 3.5)

project(boost_math VERSION "${BOOST_SUPERPROJECT_VERSION}" LANGUAGES CXX)

add_library(boost_math INTERFACE)

add_library(Boost::math ALIAS boost_math)

target_include_directories(boost_math INTERFACE include)

# Default to standalone mode if using CMake
if(EXISTS BOOST_MATH_STANDALONE)
else()
  add_compile_definitions(BOOST_MATH_STANDALONE=1)
endif()

if(!BOOST_MATH_STANDALONE)
  target_link_libraries(boost_math
    INTERFACE
      Boost::assert
      Boost::concept_check
      Boost::config
      Boost::core
      Boost::integer
      Boost::lexical_cast
      Boost::predef
      Boost::random
      Boost::static_assert
      Boost::throw_exception
  )
endif()

if(BUILD_TESTING AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/test/CMakeLists.txt")

  add_subdirectory(test)

endif()

if(BUILD_EXAMPLE AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/examples/CMakeLists.txt")

  add_subdirectory(example)

endif()
