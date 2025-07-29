cmake_minimum_required(VERSION 3.5)
set(CMAKE_EXPORT_COMPILE_COMMANDS YES)
set(THIS autodiff_reverse)

project(autodiff_reverse C CXX)
set(CMAKE_CXX_STANDARD 14)

enable_testing()

set(BOOST_ROOT_DIRECTORY /home/mzhelyez/Documents/code/boost/)
set(BOOST_MATH_DIRECTORY ${BOOST_ROOT_DIRECTORY}/libs/math/)
set(BOOST_MATH_TEST_DIRECTORY ${BOOST_MATH_DIRECTORY}test/)

include_directories(${BOOST_MATH_DIRECTORY}include)


set(HEADERS
    ${CMAKE_SOURCE_DIR}/autodiff_reverse.hpp
    ${BOOST_MATH_TEST_DIRECTORY}test_autodiff_reverse.hpp
)

set(TEST_MEMORY_ALLOCATORS
    ${BOOST_MATH_TEST_DIRECTORY}test_autodiff_reverse_flat_linear_allocator.cpp
)

set(TEST_CONSTRUCTORS
    ${BOOST_MATH_TEST_DIRECTORY}test_rvar_constructors.cpp
)
# add_executable(autodiff_reverse_tests ${SOURCES} ${HEADERS})
add_executable(autodiff_reverse_tests ${TEST_CONSTRUCTORS} ${HEADERS})
# add_executable(autodiff_reverse_tests ${TEST_MEMORY_ALLOCATORS} ${HEADERS})


