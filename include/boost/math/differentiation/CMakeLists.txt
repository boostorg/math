cmake_minimum_required(VERSION 3.5)
set(CMAKE_EXPORT_COMPILE_COMMANDS YES)
set(THIS autodiff_reverse)

project(autodiff_reverse C CXX)
set(CMAKE_CXX_STANDARD 14)

enable_testing()#
# add_compile_options(-fsanitize=address)
# add_link_options(-fsanitize=address)
set(BOOST_ROOT_DIRECTORY /home/mzhelyez/Documents/code/boost/)
set(BOOST_MATH_DIRECTORY ${BOOST_ROOT_DIRECTORY}/libs/math/)
set(BOOST_MATH_TEST_DIRECTORY ${BOOST_MATH_DIRECTORY}test/)

set(BOOST_MATH_EXAMPLE_DIRECTORY ${BOOST_MATH_DIRECTORY}/example)
include_directories(${BOOST_MATH_DIRECTORY}include
                    ${BOOST_ROOT_DIRECTORY}
)


set(HEADERS
    ${CMAKE_SOURCE_DIR}/autodiff_reverse.hpp
    ${CMAKE_SOURCE_DIR}/detail/reverse_mode_autodiff_basic_operator_overloads.hpp
    ${CMAKE_SOURCE_DIR}/detail/reverse_mode_autodiff_helper_functions.hpp
    ${CMAKE_SOURCE_DIR}/detail/reverse_mode_autodiff_memory_management.hpp
    ${CMAKE_SOURCE_DIR}/detail/reverse_mode_autodiff_stl_overloads.hpp
    ${CMAKE_SOURCE_DIR}/detail/reverse_mode_autodiff_expression_template_base.hpp
    ${CMAKE_SOURCE_DIR}/detail/reverse_mode_autodiff_utilities.hpp
    ${CMAKE_SOURCE_DIR}/detail/reverse_mode_autodiff_erf_overloads.hpp
    ${CMAKE_SOURCE_DIR}/detail/reverse_mode_autodiff_comparison_operator_overloads.hpp
    ${BOOST_MATH_TEST_DIRECTORY}test_autodiff_reverse.hpp
)

set(TEST_MEMORY_ALLOCATORS
    ${BOOST_MATH_TEST_DIRECTORY}test_reverse_mode_autodiff_flat_linear_allocator.cpp
)

set(TEST_CONSTRUCTORS
    ${BOOST_MATH_TEST_DIRECTORY}test_reverse_mode_autodiff_constructors.cpp
)

set(TEST_COMPARISON_OPS
    ${BOOST_MATH_TEST_DIRECTORY}test_reverse_mode_autodiff_comparison_operators.cpp
)

set(TEST_STL_SUPPORT
    ${BOOST_MATH_TEST_DIRECTORY}test_reverse_mode_autodiff_stl_support.cpp
)
SET(TEST_BASIC_OPS ${BOOST_MATH_TEST_DIRECTORY}test_reverse_mode_autodiff_basic_math_ops.cpp)
SET(TEST_ADDITIONAL_BOOST_FUNCTIONS ${BOOST_MATH_TEST_DIRECTORY}test_autodiff_reverse_additional_boost_functions.cpp)
SET(LINEAR_REGRESSION_EXAMPLE ${BOOST_MATH_EXAMPLE_DIRECTORY}/reverse_mode_linear_regression_example.cpp)
set(TEST_ERF ${BOOST_MATH_TEST_DIRECTORY}/test_reverse_mode_autodiff_error_functions.cpp)
set(BLACK_SCHOLES_EXAMPLE ${BOOST_MATH_EXAMPLE_DIRECTORY}/autodiff_reverse_black_scholes.cpp)
# add_executable(autodiff_reverse_tests ${TEST_CONSTRUCTORS} ${HEADERS})
# add_executable(autodiff_reverse_tests ${TEST_MEMORY_ALLOCATORS} ${HEADERS})
# add_executable(autodiff_reverse_tests ${TEST_BASIC_OPS} ${HEADERS})
# add_executable(autodiff_reverse_tests ${TEST_COMPARISON_OPS} ${HEADERS})
add_executable(autodiff_reverse_tests ${TEST_STL_SUPPORT} ${HEADERS})
# add_executable(autodiff_reverse_tests ${LINEAR_REGRESSION_EXAMPLE} ${HEADERS})
# add_executable(autodiff_reverse_tests ${TEST_ERF} ${HEADERS})
# add_executable(autodiff_reverse_tests ${BLACK_SCHOLES_EXAMPLE} ${HEADERS})




