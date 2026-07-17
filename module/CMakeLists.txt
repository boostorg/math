# Copyright 2026 Matt Borland
# Distributed under the Boost Software License, Version 1.0.
# https://www.boost.org/LICENSE_1_0.txt
#
# Builds Boost.Math as a C++ module (boost.math) and runs the module-aware
# tests in ../test against it, consuming the standard library module
# (import std). This file is included from ../test/CMakeLists.txt when the
# CMake option BOOST_MATH_BUILD_MODULE is set and expects the Boost
# superproject (for the Boost:: dependency targets). import std is required
# (see the check below).

if(NOT TARGET tests)
    add_custom_target(tests)
endif()

# The module interface unit.
add_library(boost_math_module)
target_sources(boost_math_module
        PUBLIC
        FILE_SET math_module TYPE CXX_MODULES FILES math.cppm
)
target_link_libraries(boost_math_module PUBLIC Boost::math)
target_compile_features(boost_math_module PUBLIC cxx_std_23)

# The purview includes the library's angled <boost/math/...> headers by design.
# clang-scan-deps warns on that and ignores -Wno-include-angled-in-module-purview,
# so -w keeps the dependency-scan output clean (this glue unit's code is fully
# warned on in the header builds).
target_compile_options(boost_math_module PRIVATE $<$<CXX_COMPILER_ID:Clang>:-w>)

# math.cppm defines BOOST_MATH_BUILD_MODULE and BOOST_MATH_STANDALONE itself,
# but the tests link this target and need both macros too: BOOST_MATH_BUILD_MODULE
# so they import the module instead of including the headers directly, and
# BOOST_MATH_STANDALONE so any textually included support header (for example
# test/math_unit_test.hpp -> tools/config.hpp) follows the same configuration
# branch the module was built with.
target_compile_definitions(boost_math_module PUBLIC BOOST_MATH_BUILD_MODULE BOOST_MATH_STANDALONE=1)

# Export the detail-namespace entities that the tests exercise directly.
# Only the interface unit needs this, so keep it private to the module target.
target_compile_definitions(boost_math_module PRIVATE BOOST_MATH_EXPORT_TESTING)

# Mixing textual standard library includes in consumers with a module built from
# header units trips include-translation bugs, so the standard library module is
# required for the test build. CMAKE_CXX_COMPILER_IMPORT_STD lists the standards
# it is available for; it is only populated when CMAKE_EXPERIMENTAL_CXX_IMPORT_STD
# is set before the project() call.
if(NOT "23" IN_LIST CMAKE_CXX_COMPILER_IMPORT_STD)
    message(FATAL_ERROR
            "Building the Boost.Math module tests requires 'import std' (the C++23 "
            "standard library module). Use a toolchain that supports it and pass "
            "-DCMAKE_EXPERIMENTAL_CXX_IMPORT_STD=<uuid> before the project() call. "
            "CMAKE_CXX_COMPILER_IMPORT_STD='${CMAKE_CXX_COMPILER_IMPORT_STD}'")
endif()

set_target_properties(boost_math_module PROPERTIES CXX_MODULE_STD ON)
target_compile_definitions(boost_math_module PRIVATE BOOST_MATH_USE_STD_MODULE)
message(STATUS "Boost.Math module: using 'import std;'")

# Tests that switch to `import boost.math;` when BOOST_MATH_BUILD_MODULE is set.
# The list grows as components are added to the module interface unit.
set(BOOST_MATH_MODULE_TESTS
        test_digamma_simple
        test_expm1_simple
        test_log1p_simple
        git_issue_1194
        git_issue_1255
        git_issue_1294
        git_issue_800
        git_issue_1120
        scipy_issue_17146
        ccmath_isinf_test
        ccmath_isnan_test
        git_issue_1139
        cubic_roots_test
        quartic_roots_test
        git_issue_898
        cardinal_quadratic_b_spline_test
        whittaker_shannon_test
        test_runs_test
        linear_regression_test
        ljung_box_test
        differential_evolution_test
        random_search_test
        quaternion_constexpr_test
        octonion_test_simple
        lambert_w_module_test
        students_t_df_module_test
        ccmath_next_module_test
)

# The full consumer suite runs under both clang/libc++ and GCC/libstdc++: the
# test sources take the standard library from `import std;` in module mode (see
# math_unit_test.hpp), so they no longer mix textual std includes with the
# imported std that libstdc++ rejects. Running the whole list on GCC matters
# because GCC enforces the module exposure rule (a namespace-scope entity with
# internal linkage cannot be reached from an exported template) that clang does
# not diagnose. MSVC stays limited to the quick test until the full suite is
# validated against it.
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    message(STATUS "Boost.Math module: ${CMAKE_CXX_COMPILER_ID} consumer tests limited to module_quick_test")
    set(BOOST_MATH_MODULE_TESTS)
endif()

# octonion_test_simple additionally uses Boost.MPL textually
set(BOOST_MATH_MODULE_EXTRA_LIBS_octonion_test_simple Boost::mpl)

add_executable(module_quick_test "${CMAKE_CURRENT_SOURCE_DIR}/quick_test.cpp")
target_link_libraries(module_quick_test PRIVATE boost_math_module Boost::core)
# The superproject pins an older policy version (CMP0155), so C++ sources are not
# scanned for module imports by default. Force scanning so `import boost.math;`
# is resolved to the module built above.
set_target_properties(module_quick_test PROPERTIES CXX_SCAN_FOR_MODULES ON CXX_MODULE_STD ON)
add_test(NAME module_quick_test COMMAND module_quick_test)
add_dependencies(tests module_quick_test)

foreach(test ${BOOST_MATH_MODULE_TESTS})
    add_executable(module_${test} "${CMAKE_CURRENT_SOURCE_DIR}/../test/${test}.cpp")
    target_link_libraries(module_${test}
            PRIVATE
            boost_math_module
            Boost::core
            ${BOOST_MATH_MODULE_EXTRA_LIBS_${test}}
    )
    set_target_properties(module_${test} PROPERTIES CXX_SCAN_FOR_MODULES ON CXX_MODULE_STD ON)
    add_test(NAME module_${test} COMMAND module_${test})
    add_dependencies(tests module_${test})
endforeach()
