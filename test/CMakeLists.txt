# Copyright 2021 Matt Borland
# Distributed under the Boost Software License, Version 1.0.
# https://www.boost.org/LICENSE_1_0.txt

include(BoostTestJamfile OPTIONAL RESULT_VARIABLE HAVE_BOOST_TEST)

if(HAVE_BOOST_TEST)

    if (BOOST_MATH_ENABLE_CUDA)

        message(STATUS "Building boost.cuda_math with CUDA")

        enable_language(CUDA)
        set(CMAKE_CUDA_EXTENSIONS OFF)

        enable_testing()

        boost_test_jamfile(FILE cuda_jamfile LINK_LIBRARIES Boost::math Boost::assert Boost::concept_check Boost::config Boost::core Boost::integer Boost::lexical_cast Boost::multiprecision Boost::predef Boost::random Boost::static_assert Boost::throw_exception Boost::unit_test_framework )
    elseif (BOOST_MATH_ENABLE_SYCL)

        message(STATUS "Building boost.cuda_math with SYCL")

        set(CMAKE_CXX_COMPILER "icpx")
        set(CMAKE_C_COMPILER "icx")

        enable_testing()

        boost_test_jamfile(FILE sycl_jamfile LINK_LIBRARIES Boost::math Boost::assert Boost::concept_check Boost::config Boost::core Boost::integer Boost::lexical_cast Boost::multiprecision Boost::predef Boost::random Boost::static_assert Boost::throw_exception Boost::unit_test_framework sycl COMPILE_OPTIONS -fsycl )
    else()

        add_subdirectory(compile_test)

    endif()

endif()