
    //           Copyright Maksym Zhelyeznyakov 2025-2026
    // Distributed under the Boost Software License, Version 1.0.
    // (See accompanying file LICENSE_1_0.txt or copy at
    //           https://www.boost.org/LICENSE_1_0.txt)
    //
    // This file was generated automatically with math/tests/autogen_rvar_specfun_tests.sh
    // DO NOT EDIT MANUALLY
    
    #define BOOST_MATH_REVERSE_MODE_ET_ON
    #include <test_autodiff_reverse.hpp>
    #include <boost/math/special_functions.hpp>
    #include <boost/math/special_functions/logaddexp.hpp>
    #include <boost/math/tools/workaround.hpp>
    #include <cmath>
    
    BOOST_AUTO_TEST_SUITE(test_heuman_lambda_compiles)
    
    using namespace rdiff;
    
    BOOST_AUTO_TEST_CASE_TEMPLATE(test_heuman_lambda, T, bin_float_types)
    {
        RandomSample<T> rng{ 0.1, 1.0
 };
        T               x = rng.next();
    
        rvar<T, 1> x_ad = x;
        auto y = boost::math::heuman_lambda(0,x_ad);
        auto y_expect = boost::math::heuman_lambda(0,x);
        BOOST_CHECK_CLOSE(y.item(), y_expect, 1000*boost_close_tol<T>());
    }
    BOOST_AUTO_TEST_SUITE_END()
    
