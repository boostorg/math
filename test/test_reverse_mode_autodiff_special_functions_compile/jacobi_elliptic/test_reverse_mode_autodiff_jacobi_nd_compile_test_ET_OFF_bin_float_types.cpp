
    //           Copyright Maksym Zhelyeznyakov 2025-2026
    // Distributed under the Boost Software License, Version 1.0.
    // (See accompanying file LICENSE_1_0.txt or copy at
    //           https://www.boost.org/LICENSE_1_0.txt)
    //
    // This file was generated automatically with math/tests/autogen_rvar_specfun_tests.sh
    // DO NOT EDIT MANUALLY
    
    #define BOOST_MATH_REVERSE_MODE_ET_OFF
    #include <test_autodiff_reverse.hpp>
    #include <boost/math/special_functions.hpp>
    #include <boost/math/special_functions/logaddexp.hpp>
    #include <boost/math/tools/workaround.hpp>
    #include <cmath>
    
    namespace test_jacobi_nd_compiles {
    
    using namespace rdiff;
    
    template <typename T>
    void test_jacobi_nd()
    {
        RandomSample<T> rng{ -10, 10
 };
        T               x = rng.next();
    
        rvar<T, 1> x_ad = x;
        auto y = boost::math::jacobi_nd(0.25,x_ad);
        auto y_expect = boost::math::jacobi_nd(0.25,x);
        CHECK_RELATIVE_ERROR(y.item(), y_expect, 1000*boost_close_tol<T>() / 100);
    }
    } // namespace test_jacobi_nd_compiles
    
    int main()
    {
        mp11::mp_for_each<mp11::mp_transform<mp11::mp_identity, bin_float_types>>([](auto t) { test_jacobi_nd_compiles::test_jacobi_nd<typename decltype(t)::type>(); });
        return boost::math::test::report_errors();
    }
    
