//  (C) Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <limits>
#include <boost/math/special_functions/sqrt.hpp>
#include <boost/core/lightweight_test.hpp>

template <typename Real>
void test_float_sqrt()
{
    using std::abs;
    
    constexpr Real tol = 2*std::numeric_limits<Real>::epsilon();
    
    constexpr Real test_val = boost::math::sqrt(Real(2));
    constexpr Real dummy = 1;
    static_assert(test_val > dummy, "Not constexpr");

    Real known_val = std::sqrt(Real(2));

    BOOST_TEST(abs(test_val - known_val) < tol);
}

template <typename Z>
void test_int_sqrt()
{
    using std::abs;

    constexpr double tol = 2*std::numeric_limits<double>::epsilon();

    constexpr double test_val = boost::math::sqrt(Z(2));
    constexpr double dummy = 1;
    static_assert(test_val > dummy, "Not constexpr");

    double known_val = std::sqrt(2.0);

    BOOST_TEST(abs(test_val - known_val) < tol);
}

int main()
{
    test_float_sqrt<float>();
    test_float_sqrt<double>();
    test_float_sqrt<long double>();

    test_int_sqrt<int>();
    test_int_sqrt<unsigned>();
    
    return boost::report_errors();
}
