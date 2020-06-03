/*
 * Copyright Nick Thompson, 2020
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include "math_unit_test.hpp"
#include <boost/core/demangle.hpp>
#include <boost/math/tools/fraction.hpp>
#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
using boost::multiprecision::float128;
#endif

using boost::math::tools::to_simple_continued_fraction;

template<class Real>
void test_fraction()
{
    std::cout << "Testing conversion to continued fraction on type " << boost::core::demangle(typeid(Real).name()) << "\n";
    Real x = Real(415)/Real(93);
    std::cout << "Decimal of x = " << x << "\n"; 
    // [4; 2, 6, 7]:
    std::vector<int64_t> a = to_simple_continued_fraction(x);
    //CHECK_EQUAL(static_cast<int64_t>(a.size()), int64_t(4));
    CHECK_EQUAL(a[0], int64_t(4));
    CHECK_EQUAL(a[1], int64_t(2));
    CHECK_EQUAL(a[2], int64_t(6));
    CHECK_EQUAL(a[3], int64_t(7));

    x = 0;
    a = to_simple_continued_fraction(x);
    CHECK_EQUAL(static_cast<int64_t>(a.size()), int64_t(1));
    CHECK_EQUAL(a[0], int64_t(0));

    x = 1;
    a = to_simple_continued_fraction(x);
    CHECK_EQUAL(static_cast<int64_t>(a.size()), int64_t(1));
    CHECK_EQUAL(a[0], int64_t(1));

    // need a negative number to test:

}




int main()
{
    test_fraction<float>();
    test_fraction<double>();
    test_fraction<long double>();
    #ifdef BOOST_HAS_FLOAT128
    test_fraction<float128>();
    #endif
    return boost::math::test::report_errors();
}
