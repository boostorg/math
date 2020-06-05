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

    std::cout << "Testing x = 0\n";
    x = 0;
    a = to_simple_continued_fraction(x);
    CHECK_EQUAL(static_cast<int64_t>(a.size()), int64_t(1));
    CHECK_EQUAL(a[0], int64_t(0));

    std::cout << "Testing x = 1\n";
    x = 1;
    a = to_simple_continued_fraction(x);
    CHECK_EQUAL(static_cast<int64_t>(a.size()), int64_t(1));
    CHECK_EQUAL(a[0], int64_t(1));

    // need to test some negative numbers:
    x = -1;
    a = to_simple_continued_fraction(x);
    CHECK_EQUAL(static_cast<int64_t>(a.size()), int64_t(1));
    CHECK_EQUAL(a[0], int64_t(-1));

    // -1.5 = -2 + 1/2
    x = -3/Real(2);
    a = to_simple_continued_fraction(x);
    CHECK_EQUAL(static_cast<int64_t>(a.size()), int64_t(2));
    CHECK_EQUAL(a[0], int64_t(-2));
    CHECK_EQUAL(a[1], int64_t(2));

    // First integer is zero:
    x = 1/Real(2);
    a = to_simple_continued_fraction(x);
    CHECK_EQUAL(static_cast<int64_t>(a.size()), int64_t(2));
    CHECK_EQUAL(a[0], int64_t(0));
    CHECK_EQUAL(a[1], int64_t(2));


    // Less trivial:

    x = Real(649)/200;
    a = to_simple_continued_fraction(x);
    CHECK_EQUAL(static_cast<int64_t>(a.size()), int64_t(4));
    CHECK_EQUAL(a[0], int64_t(3));
    CHECK_EQUAL(a[1], int64_t(4));
    CHECK_EQUAL(a[2], int64_t(12));
    CHECK_EQUAL(a[3], int64_t(4));


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
