//  (C) Copyright Jeremy Murphy 2015.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/config.hpp>
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#include <utility>

#include <boost/array.hpp>
#include <boost/math/tools/polynomial.hpp>
#include <boost/math/common_factor_rt.hpp>

using namespace boost::math;
using namespace boost::math::tools;
using namespace std;

template <typename T>
struct answer
{
    answer(std::pair< polynomial<T>, polynomial<T> > const &x) :
    quotient(x.first), remainder(x.second) {}
    
    polynomial<T> quotient;
    polynomial<T> remainder;
};

typedef polynomial<double> PR;

boost::array<double, 4> const d3a = {{10, -6, -4, 3}};
boost::array<double, 4> const d3b = {{-7, 5, 6, 1}};
boost::array<double, 4> const d3c = {{10.0/3.0, -2.0, -4.0/3.0, 1.0}};
boost::array<double, 2> const d1a = {{-2, 1}};
boost::array<double, 3> const d2a = {{-2, 2, 3}};
boost::array<double, 3> const d2b = {{-7, 5, 6}};
boost::array<double, 3> const d2c = {{31, -21, -22}};
boost::array<double, 1> const d0a = {{6}};
boost::array<double, 1> const d0b = {{3}};

PR const a(d3a.begin(), d3a.end());
PR const b(d1a.begin(), d1a.end());
PR const q(d2a.begin(), d2a.end());
PR const r(d0a.begin(), d0a.end());
PR const c(d3b.begin(), d3b.end());
PR const d(d2b.begin(), d2b.end());
PR const e(d2c.begin(), d2c.end());
PR const f(d0b.begin(), d0b.end());
PR const g(d3c.begin(), d3c.end());
PR const zero = zero_element(std::multiplies< polynomial<double> >());

BOOST_AUTO_TEST_CASE( test_division )
{
    BOOST_CHECK_THROW(quotient_remainder(a, zero), std::domain_error);
    
    answer<double> result = quotient_remainder(a, b);
    BOOST_CHECK_EQUAL(result.quotient, q);
    BOOST_CHECK_EQUAL(result.remainder, r);
    BOOST_CHECK_EQUAL(a, q * b + r); // Sanity check.
    
    result = quotient_remainder(a, c);
    BOOST_CHECK_EQUAL(result.quotient, f);
    BOOST_CHECK_EQUAL(result.remainder, e);
    BOOST_CHECK_EQUAL(a, f * c + e); // Sanity check.
    
    result = quotient_remainder(a, f);
    BOOST_CHECK_EQUAL(result.quotient, g);
    BOOST_CHECK_EQUAL(result.remainder, zero);
    BOOST_CHECK_EQUAL(a, g * f + zero); // Sanity check.
}

BOOST_AUTO_TEST_CASE( test_gcd )
{
    boost::array<double, 9> const d8 = {{0, 0, 0, 0, 0, 0, 0, 0, 1}};
    boost::array<double, 7> const d6 = {{0, 0, 0, 0, 0, 0, 1}};
    PR aa(d8.begin(), d8.end());
    PR bb(d6.begin(), d6.end());
    PR foo = gcd(aa, bb);
    BOOST_CHECK_EQUAL(foo, bb);
}
