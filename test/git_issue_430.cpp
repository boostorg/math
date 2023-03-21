// Copyright Matt Borland 2023
// Copyright John Maddock 2023
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0. (See accompanying file
// LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// See: https://godbolt.org/z/Ev4ManrsW

#include <boost/math/special_functions/round.hpp>
#include <boost/math/special_functions/next.hpp>
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <cstdint>
#include <cmath>

template <typename Real>
void test_llround_near_boundary()
{
   using std::ldexp;
   Real boundary = ldexp(static_cast<Real>(1), std::numeric_limits<long long>::digits);

   Real value;
   int i;

   for (value = boundary, i = 0; i < 100; value = boost::math::float_next(value), ++i)
   {
      BOOST_CHECK_THROW(boost::math::llround(value), boost::math::rounding_error);
   }
   for (value = boost::math::float_prior(boundary), i = 0; i < 1000; value = boost::math::float_prior(value), ++i)
   {
      BOOST_CHECK_EQUAL(static_cast<Real>(boost::math::llround(value)), boost::math::round(value));
   }
   for (value = boost::math::float_prior(-boundary), i = 0; i < 100; value = boost::math::float_prior(value), ++i)
   {
      BOOST_CHECK_THROW(boost::math::llround(value), boost::math::rounding_error);
   }
   for (value = -boundary, i = 0; i < 1000; value = boost::math::float_next(value), ++i)
   {
      BOOST_CHECK_EQUAL(static_cast<Real>(boost::math::llround(value)), boost::math::round(value));
   }
}

BOOST_AUTO_TEST_CASE( test_main )
{
    test_llround_near_boundary<float>();
    test_llround_near_boundary<double>();
    test_llround_near_boundary<long double>();
}
