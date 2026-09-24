//  Copyright John Maddock 2025.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_TEST_MAIN
#include <boost/math/distributions.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <iostream>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

typedef boost::math::policies::policy<
   boost::math::policies::promote_double<false>
> SciPyPolicy;

using mp_t = boost::multiprecision::cpp_bin_float_quad;

BOOST_AUTO_TEST_CASE(test_main)
{
   double v = 980.0;
   double l = 38.0;
   double x = 1.5;

   // With policy
   double y_policy = boost::math::cdf(
      boost::math::non_central_t_distribution<double, SciPyPolicy>(v, l), x);
   double reference = 1.1824454111413493e-291;

   std::cout << std::setprecision(17);
   std::cout << "With SciPyPolicy:    cdf(noncentral_t; x=" << x << ", v=" << v << ", l=" << l << ") = " << y_policy << std::endl;

   double tolerance = 1e-8;  // precision is limited due to cancellation
   BOOST_CHECK_CLOSE_FRACTION(y_policy, reference, tolerance);
}

