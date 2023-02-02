// Copyright John Maddock, 2023
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <limits>
#include <boost/math/distributions/negative_binomial.hpp>
#include "math_unit_test.hpp"

int main() 
{
   using namespace boost::math;

   negative_binomial n(5.0, 0.5);

   for (double k = 0; k < 15; ++k)
   {
      auto c = cdf(n, k);
      auto q = quantile(n, c);
      CHECK_EQUAL(q, k);
   }

   return boost::math::test::report_errors();
}
