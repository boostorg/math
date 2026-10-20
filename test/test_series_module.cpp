//  Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Module build check for the series-summation tools in
// <boost/math/tools/series.hpp>. The component is consumed through
// import boost.math and the public sum_series, checked_sum_series and
// kahan_sum_series entry points are exercised. There is no runtime non-module
// test for this component (only a compile-only test), so the expected values
// are mathematically-exact closed forms: the geometric series sum over n >= 0
// of r^n equals 1 / (1 - r) for |r| < 1, and the exponential series sum over
// n >= 0 of 1/n! equals e.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/tools/series.hpp>
#include <boost/math/constants/constants.hpp>
#else
import boost.math;
#endif

#include "math_unit_test.hpp"

// Yields the geometric sequence r^n for n = 0, 1, 2, ...; the partial sums
// converge to 1 / (1 - r) when |r| < 1.
template <class Real>
struct geometric_term
{
   using result_type = Real;

   result_type ratio;
   result_type current;

   geometric_term(result_type r, result_type first = static_cast<result_type>(1))
       : ratio{r}, current{first} {}

   result_type operator()()
   {
      result_type value {current};
      current *= ratio;
      return value;
   }
};

// Yields the terms 1/n! for n = 0, 1, 2, ...; the partial sums converge to e.
template <class Real>
struct exp_term
{
   using result_type = Real;

   result_type term;
   int n;

   exp_term() : term{static_cast<result_type>(1)}, n{0} {}

   result_type operator()()
   {
      result_type value {term};
      ++n;
      term /= static_cast<result_type>(n);
      return value;
   }
};

template <class Real>
void test_spots(const char* type_name)
{
   std::cout << "Testing series summation for type " << type_name << std::endl;

   const int bits {boost::math::numeric_limits<Real>::digits};

   // Geometric series with r = 1/2: sum over n >= 0 of (1/2)^n is 2.
   {
      geometric_term<Real> func {static_cast<Real>(0.5)};
      Real Q {boost::math::tools::sum_series(func, bits)};
      CHECK_ULP_CLOSE(static_cast<Real>(2), Q, 4);
   }

   // Geometric series with r = 1/4: sum over n >= 0 of (1/4)^n is 4/3.
   {
      geometric_term<Real> func {static_cast<Real>(0.25)};
      Real Q {boost::math::tools::sum_series(func, bits)};
      CHECK_ULP_CLOSE(static_cast<Real>(4) / 3, Q, 4);
   }

   // Alternating geometric series with r = -1/2: sum over n >= 0 of (-1/2)^n
   // is 2/3.
   {
      geometric_term<Real> func {static_cast<Real>(-0.5)};
      Real Q {boost::math::tools::sum_series(func, bits)};
      CHECK_ULP_CLOSE(static_cast<Real>(2) / 3, Q, 8);
   }

   // max_terms is both the cap on the number of iterations and, on return, the
   // count of terms actually evaluated. A cap of 3 stops this slowly converging
   // series early, so exactly 3 terms are summed and max_terms comes back as 3.
   {
      geometric_term<Real> func {static_cast<Real>(0.5)};
      const Real factor {std::ldexp(static_cast<Real>(1), 1 - bits)};
      boost::math::uintmax_t max_terms {3};
      boost::math::tools::sum_series(func, factor, max_terms);
      CHECK_EQUAL(static_cast<boost::math::uintmax_t>(3), max_terms);
   }

   // checked_sum_series returns the same value as sum_series and additionally
   // accumulates the L1 norm sum of |term|. For the alternating series r = -1/2
   // the value is 2/3 while the norm is the non-alternating geometric sum, 2.
   {
      geometric_term<Real> func {static_cast<Real>(-0.5)};
      const Real factor {std::ldexp(static_cast<Real>(1), 1 - bits)};
      boost::math::uintmax_t max_terms {1000};
      Real norm {0};
      Real Q {boost::math::tools::checked_sum_series(func, factor, max_terms, static_cast<Real>(0), norm)};
      CHECK_ULP_CLOSE(static_cast<Real>(2) / 3, Q, 8);
      CHECK_ULP_CLOSE(static_cast<Real>(2), norm, 8);
   }

   // Kahan-compensated summation of the exponential series sum over n >= 0 of
   // 1/n!, which is e.
   {
      exp_term<Real> func {};
      Real Q {boost::math::tools::kahan_sum_series(func, bits)};
      CHECK_ULP_CLOSE(boost::math::constants::e<Real>(), Q, 16);
   }
}

int main()
{
   test_spots<double>("double");
   test_spots<float>("float");

   return boost::math::test::report_errors();
}
