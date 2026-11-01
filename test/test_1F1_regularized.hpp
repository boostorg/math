// Copyright John Maddock 2006.
// Copyright Paul A. Bristow 2007, 2009
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_MATH_OVERFLOW_ERROR_POLICY ignore_error

#include <boost/math/concepts/real_concept.hpp>
#include <boost/math/special_functions/math_fwd.hpp>
#include <boost/math/tools/big_constant.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/relative_difference.hpp>
#include "math_unit_test.hpp"

#include <boost/math/special_functions/hypergeometric_1F1.hpp>
#include <boost/math/quadrature/exp_sinh.hpp>

#ifdef _MSC_VER
#pragma warning(disable:4127)
#endif

// The largest relative error, in epsilon, allowed for each type and group of test data:
std::uintmax_t max_error(const std::string& type_name, const std::string& test_name);

// Checks every row {a, b, z, f(a, b, z)} of the data, and prints each row that fails:
template <class Real, class T, class F>
void test_1F1_data(const T& data, const char* type_name, const char* function_name, const char* test_name, F f)
{
   std::cout << "Testing " << test_name << " with type " << type_name << "\n";
   const std::uintmax_t bound = max_error(type_name, test_name);
   for (unsigned i = 0; i < data.size(); ++i)
   {
      const Real a = Real(data[i][0]), b = Real(data[i][1]), z = Real(data[i][2]), expected = Real(data[i][3]);
      Real computed = 0, err = 0;
      try
      {
         computed = f(a, b, z);
         err = boost::math::relative_difference(computed, expected) / boost::math::tools::epsilon<Real>();
      }
      catch (const std::exception& e)
      {
         std::cerr << e.what() << "\n";
         err = boost::math::tools::max_value<Real>();
      }
      if (!(err <= bound))
      {
         std::cerr << std::setprecision(std::numeric_limits<Real>::max_digits10) << "\033[0;31mError in " << test_name << " at row " << i << ":\033[0m\n  "
                   << function_name << "(" << a << ", " << b << ", " << z << ")\n  Expected: " << expected << "\n  Computed: " << computed
                   << std::setprecision(4) << "\n  Relative error: " << err << " epsilon, which exceeds " << bound << "\n";
         ++boost::math::test::detail::global_error_count();
      }
   }
}

template <class Real, class T>
void do_test_1F1(const T& data, const char* type_name, const char* test_name)
{
   test_1F1_data<Real>(data, type_name, "hypergeometric_1F1_regularized", test_name, [](Real a, Real b, Real z) { return boost::math::hypergeometric_1F1_regularized(a, b, z); });
}

#ifndef SC_
#define SC_(x) BOOST_MATH_BIG_CONSTANT(T, 1000000, x)
#endif

template <class T>
void test_spots1(T, const char* type_name)
{
#ifdef TEST_MP
#include "hypergeometric_1f1_large_regularized_mp.ipp"
#else
#include "hypergeometric_1f1_large_regularized.ipp"
#endif

   do_test_1F1<T>(hypergeometric_1f1_large_regularized, type_name, "Large random values - regularized");
}

template <class T>
void test_spots(T z, const char* type_name)
{
   test_spots1(z, type_name);

   //
   // Special cases for coverage:
   //
   // Huge b: the result underflows to zero long before the rescaling would finish:
   CHECK_EQUAL(boost::math::hypergeometric_1F1_regularized(T(0.5), T(1e30), T(0.5)), T(0));
   BOOST_MATH_IF_CONSTEXPR(std::numeric_limits<T>::has_infinity)
   {
      CHECK_EQUAL(boost::math::hypergeometric_1F1_regularized(600.25, 102.75, 11512.0), std::numeric_limits<double>::infinity());
      BOOST_MATH_IF_CONSTEXPR(std::numeric_limits<long double>::max_exponent > 16000)
      {
         CHECK_EQUAL(boost::math::hypergeometric_1F1_regularized(600.25, 102.75, 9985.0), std::numeric_limits<double>::infinity());
      }
      else BOOST_MATH_IF_CONSTEXPR(std::numeric_limits<T>::max_exponent <= std::numeric_limits<double>::max_exponent)
      {
         CHECK_EQUAL(boost::math::hypergeometric_1F1_regularized(600.25, 102.75, 514.0), std::numeric_limits<double>::infinity());
      }
   }
}

