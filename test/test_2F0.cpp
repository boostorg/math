//  (C) Copyright John Maddock 2006.
//  Copyright Paul A. Bristow 2007, 2009
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_MATH_OVERFLOW_ERROR_POLICY ignore_error

#include "math_unit_test.hpp"
#include <boost/math/concepts/real_concept.hpp>
#include <boost/math/special_functions/hypergeometric_2F0.hpp>
#include <boost/math/special_functions/relative_difference.hpp>
#include <boost/math/tools/stats.hpp>
#if !defined(TEST) || (TEST == 3)
#include <boost/multiprecision/cpp_bin_float.hpp>
#endif
#if !defined(TEST) || (TEST == 4)
#include <boost/multiprecision/cpp_dec_float.hpp>
#endif
#include <array>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <type_traits>

//
// The tables are arrays of strings, converted to the type under test at run time:
// some 16000 multiprecision constants in one function cost the compiler gigabytes.
//
namespace tables {
using T = const char*;
#define SC_(x) #x
#include "hypergeometric_2F0.ipp"
#include "hypergeometric_2F0_large_z.ipp"
#include "hypergeometric_2F0_integer_a2.ipp"
#include "hypergeometric_2F0_half.ipp"
#undef SC_
}

template <class T>
T from_string(const char* s) { return T(s); }
template <>
float from_string<float>(const char* s) { return std::strtof(s, nullptr); }
template <>
double from_string<double>(const char* s) { return std::strtod(s, nullptr); }
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
template <>
long double from_string<long double>(const char* s) { return std::strtold(s, nullptr); }
template <>
boost::math::concepts::real_concept from_string<boost::math::concepts::real_concept>(const char* s) { return std::strtold(s, nullptr); }
#endif

// Peak and root mean square relative error, in units of epsilon.
struct error_bound { double peak, rms; };

template <class T, std::size_t N>
void test_table(const std::array<std::array<const char*, 4>, N>& table, const char* type_name, const char* test_name, error_bound bound)
{
   if (std::is_same<T, float>::value)
   {
      bound = { 1, 1 };  // float is evaluated in double.
   }

   boost::math::tools::stats<T> errors;
   std::size_t worst = 0;
   for (std::size_t i = 0; i < N; ++i)
   {
      const T expected = from_string<T>(table[i][3]);
      T computed = std::numeric_limits<T>::quiet_NaN();
      try
      {
         computed = boost::math::hypergeometric_2F0(from_string<T>(table[i][0]), from_string<T>(table[i][1]), from_string<T>(table[i][2]));
      }
      catch (const std::exception& e)
      {
         std::cerr << "  " << e.what() << "\n";
      }
      const T error = boost::math::relative_difference(computed, expected);
      if (!(error <= (errors.max)()))
      {
         worst = i;
      }
      // A non-finite result where a finite one is expected, or a gross error, fails by itself:
      if (!CHECK_TRUE(error <= 0.5f))
      {
         std::cerr << "  " << test_name << "<" << type_name << "> row " << i << ": found " << std::setprecision(std::numeric_limits<T>::max_digits10) << computed << ", expected " << table[i][3] << "\n";
      }
      errors.add(error);
   }
   const T eps = boost::math::tools::epsilon<T>();
   const T peak = (errors.max)() / eps;
   const T rms = errors.rms() / eps;
   std::cout << std::setprecision(4) << test_name << "<" << type_name << "> Max = " << peak << " RMS Mean=" << rms << ", worst at 2F0(" << table[worst][0] << ", " << table[worst][1] << "; " << table[worst][2] << ")" << std::endl;
   CHECK_LE(peak, T(bound.peak));
   CHECK_LE(rms, T(bound.rms));
}

template <class T>
void test_spots(const char* type_name)
{
   test_table<T>(tables::hypergeometric_2F0, type_name, "Random non-integer a2, |z| < 1", { 9000, 3000 });
   test_table<T>(tables::hypergeometric_2F0_large_z, type_name, "Random non-integer a2, |z| > 1", { 9000, 3000 });
   test_table<T>(tables::hypergeometric_2F0_integer_a2, type_name, "Integer a2, |z| > 1", { 300, 100 });
   test_table<T>(tables::hypergeometric_2F0_half, type_name, "a1 = a2 + 0.5", { 2000, 200 });

   T a1 = 1;
   T a2 = 0.1f;
   T z = -1;
   if (std::numeric_limits<T>::has_infinity)
   {
      CHECK_THROW(boost::math::hypergeometric_2F0(-std::numeric_limits<T>::infinity(), a2, T(0.5), boost::math::policies::make_policy(boost::math::policies::rounding_error<boost::math::policies::ignore_error>())), std::domain_error);
      CHECK_TRUE((boost::math::isinf)(boost::math::hypergeometric_2F0(a1, a2, z)));
      z = 1;
      CHECK_TRUE((boost::math::isinf)(boost::math::hypergeometric_2F0(a1, a2, z)));
      a2 = -2.4f;
      CHECK_TRUE((boost::math::isinf)(boost::math::hypergeometric_2F0(a1, a2, z)));
      a1 = -0.5f;
      CHECK_TRUE((boost::math::isinf)(boost::math::hypergeometric_2F0(a1, a2, z)));
   }
   //
   // Coverage:
   //
   CHECK_EQUAL(boost::math::hypergeometric_2F0(T(0), T(20), T(2)), T(1));
   CHECK_EQUAL(boost::math::hypergeometric_2F0(T(20), T(0), T(2)), T(1));
   CHECK_EQUAL(boost::math::hypergeometric_2F0(T(20), T(10), T(0)), T(1));
   //
   // Negative integer a1 with a2 > 0 and z > 0, where the series cancels: exact values, from the series
   // in rational arithmetic.  The last case may raise an evaluation_error rather than lose digits:
   //
   T tolerance = boost::math::tools::epsilon<T>() * 4;
   CHECK_RELATIVE_ERROR(from_string<T>("0.0009058711400604743835444245801821637120072"), boost::math::hypergeometric_2F0(T(-20), T(7.5), T(0.0703125)), tolerance);
   CHECK_RELATIVE_ERROR(from_string<T>("0.0009058711400604743835444245801821637120072"), boost::math::hypergeometric_2F0(T(7.5), T(-20), T(0.0703125)), tolerance);
   CHECK_RELATIVE_ERROR(from_string<T>("0.3052122739096701050752227260666316776838"), boost::math::hypergeometric_2F0(T(-50), T(1.25), T(0.03125)), tolerance);
   CHECK_RELATIVE_ERROR(from_string<T>("0.5866164156536446434351584193455355526903"), boost::math::hypergeometric_2F0(T(-30), T(0.5), T(0.0625)), tolerance);
   CHECK_RELATIVE_ERROR(from_string<T>("0.2278524180520048379881596526813574062803"), boost::math::hypergeometric_2F0(T(-40), T(3), T(0.015625)), tolerance);
#ifndef BOOST_NO_EXCEPTIONS
   try
   {
      CHECK_RELATIVE_ERROR(from_string<T>("0.08262089934873659278920475540110746903509"), boost::math::hypergeometric_2F0(T(-100), T(1.75), T(0.03125)), tolerance);
   }
   catch (const boost::math::evaluation_error&) {}
#endif
}

int main()
{
   BOOST_MATH_CONTROL_FP;
#if !defined(TEST) || (TEST == 1)
   test_spots<float>("float");
   test_spots<double>("double");
#endif
#if !defined(TEST) || (TEST == 2)
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
   test_spots<long double>("long double");
#ifndef BOOST_MATH_NO_REAL_CONCEPT_TESTS
   test_spots<boost::math::concepts::real_concept>("real_concept");
#endif
#endif
#endif

#ifndef BOOST_MATH_NO_MP_TESTS
#if !defined(TEST) || (TEST == 3)
   #if defined(BOOST_MATH_RUN_MP_TESTS)
   test_spots<boost::multiprecision::cpp_bin_float_quad>("cpp_bin_float_quad");
   #endif
#endif
#if !defined(TEST) || (TEST == 4)
   #if defined(BOOST_MATH_RUN_MP_TESTS)
   using dec_40 = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<40>>;
   test_spots<dec_40>("dec_40");
   #endif
#endif
#endif
   return boost::math::test::report_errors();
}
