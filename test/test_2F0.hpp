// Copyright John Maddock 2006.
// Copyright Paul A. Bristow 2007, 2009
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_MATH_OVERFLOW_ERROR_POLICY ignore_error

#include <boost/math/concepts/real_concept.hpp>
#include <boost/math/special_functions/math_fwd.hpp>
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/math/tools/stats.hpp>
#include <boost/math/tools/test.hpp>
#include <boost/math/tools/big_constant.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/array.hpp>
#include "functor.hpp"

#include "handle_test_result.hpp"
#include "table_type.hpp"

#include <boost/math/special_functions/hypergeometric_2F0.hpp>

template <class Real, class T>
void do_test_2F0(const T& data, const char* type_name, const char* test_name)
{
   typedef Real                   value_type;

   typedef value_type(*pg)(value_type, value_type, value_type);
#if defined(BOOST_MATH_NO_DEDUCED_FUNCTION_POINTERS)
   pg funcp = boost::math::hypergeometric_0F1<value_type, value_type>;
#else
   pg funcp = boost::math::hypergeometric_2F0;
#endif

   boost::math::tools::test_result<value_type> result;

   std::cout << "Testing " << test_name << " with type " << type_name
      << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

   //
   // test hypergeometric_2F0 against data:
   //
   result = boost::math::tools::test_hetero<Real>(
      data,
      bind_func<Real>(funcp, 0, 1, 2),
      extract_result<Real>(3));
   handle_test_result(result, data[result.worst()], result.worst(), type_name, "hypergeometric_0F1", test_name);
   std::cout << std::endl;
}

#ifndef SC_
#define SC_(x) BOOST_MATH_BIG_CONSTANT(T, 1000000, x)
#endif

template <class T>
void test_spots1(T, const char* type_name)
{
#include "hypergeometric_2F0.ipp"

   do_test_2F0<T>(hypergeometric_2F0, type_name, "Random non-integer a2, |z| < 1");

#include "hypergeometric_2F0_large_z.ipp"

   do_test_2F0<T>(hypergeometric_2F0_large_z, type_name, "Random non-integer a2, |z| > 1");
}

template <class T>
void test_spots2(T, const char* type_name)
{
#include "hypergeometric_2F0_integer_a2.ipp"

   do_test_2F0<T>(hypergeometric_2F0_integer_a2, type_name, "Integer a2, |z| > 1");

#include "hypergeometric_2F0_half.ipp"

   do_test_2F0<T>(hypergeometric_2F0_half, type_name, "a1 = a2 + 0.5");
}

template <class T>
void test_spots(T z, const char* type_name)
{
   test_spots1(z, type_name);
   test_spots2(z, type_name);

   T a1 = 1;
   T a2 = 0.1f;
   z = -1;
   if (std::numeric_limits<T>::has_infinity)
   {
      BOOST_CHECK((boost::math::isinf)(boost::math::hypergeometric_2F0(a1, a2, z)));
      z = 1;
      BOOST_CHECK((boost::math::isinf)(boost::math::hypergeometric_2F0(a1, a2, z)));
      a2 = -2.4f;
      BOOST_CHECK((boost::math::isinf)(boost::math::hypergeometric_2F0(a1, a2, z)));
      a1 = -0.5f;
      BOOST_CHECK((boost::math::isinf)(boost::math::hypergeometric_2F0(a1, a2, z)));
   }
   //
   // Coverage:
   //
   BOOST_CHECK_EQUAL(boost::math::hypergeometric_2F0(T(0), T(20), T(2)), T(1));
   BOOST_CHECK_EQUAL(boost::math::hypergeometric_2F0(T(20), T(0), T(2)), T(1));
   BOOST_CHECK_EQUAL(boost::math::hypergeometric_2F0(T(20), T(10), T(0)), T(1));
   //
   // Negative integer a1 with a2 > 0 and z > 0, where the series cancels: exact values, from the series
   // in rational arithmetic.  The last case may raise an evaluation_error rather than lose digits:
   //
   T tolerance = boost::math::tools::epsilon<T>() * 4;
   BOOST_CHECK_CLOSE_FRACTION(boost::math::hypergeometric_2F0(T(-20), T(7.5), T(0.0703125)), T(0.0009058711400604743835444245801821637120072L), tolerance);
   BOOST_CHECK_CLOSE_FRACTION(boost::math::hypergeometric_2F0(T(7.5), T(-20), T(0.0703125)), T(0.0009058711400604743835444245801821637120072L), tolerance);
   BOOST_CHECK_CLOSE_FRACTION(boost::math::hypergeometric_2F0(T(-50), T(1.25), T(0.03125)), T(0.3052122739096701050752227260666316776838L), tolerance);
   BOOST_CHECK_CLOSE_FRACTION(boost::math::hypergeometric_2F0(T(-30), T(0.5), T(0.0625)), T(0.5866164156536446434351584193455355526903L), tolerance);
   BOOST_CHECK_CLOSE_FRACTION(boost::math::hypergeometric_2F0(T(-40), T(3), T(0.015625)), T(0.2278524180520048379881596526813574062803L), tolerance);
#ifndef BOOST_NO_EXCEPTIONS
   try
   {
      BOOST_CHECK_CLOSE_FRACTION(boost::math::hypergeometric_2F0(T(-100), T(1.75), T(0.03125)), T(0.08262089934873659278920475540110746903509L), tolerance);
   }
   catch (const boost::math::evaluation_error&) {}
#endif
}
