//  (C) Copyright John Maddock 2005.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/test/test_tools.hpp>
#include <boost/test/included/test_exec_monitor.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/mpl/if.hpp>
#include <boost/static_assert.hpp>
#include <boost/math/complex.hpp>

#include <iostream>
#include <iomanip>
#include <cmath>

#ifdef BOOST_NO_STDC_NAMESPACE
namespace std{ using ::sqrt; using ::tan; using ::tanh; }
#endif

#ifndef VERBOSE
#undef BOOST_MESSAGE
#define BOOST_MESSAGE(x)
#endif

//
// check_complex:
// Verifies that expected value "a" and found value "b" have a relative error
// less than "max_error" epsilons.  Note that relative error is calculated for
// the complex number as a whole; this means that the error in the real or 
// imaginary parts alone can be much higher than max_error when the real and 
// imaginary parts are of very different magnitudes.  This is important, because
// the Hull et al analysis of the acos and asin algorithms requires that very small
// real/imaginary components can be safely ignored if they are negligible compared
// to the other component.
//
template <class T>
void check_complex(const std::complex<T>& a, const std::complex<T>& b, int max_error, int line_num = 0)
{
   //
   // a is the expected value, b is what was actually found,
   // compute | (a-b)/b | and compare with max_error which is the 
   // multiple of E to permit:
   //
   static const std::complex<T> zero(0);
   if(a == zero)
   {
      if(b != zero)
      {
         if(boost::math::fabs(b) > std::numeric_limits<T>::epsilon())
         {
            BOOST_MESSAGE("Error from complex_test.cpp, line " << line_num << ", with type: " << typeid(T).name());
            BOOST_ERROR("Expected {0,0} but got: " << b);
         }
         else
         {
            BOOST_MESSAGE("Error from complex_test.cpp, line " << line_num << ", with type: " << typeid(T).name());
            BOOST_MESSAGE("Expected {0,0} but got: " << b);
         }
      }
      return;
   }
   else if(b == zero)
   {
      if(boost::math::fabs(a) > std::numeric_limits<T>::epsilon())
      {
         BOOST_MESSAGE("Error from complex_test.cpp, line " << line_num << ", with type: " << typeid(T).name());
         BOOST_ERROR("Found {0,0} but expected: " << a);
         return;
      }
      else
      {
         BOOST_MESSAGE("Error from complex_test.cpp, line " << line_num << ", with type: " << typeid(T).name());
         BOOST_MESSAGE("Found {0,0} but expected: " << a);
      }
   }

   T rel = boost::math::fabs((b-a)/b) / std::numeric_limits<T>::epsilon();
   if( rel > max_error)
   {
      BOOST_MESSAGE("Error from complex_test.cpp, line " << line_num << ", with type: " << typeid(T).name());
      BOOST_ERROR("Error in result exceeded permitted limit of " << max_error << " (actual relative error was " << rel << "e).  Found " << b << " expected " << a);
   }
}

//
// test_inverse_trig:
// This is nothing more than a sanity check, computes trig(atrig(z)) 
// and compare the result to z.  Note that:
//
// atrig(trig(z)) != z
//
// for certain z because the inverse trig functions are multi-valued, this 
// essentially rules this out as a testing method.  On the other hand:
//
// trig(atrig(z))
//
// can vary compare to z by an arbitrarily large amount.  For one thing we 
// have no control over the implementation of the trig functions, for another
// even if both functions were accurate to 1ulp (as accurate as transcendental
// number can get, thanks to the "table makers dilemma"), the errors can still
// be arbitrarily large - often the inverse trig functions will map a very large
// part of the complex domain into a small output domain, so you can never get
// back exactly where you started from.  Consequently these tests are no more than
// sanity checks (just verifies that signs are correct and so on).
//
template <class T>
void test_inverse_trig(T)
{
   static const T interval = static_cast<T>(2.0L/128.0L);

   T x, y;

   for(x = -1; x <= 1; x += interval)
   {
      for(y = -1; y <= 1; y += interval)
      {
         // acos:
         std::complex<T> val(x, y);
         std::complex<T> result = std::cos(boost::math::acos(val));
         BOOST_MESSAGE("Testing cos(acos(z)), z = [" << x << "," << y << "]");
         check_complex(val, result, 50, __LINE__);
         // asin:
         result = std::sin(boost::math::asin(val));
         BOOST_MESSAGE("Testing sin(asin(z)), z = [" << x << "," << y << "]");
         check_complex(val, result, 5, __LINE__);
      }
   }

   static const T interval2 = static_cast<T>(3.0L/256.0L);
   for(x = -3; x <= 3; x += interval2)
   {
      for(y = -3; y <= 3; y += interval2)
      {
         // asinh:
         std::complex<T> val(x, y);
         std::complex<T> result = std::sinh(boost::math::asinh(val));
         BOOST_MESSAGE("Testing sinh(asinh(z)), z = [" << x << "," << y << "]");
         check_complex(val, result, 5, __LINE__);
         // acosh:
         if(!((y == 0) && (x <= 1))) // can't test along the branch cut
         {
            result = std::cosh(boost::math::acosh(val));
            BOOST_MESSAGE("Testing cosh(acosh(z)), z = [" << x << "," << y << "]");
            check_complex(val, result, 60, __LINE__);
         }
         //
         // There is a problem in testing atan and atanh:
         // The inverse functions map a large input range to a much
         // smaller output range, so at the extremes too rather different
         // inputs may map to the same output value once rounded to N places.
         // Consequently tan(atan(z)) can suffer from arbitrarily large errors
         // even if individually they each have a small error bound.  On the other
         // hand we can't test atan(tan(z)) either because atan is multi-valued, so
         // round-tripping in this direction isn't always possible.
         // The following heuristic is designed to make the best of a bad job,
         // using atan(tan(z)) where possible and tan(atan(z)) when it's not.
         //
         static const int tanh_error = 20;
         if((0 != x) && (0 != y) && ((std::fabs(y) < 1) || (std::fabs(x) < 1)))
         {
            // atanh:
            val = boost::math::atanh(val);
            result = boost::math::atanh(std::tanh(val));
            BOOST_MESSAGE("Testing atanh(tanh(z)), z = [" << val.real()< "," << val.imag() << "]");
            check_complex(val, result, tanh_error, __LINE__);
            // atan:
            if(!((x == 0) && (std::fabs(y) == 1))) // we can't test infinities here
            {
               val = boost::math::atan(std::complex<T>(x, y));
               result = boost::math::atan(std::tan(val));
               BOOST_MESSAGE("Testing atan(tan(z)), z = [" << val.real() << "," << val.imag() << "]");
               check_complex(val, result, tanh_error, __LINE__);
            }
         }
         else
         {
            // atanh:
            result = std::tanh(boost::math::atanh(val));
            BOOST_MESSAGE("Testing tanh(atanh(z)), z = [" << x << "," << y << "]");
            check_complex(val, result, tanh_error, __LINE__);
            // atan:
            if(!((x == 0) && (std::fabs(y) == 1))) // we can't test infinities here
            {
               result = std::tan(boost::math::atan(val));
               BOOST_MESSAGE("Testing tan(atan(z)), z = [" << x << "," << y << "]");
               check_complex(val, result, tanh_error, __LINE__);
            }
         }
      }
   }
}

//
// check_spots:
// Various spot values, mostly the C99 special cases (infinites and NAN's).
// TODO: add spot checks for the Wolfram spot values.
//
template <class T>
void check_spots(const T&)
{
   typedef std::complex<T> ct;
   ct result;
   T eps = std::numeric_limits<T>::epsilon();
   static const T zero = 0;
   static const T mzero = -zero;
   static const T one = 1;
   static const T pi = static_cast<T>(3.141592653589793238462643383279502884197L);
   static const T half_pi = static_cast<T>(1.57079632679489661923132169163975144L);
   static const T quarter_pi = static_cast<T>(0.78539816339744830961566084581987572L);
   static const T three_quarter_pi = static_cast<T>(2.35619449019234492884698253745962716L);
   //static const T log_two = static_cast<T>(0.69314718055994530941723212145817657L);
   T infinity = std::numeric_limits<T>::infinity();
   T nan = 0;
   if(std::numeric_limits<T>::has_quiet_NaN)
      nan = std::numeric_limits<T>::quiet_NaN();
   bool test_nan = false;
   if(boost::math::detail::test_is_nan(nan))
      test_nan = true;

   //
   // C99 spot tests for acos:
   //
   result = boost::math::acos(ct(zero));
   check_complex(ct(half_pi), result, 2, __LINE__);
   
   result = boost::math::acos(ct(mzero));
   check_complex(ct(half_pi), result, 2, __LINE__);
   
   result = boost::math::acos(ct(zero, mzero));
   check_complex(ct(half_pi), result, 2, __LINE__);
   
   result = boost::math::acos(ct(mzero, mzero));
   check_complex(ct(half_pi), result, 2, __LINE__);
   
   if(test_nan)
   {
      result = boost::math::acos(ct(zero,nan));
      BOOST_CHECK_CLOSE(result.real(), half_pi, eps*200);
      BOOST_CHECK(boost::math::detail::test_is_nan(result.imag()));
   
      result = boost::math::acos(ct(mzero,nan));
      BOOST_CHECK_CLOSE(result.real(), half_pi, eps*200);
      BOOST_CHECK(boost::math::detail::test_is_nan(result.imag()));
   }
   result = boost::math::acos(ct(zero, infinity));
   BOOST_CHECK_CLOSE(result.real(), half_pi, eps*200);
   BOOST_CHECK(result.imag() == -infinity);

   result = boost::math::acos(ct(zero, -infinity));
   BOOST_CHECK_CLOSE(result.real(), half_pi, eps*200);
   BOOST_CHECK(result.imag() == infinity);

   if(test_nan)
   {
      result = boost::math::acos(ct(one, nan));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.real()));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.imag()));
   }
   result = boost::math::acos(ct(-infinity, one));
   BOOST_CHECK_CLOSE(result.real(), pi, eps*200);
   BOOST_CHECK(result.imag() == -infinity);

   result = boost::math::acos(ct(infinity, one));
   BOOST_CHECK(result.real() == 0);
   BOOST_CHECK(result.imag() == -infinity);

   result = boost::math::acos(ct(-infinity, -one));
   BOOST_CHECK_CLOSE(result.real(), pi, eps*200);
   BOOST_CHECK(result.imag() == infinity);

   result = boost::math::acos(ct(infinity, -one));
   BOOST_CHECK(result.real() == 0);
   BOOST_CHECK(result.imag() == infinity);

   result = boost::math::acos(ct(-infinity, infinity));
   BOOST_CHECK_CLOSE(result.real(), three_quarter_pi, eps*200);
   BOOST_CHECK(result.imag() == -infinity);

   result = boost::math::acos(ct(infinity, infinity));
   BOOST_CHECK_CLOSE(result.real(), quarter_pi, eps*200);
   BOOST_CHECK(result.imag() == -infinity);

   result = boost::math::acos(ct(-infinity, -infinity));
   BOOST_CHECK_CLOSE(result.real(), three_quarter_pi, eps*200);
   BOOST_CHECK(result.imag() == infinity);

   result = boost::math::acos(ct(infinity, -infinity));
   BOOST_CHECK_CLOSE(result.real(), quarter_pi, eps*200);
   BOOST_CHECK(result.imag() == infinity);
   if(test_nan)
   {
      result = boost::math::acos(ct(infinity, nan));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.real()));
      BOOST_CHECK(std::fabs(result.imag()) == infinity);

      result = boost::math::acos(ct(-infinity, nan));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.real()));
      BOOST_CHECK(std::fabs(result.imag()) == infinity);

      result = boost::math::acos(ct(nan, zero));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.real()));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.imag()));

      result = boost::math::acos(ct(nan, -zero));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.real()));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.imag()));

      result = boost::math::acos(ct(nan, one));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.real()));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.imag()));

      result = boost::math::acos(ct(nan, -one));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.real()));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.imag()));

      result = boost::math::acos(ct(nan, nan));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.real()));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.imag()));

      result = boost::math::acos(ct(nan, infinity));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.real()));
      BOOST_CHECK(result.imag() == -infinity);
      
      result = boost::math::acos(ct(nan, -infinity));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.real()));
      BOOST_CHECK(result.imag() == infinity);
   }

   //
   // C99 spot tests for acosh:
   //
   result = boost::math::acosh(ct(zero, zero));
   BOOST_CHECK(result.real() == 0);
   BOOST_CHECK_CLOSE(result.imag(), half_pi, eps*200);

   result = boost::math::acosh(ct(zero, mzero));
   BOOST_CHECK(result.real() == 0);
   BOOST_CHECK_CLOSE(result.imag(), half_pi, eps*200);

   result = boost::math::acosh(ct(mzero, zero));
   BOOST_CHECK(result.real() == 0);
   BOOST_CHECK_CLOSE(result.imag(), half_pi, eps*200);
   
   result = boost::math::acosh(ct(mzero, mzero));
   BOOST_CHECK(result.real() == 0);
   BOOST_CHECK_CLOSE(result.imag(), half_pi, eps*200);
   
   result = boost::math::acosh(ct(one, infinity));
   BOOST_CHECK(result.real() == infinity);
   BOOST_CHECK_CLOSE(result.imag(), half_pi, eps*200);

   result = boost::math::acosh(ct(one, -infinity));
   BOOST_CHECK(result.real() == infinity);
   BOOST_CHECK_CLOSE(result.imag(), -half_pi, eps*200);

   if(test_nan)
   {
      result = boost::math::acosh(ct(one, nan));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.real()));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.imag()));
   }
   result = boost::math::acosh(ct(-infinity, one));
   BOOST_CHECK(result.real() == infinity);
   BOOST_CHECK_CLOSE(result.imag(), pi, eps*200);
   
   result = boost::math::acosh(ct(infinity, one));
   BOOST_CHECK(result.real() == infinity);
   BOOST_CHECK(result.imag() == 0);
   
   result = boost::math::acosh(ct(-infinity, -one));
   BOOST_CHECK(result.real() == infinity);
   BOOST_CHECK_CLOSE(result.imag(), -pi, eps*200);
   
   result = boost::math::acosh(ct(infinity, -one));
   BOOST_CHECK(result.real() == infinity);
   BOOST_CHECK(result.imag() == 0);
   
   result = boost::math::acosh(ct(-infinity, infinity));
   BOOST_CHECK(result.real() == infinity);
   BOOST_CHECK_CLOSE(result.imag(), three_quarter_pi, eps*200);
   
   result = boost::math::acosh(ct(infinity, infinity));
   BOOST_CHECK(result.real() == infinity);
   BOOST_CHECK_CLOSE(result.imag(), quarter_pi, eps*200);
   
   result = boost::math::acosh(ct(-infinity, -infinity));
   BOOST_CHECK(result.real() == infinity);
   BOOST_CHECK_CLOSE(result.imag(), -three_quarter_pi, eps*200);
   
   result = boost::math::acosh(ct(infinity, -infinity));
   BOOST_CHECK(result.real() == infinity);
   BOOST_CHECK_CLOSE(result.imag(), -quarter_pi, eps*200);
   
   if(test_nan)
   {
      result = boost::math::acosh(ct(infinity, nan));
      BOOST_CHECK(result.real() == infinity);
      BOOST_CHECK(boost::math::detail::test_is_nan(result.imag()));
      
      result = boost::math::acosh(ct(-infinity, nan));
      BOOST_CHECK(result.real() == infinity);
      BOOST_CHECK(boost::math::detail::test_is_nan(result.imag()));
      
      result = boost::math::acosh(ct(nan, one));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.real()));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.imag()));
      
      result = boost::math::acosh(ct(nan, infinity));
      BOOST_CHECK(result.real() == infinity);
      BOOST_CHECK(boost::math::detail::test_is_nan(result.imag()));
      
      result = boost::math::acosh(ct(nan, -one));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.real()));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.imag()));
      
      result = boost::math::acosh(ct(nan, -infinity));
      BOOST_CHECK(result.real() == infinity);
      BOOST_CHECK(boost::math::detail::test_is_nan(result.imag()));
      
      result = boost::math::acosh(ct(nan, nan));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.real()));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.imag()));
   }
   //
   // C99 spot checks for asinh:
   //
   result = boost::math::asinh(ct(zero, zero));
   BOOST_CHECK(result.real() == 0);
   BOOST_CHECK(result.imag() == 0);

   result = boost::math::asinh(ct(mzero, zero));
   BOOST_CHECK(result.real() == 0);
   BOOST_CHECK(result.imag() == 0);

   result = boost::math::asinh(ct(zero, mzero));
   BOOST_CHECK(result.real() == 0);
   BOOST_CHECK(result.imag() == 0);

   result = boost::math::asinh(ct(mzero, mzero));
   BOOST_CHECK(result.real() == 0);
   BOOST_CHECK(result.imag() == 0);

   result = boost::math::asinh(ct(one, infinity));
   BOOST_CHECK(result.real() == infinity);
   BOOST_CHECK_CLOSE(result.imag(), half_pi, eps*200);
   
   result = boost::math::asinh(ct(one, -infinity));
   BOOST_CHECK(result.real() == infinity);
   BOOST_CHECK_CLOSE(result.imag(), -half_pi, eps*200);
   
   result = boost::math::asinh(ct(-one, -infinity));
   BOOST_CHECK(result.real() == -infinity);
   BOOST_CHECK_CLOSE(result.imag(), -half_pi, eps*200);
   
   result = boost::math::asinh(ct(-one, infinity));
   BOOST_CHECK(result.real() == -infinity);
   BOOST_CHECK_CLOSE(result.imag(), half_pi, eps*200);

   if(test_nan)
   {
      result = boost::math::asinh(ct(one, nan));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.real()));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.imag()));

      result = boost::math::asinh(ct(-one, nan));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.real()));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.imag()));

      result = boost::math::asinh(ct(zero, nan));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.real()));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.imag()));
   }

   result = boost::math::asinh(ct(infinity, one));
   BOOST_CHECK(result.real() == infinity);
   BOOST_CHECK(result.imag() == 0);
   
   result = boost::math::asinh(ct(infinity, -one));
   BOOST_CHECK(result.real() == infinity);
   BOOST_CHECK(result.imag() == 0);
   
   result = boost::math::asinh(ct(-infinity, -one));
   BOOST_CHECK(result.real() == -infinity);
   BOOST_CHECK(result.imag() == 0);
   
   result = boost::math::asinh(ct(-infinity, one));
   BOOST_CHECK(result.real() == -infinity);
   BOOST_CHECK(result.imag() == 0);
   
   result = boost::math::asinh(ct(infinity, infinity));
   BOOST_CHECK(result.real() == infinity);
   BOOST_CHECK_CLOSE(result.imag(), quarter_pi, eps*200);
   
   result = boost::math::asinh(ct(infinity, -infinity));
   BOOST_CHECK(result.real() == infinity);
   BOOST_CHECK_CLOSE(result.imag(), -quarter_pi, eps*200);
   
   result = boost::math::asinh(ct(-infinity, -infinity));
   BOOST_CHECK(result.real() == -infinity);
   BOOST_CHECK_CLOSE(result.imag(), -quarter_pi, eps*200);
   
   result = boost::math::asinh(ct(-infinity, infinity));
   BOOST_CHECK(result.real() == -infinity);
   BOOST_CHECK_CLOSE(result.imag(), quarter_pi, eps*200);

   if(test_nan)
   {
      result = boost::math::asinh(ct(infinity, nan));
      BOOST_CHECK(result.real() == infinity);
      BOOST_CHECK(boost::math::detail::test_is_nan(result.imag()));

      result = boost::math::asinh(ct(-infinity, nan));
      BOOST_CHECK(result.real() == -infinity);
      BOOST_CHECK(boost::math::detail::test_is_nan(result.imag()));

      result = boost::math::asinh(ct(nan, zero));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.real()));
      BOOST_CHECK(result.imag() == 0);

      result = boost::math::asinh(ct(nan,  mzero));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.real()));
      BOOST_CHECK(result.imag() == 0);

      result = boost::math::asinh(ct(nan, one));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.real()));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.imag()));

      result = boost::math::asinh(ct(nan,  -one));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.real()));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.imag()));

      result = boost::math::asinh(ct(nan,  nan));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.real()));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.imag()));

      result = boost::math::asinh(ct(nan, infinity));
      BOOST_CHECK(std::fabs(result.real()) == infinity);
      BOOST_CHECK(boost::math::detail::test_is_nan(result.imag()));

      result = boost::math::asinh(ct(nan,  -infinity));
      BOOST_CHECK(std::fabs(result.real()) == infinity);
      BOOST_CHECK(boost::math::detail::test_is_nan(result.imag()));
   }
   
   //
   // C99 special cases for atanh:
   //
   result = boost::math::atanh(ct(zero, zero));
   BOOST_CHECK(result.real() == zero);
   BOOST_CHECK(result.imag() == zero);

   result = boost::math::atanh(ct(mzero, zero));
   BOOST_CHECK(result.real() == zero);
   BOOST_CHECK(result.imag() == zero);

   result = boost::math::atanh(ct(zero, mzero));
   BOOST_CHECK(result.real() == zero);
   BOOST_CHECK(result.imag() == zero);

   result = boost::math::atanh(ct(mzero, mzero));
   BOOST_CHECK(result.real() == zero);
   BOOST_CHECK(result.imag() == zero);

   if(test_nan)
   {
      result = boost::math::atanh(ct(zero, nan));
      BOOST_CHECK(result.real() == zero);
      BOOST_CHECK(boost::math::detail::test_is_nan(result.imag()));

      result = boost::math::atanh(ct(-zero, nan));
      BOOST_CHECK(result.real() == zero);
      BOOST_CHECK(boost::math::detail::test_is_nan(result.imag()));
   }

   result = boost::math::atanh(ct(one, zero));
   BOOST_CHECK_EQUAL(result.real(), infinity);
   BOOST_CHECK_EQUAL(result.imag(), zero);

   result = boost::math::atanh(ct(-one, zero));
   BOOST_CHECK_EQUAL(result.real(), -infinity);
   BOOST_CHECK_EQUAL(result.imag(), zero);

   result = boost::math::atanh(ct(-one, -zero));
   BOOST_CHECK_EQUAL(result.real(), -infinity);
   BOOST_CHECK_EQUAL(result.imag(), zero);

   result = boost::math::atanh(ct(one, -zero));
   BOOST_CHECK_EQUAL(result.real(), infinity);
   BOOST_CHECK_EQUAL(result.imag(), zero);

   result = boost::math::atanh(ct(pi, infinity));
   BOOST_CHECK_EQUAL(result.real(), zero);
   BOOST_CHECK_CLOSE(result.imag(), half_pi, eps*200);

   result = boost::math::atanh(ct(pi, -infinity));
   BOOST_CHECK_EQUAL(result.real(), zero);
   BOOST_CHECK_CLOSE(result.imag(), -half_pi, eps*200);

   result = boost::math::atanh(ct(-pi, -infinity));
   BOOST_CHECK_EQUAL(result.real(), zero);
   BOOST_CHECK_CLOSE(result.imag(), -half_pi, eps*200);

   result = boost::math::atanh(ct(-pi, infinity));
   BOOST_CHECK_EQUAL(result.real(), zero);
   BOOST_CHECK_CLOSE(result.imag(), half_pi, eps*200);

   if(test_nan)
   {
      result = boost::math::atanh(ct(pi, nan));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.real()));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.imag()));

      result = boost::math::atanh(ct(-pi, nan));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.real()));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.imag()));
   }

   result = boost::math::atanh(ct(infinity, pi));
   BOOST_CHECK(result.real() == zero);
   BOOST_CHECK_CLOSE(result.imag(), half_pi, eps*200);

   result = boost::math::atanh(ct(infinity, -pi));
   BOOST_CHECK(result.real() == zero);
   BOOST_CHECK_CLOSE(result.imag(), -half_pi, eps*200);

   result = boost::math::atanh(ct(-infinity, -pi));
   BOOST_CHECK(result.real() == zero);
   BOOST_CHECK_CLOSE(result.imag(), -half_pi, eps*200);

   result = boost::math::atanh(ct(-infinity, pi));
   BOOST_CHECK(result.real() == zero);
   BOOST_CHECK_CLOSE(result.imag(), half_pi, eps*200);

   result = boost::math::atanh(ct(infinity, infinity));
   BOOST_CHECK(result.real() == zero);
   BOOST_CHECK_CLOSE(result.imag(), half_pi, eps*200);

   result = boost::math::atanh(ct(infinity, -infinity));
   BOOST_CHECK(result.real() == zero);
   BOOST_CHECK_CLOSE(result.imag(), -half_pi, eps*200);

   result = boost::math::atanh(ct(-infinity, -infinity));
   BOOST_CHECK(result.real() == zero);
   BOOST_CHECK_CLOSE(result.imag(), -half_pi, eps*200);

   result = boost::math::atanh(ct(-infinity, infinity));
   BOOST_CHECK(result.real() == zero);
   BOOST_CHECK_CLOSE(result.imag(), half_pi, eps*200);

   if(test_nan)
   {
      result = boost::math::atanh(ct(infinity, nan));
      BOOST_CHECK(result.real() == 0);
      BOOST_CHECK(boost::math::detail::test_is_nan(result.imag()));

      result = boost::math::atanh(ct(-infinity, nan));
      BOOST_CHECK(result.real() == 0);
      BOOST_CHECK(boost::math::detail::test_is_nan(result.imag()));

      result = boost::math::atanh(ct(nan, pi));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.real()));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.imag()));

      result = boost::math::atanh(ct(nan, -pi));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.real()));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.imag()));

      result = boost::math::atanh(ct(nan, infinity));
      BOOST_CHECK(result.real() == 0);
      BOOST_CHECK_CLOSE(result.imag(), half_pi, eps*200);

      result = boost::math::atanh(ct(nan, -infinity));
      BOOST_CHECK(result.real() == 0);
      BOOST_CHECK_CLOSE(result.imag(), -half_pi, eps*200);

      result = boost::math::atanh(ct(nan, nan));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.real()));
      BOOST_CHECK(boost::math::detail::test_is_nan(result.imag()));

   }
}

//
// test_boundaries:
// This is an accuracy test, sets the real and imaginary components
// of the input argument to various "boundary conditions" that exist
// inside the implementation.  Then computes the result at double precision
// and again at float precision.  The double precision result will be
// computed using the "regular" code, where as the float precision versions
// will calculate the result using the "exceptional value" handlers, so
// we end up comparing the values calculated by two different methods.
//
const float boundaries[] = {
   0,
   1,
   2,
   (std::numeric_limits<float>::max)(),
   (std::numeric_limits<float>::min)(),
   std::numeric_limits<float>::epsilon(),
   std::sqrt((std::numeric_limits<float>::max)()) / 8,
   static_cast<float>(4) * std::sqrt((std::numeric_limits<float>::min)()),
   0.6417F,
   1.5F,
   std::sqrt((std::numeric_limits<float>::max)()) / 2,
   std::sqrt((std::numeric_limits<float>::min)()),
   1.0F / 0.3F,
};

void do_test_boundaries(float x, float y)
{
   std::complex<float> r1 = boost::math::asin(std::complex<float>(x, y));
   std::complex<double> dr = boost::math::asin(std::complex<double>(x, y));
   std::complex<float> r2(static_cast<float>(dr.real()), static_cast<float>(dr.imag()));
   check_complex(r2, r1, 5, __LINE__);
   r1 = boost::math::acos(std::complex<float>(x, y));
   dr = boost::math::acos(std::complex<double>(x, y));
   r2 = std::complex<float>(std::complex<double>(dr.real(), dr.imag()));
   check_complex(r2, r1, 5, __LINE__);
   r1 = boost::math::atanh(std::complex<float>(x, y));
   dr = boost::math::atanh(std::complex<double>(x, y));
   r2 = std::complex<float>(std::complex<double>(dr.real(), dr.imag()));
   check_complex(r2, r1, 5, __LINE__);
}

void test_boundaries(float x, float y)
{
   do_test_boundaries(x, y);
   do_test_boundaries(-x, y); 
   do_test_boundaries(-x, -y);
   do_test_boundaries(x, -y);
}

void test_boundaries(float x)
{
   for(unsigned i = 0; i < sizeof(boundaries)/sizeof(float); ++i)
   {
      test_boundaries(x, boundaries[i]);
      test_boundaries(x, boundaries[i] + std::numeric_limits<float>::epsilon()*boundaries[i]);
      test_boundaries(x, boundaries[i] - std::numeric_limits<float>::epsilon()*boundaries[i]);
   }
}

void test_boundaries()
{
   for(unsigned i = 0; i < sizeof(boundaries)/sizeof(float); ++i)
   {
      test_boundaries(boundaries[i]);
      test_boundaries(boundaries[i] + std::numeric_limits<float>::epsilon()*boundaries[i]);
      test_boundaries(boundaries[i] - std::numeric_limits<float>::epsilon()*boundaries[i]);//here
   }
}


int test_main(int, char*[])
{
   std::cout << "Running complex trig sanity checks for type float." << std::endl;
   test_inverse_trig(float(0));
   std::cout << "Running complex trig sanity checks for type double." << std::endl;
   test_inverse_trig(double(0));
   //test_inverse_trig((long double)(0));

   std::cout << "Running complex trig spot checks for type float." << std::endl;
   check_spots(float(0));
   std::cout << "Running complex trig spot checks for type double." << std::endl;
   check_spots(double(0));
   std::cout << "Running complex trig spot checks for type long double." << std::endl;
   check_spots((long double)(0));
   
   std::cout << "Running complex trig boundary and accuracy tests." << std::endl;
   test_boundaries();
   return 0;
}

