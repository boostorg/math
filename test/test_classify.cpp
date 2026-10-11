//  Copyright John Maddock 2006.
//  Copyright Paul A. Bristow 2007
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <math.h>
#include <boost/limits.hpp>
#ifndef BOOST_MATH_STANDALONE
#include <boost/math/concepts/real_concept.hpp>
#endif
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/special_functions/next.hpp>  // for has_denorm_now
#include "math_unit_test.hpp"
#include <iostream>
#include <iomanip>

#include "test_autodiff.hpp"

#ifdef _MSC_VER
#pragma warning(disable: 4127 4146) //  conditional expression is constant
#endif

const char* method_name(const boost::math::detail::native_tag&)
{
   return "Native";
}

const char* method_name(const boost::math::detail::generic_tag<true>&)
{
   return "Generic (with numeric limits)";
}

const char* method_name(const boost::math::detail::generic_tag<false>&)
{
   return "Generic (without numeric limits)";
}

const char* method_name(const boost::math::detail::ieee_tag&)
{
   return "IEEE std";
}

const char* method_name(const boost::math::detail::ieee_copy_all_bits_tag&)
{
   return "IEEE std, copy all bits";
}

const char* method_name(const boost::math::detail::ieee_copy_leading_bits_tag&)
{
   return "IEEE std, copy leading bits";
}

template <class T>
void test_classify(T t, const char* type)
{
   std::cout << "Testing type " << type << std::endl;

   typedef typename boost::math::detail::fp_traits<T>::type traits;
   typedef typename traits::method method;

   std::cout << "Evaluation method = " << method_name(method()) << std::endl;   

   t = 2;
   T u = 2;
   CHECK_EQUAL((::boost::math::fpclassify)(t), (int)FP_NORMAL);
   CHECK_EQUAL((::boost::math::fpclassify)(-t), (int)FP_NORMAL);
   CHECK_EQUAL((::boost::math::isfinite)(t), true);
   CHECK_EQUAL((::boost::math::isfinite)(-t), true);
   CHECK_EQUAL((::boost::math::isinf)(t), false);
   CHECK_EQUAL((::boost::math::isinf)(-t), false);
   CHECK_EQUAL((::boost::math::isnan)(t), false);
   CHECK_EQUAL((::boost::math::isnan)(-t), false);
   CHECK_EQUAL((::boost::math::isnormal)(t), true);
   CHECK_EQUAL((::boost::math::isnormal)(-t), true);
   if(std::numeric_limits<T>::is_specialized)
   {
      t = (std::numeric_limits<T>::max)();
      CHECK_EQUAL((::boost::math::fpclassify)(t), (int)FP_NORMAL);
      CHECK_EQUAL((::boost::math::fpclassify)(-t), (int)FP_NORMAL);
      CHECK_EQUAL((::boost::math::isfinite)(t), true);
      CHECK_EQUAL((::boost::math::isfinite)(-t), true);
      CHECK_EQUAL((::boost::math::isinf)(t), false);
      CHECK_EQUAL((::boost::math::isinf)(-t), false);
      CHECK_EQUAL((::boost::math::isnan)(t), false);
      CHECK_EQUAL((::boost::math::isnan)(-t), false);
      CHECK_EQUAL((::boost::math::isnormal)(t), true);
      CHECK_EQUAL((::boost::math::isnormal)(-t), true);
      t = (std::numeric_limits<T>::min)();
      if(t != 0)
      {
         CHECK_EQUAL((::boost::math::fpclassify)(t), (int)FP_NORMAL);
         CHECK_EQUAL((::boost::math::isfinite)(t), true);
         CHECK_EQUAL((::boost::math::isinf)(t), false);
         CHECK_EQUAL((::boost::math::isnan)(t), false);
         CHECK_EQUAL((::boost::math::isnormal)(t), true);
         if(!std::numeric_limits<T>::is_integer)
         {
            CHECK_EQUAL((::boost::math::fpclassify)(-t), (int)FP_NORMAL);
            CHECK_EQUAL((::boost::math::isfinite)(-t), true);
            CHECK_EQUAL((::boost::math::isinf)(-t), false);
            CHECK_EQUAL((::boost::math::isnormal)(-t), true);
            CHECK_EQUAL((::boost::math::isnan)(-t), false);
         }
      }
   }
   if(boost::math::detail::has_denorm_now<T>())
   {
      t = (std::numeric_limits<T>::min)();
      t /= 2;
      if(t != 0)
      {
         CHECK_EQUAL((::boost::math::fpclassify)(t), (int)FP_SUBNORMAL);
         CHECK_EQUAL((::boost::math::fpclassify)(-t), (int)FP_SUBNORMAL);
         CHECK_EQUAL((::boost::math::isfinite)(t), true);
         CHECK_EQUAL((::boost::math::isfinite)(-t), true);
         CHECK_EQUAL((::boost::math::isinf)(t), false);
         CHECK_EQUAL((::boost::math::isinf)(-t), false);
         CHECK_EQUAL((::boost::math::isnan)(t), false);
         CHECK_EQUAL((::boost::math::isnan)(-t), false);
         CHECK_EQUAL((::boost::math::isnormal)(t), false);
         CHECK_EQUAL((::boost::math::isnormal)(-t), false);
      }
      t = std::numeric_limits<T>::denorm_min();
      if((t != 0) && (t < (std::numeric_limits<T>::min)()))
      {
         CHECK_EQUAL((::boost::math::fpclassify)(t), (int)FP_SUBNORMAL);
         CHECK_EQUAL((::boost::math::fpclassify)(-t), (int)FP_SUBNORMAL);
         CHECK_EQUAL((::boost::math::isfinite)(t), true);
         CHECK_EQUAL((::boost::math::isfinite)(-t), true);
         CHECK_EQUAL((::boost::math::isinf)(t), false);
         CHECK_EQUAL((::boost::math::isinf)(-t), false);
         CHECK_EQUAL((::boost::math::isnan)(t), false);
         CHECK_EQUAL((::boost::math::isnan)(-t), false);
         CHECK_EQUAL((::boost::math::isnormal)(t), false);
         CHECK_EQUAL((::boost::math::isnormal)(-t), false);
      }
   }
   else
   {
      std::cout << "Denormalised forms not tested" << std::endl;
   }
   t = 0;
   CHECK_EQUAL((::boost::math::fpclassify)(t), (int)FP_ZERO);
   CHECK_EQUAL((::boost::math::fpclassify)(-t), (int)FP_ZERO);
   CHECK_EQUAL((::boost::math::isfinite)(t), true);
   CHECK_EQUAL((::boost::math::isfinite)(-t), true);
   CHECK_EQUAL((::boost::math::isinf)(t), false);
   CHECK_EQUAL((::boost::math::isinf)(-t), false);
   CHECK_EQUAL((::boost::math::isnan)(t), false);
   CHECK_EQUAL((::boost::math::isnan)(-t), false);
   CHECK_EQUAL((::boost::math::isnormal)(t), false);
   CHECK_EQUAL((::boost::math::isnormal)(-t), false);
   t /= -u; // create minus zero if it exists
   CHECK_EQUAL((::boost::math::fpclassify)(t), (int)FP_ZERO);
   CHECK_EQUAL((::boost::math::fpclassify)(-t), (int)FP_ZERO);
   CHECK_EQUAL((::boost::math::isfinite)(t), true);
   CHECK_EQUAL((::boost::math::isfinite)(-t), true);
   CHECK_EQUAL((::boost::math::isinf)(t), false);
   CHECK_EQUAL((::boost::math::isinf)(-t), false);
   CHECK_EQUAL((::boost::math::isnan)(t), false);
   CHECK_EQUAL((::boost::math::isnan)(-t), false);
   CHECK_EQUAL((::boost::math::isnormal)(t), false);
   CHECK_EQUAL((::boost::math::isnormal)(-t), false);
   // infinity:
   if(std::numeric_limits<T>::has_infinity) 
   {
      // At least one std::numeric_limits<T>::infinity)() returns zero 
      // (Compaq true64 cxx), hence the check.
      t = (std::numeric_limits<T>::infinity)();
      CHECK_EQUAL((::boost::math::fpclassify)(t), (int)FP_INFINITE);
      CHECK_EQUAL((::boost::math::fpclassify)(-t), (int)FP_INFINITE);
      CHECK_EQUAL((::boost::math::isfinite)(t), false);
      CHECK_EQUAL((::boost::math::isfinite)(-t), false);
      CHECK_EQUAL((::boost::math::isinf)(t), true);
      CHECK_EQUAL((::boost::math::isinf)(-t), true);
      CHECK_EQUAL((::boost::math::isnan)(t), false);
      CHECK_EQUAL((::boost::math::isnan)(-t), false);
      CHECK_EQUAL((::boost::math::isnormal)(t), false);
      CHECK_EQUAL((::boost::math::isnormal)(-t), false);
#if !defined(BOOST_BORLANDC) && !(defined(__DECCXX) && !defined(_IEEE_FP))
      // divide by zero on Borland triggers a C++ exception :-(
      // divide by zero on Compaq CXX triggers a C style signal :-(
      t = 2;
      u = 0;
      t /= u;
      CHECK_EQUAL((::boost::math::fpclassify)(t), (int)FP_INFINITE);
      CHECK_EQUAL((::boost::math::fpclassify)(-t), (int)FP_INFINITE);
      CHECK_EQUAL((::boost::math::isfinite)(t), false);
      CHECK_EQUAL((::boost::math::isfinite)(-t), false);
      CHECK_EQUAL((::boost::math::isinf)(t), true);
      CHECK_EQUAL((::boost::math::isinf)(-t), true);
      CHECK_EQUAL((::boost::math::isnan)(t), false);
      CHECK_EQUAL((::boost::math::isnan)(-t), false);
      CHECK_EQUAL((::boost::math::isnormal)(t), false);
      CHECK_EQUAL((::boost::math::isnormal)(-t), false);
      t = -2;
      t /= u;
      CHECK_EQUAL((::boost::math::fpclassify)(t), (int)FP_INFINITE);
      CHECK_EQUAL((::boost::math::fpclassify)(-t), (int)FP_INFINITE);
      CHECK_EQUAL((::boost::math::isfinite)(t), false);
      CHECK_EQUAL((::boost::math::isfinite)(-t), false);
      CHECK_EQUAL((::boost::math::isinf)(t), true);
      CHECK_EQUAL((::boost::math::isinf)(-t), true);
      CHECK_EQUAL((::boost::math::isnan)(t), false);
      CHECK_EQUAL((::boost::math::isnan)(-t), false);
      CHECK_EQUAL((::boost::math::isnormal)(t), false);
      CHECK_EQUAL((::boost::math::isnormal)(-t), false);
#else
      std::cout << "Infinities from divide by zero not tested" << std::endl;
#endif
   }
   else
   {
      std::cout << "Infinity not tested" << std::endl;
   }
#ifndef BOOST_BORLANDC
   // NaN's:
   // Note that Borland throws an exception if we even try to obtain a Nan
   // by calling std::numeric_limits<T>::quiet_NaN() !!!!!!!
   if(std::numeric_limits<T>::has_quiet_NaN)
   {
      t = std::numeric_limits<T>::quiet_NaN();
      CHECK_EQUAL((::boost::math::fpclassify)(t), (int)FP_NAN);
      CHECK_EQUAL((::boost::math::fpclassify)(-t), (int)FP_NAN);
      CHECK_EQUAL((::boost::math::isfinite)(t), false);
      CHECK_EQUAL((::boost::math::isfinite)(-t), false);
      CHECK_EQUAL((::boost::math::isinf)(t), false);
      CHECK_EQUAL((::boost::math::isinf)(-t), false);
      CHECK_EQUAL((::boost::math::isnan)(t), true);
      CHECK_EQUAL((::boost::math::isnan)(-t), true);
      CHECK_EQUAL((::boost::math::isnormal)(t), false);
      CHECK_EQUAL((::boost::math::isnormal)(-t), false);
   }
   else
   {
      std::cout << "Quiet NaN's not tested" << std::endl;
   }
   if(std::numeric_limits<T>::has_signaling_NaN)
   {
      t = std::numeric_limits<T>::signaling_NaN();
      CHECK_EQUAL((::boost::math::fpclassify)(t), (int)FP_NAN);
      CHECK_EQUAL((::boost::math::fpclassify)(-t), (int)FP_NAN);
      CHECK_EQUAL((::boost::math::isfinite)(t), false);
      CHECK_EQUAL((::boost::math::isfinite)(-t), false);
      CHECK_EQUAL((::boost::math::isinf)(t), false);
      CHECK_EQUAL((::boost::math::isinf)(-t), false);
      CHECK_EQUAL((::boost::math::isnan)(t), true);
      CHECK_EQUAL((::boost::math::isnan)(-t), true);
      CHECK_EQUAL((::boost::math::isnormal)(t), false);
      CHECK_EQUAL((::boost::math::isnormal)(-t), false);
   }
   else
   {
      std::cout << "Signaling NaN's not tested" << std::endl;
   }
#endif
}



void test_main()
{
   BOOST_MATH_CONTROL_FP;
   // start by printing some information:
#ifdef isnan
   std::cout << "Platform has isnan macro." << std::endl;
#endif
#ifdef fpclassify
   std::cout << "Platform has fpclassify macro." << std::endl;
#endif
#ifdef BOOST_HAS_FPCLASSIFY
   std::cout << "Platform has FP_NORMAL macro." << std::endl;
#endif
   std::cout << "FP_ZERO: " << (int)FP_ZERO << std::endl;
   std::cout << "FP_NORMAL: " << (int)FP_NORMAL << std::endl;
   std::cout << "FP_INFINITE: " << (int)FP_INFINITE << std::endl;
   std::cout << "FP_NAN: " << (int)FP_NAN << std::endl;
   std::cout << "FP_SUBNORMAL: " << (int)FP_SUBNORMAL << std::endl;

   // then run the tests:
   test_classify(float(0), "float");
   test_classify(double(0), "double");
   // long double support for fpclassify is considered "core" so we always test it
   // even when long double support is turned off via BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
   test_classify((long double)(0), "long double");
#ifndef BOOST_MATH_STANDALONE
   test_classify((boost::math::concepts::real_concept)(0), "real_concept");
#endif

   // We should test with integer types as well:
   test_classify(int(0), "int");
   test_classify(unsigned(0), "unsigned");
}

template <class T>
void fpclassify_autodiff() {
   test_classify(boost::math::differentiation::make_fvar<T, 1>(0), "autodiff float");
   test_classify(boost::math::differentiation::make_fvar<T, 2>(0), "autodiff float");
   test_classify(boost::math::differentiation::make_fvar<T, 3>(0), "autodiff float");
   test_classify(boost::math::differentiation::make_fvar<T, 7>(0), "autodiff float");
   test_classify(boost::math::differentiation::make_fvar<T, 12>(0), "autodiff float");
}


/*
Autorun "i:\Boost-sandbox\math_toolkit\libs\math\test\MSVC80\debug\test_classify.exe"
Running 1 test case...
FP_ZERO: 0
FP_NORMAL: 1
FP_INFINITE: 2
FP_NAN: 3
FP_SUBNORMAL: 4
Testing type float
Testing type double
Testing type long double
Testing type real_concept
Denormalised forms not tested
Infinity not tested
Quiet NaN's not tested
Signaling NaN's not tested
Test suite "Test Program" passed with:
  79 assertions out of 79 passed
  1 test case out of 1 passed
  Test case "test_main_caller( argc, argv )" passed with:
    79 assertions out of 79 passed

*/

int main()
{
   test_main();
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { fpclassify_autodiff<typename decltype(t)::type>(); });
   return boost::math::test::report_errors();
}
