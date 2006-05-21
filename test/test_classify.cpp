//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <math.h>
#include <boost/limits.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/test/included/test_exec_monitor.hpp>

template <class T>
void test_classify(T t, const char* type)
{
   std::cout << "Testing type " << type << std::endl;
   t = 2;
   T u = 2;
   BOOST_CHECK_EQUAL((::boost::math::fpclassify)(t), (int)FP_NORMAL);
   BOOST_CHECK_EQUAL((::boost::math::fpclassify)(-t), (int)FP_NORMAL);
   t = (std::numeric_limits<T>::max)();
   BOOST_CHECK_EQUAL((::boost::math::fpclassify)(t), (int)FP_NORMAL);
   BOOST_CHECK_EQUAL((::boost::math::fpclassify)(-t), (int)FP_NORMAL);
   t = (std::numeric_limits<T>::min)();
   BOOST_CHECK_EQUAL((::boost::math::fpclassify)(t), (int)FP_NORMAL);
   BOOST_CHECK_EQUAL((::boost::math::fpclassify)(-t), (int)FP_NORMAL);
   if(std::numeric_limits<T>::has_denorm)
   {
      t /= 2;
      if(t != 0)
      {
         BOOST_CHECK_EQUAL((::boost::math::fpclassify)(t), (int)FP_SUBNORMAL);
         BOOST_CHECK_EQUAL((::boost::math::fpclassify)(-t), (int)FP_SUBNORMAL);
      }
      t = std::numeric_limits<T>::denorm_min();
      if(t != 0)
      {
         BOOST_CHECK_EQUAL((::boost::math::fpclassify)(t), (int)FP_SUBNORMAL);
         BOOST_CHECK_EQUAL((::boost::math::fpclassify)(-t), (int)FP_SUBNORMAL);
      }
   }
   else
   {
      std::cout << "Denormalised forms not tested" << std::endl;
   }
   t = 0;
   BOOST_CHECK_EQUAL((::boost::math::fpclassify)(t), (int)FP_ZERO);
   BOOST_CHECK_EQUAL((::boost::math::fpclassify)(-t), (int)FP_ZERO);
   t /= -u; // create minus zero if it exists
   BOOST_CHECK_EQUAL((::boost::math::fpclassify)(t), (int)FP_ZERO);
   BOOST_CHECK_EQUAL((::boost::math::fpclassify)(-t), (int)FP_ZERO);
   // inifinity:
   if(std::numeric_limits<T>::has_infinity) 
   {
      // At least one std::numeric_limits<T>::infinity)() returns zero 
      // (Compaq true64 cxx), hence the check.
      t = (std::numeric_limits<T>::infinity)();
      BOOST_CHECK_EQUAL((::boost::math::fpclassify)(t), (int)FP_INFINITE);
      BOOST_CHECK_EQUAL((::boost::math::fpclassify)(-t), (int)FP_INFINITE);
   }
   else
   {
      std::cout << "Infinity not tested" << std::endl;
   }
#if !defined(__BORLANDC__) && !(defined(__DECCXX) && !defined(_IEEE_FP))
   // divide by zero on Borland triggers a C++ exception :-(
   // divide by zero on Compaq CXX triggers a C style signal :-(
   t = 2;
   u = 0;
   t /= u;
   BOOST_CHECK_EQUAL((::boost::math::fpclassify)(t), (int)FP_INFINITE);
   BOOST_CHECK_EQUAL((::boost::math::fpclassify)(-t), (int)FP_INFINITE);
   t = -2;
   t /= u;
   BOOST_CHECK_EQUAL((::boost::math::fpclassify)(t), (int)FP_INFINITE);
   BOOST_CHECK_EQUAL((::boost::math::fpclassify)(-t), (int)FP_INFINITE);
#else
   std::cout << "Infinities from divide by zero not tested" << std::endl;
#endif
#ifndef __BORLANDC__
   // NaN's:
   // Note that Borland throws an exception if we even try to obtain a Nan
   // by calling std::numeric_limits<T>::quiet_NaN() !!!!!!!
   if(std::numeric_limits<T>::has_quiet_NaN)
   {
      t = std::numeric_limits<T>::quiet_NaN();
      BOOST_CHECK_EQUAL((::boost::math::fpclassify)(t), (int)FP_NAN);
      BOOST_CHECK_EQUAL((::boost::math::fpclassify)(-t), (int)FP_NAN);
   }
   else
   {
      std::cout << "Quite NaN's not tested" << std::endl;
   }
   if(std::numeric_limits<T>::has_signaling_NaN)
   {
      t = std::numeric_limits<T>::signaling_NaN();
      BOOST_CHECK_EQUAL((::boost::math::fpclassify)(t), (int)FP_NAN);
      BOOST_CHECK_EQUAL((::boost::math::fpclassify)(-t), (int)FP_NAN);
   }
   else
   {
      std::cout << "Signaling NaN's not tested" << std::endl;
   }
#endif
}

int test_main(int, char* [] )
{
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
   test_classify((long double)(0), "long double");
   return 0;
}







