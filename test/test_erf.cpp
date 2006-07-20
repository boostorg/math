//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/concepts/real_concept.hpp>
#include <boost/test/included/test_exec_monitor.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/tools/stats.hpp>
#include <boost/math/tools/test.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/array.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include "test_erf_hooks.hpp"


template <class T>
void print_test_result(const boost::math::tools::test_result<T>& result, 
                       T worst, const char* name, const char* test)
{
   using namespace std;
   T eps = pow(T(2), 1-boost::math::tools::digits<T>());
   std::cout << setprecision(4);
   std::cout << test << "(" << name << ") Max = " << (result.stat.max)()/eps
      << " RMS Mean=" << result.stat.rms()/eps << " worst case at point: " << worst << std::endl;
}

template <class T>
void do_test_erf(const T& data, const char* type_name, const char* test_name)
{
   typedef typename T::value_type row_type;
   typedef typename row_type::value_type value_type;

   typedef value_type (*pg)(value_type);
   pg funcp = boost::math::erf;

   boost::math::tools::test_result<value_type> result;

   std::cout << "Testing " << test_name << " with type " << type_name
      << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

   //
   // test erf against data:
   //
   result = boost::math::tools::test(
      data, 
      boost::lambda::bind(funcp, 
         boost::lambda::ret<value_type>(boost::lambda::_1[0])), 
      boost::lambda::ret<value_type>(boost::lambda::_1[1]));
   print_test_result(result, data[result.worst_case][0], type_name, "boost::math::erf");
#ifdef TEST_OTHER
   if(::boost::is_floating_point<value_type>::value){
      funcp = other::erf;
      result = boost::math::tools::test(
         data, 
         boost::lambda::bind(funcp, 
            boost::lambda::ret<value_type>(boost::lambda::_1[0])), 
         boost::lambda::ret<value_type>(boost::lambda::_1[1]));
      print_test_result(result, data[result.worst_case][0], type_name, "other::erf");
   }
#endif
   //
   // test erfc against data:
   //
   funcp = boost::math::erfc;
   result = boost::math::tools::test(
      data, 
      boost::lambda::bind(funcp, 
         boost::lambda::ret<value_type>(boost::lambda::_1[0])), 
      boost::lambda::ret<value_type>(boost::lambda::_1[2]));
   print_test_result(result, data[result.worst_case][0], type_name, "boost::math::erfc");
#ifdef TEST_OTHER
   if(::boost::is_floating_point<value_type>::value){
      funcp = other::erfc;
      result = boost::math::tools::test(
         data, 
         boost::lambda::bind(funcp, 
            boost::lambda::ret<value_type>(boost::lambda::_1[0])), 
         boost::lambda::ret<value_type>(boost::lambda::_1[2]));
      print_test_result(result, data[result.worst_case][0], type_name, "other::erfc");
   }
#endif
   std::cout << std::endl;
}

template <class T>
void do_test_erf_inv(const T& data, const char* type_name, const char* test_name)
{
   typedef typename T::value_type row_type;
   typedef typename row_type::value_type value_type;

   typedef value_type (*pg)(value_type);
   pg funcp = boost::math::erf;

   boost::math::tools::test_result<value_type> result;
   std::cout << "Testing " << test_name << " with type " << type_name
      << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
   //
   // test erf_inv against data:
   //
   funcp = boost::math::erf_inv;
   result = boost::math::tools::test(
      data, 
      boost::lambda::bind(funcp, 
         boost::lambda::ret<value_type>(boost::lambda::_1[0])), 
      boost::lambda::ret<value_type>(boost::lambda::_1[1]));
   print_test_result(result, data[result.worst_case][0], type_name, "boost::math::erf");
   std::cout << std::endl;
}

template <class T>
void do_test_erfc_inv(const T& data, const char* type_name, const char* test_name)
{
   typedef typename T::value_type row_type;
   typedef typename row_type::value_type value_type;

   typedef value_type (*pg)(value_type);
   pg funcp = boost::math::erf;

   boost::math::tools::test_result<value_type> result;
   std::cout << "Testing " << test_name << " with type " << type_name
      << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
   //
   // test erfc_inv against data:
   //
   funcp = boost::math::erfc_inv;
   result = boost::math::tools::test(
      data, 
      boost::lambda::bind(funcp, 
         boost::lambda::ret<value_type>(boost::lambda::_1[0])), 
      boost::lambda::ret<value_type>(boost::lambda::_1[1]));
   print_test_result(result, data[result.worst_case][0], type_name, "boost::math::erfc_inv");
   std::cout << std::endl;
}

template <class T>
void test_erf(T, const char* name)
{
   //
   // The actual test data is rather verbose, so it's in a separate file
   //
   // The contents are as follows, each row of data contains
   // three items, input value a, input value b and erf(a, b):
   // 
#  include "erf_small_data.ipp"

   do_test_erf(erf_small_data, name, "Erf Function: Small Values");

#  include "erf_data.ipp"

   do_test_erf(erf_data, name, "Erf Function: Medium Values");

#  include "erf_large_data.ipp"

   do_test_erf(erf_large_data, name, "Erf Function: Large Values");

#  include "erf_inv_data.ipp"

   do_test_erf_inv(erf_inv_data, name, "Inverse Erf Function");

#  include "erfc_inv_data.ipp"

   do_test_erfc_inv(erfc_inv_data, name, "Inverse Erfc Function");

}

template <class T>
void test_spots(T, const char* t)
{
   std::cout << "Testing basic sanity checks for type " << t << std::endl;
   //
   // basic sanity checks, tolerance is 10 epsilon expressed as a percentage:
   //
   T tolerance = boost::math::tools::epsilon<T>() * 1000;
   BOOST_CHECK_CLOSE(::boost::math::erfc(static_cast<T>(0.125)), static_cast<T>(0.859683795198666182606970553478), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::erfc(static_cast<T>(0.5)), static_cast<T>(0.479500122186953462317253346108), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::erfc(static_cast<T>(1)), static_cast<T>(0.157299207050285130658779364917), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::erfc(static_cast<T>(5)), static_cast<T>(1.53745979442803485018834348538e-12), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::erfc(static_cast<T>(-0.125)), static_cast<T>(1.14031620480133381739302944652), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::erfc(static_cast<T>(-0.5)), static_cast<T>(1.52049987781304653768274665389), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::erfc(static_cast<T>(0)), static_cast<T>(1), tolerance);

   BOOST_CHECK_CLOSE(::boost::math::erf(static_cast<T>(0.125)), static_cast<T>(0.140316204801333817393029446522), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::erf(static_cast<T>(0.5)), static_cast<T>(0.520499877813046537682746653892), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::erf(static_cast<T>(1)), static_cast<T>(0.842700792949714869341220635083), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::erf(static_cast<T>(5)), static_cast<T>(0.99999999999846254020557196515), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::erf(static_cast<T>(-0.125)), static_cast<T>(-0.140316204801333817393029446522), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::erf(static_cast<T>(-0.5)), static_cast<T>(-0.520499877813046537682746653892), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::erf(static_cast<T>(0)), static_cast<T>(0), tolerance);

   tolerance = boost::math::tools::epsilon<T>() * 100 * 200; // 200 eps %.
#if defined(__CYGWIN__)
   // some platforms long double is only reliably accurate to double precision:
   if(sizeof(T) == sizeof(long double))
      tolerance = boost::math::tools::epsilon<double>() * 100 * 200; // 200 eps %.
#endif
  
   for(T i = -0.95; i < 1; i += 0.125)
   {
      T inv = boost::math::erf_inv(i);
      T b = boost::math::erf(inv);
      BOOST_CHECK_CLOSE(b, i, tolerance);
   }
   for(T j = 0.125; j < 2; j += 0.125)
   {
      T inv = boost::math::erfc_inv(j);
      T b = boost::math::erfc(inv);
      BOOST_CHECK_CLOSE(b, j, tolerance);
   }
}

int test_main(int, char* [])
{
   //test_spots(0.0F);
   test_spots(0.0, "double");
   test_spots(0.0L, "long double");
   test_spots(boost::math::concepts::real_concept(0.1), "real_concept");

   test_erf(0.1F, "float");
   test_erf(0.1, "double");
   test_erf(0.1L, "long double");
   test_erf(boost::math::concepts::real_concept(0.1), "real_concept");
   return 0;
}


