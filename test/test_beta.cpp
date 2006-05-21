//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/concepts/real_concept.hpp>
#include <boost/test/included/test_exec_monitor.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/tools/stats.hpp>
#include <boost/math/tools/test.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/array.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include "test_beta_hooks.hpp"


template <class T>
void print_test_result(const boost::math::tools::test_result<T>& result, 
                       T worst, const char* name, const char* test)
{
   using namespace std;
   T eps = pow(T(2), 1-boost::math::tools::digits(worst));
   std::cout << setprecision(4);
   std::cout << test << "(" << name << ") Max = " << (result.stat.max)()/eps
      << " RMS Mean=" << result.stat.rms()/eps << " worst case at point: " << worst << std::endl;
}

template <class T>
void do_test_beta(const T& data, const char* type_name, const char* test_name)
{
   typedef typename T::value_type row_type;
   typedef typename row_type::value_type value_type;

   typedef value_type (*pg)(value_type, value_type);
   pg funcp = boost::math::beta;

   boost::math::tools::test_result<value_type> result;

   std::cout << "Testing " << test_name << " with type " << type_name
      << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

   //
   // test beta against data:
   //
   result = boost::math::tools::test(
      data, 
      boost::lambda::bind(funcp, 
         boost::lambda::ret<value_type>(boost::lambda::_1[0]),
         boost::lambda::ret<value_type>(boost::lambda::_1[1])), 
      boost::lambda::ret<value_type>(boost::lambda::_1[2]));
   print_test_result(result, data[result.worst_case][0], type_name, "boost::math::beta");
#ifdef TEST_OTHER
   if(::boost::is_floating_point<value_type>::value){
      funcp = other::beta;
      result = boost::math::tools::test(
         data, 
         boost::lambda::bind(funcp, 
            boost::lambda::ret<value_type>(boost::lambda::_1[0]),
            boost::lambda::ret<value_type>(boost::lambda::_1[1])), 
         boost::lambda::ret<value_type>(boost::lambda::_1[2]));
      print_test_result(result, data[result.worst_case][0], type_name, "other::beta");
   }
#endif
   std::cout << std::endl;
}
template <class T>
void test_beta(T, const char* name)
{
   //
   // The actual test data is rather verbose, so it's in a separate file
   //
   // The contents are as follows, each row of data contains
   // three items, input value a, input value b and beta(a, b):
   // 
#  include "beta_small_data.ipp"

   do_test_beta(beta_small_data, name, "Beta Function: Small Values");

#  include "beta_med_data.ipp"

   do_test_beta(beta_med_data, name, "Beta Function: Medium Values");

#  include "beta_exp_data.ipp"

   do_test_beta(beta_exp_data, name, "Beta Function: Divergent Values");
}

template <class T>
void test_spots(T)
{
   //
   // basic sanity checks, tolerance is 10 decimal places expressed as a persentage:
   //
   T tolerance = std::pow(10.0, -8);
   BOOST_CHECK_CLOSE(::boost::math::beta(static_cast<T>(1), static_cast<T>(1)), static_cast<T>(1), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::beta(static_cast<T>(1), static_cast<T>(4)), static_cast<T>(0.25), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::beta(static_cast<T>(4), static_cast<T>(1)), static_cast<T>(0.25), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::beta(static_cast<T>(1e-50), static_cast<T>(4)), static_cast<T>(1e50), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::beta(static_cast<T>(4), static_cast<T>(1e-50)), static_cast<T>(1e50), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::beta(static_cast<T>(4), static_cast<T>(20)), static_cast<T>(0.00002823263692828910220214568040654997176736), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::beta(static_cast<T>(0.0125), static_cast<T>(0.000023)), static_cast<T>(43558.24045647538375006349016083320744662), tolerance);
}

int test_main(int, char* [])
{
   //test_spots(0.0F);
   test_spots(0.0);
   test_spots(0.0L);
   test_spots(boost::math::concepts::real_concept(0.1));
#ifndef __HP_aCC
   test_beta(0.1F, "float");
   test_beta(0.1, "double");
#endif
   test_beta(0.1L, "long double");
   test_beta(boost::math::concepts::real_concept(0.1), "real_concept");
   return 0;
}


