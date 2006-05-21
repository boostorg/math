//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/concepts/real_concept.hpp>
#include <boost/test/included/test_exec_monitor.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/math/tools/stats.hpp>
#include <boost/math/tools/test.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/array.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/math/special_functions/cbrt.hpp>

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
void do_test_cbrt(const T& data, const char* type_name, const char* test_name)
{
   typedef typename T::value_type row_type;
   typedef typename row_type::value_type value_type;

   typedef value_type (*pg)(value_type);
   pg funcp = boost::math::cbrt;

   boost::math::tools::test_result<value_type> result;

   std::cout << "Testing " << test_name << " with type " << type_name
      << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

   //
   // test erf against data:
   //
   result = boost::math::tools::test(
      data, 
      boost::lambda::bind(funcp, 
         boost::lambda::ret<value_type>(boost::lambda::_1[1])), 
      boost::lambda::ret<value_type>(boost::lambda::_1[0]));
   result += boost::math::tools::test(
      data, 
      boost::lambda::bind(funcp, 
         -boost::lambda::ret<value_type>(boost::lambda::_1[1])), 
      -boost::lambda::ret<value_type>(boost::lambda::_1[0]));
   print_test_result(result, data[result.worst_case][1], type_name, "boost::math::cbrt");
   std::cout << std::endl;
}
template <class T>
void test_cbrt(T, const char* name)
{
   //
   // The actual test data is rather verbose, so it's in a separate file
   //
   // The contents are as follows, each row of data contains
   // three items, input value a, input value b and erf(a, b):
   // 
#  include "cbrt_data.ipp"

   do_test_cbrt(cbrt_data, name, "cbrt Function");

}

template <class T>
void test_spots(T)
{
   //
   // basic sanity checks, tolerance is 10 decimal places expressed as a persentage:
   //
   T tolerance = std::pow(10.0, -8);
   //BOOST_CHECK_CLOSE(::boost::math::erfc(static_cast<T>(0.125)), static_cast<T>(0.859683795198666182606970553478), tolerance);
}

int test_main(int, char* [])
{
   //test_spots(0.0F);
   test_spots(0.0);
   test_spots(0.0L);
   test_spots(boost::math::concepts::real_concept(0.1));

   test_cbrt(0.1F, "float");
   test_cbrt(0.1, "double");
   test_cbrt(0.1L, "long double");
   test_cbrt(boost::math::concepts::real_concept(0.1), "real_concept");
   return 0;
}


