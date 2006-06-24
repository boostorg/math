//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#include <boost/math/concepts/real_concept.hpp>
#include <boost/test/included/test_exec_monitor.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/tools/stats.hpp>
#include <boost/math/tools/test.hpp>
#include <boost/array.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

template <class T, class Seq>
void print_test_result(const boost::math::tools::test_result<T>& result, 
                       const Seq& worst, int row, const char* name, const char* test)
{
   using namespace std;
   T eps = pow(T(2), 1-boost::math::tools::digits<T>());
   std::cout << setprecision(4);
   std::cout << test << "(" << name << ") Max = " << (result.stat.max)()/eps
      << " RMS Mean=" << result.stat.rms()/eps 
      << "\n    worst case at row: " 
      << row << "\n    { ";
   for(unsigned i = 0; i < worst.size(); ++i)
   {
      if(i)
         std::cout << ", ";
      std::cout << worst[i];
   }
   std::cout << " }" << std::endl;
}

template <class T>
void do_test_tgamma_delta_ratio(const T& data, const char* type_name, const char* test_name)
{
   typedef typename T::value_type row_type;
   typedef typename row_type::value_type value_type;

   typedef value_type (*pg)(value_type, value_type);
   pg funcp = boost::math::tgamma_delta_ratio;

   boost::math::tools::test_result<value_type> result;

   std::cout << "Testing " << test_name << " with type " << type_name
      << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

   //
   // test tgamma_delta_ratio against data:
   //
   result = boost::math::tools::test(
      data, 
      boost::lambda::bind(funcp, 
         boost::lambda::ret<value_type>(boost::lambda::_1[0]),
         boost::lambda::ret<value_type>(boost::lambda::_1[1])), 
      boost::lambda::ret<value_type>(boost::lambda::_1[2]));
   print_test_result(result, data[result.worst_case], result.worst_case, type_name, "boost::math::tgamma_delta_ratio(a, delta)");
   result = boost::math::tools::test(
      data, 
      boost::lambda::bind(funcp, 
         boost::lambda::ret<value_type>(boost::lambda::_1[0]),
         -boost::lambda::ret<value_type>(boost::lambda::_1[1])), 
      boost::lambda::ret<value_type>(boost::lambda::_1[3]));
   print_test_result(result, data[result.worst_case], result.worst_case, type_name, "boost::math::tgamma_delta_ratio(a -delta)");
}

template <class T>
void do_test_tgamma_ratio(const T& data, const char* type_name, const char* test_name)
{
   typedef typename T::value_type row_type;
   typedef typename row_type::value_type value_type;

   typedef value_type (*pg)(value_type, value_type);
   pg funcp = boost::math::tgamma_ratio;

   boost::math::tools::test_result<value_type> result;

   std::cout << "Testing " << test_name << " with type " << type_name
      << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

   //
   // test tgamma_ratio against data:
   //
   result = boost::math::tools::test(
      data, 
      boost::lambda::bind(funcp, 
         boost::lambda::ret<value_type>(boost::lambda::_1[0]),
         boost::lambda::ret<value_type>(boost::lambda::_1[1])), 
      boost::lambda::ret<value_type>(boost::lambda::_1[2]));
   print_test_result(result, data[result.worst_case], result.worst_case, type_name, "boost::math::tgamma_ratio(a, b)");
}

template <class T>
void test_tgamma_ratio(T, const char* name)
{
   //
   // The actual test data is rather verbose, so it's in a separate file
   //
#  include "tgamma_delta_ratio_data.ipp"

   do_test_tgamma_delta_ratio(tgamma_delta_ratio_data, name, "tgamma ratios");

#  include "tgamma_ratio_data.ipp"

   do_test_tgamma_ratio(tgamma_ratio_data, name, "tgamma ratios");

}

int test_main(int, char* [])
{
   test_tgamma_ratio(0.1F, "float");
   test_tgamma_ratio(0.1, "double");
   test_tgamma_ratio(0.1L, "long double");
   test_tgamma_ratio(boost::math::concepts::real_concept(0.1), "real_concept");
   return 0;
}

