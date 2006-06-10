//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/concepts/real_concept.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/test/included/test_exec_monitor.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/math/tools/stats.hpp>
#include <boost/math/tools/test.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/array.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include "test_gamma_hooks.hpp"


template <class T, class Seq>
void print_test_result(const boost::math::tools::test_result<T>& result, 
                       const Seq& worst, int row, const char* name, const char* test)
{
   using namespace std;
   T eps = pow(T(2), 1-boost::math::tools::digits(worst[0]));
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
void do_test_gamma(const T& data, const char* type_name, const char* test_name)
{
   typedef typename T::value_type row_type;
   typedef typename row_type::value_type value_type;

   typedef value_type (*pg)(value_type);
   pg funcp = boost::math::tgamma;

   boost::math::tools::test_result<value_type> result;

   std::cout << "Testing " << test_name << " with type " << type_name
      << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

   //
   // test tgamma against data:
   //
   result = boost::math::tools::test(
      data, 
      boost::lambda::bind(funcp, boost::lambda::ret<value_type>(boost::lambda::_1[0])), 
      boost::lambda::ret<value_type>(boost::lambda::_1[1]));
   print_test_result(result, data[result.worst_case], result.worst_case, type_name, "boost::math::tgamma");
#ifdef TEST_OTHER
   if(::boost::is_floating_point<value_type>::value){
      funcp = other::tgamma;
      result = boost::math::tools::test(
         data, 
         boost::lambda::bind(funcp, boost::lambda::ret<value_type>(boost::lambda::_1[0])), 
         boost::lambda::ret<value_type>(boost::lambda::_1[1]));
      print_test_result(result, data[result.worst_case], result.worst_case, type_name, "other::tgamma");
   }
#endif
   //
   // test lgamma against data:
   //
   funcp = boost::math::lgamma;
   result = boost::math::tools::test(
      data, 
      boost::lambda::bind(funcp, boost::lambda::ret<value_type>(boost::lambda::_1[0])), 
      boost::lambda::ret<value_type>(boost::lambda::_1[2]));
   print_test_result(result, data[result.worst_case], result.worst_case, type_name, "boost::math::lgamma");
#ifdef TEST_OTHER
   if(::boost::is_floating_point<value_type>::value){
      funcp = other::lgamma;
      result = boost::math::tools::test(
         data, 
         boost::lambda::bind(funcp, boost::lambda::ret<value_type>(boost::lambda::_1[0])), 
         boost::lambda::ret<value_type>(boost::lambda::_1[2]));
      print_test_result(result, data[result.worst_case], result.worst_case, type_name, "other::lgamma");
   }
#endif

   std::cout << std::endl;
}

template <class T>
void do_test_gammap1m1(const T& data, const char* type_name, const char* test_name)
{
   typedef typename T::value_type row_type;
   typedef typename row_type::value_type value_type;

   typedef value_type (*pg)(value_type);
   pg funcp = boost::math::tgammap1m1;

   boost::math::tools::test_result<value_type> result;

   std::cout << "Testing " << test_name << " with type " << type_name
      << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

   //
   // test tgammap1m1 against data:
   //
   result = boost::math::tools::test(
      data, 
      boost::lambda::bind(funcp, boost::lambda::ret<value_type>(boost::lambda::_1[0])), 
      boost::lambda::ret<value_type>(boost::lambda::_1[1]));
   print_test_result(result, data[result.worst_case], result.worst_case, type_name, "boost::math::tgammap1m1");
   std::cout << std::endl;
}

template <class T>
void do_test_gamma_2(const T& data, const char* type_name, const char* test_name)
{
   typedef typename T::value_type row_type;
   typedef typename row_type::value_type value_type;

   typedef value_type (*pg)(value_type, value_type);
   pg funcp = boost::math::tgamma;

   using namespace boost::lambda;

   boost::math::tools::test_result<value_type> result;

   std::cout << "Testing " << test_name << " with type " << type_name
      << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

   //
   // test tgamma(T, T) against data:
   //
   if(data[0][2] > 0)
   {
      result = boost::math::tools::test(
         data, 
         bind(funcp, ret<value_type>(_1[0]), ret<value_type>(_1[1])), 
         ret<value_type>(_1[2]));
      print_test_result(result, data[result.worst_case], result.worst_case, type_name, "boost::math::tgamma");
      //
      // test tgamma_lower(T, T) against data:
      //
      funcp = boost::math::tgamma_lower;
      result = boost::math::tools::test(
         data, 
         bind(funcp, ret<value_type>(_1[0]), ret<value_type>(_1[1])), 
         ret<value_type>(_1[4]));
      print_test_result(result, data[result.worst_case], result.worst_case, type_name, "boost::math::tgamma_lower");
   }
   //
   // test gamma_Q(T, T) against data:
   //
   funcp = boost::math::gamma_Q;
   result = boost::math::tools::test(
      data, 
      bind(funcp, ret<value_type>(_1[0]), ret<value_type>(_1[1])), 
      ret<value_type>(_1[3]));
   print_test_result(result, data[result.worst_case], result.worst_case, type_name, "boost::math::gamma_Q");
#if defined(TEST_CEPHES) || defined(TEST_GSL)
   //
   // test other gamma_Q(T, T) against data:
   //
   if(boost::is_floating_point<value_type>::value)
   {
      funcp = other::gamma_Q;
      result = boost::math::tools::test(
         data, 
         bind(funcp, ret<value_type>(_1[0]), ret<value_type>(_1[1])), 
         ret<value_type>(_1[3]));
      print_test_result(result, data[result.worst_case], result.worst_case, type_name, "other::gamma_Q");
   }
#endif
   //
   // test gamma_P(T, T) against data:
   //
   funcp = boost::math::gamma_P;
   result = boost::math::tools::test(
      data, 
      bind(funcp, ret<value_type>(_1[0]), ret<value_type>(_1[1])), 
      ret<value_type>(_1[5]));
   print_test_result(result, data[result.worst_case], result.worst_case, type_name, "boost::math::gamma_P");
#if defined(TEST_CEPHES) || defined(TEST_GSL)
   //
   // test other gamma_P(T, T) against data:
   //
   if(boost::is_floating_point<value_type>::value)
   {
      funcp = other::gamma_P;
      result = boost::math::tools::test(
         data, 
         bind(funcp, ret<value_type>(_1[0]), ret<value_type>(_1[1])), 
         ret<value_type>(_1[5]));
      print_test_result(result, data[result.worst_case], result.worst_case, type_name, "other::gamma_P");
   }
#endif
   //
   // test gamma_P_inv(T, T) against data:
   //
   using namespace std;
   typedef typename T::value_type row_type;
   typedef typename row_type::value_type value_type;

   value_type precision = static_cast<value_type>(ldexp(1.0, 1-boost::math::tools::digits(value_type(0))/2)) * 100;
   if(boost::math::tools::digits(value_type(0)) < 50)
      precision = 1;   // 1% or two decimal digits, all we can hope for when the input is truncated

   for(unsigned i = 0; i < data.size(); ++i)
   {
      //
      // These inverse tests are thrown off if the output of the 
      // incomplete gamma is too close to 1: basically there is insuffient 
      // information left in the value we're using as input to the inverse
      // to be able to get back to the original value.
      //
      if(data[i][5] == 0)
         BOOST_CHECK_EQUAL(boost::math::gamma_P_inv(data[i][0], data[i][5]), value_type(0));
      else if((1 - data[i][5] > 0.001) && (fabs(data[i][5]) >= boost::math::tools::min_value(data[i][5])))
      {
         value_type inv = boost::math::gamma_P_inv(data[i][0], data[i][5]);
         BOOST_CHECK_CLOSE(data[i][1], inv, precision);
      }
      else if(1 == data[i][5])
         BOOST_CHECK_EQUAL(boost::math::gamma_P_inv(data[i][0], data[i][5]), boost::math::tools::max_value(data[i][5]));
      else
      {
         // not enough bits in our input to get back to x, but we should be in
         // the same ball park:
         value_type inv = boost::math::gamma_P_inv(data[i][0], data[i][5]);
         BOOST_CHECK_CLOSE(data[i][1], inv, 100000);
      }

      if(data[i][3] == 0)
         BOOST_CHECK_EQUAL(boost::math::gamma_Q_inv(data[i][0], data[i][3]), boost::math::tools::max_value(data[i][3]));
      else if((1 - data[i][3] > 0.001) && (fabs(data[i][3]) >= boost::math::tools::min_value(data[i][3])))
      {
         value_type inv = boost::math::gamma_Q_inv(data[i][0], data[i][3]);
         BOOST_CHECK_CLOSE(data[i][1], inv, precision);
      }
      else if(1 == data[i][3])
         BOOST_CHECK_EQUAL(boost::math::gamma_Q_inv(data[i][0], data[i][3]), value_type(0));
      else
      {
         // not enough bits in our input to get back to x, but we should be in
         // the same ball park:
         value_type inv = boost::math::gamma_Q_inv(data[i][0], data[i][3]);
         BOOST_CHECK_CLOSE(data[i][1], inv, 100);
      }
   }
   std::cout << std::endl;
}

template <class T>
void do_test_gamma_inv(const T& data, const char* type_name, const char* test_name)
{
   typedef typename T::value_type row_type;
   typedef typename row_type::value_type value_type;

   typedef value_type (*pg)(value_type, value_type);
   pg funcp = boost::math::gamma_P_inv;

   using namespace boost::lambda;

   boost::math::tools::test_result<value_type> result;

   std::cout << "Testing " << test_name << " with type " << type_name
      << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

   //
   // test gamma_P_inv(T, T) against data:
   //
   result = boost::math::tools::test(
      data, 
      bind(funcp, ret<value_type>(_1[0]), ret<value_type>(_1[1])), 
      ret<value_type>(_1[2]));
   print_test_result(result, data[result.worst_case], result.worst_case, type_name, "boost::math::gamma_P_inv");
   //
   // test gamma_Q_inv(T, T) against data:
   //
   funcp = boost::math::gamma_Q_inv;
   result = boost::math::tools::test(
      data, 
      bind(funcp, ret<value_type>(_1[0]), ret<value_type>(_1[1])), 
      ret<value_type>(_1[3]));
   print_test_result(result, data[result.worst_case], result.worst_case, type_name, "boost::math::gamma_Q_inv");
}
   
template <class T>
void test_gamma(T, const char* name)
{
   //
   // The actual test data is rather verbose, so it's in a separate file
   //
   // The contents are as follows, each row of data contains
   // three items, input value, gamma and lgamma:
   // 
   // gamma and lgamma at integer and half integer values:
   // boost::array<boost::array<T, 3>, N> factorials;
   //
   // gamma and lgamma for z near 0:
   // boost::array<boost::array<T, 3>, N> near_0;
   //
   // gamma and lgamma for z near 1:
   // boost::array<boost::array<T, 3>, N> near_1;
   //
   // gamma and lgamma for z near 2:
   // boost::array<boost::array<T, 3>, N> near_2;
   //
   // gamma and lgamma for z near -10:
   // boost::array<boost::array<T, 3>, N> near_m10;
   //
   // gamma and lgamma for z near -55:
   // boost::array<boost::array<T, 3>, N> near_m55;
   //
   // The last two cases are chosen more or less at random,
   // except that one is even and the other odd, and both are
   // at negative poles.  The data near zero also tests near
   // a pole, the data near 1 and 2 are to probe lgamma as
   // the result -> 0.
   //
#  include "test_gamma_data.ipp"

   do_test_gamma(factorials, name, "factorials");
   do_test_gamma(near_0, name, "near 0");
   do_test_gamma(near_1, name, "near 1");
   do_test_gamma(near_2, name, "near 2");
   do_test_gamma(near_m10, name, "near -10");
   do_test_gamma(near_m55, name, "near -55"); 

   //
   // And now tgammap1m1 which computes gamma(1+dz)-1:
   //
   do_test_gammap1m1(gammap1m1_data, name, "tgammap1m1(dz)");

   //
   // Now the data for the incomplete gamma function, each
   // row has the following entries:
   // Parameter a, parameter z,
   // Expected tgamma(a, z), Expected gamma_Q(a, z)
   // Expected tgamma_lower(a, z), Expected gamma_P(a, z)
   //
#  include "igamma_med_data.ipp"

   do_test_gamma_2(igamma_med_data, name, "tgamma(a, z) medium values");

#  include "igamma_small_data.ipp"

   do_test_gamma_2(igamma_small_data, name, "tgamma(a, z) small values");

#  include "igamma_big_data.ipp"

   do_test_gamma_2(igamma_big_data, name, "tgamma(a, z) large values");

#  include "gamma_inv_data.ipp"

   do_test_gamma_inv(gamma_inv_data, name, "incomplete gamma inverse(a, z) medium values");

#  include "gamma_inv_big_data.ipp"

   do_test_gamma_inv(gamma_inv_big_data, name, "incomplete gamma inverse(a, z) large values");

#  include "gamma_inv_small_data.ipp"

   do_test_gamma_inv(gamma_inv_small_data, name, "incomplete gamma inverse(a, z) small values");
}

template <class T>
void test_spots(T)
{
   //
   // basic sanity checks, tolerance is 10 decimal places expressed as a percentage:
   //
   T tolerance = std::pow(10.0, -8);
   BOOST_CHECK_CLOSE(::boost::math::tgamma(static_cast<T>(3.5)), static_cast<T>(3.3233509704478425512), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::tgamma(static_cast<T>(0.1)), static_cast<T>(9.5135076986687318363), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::tgamma(static_cast<T>(-0.1)), static_cast<T>(-10.686287021193193549), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::tgamma(static_cast<T>(-3.2)), static_cast<T>(0.68905641200597974292), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::tgamma(static_cast<T>(-52.001)), static_cast<T>(-1.2349010280547181826e-65), tolerance);

   int sign = 1;
   BOOST_CHECK_CLOSE(::boost::math::lgamma(static_cast<T>(3.5), &sign), static_cast<T>(1.2009736023470742248208170321516), tolerance);
   BOOST_CHECK(sign == 1);
   BOOST_CHECK_CLOSE(::boost::math::lgamma(static_cast<T>(0.1), &sign), static_cast<T>(2.2527126517342059598704913469682), tolerance);
   BOOST_CHECK(sign == 1);
   BOOST_CHECK_CLOSE(::boost::math::lgamma(static_cast<T>(-0.1), &sign), static_cast<T>(2.3689613327287886552092298040686), tolerance);
   BOOST_CHECK(sign == -1);
   BOOST_CHECK_CLOSE(::boost::math::lgamma(static_cast<T>(-3.2), &sign), static_cast<T>(-0.37243213612996769536089336752751), tolerance);
   BOOST_CHECK(sign == 1);
   BOOST_CHECK_CLOSE(::boost::math::lgamma(static_cast<T>(-52.001), &sign), static_cast<T>(-149.4570402169712972001218882919), tolerance);
   BOOST_CHECK(sign == -1);

   BOOST_CHECK_CLOSE(::boost::math::tgamma(static_cast<T>(5), static_cast<T>(1)), static_cast<T>(23.912163676143750904), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::tgamma(static_cast<T>(5), static_cast<T>(5)), static_cast<T>(10.571838841565097875), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::tgamma(static_cast<T>(5), static_cast<T>(10)), static_cast<T>(0.70206451384706574415), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::tgamma(static_cast<T>(5), static_cast<T>(100)), static_cast<T>(3.8734332808745531497e-36), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::tgamma(static_cast<T>(0.5), static_cast<T>(0.5)), static_cast<T>(0.562418231594407124279494957302), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::tgamma(static_cast<T>(0.5), static_cast<T>(0.9)), static_cast<T>(0.318532103604121098738593603904), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::tgamma(static_cast<T>(0.5), static_cast<T>(5)), static_cast<T>(0.00277460326041280931949083572726), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::tgamma(static_cast<T>(0.5), static_cast<T>(100)), static_cast<T>(3.70174786040827892025356644813e-45), tolerance);

   BOOST_CHECK_CLOSE(::boost::math::tgamma_lower(static_cast<T>(5), static_cast<T>(1)), static_cast<T>(0.087836323856249096291), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::tgamma_lower(static_cast<T>(5), static_cast<T>(5)), static_cast<T>(13.428161158434902125), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::tgamma_lower(static_cast<T>(5), static_cast<T>(10)), static_cast<T>(23.297935486152934256), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::tgamma_lower(static_cast<T>(5), static_cast<T>(100)), static_cast<T>(24.000000000000000000), tolerance);

   BOOST_CHECK_CLOSE(::boost::math::gamma_Q(static_cast<T>(5), static_cast<T>(1)), static_cast<T>(0.996340153172656287654543544187), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::gamma_Q(static_cast<T>(5), static_cast<T>(5)), static_cast<T>(0.440493285065212411442581665663), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::gamma_Q(static_cast<T>(5), static_cast<T>(10)), static_cast<T>(0.0292526880769610726727661331928), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::gamma_Q(static_cast<T>(5), static_cast<T>(100)), static_cast<T>(1.61393053369773047904057392250e-37), tolerance);

   BOOST_CHECK_CLOSE(::boost::math::gamma_P(static_cast<T>(5), static_cast<T>(1)), static_cast<T>(0.00365984682734371234545645581271), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::gamma_P(static_cast<T>(5), static_cast<T>(5)), static_cast<T>(0.559506714934787588557418334337), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::gamma_P(static_cast<T>(5), static_cast<T>(10)), static_cast<T>(0.970747311923038927327233866807), tolerance);
   BOOST_CHECK_CLOSE(::boost::math::gamma_P(static_cast<T>(5), static_cast<T>(100)), static_cast<T>(1.00000000000000000000000000000), tolerance);

   //BOOST_CHECK_CLOSE(::boost::math::tgamma(static_cast<T>(5), static_cast<T>(1), static_cast<T>(1.5)), static_cast<T>(0.35798614547512708694), tolerance);
   //BOOST_CHECK_CLOSE(::boost::math::tgamma(static_cast<T>(5), static_cast<T>(3), static_cast<T>(6)), static_cast<T>(12.724961860971380343), tolerance);
   //BOOST_CHECK_CLOSE(::boost::math::tgamma(static_cast<T>(5), static_cast<T>(7), static_cast<T>(9)), static_cast<T>(2.8326711932871945905), tolerance);

   tolerance = boost::math::tools::epsilon(tolerance) * 100 * 200; // 200 eps %.
#if defined(__CYGWIN__)
   // some platforms long double is only reliably accurate to double precision:
   if(sizeof(T) == sizeof(long double))
      tolerance = boost::math::tools::epsilon(double(0)) * 100 * 200; // 200 eps %.
#endif
#if 0
   // this is way too slow...
   for(T a = 0.5; a < 500; a += 2.5)
   {
      for(T p = 0.05; p < 1; p += 0.04)
      {
         T x = boost::math::gamma_P_inv(a, p);
         T p2 = boost::math::gamma_P(a, x);
         BOOST_CHECK_CLOSE(p2, p, tolerance);

         x = boost::math::gamma_Q_inv(a, p);
         p2 = boost::math::gamma_Q(a, x);
         BOOST_CHECK_CLOSE(p2, p, tolerance);
      }
   }
#endif
}

int test_main(int, char* [])
{
   //test_spots(0.0F);
   test_spots(0.0);
   test_spots(0.0L);
   test_spots(boost::math::concepts::real_concept(0.1));

   test_gamma(0.1F, "float");
   test_gamma(0.1, "double");
   test_gamma(0.1L, "long double");
   test_gamma(boost::math::concepts::real_concept(0.1), "real_concept");
   return 0;
}


