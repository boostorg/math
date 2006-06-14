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

template <class T, class Seq>
void print_test_result(const boost::math::tools::test_result<T>& result, 
                       const Seq& worst, int row, const char* name, const char* test)
{
   using namespace std;
   T eps = pow(T(2), 1-boost::math::tools::digits(worst[0]));
   std::cout << setprecision(4);
   std::cout << test << "(" << name << ") Max = " << (result.stat.max)()/eps
      << " RMS Mean=" << result.stat.rms()/eps;
   if((result.stat.max)() != 0)
   {
      std::cout << "\n    worst case at row: " 
         << row << "\n    { ";
      for(unsigned i = 0; i < worst.size(); ++i)
      {
         if(i)
            std::cout << ", ";
         std::cout << worst[i];
      }
      std::cout << " }";
   }
   std::cout << std::endl;
}


template <class T>
void do_test_beta(const T& data, const char* type_name, const char* test_name)
{
   typedef typename T::value_type row_type;
   typedef typename row_type::value_type value_type;

   typedef value_type (*pg)(value_type, value_type, value_type);
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
         boost::lambda::ret<value_type>(boost::lambda::_1[1]),
         boost::lambda::ret<value_type>(boost::lambda::_1[2])), 
      boost::lambda::ret<value_type>(boost::lambda::_1[3]));
   print_test_result(result, data[result.worst_case], result.worst_case, type_name, "boost::math::beta");
   
   funcp = boost::math::betac;
   result = boost::math::tools::test(
      data, 
      boost::lambda::bind(funcp, 
         boost::lambda::ret<value_type>(boost::lambda::_1[0]),
         boost::lambda::ret<value_type>(boost::lambda::_1[1]),
         boost::lambda::ret<value_type>(boost::lambda::_1[2])), 
      boost::lambda::ret<value_type>(boost::lambda::_1[4]));
   print_test_result(result, data[result.worst_case], result.worst_case, type_name, "boost::math::betac");

   funcp = boost::math::ibeta;
   result = boost::math::tools::test(
      data, 
      boost::lambda::bind(funcp, 
         boost::lambda::ret<value_type>(boost::lambda::_1[0]),
         boost::lambda::ret<value_type>(boost::lambda::_1[1]),
         boost::lambda::ret<value_type>(boost::lambda::_1[2])), 
      boost::lambda::ret<value_type>(boost::lambda::_1[5]));
   print_test_result(result, data[result.worst_case], result.worst_case, type_name, "boost::math::ibeta");
   
   funcp = boost::math::ibetac;
   result = boost::math::tools::test(
      data, 
      boost::lambda::bind(funcp, 
         boost::lambda::ret<value_type>(boost::lambda::_1[0]),
         boost::lambda::ret<value_type>(boost::lambda::_1[1]),
         boost::lambda::ret<value_type>(boost::lambda::_1[2])), 
      boost::lambda::ret<value_type>(boost::lambda::_1[6]));
   print_test_result(result, data[result.worst_case], result.worst_case, type_name, "boost::math::ibetac");
#ifdef TEST_OTHER
   if(::boost::is_floating_point<value_type>::value){
      funcp = other::ibeta;
      result = boost::math::tools::test(
         data, 
         boost::lambda::bind(funcp, 
            boost::lambda::ret<value_type>(boost::lambda::_1[0]),
            boost::lambda::ret<value_type>(boost::lambda::_1[1]),
            boost::lambda::ret<value_type>(boost::lambda::_1[2])), 
         boost::lambda::ret<value_type>(boost::lambda::_1[5]));
      print_test_result(result, data[result.worst_case], result.worst_case, type_name, "other::ibeta");
   }
#endif
   std::cout << std::endl;
}

template <class T>
void test_inverses(const T& data)
{
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
      // incomplete beta is too close to 1: basically there is insuffient 
      // information left in the value we're using as input to the inverse
      // to be able to get back to the original value.
      //
      if(data[i][5] == 0)
         BOOST_CHECK_EQUAL(boost::math::ibeta_inv(data[i][0], data[i][1], data[i][5]), value_type(0));
      else if((1 - data[i][5] > 0.001) && (fabs(data[i][5]) >= boost::math::tools::min_value(data[i][5])))
      {
         value_type inv = boost::math::ibeta_inv(data[i][0], data[i][1], data[i][5]);
         BOOST_CHECK_CLOSE(data[i][2], inv, precision);
      }
      else if(1 == data[i][5])
         BOOST_CHECK_EQUAL(boost::math::ibeta_inv(data[i][0], data[i][1], data[i][5]), value_type(1));

      if(data[i][6] == 0)
         BOOST_CHECK_EQUAL(boost::math::ibetac_inv(data[i][0], data[i][1], data[i][6]), value_type(1));
      else if((1 - data[i][6] > 0.001) && (fabs(data[i][6]) >= boost::math::tools::min_value(data[i][6])))
      {
         value_type inv = boost::math::ibetac_inv(data[i][0], data[i][1], data[i][6]);
         BOOST_CHECK_CLOSE(data[i][2], inv, precision);
      }
      else if(data[i][6] == 1)
         BOOST_CHECK_EQUAL(boost::math::ibetac_inv(data[i][0], data[i][1], data[i][6]), value_type(0));
   }
}

template <class T>
void test_beta(T, const char* name)
{
   //
   // The actual test data is rather verbose, so it's in a separate file
   //
   // The contents are as follows, each row of data contains
   // five items, input value a, input value b, integration limits x, beta(a, b, x) and ibeta(a, b, x):
   // 
#  include "ibeta_small_data.ipp"

   do_test_beta(ibeta_small_data, name, "Incomplete Beta Function: Small Values");
   test_inverses(ibeta_small_data);

#  include "ibeta_data.ipp"

   do_test_beta(ibeta_data, name, "Incomplete Beta Function: Medium Values");
   test_inverses(ibeta_data);

#  include "ibeta_large_data.ipp"

   do_test_beta(ibeta_large_data, name, "Incomplete Beta Function: Large and Diverse Values");
   test_inverses(ibeta_large_data);
}

template <class T>
void test_spots(T)
{
   //
   // basic sanity checks, tolerance is 10 decimal places expressed as a percentage,
   // One check per domain of the implementation:
   //
   T tolerance = std::pow(10.0, -8);
   BOOST_CHECK_CLOSE(
      ::boost::math::ibeta(
         static_cast<T>(0.015964560210704803), 
         static_cast<T>(1.1846856068586931e-005), 
         static_cast<T>(0.69176378846168518)), 
      static_cast<T>(0.0007508604820642986204162462167319506309750), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::ibeta(
         static_cast<T>(42.434902191162109), 
         static_cast<T>(0.30012050271034241), 
         static_cast<T>(0.91574394702911377)), 
      static_cast<T>(0.002844243156314242058287766324242276991912), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::ibeta(
         static_cast<T>(9.7131776809692383), 
         static_cast<T>(99.406852722167969), 
         static_cast<T>(0.083912998437881470)), 
      static_cast<T>(0.4612716118626361034813232775095335302198), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::ibeta(
         static_cast<T>(72.695472717285156), 
         static_cast<T>(1.1902070045471191), 
         static_cast<T>(0.80036874115467072)), 
      static_cast<T>(1.703685144285803673344984949797496197040e-7), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::ibeta(
         static_cast<T>(4.9854421615600586), 
         static_cast<T>(1.0665277242660522), 
         static_cast<T>(0.75997146964073181)), 
      static_cast<T>(0.2755954254731642667260071858810487404614), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::ibeta(
         static_cast<T>(6.8127136230468750), 
         static_cast<T>(1.0562920570373535), 
         static_cast<T>(0.17416560649871826)), 
      static_cast<T>(7.702362015088558153029455563361002570531e-6), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::ibeta(
         static_cast<T>(0.48983201384544373), 
         static_cast<T>(0.22512593865394592), 
         static_cast<T>(0.20032680034637451)), 
      static_cast<T>(0.170905142698145967653807992508983970176), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::ibeta(
         static_cast<T>(4.0498137474060059), 
         static_cast<T>(0.15403440594673157), 
         static_cast<T>(0.65370121598243713)), 
      static_cast<T>(0.0172702040689452906446803217247250156007), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::ibeta(
         static_cast<T>(7.2695474624633789), 
         static_cast<T>(0.11902070045471191), 
         static_cast<T>(0.80036874115467072)), 
      static_cast<T>(0.013346136714187857821168127038816508028), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::ibeta(
         static_cast<T>(2.7266697883605957), 
         static_cast<T>(0.011510574258863926), 
         static_cast<T>(0.086654007434844971)), 
      static_cast<T>(5.812020420972734916187451486321162137375e-6), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::ibeta(
         static_cast<T>(0.34317314624786377), 
         static_cast<T>(0.046342257410287857), 
         static_cast<T>(0.75823287665843964)), 
      static_cast<T>(0.151317265120184850887504097401768195067), tolerance);

   BOOST_CHECK_CLOSE(
      ::boost::math::ibeta(
         static_cast<T>(0.34317314624786377), 
         static_cast<T>(0.046342257410287857), 
         static_cast<T>(0)), 
      static_cast<T>(0), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::ibetac(
         static_cast<T>(0.34317314624786377), 
         static_cast<T>(0.046342257410287857), 
         static_cast<T>(0)), 
      static_cast<T>(1), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::ibeta(
         static_cast<T>(0.34317314624786377), 
         static_cast<T>(0.046342257410287857), 
         static_cast<T>(1)), 
      static_cast<T>(1), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::ibetac(
         static_cast<T>(0.34317314624786377), 
         static_cast<T>(0.046342257410287857), 
         static_cast<T>(1)), 
      static_cast<T>(0), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::ibeta(
         static_cast<T>(1), 
         static_cast<T>(0.046342257410287857), 
         static_cast<T>(0.32)), 
      static_cast<T>(0.0177137046180187568703202426065033413304), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::ibeta(
         static_cast<T>(0.046342257410287857), 
         static_cast<T>(1), 
         static_cast<T>(0.32)), 
      static_cast<T>(0.948565954109602496577407403168592262389), tolerance);
}

int test_main(int, char* [])
{
#ifdef TEST_GSL
   gsl_set_error_handler_off();
#endif
   //test_spots(0.0F);
   test_spots(0.0);
   test_spots(0.0L);
   test_spots(boost::math::concepts::real_concept(0.1));

   test_beta(0.1F, "float");
   test_beta(0.1, "double");
   test_beta(0.1L, "long double");
   test_beta(boost::math::concepts::real_concept(0.1), "real_concept");
   return 0;
}



