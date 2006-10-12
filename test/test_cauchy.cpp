// Copyright John Maddock 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// test_students_t.cpp

#define BOOST_MATH_THROW_ON_DOMAIN_ERROR
#define BOOST_MATH_THROW_ON_OVERFLOW_ERROR

#ifdef _MSC_VER
#  pragma warning(disable: 4127) // conditional expression is constant.
#  pragma warning(disable: 4100) // unreferenced formal parameter.
#  pragma warning(disable: 4512) // assignment operator could not be generated.
#  pragma warning(disable: 4510) // default constructor could not be generated.
#  pragma warning(disable: 4610) // can never be instantiated - user defined constructor required.
//#  pragma warning(disable: 4535) // calling _set_se_translator() requires /EHa (in Boost.test)
// Enable C++ Exceptions Yes With SEH Exceptions (/EHa) prevents warning 4535.
#endif

#include <boost/math/concepts/real_concept.hpp> // for real_concept
#include <boost/math/distributions/cauchy.hpp>
	 using boost::math::cauchy_distribution;

#include <boost/test/included/test_exec_monitor.hpp> // Boost.Test
#include <boost/test/floating_point_comparison.hpp>


#include <iostream>
	using std::cout;
	using std::endl;
	using std::setprecision;

template <class RealType>
void test_spots(RealType T)
{
   // Basic sanity checks.
   // 50eps as a persentage, up to a maximum of double precision
   // (that's the limit of our test data).
   RealType tolerance = (std::max)(
      static_cast<RealType>(boost::math::tools::epsilon<double>()),
      boost::math::tools::epsilon<RealType>());
   tolerance *= 50 * 100;  

	cout << "Tolerance for type " << typeid(T).name()  << " is " << tolerance << " %" << endl;

   //
   // CDF:
   //
   BOOST_CHECK_CLOSE(
      ::boost::math::cdf(
         cauchy_distribution<RealType>(),      
         static_cast<RealType>(0.125)),              // x
         static_cast<RealType>(0.53958342416056554201085167134004L),                // probability.
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::cdf(
         cauchy_distribution<RealType>(),      
         static_cast<RealType>(-0.125)),              // x
         static_cast<RealType>(0.46041657583943445798914832865996L),                // probability.
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::cdf(
         cauchy_distribution<RealType>(),      
         static_cast<RealType>(0.5)),              // x
         static_cast<RealType>(0.64758361765043327417540107622474L),                // probability.
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::cdf(
         cauchy_distribution<RealType>(),      
         static_cast<RealType>(-0.5)),              // x
         static_cast<RealType>(0.35241638234956672582459892377526),                // probability.
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::cdf(
         cauchy_distribution<RealType>(),      
         static_cast<RealType>(1.0)),              // x
         static_cast<RealType>(0.75),                // probability.
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::cdf(
         cauchy_distribution<RealType>(),      
         static_cast<RealType>(-1.0)),              // x
         static_cast<RealType>(0.25),                // probability.
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::cdf(
         cauchy_distribution<RealType>(),      
         static_cast<RealType>(2.0)),              // x
         static_cast<RealType>(0.85241638234956672582459892377526L),                // probability.
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::cdf(
         cauchy_distribution<RealType>(),      
         static_cast<RealType>(-2.0)),              // x
         static_cast<RealType>(0.14758361765043327417540107622474L),                // probability.
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::cdf(
         cauchy_distribution<RealType>(),      
         static_cast<RealType>(10.0)),              // x
         static_cast<RealType>(0.9682744825694464304850228813987L),                // probability.
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::cdf(
         cauchy_distribution<RealType>(),      
         static_cast<RealType>(-10.0)),              // x
         static_cast<RealType>(0.031725517430553569514977118601302L),                // probability.
			tolerance); // %

   //
   // Complements:
   //
   BOOST_CHECK_CLOSE(
      ::boost::math::cdf(
         complement(cauchy_distribution<RealType>(),      
         static_cast<RealType>(0.125))),              // x
         static_cast<RealType>(0.46041657583943445798914832865996L),                // probability.
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::cdf(
         complement(cauchy_distribution<RealType>(),      
         static_cast<RealType>(-0.125))),              // x
         static_cast<RealType>(0.53958342416056554201085167134004L),                // probability.
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::cdf(
         complement(cauchy_distribution<RealType>(),      
         static_cast<RealType>(0.5))),              // x
         static_cast<RealType>(0.35241638234956672582459892377526L),                // probability.
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::cdf(
         complement(cauchy_distribution<RealType>(),      
         static_cast<RealType>(-0.5))),              // x
         static_cast<RealType>(0.64758361765043327417540107622474L),                // probability.
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::cdf(
         complement(cauchy_distribution<RealType>(),      
         static_cast<RealType>(1.0))),              // x
         static_cast<RealType>(0.25),                // probability.
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::cdf(
         complement(cauchy_distribution<RealType>(),      
         static_cast<RealType>(-1.0))),              // x
         static_cast<RealType>(0.75),                // probability.
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::cdf(
         complement(cauchy_distribution<RealType>(),      
         static_cast<RealType>(2.0))),              // x
         static_cast<RealType>(0.14758361765043327417540107622474L),                // probability.
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::cdf(
         complement(cauchy_distribution<RealType>(),      
         static_cast<RealType>(-2.0))),              // x
         static_cast<RealType>(0.85241638234956672582459892377526L),                // probability.
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::cdf(
         complement(cauchy_distribution<RealType>(),      
         static_cast<RealType>(10.0))),              // x
         static_cast<RealType>(0.031725517430553569514977118601302L),                // probability.
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::cdf(
         complement(cauchy_distribution<RealType>(),      
         static_cast<RealType>(-10.0))),              // x
         static_cast<RealType>(0.9682744825694464304850228813987L),                // probability.
			tolerance); // %

   //
   // Quantiles:
   //
   BOOST_CHECK_CLOSE(
      ::boost::math::quantile(
         cauchy_distribution<RealType>(),      
         static_cast<RealType>(0.53958342416056554201085167134004L)),              
         static_cast<RealType>(0.125),              
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::quantile(
         cauchy_distribution<RealType>(),      
         static_cast<RealType>(0.46041657583943445798914832865996L)),              
         static_cast<RealType>(-0.125),      
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::quantile(
         cauchy_distribution<RealType>(),      
         static_cast<RealType>(0.64758361765043327417540107622474L)),              
         static_cast<RealType>(0.5),       
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::quantile(
         cauchy_distribution<RealType>(),      
         static_cast<RealType>(0.35241638234956672582459892377526)),              
         static_cast<RealType>(-0.5),        
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::quantile(
         cauchy_distribution<RealType>(),      
         static_cast<RealType>(0.75)),              
         static_cast<RealType>(1.0),      
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::quantile(
         cauchy_distribution<RealType>(),      
         static_cast<RealType>(0.25)),              
         static_cast<RealType>(-1.0),         
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::quantile(
         cauchy_distribution<RealType>(),      
         static_cast<RealType>(0.85241638234956672582459892377526L)),              
         static_cast<RealType>(2.0),       
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::quantile(
         cauchy_distribution<RealType>(),      
         static_cast<RealType>(0.14758361765043327417540107622474L)),              
         static_cast<RealType>(-2.0),      
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::quantile(
         cauchy_distribution<RealType>(),      
         static_cast<RealType>(0.9682744825694464304850228813987L)),              
         static_cast<RealType>(10.0),          
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::quantile(
         cauchy_distribution<RealType>(),      
         static_cast<RealType>(0.031725517430553569514977118601302L)),              
         static_cast<RealType>(-10.0),           
			tolerance); // %

   //
   // Quantile from complement:
   //
   BOOST_CHECK_CLOSE(
      ::boost::math::quantile(
         complement(cauchy_distribution<RealType>(),      
         static_cast<RealType>(0.46041657583943445798914832865996L))),              
         static_cast<RealType>(0.125),      
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::quantile(
         complement(cauchy_distribution<RealType>(),      
         static_cast<RealType>(0.53958342416056554201085167134004L))),              
         static_cast<RealType>(-0.125),     
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::quantile(
         complement(cauchy_distribution<RealType>(),      
         static_cast<RealType>(0.35241638234956672582459892377526L))),              
         static_cast<RealType>(0.5),   
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::quantile(
         complement(cauchy_distribution<RealType>(),      
         static_cast<RealType>(0.64758361765043327417540107622474L))),              
         static_cast<RealType>(-0.5),    
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::quantile(
         complement(cauchy_distribution<RealType>(),      
         static_cast<RealType>(0.25))),              
         static_cast<RealType>(1.0), 
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::quantile(
         complement(cauchy_distribution<RealType>(),      
         static_cast<RealType>(0.75))),              
         static_cast<RealType>(-1.0), 
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::quantile(
         complement(cauchy_distribution<RealType>(),      
         static_cast<RealType>(0.14758361765043327417540107622474L))),              
         static_cast<RealType>(2.0),
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::quantile(
         complement(cauchy_distribution<RealType>(),      
         static_cast<RealType>(0.85241638234956672582459892377526L))),              
         static_cast<RealType>(-2.0), 
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::quantile(
         complement(cauchy_distribution<RealType>(),      
         static_cast<RealType>(0.031725517430553569514977118601302L))),              
         static_cast<RealType>(10.0), 
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::quantile(
         complement(cauchy_distribution<RealType>(),      
         static_cast<RealType>(0.9682744825694464304850228813987L))),              
         static_cast<RealType>(-10.0),
			tolerance); // %

   //
   // PDF
   //
   BOOST_CHECK_CLOSE(
      ::boost::math::pdf(
         cauchy_distribution<RealType>(),      
         static_cast<RealType>(0.125)),              // x
         static_cast<RealType>(0.31341281101173235351410956479511L),                // probability.
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::pdf(
         cauchy_distribution<RealType>(),      
         static_cast<RealType>(-0.125)),              // x
         static_cast<RealType>(0.31341281101173235351410956479511L),                // probability.
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::pdf(
         cauchy_distribution<RealType>(),      
         static_cast<RealType>(0.5)),              // x
         static_cast<RealType>(0.25464790894703253723021402139602L),                // probability.
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::pdf(
         cauchy_distribution<RealType>(),      
         static_cast<RealType>(-0.5)),              // x
         static_cast<RealType>(0.25464790894703253723021402139602L),                // probability.
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::pdf(
         cauchy_distribution<RealType>(),      
         static_cast<RealType>(1.0)),              // x
         static_cast<RealType>(0.15915494309189533576888376337251L),                // probability.
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::pdf(
         cauchy_distribution<RealType>(),      
         static_cast<RealType>(-1.0)),              // x
         static_cast<RealType>(0.15915494309189533576888376337251L),                // probability.
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::pdf(
         cauchy_distribution<RealType>(),      
         static_cast<RealType>(2.0)),              // x
         static_cast<RealType>(0.063661977236758134307553505349006L),                // probability.
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::pdf(
         cauchy_distribution<RealType>(),      
         static_cast<RealType>(-2.0)),              // x
         static_cast<RealType>(0.063661977236758134307553505349006L),                // probability.
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::pdf(
         cauchy_distribution<RealType>(),      
         static_cast<RealType>(10.0)),              // x
         static_cast<RealType>(0.0031515830315226799162155200667825L),                // probability.
			tolerance); // %
   BOOST_CHECK_CLOSE(
      ::boost::math::pdf(
         cauchy_distribution<RealType>(),      
         static_cast<RealType>(-10.0)),              // x
         static_cast<RealType>(0.0031515830315226799162155200667825L),                // probability.
			tolerance); // %


   //
   // Things that are errors:
   //
   cauchy_distribution<RealType> dist;
   BOOST_CHECK_THROW(
       mean(dist),
       std::domain_error);
   BOOST_CHECK_THROW(
       variance(dist),
       std::domain_error);
   BOOST_CHECK_THROW(
       standard_deviation(dist),
       std::domain_error);
   BOOST_CHECK_THROW(
       quantile(dist, RealType(0.0)),
       std::overflow_error);
   BOOST_CHECK_THROW(
       quantile(dist, RealType(1.0)),
       std::overflow_error);
   BOOST_CHECK_THROW(
       quantile(complement(dist, RealType(0.0))),
       std::overflow_error);
   BOOST_CHECK_THROW(
       quantile(complement(dist, RealType(1.0))),
       std::overflow_error);


} // template <class RealType>void test_spots(RealType)

int test_main(int, char* [])
{
	 // Basic sanity-check spot values.
	// (Parameter value, arbitrarily zero, only communicates the floating point type).
  test_spots(0.0F); // Test float. OK at decdigits = 0 tolerance = 0.0001 %
  test_spots(0.0); // Test double. OK at decdigits 7, tolerance = 1e07 %
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
  test_spots(0.0L); // Test long double.
#if !BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x582))
  test_spots(boost::math::concepts::real_concept(0.)); // Test real concept.
#endif
#else
   std::cout << "<note>The long double tests have been disabled on this platform "
      "either because the long double overloads of the usual math functions are "
      "not available at all, or because they are too inaccurate for these tests "
      "to pass.</note>" << std::cout;
#endif

   return 0;
} // int test_main(int, char* [])

