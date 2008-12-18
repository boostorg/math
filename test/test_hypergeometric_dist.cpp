// Copyright John Maddock 
// Copyright Paul A. Bristow 
// Copyright Gautam Sewani

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_MATH_OVERFLOW_ERROR_POLICY throw_on_error
#include <boost/math/concepts/real_concept.hpp> // for real_concept
#include <boost/math/distributions/hypergeometric.hpp>

#include <boost/test/included/test_exec_monitor.hpp> // Boost.Test
#include <boost/test/floating_point_comparison.hpp>

#include <iostream>
   using std::cout;
   using std::endl;
   using std::setprecision;

#include <boost/array.hpp>
#include "functor.hpp"
#include "handle_test_result.hpp"

#define BOOST_CHECK_EX(a) \
   {\
      unsigned int failures = boost::unit_test::results_collector.results( boost::unit_test::framework::current_test_case().p_id ).p_assertions_failed;\
      BOOST_CHECK(a); \
      if(failures != boost::unit_test::results_collector.results( boost::unit_test::framework::current_test_case().p_id ).p_assertions_failed)\
      {\
         std::cerr << "Failure was with data ";\
         std::cerr << std::setprecision(35); \
         std::cerr << "x = " << x << ", r = " << r << ", n = " << n\
         << ", N = " << N << ", p = " << cp << ", q = " << ccp << std::endl;\
      }\
   }

   void expected_results()
{
   //
   // Define the max and mean errors expected for
   // various compilers and platforms.
   //
   const char* largest_type;
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
   if(boost::math::policies::digits<double, boost::math::policies::policy<> >() == boost::math::policies::digits<long double, boost::math::policies::policy<> >())
   {
      largest_type = "(long\\s+)?double|real_concept";
   }
   else
   {
      largest_type = "long double|real_concept";
   }
#else
   largest_type = "(long\\s+)?double";
#endif

   add_expected_result(
      ".*",                          // compiler
      ".*",                          // stdlib
      ".*",                          // platform
      largest_type,                  // test type(s)
      ".*",                          // test data group
      ".*", 50, 20);                 // test function
}

template <class T>
inline unsigned make_unsigned(T x)
{
   return static_cast<unsigned>(x);
}
template<>
inline unsigned make_unsigned(boost::math::concepts::real_concept x)
{
   return static_cast<unsigned>(x.value());
}

template <class T>
T pdf_tester(T r, T n, T N, T x)
{
   boost::math::hypergeometric_distribution<T> d(make_unsigned(r), make_unsigned(n), make_unsigned(N));
   return pdf(d, x);
}

template <class T>
void do_test_hypergeometric(const T& data, const char* type_name, const char* test_name)
{
   typedef typename T::value_type row_type;
   typedef typename row_type::value_type value_type;

   typedef value_type (*pg)(value_type, value_type, value_type, value_type);
#if defined(BOOST_MATH_NO_DEDUCED_FUNCTION_POINTERS)
   pg funcp = pdf_tester<value_type>;
#else
   pg funcp = pdf_tester;
#endif

   boost::math::tools::test_result<value_type> result;

   std::cout << "Testing " << test_name << " with type " << type_name
      << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

   //
   // test hypergeometric against data:
   //
   result = boost::math::tools::test(
      data, 
      bind_func(funcp, 0, 1, 2, 3), 
      extract_result(4));
   handle_test_result(result, data[result.worst()], result.worst(), type_name, "hypergeometric PDF", test_name);
   std::cout << std::endl;
}

template <class RealType>
void test_spot(unsigned x, unsigned n, unsigned r, unsigned N, 
               RealType p, RealType cp, RealType ccp, RealType tol)
{
   boost::math::hypergeometric_distribution<RealType> d(r, n, N);

   std::pair<unsigned, unsigned> extent = range(d);

   BOOST_CHECK_CLOSE(pdf(d, x), p, tol);
   BOOST_CHECK_CLOSE(cdf(d, x), cp, tol);
   BOOST_CHECK_CLOSE(cdf(complement(d, x)), ccp, tol);
   // Again with real-value arguments:
   BOOST_CHECK_CLOSE(pdf(d, static_cast<RealType>(x)), p, tol);
   BOOST_CHECK_CLOSE(cdf(d, static_cast<RealType>(x)), cp, tol);
   BOOST_CHECK_CLOSE(cdf(complement(d, static_cast<RealType>(x))), ccp, tol);
   // quantiles:
   if(boost::math::tools::digits<RealType>() > 50)
   {
      using namespace boost::math::policies;
      boost::math::hypergeometric_distribution<RealType, 
         policy<discrete_quantile<integer_round_up> > > du(r, n, N);
      BOOST_CHECK_EX(quantile(du, cp) >= x);
      BOOST_CHECK_EX(quantile(complement(du, ccp)) >= x);

      boost::math::hypergeometric_distribution<RealType, 
         policy<discrete_quantile<integer_round_down> > > dl(r, n, N);
      BOOST_CHECK_EX(quantile(dl, cp) <= x);
      BOOST_CHECK_EX(quantile(complement(dl, ccp)) <= x);
   }
}

template <class RealType>
void test_spots(RealType /*T*/, const char* type_name)
{
   // Basic sanity checks.
   // 50 eps as a percentage, up to a maximum of double precision
   // Test data taken from Mathematica 6
#define T RealType
#include "hypergeometric_test_data.ipp"
   do_test_hypergeometric(hypergeometric_test_data, type_name, "Mathematica data");

   RealType tolerance = (std::max)(
      static_cast<RealType>(2e-16L),    // limit of test data
      boost::math::tools::epsilon<RealType>());
   cout<<"Absolute tolerance:"<<tolerance<<endl;
   
   tolerance *= 50 * 100; // 50eps as a persentage
   cout << "Tolerance for type " << typeid(RealType).name()  << " is " << tolerance << " %" << endl;

   test_spot(20, 200, 50, 500, static_cast<T>(0.120748236361163), static_cast<T>(0.563532430195156), static_cast<T>(1 - 0.563532430195156), tolerance);
   test_spot(53, 452, 64, 500, static_cast<T>(0.0184749573044286), static_cast<T>(0.0299118078796907), static_cast<T>(1 - 0.0299118078796907), tolerance);
   test_spot(32, 1287, 128, 5000, static_cast<T>(0.0807012167418264), static_cast<T>(0.469768774237742), static_cast<T>(1 - 0.469768774237742), tolerance);
   test_spot(1, 13, 4, 26, static_cast<T>(0.248695652173913), static_cast<T>(0.296521739130435), static_cast<T>(1 - 0.296521739130435), tolerance);
   test_spot(2, 13, 4, 26, static_cast<T>(0.40695652173913), static_cast<T>(0.703478260869565), static_cast<T>(1 - 0.703478260869565), tolerance);
   test_spot(3, 13, 4, 26, static_cast<T>(0.248695652173913), static_cast<T>(0.952173913043478), static_cast<T>(1 - 0.952173913043478), tolerance);

   boost::math::hypergeometric_distribution<RealType> d(50, 200, 500);
   BOOST_CHECK_EQUAL(range(d).first, 0u);
   BOOST_CHECK_EQUAL(range(d).second, 50u);
   BOOST_CHECK_CLOSE(mean(d), static_cast<RealType>(20), tolerance);
   BOOST_CHECK_CLOSE(mode(d), static_cast<RealType>(20), tolerance);
   BOOST_CHECK_CLOSE(variance(d), static_cast<RealType>(10.821643286573146292585170340681L), tolerance);
   BOOST_CHECK_CLOSE(skewness(d), static_cast<RealType>(0.048833071022952084732902910189366L), tolerance);
   BOOST_CHECK_CLOSE(kurtosis_excess(d), static_cast<RealType>(2.5155486690782804816404001878293L), tolerance);
   BOOST_CHECK_CLOSE(kurtosis(d), kurtosis_excess(d) + 3, tolerance);
}


int test_main(int, char* [])
{
   expected_results();
   // Basic sanity-check spot values.
   // (Parameter value, arbitrarily zero, only communicates the floating point type).
   test_spots(0.0F, "float"); // Test float. OK at decdigits = 0 tolerance = 0.0001 %
   test_spots(0.0, "double"); // Test double. OK at decdigits 7, tolerance = 1e07 %
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
   test_spots(0.0L, "long double"); // Test long double.
   test_spots(boost::math::concepts::real_concept(0), "real_concept"); // Test real_concept.
#else
   std::cout << "<note>The long double tests have been disabled on this platform "
      "either because the long double overloads of the usual math functions are "
      "not available at all, or because they are too inaccurate for these tests "
      "to pass.</note>" << std::cout;
#endif

   return 0;
} // int test_main(int, char* [])

