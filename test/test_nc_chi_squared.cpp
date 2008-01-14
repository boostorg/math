// test_nc_chi_squared.cpp

// Copyright John Maddock 2008.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/concepts/real_concept.hpp> // for real_concept
#include <boost/math/distributions/non_central_chi_squared.hpp> // for chi_squared_distribution
#include <boost/test/included/test_exec_monitor.hpp> // for test_main
#include <boost/test/floating_point_comparison.hpp> // for BOOST_CHECK_CLOSE

#include "functor.hpp"
#include "handle_test_result.hpp"
#include "test_nccs_hooks.hpp"

#include <iostream>
using std::cout;
using std::endl;
#include <limits>
using std::numeric_limits;

#define BOOST_CHECK_CLOSE_EX(a, b, prec, i) \
   {\
      unsigned int failures = boost::unit_test::results_collector.results( boost::unit_test::framework::current_test_case().p_id ).p_assertions_failed;\
      BOOST_CHECK_CLOSE(a, b, prec); \
      if(failures != boost::unit_test::results_collector.results( boost::unit_test::framework::current_test_case().p_id ).p_assertions_failed)\
      {\
         std::cerr << "Failure was at row " << i << std::endl;\
         std::cerr << std::setprecision(35); \
         std::cerr << "{ " << data[i][0] << " , " << data[i][1] << " , " << data[i][2];\
         std::cerr << " , " << data[i][3] << " , " << data[i][4] << " , " << data[i][5] << " } " << std::endl;\
      }\
   }

template <class RealType>
RealType naive_pdf(RealType v, RealType lam, RealType x)
{
   // Formula direct from 
   // http://mathworld.wolfram.com/NoncentralChi-SquaredDistribution.html
   // with no simplification:
   RealType r = -(x+lam)/2 + log(x) * (v-1)/2 + log(lam) / 2;
   r -= log(lam * x) * v/4;
   r = exp(r) / 2;
   r *= boost::math::cyl_bessel_i(v/2 - 1, sqrt(lam * x));
   return r;
}

template <class RealType>
void test_spot(
     RealType df,    // Degrees of freedom
     RealType ncp,   // non-centrality param
     RealType cs,    // Chi Square statistic
     RealType P,     // CDF
     RealType Q,     // Complement of CDF
     RealType tol)   // Test tolerance
{
   boost::math::non_central_chi_squared_distribution<RealType> dist(df, ncp);
   BOOST_CHECK_CLOSE(
      cdf(dist, cs), P, tol);
   BOOST_CHECK_CLOSE(
      pdf(dist, cs), naive_pdf(dist.degrees_of_freedom(), ncp, cs), tol);
   if((P < 0.99) && (Q < 0.99))
   {
      //
      // We can only check this if P is not too close to 1,
      // so that we can guarentee Q is free of error:
      //
      BOOST_CHECK_CLOSE(
         cdf(complement(dist, cs)), Q, tol);
      BOOST_CHECK_CLOSE(
            quantile(dist, P), cs, tol);
      BOOST_CHECK_CLOSE(
            quantile(complement(dist, Q)), cs, tol);
   }
}

template <class RealType> // Any floating-point type RealType.
void test_spots(RealType)
{
   RealType tolerance = (std::max)(
      boost::math::tools::epsilon<RealType>(),
      (RealType)boost::math::tools::epsilon<double>() * 5) * 100;
   //
   // At float precision we need to up the tolerance, since 
   // the input values are rounded off to inexact quantities
   // the results get thrown off by a noticeable amount.
   //
   if(boost::math::tools::digits<RealType>() < 50)
      tolerance *= 50;

   cout << "Tolerance = " << tolerance << "%." << endl;

   using boost::math::chi_squared_distribution;
   using  ::boost::math::chi_squared;
   using  ::boost::math::cdf;
   using  ::boost::math::pdf;
   //
   // Test against the data from Table 6 of:
   //
   // "Self-Validating Computations of Probabilities for Selected 
   // Central and Noncentral Univariate Probability Functions."
   // Morgan C. Wang; William J. Kennedy
   // Journal of the American Statistical Association, 
   // Vol. 89, No. 427. (Sep., 1994), pp. 878-887.
   //
   test_spot(
      static_cast<RealType>(1),   // degrees of freedom
      static_cast<RealType>(6),   // non centrality
      static_cast<RealType>(0.00393),   // Chi Squared statistic
      static_cast<RealType>(0.2498463724258039e-2),       // Probability of result (CDF), P
      static_cast<RealType>(1-0.2498463724258039e-2),           // Q = 1 - P
      tolerance);
   test_spot(
      static_cast<RealType>(5),   // degrees of freedom
      static_cast<RealType>(1),   // non centrality
      static_cast<RealType>(9.23636),   // Chi Squared statistic
      static_cast<RealType>(0.8272918751175548),       // Probability of result (CDF), P
      static_cast<RealType>(1-0.8272918751175548),           // Q = 1 - P
      tolerance);
   test_spot(
      static_cast<RealType>(11),   // degrees of freedom
      static_cast<RealType>(21),   // non centrality
      static_cast<RealType>(24.72497),   // Chi Squared statistic
      static_cast<RealType>(0.2539481822183126),       // Probability of result (CDF), P
      static_cast<RealType>(1-0.2539481822183126),           // Q = 1 - P
      tolerance);
   test_spot(
      static_cast<RealType>(31),   // degrees of freedom
      static_cast<RealType>(6),   // non centrality
      static_cast<RealType>(44.98534),   // Chi Squared statistic
      static_cast<RealType>(0.8125198785064969),       // Probability of result (CDF), P
      static_cast<RealType>(1-0.8125198785064969),           // Q = 1 - P
      tolerance);
   test_spot(
      static_cast<RealType>(51),   // degrees of freedom
      static_cast<RealType>(1),   // non centrality
      static_cast<RealType>(38.56038),   // Chi Squared statistic
      static_cast<RealType>(0.8519497361859118e-1),       // Probability of result (CDF), P
      static_cast<RealType>(1-0.8519497361859118e-1),           // Q = 1 - P
      tolerance);
   test_spot(
      static_cast<RealType>(100),   // degrees of freedom
      static_cast<RealType>(16),   // non centrality
      static_cast<RealType>(82.35814),   // Chi Squared statistic
      static_cast<RealType>(0.1184348822747824e-1),       // Probability of result (CDF), P
      static_cast<RealType>(1-0.1184348822747824e-1),           // Q = 1 - P
      tolerance);
   test_spot(
      static_cast<RealType>(300),   // degrees of freedom
      static_cast<RealType>(16),   // non centrality
      static_cast<RealType>(331.78852),   // Chi Squared statistic
      static_cast<RealType>(0.7355956710306709),       // Probability of result (CDF), P
      static_cast<RealType>(1-0.7355956710306709),           // Q = 1 - P
      tolerance);
   test_spot(
      static_cast<RealType>(500),   // degrees of freedom
      static_cast<RealType>(21),   // non centrality
      static_cast<RealType>(459.92612),   // Chi Squared statistic
      static_cast<RealType>(0.2797023600800060e-1),       // Probability of result (CDF), P
      static_cast<RealType>(1-0.2797023600800060e-1),           // Q = 1 - P
      tolerance);
   test_spot(
      static_cast<RealType>(1),   // degrees of freedom
      static_cast<RealType>(1),   // non centrality
      static_cast<RealType>(0.00016),   // Chi Squared statistic
      static_cast<RealType>(0.6121428929881423e-2),       // Probability of result (CDF), P
      static_cast<RealType>(1-0.6121428929881423e-2),           // Q = 1 - P
      tolerance);
   test_spot(
      static_cast<RealType>(1),   // degrees of freedom
      static_cast<RealType>(1),   // non centrality
      static_cast<RealType>(0.00393),   // Chi Squared statistic
      static_cast<RealType>(0.3033814229753780e-1),       // Probability of result (CDF), P
      static_cast<RealType>(1-0.3033814229753780e-1),           // Q = 1 - P
      tolerance);

   RealType tol2 = boost::math::tools::epsilon<RealType>() * 5 * 100; // 5 eps as a percentage
   boost::math::non_central_chi_squared_distribution<RealType> dist(static_cast<RealType>(8), static_cast<RealType>(12));
   RealType x = 7;
   using namespace std; // ADL of std names.
   // mean:
   BOOST_CHECK_CLOSE(
      mean(dist)
      , static_cast<RealType>(8+12), tol2);
   // variance:
   BOOST_CHECK_CLOSE(
      variance(dist)
      , static_cast<RealType>(64), tol2);
   // std deviation:
   BOOST_CHECK_CLOSE(
      standard_deviation(dist)
      , static_cast<RealType>(8), tol2);
   // hazard:
   BOOST_CHECK_CLOSE(
      hazard(dist, x)
      , pdf(dist, x) / cdf(complement(dist, x)), tol2);
   // cumulative hazard:
   BOOST_CHECK_CLOSE(
      chf(dist, x)
      , -log(cdf(complement(dist, x))), tol2);
   // coefficient_of_variation:
   BOOST_CHECK_CLOSE(
      coefficient_of_variation(dist)
      , standard_deviation(dist) / mean(dist), tol2);
#if 0
   // mode:
   BOOST_CHECK_CLOSE(
      mode(dist)
      , static_cast<RealType>(6), tol2);
#endif
   BOOST_CHECK_CLOSE(
      median(dist), 
      quantile(
      non_central_chi_squared_distribution<RealType>(
         static_cast<RealType>(8),
         static_cast<RealType>(12)),
      static_cast<RealType>(0.5)), static_cast<RealType>(tol2));
   // skewness:
   BOOST_CHECK_CLOSE(
      skewness(dist)
      , static_cast<RealType>(0.6875), tol2);
   // kurtosis:
   BOOST_CHECK_CLOSE(
      kurtosis(dist)
      , static_cast<RealType>(3.65625), tol2);
   // kurtosis excess:
   BOOST_CHECK_CLOSE(
      kurtosis_excess(dist)
      , static_cast<RealType>(0.65625), tol2);
} // template <class RealType>void test_spots(RealType)

template <class T>
T nccs_cdf(T df, T nc, T x)
{
   return cdf(boost::math::non_central_chi_squared_distribution<T>(df, nc), x);
}

template <class T>
T nccs_ccdf(T df, T nc, T x)
{
   return cdf(complement(boost::math::non_central_chi_squared_distribution<T>(df, nc), x));
}

template <typename T>
void do_test_nc_chi_squared(T& data, const char* type_name, const char* test)
{
   typedef typename T::value_type row_type;
   typedef typename row_type::value_type value_type;

   std::cout << "Testing: " << test << std::endl;

   value_type (*fp1)(value_type, value_type, value_type) = nccs_cdf;
   boost::math::tools::test_result<value_type> result;

   result = boost::math::tools::test(
      data,
      bind_func(fp1, 0, 1, 2),
      extract_result(3));
   handle_test_result(result, data[result.worst()], result.worst(),
      type_name, "CDF", test);

   fp1 = nccs_ccdf;
   result = boost::math::tools::test(
      data,
      bind_func(fp1, 0, 1, 2),
      extract_result(4));
   handle_test_result(result, data[result.worst()], result.worst(),
      type_name, "CCDF", test);

#ifdef TEST_OTHER
   fp1 = other::nccs_cdf;
   result = boost::math::tools::test(
      data,
      bind_func(fp1, 0, 1, 2),
      extract_result(3));
   handle_test_result(result, data[result.worst()], result.worst(),
      type_name, "other::CDF", test);
#endif

   std::cout << std::endl;

}

template <typename T>
void quantile_sanity_check(T& data, const char* type_name, const char* test)
{
   typedef typename T::value_type row_type;
   typedef typename row_type::value_type value_type;

   std::cout << "Testing: " << type_name << " quantile sanity check, with tests " << test << std::endl;

   //
   // These sanity checks test for a round trip accuracy of one half
   // of the bits in T, unless T is type float, in which case we check
   // for just one decimal digit.  The problem here is the sensitivity
   // of the functions, not their accuracy.  This test data was generated
   // for the forward functions, which means that when it is used as
   // the input to the inverses then it is necessarily inexact.  This rounding
   // of the input is what makes the data unsuitable for use as an accuracy check,
   // and also demonstrates that you can't in general round-trip these functions.
   // It is however a useful sanity check.
   //
   value_type precision = static_cast<value_type>(ldexp(1.0, 1-boost::math::policies::digits<value_type, boost::math::policies::policy<> >()/2)) * 100;
   if(boost::math::policies::digits<value_type, boost::math::policies::policy<> >() < 50)
      precision = 1;   // 1% or two decimal digits, all we can hope for when the input is truncated to float

   for(unsigned i = 0; i < data.size(); ++i)
   {
      if(data[i][3] == 0)
      {
         BOOST_CHECK(0 == quantile(boost::math::non_central_chi_squared_distribution<value_type>(data[i][0], data[i][1]), data[i][3]));
      }
      else if(data[i][3] < 0.9999f)
      {
         value_type p = quantile(boost::math::non_central_chi_squared_distribution<value_type>(data[i][0], data[i][1]), data[i][3]);
         value_type pt = data[i][2];
         BOOST_CHECK_CLOSE_EX(pt, p, precision, i);
      }
      if(data[i][4] == 0)
      {
         BOOST_CHECK(0 == quantile(complement(boost::math::non_central_chi_squared_distribution<value_type>(data[i][0], data[i][1]), data[i][3])));
      }
      else if(data[i][4] < 0.9999f)
      {
         value_type p = quantile(complement(boost::math::non_central_chi_squared_distribution<value_type>(data[i][0], data[i][1]), data[i][4]));
         value_type pt = data[i][2];
         BOOST_CHECK_CLOSE_EX(pt, p, precision, i);
      }
   }
}

template <typename T>
void test_accuracy(T, const char* type_name)
{
#if 0
#include "nccs.ipp"
    do_test_nc_chi_squared(nccs, type_name, "Non Central Chi Squared, medium parameters");
    quantile_sanity_check(nccs, type_name, "Non Central Chi Squared, medium parameters");
#endif
#include "nccs_big.ipp"
    do_test_nc_chi_squared(nccs_big, type_name, "Non Central Chi Squared, large parameters");
    quantile_sanity_check(nccs_big, type_name, "Non Central Chi Squared, large parameters");
}

int test_main(int, char* [])
{
   BOOST_MATH_CONTROL_FP;
   // Basic sanity-check spot values.

   // (Parameter value, arbitrarily zero, only communicates the floating point type).
   test_spots(0.0F); // Test float.
   test_spots(0.0); // Test double.
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
   test_spots(0.0L); // Test long double.
#if !BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x582))
   test_spots(boost::math::concepts::real_concept(0.)); // Test real concept.
#endif
#endif
   test_accuracy(0.0F, "float"); // Test float.
   test_accuracy(0.0, "double"); // Test double.
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
   test_accuracy(0.0L, "long double"); // Test long double.
#if !BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x582))
   test_accuracy(boost::math::concepts::real_concept(0.), "real_concept"); // Test real concept.
#endif
#endif
   return 0;
} // int test_main(int, char* [])

/*


*/



