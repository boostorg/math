// test_binomial.cpp

// Copyright Paul A. Bristow 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// Basic sanity test for Binomial Cumulative Distribution Function.

#define BOOST_MATH_THROW_ON_DOMAIN_ERROR

#ifdef _MSC_VER
#  pragma warning(disable: 4127) // conditional expression is constant.
#  pragma warning(disable: 4100) // unreferenced formal parameter.
#  pragma warning(disable: 4512) // assignment operator could not be generated.
#endif

#include <boost/math/distributions/binomial.hpp> // for binomial_distribution
using boost::math::binomial_distribution;

#include <boost/math/concepts/real_concept.hpp> // for real_concept
using ::boost::math::concepts::real_concept;

#include <boost/test/included/test_exec_monitor.hpp> // for test_main
#include <boost/test/floating_point_comparison.hpp> // for BOOST_CHECK_CLOSE

#include <iostream>
using std::cout;
using std::endl;
#include <limits>
using std::numeric_limits;

template <class RealType>
void test_spot(
     RealType N,    // Number of trials
     RealType k,    // Number of successes
     RealType p,    // Probability of success
     RealType P,    // CDF
     RealType Q,    // Complement of CDF
     RealType tol)  // Test tolerance
{
   boost::math::binomial_distribution<RealType> bn(N, p);
   BOOST_CHECK_CLOSE(
      cdf(bn, k), P, tol);
   if((P < 0.99) && (Q < 0.99))
   {
      //
      // We can only check this if P is not too close to 1,
      // so that we can guarentee Q is free of error:
      //
      BOOST_CHECK_CLOSE(
         cdf(complement(bn, k)), Q, tol);
      if(k != 0)
      {
         BOOST_CHECK_CLOSE(
            quantile(bn, P), k, tol);
      }
      else
      {
         // Just check quantile is very small:
         BOOST_CHECK(quantile(bn, P) < boost::math::tools::epsilon<RealType>() * 10);
      }
      if(k != 0)
      {
         BOOST_CHECK_CLOSE(
            quantile(complement(bn, Q)), k, tol);
      }
      else
      {
         // Just check quantile is very small:
         BOOST_CHECK(quantile(complement(bn, Q)) < boost::math::tools::epsilon<RealType>() * 10);
      }
      // estimate success ratio:
      BOOST_CHECK_CLOSE(
         binomial_distribution<RealType>::estimate_lower_bound_on_p(
            N, k, Q), 
         p, tol);
      BOOST_CHECK_CLOSE(
         binomial_distribution<RealType>::estimate_upper_bound_on_p(
            N, k, P), 
         p, tol);
      
      if(Q < P)
      {
         BOOST_CHECK(
            binomial_distribution<RealType>::estimate_lower_bound_on_p(
               N, k, Q)
               < 
            binomial_distribution<RealType>::estimate_upper_bound_on_p(
               N, k, Q)
               );
      }
      else
      {
         BOOST_CHECK(
            binomial_distribution<RealType>::estimate_lower_bound_on_p(
               N, k, P)
               < 
            binomial_distribution<RealType>::estimate_upper_bound_on_p(
               N, k, P)
               );
      }
      //
      // estimate sample size:
      //
      BOOST_CHECK_CLOSE(
         binomial_distribution<RealType>::estimate_number_of_trials(
            k, p, P), 
         N, tol);
      BOOST_CHECK_CLOSE(
         binomial_distribution<RealType>::estimate_number_of_trials(
            boost::math::complement(k, p, Q)), 
         N, tol);
   }

   // Double check consistency of CDF and PDF by computing 
   // the finite sum:
   RealType sum = 0;
   for(unsigned i = 0; i <= k; ++i)
      sum += pdf(bn, RealType(i));
   BOOST_CHECK_CLOSE(
      sum, P, tol);
   // And complement as well:
   sum = 0;
   for(RealType i = N; i > k; i -= 1)
      sum += pdf(bn, i);
   if(P < 0.99)
   {
      BOOST_CHECK_CLOSE(
         sum, Q, tol);
   }
   else
   {
      // Not enough information content in P for Q to be meaningful
      RealType tol = (std::max)(2 * Q, boost::math::tools::epsilon<RealType>());
      BOOST_CHECK(sum < tol);
   }
}

template <class RealType> // Any floating-point type RealType.
void test_spots(RealType)
{
  // Basic sanity checks, test data is to double precision only
  // so set tolerance to 100eps expressed as a persent, or
  // 100eps of type double expressed as a persent, whichever
  // is the larger.

  RealType tolerance = (std::max)
      (boost::math::tools::epsilon<RealType>(),
      static_cast<RealType>(std::numeric_limits<double>::epsilon()));
   tolerance *= 100 * 1000;

  cout << "Tolerance = " << tolerance << "%." << endl;

  // Sources of spot test values:

  // MathCAD defines pbinom(k, n, p)
  // returns pr(X ,=k) when random variable X has the binomial distribution with parameters n and p.
  // 0 <= k ,= n
  // 0 <= p <= 1
  // P = pbinom(30, 500, 0.05) = 0.869147702104609

  using boost::math::binomial_distribution;
  using  ::boost::math::binomial;
  using  ::boost::math::cdf;
  using  ::boost::math::pdf;

  // Test binomial using cdf spot values from MathCAD.
  // These test quantiles and complements as well.
  test_spot(
     static_cast<RealType>(500),                     // Sample size, N
     static_cast<RealType>(30),                      // Number of successes, k
     static_cast<RealType>(0.05),                    // Probability of success, p
     static_cast<RealType>(0.869147702104609),       // Probability of result (CDF), P
     static_cast<RealType>(1 - 0.869147702104609),   // Q = 1 - P
     tolerance);

  test_spot(
     static_cast<RealType>(500),                     // Sample size, N
     static_cast<RealType>(250),                     // Number of successes, k
     static_cast<RealType>(0.05),                    // Probability of success, p
     static_cast<RealType>(1),                       // Probability of result (CDF), P
     static_cast<RealType>(0),   // Q = 1 - P
     tolerance);

  test_spot(
     static_cast<RealType>(500),                     // Sample size, N
     static_cast<RealType>(470),                     // Number of successes, k
     static_cast<RealType>(0.95),                    // Probability of success, p
     static_cast<RealType>(0.176470742656766),       // Probability of result (CDF), P
     static_cast<RealType>(1 - 0.176470742656766),   // Q = 1 - P
     tolerance * 10);                                // Note higher tolerance on this test!

  test_spot(
     static_cast<RealType>(500),                       // Sample size, N
     static_cast<RealType>(400),                       // Number of successes, k
     static_cast<RealType>(0.05),                      // Probability of success, p
     static_cast<RealType>(1),                         // Probability of result (CDF), P
     static_cast<RealType>(0),                         // Q = 1 - P
     tolerance);

  test_spot(
     static_cast<RealType>(500),                       // Sample size, N
     static_cast<RealType>(400),                       // Number of successes, k
     static_cast<RealType>(0.9),                       // Probability of success, p
     static_cast<RealType>(1.80180425681923E-11),      // Probability of result (CDF), P
     static_cast<RealType>(1 - 1.80180425681923E-11),  // Q = 1 - P
     tolerance);

  test_spot(
     static_cast<RealType>(500),                       // Sample size, N
     static_cast<RealType>(5),                         // Number of successes, k
     static_cast<RealType>(0.05),                      // Probability of success, p
     static_cast<RealType>(9.181808267643E-7),         // Probability of result (CDF), P
     static_cast<RealType>(1 - 9.181808267643E-7),     // Q = 1 - P
     tolerance);

  test_spot(
     static_cast<RealType>(2),                       // Sample size, N
     static_cast<RealType>(1),                       // Number of successes, k
     static_cast<RealType>(0.5),                     // Probability of success, p
     static_cast<RealType>(0.75),                    // Probability of result (CDF), P
     static_cast<RealType>(0.25),                    // Q = 1 - P
     tolerance);
  
  test_spot(
     static_cast<RealType>(8),                       // Sample size, N
     static_cast<RealType>(3),                       // Number of successes, k
     static_cast<RealType>(0.25),                    // Probability of success, p
     static_cast<RealType>(0.8861846923828125),      // Probability of result (CDF), P
     static_cast<RealType>(1 - 0.8861846923828125),  // Q = 1 - P
     tolerance);
  
  test_spot(
     static_cast<RealType>(8),                       // Sample size, N
     static_cast<RealType>(0),                       // Number of successes, k
     static_cast<RealType>(0.25),                    // Probability of success, p
     static_cast<RealType>(0.1001129150390625),      // Probability of result (CDF), P
     static_cast<RealType>(1 - 0.1001129150390625),  // Q = 1 - P
     tolerance);
  
  test_spot(
     static_cast<RealType>(8),                       // Sample size, N
     static_cast<RealType>(1),                       // Number of successes, k
     static_cast<RealType>(0.25),                    // Probability of success, p
     static_cast<RealType>(0.36708068847656244),     // Probability of result (CDF), P
     static_cast<RealType>(1 - 0.36708068847656244), // Q = 1 - P
     tolerance);
  
  test_spot(
     static_cast<RealType>(8),                       // Sample size, N
     static_cast<RealType>(4),                       // Number of successes, k
     static_cast<RealType>(0.25),                    // Probability of success, p
     static_cast<RealType>(0.9727020263671875),      // Probability of result (CDF), P
     static_cast<RealType>(1 - 0.9727020263671875),  // Q = 1 - P
     tolerance);
  
  test_spot(
     static_cast<RealType>(8),                       // Sample size, N
     static_cast<RealType>(7),                       // Number of successes, k
     static_cast<RealType>(0.25),                    // Probability of success, p
     static_cast<RealType>(0.9999847412109375),      // Probability of result (CDF), P
     static_cast<RealType>(1 - 0.9999847412109375),  // Q = 1 - P
     tolerance);
  
  // Tests on PDF follow:
  BOOST_CHECK_CLOSE(
     pdf(binomial_distribution<RealType>(static_cast<RealType>(20), static_cast<RealType>(0.75)),
     static_cast<RealType>(10)),  // k.
     static_cast<RealType>(0.00992227527967770583927631378173), // 0.00992227527967770583927631378173
     tolerance);

  BOOST_CHECK_CLOSE(
    pdf(binomial_distribution<RealType>(static_cast<RealType>(20), static_cast<RealType>(0.5)), 
    static_cast<RealType>(10)),  // k.
    static_cast<RealType>(0.17619705200195312500000000000000000000), // get k=10 0.049611376398388612 p = 0.25
    tolerance);
  
  // Binomial pdf Test values from
  // http://www.adsciengineering.com/bpdcalc/index.php  for example
  // http://www.adsciengineering.com/bpdcalc/index.php?n=20&p=0.25&start=0&stop=20&Submit=Generate
  // Appears to use at least 80-bit long double for 32 decimal digits accuracy,
  // (if trailings zero then are exact values?)
  // so useful for testing 64-bit double accuracy.
  // P = 0.25, n = 20, k = 0 to 20

  //0	C(20,0) * 0.25^0 * 0.75^20	0.00317121193893399322405457496643
  //1	C(20,1) * 0.25^1 * 0.75^19	0.02114141292622662149369716644287
  //2	C(20,2) * 0.25^2 * 0.75^18	0.06694780759971763473004102706909
  //3	C(20,3) * 0.25^3 * 0.75^17	0.13389561519943526946008205413818
  //4	C(20,4) * 0.25^4 * 0.75^16	0.18968545486586663173511624336242
  //5	C(20,5) * 0.25^5 * 0.75^15	0.20233115185692440718412399291992
  //6	C(20,6) * 0.25^6 * 0.75^14	0.16860929321410367265343666076660
  //7	C(20,7) * 0.25^7 * 0.75^13	0.11240619547606911510229110717773
  //8	C(20,8) * 0.25^8 * 0.75^12	0.06088668921620410401374101638793
  //9	C(20,9) * 0.25^9 * 0.75^11	0.02706075076275737956166267395019
  //10	C(20,10) * 0.25^10 * 0.75^10	0.00992227527967770583927631378173
  //11	C(20,11) * 0.25^11 * 0.75^9	0.00300675008475081995129585266113
  //12	C(20,12) * 0.25^12 * 0.75^8	0.00075168752118770498782396316528
  //13	C(20,13) * 0.25^13 * 0.75^7	0.00015419231203850358724594116210
  //14	C(20,14) * 0.25^14 * 0.75^6	0.00002569871867308393120765686035
  //15	C(20,15) * 0.25^15 * 0.75^5	0.00000342649582307785749435424804
  //16	C(20,16) * 0.25^16 * 0.75^4	0.00000035692664823727682232856750
  //17	C(20,17) * 0.25^17 * 0.75^3	0.00000002799424692057073116302490
  //18	C(20,18) * 0.25^18 * 0.75^2	0.00000000155523594003170728683471
  //19	C(20,19) * 0.25^19 * 0.75^1	0.00000000005456968210637569427490
  //20	C(20,20) * 0.25^20 * 0.75^0	0.00000000000090949470177292823791


    BOOST_CHECK_CLOSE(
    pdf(binomial_distribution<RealType>(static_cast<RealType>(20), static_cast<RealType>(0.25)), 
    static_cast<RealType>(10)),  // k.
    static_cast<RealType>(0.00992227527967770583927631378173), // k=10  p = 0.25
    tolerance); 
  
    BOOST_CHECK_CLOSE( // k = 0 use different formula - only exp so more accurate.
    pdf(binomial_distribution<RealType>(static_cast<RealType>(20), static_cast<RealType>(0.25)), 
    static_cast<RealType>(0)),  // k.
    static_cast<RealType>(0.00317121193893399322405457496643), // k=0  p = 0.25
    tolerance); 

    BOOST_CHECK_CLOSE( // k = 20 use different formula - only exp so more accurate.
    pdf(binomial_distribution<RealType>(static_cast<RealType>(20), static_cast<RealType>(0.25)), 
    static_cast<RealType>(20)),  // k == n.
    static_cast<RealType>(0.00000000000090949470177292823791), // k=20  p = 0.25
    tolerance); 

    BOOST_CHECK_CLOSE( // k = 1.
    pdf(binomial_distribution<RealType>(static_cast<RealType>(20), static_cast<RealType>(0.25)), 
    static_cast<RealType>(1)),  // k.
    static_cast<RealType>(0.02114141292622662149369716644287), // k=1  p = 0.25
    tolerance); 

    // Some exact (probably) values.
    BOOST_CHECK_CLOSE( 
    pdf(binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)), 
    static_cast<RealType>(0)),  // k.
    static_cast<RealType>(0.10011291503906250000000000000000), // k=0  p = 0.25
    tolerance); 

    BOOST_CHECK_CLOSE( // k = 1.
    pdf(binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)), 
    static_cast<RealType>(1)),  // k.
    static_cast<RealType>(0.26696777343750000000000000000000), // k=1  p = 0.25
    tolerance); 

    BOOST_CHECK_CLOSE( // k = 2.
    pdf(binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)), 
    static_cast<RealType>(2)),  // k.
    static_cast<RealType>(0.31146240234375000000000000000000), // k=2  p = 0.25
    tolerance); 

    BOOST_CHECK_CLOSE( // k = 3.
    pdf(binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)), 
    static_cast<RealType>(3)),  // k.
    static_cast<RealType>(0.20764160156250000000000000000000), // k=3  p = 0.25
    tolerance); 

    BOOST_CHECK_CLOSE( // k = 7.
    pdf(binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)), 
    static_cast<RealType>(7)),  // k.
    static_cast<RealType>(0.00036621093750000000000000000000), // k=7  p = 0.25
    tolerance); 

    BOOST_CHECK_CLOSE( // k = 8.
    pdf(binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)), 
    static_cast<RealType>(8)),  // k = n.
    static_cast<RealType>(0.00001525878906250000000000000000), // k=8  p = 0.25
    tolerance); 

    RealType tol2 = boost::math::tools::epsilon<RealType>() * 5;
    binomial_distribution<RealType> dist(static_cast<RealType>(8), static_cast<RealType>(0.25));
    RealType x = static_cast<RealType>(0.125);
    using namespace std; // ADL of std names.
    // mean:
    BOOST_CHECK_CLOSE(
       mean(dist)
       , static_cast<RealType>(8 * 0.25), tol2);
    // variance:
    BOOST_CHECK_CLOSE(
       variance(dist)
       , static_cast<RealType>(8 * 0.25 * 0.75), tol2);
    // std deviation:
    BOOST_CHECK_CLOSE(
       standard_deviation(dist)
       , static_cast<RealType>(sqrt(8 * 0.25 * 0.75)), tol2);
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

    // special cases for PDF:
    BOOST_CHECK_EQUAL(
       pdf(
          binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0)),
          static_cast<RealType>(0)), static_cast<RealType>(1)
       );
    BOOST_CHECK_EQUAL(
       pdf(
          binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0)),
          static_cast<RealType>(0.0001)), static_cast<RealType>(0)
       );
    BOOST_CHECK_EQUAL(
       pdf(
          binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(1)),
          static_cast<RealType>(0.001)), static_cast<RealType>(0)
       );
    BOOST_CHECK_EQUAL(
       pdf(
          binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(1)),
          static_cast<RealType>(8)), static_cast<RealType>(1)
       );
    BOOST_CHECK_EQUAL(
       pdf(
          binomial_distribution<RealType>(static_cast<RealType>(0), static_cast<RealType>(0.25)),
          static_cast<RealType>(0)), static_cast<RealType>(1)
       );
    BOOST_CHECK_THROW(
       pdf(
          binomial_distribution<RealType>(static_cast<RealType>(-1), static_cast<RealType>(0.25)),
          static_cast<RealType>(0)), std::domain_error
       );
    BOOST_CHECK_THROW(
       pdf(
          binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(-0.25)),
          static_cast<RealType>(0)), std::domain_error
       );
    BOOST_CHECK_THROW(
       pdf(
          binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(1.25)),
          static_cast<RealType>(0)), std::domain_error
       );
    BOOST_CHECK_THROW(
       pdf(
          binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)),
          static_cast<RealType>(-1)), std::domain_error
       );
    BOOST_CHECK_THROW(
       pdf(
          binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)),
          static_cast<RealType>(9)), std::domain_error
       );
    BOOST_CHECK_THROW(
       cdf(
          binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)),
          static_cast<RealType>(-1)), std::domain_error
       );
    BOOST_CHECK_THROW(
       cdf(
          binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)),
          static_cast<RealType>(9)), std::domain_error
       );
    BOOST_CHECK_THROW(
       cdf(
          binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(-0.25)),
          static_cast<RealType>(0)), std::domain_error
       );
    BOOST_CHECK_THROW(
       cdf(
          binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(1.25)),
          static_cast<RealType>(0)), std::domain_error
       );
    BOOST_CHECK_THROW(
       quantile(
          binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(-0.25)),
          static_cast<RealType>(0)), std::domain_error
       );
    BOOST_CHECK_THROW(
       quantile(
          binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(1.25)),
          static_cast<RealType>(0)), std::domain_error
       );
    BOOST_CHECK_EQUAL(
       cdf(
          binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0.25)),
          static_cast<RealType>(8)), static_cast<RealType>(1)
       );
    BOOST_CHECK_EQUAL(
       cdf(
          binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(0)),
          static_cast<RealType>(7)), static_cast<RealType>(1)
       );
    BOOST_CHECK_EQUAL(
       cdf(
          binomial_distribution<RealType>(static_cast<RealType>(8), static_cast<RealType>(1)),
          static_cast<RealType>(7)), static_cast<RealType>(0)
       );

   if(boost::is_floating_point<RealType>::value)
   {
     // This is an extreme test, without Lanczos approximation support
     // we get pretty inaccurate results from real_concept (about 5 sig digits only).
     // So only test for accuracy with "real" floating point types:
     BOOST_CHECK_CLOSE(::boost::math::pdf( // Really big number of trials - max()
       binomial_distribution<RealType>(static_cast<RealType>((numeric_limits<int>::max)()), static_cast<RealType>(0.5)),
       static_cast<RealType>((numeric_limits<int>::max)()/2)),  // half way to max
       static_cast<RealType>(1.7217698023561394e-005), // probability.
       tolerance);
   }

  {
    // This is a visual sanity check that everything is OK:
    binomial_distribution<RealType> my8dist(8., 0.25); // Note: double values (matching the distribution definition) avoid the need for any casting.
    cout << "mean(my8dist) = " << boost::math::mean(my8dist) << endl; // mean(my8dist) = 2
    cout << "my8dist.trials() = " << my8dist.trials()  << endl; // my8dist.trials() = 8
    cout << "my8dist.success_fraction() = " << my8dist.success_fraction()  << endl; // my8dist.success_fraction() = 0.25
    BOOST_CHECK_CLOSE(my8dist.trials(), static_cast<RealType>(8), tol2);
    BOOST_CHECK_CLOSE(my8dist.success_fraction(), static_cast<RealType>(0.25), tol2);

   {
      int n = static_cast<int>(boost::math::tools::real_cast<double>(my8dist.trials()));
      RealType sumcdf = 0.;
      for (int k = 0; k <= n; k++)
      {
        cout << k << ' ' << pdf(my8dist, static_cast<RealType>(k));
        sumcdf += pdf(my8dist, static_cast<RealType>(k));
        cout  << ' '  << sumcdf;
        cout << ' ' << cdf(my8dist, static_cast<RealType>(k));
        cout << ' ' << sumcdf - cdf(my8dist, static_cast<RealType>(k)) << endl;
      } // for k
    }
    // n = 8, p =0.25
    //k         pdf              cdf
    //0 0.1001129150390625 0.1001129150390625 
    //1 0.26696777343749994 0.36708068847656244 
    //2 0.31146240234375017 0.67854309082031261 
    //3 0.20764160156249989 0.8861846923828125 
    //4 0.086517333984375 0.9727020263671875 
    //5 0.023071289062499997 0.9957733154296875 
    //6 0.0038452148437500009 0.9996185302734375 
    //7 0.00036621093749999984 0.9999847412109375 
    //8 1.52587890625e-005 1 1 0
  }

} // template <class RealType>void test_spots(RealType)

int test_main(int, char* [])
{
  // Basic sanity-check spot values.
#ifdef BOOST_MATH_THROW_ON_DOMAIN_ERROR
  cout << "BOOST_MATH_THROW_ON_DOMAIN_ERROR" << " is defined to throw on domain error." << endl;
#else
  cout << "BOOST_MATH_THROW_ON_DOMAIN_ERROR" << " is NOT defined, so NO throw on domain error." << endl;
#endif

  // (Parameter value, arbitrarily zero, only communicates the floating point type).
  test_spots(0.0F); // Test float.
  test_spots(0.0); // Test double.
  test_spots(0.0L); // Test long double.
  test_spots(boost::math::concepts::real_concept(0.)); // Test real concept.

  return 0;
} // int test_main(int, char* [])

/*
------ Build started: Project: test_binomial, Configuration: Debug Win32 ------
Compiling...
test_binomial.cpp
Linking...
Autorun "i:\boost-06-05-03-1300\libs\math\test\Math_test\debug\test_binomial.exe"
Running 1 test case...
BOOST_MATH_THROW_ON_DOMAIN_ERROR is defined to throw on domain error.
Tolerance = 1e-011%.
mydist.trials() = 20
mydist.success_fraction() = 0.25
pdf(mydist, 10.) = 0.0099222752796776954
0 0.0031712119389339932
1 0.021141412926226611
2 0.066947807599717579
3 0.13389561519943508
4 0.18968545486586622
5 0.20233115185692488
6 0.1686092932141042
7 0.11240619547606928
8 0.060886689216204264
9 0.027060750762757328
10 0.0099222752796776954
11 0.0030067500847508143
12 0.00075168752118770694
13 0.00015419231203850383
14 2.5698718673084013e-005
15 3.4264958230778655e-006
16 3.5692664823727603e-007
17 2.7994246920570691e-008
18 1.5552359400317058e-009
19 5.4569682106375668e-011
20 9.0949470177292824e-013
0 0.100113 0.100113 0.100113 0
1 0.266968 0.367081 0.367081 -2.77556e-016
2 0.311462 0.678543 0.678543 -8.88178e-016
3 0.207642 0.886185 0.886185 -1.11022e-015
4 0.0865173 0.972702 0.972702 -1.22125e-015
5 0.0230713 0.995773 0.995773 -1.22125e-015
6 0.00384521 0.999619 0.999619 -1.22125e-015
7 0.000366211 0.999985 0.999985 -1.22125e-015
8 1.52588e-005 1 1 -1.22125e-015
1
1
*** No errors detected
Build Time 0:05
Build log was saved at "file://i:\boost-06-05-03-1300\libs\math\test\Math_test\test_binomial\Debug\BuildLog.htm"
test_binomial - 0 error(s), 0 warning(s)
========== Build: 1 succeeded, 0 failed, 0 up-to-date, 0 skipped ==========

*/