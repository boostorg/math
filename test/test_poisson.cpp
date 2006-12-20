// test_poisson.cpp

// Copyright Paul A. Bristow 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// Basic sanity test for Poisson Cumulative Distribution Function.

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

#include <boost/test/included/test_exec_monitor.hpp> // Boost.Test
#include <boost/test/floating_point_comparison.hpp>

#include <boost/math/distributions/poisson.hpp>
	 using boost::math::poisson_distribution;
#include <boost/math/concepts/real_concept.hpp> // for real_concept
#include <boost/math/tools/test.hpp> // for real_concept

#include <boost/math/special_functions/gamma.hpp> // for (incomplete) gamma.
//   using boost::math::qamma_Q;

#include <iostream>
	using std::cout;
	using std::endl;
	using std::setprecision;
	using std::showpoint;
	using std::ios;
#include <limits>
  using std::numeric_limits;

template <class RealType> // Any floating-point type RealType.
void test_spots(RealType)
{
  // Basic sanity checks, tolerance is about numeric_limits<RealType>::digits10 decimal places,
	// guaranteed for type RealType, eg 6 for float, 15 for double,
	// expressed as a percentage (so -2) for BOOST_CHECK_CLOSE,

	int decdigits = numeric_limits<RealType>::digits10;
  // May eb >15 for 80 and 128-bit FP typtes.
  if (decdigits <= 0)
  { // decdigits is not defined, for example real concept,
    // so assume precision of most test data is double (for example, MathCAD).
     decdigits = numeric_limits<double>::digits10; // == 15 for 64-bit
  }
  if (decdigits > 15 ) // numeric_limits<double>::digits10)
  { // 15 is the accuracy of the MathCAD test data.
    decdigits = 15; // numeric_limits<double>::digits10;
  }

	decdigits -= 1; // Perhaps allow some decimal digit(s) margin of numerical error.
	RealType tolerance = static_cast<RealType>(std::pow(10., -(decdigits-2))); // 1e-6 (-2 so as %)
	tolerance *= 2; // Allow some bit(s) small margin (2 means + or - 1 bit) of numerical error.
	// Typically 2e-13% = 2e-15 as fraction for double.

	// Sources of spot test values:

  // Many be some combinations for which the result is 'exact',
  // or at least is good to 40 decimal digits.
	// 40 decimal digits includes 128-bit significand User Defined Floating-Point types,
	
	// Best source of accurate values is:
	// Mathworld online calculator (40 decimal digits precision, suitable for up to 128-bit significands)
	// http://functions.wolfram.com/webMathematica/FunctionEvaluation.jsp?name=GammaRegularized
	// GammaRegularized is same as gamma incomplete, gamma or gamma_Q(a, x) or Q(a, z).

  // http://documents.wolfram.com/calculationcenter/v2/Functions/ListsMatrices/Statistics/PoissonDistribution.html

  // MathCAD defines ppois(k, lambda== mean) as k integer, k >=0.
  // ppois(0, 5) =  6.73794699908547e-3
  // ppois(1, 5) = 0.040427681994513;
  // ppois(10, 10) = 5.830397501929850E-001
  // ppois(10, 1) = 9.999999899522340E-001
  // ppois(5,5) = 0.615960654833065

  // qpois returns inverse Poission distribution, that is the smallest (floor) k so that ppois(k, lambda) >= p
  // p is real number, real mean lambda > 0
  // k is approximately the integer for which probability(X <= k) = p
  // when random variable X has the Poisson distribution with parameters lambda.
  // Uses discrete bisection.
  // qpois(6.73794699908547e-3, 5) = 1
  // qpois(0.040427681994513, 5) = 

  // Test Poisson with spot values from MathCAD 'known good'.

  using boost::math::poisson_distribution;
  using  ::boost::math::poisson;
  using  ::boost::math::cdf;
  using  ::boost::math::pdf;

   // Check that bad arguments throw.
   BOOST_CHECK_THROW(
   cdf(poisson_distribution<RealType>(static_cast<RealType>(0)), // mean zero is bad.
      static_cast<RealType>(0)),  // even for a good k.
      std::domain_error); // Expected error to be thrown.

    BOOST_CHECK_THROW(
   cdf(poisson_distribution<RealType>(static_cast<RealType>(-1)), // mean negative is bad.
      static_cast<RealType>(0)),
      std::domain_error);

   BOOST_CHECK_THROW(
   cdf(poisson_distribution<RealType>(static_cast<RealType>(1)), // mean unit OK,
      static_cast<RealType>(-1)),  // but negative events is bad.
      std::domain_error);

  BOOST_CHECK_THROW(
     cdf(poisson_distribution<RealType>(static_cast<RealType>(0)), // mean zero is bad.
      static_cast<RealType>(99999)),  // for any k events. 
      std::domain_error);
  
  BOOST_CHECK_THROW(
     cdf(poisson_distribution<RealType>(static_cast<RealType>(0)), // mean zero is bad.
      static_cast<RealType>(99999)),  // for any k events. 
      std::domain_error);

  BOOST_CHECK_THROW(
     quantile(poisson_distribution<RealType>(static_cast<RealType>(0)), // mean zero.
      static_cast<RealType>(0.5)),  // probability OK. 
      std::domain_error);

  BOOST_CHECK_THROW(
     quantile(poisson_distribution<RealType>(static_cast<RealType>(-1)), 
      static_cast<RealType>(-1)),  // bad probability. 
      std::domain_error);

  BOOST_CHECK_THROW(
     quantile(poisson_distribution<RealType>(static_cast<RealType>(1)), 
      static_cast<RealType>(-1)),  // bad probability. 
      std::domain_error);

  // Check some test values.

  BOOST_CHECK_CLOSE( // mode
     mode(poisson_distribution<RealType>(static_cast<RealType>(4))), // mode = mean = 4.
      static_cast<RealType>(4), // mode.
			tolerance);

  //BOOST_CHECK_CLOSE( // mode
  //   median(poisson_distribution<RealType>(static_cast<RealType>(4))), // mode = mean = 4.
  //    static_cast<RealType>(4), // mode.
		//	tolerance);
  poisson_distribution<RealType> dist4(static_cast<RealType>(40));

  BOOST_CHECK_CLOSE( // median
     median(dist4), // mode = mean = 4. median = 40.328333333333333 
      quantile(dist4, static_cast<RealType>(0.5)), // 39.332839138842637
			tolerance);

  // PDF
  BOOST_CHECK_CLOSE(
     pdf(poisson_distribution<RealType>(static_cast<RealType>(4)), // mean 4.
      static_cast<RealType>(0)),   
      static_cast<RealType>(1.831563888873410E-002), // probability.
			tolerance);

  BOOST_CHECK_CLOSE(
     pdf(poisson_distribution<RealType>(static_cast<RealType>(4)), // mean 4.
      static_cast<RealType>(2)),   
      static_cast<RealType>(1.465251111098740E-001), // probability.
			tolerance);

  BOOST_CHECK_CLOSE(
     pdf(poisson_distribution<RealType>(static_cast<RealType>(20)), // mean big.
      static_cast<RealType>(1)),   //  k small
      static_cast<RealType>(4.122307244877130E-008), // probability.
			tolerance);

  BOOST_CHECK_CLOSE(
     pdf(poisson_distribution<RealType>(static_cast<RealType>(4)), // mean 4.
      static_cast<RealType>(20)),   //  K>> mean 
      static_cast<RealType>(8.277463646553730E-009), // probability.
			tolerance);
  
  // CDF
  BOOST_CHECK_CLOSE(
     cdf(poisson_distribution<RealType>(static_cast<RealType>(1)), // mean unity.
      static_cast<RealType>(0)),  // zero k events. 
      static_cast<RealType>(3.678794411714420E-1), // probability.
			tolerance);

  BOOST_CHECK_CLOSE(
     cdf(poisson_distribution<RealType>(static_cast<RealType>(1)), // mean unity.
      static_cast<RealType>(1)),  // one k event. 
      static_cast<RealType>(7.357588823428830E-1), // probability.
			tolerance);

  BOOST_CHECK_CLOSE(
     cdf(poisson_distribution<RealType>(static_cast<RealType>(1)), // mean unity.
      static_cast<RealType>(2)),  // two k events. 
      static_cast<RealType>(9.196986029286060E-1), // probability.
			tolerance);

  BOOST_CHECK_CLOSE(
     cdf(poisson_distribution<RealType>(static_cast<RealType>(1)), // mean unity.
      static_cast<RealType>(10)),  // two k events. 
      static_cast<RealType>(9.999999899522340E-1), // probability.
			tolerance);

  BOOST_CHECK_CLOSE(
     cdf(poisson_distribution<RealType>(static_cast<RealType>(1)), // mean unity.
      static_cast<RealType>(15)),  // two k events. 
      static_cast<RealType>(9.999999999999810E-1), // probability.
			tolerance);

  BOOST_CHECK_CLOSE(
     cdf(poisson_distribution<RealType>(static_cast<RealType>(1)), // mean unity.
      static_cast<RealType>(16)),  // two k events. 
      static_cast<RealType>(9.999999999999990E-1), // probability.
			tolerance);

  BOOST_CHECK_CLOSE(
     cdf(poisson_distribution<RealType>(static_cast<RealType>(1)), // mean unity.
      static_cast<RealType>(17)),  // two k events. 
      static_cast<RealType>(1.), // probability unity for double.
			tolerance);

  BOOST_CHECK_CLOSE(
     cdf(poisson_distribution<RealType>(static_cast<RealType>(1)), // mean unity.
      static_cast<RealType>(33)),  // k events at limit for float unchecked_factorial table. 
      static_cast<RealType>(1.), // probability.
			tolerance);

  BOOST_CHECK_CLOSE(
     cdf(poisson_distribution<RealType>(static_cast<RealType>(100)), // mean 100.
      static_cast<RealType>(33)),  // k events at limit for float unchecked_factorial table. 
      static_cast<RealType>(6.328271240363390E-15), // probability is tiny.
			tolerance * static_cast<RealType>(2e11)); // 6.3495253382825722e-015 MathCAD
      // Note that there two tiny probability are much more different.

   BOOST_CHECK_CLOSE(
     cdf(poisson_distribution<RealType>(static_cast<RealType>(100)), // mean 100.
      static_cast<RealType>(34)),  // k events at limit for float unchecked_factorial table. 
      static_cast<RealType>(1.898481372109020E-14), // probability is tiny.
			tolerance*static_cast<RealType>(2e11)); //         1.8984813721090199e-014 MathCAD


 BOOST_CHECK_CLOSE(
     cdf(poisson_distribution<RealType>(static_cast<RealType>(33)), // mean = k
      static_cast<RealType>(33)),  // k events above limit for float unchecked_factorial table. 
      static_cast<RealType>(5.461191812386560E-1), // probability.
			tolerance);

 BOOST_CHECK_CLOSE(
     cdf(poisson_distribution<RealType>(static_cast<RealType>(33)), // mean = k-1
      static_cast<RealType>(34)),  // k events above limit for float unchecked_factorial table. 
      static_cast<RealType>(6.133535681502950E-1), // probability.
			tolerance);

 BOOST_CHECK_CLOSE(
     cdf(poisson_distribution<RealType>(static_cast<RealType>(1)), // mean unity.
      static_cast<RealType>(34)),  // k events above limit for float unchecked_factorial table. 
      static_cast<RealType>(1.), // probability.
			tolerance);

  BOOST_CHECK_CLOSE(
     cdf(poisson_distribution<RealType>(static_cast<RealType>(5.)), // mean
      static_cast<RealType>(5)),  // k events. 
      static_cast<RealType>(0.615960654833065), // probability.
			tolerance);
  BOOST_CHECK_CLOSE(
     cdf(poisson_distribution<RealType>(static_cast<RealType>(5.)), // mean
      static_cast<RealType>(1)),  // k events. 
      static_cast<RealType>(0.040427681994512805), // probability.
			tolerance);

  BOOST_CHECK_CLOSE(
     cdf(poisson_distribution<RealType>(static_cast<RealType>(5.)), // mean
      static_cast<RealType>(0)),  // k events (uses special case formula, not gamma). 
      static_cast<RealType>(0.006737946999085467), // probability.
			tolerance);

  BOOST_CHECK_CLOSE(
     cdf(poisson_distribution<RealType>(static_cast<RealType>(1.)), // mean
      static_cast<RealType>(0)),  // k events (uses special case formula, not gamma). 
      static_cast<RealType>(0.36787944117144233), // probability.
			tolerance);

  BOOST_CHECK_CLOSE(
     cdf(poisson_distribution<RealType>(static_cast<RealType>(10.)), // mean
      static_cast<RealType>(10)),  // k events. 
      static_cast<RealType>(0.5830397501929856), // probability.
			tolerance);

  BOOST_CHECK_CLOSE(
     cdf(poisson_distribution<RealType>(static_cast<RealType>(4.)), // mean
      static_cast<RealType>(5)),  // k events. 
      static_cast<RealType>(0.785130387030406), // probability.
			tolerance);

  // complement CDF
  BOOST_CHECK_CLOSE( // Complement CDF
     cdf(complement(poisson_distribution<RealType>(static_cast<RealType>(4.)), // mean
      static_cast<RealType>(5))),  // k events. 
      static_cast<RealType>(1 - 0.785130387030406), // probability.
			tolerance);

  BOOST_CHECK_CLOSE( // Complement CDF
     cdf(complement(poisson_distribution<RealType>(static_cast<RealType>(4.)), // mean
      static_cast<RealType>(0))),  // Zero k events (uses special case formula, not gamma).
      static_cast<RealType>(0.98168436111126578), // probability.
			tolerance);
  BOOST_CHECK_CLOSE( // Complement CDF
     cdf(complement(poisson_distribution<RealType>(static_cast<RealType>(1.)), // mean
      static_cast<RealType>(0))),  // Zero k events (uses special case formula, not gamma).
      static_cast<RealType>(0.63212055882855767), // probability.
			tolerance);

  // Example where k is bigger than max_factorial (>34 for float)
  // (therefore using log gamma so perhaps less accurate).
  BOOST_CHECK_CLOSE(
     cdf(poisson_distribution<RealType>(static_cast<RealType>(40.)), // mean
      static_cast<RealType>(40)),  // k events. 
      static_cast<RealType>(0.5419181783625430), // probability.
			tolerance);

   // Quantile & complement.
  BOOST_CHECK_CLOSE(
    boost::math::quantile(
         poisson_distribution<RealType>(5),  // mean.
         static_cast<RealType>(0.615960654833065)),  //  probability.
         static_cast<RealType>(5.), // Expect k = 5
         tolerance/5); // 

  // EQUAL is too optimistic - fails [5.0000000000000124 != 5]
  // BOOST_CHECK_EQUAL(boost::math::quantile( // 
  //       poisson_distribution<RealType>(5.),  // mean.
  //       static_cast<RealType>(0.615960654833065)),  //  probability.
  //       static_cast<RealType>(5.)); // Expect k = 5 events.
 
  BOOST_CHECK_CLOSE(boost::math::quantile(
         poisson_distribution<RealType>(4),  // mean.
         static_cast<RealType>(0.785130387030406)),  //  probability.
         static_cast<RealType>(5.), // Expect k = 5 events.
         tolerance/5); 

  // Check on quantile of other examples of inverse of cdf.
  BOOST_CHECK_CLOSE( 
     cdf(poisson_distribution<RealType>(static_cast<RealType>(10.)), // mean
      static_cast<RealType>(10)),  // k events. 
      static_cast<RealType>(0.5830397501929856), // probability.
			tolerance);

  BOOST_CHECK_CLOSE(boost::math::quantile( // inverse of cdf above.
         poisson_distribution<RealType>(10.),  // mean.
         static_cast<RealType>(0.5830397501929856)),  //  probability.
         static_cast<RealType>(10.), // Expect k = 10 events.
         tolerance/5); 


  BOOST_CHECK_CLOSE(
     cdf(poisson_distribution<RealType>(static_cast<RealType>(4.)), // mean
      static_cast<RealType>(5)),  // k events. 
      static_cast<RealType>(0.785130387030406), // probability.
			tolerance);

  BOOST_CHECK_CLOSE(boost::math::quantile( // inverse of cdf above.
         poisson_distribution<RealType>(4.),  // mean.
         static_cast<RealType>(0.785130387030406)),  //  probability.
         static_cast<RealType>(5.), // Expect k = 10 events.
         tolerance/5); 



  //BOOST_CHECK_CLOSE(boost::math::quantile(
  //       poisson_distribution<RealType>(5),  // mean.
  //       static_cast<RealType>(0.785130387030406)),  //  probability.
  //        // 6.1882832344329559 result but MathCAD givest smallest integer ppois(k, mean) >= prob
  //       static_cast<RealType>(6.), // Expect k = 6 events. 
  //       tolerance/5); 

  //BOOST_CHECK_CLOSE(boost::math::quantile(
  //       poisson_distribution<RealType>(5),  // mean.
  //       static_cast<RealType>(0.77)),  //  probability.
  //        // 6.1882832344329559 result but MathCAD givest smallest integer ppois(k, mean) >= prob
  //       static_cast<RealType>(7.), // Expect k = 6 events. 
  //       tolerance/5); 

  //BOOST_CHECK_CLOSE(boost::math::quantile(
  //       poisson_distribution<RealType>(5),  // mean.
  //       static_cast<RealType>(0.75)),  //  probability.
  //        // 6.1882832344329559 result but MathCAD givest smallest integer ppois(k, mean) >= prob
  //       static_cast<RealType>(6.), // Expect k = 6 events. 
  //       tolerance/5); 

  BOOST_CHECK_CLOSE(
    boost::math::quantile(
         complement(
           poisson_distribution<RealType>(4),
           static_cast<RealType>(1 - 0.785130387030406))),  // complement.
           static_cast<RealType>(5), // Expect k = 5 events.
         tolerance/5);

  BOOST_CHECK_EQUAL(boost::math::quantile( // Check case when probability < cdf(0) (== pdf(0))
         poisson_distribution<RealType>(1),  // mean is small, so cdf and pdf(0) are about 0.35.
         static_cast<RealType>(0.0001)),  //  probability < cdf(0).
         static_cast<RealType>(0)); // Expect k = 0 events exactly.
          
  BOOST_CHECK_EQUAL(
    boost::math::quantile(
         complement(
           poisson_distribution<RealType>(1),
           static_cast<RealType>(0.9999))),  // complement, so 1-probability < cdf(0)
           static_cast<RealType>(0)); // Expect k = 0 events exactly.


} // template <class RealType>void test_spots(RealType)

//

int test_main(int, char* [])
{
  // Check that can construct normal distribution using the two convenience methods:
  using namespace boost::math;
  poisson myp1(2); // Using typedef
	poisson_distribution<> myp2(2); // Using default RealType double.

	// Basic sanity-check spot values.
#ifdef BOOST_MATH_THROW_ON_DOMAIN_ERROR
  cout << "BOOST_MATH_THROW_ON_DOMAIN_ERROR" << " is defined to throw on domain error." << endl;
#else
  cout << "BOOST_MATH_THROW_ON_DOMAIN_ERROR" << " is NOT defined, so NO throw on domain error." << endl;
#endif

// Some plain double examples & tests:
  cout.precision(17); // double max_digits10
  cout.setf(ios::showpoint);
  
  poisson mypoisson(4.); // // mean = 4, default FP type is double.
  cout << "mean(mypoisson, 4.) == " << mean(mypoisson) << endl;
  cout << "mean(mypoisson, 0.) == " << mean(mypoisson) << endl;
  cout << "cdf(mypoisson, 2.) == " << cdf(mypoisson, 2.) << endl;
  cout << "pdf(mypoisson, 2.) == " << pdf(mypoisson, 2.) << endl;
  
  //poisson mydudpoisson(0.);
  // throws (if BOOST_MATH_THROW_ON_DOMAIN_ERROR defined to enable).

#define BOOST_MATH_THROW_ON_LOGIC_ERROR
  
  BOOST_CHECK_THROW(poisson mydudpoisson(-1), std::domain_error);// Mean must be > 0.
  BOOST_CHECK_THROW(poisson mydudpoisson(-1), std::logic_error);// Mean must be > 0.
  // Passes the check because logic_error is a parent????
  // BOOST_CHECK_THROW(poisson mydudpoisson(-1), std::overflow_error); // fails the check
  // because overflow_error is unrelated - except from std::exception
  BOOST_CHECK_THROW(cdf(mypoisson, -1), std::domain_error); // k must be >= 0

  BOOST_CHECK_EQUAL(mean(mypoisson), 4.);
  BOOST_CHECK_CLOSE(
  pdf(mypoisson, 2.),  // k events = 2. 
    1.465251111098740E-001, // probability.
		5e-13);

  BOOST_CHECK_CLOSE(
  cdf(mypoisson, 2.),  // k events = 2. 
    0.238103305553545, // probability.
		5e-13);


#if 0
  // Compare cdf from finite sum of pdf and gamma_Q.
  using boost::math::cdf;
  using boost::math::pdf;

  double mean = 4.;
  cout.precision(17); // double max_digits10
  cout.setf(ios::showpoint);
  cout << showpoint << endl;  // Ensure trailing zeros are shown.
  // This also helps show the expected precision max_digits10
  //cout.unsetf(ios::showpoint); // No trailing zeros are shown.

  cout << "k          pdf                     sum                  cdf                   diff" << endl;
  double sum = 0.;
  for (int i = 0; i <= 50; i++)
  {
   cout << i << ' ' ;
   double p =  pdf(poisson_distribution<double>(mean), static_cast<double>(i));
   sum += p;

   cout << p << ' ' << sum << ' ' 
   << cdf(poisson_distribution<double>(mean), static_cast<double>(i)) << ' ';
     {
       cout << boost::math::gamma_Q<double>(i+1, mean); // cdf
       double diff = boost::math::gamma_Q<double>(i+1, mean) - sum; // cdf -sum
       cout << setprecision (2) << ' ' << diff; // 0 0 to 4, 1 eps 5 to 9, 10 to 20 2 eps, 21 upwards 3 eps
      
     }
    BOOST_CHECK_CLOSE(
    cdf(mypoisson, static_cast<double>(i)),
      sum, // of pdfs.
	   4e-14); // Fails at 2e-14
   // This call puts the precision etc back to default 6 !!!
   cout << setprecision(17) << showpoint;


     cout << endl;
  }

   cout << cdf(poisson_distribution<double>(5), static_cast<double>(0)) << ' ' << endl; // 0.006737946999085467
   cout << cdf(poisson_distribution<double>(5), static_cast<double>(1)) << ' ' << endl; // 0.040427681994512805
   cout << cdf(poisson_distribution<double>(2), static_cast<double>(3)) << ' ' << endl; // 0.85712346049854715 
#endif

   {
     for (int i = 1; i < 100; i++)
     {
       poisson_distribution<double> distn(static_cast<double>(i));
       cout << i << ' ' << median(distn) << ' ' << quantile(distn, 0.5) << ' ' 
         << median(distn) - quantile(distn, 0.5) << endl;
     }
   }

	// (Parameter value, arbitrarily zero, only communicates the floating-point type).
  test_spots(0.0F); // Test float.
	test_spots(0.0); // Test double.
  if (numeric_limits<long double>::digits10 > numeric_limits<double>::digits10)
  { // long double is better than double (so not MSVC where they are same).
	  test_spots(0.0L); // Test long double.
  }

  #if !BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x582))
  test_spots(boost::math::concepts::real_concept(0.)); // Test real concept.
  #endif

	return 0;
} // int test_main(int, char* [])

/*

Output:

------ Build started: Project: test_poisson, Configuration: Debug Win32 ------
Compiling...
test_poisson.cpp
Linking...
Autorun "i:\boost-06-05-03-1300\libs\math\test\Math_test\debug\test_poisson.exe"
Running 1 test case...
BOOST_MATH_THROW_ON_DOMAIN_ERROR is defined to throw on domain error.
mean(mypoisson, 4.) == 4.0000000000000000
mean(mypoisson, 0.) == 4.0000000000000000
cdf(mypoisson, 2.) == 0.23810330555354431
pdf(mypoisson, 2.) == 0.14652511110987343
*** No errors detected
Build Time 0:06
Build log was saved at "file://i:\boost-06-05-03-1300\libs\math\test\Math_test\test_poisson\Debug\BuildLog.htm"
test_poisson - 0 error(s), 0 warning(s)
========== Build: 1 succeeded, 0 failed, 0 up-to-date, 0 skipped ==========


------ Build started: Project: test_poisson, Configuration: Debug Win32 ------
Compiling...
test_poisson.cpp
Linking...
Autorun "i:\boost-06-05-03-1300\libs\math\test\Math_test\debug\test_poisson.exe"
Running 1 test case...
BOOST_MATH_THROW_ON_DOMAIN_ERROR is defined to throw on domain error.
Mean = 4
k          pdf                     sum                  cdf                   diff
0 0.018315638888734179 0.018315638888734179 0.018315638888734179 0.018315638888734179
1 0.073262555554936715 0.091578194443670893 0.091578194443670893 0.091578194443670893 0
2 0.14652511110987343 0.23810330555354431 0.23810330555354431 0.23810330555354431 0
3 0.19536681481316456 0.43347012036670884 0.43347012036670884 0.43347012036670884 0
4 0.19536681481316456 0.62883693517987338 0.62883693517987338 0.62883693517987338 0
5 0.15629345185053165 0.78513038703040505 0.78513038703040505 0.78513038703040516 1.1102230246251565e-016
6 0.10419563456702111 0.88932602159742613 0.88932602159742613 0.88932602159742624 1.1102230246251565e-016
7 0.059540362609726345 0.94886638420715252 0.94886638420715252 0.94886638420715264 1.1102230246251565e-016
8 0.029770181304863173 0.97863656551201572 0.97863656551201572 0.97863656551201583 1.1102230246251565e-016
9 0.013231191691050298 0.99186775720306597 0.99186775720306597 0.99186775720306608 1.1102230246251565e-016
10 0.0052924766764201195 0.99716023387948605 0.99716023387948605 0.99716023387948627 2.2204460492503131e-016
11 0.0019245369732436798 0.99908477085272973 0.99908477085272973 0.99908477085272995 2.2204460492503131e-016
12 0.00064151232441456 0.99972628317714429 0.99972628317714429 0.99972628317714451 2.2204460492503131e-016
13 0.00019738840751217228 0.99992367158465645 0.99992367158465645 0.99992367158465667 2.2204460492503131e-016
14 5.6396687860620656e-005 0.99998006827251706 0.99998006827251706 0.99998006827251729 2.2204460492503131e-016
15 1.5039116762832175e-005 0.99999510738927988 0.99999510738927988 0.9999951073892801 2.2204460492503131e-016
16 3.7597791907080438e-006 0.99999886716847064 0.99999886716847064 0.99999886716847086 2.2204460492503131e-016
17 8.8465392722542207e-007 0.99999975182239786 0.99999975182239786 0.99999975182239809 2.2204460492503131e-016
18 1.9658976160564933e-007 0.99999994841215945 0.99999994841215945 0.99999994841215967 2.2204460492503131e-016
19 4.1387318232768281e-008 0.99999998979947768 0.99999998979947768 0.99999998979947791 2.2204460492503131e-016
20 8.2774636465536562e-009 0.99999999807694129 0.99999999807694151 0.99999999807694151 2.2204460492503131e-016
21 1.5766597422006965e-009 0.99999999965360098 0.99999999965360131 0.99999999965360131 3.3306690738754696e-016
22 2.8666540767285388e-010 0.99999999994026634 0.99999999994026667 0.99999999994026667 3.3306690738754696e-016
23 4.9854853508322414e-011 0.99999999999012124 0.99999999999012157 0.99999999999012157 3.3306690738754696e-016
24 8.3091422513870696e-012 0.99999999999843037 0.9999999999984307 0.9999999999984307 3.3306690738754696e-016
25 1.3294627602219311e-012 0.99999999999975986 0.99999999999976019 0.99999999999976019 3.3306690738754696e-016
26 2.0453273234183554e-013 0.99999999999996436 0.99999999999996469 0.99999999999996469 3.3306690738754696e-016
27 3.0301145532123785e-014 0.99999999999999467 0.999999999999995 0.999999999999995 3.3306690738754696e-016
28 4.3287350760176837e-015 0.999999999999999 0.99999999999999933 0.99999999999999933 3.3306690738754696e-016
29 5.9706690703692191e-016 0.99999999999999956 0.99999999999999989 0.99999999999999989 3.3306690738754696e-016
30 7.9608920938256257e-017 0.99999999999999967 1 1 3.3306690738754696e-016
31 1.0272118830742743e-017 0.99999999999999967 1 1 3.3306690738754696e-016
32 1.2840148538428429e-018 0.99999999999999967 1 1 3.3306690738754696e-016
33 1.5563816410216277e-019 0.99999999999999967 1 1 3.3306690738754696e-016
34 1.8310372247313267e-020 0.99999999999999967 1 1 3.3306690738754696e-016
35 2.092613971121516e-021 0.99999999999999967 1 1 3.3306690738754696e-016
36 2.3251266345794622e-022 0.99999999999999967 1 1 3.3306690738754696e-016
37 2.5136504157615809e-023 0.99999999999999967 1 1 3.3306690738754696e-016
38 2.645947806064822e-024 0.99999999999999967 1 1 3.3306690738754696e-016
39 2.7137926216049458e-025 0.99999999999999967 1 1 3.3306690738754696e-016
40 2.713792621604946e-026 0.99999999999999967 1 1 3.3306690738754696e-016
41 2.6476025576633616e-027 0.99999999999999967 1 1 3.3306690738754696e-016
42 2.5215262453936778e-028 0.99999999999999967 1 1 3.3306690738754696e-016
43 2.3456058096685376e-029 0.99999999999999967 1 1 3.3306690738754696e-016
44 2.1323689178804887e-030 0.99999999999999967 1 1 3.3306690738754696e-016
45 1.89543903811599e-031 0.99999999999999967 1 1 3.3306690738754696e-016
46 1.6482078592312956e-032 0.99999999999999967 1 1 3.3306690738754696e-016
47 1.4027300929628047e-033 0.99999999999999967 1 1 3.3306690738754696e-016
48 1.1689417441356707e-034 0.99999999999999967 1 1 3.3306690738754696e-016
49 9.5423815847809852e-036 0.99999999999999967 1 1 3.3306690738754696e-016
50 7.633905267824789e-037 0.99999999999999967 1 1 3.3306690738754696e-016
0.006737946999085467 
0.040427681994512805 
0.85712346049854715 
*** No errors detected
Build Time 0:06
Build log was saved at "file://i:\boost-06-05-03-1300\libs\math\test\Math_test\test_poisson\Debug\BuildLog.htm"
test_poisson - 0 error(s), 0 warning(s)
========== Build: 1 succeeded, 0 failed, 0 up-to-date, 0 skipped ==========


*/
