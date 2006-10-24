// negative_binomial_example2.cpp

// Copyright Paul A. Bristow 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// Simple examples for Negative Binomial Distribution.

#define BOOST_MATH_THROW_ON_DOMAIN_ERROR

#ifdef _MSC_VER
#  pragma warning(disable: 4127) // conditional expression is constant.
#  pragma warning(disable: 4100) // unreferenced formal parameter.
#  pragma warning(disable: 4512) // assignment operator could not be generated.
#endif

#include <boost/math/distributions/negative_binomial.hpp> // for negative_binomial_distribution
using boost::math::negative_binomial_distribution;

// In a sequence of trials or events
// (Bernoulli, independent, yes or no, succeed or fail)
// with success_fraction probability p,
// negative_binomial is the probability that k or fewer failures
// preceed the r th trial's success.

#include <boost/math/special_functions/gamma.hpp>
  using boost::math::lgamma;  // log gamma

//#include <boost/math/concepts/real_concept.hpp> // for real_concept
//using ::boost::math::concepts::real_concept;

//#include <boost/test/included/test_exec_monitor.hpp> // for test_main
//#include <boost/test/floating_point_comparison.hpp> // for BOOST_CHECK_CLOSE

#include <iostream>
using std::cout;
using std::endl;
using std::setprecision;
using std::showpoint;
using std::setw;
using std::left;
#include <limits>
using std::numeric_limits;

int main(int, char* [])
{
#ifdef BOOST_MATH_THROW_ON_DOMAIN_ERROR
  cout << "BOOST_MATH_THROW_ON_DOMAIN_ERROR" << " is defined to throw on domain error." << endl;
#else
  cout << "BOOST_MATH_THROW_ON_DOMAIN_ERROR" << " is NOT defined, so NO throw on domain error." << endl;
#endif
  cout << "negative_binomial distribution - simples examples" << endl;
  using namespace boost::math;
  cout << setprecision(17) << showpoint << endl; // max_digits10 precision including trailing zeros.

  negative_binomial_distribution<double> my8dist(8., 0.25);
  // 8 successes (r), 0.25 success fraction = 25% or 1 in 4 successes.
  // Note: double values (matching the distribution definition) avoid the need for any casting.

  cout << "mean(my8dist) = " << mean(my8dist) << endl; // 
  cout << "my8dist.successes() = " << my8dist.successes()  << endl;
  // r th trial is successful, after r-1 = k failures.
  cout << "my8dist.success_fraction() = " << my8dist.success_fraction()  << endl;
  // failures/successes. 
  cout << "cdf(my8dist, 2.) = " << cdf(my8dist, 2.) << endl; // 4.1580200195313E-4

  cout << "cdf(my8dist, 8.) = " << cdf(my8dist, 8.) << endl;
  cout << "cdf(complement(my8dist, 8.)) = " << cdf(complement(my8dist, 8.)) << endl;
  cout << "cdf + complement = " << cdf(my8dist, 8.) + cdf(complement(my8dist, 8.))  << endl;
  // Expect unity.
  // cout << "cdf(complement(my8dist, 8.)) = " << pdf(complement(my8dist, 8.)) << endl;
  // Note: No complement for pdf! 
  double sum = 0.;
  int k = 20;

  for(signed i = 0; i <= k; ++i)
  {
    sum += pdf(my8dist, double(i));
  }
  // Compare with 
  double cdf8 = cdf(my8dist, static_cast<double>(k));
  double diff = sum - cdf8;

  cout << setprecision(17) << showpoint << "Sum pdfs = " << sum <<' '  // 0.40025683281803714 
  << ", cdf = " << cdf(my8dist, static_cast<double>(k))
  << ", difference = " 
  << diff/ std::numeric_limits<double>::epsilon()
  << " in epsilon units." << endl;
  // 0.40025683281803681
  //double tol = boost::math::tools::epsilon<double>() * 100 * 10;  // 100 for %, so 10 eps
  // Use boost::math::tools::epsilon rather than std::numeric_limits to cover
  // RealTypes that do not specialize numeric_limits.
  //BOOST_CHECK_CLOSE(sum, cdf(my8dist, static_cast<double>(20)), tol);


  // Print a list of values that can be used to plot
  // using Excel, or some superior graphical display tool.

  int maxk = static_cast<int>(2. * my8dist.successes() /  my8dist.success_fraction());
  // This shows most of the range of interest.
  for (int k = 0; k < maxk; k++)
  {
    cout << setprecision(17) << showpoint
      << left << setw(3) << k  << ", "
      << left << setw(23) << pdf(my8dist, static_cast<double>(k)) << ", "
      << left << setw(23) << cdf(my8dist, static_cast<double>(k))
      << endl;
  }


  return 0;
} // int main(int, char* [])

/*

Output is

------ Build started: Project: negative_binomial_example2, Configuration: Debug Win32 ------
Compiling...
negative_binomial_example2.cpp
Linking...
Autorun "i:\boost-06-05-03-1300\libs\math\test\Math_test\debug\negative_binomial_example2.exe"
BOOST_MATH_THROW_ON_DOMAIN_ERROR is defined to throw on domain error.
negative_binomial distribution - simples examples
mean(my8dist) = 24.000000000000000
my8dist.successes() = 8.0000000000000000
my8dist.success_fraction() = 0.25000000000000000
cdf(my8dist, 2.) = 0.00041580200195312495
cdf(my8dist, 8.) = 0.027129956288263198
cdf(complement(my8dist, 8.)) = 0.97287004371173680
cdf + complement = 1.0000000000000000
Sum pdfs = 0.40025683281803714 , cdf = 0.40025683281803681, difference = 1.5000000000000000 in epsilon units.
0  , 1.5258789062499998e-005, 1.5258789062500003e-005
1  , 9.1552734374999959e-005, 0.00010681152343750000 
2  , 0.00030899047851562522 , 0.00041580200195312495 
3  , 0.00077247619628906239 , 0.0011882781982421873  
4  , 0.0015932321548461931  , 0.0027815103530883785  
5  , 0.0028678178787231463  , 0.0056493282318115226  
6  , 0.0046602040529251142  , 0.010309532284736632   
7  , 0.0069903060793876605  , 0.017299838364124295   
8  , 0.0098301179241389071  , 0.027129956288263198   
9  , 0.013106823898851858   , 0.040236780187115066   
10 , 0.016711200471036133   , 0.056947980658151202   
11 , 0.020509200578089803   , 0.077457181236240999   
12 , 0.024354675686481628   , 0.10181185692272264    
13 , 0.028101548869017345   , 0.12991340579173991    
14 , 0.031614242477644432   , 0.16152764826938434    
15 , 0.034775666725408945   , 0.19630331499479323    
16 , 0.037492515688331500   , 0.23379583068312468    
17 , 0.039697957787645122   , 0.27349378847076972    
18 , 0.041352039362130291   , 0.31484582783289999    
19 , 0.042440250924291600   , 0.35728607875719159    
20 , 0.042970754060845266   , 0.40025683281803681    
21 , 0.042970754060845252   , 0.44322758687888208    
22 , 0.042482450037426567   , 0.48571003691630860    
23 , 0.041558918514873790   , 0.52726895543118235    
24 , 0.040260202311284035   , 0.56752915774246648    
25 , 0.038649794218832613   , 0.60617895196129901    
26 , 0.036791631035234952   , 0.64297058299653465    
27 , 0.034747651533277447   , 0.67771823452981150    
28 , 0.032575923312447616   , 0.71029415784225836    
29 , 0.030329307911589130   , 0.74062346575384863    
30 , 0.028054609818219937   , 0.76867807557206813    
31 , 0.025792141284492535   , 0.79447021685656094    
32 , 0.023575629142856481   , 0.81804584599941688    
33 , 0.021432390129869510   , 0.83947823612928674    
34 , 0.019383705779220179   , 0.85886194190850684    
35 , 0.017445335201298227   , 0.87630727710980472    
36 , 0.015628112784496315   , 0.89193538989430132    
37 , 0.013938587078064250   , 0.90587397697236538    
38 , 0.012379666154859706   , 0.91825364312722502    
39 , 0.010951243136991251   , 0.92920488626421649    
40 , 0.0096507830144735348  , 0.93885566927868991    
41 , 0.0084738582566109191  , 0.94732952753530109    
42 , 0.0074146259745345999  , 0.95474415350983566    
43 , 0.0064662435824429246  , 0.96121039709227840    
44 , 0.0056212231142827853  , 0.96683162020656122    
45 , 0.0048717266990450708  , 0.97170334690560634    
46 , 0.0042098073105878457  , 0.97591315421619418    
47 , 0.0036275999165703964  , 0.97954075413276454    
48 , 0.0031174686783026818  , 0.98265822281106729    
49 , 0.0026721160099737293  , 0.98533033882104104    
50 , 0.0022846591885275322  , 0.98761499800956853    
51 , 0.0019486798960970076  , 0.98956367790566557    
52 , 0.0016582516423517852  , 0.99122192954801736    
53 , 0.0014079495076571762  , 0.99262987905567457    
54 , 0.0011928461106539983  , 0.99382272516632852    
55 , 0.0010084971662801955  , 0.99483122233260879    
56 , 0.00085091948404891532 , 0.99568214181665760    
57 , 0.00071656377604119553 , 0.99639870559269883    
58 , 0.00060228420831048911 , 0.99700098980100937    
59 , 0.00050530624256557902 , 0.99750629604357488    
60 , 0.00042319397814867321 , 0.99792949002172360    
61 , 0.00035381791615708398 , 0.99828330793788067    
62 , 0.00029532382517950264 , 0.99857863176306016    
63 , 0.00024610318764958615 , 0.99882473495070978    
Build Time 0:03
Build log was saved at "file://i:\boost-06-05-03-1300\libs\math\test\Math_test\neative_binomial_example2\Debug\BuildLog.htm"
negative_binomial_example2 - 0 error(s), 0 warning(s)
========== Build: 1 succeeded, 0 failed, 0 up-to-date, 0 skipped ==========


*/
