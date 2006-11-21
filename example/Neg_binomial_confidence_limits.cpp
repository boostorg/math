// Copyright John Maddock 2006
// Copyright Paul A. Bristow 2006
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifdef _MSC_VER
#  pragma warning(disable: 4512) // assignment operator could not be generated.
#  pragma warning(disable: 4510) // default constructor could not be generated.
#  pragma warning(disable: 4610) // can never be instantiated - user defined constructor required.
#endif

#include <iostream>
using std::cout;
using std::endl;
#include <iomanip>
#include <boost/math/distributions/negative_binomial.hpp>

void confidence_limits_on_frequency(unsigned trials, unsigned successes)
{
   // trials = Total number of trials.
   // successes = Total number of observed successes.
   // failures = trials - successes.
   //
   // Calculate confidence limits for an observed
   // frequency of occurance that follows a negative binomial distribution.

   using namespace std;
   using namespace boost::math;

   // Print out general info:
   cout <<
      "___________________________________________\n"
      "2-Sided Confidence Limits For Success Ratio\n"
      "___________________________________________\n\n";
   cout << setprecision(7);
   cout << setw(40) << left << "Number of trials" << " =  " << trials << "\n";
   cout << setw(40) << left << "Number of successes" << " =  " << successes << "\n";
   cout << setw(40) << left << "Number of failures" << " =  " << trials - successes << "\n";
   cout << setw(40) << left << "Observed probability of occurrence" << " =  " << double(successes) / trials << "\n";

   // Define a table of significance levels:
   double alpha[] = { 0.5, 0.25, 0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001 };

   // Print table header:

   cout << "\n\n"
           "___________________________________________\n"
           "Confidence        Lower          Upper\n"
           " Value (%)        Limit          Limit\n"
           "___________________________________________\n";

   // Now print out the data for the alpha table rows.
   for(unsigned i = 0; i < sizeof(alpha)/sizeof(alpha[0]); ++i)
   {
      // Confidence value:
      cout << fixed << setprecision(3) << setw(10) << right << 100 * (1-alpha[i]);
      // calculate bounds:
      double l = negative_binomial::estimate_lower_bound_on_p(trials, successes, alpha[i]/2);
      double u = negative_binomial::estimate_upper_bound_on_p(trials, successes, alpha[i]/2);
      // Print Limits:
      cout << fixed << setprecision(5) << setw(15) << right << l;
      cout << fixed << setprecision(5) << setw(15) << right << u << endl;
   }
   cout << endl;
} // void confidence_limits_on_frequency(unsigned trials, unsigned successes)

void estimate_number_of_trials(double failures, double p)
{
   // Define a table of significance levels:
   double alpha[] = { 0.5, 0.25, 0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001 };
   using namespace boost::math;
    // static RealType estimate_number_of_trials(
  //RealType k,     // number of failures, k >= 0.
  //RealType p,     // success fraction 0 <= p <= 1.
  //RealType probability) // probability threshold 0 <= p <= 0.
  using namespace std;
  cout << "\n""Number of failures = " << failures;
  cout << ",   Success ratio = " << 100 * p << "%" << endl;
  cout << "Confidence %     trials   complement" << endl;
  for(unsigned i = 0; i < sizeof(alpha)/sizeof(alpha[0]); ++i)
   { // Confidence values %:
      cout << fixed << setprecision(3) << setw(10) << right << 100 * (1-alpha[i]) << "      "
      // estimate_number_of_trials
      << setw(6) << right<< int(negative_binomial::estimate_number_of_trials(failures, p, alpha[i]/2)) << "  " 
      << setw(6) << right<< int(negative_binomial::estimate_number_of_trials(boost::math::complement(failures, p, alpha[i]/2)))
      << endl;
   }
   cout << endl;
} // void estimate_number_of_trials(double fails, double p)



int main()
{
  confidence_limits_on_frequency(20, 2);
  confidence_limits_on_frequency(200, 20);
  confidence_limits_on_frequency(2000, 200);

  estimate_number_of_trials(5, 0.5);
  estimate_number_of_trials(50, 0.5);
  estimate_number_of_trials(500, 0.5);
  estimate_number_of_trials(50, 0.1);
  estimate_number_of_trials(500, 0.1);
  estimate_number_of_trials(5, 0.9);

  return 0;
} // int main()

/*

Output:

------ Build started: Project: neg_binomial_confidence_levels, Configuration: Debug Win32 ------
Compiling...
Neg_binomial_confidence_limits.cpp
Linking...
Autorun "i:\boost-06-05-03-1300\libs\math\test\Math_test\debug\neg_binomial_confidence_levels.exe"
___________________________________________
2-Sided Confidence Limits For Success Ratio
___________________________________________
Number of trials                         =  20
Number of successes                      =  2
Number of failures                       =  18
Observed probability of occurrence       =  0.1
___________________________________________
Confidence        Lower          Upper
 Value (%)        Limit          Limit
___________________________________________
    50.000        0.08701        0.18675
    75.000        0.06229        0.23163
    90.000        0.04217        0.28262
    95.000        0.03207        0.31698
    99.000        0.01764        0.38713
    99.900        0.00786        0.47093
    99.990        0.00358        0.54084
    99.999        0.00165        0.60020
___________________________________________
2-Sided Confidence Limits For Success Ratio
___________________________________________
Number of trials                         =  200
Number of successes                      =  20
Number of failures                       =  180
Observed probability of occurrence       =  0.1000000
___________________________________________
Confidence        Lower          Upper
 Value (%)        Limit          Limit
___________________________________________
    50.000        0.08929        0.11824
    75.000        0.08023        0.12959
    90.000        0.07144        0.14199
    95.000        0.06618        0.15021
    99.000        0.05664        0.16698
    99.900        0.04676        0.18756
    99.990        0.03944        0.20571
    99.999        0.03371        0.22226
___________________________________________
2-Sided Confidence Limits For Success Ratio
___________________________________________
Number of trials                         =  2000
Number of successes                      =  200
Number of failures                       =  1800
Observed probability of occurrence       =  0.1000000
___________________________________________
Confidence        Lower          Upper
 Value (%)        Limit          Limit
___________________________________________
    50.000        0.09585        0.10491
    75.000        0.09277        0.10822
    90.000        0.08963        0.11172
    95.000        0.08767        0.11399
    99.000        0.08390        0.11850
    99.900        0.07966        0.12385
    99.990        0.07621        0.12845
    99.999        0.07325        0.13256
Number of failures = 5.00000,   Success ratio = 50.00000%
Confidence %     trials   complement
    50.000          13       8
    75.000          15       7
    90.000          17       6
    95.000          19       6
    99.000          23       5
    99.900          27       5
    99.990          32       5
    99.999          36       5
Number of failures = 50.000,   Success ratio = 50.000%
Confidence %     trials   complement
    50.000         108      94
    75.000         113      90
    90.000         118      85
    95.000         122      83
    99.000         130      78
    99.900         139      73
    99.990         147      69
    99.999         155      65
Number of failures = 500.000,   Success ratio = 50.000%
Confidence %     trials   complement
    50.000        1022     979
    75.000        1038     965
    90.000        1054     950
    95.000        1064     940
    99.000        1085     922
    99.900        1110     902
    99.990        1131     885
    99.999        1150     870
Number of failures = 50.000,   Success ratio = 10.000%
Confidence %     trials   complement
    50.000         553     462
    75.000         588     432
    90.000         626     403
    95.000         651     385
    99.000         701     352
    99.900         763     317
    99.990         818     289
    99.999         869     266
Number of failures = 500.000,   Success ratio = 10.000%
Confidence %     trials   complement
    50.000        5150    4864
    75.000        5254    4766
    90.000        5364    4665
    95.000        5434    4602
    99.000        5574    4480
    99.900        5739    4341
    99.990        5880    4227
    99.999        6006    4129
Number of failures = 5.000,   Success ratio = 90.000%
Confidence %     trials   complement
    50.000           6       5
    75.000           7       5
    90.000           7       5
    95.000           8       5
    99.000           9       5
    99.900          10       5
    99.990          12       5
    99.999          13       5
Build Time 0:03
Build log was saved at "file://i:\boost-06-05-03-1300\libs\math\test\Math_test\neg_binomial_confidence_levels\Debug\BuildLog.htm"
neg_binomial_confidence_levels - 0 error(s), 0 warning(s)
========== Build: 1 succeeded, 0 failed, 0 up-to-date, 0 skipped ==========

*/

