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
   // frequency of occurrence that follows a negative binomial distribution.

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



int main()
{
  confidence_limits_on_frequency(20, 2);
  confidence_limits_on_frequency(200, 20);
  confidence_limits_on_frequency(2000, 200);

  return 0;
} // int main()

/*

Output:

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

*/

