// Copyright John Maddock 2006
// Copyright Paul A. Bristow 2007
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/distributions/negative_binomial.hpp>

#include <iostream>
using std::cout;
using std::endl;
#include <iomanip>

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
      double l = negative_binomial::find_lower_bound_on_p(trials, successes, alpha[i]/2);
      double u = negative_binomial::find_upper_bound_on_p(trials, successes, alpha[i]/2);
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
    50.000        0.04812        0.13554
    75.000        0.03078        0.17727
    90.000        0.01807        0.22637
    95.000        0.01235        0.26028
    99.000        0.00530        0.33111
    99.900        0.00164        0.41802
    99.990        0.00051        0.49202
    99.999        0.00016        0.55574
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
    50.000        0.08462        0.11350
    75.000        0.07580        0.12469
    90.000        0.06726        0.13695
    95.000        0.06216        0.14508
    99.000        0.05293        0.16170
    99.900        0.04343        0.18212
    99.990        0.03641        0.20017
    99.999        0.03095        0.21664
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
    50.000        0.09536        0.10445
    75.000        0.09228        0.10776
    90.000        0.08916        0.11125
    95.000        0.08720        0.11352
    99.000        0.08344        0.11802
    99.900        0.07921        0.12336
    99.990        0.07577        0.12795
    99.999        0.07282        0.13206
*/


