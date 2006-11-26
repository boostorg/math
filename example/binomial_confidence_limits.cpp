// Copyright John Maddock 2006
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
#include <iomanip>
#include <boost/math/distributions/binomial.hpp>

void confidence_limits_on_frequency(unsigned trials, unsigned successes)
{
   //
   // trials = Total number of trials.
   // successes = Total number of observed successes.
   //
   // Calculate confidence limits for an observed
   // frequency of occurance that follows a binomial
   // distribution.
   //
   using namespace std;
   using namespace boost::math;

   // Print out general info:
   cout <<
      "___________________________________________\n"
      "2-Sided Confidence Limits For Success Ratio\n"
      "___________________________________________\n\n";
   cout << setprecision(7);
   cout << setw(40) << left << "Number of Observations" << "=  " << trials << "\n";
   cout << setw(40) << left << "Number of successes" << "=  " << successes << "\n";
   cout << setw(40) << left << "Sample frequency of occurrence" << "=  " << double(successes) / trials << "\n";
   //
   // Define a table of significance levels:
   //
   double alpha[] = { 0.5, 0.25, 0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001 };
   //
   // Print table header:
   //
   cout << "\n\n"
           "___________________________________________\n"
           "Confidence        Lower          Upper\n"
           " Value (%)        Limit          Limit\n"
           "___________________________________________\n";
   //
   // Now print out the data for the table rows.
   //
   for(unsigned i = 0; i < sizeof(alpha)/sizeof(alpha[0]); ++i)
   {
      // Confidence value:
      cout << fixed << setprecision(3) << setw(10) << right << 100 * (1-alpha[i]);
      // calculate bounds:
      double l = binomial_distribution<>::estimate_lower_bound_on_p(trials, successes, alpha[i]/2);
      double u = binomial_distribution<>::estimate_upper_bound_on_p(trials, successes, alpha[i]/2);
      // Print Limits:
      cout << fixed << setprecision(5) << setw(15) << right << l;
      cout << fixed << setprecision(5) << setw(15) << right << u << endl;
   }
   cout << endl;
} // void confidence_limits_on_frequency()

int main()
{
   confidence_limits_on_frequency(20, 4);
   confidence_limits_on_frequency(200, 40);
   confidence_limits_on_frequency(2000, 400);

   return 0;
} // int main()

/*

------ Build started: Project: binomial_confidence_limits, Configuration: Debug Win32 ------
Compiling...
binomial_confidence_limits.cpp
Linking...
Autorun "i:\boost-06-05-03-1300\libs\math\test\Math_test\debug\binomial_confidence_limits.exe"
___________________________________________
2-Sided Confidence Limits For Success Ratio
___________________________________________

Number of Observations                  =  20
Number of successes                     =  4
Sample frequency of occurrence          =  0.2


___________________________________________
Confidence        Lower          Upper
 Value (%)        Limit          Limit
___________________________________________
    50.000        0.12840        0.29588
    75.000        0.09775        0.34633
    90.000        0.07135        0.40103
    95.000        0.05733        0.43661
    99.000        0.03576        0.50661
    99.900        0.01905        0.58632
    99.990        0.01042        0.64997
    99.999        0.00577        0.70216

___________________________________________
2-Sided Confidence Limits For Success Ratio
___________________________________________

Number of Observations                  =  200
Number of successes                     =  40
Sample frequency of occurrence          =  0.2000000


___________________________________________
Confidence        Lower          Upper
 Value (%)        Limit          Limit
___________________________________________
    50.000        0.17949        0.22259
    75.000        0.16701        0.23693
    90.000        0.15455        0.25225
    95.000        0.14689        0.26223
    99.000        0.13257        0.28218
    99.900        0.11703        0.30601
    99.990        0.10489        0.32652
    99.999        0.09492        0.34485

___________________________________________
2-Sided Confidence Limits For Success Ratio
___________________________________________

Number of Observations                  =  2000
Number of successes                     =  400
Sample frequency of occurrence          =  0.2000000


___________________________________________
Confidence        Lower          Upper
 Value (%)        Limit          Limit
___________________________________________
    50.000        0.19382        0.20638
    75.000        0.18965        0.21072
    90.000        0.18537        0.21528
    95.000        0.18267        0.21821
    99.000        0.17745        0.22400
    99.900        0.17150        0.23079
    99.990        0.16658        0.23657
    99.999        0.16233        0.24169

*/



