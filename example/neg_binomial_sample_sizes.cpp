// Copyright Paul A. Bristow 2007
// Copyright John Maddock 2006
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/distributions/negative_binomial.hpp>
using boost::math::negative_binomial;
  // RealType find_number_of_trials(
  // RealType k,     // number of failures, k >= 0.
  // RealType p,     // success fraction 0 <= p <= 1.
  // RealType probability) // probability threshold 0 <= p <= 0.

#include <iostream>
using std::cout;
using std::endl;
using std::fixed;
using std::right;
#include <iomanip>
using std::setprecision;
using std::setw; 

void find_number_of_trials(double failures, double p)
{
   // trials = number of trials
   // failures = number of failures before achieving required success(es).
   // p         = success ratio.
   //
   // Calculate how many trials we need to ensure the
   // required number of failures DOES exceed "failures".

   // Define a table of significance levels:
   double alpha[] = { 0.5, 0.25, 0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001 };

  cout << "\n""Target number of failures = " << failures;
  cout << ",   Success fraction = " << 100 * p << "%" << endl;
   // Print table header:
   cout << "\n\n"
           "____________________________\n"
           "Confidence        Min Number\n"
           " Value (%)        Of Trials \n"
           "____________________________\n";
   // Now print out the data for the table rows.
  for(unsigned i = 0; i < sizeof(alpha)/sizeof(alpha[0]); ++i)
   { // Confidence values %:
      cout << fixed << setprecision(3) << setw(10) << right << 100 * (1-alpha[i]) << "      "
      // find_minimum_number_of_trials
      << setw(6) << right << (int)ceil(negative_binomial::find_minimum_number_of_trials(failures, p, alpha[i])) << endl;
   }
   cout << endl;
} // void find_number_of_trials(double fails, double p)


int main()
{
	 find_number_of_trials(5, 0.5);
	 find_number_of_trials(50, 0.5);
	 find_number_of_trials(500, 0.5);
	 find_number_of_trials(50, 0.1);
	 find_number_of_trials(500, 0.1);
	 find_number_of_trials(5, 0.9);
	 find_number_of_trials(10-5, 0.4); // See Evans example in Wikipedia.
    return 0;
} // int main()


/*

Output is:

Autorun "i:\boost-06-05-03-1300\libs\math\test\Math_test\debug\neg_binomial_sample_sizes.exe"
Target number of failures = 5,   Success fraction = 50%
____________________________
Confidence        Min Number
 Value (%)        Of Trials 
____________________________
    50.000          11
    75.000          14
    90.000          17
    95.000          18
    99.000          22
    99.900          27
    99.990          31
    99.999          36
Target number of failures = 50.000,   Success fraction = 50.000%
____________________________
Confidence        Min Number
 Value (%)        Of Trials 
____________________________
    50.000         101
    75.000         109
    90.000         115
    95.000         119
    99.000         128
    99.900         137
    99.990         146
    99.999         154
Target number of failures = 500.000,   Success fraction = 50.000%
____________________________
Confidence        Min Number
 Value (%)        Of Trials 
____________________________
    50.000        1001
    75.000        1023
    90.000        1043
    95.000        1055
    99.000        1078
    99.900        1104
    99.990        1126
    99.999        1146
Target number of failures = 50.000,   Success fraction = 10.000%
____________________________
Confidence        Min Number
 Value (%)        Of Trials 
____________________________
    50.000          56
    75.000          58
    90.000          60
    95.000          61
    99.000          63
    99.900          66
    99.990          68
    99.999          71
Target number of failures = 500.000,   Success fraction = 10.000%
____________________________
Confidence        Min Number
 Value (%)        Of Trials 
____________________________
    50.000         556
    75.000         562
    90.000         567
    95.000         570
    99.000         576
    99.900         583
    99.990         588
    99.999         594
Target number of failures = 5.000,   Success fraction = 90.000%
____________________________
Confidence        Min Number
 Value (%)        Of Trials 
____________________________
    50.000          57
    75.000          73
    90.000          91
    95.000         103
    99.000         127
    99.900         159
    99.990         189
    99.999         217
Target number of failures = 5.000,   Success fraction = 40.000%
____________________________
Confidence        Min Number
 Value (%)        Of Trials 
____________________________
    50.000          10
    75.000          11
    90.000          13
    95.000          15
    99.000          18
    99.900          21
    99.990          25
    99.999          28


*/
