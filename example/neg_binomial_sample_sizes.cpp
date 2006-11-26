// Copyright Paul A. Bristow 2006
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
using std::cout;
using std::endl;
using std::fixed;
using std::right;
#include <iomanip>
using std::setprecision;
using std::setw;

#include <boost/math/distributions/negative_binomial.hpp>
using boost::math::negative_binomial;
  // RealType estimate_number_of_trials(
  // RealType k,     // number of failures, k >= 0.
  // RealType p,     // success fraction 0 <= p <= 1.
  // RealType probability) // probability threshold 0 <= p <= 0.
 
void estimate_number_of_trials(double failures, double p)
{
   // trials = number of trials
   // failures = number of failures before achieving required success(es).
   // p         = success ratio.
   // successes = required number of successes.
   //
   // Calculate how many trials we can have to ensure the
   // required number of successes DOES exceed "successes".

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
      // estimate_number_of_trials
      << setw(6) << right<< int(negative_binomial::estimate_number_of_trials(failures, p, alpha[i]/2))  << endl;
   }
   cout << endl;
} // void estimate_number_of_trials(double fails, double p)


int main()
{
	 estimate_number_of_trials(5, 0.5);
	 estimate_number_of_trials(50, 0.5);
	 estimate_number_of_trials(500, 0.5);
	 estimate_number_of_trials(50, 0.1);
	 estimate_number_of_trials(500, 0.1);
	 estimate_number_of_trials(5, 0.9);
	 estimate_number_of_trials(10-5, 0.4); // See Evans example in Wikipedia.
   return 0;
} // int main()


/*

------ Build started: Project: neg_binomial_sample_sizes, Configuration: Debug Win32 ------
Compiling...
neg_binomial_sample_sizes.cpp
Linking...
Autorun "i:\boost-06-05-03-1300\libs\math\test\Math_test\debug\neg_binomial_sample_sizes.exe"
Number of failures = 5,   Success fraction = 50%
____________________________
Confidence        Max Number
 Value (%)        Of Trials 
____________________________
    50.000          13
    75.000          15
    90.000          17
    95.000          19
    99.000          23
    99.900          27
    99.990          32
    99.999          36
Number of failures = 50.000,   Success fraction = 50.000%
____________________________
Confidence        Max Number
 Value (%)        Of Trials 
____________________________
    50.000         108
    75.000         113
    90.000         118
    95.000         122
    99.000         130
    99.900         139
    99.990         147
    99.999         155
Number of failures = 500.000,   Success fraction = 50.000%
____________________________
Confidence        Max Number
 Value (%)        Of Trials 
____________________________
    50.000        1022
    75.000        1038
    90.000        1054
    95.000        1064
    99.000        1085
    99.900        1110
    99.990        1131
    99.999        1150
Number of failures = 50.000,   Success fraction = 10.000%
____________________________
Confidence        Max Number
 Value (%)        Of Trials 
____________________________
    50.000         553
    75.000         588
    90.000         626
    95.000         651
    99.000         701
    99.900         763
    99.990         818
    99.999         869
Number of failures = 500.000,   Success fraction = 10.000%
____________________________
Confidence        Max Number
 Value (%)        Of Trials 
____________________________
    50.000        5150
    75.000        5254
    90.000        5364
    95.000        5434
    99.000        5574
    99.900        5739
    99.990        5880
    99.999        6006
Number of failures = 5.000,   Success fraction = 90.000%
____________________________
Confidence        Max Number
 Value (%)        Of Trials 
____________________________
    50.000           6
    75.000           7
    90.000           7
    95.000           8
    99.000           9
    99.900          10
    99.990          12
    99.999          13
Number of failures = 5.000,   Success fraction = 40.000%
____________________________
Confidence        Max Number
 Value (%)        Of Trials 
____________________________
    50.000          17
    75.000          20
    90.000          23
    95.000          25
    99.000          30
    99.900          36
    99.990          42
    99.999          48
Build Time 0:03
Build log was saved at "file://i:\boost-06-05-03-1300\libs\math\test\Math_test\neg_binomial_sample_sizes\Debug\BuildLog.htm"
neg_binomial_sample_sizes - 0 error(s), 0 warning(s)
========== Build: 1 succeeded, 0 failed, 0 up-to-date, 0 skipped ==========


*/
