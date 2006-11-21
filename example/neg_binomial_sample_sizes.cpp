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
#include <iomanip>
#include <boost/math/distributions/negative_binomial.hpp>

void estimate_max_sample_size(double p, unsigned successes)
{
   // trials = number of trials
   // p         = success ratio.
   // successes = Total number of observed successes.
   //
   // Calculate how many trials we can have to ensure the
   // maximum number of successes does not exceed "successes".
   // A typical use would be failure analysis, where you want
   // zero or fewer "successes" with some probability.
   //
   using namespace std;
   using namespace boost::math;

   // Print out general info:
   cout <<
      "________________________\n"
      "Maximum Number of Trials\n"
      "________________________\n\n";
   cout << setprecision(7);
   cout << setw(40) << left << "Success ratio" << "=  " << p << "\n";
   cout << setw(40) << left << "Maximum Number of \"successes\" permitted" << "=  " << successes << "\n";
   //
   // Define a table of confidence intervals:
   //
   double alpha[] = { 0.5, 0.25, 0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001 };
   //
   // Print table header:
   //
   cout << "\n\n"
           "____________________________\n"
           "Confidence        Max Number\n"
           " Value (%)        Of Trials \n"
           "____________________________\n";
   //
   // Now print out the data for the table rows.
   //
   for(unsigned i = 0; i < sizeof(alpha)/sizeof(alpha[0]); ++i)
   {
      // Confidence value:
      cout << fixed << setprecision(3) << setw(10) << right << 100 * (1-alpha[i]);
      // calculate trials:
      double t = negative_binomial::estimate_number_of_trials(complement(successes, p, alpha[i]));
      t = floor(t);
      // Print Trials:
      cout << fixed << setprecision(0) << setw(15) << right << t << endl;
   }
   cout << endl;
} // void estimate_max_sample_size(double p, unsigned successes)

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
   estimate_max_sample_size(1.0/1000, 0);
   estimate_max_sample_size(1.0/10000, 0);
   estimate_max_sample_size(1.0/100000, 0);
   estimate_max_sample_size(1.0/1000000, 0);

	 estimate_number_of_trials(5, 0.5);
	 estimate_number_of_trials(50, 0.5);
	 estimate_number_of_trials(500, 0.5);
	 estimate_number_of_trials(50, 0.1);
	 estimate_number_of_trials(500, 0.1);
	 estimate_number_of_trials(5, 0.9);

   return 0;
} // int main()


/*

Autorun "i:\boost-06-05-03-1300\libs\math\test\Math_test\debug\neg_binomial_sample_sizes.exe"
________________________
Maximum Number of Trials
________________________
Success ratio                           =  0.001
Maximum Number of "successes" permitted =  0
____________________________
Confidence        Max Number
 Value (%)        Of Trials 
____________________________
    50.000              1
    75.000              1
    90.000              1
    95.000              1
    99.000              1
    99.900              1
    99.990              1
    99.999              1
________________________
Maximum Number of Trials
________________________
Success ratio                           =  0.0001000
Maximum Number of "successes" permitted =  0
____________________________
Confidence        Max Number
 Value (%)        Of Trials 
____________________________
    50.000              1
    75.000              1
    90.000              1
    95.000              1
    99.000              1
    99.900              1
    99.990              1
    99.999              1
________________________
Maximum Number of Trials
________________________
Success ratio                           =  0.0000100
Maximum Number of "successes" permitted =  0
____________________________
Confidence        Max Number
 Value (%)        Of Trials 
____________________________
    50.000              1
    75.000              1
    90.000              1
    95.000              1
    99.000              1
    99.900              1
    99.990              1
    99.999              1
________________________
Maximum Number of Trials
________________________
Success ratio                           =  0.0000010
Maximum Number of "successes" permitted =  0
____________________________
Confidence        Max Number
 Value (%)        Of Trials 
____________________________
    50.000              1
    75.000              1
    90.000              1
    95.000              1
    99.000              1
    99.900              1
    99.990              1
    99.999              1
Number of failures = 5,   Success ratio = 50%
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

*/
