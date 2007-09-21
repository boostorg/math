// Copyright John Maddock 2006, 2007
// Copyright Paul A. Bristow 2007

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <iomanip>
#include <boost/math/distributions/chi_squared.hpp>

void confidence_limits_on_std_deviation(
        double Sd,    // Sample Standard Deviation
        unsigned N)   // Sample size
{
   //
   // Calculate confidence intervals for the standard deviation.
   // For example if we set the confidence limit to
   // 0.95, we know that if we repeat the sampling
   // 100 times, then we expect that the true standard deviation
   // will be between out limits on 95 occations.
   // Note: this is not the same as saying a 95%
   // confidence interval means that there is a 95%
   // probability that the interval contains the true standard deviation.
   // The interval computed from a given sample either
   // contains the true standard deviation or it does not.
   // See http://www.itl.nist.gov/div898/handbook/eda/section3/eda358.htm

   using namespace std;
   using namespace boost::math;

   // Print out general info:
   cout <<
      "________________________________________________\n"
      "2-Sided Confidence Limits For Standard Deviation\n"
      "________________________________________________\n\n";
   cout << setprecision(7);
   cout << setw(40) << left << "Number of Observations" << "=  " << N << "\n";
   cout << setw(40) << left << "Standard Deviation" << "=  " << Sd << "\n";
   //
   // Define a table of significance/risk levels:
   //
   double alpha[] = { 0.5, 0.25, 0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001 };
   //
   // Start by declaring the distribution we'll need:
   //
   chi_squared dist(N - 1);
   //
   // Print table header:
   //
   cout << "\n\n"
           "_____________________________________________\n"
           "Confidence          Lower          Upper\n"
           " Value (%)          Limit          Limit\n"
           "_____________________________________________\n";
   //
   // Now print out the data for the table rows.
   //
   for(unsigned i = 0; i < sizeof(alpha)/sizeof(alpha[0]); ++i)
   {
      // Confidence value:
      cout << fixed << setprecision(3) << setw(10) << right << 100 * (1-alpha[i]);
      // Calculate limits:
      double lower_limit = sqrt((N - 1) * Sd * Sd / quantile(complement(dist, alpha[i] / 2)));
      double upper_limit = sqrt((N - 1) * Sd * Sd / quantile(dist, alpha[i] / 2));
      // Print Limits:
      cout << fixed << setprecision(5) << setw(15) << right << lower_limit;
      cout << fixed << setprecision(5) << setw(15) << right << upper_limit << endl;
   }
   cout << endl;
}

void chi_squared_test(
       double Sd,     // Sample std deviation
       double D,      // True std deviation
       unsigned N,    // Sample size
       double alpha)  // Significance level
{
   //
   // A Chi Squared test applied to a single set of data.
   // We are testing the null hypothesis that the true
   // standard deviation of the sample is D, and that any variation is down
   // to chance.  We can also test the alternative hypothesis
   // that any difference is not down to chance.
   // See http://www.itl.nist.gov/div898/handbook/eda/section3/eda358.htm
   //
   using namespace std;
   using namespace boost::math;

   // Print header:
   cout <<
      "______________________________________________\n"
      "Chi Squared test for sample standard deviation\n"
      "______________________________________________\n\n";
   cout << setprecision(5);
   cout << setw(55) << left << "Number of Observations" << "=  " << N << "\n";
   cout << setw(55) << left << "Sample Standard Deviation" << "=  " << Sd << "\n";
   cout << setw(55) << left << "Expected True Standard Deviation" << "=  " << D << "\n\n";
   //
   // Now we can calculate and output some stats:
   //
   // test-statistic:
   double t_stat = (N - 1) * (Sd / D) * (Sd / D);
   cout << setw(55) << left << "Test Statistic" << "=  " << t_stat << "\n";
   //
   // Finally define our distribution, and get the probability:
   //
   chi_squared dist(N - 1);
   double p = cdf(dist, t_stat);
   cout << setw(55) << left << "CDF of test statistic: " << "=  "
      << setprecision(3) << scientific << p << "\n";
   double ucv = quantile(complement(dist, alpha));
   double ucv2 = quantile(complement(dist, alpha / 2));
   double lcv = quantile(dist, alpha);
   double lcv2 = quantile(dist, alpha / 2);
   cout << setw(55) << left << "Upper Critical Value at alpha: " << "=  "
      << setprecision(3) << scientific << ucv << "\n";
   cout << setw(55) << left << "Upper Critical Value at alpha/2: " << "=  "
      << setprecision(3) << scientific << ucv2 << "\n";
   cout << setw(55) << left << "Lower Critical Value at alpha: " << "=  "
      << setprecision(3) << scientific << lcv << "\n";
   cout << setw(55) << left << "Lower Critical Value at alpha/2: " << "=  "
      << setprecision(3) << scientific << lcv2 << "\n\n";
   //
   // Finally print out results of alternative hypothesis:
   //
   cout << setw(55) << left <<
      "Results for Alternative Hypothesis and alpha" << "=  "
      << setprecision(4) << fixed << alpha << "\n\n";
   cout << "Alternative Hypothesis              Conclusion\n";
   cout << "Standard Deviation != " << setprecision(3) << fixed << D << "            ";
   if((ucv2 < t_stat) || (lcv2 > t_stat))
      cout << "NOT REJECTED\n";
   else
      cout << "REJECTED\n";
   cout << "Standard Deviation  < " << setprecision(3) << fixed << D << "            ";
   if(lcv > t_stat)
      cout << "NOT REJECTED\n";
   else
      cout << "REJECTED\n";
   cout << "Standard Deviation  > " << setprecision(3) << fixed << D << "            ";
   if(ucv < t_stat)
      cout << "NOT REJECTED\n";
   else
      cout << "REJECTED\n";
   cout << endl << endl;
}

void chi_squared_sample_sized(
        double diff,      // difference from variance to detect
        double variance)  // true variance
{
   using namespace std;
   using namespace boost::math;

   try
   {

   // Print out general info:
   cout <<
      "_____________________________________________________________\n"
      "Estimated sample sizes required for various confidence levels\n"
      "_____________________________________________________________\n\n";
   cout << setprecision(5);
   cout << setw(40) << left << "True Variance" << "=  " << variance << "\n";
   cout << setw(40) << left << "Difference to detect" << "=  " << diff << "\n";
   //
   // Define a table of significance levels:
   //
   double alpha[] = { 0.5, 0.33333333333333333333333, 0.25, 0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001 };
   //
   // Print table header:
   //
   cout << "\n\n"
           "_______________________________________________________________\n"
           "Confidence       Estimated          Estimated\n"
           " Value (%)      Sample Size        Sample Size\n"
           "                (lower one-         (upper one-\n"
           "                 sided test)        sided test)\n"
           "_______________________________________________________________\n";
   //
   // Now print out the data for the table rows.
   //
   for(unsigned i = 0; i < sizeof(alpha)/sizeof(alpha[0]); ++i)
   {
      // Confidence value:
      cout << fixed << setprecision(3) << setw(10) << right << 100 * (1-alpha[i]);
      // calculate df for a lower single-sided test:
      double df = chi_squared::find_degrees_of_freedom(
         -diff, alpha[i], alpha[i], variance);
      // convert to sample size:
      double size = ceil(df) + 1;
      // Print size:
      cout << fixed << setprecision(0) << setw(16) << right << size;
      // calculate df for an upper single-sided test:
      df = chi_squared::find_degrees_of_freedom(
         diff, alpha[i], alpha[i], variance);
      // convert to sample size:
      size = ceil(df) + 1;
      // Print size:
      cout << fixed << setprecision(0) << setw(16) << right << size << endl;
   }
   cout << endl;
   }
  catch(const std::exception& e)
  { // Always useful to include try & catch blocks because default policies
    // are to throw exceptions on arguments that cause errors like underflow, overflow.
    // Lacking try & catch blocks, the program will abort without a message below,
    // which may give some helpful clues as to the cause of the exception.
    std::cout <<
      "\n""Message from thrown exception was:\n   " << e.what() << std::endl;
  }
}

//Message from thrown exception was:  for alpha 0.5
//   Error in function boost::math::tools::bracket_and_solve_root<double>:
//Unable to bracket root, last nearest value was 1.2446030555722283e-058


int main()
{
   //
   // Run tests for Gear data
   // see http://www.itl.nist.gov/div898/handbook/eda/section3/eda3581.htm
   // Tests measurements of gear diameter.
   //
   confidence_limits_on_std_deviation(0.6278908E-02, 100);
   chi_squared_test(0.6278908E-02, 0.1, 100, 0.05);
   chi_squared_sample_sized(0.1 - 0.6278908E-02, 0.1);
   //
   // Run tests for silicon wafer fabrication data.
   // see http://www.itl.nist.gov/div898/handbook/prc/section2/prc23.htm
   // A supplier of 100 ohm.cm silicon wafers claims that his fabrication
   // process can produce wafers with sufficient consistency so that the
   // standard deviation of resistivity for the lot does not exceed
   // 10 ohm.cm. A sample of N = 10 wafers taken from the lot has a
   // standard deviation of 13.97 ohm.cm
   //
   confidence_limits_on_std_deviation(13.97, 10);
   chi_squared_test(13.97, 10.0, 10, 0.05);
   chi_squared_sample_sized(13.97 * 13.97 - 100, 100);
   chi_squared_sample_sized(55, 100);
   chi_squared_sample_sized(1, 100);
   return 0;
}

/*

Output:

________________________________________________
2-Sided Confidence Limits For Standard Deviation
________________________________________________

Number of Observations                  =  100
Standard Deviation                      =  0.006278908


_____________________________________________
Confidence          Lower          Upper
 Value (%)          Limit          Limit
_____________________________________________
    50.000        0.00601        0.00662
    75.000        0.00582        0.00685
    90.000        0.00563        0.00712
    95.000        0.00551        0.00729
    99.000        0.00530        0.00766
    99.900        0.00507        0.00812
    99.990        0.00489        0.00855
    99.999        0.00474        0.00895

______________________________________________
Chi Squared test for sample standard deviation
______________________________________________

Number of Observations                                 =  100
Sample Standard Deviation                              =  0.00628
Expected True Standard Deviation                       =  0.10000

Test Statistic                                         =  0.39030
CDF of test statistic:                                 =  1.438e-099
Upper Critical Value at alpha:                         =  1.232e+002
Upper Critical Value at alpha/2:                       =  1.284e+002
Lower Critical Value at alpha:                         =  7.705e+001
Lower Critical Value at alpha/2:                       =  7.336e+001

Results for Alternative Hypothesis and alpha           =  0.0500

Alternative Hypothesis              Conclusion
Standard Deviation != 0.100            NOT REJECTED
Standard Deviation  < 0.100            NOT REJECTED
Standard Deviation  > 0.100            REJECTED


_____________________________________________________________
Estimated sample sizes required for various confidence levels
_____________________________________________________________

True Variance                           =  0.10000
Difference to detect                    =  0.09372


_______________________________________________________________
Confidence       Estimated          Estimated
 Value (%)      Sample Size        Sample Size
                (lower one         (upper one
                 sided test)        sided test)
_______________________________________________________________
    50.000               2               2
    75.000               2              10
    90.000               4              32
    95.000               5              52
    99.000               8             102
    99.900              13             178
    99.990              18             257
    99.999              23             337

________________________________________________
2-Sided Confidence Limits For Standard Deviation
________________________________________________

Number of Observations                  =  10
Standard Deviation                      =  13.9700000


_____________________________________________
Confidence          Lower          Upper
 Value (%)          Limit          Limit
_____________________________________________
    50.000       12.41880       17.25579
    75.000       11.23084       19.74131
    90.000       10.18898       22.98341
    95.000        9.60906       25.50377
    99.000        8.62898       31.81825
    99.900        7.69466       42.51593
    99.990        7.04085       55.93352
    99.999        6.54517       73.00132

______________________________________________
Chi Squared test for sample standard deviation
______________________________________________

Number of Observations                                 =  10
Sample Standard Deviation                              =  13.97000
Expected True Standard Deviation                       =  10.00000

Test Statistic                                         =  17.56448
CDF of test statistic:                                 =  9.594e-001
Upper Critical Value at alpha:                         =  1.692e+001
Upper Critical Value at alpha/2:                       =  1.902e+001
Lower Critical Value at alpha:                         =  3.325e+000
Lower Critical Value at alpha/2:                       =  2.700e+000

Results for Alternative Hypothesis and alpha           =  0.0500

Alternative Hypothesis              Conclusion
Standard Deviation != 10.000            REJECTED
Standard Deviation  < 10.000            REJECTED
Standard Deviation  > 10.000            NOT REJECTED


_____________________________________________________________
Estimated sample sizes required for various confidence levels
_____________________________________________________________

True Variance                           =  100.00000
Difference to detect                    =  95.16090


_______________________________________________________________
Confidence       Estimated          Estimated
 Value (%)      Sample Size        Sample Size
                (lower one         (upper one
                 sided test)        sided test)
_______________________________________________________________
    50.000               2               2
    75.000               2              10
    90.000               4              32
    95.000               5              51
    99.000               7              99
    99.900              11             174
    99.990              15             251
    99.999              20             330

_____________________________________________________________
Estimated sample sizes required for various confidence levels
_____________________________________________________________

True Variance                           =  100.00000
Difference to detect                    =  55.00000


_______________________________________________________________
Confidence       Estimated          Estimated
 Value (%)      Sample Size        Sample Size
                (lower one         (upper one
                 sided test)        sided test)
_______________________________________________________________
    50.000               2               2
    75.000               8              21
    90.000              23              71
    95.000              36             115
    99.000              71             228
    99.900             123             401
    99.990             177             580
    99.999             232             762

_____________________________________________________________
Estimated sample sizes required for various confidence levels
_____________________________________________________________

True Variance                           =  100.00000
Difference to detect                    =  1.00000


_______________________________________________________________
Confidence       Estimated          Estimated
 Value (%)      Sample Size        Sample Size
                (lower one         (upper one
                 sided test)        sided test)
_______________________________________________________________
    50.000               2               2
    75.000           36033           36761
    90.000          130079          132707
    95.000          214283          218612
    99.000          428628          437287
    99.900          756333          771612
    99.990         1095435         1117564
    99.999         1440608         1469711

*/

