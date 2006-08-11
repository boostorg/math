// (C) Copyright John Maddock 2006
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <iomanip>
#include <boost/math/distributions/students_t.hpp>

void confidence_limits_on_mean(double Sm, double Sd, unsigned Sn)
{
   //
   // Sm = Sample Mean.
   // Sd = Sample Standard Deviation.
   // Sn = Sample Size.
   //
   // Calculate confidence intervals for the mean.
   // For example if we set the confidence limit to
   // 0.95, we know that if we repeat the sampling
   // 100 times, then we expect that the true mean
   // will be between out limits on 95 occations.
   // Note: this is not the same as saying a 95%
   // confidence interval means that there is a 95%
   // probability that the interval contains the true mean.
   // The interval computed from a given sample either
   // contains the true mean or it does not.
   // See http://www.itl.nist.gov/div898/handbook/eda/section3/eda352.htm

   using namespace std;
   using namespace boost::math;

   // Print out general info:
   cout <<
      "__________________________________\n"
      "2-Sided Confidence Limits For Mean\n"
      "__________________________________\n\n";
   cout << setprecision(7);
   cout << setw(40) << left << "Number of Observations" << "=  " << Sn << "\n";
   cout << setw(40) << left << "Mean" << "=  " << Sm << "\n";
   cout << setw(40) << left << "Standard Deviation" << "=  " << Sd << "\n";
   //
   // Define a table of confidence intervals:
   //
   double alpha[] = { 0.5, 0.25, 0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001 };
   //
   // Start by declaring the distribution we'll need:
   //
   students_t dist(Sn - 1);
   //
   // Print table header:
   //
   cout << "\n\n"
           "_______________________________________________________________\n"
           "Confidence       T           Interval          Lower          Upper\n"
           " Value (%)     Value          Width            Limit          Limit\n"
           "_______________________________________________________________\n";
   //
   // Now print out the data for the table rows.
   //
   for(unsigned i = 0; i < sizeof(alpha)/sizeof(alpha[0]); ++i)
   {
      // Confidence value:
      cout << fixed << setprecision(3) << setw(10) << right << 100 * (1-alpha[i]);
      // calculate T:
      double T = quantile(complement(dist, alpha[i] / 2));
      // Print T:
      cout << fixed << setprecision(3) << setw(10) << right << T;
      // Calculate width of interval (one sided):
      double w = T * Sd / sqrt(double(Sn));
      // Print width:
      if(w < 0.01)
         cout << scientific << setprecision(3) << setw(17) << right << w;
      else
         cout << fixed << setprecision(3) << setw(17) << right << w;
      // Print Limits:
      cout << fixed << setprecision(5) << setw(15) << right << Sm - w;
      cout << fixed << setprecision(5) << setw(15) << right << Sm + w << endl;
   }
   cout << endl;
}

void single_sample_t_test(double M, double Sm, double Sd, unsigned Sn, double alpha)
{
   //
   // M = true mean.
   // Sm = Sample Mean.
   // Sd = Sample Standard Deviation.
   // Sn = Sample Size.
   // alpha = Confidence Level.
   //
   // A Students t test applied to a single set of data.
   // We are testing the null hypothesis that the true
   // mean of the sample is M, and that any variation is down
   // to chance.  We can also test the alternative hypothesis
   // that any difference is not down to chance.
   // See http://www.itl.nist.gov/div898/handbook/eda/section3/eda352.htm
   //
   using namespace std;
   using namespace boost::math;

   // Print header:
   cout <<
      "__________________________________\n"
      "Student t test for a single sample\n"
      "__________________________________\n\n";
   cout << setprecision(5);
   cout << setw(55) << left << "Number of Observations" << "=  " << Sn << "\n";
   cout << setw(55) << left << "Sample Mean" << "=  " << Sm << "\n";
   cout << setw(55) << left << "Sample Standard Deviation" << "=  " << Sd << "\n";
   cout << setw(55) << left << "Expected True Mean" << "=  " << M << "\n\n";
   //
   // Now we can calculate and output some stats:
   //
   // Difference in means:
   double diff = Sm - M;
   cout << setw(55) << left << "Sample Mean - Expected Test Mean" << "=  " << diff << "\n";
   // Degrees of freedom:
   unsigned v = Sn - 1;
   cout << setw(55) << left << "Degrees of Freedom" << "=  " << v << "\n";
   // t-statistic:
   double t_stat = diff * sqrt(double(Sn)) / Sd;
   cout << setw(55) << left << "T Statistic" << "=  " << t_stat << "\n";
   //
   // Finally define our distribution, and get the probability:
   //
   students_t dist(v);
   double q = cdf(complement(dist, fabs(t_stat)));
   cout << setw(55) << left << "Probability that difference is due to chance" << "=  "
      << setprecision(3) << scientific << q << "\n\n";
   //
   // Finally print out results of alternative hypothesis:
   //
   cout << setw(55) << left <<
      "Results for Alternative Hypothesis and alpha" << "=  "
      << setprecision(4) << fixed << alpha << "\n\n";
   cout << "Alternative Hypothesis     Conclusion\n";
   cout << "Mean != " << setprecision(3) << fixed << M << "            ";
   if(q < alpha)
      cout << "ACCEPTED\n";
   else
      cout << "REJECTED\n";
   cout << "Mean  < " << setprecision(3) << fixed << M << "            ";
   if(cdf(dist, t_stat) < alpha)
      cout << "ACCEPTED\n";
   else
      cout << "REJECTED\n";
   cout << "Mean  > " << setprecision(3) << fixed << M << "            ";
   if(cdf(complement(dist, t_stat)) < alpha)
      cout << "ACCEPTED\n";
   else
      cout << "REJECTED\n";
   cout << endl << endl;
}

void single_sample_estimate_df(double M, double Sm, double Sd)
{
   //
   // M = true mean.
   // Sm = Sample Mean.
   // Sd = Sample Standard Deviation.
   //
   using namespace std;
   using namespace boost::math;

   // Print out general info:
   cout <<
      "_____________________________________________________________\n"
      "Estimated sample sizes required for various confidence levels\n"
      "_____________________________________________________________\n\n";
   cout << setprecision(5);
   cout << setw(40) << left << "True Mean" << "=  " << M << "\n";
   cout << setw(40) << left << "Sample Mean" << "=  " << Sm << "\n";
   cout << setw(40) << left << "Sample Standard Deviation" << "=  " << Sd << "\n";
   //
   // Define a table of confidence intervals:
   //
   double alpha[] = { 0.5, 0.25, 0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001 };
   //
   // Print table header:
   //
   cout << "\n\n"
           "_______________________________________________________________\n"
           "Confidence      Estimated\n"
           " Value (%)     Sample Size\n"
           "_______________________________________________________________\n";
   //
   // Now print out the data for the table rows.
   //
   for(unsigned i = 0; i < sizeof(alpha)/sizeof(alpha[0]); ++i)
   {
      // Confidence value:
      cout << fixed << setprecision(3) << setw(10) << right << 100 * (1-alpha[i]);
      // calculate df:
      double df = students_t::estimate_degrees_of_freedom(complement(M, Sm, Sd, alpha[i]));
      // convert to sample size:
      double size = ceil(df) + 1;
      // Print size:
      cout << fixed << setprecision(0) << setw(16) << right << size << endl;
   }
   cout << endl;
}

int main()
{
   //
   // Run tests for Heat Flow Meter data
   // see http://www.itl.nist.gov/div898/handbook/eda/section4/eda428.htm
   // The data was collected while calibrating a heat flow meter
   // against a known value.
   //
   confidence_limits_on_mean(9.261460, 0.2278881e-01, 195);
   single_sample_t_test(5, 9.261460, 0.2278881e-01, 195, 0.05);
   single_sample_estimate_df(5, 9.261460, 0.2278881e-01);

   //
   // Data for this example from:
   // P.K.Hou, O. W. Lau & M.C. Wong, Analyst (1983) vol. 108, p 64.
   // from Statistics for Analytical Chemistry, 3rd ed. (1994), pp 54-55
   // J. C. Miller and J. N. Miller, Ellis Horwood ISBN 0 13 0309907
   //
   // Determination of mercury by cold-vapour atomic absorption,
   // the following values were obtained fusing a trusted
   // Standard Reference Material containing 38.9% mercury,
   // which we assume is correct or 'true'.
   //
   confidence_limits_on_mean(37.8, 0.964365, 3);
   // 95% test:
   single_sample_t_test(38.9, 37.8, 0.964365, 3, 0.05);
   // 90% test:
   single_sample_t_test(38.9, 37.8, 0.964365, 3, 0.1);
   // parameter estimate:
   single_sample_estimate_df(38.9, 37.8, 0.964365);
}
