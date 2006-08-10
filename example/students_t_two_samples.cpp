// (C) Copyright John Maddock 2006
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <iomanip>
#include <boost/math/dist/students_t.hpp>

void two_samples_t_test(
        double Sm1, 
        double Sd1, 
        unsigned Sn1, 
        double Sm2,
        double Sd2,
        unsigned Sn2,
        double alpha)
{
   //
   // Sm1 = Sample Mean 1.
   // Sd1 = Sample Standard Deviation 1.
   // Sn1 = Sample Size 1.
   // Sm2 = Sample Mean 2.
   // Sd2 = Sample Standard Deviation 2.
   // Sn2 = Sample Size 2.
   // alpha = Confidence Level.
   //
   // A Students t test applied to two sets of data.
   // We are testing the null hypothesis that the two
   // samples have the same mean and that any difference
   // if due to chance.  
   // See http://www.itl.nist.gov/div898/handbook/eda/section3/eda353.htm
   //
   using namespace std;
   using namespace boost::math;

   // Print header:
   cout << 
      "_________________________________________________\n"
      "Student t test for two samples (unequal variances)\n"
      "_________________________________________________\n\n";
   cout << setprecision(5);
   cout << setw(55) << left << "Number of Observations (Sample 1)" << "=  " << Sn1 << "\n";
   cout << setw(55) << left << "Sample 1 Mean" << "=  " << Sm1 << "\n";
   cout << setw(55) << left << "Sample 1 Standard Deviation" << "=  " << Sd1 << "\n";
   cout << setw(55) << left << "Number of Observations (Sample 2)" << "=  " << Sn2 << "\n";
   cout << setw(55) << left << "Sample 2 Mean" << "=  " << Sm2 << "\n";
   cout << setw(55) << left << "Sample 2 Standard Deviation" << "=  " << Sd2 << "\n";
   //
   // Now we can calculate and output some stats:
   //
   // Degrees of freedom:
   double v = Sd1 * Sd1 / Sn1 + Sd2 * Sd2 / Sn2;
   v *= v;
   double t1 = Sd1 * Sd1 / Sn1;
   t1 *= t1;
   t1 /=  (Sn1 - 1);
   double t2 = Sd2 * Sd2 / Sn2;
   t2 *= t2;
   t2 /= (Sn2 - 1);
   v /= (t1 + t2);
   cout << setw(55) << left << "Degrees of Freedom" << "=  " << v << "\n";
   // t-statistic:
   double t_stat = (Sm1 - Sm2) / sqrt(Sd1 * Sd1 / Sn1 + Sd2 * Sd2 / Sn2);
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
   cout << "Sample 1 Mean != Sample 2 Mean       " ;
   if(q < alpha)
      cout << "ACCEPTED\n";
   else
      cout << "REJECTED\n";
   cout << "Sample 1 Mean <  Sample 2 Mean       ";
   if(cdf(dist, t_stat) < alpha)
      cout << "ACCEPTED\n";
   else
      cout << "REJECTED\n";
   cout << "Sample 1 Mean >  Sample 2 Mean       ";
   if(cdf(complement(dist, t_stat)) < alpha)
      cout << "ACCEPTED\n";
   else
      cout << "REJECTED\n";
   cout << endl << endl;
}

void two_samples_t_test_equal_sd(
        double Sm1, 
        double Sd1, 
        unsigned Sn1, 
        double Sm2,
        double Sd2,
        unsigned Sn2,
        double alpha)
{
   //
   // Sm1 = Sample Mean 1.
   // Sd1 = Sample Standard Deviation 1.
   // Sn1 = Sample Size 1.
   // Sm2 = Sample Mean 2.
   // Sd2 = Sample Standard Deviation 2.
   // Sn2 = Sample Size 2.
   // alpha = Confidence Level.
   //
   // A Students t test applied to two sets of data.
   // We are testing the null hypothesis that the two
   // samples have the same mean and that any difference
   // if due to chance.  
   // See http://www.itl.nist.gov/div898/handbook/eda/section3/eda353.htm
   //
   using namespace std;
   using namespace boost::math;

   // Print header:
   cout << 
      "_______________________________________________\n"
      "Student t test for two samples (equal variances)\n"
      "_______________________________________________\n\n";
   cout << setprecision(5);
   cout << setw(55) << left << "Number of Observations (Sample 1)" << "=  " << Sn1 << "\n";
   cout << setw(55) << left << "Sample 1 Mean" << "=  " << Sm1 << "\n";
   cout << setw(55) << left << "Sample 1 Standard Deviation" << "=  " << Sd1 << "\n";
   cout << setw(55) << left << "Number of Observations (Sample 2)" << "=  " << Sn2 << "\n";
   cout << setw(55) << left << "Sample 2 Mean" << "=  " << Sm2 << "\n";
   cout << setw(55) << left << "Sample 2 Standard Deviation" << "=  " << Sd2 << "\n";
   //
   // Now we can calculate and output some stats:
   //
   // Degrees of freedom:
   double v = Sn1 + Sn2 - 2;
   cout << setw(55) << left << "Degrees of Freedom" << "=  " << v << "\n";
   // Pooled variance:
   double sp = sqrt(((Sn1-1) * Sd1 * Sd1 + (Sn2-1) * Sd2 * Sd2) / v);
   cout << setw(55) << left << "Pooled Standard Deviation" << "=  " << v << "\n";
   // t-statistic:
   double t_stat = (Sm1 - Sm2) / (sp * sqrt(1.0 / Sn1 + 1.0 / Sn2));
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
   cout << "Sample 1 Mean != Sample 2 Mean       " ;
   if(q < alpha)
      cout << "ACCEPTED\n";
   else
      cout << "REJECTED\n";
   cout << "Sample 1 Mean <  Sample 2 Mean       ";
   if(cdf(dist, t_stat) < alpha)
      cout << "ACCEPTED\n";
   else
      cout << "REJECTED\n";
   cout << "Sample 1 Mean >  Sample 2 Mean       ";
   if(cdf(complement(dist, t_stat)) < alpha)
      cout << "ACCEPTED\n";
   else
      cout << "REJECTED\n";
   cout << endl << endl;
}

void two_samples_estimate_df(double m1, double s1, unsigned n1, double m2, double s2)
{
   //
   // m1 = Sample 1 Mean.
   // s1 = Sample 1 Standard Deviation.
   // n1 = Sample 1 Size.
   // m2 = Sample 2 Mean.
   // s2 = Sample 2 Standard Deviation.
   // alpha = confidence level
   //
   using namespace std;
   using namespace boost::math;

   // Print out general info:
   cout << 
      "_____________________________________________________________\n"
      "Estimated sample sizes required for various confidence levels\n"
      "_____________________________________________________________\n\n";
   cout << setprecision(5);
   cout << setw(40) << left << "Sample 1 Mean" << "=  " << m1 << "\n";
   cout << setw(40) << left << "Sample 1 Standard Deviation" << "=  " << s1 << "\n";
   cout << setw(40) << left << "Sample 1 Size" << "=  " << n1 << "\n";
   cout << setw(40) << left << "Sample 2 Mean" << "=  " << m2 << "\n";
   cout << setw(40) << left << "Sample 2 Standard Deviation" << "=  " << s2 << "\n";
   //
   // Define a table of confidence intervals:
   //
   double alpha[] = { 0.5, 0.25, 0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001 };
   //
   // Print table header:
   //
   cout << "\n\n"
           "_______________________________________________________________________\n"
           "Confidence     Estimated Sample Size          Estimated Sample 2 Size\n"
           " Value (%)     (With Two Equal Sizes)        (With Fixed Sample 1 Size)\n"
           "_______________________________________________________________________\n";
   //
   // Now print out the data for the table rows.
   //
   for(unsigned i = 0; i < sizeof(alpha)/sizeof(alpha[0]); ++i)
   {
      // Confidence value:
      cout << fixed << setprecision(3) << setw(10) << right << 100 * (1-alpha[i]);
      // calculate df assuming equal sample sizes:
      double df = students_t::estimate_two_equal_degrees_of_freedom(
         complement(m1, s1, m2, s2, alpha[i]));
      // convert to sample size:
      double size = (ceil(df) + 2) / 2; 
      // Print size:
      cout << fixed << setprecision(0) << setw(28) << right << size;
      // calculate df with sample 1 fixed:
      df = students_t::estimate_two_unequal_degrees_of_freedom(
         complement(m1, s1, n1, m2, s2, alpha[i]));
      // convert to sample size:
      size = (ceil(df) + 2) - n1; 
      // Print size:
      cout << fixed << setprecision(0) << setw(28) << right << size << endl;
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
   two_samples_t_test(20.14458, 6.414700, 249, 30.48101, 6.107710, 79, 0.05);
   two_samples_t_test_equal_sd(20.14458, 6.414700, 249, 30.48101, 6.107710, 79, 0.05);
   two_samples_estimate_df(20.14458, 6.414700, 249, 30.48101, 6.107710);

   return 0;
}
