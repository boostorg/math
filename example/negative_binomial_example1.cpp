// negative_binomial_example1.cpp

// Copyright Paul A. Bristow 2007.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// Example 1 of using negative_binomial distribution.

//[negative_binomial_eg1_1

/*`
First we need to #define macros to control the error and discrete handling policies.
For this simple example, we want to avoid throwing
an exception (the default policy) and just return infinity.
We want to treat the distribution as if it was continuous, so we choose a policy of real,
rather than the default policy integer_round_outwards.
*/
#define BOOST_MATH_OVERFLOW_ERROR_POLICY ignore_error
#define BOOST_MATH_DISCRETE_QUANTILE_POLICY real

/*`After that (vital to #include *after* #defines), we need some includes
to provide easy access to the negative binomial distribution,
and some std library iostream, of course.
*/

#include <boost/math/distributions/negative_binomial.hpp>
// for negative_binomial_distribution
	using boost::math::negative_binomial; // typedef provides default type is double.
  using  ::boost::math::cdf;
  using  ::boost::math::pdf; // Probability mass of negative_binomial.
  using  ::boost::math::quantile;

#include <iostream>
  using std::cout;
  using std::endl;
  using std::noshowpoint;
  using std::fixed;
  using std::right;
  using std::left;
#include <iomanip>
  using std::setprecision;
  using std::setw; 

#include <limits>
  using std::numeric_limits;

int main()
{
  cout << "Example 1 using the Negative Binomial Distribution.";
  #if defined(__FILE__) && defined(__TIMESTAMP__)
	  cout << "  " << __FILE__ << ' ' << __TIMESTAMP__ << "\n";
  #endif
  cout << endl;
  cout.precision(5); 
  // None of the values calculated have a useful accuracy as great this, but
  // INF shows wrongly with < 5 !
  // https://connect.microsoft.com/VisualStudio/feedback/ViewFeedback.aspx?FeedbackID=240227
  try
  {
    double sales_quota = 5; // Pat's sales quota - successes (r).
    double success_fraction = 0.4; // success_fraction (p) - so failure_fraction is 0.6.
    negative_binomial nb(sales_quota, success_fraction); // double by default.
    int all_houses = 30; // The number of houses on the estate.

    cout <<"Selling candy bars - using the negative binomial distribution." 
      << "\nby Dr. Diane Evans,"
      "\nProfessor of Mathematics at Rose-Hulman Institute of Technology,"
      << "\nsee http://en.wikipedia.org/wiki/Negative_binomial_distribution\n"
      << endl;
    cout << "Pat has a sales per house success rate of " << success_fraction
      << ".\nTherefore he would, on average, sell " << nb.success_fraction() * 100
      << " bars after trying 100 houses." << endl;

    cout << "With a success rate of " << nb.success_fraction() 
      << ", he might expect, on average,\n"
        " to need to visit about " << success_fraction * all_houses
        << " houses in order to sell all " << nb.successes() << " candy bars. " << endl;
/*`
To finish on or before the 8th house, Pat must finish at the 5th, 6th, 7th or 8th house.
(Obviously he could not finish on fewer than 5 houses because he must sell 5 candy bars.
So the 5th house is the first that he could possibly finish on).
The probability that he will finish on EXACTLY on any house
is the Probability Density Function (pdf).
*/
    cout << "Probability that Pat finishes on the " << sales_quota << "th house is "
      << pdf(nb, 5 - sales_quota) << endl;
    cout << "Probability that Pat finishes on the 6th house is "
      << pdf(nb, 6 - sales_quota) << endl;
    cout << "Probability that Pat finishes on the 7th house is "
      << pdf(nb, 7 - sales_quota) << endl;
    cout << "Probability that Pat finishes on the 8th house is "
      << pdf(nb, 8 - sales_quota) << endl;

    // Sum of the probabilities for these houses is the Cumulative Distribution Function (cdf).
    cout << "Probability that Pat finishes on or before the 8th house is sum "
      "\n" << "pdf(sales_quota) + pdf(6) + pdf(7) + pdf(8) = "
      // Sum each of the mass/density probabilities for houses sales_quota = 5, 6, 7, & 8.
      << pdf(nb, 5 - sales_quota) // 0
      + pdf(nb, 6 - sales_quota) // 1
      + pdf(nb, 7 - sales_quota) // 2
      + pdf(nb, 8 - sales_quota) // 3
      << endl;

    // Or using the negative binomial *cumulative* distribution function
    // (cdf instead sum of the pdfs):
    cout << "\nProbability of selling his quota of " << sales_quota
      << " candy bars\non or before the " << 8 << "th house is "
      << cdf(nb, 8 - sales_quota) << endl;

    cout << "\nProbability that Pat finishes exactly on the 10th house is "
      << pdf(nb, 10 - sales_quota) << endl;
    cout << "\nProbability of selling his quota of " << sales_quota
      << " candy bars\non or before the " << 10 << "th house is "
      << cdf(nb, 10 - sales_quota) << endl;

    cout << "Probability that Pat finishes on the 11th house is "
      << pdf(nb, 11 - sales_quota) << endl;
    cout << "\nProbability of selling his quota of " << sales_quota
      << " candy bars\non or before the " << 11 << "th house is "
      << cdf(nb, 11 - sales_quota) << endl;

    cout << "Probability that Pat finishes on the 12th house is "
      << pdf(nb, 12 - sales_quota) << endl;
    cout << "\nProbability of selling his quota of " << sales_quota
      << " candy bars\non or before the " << 12 << "th house is "
      << cdf(nb, 12 - sales_quota) << endl;

/*`
Finally consider the risk of Pat not selling his quota of 5 bars
even after visiting all the houses.
Calculate the probability that he would sell on the (non-existent) last-plus-1 house.
*/
    cout << "\nProbability of selling his quota of " << sales_quota
      << " candy bars\non or before the " << all_houses << "th house is "
      << cdf(nb, all_houses - sales_quota) << endl;

/*`So the risk of failing even at the 31th (non-existent) house is 1 - this probability,
  ``1 - cdf(nb, all_houses - sales_quota``
But using this expression may cause serious inaccuracy,
so it would be much better to use the complement of the cdf:
*/
    cout << "\nProbability of failing to sell his quota of " << sales_quota
      << " candy bars\neven after visiting all " << all_houses << " houses is "
      << cdf(complement(nb, all_houses - sales_quota)) << endl;

    double p = cdf(nb, (8 - sales_quota)); 
    cout << "Probability of meeting sales quota on or before 8th house is "<< p << endl;
    // Probability of meeting sales quota on or before 8th house is 0.174
    cout << "If the confidence of meeting sales quota is " << p
        << ", then the finishing house is " << quantile(nb, p) + sales_quota << endl;
/*`
Also demanding absolute certainty that all 5 will be sold,
which implies an infinite number of trials.
(Of course, there are only 30 houses on the estate,
so he can't even be certain of selling his quota).
*/
    cout << "If the confidence of meeting sales quota is " << 1.
        << ", then the finishing house is " << quantile(nb, 1) + sales_quota << endl;
    //  1.#INF == infinity.

    cout << "If the confidence of meeting sales quota is " << 0.
        << ", then the finishing house is " << quantile(nb, 0.) + sales_quota << endl;

    cout << "If the confidence of meeting sales quota is " << 0.5
        << ", then the finishing house is " << quantile(nb, 0.5) + sales_quota << endl;

    cout << "If the confidence of meeting sales quota is " << 1 - 0.00151 // 30 th
        << ", then the finishing house is " << quantile(nb, 1 - 0.00151) + sales_quota << endl;

/*
If the opposite is true, we don't want to assume any confidence, then this is tantamount
to assuming that all the first sales_quota trials will be successful sales.
*/
    cout << "If confidence of meeting quota is zero (we assume all houses are successful sales)" 
      ", then finishing house is " << sales_quota << endl;

    double ps[] = {0., 0.001, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 0.999, 1.};
    // Confidence as fraction = 1-alpha, as percent =  100 * (1-alpha[i]) %
    for (int i = 0; i < sizeof(ps)/sizeof(ps[0]); i++)
    {
      cout << "If confidence of meeting quota is " << ps[i]
        << ", then finishing house is " << ceil(quantile(nb, ps[i])) + sales_quota
        << endl;
   }
   cout << "If we demand a confidence of meeting sales quota of unity"
     ", then we can never be certain of selling 5 bars, so the finishing house is infinite!"  << endl;

   cout << "\nProbability (%)    House for (last) " << sales_quota << " th sale." << endl;
   for (int i = (int)sales_quota; i < all_houses + 5; i++)
   {
     cout << left << setw(6) << cdf(nb, i - sales_quota) << "    " << setw(3) << i<< endl;
   }
   cout << endl;
  }
  catch(const std::exception& e)
  { // Since we have set an overflow policy of ignore_error,
    // an exception should never be thrown.
     std::cout << "\nMessage from thrown exception was:\n " << e.what() << std::endl;
  }
	return 0;
}  // int main()
//] [/ negative_binomial_eg1_1]


/*

Output is:

Example 1 using the Negative Binomial Distribution.  ..\..\..\..\..\..\boost-sandbox\math_toolkit\libs\math\example\negative_binomial_example1.cpp Sun Sep 16 11:52:52 2007
Selling candy bars - using the negative binomial distribution.
by Dr. Diane Evans,
Professor of Mathematics at Rose-Hulman Institute of Technology,
see http://en.wikipedia.org/wiki/Negative_binomial_distribution
Pat has a sales per house success rate of 0.4.
Therefore he would, on average, sell 40 bars after trying 100 houses.
With a success rate of 0.4, he might expect, on average,
 to need to visit about 12 houses in order to sell all 5 candy bars. 
Probability that Pat finishes on the 5th house is 0.01024
Probability that Pat finishes on the 6th house is 0.03072
Probability that Pat finishes on the 7th house is 0.055296
Probability that Pat finishes on the 8th house is 0.077414
Probability that Pat finishes on or before the 8th house is sum 
pdf(sales_quota) + pdf(6) + pdf(7) + pdf(8) = 0.17367
Probability of selling his quota of 5 candy bars
on or before the 8th house is 0.17367
Probability that Pat finishes exactly on the 10th house is 0.10033
Probability of selling his quota of 5 candy bars
on or before the 10th house is 0.3669
Probability that Pat finishes on the 11th house is 0.10033
Probability of selling his quota of 5 candy bars
on or before the 11th house is 0.46723
Probability that Pat finishes on the 12th house is 0.094596
Probability of selling his quota of 5 candy bars
on or before the 12th house is 0.56182
Probability of selling his quota of 5 candy bars
on or before the 30th house is 0.99849
Probability of failing to sell his quota of 5 candy bars
even after visiting all 30 houses is 0.0015101
Probability of failing to sell his quota of 5 candy bars
even after visiting all 30 houses is 0.0015101
Probability of meeting sales quota on or before 8th house is 0.17367
If the confidence of meeting sales quota is 0.17367, then the finishing house is 8
If the confidence of meeting sales quota is 1, then the finishing house is 1.#INF
If the confidence of meeting sales quota is 0, then the finishing house is 5
If the confidence of meeting sales quota is 0.5, then the finishing house is 11.337
If the confidence of meeting sales quota is 0.99849, then the finishing house is 30
If confidence of meeting quota is zero (we assume all houses are successful sales), then finishing house is 5
If confidence of meeting quota is 0, then finishing house is 5
If confidence of meeting quota is 0.001, then finishing house is 5
If confidence of meeting quota is 0.01, then finishing house is 5
If confidence of meeting quota is 0.05, then finishing house is 7
If confidence of meeting quota is 0.1, then finishing house is 8
If confidence of meeting quota is 0.5, then finishing house is 12
If confidence of meeting quota is 0.9, then finishing house is 18
If confidence of meeting quota is 0.95, then finishing house is 21
If confidence of meeting quota is 0.99, then finishing house is 25
If confidence of meeting quota is 0.999, then finishing house is 32
If confidence of meeting quota is 1, then finishing house is 1.#INF
If we demand a confidence of meeting sales quota of unity, then we can never be certain of selling 5 bars, so the finishing house is infinite!
Probability (%)    House for (last) 5 th sale.
0.01024    5  
0.04096    6  
0.096256    7  
0.17367    8  
0.26657    9  
0.3669    10 
0.46723    11 
0.56182    12 
0.64696    13 
0.72074    14 
0.78272    15 
0.83343    16 
0.874     17 
0.90583    18 
0.93039    19 
0.94905    20 
0.96304    21 
0.97342    22 
0.98103    23 
0.98655    24 
0.99053    25 
0.99337    26 
0.99539    27 
0.99681    28 
0.9978    29 
0.99849    30 
0.99897    31 
0.9993    32 
0.99952    33 
0.99968    34 


*/






