// negative_binomial_example1.cpp

// Copyright Paul A. Bristow 2007.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// Example 1 of using negative_binomial distribution.

// http://en.wikipedia.org/wiki/Negative_binomial_distribution

// (After a problem by Dr. Diane Evans,
// Professor of mathematics at Rose-Hulman Institute of Technology)

// Pat is required to sell candy bars to raise money for the 6th grade field trip.
// There are thirty houses in the neighborhood,
// and Pat is not supposed to return home until five candy bars have been sold.
// So the child goes door to door, selling candy bars.
// At each house, there is a 0.4 probability (40%) of selling one candy bar
// and a 0.6 probability (60%) of selling nothing.

// What is the probability mass/density function for selling the last (fifth) candy bar at the nth house?

// The Negative Binomial(r, p) distribution describes the probability of k failures
// and r successes in k+r Bernoulli(p) trials with success on the last trial.
// Selling five candy bars means getting five successes, so successes r = 5.
// The total number of trials n (in this case, houses) this takes is therefore
//    = sucesses + failures or k + r = k + 5.
// The random variable we are interested in is the number of houses k
// that must be visited to sell five candy bars,
// so we substitute k = n 5 into a NegBin(5, 0.4) mass/density function
// and obtain the following mass/density function of the distribution of houses (for n >= 5):
// Obviously, the best case is that Pat makes sales on all the first five houses.

// What is the probability that Pat finishes ON the tenth house?

//    f(10) = 0.1003290624, or about 1 in 10

// What is the probability that Pat finishes ON OR BEFORE reaching the eighth house?

// To finish on or before the eighth house,
// Pat must finish at the fifth, sixth, seventh, or eighth house.
// Sum those probabilities:

    // f(5) = 0.01024
    // f(6) = 0.03072
    // f(7) = 0.055296
    // f(8) = 0.0774144
    // sum {j=5 to 8} f(j) = 0.17367

// What is the probability that Pat exhausts all 30 houses in the neighborhood,
// and still doesn't sell the required 5 candy bars?

// 1 - sum{j=5 to 30} f(j) = 1 - incomplete beta (p = 0.4)(5, 30-5+1) =~ 1 - 0.99849 = 0.00151 = 0.15%.

// see also http://www.math.uah.edu/stat/bernoulli/Introduction.xhtml
// http://www.codecogs.com/pages/catgen.php?category=stats/dists/discrete

#include <boost/math/distributions/negative_binomial.hpp> // for negative_binomial_distribution
  using boost::math::negative_binomial_distribution;
	using boost::math::negative_binomial; // typedef provides default type is double.
  using  ::boost::math::cdf;
  using  ::boost::math::pdf; // Probability of negative_binomial.
  using  ::boost::math::quantile;

#include <iostream>
	using std::cout;
	using std::endl;
	using std::noshowpoint;
#include <iomanip>
	using std::setprecision;

int main()
{
	cout << "Example 1 using the Negative Binomial Distribution.";
  #if defined(__FILE__) && defined(__TIMESTAMP__)
  	cout << "  " << __FILE__ << ' ' << __TIMESTAMP__ << ' '<< _MSC_FULL_VER << "\n";
  #endif
	cout << endl;
  cout.precision(5); // NB INF shows wrongly with < 5 !
  // https://connect.microsoft.com/VisualStudio/feedback/ViewFeedback.aspx?FeedbackID=240227

  try
  {
    double sales_quota = 5; // Pat's sales quota - successes (r).
    double success_fraction = 0.4; // success_fraction (p) - so fail_fraction is 0.6.
    negative_binomial nb(sales_quota, success_fraction); // double by default.
    int all_houses = 30; // The number of houses on the estate.

    cout <<"Selling candy bars - an example of using the negative binomial distribution. " 
      << "\n""An example by Dr. Diane Evans,"
      "\n""Professor of Mathematics at Rose-Hulman Institute of Technology,"
      << "\n""see http://en.wikipedia.org/wiki/Negative_binomial_distribution""\n"
      << endl;
    cout << "Pat has a sales per house success rate of " << success_fraction
      << ".""\n""Therefore he would, on average, sell " << nb.success_fraction() * 100 << " bars after trying 100 houses." << endl;

    cout << "With a success rate of " << nb.success_fraction()  << ", he might expect, on average,""\n"
      " to need to visit about " << success_fraction * all_houses << " houses in order to sell all " << nb.successes() << " candy bars. " << endl;

    // To finish on or before the 8th house, Pat must finish at the 5th, 6th, 7th or 8th house.
    // (Obviously he could not finish on fewer than 5 houses because he must sell 5 candy bars.
    // so the 5th house is the first  that he could possibly finish on).
    // The probability that he will finish on EXACTLY on any house is the Probability Density Function (pdf).
    cout << "Probability that Pat finishes on the " << sales_quota << "th house is " << "f(5) = " << pdf(nb, nb.successes()) << endl;
    cout << "Probability that Pat finishes on the 6th house is " << pdf(nb, 6 - sales_quota) << endl;
    cout << "Probability that Pat finishes on the 7th house is " << pdf(nb, 7 - sales_quota) << endl;
    cout << "Probability that Pat finishes on the 8th house is " << pdf(nb, 8 - sales_quota) << endl;

    // The sum of the probabilities for these houses is the Cumulative Distribution Function (cdf).
    cout << "Probability that Pat finishes on or before the 8th house is sum "
      "\n" << "pdf(sales_quota) + pdf(6) + pdf(7) + pdf(8) = "
      // Sum each of the mass/density probabilities for houses sales_quota = 5, 6, 7, & 8.
      << pdf(nb, sales_quota - sales_quota) // 0
      + pdf(nb, 6 - sales_quota) // 1
      + pdf(nb, 7 - sales_quota) // 2
      + pdf(nb, 8 - sales_quota) // 3
      << endl;

    // Or using the negative binomial **cumulative** distribution function (cdf instead sum of the pdfs):
    cout << "\n""Probability of selling his quota of " << sales_quota
      << " candy bars""\n""on or before the " << 8 << "th house is "
      << cdf(nb, 8 - sales_quota) << endl;

    cout << "\n""Probability that Pat finishes exactly on the 10th house is " << pdf(nb, 10 - sales_quota) << endl;
    cout << "\n""Probability of selling his quota of " << sales_quota
      << " candy bars""\n""on or before the " << 10 << "th house is "
      << cdf(nb, 10 - sales_quota) << endl;

    cout << "Probability that Pat finishes on the 11th house is " << pdf(nb, 11 - sales_quota) << endl;
    cout << "\n""Probability of selling his quota of " << sales_quota
      << " candy bars""\n""on or before the " << 11 << "th house is "
      << cdf(nb, 11 - sales_quota) << endl;

    cout << "Probability that Pat finishes on the 12th house is " << pdf(nb, 12 - sales_quota) << endl;
    cout << "\n""Probability of selling his quota of " << sales_quota
      << " candy bars""\n""on or before the " << 12 << "th house is "
      << cdf(nb, 12 - sales_quota) << endl;

    // Finally consider the risk of Pat not setting his quota of 5 bars even after visiting all the houses.
    // Calculate the probability that he would sell on the (non-existent) last-plus-1 house.
    cout << "\n""Probability of selling his quota of " << sales_quota
      << " candy bars""\n""on or before the " << all_houses + 1 << "th house is "
      << cdf(nb, all_houses + 1 - sales_quota) << endl;
    // So the risk of failing even at the 31th (non-existent) house is 1 - this probability.
    cout << "\n""Probability of failing to sell his quota of " << sales_quota
      << " candy bars""\n""even after visiting all " << all_houses << "  houses is "
      << 1 - cdf(nb, all_houses - sales_quota + 1) << endl;

    double p = cdf(nb, (8 - sales_quota)); 
    cout << "Probability of meeting sales quota on or before 8th house is "<< p << endl;
    // Probability of meeting sales quota on or before 8th house is 0.174
    cout << "If the confidence of meeting sales quota is " << p
        << ", then the finishing house is " << quantile(nb, p) + sales_quota << endl;
    // Also try wanting absolute certainty that all 5 will be sold
    // which implies an infinite number of sales.
    // (Of course, there are only 30 houses on the estate, so he can't even be certain of selling his quota.
    cout << "If the confidence of meeting sales quota is " << 1.
        << ", then the finishing house is " << quantile(nb, 1) + sales_quota << endl; //  1.#INF == infinity.

    cout << "If the confidence of meeting sales quota is " << 0.
        << ", then the finishing house is " << quantile(nb, 0.) + sales_quota << endl;

    cout << "If the confidence of meeting sales quota is " << 1 - 0.00151 // 30 th
        << ", then the finishing house is " << quantile(nb, 1 - 0.00151) + sales_quota << endl;

    // If the opposite is true, we don't want to assume any confidence, then
    // this is tantamount to assuming that the first sales_quota will be successful.
    cout << "If confidence of meeting quota is zero (we assume all houses are successful sales)" 
      ", then finishing house is " << sales_quota << endl;


    int const pssize  = 11;
    double ps[pssize] = {0., 0.001, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 0.999, 1.};
    for (int i = 0; i < pssize; i++)
    {
      cout << "If confidence of meeting quota is " << ps[i]
        << ", then finishing house is " << ceil(quantile(nb, ps[i])) + sales_quota << endl;
    }

    cout << "If we demand a confidence of meeting sales quota of unity"
      ", then we can never be certain of selling 5 bars, so the finishing house is infinite!"  << endl;
  }
  catch(const std::exception& e)
   {
      std::cout <<
          "\n""Message from thrown exception was:\n   " << e.what() << std::endl;
   }

	return 0;
}  // int main()

/*

Output is:

Example 1 using the Negative Binomial Distribution.  ..\..\..\..\..\..\Boost-sandbox\math_toolkit\libs\math\example\negative_binomial_example1.cpp Wed Aug  1 14:12:58 2007 140050727
Selling candy bars - an example of using the negative binomial distribution. 
An example by Dr. Diane Evans,
Professor of Mathematics at Rose-Hulman Institute of Technology,
see http://en.wikipedia.org/wiki/Negative_binomial_distribution
Pat has a sales per house success rate of 0.4.
Therefore he would, on average, sell 40 bars after trying 100 houses.
With a success rate of 0.4, he might expect, on average,
 to need to visit about 12 houses in order to sell all 5 candy bars. 
Probability that Pat finishes on the 5th house is f(5) = 0.10033
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
on or before the 31th house is 0.99897
Probability of failing to sell his quota of 5 candy bars
even after visiting all 30  houses is 0.0010314
Probability of meeting sales quota on or before 8th house is 0.17367
If the confidence of meeting sales quota is 0.17367, then the finishing house is 8
Message from thrown exception was:
   Error in function boost::math::quantile(const negative_binomial_distribution<double>&, double): Probability argument is 1, which implies infinite failures !

*/






