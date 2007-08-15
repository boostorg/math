// negative_binomial_example3.cpp

// Copyright Paul A. Bristow 2007.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// Example 3 of using constructing distributions, mainly negative_binomial.

//[neg_binomial_example3
/*`
First we need some includes to access the negative binomial distribution
(and some basic std output of course).
*/

#include <boost/math/distributions/negative_binomial.hpp> // for negative_binomial_distribution
	using boost::math::negative_binomial; // typedef provides default type is double.

#include <iostream>
	using std::cout;
	using std::endl;
  using std::setw;
#include <limits>
  using std::numeric_limits;

/*`
A number of the application examples from K. Krishnamoorthy, Handbook of Statistical Distributions with Applications,
ISBN 1 58488 635 8, page 100... are implemented using the Math Toolkit library.

A sample example is shown here:
*/

//] [/ neg_binomial_example3]

int main()
{
  cout << "Example 3 Negative_binomial, Krishnamoorthy applications.";
#if defined(__FILE__) && defined(__TIMESTAMP__)
  cout << "  " << __FILE__ << ' ' << __TIMESTAMP__ << "\n";
#endif
  cout << endl;
  { 
//[neg_binomial_example3_1
    // K. Krishnamoorthy, Handbook of Statistical Distributions with Applications,
    // ISBN 1 58488 635 8, page 100, example 7.3.1
    // r successes = 20, k failures = 18, and success probability (fraction) = 0.6
    negative_binomial nbk(20, 0.6); // successes, success_fraction.
    cout.precision(6); // default precision.
    cout << "Probability of <= 18, P(X <= 18) == cdf(18) = " <<  cdf(nbk, 18) << endl; // = 0.862419
    cout << "Probability of > 18, P(X < 18) ==  1 - cdf(18) = " <<  1 - cdf(nbk, 18) << endl; // = 0.137581
    // But this may suffer from inaccuracy, so much better to use the complement.
    cout << "Probability of > 18, P(X < 18) ==  1 - cdf(18) = " <<  cdf(complement(nbk, 18)) << endl; // 0.137581
    // And, of course, the sum of probability of X <= 18 and X > 18 really should be unity!
    BOOST_ASSERT(cdf(nbk, 18) + cdf(complement(nbk, 18))); //  = 0.862419 + 0.137581 == 1
    cout << "Probability of < 18 == <= 17 == cdf(17) = " << cdf(nbk, 17) << endl; //
    cout << "Probability of >= 18 == > 17 == 1 - cdf(17) = " <<  1 - cdf(nbk, 17) << endl; // 0.181983
    // But this may suffer from inaccuracy, so much better to use the complement.
    cout << "Probability of >= 18 == > 17 == cdf(complement(nbk, 17))  = "
      << cdf(complement(nbk, 17)) << endl; //  = 0.181983
    cout << "Probability of exactly == 18, P(X =18) = pdf(18) = " <<  pdf(nbk, 18) << endl; // 0.044402
    cout << endl;
//] [/negative_binomial_example3_1 Quickbook end]
  }
  { // Example 7.3.2 to find the success probability when r = 4, k = 5, P(X <= k) = 0.56 and
    // calculates success probability = 0.417137.
    negative_binomial nb(4, 0.417137);
    cout <<  "P(X <= k) = " << cdf(nb, 5) << endl; // expect P(X <= k) = 0.56 got 0.560001
    // Now check the inverse by calculating the k failures.
    cout << " k = " << quantile(nb, 0.56) << endl; // expect to get k = 5. and got 5.000000
    // P(X <= k) = 0.56001
    // P(X >= k) = 0.55406
    // P(X == k) = 0.114061
    // Compute moments:
    // mean 5.58918, sd 3.66045, mode = 4, cv = 0.654918, skew 1.03665, Mean dev 2.86877, kurtosis 4.57463
    cout << "Mean = " << mean(nb) //  5.589176
      << ", sd = " << standard_deviation(nb) 
      << ", Coefficient of variation (sd/mean) = " << standard_deviation(nb) / mean(nb) 
      << ", mode = " << mode(nb) << ", skew " << skewness(nb) 
     // << ", mean deviation = " << "??" // perhaps todo?
      << ", kurt = " << kurtosis(nb) << endl; 
  }
  { // 7.3.3 Coin flipping.  What are chances that 10th head will occur at the 12th flip.
    negative_binomial nb(10, 0.5);

    cout << "Probability that 10th head will occur at the 12th flip is " << pdf(nb, 12-10) << endl; // 0.013428
    cout << "P(X == k) " << pdf(nb, 12-10) << endl; // 0.013428

    cout << "Probability that 10th head will occur before the 12th flip is " << cdf(nb, 1) << endl; //
    cout << "P(X < k) is " << cdf(nb, 1) << endl; // 1 = 12 - 11, // 0.005859

    cout << "Probability that 10th head will occur on or before the 12th flip is " << cdf(nb, 2) << endl; //
    cout << "P(X <= k) is " << cdf(nb, 2) << endl; //  P(X <= k) = 0.01928711

    cout << "Probability that 10th head will occur on or after the 12th flip is " << 1 - cdf(nb, 1) << endl; //
    cout << "P(X >= k) is " << cdf(complement(nb, 12-11)) << endl; // P(X >= k) =  0.994140625

    cout << "Probability that 10th head will occur after the 12th flip is " << 1 - cdf(nb, 2) << endl; //
    cout << "P(X > k) is " << cdf(complement(nb, 2)) << endl; //  P(X > k) 0.980713
    /*    
    Probability that 10th head will occur at the 12th flip is 0.013428
    P(X == k) 0.013428
    Probability that 10th head will occur before the 12th flip is 0.005859
    P(X < k) is 0.005859
    Probability that 10th head will occur on or before the 12th flip is 0.019287
    P(X <= k) is 0.019287
    Probability that 10th head will occur on or after the 12th flip is 0.994141
    P(X >= k) is 0.994141
    Probability that 10th head will occur after the 12th flip is 0.980713
    P(X > k) is 0.980713
    */
  }
  { // Example 7.3.4, Sample up to 30 items from each lot.
    // Chances of rejecting the lot if it indeed contains 15% defectives?
    // stop at the 3rd defective, or go on to sample all 30.
    negative_binomial nb(3, 0.15); // 3rd defective.
    // probability of observing 27 or less OK items to get the 3rd defective.
    // P(X <= k) = 0.8485994
    cout << "Chance of rejecting the lot is " << cdf(nb, 27) << " if the lot actually contains 15% defectives." << endl; // 0.848599
  }
  { // K. Krishnamoorthy, ISBN 1 58488 635 8, page 101, section 7.46 confidence Intervals for the Proportion
    // 7.6.1 Shipment inspected on-by-one randomly and found 6th defective at 30th inspection.
    // 2-sided 95% confidence 0.0771355 to 0.357748
    double alpha = 0.05; // Note divide by 2 if two-sided test.
    double successes = 6;
    double failures = 24;
    double trials = successes + failures;
    double l = negative_binomial::find_lower_bound_on_p(trials, successes, alpha/2); // 0.077136
    double u = negative_binomial::find_upper_bound_on_p(trials, successes, alpha/2); // 0.357748
    // These find_ use the Clopper-Pearson approach.
    cout << (1 - alpha) * 100 << "% confidence interval low " << l << ", up " << u << endl;
    cout << "So we conclude that the true percentage of defective items is between "
      << static_cast<int>(l * 100) << " and " << setw(2) << u * 100 << "%." << endl;

    cout << "The rounding style for a double type is: " 
      << numeric_limits<double>::round_style << endl;
  }
  { // K. Krishnamoorthy, ISBN 1 58488 635 8, page 102, section 7.4 & 7.5
    // Suppose required k + r trials to get the rth success.
    double r = 5;     double k = 25; // Example 7.5.1 values
    double rm1 = r -1; // r-1 is the penultimate success.
    double p = rm1/ (rm1 + k); 
    // A point estimate phat of the actual proportion of defective items.
    //  0.137931
    // 'True' estimate, if this was all the items ever available, is successes/failures = r/k = 5/25 = 0.2
    // so the point estimate 'how we are doing from the info so far' is rather less at 0.138.
    cout << "Uniformly minimum variance unbiased estimator of success probability is " << p << endl; 
    double v = 0.5 * p * p * (1 - p) * (k + k + 2 - p) / (k * (k - p + 2));
    cout << "Variance of estimator is " << v << endl;  // 0.000633
    cout << "Standard deviation of estimator is " << sqrt(v) << endl;  // 0.025 - seems smallish?
    // Would expect point estimate 0.14 + a couple of sd = 0.25 * 2 = 0.05 = 1.9 ~= 0.2?
    // So perhaps near enough?

    // Section 7.5 testing the true success probability.
    // 7.5.1  A lot is inspected randomly.  5th defective found at 30th inspection.
    // so 25 are OK and 5 are duds.
    // Using StatCalc r = 5, k = 25, so r + r = 30, value for p0 0.3
    negative_binomial nb(5, 0.2); // 1/3rd are defective?
    cout << "P(X <= 25) = " << cdf(nb, 25) << endl; //  0.744767
    cout << "quantile(nb, 0.3) = " << quantile(nb, 0.3) << endl; // 13
    // Don't understand what the p-value calculation does.
    // TODO????
  }
  return 0;
}  // int main()


/*

Output is:

Example 3 Negative_binomial, Krishnamoorthy applications.  ..\..\..\..\..\..\boost-sandbox\math_toolkit\libs\math\example\negative_binomial_example3.cpp Wed Aug 15 11:08:14 2007
Probability of <= 18, P(X <= 18) == cdf(18) = 0.862419
Probability of > 18, P(X < 18) ==  1 - cdf(18) = 0.137581
Probability of > 18, P(X < 18) ==  1 - cdf(18) = 0.137581
Probability of < 18 == <= 17 == cdf(17) = 0.818017
Probability of >= 18 == > 17 == 1 - cdf(17) = 0.181983
Probability of >= 18 == > 17 == cdf(complement(nbk, 17))  = 0.181983
Probability of exactly == 18, P(X =18) = pdf(18) = 0.0444024
P(X <= k) = 0.560001
 k = 5
Mean = 5.58918, sd = 3.66045, Coefficient of variation (sd/mean) = 0.654918, mode = 4, skew 1.03665, kurt = 4.57463
Probability that 10th head will occur at the 12th flip is 0.0134277
P(X == k) 0.0134277
Probability that 10th head will occur before the 12th flip is 0.00585938
P(X < k) is 0.00585938
Probability that 10th head will occur on or before the 12th flip is 0.0192871
P(X <= k) is 0.0192871
Probability that 10th head will occur on or after the 12th flip is 0.994141
P(X >= k) is 0.994141
Probability that 10th head will occur after the 12th flip is 0.980713
P(X > k) is 0.980713
Chance of rejecting the lot is 0.848599 if the lot actually contains 15% defectives.
95% confidence interval low 0.0771355, up 0.357748
So we conclude that the true percentage of defective items is between 7 and 35.7748%.
The rounding style for a double type is: 1
Uniformly minimum variance unbiased estimator of success probability is 0.137931
Variance of estimator is 0.000633295
Standard deviation of estimator is 0.0251654
P(X <= 25) = 0.744767
13


*/



