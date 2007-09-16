// Copyright Paul A. Bristow 2007
// Copyright John Maddock 2006

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// binomial_examples_quiz.cpp

// Simple example of computing probabilities and quantiles for a binomial random variable
// representing the correct guesses on a multiple-choice test.

// source http://www.stat.wvu.edu/SRS/Modules/Binomial/test.html

//[binomial_quiz_example1
/*`
A multiple choice test has four possible answers to each of 16 questions.
A student guesses the answer to each question,
so the probability of getting a correct answer on any given question is 1/4 = 0.25.
The conditions of the binomial experiment are assumed to be met:
n = 16 questions constitute the trials;
each question results in one of two possible outcomes (correct or incorrect);
the probability of being correct is 0.25 and is constant if no knowledge about the subject is assumed;
the questions are answered independently if the student's answer to a question
in no way influences his/her answer to another question.

The number of correct answers, X, is distributed as a binomial random variable
with binomial distribution parameters n = 16 and p = 0.25.
The program below displays the probabilities for each of the 17 possible outcomes,
i.e., for X = 0, 1, ..., 16, in a line chart.

First we need to be able to use the binomial distribution constructor
(and some std input/output, of course)
*/

#include <boost/math/distributions/binomial.hpp>
  using boost::math::binomial;

#include <iostream>
using std::cout;
using std::endl;
using std::ios;
using std::flush;
using std::left;
using std::right;
using std::fixed;
#include <iomanip>
using std::setw;
using std::setprecision;
//][/binomial_quiz_example1]

//[binomial_quiz_example2
int main()
{
  try
  {
  cout << "Binomial distribution example - guessing in a quiz." << endl;
  /*`
  The number of correct answers, X, is distributed as a binomial random variable
  with binomial distribution parameters n = 16 and p = 0.25.
  */
  cout.precision(5); // Might be able to calculate a best value for this?
  int questions = 16;
  int answers = 4; // possible answers to each question.
  double success_fraction = (double)answers / (double)questions; // If a random guess.
  // Caution:  = answers / questions would be zero (because they are integers)!
  int pass_score = 11;

  /*`
  Construct our Binomial distribution.
  */
  binomial quiz(questions, success_fraction);
  /*`
  and display the parameters we used.
  */
  cout << "In a quiz with " << quiz.trials()
    << " questions and with a probability of guessing right of "
    << quiz.success_fraction() * 100 << " %" 
    << " or 1 in " << static_cast<int>(1. / quiz.success_fraction()) << endl;
  /*`
  Show some probabilities of just guessing: these don't give any 
  encouragement to guessers!
  */
  cout << "Probability of getting none right is " << pdf(quiz, 0) << endl; // 0.010023
  cout << "Probability of getting at least one right is " << 1 - pdf(quiz, 1) << endl; // 0.94655
  cout << "Probability of getting none or one right is " << pdf(quiz, 0) + pdf(quiz, 1) << endl; // 0.063476
  cout << "Probability of getting exactly one right is " << pdf(quiz, 1) << endl;
  cout << "Probability of getting exactly 11 right is " << pdf(quiz, 11) << endl;
  cout << "Probability of getting > 10 right (to pass) is " << cdf(complement(quiz, 10)) << endl;

  // Using Binomial probabilities.
  cout << "The probability of getting all the answers wrong by chance is "
    << pdf(quiz, 0) << endl;
  cout << "The probability of getting all the answers right by chance is " 
    << pdf(quiz, questions) << endl;
  cout << "The probability of getting exactly " << pass_score
    << " answers right by guessing is " << pdf(quiz, pass_score) << endl << endl;

  cout << "The probability of getting less then " << pass_score
    << "(< " << pass_score << ") answers right by guessing is "
    << cdf(quiz, pass_score) << endl;

  cout << "The probability of getting at least " << pass_score 
    << "(>= " << pass_score << ") answers right by guessing is "
    << cdf(complement(quiz, pass_score-1))
    << " only 1 in " << 1/cdf(complement(quiz, pass_score-1)) << endl;

  /*`
  Tabulate probability versus number right.
  */
  cout << "\n" "Guessed right Probability" << right << endl;
  for (int successes = 0; successes <= questions; successes++)
  {
    double probability = pdf(quiz, successes);
    cout << setw(2) << successes << "             " << probability << endl;
  }
  cout << endl;


  cout << "\n" "At most (<=)""\n""Guessed right   Probability" << right << endl;
  for (int score = 0; score <= questions; score++)
  {
    cout << setw(2) << score << "                " << cdf(quiz, score) << endl;
  }
  cout << endl;

  cout << "\n" "At least (>=)""\n""Guessed right   Probability" << right << endl;
  for (int score = 0; score <= questions; score++)
  {
    cout << setw(2) << score << "                " << cdf(complement(quiz, score)) << endl;
  }
  /*`
  Calculate the probability of getting a range of guesses right,
  first by adding the exact probabilities of each of low ... high.
  */
  int low = 3;
  int high = 5;
  double sum = 0.;
  for (int i = low; i <= high; i++)
  {
    sum += pdf(quiz, i);
  }
  cout << "The probability of getting between "
    << low << " and " << high << " answers right by guessing is "
    << sum  << endl; // 0.61323

  /*`
  Or, better, we can use the difference of cdfs instead:
  */
  cout << "The probability of getting between " << low << " and " << high << " answers right by guessing is "
    <<  cdf(quiz, high) - cdf(quiz, low - 1) << endl; // 0.61323
  // And a few more combinations of high and low choices:
  low = 1; high = 6; 
  cout << "The probability of getting between " << low << " and " << high << " answers right by guessing is "
    <<  cdf(quiz, high) - cdf(quiz, low - 1) << endl; // 1 and 6 P= 0.91042
  low = 1; high = 8; 
  cout << "The probability of getting between " << low << " and " << high << " answers right by guessing is "
    <<  cdf(quiz, high) - cdf(quiz, low - 1) << endl; // 1 <= x 8 P = 0.9825
  low = 4; high = 4; 
  cout << "The probability of getting between " << low << " and " << high << " answers right by guessing is "
    <<  cdf(quiz, high) - cdf(quiz, low - 1) << endl; // 4 <= x 4 P = 0.22520
  low = 3; high = 5; 
  cout << "The probability of getting between " << low << " and " << high << " answers right by guessing is "
    <<  cdf(quiz, high) - cdf(quiz, low - 1) << endl; // P 3 to 5 right

  /*`
  Using Binomial distribution moments,
  we can say more about the spread of results from guessing.
  */
  cout << "By guessing, on average, one can expect to get " << mean(quiz) << " correct answers." << endl;
  cout << "Standard deviation is " << standard_deviation(quiz) << endl;
  cout << "So about 2/3 will lie within 1 standard deviation and get between "
    <<  ceil(mean(quiz) - standard_deviation(quiz))  << " and "
    << floor(mean(quiz) + standard_deviation(quiz)) << " correct." << endl; 
  cout << "Mode (the most frequent) is " << mode(quiz) << endl;
  cout << "Skewness is " << skewness(quiz) << endl;

  /*`
  Show the use of quantiles (percentiles or percentage points) for a 
  few probability levels:
  */
  cout << "Quantiles" << endl;
  cout << "Quartiles " << quantile(quiz, 0.25) << " to "
    << quantile(complement(quiz, 0.25)) << endl; // Quartiles 2.2821 4.6212
  cout << "1 sd " << quantile(quiz, 0.33) << " to " 
    << quantile(quiz, 0.67) << endl; // 1 sd 2.6654 4.1935
  cout << "Deciles " << quantile(quiz, 0.1)  << " to "
    << quantile(complement(quiz, 0.1))<< endl; // Deciles 1.3487 5.7583
  cout << "5 to 95% " << quantile(quiz, 0.05)  << " to "
    << quantile(complement(quiz, 0.05))<< endl; // 5 to 95% 0.83739 6.4559
  cout << "2.5 to 97.5% " << quantile(quiz, 0.025) << " to "
    <<  quantile(complement(quiz, 0.025)) << endl; // 2.5 to 97.5% 0.42806 7.0688
  cout << "2 to 98% " << quantile(quiz, 0.02)  << " to "
    << quantile(complement(quiz, 0.02)) << endl; //  2 to 98% 0.31311 7.7880

  cout << "If guessing then percentiles 1 to 99% will get " << quantile(quiz, 0.01) 
    << " to " << quantile(complement(quiz, 0.01)) << " right." << endl; 

//] [/binomial_quiz_example2]

//[discrete_quantile_real
/*`
The quantiles values are controlled by the discrete quantile policy chosen.
The default is `integer_round_outwards`,
so the lower quantile is rounded down, and the upper quantile is rounded up.

We can control the policy for all distributions by
  #define BOOST_MATH_DISCRETE_QUANTILE_POLICY real
at the head of the program would make this policy apply
to this *one, and only*, translation unit.

Or we can create a (typedef for) policy that has discrete quantiles real.
*/
  using namespace boost::math::policies;
/*`
Convenient for all policy and typelist values like discrete_quantile.
*/
  using namespace boost::math;
/*`
for binomial_distribution

Or to be more specific, to avoid 'using namespaces ...' statements:
*/
  using boost::math::policies::policy;
  using boost::math::policies::discrete_quantile;
  using boost::math::policies::real;
  using boost::math::policies::integer_round_outwards; // Default.
  typedef boost::math::policies::policy<discrete_quantile<real> > real_quantile_policy;
  /*`
  Add a binomial distribution called real_quantile_binomial that uses real_quantile_policy.
  */
  using boost::math::binomial_distribution;
  typedef binomial_distribution<double, real_quantile_policy> real_quantile_binomial;
  /*`
  Construct a distribution of this custom real_quantile_binomial distribution;
  */
  real_quantile_binomial quiz_real(questions, success_fraction);
  /*`
  And use this to show some quantiles - that now have real rather than integer values.
  */

  cout << "Real Quartiles " << quantile(quiz_real, 0.25) << " to "
    << quantile(complement(quiz_real, 0.25)) << endl; // Real Quartiles 2.2821 to 4.6212

//] [/discrete_quantile_real]
  }
  catch(const std::exception& e)
  { // Always useful to include try & catch blocks because
    // default policies are to throw exceptions on arguments that cause
    // errors like underflow, overflow. 
    // Lacking try & catch blocks, the program will abort without a message below,
    // which may give some helpful clues as to the cause of the exception.
    std::cout <<
      "\n""Message from thrown exception was:\n   " << e.what() << std::endl;
  }
  return 0;
} // int main()



/*

Output is:

Binomial distribution example - guessing in a quiz.
In a quiz with 16 questions and with a probability of guessing right of 25 % or 1 in 4
Probability of getting none right is 0.010023
Probability of getting at least one right is 0.94655
Probability of getting none or one right is 0.063476
Probability of getting exactly one right is 0.053454
Probability of getting exactly 11 right is 0.00024713
Probability of getting > 10 right (to pass) is 0.00028524
The probability of getting all the answers wrong by chance is 0.010023
The probability of getting all the answers right by chance is 2.3283e-010
The probability of getting exactly 11 answers right by guessing is 0.00024713
The probability of getting less then 11(< 11) answers right by guessing is 0.99996
The probability of getting at least 11(>= 11) answers right by guessing is 0.00028524 only 1 in 3505.8
Guessed right Probability
 0             0.010023
 1             0.053454
 2             0.13363
 3             0.20788
 4             0.2252
 5             0.18016
 6             0.1101
 7             0.052427
 8             0.01966
 9             0.0058253
10             0.0013592
11             0.00024713
12             3.4324e-005
13             3.5204e-006
14             2.5146e-007
15             1.1176e-008
16             2.3283e-010
At most (<=)
Guessed right   Probability
 0                0.010023
 1                0.063476
 2                0.19711
 3                0.40499
 4                0.63019
 5                0.81035
 6                0.92044
 7                0.97287
 8                0.99253
 9                0.99836
10                0.99971
11                0.99996
12                1
13                1
14                1
15                1
16                1
At least (>=)
Guessed right   Probability
 0                0.98998
 1                0.93652
 2                0.80289
 3                0.59501
 4                0.36981
 5                0.18965
 6                0.079557
 7                0.02713
 8                0.0074697
 9                0.0016445
10                0.00028524
11                3.8107e-005
12                3.7833e-006
13                2.6287e-007
14                1.1409e-008
15                2.3283e-010
16                0
The probability of getting between 3 and 5 answers right by guessing is 0.61323
The probability of getting between 3 and 5 answers right by guessing is 0.61323
The probability of getting between 1 and 6 answers right by guessing is 0.91042
The probability of getting between 1 and 8 answers right by guessing is 0.98251
The probability of getting between 4 and 4 answers right by guessing is 0.2252
The probability of getting between 3 and 5 answers right by guessing is 0.61323
By guessing, on average, one can expect to get 4 correct answers.
Standard deviation is 1.7321
So about 2/3 will lie within 1 standard deviation and get between 3 and 5 correct.
Mode (the most frequent) is 4
Skewness is 0.28868
Quantiles
Quartiles 2 to 5
1 sd 2 to 5
Deciles 1 to 6
5 to 95% 0 to 7
2.5 to 97.5% 0 to 8
2 to 98% 0 to 8
If guessing then percentiles 1 to 99% will get 0 to 8 right.
Real Quartiles 2.2821 to 4.6212

*/

