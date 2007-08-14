// distribution_construction.cpp

// Copyright Paul A. Bristow 2007.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// Caution: this file contains Quickbook markup as well as code
// and comments, don't change any of the special comment markups!

//[distribution_construction1

/*`

The structure of distributions is rather different from some other statistical libraries,
for example in less object-oriented language like FORTRAN and C,
that provide a few arguments to each free function.
This library provides each distribution as a template C++ class.
A distribution is constructed with a few arguments, and then
member and non-member functions are used to find values of the
distribution, often a function of a random variate.

First we need some includes to access the negative binomial distribution
(and the binomial, beta and gamma too).

*/

#include <boost/math/distributions/negative_binomial.hpp> // for negative_binomial_distribution
  using boost::math::negative_binomial_distribution; // default type is double.
  using boost::math::negative_binomial; // typedef provides default type is double.
#include <boost/math/distributions/binomial.hpp> // for binomial_distribution.
#include <boost/math/distributions/beta.hpp> // for beta_distribution.
#include <boost/math/distributions/gamma.hpp> // for gamma_distribution.

//] [/distribution_construction1 end of Quickbook in C++ markup]

//[distribution_construction2

/*`
Several examples of constructing distributions follow:
*/

int main()
{
  // First, a negative binomial distribution with 8 successes
  // and a success fraction 0.25, 25% or 1 in 4, is constructed like this:

  boost::math::negative_binomial_distribution<double> mydist0(8., 0.25);
  // But this is inconveniently long, so by writing
  using namespace boost::math;
  // or
  using boost::math::negative_binomial_distribution;
  // we can reduce typing.

  // Since the vast majority of applications use double,
  // the RealType default is chosen to be double, so we can also write:

  negative_binomial_distribution<> mydist9(8., 0.25); // Uses default RealType = double.

  // But the name "negative_binomial_distribution" is still inconveniently long,
  // so for most distributions, a typedef is provided, for example:
  // typedef negative_binomial_distribution<double> negative_binomial; // Reserved name of type double.
  // (Unless a  clash would occur with the name of a function, for example gamma & beta)
  // So, after a using statement,
  using boost::math::negative_binomial;
  // we have a convenient reference to negative_binomial_distribution<double> thus:
  negative_binomial mydist(8., 0.25); //

  // Some more examples using the provided convenience typedef:
  negative_binomial mydist10(5., 0.4); // Both arguments double.
  // And automatic conversion takes place, so you can use integers and floats:
  negative_binomial mydist11(5, 0.4); // Using provided typedef double, int and double arguments.
  // This is probably the most common usage.
  negative_binomial mydist12(5., 0.4F); // Double and float arguments.
  negative_binomial mydist13(5, 1); // Both arguments integer.

  // Similarly for most other distributions like the binomial.
  binomial mybinomial(1, 0.5); // is more concise than
  binomial_distribution<> mybinomd1(1, 0.5);

  // For cases when the typdef distribution name would clash with a math special function
  // (for example: beta and gamma)
  // the typedef is deliberately not provided, and the longer version(s) must be used.
  // For example:
  using boost::math::beta;
  // NOT beta mybetad0(1, 0.5); because beta is a math FUNCTION!
  // error C2146: syntax error : missing ';' before identifier 'mybetad0'
  // warning C4551: function call missing argument list
  // error C3861: 'mybetad0': identifier not found
  // Instead use:
  using boost::math::beta_distribution;
  beta_distribution<> mybetad1(1, 0.5);
  // or for the gamm distribution:
  gamma_distribution<> mygammad1(1, 0.5);

  // We can, of course, still provide the type explicitly thus:
  negative_binomial_distribution<double> mydist1(8., 0.25); // Explicit double.
  negative_binomial_distribution<float>  mydist2(8., 0.25); // Explicit float, double arguments -> float.
  negative_binomial_distribution<float>  mydist3(8, 0.25); // Explicit float, integer & double arguments -> float.
  negative_binomial_distribution<float>  mydist4(8.F, 0.25F); // Explicit float, float arguments, no conversion.
  negative_binomial_distribution<float>  mydist5(8, 1); // Explicit integer, integer arguments -> float.
  negative_binomial_distribution<double> mydist6(8., 0.25); // Explicit double.
  negative_binomial_distribution<long double> mydist7(8., 0.25); // Explicit long double.
  // And if you have your own RealType called MyFPType,
  // for example NTL quad_float (128-bit floating-point), then:
  // negative_binomial_distribution<MyFPType>  mydist6(8, 1); // Integer arguments -> MyFPType.

  /*`
  Default arguments to distribution constructors.
  */
  // Note that default constructor arguments are only provided for some distributions.
  // So if you wrongly assume a default argument you will get an error message, for example:
  //   negative_binomial_distribution<> mydist8;
  //   error C2512 no appropriate default constructor available.

  // No default constructors are provided for the negative binomial,
  // because it is difficult to chose any sensible default values for this distribution.
  // For other distributions, like the normal distribution,
  // it is obviously very useful to provide 'standard'
  // defaults for the mean and standard deviation thus:
  // normal_distribution(RealType mean = 0, RealType sd = 1)

  return 0;
}  // int main()

/*`There is no useful output from this program, of course. */

//] [/end of distribution_construction2]

