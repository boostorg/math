// distribution_construction.cpp

// Copyright Paul A. Bristow 2007.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// Examples of using constructing distributions, mainly negative_binomial.
// The structure of distributions is rather different from some other
// statistical libraries, in less object oriented language like FOTRTRAN and C,
// that provide a few arguments to each free function.
// This library provides each distribution as a template C++ class.
// A distribution is constructed with a few arguments, and then
// member and non-member functions are used to find values of the
// distribution, often a functions of a random variate.

#include <boost/math/distributions/negative_binomial.hpp> // for negative_binomial_distribution
  using boost::math::negative_binomial_distribution; // default type is double.
	using boost::math::negative_binomial; // typedef provides default type is double.
#include <boost/math/distributions/binomial.hpp> // for binomial_distribution.
#include <boost/math/distributions/beta.hpp> // for beta_distribution.
#include <boost/math/distributions/gamma.hpp> // for gamma_distribution.

#include <iostream>
	using std::cout;
	using std::endl;

int main()
{
	cout << "Examples of constructing Distributions (negative_binomial).";
  #if defined(__FILE__) && defined(__TIMESTAMP__)
  	cout << "  " << __FILE__ << ' ' << __TIMESTAMP__ << ' '<< _MSC_FULL_VER << "\n";
  #endif
	cout << endl;

  // Several examples of constructing distributions, for example, negative binomial:

  // Distributions are fundamentally constructed like this:
  boost::math::negative_binomial_distribution<double> mydist0(8., 0.25);
  // But this is inconveniently long, so by writing
  using namespace boost::math;
  // or 
  using boost::math::negative_binomial_distribution;
  // can reduce typing.

  // Since the vast majority of applications use double,
  // the RealType default is made double, so one can write:

  negative_binomial_distribution<>  mydist9(8., 0.25); // Uses default RealType = double.

  // But the name "negative_binomial_distribution" is still inconveniently long,
  // so for most distributions, a typedef is provided, for example:
  // typedef negative_binomial_distribution<double> negative_binomial; // Reserved name of type double.
  // (Unless a  clash would occur with the name of a function, for example gamma, beta)

  using boost::math::negative_binomial;
  // Allows convenient reference to negative_binomial_distribution<double> thus:
  negative_binomial mydist(8., 0.25); //

  // You can, of course, still provide the type explicitly thus:
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

  // Deafult arguments to distribution constructors.
  // Note that default constructor arguments are only provided for some distributions.
  //   negative_binomial_distribution<> mydist8;
  //   error C2512 no appropriate default constructor available.
  // No default constructor are provided,
  // because it is difficult to chose any sensible default values for this distribution.
  // For other distribution, like the normal distribution,
  // it is obviously very useful to provide 'standard'
  // defaults for the mean and standard deviation thus:
  // normal_distribution(RealType mean = 0, RealType sd = 1)

  // Some examples using the provided typedef:
  using boost::math::negative_binomial; // Convenient access to the name.
  // Allows convenient reference to negative_binomial of default type double.
  negative_binomial mydist10(5., 0.4); // Both arguments double.
  // And automatic conversion takes place, so you can use integers and floats:
  negative_binomial mydist11(5, 0.4); // Using provided typedef double, int and double arguments.
  // This is probably the most common usage.
  negative_binomial mydist12(5., 0.4F); // Double and float arguments.
  negative_binomial mydist13(5, 1); // Both arguments integer.

  // For cases when the typdef distribution name
  // would clash with a math special function
  // (for example: binomial, beta and gamma)
  // the typedef is deliberately not provided, and the longer version(s) must be used.
  // For example:
  using namespace boost::math;
  // NOT binomial mybd0(1, 0.5); // error C2065: 'binomial' : undeclared identifier

  using boost::math::beta;
  // NOT beta mybetad0(1, 0.5); because beta is a math FUNCTION!
  // error C2146: syntax error : missing ';' before identifier 'mybetad0'
  // warning C4551: function call missing argument list
  // error C3861: 'mybetad0': identifier not found
  // Instead use:
  using boost::math::beta_distribution;
  beta_distribution<> mybetad1(1, 0.5);
  // Other examples are:
  binomial_distribution<> mybinomd1(1, 0.5);
  gamma_distribution<> mygammad1(1, 0.5);

  return 0;
}  // int main()

/*

Output is:

Examples of constructing Distributions (Negative_binomial).  ..\..\..\..\..\..\Boost-sandbox\

*/






