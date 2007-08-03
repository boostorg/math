// negative_binomial_example2.cpp

// Copyright Paul A. Bristow 2007.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// Example 2 of using constructing distributions, mainly negative_binomial.

#include <boost/math/distributions/negative_binomial.hpp> // for negative_binomial_distribution
  using boost::math::negative_binomial_distribution; // default type is double.
	using boost::math::negative_binomial; // typedef provides default type is double.
#include <boost/math/distributions/binomial.hpp> // for negative_binomial_distribution

#include <iostream>
	using std::cout;
	using std::endl;

int main()
{
	cout << "Example 2 constructing Distributions (Negative_binomial).";
  #if defined(__FILE__) && defined(__TIMESTAMP__)
  	cout << "  " << __FILE__ << ' ' << __TIMESTAMP__ << ' '<< _MSC_FULL_VER << "\n";
  #endif
	cout << endl;

  // Several examples of constructing distributions, for example, negative binomial:

  // Fundamentally constructed like this:
  boost::math::negative_binomial_distribution<double> mydist0(8., 0.25);
  // But this is inconveniently long.

  using boost::math::negative_binomial_distribution;
  // Allows convenient reference to negative_binomial_distribution.

  // You can provide the type explicitly thus:
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

  // Note that default constructor arguments are only provided for some distributions. 
  //   negative_binomial_distribution<> mydist;
  //   error C2512 no appropriate default constructor available.
  // Since there are no accessor functions, no default constructor are provided,
  // because it is difficult to chose any sensible default values for this distribution.
  // For other distribution, like the normal distribution,
  // it is obviously very useful to provide
  // defaults for the mean and standard deviation thus:
  // normal_distribution(RealType mean = 0, RealType sd = 1)

  negative_binomial_distribution<>  mydist9(8., 0.25); // Uses default RealType = double.
  // But the name "negative_binomial_distribution" is still inconveniently long,
  // so for most distributions, a typedef is provided, for example:

  // typedef negative_binomial_distribution<double> negative_binomial; // Reserved name of type double.

  // Some examples using the provided typedef:
  using boost::math::negative_binomial; // Convenient access to the name.
  // Allows convenient reference to negative_binomial of default type double.
  negative_binomial mydist10(5., 0.4); // Both arguments double.
  // And automatic conversion takes place, so you can use integers and floats:
  negative_binomial mydist11(5, 0.4); // Using provided typedef double, int and double arguments.
  // This is probably the most common usage.
  negative_binomial mydist12(5., 0.4F); // Double and float arguments. 
  negative_binomial mydist13(5, 1); // Both arguments integer.

  // But for cases when the typdef distribution name
  // would clash with a math special function
  // (for example binomial, beta and gamma)
  // the typedef is deliberately not provided, and
  // the longer version(s) must be used.
  // For example:
  using namespace boost::math;
  // NOT binomial myb010(1, 0.5); but
  binomial_distribution<> myb1(1, 0.5);

 	return 0;
}  // int main()

/*

Output is:

math_toolkit\libs\math\example\negative_binomial_construction_examples.cpp Wed Aug  1 13:59:34 2007 140050727


*/






