// negative_binomial_example3.cpp

// Copyright Paul A. Bristow 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_MATH_THROW_ON_DOMAIN_ERROR

#ifdef _MSC_VER
#  pragma warning(disable: 4127) // conditional expression is constant.
#  pragma warning(disable: 4512) // assignment operator could not be generated.
#  pragma warning(disable: 4996) // 'std::char_traits<char>::copy' was declared deprecated.
#endif

// Example 3 of using constructing distributions, mainly negative_binomial.

#include <boost/math/distributions/negative_binomial.hpp> // for negative_binomial_distribution
  using boost::math::negative_binomial_distribution; // default type is double.
	using boost::math::negative_binomial; // typedef provides default type is double.
#include <boost/math/distributions/binomial.hpp> // for negative_binomial_distribution

#include <iostream>
	using std::cout;
	using std::endl;

int main()
{
	cout << "Example 3 constructing Distributions (Negative_binomial).";
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
  // negative_binomial_distribution<MyFPType>  mydist6(8, 1);
  // Integer arguments -> MyFPType.
  // (assuming integers are converted automatically to MyFPType).

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

------ Rebuild All started: Project: negative_binomial_example3.cpp, Configuration: Debug Win32 ------
Deleting intermediate and output files for project 'negative_binomial_example3.cpp', configuration 'Debug|Win32'
Compiling...
negative_binomial_example3.cpp
Linking...
Autorun "i:\boost-06-05-03-1300\libs\math\test\Math_test\debug\negative_binomial_example3.cpp.exe"
Example 3 constructing Distributions (Negative_binomial).  ..\..\..\..\..\..\boost-sandbox\libs\math_functions\example\negative_binomial_example3.cpp Wed Nov  8 18:17:49 2006 140050727
Build Time 0:03
Build log was saved at "file://i:\boost-06-05-03-1300\libs\math\test\Math_test\negative_binomial_example3.cpp\Debug\BuildLog.htm"
negative_binomial_example3.cpp - 0 error(s), 0 warning(s)
========== Rebuild All: 1 succeeded, 0 failed, 0 skipped ==========


*/





