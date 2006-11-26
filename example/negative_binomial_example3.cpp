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


 	return 0;
}  // int main()

/*

Output is:



*/





