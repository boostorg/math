// example_error_handling.cpp

// Copyright Paul A. Bristow 2007.
// Copyright John Maddock 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// Shows use of macro definition to change policy for
// domain_error - negative degrees of freedom argument
// for student's t distribution CDF.

// Uncomment this line to see the effect of changing policy,
// to ignore the error & return NaN instead of throwing an exception.
// #define BOOST_MATH_DOMAIN_ERROR_POLICY ignore_error
// Note that these policy #defines MUST preceed the #include of
// any boost/math #includes.
// If placed after, they will have no effect!
// warning C4005: 'BOOST_MATH_OVERFLOW_ERROR_POLICY' : macro redefinition
// is a warning that it will NOT have the desired effect.

// Boost
#include <boost/math/distributions/students_t.hpp>
	using boost::math::students_t;  // Probability of students_t(df, t).

// std
#include <iostream>
	using std::cout;
	using std::endl;

#include <stdexcept>
	using std::exception;

int main()
{  // Example of error handling of bad argument(s) to a distribution.
	cout << "Example error handling using Student's t function. " << endl;

  double degrees_of_freedom = -1; double t = -1.; // Bad arguments!
  // If we use
  // cout << "Probability of Student's t is " << cdf(students_t(-1), -1) << endl; 
  // Will terminate/abort (without try & catch blocks)
  // if BOOST_MATH_DOMAIN_ERROR_POLICY has the default value throw_on_error.

  try
  {
    cout << "Probability of Student's t is " << cdf(students_t(degrees_of_freedom), t) << endl; 
  }
  catch(const std::exception& e)
  {
    std::cout <<
      "\n""Message from thrown exception was:\n   " << e.what() << std::endl;
  }

	return 0;
} // int main()

/*

Output:

Example error handling using Student's t function.

With

#define BOOST_MATH_DOMAIN_ERROR_POLICY ignore_error

Example error handling using Student's t function. 
Probability of Student's t is 1.#QNAN

Default behaviour without:

Example error handling using Student's t function. 
Message from thrown exception was:
   Error in function boost::math::students_t_distribution<double>::students_t_distribution: Degrees of freedom argument is -1, but must be > 0 !


*/
