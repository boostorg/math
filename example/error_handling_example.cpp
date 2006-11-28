// example_error_handling.cpp

// Copyright Paul A. Bristow 2006.
// Copyright John Maddock 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// Uncomment these line(s) to enable throw on errors.
#define BOOST_MATH_THROW_ON_DOMAIN_ERROR
//#define BOOST_MATH_THROW_ON_OVERFLOW_ERROR
//#define BOOST_MATH_THROW_ON_OVERFLOW_ERROR
//#define BOOST_MATH_THROW_ON_UNDERFLOW_ERROR
//#define BOOST_MATH_THROW_ON_DENORM_ERROR

// Boost
#include <boost/math/distributions/students_t.hpp>
	using boost::math::students_t;  // Probability of students_t(df, t).

#include <boost/math/tools/error_handling.hpp> // for domain_error.
	using ::boost::math::tools::domain_error;
	using ::boost::math::tools::pole_error;
	using ::boost::math::tools::overflow_error;
	using ::boost::math::tools::underflow_error;
	using ::boost::math::tools::denorm_error;
	using ::boost::math::tools::logic_error;

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
  // cout << "Probability of Student's t is " << cdf(students_t(-1), -1) << endl; 
  // Will terminate/abort if #define BOOST_MATH_THROW_ON_DOMAIN_ERROR without try & catch.
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

Without

#define BOOST_MATH_THROW_ON_DOMAIN_ERROR 

Probability of Student's t is 1.#QNAN

With

#define BOOST_MATH_THROW_ON_DOMAIN_ERROR 

Message from thrown exception was:
   Error in function __thiscall
   boost::math::students_t_distribution<double>::students_t_distribution(double):
   Degrees of freedom argument is -1, but must be > 0 !

*/
