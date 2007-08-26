// test_find_location.cpp

// Copyright John Maddock 2006.
// Copyright  Paul A. Bristow 2007.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// Basic sanity test for find_location Function.

// Default domain error policy is
// #define BOOST_MATH_DOMAIN_ERROR_POLICY throw_on_error

#include <boost/math/distributions/normal.hpp> // for normal_distribution
  using boost::math::normal; 
  using boost::math::normal_distribution; 
#include <boost/math/distributions/cauchy.hpp> // for cauchy_distribution
  using boost::math::cauchy;
#include <boost/math/distributions/pareto.hpp> // for cauchy_distribution
  using boost::math::pareto;
#include <boost/math/distributions/find_location.hpp>
  using boost::math::find_location;

  using boost::math::policies::policy;

#include <boost/test/included/test_exec_monitor.hpp> // for test_main
#include <boost/test/floating_point_comparison.hpp> // for BOOST_CHECK_CLOSE_FRACTION, BOOST_CHECK_EQUAL...

#include <iostream>
using std::cout;
using std::endl;
using std::fixed;
using std::right;
using std::left;
using std::showpoint;
using std::showpos;
using std::setw;
using std::setprecision;

#include <limits>
using std::numeric_limits;

template <class RealType> // Any floating-point type RealType.
void test_spots(RealType)
{ // Parameter only provides the type, float, double... value ignored.

  // Basic sanity checks, test data may be to double precision only
  // so set tolerance to 100 eps expressed as a fraction,
  // or 100 eps of type double expressed as a fraction,
  // whichever is the larger.

  RealType tolerance = (std::max)
      (boost::math::tools::epsilon<RealType>(),
      static_cast<RealType>(std::numeric_limits<double>::epsilon()));
   tolerance *= 100;

  cout << "Tolerance for type " << typeid(RealType).name()  << " is "
    << setprecision(3) << tolerance  << " (or " << tolerance * 100 << "%)." << endl;

  BOOST_CHECK_THROW( //  probability outside 0 to 1.
       find_location<normal_distribution<RealType> >(static_cast<RealType>(0.), static_cast<RealType>(-1.), static_cast<RealType>(0.) ), std::domain_error);
  
} // template <class RealType>void test_spots(RealType)

int test_main(int, char* [])
{
	// Check for 'bad' arguments.
  BOOST_CHECK_THROW(find_location<normal>(0., -1., 0.), std::domain_error); // p outside 0 to 1.
  BOOST_CHECK_THROW(find_location<normal>(0., 2., 0.), std::domain_error); // p outside 0 to 1.
  BOOST_CHECK_THROW(find_location<normal>(numeric_limits<double>::infinity(), 0.5, 0.), std::domain_error); // z not finite.
  BOOST_CHECK_THROW(find_location<normal>(numeric_limits<double>::quiet_NaN(), -1., 0.), std::domain_error); // z not finite
  BOOST_CHECK_THROW(find_location<normal>(0., -1., numeric_limits<double>::quiet_NaN()), std::domain_error); // scale not finite

  //// Check for ab-use with unsuitable distribution(s).
  //BOOST_CHECK_THROW(find_location<pareto>(0., 0.5, 0.), std::domain_error); // 

  double tol5eps = std::numeric_limits<double>::epsilon() * 5; // 5 eps as a fraction.
  double tol100eps = std::numeric_limits<double>::epsilon() * 100; // 100 eps as a fraction.

  double mean = 3.; // kg
  double standard_deviation = 0.1; // kg
  normal n(mean, standard_deviation);
  double z = 2.9;
  double p = 0.05;
  double l = find_location<normal>(z, p, n.scale());
  BOOST_CHECK_CLOSE_FRACTION(l, 3.0644853626951472, tol100eps);

  //l = find_location<normal>(complement(z, 1 - p, n.scale()));
  BOOST_CHECK_CLOSE_FRACTION(l, 3.0644853626951472, tol100eps);


  ////using namespace boost::math; or 
  //using boost::math::normal;


  // Basic sanity-check spot values.

  // (Parameter value, arbitrarily zero, only communicates the floating point type).
  test_spots(0.0F); // Test float.
  test_spots(0.0); // Test double.
  //test_spots(0.0L); // Test long double.
//#if !BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x582))
//  test_spots(boost::math::concepts::real_concept(0.)); // Test real concept.
//#endif

  return 0;
} // int test_main(int, char* [])

/*

Output is:


*/


