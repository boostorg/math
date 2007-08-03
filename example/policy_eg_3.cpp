//  Copyright John Maddock 2007.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>

//[policy_eg_3

#include <boost/math/distributions/binomial.hpp>

//
// Begin by defining a policy type, that gives the
// behaviour we want:
//
using namespace boost::math::policies;
typedef policy<
   promote_float<false>, 
   discrete_quantile<integer_nearest> 
> mypolicy;
//
// Then define a distribution that uses it:
//
typedef boost::math::binomial_distribution<float, mypolicy> mybinom;
//
//  And now use it to get the quantile:
//
int main()
{
   std::cout << "quantile is: " <<
      quantile(mybinom(200, 0.25), 0.05) << std::endl;
}

//]

