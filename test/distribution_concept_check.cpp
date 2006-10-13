//  Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/concepts/distributions.hpp>

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/cauchy.hpp>
#include <boost/math/distributions/chi_squared.hpp>


int main()
{
   using namespace boost;
   using namespace boost::math;
   using namespace boost::math::concepts;

   function_requires<DistributionConcept<normal> >();
   function_requires<DistributionConcept<students_t> >();
   function_requires<DistributionConcept<binomial> >();
   function_requires<DistributionConcept<cauchy> >();
   function_requires<DistributionConcept<chi_squared> >();
}

