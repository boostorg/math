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
#include <boost/math/distributions/exponential.hpp>
#include <boost/math/distributions/extreme_value.hpp>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/weibull.hpp>
#include <boost/math/distributions/lognormal.hpp>

template <class RealType>
void instantiate(RealType)
{
   using namespace boost;
   using namespace boost::math;
   using namespace boost::math::concepts;

   function_requires<DistributionConcept<normal_distribution<RealType> > >();
   function_requires<DistributionConcept<binomial_distribution<RealType> > >();
   function_requires<DistributionConcept<cauchy_distribution<RealType> > >();
   function_requires<DistributionConcept<chi_squared_distribution<RealType> > >();
   function_requires<DistributionConcept<exponential_distribution<RealType> > >();
   function_requires<DistributionConcept<extreme_value_distribution<RealType> > >();
   function_requires<DistributionConcept<fisher_f_distribution<RealType> > >();
   function_requires<DistributionConcept<students_t_distribution<RealType> > >();
   function_requires<DistributionConcept<weibull_distribution<RealType> > >();
   function_requires<DistributionConcept<lognormal_distribution<RealType> > >();
}


int main()
{
   instantiate(float(0));
   instantiate(double(0));
   instantiate((long double)(0));
}

