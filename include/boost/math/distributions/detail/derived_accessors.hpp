//  Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_STATS_DERIVED_HPP
#define BOOST_STATS_DERIVED_HPP

//
// This file implements various common properties of distributions
// that can be implemented in terms of other properties.
//
// Of course there may be more efficient means of implementing these 
// that are specific to a particular distribution, but these generic
// versions give these properties "for free" with most distributions.
//
// In order to make use of this header, it must be included AT THE END
// of the distribution header, AFTER the distribution and it's core
// property accessors have been defined: this is so that compilers
// that implement 2-phase lookup and early-type-checking of templates
// can find the definitions refered to herein.
//
// Note that while both variance and standard_deiviation are provided
// here, each distribution must specialise at least one of these
// otherwise these two versions will just call each other over and over
// until stack space runs out....
//
namespace boost{ namespace math{

template <class Distribution>
typename Distribution::value_type variance(const Distribution& dist);

template <class Distribution>
inline typename Distribution::value_type standard_deviation(const Distribution& dist)
{
   using namespace std;  // ADL of sqrt.
   return sqrt(variance(dist));
}

template <class Distribution>
inline typename Distribution::value_type variance(const Distribution& dist)
{
   typename Distribution::value_type result = standard_deviation(dist);
   return result * result;
}

template <class Distribution, class RealType>
typename Distribution::value_type hazard(const Distribution& dist, const RealType& x)
{
   typedef typename Distribution::value_type value_type;
   value_type p = cdf(complement(dist, x));
   value_type d = pdf(dist, x);
   if(d > p * tools::max_value<value_type>())
      return tools::overflow_error<value_type>(
         BOOST_CURRENT_FUNCTION, 0);
   if(d == 0)
   {
      // This protects against 0/0, but is it the right thing to do?
      return 0;  
   }
   return d / p;
}

template <class Distribution, class RealType>
inline typename Distribution::value_type chf(const Distribution& dist, const RealType& x)
{
   using namespace std;
   return -log(cdf(complement(dist, x)));
}

template <class Distribution>
typename Distribution::value_type coefficient_of_variation(const Distribution& dist)
{
   typedef typename Distribution::value_type value_type;
   value_type m = mean(dist);
   value_type d = standard_deviation(dist);
   if(d > m * tools::max_value<value_type>())
      return tools::overflow_error<value_type>(
         BOOST_CURRENT_FUNCTION, 0);
   return d / m;
}


} // namespace math
} // namespace boost

#endif // BOOST_STATS_DERIVED_HPP
