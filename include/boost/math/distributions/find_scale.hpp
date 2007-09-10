//  Copyright John Maddock 2007.
//  Copyright Paul A. Bristow 2007.

//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_STATS_FIND_SCALE_HPP
#define BOOST_STATS_FIND_SCALE_HPP

#include <boost/math/distributions/fwd.hpp> // for all distribution signatures.
#include <boost/math/distributions/complement.hpp>
#include <boost/math/policies/policy.hpp>
#include <boost/math/tools/traits.hpp>
#include <boost/static_assert.hpp>
// using boost::math::policies::policy;
// using boost::math::complement; // will be needed by users who want complement,
// but NOT placed here to avoid putting it in global scope.

namespace boost
{
  namespace math
  {
  // Function to find location of random variable z
  // to give probability p (given scale)
  // Apply to normal, lognormal, extreme value, Cauchy, (and symmetrical triangular).
  // BOOST_STATIC_ASSERTs are used to enforce this.

    template <class Dist, class Policy>
    inline
      typename Dist::value_type find_scale( // For example, normal mean.
      typename Dist::value_type z, // location of random variable z to give probability, P(X > z) == p.
      // For example, a nominal minimum acceptable z, so that p * 100 % are > z
      typename Dist::value_type p, // probability value desired at x, say 0.95 for 95% > z.
      typename Dist::value_type location, // location parameter, for example, normal mean.
      const Policy& pol 
      )
    {
      BOOST_STATIC_ASSERT(::boost::math::tools::is_distribution<Dist>::value); 
      BOOST_STATIC_ASSERT(::boost::math::tools::is_scaled_distribution<Dist>::value); 
      static const char* function = "boost::math::find_scale<%1%>&, %1%)";

      if(!(boost::math::isfinite)(p) || (p < 0) || (p > 1))
      {
       return policies::raise_domain_error<typename Dist::value_type>(
           function, "Probability parameter was %1%, but must be >= 0 and <= 1!", p, pol);
      }
      if(!(boost::math::isfinite)(z))
      {
       return policies::raise_domain_error<typename Dist::value_type>(
           function, "z parameter was %1%, but must be finite!", z, pol);
      }
      if(!(boost::math::isfinite)(location))
      {
       return policies::raise_domain_error<typename Dist::value_type>(
           function, "location parameter was %1%, but must be finite!", location, pol);
      }
        
      //cout << "z " << z << ", p " << p << ",  quantile(Dist(), p) "
      //<< quantile(Dist(), p) << ", x - mean " << z - location 
      //<<", sd " << (z - location)  / quantile(Dist(), p) << endl;

       return (z - location)  // difference between desired x and current.
         / quantile(Dist(), p);

    } // find_scale

    template <class Dist>
    inline // with default policy.
      typename Dist::value_type find_scale( // For example, normal mean.
      typename Dist::value_type z, // location of random variable z to give probability, P(X > z) == p.
      // For example, a nominal minimum acceptable z, so that p * 100 % are > z
      typename Dist::value_type p, // probability value desired at x, say 0.95 for 95% > z.
      typename Dist::value_type location) // location parameter, for example, mean.
    { // Forward to find_scale with default policy.
       return (find_scale<Dist>(z, p, location, policies::policy<>()));
    } // find_scale

    // So the user can start from the complement q = (1 - p) of the probability p,
    // for example, s = find_scale<normal>(complement(z, q, l));

    template <class Dist, class Real1, class Real2, class Real3>
    inline typename Dist::value_type find_scale(
      complemented3_type<Real1, Real2, Real3> const& c)
    {
      //cout << "cparam1 q " << c.param1 // q
      //  << ", c.dist z " << c.dist // z
      //  << ", c.param2 l " << c.param2 // l
      //  << ", quantile (Dist(), c.param1 = q) "
      //  << quantile(Dist(), c.param1) //q
      //  << endl;

      return -(c.dist - c.param2) / quantile(Dist(), c.param1);
      //     (  z    - location) / (quantile(Dist(),  p) 
    }

  } // namespace boost
} // namespace math

#endif // BOOST_STATS_FIND_SCALE_HPP
