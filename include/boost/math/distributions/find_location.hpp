//  Copyright John Maddock 2007.
//  Copyright Paul A. Bristow 2007.

//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_STATS_FIND_LOCATION_HPP
#define BOOST_STATS_FIND_LOCATION_HPP

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
  // TODO use concepts to enforce this.

    template <class Dist, class Policy>
    inline
      typename Dist::value_type find_location( // For example, normal mean.
      typename Dist::value_type z, // location of random variable z to give probability, P(X > z) == p.
      // For example, a nominal minimum acceptable z, so that p * 100 % are > z
      typename Dist::value_type p, // probability value desired at x, say 0.95 for 95% > z.
      typename Dist::value_type scale, // scale parameter, for example, normal standard deviation.
      const Policy& pol 
      )
    {
      BOOST_STATIC_ASSERT(::boost::math::tools::is_distribution<Dist>::value); 
      BOOST_STATIC_ASSERT(::boost::math::tools::is_scaled_distribution<Dist>::value); 
      static const char* function = "boost::math::find_location<%1%>&, %1%)";

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
      if(!(boost::math::isfinite)(scale))
      {
       return policies::raise_domain_error<typename Dist::value_type>(
           function, "scale parameter was %1%, but must be finite!", scale, pol);
      }
        
      //cout << "z " << z << ", p " << p << ",  quantile(Dist(), p) "
      //  << quantile(Dist(), p) << ", quan * scale " << quantile(Dist(), p) * scale << endl;
      return z - (quantile(Dist(), p) * scale);
    } // find_location

    template <class Dist>
    inline // with default policy.
      typename Dist::value_type find_location( // For example, normal mean.
      typename Dist::value_type z, // location of random variable z to give probability, P(X > z) == p.
      // For example, a nominal minimum acceptable z, so that p * 100 % are > z
      typename Dist::value_type p, // probability value desired at x, say 0.95 for 95% > z.
      typename Dist::value_type scale) // scale parameter, for example, normal standard deviation.
    { // Forward to find_location with default policy.
       return (find_location<Dist>(z, p, scale, policies::policy<>()));
    } // find_location

    // So the user can start from the complement q = (1 - p) of the probability p,
    // for example, l = find_location<normal>(complement(z, q, sd));

    template <class Dist, class Real1, class Real2, class Real3>
    inline typename Dist::value_type find_location(
      complemented3_type<Real1, Real2, Real3> const& c)
    {
       // cout << "z " << c.dist << ", quantile (Dist(), " << c.param1 << ") * scale " << c.param2 << endl;
       return c.dist - quantile(Dist(), c.param1) * c.param2;
       //       z    - (quantile(Dist(),  p)      *   scale
    }

  } // namespace boost
} // namespace math

#endif // BOOST_STATS_FIND_LOCATION_HPP

