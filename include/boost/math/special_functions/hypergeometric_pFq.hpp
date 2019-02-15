
///////////////////////////////////////////////////////////////////////////////
//  Copyright 2018 John Maddock
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_HYPERGEOMETRIC_PFQ_HPP
#define BOOST_MATH_HYPERGEOMETRIC_PFQ_HPP

#include <boost/math/special_functions/detail/hypergeometric_pFq_checked_series.hpp>

namespace boost {
   namespace math {

      template <class Seq, class Real, class Policy>
      inline Real hypergeometric_pFq(const Seq& aj, const Seq& bj, const Real& z, Real* pNorm, const Policy& pol)
      {
         int scale = 0;
         std::pair<Real, Real> r = boost::math::detail::hypergeometric_pFq_checked_series_impl(aj, bj, z, pol, boost::math::detail::iteration_terminator(boost::math::policies::get_max_series_iterations<Policy>()), scale);
         r.first *= exp(Real(scale));
         r.second *= exp(Real(scale));
         if (pNorm)
            *pNorm = r.second;
         return r.first;
      }

      template <class Seq, class Real>
      inline std::pair<Real, Real> hypergeometric_pFq(const Seq& aj, const Seq& bj, const Real& z, Real* pNorm = 0)
      {
         return hypergeometric_pFq(aj, bj, z, pNorm, boost::math::policies::policy<>());
      }

      template <class Seq, class Real>
      Real hypergeometric_pFq_precision(const Seq& aj, const Seq& bj, const Real& z, unsigned digits10)
      {
         unsigned current_precision = digits10 + 5;

         for (auto ai = aj.begin(); ai != aj.end(); ++ai)
         {
            current_precision = std::max(current_precision, ai->precision());
         }
         for (auto bi = bj.begin(); bi != bj.end(); ++bi)
         {
            current_precision = std::max(current_precision, bi->precision());
         }
         current_precision = std::max(current_precision, z.precision());

         Real r, norm;
         std::vector<Real> aa(aj), bb(bj);
         do
         {
            Real::default_precision(current_precision);
            for (auto ai = aa.begin(); ai != aa.end(); ++ai)
               ai->precision(current_precision);
            for (auto bi = bb.begin(); bi != bb.end(); ++bi)
               bi->precision(current_precision);
            try
            {
               r = hypergeometric_pFq(aa, bb, z, &norm, boost::math::policies::policy<>());

               unsigned precision_obtained = current_precision - 1 - itrunc(log10(abs(norm / r)));
               if (precision_obtained < digits10)
               {
                  current_precision += digits10 - precision_obtained + 5;
               }
               else
                  break;
            }
            catch (const boost::math::evaluation_error&)
            {
               current_precision *= 2;
            }
         } while (true);

         return r;
      }

      template <class Real>
      Real hypergeometric_pFq_precision(const std::initializer_list<Real>& aj, const std::initializer_list<Real>& bj, const Real& z, unsigned digits10)
      {
         return hypergeometric_pFq_precision< std::initializer_list<Real>, Real>(aj, bj, z, digits10);
      }

   }
} // namespaces

#endif // BOOST_MATH_BESSEL_ITERATORS_HPP
