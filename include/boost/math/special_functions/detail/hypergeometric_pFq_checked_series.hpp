///////////////////////////////////////////////////////////////////////////////
//  Copyright 2018 John Maddock
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_HYPERGEOMETRIC_PFQ_SERIES_HPP_
#define BOOST_HYPERGEOMETRIC_PFQ_SERIES_HPP_

#include <boost/array.hpp>

  namespace boost { namespace math { namespace detail {

     template <class Seq, class Real, class Policy, class Terminal>
     std::pair<Real, Real> hypergeometric_pFq_checked_series_impl(const Seq& aj, const Seq& bj, const Real& z, const Policy& pol, const Terminal& termination, int& log_scale)
     {
        using std::abs;
        Real result = 1;
        Real abs_result = 1;
        Real term = 1;
        Real tol = boost::math::policies::get_epsilon<Real, Policy>();
        boost::uintmax_t k = 0;
        int log_scaling_factor = boost::math::itrunc(boost::math::tools::log_max_value<Real>()) - 2;
        Real scaling_factor = exp(log_scaling_factor);
        Real upper_limit(sqrt(boost::math::tools::max_value<Real>()));
        Real lower_limit(1 / upper_limit);

        while (!termination(k))
        {
           for (auto ai = aj.begin(); ai != aj.end(); ++ai)
           {
              term *= *ai + k;
           }
           if (term == 0)
           {
              // There is a negative integer in the aj's:
              return std::make_pair(result, abs_result);
           }
           for (auto bi = bj.begin(); bi != bj.end(); ++bi)
           {
              if (*bi + k == 0)
              {
                 // The series is undefined:
                 result = boost::math::policies::raise_domain_error("boost::math::hypergeometric_pFq<%1%>", "One of the b values was the negative integer %1%", *bi, pol);
                 return std::make_pair(result, result);
              }
              term /= *bi + k;
           }
           term *= z;
           ++k;
           term /= k;
           result += term;
           abs_result += abs(term);

           if (fabs(result) >= upper_limit)
           {
              result /= scaling_factor;
              abs_result /= scaling_factor;
              term /= scaling_factor;
              log_scale += log_scaling_factor;
           }
           if (fabs(result) < lower_limit)
           {
              result *= scaling_factor;
              abs_result *= scaling_factor;
              term *= scaling_factor;
              log_scale -= log_scaling_factor;
           }

           if (abs(result * tol) > abs(term))
              break;
           if (abs_result * tol > abs(result))
           {
              // We have no correct bits in the result... just give up!
              result = boost::math::policies::raise_evaluation_error("boost::math::hypergeometric_pFq<%1%>", "Cancellation is so severe that no bits in the result are correct, last result was %1%", result, pol);
              return std::make_pair(result, result);
           }
        }
        return std::make_pair(result, abs_result);
     }

     struct iteration_terminator
     {
        iteration_terminator(boost::uintmax_t i) : m(i) {}

        bool operator()(boost::uintmax_t v) const { return v >= m; }

        boost::uintmax_t m;
     };

     template <class Seq, class Real, class Policy>
     Real hypergeometric_pFq_checked_series_impl(const Seq& aj, const Seq& bj, const Real& z, const Policy& pol, int& log_scale)
     {
        using std::abs;
        iteration_terminator term(boost::math::policies::get_max_series_iterations<Policy>());
        std::pair<Real, Real> result = hypergeometric_pFq_checked_series_impl(aj, bj, z, pol, term, log_scale);
        //
        // Check to see how many digits we've lost, if it's more than half, raise an evaluation error -
        // this is an entirely arbitrary cut off, but not unreasonable.
        //
        if (result.second * boost::math::policies::get_epsilon<Real, Policy>() > abs(result.first))
        {
           return boost::math::policies::raise_evaluation_error("boost::math::hypergeometric_pFq<%1%>", "Cancellation is so severe that fewer than half the bits in the result are correct, last result was %1%", result.first, pol);
        }
        return result.first;
     }

     template <class Real, class Policy>
     inline Real hypergeometric_1F1_checked_series_impl(const Real& a, const Real& b, const Real& z, const Policy& pol, int& log_scale)
     {
        boost::array<Real, 1> aj = { a };
        boost::array<Real, 1> bj = { b };
        return hypergeometric_pFq_checked_series_impl(aj, bj, z, pol, log_scale);
     }

  } } } // namespaces

#endif // BOOST_HYPERGEOMETRIC_PFQ_SERIES_HPP_
