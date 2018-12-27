///////////////////////////////////////////////////////////////////////////////
//  Copyright 2018 John Maddock
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_HYPERGEOMETRIC_PFQ_SERIES_HPP_
#define BOOST_HYPERGEOMETRIC_PFQ_SERIES_HPP_

#ifndef BOOST_MATH_PFQ_MAX_B_TERMS
#  define BOOST_MATH_PFQ_MAX_B_TERMS 5
#endif

#include <boost/array.hpp>
#include <boost/math/special_functions/detail/hypergeometric_series.hpp>

  namespace boost { namespace math { namespace detail {

     template <class Seq, class Real, class Policy, class Terminal>
     std::pair<Real, Real> hypergeometric_pFq_checked_series_impl(const Seq& aj, const Seq& bj, const Real& z, const Policy& pol, const Terminal& termination, int& log_scale)
     {
        using std::abs;
        Real result = 1;
        Real abs_result = 1;
        Real term = 1;
        Real term0 = 0;
        Real tol = boost::math::policies::get_epsilon<Real, Policy>();
        boost::uintmax_t k = 0;
        Real upper_limit(sqrt(boost::math::tools::max_value<Real>())), diff;
        Real lower_limit(1 / upper_limit);
        int log_scaling_factor = boost::math::itrunc(boost::math::tools::log_max_value<Real>()) - 2;
        Real scaling_factor = exp(log_scaling_factor);
        Real term_m1 = 0;
        int local_scaling = 0;

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
           //std::cout << k << " " << *bj.begin() + k << " " << result << " " << term << /*" " << term_at_k(*aj.begin(), *bj.begin(), z, k, pol) <<*/ std::endl;
           result += term;
           abs_result += abs(term);
           //std::cout << result << std::endl;
           //std::cout << abs_result << std::endl;

           //
           // Rescaling:
           //
           if (fabs(abs_result) >= upper_limit)
           {
              abs_result /= scaling_factor;
              result /= scaling_factor;
              term /= scaling_factor;
              log_scale += log_scaling_factor;
              local_scaling += log_scaling_factor;
           }
           if (fabs(abs_result) < lower_limit)
           {
              abs_result *= scaling_factor;
              result *= scaling_factor;
              term *= scaling_factor;
              log_scale -= log_scaling_factor;
              local_scaling -= log_scaling_factor;
           }

           if ((abs(result * tol) > abs(term)) && (abs(term0) > abs(term)))
              break;
           if (abs_result * tol > abs(result))
           {
              // We have no correct bits in the result... just give up!
              result = boost::math::policies::raise_evaluation_error("boost::math::hypergeometric_pFq<%1%>", "Cancellation is so severe that no bits in the reuslt are correct, last result was %1%", Real(result * exp(Real(log_scale))), pol);
              return std::make_pair(result, result);
           }
           term0 = term;
        }
        //
        // We have to be careful when one of the b's crosses the origin:
        //
        if(bj.size() > BOOST_MATH_PFQ_MAX_B_TERMS)
           policies::raise_domain_error<Real>("boost::math::hypergeometric_pFq<%1%>(Seq, Seq, %1%)", 
              "The number of b terms must be less than the value of BOOST_MATH_PFQ_MAX_B_TERMS (" BOOST_STRINGIZE(BOOST_MATH_PFQ_MAX_B_TERMS)  "), but got %1%.",
              Real(bj.size()), pol);

        int crossover_locations[BOOST_MATH_PFQ_MAX_B_TERMS];
        int n = 0;
        for (auto bi = bj.begin(); bi != bj.end(); ++bi, ++n)
        {
           crossover_locations[n] = *bi >= 0 ? 0 : itrunc(-*bi) + 1;
        }
        std::sort(crossover_locations, crossover_locations + bj.size(), std::greater<Real>());

        for (n = 0; n < bj.size(); ++n)
        {
           if (k < crossover_locations[n])
           {
              for(auto ai = aj.begin(); ai != aj.end(); ++ai)
                 if ((*ai < 0) && (floor(*ai) == *ai) && (*ai > crossover_locations[n]))
                    return std::make_pair(result, abs_result);  // b's will never cross the origin!
              //
              // b hasn't crossed the origin yet and the series may spring back into life at that point
              // so we need to jump forward to that term and then evaluate forwards and backwards from there:
              //
              unsigned s = crossover_locations[n];
              boost::uintmax_t backstop = k;
              int s1(1), s2(1);
              term = 0;
              for(auto ai = aj.begin(); ai != aj.end(); ++ai)
                 term += log_pochhammer(*ai, s, pol, &s1);
              for(auto bi = bj.begin(); bi != bj.end(); ++bi)
                 term -= log_pochhammer(*bi, s, pol, &s2);
              term -= lgamma(Real(s + 1), pol);
              term += s * log(fabs(z));
              if (z < 0)
                 s1 *= (s & 1 ? -1 : 1);
              term -= local_scaling;
              if (term > -tools::log_max_value<Real>())
              {
                 if (term > 10)
                 {
                    int scale = itrunc(floor(term));
                    term -= scale;
                    log_scale += scale;
                    result *= exp(Real(-scale));
                 }
                 term = s1 * s2 * exp(term);
                 k = s;
                 term0 = term;
                 do
                 {
                    result += term;
                    abs_result += fabs(term);
                    if (abs_result * tol > abs(result))
                    {
                       // We have no correct bits in the result... just give up!
                       result = boost::math::policies::raise_evaluation_error("boost::math::hypergeometric_pFq<%1%>", "Cancellation is so severe that no bits in the reuslt are correct, last result was %1%", Real(result * exp(Real(log_scale))), pol);
                       return std::make_pair(result, result);
                    }
                    if (fabs(result) >= upper_limit)
                    {
                       result /= scaling_factor;
                       abs_result /= scaling_factor;
                       term /= scaling_factor;
                       log_scale += log_scaling_factor;
                       term0 /= scaling_factor;
                    }
                    if (fabs(result) < lower_limit)
                    {
                       result *= scaling_factor;
                       abs_result *= scaling_factor;
                       term *= scaling_factor;
                       log_scale -= log_scaling_factor;
                       term0 *= scaling_factor;
                    }
                    term_m1 = term;
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
                    term *= z / (k + 1);

                    ++k;
                    diff = fabs(term / result);
                 } while ((diff > boost::math::policies::get_epsilon<Real, Policy>()) || (fabs(term) > fabs(term_m1)));
                 //
                 // Now go backwards as well:
                 //
                 k = s;
                 term = term0;
                 do
                 {
                    --k;
                    if (k == backstop)
                       break;
                    term_m1 = term;
                    for (auto ai = aj.begin(); ai != aj.end(); ++ai)
                    {
                       term /= *ai + k;
                    }
                    for (auto bi = bj.begin(); bi != bj.end(); ++bi)
                    {
                       if (*bi + k == 0)
                       {
                          // The series is undefined:
                          result = boost::math::policies::raise_domain_error("boost::math::hypergeometric_pFq<%1%>", "One of the b values was the negative integer %1%", *bi, pol);
                          return std::make_pair(result, result);
                       }
                       term *= *bi + k;
                    }
                    term *= (k + 1) / z;
                    result += term;
                    abs_result += fabs(term);
                    if (abs_result * tol > abs(result))
                    {
                       // We have no correct bits in the result... just give up!
                       result = boost::math::policies::raise_evaluation_error("boost::math::hypergeometric_pFq<%1%>", "Cancellation is so severe that no bits in the reuslt are correct, last result was %1%", Real(result * exp(Real(log_scale))), pol);
                       return std::make_pair(result, result);
                    }
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
                    diff = fabs(term / result);
                 } while ((diff > boost::math::policies::get_epsilon<Real, Policy>()) || (fabs(term) > fabs(term_m1)));
              }
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
           return boost::math::policies::raise_evaluation_error("boost::math::hypergeometric_pFq<%1%>", "Cancellation is so severe that fewer than half the bits in the result are correct, last result was %1%", result.first * exp(log_scale), pol);
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
