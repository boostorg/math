//  Copyright John Maddock 2015.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_QUADRATURE_GAUSS_HPP
#define BOOST_MATH_QUADRATURE_GAUSS_HPP

#ifdef _MSC_VER
#pragma once
#endif

#include <vector>
#include <boost/math/special_functions/legendre.hpp>

namespace boost { namespace math{ namespace quadrature{ namespace detail{

template <class Real, unsigned N>
class gauss_detail
{
   static std::vector<Real> calculate_weights()
   {
      std::vector<Real> result(abscissa().size(), 0);
      for (unsigned i = 0; i < abscissa().size(); ++i)
      {
         Real x = abscissa()[i];
         Real p = boost::math::legendre_p_prime(N, x);
         result[i] = 2 / ((1 - x * x) * p * p);
      }
      return result;
   }
public:
   static const std::vector<Real>& abscissa()
   {
      static std::vector<Real> data = boost::math::legendre_p_zeros<Real>(N);
      return data;
   }
   static const std::vector<Real>& weights()
   {
      static std::vector<Real> data = calculate_weights();
      return data;
   }
};

template <>
class gauss_detail<double, 7>
{
public:
   static constexpr std::array<double, 4> const & abscissa()
   {
      static constexpr std::array<double, 4> data = {
         0.000000000000000000000000000000000e+00,
         4.058451513773971669066064120769615e-01,
         7.415311855993944398638647732807884e-01,
         9.491079123427585245261896840478513e-01,
      };
      return data;
   }
   static constexpr std::array<double, 4> const & weights()
   {
      static constexpr std::array<double, 4> data = {
         4.179591836734693877551020408163265e-01,
         3.818300505051189449503697754889751e-01,
         2.797053914892766679014677714237796e-01,
         1.294849661688696932706114326790820e-01,
      };
      return data;
   }
};

}

template <class Real, unsigned N, class Policy = boost::math::policies::policy<> >
class gauss : public detail::gauss_detail<Real, N>
{
public:
   typedef Real value_type;

   template <class F>
   static value_type integrate(F f, Real* pL1 = nullptr)
   {
      using std::fabs;
      unsigned non_zero_start = 1;
      value_type result = 0;
      if (N & 1)
         result = f(value_type(0)) * weights()[0];
      else
         non_zero_start = 0;
      value_type L1 = fabs(result);
      for (unsigned i = non_zero_start; i < abscissa().size(); ++i)
      {
         value_type fp = f(abscissa()[i]);
         value_type fm = f(-abscissa()[i]);
         result += (fp + fm) * weights()[i];
         L1 += (fabs(fp) + fabs(fm)) *  weights()[i];
      }
      if (pL1)
         *pL1 = L1;
      return result;
   }
   template <class F>
   static value_type integrate(F f, Real a, Real b, Real* pL1 = nullptr)
   {
      static const char* function = "boost::math::quadrature::gauss<%1%>::integrate(f, %1%, %1%)";
      if (!(boost::math::isnan)(a) && !(boost::math::isnan)(b))
      {
         // Infinite limits:
         if ((a <= -tools::max_value<Real>()) && (b >= tools::max_value<Real>()))
         {
            auto u = [&](const Real& t)->Real
            {
               Real t_sq = t*t;
               Real inv = 1 / (1 - t_sq);
               return f(t*inv)*(1 + t_sq)*inv*inv;
            };
            return integrate(u, pL1);
         }

         // Right limit is infinite:
         if ((boost::math::isfinite)(a) && (b >= tools::max_value<Real>()))
         {
            auto u = [&](const Real& t)->Real
            {
               Real z = 1 / (t + 1);
               Real arg = 2 * z + a - 1;
               return f(arg)*z*z;
            };
            Real Q = 2 * integrate(u, pL1);
            if (pL1)
            {
               *pL1 *= 2;
            }
            return Q;
         }

         if ((boost::math::isfinite)(b) && (a <= -tools::max_value<Real>()))
         {
            auto v = [&](const Real& t)->Real
            {
               Real z = 1 / (t + 1);
               Real arg = 2 * z - 1;
               return f(b - arg) * z * z;
            };
            Real Q = 2 * integrate(v, pL1);
            if (pL1)
            {
               *pL1 *= 2;
            }
            return Q;
         }

         if ((boost::math::isfinite)(a) && (boost::math::isfinite)(b))
         {
            if (b <= a)
            {
               return policies::raise_domain_error(function, "Arguments to integrate are in wrong order; integration over [a,b] must have b > a.", a, Policy());
            }
            Real avg = (a + b)*half<Real>();
            Real scale = (b - a)*half<Real>();

            auto u = [&](Real z)->Real
            {
               return f(avg + scale*z);
            };
            Real Q = scale*integrate(u, pL1);

            if (pL1)
            {
               *pL1 *= scale;
            }
            return Q;
         }
      }
      return policies::raise_domain_error(function, "The domain of integration is not sensible; please check the bounds.", a, Policy());
   }
};

} // namespace quadrature
} // namespace math
} // namespace boost

#endif // BOOST_MATH_QUADRATURE_GAUSS_HPP

