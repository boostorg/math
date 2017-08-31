//  Copyright John Maddock 2017.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_QUADRATURE_GAUSS_KRONROD_HPP
#define BOOST_MATH_QUADRATURE_GAUSS_KRONROD_HPP

#ifdef _MSC_VER
#pragma once
#pragma warning(push)
#pragma warning(disable: 4127)
#endif

#include <vector>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/legendre_stieltjes.hpp>
#include <boost/math/quadrature/gauss.hpp>

namespace boost { namespace math{ namespace quadrature{ namespace detail{

template <class Real, unsigned N>
class gauss_kronrod_detail
{
   static legendre_stieltjes<Real> const& get_legendre_stieltjes()
   {
      static const legendre_stieltjes<Real> data((N - 1) / 2 + 1);
      return data;
   }
   static std::vector<Real> calculate_abscissa()
   {
      static std::vector<Real> result = boost::math::legendre_p_zeros<Real>((N - 1) / 2);
      const legendre_stieltjes<Real> E = get_legendre_stieltjes();
      std::vector<Real> ls_zeros = E.zeros();
      result.insert(result.end(), ls_zeros.begin(), ls_zeros.end());
      std::sort(result.begin(), result.end());
      return result;
   }
   static std::vector<Real> calculate_weights()
   {
      std::vector<Real> result(abscissa().size(), 0);
      unsigned gauss_order = (N - 1) / 2;
      unsigned gauss_start = gauss_order & 1 ? 0 : 1;
      const legendre_stieltjes<Real>& E = get_legendre_stieltjes();

      for (unsigned i = gauss_start; i < abscissa().size(); i += 2)
      {
         Real x = abscissa()[i];
         Real p = boost::math::legendre_p_prime(gauss_order, x);
         Real gauss_weight = 2 / ((1 - x * x) * p * p);
         result[i] = gauss_weight + static_cast<Real>(2) / (static_cast<Real>(gauss_order + 1) * legendre_p_prime(gauss_order, x) * E(x));
      }
      for (unsigned i = gauss_start ? 0 : 1; i < abscissa().size(); i += 2)
      {
         Real x = abscissa()[i];
         result[i] = static_cast<Real>(2) / (static_cast<Real>(gauss_order + 1) * legendre_p(gauss_order, x) * E.prime(x));
      }
      return result;
   }
public:
   static const std::vector<Real>& abscissa()
   {
      static std::vector<Real> data = calculate_abscissa();
      return data;
   }
   static const std::vector<Real>& weights()
   {
      static std::vector<Real> data = calculate_weights();
      return data;
   }
};

template <>
class gauss_kronrod_detail<double, 15>
{
public:
   static constexpr std::array<double, 8> const & abscissa()
   {
      static constexpr std::array<double, 8> data = {
         0.000000000000000000000000000000000e+00,
         2.077849550078984676006894037732449e-01,
         4.058451513773971669066064120769615e-01,
         5.860872354676911302941448382587296e-01,
         7.415311855993944398638647732807884e-01,
         8.648644233597690727897127886409262e-01,
         9.491079123427585245261896840478513e-01,
         9.914553711208126392068546975263285e-01,
      };
      return data;
   }
   static constexpr std::array<double, 8> const & weights()
   {
      static constexpr std::array<double, 8> data = {
         2.094821410847278280129991748917143e-01,
         2.044329400752988924141619992346491e-01,
         1.903505780647854099132564024210137e-01,
         1.690047266392679028265834265985503e-01,
         1.406532597155259187451895905102379e-01,
         1.047900103222501838398763225415180e-01,
         6.309209262997855329070066318920429e-02,
         2.293532201052922496373200805896959e-02,
      };
      return data;
   }
};

}

template <class Real, unsigned N, class Policy = boost::math::policies::policy<> >
class gauss_kronrod : public detail::gauss_kronrod_detail<Real, N>
{
public:
   typedef Real value_type;
private:
   template <class F>
   static value_type integrate_non_adaptive_m1_1(F f, Real* error = nullptr, Real* pL1 = nullptr)
   {
      using std::fabs;
      unsigned gauss_start = 2;
      unsigned kronrod_start = 1;
      unsigned gauss_order = (N - 1) / 2;
      value_type kronrod_result = 0;
      value_type gauss_result = 0;
      value_type fp, fm;
      if (gauss_order & 1)
      {
         fp = f(value_type(0));
         kronrod_result = fp * weights()[0];
         gauss_result += fp * gauss<Real, (N - 1) / 2>::weights()[0];
      }
      else
      {
         fp = f(value_type(0));
         kronrod_result = fp * weights()[0];
         gauss_start = 1;
         kronrod_start = 2;
      }
      value_type L1 = fabs(kronrod_result);
      for (unsigned i = gauss_start; i < abscissa().size(); i += 2)
      {
         fp = f(abscissa()[i]);
         fm = f(-abscissa()[i]);
         kronrod_result += (fp + fm) * weights()[i];
         L1 += (fabs(fp) + fabs(fm)) *  weights()[i];
         gauss_result += (fp + fm) * gauss<Real, (N - 1) / 2>::weights()[i / 2];
      }
      for (unsigned i = kronrod_start; i < abscissa().size(); i += 2)
      {
         fp = f(abscissa()[i]);
         fm = f(-abscissa()[i]);
         kronrod_result += (fp + fm) * weights()[i];
         L1 += (fabs(fp) + fabs(fm)) *  weights()[i];
      }
      if (pL1)
         *pL1 = L1;
      if (error)
         *error = fabs(kronrod_result - gauss_result);
      return kronrod_result;
   }

   template <class F>
   struct recursive_info
   {
      F f;
      Real tol;
   };

   template <class F>
   static value_type recursive_adaptive_integrate(const recursive_info<F>* info, Real a, Real b, unsigned max_levels, Real abs_tol, Real* error, Real* L1)
   {
      Real error_local;
      Real mean = (b + a) / 2;
      Real scale = (b - a) / 2;
      auto ff = [&](const Real& x)
      {
         return info->f(scale * x + mean);
      };
      Real estimate = scale * integrate_non_adaptive_m1_1(ff, &error_local, L1);

      Real abs_tol1 = fabs(estimate * info->tol);
      if (abs_tol == 0)
         abs_tol = abs_tol1;

      if (max_levels && (abs_tol1 < error_local) && (abs_tol < error_local))
      {
         Real mid = (a + b) / 2;
         Real L1_local;
         estimate = recursive_adaptive_integrate(info, a, mid, max_levels - 1, abs_tol / 2, error, L1);
         estimate += recursive_adaptive_integrate(info, mid, b, max_levels - 1, abs_tol / 2, &error_local, &L1_local);
         if (error)
            *error += error_local;
         if (L1)
            *L1 += L1_local;
         return estimate;
      }
      if(L1)
         *L1 *= scale;
      if (error)
         *error = error_local;
      return estimate;
   }

public:
   template <class F>
   static value_type integrate(F f, Real a, Real b, unsigned max_depth = 15, Real tol = tools::root_epsilon<Real>(), Real* error = nullptr, Real* pL1 = nullptr)
   {
      static const char* function = "boost::math::quadrature::gauss_kronrod<%1%>::integrate(f, %1%, %1%)";
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
            recursive_info<decltype(u)> info = { u, tol };
            return recursive_adaptive_integrate(&info, Real(-1), Real(1), max_depth, Real(0), error, pL1);
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
            recursive_info<decltype(u)> info = { u, tol };
            Real Q = 2 * recursive_adaptive_integrate(&info, Real(-1), Real(1), max_depth, Real(0), error, pL1);
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
            recursive_info<decltype(v)> info = { v, tol };
            Real Q = 2 * recursive_adaptive_integrate(&info, Real(-1), Real(1), max_depth, Real(0), error, pL1);
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
            recursive_info<F> info = { f, tol };
            return recursive_adaptive_integrate(&info, a, b, max_depth, Real(0), error, pL1);
         }
      }
      return policies::raise_domain_error(function, "The domain of integration is not sensible; please check the bounds.", a, Policy());
   }
};

} // namespace quadrature
} // namespace math
} // namespace boost

#ifdef _MSC_VER
#pragma warning(pop)
#endif

#endif // BOOST_MATH_QUADRATURE_GAUSS_KRONROD_HPP

