// Copyright Nick Thompson, 2017
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_QUADRATURE_DETAIL_TANH_SINH_DETAIL_HPP
#define BOOST_MATH_QUADRATURE_DETAIL_TANH_SINH_DETAIL_HPP

#include <cmath>
#include <vector>
#include <typeinfo>
#include <boost/math/constants/constants.hpp>
#include <boost/math/quadrature/detail/tanh_sinh_constants.hpp>
#include <boost/math/special_functions/next.hpp>


namespace boost{ namespace math{ namespace quadrature { namespace detail{


// Returns the tanh-sinh quadrature of a function f over the open interval (-1, 1)

template<class Real, class Policy>
class tanh_sinh_detail 
   : public tanh_sinh_detail_constants<Real, 
       !std::numeric_limits<Real>::is_specialized || (std::numeric_limits<Real>::radix != 2) ? 0 :
       (std::numeric_limits<Real>::digits <= 24) ? 24 :
       (std::numeric_limits<Real>::digits <= 64) ? 64 : 
#ifdef BOOST_HAS_FLOAT128
       (std::numeric_limits<Real>::digits <= 113) ? 113 : 
#endif
       (std::numeric_limits<Real>::digits <= 334) ? 334 :
       0>
{
   typedef tanh_sinh_detail_constants<Real, !std::numeric_limits<Real>::is_specialized || (std::numeric_limits<Real>::radix != 2) ? 0 :
      std::numeric_limits<Real>::digits <= std::numeric_limits<float>::digits ? 24 : 
      (std::numeric_limits<Real>::digits <= 64) ? 64 : 
#ifdef BOOST_HAS_FLOAT128
      (std::numeric_limits<Real>::digits <= 113) ? 113 :
#endif
      (std::numeric_limits<Real>::digits <= 334) ? 334 :
      0> base_type;
#if !BOOST_WORKAROUND(BOOST_MSVC, < 1900)
   using base_type::m_abscissas;
   using base_type::m_weights;
#endif
public:
    tanh_sinh_detail(Real tol, size_t max_refinements);

    template<class F>
    Real integrate(const F f, Real* error, Real* L1, const char* function) const;

private:
    Real m_tol;
    Real m_t_max;
    size_t m_max_refinements;
};

template<class Real, class Policy>
tanh_sinh_detail<Real, Policy>::tanh_sinh_detail(Real tol, size_t max_refinements)
{
    m_tol = tol;
    m_max_refinements = max_refinements;
    /*
     * Our goal is to calculate t_max such that tanh(pi/2 sinh(t_max)) < 1 in the requested precision.
     * What follows is a good estimate for t_max, but in fact we can get closer by starting with this estimate
     * and then calculating tanh(pi/2 sinh(t_max + eps)) until it = 1 (the second to last is t_max).
     * However, this is computationally expensive, so we can't do it.
     * An alternative is to cache the value of t_max for various types (float, double, long double, float128, cpp_bin_float_50, cpp_bin_float_100)
     * and then simply use them, but unfortunately the values will be platform dependent.
     * As such we are then susceptible to the catastrophe where we evaluate the function at x = 1, when we have promised we wouldn't do that.
     * Other authors solve this problem by computing the abscissas in double the requested precision, and then returning the result at the request precision.
     * This seems to be overkill to me, but presumably it's an option if we find integrals on which this method struggles.
     */

     using std::tanh;
     using std::sinh;
     using std::asinh;
     using std::atanh;
     using boost::math::constants::half_pi;
     using boost::math::constants::two_div_pi;

     auto g = [](Real t) { return tanh(half_pi<Real>()*sinh(t)); };
     auto g_inv = [](Real x) { return asinh(two_div_pi<Real>()*atanh(x)); };

     Real x = float_prior((Real) 1);
     m_t_max = g_inv(x);
     while(!(boost::math::isnormal)(m_t_max))
     {
         // Although you might suspect that we could be in this loop essentially for ever, in point of fact it is only called once
         // even for 100 digit accuracy, and isn't called at all up to float128.
         x = float_prior(x);
         m_t_max = g_inv(x);
     }
     // This occurs once on 100 digit arithmetic:
     while(!(g(m_t_max) < (Real) 1))
     {
         x = float_prior(x);
         m_t_max = g_inv(x);
     }
}


template<class Real, class Policy>
template<class F>
Real tanh_sinh_detail<Real, Policy>::integrate(const F f, Real* error, Real* L1, const char* function) const
{
    using std::abs;
    using std::floor;
    using std::tanh;
    using std::sinh;
    using std::sqrt;
    using boost::math::constants::half;
    using boost::math::constants::half_pi;
    Real h = 1;
    Real I0 = half_pi<Real>()*f(0);
    Real L1_I0 = abs(I0);
    for(size_t i = 1; i <= m_t_max; ++i)
    {
        Real x = m_abscissas[0][i];
        Real w = m_weights[0][i];
        Real yp = f(x);
        Real ym = f(-x);
        I0 += (yp + ym)*w;
        L1_I0 += (abs(yp) + abs(ym))*w;
    }
    // Uncomment to understand the convergence rate:
    // std::cout << std::setprecision(std::numeric_limits<Real>::digits10);
    // std::cout << "First estimate : " << I0 << std::endl;
    Real I1 = half<Real>()*I0;
    Real L1_I1 = half<Real>()*L1_I0;
    h /= 2;
    Real sum = 0;
    Real absum = 0;
    for(size_t i = 0; i < m_weights[1].size() && h + i <= m_t_max; ++i)
    {
        Real x = m_abscissas[1][i];
        Real w = m_weights[1][i];

        Real yp = f(x);
        Real ym = f(-x);
        sum += (yp + ym)*w;
        absum += (abs(yp) + abs(ym))*w;
    }
    I1 += sum*h;
    L1_I1 += absum*h;
    Real err = abs(I0 - I1);

    // std::cout << "Second estimate: " << I1 << " Error estimate at level " << 1 << " = " << err << std::endl;

    size_t k = 2;

    while (k < 4 || (k < m_weights.size() && k < m_max_refinements) )
    {
        I0 = I1;
        L1_I0 = L1_I1;

        I1 = half<Real>()*I0;
        L1_I1 = half<Real>()*L1_I0;
        h *= half<Real>();
        sum = 0;
        absum = 0;
        auto const& abscissa_row = m_abscissas[k];
        auto const& weight_row = m_weights[k];

        // We have no guarantee that round-off error won't cause the function to be evaluated strictly at the endpoints.
        // However, making sure x < 1 - eps is a reasonable compromise between accuracy and usability..
        for(size_t j = 0; j < weight_row.size() && abscissa_row[j] < (Real) 1 - tools::epsilon<Real>(); ++j)
        {
            Real x = abscissa_row[j];
            Real w = weight_row[j];

            Real yp = f(x);
            Real ym = f(-x);
            Real term = (yp + ym)*w;
            sum += term;

            // A question arises as to how accurately we actually need to estimate the L1 integral.
            // For simple integrands, computing the L1 norm makes the integration 20% slower,
            // but for more complicated integrands, this calculation is not noticeable.
            Real abterm = (abs(yp) + abs(ym))*w;
            absum += abterm;
        }

        I1 += sum*h;
        L1_I1 += absum*h;
        ++k;
        err = abs(I0 - I1);
        // std::cout << "Estimate:        " << I1 << " Error estimate at level " << k  << " = " << err << std::endl;

        if (!(boost::math::isfinite)(I1))
        {
            return policies::raise_evaluation_error(function, "The tanh_sinh quadrature evaluated your function at a singular point at got %1%. Please narrow the bounds of integration or check your function for singularities.", I1, Policy());
        }

        if (err <= m_tol*L1_I1)
        {
            break;
        }

    }

    // Since we are out of precomputed abscissas and weights, switch to computing abscissas and weights on the fly.
    while (k < m_max_refinements && err > m_tol*L1_I1)
    {
        I0 = I1;
        L1_I0 = L1_I1;

        I1 = half<Real>()*I0;
        L1_I1 = half<Real>()*L1_I0;
        h *= half<Real>();
        sum = 0;
        absum = 0;
        for(Real t = h; t < m_t_max - tools::epsilon<Real>(); t += 2*h)
        {
            Real s = sinh(t);
            Real c = sqrt(1+s*s);
            Real x = tanh(half_pi<Real>()*s);
            Real w = half_pi<Real>()*c*(1-x*x);

            Real yp = f(x);
            Real ym = f(-x);
            Real term = (yp + ym)*w;
            sum += term;
            Real abterm = (abs(yp) + abs(ym))*w;
            absum += abterm;
            // There are claims that this test improves performance,
            // however my benchmarks show that it's slower!
            // However, I leave this comment here because it totally stands to reason that this *should* help:
            //if (abterm < std::numeric_limits<Real>::epsilon()) { break; }
        }

        I1 += sum*h;
        L1_I1 += absum*h;
        ++k;
        err = abs(I0 - I1);

        // std::cout << "Estimate:        " << I1 << " Error estimate at level " << k  << " = " << err << std::endl;
        if (!(boost::math::isfinite)(I1))
        {
            return policies::raise_evaluation_error(function, "The tanh_sinh quadrature evaluated your function at a singular point at got %1%. Please narrow the bounds of integration or check your function for singularities.", I1, Policy());
        }

        if (err <= m_tol*L1_I1)
        {
            break;
        }
    }

    if (error)
    {
        *error = err;
    }

    if (L1)
    {
        *L1 = L1_I1;
    }

    return I1;
}

}}}}
#endif
