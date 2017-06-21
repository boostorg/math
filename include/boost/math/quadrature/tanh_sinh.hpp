// Copyright Nick Thompson, 2017
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

/*
 * This class performs tanh-sinh quadrature on the real line.
 * Tanh-sinh quadrature is exponentially convergent for integrands in Hardy spaces,
 * (see https://en.wikipedia.org/wiki/Hardy_space for a formal definition), and is optimal for a random function from that class.
 *
 * The tanh-sinh quadrature is one of a class of so called "double exponential quadratures"-there is a large family of them,
 * but this one seems to be the most commonly used.
 *
 * As always, there are caveats: For instance, if the function you want to integrate is not holomorphic on the unit disk,
 * then the rapid convergence will be spoiled. In this case, a more appropriate quadrature is (say) Romberg, which does not
 * require the function to be holomorphic, only differentiable up to some order.
 *
 * In addition, if you are integrating a periodic function over a period, the trapezoidal rule is better.
 *
 * References:
 *
 * 1) Mori, Masatake. "Quadrature formulas obtained by variable transformation and the DE-rule." Journal of Computational and Applied Mathematics 12 (1985): 119-130.
 * 2) Bailey, David H., Karthik Jeyabalan, and Xiaoye S. Li. "A comparison of three high-precision quadrature schemes." Experimental Mathematics 14.3 (2005): 317-329.
 * 3) Press, William H., et al. "Numerical recipes third edition: the art of scientific computing." Cambridge University Press 32 (2007): 10013-2473.
 *
 */

#ifndef BOOST_MATH_QUADRATURE_TANH_SINH_HPP
#define BOOST_MATH_QUADRATURE_TANH_SINH_HPP

#include <cmath>
#include <limits>
#include <memory>
#include <boost/math/quadrature/detail/tanh_sinh_detail.hpp>

namespace boost{ namespace math{ namespace quadrature {

template<class Real, class Policy = policies::policy<> >
class tanh_sinh
{
public:
    tanh_sinh(Real tol = tools::root_epsilon<Real>(), size_t max_refinements = 15)
    : m_imp(std::make_shared<detail::tanh_sinh_detail<Real, Policy>>(tol, max_refinements)) {}

    template<class F>
    Real integrate(const F f, Real a, Real b, Real* error = nullptr, Real* L1 = nullptr) const;

private:
    std::shared_ptr<detail::tanh_sinh_detail<Real, Policy>> m_imp;
};

template<class Real, class Policy>
template<class F>
Real tanh_sinh<Real, Policy>::integrate(const F f, Real a, Real b, Real* error, Real* L1) const
{
    using boost::math::constants::half;
    using boost::math::quadrature::detail::tanh_sinh_detail;

    static const char* function = "tanh_sinh<%1%>::integrate";

    if ((boost::math::isfinite)(a) && (boost::math::isfinite)(b))
    {
        if (b <= a)
        {
           return policies::raise_domain_error(function, "Arguments to integrate are in wrong order; integration over [a,b] must have b > a.", a, Policy());
        }
        Real avg = (a+b)*half<Real>();
        Real diff = (b-a)*half<Real>();
        auto u = [&](Real z) { return f(avg + diff*z); };
        Real Q = diff*m_imp->integrate(u, error, L1, function);

        if(L1)
        {
            *L1 *= diff;
        }
        return Q;
    }

    // Infinite limits:
    if ((a <= -tools::max_value<Real>()) && (b >= tools::max_value<Real>()))
    {
        auto u = [&](Real t) { auto t_sq = t*t; auto inv = 1/(1 - t_sq); return f(t*inv)*(1+t_sq)*inv*inv; };
        return m_imp->integrate(u, error, L1, function);
    }

    // Right limit is infinite:
    if ((boost::math::isfinite)(a) && (b >= tools::max_value<Real>()))
    {
        auto u = [&](Real t) { auto z = 1/(t+1); auto arg = 2*z + a - 1; return f(arg)*z*z; };
        Real Q = 2*m_imp->integrate(u, error, L1, function);
        if(L1)
        {
            *L1 *= 2;
        }

        return Q;
    }

    if ((boost::math::isfinite)(b) && (a <= -tools::max_value<Real>()))
    {
        auto u = [&](Real t) { return f(b-t);};
        auto v = [&](Real t) { auto z = 1/(t+1); auto arg = 2*z - 1; return u(arg)*z*z; };

        Real Q = 2*m_imp->integrate(v, error, L1, function);
        if (L1)
        {
            *L1 *= 2;
        }
        return Q;
    }

    return policies::raise_domain_error(function, "The domain of integration is not sensible; please check the bounds.", a, Policy());
}


}}}
#endif
