// Copyright Nick Thompson, 2017
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

/*
 * This class performs exp-sinh quadrature on half infinite intervals.
 *
 * References:
 *
 * 1) Tanaka, Ken’ichiro, et al. "Function classes for double exponential integration formulas." Numerische Mathematik 111.4 (2009): 631-655.
 */

#ifndef BOOST_MATH_QUADRATURE_EXP_SINH_HPP
#define BOOST_MATH_QUADRATURE_EXP_SINH_HPP

#include <cmath>
#include <limits>
#include <memory>
#include <boost/math/quadrature/detail/exp_sinh_detail.hpp>

namespace boost{ namespace math{ namespace quadrature {


template<class Real>
class exp_sinh
{
public:
    exp_sinh(Real tol = sqrt(std::numeric_limits<Real>::epsilon()), size_t max_refinements = 9);

    template<class F>
    Real integrate(const F& f, Real a = 0, Real b = std::numeric_limits<Real>::infinity(), Real* error = nullptr, Real* L1 = nullptr) const;

private:
    std::shared_ptr<detail::exp_sinh_detail<Real>> m_imp;
};

template<class Real>
exp_sinh<Real>::exp_sinh(Real tol, size_t max_refinements) : m_imp(std::make_shared<detail::exp_sinh_detail<Real>>(tol, max_refinements))
{
    return;
}

template<class Real>
template<class F>
Real exp_sinh<Real>::integrate(const F& f, Real a, Real b, Real* error, Real* L1) const
{
    using std::isfinite;
    using std::abs;
    using boost::math::constants::half;
    using boost::math::quadrature::detail::exp_sinh_detail;

    // Right limit is infinite:
    if (isfinite(a) && b >= std::numeric_limits<Real>::max())
    {
        // If a = 0, don't use an additional level of indirection:
        if (a == (Real) 0)
        {
            return m_imp->integrate(f, error, L1);
        }
        const auto u = [&](Real t) { return f(t + a); };
        return m_imp->integrate(u, error, L1);
    }

    if (isfinite(b) && a <= std::numeric_limits<Real>::lowest())
    {
        const auto u = [&](Real t) { return f(b-t);};
        return m_imp->integrate(u, error, L1);
    }

    if (isfinite(a) && isfinite(b))
    {
        throw std::domain_error("Use tanh_sinh quadrature for integration over finite domains; exp_sinh is for half infinite integrals.\n");
    }

    // Infinite limits:
    if (a <= std::numeric_limits<Real>::lowest() && b >= std::numeric_limits<Real>::max())
    {
        throw std::domain_error("Use sinh_sinh quadrature for integration over the whole real line; exp_sinh is for half infinite integrals.\n");
    }

    throw std::domain_error("The domain of integration is not sensible; please check the bounds.\n");
}


}}}
#endif
