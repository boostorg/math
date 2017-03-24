// Copyright Nick Thompson, 2017
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

/*
 * This class performs sinh-sinh quadrature over the entire real line.
 *
 * References:
 *
 * 1) Tanaka, Ken’ichiro, et al. "Function classes for double exponential integration formulas." Numerische Mathematik 111.4 (2009): 631-655.
 */

#ifndef BOOST_MATH_QUADRATURE_SINH_SINH_HPP
#define BOOST_MATH_QUADRATURE_SINH_SINH_HPP

#include <cmath>
#include <limits>
#include <memory>
#include <boost/math/quadrature/detail/sinh_sinh_detail.hpp>

namespace boost{ namespace math{


template<class Real>
class sinh_sinh
{
public:
    sinh_sinh(Real tol = sqrt(std::numeric_limits<Real>::epsilon()), size_t max_refinements = 9);

    template<class F>
    Real integrate(const F f, Real* error = nullptr) const;

private:
    std::shared_ptr<detail::sinh_sinh_detail<Real>> m_imp;
};

template<class Real>
sinh_sinh<Real>::sinh_sinh(Real tol, size_t max_refinements) : m_imp(std::make_shared<detail::sinh_sinh_detail<Real>>(tol, max_refinements))
{
    return;
}

template<class Real>
template<class F>
Real sinh_sinh<Real>::integrate(const F f, Real* error) const
{
    using std::abs;
    using boost::math::detail::sinh_sinh_detail;
    if(abs(f(std::numeric_limits<Real>::max())) > std::numeric_limits<Real>::epsilon() || abs(f(std::numeric_limits<Real>::lowest())) > std::numeric_limits<Real>::epsilon())
    {
        throw std::domain_error("The function you are trying to integrate does not go to zero at infinity.\n");
    }
    return m_imp->integrate(f, error);
}

}}
#endif
