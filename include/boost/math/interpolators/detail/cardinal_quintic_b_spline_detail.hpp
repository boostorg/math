// Copyright Nick Thompson, 2019
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_INTERPOLATORS_CARDINAL_QUINTIC_B_SPLINE_DETAIL_HPP
#define BOOST_MATH_INTERPOLATORS_CARDINAL_QUINTIC_B_SPLINE_DETAIL_HPP
#include <vector>
#include <utility>

namespace boost{ namespace math{ namespace interpolators{ namespace detail{


template <class Real>
class cardinal_quintic_b_spline_detail
{
public:
    // If you don't know the value of the derivative at the endpoints, leave them as nans and the routine will estimate them.
    // y[0] = y(a), y[n -1] = y(b), step_size = (b - a)/(n -1).

    cardinal_quintic_b_spline_detail(const Real* const y,
                                size_t n,
                                Real t0 /* initial time, left endpoint */,
                                Real h  /*spacing, stepsize*/,
                                std::pair<Real, Real> left_endpoint_derivatives,
                                std::pair<Real, Real> right_endpoint_derivatives)
    {
        if (h <= 0) {
            throw std::logic_error("Spacing must be > 0.");
        }
        m_inv_h = 1/h;
        m_t0 = t0;

        if (n < 3) {
            throw std::logic_error("The interpolator requires at least 3 points.");
        }

        m_alpha.resize(n + 4);
        std::vector<Real> rhs(n+4);
        rhs[0] = 6*h*h*left_endpoint_derivatives.second;
        rhs[1] = 24*h*left_endpoint_derivatives.first;
        for (size_t i = 2; i < n + 2; ++i) {
            rhs[i] = 120*y[i-2];
        }
        rhs[n+2] = 24*h*right_endpoint_derivatives.first;
        rhs[n+3] = 6*h*h*right_endpoint_derivatives.second;

        std::vector<Real> diagonal(n+4);

    }

    Real operator()(Real t) const {
        if (t < m_t0 || t > m_t0 + (m_alpha.size()-2)/m_inv_h) {
            const char* err_msg = "Tried to evaluate the cardinal quadratic b-spline outside the domain of of interpolation; extrapolation does not work.";
            throw std::domain_error(err_msg);
        }
        return std::numeric_limits<Real>::quiet_NaN();
    }

    Real prime(Real t) const {
        if (t < m_t0 || t > m_t0 + (m_alpha.size()-2)/m_inv_h) {
            const char* err_msg = "Tried to evaluate the cardinal quadratic b-spline outside the domain of of interpolation; extrapolation does not work.";
            throw std::domain_error(err_msg);
        }
        return std::numeric_limits<Real>::quiet_NaN();
    }

    Real double_prime(Real t) const {
        if (t < m_t0 || t > m_t0 + (m_alpha.size()-2)/m_inv_h) {
            const char* err_msg = "Tried to evaluate the cardinal quadratic b-spline outside the domain of of interpolation; extrapolation does not work.";
            throw std::domain_error(err_msg);
        }
        return std::numeric_limits<Real>::quiet_NaN();
    }


    Real t_max() const {
        return m_t0 + (m_alpha.size()-3)/m_inv_h;
    }

private:
    std::vector<Real> m_alpha;
    Real m_inv_h;
    Real m_t0;
};

}}}}
#endif
