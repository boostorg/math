//  (C) Copyright Nick Thompson 2019.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_CARDINAL_B_SPLINE_HPP
#define BOOST_MATH_SPECIAL_CARDINAL_B_SPLINE_HPP

#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/binomial.hpp>

namespace boost { namespace math {

namespace detail {

template<class Real>
Real cardinal_b_spline_summation_impl(unsigned n, Real x)
{
    using std::floor;
    using std::pow;
    if (x < 0) {
        return cardinal_b_spline_summation_impl(n, -x);
    }
    Real z = x + Real(n+1)/2;
    unsigned kmax = std::min(static_cast<unsigned>(floor(z)), n+1);

    Real result = 0;
    for (unsigned k = 0; k <= kmax; ++k) {
        Real term = pow(z,n)*boost::math::binomial_coefficient<Real>(n+1,k);
        if (k&1) {
            result -= term;
        } else {
            result += term;
        }
        z -= 1;
    }
    return result/factorial<Real>(n);
}

template<class Real>
Real cardinal_b_spline_recursive_impl(unsigned n, Real x) {
    if (x < 0)
    {
        return detail::cardinal_b_spline_recursive_impl(n, -x);
    }
    if (n == 0) {
        if (x < Real(1)/Real(2)) {
            return Real(1);
        }
        else if (x > Real(1)/Real(2)) {
            return Real(0);
        }
        else if (x == Real(1)/Real(2)) {
            // Value of the Fourier series at the discontinuity:
            return Real(1)/Real(2);
        }
        else {
            // Who knows what happened?
            return std::numeric_limits<Real>::quiet_NaN();
        }
    }

    Real n_ = n;
    Real supp_max = (n_ + 1)/Real(2);
    if (x >= supp_max) {
        return Real(0);
    }
    Real a = (supp_max + x)/n_;
    Real b = (supp_max - x)/n_;
    return a*cardinal_b_spline_recursive_impl(n-1, x+1/Real(2)) + b*cardinal_b_spline_recursive_impl(n-1, x - 1/Real(2));
}
}

template<class Real>
Real cardinal_b_spline(unsigned n, Real x) {
    if (n <= 8) {
        return detail::cardinal_b_spline_recursive_impl(n, x);
    }
    return detail::cardinal_b_spline_summation_impl(n,x);
}


template<class Real>
Real cardinal_b_spline_prime(unsigned n, Real x) {
    if (x < 0) {
        return -cardinal_b_spline_prime(n, -x);
    }

    if (n==0) {
        // What is the sensible thing to do for the weak derivative?
        // Any better ideas are welcome.
        if (x == Real(1)/Real(2)) {
            return std::numeric_limits<Real>::infinity();
        }
        return Real(0);
    }
    return cardinal_b_spline(n-1, x+1/Real(2)) - cardinal_b_spline(n-1, x-1/Real(2));
}

template<class Real>
Real cardinal_b_spline_double_prime(unsigned n, Real x) {
    if (x < 0) {
        return cardinal_b_spline_double_prime(n, -x);
    }

    if (n==0) {
        throw std::domain_error("The second derivative of a discontinuous function does not exist.");
    }
    if (n==1) {
        if (x==1) {
            return std::numeric_limits<Real>::infinity();
        }
        return Real(0);
    }
    return cardinal_b_spline(n-2, x+1) - 2*cardinal_b_spline(n-2, x) + cardinal_b_spline(n-2, x-1);
}


template<class Real>
Real forward_cardinal_b_spline(unsigned n, Real x) {
    return cardinal_b_spline(n, x - (n+1)/Real(2));
}

}}
#endif
