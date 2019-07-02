//  (C) Copyright Nick Thompson 2019.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_CARDINAL_B_SPLINE_HPP
#define BOOST_MATH_SPECIAL_CARDINAL_B_SPLINE_HPP

#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/binomial.hpp>

namespace boost { namespace math {

namespace detail {

template<class Real>
Real cardinal_b_spline_summation_impl(unsigned n, Real x)
{
    using std::floor;
    using std::pow;
    using std::ceil;
    if (x < 0) {
        return cardinal_b_spline_summation_impl(n, -x);
    }
    if (x >= Real(n+1)/2) {
        return Real(0);
    }
    if (x <= Real(n+1)/4) {
        Real z = x + Real(n+1)/2;
        int kmax = floor(z);

        Real result = pow(z,n);
        z -= 1;
        Real kfact = 1;
        Real num = n+1;
        for (int k = 1; k <= kmax; ++k)
        {
            kfact *= k;
            Real term = pow(z,n)*num/kfact;
            num *= n-k+1;
            if (k&1) {
                result -= term;
            } else {
                result += term;
            }
            z -= 1;
        }
        return result/factorial<Real>(n);
    }
    else {
        // Schoenberg, Cardinal Spline interpolation, equation 1.4:
        using std::ceil;
        int kmin = ceil(x + Real(n+1)/2);
        Real z = kmin - (x + Real(n+1)/2);
        Real result = 0;
        for (int k = kmin; k <= n+1; ++k)
        {
            Real term = binomial_coefficient<Real>(n+1, k)*pow(z, n);
            if (k&1)
            {
                result -= term;
            } else {
                result += term;
            }
            z += 1;
        }
        if (n & 1)
        {
            return result/factorial<Real>(n);
        }
        else
        {
            return -result/factorial<Real>(n);
        }

    }
}

template<class Real>
Real cardinal_b_spline_recursive_impl(unsigned n, Real x)
{
    if (x < 0)
    {
        return detail::cardinal_b_spline_recursive_impl(n, -x);
    }
    if (n == 0)
    {
        if (x < Real(1)/Real(2))
        {
            return Real(1);
        }
        else if (x > Real(1)/Real(2))
        {
            return Real(0);
        }
        else if (x == Real(1)/Real(2))
        {
            // Value of the Fourier series at the discontinuity:
            return Real(1)/Real(2);
        }
        else
        {
            // Who knows what happened?
            return std::numeric_limits<Real>::quiet_NaN();
        }
    }

    Real n_ = n;
    Real supp_max = (n_ + 1)/Real(2);
    if (x >= supp_max)
    {
        return Real(0);
    }
    Real a = (supp_max + x)/n_;
    Real b = (supp_max - x)/n_;
    return a*cardinal_b_spline_recursive_impl(n-1, x+1/Real(2)) + b*cardinal_b_spline_recursive_impl(n-1, x - 1/Real(2));
}
}

template<class Real>
Real cardinal_b_spline(unsigned n, Real x)
{
    if (n <= 5)
    {
        return detail::cardinal_b_spline_recursive_impl(n, x);
    }
    return detail::cardinal_b_spline_summation_impl(n,x);
}


template<class Real>
Real cardinal_b_spline_prime(unsigned n, Real x)
{
    if (x < 0)
    {
        return -cardinal_b_spline_prime(n, -x);
    }

    if (n==0)
    {
        // What is the sensible thing to do for the weak derivative?
        // Any better ideas are welcome.
        if (x == Real(1)/Real(2))
        {
            return std::numeric_limits<Real>::infinity();
        }
        return Real(0);
    }
    return cardinal_b_spline(n-1, x+1/Real(2)) - cardinal_b_spline(n-1, x-1/Real(2));
}

template<class Real>
Real cardinal_b_spline_double_prime(unsigned n, Real x) {
    if (x < 0)
    {
        return cardinal_b_spline_double_prime(n, -x);
    }

    if (n==0)
    {
        throw std::domain_error("The second derivative of a discontinuous function does not exist.");
    }
    if (n==1)
    {
        if (x==1)
        {
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

// Fractional Splines and Wavelets, Unser and Blu, Siam Review, 2000:
// https://infoscience.epfl.ch/record/63073/files/unser9901.pdf, equation 2.2:
template<class Real>
Real forward_fractional_cardinal_b_spline(Real alpha, Real x) {
    using std::pow;
    using std::floor;

    if (alpha <= -1)
    {
        throw std::domain_error("Fractional B-splines are not defined for alpha < -1.");
    }
    if (x <= 0)
    {
        return Real(0);
    }
    using std::floor;
    int kmax = static_cast<int>(floor(x));
    Real kfact = 1;
    Real numerator = 1;
    Real result = pow(x, alpha);
    for (int k = 1; k <= kmax; ++k)
    {
        x -= 1;
        kfact *= -k;
        numerator *=  (alpha + 2 - k);
        result += numerator*pow(x, alpha)/kfact;
    }
    return result/boost::math::tgamma(alpha+1);
}

}}
#endif
