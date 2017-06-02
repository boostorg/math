//  (C) Copyright Nick Thompson 2017.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_CHEBYSHEV_HPP
#define BOOST_MATH_SPECIAL_CHEBYSHEV_HPP
#include <cmath>

namespace boost { namespace math {

template<class Real>
inline Real chebyshev_next(Real const & x, Real const & Tn, Real const & Tn_1)
{
    return 2*x*Tn - Tn_1;
}

namespace detail {

template<class Real, bool second=false>
inline Real chebyshev_imp(unsigned n, Real const & x)
{
    using std::cosh;
    using std::acosh;
    using std::pow;
    Real T0 = 1;
    Real T1;
    if (second)
    {
        if (x > 1 || x < -1)
        {
            Real t = sqrt(x*x -1);
            return (pow(x+t, n+1) - pow(x-t, n+1))/(2*t);
        }
        T1 = 2*x;
    }
    else
    {
        if (x > 1)
        {
            return cosh(n*acosh(x));
        }
        if (x < -1)
        {
            if (n & 1)
            {
                return -cosh(n*acosh(-x));
            }
            else
            {
                return cosh(n*acosh(-x));
            }
        }
        T1 = x;
    }

    if (n == 0)
    {
        return T0;
    }

    unsigned l = 1;
    while(l < n)
    {
       std::swap(T0, T1);
       T1 = boost::math::chebyshev_next(x, T0, T1);
       ++l;
    }
    return T1;
}
} // namespace detail

template<class Real>
Real chebyshev_t(unsigned n, Real const & x)
{
    return detail::chebyshev_imp<Real, false>(n, x);
}

template<class Real>
Real chebyshev_u(unsigned n, Real const & x)
{
    return detail::chebyshev_imp<Real, true>(n, x);
}

template<class Real>
Real chebyshev_t_prime(unsigned n, Real const & x)
{
    if (n == 0)
    {
        return 0;
    }
    return n*detail::chebyshev_imp<Real, true>(n - 1, x);
}

}}
#endif
