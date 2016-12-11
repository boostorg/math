//  (C) Copyright John Maddock and Jeremy William Murphy 2016.

//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_POLYNOMIAL_GCD_HPP
#define BOOST_MATH_TOOLS_POLYNOMIAL_GCD_HPP

#ifdef _MSC_VER
#pragma once
#endif

#include <boost/math/tools/polynomial.hpp>
#include <boost/math/common_factor_rt.hpp>

namespace boost { namespace math {

// Common gcd traits for polynomials.
template <typename T>
struct gcd_traits_polynomial_defaults : public gcd_traits_defaults< boost::math::tools::polynomial<T> >
{
    typedef boost::math::tools::polynomial<T> polynomial_type;

    static polynomial_type
    abs(const polynomial_type &val)
    {
        return leading_coefficient(val) < T(0) ? -val : val;
    }
    
    inline static int make_odd(polynomial_type &x)
    {
        unsigned r = 0;
        while (even(x))
        {
            x >>= 1;
            r++;
        }
        return r;
    }
    
    inline static bool
    less(polynomial_type const &a, polynomial_type const &b)
    {
        return a.size() < b.size();
    }

    inline static 
    void normalize(polynomial_type &x)
    {
        using boost::lambda::_1;
        
        // This does assume that the coefficients are totally ordered.
        if (x)
        {
            if (leading_coefficient(x) < T(0))
                x.negate();
            // Find the first non-zero, because we can't do gcd(0, 0).
            T const d = gcd_range(find_if(x.data().begin(), x.data().end(), _1 != T(0)), x.data().end()).first;
            x /= d;
        }
    }

    inline static void
    subtract(polynomial_type &a, polynomial_type const &b)
    {
        // We use Stepanov's implementation when the constant coefficients
        // divide evenly, and Joux's otherwise.
        T m = constant_coefficient(a);
        gcd_traits<T>::modulo(m, constant_coefficient(b));
        if (!m)
        {
            T const r = constant_coefficient(a) / constant_coefficient(b);
            a -= r * b;
        }
        else
        {
            // Antoine Joux's implementation tempered by coefficient gcd.
            T const a0 = constant_coefficient(a);
            T const b0 = constant_coefficient(b);
            T const gcd_a0b0 = gcd(a0, b0);
            a *= b0 / gcd_a0b0;
            a -= (a0 / gcd_a0b0) * b;
        }
    }
};


//
// Special handling for polynomials:
// Note that gcd_traits_polynomial is templated on the coefficient type T,
// not on polynomial<T>.
//
template <typename T, typename Enabled = void>
struct gcd_traits_polynomial 
{};


// Note: Use these trait classes for Z[x] and R[x], etc to customize for operations
// on different kinds of polynomials.

// gcd_traits for Z[x].
template <typename T>
struct gcd_traits_polynomial<T, typename enable_if_c<std::numeric_limits<T>::is_integer>::type> : public gcd_traits_polynomial_defaults<T>
{};


// gcd_traits for R[x].
template <typename T>
struct gcd_traits_polynomial<T, typename enable_if_c<is_floating_point<T>::value || (std::numeric_limits<T>::is_specialized && !std::numeric_limits<T>::is_integer) >::type> : public gcd_traits_polynomial_defaults<T>
{};


template <typename T>
struct gcd_traits< boost::math::tools::polynomial<T> > : public gcd_traits_polynomial<T>
{};

} // namespace math
} // namespace boost

#endif // BOOST_MATH_TOOLS_POLYNOMIAL_GCD_HPP
