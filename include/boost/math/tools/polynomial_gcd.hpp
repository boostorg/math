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
};




//
// Special handling for polynomials:
// Note that gcd_traits_polynomial is templated on the coefficient type T,
// not on polynomial<T>.
//
template <typename T, typename Enabled = void>
struct gcd_traits_polynomial 
{};


// gcd_traits for Z[x].
template <typename T>
struct gcd_traits_polynomial<T, typename enable_if_c< std::numeric_limits<T>::is_integer >::type > : public gcd_traits_polynomial_defaults<T>
{
    typedef boost::math::tools::polynomial<T> polynomial_type;
    using gcd_traits_polynomial_defaults<T>::normalize;
    
    inline static void
    subtract(polynomial_type &a, polynomial_type const &b)
    {
        // We want to use Stepanov's implementation as often as possible because
        // it results in the smallest coefficients, however, whole numbers take
        // precedence.
        T m = constant_coefficient(a);
        gcd_traits<T>::modulo(m, constant_coefficient(b));
        if (m == 0)
        {
            T const r = constant_coefficient(a) / constant_coefficient(b);
            // Stepanov's implementation; suffers from floating point inaccuracy
            // when r does not divide ___ evenly.
            a -= r * b;
        }
        else
        {
            // Antoine Joux's implementation: produces huge coefficients.
            T const tmp = constant_coefficient(a);
            a *= constant_coefficient(b);
            a -= tmp * b;
            normalize(a);
        }
    }
};


// TODO: Enable for multiprecision floating point types.
template <typename T>
struct gcd_traits_polynomial<T, typename disable_if_c< std::numeric_limits<T>::is_integer >::type > : public gcd_traits_polynomial_defaults<T>
{
    typedef boost::math::tools::polynomial<T> polynomial_type;
    typedef gcd_traits_defaults< polynomial_type > parent_class;
    using typename parent_class::method_type;
    // static const method_type method = method_type::method_euclid;

    inline static void
    subtract(polynomial_type &a, polynomial_type const &b)
    {
        using std::modf;
        
        // We want to use Stepanov's implementation as often as possible because
        // it results in the smallest coefficients, however, whole numbers take
        // precedence.
        T const r = constant_coefficient(a) / constant_coefficient(b);
        T r_int;
        T const r_frac = modf(r, &r_int);
        // Should we also consider 0.25, 0.125, etc?
        if (r_frac == T(0) || r_frac == T(0.5))
        {
            // Stepanov's implementation; suffers from floating point inaccuracy
            // when r does not divide ___ evenly.
            a -= r * b;
            // normalize coefficients so that leading coefficient is whole
            if (a && modf(a.data().back(), &r_int) != T(0))
                a /= a.data().back();
        }
        else
        {
            // Antoine Joux's implementation: produces huge coefficients.
            T const tmp = constant_coefficient(a);
            a *= constant_coefficient(b);
            a -= tmp * b;
#ifndef NDEBUG
            using boost::lambda::_1;
            if (std::numeric_limits<T>::has_infinity)
            {
                T const &inf(std::numeric_limits<T>::infinity());
                if (std::find_if(a.data().begin(), a.data().end(), _1 == inf || _1 == -inf) != a.data().end())
                    throw std::domain_error("Floating point overflow.");
            }
#endif
        }
    }
};


template <typename T>
struct gcd_traits< boost::math::tools::polynomial<T> > : public gcd_traits_polynomial<T>
{
    
};

} // namespace math
} // namespace boost

#endif // BOOST_MATH_TOOLS_POLYNOMIAL_GCD_HPP
