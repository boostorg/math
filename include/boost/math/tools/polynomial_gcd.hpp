//  (C) Copyright Jeremy William Murphy 2016.

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


namespace boost{ namespace math{ namespace tools{
    
/* From Knuth, 4.6.1:
* 
* We may write any nonzero polynomial u(x) from R[x] where R is a UFD as
*
*      u(x) = cont(u) · pp(u(x))
*
* where cont(u), the content of u, is an element of S, and pp(u(x)), the primitive
* part of u(x), is a primitive polynomial over S. 
* When u(x) = 0, it is convenient to define cont(u) = pp(u(x)) = O.
*/

template <class T>
T content(polynomial<T> const &x)
{
    return x ? gcd_range(x.data().begin(), x.data().end()).first : T(0);
}

// Knuth, 4.6.1
template <class T>
polynomial<T> primitive_part(polynomial<T> const &x, T const &cont)
{
    return x ? x / cont : polynomial<T>();
}


template <class T>
polynomial<T> primitive_part(polynomial<T> const &x)
{
    return primitive_part(x, content(x));
}


// Trivial but useful convenience function referred to simply as l() in Knuth.
template <class T>
T leading_coefficient(polynomial<T> const &x)
{
    return x.data().back();
}

namespace detail
{
    /* Reduce u and v to their primitive parts and return the gcd of their 
    * contents. Used in a couple of gcd algorithms.
    */
    template <class T>
    T reduce_to_primitive(polynomial<T> &u, polynomial<T> &v)
    {
        T const u_cont = content(u), v_cont = content(v);
        u /= u_cont;
        v /= v_cont;
        return gcd(u_cont, v_cont);
    }
}


/**
* Knuth, The Art of Computer Programming: Volume 2, Third edition, 1998
* Algorithm 4.6.1C: Greatest common divisor over a unique factorization domain.
* 
* The subresultant algorithm by George E. Collins [JACM 14 (1967), 128-142], 
* later improved by W. S. Brown and J. F. Traub [JACM 18 (1971), 505-514].
* 
* Although step C3 keeps the coefficients to a "reasonable" size, they are
* still potentially several binary orders of magnitude larger than the inputs.
* Thus, this algorithm should only be used where T is a multi-precision type.
* 
* @tparam  T   Polynomial coefficient type.
* @return      Greatest common divisor of polynomials u and v.
*/
template <class T>
typename enable_if_c< std::numeric_limits<T>::is_integer, polynomial<T> >::type
subresultant_gcd(polynomial<T> u, polynomial<T> v)
{
    BOOST_ASSERT(u || v);
    
    if (!u)
        return v;
    if (!v)
        return u;
    
    typedef typename polynomial<T>::size_type N;
    
    T const d = detail::reduce_to_primitive(u, v);
    T g = 1, h = 1;
    polynomial<T> r;
    while (true)
    {
        // Pseudo-division.
        r = u % v;
        if (!r)
            return d * primitive_part(v); // Attach the content.
        if (r.degree() == 0)
        {
            v = polynomial<T>(T(1));
            return d * primitive_part(v); // Attach the content.
        }
        N const delta = u.degree() - v.degree();
        // Adjust remainder.
        u = v;
        v = r / (g * detail::integer_power(h, delta));
        g = leading_coefficient(u);
        T const tmp = detail::integer_power(g, delta);
        if (delta <= N(1))
            h = tmp * detail::integer_power(h, N(1) - delta);
        else
            h = tmp / detail::integer_power(h, delta - N(1));
    }
}


template <class T>
typename enable_if_c< std::numeric_limits<T>::is_integer, polynomial<T> >::type
gcd(polynomial<T> const &u, polynomial<T> const &v)
{
    return subresultant_gcd(u, v);
}
    
}}} // namespace boost::math::tools

#endif
