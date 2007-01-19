//  Copyright (c) 2006 Xiaogang Zhang
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_BESSEL_KN_HPP
#define BOOST_MATH_BESSEL_KN_HPP

#include <boost/math/special_functions/detail/bessel_k0.hpp>
#include <boost/math/special_functions/detail/bessel_k1.hpp>
#include <boost/math/tools/error_handling.hpp>

// Modified Bessel function of the second kind of integer order
// K_n(z) is the dominant solution, forward recurrence always OK (though unstable)

namespace boost { namespace math { namespace detail{

template <typename T>
T bessel_kn(int n, T x)
{
    T value, current, prev;

    using namespace boost::math::tools;

    if (x < 0)
    {
        return domain_error<T>(BOOST_CURRENT_FUNCTION,
            "Got x = %1%, but argument x must be non-negative, complex number result not supported.", x);
    }
    if (x == 0)
    {
        return overflow_error<T>(BOOST_CURRENT_FUNCTION);
    }

    if (n == 0)
    {
        return bessel_k0(x);
    }
    if (n == 1)
    {
        return bessel_k1(x);
    }
    if (n < 0)
    {
        n = -n;                             // K_{-n}(z) = K_n(z)
    }

    prev = bessel_k0(x);
    current = bessel_k1(x);
    for (int k = 1; k < n; k++)            // n >= 2
    {
        value = 2 * k * current / x + prev;
        prev = current;
        current = value;
    }

    return value;
}

}}} // namespaces

#endif // BOOST_MATH_BESSEL_KN_HPP
