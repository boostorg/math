/*
 *  (C) Copyright Nick Thompson 2018.
 *  Use, modification and distribution are subject to the
 *  Boost Software License, Version 1.0. (See accompanying file
 *  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 *  The Halton sequence for quasi-random number generation.
 */
#ifndef BOOST_MATH_TOOLS_VAN_DER_CORPUT_HPP
#define BOOST_MATH_TOOLS_VAN_DER_CORPUT_HPP
#include <boost/math/policies/error_handling.hpp>

namespace boost { namespace math {

template<class Real, class Z, class Policy = policies::policy<>>
Real van_der_corput(Z x, Z base)
{
    static_assert(std::numeric_limits<Z>::is_integer,
                  "The van der Corput function takes integers as arguments and returns real values.\n");
    if (base < 2 || x < 0)
    {
        static const char* function = "boost::math::van_der_corput<%1%>";
        boost::math::policies::raise_domain_error<Real>(
              function, "Base must be >= 2 and argument >= 0, but got base %1%", base, Policy());
    }
    Real r = 0;
    Real v = 1;
    Real inv_b = static_cast<Real>(1) / static_cast<Real>(base);

    while (x > 0)
    {
        v *= inv_b;
        r += v*static_cast<Real>(x % base);
        x /= base;
    }
    return r;
}

}}
#endif
