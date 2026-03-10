// Copyright 2026 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/math/tr1.hpp>
#include <cstdlib>

int main()
{
    // boost_cbrt from boost_math_c99
    auto cbrt_result = boost::math::tr1::boost_cbrt(8.0);
    if (cbrt_result < 1.99 || cbrt_result > 2.01)
    {
        return EXIT_FAILURE;
    }

    // boost_riemann_zeta from boost_math_tr1
    // zeta(2) == pi^2/6 ~= 1.6449
    auto zeta_result = boost::math::tr1::boost_riemann_zeta(2.0);
    if (zeta_result < 1.64 || zeta_result > 1.65)
    {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
