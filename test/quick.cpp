// Copyright 2026 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/math/special_functions/gamma.hpp>
#include <cstdlib>

int main()
{
    auto result = boost::math::tgamma(2.0);
    // tgamma(2) == 1! == 1.0
    if (result < 0.99 || result > 1.01)
    {
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
