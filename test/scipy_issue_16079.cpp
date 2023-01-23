// Copyright Matt Borland, 2022
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0. (See accompanying fi
//
// See: https://github.com/scipy/scipy/issues/16079

#include <boost/math/distributions/hypergeometric.hpp>

using boost::math::cdf;
using boost::math::hypergeometric;

int main(void) 
{
    for (int N = 48000; N < 50000; N++) 
    {
        for (int r = 55; r < 59; r++) 
        {
            hypergeometric hg(r, 35775, N);
            cdf(hg, 0);
        }
    }

    return 0;
}
