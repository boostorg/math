/*
 * Copyright Nick Thompson, 2018
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef BOOST_MATH_SPECIAL_DAUBECHIES_SCALING_HPP
#define BOOST_MATH_SPECIAL_DAUBECHIES_SCALING_HPP
#include <vector>
#include <array>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/detail/daubechies_scaling_integer_grid.hpp>
#include <boost/math/special_functions/daubechies_filters.hpp>


namespace boost::math {

namespace detail {

template<class Real, int p>
std::vector<Real> dyadic_grid(int j_max)
{
    std::array<Real, 2*p> c = boost::math::daubechies_scaling_filter<Real, p>();
    for (auto & x : c)
    {
        x *= boost::math::constants::root_two<Real>();
    }

    auto phik = daubechies_scaling_integer_grid<Real, p, 0>();

    // Maximum sensible j for 32 bit floats is j_max = 22:
    std::vector<Real> v(2*p + (2*p-1)*( (1<<j_max) -1), std::numeric_limits<Real>::quiet_NaN());
    v.resize(2*p + (2*p-1)*( (1<<j_max) -1), std::numeric_limits<Real>::quiet_NaN());
    v[0] = 0;
    v[v.size()-1] = 0;
    for (size_t i = 0; i < 2*p - 2; ++i) {
        v[(i+1)*(1<<(j_max))] = phik[i];
    }

    for (int j = 1; j <= j_max; ++j)
    {
        int k_max = v.size()/(1 << (j_max-j));
        for (int k = 1; k < k_max;  k += 2)
        {
            // Where this value will go:
            int delivery_idx = k*(1 << (j_max-j));
            if (delivery_idx >= v.size())
            {
                std::cout << "Delivery index out of range!\n";
                continue;
            }
            Real term = 0;
            for (int l = 0; l < c.size(); ++l) {
                int idx = k*(1 << (j_max - j + 1)) - l*(1 << j_max);
                if (idx >= 0 && idx < v.size()) {
                    term += c[l]*v[idx];
                }
            }
            if (!isnan(v[delivery_idx])) {
                std::cout << "Delivery index already populated!, = " << v[delivery_idx] << "\n";
                std::cout << "would overwrite with " << term << "\n";
            }
            v[delivery_idx] = term;
        }
    }
    return v;
}

}

template<class Real, int p>
class daubechies_scaling {
public:
    daubechies_scaling(int levels = -1) {
        std::vector<Real> v;
        if (levels <= 0) {
            v = detail::dyadic_grid<Real, p>(20);
        }
        else {
            v = detail::dyadic_grid<Real, p>(levels);
        }

    }

    Real operator()(Real x) const { return x; }
};

}
#endif
