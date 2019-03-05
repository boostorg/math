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
#include <cmath>
#include <boost/multiprecision/float128.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/detail/daubechies_scaling_integer_grid.hpp>
#include <boost/math/special_functions/daubechies_filters.hpp>


namespace boost::math {

namespace detail {

template<class Real, int p, int order>
std::vector<Real> dyadic_grid(size_t j_max)
{
    auto c = boost::math::daubechies_scaling_filter<Real, p>();
    Real scale = boost::math::constants::root_two<Real>()*(1 << order);
    for (auto & x : c)
    {
        x *= scale;
    }

    auto phik = daubechies_scaling_integer_grid<Real, p, order>();

    // Maximum sensible j for 32 bit floats is j_max = 22:
    std::vector<Real> v(2*p + (2*p-1)*((1<<j_max) -1), std::numeric_limits<Real>::quiet_NaN());
    v.resize(2*p + (2*p-1)*((1<<j_max) -1), std::numeric_limits<Real>::quiet_NaN());
    v[0] = 0;
    v[v.size()-1] = 0;
    for (size_t i = 0; i < phik.size(); ++i) {
        v[i*(1<<(j_max))] = phik[i];
    }

    for (size_t j = 1; j <= j_max; ++j)
    {
        size_t k_max = v.size()/(1 << (j_max-j));
        for (size_t k = 1; k < k_max;  k += 2)
        {
            // Where this value will go:
            size_t delivery_idx = k*(1 << (j_max-j));
            if (delivery_idx >= v.size())
            {
                std::cout << "Delivery index out of range!\n";
                continue;
            }
            Real term = 0;
            for (size_t l = 0; l < c.size(); ++l) {
                size_t idx = k*(1 << (j_max - j + 1)) - l*(1 << j_max);
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
    daubechies_scaling(int levels = -1)
    {
        using boost::multiprecision::float128;
        if (levels < 0)
        {
            m_levels = 22;
            if constexpr (std::is_same_v<Real, float128>)
            {
                m_v = detail::dyadic_grid<Real, p, 0>(m_levels);
            }
            else
            {
                auto v = detail::dyadic_grid<float128, p, 0>(m_levels);
                m_v.resize(v.size());
                for (size_t i = 0; i < v.size(); ++i)
                {
                    m_v[i] = static_cast<Real>(v[i]);
                }
            }
        }
        else
        {
            m_levels = levels;
            m_v = detail::dyadic_grid<Real, p, 0>(m_levels);
        }

        //float128 inv_spacing = float128(m_v.size()-1)/float128(2*p-1);
        //m_inv_spacing = static_cast<Real>(inv_spacing);
        m_inv_spacing = (1 << m_levels);
    }

    Real operator()(Real x) const { return x; }

    std::pair<int, int> support() const {
        return {0, 2*p-1};
    }

    Real constant_interpolation(Real x) const {
        if (x <= 0 || x >= 2*p-1) {
            return Real(0);
        }

        using std::ceil;
        using std::floor;
        Real rescale = (1<<m_levels)*x;
        Real idx1 = floor((1<<m_levels)*x);
        Real idx2 = ceil((1<<m_levels)*x);

        size_t idx = 0;
        if (idx2 - rescale < rescale - idx1) {
            idx = static_cast<size_t>(idx2);
        } else {
            idx = static_cast<size_t>(idx1);
        }
        //std::cout << "Constant idx1 = " << idx1 << ", idx2 = " << idx2 << "\n";
        return m_v[idx];
    }

    Real linear_interpolation(Real x) const {
        if (x <= 0 || x >= 2*p-1) {
            return 0;
        }
        using std::ceil;
        using std::floor;

        Real y = (1<<m_levels)*x;
        Real idx1 = floor(y);
        Real idx2 = ceil(y);

        auto v1 = m_v[static_cast<size_t>(idx1)];
        auto v2 = m_v[static_cast<size_t>(idx2)];

        Real t = y - idx1;
        return (1-t)*v1 + t*v2;
    }

    size_t bytes() const
    {
        size_t s = sizeof(*this);
        s += m_v.size()*sizeof(Real);
        return s;
    }

    size_t size() const {
        return m_v.size();
    }

    Real index_to_abscissa(size_t i) {
        return i/m_inv_spacing;
    }

    // This is for debugging only; use the .at()
    Real operator[](size_t i) const {
        return m_v.at(i);
    }

    Real spacing() const {
        return 1/m_inv_spacing;
    }

private:
    size_t m_levels;
    Real m_inv_spacing;
    //std::array<Real, 2*p-1> m_c;
    std::vector<Real> m_v;
};

}
#endif
