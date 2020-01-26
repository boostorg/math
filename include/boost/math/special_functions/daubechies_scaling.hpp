/*
 * Copyright Nick Thompson, 2020
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef BOOST_MATH_SPECIAL_DAUBECHIES_SCALING_HPP
#define BOOST_MATH_SPECIAL_DAUBECHIES_SCALING_HPP
#include <vector>
#include <array>
#include <cmath>
#include <thread>
#include <boost/multiprecision/float128.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/detail/daubechies_scaling_integer_grid.hpp>
#include <boost/math/filters/daubechies.hpp>


namespace boost::math {

namespace detail {

template<class Real, int p, int order>
std::vector<Real> dyadic_grid(size_t j_max)
{
    auto c = boost::math::filters::daubechies_scaling_filter<Real, p>();
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
                std::cerr << "Delivery index out of range!\n";
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
                std::cerr << "Delivery index already populated!, = " << v[delivery_idx] << "\n";
                std::cerr << "would overwrite with " << term << "\n";
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
        }
        else {
            m_levels = levels;
        }

        auto f1 = [this] {
            auto v = detail::dyadic_grid<float128, p, 0>(this->m_levels);
            this->m_v.resize(v.size());
            for (size_t i = 0; i < v.size(); ++i) {
                this->m_v[i] = static_cast<Real>(v[i]);
            }
        };

        auto f2 = [this] {
            auto v_prime = detail::dyadic_grid<float128, p, 1>(this->m_levels);
            this->m_v_prime.resize(v_prime.size());
            for (size_t i = 0; i < v_prime.size(); ++i) {
                this->m_v_prime[i] = static_cast<Real>(v_prime[i]);
            }
        };


        auto f3 = [this] {
            auto v_dbl_prime = detail::dyadic_grid<float128, p, 2>(this->m_levels);
            this->m_v_dbl_prime.resize(v_dbl_prime.size());
            for (size_t i = 0; i < v_dbl_prime.size(); ++i) {
                this->m_v_dbl_prime[i] = static_cast<Real>(v_dbl_prime[i]);
            }
        };

        std::thread t1(f1);
        std::thread t2(f2);
        std::thread t3(f3);

        t1.join();
        t2.join();
        t3.join();

        m_inv_spacing = (1 << m_levels);

        m_c = boost::math::filters::daubechies_scaling_filter<Real, p>();
        Real scale = boost::math::constants::root_two<Real>();
        for (auto & x : m_c)
        {
            x *= scale;
        }
    }

    Real operator()(Real x) const { return this->linear_interpolation(x); }

    std::pair<int, int> support() const {
        return {0, 2*p-1};
    }

    Real constant_interpolation(Real x) const {
        if (x <= 0 || x >= 2*p-1) {
            return Real(0);
        }
        using std::floor;
        Real y = (1 << m_levels)*x;
        Real k = floor(y);

        if (y - k < k + 1 - y)
        {
            return m_v[static_cast<size_t>(k)];
        }
        return m_v[static_cast<size_t>(k)+1];
    }

    Real linear_interpolation(Real x) const {
        if (x <= 0 || x >= 2*p-1) {
            return Real(0);
        }
        using std::floor;

        Real y = (1<<m_levels)*x;
        Real k = floor(y);

        size_t kk = static_cast<size_t>(k);

        Real t = y - k;
        return (1-t)*m_v[kk] + t*m_v[kk+1];
    }

    Real single_crank_linear(Real x) const {
        if (x <= 0 || x >= 2*p-1) {
            return Real(0);
        }
        using std::floor;
        Real y = (1<<m_levels)*x;
        Real idx1 = floor(y);

        Real term = 0;
        size_t k = 2*idx1 + 1;
        for (size_t l = 0; l < m_c.size(); ++l) {
            uint64_t idx = k - l*(1 << m_levels);
            if (idx >= 0 && idx < m_v.size()) {
                term += m_c[l]*m_v[idx];
            }
        }

        if (y - idx1 < idx1 + 1 - y)
        {
            Real t = 2*(y - idx1);
            return (1-t)*m_v[static_cast<size_t>(idx1)] + t*term;
        }
        else {
            Real t = 2*(idx1 + 1 - y);
            return t*term + (1-t)*m_v[static_cast<size_t>(idx1)+1];
        }
    }

    Real single_crank_quadratic(Real x) const {
        if (x <= 0 || x >= 2*p-1) {
            return Real(0);
        }
        using std::floor;
        Real y = (1<<m_levels)*x;
        Real idx1 = floor(y);

        Real term = 0;
        size_t k = 2*idx1 + 1;
        for (size_t l = 0; l < m_c.size(); ++l) {
            uint64_t idx = k - l*(1 << m_levels);
            if (idx >= 0 && idx < m_v.size()) {
                term += m_c[l]*m_v[idx];
            }
        }

        Real y0 = m_v[static_cast<size_t>(idx1)];
        Real y1 = term;
        Real y2 = m_v[static_cast<size_t>(idx1)+1];

        Real a = (y2+y0-2*y1)/2;
        Real b = (4*y1-3*y0 - y2)/2;
        Real t = 2*(y - idx1);
        return a*t*t + b*t + y0;
    }


    Real double_crank_linear(Real x) const {
        return std::numeric_limits<Real>::quiet_NaN();
    }


    Real first_order_taylor(Real x) const {
        if (x <= 0 || x >= 2*p-1) {
            return 0;
        }
        using std::floor;

        Real y = (1<<m_levels)*x;
        Real k = floor(y);

        size_t kk = static_cast<size_t>(k);
        if (y - k < k + 1 - y)
        {
            Real eps = (y-k)/(1<<m_levels);
            return m_v[kk] + eps*m_v_prime[kk];
        }
        else {
            Real eps = (y-k-1)/(1<<m_levels);
            return m_v[kk+1] + eps*m_v_prime[kk+1];
        }
    }

    Real second_order_taylor(Real x) const {
        if (x <= 0 || x >= 2*p-1) {
            return 0;
        }
        using std::floor;

        Real y = (1<<m_levels)*x;
        Real k = floor(y);

        size_t kk = static_cast<size_t>(k);
        if (y - k < k + 1 - y)
        {
            Real eps = (y-k)/(1<<m_levels);
            return m_v[kk] + eps*m_v_prime[kk] + eps*eps*m_v_dbl_prime[kk]/2;
        }
        else {
            Real eps = (y-k-1)/(1<<m_levels);
            return m_v[kk+1] + eps*m_v_prime[kk+1] + eps*eps*m_v_dbl_prime[kk+1]/2;
        }
    }

    Real third_order_taylor(Real x) const {
        return std::numeric_limits<Real>::quiet_NaN();
    }

    Real single_crank_first_order_taylor(Real x) const {
        return std::numeric_limits<Real>::quiet_NaN();
    }

    Real double_crank_first_order_taylor(Real x) const {
        return std::numeric_limits<Real>::quiet_NaN();
    }

    Real single_crank_second_order_taylor(Real x) const {
        return std::numeric_limits<Real>::quiet_NaN();
    }

    Real double_crank_second_order_taylor(Real x) const {
        return std::numeric_limits<Real>::quiet_NaN();
    }

    Real single_crank_third_order_taylor(Real x) const {
        return std::numeric_limits<Real>::quiet_NaN();
    }

    Real double_crank_third_order_taylor(Real x) const {
        return std::numeric_limits<Real>::quiet_NaN();
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

    Real prime(size_t i) const {
        return m_v_prime.at(i);
    }

    Real spacing() const {
        return 1/m_inv_spacing;
    }

    auto begin() const {
        return m_v.begin();
    }

    auto end() const {
        return m_v.end();
    }

    auto data() const {
        return m_v.data();
    }

private:
    size_t m_levels;
    Real m_inv_spacing;
    std::array<Real, 2*p> m_c;
    std::vector<Real> m_v;
    std::vector<Real> m_v_prime;
    std::vector<Real> m_v_dbl_prime;
};

}
#endif
