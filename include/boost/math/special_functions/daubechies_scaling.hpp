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
#include <future>
#include <iostream>
#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
#endif
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/detail/daubechies_scaling_integer_grid.hpp>
#include <boost/math/filters/daubechies.hpp>
#include <boost/math/interpolators/detail/cubic_hermite_detail.hpp>
#include <boost/math/interpolators/detail/quintic_hermite_detail.hpp>
#include <boost/math/interpolators/detail/septic_hermite_detail.hpp>


namespace boost::math {

namespace detail {

template<class Real, int p, int order>
std::vector<Real> dyadic_grid(int64_t j_max)
{
    using std::isnan;
    auto c = boost::math::filters::daubechies_scaling_filter<Real, p>();
    Real scale = boost::math::constants::root_two<Real>()*(1 << order);
    for (auto & x : c)
    {
        x *= scale;
    }

    auto phik = daubechies_scaling_integer_grid<Real, p, order>();

    // Maximum sensible j for 32 bit floats is j_max = 22:
    std::vector<Real> v(2*p + (2*p-1)*((1<<j_max) -1), std::numeric_limits<Real>::quiet_NaN());
    v[0] = 0;
    v[v.size()-1] = 0;
    for (int64_t i = 0; i < (int64_t) phik.size(); ++i) {
        v[i*(1<<j_max)] = phik[i];
    }

    for (int64_t j = 1; j <= j_max; ++j)
    {
        int64_t k_max = v.size()/(1 << (j_max-j));
        for (int64_t k = 1; k < k_max;  k += 2)
        {
            // Where this value will go:
            int64_t delivery_idx = k*(1 << (j_max-j));
            // This is a nice check, but we've tested this exhaustively, and it's expensive:
            //if (delivery_idx >= (int64_t) v.size()) {
            //    std::cerr << "Delivery index out of range!\n";
            //    continue;
            //}
            Real term = 0;
            for (int64_t l = 0; l < (int64_t) c.size(); ++l)
            {
                int64_t idx = k*(1 << (j_max - j + 1)) - l*(1 << j_max);
                if (idx < 0)
                {
                    break;
                }
                if (idx < (int64_t) v.size())
                {
                    term += c[l]*v[idx];
                }
            }
            // Again, another nice check:
            //if (!isnan(v[delivery_idx])) {
            //    std::cerr << "Delivery index already populated!, = " << v[delivery_idx] << "\n";
            //    std::cerr << "would overwrite with " << term << "\n";
            //}
            v[delivery_idx] = term;
        }
    }


    return v;
}

template<class RandomAccessContainer>
class matched_holder {
public:
    using Real = typename RandomAccessContainer::value_type;

    matched_holder(RandomAccessContainer && y, RandomAccessContainer && dydx, int grid_refinements) : y_{std::move(y)}, dy_{std::move(dydx)}
    {
        inv_h_ = (1 << grid_refinements);
        Real h = 1/inv_h_;
        for (auto & dy : dy_)
        {
            dy *= h;
        }
    }

    inline Real operator()(Real x) const {
        using std::floor;
        using std::sqrt;
        // This is the exact Holder exponent, but it's pessimistic almost everywhere!
        // It's only exactly right at dyadic rationals.
        //Real const alpha = 2 - log(1+sqrt(Real(3)))/log(Real(2));
        // We're gonna use alpha = 1/2, rather than 0.5500...
        Real s = x*inv_h_;
        Real ii = floor(s);
        auto i = static_cast<decltype(y_.size())>(ii);
        Real t = s - ii;
        Real dphi = dy_[i+1];
        Real diff = y_[i+1] - y_[i];
        return y_[i] + (2*dphi - diff)*t + 2*sqrt(t)*(diff-dphi);
    }

private:
    Real inv_h_;
    RandomAccessContainer y_;
    RandomAccessContainer dy_;
};


template<class RandomAccessContainer>
class linear_interpolation {
public:
    using Real = typename RandomAccessContainer::value_type;

    linear_interpolation(RandomAccessContainer && y,  int grid_refinements) : y_{std::move(y)}
    {
        s_ = (1 << grid_refinements);
    }

    inline Real operator()(Real x) const {
        using std::floor;
        Real y = x*s_;
        Real k = floor(y);

        int64_t kk = static_cast<int64_t>(k);
        Real t = y - k;
        return (1-t)*y_[kk] + t*y_[kk+1];
    }

private:
    Real s_;
    RandomAccessContainer y_;
};

}

template<class Real, int p>
class daubechies_scaling {
public:
    daubechies_scaling(int grid_refinements = -1)
    {
        static_assert(p <= 15, "Scaling functions only implements up to p = 15");
        #ifdef BOOST_HAS_FLOAT128
        using boost::multiprecision::float128;
        #endif
        if (grid_refinements < 0)
        {
            if (std::is_same_v<Real, float>)
            {
                //                          p= 2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15
                std::array<int, 16> r{-1, -1, 22, 21, 19, 17, 16, 15, 14, 13, 12, 11, 11, 11, 11, 11};
                grid_refinements = r[p];
            }
            else if (std::is_same_v<Real, double>)
            {
                //                          p= 2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15
                std::array<int, 16> r{-1, -1, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 20, 19, 18, 18};
                grid_refinements = r[p];
            }
            else
            {
                grid_refinements = 21;
            }
        }

        // Compute the refined grid:
        // In fact for float precision I know the grid must be computed in double precision and then cast back down, or else parts of the support are systematically inaccurate.
        std::future<std::vector<Real>> t0 = std::async(std::launch::async, [&grid_refinements]() {
            // Computing in higher precision and downcasting is essential for 1ULP evaluation in float precision:
            if constexpr (std::is_same_v<Real, float>) {
                auto v = detail::dyadic_grid<double, p, 0>(grid_refinements);
                std::vector<float> w(v.size());
                for (size_t i = 0; i < v.size(); ++i) {
                    w[i] = v[i];
                }
                return w;
            }
            else if constexpr (std::is_same_v<Real, double>) {
                auto v = detail::dyadic_grid<long double, p, 0>(grid_refinements);
                std::vector<double> w(v.size());
                for (size_t i = 0; i < v.size(); ++i) {
                    w[i] = v[i];
                }
                return w;
            }

            return detail::dyadic_grid<Real, p, 0>(grid_refinements);
        });
        // Compute the derivative of the refined grid:
        std::future<std::vector<Real>> t1 = std::async(std::launch::async, [&grid_refinements]() {
            if constexpr (std::is_same_v<Real, float>) {
                auto v = detail::dyadic_grid<double, p, 1>(grid_refinements);
                std::vector<float> w(v.size());
                for (size_t i = 0; i < v.size(); ++i) {
                    w[i] = v[i];
                }
                return w;
            }
            else if constexpr (std::is_same_v<Real, double>) {
                auto v = detail::dyadic_grid<long double, p, 1>(grid_refinements);
                std::vector<double> w(v.size());
                for (size_t i = 0; i < v.size(); ++i) {
                    w[i] = v[i];
                }
                return w;
            }

            return detail::dyadic_grid<Real, p, 1>(grid_refinements);
        });

        // if necessary, compute the second and third derivative:
        std::vector<Real> d2ydx2;
        std::vector<Real> d3ydx3;
        if constexpr (p >= 6) {
            std::future<std::vector<Real>> t3 = std::async(std::launch::async, [&grid_refinements]() {
                if constexpr (std::is_same_v<Real, float>) {
                    auto v = detail::dyadic_grid<double, p, 2>(grid_refinements);
                    std::vector<float> w(v.size());
                    for (size_t i = 0; i < v.size(); ++i) {
                        w[i] = v[i];
                    }
                    return w;
                }
                else if constexpr (std::is_same_v<Real, double>) {
                    auto v = detail::dyadic_grid<long double, p, 2>(grid_refinements);
                    std::vector<double> w(v.size());
                    for (size_t i = 0; i < v.size(); ++i) {
                        w[i] = v[i];
                    }
                    return w;
                }

                return detail::dyadic_grid<Real, p, 2>(grid_refinements);
            });

            if constexpr (p >= 10) {
                std::future<std::vector<Real>> t4 = std::async(std::launch::async, [&grid_refinements]() {
                    if constexpr (std::is_same_v<Real, float>) {
                        auto v = detail::dyadic_grid<double, p, 3>(grid_refinements);
                        std::vector<float> w(v.size());
                        for (size_t i = 0; i < v.size(); ++i) {
                            w[i] = v[i];
                        }
                        return w;
                    }
                    else if constexpr (std::is_same_v<Real, double>) {
                        auto v = detail::dyadic_grid<long double, p, 3>(grid_refinements);
                        std::vector<double> w(v.size());
                        for (size_t i = 0; i < v.size(); ++i) {
                            w[i] = v[i];
                        }
                        return w;
                    }

                    return detail::dyadic_grid<Real, p, 3>(grid_refinements); });
                d3ydx3 = t4.get();
            }
            d2ydx2 = t3.get();
        }


        auto y = t0.get();
        auto dydx = t1.get();

        if constexpr (p==2) {
            m_mh = std::make_shared<detail::matched_holder<std::vector<Real>>>(std::move(y), std::move(dydx), grid_refinements);
        }
        if constexpr (p==3) {
            m_lin = std::make_shared<detail::linear_interpolation<std::vector<Real>>>(std::move(y), grid_refinements);
        }
        if constexpr (p == 4 || p == 5) {
            Real dx = Real(1)/(1 << grid_refinements);
            m_cbh = std::make_shared<interpolators::detail::cardinal_cubic_hermite_detail<std::vector<Real>>>(std::move(y), std::move(dydx), Real(0), dx);
        }
        if constexpr (p >= 6 && p <= 9) {
            Real dx = Real(1)/(1 << grid_refinements);
            m_qh = std::make_shared<interpolators::detail::cardinal_quintic_hermite_detail<std::vector<Real>>>(std::move(y), std::move(dydx), std::move(d2ydx2), Real(0), dx);
        }
        if constexpr (p >= 10) {
            Real dx = Real(1)/(1 << grid_refinements);
            m_sh = std::make_shared<interpolators::detail::cardinal_septic_hermite_detail<std::vector<Real>>>(std::move(y), std::move(dydx), std::move(d2ydx2), std::move(d3ydx3), Real(0), dx);
        }
     }


    inline Real operator()(Real x) const {
        if (x <= 0 || x >= 2*p-1) {
            return 0;
        }
        if constexpr (p==2) {
            return m_mh->operator()(x);
        }
        if constexpr (p==3) {
            return m_lin->operator()(x);
        }
        if constexpr (p==4 || p ==5) {
            return m_cbh->unchecked_evaluation(x);
        }
        if constexpr (p >= 6 && p <= 9) {
            return m_qh->unchecked_evaluation(x);
        }
        if constexpr (p >= 10) {
            return m_sh->unchecked_evaluation(x);
        }
    }

    inline Real prime(Real x) const {
        if (x <= 0 || x >= 2*p-1)
        {
            return 0;
        }
        if constexpr (p == 2 || p == 3) {
            throw std::domain_error("The 2 and 3-vanishing moment Daubechies scaling function is not continuously differentiable.");
        }
        if constexpr (p == 4 || p == 5) {
            return m_cbh->unchecked_prime(x);
        }
        if constexpr (p >= 6 && p <= 9) {
            return m_qh->unchecked_prime(x);
        }
        if constexpr (p >= 10) {
            return m_sh->unchecked_prime(x);
        }
    }

    inline Real double_prime(Real x) const {
        if (x <= 0 || x >= 2*p - 1) {
            return Real(0);
        }
        if constexpr (p >= 6 && p <= 9) {
            return m_qh->unchecked_double_prime(x);
        }
    }
    std::pair<int, int> support() const
    {
        return {0, 2*p-1};
    }

private:
    // Need this for p = 2:
    std::shared_ptr<detail::matched_holder<std::vector<Real>>> m_mh;
    // Need this for p = 3:
    std::shared_ptr<detail::linear_interpolation<std::vector<Real>>> m_lin;
    // Need this for p = 4,5:
    std::shared_ptr<interpolators::detail::cardinal_cubic_hermite_detail<std::vector<Real>>> m_cbh;
    // Need this for p = 6,7,8,9:
    std::shared_ptr<interpolators::detail::cardinal_quintic_hermite_detail<std::vector<Real>>> m_qh;
    // Need this for p >= 10:
    std::shared_ptr<interpolators::detail::cardinal_septic_hermite_detail<std::vector<Real>>> m_sh;

    /*
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

    size_t bytes() const
    {
        size_t s = sizeof(*this);
        s += m_v.size()*sizeof(Real);
        return s;
    }
*/
};

}
#endif
