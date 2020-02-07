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
#include <boost/multiprecision/float128.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/detail/daubechies_scaling_integer_grid.hpp>
#include <boost/math/filters/daubechies.hpp>
#include <boost/math/interpolators/detail/cubic_hermite_detail.hpp>
#include <boost/math/interpolators/detail/quintic_hermite_detail.hpp>


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
            if (delivery_idx >= (int64_t) v.size())
            {
                std::cerr << "Delivery index out of range!\n";
                continue;
            }
            Real term = 0;
            for (int64_t l = 0; l < (int64_t) c.size(); ++l) {
                int64_t idx = k*(1 << (j_max - j + 1)) - l*(1 << j_max);
                if (idx < 0) {
                    break;
                }
                if (idx < (int64_t) v.size()) {
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

template<class RandomAccessContainer>
class matched_holder {
public:
    using Real = typename RandomAccessContainer::value_type;

    matched_holder(RandomAccessContainer && y, RandomAccessContainer && dydx, int grid_refinements) : y_{std::move(y)}, dydx_{std::move(dydx)}
    {
        h_ = Real(1)/(1<< grid_refinements);
    }
    
    Real operator()(Real x) const {
        using std::log;
        using std::floor;
        using std::sqrt;
        using std::pow;
        if (x <= 0 || x >= 3) {
            return 0;
        }
        // This is the exact Holder exponent, but it's pessimistic almost everywhere!
        // It's only exactly right at dyadic rationals.
        Real const alpha = 2 - log(1+sqrt(Real(3)))/log(Real(2));
        // So we're gonna make the graph dip a little harder; this will capture more of the self-similar behavior:
        //Real constexpr const alpha = Real(3)/Real(10);
        int64_t i = static_cast<int64_t>(floor(x/h_));
        Real t = (x- i*h_)/h_;
        Real v = y_[i];
        Real dphi = dydx_[i+1]*h_;
        v += (dphi - alpha*(y_[i+1] - y_[i]))*t/(1-alpha);
        v += (-dphi + y_[i+1] - y_[i])*pow(t, alpha)/(1-alpha);
        return v;
    }

private:
    Real h_;
    RandomAccessContainer y_;
    RandomAccessContainer dydx_;
};


template<class RandomAccessContainer>
class linear_interpolation {
public:
    using Real = typename RandomAccessContainer::value_type;

    linear_interpolation(RandomAccessContainer && y,  int grid_refinements) : y_{std::move(y)}
    {
        grid_refinements_ = grid_refinements;
    }

    Real operator()(Real x) const {
        using std::floor;
        if (x <= 0 || x >= 5) {
            return 0;
        }
        Real y = x*(1<<grid_refinements_);
        Real k = floor(y);

        size_t kk = static_cast<size_t>(k);

        Real t = y - k;
        return (1-t)*y_[kk] + t*y_[kk+1];
    }

private:
    int grid_refinements_;
    RandomAccessContainer y_;
};

}

template<class Real, int p>
class daubechies_scaling {
public:
    daubechies_scaling(int grid_refinements = -1)
    {
        static_assert(p <= 15, "Scaling functions only implements up to p = 15");
        using boost::multiprecision::float128;
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
                std::array<int, 16> r{-1, -1, 22, 22, 22, 22, 22, 22, 22, 22, 21, 21, 20, 19, 18, 18};
                grid_refinements = r[p];
            }
            else
            {
                grid_refinements = 22;
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

        // if necessary, compute the second derivative:
        std::vector<Real> d2ydx2;
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

                return detail::dyadic_grid<Real, p, 2>(grid_refinements); });
            d2ydx2 = t3.get();
        }

        auto y = t0.get();
        auto dydx = t1.get();



        // Storing the vector of x's is unnecessary; it's only because I haven't implemented an equispaced cubic Hermite interpolator:
        std::vector<Real> x(y.size());
        Real h = Real(2*p-1)/(x.size()-1);
        for (size_t i = 0; i < x.size(); ++i) {
            x[i] = i*h;
        }

        if constexpr (p==2) {
            m_mh = std::make_shared<detail::matched_holder<std::vector<Real>>>(std::move(y), std::move(dydx), grid_refinements);
        }
        if constexpr (p==3) {
            m_lin = std::make_shared<detail::linear_interpolation<std::vector<Real>>>(std::move(y), grid_refinements);
        }
        if constexpr (p == 4 || p == 5) {
            m_cbh = std::make_shared<interpolators::detail::cubic_hermite_detail<std::vector<Real>>>(std::move(x), std::move(y), std::move(dydx));
        }
        else if constexpr (p >= 6) {
            m_qh = std::make_shared<interpolators::detail::quintic_hermite_detail<std::vector<Real>>>(std::move(x), std::move(y), std::move(dydx), std::move(d2ydx2));
        }
     }


    Real operator()(Real x) const {
        if constexpr (p==2) {
            return m_mh->operator()(x);
        }
        if constexpr (p==3) {
            return m_lin->operator()(x);
        }
        if constexpr (p==4 || p ==5) {
            return m_cbh->operator()(x);
        }
        else if constexpr (p >= 6) {
            return m_qh->operator()(x);
        }
    }

    Real prime(Real x) const {
        if constexpr (p==2 || p == 3) {
            throw std::domain_error("The 2 and 3-vanishing moment Daubechies scaling function is not continuously differentiable.");
        }
        if constexpr (p==4 || p ==5) {
            return m_cbh->prime(x);
        }
        else if constexpr (p >= 6) {
            return m_qh->prime(x);
        }
    }

    std::pair<int, int> support() const {
        return {0, 2*p-1};
    }

private:
    size_t m_levels;
    // Need this for p = 2:
    std::shared_ptr<detail::matched_holder<std::vector<Real>>> m_mh;
    // Need this for p = 3:
    std::shared_ptr<detail::linear_interpolation<std::vector<Real>>> m_lin;
    // need this for p = 4,5:
    std::shared_ptr<interpolators::detail::cubic_hermite_detail<std::vector<Real>>> m_cbh;
    // need this for p >= 6:
    std::shared_ptr<interpolators::detail::quintic_hermite_detail<std::vector<Real>>> m_qh;

    /*Real constant_interpolation(Real x) const {
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
    std::vector<Real> m_v_dbl_prime;*/
};

}
#endif
