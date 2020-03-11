/*
 * Copyright Nick Thompson, 2020
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef BOOST_MATH_SPECIAL_DAUBECHIES_WAVELET_HPP
#define BOOST_MATH_SPECIAL_DAUBECHIES_WAVELET_HPP
#include <vector>
#include <array>
#include <cmath>
#include <thread>
#include <future>
#include <iostream>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/detail/daubechies_scaling_integer_grid.hpp>
#include <boost/math/special_functions/daubechies_scaling.hpp>
#include <boost/math/filters/daubechies.hpp>
#include <boost/math/interpolators/detail/cubic_hermite_detail.hpp>
#include <boost/math/interpolators/detail/quintic_hermite_detail.hpp>
#include <boost/math/interpolators/detail/septic_hermite_detail.hpp>

namespace boost::math {

template<class Real, int p, int order>
std::vector<Real> wavelet_dyadic_grid(int64_t j_max)
{
    auto phijk = dyadic_grid<Real, p, order>(j_max);
}


template<class Real, int p>
class daubechies_wavelet {
public:
    daubechies_wavelet(int grid_refinements = -1)
    {
        static_assert(p < 20, "Daubechies wavelets are only implemented for p < 20.");
        static_assert(p > 0, "Daubechies wavelets must have at least 1 vanishing moment.");
        if constexpr (p == 1)
        {
            return;
        }
        else
        {
        if (grid_refinements < 0)
        {
            if (std::is_same_v<Real, float>)
            {
                if (grid_refinements == -2)
                {
                    // Control absolute error:
                    //                          p= 2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19
                    std::array<int, 20> r{-1, -1, 18, 19, 16, 11,  8,  7,  7,  7,  5,  5,  4,  4,  4,  4,  3,  3,  3,  3};
                    grid_refinements = r[p];
                }
                else
                {
                    // Control relative error:
                    //                          p= 2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19
                    std::array<int, 20> r{-1, -1, 21, 21, 21, 17, 16, 15, 14, 13, 12, 11, 11, 11, 11, 11, 11, 11, 11, 11};
                    grid_refinements = r[p];
                }
            }
            else if (std::is_same_v<Real, double>)
            {
                //                          p= 2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19
                std::array<int, 20> r{-1, -1, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 20, 19, 18, 18, 18, 18, 18, 18};
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
            if constexpr (std::is_same_v<Real, float>)
            {
                auto v = dyadic_grid<double, p, 0>(grid_refinements);
                std::vector<float> w(v.size());
                for (size_t i = 0; i < v.size(); ++i)
                {
                    w[i] = static_cast<float>(v[i]);
                }
                return w;
            }
            else if constexpr (std::is_same_v<Real, double>)
            {
                auto v = dyadic_grid<long double, p, 0>(grid_refinements);
                std::vector<double> w(v.size());
                for (size_t i = 0; i < v.size(); ++i)
                {
                    w[i] = static_cast<double>(v[i]);
                }
                return w;
            }

            return dyadic_grid<Real, p, 0>(grid_refinements);
        });
        // Compute the derivative of the refined grid:
        std::future<std::vector<Real>> t1 = std::async(std::launch::async, [&grid_refinements]() {
            if constexpr (std::is_same_v<Real, float>)
            {
                auto v = dyadic_grid<double, p, 1>(grid_refinements);
                std::vector<float> w(v.size());
                for (size_t i = 0; i < v.size(); ++i)
                {
                    w[i] = static_cast<float>(v[i]);
                }
                return w;
            }
            else if constexpr (std::is_same_v<Real, double>)
            {
                auto v = dyadic_grid<long double, p, 1>(grid_refinements);
                std::vector<double> w(v.size());
                for (size_t i = 0; i < v.size(); ++i)
                {
                    w[i] = static_cast<double>(v[i]);
                }
                return w;
            }

            return dyadic_grid<Real, p, 1>(grid_refinements);
        });

        // if necessary, compute the second and third derivative:
        std::vector<Real> d2ydx2;
        std::vector<Real> d3ydx3;
        if constexpr (p >= 6) {
            std::future<std::vector<Real>> t3 = std::async(std::launch::async, [&grid_refinements]() {
                if constexpr (std::is_same_v<Real, float>)
                {
                    auto v = dyadic_grid<double, p, 2>(grid_refinements);
                    std::vector<float> w(v.size());
                    for (size_t i = 0; i < v.size(); ++i)
                    {
                        w[i] = static_cast<float>(v[i]);
                    }
                    return w;
                }
                else if constexpr (std::is_same_v<Real, double>)
                {
                    auto v = dyadic_grid<long double, p, 2>(grid_refinements);
                    std::vector<double> w(v.size());
                    for (size_t i = 0; i < v.size(); ++i)
                    {
                        w[i] = static_cast<double>(v[i]);
                    }
                    return w;
                }

                return dyadic_grid<Real, p, 2>(grid_refinements);
            });

            if constexpr (p >= 10) {
                std::future<std::vector<Real>> t4 = std::async(std::launch::async, [&grid_refinements]() {
                    if constexpr (std::is_same_v<Real, float>)
                    {
                        auto v = dyadic_grid<double, p, 3>(grid_refinements);
                        std::vector<float> w(v.size());
                        for (size_t i = 0; i < v.size(); ++i)
                        {
                            w[i] = static_cast<float>(v[i]);
                        }
                        return w;
                    }
                    else if constexpr (std::is_same_v<Real, double>)
                    {
                        auto v = dyadic_grid<long double, p, 3>(grid_refinements);
                        std::vector<double> w(v.size());
                        for (size_t i = 0; i < v.size(); ++i)
                        {
                            w[i] = static_cast<double>(v[i]);
                        }
                        return w;
                    }

                    return dyadic_grid<Real, p, 3>(grid_refinements);
                });
                d3ydx3 = t4.get();
            }
            d2ydx2 = t3.get();
        }


        auto y = t0.get();
        auto dydx = t1.get();

        if constexpr (p==2)
        {
            std::vector<std::array<Real, 2>> data(y.size());
            for (size_t i = 0; i < y.size(); ++i)
            {
                data[i][0] = y[i];
                data[i][1] = dydx[i];
            }
            m_mh = std::make_shared<detail::matched_holder_aos<std::vector<std::array<Real,2>>>>(std::move(data), grid_refinements, Real(-p+1));
        }
        if constexpr (p==3)
        {
            std::vector<std::array<Real, 2>> data(y.size());
            for (size_t i = 0; i < y.size(); ++i)
            {
                data[i][0] = y[i];
                data[i][1] = dydx[i];
            }
            m_lin = std::make_shared<detail::linear_interpolation_aos<std::vector<std::array<Real, 2>>>>(std::move(data), grid_refinements, Real(-p+1));
        }
        if constexpr (p == 4 || p == 5)
        {
            Real dx = Real(1)/(1 << grid_refinements);
            std::vector<std::array<Real, 2>> data(y.size());
            for (size_t i = 0; i < y.size(); ++i)
            {
                data[i][0] = y[i];
                data[i][1] = dydx[i];
            }
            m_cbh = std::make_shared<interpolators::detail::cardinal_cubic_hermite_detail_aos<std::vector<std::array<Real,2>>>>(std::move(data), Real(-p+1), dx);
        }
        if constexpr (p >= 6 && p <= 9)
        {
            Real dx = Real(1)/(1 << grid_refinements);
            std::vector<std::array<Real, 3>> data(y.size());
            for (size_t i = 0; i < y.size(); ++i)
            {
                data[i][0] = y[i];
                data[i][1] = dydx[i];
                data[i][2] = d2ydx2[i];
            }

            m_qh = std::make_shared<interpolators::detail::cardinal_quintic_hermite_detail_aos<std::vector<std::array<Real,3>>>>(std::move(data), Real(-p+1), dx);
        }
        if constexpr (p >= 10)
        {
            Real dx = Real(1)/(1 << grid_refinements);
            std::vector<std::array<Real, 4>> data(y.size());
            for (size_t i = 0; i < y.size(); ++i)
            {
                data[i][0] = y[i];
                data[i][1] = dydx[i];
                data[i][2] = d2ydx2[i];
                data[i][3] = d3ydx3[i];
            }
            m_sh = std::make_shared<interpolators::detail::cardinal_septic_hermite_detail_aos<std::vector<std::array<Real, 4>>>>(std::move(data), Real(-p+1), dx);
        }
        }
     }


    inline Real operator()(Real x) const
    {
        if (x <= -p + 1 || x >= p)
        {
            return 0;
        }
        if constexpr (p==1)
        {
            if (x < Real(1)/Real(2))
            {
                return 1;
            }
            else if (x == Real(1)/Real(2))
            {
                return 0;
            }
            return -1;
        }
        if constexpr (p==2)
        {
            return m_mh->operator()(x);
        }
        if constexpr (p==3)
        {
            return m_lin->operator()(x);
        }
        if constexpr (p==4 || p ==5)
        {
            return m_cbh->unchecked_evaluation(x);
        }
        if constexpr (p >= 6 && p <= 9)
        {
            return m_qh->unchecked_evaluation(x);
        }
        if constexpr (p >= 10)
        {
            return m_sh->unchecked_evaluation(x);
        }
    }

    inline Real prime(Real x) const
    {
        static_assert(p > 2, "The 3-vanishing moment Daubechies wavelet is the first which is continuously differentiable.");
        if (x <= -p + 1 || x >= p)
        {
            return 0;
        }
        if constexpr (p == 3)
        {
            return m_lin->prime(x);
        }
        if constexpr (p == 4 || p == 5)
        {
            return m_cbh->unchecked_prime(x);
        }
        if constexpr (p >= 6 && p <= 9)
        {
            return m_qh->unchecked_prime(x);
        }
        if constexpr (p >= 10)
        {
            return m_sh->unchecked_prime(x);
        }
    }

    inline Real double_prime(Real x) const
    {
        static_assert(p >= 6, "Second derivatives of Daubechies wavelets require at least 6 vanishing moments.");
        if (x <= -p + 1 || x >= p)
        {
            return Real(0);
        }
        if constexpr (p >= 6 && p <= 9)
        {
            return m_qh->unchecked_double_prime(x);
        }
        if constexpr (p >= 10)
        {
            return m_sh->unchecked_double_prime(x);
        }
    }

    std::pair<Real, Real> support() const
    {
        return {Real(-p+1), Real(p)};
    }

    int64_t bytes() const
    {
        if constexpr (p==1)
        {
            return sizeof(this);
        }
        if constexpr (p==2)
        {
            return m_mh->bytes() + sizeof(m_mh);
        }
        if constexpr (p == 3)
        {
            return m_lin->bytes() + sizeof(m_lin);
        }
        if constexpr (p == 4 || p == 5)
        {
            return m_cbh->bytes() + sizeof(m_cbh);
        }
        if constexpr (p >= 6 && p <= 9)
        {
            return m_qh->bytes() + sizeof(m_qh);
        }
        if constexpr (p >= 10)
        {
            return m_sh->bytes() + sizeof(m_sh);
        }

        return -1;
    }

private:
    // Need this for p = 2:
    std::shared_ptr<detail::matched_holder_aos<std::vector<std::array<Real,2>>>> m_mh;
    // Need this for p = 3:
    std::shared_ptr<detail::linear_interpolation_aos<std::vector<std::array<Real, 2>>>> m_lin;
    // Need this for p = 4,5:
    std::shared_ptr<interpolators::detail::cardinal_cubic_hermite_detail_aos<std::vector<std::array<Real, 2>>>> m_cbh;
    // Need this for p = 6,7,8,9:
    std::shared_ptr<interpolators::detail::cardinal_quintic_hermite_detail_aos<std::vector<std::array<Real, 3>>>> m_qh;
    // Need this for p >= 10:
    std::shared_ptr<interpolators::detail::cardinal_septic_hermite_detail_aos<std::vector<std::array<Real, 4>>>> m_sh;
};

}
#endif
