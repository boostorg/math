/*
 *  Copyright Nick Thompson, 2019
 *  Use, modification and distribution are subject to the
 *  Boost Software License, Version 1.0. (See accompanying file
 *  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef BOOST_MATH_INTERPOLATORS_VECTOR_BARYCENTRIC_RATIONAL_DETAIL_HPP
#define BOOST_MATH_INTERPOLATORS_VECTOR_BARYCENTRIC_RATIONAL_DETAIL_HPP

#ifndef BOOST_MATH_BUILD_MODULE
#include <cstdint>
#include <cmath>
#include <vector>
#include <utility> // for std::move
#include <limits>
#include <algorithm>
#include <stdexcept>
#endif
#include <boost/math/tools/assert.hpp>

BOOST_MATH_NAMESPACE_BEGIN namespace interpolators{ namespace detail{

template <class TimeContainer, class SpaceContainer>
class vector_barycentric_rational_imp
{
public:
    using Real = typename TimeContainer::value_type;
    using Point = typename SpaceContainer::value_type;

    vector_barycentric_rational_imp(TimeContainer&& t, SpaceContainer&& y, size_t approximation_order);

    void operator()(Point& p, Real t) const;

    void eval_with_prime(Point& x, Point& dxdt, Real t) const;

    // The barycentric weights are only interesting to the unit tests:
    Real weight(size_t i) const { return w_[i]; }

private:

    void calculate_weights(size_t approximation_order);

    TimeContainer t_;
    SpaceContainer y_;
    TimeContainer w_;
};

template <class TimeContainer, class SpaceContainer>
vector_barycentric_rational_imp<TimeContainer, SpaceContainer>::vector_barycentric_rational_imp(TimeContainer&& t, SpaceContainer&& y, size_t approximation_order)
{
    using std::numeric_limits;
    t_ = std::move(t);
    y_ = std::move(y);

    if (t_.size() != y_.size())
    {
        throw std::domain_error("There must be the same number of time points as space points.");
    }
    if (approximation_order >= y_.size())
    {
        throw std::domain_error("Approximation order must be < data length.");
    }
    for (size_t i = 1; i < t_.size(); ++i)
    {
        BOOST_MATH_ASSERT_MSG(t_[i] - t_[i-1] >  (numeric_limits<typename TimeContainer::value_type>::min)(), "The abscissas must be listed in strictly increasing order t[0] < t[1] < ... < t[n-1].");
    }
    calculate_weights(approximation_order);
}


template<class TimeContainer, class SpaceContainer>
void vector_barycentric_rational_imp<TimeContainer, SpaceContainer>::calculate_weights(size_t approximation_order)
{
    using std::abs;
    std::int64_t n = t_.size();
    w_.resize(n, Real(0));
    for(std::int64_t k = 0; k < n; ++k)
    {
        std::int64_t i_min = (std::max)(k - static_cast<std::int64_t>(approximation_order), static_cast<std::int64_t>(0));
        std::int64_t i_max = k;
        if (k >= n - (std::ptrdiff_t)approximation_order)
        {
            i_max = n - approximation_order - 1;
        }

        for(std::int64_t i = i_min; i <= i_max; ++i)
        {
            Real inv_product = 1;
            std::int64_t j_max = (std::min)(static_cast<std::int64_t>(i + approximation_order), static_cast<std::int64_t>(n - 1));
            for(std::int64_t j = i; j <= j_max; ++j)
            {
                if (j == k)
                {
                    continue;
                }
                Real diff = t_[k] - t_[j];
                inv_product *= diff;
            }
            if (i % 2 == 0)
            {
                w_[k] += 1/inv_product;
            }
            else
            {
                w_[k] -= 1/inv_product;
            }
        }
    }
}


template<class TimeContainer, class SpaceContainer>
void vector_barycentric_rational_imp<TimeContainer, SpaceContainer>::operator()(typename SpaceContainer::value_type& p, typename TimeContainer::value_type t) const
{
    for (auto & x : p)
    {
        x = Real(0);
    }
    Real denominator = 0;
    for(size_t i = 0; i < t_.size(); ++i)
    {
        // See associated commentary in the scalar version of this function.
        if (t == t_[i])
        {
            p = y_[i];
            return;
        }
        Real x = w_[i]/(t - t_[i]);
        for (decltype(p.size()) j = 0; j < p.size(); ++j)
        {
            p[j] += x*y_[i][j];
        }
        denominator += x;
    }
    for (decltype(p.size()) j = 0; j < p.size(); ++j)
    {
        p[j] /= denominator;
    }
    return;
}

template<class TimeContainer, class SpaceContainer>
void vector_barycentric_rational_imp<TimeContainer, SpaceContainer>::eval_with_prime(typename SpaceContainer::value_type& x, typename SpaceContainer::value_type& dxdt, typename TimeContainer::value_type t) const
{
    this->operator()(x, t);
    // rk = r(t) - y_k, for the node k nearest t, in a form that does not cancel near t_k:
    size_t nearest = std::lower_bound(t_.begin(), t_.end(), t) - t_.begin();
    if (nearest == t_.size() || (nearest > 0 && t - t_[nearest-1] < t_[nearest] - t))
    {
        --nearest;
    }
    Point numerator;
    Point rk;
    for (decltype(x.size()) i = 0; i < x.size(); ++i)
    {
        numerator[i] = 0;
        rk[i] = 0;
    }
    Real denominator = 0;
    if (t != t_[nearest])
    {
        for (decltype(t_.size()) i = 0; i < t_.size(); ++i)
        {
            Real tw = w_[i]/(t - t_[i]);
            for (decltype(rk.size()) j = 0; j < rk.size(); ++j)
            {
                rk[j] += tw*(y_[i][j] - y_[nearest][j]);
            }
            denominator += tw;
        }
        for (decltype(rk.size()) j = 0; j < rk.size(); ++j)
        {
            rk[j] /= denominator;
        }
    }
    for(decltype(t_.size()) i = 0; i < t_.size(); ++i)
    {
        if (t == t_[i])
        {
            Point sum;
            for (decltype(x.size()) m = 0; m < x.size(); ++m)
            {
                sum[m] = 0;
            }

            for (decltype(t_.size()) j = 0; j < t_.size(); ++j)
            {
                if (j == i)
                {
                    continue;
                }
                for (decltype(sum.size()) k = 0; k < sum.size(); ++k)
                {
                    sum[k] += w_[j]*(y_[i][k] - y_[j][k])/(t_[i] - t_[j]);
                }
            }
            for (decltype(sum.size()) k = 0; k < sum.size(); ++k)
            {
                dxdt[k] = -sum[k]/w_[i];
            }
            return;
        }
        Real tw = w_[i]/(t - t_[i]);
        Point diff;
        for (decltype(diff.size()) j = 0; j < diff.size(); ++j)
        {
            diff[j] = (rk[j] + (y_[nearest][j] - y_[i][j]))/(t-t_[i]);
        }
        for (decltype(diff.size()) j = 0; j < diff.size(); ++j)
        {
            numerator[j] += tw*diff[j];
        }
    }

    for (decltype(dxdt.size()) j = 0; j < dxdt.size(); ++j)
    {
        dxdt[j] = numerator[j]/denominator;
    }
    return;
}

}} BOOST_MATH_NAMESPACE_END
#endif
