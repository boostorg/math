/*
 *  Copyright Nick Thompson, 2019
 *  Use, modification and distribution are subject to the
 *  Boost Software License, Version 1.0. (See accompanying file
 *  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef BOOST_MATH_INTERPOLATORS_VECTOR_BARYCENTRIC_RATIONAL_DETAIL_HPP
#define BOOST_MATH_INTERPOLATORS_VECTOR_BARYCENTRIC_RATIONAL_DETAIL_HPP

#include <vector>
#include <utility> // for std::move
#include <boost/assert.hpp>

namespace boost{ namespace math{ namespace detail{

template <class TimeContainer, class SpaceContainer>
class vector_barycentric_rational_imp
{
public:
    using Real = typename TimeContainer::value_type;
    using Point = typename SpaceContainer::value_type;

    vector_barycentric_rational_imp(TimeContainer&& t, SpaceContainer&& y, size_t approximation_order);

    void operator()(Point& p, Real t) const;

    /*Real prime(Real x) const;

    // The barycentric weights are not really that interesting; except to the unit tests!
    Real weight(size_t i) const { return m_w[i]; }

    std::vector<Real>&& return_x()
    {
        return std::move(m_x);
    }

    std::vector<Real>&& return_y()
    {
        return std::move(m_y);
    }*/

private:

    void calculate_weights(size_t approximation_order);

    TimeContainer t_;
    SpaceContainer y_;
    TimeContainer w_;
};

template <class TimeContainer, class SpaceContainer>
vector_barycentric_rational_imp<TimeContainer, SpaceContainer>::vector_barycentric_rational_imp(TimeContainer&& t, SpaceContainer&& y, size_t approximation_order)
{
    using Real = typename TimeContainer::value_type;
    using std::numeric_limits;
    t_ = std::move(t);
    y_ = std::move(y);
    BOOST_ASSERT_MSG(t_.size() == y_.size(), "There must be the same number of time points as space points.");
    BOOST_ASSERT_MSG(approximation_order < y_.size(), "Approximation order must be < data length.");
    for (size_t i = 1; i < t_.size(); ++i) {
        BOOST_ASSERT_MSG(t_[i] - t_[i-1] >  (numeric_limits<Real>::min)(), "The abscissas must be listed in strictly increasing order t[0] < t[1] < ... < t[n-1].");
    }
    calculate_weights(approximation_order);
}


template<class TimeContainer, class SpaceContainer>
void vector_barycentric_rational_imp<TimeContainer, SpaceContainer>::calculate_weights(size_t approximation_order)
{
    using Real = typename TimeContainer::value_type;
    using std::abs;
    int64_t n = t_.size();
    w_.resize(n, Real(0));
    for(int64_t k = 0; k < n; ++k)
    {
        int64_t i_min = (std::max)(k - (int64_t) approximation_order, (int64_t) 0);
        int64_t i_max = k;
        if (k >= n - (std::ptrdiff_t)approximation_order)
        {
            i_max = n - approximation_order - 1;
        }

        for(int64_t i = i_min; i <= i_max; ++i)
        {
            Real inv_product = 1;
            int64_t j_max = (std::min)(static_cast<int64_t>(i + approximation_order), static_cast<int64_t>(n - 1));
            for(int64_t j = i; j <= j_max; ++j)
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
                t_[k] += 1/inv_product;
            }
            else
            {
                t_[k] -= 1/inv_product;
            }
        }
    }
}


template<class TimeContainer, class SpaceContainer>
void vector_barycentric_rational_imp<TimeContainer, SpaceContainer>::operator()(typename SpaceContainer::value_type& p, typename TimeContainer::value_type t) const
{
    using Real = typename TimeContainer::value_type;
    //using Point = typename SpaceContainer::value_type;
    for (auto & x : p) {
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
        p += x*y_[i];
        denominator += x;
    }
    p /= denominator;
    return;
}


template<class TimeContainer, class SpaceContainer>
void vector_barycentric_rational_imp<TimeContainer, SpaceContainer>::prime(typename SpaceContainer::value_type& p, typename TimeContainer::value_type t) const
{
    Real rx = this->operator()(x);
    Real numerator = 0;
    Real denominator = 0;
    for(size_t i = 0; i < m_x.size(); ++i)
    {
        if (x == m_x[i])
        {
            Real sum = 0;
            for (size_t j = 0; j < m_x.size(); ++j)
            {
                if (j == i)
                {
                    continue;
                }
                sum += m_w[j]*(m_y[i] - m_y[j])/(m_x[i] - m_x[j]);
            }
            return -sum/m_w[i];
        }
        Real t = m_w[i]/(x - m_x[i]);
        Real diff = (rx - m_y[i])/(x-m_x[i]);
        numerator += t*diff;
        denominator += t;
    }

    return numerator/denominator;
}

}}}
#endif
