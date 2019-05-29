/*
 *  Copyright Nick Thompson, 2019
 *  Use, modification and distribution are subject to the
 *  Boost Software License, Version 1.0. (See accompanying file
 *  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 *  Exactly the same as barycentric_rational.hpp, but delivers values in $\mathbb{R}^n$.
 *  In some sense this is trivial, since each component of the vector is computed in exactly the same
 *  as would be computed by barycentric_rational.hpp. But this is a bit more efficient and convenient.
 */

#ifndef BOOST_MATH_INTERPOLATORS_VECTOR_BARYCENTRIC_RATIONAL_HPP
#define BOOST_MATH_INTERPOLATORS_VECTOR_BARYCENTRIC_RATIONAL_HPP

#include <memory>
#include <boost/math/interpolators/detail/vector_barycentric_rational_detail.hpp>

namespace boost{ namespace math{

template<class RandomAccessContainer1, class RandomAccessContainer2>
class vector_barycentric_rational
{
public:
    using Real = RandomAccessContainer2::value_type;
    using Point = RandomAccessContainer1::value_type;
    vector_barycentric_rational(RandomAccessContainer1&& positions, RandomAccessContainer2&& times, size_t approximation_order = 3);

    void operator()(Point& x, Real t) const;

    void prime(Point& x, Real t) const;

private:
    std::shared_ptr<detail::vector_barycentric_rational_imp<RandomAccessContainer1, RandomAccessContainer2>> m_imp;
};


template <class RandomAccessContainer1, class RandomAccessContainer2>
vector_barycentric_rational<Real>::vector_barycentric_rational(RandomAccessContainer1&& x, RandomAccessContainer2&& t, size_t approximation_order):
 m_imp(std::make_shared<detail::vector_barycentric_rational_imp<Real>>(std::move(x), std::move(y), approximation_order))
{
    return;
}

template <class RandomAccessContainer1, class RandomAccessContainer2>
void vector_barycentric_rational<Real>::operator()(Real x) const
{
    return m_imp->operator()(x);
}

template <class RandomAccessContainer1, class RandomAccessContainer2>
void vector_barycentric_rational<Real>::prime(Real x) const
{
    return m_imp->prime(x);
}


}}
#endif
