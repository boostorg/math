// Copyright Nick Thompson, 2021
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
#ifndef BOOST_MATH_INTERPOLATORS_BEZIER_POLYNOMIAL_HPP
#define BOOST_MATH_INTERPOLATORS_BEZIER_POLYNOMIAL_HPP
#include <memory>
#include <boost/math/interpolators/detail/bezier_polynomial_detail.hpp>

namespace boost::math::interpolators {

template <class RandomAccessContainer>
class bezier_polynomial
{
public:
    using Point = typename RandomAccessContainer::value_type;
    using Real = typename Point::value_type;
    using Z = typename RandomAccessContainer::size_type;

    bezier_polynomial(RandomAccessContainer && control_points)
    : m_imp(std::make_shared<detail::bezier_polynomial_imp<RandomAccessContainer>>(std::move(control_points)))
    {
    }

    inline Point operator()(Real t) const
    {
        return (*m_imp)(t);
    }

    void edit_control_point(Point const & p, Z index)
    {
        m_imp->edit_control_point(p, index);
    }

private:
    std::shared_ptr<detail::bezier_polynomial_imp<RandomAccessContainer>> m_imp;
};

}
#endif
