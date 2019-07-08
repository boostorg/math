// Copyright Nick Thompson, 2019
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
#ifndef BOOST_MATH_INTERPOLATORS_THEIS_HPP
#define BOOST_MATH_INTERPOLATORS_THEIS_HPP
#include <memory>
#include <boost/math/interpolators/detail/theis_detail.hpp>

namespace boost { namespace math { namespace interpolators {

template<class RandomAccessContainer>
class theis {
public:

    using Real = typename RandomAccessContainer::value_type;
    theis(RandomAccessContainer&& y, Real const & t0, Real const & h)
     : m_impl(std::make_shared<detail::theis_detail<RandomAccessContainer>>(std::move(y), t0, h))
    {}

    inline Real operator()(Real t) const
    {
        return m_impl->operator()(t);
    }

    inline Real operator[](size_t i) const
    {
        return m_impl->operator[](i);
    }

    RandomAccessContainer&& return_data()
    {
        return m_impl->return_data();
    }


private:
    std::shared_ptr<detail::theis_detail<RandomAccessContainer>> m_impl;
};
}}}
#endif
