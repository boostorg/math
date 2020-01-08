// Copyright Nick Thompson, 2020
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// See: https://blogs.mathworks.com/cleve/2019/04/29/makima-piecewise-cubic-interpolation/
// And: https://doi.org/10.1145/321607.321609

#ifndef BOOST_MATH_INTERPOLATORS_MAKIMA_HPP
#define BOOST_MATH_INTERPOLATORS_MAKIMA_HPP
#include <memory>
#include <boost/math/interpolators/detail/makima_detail.hpp>

namespace boost::math::interpolators {

template<class RandomAccessContainer>
class makima {
public:
    using Real = typename RandomAccessContainer::value_type;

    makima(RandomAccessContainer && x, RandomAccessContainer && y,
           Real left_endpoint_derivative = std::numeric_limits<Real>::quiet_NaN(),
           Real right_endpoint_derivative = std::numeric_limits<Real>::quiet_NaN()) : impl_(std::make_shared<detail::makima_detail<RandomAccessContainer>>(std::move(x), std::move(y), left_endpoint_derivative, right_endpoint_derivative))
    {}

    Real operator()(Real x) const {
        return impl_->operator()(x);
    }

    Real prime(Real x) const {
        return impl_->prime(x);
    }

    friend std::ostream& operator<<(std::ostream & os, const makima & m)
    {
        os << *m.impl_;
        return os;
    }

    void push_back(Real x, Real y) {
        impl_->push_back(x, y);
    }

private:
    std::shared_ptr<detail::makima_detail<RandomAccessContainer>> impl_;
};

}
#endif