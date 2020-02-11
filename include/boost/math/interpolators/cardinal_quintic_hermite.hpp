#ifndef BOOST_MATH_INTERPOLATORS_CARDINAL_QUINTIC_HERMITE_HPP
#define BOOST_MATH_INTERPOLATORS_CARDINAL_QUINTIC_HERMITE_HPP
#include <algorithm>
#include <stdexcept>
#include <memory>
#include <boost/math/interpolators/detail/cardinal_quintic_hermite_detail.hpp>

namespace boost::math::interpolators {

template<class RandomAccessContainer>
class cardinal_quintic_hermite {
public:
    using Real = typename RandomAccessContainer::value_type;
    cardinal_quintic_hermite(RandomAccessContainer && y, RandomAccessContainer && dydx, RandomAccessContainer && d2ydx2, Real x0, Real dx)
     : impl_(std::make_shared<detail::cardinal_quintic_hermite_detail<RandomAccessContainer>>(std::move(y), std::move(dydx), std::move(d2ydx2), x0, dx))
    {}

    inline Real operator()(Real x) const {
        return impl_->operator()(x);
    }

    inline Real prime(Real x) const {
        return impl_->prime(x);
    }

    inline Real double_prime(Real x) const {
        return impl_->double_prime(x);
    }

private:
    std::shared_ptr<detail::cardinal_quintic_hermite_detail<RandomAccessContainer>> impl_;
};

template<class RandomAccessContainer>
class cardinal_quintic_hermite_aos {
public:
    using Point = typename RandomAccessContainer::value_type;
    using Real = typename Point::value_type;
    cardinal_quintic_hermite_aos(RandomAccessContainer && data, Real x0, Real dx)
     : impl_(std::make_shared<detail::cardinal_quintic_hermite_detail_aos<RandomAccessContainer>>(std::move(data), x0, dx))
    {}

    inline Real operator()(Real x) const {
        return impl_->operator()(x);
    }

    inline Real prime(Real x) const {
        return impl_->prime(x);
    }

private:
    std::shared_ptr<detail::cardinal_quintic_hermite_detail_aos<RandomAccessContainer>> impl_;
};

}
#endif 