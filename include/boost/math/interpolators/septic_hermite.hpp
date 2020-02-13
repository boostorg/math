#ifndef BOOST_MATH_INTERPOLATORS_SEPTIC_HERMITE_HPP
#define BOOST_MATH_INTERPOLATORS_SEPTIC_HERMITE_HPP
#include <algorithm>
#include <stdexcept>
#include <memory>
#include <boost/math/interpolators/detail/septic_hermite_detail.hpp>

namespace boost::math::interpolators {

template<class RandomAccessContainer>
class septic_hermite {
public:
    using Real = typename RandomAccessContainer::value_type;
    septic_hermite(RandomAccessContainer && x, RandomAccessContainer && y, RandomAccessContainer && dydx, 
                    RandomAccessContainer && d2ydx2, RandomAccessContainer && d3ydx3)
     : impl_(std::make_shared<detail::septic_hermite_detail<RandomAccessContainer>>(std::move(x), 
     std::move(y), std::move(dydx), std::move(d2ydx2), std::move(d3ydx3)))
    {}

    Real operator()(Real x) const {
        return impl_->operator()(x);
    }

    Real prime(Real x) const {
        return impl_->prime(x);
    }

    friend std::ostream& operator<<(std::ostream & os, const septic_hermite & m)
    {
        os << *m.impl_;
        return os;
    }

private:
    std::shared_ptr<detail::septic_hermite_detail<RandomAccessContainer>> impl_;
};

template<class RandomAccessContainer>
class cardinal_septic_hermite {
public:
    using Real = typename RandomAccessContainer::value_type;
    cardinal_septic_hermite(RandomAccessContainer && y, RandomAccessContainer && dydx,
                    RandomAccessContainer && d2ydx2, RandomAccessContainer && d3ydx3, Real x0, Real dx)
     : impl_(std::make_shared<detail::cardinal_septic_hermite_detail<RandomAccessContainer>>(
     std::move(y), std::move(dydx), std::move(d2ydx2), std::move(d3ydx3), x0, dx))
    {}

    Real operator()(Real x) const {
        return impl_->operator()(x);
    }

    Real prime(Real x) const {
        return impl_->prime(x);
    }

private:
    std::shared_ptr<detail::cardinal_septic_hermite_detail<RandomAccessContainer>> impl_;
};

}
#endif 