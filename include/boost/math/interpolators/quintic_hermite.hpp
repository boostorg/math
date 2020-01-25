#ifndef BOOST_MATH_INTERPOLATORS_QUINTIC_HERMITE_HPP
#define BOOST_MATH_INTERPOLATORS_QUINTIC_HERMITE_HPP
#include <algorithm>
#include <stdexcept>
#include <memory>
#include <boost/math/interpolators/detail/quintic_hermite_detail.hpp>

namespace boost::math::interpolators {

template<class RandomAccessContainer>
class quintic_hermite {
public:
    using Real = typename RandomAccessContainer::value_type;
    quintic_hermite(RandomAccessContainer && x, RandomAccessContainer && y, RandomAccessContainer && dydx, RandomAccessContainer && d2ydx2)
     : impl_(std::make_shared<detail::quintic_hermite_detail<RandomAccessContainer>>(std::move(x), std::move(y), std::move(dydx), std::move(d2ydx2)))
    {}

    Real operator()(Real x) const {
        return impl_->operator()(x);
    }

    Real prime(Real x) const {
        return impl_->prime(x);
    }

    friend std::ostream& operator<<(std::ostream & os, const quintic_hermite & m)
    {
        os << *m.impl_;
        return os;
    }

    void push_back(Real x, Real y, Real dydx, Real d2ydx2) {
        impl_->push_back(x, y, dydx, d2ydx2);
    }


private:
    std::shared_ptr<detail::quintic_hermite_detail<RandomAccessContainer>> impl_;
};
}
#endif 