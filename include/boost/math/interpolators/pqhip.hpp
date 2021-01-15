// Copyright Nick Thompson, 2021
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// This is a scalar-valued PQUIP interpolator as discussed in
// Corless, A Graduate Introduction to Numerical Methods, problem 8.40.
// It's quite similar to PCHIP, but the PCHIP procedure is applied to the derivatives,
// not to the function values.
// An implementation in Matlab is found here:
// https://web.archive.org/web/20210115165752/http://www.nfillion.com/coderepository/Graduate_Introduction_to_Numerical_Methods/pqhip.m

// The generalized barycentric weights of problem 8.40 of the reference are easily found in Mathematica:
// Apart[1/(z^3*(z - 1)^3)], giving:
// 1/(z-1)^3 - 3/(z-1)^2 + 6/(z-1) - 1/z^3 - 3/z^2 - 6/z,
// which gives
// beta00 = -6, beta01 = -3, beta01 = -1
// beta10 = 6, beta11 = -3, beta11 = 1

#ifndef BOOST_MATH_INTERPOLATORS_PQHIP_HPP
#define BOOST_MATH_INTERPOLATORS_PQHIP_HPP
#include <memory>
#include <boost/math/interpolators/detail/quintic_hermite_detail.hpp>

namespace boost {
namespace math {
namespace interpolators {

template<class RandomAccessContainer>
class pqhip {
public:
    using Real = typename RandomAccessContainer::value_type;

    pqhip(RandomAccessContainer && times, RandomAccessContainer && y, RandomAccessContainer && dydt,
          Real left_endpoint_second_derivative = std::numeric_limits<Real>::quiet_NaN(),
          Real right_endpoint_second_derivative = std::numeric_limits<Real>::quiet_NaN())
    {
        using std::isnan;
        // Validation of inputs is done in the detail class.
        RandomAccessContainer d2ydt2(times.size(), std::numeric_limits<Real>::quiet_NaN());
        if (isnan(left_endpoint_second_derivative))
        {
          //
        }
        else
        {
            d2ydt2[0] = left_endpoint_second_derivative;
        }

        for (decltype(d2ydt2.size()) k = 1; k < d2ydt2.size()-1; ++k) {
            Real hkm1 = times[k] - times[k-1];
            Real dkm1 = (dydt[k] - dydt[k-1])/hkm1;

            Real hk = times[k+1] - times[k];
            Real dk = (dydt[k+1] - dydt[k])/hk;
            Real w1 = 2*hk + hkm1;
            Real w2 = hk + 2*hkm1;
            if ( (dk > 0 && dkm1 < 0) || (dk < 0 && dkm1 > 0) || dk == 0 || dkm1 == 0)
            {
                d2ydt2[k] = 0;
            }
            else
            {
                d2ydt2[k] = (w1+w2)/(w1/dkm1 + w2/dk);
            }

        }
        // Quadratic extrapolation at the other end:
        auto n = d2ydt2.size();
        if (isnan(right_endpoint_second_derivative))
        {
            //
        }
        else
        {
            d2ydt2[n-1] = right_endpoint_second_derivative;
        }
        impl_ = std::make_shared<detail::quintic_hermite_detail<RandomAccessContainer>>(std::move(times), std::move(y), std::move(dydt), std::move(d2ydt2));
    }

    Real operator()(Real x) const {
        return impl_->operator()(x);
    }

    Real prime(Real x) const {
        return impl_->prime(x);
    }

    friend std::ostream& operator<<(std::ostream & os, const pqhip & m)
    {
        os << *m.impl_;
        return os;
    }

    void push_back(Real x, Real y) {
        using std::abs;
        using std::isnan;
        if (x <= impl_->x_.back()) {
             throw std::domain_error("Calling push_back must preserve the monotonicity of the x's");
        }
        impl_->x_.push_back(x);
        impl_->y_.push_back(y);
        impl_->dydx_.push_back(std::numeric_limits<Real>::quiet_NaN());
        auto n = impl_->size();
        impl_->dydx_[n-1] = (impl_->y_[n-1]-impl_->y_[n-2])/(impl_->x_[n-1] - impl_->x_[n-2]);
        // Now fix s_[n-2]:
        auto k = n-2;
        Real hkm1 = impl_->x_[k] - impl_->x_[k-1];
        Real dkm1 = (impl_->y_[k] - impl_->y_[k-1])/hkm1;

        Real hk = impl_->x_[k+1] - impl_->x_[k];
        Real dk = (impl_->y_[k+1] - impl_->y_[k])/hk;
        Real w1 = 2*hk + hkm1;
        Real w2 = hk + 2*hkm1;
        if ( (dk > 0 && dkm1 < 0) || (dk < 0 && dkm1 > 0) || dk == 0 || dkm1 == 0)
        {
            impl_->dydx_[k] = 0;
        }
        else
        {
            impl_->dydx_[k] = (w1+w2)/(w1/dkm1 + w2/dk);
        }
    }

private:
    std::shared_ptr<detail::quintic_hermite_detail<RandomAccessContainer>> impl_;
};

}
}
}
#endif
