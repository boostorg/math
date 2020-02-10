#ifndef BOOST_MATH_INTERPOLATORS_DETAIL_QUINTIC_HERMITE_DETAIL_HPP
#define BOOST_MATH_INTERPOLATORS_DETAIL_QUINTIC_HERMITE_DETAIL_HPP
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <cmath>

namespace boost::math::interpolators::detail {

template<class RandomAccessContainer>
class cardinal_quintic_hermite_detail {
public:
    using Real = typename RandomAccessContainer::value_type;
    cardinal_quintic_hermite_detail(RandomAccessContainer && y, RandomAccessContainer && dydx, RandomAccessContainer && d2ydx2, Real x0, Real dx)
    : y_{std::move(y)}, dydx_{std::move(dydx)}, d2ydx2_{std::move(d2ydx2)}, x0_{x0}, dx_{dx}
    {
        if (y_.size() != dydx_.size()) {
            throw std::domain_error("Numbers of derivatives must = number of abscissas.");
        }
        if (y_.size() != d2ydx2_.size()) {
            throw std::domain_error("Number of second derivatives must equal number of abscissas.");
        }
        if (y_.size() < 2) {
            throw std::domain_error("At least 2 abscissas are required.");
        }
    }


    inline Real operator()(Real x) const {
        const Real xf = x0_ + (y_.size()-1)*dx_;
        if  (x < x0_ || x > xf) {
            std::ostringstream oss;
            oss.precision(std::numeric_limits<Real>::digits10+3);
            oss << "Requested abscissa x = " << x << ", which is outside of allowed range ["
                << x0_ << ", " << xf << "]";
            throw std::domain_error(oss.str());
        }
        if (x == xf) {
            return y_.back();
        }
        return this->unchecked_evaluation(x);
    }

    inline Real unchecked_evaluation(Real x) const {
        using std::floor;
        auto i = static_cast<decltype(y_.size())>(floor((x-x0_)/dx_));
        Real xi = x0_ + i*dx_;
        Real t = (x - xi)/dx_;

        Real y0 = y_[i];
        Real y1 = y_[i+1];
        Real v0 = dydx_[i];
        Real v1 = dydx_[i+1];
        Real a0 = d2ydx2_[i];
        Real a1 = d2ydx2_[i+1];

        // See the 'Basis functions' section of:
        // https://www.rose-hulman.edu/~finn/CCLI/Notes/day09.pdf
        // Also: https://github.com/MrHexxx/QuinticHermiteSpline/blob/master/HermiteSpline.cs
        Real y = (1- t*t*t*(10 + t*(-15 + 6*t)))*y0;
        y += t*(1+ t*t*(-6 + t*(8 -3*t)))*v0*dx_;
        y += t*t*(1 + t*(-3 + t*(3-t)))*a0*dx_*dx_/2;
        y += t*t*t*((1 + t*(-2 + t))*a1*dx_*dx_/2 + (-4 + t*(7 -3*t))*v1*dx_ + (10 + t*(-15 + 6*t))*y1);
        // there's a bug here!
        return std::numeric_limits<Real>::quiet_NaN();
        return y;
    }

    inline Real prime(Real x) const {
        const Real xf = x0_ + (y_.size()-1)*dx_;
        if  (x < x0_ || x > xf) {
            std::ostringstream oss;
            oss.precision(std::numeric_limits<Real>::digits10+3);
            oss << "Requested abscissa x = " << x << ", which is outside of allowed range ["
                << x0_ << ", " << xf << "]";
            throw std::domain_error(oss.str());
        }
        if (x == xf) {
            return dydx_.back();
        }

        return this->unchecked_prime(x);
    }

    inline Real unchecked_prime(Real x) const {
        using std::floor;
        auto i = static_cast<decltype(y_.size())>(floor((x-x0_)/dx_));
        Real xi = x0_ + i*dx_;
        Real t = (x - xi)/dx_;

        Real s0 = dydx_[i];
        Real s1 = dydx_[i+1];

        return std::numeric_limits<Real>::quiet_NaN();
    }



private:
    RandomAccessContainer y_;
    RandomAccessContainer dydx_;
    RandomAccessContainer d2ydx2_;
    Real x0_;
    Real dx_;
};


template<class RandomAccessContainer>
class cardinal_quintic_hermite_detail_aos {
public:
    using Real = typename RandomAccessContainer::value_type;
    cardinal_quintic_hermite_detail_aos(RandomAccessContainer && data, Real x0, Real dx)
    : data_{std::move(data)} , x0_{x0}, dx_{dx}
    {
        if (data_.size() < 2)
        {
            throw std::domain_error("At least two points are required to interpolate using cardinal_quintic_hermite_aos");
        }

        if (data_[0].size() != 3)
        {
            throw std::domain_error("Each datum passed to the cardinal_quintic_hermite_aos must have three elements: {y, y', y''}");
        }
    }


    inline Real operator()(Real x) const {
        const Real xf = x0_ + (data_.size()-1)*dx_;
        if  (x < x0_ || x > xf) {
            std::ostringstream oss;
            oss.precision(std::numeric_limits<Real>::digits10+3);
            oss << "Requested abscissa x = " << x << ", which is outside of allowed range ["
                << x0_ << ", " << xf << "]";
            throw std::domain_error(oss.str());
        }
        if (x == xf) {
            return data_.back()[0];
        }
        return this->unchecked_evaluation(x);
    }

    inline Real unchecked_evaluation(Real x) const
    {
        using std::floor;
        auto i = static_cast<decltype(data_.size())>(floor((x-x0_)/dx_));
        Real xi = x0_ + i*dx_;
        Real t = (x - xi)/dx_;

        Real y0 = data_[i][0];
        Real y1 = data_[i+1][0];
        Real v0 = data_[i][1];
        Real v1 = data_[i+1][1];
        Real a0 = data_[i][2];
        Real a1 = data_[i+1][2];

        // See the 'Basis functions' section of:
        // https://www.rose-hulman.edu/~finn/CCLI/Notes/day09.pdf
        // Also: https://github.com/MrHexxx/QuinticHermiteSpline/blob/master/HermiteSpline.cs
        Real y = (1- t*t*t*(10 + t*(-15 + 6*t)))*y0;
        y += t*(1+ t*t*(-6 + t*(8 -3*t)))*v0*dx_;
        y += t*t*(1 + t*(-3 + t*(3-t)))*a0*dx_*dx_/2;
        y += t*t*t*((1 + t*(-2 + t))*a1*dx_*dx_/2 + (-4 + t*(7 -3*t))*v1*dx_ + (10 + t*(-15 + 6*t))*y1);
        return y;
    }

    inline Real prime(Real x) const {
        const Real xf = x0_ + (data_.size()-1)*dx_;
        if  (x < x0_ || x > xf) {
            std::ostringstream oss;
            oss.precision(std::numeric_limits<Real>::digits10+3);
            oss << "Requested abscissa x = " << x << ", which is outside of allowed range ["
                << x0_ << ", " << xf << "]";
            throw std::domain_error(oss.str());
        }
        if (x == xf) {
            return data_.back()[1];
        }

        return this->unchecked_prime(x);
    }

    inline Real unchecked_prime(Real x) const {
        using std::floor;
        auto i = static_cast<decltype(data_.size())>(floor((x-x0_)/dx_));
        Real xi = x0_ + i*dx_;
        Real t = (x - xi)/dx_;

        Real v0 = data_[i][1];
        Real v1 = data_[i+1][1];

        return std::numeric_limits<Real>::quiet_NaN();
    }



private:
    RandomAccessContainer data_;
    Real x0_;
    Real dx_;
};

}
#endif 