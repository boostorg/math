#ifndef BOOST_MATH_INTERPOLATORS_DETAIL_CARDINAL_CUBIC_HERMITE_DETAIL_HPP
#define BOOST_MATH_INTERPOLATORS_DETAIL_CARDINAL_CUBIC_HERMITE_DETAIL_HPP
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <limits>

namespace boost::math::interpolators::detail {

template<class RandomAccessContainer>
class cardinal_cubic_hermite_detail {
public:
    using Real = typename RandomAccessContainer::value_type;

    cardinal_cubic_hermite_detail(RandomAccessContainer && y, RandomAccessContainer dydx, Real x0, Real dx)
    : y_{std::move(y)}, dy_{std::move(dydx)}, x0_{x0}, inv_dx_{1/dx}
    {
        using std::abs;
        using std::isnan;
        if (y_.size() != dy_.size())
        {
            throw std::domain_error("There must be the same number of derivatives as ordinates.");
        }
        if (y_.size() < 2)
        {
            throw std::domain_error("Must be at least two data points.");
        }
        if (dx <= 0)
        {
            throw std::domain_error("dx > 0 is required.");
        }

        for (auto & dy : dy_)
        {
            dy *= dx;
        }
    }

    // Why not implement push_back? It's awkward: If the buffer is circular, x0_ += dx_.
    // If the buffer is not circular, x0_ is unchanged.
    // We need a concept for circular_buffer!

    inline Real operator()(Real x) const {
        const Real xf = x0_ + (y_.size()-1)/inv_dx_;
        if  (x < x0_ || x > xf)
        {
            std::ostringstream oss;
            oss.precision(std::numeric_limits<Real>::digits10+3);
            oss << "Requested abscissa x = " << x << ", which is outside of allowed range ["
                << x0_ << ", " << xf << "]";
            throw std::domain_error(oss.str());
        }
        if (x == xf)
        {
            return y_.back();
        }
        return this->unchecked_evaluation(x);
    }

    inline Real unchecked_evaluation(Real x) const
    {
        using std::floor;
        Real s = (x-x0_)*inv_dx_;
        Real ii = floor(s);
        auto i = static_cast<decltype(y_.size())>(ii);
        Real t = s - ii;
        Real y0 = y_[i];
        Real y1 = y_[i+1];
        Real dy0 = dy_[i];
        Real dy1 = dy_[i+1];

        Real r = 1-t;
        return r*r*(y0*(1+2*t) + dy0*t)
              + t*t*(y1*(3-2*t) - dy1*r);
    }

    inline Real prime(Real x) const {
        const Real xf = x0_ + (y_.size()-1)/inv_dx_;
        if  (x < x0_ || x > xf)
        {
            std::ostringstream oss;
            oss.precision(std::numeric_limits<Real>::digits10+3);
            oss << "Requested abscissa x = " << x << ", which is outside of allowed range ["
                << x0_ << ", " << xf << "]";
            throw std::domain_error(oss.str());
        }
        if (x == xf)
        {
            return dy_.back()*inv_dx_;
        }
        return this->unchecked_prime(x);
    }

    inline Real unchecked_prime(Real x) const {
        using std::floor;
        Real s = (x-x0_)*inv_dx_;
        Real ii = floor(s);
        auto i = static_cast<decltype(y_.size())>(ii);
        Real t = s - ii;
        Real y0 = y_[i];
        Real y1 = y_[i+1];
        Real dy0 = dy_[i];
        Real dy1 = dy_[i+1];

        Real dy = 6*t*(1-t)*(y1 - y0)  + (3*t*t - 4*t+1)*dy0 + t*(3*t-2)*dy1;
        return dy*inv_dx_;
    }


    auto size() const
    {
        return y_.size();
    }

    RandomAccessContainer y_;
    RandomAccessContainer dy_;
    Real x0_;
    Real inv_dx_;
};


template<class RandomAccessContainer>
class cardinal_cubic_hermite_detail_aos {
public:
    using Point = typename RandomAccessContainer::value_type;
    using Real = typename Point::value_type;

    cardinal_cubic_hermite_detail_aos(RandomAccessContainer && dat, Real x0, Real dx)
    : dat_{std::move(dat)}, x0_{x0}, inv_dx_{1/dx}
    {
        if (dat_.size() < 2)
        {
            throw std::domain_error("Must be at least two data points.");
        }
        if (dat_[0].size() != 2)
        {
            throw std::domain_error("Each datum must contain (y, y'), and nothing else.");
        }
        if (dx <= 0)
        {
            throw std::domain_error("dx > 0 is required.");
        }

        for (auto & d : dat_)
        {
            d[1] *= dx;
        }
    }

    inline Real operator()(Real x) const {
        const Real xf = x0_ + (dat_.size()-1)/inv_dx_;
        if  (x < x0_ || x > xf)
        {
            std::ostringstream oss;
            oss.precision(std::numeric_limits<Real>::digits10+3);
            oss << "Requested abscissa x = " << x << ", which is outside of allowed range ["
                << x0_ << ", " << xf << "]";
            throw std::domain_error(oss.str());
        }
        if (x == xf)
        {
            return dat_.back()[0];
        }
        return this->unchecked_evaluation(x);
    }

    inline Real unchecked_evaluation(Real x) const
    {
        using std::floor;
        Real s = (x-x0_)*inv_dx_;
        Real ii = floor(s);
        auto i = static_cast<decltype(dat_.size())>(ii);

        Real t = s - ii;
        Real y0 = dat_[i][0];
        Real y1 = dat_[i+1][0];
        Real dy0 = dat_[i][1];
        Real dy1 = dat_[i+1][1];

        Real r = 1-t;
        return r*r*(y0*(1+2*t) + dy0*t)
              + t*t*(y1*(3-2*t) - dy1*r);
    }

    inline Real prime(Real x) const {
        const Real xf = x0_ + (dat_.size()-1)/inv_dx_;
        if  (x < x0_ || x > xf)
        {
            std::ostringstream oss;
            oss.precision(std::numeric_limits<Real>::digits10+3);
            oss << "Requested abscissa x = " << x << ", which is outside of allowed range ["
                << x0_ << ", " << xf << "]";
            throw std::domain_error(oss.str());
        }
        if (x == xf)
        {
            return dat_.back()[1];
        }
        return this->unchecked_prime(x);
    }

    inline Real unchecked_prime(Real x) const
    {
        using std::floor;
        Real s = (x-x0_)*inv_dx_;
        Real ii = floor(s);
        auto i = static_cast<decltype(dat_.size())>(ii);
        Real t = s - ii;
        Real y0 = dat_[i][0];
        Real dy0 = dat_[i][1];
        Real y1 = dat_[i+1][0];
        Real dy1 = dat_[i+1][1];

        Real dy = 6*t*(1-t)*(y1 - y0)  + (3*t*t - 4*t+1)*dy0 + t*(3*t-2)*dy1;
        return dy*inv_dx_;
    }


    auto size() const
    {
        return dat_.size();
    }

    RandomAccessContainer dat_;
    Real x0_;
    Real inv_dx_;
};


}
#endif