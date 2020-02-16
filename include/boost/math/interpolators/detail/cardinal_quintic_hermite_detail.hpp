#ifndef BOOST_MATH_INTERPOLATORS_DETAIL_CARDINAL_QUINTIC_HERMITE_DETAIL_HPP
#define BOOST_MATH_INTERPOLATORS_DETAIL_CARDINAL_QUINTIC_HERMITE_DETAIL_HPP
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
    : y_{std::move(y)}, dy_{std::move(dydx)}, d2y_{std::move(d2ydx2)}, x0_{x0}, inv_dx_{1/dx}
    {
        if (y_.size() != dy_.size())
        {
            throw std::domain_error("Numbers of derivatives must = number of abscissas.");
        }
        if (y_.size() != d2y_.size())
        {
            throw std::domain_error("Number of second derivatives must equal number of abscissas.");
        }
        if (y_.size() < 2)
        {
            throw std::domain_error("At least 2 abscissas are required.");
        }
        if (dx <= 0)
        {
            throw std::domain_error("dx > 0 is required.");
        }

        for (auto & dy : dy_)
        {
            dy *= dx;
        }

        for (auto & d2y : d2y_)
        {
            d2y *= (dx*dx)/2;
        }
    }


    inline Real operator()(Real x) const
    {
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
        Real d2y0 = d2y_[i];
        Real d2y1 = d2y_[i+1];

        // See the 'Basis functions' section of:
        // https://www.rose-hulman.edu/~finn/CCLI/Notes/day09.pdf
        // Also: https://github.com/MrHexxx/QuinticHermiteSpline/blob/master/HermiteSpline.cs
        Real y = (1- t*t*t*(10 + t*(-15 + 6*t)))*y0;
        y += t*(1+ t*t*(-6 + t*(8 -3*t)))*dy0;
        y += t*t*(1 + t*(-3 + t*(3-t)))*d2y0;
        y += t*t*t*((1 + t*(-2 + t))*d2y1 + (-4 + t*(7 -3*t))*dy1 + (10 + t*(-15 + 6*t))*y1);
        return y;
    }

    inline Real prime(Real x) const
    {
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

    inline Real unchecked_prime(Real x) const
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
        Real d2y0 = d2y_[i];
        Real d2y1 = d2y_[i+1];

        Real dydx = 30*t*t*(1 - 2*t + t*t)*(y1-y0);
        dydx += (1-18*t*t + 32*t*t*t - 15*t*t*t*t)*dy0 - t*t*(12 - 28*t + 15*t*t)*dy1;
        dydx += t*((2 - 9*t + 12*t*t - 5*t*t*t)*d2y0 + t*(3 - 8*t + 5*t*t)*d2y1);
        dydx *= inv_dx_;
        return dydx;
    }

    inline Real double_prime(Real x) const {
        const Real xf = x0_ + (y_.size()-1)/inv_dx_;
        if  (x < x0_ || x > xf) {
            std::ostringstream oss;
            oss.precision(std::numeric_limits<Real>::digits10+3);
            oss << "Requested abscissa x = " << x << ", which is outside of allowed range ["
                << x0_ << ", " << xf << "]";
            throw std::domain_error(oss.str());
        }
        if (x == xf)
        {
            return d2y_.back()*2*inv_dx_*inv_dx_;
        }

        return this->unchecked_double_prime(x);
    }

    inline Real unchecked_double_prime(Real x) const {
        using std::floor;
        Real s = (x-x0_)*inv_dx_;
        Real ii = floor(s);
        auto i = static_cast<decltype(y_.size())>(ii);
        Real t = s - ii;


        Real y0 = y_[i];
        Real y1 = y_[i+1];
        Real dy0 = dy_[i];
        Real dy1 = dy_[i+1];
        Real d2y0 = d2y_[i];
        Real d2y1 = d2y_[i+1];

        Real d2ydx2 = 60*t*(1 - 3*t + 2*t*t)*(y1 - y0)*inv_dx_*inv_dx_;
        d2ydx2 += (12*t)*((-3 + 8*t - 5*t*t)*dy0 - (2 - 7*t + 5*t*t)*dy1);
        d2ydx2 += (1 - 9*t + 18*t*t - 10*t*t*t)*d2y0*(2*inv_dx_*inv_dx_) + t*(3 - 12*t + 10*t*t)*d2y1*(2*inv_dx_*inv_dx_);
        return d2ydx2;
    }

private:
    RandomAccessContainer y_;
    RandomAccessContainer dy_;
    RandomAccessContainer d2y_;
    Real x0_;
    Real inv_dx_;
};


template<class RandomAccessContainer>
class cardinal_quintic_hermite_detail_aos {
public:
    using Point = typename RandomAccessContainer::value_type;
    using Real = typename Point::value_type;
    cardinal_quintic_hermite_detail_aos(RandomAccessContainer && data, Real x0, Real dx)
    : data_{std::move(data)} , x0_{x0}, inv_dx_{1/dx}
    {
        if (data_.size() < 2)
        {
            throw std::domain_error("At least two points are required to interpolate using cardinal_quintic_hermite_aos");
        }

        if (data_[0].size() != 3)
        {
            throw std::domain_error("Each datum passed to the cardinal_quintic_hermite_aos must have three elements: {y, y', y''}");
        }
        if (dx <= 0)
        {
            throw std::domain_error("dx > 0 is required.");
        }

        for (auto & datum : data_)
        {
            datum[1] *= dx;
            datum[2] *= (dx*dx/2);
        }
    }


    inline Real operator()(Real x) const
    {
        const Real xf = x0_ + (data_.size()-1)/inv_dx_;
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
            return data_.back()[0];
        }
        return this->unchecked_evaluation(x);
    }

    inline Real unchecked_evaluation(Real x) const
    {
        using std::floor;
        Real s = (x-x0_)*inv_dx_;
        Real ii = floor(s);
        auto i = static_cast<decltype(data_.size())>(ii);
        Real t = s - ii;


        Real y0 = data_[i][0];
        Real dy0 = data_[i][1];
        Real d2y0 = data_[i][2];
        Real y1 = data_[i+1][0];
        Real dy1 = data_[i+1][1];
        Real d2y1 = data_[i+1][2];

        Real y = (1- t*t*t*(10 + t*(-15 + 6*t)))*y0;
        y += t*(1+ t*t*(-6 + t*(8 -3*t)))*dy0;
        y += t*t*(1 + t*(-3 + t*(3-t)))*d2y0;
        y += t*t*t*((1 + t*(-2 + t))*d2y1 + (-4 + t*(7 -3*t))*dy1 + (10 + t*(-15 + 6*t))*y1);
        return y;
    }

    inline Real prime(Real x) const
    {
        const Real xf = x0_ + (data_.size()-1)/inv_dx_;
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
            return data_.back()[1]*inv_dx_;
        }

        return this->unchecked_prime(x);
    }

    inline Real unchecked_prime(Real x) const {
        using std::floor;
        Real s = (x-x0_)*inv_dx_;
        Real ii = floor(s);
        auto i = static_cast<decltype(data_.size())>(ii);
        Real t = s - ii;


        Real y0 = data_[i][0];
        Real y1 = data_[i+1][0];
        Real v0 = data_[i][1];
        Real v1 = data_[i+1][1];
        Real a0 = data_[i][2];
        Real a1 = data_[i+1][2];

        Real dy = 30*t*t*(1 - 2*t + t*t)*(y1-y0);
        dy += (1-18*t*t + 32*t*t*t - 15*t*t*t*t)*v0 - t*t*(12 - 28*t + 15*t*t)*v1;
        dy += t*((2 - 9*t + 12*t*t - 5*t*t*t)*a0 + t*(3 - 8*t + 5*t*t)*a1);
        return dy*inv_dx_;
    }



private:
    RandomAccessContainer data_;
    Real x0_;
    Real inv_dx_;
};

}
#endif
