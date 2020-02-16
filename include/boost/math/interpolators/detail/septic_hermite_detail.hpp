#ifndef BOOST_MATH_INTERPOLATORS_DETAIL_SEPTIC_HERMITE_DETAIL_HPP
#define BOOST_MATH_INTERPOLATORS_DETAIL_SEPTIC_HERMITE_DETAIL_HPP
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <cmath>

namespace boost::math::interpolators::detail {

template<class RandomAccessContainer>
class septic_hermite_detail {
public:
    using Real = typename RandomAccessContainer::value_type;
    septic_hermite_detail(RandomAccessContainer && x, RandomAccessContainer && y, RandomAccessContainer && dydx, RandomAccessContainer && d2ydx2, RandomAccessContainer && d3ydx3) 
    : x_{std::move(x)}, y_{std::move(y)}, dydx_{std::move(dydx)}, d2ydx2_{std::move(d2ydx2)}, d3ydx3_{std::move(d3ydx3)}
    {
        if (x_.size() != y_.size()) {
            throw std::domain_error("Number of abscissas must = number of ordinates.");
        }
        if (x_.size() != dydx_.size()) {
            throw std::domain_error("Numbers of derivatives must = number of abscissas.");
        }
        if (x_.size() != d2ydx2_.size()) {
            throw std::domain_error("Number of second derivatives must equal number of abscissas.");
        }
        if (x_.size() != d3ydx3_.size()) {
            throw std::domain_error("Number of third derivatives must equal number of abscissas.");
        }

        if (x_.size() < 2) {
            throw std::domain_error("At least 2 abscissas are required.");
        }
        Real x0 = x_[0];
        for (decltype(x_.size()) i = 1; i < x_.size(); ++i) {
            Real x1 = x_[i];
            if (x1 <= x0)
            {
                throw std::domain_error("Abscissas must be sorted in strictly increasing order x0 < x1 < ... < x_{n-1}");
            }
            x0 = x1;
        }
    }

    void push_back(Real x, Real y, Real dydx, Real d2ydx2, Real d3ydx3) {
        using std::abs;
        using std::isnan;
        if (x <= x_.back()) {
             throw std::domain_error("Calling push_back must preserve the monotonicity of the x's");
        }
        x_.push_back(x);
        y_.push_back(y);
        dydx_.push_back(dydx);
        d2ydx2_.push_back(d2ydx2);
        d3ydx3_.push_back(d3ydx3);
    }

    Real operator()(Real x) const {
        if  (x < x_[0] || x > x_.back()) {
            std::ostringstream oss;
            oss.precision(std::numeric_limits<Real>::digits10+3);
            oss << "Requested abscissa x = " << x << ", which is outside of allowed range ["
                << x_[0] << ", " << x_.back() << "]";
            throw std::domain_error(oss.str());
        }
        // We need t := (x-x_k)/(x_{k+1}-x_k) \in [0,1) for this to work.
        // Sadly this neccessitates this loathesome check, otherwise we get t = 1 at x = xf.
        if (x == x_.back()) {
            return y_.back();
        }

        auto it = std::upper_bound(x_.begin(), x_.end(), x);
        auto i = std::distance(x_.begin(), it) -1;
        Real x0 = *(it-1);
        Real x1 = *it;
        Real y0 = y_[i];
        Real y1 = y_[i+1];
        // Velocity:
        Real v0 = dydx_[i];
        Real v1 = dydx_[i+1];
        // Acceleration:
        Real a0 = d2ydx2_[i];
        Real a1 = d2ydx2_[i+1];
        // Jerk:
        Real j0 = d3ydx3_[i];
        Real j1 = d3ydx3_[i+1];

        Real dx = (x1-x0);
        Real t = (x-x0)/dx;

        // See: 
        // http://seisweb.usask.ca/classes/GEOL481/2017/Labs/interpolation_utilities_matlab/shermite.m
        Real t2 = t*t;
        Real t3 = t2*t;
        Real t4 = t3*t;
        Real t5 = t4*t;
        Real t6 = t5*t;
        Real t7 = t6*t;
        Real dx2 = dx*dx;
        Real dx3 = dx2*dx;


        Real s = t4*(-35 + t*(84 + t*(-70 + 20*t)));
        Real z4 = -s;
        Real z0 = s + 1;
        Real z1 = 10*t7 - 36*t6 + 45*t5 - 20*t4 + t;
        Real z2 = 2*t7 - 15*t6/2 + 10*t5 - 5*t4 + t2/2;
        Real z3 = t7/6 - 2*t6/3 + t5 - 2*t4/3 + t3/6;
        
        Real z5 = 10*t7 - 34*t6 + 39*t5 - 15*t4;
        Real z6 = -2*t7 + 13*t6/2 - 7*t5 + 5*t4/2;
        Real z7 = t7/6 - t6/2 + t5/2 - t4/6;

        Real y = z0*y0 + z1*dx*v0 + z2*dx2*a0 + z3*dx3*j0 + z4*y1 + z5*dx*v1 + z6*dx2*a1 + z7*dx3*j1;
        return y;
    }

    Real prime(Real x) const {
        if  (x < x_[0] || x > x_.back()) {
            std::ostringstream oss;
            oss.precision(std::numeric_limits<Real>::digits10+3);
            oss << "Requested abscissa x = " << x << ", which is outside of allowed range ["
                << x_[0] << ", " << x_.back() << "]";
            throw std::domain_error(oss.str());
        }
        if (x == x_.back()) {
            return dydx_.back();
        }

        auto it = std::upper_bound(x_.begin(), x_.end(), x);
        auto i = std::distance(x_.begin(), it) -1;
        Real x0 = *(it-1);
        Real x1 = *it;
        Real s0 = dydx_[i];
        Real s1 = dydx_[i+1];

        // Ridiculous linear interpolation. Fine for now:
        Real numerator = s0*(x1-x) + s1*(x-x0);
        Real denominator = x1 - x0;
        return numerator/denominator;
    }


    friend std::ostream& operator<<(std::ostream & os, const septic_hermite_detail & m)
    {
        os << "(x,y,y') = {";
        for (size_t i = 0; i < m.x_.size() - 1; ++i) {
            os << "(" << m.x_[i] << ", " << m.y_[i] << ", " << m.dydx_[i] << ", " << m.d2ydx2_[i] <<  ", " << m.d3ydx3_[i] << "),  ";
        }
        auto n = m.x_.size()-1;
        os << "(" << m.x_[n] << ", " << m.y_[n] << ", " << m.dydx_[n] << ", " << m.d2ydx2_[n] << m.d3ydx3_[n] << ")}";
        return os;
    }



private:
    RandomAccessContainer x_;
    RandomAccessContainer y_;
    RandomAccessContainer dydx_;
    RandomAccessContainer d2ydx2_;
    RandomAccessContainer d3ydx3_;
};

template<class RandomAccessContainer>
class cardinal_septic_hermite_detail {
public:
    using Real = typename RandomAccessContainer::value_type;
    cardinal_septic_hermite_detail(RandomAccessContainer && y, RandomAccessContainer && dydx, RandomAccessContainer && d2ydx2, RandomAccessContainer && d3ydx3, Real x0, Real dx)
    : y_{std::move(y)}, dy_{std::move(dydx)}, d2y_{std::move(d2ydx2)}, d3y_{std::move(d3ydx3)}, x0_{x0}, inv_dx_{1/dx}
    {
        if (y_.size() != dy_.size())
        {
            throw std::domain_error("Numbers of derivatives must = number of ordinates.");
        }
        if (y_.size() != d2y_.size())
        {
            throw std::domain_error("Number of second derivatives must equal number of ordinates.");
        }
        if (y_.size() != d3y_.size())
        {
            throw std::domain_error("Number of third derivatives must equal number of ordinates.");
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
            d2y *= (dx*dx/2);
        }

        for (auto & d3y : d3y_)
        {
            d3y *= (dx*dx*dx/6);
        }

    }

    inline Real operator()(Real x) const
    {
        Real xf = x0_ + (y_.size()-1)/inv_dx_;
        if  (x < x0_ || x > xf)
        {
            std::ostringstream oss;
            oss.precision(std::numeric_limits<Real>::digits10+3);
            oss << "Requested abscissa x = " << x << ", which is outside of allowed range ["
                << x0_ << ", " << xf << "]";
            throw std::domain_error(oss.str());
        }
        // We need t := (x-x_k)/(x_{k+1}-x_k) \in [0,1) for this to work.
        // Sadly this neccessitates this loathesome check, otherwise we get t = 1 at x = xf.
        if (x == xf)
        {
            return y_.back();
        }
        return this->unchecked_evaluation(x);
    }

    inline Real unchecked_evaluation(Real x) const {
        using std::floor;
        Real s3 = (x-x0_)*inv_dx_;
        Real ii = floor(s3);
        auto i = static_cast<decltype(y_.size())>(ii);
        Real t = s3 - ii;

        Real y0 = y_[i];
        Real y1 = y_[i+1];
        // Velocity:
        Real dy0 = dy_[i];
        Real dy1 = dy_[i+1];
        // Acceleration:
        Real a0 = d2y_[i];
        Real a1 = d2y_[i+1];
        // Jerk:
        Real j0 = d3y_[i];
        Real j1 = d3y_[i+1];

        // See:
        // http://seisweb.usask.ca/classes/GEOL481/2017/Labs/interpolation_utilities_matlab/shermite.m
        Real t2 = t*t;
        Real t3 = t2*t;
        Real t4 = t3*t;

        Real s = t4*(-35 + t*(84 + t*(-70 + 20*t)));
        Real z4 = -s;
        Real z0 = s + 1;
        //Real z1 = t*(10*t6 - 36*t5 + 45*t4 - 20*t3 + 1);
        Real z1 = t*(1 + t3*(-20 + t*(45 + t*(-36+10*t))));
        //Real z2 = t2*(4*t5 - 15*t4 + 20*t3 - 10*t2 + 1);
        Real z2 = t2*(1+ t2*(4*t3 - 15*t + 20*t - 10));
        //Real z3 = t3*(t4 - 4*t3 + 6*t2 - 4*t + 1);
        Real z3 = t3*(1 + t*(-4+t*(6+t*(-4+t))));

        //Real z5 = t4*(10*t3 - 34*t2 + 39*t - 15);
        Real z5 = t4*(-15 + t*(39 + t*(-34 + 10*t)));
        //Real z6 = t4*(-4*t3 + 13*t2 - 14*t + 5);
        Real z6 = t4*(5 + t*(-14 + t*(13-4*t)));
        //Real z7 = t4*(t3 - 3*t2 + 3*t - 1);
        Real z7 = t4*(-1 + t*(3+t*(-3+t)));

        return z0*y0 + z1*dy0 + z2*a0 + z3*j0 + z4*y1 + z5*dy1 + z6*a1 + z7*j1;
    }

    inline Real prime(Real x) const {
        Real xf = x0_ + (y_.size()-1)/inv_dx_;
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
            return dy_.back()/inv_dx_;
        }

        return this->unchecked_prime(x);
    }

    inline Real unchecked_prime(Real x) const {
        //TODO: Get the high accuracy approximation by differentiating the interpolant!
        using std::floor;
        Real s3 = (x-x0_)*inv_dx_;
        Real ii = floor(s3);
        auto i = static_cast<decltype(y_.size())>(ii);
        Real t = s3 - ii;

        // Velocity:
        Real v0 = dy_[i];
        Real v1 = dy_[i+1];
        return (v0+v1)*inv_dx_/2;
    }

private:
    RandomAccessContainer y_;
    RandomAccessContainer dy_;
    RandomAccessContainer d2y_;
    RandomAccessContainer d3y_;
    Real x0_;
    Real inv_dx_;
};


template<class RandomAccessContainer>
class cardinal_septic_hermite_detail_aos {
public:
    using Point = typename RandomAccessContainer::value_type;
    using Real = typename Point::value_type;
    cardinal_septic_hermite_detail_aos(RandomAccessContainer && dat, Real x0, Real dx)
    : data_{std::move(dat)}, x0_{x0}, inv_dx_{1/dx}
    {
        if (data_.size() < 2) {
            throw std::domain_error("At least 2 abscissas are required.");
        }
        if (data_[0].size() != 4) {
            throw std::domain_error("There must be 4 data items per struct.");
        }

        for (auto & datum : data_)
        {
            datum[1] *= dx;
            datum[2] *= (dx*dx/2);
            datum[3] *= (dx*dx*dx/6);
        }
    }

    inline Real operator()(Real x) const
    {
        Real xf = x0_ + (data_.size()-1)/inv_dx_;
        if  (x < x0_ || x > xf)
        {
            std::ostringstream oss;
            oss.precision(std::numeric_limits<Real>::digits10+3);
            oss << "Requested abscissa x = " << x << ", which is outside of allowed range ["
                << x0_ << ", " << xf << "]";
            throw std::domain_error(oss.str());
        }
        // We need t := (x-x_k)/(x_{k+1}-x_k) \in [0,1) for this to work.
        // Sadly this neccessitates this loathesome check, otherwise we get t = 1 at x = xf.
        if (x == xf)
        {
            return data_.back()[0];
        }
        return this->unchecked_evaluation(x);
    }

    inline Real unchecked_evaluation(Real x) const {
        using std::floor;
        Real s3 = (x-x0_)*inv_dx_;
        Real ii = floor(s3);
        auto i = static_cast<decltype(data_.size())>(ii);
        Real t = s3 - ii;


        Real y0 = data_[i][0];
        Real y1 = data_[i+1][0];
        // Velocity:
        Real dy0 = data_[i][1];
        Real dy1 = data_[i+1][1];
        // Acceleration:
        Real a0 = data_[i][2];
        Real a1 = data_[i+1][2];
        // Jerk:
        Real j0 = data_[i][3];
        Real j1 = data_[i+1][3];

        Real t2 = t*t;
        Real t3 = t2*t;
        Real t4 = t3*t;

        Real s = t4*(-35 + t*(84 + t*(-70 + 20*t)));
        Real z4 = -s;
        Real z0 = s + 1;
        //Real z1 = t*(10*t6 - 36*t5 + 45*t4 - 20*t3 + 1);
        Real z1 = t*(1 + t3*(-20 + t*(45 + t*(-36+10*t))));
        //Real z2 = t2*(4*t5 - 15*t4 + 20*t3 - 10*t2 + 1);
        Real z2 = t2*(1+ t2*(4*t3 - 15*t + 20*t - 10));
        //Real z3 = t3*(t4 - 4*t3 + 6*t2 - 4*t + 1);
        Real z3 = t3*(1 + t*(-4+t*(6+t*(-4+t))));

        //Real z5 = t4*(10*t3 - 34*t2 + 39*t - 15);
        Real z5 = t4*(-15 + t*(39 + t*(-34 + 10*t)));
        //Real z6 = t4*(-4*t3 + 13*t2 - 14*t + 5);
        Real z6 = t4*(5 + t*(-14 + t*(13-4*t)));
        //Real z7 = t4*(t3 - 3*t2 + 3*t - 1);
        Real z7 = t4*(-1 + t*(3+t*(-3+t)));

        return z0*y0 + z1*dy0 + z2*a0 + z3*j0 + z4*y1 + z5*dy1 + z6*a1 + z7*j1;

    }

    inline Real prime(Real x) const {
        Real xf = x0_ + (data_.size()-1)/inv_dx_;
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
            return data_.back()[1];
        }

        return this->unchecked_prime(x);
    }

    inline Real unchecked_prime(Real x) const {
        return std::numeric_limits<Real>::quiet_NaN();
    }

private:
    RandomAccessContainer data_;
    Real x0_;
    Real inv_dx_;
};


}
#endif
