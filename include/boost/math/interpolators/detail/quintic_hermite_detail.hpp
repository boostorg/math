#ifndef BOOST_MATH_INTERPOLATORS_DETAIL_QUINTIC_HERMITE_DETAIL_HPP
#define BOOST_MATH_INTERPOLATORS_DETAIL_QUINTIC_HERMITE_DETAIL_HPP
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <cmath>

namespace boost::math::interpolators::detail {

template<class RandomAccessContainer>
class quintic_hermite_detail {
public:
    using Real = typename RandomAccessContainer::value_type;
    quintic_hermite_detail(RandomAccessContainer && x, RandomAccessContainer && y, RandomAccessContainer && dydx, RandomAccessContainer && d2ydx2) : x_{std::move(x)}, y_{std::move(y)}, dydx_{std::move(dydx)}, d2ydx2_{std::move(d2ydx2)}
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

    void push_back(Real x, Real y, Real dydx, Real d2ydx2) {
        using std::abs;
        using std::isnan;
        if (x <= x_.back()) {
             throw std::domain_error("Calling push_back must preserve the monotonicity of the x's");
        }
        x_.push_back(x);
        y_.push_back(y);
        dydx_.push_back(dydx);
        d2ydx2_.push_back(d2ydx2);
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
        Real v0 = dydx_[i];
        Real v1 = dydx_[i+1];
        Real a0 = d2ydx2_[i];
        Real a1 = d2ydx2_[i+1];

        Real dx = (x1-x0);
        Real t = (x-x0)/dx;

        // See the 'Basis functions' section of:
        // https://www.rose-hulman.edu/~finn/CCLI/Notes/day09.pdf
        // Also: https://github.com/MrHexxx/QuinticHermiteSpline/blob/master/HermiteSpline.cs
        Real y = (1- t*t*t*(10 + t*(-15 + 6*t)))*y0;
        y += t*(1+ t*t*(-6 + t*(8 -3*t)))*v0*dx;
        y += t*t*(1 + t*(-3 + t*(3-t)))*a0*dx*dx/2;
        y += t*t*t*((1 + t*(-2 + t))*a1*dx*dx/2 + (-4 + t*(7 -3*t))*v1*dx + (10 + t*(-15 + 6*t))*y1);
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


    friend std::ostream& operator<<(std::ostream & os, const quintic_hermite_detail & m)
    {
        os << "(x,y,y') = {";
        for (size_t i = 0; i < m.x_.size() - 1; ++i) {
            os << "(" << m.x_[i] << ", " << m.y_[i] << ", " << m.dydx_[i] << ", " << m.d2ydx2_[i] << "),  ";
        }
        auto n = m.x_.size()-1;
        os << "(" << m.x_[n] << ", " << m.y_[n] << ", " << m.dydx_[n] << ", " << m.d2ydx2_[n] << ")}";
        return os;
    }



private:
    RandomAccessContainer x_;
    RandomAccessContainer y_;
    RandomAccessContainer dydx_;
    RandomAccessContainer d2ydx2_;
};
}
#endif 