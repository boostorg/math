// Copyright Nick Thompson, 2020
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// See Fritsch and Carlson: https://doi.org/10.1137/0717021
#ifndef BOOST_MATH_INTERPOLATORS_DETAIL_PCHIP_DETAIL_HPP
#define BOOST_MATH_INTERPOLATORS_DETAIL_PCHIP_DETAIL_HPP
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <limits>

namespace boost::math::interpolators::detail {

template<class RandomAccessContainer>
class pchip_detail {
public:
    using Real = typename RandomAccessContainer::value_type;

    pchip_detail(RandomAccessContainer && x, RandomAccessContainer && y,
                 Real left_endpoint_derivative = std::numeric_limits<Real>::quiet_NaN(),
                 Real right_endpoint_derivative = std::numeric_limits<Real>::quiet_NaN()) : x_{std::move(x)}, y_{std::move(y)}
    {
        using std::abs;
        using std::isnan;
        if (x_.size() != y_.size())
        {
            throw std::domain_error("There must be the same number of ordinates as abscissas.");
        }
        if (x_.size() < 4)
        {
            throw std::domain_error("Must be at least four data points.");
        }
        Real x0 = x_[0];
        for (size_t i = 1; i < x_.size(); ++i) {
            Real x1 = x_[i];
            if (x1 <= x0) {
                throw std::domain_error("Abscissas must be listed in strictly increasing order x0 < x1 < ... < x_{n-1}");
            }
            x0 = x1;
        }

        s_.resize(x_.size(), std::numeric_limits<Real>::quiet_NaN());
        if (isnan(left_endpoint_derivative))
        {
            // O(h) finite difference derivative:
            s_[0] = (y_[1]-y_[0])/(x_[1]-x_[0]);
        }
        else
        {
            s_[0] = left_endpoint_derivative;
        }

        for (decltype(s_.size()) k = 1; k < s_.size()-1; ++k) {
            Real hkm1 = x_[k] - x_[k-1];
            Real dkm1 = (y_[k] - y_[k-1])/hkm1;

            Real hk = x_[k+1] - x_[k];
            Real dk = (y_[k+1] - y_[k])/hk;
            Real w1 = 2*hk + hkm1;
            Real w2 = hk + 2*hkm1;
            if ( (dk > 0 && dkm1 < 0) || (dk < 0 && dkm1 > 0) || dk == 0 || dkm1 == 0)
            {
                s_[k] = 0;
            }
            else
            {
                s_[k] = (w1+w2)/(w1/dkm1 + w2/dk);
            }

        }
        // Quadratic extrapolation at the other end:
        
        decltype(s_.size()) n = s_.size();
        if (isnan(right_endpoint_derivative))
        {
            s_[n-1] = (y_[n-1]-y_[n-2])/(x_[n-1] - x_[n-2]);
        }
        else
        {
            s_[n-1] = right_endpoint_derivative;
        }
    }

    void push_back(Real x, Real y) {
        using std::abs;
        using std::isnan;
        if (x <= x_.back()) {
             throw std::domain_error("Calling push_back must preserve the monotonicity of the x's");
        }
        x_.push_back(x);
        y_.push_back(y);
        s_.push_back(std::numeric_limits<Real>::quiet_NaN());
        decltype(s_.size()) n = s_.size();
        s_[n-1] = (y_[n-1]-y_[n-2])/(x_[n-1] - x_[n-2]);
        // Now fix s_[n-2]:
        auto k = n-2;
        Real hkm1 = x_[k] - x_[k-1];
        Real dkm1 = (y_[k] - y_[k-1])/hkm1;

        Real hk = x_[k+1] - x_[k];
            Real dk = (y_[k+1] - y_[k])/hk;
            Real w1 = 2*hk + hkm1;
            Real w2 = hk + 2*hkm1;
            if ( (dk > 0 && dkm1 < 0) || (dk < 0 && dkm1 > 0) || dk == 0 || dkm1 == 0)
            {
                s_[k] = 0;
            }
            else
            {
                s_[k] = (w1+w2)/(w1/dkm1 + w2/dk);
            }

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
        Real s0 = s_[i];
        Real s1 = s_[i+1];
        Real dx = (x1-x0);
        Real t = (x-x0)/dx;

        // See the section 'Representations' in the page
        // https://en.wikipedia.org/wiki/Cubic_Hermite_spline
        // This uses the factorized form:
        //Real y = y0*(1+2*t)*(1-t)*(1-t) + dx*s0*t*(1-t)*(1-t)
        //       + y1*t*t*(3-2*t) + dx*s1*t*t*(t-1);
        // And then factorized further:
        Real y = (1-t)*(1-t)*(y0*(1+2*t) + s0*(x-x0))
              + t*t*(y1*(3-2*t) + dx*s1*(t-1));
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
            return s_.back();
        }

        auto it = std::upper_bound(x_.begin(), x_.end(), x);
        auto i = std::distance(x_.begin(), it) -1;
        Real x0 = *(it-1);
        Real x1 = *it;
        Real s0 = s_[i];
        Real s1 = s_[i+1];

        // Ridiculous linear interpolation. Fine for now:
        Real numerator = s0*(x1-x) + s1*(x-x0);
        Real denominator = x1 - x0;
        return numerator/denominator;
    }


    friend std::ostream& operator<<(std::ostream & os, const pchip_detail & m)
    {
        os << "(x,y,y') = {";
        for (size_t i = 0; i < m.x_.size() - 1; ++i) {
            os << "(" << m.x_[i] << ", " << m.y_[i] << ", " << m.s_[i] << "),  ";
        }
        auto n = m.x_.size()-1;
        os << "(" << m.x_[n] << ", " << m.y_[n] << ", " << m.s_[n] << ")}";
        return os;
    }

private:
    RandomAccessContainer x_;
    RandomAccessContainer y_;
    RandomAccessContainer s_;
};
}
#endif
