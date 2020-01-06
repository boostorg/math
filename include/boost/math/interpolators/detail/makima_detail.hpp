// Copyright Nick Thompson, 2020
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// See: https://blogs.mathworks.com/cleve/2019/04/29/makima-piecewise-cubic-interpolation/
// And: https://doi.org/10.1145/321607.321609

#ifndef BOOST_MATH_INTERPOLATORS_DETAIL_MAKIMA_DETAIL_HPP
#define BOOST_MATH_INTERPOLATORS_DETAIL_MAKIMA_DETAIL_HPP
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>

namespace boost::math::interpolators::detail {

template<class RandomAccessContainer>
class makima_detail {
public:
    using Real = typename RandomAccessContainer::value_type;

    makima_detail(RandomAccessContainer && x, RandomAccessContainer && y) : x_{std::move(x)}, y_{std::move(y)}
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
        Real m2 = (y_[3]-y_[2])/(x_[3]-x_[2]);
        Real m1 = (y_[2]-y_[1])/(x_[2]-x_[1]);
        Real m0 = (y_[1]-y_[0])/(x_[1]-x_[0]);
        // Quadratic extrapolation: m_{-1} = 2m_0 - m_1:
        Real mm1 = 2*m0 - m1;
        // Quadratic extrapolation: m_{-2} = 2*m_{-1}-m_0:
        Real mm2 = 2*mm1 - m0;
        Real w1 = abs(m1-m0) + abs(m1+m0)/2;
        Real w2 = abs(mm1-mm2) + abs(mm1+mm2)/2;
        s_[0] = (w1*mm1 + w2*m0)/(w1+w2);
        if (isnan(s_[0])) {
            s_[0] = 0;
        }

        w1 = abs(m2-m1) + abs(m2+m1)/2;
        w2 = abs(m0-mm1) + abs(m0+mm1)/2;
        s_[1] = (w1*m0 + w2*m1)/(w1+w2);
        if (isnan(s_[1])) {
            s_[1] = 0;
        }

        for (decltype(s_.size()) i = 2; i < s_.size()-2; ++i) {
            Real mim2 = (y_[i-1]-y_[i-2])/(x_[i-1]-x_[i-2]);
            Real mim1 = (y_[i  ]-y_[i-1])/(x_[i  ]-x_[i-1]);
            Real mi   = (y_[i+1]-y_[i  ])/(x_[i+1]-x_[i  ]);
            Real mip1 = (y_[i+2]-y_[i+1])/(x_[i+2]-x_[i+1]);
            w1 = abs(mip1-mi) + abs(mip1+mi)/2;
            w2 = abs(mim1-mim2) + abs(mim1+mim2)/2;
            s_[i] = (w1*mim1 + w2*mi)/(w1+w2);
            if (isnan(s_[i])) {
                s_[i] = 0;
            }
        }
        // Quadratic extrapolation at the other end:
        
        decltype(s_.size()) n = s_.size();
        Real mnm4 = (y_[n-3]-y_[n-4])/(x_[n-3]-x_[n-4]);
        Real mnm3 = (y_[n-2]-y_[n-3])/(x_[n-2]-x_[n-3]);
        Real mnm2 = (y_[n-1]-y_[n-2])/(x_[n-1]-x_[n-2]);
        Real mnm1 = 2*mnm2 - mnm3;
        Real mn = 2*mnm1 - mnm2;
        w1 = abs(mnm1 - mnm2) + abs(mnm1+mnm2)/2;
        w2 = abs(mnm3 - mnm4) + abs(mnm3+mnm4)/2;

        s_[n-2] = (w1*mnm3 + w2*mnm2)/(w1 + w2);
        if (isnan(s_[n-2])) {
            s_[n-2] = 0;
        }

        w1 = abs(mn - mnm1) + abs(mn+mnm1)/2;
        w2 = abs(mnm2 - mnm3) + abs(mnm2+mnm3)/2;

        s_[n-1] = (w1*mnm2 + w2*mnm1)/(w1+w2);
        if (isnan(s_[n-1])) {
            s_[n-1] = 0;
        }

    }

    Real operator()(Real x) const {
        if  (x < x_[0] || x > x_.back()) {
            std::ostringstream oss;
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

    friend std::ostream& operator<<(std::ostream & os, const makima_detail & m)
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