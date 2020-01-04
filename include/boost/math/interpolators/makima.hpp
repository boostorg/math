// Copyright Nick Thompson, 2019
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// See: https://blogs.mathworks.com/cleve/2019/04/29/makima-piecewise-cubic-interpolation/
// And: https://doi.org/10.1145/321607.321609

#ifndef BOOST_MATH_INTERPOLATORS_MAKIMA_HPP
#define BOOST_MATH_INTERPOLATORS_MAKIMA_HPP
#include <stdexcept>
#include <algorithm>

namespace boost::math::interpolators {

template<class RandomAccessContainer>
class makima {
public:
    using Real = typename RandomAccessContainer::value_type;

    makima(RandomAccessContainer && x, RandomAccessContainer && y) : x_{std::move(x)}, y_{std::move(y)}
    {
        using std::abs;
        if (x_.size() != y_.size())
        {
            throw std::domain_error("There must be the same number of ordinates as abscissas.");
        }
        if (x_.size() < 2)
        {
            throw std::domain_error("Must be at least two data points.");
        }
        Real x0 = x_[0];
        for (size_t i = 1; i < x_.size(); ++i) {
            Real x1 = x_[i];
            if (x1 <= x0) {
                throw std::domain_error("Abscissas must be listed in strictly increasing order x0 < x1 < ... < x_{n-1}");
            }
            x0 = x1;
        }

        s_.resize(x_.size());
        for (decltype(s_.size()) i = 2; i < s_.size()-1; ++i) {
            Real mim2 = (y_[i-1]-y_[i-2])/(x_[i-1]-x_[i-2]);
            Real mim1 = (y_[i]-y_[i-1])/(x_[i]-x_[i-1]);
            Real mi = (y_[i+1]-y_[i])/(x_[i+1]-x_[i]);
            Real mip1 = (y_[i+2]-y_[i+1])/(x_[i+2]-x_[i+1]);
            Real numerator = (abs(mip1-mi)*mim1 + abs(mim1-mim2)*mi);
            Real denominator = abs(mip1-mi) + abs(mim1-mim2);
            s_[i] = numerator/denominator;
        }
    }

    Real operator()(Real x) const {
        if  (x < x_[0] || x > x_[x_.size()-1]) {
            std::string err = "Requested abscissa outside of allowed range [" 
                             + std::to_string(x_[0]) + ", "  + std::to_string(x_[x_.size()-1]) + "]";
            throw std::domain_error(err);
        }

        auto it = std::lower_bound(x_.begin(), x_.end(), x);
        auto i = std::distance(x_.begin(), it);
        Real x1 = *it;
        Real x2 = *(it+1);
        Real y1 = y_[i];
        Real y2 = y_[i+1];
        Real s1 = s_[i];
        Real s2 = s_[i+1];
        return std::numeric_limits<Real>::quiet_NaN();
    }    

private:
    RandomAccessContainer x_;
    RandomAccessContainer y_;
    RandomAccessContainer s_;
};
}
#endif