// Copyright Nick Thompson, 2021
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_INTERPOLATORS_BEZIER_POLYNOMIAL_DETAIL_HPP
#define BOOST_MATH_INTERPOLATORS_BEZIER_POLYNOMIAL_DETAIL_HPP
#include <stdexcept>
#include <iostream>
#include <string>

namespace boost::math::interpolators::detail {

template <class RandomAccessContainer>
class bezier_polynomial_imp
{
public:
    using Point = typename RandomAccessContainer::value_type;
    using Real = typename Point::value_type;
    using Z = typename RandomAccessContainer::size_type;

    bezier_polynomial_imp(RandomAccessContainer && control_points)
    {
        using std::to_string;
        if (control_points.size() < 2) {
            std::string err = std::string(__FILE__) + ":" + to_string(__LINE__)
               + " At least two points are required to form a Bezier curve. Only " + to_string(control_points.size())  + " points have been provided.";
            throw std::logic_error(err);
        }
        control_points_ = std::move(control_points);
    }

    Point operator()(Real t) const
    {
        if (t < 0 || t > 1) {
            std::cerr << __FILE__ << ":" << __LINE__ << ":" << __func__ << "\n";
            std::cerr << "Querying the Bezier curve interpolator at t = " << t << " is not allowed; t in [0,1] is required.\n";
            Point p;
            for (Z i = 0; i < p.size(); ++i) {
                p[i] = std::numeric_limits<Real>::quiet_NaN();
            }
            return p;
        }

        auto P = control_points_[0];

        return P;
    }

    friend std::ostream& operator<<(std::ostream& out, bezier_polynomial_imp<RandomAccessContainer> const & bc) {
        return out;
    }

private:
    RandomAccessContainer control_points_;
};


}
#endif
