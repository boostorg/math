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

namespace boost::math::interpolators {

template<class RandomAccessContainer>
class makima {
public:
    using Real = RandomAccessContainer::value_type;

    makima(RandomAccessContainer const & x, RandomAccessContainer const & y) {
        if (x.size() != y.size()) {
            throw std::domain_error("There must be the same number of ordinates as abscissas.");
        }

    }

    Real operator()(Real t) const {
        return std::numeric_limits<Real>::quiet_NaN();
    }    

private:

};
}
#endif