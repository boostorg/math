//  (C) Copyright Nick Thompson 2019.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_CARDINAL_B_SPLINE_HPP
#define BOOST_MATH_SPECIAL_CARDINAL_B_SPLINE_HPP

#include <array>
#include <cmath>
#include <boost/math/tools/recurrence.hpp>

namespace boost { namespace math {


template<typename Real1, typename Real2>
Real2 gegenbauer(unsigned n, Real1 alpha, Real2 x)
{
    if (x < 0) {
        if (n&1) {
            return -gegenbauer(n, alpha, -x);
        }
        return gegenbauer(n, alpha, -x);
    }

    if (x > 1) {
        throw std::domain_error("Recurrence implementations for the Gegenbauer polynomials have unknown stability for x > 1.");
    }

    // (n+1)C_{n+1}(alpha, x) - (2x(n+alpha))C_{n}(alpha, x) + (n-1+2alpha)C_{n-1}(alpha, x) = 0
    auto Recurrence = [&](unsigned j) { return std::make_tuple<Real2, Real2, Real2>( n-1+2*alpha, -2*x*(n+alpha), n+1); };

    Real2 factor = 1;
    boost::uintmax_t max_iter = n + 1;
    Real2 Calphax = boost::math::tools::function_ratio_from_backwards_recurrence(Recurrence, factor, max_iter);

    return Calphax;
}

}}
#endif
