//  (C) Copyright Nick Thompson 2019.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_GEGENBAUER_HPP
#define BOOST_MATH_SPECIAL_GEGENBAUER_HPP

#include <stdexcept>
//#include <cmath>
//#include <boost/math/tools/recurrence.hpp>

namespace boost { namespace math {


template<typename Real>
Real gegenbauer(unsigned n, Real lambda, Real x)
{
    if (lambda <= -1/Real(2)) {
        throw std::domain_error("lambda > -1/2 is required.");
    }
    if (x < 0) {
        if (n&1) {
            return -gegenbauer(n, lambda, -x);
        }
        return gegenbauer(n, lambda, -x);
    }

    if (x > 1) {
        throw std::domain_error("Recurrence implementations for the Gegenbauer polynomials have unknown stability for x > 1.");
    }

    if (n == 0) {
        return Real(1);
    }
    if (n == 1) {
        return 2*lambda*x;
    }
    Real y0 = 1;
    Real y1 = 2*lambda*x;

    Real yk;
    for (unsigned k = 2; k <= n; ++k) {
        yk = ( 2*(k+lambda-1)*x*y1 - (k+2*lambda-2)*y0)/Real(k);
        y0 = y1;
        y1 = yk;
    }
    return yk;

    /*// (n+1)C_{n+1}(lambda, x) - (2x(n+lambda))C_{n}(lambda, x) + (n-1+2lambda)C_{n-1}(lambda, x) = 0
    auto Recurrence = [&](unsigned j) { return std::make_tuple<Real2, Real2, Real2>( n-1+2*lambda, -2*x*(n+lambda), n+1); };

    Real2 factor = 1;
    boost::uintmax_t max_iter = n + 1;
    Real2 Clambdax = boost::math::tools::function_ratio_from_backwards_recurrence(Recurrence, factor, max_iter);

    return Clambdax;*/
}

}}
#endif
