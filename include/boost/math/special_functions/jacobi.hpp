//  (C) Copyright Nick Thompson 2019.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_JACOBI_HPP
#define BOOST_MATH_SPECIAL_JACOBI_HPP

#include <limits>
#include <stdexcept>

namespace boost { namespace math {

template<typename Real>
Real jacobi(unsigned n, Real alpha, Real beta, Real x)
{
    if (n == 0) {
        return Real(1);
    }
    Real y0 = 1;
    Real y1 = (alpha+1) + (alpha+beta+2)*(x-1)/Real(2);

    Real yk = y1;
    Real k = 2;
    Real k_max = n*(1+std::numeric_limits<Real>::epsilon());
    while(k < k_max)
    {
        // Hoping for lots of common subexpression elimination by the compiler:
        Real denom = 2*k*(k+alpha+beta)*(2*k+alpha+beta-2);
        Real gamma1 = (2*k+alpha+beta-1)*( (2*k+alpha+beta)*(2*k+alpha+beta-2)*x + alpha*alpha -beta*beta);
        Real gamma0 = -2*(k+alpha-1)*(k+beta-1)*(2*k+alpha+beta);
        yk = (gamma1*y1 + gamma0*y0)/denom;
        y0 = y1;
        y1 = yk;
        k += 1;
    }
    return yk;
}

}}
#endif
