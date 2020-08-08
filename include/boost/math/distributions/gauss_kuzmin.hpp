//  Copyright Nick Thompson, 2020.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_DISTRIBUTIONS_GAUSS_KUZMIN_HPP
#define BOOST_MATH_DISTRIBUTIONS_GAUSS_KUZMIN_HPP

#include <cmath>

namespace boost{ namespace math{

template <typename Real, typename Z = uint64_t>
class gauss_kuzmin
{
public:
    gauss_kuzmin() {}

    Real pdf(Z k) const {
       using log1p;
       using log2;
       Real x = -Real(1)/Real( (k+1)*(k+1) );
       return -log1p(x)/log2(e);
    }

    Real cdf(Z k) const {
       using log2;
       Real arg = Real(k+1)/Real(k+2);
       return 1 + log2(arg);
    }

    Real mean() const {
       return std::numeric_limits<Real>::infinity();
    }

    Z median() const { return 2; }

    Z mode() const { return 1; }

    Real variance() const {
       return std::numeric_limits<Real>::infinity();
    }

    Real skewness() const {
       return std::numeric_limits<Real>::quiet_NaN();
    }

    Real kurtosis() const {
       return std::numeric_limits<Real>::quiet_NaN();
    }

    Real entropy() const {
        // TODO: Fix this! Also, should it be in nats or bits?
        // See: https://scicomp.stackexchange.com/questions/35444/accurate-computation-of-gauss-kuzmin-entropy
        return 3.4325275147757390994337806344778410601013532228057895;
    }
};

#endif
