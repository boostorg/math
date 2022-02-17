//  (C) Copyright Matt Borland 2022.
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

#include "math_unit_test.hpp"
#include <cmath>
#include <vector>
#include <boost/math/special_functions/logsumexp.hpp>

template <typename Real>
void test()
{
    using boost::math::logsumexp;
    using std::log;
    using std::exp;

    // Spot check 2 values
    // Also validate that 2 values does not attempt to instantiate the iterator version
    // https://numpy.org/doc/stable/reference/generated/numpy.logaddexp.html
    // Calculated at higher precision using wolfram alpha
    Real spot1 = static_cast<Real>(log(1e-50l));
    Real spot2 = static_cast<Real>(log(2.5e-50l));
    Real spot12 = logsumexp(spot1, spot2);

    CHECK_ULP_CLOSE(Real(-113.87649168120691620521145211223320721849348983622l), spot12, 1);

    // Spot check 3 values and compare result of each different interface
    Real spot3 = static_cast<Real>(log(5e-50l));
    std::vector<Real> spots {spot1, spot2, spot3};

    Real spot123 = logsumexp(spot1, spot2, spot3);
    Real spot123_container = logsumexp(spots);
    Real spot123_iter = logsumexp(spots.begin(), spots.end());

    CHECK_EQUAL(spot123, spot123_container);
    CHECK_EQUAL(spot123_container, spot123_iter);
    CHECK_ULP_CLOSE(Real(-112.98918848620601343006727023780326041254237155321l), spot123, 1);

    // Spot check 4 values with repeated largest value
    Real spot4 = spot3;

    Real spot1234 = logsumexp(spot1, spot2, spot3, spot4);
    spots.emplace_back(spot4);
    Real spot1234_container = logsumexp(spots);

    CHECK_EQUAL(spot1234, spot1234_container);
    CHECK_ULP_CLOSE(Real(-112.52656496425790043613106914490880983418810289233l), spot1234, 1);
}

int main (void)
{
    test<float>();
    test<double>();
    test<long double>();
}
