//  (C) Copyright Matt Borland 2022.
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

#include "math_unit_test.hpp"
#include <limits>
#include <boost/math/special_functions/logaddexp.hpp>
#include <boost/math/constants/constants.hpp>

template <typename Real>
void test()
{
    using boost::math::logaddexp;
    using std::log;
    using std::exp;

    constexpr Real nan_val = std::numeric_limits<Real>::quiet_NaN();
    constexpr Real inf_val = std::numeric_limits<Real>::infinity();

    // NAN
    CHECK_NAN(logaddexp(nan_val, Real(1)));
    CHECK_NAN(logaddexp(Real(1), nan_val));
    CHECK_NAN(logaddexp(nan_val, nan_val));

    // INF
    CHECK_EQUAL(logaddexp(inf_val, Real(1)), inf_val);
    CHECK_EQUAL(logaddexp(Real(1), inf_val), inf_val);
    CHECK_EQUAL(logaddexp(inf_val, inf_val), inf_val);

    // Equal values
    constexpr Real ln2 = boost::math::constants::ln_two<Real>();
    CHECK_ULP_CLOSE(Real(2) + ln2, logaddexp(Real(2), Real(2)), 1);
    CHECK_ULP_CLOSE(Real(1e-50) + ln2, logaddexp(Real(1e-50), Real(1e-50)), 1);

    // Spot check
    // https://numpy.org/doc/stable/reference/generated/numpy.logaddexp.html
    // Calculated at higher precision using wolfram alpha
    Real spot1 = static_cast<Real>(log(1e-50l));
    Real spot2 = static_cast<Real>(log(2.5e-50l));
    Real spot12 = logaddexp(spot1, spot2);

    CHECK_ULP_CLOSE(Real(-113.87649168120691620521145211223320721849348983622l), spot12, 1);

    // GCCs ASAN does not like this test so suppress when active
    // See: https://drone.cpp.al/boostorg/math/607/1/2 line 1508
    #ifndef __SANITIZE_ADDRESS__
    CHECK_ULP_CLOSE(Real(3.5e-50l), exp(spot12), 1);
    #endif
}

int main (void)
{
    test<float>();
    test<double>();
    test<long double>();
}
