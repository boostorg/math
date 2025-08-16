//           Copyright Maksym Zhelyenzyakov 2025-2026.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)
#include "test_autodiff_reverse.hpp"
#include <boost/math/tools/test_value.hpp>
#include <boost/utility/identity_type.hpp>
#include <cmath>
BOOST_AUTO_TEST_SUITE(erf_support)

using namespace rdiff;
BOOST_AUTO_TEST_CASE_TEMPLATE(test_erf, T, all_float_types)
{
    std::array<T, 5> answers{
        {BOOST_MATH_TEST_VALUE(
             T, 0.99997790950300141455862722387041767962015229291260075034275901649691099714918),
         BOOST_MATH_TEST_VALUE(T, 0.0001392530519467478538904180310183558561119961779459213),
         BOOST_MATH_TEST_VALUE(T, -0.0008355183116804871233425081861101351366719770676755276),
         BOOST_MATH_TEST_VALUE(T, 0.004734603766189427032274213054624099107807870050161323),
         BOOST_MATH_TEST_VALUE(T, -0.02506554935041461370027524558330405410015931203026583)}};
    const T       eps = 0.01;
    constexpr int ord = 4;
    T             cx  = 3.0;

    rvar<T, ord> x  = cx;
    rvar<T, ord> y  = rdiff::erf(x);
    auto         g1 = grad_nd<1>(y, &x);
    auto         g2 = grad_nd<2>(y, &x);
    auto         g3 = grad_nd<3>(y, &x);
    auto         g4 = grad_nd<4>(y, &x);

    BOOST_CHECK_CLOSE(y.item(), answers[0], eps);
    BOOST_CHECK_CLOSE(g1[0], answers[1], eps);
    BOOST_CHECK_CLOSE(g2[0][0], answers[2], eps);
    BOOST_CHECK_CLOSE(g3[0][0][0], answers[3], eps);
    BOOST_CHECK_CLOSE(g4[0][0][0][0], answers[4], eps);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_erfc, T, all_float_types)
{
    std::array<T, 5> answers{
        {BOOST_MATH_TEST_VALUE(T, 0.00002209049699858544137277612958232037984770708739924966),
         BOOST_MATH_TEST_VALUE(T, -0.0001392530519467478538904180310183558561119961779459213),
         BOOST_MATH_TEST_VALUE(T, 0.0008355183116804871233425081861101351366719770676755276),
         BOOST_MATH_TEST_VALUE(T, -0.004734603766189427032274213054624099107807870050161323),
         BOOST_MATH_TEST_VALUE(T, 0.02506554935041461370027524558330405410015931203026583)}};
    const T       eps = 0.01;
    constexpr int ord = 4;
    T             cx  = 3.0;

    rvar<T, ord> x  = cx;
    rvar<T, ord> y  = rdiff::erfc(x);
    auto         g1 = grad_nd<1>(y, &x);
    auto         g2 = grad_nd<2>(y, &x);
    auto         g3 = grad_nd<3>(y, &x);
    auto         g4 = grad_nd<4>(y, &x);

    BOOST_CHECK_CLOSE(y.item(), answers[0], eps);
    BOOST_CHECK_CLOSE(g1[0], answers[1], eps);
    BOOST_CHECK_CLOSE(g2[0][0], answers[2], eps);
    BOOST_CHECK_CLOSE(g3[0][0][0], answers[3], eps);
    BOOST_CHECK_CLOSE(g4[0][0][0][0], answers[4], eps);
}

using namespace rdiff;
BOOST_AUTO_TEST_CASE_TEMPLATE(test_erf_inv, T, all_float_types)
{
    std::array<T, 5> answers{{BOOST_MATH_TEST_VALUE(T, 0.4769362762044699),
                              BOOST_MATH_TEST_VALUE(T, 1.11258481897195),
                              BOOST_MATH_TEST_VALUE(T, 1.1807463499934),
                              BOOST_MATH_TEST_VALUE(T, 5.26058253926122),
                              BOOST_MATH_TEST_VALUE(T, 28.44125004446708)}};
    const T          eps = 0.01;
    constexpr int    ord = 4;
    T                cx  = 0.5;

    rvar<T, ord> x  = cx;
    rvar<T, ord> y  = rdiff::erf_inv(x);
    auto         g1 = grad_nd<1>(y, &x);
    auto         g2 = grad_nd<2>(y, &x);
    auto         g3 = grad_nd<3>(y, &x);
    auto         g4 = grad_nd<4>(y, &x);

    BOOST_CHECK_CLOSE(y.item(), answers[0], eps);
    BOOST_CHECK_CLOSE(g1[0], answers[1], eps);
    BOOST_CHECK_CLOSE(g2[0][0], answers[2], eps);
    BOOST_CHECK_CLOSE(g3[0][0][0], answers[3], eps);
    BOOST_CHECK_CLOSE(g4[0][0][0][0], answers[4], eps);
}
BOOST_AUTO_TEST_CASE_TEMPLATE(test_erfc_inv, T, all_float_types)
{
    std::array<T, 5> answers{{BOOST_MATH_TEST_VALUE(T, 0.4769362762044699),
                              BOOST_MATH_TEST_VALUE(T, -1.11258481897195),
                              BOOST_MATH_TEST_VALUE(T, 1.1807463499934),
                              BOOST_MATH_TEST_VALUE(T, -5.260582539261222),
                              BOOST_MATH_TEST_VALUE(T, 28.44125004446708)}};
    const T          eps = 0.01;
    constexpr int    ord = 4;
    T                cx  = 0.5;

    rvar<T, ord> x  = cx;
    rvar<T, ord> y  = rdiff::erfc_inv(x);
    auto         g1 = grad_nd<1>(y, &x);
    auto         g2 = grad_nd<2>(y, &x);
    auto         g3 = grad_nd<3>(y, &x);
    auto         g4 = grad_nd<4>(y, &x);

    BOOST_CHECK_CLOSE(y.item(), answers[0], eps);
    BOOST_CHECK_CLOSE(g1[0], answers[1], eps);
    BOOST_CHECK_CLOSE(g2[0][0], answers[2], eps);
    BOOST_CHECK_CLOSE(g3[0][0][0], answers[3], eps);
    BOOST_CHECK_CLOSE(g4[0][0][0][0], answers[4], eps);
}
BOOST_AUTO_TEST_SUITE_END()
