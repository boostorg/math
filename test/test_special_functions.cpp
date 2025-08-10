#include "/home/mzhelyez/Documents/code/boost/libs/math/include/boost/math/differentiation/Eigen/Eigen"
#include "test_autodiff_reverse.hpp"
#include <boost/math/tools/test_value.hpp>
#include <boost/utility/identity_type.hpp>
BOOST_AUTO_TEST_SUITE(test_special_functions)

using namespace rdiff;
BOOST_AUTO_TEST_CASE_TEMPLATE(test_tgamma, T, all_float_types)
{
    const T              eps = 1000 * boost::math::tools::epsilon<T>(); // percent
    constexpr unsigned   m   = 4; // 10 is too many derivatives for reverse mode
    const T              cx  = 3;
    // Mathematica: N[Table[D[Gamma[x],{x,n}] /. x->3, {n, 0, 10}], 52]
    std::array<T, m + 1> answers{
        {BOOST_MATH_TEST_VALUE(T, 2.0),
         BOOST_MATH_TEST_VALUE(T, 1.845568670196934278786975819835195137915681328120153),
         BOOST_MATH_TEST_VALUE(T, 2.492929991902693057942510065508124245503778067273315),
         BOOST_MATH_TEST_VALUE(T, 3.449965013523673365279327178241708777509009968597547),
         BOOST_MATH_TEST_VALUE(T, 5.521798578098737512443417699412265532987916790978887)}};
    rvar<T, m> x = make_rvar<T, m>(cx);
    rvar<T, m> y = tgamma(x);

    auto g1 = grad_nd<1>(y, &x);
    auto g2 = grad_nd<2>(y, &x);
    auto g3 = grad_nd<3>(y, &x);
    auto g4 = grad_nd<4>(y, &x);

    BOOST_CHECK_CLOSE(y.item(), answers[0], eps);
    BOOST_CHECK_CLOSE(g1[0]->item(), answers[1], eps);
    BOOST_CHECK_CLOSE(g2[0][0]->item(), answers[2], eps);
    BOOST_CHECK_CLOSE(g3[0][0][0]->item(), answers[3], eps);
    BOOST_CHECK_CLOSE(g4[0][0][0][0], answers[4], eps);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(tgamma2_test, T, all_float_types)
{
    //const T eps = 5000 * boost::math::tools::epsilon<T>(); // ok for non-multiprecision
    const T              eps = 500000 * boost::math::tools::epsilon<T>(); // percent
    constexpr unsigned   m   = 4; // 10 is too many derivatives for reverse mode
    const T              cx  = T(-1.5);
    // Mathematica: N[Table[D[Gamma[x],{x,n}] /. x->-3/2, {n, 0, 10}], 52]
    std::array<T, m + 1> answers{{
        BOOST_MATH_TEST_VALUE(T, 2.363271801207354703064223311121526910396732608163183),
        BOOST_MATH_TEST_VALUE(T, 1.661750260668596505586468565464938761014714509096807),
        BOOST_MATH_TEST_VALUE(T, 23.33417984355457252918927856618603412638766668207679),
        BOOST_MATH_TEST_VALUE(T, 47.02130025080143055642555842335081335790754507072526),
        BOOST_MATH_TEST_VALUE(T, 1148.336052788822231948472800239024335856568111484074),
    }};
    rvar<T, m>           x = make_rvar<T, m>(cx);
    rvar<T, m>           y = tgamma(x);

    auto g1 = grad_nd<1>(y, &x);
    auto g2 = grad_nd<2>(y, &x);
    auto g3 = grad_nd<3>(y, &x);
    auto g4 = grad_nd<4>(y, &x);

    BOOST_CHECK_CLOSE(y.item(), answers[0], eps);
    BOOST_CHECK_CLOSE(g1[0]->item(), answers[1], eps);
    BOOST_CHECK_CLOSE(g2[0][0]->item(), answers[2], eps);
    BOOST_CHECK_CLOSE(g3[0][0][0]->item(), answers[3], eps);
    BOOST_CHECK_CLOSE(g4[0][0][0][0], answers[4], eps);
}
BOOST_AUTO_TEST_CASE_TEMPLATE(t_mm, T, all_float_types)
{
    Eigen::Matrix<rvar<T, 1>, 2, 2> A, B, C;
    rvar<T, 1>                      a00, a01, a10, a11;
    rvar<T, 1>                      b00, b01, b10, b11;

    a00 = 1.0;
    a01 = 0.0;
    a10 = 1.0;
    a11 = 1.0;
    b00 = 2.0;
    b01 = 3.0;
    b10 = 1.0;
    b11 = std::sqrt(2);

    A(0, 0) = a00;
    A(0, 1) = a01;
    A(1, 0) = a10;
    A(1, 1) = a11;

    B(0, 0) = b00;
    B(0, 1) = b01;
    B(1, 0) = b10;
    B(1, 1) = b11;

    C = A * B;

    rvar<T, 1> z;
    z = C(0, 0) + C(1, 0) + C(0, 1) + C(1, 1);
    z.backward();
    std::cout << a00.adjoint() << std::endl;
}
BOOST_AUTO_TEST_SUITE_END()
