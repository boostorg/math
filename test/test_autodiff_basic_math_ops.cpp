/*#include <boost/math/differentiation/autodiff_reverse.hpp>
using namespace boost::math::differentiation::reverse_mode;
int main()
{
    rvar<double,1> x = 1.5;
    auto z = x*x*x*x*x;
    rvar<double,1>zz(z);
    zz.backward();
    std::cout<<x.adjoint()<<std::endl;
}*/

#include "test_autodiff_reverse.hpp"
#include <vector>
BOOST_AUTO_TEST_SUITE(test_basic_math_ops)

using namespace rdiff;
template<typename T>
T test_add(T x[5])
{
    return x[0] + x[1] + x[2] + x[3] + x[4];
}
BOOST_AUTO_TEST_CASE_TEMPLATE(addition, T, all_float_types)
{
    RandomSample<T> rng{-100, 100};
    T               x1_v               = rng.next();
    T               x2_v               = rng.next();
    T               test_rvar_p_rvar_v = x1_v + x2_v;

    rvar<T, 1>      x1                 = x1_v;
    rvar<T, 1>      x2                 = x2_v;

    rvar<T, 1>      test_rvar_p_rvar   = x1 + x2;
    rvar<T, 1>      test_rvar_p_float  = x1 + x2_v;
    rvar<T, 1>      test_float_p_rvar  = x1_v + x2;
    BOOST_REQUIRE_EQUAL(test_rvar_p_rvar_v, test_rvar_p_rvar.item());
    BOOST_REQUIRE_EQUAL(test_rvar_p_rvar_v, test_rvar_p_float.item());
    BOOST_REQUIRE_EQUAL(test_float_p_rvar.item(), test_rvar_p_rvar_v);

    gradient_tape<T, 1>& tape = get_active_tape<T, 1>();
    tape.zero_grad();
    test_rvar_p_rvar.backward();
    BOOST_REQUIRE_EQUAL(x1.adjoint(), 1.0);
    BOOST_REQUIRE_EQUAL(x2.adjoint(), 1.0);

    tape.zero_grad();
    rvar<T, 1> z = x1 + x1 + x1 + x1 + x1;
    z.backward();
    BOOST_REQUIRE_EQUAL(x1.adjoint(), 5.0);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(multiplication, T, all_float_types)
{
    RandomSample<T> rng{-100, 100};
    T               x1_v               = rng.next();
    T               x2_v               = rng.next();
    T               test_rvar_p_rvar_v = x1_v * x2_v;

    rvar<T, 1>      x1                 = x1_v;
    rvar<T, 1>      x2                 = x2_v;

    rvar<T, 1>      test_rvar_p_rvar   = x1 * x2;
    rvar<T, 1>      test_rvar_p_float  = x1 * x2_v;
    rvar<T, 1>      test_float_p_rvar  = x1_v * x2;

    BOOST_REQUIRE_EQUAL(test_rvar_p_rvar_v, test_rvar_p_rvar.item());
    BOOST_REQUIRE_EQUAL(test_rvar_p_rvar_v, test_rvar_p_float.item());
    BOOST_REQUIRE_EQUAL(test_float_p_rvar.item(), test_rvar_p_rvar_v);

    gradient_tape<T, 1>& tape = get_active_tape<T, 1>();
    tape.zero_grad();
    test_rvar_p_rvar.backward();
    BOOST_REQUIRE_EQUAL(x1.adjoint(), x2_v);
    BOOST_REQUIRE_EQUAL(x2.adjoint(), x1_v);

    tape.zero_grad();
    auto       z = x1 * x1 * x1 * x1 * x1;
    rvar<T, 1> zz(z);
    zz.backward();
    BOOST_REQUIRE_CLOSE(x1.adjoint(), 5 * x1_v * x1_v * x1_v * x1_v, boost_close_tol<T>());
    tape.clear();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(minus, T, all_float_types)
{
    RandomSample<T> rng{-100, 100};
    T               x1_v               = rng.next();
    T               x2_v               = rng.next();
    T               test_rvar_p_rvar_v = x1_v - x2_v;

    rvar<T, 1>      x1                 = x1_v;
    rvar<T, 1>      x2                 = x2_v;

    rvar<T, 1>      test_rvar_p_rvar   = x1 - x2;
    rvar<T, 1>      test_rvar_p_float  = x1 - x2_v;
    rvar<T, 1>      test_float_p_rvar  = x1_v - x2;

    BOOST_REQUIRE_EQUAL(test_rvar_p_rvar_v, test_rvar_p_rvar.item());
    BOOST_REQUIRE_EQUAL(test_rvar_p_rvar_v, test_rvar_p_float.item());
    BOOST_REQUIRE_EQUAL(test_float_p_rvar.item(), test_rvar_p_rvar_v);

    gradient_tape<T, 1>& tape = get_active_tape<T, 1>();
    tape.zero_grad();
    test_rvar_p_rvar.backward();
    BOOST_REQUIRE_EQUAL(x1.adjoint(), 1.0);
    BOOST_REQUIRE_EQUAL(x2.adjoint(), -1.0);

    tape.zero_grad();
    auto       z = -x1 - x1 - x1 - x1 - x1;
    rvar<T, 1> zz(z);
    zz.backward();
    BOOST_REQUIRE_CLOSE(x1.adjoint(), -5, boost_close_tol<T>());
    tape.clear();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(division, T, all_float_types)
{
    RandomSample<T> rng{-1, 1};
    T               x1_v               = rng.next();
    T               x2_v               = rng.next();
    T               test_rvar_p_rvar_v = x1_v / x2_v;

    rvar<T, 1>      x1                 = x1_v;
    rvar<T, 1>      x2                 = x2_v;

    rvar<T, 1>      test_rvar_p_rvar   = x1 / x2;
    rvar<T, 1>      test_rvar_p_float  = x1 / x2_v;
    rvar<T, 1>      test_float_p_rvar  = x1_v / x2;

    BOOST_REQUIRE_EQUAL(test_rvar_p_rvar_v, test_rvar_p_rvar.item());
    BOOST_REQUIRE_EQUAL(test_rvar_p_rvar_v, test_rvar_p_float.item());
    BOOST_REQUIRE_EQUAL(test_float_p_rvar.item(), test_rvar_p_rvar_v);

    gradient_tape<T, 1>& tape = get_active_tape<T, 1>();
    tape.zero_grad();
    test_rvar_p_rvar.backward();
    BOOST_REQUIRE_CLOSE(x1.adjoint(), 1.0 / x2_v, boost_close_tol<T>());
    BOOST_REQUIRE_CLOSE(x2.adjoint(), -x1_v / (x2_v * x2_v), boost_close_tol<T>());

    tape.zero_grad();
    auto       z = 2.0 / x1 / x1 / x1 / x1 / x1;
    rvar<T, 1> zz(z);
    zz.backward();
    BOOST_REQUIRE_CLOSE(x1.adjoint(),
                        -10 / (x1_v * x1_v * x1_v * x1_v * x1_v * x1_v),
                        boost_close_tol<T>());
    tape.clear();
}

/* gradient test function */
template<typename T>
T f(T x, T y)
{
    //auto z1 = x / (y + x);
    //auto z2 = y / (x - y);
    //T    f  = z1 * z2;
    return (x * y) / ((x + y) * (x - y));
}

template<typename T>
std::vector<T> grad_f_a(T x, T y)
{
    T f_x = static_cast<T>(y * (-pow(x, 2) - pow(y, 2))
                           / (pow(x, 4) - 2 * pow(x, 2) * pow(y, 2) + pow(y, 4)));
    T f_y = static_cast<T>(x * (pow(x, 2) + pow(y, 2))
                           / (pow(x, 4) - 2 * pow(x, 2) * pow(y, 2) + pow(y, 4)));
    return {f_x, f_y};
}

template<typename T>
std::vector<std::vector<T>> hess_f_a(T x, T y)
{
    T f_xx = static_cast<T>(
        2 * x * y * (pow(x, 2) + 3 * pow(y, 2))
        / (pow(x, 6) - 3 * pow(x, 4) * pow(y, 2) + 3 * pow(x, 2) * pow(y, 4) - pow(y, 6)));
    T f_xy = static_cast<T>(
        (-pow(x, 4) - 6 * pow(x, 2) * pow(y, 2) - pow(y, 4))
        / (pow(x, 6) - 3 * pow(x, 4) * pow(y, 2) + 3 * pow(x, 2) * pow(y, 4) - pow(y, 6)));
    T f_yx = static_cast<T>(
        (-pow(x, 4) - 6 * pow(x, 2) * pow(y, 2) - pow(y, 4))
        / (pow(x, 6) - 3 * pow(x, 4) * pow(y, 2) + 3 * pow(x, 2) * pow(y, 4) - pow(y, 6)));
    T f_yy = static_cast<T>(
        2 * x * y * (3 * pow(x, 2) + pow(y, 2))
        / (pow(x, 6) - 3 * pow(x, 4) * pow(y, 2) + 3 * pow(x, 2) * pow(y, 4) - pow(y, 6)));
    return {{f_xx, f_xy}, {f_yx, f_yy}};
}

template<typename T>
std::vector<std::vector<std::vector<T>>> df_3_a(T x, T y)
{
    T f_xxx = static_cast<T>(6 * y * (-pow(x, 4) - 6 * pow(x, 2) * pow(y, 2) - pow(y, 4))
                             / (pow(x, 8) - 4 * pow(x, 6) * pow(y, 2) + 6 * pow(x, 4) * pow(y, 4)
                                - 4 * pow(x, 2) * pow(y, 6) + pow(y, 8)));
    T f_xxy = static_cast<T>(2 * x * (pow(x, 4) + 14 * pow(x, 2) * pow(y, 2) + 9 * pow(y, 4))
                             / (pow(x, 8) - 4 * pow(x, 6) * pow(y, 2) + 6 * pow(x, 4) * pow(y, 4)
                                - 4 * pow(x, 2) * pow(y, 6) + pow(y, 8)));
    T f_xyx = static_cast<T>(2 * x * (pow(x, 4) + 14 * pow(x, 2) * pow(y, 2) + 9 * pow(y, 4))
                             / (pow(x, 8) - 4 * pow(x, 6) * pow(y, 2) + 6 * pow(x, 4) * pow(y, 4)
                                - 4 * pow(x, 2) * pow(y, 6) + pow(y, 8)));
    T f_xyy = static_cast<T>(2 * y * (-9 * pow(x, 4) - 14 * pow(x, 2) * pow(y, 2) - pow(y, 4))
                             / (pow(x, 8) - 4 * pow(x, 6) * pow(y, 2) + 6 * pow(x, 4) * pow(y, 4)
                                - 4 * pow(x, 2) * pow(y, 6) + pow(y, 8)));
    T f_yxx = static_cast<T>(2 * x * (pow(x, 4) + 14 * pow(x, 2) * pow(y, 2) + 9 * pow(y, 4))
                             / (pow(x, 8) - 4 * pow(x, 6) * pow(y, 2) + 6 * pow(x, 4) * pow(y, 4)
                                - 4 * pow(x, 2) * pow(y, 6) + pow(y, 8)));
    T f_yxy = static_cast<T>(2 * y * (-9 * pow(x, 4) - 14 * pow(x, 2) * pow(y, 2) - pow(y, 4))
                             / (pow(x, 8) - 4 * pow(x, 6) * pow(y, 2) + 6 * pow(x, 4) * pow(y, 4)
                                - 4 * pow(x, 2) * pow(y, 6) + pow(y, 8)));
    T f_yyx = static_cast<T>(2 * y * (-9 * pow(x, 4) - 14 * pow(x, 2) * pow(y, 2) - pow(y, 4))
                             / (pow(x, 8) - 4 * pow(x, 6) * pow(y, 2) + 6 * pow(x, 4) * pow(y, 4)
                                - 4 * pow(x, 2) * pow(y, 6) + pow(y, 8)));
    T f_yyy = static_cast<T>(6 * x * (pow(x, 4) + 6 * pow(x, 2) * pow(y, 2) + pow(y, 4))
                             / (pow(x, 8) - 4 * pow(x, 6) * pow(y, 2) + 6 * pow(x, 4) * pow(y, 4)
                                - 4 * pow(x, 2) * pow(y, 6) + pow(y, 8)));
    return {{{f_xxx, f_xxy}, {f_xyx, f_xyy}}, {{f_yxx, f_yxy}, {f_yyx, f_yyy}}};
}
template<typename T>
std::vector<std::vector<std::vector<std::vector<T>>>> df_4_a(T x, T y)
{
    T f_xxxx = static_cast<T>(24 * x * y * (pow(x, 4) + 10 * pow(x, 2) * pow(y, 2) + 5 * pow(y, 4))
                              / (pow(x, 10) - 5 * pow(x, 8) * pow(y, 2) + 10 * pow(x, 6) * pow(y, 4)
                                 - 10 * pow(x, 4) * pow(y, 6) + 5 * pow(x, 2) * pow(y, 8)
                                 - pow(y, 10)));
    T f_xxxy = static_cast<T>(
        6 * (-pow(x, 6) - 25 * pow(x, 4) * pow(y, 2) - 35 * pow(x, 2) * pow(y, 4) - 3 * pow(y, 6))
        / (pow(x, 10) - 5 * pow(x, 8) * pow(y, 2) + 10 * pow(x, 6) * pow(y, 4)
           - 10 * pow(x, 4) * pow(y, 6) + 5 * pow(x, 2) * pow(y, 8) - pow(y, 10)));
    T f_xxyx = static_cast<T>(
        6 * (-pow(x, 6) - 25 * pow(x, 4) * pow(y, 2) - 35 * pow(x, 2) * pow(y, 4) - 3 * pow(y, 6))
        / (pow(x, 10) - 5 * pow(x, 8) * pow(y, 2) + 10 * pow(x, 6) * pow(y, 4)
           - 10 * pow(x, 4) * pow(y, 6) + 5 * pow(x, 2) * pow(y, 8) - pow(y, 10)));
    T f_xxyy = static_cast<T>(
        (72 * pow(x, 5) * y + 240 * pow(x, 3) * pow(y, 3) + 72 * x * pow(y, 5))
        / (pow(x, 10) - 5 * pow(x, 8) * pow(y, 2) + 10 * pow(x, 6) * pow(y, 4)
           - 10 * pow(x, 4) * pow(y, 6) + 5 * pow(x, 2) * pow(y, 8) - pow(y, 10)));
    T f_xyxx = static_cast<T>(
        6 * (-pow(x, 6) - 25 * pow(x, 4) * pow(y, 2) - 35 * pow(x, 2) * pow(y, 4) - 3 * pow(y, 6))
        / (pow(x, 10) - 5 * pow(x, 8) * pow(y, 2) + 10 * pow(x, 6) * pow(y, 4)
           - 10 * pow(x, 4) * pow(y, 6) + 5 * pow(x, 2) * pow(y, 8) - pow(y, 10)));
    T f_xyxy = static_cast<T>(
        (72 * pow(x, 5) * y + 240 * pow(x, 3) * pow(y, 3) + 72 * x * pow(y, 5))
        / (pow(x, 10) - 5 * pow(x, 8) * pow(y, 2) + 10 * pow(x, 6) * pow(y, 4)
           - 10 * pow(x, 4) * pow(y, 6) + 5 * pow(x, 2) * pow(y, 8) - pow(y, 10)));
    T f_xyyx = static_cast<T>(
        (72 * pow(x, 5) * y + 240 * pow(x, 3) * pow(y, 3) + 72 * x * pow(y, 5))
        / (pow(x, 10) - 5 * pow(x, 8) * pow(y, 2) + 10 * pow(x, 6) * pow(y, 4)
           - 10 * pow(x, 4) * pow(y, 6) + 5 * pow(x, 2) * pow(y, 8) - pow(y, 10)));
    T f_xyyy = static_cast<T>(
        6 * (-3 * pow(x, 6) - 35 * pow(x, 4) * pow(y, 2) - 25 * pow(x, 2) * pow(y, 4) - pow(y, 6))
        / (pow(x, 10) - 5 * pow(x, 8) * pow(y, 2) + 10 * pow(x, 6) * pow(y, 4)
           - 10 * pow(x, 4) * pow(y, 6) + 5 * pow(x, 2) * pow(y, 8) - pow(y, 10)));
    T f_yxxx = static_cast<T>(
        6 * (-pow(x, 6) - 25 * pow(x, 4) * pow(y, 2) - 35 * pow(x, 2) * pow(y, 4) - 3 * pow(y, 6))
        / (pow(x, 10) - 5 * pow(x, 8) * pow(y, 2) + 10 * pow(x, 6) * pow(y, 4)
           - 10 * pow(x, 4) * pow(y, 6) + 5 * pow(x, 2) * pow(y, 8) - pow(y, 10)));
    T f_yxxy = static_cast<T>(
        (72 * pow(x, 5) * y + 240 * pow(x, 3) * pow(y, 3) + 72 * x * pow(y, 5))
        / (pow(x, 10) - 5 * pow(x, 8) * pow(y, 2) + 10 * pow(x, 6) * pow(y, 4)
           - 10 * pow(x, 4) * pow(y, 6) + 5 * pow(x, 2) * pow(y, 8) - pow(y, 10)));
    T f_yxyx = static_cast<T>(
        (72 * pow(x, 5) * y + 240 * pow(x, 3) * pow(y, 3) + 72 * x * pow(y, 5))
        / (pow(x, 10) - 5 * pow(x, 8) * pow(y, 2) + 10 * pow(x, 6) * pow(y, 4)
           - 10 * pow(x, 4) * pow(y, 6) + 5 * pow(x, 2) * pow(y, 8) - pow(y, 10)));
    T f_yxyy = static_cast<T>(
        6 * (-3 * pow(x, 6) - 35 * pow(x, 4) * pow(y, 2) - 25 * pow(x, 2) * pow(y, 4) - pow(y, 6))
        / (pow(x, 10) - 5 * pow(x, 8) * pow(y, 2) + 10 * pow(x, 6) * pow(y, 4)
           - 10 * pow(x, 4) * pow(y, 6) + 5 * pow(x, 2) * pow(y, 8) - pow(y, 10)));
    T f_yyxx = static_cast<T>(
        (72 * pow(x, 5) * y + 240 * pow(x, 3) * pow(y, 3) + 72 * x * pow(y, 5))
        / (pow(x, 10) - 5 * pow(x, 8) * pow(y, 2) + 10 * pow(x, 6) * pow(y, 4)
           - 10 * pow(x, 4) * pow(y, 6) + 5 * pow(x, 2) * pow(y, 8) - pow(y, 10)));
    T f_yyxy = static_cast<T>(
        6 * (-3 * pow(x, 6) - 35 * pow(x, 4) * pow(y, 2) - 25 * pow(x, 2) * pow(y, 4) - pow(y, 6))
        / (pow(x, 10) - 5 * pow(x, 8) * pow(y, 2) + 10 * pow(x, 6) * pow(y, 4)
           - 10 * pow(x, 4) * pow(y, 6) + 5 * pow(x, 2) * pow(y, 8) - pow(y, 10)));
    T f_yyyx = static_cast<T>(
        6 * (-3 * pow(x, 6) - 35 * pow(x, 4) * pow(y, 2) - 25 * pow(x, 2) * pow(y, 4) - pow(y, 6))
        / (pow(x, 10) - 5 * pow(x, 8) * pow(y, 2) + 10 * pow(x, 6) * pow(y, 4)
           - 10 * pow(x, 4) * pow(y, 6) + 5 * pow(x, 2) * pow(y, 8) - pow(y, 10)));
    T f_yyyy = static_cast<T>(24 * x * y * (5 * pow(x, 4) + 10 * pow(x, 2) * pow(y, 2) + pow(y, 4))
                              / (pow(x, 10) - 5 * pow(x, 8) * pow(y, 2) + 10 * pow(x, 6) * pow(y, 4)
                                 - 10 * pow(x, 4) * pow(y, 6) + 5 * pow(x, 2) * pow(y, 8)
                                 - pow(y, 10)));
    return {{{{f_xxxx, f_xxxy}, {f_xxyx, f_xxyy}}, {{f_xyxx, f_xyxy}, {f_xyyx, f_xyyy}}},
            {{{f_yxxx, f_yxxy}, {f_yxyx, f_yxyy}}, {{f_yyxx, f_yyxy}, {f_yyyx, f_yyyy}}}};
}

BOOST_AUTO_TEST_CASE_TEMPLATE(first_derivative, T, all_float_types)
{
    //RandomSample<T> rng{-1, 1};
    T x = 1.1224;   //rng.next();
    T y = 5.213414; //rng.next();

    rvar<T, 1>      x_ad              = x;
    rvar<T, 1>      y_ad              = y;

    T               fv                = f(x, y);
    rvar<T, 1>      f_ad              = f(x_ad, y_ad);
    auto            grad_f_analytical = grad_f_a(x, y);
    /* intended use case */
    f_ad.backward();
    BOOST_REQUIRE_CLOSE(x_ad.adjoint(), grad_f_analytical[0], boost_close_tol<T>());
    BOOST_REQUIRE_CLOSE(y_ad.adjoint(), grad_f_analytical[1], boost_close_tol<T>());

    gradient_tape<T, 1, BOOST_MATH_BUFFER_SIZE>& tape = get_active_tape<T, 1>();
    tape.zero_grad();

    /* grad test */
    auto grad_func_test_grad    = grad(f_ad, &x_ad, &y_ad);
    /* grad_nd test */
    auto grad_nd_func_test_grad = grad_nd<1>(f_ad, &x_ad, &y_ad);

    for (int i = 0; i < 2; i++) {
        BOOST_REQUIRE_CLOSE(grad_func_test_grad[i], grad_f_analytical[i], boost_close_tol<T>());
        BOOST_REQUIRE_CLOSE(grad_nd_func_test_grad[i], grad_f_analytical[i], boost_close_tol<T>());
    }
}
BOOST_AUTO_TEST_CASE_TEMPLATE(second_derivative_and_hessian, T, all_float_types)
{
    //RandomSample<T>                   rng{1, 1};
    T x = 1.0; //; rng.next();
    T y = 1.5; //rng.next();

    rvar<T, 2>                        x_ad = x;
    rvar<T, 2>                        y_ad = y;

    T                                 fv   = f(x, y);
    rvar<T, 2>                        f_ad = f(x_ad, y_ad);
    gradient_tape<T, 2, BOOST_MATH_BUFFER_SIZE>& tape = get_active_tape<T, 2>();
    tape.zero_grad();

    auto hess_analytical   = hess_f_a(x, y);
    auto hess_func_test    = hess(f_ad, &x_ad, &y_ad);
    auto grad_nd_func_test = grad_nd<2>(f_ad, &x_ad, &y_ad);

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            BOOST_REQUIRE_CLOSE(hess_func_test[i][j], hess_analytical[i][j], boost_close_tol<T>());
            BOOST_REQUIRE_CLOSE(hess_func_test[i][j], grad_nd_func_test[i][j], boost_close_tol<T>());
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(order_3_der, T, all_float_types)
{
    //RandomSample<T>                   rng{1, 2};
    T x = 1; //rng.next();
    T y = 2; //rng.next();

    rvar<T, 3>                        x_ad = x;
    rvar<T, 3>                        y_ad = y;

    T                                 fv   = f(x, y);
    rvar<T, 3>                        f_ad = f(x_ad, y_ad);
    gradient_tape<T, 3, BOOST_MATH_BUFFER_SIZE>& tape = get_active_tape<T, 3>();
    tape.zero_grad();

    auto                                     df3     = df_3_a(x, y);
    auto                                     grad_ad = grad(f_ad, &x_ad, &y_ad);

    std::vector<std::vector<std::vector<T>>> grad_tensor;
    for (int i = 0; i < 2; i++) {
        auto df_hess = hess(grad_ad[i], &x_ad, &y_ad);
        grad_tensor.push_back(df_hess);
    }
    auto grad_nd_func_test = grad_nd<3>(f_ad, &x_ad, &y_ad);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                BOOST_REQUIRE_CLOSE(grad_tensor[i][j][k], df3[i][j][k], boost_close_tol<T>());
                BOOST_REQUIRE_CLOSE(grad_nd_func_test[i][j][k], df3[i][j][k], boost_close_tol<T>());
            }
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(fourth_derivative, T, all_float_types)
{
    //RandomSample<T> rng{1, 2};
    // this derivative is a bit unstable and so the tests
    // are difficult to get to pass for floats and doubles (need long doubles)
    // but its stable enough at x=1 y=2
    T x_val = 1; //rng.next();
    T y_val = 2; //rng.next();

    rvar<T, 4>                                                    x_ad  = x_val;
    rvar<T, 4>                                                    y_ad  = y_val;
    rvar<T, 4>                                                    f_ad  = f(x_ad, y_ad);

    auto                                                          df4   = df_4_a(x_val, y_val);
    auto                                                          gf    = grad(f_ad, &x_ad, &y_ad);
    std::array<std::array<std::array<std::array<T, 2>, 2>, 2>, 2> ggggf;
    for (int i = 0; i < 2; ++i) {
        auto hess1 = grad(gf[i], &x_ad, &y_ad);
        for (int j = 0; j < 2; ++j) {
            auto hess2 = grad(hess1[j], &x_ad, &y_ad);
            for (int k = 0; k < 2; ++k) {
                auto hess3 = grad(hess2[k], &x_ad, &y_ad);
                for (int l = 0; l < 2; ++l) {
                    ggggf[i][j][k][l] = hess3[l];
                }
            }
        }
    }
    auto grad_nd_func_test = grad_nd<4>(f_ad, &x_ad, &y_ad);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                for (int l = 0; l < 2; l++) {
                    BOOST_REQUIRE_CLOSE(ggggf[i][j][k][l], df4[i][j][k][l], boost_close_tol<T>());
                    BOOST_REQUIRE_CLOSE(grad_nd_func_test[i][j][k][l],
                                        df4[i][j][k][l],
                                        boost_close_tol<T>());
                }
            }
        }
    }
}
BOOST_AUTO_TEST_SUITE_END()
