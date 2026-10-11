//           Copyright Maksym Zhelyeznyakov 2025-2026.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)
#include "test_autodiff_reverse.hpp"
#include <boost/math/tools/workaround.hpp>
#include <cmath>
#include <vector>


using namespace rdiff;
template <class T>
void test_fabs()
{
    RandomSample<T> rng_neg{-10, 0};
    T               x_neg_v         = rng_neg.next();
    T               test_fabs_neg_v = fabs(x_neg_v);

    RandomSample<T> rng_pos{0, 10};
    T               x_pos_v         = rng_pos.next();
    T               test_fabs_pos_v = fabs(x_pos_v);

    rvar<T, 1>      x_neg           = x_neg_v;
    rvar<T, 1>      x_pos           = x_pos_v;

    rvar<T, 1>      fabs_neg        = fabs(x_neg);
    rvar<T, 1>      fabs_pos        = fabs(x_pos);

    CHECK_RELATIVE_ERROR(fabs_neg.item(), test_fabs_neg_v, (boost_close_tol<T>()) / 100);
    CHECK_RELATIVE_ERROR(fabs_pos.item(), test_fabs_pos_v, (boost_close_tol<T>()) / 100);

    gradient_tape<T, 1>& tape = get_active_tape<T, 1>();

    tape.zero_grad();
    fabs_neg.backward();
    CHECK_RELATIVE_ERROR(x_neg.adjoint(), -1, (boost_close_tol<T>()) / 100);
    ;

    tape.zero_grad();
    fabs_pos.backward();
    CHECK_RELATIVE_ERROR(x_pos.adjoint(), 1, (boost_close_tol<T>()) / 100);

    rvar<T, 1> z = fabs(x_neg * x_neg * x_neg);
    tape.zero_grad();
    z.backward();
    T expected_grad = 3 * x_neg_v * x_neg_v * (-1);
    CHECK_RELATIVE_ERROR(x_neg.adjoint(), expected_grad, (boost_close_tol<T>()) / 100);

    tape.clear();
}

template <class T>
void test_abs()
{
    using std::abs;
    RandomSample<T> rng_neg{-1, 0};
    T               x_neg_v        = rng_neg.next();
    T               test_abs_neg_v = fabs(x_neg_v);

    RandomSample<T> rng_pos{0, 1};
    T               x_pos_v        = rng_pos.next();
    T               test_abs_pos_v = abs(x_pos_v);

    rvar<T, 1>      x_neg          = x_neg_v;
    rvar<T, 1>      x_pos          = x_pos_v;

    rvar<T, 1>      abs_neg        = abs(x_neg);
    rvar<T, 1>      abs_pos        = abs(x_pos);

    CHECK_RELATIVE_ERROR(abs_neg.item(), test_abs_neg_v, (boost_close_tol<T>()) / 100);
    CHECK_RELATIVE_ERROR(abs_pos.item(), test_abs_pos_v, (boost_close_tol<T>()) / 100);

    gradient_tape<T, 1>& tape = get_active_tape<T, 1>();

    tape.zero_grad();
    abs_neg.backward();
    CHECK_RELATIVE_ERROR(x_neg.adjoint(), -1, (boost_close_tol<T>()) / 100);

    tape.zero_grad();
    abs_pos.backward();
    CHECK_RELATIVE_ERROR(x_pos.adjoint(), 1, (boost_close_tol<T>()) / 100);

    tape.zero_grad();
    rvar<T, 1> z = abs(x_neg * x_neg * x_neg);
    z.backward();

    T expected_grad = 3 * x_neg_v * x_neg_v * -1;
    CHECK_RELATIVE_ERROR(x_neg.adjoint(), expected_grad, (boost_close_tol<T>()) / 100);

    tape.clear();
}

template <class T>
void test_ceil()
{
    RandomSample<T> rng{-1, 1};
    T               x_v         = rng.next();
    T               test_ceil_v = ceil(x_v);

    rvar<T, 1>      x_rvar      = x_v;

    rvar<T, 1>      x_ceil      = ceil(x_rvar);

    CHECK_RELATIVE_ERROR(x_ceil.item(), test_ceil_v, (boost_close_tol<T>()) / 100);

    gradient_tape<T, 1>& tape = get_active_tape<T, 1>();

    tape.zero_grad();
    x_ceil.backward();
    CHECK_RELATIVE_ERROR(x_rvar.adjoint(), 0.0, (boost_close_tol<T>()) / 100);

    tape.clear();
}
template <class T>
void test_floor()
{
    RandomSample<T> rng{-1, 1};
    T               x_v         = rng.next();
    T               test_ceil_v = floor(x_v);

    rvar<T, 1>      x_rvar      = x_v;

    rvar<T, 1>      x_floor     = floor(x_rvar);

    CHECK_RELATIVE_ERROR(x_floor.item(), test_ceil_v, (boost_close_tol<T>()) / 100);

    gradient_tape<T, 1>& tape = get_active_tape<T, 1>();

    tape.zero_grad();
    x_floor.backward();
    CHECK_RELATIVE_ERROR(x_rvar.adjoint(), 0.0, (boost_close_tol<T>()) / 100);

    tape.clear();
}

template <class T>
void test_exp()
{
    RandomSample<T> rng{-1, 1};
    T               x_v = rng.next();
    rng.next();
    T          test_func_v = exp(x_v);

    rvar<T, 1> x_rvar      = x_v;

    rvar<T, 1> x_func      = exp(x_rvar);

    CHECK_RELATIVE_ERROR(x_func.item(), test_func_v, (boost_close_tol<T>()) / 100);

    gradient_tape<T, 1>& tape = get_active_tape<T, 1>();

    tape.zero_grad();
    rvar<T, 1> test_func_2    = exp(-x_rvar * x_rvar);
    T          expected_deriv = -2.0 * x_v * exp(-x_v * x_v);
    test_func_2.backward();

    CHECK_RELATIVE_ERROR(x_rvar.adjoint(), expected_deriv, (boost_close_tol<T>()) / 100);

    tape.clear();
}

template <class T>
void test_pow()
{
    RandomSample<T> base_rng{0, 2};
    RandomSample<T> exp_rng{-5, 5};
    T               x_v = 2.0;
    base_rng.next();
    T y_v = 0.5;
    exp_rng.next();

    rvar<T, 1> x_rvar     = x_v;
    rvar<T, 1> y_rvar     = y_v;

    T          expected_v = pow(x_v, y_v);
    rvar<T, 1> v_func     = pow(x_rvar, y_rvar);

    CHECK_RELATIVE_ERROR(expected_v, v_func.item(), (boost_close_tol<T>()) / 100);

    gradient_tape<T, 1>& tape = get_active_tape<T, 1>();
    tape.zero_grad();
    v_func.backward();

    T df_dx_expected = y_v * pow(x_v, y_v - 1);
    T df_dy_expected = pow(x_v, y_v) * log(x_v);
    CHECK_RELATIVE_ERROR(df_dx_expected, x_rvar.adjoint(), (boost_close_tol<T>()) / 100);
    CHECK_RELATIVE_ERROR(df_dy_expected, y_rvar.adjoint(), (boost_close_tol<T>()) / 100);

    rvar<T, 1> v_rvar_by_float = pow(x_rvar, y_v);
    rvar<T, 1> v_float_by_rvar = pow(x_v, y_rvar);

    CHECK_RELATIVE_ERROR(expected_v, v_rvar_by_float.item(), (boost_close_tol<T>()) / 100);
    CHECK_RELATIVE_ERROR(expected_v, v_float_by_rvar.item(), (boost_close_tol<T>()) / 100);

    tape.zero_grad();
    v_rvar_by_float.backward();
    CHECK_RELATIVE_ERROR(df_dx_expected, x_rvar.adjoint(), (boost_close_tol<T>()) / 100);
    CHECK_RELATIVE_ERROR(0.0, y_rvar.adjoint(), (boost_close_tol<T>()) / 100);

    tape.zero_grad();
    v_float_by_rvar.backward();
    CHECK_RELATIVE_ERROR(0.0, x_rvar.adjoint(), (boost_close_tol<T>()) / 100);
    CHECK_RELATIVE_ERROR(df_dy_expected, y_rvar.adjoint(), (boost_close_tol<T>()) / 100);
}

template <class T>
void test_log()
{
    RandomSample<T> rng{0.0001, 2};
    T               x_v         = rng.next();
    T               test_func_v = log(x_v);

    rvar<T, 1>      x_rvar      = x_v;

    rvar<T, 1>      x_func      = log(x_rvar);

    CHECK_RELATIVE_ERROR(x_func.item(), test_func_v, (boost_close_tol<T>()) / 100);

    gradient_tape<T, 1>& tape = get_active_tape<T, 1>();

    tape.zero_grad();
    rvar<T, 1> test_func_2    = log(x_rvar * x_rvar);
    T          expected_deriv = 2.0 / x_v;
    test_func_2.backward();

    CHECK_RELATIVE_ERROR(x_rvar.adjoint(), expected_deriv, (boost_close_tol<T>()) / 100);

    tape.clear();
}
template <class T>
void test_sqrt()
{
    RandomSample<T> rng{0.0001, 2};
    T               x_v         = rng.next();
    T               test_func_v = sqrt(x_v);

    rvar<T, 1>      x_rvar      = x_v;
    rvar<T, 1>      x_func      = sqrt(x_rvar);

    CHECK_RELATIVE_ERROR(x_func.item(), test_func_v, (boost_close_tol<T>()) / 100);

    gradient_tape<T, 1>& tape = get_active_tape<T, 1>();

    tape.zero_grad();
    rvar<T, 1> test_func_2    = sqrt(x_rvar * x_rvar * x_rvar);
    T          expected_deriv = 3.0 * sqrt(x_v) / 2.0;
    test_func_2.backward();

    CHECK_RELATIVE_ERROR(x_rvar.adjoint(), expected_deriv, (boost_close_tol<T>()) / 100);

    tape.clear();
}

template <class T>
void test_frexp()
{
    /** @brief frexp(x) returns mantissa in
     *  x = m * 2 ^ e
     *
     *  its implemented as x/2^e, and technically has
     *  a derivative except at finitely many
     *  discontinuous points. In reverse mode autodiff
     *  its treated as a piecewise linear funcfion.
     */
    RandomSample<T> rng{-1, 1};
    T               x_v = rng.next();

    int             i1, i2;
    T               test_func_v = frexp(x_v, &i1);

    rvar<T, 1>      x_rvar      = x_v;
    rvar<T, 1>      x_func      = frexp(x_rvar, &i2);

    CHECK_RELATIVE_ERROR(x_func.item(), test_func_v, (boost_close_tol<T>()) / 100);

    gradient_tape<T, 1>& tape = get_active_tape<T, 1>();

    tape.zero_grad();
    x_func.backward();

    // int i3, i4;
    int i4;
    frexp(x_v, &i4);
    // T   h              = 0.001;
    // T   f1             = frexp(x_v + h, &i3);
    // T   f2             = frexp(x_v, &i4);
    T   expected_deriv = 1.0 / pow(2.0, i4);
    CHECK_RELATIVE_ERROR(x_rvar.adjoint(), expected_deriv, (boost_close_tol<T>()) / 100);

    tape.clear();
}

template <class T>
void test_cos()
{
    RandomSample<T> rng{-10, 10};
    T               x_v         = rng.next();
    T               test_func_v = cos(x_v);

    rvar<T, 1>      x_rvar      = x_v;

    rvar<T, 1>      x_func      = cos(x_rvar);

    CHECK_RELATIVE_ERROR(x_func.item(), test_func_v, (boost_close_tol<T>()) / 100);

    gradient_tape<T, 1>& tape = get_active_tape<T, 1>();

    tape.zero_grad();
    rvar<T, 1> test_func_2    = cos(x_rvar * x_rvar);
    T          expected_deriv = -2 * x_v * sin(x_v * x_v);
    test_func_2.backward();

    CHECK_RELATIVE_ERROR(x_rvar.adjoint(), expected_deriv, (boost_close_tol<T>()) / 100);

    tape.clear();
}

template <class T>
void test_sin()
{
    RandomSample<T> rng{-10, 10};
    T               x_v         = rng.next();
    T               test_func_v = sin(x_v);

    rvar<T, 1>      x_rvar      = x_v;

    rvar<T, 1>      x_func      = sin(x_rvar);

    CHECK_RELATIVE_ERROR(x_func.item(), test_func_v, (boost_close_tol<T>()) / 100);

    gradient_tape<T, 1>& tape = get_active_tape<T, 1>();

    tape.zero_grad();
    rvar<T, 1> test_func_2    = sin(x_rvar * x_rvar);
    T          expected_deriv = 2 * x_v * cos(x_v * x_v);
    test_func_2.backward();

    CHECK_RELATIVE_ERROR(x_rvar.adjoint(), expected_deriv, (boost_close_tol<T>()) / 100);

    tape.clear();
}

template <class T>
void test_tan()
{
    RandomSample<T> rng{-10, 10};
    T               x_v         = rng.next();
    T               test_func_v = tan(x_v);

    rvar<T, 1>      x_rvar      = x_v;

    rvar<T, 1>      x_func      = tan(x_rvar);

    CHECK_RELATIVE_ERROR(x_func.item(), test_func_v, (boost_close_tol<T>()) / 100);

    gradient_tape<T, 1>& tape = get_active_tape<T, 1>();

    tape.zero_grad();
    rvar<T, 1> test_func_2    = tan(x_rvar * x_rvar);
    T          expected_deriv = 2 * x_v * 1 / (cos(x_v * x_v) * cos(x_v * x_v));
    test_func_2.backward();

    CHECK_RELATIVE_ERROR(x_rvar.adjoint(), expected_deriv, (boost_close_tol<T>()) / 100);

    tape.clear();
}

template <class T>
void test_acos()
{
    RandomSample<T> rng{-1, 1};
    T               x_v         = rng.next();
    T               test_func_v = acos(x_v);

    rvar<T, 1>      x_rvar      = x_v;

    rvar<T, 1>      x_func      = acos(x_rvar);

    CHECK_RELATIVE_ERROR(x_func.item(), test_func_v, (boost_close_tol<T>()) / 100);

    gradient_tape<T, 1>& tape = get_active_tape<T, 1>();

    tape.zero_grad();
    rvar<T, 1> test_func_2    = acos(x_rvar * x_rvar);
    T          expected_deriv = -2 * x_v / sqrt(1 - x_v * x_v * x_v * x_v);
    test_func_2.backward();

    CHECK_RELATIVE_ERROR(x_rvar.adjoint(), expected_deriv, (boost_close_tol<T>()) / 100);

    tape.clear();
}

template <class T>
void test_asin()
{
    RandomSample<T> rng{-1, 1};
    T               x_v         = rng.next();
    T               test_func_v = asin(x_v);

    rvar<T, 1>      x_rvar      = x_v;

    rvar<T, 1>      x_func      = asin(x_rvar);

    CHECK_RELATIVE_ERROR(x_func.item(), test_func_v, (boost_close_tol<T>()) / 100);

    gradient_tape<T, 1>& tape = get_active_tape<T, 1>();

    tape.zero_grad();
    rvar<T, 1> test_func_2    = asin(x_rvar * x_rvar);
    T          expected_deriv = 2 * x_v / sqrt(1 - x_v * x_v * x_v * x_v);
    test_func_2.backward();

    CHECK_RELATIVE_ERROR(x_rvar.adjoint(), expected_deriv, (boost_close_tol<T>()) / 100);

    tape.clear();
}

template <class T>
void test_atan()
{
    RandomSample<T> rng{-1, 1};
    T               x_v         = rng.next();
    T               test_func_v = atan(x_v);

    rvar<T, 1>      x_rvar      = x_v;

    rvar<T, 1>      x_func      = atan(x_rvar);

    CHECK_RELATIVE_ERROR(x_func.item(), test_func_v, (boost_close_tol<T>()) / 100);

    gradient_tape<T, 1>& tape = get_active_tape<T, 1>();

    tape.zero_grad();
    rvar<T, 1> test_func_2    = atan(x_rvar * x_rvar);
    T          expected_deriv = 2 * x_v / (1 + x_v * x_v * x_v * x_v);
    test_func_2.backward();

    CHECK_RELATIVE_ERROR(x_rvar.adjoint(), expected_deriv, (boost_close_tol<T>()) / 100);

    tape.clear();
}
template <class T>
void test_atan2()
{
    RandomSample<T> rng{-10, 10};

    T               y_v      = rng.next();
    T               x_v      = rng.next();

    T               expected = atan2(y_v, x_v);

    rvar<T, 1>      y_rvar   = y_v;
    rvar<T, 1>      x_rvar   = x_v;

    {
        rvar<T, 1> z = atan2(y_rvar, x_rvar);
        CHECK_RELATIVE_ERROR(z.item(), expected, (boost_close_tol<T>()) / 100);

        gradient_tape<T, 1>& tape = get_active_tape<T, 1>();
        tape.zero_grad();
        z.backward();

        T denom       = x_v * x_v + y_v * y_v;
        T dy_expected = x_v / denom;
        T dx_expected = -y_v / denom;

        CHECK_RELATIVE_ERROR(y_rvar.adjoint(), dy_expected, (boost_close_tol<T>()) / 100);
        CHECK_RELATIVE_ERROR(x_rvar.adjoint(), dx_expected, (boost_close_tol<T>()) / 100);
    }

    {
        rvar<T, 1> y2 = y_v;
        rvar<T, 1> z  = atan2(y2, x_v);
        CHECK_RELATIVE_ERROR(z.item(), expected, (boost_close_tol<T>()) / 100);

        gradient_tape<T, 1>& tape = get_active_tape<T, 1>();
        tape.zero_grad();
        z.backward();

        T dy_expected = x_v / (x_v * x_v + y_v * y_v);

        CHECK_RELATIVE_ERROR(y2.adjoint(), dy_expected, (boost_close_tol<T>()) / 100);
    }

    {
        rvar<T, 1> x2 = x_v;
        rvar<T, 1> z  = atan2(y_v, x2);
        CHECK_RELATIVE_ERROR(z.item(), expected, (boost_close_tol<T>()) / 100);

        gradient_tape<T, 1>& tape = get_active_tape<T, 1>();
        tape.zero_grad();
        z.backward();

        T dx_expected = -y_v / (x_v * x_v + y_v * y_v);

        CHECK_RELATIVE_ERROR(x2.adjoint(), dx_expected, (boost_close_tol<T>()) / 100);
    }
}

template <class T>
void test_trunc()
{
    RandomSample<T> rng{-1, 1};
    T               x_v          = rng.next();
    T               test_trunc_v = trunc(x_v);

    rvar<T, 1>      x_rvar       = x_v;

    rvar<T, 1>      x_trunc      = trunc(x_rvar);

    CHECK_RELATIVE_ERROR(x_trunc.item(), test_trunc_v, (boost_close_tol<T>()) / 100);

    gradient_tape<T, 1>& tape = get_active_tape<T, 1>();

    tape.zero_grad();
    x_trunc.backward();
    CHECK_RELATIVE_ERROR(x_rvar.adjoint(), 0.0, (boost_close_tol<T>()) / 100);

    tape.clear();
}

template <class T>
void test_fmod()
{
    using boost::math::tools::fmod_workaround;
    RandomSample<T> rng_x{-1, 1};
    RandomSample<T> rng_y{-1, 1};

    T               x_v = rng_x.next();
    T               y_v = rng_y.next();
    while (y_v == 0.0) {
        y_v = rng_y.next();
    }

    T expected_fmod = fmod(x_v, y_v);

    rvar<T, 1> x_rvar        = x_v;
    rvar<T, 1> y_rvar        = y_v;

    rvar<T, 1> fmod_rvar     = fmod(x_rvar, y_rvar);
    CHECK_RELATIVE_ERROR(fmod_rvar.item(), expected_fmod, (boost_close_tol<T>()) / 100);

    gradient_tape<T, 1>& tape = get_active_tape<T, 1>();
    tape.zero_grad();
    fmod_rvar.backward();

    T dx_expected = 1.0;
    T dy_expected = -1.0 * trunc(x_v / y_v);

    CHECK_RELATIVE_ERROR(x_rvar.adjoint(), dx_expected, (boost_close_tol<T>()) / 100);
    CHECK_RELATIVE_ERROR(y_rvar.adjoint(), dy_expected, (boost_close_tol<T>()) / 100);

    rvar<T, 1> fmod_rvar_float = fmod(x_rvar, y_v);
    CHECK_RELATIVE_ERROR(fmod_rvar_float.item(), expected_fmod, (boost_close_tol<T>()) / 100);

    tape.zero_grad();
    fmod_rvar_float.backward();
    CHECK_RELATIVE_ERROR(x_rvar.adjoint(), dx_expected, (boost_close_tol<T>()) / 100);

    rvar<T, 1> fmod_float_rvar = fmod(x_v, y_rvar);
    CHECK_RELATIVE_ERROR(fmod_float_rvar.item(), expected_fmod, (boost_close_tol<T>()) / 100);

    tape.zero_grad();
    fmod_float_rvar.backward();
    CHECK_RELATIVE_ERROR(y_rvar.adjoint(), dy_expected, (boost_close_tol<T>()) / 100);

    tape.clear();
}

template <class T>
void test_round()
{
    RandomSample<T> rng{-1, 1};
    T               x_v          = rng.next();
    T               test_round_v = round(x_v);

    rvar<T, 1>      x_rvar       = x_v;

    rvar<T, 1>      x_round      = round(x_rvar);

    CHECK_RELATIVE_ERROR(x_round.item(), test_round_v, (boost_close_tol<T>()) / 100);

    gradient_tape<T, 1>& tape = get_active_tape<T, 1>();

    tape.zero_grad();
    x_round.backward();
    CHECK_RELATIVE_ERROR(x_rvar.adjoint(), 0.0, (boost_close_tol<T>()) / 100);

    tape.clear();
}
template <class T>
void test_iround()
{
    using boost::math::iround;
    RandomSample<T> rng{-1, 1};
    T               x_v          = rng.next();
    int             test_round_v = iround(x_v);
    rvar<T, 1>      x_rvar       = x_v;
    int             x_round      = iround(x_rvar);
    CHECK_EQUAL(x_round, test_round_v);
}

template <class T>
void test_lround()
{
    using boost::math::lround;
    RandomSample<T> rng{-1, 1};
    T               x_v          = rng.next();
    long            test_round_v = lround(x_v);
    rvar<T, 1>      x_rvar       = x_v;
    long            x_round      = lround(x_rvar);
    CHECK_EQUAL(x_round, test_round_v);
}

template <class T>
void test_llround()
{
    using boost::math::llround;
    RandomSample<T> rng{-1, 1};
    T               x_v          = rng.next();
    long long       test_round_v = llround(x_v);
    rvar<T, 1>      x_rvar       = x_v;
    long long       x_round      = llround(x_rvar);
    CHECK_EQUAL(x_round, test_round_v);
}

template <class T>
void test_itrunc()
{
    using boost::math::itrunc;
    RandomSample<T> rng{-1, 1};
    T               x_v          = rng.next();
    int             test_trunc_v = itrunc(x_v);
    rvar<T, 1>      x_rvar       = x_v;
    int             x_trunc      = itrunc(x_rvar);
    CHECK_EQUAL(x_trunc, test_trunc_v);
}
template <class T>
void test_ltrunc()
{
    using boost::math::ltrunc;
    RandomSample<T> rng{-1, 1};
    T               x_v          = rng.next();
    long            test_trunc_v = ltrunc(x_v);
    rvar<T, 1>      x_rvar       = x_v;
    long            x_trunc      = ltrunc(x_rvar);
    CHECK_EQUAL(x_trunc, test_trunc_v);
}
template <class T>
void test_lltrunc()
{
    using boost::math::lltrunc;
    RandomSample<T> rng{-1, 1};
    T               x_v          = rng.next();
    int             test_trunc_v = lltrunc(x_v);
    rvar<T, 1>      x_rvar       = x_v;
    int             x_trunc      = lltrunc(x_rvar);
    CHECK_EQUAL(x_trunc, test_trunc_v);
}

template <class T>
void test_sinh()
{
    RandomSample<T> rng{-1, 1};
    T               x_v         = rng.next();
    T               test_func_v = sinh(x_v);

    rvar<T, 1>      x_rvar      = x_v;

    rvar<T, 1>      x_func      = sinh(x_rvar);

    CHECK_RELATIVE_ERROR(x_func.item(), test_func_v, (boost_close_tol<T>()) / 100);

    gradient_tape<T, 1>& tape = get_active_tape<T, 1>();

    tape.zero_grad();
    rvar<T, 1> test_func_2    = sinh(x_rvar * x_rvar);
    T          expected_deriv = 2 * x_v * cosh(x_v * x_v);
    test_func_2.backward();

    CHECK_RELATIVE_ERROR(x_rvar.adjoint(), expected_deriv, (boost_close_tol<T>()) / 100);

    tape.clear();
}

template <class T>
void test_cosh()
{
    RandomSample<T> rng{-1, 1};
    T               x_v         = rng.next();
    T               test_func_v = cosh(x_v);

    rvar<T, 1>      x_rvar      = x_v;

    rvar<T, 1>      x_func      = cosh(x_rvar);

    CHECK_RELATIVE_ERROR(x_func.item(), test_func_v, (boost_close_tol<T>()) / 100);

    gradient_tape<T, 1>& tape = get_active_tape<T, 1>();

    tape.zero_grad();
    rvar<T, 1> test_func_2    = cosh(x_rvar * x_rvar);
    T          expected_deriv = 2 * x_v * sinh(x_v * x_v);
    test_func_2.backward();

    CHECK_RELATIVE_ERROR(x_rvar.adjoint(), expected_deriv, (boost_close_tol<T>()) / 100);

    tape.clear();
}

template <class T>
void test_tanh()
{
    RandomSample<T> rng{-1, 1};
    T               x_v         = rng.next();
    T               test_func_v = tanh(x_v);
    rvar<T, 1>      x_rvar      = x_v;
    rvar<T, 1>      x_func      = tanh(x_rvar);
    CHECK_RELATIVE_ERROR(x_func.item(), test_func_v, (boost_close_tol<T>()) / 100);
    gradient_tape<T, 1>& tape = get_active_tape<T, 1>();
    tape.zero_grad();
    rvar<T, 1> test_func_2    = tanh(x_rvar * x_rvar);
    T          expected_deriv = 2 * x_v / (cosh(x_v * x_v) * cosh(x_v * x_v));
    test_func_2.backward();
    CHECK_RELATIVE_ERROR(x_rvar.adjoint(), expected_deriv, (boost_close_tol<T>()) / 100);
    tape.clear();
}
template <class T>
void test_ldexp()
{
    RandomSample<T> rng{-2, 2};
    T               x_v         = rng.next();
    int             exp_i       = 3;
    T               test_func_v = ldexp(x_v, exp_i);

    rvar<T, 1> x_rvar = x_v;

    rvar<T, 1> x_func = ldexp(x_rvar, exp_i);

    CHECK_RELATIVE_ERROR(x_func.item(), test_func_v, (boost_close_tol<T>()) / 100);

    gradient_tape<T, 1>& tape = get_active_tape<T, 1>();

    tape.zero_grad();

    // f(x) = ldexp(x^2, exp_i) = x^2 * 2^exp_i
    rvar<T, 1> test_func_2 = ldexp(x_rvar * x_rvar, exp_i);

    // Expected derivative: d/dx [x^2 * 2^exp_i] = 2x * 2^exp_i
    T expected_deriv = static_cast<T>(2) * x_v * pow(static_cast<T>(2.0), exp_i);

    test_func_2.backward();

    CHECK_RELATIVE_ERROR(x_rvar.adjoint(), expected_deriv, (boost_close_tol<T>()) / 100);

    tape.clear();
}
template <class T>
void test_ldexp_huge_val()
{
    const T   huge_val    = std::numeric_limits<T>::infinity();
    const int exp         = 1;
    T         result_huge = ldexp(huge_val, exp);
    CHECK_EQUAL(result_huge, huge_val);
    rvar<T, 1> x_huge    = huge_val;
    rvar<T, 1> func_huge = ldexp(x_huge, exp);
    CHECK_EQUAL(func_huge.item(), huge_val);
    func_huge.backward();
    T expected_adjoint = static_cast<T>(2);
    CHECK_EQUAL(x_huge.adjoint(), expected_adjoint);
    get_active_tape<T, 1>().clear();
}

template <class T>
void test_log10()
{
    RandomSample<T> rng{0.0001, 10};
    T               x_v         = rng.next();
    T               test_func_v = log10(x_v);
    rvar<T, 1>      x_rvar      = x_v;
    rvar<T, 1>      x_func      = log10(x_rvar);
    CHECK_RELATIVE_ERROR(x_func.item(), test_func_v, (boost_close_tol<T>()) / 100);
    gradient_tape<T, 1>& tape = get_active_tape<T, 1>();
    tape.zero_grad();
    rvar<T, 1> test_func_2    = log10(x_rvar * x_rvar);
    T          expected_deriv = T{2.0} / (x_v * log(T{10.0}));
    test_func_2.backward();
    CHECK_RELATIVE_ERROR(x_rvar.adjoint(), expected_deriv, (boost_close_tol<T>()) / 100);
    tape.clear();
}

template<typename T>
T test_function_1(const T& x, const T& y)
{
    /** @brief test function for
     *  log
     *  abs
     *  exp
     *  pow
     *  sqrt
     */
    auto y1 = log(1 + abs(x)) * exp(y);
    auto y2 = pow(x + y, 2.5);
    auto y3 = sqrt(1 + x * y);
    return y1 + y2 - y3;
}

template<typename T>
std::vector<T> grad_test_func_1_analytical(T x, T y)
{
    T f_x = static_cast<T>(-1.0 / 2.0 * y / sqrt(x * y + 1) + 2.5 * pow(x + y, 1.5)
                           + sqrt(pow(x, 2)) * exp(y) / (x * (sqrt(pow(x, 2)) + 1)));
    T f_y = static_cast<T>(-1.0 / 2.0 * x / sqrt(x * y + 1) + 2.5 * pow(x + y, 1.5)
                           + exp(y) * log(sqrt(pow(x, 2)) + 1));
    return {f_x, f_y};
}

template<typename T>
std::vector<std::vector<T>> grad2_test_func_1_analytical(T x, T y)
{
    T f_xx = static_cast<T>((1.0 / 4.0) * pow(y, 2) / pow(x * y + 1, 3.0 / 2.0) + 3.75 * sqrt(x + y)
                            - exp(y) / pow(sqrt(pow(x, 2)) + 1, 2));
    T f_xy = static_cast<T>((1.0 / 4.0) * x * y / pow(x * y + 1, 3.0 / 2.0) + 3.75 * sqrt(x + y)
                            - (1.0 / 2.0) / sqrt(x * y + 1)
                            + sqrt(pow(x, 2)) * exp(y) / (x * (sqrt(pow(x, 2)) + 1)));
    T f_yx = static_cast<T>((1.0 / 4.0) * x * y / pow(x * y + 1, 3.0 / 2.0) + 3.75 * sqrt(x + y)
                            - (1.0 / 2.0) / sqrt(x * y + 1)
                            + sqrt(pow(x, 2)) * exp(y) / (x * (sqrt(pow(x, 2)) + 1)));
    T f_yy = static_cast<T>((1.0 / 4.0) * pow(x, 2) / pow(x * y + 1, 3.0 / 2.0) + 3.75 * sqrt(x + y)
                            + exp(y) * log(sqrt(pow(x, 2)) + 1));
    return {{f_xx, f_xy}, {f_yx, f_yy}};
}

template<typename T>
std::vector<std::vector<std::vector<T>>> grad3_test_func_1_analytical(T x, T y)
{
    T f_xxx = static_cast<T>(-3.0 / 8.0 * pow(y, 3) / pow(x * y + 1, 5.0 / 2.0)
                             + 1.875 * pow(x + y, -0.5)
                             + 2 * sqrt(pow(x, 2)) * exp(y) / (x * pow(sqrt(pow(x, 2)) + 1, 3)));
    T f_xxy = static_cast<T>(-3.0 / 8.0 * x * pow(y, 2) / pow(x * y + 1, 5.0 / 2.0)
                             + (1.0 / 2.0) * y / pow(x * y + 1, 3.0 / 2.0)
                             + 1.875 * pow(x + y, -0.5) - exp(y) / pow(sqrt(pow(x, 2)) + 1, 2));
    T f_xyx = static_cast<T>(-3.0 / 8.0 * x * pow(y, 2) / pow(x * y + 1, 5.0 / 2.0)
                             + (1.0 / 2.0) * y / pow(x * y + 1, 3.0 / 2.0)
                             + 1.875 * pow(x + y, -0.5) - exp(y) / pow(sqrt(pow(x, 2)) + 1, 2));
    T f_xyy = static_cast<T>(-3.0 / 8.0 * pow(x, 2) * y / pow(x * y + 1, 5.0 / 2.0)
                             + (1.0 / 2.0) * x / pow(x * y + 1, 3.0 / 2.0) + 1.875 * pow(x + y, -0.5)
                             + sqrt(pow(x, 2)) * exp(y) / (x * (sqrt(pow(x, 2)) + 1)));
    T f_yxx = static_cast<T>(-3.0 / 8.0 * x * pow(y, 2) / pow(x * y + 1, 5.0 / 2.0)
                             + (1.0 / 2.0) * y / pow(x * y + 1, 3.0 / 2.0)
                             + 1.875 * pow(x + y, -0.5) - exp(y) / pow(sqrt(pow(x, 2)) + 1, 2));
    T f_yxy = static_cast<T>(-3.0 / 8.0 * pow(x, 2) * y / pow(x * y + 1, 5.0 / 2.0)
                             + (1.0 / 2.0) * x / pow(x * y + 1, 3.0 / 2.0) + 1.875 * pow(x + y, -0.5)
                             + sqrt(pow(x, 2)) * exp(y) / (x * (sqrt(pow(x, 2)) + 1)));
    T f_yyx = static_cast<T>(-3.0 / 8.0 * pow(x, 2) * y / pow(x * y + 1, 5.0 / 2.0)
                             + (1.0 / 2.0) * x / pow(x * y + 1, 3.0 / 2.0) + 1.875 * pow(x + y, -0.5)
                             + sqrt(pow(x, 2)) * exp(y) / (x * (sqrt(pow(x, 2)) + 1)));
    T f_yyy = static_cast<T>(-3.0 / 8.0 * pow(x, 3) / pow(x * y + 1, 5.0 / 2.0)
                             + 1.875 * pow(x + y, -0.5) + exp(y) * log(sqrt(pow(x, 2)) + 1));
    return {{{f_xxx, f_xxy}, {f_xyx, f_xyy}}, {{f_yxx, f_yxy}, {f_yyx, f_yyy}}};
}

template<typename T>
std::vector<std::vector<std::vector<std::vector<T>>>> grad4_test_func_1_analytical(T x, T y)
{
    T f_xxxx = static_cast<T>((15.0 / 16.0) * pow(y, 4) / pow(x * y + 1, 7.0 / 2.0)
                              - 0.9375 * pow(x + y, -1.5)
                              - 6 * exp(y) / pow(sqrt(pow(x, 2)) + 1, 4));
    T f_xxxy = static_cast<T>((15.0 / 16.0) * x * pow(y, 3) / pow(x * y + 1, 7.0 / 2.0)
                              - 9.0 / 8.0 * pow(y, 2) / pow(x * y + 1, 5.0 / 2.0)
                              - 0.9375 * pow(x + y, -1.5)
                              + 2 * sqrt(pow(x, 2)) * exp(y) / (x * pow(sqrt(pow(x, 2)) + 1, 3)));
    T f_xxyx = static_cast<T>((15.0 / 16.0) * x * pow(y, 3) / pow(x * y + 1, 7.0 / 2.0)
                              - 9.0 / 8.0 * pow(y, 2) / pow(x * y + 1, 5.0 / 2.0)
                              - 0.9375 * pow(x + y, -1.5)
                              + 2 * sqrt(pow(x, 2)) * exp(y) / (x * pow(sqrt(pow(x, 2)) + 1, 3)));
    T f_xxyy = static_cast<T>((15.0 / 16.0) * pow(x, 2) * pow(y, 2) / pow(x * y + 1, 7.0 / 2.0)
                              - 3.0 / 2.0 * x * y / pow(x * y + 1, 5.0 / 2.0)
                              - 0.9375 * pow(x + y, -1.5) - exp(y) / pow(sqrt(pow(x, 2)) + 1, 2)
                              + (1.0 / 2.0) / pow(x * y + 1, 3.0 / 2.0));
    T f_xyxx = static_cast<T>((15.0 / 16.0) * x * pow(y, 3) / pow(x * y + 1, 7.0 / 2.0)
                              - 9.0 / 8.0 * pow(y, 2) / pow(x * y + 1, 5.0 / 2.0)
                              - 0.9375 * pow(x + y, -1.5)
                              + 2 * sqrt(pow(x, 2)) * exp(y) / (x * pow(sqrt(pow(x, 2)) + 1, 3)));
    T f_xyxy = static_cast<T>((15.0 / 16.0) * pow(x, 2) * pow(y, 2) / pow(x * y + 1, 7.0 / 2.0)
                              - 3.0 / 2.0 * x * y / pow(x * y + 1, 5.0 / 2.0)
                              - 0.9375 * pow(x + y, -1.5) - exp(y) / pow(sqrt(pow(x, 2)) + 1, 2)
                              + (1.0 / 2.0) / pow(x * y + 1, 3.0 / 2.0));
    T f_xyyx = static_cast<T>((15.0 / 16.0) * pow(x, 2) * pow(y, 2) / pow(x * y + 1, 7.0 / 2.0)
                              - 3.0 / 2.0 * x * y / pow(x * y + 1, 5.0 / 2.0)
                              - 0.9375 * pow(x + y, -1.5) - exp(y) / pow(sqrt(pow(x, 2)) + 1, 2)
                              + (1.0 / 2.0) / pow(x * y + 1, 3.0 / 2.0));
    T f_xyyy = static_cast<T>((15.0 / 16.0) * pow(x, 3) * y / pow(x * y + 1, 7.0 / 2.0)
                              - 9.0 / 8.0 * pow(x, 2) / pow(x * y + 1, 5.0 / 2.0)
                              - 0.9375 * pow(x + y, -1.5)
                              + sqrt(pow(x, 2)) * exp(y) / (x * (sqrt(pow(x, 2)) + 1)));
    T f_yxxx = static_cast<T>((15.0 / 16.0) * x * pow(y, 3) / pow(x * y + 1, 7.0 / 2.0)
                              - 9.0 / 8.0 * pow(y, 2) / pow(x * y + 1, 5.0 / 2.0)
                              - 0.9375 * pow(x + y, -1.5)
                              + 2 * sqrt(pow(x, 2)) * exp(y) / (x * pow(sqrt(pow(x, 2)) + 1, 3)));
    T f_yxxy = static_cast<T>((15.0 / 16.0) * pow(x, 2) * pow(y, 2) / pow(x * y + 1, 7.0 / 2.0)
                              - 3.0 / 2.0 * x * y / pow(x * y + 1, 5.0 / 2.0)
                              - 0.9375 * pow(x + y, -1.5) - exp(y) / pow(sqrt(pow(x, 2)) + 1, 2)
                              + (1.0 / 2.0) / pow(x * y + 1, 3.0 / 2.0));
    T f_yxyx = static_cast<T>((15.0 / 16.0) * pow(x, 2) * pow(y, 2) / pow(x * y + 1, 7.0 / 2.0)
                              - 3.0 / 2.0 * x * y / pow(x * y + 1, 5.0 / 2.0)
                              - 0.9375 * pow(x + y, -1.5) - exp(y) / pow(sqrt(pow(x, 2)) + 1, 2)
                              + (1.0 / 2.0) / pow(x * y + 1, 3.0 / 2.0));
    T f_yxyy = static_cast<T>((15.0 / 16.0) * pow(x, 3) * y / pow(x * y + 1, 7.0 / 2.0)
                              - 9.0 / 8.0 * pow(x, 2) / pow(x * y + 1, 5.0 / 2.0)
                              - 0.9375 * pow(x + y, -1.5)
                              + sqrt(pow(x, 2)) * exp(y) / (x * (sqrt(pow(x, 2)) + 1)));
    T f_yyxx = static_cast<T>((15.0 / 16.0) * pow(x, 2) * pow(y, 2) / pow(x * y + 1, 7.0 / 2.0)
                              - 3.0 / 2.0 * x * y / pow(x * y + 1, 5.0 / 2.0)
                              - 0.9375 * pow(x + y, -1.5) - exp(y) / pow(sqrt(pow(x, 2)) + 1, 2)
                              + (1.0 / 2.0) / pow(x * y + 1, 3.0 / 2.0));
    T f_yyxy = static_cast<T>((15.0 / 16.0) * pow(x, 3) * y / pow(x * y + 1, 7.0 / 2.0)
                              - 9.0 / 8.0 * pow(x, 2) / pow(x * y + 1, 5.0 / 2.0)
                              - 0.9375 * pow(x + y, -1.5)
                              + sqrt(pow(x, 2)) * exp(y) / (x * (sqrt(pow(x, 2)) + 1)));
    T f_yyyx = static_cast<T>((15.0 / 16.0) * pow(x, 3) * y / pow(x * y + 1, 7.0 / 2.0)
                              - 9.0 / 8.0 * pow(x, 2) / pow(x * y + 1, 5.0 / 2.0)
                              - 0.9375 * pow(x + y, -1.5)
                              + sqrt(pow(x, 2)) * exp(y) / (x * (sqrt(pow(x, 2)) + 1)));
    T f_yyyy = static_cast<T>((15.0 / 16.0) * pow(x, 4) / pow(x * y + 1, 7.0 / 2.0)
                              - 0.9375 * pow(x + y, -1.5) + exp(y) * log(sqrt(pow(x, 2)) + 1));
    return {{{{f_xxxx, f_xxxy}, {f_xxyx, f_xxyy}}, {{f_xyxx, f_xyxy}, {f_xyyx, f_xyyy}}},
            {{{f_yxxx, f_yxxy}, {f_yxyx, f_yxyy}}, {{f_yyxx, f_yyxy}, {f_yyyx, f_yyyy}}}};
}

template <class T>
void test_func_1_first_derivative()
{
    RandomSample<T> rng{0.0001, 10};
    T               x                 = rng.next();
    T               y                 = rng.next();

    rvar<T, 1>      x_ad              = x;
    rvar<T, 1>      y_ad              = y;

    //T               fv                = test_function_1(x, y);
    rvar<T, 1>      f_ad              = test_function_1(x_ad, y_ad);
    auto            grad_f_analytical = grad_test_func_1_analytical(x, y);
    /* intended use case */
    f_ad.backward();
    CHECK_RELATIVE_ERROR(x_ad.adjoint(), grad_f_analytical[0], (boost_close_tol<T>()) / 100);
    CHECK_RELATIVE_ERROR(y_ad.adjoint(), grad_f_analytical[1], (boost_close_tol<T>()) / 100);

    gradient_tape<T, 1>& tape = get_active_tape<T, 1>();
    tape.zero_grad();

    /* grad test */
    auto grad_func_test_grad    = grad(f_ad, &x_ad, &y_ad);
    /* grad_nd test */
    auto grad_nd_func_test_grad = grad_nd<1>(f_ad, &x_ad, &y_ad);

    for (int i = 0; i < 2; i++) {
        CHECK_RELATIVE_ERROR(grad_func_test_grad[i], grad_f_analytical[i], (boost_close_tol<T>()) / 100);
        CHECK_RELATIVE_ERROR(grad_nd_func_test_grad[i], grad_f_analytical[i], (boost_close_tol<T>()) / 100);
    }
}
template <class T>
void test_func_1_second_derivative_and_hessian()
{
    RandomSample<T>                   rng{0.0001, 10};
    T                                 x    = rng.next();
    T                                 y    = rng.next();

    rvar<T, 2>                        x_ad = x;
    rvar<T, 2>                        y_ad = y;

    //T                                 fv   = test_function_1(x, y);
    rvar<T, 2>                        f_ad = test_function_1(x_ad, y_ad);
    gradient_tape<T, 2>&              tape = get_active_tape<T, 2>();
    tape.zero_grad();

    auto hess_analytical   = grad2_test_func_1_analytical(x, y);
    auto hess_func_test    = hess(f_ad, &x_ad, &y_ad);
    auto grad_nd_func_test = grad_nd<2>(f_ad, &x_ad, &y_ad);

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            CHECK_RELATIVE_ERROR(hess_func_test[i][j], hess_analytical[i][j], (boost_close_tol<T>()) / 100);
            CHECK_RELATIVE_ERROR(hess_func_test[i][j], grad_nd_func_test[i][j], (boost_close_tol<T>()) / 100);
        }
    }
}

template <class T>
void test_func_1_order_3_der()
{
    RandomSample<T>                   rng{0.0001, 10};
    T                                 x    = rng.next();
    T                                 y    = rng.next();

    rvar<T, 3>                        x_ad = x;
    rvar<T, 3>                        y_ad = y;

    //T                                 fv   = test_function_1(x, y);
    rvar<T, 3>                        f_ad = test_function_1(x_ad, y_ad);
    gradient_tape<T, 3>&              tape = get_active_tape<T, 3>();
    tape.zero_grad();

    auto                                     df3     = grad3_test_func_1_analytical(x, y);
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
                CHECK_RELATIVE_ERROR(grad_tensor[i][j][k], df3[i][j][k], (boost_close_tol<T>()) / 100);
                CHECK_RELATIVE_ERROR(grad_nd_func_test[i][j][k], df3[i][j][k], (boost_close_tol<T>()) / 100);
            }
        }
    }
}

template <class T>
void test_func_1_fourth_derivative()
{
    RandomSample<T> rng{0.0001, 10};

    T               x_val = rng.next();
    T               y_val = rng.next();

    rvar<T, 4>      x_ad  = x_val;
    rvar<T, 4>      y_ad  = y_val;
    rvar<T, 4>      f_ad  = test_function_1(x_ad, y_ad);

    auto            df4   = grad4_test_func_1_analytical(x_val, y_val);
    auto            gf    = grad(f_ad, &x_ad, &y_ad);
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
                    CHECK_RELATIVE_ERROR(ggggf[i][j][k][l],
                                        df4[i][j][k][l], (1 * boost_close_tol<T>()) / 100);
                    CHECK_RELATIVE_ERROR(grad_nd_func_test[i][j][k][l],
                                        df4[i][j][k][l], (1 * boost_close_tol<T>()) / 100);
                }
            }
        }
    }
}

int main()
{
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { test_fabs<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { test_abs<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { test_ceil<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { test_floor<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { test_exp<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { test_pow<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { test_log<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { test_sqrt<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { test_frexp<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { test_cos<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { test_sin<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { test_tan<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { test_acos<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { test_asin<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { test_atan<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { test_atan2<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { test_trunc<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { test_fmod<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { test_round<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { test_iround<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { test_lround<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { test_llround<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { test_itrunc<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { test_ltrunc<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { test_lltrunc<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { test_sinh<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { test_cosh<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { test_tanh<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { test_ldexp<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { test_ldexp_huge_val<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { test_log10<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { test_func_1_first_derivative<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { test_func_1_second_derivative_and_hessian<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { test_func_1_order_3_der<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { test_func_1_fourth_derivative<typename decltype(t)::type>(); });
   return boost::math::test::report_errors();
}
