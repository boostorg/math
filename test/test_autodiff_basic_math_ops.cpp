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
    return x[0]+x[1]+x[2]+x[3]+x[4];
}
BOOST_AUTO_TEST_CASE_TEMPLATE(addition, T, all_float_types)
{
    RandomSample<T> rng{-100, 100};
    T x1_v = rng.next();
    T x2_v = rng.next();
    T test_rvar_p_rvar_v = x1_v+x2_v;

    rvar<T,1> x1 = x1_v;
    rvar<T,1> x2 = x2_v;

    rvar<T,1> test_rvar_p_rvar = x1+x2;
    rvar<T,1> test_rvar_p_float = x1 + x2_v;
    rvar<T,1> test_float_p_rvar = x1_v + x2;
    BOOST_REQUIRE_EQUAL(test_rvar_p_rvar_v, test_rvar_p_rvar.item());
    BOOST_REQUIRE_EQUAL(test_rvar_p_rvar_v, test_rvar_p_float.item());
    BOOST_REQUIRE_EQUAL(test_float_p_rvar.item(), test_rvar_p_rvar_v);

    gradient_tape<T,1>& tape = get_active_tape<T,1>();
    tape.zero_grad();
    test_rvar_p_rvar.backward();
    BOOST_REQUIRE_EQUAL(x1.adjoint(), 1.0);
    BOOST_REQUIRE_EQUAL(x2.adjoint(), 1.0);

    tape.zero_grad();
    rvar<T,1> z = x1+x1+x1+x1+x1;
    z.backward();
    BOOST_REQUIRE_EQUAL(x1.adjoint(), 5.0);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(multiplication, T, all_float_types)
{
    RandomSample<T> rng{-100, 100};
    T x1_v = rng.next();
    T x2_v = rng.next();
    T test_rvar_p_rvar_v = x1_v*x2_v;

    rvar<T,1> x1 = x1_v;
    rvar<T,1> x2 = x2_v;

    rvar<T,1> test_rvar_p_rvar = x1*x2;
    rvar<T,1> test_rvar_p_float = x1 * x2_v;
    rvar<T,1> test_float_p_rvar = x1_v * x2;

    BOOST_REQUIRE_EQUAL(test_rvar_p_rvar_v, test_rvar_p_rvar.item());
    BOOST_REQUIRE_EQUAL(test_rvar_p_rvar_v, test_rvar_p_float.item());
    BOOST_REQUIRE_EQUAL(test_float_p_rvar.item(), test_rvar_p_rvar_v);

    gradient_tape<T,1>& tape = get_active_tape<T,1>();
    tape.zero_grad();
    test_rvar_p_rvar.backward();
    BOOST_REQUIRE_EQUAL(x1.adjoint(), x2_v);
    BOOST_REQUIRE_EQUAL(x2.adjoint(), x1_v);

    tape.zero_grad();
    auto       z = x1 * x1 * x1 * x1 * x1;
    rvar<T, 1> zz(z);
    zz.backward();
    BOOST_REQUIRE_CLOSE(x1.adjoint(),
                        5 * x1_v * x1_v * x1_v * x1_v,
                        300 * std::numeric_limits<T>::epsilon());
    tape.clear();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(minus, T, all_float_types)
{
    RandomSample<T> rng{-100, 100};
    T               x1_v               = rng.next();
    T               x2_v               = rng.next();
    T               test_rvar_p_rvar_v = x1_v - x2_v;

    rvar<T, 1> x1 = x1_v;
    rvar<T, 1> x2 = x2_v;

    rvar<T, 1> test_rvar_p_rvar  = x1 - x2;
    rvar<T, 1> test_rvar_p_float = x1 - x2_v;
    rvar<T, 1> test_float_p_rvar = x1_v - x2;

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
    BOOST_REQUIRE_CLOSE(x1.adjoint(), -5, 300 * std::numeric_limits<T>::epsilon());
    tape.clear();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(division, T, all_float_types)
{
    RandomSample<T> rng{-1, 1};
    T               x1_v               = rng.next();
    T               x2_v               = rng.next();
    T               test_rvar_p_rvar_v = x1_v / x2_v;

    rvar<T, 1> x1 = x1_v;
    rvar<T, 1> x2 = x2_v;

    rvar<T, 1> test_rvar_p_rvar  = x1 / x2;
    rvar<T, 1> test_rvar_p_float = x1 / x2_v;
    rvar<T, 1> test_float_p_rvar = x1_v / x2;

    BOOST_REQUIRE_EQUAL(test_rvar_p_rvar_v, test_rvar_p_rvar.item());
    BOOST_REQUIRE_EQUAL(test_rvar_p_rvar_v, test_rvar_p_float.item());
    BOOST_REQUIRE_EQUAL(test_float_p_rvar.item(), test_rvar_p_rvar_v);

    gradient_tape<T, 1>& tape = get_active_tape<T, 1>();
    tape.zero_grad();
    test_rvar_p_rvar.backward();
    BOOST_REQUIRE_CLOSE(x1.adjoint(), 1.0 / x2_v, 300 * std::numeric_limits<T>::epsilon());
    BOOST_REQUIRE_CLOSE(x2.adjoint(),
                        -x1_v / (x2_v * x2_v),
                        300 * std::numeric_limits<T>::epsilon());

    tape.zero_grad();
    auto       z = T(2.0) / x1 / x1 / x1 / x1 / x1;
    rvar<T, 1> zz(z);
    zz.backward();
    BOOST_REQUIRE_CLOSE(x1.adjoint(),
                        -10 / (x1_v * x1_v * x1_v * x1_v * x1_v * x1_v),
                        300 * std::numeric_limits<T>::epsilon());
    tape.clear();
}

template<typename T>
T derivative_test_function(T x[6])
{
    auto y1 = 1.0 / (x[0] + x[1] * x[2]);
    auto y2 = 2.0 / (x[3] - x[4] / x[5] * x[5]);
    T    z  = y1 * y2;
    return z;
}
BOOST_AUTO_TEST_SUITE_END()

