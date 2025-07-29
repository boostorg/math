#include "test_autodiff_reverse.hpp"

BOOST_AUTO_TEST_SUITE(explicit_rvar_constructors)

BOOST_AUTO_TEST_CASE_TEMPLATE(rvar_constructors_and_utils, T, all_float_types)
{
    RandomSample<T>             rng{-1, 1};
    /* raw constructors */
    T                                        x_value = rng.next();
    rdiff::rvar<T, 1>                        x0(x_value);
    rdiff::rvar<T, 2>                        x1(x_value);
    rdiff::rvar<T, 3>                        x2(x_value);

    /* raw eval */
    BOOST_CHECK_EQUAL(x_value, x0.evaluate());
    BOOST_CHECK_EQUAL(x_value, x1.evaluate().evaluate());
    BOOST_CHECK_EQUAL(x_value, x2.evaluate().evaluate().evaluate());

    auto v = x1.item();
    /* get item helper */
    BOOST_CHECK_EQUAL(x_value, x0.item());
    BOOST_CHECK_EQUAL(x_value, x1.item());
    BOOST_CHECK_EQUAL(x_value, x2.item());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(make_rvar_constructors, T, all_float_types)
{
    RandomSample<T> rng{-1, 1};
    T               x_value = rng.next();

    rdiff::rvar<T,1> x1 = rdiff::make_rvar<T,1>(x_value);
    rdiff::rvar<T,2> x2 = rdiff::make_rvar<T,2>(x_value);
    rdiff::rvar<T,3> x3 = rdiff::make_rvar<T,3>(x_value);
    BOOST_CHECK_EQUAL(x_value, x1.item());
    BOOST_CHECK_EQUAL(x_value, x2.item());
    BOOST_CHECK_EQUAL(x_value, x3.item());

    T result = x_value+x_value;
    auto z1 = x1+x1;
    auto z2 = x2+x2;
    auto z3 = x3+x3;

    rdiff::rvar<T,1> y1 = rdiff::make_rvar(z1);
    rdiff::rvar<T,2> y2 = rdiff::make_rvar(z2);
    rdiff::rvar<T,3> y3 = rdiff::make_rvar(z3);

    BOOST_CHECK_EQUAL(y1.item(),result);
    BOOST_CHECK_EQUAL(y2.item(),result);
    BOOST_CHECK_EQUAL(y3.item(),result);
}

BOOST_AUTO_TEST_SUITE_END()
