#include "test_autodiff_reverse.hpp"

BOOST_AUTO_TEST_SUITE(explicit_rvar_constructors)

BOOST_AUTO_TEST_CASE_TEMPLATE(rvar_constructors_and_utils, T, all_float_types)
{
    RandomSample<T>                          rng{-100, 100};
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
    RandomSample<T> rng{-100, 100};
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

BOOST_AUTO_TEST_CASE_TEMPLATE(copy_construction, T, all_float_types)
{
    RandomSample<T> rng{-100, 100};
    T               original_value_1 = rng.next();

    /* copy constructor checks */
    rdiff::rvar<T, 1> original_rvar_1(original_value_1);
    rdiff::rvar<T, 1> copied_rvar_1 = original_rvar_1;

    BOOST_CHECK_EQUAL(copied_rvar_1.item(), original_value_1);

    /* change original, copy stays same */
    T new_value_1   = rng.next();
    original_rvar_1 = new_value_1;
    BOOST_CHECK_EQUAL(copied_rvar_1.item(), original_value_1);
    BOOST_CHECK_EQUAL(original_rvar_1.item(), new_value_1);

    T                 original_value_2 = rng.next();
    rdiff::rvar<T, 2> original_rvar_2(original_value_2);
    rdiff::rvar<T, 2> copied_rvar_2 = original_rvar_2;

    BOOST_CHECK_EQUAL(copied_rvar_2.item(), original_value_2);

    T new_value_2   = rng.next();
    original_rvar_2 = new_value_2;
    BOOST_CHECK_EQUAL(copied_rvar_2.item(), original_value_2);
    BOOST_CHECK_EQUAL(original_rvar_2.item(), new_value_2);

    T                 original_value_3 = rng.next();
    rdiff::rvar<T, 3> original_rvar_3(original_value_3);
    rdiff::rvar<T, 3> copied_rvar_3 = original_rvar_3;

    BOOST_CHECK_EQUAL(copied_rvar_3.item(), original_value_3);

    T new_value_3   = rng.next();
    original_rvar_3 = new_value_3;
    BOOST_CHECK_EQUAL(copied_rvar_3.item(), original_value_3);
    BOOST_CHECK_EQUAL(original_rvar_3.item(), new_value_3);
}
BOOST_AUTO_TEST_CASE_TEMPLATE(copy_assignment, T, all_float_types)
{
    RandomSample<T>   rng{-100, 100};
    T                 source_value_1         = rng.next();
    T                 initial_target_value_1 = rng.next();
    rdiff::rvar<T, 1> source_rvar_1(source_value_1);
    rdiff::rvar<T, 1> target_rvar_1(initial_target_value_1);

    target_rvar_1 = source_rvar_1;

    BOOST_CHECK_EQUAL(target_rvar_1.item(), source_value_1);

    T new_source_value_1 = rng.next();
    source_rvar_1        = new_source_value_1;
    BOOST_CHECK_EQUAL(target_rvar_1.item(), source_value_1);
    BOOST_CHECK_EQUAL(source_rvar_1.item(), new_source_value_1);

    rdiff::rvar<T, 1> self_assign_rvar_1(rng.next());
    T                 value_before_self_assign_1 = self_assign_rvar_1.item();
    self_assign_rvar_1                           = self_assign_rvar_1;
    BOOST_CHECK_EQUAL(self_assign_rvar_1.item(), value_before_self_assign_1);

    T                 source_value_2         = rng.next();
    T                 initial_target_value_2 = rng.next();
    rdiff::rvar<T, 2> source_rvar_2(source_value_2);
    rdiff::rvar<T, 2> target_rvar_2(initial_target_value_2);

    target_rvar_2 = source_rvar_2;

    BOOST_CHECK_EQUAL(target_rvar_2.item(), source_value_2);

    T new_source_value_2 = rng.next();
    source_rvar_2        = new_source_value_2;
    BOOST_CHECK_EQUAL(target_rvar_2.item(), source_value_2);
    BOOST_CHECK_EQUAL(source_rvar_2.item(), new_source_value_2);

    T                 source_value_3         = rng.next();
    T                 initial_target_value_3 = rng.next();
    rdiff::rvar<T, 3> source_rvar_3(source_value_3);
    rdiff::rvar<T, 3> target_rvar_3(initial_target_value_3);

    target_rvar_3 = source_rvar_3;

    BOOST_CHECK_EQUAL(target_rvar_3.item(), source_value_3);

    T new_source_value_3 = rng.next();
    source_rvar_3        = new_source_value_3;
    BOOST_CHECK_EQUAL(target_rvar_3.item(), source_value_3);
    BOOST_CHECK_EQUAL(source_rvar_3.item(), new_source_value_3);
}
BOOST_AUTO_TEST_CASE_TEMPLATE(inplace_addition, T, all_float_types)
{
    using namespace rdiff;
    RandomSample<T> rng{-1, 1};
    T               x1_v     = rng.next();
    T               x2_v     = rng.next();
    T               expected = x1_v + x2_v;

    rvar<T, 1> x1 = x1_v;
    rvar<T, 1> x2 = x2_v;

    rvar<T, 1> orig_x1 = x1;
    rvar<T, 1> orig_x2 = x2;

    x1 += x2;

    BOOST_REQUIRE_EQUAL(x1.item(), expected);

    gradient_tape<T, 1>& tape = get_active_tape<T, 1>();
    tape.zero_grad();
    x1.backward();

    BOOST_REQUIRE_CLOSE(orig_x1.adjoint(), 1.0, 300 * std::numeric_limits<T>::epsilon());
    BOOST_REQUIRE_CLOSE(orig_x2.adjoint(), 1.0, 300 * std::numeric_limits<T>::epsilon());
    tape.clear();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(inplace_subtraction, T, all_float_types)
{
    using namespace rdiff;
    RandomSample<T> rng{-1, 1};
    T               x1_v     = rng.next();
    T               x2_v     = rng.next();
    T               expected = x1_v - x2_v;

    rvar<T, 1> x1 = x1_v;
    rvar<T, 1> x2 = x2_v;

    rvar<T, 1> orig_x1 = x1;
    rvar<T, 1> orig_x2 = x2;

    x1 -= x2;

    BOOST_REQUIRE_EQUAL(x1.item(), expected);

    gradient_tape<T, 1>& tape = get_active_tape<T, 1>();
    tape.zero_grad();
    x1.backward();

    BOOST_REQUIRE_CLOSE(orig_x1.adjoint(), 1.0, 300 * std::numeric_limits<T>::epsilon());
    BOOST_REQUIRE_CLOSE(orig_x2.adjoint(), -1.0, 300 * std::numeric_limits<T>::epsilon());
    tape.clear();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(inplace_multiplication, T, all_float_types)
{
    using namespace rdiff;
    RandomSample<T> rng{-1, 1};
    T               x1_v     = rng.next();
    T               x2_v     = rng.next();
    T               expected = x1_v * x2_v;

    rvar<T, 1> x1 = x1_v;
    rvar<T, 1> x2 = x2_v;

    rvar<T, 1> orig_x1 = x1;
    rvar<T, 1> orig_x2 = x2;

    x1 *= x2;

    BOOST_REQUIRE_EQUAL(x1.item(), expected);

    gradient_tape<T, 1>& tape = get_active_tape<T, 1>();
    tape.zero_grad();
    x1.backward();

    BOOST_REQUIRE_CLOSE(orig_x1.adjoint(), x2_v, 300 * std::numeric_limits<T>::epsilon());
    BOOST_REQUIRE_CLOSE(orig_x2.adjoint(), x1_v, 300 * std::numeric_limits<T>::epsilon());
    tape.clear();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(inplace_division, T, all_float_types)
{
    using namespace rdiff;
    RandomSample<T> rng{-1, 1};
    T               x1_v = rng.next();
    T               x2_v = rng.next() + 1e-2;

    T expected = x1_v / x2_v;

    rvar<T, 1> x1 = x1_v;
    rvar<T, 1> x2 = x2_v;

    rvar<T, 1> orig_x1 = x1;
    rvar<T, 1> orig_x2 = x2;

    x1 /= x2;

    BOOST_REQUIRE_EQUAL(x1.item(), expected);

    gradient_tape<T, 1>& tape = get_active_tape<T, 1>();
    tape.zero_grad();
    x1.backward();

    BOOST_REQUIRE_CLOSE(orig_x1.adjoint(), 1.0 / x2_v, 300 * std::numeric_limits<T>::epsilon());
    BOOST_REQUIRE_CLOSE(orig_x2.adjoint(),
                        -x1_v / (x2_v * x2_v),
                        300 * std::numeric_limits<T>::epsilon());
    tape.clear();
}
BOOST_AUTO_TEST_CASE_TEMPLATE(test_rvar_ostream_output, T, all_float_types)
{
    using namespace rdiff;
    rvar<T, 1> x = 2.0;
    rvar<T, 1> y = 3.0;

    rvar<T, 1> z = x * y;
    z.backward();

    std::ostringstream oss_x, oss_y, oss_z;
    oss_x << x;
    oss_y << y;
    oss_z << z;

    BOOST_CHECK_EQUAL(oss_x.str(), "rvar<1>(2,3)");
    BOOST_CHECK_EQUAL(oss_y.str(), "rvar<1>(3,2)");
    BOOST_CHECK_EQUAL(oss_z.str(), "rvar<1>(6,1)");
}

BOOST_AUTO_TEST_CASE_TEMPLATE(boost_math_tools_numeric_traits, T, all_float_types)
{
    using namespace rdiff;
    using rvar_t = rvar<T, 1>;

    BOOST_CHECK_EQUAL(boost::math::tools::digits<rvar_t>(), boost::math::tools::digits<T>());
    BOOST_CHECK_EQUAL(boost::math::tools::max_value<rvar_t>(), boost::math::tools::max_value<T>());
    BOOST_CHECK_EQUAL(boost::math::tools::min_value<rvar_t>(), boost::math::tools::min_value<T>());
    BOOST_CHECK_EQUAL(boost::math::tools::log_max_value<rvar_t>(),
                      boost::math::tools::log_max_value<T>());
    BOOST_CHECK_EQUAL(boost::math::tools::log_min_value<rvar_t>(),
                      boost::math::tools::log_min_value<T>());
    BOOST_CHECK_EQUAL(boost::math::tools::epsilon<rvar_t>(), boost::math::tools::epsilon<T>());
}
BOOST_AUTO_TEST_SUITE_END()
