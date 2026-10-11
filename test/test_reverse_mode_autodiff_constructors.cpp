//           Copyright Maksym Zhelyeznyakov 2025-2026.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)

#include "test_autodiff_reverse.hpp"


template <class T>
void rvar_constructors_and_utils()
{
    RandomSample<T>                          rng{T(-100), T(100)};
    /* raw constructors */
    T                                        x_value = rng.next();
    rdiff::rvar<T, 1>                        x0(x_value);
    rdiff::rvar<T, 2>                        x1(x_value);
    rdiff::rvar<T, 3>                        x2(x_value);

    /* raw eval */
    CHECK_EQUAL(x_value, x0.evaluate());
    CHECK_EQUAL(x_value, x1.evaluate().evaluate());
    CHECK_EQUAL(x_value, x2.evaluate().evaluate().evaluate());

    /* get item helper */
    CHECK_EQUAL(x_value, x0.item());
    CHECK_EQUAL(x_value, x1.item());
    CHECK_EQUAL(x_value, x2.item());
}

template <class T>
void make_rvar_constructors()
{
    RandomSample<T> rng{T(-100), T(100)};
    T               x_value = rng.next();

    rdiff::rvar<T,1> x1 = rdiff::make_rvar<T,1>(x_value);
    rdiff::rvar<T,2> x2 = rdiff::make_rvar<T,2>(x_value);
    rdiff::rvar<T,3> x3 = rdiff::make_rvar<T,3>(x_value);
    CHECK_EQUAL(x_value, x1.item());
    CHECK_EQUAL(x_value, x2.item());
    CHECK_EQUAL(x_value, x3.item());

    T result = x_value+x_value;
    auto z1 = x1+x1;
    auto z2 = x2+x2;
    auto z3 = x3+x3;

    rdiff::rvar<T,1> y1 = rdiff::make_rvar(z1);
    rdiff::rvar<T,2> y2 = rdiff::make_rvar(z2);
    rdiff::rvar<T,3> y3 = rdiff::make_rvar(z3);

    CHECK_EQUAL(y1.item(),result);
    CHECK_EQUAL(y2.item(),result);
    CHECK_EQUAL(y3.item(),result);
}

template <class T>
void copy_construction()
{
    RandomSample<T> rng{-100, 100};
    T               original_value_1 = rng.next();

    /* copy constructor checks */
    rdiff::rvar<T, 1> original_rvar_1(original_value_1);
    rdiff::rvar<T, 1> copied_rvar_1 = original_rvar_1;

    CHECK_EQUAL(copied_rvar_1.item(), original_value_1);

    /* change original, copy stays same */
    T new_value_1   = rng.next();
    original_rvar_1 = new_value_1;
    CHECK_EQUAL(copied_rvar_1.item(), original_value_1);
    CHECK_EQUAL(original_rvar_1.item(), new_value_1);

    T                 original_value_2 = rng.next();
    rdiff::rvar<T, 2> original_rvar_2(original_value_2);
    rdiff::rvar<T, 2> copied_rvar_2 = original_rvar_2;

    CHECK_EQUAL(copied_rvar_2.item(), original_value_2);

    T new_value_2   = rng.next();
    original_rvar_2 = new_value_2;
    CHECK_EQUAL(copied_rvar_2.item(), original_value_2);
    CHECK_EQUAL(original_rvar_2.item(), new_value_2);

    T                 original_value_3 = rng.next();
    rdiff::rvar<T, 3> original_rvar_3(original_value_3);
    rdiff::rvar<T, 3> copied_rvar_3 = original_rvar_3;

    CHECK_EQUAL(copied_rvar_3.item(), original_value_3);

    T new_value_3   = rng.next();
    original_rvar_3 = new_value_3;
    CHECK_EQUAL(copied_rvar_3.item(), original_value_3);
    CHECK_EQUAL(original_rvar_3.item(), new_value_3);
}
template <class T>
void copy_assignment()
{
    RandomSample<T>   rng{-100, 100};
    T                 source_value_1         = rng.next();
    T                 initial_target_value_1 = rng.next();
    rdiff::rvar<T, 1> source_rvar_1(source_value_1);
    rdiff::rvar<T, 1> target_rvar_1(initial_target_value_1);

    target_rvar_1 = source_rvar_1;

    CHECK_EQUAL(target_rvar_1.item(), source_value_1);

    T new_source_value_1 = rng.next();
    source_rvar_1        = new_source_value_1;
    CHECK_EQUAL(target_rvar_1.item(), source_value_1);
    CHECK_EQUAL(source_rvar_1.item(), new_source_value_1);

    rdiff::rvar<T, 1> self_assign_rvar_1(rng.next());
    T                 value_before_self_assign_1 = self_assign_rvar_1.item();
    self_assign_rvar_1                           = self_assign_rvar_1;
    CHECK_EQUAL(self_assign_rvar_1.item(), value_before_self_assign_1);

    T                 source_value_2         = rng.next();
    T                 initial_target_value_2 = rng.next();
    rdiff::rvar<T, 2> source_rvar_2(source_value_2);
    rdiff::rvar<T, 2> target_rvar_2(initial_target_value_2);

    target_rvar_2 = source_rvar_2;

    CHECK_EQUAL(target_rvar_2.item(), source_value_2);

    T new_source_value_2 = rng.next();
    source_rvar_2        = new_source_value_2;
    CHECK_EQUAL(target_rvar_2.item(), source_value_2);
    CHECK_EQUAL(source_rvar_2.item(), new_source_value_2);

    T                 source_value_3         = rng.next();
    T                 initial_target_value_3 = rng.next();
    rdiff::rvar<T, 3> source_rvar_3(source_value_3);
    rdiff::rvar<T, 3> target_rvar_3(initial_target_value_3);

    target_rvar_3 = source_rvar_3;

    CHECK_EQUAL(target_rvar_3.item(), source_value_3);

    T new_source_value_3 = rng.next();
    source_rvar_3        = new_source_value_3;
    CHECK_EQUAL(target_rvar_3.item(), source_value_3);
    CHECK_EQUAL(source_rvar_3.item(), new_source_value_3);
}
template <class T>
void inplace_addition()
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

    if (!CHECK_EQUAL(x1.item(), expected)) { return; }

    gradient_tape<T, 1>& tape = get_active_tape<T, 1>();
    tape.zero_grad();
    x1.backward();

    if (!CHECK_RELATIVE_ERROR(orig_x1.adjoint(), 1.0, (300 * std::numeric_limits<T>::epsilon()) / 100)) { return; }
    if (!CHECK_RELATIVE_ERROR(orig_x2.adjoint(), 1.0, (300 * std::numeric_limits<T>::epsilon()) / 100)) { return; }
    tape.clear();
}

template <class T>
void inplace_subtraction()
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

    if (!CHECK_EQUAL(x1.item(), expected)) { return; }

    gradient_tape<T, 1>& tape = get_active_tape<T, 1>();
    tape.zero_grad();
    x1.backward();

    if (!CHECK_RELATIVE_ERROR(orig_x1.adjoint(), 1.0, (300 * std::numeric_limits<T>::epsilon()) / 100)) { return; }
    if (!CHECK_RELATIVE_ERROR(orig_x2.adjoint(), -1.0, (300 * std::numeric_limits<T>::epsilon()) / 100)) { return; }
    tape.clear();
}

template <class T>
void inplace_multiplication()
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

    if (!CHECK_EQUAL(x1.item(), expected)) { return; }

    gradient_tape<T, 1>& tape = get_active_tape<T, 1>();
    tape.zero_grad();
    x1.backward();

    if (!CHECK_RELATIVE_ERROR(orig_x1.adjoint(), x2_v, (300 * std::numeric_limits<T>::epsilon()) / 100)) { return; }
    if (!CHECK_RELATIVE_ERROR(orig_x2.adjoint(), x1_v, (300 * std::numeric_limits<T>::epsilon()) / 100)) { return; }
    tape.clear();
}

template <class T>
void inplace_division()
{
    using namespace rdiff;
    RandomSample<T> rng{T(-1), T(1)};
    T               x1_v = rng.next();
    T               x2_v = rng.next() + T(1e-2);

    T expected = x1_v / x2_v;

    rvar<T, 1> x1 = x1_v;
    rvar<T, 1> x2 = x2_v;

    rvar<T, 1> orig_x1 = x1;
    rvar<T, 1> orig_x2 = x2;

    x1 /= x2;

    if (!CHECK_EQUAL(x1.item(), expected)) { return; }

    gradient_tape<T, 1>& tape = get_active_tape<T, 1>();
    tape.zero_grad();
    x1.backward();

    if (!CHECK_RELATIVE_ERROR(orig_x1.adjoint(), 1.0 / x2_v, (300 * std::numeric_limits<T>::epsilon()) / 100)) { return; }
    if (!CHECK_RELATIVE_ERROR(orig_x2.adjoint(),
                        -x1_v / (x2_v * x2_v), (300 * std::numeric_limits<T>::epsilon()) / 100)) { return; }
    tape.clear();
}
template <class T>
void test_rvar_ostream_output()
{
    using namespace rdiff;
    rvar<T, 1> x = static_cast<T>(2.0);
    rvar<T, 1> y = static_cast<T>(3.0);

    rvar<T, 1> z = x * y;
    z.backward();

    std::ostringstream oss_x, oss_y, oss_z;
    oss_x << x;
    oss_y << y;
    oss_z << z;

    CHECK_EQUAL(oss_x.str(), std::string("rvar<1>(2,3)"));
    CHECK_EQUAL(oss_y.str(), std::string("rvar<1>(3,2)"));
    CHECK_EQUAL(oss_z.str(), std::string("rvar<1>(6,1)"));
}

template <class T>
void boost_math_tools_numeric_traits()
{
    using namespace rdiff;
    using rvar_t = rvar<T, 1>;

    CHECK_EQUAL(boost::math::tools::digits<rvar_t>(), boost::math::tools::digits<T>());
    CHECK_EQUAL(static_cast<T>(boost::math::tools::max_value<rvar_t>()), boost::math::tools::max_value<T>());
    CHECK_EQUAL(static_cast<T>(boost::math::tools::min_value<rvar_t>()), boost::math::tools::min_value<T>());
    CHECK_EQUAL(static_cast<T>(boost::math::tools::log_max_value<rvar_t>()),
                      boost::math::tools::log_max_value<T>());
    CHECK_EQUAL(static_cast<T>(boost::math::tools::log_min_value<rvar_t>()),
                      boost::math::tools::log_min_value<T>());
    CHECK_EQUAL(static_cast<T>(boost::math::tools::epsilon<rvar_t>()), boost::math::tools::epsilon<T>());
}

int main()
{
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { rvar_constructors_and_utils<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { make_rvar_constructors<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { copy_construction<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { copy_assignment<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { inplace_addition<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { inplace_subtraction<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { inplace_multiplication<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { inplace_division<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { test_rvar_ostream_output<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { boost_math_tools_numeric_traits<typename decltype(t)::type>(); });
   return boost::math::test::report_errors();
}
