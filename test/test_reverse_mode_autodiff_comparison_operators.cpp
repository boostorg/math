//           Copyright Maksym Zhelyeznyakov 2025-2026.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)
#include "test_autodiff_reverse.hpp"
using namespace rdiff;

template <class T>
void equals_rvar_and_rvar()
{
    T          a    = T(3.14);
    T          b    = T(3.14);
    T          c    = T(2.71);

    rvar<T, 1> var1 = rvar<T, 1>(a);
    rvar<T, 2> var2 = rvar<T, 2>(b);
    rvar<T, 3> var3 = rvar<T, 3>(a);
    rvar<T, 4> var4 = rvar<T, 4>(c);

    CHECK_TRUE(var1 == var2);
    CHECK_TRUE(!(var1 == var4));
    CHECK_TRUE(var2 == var3);
    CHECK_TRUE(!(var3 == var4));
}

template <class T>
void equals_rvar_and_scalar()
{
    T          a    = T(3.14);
    T          b    = T(3.14);
    T          c    = T(2.71);

    rvar<T, 1> var1 = rvar<T, 1>(a);
    rvar<T, 2> var2 = rvar<T, 2>(a);
    rvar<T, 3> var3 = rvar<T, 3>(a);
    rvar<T, 4> var4 = rvar<T, 4>(a);

    CHECK_TRUE(var1 == b);
    CHECK_TRUE(!(var1 == c));
    CHECK_TRUE(b == var1);
    CHECK_TRUE(!(c == var1));

    CHECK_TRUE(var2 == b);
    CHECK_TRUE(!(var2 == c));
    CHECK_TRUE(b == var2);
    CHECK_TRUE(!(c == var2));

    CHECK_TRUE(var3 == b);
    CHECK_TRUE(!(var3 == c));
    CHECK_TRUE(b == var3);
    CHECK_TRUE(!(c == var3));

    CHECK_TRUE(var4 == b);
    CHECK_TRUE(!(var4 == c));
    CHECK_TRUE(b == var4);
    CHECK_TRUE(!(c == var4));
}

template <class T>
void inequality_rvar_and_rvar()
{
    T          a    = T(3.14);
    T          b    = T(3.14);
    T          c    = T(2.71);

    rvar<T, 1> var1 = rvar<T, 1>(a);
    rvar<T, 2> var2 = rvar<T, 2>(b);
    rvar<T, 3> var3 = rvar<T, 3>(a);
    rvar<T, 4> var4 = rvar<T, 4>(c);

    CHECK_TRUE(var1 != var4);
    CHECK_TRUE(!(var1 != var2));
    CHECK_TRUE(var2 != var4);
    CHECK_TRUE(!(var3 != var2));
}

template <class T>
void inequality_rvar_and_scalar()
{
    T          a    = T(3.14);
    T          b    = T(3.14);
    T          c    = T(2.71);

    rvar<T, 1> var1 = rvar<T, 1>(a);
    rvar<T, 2> var2 = rvar<T, 2>(a);
    rvar<T, 3> var3 = rvar<T, 3>(a);
    rvar<T, 4> var4 = rvar<T, 4>(a);

    CHECK_TRUE(var1 != c);
    CHECK_TRUE(!(var1 != b));
    CHECK_TRUE(c != var1);
    CHECK_TRUE(!(b != var1));

    CHECK_TRUE(var2 != c);
    CHECK_TRUE(!(var2 != b));
    CHECK_TRUE(c != var2);
    CHECK_TRUE(!(b != var2));

    CHECK_TRUE(var3 != c);
    CHECK_TRUE(!(var3 != b));
    CHECK_TRUE(c != var3);
    CHECK_TRUE(!(b != var3));

    CHECK_TRUE(var4 != c);
    CHECK_TRUE(!(var4 != b));
    CHECK_TRUE(c != var4);
    CHECK_TRUE(!(b != var4));
}

template <class T>
void less_than_rvar_and_rvar()
{
    T          a    = T(3.14);
    T          b    = T(3.14);
    T          c    = T(2.71);

    rvar<T, 1> var1 = rvar<T, 1>(a);
    rvar<T, 2> var2 = rvar<T, 2>(b);
    rvar<T, 3> var3 = rvar<T, 3>(c);
    rvar<T, 4> var4 = rvar<T, 4>(c);

    CHECK_TRUE(var3 < var1);
    CHECK_TRUE(!(var1 < var2));
    CHECK_TRUE(!(var1 < var3));
    CHECK_TRUE(!(var3 < var4));
}

template <class T>
void greater_than_rvar_and_rvar()
{
    T          a    = T(3.14);
    T          b    = T(3.14);
    T          c    = T(2.71);

    rvar<T, 1> var1 = rvar<T, 1>(a);
    rvar<T, 2> var2 = rvar<T, 2>(b);
    rvar<T, 3> var3 = rvar<T, 3>(c);
    rvar<T, 4> var4 = rvar<T, 4>(c);

    CHECK_TRUE(var1 > var3);
    CHECK_TRUE(!(var1 > var2));
    CHECK_TRUE(!(var3 > var1));
    CHECK_TRUE(!(var3 > var4));
}
template <class T>
void less_than_rvar_and_scalar()
{
    T          a    = T(3.14);
    T          b    = T(3.14);
    T          c    = T(2.71);
    T          d    = T(4.0);

    rvar<T, 1> var1 = rvar<T, 1>(a);
    rvar<T, 2> var2 = rvar<T, 2>(a);
    rvar<T, 3> var3 = rvar<T, 3>(c);
    rvar<T, 4> var4 = rvar<T, 4>(c);

    CHECK_TRUE(var3 < a);
    CHECK_TRUE(!(var1 < b));
    CHECK_TRUE(c < var1);
    CHECK_TRUE(!(d < var2));

    CHECK_TRUE(var4 < a);
    CHECK_TRUE(!(var2 < b));
    CHECK_TRUE(c < var2);
    CHECK_TRUE(!(d < var1));

    CHECK_TRUE(var3 < b);
    CHECK_TRUE(var4 < b);
    CHECK_TRUE(c < var1);
    CHECK_TRUE(!(d < var3));

    CHECK_TRUE(var3 < a);
    CHECK_TRUE(!(var2 < a));
    CHECK_TRUE(!(c < var4));
    CHECK_TRUE(!(d < var1));
}

template <class T>
void greater_than_rvar_and_scalar()
{
    T          a    = T(3.14);
    T          b    = T(3.14);
    T          c    = T(2.71);

    rvar<T, 1> var1 = rvar<T, 1>(a);
    rvar<T, 2> var2 = rvar<T, 2>(a);
    rvar<T, 3> var3 = rvar<T, 3>(c);
    rvar<T, 4> var4 = rvar<T, 4>(c);

    CHECK_TRUE(var1 > c);
    CHECK_TRUE(!(var1 > b));
    CHECK_TRUE(b > var3);
    CHECK_TRUE(!(c > var4));

    CHECK_TRUE(var2 > c);
    CHECK_TRUE(!(var2 > b));
    CHECK_TRUE(b > var4);
    CHECK_TRUE(!(c > var3));

    CHECK_TRUE(var1 > c);
    CHECK_TRUE(!(var1 > b));
    CHECK_TRUE(a > var4);
    CHECK_TRUE(!(c > var3));

    CHECK_TRUE(var2 > c);
    CHECK_TRUE(!(var2 > b));
    CHECK_TRUE(a > var4);
    CHECK_TRUE(!(c > var3));
}

template <class T>
void less_than_or_equal_rvar_and_rvar()
{
    T          a    = T(3.14);
    T          b    = T(3.14);
    T          c    = T(2.71);

    rvar<T, 1> var1 = rvar<T, 1>(a);
    rvar<T, 2> var2 = rvar<T, 2>(b);
    rvar<T, 3> var3 = rvar<T, 3>(c);
    rvar<T, 4> var4 = rvar<T, 4>(c);

    CHECK_TRUE(var3 <= var1);
    CHECK_TRUE(var1 <= var2);
    CHECK_TRUE(!(var1 <= var3));
    CHECK_TRUE(var3 <= var4);
}

template <class T>
void less_than_or_equal_rvar_and_scalar()
{
    T          a    = T(3.14);
    T          b    = T(3.14);
    T          c    = T(2.71);
    T          d    = T(4.0);

    rvar<T, 1> var1 = rvar<T, 1>(a);
    rvar<T, 2> var2 = rvar<T, 2>(a);
    rvar<T, 3> var3 = rvar<T, 3>(c);
    rvar<T, 4> var4 = rvar<T, 4>(c);

    CHECK_TRUE(var3 <= a);
    CHECK_TRUE(var1 <= b);
    CHECK_TRUE(c <= var1);
    CHECK_TRUE(!(d <= var2));

    CHECK_TRUE(var4 <= a);
    CHECK_TRUE(var2 <= b);
    CHECK_TRUE(c <= var2);
    CHECK_TRUE(!(d <= var1));

    CHECK_TRUE(var3 <= b);
    CHECK_TRUE(var4 <= b);
    CHECK_TRUE(c <= var1);
    CHECK_TRUE(!(d <= var3));

    CHECK_TRUE(var3 <= a);
    CHECK_TRUE(var2 <= a);
    CHECK_TRUE(c <= var4);
    CHECK_TRUE(!(d <= var1));
}

template <class T>
void greater_than_or_equal_rvar_and_rvar()
{
    T          a    = T(3.14);
    T          b    = T(3.14);
    T          c    = T(2.71);

    rvar<T, 1> var1 = rvar<T, 1>(a);
    rvar<T, 2> var2 = rvar<T, 2>(b);
    rvar<T, 3> var3 = rvar<T, 3>(c);
    rvar<T, 4> var4 = rvar<T, 4>(c);

    CHECK_TRUE(var1 >= var3);
    CHECK_TRUE(var1 >= var2);
    CHECK_TRUE(!(var3 >= var1));
    CHECK_TRUE(var3 >= var4);
}

template <class T>
void greater_than_or_equal_rvar_and_scalar()
{
    T          a    = T(3.14);
    T          b    = T(3.14);
    T          c    = T(2.71);

    rvar<T, 1> var1 = rvar<T, 1>(a);
    rvar<T, 2> var2 = rvar<T, 2>(a);
    rvar<T, 3> var3 = rvar<T, 3>(c);
    rvar<T, 4> var4 = rvar<T, 4>(c);

    CHECK_TRUE(var1 >= c);
    CHECK_TRUE(var1 >= b);
    CHECK_TRUE(b >= var3);
    CHECK_TRUE(c >= var4);

    CHECK_TRUE(var2 >= c);
    CHECK_TRUE(var2 >= b);
    CHECK_TRUE(b >= var4);
    CHECK_TRUE(c >= var3);

    CHECK_TRUE(var1 >= c);
    CHECK_TRUE(var1 >= b);
    CHECK_TRUE(a >= var4);
    CHECK_TRUE(c >= var3);

    CHECK_TRUE(var2 >= c);
    CHECK_TRUE(var2 >= b);
    CHECK_TRUE(a >= var4);
    CHECK_TRUE(c >= var3);
}

int main()
{
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { equals_rvar_and_rvar<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { equals_rvar_and_scalar<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { inequality_rvar_and_rvar<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { inequality_rvar_and_scalar<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { less_than_rvar_and_rvar<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { greater_than_rvar_and_rvar<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { less_than_rvar_and_scalar<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { greater_than_rvar_and_scalar<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { less_than_or_equal_rvar_and_rvar<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { less_than_or_equal_rvar_and_scalar<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { greater_than_or_equal_rvar_and_rvar<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { greater_than_or_equal_rvar_and_scalar<typename decltype(t)::type>(); });
   return boost::math::test::report_errors();
}
