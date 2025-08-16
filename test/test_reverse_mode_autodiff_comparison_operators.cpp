//           Copyright Maksym Zhelyenzyakov 2025-2026.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)
#include "test_autodiff_reverse.hpp"
BOOST_AUTO_TEST_SUITE(test_comparison_operators)
using namespace rdiff;

BOOST_AUTO_TEST_CASE_TEMPLATE(equals_rvar_and_rvar, T, all_float_types)
{
    T          a    = T(3.14);
    T          b    = T(3.14);
    T          c    = T(2.71);

    rvar<T, 1> var1 = rvar<T, 1>(a);
    rvar<T, 2> var2 = rvar<T, 2>(b);
    rvar<T, 3> var3 = rvar<T, 3>(a);
    rvar<T, 4> var4 = rvar<T, 4>(c);

    BOOST_CHECK(var1 == var2);
    BOOST_CHECK(!(var1 == var4));
    BOOST_CHECK(var2 == var3);
    BOOST_CHECK(!(var3 == var4));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(equals_rvar_and_scalar, T, all_float_types)
{
    T          a    = T(3.14);
    T          b    = T(3.14);
    T          c    = T(2.71);

    rvar<T, 1> var1 = rvar<T, 1>(a);
    rvar<T, 2> var2 = rvar<T, 2>(a);
    rvar<T, 3> var3 = rvar<T, 3>(a);
    rvar<T, 4> var4 = rvar<T, 4>(a);

    BOOST_CHECK(var1 == b);
    BOOST_CHECK(!(var1 == c));
    BOOST_CHECK(b == var1);
    BOOST_CHECK(!(c == var1));

    BOOST_CHECK(var2 == b);
    BOOST_CHECK(!(var2 == c));
    BOOST_CHECK(b == var2);
    BOOST_CHECK(!(c == var2));

    BOOST_CHECK(var3 == b);
    BOOST_CHECK(!(var3 == c));
    BOOST_CHECK(b == var3);
    BOOST_CHECK(!(c == var3));

    BOOST_CHECK(var4 == b);
    BOOST_CHECK(!(var4 == c));
    BOOST_CHECK(b == var4);
    BOOST_CHECK(!(c == var4));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(inequality_rvar_and_rvar, T, all_float_types)
{
    T          a    = T(3.14);
    T          b    = T(3.14);
    T          c    = T(2.71);

    rvar<T, 1> var1 = rvar<T, 1>(a);
    rvar<T, 2> var2 = rvar<T, 2>(b);
    rvar<T, 3> var3 = rvar<T, 3>(a);
    rvar<T, 4> var4 = rvar<T, 4>(c);

    BOOST_CHECK(var1 != var4);
    BOOST_CHECK(!(var1 != var2));
    BOOST_CHECK(var2 != var4);
    BOOST_CHECK(!(var3 != var2));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(inequality_rvar_and_scalar, T, all_float_types)
{
    T          a    = T(3.14);
    T          b    = T(3.14);
    T          c    = T(2.71);

    rvar<T, 1> var1 = rvar<T, 1>(a);
    rvar<T, 2> var2 = rvar<T, 2>(a);
    rvar<T, 3> var3 = rvar<T, 3>(a);
    rvar<T, 4> var4 = rvar<T, 4>(a);

    BOOST_CHECK(var1 != c);
    BOOST_CHECK(!(var1 != b));
    BOOST_CHECK(c != var1);
    BOOST_CHECK(!(b != var1));

    BOOST_CHECK(var2 != c);
    BOOST_CHECK(!(var2 != b));
    BOOST_CHECK(c != var2);
    BOOST_CHECK(!(b != var2));

    BOOST_CHECK(var3 != c);
    BOOST_CHECK(!(var3 != b));
    BOOST_CHECK(c != var3);
    BOOST_CHECK(!(b != var3));

    BOOST_CHECK(var4 != c);
    BOOST_CHECK(!(var4 != b));
    BOOST_CHECK(c != var4);
    BOOST_CHECK(!(b != var4));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(less_than_rvar_and_rvar, T, all_float_types)
{
    T          a    = T(3.14);
    T          b    = T(3.14);
    T          c    = T(2.71);

    rvar<T, 1> var1 = rvar<T, 1>(a);
    rvar<T, 2> var2 = rvar<T, 2>(b);
    rvar<T, 3> var3 = rvar<T, 3>(c);
    rvar<T, 4> var4 = rvar<T, 4>(c);

    BOOST_CHECK(var3 < var1);
    BOOST_CHECK(!(var1 < var2));
    BOOST_CHECK(!(var1 < var3));
    BOOST_CHECK(!(var3 < var4));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(greater_than_rvar_and_rvar, T, all_float_types)
{
    T          a    = T(3.14);
    T          b    = T(3.14);
    T          c    = T(2.71);

    rvar<T, 1> var1 = rvar<T, 1>(a);
    rvar<T, 2> var2 = rvar<T, 2>(b);
    rvar<T, 3> var3 = rvar<T, 3>(c);
    rvar<T, 4> var4 = rvar<T, 4>(c);

    BOOST_CHECK(var1 > var3);
    BOOST_CHECK(!(var1 > var2));
    BOOST_CHECK(!(var3 > var1));
    BOOST_CHECK(!(var3 > var4));
}
BOOST_AUTO_TEST_CASE_TEMPLATE(less_than_rvar_and_scalar, T, all_float_types)
{
    T          a    = T(3.14);
    T          b    = T(3.14);
    T          c    = T(2.71);
    T          d    = T(4.0);

    rvar<T, 1> var1 = rvar<T, 1>(a);
    rvar<T, 2> var2 = rvar<T, 2>(a);
    rvar<T, 3> var3 = rvar<T, 3>(c);
    rvar<T, 4> var4 = rvar<T, 4>(c);

    BOOST_CHECK(var3 < a);
    BOOST_CHECK(!(var1 < b));
    BOOST_CHECK(c < var1);
    BOOST_CHECK(!(d < var2));

    BOOST_CHECK(var4 < a);
    BOOST_CHECK(!(var2 < b));
    BOOST_CHECK(c < var2);
    BOOST_CHECK(!(d < var1));

    BOOST_CHECK(var3 < b);
    BOOST_CHECK(var4 < b);
    BOOST_CHECK(c < var1);
    BOOST_CHECK(!(d < var3));

    BOOST_CHECK(var3 < a);
    BOOST_CHECK(!(var2 < a));
    BOOST_CHECK(!(c < var4));
    BOOST_CHECK(!(d < var1));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(greater_than_rvar_and_scalar, T, all_float_types)
{
    T          a    = T(3.14);
    T          b    = T(3.14);
    T          c    = T(2.71);

    rvar<T, 1> var1 = rvar<T, 1>(a);
    rvar<T, 2> var2 = rvar<T, 2>(a);
    rvar<T, 3> var3 = rvar<T, 3>(c);
    rvar<T, 4> var4 = rvar<T, 4>(c);

    BOOST_CHECK(var1 > c);
    BOOST_CHECK(!(var1 > b));
    BOOST_CHECK(b > var3);
    BOOST_CHECK(!(c > var4));

    BOOST_CHECK(var2 > c);
    BOOST_CHECK(!(var2 > b));
    BOOST_CHECK(b > var4);
    BOOST_CHECK(!(c > var3));

    BOOST_CHECK(var1 > c);
    BOOST_CHECK(!(var1 > b));
    BOOST_CHECK(a > var4);
    BOOST_CHECK(!(c > var3));

    BOOST_CHECK(var2 > c);
    BOOST_CHECK(!(var2 > b));
    BOOST_CHECK(a > var4);
    BOOST_CHECK(!(c > var3));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(less_than_or_equal_rvar_and_rvar, T, all_float_types)
{
    T          a    = T(3.14);
    T          b    = T(3.14);
    T          c    = T(2.71);

    rvar<T, 1> var1 = rvar<T, 1>(a);
    rvar<T, 2> var2 = rvar<T, 2>(b);
    rvar<T, 3> var3 = rvar<T, 3>(c);
    rvar<T, 4> var4 = rvar<T, 4>(c);

    BOOST_CHECK(var3 <= var1);
    BOOST_CHECK(var1 <= var2);
    BOOST_CHECK(!(var1 <= var3));
    BOOST_CHECK(var3 <= var4);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(less_than_or_equal_rvar_and_scalar, T, all_float_types)
{
    T          a    = T(3.14);
    T          b    = T(3.14);
    T          c    = T(2.71);
    T          d    = T(4.0);

    rvar<T, 1> var1 = rvar<T, 1>(a);
    rvar<T, 2> var2 = rvar<T, 2>(a);
    rvar<T, 3> var3 = rvar<T, 3>(c);
    rvar<T, 4> var4 = rvar<T, 4>(c);

    BOOST_CHECK(var3 <= a);
    BOOST_CHECK(var1 <= b);
    BOOST_CHECK(c <= var1);
    BOOST_CHECK(!(d <= var2));

    BOOST_CHECK(var4 <= a);
    BOOST_CHECK(var2 <= b);
    BOOST_CHECK(c <= var2);
    BOOST_CHECK(!(d <= var1));

    BOOST_CHECK(var3 <= b);
    BOOST_CHECK(var4 <= b);
    BOOST_CHECK(c <= var1);
    BOOST_CHECK(!(d <= var3));

    BOOST_CHECK(var3 <= a);
    BOOST_CHECK(var2 <= a);
    BOOST_CHECK(c <= var4);
    BOOST_CHECK(!(d <= var1));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(greater_than_or_equal_rvar_and_rvar, T, all_float_types)
{
    T          a    = T(3.14);
    T          b    = T(3.14);
    T          c    = T(2.71);

    rvar<T, 1> var1 = rvar<T, 1>(a);
    rvar<T, 2> var2 = rvar<T, 2>(b);
    rvar<T, 3> var3 = rvar<T, 3>(c);
    rvar<T, 4> var4 = rvar<T, 4>(c);

    BOOST_CHECK(var1 >= var3);
    BOOST_CHECK(var1 >= var2);
    BOOST_CHECK(!(var3 >= var1));
    BOOST_CHECK(var3 >= var4);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(greater_than_or_equal_rvar_and_scalar, T, all_float_types)
{
    T          a    = T(3.14);
    T          b    = T(3.14);
    T          c    = T(2.71);

    rvar<T, 1> var1 = rvar<T, 1>(a);
    rvar<T, 2> var2 = rvar<T, 2>(a);
    rvar<T, 3> var3 = rvar<T, 3>(c);
    rvar<T, 4> var4 = rvar<T, 4>(c);

    BOOST_CHECK(var1 >= c);
    BOOST_CHECK(var1 >= b);
    BOOST_CHECK(b >= var3);
    BOOST_CHECK(c >= var4);

    BOOST_CHECK(var2 >= c);
    BOOST_CHECK(var2 >= b);
    BOOST_CHECK(b >= var4);
    BOOST_CHECK(c >= var3);

    BOOST_CHECK(var1 >= c);
    BOOST_CHECK(var1 >= b);
    BOOST_CHECK(a >= var4);
    BOOST_CHECK(c >= var3);

    BOOST_CHECK(var2 >= c);
    BOOST_CHECK(var2 >= b);
    BOOST_CHECK(a >= var4);
    BOOST_CHECK(c >= var3);
}
BOOST_AUTO_TEST_SUITE_END()
