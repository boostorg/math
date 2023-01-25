//  (C) Copyright Matt Borland 2023.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "math_unit_test.hpp"
#include <boost/math/special_functions/powm1.hpp>
#include <limits>

template <typename T>
void test()
{
    using namespace boost::math::policies;
    
    constexpr T x = -std::numeric_limits<T>::infinity();
    constexpr T y = T(3);

    T p = boost::math::powm1(x, y, make_policy(overflow_error<errno_on_error>()));

    CHECK_EQUAL(p, -std::numeric_limits<T>::infinity());
}

int main(void)
{
    test<float>();
    test<double>();
    test<long double>();

    return boost::math::test::report_errors();
}
