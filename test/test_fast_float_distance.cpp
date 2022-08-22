//  (C) Copyright John Maddock 2008.
//  (C) Copyright Matt Borland 2022.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/concepts/real_concept.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/math/special_functions/next.hpp>
#include <boost/math/special_functions/ulp.hpp>
#include <boost/math/special_functions/fast_float_distance.hpp>
#include <boost/multiprecision/float128.hpp>

#include "math_unit_test.hpp"


template <class T>
void test_value(const T& val)
{
    using namespace boost::math;

    assert(fast_float_distance(float_next(val), val) == -1);
    assert(float_next(val) > val);
    assert(float_next(float_prior(val)) == val);

    assert(fast_float_distance(float_advance(val, 4), val) == -4);
    assert(fast_float_distance(float_advance(val, -4), val) == 4);
    if(std::numeric_limits<T>::is_specialized && (std::numeric_limits<T>::has_denorm == std::denorm_present))
    {
        assert(fast_float_distance(float_advance(float_next(float_next(val)), 4), float_next(float_next(val))) == -4);
        assert(fast_float_distance(float_advance(float_next(float_next(val)), -4), float_next(float_next(val))) == 4);
    }
}

int main(void)
{
    test_value(1.0f);
    test_value(1.0);

    #ifdef BOOST_MATH_USE_FAST_FLOAT128
    test_value(boost::multiprecision::float128_type(0));
    test_value(__float128(0));
    #endif
}
