//  (C) Copyright Matt Borland 2023.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/concepts/real_concept.hpp>

int main(void)
{
    using namespace boost::math::concepts;
    
    constexpr real_concept val1(1.0);
    constexpr real_concept val2(2.0);

    // Member functions
    static_assert(val1 < val2, "val1 should be less than val2");
    static_assert(val1 + val2 >= 2.999L, "Addition is 3");
    static_assert(val1 * val2 >= 1.999L, "Multiplication is 2");
    static_assert(val2 / val1 >= 1.999L, "Division is 2");

    // Tools
    static_assert(boost::math::tools::max_value<real_concept>() > 0, "max value");
    static_assert(boost::math::tools::min_value<real_concept>() > 0, "min value");
    static_assert(boost::math::tools::epsilon<real_concept>() > 0, "Epsilon");
    static_assert(boost::math::tools::digits<real_concept>() > 0, "Digits");
}
