//  (C) Copyright Matt Borland 2022.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <limits>
#include <boost/math/constants/constants.hpp>
#include <boost/math/ccmath/sin.hpp>
#include <boost/math/ccmath/isinf.hpp>
#include <boost/math/ccmath/isnan.hpp>
#include <boost/math/ccmath/abs.hpp>

#if !defined(BOOST_MATH_NO_CONSTEXPR_DETECTION) && !defined(BOOST_MATH_USING_BUILTIN_CONSTANT_P)
template <typename T>
void test()
{
    constexpr T half_pi = boost::math::constants::half_pi<T>();
    constexpr T quarter_pi = boost::math::constants::quarter_pi<T>();
    constexpr T one_sqrt2 = boost::math::constants::one_div_root_two<T>();
    constexpr T tol = std::is_same_v<T, long double> ? 2e-16L : std::numeric_limits<T>::epsilon();
    
    // Edge cases
    static_assert(boost::math::ccmath::sin(T(0)) == 0);
    static_assert(boost::math::ccmath::isnan(boost::math::ccmath::sin(std::numeric_limits<T>::quiet_NaN())));
    static_assert(boost::math::ccmath::isinf(boost::math::ccmath::sin(std::numeric_limits<T>::infinity())));

    // sin(pi/2) == 1
    static_assert(boost::math::ccmath::abs(boost::math::ccmath::sin(half_pi) - T(1)) < tol);
    
    // sin(pi/4) = 1/sqrt(2)
    static_assert(boost::math::ccmath::abs(boost::math::ccmath::sin(quarter_pi) - one_sqrt2) < tol);

    // sin(-pi/2) == -1
    static_assert(boost::math::ccmath::abs(boost::math::ccmath::sin(-half_pi) + T(1)) < tol);

    // sin(-pi/4) = -1/sqrt(2)
    static_assert(boost::math::ccmath::abs(boost::math::ccmath::sin(-quarter_pi) + one_sqrt2) < tol);
    
}

int main(void)
{
    test<float>();
    test<double>();
    test<long double>();
}
#else
int main(void)
{
    return 0;
}
#endif
