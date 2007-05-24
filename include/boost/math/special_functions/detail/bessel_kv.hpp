//  Copyright (c) 2006 Xiaogang Zhang
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_BESSEL_KV_HPP
#define BOOST_MATH_BESSEL_KV_HPP

#include <boost/math/special_functions/bessel_kn.hpp>
#include <boost/math/special_functions/bessel_ik.hpp>
#include <boost/math/constants/constants.hpp>

// Modified Bessel and Spherical Bessel functions of the second kind

namespace boost { namespace math {

// Modified Bessel function of the second kind of any order
template <typename T>
inline T bessel_kv(T v, T x)
{
    int n = static_cast<int>(v);
    if (v == n)
    {
        return bessel_kn(n, x);             // v is integer
    }
    else
    {
        T I, K;
        bessel_ik(v, x, &I, &K);
        return K;
    }
}

// Modified Spherical Bessel function of the second kind of non-negative order
template <typename T>
inline T bessel_skv(unsigned n, T x)
{
    using namespace std;
    using namespace boost::math::tools;
    using namespace boost::math::constants;

    T v = n + 0.5L;
    if (x == 0)
    {
        return overflow_error<T>("boost::math::bessel_skv(n, x)",
            "infinity occurred but not supported");
    }
    else
    {
        return sqrt(0.5L * pi<T>() / x) * bessel_kv(v, x);
    }
}

// -------------------- TR1 functions --------------------

inline float cyl_bessel_kf(float nu, float x)
{
    if (nu >= 128)
    {
        std::cout << "Warning: cyl_bessel_kf(nu, x), nu >= 128, "
                  << "result is implementation defined according to C++ TR1"
                  << std::endl;
    }
    return bessel_kv<float>(nu, x);
}

inline double cyl_bessel_k(double nu, double x)
{
    if (nu >= 128)
    {
        std::cout << "Warning: cyl_bessel_k(nu, x), nu >= 128, "
                  << "result is implementation defined according to C++ TR1"
                  << std::endl;
    }
    return bessel_kv<double>(nu, x);
}

inline long double cyl_bessel_kl(long double nu, long double x)
{
    if (nu >= 128)
    {
        std::cout << "Warning: cyl_bessel_kl(nu, x), nu >= 128, "
                  << "result is implementation defined according to C++ TR1"
                  << std::endl;
    }
    return bessel_kv<long double>(nu, x);
}

}} // namespaces

#endif // BOOST_MATH_BESSEL_KV_HPP
