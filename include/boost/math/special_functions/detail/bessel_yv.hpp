//  Copyright (c) 2006 Xiaogang Zhang
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_BESSEL_YV_HPP
#define BOOST_MATH_BESSEL_YV_HPP

#include <boost/math/special_functions/bessel_yn.hpp>
#include <boost/math/special_functions/bessel_jy.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/error_handling.hpp>

// Bessel and Spherical Bessel functions of the second kind

namespace boost { namespace math {

// Bessel function of the second kind of any order
template <typename T>
inline T bessel_yv(T v, T x)
{
    int n = static_cast<int>(v);
    if (v == n)
    {
        return bessel_yn(n, x);             // v is integer
    }
    else
    {
        T J, Y;
        bessel_jy(v, x, &J, &Y);
        return Y;
    }
}

// Spherical Bessel function of the second kind of non-negative order
template <typename T>
inline T bessel_syv(unsigned n, T x)
{
    using namespace std;
    using namespace boost::math::tools;
    using namespace boost::math::constants;

    T v = n + 0.5L;
    if (x == 0)
    {
        return -overflow_error<T>("boost::math::bessel_syv(n, x)",
            "infinity occurred but not supported");
    }
    else
    {
        return sqrt(0.5L * pi<T>() / x) * bessel_yv(v, x);
    }
}

// -------------------- TR1 functions --------------------

inline float cyl_neumannf(float nu, float x)
{
    if (nu >= 128)
    {
        std::cout << "Warning: cyl_neumannf(nu, x), nu >= 128, "
                  << "result is implementation defined according to C++ TR1"
                  << std::endl;
    }
    return bessel_yv<float>(nu, x);
}

inline double cyl_neumann(double nu, double x)
{
    if (nu >= 128)
    {
        std::cout << "Warning: cyl_neumann(nu, x), nu >= 128, "
                  << "result is implementation defined according to C++ TR1"
                  << std::endl;
    }
    return bessel_yv<double>(nu, x);
}

inline long double cyl_neumannl(long double nu, long double x)
{
    if (nu >= 128)
    {
        std::cout << "Warning: cyl_neumannl(nu, x), nu >= 128, "
                  << "result is implementation defined according to C++ TR1"
                  << std::endl;
    }
    return bessel_yv<long double>(nu, x);
}

inline float sph_neumannf(unsigned n, float x)
{
    if (n >= 128)
    {
        std::cout << "Warning: sph_neumannf(n, x), n >= 128, "
                  << "result is implementation defined according to C++ TR1"
                  << std::endl;
    }
    return bessel_syv<float>(n, x);
}

inline double sph_neumann(unsigned n, double x)
{
    if (n >= 128)
    {
        std::cout << "Warning: sph_neumann(n, x), n >= 128, "
                  << "result is implementation defined according to C++ TR1"
                  << std::endl;
    }
    return bessel_syv<double>(n, x);
}

inline long double sph_neumannl(unsigned n, long double x)
{
    if (n >= 128)
    {
        std::cout << "Warning: sph_neumannl(n, x), n >= 128, "
                  << "result is implementation defined according to C++ TR1"
                  << std::endl;
    }
    return bessel_syv<long double>(n, x);
}

}} // namespaces

#endif // BOOST_MATH_BESSEL_YV_HPP
