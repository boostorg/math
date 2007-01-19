//  Copyright (c) 2006 Xiaogang Zhang
//
//  This code may be used under either of the following two licences:
//
//  Permission is hereby granted, free of charge, to any person obtaining a copy
//  of this software and associated documentation files (the "Software"), to deal
//  in the Software without restriction, including without limitation the rights
//  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//  copies of the Software, and to permit persons to whom the Software is
//  furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included in
//  all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
//  THE SOFTWARE. OF SUCH DAMAGE.
//
//  Or:
//
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_BESSEL_JV_HPP
#define BOOST_MATH_BESSEL_JV_HPP

#include <boost/math/special_functions/bessel_jn.hpp>
#include <boost/math/special_functions/bessel_jy.hpp>
#include <boost/math/constants/constants.hpp>

// Bessel and Spherical Bessel functions of the first kind

namespace boost { namespace math {

// Bessel function of the first kind of any order
template <typename T>
inline T bessel_jv(T v, T x)
{
    int n = static_cast<int>(v);
    if (v == n)
    {
        return bessel_jn(n, x);             // v is integer
    }
    else
    {
        T J, Y;
        bessel_jy(v, x, &J, &Y);
        return J;
    }
}

// Spherical Bessel function of the first kind of non-negative order
template <typename T>
inline T bessel_sjv(unsigned n, T x)
{
    using namespace std;
    using namespace boost::math::constants;

    T v = n + 0.5L;
    if (x == 0)
    {
        return (n == 0) ? static_cast<T>(1) : static_cast<T>(0);
    }
    else
    {
        return sqrt(0.5L * pi<T>() / x) * bessel_jv(v, x);
    }
}

// -------------------- TR1 functions --------------------

inline float cyl_bessel_jf(float nu, float x)
{
    if (nu >= 128)
    {
        std::cout << "Warning: cyl_bessel_jf(nu, x), nu >= 128, "
                  << "result is implementation defined according to C++ TR1"
                  << std::endl;
    }
    return bessel_jv<float>(nu, x);
}

inline double cyl_bessel_j(double nu, double x)
{
    if (nu >= 128)
    {
        std::cout << "Warning: cyl_bessel_j(nu, x), nu >= 128, "
                  << "result is implementation defined according to C++ TR1"
                  << std::endl;
    }
    return bessel_jv<double>(nu, x);
}

inline long double cyl_bessel_jl(long double nu, long double x)
{
    if (nu >= 128)
    {
        std::cout << "Warning: cyl_bessel_jl(nu, x), nu >= 128, "
                  << "result is implementation defined according to C++ TR1"
                  << std::endl;
    }
    return bessel_jv<long double>(nu, x);
}

inline float sph_besself(unsigned n, float x)
{
    if (n >= 128)
    {
        std::cout << "Warning: sph_besself(n, x), n >= 128, "
                  << "result is implementation defined according to C++ TR1"
                  << std::endl;
    }
    return bessel_sjv<float>(n, x);
}

inline double sph_bessel(unsigned n, double x)
{
    if (n >= 128)
    {
        std::cout << "Warning: sph_bessel(n, x), n >= 128, "
                  << "result is implementation defined according to C++ TR1"
                  << std::endl;
    }
    return bessel_sjv<double>(n, x);
}

inline long double sph_bessell(unsigned n, long double x)
{
    if (n >= 128)
    {
        std::cout << "Warning: sph_bessell(n, x), n >= 128, "
                  << "result is implementation defined according to C++ TR1"
                  << std::endl;
    }
    return bessel_sjv<long double>(n, x);
}

}} // namespaces

#endif // BOOST_MATH_BESSEL_JV_HPP
