//  (C) Copyright Nick Thompson 2018.
//  (C) Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_RANDOM_VECTOR
#define BOOST_MATH_TOOLS_RANDOM_VECTOR

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
#include <vector>
#include <random>
#include <type_traits>
#include <limits>
#include <cstddef>

using boost::multiprecision::cpp_bin_float_50;
using boost::multiprecision::cpp_complex_50;

// To stress test, set global_seed = 0, global_size = huge.
static const std::size_t global_seed = 0;
static const std::size_t global_size = 128;

template<bool B, class T = void>
using enable_if_t = typename std::enable_if<B, T>::type;

template<class T, enable_if_t<std::is_floating_point<T>::value, bool> = true>
std::vector<T> generate_random_vector(std::size_t size, std::size_t seed)
{
    if (seed == 0)
    {
        std::random_device rd;
        seed = rd();
    }
    std::vector<T> v(size);

    std::mt19937 gen(seed);

    std::normal_distribution<T> dis(0, 1);
    for(std::size_t i = 0; i < v.size(); ++i)
    {
        v[i] = dis(gen);
    }
    return v;
}

template<class T, enable_if_t<std::is_integral<T>::value, bool> = true>
std::vector<T> generate_random_vector(std::size_t size, std::size_t seed)
{
    if (seed == 0)
    {
        std::random_device rd;
        seed = rd();
    }
    std::vector<T> v(size);

    std::mt19937 gen(seed);

    // Rescaling by larger than 2 is UB!
    std::uniform_int_distribution<T> dis(std::numeric_limits<T>::lowest()/2, (std::numeric_limits<T>::max)()/2);
    for (std::size_t i = 0; i < v.size(); ++i)
    {
        v[i] = dis(gen);
    }
    return v;
}

template<class T, enable_if_t<boost::is_complex<T>::value, bool> = true>
std::vector<T> generate_random_vector(std::size_t size, std::size_t seed)
{
    if (seed == 0)
    {
        std::random_device rd;
        seed = rd();
    }
    std::vector<T> v(size);

    std::mt19937 gen(seed);
    
    std::normal_distribution<typename T::value_type> dis(0, 1);
    for (std::size_t i = 0; i < v.size(); ++i)
    {
        v[i] = {dis(gen), dis(gen)};
    }
    return v;  
}

template<class T, enable_if_t<std::is_same<T, cpp_complex_50>::value , bool> = true>
std::vector<T> generate_random_vector(std::size_t size, std::size_t seed)
{
    if (seed == 0)
    {
        std::random_device rd;
        seed = rd();
    }
    std::vector<T> v(size);

    std::mt19937 gen(seed);
    
    std::normal_distribution<long double> dis(0, 1);
    for (std::size_t i = 0; i < v.size(); ++i)
    {
        v[i] = {dis(gen), dis(gen)};
    }
    return v;
}

template<class T, enable_if_t<std::is_same<T, cpp_bin_float_50>::value , bool> = true>
std::vector<T> generate_random_vector(std::size_t size, std::size_t seed)
{
    if (seed == 0)
    {
        std::random_device rd;
        seed = rd();
    }
    std::vector<T> v(size);

    std::mt19937 gen(seed);

    std::normal_distribution<long double> dis(0, 1);
    for (std::size_t i = 0; i < v.size(); ++i)
    {
        v[i] = dis(gen);
    }
    return v;
}

#endif // BOOST_MATH_TOOLS_RANDOM_VECTOR
