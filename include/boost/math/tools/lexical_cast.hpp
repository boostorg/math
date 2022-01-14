//  (C) Copyright Matt Borland 2021 - 2022
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_LEXICAL_CAST
#define BOOST_MATH_TOOLS_LEXICAL_CAST

#include <cstdlib>
#include <string>
#include <type_traits>
#include <boost/math/tools/is_standalone.hpp>

#ifndef BOOST_MATH_STANDALONE
#include <boost/lexical_cast.hpp>
#else

#ifndef BOOST_MATH_NO_LEXICAL_CAST
#define BOOST_MATH_NO_LEXICAL_CAST
#endif

namespace boost {

    template <typename T, typename std::enable_if<std::is_same<T, float>::value, bool>::type = true>
    inline T lexical_cast(const char* s)
    {
        return std::strtof(s, nullptr);
    }

    template <typename T, typename std::enable_if<std::is_same<T, double>::value, bool>::type = true>
    inline T lexical_cast(const char* s)
    {
        return std::strtod(s, nullptr);
    }

    template <typename T, typename std::enable_if<std::is_same<T, long double>::value, bool>::type = true>
    inline T lexical_cast(const char* s)
    {
        return std::strtold(s, nullptr);
    }

    template <typename T, typename std::enable_if<std::is_same<T, float>::value, bool>::type = true>
    inline T lexical_cast(const std::string& s)
    {
        return std::strtof(s.c_str(), nullptr);
    }

    template <typename T, typename std::enable_if<std::is_same<T, double>::value, bool>::type = true>
    inline T lexical_cast(const std::string& s)
    {
        return std::strtod(s.c_str(), nullptr);
    }

    template <typename T, typename std::enable_if<std::is_same<T, long double>::value, bool>::type = true>
    inline T lexical_cast(const std::string& s)
    {
        return std::strtold(s.c_str(), nullptr);
    }

    template <typename T, typename T2, typename std::enable_if<(std::is_same<T, std::string>::value && std::is_arithmetic<T2>::value), bool>::type = true>
    inline T lexical_cast(const T2 v)
    {
        return std::to_string(v);
    }

    template <typename T, typename T2, typename std::enable_if<(std::is_same<T, const char*>::value && std::is_arithmetic<T2>::value), bool>::type = true>
    inline T lexical_cast(const T2 v)
    {
        const std::string temp = std::to_string(v);
        return temp.c_str();
    }
}
#endif // BOOST_MATH_STANDALONE
#endif // BOOST_MATH_TOOLS_LEXICAL_CAST
