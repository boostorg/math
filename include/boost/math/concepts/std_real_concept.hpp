//  Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// std_real_concept is an archetype for built-in Real types.

// The main purpose in providing this type is to verify
// that std lib functions are found via a using declaration
// bringing those functions into the current scope, and not
// just because they happen to be in global scope.
//
// If ::pow is found rather than std::pow say, then the code
// will silently compile, but truncation of long doubles to
// double will cause a significant loss of precision.
// A template instantiated with std_real_concept will *only*
// compile if it std::whatever is in scope.

#include <boost/config.hpp>
#include <boost/limits.hpp>
#include <boost/math/tools/real_cast.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/math/policies/policy.hpp>

#include <ostream>
#include <istream>
#include <cmath>
#include <math.h> // fmodl

#ifndef BOOST_MATH_STD_REAL_CONCEPT_HPP
#define BOOST_MATH_STD_REAL_CONCEPT_HPP

namespace boost{ namespace math{

namespace concepts
{

class std_real_concept
{
public:
   // Constructors:
   std_real_concept() : m_value(0){}
   std_real_concept(char c) : m_value(c){}
#ifndef BOOST_NO_INTRINSIC_WCHAR_T
   std_real_concept(wchar_t c) : m_value(c){}
#endif
   std_real_concept(unsigned char c) : m_value(c){}
   std_real_concept(signed char c) : m_value(c){}
   std_real_concept(unsigned short c) : m_value(c){}
   std_real_concept(short c) : m_value(c){}
   std_real_concept(unsigned int c) : m_value(c){}
   std_real_concept(int c) : m_value(c){}
   std_real_concept(unsigned long c) : m_value(c){}
   std_real_concept(long c) : m_value(c){}
#if defined(BOOST_HAS_LONG_LONG) || defined(__DECCXX) || defined(__SUNPRO_CC)
   std_real_concept(unsigned long long c) : m_value(static_cast<long double>(c)){}
   std_real_concept(long long c) : m_value(static_cast<long double>(c)){}
#endif
   std_real_concept(float c) : m_value(c){}
   std_real_concept(double c) : m_value(c){}
   std_real_concept(long double c) : m_value(c){}

   // Assignment:
   std_real_concept& operator=(char c) { m_value = c; return *this; }
   std_real_concept& operator=(unsigned char c) { m_value = c; return *this; }
   std_real_concept& operator=(signed char c) { m_value = c; return *this; }
#ifndef BOOST_NO_INTRINSIC_WCHAR_T
   std_real_concept& operator=(wchar_t c) { m_value = c; return *this; }
#endif
   std_real_concept& operator=(short c) { m_value = c; return *this; }
   std_real_concept& operator=(unsigned short c) { m_value = c; return *this; }
   std_real_concept& operator=(int c) { m_value = c; return *this; }
   std_real_concept& operator=(unsigned int c) { m_value = c; return *this; }
   std_real_concept& operator=(long c) { m_value = c; return *this; }
   std_real_concept& operator=(unsigned long c) { m_value = c; return *this; }
#if defined(BOOST_HAS_LONG_LONG) || defined(__DECCXX) || defined(__SUNPRO_CC)
   std_real_concept& operator=(long long c) { m_value = static_cast<long double>(c); return *this; }
   std_real_concept& operator=(unsigned long long c) { m_value = static_cast<long double>(c); return *this; }
#endif
   std_real_concept& operator=(float c) { m_value = c; return *this; }
   std_real_concept& operator=(double c) { m_value = c; return *this; }
   std_real_concept& operator=(long double c) { m_value = c; return *this; }

   // Access:
   long double value()const{ return m_value; }

   // Member arithmetic:
   std_real_concept& operator+=(const std_real_concept& other)
   { m_value += other.value(); return *this; }
   std_real_concept& operator-=(const std_real_concept& other)
   { m_value -= other.value(); return *this; }
   std_real_concept& operator*=(const std_real_concept& other)
   { m_value *= other.value(); return *this; }
   std_real_concept& operator/=(const std_real_concept& other)
   { m_value /= other.value(); return *this; }
   std_real_concept operator-()const
   { return -m_value; }
   std_real_concept const& operator+()const
   { return *this; }

private:
   long double m_value;
};

// Non-member arithmetic:
inline std_real_concept operator+(const std_real_concept& a, const std_real_concept& b)
{
   std_real_concept result(a);
   result += b;
   return result;
}
inline std_real_concept operator-(const std_real_concept& a, const std_real_concept& b)
{
   std_real_concept result(a);
   result -= b;
   return result;
}
inline std_real_concept operator*(const std_real_concept& a, const std_real_concept& b)
{
   std_real_concept result(a);
   result *= b;
   return result;
}
inline std_real_concept operator/(const std_real_concept& a, const std_real_concept& b)
{
   std_real_concept result(a);
   result /= b;
   return result;
}

// Comparison:
inline bool operator == (const std_real_concept& a, const std_real_concept& b)
{ return a.value() == b.value(); }
inline bool operator != (const std_real_concept& a, const std_real_concept& b)
{ return a.value() != b.value();}
inline bool operator < (const std_real_concept& a, const std_real_concept& b)
{ return a.value() < b.value(); }
inline bool operator <= (const std_real_concept& a, const std_real_concept& b)
{ return a.value() <= b.value(); }
inline bool operator > (const std_real_concept& a, const std_real_concept& b)
{ return a.value() > b.value(); }
inline bool operator >= (const std_real_concept& a, const std_real_concept& b)
{ return a.value() >= b.value(); }

#if 0
// Non-member mixed compare:
template <class T>
inline bool operator == (const T& a, const std_real_concept& b)
{
   return a == b.value();
}
template <class T>
inline bool operator != (const T& a, const std_real_concept& b)
{
   return a != b.value();
}
template <class T>
inline bool operator < (const T& a, const std_real_concept& b)
{
   return a < b.value();
}
template <class T>
inline bool operator > (const T& a, const std_real_concept& b)
{
   return a > b.value();
}
template <class T>
inline bool operator <= (const T& a, const std_real_concept& b)
{
   return a <= b.value();
}
template <class T>
inline bool operator >= (const T& a, const std_real_concept& b)
{
   return a >= b.value();
}
#endif  // Non-member mixed compare:

} // namespace concepts
} // namespace math
} // namespace boost

namespace std{

// Non-member functions:
inline boost::math::concepts::std_real_concept acos(boost::math::concepts::std_real_concept a)
{ return std::acos(a.value()); }
inline boost::math::concepts::std_real_concept cos(boost::math::concepts::std_real_concept a)
{ return std::cos(a.value()); }
inline boost::math::concepts::std_real_concept asin(boost::math::concepts::std_real_concept a)
{ return std::asin(a.value()); }
inline boost::math::concepts::std_real_concept atan(boost::math::concepts::std_real_concept a)
{ return std::atan(a.value()); }
inline boost::math::concepts::std_real_concept atan2(boost::math::concepts::std_real_concept a, boost::math::concepts::std_real_concept b)
{ return std::atan2(a.value(), b.value()); }
inline boost::math::concepts::std_real_concept ceil(boost::math::concepts::std_real_concept a)
{ return std::ceil(a.value()); }
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
inline boost::math::concepts::std_real_concept fmod(boost::math::concepts::std_real_concept a, boost::math::concepts::std_real_concept b)
{ return fmodl(a.value(), b.value()); }
#else
inline boost::math::concepts::std_real_concept fmod(boost::math::concepts::std_real_concept a, boost::math::concepts::std_real_concept b)
{ return std::fmod(a.value(), b.value()); }
#endif
inline boost::math::concepts::std_real_concept cosh(boost::math::concepts::std_real_concept a)
{ return std::cosh(a.value()); }
inline boost::math::concepts::std_real_concept exp(boost::math::concepts::std_real_concept a)
{ return std::exp(a.value()); }
inline boost::math::concepts::std_real_concept fabs(boost::math::concepts::std_real_concept a)
{ return std::fabs(a.value()); }
inline boost::math::concepts::std_real_concept abs(boost::math::concepts::std_real_concept a)
{ return std::abs(a.value()); }
inline boost::math::concepts::std_real_concept floor(boost::math::concepts::std_real_concept a)
{ return std::floor(a.value()); }
inline boost::math::concepts::std_real_concept modf(boost::math::concepts::std_real_concept a, boost::math::concepts::std_real_concept* ipart)
{
   long double ip;
   long double result = std::modf(a.value(), &ip);
   *ipart = ip;
   return result;
}
inline boost::math::concepts::std_real_concept frexp(boost::math::concepts::std_real_concept a, int* expon)
{ return std::frexp(a.value(), expon); }
inline boost::math::concepts::std_real_concept ldexp(boost::math::concepts::std_real_concept a, int expon)
{ return std::ldexp(a.value(), expon); }
inline boost::math::concepts::std_real_concept log(boost::math::concepts::std_real_concept a)
{ return std::log(a.value()); }
inline boost::math::concepts::std_real_concept log10(boost::math::concepts::std_real_concept a)
{ return std::log10(a.value()); }
inline boost::math::concepts::std_real_concept tan(boost::math::concepts::std_real_concept a)
{ return std::tan(a.value()); }
inline boost::math::concepts::std_real_concept pow(boost::math::concepts::std_real_concept a, boost::math::concepts::std_real_concept b)
{ return std::pow(a.value(), b.value()); }
#if !defined(__SUNPRO_CC)
inline boost::math::concepts::std_real_concept pow(boost::math::concepts::std_real_concept a, int b)
{ return std::pow(a.value(), b); }
#else
inline boost::math::concepts::std_real_concept pow(boost::math::concepts::std_real_concept a, int b)
{ return std::pow(a.value(), static_cast<long double>(b)); }
#endif
inline boost::math::concepts::std_real_concept sin(boost::math::concepts::std_real_concept a)
{ return std::sin(a.value()); }
inline boost::math::concepts::std_real_concept sinh(boost::math::concepts::std_real_concept a)
{ return std::sinh(a.value()); }
inline boost::math::concepts::std_real_concept sqrt(boost::math::concepts::std_real_concept a)
{ return std::sqrt(a.value()); }
inline boost::math::concepts::std_real_concept tanh(boost::math::concepts::std_real_concept a)
{ return std::tanh(a.value()); }

} // namespace std

namespace boost{ namespace math{ namespace concepts{

// Streaming:
template <class charT, class traits>
inline std::basic_ostream<charT, traits>& operator<<(std::basic_ostream<charT, traits>& os, const std_real_concept& a)
{
   return os << a.value();
}
template <class charT, class traits>
inline std::basic_istream<charT, traits>& operator>>(std::basic_istream<charT, traits>& is, std_real_concept& a)
{
   long double v;
   is >> v;
   a = v;
   return is;
}

} // namespace concepts

namespace tools
{
// real_cast converts from T to integer and narrower floating-point types.

// Convert from T to integer types.

template <>
inline unsigned int real_cast<unsigned int, concepts::std_real_concept>(concepts::std_real_concept r)
{
   return static_cast<unsigned int>(r.value());
}

template <>
inline int real_cast<int, concepts::std_real_concept>(concepts::std_real_concept r)
{
   return static_cast<int>(r.value());
}

template <>
inline long real_cast<long, concepts::std_real_concept>(concepts::std_real_concept r)
{
   return static_cast<long>(r.value());
}

// Converts from T to narrower floating-point types, float, double & long double.

template <>
inline float real_cast<float, concepts::std_real_concept>(concepts::std_real_concept r)
{
   return static_cast<float>(r.value());
}
template <>
inline double real_cast<double, concepts::std_real_concept>(concepts::std_real_concept r)
{
   return static_cast<double>(r.value());
}
template <>
inline long double real_cast<long double, concepts::std_real_concept>(concepts::std_real_concept r)
{
   return r.value();
}

template <>
inline concepts::std_real_concept max_value<concepts::std_real_concept>(BOOST_MATH_EXPLICIT_TEMPLATE_TYPE_SPEC(concepts::std_real_concept))
{
   return max_value<long double>();
}

template <>
inline concepts::std_real_concept min_value<concepts::std_real_concept>(BOOST_MATH_EXPLICIT_TEMPLATE_TYPE_SPEC(concepts::std_real_concept))
{
   return min_value<long double>();
}

template <>
inline concepts::std_real_concept log_max_value<concepts::std_real_concept>(BOOST_MATH_EXPLICIT_TEMPLATE_TYPE_SPEC(concepts::std_real_concept))
{
   return log_max_value<long double>();
}

template <>
inline concepts::std_real_concept log_min_value<concepts::std_real_concept>(BOOST_MATH_EXPLICIT_TEMPLATE_TYPE_SPEC(concepts::std_real_concept))
{
   return log_min_value<long double>();
}

template <>
inline concepts::std_real_concept epsilon<concepts::std_real_concept>(BOOST_MATH_EXPLICIT_TEMPLATE_TYPE_SPEC(concepts::std_real_concept))
{
   return tools::epsilon<long double>();
}

template <>
inline int digits<concepts::std_real_concept>(BOOST_MATH_EXPLICIT_TEMPLATE_TYPE_SPEC(concepts::std_real_concept))
{ // Assume number of significand bits is same as long double,
  // unless std::numeric_limits<T>::is_specialized to provide digits.
   return digits<long double>();
}

} // namespace tools
} // namespace math
} // namespace boost

#endif // BOOST_MATH_STD_REAL_CONCEPT_HPP


