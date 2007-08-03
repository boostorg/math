//  Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Test real concept.

// real_concept is an archetype for User defined Real types.

// This file defines the features, constructors, operators, functions...
// that are essential to use mathematical and statistical functions.
// The template typename "RealType" is used where this type
// (as well as the normal built-in types, float, double & long double)
// can be used.
// That this is the minimum set is confirmed by use as a type
// in tests of all functions & distributions, for example:
//   test_spots(0.F); & test_spots(0.);  for float and double, but also
//   test_spots(boost::math::concepts::real_concept(0.));
// NTL quad_float type is an example of a type meeting the requirements,
// but note minor additions are needed - see ntl.diff and documentation
// "Using With NTL - a High-Precision Floating-Point Library".

#include <boost/config.hpp>
#include <boost/limits.hpp>
#include <boost/math/tools/real_cast.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/math/policy/policy.hpp>

#include <ostream>
#include <istream>
#include <cmath>

#ifndef BOOST_MATH_REAL_CONCEPT_HPP
#define BOOST_MATH_REAL_CONCEPT_HPP

namespace boost{ namespace math{

namespace concepts
{

class real_concept
{
public:
   // Constructors:
   real_concept() : m_value(0){}
   real_concept(char c) : m_value(c){}
#ifndef BOOST_NO_INTRINSIC_WCHAR_T
   real_concept(wchar_t c) : m_value(c){}
#endif
   real_concept(unsigned char c) : m_value(c){}
   real_concept(signed char c) : m_value(c){}
   real_concept(unsigned short c) : m_value(c){}
   real_concept(short c) : m_value(c){}
   real_concept(unsigned int c) : m_value(c){}
   real_concept(int c) : m_value(c){}
   real_concept(unsigned long c) : m_value(c){}
   real_concept(long c) : m_value(c){}
#ifdef BOOST_HAS_LONG_LONG
   real_concept(unsigned long long c) : m_value(static_cast<long double>(c)){}
   real_concept(long long c) : m_value(static_cast<long double>(c)){}
#endif
   real_concept(float c) : m_value(c){}
   real_concept(double c) : m_value(c){}
   real_concept(long double c) : m_value(c){}

   // Assignment:
   real_concept& operator=(char c) { m_value = c; return *this; }
   real_concept& operator=(unsigned char c) { m_value = c; return *this; }
   real_concept& operator=(signed char c) { m_value = c; return *this; }
#ifndef BOOST_NO_INTRINSIC_WCHAR_T
   real_concept& operator=(wchar_t c) { m_value = c; return *this; }
#endif
   real_concept& operator=(short c) { m_value = c; return *this; }
   real_concept& operator=(unsigned short c) { m_value = c; return *this; }
   real_concept& operator=(int c) { m_value = c; return *this; }
   real_concept& operator=(unsigned int c) { m_value = c; return *this; }
   real_concept& operator=(long c) { m_value = c; return *this; }
   real_concept& operator=(unsigned long c) { m_value = c; return *this; }
#ifdef BOOST_HAS_LONG_LONG
   real_concept& operator=(long long c) { m_value = static_cast<long double>(c); return *this; }
   real_concept& operator=(unsigned long long c) { m_value = static_cast<long double>(c); return *this; }
#endif
   real_concept& operator=(float c) { m_value = c; return *this; }
   real_concept& operator=(double c) { m_value = c; return *this; }
   real_concept& operator=(long double c) { m_value = c; return *this; }

   // Access:
   long double value()const{ return m_value; }

   // Member arithmetic:
   real_concept& operator+=(const real_concept& other)
   { m_value += other.value(); return *this; }
   real_concept& operator-=(const real_concept& other)
   { m_value -= other.value(); return *this; }
   real_concept& operator*=(const real_concept& other)
   { m_value *= other.value(); return *this; }
   real_concept& operator/=(const real_concept& other)
   { m_value /= other.value(); return *this; }
   real_concept operator-()const
   { return -m_value; }
   real_concept const& operator+()const
   { return *this; }

private:
   long double m_value;
};

// Non-member arithmetic:
inline real_concept operator+(const real_concept& a, const real_concept& b)
{
   real_concept result(a);
   result += b;
   return result;
}
inline real_concept operator-(const real_concept& a, const real_concept& b)
{
   real_concept result(a);
   result -= b;
   return result;
}
inline real_concept operator*(const real_concept& a, const real_concept& b)
{
   real_concept result(a);
   result *= b;
   return result;
}
inline real_concept operator/(const real_concept& a, const real_concept& b)
{
   real_concept result(a);
   result /= b;
   return result;
}

// Comparison:
inline bool operator == (const real_concept& a, const real_concept& b)
{ return a.value() == b.value(); }
inline bool operator != (const real_concept& a, const real_concept& b)
{ return a.value() != b.value();}
inline bool operator < (const real_concept& a, const real_concept& b)
{ return a.value() < b.value(); }
inline bool operator <= (const real_concept& a, const real_concept& b)
{ return a.value() <= b.value(); }
inline bool operator > (const real_concept& a, const real_concept& b)
{ return a.value() > b.value(); }
inline bool operator >= (const real_concept& a, const real_concept& b)
{ return a.value() >= b.value(); }

#if 0
// Non-member mixed compare:
template <class T>
inline bool operator == (const T& a, const real_concept& b)
{
   return a == b.value();
}
template <class T>
inline bool operator != (const T& a, const real_concept& b)
{
   return a != b.value();
}
template <class T>
inline bool operator < (const T& a, const real_concept& b)
{
   return a < b.value();
}
template <class T>
inline bool operator > (const T& a, const real_concept& b)
{
   return a > b.value();
}
template <class T>
inline bool operator <= (const T& a, const real_concept& b)
{
   return a <= b.value();
}
template <class T>
inline bool operator >= (const T& a, const real_concept& b)
{
   return a >= b.value();
}
#endif  // Non-member mixed compare:

// Non-member functions:
inline real_concept acos(real_concept a)
{ return std::acos(a.value()); }
inline real_concept cos(real_concept a)
{ return std::cos(a.value()); }
inline real_concept asin(real_concept a)
{ return std::asin(a.value()); }
inline real_concept atan(real_concept a)
{ return std::atan(a.value()); }
inline real_concept atan2(real_concept a, real_concept b)
{ return std::atan2(a.value(), b.value()); }
inline real_concept ceil(real_concept a)
{ return std::ceil(a.value()); }
inline real_concept fmod(real_concept a, real_concept b)
{ return std::fmod(a.value(), b.value()); }
inline real_concept cosh(real_concept a)
{ return std::cosh(a.value()); }
inline real_concept exp(real_concept a)
{ return std::exp(a.value()); }
inline real_concept fabs(real_concept a)
{ return std::fabs(a.value()); }
inline real_concept abs(real_concept a)
{ return std::abs(a.value()); }
inline real_concept floor(real_concept a)
{ return std::floor(a.value()); }
inline real_concept modf(real_concept a, real_concept* ipart)
{
   long double ip;
   long double result = std::modf(a.value(), &ip);
   *ipart = ip;
   return result;
}
inline real_concept frexp(real_concept a, int* expon)
{ return std::frexp(a.value(), expon); }
inline real_concept ldexp(real_concept a, int expon)
{ return std::ldexp(a.value(), expon); }
inline real_concept log(real_concept a)
{ return std::log(a.value()); }
inline real_concept log10(real_concept a)
{ return std::log10(a.value()); }
inline real_concept tan(real_concept a)
{ return std::tan(a.value()); }
inline real_concept pow(real_concept a, real_concept b)
{ return std::pow(a.value(), b.value()); }
inline real_concept pow(real_concept a, int b)
{ return std::pow(a.value(), b); }
inline real_concept sin(real_concept a)
{ return std::sin(a.value()); }
inline real_concept sinh(real_concept a)
{ return std::sinh(a.value()); }
inline real_concept sqrt(real_concept a)
{ return std::sqrt(a.value()); }
inline real_concept tanh(real_concept a)
{ return std::tanh(a.value()); }

// Streaming:
template <class charT, class traits>
inline std::basic_ostream<charT, traits>& operator<<(std::basic_ostream<charT, traits>& os, const real_concept& a)
{
   return os << a.value();
}
template <class charT, class traits>
inline std::basic_istream<charT, traits>& operator>>(std::basic_istream<charT, traits>& is, real_concept& a)
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
inline unsigned int real_cast<unsigned int, concepts::real_concept>(concepts::real_concept r)
{
   return static_cast<unsigned int>(r.value());
}

template <>
inline int real_cast<int, concepts::real_concept>(concepts::real_concept r)
{
   return static_cast<int>(r.value());
}

template <>
inline long real_cast<long, concepts::real_concept>(concepts::real_concept r)
{
   return static_cast<long>(r.value());
}

// Converts from T to narrower floating-point types, float, double & long double.

template <>
inline float real_cast<float, concepts::real_concept>(concepts::real_concept r)
{
   return static_cast<float>(r.value());
}
template <>
inline double real_cast<double, concepts::real_concept>(concepts::real_concept r)
{
   return static_cast<double>(r.value());
}
template <>
inline long double real_cast<long double, concepts::real_concept>(concepts::real_concept r)
{
   return r.value();
}

template <>
inline concepts::real_concept max_value<concepts::real_concept>(BOOST_EXPLICIT_TEMPLATE_TYPE_SPEC(concepts::real_concept))
{
   return max_value<long double>();
}

template <>
inline concepts::real_concept min_value<concepts::real_concept>(BOOST_EXPLICIT_TEMPLATE_TYPE_SPEC(concepts::real_concept))
{
   return min_value<long double>();
}

template <>
inline concepts::real_concept log_max_value<concepts::real_concept>(BOOST_EXPLICIT_TEMPLATE_TYPE_SPEC(concepts::real_concept))
{
   return log_max_value<long double>();
}

template <>
inline concepts::real_concept log_min_value<concepts::real_concept>(BOOST_EXPLICIT_TEMPLATE_TYPE_SPEC(concepts::real_concept))
{
   return log_min_value<long double>();
}

template <>
inline concepts::real_concept epsilon(BOOST_EXPLICIT_TEMPLATE_TYPE_SPEC(concepts::real_concept))
{
   return tools::epsilon<long double>();
}

} // namespace tools

namespace policies{

template <>
inline int digits<concepts::real_concept, policy<> >(BOOST_EXPLICIT_TEMPLATE_TYPE_SPEC(concepts::real_concept))
{ // Assume number of significand bits is same as long double,
  // unless std::numeric_limits<T>::is_specialized to provide digits.
   return boost::math::policies::digits<long double, boost::math::policies::policy<> >();
   // Note that if numeric_limits real concept is NOT specialized to provide digits10
   // (or max_digits10) then the default precision of 6 decimal digits will be used
   // by Boost test (giving misleading error messages like
   // "difference between {9.79796} and {9.79796} exceeds 5.42101e-19%"
   // and by Boost lexical cast and serialization causing loss of accuracy.
}

template <>
inline int digits<concepts::real_concept, policy<detail::forwarding_arg1, detail::forwarding_arg2 > >(BOOST_EXPLICIT_TEMPLATE_TYPE_SPEC(concepts::real_concept))
{ return digits<concepts::real_concept, policy<> >(); }

template <>
inline int digits<concepts::real_concept, policy<discrete_quantile<real> > >(BOOST_EXPLICIT_TEMPLATE_TYPE_SPEC(concepts::real_concept))
{ return digits<concepts::real_concept, policy<> >(); }

template <>
inline int digits<concepts::real_concept, policy<discrete_quantile<integer_below> > >(BOOST_EXPLICIT_TEMPLATE_TYPE_SPEC(concepts::real_concept))
{ return digits<concepts::real_concept, policy<> >(); }

template <>
inline int digits<concepts::real_concept, policy<discrete_quantile<integer_above> > >(BOOST_EXPLICIT_TEMPLATE_TYPE_SPEC(concepts::real_concept))
{ return digits<concepts::real_concept, policy<> >(); }

template <>
inline int digits<concepts::real_concept, policy<discrete_quantile<integer_outside> > >(BOOST_EXPLICIT_TEMPLATE_TYPE_SPEC(concepts::real_concept))
{ return digits<concepts::real_concept, policy<> >(); }

template <>
inline int digits<concepts::real_concept, policy<discrete_quantile<integer_inside> > >(BOOST_EXPLICIT_TEMPLATE_TYPE_SPEC(concepts::real_concept))
{ return digits<concepts::real_concept, policy<> >(); }

template <>
inline int digits<concepts::real_concept, policy<discrete_quantile<integer_nearest> > >(BOOST_EXPLICIT_TEMPLATE_TYPE_SPEC(concepts::real_concept))
{ return digits<concepts::real_concept, policy<> >(); }

}

} // namespace math
} // namespace boost

#endif // BOOST_MATH_REAL_CONCEPT_HPP


