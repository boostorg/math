//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <ostream>
#include <istream>
#include <cmath>
#include <boost/config.hpp>
#include <boost/limits.hpp>
#include <boost/math/tools/real_cast.hpp>
#include <boost/math/tools/precision.hpp>

#ifndef BOOST_MATH_REAL_CONCEPT_HPP
#define BOOST_MATH_REAL_CONCEPT_HPP

namespace boost{ namespace math{ namespace concepts{

class real_concept
{
public:
   // constructors:
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

   // assignment:
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

   // access:
   long double value()const{ return m_value; }

   // member arithmetic:
   real_concept& operator+=(const real_concept& other)
   { m_value += other.value(); return *this; }
   real_concept& operator-=(const real_concept& other)
   { m_value -= other.value(); return *this; }
   real_concept& operator*=(const real_concept& other)
   { m_value *= other.value(); return *this; }
   real_concept& operator/=(const real_concept& other)
   { m_value /= other.value(); return *this; }
   real_concept operator-()
   { return -m_value; }
   real_concept& operator+()
   { return *this; }

private:
   long double m_value;
};

// Non-member arithmetic:
real_concept operator+(const real_concept& a, const real_concept& b)
{
   real_concept result(a);
   result += b;
   return result;
}
real_concept operator-(const real_concept& a, const real_concept& b)
{
   real_concept result(a);
   result -= b;
   return result;
}
real_concept operator*(const real_concept& a, const real_concept& b)
{
   real_concept result(a);
   result *= b;
   return result;
}
real_concept operator/(const real_concept& a, const real_concept& b)
{
   real_concept result(a);
   result /= b;
   return result;
}

// comparison:
bool operator == (const real_concept& a, const real_concept& b)
{ return a.value() == b.value(); }
bool operator != (const real_concept& a, const real_concept& b)
{ return a.value() != b.value();}
bool operator < (const real_concept& a, const real_concept& b)
{ return a.value() < b.value(); }
bool operator <= (const real_concept& a, const real_concept& b)
{ return a.value() <= b.value(); }
bool operator > (const real_concept& a, const real_concept& b)
{ return a.value() > b.value(); }
bool operator >= (const real_concept& a, const real_concept& b)
{ return a.value() >= b.value(); }
#if 0
// non-member mixed compare:
template <class T>
bool operator == (const T& a, const real_concept& b)
{
   return a == b.value();
}
template <class T>
bool operator != (const T& a, const real_concept& b)
{
   return a != b.value();
}
template <class T>
bool operator < (const T& a, const real_concept& b)
{
   return a < b.value();
}
template <class T>
bool operator > (const T& a, const real_concept& b)
{
   return a > b.value();
}
template <class T>
bool operator <= (const T& a, const real_concept& b)
{
   return a <= b.value();
}
template <class T>
bool operator >= (const T& a, const real_concept& b)
{
   return a >= b.value();
}
#endif
// non-member functions:
real_concept acos(real_concept a)
{ return std::acos(a.value()); }
real_concept cos(real_concept a)
{ return std::cos(a.value()); }
real_concept asin(real_concept a)
{ return std::asin(a.value()); }
real_concept atan(real_concept a)
{ return std::atan(a.value()); }
real_concept atan2(real_concept a, real_concept b)
{ return std::atan2(a.value(), b.value()); }
real_concept ceil(real_concept a)
{ return std::ceil(a.value()); }
real_concept fmod(real_concept a, real_concept b)
{ return std::fmod(a.value(), b.value()); }
real_concept cosh(real_concept a)
{ return std::cosh(a.value()); }
real_concept exp(real_concept a)
{ return std::exp(a.value()); }
real_concept fabs(real_concept a)
{ return std::fabs(a.value()); }
real_concept floor(real_concept a)
{ return std::floor(a.value()); }
real_concept modf(real_concept a, real_concept* ipart)
{ 
   long double ip;
   long double result = std::modf(a.value(), &ip); 
   *ipart = ip;
   return result;
}
real_concept frexp(real_concept a, int* expon)
{ return std::frexp(a.value(), expon); }
real_concept ldexp(real_concept a, int expon)
{ return std::ldexp(a.value(), expon); }
real_concept log(real_concept a)
{ return std::log(a.value()); }
real_concept log10(real_concept a)
{ return std::log10(a.value()); }
real_concept tan(real_concept a)
{ return std::tan(a.value()); }
real_concept pow(real_concept a, real_concept b)
{ return std::pow(a.value(), b.value()); }
real_concept pow(real_concept a, int b)
{ return std::pow(a.value(), b); }
real_concept sin(real_concept a)
{ return std::sin(a.value()); }
real_concept sinh(real_concept a)
{ return std::sinh(a.value()); }
real_concept sqrt(real_concept a)
{ return std::sqrt(a.value()); }
real_concept tanh(real_concept a)
{ return std::tanh(a.value()); }

// streaming:
template <class charT, class traits>
std::basic_ostream<charT, traits>& operator<<(std::basic_ostream<charT, traits>& os, const real_concept& a)
{
   return os << a.value();
}
template <class charT, class traits>
std::basic_istream<charT, traits>& operator>>(std::basic_istream<charT, traits>& is, real_concept& a)
{
   long double v;
   is >> v;
   a = v;
   return is;
}

} // namepace concepts

namespace tools{

int digits(const concepts::real_concept &)
{
   return std::numeric_limits<long double>::digits;
}

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
concepts::real_concept max_value<concepts::real_concept>(concepts::real_concept const& arg)
{
   return max_value<long double>(arg.value());
}

template <>
concepts::real_concept min_value<concepts::real_concept>(concepts::real_concept const& arg)
{
   return min_value<long double>(arg.value());
}

template <>
concepts::real_concept log_max_value<concepts::real_concept>(concepts::real_concept const& arg)
{
   return log_max_value<long double>(arg.value());
}

template <>
concepts::real_concept log_min_value<concepts::real_concept>(concepts::real_concept const& arg)
{
   return log_min_value<long double>(arg.value());
}

template <>
concepts::real_concept epsilon(concepts::real_concept const&)
{
   return std::numeric_limits<long double>::epsilon();
}

}

} } // namespaces

#endif // BOOST_MATH_REAL_CONCEPT_HPP


