//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_NTL_HPP
#define BOOST_MATH_TOOLS_NTL_HPP

#include <NTL/RR.h>
#include <boost/math/tools/real_cast.hpp>
#include <boost/math/tools/precision.hpp>

namespace NTL{

   inline quad_float pow(quad_float a, quad_float b)
   {
      return to_quad_float(pow(RR(a), RR(b)));
   }

   inline quad_float pow(quad_float a, long b)
   {
      return to_quad_float(pow(RR(a), RR(b)));
   }

   inline RR pow(const RR& r, long l)
   {
      return power(r, l);
   }
   inline RR tan(const RR& a)
   {
      return sin(a)/cos(a);
   }
   inline RR frexp(RR r, int* exp)
   {
      *exp = r.e;
      r.e = 0;
      while(r >= 1)
      {
         *exp += 1;
         r.e -= 1;
      }
      while(r < 0.5)
      {
         *exp -= 1;
         r.e += 1;
      }
      BOOST_ASSERT(r < 1);
      BOOST_ASSERT(r >= 0.5);
      return r;
   }
   inline RR ldexp(RR r, int exp)
   {
      r.e += exp;
      return r;
   }
}

namespace boost{ namespace math{ namespace tools{

template <>
inline float real_cast<float, NTL::RR>(NTL::RR t)
{
   double r;
   conv(r, t);
   return static_cast<float>(r);
}
template <>
inline double real_cast<double, NTL::RR>(NTL::RR t)
{
   double r;
   conv(r, t);
   return r;
}
template <>
inline long double real_cast<long double, NTL::RR>(NTL::RR t)
{
   long double result(0), last_result(0);
   double term;
   do
   {
      term = real_cast<double>(t);
      last_result = result;
      result += term;
      t -= term;
   }while(result != last_result);
   return result;
}
template <>
inline NTL::RR real_cast<NTL::RR, NTL::RR>(NTL::RR t)
{
   return t;
}

int digits(NTL::RR const &)
{
   return NTL::RR::precision();
}

template <>
NTL::RR max_value<NTL::RR>(NTL::RR const&)
{
   static bool has_init = false;
   static NTL::RR val;
   if(!has_init)
   {
      val = 1;
      val.e = NTL_OVFBND-2;
      has_init = true;
   }
   return val;
}

template <>
NTL::RR min_value<NTL::RR>(NTL::RR const&)
{
   static bool has_init = false;
   static NTL::RR val;
   if(!has_init)
   {
      val = 1;
      val.e = -NTL_OVFBND+2;
      has_init = true;
   }
   return val;
}

template <>
NTL::RR log_max_value<NTL::RR>(NTL::RR const&)
{
   static bool has_init = false;
   static NTL::RR val;
   if(!has_init)
   {
      val = 1;
      val.e = NTL_OVFBND-2;
      val = log(val);
      has_init = true;
   }
   return val;
}

template <>
NTL::RR log_min_value<NTL::RR>(NTL::RR const&)
{
   static bool has_init = false;
   static NTL::RR val;
   if(!has_init)
   {
      val = 1;
      val.e = -NTL_OVFBND+2;
      val = log(val);
      has_init = true;
   }
   return val;
}

template <>
NTL::RR epsilon<NTL::RR>(NTL::RR const& x)
{
   static const NTL::RR val = pow(NTL::RR(2.0), NTL::RR(1-digits(x)));
   return val;
}

void setprecision(std::ostream& os, NTL::RR, int p)
{
   NTL::RR::SetOutputPrecision(p);
}

template <>
inline float real_cast<float, NTL::quad_float>(NTL::quad_float t)
{
   return to_float(t);
}
template <>
inline double real_cast<double, NTL::quad_float>(NTL::quad_float t)
{
   return to_double(t);
}
template <>
inline long double real_cast<long double, NTL::quad_float>(NTL::quad_float x)
{
   long double result = x.hi;
   result += x.lo;
   return result;
}
template <>
inline NTL::quad_float real_cast<NTL::quad_float, NTL::quad_float>(NTL::quad_float t)
{
   return t;
}

template <>
inline NTL::quad_float real_cast<NTL::quad_float, NTL::RR>(NTL::RR t)
{
   return to_quad_float(t);
}

int digits(NTL::quad_float const &)
{
   return 106;
}

template <>
NTL::quad_float max_value<NTL::quad_float>(NTL::quad_float const& f)
{
   return max_value(f.hi);
}

template <>
NTL::quad_float min_value<NTL::quad_float>(NTL::quad_float const& f)
{
   return min_value(f.hi);
}

template <>
NTL::quad_float log_max_value<NTL::quad_float>(NTL::quad_float const& f)
{
   return log_max_value(f.hi);
}

template <>
NTL::quad_float log_min_value<NTL::quad_float>(NTL::quad_float const& f)
{
   return log_min_value(f.hi);
}

template <>
NTL::quad_float epsilon<NTL::quad_float>(NTL::quad_float const& x)
{
   static const NTL::quad_float val = pow(NTL::quad_float(2.0), NTL::quad_float(1-digits(x)));
   return val;
}

void setprecision(std::ostream& os, NTL::quad_float, int p)
{
   NTL::quad_float::SetOutputPrecision(p);
}

}}} // namespaces


#endif // BOOST_MATH_TOOLS_NTL_HPP



