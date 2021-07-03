///////////////////////////////////////////////////////////////////
//  Copyright Eduardo Quintana 2021
//  Copyright Janek Kozicki 2021
//  Copyright Christopher Kormanyos 2021
//  Distributed under the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt
//  or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>
#include <limits>
#include <type_traits>

#ifndef BOOST_MATH_FFT_DISCRETE_HPP
  #define BOOST_MATH_FFT_DISCRETE_HPP
  
  namespace boost { namespace math {  namespace fft {
  namespace detail {
  
  template <typename integer>
  integer modulo(integer a, integer b)
  // it computes 'a mod b' in the mathematical sense, not 'a%b'
  // it can handle a<0 values
  // precondition: b>0
  {
    return ( a%b + b ) % b;
  }
  
  template <typename integer>
  struct euclid_return
  // helper class for the extended euclidean algorithm
  {
    integer gcd, x, y;
  };
  
  template <typename integer>
  euclid_return<integer> extended_euclid(integer a, integer b)
  // extended euclidean algorithm. It solves Bezout equation for ANY integer
  // pair (a,b):
  // g = a x  + b y
  {
    using std::swap;
    static_assert(std::numeric_limits<integer>::is_signed,
                  "Signed integers are needed");
    integer x = 1, y = 0, g = a;
    for (integer x1 = 0, y1 = 1, g1 = b; g1;)
    {
      integer q = g / g1;
      x = x - x1 * q;
      swap(x, x1);
      y = y - y1 * q;
      swap(y, y1);
      g = g - g1 * q;
      swap(g, g1);
    }
    return {g, x, y};
  }
  
  template <typename integer>
  integer gcd(integer a, integer b)
  // Greatest Commond Divisor
  {
    euclid_return<integer> ret = extended_euclid(a,b);
    return ret.gcd;
  }
  
  template <typename T,typename integer>
  T power(const T& x, integer n)
  // precondition: n>0
  {
    /* WARNING: the case n=0 works only if T{1} means the multiplicative neutral
     * element in T */
    if(n<=0)
      return T{1};
      
    bool identity = true;
    T r{};

    for (T aux{x}; n; n /= 2)
    {
      if (n % 2)
      {
        r = identity ? aux : r * aux;
        identity = false;
      }
      aux *= aux;
    }
    return r;
  }
  
  template<typename integer>
  constexpr 
  typename std::enable_if<std::is_integral<integer>::value,bool>::type 
  // enabled only for builtin integer types
  is_power2(integer x) 
  { 
    return (x>0) && (x == (x & -x));
  }
  
  template<typename integer>
  constexpr 
  typename std::enable_if<std::is_integral<integer>::value,integer>::type 
  // enabled only for builtin integer types
  lower_bound_power2(integer x)
  // Returns the biggest power of two that is smaller than or equal to x.
  // returns 1 if no such power of two exists, ie. for x<=0
  {
    if(x<=0)
      return 0;
      
    for(long lsb = x & -x; lsb!=x;)
    {
      x-=lsb;
      lsb = x&-x;
    }
    return x;
  }
  
  
  template<typename integer>
  constexpr 
  typename std::enable_if<std::is_integral<integer>::value,integer>::type 
  // enabled only for builtin integer types
  upper_bound_power2(integer x)
  // Returns the smallest power of two that is greater or equal to x.
  // returns 0 if no such power of two exits, ie. when the value overflows
  // the integer
  {
    if(x<=1)
      return 1;
    
    constexpr integer up_max = lower_bound_power2(
      std::numeric_limits<integer>::max() );
    if(up_max<x)
      return 0;
    
    integer up=1;
    for(;up<x;up<<=1);
    return up;
  }

  } // namespace detail
  } } } // namespace boost::math::fft
#endif
