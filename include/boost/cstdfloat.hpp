///////////////////////////////////////////////////////////////////////////////
// Copyright Christopher Kormanyos 2014.
// Copyright John Maddock 2014.
// Copyright Paul Bristow 2014.
// Distributed under the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
//

// <boost/cstdfloat.hpp> implements floating-point typedefs having
// specified widths, as described in N3626 (proposed for C++14).
// See: http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2013/n3626.pdf

#ifndef _BOOST_CSTDFLOAT_2014_01_09_HPP_
  #define _BOOST_CSTDFLOAT_2014_01_09_HPP_

  #include <boost/detail/cstdfloat_types.hpp>

  #if !defined(BOOST_CSTDFLOAT_NO_LIBQUADMATH_LIMITS)
  #include <boost/detail/cstdfloat_limits.hpp>
  #endif

  // Here is a helper function used for raising a value of a given
  // floating-point type to the power of n, where n has integral type.
  namespace boost { namespace cstdfloat { namespace detail {
  template<class float_type, class type_n> inline float_type pown(const float_type& cb, const type_n p)
  {
    if     (p <  static_cast<type_n>(0)) { return 1 / pown(cb, static_cast<type_n>(-p)); }
    else if(p == static_cast<type_n>(0)) { return float_type(1); }
    else if(p == static_cast<type_n>(1)) { return  cb; }
    else if(p == static_cast<type_n>(2)) { return  cb * cb; }
    else if(p == static_cast<type_n>(3)) { return (cb * cb) * cb; }
    else
    {
      float_type value = cb;

      type_n n;

      for(n = static_cast<type_n>(1); n <= static_cast<type_n>(p / 2); n *= 2)
      {
        value *= value;
      }

      const type_n p_minus_n = static_cast<type_n>(p - n);

      // Call the function recursively for computing the remaining power of n.
      return ((p_minus_n == static_cast<type_n>(0)) ? value : (value * pown(cb, p_minus_n)));
    }
  }
  } } } // boost::cstdfloat::detail

  #if !defined(BOOST_CSTDFLOAT_NO_LIBQUADMATH_CMATH)
  #include <boost/detail/cstdfloat_cmath.hpp>
  #endif
  #if !defined(BOOST_CSTDFLOAT_NO_LIBQUADMATH_IOSTREAM)
  #include <boost/detail/cstdfloat_iostream.hpp>
  #endif
  #if !defined(BOOST_CSTDFLOAT_NO_LIBQUADMATH_COMPLEX)
  #include <boost/detail/cstdfloat_complex.hpp>
  #endif

#endif // _BOOST_CSTDFLOAT_2014_01_09_HPP_
