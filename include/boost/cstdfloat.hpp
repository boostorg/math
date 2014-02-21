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

  #if !defined(BOOST_CSTDFLOAT_NO_LIBQUADMATH_CMATH)
    #include <boost/detail/cstdfloat_cmath.hpp>
  #endif

  #if !defined(BOOST_CSTDFLOAT_NO_LIBQUADMATH_IOSTREAM)
    #if defined(BOOST_CSTDFLOAT_NO_LIBQUADMATH_CMATH)
    #error You can not use <boost/detail/cstdfloat_iostream.hpp> with BOOST_CSTDFLOAT_NO_LIBQUADMATH_CMATH defined.
    #endif
    #include <boost/detail/cstdfloat_iostream.hpp>
  #endif

  #if !defined(BOOST_CSTDFLOAT_NO_LIBQUADMATH_COMPLEX)
    #if defined(BOOST_CSTDFLOAT_NO_LIBQUADMATH_CMATH)
    #error You can not use <boost/detail/cstdfloat_complex.hpp> with BOOST_CSTDFLOAT_NO_LIBQUADMATH_CMATH defined.
    #endif
    #if defined(BOOST_CSTDFLOAT_NO_LIBQUADMATH_IOSTREAM)
    #error You can not use <boost/detail/cstdfloat_complex.hpp> with BOOST_CSTDFLOAT_NO_LIBQUADMATH_IOSTREAM defined.
    #endif
    #include <boost/detail/cstdfloat_complex.hpp>
  #endif

#endif // _BOOST_CSTDFLOAT_2014_01_09_HPP_
