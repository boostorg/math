///////////////////////////////////////////////////////////////////////////////
// Copyright Christopher Kormanyos 2014.
// Copyright John Maddock 2014.
// Copyright Paul Bristow 2014.
// Distributed under the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
//

// Implement quadruple-precision (and extended) support for <complex>.

#ifndef _BOOST_CSTDFLOAT_COMPLEX_2014_02_15_HPP_
  #define _BOOST_CSTDFLOAT_COMPLEX_2014_02_15_HPP_

  #include <boost/detail/cstdfloat_cmath.hpp>
  #include <boost/detail/cstdfloat_iostream.hpp>

  #if defined(BOOST_CSTDFLOAT_NO_LIBQUADMATH_IOSTREAM)
  #error Quadruple-precision complex support in <boost/cstdfloat.hpp> is incompatible with BOOST_CSTDFLOAT_NO_LIBQUADMATH_IOSTREAM
  #endif

  #if defined(BOOST_CSTDFLOAT_NO_LIBQUADMATH_CMATH)
  #error Quadruple-precision complex support in <boost/cstdfloat.hpp> is incompatible with BOOST_CSTDFLOAT_NO_LIBQUADMATH_CMATH
  #endif

  #if !defined(BOOST_NO_FLOAT128_T) && defined(BOOST_MATH_USE_FLOAT128) && !defined(BOOST_CSTDFLOAT_NO_LIBQUADMATH_SUPPORT)

  #include <boost/detail/cstdfloat_cmath.hpp>

  #define BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE boost::cstdfloat::detail::float_internal128_t

  #include <boost/detail/cstdfloat_complex_std.hpp>

  #undef BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE

  #endif // Not BOOST_CSTDFLOAT_NO_LIBQUADMATH_SUPPORT (i.e., the user would like to have libquadmath support)

#endif // _BOOST_CSTDFLOAT_COMPLEX_2014_02_15_HPP_
