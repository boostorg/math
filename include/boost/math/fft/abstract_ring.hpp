///////////////////////////////////////////////////////////////////
//  Copyright Eduardo Quintana 2021
//  Copyright Janek Kozicki 2021
//  Copyright Christopher Kormanyos 2021
//  Distributed under the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt
//  or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_FFT_ABSTRACT_RING_HPP
  #define BOOST_MATH_FFT_ABSTRACT_RING_HPP
  // FIXME: only the complex rings are supported at the moment.
  #include <boost/config.hpp>
  #include <boost/cstdfloat.hpp>

  // TODO: discuss following problem: I include here a bunch of files, each
  // providing a different complex type. The preprocessor chews a lot, most of
  // the time unnecessarily.
  // possible solution: an #ifdef to check if this file was already included by the user,
  // then write in the documentation, that the type
  // ( e.g. <boost/multiprecision/cpp_complex.hpp>)
  // has to be included before <boost/math/fft.hpp>
  //
  // Another solution: specialise the templates on ComplexType, not on the RealType.
  // Then ring axioms will be easier and select_complex will become unnecessary.
  //
  #ifdef BOOST_MATH_USE_FLOAT128
    #include <boost/multiprecision/complex128.hpp>
  #endif
  #include <boost/multiprecision/cpp_complex.hpp>

  /*
    RingType axioms:
    1. Abelian group addition (operator+)
      -> closure
      -> associativity
      -> neutral element (0)
      -> inverse (operator-)
      -> commutativity
    2. Monoid multiplication (operator*)
      -> closure
      -> associativity
      -> neutral element (1)
    3. addition and multiplication compatibility
      -> left distributivity, ie. a*(b+c) == a*b + a*c
      -> right distributivity, ie. (b+c)*a == b*a + c*a

  TODO: add suppport for non - complex rings.
        possibly with verification whether the user-provided type satisfies the requested axioms.
        (best would be to verify this at compile-time, but it's unclear about whether it's possible)
  */

  namespace boost { namespace math {  namespace fft {  namespace detail {

  template<typename RealType> struct select_complex {};
  template<> struct select_complex<float> { using type = std::complex<float>; };
  template<> struct select_complex<double> { using type = std::complex<double>; };
  template<> struct select_complex<long double> { using type = std::complex<long double>; };
  #ifdef BOOST_MATH_USE_FLOAT128
  template<> struct select_complex<boost::multiprecision::float128> { using type = boost::multiprecision::complex128; };
  #endif
  // TODO: discuss - use traits or something else to recognize all boost::multiprecision complex types
  template<> struct select_complex<typename boost::multiprecision::cpp_complex_50::value_type> { using type = boost::multiprecision::cpp_complex_50; };

  } } } } // namespace boost::math::fft::detail

#endif

