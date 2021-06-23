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
  #include <boost/multiprecision/complex_adaptor.hpp>

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

  // Recognizes:
  // → fundamental types
  // → optionally float128
  // → boost::multiprecision::number<boost::multiprecision::complex_adaptor             <…>,…>
  // TODO:
  // → boost::multiprecision::number<boost::multiprecision::backends::mpfr_float_backend<…>,…>

  // FIXME: this gets unwieldy long. I feel it can be much shorter.
  template<typename RealType> struct select_complex {
    using R    = typename std::decay<RealType>::type;
    static constexpr bool is_fundamental_fpnumber = std::is_same<R,float>::value or std::is_same<R,double>::value or std::is_same<R,long double>::value;
  #ifdef BOOST_MATH_USE_FLOAT128
    static constexpr bool is_float128             = std::is_same<R,boost::multiprecision::float128>::value;
  #else
    static constexpr bool is_float128             = false;
  #endif

    template<typename Real>
    struct is_multiprecision_number {
      static constexpr bool value = false;
      using make_complex = void;
    };
    template<typename Real, boost::multiprecision::expression_template_option Et>
    struct is_multiprecision_number <boost::multiprecision::number<Real,Et>> {
      static constexpr bool value = true;
      using make_complex = boost::multiprecision::number<boost::multiprecision::complex_adaptor<Real>,Et>;
    };
    /* TODO  boost::multiprecision::mpfr_float_100
     *       or maybe put it inside the make_complex defined above.
    template<typename Real, class Alloc, boost::multiprecision::expression_template_option Et>
    struct is_multiprecision_number <boost::multiprecision::number<boost::multiprecision::backends::mpfr_float_backend<Real,Alloc>,Et>> {
      static constexpr bool value = true;
      using make_complex = boost::multiprecision::number<boost::multiprecision::mpc_complex_backend<Real>,Et>;
    };
    */

    static constexpr bool is_acceptable_number = is_fundamental_fpnumber or is_float128 or is_multiprecision_number<R>::value;
    static_assert(is_acceptable_number , "Error: cannot create complex for given real type.");

    using type = typename std::conditional< is_fundamental_fpnumber,
        std::complex<R>,
  #ifdef BOOST_MATH_USE_FLOAT128
        typename std::conditional< is_float128,
          boost::multiprecision::complex128,
  #endif
          // only the boost::multiprecision::number remain, thanks to static_assert above.
          typename is_multiprecision_number<R>::make_complex
  #ifdef BOOST_MATH_USE_FLOAT128
        >::type
  #endif
      >::type;
  };
  } } } } // namespace boost::math::fft::detail

#endif

