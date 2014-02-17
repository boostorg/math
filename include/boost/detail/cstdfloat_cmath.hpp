///////////////////////////////////////////////////////////////////////////////
// Copyright Christopher Kormanyos 2014.
// Copyright John Maddock 2014.
// Copyright Paul Bristow 2014.
// Distributed under the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
//

// Implement quadruple-precision <cmath> support.

#ifndef _BOOST_CSTDFLOAT_CMATH_2014_02_15_HPP_
  #define _BOOST_CSTDFLOAT_CMATH_2014_02_15_HPP_

  #include <boost/detail/cstdfloat_types.hpp>

  #if !defined(BOOST_NO_FLOAT128_T) && defined(BOOST_MATH_USE_FLOAT128) && !defined(BOOST_CSTDFLOAT_NO_LIBQUADMATH_SUPPORT)

  #include <cmath>
  #include <stdexcept>
  #include <boost/cstdint.hpp>
  #include <boost/math/constants/constants.hpp>
  #include <boost/static_assert.hpp>
  #include <boost/throw_exception.hpp>

  #if defined(BOOST_INTEL)
    #define BOOST_CSTDFLOAT_FLOAT128_LDEXP  __ldexpq
    #define BOOST_CSTDFLOAT_FLOAT128_FREXP  __frexpq
    #define BOOST_CSTDFLOAT_FLOAT128_FABS   __fabsq
    #define BOOST_CSTDFLOAT_FLOAT128_FLOOR  __floorq
    #define BOOST_CSTDFLOAT_FLOAT128_CEIL   __ceilq
    #if !defined(BOOST_CSTDFLOAT_FLOAT128_SQRT)
    #define BOOST_CSTDFLOAT_FLOAT128_SQRT   __sqrtq
    #endif
    #define BOOST_CSTDFLOAT_FLOAT128_TRUNC  __truncq
    #define BOOST_CSTDFLOAT_FLOAT128_EXP    __expq_patch
    #define BOOST_CSTDFLOAT_FLOAT128_POW    __powq
    #define BOOST_CSTDFLOAT_FLOAT128_LOG    __logq
    #define BOOST_CSTDFLOAT_FLOAT128_LOG10  __log10q
    #define BOOST_CSTDFLOAT_FLOAT128_SIN    __sinq
    #define BOOST_CSTDFLOAT_FLOAT128_COS    __cosq
    #define BOOST_CSTDFLOAT_FLOAT128_TAN    __tanq
    #define BOOST_CSTDFLOAT_FLOAT128_ASIN   __asinq
    #define BOOST_CSTDFLOAT_FLOAT128_ACOS   __acosq
    #define BOOST_CSTDFLOAT_FLOAT128_ATAN   __atanq
    #define BOOST_CSTDFLOAT_FLOAT128_SINH   __sinhq_patch
    #define BOOST_CSTDFLOAT_FLOAT128_COSH   __coshq_patch
    #define BOOST_CSTDFLOAT_FLOAT128_TANH   __tanhq_patch
    #define BOOST_CSTDFLOAT_FLOAT128_FMOD   __fmodq
    #define BOOST_CSTDFLOAT_FLOAT128_ATAN2  __atan2q
    #define BOOST_CSTDFLOAT_FLOAT128_LGAMMA __lgammaq
    #define BOOST_CSTDFLOAT_FLOAT128_TGAMMA __tgammaq
  #elif defined(__GNUC__)
    #define BOOST_CSTDFLOAT_FLOAT128_LDEXP  ldexpq
    #define BOOST_CSTDFLOAT_FLOAT128_FREXP  frexpq
    #define BOOST_CSTDFLOAT_FLOAT128_FABS   fabsq
    #define BOOST_CSTDFLOAT_FLOAT128_FLOOR  floorq
    #define BOOST_CSTDFLOAT_FLOAT128_CEIL   ceilq
    #if !defined(BOOST_CSTDFLOAT_FLOAT128_SQRT)
    #define BOOST_CSTDFLOAT_FLOAT128_SQRT   sqrtq
    #endif
    #define BOOST_CSTDFLOAT_FLOAT128_TRUNC  truncq
    #define BOOST_CSTDFLOAT_FLOAT128_EXP    expq_patch
    #define BOOST_CSTDFLOAT_FLOAT128_POW    powq
    #define BOOST_CSTDFLOAT_FLOAT128_LOG    logq
    #define BOOST_CSTDFLOAT_FLOAT128_LOG10  log10q
    #define BOOST_CSTDFLOAT_FLOAT128_SIN    sinq
    #define BOOST_CSTDFLOAT_FLOAT128_COS    cosq
    #define BOOST_CSTDFLOAT_FLOAT128_TAN    tanq
    #define BOOST_CSTDFLOAT_FLOAT128_ASIN   asinq
    #define BOOST_CSTDFLOAT_FLOAT128_ACOS   acosq
    #define BOOST_CSTDFLOAT_FLOAT128_ATAN   atanq
    #define BOOST_CSTDFLOAT_FLOAT128_SINH   sinhq_patch
    #define BOOST_CSTDFLOAT_FLOAT128_COSH   coshq_patch
    #define BOOST_CSTDFLOAT_FLOAT128_TANH   tanhq_patch
    #define BOOST_CSTDFLOAT_FLOAT128_FMOD   fmodq
    #define BOOST_CSTDFLOAT_FLOAT128_ATAN2  atan2q
    #define BOOST_CSTDFLOAT_FLOAT128_LGAMMA lgammaq
    #define BOOST_CSTDFLOAT_FLOAT128_TGAMMA tgammaq
  #endif

  // Implement quadruple-precision <cmath> functions in the namespace
  // boost::cstdfloat::detail. Subsequently inject these into the
  // std namespace via *using* directive.

  // Begin with some forward function declarations.

  // Forward declarations of quadruple-precision elementary functions.
  extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_LDEXP (boost::cstdfloat::detail::float_internal128_t, int)  throw();
  extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_FREXP (boost::cstdfloat::detail::float_internal128_t, int*) throw();
  extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_FABS  (boost::cstdfloat::detail::float_internal128_t) throw();
  extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_FLOOR (boost::cstdfloat::detail::float_internal128_t) throw();
  extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_CEIL  (boost::cstdfloat::detail::float_internal128_t) throw();
  extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_SQRT  (boost::cstdfloat::detail::float_internal128_t) throw();
  extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_TRUNC (boost::cstdfloat::detail::float_internal128_t) throw();
  inline     boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_EXP   (boost::cstdfloat::detail::float_internal128_t x)
  {
    // Patch the expq() function for a subset of broken GCC compilers
    // like GCC 4.7, 4.8 on MinGW.
    typedef boost::cstdfloat::detail::float_internal128_t float_type;

    // Use an order-36 polynomial approximation of the exponential function
    // in the range of (-ln2 < x < ln2). Scale the argument to this range
    // and multiply the result by 2^n accordingly.

    // Derive the polynomial coefficients with Mathematica(R) by generating
    // a table of high-precision values in the range of (-ln2 < x < ln2)
    // and subsequently applying the *Fit* function.

    // Table[{x, Exp[x] - 1}, {x, -Log[2], Log[2], 1/180}]
    // N[%, 120]
    // Fit[%, {x, x^2, x^3, x^4, x^5, x^6, x^7, x^8, x^9, x^10, x^11, x^12,
    //         x^13, x^14, x^15, x^16, x^17, x^18, x^19, x^20, x^21, x^22,
    //         x^23, x^24, x^25, x^26, x^27, x^28, x^29, x^30, x^31, x^32,
    //         x^33, x^34, x^35, x^36}, x]

    // Scale the argument x to the range (-ln2 < x < ln2).
    BOOST_CONSTEXPR_OR_CONST float_type one_over_ln2 = float_type(BOOST_FLOAT128_C(1.44269504088896340735992468100189213742664595415299));
    const float_type x_over_ln2   = x * one_over_ln2;

    boost::int_fast32_t n;

    if     (x < -1) { n = static_cast<boost::int_fast32_t>(BOOST_CSTDFLOAT_FLOAT128_CEIL (x_over_ln2)); }
    else if(x > +1) { n = static_cast<boost::int_fast32_t>(BOOST_CSTDFLOAT_FLOAT128_FLOOR(x_over_ln2)); }
    else            { n = static_cast<boost::int_fast32_t>(0); }

    // Check if the argument is very near an integer.
    const boost::int_fast32_t nn = ((n < static_cast<boost::int_fast32_t>(0)) ? -n : n);

    if(BOOST_CSTDFLOAT_FLOAT128_FABS(x - n) < float_type(BOOST_CSTDFLOAT_FLOAT128_EPS * nn))
    {
      // Return e**n for arguments very near an integer.
      return boost::cstdfloat::detail::pown(boost::math::constants::e<float_type>(), n);
    }

    // Here, alpha is the scaled argument.
    const float_type alpha = x - (n * boost::math::constants::ln_two<float_type>());

    // Compute the polynomial approximation of the scaled-argument exponential function.
    const float_type sum =
      ((((((((((((((((((((((((((((((((((((  float_type(BOOST_FLOAT128_C(2.69291698127774166063293705964720493864630783729857438187365E-42))  * alpha
                                          + float_type(BOOST_FLOAT128_C(9.70937085471487654794114679403710456028986572118859594614033E-41))) * alpha
                                          + float_type(BOOST_FLOAT128_C(3.38715585158055097155585505318085512156885389014410753080500E-39))) * alpha
                                          + float_type(BOOST_FLOAT128_C(1.15162718532861050809222658798662695267019717760563645440433E-37))) * alpha
                                          + float_type(BOOST_FLOAT128_C(3.80039074689434663295873584133017767349635602413675471702393E-36))) * alpha
                                          + float_type(BOOST_FLOAT128_C(1.21612504934087520075905434734158045947460467096773246215239E-34))) * alpha
                                          + float_type(BOOST_FLOAT128_C(3.76998762883139753126119821241037824830069851253295480396224E-33))) * alpha
                                          + float_type(BOOST_FLOAT128_C(1.13099628863830344684998293828608215735777107850991029729440E-31))) * alpha
                                          + float_type(BOOST_FLOAT128_C(3.27988923706982293204067897468714277771890104022419696770352E-30))) * alpha
                                          + float_type(BOOST_FLOAT128_C(9.18368986379558482800593745627556950089950023355628325088207E-29))) * alpha
                                          + float_type(BOOST_FLOAT128_C(2.47959626322479746949155352659617642905315302382639380521497E-27))) * alpha
                                          + float_type(BOOST_FLOAT128_C(6.44695028438447337900255966737803112935639344283098705091949E-26))) * alpha
                                          + float_type(BOOST_FLOAT128_C(1.61173757109611834904452725462599961406036904573072897122957E-24))) * alpha
                                          + float_type(BOOST_FLOAT128_C(3.86817017063068403772269360016918092488847584660382953555804E-23))) * alpha
                                          + float_type(BOOST_FLOAT128_C(8.89679139245057328674891109315654704307721758924206107351744E-22))) * alpha
                                          + float_type(BOOST_FLOAT128_C(1.95729410633912612308475595397946731738088422488032228717097E-20))) * alpha
                                          + float_type(BOOST_FLOAT128_C(4.11031762331216485847799061511674191805055663711439605760231E-19))) * alpha
                                          + float_type(BOOST_FLOAT128_C(8.22063524662432971695598123977873600603370758794431071426640E-18))) * alpha
                                          + float_type(BOOST_FLOAT128_C(1.56192069685862264622163643500633782667263448653185159383285E-16))) * alpha
                                          + float_type(BOOST_FLOAT128_C(2.81145725434552076319894558300988749849555291507956994126835E-15))) * alpha
                                          + float_type(BOOST_FLOAT128_C(4.77947733238738529743820749111754320727153728139716409114011E-14))) * alpha
                                          + float_type(BOOST_FLOAT128_C(7.64716373181981647590113198578807092707697416852226691068627E-13))) * alpha
                                          + float_type(BOOST_FLOAT128_C(1.14707455977297247138516979786821056670509688396295740818677E-11))) * alpha
                                          + float_type(BOOST_FLOAT128_C(1.60590438368216145993923771701549479323291461578567184216302E-10))) * alpha
                                          + float_type(BOOST_FLOAT128_C(2.08767569878680989792100903212014323125428376052986408239620E-09))) * alpha
                                          + float_type(BOOST_FLOAT128_C(2.50521083854417187750521083854417187750523408006206780016659E-08))) * alpha
                                          + float_type(BOOST_FLOAT128_C(2.75573192239858906525573192239858906525573195144226062684604E-07))) * alpha
                                          + float_type(BOOST_FLOAT128_C(2.75573192239858906525573192239858906525573191310049321957902E-06))) * alpha
                                          + float_type(BOOST_FLOAT128_C(0.00002480158730158730158730158730158730158730158730149317774)))     * alpha
                                          + float_type(BOOST_FLOAT128_C(0.00019841269841269841269841269841269841269841269841293575920)))     * alpha
                                          + float_type(BOOST_FLOAT128_C(0.00138888888888888888888888888888888888888888888888889071045)))     * alpha
                                          + float_type(BOOST_FLOAT128_C(0.00833333333333333333333333333333333333333333333333332986595)))     * alpha
                                          + float_type(BOOST_FLOAT128_C(0.04166666666666666666666666666666666666666666666666666664876)))     * alpha
                                          + float_type(BOOST_FLOAT128_C(0.16666666666666666666666666666666666666666666666666666669048)))     * alpha
                                          + float_type(BOOST_FLOAT128_C(0.50000000000000000000000000000000000000000000000000000000006)))     * alpha
                                          + float_type(BOOST_FLOAT128_C(0.99999999999999999999999999999999999999999999999999999999995)))     * alpha
                                          + float_type(1));

    // Rescale the result and return it.
    return sum * boost::cstdfloat::detail::pown(float_type(2), n);
  }
  extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_POW   (boost::cstdfloat::detail::float_internal128_t, boost::cstdfloat::detail::float_internal128_t) throw();
  extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_LOG   (boost::cstdfloat::detail::float_internal128_t) throw();
  extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_LOG10 (boost::cstdfloat::detail::float_internal128_t) throw();
  extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_SIN   (boost::cstdfloat::detail::float_internal128_t) throw();
  extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_COS   (boost::cstdfloat::detail::float_internal128_t) throw();
  extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_TAN   (boost::cstdfloat::detail::float_internal128_t) throw();
  extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_ASIN  (boost::cstdfloat::detail::float_internal128_t) throw();
  extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_ACOS  (boost::cstdfloat::detail::float_internal128_t) throw();
  extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_ATAN  (boost::cstdfloat::detail::float_internal128_t) throw();
  inline     boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_SINH  (boost::cstdfloat::detail::float_internal128_t x)
  {
    // Patch the sinhq() function for a subset of broken GCC compilers
    // like GCC 4.7, 4.8 on MinGW.
    const boost::cstdfloat::detail::float_internal128_t ex = ::BOOST_CSTDFLOAT_FLOAT128_EXP(x);
    return (ex - (boost::cstdfloat::detail::float_internal128_t(1) / ex)) / 2;
  }
  inline     boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_COSH  (boost::cstdfloat::detail::float_internal128_t x)
  {
    // Patch the coshq() function for a subset of broken GCC compilers
    // like GCC 4.7, 4.8 on MinGW.
    const boost::cstdfloat::detail::float_internal128_t ex = ::BOOST_CSTDFLOAT_FLOAT128_EXP(x);
    return (ex + (boost::cstdfloat::detail::float_internal128_t(1) / ex)) / 2;
  }
  inline     boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_TANH  (boost::cstdfloat::detail::float_internal128_t x)
  {
    // Patch the tanhq() function for a subset of broken GCC compilers
    // like GCC 4.7, 4.8 on MinGW.
    const boost::cstdfloat::detail::float_internal128_t ex_plus  = ::BOOST_CSTDFLOAT_FLOAT128_EXP(x);
    const boost::cstdfloat::detail::float_internal128_t ex_minus = (boost::cstdfloat::detail::float_internal128_t(1) / ex_plus);
    return (ex_plus - ex_minus) / (ex_plus + ex_minus);
  }
  extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_FMOD  (boost::cstdfloat::detail::float_internal128_t, boost::cstdfloat::detail::float_internal128_t) throw();
  extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_ATAN2 (boost::cstdfloat::detail::float_internal128_t, boost::cstdfloat::detail::float_internal128_t) throw();
  extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_LGAMMA(boost::cstdfloat::detail::float_internal128_t) throw();
  extern "C" boost::cstdfloat::detail::float_internal128_t BOOST_CSTDFLOAT_FLOAT128_TGAMMA(boost::cstdfloat::detail::float_internal128_t) throw();

  // Define the quadruple-precision <cmath> functions in the namespace
  // boost::cstdfloat::detail.

  namespace boost { namespace cstdfloat { namespace detail {
  inline   boost::cstdfloat::detail::float_internal128_t ldexp (boost::cstdfloat::detail::float_internal128_t x, int n)                                           { return ::BOOST_CSTDFLOAT_FLOAT128_LDEXP (x, n); }
  inline   boost::cstdfloat::detail::float_internal128_t frexp (boost::cstdfloat::detail::float_internal128_t x, int* pn)                                         { return ::BOOST_CSTDFLOAT_FLOAT128_FREXP (x, pn); }
  inline   boost::cstdfloat::detail::float_internal128_t fabs  (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_FABS  (x); }
  inline   boost::cstdfloat::detail::float_internal128_t abs   (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_FABS  (x); }
  inline   boost::cstdfloat::detail::float_internal128_t floor (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_FLOOR (x); }
  inline   boost::cstdfloat::detail::float_internal128_t ceil  (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_CEIL  (x); }
  inline   boost::cstdfloat::detail::float_internal128_t sqrt  (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_SQRT  (x); }
  inline   boost::cstdfloat::detail::float_internal128_t trunc (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_TRUNC (x); }
  inline   boost::cstdfloat::detail::float_internal128_t exp   (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_EXP   (x); }
  inline   boost::cstdfloat::detail::float_internal128_t pow   (boost::cstdfloat::detail::float_internal128_t x, boost::cstdfloat::detail::float_internal128_t a) { return ::BOOST_CSTDFLOAT_FLOAT128_POW   (x, a); }
  inline   boost::cstdfloat::detail::float_internal128_t log   (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_LOG   (x); }
  inline   boost::cstdfloat::detail::float_internal128_t log10 (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_LOG10 (x); }
  inline   boost::cstdfloat::detail::float_internal128_t sin   (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_SIN   (x); }
  inline   boost::cstdfloat::detail::float_internal128_t cos   (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_COS   (x); }
  inline   boost::cstdfloat::detail::float_internal128_t tan   (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_TAN   (x); }
  inline   boost::cstdfloat::detail::float_internal128_t asin  (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_ASIN  (x); }
  inline   boost::cstdfloat::detail::float_internal128_t acos  (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_ACOS  (x); }
  inline   boost::cstdfloat::detail::float_internal128_t atan  (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_ATAN  (x); }
  inline   boost::cstdfloat::detail::float_internal128_t sinh  (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_SINH  (x); }
  inline   boost::cstdfloat::detail::float_internal128_t cosh  (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_COSH  (x); }
  inline   boost::cstdfloat::detail::float_internal128_t tanh  (boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_TANH  (x); }
  inline   boost::cstdfloat::detail::float_internal128_t fmod  (boost::cstdfloat::detail::float_internal128_t a, boost::cstdfloat::detail::float_internal128_t b) { return ::BOOST_CSTDFLOAT_FLOAT128_FMOD  (a, b); }
  inline   boost::cstdfloat::detail::float_internal128_t atan2 (boost::cstdfloat::detail::float_internal128_t y, boost::cstdfloat::detail::float_internal128_t x) { return ::BOOST_CSTDFLOAT_FLOAT128_ATAN2 (y, x); }
  inline   boost::cstdfloat::detail::float_internal128_t lgamma(boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_LGAMMA(x); }
  inline   boost::cstdfloat::detail::float_internal128_t tgamma(boost::cstdfloat::detail::float_internal128_t x)                                                  { return ::BOOST_CSTDFLOAT_FLOAT128_TGAMMA(x); }
  } } }  // boost::cstdfloat::detail

  // Now inject the quadruple-precision <cmath> functions into the std namespace.
  namespace std
  {
    using boost::cstdfloat::detail::ldexp;
    using boost::cstdfloat::detail::frexp;
    using boost::cstdfloat::detail::fabs;
    using boost::cstdfloat::detail::abs;
    using boost::cstdfloat::detail::floor;
    using boost::cstdfloat::detail::ceil;
    using boost::cstdfloat::detail::sqrt;
    using boost::cstdfloat::detail::trunc;
    using boost::cstdfloat::detail::exp;
    using boost::cstdfloat::detail::pow;
    using boost::cstdfloat::detail::log;
    using boost::cstdfloat::detail::log10;
    using boost::cstdfloat::detail::sin;
    using boost::cstdfloat::detail::cos;
    using boost::cstdfloat::detail::tan;
    using boost::cstdfloat::detail::asin;
    using boost::cstdfloat::detail::acos;
    using boost::cstdfloat::detail::atan;
    using boost::cstdfloat::detail::sinh;
    using boost::cstdfloat::detail::cosh;
    using boost::cstdfloat::detail::tanh;
    using boost::cstdfloat::detail::fmod;
    using boost::cstdfloat::detail::atan2;
    using boost::cstdfloat::detail::lgamma;
    using boost::cstdfloat::detail::tgamma;
  } // namespace std

  #endif // Not BOOST_CSTDFLOAT_NO_LIBQUADMATH_SUPPORT (i.e., the user would like to have libquadmath support)

#endif // _BOOST_CSTDFLOAT_CMATH_2014_02_15_HPP_
