// Copyright John Maddock 2006.
// Copyright Paul A. Bristow 2007, 2009
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_MATH_OVERFLOW_ERROR_POLICY ignore_error

#include <boost/math/concepts/real_concept.hpp>
#include <boost/math/special_functions/math_fwd.hpp>
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/math/tools/stats.hpp>
#include <boost/math/tools/test.hpp>
#include <boost/math/tools/big_constant.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/array.hpp>
#include "functor.hpp"

#include "handle_test_result.hpp"
#include "table_type.hpp"

#include <boost/math/special_functions/hypergeometric_1F1.hpp>
#include <boost/math/quadrature/exp_sinh.hpp>

#ifdef BOOST_MSVC
#pragma warning(disable:4127)
#endif

template <class Real, class T>
void do_test_1F1(const T& data, const char* type_name, const char* test_name)
{
   typedef Real                   value_type;

   typedef value_type(*pg)(value_type, value_type, value_type);
#if defined(BOOST_MATH_NO_DEDUCED_FUNCTION_POINTERS)
   pg funcp = boost::math::hypergeometric_0F1<value_type, value_type>;
#else
   pg funcp = boost::math::hypergeometric_1F1;
#endif

   boost::math::tools::test_result<value_type> result;

   std::cout << "Testing " << test_name << " with type " << type_name
      << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

   //
   // test hypergeometric_2F0 against data:
   //
   result = boost::math::tools::test_hetero<Real>(
      data,
      bind_func<Real>(funcp, 0, 1, 2),
      extract_result<Real>(3));
   handle_test_result(result, data[result.worst()], result.worst(), type_name, "hypergeometric_1F1", test_name);
   std::cout << std::endl;
}

#ifndef SC_
#define SC_(x) BOOST_MATH_BIG_CONSTANT(T, 1000000, x)
#endif

template <class T>
void test_spots1(T, const char* type_name)
{
#include "hypergeometric_1F1.ipp"

   do_test_1F1<T>(hypergeometric_1F1, type_name, "Integer a values");

#include "hypergeometric_1F1_small_random.ipp"

   do_test_1F1<T>(hypergeometric_1F1_small_random, type_name, "Small random values");
}

template <class T>
void test_spots2(T, const char* type_name)
{
#include "hypergeometric_1F1_big.ipp"

   do_test_1F1<T>(hypergeometric_1F1_big, type_name, "Large random values");
}

template <class T>
void test_spots3(T, const char* type_name)
{
#include "hypergeometric_1F1_big_double_limited.ipp"

   do_test_1F1<T>(hypergeometric_1F1_big_double_limited, type_name, "Large random values - double limited precision");
}

template <class T>
void test_spots4(T, const char* type_name)
{
#include "hypergeometric_1F1_big_unsolved.ipp"

   do_test_1F1<T>(hypergeometric_1F1_big, type_name, "Large random values - unsolved domains");
}

template <class T>
void test_spots5(T, const char* type_name)
{
   std::cout << "Testing special cases for type " << type_name << std::endl;
   BOOST_MATH_STD_USING
   //
   // Special cases:
   //
   using boost::math::hypergeometric_1F1;
   T tol = boost::math::tools::epsilon<T>() * 200;
   if (std::numeric_limits<T>::digits > std::numeric_limits<double>::digits)
      tol *= 2;
   if (boost::is_class<T>::value)
      tol *= 4;
   // b = 2a
   T computed = hypergeometric_1F1(T(-12.25), T(2 * -12.25), T(6.75));
   T expected = boost::lexical_cast<T>("22.995348157760091167706081204212893687052775606591209203948675272473773725021024450870565197330528784707135828761");
   BOOST_CHECK_CLOSE(computed, expected, tol);
   computed = hypergeometric_1F1(T(12.25), T(2 * 12.25), T(6.75));
   expected = boost::lexical_cast<T>("36.47281964229300610642392880149257389834650024065756742702265701321933782423217084029882132197130099355867287657");
   BOOST_CHECK_CLOSE(computed, expected, tol);
   computed = hypergeometric_1F1(T(-11), T(-12), T(6.75));
   expected = boost::lexical_cast<T>("376.3166426246459656334542608880377435064935064935064935064935064935064935064935064935064935064935064935064935064");
   BOOST_CHECK_CLOSE(computed, expected, tol);
   computed = hypergeometric_1F1(T(-2), T(-12), T(6.75));
   expected = boost::lexical_cast<T>("2.470170454545454545454545454545454545454545454545454545454545454545454545454545454545454545454545454545454545");
   BOOST_CHECK_CLOSE(computed, expected, tol);
   computed = hypergeometric_1F1(T(-224), T(-1205), T(6.75));
   expected = boost::lexical_cast<T>("3.497033449657595724636676193024114597507981035316405619832857546161530808157860391434240068189887198094611519953");
   BOOST_CHECK_CLOSE(computed, expected, tol);
   computed = hypergeometric_1F1(T(0.5), T(-1205.5), T(-6.75));
   expected = boost::lexical_cast<T>("1.00281149043026925155096279505879868076290060374397866773878698584557482321961231721407215665017657501846692575");
   BOOST_CHECK_CLOSE(computed, expected, tol);
   computed = hypergeometric_1F1(T(-0.5), T(-1205.5), T(-6.75));
   expected = boost::lexical_cast<T>("0.99719639844965644594352920596780535220516138060108955206195178371227403775248888108818326220977962797312690");
   BOOST_CHECK_CLOSE(computed, expected, tol);
   computed = hypergeometric_1F1(T(-12), T(16.25), T(1043.75));
   expected = boost::lexical_cast<T>("1.26527673505477678311707565502355407505496430400394171269315320194708537626079491650410923064978320042481912e20");
   BOOST_CHECK_CLOSE(computed, expected, tol * 3);
   
   computed = hypergeometric_1F1(T(3.5), T(3.5), T(36.25));
   expected = exp(T(36.25));
   BOOST_CHECK_CLOSE(computed, expected, tol);
   computed = hypergeometric_1F1(T(-3.5), T(-3.5), T(36.25));
   expected = exp(T(36.25));
   BOOST_CHECK_CLOSE(computed, expected, tol);
   computed = hypergeometric_1F1(T(1), T(2), T(36.25));
   expected = boost::math::expm1(T(36.25)) / T(36.25);
   BOOST_CHECK_CLOSE(computed, expected, tol);
   computed = hypergeometric_1F1(T(10.25), T(9.25), T(36.25));
   expected = exp(T(36.25)) * (T(9.25) + T(36.25)) / T(9.25);
   BOOST_CHECK_CLOSE(computed, expected, tol);
   computed = hypergeometric_1F1(T(-10.25), T(-11.25), T(36.25));
   expected = exp(T(36.25)) * (T(-11.25) + T(36.25)) / T(-11.25);
   BOOST_CHECK_CLOSE(computed, expected, tol);
   computed = hypergeometric_1F1(T(-10.25), T(-11.25), T(-36.25));
   expected = exp(T(-36.25)) * (T(-11.25) + T(-36.25)) / T(-11.25);
   BOOST_CHECK_CLOSE(computed, expected, tol);
}

template <class T>
void test_spots6(T, const char* type_name)
{
#ifndef SC_
#  define SC_(x) static_cast<T>(BOOST_JOIN(x, L))
#endif
   static const boost::array<boost::array<T, 4>, 35> hypergeometric_1F1_bugs = { {
  { { SC_(17955.561660766602), SC_(9.6968994205831605e-09), SC_(-82.406154185533524), SC_(6.98056008378736714088730927132364938220428678e-11) }},
  { { SC_(17955.561660766602), SC_(-9.6968994205831605e-09), SC_(-82.406154185533524), SC_(-6.98055306629610746072607353939306734740549551e-11) }},
  { { SC_(-17955.561660766602), SC_(-9.6968994205831605e-09), SC_(82.406154185533524), SC_(-42897094853118832762870100.8669248353530950866) }} ,
  { { SC_(17955.561660766602), SC_(17956.061660766602), SC_(82.406154185533524), SC_(613117565438499794408370861624072730.553215432) }},
  { { SC_(2.9127331452327709e-07), SC_(-0.99999970872668542), SC_(0.15018942760070786), SC_(0.987526018990506843793601092932108059727149508) }},
  { { SC_(-2.9127331452327709e-07), SC_(-1.0000002912733146), SC_(0.15018942760070786), SC_(0.987526120661366412484942089372497015837368389) }},
  { { SC_(6.7191087900739423e-13), SC_(-0.99999999999932809), SC_(0.0011913633891253994), SC_(0.999999289758605006762757201699750974296453229) }},
  { { SC_(6.7191087900739423e-13), SC_(-0.99999999999932809), SC_(-0.0011913633891253994), SC_(0.999999290885918468326416221021126912154021802) }},
  { { SC_(-6.7191087900739423e-13), SC_(-1.0000000000006719), SC_(0.0011913633891253994), SC_(0.999999289758606609651292394510404091049823243) }},
  { { SC_(-6.7191087900739423e-13), SC_(-1.0000000000006719), SC_(-0.0011913633891253994), SC_(0.999999290885916869252591036674587894145399498) }},
  { { SC_(1.2860067365774887e-17), SC_(6.2442285664031425e-16), SC_(-2539.60133934021), SC_(0.979404874070484696999110600576068012417904384) }},
  { { SC_(1.2860067365774887e-17), SC_(-6.2442285664031425e-16), SC_(-2539.60133934021), SC_(1.0205951259295150865252112924093487321207727) }},
  { { SC_(-1.2860067365774887e-17), SC_(6.2442285664031425e-16), SC_(-2539.60133934021), SC_(1.02059512592951530745923325071510441026202975) }},
  { { SC_(-1.2860067365774887e-17), SC_(-6.2442285664031425e-16), SC_(-2539.60133934021), SC_(0.979404874070484909016444856299500644331897735) }},
  { { SC_(1.2860067365774887e-17), SC_(1), SC_(-2539.60133934021), SC_(0.999999999999999891757095137551552220860540801) }},
  { { SC_(-1.2860067365774887e-17), SC_(1), SC_(-2539.60133934021), SC_(1.00000000000000010824290486244845922375479178) }},
  { { SC_(1.2860067365774887e-17), SC_(0.5), SC_(-2539.60133934021), SC_(0.999999999999999873931788919689096760455570214) }},
  { { SC_(-1.2860067365774887e-17), SC_(0.5), SC_(-2539.60133934021), SC_(1.0000000000000001260682110803109183167444166) }},
  { { SC_(1.2860067365774887e-17), SC_(-0.5), SC_(-2539.60133934021), SC_(0.999999999999999899656990458526368219886894767) }},
  { { SC_(-1.2860067365774887e-17), SC_(-0.5), SC_(-2539.60133934021), SC_(1.00000000000000010034300954147364037131355735) }},
  { { SC_(1.9561377367172441e-13), SC_(-0.99999999999980438), SC_(0.53720525559037924), SC_(0.791950585963666119273677451162365759080483409) }},
  { { SC_(1.9561377367172441e-13), SC_(-0.99999999999980438), SC_(-0.53720525559037924), SC_(0.898314630992769591673208399706587643905527327) }},
  { { SC_(-1.9561377367172441e-13), SC_(-1.0000000000001956), SC_(0.53720525559037924), SC_(0.791950585964025761367113514279915403442035074) }},
  { { SC_(-1.9561377367172441e-13), SC_(-1.0000000000001956), SC_(-0.53720525559037924), SC_(0.898314630992646771749564140770704893561753597) }},
  { { SC_(5.1851756946064858e-12), SC_(-0.99999999999481481), SC_(-774.06985878944397), SC_(1.91306610467163858324476828831735612399803649e-06) }},
  { { SC_(-5.1851756946064858e-12), SC_(-1.0000000000051852), SC_(-774.06985878944397), SC_(1.91306610479516297551035931150910859922270467e-06) }},

  {{ SC_(4.782769898853794e-15), SC_(1.0000000000000049), SC_(43.289540141820908), SC_(715.678254892476818206948251991084031658534788) }},
  { { SC_(-4.782769898853794e-15), SC_(0.99999999999999523), SC_(43.289540141820908), SC_(-713.67825489247727251051792450091274703212426) }},
  { { SC_(4.782769898853794e-15), SC_(0.50000000000000477), SC_(43.289540141820908), SC_(8235.578376364917373771471380274179857713986) }},
  { { SC_(-4.782769898853794e-15), SC_(0.49999999999999523), SC_(43.289540141820908), SC_(-8233.57837636502669085205930058992320862281194) }},
  { { SC_(4.782769898853794e-15), SC_(-0.49999999999999523), SC_(43.289540141820908), SC_(-696269.800378137841948029488304613132151506346) }},
  { { SC_(-4.782769898853794e-15), SC_(-0.50000000000000477), SC_(43.289540141820908), SC_(696271.8003781336298001417038674968935893361) }},
  { { SC_(8.1104991963343309e-05), SC_(-0.99991889500803666), SC_(-289.12455415725708), SC_(7.89625448009377635153307897651433007437615965e-124) }},
  { { SC_(-8.1104991963343309e-05), SC_(-1.0000811049919633), SC_(-289.12455415725708), SC_(7.8949781467741574268884621364833028722017032e-124) }},

   {{ SC_(-1.98018241448205767), SC_(1.98450573845762079), SC_(54.4977916804564302), SC_(2972026581564772.790187123046255523239732028) }},
  } };
  //#undef SC_

   do_test_1F1<T>(hypergeometric_1F1_bugs, type_name, "Bug cases");
}

template <class T>
void test_spots(T z, const char* type_name)
{
   test_spots1(z, type_name);
   test_spots2(z, type_name);
   if (std::numeric_limits<T>::digits10 < 20)
      test_spots3(z, type_name);
#ifdef TEST_UNSOLVED
   test_spots4(z, type_name);
#endif
   test_spots5(z, type_name);
   test_spots6(z, type_name);
}


// Tests the Mellin transform formula given here: https://dlmf.nist.gov/13.10, Equation 13.10.10
template <class Real>
void test_hypergeometric_mellin_transform()
{
    using boost::math::hypergeometric_1F1;
    using boost::math::quadrature::exp_sinh;
    using boost::math::tgamma;
    using std::pow;

    // Constraint: 0 < lambda < a.
    Real lambda = 0.5;
    Real a = 1;
    Real b = 3;
    auto f = [&](Real t)->Real { return pow(t, lambda - 1)*hypergeometric_1F1(a, b, -t); };

    auto integrator = exp_sinh<double>();
    Real computed = integrator.integrate(f, boost::math::tools::epsilon<Real>());
    Real expected = tgamma(b)*tgamma(lambda)*tgamma(a-lambda)/(tgamma(a)*tgamma(b-lambda));

    Real tol = boost::math::tools::epsilon<Real>() * 5;
    BOOST_CHECK_CLOSE_FRACTION(computed, expected, tol);
}


// Tests the Laplace transform formula given here: https://dlmf.nist.gov/13.10, Equation 13.10.4
template <class Real>
void test_hypergeometric_laplace_transform()
{
    using boost::math::hypergeometric_1F1;
    using boost::math::quadrature::exp_sinh;
    using boost::math::tgamma;
    using std::pow;
    using std::exp;

    // Set a = 1 blows up for some reason . . .
    Real a = -1;
    Real b = 3;
    Real z = 1.5;
    auto f = [&](Real t)->Real { return exp(-z*t)*pow(t, b - 1)*hypergeometric_1F1(a, b, t); };

    auto integrator = exp_sinh<double>();
    Real computed = integrator.integrate(f, boost::math::tools::epsilon<Real>());
    Real expected = tgamma(b)/(pow(z,b)*pow(1-1/z, a));

    Real tol = boost::math::tools::epsilon<Real>() * 200;
    BOOST_CHECK_CLOSE(computed, expected, tol);
}
