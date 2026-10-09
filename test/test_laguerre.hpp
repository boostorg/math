// Copyright John Maddock 2006.
// Copyright Paul A. Bristow 2007, 2009
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/concepts/real_concept.hpp>
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/math/special_functions/math_fwd.hpp>
#include <boost/math/special_functions/laguerre.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/test_value.hpp>
#include <boost/array.hpp>
#include "functor.hpp"

#include "math_unit_test.hpp"
#include "handle_test_result.hpp"
#include "table_type.hpp"

#ifndef SC_
#define SC_(x) static_cast<typename table_type<T>::type>(BOOST_JOIN(x, L))
#endif

template <class Real, class T>
void do_test_laguerre2(const T& data, const char* type_name, const char* test_name)
{
#if !(defined(ERROR_REPORTING_MODE) && !defined(LAGUERRE_FUNCTION_TO_TEST))
   typedef Real                   value_type;

   typedef value_type (*pg)(unsigned, value_type);
#ifdef LAGUERRE_FUNCTION_TO_TEST
   pg funcp = LAGUERRE_FUNCTION_TO_TEST;
#elif defined(BOOST_MATH_NO_DEDUCED_FUNCTION_POINTERS)
   pg funcp = boost::math::laguerre<value_type>;
#else
   pg funcp = boost::math::laguerre;
#endif

   boost::math::tools::test_result<value_type> result;

   std::cout << "Testing " << test_name << " with type " << type_name
      << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

   //
   // test laguerre against data:
   //
   result = boost::math::tools::test_hetero<Real>(
      data,
      bind_func_int1<Real>(funcp, 0, 1),
      extract_result<Real>(2));
   handle_test_result(result, data[result.worst()], result.worst(), type_name, "laguerre(n, x)", test_name);

   std::cout << std::endl;
#endif
}

template <class Real, class T>
void do_test_laguerre3(const T& data, const char* type_name, const char* test_name)
{
#if !(defined(ERROR_REPORTING_MODE) && !defined(ASSOC_LAGUERRE_FUNCTION_TO_TEST))
   typedef Real                   value_type;

   typedef value_type (*pg)(unsigned, unsigned, value_type);
#ifdef ASSOC_LAGUERRE_FUNCTION_TO_TEST
   pg funcp = ASSOC_LAGUERRE_FUNCTION_TO_TEST;
#elif defined(BOOST_MATH_NO_DEDUCED_FUNCTION_POINTERS)
   pg funcp = boost::math::laguerre<unsigned, value_type>;
#else
   pg funcp = boost::math::laguerre;
#endif

   boost::math::tools::test_result<value_type> result;

   std::cout << "Testing " << test_name << " with type " << type_name
      << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

   //
   // test laguerre against data:
   //
   result = boost::math::tools::test_hetero<Real>(
      data,
      bind_func_int2<Real>(funcp, 0, 1, 2),
      extract_result<Real>(3));
   handle_test_result(result, data[result.worst()], result.worst(), type_name, "laguerre(n, m, x)", test_name);
   std::cout << std::endl;
#endif
}

template <class T>
void test_laguerre(T, const char* name)
{
   //
   // The actual test data is rather verbose, so it's in a separate file
   //
   // The contents are as follows, each row of data contains
   // three items, input value a, input value b and erf(a, b):
   //
#  include "laguerre2.ipp"

   do_test_laguerre2<T>(laguerre2, name, "Laguerre Polynomials");

#  include "laguerre3.ipp"

   do_test_laguerre3<T>(laguerre3, name, "Associated Laguerre Polynomials");
}

template <class T>
void test_spots(T, const char* t)
{
   std::cout << "Testing basic sanity checks for type " << t << std::endl;
   //
   // basic sanity checks, tolerance is 100 epsilon:
   //
   T tolerance = boost::math::tools::epsilon<T>() * 100;
   BOOST_CHECK_CLOSE_FRACTION(::boost::math::laguerre(1, static_cast<T>(0.5L)), static_cast<T>(0.5L), tolerance);
   BOOST_CHECK_CLOSE_FRACTION(::boost::math::laguerre(4, static_cast<T>(0.5L)), static_cast<T>(-0.3307291666666666666666666666666666666667L), tolerance);
   BOOST_CHECK_CLOSE_FRACTION(::boost::math::laguerre(7, static_cast<T>(0.5L)), static_cast<T>(-0.5183392237103174603174603174603174603175L), tolerance);
   BOOST_CHECK_CLOSE_FRACTION(::boost::math::laguerre(20, static_cast<T>(0.5L)), static_cast<T>(0.3120174870800154148915399248893113634676L), tolerance);
   BOOST_CHECK_CLOSE_FRACTION(::boost::math::laguerre(50, static_cast<T>(0.5L)), static_cast<T>(-0.3181388060269979064951118308575628226834L), tolerance);

   BOOST_CHECK_CLOSE_FRACTION(::boost::math::laguerre(1, static_cast<T>(-0.5L)), static_cast<T>(1.5L), tolerance);
   BOOST_CHECK_CLOSE_FRACTION(::boost::math::laguerre(4, static_cast<T>(-0.5L)), static_cast<T>(3.835937500000000000000000000000000000000L), tolerance);
   BOOST_CHECK_CLOSE_FRACTION(::boost::math::laguerre(7, static_cast<T>(-0.5L)), static_cast<T>(7.950934709821428571428571428571428571429L), tolerance);
   BOOST_CHECK_CLOSE_FRACTION(::boost::math::laguerre(20, static_cast<T>(-0.5L)), static_cast<T>(76.12915699869631476833699787070874048223L), tolerance);
   BOOST_CHECK_CLOSE_FRACTION(::boost::math::laguerre(50, static_cast<T>(-0.5L)), static_cast<T>(2307.428631277506570629232863491518399720L), tolerance);

   BOOST_CHECK_CLOSE_FRACTION(::boost::math::laguerre(1, static_cast<T>(4.5L)), static_cast<T>(-3.500000000000000000000000000000000000000L), tolerance);
   BOOST_CHECK_CLOSE_FRACTION(::boost::math::laguerre(4, static_cast<T>(4.5L)), static_cast<T>(0.08593750000000000000000000000000000000000L), tolerance);
   BOOST_CHECK_CLOSE_FRACTION(::boost::math::laguerre(7, static_cast<T>(4.5L)), static_cast<T>(-1.036928013392857142857142857142857142857L), tolerance);
   BOOST_CHECK_CLOSE_FRACTION(::boost::math::laguerre(20, static_cast<T>(4.5L)), static_cast<T>(1.437239150257817378525582974722170737587L), tolerance);
   BOOST_CHECK_CLOSE_FRACTION(::boost::math::laguerre(50, static_cast<T>(4.5L)), static_cast<T>(-0.7795068145562651416494321484050019245248L), tolerance);

   BOOST_CHECK_CLOSE_FRACTION(::boost::math::laguerre(4, 5, static_cast<T>(0.5L)), static_cast<T>(88.31510416666666666666666666666666666667L), tolerance);
   BOOST_CHECK_CLOSE_FRACTION(::boost::math::laguerre(10, 0, static_cast<T>(2.5L)), static_cast<T>(-0.8802526766660982969576719576719576719577L), tolerance);
   BOOST_CHECK_CLOSE_FRACTION(::boost::math::laguerre(10, 1, static_cast<T>(4.5L)), static_cast<T>(1.564311458042689732142857142857142857143L), tolerance);
   BOOST_CHECK_CLOSE_FRACTION(::boost::math::laguerre(10, 6, static_cast<T>(8.5L)), static_cast<T>(20.51596541066649098875661375661375661376L), tolerance);
   BOOST_CHECK_CLOSE_FRACTION(::boost::math::laguerre(10, 12, static_cast<T>(12.5L)), static_cast<T>(-199.5560968456234671241181657848324514991L), tolerance);
   BOOST_CHECK_CLOSE_FRACTION(::boost::math::laguerre(50, 40, static_cast<T>(12.5L)), static_cast<T>(-4.996769495006119488583146995907246595400e16L), tolerance);

   BOOST_CHECK_EQUAL(::boost::math::laguerre(0, T(40)), T(1));
   BOOST_CHECK_EQUAL(::boost::math::laguerre(0, T(400)), T(1));
}

template <class T>
void test_zeros(T tol, const char* t)
{
   std::cout << "Testing zeros of Laguerre polynomials for type " << t << std::endl;

   unsigned N = 5;
   std::vector<T> roots = boost::math::laguerre_zeros<T>(N);
   // Values calculated via Mathematica to 115 digits of precision.
   BOOST_CHECK_CLOSE_FRACTION(roots[0], BOOST_MATH_TEST_VALUE(T, 0.263560319718140910203061943360833334689007569905516927862615028383113143804688485388801428636768142471684891181097), tol);
   BOOST_CHECK_CLOSE_FRACTION(roots[1], BOOST_MATH_TEST_VALUE(T, 1.413403059106516792218407980187557749539096004380880280554408297900742377592354294632071875549162504722547195499148), tol);
   BOOST_CHECK_CLOSE_FRACTION(roots[2], BOOST_MATH_TEST_VALUE(T, 3.596425771040722081223186588782971665671150940710574538253074419444595035340244283770723676354534463586988774883226), tol);
   BOOST_CHECK_CLOSE_FRACTION(roots[3], BOOST_MATH_TEST_VALUE(T, 7.085810005858837556922124181108086000385935670784940882578921655737736629317011861316381315331093202242780296517779), tol);
   BOOST_CHECK_CLOSE_FRACTION(roots[4], BOOST_MATH_TEST_VALUE(T, 12.64080084427578265943321930656055124971480981421808737075098059853381281394570107489202170412844168697599884191874), tol);

   N = 2;
   roots = boost::math::laguerre_zeros<T>(N);
   BOOST_CHECK_CLOSE_FRACTION(roots[0], BOOST_MATH_TEST_VALUE(T, 0.585786437626904951198311275790301921430328124623051926823320262009267521537892961149612465672358427264986153769087), tol);
   BOOST_CHECK_CLOSE_FRACTION(roots[1], BOOST_MATH_TEST_VALUE(T, 3.414213562373095048801688724209698078569671875376948073176679737990732478462107038850387534327641572735013846230912), tol);
}

template <class T>
void test_zeros_large(T tol, const char* t)
{
   std::cout << "Testing large zeros of Laguerre polynomials for type " << t << std::endl;

   unsigned N = 30;
   std::vector<T> roots = boost::math::laguerre_zeros<T>(N);
   // Values calculated via Mathematica to 115 digits of precision.
   BOOST_CHECK_CLOSE_FRACTION(roots[0], BOOST_MATH_TEST_VALUE(T, 0.047407180540804851461661096626938515574196972256004488033124253892354077063440713062711305445935060296219853729199), tol);
   BOOST_CHECK_CLOSE_FRACTION(roots[1], BOOST_MATH_TEST_VALUE(T, 0.249923916753160223993729741481384213586009854858782922019446030331750308159544972978742717449917581786215730371944), tol);
   BOOST_CHECK_CLOSE_FRACTION(roots[2], BOOST_MATH_TEST_VALUE(T, 0.614833454392768284612976751633664856300777548726368191642899058361772180963140823942717398338905410794105984943534), tol);
   BOOST_CHECK_CLOSE_FRACTION(roots[3], BOOST_MATH_TEST_VALUE(T, 1.143195825666100798284410186232713050561122938419877909989949066107255568412352942668557861479865752638218002014316), tol);
   BOOST_CHECK_CLOSE_FRACTION(roots[4], BOOST_MATH_TEST_VALUE(T, 1.836454554622572291486188518379157894668816037615102323082905276230690030155087909682240820321762881725087181348823), tol);
}

template <class T>
void test_zeros_accuracy(T tol, const char* t)
{
   std::cout << "Testing zeros accuracy for type " << t << std::endl;

   // The machine epsilon for double precisionis about 1e-16, and the polynomial is very
   // sensitive to this error for large roots and large n. Thus, the accuracy is quite
   // poor for large n simply due to precision.
   for (unsigned n = 6; n < 8; ++n)
   {
      std::vector<T> zeros = boost::math::laguerre_zeros<T>(n);
      BOOST_CHECK(zeros.size() == n);
      BOOST_CHECK(zeros[0] > 0);

      for (unsigned k=0; k < zeros.size(); ++k)
      {
         BOOST_CHECK_SMALL(boost::math::laguerre(n, zeros[k]), 5000 * tol);
      }
   }
   return;
}

template <class T>
void test_zero_special_case()
{
   std::cout << "Testing special cases for zeros" << std::endl;

   BOOST_CHECK_THROW(boost::math::laguerre_zeros<T>(0), std::domain_error);
   BOOST_CHECK_THROW(boost::math::laguerre_zeros<T>(std::numeric_limits<unsigned>::quiet_NaN()), std::domain_error);

   using ignore_policy = boost::math::policies::policy<boost::math::policies::domain_error<boost::math::policies::ignore_error>,
                                                       boost::math::policies::evaluation_error<boost::math::policies::ignore_error> >;
   CHECK_NAN((boost::math::laguerre_zeros<T, ignore_policy>(0)[0]));
   CHECK_NAN((boost::math::laguerre_zeros<T, ignore_policy>(std::numeric_limits<unsigned>::quiet_NaN())[0]));

   std::vector<T> zeros = boost::math::laguerre_zeros<T>(1);
   BOOST_CHECK(zeros.size() == 1);
   BOOST_CHECK(zeros[0] == static_cast<T>(1));

   unsigned N = 500;
   BOOST_CHECK_THROW(boost::math::laguerre_zeros<T>(N), boost::math::evaluation_error);

   std::vector<T> roots = boost::math::laguerre_zeros<T, ignore_policy>(N);
   BOOST_CHECK(roots.size() == N);
   CHECK_NAN(roots[0]);

}
