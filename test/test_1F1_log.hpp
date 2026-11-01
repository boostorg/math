// Copyright John Maddock 2006.
// Copyright Paul A. Bristow 2007, 2009
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_MATH_OVERFLOW_ERROR_POLICY ignore_error

#include <boost/math/concepts/real_concept.hpp>
#include <boost/math/special_functions/math_fwd.hpp>
#include <boost/math/tools/big_constant.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/relative_difference.hpp>
#include "math_unit_test.hpp"

#include <boost/math/special_functions/hypergeometric_1F1.hpp>
#include <boost/math/quadrature/exp_sinh.hpp>

#ifdef _MSC_VER
#pragma warning(disable:4127)
#endif

// The largest relative error, in epsilon, allowed for each type and group of test data:
std::uintmax_t max_error(const std::string& type_name, const std::string& test_name);

// Checks every row {a, b, z, f(a, b, z)} of the data, and prints each row that fails:
template <class Real, class T, class F>
void test_1F1_data(const T& data, const char* type_name, const char* function_name, const char* test_name, F f)
{
   std::cout << "Testing " << test_name << " with type " << type_name << "\n";
   const std::uintmax_t bound = max_error(type_name, test_name);
   for (unsigned i = 0; i < data.size(); ++i)
   {
      const Real a = Real(data[i][0]), b = Real(data[i][1]), z = Real(data[i][2]), expected = Real(data[i][3]);
      Real computed = 0, err = 0;
      try
      {
         computed = f(a, b, z);
         err = boost::math::relative_difference(computed, expected) / boost::math::tools::epsilon<Real>();
      }
      catch (const std::exception& e)
      {
         std::cerr << e.what() << "\n";
         err = boost::math::tools::max_value<Real>();
      }
      if (!(err <= bound))
      {
         std::cerr << std::setprecision(std::numeric_limits<Real>::max_digits10) << "\033[0;31mError in " << test_name << " at row " << i << ":\033[0m\n  "
                   << function_name << "(" << a << ", " << b << ", " << z << ")\n  Expected: " << expected << "\n  Computed: " << computed
                   << std::setprecision(4) << "\n  Relative error: " << err << " epsilon, which exceeds " << bound << "\n";
         ++boost::math::test::detail::global_error_count();
      }
   }
}

template <class Real, class T>
void do_test_1F1(const T& data, const char* type_name, const char* test_name)
{
   test_1F1_data<Real>(data, type_name, "log_hypergeometric_1F1", test_name, [](Real a, Real b, Real z) { return boost::math::log_hypergeometric_1F1(a, b, z); });
}

#ifndef SC_
#define SC_(x) BOOST_MATH_BIG_CONSTANT(T, 1000000, x)
#endif

template <class T>
void test_spots1(T, const char* type_name)
{
#include "hypergeometric_1f1_log_large.ipp"

   do_test_1F1<T>(hypergeometric_1f1_log_large, type_name, "Large random values - log");
}

template <class T>
void test_spots2(T, const char* type_name)
{
#include "hypergeometric_1f1_log_large_unsolved.ipp"

   do_test_1F1<T>(hypergeometric_1f1_log_large_unsolved, type_name, "Large random values - log - unsolved");
}

template <class T>
void test_sign(T, const char* /*type_name*/)
{
#include "hypergeometric_1F1_small_random.ipp"

   for (unsigned i = 0; i < hypergeometric_1F1_small_random.size(); ++i)
   {
      int s;
      boost::math::log_hypergeometric_1F1(hypergeometric_1F1_small_random[i][0], hypergeometric_1F1_small_random[i][1], hypergeometric_1F1_small_random[i][2], &s);
      CHECK_EQUAL(s, boost::math::sign(hypergeometric_1F1_small_random[i][3]));
   }
}

template <class T>
void test_spots_bugs(T, const char* type_name)
{
   static const std::array<std::array<T, 4>, 7> hypergeometric_1F1_bugs = { {
    // Found while investigating https://github.com/boostorg/math/issues/1034
    {{ 21156.0f, 21156.0f, 11322, SC_(11322.0)}},
    {{ 21156.0f, 21155.0f, 11322, SC_(11322.428655862323560632951631114666466652986288119296800531328684)}},
    {{ 21156.0f, 21154.0f, 11322, SC_(11322.857338938931770780542471014439235046236959098048505808665509)}},
    {{ 21156.0f, 21154.5f, 11322, SC_(11322.642993998700652342915766513423502332460377941025163971463001)}},
    {{ 21156.0f, 21155.0f - 1.0f / 128, 11322, SC_(11322.43200484338133063401686149227756707)}},
    {{ 21156.0f, 21154.5f - 1.0f / 128, 11322, SC_(11322.646343086066466097278446364687150256282052775170519455915995)}},
    {{ 21156.0f, 21154.0f + 1.0f / 128, 11322, SC_(11322.85398974691465958225700429672975704)}},
   } };

   do_test_1F1<T>(hypergeometric_1F1_bugs, type_name, "Large random values - log - bug cases");
}

template <class T>
void test_spots(T z, const char* type_name)
{
   test_spots1(z, type_name);
   test_spots_bugs(z, type_name);
   test_sign(z, type_name);
#ifdef TEST_UNSOLVED
   test_spots2(z, type_name);
#endif
}

