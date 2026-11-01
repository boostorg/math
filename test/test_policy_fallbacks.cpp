//  (C) Copyright Nick Thompson 2026.
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//
// Where an internal evaluation fails and Boost falls back to another method,
// every error handling policy must give the default policy's result, and the
// error that was recovered from must not be reported.
//
#include "math_unit_test.hpp"
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/hypergeometric_1F1.hpp>
#include <boost/math/special_functions/owens_t.hpp>

static int user_error_calls = 0;

namespace boost { namespace math { namespace policies {
template <class T>
T user_overflow_error(const char*, const char*, const T& val)
{ ++user_error_calls; return val; }
template <class T>
T user_evaluation_error(const char*, const char*, const T& val)
{ ++user_error_calls; return val; }
}}}

using namespace boost::math::policies;
using errno_policy = policy<overflow_error<errno_on_error>, evaluation_error<errno_on_error>>;
using ignore_policy = policy<overflow_error<ignore_error>, evaluation_error<ignore_error>>;
using user_policy = policy<overflow_error<user_error>, evaluation_error<user_error>>;

template <class F>
void check_policies(F f)
{
   auto expected = f(policy<>());
   // errno itself isn't checked: glibc's and MSVC's exp, pow, ... set ERANGE on intermediate underflow.
   // The user_error count below is what shows that nothing was reported.
   CHECK_EQUAL(expected, f(errno_policy()));
   CHECK_EQUAL(expected, f(ignore_policy()));
   user_error_calls = 0;
   CHECK_EQUAL(expected, f(user_policy()));
   CHECK_EQUAL(0, user_error_calls);
}

int main()
{
   using boost::math::hypergeometric_1F1;
   // tgamma(200.5) overflows, so tgamma(-200.5) underflows:
   check_policies([](auto pol) { return boost::math::tgamma(-200.5, pol); });
   // tgamma(a + 1) overflows in the small x series:
   check_policies([](auto pol) { return boost::math::gamma_p(200.5, 1e-10, pol); });
   check_policies([](auto pol) { return boost::math::gamma_q(200.5, 1e-10, pol); });
   // beta(a, b) overflows in the initial estimates:
   check_policies([](auto pol) { return boost::math::ibeta_inv(1e-310, 20.0, 1e-100, pol); });
   check_policies([](auto pol) { return boost::math::ibeta_inv(20.0, 1e-310, 1e-100, pol); });
   check_policies([](auto pol) { return boost::math::ibetac_inv(1e-314, 1e-319, 1e-95, pol); });
   // The large z asymptotic series fails:
   check_policies([](auto pol) { return hypergeometric_1F1(-100.5, 1.5, 1000.0, pol); });
   check_policies([](auto pol) { return hypergeometric_1F1(100.5, 0.25, -1000.0, pol); });
   // tgamma(b) overflows, or the Bessel functions fail, in Tricomi's approximation:
   check_policies([](auto pol) { return hypergeometric_1F1(-250.5, 224.5, 38.0, pol); });
   check_policies([](auto pol) { return hypergeometric_1F1(-204.5, 715.5, 17.0, pol); });
   check_policies([](auto pol) { return hypergeometric_1F1(600.0, -500.5, -5.0, pol); });
   // The T4 series for Owen's T fails to converge:
   using boost::multiprecision::cpp_bin_float_50;
   check_policies([](auto pol) { return boost::math::owens_t(cpp_bin_float_50(7.25), cpp_bin_float_50(1.03125), pol); });
   return boost::math::test::report_errors();
}
