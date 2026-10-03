// Copyright 2026 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// Built with BOOST_MATH_NAMESPACE_COMPONENTS=my_lib,math (https://github.com/boostorg/math/issues/769):
// every entity must be reachable as both BOOST_MATH_NAMESPACE::name and my_lib::math::name,
// and the two spellings must name the same entity.

#include <boost/math/special_functions.hpp>
#include <boost/math/distributions.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/policies/policy.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <cmath>
#include <iostream>
#include <limits>
#include <type_traits>

#if (defined(_MSVC_LANG) && _MSVC_LANG >= 201703L) || __cplusplus >= 201703L
#  include <boost/math/ccmath/ccmath.hpp>
#  include <boost/math/statistics/univariate_statistics.hpp>
#  include <vector>
#  define BOOST_MATH_QUALIFIED_NAME_TEST_CXX17
#endif

// Types are identical through either spelling, with or without the leading ::
static_assert(std::is_same<BOOST_MATH_NAMESPACE::normal, my_lib::math::normal>::value, "normal");
static_assert(std::is_same<::BOOST_MATH_NAMESPACE::students_t_distribution<float>, ::my_lib::math::students_t_distribution<float>>::value, "students_t");
static_assert(std::is_same<BOOST_MATH_NAMESPACE::policies::policy<>, my_lib::math::policies::policy<>>::value, "policy");
static_assert(std::is_same<BOOST_MATH_NAMESPACE::tools::eps_tolerance<double>, my_lib::math::tools::eps_tolerance<double>>::value, "eps_tolerance");

// A namespace alias works with the macro
namespace bm = BOOST_MATH_NAMESPACE;
static_assert(std::is_same<bm::normal, my_lib::math::normal>::value, "alias");

// The unqualified macro still resolves from inside a nested user namespace
namespace user
{
namespace nested
{

inline double gamma_via_macro(double x)
{
    return BOOST_MATH_NAMESPACE::tgamma(x);
}

inline double gamma_via_namespace(double x)
{
    return my_lib::math::tgamma(x);
}

} // namespace nested
} // namespace user

namespace
{

int failures {0};

// Both spellings must give identical results that match the reference value
void check_same(double via_macro, double via_namespace, double expected, const char* what)
{
    const auto tolerance {1e-13 * (std::fabs(expected) > 1 ? std::fabs(expected) : 1.0)};
    if (!(via_macro == via_namespace) || !(std::fabs(via_macro - expected) <= tolerance))
    {
        std::cerr << "FAIL " << what << ": BOOST_MATH_NAMESPACE gave " << via_macro
                  << ", my_lib::math gave " << via_namespace << ", expected " << expected << '\n';
        ++failures;
    }
}

void check_true(bool condition, const char* what)
{
    if (!condition)
    {
        std::cerr << "FAIL " << what << '\n';
        ++failures;
    }
}

using unary_function = double (*)(double);

// The same function, not two functions that happen to agree
constexpr unary_function tgamma_via_macro {&BOOST_MATH_NAMESPACE::tgamma};
constexpr unary_function tgamma_via_namespace {&my_lib::math::tgamma};
static_assert(tgamma_via_macro == tgamma_via_namespace, "tgamma is one function");

constexpr unary_function erf_via_macro {&::BOOST_MATH_NAMESPACE::erf};
constexpr unary_function erf_via_namespace {&::my_lib::math::erf};
static_assert(erf_via_macro == erf_via_namespace, "erf is one function");

} // namespace

int main()
{
    // Special functions
    check_same(BOOST_MATH_NAMESPACE::tgamma(5.0), my_lib::math::tgamma(5.0), 24.0, "tgamma");
    check_same(::BOOST_MATH_NAMESPACE::lgamma(10.0), ::my_lib::math::lgamma(10.0), 12.801827480081469611, "lgamma");
    check_same(BOOST_MATH_NAMESPACE::erf(0.5), my_lib::math::erf(0.5), 0.52049987781304653768, "erf");
    check_same(BOOST_MATH_NAMESPACE::cyl_bessel_j(0, 1.0), my_lib::math::cyl_bessel_j(0, 1.0), 0.76519768655796655145, "cyl_bessel_j");
    check_same(BOOST_MATH_NAMESPACE::beta(2.0, 3.0), my_lib::math::beta(2.0, 3.0), 1.0 / 12, "beta");
    check_same(BOOST_MATH_NAMESPACE::zeta(2.0), my_lib::math::zeta(2.0), 1.6449340668482264365, "zeta");
    check_same(bm::tgamma(5.0), my_lib::math::tgamma(5.0), 24.0, "tgamma through alias");
    check_same(user::nested::gamma_via_macro(5.0), user::nested::gamma_via_namespace(5.0), 24.0, "tgamma from nested namespace");

    const auto nan {std::numeric_limits<double>::quiet_NaN()};
    check_true(BOOST_MATH_NAMESPACE::isnan(nan) && my_lib::math::isnan(nan), "isnan");

    // Distributions, with the non-member functions called by qualified name
    check_same(BOOST_MATH_NAMESPACE::cdf(BOOST_MATH_NAMESPACE::normal {}, 0.0), my_lib::math::cdf(my_lib::math::normal {}, 0.0), 0.5, "normal cdf");
    check_same(BOOST_MATH_NAMESPACE::quantile(BOOST_MATH_NAMESPACE::students_t {5}, 0.975),
               my_lib::math::quantile(my_lib::math::students_t {5}, 0.975), 2.5705818356363146, "students_t quantile");

    // An object made through one spelling works with functions named through the other
    const BOOST_MATH_NAMESPACE::normal macro_normal {1.0, 2.0};
    check_same(my_lib::math::mean(macro_normal), BOOST_MATH_NAMESPACE::mean(macro_normal), 1.0, "mixed spelling mean");

    // Constants
    check_same(BOOST_MATH_NAMESPACE::constants::pi<double>(), my_lib::math::constants::pi<double>(), 3.141592653589793238, "pi");

    // Policies spelled one way passed to functions spelled the other way
    using macro_policy = BOOST_MATH_NAMESPACE::policies::policy<
        BOOST_MATH_NAMESPACE::policies::domain_error<BOOST_MATH_NAMESPACE::policies::ignore_error>,
        BOOST_MATH_NAMESPACE::policies::pole_error<BOOST_MATH_NAMESPACE::policies::ignore_error>>;
    using namespace_policy = my_lib::math::policies::policy<
        my_lib::math::policies::domain_error<my_lib::math::policies::ignore_error>,
        my_lib::math::policies::pole_error<my_lib::math::policies::ignore_error>>;
    static_assert(std::is_same<macro_policy, namespace_policy>::value, "policies");
    check_true(std::isnan(my_lib::math::tgamma(-1.0, macro_policy {})), "my_lib::math::tgamma with BOOST_MATH_NAMESPACE policy");
    check_true(std::isnan(BOOST_MATH_NAMESPACE::tgamma(-1.0, namespace_policy {})), "BOOST_MATH_NAMESPACE::tgamma with my_lib::math policy");

    // Tools
    const auto f {[](double x) { return x * x - 2; }};
    const auto macro_root {BOOST_MATH_NAMESPACE::tools::bisect(f, 0.0, 2.0, BOOST_MATH_NAMESPACE::tools::eps_tolerance<double> {50})};
    const auto namespace_root {my_lib::math::tools::bisect(f, 0.0, 2.0, my_lib::math::tools::eps_tolerance<double> {50})};
    check_same((macro_root.first + macro_root.second) / 2, (namespace_root.first + namespace_root.second) / 2, 1.4142135623730950488, "bisect");

    // Quadrature
    const auto g {[](double x) { return x * x; }};
    check_same(BOOST_MATH_NAMESPACE::quadrature::gauss_kronrod<double, 15>::integrate(g, 0.0, 1.0),
               my_lib::math::quadrature::gauss_kronrod<double, 15>::integrate(g, 0.0, 1.0), 1.0 / 3, "gauss_kronrod");

    #ifdef BOOST_MATH_QUALIFIED_NAME_TEST_CXX17
    check_same(BOOST_MATH_NAMESPACE::ccmath::sqrt(4.0), my_lib::math::ccmath::sqrt(4.0), 2.0, "ccmath sqrt");
    const std::vector<double> v {1.0, 2.0, 3.0, 4.0};
    check_same(BOOST_MATH_NAMESPACE::statistics::mean(v), my_lib::math::statistics::mean(v), 2.5, "statistics mean");
    #endif

    if (failures == 0)
    {
        std::cout << "BOOST_MATH_NAMESPACE and my_lib::math name the same entities\n";
    }
    return failures == 0 ? 0 : 1;
}
