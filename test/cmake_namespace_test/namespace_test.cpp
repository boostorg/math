// Copyright 2026 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// Built with BOOST_MATH_NAMESPACE_COMPONENTS=my_lib,math (https://github.com/boostorg/math/issues/769):
// the whole library must be usable through ::my_lib::math and nothing may be left in ::boost::math.

#include <boost/math/special_functions.hpp>
#include <boost/math/distributions.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/policies/policy.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/polynomial.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <cerrno>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>

#if (defined(_MSVC_LANG) && _MSVC_LANG >= 201703L) || __cplusplus >= 201703L
#  include <boost/math/ccmath/ccmath.hpp>
#  include <boost/math/statistics/univariate_statistics.hpp>
#  include <vector>
#  define BOOST_MATH_NAMESPACE_TEST_CXX17
#endif

static_assert(BOOST_MATH_NAMESPACE_LEVELS == 2, "my_lib,math has two components");
static_assert(std::is_same<::BOOST_MATH_NAMESPACE::normal, ::my_lib::math::normal>::value, "BOOST_MATH_NAMESPACE names my_lib::math");

// If anything still declares ::boost::math (or ::boost::math_detail) the lookups below are ambiguous
namespace boost {}
namespace probe
{
    namespace math { constexpr bool absent {true}; }
    namespace math_detail { constexpr bool absent {true}; }
}
namespace check { using namespace ::probe; using namespace ::boost; }
static_assert(check::math::absent, "::boost::math is still declared");
static_assert(check::math_detail::absent, "::boost::math_detail is still declared");

// The forwarding macros must work from a user namespace
namespace user_functions
{
    using c_policy = ::my_lib::math::policies::policy<
        ::my_lib::math::policies::domain_error<::my_lib::math::policies::errno_on_error>,
        ::my_lib::math::policies::pole_error<::my_lib::math::policies::errno_on_error>>;

    BOOST_MATH_DECLARE_SPECIAL_FUNCTIONS(c_policy)
}

namespace user_distributions
{
    using my_policy = ::my_lib::math::policies::policy<::my_lib::math::policies::promote_double<false>>;

    BOOST_MATH_DECLARE_DISTRIBUTIONS(double, my_policy)
}

namespace
{

int failures {0};

void check_close(double result, double expected, const char* what)
{
    const auto tolerance {1e-13 * (std::fabs(expected) > 1 ? std::fabs(expected) : 1.0)};
    if (!(std::fabs(result - expected) <= tolerance))
    {
        std::cerr << "FAIL " << what << ": got " << result << ", expected " << expected << '\n';
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

} // namespace

int main()
{
    // Special functions
    check_close(::my_lib::math::tgamma(5.0), 24.0, "tgamma");
    check_close(::my_lib::math::lgamma(10.0), 12.801827480081469611, "lgamma");
    check_close(::my_lib::math::erf(0.5), 0.52049987781304653768, "erf");
    check_close(::my_lib::math::cyl_bessel_j(0, 1.0), 0.76519768655796655145, "cyl_bessel_j");
    check_close(::my_lib::math::beta(2.0, 3.0), 1.0 / 12, "beta");
    check_close(::my_lib::math::zeta(2.0), 1.6449340668482264365, "zeta");
    check_true(::my_lib::math::isnan(std::numeric_limits<double>::quiet_NaN()), "isnan");
    check_true(::my_lib::math::isfinite(1.0), "isfinite");

    // Distributions, with cdf and quantile found by ADL
    const ::my_lib::math::normal standard_normal {};
    check_close(cdf(standard_normal, 0.0), 0.5, "normal cdf");
    check_close(quantile(::my_lib::math::students_t {5}, 0.975), 2.5705818356363146, "students_t quantile");

    // Constants
    check_close(::my_lib::math::constants::pi<double>(), 3.141592653589793238, "pi");

    // Tools
    const auto root {::my_lib::math::tools::bisect([](double x) { return x * x - 2; }, 0.0, 2.0,
                                                    ::my_lib::math::tools::eps_tolerance<double> {50})};
    check_close((root.first + root.second) / 2, 1.4142135623730950488, "bisect");
    const ::my_lib::math::tools::polynomial<double> p {{1.0, 2.0, 3.0}};
    check_close(p(2.0), 17.0, "polynomial");

    // Quadrature
    const auto area {::my_lib::math::quadrature::gauss_kronrod<double, 15>::integrate([](double x) { return x * x; }, 0.0, 1.0)};
    check_close(area, 1.0 / 3, "gauss_kronrod");

    // Default policy throws, from the relocated error handling
    bool threw {false};
    try
    {
        static_cast<void>(::my_lib::math::tgamma(0.0));
    }
    catch (const std::domain_error&)
    {
        threw = true;
    }
    check_true(threw, "tgamma(0) throws std::domain_error");

    // Policies through the forwarding macros
    errno = 0;
    check_true(std::isnan(user_functions::tgamma(-1.0)), "user policy tgamma(-1) is NaN");
    check_true(errno == EDOM, "user policy tgamma(-1) sets EDOM");
    check_close(user_functions::tgamma(5.0), 24.0, "user policy tgamma");
    check_close(cdf(user_distributions::normal {}, 0.0), 0.5, "user policy normal cdf");

    #ifdef BOOST_MATH_NAMESPACE_TEST_CXX17
    check_close(::my_lib::math::ccmath::sqrt(4.0), 2.0, "ccmath sqrt");
    const std::vector<double> v {1.0, 2.0, 3.0, 4.0};
    check_close(::my_lib::math::statistics::mean(v), 2.5, "statistics mean");
    #endif

    if (failures == 0)
    {
        std::cout << "All checks passed in namespace my_lib::math\n";
    }
    return failures == 0 ? 0 : 1;
}
