/*
 * Copyright Nick Thompson, 2020
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include "math_unit_test.hpp"
#include <boost/math/tools/simple_continued_fraction.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/next.hpp>
#include <boost/core/demangle.hpp>
#include <sstream>
#include <type_traits>
#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
using boost::multiprecision::float128;
#endif
#include <boost/multiprecision/cpp_bin_float.hpp>

using boost::math::tools::simple_continued_fraction;
using boost::multiprecision::cpp_bin_float_100;
using boost::math::constants::pi;

template<class Real>
void test_integral()
{
    for (int64_t i = -20; i < 20; ++i) {
        Real ii = i;
        auto cfrac = simple_continued_fraction<Real>(ii);
        auto const & a = cfrac.partial_denominators();
        CHECK_EQUAL(size_t(1), a.size());
        CHECK_EQUAL(i, a.front());
    }
}

template<class Real>
void test_halves()
{
    for (int64_t i = -20; i < 20; ++i) {
        Real x = i + Real(1)/Real(2);
        auto cfrac = simple_continued_fraction<Real>(x);
        auto const & a = cfrac.partial_denominators();
        CHECK_EQUAL(size_t(2), a.size());
        CHECK_EQUAL(i, a.front());
        CHECK_EQUAL(int64_t(2), a.back());
    }

    // We'll also test quarters; why not?
    for (int64_t i = -20; i < 20; ++i) {
        Real x = i + Real(1)/Real(4);
        auto cfrac = simple_continued_fraction<Real>(x);
        auto const & a = cfrac.partial_denominators();
        CHECK_EQUAL(size_t(2), a.size());
        CHECK_EQUAL(i, a.front());
        CHECK_EQUAL(int64_t(4), a.back());
    }

    for (int64_t i = -20; i < 20; ++i) {
        Real x = i + Real(1)/Real(8);
        auto cfrac = simple_continued_fraction<Real>(x);
        auto const & a = cfrac.partial_denominators();
        CHECK_EQUAL(size_t(2), a.size());
        CHECK_EQUAL(i, a.front());
        CHECK_EQUAL(int64_t(8), a.back());
    }

    for (int64_t i = -20; i < 20; ++i) {
        Real x = i + Real(3)/Real(4);
        auto cfrac = simple_continued_fraction<Real>(x);
        auto const & a = cfrac.partial_denominators();
        CHECK_EQUAL(size_t(3), a.size());
        CHECK_EQUAL(i, a.front());
        CHECK_EQUAL(int64_t(1), a[1]);
        CHECK_EQUAL(int64_t(3), a.back());
    }

    for (int64_t i = -20; i < 20; ++i) {
        Real x = i + Real(7)/Real(8);
        auto cfrac = simple_continued_fraction<Real>(x);
        auto const & a = cfrac.partial_denominators();
        CHECK_EQUAL(size_t(3), a.size());
        CHECK_EQUAL(i, a.front());
        CHECK_EQUAL(int64_t(1), a[1]);
        CHECK_EQUAL(int64_t(7), a.back());
    }
}

template<typename Real>
void test_simple()
{
    std::cout << "Testing rational numbers on type " << boost::core::demangle(typeid(Real).name()) << "\n";
    {
        Real x = Real(649)/200;
        // ContinuedFraction[649/200] = [3; 4, 12, 4]
        auto cfrac = simple_continued_fraction(x);
        auto const & a = cfrac.partial_denominators();
        CHECK_EQUAL(size_t(4), a.size());
        CHECK_EQUAL(int64_t(3), a[0]);
        CHECK_EQUAL(int64_t(4), a[1]);
        CHECK_EQUAL(int64_t(12), a[2]);
        CHECK_EQUAL(int64_t(4), a[3]);
    }

    {
        Real x = Real(415)/Real(93);
        // [4; 2, 6, 7]:
        auto cfrac = simple_continued_fraction(x);
        auto const & a = cfrac.partial_denominators();
        CHECK_EQUAL(size_t(4), a.size());
        CHECK_EQUAL(int64_t(4), a[0]);
        CHECK_EQUAL(int64_t(2), a[1]);
        CHECK_EQUAL(int64_t(6), a[2]);
        CHECK_EQUAL(int64_t(7), a[3]);
    }

}

template<typename Real>
void test_khinchin()
{
    // These are simply sanity checks; the convergence is too slow otherwise:
    auto cfrac = simple_continued_fraction(pi<Real>());
    auto K0 = cfrac.khinchin_geometric_mean();
    CHECK_MOLLIFIED_CLOSE(Real(2.6854520010), K0, 0.1);
    auto Km1 = cfrac.khinchin_harmonic_mean();
    CHECK_MOLLIFIED_CLOSE(Real(1.74540566240), Km1, 0.1);
    
    using std::sqrt;
    auto rt_cfrac = simple_continued_fraction(sqrt(static_cast<Real>(2)));
    K0 = rt_cfrac.khinchin_geometric_mean();
    CHECK_ULP_CLOSE(Real(2), K0, 10);
    Km1 = rt_cfrac.khinchin_harmonic_mean();
    CHECK_ULP_CLOSE(Real(2), Km1, 10);
}

template<typename Real>
void test_just_below_integer()
{
    // x = n - ulp is within tolerance of the convergent n after one step, which
    // gives [n - 1; 1].  The canonical form of that is [n]:
    for (int64_t n : {2, 7, 20, -3, -19}) {
        Real x = boost::math::float_prior(static_cast<Real>(n));
        auto cfrac = simple_continued_fraction<Real>(x);
        auto const & a = cfrac.partial_denominators();
        CHECK_EQUAL(size_t(1), a.size());
        CHECK_EQUAL(n, a.front());
    }
}

void test_canonical_form_overflow()
{
    // Just below 2^31 the expansion is [INT32_MAX; 1].  Canonicalizing it would
    // overflow int32_t, so the longer (still valid) form must be kept:
    double x = std::nextafter(2147483648.0, 0.0);
    auto cfrac = simple_continued_fraction<double, int32_t>(x);
    auto const & a = cfrac.partial_denominators();
    CHECK_EQUAL(size_t(2), a.size());
    CHECK_EQUAL((std::numeric_limits<int32_t>::max)(), a[0]);
    CHECK_EQUAL(int32_t(1), a[1]);
    // With a wider integer type it is canonicalized as usual:
    auto wide = simple_continued_fraction<double, int64_t>(x);
    CHECK_EQUAL(size_t(1), wide.partial_denominators().size());
    CHECK_EQUAL(int64_t(2147483648), wide.partial_denominators().front());
}

template<typename Real>
void test_value_semantics()
{
    static_assert(std::is_copy_assignable_v<simple_continued_fraction<Real>>);
    static_assert(std::is_move_assignable_v<simple_continued_fraction<Real>>);
    auto cfrac = simple_continued_fraction<Real>(Real(649)/200);
    cfrac = simple_continued_fraction<Real>(Real(415)/93);
    auto const & a = cfrac.partial_denominators();
    CHECK_EQUAL(size_t(4), a.size());
    CHECK_EQUAL(int64_t(4), a[0]);
    CHECK_EQUAL(int64_t(7), a[3]);
}

template<typename Real>
void test_output()
{
    // Printing a const object must compile, and must not change the stream's state:
    const auto cfrac = simple_continued_fraction<Real>(Real(649)/200);
    std::ostringstream oss;
    oss.precision(3);
    oss << cfrac;
    CHECK_EQUAL(std::string("[3; 4, 12, 4]"), oss.str());
    CHECK_EQUAL(std::streamsize(3), oss.precision());

    std::ostringstream oss2;
    oss2 << simple_continued_fraction<Real>(Real(2));
    CHECK_EQUAL(std::string("[2]"), oss2.str());
}

int main()
{
    test_integral<float>();
    test_integral<double>();
    test_integral<long double>();
    test_integral<cpp_bin_float_100>();

    test_halves<float>();
    test_halves<double>();
    test_halves<long double>();
    test_halves<cpp_bin_float_100>();

    test_simple<float>();
    test_simple<double>();
    test_simple<long double>();
    test_simple<cpp_bin_float_100>();
    
    test_khinchin<cpp_bin_float_100>();

    test_just_below_integer<float>();
    test_just_below_integer<double>();
    test_just_below_integer<long double>();
    test_just_below_integer<cpp_bin_float_100>();

    test_canonical_form_overflow();

    test_value_semantics<double>();
    test_value_semantics<cpp_bin_float_100>();

    test_output<double>();
    test_output<cpp_bin_float_100>();
    
    #ifdef BOOST_HAS_FLOAT128
    test_integral<float128>();
    test_halves<float128>();
    test_simple<float128>();
    test_khinchin<float128>();
    #endif
    return boost::math::test::report_errors();
}
