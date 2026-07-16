// Copyright 2026 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

module;

// The module is always built in standalone mode so that no entity from another
// Boost library can attach itself to the boost.math module.
#ifndef BOOST_MATH_STANDALONE
#  define BOOST_MATH_STANDALONE
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#  define BOOST_MATH_BUILD_MODULE
#endif

// Marks this translation unit as the module interface unit (as opposed to a
// consumer that imports the module). Some headers key off this to declare
// entities here that consumers instead receive through the import.
#define BOOST_MATH_INTERFACE_UNIT

// Platform intrinsic headers are not part of the standard library module, so
// they are always brought in here. The conditions mirror the include sites in
// special_functions/next.hpp and special_functions/detail/lanczos_sse2.hpp.
#if !defined(_CRAYC) && (!defined(__GNUC__) || (__GNUC__ > 3) || ((__GNUC__ == 3) && (__GNUC_MINOR__ > 3)))
#  if (defined(_M_IX86_FP) && (_M_IX86_FP >= 2)) || defined(__SSE2__)
#    include <xmmintrin.h>
#    include <emmintrin.h>
#  endif
#endif

// import std exports declarations but not macros. The headers below provide the
// feature-test macros and the object-like macros the library uses (assert,
// errno, FP_ classification, HUGE_VAL, FLT_/DBL_ limits, INT64_C, ...) so they
// are included regardless of whether import std is used.
#include <version>
#include <cassert>
#include <cerrno>
#include <cfenv>
#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdint>

// The macro hub and its helper entities are deliberately kept in the global
// module so that consumer translation units may also include them textually
// (for example through test/math_unit_test.hpp).
#include <boost/math/tools/config.hpp>
#include <boost/math/tools/assert.hpp>

// When the standard library module is available these are provided by import
// std, otherwise they are supplied here in the global module fragment.
#ifndef BOOST_MATH_USE_STD_MODULE

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cctype>
#include <charconv>
#include <chrono>
#include <cinttypes>
#include <complex>
#include <concepts>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <exception>
#include <fstream>
#include <functional>
#include <future>
#include <initializer_list>
#include <iomanip>
#include <ios>
#include <iosfwd>
#include <iostream>
#include <istream>
#include <iterator>
#include <limits>
#include <list>
#include <locale>
#include <map>
#include <memory>
#include <mutex>
#include <numeric>
#include <ostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <type_traits>
#include <typeinfo>
#include <utility>
#include <valarray>
#include <vector>

#ifndef BOOST_MATH_NO_CXX17_HDR_EXECUTION
#  include <execution>
#endif

#endif // BOOST_MATH_USE_STD_MODULE

export module boost.math;

#ifdef BOOST_MATH_USE_STD_MODULE
import std;
#endif

#ifdef _MSC_VER
#  pragma warning( push )
#  pragma warning( disable : 5244 )
#elif defined(__clang__)
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Winclude-angled-in-module-purview"
#endif

extern "C++" {

// Foundation components
#include <boost/math/constants/constants.hpp>
#include <boost/math/policies/policy.hpp>
#include <boost/math/policies/error_handling.hpp>

// Special functions (umbrella)
#include <boost/math/special_functions.hpp>

// Statistical distributions (umbrella)
#include <boost/math/distributions.hpp>

// Constexpr cmath (umbrella)
#include <boost/math/ccmath/ccmath.hpp>

// Tools (curated public set)
#include <boost/math/tools/agm.hpp>
#include <boost/math/tools/centered_continued_fraction.hpp>
#include <boost/math/tools/cohen_acceleration.hpp>
#include <boost/math/tools/color_maps.hpp>
#include <boost/math/tools/condition_numbers.hpp>
#include <boost/math/tools/cubic_roots.hpp>
#include <boost/math/tools/engel_expansion.hpp>
#include <boost/math/tools/estrin.hpp>
#include <boost/math/tools/fraction.hpp>
#include <boost/math/tools/luroth_expansion.hpp>
#include <boost/math/tools/minima.hpp>
#include <boost/math/tools/norms.hpp>
#include <boost/math/tools/polynomial.hpp>
#include <boost/math/tools/quartic_roots.hpp>
#include <boost/math/tools/rational.hpp>
#include <boost/math/tools/recurrence.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/series.hpp>
#include <boost/math/tools/simple_continued_fraction.hpp>
#include <boost/math/tools/toms748_solve.hpp>

// Quadrature
#include <boost/math/quadrature/exp_sinh.hpp>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/quadrature/naive_monte_carlo.hpp>
#include <boost/math/quadrature/ooura_fourier_integrals.hpp>
#include <boost/math/quadrature/sinh_sinh.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/quadrature/wavelet_transforms.hpp>

// Interpolators (the deprecated boost::math namespace shims and the FFTW
// dependent cardinal_trigonometric interpolator are intentionally excluded)
#include <boost/math/interpolators/barycentric_rational.hpp>
#include <boost/math/interpolators/bezier_polynomial.hpp>
#include <boost/math/interpolators/bilinear_uniform.hpp>
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
#include <boost/math/interpolators/cardinal_quadratic_b_spline.hpp>
#include <boost/math/interpolators/cardinal_quintic_b_spline.hpp>
#include <boost/math/interpolators/catmull_rom.hpp>
#include <boost/math/interpolators/cubic_hermite.hpp>
#include <boost/math/interpolators/makima.hpp>
#include <boost/math/interpolators/pchip.hpp>
#include <boost/math/interpolators/quintic_hermite.hpp>
#include <boost/math/interpolators/septic_hermite.hpp>
#include <boost/math/interpolators/vector_barycentric_rational.hpp>
#include <boost/math/interpolators/whittaker_shannon.hpp>

// Statistics
#include <boost/math/statistics/anderson_darling.hpp>
#include <boost/math/statistics/bivariate_statistics.hpp>
#include <boost/math/statistics/chatterjee_correlation.hpp>
#include <boost/math/statistics/linear_regression.hpp>
#include <boost/math/statistics/ljung_box.hpp>
#include <boost/math/statistics/runs_test.hpp>
#include <boost/math/statistics/signal_statistics.hpp>
#include <boost/math/statistics/t_test.hpp>
#include <boost/math/statistics/univariate_statistics.hpp>
#include <boost/math/statistics/z_test.hpp>

// Optimization (cma_es is excluded: it requires Eigen)
#include <boost/math/optimization/differential_evolution.hpp>
#include <boost/math/optimization/gradient_descent.hpp>
#include <boost/math/optimization/gradient_optimizers.hpp>
#include <boost/math/optimization/jso.hpp>
#include <boost/math/optimization/lbfgs.hpp>
#include <boost/math/optimization/minimizer.hpp>
#include <boost/math/optimization/nesterov.hpp>
#include <boost/math/optimization/random_search.hpp>

// Differentiation. finite_difference and lanczos_smoothing must precede
// autodiff: their unqualified detail:: references become ambiguous once
// autodiff's inline namespace introduces a second differentiation::detail.
#include <boost/math/differentiation/finite_difference.hpp>
#include <boost/math/differentiation/lanczos_smoothing.hpp>
#include <boost/math/differentiation/autodiff.hpp>

// Algebra types and complex inverse trigonometric functions
#include <boost/math/complex.hpp>
#include <boost/math/quaternion.hpp>
#include <boost/math/octonion.hpp>

} // extern "C++"

#ifdef _MSC_VER
#  pragma warning( pop )
#elif defined(__clang__)
#  pragma clang diagnostic pop
#endif
