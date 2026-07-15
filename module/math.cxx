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

#ifdef _MSC_VER
#  pragma warning( pop )
#elif defined(__clang__)
#  pragma clang diagnostic pop
#endif
