//  (C) Copyright John Maddock 2022.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

module;

#define BOOST_MATH_POLICIES_AS_MODULE
#define BOOST_MATH_AS_MODULE
#define BOOST_MATH_STANDALONE

// For config:
#include <algorithm>  // for min and max
#include <limits>
#include <cmath>
#include <climits>
#include <cfloat>
#if (defined(macintosh) || defined(__APPLE__) || defined(__APPLE_CC__))
#  include <math.h>
#endif
#include <type_traits>
#include <cfenv>

// For error_handling.hpp and policies.hpp:
#include <iomanip>
#include <string>
#include <cstring>
#include <typeinfo>
#include <cerrno>
#include <complex>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <boost/math/tools/mp.hpp>
#include <boost/math/tools/config.hpp>
#include <boost/math/tools/throw_exception.hpp>
#include <sstream>
#include <limits>
#include <cstdint>
#include <cstddef>
#include <boost/math/tools/assert.hpp>
#include <limits>
#include <climits>
#include <cstdint>
#include <cfloat> // LDBL_MANT_DIG

// for fptraits.hpp:
#if (__cplusplus >= 202002L || _MSVC_LANG >= 202002L)
#if __has_include(<bit>)
#include <bit>
#endif
#endif

// for fpclassify.hpp:
#include <float.h>
#ifdef BOOST_MATH_USE_FLOAT128
#ifdef __has_include
#if  __has_include("quadmath.h")
#include "quadmath.h"
#define BOOST_MATH_HAS_QUADMATH_H
#endif
#endif
#endif

// for next.hpp:
#if !defined(_CRAYC) && !defined(__CUDACC__) && (!defined(__GNUC__) || (__GNUC__ > 3) || ((__GNUC__ == 3) && (__GNUC_MINOR__ > 3)))
#if (defined(_M_IX86_FP) && (_M_IX86_FP >= 2)) || defined(__SSE2__)
#include "xmmintrin.h"
#define BOOST_MATH_CHECK_SSE2
#endif
#endif

// for roots:
#include <tuple>

// for rational.hpp:
#include <array>



export module boost.math.core;

#include <boost/math/policies/policy.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/tools/promotion.hpp> // for argument promotion.
#include <boost/math/special_functions/sign.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/special_functions/trunc.hpp>
#include <boost/math/special_functions/next.hpp>
#include <boost/math/tools/toms748_solve.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/rational.hpp>

