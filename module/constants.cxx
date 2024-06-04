//  (C) Copyright John Maddock 2022.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

module;

#define BOOST_MATH_CONSTANTS_AS_MODULE
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
#include <utility>
#include <cfenv>
#include <boost/math/tools/config.hpp>

export module boost.math.constants;

import boost.math.core;

#include <boost/math/constants/constants.hpp>
