/*
 * Copyright Thomas Dybdahl Ahle, Nick Thompson, Matt Borland, 2023
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef BOOST_MATH_TOOLS_ESTRIN_HPP
#define BOOST_MATH_TOOLS_ESTRIN_HPP

#include <array>
#include <vector>

namespace boost::math::tools {

template<typename RandomAccessContainer1, typename RandomAccessContainer2, typename RealOrComplex>
inline auto estrin(RandomAccessContainer1 const & coeffs, RandomAccessContainer2& scratch, RealOrComplex z) {
  static_assert(std::is_same_v<typename RandomAccessContainer2::value_type, RealOrComplex>, "The value type of the scratch space must be the same as the abscissa.");
  auto n = coeffs.size();
  for (decltype(n) i = 0; i < n / 2; i++) {
    scratch[i] = coeffs[2 * i] + coeffs[2 * i + 1] * z;
  }
  if (n & 1) {
    scratch[n / 2] = coeffs[n - 1];
  }
  auto m = (n + 1) / 2;

  while (m != 1) {
    z = z * z;
    for (decltype(n) i = 0; i < m / 2; i++) {
      scratch[i] = scratch[2 * i] + scratch[2 * i + 1] * z;
    }
    if (m & 1) {
      scratch[m / 2] = scratch[m - 1];
    }
    m = (m + 1) / 2;
  }
  return scratch[0];
}


// The std::array template specialization doesn't need to allocate:
template <typename RealOrComplex1, size_t n, typename RealOrComplex2>
inline RealOrComplex2 estrin(const std::array<RealOrComplex1, n> &coeffs, RealOrComplex2 z) {
  std::array<RealOrComplex2, (n + 1) / 2> ds;
  return estrin(coeffs, ds, z);
}

template <typename RandomAccessContainer, typename RealOrComplex>
inline RealOrComplex estrin(const RandomAccessContainer &coeffs, RealOrComplex z) {
  auto n = coeffs.size();
  // Normally, I'd make `ds` a RandomAccessContainer, but its value type needs to be RealOrComplex,
  // and the value_type of the passed RandomAccessContainer can just be Real.
  // Allocation of the std::vector is not ideal, but I have no other ideas at the moment:
  std::vector<RealOrComplex> ds((n+1)/2);
  return estrin(coeffs, ds, z);
}
}

#endif
