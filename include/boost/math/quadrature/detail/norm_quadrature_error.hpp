// Copyright Nick Thompson, 2026
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
#ifndef BOOST_MATH_QUADRATURE_DETAIL_NORM_QUADRATURE_ERROR_HPP
#define BOOST_MATH_QUADRATURE_DETAIL_NORM_QUADRATURE_ERROR_HPP
#ifndef BOOST_MATH_BUILD_MODULE
#include <cstddef>
#endif
#include <boost/math/tools/config.hpp>

BOOST_MATH_NAMESPACE_BEGIN namespace quadrature { namespace detail {
template<class K, class Real>
K norm_quadrature_error(const K& zero, Real invalid, Real* error, Real* l1, std::size_t* levels = nullptr)
{
    if (error) *error = invalid;
    if (l1) *l1 = invalid;
    if (levels) *levels = 0;
    K result = zero * invalid;
    return result;
}
}} BOOST_MATH_NAMESPACE_END
#endif
