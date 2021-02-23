// Copyright Matt Borland, 2021
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Source: https://www.its.bldrdoc.gov/DownloadPublications.ashx?id=OT-75-70.pdf
//
// Description:
//
// In this method the x-y plane is divided into a number of triangular cells;
// each having projections of three data points in the in the plane as its vertexes,
// and a bivariate fifth-degree polynomial in x and y is applied to each triangular cell.

#ifndef BOOST_MATH_INTERPOLATORS_2D_AKIMA_HPP
#define BOOST_MATH_INTERPOLATORS_2D_AKIMA_HPP

#include <memory>
#include <utility>
#include <cstddef>
#include <type_traits>
#include <boost/math/interpolators/detail/2d_akima_detail.hpp>

namespace boost { namespace math { namespace interpolators { 

namespace detail {
// std::make_unique not implemented until C++14
template <typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

} // Namespace detail

template <typename Matrix, typename VectorType = typename Matrix::value_type, typename Real = typename VectorType::value_type>
class akima_2d_uniform
{
private:
    std::unique_ptr<detail::akima_2d_uniform_impl<Matrix, VectorType, Real>> impl_;

public:
    akima_2d_uniform(Matrix&& data, const std::size_t n) : 
        impl_ {detail::make_unique<detail::akima_2d_uniform_impl<Matrix, VectorType, Real>>(std::move(data), n)} {}
};

}}} // namespaces

#endif // BOOST_MATH_INTERPOLATORS_2D_AKIMA_HPP
