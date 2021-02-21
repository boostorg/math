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

#include <memory>
#include <execution>
#include <utility>
#include <cstddef>
#include <boost/math/interpolators/detail/2d_akima_detail.hpp>

namespace boost::math::interpolators { 

template <typename ExecutionPolicy = decltype(std::execution::seq), typename ForwardContainer = void>
class akima
{
private:
    std::unique_ptr<detail::akima_impl<ExecutionPolicy, ForwardContainer>> impl_;

public:
    akima(ExecutionPolicy&& exec, const ForwardContainer& x, const ForwardContainer& y, const ForwardContainer& z, const std::size_t n) :
        impl_ {std::make_unique<detail::akima_impl<ExecutionPolicy, ForwardContainer>>(std::move(exec), x, y, z, n)} {}

    akima(const ForwardContainer& x, const ForwardContainer& y, const ForwardContainer& z, const std::size_t n) : 
        impl_ {std::make_unique<detail::akima_impl<ExecutionPolicy, ForwardContainer>>(std::move(std::execution::seq), x, y, z, n)} {}
};

} // namespace boost::math::interpolators
