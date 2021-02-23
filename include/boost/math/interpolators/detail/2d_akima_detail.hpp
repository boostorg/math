// Copyright Matt Borland, 2021
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Source: https://www.its.bldrdoc.gov/DownloadPublications.ashx?id=OT-75-70.pdf

#ifndef BOOST_MATH_INTERPOLATORS_2D_AKIMA_DETIAL_HPP
#define BOOST_MATH_INTERPOLATORS_2D_AKIMA_DETIAL_HPP

#include <type_traits>
#include <stdexcept>
#include <utility>
#include <algorithm>
#include <iostream>
#include <cstddef>
#include <cmath>

namespace boost { namespace math { namespace interpolators { namespace detail {

template <typename Matrix, typename VectorType = typename Matrix::value_type, typename Real = typename VectorType::value_type>
class akima_2d_uniform_impl
{
private:    
    const Matrix data_;
    const std::size_t n_interpolation_points_;
    VectorType x_interpolation_points_;
    VectorType y_interpolation_points_;
    VectorType z_interpolation_points_;

    void validate_inputs() const;
    VectorType apply_uniform_grid(const VectorType& vals);

public:
    akima_2d_uniform_impl(Matrix&& data, std::size_t n) : data_{std::move(data)}, n_interpolation_points_{n}
    {
        validate_inputs();
        x_interpolation_points_ = apply_uniform_grid(data_[0]);
        y_interpolation_points_ = apply_uniform_grid(data_[1]);
        z_interpolation_points_.reserve(n_interpolation_points_);
    }
};

template <typename Matrix, typename VectorType, typename Real>
void akima_2d_uniform_impl<Matrix, VectorType, Real>::validate_inputs() const
{
    const auto data_x_size = std::size(data_[0]);
    const auto data_y_size = std::size(data_[1]);
    const auto data_z_size = std::size(data_[2]);

    if(data_x_size != data_y_size || data_x_size != data_z_size || data_y_size != data_z_size)
    {
        throw std::domain_error("X, Y, and Z must all be of the same size");
    }

    if(data_x_size < 4)
    {
        throw std::domain_error("X, Y, and Z must all have at least 4 values");
    }

    if(n_interpolation_points_ == 0)
    {
        throw std::domain_error("N must be at least 1");
    }
}

// Generates a uniform grid in the x-y plane of the dataset
template <typename Matrix, typename VectorType, typename Real>
VectorType akima_2d_uniform_impl<Matrix, VectorType, Real>::apply_uniform_grid(const VectorType& vals)
{
    using ForwardIt = decltype(std::begin(vals));
    
    const std::pair<ForwardIt, ForwardIt> extremes = std::minmax_element(std::begin(vals), std::end(vals));
    const Real step = (*extremes.second - *extremes.first) / (n_interpolation_points_ - 1);

    VectorType grid(n_interpolation_points_);

    for(std::size_t i = 0; i < n_interpolation_points_; ++i)
    {
        grid[i] = *extremes.first + step * i;
    }

    return grid;
}

// EQN A-4
template <typename Real>
std::tuple<Real, Real, Real, Real> coordinate_transform(Real x1, Real y1, Real x2, Real y2, Real x3, Real y3)
{
    const Real a = x2 - x1;
    const Real b = x3 - x1;
    const Real c = y2 - y1;
    const Real d = y3 - y1;

    return std::make_tuple(a, b, c, d);
}

// EQN A-6
template <typename Real>
std::tuple<Real, Real, Real, Real, Real> parital_derivative_transform(std::tuple<Real, Real, Real, Real> coefficients, 
                                                                      Real zx, Real zy, Real zxx, Real zxy, Real zyy)
{
    const Real a = std::get<0>(coefficients);
    const Real b = std::get<1>(coefficients);
    const Real c = std::get<2>(coefficients);
    const Real d = std::get<3>(coefficients);

    const Real zu = a*zx + c*zy;
    const Real zv = b*zx + d*zy;
    const Real zuu = a*a*zxx + 2*a*c*zxy+c*c*zyy;
    const Real zuv = a*b*zxx + (a*d + b*c)*zxy + c*d*zyy;
    const Real zvv = b*b*zxx + 2*b*d*zxy + d*d*zyy;

    return std::make_tuple(zu, zv, zuu, zuv, zvv);
}

// EQN A-9
template <typename Real>
std::tuple<Real, Real, Real> unit_vectors(std::tuple<Real, Real, Real, Real> coefficients)
{
    using std::atan;
    
    const Real a = std::get<0>(coefficients);
    const Real b = std::get<1>(coefficients);
    const Real c = std::get<2>(coefficients);
    const Real d = std::get<3>(coefficients);

    const Real lu = a*a + c*c;
    const Real lv = b*b + d*d;
    const Real theta_uv = atan(d/b) - atan(c/a);

    return std::make_tuple(lu, lv, theta_uv);
}

// EQN A-17
template <typename Real>
std::tuple<Real, Real, Real, Real> transform_coefficients(std::tuple<Real, Real, Real, Real> uv)
{
    using std::sin;
    using std::cos;

    const Real lu = std::get<0>(uv);
    const Real lv = std::get<1>(uv);
    const Real theta_uv = std::get<2>(uv);
    const Real theta_us = std::get<3>(uv);
    const Real sin_theta_uv = sin(theta_uv);

    const Real A = sin(theta_uv - theta_us) / (lu * sin_theta_uv);
    const Real B = -cos(theta_uv - theta_us) / (lu * sin_theta_uv);
    const Real C = sin(theta_us) / (lv * sin_theta_uv);
    const Real D = cos(theta_us) / (lv * sin_theta_uv);

    return std::make_tuple(A, B, C, D);
}

}}}} // namespaces

#endif // BOOST_MATH_INTERPOLATORS_2D_AKIMA_DETIAL_HPP
