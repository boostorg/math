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

// The 18 polynomial Coefficients
template <typename Real>
struct coefficients
{
    Real p00;
    Real p10;
    Real p01;
    Real p20;
    Real p11;
    Real p02;
    Real p30;
    Real p40;
    Real p50;
    Real p03;
    Real p04;
    Real p05;
    Real p41;
    Real p14;
    Real p21;
    Real p31;
    Real p12;
    Real p13;
};

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

// Determine the Coefficients of the Polynomial. EQNs A-20 through 25 calculated in sequential order
template <typename Real>
coefficients polynomial_coefficient(Real z00, std::tuple<Real, Real, Real, Real, Real> z00_partial_derivatives,
                                    Real z10, std::tuple<Real, Real, Real, Real, Real> z10_partial_derivatives,
                                    Real z01, std::tuple<Real, Real, Real, Real, Real> z01_partial_derivatives,
                                    std::tuple<Real, Real, Real> vectors)
{
    using std::cos;
    
    coefficients<Real> c;

    // EQN A-20
    const Real zu00  = std::get<0>(z00_partial_derivatives);
    const Real zv00  = std::get<1>(z00_partial_derivatives);
    const Real zuu00 = std::get<2>(z00_partial_derivatives);
    const Real zuv00 = std::get<3>(z00_partial_derivatives);
    const Real zvv00 = std::get<4>(z00_partial_derivatives);

    c.p00 = z00;
    c.p10 = zu00;
    c.p01 = zv00;
    c.p20 = zuu00 / 2;
    c.p11 = zuv00;
    c.p02 = zvv00 / 2;

    // EQN A-21
    const Real zu10  = std::get<0>(z10_partial_derivatives);
    const Real zv10  = std::get<1>(z10_partial_derivatives);
    const Real zuu10 = std::get<2>(z10_partial_derivatives);
    const Real zuv10 = std::get<3>(z10_partial_derivatives);
    const Real zvv10 = std::get<4>(z10_partial_derivatives);

    c.p30 = (10*z10 - 8*zu10 + zuu10 - 20*c.p00 - 12*c.p10 - 6*c.p20) / 2;
    c.p40 = -15*z10 + 7*zu10 - zuu10 + 15*c.p00 + 8*c.p10 + 3*c.p20;
    c.p50 = (12*z10 - 6*zu10 + zuu10 - 12*c.p00 - 6*c.p10 - 2*c.p20) / 2;

    // EQN A-22
    const Real zu01  = std::get<0>(z01_partial_derivatives);
    const Real zv01  = std::get<1>(z01_partial_derivatives);
    const Real zuu01 = std::get<2>(z01_partial_derivatives);
    const Real zuv01 = std::get<3>(z01_partial_derivatives);
    const Real zvv01 = std::get<4>(z01_partial_derivatives);

    c.p03 = (20*z01 - 8*zv01 + zvv01 - 20*c.p00 - 12*c.p01 - 6*c.p02) / 2;
    c.p04 = -15*z01 + 7*zv01 - zvv01 + 15*c.p00 + 8*c.p01 + 3*c.p02;
    c.p05 = (12*z01 - 6*zv01 + zvv01 -12*c.p00 - 6*c.p01 - 2*c.p02) / 2;

    // EQN A-23
    const Real lu = std::get<0>(vectors);
    const Real lv = std::get<1>(vectors);
    const Real theta_uv = std::get<2>(vectors);

    c.p41 = (5*lv*cos(theta_uv)/lu)*c.p50;
    c.p14 = (5*lu*cos(theta_uv)/lv)*c.p05;

    // EQN A-24
    c.p21 = 3*zv10 - zuv10 - 3*c.p01 - 2*c.p11 + c.p41;
    c.p31 = -2*zv10 + zuv10 + 2*c.p01 + c.p11 - 2*c.p41;

    // EQN A-25
    c.p12 = 3*zu01 - zuv01 - 3*c.p10 - 2*c.p11 + c.p14;
    c.p13 = -2*zu01 + zuv01 + 2*c.p10 + c.p11 - 2*c.p14;

    return c;
}

}}}} // namespaces

#endif // BOOST_MATH_INTERPOLATORS_2D_AKIMA_DETIAL_HPP
