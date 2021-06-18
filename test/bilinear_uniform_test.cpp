/*
 * Copyright Nick Thompson, 2021
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include "math_unit_test.hpp"
#include <numeric>
#include <random>
#include <array>
#include <boost/core/demangle.hpp>
#include <boost/math/interpolators/bilinear_uniform.hpp>
#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
using boost::multiprecision::float128;
#endif

using boost::math::interpolators::bilinear_uniform;

template<class Real>
void test_four_values()
{
    Real x0 = 0;
    Real y0 = 0;
    Real dx = 1;
    Real dy = 1;
    Real f = 1.5;
    std::vector<Real> v(2*2, f);
    auto v_copy = v;
    auto ub = bilinear_uniform<decltype(v)>(std::move(v_copy), 2, 2, dx, dy, x0, y0);
    for (Real x = x0; x <= x0 + dx; x += dx/8) {
        for (Real y = y0; y <= y0 + dx; y += dy/8) {
            CHECK_ULP_CLOSE(f, ub(x, y), 1);
        }
    }
    std::cout << ub << "\n";

    // Now we test the unit square:
    std::random_device rd;
    std::uniform_real_distribution<Real> dis(1,2);

    int i = 0;
    while (i++ < 300) {
        v[0] = dis(rd);
        v[1] = dis(rd);
        v[2] = dis(rd);
        v[3] = dis(rd);

        // See https://en.wikipedia.org/wiki/Bilinear_interpolation, section: Unit square
        auto f = [&v](Real x, Real y) {
            return v[0]*(1-x)*(1-y) + v[1]*x*(1-y) + v[2]*(1-x)*y + v[3]*x*y;
        };

        auto v_copy = v;
        ub = bilinear_uniform<decltype(v_copy)>(std::move(v_copy), 2, 2, dx, dy, x0, y0);
        for (Real x = x0; x <= x0 + dx; x += dx/16) {
            for (Real y = y0; y <= y0 + dx; y += dy/16) {
                CHECK_ULP_CLOSE(f(x,y), ub(x, y), 3);
            }
        }
        std::cout << ub << "\n";
    }
}

template<typename Real>
void test_linear()
{
    std::random_device rd;
    std::uniform_real_distribution<Real> dis(1,2);
    std::array<Real, 4> a{dis(rd), dis(rd), dis(rd), dis(rd)};
    auto f = [&a](Real x, Real y) {
        return a[0] + a[1]*x + a[2]*y + a[3]*x*y;
    };

    Real dx = dis(rd);
    Real dy = dis(rd);
    Real x0 = dis(rd);
    Real y0 = dis(rd);
    for (int rows = 2; rows < 10; ++rows) {
        for (int cols = 2; cols < 10; ++cols) {
            std::vector<Real> v(rows*cols);
            for (int i = 0; i < cols; ++i) {
                for (int j = 0; j < rows; ++j) {
                    v[j*cols + i] = f(x0 + i*dx, y0 + j*dy);
                }
            }
            auto ub = bilinear_uniform<decltype(v)>(std::move(v), rows, cols, dx, dy, x0, y0);

            for (Real x = x0; x <= x0 + cols*dx; x += dx/8) {
                for (Real y = y0; y <= y0 + rows*dy; y += dy/8) {
                    CHECK_ULP_CLOSE(f(x,y), ub(x, y), 3);
                }
            }
        }
    }
}

int main()
{
    test_four_values<float>();
    test_four_values<double>();
    test_four_values<long double>();
    test_linear<long double>();
    return boost::math::test::report_errors();
}
