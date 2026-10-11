/*
 * Copyright Nick Thompson, 2017
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
#include <array>
#include <random>
#include <boost/cstdfloat.hpp>
#include <boost/type_index.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/interpolators/catmull_rom.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include "math_unit_test.hpp"

#if __has_include(<stdfloat>)
#  include <stdfloat>
#endif

using std::abs;
using boost::multiprecision::cpp_bin_float_50;
using boost::math::catmull_rom;

template<class Real>
void test_alpha_distance()
{
    std::array<Real, 3> v1 = {0,0,0};
    std::array<Real, 3> v2 = {1,0,0};
    Real alpha = 0.5;
    Real d = boost::math::detail::alpha_distance<std::array<Real, 3>>(v1, v2, alpha);
    CHECK_ULP_CLOSE(1, d, 2);

    d = boost::math::detail::alpha_distance<std::array<Real, 3>>(v1, v2, 0.0);
    CHECK_ULP_CLOSE(1, d, 2);

    d = boost::math::detail::alpha_distance<std::array<Real, 3>>(v1, v2, 1.0);
    CHECK_ULP_CLOSE(1, d, 2);

    v2[0] = 2;
    d = boost::math::detail::alpha_distance<std::array<Real, 3>>(v1, v2, alpha);
    CHECK_ULP_CLOSE(pow(2, (Real)1/ (Real) 2), d, 2);

    d = boost::math::detail::alpha_distance<std::array<Real, 3>>(v1, v2, 0.0);
    CHECK_ULP_CLOSE(1, d, 2);

    d = boost::math::detail::alpha_distance<std::array<Real, 3>>(v1, v2, 1.0);
    CHECK_ULP_CLOSE(2, d, 2);
}


template<class Real>
void test_linear()
{
    std::cout << "Testing that the Catmull-Rom spline interpolates linear functions correctly on type "
              << boost::typeindex::type_id<Real>().pretty_name() << "\n";

    Real tol = 10*std::numeric_limits<Real>::epsilon();
    std::vector<std::array<Real, 3>> v(4);
    v[0] = {0,0,0};
    v[1] = {1,0,0};
    v[2] = {2,0,0};
    v[3] = {3,0,0};
    catmull_rom<std::array<Real, 3>> cr(std::move(v));

    // Test that the interpolation condition is obeyed:
    CHECK_ULP_CLOSE(3, cr.max_parameter(), 4);
    auto p0 = cr(0.0);
    CHECK_ABSOLUTE_ERROR(0, p0[0], tol);
    CHECK_ABSOLUTE_ERROR(0, p0[1], tol);
    CHECK_ABSOLUTE_ERROR(0, p0[2], tol);
    auto p1 = cr(1.0);
    CHECK_ULP_CLOSE(1, p1[0], 4);
    CHECK_ABSOLUTE_ERROR(0, p1[1], tol);
    CHECK_ABSOLUTE_ERROR(0, p1[2], tol);

    auto p2 = cr(2.0);
    CHECK_ULP_CLOSE(2, p2[0], 4);
    CHECK_ABSOLUTE_ERROR(0, p2[1], tol);
    CHECK_ABSOLUTE_ERROR(0, p2[2], tol);


    auto p3 = cr(3.0);
    CHECK_ULP_CLOSE(3, p3[0], 4);
    CHECK_ABSOLUTE_ERROR(0, p3[1], tol);
    CHECK_ABSOLUTE_ERROR(0, p3[2], tol);

    Real s = cr.parameter_at_point(0);
    CHECK_ABSOLUTE_ERROR(0, s, tol);

    s = cr.parameter_at_point(1);
    CHECK_ULP_CLOSE(1, s, 4);

    s = cr.parameter_at_point(2);
    CHECK_ULP_CLOSE(2, s, 4);

    s = cr.parameter_at_point(3);
    CHECK_ULP_CLOSE(3, s, 4);

    // Test that the function is linear on the interval [1,2]:
    for (double s = 1; s < 2; s += 0.01)
    {
        auto p = cr(s);
        CHECK_ULP_CLOSE(Real(s), p[0], 4);
        CHECK_ABSOLUTE_ERROR(0, p[1], tol);
        CHECK_ABSOLUTE_ERROR(0, p[2], tol);

        auto tangent = cr.prime(s);
        CHECK_ULP_CLOSE(1, tangent[0], 8);
        CHECK_ABSOLUTE_ERROR(0, tangent[1], tol);
        CHECK_ABSOLUTE_ERROR(0, tangent[2], tol);
    }

}

template<class Real>
void test_circle()
{
    using boost::math::constants::pi;
    using std::cos;
    using std::sin;

    std::cout << "Testing that the Catmull-Rom spline interpolates circles correctly on type "
              << boost::typeindex::type_id<Real>().pretty_name() << "\n";

    Real tol = 10*std::numeric_limits<Real>::epsilon();
    std::vector<std::array<Real, 2>> v(20*sizeof(Real));
    std::vector<std::array<Real, 2>> u(20*sizeof(Real));
    for (size_t i = 0; i < v.size(); ++i)
    {
        Real theta = ((Real) i/ (Real) v.size())*2*pi<Real>();
        v[i] = {cos(theta), sin(theta)};
        u[i] = v[i];
    }
    catmull_rom<std::array<Real, 2>> circle(std::move(v), true);

    // Interpolation condition (v was moved from, so iterate over its copy u):
    for (size_t i = 0; i < u.size(); ++i)
    {
        Real s = circle.parameter_at_point(i);
        auto p = circle(s);
        Real x = p[0];
        Real y = p[1];
        if (abs(x) < std::numeric_limits<Real>::epsilon())
        {
            CHECK_ABSOLUTE_ERROR(0, u[i][0], tol);
        }
        if (abs(y) < std::numeric_limits<Real>::epsilon())
        {
            CHECK_ABSOLUTE_ERROR(0, u[i][1], tol);
        }
        else
        {
            CHECK_ULP_CLOSE(u[i][0], x, 12);
            CHECK_ULP_CLOSE(u[i][1], y, 12);
        }
    }

    Real max_s = circle.max_parameter();
    for(Real s = 0; s < max_s; s += Real(0.01))
    {
        auto p = circle(s);
        Real x = p[0];
        Real y = p[1];
        CHECK_ABSOLUTE_ERROR(1, x*x+y*y, Real(0.001));
    }
}


template<class Real, size_t dimension>
void test_affine_invariance()
{
    std::cout << "Testing that the Catmull-Rom spline is affine invariant in dimension "
              << dimension << " on type "
              << boost::typeindex::type_id<Real>().pretty_name() << "\n";

    Real tol = 1000*std::numeric_limits<Real>::epsilon();
    std::vector<std::array<Real, dimension>> v(100);
    std::vector<std::array<Real, dimension>> u(100);
    std::mt19937_64 gen(438232);
    Real inv_denom = (Real) 100/( (Real) (gen.max)() + (Real) 2);
    for(size_t j = 0; j < dimension; ++j)
    {
        v[0][j] = gen()*inv_denom;
        u[0][j] = v[0][j];
    }

    for (size_t i = 1; i < v.size(); ++i)
    {
        for(size_t j = 0; j < dimension; ++j)
        {
            v[i][j] = v[i-1][j] + gen()*inv_denom;
            u[i][j] = v[i][j];
        }
    }
    std::array<Real, dimension> affine_shift;
    for (size_t j = 0; j < dimension; ++j)
    {
        affine_shift[j] = gen()*inv_denom;
    }

    catmull_rom<std::array<Real, dimension>> cr1(std::move(v));

    for(size_t i = 0; i< u.size(); ++i)
    {
        for(size_t j = 0; j < dimension; ++j)
        {
            u[i][j] += affine_shift[j];
        }
    }

    catmull_rom<std::array<Real, dimension>> cr2(std::move(u));

    CHECK_ULP_CLOSE(cr2.max_parameter(), cr1.max_parameter(), 4);

    Real ds = cr1.max_parameter()/1024;
    for (Real s = 0; s < cr1.max_parameter(); s += ds)
    {
        auto p0 = cr1(s);
        auto p1 = cr2(s);
        auto tangent0 = cr1.prime(s);
        auto tangent1 = cr2.prime(s);
        for (size_t j = 0; j < dimension; ++j)
        {
            CHECK_ULP_CLOSE(p1[j], p0[j] + affine_shift[j], 28);
            if (abs(tangent0[j]) > 5000*tol)
            {
                CHECK_ULP_CLOSE(tangent1[j], tangent0[j], 1600868);
            }
        }
    }
}

template<class Real>
void test_helix()
{
    using boost::math::constants::pi;
    std::cout << "Testing that the Catmull-Rom spline interpolates helices correctly on type "
              << boost::typeindex::type_id<Real>().pretty_name() << "\n";

    Real tol = 0.001;
    std::vector<std::array<Real, 3>> v(400*sizeof(Real));
    for (size_t i = 0; i < v.size(); ++i)
    {
        Real theta = ((Real) i/ (Real) v.size())*2*pi<Real>();
        v[i] = {cos(theta), sin(theta), theta};
    }
    const size_t num_points = v.size();
    catmull_rom<std::array<Real, 3>> helix(std::move(v));

    // Interpolation condition (v was moved from, so use the saved size):
    for (size_t i = 0; i < num_points; ++i)
    {
        Real s = helix.parameter_at_point(i);
        auto p = helix(s);
        Real t = p[2];

        Real x = p[0];
        Real y = p[1];
        if (abs(x) < tol)
        {
            CHECK_ABSOLUTE_ERROR(0, Real(cos(t)), tol);
        }
        else
        {
            CHECK_ABSOLUTE_ERROR(cos(t), x, Real(tol*abs(cos(t))));
        }
        if (abs(y) < tol)
        {
            CHECK_ABSOLUTE_ERROR(0, Real(sin(t)), tol);
        }
        else
        {
            CHECK_ABSOLUTE_ERROR(sin(t), y, Real(tol*abs(sin(t))));
        }
    }

    Real max_s = helix.max_parameter();
    for(Real s = helix.parameter_at_point(1); s < max_s; s += 0.01)
    {
        auto p = helix(s);
        Real x = p[0];
        Real y = p[1];
        Real t = p[2];
        CHECK_ABSOLUTE_ERROR(Real(1), x*x+y*y, Real(0.01));
        if (abs(x) < 0.01)
        {
            CHECK_ABSOLUTE_ERROR(0, Real(cos(t)), (Real) 0.05);
        }
        if (abs(y) < 0.01)
        {
            CHECK_ABSOLUTE_ERROR(0, Real(sin(t)), (Real) 0.05);
        }
        else
        {
            CHECK_ABSOLUTE_ERROR(Real(cos(t)), x, Real(((Real) 0.05)*abs(Real(cos(t)))));
            CHECK_ABSOLUTE_ERROR(Real(sin(t)), y, Real(((Real) 0.05)*abs(Real(sin(t)))));
        }
    }
}


template<class Real>
class mypoint3d
{
public:
    // Must define a value_type:
    typedef Real value_type;

    // Regular constructor:
    mypoint3d(Real x, Real y, Real z)
    {
        m_vec[0] = x;
        m_vec[1] = y;
        m_vec[2] = z;
    }

    // Must define a default constructor:
    mypoint3d() {}

    // Must define array access:
    Real operator[](size_t i) const
    {
        return m_vec[i];
    }

    // Array element assignment:
    Real& operator[](size_t i)
    {
        return m_vec[i];
    }


private:
    std::array<Real, 3>  m_vec;
};


// Must define the free function "size()":
template<class Real>
constexpr std::size_t size(const mypoint3d<Real>&)
{
    return 3;
}

template<class Real>
void test_data_representations()
{
    std::cout << "Testing that the Catmull-Rom spline works with multiple data representations.\n";
    mypoint3d<Real> p0(Real(0.1), Real(0.2), Real(0.3));
    mypoint3d<Real> p1(Real(0.2), Real(0.3), Real(0.4));
    mypoint3d<Real> p2(Real(0.3), Real(0.4), Real(0.5));
    mypoint3d<Real> p3(Real(0.4), Real(0.5), Real(0.6));
    mypoint3d<Real> p4(Real(0.5), Real(0.6), Real(0.7));
    mypoint3d<Real> p5(Real(0.6), Real(0.7), Real(0.8));


    // Tests initializer_list:
    catmull_rom<mypoint3d<Real>> cat({p0, p1, p2, p3, p4, p5});

    Real tol = Real(0.001);
    auto p = cat(cat.parameter_at_point(0));
    CHECK_ABSOLUTE_ERROR(p0[0], p[0], Real(tol*abs(p0[0])));
    CHECK_ABSOLUTE_ERROR(p0[1], p[1], Real(tol*abs(p0[1])));
    CHECK_ABSOLUTE_ERROR(p0[2], p[2], Real(tol*abs(p0[2])));
    p = cat(cat.parameter_at_point(1));
    CHECK_ABSOLUTE_ERROR(p1[0], p[0], Real(tol*abs(p1[0])));
    CHECK_ABSOLUTE_ERROR(p1[1], p[1], Real(tol*abs(p1[1])));
    CHECK_ABSOLUTE_ERROR(p1[2], p[2], Real(tol*abs(p1[2])));
}

template<class Real>
void test_random_access_container()
{
    std::cout << "Testing that the Catmull-Rom spline works with multiple data representations.\n";
    mypoint3d<Real> p0(0.1, 0.2, 0.3);
    mypoint3d<Real> p1(0.2, 0.3, 0.4);
    mypoint3d<Real> p2(0.3, 0.4, 0.5);
    mypoint3d<Real> p3(0.4, 0.5, 0.6);
    mypoint3d<Real> p4(0.5, 0.6, 0.7);
    mypoint3d<Real> p5(0.6, 0.7, 0.8);

    boost::numeric::ublas::vector<mypoint3d<Real>> u(6);
    u[0] = p0;
    u[1] = p1;
    u[2] = p2;
    u[3] = p3;
    u[4] = p4;
    u[5] = p5;

    // Tests initializer_list:
    catmull_rom<mypoint3d<Real>, decltype(u)> cat(std::move(u));

    Real tol = Real(0.001);
    auto p = cat(cat.parameter_at_point(0));
    CHECK_ABSOLUTE_ERROR(p0[0], p[0], Real(tol*abs(p0[0])));
    CHECK_ABSOLUTE_ERROR(p0[1], p[1], Real(tol*abs(p0[1])));
    CHECK_ABSOLUTE_ERROR(p0[2], p[2], Real(tol*abs(p0[2])));
    p = cat(cat.parameter_at_point(1));
    CHECK_ABSOLUTE_ERROR(p1[0], p[0], Real(tol*abs(p1[0])));
    CHECK_ABSOLUTE_ERROR(p1[1], p[1], Real(tol*abs(p1[1])));
    CHECK_ABSOLUTE_ERROR(p1[2], p[2], Real(tol*abs(p1[2])));
}

void catmull_rom_test()
{
#if !defined(TEST) || (TEST == 1)

    #ifdef __STDCPP_FLOAT32_T__
    test_circle<std::float32_t>();
    test_data_representations<std::float32_t>();
    #else
    test_circle<float>();
    test_data_representations<float>();
    #endif

    #ifdef __STDCPP_FLOAT64_T__
    test_alpha_distance<std::float64_t>();
    test_linear<std::float64_t>();
    test_circle<std::float64_t>();
    #else
    test_alpha_distance<double>();
    test_linear<double>();
    test_circle<double>();
    #endif
    
    #ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
    test_linear<long double>();
    #endif

#endif

#if !defined(TEST) || (TEST == 2)

    #ifdef __STDCPP_FLOAT64_T__
    test_helix<std::float64_t>();
    test_affine_invariance<std::float64_t, 1>();
    test_affine_invariance<std::float64_t, 2>();
    test_affine_invariance<std::float64_t, 3>();
    test_affine_invariance<std::float64_t, 4>();
    test_random_access_container<std::float64_t>();
    #else
    test_helix<double>();
    test_affine_invariance<double, 1>();
    test_affine_invariance<double, 2>();
    test_affine_invariance<double, 3>();
    test_affine_invariance<double, 4>();
    test_random_access_container<double>();
    #endif

#endif

#if !defined(TEST) || (TEST == 3)
    test_affine_invariance<cpp_bin_float_50, 4>();
#endif
}

int main()
{
    catmull_rom_test();
    return boost::math::test::report_errors();
}
