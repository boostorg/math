#define BOOST_TEST_MODULE test_cubic_b_spline

#include <random>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/type_index.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/math/interpolators/cubic_b_spline.hpp>



using boost::math::constants::third;
using boost::math::constants::half;

template<class Real>
void test_b3_spline()
{
    std::cout << "Testing evaluation of spline basis functions on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    // Outside the support:
    auto eps = std::numeric_limits<Real>::epsilon();
    BOOST_CHECK_SMALL(boost::math::b3_spline<Real>(2.5), (Real) 0);
    BOOST_CHECK_SMALL(boost::math::b3_spline<Real>(-2.5), (Real) 0);
    BOOST_CHECK_SMALL(boost::math::b3_spline_prime<Real>(2.5), (Real) 0);
    BOOST_CHECK_SMALL(boost::math::b3_spline_prime<Real>(-2.5), (Real) 0);

    // On the boundary of support:
    BOOST_CHECK_SMALL(boost::math::b3_spline<Real>(2), (Real) 0);
    BOOST_CHECK_SMALL(boost::math::b3_spline<Real>(-2), (Real) 0);
    BOOST_CHECK_SMALL(boost::math::b3_spline_prime<Real>(2), (Real) 0);
    BOOST_CHECK_SMALL(boost::math::b3_spline_prime<Real>(-2), (Real) 0);

    // Special values:
    BOOST_CHECK_CLOSE(boost::math::b3_spline<Real>(-1), third<Real>()*half<Real>(), eps);
    BOOST_CHECK_CLOSE(boost::math::b3_spline<Real>( 1), third<Real>()*half<Real>(), eps);
    BOOST_CHECK_CLOSE(boost::math::b3_spline<Real>(0), 2*third<Real>(), eps);

    BOOST_CHECK_CLOSE(boost::math::b3_spline_prime<Real>(-1), half<Real>(), eps);
    BOOST_CHECK_CLOSE(boost::math::b3_spline_prime<Real>( 1), -half<Real>(), eps);
    BOOST_CHECK_SMALL(boost::math::b3_spline_prime<Real>(0), eps);

    // Properties: B3 is an even function, B3' is an odd function.
    for (size_t i = 1; i < 200; ++i)
    {
        auto arg = i*0.01;
        BOOST_CHECK_CLOSE(boost::math::b3_spline<Real>(arg), boost::math::b3_spline<Real>(arg), eps);
        BOOST_CHECK_CLOSE(boost::math::b3_spline_prime<Real>(-arg), -boost::math::b3_spline_prime<Real>(arg), eps);
    }

}
/*
 * This test ensures that the interpolant s(x_j) = f(x_j) at all grid points.
 */
template<class Real>
void test_interpolation_condition()
{
    std::cout << "Testing interpolation condition for cubic b splines on type " << boost::typeindex::type_id<Real>().pretty_name()  << "\n";
    std::random_device rd;
    std::mt19937 gen(rd());
    boost::random::uniform_real_distribution<Real> dis(1, 10);
    boost::random::uniform_real_distribution<Real> step_distribution(0.001, 0.01);
    std::vector<Real> v(5000);
    for (size_t i = 0; i < v.size(); ++i)
    {
        v[i] = dis(gen);
    }

    Real step = step_distribution(gen);
    Real a = dis(gen);
    boost::math::cubic_b_spline<Real> spline(v.data(), v.size(), a, step);

    for (size_t i = 0; i < v.size(); ++i)
    {
        auto y = spline.interpolate_at(i*step + a);
        // This seems like a very large tolerance, but I don't know of any other interpolators
        // that will be able to do much better on random data.
        BOOST_CHECK_CLOSE(y, v[i], 10000000*std::numeric_limits<Real>::epsilon());
    }

}

template<class Real>
void test_constant_function()
{
    std::cout << "Testing that constants are interpolated correctly by cubic b splines on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    std::random_device rd;
    std::mt19937 gen(rd());
    boost::random::uniform_real_distribution<Real> dis(-100, 100);
    boost::random::uniform_real_distribution<Real> step_distribution(0.0001, 0.1);
    std::vector<Real> v(500);
    Real constant = dis(gen);
    for (size_t i = 0; i < v.size(); ++i)
    {
        v[i] = constant;
    }

    Real step = step_distribution(gen);
    Real a = dis(gen);
    boost::math::cubic_b_spline<Real> spline(v.data(), v.size(), a, step);

    for (size_t i = 0; i < v.size(); ++i)
    {
        auto y = spline.interpolate_at(i*step + a);
        BOOST_CHECK_CLOSE(y, v[i], 10*std::numeric_limits<Real>::epsilon());
        auto y_prime = spline.interpolate_derivative(i*step + a);
        BOOST_CHECK_SMALL(y_prime, 5000*std::numeric_limits<Real>::epsilon());
    }

    // Test that correctly specified left and right-derivatives work properly:
    spline = boost::math::cubic_b_spline<Real>(v.data(), v.size(), a, step, 0, 0);

    for (size_t i = 0; i < v.size(); ++i)
    {
        auto y = spline.interpolate_at(i*step + a);
        BOOST_CHECK_CLOSE(y, v[i], std::numeric_limits<Real>::epsilon());
        auto y_prime = spline.interpolate_derivative(i*step + a);
        BOOST_CHECK_SMALL(y_prime, std::numeric_limits<Real>::epsilon());
    }
}


template<class Real>
void test_affine_function()
{
    std::cout << "Testing that affine functions are interpolated correctly by cubic b splines on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    std::random_device rd;
    std::mt19937 gen(rd());
    boost::random::uniform_real_distribution<Real> dis(-100, 100);
    boost::random::uniform_real_distribution<Real> step_distribution(0.0001, 0.01);
    std::vector<Real> v(500);
    Real a = dis(gen);
    Real b = dis(gen);
    Real step = step_distribution(gen);
    Real x0 = dis(gen);

    for (size_t i = 0; i < v.size(); ++i)
    {
        v[i] = a*(i*step - x0) + b;
    }


    boost::math::cubic_b_spline<Real> spline(v.data(), v.size(), x0, step);

    for (size_t i = 0; i < v.size(); ++i)
    {
        auto y = spline.interpolate_at(i*step + x0);
        BOOST_CHECK_CLOSE(y, v[i], 10000*std::numeric_limits<Real>::epsilon());
        auto y_prime = spline.interpolate_derivative(i*step + x0);
        BOOST_CHECK_CLOSE(y_prime, a, 10000000*std::numeric_limits<Real>::epsilon());
    }

    // Test that correctly specified left and right-derivatives work properly:
    spline = boost::math::cubic_b_spline<Real>(v.data(), v.size(), x0, step, a, a);

    for (size_t i = 0; i < v.size(); ++i)
    {
        auto y = spline.interpolate_at(i*step + x0);
        BOOST_CHECK_CLOSE(y, v[i], 10000*std::numeric_limits<Real>::epsilon());
        auto y_prime = spline.interpolate_derivative(i*step + x0);
        BOOST_CHECK_CLOSE(y_prime, a, 10000000*std::numeric_limits<Real>::epsilon());
    }
}

template<class Real>
void test_quadratic_function()
{
    std::cout << "Testing that quadratic functions are interpolated correctly by cubic b splines on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    std::random_device rd;
    std::mt19937 gen(rd());
    boost::random::uniform_real_distribution<Real> dis(1, 5);
    boost::random::uniform_real_distribution<Real> step_distribution(0.0001, 1);
    std::vector<Real> v(500);
    Real a = dis(gen);
    Real b = dis(gen);
    Real c = dis(gen);
    Real step = step_distribution(gen);
    Real x0 = dis(gen);

    for (size_t i = 0; i < v.size(); ++i)
    {
        v[i] = a*(i*step - x0)*(i*step - x0) + b*(i*step - x0) + c;
    }

    boost::math::cubic_b_spline<Real> spline(v.data(), v.size(), x0, step);

    for (size_t i = 0; i < v.size(); ++i)
    {
        auto y = spline.interpolate_at(i*step + x0);
        BOOST_CHECK_CLOSE(y, v[i], 100000000*std::numeric_limits<Real>::epsilon());
        auto y_prime = spline.interpolate_derivative(i*step + x0);
        BOOST_CHECK_CLOSE(y_prime, 2*a*(i*step-x0) + b, 2.0);
    }
}


BOOST_AUTO_TEST_CASE(test_cubic_b_spline)
{
    test_b3_spline<float>();
    test_b3_spline<double>();
    test_b3_spline<long double>();

    test_interpolation_condition<float>();
    test_interpolation_condition<double>();
    test_interpolation_condition<long double>();

    test_constant_function<float>();
    test_constant_function<double>();
    test_constant_function<long double>();

    test_affine_function<float>();
    test_affine_function<double>();
    test_affine_function<long double>();

    test_quadratic_function<float>();
    test_quadratic_function<double>();
    test_quadratic_function<long double>();

}
