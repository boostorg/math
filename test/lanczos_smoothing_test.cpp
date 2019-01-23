/*
 * Copyright Nick Thompson, 2019
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
#define BOOST_TEST_MODULE lanczos_smoothing_test

#include <random>
#include <boost/math/constants/constants.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/math/differentiation/lanczos_smoothing.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/math/special_functions/next.hpp> // for float_distance

using std::abs;
using std::pow;
using std::sqrt;
using std::sin;
using boost::math::constants::two_pi;
using boost::multiprecision::cpp_bin_float_50;
using boost::multiprecision::cpp_bin_float_100;
using boost::math::differentiation::discrete_lanczos_derivative;
using boost::math::differentiation::detail::discrete_legendre;
using boost::math::differentiation::detail::interior_velocity_filter;
using boost::math::differentiation::detail::boundary_velocity_filter;

template<class Real>
void test_dlp_norms()
{
    std::cout << "Testing Discrete Legendre Polynomial norms on type " << typeid(Real).name() << "\n";
    Real tol = std::numeric_limits<Real>::epsilon();
    auto dlp = discrete_legendre<Real>(1);
    BOOST_CHECK_CLOSE_FRACTION(dlp.norm_sq(0), 3, tol);
    BOOST_CHECK_CLOSE_FRACTION(dlp.norm_sq(1), 2, tol);
    dlp = discrete_legendre<Real>(2);
    BOOST_CHECK_CLOSE_FRACTION(dlp.norm_sq(0), Real(5)/Real(2), tol);
    BOOST_CHECK_CLOSE_FRACTION(dlp.norm_sq(1), Real(5)/Real(4), tol);
    BOOST_CHECK_CLOSE_FRACTION(dlp.norm_sq(2), Real(3*3*7)/Real(pow(2,6)), 2*tol);
    dlp = discrete_legendre<Real>(200);
    for(size_t r = 0; r < 10; ++r)
    {
        Real calc = dlp.norm_sq(r);
        Real expected = Real(2)/Real(2*r+1);
        // As long as r << n, ||q_r||^2 -> 2/(2r+1) as n->infty
        BOOST_CHECK_CLOSE_FRACTION(calc, expected, 0.05);
    }

}

template<class Real>
void test_dlp_evaluation()
{
    std::cout << "Testing evaluation of Discrete Legendre polynomials on type " << typeid(Real).name() << "\n";
    Real tol = std::numeric_limits<Real>::epsilon();
    size_t n = 25;
    auto dlp = discrete_legendre<Real>(n);
    Real x = 0.72;
    Real q0 = dlp(x, 0);
    BOOST_TEST(q0 == 1);
    Real q1 = dlp(x, 1);
    BOOST_TEST(q1 == x);
    Real q2 = dlp(x, 2);
    int N = 2*n+1;
    Real expected = 0.5*(3*x*x - Real(N*N - 1)/Real(4*n*n));
    BOOST_CHECK_CLOSE_FRACTION(q2, expected, tol);
    Real q3 = dlp(x, 3);
    expected = (x/3)*(5*expected - (Real(N*N - 4))/(2*n*n));
    BOOST_CHECK_CLOSE_FRACTION(q3, expected, 2*tol);

    // q_r(x) is even for even r, and odd for odd r:
    for (size_t n = 8; n < 22; ++n)
    {
        dlp = discrete_legendre<Real>(n);
        for(size_t r = 2; r <= n; ++r)
        {
            if (r & 1)
            {
                Real q1 = dlp(x, r);
                Real q2 = -dlp(-x, r);
                BOOST_CHECK_CLOSE_FRACTION(q1, q2, tol);
            }
            else
            {
                Real q1 = dlp(x, r);
                Real q2 = dlp(-x, r);
                BOOST_CHECK_CLOSE_FRACTION(q1, q2, tol);
            }

            Real l2_sq = 0;
            for (int j = -(int)n; j <= (int) n; ++j)
            {
                Real y = Real(j)/Real(n);
                Real term = dlp(y, r);
                l2_sq += term*term;
            }
            l2_sq /= n;
            Real l2_sq_expected = dlp.norm_sq(r);
            BOOST_CHECK_CLOSE_FRACTION(l2_sq, l2_sq_expected, 20*tol);
        }
    }
}

template<class Real>
void test_dlp_next()
{
    std::cout << "Testing Discrete Legendre polynomial 'next' function on type " << typeid(Real).name() << "\n";
    Real tol = std::numeric_limits<Real>::epsilon();

    for(size_t n = 2; n < 20; ++n)
    {
        auto dlp = discrete_legendre<Real>(n);
        for(Real x = -1; x <= 1; x += 0.1)
        {
            dlp.initialize_recursion(x);
            for (size_t k = 2; k < n; ++k)
            {
                BOOST_CHECK_CLOSE(dlp.next(), dlp(x, k), tol);
            }

            dlp.initialize_recursion(x);
            for (size_t k = 2; k < n; ++k)
            {
                BOOST_CHECK_CLOSE(dlp.next_prime(), dlp.prime(x, k), tol);
            }
        }
    }
}


template<class Real>
void test_dlp_derivatives()
{
    std::cout << "Testing Discrete Legendre polynomial derivatives on type " << typeid(Real).name() << "\n";
    Real tol = 10*std::numeric_limits<Real>::epsilon();
    int n = 25;
    auto dlp = discrete_legendre<Real>(n);
    Real x = 0.72;
    Real q0p = dlp.prime(x, 0);
    BOOST_TEST(q0p == 0);
    Real q1p = dlp.prime(x, 1);
    BOOST_TEST(q1p == 1);
    Real q2p = dlp.prime(x, 2);
    Real expected = 3*x;
    BOOST_CHECK_CLOSE_FRACTION(q2p, expected, tol);
}

template<class Real>
void test_dlp_second_derivative()
{
    std::cout << "Testing Discrete Legendre polynomial derivatives on type " << typeid(Real).name() << "\n";
    int n = 25;
    auto dlp = discrete_legendre<Real>(n);
    Real x = Real(1)/Real(3);
    dlp.initialize_recursion(x);
    Real q2pp = dlp.next_dbl_prime();
    BOOST_TEST(q2pp == 3);
}


template<class Real>
void test_interior_velocity_filter()
{
    std::cout << "Testing interior filter on type " << typeid(Real).name() << "\n";
    Real tol = std::numeric_limits<Real>::epsilon();
    for(int n = 1; n < 10; ++n)
    {
        for (int p = 1; p < n; p += 2)
        {
            auto f = interior_velocity_filter<Real>(n,p);
            // Since we only store half the filter coefficients,
            // we need to reindex the moment sums:
            Real sum = 0;
            for (size_t j = 0; j < f.size(); ++j)
            {
                sum += j*f[j];
            }
            BOOST_CHECK_CLOSE_FRACTION(2*sum, 1, 1000*tol);

            for (int l = 3; l <= p; l += 2)
            {
                sum = 0;
                for (size_t j = 0; j < f.size(); ++j)
                {
                    // The condition number of this sum is infinite!
                    // No need to get to worked up about the tolerance.
                    sum += pow(Real(j), l)*f[j];
                }
                BOOST_CHECK_SMALL(sum, sqrt(tol)/100);
            }
            //std::cout << "(n,p) = (" << n  << "," << p << ") = {";
            //for (auto & x : f)
            //{
            //    std::cout << x << ", ";
            //}
            //std::cout << "}\n";
        }
    }
}

template<class Real>
void test_interior_lanczos()
{
    std::cout << "Testing interior Lanczos on type " << typeid(Real).name() << "\n";
    Real tol = std::numeric_limits<Real>::epsilon();
    std::vector<Real> v(500);
    std::fill(v.begin(), v.end(), 7);

    for (size_t n = 1; n < 10; ++n)
    {
        for (size_t p = 2; p < 2*n; p += 2)
        {
            auto dld = discrete_lanczos_derivative(Real(0.1), n, p);
            for (size_t m = n; m < v.size() - n; ++m)
            {
                Real dvdt = dld(v, m);
                BOOST_CHECK_SMALL(dvdt, tol);
            }
            auto dvdt = dld(v);
            for (size_t m = n; m < v.size() - n; ++m)
            {
                BOOST_CHECK_SMALL(dvdt[m], tol);
            }
        }
    }


    for(size_t i = 0; i < v.size(); ++i)
    {
        v[i] = 7*i+8;
    }

    for (size_t n = 1; n < 10; ++n)
    {
        for (size_t p = 2; p < 2*n; p += 2)
        {
            auto dld = discrete_lanczos_derivative(Real(1), n, p);
            for (size_t m = n; m < v.size() - n; ++m)
            {
                Real dvdt = dld(v, m);
                BOOST_CHECK_CLOSE_FRACTION(dvdt, 7, 2000*tol);
            }
            auto dvdt = dld(v);
            for (size_t m = n; m < v.size() - n; ++m)
            {
                BOOST_CHECK_CLOSE_FRACTION(dvdt[m], 7, 2000*tol);
            }
        }
    }

    //std::random_device rd{};
    //auto seed = rd();
    //std::cout << "Seed = " << seed << "\n";
    std::mt19937 gen(4172378669);
    std::normal_distribution<> dis{0, 0.01};
    for (size_t i = 0; i < v.size(); ++i)
    {
        v[i] = 7*i+8 + dis(gen);
    }

    for (size_t n = 1; n < 10; ++n)
    {
        for (size_t p = 2; p < 2*n; p += 2)
        {
            auto dld = discrete_lanczos_derivative(Real(1), n, p);
            for (size_t m = n; m < v.size() - n; ++m)
            {
                BOOST_CHECK_CLOSE_FRACTION(dld(v, m), Real(7), Real(0.0042));
            }
        }
    }


    for (size_t i = 0; i < v.size(); ++i)
    {
        v[i] = 15*i*i + 7*i+8 + dis(gen);
    }

    for (size_t n = 1; n < 10; ++n)
    {
        for (size_t p = 2; p < 2*n; p += 2)
        {
            auto dld = discrete_lanczos_derivative(Real(1), n, p);
            for (size_t m = n; m < v.size() - n; ++m)
            {
                BOOST_CHECK_CLOSE_FRACTION(dld(v,m), Real(30*m + 7), Real(0.00008));
            }
        }
    }

    std::normal_distribution<> dis1{0, 0.0001};
    Real omega = Real(1)/Real(16);
    for (size_t i = 0; i < v.size(); ++i)
    {
        v[i] = sin(i*omega) + dis1(gen);
    }

    for (size_t n = 10; n < 20; ++n)
    {
        for (size_t p = 3; p < 100 && p < n/2; p += 2)
        {
            auto dld = discrete_lanczos_derivative(Real(1), n, p);

            for (size_t m = n; m < v.size() - n && m < n + 10; ++m)
            {
                BOOST_CHECK_CLOSE_FRACTION(dld(v,m), omega*cos(omega*m), Real(0.03));
            }
        }
    }
}

template<class Real>
void test_boundary_velocity_filters()
{
    std::cout << "Testing boundary filters on type " << typeid(Real).name() << "\n";
    Real tol = std::numeric_limits<Real>::epsilon();
    for(int n = 1; n < 5; ++n)
    {
        for (int p = 1; p < 2*n+1; ++p)
        {
            for (int s = -n; s <= n; ++s)
            {
                auto f = boundary_velocity_filter<Real>(n, p, s);
                // Sum is zero:
                Real sum = 0;
                Real c = 0;
                for (auto & x : f)
                {
                    Real y = x - c;
                    Real t = sum + y;
                    c = (t-sum) -y;
                    sum = t;
                }
                BOOST_CHECK_SMALL(sum, 200*tol);

                sum = 0;
                c = 0;
                for (size_t k = 0; k < f.size(); ++k)
                {
                    Real j = Real(k) - Real(n);
                    // note the shifted index here:
                    Real x = (j-s)*f[k];
                    Real y = x - c;
                    Real t = sum + y;
                    c = (t-sum) -y;
                    sum = t;
                }
                BOOST_CHECK_CLOSE_FRACTION(sum, 1, 350*tol);


                for (int l = 2; l <= p; ++l)
                {
                    sum = 0;
                    c = 0;
                    for (size_t k = 0; k < f.size(); ++k)
                    {
                        Real j = Real(k) - Real(n);
                        // The condition number of this sum is infinite!
                        // No need to get to worked up about the tolerance.
                        Real x = pow(j-s, l)*f[k];
                        Real y = x - c;
                        Real t = sum + y;
                        c = (t-sum) -y;
                        sum = t;
                    }
                    BOOST_CHECK_SMALL(sum, sqrt(tol)/10);
                }

                //std::cout << "(n,p,s) = ("<< n << ", " << p << "," << s << ") = {";
                //for (auto & x : f)
                //{
                //    std::cout << x << ", ";
                //}
                //std::cout << "}\n";*/
            }
        }
    }
}

template<class Real>
void test_boundary_lanczos()
{
    std::cout << "Testing Lanczos boundary on type " << typeid(Real).name() << "\n";
    Real tol = std::numeric_limits<Real>::epsilon();
    std::vector<Real> v(500, 7);

    for (size_t n = 1; n < 10; ++n)
    {
        for (size_t p = 2; p < 2*n; ++p)
        {
            auto lsd = discrete_lanczos_derivative(Real(0.0125), n, p);
            for (size_t m = 0; m < n; ++m)
            {
                Real dvdt = lsd(v,m);
                BOOST_CHECK_SMALL(dvdt, 4*sqrt(tol));
            }
            for (size_t m = v.size() - n; m < v.size(); ++m)
            {
                Real dvdt = lsd(v,m);
                BOOST_CHECK_SMALL(dvdt, 4*sqrt(tol));
            }
        }
    }

    for(size_t i = 0; i < v.size(); ++i)
    {
        v[i] = 7*i+8;
    }

    for (size_t n = 3; n < 10; ++n)
    {
        for (size_t p = 2; p < 2*n; ++p)
        {
            auto lsd = discrete_lanczos_derivative(Real(1), n, p);
            for (size_t m = 0; m < n; ++m)
            {
                Real dvdt = lsd(v,m);
                BOOST_CHECK_CLOSE_FRACTION(dvdt, 7, sqrt(tol));
            }

            for (size_t m = v.size() - n; m < v.size(); ++m)
            {
                Real dvdt = lsd(v,m);
                BOOST_CHECK_CLOSE_FRACTION(dvdt, 7, 4*sqrt(tol));
            }
        }
    }

    for (size_t i = 0; i < v.size(); ++i)
    {
        v[i] = 15*i*i + 7*i+8;
    }

    for (size_t n = 1; n < 10; ++n)
    {
        for (size_t p = 2; p < 2*n; ++p)
        {
            auto lsd = discrete_lanczos_derivative(Real(1), n, p);
            for (size_t m = 0; m < v.size(); ++m)
            {
                BOOST_CHECK_CLOSE_FRACTION(lsd(v,m), 30*m+7, 30*sqrt(tol));
            }
        }
    }

    // Demonstrate that the boundary filters are also denoising:
    //std::random_device rd{};
    //auto seed = rd();
    //std::cout << "seed = " << seed << "\n";
    std::mt19937 gen(311354333);
    std::normal_distribution<> dis{0, 0.01};
    for (size_t i = 0; i < v.size(); ++i)
    {
        v[i] += dis(gen);
    }

    for (size_t n = 1; n < 10; ++n)
    {
        for (size_t p = 2; p < n; ++p)
        {
            auto lsd = discrete_lanczos_derivative(Real(1), n, p);
            for (size_t m = 0; m < v.size(); ++m)
            {
                BOOST_CHECK_CLOSE_FRACTION(lsd(v,m), 30*m+7, 0.005);
            }
            auto dvdt = lsd(v);
            for (size_t m = 0; m < v.size(); ++m)
            {
                BOOST_CHECK_CLOSE_FRACTION(dvdt[m], 30*m+7, 0.005);
            }
        }
    }
}

template<class Real>
void test_acceleration_filters()
{
    Real eps = std::numeric_limits<Real>::epsilon();
    for (size_t n = 1; n < 5; ++n)
    {
        for(size_t p = 3; p <= 2*n; ++p)
        {
            for(int64_t s = -int64_t(n); s <= 0; ++s)
            {
                auto g = boost::math::differentiation::detail::acceleration_filter<long double>(n,p,s);

                std::vector<Real> f(g.size());
                for (size_t i = 0; i < g.size(); ++i)
                {
                    f[i] = static_cast<Real>(g[i]);
                }

                Real sum = 0;
                Real c = 0;
                for (auto & x : f)
                {
                    Real y = x - c;
                    Real t = sum + y;
                    c = (t - sum) - y;
                    sum = t;
                }
                BOOST_CHECK_SMALL(abs(sum), 1000*eps);

                sum = 0;
                for (size_t k = 0; k < f.size(); ++k)
                {
                    Real j = Real(k) - Real(n);
                    sum += (j-s)*f[k];
                }
                BOOST_CHECK_SMALL(sum, sqrt(eps));

                sum = 0;
                for (size_t k = 0; k < f.size(); ++k)
                {
                    Real j = Real(k) - Real(n);
                    sum += (j-s)*(j-s)*f[k];
                }
                BOOST_CHECK_CLOSE_FRACTION(sum, 2, 4*sqrt(eps));
                // See unlabelled equation in McDevitt, 2012, just after equation 26:
                // It appears that there is an off-by-one error in that equation, since p + 1 moments don't vanish, only p.
                // This test is itself suspect; the condition number of the moment sum is infinite.
                // So the *slightest* error in the filter gets amplified by the test; in terms of the
                // behavior of the actual filter, it's not a big deal.
                for (size_t l = 3; l <= p; ++l)
                {
                    sum = 0;
                    Real c = 0;
                    for (size_t k = 0; k < f.size(); ++k)
                    {
                        Real j = Real(k) - Real(n);
                        Real term = pow((j-s), l)*f[k];
                        Real y = term - c;
                        Real t = sum + y;
                        c = (t - sum) - y;
                        sum = t;
                    }
                    BOOST_CHECK_SMALL(abs(sum), 500*sqrt(eps));
                }
            }
        }
    }
}

template<class Real>
void test_lanczos_acceleration()
{
    Real eps = std::numeric_limits<Real>::epsilon();
    std::vector<Real> v(100, 7);
    auto lanczos = discrete_lanczos_derivative<Real, 2>(Real(1), 4, 3);
    for (size_t i = 0; i < v.size(); ++i)
    {
        BOOST_CHECK_SMALL(lanczos(v, i), eps);
    }

    for(size_t i = 0; i < v.size(); ++i)
    {
        v[i] = 7*i + 6;
    }
    for (size_t i = 0; i < v.size(); ++i)
    {
        BOOST_CHECK_SMALL(lanczos(v,i), 200*eps);
    }

    for(size_t i = 0; i < v.size(); ++i)
    {
        v[i] = 7*i*i + 9*i + 6;
    }
    for (size_t i = 0; i < v.size(); ++i)
    {
        BOOST_CHECK_CLOSE_FRACTION(lanczos(v, i), 14, 1500*eps);
    }

    // Now add noise, and kick up the smoothing of the Lanzcos derivative (increase n):
    //std::random_device rd{};
    //auto seed = rd();
    //std::cout << "seed = " << seed << "\n";
    size_t seed = 2507134629;
    std::mt19937 gen(seed);
    Real std_dev = 0.1;
    std::normal_distribution<Real> dis{0, std_dev};
    for (size_t i = 0; i < v.size(); ++i)
    {
        v[i] += dis(gen);
    }
    lanczos = discrete_lanczos_derivative<Real, 2>(Real(1), 18, 3);
    auto w = lanczos(v);
    for (size_t i = 0; i < v.size(); ++i)
    {
        BOOST_CHECK_CLOSE_FRACTION(w[i], 14, std_dev/200);
    }
}

BOOST_AUTO_TEST_CASE(lanczos_smoothing_test)
{
    test_dlp_second_derivative<double>();
    test_dlp_norms<double>();
    test_dlp_evaluation<double>();
    test_dlp_derivatives<double>();
    test_dlp_next<double>();
    test_dlp_norms<cpp_bin_float_50>();
    test_boundary_velocity_filters<double>();
    test_boundary_velocity_filters<long double>();
    test_boundary_velocity_filters<cpp_bin_float_50>();
    test_boundary_lanczos<double>();
    test_boundary_lanczos<long double>();
    // Takes too long!
    //test_boundary_lanczos<cpp_bin_float_50>();

    test_interior_velocity_filter<double>();
    test_interior_velocity_filter<long double>();
    test_interior_lanczos<double>();

    test_acceleration_filters<double>();

    test_lanczos_acceleration<float>();
    test_lanczos_acceleration<double>();
}
