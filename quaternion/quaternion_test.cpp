// test file for quaternion.hpp

//  (C) Copyright Hubert Holin 2001. Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.


#define    BOOST_INTERACTIVE_TEST_INPUT_ITERATOR    0

#include <iostream>
#include <sstream>

#include <boost/math/quaternion.hpp>

#include <boost/config.hpp>
#include <boost/cstdlib.hpp>  // for exit_success
#ifdef BOOST_NO_STDC_NAMESPACE
namespace std {
    using ::sqrt;
    using ::atan;
    using ::log;
    using ::exp;
    using ::cos;
    using ::sin;
    using ::tan;
    using ::cosh;
    using ::sinh;
    using ::tanh;
}
#endif

// explicit (if ludicrous) instanciation
#ifndef __GNUC__
template    class ::boost::math::quaternion<int>;
#else
// gcc doesn't like the absolutely-qualified namespace
template class boost::math::quaternion<int>;
#endif


#define BOOST_INCLUDE_MAIN  // for testing, include rather than link
#include <boost/test/test_tools.hpp>


int    test_main(int, char *[])

{
    // tests for evaluation by humans
    
    
    // using default constructor
    ::boost::math::quaternion<float>        q0;
    
    ::boost::math::quaternion<float>        qa[2];
    
    // using constructor "H seen as R^4"
    ::boost::math::quaternion<double>       q1(1,2,3,4);
    
    ::std::complex<float>                   c0(5,6);
    
    // using constructor "H seen as C^2"
    ::boost::math::quaternion<float>        q2(c0);
    
    // using UNtemplated copy constructor
    ::boost::math::quaternion<float>        q3(q2);
    
    // using templated copy constructor
    ::boost::math::quaternion<long double>  q4(q3);
    
    // using UNtemplated assignment operator
    q3 = q0;
    qa[0] = q0;
    
    // using templated assignment operator
    q4 = q0;
    qa[1] = q1;
    
    float                                   f0(7);
    
    // using converting assignment operator
    q2 = f0;
    
    // using converting assignment operator
    q3 = c0;
    
    // using += (const T &)
    q2 += f0;
    
    // using += (const ::std::complex<T> &)
    q2 += c0;
    
    // using += (const quaternion<X> &)
    q2 += q3;
    
    // using -= (const T &)
    q3 -= f0;
    
    // using -= (const ::std::complex<T> &)
    q3 -= c0;
    
    // using -= (const quaternion<X> &)
    q3 -= q2;
    
    double                                  d0(8);
    ::std::complex<double>                  c1(9,10);
    
    // using *= (const T &)
    q1 *= d0;
    
    // using *= (const ::std::complex<T> &)
    q1 *= c1;
    
    // using *= (const quaternion<X> &)
    q1 *= q1;
    
    long double                             l0(11);
    ::std::complex<long double>             c2(12,13);
    
    // using /= (const T &)
    q4 /= l0;
    
    // using /= (const ::std::complex<T> &)
    q4 /= c2;
    
    // using /= (const quaternion<X> &)
    q4 /= q1;
    
    // using + (const T &, const quaternion<T> &)
    ::boost::math::quaternion<float>        q5 = f0+q2;
    
    // using + (const quaternion<T> &, const T &)
    ::boost::math::quaternion<float>        q6 = q2+f0;
    
    // using + (const ::std::complex<T> &, const quaternion<T> &)
    ::boost::math::quaternion<float>        q7 = c0+q2;
    
    // using + (const quaternion<T> &, const ::std::complex<T> &)
    ::boost::math::quaternion<float>        q8 = q2+c0;
    
    // using + (const quaternion<T> &,const quaternion<T> &)
    ::boost::math::quaternion<float>        q9 = q2+q3;
    
    // using - (const T &, const quaternion<T> &)
    q5 = f0-q2;
    
    // using - (const quaternion<T> &, const T &)
    q6 = q2-f0;
    
    // using - (const ::std::complex<T> &, const quaternion<T> &)
    q7 = c0-q2;
    
    // using - (const quaternion<T> &, const ::std::complex<T> &)
    q8 = q2-c0;
    
    // using - (const quaternion<T> &,const quaternion<T> &)
    q9 = q2-q3;
    
    // using * (const T &, const quaternion<T> &)
    q5 = f0*q2;
    
    // using * (const quaternion<T> &, const T &)
    q6 = q2*f0;
    
    // using * (const ::std::complex<T> &, const quaternion<T> &)
    q7 = c0*q2;
    
    // using * (const quaternion<T> &, const ::std::complex<T> &)
    q8 = q2*c0;
    
    // using * (const quaternion<T> &,const quaternion<T> &)
    q9 = q2*q3;
    
    // using / (const T &, const quaternion<T> &)
    q5 = f0/q2;
    
    // using / (const quaternion<T> &, const T &)
    q6 = q2/f0;
    
    // using / (const ::std::complex<T> &, const quaternion<T> &)
    q7 = c0/q2;
    
    // using / (const quaternion<T> &, const ::std::complex<T> &)
    q8 = q2/c0;
    
    // using / (const quaternion<T> &,const quaternion<T> &)
    q9 = q2/q3;
    
    // using + (const quaternion<T> &)
    q2 = +q0;
    
    // using - (const quaternion<T> &)
    q2 = -q3;
    
    // using == (const T &, const quaternion<T> &)
    f0 == q2;
    
    // using == (const quaternion<T> &, const T &)
    q2 == f0;
    
    // using == (const ::std::complex<T> &, const quaternion<T> &)
    c0 == q2;
    
    // using == (const quaternion<T> &, const ::std::complex<T> &)
    q2 == c0;
    
    // using == (const quaternion<T> &,const quaternion<T> &)
    q2 == q3;
    
    // using != (const T &, const quaternion<T> &)
    f0 != q2;
    
    // using != (const quaternion<T> &, const T &)
    q2 != f0;
    
    // using != (const ::std::complex<T> &, const quaternion<T> &)
    c0 != q2;
    
    // using != (const quaternion<T> &, const ::std::complex<T> &)
    q2 != c0;
    
    // using != (const quaternion<T> &,const quaternion<T> &)
    q2 != q3;
    
    ::std::cout << "Please input a quaternion..." << ::std::endl;
    
#if BOOST_INTERACTIVE_TEST_INPUT_ITERATOR
    ::std::cin >> q0;
    
    if    (::std::cin.fail())
    {
        ::std::cout << "You have entered nonsense!" << ::std::endl;
    }
    else
    {
        ::std::cout << "You have entered the quaternion "
                    << q0 << " ." << ::std::endl;
    }
#else
    ::std::istringstream                    bogus("(1,2,3,4)");
    
    bogus >> q0;
    
    ::std::cout << "You have entered the quaternion "
                << q0 << " ." << ::std::endl;
#endif
    
    ::std::cout << "For this quaternion:" << ::std::endl;
    
    ::std::cout << "the value of "
                << "the real part is "
                << real(q0) << ::std::endl;
    
    ::std::cout << "the value of "
                << "the unreal part is "
                << unreal(q0) << ::std::endl;
    
    ::std::cout << "the value of "
                << "the sup norm is "
                << sup(q0) << ::std::endl;
    
    ::std::cout << "the value of "
                << "the l1 norm is "
                << l1(q0) << ::std::endl;
    
    ::std::cout << "the value of "
                << "the magnitude (euclidian norm) is "
                << abs(q0) << ::std::endl;
    
    ::std::cout << "the value of "
                << "the (Cayley) norm is "
                << norm(q0) << ::std::endl;
    
    ::std::cout << "the value of "
                << "the conjugate is "
                << conj(q0) << ::std::endl;
    
    ::std::cout << "the value of "
                << "the exponential is "
                << exp(q0) << ::std::endl;
    
    ::std::cout << "the value of "
                << "the cube is "
                << pow(q0,3) << ::std::endl;
    
    ::std::cout << "the value of "
                << "the cosinus is "
                << cos(q0) << ::std::endl;
    
    ::std::cout << "the value of "
                << "the sinus is "
                << sin(q0) << ::std::endl;
    
    ::std::cout << "the value of "
                << "the tangent is "
                << tan(q0) << ::std::endl;
    
    ::std::cout << "the value of "
                << "the hyperbolic cosinus is "
                << cosh(q0) << ::std::endl;
    
    ::std::cout << "the value of "
                << "the hyperbolic sinus is "
                << sinh(q0) << ::std::endl;
    
    ::std::cout << "the value of "
                << "the hyperbolic tangent is "
                << tanh(q0) << ::std::endl;
    
#if !defined(__GNUC__) && !defined(__COMO__)
    // somehow, this does not work with either gcc or Comeau :-(
    ::std::cout << "the value of "
                << "the Sinus Cardinal (of index pi) is "
                << sinc_pi(q0) << ::std::endl;
    
    ::std::cout << "the value of "
                << "the Hyperbolic Sinus Cardinal (of index pi) is "
                << sinhc_pi(q0) << ::std::endl;
#endif
    
    ::std::cout << ::std::endl;
    
    float rho = ::std::sqrt(8.0f);
    float theta = ::std::atan(1.0f);
    float phi1 = ::std::atan(1.0f);
    float phi2 = ::std::atan(1.0f);
    
    ::std::cout << "The value of the quaternion represented "
                << "in spherical form by "
                << "rho = " << rho << " , theta = " << theta
                << " , phi1 = " << phi1 << " , phi2 = " << phi2
                << " is "
                << ::boost::math::spherical(rho, theta, phi1, phi2)
                << ::std::endl;
    
    float                            alpha = ::std::atan(1.0f);
    
    ::std::cout << "The value of the quaternion represented "
                << "in semipolar form by "
                << "rho = " << rho << " , alpha = " << alpha
                << " , phi1 = " << phi1 << " , phi2 = " << phi2
                << " is "
                << ::boost::math::semipolar(rho, alpha, phi1, phi2)
                << ::std::endl;
    
    float rho1 = 1;
    float rho2 = 2;
    float theta1 = 0;
    float theta2 = ::std::atan(1.0f)*2;
    
    ::std::cout << "The value of the quaternion represented "
                << "in multipolar form by "
                << "rho1 = " << rho1 << " , theta1 = " << theta1
                << " , rho2 = " << rho2 << " , theta2 = " << theta2
                << " is "
                << ::boost::math::multipolar(rho1, theta1, rho2, theta2)
                << ::std::endl;
    
    float t = 5;
    float radius = ::std::sqrt(2.0f);
    float longitude = ::std::atan(1.0f);
    float lattitude = ::std::atan(::std::sqrt(3.0f));
    
    ::std::cout << "The value of the quaternion represented "
                << "in cylindrospherical form by "
                << "t = " << t << " , radius = " << radius
                << " , longitude = " << longitude << " , latitude = "
                << lattitude << " is "
                << ::boost::math::cylindrospherical(t, radius,
                                                    longitude, lattitude)
                << ::std::endl;
    
    float r = ::std::sqrt(2.0f);
    float angle = ::std::atan(1.0f);
    float h1 = 3;
    float h2 = 4;
    
    ::std::cout << "The value of the quaternion represented "
                << "in cylindrical form by "
                << "r = " << r << " , angle = " << angle
                << " , h1 = " << h1 << " , h2 = " << h2
                << " is "
                << ::boost::math::cylindrical(r, angle, h1, h2)
                << ::std::endl;
    
    double                                   real_1(1);
    ::std::complex<double>                   complex_1(1);
    ::std::complex<double>                   complex_i(0,1);
    ::boost::math::quaternion<double>        quaternion_1(1);
    ::boost::math::quaternion<double>        quaternion_i(0,1);
    ::boost::math::quaternion<double>        quaternion_j(0,0,1);
    ::boost::math::quaternion<double>        quaternion_k(0,0,0,1);
    
    ::std::cout << ::std::endl;
    
    ::std::cout << "Real 1: " << real_1
                << " ; Complex 1: " << complex_1
                << " ; Quaternion 1: " << quaternion_1
                << " ." << ::std::endl;
    
    ::std::cout << "Complex i: " << complex_i
                << " ; Quaternion i: "
                << quaternion_i << " ."
                << ::std::endl;
    
    ::std::cout << "Quaternion j: " << quaternion_j
                << " ." << ::std::endl;
    
    ::std::cout << "Quaternion k: " << quaternion_k
                << " ." << ::std::endl;
    
    
    ::std::cout << ::std::endl;
    
    
    ::std::cout << "i*i: " << quaternion_i*quaternion_i << " ; ";
    ::std::cout << "j*j: " << quaternion_j*quaternion_j << " ; ";
    ::std::cout << "k*k: " << quaternion_k*quaternion_k << " ." << ::std::endl;
    ::std::cout << "i*j: " << quaternion_i*quaternion_j << " ; ";
    ::std::cout << "j*i: " << quaternion_j*quaternion_i << " ." << ::std::endl;
    ::std::cout << "j*k: " << quaternion_j*quaternion_k << " ; ";
    ::std::cout << "k*j: " << quaternion_k*quaternion_j << " ." << ::std::endl;
    ::std::cout << "k*i: " << quaternion_k*quaternion_i << " ; ";
    ::std::cout << "i*k: " << quaternion_i*quaternion_k << " ." << ::std::endl;
    
    ::std::cout << ::std::endl;
    
    
    // tests for evaluation by scripts
        
    using ::std::numeric_limits;
    
    using ::boost::math::abs;
    
#define    BOOST_QUATERNION_MULTIPLICATION_TEST(type)                   \
                                                                        \
    ::std::cout << "Testing multiplication." << std::endl;              \
                                                                        \
    BOOST_CRITICAL_TEST(abs(::boost::math::quaternion<type>(1,0,0,0)*   \
        ::boost::math::quaternion<type>(1,0,0,0)-                       \
        static_cast<type>(1)) <=                                        \
        numeric_limits<type>::epsilon());                               \
                                                                        \
    BOOST_CRITICAL_TEST(abs(::boost::math::quaternion<type>(0,1,0,0)*   \
        ::boost::math::quaternion<type>(0,1,0,0)+                       \
        static_cast<type>(1)) <=                                        \
        numeric_limits<type>::epsilon());                               \
                                                                        \
    BOOST_CRITICAL_TEST(abs(::boost::math::quaternion<type>(0,0,1,0)*   \
        ::boost::math::quaternion<type>(0,0,1,0)+                       \
        static_cast<type>(1)) <=                                        \
        numeric_limits<type>::epsilon());                               \
                                                                        \
    BOOST_CRITICAL_TEST(abs(::boost::math::quaternion<type>(0,0,0,1)*   \
        ::boost::math::quaternion<type>(0,0,0,1)+                       \
        static_cast<type>(1)) <=                                        \
        numeric_limits<type>::epsilon());                               \
                                                                        \
    BOOST_CRITICAL_TEST(abs(::boost::math::quaternion<type>(0,1,0,0)*   \
        ::boost::math::quaternion<type>(0,0,1,0)-                       \
        ::boost::math::quaternion<type>(0,0,0,1)) <=                    \
        numeric_limits<type>::epsilon());                               \
                                                                        \
    BOOST_CRITICAL_TEST(abs(::boost::math::quaternion<type>(0,0,1,0)*   \
        ::boost::math::quaternion<type>(0,1,0,0)+                       \
        ::boost::math::quaternion<type>(0,0,0,1)) <=                    \
        numeric_limits<type>::epsilon());                               \
                                                                        \
    BOOST_CRITICAL_TEST(abs(::boost::math::quaternion<type>(0,0,1,0)*   \
        ::boost::math::quaternion<type>(0,0,0,1)-                       \
        ::boost::math::quaternion<type>(0,1,0,0)) <=                    \
        numeric_limits<type>::epsilon());                               \
                                                                        \
    BOOST_CRITICAL_TEST(abs(::boost::math::quaternion<type>(0,0,0,1)*   \
        ::boost::math::quaternion<type>(0,0,1,0)+                       \
        ::boost::math::quaternion<type>(0,1,0,0)) <=                    \
        numeric_limits<type>::epsilon());                               \
                                                                        \
    BOOST_CRITICAL_TEST(abs(::boost::math::quaternion<type>(0,0,0,1)*   \
        ::boost::math::quaternion<type>(0,1,0,0)-                       \
        ::boost::math::quaternion<type>(0,0,1,0)) <=                    \
        numeric_limits<type>::epsilon());                               \
                                                                        \
    BOOST_CRITICAL_TEST(abs(::boost::math::quaternion<type>(0,1,0,0)*   \
        ::boost::math::quaternion<type>(0,0,0,1)+                       \
        ::boost::math::quaternion<type>(0,0,1,0)) <=                    \
        numeric_limits<type>::epsilon());
    
    
#define    BOOST_QUATERNION_EXP_TEST(type)                              \
                                                                        \
    ::std::cout << "Testing exp." << std::endl;                         \
                                                                        \
    BOOST_TEST(abs(exp(::boost::math::quaternion<type>                  \
            (0,4*::std::atan(static_cast<type>(1)),0,0)                 \
        )+static_cast<type>(1)) <= 2*numeric_limits<type>::epsilon());  \
                                                                        \
    BOOST_TEST(abs(exp(::boost::math::quaternion<type>                  \
            (0,0,4*::std::atan(static_cast<type>(1)),0)                 \
        )+static_cast<type>(1)) <= 2*numeric_limits<type>::epsilon());  \
                                                                        \
    BOOST_TEST(abs(exp(::boost::math::quaternion<type>                  \
            (0,0,0,4*::std::atan(static_cast<type>(1)))                 \
        )+static_cast<type>(1)) <= 2*numeric_limits<type>::epsilon());
    
    
#define    BOOST_QUATERNION_COS_TEST(type)                              \
                                                                        \
    ::std::cout << "Testing cos." << std::endl;                         \
                                                                        \
    BOOST_TEST(abs(static_cast<type>(4)*                                \
        cos(::boost::math::quaternion<type>                             \
            (0,::std::log(static_cast<type>(2)),0,0)                    \
        )-static_cast<type>(5)) <= 4*numeric_limits<type>::epsilon());  \
                                                                        \
    BOOST_TEST(abs(static_cast<type>(4)*                                \
        cos(::boost::math::quaternion<type>                             \
            (0,0,::std::log(static_cast<type>(2)),0)                    \
        )-static_cast<type>(5)) <= 4*numeric_limits<type>::epsilon());  \
                                                                        \
    BOOST_TEST(abs(static_cast<type>(4)*                                \
        cos(::boost::math::quaternion<type>                             \
            (0,0,0,::std::log(static_cast<type>(2)))                    \
        )-static_cast<type>(5)) <= 4*numeric_limits<type>::epsilon());
    
    
#define    BOOST_QUATERNION_SIN_TEST(type)                              \
                                                                        \
    ::std::cout << "Testing sin." << std::endl;                         \
                                                                        \
    BOOST_TEST(abs(static_cast<type>(4)*                                \
        sin(::boost::math::quaternion<type>                             \
            (0,::std::log(static_cast<type>(2)),0,0)                    \
        )-::boost::math::quaternion<type>(0,3,0,0)) <=                  \
        4*numeric_limits<type>::epsilon());                             \
                                                                        \
    BOOST_TEST(abs(static_cast<type>(4)*                                \
        sin(::boost::math::quaternion<type>                             \
            (0,0,::std::log(static_cast<type>(2)),0)                    \
        )-::boost::math::quaternion<type>(0,0,3,0)) <=                  \
        4*numeric_limits<type>::epsilon());                             \
                                                                        \
    BOOST_TEST(abs(static_cast<type>(4)*                                \
        sin(::boost::math::quaternion<type>                             \
            (0,0,0,::std::log(static_cast<type>(2)))                    \
        )-::boost::math::quaternion<type>(0,0,0,3)) <=                  \
        4*numeric_limits<type>::epsilon());
    
    
#define    BOOST_QUATERNION_COSH_TEST(type)                             \
                                                                        \
    ::std::cout << "Testing cosh." << std::endl;                        \
                                                                        \
    BOOST_TEST(abs(                                                     \
        cosh(::boost::math::quaternion<type>                            \
            (0,4*::std::atan(static_cast<type>(1)),0,0)                 \
        )+static_cast<type>(1)) <= 4*numeric_limits<type>::epsilon());  \
                                                                        \
    BOOST_TEST(abs(                                                     \
        cosh(::boost::math::quaternion<type>                            \
            (0,0,4*::std::atan(static_cast<type>(1)),0)                 \
        )+static_cast<type>(1)) <= 4*numeric_limits<type>::epsilon());  \
                                                                        \
    BOOST_TEST(abs(                                                     \
        cosh(::boost::math::quaternion<type>                            \
            (0,0,0,4*::std::atan(static_cast<type>(1)))                 \
        )+static_cast<type>(1)) <= 4*numeric_limits<type>::epsilon());
    
    
#define    BOOST_QUATERNION_SINH_TEST(type)                             \
                                                                        \
    ::std::cout << "Testing sinh." << std::endl;                        \
                                                                        \
    BOOST_TEST(abs(                                                     \
        sinh(::boost::math::quaternion<type>                            \
            (0,2*::std::atan(static_cast<type>(1)),0,0)                 \
        )-::boost::math::quaternion<type>(0,1,0,0)) <=                  \
        4*numeric_limits<type>::epsilon());                             \
                                                                        \
    BOOST_TEST(abs(                                                     \
        sinh(::boost::math::quaternion<type>                            \
            (0,0,2*::std::atan(static_cast<type>(1)),0)                 \
        )-::boost::math::quaternion<type>(0,0,1,0)) <=                  \
        4*numeric_limits<type>::epsilon());                             \
                                                                        \
    BOOST_TEST(abs(                                                     \
        sinh(::boost::math::quaternion<type>                            \
            (0,0,0,2*::std::atan(static_cast<type>(1)))                 \
        )-::boost::math::quaternion<type>(0,0,0,1)) <=                  \
        4*numeric_limits<type>::epsilon());
    
    
#define    BOOST_QUATERNION_TRENSCENDENTALS_TEST(type)                  \
                                                                        \
    BOOST_QUATERNION_EXP_TEST(type)                                     \
    BOOST_QUATERNION_COS_TEST(type)                                     \
    BOOST_QUATERNION_SIN_TEST(type)                                     \
    BOOST_QUATERNION_COSH_TEST(type)                                    \
    BOOST_QUATERNION_SINH_TEST(type)
    
    
#if defined(__GNUC__) || defined(__COMO__) || \
    (defined(__MWERKS__) && (__MWERKS__ <= 0x2301))
#define    BOOST_QUATERNION_TEST(type)                      \
                                                            \
    ::std::cout << "Testing " << #type << "." << std::endl; \
                                                            \
    BOOST_QUATERNION_MULTIPLICATION_TEST(type)
#else
#define    BOOST_QUATERNION_TEST(type)                      \
                                                            \
    ::std::cout << "Testing " << #type << "." << std::endl; \
                                                            \
    BOOST_QUATERNION_MULTIPLICATION_TEST(type)              \
    BOOST_QUATERNION_TRENSCENDENTALS_TEST(type)
#endif
    
    
    BOOST_QUATERNION_TEST(float)
    BOOST_QUATERNION_TEST(double)
    BOOST_QUATERNION_TEST(long double)
    
    
#undef    BOOST_QUATERNION_TEST
    
#undef    BOOST_QUATERNION_MULTIPLICATION_TEST
#undef    BOOST_QUATERNION_TRENSCENDENTALS_TEST
    
#undef    BOOST_QUATERNION_EXP_TEST
#undef    BOOST_QUATERNION_COS_TEST
#undef    BOOST_QUATERNION_SIN_TEST
#undef    BOOST_QUATERNION_COSH_TEST
#undef    BOOST_QUATERNION_SINH_TEST
    
    
    return(::boost::exit_success);
}
