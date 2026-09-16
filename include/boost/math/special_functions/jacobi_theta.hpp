// Jacobi theta functions
// Copyright Evan Miller 2020
//
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0. (See accompanying file
// LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Four main theta functions with various flavors of parameterization,
// floating-point policies, and bonus "minus 1" versions of functions 3 and 4
// designed to preserve accuracy for small q. Twenty-four C++ functions are
// provided in all.
//
// The functions take a real argument z and a parameter known as q, or its close
// relative tau.
//
// The mathematical functions are best understood in terms of their Fourier
// series. Using the q parameterization, and summing from n = 0 to INF:
//
// theta_1(z,q) = 2 SUM (-1)^n * q^(n+1/2)^2 * sin((2n+1)z)
// theta_2(z,q) = 2 SUM q^(n+1/2)^2 * cos((2n+1)z)
// theta_3(z,q) = 1 + 2 SUM q^n^2 * cos(2nz)
// theta_4(z,q) = 1 + 2 SUM (-1)^n * q^n^2 * cos(2nz)
//
// Appropriately multiplied and divided, these four theta functions can be used
// to implement the famous Jacabi elliptic functions - but this is not really
// recommended, as the existing Boost implementations are likely faster and
// more accurate.  More saliently, setting z = 0 on the fourth theta function
// will produce the limiting CDF of the Kolmogorov-Smirnov distribution, which
// is this particular implementation's raison d'etre.
//
// Separate C++ functions are provided for q and for tau. The main q functions are:
//
// template <class T> inline T jacobi_theta1(T z, T q);
// template <class T> inline T jacobi_theta2(T z, T q);
// template <class T> inline T jacobi_theta3(T z, T q);
// template <class T> inline T jacobi_theta4(T z, T q);
//
// The parameter q, also known as the nome, is restricted to the domain (0, 1),
// and will throw a domain error otherwise.
//
// The equivalent functions that use tau instead of q are:
//
// template <class T> inline T jacobi_theta1tau(T z, T tau);
// template <class T> inline T jacobi_theta2tau(T z, T tau);
// template <class T> inline T jacobi_theta3tau(T z, T tau);
// template <class T> inline T jacobi_theta4tau(T z, T tau);
//
// Mathematically, q and tau are related by:
//
// q = exp(i PI*Tau)
//
// However, the tau in the equation above is *not* identical to the tau in the function
// signature. Instead, `tau` is the imaginary component of tau. Mathematically, tau can
// be complex - but practically, most applications call for a purely imaginary tau.
// Rather than provide a full complex-number API, the author decided to treat the
// parameter `tau` as an imaginary number. So in computational terms, the
// relationship between `q` and `tau` is given by:
//
// q = exp(-constants::pi<T>() * tau)
//
// The tau versions are provided for the sake of accuracy, as well as conformance
// with common notation. If your q is an exponential, you are better off using
// the tau versions, e.g.
//
// jacobi_theta1(z, exp(-a)); // rather poor accuracy
// jacobi_theta1tau(z, a / constants::pi<T>()); // better accuracy
//
// Similarly, if you have a precise (small positive) value for the complement
// of q, you can obtain a more precise answer overall by passing the result of
// `log1p` to the tau parameter:
//
// jacobi_theta1(z, 1-q_complement); // precision lost in subtraction
// jacobi_theta1tau(z, -log1p(-q_complement) / constants::pi<T>()); // better!
//
// A third quartet of functions are provided for improving accuracy in cases
// where q is small, specifically |q| < exp(-PI) = 0.04 (or, equivalently, tau
// greater than unity). In this domain of q values, the third and fourth theta
// functions always return values close to 1. So the following "m1" functions
// are provided, similar in spirit to `expm1`, which return one less than their
// regular counterparts:
//
// template <class T> inline T jacobi_theta3m1(T z, T q);
// template <class T> inline T jacobi_theta4m1(T z, T q);
// template <class T> inline T jacobi_theta3m1tau(T z, T tau);
// template <class T> inline T jacobi_theta4m1tau(T z, T tau);
//
// Note that "m1" versions of the first and second theta would not be useful,
// as their ranges are not confined to a neighborhood around 1 (see the Fourier
// transform representations above).
//
// Finally, the twelve functions above are each available with a third Policy
// argument, which can be used to define a custom epsilon value. These Policy
// versions bring the total number of functions provided by jacobi_theta.hpp
// to twenty-four.
//
// See:
// https://mathworld.wolfram.com/JacobiThetaFunctions.html
// https://dlmf.nist.gov/20

#ifndef BOOST_MATH_JACOBI_THETA_HPP
#define BOOST_MATH_JACOBI_THETA_HPP

#include <boost/math/tools/complex.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/math/tools/promotion.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/expm1.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

namespace boost{ namespace math{

// Simple functions - parameterized by q
BOOST_MATH_EXPORT template <class T, class U>
inline typename tools::promote_args<T, U>::type jacobi_theta1(T z, U q);
BOOST_MATH_EXPORT template <class T, class U>
inline typename tools::promote_args<T, U>::type jacobi_theta2(T z, U q);
BOOST_MATH_EXPORT template <class T, class U>
inline typename tools::promote_args<T, U>::type jacobi_theta3(T z, U q);
BOOST_MATH_EXPORT template <class T, class U>
inline typename tools::promote_args<T, U>::type jacobi_theta4(T z, U q);

// Simple functions - parameterized by tau (assumed imaginary)
// q = exp(i*PI*TAU)
// tau = -log(q)/PI
BOOST_MATH_EXPORT template <class T, class U>
inline typename tools::promote_args<T, U>::type jacobi_theta1tau(T z, U tau);
BOOST_MATH_EXPORT template <class T, class U>
inline typename tools::promote_args<T, U>::type jacobi_theta2tau(T z, U tau);
BOOST_MATH_EXPORT template <class T, class U>
inline typename tools::promote_args<T, U>::type jacobi_theta3tau(T z, U tau);
BOOST_MATH_EXPORT template <class T, class U>
inline typename tools::promote_args<T, U>::type jacobi_theta4tau(T z, U tau);

// Minus one versions for small q / large tau
BOOST_MATH_EXPORT template <class T, class U>
inline typename tools::promote_args<T, U>::type jacobi_theta3m1(T z, U q);
BOOST_MATH_EXPORT template <class T, class U>
inline typename tools::promote_args<T, U>::type jacobi_theta4m1(T z, U q);
BOOST_MATH_EXPORT template <class T, class U>
inline typename tools::promote_args<T, U>::type jacobi_theta3m1tau(T z, U tau);
BOOST_MATH_EXPORT template <class T, class U>
inline typename tools::promote_args<T, U>::type jacobi_theta4m1tau(T z, U tau);

// Policied versions - parameterized by q
BOOST_MATH_EXPORT template <class T, class U, class Policy>
inline typename tools::promote_args<T, U>::type jacobi_theta1(T z, U q, const Policy& pol);
BOOST_MATH_EXPORT template <class T, class U, class Policy>
inline typename tools::promote_args<T, U>::type jacobi_theta2(T z, U q, const Policy& pol);
BOOST_MATH_EXPORT template <class T, class U, class Policy>
inline typename tools::promote_args<T, U>::type jacobi_theta3(T z, U q, const Policy& pol);
BOOST_MATH_EXPORT template <class T, class U, class Policy>
inline typename tools::promote_args<T, U>::type jacobi_theta4(T z, U q, const Policy& pol);

// Policied versions - parameterized by tau
BOOST_MATH_EXPORT template <class T, class U, class Policy>
inline typename tools::promote_args<T, U>::type jacobi_theta1tau(T z, U tau, const Policy& pol);
BOOST_MATH_EXPORT template <class T, class U, class Policy>
inline typename tools::promote_args<T, U>::type jacobi_theta2tau(T z, U tau, const Policy& pol);
BOOST_MATH_EXPORT template <class T, class U, class Policy>
inline typename tools::promote_args<T, U>::type jacobi_theta3tau(T z, U tau, const Policy& pol);
BOOST_MATH_EXPORT template <class T, class U, class Policy>
inline typename tools::promote_args<T, U>::type jacobi_theta4tau(T z, U tau, const Policy& pol);

// Policied m1 functions
BOOST_MATH_EXPORT template <class T, class U, class Policy>
inline typename tools::promote_args<T, U>::type jacobi_theta3m1(T z, U q, const Policy& pol);
BOOST_MATH_EXPORT template <class T, class U, class Policy>
inline typename tools::promote_args<T, U>::type jacobi_theta4m1(T z, U q, const Policy& pol);
BOOST_MATH_EXPORT template <class T, class U, class Policy>
inline typename tools::promote_args<T, U>::type jacobi_theta3m1tau(T z, U tau, const Policy& pol);
BOOST_MATH_EXPORT template <class T, class U, class Policy>
inline typename tools::promote_args<T, U>::type jacobi_theta4m1tau(T z, U tau, const Policy& pol);

// Compare the non-oscillating component of the delta to the previous delta.
// Both are assumed to be non-negative. Written so that a NaN delta counts as
// converged: otherwise a NaN would never satisfy the test and the summation
// loops below would never terminate.
template <class RealType>
inline bool
_jacobi_theta_converged(RealType last_delta, RealType delta, RealType eps) {
    return !(delta > eps*last_delta);
}

template <class RealType, class Policy>
inline bool
_jacobi_theta_check_z(RealType z, const Policy& pol, const char* function, RealType* result) {
    if (!(boost::math::isfinite)(z)) {
        *result = policies::raise_domain_error<RealType>(function, "z must be finite but got %1%.", z, pol);
        return false;
    }
    return true;
}

template <class RealType, class Policy>
inline bool
_jacobi_theta_check_tau(RealType tau, const Policy& pol, const char* function, RealType* result) {
    // The negated comparison also rejects NaN.
    if (!(tau > 0)) {
        *result = policies::raise_domain_error<RealType>(function, "tau must be greater than 0 but got %1%.", tau, pol);
        return false;
    }
    return true;
}

template <class RealType, class Policy>
inline bool
_jacobi_theta_check_q(RealType q, const Policy& pol, const char* function, RealType* result) {
    // The negated comparison also rejects NaN.
    if (!(q > 0 && q < 1)) {
        *result = policies::raise_domain_error<RealType>(function, "q must be greater than 0 and less than 1 but got %1%.", q, pol);
        return false;
    }
    return true;
}

// Powers of the nome for the direct Fourier series below. When the caller
// supplies q directly we raise q to the power with pow(), which is accurate to
// about an ulp for any exponent. Going through tau = -log(q)/pi and back via
// exp() would multiply the rounding error by |log q|, which is exactly the
// small-q regime the "m1" functions exist to serve.
template <class RealType>
struct _jacobi_theta_q_power {
    RealType q;
    RealType operator()(RealType exponent) const {
        BOOST_MATH_STD_USING
        return pow(q, exponent);
    }
};

template <class RealType>
struct _jacobi_theta_tau_power {
    RealType tau;
    RealType operator()(RealType exponent) const {
        BOOST_MATH_STD_USING
        return exp(-tau * constants::pi<RealType>() * exponent);
    }
};

// Direct Fourier series (DLMF 20.2.1 - 20.2.4). These converge quickly when
// q < exp(-Pi), i.e. tau > 1, and are used in that regime by both the q and
// the tau parameterizations. NomePower(e) must return q^e.

// = 2 * SUM (-1)^n * q^(n+1/2)^2 * sin((2n+1)z)
template <class RealType, class NomePower, class Policy>
inline RealType
_jacobi_theta1_series(RealType z, const NomePower& q_pow, const Policy&) {
    BOOST_MATH_STD_USING
    unsigned n = 0;
    RealType eps = policies::get_epsilon<RealType, Policy>();
    RealType q_n = 0, last_q_n, delta, result = 0;

    do {
        last_q_n = q_n;
        q_n = q_pow(RealType(n + 0.5)*RealType(n + 0.5));
        delta = q_n * sin(RealType(2*n+1)*z);
        if (n%2)
            delta = -delta;

        result += delta + delta;
        n++;
    } while (!_jacobi_theta_converged(last_q_n, q_n, eps));

    return result;
}

// = 2 * SUM q^(n+1/2)^2 * cos((2n+1)z)
template <class RealType, class NomePower, class Policy>
inline RealType
_jacobi_theta2_series(RealType z, const NomePower& q_pow, const Policy&) {
    BOOST_MATH_STD_USING
    unsigned n = 0;
    RealType eps = policies::get_epsilon<RealType, Policy>();
    RealType q_n = 0, last_q_n, delta, result = 0;

    do {
        last_q_n = q_n;
        q_n = q_pow(RealType(n + 0.5)*RealType(n + 0.5));
        delta = q_n * cos(RealType(2*n+1)*z);
        result += delta + delta;
        n++;
    } while (!_jacobi_theta_converged(last_q_n, q_n, eps));

    return result;
}

// = 2 * SUM q^n^2 * cos(2nz), n >= 1 (i.e. theta3 minus one)
template <class RealType, class NomePower, class Policy>
inline RealType
_jacobi_theta3m1_series(RealType z, const NomePower& q_pow, const Policy&) {
    BOOST_MATH_STD_USING
    unsigned n = 1;
    RealType eps = policies::get_epsilon<RealType, Policy>();
    RealType q_n = 0, last_q_n, delta, result = 0;

    do {
        last_q_n = q_n;
        q_n = q_pow(RealType(n)*RealType(n));
        delta = q_n * cos(RealType(2*n)*z);
        result += delta + delta;
        n++;
    } while (!_jacobi_theta_converged(last_q_n, q_n, eps));

    return result;
}

// = 2 * SUM (-1)^n q^n^2 * cos(2nz), n >= 1 (i.e. theta4 minus one)
template <class RealType, class NomePower, class Policy>
inline RealType
_jacobi_theta4m1_series(RealType z, const NomePower& q_pow, const Policy&) {
    BOOST_MATH_STD_USING
    unsigned n = 1;
    RealType eps = policies::get_epsilon<RealType, Policy>();
    RealType q_n = 0, last_q_n, delta, result = 0;

    do {
        last_q_n = q_n;
        q_n = q_pow(RealType(n)*RealType(n));
        delta = q_n * cos(RealType(2*n)*z);
        if (n%2)
            delta = -delta;

        result += delta + delta;
        n++;
    } while (!_jacobi_theta_converged(last_q_n, q_n, eps));

    return result;
}

template <class RealType>
inline RealType
_jacobi_theta_sum(RealType tau, RealType z_n, RealType z_increment, RealType eps) {
    BOOST_MATH_STD_USING
    RealType delta = 0, partial_result = 0;
    RealType last_delta = 0;

    do {
        last_delta = delta;
        delta = exp(-tau*z_n*z_n/constants::pi<RealType>());
        partial_result += delta;
        z_n += z_increment;
    } while (!_jacobi_theta_converged(last_delta, delta, eps));

    return partial_result;
}

// The following _IMAGINARY theta functions assume imaginary z and are for
// internal use only. They are designed to increase accuracy and reduce the
// number of iterations required for convergence for large |q|. The z argument
// is scaled by tau, and the summations are rewritten to be double-sided
// following DLMF 20.13.4 and 20.13.5. Each term is a Gaussian
// exp(-tau*(z - c)^2/Pi) centered at a multiple of Pi or Pi/2, and the
// results are scaled by sqrt(tau).
//
// These functions are triggered when tau < 1, i.e. |q| > exp(-Pi) = 0.043
//
// Note that jacobi_theta4 uses the imaginary version of jacobi_theta2 (and
// vice-versa). jacobi_theta1 and jacobi_theta3 use the imaginary versions of
// themselves, following DLMF 20.7.30 - 20.7.33.

// theta1(z|i/tau) = sqrt(tau) * SUM_{n>=0} (-1)^n [G(z - c_n) - G(z + c_n)]
// with c_n = Pi*(n+1/2) and G(x) = exp(-tau*x^2/Pi).
//
// Each bracket is a difference of two Gaussians which nearly cancel when z is
// small, so it is evaluated instead through the exact identity
//   G(z - c) - G(z + c) = -G(z - c) * expm1(-2*tau*z*(2n+1)),
// which keeps full relative precision all the way down to z -> 0.
// Requires 0 <= z <= Pi/2; the caller reduces z into this range.
template <class RealType, class Policy>
inline RealType
_IMAGINARY_jacobi_theta1tau(RealType z, RealType tau, const Policy& pol) {
    BOOST_MATH_STD_USING
    RealType eps = policies::get_epsilon<RealType, Policy>();
    RealType result = 0, g = 0, last_g, c, pair;
    unsigned n = 0;

    do {
        last_g = g;
        c = constants::pi<RealType>() * RealType(n + 0.5);
        g = exp(-tau * (z - c) * (z - c) / constants::pi<RealType>());
        pair = -g * boost::math::expm1(RealType(-2 * tau * z * RealType(2*n + 1)), pol);
        if (n%2)
            pair = -pair;

        result += pair;
        n++;
    } while (!_jacobi_theta_converged(last_g, g, eps));

    return result * sqrt(tau);
}

template <class RealType, class Policy>
inline RealType
_IMAGINARY_jacobi_theta2tau(RealType z, RealType tau, const Policy&) {
    BOOST_MATH_STD_USING
    RealType eps = policies::get_epsilon<RealType, Policy>();
    RealType result = RealType(0);

    // n>=0
    result += _jacobi_theta_sum(tau, RealType(z + constants::half_pi<RealType>()), constants::pi<RealType>(), eps);
    // n<0
    result += _jacobi_theta_sum(tau, RealType(z - constants::half_pi<RealType>()), RealType (-constants::pi<RealType>()), eps);

    return result * sqrt(tau);
}

template <class RealType, class Policy>
inline RealType
_IMAGINARY_jacobi_theta3tau(RealType z, RealType tau, const Policy&) {
    BOOST_MATH_STD_USING
    RealType eps = policies::get_epsilon<RealType, Policy>();
    RealType result = 0;

    // n=0
    result += exp(-z*z*tau/constants::pi<RealType>());
    // n>0
    result += _jacobi_theta_sum(tau, RealType(z + constants::pi<RealType>()), constants::pi<RealType>(), eps);
    // n<0
    result += _jacobi_theta_sum(tau, RealType(z - constants::pi<RealType>()), RealType(-constants::pi<RealType>()), eps);

    return result * sqrt(tau);
}

template <class RealType, class Policy>
inline RealType
_IMAGINARY_jacobi_theta4tau(RealType z, RealType tau, const Policy&) {
    BOOST_MATH_STD_USING
    RealType eps = policies::get_epsilon<RealType, Policy>();
    RealType result = 0;

    // n = 0
    result += exp(-z*z*tau/constants::pi<RealType>());

    // n > 0 odd
    result -= _jacobi_theta_sum(tau, RealType(z + constants::pi<RealType>()), constants::two_pi<RealType>(), eps);
    // n < 0 odd
    result -= _jacobi_theta_sum(tau, RealType(z - constants::pi<RealType>()), RealType (-constants::two_pi<RealType>()), eps);
    // n > 0 even
    result += _jacobi_theta_sum(tau, RealType(z + constants::two_pi<RealType>()), constants::two_pi<RealType>(), eps);
    // n < 0 even
    result += _jacobi_theta_sum(tau, RealType(z - constants::two_pi<RealType>()), RealType (-constants::two_pi<RealType>()), eps);

    return result * sqrt(tau);
}

// First Jacobi theta function (Parameterized by tau - assumed imaginary)
// = 2 * SUM (-1)^n * exp(i*Pi*Tau*(n+1/2)^2) * sin((2n+1)z)
template <class RealType, class Policy>
inline RealType
jacobi_theta1tau_imp(RealType z, RealType tau, const Policy& pol, const char *function)
{
    BOOST_MATH_STD_USING
    RealType result = 0;

    if (!_jacobi_theta_check_tau(tau, pol, function, &result))
        return result;
    if (!_jacobi_theta_check_z(z, pol, function, &result))
        return result;

    if (abs(z) == 0.0)
        return result;

    if (tau < 1.0) {
        // Reduce to -Pi <= z <= Pi (theta1 has period 2*Pi)...
        z = fmod(z, constants::two_pi<RealType>());
        while (z > constants::pi<RealType>()) {
            z -= constants::two_pi<RealType>();
        }
        while (z < -constants::pi<RealType>()) {
            z += constants::two_pi<RealType>();
        }
        // ...then to -Pi/2 <= z <= Pi/2 using theta1(z + Pi) = -theta1(z)...
        RealType sign = 1;
        if (z > constants::half_pi<RealType>()) {
            z -= constants::pi<RealType>();
            sign = -sign;
        } else if (z < -constants::half_pi<RealType>()) {
            z += constants::pi<RealType>();
            sign = -sign;
        }
        // ...and finally to 0 <= z <= Pi/2 since theta1 is odd.
        if (z < 0) {
            z = -z;
            sign = -sign;
        }

        return sign * _IMAGINARY_jacobi_theta1tau(z, RealType(1/tau), pol);
    }

    return _jacobi_theta1_series(z, _jacobi_theta_tau_power<RealType>{tau}, pol);
}

// First Jacobi theta function (Parameterized by q)
// = 2 * SUM (-1)^n * q^(n+1/2)^2 * sin((2n+1)z)
template <class RealType, class Policy>
inline RealType
jacobi_theta1_imp(RealType z, RealType q, const Policy& pol, const char *function) {
    BOOST_MATH_STD_USING
    RealType result = 0;

    if (!_jacobi_theta_check_q(q, pol, function, &result))
        return result;
    if (!_jacobi_theta_check_z(z, pol, function, &result))
        return result;

    if (q < exp(-constants::pi<RealType>()))
        return _jacobi_theta1_series(z, _jacobi_theta_q_power<RealType>{q}, pol);

    return jacobi_theta1tau_imp(z, RealType (-log(q)/constants::pi<RealType>()), pol, function);
}

// Second Jacobi theta function (Parameterized by tau - assumed imaginary)
// = 2 * SUM exp(i*Pi*Tau*(n+1/2)^2) * cos((2n+1)z)
template <class RealType, class Policy>
inline RealType
jacobi_theta2tau_imp(RealType z, RealType tau, const Policy& pol, const char *function)
{
    BOOST_MATH_STD_USING
    RealType result = 0;

    if (!_jacobi_theta_check_tau(tau, pol, function, &result))
        return result;
    if (!_jacobi_theta_check_z(z, pol, function, &result))
        return result;

    if (tau < 1.0 && abs(z) == 0.0) {
        return jacobi_theta4tau(z, RealType(1/tau), pol) / sqrt(tau);
    } else if (tau < 1.0) { // DLMF 20.7.31
        z = fmod(z, constants::two_pi<RealType>());
        while (z > constants::pi<RealType>()) {
            z -= constants::two_pi<RealType>();
        }
        while (z < -constants::pi<RealType>()) {
            z += constants::two_pi<RealType>();
        }

        return _IMAGINARY_jacobi_theta4tau(z, RealType(1/tau), pol);
    }

    return _jacobi_theta2_series(z, _jacobi_theta_tau_power<RealType>{tau}, pol);
}

// Second Jacobi theta function, parameterized by q
// = 2 * SUM q^(n+1/2)^2 * cos((2n+1)z)
template <class RealType, class Policy>
inline RealType
jacobi_theta2_imp(RealType z, RealType q, const Policy& pol, const char *function) {
    BOOST_MATH_STD_USING
    RealType result = 0;

    if (!_jacobi_theta_check_q(q, pol, function, &result))
        return result;
    if (!_jacobi_theta_check_z(z, pol, function, &result))
        return result;

    if (q < exp(-constants::pi<RealType>()))
        return _jacobi_theta2_series(z, _jacobi_theta_q_power<RealType>{q}, pol);

    return jacobi_theta2tau_imp(z, RealType (-log(q)/constants::pi<RealType>()), pol, function);
}

// Third Jacobi theta function, parameterized by tau
// = 1 + 2 * SUM exp(i*Pi*Tau*(n)^2) * cos(2nz)
template <class RealType, class Policy>
inline RealType
jacobi_theta3tau_imp(RealType z, RealType tau, const Policy& pol, const char *function)
{
    BOOST_MATH_STD_USING
    RealType result = 0;

    if (!_jacobi_theta_check_tau(tau, pol, function, &result))
        return result;
    if (!_jacobi_theta_check_z(z, pol, function, &result))
        return result;

    if (tau < 1.0 && abs(z) == 0.0) {
        return jacobi_theta3tau(z, RealType(1/tau), pol) / sqrt(tau);
    } else if (tau < 1.0) { // DLMF 20.7.32
        z = fmod(z, constants::pi<RealType>());
        while (z > constants::half_pi<RealType>()) {
            z -= constants::pi<RealType>();
        }
        while (z < -constants::half_pi<RealType>()) {
            z += constants::pi<RealType>();
        }
        return _IMAGINARY_jacobi_theta3tau(z, RealType(1/tau), pol);
    }
    return RealType(1) + _jacobi_theta3m1_series(z, _jacobi_theta_tau_power<RealType>{tau}, pol);
}

// Third Jacobi theta function, minus one (Parameterized by tau - assumed imaginary)
// This function preserves accuracy for small values of q (i.e. |q| < exp(-Pi) = 0.043)
// For larger values of q, the minus one version usually won't help.
// = 2 * SUM exp(i*Pi*Tau*(n)^2) * cos(2nz)
template <class RealType, class Policy>
inline RealType
jacobi_theta3m1tau_imp(RealType z, RealType tau, const Policy& pol, const char *function)
{
    BOOST_MATH_STD_USING
    RealType result = 0;

    if (!_jacobi_theta_check_tau(tau, pol, function, &result))
        return result;
    if (!_jacobi_theta_check_z(z, pol, function, &result))
        return result;

    if (tau < 1.0)
        return jacobi_theta3tau_imp(z, tau, pol, function) - RealType(1);

    return _jacobi_theta3m1_series(z, _jacobi_theta_tau_power<RealType>{tau}, pol);
}

// Third Jacobi theta function, minus one (parameterized by q)
// = 2 * SUM q^n^2 * cos(2nz)
template <class RealType, class Policy>
inline RealType
jacobi_theta3m1_imp(RealType z, RealType q, const Policy& pol, const char *function) {
    BOOST_MATH_STD_USING
    RealType result = 0;

    if (!_jacobi_theta_check_q(q, pol, function, &result))
        return result;
    if (!_jacobi_theta_check_z(z, pol, function, &result))
        return result;

    if (q < exp(-constants::pi<RealType>()))
        return _jacobi_theta3m1_series(z, _jacobi_theta_q_power<RealType>{q}, pol);

    return jacobi_theta3m1tau_imp(z, RealType (-log(q)/constants::pi<RealType>()), pol, function);
}

// Third Jacobi theta function (parameterized by q)
// = 1 + 2 * SUM q^n^2 * cos(2nz)
template <class RealType, class Policy>
inline RealType
jacobi_theta3_imp(RealType z, RealType q, const Policy& pol, const char *function) {
    BOOST_MATH_STD_USING
    RealType result = 0;

    if (!_jacobi_theta_check_q(q, pol, function, &result))
        return result;
    if (!_jacobi_theta_check_z(z, pol, function, &result))
        return result;

    if (q < exp(-constants::pi<RealType>()))
        return RealType(1) + _jacobi_theta3m1_series(z, _jacobi_theta_q_power<RealType>{q}, pol);

    return jacobi_theta3tau_imp(z, RealType (-log(q)/constants::pi<RealType>()), pol, function);
}

// Fourth Jacobi theta function (Parameterized by tau)
// = 1 + 2 * SUM (-1)^n exp(i*Pi*Tau*(n)^2) * cos(2nz)
template <class RealType, class Policy>
inline RealType
jacobi_theta4tau_imp(RealType z, RealType tau, const Policy& pol, const char *function)
{
    BOOST_MATH_STD_USING
    RealType result = 0;

    if (!_jacobi_theta_check_tau(tau, pol, function, &result))
        return result;
    if (!_jacobi_theta_check_z(z, pol, function, &result))
        return result;

    if (tau < 1.0 && abs(z) == 0.0) {
        return jacobi_theta2tau(z, RealType(1/tau), pol) / sqrt(tau);
    } else if (tau < 1.0) { // DLMF 20.7.33
        z = fmod(z, constants::pi<RealType>());
        while (z > constants::half_pi<RealType>()) {
            z -= constants::pi<RealType>();
        }
        while (z < -constants::half_pi<RealType>()) {
            z += constants::pi<RealType>();
        }
        return _IMAGINARY_jacobi_theta2tau(z, RealType(1/tau), pol);
    }

    return RealType(1) + _jacobi_theta4m1_series(z, _jacobi_theta_tau_power<RealType>{tau}, pol);
}

// Fourth Jacobi theta function, minus one (Parameterized by tau)
// This function preserves accuracy for small values of q (i.e. tau > 1)
// = 2 * SUM (-1)^n exp(i*Pi*Tau*(n)^2) * cos(2nz)
template <class RealType, class Policy>
inline RealType
jacobi_theta4m1tau_imp(RealType z, RealType tau, const Policy& pol, const char *function)
{
    BOOST_MATH_STD_USING
    RealType result = 0;

    if (!_jacobi_theta_check_tau(tau, pol, function, &result))
        return result;
    if (!_jacobi_theta_check_z(z, pol, function, &result))
        return result;

    if (tau < 1.0)
        return jacobi_theta4tau_imp(z, tau, pol, function) - RealType(1);

    return _jacobi_theta4m1_series(z, _jacobi_theta_tau_power<RealType>{tau}, pol);
}

// Fourth Jacobi theta function, minus one (Parameterized by q)
// This function preserves accuracy for small values of q
// = 2 * SUM (-1)^n q^n^2 * cos(2nz)
template <class RealType, class Policy>
inline RealType
jacobi_theta4m1_imp(RealType z, RealType q, const Policy& pol, const char *function) {
    BOOST_MATH_STD_USING
    RealType result = 0;

    if (!_jacobi_theta_check_q(q, pol, function, &result))
        return result;
    if (!_jacobi_theta_check_z(z, pol, function, &result))
        return result;

    if (q < exp(-constants::pi<RealType>()))
        return _jacobi_theta4m1_series(z, _jacobi_theta_q_power<RealType>{q}, pol);

    return jacobi_theta4m1tau_imp(z, RealType (-log(q)/constants::pi<RealType>()), pol, function);
}

// Fourth Jacobi theta function, parameterized by q
// = 1 + 2 * SUM (-1)^n q^n^2 * cos(2nz)
template <class RealType, class Policy>
inline RealType
jacobi_theta4_imp(RealType z, RealType q, const Policy& pol, const char *function) {
    BOOST_MATH_STD_USING
    RealType result = 0;

    if (!_jacobi_theta_check_q(q, pol, function, &result))
        return result;
    if (!_jacobi_theta_check_z(z, pol, function, &result))
        return result;

    if (q < exp(-constants::pi<RealType>()))
        return RealType(1) + _jacobi_theta4m1_series(z, _jacobi_theta_q_power<RealType>{q}, pol);

    return jacobi_theta4tau_imp(z, RealType(-log(q)/constants::pi<RealType>()), pol, function);
}

// Begin public API

BOOST_MATH_EXPORT template <class T, class U, class Policy>
inline typename tools::promote_args<T, U>::type jacobi_theta1tau(T z, U tau, const Policy&) {
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename tools::promote_args<T, U>::type result_type;
   typedef typename policies::normalise<
      Policy,
      policies::promote_float<false>,
      policies::promote_double<false>,
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;

   static const char* function = "boost::math::jacobi_theta1tau<%1%>(%1%)";

   return policies::checked_narrowing_cast<result_type, Policy>(jacobi_theta1tau_imp(static_cast<result_type>(z), static_cast<result_type>(tau), forwarding_policy(), function), function);
}

BOOST_MATH_EXPORT template <class T, class U>
inline typename tools::promote_args<T, U>::type jacobi_theta1tau(T z, U tau) {
    return jacobi_theta1tau(z, tau, policies::policy<>());
}

BOOST_MATH_EXPORT template <class T, class U, class Policy>
inline typename tools::promote_args<T, U>::type jacobi_theta1(T z, U q, const Policy&) {
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename tools::promote_args<T, U>::type result_type;
   typedef typename policies::normalise<
      Policy,
      policies::promote_float<false>,
      policies::promote_double<false>,
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;

   static const char* function = "boost::math::jacobi_theta1<%1%>(%1%)";

   return policies::checked_narrowing_cast<result_type, Policy>(jacobi_theta1_imp(static_cast<result_type>(z), static_cast<result_type>(q), forwarding_policy(), function), function);
}

BOOST_MATH_EXPORT template <class T, class U>
inline typename tools::promote_args<T, U>::type jacobi_theta1(T z, U q) {
    return jacobi_theta1(z, q, policies::policy<>());
}

BOOST_MATH_EXPORT template <class T, class U, class Policy>
inline typename tools::promote_args<T, U>::type jacobi_theta2tau(T z, U tau, const Policy&) {
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename tools::promote_args<T, U>::type result_type;
   typedef typename policies::normalise<
      Policy,
      policies::promote_float<false>,
      policies::promote_double<false>,
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;

   static const char* function = "boost::math::jacobi_theta2tau<%1%>(%1%)";

   return policies::checked_narrowing_cast<result_type, Policy>(jacobi_theta2tau_imp(static_cast<result_type>(z), static_cast<result_type>(tau), forwarding_policy(), function), function);
}

BOOST_MATH_EXPORT template <class T, class U>
inline typename tools::promote_args<T, U>::type jacobi_theta2tau(T z, U tau) {
    return jacobi_theta2tau(z, tau, policies::policy<>());
}

BOOST_MATH_EXPORT template <class T, class U, class Policy>
inline typename tools::promote_args<T, U>::type jacobi_theta2(T z, U q, const Policy&) {
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename tools::promote_args<T, U>::type result_type;
   typedef typename policies::normalise<
      Policy,
      policies::promote_float<false>,
      policies::promote_double<false>,
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;

   static const char* function = "boost::math::jacobi_theta2<%1%>(%1%)";

   return policies::checked_narrowing_cast<result_type, Policy>(jacobi_theta2_imp(static_cast<result_type>(z), static_cast<result_type>(q), forwarding_policy(), function), function);
}

BOOST_MATH_EXPORT template <class T, class U>
inline typename tools::promote_args<T, U>::type jacobi_theta2(T z, U q) {
    return jacobi_theta2(z, q, policies::policy<>());
}

BOOST_MATH_EXPORT template <class T, class U, class Policy>
inline typename tools::promote_args<T, U>::type jacobi_theta3m1tau(T z, U tau, const Policy&) {
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename tools::promote_args<T, U>::type result_type;
   typedef typename policies::normalise<
      Policy,
      policies::promote_float<false>,
      policies::promote_double<false>,
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;

   static const char* function = "boost::math::jacobi_theta3m1tau<%1%>(%1%)";

   return policies::checked_narrowing_cast<result_type, Policy>(
           jacobi_theta3m1tau_imp(static_cast<result_type>(z), static_cast<result_type>(tau), forwarding_policy(), function), function);
}

BOOST_MATH_EXPORT template <class T, class U>
inline typename tools::promote_args<T, U>::type jacobi_theta3m1tau(T z, U tau) {
    return jacobi_theta3m1tau(z, tau, policies::policy<>());
}

BOOST_MATH_EXPORT template <class T, class U, class Policy>
inline typename tools::promote_args<T, U>::type jacobi_theta3tau(T z, U tau, const Policy&) {
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename tools::promote_args<T, U>::type result_type;
   typedef typename policies::normalise<
      Policy,
      policies::promote_float<false>,
      policies::promote_double<false>,
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;

   static const char* function = "boost::math::jacobi_theta3tau<%1%>(%1%)";

   return policies::checked_narrowing_cast<result_type, Policy>(jacobi_theta3tau_imp(static_cast<result_type>(z), static_cast<result_type>(tau), forwarding_policy(), function), function);
}

BOOST_MATH_EXPORT template <class T, class U>
inline typename tools::promote_args<T, U>::type jacobi_theta3tau(T z, U tau) {
    return jacobi_theta3tau(z, tau, policies::policy<>());
}


BOOST_MATH_EXPORT template <class T, class U, class Policy>
inline typename tools::promote_args<T, U>::type jacobi_theta3m1(T z, U q, const Policy&) {
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename tools::promote_args<T, U>::type result_type;
   typedef typename policies::normalise<
      Policy,
      policies::promote_float<false>,
      policies::promote_double<false>,
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;

   static const char* function = "boost::math::jacobi_theta3m1<%1%>(%1%)";

   return policies::checked_narrowing_cast<result_type, Policy>(jacobi_theta3m1_imp(static_cast<result_type>(z), static_cast<result_type>(q), forwarding_policy(), function), function);
}

BOOST_MATH_EXPORT template <class T, class U>
inline typename tools::promote_args<T, U>::type jacobi_theta3m1(T z, U q) {
    return jacobi_theta3m1(z, q, policies::policy<>());
}

BOOST_MATH_EXPORT template <class T, class U, class Policy>
inline typename tools::promote_args<T, U>::type jacobi_theta3(T z, U q, const Policy&) {
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename tools::promote_args<T, U>::type result_type;
   typedef typename policies::normalise<
      Policy,
      policies::promote_float<false>,
      policies::promote_double<false>,
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;

   static const char* function = "boost::math::jacobi_theta3<%1%>(%1%)";

   return policies::checked_narrowing_cast<result_type, Policy>(jacobi_theta3_imp(static_cast<result_type>(z), static_cast<result_type>(q), forwarding_policy(), function), function);
}

BOOST_MATH_EXPORT template <class T, class U>
inline typename tools::promote_args<T, U>::type jacobi_theta3(T z, U q) {
    return jacobi_theta3(z, q, policies::policy<>());
}

BOOST_MATH_EXPORT template <class T, class U, class Policy>
inline typename tools::promote_args<T, U>::type jacobi_theta4m1tau(T z, U tau, const Policy&) {
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename tools::promote_args<T, U>::type result_type;
   typedef typename policies::normalise<
      Policy,
      policies::promote_float<false>,
      policies::promote_double<false>,
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;

   static const char* function = "boost::math::jacobi_theta4m1tau<%1%>(%1%)";

   return policies::checked_narrowing_cast<result_type, Policy>(jacobi_theta4m1tau_imp(static_cast<result_type>(z), static_cast<result_type>(tau), forwarding_policy(), function), function);
}

BOOST_MATH_EXPORT template <class T, class U>
inline typename tools::promote_args<T, U>::type jacobi_theta4m1tau(T z, U tau) {
    return jacobi_theta4m1tau(z, tau, policies::policy<>());
}

BOOST_MATH_EXPORT template <class T, class U, class Policy>
inline typename tools::promote_args<T, U>::type jacobi_theta4tau(T z, U tau, const Policy&) {
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename tools::promote_args<T, U>::type result_type;
   typedef typename policies::normalise<
      Policy,
      policies::promote_float<false>,
      policies::promote_double<false>,
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;

   static const char* function = "boost::math::jacobi_theta4tau<%1%>(%1%)";

   return policies::checked_narrowing_cast<result_type, Policy>(jacobi_theta4tau_imp(static_cast<result_type>(z), static_cast<result_type>(tau), forwarding_policy(), function), function);
}

BOOST_MATH_EXPORT template <class T, class U>
inline typename tools::promote_args<T, U>::type jacobi_theta4tau(T z, U tau) {
    return jacobi_theta4tau(z, tau, policies::policy<>());
}

BOOST_MATH_EXPORT template <class T, class U, class Policy>
inline typename tools::promote_args<T, U>::type jacobi_theta4m1(T z, U q, const Policy&) {
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename tools::promote_args<T, U>::type result_type;
   typedef typename policies::normalise<
      Policy,
      policies::promote_float<false>,
      policies::promote_double<false>,
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;

   static const char* function = "boost::math::jacobi_theta4m1<%1%>(%1%)";

   return policies::checked_narrowing_cast<result_type, Policy>(jacobi_theta4m1_imp(static_cast<result_type>(z), static_cast<result_type>(q), forwarding_policy(), function), function);
}

BOOST_MATH_EXPORT template <class T, class U>
inline typename tools::promote_args<T, U>::type jacobi_theta4m1(T z, U q) {
    return jacobi_theta4m1(z, q, policies::policy<>());
}

BOOST_MATH_EXPORT template <class T, class U, class Policy>
inline typename tools::promote_args<T, U>::type jacobi_theta4(T z, U q, const Policy&) {
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename tools::promote_args<T, U>::type result_type;
   typedef typename policies::normalise<
      Policy,
      policies::promote_float<false>,
      policies::promote_double<false>,
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;

   static const char* function = "boost::math::jacobi_theta4<%1%>(%1%)";

   return policies::checked_narrowing_cast<result_type, Policy>(jacobi_theta4_imp(static_cast<result_type>(z), static_cast<result_type>(q), forwarding_policy(), function), function);
}

BOOST_MATH_EXPORT template <class T, class U>
inline typename tools::promote_args<T, U>::type jacobi_theta4(T z, U q) {
    return jacobi_theta4(z, q, policies::policy<>());
}

}}

#endif
