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

#include <cmath>
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

// When the caller supplies tau, every term is an exponential exp(-E) whose
// argument E is a product or quotient of tau, pi and a small integer (or, for
// tau < 1, the square of z plus a multiple of pi/2). Rounding E to working
// precision costs |E| ulps in the exponential, so instead the rounding
// errors of pi, of the products and of the division are tracked exactly with
// Dekker's product and Knuth's sum (which work for any binary floating-point
// type) and applied as a correction factor exp(-dE) = 1 - dE. What remains
// is the rounding of tau itself, which is the caller's.
// (Cached per type: for a type whose precision can change at run time the
// cached value may split at the wrong place, in which case the corrections
// below are merely inexact, i.e. no worse than not applying them.)
template <class RealType>
inline RealType _jacobi_theta_splitter() {
    BOOST_MATH_STD_USING
    static const RealType splitter = ldexp(RealType(1), (tools::digits<RealType>() + 1) / 2) + 1;
    return splitter;
}

// Veltkamp's splitting: t = hi + lo where hi holds only the leading half of
// the significand, so that products of hi parts are exact.
template <class RealType>
inline void _jacobi_theta_split(RealType t, RealType splitter, RealType& hi, RealType& lo) {
    RealType g = splitter * t;
    hi = g - (g - t);
    lo = t - hi;
}

// The exact product a*b is ab + (return value), where ab = fl(a*b). Where
// the standard library advertises a fast fused multiply-add it is used
// instead, since fma(a, b, -ab) is exactly this quantity.
template <class RealType>
inline RealType _jacobi_theta_product_error(RealType a, RealType b, RealType ab, RealType splitter) {
    RealType ah, al, bh, bl;
    _jacobi_theta_split(a, splitter, ah, al);
    _jacobi_theta_split(b, splitter, bh, bl);
    return (((ah * bh - ab) + ah * bl) + al * bh) + al * bl;
}
#ifdef FP_FAST_FMAF
inline float _jacobi_theta_product_error(float a, float b, float ab, float) {
    return std::fma(a, b, -ab);
}
#endif
#ifdef FP_FAST_FMA
inline double _jacobi_theta_product_error(double a, double b, double ab, double) {
    return std::fma(a, b, -ab);
}
#endif
#ifdef FP_FAST_FMAL
inline long double _jacobi_theta_product_error(long double a, long double b, long double ab, long double) {
    return std::fma(a, b, -ab);
}
#endif

// The exact sum a+b is s + (return value), where s = fl(a+b).
template <class RealType>
inline RealType _jacobi_theta_sum_error(RealType a, RealType b, RealType s) {
    RealType bb = s - a;
    return (a - (s - bb)) + (b - bb);
}

// Pi as hi + lo: hi is exactly representable and lo completes it to about
// 120 bits, as five pieces of 24 significant bits each. Types with more
// precision than that use their own rounding of pi with lo = 0.
template <class RealType>
inline RealType _jacobi_theta_pi(RealType& lo) {
    BOOST_MATH_STD_USING
    if (tools::digits<RealType>() > 116) {
        lo = 0;
        return constants::pi<RealType>();
    }
    static const RealType pi_lo = ldexp(RealType(10625384), -46) + ldexp(RealType(12727492), -70) + ldexp(RealType(13001355), -94) + ldexp(RealType(8444956), -118);
    static const RealType pi_hi = ldexp(RealType(13176794), -22);
    lo = pi_lo;
    return pi_hi;
}

// The quantity everything is expressed in is a = pi*tau = -ln(q), held as
// a_hi + a_lo. Built from tau it is exact (to second order); built from q it
// carries the rounding of the logarithm, which is inherent to that
// parameterization, but nothing else.
template <class RealType>
struct _jacobi_theta_exponents {
    RealType tau;           // only used for the scale factor 1/sqrt(tau)
    RealType splitter;
    RealType pi_hi, pi_lo;  // pi = pi_hi + pi_lo
    RealType a_hi, a_lo;    // pi * tau = -ln(q) = a_hi + a_lo
    RealType inv_a_hi;      // 1 / a_hi, only needed when tau < 1

    static _jacobi_theta_exponents from_tau(RealType tau) {
        _jacobi_theta_exponents x;
        x.tau = tau;
        x.a_hi = x.pi_times(tau, x.a_lo);
        x.inv_a_hi = (tau < 1) ? 1 / x.a_hi : RealType(0);
        return x;
    }

    static _jacobi_theta_exponents from_nome(RealType q) {
        BOOST_MATH_STD_USING
        _jacobi_theta_exponents x;
        x.a_hi = -log(q);
        x.a_lo = 0;
        x.tau = x.a_hi / constants::pi<RealType>();
        x.inv_a_hi = (x.tau < 1) ? 1 / x.a_hi : RealType(0);
        return x;
    }

    // pi * e as hi + lo for a small exact e
    RealType pi_times(RealType e, RealType& lo) const {
        RealType p = pi_hi * e;
        RealType q = pi_lo * e;
        RealType hi = p + q;
        // |q| << |p|, so the sum's error is simply q - (hi - p)
        lo = (q - (hi - p)) + _jacobi_theta_product_error(pi_hi, e, p, splitter);
        return hi;
    }

    // (N + dN) / (a_hi + a_lo) as E + dE, to second order
    RealType divide_by_a(RealType N, RealType dN, RealType& dE) const {
        RealType E = N * inv_a_hi;
        RealType Ea = E * a_hi;
        RealType r = (N - Ea) - _jacobi_theta_product_error(E, a_hi, Ea, splitter);
        dE = (r + dN - E * a_lo) * inv_a_hi;
        return E;
    }

    // exp(-a * e) = q^e for a small exact e such as n^2 or (n+1/2)^2
    RealType direct(RealType e) const {
        BOOST_MATH_STD_USING
        RealType E = a_hi * e;
        RealType result = exp(-E);
        if (result == 0)
            return result;
        RealType dE = _jacobi_theta_product_error(a_hi, e, E, splitter) + a_lo * e;
        return result * (1 - dE);
    }

    // exp(-pi * e / tau) = exp(-pi^2 * e / a), the nome of 1/tau raised to e
    RealType inverse(RealType e) const {
        BOOST_MATH_STD_USING
        RealType dN, dP;
        RealType N = pi_times(e, dN);
        RealType P = pi_times(N, dP);
        dP += pi_hi * dN;
        RealType dE;
        RealType E = divide_by_a(P, dP, dE);
        RealType result = exp(-E);
        if (result == 0)
            return result;
        return result * (1 - dE);
    }

    // Reduces z by the nearest multiple of the period (half_pis * pi/2),
    // returning the remainder along with its error dz (the remainder is only
    // exact with respect to the rounded pi, and the Gaussians below amplify
    // that error by 2z/(pi*tau)), and the multiple k for the caller's sign.
    RealType reduce(RealType z, int half_pis, RealType& dz, RealType& k) const {
        BOOST_MATH_STD_USING
        RealType period = pi_hi * RealType(half_pis) / 2;
        k = floor(z / period + RealType(0.5));
        RealType dc;
        RealType c = pi_times(k * RealType(half_pis) / 2, dc);
        RealType r = z - c;
        dz = _jacobi_theta_sum_error(z, RealType(-c), r) - dc;
        return r;
    }

    // exp(-(z + m*pi/2)^2 / (pi * tau)) for an integer m, where dz is the
    // error in z from the argument reduction
    RealType gaussian(RealType z, RealType dz, int m) const {
        BOOST_MATH_STD_USING
        RealType z_n = z;
        if (m != 0) {
            RealType dc;
            RealType c = pi_times(RealType(m), dc) / 2;
            z_n = z + c;
            dz += _jacobi_theta_sum_error(z, c, z_n) + dc / 2;
        }
        RealType s = z_n * z_n;
        RealType ds = _jacobi_theta_product_error(z_n, z_n, s, splitter) + 2 * z_n * dz;
        RealType dE;
        RealType E = divide_by_a(s, ds, dE);
        RealType result = exp(-E);
        if (result == 0)
            return result;
        return result * (1 - dE);
    }

    // expm1(-2 * z * k / tau) = expm1(-2 * pi * z * k / a) for a small exact
    // k, where dz is the error in z
    RealType expm1_scaled(RealType z, RealType dz, RealType k) const {
        BOOST_MATH_STD_USING
        RealType P = z * k;
        RealType dP = _jacobi_theta_product_error(z, k, P, splitter) + dz * k;
        RealType dN;
        RealType N = pi_times(P, dN);
        dN += pi_hi * dP;
        RealType dQ;
        RealType Q = divide_by_a(N, dN, dQ);
        RealType result = boost::math::expm1(RealType(-2 * Q));
        if (result == -1)
            return result;
        return result - 2 * dQ * (1 + result);
    }

private:
    _jacobi_theta_exponents() : splitter(_jacobi_theta_splitter<RealType>()) {
        pi_hi = _jacobi_theta_pi(pi_lo);
    }
};

template <class RealType>
struct _jacobi_theta_tau_power {
    const _jacobi_theta_exponents<RealType>& x;
    RealType operator()(RealType exponent) const {
        return x.direct(exponent);
    }
};

// exp(-pi * e / tau): the nome of 1/tau, used by the modular transformation
// at z = 0 without forming the rounded reciprocal.
template <class RealType>
struct _jacobi_theta_inverse_tau_power {
    const _jacobi_theta_exponents<RealType>& x;
    RealType operator()(RealType exponent) const {
        return x.inverse(exponent);
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

// SUM exp(-(z + m*pi/2)^2 / (pi*tau)) over m = m0, m0 + m_step, ... until the
// terms are negligible.
template <class RealType>
inline RealType
_jacobi_theta_sum(const _jacobi_theta_exponents<RealType>& x, RealType z, RealType dz, int m, int m_step, RealType eps) {
    RealType delta = 0, partial_result = 0;
    RealType last_delta = 0;

    do {
        last_delta = delta;
        delta = x.gaussian(z, dz, m);
        partial_result += delta;
        m += m_step;
    } while (!_jacobi_theta_converged(last_delta, delta, eps));

    return partial_result;
}

// The following _IMAGINARY theta functions assume imaginary z and are for
// internal use only. They are designed to increase accuracy and reduce the
// number of iterations required for convergence for large |q|. The z argument
// is scaled by 1/tau, and the summations are rewritten to be double-sided
// following DLMF 20.13.4 and 20.13.5. Each term is a Gaussian
// exp(-(z - c)^2/(Pi*tau)) centered at a multiple of Pi or Pi/2, and the
// results are scaled by 1/sqrt(tau).
//
// These functions are triggered when tau < 1, i.e. |q| > exp(-Pi) = 0.043
//
// Note that jacobi_theta4 uses the imaginary version of jacobi_theta2 (and
// vice-versa). jacobi_theta1 and jacobi_theta3 use the imaginary versions of
// themselves, following DLMF 20.7.30 - 20.7.33.

// theta1(z|tau) = 1/sqrt(tau) * SUM_{n>=0} (-1)^n [G(z - c_n) - G(z + c_n)]
// with c_n = Pi*(n+1/2) and G(x) = exp(-x^2/(Pi*tau)).
//
// Each bracket is a difference of two Gaussians which nearly cancel when z is
// small, so it is evaluated instead through the exact identity
//   G(z - c) - G(z + c) = -G(z - c) * expm1(-2*z*(2n+1)/tau),
// which keeps full relative precision all the way down to z -> 0.
// Requires 0 <= z <= Pi/2; the caller reduces z into this range, and
// passes the error dz of the reduced z.
template <class RealType, class Policy>
inline RealType
_IMAGINARY_jacobi_theta1tau(RealType z, RealType dz, const _jacobi_theta_exponents<RealType>& x, const Policy&) {
    BOOST_MATH_STD_USING
    RealType eps = policies::get_epsilon<RealType, Policy>();
    RealType result = 0, g = 0, last_g, pair;
    unsigned n = 0;

    do {
        last_g = g;
        g = x.gaussian(z, dz, -static_cast<int>(2*n + 1));
        pair = -g * x.expm1_scaled(z, dz, RealType(2*n + 1));
        if (n%2)
            pair = -pair;

        result += pair;
        n++;
    } while (!_jacobi_theta_converged(last_g, g, eps));

    return result / sqrt(x.tau);
}

template <class RealType, class Policy>
inline RealType
_IMAGINARY_jacobi_theta2tau(RealType z, RealType dz, const _jacobi_theta_exponents<RealType>& x, const Policy&) {
    BOOST_MATH_STD_USING
    RealType eps = policies::get_epsilon<RealType, Policy>();
    RealType result = RealType(0);

    // n>=0: centers at z + Pi/2 + n*Pi
    result += _jacobi_theta_sum(x, z, dz, 1, 2, eps);
    // n<0
    result += _jacobi_theta_sum(x, z, dz, -1, -2, eps);

    return result / sqrt(x.tau);
}

template <class RealType, class Policy>
inline RealType
_IMAGINARY_jacobi_theta3tau(RealType z, RealType dz, const _jacobi_theta_exponents<RealType>& x, const Policy&) {
    BOOST_MATH_STD_USING
    RealType eps = policies::get_epsilon<RealType, Policy>();
    RealType result = 0;

    // n=0
    result += x.gaussian(z, dz, 0);
    // n>0: centers at z + n*Pi
    result += _jacobi_theta_sum(x, z, dz, 2, 2, eps);
    // n<0
    result += _jacobi_theta_sum(x, z, dz, -2, -2, eps);

    return result / sqrt(x.tau);
}

template <class RealType, class Policy>
inline RealType
_IMAGINARY_jacobi_theta4tau(RealType z, RealType dz, const _jacobi_theta_exponents<RealType>& x, const Policy&) {
    BOOST_MATH_STD_USING
    RealType eps = policies::get_epsilon<RealType, Policy>();
    RealType result = 0;

    // n = 0
    result += x.gaussian(z, dz, 0);

    // n > 0 odd: centers at z + Pi + 2n*Pi
    result -= _jacobi_theta_sum(x, z, dz, 2, 4, eps);
    // n < 0 odd
    result -= _jacobi_theta_sum(x, z, dz, -2, -4, eps);
    // n > 0 even: centers at z + 2*Pi + 2n*Pi
    result += _jacobi_theta_sum(x, z, dz, 4, 4, eps);
    // n < 0 even
    result += _jacobi_theta_sum(x, z, dz, -4, -4, eps);

    return result / sqrt(x.tau);
}

// Dispatch on the size of tau (i.e. of the nome): the direct Fourier series
// for tau >= 1, otherwise the modular transformation to 1/tau. At z = 0 the
// transformed series is evaluated directly (single-sided, with the nome of
// 1/tau formed without rounding the reciprocal); otherwise as double-sided
// Gaussian sums.

// = 2 * SUM (-1)^n * exp(i*Pi*Tau*(n+1/2)^2) * sin((2n+1)z)
template <class RealType, class Policy>
inline RealType
_jacobi_theta1_dispatch(RealType z, const _jacobi_theta_exponents<RealType>& x, const Policy& pol) {
    BOOST_MATH_STD_USING
    if (x.tau < 1.0) {
        // Reduce to -Pi/2 <= z <= Pi/2 using theta1(z + Pi) = -theta1(z)...
        RealType dz, k;
        z = x.reduce(z, 2, dz, k);
        RealType sign = (fmod(k, RealType(2)) == 0) ? 1 : -1;
        // ...and then to 0 <= z <= Pi/2 since theta1 is odd.
        if (z < 0) {
            z = -z;
            dz = -dz;
            sign = -sign;
        }
        return sign * _IMAGINARY_jacobi_theta1tau(z, dz, x, pol);
    }
    return _jacobi_theta1_series(z, _jacobi_theta_tau_power<RealType>{x}, pol);
}

// = 2 * SUM exp(i*Pi*Tau*(n+1/2)^2) * cos((2n+1)z)
template <class RealType, class Policy>
inline RealType
_jacobi_theta2_dispatch(RealType z, const _jacobi_theta_exponents<RealType>& x, const Policy& pol) {
    BOOST_MATH_STD_USING
    if (x.tau < 1.0 && abs(z) == 0.0) { // theta4(0|1/tau)/sqrt(tau)
        return (RealType(1) + _jacobi_theta4m1_series(z, _jacobi_theta_inverse_tau_power<RealType>{x}, pol)) / sqrt(x.tau);
    } else if (x.tau < 1.0) { // DLMF 20.7.31
        // Reduce to -Pi <= z <= Pi (theta2 has period 2*Pi)
        RealType dz, k;
        z = x.reduce(z, 4, dz, k);
        return _IMAGINARY_jacobi_theta4tau(z, dz, x, pol);
    }
    return _jacobi_theta2_series(z, _jacobi_theta_tau_power<RealType>{x}, pol);
}

// = 1 + 2 * SUM exp(i*Pi*Tau*(n)^2) * cos(2nz)
template <class RealType, class Policy>
inline RealType
_jacobi_theta3_dispatch(RealType z, const _jacobi_theta_exponents<RealType>& x, const Policy& pol) {
    BOOST_MATH_STD_USING
    if (x.tau < 1.0 && abs(z) == 0.0) { // theta3(0|1/tau)/sqrt(tau)
        return (RealType(1) + _jacobi_theta3m1_series(z, _jacobi_theta_inverse_tau_power<RealType>{x}, pol)) / sqrt(x.tau);
    } else if (x.tau < 1.0) { // DLMF 20.7.32
        // Reduce to -Pi/2 <= z <= Pi/2 (theta3 has period Pi)
        RealType dz, k;
        z = x.reduce(z, 2, dz, k);
        return _IMAGINARY_jacobi_theta3tau(z, dz, x, pol);
    }
    return RealType(1) + _jacobi_theta3m1_series(z, _jacobi_theta_tau_power<RealType>{x}, pol);
}

// = 2 * SUM exp(i*Pi*Tau*(n)^2) * cos(2nz), n >= 1 (theta3 minus one)
// This preserves accuracy for small values of q (i.e. tau > 1). For larger
// values of q, the minus one version usually won't help.
template <class RealType, class Policy>
inline RealType
_jacobi_theta3m1_dispatch(RealType z, const _jacobi_theta_exponents<RealType>& x, const Policy& pol) {
    if (x.tau < 1.0)
        return _jacobi_theta3_dispatch(z, x, pol) - RealType(1);
    return _jacobi_theta3m1_series(z, _jacobi_theta_tau_power<RealType>{x}, pol);
}

// = 1 + 2 * SUM (-1)^n exp(i*Pi*Tau*(n)^2) * cos(2nz)
template <class RealType, class Policy>
inline RealType
_jacobi_theta4_dispatch(RealType z, const _jacobi_theta_exponents<RealType>& x, const Policy& pol) {
    BOOST_MATH_STD_USING
    if (x.tau < 1.0 && abs(z) == 0.0) { // theta2(0|1/tau)/sqrt(tau)
        return _jacobi_theta2_series(z, _jacobi_theta_inverse_tau_power<RealType>{x}, pol) / sqrt(x.tau);
    } else if (x.tau < 1.0) { // DLMF 20.7.33
        // Reduce to -Pi/2 <= z <= Pi/2 (theta4 has period Pi)
        RealType dz, k;
        z = x.reduce(z, 2, dz, k);
        return _IMAGINARY_jacobi_theta2tau(z, dz, x, pol);
    }
    return RealType(1) + _jacobi_theta4m1_series(z, _jacobi_theta_tau_power<RealType>{x}, pol);
}

// = 2 * SUM (-1)^n exp(i*Pi*Tau*(n)^2) * cos(2nz), n >= 1 (theta4 minus one)
// This preserves accuracy for small values of q (i.e. tau > 1).
template <class RealType, class Policy>
inline RealType
_jacobi_theta4m1_dispatch(RealType z, const _jacobi_theta_exponents<RealType>& x, const Policy& pol) {
    if (x.tau < 1.0)
        return _jacobi_theta4_dispatch(z, x, pol) - RealType(1);
    return _jacobi_theta4m1_series(z, _jacobi_theta_tau_power<RealType>{x}, pol);
}

// The twelve _imp functions below validate their arguments and then hand
// over to the dispatchers above. The q versions use the direct series with
// pow() when q < exp(-Pi), and otherwise go through a = -ln(q).

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

    return _jacobi_theta1_dispatch(z, _jacobi_theta_exponents<RealType>::from_tau(tau), pol);
}

template <class RealType, class Policy>
inline RealType
jacobi_theta1_imp(RealType z, RealType q, const Policy& pol, const char *function) {
    BOOST_MATH_STD_USING
    RealType result = 0;

    if (!_jacobi_theta_check_q(q, pol, function, &result))
        return result;
    if (!_jacobi_theta_check_z(z, pol, function, &result))
        return result;
    if (abs(z) == 0.0)
        return result;

    if (q < exp(-constants::pi<RealType>()))
        return _jacobi_theta1_series(z, _jacobi_theta_q_power<RealType>{q}, pol);

    return _jacobi_theta1_dispatch(z, _jacobi_theta_exponents<RealType>::from_nome(q), pol);
}

template <class RealType, class Policy>
inline RealType
jacobi_theta2tau_imp(RealType z, RealType tau, const Policy& pol, const char *function)
{
    RealType result = 0;

    if (!_jacobi_theta_check_tau(tau, pol, function, &result))
        return result;
    if (!_jacobi_theta_check_z(z, pol, function, &result))
        return result;

    return _jacobi_theta2_dispatch(z, _jacobi_theta_exponents<RealType>::from_tau(tau), pol);
}

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

    return _jacobi_theta2_dispatch(z, _jacobi_theta_exponents<RealType>::from_nome(q), pol);
}

template <class RealType, class Policy>
inline RealType
jacobi_theta3tau_imp(RealType z, RealType tau, const Policy& pol, const char *function)
{
    RealType result = 0;

    if (!_jacobi_theta_check_tau(tau, pol, function, &result))
        return result;
    if (!_jacobi_theta_check_z(z, pol, function, &result))
        return result;

    return _jacobi_theta3_dispatch(z, _jacobi_theta_exponents<RealType>::from_tau(tau), pol);
}

template <class RealType, class Policy>
inline RealType
jacobi_theta3m1tau_imp(RealType z, RealType tau, const Policy& pol, const char *function)
{
    RealType result = 0;

    if (!_jacobi_theta_check_tau(tau, pol, function, &result))
        return result;
    if (!_jacobi_theta_check_z(z, pol, function, &result))
        return result;

    return _jacobi_theta3m1_dispatch(z, _jacobi_theta_exponents<RealType>::from_tau(tau), pol);
}

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

    return _jacobi_theta3m1_dispatch(z, _jacobi_theta_exponents<RealType>::from_nome(q), pol);
}

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

    return _jacobi_theta3_dispatch(z, _jacobi_theta_exponents<RealType>::from_nome(q), pol);
}

template <class RealType, class Policy>
inline RealType
jacobi_theta4tau_imp(RealType z, RealType tau, const Policy& pol, const char *function)
{
    RealType result = 0;

    if (!_jacobi_theta_check_tau(tau, pol, function, &result))
        return result;
    if (!_jacobi_theta_check_z(z, pol, function, &result))
        return result;

    return _jacobi_theta4_dispatch(z, _jacobi_theta_exponents<RealType>::from_tau(tau), pol);
}

template <class RealType, class Policy>
inline RealType
jacobi_theta4m1tau_imp(RealType z, RealType tau, const Policy& pol, const char *function)
{
    RealType result = 0;

    if (!_jacobi_theta_check_tau(tau, pol, function, &result))
        return result;
    if (!_jacobi_theta_check_z(z, pol, function, &result))
        return result;

    return _jacobi_theta4m1_dispatch(z, _jacobi_theta_exponents<RealType>::from_tau(tau), pol);
}

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

    return _jacobi_theta4m1_dispatch(z, _jacobi_theta_exponents<RealType>::from_nome(q), pol);
}

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

    return _jacobi_theta4_dispatch(z, _jacobi_theta_exponents<RealType>::from_nome(q), pol);
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
