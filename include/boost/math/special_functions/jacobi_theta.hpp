// Jacobi theta functions
// Copyright Evan Miller 2020
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
// series. Using the q parameterization, and summing from n = 0 to ∞:
//
// θ₁(z,q) = 2 Σ (-1)ⁿ * q^(n+1/2)² * sin((2n+1)z)
// θ₂(z,q) = 2 Σ q^(n+1/2)² * cos((2n+1)z)
// θ₃(z,q) = 1 + 2 Σ q^n² * cos(2nz)
// θ₄(z,q) = 1 + 2 Σ q^n² * cos(2nz)
//
// Appropriately multiplied and divided, these four theta functions can be used
// to implement the famous Jacabi elliptic functions - but this is not really
// recommended, as the existing Boost implementations are likely faster and
// more accurate.  More saliently, setting z = 0 on the fourth theta function
// will produce the limiting CDF of the Kolmogorov-Smirnov distribution, which
// is this particular implementation’s raison d'être.
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
// Mathematically, q and τ are related by:
//
// q = exp(iπτ)
//
// However, the τ in the equation above is *not* identical to the tau in the function
// signature. Instead, `tau` is the imaginary component of τ. Mathematically, τ can
// be complex - but practically, most applications call for a purely imaginary τ.
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
// A third quartet of functions are provided for improving accuracy in cases
// where q is small, e.g. |q| < exp(-π) ≅ 0.04. In this domain of q values,
// the third and fourth theta functions always return values close to 1. So
// the following "m1" functions are provided, similar in spirit to `expm1`,
// which return one less than their regular counterparts:
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

// In order to preserve accuracy with large q, __JACOBI_THETA_USE_IMAGINARY
// will switch over to the imaginary version of functions for |q| > 0.85. This
// cuts down on the number of required iterations, and improves accuracy for
// large q, but comes at the cost of precision. Maybe someone smarter than me
// can fix the precision issues with the orgy of exponentials.
#define __JACOBI_THETA_USE_IMAGINARY

namespace boost{ namespace math{

// Simple functions - parameterized by q
template <class T>
inline T jacobi_theta1(T z, T q);
template <class T>
inline T jacobi_theta2(T z, T q);
template <class T>
inline T jacobi_theta3(T z, T q);
template <class T>
inline T jacobi_theta4(T z, T q);

// Simple functions - parameterized by tau (assumed imaginary)
// q = exp(iπτ)
// tau = -log(q)/π
template <class T>
inline T jacobi_theta1tau(T z, T tau);
template <class T>
inline T jacobi_theta2tau(T z, T tau);
template <class T>
inline T jacobi_theta3tau(T z, T tau);
template <class T>
inline T jacobi_theta4tau(T z, T tau);

// Minus one versions for small q / large tau
template <class T>
inline T jacobi_theta3m1(T z, T q);
template <class T>
inline T jacobi_theta4m1(T z, T q);
template <class T>
inline T jacobi_theta3m1tau(T z, T tau);
template <class T>
inline T jacobi_theta4m1tau(T z, T tau);

// Policied versions - parameterized by q
template <class RealType, class Policy>
inline RealType jacobi_theta1(RealType z, RealType q, const Policy& pol);
template <class RealType, class Policy>
inline RealType jacobi_theta2(RealType z, RealType q, const Policy& pol);
template <class RealType, class Policy>
inline RealType jacobi_theta3(RealType z, RealType q, const Policy& pol);
template <class RealType, class Policy>
inline RealType jacobi_theta4(RealType z, RealType q, const Policy& pol);

// Policied versions - parameterized by tau
template <class RealType, class Policy>
inline RealType jacobi_theta1tau(RealType z, RealType tau, const Policy& pol);
template <class RealType, class Policy>
inline RealType jacobi_theta2tau(RealType z, RealType tau, const Policy& pol);
template <class RealType, class Policy>
inline RealType jacobi_theta3tau(RealType z, RealType tau, const Policy& pol);
template <class RealType, class Policy>
inline RealType jacobi_theta4tau(RealType z, RealType tau, const Policy& pol);

// Policied m1 functions
template <class RealType, class Policy>
inline RealType jacobi_theta3m1(RealType z, RealType q, const Policy& pol);
template <class RealType, class Policy>
inline RealType jacobi_theta4m1(RealType z, RealType q, const Policy& pol);
template <class RealType, class Policy>
inline RealType jacobi_theta3m1tau(RealType z, RealType tau, const Policy& pol);
template <class RealType, class Policy>
inline RealType jacobi_theta4m1tau(RealType z, RealType tau, const Policy& pol);

#ifdef __JACOBI_THETA_USE_IMAGINARY
// The following _IMAGINARY theta functions assume imaginary z and are for
// internal use only. The z argument is scaled by tau, and sines and cosines
// are replaced with hyperbolic sines and cosines to accommodate the imaginary
// argument. (Recall sinh=(exp(x)-exp(-x))/2 and cosh=(exp(x)+exp(-x))/2)
//
// The return values are scaled by exp(-tau*z²/π)/sqrt(tau).
//
// Using the convention τ'=-1/τ and noting that both τ and τ' are always
// imaginary (at least in our applications), the _IMAGINARY functions are used
// to implement the following four relations:
//
// [20.7.30] sqrt(-iτ)θ₁(z|τ) = -i*exp(iτ'z²/π)*θ₁(zτ'|τ')
// [20.7.31] sqrt(-iτ)θ₂(z|τ) =    exp(iτ'z²/π)*θ₄(zτ'|τ')
// [20.7.32] sqrt(-iτ)θ₃(z|τ) =    exp(iτ'z²/π)*θ₃(zτ'|τ')
// [20.7.33] sqrt(-iτ)θ₄(z|τ) =    exp(iτ'z²/π)*θ₂(zτ'|τ')
template <class RealType, class Policy>
inline RealType
_IMAGINARY_jacobi_theta1tau(RealType z, RealType tau, const Policy& pol) {
    BOOST_MATH_STD_USING
    unsigned n = 0;
    RealType eps = policies::get_epsilon<RealType, Policy>();
    RealType /* q_n, */ delta, delta1, delta2, result = RealType(0);

    do {
        delta1 = exp(tau*(z*(RealType(2*n+1) - z/constants::pi<RealType>())
                    -constants::pi<RealType>() * RealType(n + 0.5)*RealType(n + 0.5)));
        // q_n = exp(-tau * constants::pi<RealType>() * RealType(n + 0.5)*RealType(n + 0.5) );
        delta2 = exp(tau*(z*(-RealType(2*n+1) - z/constants::pi<RealType>())
                    -constants::pi<RealType>() * RealType(n + 0.5)*RealType(n + 0.5)));

        if (n%2) {
            delta = delta2 - delta1;
        } else {
            delta = delta1 - delta2;
        }

        result += delta;
        n++;
    } while (abs(delta * sqrt(tau)) > eps || (abs(result * sqrt(tau)) > eps * eps && abs(delta/result) > eps));
    
    if (abs(result) < eps * eps)
        return RealType(0);

    return result * sqrt(tau);
}

template <class RealType, class Policy>
inline RealType
_IMAGINARY_jacobi_theta2tau(RealType z, RealType tau, const Policy& pol) {
    BOOST_MATH_STD_USING
    unsigned n = 0;
    RealType eps = policies::get_epsilon<RealType, Policy>();
    RealType /* q_n, */ delta, delta1, delta2, result = RealType(0);

    do {
        delta1 = exp(tau*(z*(RealType(2*n+1) - z/constants::pi<RealType>())
                    -constants::pi<RealType>() * RealType(n + 0.5)*RealType(n + 0.5)));
        // q_n = exp(-tau * constants::pi<RealType>() * RealType(n + 0.5)*RealType(n + 0.5) );
        delta2 = exp(tau*(z*(-RealType(2*n+1) - z/constants::pi<RealType>())
                    -constants::pi<RealType>() * RealType(n + 0.5)*RealType(n + 0.5)));

        delta = delta1 + delta2;

        result += delta;
        n++;
    } while (abs(delta * sqrt(tau)) > eps || (abs(result * sqrt(tau)) > eps * eps && abs(delta/result) > eps));
    
    if (abs(result) < eps * eps)
        return RealType(0);

    return result * sqrt(tau);
}

template <class RealType, class Policy>
inline RealType
_IMAGINARY_jacobi_theta3tau(RealType z, RealType tau, const Policy& pol) {
    BOOST_MATH_STD_USING
    unsigned n = 1;
    RealType eps = policies::get_epsilon<RealType, Policy>();
    RealType /* q_n, */ delta, delta1, delta2, result = exp(-z*z*tau/constants::pi<RealType>());

    do {
        delta1 = exp(tau*(z*(RealType(2*n) - z/constants::pi<RealType>())
                    -constants::pi<RealType>() * RealType(n)*RealType(n)));
        // q_n = exp(-tau * constants::pi<RealType>() * RealType(n + 0.5)*RealType(n + 0.5) );
        delta2 = exp(tau*(z*(-RealType(2*n) - z/constants::pi<RealType>())
                    -constants::pi<RealType>() * RealType(n)*RealType(n)));

        delta = delta1 + delta2;

        result += delta;
        n++;
    } while (abs(delta * sqrt(tau)) > eps || (abs(result * sqrt(tau)) > eps * eps && abs(delta/result) > eps));
    
    if (abs(result) < eps * eps)
        return RealType(0);

    return result * sqrt(tau);
}

template <class RealType, class Policy>
inline RealType
_IMAGINARY_jacobi_theta4tau(RealType z, RealType tau, const Policy& pol) {
    BOOST_MATH_STD_USING
    unsigned n = 1;
    RealType eps = policies::get_epsilon<RealType, Policy>();
    RealType /* q_n, */ delta, delta1, delta2, result = exp(-z*z*tau/constants::pi<RealType>());

    do {
        delta1 = exp(tau*(z*(RealType(2*n) - z/constants::pi<RealType>())
                -constants::pi<RealType>() * RealType(n)*RealType(n)));
        // q_n = exp(-tau * constants::pi<RealType>() * RealType(n + 0.5)*RealType(n + 0.5) );
        delta2 = exp(tau*(z*(-RealType(2*n) - z/constants::pi<RealType>())
                -constants::pi<RealType>() * RealType(n)*RealType(n)));

        if (n%2) {
            delta = -(delta1 + delta2);
        } else {
            delta = delta1 + delta2;
        }

        result += delta;
        n++;
    } while (abs(delta * sqrt(tau)) > eps || (abs(result * sqrt(tau)) > eps * eps && abs(delta/result) > eps));
    
    if (abs(result) < eps * eps)
        return RealType(0);

    return result * sqrt(tau);
}
#endif

// First Jacobi theta function (Parameterized by tau - assumed imaginary)
// = 2 * Σ (-1)^n * exp(iπτ*(n+1/2)^2) * sin((2n+1)z)
template <class RealType, class Policy>
inline RealType
jacobi_theta1tau(RealType z, RealType tau, const Policy& pol)
{
    BOOST_MATH_STD_USING
    unsigned n = 0;
    RealType eps = policies::get_epsilon<RealType, Policy>();
    RealType q_n, delta, result = RealType(0);

    if (tau <= 0.0)
        return policies::raise_domain_error<RealType>("boost::math::jacobi_theta1tau<%1%>(%1%)",
                "tau must be greater than 0 but got %1%.", tau, pol);

    if (abs(z) == 0.0)
        return result;

#ifdef __JACOBI_THETA_USE_IMAGINARY
    if (exp(-tau*constants::pi<RealType>()) > 0.85) {
        z = fmod(z, constants::two_pi<RealType>());
        while (z > constants::pi<RealType>()) {
            z -= constants::two_pi<RealType>();
        }
        while (z < -constants::pi<RealType>()) {
            z += constants::two_pi<RealType>();
        }

        return _IMAGINARY_jacobi_theta1tau(z, 1/tau, pol);
    }
#endif

    do {
        q_n = exp(-tau * constants::pi<RealType>() * RealType(n + 0.5)*RealType(n + 0.5) );
        delta = q_n * sin(RealType(2*n+1)*z);
        if (n%2)
            delta = -delta;

        result += RealType(2) * delta;
        n++;
    } while (abs(q_n) > eps || (abs(result) > eps * eps && abs(q_n/result) > eps));

    if (abs(result) < eps * eps)
        return RealType(0);

    return result;
}

template <class T>
inline T jacobi_theta1tau(T z, T tau) {
    return jacobi_theta1tau(z, tau, policies::policy<>());
}

// First Jacobi theta function (Parameterized by q)
// = 2 * Σ (-1)^n * q^(n+1/2)^2 * sin((2n+1)z)
template <class RealType, class Policy>
inline RealType
jacobi_theta1(RealType z, RealType q, const Policy& pol) {
    BOOST_MATH_STD_USING
    if (q <= 0.0 || q >= 1.0) {
        return policies::raise_domain_error<RealType>("boost::math::jacobi_theta1<%1%>(%1%)",
                "q must be greater than 0 and less than 1 but got %1%.", q, pol);
    }
    return jacobi_theta1tau(z, -log(q)/constants::pi<RealType>(), pol);
}

template <class T>
inline T jacobi_theta1(T z, T q) {
    return jacobi_theta1(z, q, policies::policy<>());
}

// Second Jacobi theta function (Parameterized by tau - assumed imaginary)
// = 2 * Σ exp(iπτ*(n+1/2)^2) * cos((2n+1)z)
template <class RealType, class Policy>
inline RealType
jacobi_theta2tau(RealType z, RealType tau, const Policy& pol)
{
    BOOST_MATH_STD_USING
    unsigned n = 0;
    RealType eps = policies::get_epsilon<RealType, Policy>();
    RealType q_n, delta, result = RealType(0);

    if (tau <= 0.0) {
        return policies::raise_domain_error<RealType>("boost::math::jacobi_theta2tau<%1%>(%1%)",
                "tau must be greater than 0 but got %1%.", tau, pol);
    } else if (tau < 1.0 && abs(z) == 0.0) {
        return jacobi_theta4tau(z, 1/tau, pol) / sqrt(tau);
#ifdef __JACOBI_THETA_USE_IMAGINARY // DLMF 20.7.31
    } else if (exp(-tau*constants::pi<RealType>()) > 0.85) {
        z = fmod(z, constants::two_pi<RealType>());
        while (z > constants::pi<RealType>()) {
            z -= constants::two_pi<RealType>();
        }
        while (z < -constants::pi<RealType>()) {
            z += constants::two_pi<RealType>();
        }

        return _IMAGINARY_jacobi_theta4tau(z, 1/tau, pol);
#endif
    }

    do {
        q_n = exp(-tau * constants::pi<RealType>() * RealType(n + 0.5)*RealType(n + 0.5));
        delta = q_n * cos(RealType(2*n+1)*z);
        result += RealType(2) * delta;
        n++;
    //} while (abs(q_n) > eps);
    } while (abs(q_n) > eps || (abs(result) > eps * eps && abs(q_n/result) > eps));

    if (abs(result) < eps * eps)
        return RealType(0);

    return result;
}

template <class T>
inline T jacobi_theta2tau(T z, T tau) {
    return jacobi_theta2tau(z, tau, policies::policy<>());
}

// Second Jacobi theta function, parameterized by q
// = 2 * Σ q^(n+1/2)^2 * cos((2n+1)z)
template <class RealType, class Policy>
inline RealType
jacobi_theta2(RealType z, RealType q, const Policy& pol) {
    BOOST_MATH_STD_USING
    if (q <= 0.0 || q >= 1.0) {
        return policies::raise_domain_error<RealType>("boost::math::jacobi_theta2<%1%>(%1%)",
                "q must be greater than 0 and less than 1 but got %1%.", q, pol);
    }
    return jacobi_theta2tau(z, -log(q)/constants::pi<RealType>(), pol);
}

template <class T>
inline T jacobi_theta2(T z, T q) {
    return jacobi_theta2(z, q, policies::policy<>());
}

// Third Jacobi theta function, minus one (Parameterized by tau - assumed imaginary)
// This function preserves accuracy for small values of q (i.e. |q| < exp(-π) ≅ 0.043)
// For larger values of q, the minus one version usually won't help.
// = 2 * Σ exp(iπτ*(n)^2) * cos(2nz)
template <class RealType, class Policy>
inline RealType
jacobi_theta3m1tau(RealType z, RealType tau, const Policy& pol)
{
    BOOST_MATH_STD_USING

    RealType eps = policies::get_epsilon<RealType, Policy>();
    RealType q_n, delta, result = RealType(0);
    unsigned n = 1;

    do {
        q_n = exp(-tau * constants::pi<RealType>() * RealType(n)*RealType(n));
        delta = q_n * cos(RealType(2*n)*z);
        result += RealType(2) * delta;
        n++;
    // } while (abs(q_n) > eps);
    } while (abs(q_n) > eps || (abs(result) > eps * eps && abs(q_n/result) > eps));

    if (abs(result) < eps * eps)
        return RealType(0);

    return result;
}

template <class T>
inline T jacobi_theta3m1tau(T z, T tau) {
    return jacobi_theta3m1tau(z, tau, policies::policy<>());
}

// Third Jacobi theta function, parameterized by tau
// = 1 + 2 * Σ exp(iπτ*(n)^2) * cos(2nz)
template <class RealType, class Policy>
inline RealType
jacobi_theta3tau(RealType z, RealType tau, const Policy& pol)
{
    BOOST_MATH_STD_USING
    if (tau <= 0.0) {
        return policies::raise_domain_error<RealType>("boost::math::jacobi_theta3tau<%1%>(%1%)",
                "tau must be greater than 0 but got %1%.", tau, pol);
    } else if (tau < 1.0 && abs(z) == 0.0) {
        return jacobi_theta3tau(z, 1/tau, pol) / sqrt(tau);
#ifdef __JACOBI_THETA_USE_IMAGINARY // DLMF 20.7.32
    } else if (exp(-tau*constants::pi<RealType>()) > 0.85) {
        z = fmod(z, constants::pi<RealType>());
        while (z > constants::half_pi<RealType>()) {
            z -= constants::pi<RealType>();
        }
        while (z < -constants::half_pi<RealType>()) {
            z += constants::pi<RealType>();
        }
        return _IMAGINARY_jacobi_theta3tau(z, 1/tau, pol);
#endif
    }
    return RealType(1) + jacobi_theta3m1tau(z, tau, pol);
}

template <class T>
inline T jacobi_theta3tau(T z, T tau) {
    return jacobi_theta3tau(z, tau, policies::policy<>());
}

// Third Jacobi theta function, minus one (parameterized by q)
// = 2 * Σ q^n^2 * cos(2nz)
template <class RealType, class Policy>
inline RealType
jacobi_theta3m1(RealType z, RealType q, const Policy& pol) {
    BOOST_MATH_STD_USING
    if (q <= 0.0 || q >= 1.0) {
        return policies::raise_domain_error<RealType>("boost::math::jacobi_theta3m1<%1%>(%1%)",
                "q must be greater than 0 and less than 1 but got %1%.", q, pol);
    }
    return jacobi_theta3m1tau(z, -log(q)/constants::pi<RealType>(), pol);
}

template <class T>
inline T jacobi_theta3m1(T z, T q) {
    return jacobi_theta3m1(z, q, policies::policy<>());
}

// Third Jacobi theta function (parameterized by q)
// = 1 + 2 * Σ q^n^2 * cos(2nz)
template <class RealType, class Policy>
inline RealType
jacobi_theta3(RealType z, RealType q, const Policy& pol) {
    BOOST_MATH_STD_USING
    if (q <= 0.0 || q >= 1.0) {
        return policies::raise_domain_error<RealType>("boost::math::jacobi_theta3<%1%>(%1%)",
                "q must be greater than 0 and less than 1 but got %1%.", q, pol);
    }
    return jacobi_theta3tau(z, -log(q)/constants::pi<RealType>(), pol);
}

template <class T>
inline T jacobi_theta3(T z, T q) {
    return jacobi_theta3(z, q, policies::policy<>());
}

// Fourth Jacobi theta function, minus one (Parameterized by tau)
// This function preserves accuracy for small values of q (i.e. tau > 1)
// = 2 * Σ (-1)^n exp(iπτ*(n)^2) * cos(2nz)
template <class RealType, class Policy>
inline RealType
jacobi_theta4m1tau(RealType z, RealType tau, const Policy& pol)
{
    BOOST_MATH_STD_USING

    RealType eps = policies::get_epsilon<RealType, Policy>();
    RealType q_n, delta, result = RealType(0);
    unsigned n = 1;

    do {
        q_n = exp(-tau * constants::pi<RealType>() * RealType(n)*RealType(n));
        delta = q_n * cos(RealType(2*n)*z);
        if (n%2)
            delta = -delta;

        result += RealType(2) * delta;
        n++;
    } while (abs(q_n) > eps || (abs(result) > eps * eps && abs(q_n/result) > eps));

    if (abs(result) < eps * eps)
        return RealType(0);

    return result;
}

template <class T>
inline T jacobi_theta4m1tau(T z, T tau) {
    return jacobi_theta4m1tau(z, tau, policies::policy<>());
}

// Fourth Jacobi theta function (Parameterized by tau)
// = 1 + 2 * Σ (-1)^n exp(iπτ*(n)^2) * cos(2nz)
template <class RealType, class Policy>
inline RealType
jacobi_theta4tau(RealType z, RealType tau, const Policy& pol)
{
    BOOST_MATH_STD_USING
    if (tau <= 0.0) {
        return policies::raise_domain_error<RealType>("boost::math::jacobi_theta4tau<%1%>(%1%)",
                "tau must be greater than 0 but got %1%.", tau, pol);
    } else if (tau < 1.0 && abs(z) == 0.0) {
        return jacobi_theta2tau(z, 1/tau, pol) / sqrt(tau);
#ifdef __JACOBI_THETA_USE_IMAGINARY // DLMF 20.7.33
    } else if (exp(-tau*constants::pi<RealType>()) > 0.85) {
        z = fmod(z, constants::pi<RealType>());
        while (z > constants::half_pi<RealType>()) {
            z -= constants::pi<RealType>();
        }
        while (z < -constants::half_pi<RealType>()) {
            z += constants::pi<RealType>();
        }
        return _IMAGINARY_jacobi_theta2tau(z, 1/tau, pol);
#endif
    }

    return RealType(1) + jacobi_theta4m1tau(z, tau, pol);
}

template <class T>
inline T jacobi_theta4tau(T z, T tau) {
    return jacobi_theta4tau(z, tau, policies::policy<>());
}

// Fourth Jacobi theta function, minus one (Parameterized by q)
// This function preserves accuracy for small values of q
// = 2 * Σ q^n^2 * cos(2nz)
template <class RealType, class Policy>
inline RealType
jacobi_theta4m1(RealType z, RealType q, const Policy& pol) {
    BOOST_MATH_STD_USING
    if (q <= 0.0 || q >= 1.0) {
        return policies::raise_domain_error<RealType>("boost::math::jacobi_theta4m1<%1%>(%1%)",
                "q must be greater than 0 and less than 1 but got %1%.", q, pol);
    }
    return jacobi_theta4m1tau(z, -log(q)/constants::pi<RealType>(), pol);
}

template <class T>
inline T jacobi_theta4m1(T z, T q) {
    return jacobi_theta4m1(z, q, policies::policy<>());
}

// Fourth Jacobi theta function, parameterized by q
// = 1 + 2 * Σ q^n^2 * cos(2nz)
template <class RealType, class Policy>
inline RealType
jacobi_theta4(RealType z, RealType q, const Policy& pol) {
    BOOST_MATH_STD_USING
    if (abs(q) >= 1.0 || abs(q) == 0.0) {
        return policies::raise_domain_error<RealType>("boost::math::jacobi_theta4<%1%>(%1%)",
            "|q| must be greater than zero and less than 1, but got %1%.", q, pol);
    }
    return jacobi_theta4tau(z, -log(q)/constants::pi<RealType>(), pol);
}

template <class T>
inline T jacobi_theta4(T z, T q) {
    return jacobi_theta4(z, q, policies::policy<>());
}

}}

#endif
