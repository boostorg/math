// Kolmogorov-Smirnov 1st order asymptotic distribution
// Copyright Evan Miller 2020
//
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0. (See accompanying file
// LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// The Kolmogorov-Smirnov test in statistics compares two empirical distributions,
// or an empirical distribution against any theoretical distribution. It makes
// use of a specific distribution which doesn't have a formal name, but which
// is often called the Kolmogorv-Smirnov distribution for lack of anything
// better. This file implements the limiting form of this distribution, first
// identified by Andrey Kolmogorov in
//
// Kolmogorov, A. (1933) "Sulla Determinazione Empirica di una Legge di
// Distribuzione." Giornale dell' Istituto Italiano degli Attuari
//
// This limiting form of the CDF is a first-order Taylor expansion that is
// easily implemented by the fourth Jacobi Theta function (setting z=0). The
// PDF is then implemented here as a derivative of the Theta function. Note
// that this derivative is with respect to x, which enters into \tau, and not
// with respect to the z argument, which is always zero, and so the derivative
// identities in DLMF 20.4 do not apply here.
//
// A higher order order expansion is possible, and was first outlined by
//
// Pelz W, Good IJ (1976). "Approximating the Lower Tail-Areas of the
// Kolmogorov-Smirnov One-sample Statistic." Journal of the Royal Statistical
// Society B.
//
// The terms in this expansion get fairly complicated, and as far as I know the
// Pelz-Good expansion is not used in any statistics software. Someone could
// consider updating this implementation to use the Pelz-Good expansion in the
// future, but the math gets considerably hairier with each additional term.
//
// A formula for an exact version of the Kolmogorov-Smirnov test is laid out in
// Equation 2.4.4 of
//
// Durbin J (1973). "Distribution Theory for Tests Based on the Sample
// Distribution Func- tion." In SIAM CBMS-NSF Regional Conference Series in
// Applied Mathematics. SIAM, Philadelphia, PA.
//
// which is available in book form from Amazon and others. This exact version
// involves taking powers of large matrices. To do that right you need to
// compute eigenvalues and eigenvectors, which are beyond the scope of Boost.
// (Some recent work indicates the exact form can also be computed via FFT, see
// https://cran.r-project.org/web/packages/KSgeneral/KSgeneral.pdf).
//
// Even if the CDF of the exact distribution could be computed using Boost
// libraries (which would be cumbersome), the PDF would present another
// difficulty. Therefore I am limiting this implementation to the asymptotic
// form, even though the exact form has trivial values for certain specific
// values of x and n. For more on trivial values see
//
// Ruben H, Gambino J (1982). "The Exact Distribution of Kolmogorov's Statistic
// Dn for n <= 10." Annals of the Institute of Statistical Mathematics.
// 
// For a good bibliography and overview of the various algorithms, including
// both exact and asymptotic forms, see
// https://www.jstatsoft.org/article/view/v039i11
//
// As for this implementation: the distribution is parameterized by n (number
// of observations) in the spirit of chi-squared's degrees of freedom. It then
// takes a single argument x. In terms of the Kolmogorov-Smirnov statistical
// test, x represents the distribution of D_n, where D_n is the maximum
// difference between the CDFs being compared, that is,
//
//   D_n = sup|F_n(x) - G(x)|
//
// In the exact distribution, x is confined to the support [0, 1], but in this
// limiting approximation, we allow x to exceed unity (similar to how a normal
// approximation always spills over any boundaries).
//
// As mentioned previously, the CDF is implemented using the \tau
// parameterization of the fourth Jacobi Theta function as
//
// CDF=theta_4(0|2*x*x*n/pi)
//
// The PDF is a hand-coded derivative of that function. Actually, there are two
// (independent) derivatives, as separate code paths are used for "small x"
// (2*x*x*n < pi) and "large x", mirroring the separate code paths in the
// Jacobi Theta implementation to achieve fast convergence. Each path evaluates
// a single exponential and obtains the remaining terms of the series by
// multiplication; the rounding error in the exponent x*x*n is tracked
// explicitly so that it is not amplified by the size of the exponent.
//
// Quantiles are computed by a Newton-Raphson iteration, but not directly in x:
// the lower tail is solved in terms of W = pi^2/(8*x*x*n) and the upper tail
// in terms of V = x*x*n, using the logarithm of the CDF (or its complement) as
// the residual. In these variables the residual is very nearly linear, so the
// iteration converges from anywhere, and the starting point (an inversion of
// the leading terms of the series) is already accurate to a fraction of a
// percent, and to working precision in the tails.
//
// The mean and variance are implemented using simple closed-form expressions.
// Skewness and kurtosis use slightly more complicated closed-form expressions
// that involve the zeta function; since these (like the mode and the median)
// do not depend on n except through a factor of sqrt(n), they are stored as
// constants of the standardized (n=1) distribution accurate to 100 decimal
// digits, and only computed at run-time for types with more precision than
// that.
//
// The CDF and PDF could almost certainly be re-implemented and sped up using a
// polynomial or rational approximation, since the only meaningful argument is
// x * sqrt(n). But that is left as an exercise for the next maintainer.
//
// In the future, the Pelz-Good approximation could be added. I suggest adding
// a second parameter representing the order, e.g.
//
// kolmogorov_smirnov_dist<>(100) // N=100, order=1
// kolmogorov_smirnov_dist<>(100, 1) // N=100, order=1, i.e. Kolmogorov's formula
// kolmogorov_smirnov_dist<>(100, 4) // N=100, order=4, i.e. Pelz-Good formula
//
// The exact distribution could be added to the API with a special order
// parameter (e.g. 0 or infinity), or a separate distribution type altogether
// (e.g. kolmogorov_smirnov_exact_distribution).
//
#ifndef BOOST_MATH_DISTRIBUTIONS_KOLMOGOROV_SMIRNOV_HPP
#define BOOST_MATH_DISTRIBUTIONS_KOLMOGOROV_SMIRNOV_HPP

#include <boost/math/distributions/fwd.hpp>
#include <boost/math/distributions/complement.hpp>
#include <boost/math/distributions/detail/common_error_handling.hpp>
#include <boost/math/special_functions/jacobi_theta.hpp>
#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/tools/big_constant.hpp>
#include <boost/math/tools/tuple.hpp>
#include <boost/math/tools/roots.hpp> // Newton-Raphson

namespace boost { namespace math {

namespace detail {

// Constants of the standardized (n=1) distribution, accurate to 100 decimal
// digits. Each accessor scales by the appropriate power of n. Types with more
// precision than the constant (as requested by the policy) fall back to a
// run-time computation.
template <class RealType, class Policy>
inline bool kolmogorov_smirnov_use_constants(const Policy&) {
    return policies::digits<RealType, Policy>() <= 330;
}

template <class RealType>
inline RealType kolmogorov_smirnov_mode_constant() {
    return BOOST_MATH_BIG_CONSTANT(RealType, 1000, 0.7354679079165719820624448513051825390913503143340225122926732679283442751250018638342479353073554823);
}

template <class RealType>
inline RealType kolmogorov_smirnov_median_constant() {
    return BOOST_MATH_BIG_CONSTANT(RealType, 1000, 0.8275735551899076901138270828889768075843727832232029452002281054669822747682290224084911530683945353);
}

template <class RealType>
inline RealType kolmogorov_smirnov_skewness_constant() {
    return BOOST_MATH_BIG_CONSTANT(RealType, 1000, 0.8604261371436682558667183685173007452393495203060929924698665373211030237490355834020123201092217303);
}

template <class RealType>
inline RealType kolmogorov_smirnov_kurtosis_excess_constant() {
    return BOOST_MATH_BIG_CONSTANT(RealType, 1000, 0.8816189679105236704015538304328295408355055629289034536289160037385510119639244078254098968815665536);
}

// Splits t = hi + lo so that hi carries only the leading half of the
// significand: products of two such hi parts are then exact (Veltkamp's
// splitting, with splitter = 2^ceil(digits/2) + 1).
template <class RealType>
inline RealType kolmogorov_smirnov_splitter() {
    BOOST_MATH_STD_USING
    return ldexp(RealType(1), (tools::digits<RealType>() + 1) / 2) + 1;
}

template <class RealType>
inline void kolmogorov_smirnov_split(RealType t, RealType splitter, RealType& hi, RealType& lo) {
    RealType g = splitter * t;
    hi = g - (g - t);
    lo = t - hi;
}

// Returns the rounding error of the product ab = a*b, i.e. the exact product
// is ab + result (Dekker's algorithm, which works for any floating-point
// type, unlike fma).
template <class RealType>
inline RealType kolmogorov_smirnov_product_error(RealType a, RealType b, RealType ab, RealType splitter) {
    RealType ah, al, bh, bl;
    kolmogorov_smirnov_split(a, splitter, ah, al);
    kolmogorov_smirnov_split(b, splitter, bh, bl);
    return (((ah * bh - ab) + ah * bl) + al * bh) + al * bl;
}

// Given u = x*x*n as computed in floating point, returns the rounding error
// in u, so that exp(-c*u) can be evaluated as exp(-c*u) * (1 - c*err)
// without the rounding error of u being amplified by the size of the
// exponent c*u.
template <class RealType>
inline RealType kolmogorov_smirnov_x2n_error(RealType x, RealType n, RealType u, RealType splitter) {
    RealType x2 = x * x;
    RealType e1 = kolmogorov_smirnov_product_error(x, x, x2, splitter);
    RealType e2 = kolmogorov_smirnov_product_error(x2, n, u, splitter);
    return e2 + e1 * n; // x*x*n == u + err
}

// The nome q = exp(-2*x*x*n) of the theta function in the CDF, with the
// rounding of x*x*n compensated. Passing q rather than tau = 2*x*x*n/pi
// matters when 2*x*x*n > pi: the theta function then raises q to integer
// powers, so nothing amplifies the rounding, whereas any rounding of tau is
// multiplied back by the size of the exponent.
template <class RealType>
inline RealType kolmogorov_smirnov_nome(RealType x, RealType n) {
    BOOST_MATH_STD_USING
    RealType u = x * x * n;
    RealType q = exp(-2 * u);
    if (q == 0)
        return q;
    RealType err = kolmogorov_smirnov_x2n_error(x, n, u, kolmogorov_smirnov_splitter<RealType>());
    return q * (1 - 2 * err); // exp(-2 (u + err))
}

// Returns pi^2/8 as hi + lo, where hi is exactly representable and lo is
// its complement to about 120 bits, so that W = (pi^2/8)/u can be computed
// without the rounding error of the constant being amplified by W. Types
// with more than 116 bits of precision simply use their own rounding of the
// constant (lo = 0).
template <class RealType>
inline RealType kolmogorov_smirnov_pi_sqr_div_eight(RealType& lo) {
    BOOST_MATH_STD_USING
    if (tools::digits<RealType>() > 116) {
        lo = 0;
        return constants::pi_sqr<RealType>() / 8;
    }
    // Five pieces of 24 significant bits each; every one is exact in any
    // binary floating-point type with at least 24 bits.
    lo = ldexp(RealType(5108270), -47) + ldexp(RealType(15913558), -71) + ldexp(RealType(14839001), -95) + ldexp(RealType(8424474), -119);
    return ldexp(RealType(10349030), -23);
}

// d/dx (theta2(0, pi/(2*x*x*n))/sqrt(2*x*x*n)) - valid for all x but converges
// quickly when 2*x*x*n < pi. With u = x*x*n and W = pi^2/(8u):
//
//  pdf = sqrt(2*pi*n)/u^2 * SUM_{i>=0} exp(-(2i+1)^2 W) * ((i+1/2)^2 pi^2 - u)
template <class RealType, class Policy>
RealType kolmogorov_smirnov_pdf_small_x(RealType x, RealType n, const Policy&) {
    BOOST_MATH_STD_USING
    RealType eps = policies::get_epsilon<RealType, Policy>();
    RealType pi2 = constants::pi_sqr<RealType>();
    RealType u = x * x * n;
    if (u == 0)
        return static_cast<RealType>(0);
    RealType c_lo;
    RealType c_hi = kolmogorov_smirnov_pi_sqr_div_eight(c_lo);
    RealType W = (c_hi + c_lo) / u;
    RealType r = exp(-W);
    RealType half = 0;
    if (r < tools::min_value<RealType>()) {
        // exp(-W) is subnormal or zero although the result, which is larger
        // by a factor of order 1/u^2, may not be: only the leading term
        // contributes, and it is evaluated below as a square of exp(-W/2)
        // to avoid the intermediate underflow.
        half = exp(-W / 2);
        if (half == 0)
            return static_cast<RealType>(0);
    }
    // The rounding errors in the constant, in u and in the division, to
    // second order, so that exp(-W) can be corrected by a factor (1 - delta):
    //   W_exact = W + delta,  delta = ((c_hi + c_lo) - W*u - W*err) / u,  u_exact = u + err
    RealType splitter = kolmogorov_smirnov_splitter<RealType>();
    RealType err = kolmogorov_smirnov_x2n_error(x, n, u, splitter);
    RealType Wu = W * u;
    RealType delta = ((((c_hi - Wu) - kolmogorov_smirnov_product_error(W, u, Wu, splitter)) + c_lo) - W * err) / u;
    if (half != 0) {
        RealType leading = half * (1 - delta / 2) * sqrt(constants::root_two_pi<RealType>() * sqrt(n) * (pi2 / 4 - u)) / u;
        return leading * leading;
    }
    r *= (1 - delta); // exp(-(W + delta))
    RealType r8 = r * r;
    r8 *= r8;
    r8 *= r8;

    // pw = r^((2i+1)^2), mult = r^(8(i+1)): pw_{i+1} = pw_i * mult_i
    RealType sum = 0, term, pw = r, mult = r8;
    int i = 0;
    do {
        term = pw * (RealType(i + 0.5) * RealType(i + 0.5) * pi2 - u);
        sum += term;
        pw *= mult;
        mult *= r8;
        i++;
    } while (term > eps * sum);

    return sum * constants::root_two_pi<RealType>() * sqrt(n) / (u * u);
}

// d/dx (theta4(0, 2*x*x*n/pi)) - valid for all x but converges quickly when
// 2*x*x*n > pi. With u = x*x*n:
//
//  pdf = 8*x*n * SUM_{i>=1} (-1)^(i+1) i^2 exp(-2 i^2 u)
template <class RealType, class Policy>
inline RealType kolmogorov_smirnov_pdf_large_x(RealType x, RealType n, const Policy&) {
    BOOST_MATH_STD_USING
    RealType eps = policies::get_epsilon<RealType, Policy>();
    RealType u = x * x * n;
    RealType r = exp(-2 * u);
    if (r < tools::min_value<RealType>()) {
        // As above: evaluate the leading term as a square so that it does
        // not lose precision to an intermediate underflow.
        RealType half = exp(-u);
        if (half == 0)
            return static_cast<RealType>(0);
        RealType err = kolmogorov_smirnov_x2n_error(x, n, u, kolmogorov_smirnov_splitter<RealType>());
        RealType leading = half * (1 - err) * sqrt(8 * x * n);
        return leading * leading;
    }
    RealType err = kolmogorov_smirnov_x2n_error(x, n, u, kolmogorov_smirnov_splitter<RealType>());
    r *= (1 - 2 * err); // exp(-2 (u + err))
    RealType r2 = r * r;

    // pw = r^(i^2), odd = r^(2i+1): pw_{i+1} = pw_i * odd_i
    RealType sum = 0, term, pw = r, odd = r2 * r;
    int i = 1;
    do {
        term = RealType(i) * RealType(i) * pw;
        if (i % 2 == 0)
            sum -= term;
        else
            sum += term;
        pw *= odd;
        odd *= r2;
        i++;
    } while (term > eps * sum);

    return 8 * x * n * sum;
}

} // detail

BOOST_MATH_EXPORT template <class RealType = double, class Policy = policies::policy<> >
    class kolmogorov_smirnov_distribution
{
    public:
        typedef RealType value_type;
        typedef Policy policy_type;

        // Constructor
    kolmogorov_smirnov_distribution( RealType n ) : n_obs_(n)
    {
        RealType result;
        detail::check_df(
                "boost::math::kolmogorov_smirnov_distribution<%1%>::kolmogorov_smirnov_distribution", n_obs_, &result, Policy());
    }

    RealType number_of_observations()const
    {
        return n_obs_;
    }

    private:

    RealType n_obs_; // positive integer
};

BOOST_MATH_EXPORT typedef kolmogorov_smirnov_distribution<double> kolmogorov_k; // Convenience typedef for double version.

#ifdef __cpp_deduction_guides
BOOST_MATH_EXPORT template <class RealType>
kolmogorov_smirnov_distribution(RealType)->kolmogorov_smirnov_distribution<typename boost::math::tools::promote_args<RealType>::type>;
#endif

namespace detail {

// The quantile is found by Newton-Raphson iteration in a transformed variable
// in which the (logarithm of the) CDF is very nearly linear:
//
// Lower tail, W = pi^2 / (8*x*x*n):
//   ln cdf(x) = ln(4/sqrt(pi)) + ln(W)/2 - W + ln(1 + exp(-8W) + exp(-24W) + ...)
//
// Upper tail, V = x*x*n:
//   ln(1 - cdf(x)) = ln(2) - 2V + ln(1 - exp(-6V) + exp(-16V) - ...)
//
// The lower form is used when the lower tail probability p is less than 0.6,
// and the upper form otherwise (with q = 1 - p). The starting values below
// invert the first two terms of each expansion; the neglected terms are
// smaller than the working precision for most of each region, in which case
// the iteration returns after a single evaluation.

// Solves W - ln(W)/2 - ln(1 + exp(-8W)) = L for the starting value of W,
// L = ln(4/sqrt(pi)) - ln(p). The iteration converges quadratically from
// the asymptotic starting point L + ln(L)/2.
template <class RealType>
inline RealType kolmogorov_smirnov_lower_guess(RealType p) {
    BOOST_MATH_STD_USING
    RealType L = 2 * constants::ln_two<RealType>() - constants::log_pi<RealType>() / 2 - log(p);
    RealType W = (L > 1) ? L + log(L) / 2 : RealType(1);
    for (int i = 0; i < 3; i++) {
        RealType e = exp(-8 * W);
        RealType h = W - log(W) / 2 - boost::math::log1p(e) - L;
        RealType dh = 1 - 1 / (2 * W) + 8 * e / (1 + e);
        W -= h / dh;
    }
    return W;
}

// Solves 2 * (t - t^4 + t^9 - t^16) = q for t = exp(-2V), a polynomial
// equation that converges in a few steps from t = q/2 when q <= 0.4, and
// returns the starting value V = -ln(t)/2. The iteration is carried out on
// the ratio r = t / (q/2) so that q may be arbitrarily small.
template <class RealType>
inline RealType kolmogorov_smirnov_upper_guess(RealType q) {
    BOOST_MATH_STD_USING
    RealType t0 = q / 2;
    RealType a = t0 * t0 * t0;
    RealType b = a * a * t0 * t0;
    RealType c = b * a * a * t0;
    RealType r = 1;
    for (int i = 0; i < 3; i++) {
        RealType r3 = r * r * r;
        RealType r8 = r3 * r3 * r * r;
        RealType r15 = r8 * r3 * r3 * r;
        RealType phi = (r - 1) - a * r3 * r + b * r8 * r - c * r15 * r;
        RealType dphi = 1 - 4 * a * r3 + 9 * b * r8 - 16 * c * r15;
        r -= phi / dphi;
    }
    return -(log(q) - constants::ln_two<RealType>() + boost::math::log1p(r - 1)) / 2;
}

template <class RealType, class Policy>
struct kolmogorov_smirnov_lower_quantile_functor
{
    kolmogorov_smirnov_lower_quantile_functor(const kolmogorov_smirnov_distribution<RealType, Policy>& dist, RealType const& p)
        : distribution(dist), scale(constants::pi_sqr<RealType>() / (8 * dist.number_of_observations()))
    {
        BOOST_MATH_STD_USING // not in scope in a member initializer
        log_prob = log(p);
    }

    boost::math::tuple<RealType, RealType> operator()(RealType const& W)
    {
        BOOST_MATH_STD_USING
        RealType x = sqrt(scale / W);
        RealType F = cdf(distribution, x);
        RealType f = pdf(distribution, x);
        // g(W) = ln F(x(W)) - ln p ; dg/dW = (f/F) * dx/dW, dx/dW = -x/(2W)
        return boost::math::make_tuple(log(F) - log_prob, -f * x / (2 * W * F));
    }
private:
    const kolmogorov_smirnov_distribution<RealType, Policy>& distribution;
    RealType log_prob;
    RealType scale;
};

template <class RealType, class Policy>
struct kolmogorov_smirnov_upper_quantile_functor
{
    kolmogorov_smirnov_upper_quantile_functor(const kolmogorov_smirnov_distribution<RealType, Policy>& dist, RealType const& q)
        : distribution(dist), n(dist.number_of_observations())
    {
        BOOST_MATH_STD_USING // not in scope in a member initializer
        log_prob = log(q);
    }

    boost::math::tuple<RealType, RealType> operator()(RealType const& V)
    {
        BOOST_MATH_STD_USING
        RealType x = sqrt(V / n);
        RealType Q = cdf(complement(distribution, x));
        RealType f = pdf(distribution, x);
        // g(V) = ln Q(x(V)) - ln q ; dg/dV = -(f/Q) * dx/dV, dx/dV = x/(2V)
        return boost::math::make_tuple(log(Q) - log_prob, -f * x / (2 * V * Q));
    }
private:
    const kolmogorov_smirnov_distribution<RealType, Policy>& distribution;
    RealType log_prob;
    RealType n;
};

// Common implementation of quantile() and quantile(complement()): p is the
// lower tail probability and q = 1 - p the upper tail probability, whichever
// was supplied by the caller having been computed from the other.
template <class RealType, class Policy>
RealType kolmogorov_smirnov_quantile_imp(const kolmogorov_smirnov_distribution<RealType, Policy>& dist, RealType p, RealType q, const char* function)
{
    BOOST_MATH_STD_USING
    RealType n = dist.number_of_observations();
    if (p == 0)
        return 0;
    if (q == 0)
        return policies::raise_overflow_error<RealType>(function, 0, Policy());

    // Newton-Raphson in W or V is quadratically convergent from a starting
    // point that is accurate to at least 1e-5, so the last step overestimates
    // the remaining error by many orders of magnitude and a couple of bits
    // can be dropped from the convergence criterion.
    const int digits = policies::digits<RealType, Policy>() - 2;
    const std::uintmax_t max_root_iterations = policies::get_max_root_iterations<Policy>();
    std::uintmax_t max_iter = max_root_iterations;
    RealType result;
    if (p < RealType(0.6)) {
        RealType W = kolmogorov_smirnov_lower_guess(p);
        // Below a thousand times the smallest normalized number, the CDF can
        // no longer be evaluated to full precision; the starting value is
        // exact to working precision there anyway.
        if (p >= 1000 * tools::min_value<RealType>()) {
            W = tools::newton_raphson_iterate(
                kolmogorov_smirnov_lower_quantile_functor<RealType, Policy>(dist, p),
                W, RealType(0), tools::max_value<RealType>(), digits, max_iter);
        } else {
            max_iter = 0;
        }
        result = sqrt(constants::pi_sqr<RealType>() / (8 * n * W));
    } else {
        RealType V = kolmogorov_smirnov_upper_guess(q);
        if (q >= 1000 * tools::min_value<RealType>()) {
            V = tools::newton_raphson_iterate(
                kolmogorov_smirnov_upper_quantile_functor<RealType, Policy>(dist, q),
                V, RealType(0), tools::max_value<RealType>(), digits, max_iter);
        } else {
            max_iter = 0;
        }
        result = sqrt(V / n);
    }
    if (max_iter >= max_root_iterations)
    {
        return policies::raise_evaluation_error<RealType>(function, "Unable to locate solution in a reasonable time:" // LCOV_EXCL_LINE
            " either there is no answer to quantile or the answer is infinite.  Current best guess is %1%", result, Policy()); // LCOV_EXCL_LINE
    }
    return result;
}

} // namespace detail

BOOST_MATH_EXPORT template <class RealType, class Policy>
inline const std::pair<RealType, RealType> range(const kolmogorov_smirnov_distribution<RealType, Policy>& /*dist*/)
{ // Range of permissible values for random variable x.
   using boost::math::tools::max_value;
   return std::pair<RealType, RealType>(static_cast<RealType>(0), max_value<RealType>());
}

BOOST_MATH_EXPORT template <class RealType, class Policy>
inline const std::pair<RealType, RealType> support(const kolmogorov_smirnov_distribution<RealType, Policy>& /*dist*/)
{ // Range of supported values for random variable x.
   // This is range where cdf rises from 0 to 1, and outside it, the pdf is zero.
   // In the exact distribution, the upper limit would be 1.
   using boost::math::tools::max_value;
   return std::pair<RealType, RealType>(static_cast<RealType>(0), max_value<RealType>());
}

BOOST_MATH_EXPORT template <class RealType, class Policy>
inline RealType pdf(const kolmogorov_smirnov_distribution<RealType, Policy>& dist, const RealType& x)
{
   BOOST_FPU_EXCEPTION_GUARD
   BOOST_MATH_STD_USING  // for ADL of std functions.

   RealType n = dist.number_of_observations();
   RealType error_result;
   static const char* function = "boost::math::pdf(const kolmogorov_smirnov_distribution<%1%>&, %1%)";
   if(false == detail::check_x_not_NaN(function, x, &error_result, Policy()))
      return error_result;

   if(false == detail::check_df(function, n, &error_result, Policy()))
      return error_result;

   if (x < 0 || !(boost::math::isfinite)(x))
   {
      return policies::raise_domain_error<RealType>(
         function, "Kolmogorov-Smirnov parameter was %1%, but must be > 0 !", x, Policy());
   }

   if (2*x*x*n < constants::pi<RealType>()) {
       return detail::kolmogorov_smirnov_pdf_small_x(x, n, Policy());
   }

   return detail::kolmogorov_smirnov_pdf_large_x(x, n, Policy());
} // pdf

BOOST_MATH_EXPORT template <class RealType, class Policy>
inline RealType cdf(const kolmogorov_smirnov_distribution<RealType, Policy>& dist, const RealType& x)
{
    BOOST_MATH_STD_USING // for ADL of std function exp.
   static const char* function = "boost::math::cdf(const kolmogorov_smirnov_distribution<%1%>&, %1%)";
   RealType error_result;
   RealType n = dist.number_of_observations();
   if(false == detail::check_x_not_NaN(function, x, &error_result, Policy()))
      return error_result;
   if(false == detail::check_df(function, n, &error_result, Policy()))
      return error_result;
   if((x < 0) || !(boost::math::isfinite)(x)) {
      return policies::raise_domain_error<RealType>(
         function, "Random variable parameter was %1%, but must be between > 0 !", x, Policy());
   }

   if (x*x*n == 0)
       return 0;

   if (2*x*x*n > constants::pi<RealType>()) {
       RealType q = detail::kolmogorov_smirnov_nome(x, n);
       if (q == 0)
           return 1;
       return RealType(1) + jacobi_theta4m1(RealType(0), q, Policy());
   }

   return jacobi_theta4tau(RealType(0), 2*x*x*n/constants::pi<RealType>(), Policy());
} // cdf

BOOST_MATH_EXPORT template <class RealType, class Policy>
inline RealType cdf(const complemented2_type<kolmogorov_smirnov_distribution<RealType, Policy>, RealType>& c) {
    BOOST_MATH_STD_USING // for ADL of std function exp.
    RealType x = c.param;
   static const char* function = "boost::math::cdf(const complemented2_type<const kolmogorov_smirnov_distribution<%1%>&, %1%>)";
   RealType error_result;
   kolmogorov_smirnov_distribution<RealType, Policy> const& dist = c.dist;
   RealType n = dist.number_of_observations();

   if(false == detail::check_x_not_NaN(function, x, &error_result, Policy()))
      return error_result;
   if(false == detail::check_df(function, n, &error_result, Policy()))
      return error_result;

   if((x < 0) || !(boost::math::isfinite)(x))
      return policies::raise_domain_error<RealType>(
         function, "Random variable parameter was %1%, but must be between > 0 !", x, Policy());

   if (x*x*n == 0)
       return 1;

   if (2*x*x*n > constants::pi<RealType>()) {
       RealType q = detail::kolmogorov_smirnov_nome(x, n);
       if (q == 0)
           return 0;
       return -jacobi_theta4m1(RealType(0), q, Policy());
   }

   return RealType(1) - jacobi_theta4tau(RealType(0), 2*x*x*n/constants::pi<RealType>(), Policy());
} // cdf (complemented)

BOOST_MATH_EXPORT template <class RealType, class Policy>
inline RealType quantile(const kolmogorov_smirnov_distribution<RealType, Policy>& dist, const RealType& p)
{
   static const char* function = "boost::math::quantile(const kolmogorov_smirnov_distribution<%1%>&, %1%)";
   // Error check:
   RealType error_result;
   RealType n = dist.number_of_observations();
   if(false == detail::check_probability(function, p, &error_result, Policy()))
      return error_result;
   if(false == detail::check_df(function, n, &error_result, Policy()))
      return error_result;

   return detail::kolmogorov_smirnov_quantile_imp(dist, p, RealType(1 - p), function);
} // quantile

BOOST_MATH_EXPORT template <class RealType, class Policy>
inline RealType quantile(const complemented2_type<kolmogorov_smirnov_distribution<RealType, Policy>, RealType>& c) {
   static const char* function = "boost::math::quantile(const complemented2_type<const kolmogorov_smirnov_distribution<%1%>&, %1%>)";
   kolmogorov_smirnov_distribution<RealType, Policy> const& dist = c.dist;
   RealType n = dist.number_of_observations();
   // Error check:
   RealType error_result;
   RealType q = c.param;
   if(false == detail::check_probability(function, q, &error_result, Policy()))
      return error_result;
   if(false == detail::check_df(function, n, &error_result, Policy()))
      return error_result;

   return detail::kolmogorov_smirnov_quantile_imp(dist, RealType(1 - q), q, function);
} // quantile (complemented)

BOOST_MATH_EXPORT template <class RealType, class Policy>
inline RealType mode(const kolmogorov_smirnov_distribution<RealType, Policy>& dist)
{
    BOOST_MATH_STD_USING
   static const char* function = "boost::math::mode(const kolmogorov_smirnov_distribution<%1%>&)";
   RealType n = dist.number_of_observations();
   RealType error_result;
   if(false == detail::check_df(function, n, &error_result, Policy()))
      return error_result;

   RealType k = detail::kolmogorov_smirnov_mode_constant<RealType>();
   if (!detail::kolmogorov_smirnov_use_constants<RealType>(Policy())) {
       // Refine by Newton-Raphson on the stationarity condition of the PDF.
       // With u = k*k and a_i = (i+1/2)^2 pi^2, the derivative of the small-x
       // PDF series with respect to u vanishes when
       //
       //   G(u) = SUM_{i>=0} exp(-a_i/(2u)) (a_i^2/(2u) - 5 a_i/2 + u) = 0
       RealType eps = policies::get_epsilon<RealType, Policy>();
       RealType pi2 = constants::pi_sqr<RealType>();
       RealType u = k * k;
       for (int iter = 0; iter < 100; iter++) {
           RealType G = 0, dG = 0, term;
           int i = 0;
           do {
               RealType a = RealType(i + 0.5) * RealType(i + 0.5) * pi2;
               RealType e = exp(-a / (2 * u));
               RealType poly = a * a / (2 * u) - 5 * a / 2 + u;
               term = e * poly;
               G += term;
               dG += e * (a / (2 * u * u) * poly - a * a / (2 * u * u) + 1);
               i++;
           } while (!(fabs(term) <= eps * fabs(G)));
           RealType step = G / dG;
           u -= step;
           if (fabs(step) <= eps * u)
               break;
       }
       k = sqrt(u);
   }
   return k / sqrt(n);
}

BOOST_MATH_EXPORT template <class RealType, class Policy>
inline RealType median(const kolmogorov_smirnov_distribution<RealType, Policy>& dist)
{
    BOOST_MATH_STD_USING
   static const char* function = "boost::math::median(const kolmogorov_smirnov_distribution<%1%>&)";
   RealType n = dist.number_of_observations();
   RealType error_result;
   if(false == detail::check_df(function, n, &error_result, Policy()))
      return error_result;
   if (!detail::kolmogorov_smirnov_use_constants<RealType>(Policy()))
       return quantile(dist, RealType(0.5));
   return detail::kolmogorov_smirnov_median_constant<RealType>() / sqrt(n);
}

// Mean and variance come directly from
// https://www.jstatsoft.org/article/view/v008i18 Section 3
BOOST_MATH_EXPORT template <class RealType, class Policy>
inline RealType mean(const kolmogorov_smirnov_distribution<RealType, Policy>& dist)
{
    BOOST_MATH_STD_USING
   static const char* function = "boost::math::mean(const kolmogorov_smirnov_distribution<%1%>&)";
    RealType n = dist.number_of_observations();
    RealType error_result;
    if(false == detail::check_df(function, n, &error_result, Policy()))
        return error_result;
    return constants::root_half_pi<RealType>() * constants::ln_two<RealType>() / sqrt(n);
}

BOOST_MATH_EXPORT template <class RealType, class Policy>
inline RealType variance(const kolmogorov_smirnov_distribution<RealType, Policy>& dist)
{
   static const char* function = "boost::math::variance(const kolmogorov_smirnov_distribution<%1%>&)";
    RealType n = dist.number_of_observations();
    RealType error_result;
    if(false == detail::check_df(function, n, &error_result, Policy()))
        return error_result;
    return (constants::pi_sqr_div_six<RealType>()
            - constants::pi<RealType>() * constants::ln_two<RealType>() * constants::ln_two<RealType>()) / (2*n);
}

// Skewness and kurtosis come from integrating the PDF
// The alternating series pops out a Dirichlet eta function which is related to the zeta function
BOOST_MATH_EXPORT template <class RealType, class Policy>
inline RealType skewness(const kolmogorov_smirnov_distribution<RealType, Policy>& dist)
{
    BOOST_MATH_STD_USING
   static const char* function = "boost::math::skewness(const kolmogorov_smirnov_distribution<%1%>&)";
    RealType n = dist.number_of_observations();
    RealType error_result;
    if(false == detail::check_df(function, n, &error_result, Policy()))
        return error_result;
    if (detail::kolmogorov_smirnov_use_constants<RealType>(Policy()))
        return detail::kolmogorov_smirnov_skewness_constant<RealType>();
    RealType ex3 = RealType(0.5625) * constants::root_half_pi<RealType>() * constants::zeta_three<RealType>() / n / sqrt(n);
    RealType mean = boost::math::mean(dist);
    RealType var = boost::math::variance(dist);
    return (ex3 - 3 * mean * var - mean * mean * mean) / var / sqrt(var);
}

BOOST_MATH_EXPORT template <class RealType, class Policy>
inline RealType kurtosis_excess(const kolmogorov_smirnov_distribution<RealType, Policy>& dist)
{
    BOOST_MATH_STD_USING
   static const char* function = "boost::math::kurtosis_excess(const kolmogorov_smirnov_distribution<%1%>&)";
    RealType n = dist.number_of_observations();
    RealType error_result;
    if(false == detail::check_df(function, n, &error_result, Policy()))
        return error_result;
    if (detail::kolmogorov_smirnov_use_constants<RealType>(Policy()))
        return detail::kolmogorov_smirnov_kurtosis_excess_constant<RealType>();
    RealType ex4 = 7 * constants::pi_sqr_div_six<RealType>() * constants::pi_sqr_div_six<RealType>() / 20 / n / n;
    RealType mean = boost::math::mean(dist);
    RealType var = boost::math::variance(dist);
    RealType skew = boost::math::skewness(dist);
    return (ex4 - 4 * mean * skew * var * sqrt(var) - 6 * mean * mean * var - mean * mean * mean * mean) / var / var - 3;
}

BOOST_MATH_EXPORT template <class RealType, class Policy>
inline RealType kurtosis(const kolmogorov_smirnov_distribution<RealType, Policy>& dist)
{
    return kurtosis_excess(dist) + 3;
}
}}
#endif
