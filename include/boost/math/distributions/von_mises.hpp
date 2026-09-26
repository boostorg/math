//    Copyright John Maddock 2006, 2007.
//    Copyright Paul A. Bristow 2006, 2007.
//    Copyright Philipp C. J. Muenster, 2020.
//    Copyright Matt Borland, 2022.
//
//    Use, modification and distribution are subject to the
//    Boost Software License, Version 1.0. (See accompanying file
//    LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_DISTRIBUTIONS_VON_MISES_HPP
#define BOOST_MATH_DISTRIBUTIONS_VON_MISES_HPP

// https://en.wikipedia.org/wiki/Von_Mises_distribution
// From MathWorld--A Wolfram Web Resource.
// http://mathworld.wolfram.com/VonMisesDistribution.html

#include <boost/math/distributions/fwd.hpp>
#include <boost/math/distributions/complement.hpp>
#include <boost/math/distributions/detail/common_error_handling.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/promotion.hpp>
#include <boost/math/tools/config.hpp>
#include <boost/math/constants/constants.hpp>

#ifndef BOOST_MATH_BUILD_MODULE
#include <cstdint>
#include <utility>
#include <limits>
#endif

namespace boost { namespace math {

namespace detail {

template <class RealType, class Policy>
inline bool check_von_mises_concentration(const char* function, RealType concentration, RealType* result, const Policy& pol)
{
   if (!(boost::math::isfinite)(concentration) || (concentration < 0))
   {
      *result = policies::raise_domain_error<RealType>(
         function, "Concentration parameter is %1%, but must be finite and >= 0!", concentration, pol);
      return false;
   }
   return true;
}

template <class RealType, class Policy>
inline bool check_von_mises(const char* function, RealType mean, RealType concentration, RealType* result, const Policy& pol)
{
   return check_location(function, mean, result, pol)
      && check_von_mises_concentration(function, concentration, result, pol);
}

} // namespace detail

BOOST_MATH_EXPORT template <typename RealType = double, typename Policy = policies::policy<> >
class von_mises_distribution
{
public:
    using value_type = RealType;
    using policy_type = Policy;

    explicit von_mises_distribution(RealType l_mean = 0, RealType concentration = 1)
        : m_mean {l_mean}, m_concentration {concentration}
    { // Default is a 'standard' von Mises distribution vM01.
        RealType result;
        detail::check_von_mises("boost::math::von_mises_distribution<%1%>::von_mises_distribution",
                                l_mean, concentration, &result, Policy());
    }

    RealType mean() const
    { // alias for location.
        return m_mean;
    }

    RealType concentration() const
    { // alias for scale.
        return m_concentration;
    }

    // Synonyms, provided to allow generic use of find_location and find_scale.
    RealType location() const
    {
        return m_mean;
    }
    RealType scale() const
    {
        return m_concentration;
    }

private:
    RealType m_mean;
    RealType m_concentration;
}; // class von_mises_distribution

BOOST_MATH_EXPORT using von_mises = von_mises_distribution<double>;

#ifdef __cpp_deduction_guides
BOOST_MATH_EXPORT template <typename RealType>
von_mises_distribution(RealType)->von_mises_distribution<boost::math::tools::promote_args_t<RealType>>;
BOOST_MATH_EXPORT template <typename RealType>
von_mises_distribution(RealType,RealType)->von_mises_distribution<boost::math::tools::promote_args_t<RealType>>;
#endif

BOOST_MATH_EXPORT template <typename RealType, typename Policy>
inline std::pair<RealType, RealType> range(const von_mises_distribution<RealType, Policy>& /*dist*/)
{ // Range of permissible values for random variable x.
    using boost::math::tools::max_value;
    return std::pair<RealType, RealType>(-max_value<RealType>(), max_value<RealType>());
}

BOOST_MATH_EXPORT template <typename RealType, typename Policy>
inline std::pair<RealType, RealType> support(const von_mises_distribution<RealType, Policy>& dist)
{ // Range of x where the pdf is non-zero: one full turn centred on the mean.
    const RealType pi = boost::math::constants::pi<RealType>();
    return std::pair<RealType, RealType>(dist.mean() - pi, dist.mean() + pi);
}

namespace detail {

// e^-k I0(k), which stays finite for all k.
template <typename RealType, typename Policy>
RealType von_mises_scaled_i0(RealType k, const Policy& pol)
{
    BOOST_MATH_STD_USING
    if (k < tools::log_max_value<RealType>())
    {
        return cyl_bessel_i(0, k, pol) * exp(-k);
    }
    // Here k is far beyond the precision, so the asymptotic series converges before its terms grow.
    RealType term = 1;
    RealType sum = 1;
    for (unsigned j = 1; term > tools::epsilon<RealType>() * sum; ++j)
    {
        term *= RealType((2 * j - 1) * (2 * j - 1)) / (8 * j * k);
        sum += term;
    }
    return sum / sqrt(constants::two_pi<RealType>() * k);
}

// The circular variance 1 - I1(k)/I0(k) and ln(2 pi e^-k I0(k)), both without underflow and
// without the cancellation of forming I0 - I1: every sum below has positive terms.
template <typename RealType, typename Policy>
void von_mises_moments(RealType k, RealType* circular_variance, RealType* log_two_pi_scaled_i0, const Policy&)
{
    BOOST_MATH_STD_USING
    const RealType eps = tools::epsilon<RealType>();
    // Beyond this the smallest term of the asymptotic series, about e^(-2k), is below eps.
    const RealType asymptotic_limit = tools::digits<RealType>() * constants::ln_two<RealType>();
    if (k >= asymptotic_limit)
    {
        // e^-k I_nu(k) ~ (2 pi k)^(-1/2) sum_j (-1)^j a_j(nu) k^-j; every term of the I0 - I1 series is positive.
        RealType t0 = 1;
        RealType t1 = 1;
        RealType sum0 = 1;
        RealType sumd = 0;
        for (unsigned j = 1; ; ++j)
        {
            t0 *= RealType((2 * j - 1) * (2 * j - 1)) / (8 * j * k);
            t1 *= (RealType(2 * j) - 3) * (2 * j + 1) / (8 * j * k);
            sum0 += t0;
            sumd += t0 - t1;
            if (t0 - t1 <= eps * sumd)
            {
                break;
            }
        }
        *circular_variance = sumd / sum0;
        *log_two_pi_scaled_i0 = log(sum0 * sqrt(constants::two_pi<RealType>() / k));
        return;
    }
    // (1/pi) int_0^pi e^(-k(1 - cos t)) {1, 1 - cos t} dt by the trapezoidal rule, which for these
    // periodic integrands at least squares the error with each halving of the step.
    const RealType tol = sqrt(eps);
    RealType f_end = exp(-2 * k);
    RealType sum0 = (1 + f_end) / 2;
    RealType sumd = f_end;
    unsigned n = 1;
    for (unsigned level = 0; level < 30; ++level)
    {
        RealType mid0 = 0;
        RealType midd = 0;
        for (unsigned j = 1; j < 2 * n; j += 2)
        {
            RealType s = sin(constants::pi<RealType>() * j / (4 * n));
            RealType g = 2 * s * s;
            RealType f = exp(-k * g);
            mid0 += f;
            midd += f * g;
        }
        // Compare the midpoint rule with the trapezoidal rule of the previous level.
        const bool converged = (level > 1) && (fabs(mid0 - sum0) <= tol * sum0) && (fabs(midd - sumd) <= tol * sumd);
        sum0 += mid0;
        sumd += midd;
        n *= 2;
        if (converged)
        {
            break;
        }
    }
    *circular_variance = sumd / sum0;
    *log_two_pi_scaled_i0 = log(constants::two_pi<RealType>() * sum0 / n);
}

// int_0^width f(s) ds for the positive, smooth integrands below, which fall by at most the cutoff's
// number of e-folds over the interval. That interval is not a full period, so the trapezoidal rule
// would converge only as h^2: use Gauss-Legendre up to double precision, and tanh-sinh beyond.
template <typename RealType, typename Policy, typename F>
RealType von_mises_integrate(F f, RealType width, const Policy&)
{
    BOOST_MATH_IF_CONSTEXPR (std::numeric_limits<RealType>::is_specialized && (std::numeric_limits<RealType>::digits <= 53))
    {
        return quadrature::gauss<RealType, 30>::integrate(f, RealType(0), width);
    }
    else
    {
        static const quadrature::tanh_sinh<RealType, Policy> integrator;
        return integrator.integrate(f, RealType(0), width, policies::get_epsilon<RealType, Policy>());
    }
}

// Beyond this many e-folds the integrands are negligible: the part cut off is at most about
// e^-L sqrt(2L) of what is kept.
template <typename RealType, typename Policy>
RealType von_mises_cutoff()
{
    BOOST_MATH_STD_USING
    const RealType l0 = 1 - log(policies::get_epsilon<RealType, Policy>());
    return l0 + log(l0);
}

// The cdf of the standard distribution at -a, 0 <= a <= pi, which is at most 1/2:
//   F(-a) = (1/(2 pi e^-k I0(k))) int_a^pi e^(-k(1 - cos t)) dt.
// The integrand is positive, so the left tail keeps its relative accuracy.
// It falls off like e^(-k(cos a - cos t)), so the range is cut where it drops below eps.
template <typename RealType, typename Policy>
RealType von_mises_cdf_left(RealType k, RealType a, RealType s0, const Policy& pol)
{
    BOOST_MATH_STD_USING
    const RealType pi = constants::pi<RealType>();
    if (a >= pi)
    {
        return 0;
    }
    if (k == 0)
    {
        return (pi - a) / constants::two_pi<RealType>();
    }
    const RealType sin_half_a = sin(a / 2);
    const RealType s2 = sin_half_a * sin_half_a + von_mises_cutoff<RealType, Policy>() / (2 * k);
    const RealType width = (s2 >= 1 ? pi : RealType(2 * asin(sqrt(s2)))) - a;
    // In terms of the offset s = t - a, so that rounding t does not perturb the steep exponent:
    // cos a - cos t = 2 sin(a + s/2) sin(s/2).
    auto integrand = [&](RealType s) -> RealType
    {
        return exp(-2 * k * sin(a + s / 2) * sin(s / 2));
    };
    return exp(-2 * k * sin_half_a * sin_half_a) * von_mises_integrate(integrand, width, pol) / (constants::two_pi<RealType>() * s0);
}

// The probability between the mean and a, 0 <= a <= pi, with relative accuracy even as a -> 0:
//   F(a) - 1/2 = (1/(2 pi e^-k I0(k))) int_0^a e^(-k(1 - cos t)) dt.
template <typename RealType, typename Policy>
RealType von_mises_cdf_centre(RealType k, RealType a, RealType s0, const Policy& pol)
{
    BOOST_MATH_STD_USING
    RealType width = a;
    if (k > 0)
    {
        const RealType s2 = von_mises_cutoff<RealType, Policy>() / (2 * k);
        if (s2 < 1)
        {
            width = (std::min)(a, RealType(2 * asin(sqrt(s2))));
        }
    }
    auto integrand = [&](RealType t) -> RealType
    {
        RealType s = sin(t / 2);
        return exp(-2 * k * s * s);
    };
    return von_mises_integrate(integrand, width, pol) / (constants::two_pi<RealType>() * s0);
}

template <typename RealType, typename Policy>
RealType von_mises_pdf_imp(RealType k, RealType u, RealType s0)
{
    BOOST_MATH_STD_USING
    RealType s = sin(u / 2);
    return exp(-2 * k * s * s) / (constants::two_pi<RealType>() * s0);
}

// Solves F(u) = p for u in [-pi, 0], where p <= 1/2.
template <typename RealType, typename Policy>
RealType von_mises_quantile_left(RealType k, RealType p, const Policy& pol)
{
    BOOST_MATH_STD_USING
    static const char* function = "boost::math::quantile(const von_mises_distribution<%1%>&, %1%)";
    const RealType pi = constants::pi<RealType>();
    if (p <= 0)
    {
        return -pi;
    }
    const RealType s0 = von_mises_scaled_i0(k, pol);
    std::uintmax_t max_iter = policies::get_max_root_iterations<Policy>();
    const int digits = policies::digits<RealType, Policy>();
    RealType result;
    if (p >= 0.25f)
    {
        // Near the median solve for the distance a below it, so that the root is not
        // lost against the 1/2 it sits beside: q = 1/2 - p is exact here.
        const RealType q = 0.5f - p;
        if (q == 0)
        {
            return 0;
        }
        RealType guess;
        if (k < 1)
        {
            guess = constants::two_pi<RealType>() * q;
        }
        else
        {
            // Wrapped normal with variance 1/k: 2 sqrt(k) sin(u/2) is close to standard normal.
            RealType z = constants::root_two<RealType>() * erf_inv(2 * q, pol) / (2 * sqrt(k));
            guess = z >= 1 ? RealType(pi / 2) : RealType(2 * asin(z));
        }
        auto f = [&](RealType a) -> std::pair<RealType, RealType>
        {
            return std::make_pair(von_mises_cdf_centre(k, a, s0, pol) - q, von_mises_pdf_imp<RealType, Policy>(k, a, s0));
        };
        result = -tools::newton_raphson_iterate(f, guess, RealType(0), pi, digits, max_iter);
    }
    else
    {
        RealType guess;
        if (k < 1)
        {
            guess = (2 * p - 1) * pi;
        }
        else
        {
            RealType z = -constants::root_two<RealType>() * erfc_inv(2 * p, pol) / (2 * sqrt(k));
            guess = z <= -1 ? RealType(-pi / 2) : RealType(2 * asin(z));
        }
        auto f = [&](RealType u) -> std::pair<RealType, RealType>
        {
            return std::make_pair(von_mises_cdf_left(k, -u, s0, pol) - p, von_mises_pdf_imp<RealType, Policy>(k, u, s0));
        };
        result = tools::newton_raphson_iterate(f, guess, -pi, RealType(0), digits, max_iter);
    }
    if (max_iter >= policies::get_max_root_iterations<Policy>())
    {
        return policies::raise_evaluation_error<RealType>(function,
            "Unable to locate solution in a reasonable time: either there is no answer to quantile or the answer is infinite.  Current best guess is %1%", result, pol);
    }
    return result;
}

} // namespace detail

BOOST_MATH_EXPORT template <typename RealType, typename Policy>
inline RealType pdf(const von_mises_distribution<RealType, Policy>& dist, const RealType& x)
{
    BOOST_MATH_STD_USING
    static const char* function = "boost::math::pdf(const von_mises_distribution<%1%>&, %1%)";
    const RealType k = dist.concentration();
    const RealType mean = dist.mean();
    RealType result = 0;
    if (!detail::check_von_mises(function, mean, k, &result, Policy())
        || !detail::check_x(function, x, &result, Policy()))
    {
        return result;
    }
    const RealType u = x - mean;
    if (fabs(u) > constants::pi<RealType>())
    {
        return 0;
    }
    return detail::von_mises_pdf_imp<RealType, Policy>(k, u, detail::von_mises_scaled_i0(k, Policy()));
}

BOOST_MATH_EXPORT template <typename RealType, typename Policy>
inline RealType cdf(const von_mises_distribution<RealType, Policy>& dist, const RealType& x)
{
    static const char* function = "boost::math::cdf(const von_mises_distribution<%1%>&, %1%)";
    const RealType k = dist.concentration();
    const RealType mean = dist.mean();
    RealType result = 0;
    if (!detail::check_von_mises(function, mean, k, &result, Policy())
        || !detail::check_x(function, x, &result, Policy()))
    {
        return result;
    }
    const RealType u = x - mean;
    const RealType s0 = detail::von_mises_scaled_i0(k, Policy());
    return u <= 0 ? detail::von_mises_cdf_left(k, -u, s0, Policy())
                  : 1 - detail::von_mises_cdf_left(k, u, s0, Policy());
}

BOOST_MATH_EXPORT template <typename RealType, typename Policy>
inline RealType cdf(const complemented2_type<von_mises_distribution<RealType, Policy>, RealType>& c)
{
    static const char* function = "boost::math::cdf(const complement(von_mises_distribution<%1%>&), %1%)";
    const RealType k = c.dist.concentration();
    const RealType mean = c.dist.mean();
    const RealType x = c.param;
    RealType result = 0;
    if (!detail::check_von_mises(function, mean, k, &result, Policy())
        || !detail::check_x(function, x, &result, Policy()))
    {
        return result;
    }
    const RealType u = x - mean;
    const RealType s0 = detail::von_mises_scaled_i0(k, Policy());
    return u >= 0 ? detail::von_mises_cdf_left(k, u, s0, Policy())
                  : 1 - detail::von_mises_cdf_left(k, -u, s0, Policy());
}

BOOST_MATH_EXPORT template <typename RealType, typename Policy>
inline RealType quantile(const von_mises_distribution<RealType, Policy>& dist, const RealType& p)
{
    static const char* function = "boost::math::quantile(const von_mises_distribution<%1%>&, %1%)";
    const RealType k = dist.concentration();
    const RealType mean = dist.mean();
    RealType result = 0;
    if (!detail::check_von_mises(function, mean, k, &result, Policy())
        || !detail::check_probability(function, p, &result, Policy()))
    {
        return result;
    }
    // By symmetry F(-u) = 1 - F(u).
    return p <= 0.5f ? mean + detail::von_mises_quantile_left(k, p, Policy())
                     : mean - detail::von_mises_quantile_left(k, RealType(1 - p), Policy());
}

BOOST_MATH_EXPORT template <typename RealType, typename Policy>
inline RealType quantile(const complemented2_type<von_mises_distribution<RealType, Policy>, RealType>& c)
{
    static const char* function = "boost::math::quantile(const complement(von_mises_distribution<%1%>&), %1%)";
    const RealType k = c.dist.concentration();
    const RealType mean = c.dist.mean();
    const RealType q = c.param;
    RealType result = 0;
    if (!detail::check_von_mises(function, mean, k, &result, Policy())
        || !detail::check_probability(function, q, &result, Policy()))
    {
        return result;
    }
    return q <= 0.5f ? mean - detail::von_mises_quantile_left(k, q, Policy())
                     : mean + detail::von_mises_quantile_left(k, RealType(1 - q), Policy());
}

BOOST_MATH_EXPORT template <typename RealType, typename Policy>
inline RealType mean(const von_mises_distribution<RealType, Policy>& dist)
{
    return dist.mean();
}

BOOST_MATH_EXPORT template <typename RealType, typename Policy>
inline RealType mode(const von_mises_distribution<RealType, Policy>& dist)
{
    return dist.mean();
}

BOOST_MATH_EXPORT template <typename RealType, typename Policy>
inline RealType median(const von_mises_distribution<RealType, Policy>& dist)
{
    return dist.mean();
}

// Circular variance, 1 - I1(k)/I0(k).
BOOST_MATH_EXPORT template <typename RealType, typename Policy>
inline RealType variance(const von_mises_distribution<RealType, Policy>& dist)
{
    RealType v, log_two_pi_s0;
    detail::von_mises_moments(dist.concentration(), &v, &log_two_pi_s0, Policy());
    return v;
}

// Circular standard deviation, sqrt(-2 ln(I1(k)/I0(k))).
BOOST_MATH_EXPORT template <typename RealType, typename Policy>
inline RealType standard_deviation(const von_mises_distribution<RealType, Policy>& dist)
{
    BOOST_MATH_STD_USING
    if (dist.concentration() == 0)
    {
        return policies::raise_overflow_error<RealType>("boost::math::standard_deviation(const von_mises_distribution<%1%>&)",
            "The circular standard deviation is infinite when the concentration is zero.", Policy());
    }
    RealType v, log_two_pi_s0;
    detail::von_mises_moments(dist.concentration(), &v, &log_two_pi_s0, Policy());
    return sqrt(-2 * log1p(-v, Policy()));
}

BOOST_MATH_EXPORT template <typename RealType, typename Policy>
inline RealType skewness(const von_mises_distribution<RealType, Policy>& /*dist*/)
{
    return 0;
}

BOOST_MATH_EXPORT template <typename RealType, typename Policy>
inline RealType entropy(const von_mises_distribution<RealType, Policy>& dist)
{
    BOOST_MATH_STD_USING
    // ln(2 pi I0(k)) - k I1(k)/I0(k), with the e^k factors cancelled.
    const RealType k = dist.concentration();
    RealType v, log_two_pi_s0;
    detail::von_mises_moments(k, &v, &log_two_pi_s0, Policy());
    return log_two_pi_s0 + k * v;
}

} // namespace math
} // namespace boost

// This include must be at the end, *after* the accessors
// for this distribution have been defined, in order to
// keep compilers that support two-phase lookup happy.
#include <boost/math/distributions/detail/derived_accessors.hpp>

#endif // BOOST_MATH_DISTRIBUTIONS_VON_MISES_HPP
