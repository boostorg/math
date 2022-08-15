//    Copyright John Maddock 2006, 2007.
//    Copyright Paul A. Bristow 2006, 2007.
//    Copyright Philipp C. J. Muenster, 2020.
//    Copyright Matt Borland, 2022.

//    Use, modification and distribution are subject to the
//    Boost Software License, Version 1.0. (See accompanying file
//    LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_DISTRIBUTIONS_VON_MISES_HPP
#define BOOST_MATH_DISTRIBUTIONS_VON_MISES_HPP

// https://en.wikipedia.org/wiki/Von_Mises_distribution
// From MathWorld--A Wolfram Web Resource.
// http://mathworld.wolfram.com/VonMisesDistribution.html

#include <boost/math/distributions/complement.hpp>
#include <boost/math/distributions/detail/common_error_handling.hpp>
#include <boost/math/special_functions/bessel.hpp>  // for besseli0 and besseli1
#include <boost/math/special_functions/erf.hpp>     // for erf
#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/promotion.hpp>
#include <boost/math/tools/config.hpp>

#include <utility>
#include <limits>
#include <array>
#include <type_traits>

namespace boost { namespace math {

template <typename RealType = double, typename Policy = policies::policy<> >
class von_mises_distribution
{
public:
        using value_type = RealType;
        using policy_type = Policy;

    von_mises_distribution(RealType l_mean = 0, RealType concentration = 1)
        : m_mean {l_mean}, m_concentration {concentration}
    { // Default is a 'standard' von Mises distribution vM01.
        static const char* function = "boost::math::von_mises_distribution<%1%>::von_mises_distribution";

        RealType result;
        detail::check_positive_x(function, concentration, &result, Policy());
        detail::check_angle(function, l_mean, &result, Policy());
    }

    inline RealType mean() const
    { // alias for location.
        return m_mean;
    }

    inline RealType concentration() const
    { // alias for scale.
        return m_concentration;
    }

    // Synonyms, provided to allow generic use of find_location and find_scale.
    inline RealType location() const
    { // location.
        return m_mean;
    }
    inline RealType scale() const
    { // scale.
        return m_concentration;
    }

private:
    //
    // Data members:
    //
    RealType m_mean;                 // distribution mean or location.
    RealType m_concentration;        // distribution standard deviation or scale.
}; // class von_mises_distribution

using von_mises = von_mises_distribution<double>;

#ifdef __cpp_deduction_guides
template <typename RealType>
von_mises_distribution(RealType)->von_mises_distribution<boost::math::tools::promote_args_t<RealType>>;
template <typename RealType>
von_mises_distribution(RealType,RealType)->von_mises_distribution<boost::math::tools::promote_args_t<RealType>>;
#endif

#ifdef BOOST_MSVC
#pragma warning(push)
#pragma warning(disable:4127)
#endif

template <typename RealType, typename Policy>
inline std::pair<RealType, RealType> range(const von_mises_distribution<RealType, Policy>& /*dist*/)
{ // Range of permissible values for random variable x.
    BOOST_IF_CONSTEXPR (std::numeric_limits<RealType>::has_infinity)
    {
        return std::pair<RealType, RealType>(-std::numeric_limits<RealType>::infinity(),
                                             +std::numeric_limits<RealType>::infinity()); // - to + infinity.
    }
    else
    { // Can only use max_value.
        using boost::math::tools::max_value;
        return std::pair<RealType, RealType>(-max_value<RealType>(), +max_value<RealType>()); // - to + max value.
    }
}

template <typename RealType, typename Policy>
inline std::pair<RealType, RealType> support(const von_mises_distribution<RealType, Policy>& dist)
{ // This is range values for random variable x where cdf rises from 0 to 1, and outside it, the pdf is zero.
    constexpr RealType pi = boost::math::constants::pi<RealType>();
    return std::pair<RealType, RealType>(dist.mean() - pi, dist.mean() + pi); //    [µ-π, µ+π)
}

#ifdef BOOST_MSVC
#pragma warning(pop)
#endif

namespace detail {
// float version of pdf_impl
template <typename RealType, typename Policy>
inline RealType pdf_impl(const von_mises_distribution<RealType, Policy>& dist, RealType x,
                         const std::integral_constant<int, 24>&)
{
    const RealType mean = dist.mean();
    const RealType conc = dist.concentration();

    BOOST_MATH_STD_USING

    if(conc < 87)
    {
        RealType bessel_i0 = cyl_bessel_i(0, conc, Policy());
        return exp(conc * cos(x - mean)) / bessel_i0 / boost::math::constants::two_pi<RealType>();
    }
    else    // exp(88) > MAX_FLOAT
    {
        // we make use of I0(conc) = exp(conc) * P(conc), with polynomial P
        // polynomial coefficients from boost/math/special_functions/detail/bessel_i0.hpp
        static constexpr std::array<float, 5> P
        {
            3.98942280401432677e-01f,
            4.98677850501790847e-02f,
            2.80506290907257351e-02f,
            2.92194053028393074e-02f,
            4.47422143699726895e-02f
        };

        RealType result = exp(conc * (cos(x - mean) - 1.f));
        result /= boost::math::tools::evaluate_polynomial(P, RealType(1.f / conc)) / sqrt(conc)
                            * boost::math::constants::two_pi<RealType>();
        return result;
    }
}

// double version of pdf_impl
template <typename RealType, typename Policy>
inline RealType pdf_impl(const von_mises_distribution<RealType, Policy>& dist, RealType x,
                         const std::integral_constant<int, 53>&)
{
    RealType mean = dist.mean();
    RealType conc = dist.concentration();

    BOOST_MATH_STD_USING
    if(conc < 709)
    {
        RealType bessel_i0 = cyl_bessel_i(0, conc, Policy());
        return exp(conc * cos(x - mean)) / bessel_i0 / boost::math::constants::two_pi<RealType>();
    }
    else // exp(709) > MAX_DOUBLE
    {
        // we make use of I0(conc) = exp(conc) * P(conc), with polynomial P
        // polynomial coefficients from boost/math/special_functions/detail/bessel_i0.hpp
        static constexpr std::array<double, 5> P
        {
                3.98942280401432905e-01,
                4.98677850491434560e-02,
                2.80506308916506102e-02,
                2.92179096853915176e-02,
                4.53371208762579442e-02
        };

        RealType result = exp(conc * (cos(x - mean) - 1.0));
        result /= boost::math::tools::evaluate_polynomial(P, RealType(1.0 / conc)) / sqrt(conc)
                            * boost::math::constants::two_pi<RealType>();
        return result;
    }
}

// long double version of pdf_impl
template <typename RealType, typename Policy>
inline RealType pdf_impl(const von_mises_distribution<RealType, Policy>& dist, RealType x,
                         const std::integral_constant<int, 64>&)
{
    const RealType mean = dist.mean();
    const RealType conc = dist.concentration();

    BOOST_MATH_STD_USING
    if (conc < 1000)
    {
        RealType bessel_i0 = cyl_bessel_i(0, conc, Policy());
        return exp(conc * cos(x - mean)) / bessel_i0 / boost::math::constants::two_pi<RealType>();
    }
    else
    {
        // Bessel I0 over[50, INF]
        // Max error in interpolated form : 5.587e-20
        // Max Error found at float80 precision = Poly : 8.776852e-20
        static constexpr std::array<RealType, 18>  P
        {
            BOOST_MATH_BIG_CONSTANT(RealType, 64, +3.98942280401432677955074061e-01),
            BOOST_MATH_BIG_CONSTANT(RealType, 64, +4.98677850501789875615574058e-02),
            BOOST_MATH_BIG_CONSTANT(RealType, 64, +2.80506290908675604202206833e-02),
            BOOST_MATH_BIG_CONSTANT(RealType, 64, +2.92194052159035901631494784e-02),
            BOOST_MATH_BIG_CONSTANT(RealType, 64, +4.47422430732256364094681137e-02),
            BOOST_MATH_BIG_CONSTANT(RealType, 64, +9.05971614435738691235525172e-02),
            BOOST_MATH_BIG_CONSTANT(RealType, 64, +2.29180522595459823234266708e-01),
            BOOST_MATH_BIG_CONSTANT(RealType, 64, +6.15122547776140254569073131e-01),
            BOOST_MATH_BIG_CONSTANT(RealType, 64, +7.48491812136365376477357324e+00),
            BOOST_MATH_BIG_CONSTANT(RealType, 64, -2.45569740166506688169730713e+02),
            BOOST_MATH_BIG_CONSTANT(RealType, 64, +9.66857566379480730407063170e+03),
            BOOST_MATH_BIG_CONSTANT(RealType, 64, -2.71924083955641197750323901e+05),
            BOOST_MATH_BIG_CONSTANT(RealType, 64, +5.74276685704579268845870586e+06),
            BOOST_MATH_BIG_CONSTANT(RealType, 64, -8.89753803265734681907148778e+07),
            BOOST_MATH_BIG_CONSTANT(RealType, 64, +9.82590905134996782086242180e+08),
            BOOST_MATH_BIG_CONSTANT(RealType, 64, -7.30623197145529889358596301e+09),
            BOOST_MATH_BIG_CONSTANT(RealType, 64, +3.27310000726207055200805893e+10),
            BOOST_MATH_BIG_CONSTANT(RealType, 64, -6.64365417189215599168817064e+10)
        };

        RealType result = exp(conc * (cos(x - mean) - 1.0));
        result /= boost::math::tools::evaluate_polynomial(P, RealType(1.0 / conc)) / sqrt(conc)
                            * boost::math::constants::two_pi<RealType>();
        return result;
    }
}
template <typename RealType, typename Policy>
inline RealType pdf_impl(const von_mises_distribution<RealType, Policy>& dist, RealType x,
                         const std::integral_constant<int, 113>&)
{
    BOOST_MATH_STD_USING    // for ADL of std functions

    const RealType mean = dist.mean();
    const RealType conc = dist.concentration();
    const RealType bessel_i0 = cyl_bessel_i(0, conc, Policy());

    return exp(conc * cos(x - mean)) / bessel_i0 / boost::math::constants::two_pi<RealType>();
}
} // namespace detail

template <typename RealType, typename Policy>
inline RealType pdf(const von_mises_distribution<RealType, Policy>& dist, const RealType& x)
{
    BOOST_MATH_STD_USING    // for ADL of std functions

    const RealType conc = dist.concentration();
    const RealType mean = dist.mean();

    static const char* function = "boost::math::pdf(const von_mises_distribution<%1%>&, %1%)";

    RealType result = 0;
    if (!detail::check_positive_x(function, conc, &result, Policy()))
    {
        return result;
    }
    else if (!detail::check_angle(function, mean, &result, Policy()))
    {
        return result;
    }
    else if (!detail::check_angle(function, x - mean, &result, Policy()))
    {
        return result;
    }

    // Below produces MSVC 4127 warnings, so the above used instead.
    //if(std::numeric_limits<RealType>::has_infinity && abs(x) == std::numeric_limits<RealType>::infinity())
    //{ // pdf + and - infinity is zero.
    //    return 0;
    //}
    using tag_type = std::integral_constant<int,
         ((std::numeric_limits<RealType>::digits == 0) || (std::numeric_limits<RealType>::radix != 2)) ?
         0 :
         std::numeric_limits<RealType>::digits <= 24 ?
         24 :
         std::numeric_limits<RealType>::digits <= 53 ?
         53 :
         std::numeric_limits<RealType>::digits <= 64 ?
         64 :
         std::numeric_limits<RealType>::digits <= 113 ?
         113 : -1
         >;

    return detail::pdf_impl(dist, x, tag_type());
} // pdf

namespace detail {

template <typename RealType, typename Policy>
inline RealType cdf_impl(const von_mises_distribution<RealType, Policy>& dist, const RealType& x)
{
    BOOST_MATH_STD_USING    // for ADL of std functions

    RealType conc = dist.concentration();
    RealType mean = dist.mean();

    RealType u = x - mean;

    if (u <= -boost::math::constants::pi<RealType>())
    {
        return 0;
    }
    if (u >= +boost::math::constants::pi<RealType>())
    {
        return 1;
    }

    // We use the Fortran algorithm designed by Geoffrey W. Hill in
    // "Algorithm 518: Incomplete Bessel Function I0. The Von Mises Distribution", 1977, ACM
    // DOI: 10.1145/355744.355753
    RealType result = 0;

    int digits = std::numeric_limits<RealType>::max_digits10 - 1;
    RealType ck = ((0.1611*digits - 2.8778)*digits + 18.45)*digits - 35.4;
    if (conc > ck) 
    {
        RealType c = 24.0 * conc;
        RealType v = c - 56;
        RealType r = sqrt((54.0 / (347.0 / v + 26.0 - c) - 6.0 + c) / 12.0);
        RealType z = sin(u / 2.0) * r;
        RealType s = z * z * 2;
        v = v - s + 3;
        RealType y = (c - s - s - 16.0) / 3.0;
        y = ((s + 1.75) * s + 83.5) / v - y;
        result = boost::math::erf(z - s / (y * y) * z) / 2 + 0.5;
    }
    else 
    {
        RealType v = 0;
        if(conc > 0) {
            // extrapolation of the tables given in the paper
            RealType a1 = (0.33 * digits - 2.6666) * digits + 12;
            RealType a2 = (std::max)(0.5, (std::min)(1.5 - digits / 12, 1.0));
            RealType a3 = 8; //digits <= 6 ? 3 : (1 << (digits - 5));
            RealType a4 = digits <= 6 ? 1 : std::pow(1.5, digits - 8);

            auto iterations = static_cast<int>(ceil(a1 + conc * a2 - a3 / (conc + a4)));
            RealType r = 0;
            RealType z = 2 / conc;
            for (int j = iterations - 1; j > 0; --j) {
                RealType sj = sin(j * u);
                r = 1 / (j * z + r);
                v = (sj / j + v) * r;
            }
        }
        result = (x - mean + boost::math::constants::pi<RealType>()) / 2;
        result = (result + v) / boost::math::constants::pi<RealType>();
    }

    return result <= 0 ? 0 : (1 <= result ? 1 : result);
}
} // namespace detail

template <typename RealType, typename Policy>
inline RealType cdf(const von_mises_distribution<RealType, Policy>& dist, const RealType& x)
{
    RealType conc = dist.concentration();
    RealType mean = dist.mean();

    static const char* function = "boost::math::cdf(const von_mises_distribution<%1%>&, %1%)";
    
    RealType result = 0;
    if (!detail::check_positive_x(function, conc, &result, Policy()))
    {
        return result;
    }
    else if (!detail::check_angle(function, mean, &result, Policy()))
    {
        return result;
    }
    if (!detail::check_angle(function, x - mean, &result, Policy()))
    {
        return result;
    }

    return detail::cdf_impl(dist, x);
} // cdf

template <typename RealType, typename Policy>
inline RealType quantile(const von_mises_distribution<RealType, Policy>& dist, const RealType& p)
{
    BOOST_MATH_STD_USING    // for ADL of std functions

    RealType conc = dist.concentration();
    RealType mean = dist.mean();

    static const char* function
            = "boost::math::quantile(const von_mises_distribution<%1%>&, %1%)";

    RealType result = 0;
    if (!detail::check_positive_x(function, conc, &result, Policy()))
    {
        return result;
    }
    else if (!detail::check_angle(function, mean, &result, Policy()))
    {
        return result;
    }
    else if (!detail::check_probability(function, p, &result, Policy()))
    {
        return result;
    }

    if (p <= 0)
        return -boost::math::constants::pi<RealType>();
    if (p >= 1)
        return +boost::math::constants::pi<RealType>();

    using tag_type = std::integral_constant<int,
            ((std::numeric_limits<RealType>::digits == 0)
                    || (std::numeric_limits<RealType>::radix != 2)) ? 0 :
            std::numeric_limits<RealType>::digits <= 24 ? 24 :
            std::numeric_limits<RealType>::digits <= 53 ? 53 :
            std::numeric_limits<RealType>::digits <= 64 ? 64 :
            std::numeric_limits<RealType>::digits <= 113 ? 113 :
            -1
            >;

    struct step_func 
    {
        const von_mises_distribution<RealType, Policy>& dist;
        const RealType p;
        std::pair<RealType, RealType> operator()(RealType x) {
            return std::make_pair(detail::cdf_impl(dist, x) - p,                        // f(x)
                                                        detail::pdf_impl(dist, x, tag_type()));     // f'(x)
        }
    };

    RealType lower = mean - boost::math::constants::pi<RealType>();
    RealType upper = mean + boost::math::constants::pi<RealType>();
    RealType zero = boost::math::tools::newton_raphson_iterate(
            step_func{dist, p}, mean, lower, upper, 15 /* digits */);

    return zero;
} // quantile

template <typename RealType, typename Policy>
inline RealType cdf(const complemented2_type<von_mises_distribution<RealType, Policy>, RealType>& c)
{
    RealType conc = c.dist.concentration();
    RealType mean = c.dist.mean();
    RealType x = c.param;

    static const char* function
            = "boost::math::cdf(const complement(von_mises_distribution<%1%>&), %1%)";

    RealType result = 0;
    if (!detail::check_positive_x(function, conc, &result, Policy()))
    {
        return result;
    }
    if (!detail::check_angle(function, mean, &result, Policy()))
    {
        return result;
    }
    if (!detail::check_angle(function, x - mean, &result, Policy()))
    {
        return result;
    }

    return detail::cdf_impl(c.dist, 2 * mean - x);
} // cdf complement

template <typename RealType, typename Policy>
inline RealType quantile(const complemented2_type<von_mises_distribution<RealType, Policy>, RealType>& c)
{
    BOOST_MATH_STD_USING    // for ADL of std functions

    RealType conc = c.dist.concentration();
    RealType mean = c.dist.mean();

    static const char* function
            = "boost::math::quantile(const complement(von_mises_distribution<%1%>&), %1%)";

    RealType result = 0;
    if (!detail::check_positive_x(function, conc, &result, Policy()))
    {
         return result;
    }
    else if (!detail::check_angle(function, mean, &result, Policy()))
    {
         return result;
    }

    RealType q = c.param;
    if (!detail::check_probability(function, q, &result, Policy()))
    {
         return result;
    }

    if (q <= 0)
    {
        return +boost::math::constants::pi<RealType>();
    }
    else if (q >= 1)
    {
        return -boost::math::constants::pi<RealType>();
    }

    using tag_type = std::integral_constant<int,
            ((std::numeric_limits<RealType>::digits == 0)
                    || (std::numeric_limits<RealType>::radix != 2)) ? 0 :
            std::numeric_limits<RealType>::digits <= 24 ? 24 :
            std::numeric_limits<RealType>::digits <= 53 ? 53 :
            std::numeric_limits<RealType>::digits <= 64 ? 64 :
            std::numeric_limits<RealType>::digits <= 113 ? 113 :
            -1
            >;

    struct step_func 
    {
        const complemented2_type<von_mises_distribution<RealType, Policy>, RealType>& c;
        std::pair<RealType, RealType> operator()(RealType x) {
            RealType xc = 2 * c.dist.mean() - x;
            return std::make_pair(detail::cdf_impl(c.dist, xc) - c.param,        // f(x)
                                                     -detail::pdf_impl(c.dist, xc, tag_type())); // f'(x)
        }
    };

    RealType lower = mean - boost::math::constants::pi<RealType>();
    RealType upper = mean + boost::math::constants::pi<RealType>();
    RealType zero = boost::math::tools::newton_raphson_iterate(
            step_func{c}, mean, lower, upper, 15 /* digits */);

    return zero;
} // quantile

template <typename RealType, typename Policy>
inline RealType mean(const von_mises_distribution<RealType, Policy>& dist)
{
    return dist.mean();
}

template <typename RealType, typename Policy>
inline RealType standard_deviation(const von_mises_distribution<RealType, Policy>& dist)
{
    BOOST_MATH_STD_USING
    RealType bessel_quot = cyl_bessel_i(1, dist.concentration(), Policy())
                                             / cyl_bessel_i(0, dist.concentration(), Policy());
    return sqrt(-2 * log(bessel_quot));
}

template <typename RealType, typename Policy>
inline RealType mode(const von_mises_distribution<RealType, Policy>& dist)
{
    return dist.mean();
}

template <typename RealType, typename Policy>
inline RealType median(const von_mises_distribution<RealType, Policy>& dist)
{
    return dist.mean();
}

namespace detail {
// float version of variance_impl
template <typename RealType, typename Policy>
inline RealType variance_impl(const von_mises_distribution<RealType, Policy>& dist,
                              const std::integral_constant<int, 24>&)
{
    RealType conc = dist.concentration();
    BOOST_MATH_STD_USING

    if (conc < 7.75)
    {
        RealType bessel_i0 = cyl_bessel_i(0, conc, Policy());
        RealType bessel_i1 = cyl_bessel_i(1, conc, Policy());
        return 1 - bessel_i1 / bessel_i0;
    }
    else if (conc < 50)
    {
        // Polynomial coefficients from
        // boost/math/special_functions/detail/bessel_i0.hpp and
        // boost/math/special_functions/detail/bessel_i1.hpp
        // compute numerator as I0(conc) - I1(conc) from single polynomial
        static constexpr std::array<float, 5> P_numer 
        {
            +5.356107887570000000e-07f,
            +1.994139882543095464e-01f,
            +7.683426463016022940e-02f,
            +4.007722563185265850e-02f,
            +2.785578524715388070e-01f
        };
        static constexpr std::array<float, 5> P_denom
        {
            +3.98942651588301770e-01f,
            +4.98327234176892844e-02f,
            +2.91866904423115499e-02f,
            +1.35614940793742178e-02f,
            +1.31409251787866793e-01f
        };

        RealType x = 1 / conc;
        RealType numer = boost::math::tools::evaluate_polynomial(P_numer, x);
        RealType denom = boost::math::tools::evaluate_polynomial(P_denom, x);

        return numer / denom;
    }
    else
    {
        static constexpr std::array<float, 5> P_numer 
        {
            1.644239196640000000e-07f,
            1.994490498867993467e-01f,
            7.569820327857441460e-02f,
            5.573513685531774810e-02f,
            1.918908150536447035e-01f
        };
        static constexpr std::array<float, 5> P_denom
        {
            3.98942280401432677e-01f,
            4.98677850501790847e-02f,
            2.80506290907257351e-02f,
            2.92194053028393074e-02f,
            4.47422143699726895e-02f
        };

        RealType x = 1 / conc;
        RealType numer = boost::math::tools::evaluate_polynomial(P_numer, x);
        RealType denom = boost::math::tools::evaluate_polynomial(P_denom, x);

        return numer / denom;
    }
}

// double version of variance_impl
template <typename RealType, typename Policy>
inline RealType variance_impl(const von_mises_distribution<RealType, Policy>& dist,
                              const std::integral_constant<int, 53>&)
{
    RealType conc = dist.concentration();
    BOOST_MATH_STD_USING
    if (conc < 7.75)
    {
        RealType bessel_i0 = cyl_bessel_i(0, conc, Policy());
        RealType bessel_i1 = cyl_bessel_i(1, conc, Policy());
        return 1 - bessel_i1 / bessel_i0;
    }
    else if (conc < 500)
    {
        // Polynomial coefficients from
        // boost/math/special_functions/detail/bessel_i0.hpp and
        // boost/math/special_functions/detail/bessel_i1.hpp
        // compute numerator as I0(conc) - I1(conc) from single polynomial
        static constexpr std::array<double, 22> P_numer
        {
            -1.55174e-14,
            +1.994711402218073518e-01,
            +7.480166592881663552e-02,
            +7.013008203242116521e-02,
            +1.016110940936680100e-02,
            +2.837895300433059925e-01,
            -6.808807273294442296e+00
            +4.756438487430168291e+02
            -2.312449372965164219e+04,
            +8.645232325622160644e+05,
            -2.509248065298433807e+07,
            +5.705709779789331679e+08,
            -1.021122689535998576e+10,
            +1.439407686911892518e+11,
            -1.593043530624297061e+12,
            +1.373585989454675265e+13,
            -9.104389306577963414e+13,
            +4.540402039126762439e+14,
            -1.645981875199121119e+15,
            +4.091196019695783875e+15,
            -6.233158807315033853e+15,
            +4.389193640817412685e+15
        };
        static constexpr std::array<double, 22> P_denom
        {
            3.98942280401425088e-01,
            4.98677850604961985e-02,
            2.80506233928312623e-02,
            2.92211225166047873e-02,
            4.44207299493659561e-02,
            1.30970574605856719e-01,
            -3.35052280231727022e+00,
            2.33025711583514727e+02,
            -1.13366350697172355e+04,
            4.24057674317867331e+05,
            -1.23157028595698731e+07,
            2.80231938155267516e+08,
            -5.01883999713777929e+09,
            7.08029243015109113e+10,
            -7.84261082124811106e+11,
            6.76825737854096565e+12,
            -4.49034849696138065e+13,
            2.24155239966958995e+14,
            -8.13426467865659318e+14,
            2.02391097391687777e+15,
            -3.08675715295370878e+15,
            2.17587543863819074e+15
        };

        RealType x = 1 / conc;
        RealType numer = boost::math::tools::evaluate_polynomial(P_numer, x);
        RealType denom = boost::math::tools::evaluate_polynomial(P_denom, x);

        return numer / denom;
    }
    else
    {
        static constexpr std::array<double, 5> P_numer
        {
            +1.423e-15,
            +1.994711401959018717e-01,
            +7.480168411736836931e-02,
            +7.012212565916144652e-02,
            +1.037734243240472200e-01
        };
        static constexpr std::array<double, 5> P_denom 
        {
            3.98942280401432905e-01,
            4.98677850491434560e-02,
            2.80506308916506102e-02,
            2.92179096853915176e-02,
            4.53371208762579442e-02
        };
        RealType x = 1 / conc;
        RealType numer = boost::math::tools::evaluate_polynomial(P_numer, x);
        RealType denom = boost::math::tools::evaluate_polynomial(P_denom, x);

        return numer / denom;
    }
}

// long double version of variance_impl
template <typename RealType, typename Policy>
inline RealType variance_impl(const von_mises_distribution<RealType, Policy>& dist,
                              const std::integral_constant<int, 64>&)
{
    BOOST_MATH_STD_USING
    RealType conc = dist.concentration();
    if (conc < 50)
    {
        RealType bessel_i0 = cyl_bessel_i(0, conc, Policy());
        RealType bessel_i1 = cyl_bessel_i(1, conc, Policy());
        return 1 - bessel_i1 / bessel_i0;
    }
    else
    {
        // Polynomial for Bessel I0 / exp(conc) / sqrt(conc)
        static constexpr std::array<RealType, 16> P_denom 
        {
            BOOST_MATH_BIG_CONSTANT(RealType, 64, 3.9894228040143267793994605993438166526772e-01),
            BOOST_MATH_BIG_CONSTANT(RealType, 64, 4.9867785050179084742493257495245185241487e-02),
            BOOST_MATH_BIG_CONSTANT(RealType, 64, 2.8050629090725735167652437695397756897920e-02),
            BOOST_MATH_BIG_CONSTANT(RealType, 64, 2.9219405302839307466358297347675795965363e-02),
            BOOST_MATH_BIG_CONSTANT(RealType, 64, 4.4742214369972689474366968442268908028204e-02),
            BOOST_MATH_BIG_CONSTANT(RealType, 64, 9.0602984099194778006610058410222616383078e-02),
            BOOST_MATH_BIG_CONSTANT(RealType, 64, 2.2839502241666629677015839125593079416327e-01),
            BOOST_MATH_BIG_CONSTANT(RealType, 64, 6.8926354981801627920292655818232972385750e-01),
            BOOST_MATH_BIG_CONSTANT(RealType, 64, 2.4231921590621824187100989532173995000655e+00),
            BOOST_MATH_BIG_CONSTANT(RealType, 64, 9.7264260959693775207585700654645245723497e+00),
            BOOST_MATH_BIG_CONSTANT(RealType, 64, 4.3890136225398811195878046856373030127018e+01),
            BOOST_MATH_BIG_CONSTANT(RealType, 64, 2.1999720924619285464910452647408431234369e+02),
            BOOST_MATH_BIG_CONSTANT(RealType, 64, 1.2076909538525038580501368530598517194748e+03),
            BOOST_MATH_BIG_CONSTANT(RealType, 64, 7.5684635141332367730007149159063086133399e+03),
            BOOST_MATH_BIG_CONSTANT(RealType, 64, 3.5178192543258299267923025833141286569141e+04),
            BOOST_MATH_BIG_CONSTANT(RealType, 64, 6.2966297919851965784482163987240461837728e+05)
        };

        // Polynomial for (Bessel I0 - Bessel I1) / exp(conc) / sqrt(x)
        static constexpr std::array<RealType, 16> P_numer
        {
             BOOST_MATH_BIG_CONSTANT(RealType, 64, -3.368165e-34),
             BOOST_MATH_BIG_CONSTANT(RealType, 64, +0.1994711402007163389699730299733589146073),
             BOOST_MATH_BIG_CONSTANT(RealType, 64, +0.07480167757526862711373984952175456888896),
             BOOST_MATH_BIG_CONSTANT(RealType, 64, +0.07012657272681433791923412617430580227103),
             BOOST_MATH_BIG_CONSTANT(RealType, 64, +0.10226791855993757596915805134093796837369),
             BOOST_MATH_BIG_CONSTANT(RealType, 64, +0.2013399646648772644282448264813045181469),
             BOOST_MATH_BIG_CONSTANT(RealType, 64, +0.4983164125454637874163933874417571110056),
             BOOST_MATH_BIG_CONSTANT(RealType, 64, +1.4845676457582822590839158984567308297366),
             BOOST_MATH_BIG_CONSTANT(RealType, 64, +5.169476606935535670414552625141409318434),
             BOOST_MATH_BIG_CONSTANT(RealType, 64, +20.59713743665130419014001937211862931161),
             BOOST_MATH_BIG_CONSTANT(RealType, 64, +92.40031163861578044111910646492625263225),
             BOOST_MATH_BIG_CONSTANT(RealType, 64, +460.9440321063085921197536056693132597443),
             BOOST_MATH_BIG_CONSTANT(RealType, 64, +2.520575547528944544570101030955801999019e+03),
             BOOST_MATH_BIG_CONSTANT(RealType, 64, +1.5734053646329490893326466550032992462076e+04),
             BOOST_MATH_BIG_CONSTANT(RealType, 64, +7.319778356894459435808347175389511056370e+04),
             BOOST_MATH_BIG_CONSTANT(RealType, 64, +1.2997438696903014448182029282439919466883e+06)
        };

        RealType x = 1 / conc;
        RealType numer = boost::math::tools::evaluate_polynomial(P_numer, x);
        RealType denom = boost::math::tools::evaluate_polynomial(P_denom, x);

        return numer / denom;
        
        return 0;
    }
}

// quad version of variance_impl
template <typename RealType, typename Policy>
inline RealType variance_impl(const von_mises_distribution<RealType, Policy>& dist,
                              const std::integral_constant<int, 113>&)
{
    BOOST_MATH_STD_USING
    RealType conc = dist.concentration();
    if (conc < 100)
    {
        RealType bessel_i0 = cyl_bessel_i(0, conc, Policy());
        RealType bessel_i1 = cyl_bessel_i(1, conc, Policy());
        return 1 - bessel_i1 / bessel_i0;
    }
    else
    {
        // Polynomial for Bessel I0 / exp(conc) / sqrt(conc)
        static const std::array<RealType, 16> P_denom 
        {
            BOOST_MATH_BIG_CONSTANT(RealType, 113, 3.9894228040143267793994605993438166526772e-01),
            BOOST_MATH_BIG_CONSTANT(RealType, 113, 4.9867785050179084742493257495245185241487e-02),
            BOOST_MATH_BIG_CONSTANT(RealType, 113, 2.8050629090725735167652437695397756897920e-02),
            BOOST_MATH_BIG_CONSTANT(RealType, 113, 2.9219405302839307466358297347675795965363e-02),
            BOOST_MATH_BIG_CONSTANT(RealType, 113, 4.4742214369972689474366968442268908028204e-02),
            BOOST_MATH_BIG_CONSTANT(RealType, 113, 9.0602984099194778006610058410222616383078e-02),
            BOOST_MATH_BIG_CONSTANT(RealType, 113, 2.2839502241666629677015839125593079416327e-01),
            BOOST_MATH_BIG_CONSTANT(RealType, 113, 6.8926354981801627920292655818232972385750e-01),
            BOOST_MATH_BIG_CONSTANT(RealType, 113, 2.4231921590621824187100989532173995000655e+00),
            BOOST_MATH_BIG_CONSTANT(RealType, 113, 9.7264260959693775207585700654645245723497e+00),
            BOOST_MATH_BIG_CONSTANT(RealType, 113, 4.3890136225398811195878046856373030127018e+01),
            BOOST_MATH_BIG_CONSTANT(RealType, 113, 2.1999720924619285464910452647408431234369e+02),
            BOOST_MATH_BIG_CONSTANT(RealType, 113, 1.2076909538525038580501368530598517194748e+03),
            BOOST_MATH_BIG_CONSTANT(RealType, 113, 7.5684635141332367730007149159063086133399e+03),
            BOOST_MATH_BIG_CONSTANT(RealType, 113, 3.5178192543258299267923025833141286569141e+04),
            BOOST_MATH_BIG_CONSTANT(RealType, 113, 6.2966297919851965784482163987240461837728e+05)
        };

        // Polynomial for (Bessel I0 - Bessel I1) / exp(conc) / sqrt(x)
        static const std::array<RealType, 16> P_numer
        {
             BOOST_MATH_BIG_CONSTANT(RealType, 113, -3.368165e-34),
             BOOST_MATH_BIG_CONSTANT(RealType, 113, +0.1994711402007163389699730299733589146073),
             BOOST_MATH_BIG_CONSTANT(RealType, 113, +0.07480167757526862711373984952175456888896),
             BOOST_MATH_BIG_CONSTANT(RealType, 113, +0.07012657272681433791923412617430580227103),
             BOOST_MATH_BIG_CONSTANT(RealType, 113, +0.10226791855993757596915805134093796837369),
             BOOST_MATH_BIG_CONSTANT(RealType, 113, +0.2013399646648772644282448264813045181469),
             BOOST_MATH_BIG_CONSTANT(RealType, 113, +0.4983164125454637874163933874417571110056),
             BOOST_MATH_BIG_CONSTANT(RealType, 113, +1.4845676457582822590839158984567308297366),
             BOOST_MATH_BIG_CONSTANT(RealType, 113, +5.169476606935535670414552625141409318434),
             BOOST_MATH_BIG_CONSTANT(RealType, 113, +20.59713743665130419014001937211862931161),
             BOOST_MATH_BIG_CONSTANT(RealType, 113, +92.40031163861578044111910646492625263225),
             BOOST_MATH_BIG_CONSTANT(RealType, 113, +460.9440321063085921197536056693132597443),
             BOOST_MATH_BIG_CONSTANT(RealType, 113, +2.520575547528944544570101030955801999019e+03),
             BOOST_MATH_BIG_CONSTANT(RealType, 113, +1.5734053646329490893326466550032992462076e+04),
             BOOST_MATH_BIG_CONSTANT(RealType, 113, +7.319778356894459435808347175389511056370e+04),
             BOOST_MATH_BIG_CONSTANT(RealType, 113, +1.2997438696903014448182029282439919466883e+06)
        };
        RealType x = 1 / conc;
        RealType numer = boost::math::tools::evaluate_polynomial(P_numer, x);
        RealType denom = boost::math::tools::evaluate_polynomial(P_denom, x);

        return numer / denom;
    }
}
} // namespace detail

template <typename RealType, typename Policy>
inline RealType variance(const von_mises_distribution<RealType, Policy>& dist)
{

    using tag_type = std::integral_constant<int,
            ((std::numeric_limits<RealType>::digits == 0)
                    || (std::numeric_limits<RealType>::radix != 2)) ? 0 :
            std::numeric_limits<RealType>::digits <= 24 ? 24 :
            std::numeric_limits<RealType>::digits <= 53 ? 53 :
            std::numeric_limits<RealType>::digits <= 64 ? 64 :
            std::numeric_limits<RealType>::digits <= 113 ? 113 :
            -1
            >;

    return detail::variance_impl(dist, tag_type());
}

template <typename RealType, typename Policy>
inline RealType skewness(const von_mises_distribution<RealType, Policy>& /*dist*/)
{
    return 0;
}

template <typename RealType, typename Policy>
inline RealType entropy(const von_mises_distribution<RealType, Policy> & dist)
{
    BOOST_MATH_STD_USING
    RealType arg = constants::two_pi<RealType>() * cyl_bessel_i(0, dist.concentration(), Policy());
    RealType bessel_quot = cyl_bessel_i(1, dist.concentration(), Policy())
                                            / cyl_bessel_i(0, dist.concentration(), Policy());
    return log(arg) - dist.concentration() * bessel_quot;
}

} // namespace math
} // namespace boost

// This include must be at the end, *after* the accessors
// for this distribution have been defined, in order to
// keep compilers that support two-phase lookup happy.
#include <boost/math/distributions/detail/derived_accessors.hpp>

#endif // BOOST_MATH_DISTRIBUTIONS_VON_MISES_HPP