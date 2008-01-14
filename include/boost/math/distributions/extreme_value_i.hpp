// Copyright John Maddock 2006
//
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_STATS_EXTREME_VALUE_I_HPP
#define BOOST_STATS_EXTREME_VALUE_I_HPP

#include <boost/math/distributions/fwd.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/special_functions/expm1.hpp>

#include <boost/math/distributions/complement.hpp>
#include <boost/math/distributions/extreme_value.hpp>
#include <boost/math/distributions/detail/common_error_handling.hpp>
#include <cmath>

#include <utility>

#ifdef BOOST_MSVC
# pragma warning(push)
# pragma warning(disable: 4702) // unreachable code (return after domain_error throw).
#endif

namespace boost {
namespace math {
namespace detail{

    template <class RealType, class Policy>
    inline bool verify_sign(const char* function, RealType c, RealType* presult, const Policy& pol)
    {
        if(c != -1 && c != 1) {
            *presult = policies::raise_domain_error<RealType>(
                function,
                "The sign parameter \"c\" must be -1 or 1, but was: %1%.", c, pol);
            return false;
        }
        return true;
    }
} // namespace detail

    // This is the maximum extreme value distribution, see
    // http://www.itl.nist.gov/div898/handbook/eda/section3/eda366g.htm
    // and http://mathworld.wolfram.com/ExtremeValueDistribution.html
    // Also known as a Fisher-Tippett distribution, a log-Weibull
    // distribution or a Gumbel distribution.
    //
    // This version of the extreme_value_i_distribution is modified to take a
    // third parameter, the sign (Note not a shape parameter) that causes the
    // distribution to represent either the maximum case (sign == 1) or the
    // minimum case (sign == -1). This defaults to the maximum case.
    //
    // This distribution actually implements a class of several different
    // distributions. The extreme value distribution deals with the median
    // of different distributions. There are three distinct flavors of extreme
    // value distributions (I, II, and III). For type I distributions, there
    // are two cases: minimum and maximum. Generally when authors refer to the
    // type I distribution, they refer to the maximum case. The different cases
    // only differ by a factor of -1 in most cases.
    //
    // Note that the type II and III are closely related to the minimum and
    // maximum cases of the type I.
    template <class RealType = double, class Policy = policies::policy<> >
    class extreme_value_i_distribution
    {
    public:
        typedef RealType value_type;
        typedef Policy policy_type;

        extreme_value_i_distribution(RealType a = RealType(0),
                                RealType b = RealType(1),
                                RealType s = RealType(1))
            : m_a(a)
            , m_b(b)
            , m_sign(s)
        {
            static const char* func = "boost::math::extreme_value_i_distribution<%1%>::extreme_value_i_distribution";

            RealType err;
            detail::verify_scale_b(func, b, &err, Policy());
            detail::verify_sign(func, s, &err, Policy());
        }

        RealType location() const
        { return m_a; }

        RealType scale() const
        { return m_b; }

        RealType sign() const
        { return m_sign; }

    private:
        RealType m_a;
        RealType m_b;
        RealType m_sign;
    };

    typedef extreme_value_i_distribution<double> extreme_value_i;

    // Range of permissible values for random variable x.
    template <class RealType, class Policy>
    inline const std::pair<RealType, RealType>
    range(const extreme_value_i_distribution<RealType, Policy>&)
    {
        using boost::math::tools::max_value;
        return std::make_pair(-max_value<RealType>(), max_value<RealType>());
    }

    // Range of supported values for random variable x.
    // This is range where cdf rises from 0 to 1, and outside it, the pdf is zero.
    template <class RealType, class Policy>
    inline const std::pair<RealType, RealType>
    support(const extreme_value_i_distribution<RealType, Policy>&)
    {
        using boost::math::tools::max_value;
        return std::make_pair(-max_value<RealType>(),  max_value<RealType>());
    }

    template <class RealType, class Policy>
    inline RealType
    pdf(const extreme_value_i_distribution<RealType, Policy>& dist, const RealType& x)
    {
        static const char* func = "boost::math::pdf(const extreme_value_i_distribution<%1%>&, %1%)";

        BOOST_MATH_STD_USING; // for ADL of std functions

        RealType a = dist.location();
        RealType b = dist.scale();
        RealType s = dist.sign();
        RealType c = s * (a - x) / b;
        RealType result;

        if(!detail::verify_scale_b(func, b, &result, Policy())) {
            return result;
        }
        if(!detail::verify_sign(func, s, &result, Policy())) {
            return result;
        }

        result = exp(c) * exp(-exp(c)) / b;
        return result;
    }

    template <class RealType, class Policy>
    inline RealType
    cdf(const extreme_value_i_distribution<RealType, Policy>& dist, const RealType& x)
    {
        static const char* func = "boost::math::cdf(const extreme_value_i_distribution<%1%>&, %1%)";

        BOOST_MATH_STD_USING; // for ADL of std functions

        RealType a = dist.location();
        RealType b = dist.scale();
        RealType s = dist.sign();
        RealType c = s * (a - x) / b;
        RealType result;

        if(!detail::verify_scale_b(func, b, &result, Policy())) {
            return result;
        }
        if(!detail::verify_sign(func, s, &result, Policy())) {
            return result;
        }

        // The difference between these two values is interesting. Basically,
        // the cdf for the minimum case is the complement of the max case.
        // Basically, it's 1 - exp(-exp(c)). We could recursively call the
        // complement, but that's a little overkill. If you look at the
        // complement cdf function, you'll see these reversed.
        result = (s == 1)
                ? exp(-exp(c))
                : -boost::math::expm1(-exp(c), Policy());
        return result;
    }

    template <class RealType, class Policy>
    RealType
    quantile(const extreme_value_i_distribution<RealType, Policy>& dist, const RealType& p)
    {
        static const char* func = "boost::math::quantile(const extreme_value_i_distribution<%1%>&, %1%)";

        BOOST_MATH_STD_USING; // for ADL of std functions

        RealType a = dist.location();
        RealType b = dist.scale();
        RealType s = dist.sign();
        RealType result;

        if(!detail::verify_scale_b(func, b, &result, Policy())) {
            return result;
        }
        if(!detail::verify_sign(func, s, &result, Policy())) {
            return result;
        }
        if(!detail::check_probability(func, p, &result, Policy())) {
            return result;
        }

        if(p == 0) {
            return -policies::raise_overflow_error<RealType>(func, 0, Policy());
        }
        else if(p == 1) {
            return policies::raise_overflow_error<RealType>(func, 0, Policy());
        }

        // TODO: Check this out... it might be a - sign * log...
        result = a - log(-log(p)) * b;
        return result;
    }

    template <class RealType, class Policy>
    inline RealType
    cdf(const complemented2_type<extreme_value_i_distribution<RealType, Policy>, RealType>& comp)
    {
        static const char* func = "boost::math::cdf(const extreme_value_i_distribution<%1%>&, %1%)";

        BOOST_MATH_STD_USING; // for ADL of std functions

        RealType a = comp.dist.location();
        RealType b = comp.dist.scale();
        RealType s = comp.dist.sign();
        RealType x = comp.param;
        RealType c = s * (a - x) / b;
        RealType result;

        if(!detail::verify_scale_b(func, b, &result, Policy())) {
            return result;
        }
        if(!detail::verify_sign(func, s, &result, Policy())) {
            return result;
        }

        // See the non-complement CDF implementation for discussion.
        result = (s == 1)
                ? -boost::math::expm1(-exp(c), Policy())
                : exp(-exp(c));
        return result;
    }

    template <class RealType, class Policy>
    RealType
    quantile(const complemented2_type<extreme_value_i_distribution<RealType, Policy>, RealType>& c)
    {
        static const char* func = "boost::math::quantile(const extreme_value_i_distribution<%1%>&, %1%)";

        BOOST_MATH_STD_USING; // for ADL of std functions

        RealType a = c.dist.location();
        RealType b = c.dist.scale();
        RealType s = c.dist.sign();
        RealType q = c.param;
        RealType result;

        if(!detail::verify_scale_b(func, b, &result, Policy())) {
            return result;
        }
        if(!detail::check_probability(func, q, &result, Policy())) {
            return result;
        }
        if(!detail::verify_sign(func, s, &result, Policy())) {
            return result;
        }

        if(q == 0) {
            return policies::raise_overflow_error<RealType>(func, 0, Policy());
        }
        else if(q == 1) {
            return -policies::raise_overflow_error<RealType>(func, 0, Policy());
        }

        // TODO: This is probably not right... There's likely to be some sign-based
        // manipulation going on here. For example, it might be a - sign * log...
        result = a - log(-boost::math::log1p(-q, Policy())) * b;
        return result;
    }

    template <class RealType, class Policy>
    inline RealType
    mean(const extreme_value_i_distribution<RealType, Policy>& dist)
    {
        static const char* func = "boost::math::mean(const extreme_value_i_distribution<%1%>&)";

        RealType a = dist.location();
        RealType b = dist.scale();
        RealType s = dist.sign();
        RealType result;

        if(!detail::verify_scale_b(func, b, &result, Policy())) {
            return result;
        }
        if(!detail::verify_sign(func, s, &result, Policy())) {
            return result;
        }

        return a + s * constants::euler<RealType>() * b;
    }

    template <class RealType, class Policy>
    inline RealType
    standard_deviation(const extreme_value_i_distribution<RealType, Policy>& dist)
    {
        static const char* func = "boost::math::standard_deviation(const extreme_value_i_distribution<%1%>&)";

        BOOST_MATH_STD_USING; // for ADL of std functions.

        RealType b = dist.scale();
        RealType s = dist.sign();
        RealType result;

        if(!detail::verify_scale_b(func, b, &result, Policy())) {
            return result;
        }
        if(!detail::verify_sign(func, s, &result, Policy())) {
            return result;
        }

        return constants::pi<RealType>() * b / sqrt(RealType(6));
    }

    template <class RealType, class Policy>
    inline RealType mode(const extreme_value_i_distribution<RealType, Policy>& dist)
    {
        return dist.location();
    }

    template <class RealType, class Policy>
    inline RealType median(const extreme_value_i_distribution<RealType, Policy>& dist)
    {
        using constants::ln_ln_two;
        return dist.location() - dist.scale() * ln_ln_two<RealType>();
    }

    template <class RealType, class Policy>
    inline RealType skewness(const extreme_value_i_distribution<RealType, Policy>& /*dist*/)
    {
        // This is 12 * sqrt(6) * zeta(3) / pi^3:
        // See http://mathworld.wolfram.com/ExtremeValueDistribution.html
        return RealType(1.1395470994046486574927930193898461120875997958366L);
    }

    template <class RealType, class Policy>
    inline RealType kurtosis(const extreme_value_i_distribution<RealType, Policy>& /*dist*/)
    {
        // See http://mathworld.wolfram.com/ExtremeValueDistribution.html
        return RealType(27) / RealType(5);
    }

    template <class RealType, class Policy>
    inline RealType kurtosis_excess(const extreme_value_i_distribution<RealType, Policy>& /*dist*/)
    {
        // See http://mathworld.wolfram.com/ExtremeValueDistribution.html
        return RealType(12) / RealType(5);
    }


} // namespace math
} // namespace boost

#ifdef BOOST_MSVC
# pragma warning(pop)
#endif

// This include must be at the end, *after* the accessors
// for this distribution have been defined, in order to
// keep compilers that support two-phase lookup happy.
#include <boost/math/distributions/detail/derived_accessors.hpp>

#endif // BOOST_STATS_EXTREME_VALUE_HPP
