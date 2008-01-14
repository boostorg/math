// Copyright Andrew Sutton 2007
//
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0. (See accompanying file
// LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_GEOMETRIC_DISTRIBUTION_HPP
#define BOOST_MATH_GEOMETRIC_DISTRIBUTION_HPP

#include <boost/math/distributions/fwd.hpp>
#include <boost/math/distributions/complement.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/distributions/detail/common_error_handling.hpp>

namespace boost { namespace math {

    // The geometric distribution... Write docs here...
    template <typename RealType = double, typename Policy = policies::policy<> >
    class geometric_distribution
    {
    public:
        typedef RealType value_type;
        typedef Policy policy_type;

        geometric_distribution(value_type prob)
            : m_prob(prob)

        {
        }

        value_type success_probability() const
        { return m_prob; }

    private:
        value_type m_prob;
    };

    template <typename RealType, typename Policy>
    inline std::pair<RealType, RealType>
    range(const geometric_distribution<RealType, Policy>&)
    {
        return std::make_pair(RealType(), RealType());
    }

    template <typename RealType, typename Policy>
    inline std::pair<RealType, RealType>
    support(const geometric_distribution<RealType, Policy>&)
    {
        return std::make_pair(RealType(), RealType());
    }

    template <typename RealType, typename Policy>
    inline RealType
    pdf(const geometric_distribution<RealType, Policy>& dist, const RealType& x)
    {
        using std::pow;
        RealType one = 1;
        RealType p = dist.success_probability();
        return pow(one - p, x - one) * p;
    }

    template <typename RealType, typename Policy>
    inline RealType
    cdf(const geometric_distribution<RealType, Policy>& dist, const RealType& x)
    {
        using std::pow;
        RealType one = 1;
        RealType p = dist.success_probability();
        return one - pow(one - p, x);
    }

    template <typename RealType, typename Policy>
    inline RealType
    cdf(const complemented2_type<geometric_distribution<RealType, Policy>, RealType>& dist)
    {
        return RealType();
    }

    template <typename RealType, typename Policy>
    inline RealType
    quantile(const geometric_distribution<RealType, Policy>& dist, const RealType& x)
    {
        return RealType();
    }

    template <typename RealType, typename Policy>
    inline RealType
    quantile(const complemented2_type<geometric_distribution<RealType, Policy>, RealType>& dist)
    {
        return RealType();
    }

    template <typename RealType, typename Policy>
    inline RealType
    mean(const geometric_distribution<RealType, Policy>& dist)
    {
        RealType one = 1;
        RealType p = dist.success_probability();
        return one / p;
    }

    template <typename RealType, typename Policy>
    inline RealType
    variance(const geometric_distribution<RealType, Policy>& dist)
    {
        RealType one = 1;
        RealType p = dist.success_probability();
        return (one - p) / (p * p);
    }

    template <typename RealType, typename Policy>
    inline RealType
    skewness(const geometric_distribution<RealType, Policy>& dist)
    {
        using std::sqrt;
        RealType one = 1;
        RealType two = 2;
        RealType p = dist.success_probability();
        return (two - p) / sqrt(one - p);
    }

    template <typename RealType, typename Policy>
    inline RealType
    kurtosis(const geometric_distribution<RealType, Policy>& dist)
    {
        return RealType();
    }

    template <typename RealType, typename Policy>
    inline RealType
    kurtosis_excess(const geometric_distribution<RealType, Policy>& dist)
    {
        return RealType();
    }

} }

#endif
