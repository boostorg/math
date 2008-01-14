// Copyright Edward M. Morrison 2007
// Copyright Andrew Sutton 2007
//
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0. (See accompanying file
// LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SKEW_NORMAL_DISTRIBUTION_HPP
#define BOOST_MATH_SKEW_NORMAL_DISTRIBUTION_HPP

#include <boost/math/distributions/fwd.hpp>
#include <boost/math/distributions/complement.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/distributions/detail/common_error_handling.hpp>

namespace boost { namespace math {
    
namespace skew_normal_detail {
    // delta() := a / sqrt(1 + a^2)
    // expectation() := d * sqrt(2/pi)
    // variance() := 1 - 2*d^2 / pi
}
    

// The skew normal distribution is parameterized over location, shape and
// scale parameters. These are similar to mean and standard deviation in the
// normal distribution, but are scaled and "re-shaped".
//
// Note that we could forseeably cache a number of factors inside the distribution
// like we're doing with the correlation factor. However, these aren't really part
// of the public interface of the type. In fact, it may just be better to recompute
// them on the fly each time. Ask the math guys.
template <typename RealType = double, typename Policy = policies::policy<> >
class skew_normal_distribution
{
public:
    typedef RealType value_type;
    typedef Policy policy_type;
    
    skew_normal_distribution(value_type loc = 0, value_type scale = 1, value_type shape = 0)
        : m_location(loc)
        , m_scale(scale)
        , m_shape(shape)
        , m_correlation(make_correlation(shape))
    {
    }
    
    value_type location() const
    { return m_location; }
    
    value_type scale() const
    { return m_scale; }
    
    value_type shape() const
    { return m_shape; }
    
    value_type correlation() const
    { return m_correlation; }
    
private:
    value_type make_correlation(value_type shape)
    { return shape / sqrt(1.0 + shape * shape); }
    
private:
    value_type m_location;
    value_type m_scale;
    value_type m_shape;
    value_type m_correlation;
};

template <typename RealType, typename Policy>
inline std::pair<RealType, RealType>
range(const skew_normal_distribution<RealType, Policy>&)
{
    return std::make_pair(RealType(), RealType());
}

template <typename RealType, typename Policy>
inline std::pair<RealType, RealType>
support(const skew_normal_distribution<RealType, Policy>&)
{
    return std::make_pair(RealType(), RealType());
}

template <typename RealType, typename Policy>
inline RealType
pdf(const skew_normal_distribution<RealType, Policy>& dist, const RealType& x)
{
    return RealType();
}

template <typename RealType, typename Policy>
inline RealType
cdf(const skew_normal_distribution<RealType, Policy>& dist, const RealType& x)
{
    return RealType();
}

template <typename RealType, typename Policy>
inline RealType
cdf(const complemented2_type<skew_normal_distribution<RealType, Policy>, RealType>& dist)
{
    return RealType();
}

template <typename RealType, typename Policy>
inline RealType
quantile(const skew_normal_distribution<RealType, Policy>& dist, const RealType& x)
{
    return RealType();
}

template <typename RealType, typename Policy>
inline RealType
quantile(const complemented2_type<skew_normal_distribution<RealType, Policy>, RealType>& dist)
{
    return RealType();
}

template <typename RealType, typename Policy>
inline RealType 
mean(const skew_normal_distribution<RealType, Policy>& dist)
{
    // This is the complicated way of writing dist.location()
    RealType var = sqrt(dist.scale());
    RealType root = sqrt(2.0 / constants::pi<RealType>());
    return  dist.location() + var * root * dist.correlation();
}

template <typename RealType, typename Policy>
inline RealType
variance(const skew_normal_distribution<RealType, Policy>& dist)
{
    // This is apparently the same as dist.scale()^2.
    RealType scale = dist.scale() * dist.scale();
    RealType corr =  dist.correlation() * dist.correlation();
    return scale * (1.0 - 2.0 * corr / constants::pi<RealType>());
}

template <typename RealType, typename Policy>
inline RealType
skewness(const skew_normal_distribution<RealType, Policy>& dist)
{
    return RealType();
}

template <typename RealType, typename Policy>
inline RealType
kurtosis(const skew_normal_distribution<RealType, Policy>& dist)
{
    return RealType();
}

template <typename RealType, typename Policy>
inline RealType
kurtosis_excess(const skew_normal_distribution<RealType, Policy>& dist)
{
    return RealType();
}
    
} }

#endif
