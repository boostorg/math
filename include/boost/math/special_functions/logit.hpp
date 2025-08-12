// Copyright Matt Borland 2025.
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SF_LOGIT_HPP
#define BOOST_MATH_SF_LOGIT_HPP

#include <boost/math/tools/config.hpp>
#include <boost/math/policies/policy.hpp>
#include <cmath>
#include <cfenv>

namespace boost {
namespace math {

template <typename RealType, typename Policy>
RealType logit(RealType p, const Policy&)
{
    BOOST_MATH_STD_USING

    using promoted_real_type = typename policies::evaluation<RealType, Policy>::type;

    std::fexcept_t flags;
    std::fegetexceptflag(&flags, FE_ALL_EXCEPT);

    static const RealType crossover {RealType{1}/4};
    const auto promoted_p {static_cast<promoted_real_type>(p)};
    RealType result {};
    if (p > crossover)
    {
        result = 2 * atanh(2 * promoted_p - 1);
    }
    else
    {
        result = log((1 - promoted_p) / promoted_p);
    }

    std::fesetexceptflag(&flags, FE_ALL_EXCEPT);

    return result;
}

template <typename RealType>
RealType logit(RealType p)
{
    return logit(p, policies::policy<>());
}

} // namespace math
} // namespace boost

#endif // BOOST_MATH_SF_LOGIT_HPP
