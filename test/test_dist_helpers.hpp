/*
 * Copyright Jacob Hass, 2026
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef BOOST_MATH_TEST_DIST_HELPERS_HPP
#define BOOST_MATH_TEST_DIST_HELPERS_HPP

#include <utility>
#include <limits>
#include <stdexcept>
#include <boost/math/tools/precision.hpp>
#include <boost/math/special_functions/next.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/math/policies/policy.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#if defined(BOOST_CHECK_THROW) && defined(BOOST_MATH_NO_EXCEPTIONS)
#  undef BOOST_CHECK_THROW
#  define BOOST_CHECK_THROW(x, y)
#endif

template <class Dist, class Container>
Dist make_distribution(const Container& c)
{
    using value_type = typename Dist::value_type;
    BOOST_MATH_IF_CONSTEXPR(std::is_constructible<Dist, value_type, value_type, value_type>::value)
    {
        if (c.size() >= 3)
            return Dist(c.data()[0], c.data()[1], c.data()[2]);
    }
    BOOST_MATH_IF_CONSTEXPR(std::is_constructible<Dist, value_type, value_type>::value)
    {
        if (c.size() >= 2)
            return Dist(c.data()[0], c.data()[1]);
    }
    BOOST_MATH_IF_CONSTEXPR(std::is_constructible<Dist, value_type>::value)
    {
        if (c.size() >= 1)
            return Dist(c.data()[0]);
    }
    throw std::domain_error("Object not initialized!");
}

template <class Dist, class V>
Dist make_distribution(const std::initializer_list<V>& c)
{
    return make_distribution<Dist, std::initializer_list<V>>(c);
}

template <class Dist, class Ignore_Error_Dist, class Real>
void test_invalid_support(std::vector<Real> params)
{
    using namespace boost::math;
    
    Dist dist = make_distribution<Dist>(params); 
    Ignore_Error_Dist ignore_error_dist = make_distribution<Ignore_Error_Dist>(params);
    
    /* We will assume that std::numeric_limits<Distribution::value_type>::has_infinity */
    std::pair<Real, Real> sup = support(dist);

    // Test outside lower bound:
    Real invalid = sup.first;
    if (boost::math::isfinite(invalid))
    {
        if (invalid == -boost::math::tools::max_value<Real>())
            invalid = -std::numeric_limits<Real>::infinity();
        else if (invalid == boost::math::tools::min_value<Real>())
        {
            invalid = -boost::math::tools::min_value<Real>();
        }
        else
            invalid = boost::math::float_prior(invalid);

        // Test PDF/CDF
        BOOST_CHECK_THROW(pdf(dist, invalid), std::domain_error);
        BOOST_CHECK_THROW(cdf(dist, invalid), std::domain_error);
        BOOST_CHECK_THROW(cdf(complement(dist, invalid)), std::domain_error);
        BOOST_CHECK((boost::math::isnan)(pdf(ignore_error_dist, invalid)));
        BOOST_CHECK((boost::math::isnan)(cdf(ignore_error_dist, invalid)));
        BOOST_CHECK((boost::math::isnan)(cdf(complement(ignore_error_dist, invalid))));

        // Test log PDF/CDF
        BOOST_CHECK_THROW(logpdf(dist, invalid), std::domain_error);
        BOOST_CHECK_THROW(logcdf(dist, invalid), std::domain_error);
        BOOST_CHECK_THROW(logcdf(complement(dist, invalid)), std::domain_error);
        BOOST_CHECK((boost::math::isnan)(logpdf(ignore_error_dist, invalid)));
        BOOST_CHECK((boost::math::isnan)(logcdf(ignore_error_dist, invalid)));
        BOOST_CHECK((boost::math::isnan)(logcdf(complement(ignore_error_dist, invalid))));
    }
    else // lower bound is -infinity
    {
        BOOST_CHECK_THROW(quantile(dist, static_cast<Real>(0.0)), std::overflow_error);
        BOOST_CHECK_THROW(quantile(complement(dist, static_cast<Real>(1.0))), std::overflow_error);
        BOOST_CHECK_EQUAL(pdf(dist, invalid), static_cast<Real>(0));
        BOOST_CHECK_EQUAL(cdf(dist, invalid), static_cast<Real>(0));
        BOOST_CHECK_EQUAL(cdf(complement(dist, invalid)), static_cast<Real>(1));
    }
    // Test outside upper bound:
    invalid = sup.second;
    if (boost::math::isfinite(invalid))
    {
        if (invalid == boost::math::tools::max_value<Real>())
            invalid = std::numeric_limits<Real>::infinity();
        else
            invalid = boost::math::float_next(invalid);
        // Test PDF/CDF
        BOOST_CHECK_THROW(pdf(dist, invalid), std::domain_error);
        BOOST_CHECK_THROW(cdf(dist, invalid), std::domain_error);
        BOOST_CHECK_THROW(cdf(complement(dist, invalid)), std::domain_error);
        // Test NaN handling
        BOOST_CHECK((boost::math::isnan)(pdf(ignore_error_dist, invalid)));
        BOOST_CHECK((boost::math::isnan)(cdf(ignore_error_dist, invalid)));
        BOOST_CHECK((boost::math::isnan)(cdf(complement(ignore_error_dist, invalid))));

        // Test log PDF/CDF
        BOOST_CHECK_THROW(logpdf(dist, invalid), std::domain_error);
        BOOST_CHECK_THROW(logcdf(dist, invalid), std::domain_error);
        BOOST_CHECK_THROW(logcdf(complement(dist, invalid)), std::domain_error);
        BOOST_CHECK((boost::math::isnan)(logpdf(ignore_error_dist, invalid)));
        BOOST_CHECK((boost::math::isnan)(logcdf(ignore_error_dist, invalid)));
        BOOST_CHECK((boost::math::isnan)(logcdf(complement(ignore_error_dist, invalid))));
    }
    else // upper bound is +infinity
    {
        BOOST_CHECK_EQUAL(pdf(dist, invalid), static_cast<Real>(0));
        BOOST_CHECK_EQUAL(cdf(dist, invalid), static_cast<Real>(1));
        BOOST_CHECK_EQUAL(cdf(complement(dist, invalid)), static_cast<Real>(0));

        // If support is infinite quantile will throw overflow error
        BOOST_CHECK_THROW(quantile(dist, static_cast<Real>(1.0)), std::overflow_error);
        BOOST_CHECK_THROW(quantile(complement(dist, static_cast<Real>(0.0))), std::overflow_error);
    }

    // NaN checks for PDF/CDF
    invalid = std::numeric_limits<Real>::quiet_NaN();
    BOOST_CHECK_THROW(pdf(dist, invalid), std::domain_error);
    BOOST_CHECK_THROW(cdf(dist, invalid), std::domain_error);
    BOOST_CHECK_THROW(cdf(complement(dist, invalid)), std::domain_error);
    BOOST_CHECK((boost::math::isnan)(pdf(ignore_error_dist, invalid)));
    BOOST_CHECK((boost::math::isnan)(cdf(ignore_error_dist, invalid)));
    BOOST_CHECK((boost::math::isnan)(cdf(complement(ignore_error_dist, invalid))));

    // NaN checks for log PDF/CDF
    BOOST_CHECK_THROW(logpdf(dist, invalid), std::domain_error);
    BOOST_CHECK_THROW(logcdf(dist, invalid), std::domain_error);
    BOOST_CHECK_THROW(logcdf(complement(dist, invalid)), std::domain_error);
    BOOST_CHECK((boost::math::isnan)(logpdf(ignore_error_dist, invalid)));
    BOOST_CHECK((boost::math::isnan)(logcdf(ignore_error_dist, invalid)));
    BOOST_CHECK((boost::math::isnan)(logcdf(complement(ignore_error_dist, invalid))));

    // Quantile Checks
    BOOST_CHECK_THROW(quantile(dist, std::numeric_limits<Real>::infinity()), std::domain_error); // p == infinity
    BOOST_CHECK_THROW(quantile(dist, static_cast<Real>(-1)), std::domain_error); // p < 0
    BOOST_CHECK_THROW(quantile(dist, static_cast<Real>(2)), std::domain_error); // p > 1
    BOOST_CHECK_THROW(quantile(complement(dist, std::numeric_limits<Real>::infinity())), std::domain_error); // q == infinity
    BOOST_CHECK_THROW(quantile(complement(dist, static_cast<Real>(-1))), std::domain_error); // q < 0
    BOOST_CHECK_THROW(quantile(complement(dist, static_cast<Real>(2))), std::domain_error); // q > 1

    BOOST_CHECK((boost::math::isnan)(quantile(ignore_error_dist, std::numeric_limits<Real>::infinity()))); // p == infinity
    BOOST_CHECK((boost::math::isnan)(quantile(ignore_error_dist, static_cast<Real>(-1)))); // p < 0
    BOOST_CHECK((boost::math::isnan)(quantile(ignore_error_dist, static_cast<Real>(2)))); // p > 1
    BOOST_CHECK((boost::math::isnan)(quantile(complement(ignore_error_dist, std::numeric_limits<Real>::infinity())))); // q == infinity
    BOOST_CHECK((boost::math::isnan)(quantile(complement(ignore_error_dist, static_cast<Real>(-1))))); // q < 0
    BOOST_CHECK((boost::math::isnan)(quantile(complement(ignore_error_dist, static_cast<Real>(2))))); // q > 1
}

template <class Dist, class Ignore_Error_Dist, class Real>
void test_invalid_parameters(std::vector<std::vector<Real> > invalid_parameters)
{
    std::vector<Real> params;
    for (unsigned i=0; i<invalid_parameters.size(); i++)
    {
        params = invalid_parameters[i];
        // Check functions throw with bad constructors
        BOOST_CHECK_THROW(make_distribution<Dist>(params), std::domain_error);
        BOOST_CHECK_THROW(pdf(make_distribution<Dist>(params), static_cast<Real>(0)), std::domain_error);
        BOOST_CHECK_THROW(logpdf(make_distribution<Dist>(params), static_cast<Real>(0)), std::domain_error);
        BOOST_CHECK_THROW(cdf(make_distribution<Dist>(params), static_cast<Real>(0)), std::domain_error);
        BOOST_CHECK_THROW(cdf(complement(make_distribution<Dist>(params), static_cast<Real>(0))), std::domain_error);
        BOOST_CHECK_THROW(logcdf(make_distribution<Dist>(params), static_cast<Real>(0)), std::domain_error);
        BOOST_CHECK_THROW(logcdf(complement(make_distribution<Dist>(params), static_cast<Real>(0))), std::domain_error);
        BOOST_CHECK_THROW(quantile(make_distribution<Dist>(params), static_cast<Real>(0.5)), std::domain_error);
        BOOST_CHECK_THROW(quantile(complement(make_distribution<Dist>(params), static_cast<Real>(0.5))), std::domain_error);

        // Check return NaN
        if (std::numeric_limits<Real>::has_quiet_NaN)
        {
            BOOST_CHECK((boost::math::isnan)(pdf(make_distribution<Ignore_Error_Dist>(params), static_cast<Real>(0))));
            BOOST_CHECK((boost::math::isnan)(logpdf(make_distribution<Ignore_Error_Dist>(params), static_cast<Real>(0))));
            BOOST_CHECK((boost::math::isnan)(cdf(make_distribution<Ignore_Error_Dist>(params), static_cast<Real>(0))));
            BOOST_CHECK((boost::math::isnan)(cdf(complement(make_distribution<Ignore_Error_Dist>(params), static_cast<Real>(0)))));
            BOOST_CHECK((boost::math::isnan)(logcdf(make_distribution<Ignore_Error_Dist>(params), static_cast<Real>(0))));
            BOOST_CHECK((boost::math::isnan)(logcdf(complement(make_distribution<Ignore_Error_Dist>(params), static_cast<Real>(0)))));
            BOOST_CHECK((boost::math::isnan)(quantile(make_distribution<Ignore_Error_Dist>(params), static_cast<Real>(0.5))));
            BOOST_CHECK((boost::math::isnan)(quantile(complement(make_distribution<Ignore_Error_Dist>(params), static_cast<Real>(0.5)))));
        }
    }
}

#endif // BOOST_MATH_TEST_DIST_HELPERS_HPP