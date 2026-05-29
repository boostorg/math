#include <utility>
#include <limits>
#include <stdexcept>
#include <boost/math/tools/precision.hpp>
#include <boost/math/special_functions/next.hpp>
#include <boost/test/unit_test.hpp>

template <template <typename, typename> typename Distribution, class Real>
void test_invalid_support()
{
    using boost::math::policies::policy;
    typedef policy<
      boost::math::policies::domain_error<boost::math::policies::ignore_error>,
      boost::math::policies::overflow_error<boost::math::policies::ignore_error>,
      boost::math::policies::underflow_error<boost::math::policies::ignore_error>,
      boost::math::policies::denorm_error<boost::math::policies::ignore_error>,
      boost::math::policies::pole_error<boost::math::policies::ignore_error>,
      boost::math::policies::evaluation_error<boost::math::policies::ignore_error>
    > ignore_all_policy;

    Distribution<Real, ignore_all_policy> ignore_error_dist;
    Distribution<Real, policy<> > dist; 
    
    /* We will assume that std::numeric_limits<Distribution::value_type>::has_infinity */
    std::pair<Real, Real> sup = support(dist);

    // Test outside lower bound:
    Real invalid = sup.first;
    if (boost::math::isfinite(invalid))
    {
        if (invalid == -boost::math::tools::max_value<Real>())
            invalid = -std::numeric_limits<Real>::infinity();
        else
            invalid = boost::math::float_prior(invalid);

        // Test PDF/CDF
        BOOST_MATH_CHECK_THROW(pdf(dist, invalid), std::domain_error);
        BOOST_MATH_CHECK_THROW(cdf(dist, invalid), std::domain_error);
        BOOST_MATH_CHECK_THROW(cdf(complement(dist, invalid)), std::domain_error);
        BOOST_CHECK((boost::math::isnan)(pdf(ignore_error_dist, invalid)));
        BOOST_CHECK((boost::math::isnan)(cdf(ignore_error_dist, invalid)));
        BOOST_CHECK((boost::math::isnan)(cdf(complement(ignore_error_dist, invalid))));

        // Test log PDF/CDF
        BOOST_MATH_CHECK_THROW(logpdf(dist, invalid), std::domain_error);
        BOOST_MATH_CHECK_THROW(logcdf(dist, invalid), std::domain_error);
        BOOST_MATH_CHECK_THROW(logcdf(complement(dist, invalid)), std::domain_error);
        BOOST_CHECK((boost::math::isnan)(logpdf(ignore_error_dist, invalid)));
        BOOST_CHECK((boost::math::isnan)(logcdf(ignore_error_dist, invalid)));
        BOOST_CHECK((boost::math::isnan)(logcdf(complement(ignore_error_dist, invalid))));
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
        BOOST_MATH_CHECK_THROW(pdf(dist, invalid), std::domain_error);
        BOOST_MATH_CHECK_THROW(cdf(dist, invalid), std::domain_error);
        BOOST_MATH_CHECK_THROW(cdf(complement(dist, invalid)), std::domain_error);
        // Test NaN handnling
        BOOST_CHECK((boost::math::isnan)(pdf(ignore_error_dist, invalid)));
        BOOST_CHECK((boost::math::isnan)(cdf(ignore_error_dist, invalid)));
        BOOST_CHECK((boost::math::isnan)(cdf(complement(ignore_error_dist, invalid))));

        // Test log PDF/CDF
        BOOST_MATH_CHECK_THROW(logpdf(dist, invalid), std::domain_error);
        BOOST_MATH_CHECK_THROW(logcdf(dist, invalid), std::domain_error);
        BOOST_MATH_CHECK_THROW(logcdf(complement(dist, invalid)), std::domain_error);
        BOOST_CHECK((boost::math::isnan)(logpdf(ignore_error_dist, invalid)));
        BOOST_CHECK((boost::math::isnan)(logcdf(ignore_error_dist, invalid)));
        BOOST_CHECK((boost::math::isnan)(logcdf(complement(ignore_error_dist, invalid))));
    }

    // NaN checks for PDF/CDF
    invalid = std::numeric_limits<Real>::quiet_NaN();
    BOOST_MATH_CHECK_THROW(pdf(dist, invalid), std::domain_error);
    BOOST_MATH_CHECK_THROW(cdf(dist, invalid), std::domain_error);
    BOOST_MATH_CHECK_THROW(cdf(complement(dist, invalid)), std::domain_error);
    BOOST_CHECK((boost::math::isnan)(pdf(ignore_error_dist, invalid)));
    BOOST_CHECK((boost::math::isnan)(cdf(ignore_error_dist, invalid)));
    BOOST_CHECK((boost::math::isnan)(cdf(complement(ignore_error_dist, invalid))));

    // NaN checks for log PDF/CDF
    invalid = std::numeric_limits<Real>::quiet_NaN();
    BOOST_MATH_CHECK_THROW(logpdf(dist, invalid), std::domain_error);
    BOOST_MATH_CHECK_THROW(logcdf(dist, invalid), std::domain_error);
    BOOST_MATH_CHECK_THROW(logcdf(complement(dist, invalid)), std::domain_error);
    BOOST_CHECK((boost::math::isnan)(logpdf(ignore_error_dist, invalid)));
    BOOST_CHECK((boost::math::isnan)(logcdf(ignore_error_dist, invalid)));
    BOOST_CHECK((boost::math::isnan)(logcdf(complement(ignore_error_dist, invalid))));

    // Quantile Checks
    BOOST_MATH_CHECK_THROW(quantile(dist, std::numeric_limits<Real>::infinity()), std::domain_error); // p == infinity
    BOOST_MATH_CHECK_THROW(quantile(dist, static_cast<Real>(-1)), std::domain_error); // p < 0
    BOOST_MATH_CHECK_THROW(quantile(dist, static_cast<Real>(2)), std::domain_error); // p > 1
    BOOST_MATH_CHECK_THROW(quantile(complement(dist, std::numeric_limits<Real>::infinity())), std::domain_error); // q == infinity
    BOOST_MATH_CHECK_THROW(quantile(complement(dist, static_cast<Real>(-1))), std::domain_error); // q < 0
    BOOST_MATH_CHECK_THROW(quantile(complement(dist, static_cast<Real>(2))), std::domain_error); // q > 1
    BOOST_CHECK_EQUAL(quantile(dist, static_cast<Real>(1)), static_cast<Real>(1)); // p = 1
    BOOST_CHECK_EQUAL(quantile(dist, static_cast<Real>(0)), static_cast<Real>(0)); // p = 0
    BOOST_CHECK_EQUAL(quantile(complement(dist, static_cast<Real>(1))), static_cast<Real>(0)); // q = 1
    BOOST_CHECK_EQUAL(quantile(complement(dist, static_cast<Real>(0))), static_cast<Real>(1)); // q = 0

    BOOST_CHECK((boost::math::isnan)(quantile(ignore_error_dist, std::numeric_limits<Real>::infinity()))); // p == infinity
    BOOST_CHECK((boost::math::isnan)(quantile(ignore_error_dist, static_cast<Real>(-1)))); // p < 0
    BOOST_CHECK((boost::math::isnan)(quantile(ignore_error_dist, static_cast<Real>(2)))); // p > 1
    BOOST_CHECK((boost::math::isnan)(quantile(complement(ignore_error_dist, std::numeric_limits<Real>::infinity())))); // q == infinity
    BOOST_CHECK((boost::math::isnan)(quantile(complement(ignore_error_dist, static_cast<Real>(-1))))); // q < 0
    BOOST_CHECK((boost::math::isnan)(quantile(complement(ignore_error_dist, static_cast<Real>(2))))); // q > 1
}