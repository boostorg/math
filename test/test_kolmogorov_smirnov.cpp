// Copyright Evan Miller 2020
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
//
#include <pch_light.hpp>
#include <boost/math/concepts/real_concept.hpp>

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp> // for test_main
#include <boost/test/tools/floating_point_comparison.hpp> // for BOOST_CHECK_CLOSE
#include <boost/math/distributions/kolmogorov_smirnov.hpp>
#include <boost/math/quadrature/exp_sinh.hpp>
#include <boost/math/special_functions/zeta.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>

template <typename RealType> // Any floating-point type RealType.
void test_spots(RealType)
{
    using namespace boost::math;
    // Test quantiles, CDFs, and complements
    RealType eps = tools::epsilon<RealType>();
    RealType tol = tools::epsilon<RealType>() * 25;
    for (int n=10; n<100; n += 10) {
        kolmogorov_smirnov_distribution<RealType> dist(n);
        for (int i=0; i<1000; i++) {
            RealType p = 1.0 * (i+1) / 1001;
            RealType crit1 = quantile(dist, 1 - p);
            RealType crit2 = quantile(complement(dist, p));
            RealType p1 = cdf(dist, crit1);
            BOOST_CHECK_CLOSE_FRACTION(crit1, crit2, tol);
            BOOST_CHECK_CLOSE_FRACTION(1 - p, p1, tol);
        }

        for (int i=0; i<1000; i++) {
            RealType x = 1.0 * (i+1) / 1001;
            RealType p = cdf(dist, x);
            RealType p1 = cdf(complement(dist, x));
            if (p > tol && p1 > tol) { // skip the extreme tails
                RealType x1;
                if (p < 0.5)
                    x1 = quantile(dist, p);
                else
                    x1 = quantile(complement(dist, p1));
                BOOST_CHECK_CLOSE_FRACTION(x, x1, tol);
            }
        }
    }

    kolmogorov_smirnov_distribution<RealType> dist(100);

    // Basics
    BOOST_CHECK_THROW(pdf(dist, RealType(-1.0)), std::domain_error);
    BOOST_CHECK_THROW(cdf(dist, RealType(-1.0)), std::domain_error);
    BOOST_CHECK_THROW(quantile(dist, RealType(-1.0)), std::domain_error);
    BOOST_CHECK_THROW(quantile(dist, RealType(2.0)), std::domain_error);

    // Confirm mode is at least a local minimum
    RealType mode = boost::math::mode(dist);

    using std::sqrt;
    BOOST_TEST_CHECK(pdf(dist, mode) >= pdf(dist, RealType(mode - sqrt(eps))));
    BOOST_TEST_CHECK(pdf(dist, mode) >= pdf(dist, RealType(mode + sqrt(eps))));

    // Test the moments - each one integrates the entire distribution
    quadrature::exp_sinh<RealType> integrator;

    auto f_one = [&, dist](RealType t) { return pdf(dist, t); };
    BOOST_CHECK_CLOSE_FRACTION(integrator.integrate(f_one, eps), RealType(1), tol);

    RealType mean = boost::math::mean(dist);
    auto f_mean = [&, dist](RealType t) { return pdf(dist, t) * t; };
    BOOST_CHECK_CLOSE_FRACTION(integrator.integrate(f_mean, eps), mean, tol);

    RealType var = variance(dist);
    auto f_var = [&, dist, mean](RealType t) { return pdf(dist, t) * (t - mean) * (t - mean); };
    BOOST_CHECK_CLOSE_FRACTION(integrator.integrate(f_var, eps), var, tol);

    RealType skew = skewness(dist);
    auto f_skew = [&, dist, mean, var](RealType t) { return pdf(dist, t)
        * (t - mean) * (t - mean) * (t - mean) / var / sqrt(var); };
    BOOST_CHECK_CLOSE_FRACTION(integrator.integrate(f_skew, eps), skew, 10*tol);

    RealType kurt = kurtosis(dist);
    auto f_kurt= [&, dist, mean, var](RealType t) { return pdf(dist, t)
        * (t - mean) * (t - mean) * (t - mean) * (t - mean) / var / var; };
    BOOST_CHECK_CLOSE_FRACTION(integrator.integrate(f_kurt, eps), kurt, 5*tol);

    BOOST_CHECK_CLOSE_FRACTION(kurt, kurtosis_excess(dist) + 3, eps);
}

// Checks that x is the quantile of p to within a few ulps, by bracketing:
// cdf(x * (1 - margin)) <= p <= cdf(x * (1 + margin)). Unlike a round-trip
// comparison of cdf(quantile(p)) with p, this remains meaningful deep in the
// tails, where the CDF is so steep that one ulp of x corresponds to many
// ulps of p.
template <typename RealType>
void check_quantile(const boost::math::kolmogorov_smirnov_distribution<RealType>& dist, RealType p, RealType x, bool complemented)
{
    using namespace boost::math;
    RealType margin = 8 * tools::epsilon<RealType>();
    RealType lo = complemented ? cdf(complement(dist, RealType(x * (1 + margin)))) : cdf(dist, RealType(x * (1 - margin)));
    RealType hi = complemented ? cdf(complement(dist, RealType(x * (1 - margin)))) : cdf(dist, RealType(x * (1 + margin)));
    BOOST_TEST_CHECK(lo <= p);
    BOOST_TEST_CHECK(p <= hi);
}

template <typename RealType>
void test_quantile_extremes(RealType)
{
    using namespace boost::math;
    RealType eps = tools::epsilon<RealType>();

    // Small n, where the quantile exceeds unity in the upper tail
    for (RealType n : { RealType(1), RealType(2), RealType(3), RealType(1000000) }) {
        kolmogorov_smirnov_distribution<RealType> dist(n);
        for (int i = 1; i < 100; i++) {
            RealType p = RealType(i) / 100;
            check_quantile(dist, p, quantile(dist, p), false);
            check_quantile(dist, p, quantile(complement(dist, p)), true);
        }
    }

    kolmogorov_smirnov_distribution<RealType> dist(10);

    // Deep tails, down to where the CDF itself is no longer evaluated to
    // full precision because its terms become denormal
    for (RealType p = RealType(0.1); p > 10000 * tools::min_value<RealType>(); p /= 1000) {
        check_quantile(dist, p, quantile(dist, p), false);
        check_quantile(dist, p, quantile(complement(dist, p)), true);
    }

    // Limits of the probability range
    BOOST_CHECK_EQUAL(quantile(dist, RealType(0)), RealType(0));
    BOOST_CHECK_EQUAL(quantile(complement(dist, RealType(1))), RealType(0));
    BOOST_CHECK_THROW(quantile(dist, RealType(1)), std::overflow_error);
    BOOST_CHECK_THROW(quantile(complement(dist, RealType(0))), std::overflow_error);

    // The median is the quantile of one half
    BOOST_CHECK_CLOSE_FRACTION(cdf(dist, median(dist)), RealType(0.5), 10 * eps);
    BOOST_CHECK_CLOSE_FRACTION(median(dist), quantile(dist, RealType(0.5)), 10 * eps);
}

// Reference values computed with the same series at 120 decimal digits.
template <typename RealType>
void test_reference_values(RealType)
{
    using namespace boost::math;
    BOOST_MATH_STD_USING // without it sqrt(long double) can bind to ::sqrt(double)
    RealType eps = tools::epsilon<RealType>();
    static const struct { double n, x; long double pdf, cdf, ccdf; } data[] = {
        { 10, 0.125, 0.279400023454107971738809672357L, 0.00236117478750956819598968016445L, 0.997638825212490431804010319836L },
        { 10, 0.25,  5.19340083129894262893339472019L,  0.440440289804736223158552770951L,  0.559559710195263776841447229049L },
        { 10, 0.375, 1.80007916388412209267020060796L,  0.879916678784370975806015389694L,  0.120083321215629024193984610306L },
        { 10, 0.5,   0.269517550178839104001256352353L, 0.986524110124136310626593187475L,  0.0134758898758636893734068125254L },
        { 10, 0.75,  0.000780437859244050388721178490249L, 0.999973985404691864815291512551L, 2.60145953081351847084874493052e-05L },
        { 10, 1.0,   1.64892289795084626237275224637e-07L, 0.999999995877692755122884344068L, 4.12230724487711565593188072421e-09L },
        { 1,  0.25,  4.12855106473622609508308670705e-06L, 2.68238100848298275381058297367e-08L, 0.999999973176189915170172461894L },
        { 1,  0.5,   0.639582850940456634645459338412L, 0.0360547563351249056140861037179L, 0.963945243664875094385913896282L },
        { 1,  0.75,  1.6834609513049753816942330657L,   0.372832958223738358506334471397L,  0.627167041776261641493665528603L },
        { 1,  1.5,   0.133307227419880210037320465613L, 0.97778203738347487127945638539L,   0.0222179626165251287205436146105L },
        { 1,  2.5,   7.45330634415734044284982458077e-05L, 0.999992546693655842658399900267L, 7.45330634415734160009973335912e-06L },
    };
    for (const auto& d : data) {
        kolmogorov_smirnov_distribution<RealType> dist(static_cast<RealType>(d.n));
        RealType x = static_cast<RealType>(d.x);
        BOOST_CHECK_CLOSE_FRACTION(pdf(dist, x), static_cast<RealType>(d.pdf), 20 * eps);
        BOOST_CHECK_CLOSE_FRACTION(cdf(dist, x), static_cast<RealType>(d.cdf), 50 * eps);
        // Above 2*x*x*n = pi the complement is computed from the nome with a
        // compensated exponent; below it goes through the Jacobi Theta
        // function with a rounded tau, whose error is amplified by the
        // exponent.
        RealType ccdf_tol = (2 * x * x * d.n > constants::pi<double>()) ? 10 * eps : 200 * eps;
        BOOST_CHECK_CLOSE_FRACTION(cdf(complement(dist, x)), static_cast<RealType>(d.ccdf), ccdf_tol);
    }

    // The moment constants agree with the closed forms they were computed from
    kolmogorov_smirnov_distribution<RealType> dist(10);
    RealType n = 10;
    RealType mean = boost::math::mean(dist);
    RealType var = variance(dist);
    RealType ex3 = RealType(0.5625) * constants::root_half_pi<RealType>() * constants::zeta_three<RealType>() / n / sqrt(n);
    RealType ex4 = 7 * constants::pi_sqr_div_six<RealType>() * constants::pi_sqr_div_six<RealType>() / 20 / n / n;
    RealType skew = (ex3 - 3 * mean * var - mean * mean * mean) / var / sqrt(var);
    RealType kurt = (ex4 - 4 * mean * skew * var * sqrt(var) - 6 * mean * mean * var - mean * mean * mean * mean) / var / var;
    BOOST_CHECK_CLOSE_FRACTION(skewness(dist), skew, 500 * eps);
    BOOST_CHECK_CLOSE_FRACTION(kurtosis(dist), kurt, 500 * eps);
    BOOST_CHECK_CLOSE_FRACTION(kurtosis_excess(dist), kurt - 3, 2000 * eps);
}

// Types with more than 100 decimal digits compute the mode, median and
// moments at run-time; check that they agree with the tabulated constants.
void test_high_precision()
{
    using namespace boost::math;
    typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float<110> > RealType;
    kolmogorov_smirnov_distribution<RealType> dist(1);
    RealType tol = RealType("1e-99");
    BOOST_CHECK_CLOSE_FRACTION(mode(dist), RealType("0.7354679079165719820624448513051825390913503143340225122926732679283442751250018638342479353073554823"), tol);
    BOOST_CHECK_CLOSE_FRACTION(median(dist), RealType("0.8275735551899076901138270828889768075843727832232029452002281054669822747682290224084911530683945353"), tol);
    BOOST_CHECK_CLOSE_FRACTION(skewness(dist), RealType("0.8604261371436682558667183685173007452393495203060929924698665373211030237490355834020123201092217303"), tol);
    BOOST_CHECK_CLOSE_FRACTION(kurtosis_excess(dist), RealType("0.8816189679105236704015538304328295408355055629289034536289160037385510119639244078254098968815665536"), tol);
}

BOOST_AUTO_TEST_CASE( test_main )
{
  BOOST_MATH_CONTROL_FP;

  // (Parameter value, arbitrarily zero, only communicates the floating point type).
  test_spots(0.0F); // Test float.
  test_spots(0.0); // Test double.
  test_quantile_extremes(0.0F);
  test_quantile_extremes(0.0);
  test_reference_values(0.0F);
  test_reference_values(0.0);
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
  test_spots(0.0L); // Test long double.
  test_quantile_extremes(0.0L);
  test_reference_values(0.0L);
#if !defined(BOOST_MATH_NO_REAL_CONCEPT_TESTS)
  test_spots(boost::math::concepts::real_concept(0.)); // Test real concept.
  test_quantile_extremes(boost::math::concepts::real_concept(0.));
  test_reference_values(boost::math::concepts::real_concept(0.));
#endif
#endif
  test_high_precision();
}
