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

// Reference values computed from the two series at 60 decimal digits and given
// to 40, which is enough for a 113-bit long double.
template <typename RealType>
void test_reference_values(RealType)
{
    using namespace boost::math;
    BOOST_MATH_STD_USING // without it sqrt(long double) can bind to ::sqrt(double)
    RealType eps = tools::epsilon<RealType>();
    static const struct { double n, x; long double pdf, cdf, ccdf; } data[] = {
        { 10, 0.125, 0.2794000234541079717388096723565994485867L, 0.002361174787509568195989680164454239235797L, 0.9976388252124904318040103198355457607642L },
        { 10, 0.25,  5.193400831298942628933394720192144589244L, 0.4404402898047362231585527709511894227082L, 0.5595597101952637768414472290488105772918L },
        { 10, 0.375, 1.800079163884122092670200607964937329104L, 0.8799166787843709758060153896938758782439L, 0.1200833212156290241939846103061241217561L },
        { 10, 0.5,   0.2695175501788391040012563523532817140662L, 0.9865241101241363106265931874745960750216L, 0.01347588987586368937340681252540392497843L },
        { 10, 0.75,  7.804378592440503887211784902487782139377e-4L, 0.9999739854046918648152915125506947553913L, 2.601459530813518470848744930524460867172e-5L },
        { 10, 1.0,   1.648922897950846262372752246369412370047e-7L, 0.9999999958776927551228843440681192757854L, 4.122307244877115655931880724214614195843e-9L },
        { 1,  0.25,  4.128551064736226095083086707053011290815e-6L, 2.682381008482982753810582973667256706209e-8L, 0.9999999731761899151701724618941702633274L },
        { 1,  0.5,   0.6395828509404566346454593384123317079184L, 0.03605475633512490561408610371791087180628L, 0.9639452436648750943859138962820891281937L },
        { 1,  0.75,  1.683460951304975381694233065698018862373L, 0.3728329582237383585063344713966983906617L, 0.6271670417762616414936655286033016093383L },
        { 1,  1.5,   0.1333072274198802100373204656125755644301L, 0.9777820373834748712794563853895431189662L, 0.02221796261652512872054361461045688103382L },
        { 1,  2.5,   7.453306344157340442849824580766625946801e-5L, 0.9999925466936558426583999002666408827042L, 7.453306344157341600099733359117295757206e-6L },
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
