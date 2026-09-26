// Copyright Philipp C. J. Muenster 2020.
// Copyright Paul A. Bristow 2010.
// Copyright John Maddock 2007.
// Copyright Matt Borland 2022.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/distributions/von_mises.hpp>
#include <boost/math/quadrature/gauss.hpp>
#ifndef BOOST_MATH_STANDALONE
#include <boost/multiprecision/cpp_bin_float.hpp>
#endif
#include "math_unit_test.hpp"

#include <cmath>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <type_traits>

using boost::math::von_mises_distribution;

// strtold rather than lexical_cast, which rejects values that underflow T.
template <class T>
typename std::enable_if<std::is_floating_point<T>::value, T>::type from_string(const char* s)
{
    return static_cast<T>(std::strtold(s, nullptr));
}

template <class T>
typename std::enable_if<!std::is_floating_point<T>::value, T>::type from_string(const char* s)
{
    return T(s);
}

// The relative error of exp(-y) is |y| times that of y, so allow for it.
template <class T>
int exp_tolerance(T k, T u)
{
    using std::sin;
    T s = sin(u / 2);
    return 32 + static_cast<int>(8 * k * s * s);
}

// Beyond double the cdf uses tanh-sinh quadrature, which takes milliseconds, so those types get sparser loops.
template <class T>
bool is_fast()
{
    return std::is_floating_point<T>::value;
}

template <class T>
T grid_step()
{
    return is_fast<T>() ? T(0.125) : T(1);
}

template <class T>
void test_spots()
{
    // {kappa, u, pdf, cdf, complement of cdf}, from mpmath at 60 digits.
    static const char* const data[][5] = {
      {"0.5", "-3.0", "0.0912252976461840379820933577361329255606", "0.01287384404011132628315119923691059344663", "0.9871261559598886737168488007630894065534"},
      {"0.5", "-1.0", "0.196071550526524086413763451689326250825", "0.2715226591251793199448800702551483284961", "0.7284773408748206800551199297448516715039"},
      {"0.5", "-0.25", "0.2429327617033683559013882330159184120864", "0.4386341927772687107551605685702497556036", "0.5613658072227312892448394314297502443964"},
      {"0.5", "0.0", "0.2467383573941201474772991642849040692603", "0.5", "0.5"},
      {"0.5", "0.5", "0.232088734383606427838415056057829590115", "0.6208770268159540434436100108945377629265", "0.3791229731840459565563899891054622370735"},
      {"0.5", "2.0", "0.1215414157555752772890403705423004325392", "0.8848231100868918132924501451608733377203", "0.1151768899131081867075498548391266622797"},
      {"5.0", "-3.0", "0.00004138792752915659777616252209971139913074", "0.000005668659921619703012317037101563918361689", "0.9999943313400783802969876829628984360816"},
      {"5.0", "-1.0", "0.08706961456376827266327441069660051952476", "0.0190379545318531094541861520374823651998", "0.9809620454681468905458138479625176348002"},
      {"5.0", "-0.25", "0.7423037643464445689948461977563146163048", "0.2939650204800054031218938580632133253127", "0.7060349795199945968781061419367866746873"},
      {"5.0", "0.0", "0.8671365285423520025735184696977704534391", "0.5", "0.5"},
      {"5.0", "0.5", "0.4701770125291431414758802335102456690076", "0.8586624818092983238327382084308137380641", "0.1413375181907016761672617915691862619359"},
      {"5.0", "2.0", "0.0007293965387295747468161420565478164792032", "0.999810518550356705764482062195139854066", "0.0001894814496432942355179378048601459340429"},
      {"50.0", "-3.0", "1.726474728701131897274231580755440210442e-43", "1.771559480845494822618388772926421989873e-44", "1.0"},
      {"50.0", "-1.0", "0.0000000002931498292824351897359431582750173746876", "0.000000000006869235154586915768865963174674043998617", "0.9999999999931307648454130842311340368253"},
      {"50.0", "-0.25", "0.5946207433807447509051628129560149739889", "0.03931556815087087478901111487392325730048", "0.9606844318491291252109888851260767426995"},
      {"50.0", "0.0", "2.81383249608255054589931212712591310859", "0.5", "0.5"},
      {"50.0", "0.5", "0.006180695497707893157202489669503380565833", "0.9997582277598794729826240424180070740493", "0.0002417722401205270173759575819929259507402"},
      {"50.0", "2.0", "4.989538832605930159529018130048728362805e-31", "0.9999999999999999999999999999999889056431", "1.109435685709275223350169131875357208308e-32"},
      {"1000.0", "-3.0", "7.212659818487654971016715045664679601909e-864", "5.419758304748379754068015875706168003294e-866", "1.0"},
      {"1000.0", "-1.0", "2.862093693516050071781425702949356969262e-199", "3.39871343372991781227724531225219529813e-202", "1.0"},
      {"1000.0", "-0.25", "3.978249758044480108064723961370962347826e-13", "1.583689032406510412742370845548162517928e-15", "0.9999999999999984163109675934895872576292"},
      {"1000.0", "0.0", "12.61408496162744720188900353346495083738", "0.5", "0.5"},
      {"1000.0", "0.5", "8.622593097199372341764522813836560181269e-53", "1.0", "1.791743860692674111721063817439501913897e-55"},
      {"1000.0", "2.0", "1.191513794511670949059178707139141346055e-614", "1.0", "1.311029608541689282749467334862523737376e-617"},
      {"1000000.0", "-0.25", "2.736548875920210235984480318011804339029e-13499", "1.10608802206099512319217301932265410593e-13504", "1.0"},
      {"1000000.0", "0.0", "398.9422305336258105819158952373769957743", "0.5", "0.5"},
    };
    for (auto const& row : data)
    {
        const T k = from_string<T>(row[0]);
        const T u = from_string<T>(row[1]);
        const T expected_pdf = from_string<T>(row[2]);
        const T expected_cdf = from_string<T>(row[3]);
        const T expected_ccdf = from_string<T>(row[4]);
        const int tol = exp_tolerance(k, u);
        // The distribution is a translate, so a dyadic mean changes nothing.
        for (T mean : {T(0), T(0.75), T(-2)})
        {
            if (!is_fast<T>() && (mean != 0))
            {
                continue;
            }
            von_mises_distribution<T> dist(mean, k);
            const T x = mean + u;
            if (expected_pdf >= (std::numeric_limits<T>::min)())
            {
                CHECK_ULP_CLOSE(expected_pdf, pdf(dist, x), tol);
            }
            if (expected_cdf >= (std::numeric_limits<T>::min)())
            {
                CHECK_ULP_CLOSE(expected_cdf, cdf(dist, x), tol);
            }
            if (expected_ccdf >= (std::numeric_limits<T>::min)())
            {
                CHECK_ULP_CLOSE(expected_ccdf, cdf(complement(dist, x)), tol);
            }
        }
    }
}

template <class T>
void test_moments()
{
    // {kappa, circular variance, circular standard deviation, entropy}, from mpmath at 60 digits.
    static const char* const data[][4] = {
      {"0.5", "0.7575003874191980546492976464963645925877", "1.683303398802970259379764397122457347388", "1.77817697930442581482659287057323642698"},
      {"5.0", "0.1066168629559147784129949927750508272456", "0.4748468074586960312096852738531094009897", "0.6756431570114528094714655330324416977434"},
      {"50.0", "0.0100510326215022474073440705353173833048", "0.1421399682332739278516634966880758092353", "-0.531995800643737561909134538702056029257"},
      {"1000.0", "0.0005001251251957198010818256520440214055299", "0.03163068856284189459460209894809633476081", "-2.034688918525470043889750394297228431701"},
      {"1000000.0", "0.0000005000001250001250001953129062510478547813", "0.001000000250000197916929687995129242568135", "-5.488816495777276810013227453167050890645"},
    };
    for (auto const& row : data)
    {
        von_mises_distribution<T> dist(1, from_string<T>(row[0]));
        CHECK_ULP_CLOSE(from_string<T>(row[1]), variance(dist), 8);
        CHECK_ULP_CLOSE(from_string<T>(row[2]), standard_deviation(dist), 8);
        CHECK_ULP_CLOSE(from_string<T>(row[3]), entropy(dist), 8);
        CHECK_EQUAL(mean(dist), T(1));
        CHECK_EQUAL(mode(dist), T(1));
        CHECK_EQUAL(median(dist), T(1));
        CHECK_EQUAL(skewness(dist), T(0));
    }
}

// With zero concentration the distribution is uniform on the circle.
template <class T>
void test_uniform()
{
    const T pi = boost::math::constants::pi<T>();
    von_mises_distribution<T> dist(0, 0);
    for (T u = -3; u <= 3; u += T(0.25))
    {
        CHECK_ULP_CLOSE(1 / (2 * pi), pdf(dist, u), 2);
        CHECK_ULP_CLOSE((u + pi) / (2 * pi), cdf(dist, u), 4);
        CHECK_ULP_CLOSE((pi - u) / (2 * pi), cdf(complement(dist, u)), 4);
    }
    CHECK_ULP_CLOSE(T(1), variance(dist), 1);
    using std::log;
    CHECK_ULP_CLOSE(log(2 * pi), entropy(dist), 2);
}

template <class T>
void test_quantile()
{
    using std::abs;
    const T eps = std::numeric_limits<T>::epsilon();
    for (T k : {T(0), T(0.5), T(5), T(50), T(1000)})
    {
        if (!is_fast<T>() && (k != T(0.5)) && (k != T(50)))
        {
            continue;
        }
        von_mises_distribution<T> dist(T(0.5), k);
        for (T u = -3; u <= 3; u += grid_step<T>())
        {
            const T x = T(0.5) + u;
            // Invert whichever of the cdf and its complement is the small one, as that carries u exactly.
            if (u <= 0)
            {
                const T p = cdf(dist, x);
                if (p >= (std::numeric_limits<T>::min)())
                {
                    CHECK_ABSOLUTE_ERROR(x, quantile(dist, p), T(16 * eps * (1 + abs(x))));
                }
            }
            else
            {
                const T q = cdf(complement(dist, x));
                if (q >= (std::numeric_limits<T>::min)())
                {
                    CHECK_ABSOLUTE_ERROR(x, quantile(complement(dist, q)), T(16 * eps * (1 + abs(x))));
                }
            }
        }
        CHECK_EQUAL(quantile(dist, T(0.5)), T(0.5));
        CHECK_EQUAL(quantile(complement(dist, T(0.5))), T(0.5));
        CHECK_ULP_CLOSE(support(dist).first, quantile(dist, T(0)), 1);
        CHECK_ULP_CLOSE(support(dist).second, quantile(dist, T(1)), 1);
        CHECK_ULP_CLOSE(support(dist).second, quantile(complement(dist, T(0))), 1);
        CHECK_ULP_CLOSE(support(dist).first, quantile(complement(dist, T(1))), 1);
    }
}

// The cdf is the integral of the pdf, and the pdf integrates to one.
template <class T>
void test_pdf_integral()
{
    const T pi = boost::math::constants::pi<T>();
    for (T k : {T(0.25), T(3), T(30)})
    {
        von_mises_distribution<T> dist(0, k);
        auto f = [&](T x) { return pdf(dist, x); };
        // Composite 30-point Gauss-Legendre: the pdf is entire, so this is accurate to rounding.
        auto integrate = [&](T a, T b)
        {
            T sum = 0;
            for (int i = 0; i < 8; ++i)
            {
                sum += boost::math::quadrature::gauss<T, 30>::integrate(f, a + (b - a) * i / 8, a + (b - a) * (i + 1) / 8);
            }
            return sum;
        };
        CHECK_ULP_CLOSE(T(1), integrate(-pi, pi), 16);
        for (T u : {T(-2.5), T(-1), T(-0.125), T(0.5), T(2)})
        {
            CHECK_ULP_CLOSE(integrate(-pi, u), cdf(dist, u), 64);
        }
    }
}

template <class T>
void test_symmetry()
{
    for (T mean : {T(0), T(1.5), T(-3)})
    {
        for (T k : {T(0.5), T(4), T(40)})
        {
            if (!is_fast<T>() && ((mean != T(1.5)) || (k != 4)))
            {
                continue;
            }
            von_mises_distribution<T> dist(mean, k);
            for (T x = T(1) / 16; x < 3; x += grid_step<T>() / 2)
            {
                CHECK_EQUAL(pdf(dist, mean + x), pdf(dist, mean - x));
                CHECK_EQUAL(cdf(dist, mean - x), cdf(complement(dist, mean + x)));
                CHECK_EQUAL(cdf(dist, mean + x), cdf(complement(dist, mean - x)));
            }
        }
    }
}

template <class T>
void test_support()
{
    const T pi = boost::math::constants::pi<T>();
    von_mises_distribution<T> dist(1, 2);
    CHECK_ULP_CLOSE(1 - pi, support(dist).first, 1);
    CHECK_ULP_CLOSE(1 + pi, support(dist).second, 1);
    CHECK_EQUAL(pdf(dist, T(-3)), T(0));
    CHECK_EQUAL(pdf(dist, T(5)), T(0));
    CHECK_EQUAL(cdf(dist, T(-3)), T(0));
    CHECK_EQUAL(cdf(dist, T(5)), T(1));
    CHECK_EQUAL(cdf(complement(dist, T(-3))), T(1));
    CHECK_EQUAL(cdf(complement(dist, T(5))), T(0));
}

template <class T>
void test_errors()
{
    const T inf = std::numeric_limits<T>::infinity();
    const T nan = std::numeric_limits<T>::quiet_NaN();
    CHECK_THROW(von_mises_distribution<T>(0, -1), std::domain_error);
    CHECK_THROW(von_mises_distribution<T>(0, inf), std::domain_error);
    CHECK_THROW(von_mises_distribution<T>(0, nan), std::domain_error);
    CHECK_THROW(von_mises_distribution<T>(inf, 1), std::domain_error);
    CHECK_THROW(von_mises_distribution<T>(nan, 1), std::domain_error);
    von_mises_distribution<T> dist(0, 1);
    CHECK_THROW(pdf(dist, nan), std::domain_error);
    CHECK_THROW(cdf(dist, nan), std::domain_error);
    CHECK_THROW(cdf(complement(dist, nan)), std::domain_error);
    CHECK_THROW(quantile(dist, T(-1)), std::domain_error);
    CHECK_THROW(quantile(dist, T(2)), std::domain_error);
    CHECK_THROW(quantile(complement(dist, nan)), std::domain_error);
}

template <class T>
void test_all()
{
    test_spots<T>();
    test_moments<T>();
    test_uniform<T>();
    test_quantile<T>();
    if (is_fast<T>())
    {
        test_pdf_integral<T>();
    }
    test_symmetry<T>();
    test_support<T>();
    test_errors<T>();
}

int main()
{
    boost::math::von_mises dist;
    CHECK_EQUAL(dist.mean(), dist.location());
    CHECK_EQUAL(dist.concentration(), dist.scale());
    CHECK_EQUAL(dist.concentration(), 1.0);

    test_all<float>();
    test_all<double>();
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
    test_all<long double>();
#endif
#ifndef BOOST_MATH_STANDALONE
    test_all<boost::multiprecision::cpp_bin_float_quad>();
#endif
    return boost::math::test::report_errors();
}
