// Copyright Paul A. Bristow 2010.
// Copyright John Maddock 2010.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifdef _MSC_VER
#  pragma warning (disable : 4224) // nonstandard extension used : formal parameter 'type' was previously defined as a type
// in Boost.test and lexical_cast
#  pragma warning (disable : 4310) // cast truncates constant value
#  pragma warning (disable : 4512) // assignment operator could not be generated

#endif

//#include <pch.hpp> // include directory libs/math/src/tr1/ is needed.

#include <boost/math/tools/config.hpp>
#include "../include_private/boost/math/tools/test.hpp"

#ifndef BOOST_MATH_NO_REAL_CONCEPT_TESTS
#include <boost/math/concepts/real_concept.hpp> // for real_concept
#endif

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp> // Boost.Test
#include <boost/test/tools/floating_point_comparison.hpp>

#include <boost/math/distributions/inverse_gaussian.hpp>
#include <boost/math/special_functions/next.hpp>
#ifdef BOOST_MATH_TEST_MULTIPRECISION
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#endif
using boost::math::inverse_gaussian_distribution;
using boost::math::inverse_gaussian;

#include "test_out_of_range.hpp"

#include <iostream>
#include <iomanip>
using std::cout;
using std::endl;
using std::setprecision;
#include <limits>
#include <cfloat>
using std::numeric_limits;
#include <cmath>
using std::log;

template <class RealType>
void check_inverse_gaussian(RealType mean, RealType scale, RealType x, RealType p, RealType q, RealType tol)
{
 using boost::math::inverse_gaussian_distribution;

  BOOST_CHECK_CLOSE_FRACTION(
    ::boost::math::cdf(   // Check cdf
    inverse_gaussian_distribution<RealType>(mean, scale),      // distribution.
    x),    // random variable.
    p,     // probability.
    tol);   // tolerance.
  BOOST_CHECK_CLOSE_FRACTION(
    ::boost::math::cdf( // Check cdf complement
    complement( 
    inverse_gaussian_distribution<RealType>(mean, scale),   // distribution.
    x)),   // random variable.
    q,      // probability complement.
    tol);    // %tolerance.
  BOOST_CHECK_CLOSE_FRACTION(
    ::boost::math::quantile( // Check quantile
    inverse_gaussian_distribution<RealType>(mean, scale),    // distribution.
    p),   // probability.
    x,   // random variable.
    tol);   // tolerance.
  BOOST_CHECK_CLOSE_FRACTION(
    ::boost::math::quantile( // Check quantile complement
    complement(
    inverse_gaussian_distribution<RealType>(mean, scale),   // distribution.
    q)),   // probability complement.
    x,     // random variable.
    tol);  // tolerance.

   inverse_gaussian_distribution<RealType> dist (mean, scale);

   if((p < 0.999) && (q < 0.999))
   {  // We can only check this if P is not too close to 1,
      // so that we can guarantee Q is accurate:
      BOOST_CHECK_CLOSE_FRACTION(
        cdf(complement(dist, x)), q, tol); // 1 - cdf
      BOOST_CHECK_CLOSE_FRACTION(
        quantile(dist, p), x, tol); // quantile(cdf) = x
      BOOST_CHECK_CLOSE_FRACTION(
        quantile(complement(dist, q)), x, tol); // quantile(complement(1 - cdf)) = x
   }
}

template <class RealType>
void test_spots(RealType)
{
  // Basic sanity checks
  RealType tolerance = static_cast<RealType>(1e-4L); // 
  cout << "Tolerance for type " << typeid(RealType).name()  << " is " << tolerance << endl;

  // Check some bad parameters to the distribution,
#ifndef BOOST_NO_EXCEPTIONS
  BOOST_MATH_CHECK_THROW(boost::math::inverse_gaussian_distribution<RealType> nbad1(0, 0), std::domain_error); // zero scale
  BOOST_MATH_CHECK_THROW(boost::math::inverse_gaussian_distribution<RealType> nbad1(0, -1), std::domain_error); // negative scale
#else
  BOOST_MATH_CHECK_THROW(boost::math::inverse_gaussian_distribution<RealType>(0, 0), std::domain_error); // zero scale
  BOOST_MATH_CHECK_THROW(boost::math::inverse_gaussian_distribution<RealType>(0, -1), std::domain_error); // negative scale
#endif

  inverse_gaussian_distribution<RealType> w11;

  // Error tests:
  check_out_of_range<inverse_gaussian_distribution<RealType> >(0.25, 1);
  
  // Check complements.

    BOOST_CHECK_CLOSE_FRACTION(
     cdf(complement(w11, 1.)), static_cast<RealType>(1) - cdf(w11, 1.), tolerance); // cdf complement
    // cdf(complement = 1 - cdf  - but if cdf near unity, then loss of accuracy in cdf,
    // but cdf complement is near zero but more accurate.

     BOOST_CHECK_CLOSE_FRACTION( // quantile(complement p) == quantile(1 - p)
     quantile(complement(w11, static_cast<RealType>(0.5))), 
     quantile(w11, 1 - static_cast<RealType>(0.5)),
     tolerance); // cdf complement

  check_inverse_gaussian(
     static_cast<RealType>(2),
     static_cast<RealType>(3),
     static_cast<RealType>(1),
     static_cast<RealType>(0.28738674440477374),
     static_cast<RealType>(1 - 0.28738674440477374),
     tolerance);

  RealType tolfeweps = boost::math::tools::epsilon<RealType>() * 5;

  inverse_gaussian_distribution<RealType> dist(2, 3);

  using namespace std; // ADL of std names.
  // mean:
  BOOST_CHECK_CLOSE_FRACTION(mean(dist),
    static_cast<RealType>(2), tolfeweps);
  BOOST_CHECK_CLOSE_FRACTION(scale(dist),
    static_cast<RealType>(3), tolfeweps);

  // variance:
  BOOST_CHECK_CLOSE_FRACTION(variance(dist),
    static_cast<RealType>(2.6666666666666666666666666666666666666666666666666666666667L), 1000*tolfeweps);
  // std deviation:
  BOOST_CHECK_CLOSE_FRACTION(standard_deviation(dist), 
    static_cast<RealType>(1.632993L), 1000 * tolerance);
  //// hazard:
  //BOOST_CHECK_CLOSE_FRACTION(hazard(dist, x),
  //  pdf(dist, x) / cdf(complement(dist, x)), tolerance);
  //// cumulative hazard:
  //BOOST_CHECK_CLOSE_FRACTION(chf(dist, x),
  //  -log(cdf(complement(dist, x))), tolerance);
  // coefficient_of_variation:
  BOOST_CHECK_CLOSE_FRACTION(coefficient_of_variation(dist),
    standard_deviation(dist) / mean(dist), tolerance);
  // mode:
  BOOST_CHECK_CLOSE_FRACTION(mode(dist),
    static_cast<RealType>(0.8284271L), tolerance);

  // median
  BOOST_CHECK_CLOSE_FRACTION(median(dist),
    static_cast<RealType>(1.5122506636053668L), tolerance);
  // Fails for real_concept - because std::numeric_limits<RealType>::digits = 0

  // skewness:
  BOOST_CHECK_CLOSE_FRACTION(skewness(dist),
    static_cast<RealType>(2.449490L), tolerance);
  // kurtosis:
  BOOST_CHECK_CLOSE_FRACTION(kurtosis(dist),
    static_cast<RealType>(10-3), tolerance);
  BOOST_CHECK_CLOSE_FRACTION(kurtosis_excess(dist),
    static_cast<RealType>(10), tolerance);
} // template <class RealType>void test_spots(RealType)

template <class RealType>
void test_large_shape_cdf(RealType)
{
  // Independent high-precision values from the defining CDF, using Decimal
  // arithmetic, an erf power series, and a bounded erfc asymptotic expansion.
  const long double data[][4] = {
    { 100L, 1L, 0.519897615648327031592108503935940983113411L, 0.480102384351672968407891496064059016886589L },
    { 300L, 1L, 0.511506898482595696440706140009078301888645L, 0.488493101517404303559293859990921698111355L },
    { 500L, 1L, 0.508916166944271025203908756615335656358401L, 0.491083833055728974796091243384664343641599L },
    { 700L, 1L, 0.507536610711149382901383525090917369515136L, 0.492463389288850617098616474909082630484864L },
    { 1000L, 1L, 0.506306255528466690646614735956978724498991L, 0.493693744471533309353385264043021275501009L },
    { 2000L, 1L, 0.504459752960542115964542559424352811401032L, 0.495540247039457884035457440575647188598968L },
    { 5000L, 1L, 0.502820806891494716451778228503475775359353L, 0.497179193108505283548221771496524224640647L },
    { 10000L, 1L, 0.501994661537961729660689505416408573402360L, 0.498005338462038270339310494583591426597640L },
    { 500L, 0.96875L, 0.245795920321378254279271385912219740338036L, 0.754204079678621745720728614087780259661964L },
    { 500L, 1.03125L, 0.761341593959515617582627727530040166229774L, 0.238658406040484382417372272469959833770226L },
    { 1000L, 0.96875L, 0.161492549505796665675678816590357590761845L, 0.838507450494203334324321183409642409238155L },
    { 1000L, 1.03125L, 0.838681329955945237527053545092684829137562L, 0.161318670044054762472946454907315170862438L },
    { 500L, 0.75L, 6.20243988565182005344254309057667117068051e-11L, 0.999999999937975601143481799465574569094233L },
    { 500L, 1.25L, 0.999999746370348504344912460794892721683629L, 2.53629651495655087539205107278316371144419e-7L },
    { 10000L, 0.96875L, 0.000762081454762838844391132060849477208809302L, 0.999237918545237161155608867939150522791191L },
    { 10000L, 1.03125L, 0.998973049207792457428094243426255769746681L, 0.00102695079220754257190575657374423025331869L },
  };
  const RealType epsilon = boost::math::tools::epsilon<RealType>();
  using std::frexp;
  using std::ldexp;
  int max_exponent;
  frexp(boost::math::tools::max_value<RealType>(), &max_exponent);
  const RealType means[] = { RealType(1), RealType(0.00000000023283064365386962890625L),
    RealType(4294967296.0L), 16 * boost::math::tools::min_value<RealType>(),
    ldexp(RealType(1), max_exponent - 16), 32 * std::numeric_limits<RealType>::denorm_min() };
  for (const auto& row : data)
  {
    // Exponential relative sensitivity grows with a*a/2. Keep the central
    // 128-epsilon bound and use a 32-epsilon budget per exponent in the tails.
    const RealType y_squared = static_cast<RealType>(row[0] * (row[1] - 1) * (row[1] - 1) / (2 * row[1]));
    const RealType tolerance = (y_squared > 3 ? 32 * (1 + y_squared) : RealType(128)) * epsilon;
    for (unsigned i = 0; i < sizeof(means) / sizeof(means[0]); ++i)
    {
      const RealType mu = means[i];
      if (mu == 0)
        continue;
      BOOST_TEST_CONTEXT("type=" << typeid(RealType).name() << ", shape=" << row[0]
        << ", x/mean=" << row[1] << ", mean=" << mu)
      {
      const inverse_gaussian_distribution<RealType> dist(mu, mu * static_cast<RealType>(row[0]));
      const RealType x = mu * static_cast<RealType>(row[1]);
      const RealType p = cdf(dist, x);
      const RealType q = cdf(complement(dist, x));
      BOOST_CHECK((boost::math::isfinite)(p));
      BOOST_CHECK((boost::math::isfinite)(q));
      BOOST_CHECK(p >= 0 && p <= 1);
      BOOST_CHECK(q >= 0 && q <= 1);
      BOOST_CHECK_CLOSE_FRACTION(p, static_cast<RealType>(row[2]), tolerance);
      BOOST_CHECK_CLOSE_FRACTION(q, static_cast<RealType>(row[3]), tolerance);
      BOOST_CHECK_CLOSE_FRACTION(p + q, RealType(1), 128 * epsilon);
      if ((i < 3) && (row[0] == 500L) && (row[2] > 0.01L) && (row[3] > 0.01L))
      {
        BOOST_CHECK_CLOSE_FRACTION(quantile(dist, p), x, 128 * epsilon);
        BOOST_CHECK_CLOSE_FRACTION(quantile(complement(dist, q)), x, 128 * epsilon);
      }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE( test_large_shape )
{
  test_large_shape_cdf(0.0F);
  test_large_shape_cdf(0.0);
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
  test_large_shape_cdf(0.0L);
#ifndef BOOST_MATH_NO_REAL_CONCEPT_TESTS
  test_large_shape_cdf(boost::math::concepts::real_concept(0.));
#endif
#endif
}

template <class RealType>
void test_deep_tail_cdf(RealType)
{
  // Independent Decimal references use separately bounded positive erfcx tails.
  const long double data[][3] = {
    { 21L, 0.125L, 7.291459004754112603352941093615415096235e-30L },
    { 21L, 8L, 9.007363557162028734338325443837286235165e-31L },
    { 22L, 0.125L, 3.332796457964578353784931221427005072469e-31L },
    { 22L, 8L, 4.119266157712698228360653991204015847817e-32L },
    { 23L, 0.125L, 1.524901282029210902403448998373492554917e-32L },
    { 23L, 8L, 1.885650541084410056066653349065387755587e-33L },
    { 176L, 0.125L, 1.775905505448122785699924483934182217827e-236L },
    { 176L, 8L, 2.216690372486544511107099941569763090701e-237L },
    { 177L, 0.125L, 8.282559769103606776813053919204554506461e-238L },
    { 177L, 8L, 1.033839877291871420250393937221921397053e-238L },
    { 178L, 0.125L, 3.862924444359401490606536105737107118119e-239L },
    { 178L, 8L, 4.821791147898428401109417195269734536297e-240L },
#if LDBL_MAX_EXP > DBL_MAX_EXP
    { 2838L, 0.125L, 1.301682246950782425490244544260728669449e-3777L },
#endif
#if LDBL_MAX_EXP > DBL_MAX_EXP
    { 2838L, 8L, 1.626957235152347072563874816857493004111e-3778L },
#endif
#if LDBL_MAX_EXP > DBL_MAX_EXP
    { 2839L, 0.125L, 6.086976674034536842513982652262554314508e-3779L },
#endif
#if LDBL_MAX_EXP > DBL_MAX_EXP
    { 2839L, 8L, 7.608040345690346295911004284599736386307e-3780L },
#endif
#if LDBL_MAX_EXP > DBL_MAX_EXP
    { 2840L, 0.125L, 2.846415660587392915323404759687483152324e-3780L },
#endif
#if LDBL_MAX_EXP > DBL_MAX_EXP
    { 2840L, 8L, 3.557701471176804781448244810749566553326e-3781L },
#endif
    { 8.59258937835693359375L, 0.125L, 3.591174609673238026368113606805827534608e-13L },
    { 8.5925922393798828125L, 0.125L, 3.591142563890887388309235085918784231035e-13L },
    { 21.74999237060546875L, 1L, 0.4577044632823731484311723032884239273365L },
    { 21.75L, 1L, 0.4577044705390353035423688153929315176905L },
    { 19.333332061767578125L, 2L, 0.0006015374450221001062105790455698042612565L },
    { 19.3333339691162109375L, 2L, 0.0006015371347789923742069316649113152221024L },
    { 8.59259033203125L, 8L, 4.365322672757212439356936307049058419658e-14L },
    { 8.5925922393798828125L, 8L, 4.365296729096854790521451435110168247173e-14L },
    { 127.5390777587890625L, 2L, 4.630597763569726123318880453345364749019e-16L },
    { 127.53908538818359375L, 2L, 4.630588798766307767824499779348899387908e-16L },
    { 10.4113521575927734375L, 8L, 1.522152219372173594557917600251344281880e-16L },
    { 10.41135311126708984375L, 8L, 1.522147708819046419858753962289249089336e-16L },
    { 69.92592592592590960975940106436610221863L, 0.125L, 3.394078336304087657704167068057799217965e-95L },
    { 69.92592592592592382061411626636981964111L, 0.125L, 3.394078336303939601328638456009704213323e-95L },
    { 177L, 1L, 0.4850279186644429965584006013863277807435L },
    { 177.0000000000000284217094304040074348450L, 1L, 0.4850279186644429977570968689302612237533L },
    { 157.3333333333333143855270463973283767700L, 2L, 2.437400168809340833531564466307626617317e-19L },
    { 157.3333333333333428072364768013358116150L, 2L, 2.437400168809323302244586693893501812231e-19L },
    { 69.92592592592590960975940106436610221863L, 8L, 4.227330446778904380408703083920300612814e-96L },
    { 69.92592592592593803146883146837353706360L, 8L, 4.227330446778535577613184648394072842993e-96L },
    { 288.3492271129371715687739197164773941040L, 2L, 1.081367401572126457834984206768181645347e-33L },
    { 288.3492271129372284121927805244922637939L, 2L, 1.081367401572110986161455678712190078912e-33L },
    { 23.53871241738262654052959987893700599670L, 8L, 3.581735565814747988118739514752194335691e-34L },
    { 23.53871241738263009324327867943793535233L, 8L, 3.581735565814708756472372152945893559265e-34L },
#if LDBL_MAX_EXP > DBL_MAX_EXP
    { 1121.481481481481481177198133991623762995L, 0.125L, 2.146856822417751309921693777797030596742e-1494L },
#endif
#if LDBL_MAX_EXP > DBL_MAX_EXP
    { 1121.481481481481481510265041379170725122L, 0.125L, 2.146856822417749119771652306994510402888e-1494L },
#endif
    { 2838.749999999999999777955395074968691915L, 1L, 0.4962564963752824700551812971624815459326L },
    { 2838.75L, 1L, 0.4962564963752824700553276782176961266827L },
    { 2523.333333333333333037273860099958255887L, 2L, 8.061475687671561576395551403009632536180e-277L },
    { 2523.333333333333333481363069950020872056L, 2L, 8.061475687671560680684195544090990589351e-277L },
#if LDBL_MAX_EXP > DBL_MAX_EXP
    { 1121.481481481481481399242738916655071080L, 8L, 2.682963664854186619255662358356177570044e-1495L },
#endif
#if LDBL_MAX_EXP > DBL_MAX_EXP
    { 1121.481481481481481621287343841686379164L, 8L, 2.682963664854184794543822096148391707832e-1495L },
#endif
    { 349.3461790022124359156308059937146026641L, 2L, 2.346294912654197788941015666145882964918e-40L },
    { 349.3461790022124359433863816093435161747L, 2L, 2.346294912654197772568618770075048031901e-40L },
    { 28.51805542875203558668417702648412159760L, 8L, 7.779918162224078419313680298233549591858e-41L },
    { 28.51805542875203559188834745441454288084L, 8L, 7.779918162224078294628012728206757545041e-41L },
  };
  const RealType epsilon = boost::math::tools::epsilon<RealType>();
  for (const auto& row : data)
  {
    const inverse_gaussian_distribution<RealType> dist(1, static_cast<RealType>(row[0]));
    const RealType x = static_cast<RealType>(row[1]);
    const RealType p = cdf(dist, x);
    const RealType q = cdf(complement(dist, x));
    const RealType small = row[1] < 1 ? p : q;
    const RealType expected = static_cast<RealType>(row[2]);
    const RealType y_squared = static_cast<RealType>(row[0] * (row[1] - 1) * (row[1] - 1) / (2 * row[1]));
    // The relative error budget follows the sensitivity of exp(-a*a/2).
    const RealType tolerance = 32 * epsilon * (1 + y_squared);
    BOOST_TEST_CONTEXT("type=" << typeid(RealType).name() << ", shape=" << row[0] << ", x/mean=" << row[1])
    {
      BOOST_CHECK((boost::math::isfinite)(p));
      BOOST_CHECK((boost::math::isfinite)(q));
      BOOST_CHECK(p >= 0 && p <= 1);
      BOOST_CHECK(q >= 0 && q <= 1);
      if (expected == 0)
        BOOST_CHECK_EQUAL(small, RealType(0));
      else
        BOOST_CHECK_CLOSE_FRACTION(small, expected, tolerance);
      BOOST_CHECK_CLOSE_FRACTION(p + q, RealType(1), 128 * epsilon);
    }
  }
}

BOOST_AUTO_TEST_CASE( test_deep_tail )
{
  test_deep_tail_cdf(0.0F);
  test_deep_tail_cdf(0.0);
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
  test_deep_tail_cdf(0.0L);
#ifndef BOOST_MATH_NO_REAL_CONCEPT_TESTS
  test_deep_tail_cdf(boost::math::concepts::real_concept(0.));
#endif
#endif
}

BOOST_AUTO_TEST_CASE( test_deep_tail_small_shape )
{
  // Both positive erfc tails were evaluated independently at 700 decimal digits.
  const long double data[][3] = {
    { 9.999999999999999451532714542095716517295e-21L, 4e22L, 1.370012494729595008319815564602529481291e-111L },
    { 1.000000000000000019991899802602883619648e-100L, 3.999999999999999908198053060981346513787e+102L, 1.370012494729580787248448498687530730617e-191L },
    { 9.999999999999999821002623990827595960544e-201L, 3.999999999999999606989836527856692108829e+202L, 1.370012494729611982162726656501876826736e-291L },
  };
  for (const auto& row : data)
  {
    const inverse_gaussian dist(1, static_cast<double>(row[0]));
    const double x = static_cast<double>(row[1]);
    const double q = cdf(complement(dist, x));
    const double tolerance = 32 * boost::math::tools::epsilon<double>() * (1 + dist.scale() * x / 2);
    BOOST_CHECK((boost::math::isfinite)(q));
    BOOST_CHECK(q > 0);
    BOOST_CHECK_CLOSE_FRACTION(q, static_cast<double>(row[2]), tolerance);
    BOOST_CHECK_EQUAL(cdf(dist, x), 1);
  }
  typedef boost::math::policies::policy<boost::math::policies::max_series_iterations<1> > limited_policy;
  const inverse_gaussian_distribution<double, limited_policy> limited(1, 1e-20);
  BOOST_MATH_CHECK_THROW(cdf(complement(limited, 4e22)), boost::math::evaluation_error);
}

#ifdef BOOST_MATH_TEST_MULTIPRECISION
template <class RealType>
void check_multiprecision_cdf(const char* scale, const char* x, const char* expected, bool central = false)
{
  const inverse_gaussian_distribution<RealType> dist(RealType(1), RealType(scale));
  const RealType value(x);
  const RealType q = cdf(complement(dist, value));
  const RealType p = cdf(dist, value);
  const RealType ratio = q / RealType(expected);
  const RealType sum = p + q;
  const RealType central_tolerance = 128 * boost::math::tools::epsilon<RealType>();
  const RealType tolerance = central ? central_tolerance
    : 32 * boost::math::tools::epsilon<RealType>() * (1 + dist.scale() * value / 2);
  BOOST_TEST_CONTEXT("type=" << typeid(RealType).name() << ", scale=" << scale << ", x=" << x)
  {
    BOOST_CHECK((boost::math::isfinite)(p));
    BOOST_CHECK((boost::math::isfinite)(q));
    BOOST_CHECK(p >= 0 && p <= 1);
    BOOST_CHECK(q > 0 && q <= 1);
    // Comparing the ratio avoids underflow of the absolute error for types
    // whose exponent range is narrow relative to their significand precision.
    BOOST_CHECK_CLOSE_FRACTION(ratio, RealType(1), tolerance);
    BOOST_CHECK_CLOSE_FRACTION(sum, RealType(1), central_tolerance);
  }
}

BOOST_AUTO_TEST_CASE( test_small_shape_multiprecision )
{
  // Independent Decimal references from separately evaluated erfcx tails.
  const char* data[][3] = {
    { "1e-20", "4e22", "1.37001249472957994316271223215432462122454525140999551934197360681583347626422486125766894044479074637497857752133303e-111" },
    { "1e-100", "4e102", "1.37001249472957994314901210720702882179298680012311892480779017277481694439548691331310796824655220891320662224694514e-191" },
    { "1e-200", "4e202", "1.37001249472957994314901210720702882179298680012311892480779017277481694439548691331310796824655220877620537277398715e-291" },
  };
  for (const auto& row : data)
  {
    check_multiprecision_cdf<boost::multiprecision::cpp_dec_float_50>(row[0], row[1], row[2]);
    check_multiprecision_cdf<boost::multiprecision::cpp_dec_float_100>(row[0], row[1], row[2]);
    check_multiprecision_cdf<boost::multiprecision::cpp_bin_float_50>(row[0], row[1], row[2]);
  }
  typedef boost::multiprecision::number<boost::multiprecision::backends::cpp_bin_float<
    200, boost::multiprecision::backends::digit_base_2, void, int, -500, 500> > narrow_type;
  check_multiprecision_cdf<narrow_type>("1e-40", "4e42",
    "1.37001249472957994314901210720702882179312380137259188280210507398553764727766621882614176076966461312742370135403511e-131");
}

BOOST_AUTO_TEST_CASE( test_survival_series_multiprecision )
{
  // Independent Decimal references; u=3/2 enters recurrence at the first term.
  const char* multi_step_q = "0.04681207925721164097025559101417431566314958201996002592494513184620948224199524719621507929723713261278315099458727972090030446464";
  check_multiprecision_cdf<boost::multiprecision::cpp_dec_float_50>("1", "3", multi_step_q, true);
  check_multiprecision_cdf<boost::multiprecision::cpp_dec_float_100>("1", "3", multi_step_q, true);
  check_multiprecision_cdf<boost::multiprecision::cpp_bin_float_50>("1", "3", multi_step_q, true);

  typedef boost::multiprecision::number<boost::multiprecision::backends::cpp_bin_float<
    200, boost::multiprecision::backends::digit_base_2, void, int, -500, 500> > narrow_type;
  // Adjacent inputs crossing into the survival series at u=1, v=1/2.
  const char* corner[][3] = {
    { "1.414213562373095048801688724209698078569671875376948073176679622894547",
      "1.414213562373095048801688724209698078569671875376948073176678378291491",
      "0.2059605193715049402458665802326639341156629421427760809218245804609142450312212476827828034864383544431633973048825932372398297023" },
    { "1.414213562373095048801688724209698078569671875376948073176679622894547",
      "1.414213562373095048801688724209698078569671875376948073176679622894547",
      "0.2059605193715049402458665802326639341156629421427760809218242582284041260201878244396353187264626219882888244632805755587666907003" },
  };
  for (const auto& row : corner)
    check_multiprecision_cdf<narrow_type>(row[0], row[1], row[2], true);

  // Actual computed u crosses 2, 3 and 4 at these adjacent input pairs.
  const char* switches[][4] = {
    { "4.000000000000000000000000000000000000000000000000000000000004978412222",
      "0.02092363582111373141960376258709922521293656403561454031010929074798552781677954803056498954899168496393156412863584081561172721662",
      "4.000000000000000000000000000000000000000000000000000000000009956824445",
      "0.02092363582111373141960376258709922521293656403561454031010921014898786974815832729207062819185999206106395957955980385682254765521" },
    { "6",
      "0.004849882133702179712128113117466989887821635107637360357146681212152351539518460682022575155907088589143706740831147798465098445979",
      "6.000000000000000000000000000000000000000000000000000000000004978412222",
      "0.004849882133702179712128113117466989887821635107637360357146664385649098428104497840520520195873700211107482062515460028804533195565" },
    { "8",
      "0.001260116932497084584454773499297023095525081445059968066573156084981839328264112072637712608048683745075057908405117320026849427536",
      "8.000000000000000000000000000000000000000000000000000000000009956824445",
      "0.001260116932497084584454773499297023095525081445059968066573147874491690268630198966925169182285601294443961667733174494523088863216" },
  };
  const auto computed_u = [](const narrow_type& value) -> narrow_type
  {
    using std::sqrt;
    const narrow_type w = narrow_type(1) / sqrt(value);
    const narrow_type root_u = 1 / (boost::math::constants::root_two<narrow_type>() * w);
    return root_u * root_u;
  };
  unsigned threshold = 2;
  for (const auto& row : switches)
  {
    BOOST_TEST_CONTEXT("recurrence threshold=" << threshold)
    {
      BOOST_CHECK(computed_u(narrow_type(row[0])) <= threshold);
      BOOST_CHECK(computed_u(narrow_type(row[2])) > threshold);
      check_multiprecision_cdf<narrow_type>("1", row[0], row[1], true);
      check_multiprecision_cdf<narrow_type>("1", row[2], row[3], true);
    }
    ++threshold;
  }
}
#endif

BOOST_AUTO_TEST_CASE( test_survival_series_policy_and_subnormal )
{
  typedef boost::math::policies::policy<boost::math::policies::max_series_iterations<1> > limited_policy;
  const inverse_gaussian_distribution<double, limited_policy> limited(1, 2);
  BOOST_MATH_CHECK_THROW(cdf(complement(limited, 2.0)), boost::math::evaluation_error);
  const inverse_gaussian_distribution<double, limited_policy> one_term(1, 2 * std::sqrt(1.5e-8));
  // One additional term is sufficient; using the last allowed term must succeed.
  BOOST_CHECK_CLOSE_FRACTION(cdf(complement(one_term, std::sqrt(1.5e8))),
    4.78315592525023081717261419515075473466e-6, 128 * boost::math::tools::epsilon<double>());

  if (std::numeric_limits<double>::has_denorm != std::denorm_absent)
  {
    const double scale = (std::numeric_limits<double>::min)();
    const double x = 3 / scale;
    const inverse_gaussian_distribution<double, limited_policy> dist(1, scale);
    // v underflows, but the leading positive survival term remains representable.
    const double q = cdf(complement(dist, x));
    const double expected = scale * 0.019522369872295782685894384688502128270598616598695L;
    const double error_in_ulps = std::abs(q - expected) / std::numeric_limits<double>::denorm_min();
    BOOST_CHECK(q > 0);
    BOOST_CHECK(error_in_ulps <= 2);
    BOOST_CHECK_EQUAL(cdf(dist, x), 1);
  }
}

BOOST_AUTO_TEST_CASE( test_deep_tail_last_iteration_policy )
{
  typedef boost::math::policies::policy<boost::math::policies::max_series_iterations<1> > one_iteration_policy;
  const inverse_gaussian_distribution<double, one_iteration_policy> saturated(1, 1);
  double p = -1;
  double q = -1;
  BOOST_CHECK_NO_THROW(p = cdf(saturated, 1e20));
  BOOST_CHECK_NO_THROW(q = cdf(complement(saturated, 1e20)));
  BOOST_CHECK_EQUAL(p, 1);
  BOOST_CHECK_EQUAL(q, 0);

  typedef boost::math::policies::policy<boost::math::policies::max_series_iterations<10> > final_iteration_policy;
  const inverse_gaussian_distribution<double, final_iteration_policy> finite(1, 1);
  const double x = 400;
  const double y_squared = (x - 1) * (x - 1) / (2 * x);
  const double epsilon = boost::math::tools::epsilon<double>();
  // Independent Decimal reference; the tenth term satisfies the remainder bound.
  const double expected_q = 3.719450726804640441209187365648001756740571742793735086412992904402574891749769661698100891154651936937197258478273770549830930048e-91;
  p = q = -1;
  BOOST_CHECK_NO_THROW(p = cdf(finite, x));
  BOOST_CHECK_NO_THROW(q = cdf(complement(finite, x)));
  BOOST_CHECK_EQUAL(p, 1);
  BOOST_CHECK((boost::math::isfinite)(q));
  BOOST_CHECK(q > 0);
  BOOST_CHECK_CLOSE_FRACTION(q, expected_q, 32 * epsilon * (1 + y_squared));
  BOOST_CHECK_CLOSE_FRACTION(p + q, 1, 128 * epsilon);

  typedef boost::math::policies::policy<boost::math::policies::max_series_iterations<9> > exhausted_policy;
  const inverse_gaussian_distribution<double, exhausted_policy> exhausted(1, 1);
  BOOST_MATH_CHECK_THROW(cdf(exhausted, x), boost::math::evaluation_error);
  BOOST_MATH_CHECK_THROW(cdf(complement(exhausted, x)), boost::math::evaluation_error);
}

template <class RealType>
void test_finite_extremes_cdf(RealType)
{
  const RealType maximum = boost::math::tools::max_value<RealType>();
  const RealType minimum = boost::math::tools::min_value<RealType>();
  const RealType epsilon = boost::math::tools::epsilon<RealType>();
  for (const RealType& mu : { RealType(1), minimum })
  {
    const inverse_gaussian_distribution<RealType> dist(mu, maximum);
    const RealType lower = boost::math::float_prior(mu);
    const RealType upper = boost::math::float_next(mu);
    BOOST_CHECK_EQUAL(cdf(dist, lower), RealType(0));
    BOOST_CHECK_EQUAL(cdf(complement(dist, lower)), RealType(1));
    BOOST_CHECK_EQUAL(cdf(dist, mu), RealType(0.5L));
    BOOST_CHECK_EQUAL(cdf(complement(dist, mu)), RealType(0.5L));
    BOOST_CHECK_EQUAL(cdf(dist, upper), RealType(1));
    BOOST_CHECK_EQUAL(cdf(complement(dist, upper)), RealType(0));
  }
  const inverse_gaussian_distribution<RealType> concentrated(1, maximum);
  BOOST_CHECK_EQUAL(cdf(concentrated, minimum), RealType(0));
  BOOST_CHECK_EQUAL(cdf(complement(concentrated, minimum)), RealType(1));
  BOOST_CHECK_EQUAL(cdf(concentrated, maximum), RealType(1));
  BOOST_CHECK_EQUAL(cdf(complement(concentrated, maximum)), RealType(0));
  const inverse_gaussian_distribution<RealType> equal_large(maximum / 2, maximum / 2);
  const RealType at_mean = static_cast<RealType>(0.668102001223170606427149114029083698375984L);
  BOOST_CHECK_CLOSE_FRACTION(cdf(equal_large, maximum / 2), at_mean, 128 * epsilon);
  BOOST_CHECK_CLOSE_FRACTION(cdf(complement(equal_large, maximum / 2)), 1 - at_mean, 128 * epsilon);
}

BOOST_AUTO_TEST_CASE( test_finite_extremes )
{
  test_finite_extremes_cdf(0.0F);
  test_finite_extremes_cdf(0.0);
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
  test_finite_extremes_cdf(0.0L);
#ifndef BOOST_MATH_NO_REAL_CONCEPT_TESTS
  test_finite_extremes_cdf(boost::math::concepts::real_concept(0.));
#endif
#endif
}

BOOST_AUTO_TEST_CASE( test_main )
{
  using boost::math::inverse_gaussian;
  using boost::math::inverse_gaussian_distribution;

  //int precision = 17; // std::numeric_limits<double::max_digits10;
  double tolfeweps = numeric_limits<double>::epsilon() * 5;
  //double tol6decdigits = numeric_limits<float>::epsilon() * 2;
  // Check that can generate inverse_gaussian distribution using the two convenience methods:
  boost::math::inverse_gaussian w12(1., 2); // Using typedef
  inverse_gaussian_distribution<> w23(2., 3); // Using default RealType double.
  boost::math::inverse_gaussian w11; // Use default unity values for mean and scale.
  // Note NOT myn01() as the compiler will interpret as a function!
  BOOST_CHECK_EQUAL(w11.mean(), 1);
  BOOST_CHECK_EQUAL(w11.scale(), 1);
  BOOST_CHECK_EQUAL(w23.mean(), 2);
  BOOST_CHECK_EQUAL(w23.scale(), 3);
  BOOST_CHECK_EQUAL(w23.shape(), 1.5L);

  // Check the synonyms, provided to allow generic use of find_location and find_scale.
  BOOST_CHECK_EQUAL(w11.mean(), w11.location());
  BOOST_CHECK_EQUAL(w11.scale(), w11.scale());

  BOOST_CHECK_CLOSE_FRACTION(mean(w11), static_cast<double>(1), tolfeweps); // Default mean == unity
  BOOST_CHECK_CLOSE_FRACTION(scale(w11), static_cast<double>(1), tolfeweps); // Default mean == unity

  // median
  // (test double because fails for real_concept because numeric_limits<real_concept>::digits = 0)
  BOOST_CHECK_CLOSE_FRACTION(median(w11),
    static_cast<double>(0.67584130569523893), tolfeweps);
  BOOST_CHECK_CLOSE_FRACTION(median(w23),
    static_cast<double>(1.5122506636053668), tolfeweps);
  
  // Initial spot tests using double values from R.
  // library(SuppDists)
  // formatC(SuppDists::dinverse_gaussian(1, 1, 1), digits=17) ...
  BOOST_CHECK_CLOSE_FRACTION( //  x = 1
    pdf(w11, 1.), static_cast<double>(0.3989422804014327), tolfeweps); // pdf
  BOOST_CHECK_CLOSE_FRACTION( //  x = 1
    logpdf(w11, 1.), static_cast<double>(log(0.3989422804014327)), tolfeweps); // logpdf
  BOOST_CHECK_CLOSE_FRACTION(
    cdf(w11, 1.), static_cast<double>(0.66810200122317065), 10 * tolfeweps); // cdf

  BOOST_CHECK_CLOSE_FRACTION(
    pdf(w11, 0.1), static_cast<double>(0.21979480031862672), tolfeweps); // pdf
  BOOST_CHECK_CLOSE_FRACTION(
    logpdf(w11, 0.1), static_cast<double>(log(0.21979480031862672)), tolfeweps); // logpdf
  BOOST_CHECK_CLOSE_FRACTION(
    cdf(w11, 0.1), static_cast<double>(0.0040761113207110162), 10 * tolfeweps); // cdf

  BOOST_CHECK_CLOSE_FRACTION( // small x
    pdf(w11, 0.01), static_cast<double>(2.0811768202028392e-19), tolfeweps); // pdf
  BOOST_CHECK_CLOSE_FRACTION( // small x
    logpdf(w11, 0.01), static_cast<double>(log(2.0811768202028392e-19)), tolfeweps); // logpdf
  BOOST_CHECK_CLOSE_FRACTION(
    cdf(w11, 0.01), static_cast<double>(4.122313403318778e-23), 10 * tolfeweps); // cdf

  BOOST_CHECK_CLOSE_FRACTION( // smaller x
    pdf(w11, 0.001), static_cast<double>(2.4420044378793562e-213),  tolfeweps); // pdf
  BOOST_CHECK_CLOSE_FRACTION( // smaller x
    logpdf(w11, 0.001), static_cast<double>(log(2.4420044378793562e-213)),  tolfeweps); // pdf
  BOOST_CHECK_CLOSE_FRACTION(
    cdf(w11, 0.001), static_cast<double>(4.8791443010851493e-219), 1000 * tolfeweps); // cdf
  // 4.8791443010859224e-219 versus 4.8791443010851493e-219 so still 14 decimal digits.

  BOOST_CHECK_CLOSE_FRACTION(
    quantile(w11, 0.66810200122317065), static_cast<double>(1.), 1 * tolfeweps); // cdf
  BOOST_CHECK_CLOSE_FRACTION(
    quantile(w11, 0.0040761113207110162), static_cast<double>(0.1), 1 * tolfeweps); // cdf
  BOOST_CHECK_CLOSE_FRACTION(
    quantile(w11, 4.122313403318778e-23), 0.01, 1 * tolfeweps); // quantile
  BOOST_CHECK_CLOSE_FRACTION(
    quantile(w11, 2.4420044378793562e-213), 0.001, 0.03); // quantile
  // quantile 0.001026926242348481 compared to expected 0.001, so much less accurate,
  // but better than R that gives up completely!
  // R Error in SuppDists::qinverse_gaussian(4.87914430108515e-219, 1, 1) : Infinite value in NewtonRoot()

  inverse_gaussian w_big(66.99652081);
  BOOST_CHECK_CLOSE_FRACTION(
     quantile(w_big, 0.97969), 591.567880739988823, 15 * tolfeweps);
  BOOST_CHECK_CLOSE_FRACTION(
     quantile(complement(w_big, 1 - 0.97969)), 591.567880739988823, 10 * tolfeweps);


  BOOST_CHECK_CLOSE_FRACTION(
    pdf(w11, 0.5), static_cast<double>(0.87878257893544476), tolfeweps); // pdf
  BOOST_CHECK_CLOSE_FRACTION(
    logpdf(w11, 0.5), static_cast<double>(log(0.87878257893544476)), tolfeweps); // logpdf
  BOOST_CHECK_CLOSE_FRACTION(
    cdf(w11, 0.5), static_cast<double>(0.3649755481729598), tolfeweps); // cdf

  BOOST_CHECK_CLOSE_FRACTION(
    pdf(w11, 2), static_cast<double>(0.10984782236693059), tolfeweps); // pdf
  BOOST_CHECK_CLOSE_FRACTION(
    logpdf(w11, 2), static_cast<double>(log(0.10984782236693059)), tolfeweps); // logpdf
  BOOST_CHECK_CLOSE_FRACTION(
    cdf(w11, 2), static_cast<double>(.88547542598600637), tolfeweps); // cdf

  BOOST_CHECK_CLOSE_FRACTION(
    pdf(w11, 10), static_cast<double>(0.00021979480031862676), tolfeweps); // pdf
  BOOST_CHECK_CLOSE_FRACTION(
    logpdf(w11, 10), static_cast<double>(log(0.00021979480031862676)), tolfeweps); // logpdf
  BOOST_CHECK_CLOSE_FRACTION(
    cdf(w11, 10), static_cast<double>(0.99964958546279115), tolfeweps); // cdf

  BOOST_CHECK_CLOSE_FRACTION(
    pdf(w11, 100), static_cast<double>(2.0811768202028246e-25), tolfeweps); // pdf
  BOOST_CHECK_CLOSE_FRACTION(
    logpdf(w11, 100), static_cast<double>(log(2.0811768202028246e-25)), tolfeweps); // logpdf
  BOOST_CHECK_CLOSE_FRACTION(
    cdf(w11, 100), static_cast<double>(1), tolfeweps); // cdf
  BOOST_CHECK_CLOSE_FRACTION(
    pdf(w11, 1000), static_cast<double>(2.4420044378793564e-222), 10 * tolfeweps); // pdf
  BOOST_CHECK_CLOSE_FRACTION(
    logpdf(w11, 1000), static_cast<double>(log(2.4420044378793564e-222)), 10 * tolfeweps); // pdf
  BOOST_CHECK_CLOSE_FRACTION(
    cdf(w11, 1000), static_cast<double>(1.), tolfeweps); // cdf

  // A few more misc tests, probably not very useful.  
  BOOST_CHECK_CLOSE_FRACTION(
    cdf(w11, 1.), static_cast<double>(0.66810200122317065), tolfeweps); // cdf
  BOOST_CHECK_CLOSE_FRACTION(
    cdf(w11, 0.1), static_cast<double>(0.0040761113207110162), tolfeweps * 5); // cdf
  // 0.0040761113207110162   0.0040761113207110362
  BOOST_CHECK_CLOSE_FRACTION(
    cdf(w11, 0.2), static_cast<double>(0.063753567519976254), tolfeweps * 5); // cdf
  BOOST_CHECK_CLOSE_FRACTION(
    cdf(w11, 0.5), static_cast<double>(0.3649755481729598), tolfeweps); // cdf

  BOOST_CHECK_CLOSE_FRACTION(
    cdf(w11, 0.9), static_cast<double>(0.62502320258649202), tolfeweps); // cdf
  BOOST_CHECK_CLOSE_FRACTION(
    cdf(w11, 0.99), static_cast<double>(0.66408247396139031), tolfeweps); // cdf
  BOOST_CHECK_CLOSE_FRACTION(
    cdf(w11, 0.999), static_cast<double>(0.66770275955311675), tolfeweps); // cdf
  BOOST_CHECK_CLOSE_FRACTION(
    cdf(w11, 10.), static_cast<double>(0.99964958546279115), tolfeweps); // cdf
  BOOST_CHECK_CLOSE_FRACTION(
    cdf(w11, 50.), static_cast<double>(0.99999999999992029), tolfeweps); // cdf

  BOOST_CHECK_CLOSE_FRACTION(
    quantile(w11, 0.3649755481729598), static_cast<double>(0.5), tolfeweps); // quantile
  BOOST_CHECK_CLOSE_FRACTION(
    quantile(w11, 0.62502320258649202), static_cast<double>(0.9), tolfeweps); // quantile
  BOOST_CHECK_CLOSE_FRACTION(
    quantile(w11, 0.0040761113207110162), static_cast<double>(0.1), tolfeweps); // quantile

  // Wald(2,3) tests
  // ===================
  BOOST_CHECK_CLOSE_FRACTION( // formatC(SuppDists::dinvGauss(1, 2, 3), digits=17) "0.47490884963330904"
    pdf(w23, 1.), static_cast<double>(0.47490884963330904), tolfeweps ); // pdf
  BOOST_CHECK_CLOSE_FRACTION(
    logpdf(w23, 1.), static_cast<double>(log(0.47490884963330904)), tolfeweps ); // logpdf
  BOOST_CHECK_CLOSE_FRACTION(
    pdf(w23, 0.1), static_cast<double>(2.8854207087665401e-05), tolfeweps * 2); // pdf
  BOOST_CHECK_CLOSE_FRACTION(
    logpdf(w23, 0.1), static_cast<double>(log(2.8854207087665401e-05)), tolfeweps * 2); // logpdf
  //2.8854207087665452e-005 2.8854207087665401e-005
  BOOST_CHECK_CLOSE_FRACTION(
    pdf(w23, 10.), static_cast<double>(0.0019822751498574636), tolfeweps); // pdf
  BOOST_CHECK_CLOSE_FRACTION(
    logpdf(w23, 10.), static_cast<double>(log(0.0019822751498574636)), tolfeweps); // logpdf
  BOOST_CHECK_CLOSE_FRACTION(
    pdf(w23, 10.), static_cast<double>(0.0019822751498574636), tolfeweps); // pdf
  BOOST_CHECK_CLOSE_FRACTION(
    logpdf(w23, 10.), static_cast<double>(log(0.0019822751498574636)), tolfeweps); // logpdf

  // Bigger changes in mean and scale.

  inverse_gaussian w012(0.1, 2);
  BOOST_CHECK_CLOSE_FRACTION(
    pdf(w012, 1.), static_cast<double>(3.7460367141230404e-36), tolfeweps ); // pdf
  BOOST_CHECK_CLOSE_FRACTION(
    logpdf(w012, 1.), static_cast<double>(log(3.7460367141230404e-36)), tolfeweps ); // logpdf
  BOOST_CHECK_CLOSE_FRACTION(
    cdf(w012, 1.), static_cast<double>(1), tolfeweps ); // pdf

  inverse_gaussian w0110(0.1, 10);
  BOOST_CHECK_CLOSE_FRACTION(
    pdf(w0110, 1.), static_cast<double>(1.6279643678071011e-176), 100 * tolfeweps ); // pdf
  BOOST_CHECK_CLOSE_FRACTION(
    logpdf(w0110, 1.), static_cast<double>(log(1.6279643678071011e-176)), 100 * tolfeweps ); // logpdf
  BOOST_CHECK_CLOSE_FRACTION(
    cdf(w0110, 1.), static_cast<double>(1), tolfeweps ); // cdf
  BOOST_CHECK_CLOSE_FRACTION(
     cdf(complement(w0110, 1.)), static_cast<double>(3.2787685715328683e-179), 1e6 * tolfeweps ); // cdf complement
  // Differs because of loss of accuracy.

  BOOST_CHECK_CLOSE_FRACTION(
    pdf(w0110, 0.1), static_cast<double>(39.894228040143268), tolfeweps ); // pdf
  BOOST_CHECK_CLOSE_FRACTION(
    logpdf(w0110, 0.1), static_cast<double>(log(39.894228040143268)), tolfeweps ); // logpdf
  BOOST_CHECK_CLOSE_FRACTION(
    cdf(w0110, 0.1), static_cast<double>(0.51989761564832704), 10 * tolfeweps ); // cdf

    // Basic sanity-check spot values for all floating-point types..
  // (Parameter value, arbitrarily zero, only communicates the floating point type).
  test_spots(0.0F); // Test float. OK at decdigits = 0 tolerance = 0.0001 %
  test_spots(0.0); // Test double. OK at decdigits 7, tolerance = 1e07 %
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
  test_spots(0.0L); // Test long double.
#ifndef BOOST_MATH_NO_REAL_CONCEPT_TESTS
  test_spots(boost::math::concepts::real_concept(0.)); // Test real concept.
#endif
#else
  std::cout << "<note>The long double tests have been disabled on this platform "
    "either because the long double overloads of the usual math functions are "
    "not available at all, or because they are too inaccurate for these tests "
    "to pass.</note>" << std::endl;
#endif
  /*      */
  
} // BOOST_AUTO_TEST_CASE( test_main )

/*

Output:


*/
