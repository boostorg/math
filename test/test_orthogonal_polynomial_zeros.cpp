//  Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "math_unit_test.hpp"
#include <cmath>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <vector>
#include <boost/math/special_functions/hermite.hpp>
#include <boost/math/special_functions/laguerre.hpp>
#include <boost/math/special_functions/jacobi.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#ifndef BOOST_MATH_NO_CXX17_IF_CONSTEXPR
#include <boost/math/quadrature/gauss_hermite.hpp>
#include <boost/math/quadrature/gauss_laguerre.hpp>
#endif

using boost::math::hermite_zeros;
using boost::math::laguerre_zeros;
using boost::math::jacobi_zeros;
using boost::multiprecision::cpp_bin_float_50;
// The zeros are checked with Newton's method on the polynomials in this type.
using high = boost::multiprecision::cpp_bin_float_100;

void parse(const char* s, float& x) { x = std::strtof(s, nullptr); }
void parse(const char* s, double& x) { x = std::strtod(s, nullptr); }
template <class T>
void parse(const char* s, T& x) { x = T(s); }

template <class T>
T value(const char* s)
{
   T x;
   parse(s, x);
   return x;
}

// Zeros from mpmath's polyroots on the coefficients from the recurrences, in exact arithmetic, at 150 digits.
struct hermite_case { unsigned n; std::vector<const char*> zeros; };
struct laguerre_case { unsigned n; const char* alpha; std::vector<const char*> zeros; };
struct jacobi_case { unsigned n; const char* alpha; const char* beta; std::vector<const char*> zeros; };

const std::vector<hermite_case> hermite_cases = {
      { 1, { "0.0" } },
      { 6, { "0.4360774119276165086792159482506249107662566689318268059", "1.335849074013696949714895282970370673972589526125010814", "2.350604973674492222833921987060920846537863827294099092" } },
      { 7, { "0.0", "0.8162878828589646630387109590271458167428894003786361568", "1.673551628767471445031801398303594819107810057735408927", "2.651961356835233492447082006516616114438158478625529417" } },
      { 20, { "0.2453407083009012499038365306336166239661338513034857349", "0.7374737285453943587056051442521042290772162039767941873", "1.234076215395323007885818346959410229585445930069362387", "1.738537712116586206780865662136406442951409419600399894", "2.254974002089275523082333344734565128082265316026396151", "2.788806058428130480525033756403185410670698887902043974", "3.347854567383216326914924522996463698510478590293677057", "3.944764040115625210375628800524411807149768127888313021", "4.603682449550744273077675248978347585113398487761923654", "5.38748089001123286201690041068112075399628644906591489" } },
};

const std::vector<laguerre_case> laguerre_cases = {
      { 1, "0", { "1.0" } },
      { 10, "0", { "0.1377934705404924308307725056527111881079916807457832942", "0.729454549503170498160373121676078781076072733312249246", "1.808342901740316048232920075750608833283060282371448152", "3.401433697854899514482532221408390679273156614203473229", "5.552496140063803632417558486868762857974064287317813828", "8.330152746764496700238767197274522182709438972030989345", "11.84378583790006556491853891914161398580281690946533587", "16.27925783137810209953265393583362233525599560303305908", "21.99658581198076195127709019559449397680673234000188779", "29.92069701227389155990879334079919517971067057751796017" } },
      { 8, "0.5", { "0.2826336481165991273853395912392153400317888518249301018", "1.139873801581613696920923375954489917014082188982960013", "2.601524843406029388441746428168779837605977211502989905", "4.724114537527790745575326684891977598902584439580712876", "7.605256299231613702365957655278871109340770299041742325", "11.41718207654583007908074696912694469890867755572885602", "16.49941079765581552562593445152671228778584349232290796", "23.7300039959347077346040248438130092104102759610149008" } },
      { 6, "-0.75", { "0.04573846887456520038620989561492032168122358426547260214", "0.7565209910189564466591515188946810321340892868201833925", "2.332774721240584867170333042476367544765865975492675727", "4.9043829130172385179626546590353431609376933642702237", "8.763281241957814470464143886421113679626392652705657415", "14.69730166389084049735750699755757426085473513644578716" } },
      { 12, "7.25", { "2.043758150174322591200592717784207822625618308304800389", "3.64718495206882283124005315221826034011272656863142807", "5.575597754279745303521746639232348961399142881711056736", "7.872192335862900590557439817561819120335257336861408761", "10.57416475679685333249976666062364370147795985778201437", "13.72658506718814656530926123361730009366667687863342647", "17.38980836498542128124097505393375852959252362125014268", "21.64946129816789861618352101477528118586033167374316479", "26.63486956479782523857710351117562983743330213816664901", "32.55955622282857313368274043834394910550846101755621263", "39.83152790701440407236240841409148244203427925887579985", "49.49529362583508644362439134664231885995372045848389624" } },
      { 5, "2", { "1.031109144093381623645570337602608646127321629493926668", "2.837212823953821742048730848650826181089626364845069705", "5.620294272598707887045331611974069093873903434501008406", "9.682909837664027130647489393135938476572439587389522246", "15.82847392169006161661287780863655760233670898377047297" } },
};

const std::vector<jacobi_case> jacobi_cases = {
      { 1, "0.5", "-0.25", { "-0.3333333333333333333333333333333333333333333333333333333" } },
      { 7, "0.5", "-0.25", { "-0.9655920538134315882616228274345519487842059231148001697", "-0.7827112226512704387664840229775231549219070348824486075", "-0.4688025808677373384204836390266046588776741913904261808", "-0.07631191305758712368332900594327470750340950314833799915", "0.3290815045839580861730936900104315615652550129272236992", "0.6795310452513747125961114755008135964194259283087286156", "0.916384167923114742994293277239130364734094658668481695" } },
      { 10, "-0.75", "2.5", { "-0.873804943020537125196403522735511394175118098708810719", "-0.6957101323762201584473665608991765343988903535509152933", "-0.4661623745701520014145407018286872858255265259102531779", "-0.202240673289828344522110074907991089824321169434614337", "0.07612257657064914208109426007100732140603393528111644107", "0.347853249920195295634006053997591856044673681182151568", "0.5923656663132282918825213480893196820217778333046081117", "0.7911306237246019776851089162728341925726325606724894739", "0.9290770430278378419938374833651644494121350295933970337", "0.995621837263443471108450499724874090122924946651290669" } },
      { 6, "3", "3", { "-0.7755355259865794571856760657808861808477886486064745923", "-0.5009134932118855608610636756018708834606319061246200168", "-0.1731567539181085648861776049085975765417940108723871403", "0.1731567539181085648861776049085975765417940108723871403", "0.5009134932118855608610636756018708834606319061246200168", "0.7755355259865794571856760657808861808477886486064745923" } },
      { 5, "-0.5", "-0.5", { "-0.9510565162951535721164393333793821434056986341257502224", "-0.5877852522924731291687059546390727685976524376431459911", "0.0", "0.5877852522924731291687059546390727685976524376431459911", "0.9510565162951535721164393333793821434056986341257502224" } },
};

template <class T>
void check_against(const std::vector<const char*>& expected, const std::vector<T>& zeros, int ulps)
{
   CHECK_EQUAL(zeros.size(), expected.size());
   for (std::size_t i = 0; i < (std::min)(zeros.size(), expected.size()); ++i)
   {
      const T x = value<T>(expected[i]);
      if (x == 0)
         CHECK_ABSOLUTE_ERROR(x, zeros[i], ulps * std::numeric_limits<T>::epsilon());
      else
         CHECK_ULP_CLOSE(x, zeros[i], ulps);
   }
}

template <class T>
void test_reference_values()
{
   for (const auto& c : hermite_cases)
      check_against(c.zeros, hermite_zeros<T>(c.n), 1);
   for (const auto& c : laguerre_cases)
      check_against(c.zeros, laguerre_zeros<T>(c.n, value<T>(c.alpha)), 2);
   for (const auto& c : jacobi_cases)
      check_against(c.zeros, jacobi_zeros<T>(c.n, value<T>(c.alpha), value<T>(c.beta)), 2);
}

// The zeros are within `ulps` of those of p, by one step of Newton's method in high precision: p(x) / p'(x) is
// the distance to the zero to first order. The ulps are those of max(|x|, smallest): the zero of a symmetric Jacobi
// polynomial at 0 comes out as a few eps squared.
template <class T, class P, class Prime>
void check_newton_step(const std::vector<T>& zeros, P p, Prime p_prime, int ulps, double smallest = 0)
{
   using std::abs;
   for (const T& z : zeros)
   {
      const high x = static_cast<high>(z);
      const high step = p(x) / p_prime(x);
      CHECK_LE(static_cast<double>(abs(step)), ulps * static_cast<double>(std::numeric_limits<T>::epsilon()) * (std::max)(static_cast<double>(abs(z)), smallest));
   }
}

template <class T>
void test_polynomial_vanishes()
{
   for (unsigned n = 1; n <= 40; n += 3)
   {
      check_newton_step(hermite_zeros<T>(n), [n](const high& x) { return boost::math::hermite(n, x); },
                        [n](const high& x) { return 2 * n * boost::math::hermite(n - 1, x); }, 1);
      for (unsigned m : { 0u, 1u, 4u })
      {  // The derivative of L_n^m is -L_{n-1}^(m+1).
         check_newton_step(laguerre_zeros<T>(n, T(m)), [n, m](const high& x) { return boost::math::laguerre(n, m, x); },
                           [n, m](const high& x) { return -boost::math::laguerre(n - 1, m + 1, x); }, 2);
      }
      for (const auto& ab : { std::make_pair(0.5, -0.25), std::make_pair(-0.75, 2.5), std::make_pair(4.0, 4.0) })
      {
         const T alpha = static_cast<T>(ab.first), beta = static_cast<T>(ab.second);
         const high a = static_cast<high>(alpha), b = static_cast<high>(beta);
         check_newton_step(jacobi_zeros<T>(n, alpha, beta), [=](const high& x) { return boost::math::jacobi(n, a, b, x); },
                           [=](const high& x) { return boost::math::jacobi_prime(n, a, b, x); }, 2, 1);
      }
   }
}

// The count, the order, the domain, and that the zeros of p_n separate those of p_{n+1}.
template <class T>
void check_interlacing(const std::vector<T>& x, const std::vector<T>& y)
{
   CHECK_EQUAL(y.size(), x.size() + 1);
   for (std::size_t i = 0; i + 1 < y.size() && i < x.size(); ++i)
   {
      CHECK_LE(y[i], x[i]);
      CHECK_LE(x[i], y[i + 1]);
      CHECK_TRUE(y[i] < y[i + 1]);
   }
}

template <class T>
std::vector<T> mirrored(const std::vector<T>& nonnegative)
{
   std::vector<T> all;
   for (std::size_t i = nonnegative.size(); i-- > 0;)
   {
      if (nonnegative[i] != 0)
         all.push_back(-nonnegative[i]);
   }
   all.insert(all.end(), nonnegative.begin(), nonnegative.end());
   return all;
}

template <class T>
void test_structure(unsigned max_n)
{
   using std::sqrt;
   CHECK_EQUAL(hermite_zeros<T>(0).size(), std::size_t(0));
   CHECK_EQUAL(laguerre_zeros<T>(0).size(), std::size_t(0));
   CHECK_EQUAL(jacobi_zeros<T>(0, T(1), T(2)).size(), std::size_t(0));
   for (unsigned n = 1; n < max_n; ++n)
   {
      const std::vector<T> h = hermite_zeros<T>(n);
      CHECK_EQUAL(h.size(), std::size_t((n + 1) / 2));
      CHECK_TRUE(h.front() >= 0);
      CHECK_TRUE(h.back() < sqrt(T(2 * n + 1)));
      CHECK_EQUAL(h.front() == 0, n % 2 == 1);
      check_interlacing(mirrored(h), mirrored(hermite_zeros<T>(n + 1)));

      for (T alpha : { T(0), T(-0.5), T(2.75) })
      {
         const std::vector<T> l = laguerre_zeros<T>(n, alpha);
         CHECK_TRUE(l.front() > 0);
         check_interlacing(l, laguerre_zeros<T>(n + 1, alpha));
      }
      for (const auto& ab : { std::make_pair(0.0, 0.0), std::make_pair(-0.5, 0.5), std::make_pair(-0.875, 6.0) })
      {
         const T alpha = static_cast<T>(ab.first), beta = static_cast<T>(ab.second);
         const std::vector<T> j = jacobi_zeros<T>(n, alpha, beta);
         CHECK_TRUE(j.front() > -1);
         CHECK_TRUE(j.back() < 1);
         check_interlacing(j, jacobi_zeros<T>(n + 1, alpha, beta));
      }
   }
}

#ifndef BOOST_MATH_NO_CXX17_IF_CONSTEXPR
// The tabulated rules, for 10 points.
template <class T>
void test_quadrature_tables()
{
   const auto& h = boost::math::quadrature::gauss_hermite<T, 10>::abscissa();
   const boost::math::quadrature::gauss_laguerre<T, 10> laguerre_rule;
   const auto& l = laguerre_rule.abscissa();
   const std::vector<T> hz = hermite_zeros<T>(10), lz = laguerre_zeros<T>(10);
   CHECK_EQUAL(hz.size(), std::size_t(h.size()));
   CHECK_EQUAL(lz.size(), std::size_t(l.size()));
   for (std::size_t i = 0; i < (std::min)(hz.size(), std::size_t(h.size())); ++i)
      CHECK_ULP_CLOSE(static_cast<T>(h[i]), hz[i], 1);
   for (std::size_t i = 0; i < (std::min)(lz.size(), std::size_t(l.size())); ++i)
      CHECK_ULP_CLOSE(static_cast<T>(l[i]), lz[i], 1);
}
#endif

template <class T>
void test_errors()
{
   using std::numeric_limits;
   const T nan = numeric_limits<T>::quiet_NaN();
   const T inf = numeric_limits<T>::infinity();
   for (T bad : { T(-1), T(-2), nan, inf })
   {
      CHECK_THROW(laguerre_zeros<T>(3, bad), std::domain_error);
      CHECK_THROW(jacobi_zeros<T>(3, bad, T(0)), std::domain_error);
      CHECK_THROW(jacobi_zeros<T>(3, T(0), bad), std::domain_error);
   }
   using namespace boost::math::policies;
   using ignore = policy<domain_error<ignore_error>>;
   const std::vector<T> l = laguerre_zeros<T>(3, T(-1), ignore());
   CHECK_EQUAL(l.size(), std::size_t(3));
   CHECK_NAN(l[0]);
   const std::vector<T> j = jacobi_zeros<T>(4, T(0), nan, ignore());
   CHECK_EQUAL(j.size(), std::size_t(4));
   CHECK_NAN(j[3]);
}

int main()
{
   test_reference_values<float>();
   test_reference_values<double>();
   test_reference_values<cpp_bin_float_50>();

   test_polynomial_vanishes<float>();
   test_polynomial_vanishes<double>();
   test_polynomial_vanishes<cpp_bin_float_50>();

   test_structure<float>(40);
   test_structure<double>(40);
   test_structure<cpp_bin_float_50>(12);

   // Many zeros, to O(n) work: the largest zero of H_n is below sqrt(2n + 1).
   const std::vector<double> h = hermite_zeros<double>(2000);
   CHECK_EQUAL(h.size(), std::size_t(1000));
   CHECK_TRUE(h.back() < std::sqrt(4001.0));
   check_newton_step(std::vector<double>{ h.front(), h.back() }, [](const high& x) { return boost::math::hermite(2000, x); },
                     [](const high& x) { return 4000 * boost::math::hermite(1999, x); }, 4);

#ifndef BOOST_MATH_NO_CXX17_IF_CONSTEXPR
   test_quadrature_tables<float>();
   test_quadrature_tables<double>();
   test_quadrature_tables<cpp_bin_float_50>();
#endif

   test_errors<float>();
   test_errors<double>();
   test_errors<cpp_bin_float_50>();

   return boost::math::test::report_errors();
}
