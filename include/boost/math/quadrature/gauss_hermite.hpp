//  Copyright Jacob Hass 2026.
//  Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_QUADRATURE_GAUSS_HERMITE_HPP
#define BOOST_MATH_QUADRATURE_GAUSS_HERMITE_HPP

#ifdef _MSC_VER
#pragma once
#endif

#include <boost/math/tools/config.hpp>

#ifdef BOOST_MATH_NO_CXX17_IF_CONSTEXPR
#error "The header <boost/math/quadrature/gauss_hermite.hpp> requires C++17 or later."
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <array>
#include <cmath>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>
#endif
#include <boost/math/constants/constants.hpp>
#include <boost/math/policies/policy.hpp>
#include <boost/math/tools/big_constant.hpp>
#include <boost/math/tools/precision.hpp>

namespace boost { namespace math { namespace quadrature { namespace detail {

template <class T>
struct gauss_hermite_constant_category
{
   static const unsigned value =
      (std::numeric_limits<T>::is_specialized == 0) ? 999 :
      (std::numeric_limits<T>::radix == 2) ?
      (
#ifdef BOOST_HAS_FLOAT128
         (std::numeric_limits<T>::digits <= 113) && std::is_constructible<T, __float128>::value ? 0 :
#else
         (std::numeric_limits<T>::digits <= std::numeric_limits<long double>::digits) && std::is_constructible<T, long double>::value ? 0 :
#endif
         (std::numeric_limits<T>::digits10 <= 110) && std::is_constructible<T, const char*>::value ? 4 : 999
      ) : (std::numeric_limits<T>::digits10 <= 110) && std::is_constructible<T, const char*>::value ? 4 : 999;

   using storage_type =
      std::conditional_t<(std::numeric_limits<T>::is_specialized == 0), T,
         std::conditional_t<(std::numeric_limits<T>::radix == 2),
            std::conditional_t< ((std::numeric_limits<T>::digits <= std::numeric_limits<float>::digits) && std::is_constructible<T, float>::value),
               float,
               std::conditional_t<((std::numeric_limits<T>::digits <= std::numeric_limits<double>::digits) && std::is_constructible<T, double>::value),
                  double,
                  std::conditional_t<((std::numeric_limits<T>::digits <= std::numeric_limits<long double>::digits) && std::is_constructible<T, long double>::value),
                     long double,
#ifdef BOOST_HAS_FLOAT128
                     std::conditional_t<((std::numeric_limits<T>::digits <= 113) && std::is_constructible<T, __float128>::value),
                        __float128,
                        T
                     >
                  >
#else
                     T
                  >
#endif
               >
            >, T
         >
      >;
};

#ifndef BOOST_MATH_GAUSS_NO_COMPUTE_ON_DEMAND

template <class Real, unsigned N, unsigned Category>
class gauss_hermite_detail
{
   static_assert(N > 0, "Gauss-Hermite quadrature needs at least one point.");

   // The nonnegative zeros of the Hermite polynomial H_N, in increasing order, and their weights.
   // Everything is evaluated with the orthonormal recurrence
   //    p_0 = pi^(-1/4),  p_{k+1} = x sqrt(2/(k+1)) p_k - sqrt(k/(k+1)) p_{k-1},
   // whose values stay representable long after H_N itself overflows. Then p_N' = sqrt(2N) p_{N-1},
   // and the weight is 2/p_N'(x)^2. Evaluating H_N from its monomial coefficients instead
   // would be hopelessly ill-conditioned: the largest zero of H_100 has condition number ~1e13.
   //
   // By Sturm's theorem for orthogonal polynomials, the number of sign changes in p_0(x), ..., p_N(x)
   // is the number of zeros of p_N above x. Bisecting on that count isolates each zero, and Newton's
   // method then converges to it inside its bracket. Asymptotic initial guesses alone are not enough:
   // for N = 200, those of Numerical Recipes' gauher send Newton to the wrong zero.
   struct recurrence
   {
      Real p;
      Real p_prime;
      unsigned zeros_above;
   };

   static recurrence evaluate(Real x)
   {
      using std::sqrt;
      Real p = 1 / sqrt(sqrt(boost::math::constants::pi<Real>()));
      Real p_previous = 0;
      bool last_negative = false;
      unsigned sign_changes = 0;
      for (unsigned k = 0; k < N; ++k)
      {
         Real p_next = x * sqrt(Real(2) / Real(k + 1)) * p - sqrt(Real(k) / Real(k + 1)) * p_previous;
         p_previous = p;
         p = p_next;
         // A zero p_k(x) lies between values of opposite sign, so skipping it keeps the count right.
         if ((p != 0) && ((p < 0) != last_negative))
         {
            ++sign_changes;
            last_negative = !last_negative;
         }
      }
      return recurrence{ p, sqrt(Real(2 * N)) * p_previous, sign_changes };
   }

   // The weight is 2/p_N'(x)^2 at the exact zero, but d(log w)/dx is about -4x, so evaluating it at
   // the rounded zero z costs far more accuracy than the rounding itself. Since p'' = 2x p' - 2N p,
   // move p' to the exact zero, which lies a Newton step delta = -p/p' away.
   static Real weight(Real z, const recurrence& r)
   {
      Real delta = -r.p / r.p_prime;
      Real p_prime = r.p_prime + (2 * z * r.p_prime - 2 * Real(N) * r.p) * delta;
      return 2 / (p_prime * p_prime);
   }

   static std::pair<std::vector<Real>, std::vector<Real> > calculate_values()
   {
      using std::abs;
      using std::sqrt;
      const unsigned positive_zeros = N / 2;
      const unsigned offset = N & 1;
      std::vector<Real> x(positive_zeros + offset), w(positive_zeros + offset);
      if (offset)
      {
         recurrence r = evaluate(Real(0));
         w[0] = weight(Real(0), r);
      }
      // All zeros lie below sqrt(2N + 1), so none lie above `upper`. Find them from the largest down:
      // after each one, `upper` is a point with exactly t - 1 zeros above it.
      Real upper = sqrt(Real(2 * N + 2));
      for (unsigned t = 1; t <= positive_zeros; ++t)
      {
         Real lower = 0;
         Real high = upper;
         for (int i = 0; i < tools::digits<Real>(); ++i)
         {
            Real mid = (lower + high) / 2;
            unsigned above = evaluate(mid).zeros_above;
            if (above >= t)
               lower = mid;
            else
               high = mid;
            if (above == t)
               break;
         }
         // Newton's method, falling back to bisection whenever a step would leave the bracket.
         const bool lower_negative = evaluate(lower).p < 0;
         Real z = (lower + high) / 2;
         recurrence r = evaluate(z);
         for (int i = 0; (r.p != 0) && (i < 4 * tools::digits<Real>()); ++i)
         {
            if ((r.p < 0) == lower_negative)
               lower = z;
            else
               high = z;
            Real next = z - r.p / r.p_prime;
            if (!((next > lower) && (next < high)))
               next = (lower + high) / 2;
            const bool converged = abs(next - z) <= 2 * tools::epsilon<Real>() * abs(next);
            z = next;
            r = evaluate(z);
            if (converged)
               break;
         }
         x[positive_zeros + offset - t] = z;
         w[positive_zeros + offset - t] = weight(z, r);
         upper = lower;
      }
      return std::make_pair(x, w);
   }

   static const std::pair<std::vector<Real>, std::vector<Real> >& values()
   {
      static const std::pair<std::vector<Real>, std::vector<Real> > data = calculate_values();
      return data;
   }

public:
   static const std::vector<Real>& abscissa()
   {
      return values().first;
   }
   static const std::vector<Real>& weights()
   {
      return values().second;
   }
};

#else

template <class Real, unsigned N, unsigned Category>
class gauss_hermite_detail;

#endif

#ifndef BOOST_HAS_FLOAT128
template <class T>
class gauss_hermite_detail<T, 7, 0>
{
   using storage_type = typename gauss_hermite_constant_category<T>::storage_type;
   public:
      static std::array<storage_type, 4> const & abscissa()
      {
         static std::array<storage_type, 4> data = {
            static_cast<storage_type>(0.00000000000000000000000000000000000e+00L),
            static_cast<storage_type>(8.16287882858964663038710959027145817e-01L),
            static_cast<storage_type>(1.67355162876747144503180139830359482e+00L),
            static_cast<storage_type>(2.65196135683523349244708200651661611e+00L),
         };
         return data;
      }
      static std::array<storage_type, 4> const & weights()
      {
         static std::array<storage_type, 4> data = {
            static_cast<storage_type>(8.10264617556807326764876563813094941e-01L),
            static_cast<storage_type>(4.25607252610127800520317466666391035e-01L),
            static_cast<storage_type>(5.45155828191270305921785688416951260e-02L),
            static_cast<storage_type>(9.71781245099519154149424255938959644e-04L),
         };
         return data;
      }
   };

#else
template <class T>
class gauss_hermite_detail<T, 7, 0>
{
   using storage_type = typename gauss_hermite_constant_category<T>::storage_type;
   public:
      static std::array<storage_type, 4> const & abscissa()
      {
         static std::array<storage_type, 4> data = {
            static_cast<storage_type>(0.00000000000000000000000000000000000e+00Q),
            static_cast<storage_type>(8.16287882858964663038710959027145817e-01Q),
            static_cast<storage_type>(1.67355162876747144503180139830359482e+00Q),
            static_cast<storage_type>(2.65196135683523349244708200651661611e+00Q),
         };
         return data;
      }
      static std::array<storage_type, 4> const & weights()
      {
         static std::array<storage_type, 4> data = {
            static_cast<storage_type>(8.10264617556807326764876563813094941e-01Q),
            static_cast<storage_type>(4.25607252610127800520317466666391035e-01Q),
            static_cast<storage_type>(5.45155828191270305921785688416951260e-02Q),
            static_cast<storage_type>(9.71781245099519154149424255938959644e-04Q),
         };
         return data;
      }
   };

#endif
template <class T>
class gauss_hermite_detail<T, 7, 4>
{
   using storage_type = typename gauss_hermite_constant_category<T>::storage_type;
   public:
      static std::array<T, 4> const & abscissa()
      {
         static std::array<T, 4> data = {
            BOOST_MATH_HUGE_CONSTANT(T, 0, 0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000e+00),
            BOOST_MATH_HUGE_CONSTANT(T, 0, 8.1628788285896466303871095902714581674288940037863615684472203343594907048766511668519794976704116666704491757953733e-01),
            BOOST_MATH_HUGE_CONSTANT(T, 0, 1.6735516287674714450318013983035948191078100577354089269242175099937138333690347986612540761793805055242734953861690e+00),
            BOOST_MATH_HUGE_CONSTANT(T, 0, 2.6519613568352334924470820065166161144381584786255294172031003071471600949016031668052906818697878641942881494964276e+00),
         };
         return data;
      }
      static std::array<T, 4> const & weights()
      {
         static std::array<T, 4> data = {
            BOOST_MATH_HUGE_CONSTANT(T, 0, 8.1026461755680732676487656381309494070745117994166268718345498964704515867018614005712030022333470419027511937110060e-01),
            BOOST_MATH_HUGE_CONSTANT(T, 0, 4.2560725261012780052031746666639103543970682101400152587878953913358027297072602373305673647185542480994656431711885e-01),
            BOOST_MATH_HUGE_CONSTANT(T, 0, 5.4515582819127030592178568841695125960899915859346492986586846631369448562711195465546454603039043598893480308577392e-02),
            BOOST_MATH_HUGE_CONSTANT(T, 0, 9.7178124509951915414942425593895964444240121701420164980001433798334142698580146031198718271051220413580750089473744e-04),
         };
         return data;
      }
   };

#ifndef BOOST_HAS_FLOAT128
template <class T>
class gauss_hermite_detail<T, 10, 0>
{
   using storage_type = typename gauss_hermite_constant_category<T>::storage_type;
   public:
      static std::array<storage_type, 5> const & abscissa()
      {
         static std::array<storage_type, 5> data = {
            static_cast<storage_type>(3.42901327223704608789165025557258031e-01L),
            static_cast<storage_type>(1.03661082978951365417749191675920902e+00L),
            static_cast<storage_type>(1.75668364929988177345140122010615676e+00L),
            static_cast<storage_type>(2.53273167423278979640896079775479348e+00L),
            static_cast<storage_type>(3.43615911883773760332672549431912138e+00L),
         };
         return data;
      }
      static std::array<storage_type, 5> const & weights()
      {
         static std::array<storage_type, 5> data = {
            static_cast<storage_type>(6.10862633735325798783564990433419713e-01L),
            static_cast<storage_type>(2.40138611082314686416523295005861395e-01L),
            static_cast<storage_type>(3.38743944554810631361647312775859737e-02L),
            static_cast<storage_type>(1.34364574678123269220156558584591387e-03L),
            static_cast<storage_type>(7.64043285523262062915936785959522211e-06L),
         };
         return data;
      }
   };

#else
template <class T>
class gauss_hermite_detail<T, 10, 0>
{
   using storage_type = typename gauss_hermite_constant_category<T>::storage_type;
   public:
      static std::array<storage_type, 5> const & abscissa()
      {
         static std::array<storage_type, 5> data = {
            static_cast<storage_type>(3.42901327223704608789165025557258031e-01Q),
            static_cast<storage_type>(1.03661082978951365417749191675920902e+00Q),
            static_cast<storage_type>(1.75668364929988177345140122010615676e+00Q),
            static_cast<storage_type>(2.53273167423278979640896079775479348e+00Q),
            static_cast<storage_type>(3.43615911883773760332672549431912138e+00Q),
         };
         return data;
      }
      static std::array<storage_type, 5> const & weights()
      {
         static std::array<storage_type, 5> data = {
            static_cast<storage_type>(6.10862633735325798783564990433419713e-01Q),
            static_cast<storage_type>(2.40138611082314686416523295005861395e-01Q),
            static_cast<storage_type>(3.38743944554810631361647312775859737e-02Q),
            static_cast<storage_type>(1.34364574678123269220156558584591387e-03Q),
            static_cast<storage_type>(7.64043285523262062915936785959522211e-06Q),
         };
         return data;
      }
   };

#endif
template <class T>
class gauss_hermite_detail<T, 10, 4>
{
   using storage_type = typename gauss_hermite_constant_category<T>::storage_type;
   public:
      static std::array<T, 5> const & abscissa()
      {
         static std::array<T, 5> data = {
            BOOST_MATH_HUGE_CONSTANT(T, 0, 3.4290132722370460878916502555725803120830265867773062820577103907292018192312101560713520305443248063627958675235231e-01),
            BOOST_MATH_HUGE_CONSTANT(T, 0, 1.0366108297895136541774919167592090162982561106572157244372326874036404187243878992980306028243900733070850998719429e+00),
            BOOST_MATH_HUGE_CONSTANT(T, 0, 1.7566836492998817734514012201061567632954744937388471000872849560961705796098190771278266711469838314640623161517203e+00),
            BOOST_MATH_HUGE_CONSTANT(T, 0, 2.5327316742327897964089607977547934803078465081567249459332874134930859380473865340227292325941771646240824732515703e+00),
            BOOST_MATH_HUGE_CONSTANT(T, 0, 3.4361591188377376033267254943191213848406783093901772906871762944002361145649511153468627231618843292611320162500620e+00),
         };
         return data;
      }
      static std::array<T, 5> const & weights()
      {
         static std::array<T, 5> data = {
            BOOST_MATH_HUGE_CONSTANT(T, 0, 6.1086263373532579878356499043341971323859272928403754564871663738953676910082005768001121074412857938126204657377493e-01),
            BOOST_MATH_HUGE_CONSTANT(T, 0, 2.4013861108231468641652329500586139537004249390192739899066056147994315237169314350600241149865696523148172530499056e-01),
            BOOST_MATH_HUGE_CONSTANT(T, 0, 3.3874394455481063136164731277585973698161942481867016739715157699771381277813573402914401584449869233310982466662099e-02),
            BOOST_MATH_HUGE_CONSTANT(T, 0, 1.3436457467812326922015655858459138698692732786257703646162498020762483818413155201782152234375992775169973716992552e-03),
            BOOST_MATH_HUGE_CONSTANT(T, 0, 7.6404328552326206291593678595952221082891147358323631952885551280911633480005783690893185993195845416600950144409630e-06),
         };
         return data;
      }
   };

#ifndef BOOST_HAS_FLOAT128
template <class T>
class gauss_hermite_detail<T, 15, 0>
{
   using storage_type = typename gauss_hermite_constant_category<T>::storage_type;
   public:
      static std::array<storage_type, 8> const & abscissa()
      {
         static std::array<storage_type, 8> data = {
            static_cast<storage_type>(0.00000000000000000000000000000000000e+00L),
            static_cast<storage_type>(5.65069583255575748526020337198188867e-01L),
            static_cast<storage_type>(1.13611558521092066631913490555610213e+00L),
            static_cast<storage_type>(1.71999257518648893241583152515256151e+00L),
            static_cast<storage_type>(2.32573248617385774545404479448534701e+00L),
            static_cast<storage_type>(2.96716692790560324848896036354980632e+00L),
            static_cast<storage_type>(3.66995037340445253472922383311568148e+00L),
            static_cast<storage_type>(4.49999070730939155366438053053482422e+00L),
         };
         return data;
      }
      static std::array<storage_type, 8> const & weights()
      {
         static std::array<storage_type, 8> data = {
            static_cast<storage_type>(5.64100308726417532852625797339963533e-01L),
            static_cast<storage_type>(4.12028687498898627025891079567810319e-01L),
            static_cast<storage_type>(1.58488915795935746883839384959994027e-01L),
            static_cast<storage_type>(3.07800338725460822286814158757801456e-02L),
            static_cast<storage_type>(2.77806884291277589607887049229213490e-03L),
            static_cast<storage_type>(1.00004441232499868127296736176978509e-04L),
            static_cast<storage_type>(1.05911554771106663577520791055016329e-06L),
            static_cast<storage_type>(1.52247580425351702016062666964826618e-09L),
         };
         return data;
      }
   };

#else
template <class T>
class gauss_hermite_detail<T, 15, 0>
{
   using storage_type = typename gauss_hermite_constant_category<T>::storage_type;
   public:
      static std::array<storage_type, 8> const & abscissa()
      {
         static std::array<storage_type, 8> data = {
            static_cast<storage_type>(0.00000000000000000000000000000000000e+00Q),
            static_cast<storage_type>(5.65069583255575748526020337198188867e-01Q),
            static_cast<storage_type>(1.13611558521092066631913490555610213e+00Q),
            static_cast<storage_type>(1.71999257518648893241583152515256151e+00Q),
            static_cast<storage_type>(2.32573248617385774545404479448534701e+00Q),
            static_cast<storage_type>(2.96716692790560324848896036354980632e+00Q),
            static_cast<storage_type>(3.66995037340445253472922383311568148e+00Q),
            static_cast<storage_type>(4.49999070730939155366438053053482422e+00Q),
         };
         return data;
      }
      static std::array<storage_type, 8> const & weights()
      {
         static std::array<storage_type, 8> data = {
            static_cast<storage_type>(5.64100308726417532852625797339963533e-01Q),
            static_cast<storage_type>(4.12028687498898627025891079567810319e-01Q),
            static_cast<storage_type>(1.58488915795935746883839384959994027e-01Q),
            static_cast<storage_type>(3.07800338725460822286814158757801456e-02Q),
            static_cast<storage_type>(2.77806884291277589607887049229213490e-03Q),
            static_cast<storage_type>(1.00004441232499868127296736176978509e-04Q),
            static_cast<storage_type>(1.05911554771106663577520791055016329e-06Q),
            static_cast<storage_type>(1.52247580425351702016062666964826618e-09Q),
         };
         return data;
      }
   };

#endif
template <class T>
class gauss_hermite_detail<T, 15, 4>
{
   using storage_type = typename gauss_hermite_constant_category<T>::storage_type;
   public:
      static std::array<T, 8> const & abscissa()
      {
         static std::array<T, 8> data = {
            BOOST_MATH_HUGE_CONSTANT(T, 0, 0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000e+00),
            BOOST_MATH_HUGE_CONSTANT(T, 0, 5.6506958325557574852602033719818886681552078800505788680930739682745058470915375159973526283034610292989738379345761e-01),
            BOOST_MATH_HUGE_CONSTANT(T, 0, 1.1361155852109206663191349055561021252835292947846650850648392599800715967101004409631966601830482707422651868052475e+00),
            BOOST_MATH_HUGE_CONSTANT(T, 0, 1.7199925751864889324158315251525615148429195630613888498938216300095132434761832428439718888176122448729491067241320e+00),
            BOOST_MATH_HUGE_CONSTANT(T, 0, 2.3257324861738577454540447944853470123070950095634828957721625948740238048583608727061841554643103134570736658336089e+00),
            BOOST_MATH_HUGE_CONSTANT(T, 0, 2.9671669279056032484889603635498063155705164076381381507212727522118956916046294049934561649283587290854852693064039e+00),
            BOOST_MATH_HUGE_CONSTANT(T, 0, 3.6699503734044525347292238331156814846958556260108268200827250197712929612325341525137596165272567610075980577793224e+00),
            BOOST_MATH_HUGE_CONSTANT(T, 0, 4.4999907073093915536643805305348242199319830534949905938271742846019366228525630586264389556919589962713446071341276e+00),
         };
         return data;
      }
      static std::array<T, 8> const & weights()
      {
         static std::array<T, 8> data = {
            BOOST_MATH_HUGE_CONSTANT(T, 0, 5.6410030872641753285262579733996353292453477640072243023805413420648986959478382400247069852378235816199417789938317e-01),
            BOOST_MATH_HUGE_CONSTANT(T, 0, 4.1202868749889862702589107956781031893375232741692141877844273550979447720862335855236503495494156330750630603931556e-01),
            BOOST_MATH_HUGE_CONSTANT(T, 0, 1.5848891579593574688383938495999402721002703620834420484147417162556491543397074116784210832986285797844312792251258e-01),
            BOOST_MATH_HUGE_CONSTANT(T, 0, 3.0780033872546082228681415875780145569645078095357831670126204320943426296907483543543105075772562113024847790201415e-02),
            BOOST_MATH_HUGE_CONSTANT(T, 0, 2.7780688429127758960788704922921349030124251134822666685717061210098979501372205289327949738880183107056484267658214e-03),
            BOOST_MATH_HUGE_CONSTANT(T, 0, 1.0000444123249986812729673617697850851699387810876387789234107161752152460966503439378917351614005951603805751399265e-04),
            BOOST_MATH_HUGE_CONSTANT(T, 0, 1.0591155477110666357752079105501632872963524018315499158761728455091303265515893684325618169279641299843625141139696e-06),
            BOOST_MATH_HUGE_CONSTANT(T, 0, 1.5224758042535170201606266696482661827962160316014537930014349599535491582697947140375830838937903702636262186650310e-09),
         };
         return data;
      }
   };

} // namespace detail

template <class Real, unsigned N, class Policy = boost::math::policies::policy<> >
class gauss_hermite : public detail::gauss_hermite_detail<Real, N, detail::gauss_hermite_constant_category<Real>::value>
{
   using base = detail::gauss_hermite_detail<Real, N, detail::gauss_hermite_constant_category<Real>::value>;
public:
   // Integrates f(x) exp(-x^2) over the real line.
   template <class F>
   static auto integrate(F f, Real* pL1 = nullptr)->decltype(std::declval<F>()(std::declval<Real>()))
   {
      // In many math texts, K represents the field of real or complex numbers.
      // Too bad we can't put blackboard bold into C++ source!
      using K = decltype(f(Real(0)));
      static_assert(!std::is_integral<K>::value,
                    "The return type cannot be integral, it must be either a real or complex floating point type.");
      using std::abs;
      unsigned non_zero_start = 1;
      K result = Real(0);
      if (N & 1)
      {
         result = f(Real(0)) * static_cast<Real>(base::weights()[0]);
      }
      else
      {
         non_zero_start = 0;
      }
      Real L1 = abs(result);
      for (unsigned i = non_zero_start; i < base::abscissa().size(); ++i)
      {
         K fp = f(static_cast<Real>(base::abscissa()[i]));
         K fm = f(-static_cast<Real>(base::abscissa()[i]));
         result += (fp + fm) * static_cast<Real>(base::weights()[i]);
         L1 += (abs(fp) + abs(fm)) * static_cast<Real>(base::weights()[i]);
      }
      if (pL1)
         *pL1 = L1;
      return result;
   }

   // Explicit zero and norm for vector- and matrix-valued integrands.
   template <class F, class Norm>
   static auto integrate(F f, const decltype(std::declval<F>()(std::declval<Real>()))& zero, Norm norm, Real* pL1 = nullptr)
      ->decltype(static_cast<Real>(norm(f(Real(0)))), f(Real(0)))
   {
      using K = decltype(f(Real(0)));
      static_assert(!std::is_integral<K>::value, "The return type cannot be integral.");
      unsigned non_zero_start = 1;
      K result = zero;
      if (N & 1)
      {
         result = f(Real(0)) * static_cast<Real>(base::weights()[0]);
      }
      else
      {
         non_zero_start = 0;
      }
      Real L1 = static_cast<Real>(norm(result));
      for (unsigned i = non_zero_start; i < base::abscissa().size(); ++i)
      {
         K fp = f(static_cast<Real>(base::abscissa()[i]));
         K fm = f(-static_cast<Real>(base::abscissa()[i]));
         result += (fp + fm) * static_cast<Real>(base::weights()[i]);
         L1 += (static_cast<Real>(norm(fp)) + static_cast<Real>(norm(fm))) * static_cast<Real>(base::weights()[i]);
      }
      if (pL1)
         *pL1 = L1;
      return result;
   }
};

} // namespace quadrature
} // namespace math
} // namespace boost

#endif // BOOST_MATH_QUADRATURE_GAUSS_HERMITE_HPP
