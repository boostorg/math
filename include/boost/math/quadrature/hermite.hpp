//  Copyright Jacob Hass 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_QUADRATURE_HERMITE_HPP
#define BOOST_MATH_QUADRATURE_HERMITE_HPP

#include <array>
#include <cmath>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>
#include <boost/math/policies/policy.hpp>
#include <boost/math/tools/big_constant.hpp>
#include <iostream>

#if __has_include(<Eigen/Dense>)
   #include <Eigen/Dense>
   #include <Eigen/Eigenvalues>
   #include <boost/multiprecision/eigen.hpp>
   #include <boost/multiprecision/cpp_bin_float.hpp>
   #include <boost/math/constants/constants.hpp>
   #include <boost/math/special_functions/hermite.hpp>

   #define EIGEN_SUPPORT
#endif

namespace boost { namespace math{ namespace quadrature{ namespace detail {


#if defined(EIGEN_SUPPORT) && !defined(BOOST_MATH_GAUSS_NO_COMPUTE_ON_DEMAND)

template <typename Real>
Real factorial(unsigned n)
{
   Real i = 1;
   Real factorial = 1;

   while (i <= n)
   {
      factorial *= i;
      i++;
   }
   return factorial;
}

template <class Real, unsigned N, unsigned Category>
class hermite_detail
{
public:
   static std::pair<std::vector<Real>, std::vector<Real> > calculate_values()
   {
      namespace mp = boost::multiprecision;

      Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> P = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>::Zero(N, N);
      for (unsigned n = 0; n < N; ++n)
      {
         P(n, n) = Real(0);

         if (n > 0)
         {
               P(n, n-1) = sqrt(Real(n) / Real(2));
               P(n-1, n) = sqrt(Real(n) / Real(2));
         }
      }

      Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> > solver(P);
      Eigen::Vector<Real, Eigen::Dynamic> roots = solver.eigenvalues();

      // Filter out negative roots
      std::vector<Real> pos_roots;
      pos_roots.reserve(N);
      for (unsigned i=0; i < roots.size(); i++)
      {
        // Need to account for root near 0
        if (roots(i) >= -std::numeric_limits<Real>::epsilon() * 1000)
        {
            pos_roots.push_back(roots(i));
        }
      }

      std::vector<Real> abscissa_vals(pos_roots.size());
      std::vector<Real> weight_vals(pos_roots.size());

      for (unsigned i=0; i < pos_roots.size(); i++)
      {
         abscissa_vals[i] = pos_roots[i];
         Real hermite_val = boost::math::hermite<Real>(N-1, pos_roots[i]);
         Real weight = mp::pow(Real(2), N-1) * factorial<Real>(N) * boost::math::constants::root_pi<Real>() / mp::pow(Real(N), 2) / mp::pow(hermite_val, 2.0);
         weight_vals[i] = weight;
      }

      // If N is odd, first abscissa is 0
      if (N % 2 != 0)
      {
         abscissa_vals[0] = Real(0);
      }

      return std::make_pair(abscissa_vals, weight_vals);
   }

   static const std::vector<Real>& weights()
   {
      static std::pair<std::vector<Real>, std::vector<Real> > data = calculate_values();
      return data.second;
   }

   static const std::vector<Real>& abscissa()
   {
      static std::pair<std::vector<Real>, std::vector<Real> > data = calculate_values();
      return data.first;
   }
};

#else

template <class Real, unsigned N, unsigned Category>
class hermite_detail;

#endif

template <class T>
struct hermite_constant_category
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

#ifndef BOOST_HAS_FLOAT128
template <class T>
class hermite_detail<T, 7, 0>
{
   using storage_type = typename hermite_constant_category<T>::storage_type;
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
class hermite_detail<T, 7, 0>
{
   using storage_type = typename hermite_constant_category<T>::storage_type;
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
class hermite_detail<T, 7, 4>
{
   using storage_type = typename hermite_constant_category<T>::storage_type;
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
class hermite_detail<T, 10, 0>
{
   using storage_type = typename hermite_constant_category<T>::storage_type;
   public:
      static std::array<storage_type, 5> const & abscissa()
      {
         static std::array<storage_type, 5> data = {
            static_cast<storage_type>(3.42901327223704608789165025557258046e-01L),
            static_cast<storage_type>(1.03661082978951365417749191675920910e+00L),
            static_cast<storage_type>(1.75668364929988177345140122010615672e+00L),
            static_cast<storage_type>(2.53273167423278979640896079775479347e+00L),
            static_cast<storage_type>(3.43615911883773760332672549431912143e+00L),
};
         return data;
      }
      static std::array<storage_type, 5> const & weights()
      {
         static std::array<storage_type, 5> data = {
            static_cast<storage_type>(6.10862633735325798783564990433419732e-01L),
            static_cast<storage_type>(2.40138611082314686416523295005861392e-01L),
            static_cast<storage_type>(3.38743944554810631361647312775859719e-02L),
            static_cast<storage_type>(1.34364574678123269220156558584591379e-03L),
            static_cast<storage_type>(7.64043285523262062915936785959522150e-06L),
         };
         return data;
      }
   };

#else
template <class T>
class hermite_detail<T, 10, 0>
{
   using storage_type = typename hermite_constant_category<T>::storage_type;
   public:
      static std::array<storage_type, 5> const & abscissa()
      {
         static std::array<storage_type, 5> data = {
            static_cast<storage_type>(3.42901327223704608789165025557258046e-01Q),
            static_cast<storage_type>(1.03661082978951365417749191675920910e+00Q),
            static_cast<storage_type>(1.75668364929988177345140122010615672e+00Q),
            static_cast<storage_type>(2.53273167423278979640896079775479347e+00Q),
            static_cast<storage_type>(3.43615911883773760332672549431912143e+00Q),
};
         return data;
      }
      static std::array<storage_type, 5> const & weights()
      {
         static std::array<storage_type, 5> data = {
            static_cast<storage_type>(6.10862633735325798783564990433419732e-01Q),
            static_cast<storage_type>(2.40138611082314686416523295005861392e-01Q),
            static_cast<storage_type>(3.38743944554810631361647312775859719e-02Q),
            static_cast<storage_type>(1.34364574678123269220156558584591379e-03Q),
            static_cast<storage_type>(7.64043285523262062915936785959522150e-06Q),
         };
         return data;
      }
   };

#endif
template <class T>
class hermite_detail<T, 10, 4>
{
   using storage_type = typename hermite_constant_category<T>::storage_type;
   public:
      static std::array<T, 5> const & abscissa()
      {
         static std::array<T, 5> data = {
            BOOST_MATH_HUGE_CONSTANT(T, 0, 3.4290132722370460878916502555725804574577169343357617326210300863769016299631431365924072451889514923095703125000000e-01),
            BOOST_MATH_HUGE_CONSTANT(T, 0, 1.0366108297895136541774919167592091040917543343984240793575756300019937924949964269671909278258681297302246093750000e+00),
            BOOST_MATH_HUGE_CONSTANT(T, 0, 1.7566836492998817734514012201061567205847985641023178720560650082445006137099596799089340493083000183105468750000000e+00),
            BOOST_MATH_HUGE_CONSTANT(T, 0, 2.5327316742327897964089607977547934701687784149754958919647942962061508073712268185317952884361147880554199218750000e+00),
            BOOST_MATH_HUGE_CONSTANT(T, 0, 3.4361591188377376033267254943191214342191179952800433354144214837346850421884170145858661271631717681884765625000000e+00),
};
         return data;
      }
      static std::array<T, 5> const & weights()
      {
         static std::array<T, 5> data = {
            BOOST_MATH_HUGE_CONSTANT(T, 0, 6.1086263373532579878356499043341973238818374059758277936763908978968140237042483420282223960384726524353027343750000e-01),
            BOOST_MATH_HUGE_CONSTANT(T, 0, 2.4013861108231468641652329500586139164344096899966477839326207760887354125331483167826718272408470511436462402343750e-01),
            BOOST_MATH_HUGE_CONSTANT(T, 0, 3.3874394455481063136164731277585971926477457178246590999848008213597529213016723570461863346281461417675018310546875e-02),
            BOOST_MATH_HUGE_CONSTANT(T, 0, 1.3436457467812326922015655858459137918851865251396942830563543358449043010113127527560550333873834460973739624023438e-03),
            BOOST_MATH_HUGE_CONSTANT(T, 0, 7.6404328552326206291593678595952215021994824892912048011683019177853182368166190746930355182087168941507115960121155e-06),
         };
         return data;
      }
   };

#ifndef BOOST_HAS_FLOAT128
template <class T>
class hermite_detail<T, 15, 0>
{
   using storage_type = typename hermite_constant_category<T>::storage_type;
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
class hermite_detail<T, 15, 0>
{
   using storage_type = typename hermite_constant_category<T>::storage_type;
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
class hermite_detail<T, 15, 4>
{
   using storage_type = typename hermite_constant_category<T>::storage_type;
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
class hermite : public detail::hermite_detail<Real, N, detail::hermite_constant_category<Real>::value>
{
   typedef detail::hermite_detail<Real, N, detail::hermite_constant_category<Real>::value> base;
public:

   template <class F>
   static auto integrate(F f, Real* pL1 = nullptr)->decltype(std::declval<F>()(std::declval<Real>()))
   {
     // In many math texts, K represents the field of real or complex numbers.
     // Too bad we can't put blackboard bold into C++ source!
      typedef decltype(f(Real(0))) K;
      static_assert(!std::is_integral<K>::value,
                   "The return type cannot be integral, it must be either a real or complex floating point type.");
      using std::abs;
      unsigned non_zero_start = 1;
      K result = Real(0);
      if (N & 1) {
         result = f(Real(0)) * static_cast<Real>(base::weights()[0]);
      }
      else {
         result = 0;
         non_zero_start = 0;
      }
      Real L1 = abs(result);
      Real weight_total;
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
};

} // namespace quadrature
} // namespace math
} // namespace boost

#endif // BOOST_MATH_QUADRATURE_HERMITE_HPP