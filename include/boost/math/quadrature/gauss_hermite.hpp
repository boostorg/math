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
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>
#endif
#include <boost/math/constants/constants.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/policies/policy.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/tools/compensated_arithmetic.hpp>
#include <boost/math/tools/big_constant.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/math/quadrature/detail/quadrature_constant.hpp>
#include <boost/math/special_functions/detail/orthogonal_polynomial.hpp>

BOOST_MATH_NAMESPACE_BEGIN namespace detail {

template <class Real>
class hermite_family
{
public:
   hermite_family(unsigned n) : N(n), a_values(N), b_values(N)
   {
      using std::sqrt;
      for (unsigned k = 0; k < N; ++k)
      {
         a_values[k] = sqrt(Real(2) / Real(k + 1));
         b_values[k] = -sqrt(Real(k) / Real(k + 1));
      }
      p0_value = 1 / sqrt(sqrt(BOOST_MATH_NAMESPACE::constants::pi<Real>()));
      derivative_scale = sqrt(Real(2) * Real(N));
   }

   const Real& a(unsigned k) const
   {
      return a_values[k];
   }

   const Real& b(unsigned k) const
   {
      return b_values[k];
   }

   const Real& p0() const
   {
      return p0_value;
   }

   unsigned num_roots() const
   {  // Need to account for odd N having root at 0
      return N / 2 + N % 2;
   }

   Real lower_bound() const
   {
      return 0;
   }

   Real upper_bound() const
   {
      using std::sqrt;
      return sqrt(Real(2) * Real(N) + Real(2));
   }

   Real area() const
   {
      return BOOST_MATH_NAMESPACE::constants::root_pi<Real>();
   }

   unsigned weight_multiplicity(unsigned i) const
   {
      return (N % 2 && i == 0) ? 1 : 2;
   }

   Real next(const Real& x, const Real& p, const Real& p_previous, unsigned k) const
   {
      return a(k) * x * p + b(k) * p_previous;
   }

   Real slope(unsigned k) const
   {  // The coefficient of x p in next()
      return a(k);
   }

   Real derivative(const Real&, const Real&, const Real& p_previous) const
   {
      return derivative_scale * p_previous;
   }

   Real second_derivative(const Real& x, const Real& p, const Real& p_prime) const
   {
      return Real(2) * x * p_prime - Real(2) * Real(N) * p;
   }

   Real second_derivative_ratio(const Real& x, const Real& p, const Real& p_prime) const
   {
      return Real(2) * x - Real(2) * Real(N) * p / p_prime;
   }

   Real weight(const Real& z, const Real& p, const Real& p_prime) const
   {
      Real delta = -p / p_prime;
      Real p_prime_at_root = p_prime + second_derivative(z, p, p_prime) * delta;
      return Real(2) / (p_prime_at_root * p_prime_at_root);
   }

   // For glaser_liu_rokhlin: u = exp(-x^2 / 2) p_N satisfies u'' + (2N + 1 - x^2) u = 0; the coefficients of
   // s(t + h) = 1 and g(t + h) as polynomials in h. Its solutions are entire.
   static constexpr unsigned ends = 0;

   void normal_form(const double_word_arithmetic<Real>& dw, const double_word<Real>& t, double_word<Real> (&s)[5], double_word<Real> (&g)[3], unsigned m) const
   {
      s[0] = double_word<Real>{ Real(1), Real(0) };
      g[0] = dw.add(dw.negate(dw.multiply(t, t)), Real(2) * Real(m) + 1);
      g[1] = double_word<Real>{ -2 * t.hi, -2 * t.lo };
      g[2] = double_word<Real>{ Real(-1), Real(0) };
   }

   Real singular_distance(const Real&) const
   {
      return tools::max_value<Real>();
   }

   Real weight_function(const double_word<Real>& z) const
   {
      using std::exp;
      const double_word<Real> square = double_word_arithmetic<Real>().multiply(z, z);
      return exp(-square.hi) * (1 - square.lo);
   }

   Real sturm_liouville(const double_word<Real>&) const
   {
      return 1;
   }

   // p_N = (A x + B) p_{N-1} - C p_{N-2}: A = sqrt(2 / N), B = 0, and C_2 = -2 C = -2 sqrt((N - 1) / N), in double words.
   void anti_gauss_recurrence(const double_word_arithmetic<Real>& dw, double_word<Real>& A, double_word<Real>& B, double_word<Real>& C2) const
   {
      A = square_root(dw, dw.divide(double_word<Real>{ Real(2), Real(0) }, Real(N)));
      B = double_word<Real>{ Real(0), Real(0) };
      const double_word<Real> C = square_root(dw, dw.divide(double_word<Real>{ Real(N - 1), Real(0) }, Real(N)));
      C2 = double_word<Real>{ -2 * C.hi, -2 * C.lo };
   }

   // One step of Newton's method for the square root, in double words.
   static double_word<Real> square_root(const double_word_arithmetic<Real>& dw, const double_word<Real>& x)
   {
      using std::sqrt;
      const Real root = sqrt(x.hi);
      Real p, e;
      two_product<Real>()(root, root, p, e);
      return dw.add(double_word<Real>{ root, Real(0) }, (((x.hi - p) - e) + x.lo) / (2 * root));
   }

   // For the anti-Gauss polynomial: p_{N-1}(0) / p_{N-2}'(0) = -1 / sqrt(2 (N - 1)) or p_{N-1}'(0) / p_{N-2}(0) = sqrt(2 (N - 1)),
   // from p_k' = sqrt(2 k) p_{k-1} and p_{N-1}(0) = b(N - 2) p_{N-3}(0).
   double_word<Real> anti_gauss_ratio(const double_word_arithmetic<Real>& dw, unsigned) const
   {
      const double_word<Real> r = square_root(dw, double_word<Real>{ 2 * Real(N - 1), Real(0) });
      return (N - 1) % 2 ? r : dw.divide(double_word<Real>{ Real(-1), Real(0) }, r);
   }

private:
   unsigned N;
   std::vector<Real> a_values;
   std::vector<Real> b_values;
   Real p0_value;
   Real derivative_scale;
};

} // namespace detail
BOOST_MATH_NAMESPACE_END

BOOST_MATH_NAMESPACE_BEGIN namespace quadrature { namespace detail {

#ifndef BOOST_MATH_GAUSS_NO_COMPUTE_ON_DEMAND

template <class Real, unsigned N, unsigned Category>
class gauss_hermite_detail
{
   static_assert(N > 0, "Gauss-Hermite quadrature needs at least one point.");

   static const BOOST_MATH_NAMESPACE::detail::orthogonal_polynomial<Real, BOOST_MATH_NAMESPACE::detail::hermite_family<Real>>& polynomial()
   {
      static const BOOST_MATH_NAMESPACE::detail::orthogonal_polynomial<Real, BOOST_MATH_NAMESPACE::detail::hermite_family<Real>> value(N);
      return value;
   }

public:
   static const std::vector<Real>& abscissa()
   {
      return polynomial().abscissa();
   }

   static const std::vector<Real>& weights()
   {
      return polynomial().weights();
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
   using storage_type = typename quadrature_constant_category<T>::storage_type;
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
   using storage_type = typename quadrature_constant_category<T>::storage_type;
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
   using storage_type = typename quadrature_constant_category<T>::storage_type;
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
   using storage_type = typename quadrature_constant_category<T>::storage_type;
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
   using storage_type = typename quadrature_constant_category<T>::storage_type;
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
   using storage_type = typename quadrature_constant_category<T>::storage_type;
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
   using storage_type = typename quadrature_constant_category<T>::storage_type;
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
   using storage_type = typename quadrature_constant_category<T>::storage_type;
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
   using storage_type = typename quadrature_constant_category<T>::storage_type;
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

template <class Real, unsigned N, class Policy = BOOST_MATH_NAMESPACE::policies::policy<> >
class gauss_hermite : public detail::gauss_hermite_detail<Real, N, detail::quadrature_constant_category<Real>::value>
{
   using base = detail::gauss_hermite_detail<Real, N, detail::quadrature_constant_category<Real>::value>;

   // Laurie's anti-Gauss rule with N + 1 points, computed the first time integrate_averaged is called.
   static const BOOST_MATH_NAMESPACE::detail::orthogonal_polynomial<Real, BOOST_MATH_NAMESPACE::detail::hermite_family<Real>>& anti_gauss()
   {
      static const BOOST_MATH_NAMESPACE::detail::orthogonal_polynomial<Real, BOOST_MATH_NAMESPACE::detail::hermite_family<Real>> value(N + 1, BOOST_MATH_NAMESPACE::detail::hermite_family<Real>(N + 1), true);
      return value;
   }

   template <class F, class K, class Norm>
   static K anti_gauss_average(F& f, K zero, const K& result, Real& L1, Norm norm, Real* error)
   {
      return detail::anti_gauss_average(f, anti_gauss().abscissa(), anti_gauss().weights(), true, zero, result, L1, norm, error, "boost::math::quadrature::gauss_hermite<%1%>::integrate_averaged", Policy());
   }

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
      Real nan;
      if (detail::overflowed(base::abscissa(), "boost::math::quadrature::gauss_hermite<%1%>::integrate", pL1, &nan, Policy()))
      {
         return static_cast<K>(nan);
      }
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
      Real nan;
      if (detail::overflowed(base::abscissa(), "boost::math::quadrature::gauss_hermite<%1%>::integrate", pL1, &nan, Policy()))
      {
         K result = zero * nan;
         return result;
      }
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
   // Laurie's averaged rule (G + A) / 2: N + 1 more evaluations of f than integrate, exact for polynomials of degree
   // 2N + 1. The error is half the difference of G and A, a bound on the error when they bracket the integral.
   template <class F>
   static auto integrate_averaged(F f, Real* error = nullptr, Real* pL1 = nullptr)->decltype(std::declval<F>()(std::declval<Real>()))
   {
      using K = decltype(f(Real(0)));
      Real L1;
      K result = integrate(f, &L1);
      result = anti_gauss_average(f, K(Real(0)), result, L1, [](const K& z) { using std::abs; return abs(z); }, error);
      if (pL1)
         *pL1 = L1;
      return result;
   }

   template <class F, class Norm>
   static auto integrate_averaged(F f, const decltype(std::declval<F>()(std::declval<Real>()))& zero, Norm norm, Real* error = nullptr, Real* pL1 = nullptr)
      ->decltype(static_cast<Real>(norm(f(Real(0)))), f(Real(0)))
   {
      Real L1;
      auto result = integrate(f, zero, norm, &L1);
      result = anti_gauss_average(f, zero, result, L1, norm, error);
      if (pL1)
         *pL1 = L1;
      return result;
   }
};

} // namespace quadrature
BOOST_MATH_NAMESPACE_END

#endif // BOOST_MATH_QUADRATURE_GAUSS_HERMITE_HPP
