//  Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SF_GAMMA_HPP
#define BOOST_MATH_SF_GAMMA_HPP

#include <boost/config.hpp>
#ifdef BOOST_MSVC
# pragma warning(push)
# pragma warning(disable: 4127 4701)
//  // For lexical_cast, until fixed in 1.35?
//  // conditional expression is constant &
//  // Potentially uninitialized local variable 'name' used
#endif
#include <boost/lexical_cast.hpp>
#ifdef BOOST_MSVC
# pragma warning(pop)
#endif
#include <boost/math/tools/series.hpp>
#include <boost/math/tools/fraction.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/math/tools/real_cast.hpp>
#include <boost/math/tools/error_handling.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/math_fwd.hpp>
#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/special_functions/powm1.hpp>
#include <boost/math/special_functions/sqrt1pm1.hpp>
#include <boost/math/special_functions/lanczos.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/special_functions/detail/igamma_large.hpp>
#include <boost/math/special_functions/detail/unchecked_factorial.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/assert.hpp>

#include <cmath>
#include <algorithm>

#ifdef BOOST_MATH_INSTRUMENT
#include <iostream>
#include <iomanip>
#include <typeinfo>
#endif

#ifdef BOOST_MSVC
# pragma warning(push)
# pragma warning(disable: 4702) // unreachable code (return after domain_error throw).
# pragma warning(disable: 4127) // conditional expression is constant.
#endif

namespace boost{ namespace math{

namespace detail{

template <class T>
inline bool is_odd(T v, const boost::true_type&)
{
   int i = static_cast<int>(v);
   return i&1;
}
template <class T>
inline bool is_odd(T v, const boost::false_type&)
{
   // Oh dear can't cast T to int!
   using namespace std;
   T modulus = v - 2 * floor(v/2);
   return static_cast<bool>(modulus != 0);
}
template <class T>
inline bool is_odd(T v)
{
   return is_odd(v, ::boost::is_convertible<T, int>());
}

template <class T>
T sinpx(T z)
{
   // Ad hoc function calculates x * sin(pi * x),
   // taking extra care near when x is near a whole number.
   using namespace std;
   int sign = 1;
   if(z < 0)
   {
      z = -z;
   }
   else
   {
      sign = -sign;
   }
   T fl = floor(z);
   T dist;
   if(is_odd(fl))
   {
      fl += 1;
      dist = fl - z;
      sign = -sign;
   }
   else
   {
      dist = z - fl;
   }
   BOOST_ASSERT(fl >= 0);
   if(dist > 0.5)
      dist = 1 - dist;
   T result = sin(dist*boost::math::constants::pi<T>());
   return sign*z*result;
} // template <class T> T sinpx(T z)
//
// tgamma(z), with Lanczos support:
//
template <class T, class L>
T gamma_imp(T z, const L& l)
{
   using namespace std;

   T result = 1;

#ifdef BOOST_MATH_INSTRUMENT
   static bool b = false;
   if(!b)
   {
      std::cout << "tgamma_imp called with " << typeid(z).name() << " " << typeid(l).name() << std::endl;
      b = true;
   }
#endif

   if((z <= 0) && (floor(z) == z))
      return tools::pole_error<T>(BOOST_CURRENT_FUNCTION, "Evaluation of tgamma at a negative integer %1%.", z);
   if(z <= -20)
   {
      result = gamma_imp(-z, l) * sinpx(z);
      if((fabs(result) < 1) && (tools::max_value<T>() * fabs(result) < boost::math::constants::pi<T>()))
         return tools::overflow_error<T>(BOOST_CURRENT_FUNCTION, "Result of tgamma is too large to represent.");
      result = -boost::math::constants::pi<T>() / result;
      if(result == 0)
         return tools::underflow_error<T>(BOOST_CURRENT_FUNCTION, "Result of tgamma is too small to represent.");
      if(boost::math::fpclassify(result) == FP_SUBNORMAL)
         return tools::denorm_error<T>(result, BOOST_CURRENT_FUNCTION, "Result of tgamma is denormalized.");
      return result;
   }

   // shift z to > 1:
   while(z < 1)
   {
      result /= z;
      z += 1;
   }
   if((floor(z) == z) && (z < max_factorial<T>::value))
   {
      result *= unchecked_factorial<T>(tools::real_cast<unsigned>(z) - 1);
   }
   else
   {
      result *= L::lanczos_sum(z);
      if(z * log(z) > tools::log_max_value<T>())
      {
         // we're going to overflow unless this is done with care:
         T zgh = (z + L::g() - boost::math::constants::half<T>());
         T hp = pow(zgh, (z / 2) - T(0.25));
         result *= hp / exp(zgh);
         if(tools::max_value<T>() / hp < result)
            return tools::overflow_error<T>(BOOST_CURRENT_FUNCTION, "Result of tgamma is too large to represent.");
         result *= hp;
      }
      else
      {
         T zgh = (z + L::g() - boost::math::constants::half<T>());
         result *= pow(zgh, z - boost::math::constants::half<T>()) / exp(zgh);
      }
   }
   return result;
}
//
// lgamma for small arguments:
//
template <class T, class L>
T lgamma_small_imp(T z, T zm1, T zm2, const mpl::int_<64>&, const L& /* l */)
{
   // This version uses rational approximations for small
   // values of z accurate enough for 64-bit mantissas
   // (80-bit long doubles), works well for 53-bit doubles as well.
   // L is only used to select the Lanczos function.

   using namespace std;  // for ADL of std names
   T result = 0;
   if(z < tools::epsilon<T>())
   {
      result = -log(z);
   }
   else if((zm1 == 0) || (zm2 == 0))
   {
      // nothing to do, result is zero....
   }
   else if(z > 2)
   {
      //
      // Begin by performing argument reduction until
      // z is in [2,3):
      //
      if(z >= 3)
      {
         do
         {
            z -= 1;
            zm2 -= 1;
            result += log(z);
         }while(z >= 3);
         // Update zm2, we need it below:
         zm2 = z - 2;
      }

      //
      // Use the following form:
      //
      // lgamma(z) = (z-2)(z+1)(Y + R(z-2))
      //
      // where R(z-2) is a rational approximation optimised for
      // low absolute error - as long as it's absolute error
      // is small compared to the constant Y - then any rounding
      // error in it's computation will get wiped out.
      //
      // R(z-2) has the following properties:
      //
      // At double: Max error found:                    4.231e-18
      // At long double: Max error found:               1.987e-21
      // Maximum Deviation Found (approximation error): 5.900e-24
      //
      static const T P[] = {
         -0.180355685678449379109e-1L,
         0.25126649619989678683e-1L,
         0.494103151567532234274e-1L,
         0.172491608709613993966e-1L,
         -0.259453563205438108893e-3L,
         -0.541009869215204396339e-3L,
         -0.324588649825948492091e-4L
      };
      static const T Q[] = {
         0.1e1,
         0.196202987197795200688e1L,
         0.148019669424231326694e1L,
         0.541391432071720958364e0L,
         0.988504251128010129477e-1L,
         0.82130967464889339326e-2L,
         0.224936291922115757597e-3L,
         -0.223352763208617092964e-6L
      };

      static const float Y = 0.158963680267333984375e0f;

      T r = zm2 * (z + 1);
      T R = tools::evaluate_polynomial(P, zm2);
      R /= tools::evaluate_polynomial(Q, zm2);

      result +=  r * Y + r * R;
   }
   else
   {
      //
      // If z is less than 1 use recurrance to shift to
      // z in the interval [1,2]:
      //
      if(z < 1)
      {
         result += -log(z);
         zm2 = zm1;
         zm1 = z;
         z += 1;
      }
      //
      // Two approximations, on for z in [1,1.5] and
      // one for z in [1.5,2]:
      //
      if(z <= 1.5)
      {
         //
         // Use the following form:
         //
         // lgamma(z) = (z-1)(z-2)(Y + R(z-1))
         //
         // where R(z-1) is a rational approximation optimised for
         // low absolute error - as long as it's absolute error
         // is small compared to the constant Y - then any rounding
         // error in it's computation will get wiped out.
         //
         // R(z-1) has the following properties:
         //
         // At double precision: Max error found:                1.230011e-17
         // At 80-bit long double precision:   Max error found:  5.631355e-21
         // Maximum Deviation Found:                             3.139e-021
         // Expected Error Term:                                 3.139e-021

         //
         static const float Y = 0.52815341949462890625f;

         static const T P[] = {
            0.490622454069039543534e-1L,
            -0.969117530159521214579e-1L,
            -0.414983358359495381969e0L,
            -0.406567124211938417342e0L,
            -0.158413586390692192217e0L,
            -0.240149820648571559892e-1L,
            -0.100346687696279557415e-2L
         };
         static const T Q[] = {
            0.1e1L,
            0.302349829846463038743e1L,
            0.348739585360723852576e1L,
            0.191415588274426679201e1L,
            0.507137738614363510846e0L,
            0.577039722690451849648e-1L,
            0.195768102601107189171e-2L
         };

         T r = tools::evaluate_polynomial(P, zm1) / tools::evaluate_polynomial(Q, zm1);
         T prefix = zm1 * zm2;

         result += prefix * Y + prefix * r;
      }
      else
      {
         //
         // Use the following form:
         //
         // lgamma(z) = (2-z)(1-z)(Y + R(2-z))
         //
         // where R(2-z) is a rational approximation optimised for
         // low absolute error - as long as it's absolute error
         // is small compared to the constant Y - then any rounding
         // error in it's computation will get wiped out.
         //
         // R(2-z) has the following properties:
         //
         // At double precision, max error found:              1.797565e-17
         // At 80-bit long double precision, max error found:  9.306419e-21
         // Maximum Deviation Found:                           2.151e-021
         // Expected Error Term:                               2.150e-021
         //
         static const float Y = 0.452017307281494140625f;

         static const T P[] = {
            -0.292329721830270012337e-1L, 
            0.144216267757192309184e0L,
            -0.142440390738631274135e0L,
            0.542809694055053558157e-1L,
            -0.850535976868336437746e-2L,
            0.431171342679297331241e-3L
         };
         static const T Q[] = {
            0.1e1,
            -0.150169356054485044494e1L,
            0.846973248876495016101e0L,
            -0.220095151814995745555e0L,
            0.25582797155975869989e-1L,
            -0.100666795539143372762e-2L,
            -0.827193521891290553639e-6L
         };
         T r = zm2 * zm1;
         T R = tools::evaluate_polynomial(P, -zm2) / tools::evaluate_polynomial(Q, -zm2);

         result += r * Y + r * R;
      }
   }
   return result;
}
template <class T, class L>
T lgamma_small_imp(T z, T zm1, T zm2, const mpl::int_<113>&, const L& /* l */)
{
   //
   // This version uses rational approximations for small
   // values of z accurate enough for 113-bit mantissas
   // (128-bit long doubles).
   //
   using namespace std;  // for ADL of std names
   T result = 0;
   if(z < tools::epsilon<T>())
   {
      result = -log(z);
   }
   else if((zm1 == 0) || (zm2 == 0))
   {
      // nothing to do, result is zero....
   }
   else if(z > 2)
   {
      //
      // Begin by performing argument reduction until
      // z is in [2,3):
      //
      if(z >= 3)
      {
         do
         {
            z -= 1;
            result += log(z);
         }while(z >= 3);
         zm2 = z - 2;
      }

      //
      // Use the following form:
      //
      // lgamma(z) = (z-2)(z+1)(Y + R(z-2))
      //
      // where R(z-2) is a rational approximation optimised for
      // low absolute error - as long as it's absolute error
      // is small compared to the constant Y - then any rounding
      // error in it's computation will get wiped out.
      //
      // Maximum Deviation Found (approximation error)      3.73e-37

      static const T P[] = {
         -0.018035568567844937910504030027467476655L,
         0.013841458273109517271750705401202404195L,
         0.062031842739486600078866923383017722399L,
         0.052518418329052161202007865149435256093L,
         0.01881718142472784129191838493267755758L,
         0.0025104830367021839316463675028524702846L,
         -0.00021043176101831873281848891452678568311L,
         -0.00010249622350908722793327719494037981166L,
         -0.11381479670982006841716879074288176994e-4L,
         -0.49999811718089980992888533630523892389e-6L,
         -0.70529798686542184668416911331718963364e-8L
      };
      static const T Q[] = {
         1L,
         2.5877485070422317542808137697939233685L,
         2.8797959228352591788629602533153837126L,
         1.8030885955284082026405495275461180977L,
         0.69774331297747390169238306148355428436L,
         0.17261566063277623942044077039756583802L,
         0.02729301254544230229429621192443000121L,
         0.0026776425891195270663133581960016620433L,
         0.00015244249160486584591370355730402168106L,
         0.43997034032479866020546814475414346627e-5L,
         0.46295080708455613044541885534408170934e-7L,
         -0.93326638207459533682980757982834180952e-11L,
         0.42316456553164995177177407325292867513e-13L
      };

      T R = tools::evaluate_polynomial(P, zm2);
      R /= tools::evaluate_polynomial(Q, zm2);

      static const float Y = 0.158963680267333984375F;

      T r = zm2 * (z + 1);

      result +=  r * Y + r * R;
   }
   else
   {
      //
      // If z is less than 1 use recurrance to shift to
      // z in the interval [1,2]:
      //
      if(z < 1)
      {
         result += -log(z);
         zm2 = zm1;
         zm1 = z;
         z += 1;
      }
      //
      // Three approximations, on for z in [1,1.35], [1.35,1.625] and [1.625,1]
      //
      if(z <= 1.35)
      {
         //
         // Use the following form:
         //
         // lgamma(z) = (z-1)(z-2)(Y + R(z-1))
         //
         // where R(z-1) is a rational approximation optimised for
         // low absolute error - as long as it's absolute error
         // is small compared to the constant Y - then any rounding
         // error in it's computation will get wiped out.
         //
         // R(z-1) has the following properties:
         //
         // Maximum Deviation Found (approximation error)            1.659e-36
         // Expected Error Term (theoretical error)                  1.343e-36
         // Max error found at 128-bit long double precision         1.007e-35
         //
         static const float Y = 0.54076099395751953125f;

         static const T P[] = {
            0.036454670944013329356512090082402429697L,
            -0.066235835556476033710068679907798799959L,
            -0.67492399795577182387312206593595565371L,
            -1.4345555263962411429855341651960000166L,
            -1.4894319559821365820516771951249649563L,
            -0.87210277668067964629483299712322411566L,
            -0.29602090537771744401524080430529369136L,
            -0.0561832587517836908929331992218879676L,
            -0.0053236785487328044334381502530383140443L,
            -0.00018629360291358130461736386077971890789L,
            -0.10164985672213178500790406939467614498e-6L,
            0.13680157145361387405588201461036338274e-8L
         };
         static const T Q[] = {
            1,
            4.9106336261005990534095838574132225599L,
            10.258804800866438510889341082793078432L,
            11.88588976846826108836629960537466889L,
            8.3455000546999704314454891036700998428L,
            3.6428823682421746343233362007194282703L,
            0.97465989807254572142266753052776132252L,
            0.15121052897097822172763084966793352524L,
            0.012017363555383555123769849654484594893L,
            0.0003583032812720649835431669893011257277L
         };

         T r = tools::evaluate_polynomial(P, zm1) / tools::evaluate_polynomial(Q, zm1);
         T prefix = zm1 * zm2;

         result += prefix * Y + prefix * r;
      }
      else if(z <= 1.625)
      {
         //
         // Use the following form:
         //
         // lgamma(z) = (2-z)(1-z)(Y + R(2-z))
         //
         // where R(2-z) is a rational approximation optimised for
         // low absolute error - as long as it's absolute error
         // is small compared to the constant Y - then any rounding
         // error in it's computation will get wiped out.
         //
         // R(2-z) has the following properties:
         //
         // Max error found at 128-bit long double precision  9.634e-36
         // Maximum Deviation Found (approximation error)     1.538e-37
         // Expected Error Term (theoretical error)           2.350e-38
         //
         static const float Y = 0.483787059783935546875f;

         static const T P[] = {
            -0.017977422421608624353488126610933005432L,
            0.18484528905298309555089509029244135703L,
            -0.40401251514859546989565001431430884082L,
            0.40277179799147356461954182877921388182L,
            -0.21993421441282936476709677700477598816L,
            0.069595742223850248095697771331107571011L,
            -0.012681481427699686635516772923547347328L,
            0.0012489322866834830413292771335113136034L,
            -0.57058739515423112045108068834668269608e-4L,
            0.8207548771933585614380644961342925976e-6L
         };
         static const T Q[] = {
            1,
            -2.9629552288944259229543137757200262073L,
            3.7118380799042118987185957298964772755L,
            -2.5569815272165399297600586376727357187L,
            1.0546764918220835097855665680632153367L,
            -0.26574021300894401276478730940980810831L,
            0.03996289731752081380552901986471233462L,
            -0.0033398680924544836817826046380586480873L,
            0.00013288854760548251757651556792598235735L,
            -0.17194794958274081373243161848194745111e-5L
         };
         T r = zm2 * zm1;
         T R = tools::evaluate_polynomial(P, 0.625 - zm1) / tools::evaluate_polynomial(Q, 0.625 - zm1);

         result += r * Y + r * R;
      }
      else
      {
         //
         // Same form as above.
         //
         // Max error found (at 128-bit long double precision) 1.831e-35
         // Maximum Deviation Found (approximation error)      8.588e-36
         // Expected Error Term (theoretical error)            1.458e-36
         //
         static const float Y = 0.443811893463134765625f;

         static const T P[] = {
            -0.021027558364667626231512090082402429494L,
            0.15128811104498736604523586803722368377L,
            -0.26249631480066246699388544451126410278L,
            0.21148748610533489823742352180628489742L,
            -0.093964130697489071999873506148104370633L,
            0.024292059227009051652542804957550866827L,
            -0.0036284453226534839926304745756906117066L,
            0.0002939230129315195346843036254392485984L,
            -0.11088589183158123733132268042570710338e-4L,
            0.13240510580220763969511741896361984162e-6L
         };
         static const T Q[] = {
            1,
            -2.4240003754444040525462170802796471996L,
            2.4868383476933178722203278602342786002L,
            -1.4047068395206343375520721509193698547L,
            0.47583809087867443858344765659065773369L,
            -0.09865724264554556400463655444270700132L,
            0.012238223514176587501074150988445109735L,
            -0.00084625068418239194670614419707491797097L,
            0.2796574430456237061420839429225710602e-4L,
            -0.30202973883316730694433702165188835331e-6L
         };
         // (2 - x) * (1 - x) * (c + R(2 - x))
         T r = zm2 * zm1;
         T R = tools::evaluate_polynomial(P, -zm2) / tools::evaluate_polynomial(Q, -zm2);

         result += r * Y + r * R;
      }
   }
   return result;
}
template <class T, class L>
T lgamma_small_imp(T z, T zm1, T zm2, const mpl::int_<0>&, const L& l)
{
   //
   // No rational approximations are available because either
   // T has no numeric_limits support (so we can't tell how
   // many digits it has), or T has more digits than we know
   // what to do with.... we do have a Lanczos approximation
   // though, and that can be used to keep errors under control.
   //
   using namespace std;  // for ADL of std names
   T result = 0;
   if(z < tools::epsilon<T>())
   {
      result = -log(z);
   }
   else if(z < 0.5)
   {
      // taking the log of tgamma reduces the error, no danger of overflow here:
      result = log(gamma_imp(z, l));
   }
   else if(z >= 3)
   {
      // taking the log of tgamma reduces the error, no danger of overflow here:
      result = log(gamma_imp(z, l));
   }
   else if(z >= 1.5)
   {
      // special case near 2:
      T dz = zm2;
      result = dz * log((z + L::g() - T(0.5)) / boost::math::constants::e<T>());
      result += boost::math::log1p(dz / (L::g() + T(1.5))) * T(1.5);
      result += boost::math::log1p(L::lanczos_sum_near_2(dz));
   }
   else
   {
      // special case near 1:
      T dz = zm1;
      result = dz * log((z + L::g() - T(0.5)) / boost::math::constants::e<T>());
      result += boost::math::log1p(dz / (L::g() + T(0.5))) / 2;
      result += boost::math::log1p(L::lanczos_sum_near_1(dz));
   }
   return result;
}
//
// lgamma(z) with Lanczos support:
//
template <class T, class L>
T lgamma_imp(T z, const L& l, int* sign = 0)
{
#ifdef BOOST_MATH_INSTRUMENT
   static bool b = false;
   if(!b)
   {
      std::cout << "lgamma_imp called with " << typeid(z).name() << " " << typeid(l).name() << std::endl;
      b = true;
   }
#endif

   using namespace std;

   T result = 0;
   int sresult = 1;
   if(z <= 0)
   {
      // reflection formula:
      if(floor(z) == z)
         return tools::pole_error<T>(BOOST_CURRENT_FUNCTION, "Evaluation of lgamma at a negative integer %1%.", z);

      T t = sinpx(z);
      z = -z;
      if(t < 0)
      {
         t = -t;
      }
      else
      {
         sresult = -sresult;
      }
      result = log(boost::math::constants::pi<T>()) - lgamma_imp(z, l) - log(t);
   }
   else if(z < 15)
   {
      typedef typename mpl::if_c<
         (std::numeric_limits<T>::is_specialized && (std::numeric_limits<T>::digits <= 64)),
         mpl::int_<64>,
         typename mpl::if_c<
            (std::numeric_limits<T>::is_specialized && (std::numeric_limits<T>::digits <= 113)),
            mpl::int_<113>, mpl::int_<0> >::type
          >::type tag_type;
      result = lgamma_small_imp(z, z - 1, z - 2, tag_type(), l);
   }
   else if((z >= 3) && (z < 100))
   {
      // taking the log of tgamma reduces the error, no danger of overflow here:
      result = log(gamma_imp(z, l));
   }
   else
   {
      // regular evaluation:
      T zgh = (z + L::g() - boost::math::constants::half<T>());
      T l = L::lanczos_sum_expG_scaled(z);
      result = log(zgh) - 1;
      result *= z - 0.5f;
      result += log(l);
   }

   if(sign)
      *sign = sresult;
   return result;
}

//
// Incomplete gamma functions follow:
//
template <class T>
struct upper_incomplete_gamma_fract
{
private:
   T z, a;
   int k;
public:
   typedef std::pair<T,T> result_type;

   upper_incomplete_gamma_fract(T a1, T z1)
      : z(z1-a1+1), a(a1), k(0)
   {
   }

   result_type operator()()
   {
      ++k;
      z += 2;
      return result_type(k * (a - k), z);
   }
};

template <class T>
T upper_gamma_fraction(T a, T z, int bits)
{
   // Multiply result by z^a * e^-z to get the full
   // upper incomplete integral.  Divide by tgamma(z)
   // to normalise.
   upper_incomplete_gamma_fract<T> f(a, z);
   return 1 / (z - a + 1 + boost::math::tools::continued_fraction_a(f, bits));
}

template <class T>
struct lower_incomplete_gamma_series
{
private:
   T a, z, result;
public:
   typedef T result_type;
   lower_incomplete_gamma_series(T a1, T z1) : a(a1), z(z1), result(1){}

   T operator()()
   {
      T r = result;
      a += 1;
      result *= z/a;
      return r;
   }
};

template <class T>
T lower_gamma_series(T a, T z, int bits)
{
   // Multiply result by ((z^a) * (e^-z) / a) to get the full
   // lower incomplete integral. Then divide by tgamma(a)
   // to get the normalised value.
   lower_incomplete_gamma_series<T> s(a, z);
   boost::uintmax_t max_iter = BOOST_MATH_MAX_ITER;
#if BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x582))
   T zero = 0;
   T result = boost::math::tools::sum_series(s, bits, max_iter, zero);
#else
   T result = boost::math::tools::sum_series(s, bits, max_iter);
#endif
   tools::check_series_iterations(BOOST_CURRENT_FUNCTION, max_iter);
   return result;
}

//
// Fully generic tgamma and lgamma use the incomplete partial
// sums added together:
//
template <class T>
T gamma_imp(T z, const lanczos::undefined_lanczos& l)
{
   using namespace std;
   if((z <= 0) && (floor(z) == z))
      return tools::pole_error<T>(BOOST_CURRENT_FUNCTION, "Evaluation of tgamma at a negative integer %1%.", z);
   if(z <= -20)
   {
      T result = gamma_imp(-z, l) * sinpx(z);
      if((fabs(result) < 1) && (tools::max_value<T>() * fabs(result) < boost::math::constants::pi<T>()))
         return tools::overflow_error<T>(BOOST_CURRENT_FUNCTION, "Result of tgamma is too large to represent.");
      result = -boost::math::constants::pi<T>() / result;
      if(result == 0)
         return tools::underflow_error<T>(BOOST_CURRENT_FUNCTION, "Result of tgamma is too small to represent.");
      if(boost::math::fpclassify(result) == FP_SUBNORMAL)
         return tools::denorm_error<T>(result, BOOST_CURRENT_FUNCTION, "Result of tgamma is denormalized.");
      return result;
   }
   //
   // The upper gamma fraction is *very* slow for z < 6, actually it's very
   // slow to converge everywhere but recursing until z > 6 gets rid of the
   // worst of it's behaviour.
   //
   T prefix = 1;
   while(z < 6)
   {
      prefix /= z;
      z += 1;
   }
   if((floor(z) == z) && (z < max_factorial<T>::value))
   {
      prefix *= unchecked_factorial<T>(tools::real_cast<unsigned>(z) - 1);
   }
   else
   {
      prefix = prefix * pow(z / boost::math::constants::e<T>(), z);
      T sum = detail::lower_gamma_series(z, z, ::boost::math::tools::digits<T>()) / z;
      sum += detail::upper_gamma_fraction(z, z, ::boost::math::tools::digits<T>());
      if(fabs(tools::max_value<T>() / prefix) < fabs(sum))
         return tools::overflow_error<T>(BOOST_CURRENT_FUNCTION, "Result of tgamma is too large to represent.");
      return sum * prefix;
   }
   return prefix;
}

template <class T>
T lgamma_imp(T z, const lanczos::undefined_lanczos&, int*sign)
{
   using namespace std;

   T result = 0;
   int sresult = 1;
   if(z <= 0)
   {
      if(floor(z) == z)
         return tools::pole_error<T>(BOOST_CURRENT_FUNCTION, "Evaluation of tgamma at a negative integer %1%.", z);
      T t = detail::sinpx(z);
      z = -z;
      if(t < 0)
      {
         t = -t;
      }
      else
      {
         sresult = -sresult;
      }
      result = log(boost::math::constants::pi<T>()) - lgamma(z) - log(t);
   }
   else if((z != 1) && (z != 2))
   {
      T limit = (std::max)(z+1, T(10));
      T prefix = z * log(limit) - limit;
      T sum = detail::lower_gamma_series(z, limit, ::boost::math::tools::digits<T>()) / z;
      sum += detail::upper_gamma_fraction(z, limit, ::boost::math::tools::digits<T>());
      result = log(sum) + prefix;
   }
   if(sign)
      *sign = sresult;
   return result;
}
//
// This helper calculates tgamma(dz+1)-1 without cancellation errors,
// used by the upper incomplete gamma with z < 1:
//
template <class T, class Tag, class L>
T tgammap1m1_imp(T dz, Tag const& tag, const L& l)
{
   using namespace std;

   T result;
   if(dz < 0)
   {
      if(dz < -0.5)
      {
         // Best method is simply to subtract 1 from tgamma:
         result = boost::math::tgamma(1+dz) - 1;
      }
      else
      {
         // Use expm1 on lgamma:
         result = boost::math::expm1(-boost::math::log1p(dz) 
            + lgamma_small_imp(dz+2, dz + 1, dz, tag, l));
      }
   }
   else
   {
      if(dz < 2)
      {
         // Use expm1 on lgamma:
         result = boost::math::expm1(lgamma_small_imp(dz+1, dz, dz-1, tag, l));
      }
      else
      {
         // Best method is simply to subtract 1 from tgamma:
         result = boost::math::tgamma(1+dz) - 1;
      }
   }

   return result;
}

template <class T, class Tag>
T tgammap1m1_imp(T dz, Tag const& tag,
                 const ::boost::math::lanczos::undefined_lanczos& l)
{
   using namespace std; // ADL of std names
   //
   // There should be a better solution than this, but the
   // algebra isn't easy for the general case....
   // Start by subracting 1 from tgamma:
   //
   T result = gamma_imp(1 + dz, l) - 1;
   //
   // Test the level of cancellation error observed: we loose one bit
   // for each power of 2 the result is less than 1.  If we would get
   // more bits from our most precise lgamma rational approximation, 
   // then use that instead:
   //
   if((dz > -0.5) && (dz < 2) && (ldexp(1.0, boost::math::tools::digits<T>()) * fabs(result) < 1e34))
   {
      result = tgammap1m1_imp(dz, mpl::int_<113>(), boost::math::lanczos::lanczos24m113());
   }
   return result;
}

template <class T, class L>
inline T tgammap1m1_imp(T dz, const L& l)
{
   // Simple forwarding function.
   typedef typename mpl::if_c<
      (std::numeric_limits<T>::is_specialized && (std::numeric_limits<T>::digits <= 64)),
      mpl::int_<64>,
      typename mpl::if_c<
         (std::numeric_limits<T>::is_specialized && (std::numeric_limits<T>::digits <= 113)),
         mpl::int_<113>, mpl::int_<0> >::type
       >::type tag_type;
   return tgammap1m1_imp(dz, tag_type(), l);
}

//
// Series representation for upper fraction when z is small:
//
template <class T>
struct small_gamma2_series
{
   typedef T result_type;

   small_gamma2_series(T a_, T x_) : result(-x_), x(-x_), apn(a_+1), n(1){}

   T operator()()
   {
      T r = result / (apn);
      result *= x;
      result /= ++n;
      apn += 1;
      return r;
   }

private:
   T result, x, apn;
   int n;
};
//
// calculate power term prefix (z^a)(e^-z) used in the non-normalised
// incomplete gammas:
//
template <class T>
T full_igamma_prefix(T a, T z)
{
   using namespace std;

   T prefix;
   T alz = a * log(z);

   if(z >= 1)
   {
      if((alz < tools::log_max_value<T>()) && (-z > tools::log_min_value<T>()))
      {
         prefix = pow(z, a) * exp(-z);
      }
      else if(a >= 1)
      {
         prefix = pow(z / exp(z/a), a);
      }
      else
      {
         prefix = exp(alz - z);
      }
   }
   else
   {
      if(alz > tools::log_min_value<T>())
      {
         prefix = pow(z, a) * exp(-z);
      }
      else if(z/a < tools::log_max_value<T>())
      {
         prefix = pow(z / exp(z/a), a);
      }
      else
      {
         prefix = exp(alz - z);
      }
   }
   //
   // This error handling isn't very good: it happens after the fact
   // rather than before it...
   //
   if(boost::math::fpclassify(prefix) == FP_INFINITE)
      tools::overflow_error<T>(BOOST_CURRENT_FUNCTION, "Result of incomplete gamma function is too large to represent.");

   return prefix;
}
//
// Helper to compute log(1+x)-x:
//
template <class T>
T log1pmx(T x)
{
   boost::math::detail::log1p_series<T> s(x);
   s();
   boost::uintmax_t max_iter = BOOST_MATH_MAX_ITER;
#if BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x582))
   T zero = 0;
   T result = boost::math::tools::sum_series(s, tools::digits<T>() + 2, max_iter, zero);
#else
   T result = boost::math::tools::sum_series(s, tools::digits<T>() + 2, max_iter);
#endif
   tools::check_series_iterations(BOOST_CURRENT_FUNCTION, max_iter);
   return result;
}
//
// Compute (z^a)(e^-z)/tgamma(a)
// most if the error occurs in this function:
//
template <class T, class L>
T regularised_gamma_prefix(T a, T z, const L& l)
{
   using namespace std;
   T agh = a + L::g() - T(0.5);
   T prefix;
   T d = ((z - a) - L::g() + T(0.5)) / agh;

   if(a < 1)
   {
      //
      // We have to treat a < 1 as a special case because our Lanczos
      // approximations are optimised against the factorials with a > 1,
      // and for high precision types especially (128-bit reals for example)
      // very small values of a can give rather eroneous results for gamma
      // unless we do this:
      //
      // TODO: is this still required?  Lanczos approx should be better now?
      //
      if(z <= tools::log_min_value<T>())
      {
         // Oh dear, have to use logs, should be free of cancellation errors though:
         return exp(a * log(z) - z - lgamma_imp(a, l));
      }
      else
      {
         // direct calculation, no danger of overflow as gamma(a) < 1/a
         // for small a.
         return pow(z, a) * exp(-z) / gamma_imp(a, l);
      }
   }
   else if((fabs(d*d*a) <= 100) && (a > 150))
   {
      // special case for large a and a ~ z.
      prefix = a * log1pmx(d) + z * (0.5 - L::g()) / agh;
      prefix = exp(prefix);
   }
   else
   {
      //
      // general case.
      // direct computation is most accurate, but use various fallbacks
      // for different parts of the problem domain:
      //
      T alz = a * log(z / agh);
      T amz = a - z;
      if(((std::min)(alz, amz) <= tools::log_min_value<T>()) || ((std::max)(alz, amz) >= tools::log_max_value<T>()))
      {
         T amza = amz / a;
         if(((std::min)(alz, amz)/2 > tools::log_min_value<T>()) && ((std::max)(alz, amz)/2 < tools::log_max_value<T>()))
         {
            // compute square root of the result and then square it:
            T sq = pow(z / agh, a / 2) * exp(amz / 2);
            prefix = sq * sq;
         }
         else if(((std::min)(alz, amz)/4 > tools::log_min_value<T>()) && ((std::max)(alz, amz)/4 < tools::log_max_value<T>()) && (z > a))
         {
            // compute the 4th root of the result then square it twice:
            T sq = pow(z / agh, a / 4) * exp(amz / 4);
            prefix = sq * sq;
            prefix *= prefix;
         }
         else if((amza > tools::log_min_value<T>()) && (amza < tools::log_max_value<T>()))
         {
            prefix = pow((z * exp(amza)) / agh, a);
         }
         else
         {
            prefix = exp(alz + amz);
         }
      }
      else
      {
         prefix = pow(z / agh, a) * exp(amz);
      }
   }
   prefix *= sqrt(agh / boost::math::constants::e<T>()) / L::lanczos_sum_expG_scaled(a);
   return prefix;
}
//
// And again, without Lanczos support:
//
template <class T>
T regularised_gamma_prefix(T a, T z, const lanczos::undefined_lanczos&)
{
   using namespace std;

   T limit = (std::max)(T(10), a);
   T sum = detail::lower_gamma_series(a, limit, ::boost::math::tools::digits<T>()) / a;
   sum += detail::upper_gamma_fraction(a, limit, ::boost::math::tools::digits<T>());

   if(a < 10)
   {
      // special case for small a:
      T prefix = pow(z / 10, a);
      prefix *= exp(10-z);
      if(0 == prefix)
      {
         prefix = pow((z * exp((10-z)/a)) / 10, a);
      }
      prefix /= sum;
      return prefix;
   }

   T zoa = z / a;
   T amz = a - z;
   T alzoa = a * log(zoa);
   T prefix;
   if(((std::min)(alzoa, amz) <= tools::log_min_value<T>()) || ((std::max)(alzoa, amz) >= tools::log_max_value<T>()))
   {
      T amza = amz / a;
      if((amza <= tools::log_min_value<T>()) || (amza >= tools::log_max_value<T>()))
      {
         prefix = exp(alzoa + amz);
      }
      else
      {
         prefix = pow(zoa * exp(amza), a);
      }
   }
   else
   {
      prefix = pow(zoa, a) * exp(amz);
   }
   prefix /= sum;
   return prefix;
}
//
// Upper gamma fraction for very small a:
//
template <class T, class L>
T tgamma_small_upper_part(T a, T x, const L& l)
{
   using namespace std;  // ADL of std functions.
   //
   // Compute the full upper fraction (Q) when a is very small:
   //
   T result = tgammap1m1_imp(a, l) - powm1(x, a);
   result /= a;
   detail::small_gamma2_series<T> s(a, x);
   boost::uintmax_t max_iter = BOOST_MATH_MAX_ITER;
#if BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x582))
   T zero = 0;
   result -= pow(x, a) * tools::sum_series(s, boost::math::tools::digits<T>(), max_iter, zero);
#else
   result -= pow(x, a) * tools::sum_series(s, boost::math::tools::digits<T>(), max_iter);
#endif
   tools::check_series_iterations(BOOST_CURRENT_FUNCTION, max_iter);
   return result;
}
//
// Upper gamma fraction for integer a:
//
template <class T>
T finite_gamma_Q(T a, T x)
{
   //
   // Calculates normalised Q when a is an integer:
   //
   using namespace std;
   T sum = exp(-x);
   if(sum != 0)
   {
      T term = sum;
      for(unsigned n = 1; n < a; ++n)
      {
         term /= n;
         term *= x;
         sum += term;
      }
   }
   return sum;
}
//
// Upper gamma fraction for half integer a:
//
template <class T>
T finite_half_gamma_Q(T a, T x)
{
   //
   // Calculates normalised Q when a is a half-integer:
   //
   using namespace std;
   T e = boost::math::erfc(sqrt(x));
   if((e != 0) && (a > 1))
   {
      T term = exp(-x) / sqrt(constants::pi<T>() * x);
      term *= x;
      static const T half = T(1) / 2;
      term /= half;
      T sum = term;
      for(unsigned n = 2; n < a; ++n)
      {
         term /= n - half;
         term *= x;
         sum += term;
      }
      e += sum;
   }
   return e;
}
//
// Main incomplete gamma entry point, handles all four incomplete gamma's:
//
template <class T, class L>
T gamma_incomplete_imp(T a, T x, bool normalised, bool invert, const L& l)
{
   if(a <= 0)
      tools::domain_error<T>(BOOST_CURRENT_FUNCTION, "Argument a to the incomplete gamma function must be greater than zero (got a=%1%).", a);
   if(x < 0)
      tools::domain_error<T>(BOOST_CURRENT_FUNCTION, "Argument x to the incomplete gamma function must be >= 0 (got x=%1%).", x);

   using namespace std;

   T result;

   bool is_int = floor(a) == a;
   bool is_half_int = (floor(2 * a) == 2 * a) && !is_int;
   bool is_small_a = (a < 30) && (a <= x + 1);
   
   if(is_int && is_small_a && (x > 0.6))
   {
      // calculate Q via finite sum:
      invert = !invert;
      result = finite_gamma_Q(a, x);
      if(normalised == false)
         result *= gamma_imp(a, l);
   }
   else if(is_half_int && is_small_a && (x > 0.2))
   {
      // calculate Q via finite sum for half integer a:
      invert = !invert;
      result = finite_half_gamma_Q(a, x);
      if(normalised == false)
         result *= gamma_imp(a, l);
   }
   else if(x < 0.5)
   {
      //
      // Changeover criterion chosen to give a changeover at Q ~ 0.33
      //
      if(-0.4 / log(x) < a)
      {
         // Compute P:
         result = normalised ? regularised_gamma_prefix(a, x, l) : full_igamma_prefix(a, x);
         if(result != 0)
            result *= detail::lower_gamma_series(a, x, boost::math::tools::digits<T>()) / a;
      }
      else
      {
         // Compute Q:
         invert = !invert;
         result = tgamma_small_upper_part(a, x, l);
         if(normalised)
            result /= gamma_imp(a, l);
      }
   }
   else if(x < 1.1)
   {
      //
      // Changover here occurs when P ~ 0.6 or Q ~ 0.4:
      //
      if(x * 1.1f < a)
      {
         // Compute P:
         result = normalised ? regularised_gamma_prefix(a, x, l) : full_igamma_prefix(a, x);
         if(result != 0)
            result *= detail::lower_gamma_series(a, x, boost::math::tools::digits<T>()) / a;
      }
      else
      {
         // Compute Q:
         invert = !invert;
         result = tgamma_small_upper_part(a, x, l);
         if(normalised)
            result /= gamma_imp(a, l);
      }
   }
   else
   {
      //
      // Begin by testing whether we're in the "bad" zone
      // where the result will be near 0.5 and the usual
      // series and continued fractions are slow to converge:
      //
      bool use_temme = false;
      if(normalised && std::numeric_limits<T>::is_specialized && (a > 20))
      {
         T sigma = fabs((x-a)/a);
         if((a > 200) && (tools::digits<T>() <= 113))
         {
            //
            // This limit is chosen so that we use Temme's expansion
            // only if the result would be larger than about 10^-6.
            // Below that the regular series and continued fractions
            // converge OK, and if we use Temme's method we get increasing
            // errors from the dominant erfc term as it's (inexact) argument
            // increases in magnitude.
            //
            if(20 / a > sigma * sigma)
               use_temme = true;
         }
         else if(tools::digits<T>() <= 64)
         {
            // Note in this zone we can't use Temme's expansion for 
            // types longer than an 80-bit real:
            // it would require too many terms in the polynomials.
            if(sigma < 0.4)
               use_temme = true;
         }
      }
      if(use_temme)
      {
         //
         // Use compile time dispatch to the appropriate
         // Temme asymptotic expansion.  This may be dead code
         // if T does not have numeric limits support, or has
         // too many digits for the most precise version of
         // these expansions, in that case we'll be calling
         // an empty function.
         //
         typedef typename mpl::if_c<
            (std::numeric_limits<T>::digits == 0)
            ||
            (std::numeric_limits<T>::digits > 113),
            mpl::int_<0>,
            typename mpl::if_c<
               (std::numeric_limits<T>::digits <= 53),
               mpl::int_<53>,
               typename mpl::if_c<
                  (std::numeric_limits<T>::digits <= 64),
                  mpl::int_<64>,
                  mpl::int_<113>
               >::type
            >::type
         >::type tag_type;

         result = igamma_temme_large(a, x, static_cast<tag_type const*>(0));
         if(x >= a)
            invert = !invert;
      }
      else
      {
         //
         // Regular case where the result will not be too close to 0.5.
         //
         // Changeover here occurs at P ~ Q ~ 0.5
         //
         result = normalised ? regularised_gamma_prefix(a, x, l) : full_igamma_prefix(a, x);
         if(x < a)
         {
            // Compute P:
            if(result != 0)
               result *= detail::lower_gamma_series(a, x, boost::math::tools::digits<T>()) / a;
         }
         else
         {
            // Compute Q:
            invert = !invert;
            if(result != 0)
               result *= upper_gamma_fraction(a, x, tools::digits<T>());
         }
      }
   }

   if(invert)
   {
      T gam = normalised ? 1 : gamma_imp(a, l);
      result = gam - result;
   }

   return result;
}
//
// Ratios of two gamma functions:
//
template <class T, class L>
T tgamma_delta_ratio_imp_lanczos(T z, T delta, const L&)
{
   using namespace std;
   T zgh = z + L::g() - constants::half<T>();
   T result;
   if(fabs(delta) < 10)
   {
      result = exp((constants::half<T>() - z) * boost::math::log1p(delta / zgh));
   }
   else
   {
      result = pow(zgh / (zgh + delta), z - constants::half<T>());
   }
   result *= pow(constants::e<T>() / (zgh + delta), delta);
   result *= L::lanczos_sum(z) / L::lanczos_sum(z + delta);
   return result;
}
//
// And again without Lanczos support this time:
//
template <class T>
T tgamma_delta_ratio_imp_lanczos(T z, T delta, const lanczos::undefined_lanczos&)
{
   using namespace std;
   //
   // The upper gamma fraction is *very* slow for z < 6, actually it's very
   // slow to converge everywhere but recursing until z > 6 gets rid of the
   // worst of it's behaviour.
   //
   T prefix = 1;
   T zd = z + delta;
   while((zd < 6) && (z < 6))
   {
      prefix /= z;
      prefix *= zd;
      z += 1;
      zd += 1;
   }
   if(delta < 10)
   {
      prefix *= exp(-z * boost::math::log1p(delta / z));
   }
   else
   {
      prefix *= pow(z / zd, z);
   }
   prefix *= pow(constants::e<T>() / zd, delta);
   T sum = detail::lower_gamma_series(z, z, ::boost::math::tools::digits<T>()) / z;
   sum += detail::upper_gamma_fraction(z, z, ::boost::math::tools::digits<T>());
   T sumd = detail::lower_gamma_series(zd, zd, ::boost::math::tools::digits<T>()) / zd;
   sumd += detail::upper_gamma_fraction(zd, zd, ::boost::math::tools::digits<T>());
   sum /= sumd;
   if(fabs(tools::max_value<T>() / prefix) < fabs(sum))
      return tools::overflow_error<T>(BOOST_CURRENT_FUNCTION, "Result of tgamma is too large to represent.");
   return sum * prefix;
}

template <class T, class L>
T tgamma_delta_ratio_imp(T z, T delta, const L& l)
{
   using namespace std;

   if(z <= 0)
      tools::domain_error<T>(BOOST_CURRENT_FUNCTION, "Gamma function ratios only implemented for positive arguments (got a=%1%).", z);
   if(z+delta <= 0)
      tools::domain_error<T>(BOOST_CURRENT_FUNCTION, "Gamma function ratios only implemented for positive arguments (got b=%1%).", z+delta);

   if(floor(delta) == delta)
   {
      if(floor(z) == z)
      {
         //
         // Both z and delta are integers, see if we can just use table lookup
         // of the factorials to get the result:
         //
         if((z <= max_factorial<T>::value) && (z + delta <= max_factorial<T>::value))
         {
            return unchecked_factorial<T>(tools::real_cast<unsigned>(z) - 1) / unchecked_factorial<T>(tools::real_cast<unsigned>(z + delta) - 1);
         }
      }
      if(fabs(delta) < 20)
      {
         //
         // delta is a small integer, we can use a finite product:
         //
         if(delta == 0)
            return 1;
         if(delta < 0)
         {
            z -= 1;
            T result = z;
            while(0 != (delta += 1))
            {
               z -= 1;
               result *= z;
            }
            return result;
         }
         else
         {
            T result = 1 / z;
            while(0 != (delta -= 1))
            {
               z += 1;
               result /= z;
            }
            return result;
         }
      }
   }
   return tgamma_delta_ratio_imp_lanczos(z, delta, l);
}

template <class T, class L>
T gamma_P_derivative_imp(T a, T x, const L& l)
{
   //
   // Usual error checks first:
   //
   if(a <= 0)
      tools::domain_error<T>(BOOST_CURRENT_FUNCTION, "Argument a to the incomplete gamma function must be greater than zero (got a=%1%).", a);
   if(x < 0)
      tools::domain_error<T>(BOOST_CURRENT_FUNCTION, "Argument x to the incomplete gamma function must be >= 0 (got x=%1%).", x);
   //
   // Now special cases:
   //
   if(x == 0)
   {
      return (a > 1) ? 0 :
         (a == 1) ? 1 : tools::overflow_error<T>(BOOST_CURRENT_FUNCTION, 0);
   }
   //
   // Normal case:
   //
   T f1 = detail::regularised_gamma_prefix(a, x, l);
   if((x < 1) && (tools::max_value<T>() * x < f1))
   {
      // overflow:
      return tools::overflow_error<T>(BOOST_CURRENT_FUNCTION, 0);
   }

   f1 /= x;

   return f1;
}

} // namespace detail

template <class T>
inline T tgamma(T z)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename lanczos::lanczos_traits<typename remove_cv<T>::type>::value_type value_type;
   typedef typename lanczos::lanczos_traits<typename remove_cv<T>::type>::evaluation_type evaluation_type;
   return tools::checked_narrowing_cast<typename remove_cv<T>::type>(detail::gamma_imp(static_cast<value_type>(z), evaluation_type()), BOOST_CURRENT_FUNCTION);
}

template <class T>
inline T lgamma(T z, int* sign)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename lanczos::lanczos_traits<typename remove_cv<T>::type>::value_type value_type;
   typedef typename lanczos::lanczos_traits<typename remove_cv<T>::type>::evaluation_type evaluation_type;
   return tools::checked_narrowing_cast<typename remove_cv<T>::type>(detail::lgamma_imp(static_cast<value_type>(z), evaluation_type(), sign), BOOST_CURRENT_FUNCTION);
}

template <class T>
inline T lgamma(T x)
{
   BOOST_FPU_EXCEPTION_GUARD
   return ::boost::math::lgamma(x, 0);
}

template <class T>
inline T tgamma1pm1(T z)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename lanczos::lanczos_traits<typename remove_cv<T>::type>::value_type value_type;
   typedef typename lanczos::lanczos_traits<typename remove_cv<T>::type>::evaluation_type evaluation_type;
   typedef typename mpl::if_c<
      (std::numeric_limits<T>::is_specialized && (std::numeric_limits<T>::digits <= 64)),
      mpl::int_<64>,
      typename mpl::if_c<
         (std::numeric_limits<T>::is_specialized && (std::numeric_limits<T>::digits <= 113)),
         mpl::int_<113>, mpl::int_<0> >::type
       >::type tag_type;

   return tools::checked_narrowing_cast<typename remove_cv<T>::type>(detail::tgammap1m1_imp(static_cast<value_type>(z), tag_type(), evaluation_type()), BOOST_CURRENT_FUNCTION);
}

//
// Full upper incomplete gamma:
//
template <class T>
inline T tgamma(T a, T z)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename lanczos::lanczos_traits<typename remove_cv<T>::type>::value_type value_type;
   typedef typename lanczos::lanczos_traits<typename remove_cv<T>::type>::evaluation_type evaluation_type;
   return tools::checked_narrowing_cast<typename remove_cv<T>::type>(
      detail::gamma_incomplete_imp(static_cast<value_type>(a),
      static_cast<value_type>(z), false, true,
      evaluation_type()), BOOST_CURRENT_FUNCTION);
}
//
// Full lower incomplete gamma:
//
template <class T>
inline T tgamma_lower(T a, T z)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename lanczos::lanczos_traits<typename remove_cv<T>::type>::value_type value_type;
   typedef typename lanczos::lanczos_traits<typename remove_cv<T>::type>::evaluation_type evaluation_type;
   return tools::checked_narrowing_cast<typename remove_cv<T>::type>(
      detail::gamma_incomplete_imp(static_cast<value_type>(a),
      static_cast<value_type>(z), false, false,
      evaluation_type()), BOOST_CURRENT_FUNCTION);
}
//
// Regularised upper incomplete gamma:
//
template <class T>
inline T gamma_Q(T a, T z)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename lanczos::lanczos_traits<typename remove_cv<T>::type>::value_type value_type;
   typedef typename lanczos::lanczos_traits<typename remove_cv<T>::type>::evaluation_type evaluation_type;
   return tools::checked_narrowing_cast<typename remove_cv<T>::type>(
      detail::gamma_incomplete_imp(static_cast<value_type>(a),
      static_cast<value_type>(z), true, true,
      evaluation_type()), BOOST_CURRENT_FUNCTION);
}
//
// Regularised lower incomplete gamma:
//
template <class T>
inline T gamma_P(T a, T z)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename lanczos::lanczos_traits<typename remove_cv<T>::type>::value_type value_type;
   typedef typename lanczos::lanczos_traits<typename remove_cv<T>::type>::evaluation_type evaluation_type;
   return tools::checked_narrowing_cast<typename remove_cv<T>::type>(
      detail::gamma_incomplete_imp(static_cast<value_type>(a),
      static_cast<value_type>(z), true, false,
      evaluation_type()), BOOST_CURRENT_FUNCTION);
}

// ratios of gamma functions:
template <class T>
inline T tgamma_delta_ratio(T z, T delta)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename lanczos::lanczos_traits<typename remove_cv<T>::type>::value_type value_type;
   typedef typename lanczos::lanczos_traits<typename remove_cv<T>::type>::evaluation_type evaluation_type;
   return tools::checked_narrowing_cast<typename remove_cv<T>::type>(detail::tgamma_delta_ratio_imp(static_cast<value_type>(z), static_cast<value_type>(delta), evaluation_type()), BOOST_CURRENT_FUNCTION);
}
template <class T>
inline T tgamma_ratio(T a, T b)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename lanczos::lanczos_traits<typename remove_cv<T>::type>::value_type value_type;
   typedef typename lanczos::lanczos_traits<typename remove_cv<T>::type>::evaluation_type evaluation_type;
   return tools::checked_narrowing_cast<typename remove_cv<T>::type>(detail::tgamma_delta_ratio_imp(static_cast<value_type>(a), static_cast<value_type>(b - a), evaluation_type()), BOOST_CURRENT_FUNCTION);
}

template <class T>
T gamma_P_derivative(T a, T x)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename lanczos::lanczos_traits<typename remove_cv<T>::type>::value_type value_type;
   typedef typename lanczos::lanczos_traits<typename remove_cv<T>::type>::evaluation_type evaluation_type;
   return tools::checked_narrowing_cast<typename remove_cv<T>::type>(detail::gamma_P_derivative_imp(static_cast<value_type>(a), static_cast<value_type>(x), evaluation_type()), BOOST_CURRENT_FUNCTION);
}

} // namespace math
} // namespace boost

#ifdef BOOST_MSVC
# pragma warning(pop)
#endif

#include <boost/math/special_functions/detail/igamma_inverse.hpp>
#include <boost/math/special_functions/detail/gamma_inva.hpp>
#include <boost/math/special_functions/erf.hpp>

#endif // BOOST_MATH_SF_GAMMA_HPP


