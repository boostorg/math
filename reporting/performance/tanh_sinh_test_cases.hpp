//  (C) Copyright John Maddock 2021.
//  (C) Copyright Robert-van-Engelen 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/tools/config.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/pow.hpp>
#include <boost/math/special_functions/sign.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/cbrt.hpp>

template <class T>
void log_test_call(const T&);

template <class T>
inline T cot(const T& x)
{
   return 1 / tan(x);
}
template <class T>
inline T sec(const T& x)
{
   return 1 / cos(x);
}
template <class T>
inline T csc(const T& x)
{
   return 1 / sin(x);
}

template <class T>
inline T csch(const T& x)
{
   return 1 / sinh(x);
}

template <class T>
inline T sech(const T& x)
{
   return 1 / cosh(x);
}

template <class T>
inline T coth(const T& x)
{
   return 1 / tanh(x);
}

struct test_entry
{
   double (*proc)(const double&);
   double a;
   double b;
   double exact_result;
};

std::pair<const test_entry*, const test_entry*> get_tests();

template <class T>
inline T test_case_001(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / sqrt(x);
}


template <class T>
inline T test_case_002(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(4 - x * x);
}


template <class T>
inline T test_case_003(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(x);
}


template <class T>
inline T test_case_004(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * log(x);
}


template <class T>
inline T test_case_005(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(x) / sqrt(x);
}


template <class T>
inline T test_case_006(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 4 / (1 + x * x);
}


template <class T>
inline T test_case_007(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<4>(sin(x)) * boost::math::pow<2>(cos(x));
}


template <class T>
inline T test_case_008(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cos(x);
}


template <class T>
inline T test_case_009(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cos(log(x));
}


template <class T>
inline T test_case_010(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(4 * x - x * x);
}


template <class T>
inline T test_case_011(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 5 * x * x;
}


template <class T>
inline T test_case_012(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return pow(x, 0.125);
}


template <class T>
inline T test_case_013(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / x;
}


template <class T>
inline T test_case_014(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(x) / (1 - x);
}


template <class T>
inline T test_case_015(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(-1 / cos(x));
}


template <class T>
inline T test_case_016(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(x * (x + 88) * (x - 88) * (x + 47) * (x - 47) * (x + 117) * (x - 117));
}


template <class T>
inline T test_case_017(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (2 * log(1 / x) + 100);
}


template <class T>
inline T test_case_018(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 2 * x * x / (x + 1) / (x - 1) - x / log(x);
}


template <class T>
inline T test_case_019(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * log(1 + x);
}


template <class T>
inline T test_case_020(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * x * atan(x);
}


template <class T>
inline T test_case_021(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(x) * cos(x);
}


template <class T>
inline T test_case_022(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return atan(sqrt(x * x + 2)) / (1 + x * x) / sqrt(x * x + 2);
}


template <class T>
inline T test_case_023(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(x) * sqrt(x);
}


template <class T>
inline T test_case_024(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(1 - x * x);
}


template <class T>
inline T test_case_025(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(x) / sqrt(1 - x * x);
}


template <class T>
inline T test_case_026(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(log(x));
}


template <class T>
inline T test_case_027(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(cos(x));
}


template <class T>
inline T test_case_028(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(tan(x));
}


template <class T>
inline T test_case_029(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(x * x);
}


template <class T>
inline T test_case_030(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * sin(x) / (1 + boost::math::pow<2>(cos(x)));
}


template <class T>
inline T test_case_031(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (x - 2) / (pow((1 - x), 0.25)) / pow((1 + x), 0.75);
}


template <class T>
inline T test_case_032(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cos(boost::math::constants::pi<T>() * x) / sqrt(1 - x);
}


template <class T>
inline T test_case_033(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * x * log(x) / (x * x - 1) / (boost::math::pow<4>(x) + 1);
}


template <class T>
inline T test_case_034(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 - 2 * x + 2 * x * x);
}


template <class T>
inline T test_case_035(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(1 - 1 / x) / sqrt(boost::math::pow<3>(x) - boost::math::pow<4>(x));
}


template <class T>
inline T test_case_036(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(1 - 1 / x) * cos(1 / x - 1) / boost::math::pow<2>(x);
}


template <class T>
inline T test_case_037(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x / sqrt(1 - x * x);
}


template <class T>
inline T test_case_038(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<4>(1 - x) * boost::math::pow<4>(x) / (1 + x * x);
}


template <class T>
inline T test_case_039(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<4>(x) * boost::math::pow<4>((1 - x));
}


template <class T>
inline T test_case_040(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return atan(sqrt(boost::math::pow<2>(x) + 1)) / pow(boost::math::pow<2>(x) + 1, T(3) / 2);
}


template <class T>
inline T test_case_041(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 + boost::math::pow<2>(x) + boost::math::pow<4>(x) + boost::math::pow<6>(x));
}


template <class T>
inline T test_case_042(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (boost::math::constants::pi<T>() / 4 - x * tan(x)) * tan(x);
}


template <class T>
inline T test_case_043(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(x) / boost::math::pow<2>(sin(x));
}


template <class T>
inline T test_case_044(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(log(cos(x)));
}


template <class T>
inline T test_case_045(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(log(x)) / (boost::math::pow<2>(x) + x + 1);
}


template <class T>
inline T test_case_046(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(1 + boost::math::pow<2>(x)) / boost::math::pow<2>(x);
}


template <class T>
inline T test_case_047(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return atan(x) / (x * sqrt(1 - boost::math::pow<2>(x)));
}


template <class T>
inline T test_case_048(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(x) / (1 + boost::math::pow<4>(x)) / sqrt(1 - boost::math::pow<4>(x));
}


template <class T>
inline T test_case_049(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(-(boost::math::pow<2>(x))) * boost::math::pow<2>(log(x));
}


template <class T>
inline T test_case_050(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / sqrt(1 - boost::math::pow<2>(x));
}


template <class T>
inline T test_case_051(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(1 + x) * sin(2 * boost::math::constants::pi<T>() / (1 + x));
}


template <class T>
inline T test_case_052(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * pow(1 - x, 0.1);
}


template <class T>
inline T test_case_053(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(boost::math::pow<3>(sin(x))) * cos(x);
}


template <class T>
inline T test_case_054(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x / (exp(x) - 1);
}


template <class T>
inline T test_case_055(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 + exp(x));
}


template <class T>
inline T test_case_056(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(x);
}


template <class T>
inline T test_case_057(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 + x);
}


template <class T>
inline T test_case_058(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 0.92 * cosh(x) - cos(x);
}


template <class T>
inline T test_case_059(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 + boost::math::pow<4>(x));
}


template <class T>
inline T test_case_060(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (boost::math::pow<4>(x) + boost::math::pow<2>(x) + 0.9);
}


template <class T>
inline T test_case_061(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (boost::math::pow<2>(x) + 1.005);
}


template <class T>
inline T test_case_062(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 25 * exp(-25 * x);
}


template <class T>
inline T test_case_063(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(50) * exp(-50 * boost::math::constants::pi<T>() * boost::math::pow<2>(x));
}


template <class T>
inline T test_case_064(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 50 / (1 + 2500 * boost::math::pow<2>(x)) / boost::math::constants::pi<T>();
}


template <class T>
inline T test_case_065(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * sqrt(x);
}


template <class T>
inline T test_case_066(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(x);
}


template <class T>
inline T test_case_067(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cos(cos(x) + 3 * sin(x) + 2 * cos(2 * x) + 3 * sin(2 * x) + 3 * cos(3 * x));
}


template <class T>
inline T test_case_068(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(-x / 5) * (2 + sin(2 * x));
}


template <class T>
inline T test_case_069(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return pow(x, -T(1) / 3) * boost::math::pow<5>(1 - x);
}


template <class T>
inline T test_case_070(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(boost::math::constants::pi<T>() * x);
}


template <class T>
inline T test_case_071(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return pow(x, 0.25);
}


template <class T>
inline T test_case_072(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * sin(20 * x) * cos(50 * x);
}


template <class T>
inline T test_case_073(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(1 / x) * boost::math::pow<4>(x);
}


template <class T>
inline T test_case_074(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cos(8 * sin(x) - x);
}


template <class T>
inline T test_case_075(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * cos(x) * sin(30 * x);
}


template <class T>
inline T test_case_076(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (sqrt(1 - 0.81 * boost::math::pow<2>(sin(x))));
}


template <class T>
inline T test_case_077(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(x) * log(1 - x);
}


template <class T>
inline T test_case_078(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(sin(2 * x));
}


template <class T>
inline T test_case_079(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return pow(x, T(5) / 6) * pow(4 - x, T(1) / 6) / (5 - x) / (6 - x) / (7 - x);
}


template <class T>
inline T test_case_080(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 + boost::math::pow<3>(tan(x)));
}


template <class T>
inline T test_case_081(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return atan(sqrt(boost::math::pow<2>(x) + 1)) / sqrt(boost::math::pow<2>(x) + 1) / (boost::math::pow<2>(x) + 1);
}


template <class T>
inline T test_case_082(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(1 - 4 * cos(x) + 4);
}


template <class T>
inline T test_case_083(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 + pow(cos(x), x));
}


template <class T>
inline T test_case_084(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / sqrt(abs(x - 1));
}


template <class T>
inline T test_case_085(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / sqrt(x) * log(1 / x);
}


template <class T>
inline T test_case_086(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / pow(x, 0.25) * log(1 / x);
}


template <class T>
inline T test_case_087(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(x) / x;
}


template <class T>
inline T test_case_088(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(sin(x)) * boost::math::pow<4>(cos(x));
}


template <class T>
inline T test_case_089(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(1 + tan(x));
}


template <class T>
inline T test_case_090(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(sin(x));
}


template <class T>
inline T test_case_091(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * cos(boost::math::pow<2>(x) + 1);
}


template <class T>
inline T test_case_092(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(x) * (boost::math::pow<2>(x) - 2) * sin(x);
}


template <class T>
inline T test_case_093(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<3>(x) / (exp(x) - 1);
}


template <class T>
inline T test_case_094(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(x) / (1 + x);
}


template <class T>
inline T test_case_095(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(1 - 0.5 * boost::math::pow<2>(sin(x)));
}


template <class T>
inline T test_case_096(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<3>(log(x));
}


template <class T>
inline T test_case_097(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (boost::math::constants::pi<T>() * sqrt(1 - boost::math::pow<2>(x)));
}


template <class T>
inline T test_case_098(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log((1 + x) / (1 - x)) / 4 / log(2);
}


template <class T>
inline T test_case_099(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return ((T(23) / 25) * cosh(x) - cos(x));
}


template <class T>
inline T test_case_100(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 2 / (2 + sin(10 * boost::math::constants::pi<T>() * x));
}


template <class T>
inline T test_case_101(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(100 * boost::math::constants::pi<T>() * x) / boost::math::constants::pi<T>() / x;
}


template <class T>
inline T test_case_102(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 + cos(x));
}


template <class T>
inline T test_case_103(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1.005 + boost::math::pow<2>(x));
}


template <class T>
inline T test_case_104(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 4 * boost::math::pow<2>(boost::math::constants::pi<T>()) * x * sin(20 * boost::math::constants::pi<T>() * x) * cos(2 * boost::math::constants::pi<T>() * x);
}


template <class T>
inline T test_case_105(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 + boost::math::pow<2>((230 * x - 30)));
}


template <class T>
inline T test_case_106(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 50 * boost::math::pow<2>((sin(50 * boost::math::constants::pi<T>() * x) / (50 * boost::math::constants::pi<T>() * x)));
}


template <class T>
inline T test_case_107(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cos(cos(x) + 3 * sin(x) + 2 * cos(2 * x) + 3 * cos(3 * x));
}


template <class T>
inline T test_case_108(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<63>(x);
}


template <class T>
inline T test_case_109(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (x + 0.5);
}


template <class T>
inline T test_case_110(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(12.25 - boost::math::pow<2>((5 * x - 3)));
}


template <class T>
inline T test_case_111(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 10 / (1 + boost::math::pow<2>((10 * x - 4)));
}


template <class T>
inline T test_case_112(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / sqrt(x * (1 - x));
}


template <class T>
inline T test_case_113(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return pow(x, T(-3) / 4) * pow(1 - x, -0.25) / (3 - 2 * x);
}


template <class T>
inline T test_case_114(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (5 + 4 * cos(x));
}


template <class T>
inline T test_case_115(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(x) * sqrt(x / (1 - x));
}


template <class T>
inline T test_case_116(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(abs((tan(x) + sqrt(7)) / (tan(x) - sqrt(7))));
}


template <class T>
inline T test_case_117(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 + 25 * boost::math::pow<2>(x));
}


template <class T>
inline T test_case_118(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<3>(cos(x));
}


template <class T>
inline T test_case_119(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 + sin(x));
}


template <class T>
inline T test_case_120(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 - boost::math::pow<4>(x) / 2);
}


template <class T>
inline T test_case_121(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 + 100 * boost::math::pow<2>(x));
}


template <class T>
inline T test_case_122(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(x) / x;
}


template <class T>
inline T test_case_123(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (exp(x) - 1);
}


template <class T>
inline T test_case_124(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 + boost::math::pow<2>(x));
}


template <class T>
inline T test_case_125(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 + boost::math::pow<2>(sin(x)));
}


template <class T>
inline T test_case_126(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(abs(boost::math::pow<2>(x) - 0.25));
}


template <class T>
inline T test_case_127(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(sin(boost::math::constants::pi<T>() * x));
}


template <class T>
inline T test_case_128(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return abs(x - 0.4);
}


template <class T>
inline T test_case_129(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(-2 * abs(x - 0.4));
}


template <class T>
inline T test_case_130(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(abs(x - 0.3));
}


template <class T>
inline T test_case_131(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / boost::math::pow<5>((x + 0.01));
}


template <class T>
inline T test_case_132(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / sqrt(x + 0.0001);
}


template <class T>
inline T test_case_133(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (x + 0.0001);
}


template <class T>
inline T test_case_134(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (boost::math::pow<2>((230 * x - 30)) + 1);
}


template <class T>
inline T test_case_135(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (x + 0.01);
}


template <class T>
inline T test_case_136(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * sin(30 * x) * cos(x);
}


template <class T>
inline T test_case_137(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(-x) * sin(10 * x);
}


template <class T>
inline T test_case_138(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * boost::math::pow<2>((1 - x)) / boost::math::pow<3>((1 + x));
}


template <class T>
inline T test_case_139(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<12>(x) * boost::math::pow<12>((1 - x)) / 16 / (1 + boost::math::pow<2>(x));
}


template <class T>
inline T test_case_140(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<12>(x) * boost::math::pow<12>((1 - x)) / 16;
}


template <class T>
inline T test_case_141(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return atan(sqrt(boost::math::pow<2>(x) + 2)) / sqrt(boost::math::pow<2>(x) + 2) / (boost::math::pow<2>(x) + 1);
}


template <class T>
inline T test_case_142(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 + boost::math::pow<2>(x) + boost::math::pow<4>(x));
}


template <class T>
inline T test_case_143(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<5>((x + 0.1));
}


template <class T>
inline T test_case_144(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (boost::math::pow<2>(x) + 0.0001);
}


template <class T>
inline T test_case_145(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(x * (1 - x));
}


template <class T>
inline T test_case_146(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (9 + boost::math::pow<6>(x));
}


template <class T>
inline T test_case_147(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<5>((x + 0.01));
}


template <class T>
inline T test_case_148(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 - 0.99 * boost::math::pow<4>(x));
}


template <class T>
inline T test_case_149(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * exp(-3 * boost::math::pow<2>(x));
}


template <class T>
inline T test_case_150(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(x) * log(exp(1) / x);
}


template <class T>
inline T test_case_151(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (4 * sqrt(2) - 8 * boost::math::pow<3>(x) - 4 * sqrt(2) * boost::math::pow<4>(x) - 8 * boost::math::pow<5>(x)) / (1 - boost::math::pow<8>(x));
}


template <class T>
inline T test_case_152(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (3 - 2 * x) * 1 / sqrt(x * (1 - x));
}


template <class T>
inline T test_case_153(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cos(2 * boost::math::constants::pi<T>() * x) / sqrt(1 - x);
}


template <class T>
inline T test_case_154(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 2 * (1 - boost::math::pow<2>(x)) / (boost::math::pow<2>(tan(0.5)) + boost::math::pow<2>(x));
}


template <class T>
inline T test_case_155(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / sqrt(1 - 0.5 * boost::math::pow<2>(sin(x)));
}


template <class T>
inline T test_case_156(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 2 * (1 - boost::math::pow<2>(x)) / (cos(4 * atanh(x)) + cosh(2));
}


template <class T>
inline T test_case_157(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return acos(boost::math::pow<2>(cos(boost::math::constants::pi<T>() * x))) / 3 / sin(boost::math::constants::pi<T>() * x) + 0.5;
}


template <class T>
inline T test_case_158(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(x);
}


template <class T>
inline T test_case_159(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(x) * cos(x);
}


template <class T>
inline T test_case_160(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(2 * x + 7);
}


template <class T>
inline T test_case_161(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<3>(x) * exp(-x);
}


template <class T>
inline T test_case_162(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(10 * x) * exp(-x);
}


template <class T>
inline T test_case_163(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 + boost::math::pow<64>(x));
}


template <class T>
inline T test_case_164(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(x) / (x - 1) - 1 / log(x);
}


template <class T>
inline T test_case_165(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(boost::math::pow<2>(x) * boost::math::pow<2>((1 - x)));
}


template <class T>
inline T test_case_166(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(x) * sin(10 * x);
}


template <class T>
inline T test_case_167(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cos(8 * sin(x) - x) / boost::math::constants::pi<T>();
}


template <class T>
inline T test_case_168(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * sin(10 * x) * cos(x);
}


template <class T>
inline T test_case_169(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 23 * cosh(x) / 25 - cos(x);
}


template <class T>
inline T test_case_170(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<3>(x) * sin(10 * boost::math::constants::pi<T>() * x) / sqrt(1 - boost::math::pow<2>(x));
}


template <class T>
inline T test_case_171(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cos(x) * log(x);
}


template <class T>
inline T test_case_172(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(1 + x) * sin(10 * x);
}


template <class T>
inline T test_case_173(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   // x^(PI()/4)*SIN(PI()/(8-4*x))
   return pow(x, (boost::math::constants::pi<T>() / 4)) * sin(boost::math::constants::pi<T>() / (8 - 4 * x));
}


template <class T>
inline T test_case_174(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return pow(x, x);
}


template <class T>
inline T test_case_175(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return pow(x, -x);
}


template <class T>
inline T test_case_176(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (4 + cos(x)) * log(3 + cos(x));
}


template <class T>
inline T test_case_177(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<3>(x) - 3 * boost::math::pow<2>(x);
}


template <class T>
inline T test_case_178(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 2 * sqrt(1 - boost::math::pow<2>(x));
}


template <class T>
inline T test_case_179(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<3>(x) * sin(5 * boost::math::constants::pi<T>() * x);
}


template <class T>
inline T test_case_180(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / x * sin(1 / x);
}


template <class T>
inline T test_case_181(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(x) / ((1 + boost::math::pow<4>(x)) * sqrt(1 - boost::math::pow<4>(x)));
}


template <class T>
inline T test_case_182(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(x) / sqrt(1 - x);
}


template <class T>
inline T test_case_183(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(x) / boost::math::pow<2>((1 + log(x)));
}


template <class T>
inline T test_case_184(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return atan(x) / x / (boost::math::pow<2>(x) + 1);
}


template <class T>
inline T test_case_185(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(cos(3 * x)) / (5 - 4 * cos(2 * x));
}


template <class T>
inline T test_case_186(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * tan(x);
}


template <class T>
inline T test_case_187(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(x) / (1 - cos(x));
}


template <class T>
inline T test_case_188(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<4>(x) * boost::math::pow<4>((1 - x)) / (1 + boost::math::pow<2>(x));
}


template <class T>
inline T test_case_189(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   // x^4*(1-x)^4/(1+x^2)
   return boost::math::pow<2>((boost::math::pow<2>(x) - x)) / boost::math::pow<2>((boost::math::pow<3>(x) - 3 * x + 1));
}


template <class T>
inline T test_case_190(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cos(sin(x) - x);
}


template <class T>
inline T test_case_191(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / boost::math::constants::pi<T>() / sqrt(1 - boost::math::pow<2>(x));
}


template <class T>
inline T test_case_192(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / 4 / log(2) * log((1 + x) / (1 - x));
}


template <class T>
inline T test_case_193(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>((1 + x)) * exp(-x) / boost::math::pow<2>((1 + boost::math::pow<2>(x)));
}


template <class T>
inline T test_case_194(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(1 - x) / x;
}


template <class T>
inline T test_case_195(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x / (1 - x);
}


template <class T>
inline T test_case_196(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 + pow(tan(x), sqrt(2)));
}


template <class T>
inline T test_case_197(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(1 - boost::math::pow<4>(x));
}


template <class T>
inline T test_case_198(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (boost::math::pow<2>((x - 0.3)) + 0.01) + 1 / (boost::math::pow<2>((x - 0.9)) + 0.04) - 6;
}


template <class T>
inline T test_case_199(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cos(boost::math::pow<2>(x)) / sqrt(1 - boost::math::pow<2>(x));
}


template <class T>
inline T test_case_200(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (exp(2 * x));
}


template <class T>
inline T test_case_201(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return pow(x, 0.1);
}


template <class T>
inline T test_case_202(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<4>(x) * log(x + sqrt(boost::math::pow<2>(x) + 1));
}


template <class T>
inline T test_case_203(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<3>(x) * exp(x);
}


template <class T>
inline T test_case_204(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 100 / boost::math::pow<2>(x) * sin(10 / x);
}


template <class T>
inline T test_case_205(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<3>(x) / sqrt(1 - boost::math::pow<2>(x)) * sin(3 * x);
}


template <class T>
inline T test_case_206(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<15>(x) * exp(x - 1);
}


template <class T>
inline T test_case_207(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (x - 2);
}


template <class T>
inline T test_case_208(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(exp(x)) / sqrt(x);
}


template <class T>
inline T test_case_209(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (boost::math::pow<2>(x) + 2 * x + 4) / (boost::math::pow<4>(x) - 7 * boost::math::pow<2>(x) + 2 * x + 17);
}


template <class T>
inline T test_case_210(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return pow(tan(x), (1 / boost::math::constants::pi<T>()));
}


template <class T>
inline T test_case_211(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 5 * sin(x) + (9 * x - 4) * (9 * x - 8) * (3 * x - 4) * (9 * x - 10) * (boost::math::constants::pi<T>() - 2 * x) / (1 + boost::math::pow<4>((90 * x - 110)));
}


template <class T>
inline T test_case_212(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return - log(cos(boost::math::constants::pi<T>() / 2 * tanh(boost::math::constants::pi<T>() / 2 * sinh(x)))) * boost::math::constants::pi<T>() / 2 * cosh(x) / boost::math::pow<2>(cosh(boost::math::constants::pi<T>() / 2 * sinh(x)));
}


template <class T>
inline T test_case_213(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (16 * boost::math::pow<2>((x - boost::math::constants::pi<T>() / 4)) + T(1) / 16);
}


template <class T>
inline T test_case_214(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cos(64 * sin(x));
}


template <class T>
inline T test_case_215(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(20 * (x - 1)) * sin(256 * x);
}


template <class T>
inline T test_case_216(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (boost::math::constants::pi<T>() * sqrt(1 - boost::math::pow<2>(x)));
}


template <class T>
inline T test_case_217(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / 4 / log(2) * log((1 + x) / (1 - x));
}


template <class T>
inline T test_case_218(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 2 / boost::math::constants::pi<T>() * sqrt(1 - boost::math::pow<2>(x));
}


template <class T>
inline T test_case_219(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 2 / boost::math::constants::pi<T>() / (1 + boost::math::pow<2>(x));
}


template <class T>
inline T test_case_220(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(x) / (1 + 4 * x + 3 * boost::math::pow<2>(x) - 4 * boost::math::pow<3>(x) - 2 * boost::math::pow<4>(x) + 2 * boost::math::pow<5>(x) + boost::math::pow<6>(x));
}


template <class T>
inline T test_case_221(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 2 / (1 + boost::math::pow<2>(x)) + 1 / (1 + boost::math::pow<2>((x - 10)));
}


template <class T>
inline T test_case_222(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(x) / log(1 + x);
}


template <class T>
inline T test_case_223(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return tan(boost::math::constants::pi<T>() * x / 2) * sin(4 * boost::math::constants::pi<T>() * x);
}


template <class T>
inline T test_case_224(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cos(sin(x));
}


template <class T>
inline T test_case_225(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * boost::math::constants::pi<T>() * sin(30 * x) / sqrt(4 * boost::math::pow<2>(boost::math::constants::pi<T>()) - boost::math::pow<2>(x));
}


template <class T>
inline T test_case_226(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(x) * log(x) * log(x + 1) / (x + 1);
}


template <class T>
inline T test_case_227(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(abs(sin(x) / x));
}


template <class T>
inline T test_case_228(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 3 * boost::math::pow<4>(x) * log(x);
}


template <class T>
inline T test_case_229(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>((boost::math::constants::pi<T>() - x)) * log(2 * boost::math::pow<2>(sin(x / 2)));
}


template <class T>
inline T test_case_230(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(x) / boost::math::pow<2>(sin(x));
}


template <class T>
inline T test_case_231(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (16 * x - 16) / (boost::math::pow<4>(x) - 2 * boost::math::pow<3>(x) + 4 * x - 4);
}


template <class T>
inline T test_case_232(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<4>(x) * (1 - boost::math::pow<4>(x)) / (1 + boost::math::pow<2>(x));
}


template <class T>
inline T test_case_233(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return pow(T(2), x) * x / (pow(T(2), x) - 1);
}


template <class T>
inline T test_case_234(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (exp(-x) - exp(-10 * x));
}


template <class T>
inline T test_case_235(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(x / (1 + boost::math::pow<2>(x)));
}


template <class T>
inline T test_case_236(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(sin(x));
}


template <class T>
inline T test_case_237(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 2 * (exp(-9 * boost::math::pow<2>(x)) + exp(-1024 * boost::math::pow<2>((x - 0.25)))) / sqrt(boost::math::constants::pi<T>());
}


template <class T>
inline T test_case_238(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 50 / (boost::math::constants::pi<T>() * (2500 * boost::math::pow<2>(x) + 1));
}


template <class T>
inline T test_case_239(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / sqrt(abs(x));
}


template <class T>
inline T test_case_240(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 2 * (exp(-9 * boost::math::pow<2>(x)) + exp(-1024 * boost::math::pow<2>((x - 0.25)))) / sqrt(boost::math::constants::pi<T>());
}


template <class T>
inline T test_case_241(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / sqrt(1 - boost::math::pow<2>(x) * boost::math::pow<2>(sin(x)));
}


template <class T>
inline T test_case_242(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(1 + x) / x;
}


template <class T>
inline T test_case_243(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<4>(x) * log(x + sqrt(boost::math::pow<2>(x) + 1));
}


template <class T>
inline T test_case_244(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cos(x) * (boost::math::pow<2>(x) + x + 1);
}


template <class T>
inline T test_case_245(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 6 * x - boost::math::pow<4>(x) - 1;
}


template <class T>
inline T test_case_246(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(x) + x - 2;
}


template <class T>
inline T test_case_247(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(x) - 4 * x - 4 + 4 * log(4);
}


template <class T>
inline T test_case_248(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(x) - 5 * x + 3;
}


template <class T>
inline T test_case_249(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(x) - 20 * x + 90;
}


template <class T>
inline T test_case_250(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<3>(x) - x - 1;
}


template <class T>
inline T test_case_251(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<3>(x) - 3 * boost::math::pow<2>(x) + 4;
}


template <class T>
inline T test_case_252(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<9>(x) - boost::math::pow<8>(x) + boost::math::pow<7>(x) - boost::math::pow<6>(x) + boost::math::pow<5>(x) - boost::math::pow<4>(x) - boost::math::pow<3>(x) + 2 * boost::math::pow<2>(x) - x + 0.5;
}


template <class T>
inline T test_case_253(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 2 * boost::math::pow<4>(x) - 9 * exp(x) - 22.5;
}


template <class T>
inline T test_case_254(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(x) - log(x) - 0.7;
}


template <class T>
inline T test_case_255(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 7 * sin(x) - exp(-x) * cos(x) - 0.7;
}


template <class T>
inline T test_case_256(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(x) - 20;
}


template <class T>
inline T test_case_257(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / x / (1 + boost::math::pow<2>(log(x)));
}


template <class T>
inline T test_case_258(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(1 + sqrt(x));
}


template <class T>
inline T test_case_259(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / sqrt(sin(x));
}


template <class T>
inline T test_case_260(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cos(100 * sin(x));
}


template <class T>
inline T test_case_261(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (0.1 * cos(100 * sin(x)) - 1) * boost::math::sign(x - 1);
}


template <class T>
inline T test_case_262(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (boost::math::pow<6>(x) + 0.9);
}


template <class T>
inline T test_case_263(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(abs(x + 0.5));
}


template <class T>
inline T test_case_264(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>((sin(50 * boost::math::constants::pi<T>() * x)));
}


template <class T>
inline T test_case_265(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x / (exp(x) + 1);
}


template <class T>
inline T test_case_266(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>((1 / cosh(10 * (x - 0.2)))) + boost::math::pow<4>((1 / cosh(100 * (x - 0.4)))) + boost::math::pow<6>((1 / cosh(1000 * (x - 0.6))));
}


template <class T>
inline T test_case_267(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(abs(x - 0.7));
}


template <class T>
inline T test_case_268(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(cos(x));
}


template <class T>
inline T test_case_269(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (sqrt(x) + boost::math::cbrt(x));
}


template <class T>
inline T test_case_270(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(-x) * sin(50 * x);
}


template <class T>
inline T test_case_271(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cos(x) + 5 * cos(1.6 * x) - 2 * cos(2 * x) + 5 * cos(4.5 * x) + 7 * cos(9 * x);
}


template <class T>
inline T test_case_272(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(x) * cos(x);
}


template <class T>
inline T test_case_273(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(-log(x));
}


template <class T>
inline T test_case_274(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(-20 * x);
}


template <class T>
inline T test_case_275(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 2 * boost::math::pow<2>(x) / ((x - 1) * (x + 1)) - x / log(x);
}


template <class T>
inline T test_case_276(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(x) * log(1 / x);
}


template <class T>
inline T test_case_277(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return T(23) / T(50) * (exp(x) + exp(-x)) - cos(x);
}


template <class T>
inline T test_case_278(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cos(32 * sin(x));
}


template <class T>
inline T test_case_279(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>((abs(x - T(1) / 3)));
}


template <class T>
inline T test_case_280(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>((abs(x - boost::math::constants::pi<T>() / 4)));
}


template <class T>
inline T test_case_281(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(20 * (x - 1)) * sin(pow(T(2), 5) * x);
}


template <class T>
inline T test_case_282(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (1 - boost::math::pow<2>((sin(boost::math::constants::pi<T>() / 12)))) / (1 - (sin(boost::math::constants::pi<T>() / 12)) * cos(x));
}


template <class T>
inline T test_case_283(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(x) * log(x);
}


template <class T>
inline T test_case_284(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<3>(x) * log(abs((boost::math::pow<2>(x) - 1) * (boost::math::pow<2>(x) - 2)));
}


template <class T>
inline T test_case_285(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / sqrt(abs(boost::math::pow<2>(x) + 2 * x - 2));
}


template <class T>
inline T test_case_286(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return pow(x, 0.2);
}


template <class T>
inline T test_case_287(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(x) * pow(1 - x, 0.3);
}


template <class T>
inline T test_case_288(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return pow(x, T(5) / 2);
}


template <class T>
inline T test_case_289(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (boost::math::pow<2>(x) + x + 1) * cos(x);
}


template <class T>
inline T test_case_290(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 - 0.998 * boost::math::pow<2>(x));
}


template <class T>
inline T test_case_291(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 - boost::math::cbrt(boost::math::pow<2>((x - boost::math::constants::pi<T>() / (2 * exp(1)))));
}


template <class T>
inline T test_case_292(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(2 * sin(x / 2));
}


template <class T>
inline T test_case_293(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(x) / (exp(x) - 0.99);
}


template <class T>
inline T test_case_294(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (3 + (x - 1) * exp(x)) / boost::math::pow<2>((3 - exp(x)));
}


template <class T>
inline T test_case_295(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(x) / boost::math::pow<2>((3 - exp(x)));
}


template <class T>
inline T test_case_296(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 + boost::math::pow<2>(tan(x));
}


template <class T>
inline T test_case_297(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (2 * (1 - x) * sin(x) + cos(x)) / sqrt(1 - x);
}


template <class T>
inline T test_case_298(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<4>((1 - pow(x, 0.25)));
}


template <class T>
inline T test_case_299(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 5 / (1 - exp(-5)) * exp(-5 * x);
}


template <class T>
inline T test_case_300(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 100 * exp(-100 * x);
}


template <class T>
inline T test_case_301(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 - 0.98 * boost::math::pow<2>(x));
}


template <class T>
inline T test_case_302(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * (x - 2) * log(x / 2);
}


template <class T>
inline T test_case_303(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<3>((x - 15)) * log(x / 2);
}


template <class T>
inline T test_case_304(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (sqrt(x) - 1 / sqrt(x)) / (1 + x);
}


template <class T>
inline T test_case_305(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(x) / sqrt(1 - x);
}


template <class T>
inline T test_case_306(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (boost::math::pow<2>(x) + x) / (1 + boost::math::pow<5>(x));
}


template <class T>
inline T test_case_307(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>((1 - sqrt(x)));
}


template <class T>
inline T test_case_308(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(x) * pow(1 - boost::math::pow<5>(x), T(-3) / 5);
}


template <class T>
inline T test_case_309(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 - x) - 3 * boost::math::pow<2>(x) / (1 - boost::math::pow<3>(x));
}


template <class T>
inline T test_case_310(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (sqrt(x) - 1 / sqrt(x)) * x / (1 - boost::math::pow<2>(x));
}


template <class T>
inline T test_case_311(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (boost::math::pow<3>(x) + pow(x, T(7) / 2) - 2 * boost::math::pow<7>(x)) / (1 - x);
}


template <class T>
inline T test_case_312(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(x) / sqrt(1 - 4 * cos(x) + 4);
}


template <class T>
inline T test_case_313(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return tan(x) / (boost::math::pow<3>(cos(x)) + 1 / boost::math::pow<3>(cos(x)));
}


template <class T>
inline T test_case_314(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(4 * sin(x)) * sin(2 * x);
}


template <class T>
inline T test_case_315(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(4 * tan(x)) * sin(2 * x);
}


template <class T>
inline T test_case_316(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 - x) * log(x);
}


template <class T>
inline T test_case_317(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log((1 - x) / x) / (1 + boost::math::pow<2>(x));
}


template <class T>
inline T test_case_318(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<3>(x) / (1 - x) - 3 * boost::math::pow<11>(x) / (1 - boost::math::pow<3>(x));
}


template <class T>
inline T test_case_319(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * exp(-16 * boost::math::pow<2>(x));
}


template <class T>
inline T test_case_320(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * exp(x) / boost::math::pow<2>((1 + x));
}


template <class T>
inline T test_case_321(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x / (exp(x) + 2 * exp(-x) - 2);
}


template <class T>
inline T test_case_322(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (3 * exp(1 - boost::math::pow<-3>(x)) / (1 - boost::math::pow<3>(x)) - exp(1 - 1 / x) / (1 - x)) / x;
}


template <class T>
inline T test_case_323(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x / sqrt(cosh(6) - cosh(6 * x)) / sinh(3 * x);
}


template <class T>
inline T test_case_324(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(3 * x) * cos(2 * x) / sin(x);
}


template <class T>
inline T test_case_325(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(3 * x) / sin(x);
}


template <class T>
inline T test_case_326(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(5 * x) / sin(x);
}


template <class T>
inline T test_case_327(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cos(3 * x) / (1 + cos(x) / 2);
}


template <class T>
inline T test_case_328(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(3 * x) * sin(x) / (1 - 4 * cos(x) + 4);
}


template <class T>
inline T test_case_329(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (cos(3 * x) - 0.5 * cos(4 * x)) / (1 - cos(x) + 0.25);
}


template <class T>
inline T test_case_330(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 - boost::math::pow<2>(x) / 2);
}


template <class T>
inline T test_case_331(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(cos(x)) * sin(4 * x) / tan(x);
}


template <class T>
inline T test_case_332(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<3>(cos(x)) * sin(4 * x) / sin(x);
}


template <class T>
inline T test_case_333(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 + 10 * boost::math::pow<2>(x));
}


template <class T>
inline T test_case_334(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<3>(sin(x)) / (4 + 3 * cos(x));
}


template <class T>
inline T test_case_335(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(tan(x)) / (sin(x) + cos(x)) / sin(x);
}


template <class T>
inline T test_case_336(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * sin(x) / (1 - cos(x));
}


template <class T>
inline T test_case_337(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (x - sin(x)) / (1 - cos(x));
}


template <class T>
inline T test_case_338(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 - 6 * cos(x) + 9);
}


template <class T>
inline T test_case_339(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x / boost::math::pow<2>((sin(x) + 3 * cos(x)));
}


template <class T>
inline T test_case_340(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x / (9 * boost::math::pow<2>(cos(x)) + 16 * boost::math::pow<2>(sin(x)));
}


template <class T>
inline T test_case_341(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (1 - x / tan(x)) / boost::math::pow<2>(sin(x));
}


template <class T>
inline T test_case_342(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 + 5 * boost::math::pow<2>(x));
}


template <class T>
inline T test_case_343(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<3>(x) * cos(x) / boost::math::pow<3>(sin(x));
}


template <class T>
inline T test_case_344(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * sin(10 * boost::math::constants::pi<T>() * x);
}


template <class T>
inline T test_case_345(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(3 * cos(x)) * sin(x);
}


template <class T>
inline T test_case_346(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(4 * cos(x)) * sin(4 * sin(x)) / sin(x);
}


template <class T>
inline T test_case_347(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * exp(-4 * boost::math::pow<2>(tan(x))) * (4 - boost::math::pow<2>(cos(x))) / boost::math::pow<4>(cos(x)) * tan(x);
}


template <class T>
inline T test_case_348(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 + boost::math::pow<2>(x) + boost::math::pow<4>(x));
}


template <class T>
inline T test_case_349(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / sqrt(log(1 / x));
}


template <class T>
inline T test_case_350(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(sin(x));
}


template <class T>
inline T test_case_351(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / sqrt(x) * log(log(1 / x));
}


template <class T>
inline T test_case_352(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(4 + 2 * cos(x));
}


template <class T>
inline T test_case_353(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(1 - 8 * cos(x) + 16);
}


template <class T>
inline T test_case_354(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log((1 + sin(1) * boost::math::pow<2>(cos(x))) / (1 - sin(1) * boost::math::pow<2>(cos(x))));
}


template <class T>
inline T test_case_355(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(2 * tan(x));
}


template <class T>
inline T test_case_356(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(9 + 16 * boost::math::pow<2>(tan(x)));
}


template <class T>
inline T test_case_357(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(x) + 3 * boost::math::pow<3>(x);
}


template <class T>
inline T test_case_358(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (1 - x) * log(x) / (1 + x);
}


template <class T>
inline T test_case_359(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * log(x) / boost::math::pow<2>((1 + boost::math::pow<2>(x)));
}


template <class T>
inline T test_case_360(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 - x) + x * log(x) / boost::math::pow<2>((1 - x));
}


template <class T>
inline T test_case_361(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * log(x) / sqrt(1 - boost::math::pow<4>(x));
}


template <class T>
inline T test_case_362(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * log(x) / boost::math::cbrt(boost::math::pow<2>((1 - boost::math::pow<3>(x))));
}


template <class T>
inline T test_case_363(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(log(x)) / (boost::math::pow<2>(x) - x + 1);
}


template <class T>
inline T test_case_364(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (1 - x) / (1 + x) / log(x);
}


template <class T>
inline T test_case_365(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(x) / sqrt(log(1 / x));
}


template <class T>
inline T test_case_366(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / log(x) + 1 / (1 - x);
}


template <class T>
inline T test_case_367(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (boost::math::pow<2>(boost::math::constants::pi<T>()) + boost::math::pow<2>(log(x))) / (1 + boost::math::pow<2>(x));
}


template <class T>
inline T test_case_368(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (x * log(x) + 1 - x) / x / boost::math::pow<2>(log(x)) * log(1 + x);
}


template <class T>
inline T test_case_369(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * exp(x) * log(1 - x);
}


template <class T>
inline T test_case_370(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<3>(sin(x)) * log(sin(x)) / sqrt(1 + boost::math::pow<2>(sin(x)));
}


template <class T>
inline T test_case_371(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x / (1 + sin(x));
}


template <class T>
inline T test_case_372(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(1 + tan(x));
}


template <class T>
inline T test_case_373(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(1 / tan(x) - 1);
}


template <class T>
inline T test_case_374(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(log(1 / x)) / sqrt(log(1 / x));
}


template <class T>
inline T test_case_375(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(x) / (1 + boost::math::pow<2>(x));
}


template <class T>
inline T test_case_376(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(x) / sqrt(1 - boost::math::pow<2>(x));
}


template <class T>
inline T test_case_377(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(1 - boost::math::pow<2>(x)) * log(x);
}


template <class T>
inline T test_case_378(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(log(x)) * (1 + boost::math::pow<2>(x)) / (1 + boost::math::pow<4>(x));
}


template <class T>
inline T test_case_379(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<3>(log(x)) / (1 + x);
}


template <class T>
inline T test_case_380(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<4>(log(x)) / (1 + boost::math::pow<2>(x));
}


template <class T>
inline T test_case_381(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (1 - x) / (1 + x) * boost::math::pow<2>(x) / (1 + boost::math::pow<2>(x)) / log(x);
}


template <class T>
inline T test_case_382(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (boost::math::pow<3>(x) - boost::math::pow<2>(x)) / log(x);
}


template <class T>
inline T test_case_383(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (boost::math::pow<4>(x) - boost::math::pow<3>(x)) * boost::math::pow<4>(x) / log(x);
}


template <class T>
inline T test_case_384(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (boost::math::pow<5>(x) - 1) / x / boost::math::pow<2>(log(x)) - 5 / log(x);
}


template <class T>
inline T test_case_385(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log((1 + x) / 2) / (1 - x);
}


template <class T>
inline T test_case_386(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(1 + x) / boost::math::pow<2>((3 * x + 3));
}


template <class T>
inline T test_case_387(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(1 + x) * (1 + boost::math::pow<2>(x)) / boost::math::pow<4>((1 + x));
}


template <class T>
inline T test_case_388(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log((1 + x) / (1 - x)) / (1 + boost::math::pow<2>(x));
}


template <class T>
inline T test_case_389(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(log(1 / x)) / (1 + x);
}


template <class T>
inline T test_case_390(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (1 - x) * exp(-x) * log(x);
}


template <class T>
inline T test_case_391(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(sin(x)) * boost::math::pow<2>(sin(x));
}


template <class T>
inline T test_case_392(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(x) * log(cot(x / 2));
}


template <class T>
inline T test_case_393(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(1 + cos(x) / 2) / cos(x);
}


template <class T>
inline T test_case_394(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(sin(x)) * tan(x);
}


template <class T>
inline T test_case_395(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(x) * sin(2 * x);
}


template <class T>
inline T test_case_396(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<4>(x) * boost::math::pow<2>(sin(boost::math::constants::pi<T>() * x));
}


template <class T>
inline T test_case_397(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<3>(sec(x));
}


template <class T>
inline T test_case_398(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return - log(1 - x) / (1 - x) * boost::math::pow<2>(log(x));
}


template <class T>
inline T test_case_399(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / sqrt(x) * 1 / sqrt(1 - x);
}


template <class T>
inline T test_case_400(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (4 * boost::math::pow<3>(x) - 1) / (2 * boost::math::pow<2>(x) + 1);
}


template <class T>
inline T test_case_401(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 2 / (2 * boost::math::pow<2>(x) - 1);
}


template <class T>
inline T test_case_402(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return atan((x + 1) / 3) - atan((x - 1) / 3);
}


template <class T>
inline T test_case_403(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cos(20 * boost::math::pow<2>(x));
}


template <class T>
inline T test_case_404(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 - log(x)) / sqrt(-log(x));
}


template <class T>
inline T test_case_405(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x / sin(x) / cos(x);
}


template <class T>
inline T test_case_406(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x / sin(x);
}


template <class T>
inline T test_case_407(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(cot(x));
}


template <class T>
inline T test_case_408(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(log(1 / x));
}


template <class T>
inline T test_case_409(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(x) * cos(10 * x);
}


template <class T>
inline T test_case_410(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return atan(x) / x;
}


template <class T>
inline T test_case_411(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(3 * x) * cos(5 * x) / x;
}


template <class T>
inline T test_case_412(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (4 + 3 * cos(x));
}


template <class T>
inline T test_case_413(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / boost::math::pow<2>((4 + 3 * sin(x)));
}


template <class T>
inline T test_case_414(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<3>(sin(5 * x));
}


template <class T>
inline T test_case_415(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x / (1 + sin(3 * x));
}


template <class T>
inline T test_case_416(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cos(3 * x) / (1 + cos(3 * x));
}


template <class T>
inline T test_case_417(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return tan(3 * x) / (tan(3 * x) - 1);
}


template <class T>
inline T test_case_418(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(sec(x));
}


template <class T>
inline T test_case_419(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (sec(x) + 1);
}


template <class T>
inline T test_case_420(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(3 * x) * cos(3 * x);
}


template <class T>
inline T test_case_421(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sec(x) * tan(x);
}


template <class T>
inline T test_case_422(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(x) * boost::math::pow<2>(cos(3 * boost::math::constants::pi<T>() * x / 2));
}


template <class T>
inline T test_case_423(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x / (boost::math::pow<2>(x) + 9) * log(boost::math::pow<2>(x) + 9);
}


template <class T>
inline T test_case_424(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / x / log(x);
}


template <class T>
inline T test_case_425(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / exp(x) * (1 / x - log(x));
}


template <class T>
inline T test_case_426(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return T(1) / T(5) * (T(1) / 100 * (322 + 3 * x * (98 + x * (37 + x))) - 24 * x / (1 + boost::math::pow<2>(x)));
}


template <class T>
inline T test_case_427(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(x) * boost::math::pow<2>(sin(3 * x));
}


template <class T>
inline T test_case_428(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * sin(3 * x);
}


template <class T>
inline T test_case_429(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 + cos(3 * x));
}


template <class T>
inline T test_case_430(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cos(3 * x) / (1 + cos(3 * x));
}


template <class T>
inline T test_case_431(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(tan(x));
}


template <class T>
inline T test_case_432(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<5>((log(1 / x)));
}


template <class T>
inline T test_case_433(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<8>(x) * boost::math::pow<4>((1 - x));
}


template <class T>
inline T test_case_434(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return acosh(4 * x);
}


template <class T>
inline T test_case_435(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (4 * x + 3) / (6 * boost::math::pow<2>(x) + 3 * x + 8);
}


template <class T>
inline T test_case_436(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cos(boost::math::constants::pi<T>() * boost::math::pow<2>(x) / 2);
}


template <class T>
inline T test_case_437(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(boost::math::constants::pi<T>() * boost::math::pow<2>(x) / 2);
}


template <class T>
inline T test_case_438(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (3 + 4 * exp(5 * x));
}


template <class T>
inline T test_case_439(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<4>(x) * exp(3 * x);
}


template <class T>
inline T test_case_440(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * exp(3 * x) / boost::math::pow<2>((1 + 3 * x));
}


template <class T>
inline T test_case_441(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<5>(sinh(x));
}


template <class T>
inline T test_case_442(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<5>(cosh(x));
}


template <class T>
inline T test_case_443(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<4>(sinh(x)) * boost::math::pow<2>(cosh(x));
}


template <class T>
inline T test_case_444(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<3>(sinh(x)) * boost::math::pow<3>(cosh(x));
}


template <class T>
inline T test_case_445(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sinh(x) * boost::math::pow<4>(cosh(x));
}


template <class T>
inline T test_case_446(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / sinh(x);
}


template <class T>
inline T test_case_447(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / cosh(x);
}


template <class T>
inline T test_case_448(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sinh(x) / cosh(x);
}


template <class T>
inline T test_case_449(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<3>(sinh(x)) / boost::math::pow<2>(cosh(x));
}


template <class T>
inline T test_case_450(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cosh(x) / sinh(x);
}


template <class T>
inline T test_case_451(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(cosh(x)) / sinh(x);
}


template <class T>
inline T test_case_452(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / sinh(x) / cosh(x);
}


template <class T>
inline T test_case_453(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sinh(2 * x) / boost::math::pow<2>(cosh(x));
}


template <class T>
inline T test_case_454(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sinh(x) / (cosh(x) + sinh(x));
}


template <class T>
inline T test_case_455(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (cosh(x) - sinh(x));
}


template <class T>
inline T test_case_456(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 - boost::math::pow<2>(cosh(x)));
}


template <class T>
inline T test_case_457(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * sinh(x);
}


template <class T>
inline T test_case_458(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * cosh(x);
}


template <class T>
inline T test_case_459(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x / boost::math::pow<2>(cosh(x));
}


template <class T>
inline T test_case_460(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x / (1 + cosh(x));
}


template <class T>
inline T test_case_461(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(sin(x)) * boost::math::pow<2>(cos(x));
}


template <class T>
inline T test_case_462(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(x) * boost::math::pow<4>(cos(x));
}


template <class T>
inline T test_case_463(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<4>(sin(x)) * boost::math::pow<4>(cos(x));
}


template <class T>
inline T test_case_464(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<3>(sin(x)) / boost::math::pow<2>(cos(x));
}


template <class T>
inline T test_case_465(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(x) / boost::math::pow<4>(cos(x));
}


template <class T>
inline T test_case_466(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cos(x) / boost::math::pow<3>(sin(x));
}


template <class T>
inline T test_case_467(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (boost::math::pow<4>(sin(x)) * boost::math::pow<4>(cos(x)));
}


template <class T>
inline T test_case_468(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cos(2 * x) / boost::math::pow<2>(cos(x));
}


template <class T>
inline T test_case_469(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(3 * x) / boost::math::pow<3>(sin(x));
}


template <class T>
inline T test_case_470(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cos(3 * x) / boost::math::pow<3>(cos(x));
}


template <class T>
inline T test_case_471(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cos(log(x));
}


template <class T>
inline T test_case_472(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cos(x) / pow(1 - 0.01 * boost::math::pow<2>(sin(x)), T(3) / 2);
}


template <class T>
inline T test_case_473(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cosh(x) * cos(x);
}


template <class T>
inline T test_case_474(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(x) * exp(-x * log(3));
}


template <class T>
inline T test_case_475(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * exp(-x * log(3));
}


template <class T>
inline T test_case_476(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return ((sin(2 * boost::math::constants::pi<T>() * x) + cos(boost::math::constants::pi<T>() * x)) / boost::math::constants::pi<T>() / (2 * x + 1)) * 2 * (-1);
}


template <class T>
inline T test_case_477(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (sin(boost::math::constants::pi<T>() * (x - 0.5)) - sin(2 * boost::math::constants::pi<T>() * (x - 0.5))) / boost::math::constants::pi<T>() / (x - 0.5);
}


template <class T>
inline T test_case_478(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 + exp(-x) * sin(8 * pow(x, T(2) / 3));
}


template <class T>
inline T test_case_479(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 3 * exp(-x) * sin(boost::math::pow<2>(x)) + 1;
}


template <class T>
inline T test_case_480(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(-2 * x) * (14 * x - 11 * boost::math::pow<2>(x));
}


template <class T>
inline T test_case_481(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<10>(x) * exp(4 * boost::math::pow<3>(x) - 3 * boost::math::pow<4>(x));
}


template <class T>
inline T test_case_482(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 2 + cos(2 * sqrt(x));
}


template <class T>
inline T test_case_483(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(boost::math::constants::pi<T>() * x) / boost::math::constants::pi<T>() / x;
}


template <class T>
inline T test_case_484(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<4>(x) * log(x + sqrt(boost::math::pow<2>(x) + 1));
}


template <class T>
inline T test_case_485(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * exp(x);
}


template <class T>
inline T test_case_486(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return atan(x);
}


template <class T>
inline T test_case_487(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(x) * cos(x);
}


template <class T>
inline T test_case_488(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::constants::pi<T>() / 2 * cosh(x) * sin(exp(boost::math::constants::pi<T>() / 2 * sinh(x)));
}


template <class T>
inline T test_case_489(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 4 * boost::math::pow<3>(x) * sqrt(boost::math::pow<4>(x) + 7);
}


template <class T>
inline T test_case_490(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * cos(x);
}


template <class T>
inline T test_case_491(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(cos(x));
}


template <class T>
inline T test_case_492(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(tan(x)) * boost::math::pow<4>(sec(x));
}


template <class T>
inline T test_case_493(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (3 * x + 6) / (boost::math::pow<2>(x) + 5 * x + 4);
}


template <class T>
inline T test_case_494(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x / (boost::math::pow<2>(x) + 4);
}


template <class T>
inline T test_case_495(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (boost::math::pow<2>(x) + 4);
}


template <class T>
inline T test_case_496(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x / boost::math::pow<2>((boost::math::pow<2>(x) + 4));
}


template <class T>
inline T test_case_497(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return tan(x);
}


template <class T>
inline T test_case_498(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sec(x);
}


template <class T>
inline T test_case_499(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(-(boost::math::pow<2>(x)));
}


template <class T>
inline T test_case_500(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>((2 * x + 3 / x));
}


template <class T>
inline T test_case_501(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return pow(x, 0.1) * (1.2 - x) * (1 - exp(20 * (x - 1)));
}


template <class T>
inline T test_case_502(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (9 + boost::math::pow<2>(x));
}


template <class T>
inline T test_case_503(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(abs(x - 0.5));
}


template <class T>
inline T test_case_504(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(1 - 30 * boost::math::pow<2>(x));
}


template <class T>
inline T test_case_505(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (2 + cos(1 + pow(x, T(3) / 2)) * exp(0.5 * x) / sqrt(1 + sin(x) / 2)) * exp(0.5 * x);
}


template <class T>
inline T test_case_506(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 + boost::math::pow<2>((x + 3)));
}


template <class T>
inline T test_case_507(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(-x) / (1 + boost::math::pow<2>(x));
}


template <class T>
inline T test_case_508(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(1 + x) * log(1 - x);
}


template <class T>
inline T test_case_509(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return - boost::math::pow<3>(x) / 10 + 23 * x - 3.5;
}


template <class T>
inline T test_case_510(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<25>(x) * boost::math::pow<2>((1 - x));
}


template <class T>
inline T test_case_511(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<10>(x) - 10 * boost::math::pow<8>(x) + 33 * boost::math::pow<6>(x) - 40 * boost::math::pow<4>(x) + 16 * boost::math::pow<2>(x);
}


template <class T>
inline T test_case_512(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(-x) * sin(x) / x;
}


template <class T>
inline T test_case_513(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (boost::math::pow<2>(x) * log(x)) / (boost::math::pow<2>(x) - 1) / (boost::math::pow<4>(x) + 1);
}


template <class T>
inline T test_case_514(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return asin(sqrt(2) / 2 * sin(x)) * sin(x) / sqrt(4 - 2 * boost::math::pow<2>(sin(x)));
}


template <class T>
inline T test_case_515(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x / (1 - exp(-6 * sinh(x)));
}


template <class T>
inline T test_case_516(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return tanh(boost::math::constants::pi<T>() / 2 * sinh(x));
}


template <class T>
inline T test_case_517(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<3>(x) + boost::math::pow<2>(x) + 5 * x + 3;
}


template <class T>
inline T test_case_518(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<3>(sin(x)) / (boost::math::pow<3>(sin(x)) + boost::math::pow<3>(cos(x)));
}


template <class T>
inline T test_case_519(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(-x) * sin(x);
}


template <class T>
inline T test_case_520(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(x + sqrt(1 + boost::math::pow<2>(x)));
}


template <class T>
inline T test_case_521(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / sqrt(exp(2 * x) + 1);
}


template <class T>
inline T test_case_522(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<8>(sec(x)) * tan(x);
}


template <class T>
inline T test_case_523(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(-x) / (1 + boost::math::pow<4>(x));
}


template <class T>
inline T test_case_524(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 8 * exp(-8 * x);
}


template <class T>
inline T test_case_525(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 + cos(1.95 * boost::math::constants::pi<T>() * x);
}


template <class T>
inline T test_case_526(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(-x) / sqrt(x);
}


template <class T>
inline T test_case_527(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(-x) * cos(x);
}


template <class T>
inline T test_case_528(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(-x) * boost::math::pow<2>(x);
}


template <class T>
inline T test_case_529(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 100 * exp(-200 * boost::math::pow<2>(x));
}


template <class T>
inline T test_case_530(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (14 * x - 11 * boost::math::pow<2>(x)) * exp(-2 * x);
}


template <class T>
inline T test_case_531(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 2000 * log(140000 / (140000 - 2100 * x)) - 9.8 * x;
}


template <class T>
inline T test_case_532(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (boost::math::pow<2>(x) + boost::math::pow<-6>(10.0));
}


template <class T>
inline T test_case_533(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 + exp(-25 * x) / 2;
}


template <class T>
inline T test_case_534(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / boost::math::pow<3>(log(1 / x)) / x;
}


template <class T>
inline T test_case_535(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(3 * x) * sin(2 * x);
}


template <class T>
inline T test_case_536(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * exp(x) / boost::math::pow<2>((x + 1));
}


template <class T>
inline T test_case_537(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 5 / (exp(boost::math::constants::pi<T>()) - 2) * exp(2 * x) * cos(x);
}


template <class T>
inline T test_case_538(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return - (boost::math::pow<3>(x) + 23 * x - 3.5);
}


template <class T>
inline T test_case_539(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 10 * exp(-100 * boost::math::pow<2>(x));
}


template <class T>
inline T test_case_540(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<3>(x);
}


template <class T>
inline T test_case_541(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<6>(x);
}


template <class T>
inline T test_case_542(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<10>(x);
}


template <class T>
inline T test_case_543(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(boost::math::pow<7>(x));
}


template <class T>
inline T test_case_544(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 + boost::math::pow<2>((sin(x))));
}


template <class T>
inline T test_case_545(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (x - 1) * (x - 2) * (x - 3) * (x - 4) * (x - 5) / 120;
}


template <class T>
inline T test_case_546(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / boost::math::pow<5>((x + T(1) / 100));
}


template <class T>
inline T test_case_547(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 50 / (boost::math::constants::pi<T>() * (1 + 2500 * boost::math::pow<2>(x)));
}


template <class T>
inline T test_case_548(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (sin(100 * boost::math::constants::pi<T>() * x) / (boost::math::constants::pi<T>() * x));
}


template <class T>
inline T test_case_549(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * cos(x) * sin(3 * x);
}


template <class T>
inline T test_case_550(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(x) * log(1 - x);
}


template <class T>
inline T test_case_551(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * sin(30 * x) / sqrt(1 - boost::math::pow<2>(x) / (4 * boost::math::pow<2>(boost::math::constants::pi<T>())));
}


template <class T>
inline T test_case_552(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / sqrt(x) * log(log(1 / x));
}


template <class T>
inline T test_case_553(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / sqrt(x) / sqrt(1 - x);
}


template <class T>
inline T test_case_554(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 - log(x)) / sqrt(-log(x));
}


template <class T>
inline T test_case_555(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return pow(x, 0.95) * exp(x);
}


template <class T>
inline T test_case_556(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(log(x)) * 1 / (1 + boost::math::pow<2>(x));
}


template <class T>
inline T test_case_557(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(x) * exp(-x) / (1 + x);
}


template <class T>
inline T test_case_558(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(abs(3 + x));
}


template <class T>
inline T test_case_559(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<3>(x) * boost::math::pow<5>((1 - x));
}


template <class T>
inline T test_case_560(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(x) / sqrt(1 - x);
}


template <class T>
inline T test_case_561(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / sqrt(x) / sqrt(1 - x);
}


template <class T>
inline T test_case_562(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<8>(x) / sqrt(1 - x);
}


template <class T>
inline T test_case_563(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(x) / (1 + boost::math::pow<4>(x)) / sqrt(1 - boost::math::pow<4>(x));
}


template <class T>
inline T test_case_564(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * exp(-3 * x);
}


template <class T>
inline T test_case_565(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * exp(x) / boost::math::pow<2>((1 + x));
}


template <class T>
inline T test_case_566(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(-4 * x) / sqrt(x);
}


template <class T>
inline T test_case_567(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x / (exp(x) + 2 * exp(-x) - 2);
}


template <class T>
inline T test_case_568(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x / (1 - exp(-x));
}


template <class T>
inline T test_case_569(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (4 * exp(1 - boost::math::pow<-4>(x)) / (1 - boost::math::pow<4>(x)) - exp(1 - 1 / x) / (1 - x)) / x;
}


template <class T>
inline T test_case_570(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(-5 / x) / boost::math::pow<2>(x);
}


template <class T>
inline T test_case_571(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(-(boost::math::pow<2>(tan(x))));
}


template <class T>
inline T test_case_572(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x / sqrt(cosh(8) - cosh(8 * x)) / sinh(4 * x);
}


template <class T>
inline T test_case_573(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(7 * x) / sin(x);
}


template <class T>
inline T test_case_574(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(8 * x) * cos(5 * x) / sin(x);
}


template <class T>
inline T test_case_575(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(5 * x) / sin(x);
}


template <class T>
inline T test_case_576(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cos(7 * x) / cos(x);
}


template <class T>
inline T test_case_577(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cos(2 * x) / (1 - 6 * cos(x) + 9);
}


template <class T>
inline T test_case_578(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(4 * x) * sin(x) / (1 - 10 * cos(x) + 25);
}


template <class T>
inline T test_case_579(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<6>(sin(x)) / boost::math::pow<8>(cos(x));
}


template <class T>
inline T test_case_580(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(10 * x) * boost::math::pow<9>(cos(x)) / sin(x);
}


template <class T>
inline T test_case_581(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(sin(x));
}


template <class T>
inline T test_case_582(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / sqrt(sin(x));
}


template <class T>
inline T test_case_583(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<5>(cos(x)) * cos(5 * x) / (16 * boost::math::pow<2>(sin(x)) + 9 * boost::math::pow<2>(cos(x)));
}


template <class T>
inline T test_case_584(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(sin(x)) / (5 + 4 * cos(x));
}


template <class T>
inline T test_case_585(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(sin(x)) / (1 / sqrt((cos(x) - sin(x)))) / boost::math::pow<3>(cos(x));
}


template <class T>
inline T test_case_586(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(x) / sqrt(1 + 25 * boost::math::pow<2>(sin(x)));
}


template <class T>
inline T test_case_587(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(x) / sqrt(1 - 8 * cos(x) + 16);
}


template <class T>
inline T test_case_588(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt((cos(2 * x) - cos(2)) / (cos(2 * x) + 1));
}


template <class T>
inline T test_case_589(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (boost::math::pow<3>(tan(x)) + boost::math::pow<3>(cot(x))) / sin(2 * x);
}


template <class T>
inline T test_case_590(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return tan(x) / (boost::math::pow<4>(cos(x)) + boost::math::pow<4>(sec(x)));
}


template <class T>
inline T test_case_591(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(4 * sin(x)) * sin(2 * x);
}


template <class T>
inline T test_case_592(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(4 * cos(x)) * sin(2 * x);
}


template <class T>
inline T test_case_593(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(4 * tan(x)) * sin(2 * x);
}


template <class T>
inline T test_case_594(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * tan(x);
}


template <class T>
inline T test_case_595(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * cot(x);
}


template <class T>
inline T test_case_596(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x / sin(x);
}


template <class T>
inline T test_case_597(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * cot(x);
}


template <class T>
inline T test_case_598(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (boost::math::constants::pi<T>() / 2 - x) * tan(x);
}


template <class T>
inline T test_case_599(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(x) / boost::math::pow<2>(cos(x));
}


template <class T>
inline T test_case_600(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * boost::math::pow<3>(tan(x));
}


template <class T>
inline T test_case_601(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(x) * tan(x) / boost::math::pow<2>(cos(x));
}


template <class T>
inline T test_case_602(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(x) * boost::math::pow<2>(tan(x)) / boost::math::pow<2>(cos(x));
}


template <class T>
inline T test_case_603(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / sqrt(log(1 / x));
}


template <class T>
inline T test_case_604(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(x) / boost::math::pow<2>((1 + log(x)));
}


template <class T>
inline T test_case_605(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(x) * log(1 - x);
}


template <class T>
inline T test_case_606(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<4>((log(tan(x))));
}


template <class T>
inline T test_case_607(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(tan(x));
}


template <class T>
inline T test_case_608(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(1 + tan(x));
}


template <class T>
inline T test_case_609(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(sqrt(tan(x)) + sqrt(cot(x)));
}


template <class T>
inline T test_case_610(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return - log(cos(boost::math::constants::pi<T>() / 2 * x));
}


template <class T>
inline T test_case_611(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(x) / sqrt(1 + 9 * boost::math::pow<2>(sin(x)));
}


template <class T>
inline T test_case_612(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(x) / sqrt(1 - 6 * cos(x) + 9);
}


template <class T>
inline T test_case_613(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt((cos(2 * x) - cos(2)) / (cos(2 * x) + 1));
}


template <class T>
inline T test_case_614(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(tan(x)) + sqrt(cot(x));
}


template <class T>
inline T test_case_615(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (sqrt(tan(x)) - sqrt(cot(x))) * tan(x);
}


template <class T>
inline T test_case_616(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (boost::math::pow<3>(tan(x)) + boost::math::pow<3>(cot(x))) / sin(2 * x);
}


template <class T>
inline T test_case_617(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(5 * cos(x)) * sin(2 * x);
}


template <class T>
inline T test_case_618(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / sqrt(1 - 2 * x * cos(3) + boost::math::pow<2>(x));
}


template <class T>
inline T test_case_619(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (cos(x) - sin(x)) / (cos(x) + sin(x)) * x;
}


template <class T>
inline T test_case_620(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (boost::math::constants::pi<T>() / 4 - x * tan(x)) * tan(x);
}


template <class T>
inline T test_case_621(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (boost::math::constants::pi<T>() / 4 - x) * tan(x) / cos(2 * x);
}


template <class T>
inline T test_case_622(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (boost::math::constants::pi<T>() / 4 - x * tan(x)) / cos(2 * x);
}


template <class T>
inline T test_case_623(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x / sin(x) / (cos(x) + sin(x));
}


template <class T>
inline T test_case_624(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(x) * x / (sin(x) + cos(x)) / boost::math::pow<2>(cos(x));
}


template <class T>
inline T test_case_625(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / x - cot(x);
}


template <class T>
inline T test_case_626(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (4 * boost::math::pow<2>(x) * cos(x) + (boost::math::constants::pi<T>() - x) * x) / sin(x);
}


template <class T>
inline T test_case_627(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * cos(x) / (1 + sin(x));
}


template <class T>
inline T test_case_628(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * sin(x) / (1 - cos(x));
}


template <class T>
inline T test_case_629(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x / (sin(x) + cos(x)) / sin(x);
}


template <class T>
inline T test_case_630(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * sin(2 * x) / (9 * boost::math::pow<2>(cos(x)) + 16 * boost::math::pow<2>(sin(x)));
}


template <class T>
inline T test_case_631(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * sin(2 * x) / (1 + 3 * boost::math::pow<2>(sin(x))) / (1 + 4 * boost::math::pow<2>(sin(x)));
}


template <class T>
inline T test_case_632(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(x) * tan(x) / boost::math::pow<2>(cos(x));
}


template <class T>
inline T test_case_633(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<3>(x) * cos(x) / boost::math::pow<3>(sin(x));
}


template <class T>
inline T test_case_634(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(3 * cos(x)) * sin(3 * sin(x)) * cot(x / 2);
}


template <class T>
inline T test_case_635(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(cos(x) - sin(x));
}


template <class T>
inline T test_case_636(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(cos(x) + sin(x));
}


template <class T>
inline T test_case_637(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>((log(tan(x))));
}


template <class T>
inline T test_case_638(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<4>((log(tan(x))));
}


template <class T>
inline T test_case_639(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(log(cos(x)));
}


template <class T>
inline T test_case_640(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(cos(x) + sin(x));
}


template <class T>
inline T test_case_641(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * sin(30 * x) / sqrt(1 - boost::math::pow<2>((x / 2 / boost::math::constants::pi<T>())));
}


template <class T>
inline T test_case_642(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cos(sqrt(377) * x);
}


template <class T>
inline T test_case_643(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(x) * exp(-(boost::math::pow<2>(x))) * tan(x) * acos(x);
}


template <class T>
inline T test_case_644(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(-3 * x) * cos(16 * sqrt(3) * boost::math::constants::pi<T>() * x);
}


template <class T>
inline T test_case_645(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(cos(sqrt(47 * boost::math::constants::pi<T>()) * x));
}


template <class T>
inline T test_case_646(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cosh(tanh(sinh(x)));
}


template <class T>
inline T test_case_647(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return T(23) / T(25) * cosh(x) - cos(x);
}


template <class T>
inline T test_case_648(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(-2 * x) * cos(16 * sqrt(2) * x);
}


template <class T>
inline T test_case_649(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * atan(boost::math::pow<3>(x));
}


template <class T>
inline T test_case_650(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(x) * atan(boost::math::pow<3>(x));
}


template <class T>
inline T test_case_651(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * sin(30 * x) / sqrt(1 - boost::math::pow<2>(x) / 4 / boost::math::pow<2>(boost::math::constants::pi<T>()));
}


template <class T>
inline T test_case_652(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * sin(50 * x) * cos(75 * x);
}


template <class T>
inline T test_case_653(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (boost::math::pow<4>(x) + boost::math::pow<2>(x) + exp(1));
}


template <class T>
inline T test_case_654(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return tan(x) / (1 + exp(x) * sin(boost::math::constants::pi<T>() * x));
}


template <class T>
inline T test_case_655(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(-3 * x) * cos(5 * boost::math::constants::pi<T>() * x);
}


template <class T>
inline T test_case_656(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(x) * sin(3 * x) * tanh(5 * cos(30 * x));
}


template <class T>
inline T test_case_657(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return abs(x - 1 / sqrt(3)) + abs(x - 1 / sqrt(2));
}


template <class T>
inline T test_case_658(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return pow(x, T(-2) / 3);
}


template <class T>
inline T test_case_659(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(10 * x) / sqrt(x * (1 - x));
}


template <class T>
inline T test_case_660(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<3>(x) * sin(x) / sqrt(1 - boost::math::pow<2>(x));
}


template <class T>
inline T test_case_661(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<5>(x) * log(1 / x);
}


template <class T>
inline T test_case_662(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   // x^(-x)^((-3/10)*(-LN(x))^(-7/10))
   return pow(x,  pow(-x, ((T(-3) / 10) * pow(-log(x), T(-7) / 10))));
}


template <class T>
inline T test_case_663(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(x) / (100 * pow(x, 2 * x) + 1);
}


template <class T>
inline T test_case_664(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * (-log(x)) / (100 * pow(x, 2 * x) + 1);
}


template <class T>
inline T test_case_665(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(x * (-log(x)) / (100 * pow(x, -2 * x) + 1));
}


template <class T>
inline T test_case_666(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(100 * x * (-log(x))) / (100 * pow(x, -2 * x) + 1);
}


template <class T>
inline T test_case_667(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(70 * boost::math::pow<2>(x)) / (100 * pow(x, x) + 1);
}


template <class T>
inline T test_case_668(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(70 * sin(x * (-log(x)))) / (100 * pow(x, x) + 1);
}


template <class T>
inline T test_case_669(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(pow(sqrt(x) * (-log(x)), T(1) / 10)) / (pow(x, T(1) / 10) * pow(-log(x), T(1) / 5));
}


template <class T>
inline T test_case_670(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(30 * pow(x, T(1) / 20) * pow(-log(x), T(1) / 10) / (100 * pow(x, T(1) / 10) * pow(-log(x), T(1) / 5) + 1));
}


template <class T>
inline T test_case_671(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(x) * exp(x);
}


template <class T>
inline T test_case_672(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<5>((-log(x))) * boost::math::pow<4>(x);
}


template <class T>
inline T test_case_673(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<6>(sin(x)) * boost::math::pow<6>(cos(x));
}


template <class T>
inline T test_case_674(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(log(sin(x)));
}


template <class T>
inline T test_case_675(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(-log(x));
}


template <class T>
inline T test_case_676(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cos(x * exp(x));
}


template <class T>
inline T test_case_677(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<9>((1 - sqrt(x)));
}


template <class T>
inline T test_case_678(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(x) * pow(1 - x, 0.3);
}


template <class T>
inline T test_case_679(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 3 * boost::math::pow<2>(x) / (boost::math::pow<6>(x) + 1);
}


template <class T>
inline T test_case_680(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cos(4 * x) / sqrt(x);
}


template <class T>
inline T test_case_681(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (cos(x) - cos(2 * x)) / x;
}


template <class T>
inline T test_case_682(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return asin(sqrt(2) / 2 * sin(x)) * sin(x) / sqrt(4 - 2 * boost::math::pow<2>(sin(x)));
}


template <class T>
inline T test_case_683(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::constants::pi<T>() * cosh(x) * sin(exp(boost::math::constants::pi<T>() / 2 * sinh(x))) / 2;
}


template <class T>
inline T test_case_684(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<6>(log(x)) * atan(abs(x * sqrt(3) / (x - 2))) / (x + 1);
}


template <class T>
inline T test_case_685(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(1 - 1 / x) * cos(1 / x - 1) / boost::math::pow<2>(x);
}


template <class T>
inline T test_case_686(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(x) + exp(x) * sin(x);
}


template <class T>
inline T test_case_687(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 3 * boost::math::pow<2>(x) * boost::math::pow<3>(sin(100 * x)) * exp(1 / 3 * x);
}


template <class T>
inline T test_case_688(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(cos(x));
}


template <class T>
inline T test_case_689(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(1 - 0.36 * boost::math::pow<2>(sin(x)));
}


template <class T>
inline T test_case_690(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (1 + boost::math::pow<2>(tan(x / 2))) / (1 + boost::math::pow<4>(tan(x / 2)));
}


template <class T>
inline T test_case_691(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(-(boost::math::pow<2>(tan(x / 2))));
}


template <class T>
inline T test_case_692(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (1 + x) / boost::math::cbrt(boost::math::pow<2>(x) + 2 * x + 5);
}


template <class T>
inline T test_case_693(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / sqrt(x) * log(exp(1) / x);
}


template <class T>
inline T test_case_694(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * exp(boost::math::pow<2>(x)) * erf(x) / sqrt(1 - boost::math::pow<2>(x));
}


template <class T>
inline T test_case_695(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<3>(x) * boost::math::pow<4>(log(1 / x));
}


template <class T>
inline T test_case_696(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 - boost::math::cbrt(boost::math::pow<2>((x - boost::math::constants::pi<T>() / 2 / exp(1))));
}


template <class T>
inline T test_case_697(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (sqrt(x) + boost::math::cbrt(x));
}


template <class T>
inline T test_case_698(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(x) / sqrt(x + 0.01);
}


template <class T>
inline T test_case_699(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / sqrt(x) / (1 + x);
}


template <class T>
inline T test_case_700(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(-3 * boost::math::pow<2>(x)) * log(1 + x + boost::math::pow<2>(x));
}


template <class T>
inline T test_case_701(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(x) * log(x) / sqrt(1 - boost::math::pow<2>(x));
}


template <class T>
inline T test_case_702(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(4095 * boost::math::pow<2>(cos(x)) + 1);
}


template <class T>
inline T test_case_703(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<6>(x) - T(105) / T(4) * boost::math::pow<4>(x) + T(315) / T(2) * boost::math::pow<2>(x) - T(315) / T(4);
}


template <class T>
inline T test_case_704(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<8>(x) - T(104) / T(3) * boost::math::pow<6>(x) + 658 * boost::math::pow<4>(x) - 2940 * boost::math::pow<2>(x) + 1785;
}


template <class T>
inline T test_case_705(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (boost::math::pow<2>(x) + x + 1) * cos(x);
}


template <class T>
inline T test_case_706(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return pow(x, T(2) / 3) / pow(boost::math::pow<2>(x) + boost::math::pow<2>((1 - x)), T(4) / 3);
}


template <class T>
inline T test_case_707(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(x) / exp(boost::math::pow<2>(x));
}


template <class T>
inline T test_case_708(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return atan(x) / pow(x, T(3) / 2);
}


template <class T>
inline T test_case_709(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return atan(x) / sqrt(x);
}


template <class T>
inline T test_case_710(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / x / exp(boost::math::pow<2>((-log(x))));
}


template <class T>
inline T test_case_711(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * abs(x);
}


template <class T>
inline T test_case_712(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x / (5 + x) * exp(-2 * x / 30);
}


template <class T>
inline T test_case_713(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(x) / boost::math::pow<3>(x);
}


template <class T>
inline T test_case_714(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(-(boost::math::pow<2>(cos(x)))) / sqrt(1 - boost::math::pow<2>(x));
}


template <class T>
inline T test_case_715(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 12 * x * exp(-2 * x);
}


template <class T>
inline T test_case_716(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(x) * (2 + sin(1 / x));
}


template <class T>
inline T test_case_717(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(2 * sin(2 * sin(2 * sin(x))));
}


template <class T>
inline T test_case_718(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return - 2 / sqrt(boost::math::constants::pi<T>()) * log(cos(boost::math::constants::pi<T>() * tanh(x) / 2) / boost::math::pow<2>(cosh(x)));
}


template <class T>
inline T test_case_719(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(1 + 1 / x);
}


template <class T>
inline T test_case_720(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return - log(log(2 / (x + 1)));
}


template <class T>
inline T test_case_721(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * sin(x) / (1 + boost::math::pow<2>(cos(x)));
}


template <class T>
inline T test_case_722(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(x) * log(x) / (boost::math::pow<2>(x) - 1) / (boost::math::pow<4>(x) + 1);
}


template <class T>
inline T test_case_723(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin((x - 4) * 55) / (x - 4.99);
}


template <class T>
inline T test_case_724(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(-0.01 * x) * (0.01 * cos(0.3 * x) + 0.3 * sin(0.3 * x));
}


template <class T>
inline T test_case_725(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (0.2 * (1 + x) * cos(0.2 * x) - sin(0.2 * x)) / boost::math::pow<2>((1 + x));
}


template <class T>
inline T test_case_726(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 80 * sin(sqrt(1 + 80 * x)) / 2 / sqrt(1 + 80 * x);
}


template <class T>
inline T test_case_727(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(-sqrt(1 + x)) * (0.5 * cos(0.5 * x) - sin(0.5 * x) / 2 / sqrt(1 + x));
}


template <class T>
inline T test_case_728(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(-2 / (1 + x) - 2 / (1 - x));
}


template <class T>
inline T test_case_729(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / pow(1 + boost::math::pow<2>(x), T(5) / 4);
}


template <class T>
inline T test_case_730(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 3 * boost::math::pow<4>(x) * pow(boost::math::pow<6>(x) + boost::math::pow<2>((1 - boost::math::pow<3>(x))), T(-4) / 3);
}


template <class T>
inline T test_case_731(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(21 * boost::math::constants::pi<T>() * x) + sin(31 * boost::math::constants::pi<T>() * x) / 2;
}


template <class T>
inline T test_case_732(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(100 * x) / x;
}


template <class T>
inline T test_case_733(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (cos(x) - 1) / x;
}


template <class T>
inline T test_case_734(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return csch(x);
}


template <class T>
inline T test_case_735(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sech(x);
}


template <class T>
inline T test_case_736(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return coth(x);
}


template <class T>
inline T test_case_737(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return csc(x) * cot(x);
}


template <class T>
inline T test_case_738(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return tanh(x);
}


template <class T>
inline T test_case_739(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<5>(log(1 / x));
}


template <class T>
inline T test_case_740(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return pow(x, -x);
}


template <class T>
inline T test_case_741(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return pow(x, x);
}


template <class T>
inline T test_case_742(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(boost::math::pow<4>(x)) / x;
}


template <class T>
inline T test_case_743(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (x + 3) * boost::math::pow<2>((x - 1));
}


template <class T>
inline T test_case_744(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (1 - cos(x)) / boost::math::pow<2>(x);
}


template <class T>
inline T test_case_745(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(1 / x) / pow(x, 0.25);
}


template <class T>
inline T test_case_746(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (16 * boost::math::pow<2>((x - boost::math::constants::pi<T>() / 4)) + T(1) / 16);
}


template <class T>
inline T test_case_747(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cos(64 * sin(x));
}


template <class T>
inline T test_case_748(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(20 * (x - 1)) * sin(256 * x);
}


template <class T>
inline T test_case_749(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * sin(2 * exp(2 * sin(2 * exp(2 * x))));
}


template <class T>
inline T test_case_750(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return abs(sin(x) / (2 * ((pow(x, x) - (boost::math::constants::pi<T>() / 2)) / boost::math::constants::pi<T>())));
}


template <class T>
inline T test_case_751(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(2 * x) * sin(3 * x);
}


template <class T>
inline T test_case_752(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 2 * x * cos(2 * x) - boost::math::pow<2>((x - 2));
}


template <class T>
inline T test_case_753(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(x / (1 - x * x));
}


template <class T>
inline T test_case_754(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(x) * log(1 - x);
}


template <class T>
inline T test_case_755(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 + 25 * x * x);
}


template <class T>
inline T test_case_756(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 + 0.04 * x * x);
}


template <class T>
inline T test_case_757(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(x) * cos(10 * boost::math::constants::pi<T>() * x);
}


template <class T>
inline T test_case_758(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(sech(10 * (x - 0.2))) + boost::math::pow<4>(sech(100 * (x - 0.4))) + boost::math::pow<6>(sech(1000 * (x - 0.6))) + boost::math::pow<8>(sech(1000 * (x - 0.8)));
}


template <class T>
inline T test_case_759(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cos(10 * x) * boost::math::tgamma(x + 2) * erf(sqrt(1 + x));
}


template <class T>
inline T test_case_760(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * boost::math::pow<2>(cos(20 * x));
}


template <class T>
inline T test_case_761(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(abs(x + 0.5));
}


template <class T>
inline T test_case_762(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (1 + 1000 * boost::math::pow<2>(x));
}


template <class T>
inline T test_case_763(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return pow(1 - x, 0.6) * 1 / sqrt(1 + x);
}


template <class T>
inline T test_case_764(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return pow(1 - x, -0.3) * pow(1 + x, 0.2);
}


template <class T>
inline T test_case_765(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(1 / x) * asin(x);
}


template <class T>
inline T test_case_766(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (boost::math::pow<2>((x - sqrt(3))) + 0.0001);
}


template <class T>
inline T test_case_767(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sin(log(x));
}


template <class T>
inline T test_case_768(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(x) * log(x) / sqrt(1 - boost::math::pow<2>(x));
}


template <class T>
inline T test_case_769(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return asinh(x * (1 + x / 2) / (1 + x));
}


template <class T>
inline T test_case_770(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (x + boost::math::pow<-4>(10.));
}


template <class T>
inline T test_case_771(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(1 / x) * asin(x);
}


template <class T>
inline T test_case_772(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x * boost::math::pow<2>(cos(20 * x));
}


template <class T>
inline T test_case_773(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return (sin(1 / x) + 1 + x * cos(1 / x)) / x;
}


template <class T>
inline T test_case_774(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<20>(cos(x)) * cos(20 * x);
}


template <class T>
inline T test_case_775(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(10 * cos(x)) * cos(10 * sin(x));
}


template <class T>
inline T test_case_776(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(-x) / sqrt(x);
}


template <class T>
inline T test_case_777(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(1 - boost::math::pow<2>(x)) / x / (1 + x);
}


template <class T>
inline T test_case_778(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(-(boost::math::pow<2>((1 / x - 1))) / 2) / boost::math::pow<2>(x);
}


template <class T>
inline T test_case_779(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(1 / x) * pow(x, -0.25);
}


template <class T>
inline T test_case_780(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 2 * exp(-(boost::math::pow<2>(x)));
}


template <class T>
inline T test_case_781(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(1 + boost::math::pow<2>(x)) / x;
}


template <class T>
inline T test_case_782(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(-(boost::math::pow<2>(x))) / x;
}


template <class T>
inline T test_case_783(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(1 + exp(-x)) / pow(x, 1.5);
}


template <class T>
inline T test_case_784(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(-(boost::math::pow<2>(x))) / pow(sin(x), 0.7);
}


template <class T>
inline T test_case_785(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(1 + x + exp(x)) / (pow(x, 0.2) + 2);
}


template <class T>
inline T test_case_786(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return asin(sqrt(2) / 2 * sin(x)) * sin(x) / sqrt(4 - 2 * boost::math::pow<2>(sin(x)));
}


template <class T>
inline T test_case_787(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(x) * atan(x);
}


template <class T>
inline T test_case_788(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(x) * cos(x);
}


template <class T>
inline T test_case_789(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return atan(sqrt(2 + boost::math::pow<2>(x))) / (1 + boost::math::pow<2>(x)) / sqrt(2 + boost::math::pow<2>(x));
}


template <class T>
inline T test_case_790(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(x) * log(x);
}


template <class T>
inline T test_case_791(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(1 - boost::math::pow<2>(x));
}


template <class T>
inline T test_case_792(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(x) / sqrt(1 - boost::math::pow<2>(x));
}


template <class T>
inline T test_case_793(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return boost::math::pow<2>(log(x));
}


template <class T>
inline T test_case_794(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(cos(x));
}


template <class T>
inline T test_case_795(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(tan(x));
}


template <class T>
inline T test_case_796(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(x + 1) / (boost::math::pow<2>(x) + 1);
}


template <class T>
inline T test_case_797(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(x * (4 - x));
}


template <class T>
inline T test_case_798(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt((1 - boost::math::pow<2>(x)) * (2 - x));
}


template <class T>
inline T test_case_799(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(abs(x - 1));
}


template <class T>
inline T test_case_800(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return pow(x, T(-1) / 4) * log(1 / x);
}


template <class T>
inline T test_case_801(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 0.5 * 1 / sqrt(x) * exp(-(sqrt(x)));
}


template <class T>
inline T test_case_802(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1.5 * sqrt(x) * exp(-(pow(x, 1.5)));
}


template <class T>
inline T test_case_803(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / boost::math::constants::pi<T>() / sqrt(1 - boost::math::pow<2>(x));
}


template <class T>
inline T test_case_804(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 2 / boost::math::constants::pi<T>() * sqrt(1 - boost::math::pow<2>(x));
}


template <class T>
inline T test_case_805(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 2 / boost::math::constants::pi<T>() / (1 + boost::math::pow<2>(x));
}


template <class T>
inline T test_case_806(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 2 / boost::math::constants::pi<T>() / (1 + boost::math::pow<2>(x));
}


template <class T>
inline T test_case_807(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 1 / (x - 2);
}


template <class T>
inline T test_case_808(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return T(1) / T(4) / log(2) * log((1 + x) / (1 - x));
}


template <class T>
inline T test_case_809(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(1 / x) / pow(x, 0.25);
}


template <class T>
inline T test_case_810(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return sqrt(tan(x));
}


template <class T>
inline T test_case_811(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(1 + exp(-x)) / pow(x, 1.5);
}


template <class T>
inline T test_case_812(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return log(exp(x) + x + 1) / (pow(x, 0.2) + 2);
}


template <class T>
inline T test_case_813(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x / (1 + boost::math::pow<6>(x) * boost::math::pow<2>(sin(x)));
}


template <class T>
inline T test_case_814(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return pow(x, -0.12345);
}


template <class T>
inline T test_case_815(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return cos(x) * pow(x, -0.12345);
}


template <class T>
inline T test_case_816(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return exp(-1 / x - 1 / (1 - x));
}


template <class T>
inline T test_case_817(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return 2 * (1 - boost::math::pow<2>(x)) / (cos(4 * atanh(x)) + cosh(2));
}


template <class T>
inline T test_case_818(const T& x)
{
   BOOST_MATH_STD_USING
   log_test_call(x);
   return x / (1 + boost::math::pow<6>(x) * boost::math::pow<2>((sinh(x))));
}


