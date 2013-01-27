// Copyright Christopher Kormanyos 2013.
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or
// copy at http://www.boost.org/LICENSE_1_0.txt).

#ifdef _MSC_VER
#  pragma warning (disable : 4512) // assignment operator could not be generated.
#  pragma warning (disable : 4996) // assignment operator could not be generated.
#endif

#include <iostream>
#include <limits>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <iterator>

#include <boost/math/special_functions/bessel.hpp>

/*
using

template <class output_iterator, class T>
inline void cyl_bessel_j_zero(
  output_iterator out_it,
  T v,
   unsigned number_of_zeros,
   unsigned start_index)
*/

#include <boost/multiprecision/cpp_dec_float.hpp>

namespace
{
//  typedef boost::multiprecision::cpp_dec_float_50 float_type;
 // typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<50>, boost::multiprecision::et_off> float_type;
  typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<50>, boost::multiprecision::et_on> float_type;
//  typedef float float_type;
}

int main()
{
  // Generate three double roots of Jv for integral order 2.

  std::vector<double> roots_3(3U);
  boost::math::cyl_bessel_j_zero(roots_3.begin(), 2, 3U, 1U);

  std::cout << "1st root = " << roots_3[0] << std::endl;
  std::cout << "2nd root = " << roots_3[1] << std::endl;
  std::cout << "3rd root = " << roots_3[2] << std::endl;



  // Generate 20 roots of Jv for order v=71/19.
  std::vector<float_type> roots(20U);

  boost::math::cyl_bessel_j_zero(roots.begin(), float_type(71) / 19, unsigned(roots.size()), 1U);

  // Set the precision of the output stream.
  std::cout.precision(std::numeric_limits<float_type>::digits10);

  // Print the roots to the output stream.
  std::copy(roots.begin(),
            roots.end(),
            std::ostream_iterator<float_type>(std::cout, "\n"));
}

/*
Mathematica: Table[N[BesselJZero[71/19, n], 50], {n, 1, 20, 1}]

7.2731751938316489503185694262290765588963196701623, \
10.724858308883141732536172745851416647110749599085, \
14.018504599452388106120459558042660282427471931581, \
17.252498459170417182162487166549777349195903838610, \
20.456678874044517595180234083894285885460502077814, \
23.643630897142345224945514227147319599854051725040, \
26.819671140255087745421311470965019261522390519297, \
29.988343117423674742679141796661432043878868194142, \
33.151796897690520871250862469973445265444791966114, \
36.311416000216207415724354035039386081316520184200, \
39.468132467505236587945197808083337887765967032029, \
42.622597801391236474855034831297954018844433480227, \
45.775281464536847753390206207806726581495950012439, \
48.926530489173566198367766817478553992471739894799, \
52.076607045343002794279746041878924876873478063472, \
55.225712944912571393594224327817265689059002890192, \
58.374006101538886436775188150439025201735151418932, \
61.521611873000965273726742659353136266390944103571, \
64.668631053790930368346482214873660794565966287160, \
67.815145619696290925556791375555951165111460585458

7.2731751938316489503185694262290765588963196701623
10.724858308883141732536172745851416647110749599085
14.018504599452388106120459558042660282427471931581
17.25249845917041718216248716654977734919590383861
20.456678874044517595180234083894285885460502077814
23.64363089714234522494551422714731959985405172504
26.819671140255087745421311470965019261522390519297
29.988343117423674742679141796661432043878868194142
33.151796897690520871250862469973445265444791966114
36.3114160002162074157243540350393860813165201842
39.468132467505236587945197808083337887765967032029
42.622597801391236474855034831297954018844433480227
45.775281464536847753390206207806726581495950012439
48.926530489173566198367766817478553992471739894799
52.076607045343002794279746041878924876873478063472
55.225712944912571393594224327817265689059002890192
58.374006101538886436775188150439025201735151418932
61.521611873000965273726742659353136266390944103571
64.66863105379093036834648221487366079456596628716
67.815145619696290925556791375555951165111460585458
*/


/*

------ Rebuild All started: Project: bessel_zeros_example, Configuration: Debug Win32 ------
Build started 27-Jan-2013 17:55:44.
_PrepareForClean:
  Deleting file "Debug\bessel_zeros_example.lastbuildstate".
InitializeBuildStatus:
  Creating "Debug\bessel_zeros_example.unsuccessfulbuild" because "AlwaysCreate" was specified.
ClCompile:
  bessel_zeros_example.cpp
I:\boost-trunk\boost/math/special_functions/bessel.hpp(418): warning C4244: '=' : conversion from 'value_type' to 'double', possible loss of data
          I:\boost-trunk\boost/math/special_functions/bessel.hpp(610) : see reference to function template instantiation 'void boost::math::detail::cyl_bessel_j_zero_imp<output_iterator,value_type,boost::math::policies::policy<boost::math::policies::default_policy,boost::math::policies::default_policy,boost::math::policies::default_policy,boost::math::policies::default_policy,boost::math::policies::default_policy,boost::math::policies::default_policy,boost::math::policies::default_policy,boost::math::policies::default_policy,boost::math::policies::default_policy,boost::math::policies::default_policy,boost::math::policies::default_policy>>(output_iterator,T,unsigned int,unsigned int,const Policy &)' being compiled
          with
          [
              output_iterator=std::_Vector_iterator<std::_Vector_val<double,std::allocator<double>>>,
              T=value_type,
              Policy=boost::math::policies::policy<boost::math::policies::default_policy,boost::math::policies::default_policy,boost::math::policies::default_policy,boost::math::policies::default_policy,boost::math::policies::default_policy,boost::math::policies::default_policy,boost::math::policies::default_policy,boost::math::policies::default_policy,boost::math::policies::default_policy,boost::math::policies::default_policy,boost::math::policies::default_policy>
          ]
          I:\boost-trunk\libs\math\example\bessel_zeros_example.cpp(46) : see reference to function template instantiation 'void boost::math::cyl_bessel_j_zero<std::_Vector_iterator<_Myvec>,int>(output_iterator,T,unsigned int,unsigned int)' being compiled
          with
          [
              _Myvec=std::_Vector_val<double,std::allocator<double>>,
              output_iterator=std::_Vector_iterator<std::_Vector_val<double,std::allocator<double>>>,
              T=int
          ]
Manifest:
  Deleting file "Debug\bessel_zeros_example.exe.embed.manifest".
LinkEmbedManifest:
  bessel_zeros_example.vcxproj -> J:\Cpp\big_number\Debug\bessel_zeros_example.exe
CustomBuildStep:
  Description: Autorun "J:\Cpp\big_number\Debug\bessel_zeros_example.exe"
  1st root = 5.13562
  2nd root = 8.41724
  3rd root = 11.6198
  7.2731751938316489503185694262290765588963196701623
  10.724858308883141732536172745851416647110749599085
  14.018504599452388106120459558042660282427471931581
  17.25249845917041718216248716654977734919590383861
  20.456678874044517595180234083894285885460502077814
  23.64363089714234522494551422714731959985405172504
  26.819671140255087745421311470965019261522390519297
  29.988343117423674742679141796661432043878868194142
  33.151796897690520871250862469973445265444791966114
  36.3114160002162074157243540350393860813165201842
  39.468132467505236587945197808083337887765967032029
  42.622597801391236474855034831297954018844433480227
  45.775281464536847753390206207806726581495950012439
  48.926530489173566198367766817478553992471739894799
  52.076607045343002794279746041878924876873478063472
  55.225712944912571393594224327817265689059002890192
  58.374006101538886436775188150439025201735151418932
  61.521611873000965273726742659353136266390944103571
  64.66863105379093036834648221487366079456596628716
  67.815145619696290925556791375555951165111460585458
FinalizeBuildStatus:
  Deleting file "Debug\bessel_zeros_example.unsuccessfulbuild".
  Touching "Debug\bessel_zeros_example.lastbuildstate".

Build succeeded.

Time Elapsed 00:00:22.68
*/
