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

#include <boost/multiprecision/cpp_dec_float.hpp>

#include <boost/math/special_functions/math_fwd.hpp>
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

// Weisstein, Eric W. "Bessel Function Zeros." From MathWorld--A Wolfram Web Resource.
// http://mathworld.wolfram.com/BesselFunctionZeros.html 
// See also http://dlmf.nist.gov/10.21

// To use the Boost.Multiprecision 50 decimal digits type.
typedef boost::multiprecision::cpp_dec_float_50 float_type;

using namespace boost::math;

/*`This example shows obtaining both a single zero or root of the Bessel function,
and placing multiple zeros into a container like `std::vector` by providing an iterator.
The signature of the single value function is:


  template <class T>
  inline typename detail::bessel_traits<T, T, policies::policy<> >::result_type 
    cyl_bessel_j_zero(T v,  // Floating-point value for Jv.
    unsigned m); // start index.

The result type is controlled by the floating-point type of parameter `v`
(but subject to the usual policy and promotion rules).

and for multiple zeros:

  template <class T, class OutputIterator>
  inline OutputIterator cyl_bessel_j_zero(T v, // Floating-point value for Jv.
                                unsigned number_of_zeros,
                                unsigned start_index,
                                OutputIterator out_it); // 

There is also a version which allows control of the policy for error handling and precision.

  template <class T, class OutputIterator, class Policy>
  inline OutputIterator cyl_bessel_j_zero(T v, // Floating-point value for Jv.
                                unsigned number_of_zeros,
                                unsigned start_index,
                                OutputIterator out_it,
                                const Policy& pol);

This exmaple also shows how to use the output iterator to create a sum of 1/zeros^2.
*/

template <class T>
struct output_summation_iterator
{
   output_summation_iterator(T* p) : p_sum(p)
   {
   }
   output_summation_iterator& operator*()
   {
     return *this;
   }
   output_summation_iterator& operator = (T const& val)
   {
     *p_sum += 1./ (val * val); // 1/zero^2
     //std::cout << *p_sum << ", ";
     return *this;
   }
   output_summation_iterator& operator++()
   { return *this; }
   output_summation_iterator& operator++(int)
   {
     return *this;
   }
private:
   T* p_sum;
};

double s_m_nu(int m, int nu)
{
  switch(m)
  {

    case 2:
    return 1./(4 * (nu + 1));
    case 4:
      return 1./(16 * (nu+1)*(nu+1)*(nu+2));
    default:
      return 0;
  }
} // double s_m_nu(int m, int nu)

int main()
{

  { // Evaluate a single Bessel zero.

   //  T cyl_bessel_j_zero(T v, unsigned int index)
    // The precision is controlled by the float-point type of template parameter T of v.
    // so this example has double precision, at least 15 but up to 17 decimal digits.
    double root = boost::math::cyl_bessel_j_zero(0.0, 1U);
    // Using default precision of 6 decimal digits:
    std::cout << "boost::math::cyl_bessel_j_zero(0.0, 1U) " << root << std::endl; // 2.40483
    std::cout.precision(std::numeric_limits<double>::digits10);
    std::cout << "boost::math::cyl_bessel_j_zero(0.0, 1U) " << root << std::endl; // 2.40482555769577


/*`But note that because the parameter v controls the precision of the result,
it *must* be a [[floating-point type].
So if you provide an integer type, say 0, rather than 0.0, then it will fail to compile thus:
``
    root = boost::math::cyl_bessel_j_zero(0, 1U);
``
error C2338: Order must be a floating-point type.
*/
  }
  {

    // IAN N. SNEDDON, Infinite sums of Bessel Zeros.
    // page 150 equation 40.
    using boost::math::cyl_bessel_j_zero;
    std::cout.precision(std::numeric_limits<double>::digits10);
    double nu = 1.;
    double sum = 0;
    output_summation_iterator<double> it(&sum);  // sum of 1/zeros^2
    cyl_bessel_j_zero(nu, 100000U, 1U, it);

    std::cout << "Final " << sum << std::endl; //   0.0 Final 0.249999
    // 1.0  Final 0.124998986795763

    double s = 1/(4 * (nu + 1)); // 0.125 = 1/8 is exact analytical solution.
    std::cout << s << std::endl;

  }

  {
/*`The Neumann functions zeros are evaluated very similarly:
*/

    using boost::math::cyl_neumann_zero;

    double sum = 0;
    output_summation_iterator<double> it(&sum);
    cyl_neumann_zero(2.5, 1, 10, it);

    std::cout << sum << std::endl;

  }
/*`Another version allows calculation of multiple zeros with one call,
placing the results in a container, often `std::vector`.
For example, generate five double roots of Jv for integral order 2.
*/
  {
    double azero = boost::math::cyl_bessel_j_zero(0.0, 1U);



    unsigned int n_roots = 5U;
    std::vector<double> roots;

    boost::math::cyl_bessel_j_zero(0.0, n_roots, 1U, std::back_inserter(roots));

     // Note must provide an floating-point type, not an integer type, so v = 2.0, not 2.
     //boost::math::cyl_bessel_j_zero(std::back_inserter(roots_3), 2, 3U, 1U);
     // error C2338: Order must be a floating-point type.

      std::copy(roots.begin(),
              roots.end(),
              std::ostream_iterator<double>(std::cout, "\n"));

  // 5 roots v = 0.0
  // 1.#QNAN
  //2.40483
  //5.52008
  //8.65373
  //11.7915
  //14.9309

  }

/*`Or generate 20 roots of Jv for non-integral order v=71/19.
*/
  {
    // Set the precision of the output stream.
    std::cout << std::showpoint << std::endl; // Show trailing zeros.
    std::cout.precision(std::numeric_limits<float_type>::digits10);

    float_type x = float_type(71) / 19;
    float_type r = boost::math::cyl_bessel_j_zero(x, 1U);

    std::cout << "x = " << x << ", r = " << r << std::endl;

    r = boost::math::cyl_bessel_j_zero(x, 50U);

    std::cout << "x = " << x << ", r = " << r << std::endl;

    std::vector<float_type> roots(20U);

    boost::math::cyl_bessel_j_zero(float_type(71) / 19, unsigned(roots.size()), 1U, roots.begin());

    // Print the roots to the output stream.
    std::copy(roots.begin(),
              roots.end(),
              std::ostream_iterator<float_type>(std::cout, "\n"));
 }
/*
*/

/*`Test some corner cases:



*/
  try
  { 
/*
     [N[BesselJZero[0, 1000], 50]] 
     3140.8072952250786288955454534711266789940767025137
     3.140807295225079e+003
     j 1000(x = 0), r = 3.140807295225079e+003
     j 1000(x = 0.000000000000000), r = 3140.80729522508

    [N[BesselJZero[0, 1000000], 50]]
     3.1415918681916696297600539252789979000145664979511×10^6
     3.141591868191670e+006
      j 1000000(x = 0.000000000000000), r = 3141591.86819167

    [N[BesselJZero[0, 1000000000], 50]]
     3.1415926528043950751049838094467630562626412341405×10^9
     3.141592652804395e+009
       j 1000000000(x = 0), r = 3.141592652804395e+009
       j 1000000000(x = 0), r = 3141592652.8044

       */

//  [N[BesselJZero[0, 4294967295], 50]]
    
    std::cout.precision(std::numeric_limits<double>::digits10);
    double x = 0.;
    // double r = boost::math::cyl_bessel_j_zero(x, 1U); // 2.4
    // double r = boost::math::cyl_bessel_j_zero(x, 0U); // NAN

    unsigned int j = std::numeric_limits<unsigned int>::max();
    j = 1000U;
    double r = boost::math::cyl_bessel_j_zero(x, j);
   
    std::cout << "j " << j << "(x = " << x << "), r = "
      << std::scientific << r << std::endl;
    // j 4294967295(x = 0.000000000000000), r = 6746518848.33402

    // 1.3493037700595028141621005137780845320949701145378×10^10

  }
  catch (std::exception ex)
  {
    std::cout << "Throw exception " << ex.what() << std::endl;
  }

 } // int main()

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
