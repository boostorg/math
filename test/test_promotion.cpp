// boost\math\test\test_promotion.hpp

// Copyright John Maddock & Paul A. Bristow 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_MATH_THROW_ON_DOMAIN_ERROR

#ifdef _MSC_VER
#  pragma warning(disable: 4127) // conditional expression is constant.
//#  pragma warning(disable: 4100) // unreferenced formal parameter.
#  pragma warning(disable: 4512) // assignment operator could not be generated.
#  if !(defined _SCL_SECURE_NO_DEPRECATE) || (_SCL_SECURE_NO_DEPRECATE == 0)
#     pragma warning(disable: 4996) // 'std::char_traits<char>::copy' was declared deprecated.
// #define _SCL_SECURE_NO_DEPRECATE = 1 // avoid C4996 warning.
#  endif
#  pragma warning(disable: 4244) // conversion from 'double' to 'float', possible loss of data.
#endif

#include <boost/test/included/test_exec_monitor.hpp> // Boost.Test
#include <boost/test/floating_point_comparison.hpp>
#include <boost/current_function.hpp> // for BOOST_CURRENT_FUNCTION
#include <boost/type_traits/is_same.hpp>
  using boost::is_same;

#include <boost/math/concepts/real_concept.hpp> // for real_concept

#include <boost/math/tools/promotion.hpp>
  using boost::math::tools::promote_arg2;
  using boost::math::tools::promote_arg3;

#include <iostream>
  using std::cout;
  using std::endl;
  using std::setprecision;

namespace boost
{
  namespace math
  {
    namespace detail
    {
      // RealType is floating-point (built-in or User Defined), IntegerType is integer, NumericType is either.
      template <class RealType>
      RealType test2_imp(RealType arg1, RealType arg2)
      { // Arguments promoted to same RealType, if necessary.
        cout << "Current function test is " << BOOST_CURRENT_FUNCTION << endl;
        RealType result = arg1 * arg2; // A useless test function.
        cout << "Result for test_func2 is " << result << endl;
        return result;
      } // test_imp

      // RealType is floating-point (built-in or User Defined), IntegerType is integer, NumericType is either.
      template <class RealType>
      RealType test3_imp(RealType arg1, RealType arg2, RealType arg3)
      { // Arguments have been promoted to same RealType, if necessary.
        cout << "Current function test is " << BOOST_CURRENT_FUNCTION << endl;
        RealType result = arg1 * arg2 * arg3; // A useless test function.
        cout << "Result for test_func3 is " << result << endl;
        return result;
      } // test_imp
    } // namespace detail
  } // namespace math
} // namespace boost

template <class NumericType, class RealType>
inline typename promote_arg2<RealType, NumericType>::type 
  test_func2(RealType arg1, NumericType arg2)
{ // return type is the wider of the two (if necessary, promoted) floating-point types.
  typedef typename promote_arg2<RealType, NumericType>::type promote_type; // Arguments type.
  // A RealType from rules above.

  cout << "arg1 is " << typeid(arg1).name() 
    << ", arg2 is " << typeid(arg2).name() 
    << ", promoted type is " << typeid(promote_type).name() << endl;

  return boost::math::detail::test2_imp(static_cast<promote_type>(arg1), static_cast<promote_type>(arg2));
} // test_func2

template <class NumericType, class RealType>
inline typename promote_arg3<RealType, NumericType, NumericType>::type
  test_func3(RealType arg1, NumericType arg2, NumericType arg3)
{  // return type is the wider of the two (if necessary, promoted) floating-point types.
  typedef typename promote_arg3<RealType, NumericType, NumericType>::type promote_type; // Arguments type.
  // A RealType from rules above.
  cout << "arg1 is " << typeid(arg1).name() 
    << ", arg2 is " << typeid(arg2).name() 
    << ", arg3 is " << typeid(arg3).name() 
    << ", promoted type is " << typeid(promote_type).name() << endl;
  return boost::math::detail::test3_imp(static_cast<promote_type>(arg1), static_cast<promote_type>(arg2), static_cast<promote_type>(arg3));
} // test_func3

template <class RealType>
void test_spots(RealType)
{
  cout << "Current floating-point type function is " << BOOST_CURRENT_FUNCTION << endl;

 // Tests of two argument function.
  RealType x = test_func2(1.F, 1.F); // Tests of two argument function.
        x = test_func2(1.F, 2.); // warning C4244: '=' : conversion from 'double' to 'float', possible loss of data
        x = test_func2(1.F, 3.L); 
        x = test_func2(1.F, 4L);
        x = test_func2(2., 5);
        x = test_func2(2., 6L);
        x = test_func2(2., 7.);  
        x = test_func2(2., 8.L);
        x = test_func2(3.L, 9.L);
        x = test_func2(4, 10.L);
        x = test_func2(4.F, 11.L);
        x = test_func2(4.L, 12.L);

 // Tests of three argument function.
  RealType y = test_func3(1.F, 1.F, 1.F);
        y = test_func3(1.F, 2, 3); // warning C4244: '=' : conversion from 'double' to 'float', possible loss of data
        y = test_func3(1.F, 4L, 4L);
        y = test_func3(1.F, 4L, 4L);
        y = test_func3(1., 4L, 4L);
} // void test_spots(RealType)

int test_main(int, char* [])
{
  // Basic sanity-check spot values.
  // (Parameter value, arbitrarily zero, only communicates the RealType (floating_point)).
  test_spots(0.0F); // Test float. 
  // Many float from double conversions here as RealType is float.
  test_spots(0.0); // Test double.
  // Some double from long double conversions here.
  test_spots(0.0L); // Test long double.
  test_spots(boost::math::concepts::real_concept(0.)); // Test real concept.

  return 0;
} // int test_main(int, char* [])


/*

Output:

------ Build started: Project: test_promotion, Configuration: Debug Win32 ------
Compiling...
test_promotion.cpp
Linking...
Autorun "i:\boost-06-05-03-1300\libs\math\test\Math_test\debug\test_promotion.exe"
Running 1 test case...
Current floating-point type function is void __cdecl test_spots<float>(float)
arg1 is float, arg2 is float, promoted type is float
Current function test is float __cdecl boost::math::detail::test2_imp<float>(float,float)
Result for test_func2 is 1
arg1 is float, arg2 is double, promoted type is double
Current function test is double __cdecl boost::math::detail::test2_imp<double>(double,double)
Result for test_func2 is 2
arg1 is float, arg2 is long double, promoted type is long double
Current function test is long double __cdecl boost::math::detail::test2_imp<long double>(long double,long double)
Result for test_func2 is 3
arg1 is float, arg2 is long, promoted type is double
Current function test is double __cdecl boost::math::detail::test2_imp<double>(double,double)
Result for test_func2 is 4
arg1 is double, arg2 is int, promoted type is double
Current function test is double __cdecl boost::math::detail::test2_imp<double>(double,double)
Result for test_func2 is 10
arg1 is double, arg2 is long, promoted type is double
Current function test is double __cdecl boost::math::detail::test2_imp<double>(double,double)
Result for test_func2 is 12
arg1 is double, arg2 is double, promoted type is double
Current function test is double __cdecl boost::math::detail::test2_imp<double>(double,double)
Result for test_func2 is 14
arg1 is double, arg2 is long double, promoted type is long double
Current function test is long double __cdecl boost::math::detail::test2_imp<long double>(long double,long double)
Result for test_func2 is 16
arg1 is long double, arg2 is long double, promoted type is long double
Current function test is long double __cdecl boost::math::detail::test2_imp<long double>(long double,long double)
Result for test_func2 is 27
arg1 is int, arg2 is long double, promoted type is long double
Current function test is long double __cdecl boost::math::detail::test2_imp<long double>(long double,long double)
Result for test_func2 is 40
arg1 is float, arg2 is long double, promoted type is long double
Current function test is long double __cdecl boost::math::detail::test2_imp<long double>(long double,long double)
Result for test_func2 is 44
arg1 is long double, arg2 is long double, promoted type is long double
Current function test is long double __cdecl boost::math::detail::test2_imp<long double>(long double,long double)
Result for test_func2 is 48
arg1 is float, arg2 is float, arg3 is float, promoted type is float
Current function test is float __cdecl boost::math::detail::test3_imp<float>(float,float,float)
Result for test_func3 is 1
arg1 is float, arg2 is int, arg3 is int, promoted type is double
Current function test is double __cdecl boost::math::detail::test3_imp<double>(double,double,double)
Result for test_func3 is 6
arg1 is float, arg2 is long, arg3 is long, promoted type is double
Current function test is double __cdecl boost::math::detail::test3_imp<double>(double,double,double)
Result for test_func3 is 16
arg1 is float, arg2 is long, arg3 is long, promoted type is double
Current function test is double __cdecl boost::math::detail::test3_imp<double>(double,double,double)
Result for test_func3 is 16
arg1 is double, arg2 is long, arg3 is long, promoted type is double
Current function test is double __cdecl boost::math::detail::test3_imp<double>(double,double,double)
Result for test_func3 is 16
Current floating-point type function is void __cdecl test_spots<double>(double)
arg1 is float, arg2 is float, promoted type is float
Current function test is float __cdecl boost::math::detail::test2_imp<float>(float,float)
Result for test_func2 is 1
arg1 is float, arg2 is double, promoted type is double
Current function test is double __cdecl boost::math::detail::test2_imp<double>(double,double)
Result for test_func2 is 2
arg1 is float, arg2 is long double, promoted type is long double
Current function test is long double __cdecl boost::math::detail::test2_imp<long double>(long double,long double)
Result for test_func2 is 3
arg1 is float, arg2 is long, promoted type is double
Current function test is double __cdecl boost::math::detail::test2_imp<double>(double,double)
Result for test_func2 is 4
arg1 is double, arg2 is int, promoted type is double
Current function test is double __cdecl boost::math::detail::test2_imp<double>(double,double)
Result for test_func2 is 10
arg1 is double, arg2 is long, promoted type is double
Current function test is double __cdecl boost::math::detail::test2_imp<double>(double,double)
Result for test_func2 is 12
arg1 is double, arg2 is double, promoted type is double
Current function test is double __cdecl boost::math::detail::test2_imp<double>(double,double)
Result for test_func2 is 14
arg1 is double, arg2 is long double, promoted type is long double
Current function test is long double __cdecl boost::math::detail::test2_imp<long double>(long double,long double)
Result for test_func2 is 16
arg1 is long double, arg2 is long double, promoted type is long double
Current function test is long double __cdecl boost::math::detail::test2_imp<long double>(long double,long double)
Result for test_func2 is 27
arg1 is int, arg2 is long double, promoted type is long double
Current function test is long double __cdecl boost::math::detail::test2_imp<long double>(long double,long double)
Result for test_func2 is 40
arg1 is float, arg2 is long double, promoted type is long double
Current function test is long double __cdecl boost::math::detail::test2_imp<long double>(long double,long double)
Result for test_func2 is 44
arg1 is long double, arg2 is long double, promoted type is long double
Current function test is long double __cdecl boost::math::detail::test2_imp<long double>(long double,long double)
Result for test_func2 is 48
arg1 is float, arg2 is float, arg3 is float, promoted type is float
Current function test is float __cdecl boost::math::detail::test3_imp<float>(float,float,float)
Result for test_func3 is 1
arg1 is float, arg2 is int, arg3 is int, promoted type is double
Current function test is double __cdecl boost::math::detail::test3_imp<double>(double,double,double)
Result for test_func3 is 6
arg1 is float, arg2 is long, arg3 is long, promoted type is double
Current function test is double __cdecl boost::math::detail::test3_imp<double>(double,double,double)
Result for test_func3 is 16
arg1 is float, arg2 is long, arg3 is long, promoted type is double
Current function test is double __cdecl boost::math::detail::test3_imp<double>(double,double,double)
Result for test_func3 is 16
arg1 is double, arg2 is long, arg3 is long, promoted type is double
Current function test is double __cdecl boost::math::detail::test3_imp<double>(double,double,double)
Result for test_func3 is 16
Current floating-point type function is void __cdecl test_spots<long double>(long double)
arg1 is float, arg2 is float, promoted type is float
Current function test is float __cdecl boost::math::detail::test2_imp<float>(float,float)
Result for test_func2 is 1
arg1 is float, arg2 is double, promoted type is double
Current function test is double __cdecl boost::math::detail::test2_imp<double>(double,double)
Result for test_func2 is 2
arg1 is float, arg2 is long double, promoted type is long double
Current function test is long double __cdecl boost::math::detail::test2_imp<long double>(long double,long double)
Result for test_func2 is 3
arg1 is float, arg2 is long, promoted type is double
Current function test is double __cdecl boost::math::detail::test2_imp<double>(double,double)
Result for test_func2 is 4
arg1 is double, arg2 is int, promoted type is double
Current function test is double __cdecl boost::math::detail::test2_imp<double>(double,double)
Result for test_func2 is 10
arg1 is double, arg2 is long, promoted type is double
Current function test is double __cdecl boost::math::detail::test2_imp<double>(double,double)
Result for test_func2 is 12
arg1 is double, arg2 is double, promoted type is double
Current function test is double __cdecl boost::math::detail::test2_imp<double>(double,double)
Result for test_func2 is 14
arg1 is double, arg2 is long double, promoted type is long double
Current function test is long double __cdecl boost::math::detail::test2_imp<long double>(long double,long double)
Result for test_func2 is 16
arg1 is long double, arg2 is long double, promoted type is long double
Current function test is long double __cdecl boost::math::detail::test2_imp<long double>(long double,long double)
Result for test_func2 is 27
arg1 is int, arg2 is long double, promoted type is long double
Current function test is long double __cdecl boost::math::detail::test2_imp<long double>(long double,long double)
Result for test_func2 is 40
arg1 is float, arg2 is long double, promoted type is long double
Current function test is long double __cdecl boost::math::detail::test2_imp<long double>(long double,long double)
Result for test_func2 is 44
arg1 is long double, arg2 is long double, promoted type is long double
Current function test is long double __cdecl boost::math::detail::test2_imp<long double>(long double,long double)
Result for test_func2 is 48
arg1 is float, arg2 is float, arg3 is float, promoted type is float
Current function test is float __cdecl boost::math::detail::test3_imp<float>(float,float,float)
Result for test_func3 is 1
arg1 is float, arg2 is int, arg3 is int, promoted type is double
Current function test is double __cdecl boost::math::detail::test3_imp<double>(double,double,double)
Result for test_func3 is 6
arg1 is float, arg2 is long, arg3 is long, promoted type is double
Current function test is double __cdecl boost::math::detail::test3_imp<double>(double,double,double)
Result for test_func3 is 16
arg1 is float, arg2 is long, arg3 is long, promoted type is double
Current function test is double __cdecl boost::math::detail::test3_imp<double>(double,double,double)
Result for test_func3 is 16
arg1 is double, arg2 is long, arg3 is long, promoted type is double
Current function test is double __cdecl boost::math::detail::test3_imp<double>(double,double,double)
Result for test_func3 is 16
Current floating-point type function is void __cdecl test_spots<class boost::math::concepts::real_concept>(class boost::math::concepts::real_concept)
arg1 is float, arg2 is float, promoted type is float
Current function test is float __cdecl boost::math::detail::test2_imp<float>(float,float)
Result for test_func2 is 1
arg1 is float, arg2 is double, promoted type is double
Current function test is double __cdecl boost::math::detail::test2_imp<double>(double,double)
Result for test_func2 is 2
arg1 is float, arg2 is long double, promoted type is long double
Current function test is long double __cdecl boost::math::detail::test2_imp<long double>(long double,long double)
Result for test_func2 is 3
arg1 is float, arg2 is long, promoted type is double
Current function test is double __cdecl boost::math::detail::test2_imp<double>(double,double)
Result for test_func2 is 4
arg1 is double, arg2 is int, promoted type is double
Current function test is double __cdecl boost::math::detail::test2_imp<double>(double,double)
Result for test_func2 is 10
arg1 is double, arg2 is long, promoted type is double
Current function test is double __cdecl boost::math::detail::test2_imp<double>(double,double)
Result for test_func2 is 12
arg1 is double, arg2 is double, promoted type is double
Current function test is double __cdecl boost::math::detail::test2_imp<double>(double,double)
Result for test_func2 is 14
arg1 is double, arg2 is long double, promoted type is long double
Current function test is long double __cdecl boost::math::detail::test2_imp<long double>(long double,long double)
Result for test_func2 is 16
arg1 is long double, arg2 is long double, promoted type is long double
Current function test is long double __cdecl boost::math::detail::test2_imp<long double>(long double,long double)
Result for test_func2 is 27
arg1 is int, arg2 is long double, promoted type is long double
Current function test is long double __cdecl boost::math::detail::test2_imp<long double>(long double,long double)
Result for test_func2 is 40
arg1 is float, arg2 is long double, promoted type is long double
Current function test is long double __cdecl boost::math::detail::test2_imp<long double>(long double,long double)
Result for test_func2 is 44
arg1 is long double, arg2 is long double, promoted type is long double
Current function test is long double __cdecl boost::math::detail::test2_imp<long double>(long double,long double)
Result for test_func2 is 48
arg1 is float, arg2 is float, arg3 is float, promoted type is float
Current function test is float __cdecl boost::math::detail::test3_imp<float>(float,float,float)
Result for test_func3 is 1
arg1 is float, arg2 is int, arg3 is int, promoted type is double
Current function test is double __cdecl boost::math::detail::test3_imp<double>(double,double,double)
Result for test_func3 is 6
arg1 is float, arg2 is long, arg3 is long, promoted type is double
Current function test is double __cdecl boost::math::detail::test3_imp<double>(double,double,double)
Result for test_func3 is 16
arg1 is float, arg2 is long, arg3 is long, promoted type is double
Current function test is double __cdecl boost::math::detail::test3_imp<double>(double,double,double)
Result for test_func3 is 16
arg1 is double, arg2 is long, arg3 is long, promoted type is double
Current function test is double __cdecl boost::math::detail::test3_imp<double>(double,double,double)
Result for test_func3 is 16
*** No errors detected
Build Time 0:04
Build log was saved at "file://i:\boost-06-05-03-1300\libs\math\test\Math_test\promotion\Debug\BuildLog.htm"
test_promotion - 0 error(s), 0 warning(s)
========== Build: 1 succeeded, 0 failed, 0 up-to-date, 0 skipped ==========




*/

