// test_error_handling.cpp

// Test error handling.

// Copyright Paul A. Bristow 2006.
// Copyright John Maddock 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_MATH_THROW_ON_DOMAIN_ERROR

// Boost
#include <boost/math/tools/error_handling.hpp> // for domain_error.
	using ::boost::math::tools::domain_error;
	using ::boost::math::tools::pole_error;
	using ::boost::math::tools::overflow_error;
	using ::boost::math::tools::underflow_error;
	using ::boost::math::tools::denorm_error;
	using ::boost::math::tools::logic_error;
#include <boost/math/special_functions/chisqr.hpp> // chisqr
   using ::boost::math::chisqr;
#include <boost/math/special_functions/gamma.hpp> // tgamma
   using ::boost::math::tgamma;
#include <boost/math/special_functions/fpclassify.hpp>
	using boost::math::isnan;
// std
#include <iostream>
	using std::cout;
	using std::endl;
#include <limits>
  using std::numeric_limits;
#include <stdexcept>
	using std::exception;

#include <boost/test/included/test_exec_monitor.hpp> // for test_main
#include <boost/test/floating_point_comparison.hpp> // for BOOST_CHECK_CLOSE

#include <boost/math/concepts/real_concept.hpp> // for real_concept
   using ::boost::math::concepts::real_concept;

//
// The macro BOOST_CHECK_THROW_MSG is the same as BOOST_CHECK_THROW
// which is to say it verifies that the exception we expect is
// thrown from the function under test, but it also prints the message
// contained in the thrown exception so we can manually inspect the
// quality of the message.
//
// The expanded code is:
//
//   try
//   { 
//      code; 
//   } 
//   catch(const std::exception& e)
//   {
//      std::cout << 
//          "Message from thrown exception was:\n   " << e.what() << std::endl; 
//   }
//   BOOST_CHECK_THROW(code, t);
//
#define BOOST_CHECK_THROW_MSG(code,t)\
   try{ code; } catch(const std::exception& e){\
   std::cout << "Message from thrown exception was:\n   " << e.what() << std::endl; }\
   BOOST_CHECK_THROW(code, t);

template <class FPT> // Any floating-point type FPT.
void test_error(FPT)
{
	cout << "Current function is " << BOOST_CURRENT_FUNCTION << endl;
// 2 argument version now removed, so these two commented out.
//   BOOST_CHECK_THROW_MSG(domain_error<FPT>(BOOST_CURRENT_FUNCTION, 0), std::domain_error);
//   BOOST_CHECK_THROW_MSG(domain_error<FPT>(BOOST_CURRENT_FUNCTION, "Out of range argument"), std::domain_error);
   BOOST_CHECK_THROW_MSG(domain_error<FPT>(BOOST_CURRENT_FUNCTION, 0, static_cast<FPT>(3.124567890123456789012345678901L)), std::domain_error);
   BOOST_CHECK_THROW_MSG(domain_error<FPT>(BOOST_CURRENT_FUNCTION, "Out of range argument %1% in test invocation", static_cast<FPT>(3.124567890123456789012345678901L)), std::domain_error);
//   BOOST_CHECK_THROW_MSG(logic_error<FPT>(BOOST_CURRENT_FUNCTION, 0), std::logic_error);
//   BOOST_CHECK_THROW_MSG(logic_error<FPT>(BOOST_CURRENT_FUNCTION, "Internal error"), std::logic_error);
   BOOST_CHECK_THROW_MSG(logic_error<FPT>(BOOST_CURRENT_FUNCTION, 0, static_cast<FPT>(3.124567890123456789012345678901L)), std::logic_error);
   BOOST_CHECK_THROW_MSG(logic_error<FPT>(BOOST_CURRENT_FUNCTION, "Internal error, computed result was %1%, but should be in the range [0,1]", static_cast<FPT>(3.124567890123456789012345678901L)), std::logic_error);
   BOOST_CHECK_THROW_MSG(chisqr(-1, static_cast<FPT>(1)), std::domain_error);
   BOOST_CHECK_THROW_MSG(chisqr(static_cast<FPT>(1), static_cast<FPT>(-1)), std::domain_error);
   BOOST_CHECK_THROW_MSG(chisqr(static_cast<FPT>(0), static_cast<FPT>(1)), std::domain_error);
   BOOST_CHECK_THROW_MSG(tgamma(static_cast<FPT>(0)), std::domain_error);
   BOOST_CHECK_THROW_MSG(tgamma(static_cast<FPT>(-10)), std::domain_error);
   BOOST_CHECK_THROW_MSG(tgamma(static_cast<FPT>(-10123457772243.0)), std::domain_error);

   cout << endl;

} // template <class FPT>void test_error(FPT)

int test_main(int, char* [])
{
	// Test error handling.
	// (Parameter value, arbitrarily zero, only communicates the floating point type FPT).
	test_error(0.0F); // Test float.
	test_error(0.0); // Test double.
	test_error(0.0L); // Test long double.
	test_error(real_concept(0.0L)); // Test concepts.
	return 0;
} // int test_main(int, char* [])

/*

Output:

------ Build started: Project: test_error_handling, Configuration: Debug Win32 ------
Compiling...
test_error_handling.cpp
Linking...
Embedding manifest...
Autorun "i:\boost-06-05-03-1300\libs\math\test\Math_test\debug\test_error_handling.exe"
Running 1 test case...
Current function is void __cdecl test_error<float>(float)
Message from thrown exception was:
   Error in function void __cdecl test_error<float>(float): Domain Error on value 3.12456799
Message from thrown exception was:
   Error in function void __cdecl test_error<float>(float): Out of range argument 3.12456799 in test invocation
Message from thrown exception was:
   Error in function void __cdecl test_error<float>(float): Internal logic error, computed value was 3.12456799
Message from thrown exception was:
   Error in function void __cdecl test_error<float>(float): Internal error, computed result was 3.12456799, but should be in the range [0,1]
Message from thrown exception was:
   Error in function double __cdecl boost::math::detail::chisqr_imp<double>(double,double): degrees of freedom argument is -1, but must be > 0 !
Message from thrown exception was:
   Error in function float __cdecl boost::math::detail::chisqr_imp<float>(float,float): chisqr argument is -1, but must be > 0 !
Message from thrown exception was:
   Error in function float __cdecl boost::math::detail::chisqr_imp<float>(float,float): degrees of freedom argument is 0, but must be > 0 !
Message from thrown exception was:
   Error in function double __cdecl boost::math::detail::gamma_imp<double,struct boost::math::lanczos::lanczos6>(double,const struct boost::math::lanczos::lanczos6 &): Evaluation of tgamma at a negative integer 0.
Message from thrown exception was:
   Error in function double __cdecl boost::math::detail::gamma_imp<double,struct boost::math::lanczos::lanczos6>(double,const struct boost::math::lanczos::lanczos6 &): Evaluation of tgamma at a negative integer -10.
Message from thrown exception was:
   Error in function double __cdecl boost::math::detail::gamma_imp<double,struct boost::math::lanczos::lanczos6>(double,const struct boost::math::lanczos::lanczos6 &): Evaluation of tgamma at a negative integer -10123458117632.
Current function is void __cdecl test_error<double>(double)
Message from thrown exception was:
   Error in function void __cdecl test_error<double>(double): Domain Error on value 3.1245678901234566
Message from thrown exception was:
   Error in function void __cdecl test_error<double>(double): Out of range argument 3.1245678901234566 in test invocation
Message from thrown exception was:
   Error in function void __cdecl test_error<double>(double): Internal logic error, computed value was 3.1245678901234566
Message from thrown exception was:
   Error in function void __cdecl test_error<double>(double): Internal error, computed result was 3.1245678901234566, but should be in the range [0,1]
Message from thrown exception was:
   Error in function double __cdecl boost::math::detail::chisqr_imp<double>(double,double): degrees of freedom argument is -1, but must be > 0 !
Message from thrown exception was:
   Error in function double __cdecl boost::math::detail::chisqr_imp<double>(double,double): chisqr argument is -1, but must be > 0 !
Message from thrown exception was:
   Error in function double __cdecl boost::math::detail::chisqr_imp<double>(double,double): degrees of freedom argument is 0, but must be > 0 !
Message from thrown exception was:
   Error in function double __cdecl boost::math::detail::gamma_imp<double,struct boost::math::lanczos::lanczos13m53>(double,const struct boost::math::lanczos::lanczos13m53 &): Evaluation of tgamma at a negative integer 0.
Message from thrown exception was:
   Error in function double __cdecl boost::math::detail::gamma_imp<double,struct boost::math::lanczos::lanczos13m53>(double,const struct boost::math::lanczos::lanczos13m53 &): Evaluation of tgamma at a negative integer -10.
Message from thrown exception was:
   Error in function double __cdecl boost::math::detail::gamma_imp<double,struct boost::math::lanczos::lanczos13m53>(double,const struct boost::math::lanczos::lanczos13m53 &): Evaluation of tgamma at a negative integer -10123457772243.
Current function is void __cdecl test_error<long double>(long double)
Message from thrown exception was:
   Error in function void __cdecl test_error<long double>(long double): Domain Error on value 3.1245678901234566
Message from thrown exception was:
   Error in function void __cdecl test_error<long double>(long double): Out of range argument 3.1245678901234566 in test invocation
Message from thrown exception was:
   Error in function void __cdecl test_error<long double>(long double): Internal logic error, computed value was 3.1245678901234566
Message from thrown exception was:
   Error in function void __cdecl test_error<long double>(long double): Internal error, computed result was 3.1245678901234566, but should be in the range [0,1]
Message from thrown exception was:
   Error in function long double __cdecl boost::math::detail::chisqr_imp<long double>(long double,long double): degrees of freedom argument is -1, but must be > 0 !
Message from thrown exception was:
   Error in function long double __cdecl boost::math::detail::chisqr_imp<long double>(long double,long double): chisqr argument is -1, but must be > 0 !
Message from thrown exception was:
   Error in function long double __cdecl boost::math::detail::chisqr_imp<long double>(long double,long double): degrees of freedom argument is 0, but must be > 0 !
Message from thrown exception was:
   Error in function long double __cdecl boost::math::detail::gamma_imp<long double,struct boost::math::lanczos::lanczos13m53>(long double,const struct boost::math::lanczos::lanczos13m53 &): Evaluation of tgamma at a negative integer 0.
Message from thrown exception was:
   Error in function long double __cdecl boost::math::detail::gamma_imp<long double,struct boost::math::lanczos::lanczos13m53>(long double,const struct boost::math::lanczos::lanczos13m53 &): Evaluation of tgamma at a negative integer -10.
Message from thrown exception was:
   Error in function long double __cdecl boost::math::detail::gamma_imp<long double,struct boost::math::lanczos::lanczos13m53>(long double,const struct boost::math::lanczos::lanczos13m53 &): Evaluation of tgamma at a negative integer -10123457772243.
Current function is void __cdecl test_error<class boost::math::concepts::real_concept>(class boost::math::concepts::real_concept)
Message from thrown exception was:
   Error in function void __cdecl test_error<class boost::math::concepts::real_concept>(class boost::math::concepts::real_concept): Domain Error on value 3.1245678901234566
Message from thrown exception was:
   Error in function void __cdecl test_error<class boost::math::concepts::real_concept>(class boost::math::concepts::real_concept): Out of range argument 3.1245678901234566 in test invocation
Message from thrown exception was:
   Error in function void __cdecl test_error<class boost::math::concepts::real_concept>(class boost::math::concepts::real_concept): Internal logic error, computed value was 3.1245678901234566
Message from thrown exception was:
   Error in function void __cdecl test_error<class boost::math::concepts::real_concept>(class boost::math::concepts::real_concept): Internal error, computed result was 3.1245678901234566, but should be in the range [0,1]
Message from thrown exception was:
   Error in function class boost::math::concepts::real_concept __cdecl boost::math::detail::chisqr_imp<class boost::math::concepts::real_concept>(class boost::math::concepts::real_concept,class boost::math::concepts::real_concept): degrees of freedom argument is -1, but must be > 0 !
Message from thrown exception was:
   Error in function class boost::math::concepts::real_concept __cdecl boost::math::detail::chisqr_imp<class boost::math::concepts::real_concept>(class boost::math::concepts::real_concept,class boost::math::concepts::real_concept): chisqr argument is -1, but must be > 0 !
Message from thrown exception was:
   Error in function class boost::math::concepts::real_concept __cdecl boost::math::detail::chisqr_imp<class boost::math::concepts::real_concept>(class boost::math::concepts::real_concept,class boost::math::concepts::real_concept): degrees of freedom argument is 0, but must be > 0 !
Message from thrown exception was:
   Error in function class boost::math::concepts::real_concept __cdecl boost::math::detail::gamma_imp<class boost::math::concepts::real_concept>(class boost::math::concepts::real_concept,const struct boost::math::lanczos::undefined_lanczos &): Evaluation of tgamma at a negative integer 0.
Message from thrown exception was:
   Error in function class boost::math::concepts::real_concept __cdecl boost::math::detail::gamma_imp<class boost::math::concepts::real_concept>(class boost::math::concepts::real_concept,const struct boost::math::lanczos::undefined_lanczos &): Evaluation of tgamma at a negative integer -10.
Message from thrown exception was:
   Error in function class boost::math::concepts::real_concept __cdecl boost::math::detail::gamma_imp<class boost::math::concepts::real_concept>(class boost::math::concepts::real_concept,const struct boost::math::lanczos::undefined_lanczos &): Evaluation of tgamma at a negative integer -10123457772243.
*** No errors detected
Build Time 0:07
Build log was saved at "file://i:\boost-06-05-03-1300\libs\math\test\Math_test\test_error_handling\Debug\BuildLog.htm"
test_error_handling - 0 error(s), 0 warning(s)
========== Build: 1 succeeded, 0 failed, 0 up-to-date, 0 skipped ==========

*/
