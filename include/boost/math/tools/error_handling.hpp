//  Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_ERROR_HANDLING_HPP
#define BOOST_MATH_TOOLS_ERROR_HANDLING_HPP

// Boost
#include <boost/throw_exception.hpp>
#include <boost/limits.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/math/tools/precision.hpp>
#ifdef BOOST_MSVC
#  pragma warning(push) // Quiet warnings in boost/format.hpp
#  pragma warning(disable: 4996) // _SCL_SECURE_NO_DEPRECATE
#  pragma warning(disable: 4512) // assignment operator could not be generated.
#endif
#include <boost/format.hpp>
#ifdef BOOST_MSVC
#  pragma warning(pop)
#endif
#ifdef BOOST_MSVC
#  pragma warning(push) // Quiet warnings in error_handling
// May also need to quiet warnings within distributions like binomial, etc
// depending on definition, or not, of macro BOOST_MATH_THROW_ON_DOMAIN_ERROR
// that controls if error_handling throws, or not, 
// that, in turn, determines if code is reachable, or not.
#  pragma warning(disable: 4702) // unreachable code in domain_error_imp etc
#endif

// We don't use this one directly,
// but our clients will use it for BOOST_CURRENT_FUNCTION.
#include <boost/current_function.hpp>

// std
#include <limits> // for numeric_limits
#include <stdexcept>
#include <cmath>
#include <cerrno> // for errno

namespace boost{ namespace math
{
  namespace tools
  {
    namespace detail
    {
	    template <class E>
	    void raise_error(const char* function, const char* message)
	    {
		     if(function == 0)
				    function = "Unknown function";
		     if(message == 0)
				    message = "Cause unknown";

		     std::string msg("Error in function ");
		     msg += function;
		     msg += ": ";
		     msg += message;

		     E e(msg);
		     boost::throw_exception(e);
	    }

	    template <class E, class T>
	    void raise_error(const char* function, const char* message, const T& val)
	    {
		     if(function == 0)
				    function = "Unknown function";
		     if(message == 0)
				    message = "Cause unknown";

		     std::string msg("Error in function ");
		     msg += function;
		     msg += ": ";
		     msg += message;

		     int prec = 2 + (tools::digits<T>() * 30103UL) / 100000UL;
		     msg = (boost::format(msg) % boost::io::group(std::setprecision(prec), val)).str();

		     E e(msg);
		     boost::throw_exception(e);
	    }

      #ifndef BOOST_MATH_THROW_ON_DOMAIN_ERROR
        template <class T>
        inline T domain_error_imp(const char* /*function*/, const char* /*message*/, const T&, const mpl::true_*)
        { // No throw case.
           errno = ERANGE;
           return std::numeric_limits<T>::quiet_NaN();
        }
      #endif // BOOST_MATH_THROW_ON_DOMAIN_ERROR
      template <class T>
      inline T domain_error_imp(const char* function, const char* message, const T& val, const void*)
      {  // Throw case.
         errno = ERANGE;
         detail::raise_error<std::domain_error>(function, message ? message : "Domain Error on value %1%", val);
         // We don't get here: (and may get warning C4702: unreachable code)
         return 0;
      }

      #ifndef BOOST_MATH_THROW_ON_OVERFLOW_ERROR
        template <class T>
        inline T overflow_error_imp(const char* /*function*/, const char* /*message*/, mpl::bool_<true> const*)
        { // No throw case.
           errno = ERANGE;
           return std::numeric_limits<T>::infinity();
        }
      #endif // BOOST_MATH_THROW_ON_OVERFLOW_ERROR
      template <class T>
      inline T overflow_error_imp(const char* function, const char* message, const void*)
      {  // Throw case.
         errno = ERANGE;
         detail::raise_error<std::overflow_error>(function, message ? message : "Overflow");
         // We don't get here: (and may get warning C4702: unreachable code)
         return 0;
      }
    } // namespace detail

// boost::math::tools:: Error handling functions:
// domain_error, pole_error, overflow_error, underflow_error, denorm_error, logic_error
// & checked_narrowing_cast.

// An argument is outside its allowed range:
// so show message & return a NaN.
template <class T>
inline T domain_error(const char* function, const char* message, const T& val)
{ // returns std::numeric_limits<T>::quiet_NaN(), or as near as possible.
#ifdef BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS
   if(std::numeric_limits<T>::has_quiet_NaN)
      return detail::domain_error_imp<T>(function, message, val, static_cast<mpl::true_*>(0));
   return detail::domain_error_imp<T>(function, message, val, static_cast<mpl::false_*>(0));
#else
   typedef boost::mpl::bool_< ::std::numeric_limits<T>::has_quiet_NaN> tag_type;
   return detail::domain_error_imp<T>(function, message, val, static_cast<tag_type*>(0));
#endif
} // domain_error

// Evaluation at a pole, this is currently treated the same as a domain error:
template <class T>
inline T pole_error(const char* function, const char* message, const T& val)
{
   return domain_error<T>(function, message ? message : "Evaluation at pole %1%", val);
}
//
// Result too large to be represented in type T:
// so  show message & return an infinity.
template <class T>
inline T overflow_error(const char* function, const char* message)
{ // returns std::numeric_limits<T>::infinity(), or as near as possible.
#ifdef BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS
   if(std::numeric_limits<T>::has_infinity)
      return detail::overflow_error_imp<T>(function, message, static_cast<mpl::true_*>(0));
   return detail::overflow_error_imp<T>(function, message, static_cast<mpl::false_*>(0));
#else
   typedef boost::mpl::bool_< ::std::numeric_limits<T>::has_infinity> tag_type;
   return detail::overflow_error_imp<T>(function, message, static_cast<tag_type*>(0));
#endif
} // overflow_error

// Result too small to be represented in type T,
// called only when we know the result is not actually zero:
// but show message & return true zero.
template <class T>
#ifdef BOOST_MATH_THROW_ON_UNDERFLOW_ERROR
  inline T underflow_error(const char* function, const char* message)
#else
  inline T underflow_error(const char*, const char*) // Warning suppressed version.
#endif
{
   errno = ERANGE;
#ifdef BOOST_MATH_THROW_ON_UNDERFLOW_ERROR
   detail::raise_error<std::underflow_error>(function, message ? message : "Underflow");
#endif
   return 0;
} // underflow_error

// Result is denormalised:
// so use domain_error to show message & return unchanged denormalised value.
template <class T>
#ifdef BOOST_MATH_THROW_ON_DENORM_ERROR
  inline T denorm_error(T const& t, const char* function, const char* message)
#else
  inline T denorm_error(T const& t, const char*, const char*) // Warning suppressed version.
#endif
{
   // Don't set errno if denormalized.
#ifdef BOOST_MATH_THROW_ON_DENORM_ERROR
   // We don't have a std exception type for this, but this feels about right:
   detail::raise_error<std::underflow_error>(function, message ? message : "Denormalized value");
#endif
   return t; // Unchanged denormalised value.
} // denorm_error
//
// Computed result is garbage / internal error:
template <class T>
inline T logic_error(const char* function, const char* message, const T& val)
{
   errno = EDOM;
   detail::raise_error<std::logic_error>(function, message ? message : "Internal logic error, computed value was %1%", val);
   return 0; // We don't get here:
} // logic_error

namespace detail
{
	template <class T, class U>
	inline T checked_narrowing_cast(U const& val, const char* , const boost::mpl::true_*)
	{
		 return val;
	}

	template <class T, class U>
	inline T checked_narrowing_cast(U const& val, const char* function, const boost::mpl::false_*)
	{
		 if(val > tools::max_value<T>())
				return tools::overflow_error<T>(function, 0);
		 T result = static_cast<T>(val);
		 if((result == 0) && (val != 0))
				return tools::underflow_error<T>(function, 0);
		 if((fabs(result) < tools::min_value<T>()) && (fabs(val) >= tools::min_value<U>()))
				return tools::denorm_error<T>(result, function, 0);
		 return result;
	}
} // namespace detail

template <class T, class U>
inline T checked_narrowing_cast(U const& val, const char* function)
{
#ifndef BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS
   typedef boost::mpl::bool_<
      ::boost::is_same<T, U>::value
      || (std::numeric_limits<T>::digits >= std::numeric_limits<U>::digits
      && std::numeric_limits<T>::min_exponent <= std::numeric_limits<U>::min_exponent
      && std::numeric_limits<T>::max_exponent >= std::numeric_limits<U>::max_exponent) > select_type;
#else
   typedef boost::mpl::bool_<
      ::boost::is_same<T, U>::value> select_type;
#endif
   return detail::checked_narrowing_cast<T>(val, function, static_cast<const select_type*>(0));
} // checked_narrowing_cast

} // namespace tools
} // namespace math
} // namespace boost

#ifdef BOOST_MSVC
#  pragma warning(pop)
#endif

#endif // BOOST_MATH_TOOLS_ERROR_HANDLING_HPP

