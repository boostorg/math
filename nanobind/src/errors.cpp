// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// The user_error handlers of boost/math/policies/error_handling.hpp, and boost.math.errstate's state.

#include "bindings.hpp"

#include <boost/math/policies/error_handling.hpp>

#ifdef BMPY_HAVE_ULPS
#include "ulps.hpp"
#endif

#include <limits>
#include <stdexcept>
#include <string>

namespace bmpy {

namespace {

namespace bp = boost::math::policies;

// The message throw_on_error would give.
template <class T>
std::string format(char const* function, char const* message, T const& val) {
   std::string f(function ? function : "Unknown function operating on type %1%");
   std::string m(message ? message : "Cause unknown");
   bp::detail::replace_all_in_string(f, "%1%", bp::detail::name_of<T>());
   bp::detail::replace_all_in_string(m, "%1%", bp::detail::prec_format(val).c_str());
   return "Error in function " + f + ": " + m;
}

// Throws E, or records the error in category c if a recording is active.
template <class E, class T>
void handle(error_category c, char const* function, char const* message, T const& val) {
   error_record* r = active_record;
   if (r)
      r->failed[c] = true;
   if (r && !r->first[c].empty())
      return;
   std::string msg = format(function, message, val);
   if (!r)
      throw E(msg);
   r->first[c] = std::move(msg);
}

template <class T>
T nan() {
   return std::numeric_limits<T>::has_quiet_NaN ? std::numeric_limits<T>::quiet_NaN() : T(0);
}

}  // namespace

void error_record::merge(error_record const& other) {
   for (int c = 0; c < n_error_categories; ++c) {
      if (first[c].empty())
         first[c] = other.first[c];
      count[c] += other.count[c] + other.failed[c];
   }
}

void record_exception(std::exception const& e) {
   error_category c = dynamic_cast<std::overflow_error const*>(&e)                  ? overflow_errors
                      : dynamic_cast<boost::math::evaluation_error const*>(&e) ? evaluation_errors
                                                                                : domain_errors;
   active_record->failed[c] = true;
   if (active_record->first[c].empty())
      active_record->first[c] = e.what();
}

void report(char const* name, error_record record, std::size_t size, bool scalar, bool integer_result,
            error_actions const& actions) {
   static char const* what[] = {"domain error", "overflow", "evaluation error"};
   static char const* value[] = {"NaN", "an infinity", "NaN"};
   static PyObject* const* exception[] = {&PyExc_ValueError, &PyExc_OverflowError, &PyExc_RuntimeError};
   record.next();
   std::string warning, raise;
   int raised = -1;
   for (int c = 0; c < n_error_categories; ++c) {
      if (!record.count[c])
         continue;
      error_action a = actions[c];
      if (a == action_default)
         a = scalar || integer_result ? action_raise : c == overflow_errors ? action_ignore : action_warn;
      else if (a == action_warn && integer_result)
         a = action_raise;
      std::string counted = std::string(name) + ": " + what[c] + " in " + std::to_string(record.count[c]) + " of " +
                            std::to_string(size) + " evaluations";
      if (a == action_raise && raised < 0) {
         raised = c;
         raise = size == 1 ? record.first[c] : counted + "; first: " + record.first[c];
      }
      if (a == action_warn)
         warning += (warning.empty() ? "" : "; ") + counted + ", set to " + value[c] + "; first: " + record.first[c];
   }
   if (!warning.empty() && PyErr_WarnEx(PyExc_RuntimeWarning, warning.c_str(), 1) != 0)
      throw nb::python_error();
   if (raised >= 0) {
      PyErr_SetString(*exception[raised], raise.c_str());
      throw nb::python_error();
   }
}

void bind_errors(nb::module_& m) {
   m.def("_geterr", [] {
      error_actions const& a = current_error_actions;
      return nb::make_tuple(int(a[0]), int(a[1]), int(a[2]));
   });
   m.def("_seterr", [](int domain, int overflow, int evaluation) {
      current_error_actions = {error_action(domain), error_action(overflow), error_action(evaluation)};
   });
}

}  // namespace bmpy

namespace boost::math::policies {

template <class T>
T user_domain_error(char const* function, char const* message, T const& val) {
   bmpy::handle<std::domain_error>(bmpy::domain_errors, function, message, val);
   return bmpy::nan<T>();
}

template <class T>
T user_pole_error(char const* function, char const* message, T const& val) {
   bmpy::handle<std::domain_error>(bmpy::domain_errors, function, message, val);
   return bmpy::nan<T>();
}

// val is the infinity to return; callers negate it for negative overflows.
template <class T>
T user_overflow_error(char const* function, char const* message, T const& val) {
   bmpy::handle<std::overflow_error>(bmpy::overflow_errors, function, message, val);
   return val;
}

template <class T>
T user_evaluation_error(char const* function, char const* message, T const& val) {
   bmpy::handle<boost::math::evaluation_error>(bmpy::evaluation_errors, function, message, val);
   return bmpy::nan<T>();
}

template <class T>
T user_indeterminate_result_error(char const* function, char const* message, T const& val) {
   bmpy::handle<std::domain_error>(bmpy::domain_errors, function, message, val);
   return bmpy::nan<T>();
}

// As ignore_error: the largest value of the target type, with the sign of val.
template <class T, class TargetType>
TargetType user_rounding_error(char const* function, char const* message, T const& val, TargetType const&) {
   // Always thrown, as by the default policy: a value returned here would become a loop count or an index
   // (NaN to LLONG_MIN), which hangs or reads out of bounds. Inside a recording, the elementwise loop catches it.
   throw boost::math::rounding_error(bmpy::format(function, message, val));
}

#define BMPY_USER_ERRORS(T)                                                                \
   template T user_domain_error<T>(char const*, char const*, T const&);                    \
   template T user_pole_error<T>(char const*, char const*, T const&);                      \
   template T user_overflow_error<T>(char const*, char const*, T const&);                  \
   template T user_evaluation_error<T>(char const*, char const*, T const&);                \
   template T user_indeterminate_result_error<T>(char const*, char const*, T const&);
#define BMPY_USER_ROUNDING_ERRORS(T)                                                                      \
   template T user_rounding_error<T, T>(char const*, char const*, T const&, T const&);                       \
   template int user_rounding_error<T, int>(char const*, char const*, T const&, int const&);                 \
   template long user_rounding_error<T, long>(char const*, char const*, T const&, long const&);              \
   template long long user_rounding_error<T, long long>(char const*, char const*, T const&, long long const&);
// The types the bindings instantiate the error handlers with (prime and fibonacci raise in integers):
BMPY_USER_ERRORS(float)
BMPY_USER_ERRORS(double)
BMPY_USER_ERRORS(long double)
BMPY_USER_ERRORS(unsigned)
BMPY_USER_ERRORS(unsigned long)
BMPY_USER_ERRORS(unsigned long long)
BMPY_USER_ROUNDING_ERRORS(float)
BMPY_USER_ROUNDING_ERRORS(double)
BMPY_USER_ROUNDING_ERRORS(long double)
// The Hermite interpolators, makima and pchip (daubechies_scaling too) raise with a bool.
template bool user_domain_error<bool>(char const*, char const*, bool const&);
#ifdef BMPY_HAVE_ULPS
// boost.math.ulps's reference evaluations.
BMPY_USER_ERRORS(bmpy::ulps::Precise)
BMPY_USER_ROUNDING_ERRORS(bmpy::ulps::Precise)
#endif
#undef BMPY_USER_ERRORS
#undef BMPY_USER_ROUNDING_ERRORS

}  // namespace boost::math::policies
