// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// Boost.Math errors, with the user_error policy (see CMakeLists.txt). The handlers, in errors.cpp, throw
// the exceptions throw_on_error would, unless a `recording` is active on the calling thread: then they
// mark the error, keep its first message and return the natural value (NaN, or the signed infinity of an
// overflow). The elementwise functions record, and then report() applies the boost.math.errstate actions.

#ifndef BOOST_MATH_NANOBIND_ERRORS_HPP
#define BOOST_MATH_NANOBIND_ERRORS_HPP

#include <array>
#include <cstddef>
#include <exception>
#include <string>
#include <utility>

namespace bmpy {

// Pole, indeterminate-result and rounding errors count as domain errors.
enum error_category { domain_errors, overflow_errors, evaluation_errors, n_error_categories };

// What a category does; by default, raise for scalars, warn for arrays (ignore array overflows).
enum error_action : unsigned char { action_default, action_ignore, action_warn, action_raise };
using error_actions = std::array<error_action, n_error_categories>;
inline thread_local error_actions current_error_actions{};

// count is of the evaluations that failed; failed marks the categories of the current one, until next().
struct error_record {
   std::array<std::size_t, n_error_categories> count{};
   std::array<bool, n_error_categories> failed{};
   std::array<std::string, n_error_categories> first;
   void next() {
      for (int c = 0; c < n_error_categories; ++c)
         count[c] += std::exchange(failed[c], false);
   }
   bool any() const { return count[0] + count[1] + count[2] + failed[0] + failed[1] + failed[2] != 0; }
   void merge(error_record const& other);
};

inline thread_local error_record* active_record = nullptr;

// While alive, Boost.Math errors on this thread are recorded in r instead of thrown; recording(nullptr)
// suspends an outer recording, around a call of Python code.
struct recording {
   error_record* previous;
   explicit recording(error_record& r) : previous(active_record) { active_record = &r; }
   explicit recording(std::nullptr_t) : previous(active_record) { active_record = nullptr; }
   ~recording() { active_record = previous; }
   recording(recording const&) = delete;
   recording& operator=(recording const&) = delete;
};

// Records an exception that escaped the error policy (an integer argument out of range, ...).
void record_exception(std::exception const& e);

// Applies `actions` to what `record` holds for `size` evaluations of `name`: emits one RuntimeWarning for
// the categories whose action is warn, then raises the first whose action is raise. `scalar` selects the
// defaults; integer results, which have no NaN to mark a failure, raise instead of warning.
void report(char const* name, error_record record, std::size_t size, bool scalar, bool integer_result,
            error_actions const& actions);

}  // namespace bmpy

#endif
