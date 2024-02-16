//  Copyright John Maddock 2014.
//  Copyright Christopher Kormanyos 2024.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#define TEST_COMPLEX

#include <boost/core/lightweight_test.hpp>

#include "compile_test/instantiate.hpp"

#include <limits>
#include <type_traits>

namespace local
{
  namespace detail
  {
    using float64_t =
        typename std::conditional<std::numeric_limits<float>::digits == 53,
                                  float,
                                  typename std::conditional<std::numeric_limits<double>::digits == 53,
                                                            double,
                                                            long double>::type>::type;
  } // namespace detail

  auto instantiate_runner() -> void
  {
    instantiate(static_cast<detail::float64_t>(1.23L));

    const bool result_instantiate_and_run_is_ok = instantiate_runner_result<detail::float64_t>::value;

    BOOST_TEST(result_instantiate_and_run_is_ok);
  }
} // namespace local

auto main() -> int
{
  static_assert(std::numeric_limits<local::detail::float64_t>::digits == 53, "Error: rong digits in local::float64_t");

  local::instantiate_runner();

  return boost::report_errors();
}
