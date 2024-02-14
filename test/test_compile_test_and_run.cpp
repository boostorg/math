//  Copyright John Maddock 2014.
//  Copyright Christopher Kormanyos 2024.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/core/lightweight_test.hpp>
#include <boost/cstdfloat.hpp>

#if !defined(__CYGWIN__)
#include "compile_test/instantiate.hpp"
#endif

namespace local
{
  auto instantiate_runner() -> void
  {
    #if (defined(BOOST_FLOAT64_C) && !defined(__CYGWIN__))

    instantiate(static_cast<boost::float64_t>(BOOST_FLOAT64_C(1.23)));

    const bool result_instantiate_and_run_is_ok = instantiate_runner_result<boost::float64_t>::value;

    BOOST_TEST(result_instantiate_and_run_is_ok);

    #else

    BOOST_TEST(true);

    #endif
  }
}

auto main() -> int
{
  local::instantiate_runner();

  const auto result_is_ok = (boost::report_errors() == 0);

  return (result_is_ok ? 0 : -1);
}
