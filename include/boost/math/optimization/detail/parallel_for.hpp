//  Copyright Nick Thompson, 2026
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#ifndef BOOST_MATH_OPTIMIZATION_DETAIL_PARALLEL_FOR_HPP
#define BOOST_MATH_OPTIMIZATION_DETAIL_PARALLEL_FOR_HPP
#include <boost/math/tools/config.hpp>
#ifndef BOOST_MATH_BUILD_MODULE
#include <algorithm>
#include <cstddef>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>
#endif

BOOST_MATH_NAMESPACE_BEGIN namespace optimization::detail {

// Calls f(i) for i < count on up to `threads` threads, of which the calling thread is one. The first exception thrown,
// by f or by the start of a thread, is rethrown once every thread has been joined.
template <class F> void parallel_for(size_t count, size_t threads, F &&f) {
  threads = (std::max)(size_t(1), (std::min)(threads, count));
  std::mutex mutex;
  std::exception_ptr error;
  auto guarded = [&](auto &&g) {
    try {
      g();
    } catch (...) {
      std::lock_guard<std::mutex> lock(mutex);
      error = error ? error : std::current_exception();
    }
  };
  auto work = [&](size_t j) {
    guarded([&] {
      for (size_t i = j; i < count; i += threads) {
        f(i);
      }
    });
  };
  std::vector<std::thread> pool;
  guarded([&] {
    for (size_t j = 1; j < threads; ++j) {
      pool.emplace_back(work, j);
    }
  });
  work(0);
  for (auto &thread : pool) {
    thread.join();
  }
  if (error) {
    std::rethrow_exception(error);
  }
}

} BOOST_MATH_NAMESPACE_END // namespace boost::math::optimization::detail
#endif
