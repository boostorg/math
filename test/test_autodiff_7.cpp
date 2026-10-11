//           Copyright Matthew Pulver 2018 - 2019.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)

#include "test_autodiff.hpp"


template <class T>
void expm1_hpp() {
  using boost::math::differentiation::detail::log;
  using boost::multiprecision::log;
  using std::log;
  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;
  test_detail::RandomSample<T> x_sampler{-log(T(2000)), log(T(2000))};
  for (auto i : boost::irange(test_constants::n_samples)) {
    std::ignore = i;
    auto x = x_sampler.next();
    CHECK_RELATIVE_ERROR(boost::math::expm1(make_fvar<T, m>(x)).derivative(0u),
                      boost::math::expm1(x), (50 * test_constants::pct_epsilon()) / 100);
  }
}

int main()
{
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { expm1_hpp<typename decltype(t)::type>(); });
   return boost::math::test::report_errors();
}
