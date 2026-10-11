//           Copyright Matthew Pulver 2018 - 2019.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)

#include "test_autodiff.hpp"
#include <boost/math/special_functions.hpp>


template <class T>
void binomial_hpp() {
  using boost::multiprecision::min;
  using std::fabs;
  using std::min;

  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;
  test_detail::RandomSample<unsigned> n_sampler{0u, 30u};
  test_detail::RandomSample<unsigned> r_sampler{0u, 30u};

  for (auto i : boost::irange(test_constants::n_samples)) {
    std::ignore = i;
    auto n = n_sampler.next();
    auto r = n == 0 ? 0 : (min)(r_sampler.next(), n - 1);

    // This is a hard function to test for type float due to a specialization of
    // boost::math::binomial_coefficient
    auto autodiff_v =
        std::is_same<T, float>::value
            ? make_fvar<T, m>(boost::math::binomial_coefficient<T>(n, r))
            : boost::math::binomial_coefficient<T>(n, r);
    auto anchor_v = boost::math::binomial_coefficient<T>(n, r);
    CHECK_EQUAL(autodiff_v.derivative(0u), anchor_v);
  }
}

template <class T>
void cbrt_hpp() {
  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;
  test_detail::RandomSample<T> x_sampler{-2000, 2000};
  for (auto i : boost::irange(test_constants::n_samples)) {
    std::ignore = i;
    auto x = x_sampler.next();
    CHECK_RELATIVE_ERROR(boost::math::cbrt(make_fvar<T, m>(x)).derivative(0u),
                      boost::math::cbrt(x), (50 * test_constants::pct_epsilon()) / 100);
  }
}

#if !defined(__APPLE__)
template <class T>
void chebyshev_hpp() {
  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;
  {
    test_detail::RandomSample<unsigned> n_sampler{0u, 10u};
    test_detail::RandomSample<T> x_sampler{-2, 2};
    for (auto i : boost::irange(test_constants::n_samples)) {
      std::ignore = i;
      auto n = n_sampler.next();
      auto x = x_sampler.next();
      CHECK_RELATIVE_ERROR(
          boost::math::chebyshev_t(n, make_fvar<T, m>(x)).derivative(0u),
          boost::math::chebyshev_t(n, x), (40 * test_constants::pct_epsilon()) / 100);
      // Lower accuracy with clang/apple:
      CHECK_RELATIVE_ERROR(
          boost::math::chebyshev_u(n, make_fvar<T, m>(x)).derivative(0u),
          boost::math::chebyshev_u(n, x), (80 * test_constants::pct_epsilon()) / 100);

      CHECK_RELATIVE_ERROR(
          boost::math::chebyshev_t_prime(n, make_fvar<T, m>(x)).derivative(0u),
          boost::math::chebyshev_t_prime(n, x), (40 * test_constants::pct_epsilon()) / 100);

      /*/usr/include/boost/math/special_functions/chebyshev.hpp:164:40: error:
       cannot convert
       boost::math::differentiation::autodiff_v1::detail::fvar<double, 3> to
       double in return
       CHECK_EQUAL(boost::math::chebyshev_clenshaw_recurrence(c.data(),c.size(),make_fvar<T,m>(0.20))
       ,
       boost::math::chebyshev_clenshaw_recurrence(c.data(),c.size(),static_cast<T>(0.20)));*/
      /*try {
        std::array<T, 4> c0{{14.2, -13.7, 82.3, 96}};
        CHECK_RELATIVE_ERROR(boost::math::chebyshev_clenshaw_recurrence(c0.data(),
      c0.size(), make_fvar<T,m>(x)),
                                     boost::math::chebyshev_clenshaw_recurrence(c0.data(),
      c0.size(), x), (10*test_constants::pct_epsilon()) / 100); } catch (...) {
        std::rethrow_exception(std::exception_ptr(std::current_exception()));
      }*/
    }
  }
}
#endif

template <class T>
void cospi_hpp() {
  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;
  test_detail::RandomSample<T> x_sampler{-2000, 2000};
  for (auto i : boost::irange(test_constants::n_samples)) {
    std::ignore = i;
    auto x = x_sampler.next();
    CHECK_RELATIVE_ERROR(boost::math::cos_pi(make_fvar<T, m>(x)).derivative(0u),
                      boost::math::cos_pi(x), (test_constants::pct_epsilon()) / 100);
  }
}

template <class T>
void digamma_hpp() {

  using boost::math::nextafter;
  using std::nextafter;

  using test_constants = test_constants_t<T>;
  static constexpr auto m = test_constants::order;
  test_detail::RandomSample<T> x_sampler{-1, 2000};
  for (auto i : boost::irange(test_constants::n_samples)) {
    std::ignore = i;
    auto x = nextafter(x_sampler.next(), ((std::numeric_limits<T>::max))());
    auto autodiff_v = boost::math::digamma(make_fvar<T, m>(x));
    auto anchor_v = boost::math::digamma(x);
    CHECK_RELATIVE_ERROR(autodiff_v.derivative(0u), anchor_v, (1e4 * test_constants::pct_epsilon()) / 100);
  }
}

int main()
{
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { binomial_hpp<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { cbrt_hpp<typename decltype(t)::type>(); });
#if !defined(__APPLE__)
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { chebyshev_hpp<typename decltype(t)::type>(); });
#endif
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { cospi_hpp<typename decltype(t)::type>(); });
   boost::mp11::mp_for_each<boost::mp11::mp_transform<boost::mp11::mp_identity, all_float_types>>([](auto t) { digamma_hpp<typename decltype(t)::type>(); });
   return boost::math::test::report_errors();
}
